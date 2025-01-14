//
//  CPXOPF.cpp
//  Gravity
//
//  Created by Hassan Hijazi on Nov 17 2018
//
//
#include <stdio.h>
#include <PowerNet.h>
#include <gravity/solver.h>



using namespace std;
using namespace gravity;



int main (int argc, char * argv[])
{
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m", mtype = "polar";
    int output = 0;
    bool relax = false;
    double tol = 1e-6;
    string log_level="0";
    
    /** Read user options */
    auto options = readOptions(argc, argv);
    options.add_option("f", "file", "Input file name (def. ../../data_sets/Power/nesta_case5_pjm.m)", fname );
    options.add_option("m", "model", "Complex Number Model: polar/rect (def. polar)", mtype );
    /** parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = options.parse_options(argc, argv);
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    fname = options["f"];
    mtype = options["m"];
    output = op::str2int(options["l"]);
    bool has_help = op::str2bool(options["h"]);
    /** show help */
    if(has_help) {
        options.show_help();
        exit(0);
    }
    /* Start timer */
    auto total_time_start = get_wall_time();
    
    /* Construct Grid **/
    PowerNet grid;
    grid.readgrid(fname);
    
    /** Sets */
    auto nodes = indices(grid.nodes);
    auto arcs = indices(grid.arcs);
    auto gens = indices(grid.gens);
    auto gen_nodes = grid.gens_per_node();
    auto out_arcs = grid.out_arcs_per_node();
    auto in_arcs = grid.in_arcs_per_node();
    
    /* Grid Parameters */
    unsigned nb_gen = grid.get_nb_active_gens();
    unsigned nb_lines = grid.get_nb_active_arcs();
    unsigned nb_buses = grid.get_nb_active_nodes();
    /* Line properties */
    auto Y = grid.Y.in(arcs);
    /* Line charging */
    auto Ych = grid.Ych.in(arcs);
    /* Line limit */
    auto S_max = grid.S_max.in(arcs);
    /* Tap ratio */
    auto T = grid.T.in(arcs);
    /* Bus shunt */
    auto Ysh = grid.Ysh.in(nodes);
    /* Power Demand */
    auto Sd = grid.Sd.in(nodes);
    /* Generation Bounds */
    auto Sg_min = grid.Sg_min.in(gens);
    auto Sg_max = grid.Sg_max.in(gens);
    /* Voltage Bounds */
    auto V_min = grid.V_min.in(nodes);
    auto V_max = grid.V_max.in(nodes);
    /* Phase angle Bounds */
    auto th_min = grid.th_min.in(arcs);
    auto th_max = grid.th_max.in(arcs);
    /* Generators' costs */
    auto c0 = grid.c0.in(gens);
    auto c1 = grid.c1.in(gens);
    auto c2 = grid.c2.in(gens);
    /* Reference Bus */
    auto ref_bus = grid.ref_bus;
    
    
    DebugOn("nb gens = " << nb_gen << endl);
    DebugOn("nb lines = " << nb_lines << endl);
    DebugOn("nb buses = " << nb_buses << endl);
    
    PowerModelType pmt = ACPOL;
    if(mtype=="rect") pmt = ACRECT;
    bool polar = (pmt==ACPOL);
    if (polar) {
        DebugOn("Using polar model\n");
    }
    else {
        DebugOn("Using rectangular model\n");
    }
    Model<Cpx> CPXOPF("CPX-OPF Model");
    /** Variables */
    /* Power generation variables */
    var<Cpx> Sg("Sg", Sg_min, Sg_max);
    CPXOPF.add(Sg.in(gens));
    
    /* Power flow variables */
    var<Cpx> S_from("S_from", -1*S_max,S_max);
    var<Cpx> S_to("S_to", -1*S_max,S_max);
    CPXOPF.add(S_from.in(arcs),S_to.in(arcs));
    
    /** Voltage variables */
    var<Cpx> V("V", V_min, V_max);
    CPXOPF.add(V.in(nodes));
    
    auto V_from = V.from(arcs);
    auto V_to = V.to(arcs);
    
    /** Construct the objective function */
    auto obj = c1.tr()*real(Sg) + c2.tr()*pow(real(Sg),2) + sum(c0);
    CPXOPF.min(obj);
    
    /** Define constraints */
    
    /* Ref Bus */
    Constraint<Cpx> Ref_Bus("Ref_Bus");
    Ref_Bus = imag(V)(ref_bus);
    CPXOPF.add(Ref_Bus == 0);
    
    /** KCL Flow conservation */
    Constraint<Cpx> KCL("KCL");
    KCL  = sum(S_from, out_arcs) + sum(S_to, in_arcs) + Sd - sum(Sg, gen_nodes) + conj(Ysh)*sqrmag(V);
    CPXOPF.add(KCL.in(nodes) == 0);
    
    /** AC Power Flows */
    Constraint<Cpx> Flow_From("Flow_From");
    Flow_From = S_from - ((conj(Y)+conj(Ych))/sqrmag(T))*sqrmag(V_from) + (conj(Y)/T)*V_from*conj(V_to);
    CPXOPF.add(Flow_From.in(arcs)==0);
    
    Constraint<Cpx> Flow_To("Flow_To");
    Flow_To = S_to - (conj(Y)+conj(Ych))*sqrmag(V_to) + (conj(Y)/conj(T))*conj(V_from)*V_to;
    CPXOPF.add(Flow_To.in(arcs)==0);
    
    /* Phase Angle Bounds constraints */
    auto bus_pairs = grid.get_bus_pairs();
    Constraint<Cpx> PAD_UB("PAD_UB");
    PAD_UB = ang(V.from(bus_pairs)) - ang(V.to(bus_pairs));
    PAD_UB -= th_max;
    CPXOPF.add(PAD_UB.in(bus_pairs) <= 0);
    Constraint<Cpx> PAD_LB("PAD_LB");
    PAD_LB = ang(V.from(bus_pairs)) - ang(V.to(bus_pairs));
    PAD_LB -= th_min;
    CPXOPF.add(PAD_LB.in(bus_pairs) >= 0);
    
    /*  Thermal Limit Constraints */
    Constraint<Cpx> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = sqrmag(S_from) - pow(S_max, 2);
    CPXOPF.add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint<Cpx> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = sqrmag(S_to) - pow(S_max,2);
    CPXOPF.add(Thermal_Limit_to.in(arcs) <= 0);
    CPXOPF.print_symbolic();
    CPXOPF.print();
    return 0;
//    solver<Cpx> OPF(CPXOPF,ipopt);
    auto solver_time_start = get_wall_time();
//    OPF.run(output, tol = 1e-6);
    auto solver_time_end = get_wall_time();
    auto total_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    auto total_time = total_time_end - total_time_start;
    /* Uncomment line below to print expanded model */
    
    
    /* CPXOPF.print(); */
    
    /** Terminal output **/
    string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(CPXOPF.get_obj_val().real()) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
    DebugOn(out <<endl);
    
    return 0;
}
