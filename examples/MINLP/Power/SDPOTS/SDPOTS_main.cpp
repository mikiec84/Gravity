//
// Created by ksenia on 21.11.18.
//

#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <optionParser.hpp>

using namespace std;
using namespace gravity;

/* SDP OPF constraints */

void add_PowerFlow(Model& model, PowerNet& grid, var<Real>& R_Wij, var<Real>& Im_Wij, var<Real>& Wii,
        var<Real>& Pf_from, var<Real>& Pf_to, var<Real>& Qf_from, var<Real>& Qf_to, var<Real>& Pg, var<Real>& Qg){
   /* Flow conservation */
   Constraint KCL_P("KCL_P");
   KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid.pl - sum(Pg.in_gens()) + grid.gs*Wii;
   model.add(KCL_P.in(grid.nodes) == 0);

   Constraint KCL_Q("KCL_Q");
   KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens()) - grid.bs*Wii;
   model.add(KCL_Q.in(grid.nodes) == 0);

   /* AC Power Flow */
   Constraint Flow_P_From("Flow_P_From");
   Flow_P_From = Pf_from - (grid.g_ff*Wii.from() + grid.g_ft*R_Wij.in_pairs() + grid.b_ft*Im_Wij.in_pairs());
   model.add(Flow_P_From.in(grid.arcs) == 0);

   Constraint Flow_P_To("Flow_P_To");
   Flow_P_To = Pf_to - (grid.g_tt*Wii.to() + grid.g_tf*R_Wij.in_pairs() - grid.b_tf*Im_Wij.in_pairs());
   model.add(Flow_P_To.in(grid.arcs) == 0);

   Constraint Flow_Q_From("Flow_Q_From");
   Flow_Q_From = Qf_from - (grid.g_ft*Im_Wij.in_pairs() - grid.b_ff*Wii.from() - grid.b_ft*R_Wij.in_pairs());
   model.add(Flow_Q_From.in(grid.arcs) == 0);

   Constraint Flow_Q_To("Flow_Q_To");
   Flow_Q_To = Qf_to + (grid.b_tt*Wii.to() + grid.b_tf*R_Wij.in_pairs() + grid.g_tf*Im_Wij.in_pairs());
   model.add(Flow_Q_To.in(grid.arcs) == 0);
}

void add_AngleBounds(Model& model, PowerNet& grid, var<Real>& R_Wij, var<Real>& Im_Wij){
   auto bus_pairs = grid.get_bus_pairs();

   Constraint PAD_UB("PAD_UB");
   PAD_UB = Im_Wij;
   PAD_UB <= grid.tan_th_max*R_Wij;
//    SDP.add_lazy(PAD_UB.in(bus_pairs));
   model.add(PAD_UB.in(bus_pairs));

   Constraint PAD_LB("PAD_LB");
   PAD_LB =  Im_Wij;
   PAD_LB >= grid.tan_th_min*R_Wij;
//    SDP.add_lazy(PAD_LB.in(bus_pairs));
   model.add(PAD_LB.in(bus_pairs));
}

void add_Thermal(Model& model, PowerNet& grid, var<Real>& Pf_from, var<Real>& Pf_to, var<Real>& Qf_from, var<Real>& Qf_to){
   Constraint Thermal_Limit_from("Thermal_Limit_from");
   Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
   Thermal_Limit_from <= power(grid.S_max,2);
   model.add(Thermal_Limit_from.in(grid.arcs));

   Constraint Thermal_Limit_to("Thermal_Limit_to");
   Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
   Thermal_Limit_to <= power(grid.S_max,2);
   model.add(Thermal_Limit_to.in(grid.arcs));
}

void add_SOC(Model& model, PowerNet& grid, var<Real>& R_Wij, var<Real>& Im_Wij, var<Real>& Wii){
   auto bus_pairs_chord = grid.get_bus_pairs_chord();

   /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to();
    model.add(SOC.in(bus_pairs_chord) <= 0);
}

void add_LNC(Model& model, PowerNet& grid, var<Real>& R_Wij, var<Real>& Im_Wij, var<Real>& Wii){
   auto bus_pairs = grid.get_bus_pairs();

   /* Lifted Nonlinear Cuts */
   Constraint LNC1("LNC1");
   LNC1 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
   LNC1 -= grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
   LNC1 -= grid.v_max.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
   LNC1 -= grid.v_max.from()*grid.v_max.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    model.add(LNC1.in(bus_pairs) >= 0);
//    model.add_lazy(LNC1.in(bus_pairs) >= 0);

   Constraint LNC2("LNC2");
   LNC2 = (grid.v_min.from()+grid.v_max.from())*(grid.v_min.to()+grid.v_max.to())*(sin(0.5*(grid.th_max+grid.th_min))*Im_Wij + cos(0.5*(grid.th_max+grid.th_min))*R_Wij);
   LNC2 -= grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.to()+grid.v_max.to())*Wii.from();
   LNC2 -= grid.v_min.from()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()+grid.v_max.from())*Wii.to();
   LNC2 += grid.v_min.from()*grid.v_min.to()*cos(0.5*(grid.th_max-grid.th_min))*(grid.v_min.from()*grid.v_min.to() - grid.v_max.from()*grid.v_max.to());
//    model.add(LNC2.in(bus_pairs) >= 0);
//    model.add_lazy(LNC2.in(bus_pairs) >= 0);
}

void add_SDP3(Model& model, PowerNet& grid, var<Real>& R_Wij, var<Real>& Im_Wij, var<Real>& Wii){
   auto bags_copy = grid._bags;
   auto bag_size = grid._bags.size();
   DebugOn("\nNum of bags = " << bag_size << endl);
   DebugOn("Adding 3d determinant polynomial cuts\n");
   auto R_Wij_ = R_Wij.pairs_in_directed(grid, bags_copy, 3);
   auto Im_Wij_ = Im_Wij.pairs_in_directed(grid, bags_copy, 3);
   auto Wii_ = Wii.in_bags(bags_copy, 3);

   Constraint SDP3("SDP_3D");
   SDP3 = 2 * R_Wij_[0] * (R_Wij_[1] * R_Wij_[2] + Im_Wij_[1] * Im_Wij_[2]);
   SDP3 -= 2 * Im_Wij_[0] * (R_Wij_[2] * Im_Wij_[1] - Im_Wij_[2] * R_Wij_[1]);
   SDP3 -= (power(R_Wij_[0], 2) + power(Im_Wij_[0], 2)) * Wii_[2];
   SDP3 -= (power(R_Wij_[1], 2) + power(Im_Wij_[1], 2)) * Wii_[0];
   SDP3 -= (power(R_Wij_[2], 2) + power(Im_Wij_[2], 2)) * Wii_[1];
   SDP3 += Wii_[0] * Wii_[1] * Wii_[2];
   model.add(SDP3 >= 0);

   DebugOn("\nsdp 3d nb inst = " << SDP3.get_nb_instances() << endl);
   //        SDP3.print_expanded();
}

/* OTS constraints */

void add_PowerFlow_onoff(Model& model, PowerNet& grid, var<Real>& R_Wij, var<Real>& Im_Wij, var<Real>& Wii, var<Real>& W_line_ij,
                         var<Real>& W_line_ji, var<Real>& Pf_from, var<Real>& Pf_to, var<Real>& Qf_from, var<Real>& Qf_to, var<Real>& Pg, var<Real>& Qg){
   /* Flow conservation */
   Constraint KCL_P("KCL_P");
   KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid.pl - sum(Pg.in_gens()) + grid.gs*Wii;
   model.add(KCL_P.in(grid.nodes) == 0);

   Constraint KCL_Q("KCL_Q");
   KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid.ql - sum(Qg.in_gens()) - grid.bs*Wii;
   model.add(KCL_Q.in(grid.nodes) == 0);

   /* AC Power Flow */
   Constraint Flow_P_From("Flow_P_From");
   Flow_P_From = Pf_from - (grid.g_ff*W_line_ij + grid.g_ft*R_Wij.in_pairs() + grid.b_ft*Im_Wij.in_pairs());
   model.add(Flow_P_From.in(grid.arcs) == 0);
   Flow_P_From.print_expanded();

   Constraint Flow_P_To("Flow_P_To");
   Flow_P_To = Pf_to - (grid.g_tt*W_line_ji + grid.g_tf*R_Wij.in_pairs() - grid.b_tf*Im_Wij.in_pairs());
   model.add(Flow_P_To.in(grid.arcs) == 0);

   Constraint Flow_Q_From("Flow_Q_From");
   Flow_Q_From = Qf_from - (grid.g_ft*Im_Wij.in_pairs() - grid.b_ff*W_line_ij - grid.b_ft*R_Wij.in_pairs());
   model.add(Flow_Q_From.in(grid.arcs) == 0);

   Constraint Flow_Q_To("Flow_Q_To");
   Flow_Q_To = Qf_to + (grid.b_tt*W_line_ji + grid.b_tf*R_Wij.in_pairs() + grid.g_tf*Im_Wij.in_pairs());
   model.add(Flow_Q_To.in(grid.arcs) == 0);
}

void add_Thermal_onoff(Model& model, PowerNet& grid, var<Real>& Pf_from, var<Real>& Pf_to, var<Real>& Qf_from, var<Real>& Qf_to,
                 var<Binary>& z){
   Constraint Thermal_Limit_from("Thermal_Limit_from");
   Thermal_Limit_from = power(Pf_from, 2) + power(Qf_from, 2);
   Thermal_Limit_from <= power(grid.S_max,2)*z*z;
   model.add(Thermal_Limit_from.in(grid.arcs));

   Constraint Thermal_Limit_to("Thermal_Limit_to");
   Thermal_Limit_to = power(Pf_to, 2) + power(Qf_to, 2);
   Thermal_Limit_to <= power(grid.S_max,2)*z*z;
   model.add(Thermal_Limit_to.in(grid.arcs));
}

/* constraints defining on/off squared voltage variables:
 * W_line_ij = Wii if z = 1,
 * W_line_ij = 0   if z = 0
 */
void add_Wline_W(Model& model, PowerNet& grid, var<Real>& W_line_ij, var<Real>& W_line_ji, var<Real>& Wii, var<Binary>& z){
   Constraint W_line_ij_le("W_line_ij_le");
   W_line_ij_le = W_line_ij - Wii.from();
   W_line_ij_le += grid.w_min.from()*(1-z);
   model.add(W_line_ij_le.in(grid.arcs) <= 0);
   W_line_ij_le.print_expanded();

   Constraint W_line_ij_ge("W_line_ij_ge");
   W_line_ij_ge = Wii.from() - W_line_ij;
   W_line_ij_ge -= grid.w_max.from()*(1-z);
   model.add(W_line_ij_ge.in(grid.arcs) <= 0);
   W_line_ij_ge.print_expanded();

   Constraint W_line_ji_le("W_line_ji_le");
   W_line_ji_le = W_line_ji - Wii.to();
   W_line_ji_le += grid.w_min.to()*(1-z);
   model.add(W_line_ji_le.in(grid.arcs) <= 0);
   W_line_ji_le.print_expanded();

   Constraint W_line_ji_ge("W_line_ji_ge");
   W_line_ji_ge = Wii.to() - W_line_ji;
   W_line_ji_ge -= grid.w_max.to()*(1-z);
   model.add(W_line_ji_ge.in(grid.arcs) <= 0);
}

/* second-order cone constraints, on/off version */
void add_SOC_onoff(Model& model, PowerNet& grid, var<Real>& R_Wij, var<Real>& Im_Wij, var<Real>& Wii, var<Binary>& z){
   Constraint SOC("SOC");
   SOC = power(R_Wij.in_pairs(), 2) + power(Im_Wij.in_pairs(), 2) - Wii.from()*Wii.to();
   model.add(SOC.in(grid.arcs) <= 0);

   Constraint SOC_bis("SOC_bis");
   SOC_bis = power(R_Wij.in_pairs(), 2) + power(Im_Wij.in_pairs(), 2) - z*Wii.from()*grid.w_max.to();
   model.add(SOC_bis.in(grid.arcs) <= 0);

   Constraint SOC_bis2("SOC_bis2");
   SOC_bis2 = power(R_Wij.in_pairs(), 2) + power(Im_Wij.in_pairs(), 2) - z*Wii.to()*grid.w_max.from();
   model.add(SOC_bis2.in(grid.arcs) <= 0);
}

/* on/off variable bounds of the form:
 * x <= xub*z
 * x >= xlb*z
 */
void add_onoff_bnds(Model& model, PowerNet& grid, var<Real>& R_Wij, var<Real>& Im_Wij, var<Real>& W_line_ij, var<Real>& W_line_ji, var<Binary>& z){
   Constraint R_Wij_ub("R_Wij_ub");
   R_Wij_ub = R_Wij.in_pairs() - grid.wr_max.in_pairs()*z;
   model.add(R_Wij_ub.in(grid.arcs) <= 0);

   Constraint R_Wij_lb("R_Wij_lb");
   R_Wij_lb = R_Wij.in_pairs() - grid.wr_min.in_pairs()*z;
   model.add(R_Wij_lb.in(grid.arcs) >= 0);

   Constraint Im_Wij_ub("Im_Wij_ub");
   Im_Wij_ub = Im_Wij.in_pairs() - grid.wi_max.in_pairs()*z;
   model.add(Im_Wij_ub.in(grid.arcs) <= 0);

   Constraint Im_Wij_lb("Im_Wij_lb");
   Im_Wij_lb = Im_Wij.in_pairs() - grid.wi_min.in_pairs()*z;
   model.add(Im_Wij_lb.in(grid.arcs) >= 0);

   Constraint W_line_ij_ub("W_line_ij_ub");
   W_line_ij_ub = W_line_ij.in_pairs() - grid.w_max.from()*z;
   model.add(W_line_ij_ub.in(grid.arcs) <= 0);

   Constraint W_line_ij_lb("W_line_ij_lb");
   W_line_ij_lb = W_line_ij.in_pairs() - grid.w_min.from()*z;
   model.add(W_line_ij_lb.in(grid.arcs) >= 0);

   Constraint W_line_ji_ub("W_line_ji_ub");
   W_line_ji_ub = W_line_ji.in_pairs() - grid.w_max.to()*z;
   model.add(W_line_ji_ub.in(grid.arcs) <= 0);

   Constraint W_line_ji_lb("W_line_ji_lb");
   W_line_ji_lb = W_line_ji.in_pairs() - grid.w_min.to()*z;
   model.add(W_line_ji_lb.in(grid.arcs) >= 0);
}

/* main */
int main (int argc, char * argv[]) {
   int output = 0;
   bool relax = false, sdp_cuts = true;
   size_t num_bags = 0;
   string num_bags_s = "100";
   string solver_str = "ipopt";
   string sdp_cuts_s = "yes";
   SolverType solv_type = ipopt;
   double tol = 1e-6;
   double gap, solver_time_start;
   string mehrotra = "no";
   string fname = "../data_sets/Power/nesta_case5_pjm.m";

   // create a OptionParser with options
   op::OptionParser opt;
   opt.add_option("h", "help",
                  "shows option help"); // no default value means boolean options, which default value is false
   opt.add_option("f", "file", "Input file name", fname);
   opt.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
   opt.add_option("b", "numbags", "Number of bags per iteration", num_bags_s);
   opt.add_option("c", "sdpcuts", "Generate 3d SDP cuts, default = yes", sdp_cuts_s);
   // parse the options and verify that all went well. If not, errors and help will be shown
   bool correct_parsing = opt.parse_options(argc, argv);

   if (!correct_parsing) {
      return EXIT_FAILURE;
   }

   fname = opt["f"];
   bool has_help = op::str2bool(opt["h"]);
   if (has_help) {
      opt.show_help();
      exit(0);
   }
   solver_str = opt["s"];
   if (solver_str.compare("gurobi")==0) {
      solv_type = gurobi;
   }
   else if(solver_str.compare("cplex")==0) {
      solv_type = cplex;
   }else if(solver_str.compare("Mosek")==0) {
      solv_type = Mosek;
   }
   solver_str = opt["c"];
   if (solver_str.compare("no")==0) {
      sdp_cuts = false;
   }

   num_bags = atoi(opt["b"].c_str());

   cout << "\nnum bags = " << num_bags << endl;

   double total_time_start = get_wall_time();
   PowerNet grid;
   grid.readgrid(fname.c_str());

   grid.get_tree_decomp_bags();
   grid.update_net();

   // Grid Parameters
   auto bus_pairs = grid.get_bus_pairs();
   auto bus_pairs_chord = grid.get_bus_pairs_chord();
   auto nb_gen = grid.get_nb_active_gens();
   auto nb_lines = grid.get_nb_active_arcs();
   auto nb_buses = grid.get_nb_active_nodes();
   DebugOff("nb gens = " << nb_gen << endl);
   DebugOff("nb lines = " << nb_lines << endl);
   DebugOff("nb buses = " << nb_buses << endl);
   DebugOff("nb bus_pairs = " << nb_bus_pairs_chord << endl);

   double upper_bound = grid.solve_acopf();

   /** Build model */

   Model SDP("SDP Model");

   /* power generation variables */
   var<Real> Pg("Pg", grid.pg_min, grid.pg_max);
   var<Real> Qg ("Qg", grid.qg_min, grid.qg_max);
   SDP.add(Pg.in(grid.gens));
   SDP.add(Qg.in(grid.gens));
   /* power flow variables */
   var<Real> Pf_from("Pf_from", grid.S_max);
   var<Real> Qf_from("Qf_from", grid.S_max);
   var<Real> Pf_to("Pf_to", grid.S_max);
   var<Real> Qf_to("Qf_to", grid.S_max);
   SDP.add(Pf_from.in(grid.arcs));
   SDP.add(Qf_from.in(grid.arcs));
   SDP.add(Pf_to.in(grid.arcs));
   SDP.add(Qf_to.in(grid.arcs));

   /* Real part of Wij = ViVj */
   var<Real>  R_Wij("R_Wij", grid.wr_min, grid.wr_max);
   /* Imaginary part of Wij = ViVj */
   var<Real>  Im_Wij("Im_Wij", grid.wi_min, grid.wi_max);
   /* Magnitude of Wii = Vi^2 */
   var<Real>  Wii("Wii", grid.w_min, grid.w_max);

   SDP.add(Wii.in(grid.nodes));
   SDP.add(R_Wij.in(bus_pairs_chord));
   SDP.add(Im_Wij.in(bus_pairs_chord));

   /* Initialize variables */
   R_Wij.initialize_all(1.0);
   Wii.initialize_all(1.001);

   /**  Objective */
   auto obj = product(grid.c1, Pg) + product(grid.c2, power(Pg,2)) + sum(grid.c0);
   SDP.min(obj.in(grid.gens));

   /** Constraints */

   add_SOC(SDP, grid, R_Wij, Im_Wij, Wii);
   add_LNC(SDP, grid, R_Wij, Im_Wij, Wii);
   add_PowerFlow(SDP, grid, R_Wij, Im_Wij, Wii, Pf_from, Pf_to, Qf_from, Qf_to, Pg, Qg);
   add_Thermal(SDP, grid, Pf_from, Pf_to, Qf_from, Qf_to);
   add_AngleBounds(SDP, grid, R_Wij, Im_Wij);
//   add_SDP3(SDP, grid, R_Wij, Im_Wij, Wii);

   /* Solver selection */
   solver SDPOPF(SDP,solv_type);
   solver_time_start = get_wall_time();

   /* Solve the NLP */
   SDPOPF.run(output = 5, relax = false);

   /* after solving, update the variable values */
   for(auto& arc: grid.arcs){
      ((Line*)arc)->wr = R_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
      ((Line*)arc)->wi = Im_Wij(arc->_src->_name+","+arc->_dest->_name).eval();
//      if(grid.add_3d_nlin && !arc->_free) arc->_active = true;
   }
   for(auto& node: grid.nodes) ((Bus*)node)->w = Wii(node->_name).eval();

   gap = 100*(upper_bound - SDP._obj_val)/upper_bound;

   /* TODO build OTS */

   Model OTS("SDPOTS Model");

   /* power generation variables */
   var<Real> Pg1("Pg1", grid.pg_min, grid.pg_max);
   var<Real> Qg1("Qg1", grid.qg_min, grid.qg_max);
   OTS.add(Pg1.in(grid.gens));
   OTS.add(Qg1.in(grid.gens));
   /* power flow variables */
   var<Real> Pf_from1("Pf_from1", grid.S_max);
   var<Real> Qf_from1("Qf_from1", grid.S_max);
   var<Real> Pf_to1("Pf_to1", grid.S_max);
   var<Real> Qf_to1("Qf_to1", grid.S_max);
   OTS.add(Pf_from1.in(grid.arcs));
   OTS.add(Qf_from1.in(grid.arcs));
   OTS.add(Pf_to1.in(grid.arcs));
   OTS.add(Qf_to1.in(grid.arcs));

   /* Real part of Wij = ViVj */
   var<Real>  R_Wij1("R_Wij1", grid.wr_min, grid.wr_max);
   /* Imaginary part of Wij = ViVj */
   var<Real>  Im_Wij1("Im_Wij1", grid.wi_min, grid.wi_max);
   /* Magnitude of Wii = Vi^2 */
   var<Real>  Wii1("Wii1", grid.w_min, grid.w_max);
   /* On/off squared voltage magnitudes */
   var<Real>  W_line_ij("W_line_ij", 0.0, grid.w_max.from());
   var<Real>  W_line_ji("W_line_ji", 0.0, grid.w_max.to());
   /* Variales indicating whether a line is switched on */
   var<Binary> z("z", 1.0, 1.0);

   OTS.add(Wii1.in(grid.nodes));
   OTS.add(R_Wij1.in(bus_pairs_chord));
   OTS.add(Im_Wij1.in(bus_pairs_chord));
   OTS.add(W_line_ij.in(grid.arcs));
   OTS.add(W_line_ji.in(grid.arcs));
   OTS.add(z.in(grid.arcs));

   W_line_ij.print(true);
   z.print(true);

   /* Initialize variables */
   R_Wij1.initialize_all(1.0);
   Wii1.initialize_all(1.001);
   W_line_ij.initialize_all(1.001);
   z.initialize_all(1);

   /**  Objective */
   auto obj1 = product(grid.c1, Pg1) + product(grid.c2, power(Pg1,2)) + sum(grid.c0);
   OTS.min(obj1.in(grid.gens));

   /** Constraints */

   /* TODO OTS constraints */
   add_SOC_onoff(OTS, grid, R_Wij1, Im_Wij1, Wii1, z);
   add_LNC(OTS, grid, R_Wij1, Im_Wij1, Wii1);
   add_PowerFlow_onoff(OTS, grid, R_Wij1, Im_Wij1, Wii1, W_line_ij, W_line_ji, Pf_from1, Pf_to1, Qf_from1, Qf_to1, Pg1, Qg1);
   add_Thermal_onoff(OTS, grid, Pf_from1, Pf_to1, Qf_from1, Qf_to1, z);
   add_AngleBounds(OTS, grid, R_Wij1, Im_Wij1);
   add_Wline_W(OTS, grid, W_line_ij, W_line_ji, Wii1, z);
   add_onoff_bnds(OTS, grid, R_Wij1, Im_Wij1, W_line_ij, W_line_ji, z);
//   add_SDP_onoff(OTS, grid, R_Wij1, Im_Wij1, Wii1);

   /* TODO generate cuts at the optimal solution of the NLP */

   /* Solver selection */
   solver SDPOTS(OTS, gurobi);

   /* Solve OTS */
   SDPOTS.run(output = 5, relax = false);

   for(auto& arc: grid.arcs){
      cout << "\nz(" << arc->_name << ") = " << z(arc->_name).eval();
      cout << "\nw_line_ij(" << arc->_name << ") = " << W_line_ij(arc->_name).eval();
      cout << "\nR_Wij(" << arc->_name << ") = " << R_Wij1(arc->_src->_name+","+arc->_dest->_name).eval();
      cout << "\nIm_Wij(" << arc->_name << ") = " << Im_Wij1(arc->_src->_name+","+arc->_dest->_name).eval();
   }

   double solver_time_end = get_wall_time();
   double total_time_end = get_wall_time();
   auto solve_time = solver_time_end - solver_time_start;
   auto total_time = total_time_end - total_time_start;
   string out = "\nDATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(SDP._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
   DebugOn(out <<endl);
   DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
   DebugOn("Upper bound = " << to_string(upper_bound) << "."<<endl);
   DebugOn("Lower bound = " << to_string(SDP._obj_val) << "."<<endl);
   DebugOn("\nResults: " << grid._name << " " << to_string(SDP._obj_val) << " " << to_string(total_time)<<endl);
   return 0;
}
