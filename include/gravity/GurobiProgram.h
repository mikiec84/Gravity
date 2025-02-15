#ifndef __Gravity____GurobiProgram
#define __Gravity____GurobiProgram

#ifdef USE_GUROBI
#include <gurobi_c++.h>
#endif
#include <gravity/model.h>
using namespace gravity;


class GurobiProgram: public Program<>{
private:

    GRBModel* grb_mod;
    GRBEnv* grb_env;
    vector<GRBVar> _grb_vars; /** Mapping variables to Gurobi variables */
public:
    shared_ptr<Model<>>  _model;
    int     _output;
    GurobiProgram();
    GurobiProgram(const shared_ptr<Model<>>& m);
    ~GurobiProgram();
    void reset_model();

    bool solve(bool relax = false, double mipgap = 0.01);
    void prepare_model();
    void update_model();
    void update_solution();
    void relax_model();

    void fill_in_grb_vmap();
    void create_grb_constraints();
    void set_grb_objective();

    void print_constraints();
};


#endif /* defined(__Gravity____GurobiProgram) */
