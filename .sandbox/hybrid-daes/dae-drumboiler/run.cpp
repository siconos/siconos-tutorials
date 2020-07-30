#include "model_common.hpp"
#include "model.hpp"
#include "simu.hpp"

int main(void)
{
    double P0 = 1060740.0; // 106074.0; // 101325.0;      // TODO need to make the definition of Pext in model.cpp global [further refactoring with c++ object approach instead of c needed]
    double h = 1.0;          // s
    double Tfinal = 2600.0; // s

    // /* Initial Conditions */  
    // ######## MIXTURE 1 (vapor)
    double M10 = 1.0;// kg
    double x10 = 1.0; // Kg
    double T10 = 600; //298.15; // 374.156; //600 //K 
    // ######## MIXTURE 2 INITIAL CONDITIONS (liquid)
    double M20 = 0.0; // kg 
    double x20 = 1.0; // 1.0;  // Kg
    double T20 = 600; // 298.15; //373.156; //600.0;  //K 

    NumericalSimuInfo simu_info;
    simu_info.nb_eqs = 21;
    simu_info.h = h;
    simu_info.Tfinal = Tfinal;

    initialize_simu(&simu_info,P0,M10,x10,T10,M20,x20,T20); // hand writing

    Options solver_options;
    solver_options.max_iter = 50;
    solver_options.tol = 7e-7;
    solver_options.verbose = 1;
    solver_options.solverID = SICONOS_MCP_NEWTON_FB_FBLSA; 
    solver_options.output_file_name = "result-double-mix-mcp-Cond-step-1s-FB-50iterations-Tf=2600.dat";
    
    simulate(&simu_info, &solver_options);
}