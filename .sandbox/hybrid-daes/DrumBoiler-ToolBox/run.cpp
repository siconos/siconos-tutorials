#include "model_assets.hpp"
#include "SimuEngine.hpp"
#include "ModelIsochore.hpp"
#include "ModelIsobar.hpp"

#define USE_PROGRESS_BAR 0

int run_isobar_evap(){
    std::cout << "\n####### Isobar Evaporation ########\n" << std::endl;
    double P0 = 5000000.0; 
    double h = 0.1;         
    double Tfinal = 3000;

    // /* Initial Conditions */  
    // ######## MIXTURE 1 (vapor)
    double M10 = 0.0;
    double x10 = 0.0;
    double T10 = 300;
    // ######## MIXTURE 2 INITIAL CONDITIONS (liquid)
    double M20 = 1.0;
    double x20 = 0.0; 
    double T20 = 300; 

    ModelNoWallIsobar* model_isobar = new ModelNoWallIsobar();
    ThLPT_perfGaz_incompLiquid* laws = new ThLPT_perfGaz_incompLiquid();
    model_isobar->set_lawsPT(laws);
    model_isobar->initialize_model(P0,M10,x10,T10,M20,x20,T20);
    // model_isobar->display_model_state(1,0);
    model_isobar->set_ExternalHeatExchange(0,0,500,2000,0,0);

    FiniteDifference* diff_method = new FiniteDifference();
    EulerDiscretization* discr = new EulerDiscretization(model_isobar, diff_method, h);
    
    SimuEngine::Options *solver_options = new SimuEngine::Options();
    solver_options->max_iter = 30;
    solver_options->tol = 5e-7;
    solver_options->verbose = 0;
    solver_options->solverID = SICONOS_MCP_NEWTON_FB_FBLSA; 
    solver_options->output_file_name = "results-files/result-isobar-evap.dat";

    SimuEngine simu = SimuEngine(discr, solver_options, Tfinal);

    simu.simulate();
}

int run_isochore_detent(){
    std::cout << "\n####### Isochore Detente ########\n" << std::endl;
    double P0 = 5000000.0; 
    double h = 1;         
    double Tfinal = 4001;

    // /* Initial Conditions */  
    // ######## MIXTURE 1 (vapor)
    double M10 = 1.0;
    double x10 = 1.0;
    double T10 = 834.6821577647171;
    // ######## MIXTURE 2 INITIAL CONDITIONS (liquid)
    double M20 = 0.0;
    double x20 = 1.0; 
    double T20 = 834.6821577647171; 

    ModelNoWallIsochore* model_isochore = new ModelNoWallIsochore();
    ThLPT_perfGaz_incompLiquid* laws = new ThLPT_perfGaz_incompLiquid();
    model_isochore->set_lawsPT(laws);
    model_isochore->initialize_model(P0,M10,x10,T10,M20,x20,T20);
    // model_isochore->display_model_state(1,0);
    model_isochore->set_ExternalHeatExchange(0,0,-1000,0,0,0);

    FiniteDifference* diff_method = new FiniteDifference();
    EulerDiscretization* discr = new EulerDiscretization(model_isochore, diff_method, h);
    
    SimuEngine::Options *solver_options = new SimuEngine::Options();
    solver_options->max_iter = 30;
    solver_options->tol = 5e-7;
    solver_options->verbose = 0;
    solver_options->solverID = SICONOS_MCP_NEWTON_FB_FBLSA; 
    solver_options->output_file_name = "results-files/result-isochore-detent.dat";

    SimuEngine simu = SimuEngine(discr, solver_options, Tfinal);

    simu.simulate();
}

int run_isobar_cond(){
    std::cout << "\n####### Isobar Condensation ########\n" << std::endl;
    double P0 = 470622.1109526682; 
    double h = 0.1;         
    double Tfinal = 861;

    // /* Initial Conditions */  
    // ######## MIXTURE 1 (vapor)
    double M10 = 0.1935257290152718;
    double x10 = 0.9480766560879141;
    double T10 = 422.5958904656754;
    // ######## MIXTURE 2 INITIAL CONDITIONS (liquid)
    double M20 = 0.8064742709847305;
    double x20 = 0.000355928278432825; 
    double T20 = 422.5958904656754; 

    ModelNoWallIsobar* model_isobar = new ModelNoWallIsobar();
    ThLPT_perfGaz_incompLiquid* laws = new ThLPT_perfGaz_incompLiquid();
    model_isobar->set_lawsPT(laws);
    model_isobar->initialize_model(P0,M10,x10,T10,M20,x20,T20);
    // model_isobar->display_model_state(1,0);
    model_isobar->set_ExternalHeatExchange(0,0,-1000,-1000,0,0);

    FiniteDifference* diff_method = new FiniteDifference();
    EulerDiscretization* discr = new EulerDiscretization(model_isobar, diff_method, h);
    
    SimuEngine::Options *solver_options = new SimuEngine::Options();
    solver_options->max_iter = 30;
    solver_options->tol = 5e-7;
    solver_options->verbose = 0;
    solver_options->solverID = SICONOS_MCP_NEWTON_FB_FBLSA; 
    solver_options->output_file_name = "results-files/result-isobar-cond.dat";

    SimuEngine simu = SimuEngine(discr, solver_options, Tfinal);

    simu.simulate();
}

/* 
TODO need update to choose the thermo laws to be as in the report 
(currently need to change rhol() related functions in file laws.cpp  to enable liquid compressibility)
*/
int run_isochore_compress(){
    std::cout << "\n####### Isochore Compression ########\n" << std::endl;
    double P0 = 470622.1109526682; 
    double h = 0.1;         
    double Tfinal = 64.9;

    // /* Initial Conditions */  
    // ######## MIXTURE 1 (vapor)
    double M10 = 0.0003203752961351743;
    double x10 = 1.268783255539628e-13;
    double T10 = 300.0025692754845;
    // ######## MIXTURE 2 INITIAL CONDITIONS (liquid)
    double M20 = 0.9996796179171713;
    double x20 = 0.0; 
    double T20 = 300.0025707874193; 

    ModelNoWallIsochore* model_isochore = new ModelNoWallIsochore();
    ThLPT_CompressLiquid* laws = new ThLPT_CompressLiquid();
    model_isochore->set_lawsPT(laws);
    model_isochore->initialize_model(P0,M10,x10,T10,M20,x20,T20);
    // model_isochore->display_model_state(1,0);
    model_isochore->set_ExternalHeatExchange(0,0,0,1000,0,0);

    FiniteDifference* diff_method = new FiniteDifference();
    EulerDiscretization* discr = new EulerDiscretization(model_isochore, diff_method, h);
    
    SimuEngine::Options *solver_options = new SimuEngine::Options();
    solver_options->max_iter = 30;
    solver_options->tol = 5e-7;
    solver_options->verbose = 0;
    solver_options->solverID = SICONOS_MCP_NEWTON_FB_FBLSA; 
    solver_options->output_file_name = "results-files/result-isochore-compress-incorrect.dat";

    SimuEngine simu = SimuEngine(discr, solver_options, Tfinal);

    simu.simulate();
}


int main(void)
{

    run_isobar_evap();
    // run_isochore_detent();
    // run_isobar_cond();
    run_isochore_compress();


    // double P0 = 5000000.0; // 470622.1109526682;//470567.2629;//5000000.0; 
    // // 10607.40; // 647089.1472783361; // lowpressure end turbine // 5060740.0; // 1060740.0; // 106074.0; // 101325.0;      // TODO need to make the definition of Pext in model.cpp global [further refactoring with c++ object approach instead of c needed]
    // double h = 1;          // s
    // double Tfinal = 4001;//5000;//8619.0; // s (Put Tf-h)

    // // /* Initial Conditions */  
    // // ######## MIXTURE 1 (vapor)
    // double M10 = 1.0;//0.0;//0.19352572; //0.222242570108856;//0.0;
    // //0.2901454432822508;// kg
    // double x10 = 1.0;//0.0;//0.94807665; // 0.8267739; //0.0;
    // //0.8484152486771914; // 1.0; // [0.0 for liquid only | 1.0 for gaz only] 
    // double T10 = 834.6821577647171;//300;//422.59589;// 422.59;// 300;//435.0513307636953; // 300.0; //298.15; // 374.156; //600 //K 
    // // ######## MIXTURE 2 INITIAL CONDITIONS (liquid)
    // double M20 = 0.0;//1.0; //0.80647427; //0.7777574298910757;
    // // 1.0;//0.709854556717751; // kg 
    // double x20 = 1.0; //0.00035592;//0.001381661136610513; // 1.0;  // [0.0 for liquid only | 1.0 for gaz only] 
    // double T20 = 834.6821577647171;//300.0; //422.59589;//422.59;//300;//435.0513307636953; // 300.0; // 298.15; //373.156; //600.0;  //K 

    // // ModelNoWallIsobar* model_isobar = new ModelNoWallIsobar();
    // ModelNoWallIsochore* model_isochore = new ModelNoWallIsochore();
    // model_isochore->initialize_model(P0,M10,x10,T10,M20,x20,T20);
    // model_isochore->display_model_state(1,0);
    // model_isochore->set_ExternalHeatExchange(0,0,-1000,0,0,0);

    // FiniteDifference* diff_method = new FiniteDifference();
    // EulerDiscretization* discr = new EulerDiscretization(model_isochore, diff_method, h);
    
    // SimuEngine::Options *solver_options = new SimuEngine::Options();
    // solver_options->max_iter = 30;
    // solver_options->tol = 5e-7;
    // solver_options->verbose = 1;
    // solver_options->solverID = SICONOS_MCP_NEWTON_FB_FBLSA; 
    // solver_options->output_file_name = "results-files/result-test-object-version.dat";

    // SimuEngine simu = SimuEngine(discr, solver_options, Tfinal);

    // simu.simulate();

}