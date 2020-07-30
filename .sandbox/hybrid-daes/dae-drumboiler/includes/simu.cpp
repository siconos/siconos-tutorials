#include "simu.h"

// (self note) static makes function callable only in this file (can still be transferred through pointer)
static MixedComplementarityProblem * create_mcp(NumericalSimuInfo* info, ptrFunctionMCP2 F, ptrFunctionMCP_nabla NablaF)
{
    /* Create a MixedComplementarityProblem */
    MixedComplementarityProblem* problem = (MixedComplementarityProblem *)malloc(sizeof(MixedComplementarityProblem));

    int n = info->nb_eqs;
    
    problem->n1 = n-info->nb_compl_eqs; //nb equality constraint
    problem->n2 = info->nb_compl_eqs; // nb complementarity constraint
    problem->compute_Fmcp = F ;
    problem->compute_nabla_Fmcp = NablaF ;
    problem->nabla_Fmcp =  NM_create(NM_DENSE, n, n); 
    problem->env = (void*)info;
 
    return problem;
}


// solve a MCP problem in format MixedComplementarityProblem using siconos newton algorithm
static void solve_mcp_newton(MixedComplementarityProblem* problem, NumericalSimuInfo* info, Options* sopt, double* z_sol)
{
  printf("test_mcp_newton() starts for solver %s.\n", solver_options_id_to_name(sopt->solverID));
  
  // AT this point z_sol = info->z_prev  
  
  int return_info = 1 ;
  
  /* Set solver options */
  SolverOptions * options = solver_options_create(sopt->solverID);
  options->iparam[SICONOS_IPARAM_MAX_ITER] = sopt->max_iter;
  options->dparam[SICONOS_DPARAM_TOL] = sopt->tol;
  
  /* VERBOSE MODE */
  numerics_set_verbose(sopt->verbose);

  int size = problem->n1 + problem->n2 ;
  
  double * initial_guess = (double *)calloc(size, sizeof(double));

  problem->compute_Fmcp( (void*)info, size, z_sol, initial_guess ); // initial guess of w = F(z_curr)

  // display constraints values at current (soon previous) point
  std::cout << "F(z_prev) = "  << std::endl;
  for(unsigned int i=0; i<size; i++)
  {
    std::cout << initial_guess[i] << " " ;
  }
  std::cout << std::endl;

  std::cout << "initial err on liquid mass balance eq. , F[8] = "<< initial_guess[8] << std::endl;
  std::cout << "initial err on liquid energy balance  eq. , F[9] = "<< initial_guess[9] << std::endl;
  std::cout << "[solve_mcp_newton] -- z_prev = "  << std::endl;
  for(unsigned int i=0;i<size;i++)
  {
    std::cout << info->z_prev[i]<< " " ;
  }
  std::cout << std::endl;

  // computing next step -> stored in z 
  // info = mcp_driver(problem, z , w,  &options); // OLD CODE
  return_info = mcp_driver(problem, z_sol , initial_guess,  options); 

  // AT this point z_sol is now updated


  free(initial_guess); // free initial_guess memory (could be avoided if the pointer is made "global")

  //  double * nablaFz = problem->nabla_Fmcp->matrix0;
  //  double cond_num = cond(nablaFz, size, size);
  //  cout << "########## nablaF(z) conditionning: " << cond_num  << endl;
  //   NM_display(problem->nabla_Fmcp);

  /* Freeing the MCP*/
  // delete(problem->nabla_Fmcp); // matrix0 is always re-written, is it enough ?
  // free(problem); // one problem
}


int simulate(NumericalSimuInfo* simu_info, Options* solver_options){

    int N = ceil(simu_info->Tfinal/simu_info->h)+1;
    int n = simu_info->nb_eqs;
    unsigned int outputSize = n+1;
    SimpleMatrix dataPlot(N + 1, outputSize);
    dataPlot(0, 0) = 0.0;
    for(unsigned int i=0;i<n;i++)
    {
      dataPlot(0, i+1) = simu_info->z_prev[i] ;
    }

    // Solution vector
    double* z_sol = (double *)calloc(n,sizeof(double));
    // z_sol = z_prev =  initial condition
    memcpy(z_sol, simu_info->z_prev, n * sizeof(double));

    // Discrete problem to solve
    MixedComplementarityProblem* discrete_pb = create_mcp(simu_info, &implicitEulerDiscr, &implicitEulerDiscrJac);

    // --- Time loop ---
    std::cout << "====> Start computation ... " << std::endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    for(unsigned int j=0;j<N;j++)
    {
        // Solve MCP discrete problem    
        solve_mcp_newton(discrete_pb, simu_info, solver_options, z_sol);         
        
        // Previous-step = Current-step
        memcpy(simu_info->z_prev, z_sol, n * sizeof(double)); 

        display_model_state(z_sol, n, simu_info->h, k);
         
        // store result in table
        dataPlot(k, 0) = (k*simu_info->h);
        for(unsigned int i=0;i<n;i++)
        {
            dataPlot(k, i+1) = z_sol[i] ;
        }
      
      k++;
    }
    std::cout << "End of computation - Number of iterations done: " << k - 1 << std::endl;
    // --- Output files ---
    std::cout << "====> Output file writing ..." << std::endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write(solver_options->output_file_name, "ascii", dataPlot, "noDim");

    return 0;
}

