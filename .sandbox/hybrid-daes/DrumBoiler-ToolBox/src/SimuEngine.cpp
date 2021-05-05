#include "SimuEngine.hpp"

// (self note) static makes function callable only in this file (can still be transferred through pointer)
MixedComplementarityProblem * SimuEngine::create_mcp(EulerDiscretization::NumericalSimuInfo* info, ptrFunctionMCP2 F, ptrFunctionMCP_nabla NablaF)
{
    // std::cout << "\n Defining the MCP corresponding to the Discretization";
    //  std::cout << "\n In [simu.create_mcp()] : ";
    /* Create a MixedComplementarityProblem */
    MixedComplementarityProblem* problem = (MixedComplementarityProblem *)malloc(sizeof(MixedComplementarityProblem));
    Model* model = discr->get_model();
    int n = model->get_nb_eqs();

    // std::cout << "\n Problem Size : " << n ;

    problem->n1 = model->get_nb_equaeqs(); //nb equality constraint
    problem->n2 = model->get_nb_compleqs(); // nb complementarity constraint

    // std::cout << "\n Equality Pb Size : " << problem->n1 ;
    // std::cout << "\n Complementarity Pb Size : " << problem->n2 << std::endl;

    problem->compute_Fmcp = F ;
    problem->compute_nabla_Fmcp = NablaF ;
    problem->nabla_Fmcp =  NM_create(NM_DENSE, n, n); 
    problem->env = (void*)info;

    //NumericalSimuInfo* test = (NumericalSimuInfo*)problem->env;
    // std::cout << "\n TIME STEP : " <<  test->h << std::endl;

 
    return problem;
}


// solve a MCP problem in format MixedComplementarityProblem using siconos newton algorithm
void SimuEngine::solve_mcp_newton(MixedComplementarityProblem* problem, EulerDiscretization::NumericalSimuInfo* info, double* z_sol)
{
  // printf("mcp_newton() starts for solver %s.\n", solver_options_id_to_name(solver_options->solverID));
  
  // AT this point z_sol = info->z_prev  
  int return_info = 1 ;
  
  /* Set solver options */
  SolverOptions * options = solver_options_create(solver_options->solverID);
  options->iparam[SICONOS_IPARAM_MAX_ITER] = solver_options->max_iter;
  options->dparam[SICONOS_DPARAM_TOL] = solver_options->tol;
  
  /* VERBOSE MODE */
  numerics_set_verbose(solver_options->verbose);

  // Initial inputs
  int size = problem->n1 + problem->n2 ;

  // std::cout << "\nInitial input, z_sol = [";
  //   for(int i=0;i<size;i++){
  //     std::cout << z_sol[i] << ", ";
  //   }
  // std::cout << " ]" << std::endl;

  // std::cout << "size = " << size << std::endl;
  double * initial_guess = (double *)calloc(size, sizeof(double));

  // printf(" [solve_mcp_newton] simu_info = %p\n", info);
  // std::cout << "\n [solve_mcp_newton] Compute initial guess compute_Fmcp(z_sol)" << std::endl;
  problem->compute_Fmcp( (void*)info, size, z_sol, initial_guess ); // initial guess of w = F(z_curr)

  // display constraints values at current (soon previous) point
  // std::cout << "F(z_prev) = "  << std::endl;
  // for(unsigned int i=0; i<size; i++)
  // {
  //   std::cout << initial_guess[i] << " " ;
  // }
  // std::cout << std::endl;

  // std::cout << "initial err on liquid mass balance eq. , F[2] = "<< initial_guess[2] << std::endl;
  // std::cout << "initial err on liquid energy balance  eq. , F[3] = "<< initial_guess[3] << std::endl;
  // std::cout << "[solve_mcp_newton]  "  << std::endl;
  
  // computing next step -> stored in z 
  // info = mcp_driver(problem, z , w,  &options); // OLD CODE
  return_info = mcp_driver(problem, z_sol , initial_guess,  options); 

  // printf("[solve_mcp_newton] -- Done : err = %d.\n", return_info);
  
  // AT this point z_sol is now updated
  // std::cout << "\n Computed solution, z_sol = [";
  //   for(int i=0;i<size;i++){
  //     std::cout << z_sol[i] << ", ";
  //   }
  // std::cout << " ]" << std::endl;


  free(initial_guess); // free initial_guess memory (could be avoided if the pointer is made "global")

  //  double * nablaFz = problem->nabla_Fmcp->matrix0;
  //  double cond_num = cond(nablaFz, size, size);
  //  cout << "########## nablaF(z) conditionning: " << cond_num  << endl;
  //   NM_display(problem->nabla_Fmcp);

  /* Freeing the MCP*/
  // delete(problem->nabla_Fmcp); // matrix0 is always re-written, is it enough ?
  // free(problem); // one problem
}

//Options* solver_options
//EulerDiscretization* euler
int SimuEngine::simulate(){
    Model* model = discr->get_model();

    EulerDiscretization::NumericalSimuInfo* simu_info = new EulerDiscretization::NumericalSimuInfo();
    simu_info->discretization = discr;

    int N = ceil(Tfinal/discr->h)+1;
    int n = model->get_nb_eqs();
    unsigned int outputSize = n+1;
    SimpleMatrix dataPlot(N + 1, outputSize);
    dataPlot(0, 0) = 0.0;
    
    
    double* z_prev = model->get_z();

    for(unsigned int i=0;i<n;i++)
    {
      dataPlot(0, i+1) = z_prev[i] ;
    }

    // Solution vector
    double* z_sol = (double *)calloc(n,sizeof(double));

    
    // printf(" [simulate] (before update) h = %f\n", discr->h);
    // z_sol = z_prev =  initial condition
    memcpy(z_sol, z_prev, n * sizeof(double));
    // printf(" [simulate](before update) h = %f\n", simu_info->h);

    // Discrete problem to solve
    // printf(" [simulate] creation mcp \n");
    MixedComplementarityProblem* discrete_pb = create_mcp(simu_info, &EulerDiscretization::Discr, &EulerDiscretization::DiscrJac);
    
    // --- Time loop ---
    std::cout << "====> Start computation ... " << std::endl;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    for(unsigned int j=0;j<N;j++)
    {
        // Solve MCP discrete problem    
        solve_mcp_newton(discrete_pb, simu_info, z_sol);  
        // Previous-step = Current-step
        memcpy(z_prev, z_sol, n * sizeof(double)); 

        if(solver_options->verbose>0)
          model->display_model_state(discr->h, k);
         
        // store result in table
        dataPlot(k, 0) = (k*discr->h);
        for(unsigned int i=0;i<n;i++)
        {
            dataPlot(k, i+1) = z_sol[i] ;
        }
      //  if(z_sol[P]>=5000000.0)
      //    break;
      //  if(z_sol[T1]<=300.0)
      //    break;
      //  if(z_sol[P]<=470567.0)
      //    break;
      k++;
      #if USE_PROGRESS_BAR
        printProgress(((float)j/(N-1)));
      #endif  
    }
    std::cout << std::endl << "End of computation - Number of iterations done: " << k - 1 << std::endl;
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    std::cout << "Computation time : " << elapsed << " ms" << std::endl;

    model->display_model_state(discr->h, k-1);

    // --- Output files ---
    std::cout << "====> Output file writing ..." << std::endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write(solver_options->output_file_name, "ascii", dataPlot, "noDim");

    return 0;
}

