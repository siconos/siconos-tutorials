#ifndef SIMU_H
#define SIMU_H

#include "model_common.h"
#include "discrete_system.h"

#include "SiconosKernel.hpp"
#include "MCP_cst.h"
#include "NonSmoothDrivers.h"
#include "MixedComplementarityProblem.h"
#include "SolverOptions.h"
#include "MCP_Solvers.h"
#include "NumericsVerbose.h"


typedef struct
{
    int max_iter;
    double tol;
    int verbose; // 0 -- 1 -- 2 -- 3 (increasing from none to max)
    MCP_SOLVER solverID;
    std::string output_file_name;
}Options;

static MixedComplementarityProblem * create_mcp(NumericalSimuInfo* info);

static void solve_mcp_newton(MixedComplementarityProblem* problem, NumericalSimuInfo* info, Options* sopt, double* z_sol);

int simulate(NumericalSimuInfo* simulation_info, Options* solver_options);

#endif // SIMU_H