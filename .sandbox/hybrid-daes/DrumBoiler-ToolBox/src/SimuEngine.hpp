#ifndef SIMUENGINE_HPP
#define SIMUENGINE_HPP

#include "EulerDiscretization.hpp"

#include "SiconosKernel.hpp"
#include "MCP_cst.h"
#include "NonSmoothDrivers.h"
#include "MixedComplementarityProblem.h"
#include "SolverOptions.h"
#include "MCP_Solvers.h"
#include "NumericsVerbose.h"


class SimuEngine{

    public:

        typedef struct
        {
            int max_iter;
            double tol;
            int verbose; // 0 -- 1 -- 2 -- 3 (increasing from none to max)
            MCP_SOLVER solverID;
            std::string output_file_name;
        } Options;

        SimuEngine(EulerDiscretization* discr_, Options* solver_options_, double Tf)
        {
            discr = discr_;
            solver_options = solver_options_;
            Tfinal = Tf;

        };
        ~SimuEngine(){
            delete(discr);
            delete(solver_options);
        };

        MixedComplementarityProblem * create_mcp(EulerDiscretization::NumericalSimuInfo* info, ptrFunctionMCP2 F, ptrFunctionMCP_nabla NablaF);

        void solve_mcp_newton(MixedComplementarityProblem* problem, EulerDiscretization::NumericalSimuInfo* info, double* z_sol);

        int simulate();

    private:
        EulerDiscretization* discr;
        Options* solver_options;
        double Tfinal;
};

#endif // SIMULATOR_HPP