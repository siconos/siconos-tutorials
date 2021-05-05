/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*!\file absolute-value-and-sign-lcp.cpp
 
  Example of Linear switching DAE with a constraint defined by a multivalued operator
*/

#include "SiconosKernel.hpp"
#include "lcp_cst.h"

using namespace std;

int main(int argc, char* argv[])
{
    try
    {

        // ================= Creation of the model =======================

        // User-defined main parameters
        unsigned int dimX       = 3;    // Dimension of the system state variables
        unsigned int dimLambda  = 3;    // Dimension of the system lambda variables

        double t0       = 0.0;          // initial computation time
        double T        = 10.0;         // final computation time 
      
        double h        =  0.05;        // time step

        double x1_0     = -1.0;         // initial condition in state variable x1
        double x2_0     = 0.0;          // initial condition in state variable x2
        double z_0      = 0;            // initial condition in algebraic variable z


        // -------------------------
        // --- Dynamical systems ---
        // -------------------------
        cout << "###### DAE with absolute value and sign in constraint #####" <<  endl;  
        cout << "====> Model definition ..." <<  endl;


        SP::SiconosVector init(new SiconosVector({x1_0, x2_0, z_0}));

        SP::SiconosMatrix A( new SimpleMatrix(dimX,dimX) ); 
        // *** sliding
        // double B0 = 0.0;
        // double B1 = 1.0;
        // should have slide or jump (multiple solutions)
        double B0 = -1.0;
        double B1 = 0.5;
        A->setRow(0,SiconosVector({0.0, 0.0, B0}));
        A->setRow(1,SiconosVector({0.0, 0.0, B1}));
        A->setRow(2,SiconosVector({1.0, -1.0, 0.0}));

        cout << "matrix A: " << endl;
        A->display();

        SP::SimpleMatrix E(new SimpleMatrix(dimX,dimX));
        (*E)(0,0) = 1.0;
        (*E)(1,1) = 1.0;

        cout << "matrix E: " << endl;
        E->display();

        SP::SiconosVector b(new SiconosVector({1.0, 0.0, 1.0}));

        cout << "vector b: " << endl;
        b->display();

        // Siconos smooth dynamical system
        SP::FirstOrderLinearTIDS dyn(new FirstOrderLinearTIDS(init,A));
        dyn->setbPtr(b);
        dyn->setMPtr(E);

        // -------------------------
        // --- LCP Relation ---
        // -------------------------


        // IDE complains when smartpoint is SP::SiconosMatrix why ? not virtual ?
        SP::SimpleMatrix C( new SimpleMatrix(dimLambda,dimX) );
        (*C)(0,0) = 2.0;
        (*C)(1,0) = 1.0;

        SP::SimpleMatrix D( new SimpleMatrix(dimLambda,dimLambda) );
        (*D)(0,0) = 1.0;
        (*D)(1,2) = 1.0;
        (*D)(2,1) = -1.0;

        SP::SimpleMatrix R( new SimpleMatrix(dimX,dimLambda) );
        (*R)(2,0) = 1.0;
        (*R)(2,1) = -1.0;

        SP::SiconosVector e(new SiconosVector({0.0, 0.0, 2.0}));

        // Relation LCP lhs
        SP::FirstOrderLinearTIR relation(new FirstOrderLinearTIR(C, R) );
        relation->setDPtr(D);
        relation->setePtr(e);

        // NonSmooth law: LCP
        SP::NonSmoothLaw nslaw(new ComplementarityConditionNSL(dimLambda));

        // interaction 
        SP::Interaction inter(new Interaction(nslaw, relation));

        // -----------------------------
        // --- Siconos Model Entity ---
        // ----------------------------
        SP::NonSmoothDynamicalSystem switch_dae(new NonSmoothDynamicalSystem(t0, T));

        // add the dynamical system in the non smooth dynamical system
        switch_dae->insertDynamicalSystem(dyn);

        // link the interaction and the dynamical system
        switch_dae->link(inter, dyn);

        // -----------------------------
        // --- Simulation Definition ---
        // -----------------------------

        // -- (1) OneStepIntegrators --
        double theta = 1.0;
        double gamma = 1.0;
        SP::EulerMoreauOSI osi(new EulerMoreauOSI(theta,gamma));


        // -- (2) Time discretisation --
        SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));

        // -- (3) one step non smooth problem
        SP::LCP osnspb(new LCP(SICONOS_LCP_ENUM));
        // SP::LCP osnspb(new LCP());
        osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS] = 1; // 1 for multiple solutions,
                                                                                         // 0 else .
        osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_USE_DGELS] = 0; // 1 if useDGELS 
        osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_SKIP_TRIVIAL] = SICONOS_LCP_SKIP_TRIVIAL_YES;
        osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_SEED] = 4; // SEED
        osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_MIN_NORM_SOLUTION] = 1;
        osnspb->setNumericsVerboseMode(true);
        
        
        // -- (4) Simulation setup with (1) (2) (3)
        SP::TimeStepping s(new TimeStepping(switch_dae, td, osi, osnspb));

        // s->setResetAllLambda(false);

        // =========================== End of model definition ===========================
        cout <<  "====> ... End of Model definition" <<  endl;
        //   // ================================= Computation =================================

        int N = ceil((T - t0) / h)+1; // Number of time steps

        // --- Get the values to be plotted ---
        // -> saved in a matrix dataPlot
        unsigned int outputSize = 7;
        SimpleMatrix dataPlot(N + 1, outputSize);

        SP::SiconosVector x = dyn->x();
        SP::SiconosVector lambda = inter->lambda(0);

        dataPlot(0, 0) = switch_dae->t0();
        dataPlot(0, 1) = (*x)(0);   // x1
        dataPlot(0, 2) = (*x)(1);   // x2
        dataPlot(0, 3) = (*x)(2);   // z
        dataPlot(0, 4) = (*lambda)(0);
        dataPlot(0, 5) = (*lambda)(1);       // 
        dataPlot(0, 6) = (*lambda)(2);       // 
        // --- Time loop ---
        cout << "====> Start computation ... " << endl;

        cout << "sigma = h*z = " << h*dataPlot(0, 3) << endl; // value \sigma_z
        cout << "x1 = " << dataPlot(0, 1) << endl; // position in x1

        // ==== Simulation loop - Writing without explicit event handling =====
        int k = 1;
        boost::progress_display show_progress(N);
        boost::timer time;
        time.restart();

        // // element M[0,0] of the matrix M (normaly of size 1) s.t. y = M*lambda+q
        // double M_00;
        //  // element q[0] of the vector q (normaly of size 1) s.t. y = M*lambda+q
        // double q_0;

        while (s->hasNextEvent())
        {
            
            s->computeOneStep();
            osnspb->display();

            cout << "# of solutions: " 
                 << osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_NUMBER_OF_SOLUTIONS] << endl; // Number of solutions
            cout << "Mode of solutions: " 
                 << osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_CURRENT_ENUM] << endl; // Number of solutions
            
            // M_00 = osnspb->M()->defaultMatrix()->getValue(0,0);
            // q_0  = osnspb->q()->getValue(0);
            // --- Get values to be plotted ---
            dataPlot(k, 0) =  s->nextTime();
            dataPlot(k, 1) = (*x)(0);
            dataPlot(k, 2) = (*x)(1);
            dataPlot(k, 3) = (*x)(2);
            dataPlot(k, 4) = (*lambda)(0);
            dataPlot(k, 5) = (*lambda)(1);
            dataPlot(k, 6) = (*lambda)(2);
            cout << "sigma = h*z = " << h*dataPlot(k, 3) << endl; // value \sigma_z
            cout << "x1 = " << dataPlot(k, 1) << endl; // position in x1
            SiconosVector temp_res = SiconosVector(outputSize);
            dataPlot.getRow(k,temp_res);
            cout << "all results: " << temp_res << endl;
            s->nextStep();
            ++show_progress;
            k++;

        }
        cout  << "End of computation - Number of iterations done: " << k - 1 << endl;
        cout << "Computation Time " << time.elapsed()  << endl;

        // --- Output files ---
        cout << "====> Output file writing ..." << endl;
        dataPlot.resize(k, outputSize);
        ioMatrix::write("result_CONTINUE-big-h_absolute-value-and-sign-lcp.dat", "ascii", dataPlot, "noDim");
    }

    catch (SiconosException& e)    
    {
    cerr << e.report() << endl;
    return 1;
    }

    catch (...) 
    {
    cerr << "Exception caught in absolute-value-and-sign-lcp.cpp" << endl;
    return 1;
    }


}