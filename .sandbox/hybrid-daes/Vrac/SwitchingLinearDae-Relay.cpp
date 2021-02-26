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

/*!\file clutchFrictionEngaging.cpp
 
  Engaging of a clutch modeled using a friction coefficient.
  Simulation with a Time-Stepping scheme.
*/

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    try
    {

        // ================= Creation of the model =======================

        // User-defined main parameters
        unsigned int dimX       = 3;        // Dimension of the system state variables
        unsigned int dimLambda  = 1;        // Dimension of the system lambda variables

        double t0   = 0.0;                  // initial computation time
        double T    = 1.53;                  // final computation time 

        double h =  0.0005;                   // time step

        double x1_0 = -1.0;                 // initial condition in state variable x1
        double x2_0 = 0.1;                  // initial condition in state variable x2
        double z_0  = 0;                    // initial condition in algebraic variable z


        // -------------------------
        // --- Dynamical systems ---
        // -------------------------

        cout << "====> Model definition ..." <<  endl;


        SP::SiconosVector init(new SiconosVector({x1_0, x2_0, z_0}));

        SP::SiconosMatrix A( new SimpleMatrix(dimX,dimX) ); 
        double B0 = 0.1;
        double B1 = 0.9;
        A->setRow(0,SiconosVector({0.0, 0.0, B0}));
        A->setRow(1,SiconosVector({0.0, 0.0, B1}));
        A->setRow(2,SiconosVector({0.0, -1.0, 0.0}));

        cout << "matrix A: " << endl;
        A->display();

        SP::SimpleMatrix E(new SimpleMatrix(dimX,dimX));
        (*E)(0,0) = 1.0;
        (*E)(1,1) = 1.0;

        cout << "matrix E: " << endl;
        E->display();

        SP::SiconosVector b(new SiconosVector({1.0, 0.0, 0.0}));

        cout << "vector b: " << endl;
        b->display();

        // Siconos smooth dynamical system
        SP::FirstOrderLinearTIDS dyn(new FirstOrderLinearTIDS(init,A));
        dyn->setbPtr(b);
        dyn->setMPtr(E);

        // -------------------------
        // --- Relay Relation ---
        // ---  -y \in N_{[bl,bu]}(lambda) ---
        // ---  dynamic input: B.dot(lambda)
        // -------------------------

        // Relation LCP lhs
        SP::FirstOrderR relation(new FirstOrderNonLinearR() );

        //Plugin for output function h(x,lambda)
        relation->setComputehFunction("NLRelation_Plugin","Computeh");
        //Plugin for output function jacobian wrt x
        relation->setComputeJachxFunction ("NLRelation_Plugin","ComputeJxh");
        //Plugin for output function jacobian wrt lambda
        relation->setComputeJachlambdaFunction("NLRelation_Plugin","ComputeJlh");
        //Plugin for input function g(x,lambda)
        relation->setComputegFunction("NLRelation_Plugin","Computeg");
        //Plugin for input function jacobian wrt x
        relation->setComputeJacgxFunction ("NLRelation_Plugin","ComputeJxg");
        //Plugin for input function jacobian wrt lambda
        relation->setComputeJacglambdaFunction("NLRelation_Plugin","ComputeJlg");

        // NonSmooth law: relay on [-1,1] rhs
        SP::NonSmoothLaw nslaw(new RelayNSL(dimLambda, -1.0, 1.0));

        // interaction -y in relay[-1,1]
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
        SP::Relay osnspb(new Relay());

        // osnspb->setNumericsVerboseMode(true);
        

        // -- (4) Simulation setup with (1) (2) (3)
        SP::TimeStepping s(new TimeStepping(switch_dae, td, osi, osnspb));
        
        
        // Non linear flags set to true
        // Newton iteration and accuracy increased ...
        s->setComputeResiduY(true);
        s->setComputeResiduR(true);
        s->setNewtonMaxIteration(100);
        s->setNewtonTolerance(1e-10);
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
        dataPlot(0, 1) = (*x)(0);
        dataPlot(0, 2) = (*x)(1);
        dataPlot(0, 3) = (*x)(2);
        dataPlot(0, 4) = (*lambda)(0);
        dataPlot(0, 5) = 0.0;   // M[0,0] at t= 0 initialized at 0.0
        dataPlot(0, 6) = 0.0;   // q[0] at t= 0 initialized at 0.0
        // --- Time loop ---
        cout << "====> Start computation ... " << endl;
        // ==== Simulation loop - Writing without explicit event handling =====
        int k = 1;
        boost::progress_display show_progress(N);
        boost::timer time;
        time.restart();

        // element M[0,0] of the matrix M (normaly of size 1) s.t. y = M*lambda+q
        double M_00;
         // element q[0] of the vector q (normaly of size 1) s.t. y = M*lambda+q
        double q_0;
        while (s->hasNextEvent())
        {
            
            s->computeOneStep();
            // osnspb->display();
            M_00 = osnspb->M()->defaultMatrix()->getValue(0,0);
            q_0  = osnspb->q()->getValue(0);
            // --- Get values to be plotted ---
            dataPlot(k, 0) =  s->nextTime();
            dataPlot(k, 1) = (*x)(0);
            dataPlot(k, 2) = (*x)(1);
            dataPlot(k, 3) = (*x)(2);
            dataPlot(k, 4) = (*lambda)(0);
            dataPlot(k, 5) = M_00;
            dataPlot(k, 6) = q_0;
            s->nextStep();
            ++show_progress;
            k++;

        }
        cout  << "End of computation - Number of iterations done: " << k - 1 << endl;
        cout << "Computation Time " << time.elapsed()  << endl;

        // --- Output files ---
        cout << "====> Output file writing ..." << endl;
        dataPlot.resize(k, outputSize);
        ioMatrix::write("result_SwitchingLinearDae.dat", "ascii", dataPlot, "noDim");
    }

    catch (SiconosException& e)    
    {
    cerr << e.report() << endl;
    return 1;
    }

    catch (...) 
    {
    cerr << "Exception caught in SwitchingLinearDae.cpp" << endl;
    return 1;
    }


}