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
        unsigned int dimX = 3;            // Dimension of the system state variables
        unsigned int dimLambda = 3;       // Dimension of the system lambda variables

        double t0 = 0;                    // initial computation time
        double tswitch = 5;               // switching time
        double T = 10;                    // final computation time
        double h = 0.05;                   // time step

        double omega1_init = 10;           // initial rotation speed of disk 1
        double omega2_init = 5;            // initial rotation speed of disk 2

        double J1 = 1.0;                    // inertial of disk 1
        double J2 = 0.25;                   // inertial of disk 2

        double T1 = 10;                     // Torque on disk 1 when clutch is not engaged
        double T2 = -10;                    // Torque on disk 2 when clutch is not engaged

        double alpha = 3000;                // Friction coefficient

        // -------------------------
        // --- Dynamical systems ---
        // -------------------------

        cout << "====> Model definition ..." <<  endl;


        SP::SiconosVector init(new SiconosVector({tswitch, omega1_init, omega2_init}));

        SP::SiconosMatrix A(new SimpleMatrix( ZeroMat(dimX,dimX) )); 
        SP::SiconosVector b(new SiconosVector({-1.0,T1/J1,T2/J2}));

        // Siconos smooth dynamical system
        SP::FirstOrderLinearTIDS dyn(new FirstOrderLinearTIDS(init,A,b));

        // -------------------------
        // --- Complemetary Relation ---
        // ---  y _|_ lambda ---
        // ---  dynamic input: B.dot(lambda)
        // -------------------------

        // Jacobian of y wrt x
        SP::SimpleMatrix C(new SimpleMatrix(dimLambda,dimX));
        (*C)(0,0) = 1.0;
        (*C)(1,1) = 1.0;
        (*C)(1,2) = -1.0;

        // Jacobian of y wrt lambda
        SP::SimpleMatrix D(new SimpleMatrix(dimLambda,dimLambda));
        (*D)(1,2) = 1.0;
        (*D)(2,0) = 2.0;
        (*D)(2,1) = -1.0;

        // Jacobian of r wrt lambda
        SP::SimpleMatrix B(new SimpleMatrix(dimX,dimLambda));
        (*B)(0,0) = 1.0/alpha;
        B->setRow(1,SiconosVector({ -1.0/J1, 1.0/J1, 0.0}));
        B->setRow(2,SiconosVector({ 1.0/J2, -1.0/J2, 0.0}));

        // Relation LCP lhs
        SP::FirstOrderLinearTIR relation(new FirstOrderLinearTIR(C,B));
        relation->setDPtr(D);

        // NonSmooth law LCP rhs
        SP::NonSmoothLaw nslaw(new ComplementarityConditionNSL(dimLambda));

        // interaction: complete LCP
        SP::Interaction inter(new Interaction(nslaw, relation));

        // -----------------------------
        // --- Siconos Model Entity ---
        // ----------------------------
        SP::NonSmoothDynamicalSystem clutch(new NonSmoothDynamicalSystem(t0, T));

        // add the dynamical system in the non smooth dynamical system
        clutch->insertDynamicalSystem(dyn);

        // link the interaction and the dynamical system
        clutch->link(inter, dyn);

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
        SP::OneStepNSProblem osnspb(new LCP());

        // -- (4) Simulation setup with (1) (2) (3)
        SP::TimeStepping s(new TimeStepping(clutch, td, osi, osnspb));

        // =========================== End of model definition ===========================
        cout <<  "====> ... End of Model definition" <<  endl;
        //   // ================================= Computation =================================

        int N = ceil((T - t0) / h)+10000; // Number of time steps

        // --- Get the values to be plotted ---
        // -> saved in a matrix dataPlot
        unsigned int outputSize = 7;
        SimpleMatrix dataPlot(N + 1, outputSize);

        SP::SiconosVector x = dyn->x();
        SP::SiconosVector lambda = inter->lambda(0);

        dataPlot(0, 0) = clutch->t0();
        dataPlot(0, 1) = (*x)(0);
        dataPlot(0, 2) = (*x)(1);
        dataPlot(0, 3) = (*x)(2);
        dataPlot(0, 4) = (*lambda)(0);
        dataPlot(0, 5) = (*lambda)(1);
        dataPlot(0, 6) = (*lambda)(2);
        // --- Time loop ---
        cout << "====> Start computation ... " << endl;
        // ==== Simulation loop - Writing without explicit event handling =====
        int k = 1;
        boost::progress_display show_progress(N);
        boost::timer time;
        time.restart();

        while (s->hasNextEvent())
        {
            s->computeOneStep();
            // --- Get values to be plotted ---
            dataPlot(k, 0) =  s->nextTime();
            dataPlot(k, 1) = (*x)(0);
            dataPlot(k, 2) = (*x)(1);
            dataPlot(k, 3) = (*x)(2);
            dataPlot(k, 4) = (*lambda)(0);
            dataPlot(k, 5) = (*lambda)(1);
            dataPlot(k, 6) = (*lambda)(2);
            s->nextStep();
            ++show_progress;
            k++;

        }
        cout  << "End of computation - Number of iterations done: " << k - 1 << endl;
        cout << "Computation Time " << time.elapsed()  << endl;

        // --- Output files ---
        cout << "====> Output file writing ..." << endl;
        dataPlot.resize(k, outputSize);
        ioMatrix::write("result_ClutchEngaging.dat", "ascii", dataPlot, "noDim");
        //double error=0.0, eps=1e-12;
        // if ((error=ioMatrix::compareRefFile(dataPlot, "ClutchEngaging.ref", eps)) >= 0.0
        //     && error > eps)
        //     return 1;
        // }
    }

    catch (SiconosException& e)    
    {
    cerr << e.report() << endl;
    return 1;
    }

    catch (...) 
    {
    cerr << "Exception caught in ClutchFrictionEngaging.cpp" << endl;
    return 1;
    }


}