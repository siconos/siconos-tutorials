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

/*!\file ImpactingStick.cpp
*/

#include "FrictionContact.hpp"
#include "TimeStepping.hpp"
#include <SiconosKernel.hpp>
#include <chrono>
#include <cmath>



using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 0.2;                  // final computation time
    double h = 1e-4;                // time step
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    
    double l = 1.0; // length of the stick
    double m = 1.;  // mass
    double g = 10.; // Gravity

    double pi = acos(-1);
    double initial_angle=pi/4.0;

    
    
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." <<  endl;

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = m*l*l/12.;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0)(0) = 0.5 * l * sin(initial_angle);
    (*q0)(1) = 0.5 * l * cos(initial_angle)+0.01;
    (*q0)(2) = initial_angle;

    
    (*v0)(0) = -0.5;
    (*v0)(1) = 0.1;
    (*v0)(2) = 0.1;
    
    // -- The dynamical system --
    SP::LagrangianDS stick(new LagrangianDS(q0, v0, Mass));


    
    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(1) = -m * g;
    stick->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 1.0;
    double mu = 0.01;
    // Interaction stick-floor
    //

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(e,0.0,mu,2));

    SP::Relation relation(new LagrangianScleronomousR("ImpactingStickPlugin:h1", "ImpactingStickPlugin:G1"));

    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem bouncingStick(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    bouncingStick->insertDynamicalSystem(stick);

    // link the interaction and the dynamical system
    bouncingStick->link(inter, stick);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));
 

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new FrictionContact(2));

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(bouncingStick, t, OSI, osnspb));

    s->setNewtonOptions(SICONOS_TS_LINEAR);
    //OSI->setGamma(3/2.0);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================


    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 15;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = stick->q();
    SP::SiconosVector v = stick->velocity();
    SP::SiconosVector p = stick->p(1);
    SP::SiconosVector lambda = inter->lambda(1);
    SP::SiconosVector y = inter->y(0);
    SP::SiconosVector u = inter->y(1);

    dataPlot(0, 0) = bouncingStick->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*v)(1);
    dataPlot(0, 6) = (*v)(2);

    dataPlot(0, 7) = (*p)(0);
    dataPlot(0, 8) = (*y)(0);
    dataPlot(0, 9) = (*y)(1);
    dataPlot(0, 10) = (*u)(0);
    dataPlot(0, 11) = (*u)(1);
    dataPlot(0, 12) = (*lambda)(0);
    dataPlot(0, 13) = (*lambda)(1);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    while(s->hasNextEvent())
    {
      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      
      dataPlot(k, 7) = (*p)(0);
      dataPlot(k, 8) = (*y)(0);
      dataPlot(k, 9) = (*y)(1);
      dataPlot(k, 10) = (*u)(0);
      dataPlot(k, 11) = (*u)(1);
      dataPlot(k, 12) = (*lambda)(0);
      dataPlot(k, 13) = (*lambda)(1);
      
      s->nextStep();
      k++;
    }
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << endl <<  "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms" << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("ImpactingStick.dat", "ascii", dataPlot, "noDim");
    ioMatrix::write("ImpactingStick.ref", "ascii", dataPlot);
    double error=0.0, eps=1e-12;
    if((error=ioMatrix::compareRefFile(dataPlot, "ImpactingStick.ref", eps)) >= 0.0
        && error > eps)
      return 1;

  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }



}
