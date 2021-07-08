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

/*!\file BouncingBallTS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/
#include <chrono>
#include <SiconosBodies.hpp>
#include <SiconosKernel.hpp>

#include <RigidBody2dDS.hpp>
#include "Contact2dR.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10;                  // final computation time
    double h = 0.005;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    double R = 0.5; // Ball radius
    double m = 1; // Ball mass
    double g = 9.81; // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." <<  endl;

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;

    // -- The dynamical system --
    SP::RigidBody2dDS ball(new RigidBody2dDS(q0, v0, Mass));

    SP::SiconosVector q01(new SiconosVector(nDof));
    SP::SiconosVector v01(new SiconosVector(nDof));
    (*q01)(0) = position_init+2*R+0.1;
    (*v01)(0) = velocity_init;



    SP::RigidBody2dDS ball1(new RigidBody2dDS(q01, v01, Mass));

    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);
    ball1->setFExtPtr(weight);
    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // // Interaction ball-floor
    // //
    // SP::SimpleMatrix H(new SimpleMatrix(1, nDof));
    // (*H)(0, 0) = 1.0;
    // SP::Relation relation(new LagrangianLinearTIR(H));

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(e,0.0,0.1,2));

    SP::Contact2dR relation(new Contact2dR());

    SP::Interaction inter(new Interaction(nslaw, relation));

    SP::Contact2dR relation1(new Contact2dR());

    SP::Interaction inter1(new Interaction(nslaw, relation1));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);
    bouncingBall->insertDynamicalSystem(ball1);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball);
    bouncingBall->link(inter1, ball1, ball);


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
    SP::TimeStepping s(new TimeStepping(bouncingBall, t, OSI, osnspb));

    // =========================== End of model definition ===========================

    // ================================= Computation =================================


    int N = ceil((T - t0) / h)+1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(1);
    SP::SiconosVector lambda = inter->lambda(1);
    SP::SiconosVector q1 = ball1->q();
    SP::SiconosVector v1 = ball->velocity();
    SP::SiconosVector p1 = ball->p(1);
    SP::SiconosVector lambda1 = inter1->lambda(1);


    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = (*q1)(0);
    dataPlot(0, 6) = (*v1)(0);
    dataPlot(0, 7) = (*p1)(0);
    dataPlot(0, 8) = (*lambda1)(0);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    SP::SiconosVector pc = relation->pc1();
    SP::SiconosVector nnc = relation->nc();

    SP::SiconosVector pc1 = relation1->pc1();
    SP::SiconosVector pc2 = relation1->pc2();
    SP::SiconosVector nnc1 = relation1->nc();




    while(s->hasNextEvent())
    {
      // a fake contact detection
      (*pc)(0) = -R + (*q)(0);
      (*pc)(1) = 0.0 + (*q)(1);
      (*nnc)(0) = 1.0;
      (*nnc)(1) = 0.0;

      (*pc1)(0) = -R + (*q1)(0);;
      (*pc1)(1) = 0.0 + (*q1)(1);

      (*pc2)(0) = R + (*q)(0);;
      (*pc2)(1) = 0.0 + (*q)(1);;
      (*nnc1)(0) = 1.0;
      (*nnc1)(1) = 0.0;

      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q1)(0);
      dataPlot(k, 6) = (*v1)(0);
      dataPlot(k, 7) = (*p1)(0);
      dataPlot(k, 8) = (*lambda1)(0);
      s->nextStep();
      k++;

    }
    cout  << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time : "  << endl;
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << "Computation time : " << elapsed << " ms" << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("Ball2D.dat", "ascii", dataPlot, "noDim");
    double error=0.0, eps=1e-12;
    if((error=ioMatrix::compareRefFile(dataPlot, "Ball2D.ref", eps)) >= 0.0
        && error > eps)
      return 1;

  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }



}
