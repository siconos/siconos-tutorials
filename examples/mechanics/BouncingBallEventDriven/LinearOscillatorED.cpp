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

/*!\file LinearOscillator.cpp
  \brief \ref EMBouncingBall - C++ input file, Event-Driven version - V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with an Event-Driven scheme.
*/

#include <SiconosKernel.hpp>
#include <chrono>
#include <ratio>

using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;       // degrees of freedom for the ball
    double t0 = 0;               // initial computation time
    double T = 8.5;              // final computation time
    double h = 0.005;            // time step
    double position_init = 1.0;  // initial position for lowest bead.
    double velocity_init = 0.0;  // initial velocity for lowest bead.
    double R = 0.1;              // Ball radius
    double m = 1e-03;            // Ball mass
    double g = 9.81;             // Gravity
    double kx = 1e+02;           // stiffness in x-axis
    double ky = 1e+02;           // stiffness in y-axis
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    std::cout << "====> Model loading ...\n";

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 3. / 5 * m * R * R;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0)(0) = position_init;
    (*q0)(1) = 0.001;
    (*v0)(0) = velocity_init;
    // -- The dynamical system --
    SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q0, v0, Mass));

    // set stiffmatrix
    SP::SiconosMatrix K(new SimpleMatrix(nDof, nDof));
    (*K)(0, 0) = kx;
    (*K)(1, 1) = ky;
    (*K)(2, 2) = 0.;
    ball->setK(*K);

    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(1, nDof));
    (*H)(0, 0) = 1.0;
    SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
    SP::Relation relation0(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(nslaw0, relation0));

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- (1) OneStepIntegrators --
    SP::OneStepIntegrator OSI(new LsodarOSI());

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) Non smooth problem --
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem acceleration(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::EventDriven s(new EventDriven(bouncingBall, t));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_ED_IMPACT);
    s->insertNonSmoothProblem(acceleration, SICONOS_OSNSP_ED_SMOOTH_ACC);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    s->setPrintStat(true);

    int N = 1854;  // Number of saved points: depends on the number of events ...

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    SimpleMatrix dataPlot(N + 1, outputSize);
    auto q = ball->q();
    auto v = ball->velocity();
    auto p = ball->p(1);
    SP::SiconosVector f;
    //   SiconosVector * y = bouncingBall->nonSmoothDynamicalSystem()->interaction(0)->y(0);

    SP::EventsManager eventsManager = s->eventsManager();

    // For the initial time step:
    // time

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = 0.0;

    dataPlot(0, 7) = (*q)(1);
    dataPlot(0, 8) = (*v)(1);

    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    bool nonSmooth = false;
    unsigned int numberOfEvent = 0;
    int k = 0;
    int kns = 0;

    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent() && k < N) {
      s->advanceToEvent();
      f = ball->p(2);
      if (eventsManager->nextEvent()->getType() == 2) nonSmooth = true;

      s->processEvents();
      // If the treated event is non smooth, the pre-impact state has been solved in memory
      // vectors during process.
      if (nonSmooth) {
        dataPlot(k, 0) = s->startingTime();
        dataPlot(k, 1) = ball->qMemory().getSiconosVector(1)(0);
        dataPlot(k, 2) = ball->velocityMemory().getSiconosVector(1)(0);
        dataPlot(k, 3) = (*p)(0);
        dataPlot(k, 4) = (*f)(0);
        k++;
        kns++;
        nonSmooth = false;
      }
      dataPlot(k, 0) = s->startingTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*f)(0);
      dataPlot(k, 5) = (*inter->lambda(1))(0);
      dataPlot(k, 6) = (*inter->lambda(2))(0);
      dataPlot(k, 7) = (*q)(1);
      dataPlot(k, 8) = (*v)(1);

      ++k;
      ++numberOfEvent;
    }
    auto end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "===== End of Event Driven simulation. \n";
    std::cout << numberOfEvent << " events have been processed. ==== \n";
    std::cout << numberOfEvent - kns << " events are of time--discretization type  ==== \n";
    std::cout << kns << " events are of nonsmooth type  ==== \n\n";
    std::cout << "\nComputation time : " << elapsed << " ms\n";
    // --- Output files ---
    std::cout << "====> Output file writing ...\n\n";
    dataPlot.resize(k, outputSize);
    ioMatrix::write("LinearOscillatorED.dat", "ascii", dataPlot, "noDim");

    double error = 0.0, eps = 1e-10;
    if ((error = ioMatrix::compareRefFile(dataPlot, "LinearOscillatorED.ref", eps)) >= 0.0 &&
        error > eps)
      return 1;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}