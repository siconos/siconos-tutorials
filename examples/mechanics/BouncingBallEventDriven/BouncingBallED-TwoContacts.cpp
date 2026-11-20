/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
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

/*!\file BouncingBallED.cpp
  \brief \ref EMBouncingBall - C++ input file, Event-Driven version - V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with an Event-Driven scheme.
*/

#include <SiconosKernel.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <chrono>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDof = 3;                 // degrees of freedom for the ball
    double t0 = 0;                // initial computation time
    double T = 10.0;              // final computation time
    double h = 0.01;              // time step
    double position_init = 1.0;   // initial position for lowest bead.
    double velocity_init = 10.0;  // initial velocity for lowest bead.
    double Heightbox = 1.5;
    double R = 0.1;   // Ball radius
    double m = 1;     // Ball mass
    double g = 10.0;  // Gravity

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    std::cout << "====> Model loading ...\n";

    Matrix mass{nDof, nDof};
    mass.setZero();
    mass(0, 0) = m;
    mass(1, 1) = m;
    mass(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --
    Vector q0{nDof};
    q0.setZero();
    q0(0) = position_init;
    Vector v0{nDof};
    v0.setZero();
    v0(0) = velocity_init;

    // -- The dynamical system --
    auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        q0, v0, mass, siconos::algebra::alias_t);

    // -- Set external forces (weight) --
    Vector weight{nDof};
    weight.setZero();
    weight(0) = -m * g;
    ball->setConstantFext(weight, siconos::algebra::alias_t);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.8;  // Warning this example does not work with e=0.0

    // Interaction ball-floor-ceiling
    //
    auto H1 = std::make_shared<Matrix>(1, nDof);
    (*H1)(0, 0) = 1.0;
    auto E1 = std::make_shared<Vector>(1);
    (*E1)(0) = 0.0;  //-1.0*R;
    //
    auto H2 = std::make_shared<Matrix>(1, nDof);
    (*H2)(0, 0) = -1.0;
    auto E2 = std::make_shared<Vector>(1);
    (*E2)(0) = Heightbox;  //- R;
    // impact law
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    // Interaction at contact 1 (ball-floor)
    auto relation1 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H1, *E1);
    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation1);
    // Interaction at contact 2 (ball-ceiling)
    auto relation2 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H2, *E2);
    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation2);
    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->link(inter1, ball);
    bouncingBall->link(inter2, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::LsodarOSI>();

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problems
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto acceleration = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::EventDriven>(bouncingBall, t);
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_ED_IMPACT);
    s->insertNonSmoothProblem(acceleration, siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC);

    // =========================== End of model definition
    // ===========================

    // ================================= Computation
    // =================================

    s->setPrintStat(true);

    int N = 1161;  // Number of saved points: depends on the number of events ...
    int ll = 0;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    Matrix dataPlot{N, outputSize};
    auto q = ball->q();             // ball position
    auto v = ball->velocity();      // ball velocity
    std::shared_ptr<Vector> gamma;  // ball acceleration
    std::shared_ptr<Vector> f;  // resultant force deduced from the LCP at acceleration level
    std::shared_ptr<Vector> p;  // resultant force deduced from the LCP at velocity level

    auto y1 = inter1->y(0);
    auto y2 = inter2->y(0);
    //   SiconosVector * y = bouncingBall->nonSmoothDynamicalSystem()->interaction(0)->y(0);

    auto eventsManager = s->eventsManager();

    OSI->display();
    // For the initial time step:
    // time

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = 0.0;
    dataPlot(0, 4) = 0;
    dataPlot(0, 5) = 1.0;
    dataPlot(0, 6) = 0.5;
    dataPlot(0, 7) = -10.0;
    dataPlot(0, 8) = 0.0;

    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    bool nonSmooth = false;
    unsigned int numberOfEvent = 0;
    int k = 1;

    //    s->setTolerance(1e-10);
    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      s->advanceToEvent();  // run simulation from one event to the next
      f = ball->p(2);       // resultant force deduced from the LCP at acceleration level
      p = ball->p(1);
      gamma = ball->acceleration();
      y1 = inter1->y(0);
      y2 = inter2->y(0);
      if (eventsManager->nextEvent()->getType() == siconos::simulation::EventType::NS)
        nonSmooth = true;

      s->processEvents();  // process events
      // If the treated event is non smooth, the pre-impact state has been solved in memory
      // vectors during process.
      if (nonSmooth)  // if the event is nonsmooth
      {
        dataPlot(k, 0) = s->startingTime();  // get the time at nonsmooth event
        dataPlot(k, 1) = ball->qMemory().getSiconosVector(1)(0);
        dataPlot(k, 2) = ball->velocityMemory().getSiconosVector(1)(0);
        k++;
        nonSmooth = false;

        dataPlot(k, 4) = 1;
        ++ll;
        //         cout << "========================================\n";
        //         cout << "Nonsmooth event\n";
      }
      dataPlot(k, 0) = s->startingTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 5) = (*y1)(0);
      dataPlot(k, 6) = (*y2)(0);
      dataPlot(k, 4) = 0;
      dataPlot(k, 7) = (*gamma)(0);
      dataPlot(k, 8) = (*f)(0);

      // cout << "========================================\n";
      // cout << " time: " << s->startingTime() << endl;
      // cout << "ball position: " << (*q)(0) << endl;
      // cout << "ball velocity: " << (*v)(0) << endl;
      // cout << "gap at contact 1: " << (*y1)(0) << endl;
      // cout << "gap at contact 2: " << (*y2)(0) << endl;
      //
      k++;
      ++numberOfEvent;
    }

    // --- Output files ---
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "===== End of Event Driven simulation. \n";
    std::cout << numberOfEvent << " events have been processed. ==== \n";
    std::cout << numberOfEvent - ll << " events are of time--discretization type  ==== \n";
    std::cout << ll << " events are of nonsmooth type  ==== \n\n";
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    siconos::algebra::io::write("BouncingBallED-TwoContacts.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // relative tolerance: the event-driven integration (Lsodar, event location) differs
    // with the platform (libm, blas) by up to 1e-9 in absolute value
    double error = 0.0, eps = 1e-8 * dataPlot.cwiseAbs().maxCoeff();
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "BouncingBallED-TwoContacts.ref", eps)) > eps)
      return 1;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
