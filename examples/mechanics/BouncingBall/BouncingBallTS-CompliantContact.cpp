/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*
  V. Acary, F. Perignon.

  A Ball bouncing between a roof and a ground. Compliant contact with the ground.
  Simulation with a Time-Stepping scheme.
*/

#include <SiconosKernel.hpp>
#include <chrono>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDof = 3;                // degrees of freedom for the ball
    double t0 = 0.;              // initial computation time
    double T = 10.;              // final computation time
    double h = 0.0005;           // time step
    double position_init = 1.0;  // initial position for lowest bead.
    double velocity_init = 0.0;  // initial velocity for lowest bead.
    double theta = 0.5;          // theta for the integrator
    double R = 0.1;              // Ball radius
    double height = 1.0;         // height to the roof
    double m = 1;                // Ball mass
    double g = 9.81;             // Gravity

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
    double e = 1.0;

    // Interaction ball-roof with impact
    //
    Matrix H{1, nDof};
    H.setZero();
    H(0, 0) = -1.0;
    Vector b{1};
    b << height - R - .2;
    auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H, b);
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // -- nslaw --
    double compliance = 0.01;

    // Interaction ball-floor with a compliant spring -
    // y = q0 + compliance.lambda - R
    Matrix Hfloor{1, nDof};
    Hfloor.setZero();
    Hfloor(0, 0) = 1.0;
    Matrix Kfloor{1, 1};
    Kfloor << compliance;
    Vector bfloor{1};
    bfloor << -R;

    auto nslawfloor = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(1);
    auto relationfloor = std::make_shared<siconos::modeling::LagrangianCompliantLinearTIR>(
        Hfloor, Kfloor, bfloor);
    auto interfloor =
        std::make_shared<siconos::modeling::Interaction>(nslawfloor, relationfloor);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball);

    // link the interaction and the dynamical system
    bouncingBall->link(interfloor, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI, osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 6;
    Matrix dataPlot(N + 1, outputSize);

    auto q = ball->q();
    auto v = ball->velocity();
    auto p = ball->p(1);
    auto lambda = inter->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = (*interfloor->y(0))(0);
    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*interfloor->y(0))(0);
      s->nextStep();
      siconos::tools::progressBar((double)k / N);
      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("BouncingBallTS-CompliantContact.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "BouncingBallTS-CompliantContact.ref", eps)) > eps)
      return 1;
    else
      return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
