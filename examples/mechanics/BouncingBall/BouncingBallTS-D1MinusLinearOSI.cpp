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

/* D1MinusLinearOSI-Time-Stepping version
  T. Schindler, V. Acary

  A Ball bouncing on the ground.
  Direct description of the model
  Simulation with a D1MinusLinearOSI-Time-Stepping scheme.
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
    double h = 5e-4;             // time step
    double hplot = 0.005;        // plot step size (larger than time step)
    double position_init = 1.0;  // initial position for lowest bead.
    double velocity_init = 0.0;  // initial velocity for lowest bead.
    double R = 0.1;              // Ball radius
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
    double e = 0.9;

    // Interaction ball-floor
    //
    Matrix H{1, nDof};
    H.setZero();
    H(0, 0) = 1.0;
    auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H);
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::D1MinusLinearOSI>(
        siconos::integrators::D1MinusLinearOSI::Type::halfexplicit_acceleration_level);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) One step non smooth problems
    auto impact =
        std::make_shared<siconos::nonsmooth_formulations::LCP>();  // impulse right limit
                                                                   // right side
    auto force =
        std::make_shared<siconos::nonsmooth_formulations::LCP>();  // contact force right
                                                                   // limit left side

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeSteppingD1Minus>(bouncingBall, t, 2);
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
    s->insertNonSmoothProblem(force, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h);  // Number of time steps
    // int Nplot = (int)((T - t0) / hplot); // Number of plot steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 7;
    Matrix dataPlot(N + 1, outputSize);

    auto q = ball->q_read();
    auto v = ball->velocity_read();
    auto p = ball->p_read(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = q(0);
    dataPlot(0, 2) = v(0);
    dataPlot(0, 3) = p(0);
    dataPlot(0, 4) = 0;
    dataPlot(0, 5) = 0;
    dataPlot(0, 6) = 0;
    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      s->advanceToEvent();
      // siconos::algebra::print(*ball);
      // --- Get values to be plotted ---
      //  if (fmod(s->nextTime(), hplot) < h)
      {
        // std::cout << "k=" << k <<std::endl;
        dataPlot(k, 0) = s->nextTime();
        dataPlot(k, 1) = q(0);
        dataPlot(k, 2) = v(0);
        dataPlot(k, 3) = p(0);
        dataPlot(k, 4) = (*inter->lambda(1))(0);
        dataPlot(k, 5) = (*ball->p(2))(0);
        dataPlot(k, 6) = (*inter->lambda(2))(0);
        k++;
      }

      s->processEvents();
      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("BouncingBallTS-D1MinusLinearOSI.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "BouncingBallTS-D1MinusLinearOSI.ref", eps)) >= eps)
      return 1;

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
