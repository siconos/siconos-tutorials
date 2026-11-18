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

#include <SiconosKernel.hpp>
#include <SiconosMatrix.hpp>
#include <chrono>

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDof = 3;                // degrees of freedom for the ball
    double t0 = 0;               // initial computation time
    double T = 10;               // final computation time
    double h = 0.005;            // time step
    double position_init = 1.0;  // initial position for lowest bead.
    double velocity_init = 0.0;  // initial velocity for lowest bead.
    double theta = 0.5;          // theta for MoreauJeanOSI integrator
    double R = 0.5;              // Ball radius
    double m = 1;                // Ball mass
    double g = 9.81;             // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    std::cout << "====> Model loading ...\n";

    siconos::algebra::SiconosDenseMatrix mass{nDof, nDof};
    mass.setZero();
    mass(0, 0) = m;
    mass(1, 1) = m;
    mass(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --
    siconos::algebra::SiconosVector q0{nDof};
    q0.setZero();
    q0(0) = position_init;
    siconos::algebra::SiconosVector v0{nDof};
    v0.setZero();
    v0(0) = velocity_init;

    // -- The dynamical system --
    auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        q0, v0, mass, siconos::algebra::alias_t);

    siconos::algebra::SiconosVector q01{nDof};
    siconos::algebra::SiconosVector v01{nDof};
    q01(0) = position_init + 2 * R + 0.1;
    v01(0) = velocity_init;

    auto ball1 = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        q01, v01, mass, siconos::algebra::alias_t);

    // -- Set external forces (weight) --
    siconos::algebra::SiconosVector weight{nDof};
    weight.setZero();
    weight(0) = -m * g;
    ball->setConstantFext(weight, siconos::algebra::alias_t);
    ball1->setConstantFext(weight, siconos::algebra::alias_t);
    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);

    auto relation = std::make_shared<siconos::modeling::Lagrangian2dNormalVelocityR>();

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    auto relation1 = std::make_shared<siconos::modeling::Lagrangian2dNormalVelocityR>();

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation1);

    // -------------
    // --- Model ---
    // -------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

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
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI, osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h) + 1;  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    siconos::algebra::SiconosDenseMatrix dataPlot(N, outputSize);

    auto q = ball->q_read();
    auto v = ball->velocity_read();
    auto p = ball->p_read(1);
    auto lambda = inter->lambda(1);
    auto q1 = ball1->q_read();
    auto v1 = ball->velocity_read();
    auto p1 = ball->p_read(1);
    auto lambda1 = inter1->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = q(0);
    dataPlot(0, 2) = v(0);
    dataPlot(0, 3) = p(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = q1(0);
    dataPlot(0, 6) = v1(0);
    dataPlot(0, 7) = p1(0);
    dataPlot(0, 8) = (*lambda1)(0);
    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    auto start = std::chrono::system_clock::now();

    // For relation
    siconos::algebra::SiconosVector2 pc1;
    siconos::algebra::SiconosVector2 pc2;
    pc2.setZero();
    siconos::algebra::SiconosVector2 normal;
    // For relation1
    siconos::algebra::SiconosVector2 pc1_1;
    siconos::algebra::SiconosVector2 pc2_1;
    siconos::algebra::SiconosVector2 normal_1;

    while (s->hasNextEvent()) {
      // a fake contact detection
      pc1(0) = -R + q(0);
      pc1(1) = q(1);
      normal(0) = 1.0;
      normal(1) = 0.0;
      relation->updateContactPoints(pc1, pc2, normal);

      pc1_1(0) = -R + q1(0);
      pc1_1(1) = q1(1);

      pc2_1(0) = R + q(0);

      pc2_1(1) = q(1);
      normal_1(0) = 1.0;
      normal_1(1) = 0.0;
      relation1->updateContactPoints(pc1_1, pc2_1, normal_1);

      s->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(0);
      dataPlot(k, 2) = v(0);
      dataPlot(k, 3) = p(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = q1(0);
      dataPlot(k, 6) = v1(0);
      dataPlot(k, 7) = p1(0);
      dataPlot(k, 8) = (*lambda1)(0);
      s->nextStep();

      k++;
    }
    std::cout << "End of computation - Number of iterations done: " << k - 1 << std::endl;
    std::cout << "Computation Time : " << std::endl;
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("Ball2D_kernel_only.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "Ball2D_kernel_only.ref",
                                                      eps)) > eps)
      return 1;
    else
      return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
