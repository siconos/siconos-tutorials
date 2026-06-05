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

/*
  A Ball bouncing on the moving ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>  // For std::numbers::pi

// #define WITH_PROJ
#define WITH_FC3D
#include "UserRelation.hpp"

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDof = 3;                // degrees of freedom for the ball
    unsigned int qDim = 7;       // degrees of freedom for the ball
    unsigned int nDim = 6;       // degrees of freedom for the ball
    double t0 = 0;               // initial computation time
    double T = 10.0;             // final computation time
    double h = 0.005;            // time step
    double position_init = 1.0;  // initial position for lowest bead.
    double velocity_init = 2.0;  // initial velocity for lowest bead.
    double omega_initx = 0.0;
    double omega_initz = 0.0;  // initial velocity for lowest bead.
    double theta = 0.5;        // theta for MoreauJeanOSI integrator
    double m = 1;              // Ball mass
    double g = 9.81;           // Gravity
    double radius = 0.0;
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    std::cout << "====> Model loading ...\n";

    // -- Initial positions and velocities --

    siconos::algebra::SiconosVector q0{qDim};
    siconos::algebra::SiconosVector v0{nDim};
    q0.setZero();
    v0.setZero();
    siconos::algebra::SiconosDenseMatrix I = Eigen::MatrixXd::Identity(3, 3);

    q0(0) = position_init;
    /*initial quaternion equal to (1,0,0,0)*/
    q0(3) = 1.0;

    v0(0) = velocity_init;
    v0(3) = omega_initx;
    v0(5) = omega_initz;
    // -- The dynamical system --
    auto ball = std::make_shared<siconos::modeling::NewtonEulerDS>(q0, v0, m, I,
                                                                   siconos::algebra::alias_t);

    // -- Set external forces (weight) --
    siconos::algebra::SiconosVector weight{nDof};
    weight(0) = -m * g;
    ball->setConstantFext(weight, siconos::algebra::alias_t);

    // -- Moving Plane --

    // -- Initial positions and velocities --
    siconos::algebra::SiconosVector q02{qDim};
    siconos::algebra::SiconosVector v02{nDim};
    v02.setZero();
    q02.setZero();
    q02(3) = 1.0;
    // -- The dynamical system --
    auto movingplane = std::make_shared<siconos::modeling::NewtonEulerDS>(
        q02, v02, m, I, siconos::algebra::alias_t);

    // // -- Set external forces (weight) --
    movingplane->setConstantFext(weight, siconos::algebra::alias_t);

    auto bd = std::make_shared<siconos::modeling::BoundaryCondition>(
        siconos::modeling::BoundaryCondition::Indices{0});
    bd->setComputePrescribedVelocityFunction(
        [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
          result.setConstant(2. + cos(0.5 * std::numbers::pi * time));
        });
    movingplane->setBoundaryConditions(bd);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

#ifdef WITH_FC3D
    auto nslaw0 = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(e, e, 0.6, 3);
#else
    auto nslaw0 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
#endif

    auto relation0 = std::make_shared<user_defined::my_NewtonEulerR>(radius);
    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw0, relation0);

    // -------------
    // --- Model ---
    // -------------
    auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);
    bouncingBall->insertDynamicalSystem(movingplane);

    // link the interaction and the dynamical system
    // r
    bouncingBall->link(inter, ball, movingplane);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
#ifdef WITH_PROJ
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanDirectProjectionOSI>(theta);
#else
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
#endif
    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GenericMechanical>();
#ifdef WITH_PROJ
    auto osnspb_pos =
        std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>(
            SICONOS_MLCP_ENUM, 1.0);
#endif
    // -- (4) Simulation setup with (1) (2) (3)
#ifdef WITH_PROJ
    auto s = std::make_shared<siconos::simulation::TimeSteppingDirectProjection>(
        bouncingBall, t, OSI, osnspb, osnspb_pos);
    s->setProjectionMaxIteration(20);
    s->setConstraintTolUnilateral(1e-08);
    s->setConstraintTol(1e-08);
#else
    auto s = std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI, osnspb);
#endif
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(10);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    siconos::algebra::Index outputSize = 19;
    siconos::algebra::SiconosDenseMatrix dataPlot(N + 1, outputSize);

    auto q = ball->q();
    auto v = ball->twist();
    auto p = ball->p(1);
    auto qplane = movingplane->q();
    auto vplane = movingplane->twist();
    auto pplane = movingplane->p(1);

    auto lambda = inter->lambda(1);
    auto y = inter->y(0);

    auto reaction = movingplane->reactionToBoundaryConditions();

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = acos((*q)(3));
    dataPlot(0, 6) = relation0->contactForce().norm();

    dataPlot(0, 7) = (*q)(0);
    dataPlot(0, 8) = (*q)(1);
    dataPlot(0, 9) = (*q)(2);
    dataPlot(0, 10) = (*q)(3);
    dataPlot(0, 11) = (*q)(4);
    dataPlot(0, 12) = (*q)(5);
    dataPlot(0, 13) = (*q)(6);

    dataPlot(0, 14) = (*v)(1);
    dataPlot(0, 15) = (*v)(2);
    dataPlot(0, 16) = (*qplane)(2);
    dataPlot(0, 17) = (*vplane)(2);
    dataPlot(0, 18) = (*reaction)(0);

    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();
    dataPlot(k, 6) = relation0->contactForce().norm();
    while (s->hasNextEvent() && k < 5000000) {
      //      s->computeOneStep();
      s->advanceToEvent();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = acos((*q)(3));
      dataPlot(k, 6) = relation0->contactForce().norm();
      dataPlot(k, 7) = (*q)(0);
      dataPlot(k, 8) = (*q)(1);
      dataPlot(k, 9) = (*q)(2);
      dataPlot(k, 10) = (*q)(3);
      dataPlot(k, 11) = (*q)(4);
      dataPlot(k, 12) = (*q)(5);
      dataPlot(k, 13) = (*q)(6);
      dataPlot(k, 14) = (*v)(1);
      dataPlot(k, 15) = (*v)(2);

      dataPlot(k, 16) = (*qplane)(0);
      dataPlot(k, 17) = (*vplane)(0);

      dataPlot(k, 18) = (*reaction)(0);

      s->nextStep();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("BallNewtonEulerOnMovingPlane.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // Comparison with a reference file
    std::cout << "====> Comparison with a reference file ...\n";
    double error = 0.0, eps = 1e-10;
#ifdef WITH_PROJ
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "resultNETS-WITHPROJ.ref.ref",
                                                      eps)) > eps)
      return 1;
#else
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "BallNewtonEulerOnMovingPlane.ref", eps)) > eps)
      return 1;
#endif
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
