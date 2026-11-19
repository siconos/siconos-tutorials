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

// =============================== Cam Follower (1DOF Impact System)
//
// The Cam Follower system is modelled as a Generalised Langrangian System impacting against a
// fixed wall the moving constraint (i.e. a rotational cam) is modelled as an input force
//
// Direct description of the model.
//
// Keywords: LagrangianLinearDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// ======================================================================================================

#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>

#include "CamState.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

int main(int argc, char* argv[]) {
  double rpm = 358;
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int dsNumber = 1;    // the Follower and the ground
    int nDof = 1;                 // degrees of freedom for the ball
    double t0 = 0;                // initial computation time
    double T = 1;                 // final computation time
    double h = 0.0001;            // time step
    double position_init = 0.40;  // initial position for lowest bead.
    double velocity_init = 0.4;   // initial velocity for lowest bead.
    double theta = 0.5;           // theta for MoreauJeanOSI integrator
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    Matrix Mass{nDof, nDof};
    Mass.setZero();
    Mass(0, 0) = user_defined::mass;
    Vector q0{nDof};
    Vector velocity0{nDof};
    q0.setZero();
    q0(0) = position_init;
    velocity0.setZero();
    velocity0(0) = velocity_init;
    auto lds = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        q0, velocity0, Mass, siconos::algebra::alias_t);
    Matrix K{nDof, nDof};
    K.setZero();
    K(0, 0) = 1430.8;
    lds->setStiffnessMatrix(K, siconos::algebra::alias_t);
    lds->setComputeFextFunction(
        [mass = user_defined::mass, gravity = user_defined::gravity](
            double time, Eigen::Ref<siconos::algebra::MapVectorType> fext) {
          fext(0) = -mass * gravity;
        });

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.8;

    // Interaction Follower-floor
    //
    auto H = std::make_shared<Matrix>(1, nDof);
    (*H)(0, 0) = 1.0;
    auto nslaw0 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation0 = std::make_shared<siconos::modeling::LagrangianRheonomousR>();

    auto hfunc = [rpm](const siconos::algebra::BlockVector& pos, double time,
                       Eigen::Ref<siconos::algebra::MapVectorType> y) {
      double CamEqForce, CamPosition, CamVelocity, CamAcceleration;
      // CamEqForce =
      //     user_defined::CamState(time, rpm, CamPosition, CamVelocity, CamAcceleration);
      //  y[0] = q[0] - CamPosition;
      y(0) = pos(0);
    };
    relation0->setComputehFunction(hfunc);

    auto jachq = [](const siconos::algebra::BlockVector& pos, double time,
                    Eigen::Ref<siconos::algebra::MapType> result) { result(0, 0) = 1; };
    relation0->setComputeContactMatrixFunction(jachq);

    auto dhdt = [rpm](const siconos::algebra::BlockVector& pos, double time,
                      Eigen::Ref<siconos::algebra::MapVectorType> result) {
      // double CamEqForce, CamPosition, CamVelocity, CamAcceleration;

      // CamEqForce =
      //     user_defined::CamState(time, rpm, CamPosition, CamVelocity, CamAcceleration);
      // result[0] = -CamVelocity;
      result.setZero();
    };
    relation0->setComputedhdtFunction(dhdt);

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw0, relation0);

    // -------------
    // --- Model ---
    // -------------

    auto Follower = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    Follower->insertDynamicalSystem(lds);
    Follower->link(inter, lds);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- OneStepIntegrator --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    // -- OneStepNsProblem --
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_QP);

    // solver
    // osnspb->numericsSolverOptions()->solverId=SICONOS_LCP_QP;

    // max number of iterations
    osnspb->numericsSolverOptions()->iparam[0] = 101;

    // tolerance
    osnspb->numericsSolverOptions()->dparam[0] = 1e-6;

    auto S = std::make_shared<siconos::simulation::TimeStepping>(Follower, t, OSI, osnspb);
    cout << "=== End of model loading === \n";
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int k = 0;
    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 8;
    Matrix DataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    DataPlot(k, 0) = t0;
    DataPlot(k, 1) = (*lds->q())(0);
    DataPlot(k, 2) = (*lds->velocity())(0);
    DataPlot(k, 3) = (*inter->lambda(1))(0);
    DataPlot(k, 4) = lds->fext()(0);

    // State of the Cam
    //    double rpm=358;
    double CamEqForce, CamPosition, CamVelocity, CamAcceleration;

    CamEqForce = user_defined::CamState(t0, rpm, CamPosition, CamVelocity, CamAcceleration);
    // Position of the Cam
    DataPlot(k, 5) = CamPosition;
    // Velocity of the Cam
    DataPlot(k, 6) = CamVelocity;
    // Acceleration of the Cam
    DataPlot(k, 7) = CamPosition + (*lds->q())(0);
    auto start = std::chrono::system_clock::now();
    // --- Time loop ---
    cout << "Start computation ... \n";
    while (k < N) {
      // get current time step
      k++;
      // solve ...
      S->computeOneStep();

      // --- Get values to be plotted ---

      DataPlot(k, 0) = S->nextTime();
      DataPlot(k, 1) = (*lds->q())(0);
      DataPlot(k, 2) = (*lds->velocity())(0);
      DataPlot(k, 3) = (*inter->lambda(1))(0);
      DataPlot(k, 4) = lds->fext()(0);

      CamEqForce = user_defined::CamState(S->nextTime(), rpm, CamPosition, CamVelocity,
                                          CamAcceleration);
      DataPlot(k, 5) = CamPosition;
      DataPlot(k, 6) = CamVelocity;
      DataPlot(k, 7) = CamPosition + (*lds->q())(0);
      // transfer of state i+1 into state i and time incrementation
      S->nextStep();
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "\nEnd of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms\n";

    // --- Output files ---
    siconos::algebra::io::write("CamFollower-Rheonomous.dat", DataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(DataPlot, "CamFollower-Rheonomous.ref",
                                                      eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
