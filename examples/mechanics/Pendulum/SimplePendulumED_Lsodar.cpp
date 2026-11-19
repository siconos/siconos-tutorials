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

// =============================== Simple Pendulum Example ===============================
//
// Author: Vincent Acary
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

constexpr double gravity = 9.8100;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================
    // User-defined main parameters
    int nDof = 2;                               // degrees of freedom for robot arm
    double t0 = 0;                              // initial computation time
    double T = 10.0;                            // final computation time
    double h = 0.01;                            // time step
    double L = 1.0;                             // Length of the pendulum
    double InitAngle = std::numbers::pi / 3.0;  // Initial inclination angle
    double m = 1.0;                             // Mass of the pendulum
    double e = 0.9;                             // nslaw
    double _rho = 0.99;

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // --- DS: Simple Pendulum ---

    // Initial position (angles in radian)
    Vector q0{nDof};
    q0.setZero();
    Vector v0{nDof};
    v0.setZero();
    q0(0) = L * sin(InitAngle);
    q0(1) = L * cos(InitAngle);
    Matrix mass{nDof, nDof};

    mass(0, 0) = m;
    mass(1, 1) = m;
    auto simplependulum = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        q0, v0, mass, siconos::algebra::alias_t);

    Vector ForceExtern{nDof};
    ForceExtern.setZero();
    ForceExtern(1) = m * gravity;
    simplependulum->setConstantFext(ForceExtern, siconos::algebra::alias_t);

    // -------------------
    // --- Interactions---
    // -------------------
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation->setComputehFunction([L](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = pow(L, 2) - (pow(q(0), 2) + pow(q(1), 2));
    });

    relation->setComputeContactMatrixFunction(
        [](const siconos::algebra::BlockVector& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = -2.0 * q(0);
          result(0, 1) = -2.0 * q(1);
        });

    relation->setComputeContactMatrixDotFunction(
        [](const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = -2.0 * qdot(0);
          result(0, 1) = -2.0 * qdot(1);
        });

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // -------------
    // --- Model ---
    // -------------

    auto Pendulum = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    Pendulum->insertDynamicalSystem(simplependulum);
    Pendulum->link(inter, simplependulum);

    // ----------------
    // --- Simulation ---
    // ----------------

    // 1. Time discretization
    auto TimeDiscret = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // 2. Integration solver for one step
    auto OSI = std::make_shared<siconos::integrators::LsodarOSI>();
    // 3. Nonsmooth problem
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto acceleration = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    // 4. Simulation with (1), (2), (3)
    auto EDscheme = std::make_shared<siconos::simulation::EventDriven>(Pendulum, TimeDiscret);
    EDscheme->insertIntegrator(OSI);
    EDscheme->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_ED_IMPACT);
    EDscheme->insertNonSmoothProblem(acceleration,
                                     siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC);

    // =========================== End of model definition ===========================

    // auto lsodar = std::static_pointer_cast<LsodarOSI>(OSI);
    // lsodar->setMinMaxStepSizes(9.5e-4,1.0e-3);
    // lsodar->setTol(1,1.0e-3,1.0e-6);
    // lsodar->setMaxOrder(2, 2);

    // ================================= Computation =================================

    auto eventsManager =
        EDscheme->eventsManager();            // ponters point to the "eventsManager" object
    auto _q = simplependulum->q();            // pointer points to the position vector
    auto _qdot = simplependulum->velocity();  // pointer points to the velocity
    simplependulum->initRhs(t0);
    simplependulum->computeRhs(t0);
    auto _qddot = simplependulum->acceleration();
    auto _g = inter->y(0);
    auto indexSet0 = Pendulum->topology()->indexSet(0);
    std::cout << "Size of IndexSet0: " << indexSet0->size() << "\n";
    //-------------------- Save the output during simulation
    //---------------------------------------------------------
    unsigned int N = ceil(T / h) + 1;  // Number of points to be saved
    Matrix DataPlot(N, 10);
    //------------- At the initial time
    //-----------------------------------------------------------------------------
    DataPlot(0, 0) = Pendulum->t0();
    DataPlot(0, 1) = (*_q)(0);      // Position X
    DataPlot(0, 2) = (*_q)(1);      // Position Y
    DataPlot(0, 3) = (*_qdot)(0);   // Velocity Vx
    DataPlot(0, 4) = (*_qdot)(1);   // Velocity Vy
    DataPlot(0, 5) = (*_qddot)(0);  // Acceleration ax
    DataPlot(0, 6) = (*_qddot)(1);  // Acceleration ay
    DataPlot(0, 7) = (*_g)(0);      // Contraint in position
    DataPlot(0, 8) = 0.0;           // Constraint in velocity
    DataPlot(0, 9) = 0.0;           // Reaction force

    //----------------------------------- Simulation starts
    //----------------------------------------------------------
    std::cout << "====> Start computation ... \n";
    bool NSEvent = false;
    unsigned int NumberNSEvent = 0;
    unsigned int k = 0;

    auto start = std::chrono::system_clock::now();
    while ((EDscheme->hasNextEvent()) && (k < N)) {
      EDscheme->advanceToEvent();  // lead the simulation run from one event to the next
      //---------- detect the statue of the current event ------------------------------------
      if (eventsManager->nextEvent()->getType() ==
          siconos::simulation::EventType::NS)  // the current event is non-smooth
      {
        NSEvent = true;
      };
      EDscheme->processEvents();  // process the current event
      //------------------- get data at the beginning of non-smooth events
      //---------------------------
      auto _gdot = inter->y(1);
      auto _lambda = inter->lambda(2);

      if (NSEvent) {
        DataPlot(k, 0) = EDscheme->startingTime();  // instant at non-smooth event

        const auto& qMemory = simplependulum->qMemory().getSiconosVector(1);
        const auto& qdotMemory = simplependulum->velocityMemory().getSiconosVector(1);
        DataPlot(k, 1) = qMemory(0);
        DataPlot(k, 2) = qMemory(1);
        DataPlot(k, 3) = qdotMemory(0);
        DataPlot(k, 4) = qdotMemory(1);
        DataPlot(k, 5) = (*_qddot)(0);   // Acceleration ax
        DataPlot(k, 6) = (*_qddot)(1);   // Acceleration ay
        DataPlot(k, 7) = (*_g)(0);       // Contraint in position
        DataPlot(k, 8) = (*_gdot)(0);    // Constraint in velocity
        DataPlot(k, 9) = (*_lambda)(0);  // Reaction force
        k++;
        ++NumberNSEvent;

        NSEvent = false;  // The next event is maybe smooth
      } else {
        //-------------------- get data at smooth events or at the end of non-smooth events
        //---------------
        DataPlot(k, 0) = EDscheme->startingTime();
        DataPlot(k, 1) = (*_q)(0);       // Position X
        DataPlot(k, 2) = (*_q)(1);       // Position Y
        DataPlot(k, 3) = (*_qdot)(0);    // Velocity Vx
        DataPlot(k, 4) = (*_qdot)(1);    // Velocity Vy
        DataPlot(k, 5) = (*_qddot)(0);   // Acceleration ax
        DataPlot(k, 6) = (*_qddot)(1);   // Acceleration ay
        DataPlot(k, 7) = (*_g)(0);       // Contraint in position
        DataPlot(k, 8) = (*_gdot)(0);    // Constraint in velocity
        DataPlot(k, 9) = (*_lambda)(0);  // Reaction force
        k++;
      }
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("SimplePendulumED_Lsodar.dat", DataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(DataPlot, "SimplePendulumED_Lsodar.ref",
                                                      eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
