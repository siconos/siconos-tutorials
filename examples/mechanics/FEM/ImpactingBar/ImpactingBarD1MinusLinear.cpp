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
  V. Acary

  A Bar bouncing on the ground
  Simulation with a Time-Stepping scheme.
*/

#include <LagrangianSparseLinearTIDS.hpp>
#include <SiconosKernel.hpp>
#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>
#include <StorageTools.hpp>
#include <chrono>

#include "UserDefinedParameter.hpp"

#define TS_VELOCITY_LEVEL
using namespace std;
namespace user =
    user_defined;  // To choose the set of parameters used in the current simulation

// #define TS_PROJ
// #define TS_COMBINED

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    std::cout << "====> Model loading ...\n";

    auto ndof = user::nDof;
    double l = user::L / ndof;  // length of an element

    siconos::algebra::SiconosSparseMatrix massMatrix{ndof, ndof};
    siconos::algebra::SiconosSparseMatrix stiffnessMatrix{ndof, ndof};

    massMatrix.insert(0, 0) = 1.0 / 3.0;
    massMatrix.insert(0, 1) = 1.0 / 6.0;
    stiffnessMatrix.insert(0, 0) = 1.0;
    stiffnessMatrix.insert(0, 1) = -1.0;

    for (unsigned int i = 1; i < ndof - 1; i++) {
      massMatrix.insert(i, i) = 2.0 / 3.0;
      massMatrix.insert(i, i - 1) = 1.0 / 6.0;
      massMatrix.insert(i, i + 1) = 1.0 / 6.0;

      stiffnessMatrix.insert(i, i) = 2.0;
      stiffnessMatrix.insert(i, i - 1) = -1.0;
      stiffnessMatrix.insert(i, i + 1) = -1.0;
    }

    massMatrix.insert(ndof - 1, ndof - 1) = 1.0 / 3.0;
    massMatrix.insert(ndof - 1, ndof - 2) = 1.0 / 6.0;

    stiffnessMatrix.insert(ndof - 1, ndof - 2) = -1.0;
    stiffnessMatrix.insert(ndof - 1, ndof - 1) = 1.0;

    // Matrix dampingMatrix{ndof, ndof};
    // dampingMatrix.setZero();
    // auto SparseDamping = std::make_shared<Matrix>(*SparseStiffness);

    massMatrix *= user::rho * user::S * l;
    stiffnessMatrix *= user::E * user::S / l;
    massMatrix.makeCompressed();
    stiffnessMatrix.makeCompressed();

    double xsi = 1000.0;
    std::cout << "xsi:" << xsi * user::S / l << "\n";
    // dampingMatrix *= xsi * user::S / l;

    // -- Initial positions and velocities --
    siconos::algebra::SiconosVector q0{ndof};
    q0.setConstant(user::position_init);
    siconos::algebra::SiconosVector v0{ndof};
    v0.setConstant(user::velocity_init);

    Eigen::Map<siconos::algebra::SiconosSparseMatrix> mapMass{
        massMatrix.rows(),          massMatrix.cols(),          massMatrix.nonZeros(),
        massMatrix.outerIndexPtr(), massMatrix.innerIndexPtr(), massMatrix.valuePtr()};

    // -- The dynamical system --
    auto bar = std::make_shared<siconos::modeling::LagrangianSparseLinearTIDS>(
        q0, v0, mapMass, siconos::algebra::alias_t);

    // -- Set stiffness matrix (weight) --
    Eigen::Map<siconos::algebra::SiconosSparseMatrix> mapK{
        stiffnessMatrix.rows(),          stiffnessMatrix.cols(),
        stiffnessMatrix.nonZeros(),      stiffnessMatrix.outerIndexPtr(),
        stiffnessMatrix.innerIndexPtr(), stiffnessMatrix.valuePtr()};
    bar->setStiffnessMatrix(mapK, siconos::algebra::alias_t);
    //    bar->setDampingMatrix(dampingMatrix);

    // -- Set external forces (weight) --
    siconos::algebra::SiconosVector weight{ndof};
    weight.setZero();
    bar->setConstantFext(weight, siconos::algebra::alias_t);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.0;

    // Interaction bar-floor
    //
    auto H = std::make_shared<siconos::algebra::SiconosDenseMatrix>(1, ndof);
    H->setZero();  // (Eigen does not initialize)
    (*H)(0, 0) = 1.0;

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H);

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // -------------
    // --- Model ---
    // -------------
    auto impactingBar =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(user::t0, user::T);

    // add the dynamical system in the non smooth dynamical system
    impactingBar->insertDynamicalSystem(bar);

    // link the interaction and the dynamical system
    impactingBar->link(inter, bar);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --

#ifdef TS_VELOCITY_LEVEL
    auto OSI = std::make_shared<siconos::integrators::D1MinusLinearOSI>(
        siconos::integrators::D1MinusLinearOSI::Type::halfexplicit_velocity_level);

#else
    auto OSI = std::make_shared<siconos::integrators::D1MinusLinearOSI>();
#endif

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(user::t0, user::h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto force = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    auto s = std::make_shared<siconos::simulation::TimeSteppingD1Minus>(impactingBar, t, 2);
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
    s->insertNonSmoothProblem(force, siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = floor((user::T - user::t0) / user::h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 13;
    siconos::algebra::SiconosDenseMatrix dataPlot(N, outputSize);

    auto q = bar->q_read();
    auto v = bar->velocity_read();
    auto p = bar->p_read(1);
    auto lambda = inter->lambda(1);

    auto y = inter->y(0);
    int k = 0;
    dataPlot(k, 0) = impactingBar->t0();
    dataPlot(k, 1) = q(0);
    dataPlot(k, 2) = v(0);
    dataPlot(k, 3) = p(0);
    dataPlot(k, 4) = (*lambda)(0);
    dataPlot(k, 11) = 0.0; /* not yet initialized (*lambdaminus)(0); // lambda1_{k+1}^- */
    dataPlot(k, 12) = 0.0;

    dataPlot(k, 7) = q(ndof - 1);
    dataPlot(k, 8) = v(ndof - 1);
    dataPlot(k, 9) = q((ndof) / 2);
    dataPlot(k, 10) = v((ndof) / 2);

    siconos::algebra::SiconosVector tmp{ndof};
    tmp = stiffnessMatrix * q;
    double potentialEnergy = q.dot(tmp);
    tmp = massMatrix * v;
    double kineticEnergy = v.dot(tmp);

    dataPlot(k, 5) = potentialEnergy;
    dataPlot(k, 6) = kineticEnergy;

    //    std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
    //     std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;

    // --- Time loop ---
    cout << "====> Start computation ... \n\n";
    // ==== Simulation loop - Writing without explicit event handling =====

    auto start = std::chrono::system_clock::now();
    //    while (s->nextTime() < T)
    while ((s->hasNextEvent())) {
      s->advanceToEvent();
      //      std::cout << "k = "  << k << std::endl;
      //       std::cout << "position"  << std::endl;
      //       siconos::algebra::print(*q);
      //       std::cout << "velocity"  << std::endl;
      //       siconos::algebra::print(*v);

      // --- Get values to be plotted ---
      const auto& lambdaplus = inter->lambdaMemory(2).getSiconosVector(0);
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(0);
      dataPlot(k, 2) = v(0);
      dataPlot(k, 3) = p(0);
      dataPlot(k, 4) = (*lambda)(0);

      dataPlot(k, 11) = (*inter->lambda(2))(0);  // lambda1_{k+1}^-
      dataPlot(k, 12) = lambdaplus(0);
      ;

      dataPlot(k, 7) = q(ndof - 1);
      dataPlot(k, 8) = v(ndof - 1);
      dataPlot(k, 9) = q((ndof) / 2);
      dataPlot(k, 10) = v((ndof) / 2);

      tmp = stiffnessMatrix * q;
      potentialEnergy = q.dot(tmp);
      tmp = massMatrix * v;
      kineticEnergy = v.dot(tmp);
      dataPlot(k, 5) = potentialEnergy;
      dataPlot(k, 6) = kineticEnergy;

      //      std::cout << "q" << std::endl;
      //       siconos::algebra::print(*q);

      //       std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
      //       std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;

      s->processEvents();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("ImpactingBarD1MinusLinear.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-11;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "ImpactingBarD1MinusLinear.ref", eps)) >= eps)
      return 1;

    return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
