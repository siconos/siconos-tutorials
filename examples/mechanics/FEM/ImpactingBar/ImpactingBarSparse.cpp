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

using namespace std;
namespace user =
    user_defined_ref;  // To choose the set of parameters used in the current simulation

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

    massMatrix *= user::rho * user::S * l;
    stiffnessMatrix *= user::E * user::S / l;

    massMatrix.makeCompressed();
    stiffnessMatrix.makeCompressed();

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
#ifdef TS_PROJ
    auto OSI =
        std::make_shared<siconos::integrators::MoreauJeanDirectProjectionOSI>(user::theta);
    OSI->setDeactivateYPosThreshold(1e-05);
    OSI->setDeactivateYVelThreshold(0.0);
    OSI->setActivateYPosThreshold(1e-09);
    OSI->setActivateYVelThreshold(100.0);
#else
#ifdef TS_COMBINED
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanCombinedProjectionOSI>(theta);
#else
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(user::theta, 0.0);
    // OSI->setConstraintActivationThreshold(1e-05);
#endif
#endif
    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(user::t0, user::h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
#ifdef TS_PROJ
    auto position =
        std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>();
    auto s = std::make_shared<siconos::simulation::TimeSteppingDirectProjection>(
        impactingBar, t, OSI, osnspb, position, 0);
    s->setProjectionMaxIteration(10);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);

#else
#ifdef TS_COMBINED
    auto position =
        std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>(
            SICONOS_MLCP_ENUM);
    auto s = std::make_shared<siconos::simulation::TimeSteppingCombinedProjection>(
        impactingBar, t, OSI, osnspb, position, 2);
    s->setProjectionMaxIteration(500);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);
#else
    auto s = std::make_shared<siconos::simulation::TimeStepping>(impactingBar, t, OSI, osnspb);
#endif
#endif

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    int N = floor((user::T - user::t0) / user::h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 12;
    siconos::algebra::SiconosDenseMatrix dataPlot(N, outputSize);

    auto q = bar->q_read();
    auto v = bar->velocity_read();
    auto p = bar->p_read(1);
    auto lambda = inter->lambda(1);

    dataPlot(0, 0) = impactingBar->t0();
    dataPlot(0, 1) = q(0);
    dataPlot(0, 2) = v(0);
    dataPlot(0, 3) = p(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 7) = q(ndof - 1);
    dataPlot(0, 8) = v(ndof - 1);
    dataPlot(0, 9) = q((ndof) / 2);
    dataPlot(0, 10) = v((ndof) / 2);

    siconos::algebra::SiconosVector tmp{ndof};
    tmp = stiffnessMatrix * q;
    double potentialEnergy = q.dot(tmp);
    tmp = massMatrix * v;
    double kineticEnergy = v.dot(tmp);
    double impactEnergy = 0.0;

    dataPlot(0, 5) = potentialEnergy;
    dataPlot(0, 6) = kineticEnergy;
    dataPlot(0, 11) = impactEnergy;

    //    std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
    //    std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;

    // --- Time loop ---
    cout << "====> Start computation ... \n\n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();
    //    while (s->nextTime() < T)
    while (k < N) {
      s->advanceToEvent();

      //       std::cout << "position"  << std::endl;
      //       siconos::algebra::print(*q);
      //       std::cout << "velocity"  << std::endl;
      //       siconos::algebra::print(*v);

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(0);
      dataPlot(k, 2) = v(0);
      dataPlot(k, 3) = p(0) / user::h;
      dataPlot(k, 4) = (*lambda)(0);
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

      double v_p_theta = user::theta * v(0) + (1 - user::theta) * dataPlot(k - 1, 2);
      impactEnergy = v_p_theta * (*lambda)(0);
      dataPlot(k, 11) = impactEnergy;

      //      std::cout << "q" << std::endl;
      //       siconos::algebra::print(*q);

      //       std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
      //       std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;

      s->processEvents();

      k++;
    }
    double totalImpactEnergy = 0.0;
    for (int p = 0; p < k; p++) {
      totalImpactEnergy = totalImpactEnergy + dataPlot(p, 11);
    }
    printf("totalImpactEnergy = %e", totalImpactEnergy);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("ImpactingBar.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-9;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "ImpactingBar.ref", eps)) >=
        eps)
      return 1;

    return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
