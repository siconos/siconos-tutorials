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
#include <PrismaticJointR.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

#include "GeomTools.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDof = 3;
    unsigned int qDim = 7;
    unsigned int nDim = 6;
    double t0 = 0;    // initial computation time
    double T = 10.0;  // final computation time
    double h = 0.01;  // time step
    double L1 = 1.0;
    double L2 = 1.0;
    double L3 = 1.0;
    double theta = 1.0;  // theta for MoreauJeanOSI integrator
    double g = 9.81;     // Gravity
    double m = 1.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    FILE* pFile;
    pFile = fopen("data.h", "w");
    if (pFile == NULL) {
      printf("fopen exampleopen filed!\n");
      fclose(pFile);
    }

    std::cout << "====> Model loading ...\n";

    // -- Initial positions and velocities --
    Vector q03{qDim};
    Vector v03{nDim};
    Matrix I3{3, 3};
    v03.setZero();
    q03.setZero();
    I3.setIdentity();
    I3(0, 0) = 0.1;
    q03(2) = -L1 * sqrt(2.0) - L1 / 2;

    double angle = std::numbers::pi / 2;
    Vector V1{3};
    V1 << 0, 1, 0;
    q03(3) = cos(angle / 2);
    q03(4) = V1(0) * sin(angle / 2);
    q03(5) = V1(1) * sin(angle / 2);
    q03(6) = V1(2) * sin(angle / 2);

    auto bouncingbeam = std::make_shared<siconos::modeling::NewtonEulerDS>(
        q03, v03, m, I3, siconos::algebra::alias_t);
    // -- Set external forces (weight) --
    Vector weight{nDof};
    weight.setZero();
    weight(2) = -m * g;
    bouncingbeam->setConstantFext(weight, siconos::algebra::alias_t);
    bouncingbeam->setComputeFintFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& twist,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
           Eigen::Ref<siconos::algebra::MapVectorType> fint) {
          auto i = 0;
          fint(0) = 1e4 * q(0);
          fint(1) = 0.0;
          fint(2) = 1e4 * q(2);
        });

    bouncingbeam->setComputeJacobianFintOver_qFunction(
        [](const Eigen::Ref<const siconos::algebra::SiconosVector>& twist,
           const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();

          result(0, 0) = 1e4;
          result(2, 2) = 1e4;
        });

    // --------------------
    // --- Interactions ---
    // --------------------

    // Interaction with the floor
    double e = 0.9;
    Matrix H{1, qDim};
    Vector eR{1};
    eR << 2.3;
    H.setZero();
    H(0, 2) = 1.0;
    auto nslaw0 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation0 = std::make_shared<siconos::modeling::NewtonEulerR>();
    relation0->setConstantConstraintJacobian(H);
    relation0->setConstanteVector(eR);

    // Interactions
    // Building the prismatic joint for bouncingbeam
    // input  - the first concerned DS : bouncingbeam
    //        - an axis in the spatial frame (absolute frame)
    // auto H4= std::make_shared<Matrix>(PrismaticJointR::numberOfConstraints(), qDim);
    // H4->setZero();

    Vector axe1{3};
    axe1 << 0., 0, 1.;
    auto relation4 =
        std::make_shared<siconos::joints::PrismaticJointR>(axe1, false, bouncingbeam);
    auto nslaw4 = std::make_shared<siconos::modeling::EqualityConditionNSL>(
        relation4->numberOfConstraints());

    auto inter4 = std::make_shared<siconos::modeling::Interaction>(nslaw4, relation4);
    auto interFloor = std::make_shared<siconos::modeling::Interaction>(nslaw0, relation0);

    // -------------
    // --- Model ---
    // -------------
    auto myModel = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    myModel->insertDynamicalSystem(bouncingbeam);
    // link the interaction and the dynamical system

    myModel->link(inter4, bouncingbeam);
    myModel->link(interFloor, bouncingbeam);
    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --

    auto OSI3 = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::MLCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(myModel, t, OSI3, osnspb);
    s->setNewtonTolerance(1e-4);
    s->setNewtonMaxIteration(50);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 15 + 7;
    int N = 1000;
    Matrix dataPlot(N, outputSize);
    Matrix bouncingbeamPlot(2, 3 * N);

    auto q3 = bouncingbeam->q_read();
    auto y = interFloor->y(0);
    auto ydot = interFloor->y(1);

    // --- Time loop ---
    std::cout << "====> Start computation ... \n\n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 0;

    auto start = std::chrono::system_clock::now();
    fprintf(pFile, "double T[%d*%d]={", N + 1, outputSize);
    std::vector<double> beamTipTrajectories(6);

    for (k = 0; k < N; k++) {
      // solve ...
      s->advanceToEvent();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();

      dataPlot(k, 1) = q3(0);
      dataPlot(k, 2) = q3(1);
      dataPlot(k, 3) = q3(2);
      dataPlot(k, 4) = q3(3);
      dataPlot(k, 5) = q3(4);
      dataPlot(k, 6) = q3(5);
      dataPlot(k, 7) = q3(6);

      dataPlot(k, 8) = y->norm();
      dataPlot(k, 9) = ydot->norm();

      geomtools::tipTrajectories(q3, beamTipTrajectories, L3);
      bouncingbeamPlot(0, 3 * k) = beamTipTrajectories[0];
      bouncingbeamPlot(0, 3 * k + 1) = beamTipTrajectories[1];
      bouncingbeamPlot(0, 3 * k + 2) = beamTipTrajectories[2];
      bouncingbeamPlot(1, 3 * k) = beamTipTrajectories[3];
      bouncingbeamPlot(1, 3 * k + 1) = beamTipTrajectories[4];
      bouncingbeamPlot(1, 3 * k + 2) = beamTipTrajectories[5];

      // printf("reaction1:%lf \n", (*interFloor->lambda(1))(0));

      for (unsigned int jj = 0; jj < outputSize; jj++) {
        if ((k || jj)) fprintf(pFile, ",");
        fprintf(pFile, "%f", dataPlot(k, jj));
      }
      fprintf(pFile, "\n");
      s->nextStep();
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    fprintf(pFile, "};");
    fclose(pFile);
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("NE_BouncingBeam.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("NE_BouncingBeam_beam.dat", bouncingbeamPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "NE_BouncingBeam.ref", eps)) >=
        eps)
      return 1;

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
