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

#include <KneeJointR.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

#include "GeomTools.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

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

    std::cout << "====> Model loading ...\n\n";

    // -- Initial positions and velocities --

    // First DS
    Vector q10{qDim};
    Vector v10{nDim};
    Matrix I1{3, 3};
    v10.setZero();
    I1.setIdentity();
    I1.setValue(0, 0, 0.1);
    // Initial position of the center of gravity CG1
    q10.setZero();
    q10(0) = 0.5 * L1 / sqrt(2.0);
    q10(2) = -0.5 * L1 / sqrt(2.0);
    // Initial orientation (a quaternion that gives the rotation w.r.t the spatial frame)
    // angle of the rotation Pi/4
    double angle = std::numbers::pi / 4;
    // vector of the rotation (Y-axis)
    Vector V1{3};
    V1 << 0., 1., 0.;
    // construction of the quaternion
    q10(3) = cos(angle * 0.5);
    q10(4) = V1(0) * sin(angle * 0.5);
    q10(5) = V1(1) * sin(angle * 0.5);
    q10(6) = V1(2) * sin(angle * 0.5);

    // -- The dynamical system --
    auto beam1 = std::make_shared<siconos::modeling::NewtonEulerDS>(q10, v10, m, I1,
                                                                    siconos::algebra::alias_t);
    // -- Set external forces (weight) --
    Vector weight{nDof};
    weight.setZero();
    weight(2) = -m * g;
    beam1->setConstantFext(weight, siconos::algebra::alias_t);

    // --------------------
    // --- Interactions ---
    // --------------------

    Vector P{3};
    P.setZero();
    // Building the first knee joint for beam1
    // input  - the concerned DS : beam1
    //        - a point in the spatial frame (absolute frame) where the knee is defined P
    auto relation1 = std::make_shared<siconos::joints::KneeJointR>(P, true, beam1);
    auto nslaw1 = std::make_shared<siconos::modeling::EqualityConditionNSL>(
        relation1->numberOfConstraints());

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);

    // -------------
    // --- Model ---
    // -------------
    auto myModel = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    myModel->insertDynamicalSystem(beam1);
    // link the interaction and the dynamical system
    myModel->link(inter1, beam1);
    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanCombinedProjectionOSI>(theta);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::MLCP>();
    auto osnspb_pos =
        std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>(
            SICONOS_MLCP_ENUM);

    // -- (4) Simulation setup with (1) (2) (3)

    auto s = std::make_shared<siconos::simulation::TimeSteppingCombinedProjection>(
        myModel, t, OSI, osnspb, osnspb_pos);
    s->setProjectionMaxIteration(1000);
    s->setConstraintTolUnilateral(1e-08);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 15 + 7;
    int N = 999;
    Matrix dataPlot(N, outputSize);
    Matrix beam1Plot(2, 3 * N);

    auto q1 = beam1->q_read();
    auto y = inter1->y(0);
    auto ydot = inter1->y(1);

    // --- Time loop ---
    std::cout << "====> Start computation ... \n\n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 0;

    auto start = std::chrono::system_clock::now();
    fprintf(pFile, "double T[%d*%d]={", N + 1, outputSize);
    std::vector<double> beamTipTrajectories(6);

    for (k = 0; k < N; k++) {
      // solve ...
      // s->newtonSolve(1e-4, 50);
      s->advanceToEvent();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();

      dataPlot(k, 1) = q1(0);
      dataPlot(k, 2) = q1(1);
      dataPlot(k, 3) = q1(2);
      dataPlot(k, 4) = q1(3);
      dataPlot(k, 5) = q1(4);
      dataPlot(k, 6) = q1(5);
      dataPlot(k, 7) = q1(6);
      dataPlot(k, 8) = y->norm();
      dataPlot(k, 9) = ydot->norm();
      // dataPlot(k, 8) = (*y)(0);
      // dataPlot(k, 9) = (*y)(1);
      // dataPlot(k, 10) = (*y)(2);
      // dataPlot(k, 11) = (*ydot)(0);
      // dataPlot(k, 12) = (*ydot)(1);
      // dataPlot(k, 13) = (*ydot)(2);
      // dataPlot(k, 14) = y->norm();
      // dataPlot(k, 15) = ydot->norm();

      geomtools::tipTrajectories(q1, beamTipTrajectories, L1);
      beam1Plot(0, 3 * k) = beamTipTrajectories[0];
      beam1Plot(0, 3 * k + 1) = beamTipTrajectories[1];
      beam1Plot(0, 3 * k + 2) = beamTipTrajectories[2];
      beam1Plot(1, 3 * k) = beamTipTrajectories[3];
      beam1Plot(1, 3 * k + 1) = beamTipTrajectories[4];
      beam1Plot(1, 3 * k + 2) = beamTipTrajectories[5];

      for (unsigned int jj = 0; jj < outputSize; jj++) {
        if ((k || jj)) fprintf(pFile, ",");
        fprintf(pFile, "%f", dataPlot(k, jj));
      }
      fprintf(pFile, "\n");
      s->processEvents();
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";
    fprintf(pFile, "};");
    fclose(pFile);

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("NE_1DS_1Knee_MLCP_MoreauJeanCombinedProjection.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("NE_1DS_1Knee_MLCP_MoreauJeanCombinedProjection_beam1.dat",
                                beam1Plot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "NE_1DS_1Knee_MLCP_MoreauJeanCombinedProjection.ref", eps)) >= eps)
      return 1;

    return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
