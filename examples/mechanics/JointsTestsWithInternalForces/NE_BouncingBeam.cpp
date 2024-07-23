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

/*!\file NE....cpp
  \brief \ref EMNE_MULTIBODY - C++ input file, Time-Stepping version - O.B.

  A multibody example.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include <KneeJointR.hpp>
#include <PrismaticJointR.hpp>
#include <SiconosKernel.hpp>
#include <chrono>

#include "GeomTools.h"

using namespace std;

int main(int argc, char *argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;
    unsigned int qDim = 7;
    unsigned int nDim = 6;
    double t0 = 0;    // initial computation time
    double T = 10.0;  // final computation time
    double h = 0.01;  // time step
    int N = 1000;
    double L1 = 1.0;
    double L2 = 1.0;
    double L3 = 1.0;
    double theta = 1.0;  // theta for MoreauJeanOSI integrator
    double g = 9.81;     // Gravity
    double m = 1.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    FILE *pFile;
    pFile = fopen("data.h", "w");
    if (pFile == NULL) {
      printf("fopen exampleopen filed!\n");
      fclose(pFile);
    }

    cout << "====> Model loading ..." << endl << endl;
    // -- Initial positions and velocities --
    SP::SiconosVector q03(new SiconosVector(qDim));
    SP::SiconosVector v03(new SiconosVector(nDim));
    SP::SimpleMatrix I3(new SimpleMatrix(3, 3));
    v03->zero();
    I3->eye();
    I3->setValue(0, 0, 0.1);
    q03->zero();
    (*q03)(2) = -L1 * sqrt(2.0) - L1 / 2;

    double angle = M_PI / 2;
    SiconosVector V1(3);
    V1.zero();
    V1.setValue(0, 0);
    V1.setValue(1, 1);
    V1.setValue(2, 0);
    q03->setValue(3, cos(angle / 2));
    q03->setValue(4, V1.getValue(0) * sin(angle / 2));
    q03->setValue(5, V1.getValue(1) * sin(angle / 2));
    q03->setValue(6, V1.getValue(2) * sin(angle / 2));

    SP::NewtonEulerDS bouncingbeam(new NewtonEulerDS(q03, v03, m, I3));
    // -- Set external forces (weight) --
    SP::SiconosVector weight3(new SiconosVector(nDof));
    (*weight3)(2) = -m * g;
    bouncingbeam->setFExtPtr(weight3);
    bouncingbeam->setComputeFIntFunction("SimplePlugin", "fInt_beam1");
    bouncingbeam->setComputeJacobianFIntqFunction("SimplePlugin", "jacobianFIntq_beam1");
    bouncingbeam->setComputeJacobianFIntvFunction("SimplePlugin", "jacobianFIntv_beam1");

    // --------------------
    // --- Interactions ---
    // --------------------

    // Interaction with the floor
    double e = 0.9;
    SP::SimpleMatrix H(new SimpleMatrix(1, qDim));
    SP::SiconosVector eR(new SiconosVector(1));
    eR->setValue(0, 2.3);
    H->zero();
    (*H)(0, 2) = 1.0;
    SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
    SP::NewtonEulerR relation0(new NewtonEulerR());
    relation0->setJachq(H);
    relation0->setE(eR);
    cout << "main jacQH" << endl;
    relation0->jachq()->display();

    SP::SiconosVector axe1(new SiconosVector(3));
    axe1->zero();
    axe1->setValue(2, 1);

    SP::PrismaticJointR relation4(new PrismaticJointR(axe1, false, bouncingbeam));
    SP::NonSmoothLaw nslaw4(new EqualityConditionNSL(relation4->numberOfConstraints()));
    SP::Interaction inter4(new Interaction(nslaw4, relation4));
    SP::Interaction interFloor(new Interaction(nslaw0, relation0));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem myModel(new NonSmoothDynamicalSystem(t0, T));
    // add the dynamical system in the non smooth dynamical system
    myModel->insertDynamicalSystem(bouncingbeam);
    // link the interaction and the dynamical system
    myModel->link(inter4, bouncingbeam);
    myModel->link(interFloor, bouncingbeam);
    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --

    SP::MoreauJeanOSI OSI3(new MoreauJeanOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new MLCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(myModel, t, OSI3, osnspb));
    s->setNewtonTolerance(1e-4);
    s->setNewtonMaxIteration(50);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 15 + 7;
    SimpleMatrix dataPlot(N, outputSize);
    SimpleMatrix bouncingbeamPlot(2, 3 * N);

    SP::SiconosVector q3 = bouncingbeam->q();
    SP::SiconosVector y = interFloor->y(0);
    SP::SiconosVector ydot = interFloor->y(1);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 0;

    auto start = std::chrono::system_clock::now();
    fprintf(pFile, "double T[%d*%d]={", N + 1, outputSize);
    double beamTipTrajectories[6];

    for (k = 0; k < N; k++) {
      // solve ...
      s->advanceToEvent();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();

      dataPlot(k, 1) = (*q3)(0);
      dataPlot(k, 2) = (*q3)(1);
      dataPlot(k, 3) = (*q3)(2);
      dataPlot(k, 4) = (*q3)(3);
      dataPlot(k, 5) = (*q3)(4);
      dataPlot(k, 6) = (*q3)(5);
      dataPlot(k, 7) = (*q3)(6);

      dataPlot(k, 8) = y->norm2();
      dataPlot(k, 9) = ydot->norm2();

      geomtools::tipTrajectories(q3, beamTipTrajectories, L3);
      bouncingbeamPlot(0, 3 * k) = beamTipTrajectories[0];
      bouncingbeamPlot(0, 3 * k + 1) = beamTipTrajectories[1];
      bouncingbeamPlot(0, 3 * k + 2) = beamTipTrajectories[2];
      bouncingbeamPlot(1, 3 * k) = beamTipTrajectories[3];
      bouncingbeamPlot(1, 3 * k + 1) = beamTipTrajectories[4];
      bouncingbeamPlot(1, 3 * k + 2) = beamTipTrajectories[5];

      // printf("reaction1:%lf \n", interFloor->lambda(1)->getValue(0));

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
    dataPlot.resize(k, outputSize);
    ioMatrix::write("NE_BouncingBeam.dat", "ascii", dataPlot, "noDim");
    ioMatrix::write("NE_BouncingBeam_beam.dat", "ascii", bouncingbeamPlot, "noDim");

    double error = 0.0, eps = 1e-12;
    if ((error = ioMatrix::compareRefFile(dataPlot, "NE_BouncingBeam.ref", eps)) >= 0.0 &&
        error > eps)
      return 1;

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
