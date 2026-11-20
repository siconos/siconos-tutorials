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

/*!\file LuenbergerObserver.cpp
  \brief Academic example of a Luenberger Observer
  O. Huber.

  The controlled plant is a double integrator
  */

#include <math.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>
#include <string>

using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

constexpr unsigned int ndof = 4;
/* 02/2010 --> 08/2010*/
/*Author: Yamen ABDENNADHER */
/*Exemple from : Rafael Ramos, Dominigo Biel, Enric Fossas and Franisco Guinjoan. Interleaving
 * Quasi-Sliding-Mode Control of Parallel-Connected Buck-based Inverters. IEEE vol 55 N11
 * Novembre 2008. For more detail, please see Yamen's report.*/
// main program
int main(int argc, char* argv[]) {
  // Exception handling
  try {
    // == User-defined parameters ==
    // unsigned int ndof = 4;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 0.02;    // Total simulation time
    double h = 1.0e-6;  // Time step
    // double Vinit= 1.0;
    double Cp = 60e-6;
    double leq = 0.385e-3;
    double Rl = 5;
    double m = 3;

    double L1 = 1.5e-3;
    double L2 = 1.22e-3;
    double L3 = 0.9e-3;
    // leq=1/((1/L1)+(1/L2)+(1/L3));

    double rl1 = 94e-3;
    double rl2 = 116e-3;
    double rl3 = 100e-3;

    double E1 = 70;
    double E2 = 70;
    double E3 = 70;

    double k1 = 1;
    double k2 = 6e-5;

    // ================= Creation of the model =======================
    // Steps:
    // - create some Dynamical Systems
    // - create some Interactions between those Dynamical Systems
    //   Interaction = some relations (constraints) and a NonSmoothLaw
    // - create a NonSmoothDynamicalSystem with the DynamicalSystems and the Interactions
    // - add this NonSmoothDynamicalSystem into a Model
    // - add a Simulation to the model
    //  Simulation = TimeDiscretisation + OneStepIntegrator and OneStepNSProblem

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // First System:
    // dx/dt = Ax + u(t) + r
    // x(0) = x0
    // Note: r = Blambda, B defines in relation below.

    auto A = std::make_shared<Matrix>(ndof, ndof);
    A->setZero();
    (*A)(0, 1) = 1;
    (*A)(1, 0) = (-1 / (Cp * leq)) - (rl1 / (Cp * L1 * Rl));
    (*A)(1, 1) = (-rl1 / L1) - (1 / (Cp * Rl));
    (*A)(1, 2) = (rl1 / (Cp * L1)) - (rl2 / (Cp * L2));
    (*A)(1, 3) = (rl1 / (Cp * L1)) - (rl3 / (Cp * L3));
    (*A)(2, 0) = -1 / L2;
    (*A)(2, 2) = -rl2 / L2;
    (*A)(3, 0) = -1 / L3;
    (*A)(3, 3) = -rl3 / L3;

    auto x0 = std::make_shared<Vector>(ndof);
    (*x0)(0) = 50;
    (*x0)(1) = 7;
    (*x0)(2) = 4;
    (*x0)(3) = 4;

    auto process = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        *x0, siconos::algebra::alias_t);
    process->setConstantA(*A, siconos::algebra::alias_t);
    double zProc = 0.;

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 3;  // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda +eDLS
    // r = Blambda
    auto B = std::make_shared<Matrix>(ndof, ninter);
    B->setZero();
    (*B)(1, 0) = E1 / (Cp * L1);
    (*B)(1, 1) = E2 / (Cp * L2);
    (*B)(1, 2) = E3 / (Cp * L3);
    (*B)(2, 1) = E2 / L2;
    (*B)(3, 2) = E3 / L3;

    auto C = std::make_shared<Matrix>(ninter, ndof);
    C->setZero();
    (*C)(0, 0) = (Cp * k1 / k2) + ((m - 1) / Rl);
    (*C)(0, 1) = m * Cp;
    (*C)(0, 2) = -m;
    (*C)(0, 3) = -m;
    (*C)(1, 0) = (Cp * k1 / k2) - (1 / Rl);
    (*C)(1, 2) = m;
    (*C)(2, 0) = (Cp * k1 / k2) - (1 / Rl);
    (*C)(2, 3) = m;
    auto myProcessRelation = std::make_shared<siconos::modeling::FirstOrderLinearR>();
    myProcessRelation->setConstantC(*C);
    myProcessRelation->setConstantB(*B);

    myProcessRelation->setComputeeVectorFunction(
        [&zProc](double t, Eigen::Ref<siconos::algebra::MapVectorType> result) {
          zProc = 55 * sin(100 * std::numbers::pi * t);
          result.setConstant(-zProc);
        });

    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda

    // NonSmoothLaw
    unsigned int nslawSize = 3;
    auto myNslaw = std::make_shared<siconos::modeling::RelayNSL>(nslawSize);

    myNslaw->display();

    // The Interaction which involves the first DS (the process)
    auto myProcessInteraction =
        std::make_shared<siconos::modeling::Interaction>(myNslaw, myProcessRelation);

    // -------------
    // --- Model ---
    // -------------
    auto simpleExampleRelay =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    simpleExampleRelay->insertDynamicalSystem(process);
    simpleExampleRelay->link(myProcessInteraction, process);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // == Creation of the Simulation ==
    auto s = std::make_shared<siconos::simulation::TimeStepping>(simpleExampleRelay, td);
    // -- OneStepIntegrators --
    double theta = 0.5;
    auto myIntegrator = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);
    s->insertIntegrator(myIntegrator);

    // -- OneStepNsProblem --
    // -- OneStepNsProblem --

    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::Relay>(SICONOS_RELAY_ENUM);

    // osnspb->setNumericsSolverName("Lemke");
    // osnspb->numericsSolverOptions()->dparam[0]=1e-08;
    s->insertNonSmoothProblem(osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ...\n\n";

    // --- Get the values to be plotted ---
    unsigned int outputSize = 17;         // number of required data
    unsigned int N = ceil((T - t0) / h);  // Number of time steps

    Matrix dataPlot(N, outputSize);

    auto xProc = process->x();
    auto lambdaProc = myProcessInteraction->lambda(0);
    auto yProc = myProcessInteraction->y(0);
    myProcessInteraction->computeOutput(t0, 0);
    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = simpleExampleRelay->t0();  // Initial time of the model
    dataPlot(0, 1) = (*xProc)(0);
    dataPlot(0, 2) = (*xProc)(1);
    dataPlot(0, 3) = (*xProc)(2);
    dataPlot(0, 4) = (*xProc)(3);
    dataPlot(0, 5) = (*lambdaProc)(0);
    dataPlot(0, 6) = (*lambdaProc)(1);
    dataPlot(0, 7) = (*lambdaProc)(2);
    dataPlot(0, 8) = (*yProc)(0);
    dataPlot(0, 9) = (*yProc)(1);
    dataPlot(0, 10) = (*yProc)(2);
    dataPlot(0, 11) = zProc;       // v_ref
    dataPlot(0, 12) = zProc / Rl;  // i_ref
    dataPlot(0, 13) = 0;           // voltage_error
    dataPlot(0, 14) = 0;           // current_error
    dataPlot(0, 15) = 0;           // il1
    dataPlot(0, 16) = (*xProc)(1);

    // ==== Simulation loop =====
    cout << "====> Start computation ... \n\n";

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    unsigned int k = 0;  // Current step

    // Simulation loop
    auto start = std::chrono::system_clock::now();

    unsigned int i = 0;
    int j = 0;
    auto err = std::make_shared<Vector>(2);
    auto il1 = std::make_shared<Vector>(1);
    auto temps = std::make_shared<Vector>(N + 1);

    while (k < N - 1) {
      k++;
      //    osnspb->setNumericsVerboseMode(1);
      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      (*err)(0) = abs(zProc - (*xProc)(0));                // voltage error
      (*err)(1) = abs((zProc / Rl) - ((*xProc)(0) / Rl));  // current error
      (*il1)(0) = ((*xProc)(0) / Rl) - ((*xProc)(2) + (*xProc)(3)) +
                  Cp * (*xProc)(1);  // current through L_1

      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);       // output voltage v0
      dataPlot(k, 2) = (*xProc)(0) / Rl;  // output current i0
      dataPlot(k, 3) = (*xProc)(2);       // il2
      dataPlot(k, 4) = (*xProc)(3);       // il3
      dataPlot(k, 5) = (*lambdaProc)(0);  // u1
      dataPlot(k, 6) = (*lambdaProc)(1);  // u2
      dataPlot(k, 7) = (*lambdaProc)(2);  // u3
      dataPlot(k, 8) = (*yProc)(0);       // y1
      dataPlot(k, 9) = (*yProc)(1);       // y2
      dataPlot(k, 10) = (*yProc)(2);      // y3
      dataPlot(k, 11) = zProc;            // v_ref
      dataPlot(k, 12) = zProc / Rl;       // i_ref
      dataPlot(k, 13) = (*err)(0);        // voltage_error
      dataPlot(k, 14) = (*err)(1);        // current_error
      dataPlot(k, 15) = (*il1)(0);        // il1
      dataPlot(k, 16) = (*xProc)(1);      // v'0
      s->nextStep();

      //////////////////////////////////////////////////////////////////////////////////
      if (((*yProc)(0) > 1e-8) || ((*yProc)(0) < -1e-8)) {
        (*temps)(k) = (*yProc)(0);
      }
      ///////////////////////////////////////////////////////////////////////////////////
    }

    while (i < N - 1) {
      if ((*temps)(i) == 0) j = j + 1;

      i++;
    }
    cout << "The sliding mode appears at the step number : \n" << N - 1 - j << endl;

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time \n";
    ;
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("ParallelInverter.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-08;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "ParallelInverter.ref", eps)) >
        eps)
      return 1;
    else
      return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
