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

#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <cmath>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

// main program

#include <chrono>

int main(int argc, char* argv[]) {
  // Exception handling
  try {
    // == User-defined parameters ==
    unsigned int ndof = 2;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 100;     // Total simulation times
    double h = 1.0e-2;  // Time step
    double xinit = 0.0;
    double vinit = 0.0;
    char filename[50] = "simu.";
    if (argc == 1) {
      xinit = 12.0;
      vinit = 6.0;
      strncpy(&filename[5], "1.0.1.0.log", 7);
    } else if (argc == 3) {
      // printf("argv[0] %s\n", argv[0]);
      printf("xinit is set to %f\n", atof(argv[1]));
      printf("vinit is set to %f\n", atof(argv[2]));

      xinit = atof(argv[1]);
      vinit = atof(argv[2]);
      int sizeofargv1 = strlen(argv[1]);
      // printf("sizeofargv1 %i\n",sizeofargv1);
      strncpy(&filename[5], argv[1], sizeofargv1);
      int sizeofargv2 = strlen(argv[2]);
      // printf("sizeofargv2 %i\n",sizeofargv2);
      strncpy(&filename[5 + sizeofargv1], ".", 1);

      strncpy(&filename[5 + sizeofargv1 + 1], argv[2], sizeofargv2);
      strncpy(&filename[5 + sizeofargv1 + sizeofargv2 + 1], ".log", 4);

      // printf("Output is written in filename %s\n",  filename);
    } else {
      cout << "wrong  number of arguments = " << argc << endl;
    }

    double m = 1, stiffness = 1;
    double alpha = 1.0;

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
    (*A)(0, 0) = 0.0;
    (*A)(0, 1) = 1.0;
    (*A)(1, 0) = -stiffness / m;
    (*A)(1, 1) = 0.0;
    auto x0 = std::make_shared<Vector>(ndof);
    (*x0)(0) = xinit;
    (*x0)(1) = vinit;

    auto process = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        *x0, siconos::algebra::alias_t);
    process->setConstantA(*A, siconos::algebra::alias_t);

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 1;  // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    auto B = std::make_shared<Matrix>(ndof, ninter);
    (*B)(0, 0) = 0.0;
    (*B)(1, 0) = alpha;

    auto C = std::make_shared<Matrix>(ninter, ndof);
    (*C)(0, 0) = 0.0;
    (*C)(0, 1) = 1.0;

    auto myProcessRelation = std::make_shared<siconos::modeling::FirstOrderLinearR>();
    myProcessRelation->setConstantC(*C);
    myProcessRelation->setConstantB(*B);

    // NonSmoothLaw
    unsigned int nslawSize = 1;
    auto myNslaw = std::make_shared<siconos::modeling::RelayNSL>(nslawSize);

    // The Interaction which involves the first DS (the process)
    auto myProcessInteraction =
        std::make_shared<siconos::modeling::Interaction>(myNslaw, myProcessRelation);

    // -------------
    // --- Model ---
    // -------------
    auto relayOscillator =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    relayOscillator->insertDynamicalSystem(process);
    relayOscillator->link(myProcessInteraction, process);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // == Creation of the Simulation ==
    auto s = std::make_shared<siconos::simulation::TimeStepping>(relayOscillator, td);
    // -- OneStepIntegrators --
    double theta = 0.5;
    auto myIntegrator = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);
    s->insertIntegrator(myIntegrator);

    // -- OneStepNsProblem --

    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::Relay>();

    osnspb->setSolverId(SICONOS_RELAY_LEMKE);
    osnspb->numericsSolverOptions()->dparam[0] = 1e-08;
    s->insertNonSmoothProblem(osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    unsigned int outputSize = 10;         // number of required data
    unsigned int N = ceil((T - t0) / h);  // Number of time steps

    Matrix dataPlot(N, outputSize);

    auto xProc = process->x();
    auto lambdaProc = myProcessInteraction->lambda(0);
    auto yProc = myProcessInteraction->y(0);
    auto vectorfield = process->rhs();
    unsigned int k = 0;  // Current step

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = relayOscillator->t0();  // Initial time of the model
    dataPlot(k, 1) = (*xProc)(0);
    dataPlot(k, 2) = (*xProc)(1);
    dataPlot(k, 3) = (*lambdaProc)(0);
    dataPlot(k, 4) = (*yProc)(0);
    dataPlot(k, 7) = (*vectorfield)(0);
    dataPlot(k, 8) = (*vectorfield)(1);

    // ==== Simulation loop =====
    cout << "====> Start computation ... \n\n";

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    // Simulation loop
    auto start = std::chrono::system_clock::now();
    while (k < N - 1) {
      k++;

      //  osnspb->setNumericsVerboseMode(1);

      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      dataPlot(k, 3) = (*lambdaProc)(0);
      dataPlot(k, 4) = (*yProc)(0);
      process->computeRhs(s->nextTime());
      if (k == 1)  // tricks just for display to avoid the computation of the initial Rhs
      {
        dataPlot(k - 1, 7) = (*vectorfield)(0);
        dataPlot(k - 1, 8) = (*vectorfield)(1);
      }

      dataPlot(k, 7) = (*vectorfield)(0);
      dataPlot(k, 8) = (*vectorfield)(1);
      s->nextStep();
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("FrictionOscillator.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write(filename, dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "FrictionOscillator.ref",
                                                      eps)) >= eps)
      return 1;

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
