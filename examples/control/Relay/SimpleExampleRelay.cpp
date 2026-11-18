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
#include <SiconosKernel.hpp>
#include <chrono>

using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

// main program
int main(int argc, char* argv[]) {
  // Exception handling
  try {
    // == User-defined parameters ==
    unsigned int ndof = 2;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 0.8;       // Total simulation time
    double h = 1.0e-3;  // Time step
    double Vinit = 1.0;

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

    auto x0 = std::make_shared<Vector>(ndof);
    (*x0)(0) = Vinit;
    (*x0)(1) = -Vinit;
    auto process = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        *x0, siconos::algebra::SiconosDenseMatrix::Zero(x0->size(), x0->size()),
        siconos::algebra::SiconosVector::Zero(x0->size()), siconos::algebra::copy_t);

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 2;  // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    auto B = std::make_shared<Matrix>(ndof, ninter);
    B->setZero();
    (*B)(0, 0) = 2.0;
    (*B)(1, 1) = 2.0;
    auto C = std::make_shared<Matrix>(ninter, ndof);
    C->setZero();
    (*C)(0, 0) = 1.0;
    (*C)(1, 1) = 1.0;
    auto myProcessRelation = std::make_shared<siconos::modeling::FirstOrderLinearR>();
    myProcessRelation->setConstantB(*B);
    myProcessRelation->setConstantC(*C);
    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda

    // NonSmoothLaw
    unsigned int nslawSize = 2;
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
    auto osnspb1 = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::Relay>(SICONOS_RELAY_PGS);
    s->insertNonSmoothProblem(osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    unsigned int outputSize = 7;          // number of required data
    unsigned int N = ceil((T - t0) / h) + 1;  // Number of time steps

    Matrix dataPlot(N, outputSize);

    auto xProc = process->x();
    auto lambdaProc = myProcessInteraction->lambda(0);
    auto yProc = myProcessInteraction->y(0);

    myProcessInteraction->computeOutput(t0, 0);

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = simpleExampleRelay->t0();  // Initial time of the model
    dataPlot(0, 1) = (*xProc)(0);
    dataPlot(0, 2) = (*xProc)(1);
    dataPlot(0, 3) = (*lambdaProc)(0);
    dataPlot(0, 4) = (*lambdaProc)(1);
    dataPlot(0, 5) = (*yProc)(0);
    dataPlot(0, 6) = (*yProc)(1);

    // ==== Simulation loop =====
    cout << "====> Start computation ... \n\n";

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    int k = 0;  // Current step

    // Simulation loop
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      k++;
      // cout << "step --> " << k << endl;
      //     osnspb->setNumericsVerboseMode(1);
      //   *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      dataPlot(k, 3) = (*lambdaProc)(0);
      dataPlot(k, 4) = (*lambdaProc)(1);
      dataPlot(k, 5) = (*yProc)(0);
      dataPlot(k, 6) = (*yProc)(1);
      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time \n";
    ;
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("SimpleExampleRelay.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "SimpleExampleRelay.ref",
                                                      eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
