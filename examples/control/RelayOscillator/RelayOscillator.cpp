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
#include <string>
using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

// main program
/* Example of a limit cycle in Relay system predicted by the describing function approach
 * see Example 7.9 in "Nonlinear Systems" H. Khalil. Third Edition. Prentice Hall, 2002 ISBN
 * 0-13-067389-7
 *
 * The transfer function of the linear part is $G(s) 1/(s(s+1)(s+2))$ and the describing
 * function is $\Psi(a)=4/(\pi a)$ The predicted limit cycle is for the output $y(t) = a \sin
 * (\omega t)$ with $a = 2/(3\pi)$ and $\omega = \sqrt(2)$.
 */

int main(int argc, char* argv[]) {
  // Exception handling
  try {
    // == User-defined parameters ==
    unsigned int ndof = 3;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 100;     // Total simulation times
    double h = 1.0e-2;  // Time step
    double xinit = 3.0 * sqrt(2.0) / (2.0 * M_PI);

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
    (*A)(0, 1) = 1.0;
    (*A)(1, 2) = 1.0;
    (*A)(2, 1) = -3.0;
    (*A)(2, 2) = -2.0;
    auto x0 = std::make_shared<Vector>(ndof);
    x0->setZero();
    (*x0)(1) = xinit;

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
    B->setZero();
    (*B)(2, 0) = 1.0;

    auto C = std::make_shared<Matrix>(ninter, ndof);
    C->setZero();
    (*C)(0, 0) = 1.0;
    auto myProcessRelation = std::make_shared<siconos::modeling::FirstOrderLinearR>();
    myProcessRelation->setConstantB(*B);
    myProcessRelation->setConstantC(*C);
    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda

    // NonSmoothLaw
    unsigned int nslawSize = 1;
    auto myNslaw = std::make_shared<siconos::modeling::RelayNSL>(nslawSize);

    myNslaw->display();

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
    unsigned int outputSize = 7;          // number of required data
    unsigned int N = ceil((T - t0) / h);  // Number of time steps

    Matrix dataPlot{N + 1, outputSize};

    auto xProc = process->x();
    auto lambdaProc = myProcessInteraction->lambda(0);
    auto yProc = myProcessInteraction->y(0);
    unsigned int k = 0;  // Current step

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = relayOscillator->t0();  // Initial time of the model
    dataPlot(k, 1) = (*xProc)(0);
    dataPlot(k, 2) = (*xProc)(1);
    dataPlot(k, 3) = (*xProc)(2);
    dataPlot(k, 4) = (*lambdaProc)(0);
    dataPlot(k, 5) = (*yProc)(0);

    // ==== Simulation loop =====
    cout << "====> Start computation ... \n\n";

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    // Simulation loop
    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      k++;

      //  osnspb->setNumericsVerboseMode(1);

      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      dataPlot(k, 3) = (*xProc)(2);
      dataPlot(k, 4) = (*lambdaProc)(0);
      dataPlot(k, 5) = (*yProc)(0);
      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    cout << "====> Output file writing ...\n";

    siconos::algebra::io::write("RelayOscillator.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "RelayOscillator.ref", eps)) >
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
