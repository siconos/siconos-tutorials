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

/*  O. Huber.

  The controlled plant is a double integrator
  */

#include <SiconosKernel.hpp>
#include <chrono>
#include <string>

using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;
using namespace std;

namespace user {
double computeControl(double time) {
  // double u;
  // double alpha = 1.0;
  // double T = 0.1;
  // int oddoreven = int(time / T);

  // if ((oddoreven / 2) == 0)
  //   u = alpha;
  // else
  //   u = 0;
  return 30. * sin(50. * time);
  ;
}
}  // namespace user

int main(int argc, char* argv[]) {
  // Exception handling
  try {
    // == User-defined parameters ==
    unsigned int ndof = 2;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 2;       // Total simulation time
    double h = 1.0e-4;  // Time step
    double Vinit = 10.0;

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
    (*A)(0, 0) = 1.0;
    (*A)(0, 1) = 1.0;
    (*A)(1, 0) = 3.0;
    (*A)(1, 1) = 1.0;
    auto x0 = std::make_shared<Vector>(ndof);
    x0->setZero();
    (*x0)(0) = Vinit;
    auto process = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        *x0, siconos::algebra::alias_t);
    process->setConstantA(*A, siconos::algebra::alias_t);
    process->setComputebVectorFunction(
        [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
          double u = user::computeControl(time);
          result(0) = u;
          result(1) = 2.0 * u;
        });

    // Second System, the observer:
    // dx/dt = A hatx + u(t) + L(y-haty)
    // hatx(0) = x0
    // y = Gx
    // u(t) + Ly  is computed with uObserver function

    unsigned int noutput = 1;
    auto L = std::make_shared<Matrix>(ndof, noutput);
    L->setConstant(1.);
    auto G = std::make_shared<Matrix>(noutput, ndof);
    G->setConstant(2.);

    // hatA is initialized with A
    auto hatA = std::make_shared<Matrix>(ndof, ndof);
    (*hatA)(0, 0) = -1.0;
    (*hatA)(0, 1) = -1.0;
    (*hatA)(1, 0) = 1.0;
    (*hatA)(1, 1) = -1.0;

    auto obsX0 = std::make_shared<Vector>(ndof);
    obsX0->setZero();
    auto param = process->x_read();
    auto observer = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        *obsX0, siconos::algebra::alias_t);
    observer->setConstantA(*hatA, siconos::algebra::alias_t);
    observer->setComputebVectorFunction(
        [&param](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
          double u = user::computeControl(time);
          double coeff = 2 * (param(0) + param(1));
          result(0) = u + coeff;
          result(1) = 2.0 * u + coeff;
        });

    // observer->setzPtr(process->x());
    //  The set of all DynamicalSystems
    //  --------------------
    //  --- Interactions ---
    //  --------------------
    unsigned int ninter = 1;  // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    auto B = std::make_shared<Matrix>(ndof, ninter);
    (*B)(0, 0) = -1.0;
    (*B)(1, 0) = 1.0;
    auto C = std::make_shared<Matrix>(ninter, ndof);
    (*C)(0, 0) = -1.0;
    (*C)(0, 1) = 1.0;
    auto myProcessRelation = std::make_shared<siconos::modeling::FirstOrderLinearR>();
    myProcessRelation->setConstantC(*C);
    myProcessRelation->setConstantB(*B);
    auto D = std::make_shared<Matrix>(ninter, ninter);
    (*D)(0, 0) = 1.0;

    myProcessRelation->setConstantD(*D);

    myProcessRelation->setComputeeVectorFunction(
        [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
          result(0) = user::computeControl(time);
        });

    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda
    auto myObserverRelation = std::make_shared<siconos::modeling::FirstOrderLinearR>();
    myObserverRelation->setConstantC(*C);
    myObserverRelation->setConstantB(*B);

    myObserverRelation->setConstantD(*D);
    myObserverRelation->setComputeeVectorFunction(
        [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
          result(0) = user::computeControl(time);
        });

    // NonSmoothLaw
    unsigned int nslawSize = 1;
    auto myNslaw = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(nslawSize);

    // The Interaction which involves the first DS (the process)
    auto myProcessInteraction =
        std::make_shared<siconos::modeling::Interaction>(myNslaw, myProcessRelation);

    // The Interaction which involves the second DS (the observer)
    auto myObserverInteraction =
        std::make_shared<siconos::modeling::Interaction>(myNslaw, myObserverRelation);

    // -------------
    // --- Model ---
    // -------------
    auto ObserverLCS = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    ObserverLCS->insertDynamicalSystem(process);
    ObserverLCS->insertDynamicalSystem(observer);
    ObserverLCS->link(myProcessInteraction, process);
    ObserverLCS->link(myObserverInteraction, observer);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // == Creation of the Simulation ==
    auto s = std::make_shared<siconos::simulation::TimeStepping>(ObserverLCS, td);
    // -- OneStepIntegrators --
    double theta = 0.5;
    auto myIntegrator = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);
    s->insertIntegrator(myIntegrator);

    // -- OneStepNsProblem --
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    s->insertNonSmoothProblem(osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    unsigned int outputSize = 10;         // number of required data
    unsigned int N = ceil((T - t0) / h);  // Number of time steps

    Matrix dataPlot(N, outputSize);

    auto xProc = process->x();
    auto xObs = observer->x();
    auto lambdaProc = myProcessInteraction->lambda(0);
    auto lambdaObs = myObserverInteraction->lambda(0);
    auto yProc = myProcessInteraction->y(0);
    auto yObs = myObserverInteraction->y(0);
    auto z = process->x();

    myProcessInteraction->computeOutput(t0, 0);

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = ObserverLCS->t0();  // Initial time of the model
    dataPlot(0, 1) = (*xProc)(0);
    dataPlot(0, 2) = (*xProc)(1);
    dataPlot(0, 3) = (*xObs)(0);
    dataPlot(0, 4) = (*xObs)(1);
    dataPlot(0, 5) = (*lambdaProc)(0);
    dataPlot(0, 6) = (*lambdaObs)(0);
    dataPlot(0, 7) = (*yProc)(0);
    dataPlot(0, 8) = (*yObs)(0);
    dataPlot(0, 9) = (*z)(0);

    // ==== Simulation loop =====
    cout << "====> Start computation ... \n\n";

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    unsigned int k = 0;  // Current step

    // Simulation loop
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    while (k < N - 1) {
      k++;
      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      dataPlot(k, 3) = (*xObs)(0);
      dataPlot(k, 4) = (*xObs)(1);
      dataPlot(k, 5) = (*lambdaProc)(0);
      dataPlot(k, 6) = (*lambdaObs)(0);
      dataPlot(k, 7) = (*yProc)(0);
      dataPlot(k, 8) = (*yObs)(0);
      dataPlot(k, 9) = (*z)(0);
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
    siconos::algebra::io::write("ObserverLCS.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-10;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "ObserverLCS.ref", eps)) > eps)
      return 1;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
