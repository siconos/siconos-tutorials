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
#include <string>

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
    double T = 20.0;    // Total simulation times
    double h = 1.0e-1;  // Time step
    double Vinit = 10.0;

    double G = 10.0;
    double beta = .3;

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
    auto x0 = std::make_shared<Vector>(ndof);
    (*x0)(0) = Vinit;
    (*x0)(1) = Vinit;
    auto doubleIntegrator = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        *x0, siconos::algebra::alias_t);
    doubleIntegrator->setConstantA(*A, siconos::algebra::alias_t);
    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 2;  // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the doubleIntegrator
    // y = Cx + Dlambda
    // r = Blambda
    auto B = std::make_shared<Matrix>(ndof, ninter);
    B->setZero();
    (*B)(1, 0) = G;
    (*B)(1, 1) = G * beta;
    auto C = std::make_shared<Matrix>(ninter, ndof);
    C->setIdentity();
    auto twistingRelation = std::make_shared<siconos::modeling::FirstOrderLinearTIR>(*C, *B);

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    auto H = std::make_shared<Matrix>(4, 2);
    H->setZero();
    (*H)(0, 0) = 1.0;
    (*H)(1, 0) = -h / 2.0;
    (*H)(2, 0) = -1.0;
    (*H)(3, 0) = h / 2.0;
    (*H)(1, 1) = 1.0;
    (*H)(3, 1) = -1.0;

    auto K = std::make_shared<Vector>(4);
    K->setConstant(-1.);
    auto nslaw = std::make_shared<siconos::modeling::NormalConeNSL>(nslawSize, H, K);

    auto twistingInteraction =
        std::make_shared<siconos::modeling::Interaction>(nslaw, twistingRelation);

    // -------------
    // --- Model ---
    // -------------
    auto itw = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    itw->insertDynamicalSystem(doubleIntegrator);
    itw->link(twistingInteraction, doubleIntegrator);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // == Creation of the Simulation ==
    auto s = std::make_shared<siconos::simulation::TimeStepping>(itw, td);
    // -- OneStepIntegrators --
    double theta = 0.5;
    auto integrator = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);
    s->insertIntegrator(integrator);
    // -- OneStepNsProblem --

    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::AVI>();
    s->insertNonSmoothProblem(osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    unsigned outputSize = 5;              // number of required data
    unsigned N = ceil((T - t0) / h) + 1;  // Number of time steps

    auto dataPlot = std::make_shared<Matrix>(N, outputSize);

    auto& xProc = *doubleIntegrator->x();
    auto& lambdaProc = *twistingInteraction->lambda(0);

    // -> saved in a matrix dataPlot
    (*dataPlot)(0, 0) = itw->t0();  // Initial time of the model
    (*dataPlot)(0, 1) = xProc(0);
    (*dataPlot)(0, 2) = xProc(1);
    (*dataPlot)(0, 3) = -1.0;

    (*dataPlot)(0, 4) = -1.0;

    // ==== Simulation loop =====
    cout << "====> Start computation ... \n\n";

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    unsigned int k = 0;  // Current step

    // Simulation loop
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      k++;
      //  osnspb->setNumericsVerboseMode(1);

      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      (*dataPlot)(k, 0) = s->nextTime();
      (*dataPlot)(k, 1) = xProc(0);
      (*dataPlot)(k, 2) = xProc(1);
      (*dataPlot)(k, 3) = lambdaProc(0);
      (*dataPlot)(k, 4) = lambdaProc(1);
      s->nextStep();
    }

    cout << "End of computation - Number of iterations done: " << k - 1 << endl;
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms\n";
    cout << "====> Output file writing ...\n";

    siconos::algebra::io::write("Twisting.dat", *dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // We do not compare the Lagrange multiplier that are very
    // sensitive to numerical approximations
    std::vector<int> idx = {0, 1, 2};
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(*dataPlot, "Twisting.ref", eps, idx)) >
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
