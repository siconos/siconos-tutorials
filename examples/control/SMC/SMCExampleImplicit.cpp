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

/* Two independent systems of dimension one controlled to slide
  on \f$x = 0\f$. An implicit scheme is used
  O. Huber
*/

#include <SiconosControl.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <string>
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;
using namespace std;

// main program
int main(int argc, char* argv[]) {
  // User-defined parameters
  unsigned int ndof = 2;     // Number of degrees of freedom of your system
  double t0 = 0.0;           // Starting time
  double T = 1;              // Total simulation time
  double h = 1.0e-4;         // Time step for simulation
  double hControl = 1.0e-2;  // Time step for control
  double Xinit = 1.0;

  if (h > hControl) {
    THROW_EXCEPTION("hControl must be bigger than h");
  }

  // ================= Creation of the model =======================
  // Steps:
  // - create a Dynamical System
  // - add a Simulation to the model

  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  // First System:
  // dx/dt = Ax + u(t) + r
  // x(0) = x0
  // Note: r = Blambda, B defines in relation below.

  // Matrix declaration
  auto x0 = std::make_shared<Vector>(ndof);
  (*x0)(0) = Xinit;
  (*x0)(1) = -Xinit;
  auto sensorC = std::make_shared<Matrix>(2, 2);
  sensorC->setIdentity();
  auto sensorD = std::make_shared<Matrix>(2, 2);
  sensorD->setZero();
  auto Csurface = std::make_shared<Matrix>(1, 2);
  Csurface->setZero();
  (*Csurface)(0, 1) = 1;
  auto Brel = std::make_shared<Matrix>(2, 1);
  Brel->setZero();
  (*Brel)(1, 0) = 2;

  // Dynamical Systems
  auto processDS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
      *x0, siconos::algebra::SiconosDenseMatrix::Zero(x0->size(), x0->size()),
      siconos::algebra::SiconosVector::Zero(x0->size()), siconos::algebra::copy_t);
  processDS->setComputebVectorFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        auto t = sin(50 * time);
        result(0) = t;
        result(1) = -t;
      });

  // -------------
  // --- Model process ---
  // -------------
  auto sim = std::make_shared<siconos::control::ControlZOHSimulation>(t0, T, h);
  // sim->setSaveOnlyMainSimulation(true);
  sim->addDynamicalSystem(processDS);

  // ------------------
  // --- Simulation ---
  // ------------------
  // Control stuff
  // use a controlSensor
  auto sens = std::make_shared<siconos::control::LinearSensor>(processDS, sensorC, sensorD);
  sim->addSensor(sens, hControl);
  // add the sliding mode controller
  auto act = std::make_shared<siconos::control::LinearSMC>(sens);
  act->setCsurface(Csurface);
  act->setB(Brel);
  sim->addActuator(act, hControl);
  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Simulation initialisation ...\n\n";
  // initialise the process and the ControlManager
  sim->initialize();

  // ==== Simulation loop =====
  cout << "====> Start computation ... \n\n";
  sim->run();
  // --- Output files ---
  cout << "====> Output file writing ...\n";
  auto& dataPlot = *sim->data();

  siconos::algebra::io::write("SMCExampleImplicit.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);

  double error = 0.0, eps = 1e-12;
  if ((error = siconos::algebra::io::compareRefFile(dataPlot, "SMCExampleImplicit.ref", eps)) >
      eps)
    return 1;
  else
    return 0;
}
