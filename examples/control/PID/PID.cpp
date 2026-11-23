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

/*!\file PID.cpp
  \brief \ref EMPID - C++ input file -
  O. Huber.

  Simple PID example.
  The controlled plant is a double integrator
  */

#include <LinearSensor.hpp>
#include <PID.hpp>
#include <SiconosControl.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <string>

using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  // ================= Creation of the model =======================

  // User-defined main parameters
  int nDof = 2;      // degrees of freedom for the system
  double t0 = 0;     // initial computation time
  double T = 100.0;  // final computation time
  double h = 0.05;   // time step
  double hControl = h;
  double position_init = 10;   // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double xFinal = 0.0;         // final value
  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  cout << "====> Model loading ...\n\n";

  auto A = std::make_shared<Matrix>(nDof, nDof);
  A->setZero();
  (*A)(0, 1) = 1.0;

  auto B = std::make_shared<Matrix>(nDof, 1);
  (*B)(1, 0) = 1;

  // -- Initial positions and velocities --
  auto x0 = std::make_shared<Vector>(nDof);
  (*x0)(0) = position_init;
  (*x0)(1) = velocity_init;

  // -- The dynamical system --
  auto doubleIntegrator =
      std::make_shared<siconos::modeling::FirstOrderLinearDS>(*x0, siconos::algebra::alias_t);
  doubleIntegrator->setConstantA(*A, siconos::algebra::alias_t);
  // -------------
  // --- Model ---
  // -------------
  auto sim = std::make_shared<siconos::control::ControlZOHSimulation>(t0, T, h);

  // add the dynamical system in the non smooth dynamical system
  sim->addDynamicalSystem(doubleIntegrator);

  // use a controlSensor
  auto C = std::make_shared<Matrix>(1, 2);
  C->setZero();
  (*C)(0, 0) = 1;
  auto sens = std::make_shared<siconos::control::LinearSensor>(doubleIntegrator, C);
  sim->addSensor(sens, hControl);
  // add the PID controller
  auto K = std::make_shared<Vector>(3);
  (*K)(0) = .25;
  (*K)(1) = .125;
  (*K)(2) = 2;
  auto act = std::make_shared<siconos::control::PID>(sens);
  act->setB(B);
  act->setRef(xFinal);
  act->setK(K);
  act->setDeltaT(h);
  sim->addActuator(act, hControl);

  cout << "=== End of model loading === \n";
  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Initialisation ...\n\n";
  // Initialize the model and the controlManager
  sim->initialize();
  sim->run();

  // --- Output files ---
  cout << "====> Output file writing ...\n";
  auto& dataPlot = *sim->data();
  siconos::algebra::io::write("PID.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Comparison with a reference file
  double error = 0.0, eps = 1e-12;
  if ((error = siconos::algebra::io::compareRefFile(dataPlot, "PID.ref", eps)) > eps)
    return 1;
  else
    return 0;
}
