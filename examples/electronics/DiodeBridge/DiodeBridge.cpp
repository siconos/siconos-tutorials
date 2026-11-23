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
//-----------------------------------------------------------------------
//
//  DiodeBridge  : sample of an electrical circuit involving :
//  - a linear dynamical system consisting of an LC oscillator (1 microF , 10 mH)
//  - a non smooth system (a 1000 Ohm resistor supplied through a 4
//    diodes bridge) in parallel with the oscillator
//
//  Expected behavior :

//  The initial state (Vc = 10 V , IL = 0) of the oscillator provides
//  an initial energy.
//  The period is 2 Pi sqrt(LC) ~ 0,628 ms.
//      The non smooth system is a full wave rectifier :
//  each phase (positive and negative) of the oscillation allows current to flow
//  through the resistor in a constant direction, resulting in an energy loss :
//  the oscillation damps.
//
//  State variables :
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  Since there is only one dynamical system, the interaction is defined by :
//  - complementarity laws between diodes current and
//    voltage. Depending on the diode position in the bridge, y stands
//    for the reverse voltage across the diode or for the diode
//    current (see figure in the template file)
//  - a linear time invariant relation between the state variables and
//    y and lambda (derived from Kirchhoff laws)
//
//-----------------------------------------------------------------------
#include <SiconosKernel.hpp>
#include <chrono>
#include <string>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  double t0 = 0.0;
  double T = 5.0e-3;                    // Total simulation time
  auto h_step = std::string("1.0e-6");  // Time step
  double Lvalue = 1e-2;                 // inductance
  double Cvalue = 1e-6;                 // capacitance
  double Rvalue = 1e3;                  // resistance
  double Vinit = 10.0;                  // initial voltage
  std::string Modeltitle = "DiodeBridge";

  auto start = std::chrono::system_clock::now();
  try {
    // --- Dynamical system specification ---
    Vector init_state{2};
    init_state << Vinit, 0.;
    auto LSDiodeBridge = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        init_state, siconos::algebra::alias_t);

    Matrix LS_A{2, 2};
    LS_A.setValue(0, 1, -1.0 / Cvalue);
    LS_A.setValue(1, 0, 1.0 / Lvalue);
    LSDiodeBridge->setConstantA(LS_A, siconos::algebra::alias_t);

    // --- Interaction between linear system and non smooth system ---
    Matrix int_C{4, 2};
    int_C.setZero();
    int_C(2, 0) = -1.0;
    int_C(3, 0) = 1.0;

    Matrix int_D{4, 4};
    int_D.setZero();
    int_D(0, 0) = 1.0 / Rvalue;
    int_D(0, 1) = 1.0 / Rvalue;
    int_D(0, 2) = -1.0;
    int_D(1, 0) = 1.0 / Rvalue;
    int_D(1, 1) = 1.0 / Rvalue;
    int_D(1, 3) = -1.0;
    int_D(2, 0) = 1.0;
    int_D(3, 1) = 1.0;

    Matrix int_B{2, 4};
    int_B.setZero();
    int_B(0, 2) = -1.0 / Cvalue;
    int_B(0, 3) = 1.0 / Cvalue;

    auto LTIRDiodeBridge =
        std::make_shared<siconos::modeling::FirstOrderLinearTIR>(int_C, int_B);
    LTIRDiodeBridge->setConstantD(int_D);

    auto nslaw = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(4);

    auto InterDiodeBridge =
        std::make_shared<siconos::modeling::Interaction>(nslaw, LTIRDiodeBridge);

    // --- Model creation ---
    auto DiodeBridge = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    DiodeBridge->setTitle(Modeltitle);
    // add the dynamical system in the non smooth dynamical system
    DiodeBridge->insertDynamicalSystem(LSDiodeBridge);
    // link the interaction and the dynamical system
    DiodeBridge->link(InterDiodeBridge, LSDiodeBridge);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    double theta = 0.5;
    auto aOSI = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);

    // -- (2) Time discretisation --
    auto aTiDisc = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h_step);

    // -- (3) Non smooth problem
    auto aLCP = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto aTS =
        std::make_shared<siconos::simulation::TimeStepping>(DiodeBridge, aTiDisc, aOSI, aLCP);

    int k = 0;
    double h = aTS->timeStep();
    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    Matrix dataPlot(N, 8);

    auto x = LSDiodeBridge->x();
    auto y = InterDiodeBridge->y(0);
    auto lambda = InterDiodeBridge->lambda(0);

    // For the initial time step:
    // time
    dataPlot(k, 0) = t0;

    // inductor voltage
    dataPlot(k, 1) = (*x)(0);

    // inductor current
    dataPlot(k, 2) = (*x)(1);

    // diode R1 current
    dataPlot(k, 3) = (*y)(0);

    // diode R1 voltage
    dataPlot(k, 4) = -(*lambda)(0);

    // diode F2 voltage
    dataPlot(k, 5) = -(*lambda)(1);

    // diode F1 current
    dataPlot(k, 6) = (*lambda)(2);

    // resistor current
    dataPlot(k, 7) = (*y)(0) + (*lambda)(2);

    auto start = std::chrono::system_clock::now();

    // --- Time loop  ---
    for (k = 1; k < N; ++k) {
      // solve ...
      aTS->computeOneStep();
      //  siconos::algebra::print(*aLCP);
      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = aTS->nextTime();

      // inductor voltage
      dataPlot(k, 1) = (*x)(0);

      // inductor current
      dataPlot(k, 2) = (*x)(1);

      // diode R1 current
      dataPlot(k, 3) = (*y)(0);

      // diode R1 voltage
      dataPlot(k, 4) = -(*lambda)(0);

      // diode F2 voltage
      dataPlot(k, 5) = -(*lambda)(1);

      // diode F1 current
      dataPlot(k, 6) = (*lambda)(2);

      // resistor current
      dataPlot(k, 7) = (*y)(0) + (*lambda)(2);

      aTS->nextStep();
    }

    // --- elapsed time computing ---
    std::cout << "time = \n";
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";

    // Number of time iterations
    std::cout << "Number of iterations done: " << k << "\n";

    siconos::algebra::io::write("DiodeBridge.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-11;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "DiodeBridge.ref", eps)) > eps)
      return 1;

  }
  // --- Exceptions handling ---
  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
