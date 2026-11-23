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
//-----------------------------------------------------------------------
//
//  DiodeBridgeCapFilter  : sample of an electrical circuit involving :
//  - a 1st linear dynamical system LSDiodeBridge1 consisting of
//        an LC oscillator (1 microF , 10 mH)
//  - a non smooth system : a 4 diodes bridge used as a full wave rectifier
//        of the supplied voltage across the LC oscillator, providing power
//    to the resistor load of the 2nd dynamical system
//      - a 2nd linear dynamical system LSDiodeBridge2 consisting of
//        a filtering capacitor in parallel with a load resistor
//
//  Expected behavior :
//  The initial state (Vc = 10 V , IL = 0) of the oscillator provides
//      an initial energy.
//  The oscillator period is 2 Pi sqrt(LC) ~ 0,628 ms.
//      The non smooth system is a full wave rectifier :
//  each phase (positive and negative) of the oscillation allows current
//      to flow in a constant direction through the load.
//      The capacitor filter acts as a tank providing energy to the load resistor
//      when the voltage across the oscillator weakens.
//      The load resistor consumes energy : the oscillation damps.
//
//  State variables LSDiodeBridge1:
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  State variable LSDiodeBridge2:
//  - the voltage across the filtering capacitor
//
//  The interaction between the two dynamical systems is defined by :
//  - complementarity laws between diodes current and voltage. Depending on
//        the diode position in the bridge, y stands for the reverse voltage across
//    the diode or for the diode current (see figure in the template file)
//  - a linear time invariant relation between the state variables and y and
//    lambda (derived from the Kirchhoff laws)
//
//-----------------------------------------------------------------------
#include <SiconosKernel.hpp>
#include <chrono>
#include <string>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  double t0 = 0.0;
  double T = 5e-3;          // Total simulation time
  double h_step = 1e-6;     // Time step
  double Lvalue = 1e-2;     // inductance
  double Cvalue = 1e-6;     // capacitance LC oscillator
  double Rvalue = 1e3;      // load resistance
  double Cfilt = 300.0e-9;  // filtering capacitor
  double VinitLS1 = 10.0;   // initial voltage LC oscillator
  double VinitLS2 = 0.0;    // initial voltage Cfilt
  std::string Modeltitle = "DiodeBridgeCapFilter";

  try {
    // --- Linear system 1 (LC oscillator) specification ---
    Vector init_stateLS1{2};
    init_stateLS1 << VinitLS1, 0.;

    auto LS1DiodeBridgeCapFilter = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        init_stateLS1, siconos::algebra::alias_t);
    Matrix LS1_A{2, 2};
    LS1_A.setZero();
    LS1_A(0, 1) = -1.0 / Cvalue;
    LS1_A(1, 0) = 1.0 / Lvalue;

    LS1DiodeBridgeCapFilter->setConstantA(LS1_A, siconos::algebra::alias_t);
    LS1DiodeBridgeCapFilter->display();

    // --- Linear system 2 (load and filter) specification ---
    Vector init_stateLS2{1};
    init_stateLS2 << VinitLS2;
    auto LS2DiodeBridgeCapFilter = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        init_stateLS2, siconos::algebra::alias_t);
    Matrix LS2_A{1, 1};
    LS2_A(0, 0) = -1.0 / (Rvalue * Cfilt);
    LS2DiodeBridgeCapFilter->setConstantA(LS2_A, siconos::algebra::alias_t);

    // --- Interaction between linear systems and non smooth system ---
    Matrix int_C{4, 3};
    int_C.setZero();
    int_C(0, 2) = 1.0;
    int_C(2, 0) = -1.0;
    int_C(2, 2) = 1.0;
    int_C(3, 0) = 1.0;

    Matrix int_D{4, 4};
    int_D.setZero();
    int_D(0, 1) = -1.0;
    int_D(1, 0) = 1.0;
    int_D(1, 2) = 1.0;
    int_D(1, 3) = -1.0;
    int_D(2, 1) = -1.0;
    int_D(3, 1) = 1.0;

    Matrix int_B{3, 4};
    int_B.setZero();
    int_B(0, 2) = -1.0 / Cvalue;
    int_B(0, 3) = 1.0 / Cvalue;
    int_B(2, 0) = 1.0 / Cfilt;
    int_B(2, 2) = 1.0 / Cfilt;

    auto LTIRDiodeBridgeCapFilter =
        std::make_shared<siconos::modeling::FirstOrderLinearTIR>(int_C, int_B);
    LTIRDiodeBridgeCapFilter->setConstantD(int_D);
    auto nslaw = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(4);

    auto InterDiodeBridgeCapFilter =
        std::make_shared<siconos::modeling::Interaction>(nslaw, LTIRDiodeBridgeCapFilter);

    // --- Model creation ---
    auto DiodeBridgeCapFilter =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    DiodeBridgeCapFilter->setTitle(Modeltitle);
    // add the dynamical system in the non smooth dynamical system
    DiodeBridgeCapFilter->insertDynamicalSystem(LS1DiodeBridgeCapFilter);
    DiodeBridgeCapFilter->insertDynamicalSystem(LS2DiodeBridgeCapFilter);
    // link the interaction and the dynamical system
    DiodeBridgeCapFilter->link(InterDiodeBridgeCapFilter, LS1DiodeBridgeCapFilter,
                               LS2DiodeBridgeCapFilter);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    double theta = 1.0;
    auto aOSI = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);
    // -- (2) Time discretisation --
    auto aTiDisc = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h_step);

    // -- (3) Non smooth problem
    auto aLCP = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto aTS = std::make_shared<siconos::simulation::TimeStepping>(DiodeBridgeCapFilter,
                                                                   aTiDisc, aOSI, aLCP);

    int k = 0;
    double h = aTS->timeStep();
    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    Matrix dataPlot{N, 8};
    // For the initial time step:
    // time
    dataPlot(k, 0) = DiodeBridgeCapFilter->t0();
    // inductor voltage
    dataPlot(k, 1) = (*LS1DiodeBridgeCapFilter->x())(0);
    // inductor current
    dataPlot(k, 2) = (*LS1DiodeBridgeCapFilter->x())(1);

    // diode R1 current
    dataPlot(k, 3) = (InterDiodeBridgeCapFilter->getLambda(0))(0);

    // diode R1 voltage
    dataPlot(k, 4) = -(*InterDiodeBridgeCapFilter->y(0))(0);

    // diode F2 voltage
    dataPlot(k, 5) = -(InterDiodeBridgeCapFilter->getLambda(0))(1);

    // diode F1 current
    dataPlot(k, 6) = (InterDiodeBridgeCapFilter->getLambda(0))(2);

    // load voltage
    dataPlot(k, 7) = (*LS2DiodeBridgeCapFilter->x())(0);

    // --- Compute elapsed time ---
    auto start = std::chrono::system_clock::now();
    // --- Time loop  ---
    for (k = 1; k < N; ++k) {
      // solve ...
      aTS->computeOneStep();
      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = aTS->nextTime();

      // inductor voltage
      dataPlot(k, 1) = (*LS1DiodeBridgeCapFilter->x())(0);

      // inductor current
      dataPlot(k, 2) = (*LS1DiodeBridgeCapFilter->x())(1);

      // diode R1 current
      dataPlot(k, 3) = (InterDiodeBridgeCapFilter->getLambda(0))(0);

      // diode R1 voltage
      dataPlot(k, 4) = -(*InterDiodeBridgeCapFilter->y(0))(0);

      // diode F2 voltage
      dataPlot(k, 5) = -(InterDiodeBridgeCapFilter->getLambda(0))(1);

      // diode F1 current
      dataPlot(k, 6) = (InterDiodeBridgeCapFilter->getLambda(0))(2);
      // load voltage
      dataPlot(k, 7) = (*LS2DiodeBridgeCapFilter->x())(0);

      aTS->nextStep();
    }

    // --- elapsed time computing ---
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";

    // Number of time iterations
    std::cout << "Number of iterations done: " << k << "\n";

    siconos::algebra::io::write("DiodeBridgeCapFilter.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    std::vector<int> idx(4);
    for (auto i = 0; i < 4; i++) idx.push_back(i);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "DiodeBridgeCapFilter.ref",
                                                      eps, idx)) > eps)
      return 1;
  }
  // --- Exceptions handling ---
  catch (...) {
    siconos::exception::process();
    return 1;
  }
  return 0;
}
