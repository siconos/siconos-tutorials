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
//  CircuitRLCD  : sample of an electrical circuit involving :
//  - a linear dynamical system consisting of an LC oscillator (1 microF , 10 mH)
//  - a non smooth system (a 1000 Ohm resistor in series with a diode) in parallel
//    with the oscillator
//
//  Expected behavior :
//  The initial state of the oscillator provides an initial energy.
//  The period is 2 Pi sqrt(LC) ~ 0,628 ms.
//  A positive voltage across the capacitor allows current to flow
//  through the resistor-diode branch , resulting in an energy loss :
//  the oscillation damps.
//
//  State variables :
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  Since there is only one dynamical system, the interaction is defined by :
//  - a complementarity law between diode current and voltage where y stands
//    for the reverse voltage across the diode and lambda stands for the
//    the diode current
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
  double T = 5e-3;          // Total simulation time
  double h_step = 10.0e-6;  // Time step
  double Lvalue = 1e-2;     // inductance
  double Cvalue = 1e-6;     // capacitance
  double Rvalue = 1e3;      // resistance
  double Vinit = 10.0;      // initial voltage

  try {
    // ================= Creation of the model =======================
    // --- Dynamical system specification ---
    siconos::algebra::SiconosVector init_state{2};
    init_state << Vinit, 0.;

    siconos::algebra::SiconosMatrix LS_A{2, 2};
    LS_A.setZero();
    LS_A(0, 1) = -1.0 / Cvalue;
    LS_A(1, 0) = 1.0 / Lvalue;

    auto LSCircuitRLCD = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        init_state, siconos::algebra::alias_t);
    LSCircuitRLCD->setConstantA(LS_A, siconos::algebra::alias_t);
    // --- Interaction between linear system and non smooth system ---
    siconos::algebra::SiconosMatrix Int_C{1, 2};
    Int_C(0, 0) = -1.0;
    Int_C(0, 1) = 0.;

    siconos::algebra::SiconosMatrix Int_D{1, 1};
    Int_D << Rvalue;

    siconos::algebra::SiconosMatrix Int_B{2, 1};
    Int_B << -1.0 / Cvalue, 0.;

    auto LTIRCircuitRLCD =
        std::make_shared<siconos::modeling::FirstOrderLinearTIR>(Int_C, Int_B);
    auto NSLaw = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(1);

    LTIRCircuitRLCD->setConstantD(Int_D);

    auto InterCircuitRLCD =
        std::make_shared<siconos::modeling::Interaction>(NSLaw, LTIRCircuitRLCD);

    // --- Model creation ---
    auto CircuitRLCD = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    CircuitRLCD->insertDynamicalSystem(LSCircuitRLCD);

    // link the interaction and the dynamical system
    CircuitRLCD->link(InterCircuitRLCD, LSCircuitRLCD);

    InterCircuitRLCD->computeOutput(t0, 0);
    InterCircuitRLCD->computeInput(t0, 0);

    // ------------------
    // --- Simulation ---
    // ------------------
    double theta = 0.5000000000001;

    // -- (1) OneStepIntegrators --
    auto OSI_RLCD = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);

    // -- (2) Time discretisation --
    auto TiDiscRLCD = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h_step);
    // --- (3) one step non smooth problem
    auto LCP_RLCD = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto StratCircuitRLCD = std::make_shared<siconos::simulation::TimeStepping>(
        CircuitRLCD, TiDiscRLCD, OSI_RLCD, LCP_RLCD);
    double h = StratCircuitRLCD->timeStep();
    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    Matrix dataPlot(N + 1, 6);
    dataPlot.setZero();
    // For the initial time step:
    LSCircuitRLCD->display();
    auto x = LSCircuitRLCD->x_read();
    siconos::algebra::print(x);
    // time
    dataPlot(0, 0) = CircuitRLCD->t0();

    // inductor voltage
    dataPlot(0, 1) = x(0);  // x(0);

    // inductor current
    dataPlot(0, 2) = x(1);

    // diode voltage
    dataPlot(0, 3) = -(*InterCircuitRLCD->y(0))(0);

    // diode current
    dataPlot(0, 4) = (InterCircuitRLCD->getLambda(0))(0);

    dataPlot(0, 5) = (*LSCircuitRLCD->r())(0);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    // --- Time loop  ---
    int k = 1;
    while (StratCircuitRLCD->hasNextEvent()) {
      // solve ...
      StratCircuitRLCD->computeOneStep();

      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = StratCircuitRLCD->nextTime();

      // inductor voltage
      dataPlot(k, 1) = x(0);

      // inductor current
      dataPlot(k, 2) = x(1);

      // diode voltage
      dataPlot(k, 3) = -(*InterCircuitRLCD->y(0))(0);

      // diode current
      dataPlot(k, 4) = (InterCircuitRLCD->getLambda(0))(0);

      // dataPlot(k,5) = (LSCircuitRLCD->getR())(0);
      //     dataPlot(k,5) = OSI_RLCD->computeResidu();
      dataPlot(k, 5) = 0;
      // transfer of state i+1 into state i and time incrementation
      StratCircuitRLCD->nextStep();
      k++;
    }
    // Number of time iterations
    std::cout << "Number of iterations done: " << k - 1 << "\n";
    std::cout << "Computation Time \n";
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";
    std::cout << dataPlot(0, 1) << "\n";
    // dataPlot (ascii) output
    siconos::algebra::io::write("CircuitRLCD.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-11;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "CircuitRLCD.ref", eps)) > eps)
      return 1;
  }

  // --- Exceptions handling ---
  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
