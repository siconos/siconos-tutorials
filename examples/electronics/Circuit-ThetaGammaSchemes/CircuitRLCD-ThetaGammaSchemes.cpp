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
//  This circuit corresponds to the circuit (b) p 20, equation 1.39 of
//  "Nonsmooth modeling and simulation for switched circuits"
//  Acary, V., Brogliato, B. and Bonnefon, B.
//  LNEE 69, Springer Verlag 2010
//-----------------------------------------------------------------------

#include <SiconosKernel.hpp>
#include <chrono>
#include <string>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  double t0 = 0.0;
  double T = 5.0;                           // Total simulation time
  double h_step = 5e-2;                     // Time step
  double Lvalue = 1.0;                      // inductance
  double Cvalue = 1.0 / (4 * M_PI * M_PI);  // capacitance
  double Rvalue = 10.0;                     // resistance
  std::string Modeltitle = "CircuitRLCD";

  try {
    // --- Dynamical system specification ---
    auto init_state = std::make_shared<Vector>(2);
    (*init_state)(0) = 1.0;
    (*init_state)(1) = 1.0;

    auto LS_A = std::make_shared<Matrix>(2, 2);
    LS_A->setValue(0, 1, -1.0);
    LS_A->setValue(1, 0, 1.0 / (Lvalue * Cvalue));

    auto LSCircuitRLCD = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        *init_state, siconos::algebra::alias_t);
    LSCircuitRLCD->setConstantA(*LS_A, siconos::algebra::alias_t);

    // --- Interaction between linear system and non smooth system ---

    auto Int_C = std::make_shared<Matrix>(1, 2);
    Int_C->setValue(0, 0, 1.0 / Cvalue);

    auto Int_D = std::make_shared<Matrix>(1, 1);
    Int_D->setValue(0, 0, Rvalue);

    auto Int_B = std::make_shared<Matrix>(2, 1);
    Int_B->setValue(0, 0, 1.0);

    auto LTIRCircuitRLCD =
        std::make_shared<siconos::modeling::FirstOrderLinearTIR>(*Int_C, *Int_B);
    auto NSLaw = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(1);

    LTIRCircuitRLCD->setConstantD(*Int_D);

    auto InterCircuitRLCD =
        std::make_shared<siconos::modeling::Interaction>(NSLaw, LTIRCircuitRLCD);

    // --- Model creation ---
    auto CircuitRLCD = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    CircuitRLCD->setTitle(Modeltitle);
    CircuitRLCD->insertDynamicalSystem(LSCircuitRLCD);
    CircuitRLCD->link(InterCircuitRLCD, LSCircuitRLCD);

    InterCircuitRLCD->computeOutput(t0, 0);
    InterCircuitRLCD->computeInput(t0, 0);

    // ------------------
    // --- Simulation ---
    // ------------------
    double theta = 0.500000000000;
    double gamma = 0.500000000000;

    // -- (1) OneStepIntegrators --
    auto OSI_RLCD = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta, gamma);
    OSI_RLCD->setUseGammaForRelation(true);
    // -- (2) Time discretisation --
    auto TiDiscRLCD = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h_step);
    // --- (3) one step non smooth problem
    auto LCP_RLCD = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto StratCircuitRLCD = std::make_shared<siconos::simulation::TimeStepping>(
        CircuitRLCD, TiDiscRLCD, OSI_RLCD, LCP_RLCD);

    double h = StratCircuitRLCD->timeStep();
    int N = ceil((T - t0) / h);  // Number of time steps
    int k = 0;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    Matrix dataPlot(N, 9);

    // For the initial time step:

    // time
    dataPlot(k, 0) = CircuitRLCD->t0();

    // inductor voltage
    dataPlot(k, 1) = (*LSCircuitRLCD->x())(0);

    // inductor current
    dataPlot(k, 2) = (*LSCircuitRLCD->x())(1);

    // diode voltage
    dataPlot(k, 3) = -(*InterCircuitRLCD->y(0))(0);

    // diode current
    dataPlot(k, 4) = (InterCircuitRLCD->getLambda(0))(0);

    dataPlot(k, 5) = (*LSCircuitRLCD->r())(0);

    // Storage function
    dataPlot(k, 6) = 1.0 / 2.0 *
                     (1.0 / (Cvalue) * (*LSCircuitRLCD->x())(0) * (*LSCircuitRLCD->x())(0) +
                      Lvalue * (*LSCircuitRLCD->x())(1) * (*LSCircuitRLCD->x())(1));
    // Dissipation function
    // dataPlot(k,7) =Rvalue*(InterCircuitRLCD->getLambda(0))(0) ;
    dataPlot(k, 7) = 0.0;

    // Total energy
    dataPlot(k, 8) = dataPlot(k, 6) + dataPlot(k, 7);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    // --- Time loop  ---
    for (k = 1; k < N; ++k) {
      // solve ...
      StratCircuitRLCD->computeOneStep();
      // siconos::algebra::print(*LCP_RLCD);
      //  --- Get values to be plotted ---
      //  time
      dataPlot(k, 0) = StratCircuitRLCD->nextTime();

      // inductor voltage
      dataPlot(k, 1) = (*LSCircuitRLCD->x())(0);

      // inductor current
      dataPlot(k, 2) = (*LSCircuitRLCD->x())(1);

      // diode voltage
      dataPlot(k, 3) = -(*InterCircuitRLCD->y(0))(0);

      // diode current
      dataPlot(k, 4) = (InterCircuitRLCD->getLambda(0))(0);

      // dataPlot(k,5) = (LSCircuitRLCD->getR())(0);
      //     dataPlot(k,5) = OSI_RLCD->computeResidu();
      dataPlot(k, 5) = 0;

      // Storage function
      dataPlot(k, 6) = 1.0 / 2.0 *
                       (1.0 / Cvalue * (*LSCircuitRLCD->x())(0) * (*LSCircuitRLCD->x())(0) +
                        Lvalue * (*LSCircuitRLCD->x())(1) * (*LSCircuitRLCD->x())(1));

      // Dissipation function
      dataPlot(k, 7) = h * Rvalue * (InterCircuitRLCD->getLambda(0))(0) *
                           (InterCircuitRLCD->getLambda(0))(0) +
                       dataPlot(k - 1, 7);

      // Total energy function
      // Total energy
      dataPlot(k, 8) = dataPlot(k, 6) + dataPlot(k, 7);

      // transfer of state i+1 into state i and time incrementation
      StratCircuitRLCD->nextStep();
    }
    // Number of time iterations
    std::cout << "Number of iterations done: " << k - 1 << "\n";
    std::cout << "Computation Time \n";
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";

    // dataPlot (ascii) output
    siconos::algebra::io::write("CircuitRLCD.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "CircuitRLCD.ref", eps)) > eps)
      return 1;

  }

  // --- Exceptions handling ---
  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
