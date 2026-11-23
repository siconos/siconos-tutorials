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
//  DiodeBridgePowSup  : sample of an electrical circuit involving :
//  - a sinusoidal voltage source
//  - a non smooth system : a 4 diodes bridge used as a full wave rectifier
//        of the supplied voltage across the sinusoidal source, providing power
//        to the resistor load of the dynamical system
//  - a linear dynamical system LSDiodeBridgePowSup consisting of
//        a filtering capacitor in parallel with a load resistor
//
//  Expected behavior :
//      The non smooth system is a full wave rectifier :
//      each phase (positive and negative) of the supplied voltage allows current
//      to flow in a constant direction through the load.
//      The capacitor filter acts as a tank providing energy to the load resistor
//      when the voltage across the source weakens.
//      The load resistor consumes energy provided by the source.
//
//  State variable LSDiodeBridgePowSup:
//  - the voltage across the filtering capacitor
//
//  The interaction between the dynamical system and the source is defined by :
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
  double T = 5e-3;               // Total simulation time
  double h_step = 1.0e-6;        // Time step
  double Rvalue = 1e3;           // load resistance
  double Cfilt = 300.0e-9;       // filtering capacitor
  double VinitLS = 0.0;          // initial voltage Cfilt
  double DiodeThreshold = 0.21;  // Guess what ???
  double tinst;
  int k = 0;

  try {
    // --- Dynamical system creation ---
    // --- Linear system  (load and filter) specification ---
    // --- Dynamical system specification ---
    Vector init_stateLS{1};
    init_stateLS << VinitLS;
    auto LSDiodeBridgePowSup = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        init_stateLS, siconos::algebra::alias_t);

    Matrix LS_A{1, 1};
    LS_A(0, 0) = -1.0 / (Rvalue * Cfilt);

    LSDiodeBridgePowSup->setConstantA(LS_A, siconos::algebra::alias_t);

    // --- Interaction between linear system and non smooth system ---

    Matrix Int_C{4, 1};
    Int_C.setZero();
    Int_C(0, 0) = 1.0;
    Int_C(2, 0) = 1.0;

    Matrix Int_D{4, 4};
    Int_D.setZero();
    Int_D(0, 1) = -1.0;
    Int_D(1, 0) = 1.0;
    Int_D(1, 2) = 1.0;
    Int_D(1, 3) = -1.0;
    Int_D(2, 1) = -1.0;
    Int_D(3, 1) = 1.0;

    Matrix Int_B{1, 4};
    Int_B.setZero();
    Int_B(0, 0) = 1.0 / Cfilt;
    Int_B(0, 2) = 1.0 / Cfilt;

    auto LTIRDiodeBridgePowSup = std::make_shared<siconos::modeling::FirstOrderLinearR>();
    LTIRDiodeBridgePowSup->setConstantC(Int_C);
    LTIRDiodeBridgePowSup->setConstantB(Int_B);
    LTIRDiodeBridgePowSup->setConstantD(Int_D);

    Vector offset_y{4};
    offset_y.setZero();
    offset_y(0) = 1.0;
    offset_y(2) = 1.0;
    offset_y(3) = 1.0;
    offset_y *= -DiodeThreshold;

    Vector offset_lambda{4};
    offset_lambda.setZero();
    offset_lambda(1) = -DiodeThreshold;

    Vector Int_z{5};
    Int_z.setZero();
    Vector tmp{4};
    tmp = Int_D * offset_lambda - offset_y;
    Int_z.head(4) = tmp;
    Int_z(4) = 10.;

    LTIRDiodeBridgePowSup->setComputeeVectorFunction(
        [&Int_z](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
          double omega = 1e4;
          double Voffset = 0.0;
          double amplitude = 10.0;
          double phase = 0.0;
          double VSinPo;

          VSinPo = Voffset + (amplitude * cos((omega * time) + phase));

          result = Int_z.head(result.size());
          result(2) -= VSinPo;
          result(3) += VSinPo;
          Int_z(4) = VSinPo;
        });

    auto nslaw = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(4);

    auto InterDiodeBridgePowSup =
        std::make_shared<siconos::modeling::Interaction>(nslaw, LTIRDiodeBridgePowSup);

    // --- Model creation ---
    auto DiodeBridgePowSup =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    DiodeBridgePowSup->insertDynamicalSystem(LSDiodeBridgePowSup);
    // link the interaction and the dynamical system
    DiodeBridgePowSup->link(InterDiodeBridgePowSup, LSDiodeBridgePowSup);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    double theta = 0.5;
    auto aOSI = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);

    // -- (2) Time discretisation --
    auto aTiDisc = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h_step);

    // -- (3) Non smooth problem
    auto aLCP = std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_NSQP);

    // -- (4) Simulation setup with (1) (2) (3)
    auto aTS = std::make_shared<siconos::simulation::TimeStepping>(DiodeBridgePowSup, aTiDisc,
                                                                   aOSI, aLCP);

    k = 0;
    double h = aTS->timeStep();
    int N = ceil((T - t0) / h);  // Number of time steps
    std::cout << "Number of time steps = " << N << "\n";
    tinst = k * h_step;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    Matrix dataPlot(N + 1, 9);
    double i_DF1, i_DR1, i_DF2, i_DR2;
    double v_DF1, v_DR1, v_DF2, v_DR2;

    // For the initial time step:

    i_DF1 = (InterDiodeBridgePowSup->getLambda(0))(2);
    i_DR1 = (InterDiodeBridgePowSup->getLambda(0))(0);
    i_DF2 = (*InterDiodeBridgePowSup->y(0))(1);
    i_DR2 = (InterDiodeBridgePowSup->getLambda(0))(3);

    v_DF1 = -(*InterDiodeBridgePowSup->y(0))(2) + DiodeThreshold;
    v_DR1 = -(*InterDiodeBridgePowSup->y(0))(0) + DiodeThreshold;
    v_DF2 = -(InterDiodeBridgePowSup->getLambda(0))(1) + DiodeThreshold;
    v_DR2 = -(*InterDiodeBridgePowSup->y(0))(3) + DiodeThreshold;

    // time
    dataPlot(k, 0) = DiodeBridgePowSup->t0();

    // source voltage
    dataPlot(k, 1) = Int_z(4);

    // source current
    dataPlot(k, 2) = i_DF1 - i_DR2;

    // diode R1 current
    dataPlot(k, 3) = i_DR1;

    // diode R1 voltage
    dataPlot(k, 4) = v_DR1;

    // diode F2 voltage
    dataPlot(k, 5) = v_DF2;

    // diode F1 current
    dataPlot(k, 6) = i_DF1;

    // diode F2 current
    dataPlot(k, 7) = i_DF2;

    // diode R2 current
    dataPlot(k, 8) = i_DR2;

    // --- Compute elapsed time ---
    auto start = std::chrono::system_clock::now();
    // --- Time loop  ---
    while (k < N) {
      // get current time step
      k++;
      tinst = k * h_step;

      // solve ...
      aTS->computeOneStep();

      // --- Get values to be plotted ---
      i_DF1 = (InterDiodeBridgePowSup->getLambda(0))(2);
      i_DR1 = (InterDiodeBridgePowSup->getLambda(0))(0);
      i_DF2 = (*InterDiodeBridgePowSup->y(0))(1);
      i_DR2 = (InterDiodeBridgePowSup->getLambda(0))(3);

      v_DF1 = -(*InterDiodeBridgePowSup->y(0))(2) + DiodeThreshold;
      v_DR1 = -(*InterDiodeBridgePowSup->y(0))(0) + DiodeThreshold;
      v_DF2 = -(InterDiodeBridgePowSup->getLambda(0))(1) + DiodeThreshold;
      v_DR2 = -(*InterDiodeBridgePowSup->y(0))(3) + DiodeThreshold;

      // time
      dataPlot(k, 0) = aTS->nextTime();

      // source voltage
      dataPlot(k, 1) = Int_z(4);

      // source current
      dataPlot(k, 2) = i_DF1 - i_DR2;

      // diode R1 current
      dataPlot(k, 3) = i_DR1;

      // diode R1 voltage
      dataPlot(k, 4) = v_DR1;

      // diode F2 voltage
      dataPlot(k, 5) = v_DF2;

      // diode F1 current
      dataPlot(k, 6) = i_DF1;

      // diode F2 current
      dataPlot(k, 7) = i_DF2;

      // diode R2 current
      dataPlot(k, 8) = i_DR2;

      aTS->nextStep();
    }

    // --- elapsed time computing ---
    std::cout << "time = \n";
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";

    // Number of time iterations
    std::cout << "Number of iterations done: " << k << "\n";

    siconos::algebra::io::write("DiodeBridgePowSup.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "DiodeBridgePowSup.ref",
                                                      eps)) > eps)
      return 1;
  }
  // --- Exceptions handling ---
  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
