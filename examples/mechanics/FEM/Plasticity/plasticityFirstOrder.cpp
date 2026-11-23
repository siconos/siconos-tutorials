
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#include <stdio.h>
#include <stdlib.h>

#include <chrono>

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[]) {
  double t0 = 0.0;
  double T = 35e-1;        // Total simulation time
  double h_step = 1.0e-3;  // Time step

  double m = 1;   // mass
  double k = 10;  // stiffness of the spring
  double b = 1;   // B matrix, trivial here

  double vInit = 1.0;            // initial velocity
  double xInit = 0.0;            // initial position
  double sigmaInit = k * xInit;  // initial stress
  double sigmaMax = 1.5;

  string Modeltitle = "SpringMass";

  try {
    // --- Dynamical system specification ---
    SP::SiconosVector init_state(new SiconosVector(3));
    init_state->setValue(0, vInit);
    init_state->setValue(1, xInit);
    init_state->setValue(2, sigmaInit);
    SP::SimpleMatrix M(new SimpleMatrix(3, 3));
    M->setValue(0, 0, m);
    M->setValue(1, 1, 1.0);
    //    M->setValue(2, 2, 1/k);
    M->setValue(2, 2, 1.0);
    SP::SimpleMatrix A(new SimpleMatrix(3, 3));
    A->setValue(0, 2, -b);
    A->setValue(1, 0, 1);
    //    A->setValue(2, 0, b);
    A->setValue(2, 0, k * b);

    SP::SiconosVector init_state_1D(new SiconosVector(2));
    init_state_1D->setValue(0, vInit);
    init_state_1D->setValue(1, xInit);
    SP::SimpleMatrix M_1D(new SimpleMatrix(2, 2));
    M_1D->setValue(0, 0, m);
    M_1D->setValue(1, 1, 1);
    SP::SimpleMatrix A_1D(new SimpleMatrix(2, 2));
    A_1D->setValue(0, 1, -k);
    A_1D->setValue(1, 0, 1);

    SP::FirstOrderLinearTIDS plasticSpring(new FirstOrderLinearTIDS(init_state, A));
    //    SP::FirstOrderLinearTIDS plasticSpring(new FirstOrderLinearTIDS(init_state_1D,
    //    A_1D)); plasticSpring->setMPtr(M_1D);
    plasticSpring->setMPtr(M);
    plasticSpring->display();

    // --- Interaction between linear system and non smooth system ---
    SP::SimpleMatrix C(new SimpleMatrix(2, 3));
    C->setValue(0, 2, -1.0);
    C->setValue(1, 2, 1.0);
    SP::SimpleMatrix B(new SimpleMatrix(3, 2));
    B->setValue(2, 0, -1.0);
    B->setValue(2, 1, 1.0);
    SP::SiconosVector e(new SiconosVector(2));
    e->setValue(0, sigmaMax);
    e->setValue(1, sigmaMax);
    SP::SimpleMatrix C_1D(new SimpleMatrix(1, 2));
    C_1D->setValue(0, 1, 1.0);

    SP::SimpleMatrix B_1D(new SimpleMatrix(2, 1));
    B_1D->setValue(0, 0, 1.0);

    SP::SiconosVector e_1D(new SiconosVector(3));
    e_1D->setValue(0, 1);

    SP::FirstOrderLinearTIR LTIRspring(new FirstOrderLinearTIR(C, B));
    //    SP::FirstOrderLinearTIR LTIRspring(new FirstOrderLinearTIR(C_1D, B_1D));
    LTIRspring->setePtr(e);
    SP::NonSmoothLaw NSLaw(new ComplementarityConditionNSL(2));
    //    SP::NonSmoothLaw NSLaw(new RelayNSL(1));

    //    LTIRCircuitRLCD->setDPtr(Int_D);

    SP::Interaction InterSpring(new Interaction(NSLaw, LTIRspring));
    InterSpring->display();
    // --- Model creation ---
    SP::NonSmoothDynamicalSystem springDS(new NonSmoothDynamicalSystem(t0, T));
    springDS->setTitle(Modeltitle);
    // add the dynamical system in the non smooth dynamical system
    springDS->insertDynamicalSystem(plasticSpring);

    // link the interaction and the dynamical system
    springDS->link(InterSpring, plasticSpring);

    //    InterSpring->computeOutput(t0,0);
    //    InterSpring->computeInput(t0,0);

    springDS->display();

    // ------------------
    // --- Simulation ---
    // ------------------
    double theta = 0.5000000000001;

    // -- (1) OneStepIntegrators --
    SP::EulerMoreauOSI OSI(new EulerMoreauOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation TiDis(new TimeDiscretisation(t0, h_step));
    // --- (3) one step non smooth problem
    SP::LCP LCP_spring(new LCP());
    //    SP::Relay LCP_spring(new Relay());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping springMassTS(new TimeStepping(springDS, TiDis, OSI, LCP_spring));
    //    SP::TimeStepping StratCircuitRLCD(new TimeStepping(springDS, TiDiscRLCD));
    //    StratCircuitRLCD->insertIntegrator(OSI_RLCD);
    double h = springMassTS->timeStep();
    int N = ceil((T - t0) / h);  // Number of time steps
    int k = 0;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 11);

    // For the initial time step:

    // time
    dataPlot(k, 0) = springDS->t0();

    // velocity
    dataPlot(k, 1) = (*plasticSpring->x())(0);

    // position
    dataPlot(k, 2) = (*plasticSpring->x())(1);

    // stress
    dataPlot(k, 3) = (*plasticSpring->x())(2);

    dataPlot(k, 4) = (*InterSpring->y(0))(0);

    // diode current
    dataPlot(k, 5) = (InterSpring->getLambda(0))(0);
    double x, v, sigma, plasticRate, epElastic, ekinetic, plasticDeformation = 0.0,
                                                          plasticDissipation = 0.0;

    dataPlot(k, 6) = 0;
    dataPlot(k, 7) = 0;
    dataPlot(k, 8) = 0;
    dataPlot(k, 9) = 0;
    dataPlot(k, 10) = 0;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    // --- Time loop  ---
    for (k = 1; k < N; ++k) {
      // solve ...
      springMassTS->computeOneStep();
      InterSpring->display();
      std::cout << "x(0): " << (*plasticSpring->x())(0) << std::endl;
      std::cout << "x(1): " << (*plasticSpring->x())(1) << std::endl;
      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = springMassTS->nextTime();

      // velocity
      dataPlot(k, 1) = (*plasticSpring->x())(0);

      // position
      dataPlot(k, 2) = (*plasticSpring->x())(1);

      // stress
      dataPlot(k, 3) = (*plasticSpring->x())(2);

      dataPlot(k, 4) = (*InterSpring->y(0))(0);

      // lambda = \dot epsilon^p = Taux de deformation plastique
      dataPlot(k, 5) = (InterSpring->getLambda(0))(0);
      x = (*plasticSpring->x())(1);
      v = (*plasticSpring->x())(0);
      sigma = (*plasticSpring->x())(2);
      plasticRate = (InterSpring->getLambda(0))(0);

      ekinetic = 0.5 * m * v * v;
      plasticDeformation += plasticRate * h_step;
      plasticDissipation += sigma * plasticRate * h_step;
      epElastic = 0.5 * sigma * (x - plasticDeformation);

      dataPlot(k, 6) = ekinetic;
      dataPlot(k, 7) = epElastic;
      dataPlot(k, 8) = epElastic + ekinetic;
      dataPlot(k, 9) = plasticDissipation;
      dataPlot(k, 10) = epElastic + ekinetic + plasticDissipation;

      // transfer of state i+1 into state i and time incrementation
      springMassTS->nextStep();
    }
    // Number of time iterations
    cout << "Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << endl;
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms" << endl;

    // dataPlot (ascii) output
    ioMatrix::write("springMass.dat", "ascii", dataPlot, "noDim");

    //    double error=0.0, eps=1e-12;
    //    if((error=ioMatrix::compareRefFile(dataPlot, "CircuitRLCD.ref", eps)) >= 0.0
    //        && error > eps)
    //      return 1;

  }

  // --- Exceptions handling ---
  catch (...) {
    Siconos::exception::process();
    return 1;
  }
}
