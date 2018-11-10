/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
//  LCS for a CircuitRLCD  : sample of an electrical circuit involving :
//  - a linear dynamical system consisting of an LC oscillator (1 µF , 10 mH)
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

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  double t0 = 0.0;
  double T = 5e-3;        // Total simulation time
  double h_step = 10.0e-6;// Time step
  double Lvalue = 1e-2;   // inductance
  double Cvalue = 1e-6;   // capacitance
  double Rvalue = 1e3;    // resistance
  double Vinit = 10.0;    // initial voltage
  string Modeltitle = "LCS for a CircuitRLCD";

  try
  {
    // --- Dynamical system specification ---
    SP::SiconosVector x0(new SiconosVector(2));
    x0->setValue(0, Vinit);
    x0->setValue(1, 0.0);

    SP::SimpleMatrix A(new SimpleMatrix(2, 2));
    A->setValue(0 , 1, -1.0 / Cvalue);
    A->setValue(1 , 0, 1.0 / Lvalue);

    // --- Interaction between linear system and non smooth system ---
    SP::SimpleMatrix C(new SimpleMatrix(1, 2));
    C->setValue(0 , 0 , -1.0);

    SP::SimpleMatrix D(new SimpleMatrix(1, 1));
    D->setValue(0 , 0, Rvalue);

    SP::SimpleMatrix B(new SimpleMatrix(2, 1));
    B->setValue(0 , 0, -1.0 / Cvalue);

    SP::SiconosVector a;
    SP::SiconosVector b;

    // --- Model creation ---
    SP::NonSmoothDynamicalSystem CircuitRLCD(new NonSmoothDynamicalSystem(t0, T));

    SP::LinearComplementaritySystemsNSDS lcs(new LinearComplementaritySystemsNSDS(t0,T, x0, A, B, C, D, a, b));
    // assert(lcs->interaction());
    // assert(lcs->relation());
    // assert(lcs->ds());
    // assert(lcs->nslaw());
    // lcs->interaction()->display();
    lcs->interaction()->computeOutput(t0,0);
    lcs->interaction()->computeInput(t0,0);

    //lcs->display();

    // ------------------
    // --- Simulation ---
    // ------------------
    double theta = 0.5000000000001;

    // -- (1) OneStepIntegrators --
    SP::EulerMoreauOSI osi(new EulerMoreauOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation td(new TimeDiscretisation(t0, h_step));
    // --- (3) one step non smooth problem
    SP::LCP lcp(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping simulation(new TimeStepping(lcs, td, osi, lcp));
    double h = simulation->timeStep();
    int N = ceil((T - t0) / h); // Number of time steps
    int k = 0;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 6);

    // For the initial time step:

    // time
    dataPlot(k, 0) = lcs->t0();

    // inductor voltage
    dataPlot(k, 1) = (*lcs->ds()->x())(0);

    // inductor current
    dataPlot(k, 2) = (*lcs->ds()->x())(1);

    // diode voltage
    dataPlot(k, 3) = - (*lcs->interaction()->y(0))(0);

    // diode current
    dataPlot(k, 4) = (lcs->interaction()->getLambda(0))(0);

    dataPlot(k, 5) = (lcs->ds()->getR())(0);

    boost::timer t;
    t.restart();

    // --- Time loop  ---
    for (k = 1 ; k < N ; ++k)
    {
      // solve ...
      simulation->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = simulation->nextTime();;

      // inductor voltage
      dataPlot(k, 1) = (*lcs->ds()->x())(0);

      // inductor current
      dataPlot(k, 2) = (*lcs->ds()->x())(1);

      // diode voltage
      dataPlot(k, 3) = - (*lcs->interaction()->y(0))(0);

      // diode current
      dataPlot(k, 4) = (lcs->interaction()->getLambda(0))(0);

      // input in ds
      dataPlot(k, 5) = (lcs->ds()->getR())(0);

      // transfer of state i+1 into state i and time incrementation
      simulation->nextStep();

    }
    // Number of time iterations
    cout << "Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << t.elapsed()  << endl;

    // dataPlot (ascii) output
    ioMatrix::write("LinearComplementaritySystemsNSDS.dat", "ascii", dataPlot, "noDim");

    double error=0.0, eps=1e-12;
    if ((error=ioMatrix::compareRefFile(dataPlot, "LinearComplementaritySystemsNSDS.ref", eps)) >= 0.0
        && error > eps)
      return 1;

  }


  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cerr << e.report() << endl;
    return 1;
  }
  catch (...)
  {
    cerr << "Exception caught " << endl;
    return 1;
  }
}
