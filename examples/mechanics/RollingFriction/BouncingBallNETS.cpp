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

/*!\file BouncingBallNETS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, O. Bonnefon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/
#include <chrono>
#include "SiconosKernel.hpp"


using namespace std;
#define R_CLASS NewtonEuler5DR

class my_NewtonEulerR : public R_CLASS
{

  double _sBallRadius ;

public:

  my_NewtonEulerR(double radius): R_CLASS(), _sBallRadius(radius) { };

  virtual void computeOutput(double t, Interaction& inter, unsigned int derivativeNumber)
  {
    VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
    if(derivativeNumber == 0)
    {
      computeh(t, *DSlink[NewtonEulerR::q0], *inter.y(0));
    }
    else
    {
      R_CLASS::computeOutput(t, inter, derivativeNumber);
    }

  }
  void computeh(double time, BlockVector& q0, SiconosVector& y)
  {
    double height = fabs(q0.getValue(0)) - _sBallRadius;
    // std::cout <<"my_NewtonEulerR:: computeh _jachq" << std:: endl;
    // _jachq->display();
    y.setValue(0, height);
    _Nc->setValue(0, 1);
    _Nc->setValue(1, 0);
    _Nc->setValue(2, 0);
    _Pc1->setValue(0, height);
    _Pc1->setValue(1, q0.getValue(1));
    _Pc1->setValue(2, q0.getValue(2));

    //_Pc2->setValue(0,hpc);
    //_Pc2->setValue(1,data[q0]->getValue(1));
    //_Pc2->setValue(2,data[q0]->getValue(2));
    //printf("my_NewtonEulerR N, Pc\n");
    //_Nc->display();
    //_Pc1->display();
  }
  //ACCEPT_VISITORS();
};
TYPEDEF_SPTR(my_NewtonEulerR);





int main(int argc, char* argv[])
{
  try
  {


    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    unsigned int qDim = 7;           // degrees of freedom for the ball
    unsigned int nDim = 6;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 15.0;                  // final computation time
    double h = 0.005;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 2.0;      // initial velocity for lowest bead.
    // position_init = 0.0;      // initial position for lowest bead.
    // velocity_init = 0.0;      // initial velocity for lowest bead.
    double omega_initx = 0.0;
    double omega_initz = 1.0;// initial velocity for lowest bead.
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    double m = 1; // Ball mass
    double g = 9.81; // Gravity
    double radius = 0.1;
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(qDim));
    SP::SiconosVector v0(new SiconosVector(nDim));
    SP::SimpleMatrix I(new SimpleMatrix(3, 3));
    v0->zero();
    q0->zero();
    I->eye();
    (*q0)(0) = position_init;
    /*initial quaternion equal to (1,0,0,0)*/
    (*q0)(3) = 1.0;

    (*v0)(0) = velocity_init;
    (*v0)(3) = omega_initx;
    (*v0)(5) = omega_initz;
    // -- The dynamical system --
    SP::NewtonEulerDS ball(new NewtonEulerDS(q0, v0, m, I));

    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //

    SP::NonSmoothLaw nslaw0(new NewtonImpactRollingFrictionNSL(e, e, 0.6, 0.01, 5));

    // Version with my_NewtonEulerR()
    SP::NewtonEulerR relation0(new my_NewtonEulerR(radius));
    SP::Interaction inter(new Interaction(nslaw0, relation0));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));
    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new RollingFrictionContact());

    // -- (4) Simulation setup with (1) (2) (3)

    SP::TimeStepping s(new TimeStepping(bouncingBall, t, OSI, osnspb));

    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(10);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 27;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->twist();
    SP::SiconosVector p = ball->p(1);
    SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();

    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*q)(3);
    dataPlot(0, 5) = (*q)(4);
    dataPlot(0, 6) = (*q)(5);
    dataPlot(0, 7) = (*q)(6);

    dataPlot(0, 8) = (*v)(0);
    dataPlot(0, 9) = (*v)(1);
    dataPlot(0, 10) = (*v)(2);
    dataPlot(0, 11) = (*v)(3);
    dataPlot(0, 12) = (*v)(4);
    dataPlot(0, 13) = (*v)(5);

    dataPlot(0, 14) = (*p)(0);
    dataPlot(0, 15) = (*p)(1);
    dataPlot(0, 16) = (*p)(2);
    dataPlot(0, 17) = (*p)(3);
    dataPlot(0, 18) = (*p)(4);
    dataPlot(0, 19) = (*p)(5);



    dataPlot(0, 20) = (*lambda)(0);
    dataPlot(0, 21) = (*lambda)(1);
    dataPlot(0, 22) = (*lambda)(2);
    dataPlot(0, 23) = (*lambda)(3);
    dataPlot(0, 24) = (*lambda)(4);


    dataPlot(0, 25) = acos((*q)(3));
    dataPlot(0, 26) = relation0->contactForce()->norm2();

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;


    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    dataPlot(k, 6) = relation0->contactForce()->norm2();
    while(s->hasNextEvent() and k<=100000)
    {
      //      s->computeOneStep();
      s->advanceToEvent();
      // osnspb->setNumericsVerboseMode(true);
      // osnspb->display();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();

      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*q)(3);
      dataPlot(k, 5) = (*q)(4);
      dataPlot(k, 6) = (*q)(5);
      dataPlot(k, 7) = (*q)(6);

      dataPlot(k, 8) = (*v)(0);
      dataPlot(k, 9) = (*v)(1);
      dataPlot(k, 10) = (*v)(2);
      dataPlot(k, 11) = (*v)(3);
      dataPlot(k, 12) = (*v)(4);
      dataPlot(k, 13) = (*v)(5);

      dataPlot(k, 14) = (*p)(0);
      dataPlot(k, 15) = (*p)(1);
      dataPlot(k, 16) = (*p)(2);
      dataPlot(k, 17) = (*p)(3);
      dataPlot(k, 18) = (*p)(4);
      dataPlot(k, 19) = (*p)(5);



      dataPlot(k, 20) = (*lambda)(0);
      dataPlot(k, 21) = (*lambda)(1);
      dataPlot(k, 22) = (*lambda)(2);
      dataPlot(k, 23) = (*lambda)(3);
      dataPlot(k, 24) = (*lambda)(4);
      // std::cout << "angle " << (*q)(3) <<  " " <<  2.0*acos((*q)(3)) << std::endl;
      dataPlot(k, 25) = 2.0*acos((*q)(3));
      dataPlot(k, 26) = relation0->contactForce()->norm2();

      s->nextStep();

      k++;
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << endl;;
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << "Computation time : " << elapsed << " ms" << endl;
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("BouncingBallNETS.dat", "ascii", dataPlot, "noDim");
    //ioMatrix::write("BouncingBallNETS.ref", "ascii", dataPlot);

    // Comparison with a reference file

    double error=0.0, eps=1e-12;
    if((error=ioMatrix::compareRefFile(dataPlot, "BouncingBallNETS.ref", eps)) >= 0.0
        && error > eps)
      return 1;


  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }

}
