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
//#include "SphereNEDSPlanR.hpp"
#include "SiconosKernel.hpp"
#include <chrono>

#define WITH_PROJ
#define WITH_FC3D
using namespace std;
#ifdef WITH_FC3D
#define R_CLASS NewtonEuler3DR
#else
#define R_CLASS NewtonEuler1DR
#endif

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

    std::cout <<"my_NewtonEulerR:: computeh" << std:: endl;
    std::cout <<"q0.size() = " << q0.size() << std:: endl;
    double height = q0.getValue(0) - _sBallRadius - q0.getValue(7);
    // std::cout <<"my_NewtonEulerR:: computeh _jachq" << std:: endl;
    // _jachq->display();

    y.setValue(0, height);
    _Nc->setValue(0, 1);
    _Nc->setValue(1, 0);
    _Nc->setValue(2, 0);
    _Pc1->setValue(0, q0.getValue(0) - _sBallRadius);
    _Pc1->setValue(1, q0.getValue(1));
    _Pc1->setValue(2, q0.getValue(2));

    _Pc2->setValue(0,q0.getValue(7));
    _Pc2->setValue(1,q0.getValue(8));
    _Pc2->setValue(2,q0.getValue(9));
    //printf("my_NewtonEulerR N, Pc\n");
    _Nc->display();
    _Pc1->display();
    _Pc2->display();
    std::cout <<"my_NewtonEulerR:: computeh ends" << std:: endl;
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
    double T = 10.0;                 // final computation time
    double h = 0.001;                // time step
    double position_init = 1.0;      // initial position
    double velocity_init = 0.0;      // initial velocity
    double omega_initx = 0.0;        // initial angular velocity
    double omega_initz = 0.0;        // initial angular velocity
    double theta = 1.0;              // theta for MoreauJeanOSI integrator
    double m = 1; // Ball mass
    double g = 10.0; // Gravity
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

    SP::IndexInt bdindex(new IndexInt(3));
    (*bdindex)[0] = 0;
    (*bdindex)[1] = 3;
    (*bdindex)[2] = 5;

    // SP::SiconosVector bdPrescribedVelocity(new SiconosVector(1));
    // bdPrescribedVelocity->setValue(0,0.5);
    // SP::BoundaryCondition bd (new BoundaryCondition(bdindex,bdPrescribedVelocity));

    SP::BoundaryCondition bd(new BoundaryCondition(bdindex));
    bd->setComputePrescribedVelocityFunction("BallOnMovingPlanePlugin", "prescribedvelocity3");

    ball->setBoundaryConditions(bd);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

#ifdef WITH_FC3D
    SP::NonSmoothLaw nslaw0(new NewtonImpactFrictionNSL(e, e, 0.6, 3));
#else
    SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
#endif

    SP::NewtonEulerR relation0(new my_NewtonEulerR(radius));
    SP::Interaction inter(new Interaction(nslaw0, relation0));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));
    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
#ifdef WITH_PROJ
    SP::MoreauJeanDirectProjectionOSI OSI(new MoreauJeanDirectProjectionOSI(theta));
#else
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));
#endif
    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new GenericMechanical());
#ifdef WITH_PROJ
    SP::OneStepNSProblem osnspb_pos(new MLCPProjectOnConstraints(SICONOS_MLCP_ENUM, 1.0));
#endif
    // -- (4) Simulation setup with (1) (2) (3)
#ifdef WITH_PROJ
    SP::TimeSteppingDirectProjection s(new TimeSteppingDirectProjection(bouncingBall, t, OSI, osnspb, osnspb_pos));
    s->setProjectionMaxIteration(20);
    s->setConstraintTolUnilateral(1e-08);
    s->setConstraintTol(1e-08);
#else
    SP::TimeStepping s(new TimeStepping(bouncingBall, t, OSI, osnspb));
#endif
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(10);
    // =========================== End of model definition ===================

    // ================================= Computation =========================


    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 20;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(1);

    // SP::SiconosVector lambda = inter->lambda(1);
    // SP::SiconosVector y = inter->y(0);

    SP::SiconosVector reaction = ball->reactionToBoundaryConditions();


    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*reaction)(0);
    dataPlot(0, 5) = acos((*q)(3));
    //dataPlot(0, 6) = relation0->contactForce()->norm2();

    dataPlot(0, 7) = (*q)(0);
    dataPlot(0, 8) = (*q)(1);
    dataPlot(0, 9) = (*q)(2);
    dataPlot(0, 10) = (*q)(3);
    dataPlot(0, 11) = (*q)(4);
    dataPlot(0, 12) = (*q)(5);
    dataPlot(0, 13) = (*q)(6);

    dataPlot(0, 14) = (*v)(0);
    dataPlot(0, 15) = (*v)(1);
    dataPlot(0, 16) = (*v)(2);
    dataPlot(0, 17) = (*v)(3);
    dataPlot(0, 18) = (*v)(4);
    dataPlot(0, 19) = (*v)(5);



    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;


    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //dataPlot(k, 6) = relation0->contactForce()->norm2();
    while(s->hasNextEvent())
    {
      //std::cout << "step " << k << std::endl;
      //      s->computeOneStep();
      s->advanceToEvent();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*reaction)(0);

      /* fix issue with machine accuracy */
      if ((*q)(3) > 1.0)
      {
        dataPlot(k, 5) = acos(1.0);
      }
      else if  ((*q)(3) < - 1.0)
      {
        dataPlot(k, 5) = acos(-1.0);
      }
      else
      {
        dataPlot(k, 5) = acos((*q)(3));
      }
      // dataPlot(k, 6) = relation0->contactForce()->norm2();
      dataPlot(k, 7) = (*q)(0);
      dataPlot(k, 8) = (*q)(1);
      dataPlot(k, 9) = (*q)(2);
      dataPlot(k, 10) = (*q)(3);
      dataPlot(k, 11) = (*q)(4);
      dataPlot(k, 12) = (*q)(5);
      dataPlot(k, 13) = (*q)(6);

      dataPlot(k, 14) = (*v)(0);
      dataPlot(k, 15) = (*v)(1);
      dataPlot(k, 16) = (*v)(2);
      dataPlot(k, 17) = (*v)(3);
      dataPlot(k, 18) = (*v)(4);
      dataPlot(k, 19) = (*v)(5);

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
    ioMatrix::write("BallNewtonEuler.dat", "ascii", dataPlot, "noDim");

    // Comparison with a reference file
    cout << "====> Comparison with a reference file ..." << endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
#ifdef WITH_PROJ
    double error=0.0, eps=1e-10;
    if((error=ioMatrix::compareRefFile(dataPlot, "BallNewtonEuler-WITHPROJ.ref",
                                       eps)) >= 0.0
        && error > eps)
      return 1;
#else
    double error=0.0, eps=1e-10;
    if((error=ioMatrix::compareRefFile(dataPlot, "BallNewtonEuler.ref", eps)) >= 0.0
        && error > eps)
      return 1;
#endif

  }

  catch(...)
  {
    siconos::exception::process();
    return 1;
  }

}
