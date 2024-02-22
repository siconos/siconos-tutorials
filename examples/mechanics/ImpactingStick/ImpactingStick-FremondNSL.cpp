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

/*!\file ImpactingStick.cpp
*/


#include <SiconosKernel.hpp>
#include <SolverOptions.h>
#include <chrono>
#include <cmath>



using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 0.2;                  // final computation time
    double h = 1e-4;                // time step
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    
    double l = 1.0; // length of the stick
    double m = 1.;  // mass
    double inertia= m*l*l/12. ;
    double g = 10.; // Gravity

    double pi = acos(-1);
    double initial_angle=pi/4.0;

    
    
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." <<  endl;

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = inertia;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0)(0) = 0.5 * l * sin(initial_angle);
    (*q0)(1) = 0.5 * l * cos(initial_angle)+0.01;
    (*q0)(2) = initial_angle;

    
    (*v0)(0) = -0.5;
    (*v0)(1) = 0.1;
    (*v0)(2) = 0.1;
    
    // -- The dynamical system --
    SP::LagrangianDS stick(new LagrangianDS(q0, v0, Mass));


    
    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(1) = -m * g;
    stick->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 1.0;
    double mu = 0.01;
    // Interaction stick-floor
    //

    SP::NonSmoothLaw nslaw(new FremondImpactFrictionNSL(e,0.0,mu,2));

    SP::Relation relation(new LagrangianScleronomousR("ImpactingStickPlugin:h1", "ImpactingStickPlugin:G1"));

    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem bouncingStick(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    bouncingStick->insertDynamicalSystem(stick);

    // link the interaction and the dynamical system
    bouncingStick->link(inter, stick);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));
 

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::FrictionContact osnspb(new FrictionContact(2));
    osnspb->setMStorageType(NM_SPARSE);
    osnspb->setAssemblyType(REDUCED_DIRECT);
    osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-10;
    
    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(bouncingStick, t, OSI, osnspb));

    //s->setNewtonOptions(SICONOS_TS_LINEAR);
    //OSI->setGamma(3/2.0);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================


    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 22;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SiconosVector& q = *stick->q();
    SiconosVector& v = *stick->velocity();
    SiconosVector& p = *stick->p(1);
    SiconosVector& lambda = *inter->lambda(1);
    SiconosVector& y = *inter->y(0);
    SiconosVector& u = *inter->y(1);


    dataPlot(0, 0) = bouncingStick->t0();
    dataPlot(0, 1) = q(0);
    dataPlot(0, 2) = q(1);
    dataPlot(0, 3) = q(2);
    
    dataPlot(0, 4) = v(0);
    dataPlot(0, 5) = v(1);
    dataPlot(0, 6) = v(2);

    dataPlot(0, 7) = p(0);
    dataPlot(0, 8) = y(0);
    dataPlot(0, 9) = y(1);
    dataPlot(0, 10) = u(0);
    dataPlot(0, 11) = u(1);
    dataPlot(0, 12) = lambda(0);
    dataPlot(0, 13) = lambda(1);

    dataPlot(0, 14) = 0.0;
    dataPlot(0, 15) = 0.0;
    dataPlot(0, 16) = 0.0;
    dataPlot(0, 17) = 0.0;
    dataPlot(0, 18) = 0.0;
    dataPlot(0, 19) = 0.0;
    dataPlot(0, 20) = 0.0;
    dataPlot(0, 21) = 0.0;
	


    
    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    while(s->hasNextEvent())
    {
      s->computeOneStep();
      
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = q(0);
      dataPlot(k, 2) = q(1);
      dataPlot(k, 3) = q(2);
      
      dataPlot(k, 4) = v(0);
      dataPlot(k, 5) = v(1);
      dataPlot(k, 6) = v(2);
      
      dataPlot(k, 7) = p(0);
      dataPlot(k, 8) = y(0);
      dataPlot(k, 9) = y(1);
      dataPlot(k, 10) = u(0);
      dataPlot(k, 11) = u(1);
      dataPlot(k, 12) = lambda(0);
      dataPlot(k, 13) = lambda(1);

      // compute work of contact forces
      const SiconosVector& u_k = inter->y_k(1);
      double normal_work= 0.5* (u(0)+u_k(0))*lambda(0);

      double tangential_work= 0.5* (u(1)+u_k(1))*lambda(1);
      double dissipation= 0.5* mu * fabs((u(1)+u_k(1)))*lambda(0);
      dataPlot(k, 14) = normal_work;
      dataPlot(k, 15) = tangential_work;
      dataPlot(k, 16) = dissipation;
      
      dataPlot(k, 17) = normal_work + dataPlot(k-1, 17);
      dataPlot(k, 18) = tangential_work + dataPlot(k-1, 18);
      
      // compute kinetic energy and potential energy
      double kinetic_energy = 0.5*m*(v(0)*v(0)+v(1)*v(1)) + 0.5*inertia*(v(2)*v(2));  
      double potential_energy = m*g* q(1);
      dataPlot(k, 19) = kinetic_energy;
      dataPlot(k, 20) = potential_energy;

      SP::SiconosMatrix wf = OSI->computeWorkForces();
      dataPlot(k, 21) = (*wf)(0,1) + dataPlot(k-1, 21);
      
      // if (lambda(0) >0)
      // 	{
      // 	  std::cout << s->nextTime() << std::endl;
      // 	  std::cout << "OSNSP iteration : "
      // 		<< osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_ITER_DONE]
      // 		<< " error:  : "
      // 		<< osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_RESIDU]
      // 		<< std::endl;
      // 	  osnspb->display();
      // 	  std::cout << lambda(0) << std::endl;
      // 	  std::cout << u(0) << std::endl;
      // 	  std::cout << u_k(0) << std::endl;
      // 	  std::cout << u(0)+u_k(0) << std::endl;
      // 	  //getchar();
      // 	}
      
      s->nextStep();
      k++;
    }
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << endl <<  "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms" << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("ImpactingStick-FremondNSL.dat", "ascii", dataPlot, "noDim");
    //ioMatrix::write("ImpactingStick-FremondNSL.ref", "ascii", dataPlot);
    double error=0.0, eps=1e-12;
    if((error=ioMatrix::compareRefFile(dataPlot, "ImpactingStick-FremondNSL.ref", eps)) >= 0.0
        && error > eps)
      return 1;

  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }



}
