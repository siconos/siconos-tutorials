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

// =============================== Double Pendulum Example ===============================
//
// Author: Vincent Acary
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#include <SiconosKernel.hpp>
#include <chrono>

using namespace std;

double gravity = 10.0;
double m1 = 1.0;
double l1 = 1.0;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 1;  // degrees of freedom for robot arm
    double t0 = 0;          // initial computation time
    double T = 50.0;        // final computation time
    double h = 0.0005;      // time step
    double criterion = 0.00005;
    unsigned int maxIter = 2000;
    double e = 0.9;  // nslaw

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;

    // --- DS: Double Pendulum ---

    // Initial position (angles in radian)
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0).zero();
    (*v0).zero();
    (*q0)(0) = 1;

    SP::LagrangianDS simplependulum(new LagrangianDS(q0, v0));

    SP::SimpleMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m1 * l1;
    simplependulum->setMassPtr(Mass);

    // external plug-in
    // simplependulum->setComputeMassFunction("SimplePendulumPlugin","mass");

    simplependulum->setComputeFIntFunction("SimplePendulumPlugin", "FInt");
    simplependulum->setComputeJacobianFIntqDotFunction("SimplePendulumPlugin",
                                                       "jacobianVFInt");
    simplependulum->setComputeJacobianFIntqFunction("SimplePendulumPlugin", "jacobianFIntq");

    // -------------------
    // --- Interactions---
    // -------------------

    // -- relations --

    //     SimpleMatrix H(1,2);
    //     SiconosVector b(1);
    //     H.zero();
    //     H(0,0) =1.0;
    //     H(0,1) =0.0;

    //     b(0) = 0.0;

    //     NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    //     Relation relation(new LagrangianLinearTIR(H,b));
    //     Interaction inter =  new Interaction("floor-mass1", allDS,1,1, nslaw, relation);)

    string G = "SimplePendulumPlugin:G0";
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianScleronomousR("SimplePendulumPlugin:h0", G));
    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------

    SP::NonSmoothDynamicalSystem Pendulum(new NonSmoothDynamicalSystem(t0, T));
    Pendulum->insertDynamicalSystem(simplependulum);
    Pendulum->link(inter, simplependulum);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(Pendulum, t));
    s->setNewtonTolerance(criterion);
    s->setNewtonMaxIteration(maxIter);

    // -- OneStepIntegrators --

    // double theta=0.500001;
    double theta = 0.500001;

    SP::OneStepIntegrator OSI(new MoreauJeanOSI(theta));
    s->insertIntegrator(OSI);

    SP::OneStepNSProblem osnspb(new LCP());
    s->insertNonSmoothProblem(osnspb);
    cout << "=== End of model loading === " << endl;

    // --- Simulation initialization ---
    int k = 0;
    int N = ceil((T - t0) / h);
    std::cout << "Number of time step   " << N << "\n";
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    dataPlot(k, 0) = Pendulum->t0();
    dataPlot(k, 1) = (*simplependulum->q())(0);
    dataPlot(k, 2) = (*simplependulum->velocity())(0);
    dataPlot(k, 3) = l1 * sin((*simplependulum->q())(0));
    dataPlot(k, 4) = -l1 * cos((*simplependulum->q())(0));
    dataPlot(k, 5) = l1 * cos((*simplependulum->q())(0)) * ((*simplependulum->velocity())(0));
    // --- Compute elapsed time ---
    auto start = std::chrono::system_clock::now();
    //    EventsManager eventsManager = s->eventsManager();
    // --- Time loop ---
    std::cout << "Start computation ... \n";
    std::cout << "Number of time steps " << N << "\n";
    while (s->hasNextEvent()) {
      k++;
      if (!(div(k, 10000).rem)) std::cout << "Step number " << k << "\n";

      // Solve problem
      s->advanceToEvent();
      // Data Output
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*simplependulum->q())(0);
      dataPlot(k, 2) = (*simplependulum->velocity())(0);
      dataPlot(k, 3) = l1 * sin((*simplependulum->q())(0));
      dataPlot(k, 4) = -l1 * cos((*simplependulum->q())(0));
      dataPlot(k, 5) =
          l1 * cos((*simplependulum->q())(0)) * ((*simplependulum->velocity())(0));
      s->nextStep();
      progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    ioMatrix::write("SimplePendulumResult.dat", "ascii", dataPlot, "noDim");

    double error = 0.0, eps = 1e-12;
    if ((error = ioMatrix::compareRefFile(dataPlot, "SimplePendulumResult.ref", eps)) >=
            0.0 &&
        error > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
