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
#include <SiconosKernel.hpp>
#include <chrono>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDof = 3;                // degrees of freedom for the ball
    double t0 = 0;               // initial computation time
    double T = 2.0;              // final computation time
    double h = 0.0005;           // time step
    double position_init = 1.0;  // initial position for lowest bead.
    double velocity_init = 0.0;  // initial velocity for lowest bead.
    double theta = 0.5;          // theta for MoreauJeanOSI integrator
    double R = 0.1;              // Ball radius
    double m = 1;                // Ball mass
    double g = 9.81;             // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ...\n\n";

    // Number of Beads
    unsigned int nBeads = 10;
    double initialGap = 0.25;
    double alert = 0.02;

    Matrix mass{nDof, nDof};
    mass.setZero();
    mass(0, 0) = m;
    mass(1, 1) = m;
    mass(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --

    std::vector<Vector> q0(nBeads);
    std::vector<Vector> v0(nBeads);

    for (unsigned int i = 0; i < nBeads; i++) {
      q0[i] = Vector::Zero(nDof);
      v0[i] = Vector::Zero(nDof);
      q0[i](0) = position_init + i * initialGap;
      v0[i](0) = velocity_init;
    }

    // -- The dynamical system --

    Vector weight{nDof};
    weight.setZero();
    weight(0) = -m * g;

    std::vector<std::shared_ptr<siconos::modeling::LagrangianLinearTIDS>> beads(nBeads);
    for (unsigned int i = 0; i < nBeads; i++) {
      beads[i] = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
          q0[i], v0[i], mass, siconos::algebra::alias_t);
      // -- Set external forces (weight) --
      beads[i]->setConstantFext(weight, siconos::algebra::alias_t);
    }

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //
    Matrix H{1, nDof};
    H.setZero();
    H(0, 0) = 1.;
    Vector b{1};
    b << -R;
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H, b);

    std::shared_ptr<siconos::modeling::Interaction> inter{nullptr};

    // beads/beads interactions
    Matrix HOfBeads{1, 2 * nDof};
    HOfBeads.setZero();
    HOfBeads(0, 0) = -1.0;
    HOfBeads(0, 3) = 1.0;
    Vector bOfBeads{1};
    bOfBeads << -2 * R;
    std::vector<std::shared_ptr<siconos::modeling::LagrangianLinearTIR>> relationOfBeads(
        nBeads - 1);
    std::vector<std::shared_ptr<siconos::modeling::Interaction>> interOfBeads(nBeads - 1);

    // --------------------------------------
    // ---      Model and simulation      ---
    // --------------------------------------

    auto columnOfBeads = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // add the dynamical system in the non smooth dynamical system
    for (unsigned int i = 0; i < nBeads; i++) {
      columnOfBeads->insertDynamicalSystem(beads[i]);
    }

    // --  (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto s =
        std::make_shared<siconos::simulation::TimeStepping>(columnOfBeads, t, OSI, osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 1 + nBeads * 4;
    Matrix dataPlot(N + 1, outputSize);

    dataPlot(0, 0) = columnOfBeads->t0();

    for (unsigned int i = 0; i < nBeads; i++) {
      dataPlot(0, 1 + i * 2) = (*(beads[i]->q()))(0);
      dataPlot(0, 2 + i * 2) = (*(beads[i]->velocity()))(0);
      //      dataPlot(0,3+i*4) = (*(beads[i]->p(1)))(0);
    }

    // for (unsigned int i =1; i< nBeads; i++)
    // {
    // dataPlot(0,4+i*4) = (*(interOfBeads[i-1]->lambda(1)))(0);
    // }

    // --- Time loop ---
    cout << "====> Start computation ... \n\n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();
    int ncontact = 0;
    while (s->hasNextEvent()) {
      // Rough contact detection
      for (unsigned int i = 0; i < nBeads - 1; i++) {
        // Between first bead and plane
        if (abs((*(beads[i])->q())(0) - R) < alert) {
          if (!inter) {
            ncontact++;
            // std::cout << "Number of contact = " << ncontact << std::endl;

            inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);
            columnOfBeads->link(inter, beads[0]);

            assert((*inter->y(0))(0) >= 0);
          }
        }

        // Between two beads
        if (abs((*(beads[i + 1])->q())(0) - (*((beads[i])->q()))(0) - 2 * R) < alert) {
          // std::cout << "Alert distance for declaring contact = ";
          // std::cout << (*abs((*(beads[i])->q())-(0)-((beads[i+1])->q()))(0))
          // <<std::endl;
          if (!interOfBeads[i].get()) {
            ncontact++;
            // std::cout << "Number of contact = " << ncontact << std::endl;

            relationOfBeads[i] =
                std::make_shared<siconos::modeling::LagrangianLinearTIR>(HOfBeads, bOfBeads);
            interOfBeads[i] =
                std::make_shared<siconos::modeling::Interaction>(nslaw, relationOfBeads[i]);

            columnOfBeads->link(interOfBeads[i], beads[i], beads[i + 1]);

            assert((*interOfBeads[i]->y(0))(0) >= 0);
          }
        }
      }

      s->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      for (unsigned int i = 0; i < nBeads; i++) {
        dataPlot(k, 1 + i * 2) = (*(beads[i]->q()))(0);
        dataPlot(k, 2 + i * 2) = (*(beads[i]->velocity()))(0);
      }
      // for (unsigned int i =1; i< nBeads; i++)
      // {
      //   dataPlot(k,4+i*4) = (*(interOfBeads[i-1]->lambda(1)))(0);
      // }
      // for (unsigned int i =1; i< nBeads; i++)
      // {
      //   std::cout <<  (*(interOfBeads[i-1]->y(0)))(0) << std::endl ;
      // }

      s->nextStep();

      k++;
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    siconos::algebra::io::write("ColumnOfBeadsTS.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "ColumnOfBeadsTS.ref", eps)) >
        eps)
      return 1;

    return 0;

  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
