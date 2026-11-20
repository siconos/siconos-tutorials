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

/*!\file ImpactingBar.cpp
  V. Acary

  A Bar bouncing on the ground
  Simulation with a Time-Stepping scheme.
*/

#include <chrono>
#include <cmath>

#include "SiconosKernel.hpp"

using namespace std;

using Matrix = siconos::algebra::SiconosDenseMatrix;
using Vector = siconos::algebra::SiconosVector;

/** add the stiffness matrix of element elementNumber (Euler-Bernoulli beam, 2 dofs per node)
 *  to the global stiffness matrix */
static void addElementaryStiffnessMatrix(Matrix& stiffness, int elementNumber,
                                         double elementLength) {
  int dofStart = elementNumber * 2;
  int ndof_per_element = 4;

  Matrix Ke{4, 4};
  Ke.setZero();
  Ke(0, 0) = 12.;
  Ke(0, 1) = 6. * elementLength;
  Ke(0, 2) = -12.;
  Ke(0, 3) = 6. * elementLength;

  Ke(1, 1) = 4. * elementLength * elementLength;
  Ke(1, 2) = -6. * elementLength;
  Ke(1, 3) = 2. * elementLength * elementLength;

  Ke(2, 2) = 12.;
  Ke(2, 3) = -6. * elementLength;

  Ke(3, 3) = 4. * elementLength * elementLength;

  for (unsigned int i = 0; i < 4; i++) {
    for (unsigned int j = 1; j < 4; j++) Ke(j, i) = Ke(i, j);
  }

  stiffness.block(dofStart, dofStart, ndof_per_element, ndof_per_element) += Ke;
}

/** add the mass matrix (consistent or lumped) of element elementNumber to the global mass
 *  matrix */
static void addElementaryMassMatrix(Matrix& mass, int elementNumber, double elementLength,
                                    bool lumpedMass) {
  int dofStart = elementNumber * 2;
  int ndof_per_element = 4;

  Matrix Me{4, 4};
  Me.setZero();

  if (!lumpedMass) {
    Me(0, 0) = 156.;
    Me(0, 1) = 22. * elementLength;
    Me(0, 2) = 54.;
    Me(0, 3) = -13. * elementLength;

    Me(1, 1) = 4. * elementLength * elementLength;
    Me(1, 2) = 13. * elementLength;
    Me(1, 3) = -3. * elementLength * elementLength;

    Me(2, 2) = 156.;
    Me(2, 3) = -22. * elementLength;

    Me(3, 3) = 4. * elementLength * elementLength;

    for (unsigned int i = 0; i < 4; i++) {
      for (unsigned int j = 1; j < 4; j++) Me(j, i) = Me(i, j);
    }
  } else {
    double alpha = 1 / 50.;
    Me(0, 0) = 0.5 * 420.;
    Me(1, 1) = alpha * elementLength * elementLength * 420.;
    Me(2, 2) = 0.5 * 420.;
    Me(3, 3) = alpha * elementLength * elementLength * 420.;
  }

  mass.block(dofStart, dofStart, ndof_per_element, ndof_per_element) += Me;
}

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
#include "UserDefinedParameter.hpp"

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ...\n\n";

    double l = L / nElement;  // length of an element

    cout << "beam length = " << l << endl;

    int nDof = (nElement + 1) * 2;
    cout << "number of element = " << nElement << endl;
    cout << "number of dof = " << nDof << endl;

    Matrix massMatrix{nDof, nDof};
    massMatrix.setZero();
    Matrix stiffnessMatrix{nDof, nDof};
    stiffnessMatrix.setZero();
    bool lumpedMass = false;
    for (int e = 0; e < static_cast<int>(nElement); e++) {
      addElementaryStiffnessMatrix(stiffnessMatrix, e, l);
      addElementaryMassMatrix(massMatrix, e, l, lumpedMass);
    }

    massMatrix *= rho * S * l / 420.;
    stiffnessMatrix *= E * I / (l * l * l);
    cout << "massMatrix(0,0) = " << massMatrix(0, 0) << endl;
    cout << "massMatrix(0,1) = " << massMatrix(0, 1) << endl;
    cout << "massMatrix(1,0) = " << massMatrix(1, 0) << endl;
    cout << "massMatrix(1,1) = " << massMatrix(1, 1) << endl;

    // -- Initial positions and velocities --
    Vector q0{nDof};
    q0.setZero();
    Vector v0{nDof};
    v0.setZero();
    //  -- The dynamical system --
    auto beam = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        q0, v0, massMatrix, siconos::algebra::alias_t);

    // -- Set stiffness matrix --
    beam->setStiffnessMatrix(stiffnessMatrix, siconos::algebra::alias_t);

    // -- Set external forces (no weight) --
    Vector weight{nDof};
    weight.setZero();
    beam->setConstantFext(weight, siconos::algebra::alias_t);

    // -- Clamped end: zero prescribed velocity on the last two dofs --
    auto bd = std::make_shared<siconos::modeling::BoundaryCondition>(
        siconos::modeling::BoundaryCondition::Indices{
            static_cast<siconos::algebra::SiconosVector::Index>(nDof - 1),
            static_cast<siconos::algebra::SiconosVector::Index>(nDof - 2)});
    beam->setBoundaryConditions(bd);

    //  Impacting ball
    double ball_mass = 1.0;
    Matrix ballMass{1, 1};
    ballMass(0, 0) = ball_mass;

    Vector q0_ball{1};
    q0_ball(0) = position_init;
    Vector v0_ball{1};
    v0_ball(0) = velocity_init;

    // -- The dynamical system --
    auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        q0_ball, v0_ball, ballMass, siconos::algebra::alias_t);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.0;

    // Interaction beam-ball
    //
    Matrix H{1, nDof + 1};
    H.setZero();
    H(0, 0) = -1.0;
    H(0, nDof) = 1.0;

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H);

    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

    // -------------
    // --- Model ---
    // -------------
    auto impactingBeam = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    impactingBeam->insertDynamicalSystem(beam);
    impactingBeam->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    impactingBeam->link(inter, beam, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta, 0.0);
    // OSI->setIsWSymmetricDefinitePositive(true);
    // OSI->setConstraintActivationThreshold(1e-05);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto s =
        std::make_shared<siconos::simulation::TimeStepping>(impactingBeam, t, OSI, osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = floor((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 17;
    Matrix dataPlot{N, outputSize};
    dataPlot.setZero();

    auto q = beam->q_read();
    auto v = beam->velocity_read();
    auto p = beam->p_read(1);
    auto lambda = inter->lambda(1);
    auto u = inter->y(1);

    int k = 0;

    dataPlot(k, 0) = impactingBeam->t0();

    dataPlot(k, 1) = q(0);
    dataPlot(k, 2) = v(0);
    dataPlot(k, 3) = p(0);
    dataPlot(k, 4) = (*lambda)(0);
    dataPlot(k, 7) = q(nDof - 2);
    dataPlot(k, 8) = v(nDof - 2);

    dataPlot(k, 9) = q((nDof) / 2);
    dataPlot(k, 10) = v((nDof) / 2);

    Vector tmp{nDof};

    tmp = stiffnessMatrix * q;
    double potentialEnergy = 0.5 * q.dot(tmp);
    tmp = massMatrix * v;
    double kineticEnergy = 0.5 * v.dot(tmp);
    double impactEnergy = 0.0;

    dataPlot(k, 5) = potentialEnergy;
    dataPlot(k, 6) = kineticEnergy;
    dataPlot(k, 11) = impactEnergy;

    auto qBall = ball->q_read();
    auto vBall = ball->velocity_read();

    dataPlot(k, 12) = qBall(0);
    dataPlot(k, 13) = vBall(0);

    double ballKineticEnergy = 0.5 * vBall(0) * ball_mass * vBall(0);
    dataPlot(k, 16) = ballKineticEnergy;

    // --- Time loop ---
    cout << "====> Start computation ... \n\n";
    // ==== Simulation loop - Writing without explicit event handling =====
    k++;

    auto start = std::chrono::system_clock::now();
    double max_impulse = 0.0;
    while (k < N) {
      if (k % 100 == 0) std::cout << "k :" << k << "/" << N << std::endl;
      s->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(0);
      dataPlot(k, 2) = v(0);
      dataPlot(k, 3) = p(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 7) = q(nDof - 2);
      dataPlot(k, 8) = v(nDof - 2);

      dataPlot(k, 9) = q((nDof) / 2);
      dataPlot(k, 10) = v((nDof) / 2);

      tmp = stiffnessMatrix * q;
      potentialEnergy = 0.5 * q.dot(tmp);
      tmp = massMatrix * v;
      kineticEnergy = 0.5 * v.dot(tmp);

      dataPlot(k, 5) = potentialEnergy;
      dataPlot(k, 6) = kineticEnergy;

      dataPlot(k, 12) = qBall(0);
      dataPlot(k, 13) = vBall(0);
      dataPlot(k, 14) = (*u)(0);

      ballKineticEnergy = 0.5 * vBall(0) * ball_mass * vBall(0);
      dataPlot(k, 16) = ballKineticEnergy;

      double u_p_theta =
          theta * (vBall(0) - v(0)) + (1 - theta) * (dataPlot(k - 1, 13) - dataPlot(k - 1, 2));
      impactEnergy = (u_p_theta) * (*lambda)(0);
      dataPlot(k, 11) = impactEnergy;
      dataPlot(k, 15) = u_p_theta;

      if ((*lambda)(0) > 0) {
        max_impulse = fmax(max_impulse, (*lambda)(0));
      }

      s->nextStep();

      k++;
    }

    double totalImpactEnergy = 0.0;
    for (int i = 0; i < k; i++) {
      totalImpactEnergy = totalImpactEnergy + dataPlot(i, 11);
    }
    printf("totalImpactEnergy = %e\n", totalImpactEnergy);
    printf("max_impulse = %e\n", max_impulse);

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms" << endl;
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    siconos::algebra::io::write("ImpactingEBBeam.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
