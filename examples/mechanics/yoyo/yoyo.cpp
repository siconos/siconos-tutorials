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

#include <SiconosKernel.hpp>
#include <cassert>
#include <chrono>
#include <cmath>

#include "donnees.h"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;
using namespace user;

int main(int argc, char* argv[]) {
  try {
    int nDof = 3;   // nombre de degrés de liberté du yoyo
    double t0 = 0;  //  instants initial et final de la simulation
    double T = 50;
    double h = 0.001;                 // pas de discrétisation du temps
    double newton_tolerance = 1e-10;  // Newton tolerance of the time-stepping (see below)
    const double theta = 0.5;         // coefficient pour le générateur de simulation
    double e = 0;                     // coefficient de restitution
    int N = ceil((T - t0) / h) + 1;

    auto law = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation = std::make_shared<siconos::modeling::LagrangianRheonomousR>();

    relation->setComputehFunction([](const siconos::algebra::BlockVector& q, double time,
                                     Eigen::Ref<siconos::algebra::MapVectorType> y) {
      y(0) = q(1) - r * q(0) + L - q(2);
    });

    relation->setComputeContactMatrixFunction(
        [](const siconos::algebra::BlockVector& pos, double time,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(0, 0) = -r;
          result(0, 1) = 1;
          result(0, 2) = -1;
        });

    Matrix H{1, nDof};
    H.setZero();
    auto law0 = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    auto relation0 = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H);

    unsigned int outputSize = 9;
    Matrix dataPlot(N, outputSize);

    Vector q0{nDof};
    q0.setZero();
    Vector v0{nDof};
    v0.setZero();

    q0(0) = L / (2 * r);     //  vlaeur de  teta( 0 )
    q0(1) = -L + r * q0(0);  // y ( 0 )
                             // valeur initiale de h
                             // vitesses initiels de teta , y et  h
    v0(1) = r * v0(0);

    // Objectifs du joueurs avec les instants correspendants
    Matrix Controle(G + 1, 2);
    Controle(0, 1) = q0(1) - q0(2);
    for (int i = 0; i < G; i++) {
      Controle(i + 1, 0) = times[i];
      Controle(i + 1, 1) = Som[i] - L;
    }

    int k = 0;

    Matrix mass_phase1{nDof, nDof};
    mass_phase1.setIdentity();
    mass_phase1(0, 0) = I + m * r * r;
    mass_phase1(0, 2) = m * r;
    mass_phase1(1, 0) = -r;
    mass_phase1(1, 1) = 1;
    mass_phase1(1, 2) = -1;
    Matrix mass_phase2{nDof, nDof};
    mass_phase2.setIdentity();
    mass_phase2(1, 1) = m;
    mass_phase2(0, 0) = I;

    Vector fext1{nDof};
    fext1.setZero();
    fext1(0) = -m * r * g;

    Vector fext2{nDof};
    fext2.setZero();
    fext2(1) = -m * g;

    std::shared_ptr<siconos::integrators::MoreauJeanOSI> OSI = nullptr;

    std::shared_ptr<Vector> q, v;
    auto start = std::chrono::system_clock::now();
    while (k < N) {
      ///////////////////////////////////////Phase contrainte
      /////////////////////////////////////////

      if (k != 0) {
        t0 = dataPlot(k - 1, 0);
        q0(0) = (*q)(0);
        q0(2) = (*q)(2);
        q0(1) = q0(2) - L + r * q0(0);
        v0(0) = (*v)(0);
        v0(2) = (*v)(2);
        v0(1) = r * v0(0) + v0(2);
      }

      // création et insertion  du système dynamique représentant la yoyo dans le récipient
      // allDS
      auto yoyo =
          std::make_shared<siconos::modeling::LagrangianDS>(q0, v0, siconos::algebra::alias_t);
      yoyo->setConstantMass(mass_phase1, siconos::algebra::alias_t);
      yoyo->setConstantFext(fext1, siconos::algebra::alias_t);
      // yoyo->setComputeFextFunction(
      //     [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
      //       fExt(0) = 0;
      //       fExt(1) = -m * g;
      //       fExt(2) = accelerationmain (5,A,Cy,time);
      //     });

      yoyo->setComputeFintFunction(
          [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
             Eigen::Ref<siconos::algebra::MapVectorType> fint) {
            fint(0) = r * epsilon * (velocity(0));
            fint(1) = 0;
            int i = 0;
            while (user::times[i] < time) i++;
            if (velocity(0) < 0 && q(0) < thetaset(Som, i))
              fint(2) = -w * g;
            else
              fint(2) = c1 * velocity(2) + c2 * q(2);
            // fInt[2] =0;
          });

      yoyo->setComputeJacobianFintOver_qFunction(
          [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
             Eigen::Ref<siconos::algebra::MapType> jacob) {
            jacob(0, 0) = 0;
            jacob(1, 0) = 0;
            int i = 0;
            while (user::times[i] < time) i++;
            if (velocity(0) < 0 && q(0) < thetaset(Som, i))
              jacob(2, 0) = 0;
            else
              jacob(2, 0) = c2;
            // jacob[2] =0;});
          });

      yoyo->setComputeJacobianFintOver_velocityFunction(
          [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
             Eigen::Ref<siconos::algebra::MapType> jacob) {
            jacob(0, 0) = r * epsilon;
            jacob(1, 0) = 0;
            int i = 0;
            while (user::times[i] < time) i++;
            if (velocity(0) < 0 && q(0) < thetaset(Som, i))
              jacob(2, 0) = 0;
            else
              jacob(2, 0) = c1;
          });

      ////////////////  loi d'impact et relations /////////////////////////////////

      auto inter = std::make_shared<siconos::modeling::Interaction>(law0, relation0);

      /////////////////////////  MODEL //////////////////////////////////////////////////
      auto jeu = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
      jeu->insertDynamicalSystem(yoyo);
      jeu->link(inter, yoyo);
      ///////////////////// SIMULATION /////////////////////////////////

      // déscrétisation du temps
      auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

      // -- OneStepIntegrators --
      OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
      // -- OneStepNsProblem --
      auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

      auto s = std::make_shared<siconos::simulation::TimeStepping>(jeu, t, OSI, osnspb);
      // Newton tolerance of the time-stepping (default 1e-6) tightened to 1e-10.
      // The internal force switches (fint depends on the sign of v(0) and on
      // q(0) < thetaset): at the steps where the phase changes (e.g. t = 9.23, the string is
      // fully unwound and the model is re-initialized), Newton needs several iterations and,
      // with the default tolerance, stopped at a residual of ~4e-7. The solution was then
      // only accurate to that level, and a round-off difference of 1e-13 in the state
      // (macOS arm64 / clang with fused multiply-add vs Linux x86_64 / gcc) changed the
      // velocity by 1.3e-7, then the whole trajectory (2.6e-7). With 1e-10, the results no
      // longer depend on the platform at this level (reference regenerated).
      s->setNewtonTolerance(newton_tolerance);

      q = yoyo->q();
      v = yoyo->velocity();
      auto p = yoyo->p(1);
      auto lambda = inter->lambda(1);

      // --- sauver les valeurs dans une matrice dataPlot

      dataPlot(k, 0) = jeu->t0();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q)(1);
      dataPlot(k, 6) = (*v)(1);
      dataPlot(k, 7) = L + (*q)(1) - r * (*q)(0) - (*q)(2);
      dataPlot(k, 8) = (*q)(1) - (*q)(2);
      k++;

      while (s->hasNextEvent() && (*q)(0) > 0.0) {
        s->computeOneStep();
        // --- Get values to be plotted ---
        dataPlot(k, 0) = s->nextTime();
        dataPlot(k, 1) = (*q)(0);
        dataPlot(k, 2) = (*v)(0);
        dataPlot(k, 3) = (*p)(0);
        dataPlot(k, 4) = (*lambda)(0);  // multiplicateur de lagrange
        dataPlot(k, 5) = (*q)(1);
        dataPlot(k, 6) = (*v)(1);
        dataPlot(k, 7) = L + (*q)(1) - r * (*q)(0) - (*q)(2);  // contrainte géométrique
        dataPlot(k, 8) = (*q)(1) - (*q)(2);
        if (abs((*v)(0)) <= 0.05)
          std::cout << "valeur max de theta est : " << (*q)(0) << std::endl;
        k++;
        s->nextStep();
      }

      ////////////////////////// Phase libre //////////////////////

      t0 = dataPlot(k - 1, 0) + h;
      if (t0 + h < T) {
        q0(0) = 0;
        q0(2) = (*q)(2);
        q0(1) = q0(2) - L;
        v0(0) = -(*v)(0);
        v0(2) = (*v)(2);
        v0(1) = -r * v0(0) + v0(2);

        yoyo = std::make_shared<siconos::modeling::LagrangianDS>(q0, v0,
                                                                 siconos::algebra::alias_t);
        yoyo->setConstantMass(mass_phase2, siconos::algebra::alias_t);
        yoyo->setConstantFext(fext2, siconos::algebra::alias_t);

        yoyo->setComputeFintFunction(
            [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
               const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
               Eigen::Ref<siconos::algebra::MapVectorType> fint) {
              fint(0) = 0.;
              fint(1) = 0;
              int i = 0;
              while (user::times[i] < time) i++;
              if (velocity(0) < 0 && q(0) < thetaset(Som, i))
                fint(2) = -w * g;
              else
                fint(2) = c1 * velocity(2) + c2 * q(2);
              // fInt[2] =0;
            });

        yoyo->setComputeJacobianFintOver_qFunction(
            [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
               const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
               Eigen::Ref<siconos::algebra::MapType> jacob) {
              jacob(0, 0) = 0;
              jacob(1, 0) = 0;
              int i = 0;
              while (user::times[i] < time) i++;
              if (velocity(0) < 0 && q(0) < thetaset(Som, i))
                jacob(2, 0) = 0;
              else
                jacob(2, 0) = c2;
              // jacob[2] =0;});
            });

        yoyo->setComputeJacobianFintOver_velocityFunction(
            [](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
               const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time,
               Eigen::Ref<siconos::algebra::MapType> jacob) {
              jacob(0, 0) = 0.;
              jacob(1, 0) = 0.;
              int i = 0;
              while (user::times[i] < time) i++;
              if (velocity(0) < 0 && q(0) < thetaset(Som, i))
                jacob(2, 0) = 0;
              else
                jacob(2, 0) = c1;
            });

        inter = std::make_shared<siconos::modeling::Interaction>(law, relation);

        jeu = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
        jeu->insertDynamicalSystem(yoyo);
        jeu->link(inter, yoyo);

        t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
        OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
        osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();
        s = std::make_shared<siconos::simulation::TimeStepping>(jeu, t, OSI, osnspb);
        s->setNewtonTolerance(newton_tolerance);  // see above

        q = yoyo->q();
        v = yoyo->velocity();
        p = yoyo->p(1);
        lambda = inter->lambda(1);

        while (s->hasNextEvent()) {
          s->computeOneStep();
          dataPlot(k, 0) = s->nextTime();
          dataPlot(k, 1) = (*q)(0);
          dataPlot(k, 2) = (*v)(0);
          dataPlot(k, 3) = (*p)(0);
          dataPlot(k, 4) = (*lambda)(0);
          dataPlot(k, 5) = (*q)(1);
          dataPlot(k, 6) = (*v)(1);
          dataPlot(k, 7) = L + (*q)(1) - r * (*q)(0) - (*q)(2);
          dataPlot(k, 8) = (*q)(1) - (*q)(2);
          k++;
          if ((*lambda)(0) > 0 && (-r * (*v)(0) + (*v)(1) - (*v)(2)) < 10e-14) break;
          s->nextStep();
        }
      }
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done : " << k - 1 << std::endl;
    std::cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("yoyo.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("fichier.dat", Controle, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // relative tolerance: the angular velocity reaches ~80 rad/s; with the Newton tolerance
    // 1e-10, macOS and Linux results differ by up to 5e-11 (6e-13 relative).
    double error = 0.0, eps = 1e-11 * dataPlot.cwiseAbs().maxCoeff();
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "yoyo.ref", eps)) > eps)
      return 1;

    return 0;
    // --- Libérer de la mémoire
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
