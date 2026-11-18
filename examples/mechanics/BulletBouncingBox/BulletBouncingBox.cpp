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

/* A Bullet box bouncing on the ground

  A box bouncing on the ground with the use of Bullet collision
  detection.

*/

#include <SolverOptions.h>

#include <SiconosBulletCollisionManager.hpp>
#include <SiconosCollision.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <iostream>
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main() {
  // User-defined main parameters
  double t0 = 0;                // initial computation time
  double T = 20.0;              // end of computation time
  double h = 0.005;             // time step
  double position_init = 10.0;  // initial position
  double velocity_init = 0.0;   // initial velocity

  double g = 9.81;
  double theta = 0.5;  // theta for MoreauJeanOSI integrator

  // -----------------------------------------
  // --- Dynamical systems && interactions ---
  // -----------------------------------------

  try {
    // ------------
    // --- Init ---
    // ------------

    std::cout << "====> Model loading ...\n\n";

    // -- OneStepIntegrators --
    auto osi = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    // -- Model --
    auto model = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // -- Shape: cube with all dimensions=1.0
    auto box1 = std::make_shared<siconos::collision::SiconosBox>(1., 1., 1.);

    // -- Initial position and velocity
    Vector q0{7};
    Vector v0{6};
    q0.setZero();
    v0.setZero();
    q0(2) = position_init;
    q0(3) = 1.0;
    v0(2) = velocity_init;

    // -- The dynamical system --
    Matrix inertia{3, 3};
    inertia.setIdentity();
    auto body = std::make_shared<siconos::collision::RigidBodyDS>(q0, v0, 1.0, inertia);

    // -- add the box to the body's set of contactactors
    // -- by default, the contactor id is 0 with no position offset,
    //    see SiconosContactor.hpp for how to change these.
    auto bodyCtrSet = std::make_shared<siconos::collision::SiconosContactorSet>();
    bodyCtrSet->push_back(std::make_shared<siconos::collision::SiconosContactor>(box1));
    body->setContactors(bodyCtrSet);

    // -- Set external forces (weight) --
    Vector FExt{3};
    FExt.setZero();
    FExt(2) = -g * body->scalarMass();
    body->setConstantFext(FExt, siconos::algebra::alias_t);

    // -- Add the dynamical system in the non smooth dynamical system
    model->insertDynamicalSystem(body);

    auto ground = std::make_shared<siconos::collision::SiconosPlane>();

    // -- Create a Z-offset of -0.5 for the ground so that contact is at zero.
    auto groundOffset = std::make_shared<Vector>(7);
    (*groundOffset)(2) = -.5;  // translation 0,0,-0.5
    (*groundOffset)(3) = 1;    // orientation 1,0,0,0

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    auto timedisc = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- OneStepNsProblem --
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(3);

    // -- Some configuration

    // Max number of iterations
    osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    // Tolerance
    osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-5;

    osnspb->setMaxSize(16384);  // max number of interactions

    osnspb->setMStorageType(NM_SPARSE_BLOCK);  // Sparse storage

    osnspb->setNumericsVerboseMode(0);  // 0 silent, 1 verbose

    osnspb->setKeepLambdaAndYState(true);  // inject previous solution

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ...\n\n";

    int N = ceil((T - t0) / h);  // Number of time steps

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(0.8, 0., 0.0, 3);

    // some options for the Bullet collision manager:
    // -- defaults are okay, see SiconosBulletCollisionManager.hpp
    // -- in particular we want to leave multipoint iterations enabled
    //    to allow Bullet to collect more points for plane-plane collisions.
    // std::shared_ptr<siconos::collision::bullet::SiconosBulletOptions> options;

    // -- The collision manager performs broadphase collision
    //    detection, we use the Bullet implementation here.
    auto collision_manager =
        std::make_shared<siconos::collision::bullet::SiconosBulletCollisionManager>();

    // -- insert a non smooth law for contactors id 0
    collision_manager->insertNonSmoothLaw(nslaw, 0, 0);

    // -- The ground is a static object.  The collision manager
    //    maintains a list of contact sets for static objects, so we add one.
    // -- We give it a group contactor id : 0
    // -- We apply the groundOffset to the SiconosContactorSet, but
    //    equivalently it could be applied to the SiconosContactor inside the
    //    set, this is a design choice allowing for re-use for more complex
    //    compound contactor sets.
    auto staticCtrSet = std::make_shared<siconos::collision::SiconosContactorSet>();
    staticCtrSet->push_back(std::make_shared<siconos::collision::SiconosContactor>(ground));
    collision_manager->addStaticBody(staticCtrSet, *groundOffset);

    // -- MoreauJeanOSI Time Stepping with Bullet collision manager as
    // -- the interaction manager.
    auto simulation = std::make_shared<siconos::simulation::TimeStepping>(model, timedisc);
    simulation->insertInteractionManager(collision_manager);

    simulation->insertIntegrator(osi);
    simulation->insertNonSmoothProblem(osnspb);

    std::cout << "====> End of initialisation ...\n\n";

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 4;
    Matrix dataPlot(N + 1, outputSize);
    dataPlot.setZero();

    auto q = body->q_read();
    auto v = body->twist_read();

    dataPlot(0, 0) = model->t0();
    dataPlot(0, 1) = q(2);
    dataPlot(0, 2) = v(2);

    // --- Time loop ---

    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    auto start = std::chrono::system_clock::now();
    while (simulation->hasNextEvent()) {
      // while (k < 2) {

      collision_manager->resetStatistics();
      simulation->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = simulation->nextTime();
      dataPlot(k, 1) = q(2);
      dataPlot(k, 2) = v(2);

      // If broadphase collision detection shows some contacts then we may
      // display contact forces.
      if ((collision_manager->statistics().new_interactions_created +
           collision_manager->statistics().existing_interactions_processed) > 0) {
        // we *must* have an indexSet0, filled by Bullet broadphase collision
        // detection and an indexSet1, filled by TimeStepping::updateIndexSet
        // with the help of Bullet getDistance() function.
        if (model->topology()->numberOfIndexSet() == 2) {
          auto index1 = simulation->indexSet(1);

          // This is the narrow phase contact detection : if
          // TimeStepping::updateIndexSet has filled indexSet1 then we
          // have some contact forces to display
          if (index1->size() > 0) {
            // Four contact points for a cube with a side facing the
            // ground. Note : changing Bullet margin for collision
            // detection may lead this assertion to be false.
            if (index1->size() == 4) {
              auto iur = index1->begin();

              // different version of bullet may not gives the same
              // contact points! So we only keep the summation.
              dataPlot(k, 3) = index1->bundle(*iur)->lambda(1)->norm() +
                               index1->bundle(*++iur)->lambda(1)->norm() +
                               index1->bundle(*++iur)->lambda(1)->norm() +
                               index1->bundle(*++iur)->lambda(1)->norm();
            }
          }
        }
      }

      // auto &DSG =
      // simulation->nonSmoothDynamicalSystem()->dynamicalSystems();
      // std::weak_ptr<SiconosMatrix>
      // work; //{nullptr};
      // // auto
      // work{nullptr};
      // DynamicalSystemsGraph::VIterator
      // dsi, dsend; for
      // (std::tie(dsi, dsend) =
      // DSG->vertices(); dsi !=
      // dsend; ++dsi) {

      //   auto &osi2 =
      //   *DSG->properties(*dsi).osi;
      //   auto ds =
      //   DSG->bundle(*dsi); auto
      //   sods =
      //   std::static_pointer_cast<SecondOrderDS>(ds);
      //   work =
      //   static_cast<MoreauJeanOSI
      //   &>(osi2).Winverse(sods);
      //   work.lock()->display();
      //   // osi->Winverse(sods,
      //   true);
      //   //      work =
      //   static_cast<MoreauJeanGOSI
      //   &>(*osi).Winverse(sods,
      //   true);
      // }
      // static_pointer_cast<MoreauJeanOSI>(osi)->Winverse(sods,
      // true);
      simulation->nextStep();
      siconos::tools::progressBar((double)k / N);
      k++;
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1 << "\n";
    std::cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("BouncingBox.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "BouncingBox.ref", eps)) > eps)
      return 1;

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }

  return 0;
}
