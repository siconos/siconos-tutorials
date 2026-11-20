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

  This is similar to BulletBouncingBox.cpp, but adds a second box
  part-way through the simulation to demonstrate dynamic changes to
  the graph.
*/

#include <SolverOptions.h>

#include <SiconosBulletCollisionManager.hpp>
#include <SiconosCollision.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <iostream>

#include "../../toolbox/ResultChecks.hpp"
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

constexpr auto mass = 1.;

auto makeBox(Eigen::Ref<Vector> q0, Eigen::Ref<Vector> v0) {
  // -- Shape: cube with all dimensions=1.0
  auto box1 = std::make_shared<siconos::collision::SiconosBox>(1.0, 1.0, 1.0);

  // -- The dynamical system --
  auto inertia = std::make_shared<siconos::algebra::SiconosMatrix33>();
  inertia->setIdentity();
  auto body = std::make_shared<siconos::collision::RigidBodyDS>(q0, v0, mass, *inertia);

  // -- add the box to the body's set of contactactors
  // -- by default, the contactor id is 0 with no position offset,
  //    see SiconosContactor.hpp for how to change these.
  auto bodyCtrSet = std::make_shared<siconos::collision::SiconosContactorSet>();
  bodyCtrSet->push_back(std::make_shared<siconos::collision::SiconosContactor>(box1));
  body->setContactors(bodyCtrSet);
  return body;
}

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
    Vector q0Body{7};
    Vector v0Body{6};

    // -- Initial position and velocity
    // -- Initial position and velocity
    q0Body.setZero();
    v0Body.setZero();
    q0Body(2) = position_init;
    q0Body(3) = 1.0;
    v0Body(2) = velocity_init;
    // -- Set external forces (weight) --
    // -- Set external forces (weight) --
    Vector FExt{3};
    FExt.setZero();
    FExt(2) = -g * mass;

    // -- Moving object --
    auto body = makeBox(q0Body, v0Body);
    body->setConstantFext(FExt, siconos::algebra::alias_t);

    // -- Add the dynamical system in the non smooth dynamical system
    model->insertDynamicalSystem(body);

    auto ground = std::make_shared<siconos::collision::SiconosPlane>();

    // -- Create a Z-offset of -0.5 for the ground so that contact is at zero.
    Vector groundOffset{7};
    groundOffset.setZero();
    groundOffset(2) = -0.5;  // translation 0,0,-0.5
    groundOffset(3) = 1;     // orientation 1,0,0,0

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    auto timedisc = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- OneStepNsProblem --
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(3);

    // -- Some configuration

    osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;  // Max number of
    // iterations
    osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-5;  // Tolerance

    osnspb->setMaxSize(16384);  // max number of
    // interactions

    osnspb->setMStorageType(NM_SPARSE_BLOCK);  // Sparse storage

    osnspb->setNumericsVerboseMode(0);  // 0 silent, 1
    // verbose

    osnspb->setKeepLambdaAndYState(true);  // inject
    // previous
    // solution

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
    collision_manager->addStaticBody(staticCtrSet, groundOffset);

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

    // --- Data for the checks of the results (see the end of main) ---
    namespace rb = siconos::examples::rigid_body;
    std::shared_ptr<siconos::collision::RigidBodyDS> body2;  // the box added at step 100
    auto totalEnergy = [&]() {
      double E = rb::energy(body->q_read(), body->twist_read(), body->totalInertiaMatrix(), g);
      if (body2)
        E += rb::energy(body2->q_read(), body2->twist_read(), body2->totalInertiaMatrix(), g);
      return E;
    };
    const double E0 = totalEnergy();
    siconos::examples::MaxIncrease energyIncrease;
    energyIncrease.push(E0);
    double minLowestZ = rb::boxLowestZ(body->q_read(), 0.5, 0.5, 0.5);
    double minCentreDistance = 1e300;
    double addTime = -1.;
    // first impact of box 2 on the ground, first impact of box 1 (on box 2)
    double impact2Time = -1., impact1Time = -1.;
    double v1Old = body->twist_read()(2), v2Old = 0.;

    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    Vector q0{7};
    Vector v0{6};

    // -- Initial position and velocity
    // -- Initial position and velocity
    q0.setZero();
    v0.setZero();
    q0(2) = 3.;
    q0(3) = 1.0;

    auto start = std::chrono::system_clock::now();
    while (simulation->hasNextEvent()) {
      // --- Add a dynamic object at step 100 of the simulation ---
      if (k == 100) {
        auto ds = makeBox(q0, v0);
        ds->setConstantFext(FExt, siconos::algebra::alias_t);
        simulation->nonSmoothDynamicalSystem()->insertDynamicalSystem(ds);
        simulation->associate(osi, ds);
        body2 = ds;
        addTime = simulation->startingTime();
        energyIncrease.restart(totalEnergy());
      }

      simulation->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = simulation->nextTime();
      dataPlot(k, 1) = q(2);
      dataPlot(k, 2) = v(2);

      energyIncrease.push(totalEnergy());
      minLowestZ = std::min(minLowestZ, rb::boxLowestZ(body->q_read(), 0.5, 0.5, 0.5));
      const double v1 = body->twist_read()(2);
      if (body2) {
        const double v2 = body2->twist_read()(2);
        minLowestZ = std::min(minLowestZ, rb::boxLowestZ(body2->q_read(), 0.5, 0.5, 0.5));
        minCentreDistance =
            std::min(minCentreDistance, rb::centreDistance(body->q_read(), body2->q_read()));
        if (impact2Time < 0. && v2Old < 0. && v2 > 0.) {
          impact2Time = dataPlot(k, 0);
        }
        if (impact1Time < 0. && v1Old < 0. && v1 > v1Old + 1.) {
          impact1Time = dataPlot(k, 0);
        }
        v2Old = v2;
      }
      v1Old = v1;

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
    siconos::algebra::io::write("result_dynamic.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // --- Checks of the results (see BulletBouncingBox.cpp) ---
    // Box 1 falls from z = 10 at t = 0; box 2 is added at z = 3 at the beginning of step
    // 100 (t = 0.495), right below it. Box 2 lands on the ground, then box 1 on box 2.
    // The restitution of these impacts is not checked: Bullet builds the contact
    // manifold over a few steps, and at the first contact it may give the points of one
    // edge only (an edge impact: the box rotates, then the opposite edge hits).
    const double e = 0.8, groundZ = -0.5;
    const double z2 = q0(2);
    siconos::examples::Checks chk("BulletBouncingBoxDynamic");
    // - first impact of box 2 (on the ground)
    const double t2 = addTime + std::sqrt(2 * z2 / g);
    chk.near("time of the first impact of box 2", impact2Time, t2, 2 * h);
    // - first impact of box 1, on box 2 bouncing up (ideal 1D impacts): z1 - z2 = 1
    const double v2up = e * g * (t2 - addTime);
    double t1 = t2;
    while (t1 < T && (position_init - 0.5 * g * t1 * t1) -
                             (v2up * (t1 - t2) - 0.5 * g * (t1 - t2) * (t1 - t2)) >
                         1.)
      t1 += 1e-6;
    chk.near("time of the first impact of box 1 (on box 2)", impact1Time, t1, 3 * h);
    // - bounded penetration (ground and box/box), no energy increase
    const double margin = h * std::sqrt(2 * g * position_init) + 0.05;
    chk.atLeast("lowest point of the boxes", minLowestZ, groundZ - margin);
    chk.atLeast("distance between the centres", minCentreDistance, 1. - margin);
    chk.atMost("max energy increase in one step", energyIncrease.value(), 1e-5 * E0);
    // - both boxes at rest, lying on a face (they may still slide and spin: no friction)
    for (auto [name, b] : {std::pair{"box 1", body}, std::pair{"box 2", body2}}) {
      chk.atMost(std::string("final vertical velocity of ") + name,
                 std::abs(b->twist_read()(2)), 1e-2);
      chk.atLeast(std::string("final orientation of ") + name + ": on a face",
                  rb::verticalAxisCosine(b->q_read()), 0.999);
    }
    return chk.report();
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }

  return 0;
}
