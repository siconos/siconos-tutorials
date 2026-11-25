/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

#include <FENode.hpp>
#include <FemTools.hpp>
#include <FiniteElementModel.hpp>
#include <Material.hpp>
#include <Mesh.hpp>
#include <MeshUtils.hpp>
#include <SiconosKernel.hpp>
#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>
#include <cmath>  // for fabs

#include "native_fem_utils.h"

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;
    //  std::shared_ptr<Mesh> mesh = create2dMesh2x1();
    //  std::shared_ptr<Mesh> mesh = create2dMeshnxm(50, 15 , 3., Ly);
    // string gmsh_filename = "./mesh_data/triangle_felippa.msh";
    // string gmsh_filename = "./mesh_data/triangle_reference.msh";
    // string gmsh_filename = "./mesh_data/square_6.msh";
    auto gmsh_filename = "./mesh_data/square_200.msh";
    // string gmsh_filename = "./mesh_data/square_2720.msh";

    // Applied forces
    siconos::algebra::SiconosVector nodal_forces{2};
    nodal_forces << 0., -1e7;
    // Boundary Conditions
    std::vector<int> node_dof_index(2);
    node_dof_index[0] = 0;
    node_dof_index[1] = 1;

    siconos::mechanics::fem::Tags tags;
    tags[siconos::mechanics::fem::MeshTags::bulk_material] = 1;
    tags[siconos::mechanics::fem::MeshTags::boundary_conditions] = 2;
    tags[siconos::mechanics::fem::MeshTags::applied_forces] = 3;

    siconos::mechanics::fem::Material mat{7800, 210e9, 1. / 3};

    auto FEsolid = siconos::mechanics::fem::build_dynamicalsystem_from_gmsh(
        gmsh_filename, tags, mat, nodal_forces, node_dof_index);
    double t0 = 0;     // initial computation time
    double T = 1e-02;  // final computation time
    auto solid = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);
    // Contact Conditions
    double e = 0.0;
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    double initial_gap = Ly * 5e-4;
    siconos::algebra::SiconosVector normal{2};
    normal << 0., 1.;

    auto contact_condition = [](const siconos::mechanics::fem::FENode& node) {
      return (std::fabs(node.y()) <= 1e-16 and std::fabs(node.x()) >= 1e-16);
    };

    auto collision_detection = std::make_shared<siconos::mechanics::fem::ContactDetection>(
        initial_gap, normal, nslaw, FEsolid, contact_condition);

    // ------------------
    // --- Simulation ---
    // ------------------
    double h = 1e-05;    // time step
    double theta = 1.0;  // theta for MoreauJeanOSI integrator

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    OSI->setIsWSymmetricDefinitePositive(true);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto simulation =
        std::make_shared<siconos::simulation::TimeStepping>(solid, t, OSI, osnspb);
    simulation->insertInteractionManager(collision_detection);

    //  Computation
    return native_fem_examples::run_T3_simulation(
        simulation, FEsolid, "T3_collision_detection", "T3_square_200.ref");
  } catch (...) {
    std::cerr << "Exception caught in T3_collision_detection.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
