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

#include <SiconosFEM.h>

#include <FEM-GlobalFrictionContact.hpp>
#include <FEM-MoreauJeanGOSI.hpp>
#include <SiconosKernel.hpp>
#include <chrono>
#include <cmath>  // for fabs

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;

    int nb_elements = 4;
    int dim = 2;
    siconos::algebra::SiconosVector3 coords_start{0., 0., 0.};
    siconos::algebra::SiconosVector3 coords_end{0., 4., 0.};
    auto mesh =
        siconos::mechanics::fem::createBeamMesh(coords_start, coords_end, nb_elements, dim);
    mesh->display(false);

    siconos::mechanics::fem::Tags tags;
    tags[siconos::mechanics::fem::MeshTags::bulk_material] = 1;
    tags[siconos::mechanics::fem::MeshTags::boundary_conditions] = 2;
    tags[siconos::mechanics::fem::MeshTags::applied_forces] = 3;

    siconos::mechanics::fem::Material mat{7800, 210e9, 1. / 3};

    // Same material for all tags.
    std::map<int, const siconos::mechanics::fem::Material> materials = {
        {1, mat}, {2, mat}, {3, mat}};

    // Create the beam
    auto beam = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(mesh, materials);

    // Apply nodal forces
    if (tags.find(siconos::mechanics::fem::MeshTags::applied_forces) != tags.end()) {
      siconos::algebra::SiconosVector nodal_forces;
      if (dim == 2) {
        nodal_forces.resize(3);
        nodal_forces << -2e9, 0., 0.;
      } else {
        nodal_forces.resize(6);
        nodal_forces << -2e9, 0., 0., 0., 0., 0.;
      }
      beam->applyNodalForces(tags[siconos::mechanics::fem::MeshTags::applied_forces],
                             nodal_forces);
    }

    // Boundary Conditions
    if (tags.find(siconos::mechanics::fem::MeshTags::boundary_conditions) != tags.end()) {
      std::vector<int> bc_dof_index;
      if (dim == 2) {
        bc_dof_index.resize(3);
        bc_dof_index[0] = 0;
        bc_dof_index[1] = 1;
        bc_dof_index[2] = 2;
      } else {
        bc_dof_index.resize(6);
        bc_dof_index[0] = 0;
        bc_dof_index[1] = 1;
        bc_dof_index[2] = 2;
        bc_dof_index[3] = 3;
        bc_dof_index[4] = 4;
        bc_dof_index[5] = 5;
      }
      beam->applyDirichletBoundaryConditions(
          tags[siconos::mechanics::fem::MeshTags::boundary_conditions], bc_dof_index);
    }

    // ------- NSDS -------
    double t0 = 0;     // initial computation time
    double T = 1e-02;  // final computation time
    auto beamNSDS = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    beamNSDS->insertDynamicalSystem(beam);

    // ------- Simulation -------
    // -- (1) integrator --
    double theta = 1.0;  // theta for MoreauJeanOSI integrator
    auto osi_fem =
        std::make_shared<siconos::mechanics::fem::integrators::MoreauJeanGOSI>(theta);
    osi_fem->setIsWSymmetricDefinitePositive(false);

    // -- (2) Time discretisation --
    double h = 1e-05;  // time step
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem --
    int SOLVER;
    if (dim == 2)
      SOLVER = SICONOS_FRICTION_2D_NSGS;
    else
      SOLVER = SICONOS_GLOBAL_FRICTION_3D_NSGS_WR;

    auto osnspb = std::make_shared<
        siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact>(dim, SOLVER);

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(beamNSDS, t, osi_fem, osnspb);

    // ---- Outputs ----
    int N = ceil((T - t0) / h);  // Number of time steps
    std::string SOFAfilename = "beam.state";
    int k = 1;
    siconos::algebra::Index outputSize = 4;
    siconos::algebra::SiconosDenseMatrix dataPlot{N + 1, outputSize};
    dataPlot.setZero();
    auto q = beam->q_read();
    auto v = beam->velocity_read();

    auto start = std::chrono::system_clock::now();

    siconos::mechanics::fem::prepareWriteBeamPositionforSOFA(SOFAfilename);

    siconos::mechanics::fem::writeBeamPositionforSOFA(*mesh, *beam->FEModel(), q, SOFAfilename,
                                                      0);

    siconos::algebra::SiconosDenseMatrix pos{N + 1, 10};
    pos.setZero();
    int lastXposition, lastYposition, lastOrientation;
    if (dim == 2) {
      lastXposition = beam->dimension() - 3;
      lastYposition = beam->dimension() - 2;
      lastOrientation = beam->dimension() - 1;
    } else {
      lastXposition = beam->dimension() - 6;
      lastYposition = beam->dimension() - 5;
      lastOrientation = beam->dimension() - 1;
    }

    dataPlot(0, 0) = beamNSDS->t0();
    dataPlot(0, 1) = q(lastXposition);
    dataPlot(0, 2) = q(lastYposition);
    dataPlot(0, 3) = q(lastOrientation);

    while (s->hasNextEvent()) {
      std::cout << "beam->dimension(): " << beam->dimension() << std::endl;
      s->computeOneStep();
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(lastXposition);
      dataPlot(k, 2) = q(lastYposition);
      dataPlot(k, 3) = q(lastOrientation);
      siconos::algebra::print(q);
      siconos::mechanics::fem::writeBeamPositionforSOFA(*mesh, *beam->FEModel(), q,
                                                        SOFAfilename, k * 0.01);

      s->nextStep();
      k++;

      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    siconos::algebra::io::write("bernoulliBeam-GOSI.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "bernoulliBeam.ref", eps)) >=
        eps)
      return 1;

    return 0;

  } catch (...) {
    std::cerr << "Exception caught in bernoulliBeam.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
