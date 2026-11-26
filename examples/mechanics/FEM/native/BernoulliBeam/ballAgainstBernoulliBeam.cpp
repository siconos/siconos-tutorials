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
#include <LagrangianSparseLinearTIDS.hpp>
#include <SiconosKernel.hpp>
#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>
#include <StorageTools.hpp>
#include <chrono>
#include <cmath>  // for fabs

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;

    int nb_elements = 4;
    int dim = 3;
    siconos::algebra::SiconosVector3 coords_start{0., 0., 0.};
    siconos::algebra::SiconosVector3 coords_end{0., 4., 0.};
    auto mesh =
        siconos::mechanics::fem::createBeamMesh(coords_start, coords_end, nb_elements, 3);
    mesh->display(false);

    siconos::mechanics::fem::Tags tags;
    tags[siconos::mechanics::fem::MeshTags::bulk_material] = 1;
    tags[siconos::mechanics::fem::MeshTags::boundary_conditions] = 2;
    tags[siconos::mechanics::fem::MeshTags::applied_forces] = 3;

    siconos::mechanics::fem::Material mat{1080, 500e7, 1. / 3};

    // Same material for all tags.
    std::map<int, const siconos::mechanics::fem::Material> materials = {
        {1, mat}, {2, mat}, {3, mat}};

    // Create the beam
    auto beam = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(mesh, materials);

    // Boundary Conditions
    if (tags.find(siconos::mechanics::fem::MeshTags::boundary_conditions) != tags.end()) {
      std::vector<int> bc_dof_index;
      bc_dof_index.resize(6);
      bc_dof_index[0] = 0;
      bc_dof_index[1] = 1;
      bc_dof_index[2] = 2;
      bc_dof_index[3] = 3;
      bc_dof_index[4] = 4;
      bc_dof_index[5] = 5;
      beam->applyDirichletBoundaryConditions(
          tags[siconos::mechanics::fem::MeshTags::boundary_conditions], bc_dof_index);
    }

    std::cout << beam->stressDimension() << "\n";
    // siconos::mechanics::fem::writeMeshforPython(mesh);

    // ------- Now the block -------

    double R = 3;     // Ball radius
    double m1 = 300;  // Ball mass
    siconos::algebra::SiconosSparseMatrix mass{3, 3};
    mass.insert(0, 0) = m1;
    mass.insert(1, 1) = m1;
    mass.insert(2, 2) = 2. / 5 * m1 * R * R;
    mass.makeCompressed();

    // -- Initial positions and velocities --
    double position_init = 0.5;      // initial position for the block.
    double velocity_init = -1000.0;  // initial velocity for the block.
    siconos::algebra::SiconosVector3 q0;
    siconos::algebra::SiconosVector3 v0;
    q0 << position_init, 4., 0.;
    v0 << velocity_init, 0., 0.;

    auto block = std::make_shared<siconos::modeling::LagrangianSparseLinearTIDS>(
        q0, v0, mass, siconos::algebra::copy_t);

    // ------- Interactions -------

    auto ndof_beam = beam->dimension();
    auto ndof_stressBeam = beam->stressDimension();
    auto ndof_ball = block->dimension();
    std::cout << " beam->dimension(): " << beam->dimension() << std::endl;
    std::cout << " block->dimension(): " << block->dimension() << std::endl;
    // -- nslaw --
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(0, 0, 0, 3);
    // Interaction ball-floor
    siconos::algebra::SiconosDenseMatrix H_bb{ndof_ball, ndof_ball + ndof_beam};
    H_bb.setZero();
    H_bb(0, 24) = -1.0;
    H_bb(1, 25) = -1.0;
    H_bb(2, 26) = -1.0;
    H_bb(0, 30) = 1.0;
    H_bb(1, 31) = 1.0;
    H_bb(2, 32) = 1.0;

    siconos::algebra::SiconosVector initialGap{3};
    initialGap(0) = 0;
    initialGap(1) = -4.0;
    initialGap(0) = 0;
    auto relation_bb =
        std::make_shared<siconos::modeling::LagrangianLinearTIR>(H_bb, initialGap);
    auto inter_bb = std::make_shared<siconos::modeling::Interaction>(nslaw, relation_bb);

    // ------- NSDS -------
    double t0 = 0;     // initial computation time
    double T = 2e-02;  // final computation time
    auto beamNSDS = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    beamNSDS->insertDynamicalSystem(beam);
    beamNSDS->insertDynamicalSystem(block);
    beamNSDS->link(inter_bb, beam, block);

    // ------- Simulation -------
    // -- (1) integrator --
    double theta = 1.0;  // theta for MoreauJeanOSI integrator
    auto osi_fem =
        std::make_shared<siconos::mechanics::fem::integrators::MoreauJeanGOSI>(theta);
    auto osi_block = std::make_shared<siconos::integrators::MoreauJeanGOSI>(theta);
    osi_fem->setIsWSymmetricDefinitePositive(false);

    // -- (2) Time discretisation --
    double h = 1e-05;  // time step
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem --
    auto osnspb = std::make_shared<
        siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact>(
        3, SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(beamNSDS, t, osi_fem, osnspb);
    s->insertIntegrator(osi_block);
    s->associate(osi_block, block);
    s->associate(osi_fem, beam);
    // ---- Outputs ----
    int N = ceil((T - t0) / h);  // Number of time steps
    std::string SOFAfilename = "beam";
    int k = 1;

    siconos::algebra::Index outputSize = 6;
    siconos::algebra::SiconosDenseMatrix dataPlot{N + 1, outputSize};
    auto qball = block->q_read();
    auto q = beam->q_read();
    auto v = beam->velocity_read();
    auto p1 = block->p_read(1);

    auto start = std::chrono::system_clock::now();

    siconos::mechanics::fem::prepareWriteBeamPositionforSOFA(SOFAfilename + ".state");

    siconos::mechanics::fem::writeBeamPositionforSOFA(*mesh, *beam->FEModel(), q,
                                                      SOFAfilename + ".state", 0);
    siconos::mechanics::fem::prepareWriteBlockPositionforSOFA(SOFAfilename + "_block.state");
    siconos::mechanics::fem::writeBlockPositionforSOFA(q0, SOFAfilename + "_block.state", 0);

    siconos::algebra::SiconosDenseMatrix pos{N + 1, 10};
    int vcnt = 0;
    dataPlot(0, 0) = beamNSDS->t0();
    dataPlot(0, 1) = q(beam->dimension() - 3);
    dataPlot(0, 2) = q(beam->dimension() - 2);
    dataPlot(0, 3) = q(beam->dimension() - 1);
    dataPlot(0, 4) = qball(0);
    dataPlot(0, 5) = p1(0);

    while (s->hasNextEvent()) {
      std::cout << "Step number " << k << std::endl;
      s->computeOneStep();
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(beam->dimension() - 6);
      dataPlot(k, 2) = q(beam->dimension() - 5);
      dataPlot(k, 3) = q(beam->dimension() - 1);
      dataPlot(k, 4) = qball(0);
      dataPlot(k, 5) = p1(0);
      siconos::algebra::print(q);
      siconos::algebra::print(qball);
      siconos::mechanics::fem::writeBeamPositionforSOFA(*mesh, *beam->FEModel(), q,
                                                        SOFAfilename + ".state", k * 0.01);
      siconos::mechanics::fem::writeBlockPositionforSOFA(qball, SOFAfilename + "_block.state",
                                                         k * 0.01);
      s->nextStep();
      k++;

      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    siconos::algebra::io::write("ballAgainstBernoulliBeam.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // columns 0-4 are of order 1, column 5 of order 1e5: relative tolerances
    double eps = 1e-10;
    double eps_5 = eps * std::max(1.0, dataPlot.col(5).cwiseAbs().maxCoeff());
    double error = siconos::algebra::io::compareRefFile(
        dataPlot, "ballAgainstBernoulliBeam.ref", eps, {0, 1, 2, 3, 4});
    double error_5 = siconos::algebra::io::compareRefFile(
        dataPlot, "ballAgainstBernoulliBeam.ref", eps_5, {5});
    if (error < 0. || error >= eps || error_5 < 0. || error_5 >= eps_5) return 1;

    return 0;
  } catch (...) {
    std::cerr << "Exception caught in ballAgainstBernoulliBeam.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
