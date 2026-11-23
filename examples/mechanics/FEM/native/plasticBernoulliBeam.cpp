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
// #include <StressLinearTIR.hpp>
#include <chrono>
#include <cmath>  // for fabs

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;

    int nb_elements = 4;
    int dim = 3;
    siconos::algebra::SiconosVector3 coords_start{0., 0., 0.};
    siconos::algebra::SiconosVector3 coords_end{0., 4., 0.};

    int dimBeamDOFs = dim == 2 ? 3 : 6;
    auto mesh =
        siconos::mechanics::fem::createBeamMesh(coords_start, coords_end, nb_elements, dim);
    mesh->display(false);

    siconos::mechanics::fem::Tags tags;
    tags[siconos::mechanics::fem::MeshTags::bulk_material] = 1;
    tags[siconos::mechanics::fem::MeshTags::boundary_conditions] = 2;
    tags[siconos::mechanics::fem::MeshTags::applied_forces] = 3;

    siconos::mechanics::fem::Material mat{1080, 500e7, 1. / 3, 2.};

    // Same material for all tags.
    // Same material for all tags.
    std::map<int, const siconos::mechanics::fem::Material> materials = {
        {1, mat}, {2, mat}, {3, mat}};

    // Create the beam
    auto beam = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(mesh, materials);

    // Apply nodal forces
    if (tags.find(siconos::mechanics::fem::MeshTags::applied_forces) != tags.end()) {
      siconos::algebra::SiconosVector nodal_forces{6};
      nodal_forces << -1e9, 0., 0., 0., 0., 0.;
      // nodal_forces << 0., 0., 4e9, 0. , 0., 0. ;
      beam->applyNodalForces(tags[siconos::mechanics::fem::MeshTags::applied_forces],
                             nodal_forces);
    }
    // Boundary Conditions

    if (tags.find(siconos::mechanics::fem::MeshTags::boundary_conditions) != tags.end()) {
      std::vector<int> bc_dof_index(6);
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

    // ------- Interactions -------

    auto ndof_beam = beam->dimension();

    // ------- NSDS -------
    double t0 = 0;     // initial computation time
    double T = 1e-02;  // final computation time
    auto beamNSDS = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    beamNSDS->insertDynamicalSystem(beam);

    // unsigned int nDof = 3;  // degrees of freedom for the ball
    siconos::algebra::Index nslawdimension = 3;  // degrees of freedom for the NS law

    // ------------------------------
    // --- Plasticity Interaction ---
    // ------------------------------

    double c = 2.0e9;  // Used only to specify elastic bounds

    double phi = M_PI / 4;  // Not used at the moment

    auto nslawPlasticity =
        std::make_shared<siconos::modeling::MohrCoulombPlasticityNSL>(c, phi, nslawdimension);
    siconos::algebra::SiconosDenseMatrix H{nslawdimension, beam->stressDimension()};
    siconos::algebra::SiconosDenseMatrix Hminus{nslawdimension, beam->stressDimension()};
    siconos::algebra::SiconosDenseMatrix H2{nslawdimension, beam->stressDimension()};
    siconos::algebra::SiconosDenseMatrix Hminus2{nslawdimension, beam->stressDimension()};

    H.setZero();
    Hminus.setZero();
    int dimStressPerElement = 6;
    for (int elNum = 2; elNum < 3; elNum++) {
      H(0, elNum * dimStressPerElement + 5) = 1.0;
      H(1, elNum * dimStressPerElement) = 1.0;
      H(2, elNum * dimStressPerElement + 1) = 1.0;
      Hminus(0, elNum * dimStressPerElement + 5) = -1.0;
      Hminus(1, elNum * dimStressPerElement) = -1.0;
      Hminus(2, elNum * dimStressPerElement + 1) = -1.0;

      auto relation = std::make_shared<siconos::mechanics::fem::StressLinearTIR>(H);
      auto inter = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relation);
      auto relationMinus = std::make_shared<siconos::mechanics::fem::StressLinearTIR>(Hminus);
      auto interMinus =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationMinus);

      beamNSDS->link(inter, beam);
      beamNSDS->link(interMinus, beam);
    }

    // ---------------------------
    // --- Integration and solver Setup ---
    // ---------------------------
    // ------- Simulation -------
    // -- (1) integrator --
    double theta = 1.0;  // theta for MoreauJeanOSI integrator
    auto osi_fem =
        std::make_shared<siconos::mechanics::fem::integrators::MoreauJeanGOSI>(theta);
    // osi_fem->setIsWSymmetricDefinitePositive(true);

    // -- (2) Time discretisation --
    double h = 1e-05;  // time step
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem --
    auto osnspb = std::make_shared<
        siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact>(
        nslawdimension, SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);
    // osnspb->setNumericsVerboseLevel(1);
    auto s = std::make_shared<siconos::simulation::TimeStepping>(beamNSDS, t, osi_fem, osnspb);

    int N = ceil((T - t0) / h);  // Number of time steps
    std::string SOFAfilename = "beam";
    int k = 1;

    siconos::algebra::Index outputSize = 5;
    siconos::algebra::SiconosDenseMatrix dataPlot{N + 1, outputSize};
    dataPlot.setZero();
    auto start = std::chrono::system_clock::now();
    auto q = beam->q_read();
    auto v = beam->velocity_read();
    auto sigma = beam->stress_read();
    auto epsilonp = beam->plasticDeformation_read();
    auto dotEpsilonp = beam->plasticRate_read();
    auto femodel = beam->FEModel();
    siconos::mechanics::fem::prepareWriteBeamPositionforSOFA(SOFAfilename + ".state");
    siconos::mechanics::fem::writeBeamPositionforSOFA(*mesh, *femodel, q,
                                                      SOFAfilename + ".state", 0);
    dataPlot(0, 0) = beamNSDS->t0();
    dataPlot(0, 1) = sigma(beam->stressDimension() - 8);
    dataPlot(0, 2) = sigma(beam->stressDimension() - 7);
    dataPlot(0, 3) = epsilonp(beam->stressDimension() - 7);
    dataPlot(0, 4) = dotEpsilonp(beam->stressDimension() - 7);

    auto filename = siconos::mechanics::fem::prepareWriteDisplacementforPython("beam");

    while (s->hasNextEvent()) {
      std::cout << "Step number " << k << std::endl;
      s->computeOneStep();
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = sigma(beam->stressDimension() - 8);
      dataPlot(k, 2) = sigma(beam->stressDimension() - 7);
      dataPlot(k, 3) = epsilonp(beam->stressDimension() - 7);
      dataPlot(k, 4) = dotEpsilonp(beam->stressDimension() - 7);

      siconos::mechanics::fem::writeBeamPositionforSOFA(*mesh, *femodel, q,
                                                        SOFAfilename + ".state", k * 0.01);
      if (k % 1 == 0)
        siconos::mechanics::fem::writeDisplacementforPython(*mesh, *femodel, q, filename);
      s->nextStep();
      k++;
      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    siconos::algebra::io::write("plasticBernoulliBeam.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-6;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "plasticBernoulliBeam.ref",
                                                      eps)) >= eps)
      return 1;

    return 0;
  } catch (...) {
    std::cerr << "Exception caught in ballAgainstBernoulliBeam.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
