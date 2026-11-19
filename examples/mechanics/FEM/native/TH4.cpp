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

// Mesh file in C++ to read directly mesh (this file can be produced with meshio_prepost)
// #include "./mesh_data/cube_fine.cpp"
// #include "./mesh_data/beam.cpp"

int main(int argc, char* argv[]) {
  try {
    /* Mesh creation *********************************************************/
    // std::shared_ptr<Mesh> mesh (createMesh()); // to create a mesh for a cpp file
    // std::shared_ptr<Mesh> mesh (createMeshFromGMSH2("./mesh_data/tetra_simple.msh"));    //
    // simple reference tetrahedron for debugging TH4 element std::shared_ptr<Mesh> mesh
    // (createMeshFromGMSH2("./mesh_data/tetra_simple2.msh"));   // simple tetrahedron for
    // debugging TH4 element from Felippa's book std::shared_ptr<Mesh> mesh
    // (createMeshFromGMSH2("./mesh_data/cube.msh"));            // simple cube example

    auto gmsh_filename = "./mesh_data/cube_multi.msh";

    // Applied forces
    siconos::algebra::SiconosVector nodal_forces{3};
    nodal_forces << 0., 0., -1e6;
    // Boundary Conditions
    std::vector<int> node_dof_index(3);
    node_dof_index[0] = 0;
    node_dof_index[1] = 1;
    node_dof_index[2] = 2;

    siconos::mechanics::fem::Tags tags;
    tags[siconos::mechanics::fem::MeshTags::bulk_material] = 1;
    tags[siconos::mechanics::fem::MeshTags::boundary_conditions] = 2;
    tags[siconos::mechanics::fem::MeshTags::applied_forces] = 3;

    siconos::mechanics::fem::Material mat{7800, 210e9, 1. / 3};

    auto FEsolid = siconos::mechanics::fem::build_dynamicalsystem_from_gmsh(
        gmsh_filename, tags, mat, nodal_forces, node_dof_index);

    double t0 = 0;     // initial computation time
    double T = 1e-03;  // final computation time

    auto solid = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    // Contact Conditions
    auto femodel = FEsolid->FEModel();
    double e = 0.0;
    double Lz = 10.;
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    siconos::algebra::SiconosVector initial_gap{1};
    initial_gap << Lz * 1e-05;
    std::vector<siconos::algebra::SiconosDenseMatrix> Hv;
    std::cout << "contact node number : [ ";
    for (auto node : femodel->nodes()) {
      if (fabs(node->z()) <= 1e-16 and fabs(node->x()) >= 40.) {
        std::cout << " " << node->num();
        auto idx_z = node->global_dof_index()[2];
        Hv.emplace_back(1, FEsolid->dimension());
        Hv.back().setZero();
        Hv.back()(0, idx_z) = 1.0;
        auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
        auto relation =
            std::make_shared<siconos::modeling::LagrangianLinearTIR>(Hv.back(), initial_gap);
        auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);
        // link the interaction and the dynamical system
        solid->link(inter, FEsolid);
      }
    }
    std::cout << "]" << std::endl;
    //  ------------------
    //  --- Simulation ---
    //  ------------------
    double h = 1e-05;    // time step
    double theta = 1.0;  // 0.5;              // theta for MoreauJeanOSI integrator

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
    // =========================== End of model definition
    // ===========================

    // ================================= Computation
    // =================================
    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    siconos::algebra::Index outputSize = 5;
    siconos::algebra::SiconosDenseMatrix dataPlot(N + 1, outputSize);
    dataPlot.setZero();

    auto q = FEsolid->q_read();
    auto v = FEsolid->velocity_read();
    auto p = FEsolid->p_read(1);
    auto fext = FEsolid->fext();

    dataPlot(0, 0) = solid->t0();
    dataPlot(0, 1) = q(FEsolid->dimension() - 1);
    dataPlot(0, 2) = v(FEsolid->dimension() - 1);
    dataPlot(0, 3) = p(0);
    // dataPlot(0, 4) = (*lambda)(0);

    auto filename = siconos::mechanics::fem::prepareWriteDisplacementforPython("TH4");
    auto mesh = femodel->mesh();
    siconos::mechanics::fem::writeDisplacementforPython(*mesh, *femodel, q, filename);

    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    auto start = std::chrono::system_clock::now();
    while (simulation->hasNextEvent()) {
      simulation->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = simulation->nextTime();
      dataPlot(k, 1) = q(FEsolid->dimension() - 1);
      dataPlot(k, 2) = v(FEsolid->dimension() - 1);
      dataPlot(k, 3) = p(0);

      if (k % 1 == 0)
        siconos::mechanics::fem::writeDisplacementforPython(*mesh, *femodel, q, filename);

      // dataPlot(k, 4) = (*lambda)(0);
      simulation->nextStep();
      k++;
      siconos::tools::progressBar((double)k / N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("TH4.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "TH4.ref", eps)) >= eps)
      return 1;

    return 0;
  } catch (...) {
    std::cerr << "Exception caught in TH4.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
