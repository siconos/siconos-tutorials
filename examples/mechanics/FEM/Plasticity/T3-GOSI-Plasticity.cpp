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

#include <SiconosKernel.hpp>
#include <chrono>
#include <cmath>  // for fabs

#include "FENode.hpp"
#include "FiniteElement.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"
#include "MeshUtils.hpp"
#include "SolidLinearTIDS.hpp"
#include "StressLinearTIR.hpp"

using namespace std;
using namespace siconos::mechanics::fem;
using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;
    //  std::shared_ptr<Mesh> mesh = create2dMesh2x1();
    //  std::shared_ptr<Mesh> mesh = create2dMeshnxm(50, 15 , 3., Ly);
    // string gmsh_filename = "./mesh_data/triangle_felippa.msh";
    // string gmsh_filename = "./mesh_data/triangle_reference.msh";
    // string gmsh_filename = "./mesh_data/square_6.msh";
    auto gmsh_filename = "../native/mesh_data/square_200.msh";
    // string gmsh_filename = "./mesh_data/square_2720.msh";

    auto mesh = siconos::mechanics::fem::createMeshFromGMSH2(gmsh_filename);
    // mesh->display(false);

    siconos::mechanics::fem::writeMeshforPython(mesh);

    int bulk_material_tag = 1;
    int boundary_condition_tag = 2;
    int applied_force_tag = 3;

    // std::shared_ptr<Material> mat1 = std::make_shared<Material>(1, 8*36/5.,
    // 1/5.); // material for  triangle_felippa.msh
    double density = 7800.;
    auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(density, 210e9, 1 / 3.);
    std::map<unsigned int, std::shared_ptr<siconos::mechanics::fem::Material>> materials = {
        {bulk_material_tag, mat1}};

    auto start = std::chrono::system_clock::now();

    auto FEsolid = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(
        mesh, materials, siconos::algebra::UblasType::SPARSE);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Assembly time : " << elapsed << " ms\n";

    auto femodel = FEsolid->FEModel();
    // FEsolid->K()->display();

    /*------------------------------------------------- Applied forces  */

    //    (*nodal_forces)(0) = 1e6;
    auto nodal_forces = std::make_shared<Vector>(2);
    nodal_forces->zero();
    (*nodal_forces)(1) = -1e7;
    FEsolid->applyNodalForces(applied_force_tag, nodal_forces);

    /*------------------------------------------------- Boundary Conditions  */
    /* This part should be hidden in a new BC function for a node number
     * and a dof index. */

    auto node_dof_index = std::make_shared<std::vector<int>>(0);
    node_dof_index->push_back(0);
    node_dof_index->push_back(1);

    FEsolid->applyDirichletBoundaryConditions(boundary_condition_tag, node_dof_index);
    FEsolid->boundaryConditions()->display();

    // -------------
    // --- Model ---
    // -------------
    double t0 = 0;       // initial computation time
    double T = 1e-02;    // final computation time
    double h = 1e-05;    // time step
    double theta = 1.0;  // theta for MoreauJeanOSI integrator

    auto solid = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    /*------------------------------------------------- Contact Conditions  */
    double e = 0.0;

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(e, e, 0, 2);

    auto initial_gap = std::make_shared<Vector>(1, Ly * 5e-4);
    for (auto& n : femodel->nodes()) {
      if (fabs(n->y()) <= 1e-16 and fabs(n->x()) >= 1e-16) {
        std::cout << "contact node number : " << n->num() << " " << n->y() << "\n";
        auto idx_x = (*n->dofIndex())[0];
        auto idx_y = (*n->dofIndex())[1];
        auto H = std::make_shared<Matrix>(2, FEsolid->dimension());
        (*H)(0, idx_y) = 1.0;
        (*H)(1, idx_x) = 1.0;
        auto relation =
            std::make_shared<siconos::modeling::LagrangianLinearTIR>(H, initial_gap);
        auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);
        // link the interaction and the dynamical system
        solid->link(inter, FEsolid);
      }
    }

    double c = 1.0e12;
    double phi = 45;
    auto nslawPlasticity =
        std::make_shared<siconos::modeling::MohrCoulombPlasticityNSL>(c, phi, 2);

    auto initialStress_gap = std::make_shared<Vector>(2, 0);
    int elcount = 0;
    for (auto& el : femodel->elements()) {
      auto H = std::make_shared<Matrix>(2, FEsolid->stressDimension());
      el->display();
      (*H)(0, elcount * 3) = -1.0;
      (*H)(1, elcount * 3 + 1) = -1.0;
      // (*H)(2, elcount*3+2) = -1.0;
      auto relation =
          std::make_shared<siconos::modeling::StressLinearTIR>(H, initialStress_gap);
      auto inter = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relation);
      solid->link(inter, FEsolid);

      elcount++;
      // if (elcount == 2)
      //   break;
    }

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanGOSI>(theta);
    OSI->setIsWSymmetricDefinitePositive(true);

    // -- (2) Time discretisation --

    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    // auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GlobalFrictionContact>(
        2, SICONOS_FRICTION_2D_NSGS);
    // osnspb->numericsSolverOptions()->

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(solid, t, OSI, osnspb);
    // =========================== End of model definition
    // ===========================

    // ================================= Computation
    // =================================
    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 5;
    Matrix dataPlot(N + 1, outputSize);

    auto q = FEsolid->q();
    auto v = FEsolid->velocity();
    auto sigma = FEsolid->stress();
    auto epsilon_p = FEsolid->epsilon_p(1);
    auto p = FEsolid->p(1);
    // auto lambda = inter->lambda(1);
    dataPlot(0, 0) = solid->t0();
    dataPlot(0, 1) = (*q)(FEsolid->dimension() - 1);
    dataPlot(0, 2) = (*v)(FEsolid->dimension() - 1);
    dataPlot(0, 3) = (*sigma)(30);
    dataPlot(0, 4) = (*epsilon_p)(30);
    // dataPlot(0, 4) = (*lambda)(0);

    auto filename =
        siconos::mechanics::fem::prepareWriteDisplacementforPython("T3-GOSI-Plasticity");
    writeDisplacementforPython(mesh, femodel, q, filename);

    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    start = std::chrono::system_clock::now();

    while (s->hasNextEvent()) {
      s->computeOneStep();
      // osnspb->display();
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(FEsolid->dimension() - 1);
      dataPlot(k, 2) = (*v)(FEsolid->dimension() - 1);
      dataPlot(k, 3) = (*sigma)(30);
      dataPlot(k, 4) = (*epsilon_p)(30);

      if (k % 1 == 0) writeDisplacementforPython(mesh, femodel, q, filename);
      // dataPlot(k, 4) = (*lambda)(0);
      s->nextStep();
      k++;

      siconos::tools::progressBar((double)k / N);
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    dataPlot.resize(k, outputSize);

    siconos::algebra::io::write("T3-GOSI-Plasticity.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "T3_square_200.ref", eps)) >=
        eps)
      return 1;

    return 0;

  } catch (...) {
    std::cerr << "Exception caught in T3-GOSI-Plasticity.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
