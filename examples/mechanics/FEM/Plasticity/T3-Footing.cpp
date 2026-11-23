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
    // auto gmsh_filename = "mesh_data/Footing_03.msh";
    // auto gmsh_filename = "mesh_data/Half_footing_rough.msh";
    auto gmsh_filename = "mesh_data/half_footing_2x2_b02_size01.msh";

    auto mesh = siconos::mechanics::fem::createMeshFromGMSH2(gmsh_filename);
    // mesh->display(false);

    siconos::mechanics::fem::writeMeshforPython(mesh);

    int bulk_material_tag = 5;
    int boundary_condition_tag = 1;
    int applied_force_tag = 4;

    // std::shared_ptr<Material> mat1 = std::make_shared<Material>(1, 8*36/5.,
    // 1/5.); // material for  triangle_felippa.msh

    // double density = 7800.;
    double density = 1.5e3;  // Rho, material density
    double E = 1e7;          // Young's Modulus
    double nu = 0.3;         // Poisson's Ratio
    auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(density, E, nu);
    std::map<unsigned int, std::shared_ptr<siconos::mechanics::fem::Material>> materials = {
        {bulk_material_tag, mat1}};

    auto start = std::chrono::system_clock::now();
    std::cout << "About to create FESolid..." << std::endl;
    auto FEsolid = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(
        mesh, materials, siconos::algebra::UblasType::SPARSE);
    std::cout << "Created FESolid..." << std::endl;

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Assembly time : " << elapsed << " ms\n";

    auto femodel = FEsolid->FEModel();
    // FEsolid->K()->display();

    /*------------------------------------------------- Applied forces  */

    //    (*nodal_forces)(0) = 1e6;
    // auto nodal_forces = std::make_shared<Vector>(2);
    // nodal_forces->zero();
    // (*nodal_forces)(1) = -8e3;
    // FEsolid->applyNodalForces(applied_force_tag, nodal_forces);

    /*------------------------------------------------- Boundary Conditions  */
    /* This part should be hidden in a new BC function for a node number
     * and a dof index. */

    auto node_dof_index = std::make_shared<std::vector<int>>(0);
    node_dof_index->push_back(0);
    node_dof_index->push_back(1);

    FEsolid->applyDirichletBoundaryConditions(boundary_condition_tag, node_dof_index);
    FEsolid->applyDirichletBoundaryConditions(boundary_condition_tag + 1, node_dof_index);
    FEsolid->applyDirichletBoundaryConditions(boundary_condition_tag + 2, node_dof_index);

    auto node_dof_index_ImposedVelocity = std::make_shared<std::vector<int>>(0);
    node_dof_index_ImposedVelocity->push_back(1);

    FEsolid->applyUniformDirichletBoundaryConditions(applied_force_tag,
                                                     node_dof_index_ImposedVelocity, 1.0e-2);

    FEsolid->boundaryConditions()->display();

    // -------------
    // --- Model ---
    // -------------
    double t0 = 0;  // initial computation time
    // double T = 1e-02;    // final computation time
    double T = 2;  // final computation time
    // double h = 1e-05;    // time step
    double h = 1e-1;     // time step
    double theta = 1.0;  // theta for MoreauJeanOSI integrator

    auto solid = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    /*------------------------------------------------- Contact Conditions  */
    double e = 0.0;
    double phi = 25 * M_PI / 180;
    double c = 2e3 * cos(phi);
    // double c = 2e1*cos(phi);
    double mu = 0;
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(e, e, mu, 2);

    // auto initial_gap = std::make_shared<Vector>(1, Ly * 5e-4);
    // for (auto& n : femodel->nodes()) {
    //   if (fabs(n->y()) <= 1e-16 and fabs(n->x()) >= 1e-16) {
    //     std::cout << "contact node number : " << n->num() << " " << n->y()
    //     << "\n";
    //     auto idx_x = (*n->dofIndex())[0];
    //     auto idx_y = (*n->dofIndex())[1];
    //     auto H = std::make_shared<Matrix>(2, FEsolid->dimension());
    //     (*H)(0, idx_y) = 1.0;
    //     (*H)(1, idx_x) = 1.0;
    //     auto relation =
    //         std::make_shared<siconos::modeling::LagrangianLinearTIR>(
    //             H, initial_gap);
    //     auto inter =
    //         std::make_shared<siconos::modeling::Interaction>(nslaw, relation);
    //     // link the interaction and the dynamical system
    //     solid->link(inter, FEsolid);
    //   }
    // }

    auto nslawPlasticity =
        std::make_shared<siconos::modeling::MohrCoulombPlasticityNSL>(c, phi, 2);

    auto initialStress_gap = std::make_shared<Vector>(2, 0);
    int elcount = 0;
    for (auto& el : femodel->elements()) {
      auto H = std::make_shared<Matrix>(2, FEsolid->stressDimension());
      auto Hplus = std::make_shared<Matrix>(2, FEsolid->stressDimension());
      auto Hyy = std::make_shared<Matrix>(2, FEsolid->stressDimension());
      auto Hyyplus = std::make_shared<Matrix>(2, FEsolid->stressDimension());
      auto Hxy = std::make_shared<Matrix>(2, FEsolid->stressDimension());
      auto Hxyplus = std::make_shared<Matrix>(2, FEsolid->stressDimension());

      el->display();
      (*H)(0, elcount * 3) = -1.0;
      (*H)(1, elcount * 3 + 1) = -1.0;
      (*Hyy)(0, elcount * 3 + 1) = -1.0;
      (*Hyy)(1, elcount * 3 + 2) = -1.0;
      (*Hxy)(0, elcount * 3 + 2) = -1.0;
      (*Hxy)(1, elcount * 3) = -1.0;

      // (*H)(2, elcount*3+2) = -1.0;
      (*Hplus)(0, elcount * 3) = 1.0;
      (*Hplus)(1, elcount * 3 + 1) = 1.0;
      (*Hyyplus)(0, elcount * 3 + 1) = 1.0;
      (*Hyyplus)(1, elcount * 3 + 2) = 1.0;
      (*Hxyplus)(0, elcount * 3 + 2) = 1.0;
      (*Hxyplus)(1, elcount * 3) = 1.0;

      auto relation =
          std::make_shared<siconos::modeling::StressLinearTIR>(H, initialStress_gap);
      auto relationPlus =
          std::make_shared<siconos::modeling::StressLinearTIR>(Hplus, initialStress_gap);
      auto relationyy =
          std::make_shared<siconos::modeling::StressLinearTIR>(Hyy, initialStress_gap);
      auto relationyyPlus =
          std::make_shared<siconos::modeling::StressLinearTIR>(Hyyplus, initialStress_gap);
      auto relationxy =
          std::make_shared<siconos::modeling::StressLinearTIR>(Hxy, initialStress_gap);
      auto relationxyPlus =
          std::make_shared<siconos::modeling::StressLinearTIR>(Hxyplus, initialStress_gap);

      auto inter = std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relation);
      auto interPlus =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationPlus);
      auto interyy =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationyy);
      auto interyyPlus =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationyyPlus);
      auto interxy =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationxy);
      auto interxyPlus =
          std::make_shared<siconos::modeling::Interaction>(nslawPlasticity, relationxyPlus);

      // solid->link(inter, FEsolid);
      // solid->link(interPlus, FEsolid);
      // solid->link(interyy, FEsolid);
      // solid->link(interyyPlus, FEsolid);
      solid->link(interxy, FEsolid);
      solid->link(interxyPlus, FEsolid);

      elcount++;
      // if (elcount == 1)
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
    unsigned int outputSize = 9;
    Matrix dataPlot(N + 1, outputSize);

    auto q = FEsolid->q();
    auto v = FEsolid->velocity();
    auto sigma = FEsolid->stress();
    auto epsilon_p = FEsolid->epsilon_p(0);
    auto dotEpsilonp = FEsolid->epsilon_p(1);
    auto p = FEsolid->p(1);
    // auto lambda = inter->lambda(1);

    std::cout << "q size: " << q->size() << std::endl;
    std::cout << "dotEpsilonp size: " << dotEpsilonp->size() << std::endl;

    dataPlot(0, 0) = solid->t0();
    dataPlot(0, 1) = (*q)(647);
    dataPlot(0, 2) = (*v)(647);
    dataPlot(0, 3) = (*sigma)(0);
    dataPlot(0, 4) = (*dotEpsilonp)(0);
    dataPlot(0, 5) = (*sigma)(1);
    dataPlot(0, 6) = (*dotEpsilonp)(1);
    dataPlot(0, 7) = (*sigma)(2);
    dataPlot(0, 8) = (*dotEpsilonp)(2);
    // dataPlot(0, 4) = (*lambda)(0);

    auto filename = siconos::mechanics::fem::prepareWriteDisplacementforPython(
        "half_footing_2x2_b02_size01");
    auto filename2 = siconos::mechanics::fem::prepareWriteTensorforPython(
        "half_footing_2x2_b02_size01", "dotEpsilonp");

    // auto filename =
    //     siconos::mechanics::fem::prepareWriteDisplacementforPython("Half_footing_rough");
    writeDisplacementforPython(mesh, femodel, q, filename);
    writeTensorforPython(femodel, dotEpsilonp, filename2, "dotEpsilonp");

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
      dataPlot(k, 1) = (*q)(647);
      dataPlot(k, 2) = (*v)(647);
      dataPlot(k, 3) = (*sigma)(0);
      dataPlot(k, 4) = (*dotEpsilonp)(0);
      dataPlot(k, 5) = (*sigma)(1);
      dataPlot(k, 6) = (*dotEpsilonp)(1);
      dataPlot(k, 7) = (*sigma)(2);
      dataPlot(0, 8) = (*dotEpsilonp)(2);

      if (k % 1 == 0) writeDisplacementforPython(mesh, femodel, q, filename);
      writeTensorforPython(femodel, dotEpsilonp, filename2, "dotEpsilonp");

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

    siconos::algebra::io::write("T3-half_footing_2x2_b02_size01.dat", dataPlot,
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
