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

/* Strip footing on a cohesive-frictional soil: Mohr-Coulomb (Drucker-Prager) plasticity,
   plane strain.

   Half of the problem (symmetry): square soil domain [0, 2] x [0, 2], smooth rigid footing
   of half-width 0.2 on the top edge, pushed down at a constant velocity (displacement
   control); bottom edge v_y = 0, lateral edges v_x = 0. Weightless soil: the bearing
   capacity of Prandtl is q_u = c N_c (N_c = 20.72 for phi = 25 degrees), i.e. 8288 for the
   half footing; the coarse T3 mesh and the small domain overestimate it.

   - SolidLinearTIDS: T3 elements, unknowns [v; sigma], lumped mass, plane strain;
   - one StressLinearTIR per element (output: the physical stress of the element) with a
     DruckerPragerPlasticityNSL::mohrCoulombPlaneStrain(c, phi, psi): R <= -p sin(phi) +
     c cos(phi), R the radius of the Mohr circle, p the in-plane mean stress, dilatancy
     angle psi (psi = phi: associated flow);
   - fem::integrators::MoreauJeanGOSI (theta = 1) and the condensed one-step problem
     fem::nonsmooth_formulations::Plasticity (numerics PlasticityProblem, global
     semi-smooth Newton PLASTICITY_2D_NSN).

   Same data as the "Footing" example of pytripopfem (half_footing.py, native T3 backend):
   the results agree to round-off (docs/devel/FEM_PLASTICITY_SCHEME.md in siconos).

   Usage: footing [mesh [dilatancy angle in degrees]] (default: mesh_data/
   half_footing_2x2_b02_size01.msh, psi = phi; the reference comparison is done for the
   default case only).

   Outputs:
   - footing.dat (t, u_y of the footing, reaction of the footing, max of the yield function
     f / (c cos phi) (<= 0, 0 on the yield surface), number of plastic elements, max plastic
     strain), compared with footing.ref;
   - for the visualization (footing_viz.py, pyvista): footing_nodes.dat (x, y, dof of u_x,
     dof of u_y), footing_triangles.dat (node indices) and the histories, one line per step
     (initial state included): footing_u.dat (displacements, dof order), footing_sigma.dat
     and footing_eps_p.dat (physical stresses and plastic strains, 3 per element, Voigt
     order).
*/

#include <SiconosFEM.h>

#include <DruckerPragerPlasticityNSL.hpp>
#include <FEM-MoreauJeanGOSI.hpp>
#include <FEM-Plasticity.hpp>
#include <FENode.hpp>
#include <FiniteElementModel.hpp>
#include <MeshUtils.hpp>
#include <SiconosKernel.hpp>
#include <SolidLinearTIDS.hpp>
#include <StressLinearTIR.hpp>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <string>
#include <unordered_map>

namespace fem = siconos::mechanics::fem;

int main(int argc, char* argv[]) {
  try {
    // gmsh physical tags of mesh_data/half_footing_2x2_b02_size01.msh
    const int bottom_tag = 1;   // "Dirichlet BC1": v_y = 0
    const int right_tag = 2;    // "Dirichlet BC2": v_x = 0
    const int left_tag = 3;     // "Dirichlet BC3": symmetry, v_x = 0
    const int footing_tag = 4;  // "Forces": footing, prescribed v_y
    const int bulk_tag = 5;

    // soil (SI units): Mohr-Coulomb, c = 2 kPa, phi = 25 degrees
    const double E = 1e7, nu = 0.3, rho = 1.5e3;
    const double c = 2e3, phi = 25. * std::numbers::pi / 180.;
    const double v_footing = -1e-2;  // prescribed velocity of the footing
    const double h = 5e-2, t0 = 0.;
    const int n_steps = 40;
    const double T = n_steps * h;

    const std::string default_mesh = "mesh_data/half_footing_2x2_b02_size01.msh";
    const std::string mesh_file = argc > 1 ? argv[1] : default_mesh;
    const double psi = argc > 2 ? std::atof(argv[2]) * std::numbers::pi / 180. : phi;
    const bool default_case = (mesh_file == default_mesh) && (argc <= 2);

    // ------- Dynamical system -------
    auto mesh = fem::createMeshFromGMSH2(mesh_file);
    std::map<int, const fem::Material> materials = {
        {bulk_tag, fem::Material(rho, E, nu, fem::AnalysisType2D::plane_strain)}};
    auto solid =
        std::make_shared<fem::SolidLinearTIDS>(mesh, materials, fem::MassType::lumped);
    solid->applyDirichletBoundaryConditions(bottom_tag, {1});
    solid->applyDirichletBoundaryConditions(right_tag, {0});
    solid->applyDirichletBoundaryConditions(left_tag, {0});
    solid->applyDirichletBoundaryConditions(footing_tag, {1}, v_footing);
    auto femodel = solid->FEModel();
    auto nel = femodel->elements().size();
    std::cout << "Footing (" << mesh_file << "): " << femodel->nodes().size() << " nodes, "
              << nel << " T3 elements, dilatancy angle " << psi * 180. / std::numbers::pi
              << " degrees, " << n_steps << " steps" << std::endl;

    auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    nsds->insertDynamicalSystem(solid);

    // ------- Plasticity: one interaction per element, y = physical stress -------
    auto nslaw =
        siconos::modeling::DruckerPragerPlasticityNSL::mohrCoulombPlaneStrain(c, phi, psi);
    std::vector<siconos::algebra::SiconosMatrix> C(nel);
    for (size_t e = 0; e < nel; e++) {
      C[e] = solid->elementStressSelection(e);
      auto relation = std::make_shared<fem::StressLinearTIR>(C[e]);
      auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);
      nsds->link(inter, solid);
    }

    // ------- Simulation -------
    auto osi = std::make_shared<fem::integrators::MoreauJeanGOSI>(1.0);
    auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // global semi-smooth Newton: ~10 iterations per step (NSGS would need more than 2e4
    // sweeps per step on this ill-conditioned problem)
    auto osnspb = std::make_shared<fem::nonsmooth_formulations::Plasticity>(PLASTICITY_2D_NSN);
    osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-10;
    osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
    auto s = std::make_shared<siconos::simulation::TimeStepping>(nsds, td, osi, osnspb);

    // ------- Outputs -------
    // reaction of the footing: R = M (v_{k+1} - v_k) / h + B^T sigma_{k+1} (theta = 1), sum
    // of the y components of the footing nodes
    siconos::algebra::SiconosSparseMatrix M(solid->dimension(), solid->dimension());
    femodel->computeMassMatrix(M, materials, fem::MassType::lumped);
    std::vector<siconos::algebra::Index> footing_dofs;
    siconos::algebra::Index footing_dof = 0;
    for (auto& node : femodel->nodes())
      if (std::abs(node->y() - 2.0) < 1e-12 && node->x() <= 0.2 + 1e-12)
        footing_dofs.push_back(node->global_dof_index()[1]);
    footing_dof = footing_dofs.front();
    std::vector<double> scaling(nel);
    for (size_t e = 0; e < nel; e++) scaling[e] = solid->elementStressScaling(e);

    auto q = solid->q_read();
    auto v = solid->velocity_read();
    auto sigma = solid->stress_read();
    auto epsilon_p = solid->plasticDeformation_read();
    Eigen::VectorXd v_old = v;

    // mesh and histories for the visualization
    {
      std::ofstream fn("footing_nodes.dat"), ft("footing_triangles.dat");
      fn << std::setprecision(17);
      std::unordered_map<const fem::FENode*, size_t> node_index;
      size_t n = 0;
      for (auto& node : femodel->nodes()) {
        node_index[node.get()] = n++;
        fn << node->x() << " " << node->y() << " " << node->global_dof_index()[0] << " "
           << node->global_dof_index()[1] << "\n";
      }
      for (auto& element : femodel->elements()) {
        const auto& nodes = element->nodes();
        ft << node_index.at(nodes[0].get()) << " " << node_index.at(nodes[1].get()) << " "
           << node_index.at(nodes[2].get()) << "\n";
      }
    }
    std::ofstream fu("footing_u.dat"), fs("footing_sigma.dat"), fp("footing_eps_p.dat");
    fu << std::setprecision(10);
    fs << std::setprecision(10);
    fp << std::setprecision(10);
    auto write_history = [&]() {
      for (Eigen::Index i = 0; i < q.size(); i++)
        fu << q(i) << (i + 1 < q.size() ? " " : "\n");
      for (size_t e = 0; e < nel; e++)
        for (int k = 0; k < 3; k++) {
          auto i = static_cast<siconos::algebra::Index>(3 * e + k);
          const char* sep = (e + 1 < nel || k < 2) ? " " : "\n";
          fs << sigma(i) / scaling[e] << sep;
          fp << epsilon_p(i) / scaling[e] << sep;
        }
    };
    write_history();

    const siconos::algebra::Index output_size = 6;
    siconos::algebra::SiconosDenseMatrix dataPlot{n_steps + 1, output_size};
    dataPlot.setZero();
    dataPlot(0, 3) = -1.;
    const double c_cos_phi = c * std::cos(phi);

    auto start = std::chrono::system_clock::now();
    int k = 0;
    while (s->hasNextEvent()) {
      s->computeOneStep();
      k++;
      Eigen::VectorXd R = M * (Eigen::VectorXd(v) - v_old) / h +
                          solid->BMatrix().transpose() * Eigen::VectorXd(sigma);
      double R_footing = 0.;
      for (auto d : footing_dofs) R_footing += R(d);
      // Mohr-Coulomb yield function f = R + p sin(phi) - c cos(phi) of each element
      double max_f = -INFINITY, max_eps_p = 0.;
      int n_plastic = 0;
      for (size_t e = 0; e < nel; e++) {
        auto i = static_cast<siconos::algebra::Index>(3 * e);
        Eigen::Vector3d se = sigma.segment<3>(i) / scaling[e];  // physical
        double p = 0.5 * (se(0) + se(1));
        double radius = std::hypot(0.5 * (se(0) - se(1)), se(2));
        double f = (radius + p * std::sin(phi) - c_cos_phi) / c_cos_phi;
        max_f = std::max(max_f, f);
        if (f > -1e-6) n_plastic++;
        max_eps_p = std::max(max_eps_p, epsilon_p.segment<3>(i).norm() / scaling[e]);
      }
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(footing_dof);
      dataPlot(k, 2) = R_footing;
      dataPlot(k, 3) = max_f;
      dataPlot(k, 4) = n_plastic;
      dataPlot(k, 5) = max_eps_p;
      write_history();
      v_old = v;
      s->nextStep();
      siconos::tools::progressBar((double)k / n_steps);
    }
    auto end = std::chrono::system_clock::now();
    std::cout << "\nEnd of computation - Number of steps done: " << k << "\nComputation time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;
    std::cout << "Footing reaction: " << dataPlot(k, 2)
              << " (Prandtl, half footing: -8288), plastic elements: " << dataPlot(k, 4)
              << std::endl;

    siconos::algebra::io::write("footing.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // ------- Comparison with the reference (relative tolerances per column) -------
    // the footing displacement is prescribed; the reaction and the plastic strain depend on
    // the solver tolerance (1e-10); the yield function is <= 0 (tight tolerance); the
    // number of plastic elements (threshold -1e-6 on f) within one element.
    if (!default_case) return 0;
    struct Group {
      std::vector<int> columns;
      double eps;
    };
    const std::vector<Group> groups = {
        {{0, 1}, 1e-12 * std::max(T, std::abs(dataPlot(k, 1)))},
        {{2}, 1e-6 * dataPlot.col(2).cwiseAbs().maxCoeff()},
        {{3}, 1e-8},
        {{4}, 1.5},
        {{5}, 1e-5 * std::max(1e-12, dataPlot.col(5).cwiseAbs().maxCoeff())}};
    bool failed = false;
    for (const auto& g : groups) {
      double error =
          siconos::algebra::io::compareRefFile(dataPlot, "footing.ref", g.eps, g.columns);
      if (error < 0. || error > g.eps) failed = true;
    }
    return failed ? 1 : 0;

  } catch (...) {
    std::cerr << "Exception caught in footing.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
