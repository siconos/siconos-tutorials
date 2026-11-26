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

/* Strip with a circular hole in tension, Von Mises plasticity (plane stress).

   Quarter of the strip (x in [0, 0.5], y in [0, 1], hole of radius 0.25 centred at the
   origin), symmetry conditions on the bottom (v_y = 0) and left (v_x = 0) edges, top edge
   pulled at a constant velocity. Elastic, then perfectly plastic: the reaction of the top
   edge saturates near the net-section limit load sigma_y * (0.5 - 0.25).

   - SolidLinearTIDS: T3 elements, unknowns [v; sigma], lumped mass;
   - one StressLinearTIR per element (output: the physical stress of the element) with a
     VonMisesPlasticityNSL;
   - fem::integrators::MoreauJeanGOSI (theta = 1) and the condensed one-step problem
     fem::nonsmooth_formulations::Plasticity (numerics PlasticityProblem, global
     semi-smooth Newton PLASTICITY_2D_NSN).

   Same data as the "Strip_with_hole" example of pytripopfem (native T3 backend); see
   docs/devel/FEM_PLASTICITY_SCHEME.md in siconos for the scheme and the comparison.

   Outputs:
   - strip_with_hole.dat (t, u_top, R_top, max sigma_vm / sigma_y, number of plastic
     elements, max plastic strain), compared with strip_with_hole.ref;
   - strip_with_hole_elements.dat (final state per element: x, y of the centroid, sigma_xx,
     sigma_yy, sigma_xy, sigma_vm, plastic strain (3));
   - for the visualization (strip_with_hole_viz.py, pyvista): strip_with_hole_nodes.dat
     (x, y, dof of u_x, dof of u_y), strip_with_hole_triangles.dat (node indices) and the
     histories, one line per step (initial state included): strip_with_hole_u.dat
     (displacements, dof order), strip_with_hole_sigma.dat and strip_with_hole_eps_p.dat
     (physical stresses and plastic strains, 3 per element, Voigt order).
*/

#include <SiconosFEM.h>

#include <FEM-MoreauJeanGOSI.hpp>
#include <FEM-Plasticity.hpp>
#include <FENode.hpp>
#include <FiniteElementModel.hpp>
#include <MeshUtils.hpp>
#include <SiconosKernel.hpp>
#include <SolidLinearTIDS.hpp>
#include <StressLinearTIR.hpp>
#include <VonMisesPlasticityNSL.hpp>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <unordered_map>

namespace fem = siconos::mechanics::fem;

static double von_mises(double sxx, double syy, double sxy) {
  return std::sqrt(sxx * sxx + syy * syy - sxx * syy + 3. * sxy * sxy);
}

int main(int argc, char* argv[]) {
  try {
    // gmsh physical tags of mesh_data/strip_with_circular_hole_72.msh
    const int prescribed_tag = 1;  // top edge
    const int bottom_tag = 2;      // symmetry, v_y = 0
    const int left_tag = 3;        // symmetry, v_x = 0
    const int bulk_tag = 4;

    // material (MPa, mm, t/mm^3)
    const double E = 210e3, nu = 0.3, rho = 7.8e-9, sigma_y = 210.;
    const double v_top = 1e2;  // prescribed velocity of the top edge
    const double h = 1e-6, t0 = 0.;
    const int n_steps = 100;
    const double T = n_steps * h;

    // ------- Dynamical system -------
    // optional argument: another mesh of the same geometry and tags (gmsh format 2.2), e.g.
    // mesh_data/strip_with_circular_hole_graded.msh; the reference comparison is done only
    // with the default mesh
    const std::string default_mesh = "mesh_data/strip_with_circular_hole_72.msh";
    const std::string mesh_file = argc > 1 ? argv[1] : default_mesh;
    auto mesh = fem::createMeshFromGMSH2(mesh_file);
    std::map<int, const fem::Material> materials = {
        {bulk_tag, fem::Material(rho, E, nu, fem::AnalysisType2D::plane_stress)}};
    auto solid =
        std::make_shared<fem::SolidLinearTIDS>(mesh, materials, fem::MassType::lumped);
    solid->applyDirichletBoundaryConditions(bottom_tag, {1});
    solid->applyDirichletBoundaryConditions(left_tag, {0});
    solid->applyDirichletBoundaryConditions(prescribed_tag, {1}, v_top);
    auto femodel = solid->FEModel();
    auto nel = femodel->elements().size();
    std::cout << "Strip with a hole (" << mesh_file << "): " << femodel->nodes().size()
              << " nodes, " << nel << " T3 elements, " << n_steps << " steps" << std::endl;

    auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    nsds->insertDynamicalSystem(solid);

    // ------- Plasticity: one interaction per element, y = physical stress -------
    auto nslaw = std::make_shared<siconos::modeling::VonMisesPlasticityNSL>(
        sigma_y, siconos::modeling::VonMisesPlasticityNSL::Hypothesis::plane_stress);
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
    // global semi-smooth Newton: 1 iteration per elastic step, 10-30 per plastic step (the
    // condensed matrix is ill-conditioned in this quasi-static regime: NSGS would need
    // 1e4-1e5 sweeps per step, see FEM_PLASTICITY_SCHEME.md)
    auto osnspb = std::make_shared<fem::nonsmooth_formulations::Plasticity>(PLASTICITY_2D_NSN);
    osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-10;
    osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
    auto s = std::make_shared<siconos::simulation::TimeStepping>(nsds, td, osi, osnspb);

    // ------- Outputs -------
    // reaction of the top edge: R = M (v_{k+1} - v_k) / h + B^T sigma_{k+1} (theta = 1),
    // sum of the y components of the top nodes
    siconos::algebra::SiconosSparseMatrix M(solid->dimension(), solid->dimension());
    femodel->computeMassMatrix(M, materials, fem::MassType::lumped);
    std::vector<siconos::algebra::Index> top_dofs;
    siconos::algebra::Index top_dof = 0;
    for (auto& node : femodel->nodes())
      if (std::abs(node->y() - 1.0) < 1e-12) top_dofs.push_back(node->global_dof_index()[1]);
    top_dof = top_dofs.front();
    std::vector<double> scaling(nel);
    for (size_t e = 0; e < nel; e++) scaling[e] = solid->elementStressScaling(e);

    auto q = solid->q_read();
    auto v = solid->velocity_read();
    auto sigma = solid->stress_read();
    auto epsilon_p = solid->plasticDeformation_read();
    Eigen::VectorXd v_old = v;

    // mesh and histories for the visualization
    {
      std::ofstream fn("strip_with_hole_nodes.dat"), ft("strip_with_hole_triangles.dat");
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
    std::ofstream fu("strip_with_hole_u.dat"), fs("strip_with_hole_sigma.dat"),
        fp("strip_with_hole_eps_p.dat");
    fu << std::setprecision(10);
    fs << std::setprecision(10);
    fp << std::setprecision(10);
    auto write_history = [&]() {
      for (Eigen::Index i = 0; i < q.size(); i++)
        fu << q(i) << (i + 1 < q.size() ? " " : "\n");
      for (size_t e = 0; e < nel; e++)
        for (int c = 0; c < 3; c++) {
          auto i = static_cast<siconos::algebra::Index>(3 * e + c);
          const char* sep = (e + 1 < nel || c < 2) ? " " : "\n";
          fs << sigma(i) / scaling[e] << sep;
          fp << epsilon_p(i) / scaling[e] << sep;
        }
    };
    write_history();

    const siconos::algebra::Index output_size = 6;
    siconos::algebra::SiconosDenseMatrix dataPlot{n_steps + 1, output_size};
    dataPlot.setZero();

    auto start = std::chrono::system_clock::now();
    int k = 0;
    while (s->hasNextEvent()) {
      s->computeOneStep();
      k++;
      Eigen::VectorXd R = M * (Eigen::VectorXd(v) - v_old) / h +
                          solid->BMatrix().transpose() * Eigen::VectorXd(sigma);
      double R_top = 0.;
      for (auto d : top_dofs) R_top += R(d);
      double max_vm = 0., max_eps_p = 0.;
      int n_plastic = 0;
      for (size_t e = 0; e < nel; e++) {
        auto i = static_cast<siconos::algebra::Index>(3 * e);
        double vm = von_mises(sigma(i), sigma(i + 1), sigma(i + 2)) / scaling[e];  // physical
        max_vm = std::max(max_vm, vm / sigma_y);
        if (vm >= sigma_y * (1. - 1e-6)) n_plastic++;
        max_eps_p = std::max(max_eps_p, epsilon_p.segment<3>(i).norm() / scaling[e]);
      }
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(top_dof);
      dataPlot(k, 2) = R_top;
      dataPlot(k, 3) = max_vm;
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
    std::cout << "Top reaction: " << dataPlot(k, 2)
              << " (limit load sigma_y * 0.25 = " << 0.25 * sigma_y
              << "), plastic elements: " << dataPlot(k, 4) << std::endl;

    siconos::algebra::io::write("strip_with_hole.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    // final state per element
    std::ofstream fe("strip_with_hole_elements.dat");
    fe << std::setprecision(10);
    size_t e = 0;
    for (auto& element : femodel->elements()) {
      double x = 0., y = 0.;
      for (auto& node : element->nodes()) {
        x += node->x() / 3.;
        y += node->y() / 3.;
      }
      auto i = static_cast<siconos::algebra::Index>(3 * e);
      Eigen::Vector3d se = sigma.segment<3>(i) / scaling[e];
      Eigen::Vector3d ep = epsilon_p.segment<3>(i) / scaling[e];
      fe << x << " " << y << " " << se(0) << " " << se(1) << " " << se(2) << " "
         << von_mises(se(0), se(1), se(2)) << " " << ep(0) << " " << ep(1) << " " << ep(2)
         << "\n";
      e++;
    }

    // ------- Comparison with the reference (relative tolerances per column) -------
    // u_top is prescribed; R_top and the plastic strain depend on the solver tolerance
    // (1e-10); max sigma_vm / sigma_y <= 1 is compared tightly; the
    // number of plastic elements (threshold 1 - 1e-6 on sigma_vm / sigma_y) within one
    // element, since an element entering the plastic zone may cross the threshold
    // differently on another platform.
    if (mesh_file != default_mesh) return 0;
    std::string ref = "strip_with_hole.ref";
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
      double error = siconos::algebra::io::compareRefFile(dataPlot, ref, g.eps, g.columns);
      if (error < 0. || error > g.eps) failed = true;
    }
    return failed ? 1 : 0;

  } catch (...) {
    std::cerr << "Exception caught in strip_with_hole.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
