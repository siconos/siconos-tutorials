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

/*
  A Ball bouncing on the ground, with friction and the Fremond impact law.
  Simulation with an event-capturing (Time-Stepping) scheme,
  Moreau-Jean integrator, with the vectorised kernel (vkernel).

  The same model as BouncingBallTS-FremondNSL.cpp, checked against its
  reference file BouncingBallTS-FremondNSL.ref.
*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <print>
#include <vector>

#include "siconos/vkernel.hpp"
#include "siconos/vkernel/config/config_builder.hpp"

namespace siconos::vkernel::config {
struct ball : model::lagrangian_ds {};
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct nslaw : model::fremond_impact_friction {};
struct relation : model::lagrangian_r<nslaw::size> {};
struct interaction : simul::interaction<nslaw, relation> {};
struct topo : simul::topology<ball, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
struct simulation : simul::time_stepping<td, osi, osnspb> {};

struct data_t : decltype(storage()
                             .with_param<"dof", 3>()
                             .with_items<simulation>()
                             .with_time_invariant<ball, "fext">()
                             .with_diagonal<ball, "mass_matrix">()
                             .with_assembled_diagonal<osi, "mass_matrix_assembled">()
                             .build_runtime_database()){};
}  // namespace siconos::vkernel::config

// the largest difference, over all the entries, between the rows and the
// reference file (first line: number of rows and columns)
static double compare_ref_file(const std::vector<std::vector<double>>& rows,
                               const char* filename) {
  std::ifstream ref(filename);
  std::size_t nrows, ncols;
  if (!(ref >> nrows >> ncols) || nrows != rows.size()) return INFINITY;
  double error = 0;
  for (const auto& row : rows) {
    if (row.size() != ncols) return INFINITY;
    for (double x : row) {
      double r;
      if (!(ref >> r)) return INFINITY;
      error = std::max(error, std::abs(x - r));
    }
  }
  return error;
}

int main() {
  namespace storage = siconos::vkernel::storage;
  namespace config = siconos::vkernel::config;
  using siconos::vkernel::storage::handle;

  // ================= Creation of the model =======================

  // User-defined main parameters
  double t0 = 0.;      // initial computation time
  double T = 10.;      // final computation time
  double h = 0.005;    // time step
  double theta = 0.5;  // theta for MoreauJeanOSI integrator
  double R = 0.1;      // Ball radius
  double m = 1;        // Ball mass
  double g = 9.81;     // Gravity

  config::data_t data;

  // -------------------------
  // --- Dynamical systems ---
  // -------------------------
  std::print("====> Model loading ...\n");

  handle ball = storage::add<config::ball>(data);
  ball.q() = {0, 1.0, 0};
  ball.velocity() = {0, 0, 1.e-01};
  ball.mass_matrix().diagonal() << m, m, 2. / 5 * m * R * R;

  // -- Set external forces (weight) --
  ball.fext() = {0., -m * g, 0.};

  // --------------------
  // --- Interactions ---
  // --------------------

  // -- nslaw --
  double e = 0.9;
  double mu = 0.5;
  handle nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;
  nslaw.mu() = mu;

  // Interaction ball-floor
  handle relation = storage::add<config::relation>(data);
  // normal gap and tangential displacement of the contact point
  Eigen::Matrix<double, 2, 3> H;
  H << 0., 1., 0., -1., 0., -R;
  relation.h_matrix() = -H;  // the vkernel lagrangian_r takes -H

  // ------------------
  // --- Simulation ---
  // ------------------
  handle simulation = storage::add<config::simulation>(data);
  simulation.one_step_integrator().theta() = theta;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = t0;
  simulation.time_discretization().h() = h;
  simulation.time_discretization().tmax() = T;

  // -- (3) one step non smooth problem
  handle fc2d = storage::add<config::fc2d>(data);
  fc2d.create();
  handle osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = fc2d;
  handle so = storage::add<siconos::vkernel::simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_NSGS);
  osnspb.options() = so;

  handle interaction = simulation.topology().link(ball);
  simulation.topology().bind(interaction, relation);
  interaction.nslaw() = nslaw;

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  simulation.initialize();

  // time, q, v(0), p(0), y, u = H v, lambda, as BouncingBallTS-FremondNSL.cpp
  std::vector<std::vector<double>> data_plot;
  auto output = [&](double p0, double lambda0, double lambda1) {
    auto step = simulation.current_step();
    Eigen::Vector3d q = storage::attr<"q">(ball, step);
    Eigen::Vector3d v = storage::attr<"velocity">(ball, step);
    Eigen::Vector2d y = H * q, u = H * v;
    // the reference has y = u = 0 at t0: the kernel computes the outputs of
    // the interaction from the first step on
    if (step == 0) y = u = Eigen::Vector2d::Zero();
    data_plot.push_back({step * simulation.time_step(), q(0), q(1), q(2), v(0), p0, y(0), y(1),
                         u(0), u(1), lambda0, lambda1});
  };
  output(0., 0., 0.);

  std::print("====> Start computation ...\n");
  while (simulation.has_next_event()) {
    auto involved = simulation.compute_one_step();
    if (involved > 0) {
      auto osi = simulation.one_step_integrator();
      auto lambda = get_vector(osi.lambda_vector_assembled(), 0);
      output(get_vector(osi.p0_vector_assembled(), 0)(0), lambda(0), lambda(1));
    } else {
      output(0., 0., 0.);
    }
  }
  std::print("End of computation - Number of iterations done: {}\n", data_plot.size() - 1);

  // --- Output files ---
  std::print("====> Output file writing ...\n");
  std::ofstream dat("BouncingBallTS-VkernelFremondNSL.dat");
  for (const auto& row : data_plot) {
    for (double x : row) std::print(dat, "{:.15e} ", x);
    std::print(dat, "\n");
  }

  double eps = 1e-12;
  double error = compare_ref_file(data_plot, "BouncingBallTS-FremondNSL.ref");
  std::print("Comparison with BouncingBallTS-FremondNSL.ref: error = {:e}\n", error);
  return error < eps ? 0 : 1;
}
