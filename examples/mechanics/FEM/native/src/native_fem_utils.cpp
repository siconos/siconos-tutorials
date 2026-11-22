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

#include "native_fem_utils.h"

#include <FiniteElementModel.hpp>
#include <MeshUtils.hpp>
#include <SiconosKernel.hpp>
#include <Tools.hpp>
#include <memory>

int native_fem_examples::run_T3_simulation(
    std::shared_ptr<siconos::simulation::TimeStepping> simulation,
    std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS> FEsolid,
    std::string basename, std::string reference_file_name) {
  double time_step = simulation->timeStep();

  auto solid = simulation->nonSmoothDynamicalSystem();
  double T = solid->finalT();
  double t0 = solid->t0();
  int N = ceil((T - t0) / time_step);  // Number of time steps

  // --- Get the values to be plotted ---
  // -> saved in a matrix dataPlot
  siconos::algebra::Index outputSize = 6;
  siconos::algebra::SiconosDenseMatrix dataPlot(N + 1, outputSize);

  auto q = FEsolid->q_read();
  auto v = FEsolid->velocity_read();
  auto p = FEsolid->p_read(1);
  // auto lambda = inter->lambda(1);
  dataPlot(0, 0) = solid->t0();
  dataPlot(0, 1) = q(FEsolid->dimension() - 1);
  dataPlot(0, 2) = v(FEsolid->dimension() - 1);
  dataPlot(0, 3) = p(0);
  dataPlot(0, 4) = q(169);
  dataPlot(0, 5) = v(169);

  auto filename = siconos::mechanics::fem::prepareWriteDisplacementforPython(basename);
  auto femodel = FEsolid->FEModel();
  auto mesh = femodel->mesh();
  siconos::mechanics::fem::writeDisplacementforPython(*mesh, *femodel, q, filename);
  // --- Time loop ---
  std::cout << "====> Start computation ... \n";
  // ==== Simulation loop - Writing without explicit event handling =====
  int k = 1;
  auto start = std::chrono::system_clock::now();
  while (simulation->hasNextEvent()) {
    simulation->computeOneStep();

    // siconos::algebra::print(*osnspb);
    // --- Get values to be plotted ---
    dataPlot(k, 0) = simulation->nextTime();
    dataPlot(k, 1) = q(FEsolid->dimension() - 1);
    dataPlot(k, 2) = v(FEsolid->dimension() - 1);
    dataPlot(k, 3) = p(0);
    dataPlot(k, 4) = q(169);
    dataPlot(k, 5) = v(169);

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
  dataPlot.conservativeResize(k, outputSize);
  siconos::algebra::io::write(basename + ".dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  double error = 0.0, eps = 1e-11;
  if ((error = siconos::algebra::io::compareRefFile(dataPlot, reference_file_name, eps)) >=
      eps)
    return 1;

  return 0;
}
