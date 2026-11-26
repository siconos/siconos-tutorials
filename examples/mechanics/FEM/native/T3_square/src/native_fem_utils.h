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
#ifndef NFEM_UTILS
#define NFEM_UTILS

#include <FENode.hpp>
#include <FiniteElementLinearTIDS.hpp>
#include <SiconosKernel.hpp>

namespace native_fem_examples {

// Run simulation, work for all examples with T3 FEM
int run_T3_simulation(
    std::shared_ptr<siconos::simulation::TimeStepping> simulation,
    std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS> FEsolid,
    std::string outputfile_name, std::string reference_file_name);

}  // namespace native_fem_examples
#endif
