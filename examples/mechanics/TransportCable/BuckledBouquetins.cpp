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

#include <TransportCableManager.h>
#include <TransportCableModel.h>
#include <TransportCableProfile.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <fstream>
#include <string>

// Siconos names used unqualified in this example
using siconos::fem::cable::TransportCableManager;
using siconos::fem::cable::TransportCableModel;
using siconos::fem::cable::TransportCableProfile;
using siconos::fem::cable::TransportCableResult;

int main()

{
  try {
    // Reads model parameters from a json file
    // Mechanical params, geometry, supports positions ...
    std::string modelFile = "buckled_bouquetins.json";
    std::ifstream modelStream(modelFile);
    auto model = std::make_shared<TransportCableModel>(nlohmann::json::parse(modelStream));

    // Creates the result object: it handles two ropeways (sequence of ropes),
    // the supports (pylons and pulleys), and some extra variables related to
    // contacts.
    auto results = std::make_shared<TransportCableResult>();

    // Creates the profile, object which owns the model and the result and which
    // is responsible for the computation of the initial state/profile using
    // catenary and FEM.
    auto profil = std::make_shared<TransportCableProfile>(*model, *results);

    // -- Applies catenary equations to compute a first profile of the ropeways
    // --
    int nb_nodes = 50;   // Catenary, number of nodes per rope span
    double tol = 1e-10;  // Tol. used in Newton-Raphson for catenary equation
    int nmax = 20;       // Newton-Raphson, max number of iterations
    profil->computeInitialProfile(nb_nodes, tol, nmax);

    // Save ropeways variables into json file
    nlohmann::ordered_json out;
    results->to_json(out, "ropeway");
    // results->to_json(out);

    std::ofstream out_ropes("catenary.json");
    out_ropes << std::setw(4) << out << std::endl;

    // std::ifstream in("bouquetins_ref.json");
    // json reader;
    // in >> reader;
    // auto qref1 = ioVector::readVectorFromJson(reader["ropes_up"]["q"]);
    // auto qref2 = ioVector::readVectorFromJson(reader["ropes_down"]["q"]);

    // std::ifstream in2("ropes.json");
    // auto reader2 = nlohmann::json::parse()
    //  std::ifstream in("ropes.json");
    // auto q1 = ioVector::readVectorFromJson(out["ropes_up"]["q"]);
    // auto q2 = ioVector::readVectorFromJson(out["ropes_down"]["q"]);

    // std::cout << qref1->norm() << " " << qref1->size() << std::endl;

    // std::cout << ((*qref1) == (*q1)) << std::endl;

    // -- Fem part --

    int numberOfElements = 1400;  // FEM number of nodes
    double eps = 0.1;
    tol = 1e-3;  // tolerance used to activate constraints
    profil->initializeFEM(numberOfElements, eps, tol);

    // nlohmann::ordered_json out;
    results->to_json(out);
    std::ofstream out_ropes_fem("ropes_fem.json");
    out_ropes_fem << std::setw(4) << out << std::endl;

    auto positions = siconos::algebra::io::readVectorFromJson(out["q"]);
    std::cout << positions.norm() << " " << positions.size() << std::endl;

    std::cout << "PAR 2 \n";
    auto manager = std::make_shared<TransportCableManager>(modelFile);
    std::string outFile = "results.json";
    nlohmann::ordered_json out2;
    nlohmann::json args;
    manager->computeFEM(args, outFile, out2);

    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
