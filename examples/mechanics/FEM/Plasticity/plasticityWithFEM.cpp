
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
//-----------------------------------------------------------------------
//
//  CircuitRLCD  : sample of an electrical circuit involving :
//  - a linear dynamical system consisting of an LC oscillator (1 microF , 10 mH)
//  - a non smooth system (a 1000 Ohm resistor in series with a diode) in parallel
//    with the oscillator
//
//  Expected behavior :
//  The initial state of the oscillator provides an initial energy.
//  The period is 2 Pi sqrt(LC) ~ 0,628 ms.
//  A positive voltage across the capacitor allows current to flow
//  through the resistor-diode branch , resulting in an energy loss :
//  the oscillation damps.
//
//  State variables :
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  Since there is only one dynamical system, the interaction is defined by :
//  - a complementarity law between diode current and voltage where y stands
//    for the reverse voltage across the diode and lambda stands for the
//    the diode current
//  - a linear time invariant relation between the state variables and
//    y and lambda (derived from Kirchhoff laws)
//
//-----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include <chrono>

#include "FiniteElementLinearTIDS.hpp"
#include "MeshUtils.hpp"
// #include "SiconosAlgebraProd.hpp"
#include "SiconosVector.hpp"

using namespace std;

int main(int argc, char* argv[]) {
  double t0 = 0.0;
  double T = 500e-1;       // Total simulation time
  double h_step = 1.0e-2;  // Time step

  double sigmaMax = 100.0;

  // string gmsh_filename = "./mesh_data/triangle_felippa.msh";
  //    string gmsh_filename = "./mesh_data/triangle_reference.msh";
  //    string gmsh_filename = "./mesh_data/triangle.msh";
  string gmsh_filename = "./mesh_data/square_6.msh";
  //      string gmsh_filename = "./mesh_data/square_200.msh";
  // string gmsh_filename = "./mesh_data/square_2720.msh";

  std::shared_ptr<siconos::mechanics::fem::native::Mesh> mesh(
      siconos::mechanics::fem::native::createMeshFromGMSH2(gmsh_filename));
  // mesh->display(false);

  writeMeshforPython(mesh);

  int bulk_material_tag = 1;
  int boundary_condition_tag = 2;
  int applied_force_tag = 3;

  double density = 7000;
  std::shared_ptr<Material> mat1 = std::make_shared<Material>(density, 1.0e3, 0.3);
  std::map<unsigned int, std::shared_ptr<Material> > materials = {{bulk_material_tag, mat1}};

  std::shared_ptr<siconos::mechanics::fem::native::FiniteElementLinearTIDS> FEsolid =
      std::make_shared<siconos::mechanics::fem::native::FiniteElementLinearTIDS>(
          mesh, materials, siconos::algebra::UblasType::SPARSE);
  std::shared_ptr<siconos::mechanics::fem::native::FiniteElementModel> femodel =
      FEsolid->FEModel();

  string Modeltitle = "2dTriangle";
  try {
    int nElements = (femodel->elements()).size();
    double ndofs = (mesh->vertices()).size() * mesh->dim();
    unsigned dim = 2 * ndofs + 3 * nElements;
    // --- Dynamical system specification ---
    auto init_state = std::make_shared<siconos::algebra::SiconosVector>(dim);
    auto Fext = std::make_shared<siconos::algebra::SiconosVector>(dim);
    Fext->setValue(3, 100.0);
    Fext->setValue(9, 100.0);
    auto mass = std::make_shared<siconos::algebra::SimpleMatrix>(ndofs, ndofs);
    auto stiffness = std::make_shared<siconos::algebra::SimpleMatrix>(ndofs, ndofs);
    auto Bfem = std::make_shared<siconos::algebra::SimpleMatrix>(3, ndofs);
    femodel->computeMassMatrix(mass, materials);
    std::shared_ptr<siconos::algebra::SimpleMatrix> D =
        std::make_shared<siconos::algebra::SimpleMatrix>(3, 3);
    double E = mat1->elasticYoungModulus();
    double nu = mat1->poissonCoefficient();

    double coef = E / ((1 + nu) * (1 - 2. * nu));
    (*D)(0, 0) = coef * (1. - nu);
    (*D)(0, 1) = coef * nu;
    (*D)(0, 2) = 0.0;

    (*D)(1, 0) = (*D)(0, 1);
    (*D)(1, 1) = (*D)(0, 0);
    (*D)(1, 2) = 0.0;

    (*D)(2, 0) = 0.0;
    (*D)(2, 1) = 0.0;
    (*D)(2, 2) = 0.5 * coef * (1.0 - 2 * nu);

    femodel->computeStiffnessMatrix(stiffness, materials);

    std::cout << "nElements:" << nElements << std::endl;
    auto BigB = std::make_shared<siconos::algebra::SimpleMatrix>(3 * nElements, ndofs);
    std::shared_ptr<siconos::algebra::SimpleMatrix> DBigB =
        std::make_shared<siconos::algebra::SimpleMatrix>(3 * nElements, ndofs);
    std::shared_ptr<siconos::algebra::SimpleMatrix> DBfem =
        std::make_shared<siconos::algebra::SimpleMatrix>(3, ndofs);
    int elem_cnt = 0;
    for (std::shared_ptr<siconos::mechanics::fem::native::FiniteElement> fe :
         femodel->elements()) {
      femodel->computeElementaryBMatrix_direct(*fe, *Bfem);
      prod(*D, *Bfem, *DBfem, true);
      femodel->AssembleElementary_B_Matrix(BigB, *Bfem, *fe, elem_cnt);
      femodel->AssembleElementary_B_Matrix(DBigB, *DBfem, *fe, elem_cnt);
      elem_cnt++;
    }
    auto M = std::make_shared<siconos::algebra::SimpleMatrix>(dim, dim);
    for (int i = 0; i < ndofs; i++) {
      for (int j = 0; j < ndofs; j++) {
        M->setValue(i, j, mass->getValue(i, j));
      }
    }
    for (int i = 0; i < ndofs; i++) {
      M->setValue(i + ndofs, i + ndofs, 1.0);
    }

    for (int i = 0; i < 3 * nElements; i++) {
      M->setValue(i + 2 * ndofs, i + 2 * ndofs, 1.0);
    }

    std::shared_ptr<siconos::algebra::SimpleMatrix> BT =
        std::make_shared<siconos::algebra::SimpleMatrix>(ndofs, 3 * nElements);
    BT->trans(*BigB);
    std::shared_ptr<siconos::algebra::SimpleMatrix> BTDB =
        std::make_shared<siconos::algebra::SimpleMatrix>(ndofs, ndofs);

    prod(*BT, *DBigB, *BTDB, true);

    auto A = std::make_shared<siconos::algebra::SimpleMatrix>(dim, dim);
    for (int i = 0; i < ndofs; i++) {
      for (int j = 0; j < 3 * nElements; j++) {
        A->setValue(i, j + 2 * ndofs, -BigB->getValue(j, i));
      }
    }
    for (int i = 0; i < 3 * nElements; i++) {
      for (int j = 0; j < ndofs; j++) {
        A->setValue(i + 2 * ndofs, j, DBigB->getValue(i, j));
      }
    }
    for (int i = 0; i < ndofs; i++) {
      A->setValue(i + ndofs, i, 1.0);
    }

    auto triangle =
        std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(init_state, A, Fext);
    triangle->setMPtr(M);

    // --- Interaction between linear system and non smooth system ---
    int nbConstraints = 10;
    std::shared_ptr<siconos::algebra::SimpleMatrix> C =
        std::make_shared<siconos::algebra::SimpleMatrix>(nbConstraints, dim);
    // Position Constraints
    C->setValue(0, ndofs, -1.0);
    C->setValue(1, ndofs, 1.0);
    C->setValue(2, ndofs + 1, -1.0);
    C->setValue(3, ndofs + 1, 1.0);
    C->setValue(4, ndofs + 6, -1.0);
    C->setValue(5, ndofs + 6, 1.0);
    C->setValue(6, ndofs + 7, -1.0);
    C->setValue(7, ndofs + 7, 1.0);
    // Stress Constraints
    C->setValue(8, 2 * ndofs, -1.0);
    C->setValue(9, 2 * ndofs, 1.0);

    auto B = std::make_shared<siconos::algebra::SimpleMatrix>(dim, nbConstraints);
    B->setValue(0, 0, -1.0);
    B->setValue(0, 1, 1.0);
    B->setValue(0 + 1, 2, -1.0);
    B->setValue(0 + 1, 3, 1.0);
    B->setValue(0 + 6, 4, -1.0);
    B->setValue(0 + 6, 5, 1.0);
    B->setValue(0 + 7, 6, -1.0);
    B->setValue(0 + 7, 7, 1.0);
    // Stress Constraints
    B->setValue(2 * ndofs, 8, -1.0);
    B->setValue(2 * ndofs, 9, 1.0);
    auto e = std::make_shared<siconos::algebra::SiconosVector>(nbConstraints);
    e->setValue(0, 0.00001);
    e->setValue(1, 0.00001);
    e->setValue(2, 0.00001);
    e->setValue(3, 0.00001);
    e->setValue(4, 0.00001);
    e->setValue(5, 0.00001);
    e->setValue(6, 0.00001);
    e->setValue(7, 0.00001);
    // Stress Constraints
    e->setValue(8, sigmaMax);
    e->setValue(9, sigmaMax);

    auto LTIRspring = std::make_shared<siconos::modeling::FirstOrderLinearTIR>(C, B);
    LTIRspring->setePtr(e);

    auto NSLaw =
        std::make_shared<siconos::modeling::ComplementarityConditionNSL>(nbConstraints);

    auto InterTriangle = std::make_shared<siconos::modeling::Interaction>(NSLaw, LTIRspring);
    // --- Model creation ---
    auto springDS = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    springDS->setTitle(Modeltitle);
    springDS->insertDynamicalSystem(triangle);

    // link the interaction and the dynamical system
    springDS->link(InterTriangle, triangle);

    InterTriangle->computeOutput(t0, 0);
    InterTriangle->computeInput(t0, 0);

    springDS->display();

    // ------------------
    // --- Simulation ---
    // ------------------
    double theta = 0.5000000000001;
    // -- (1) OneStepIntegrators --

    auto OSI = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);
    // -- (2) Time discretisation --
    auto TiDis = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h_step);
    // --- (3) one step non smooth problem
    auto LCP_triangle = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    // -- (4) Simulation setup with (1) (2) (3)
    auto springMassTS = std::make_shared<siconos::simulation::TimeStepping>(springDS, TiDis,
                                                                            OSI, LCP_triangle);
    double h = springMassTS->timeStep();
    int N = ceil((T - t0) / h);  // Number of time steps
    int k = 0;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    siconos::algebra::SimpleMatrix dataPlot(N, 22);

    // For the initial time step:

    // time
    dataPlot(k, 0) = springMassTS->nextTime();
    dataPlot(k, 1) = (*triangle->x())(2 * ndofs);
    dataPlot(k, 2) = (*triangle->x())(2 * ndofs + 1);
    dataPlot(k, 3) = (*triangle->x())(2 * ndofs + 2);
    dataPlot(k, 4) = (*triangle->x())(2);
    dataPlot(k, 5) = (*triangle->x())(3);
    dataPlot(k, 6) = (*triangle->x())(ndofs + 2);
    dataPlot(k, 7) = (*triangle->x())(ndofs + 3);
    dataPlot(k, 8) = (InterTriangle->getLambda(0))(8);
    dataPlot(k, 9) = (InterTriangle->getLambda(0))(9);

    double x, v, sigma, plasticRate, epElastic, ekinetic, plasticDeformation = 0.0,
                                                          plasticDissipation = 0.0;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    // --- Time loop  ---
    std::string filename = "triangle.state";
    FILE* foutput = fopen(filename.c_str(), "w");
    fclose(foutput);
    for (k = 1; k < N; ++k) {
      // solve ...
      springMassTS->computeOneStep();
      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = springMassTS->nextTime();
      dataPlot(k, 1) = (*triangle->x())(2 * ndofs);
      dataPlot(k, 2) = (*triangle->x())(2 * ndofs + 1);
      dataPlot(k, 3) = (*triangle->x())(2 * ndofs + 2);
      dataPlot(k, 4) = (*triangle->x())(2);
      dataPlot(k, 5) = (*triangle->x())(3);
      dataPlot(k, 6) = (*triangle->x())(ndofs + 2);
      dataPlot(k, 7) = (*triangle->x())(ndofs + 3);
      dataPlot(k, 8) = (InterTriangle->getLambda(0))(8);
      dataPlot(k, 9) = (InterTriangle->getLambda(0))(9);

      // transfer of state i+1 into state i and time incrementation
      springMassTS->nextStep();
    }
    // Number of time iterations
    cout << "Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << endl;
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Computation time : " << elapsed << " ms" << endl;

    // dataPlot (ascii) output
    siconos::algebra::io::write("plasticTriangle.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

  }

  // --- Exceptions handling ---
  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
