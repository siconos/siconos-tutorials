/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*!\file BouncingBallTS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include <ReferenceClasses.hpp>
#include <SiconosKernel.hpp>
#include <chrono>

#include "testClass.hpp"

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

auto example() {
  siconos::algebra::SiconosVector q0{7}, q01{7}, velocity0{6};

  siconos::algebra::SiconosMatrix inertia{3, 3};
  // siconos::algebra::SiconosMatrix mass{6, 6};
  auto mass = 10.;
  q0 << 1, 2, 3, 0., 1, 0, 0;
  q01 << 1, 2, 3, 1, 0, 0, 0.;
  velocity0 << 4., 5, 6, 7, 8, 9;
  inertia(0, 0) = 1;
  inertia(1, 1) = 2;
  inertia(2, 2) = 3;

  return std::make_shared<siconos::modeling::NewtonEulerDS>(q0, velocity0, mass, inertia);

  // auto mass_func = [](Eigen::Ref<siconos::algebra::MapVectorType> pos, double time,
  //                     Eigen::Ref<siconos::algebra::MapType> result) {
  //   result.setZero();
  //   result(0, 0) = 1;
  //   result(1, 1) = 2.;
  //   result(2, 2) = 3.;
  // };
}

// void testfunc(const Eigen::Ref<Vector> vin) {
//   auto temp = vin;
//   temp(0) += 3;
//   vin(2) += 2;
//   std::cout << " aaaaaaaa " << temp << " \n";
// }

void computeMgyr(const Eigen::Ref<siconos::algebra::SiconosVector>& twist,
                 const Eigen::Ref<siconos::algebra::SiconosMatrix>& inertiaMatrix,
                 Eigen::Ref<siconos::algebra::SiconosVector> result) {
  auto omega = twist.tail<3>();
  auto inertia = inertiaMatrix.block<3, 3>(3, 3);
  result = omega.cross(inertia * omega);
}

siconos::algebra::SiconosVector setup(int size) {
  Vector vec0{size};

  for (auto& v : vec0) v = 4;
  return vec0;  // RVO
}

void test_var() {
  std::cout << "Start test_var ... \n";
  auto vec0 = setup(3);
  ClassA myclass(vec0, vec0);

  auto var2 = myclass.var_read();  // read-only - var2 is a Map
  for (auto& v : var2) {
    assert(v == 4);
  }

  var2(2) == 14;
  assert((*myclass.var())(2) == 4);

  // pointer access - var is a shared_ptr
  auto var = myclass.var();

  for (auto& v : *var) {
    assert(v == 4);
  }

  (*var)(2) = 128;
  assert((*myclass.var())(2) == 128);

  std::cout << "End test_var ... \n";
}

void test_vector1() {
  std::cout << "Start test_vector1 ... \n";
  auto vec0 = setup(3);
  ClassA myclass(vec0, vec0);

  auto var2 = myclass.vector1();  // Map, read-only
  for (auto& v : var2) {
    assert(v == 4);
  }

  std::cout << myclass.vector1() << "\n";
  ;

  var2(2) == 14;  // Forbidden
  std::cout << myclass.vector1() << "\n";

  std::cout << "End test_vector1 ... \n";
}

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int ndof = 3;  // 100000000;  // degrees of freedom for the ball
    double m = 1;           // Ball mass
    double g = 9.81;        // Gravity

    // -- Initial positions and velocities --
    //          Vector q0{ndof};
    // auto q0 = std::make_shared<Vector>(ndof);
    // q0->setZero();
    // q0(0)(0) = 1.2;

    Vector q0{ndof};
    q0.setZero();
    q0(0) = 1.;

    Vector v0{ndof};
    v0.setZero();
    v0(0) = 1.;
    Vector z0{ndof};
    z0.setZero();
    z0(0) = 1.;

    z0 = v0 + q0;

    // -- Build
    // auto ball = std::make_shared<siconos::internal::devel_model::ClassA>(q0);

    // Vector fext{ndof};
    // fext.setZero();
    // fext(1) = 112.2;
    // ball->setConstantVector2(fext);

    // //siconos::algebra::print(*ball->vector2());
    // ball->computeVector2(0.4);
    // //siconos::algebra::print(*ball->vector2());

    // //siconos::algebra::print(*ball);

    std::cout << "Second test ... \n";
    // -- Set external forces (weight) --
    auto myforces = [m, g](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
      int i = 0;
      for (auto& v : result) v = time * i++;

      //  siconos::tools::print("call plugin", result);
    };

    // auto ball2 = std::make_shared<siconos::internal::devel_model::ClassA>(q0);

    // std::cout << (*ball2->vector1())(0) << " " << (*ball2->vector3())(0) << "\n";

    // ball2->computeVector2(0.);
    // if (ball2->hasVector2()) {
    //   siconos::algebra::print(*ball2->vector2());
    // }

    // //    siconos::algebra::print(*ball2);
    // ball2->setComputeVector2Function(myforces);
    // // //siconos::algebra::print(*ball2);
    // // //siconos::algebra::print(*ball2->vector2());
    // ball2->computeVector2(0.);
    // // //siconos::algebra::print(*ball2->vector2());
    // ball2->computeVector2(1.);
    // siconos::algebra::print(*ball2->vector2());
    // if (ball2->hasVector2()) {
    //   siconos::algebra::print(*ball2->vector2());

    //   auto res = 3 * ball2->vector2_view();
    //   std::cout << res << "\n";
    // }
    // siconos::algebra::print(*ball2->vector2());

    // /// ----- with vectorDirect -----

    // ball2->computeVectorDirect(0.);
    // // if (ball2->hasVectorDirect()) {
    // //   siconos::algebra::print(*ball2->vectorDirect());
    // // }

    // //    siconos::algebra::print(*ball2);
    // ball2->setComputeVectorDirectFunction(myforces);
    // // //siconos::algebra::print(*ball2);
    // // //siconos::algebra::print(*ball2->vectorDirect());
    // ball2->computeVectorDirect(0.);
    // // //siconos::algebra::print(*ball2->vectorDirect());
    // ball2->computeVectorDirect(1.);
    // // siconos::algebra::print(*ball2->vectorDirect());
    // if (ball2->hasVectorDirect()) {
    //   //   //   siconos::algebra::print(*ball2->vectorDirect());

    //   auto res = 3 * ball2->vectorDirect_view();
    //   std::cout << res << "\n";
  }
  // // siconos::algebra::print(*ball2->vectorDirect());
  ////  --------------------------------

  // /// ----- with vectorSpan -----
  // // -- Set external forces (weight) --
  // auto myforces_span = [m, g](double time, std::span<double> result) {
  //   int i = 0;
  //   for (auto& v : result) v = time * i++;

  //   result[1] = 12;
  //   //  siconos::tools::print("call plugin", result);
  // };

  // ball2->computeVectorSpan(0.);
  // if (ball2->hasVectorSpan()) {
  //   siconos::algebra::print(*ball2->vectorSpan());
  // }

  // //    siconos::algebra::print(*ball2);
  // ball2->setComputeVectorSpanFunction(myforces_span);
  // // //siconos::algebra::print(*ball2);
  // // //siconos::algebra::print(*ball2->vectorSpan());
  // ball2->computeVectorSpan(0.);
  // // //siconos::algebra::print(*ball2->vectorSpan());
  // ball2->computeVectorSpan(1.);
  // siconos::algebra::print(*ball2->vectorSpan());
  // if (ball2->hasVectorSpan()) {
  //   siconos::algebra::print(*ball2->vectorSpan());

  //   auto res = 3 * ball2->vectorSpan_view();
  //   std::cout << res << "\n";
  // }
  // siconos::algebra::print(*ball2->vectorSpan());
  // ////  --------------------------------

  // siconos::algebra::SiconosMatrix mass{ndof, ndof};
  // mass.setZero();
  // mass(1, 2) = -m * g;
  // mass(0, 1) = 12;

  // // mass(3,4) = 12; // ça marche avec ndof = 3, pourquoi ????

  // ball2->setConstantMatrix1(mass);

  // siconos::algebra::print(*ball2->matrix1());

  // auto mass_func = [m, g](Eigen::Ref<siconos::algebra::MapVectorType> pos, double time,
  //                         Eigen::Ref<siconos::algebra::MapType> result) {
  //   int i = 0;
  //   // for (auto& v : result) v = time * i++;

  //   //      result << 1, 2, 3, 4;

  //   result(2, 2) = -m * g;

  //   //  siconos::tools::print("call plugin", result);
  // };

  // ball2->setComputeMatrix1Function(mass_func);
  // siconos::algebra::print(*ball2->matrix1());

  // Vector pos{ndof};
  // pos.setZero();
  // pos(1) = 8;
  // ball2->computeMatrix1(1., pos);
  // siconos::algebra::print(*ball2->matrix1());

  // example();

  // // int size = 700000000;
  // int size = 3;
  // Vector vec0{size};

  for (auto& v : vec0) v = 4;
  // Vector vec01{size};

  // for (auto& v : vec01) v = 4;

  // auto vec1 = std::make_shared<Vector>(size);
  // for (auto& v : *vec1) v = 4;

  // Cas const& pour param1
  // OK :
  // ClassA myclass{Eigen::Ref<Vector>(vec01), vec0};
  // ClassA myclass{*vec1, vec0};
  // ClassA myclass{Eigen::Ref<Vector>(*vec1), vec0};
  // NON

  // Cas const& pour param1, & pour param2
  // OK :

  // NON
  // ClassA myclass{vec01, vec0};
  // ClassA myclass{vec01, Eigen::Ref<Vector>{vec0}};
  // ClassA myclass{vec01, *vec1};
  // ClassA myclass{vec01, Eigen::Ref<Vector>(*vec1)};

  // OK
  // Eigen::Ref<Vector> param2_ref(vec0);
  // ClassA myclass{vec01, param2_ref};
  // Eigen::Ref<Vector> param2_ref(*vec1);
  // ClassA myclass{vec01, param2_ref};

  // Cas std (ni const ni &, juste Eigen::Ref) et avec un seul param pour constr
  // Ok :
  // //    ClassA myclass{vec0, vec0};
  // ClassA myclass{Eigen::Ref<Vector>(vec0)};
  // ClassA myclass{*vec1};
  // Non :

  // std::cout << "And the winner is ...\n"
  //           << (*myclass.vector1())(1) << " " << vec0(1) << " " << (*vec1)(1) << " "
  //           << (*myclass.constvector1())(1) << " " << vec01(1) << "\n";

  // vec0(1) = 14;
  // (*vec1)(1) = 14;
  // vec01(1) = 38;
  // std::cout << "And the winner is ...\n"
  //           << (*myclass.vector1())(1) << " " << vec0(1) << " " << (*vec1)(1) << " "
  //           << (*myclass.constvector1())(1) << " " << vec01(1) << "\n";

  // myclass.update(127, 39);
  // std::cout << "And the winner is ...\n"
  //           << (*myclass.vector1())(1) << " " << vec0(1) << " " << (*vec1)(1) << " "
  //           << (*myclass.constvector1())(1) << " " << vec01(1) << "\n";

  ClassA myclassa{vec0, vec0};
  // auto T = myclassa.T_view();

  // // // auto val = T(1) + T(2);

  // std::cout << " ################# \n";

  // std::cout << T(12) << " " << T(23) << "\n";

  // auto T2 = T.dot(T);

  // myclassa.reset();

  // std::cout << T(1) << " " << vec0(1) << " " << myclassa.vector1()(1) << "\n";

  // myclassa.update(48);

  // std::cout << T(1) << " " << vec0(1) << " " << myclassa.vector1()(1) << "\n";
  // vec0(1) = 123;
  // std::cout << T(1) << " " << vec0(1) << " " << myclassa.vector1()(1) << "\n";

  // testfunc(v0);

  siconos::algebra::SiconosMatrix toto{3, 6};
  for (auto i = 0; i < toto.rows(); i++)
    for (auto j = 0; j < toto.rows(); j++) toto(i, j) = i + j;

  siconos::algebra::SiconosVector q{6};
  q << 1, 2, 3, 4, 5, 6;

  Eigen::Map<const siconos::algebra::SiconosVector3> qvect(q.data() +
                                                           4);  // view onto  q4, 5, 6

  std::cout << toto << "\n";

  for (unsigned int j = 1; j < 2; j++) {
    Eigen::Map<Eigen::Vector3d> mcol(toto.col(j).data());
    // auto t = 2 * qvect.cross(mcol);
    mcol += qvect.cross(2 * qvect.cross(mcol));

    // m.col(j) = mcol;
  }
  std::cout << toto << "\n";

  using SiconosDiagonalMatrix = Eigen::DiagonalMatrix<double_t, Eigen::Dynamic>;
  using DiagonalMatrixMapType = Eigen::Map<SiconosDiagonalMatrix>;
  using ConstDiagonalMatrixMapType = Eigen::Map<const SiconosDiagonalMatrix>;

  SiconosDiagonalMatrix diagmat{q};

  std::cout << " STATATATA " << diagmat.diagonal() << "\n";

  std::shared_ptr<siconos::algebra::MapVectorType> stiffnessMatrix_view;
  stiffnessMatrix_view = std::make_shared<siconos::algebra::MapVectorType>(
      diagmat.diagonal().data(), diagmat.diagonal().size());

  std::cout << *stiffnessMatrix_view << "\n";

  diagmat.diagonal()(2) = 112;
  std::cout << diagmat.diagonal() << "\n";
  std::cout << *stiffnessMatrix_view << "\n";

  (*stiffnessMatrix_view)(3) = 128;
  std::cout << diagmat.diagonal() << "\n";
  std::cout << *stiffnessMatrix_view << "\n";
  // auto truc = 2 * T;
  // std::cout << truc << "\n";

  std::cout << q << "\n";

  Eigen::VectorXd x(3);
  Eigen::VectorXd y(3);

  x << 1.0, 2.0, 3.0;
  y << 4.0, 5.0, 6.0;

  // Produit élément par élément
  auto w1 = x.array() * y.array();

  auto w = stiffnessMatrix_view->array() * q.array();
  // Afficher le résultat
  std::cout << "Produit élément par élément (w) :\n" << w << std::endl;

  Matrix inertiaMatrix{6, 6};
  for (auto i = 0; i < inertiaMatrix.rows(); i++)
    for (auto j = 0; j < inertiaMatrix.rows(); j++) inertiaMatrix(i, j) = i + j;

  Vector result{3};
  result.setZero();
  computeMgyr(q, inertiaMatrix, result);

  std::cout << "I" << inertiaMatrix << "\n";
  std::cout << "q" << q << "\n";
  std::cout << "res" << result << "\n";

  // Définir deux vecteurs
  Vector v1(3);
  Vector v2(2);

  v1 << 1, 2, 3;  // Initialisation de v1
  v2 << 4, 5;     // Initialisation de v2

  // Créer une matrice A
  Matrix A(5, 5);
  A << 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      1;  // Matrice identité 5x5

  // // Créer des maps pour éviter des copies
  // Eigen::Map<Vector> xx(v1.data(), v1.size() + v2.size()); // La map représente le
  // vecteur concaténé

  // // Combiner les deux vecteurs en utilisant Map
  // xx.head(v1.size()) = v1;        // Mettre v1 dans les premières entrées
  // xx.tail(v2.size()) = v2;        // Mettre v2 dans les dernières entrées

  // // Résultat dans y
  // Vector yy = A * xx; // Multiplication

  // std::cout << "y:\n" << yy << std::endl;

  test_var();
  test_vector1();

  return 0;
}

catch (...) {
  siconos::exception::process();
  return 1;
}
}
