#ifndef ADJOINTINPUT_CPP
#define ADJOINTINPUT_CPP

#include "adjointInput.hpp"

#include "BlockVector.hpp"
// #define SICONOS_DEBUG

user_defined::adjointInput::adjointInput() : FirstOrderNonLinearR{} {
  setComputehFunction([this](const siconos::algebra::BlockVector& x, double t,
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
                             Eigen::Ref<siconos::algebra::MapVectorType> result) {
    auto betatmp = beta(t, x);

    double betap = 2.0 * ((betatmp(0)) * x(2) + (betatmp(1)) * x(3));

    result(0) = lamb(1) + betap;
    result(1) = 2.0 - lamb(0);
  });

  setComputegFunction([this](const siconos::algebra::BlockVector& x, double t,
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
                             siconos::algebra::BlockVector& result) {
    auto K2P = std::make_shared<siconos::algebra::SiconosVector>(2);
    auto P = std::make_shared<siconos::algebra::SiconosVector>(2);
    (*P)(0) = x(2);
    (*P)(1) = x(3);

    *K2P = *K2 * *P;

    auto betatmp = beta(t, x);

    result(0) = betatmp(0) * (lamb(0) - 1.0);  // R=g_barre(x,lambda_barre)
    result(1) = (betatmp(1)) * (lamb(0) - 1.0);
    result(2) = ((*K2P)(0)) * (lamb(0) - 1.0);
    result(3) = ((*K2P)(1)) * (lamb(0) - 1.0);
  });

  setComputeJacobianhOver_stateFunction(
      [this](const siconos::algebra::BlockVector& x, double t,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
             Eigen::Ref<siconos::algebra::MapType> result) {
        result(0, 0) = x(3);
        result(0, 1) = -x(2);
        result(0, 2) = (-x(1) + 1.0);
        result(0, 3) = x(0);
        result(1, 0) = 0.0;
        result(1, 1) = 0.0;
        result(1, 2) = 0.0;
        result(1, 3) = 0.0;
      });

  setComputeJacobiangOver_stateFunction(
      [this](const siconos::algebra::BlockVector& x, double t,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
             Eigen::Ref<siconos::algebra::MapType> result) {
        auto jacbetaXtmp = JacobianXbeta(t, x);

        result.setZero();
        result(0, 0) = jacbetaXtmp(0, 0) * (lamb(0) - 1.0);
        result(0, 1) = jacbetaXtmp(0, 1) * (lamb(0) - 1.0);
        result(1, 0) = jacbetaXtmp(1, 0) * (lamb(0) - 1.0);
        result(1, 1) = jacbetaXtmp(1, 1) * (lamb(0) - 1.0);
        result(2, 2) = (*K2)(0, 0) * (lamb(0) - 1.0);
        result(2, 3) = (*K2)(0, 1) * (lamb(0) - 1.0);
        result(3, 2) = (*K2)(1, 0) * (lamb(0) - 1.0);
        result(3, 3) = (*K2)(1, 1) * (lamb(0) - 1.0);
      });

  setComputeJacobiangOver_lambdaFunction(
      [this](const siconos::algebra::BlockVector& x, double time,
             const Eigen::Ref<const siconos::algebra::SiconosVector>& lamb,
             Eigen::Ref<siconos::algebra::MapType> result) {
        auto K2P = std::make_shared<siconos::algebra::SiconosVector>(2);
        auto P = std::make_shared<siconos::algebra::SiconosVector>(2);
        (*P)(0) = x(2);
        (*P)(1) = x(3);

        *K2P = *K2 * *P;

        auto betatmp = beta(time, x);

        result(0, 0) = betatmp(0);
        result(0, 1) = 0.0;
        result(1, 0) = betatmp(1);
        result(1, 1) = 0.0;
        result(2, 0) = (*K2P)(0);
        result(2, 1) = 0.0;
        result(3, 0) = (*K2P)(1);
        result(3, 1) = 0.0;
      });
}

void user_defined::adjointInput::initialize(siconos::modeling::Interaction& inter) {
  FirstOrderNonLinearR::initialize(inter);
  K2 = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  K2->setValue(0, 0, 0.0);
  K2->setValue(0, 1, -1.0 / 2.0);
  K2->setValue(1, 0, 1.0 / 2.0);
  K2->setValue(1, 1, 0.0);
}

double user_defined::adjointInput::source(double t) {
  double daux = 0;
  return daux;
}

siconos::algebra::SiconosVector user_defined::adjointInput::beta(
    double t, const siconos::algebra::BlockVector& xvalue) {
  siconos::algebra::SiconosVector res{2};
  res(0) = -1.0 / 2.0 * xvalue(1) + 1.0 / 2.0;
  res(1) = 1.0 / 2.0 * xvalue(0);
#ifdef SICONOS_DEBUG
  std::cout << "beta\n" << std::endl;
  ;
  siconos::algebra::print(*beta);
#endif
  return res;  // RVO
}

siconos::algebra::SiconosMatrix user_defined::adjointInput::JacobianXbeta(
    double t, const siconos::algebra::BlockVector& xvalue) {
  siconos::algebra::SiconosMatrix res{2, 2};

  res.setValue(0, 0, 0.0);
  res.setValue(0, 1, -1.0 / 2.0);
  res.setValue(1, 0, 1.0 / 2.0);
  res.setValue(1, 1, 0.0);
#ifdef SICONOS_DEBUG
  std::cout << "JacXbeta\n" << std::endl;
  ;
  siconos::algebra::print(*JacXbeta);
#endif
  return res;  // RVO
}

#endif
