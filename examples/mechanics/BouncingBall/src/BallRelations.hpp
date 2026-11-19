#ifndef BALL_RELATION_HPP
#define BALL_RELATION_HPP

#include <SiconosKernel.hpp>

namespace user_defined {

class BallR : public siconos::modeling::LagrangianScleronomousR {
  double alpha = 0.1;

 public:
  BallR() : LagrangianScleronomousR() {
    setComputehFunction(
        [this](const siconos::algebra::BlockVector& q,
               Eigen::Ref<siconos::algebra::MapVectorType> y) { y(0) = q(0) + alpha * q(1); });

    setComputeContactMatrixFunction([this](const siconos::algebra::BlockVector& q,
                                           Eigen::Ref<siconos::algebra::MapType> result) {
      result.setValue(0, 0, 1.0);
      result.setValue(0, 1, alpha);
    });
  };
};
}  // namespace user_defined
#endif
