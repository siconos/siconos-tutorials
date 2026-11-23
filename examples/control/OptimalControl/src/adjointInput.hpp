#ifndef ADJOINTINPUT_H
#define ADJOINTINPUT_H

#include <FirstOrderNonLinearR.hpp>
#include <memory>

namespace user_defined {

class adjointInput : public siconos::modeling::FirstOrderNonLinearR {
 protected:
  std::shared_ptr<siconos::algebra::SiconosMatrix> K2{nullptr};

 public:
  adjointInput();
  virtual ~adjointInput() noexcept = default;

  virtual void initialize(siconos::modeling::Interaction& inter) override;

  double source(double t);

  siconos::algebra::SiconosVector beta(double t, const siconos::algebra::BlockVector& xvalue);

  siconos::algebra::SiconosMatrix JacobianXbeta(double t,
                                                const siconos::algebra::BlockVector& xvalue);
};
}  // namespace user_defined

#endif
