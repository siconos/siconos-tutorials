#ifndef NONLINEARRELATION_H
#define NONLINEARRELATION_H

#include "SiconosKernel.hpp"

class NonlinearRelation : public FirstOrderType2R
{
protected:
public:
  NonlinearRelation();
  virtual ~NonlinearRelation() {};

  /** default function to compute h */
  virtual void computeh(double t, const BlockVector& x, const SiconosVector& lambda, SiconosVector& y);

  /** default function to compute g */
  virtual void computeg(double t, const SiconosVector& lambda, BlockVector& r);

  /** default function to compute jacobian of h w.r.t x and lambda */
  virtual void computeJachx(double t, const BlockVector& x, const SiconosVector& lambda, SimpleMatrix& C);
  virtual void computeJachlambda(double t, const BlockVector& x, const SiconosVector& lambda, SimpleMatrix& D);

  /** default function to compute jacobian of g  w.r.t  lambda  */
  virtual void computeJacglambda(double t, const SiconosVector& lambda, SimpleMatrix& B);

};

TYPEDEF_SPTR(NonlinearRelation);

#endif
