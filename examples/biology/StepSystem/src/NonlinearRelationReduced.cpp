#ifndef NONLINEARRELATIONREDUCED_CPP
#define NONLINEARRELATIONREDUCED_CPP

#include "NonlinearRelationReduced.hpp"

//#include "const.h"
#define SICONOS_DEBUG

NonlinearRelationReduced::NonlinearRelationReduced():
  FirstOrderType2R()
{
}

/*y = h(X)*/
void NonlinearRelationReduced::computeh(double t, const BlockVector& x, const SiconosVector& lambda, SiconosVector& y)
{

  y(0) =  8.0 - x(1);

}



/** default function to compute jacobianH
 *  \param double : current time
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */

/** default function to compute jacobianG according to lambda
 *  \param double : current time
 *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
 */


/*g=g(lambda)*/
void NonlinearRelationReduced::computeg(double t, const SiconosVector& lambda, BlockVector& r)
{
  r.setValue(1, 40.0 * (1 - lambda(0)));
  r.setValue(0, 0.0);
}

void NonlinearRelationReduced::computeJachlambda(double t, const BlockVector& x, const SiconosVector& lambda, SimpleMatrix& D)
{
  D.zero();
}

void NonlinearRelationReduced::computeJachx(double t, const BlockVector& x, const SiconosVector& lambda, SimpleMatrix& C)
{

  C(0, 0) =  0;
  C(0, 1) = -1;

}

void NonlinearRelationReduced::computeJacglambda(double t, const SiconosVector& lambda, SimpleMatrix& B)
{
  B(0, 0) =  0.0;
  B(1, 0) = -40.0;
}
#endif

