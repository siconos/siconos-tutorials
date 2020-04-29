#ifndef ROCKING_R_HPP
#define ROCKING_R_HPP

#include "LagrangianScleronomousR.hpp"

class RockingBlockR : public LagrangianScleronomousR
{
public:

  double LengthBlock = 0.2;
  double HeightBlock = 0.1;

  void computeh(const BlockVector& q, BlockVector& z, SiconosVector& y)
  {
    double q1 = q.getValue(1);
    double q2 = q.getValue(2);
    y.setValue(0, q1 - 0.5 * LengthBlock * sin(q2) - 0.5 * HeightBlock * cos(q2));
  }

  void computeJachq(const BlockVector& q, BlockVector& z)
  {
    double q2 = q.getValue(2);
    _jachq->setValue(0, 0, 0.0);
    _jachq->setValue(0, 1, 1.0);
    _jachq->setValue(0, 2, -0.5 * LengthBlock * cos(q2) + 0.5 * HeightBlock * sin(q2));
  }

  void computeDotJachq(const BlockVector& q, BlockVector& z, const BlockVector& qdot)

  {
    double q2 = q.getValue(2);
    double qdot2 = qdot.getValue(2);
    _dotjachq->setValue(0, 0, 0.0);
    _dotjachq->setValue(0, 1, 0.0);
    _dotjachq->setValue(0, 2, (0.5 * LengthBlock * sin(q2) + 0.5 * HeightBlock * cos(q2)) * qdot2);
  }

};

class RockingBlockR2 : public LagrangianScleronomousR
{
public:

  double LengthBlock = 0.2;
  double HeightBlock = 0.1;

  void computeh(const BlockVector& q, BlockVector& z, SiconosVector& y)
  {
    double q1 = q.getValue(1);
    double q2 = q.getValue(2);
    y.setValue(0, q1 + 0.5 * LengthBlock * sin(q2) - 0.5 * HeightBlock * cos(q2));
  }

  void computeJachq(const BlockVector& q, BlockVector& z)
  {
    double q2 = q.getValue(2);
    _jachq->setValue(0, 0, 0.0);
    _jachq->setValue(0, 1, 1.0);
    _jachq->setValue(0, 2, 0.5 * LengthBlock * cos(q2) + 0.5 * HeightBlock * sin(q2));
  }
  
  void computeDotJachq(const BlockVector& q, BlockVector& z, const BlockVector& qdot)
  {
    double q2 = q.getValue(2);
    double qdot2 = qdot.getValue(2);
    _dotjachq->setValue(0, 0, 0.0);
    _dotjachq->setValue(0, 1, 0.0);
    _dotjachq->setValue(0, 2, (-0.5 * LengthBlock * sin(q2) + 0.5 * HeightBlock * cos(q2)) * qdot2);
  }

};



#endif
