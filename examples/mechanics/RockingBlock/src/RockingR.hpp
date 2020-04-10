#ifndef ROCKING_R_HPP
#define ROCKING_R_HPP

#include "LagrangianScleronomousR.hpp"

class RockingBlockR : public LagrangianScleronomousR
{
public:

  double LengthBlock = 0.2;
  double HeightBlock = 0.1;

  void computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y)
  {
    y(0) = q(1) - 0.5 * LengthBlock * sin(q(2)) - 0.5 * HeightBlock * cos(q(2));
  }

  void computeJachq(SiconosVector& q, SiconosVector& z)
  {
    (*_jachq)(0,0) = 0.0;
    (*_jachq)(0,1) = 1.0;
    (*_jachq)(0,2) = -0.5 * LengthBlock * cos(q(2)) + 0.5 * HeightBlock * sin(q(2));
  }

  void computeDotJachq(SiconosVector& q, SiconosVector& z, SiconosVector& qdot)

  {
    // TMP : this stage won't be necessary when the LagrangianScleronomousR
    // class will be updated
    if(!_dotjachq)
    {
      unsigned int sizeY = 1;
      unsigned int sizeDS = 3;
      _dotjachq.reset(new SimpleMatrix(sizeY, sizeDS));
    }
    (*_dotjachq)(0,0) = 0.0;
    (*_dotjachq)(0,1) = 0.0;
    (*_dotjachq)(0,2) = (0.5 * LengthBlock * sin(q(2)) + 0.5 * HeightBlock * cos(q(2))) * qdot(2);
  }

};

class RockingBlockR2 : public LagrangianScleronomousR
{
public:

  double LengthBlock = 0.2;
  double HeightBlock = 0.1;

  void computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y)
  {
    y(0) = q(1) + 0.5 * LengthBlock * sin(q(2)) - 0.5 * HeightBlock * cos(q(2));
  }

  void computeJachq(SiconosVector& q, SiconosVector& z)
  {
  (*_jachq)(0,0) = 0.0;
  (*_jachq)(0,1) = 1.0;
  (*_jachq)(0,2) = 0.5 * LengthBlock * cos(q(2)) + 0.5 * HeightBlock * sin(q(2));
  }

  void computeDotJachq(SiconosVector& q, SiconosVector& z, SiconosVector& qdot)
  {
    // TMP : this stage won't be necessary when the LagrangianScleronomousR
    // class will be updated
    if(!_dotjachq)
    {
      unsigned int sizeY = 1;
      unsigned int sizeDS = 3;
      _dotjachq.reset(new SimpleMatrix(sizeY, sizeDS));
    }
    (*_dotjachq)(0,0) = 0.0;
    (*_dotjachq)(0,1) = 0.0;
    (*_dotjachq)(0,2) = (-0.5 * LengthBlock * sin(q(2)) + 0.5 * HeightBlock * cos(q(2))) * qdot(2);

  }

};



#endif
