#ifndef BALL_RELATION_HPP
#define BALL_RELATION_HPP

class BallR : public LagrangianScleronomousR
{
  double alpha = 0.1;
  
public:


  void computeh(const BlockVector& q, BlockVector& z, SiconosVector& y)
  {
    y.setValue(0, q.getValue(0) + alpha * q.getValue(1));
  }

  void computeJachq(const BlockVector& q, BlockVector& z)
  {
    _jachq->setValue(0, 0, 1.0);
    _jachq->setValue(0, 1, alpha);
  }
};


#endif
