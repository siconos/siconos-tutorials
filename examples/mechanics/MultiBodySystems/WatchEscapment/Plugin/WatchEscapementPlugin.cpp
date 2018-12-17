#include <math.h>
#include <stdio.h>

#include <boost/math/quaternion.hpp>

static double inertia_z = 6.55929728e-05;
static double frequency =3.0;


static double stiffness= inertia_z * pow(2*M_PI*frequency,2);
static double angle_init = -0.2300258882727588;
//static double angle_init = -0.8;


extern "C" void externalForces(double t, double *f, unsigned int size_z, double *z){
  f[0]=0;
  f[1]=0;
  f[2]=-10;
}
extern "C" void externalMoment(double t,double *m, unsigned int size_z, double *z){
    m[0]=0;
    m[1]=0;
    m[2]=0;
}

extern "C" void internalForces(double t, double *q, double *v, double *f, unsigned int size_z,double *z){
  // Simple spring in z direction.
  f[0]=0;
  f[1]=0;
  f[2]=0.;
  // printf("internalForcesB1 :\n");
  // printf("f[0] = %e\t f[1] = %e\t, f[2]=%e\n",f[0],f[1],f[2]);
}

extern "C" void internalForcesB1_Jacq(double t, double *q, double *v, double *jac, unsigned int size_z,double *z){
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j<7; j++)
      jac[i+j*3]=0.0;
  }
  jac[2+2*3]=1e4;
  // printf("internalForcesB1_Jacq :\n");
  // printf("jac[2+2*3] = %e\n", jac[2+2*3]);
}

extern "C" void internalMomentsBalanceWheel(double t, double *q, double *v, double *m, unsigned int size_z,double *z){
  //  printf("internalMomentsB1 :\n");
  // Simple torsional spring around z axis
  // printf("q[3] = %e\n", q[3]);
  // printf("q[4] = %e\n", q[4]);
  // printf("q[5] = %e\n", q[5]);
  // printf("q[6] = %e\n", q[6]);

  double angle = 2*asin(q[6]);
  //printf("angle = %e\n", angle);
  //printf("stiffness = %e \n", stiffness);
  m[0]=0.0;
  m[1]=0.;
  m[2]=stiffness * (angle-angle_init);
  //printf("m[0] = %e\t m[1] = %e\t, m[2]=%e\n",m[0],m[1],m[2]);
}

extern "C" void internalMomentsBalanceWheel_Jacq(double t, double *q, double *v, double *jac, unsigned int size_z,double *z){
  //printf("internalMomentsB1_Jacq :\n");
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j<7; j++)
      jac[i+j*3]=0.0;
  }
  // printf("q[3] = %e\n", q[3]);
  // printf("q[4] = %e\n", q[4]);
  // printf("q[5] = %e\n", q[5]);
  // printf("q[6] = %e\n", q[6]);

  //double angle = 2*asin(q[6]);
  // printf("angle = %e\n", angle);
  jac[2+6*3]=stiffness * 2.0 / sqrt(1 - q[6]*q[6]) ;
  // printf("jac[3+3*3] = %e\n", jac[3+3*3]);
}

extern "C" void externalMomentEscapeWheel(double t,double *m, unsigned int size_z, double *z){
    m[0]=0;
    m[1]=0;
    m[2]=-0.001;
}


extern "C" void prescribedvelocityB1(double time, unsigned int sizeofprescribedvelocity, double *pv)
{
  /* the plugin implements v(t) = C + A cos(omega *t) */

  double C = -150.0 ;
  // double omega = M_PI / 2.0;
  // double A = 10.0;

  //pv[0] =  A * cos(omega * time*100.0);
  pv[0] =  C;
  //printf("prescribed velocity = %e\n", pv[0]);
}
