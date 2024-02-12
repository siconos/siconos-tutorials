// Define plugin functions
//=============================================================================================================
#ifdef _WIN32
#define SICONOS_EXPORT extern "C" __declspec(dllexport)
#else
#define SICONOS_EXPORT extern "C"
#endif
#include <cmath>
#include <stdio.h>
using namespace std;
//===========================================================================================================
extern "C" double l;

//1. Plugin function to calculate the gap function h1 at contact point 1 and h2 at contact point 2
SICONOS_EXPORT void h1(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* y, unsigned int sizeOfZ, double* z)
{
  double l = 1.0;
  y[0] = q[1] - 0.5 * l * cos(q[2]);
}

//2. Plugin function to calculate the gradient of h1 (G1) and the one of h2 (G2)
SICONOS_EXPORT void G1(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* Wn, unsigned int sizeOfZ, double* z)
{
  double l = 1.0;
  Wn[0] = 0.0;
  Wn[1] = -1.0;
  
  Wn[2] = 1.0;
  Wn[3] = 0.0;
  
  Wn[4] = + 0.5 * l * sin(q[2]);
  Wn[5] = - 0.5 * l * cos(q[2]);
}
