#ifdef _WIN32
#define SICONOS_EXPORT extern "C" __declspec(dllexport)
#else
#define SICONOS_EXPORT extern "C"
#endif
#include <stdio.h>
#include <math.h>

#define CONSTRAINT_0 -1.0

// extern "C" double input(double time)
// {
//   double u;
//   u = sin(50 * time);
//   return u;
// }


// output function s.t. y = h(x,lambda)
// Currently y = C*x+D*lambda
SICONOS_EXPORT void Computeh(   double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *y,
                                unsigned z_size, double *z)
{
    y[0] = CONSTRAINT_0*x[0];
}

SICONOS_EXPORT void ComputeJxh( double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *C,
                                unsigned z_size, double *z)
{
   for(unsigned int i = 0; i < x_size*size_lambda; i++)
   {
      C[i] = 0.0;
   }
   C[0] = CONSTRAINT_0;
}

SICONOS_EXPORT void ComputeJlh( double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *D,
                                unsigned z_size, double *z)
{
     for(unsigned int i = 0; i < size_lambda*size_lambda; i++)
   {
      D[i] = 0.0;
   }  
}


// input function s.t. r = g(x,lambda)
SICONOS_EXPORT void Computeg(   double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *r,
                                unsigned z_size, double *z)
{
    //printf("Compute  g(x,lambda) -- plugin version \n");

    r[0] = 0.0;
    r[1] = 0.0;
    r[2] =lambda[0]*(x[0] + 1);
}

// Jacobian of the inputs w.r.t. x
// Note that matrices are flattened in cloumn first mode
SICONOS_EXPORT void ComputeJxg( double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *K,
                                unsigned z_size, double *z)
{
   //printf("Compute  Jacobian of g wrt x -- plugin version \n"); 
   for(unsigned int i = 0; i < x_size*x_size; i++)
   {
      K[i] = 0.0;
   }
    K[2] = lambda[0]; // K[2,0] = K[2] in flat-mode  (column first)
    //K[6] = lambda[0]; // K[2,0] = K[6] in flat-mode (row first)
}

// Jacobian of the inputs w.r.t. lambda
// Note that matrices are flattened in cloumn first mode
SICONOS_EXPORT void ComputeJlg( double time, unsigned x_size, double *x,
                                unsigned size_lambda, double* lambda, double *B,
                                unsigned z_size, double *z)
{
    //printf("Compute  Jacobian of g wrt lambda -- plugin version \n");
    for(unsigned int i = 0; i < x_size*size_lambda; i++)
    {
        B[i] = 0.0;
    }
    B[2] = x[0]; // B[2,0] = B[2] in flat-mode     
}

