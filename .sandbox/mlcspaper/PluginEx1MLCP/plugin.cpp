#include <math.h>
#include <stdio.h>

#if defined(_MSC_VER)
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

/* 2d matrix conversion in 1d array 
(Column first, indexes starting at 0) */
int M2Ac(int i, int j, int nb_line)
{
    int index = nb_line*j + i;
    return index;
}

double u1(double t)
{
    return sin(t);
}    
double u2(double t)
{
    return cos(10*t);
}   

extern "C" DLLEXPORT void computef(double time, unsigned int sizeOfX, double* x, double* f, unsigned int sizeZ, double* z);
extern "C" DLLEXPORT void computef(double time, unsigned int sizeOfX, double* x, double* f, unsigned int sizeZ, double* z)
{
    f[0] = -2.12e3*x[0] + 0.301e4*u1(time);
    f[1] = 0;
}

extern "C" DLLEXPORT void computeJacobianfx(double time, unsigned int sizeOfX, double* x, double* jacob, unsigned int sizeZ, double* z);
extern "C" DLLEXPORT void computeJacobianfx(double time, unsigned int sizeOfX, double* x, double* jacob, unsigned int sizeZ, double* z)
{
	jacob[M2Ac(0,0,sizeOfX)] = -2.12e3; 
}

//////////////// Mixed Complementarity Constraint ///////////////////////////
// relation w(x,lambda) _|_ lambda
extern "C" void computeh_mlcp(   double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *f,
                                unsigned z_size, double *z)
{  
    f[0] = x[1] + 1e-4*lambda[1] + 1e-4*u2(time);
    f[1] = -0.707*x[0] - lambda[0];  
}


extern "C" void computeJxh_mlcp( double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *jacob,
                                unsigned z_size, double *z)
{
   jacob[M2Ac(0,1,size_lambda)] = 1.0;
   jacob[M2Ac(1,0,size_lambda)] = -0.707;
}

extern "C" void computeJlh_mlcp( double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *jacob,
                                unsigned z_size, double *z)
{   
    jacob[M2Ac(0,1,size_lambda)] = 1e-4;
    jacob[M2Ac(1,0,size_lambda)] = -1.0;
}

// input function s.t. r = g(x,lambda)
extern "C" void computeg_mlcp(   double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *f,
                                unsigned z_size, double *z)
{   
    f[0] = -3.0e-3*lambda[1];
    f[1] = lambda[0];
} 

// Jacobian of the inputs w.r.t. x
extern "C" void computeJxg_mlcp( double time, unsigned x_size, double *x, 
                                unsigned size_lambda, double* lambda, double *jacob,
                                unsigned z_size, double *z)
{ 
}

// Jacobian of the inputs w.r.t. lambda
extern "C" void computeJlg_mlcp( double time, unsigned x_size, double *x,
                                unsigned size_lambda, double* lambda, double *jacob,
                                unsigned z_size, double *z)
{   
    jacob[M2Ac(0,1,x_size)] = -3.0e-3; 
    jacob[M2Ac(1,0,x_size)] = 1.0;                     
}