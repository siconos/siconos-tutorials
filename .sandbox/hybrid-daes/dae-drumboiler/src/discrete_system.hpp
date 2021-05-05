#ifndef DISCRETE_SYSTEM_HPP
#define DISCRETE_SYSTEM_HPP

#include "model_common.hpp"
#include "model.hpp"

/* Information needed for the numerical simulation*/ 
typedef struct
{
    int nb_eqs;     // square system assumed
    int nb_diff_eqs;
    int nb_compl_eqs; 
    // TODO extended to double** (NumericsMatrix**) in case of theta/multistep methods (or use c++ std::array<T>)
    double* system_rhs; 
    double* system_lhs;
    NumericsMatrix* system_rhs_jac;
    NumericsMatrix* system_lhs_jac;
    double* z_prev; 
    // time step
    double h;
    double Tfinal;
}NumericalSimuInfo;


double finiteDifference(void* env, double z_i, int i);

double finiteDifferenceDz(void* env, double z_i, int i);

void implicitEulerDiscr(void* env, int size, double *z, double * system);

void implicitEulerDiscrJac(void* env, int size, double *z, NumericsMatrix * jacobian);

void initialize_simu(NumericalSimuInfo* simu_info, const double &P0, const double &M10, const double &x10, const double &T10, const double &M20, const double &x20, const double &T20);

#endif // DISCRETE_SYSTEM_HPP