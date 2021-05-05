#ifndef EULERDISCRETIZATION_HPP
#define EULERDISCRETIZATION_HPP

#include "Model.hpp"

class EulerDiscretization
{
    public:

        typedef struct
        {
            EulerDiscretization* discretization;

        }NumericalSimuInfo; // TODO passer directement EulerDiscretization object in the env.

       EulerDiscretization(Model* model_, FiniteDifference* diff_method_, double h_);

        ~EulerDiscretization();

        static void Discr(void* env, int size, double *z, double * system); 

        static void DiscrJac(void* env, int size, double *z, NumericsMatrix * jacobian);

        Model* get_model();
        FiniteDifference* get_diff_method();

        // put in protected and add getter/setter
        double* system_rhs; 
        double* system_lhs;
        NumericsMatrix* system_rhs_jac;
        NumericsMatrix* system_lhs_jac;
        double h;

    protected:

        Model* model;    
        FiniteDifference* diff_method;
        // TODO extend to double** (NumericsMatrix**) in case of theta/multistep methods (or use c++ std::array<T>)

        // time step
};


#endif //EULERDISCRETIZATION_HPP