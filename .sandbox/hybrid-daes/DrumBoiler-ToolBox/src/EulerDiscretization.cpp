
#include "EulerDiscretization.hpp"


EulerDiscretization::EulerDiscretization(Model* model_, FiniteDifference* diff_method_, double h_){
            model = model_;
            diff_method = diff_method_;
            h = h_;

            system_lhs =  (double *)calloc(model->get_nb_eqs(),sizeof(double));
            system_rhs =  (double *)calloc(model->get_nb_eqs(),sizeof(double));
            system_lhs_jac = NM_create(NM_DENSE, model->get_nb_eqs(), model->get_nb_vars());
            system_rhs_jac = NM_create(NM_DENSE, model->get_nb_eqs(), model->get_nb_vars());  

}

EulerDiscretization::~EulerDiscretization(){
            delete(model);
            delete(diff_method);
            delete(system_rhs);
            delete(system_rhs_jac);
            delete(system_lhs_jac);
}

Model* EulerDiscretization::get_model(){
    return model;
}

FiniteDifference* EulerDiscretization::get_diff_method(){
    return diff_method;
}


void EulerDiscretization::Discr(void* env, int size, double *z, double * system){

    NumericalSimuInfo* info = (NumericalSimuInfo* )env;
    EulerDiscretization* euler  = info->discretization;
    Model* model = euler->get_model();
    FiniteDifference* diff_method = euler->get_diff_method();

    // printf("[implicitEulerDiscr] -- size = %d.\n",  size);
    // printf("[implicitEulerDiscr] -- h = %f.\n",  info->h);
    // printf("[implicitEulerDiscr] -- nbDiffEqs = %d.\n",  info->nb_diff_eqs);

    // could be hidden in function pointer if generalization is needed ... 
    model->doubleMixtureRHS(z, euler->system_rhs);

    // printf("[implicitEulerDiscr] -- rhs = ["); 
    // for(int i=0;i<size;i++){
    //   printf(" %f,",info->system_rhs[i]);
    // }
    // printf(" ]\n");
    // printf(" [implicitEulerDiscr] env = %p\n", env);
    // printf(" [implicitEulerDiscr] info = %p\n", info);

    // ptrDiffFunction diff = &DiscreteSystem::finiteDifference; // need to pass the object in the env.
    model->doubleMixtureLHS(env, diff_method, z, euler->system_lhs);

    // printf("[implicitEulerDiscr] -- lhs = ["); 
    // for(int i=0;i<size;i++){
    //   printf(" %f,",info->system_lhs[i]);
    // }
    // printf(" ]\n");


    /* /!\ Writing the discretization as lhs - "h"*rhs makes the Newton solver fail to find a descent gradient */
    int ndiff = model->get_nb_diffeqs();
    for(int i=0;i<size;i++){
        if(i<ndiff) 
        {  
            system[i] = -euler->system_lhs[i] + euler->h*euler->system_rhs[i];
        }
        else
        { 
             system[i] = -euler->system_lhs[i] + euler->system_rhs[i];
        }
            
    } 

}


void EulerDiscretization::DiscrJac(void* env, int size, double *z, NumericsMatrix * jacobian){

    NumericalSimuInfo* info = (NumericalSimuInfo*)env;

    EulerDiscretization* euler  = info->discretization;
    Model* model = euler->get_model();
    FiniteDifference* diff_method = euler->get_diff_method();

    // could be hidden in function pointer if generalization is needed ...
    model->doubleMixtureRHS_Jacobian(z, euler->system_rhs_jac);
    // printf("\n #### RHS JACOBIAN \n");
    // NM_display(info->system_rhs_jac);
    // printf("\n ####  \n");

    model->doubleMixtureLHS_Jacobian(env, diff_method, z, euler->system_lhs_jac);
    //  printf("\n #### LHS JACOBIAN \n");
    //  NM_display(info->system_lhs_jac);
    //  printf("\n ####  \n");

    // printf("ndiff_eq = %d", info->nb_diff_eqs);

    int ndiff = model->get_nb_diffeqs();
    for(int j=0;j<ndiff;j++){
        for(int i=0;i<size;i++){
            // printf("scaled_rhs(%d,%d) = % f\n",j,i,info->system_rhs_jac->matrix0[M2Ac(j,i,size)]*info->h);
            euler->system_rhs_jac->matrix0[M2Ac(j,i,size)] *= euler->h;
        }
    }
    // printf("\n #### scaled RHS JACOBIAN \n");
    // NM_display(info->system_rhs_jac);


    /* /!\ Writing the discretization as lhs - "h"*rhs makes the Newton solver fail to find a descent gradient */
    NumericsMatrix * temp = NM_add(-1.0, euler->system_lhs_jac, 1.0, euler->system_rhs_jac);
    // TODO slow copy need a by reference add function
    NM_copy(temp,jacobian);
    NM_clear(temp); 

    // NM_display(jacobian);

 }   