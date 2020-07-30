#include "discrete_system.hpp"

// finite difference discretisation of the differentiation operator (scaled by 'h')
double finiteDifference(void* env, double z_i, int i){

    NumericalSimuInfo* info = (NumericalSimuInfo*)env;
    return z_i-info->z_prev[i];

}

// partial derivartive of finite difference 
double finiteDifferenceDz(void* env, double z_i, int i){
    return 1.0;
}

void implicitEulerDiscr(void* env, int size, double *z, double * system){

    NumericalSimuInfo* info = (NumericalSimuInfo*)env;

    // printf("[implicitEulerDiscr] -- size = %d.\n",  size);
    // printf("[implicitEulerDiscr] -- h = %f.\n",  info->h);
    // printf("[implicitEulerDiscr] -- nbDiffEqs = %d.\n",  info->nb_diff_eqs);

    // could be hidden in function pointer if generalization is needed ... 
    doubleMixtureRHS(z, info->system_rhs);

    // printf("[implicitEulerDiscr] -- rhs = ["); 
    // for(int i=0;i<size;i++){
    //   printf(" %f,",info->system_rhs[i]);
    // }
    // printf(" ]\n");
    // printf(" [implicitEulerDiscr] env = %p\n", env);
    // printf(" [implicitEulerDiscr] info = %p\n", info);

    doubleMixtureLHS(env, &finiteDifference, z, info->system_lhs);

    // printf("[implicitEulerDiscr] -- lhs = ["); 
    // for(int i=0;i<size;i++){
    //   printf(" %f,",info->system_lhs[i]);
    // }
    // printf(" ]\n");


    /* /!\ Writing the discretization as lhs - "h"*rhs makes the Newton solver fail to find a descent gradient */
    for(int i=0;i<size;i++){
        if(i<info->nb_diff_eqs) 
        {  
            system[i] = -info->system_lhs[i] + info->h*info->system_rhs[i];
        }
        else
        { 
             system[i] = -info->system_lhs[i] + info->system_rhs[i];
        }
            
    } 

}

void implicitEulerDiscrJac(void* env, int size, double *z, NumericsMatrix * jacobian){
    NumericalSimuInfo* info = (NumericalSimuInfo*)env;
    // could be hidden in function pointer if generalization is needed ...
    doubleMixtureRHS_Jacobian(z, size, info->system_rhs_jac);
    // printf("\n #### RHS JACOBIAN \n");
    // NM_display(info->system_rhs_jac);
    // printf("\n ####  \n");

    doubleMixtureLHS_Jacobian(env, size, &finiteDifference, &finiteDifferenceDz, z, info->system_lhs_jac);
    //  printf("\n #### LHS JACOBIAN \n");
    //  NM_display(info->system_lhs_jac);
    //  printf("\n ####  \n");

    // printf("ndiff_eq = %d", info->nb_diff_eqs);
    for(int j=0;j<info->nb_diff_eqs;j++){
        for(int i=0;i<size;i++){
            // printf("scaled_rhs(%d,%d) = % f\n",j,i,info->system_rhs_jac->matrix0[M2Ac(j,i,size)]*info->h);
            info->system_rhs_jac->matrix0[M2Ac(j,i,size)] *= info->h;
        }
    }
    // printf("\n #### scaled RHS JACOBIAN \n");
    // NM_display(info->system_rhs_jac);


    /* /!\ Writing the discretization as lhs - "h"*rhs makes the Newton solver fail to find a descent gradient */
    NumericsMatrix * temp = NM_add(-1.0, info->system_lhs_jac, 1.0, info->system_rhs_jac);
    // TODO slow copy need a by reference add function
    NM_copy(temp,jacobian);
    NM_clear(temp); 

    // NM_display(jacobian);

 }   

 void initialize_simu(  NumericalSimuInfo* simu_info, const double &P0, const double &M10, const double &x10,
                        const double &T10, const double &M20, const double &x20, const double &T20){

    double* z0 = (double *)calloc(simu_info->nb_eqs,sizeof(double)); 
    /////////////////////////// INITIALISATION MIXTURE 1 (VAPOUR) //////////////////////
    z0[M1] = M10;
    z0[x1] = x10;
    z0[T1] = T10;
    z0[P] = P0;
    z0[rv1] = rhov(z0[P],z0[T1]);   // kg/m^3
    z0[rl1] = rhol(z0[P],z0[T1]);   // kg/m^3
    z0[hv1] = h_v(z0[P],z0[T1]);    // J/kg
    z0[hl1] = h_l(z0[P],z0[T1]);    // J/kg

    z0[V1] = z0[M1]*z0[x1]/z0[rv1] + z0[M1]*(1-z0[x1])/z0[rl1];

    if (z0[P]-psat(z0[T1]) < 0.0)
      z0[DP1] = -1.0*(z0[P]-psat(z0[T1]));  // [P-Psat]⁻ in Pa
    else   
      z0[DP1] = 0.0;  // [P-Psat]⁻ in Pa

    z0[U1] = (1.0-z0[x1])*(z0[hl1] - z0[P]/z0[rl1]) + z0[x1]*(z0[hv1] - z0[P]/z0[rv1]); // J
    

    /////////////////////////// INITIALISATION MIXTURE 2 (LIQUID) //////////////////////

    z0[M2] = M20;
    z0[x2] = x20;
    z0[T2] = T20;
    z0[rv2] = rhov(z0[P],z0[T2]);   // kg/m^3
    z0[rl2] = rhol(z0[P],z0[T2]);   // kg/m^3
    z0[hv2] = h_v(z0[P],z0[T2]);    // J/kg
    z0[hl2] = h_l(z0[P],z0[T2]);    // J/kg

    z0[V2] = z0[M2]*z0[x2]/z0[rv2] + z0[M2]*(1-z0[x2])/z0[rl2];

    if (z0[P]-psat(z0[T2]) < 0.0)
      z0[DP2] = -1.0*(z0[P]-psat(z0[T2]));  // [P-Psat]⁻ in Pa
    else   
      z0[DP2] = 0.0;  // [P-Psat]⁻ in Pa

    z0[U2] = (1.0-z0[x2])*(z0[hl2] - z0[P]/z0[rl2]) + z0[x2]*(z0[hv2] - z0[P]/z0[rv2]); // J


    simu_info->z_prev = z0;

    /// matrices initialisation
    simu_info->nb_compl_eqs = 4;
    simu_info->nb_diff_eqs = 4;   
    simu_info->system_lhs =  (double *)calloc(simu_info->nb_eqs,sizeof(double));
    simu_info->system_rhs =  (double *)calloc(simu_info->nb_eqs,sizeof(double));
    simu_info->system_lhs_jac = NM_create(NM_DENSE, simu_info->nb_eqs, simu_info->nb_eqs);
    simu_info->system_rhs_jac = NM_create(NM_DENSE, simu_info->nb_eqs, simu_info->nb_eqs);  

    std::cout << "\nInitial condition, z0 = [";
    for(int i=0;i<simu_info->nb_eqs;i++){
      std::cout << z0[i] << ", ";
    }
    std::cout << " ]" << std::endl;
}
