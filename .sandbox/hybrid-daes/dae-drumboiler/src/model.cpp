#include "model.hpp"

//TODO use sparse matrix format at some point ...

static ModelConstants CTE = {
Ccond : 0.01,
Cvap : 0.09,
Kvl : 10,
Qext : -1000,//-1000, //
Pext : 1060740.0 // 106074.0 // 101325.0// 106074.0 // 1 ATM
};

/* Model of condensation mass transfert */
double massCond(const double &M_, const double &x_){
    return CTE.Ccond*M_*(1-x_);
} 
double massConddx(const double &M_, const double &x_){
    return -1.0*CTE.Ccond*M_;
} 
double massConddM(const double &M_, const double &x_){
    return CTE.Ccond*(1-x_);
} 


/* Model of evaporatation mass transfert */
double massEvap(const double &M_, const double &x_){
    return CTE.Cvap*M_*x_;
} 
double massEvapdx(const double &M_, const double &x_){
    return CTE.Cvap*M_;
} 
double massEvapdM(const double &M_, const double &x_){
    return CTE.Cvap*x_;
} 


/* Model of condensation energy transfert */
double energCond(const double &M_, const double &x_, const double &h_){
    return CTE.Ccond*h_*M_*(1-x_);
}
double energConddM(const double &M_, const double &x_, const double &h_){
    return CTE.Ccond*h_*(1-x_);
}
double energConddx(const double &M_, const double &x_, const double &h_){
    return -1.0*CTE.Ccond*h_*M_;
}
double energConddh(const double &M_, const double &x_, const double &h_){
    return CTE.Ccond*M_*(1-x_);
}

/* Model of evaportation energy transfert */
double energEvap(const double &M_, const double &x_, const double &h_){
    return CTE.Cvap*h_*M_*x_;
}
double energEvapdM(const double &M_, const double &x_, const double &h_){
    return CTE.Cvap*h_*x_;
}
double energEvapdx(const double &M_, const double &x_, const double &h_){
    return CTE.Cvap*h_*M_;
}
double energEvapdh(const double &M_, const double &x_, const double &h_){
    return CTE.Cvap*M_*x_;
}


/* Model of external energy transfert  - Constant = Qext 
   TODO remove if not used.
   TODO Add other external interaction 
*/
double energIn(void){
    return CTE.Qext;
}

/* Model of energy conduction at mixture interface 
        FROM phase 1 TO phase 2 
>0 if energy is transfered to phase 2 (T1>T2)*/
double tempConduc_1_2(const double &T1_, const double &T2_){
    return CTE.Kvl*(T1_-T2_);
}
double tempConduc_1_2dT1(const double &T1_, const double &T2_){
    return CTE.Kvl;
}
double tempConduc_1_2dT2(const double &T1_, const double &T2_){
    return -1.0*CTE.Kvl;
}

/* Volume equation of mixture */
double volMix(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return x_*vol(rv_,M_) + (1.0-x_)*vol(rl_,M_);
}
double volMixdrl(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return (1.0-x_)*voldr(rl_,M_);
}
double volMixdrv(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return x_*voldr(rv_,M_);
}
double volMixdM(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return x_*voldM(rv_,M_) + (1.0-x_)*voldM(rl_,M_);
}
double volMixdx(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return vol(rv_,M_) - vol(rl_,M_);
}


/* specific internal energy of mixture */ 
double uMix(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return x_*u(P_,hv_,rv_) + (1.0-x_)*u(P_,hl_,rl_);
}
double uMixdhl(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return (1.0-x_)*udh(P_,hl_,rl_);
}
double uMixdhv(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return x_*udh(P_,hv_,rv_);
}
double uMixdrl(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return (1.0-x_)*udr(P_,hl_,rl_);
}
double uMixdrv(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return x_*udr(P_,hv_,rv_);
}
double uMixdP(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return x_*udP(P_,hv_,rv_) + (1.0-x_)*udP(P_,hl_,rl_);
}
double uMixdx(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return u(P_,hv_,rv_) - u(P_,hl_,rl_);
}

/*Complete double mixture mixed complementarity system (MCS) Right-Hand Side 
    ### System size. 
                - 21 Equations.
                - 21 Variables. 
                - 4 Vomplementarity Variables 

    ### System structure.   
    Each equation i is of the form: 
                    LHS[i](vars,diff(vars)) = RHS[i](vars)
    Rq:
    diff(var) is the differentiation of var wrt time
    LHS[i] = 0 for complementarity equations.                             
*/   
void doubleMixtureRHS(double* z, double* rhs){


                        /* Differential Equations */

    // Vapor Mixture Mass Balance Equation
    rhs[0] = massEvap(z[M2],z[x2]) - massCond(z[M1],z[x1]);  
    // Vapor Mixture Energy Balance Equation
    rhs[1] = energEvap(z[M2],z[x2],z[hv2]) - tempConduc_1_2(z[T1],z[T2]) - energCond(z[M1],z[x1],z[hl1]) + CTE.Qext ;    
    // Liquid Mixture Mass Balance Equation
    rhs[2] = massCond(z[M1],z[x1]) - massEvap(z[M2],z[x2]); 
    // Liquid Mixture Energy Balance Equation
    rhs[3] = energCond(z[M1],z[x1],z[hl1]) + tempConduc_1_2(z[T1],z[T2]) -energEvap(z[M2],z[x2],z[hv2]) ;   

                                /* MIXTURE 1 -- VAPOR AND DROPPLETS */
    // Vapor Mixture Volume Equation
    rhs[4] = volMix(z[rl1],z[rv1],z[M1],z[x1]);
    // Vapor Mixture Specific Internal Energy Equation
    rhs[5] = uMix(z[hl1],z[hv1],z[rl1],z[rv1],z[P],z[x1]);
    // Vapor Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor         
    rhs[6] = rhov(z[P],z[T1]);
    rhs[7] = rhol(z[P],z[T1]);
    rhs[8] = h_v(z[P],z[T1]); 
    rhs[9] = h_l(z[P],z[T1]);   

                            /* MIXTURE 2 -- LIQUID AND BUBBLES */
                  
    // Liquid Mixture Volume Equation
    rhs[10] = volMix(z[rl2],z[rv2],z[M2],z[x2]);

    // Liquid Mixture Specific Internal Energy Equation
    rhs[11] = uMix(z[hl2],z[hv2],z[rl2],z[rv2],z[P],z[x2]); 
    // Liquid Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor     
    rhs[12] = rhov(z[P],z[T2]);
    rhs[13] = rhol(z[P],z[T2]);
    rhs[14] = h_v(z[P],z[T2]); 
    rhs[15] = h_l(z[P],z[T2]);  

    /////// EXTERNAL /////////
    rhs[16] = CTE.Pext;                  // Fixed pressure

    /////////// COMPLEMENTARITY ///////////////// 
    // Vapor mixture complementarity
    rhs[17] = z[P] - psat(z[T1]) + z[DP1];
    rhs[18] = 1.0-z[x1];
    // Liquid mixture complementarity
    rhs[19] = z[P] - psat(z[T2]) + z[DP2];
    rhs[20] = 1.0-z[x2];
}

// Jacobian of the MCS right hand side
void doubleMixtureRHS_Jacobian(double* z, int sys_size, NumericsMatrix* rhs_jac){

                    /* Differential Equations */

    /* Vapor Mixture Mass Balance Equation
         massEvap(z[M2],z[x2]) - massCond(z[M1],z[x1])
    */
    rhs_jac->matrix0[M2Ac(0,M1,sys_size)] = -1.0*massConddM(z[M1],z[x1]);
    rhs_jac->matrix0[M2Ac(0,x1,sys_size)] = -1.0*massConddx(z[M1],z[x1]);
    rhs_jac->matrix0[M2Ac(0,M2,sys_size)] = massEvapdM(z[M2],z[x2]);
    rhs_jac->matrix0[M2Ac(0,x2,sys_size)] = massEvapdx(z[M2],z[x2]);

    /* Vapor Mixture Energy Balance Equation
            energEvap(z[M2],z[x2],z[hv2]) - tempConduc_1_2(z[T1],z[T2]) - energCond(z[M1],z[x1],z[hl1])
    */  
    rhs_jac->matrix0[M2Ac(1,M1,sys_size)] = -1.0*energConddM(z[M1],z[x1],z[hl1]) ; 
    rhs_jac->matrix0[M2Ac(1,M2,sys_size)] = energEvapdM(z[M2],z[x2],z[hv2]) ;    
    rhs_jac->matrix0[M2Ac(1,x1,sys_size)] = -1.0*energConddx(z[M1],z[x1],z[hl1]) ;    
    rhs_jac->matrix0[M2Ac(1,x2,sys_size)] = energEvapdx(z[M2],z[x2],z[hv2]) ;    
    rhs_jac->matrix0[M2Ac(1,hl1,sys_size)] = -1.0*energConddh(z[M1],z[x1],z[hl1]) ;    
    rhs_jac->matrix0[M2Ac(1,hv2,sys_size)] = energEvapdh(z[M2],z[x2],z[hv2]) ;    
    rhs_jac->matrix0[M2Ac(1,T1,sys_size)] = -1.0*tempConduc_1_2dT1(z[T1],z[T2]);    
    rhs_jac->matrix0[M2Ac(1,T2,sys_size)] = -1.0*tempConduc_1_2dT2(z[T1],z[T2]);    

    /* Liquid Mixture Mass Balance Equation
        massCond(z[M1],z[x1]) - massEvap(z[M2],z[x2]); 
    */
    rhs_jac->matrix0[M2Ac(2,M1,sys_size)] = massConddM(z[M1],z[x1]);
    rhs_jac->matrix0[M2Ac(2,x1,sys_size)] = massConddx(z[M1],z[x1]);
    rhs_jac->matrix0[M2Ac(2,M2,sys_size)] = -1.0*massEvapdM(z[M2],z[x2]);
    rhs_jac->matrix0[M2Ac(2,x2,sys_size)] = -1.0*massEvapdx(z[M2],z[x2]);


    /* Liquid Mixture Energy Balance Equation
        energCond(z[M1],z[x1],z[hl1]) + tempConduc_1_2(z[T1],z[T2]) -energEvap(z[M2],z[x2],z[hv2]) + CTE.Qext;     
    */
    rhs_jac->matrix0[M2Ac(3,M1,sys_size)] = energConddM(z[M1],z[x1],z[hl1]) ;     
    rhs_jac->matrix0[M2Ac(3,x1,sys_size)] = energConddx(z[M1],z[x1],z[hl1]) ;     
    rhs_jac->matrix0[M2Ac(3,hl1,sys_size)] = energConddh(z[M1],z[x1],z[hl1]) ;     
    rhs_jac->matrix0[M2Ac(3,T1,sys_size)] = tempConduc_1_2dT1(z[T1],z[T2]);     
    rhs_jac->matrix0[M2Ac(3,T2,sys_size)] = tempConduc_1_2dT2(z[T1],z[T2]);     
    rhs_jac->matrix0[M2Ac(3,M2,sys_size)] = -1.0*energEvapdM(z[M2],z[x2],z[hv2]) ;     
    rhs_jac->matrix0[M2Ac(3,x2,sys_size)] = -1.0*energEvapdx(z[M2],z[x2],z[hv2]) ;     
    rhs_jac->matrix0[M2Ac(3,hv2,sys_size)] = -1.0*energEvapdh(z[M2],z[x2],z[hv2]) ;  
 
                    /* MIXTURE 1 -- VAPOR AND DROPPLETS */

    /* Vapor Mixture Volume Equation
         volMix(z[rl1],z[rv1],z[M1],z[x1]);
    */
    rhs_jac->matrix0[M2Ac(4,M1,sys_size)] =  volMixdM(z[rl1],z[rv1],z[M1],z[x1]);
    rhs_jac->matrix0[M2Ac(4,x1,sys_size)] =  volMixdx(z[rl1],z[rv1],z[M1],z[x1]);
    rhs_jac->matrix0[M2Ac(4,rl1,sys_size)] =  volMixdrl(z[rl1],z[rv1],z[M1],z[x1]);
    rhs_jac->matrix0[M2Ac(4,rv1,sys_size)] =  volMixdrv(z[rl1],z[rv1],z[M1],z[x1]);


    /* Vapor Mixture Specific Internal Energy Equation
         uMix(z[hl1],z[hv1],z[rl1],z[rv1],z[P],z[x1]);
    */
    rhs_jac->matrix0[M2Ac(5,hl1,sys_size)] = uMixdhl(z[hl1],z[hv1],z[rl1],z[rv1],z[P],z[x1]);
    rhs_jac->matrix0[M2Ac(5,hv1,sys_size)] = uMixdhv(z[hl1],z[hv1],z[rl1],z[rv1],z[P],z[x1]);
    rhs_jac->matrix0[M2Ac(5,rl1,sys_size)] = uMixdrl(z[hl1],z[hv1],z[rl1],z[rv1],z[P],z[x1]);
    rhs_jac->matrix0[M2Ac(5,rv1,sys_size)] = uMixdrv(z[hl1],z[hv1],z[rl1],z[rv1],z[P],z[x1]);
    rhs_jac->matrix0[M2Ac(5,P,sys_size)] = uMixdP(z[hl1],z[hv1],z[rl1],z[rv1],z[P],z[x1]);
    rhs_jac->matrix0[M2Ac(5,x1,sys_size)] = uMixdx(z[hl1],z[hv1],z[rl1],z[rv1],z[P],z[x1]);

    /* Vapor Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
     one equation for each constituant of the vapor */     
    /*  rhov(z[P],z[T1]);  */   
    rhs_jac->matrix0[M2Ac(6,P,sys_size)] = rhovdP(z[P],z[T1]);
    rhs_jac->matrix0[M2Ac(6,T1,sys_size)] = rhovdT(z[P],z[T1]);

    /*  rhol(z[P],z[T1]);  */   
    rhs_jac->matrix0[M2Ac(7,P,sys_size)] = rholdP(z[P],z[T1]);
    rhs_jac->matrix0[M2Ac(7,T1,sys_size)] = rholdT(z[P],z[T1]);

    /*  h_v(z[P],z[T1]);  */ 
    rhs_jac->matrix0[M2Ac(8,P,sys_size)] = h_vdP(z[P],z[T1]); 
    rhs_jac->matrix0[M2Ac(8,T1,sys_size)] = h_vdT(z[P],z[T1]); 

    /*  h_l(z[P],z[T1]);  */
    rhs_jac->matrix0[M2Ac(9,P,sys_size)] = h_ldP(z[P],z[T1]); 
    rhs_jac->matrix0[M2Ac(9,T1,sys_size)] = h_ldT(z[P],z[T1]);  

                /* MIXTURE 2 -- LIQUID AND BUBBLES */
    
        /* Liquid Mixture Volume Equation
        volMix(z[rl2],z[rv2],z[M2],z[x2]);
    */
    rhs_jac->matrix0[M2Ac(10,rl2,sys_size)] = volMixdrl(z[rl2],z[rv2],z[M2],z[x2]);
    rhs_jac->matrix0[M2Ac(10,rv2,sys_size)] = volMixdrv(z[rl2],z[rv2],z[M2],z[x2]);
    rhs_jac->matrix0[M2Ac(10,M2,sys_size)] = volMixdM(z[rl2],z[rv2],z[M2],z[x2]);
    rhs_jac->matrix0[M2Ac(10,x2,sys_size)] = volMixdx(z[rl2],z[rv2],z[M2],z[x2]);


    /* Liquid Mixture Specific Internal Energy Equation
        uMix(z[hl2],z[hv2],z[rl2],z[rv2],z[P],z[x2]); 
    */
    rhs_jac->matrix0[M2Ac(11,hl2,sys_size)] = uMixdhl(z[hl2],z[hv2],z[rl2],z[rv2],z[P],z[x2]); 
    rhs_jac->matrix0[M2Ac(11,hv2,sys_size)] = uMixdhv(z[hl2],z[hv2],z[rl2],z[rv2],z[P],z[x2]); 
    rhs_jac->matrix0[M2Ac(11,rl2,sys_size)] = uMixdrl(z[hl2],z[hv2],z[rl2],z[rv2],z[P],z[x2]); 
    rhs_jac->matrix0[M2Ac(11,rv2,sys_size)] = uMixdrv(z[hl2],z[hv2],z[rl2],z[rv2],z[P],z[x2]); 
    rhs_jac->matrix0[M2Ac(11,P,sys_size)] = uMixdP(z[hl2],z[hv2],z[rl2],z[rv2],z[P],z[x2]); 
    rhs_jac->matrix0[M2Ac(11,x2,sys_size)] = uMixdx(z[hl2],z[hv2],z[rl2],z[rv2],z[P],z[x2]);
    
    
    /* Liquid Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
     one equation for each constituant of the vapor     */
    /*  rhov(z[P],z[T2]); */
    rhs_jac->matrix0[M2Ac(12,P,sys_size)] = rhovdP(z[P],z[T2]);
    rhs_jac->matrix0[M2Ac(12,T2,sys_size)] = rhovdT(z[P],z[T2]);

    /*  rhol(z[P],z[T2]); */
    rhs_jac->matrix0[M2Ac(13,P,sys_size)] = rholdP(z[P],z[T2]);
    rhs_jac->matrix0[M2Ac(13,T2,sys_size)] = rholdT(z[P],z[T2]);

    /*  h_v(z[P],z[T2]); */ 
    rhs_jac->matrix0[M2Ac(14,P,sys_size)] = h_vdP(z[P],z[T2]);
    rhs_jac->matrix0[M2Ac(14,T2,sys_size)] = h_vdT(z[P],z[T2]);

    /*  h_l(z[P],z[T2]);  */
    rhs_jac->matrix0[M2Ac(15,P,sys_size)] = h_ldP(z[P],z[T2]);
    rhs_jac->matrix0[M2Ac(15,T2,sys_size)] = h_ldT(z[P],z[T2]);

    /////// EXTERNAL /////////
    /* Cte --> line [16] of rhs_jac full of 0    */

    /////////// COMPLEMENTARITY ///////////////// 
    /* (Vapor) mixture 1 complementarity */
    /* z[P] - psat(z[T1]) + z[DP1]; */
    rhs_jac->matrix0[M2Ac(17,P,sys_size)] = 1.0;
    rhs_jac->matrix0[M2Ac(17,T1,sys_size)] = -1.0*psatdT(z[T1]);
    rhs_jac->matrix0[M2Ac(17,DP1,sys_size)] = 1.0;

    /* 1.0-z[x1];*/
    rhs_jac->matrix0[M2Ac(18,x1,sys_size)] = -1.0;

    /* (Liquid) mixture 2 complementarity */
    /* z[P] - psat(z[T2]) + z[DP2];*/
    rhs_jac->matrix0[M2Ac(19,P,sys_size)] = 1.0;
    rhs_jac->matrix0[M2Ac(19,T2,sys_size)] = -1.0*psatdT(z[T2]);
    rhs_jac->matrix0[M2Ac(19,DP2,sys_size)] = 1.0;

    /* 1.0-z[x2];*/
    rhs_jac->matrix0[M2Ac(20,x2,sys_size)] = -1.0;
}

 /*Complete double mixture mixed complementarity system (MCS) Left-Hand Side 
    ### System size. 
                - 21 Equations.
                - 21 Variables. 
                - 4 Vomplementarity Variables 

    ### System structure.   
    Each equation i is of the form: 
                    LHS[i](vars,diff(vars)) = RHS[i](vars)
    Rq:
        - Differential Equations always comes first (needed for numerical differentiation step)
        - Complementarity Equations always come Last (needed for siconos MCP solver)
        - diff(var) is the differentiation of var wrt time
        - LHS[i] = 0 for complementarity equations. 
        - LHS[i] may be chosen = 0 for algebraic equations  

    ### Differentiation
                - ptrDifferentiationFunction is a pointer for numerical Differentiation function
                - env : contains any information needed for the numerical differentiation (previous step, etc ...)                        
*/   
void doubleMixtureLHS(void* env, ptrDifferentiationFunction diff, double* z, double* lhs){

                      /* Differential Equations */

    // Vapor Mixture Mass Balance Equation
    lhs[0] = diff(env,z[M1],M1);  
    // Vapor Mixture Energy Balance Equation
    lhs[1] = z[M1]*diff(env,z[U1],U1) + z[U1]*diff(env,z[M1],M1); 
    // Liquid Mixture Mass Balance Equation
    lhs[2] = diff(env,z[M2],M2);
    // Liquid Mixture Energy Balance Equation
    lhs[3] = z[M2]*diff(env,z[U2],U2) + z[U2]*diff(env,z[M2],M2);       

                        /* MIXTURE 1 -- VAPOR AND DROPPLETS */

    // Vapor Mixture Volume Equation
    lhs[4] = z[V1];
    // Vapor Mixture Specific Internal Energy Equation
    lhs[5] = z[U1]; 
    // Vapor Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor         
    lhs[6] = z[rv1];
    lhs[7] = z[rl1];
    lhs[8] = z[hv1]; 
    lhs[9] = z[hl1];   

                            /* MIXTURE 2 -- LIQUID AND BUBBLES */
     
    // Liquid Mixture Volume Equation
    lhs[10] = z[V2];

    // Liquid Mixture Specific Internal Energy Equation
    lhs[11] = z[U2]; 
    // Liquid Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor     
    lhs[12] = z[rv2];
    lhs[13] = z[rl2];
    lhs[14] = z[hv2]; 
    lhs[15] = z[hl2];  

    /////// EXTERNAL /////////
    lhs[16] = z[P];  // Fixed pressure

    /////////// COMPLEMENTARITY ///////////////// 
    // Vapor mixture complementarity
    lhs[17] = 0;
    lhs[18] = 0;
    // Liquid mixture complementarity
    lhs[19] = 0;
    lhs[20] = 0;
}

/* Jacobian of the MCS Left Hand side
    ptrDifferentiationFunctionDz is a pointer to a function computing the partial derivative of a differentiated variable
*/
void doubleMixtureLHS_Jacobian(void* env, int sys_size, ptrDifferentiationFunction diff, ptrDifferentiationFunctionDz diff_dz, double* z, NumericsMatrix* lhs_jac){

                            /* Differential Equations */

    // Vapor Mixture Mass Balance Equation
    lhs_jac->matrix0[M2Ac(0,M1,sys_size)] = diff_dz(env,z[M1],M1);

    // Vapor Mixture Energy Balance Equation
    lhs_jac->matrix0[M2Ac(1,M1,sys_size)] = diff(env,z[U1],U1) + z[U1]*diff_dz(env,z[M1],M1); 
    lhs_jac->matrix0[M2Ac(1,U1,sys_size)] = z[M1]*diff_dz(env,z[U1],U1) + diff(env,z[M1],M1);

    // Liquid Mixture Mass Balance Equation
    lhs_jac->matrix0[M2Ac(2,M2,sys_size)] = diff_dz(env,z[M2],M2);

    // Liquid Mixture Energy Balance Equation
    lhs_jac->matrix0[M2Ac(3,M2,sys_size)] = diff(env,z[U2],U2) + z[U2]*diff_dz(env,z[M2],M2); 
    lhs_jac->matrix0[M2Ac(3,U2,sys_size)] = z[M2]*diff_dz(env,z[U2],U2) + diff(env,z[M2],M2); 

                        /* MIXTURE 1 -- VAPOR AND DROPPLETS */

    // Vapor Mixture Volume Equation
    lhs_jac->matrix0[M2Ac(4,V1,sys_size)] = 1.0;

    // Vapor Mixture Specific Internal Energy Equation
    lhs_jac->matrix0[M2Ac(5,U1,sys_size)] = 1.0; 

    // Vapor Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor         
    lhs_jac->matrix0[M2Ac(6,rv1,sys_size)] = 1.0;

    lhs_jac->matrix0[M2Ac(7,rl1,sys_size)] = 1.0;

    lhs_jac->matrix0[M2Ac(8,hv1,sys_size)] = 1.0; 

    lhs_jac->matrix0[M2Ac(9,hl1,sys_size)] = 1.0;   

                            /* MIXTURE 2 -- LIQUID AND BUBBLES */

    // Liquid Mixture Volume Equation
    lhs_jac->matrix0[M2Ac(10,V2,sys_size)] = 1.0;

    // Liquid Mixture Specific Internal Energy Equation
    lhs_jac->matrix0[M2Ac(11,U2,sys_size)] = 1.0; 

    // Liquid Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor     
    lhs_jac->matrix0[M2Ac(12,rv2,sys_size)] = 1.0;

    lhs_jac->matrix0[M2Ac(13,rl2,sys_size)] = 1.0;

    lhs_jac->matrix0[M2Ac(14,hv2,sys_size)] = 1.0; 

    lhs_jac->matrix0[M2Ac(15,hl2,sys_size)] = 1.0;  

    /////// EXTERNAL /////////
    lhs_jac->matrix0[M2Ac(16,P,sys_size)] = 1.0;  // Fixed pressure

    /////////// COMPLEMENTARITY ///////////////// 
    // LHS is set to 0 for complementarity equations    
}

void display_model_state(double* z, int n, double h, int iteration){

    std::cout << "Solution  = "<< std::endl;
    for(unsigned int i=0;i<n;i++)
    {
        std::cout << z[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "### general variables" << std::endl;
    std::cout << " V  = "<<  z[V1]+z[V2]<< " |  P = " <<  z[P]  << std::endl; 
    std::cout << "### MIXTURE 1 (vapour) variables" << std::endl;
    std::cout << " Mvd  = " <<  z[M1] << " |  uvd = " <<  z[U1]  <<  " |  Vvd = "<<  z[V1]  << std::endl;
    std::cout << " TV  = " <<  z[T1] <<  " |  xv = " <<  z[x1] << std::endl;
    std::cout << " rv  = " <<  z[rv1] << " |  rd = " <<  z[rl1] << std::endl;
    std::cout << " hv  = " <<  z[hv1] << " |  hd = " <<  z[hl1] << std::endl;
    std::cout << " DPV = " << z[DP1] << " |  P-PsatV = " << z[P]-psat(z[T1]) << std::endl;
    std::cout << "### MIXTURE 2 (liquid) variables" << std::endl;
    std::cout << " Mlb  = " <<  z[M2] << " |  ulb = " <<  z[U2]   << " |  V2 = " <<  z[V2] <<  std::endl;
    std::cout << " TL  = " <<  z[T2] <<  " |  xb = " <<  z[x2] << std::endl;
    std::cout << " rl  = " <<  z[rl2] << " |  rb = " <<  z[rv2] << std::endl;
    std::cout << " hl  = " <<  z[hl2] << " |  hb = " <<  z[hv2] << std::endl;
    std::cout << " DPL = " << z[DP2] << " |  P-PsatL = " << z[P]-psat(z[T2]) << std::endl;
    std::cout << std::endl;
    std::cout << "t =  " << ((iteration)*h) << " seconds " << std::endl;
    std::cout << "####\n" << std::endl;
}