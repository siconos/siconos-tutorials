#include "ModelIsochore.hpp"

ModelNoWallIsochore::ModelNoWallIsochore(){

    vars_name = {"P","M1","V1","U1","rv1","rl1","hv1","hl1","T1",
                 "M2","V2","U2","rv2","rl2","hv2","hl2","T2","x1",
                  "DP1","x2","DP2"};

    eqs_name = {"dM1","dU1","dM2","dU2","V1eq","u1eq","rv1eq",
                "rl1eq","hv1eq","hl1eq","V2eq","u2eq","rv2eq","rl2eq",
                "hv2eq","hl2eq","isoV_eq","x1_nseq","DP1_nseq","x2_ns","DP2_nseq"};              

    nb_vars = vars_name.size();
    nb_tot_eqs = eqs_name.size();
    nb_diff_eqs = 4;
    nb_compl_eqs = 4;
    nb_equa_eqs = nb_tot_eqs - nb_compl_eqs;


    constants.Ccond = 0.01; // condensation rate constant
    constants.Cvap = 0.09; 
    constants.Kvl = 1000; 
    constants.xqext1 =0 ; // specific (mass-ratio) heat energy input
    constants.xqext2 = 0;
    constants.qext1 = 0; // specific (mass) heat energy input
    constants.qext2 = 0;
    constants.Qext1 = 0; // constant heat energy input
    constants.Qext2 = 0;
    constants.maxTotVolume = 0; // Fixed Pressure (= P(0) in isobar case)

}

void ModelNoWallIsochore::display_model_state(double h, int iteration){
    Model::display_model_state(h,iteration);

    std::cout << "Using Model ISOCHORE WITH NO WALL \n" << std::endl;
}

int ModelNoWallIsochore::initialize_model(double P0, double M10, double x10, double T10,
                                        double M20, double x20, double T20){

                                          
    z = (double *)calloc(nb_tot_eqs,sizeof(double)); 
    /////////////////////////// INITIALISATION MIXTURE 1 (VAPOUR) //////////////////////
    z[M1] = M10;
    z[x1] = x10;
    z[T1] = T10;
    z[P] = P0;
    z[rv1] = laws->rhov(z[P],z[T1]);   // kg/m^3
    z[rl1] = laws->rhol(z[P],z[T1]);   // kg/m^3
    z[hv1] = laws->h_v(z[P],z[T1]);    // J/kg
    z[hl1] = laws->h_l(z[P],z[T1]);    // J/kg

    z[V1] = z[M1]*z[x1]/z[rv1] + z[M1]*(1-z[x1])/z[rl1];

    if (z[P]-laws->psat(z[T1]) < 0.0)
      z[DP1] = -1.0*(z[P]-laws->psat(z[T1]));  // [P-Psat]⁻ in Pa
    else   
      z[DP1] = 0.0;  // [P-Psat]⁻ in Pa

    z[U1] = (1.0-z[x1])*(z[hl1] - z[P]/z[rl1]) + z[x1]*(z[hv1] - z[P]/z[rv1]); // J
    

    /////////////////////////// INITIALISATION MIXTURE 2 (LIQUID) //////////////////////

    z[M2] = M20;
    z[x2] = x20;
    z[T2] = T20;
    z[rv2] = laws->rhov(z[P],z[T2]);   // kg/m^3
    z[rl2] = laws->rhol(z[P],z[T2]);   // kg/m^3
    z[hv2] = laws->h_v(z[P],z[T2]);    // J/kg
    z[hl2] = laws->h_l(z[P],z[T2]);    // J/kg

    z[V2] = z[M2]*z[x2]/z[rv2] + z[M2]*(1-z[x2])/z[rl2];

    constants.maxTotVolume = z[V1]+z[V2]; // Assumed correct at t=0

    if (z[P]-laws->psat(z[T2]) < 0.0)
      z[DP2] = -1.0*(z[P]-laws->psat(z[T2]));  // [P-Psat]⁻ in Pa
    else   
      z[DP2] = 0.0;  // [P-Psat]⁻ in Pa

    z[U2] = (1.0-z[x2])*(z[hl2] - z[P]/z[rl2]) + z[x2]*(z[hv2] - z[P]/z[rv2]); // J


    std::cout << "\nInitial condition, z = [";
    for(int i=0;i<nb_tot_eqs;i++){
      std::cout << z[i] << ", ";
    }
    std::cout << " ]" << std::endl;
}

void ModelNoWallIsochore::set_ExternalHeatExchange(double xqext1, double xqext2, double qext1,
                                                 double qext2, double Qext1, double Qext2){
    constants.xqext1=xqext1;
    constants.xqext2=xqext2;
    constants.qext1=qext1;
    constants.qext2=qext2;
    constants.Qext1=Qext1;
    constants.Qext2=Qext2;
}

/* Model of condensation mass transfert */
double ModelNoWallIsochore::massCond(const double &M_, const double &x_){
    return constants.Ccond*M_*(1-x_);
} 
double ModelNoWallIsochore::massConddx(const double &M_, const double &x_){
    return -1.0*constants.Ccond*M_;
} 
double ModelNoWallIsochore::massConddM(const double &M_, const double &x_){
    return constants.Ccond*(1-x_);
} 


/* Model of evaporatation mass transfert */
double ModelNoWallIsochore::massEvap(const double &M_, const double &x_){
    return constants.Cvap*M_*x_;
} 
double ModelNoWallIsochore::massEvapdx(const double &M_, const double &x_){
    return constants.Cvap*M_;
} 
double ModelNoWallIsochore::massEvapdM(const double &M_, const double &x_){
    return constants.Cvap*x_;
} 


/* Model of condensation energy transfert */
double ModelNoWallIsochore::energCond(const double &M_, const double &x_, const double &h_){
    return constants.Ccond*h_*M_*(1-x_);
}
double ModelNoWallIsochore::energConddM(const double &M_, const double &x_, const double &h_){
    return constants.Ccond*h_*(1-x_);
}
double ModelNoWallIsochore::energConddx(const double &M_, const double &x_, const double &h_){
    return -1.0*constants.Ccond*h_*M_;
}
double ModelNoWallIsochore::energConddh(const double &M_, const double &x_, const double &h_){
    return constants.Ccond*M_*(1-x_);
}

/* Model of evaportation energy transfert */
double ModelNoWallIsochore::energEvap(const double &M_, const double &x_, const double &h_){
    return constants.Cvap*h_*M_*x_;
}
double ModelNoWallIsochore::energEvapdM(const double &M_, const double &x_, const double &h_){
    return constants.Cvap*h_*x_;
}
double ModelNoWallIsochore::energEvapdx(const double &M_, const double &x_, const double &h_){
    return constants.Cvap*h_*M_;
}
double ModelNoWallIsochore::energEvapdh(const double &M_, const double &x_, const double &h_){
    return constants.Cvap*M_*x_;
}


/* Model of external energy transfert  - Constant = Qext 
   TODO remove if not used.
   TODO Add other external interaction 
*/
// double energIn(void){
//     return CTE.Qext;
// }


double ModelNoWallIsochore::energIn(const double & M_,const double & Q_){
    return M_*Q_;
}
double ModelNoWallIsochore::energIndM(const double & M_,const double & Q_){
    return Q_;
}

/* Model of external energy transfert  - Constant = Qext 
   TODO remove if not used.
   TODO Add other external interaction 
*/
double ModelNoWallIsochore::energIn1(const double & x_,const double & Q_){
    return x_*Q_;
}
double ModelNoWallIsochore::energIn1dx(const double & x_,const double & Q_){
    return Q_;
}

double ModelNoWallIsochore::energIn2(const double & x_,const double & Q_){
    return (1.0-x_)*Q_;
}
double ModelNoWallIsochore::energIn2dx(const double & x_,const double & Q_){
    return -1.0*Q_;
}

/* Model of energy conduction at mixture interface 
        FROM phase 1 TO phase 2 
>0 if energy is transfered to phase 2 (T1>T2)*/
double ModelNoWallIsochore::tempConduc_1_2(const double &T1_, const double &T2_){
    return constants.Kvl*(T1_-T2_);
}
double ModelNoWallIsochore::tempConduc_1_2dT1(const double &T1_, const double &T2_){
    return constants.Kvl;
}
double ModelNoWallIsochore::tempConduc_1_2dT2(const double &T1_, const double &T2_){
    return -1.0*constants.Kvl;
}

/* Volume equation of mixture */
double ModelNoWallIsochore::volMix(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return x_*laws->vol(rv_,M_) + (1.0-x_)*laws->vol(rl_,M_);
}
double ModelNoWallIsochore::volMixdrl(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return (1.0-x_)*laws->voldr(rl_,M_);
}
double ModelNoWallIsochore::volMixdrv(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return x_*laws->voldr(rv_,M_);
}
double ModelNoWallIsochore::volMixdM(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return x_*laws->voldM(rv_,M_) + (1.0-x_)*laws->voldM(rl_,M_);
}
double ModelNoWallIsochore::volMixdx(const double &rl_, const double &rv_, const double &M_, const double &x_){
    return laws->vol(rv_,M_) - laws->vol(rl_,M_);
}


/* specific internal energy of mixture */ 
double ModelNoWallIsochore::uMix(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return x_*laws->u(P_,hv_,rv_) + (1.0-x_)*laws->u(P_,hl_,rl_);
}
double ModelNoWallIsochore::uMixdhl(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return (1.0-x_)*laws->udh(P_,hl_,rl_);
}
double ModelNoWallIsochore::uMixdhv(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return x_*laws->udh(P_,hv_,rv_);
}
double ModelNoWallIsochore::uMixdrl(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return (1.0-x_)*laws->udr(P_,hl_,rl_);
}
double ModelNoWallIsochore::uMixdrv(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return x_*laws->udr(P_,hv_,rv_);
}
double ModelNoWallIsochore::uMixdP(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return x_*laws->udP(P_,hv_,rv_) + (1.0-x_)*laws->udP(P_,hl_,rl_);
}
double ModelNoWallIsochore::uMixdx(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_){
    return laws->u(P_,hv_,rv_) - laws->u(P_,hl_,rl_);
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
void ModelNoWallIsochore::doubleMixtureRHS(double* z_, double* rhs){

                        /* Differential Equations */

    // Vapor Mixture Mass Balance Equation
    rhs[dM1] = massEvap(z_[M2],z_[x2]) - massCond(z_[M1],z_[x1]);  
    // Vapor Mixture Energy Balance Equation
    rhs[dU1] = energEvap(z_[M2],z_[x2],z_[hv2]) - tempConduc_1_2(z_[T1],z_[T2]) - energCond(z_[M1],z_[x1],z_[hl1]) +
            constants.Qext1 + energIn(z_[M1],constants.qext1) + energIn1(z_[x1],constants.xqext1);    
    // Liquid Mixture Mass Balance Equation
    rhs[dM2] = massCond(z_[M1],z_[x1]) - massEvap(z_[M2],z_[x2]); 
    // Liquid Mixture Energy Balance Equation
    rhs[dU2] = energCond(z_[M1],z_[x1],z_[hl1]) + tempConduc_1_2(z_[T1],z_[T2]) -energEvap(z_[M2],z_[x2],z_[hv2]) +
             constants.Qext2 + energIn(z_[M2],constants.qext2) + energIn2(z_[x2],constants.xqext2);;  

                                /* MIXTURE 1 -- VAPOR AND DROPPLETS */
    // Vapor Mixture Volume Equation
    rhs[V1eq] = volMix(z_[rl1],z_[rv1],z_[M1],z_[x1]);
    // Vapor Mixture Specific Internal Energy Equation
    rhs[u1eq] = uMix(z_[hl1],z_[hv1],z_[rl1],z_[rv1],z_[P],z_[x1]);
    // Vapor Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor         
    rhs[rv1eq] = laws->rhov(z_[P],z_[T1]);
    rhs[rl1eq] = laws->rhol(z_[P],z_[T1]);
    rhs[hv1eq] = laws->h_v(z_[P],z_[T1]); 
    rhs[hl1eq] = laws->h_l(z_[P],z_[T1]);   

                            /* MIXTURE 2 -- LIQUID AND BUBBLES */
                  
    // Liquid Mixture Volume Equation
    rhs[V2eq] = volMix(z_[rl2],z_[rv2],z_[M2],z_[x2]);

    // Liquid Mixture Specific Internal Energy Equation
    rhs[u2eq] = uMix(z_[hl2],z_[hv2],z_[rl2],z_[rv2],z_[P],z_[x2]); 
    // Liquid Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor     
    rhs[rv2eq] = laws->rhov(z_[P],z_[T2]);
    rhs[rl2eq] = laws->rhol(z_[P],z_[T2]);
    rhs[hv2eq] = laws->h_v(z_[P],z_[T2]); 
    rhs[hl2eq] = laws->h_l(z_[P],z_[T2]);  

    /////// EXTERNAL /////////
    rhs[isoV_eq] = constants.maxTotVolume;      // Fixed total volume

    /////////// COMPLEMENTARITY ///////////////// 
    // Vapor mixture complementarity
    rhs[x1_nseq] = z_[P] - laws->psat(z_[T1]) + z_[DP1];
    rhs[DP1_nseq] = 1.0-z_[x1];
    // Liquid mixture complementarity
    rhs[x2_ns] = z_[P] - laws->psat(z_[T2]) + z_[DP2];
    rhs[DP2_nseq] = 1.0-z_[x2];
}

// Jacobian of the MCS right hand side
void ModelNoWallIsochore::doubleMixtureRHS_Jacobian(double* z_, NumericsMatrix* rhs_jac){

                    /* Differential Equations */

    /* Vapor Mixture Mass Balance Equation
         massEvap(z_[M2],z_[x2]) - massCond(z_[M1],z_[x1])
    */
    rhs_jac->matrix0[M2Ac(dM1,M1,nb_vars)] = -1.0*massConddM(z_[M1],z_[x1]);
    rhs_jac->matrix0[M2Ac(dM1,x1,nb_vars)] = -1.0*massConddx(z_[M1],z_[x1]);
    rhs_jac->matrix0[M2Ac(dM1,M2,nb_vars)] = massEvapdM(z_[M2],z_[x2]);
    rhs_jac->matrix0[M2Ac(dM1,x2,nb_vars)] = massEvapdx(z_[M2],z_[x2]);

    /* Vapor Mixture Energy Balance Equation
            energEvap(z_[M2],z_[x2],z_[hv2]) - tempConduc_1_2(z_[T1],z_[T2]) - energCond(z_[M1],z_[x1],z_[hl1])
    */  
    rhs_jac->matrix0[M2Ac(dU1,M1,nb_vars)] = -1.0*energConddM(z_[M1],z_[x1],z_[hl1]) + energIndM(z_[M1],constants.qext1) ; 
    rhs_jac->matrix0[M2Ac(dU1,M2,nb_vars)] = energEvapdM(z_[M2],z_[x2],z_[hv2]) ;    
    rhs_jac->matrix0[M2Ac(dU1,x1,nb_vars)] = -1.0*energConddx(z_[M1],z_[x1],z_[hl1]) + energIn1dx(z_[x1],constants.xqext1);    
    rhs_jac->matrix0[M2Ac(dU1,x2,nb_vars)] = energEvapdx(z_[M2],z_[x2],z_[hv2]) ;    
    rhs_jac->matrix0[M2Ac(dU1,hl1,nb_vars)] = -1.0*energConddh(z_[M1],z_[x1],z_[hl1]) ;    
    rhs_jac->matrix0[M2Ac(dU1,hv2,nb_vars)] = energEvapdh(z_[M2],z_[x2],z_[hv2]) ;    
    rhs_jac->matrix0[M2Ac(dU1,T1,nb_vars)] = -1.0*tempConduc_1_2dT1(z_[T1],z_[T2]);    
    rhs_jac->matrix0[M2Ac(dU1,T2,nb_vars)] = -1.0*tempConduc_1_2dT2(z_[T1],z_[T2]);    

    /* Liquid Mixture Mass Balance Equation
        massCond(z_[M1],z_[x1]) - massEvap(z_[M2],z_[x2]); 
    */
    rhs_jac->matrix0[M2Ac(dM2,M1,nb_vars)] = massConddM(z_[M1],z_[x1]);
    rhs_jac->matrix0[M2Ac(dM2,x1,nb_vars)] = massConddx(z_[M1],z_[x1]);
    rhs_jac->matrix0[M2Ac(dM2,M2,nb_vars)] = -1.0*massEvapdM(z_[M2],z_[x2]);
    rhs_jac->matrix0[M2Ac(dM2,x2,nb_vars)] = -1.0*massEvapdx(z_[M2],z_[x2]);


    /* Liquid Mixture Energy Balance Equation
        energCond(z_[M1],z_[x1],z_[hl1]) + tempConduc_1_2(z_[T1],z_[T2]) -energEvap(z_[M2],z_[x2],z_[hv2]) + CTE.Qext;     
    */
    rhs_jac->matrix0[M2Ac(dU2,M1,nb_vars)] = energConddM(z_[M1],z_[x1],z_[hl1]) ;     
    rhs_jac->matrix0[M2Ac(dU2,x1,nb_vars)] = energConddx(z_[M1],z_[x1],z_[hl1]) ;     
    rhs_jac->matrix0[M2Ac(dU2,hl1,nb_vars)] = energConddh(z_[M1],z_[x1],z_[hl1]) ;     
    rhs_jac->matrix0[M2Ac(dU2,T1,nb_vars)] = tempConduc_1_2dT1(z_[T1],z_[T2]);     
    rhs_jac->matrix0[M2Ac(dU2,T2,nb_vars)] = tempConduc_1_2dT2(z_[T1],z_[T2]);     
    rhs_jac->matrix0[M2Ac(dU2,M2,nb_vars)] = -1.0*energEvapdM(z_[M2],z_[x2],z_[hv2]) + energIndM(z_[M2],constants.qext2);     
    rhs_jac->matrix0[M2Ac(dU2,x2,nb_vars)] = -1.0*energEvapdx(z_[M2],z_[x2],z_[hv2]) + energIn2dx(z_[x2],constants.xqext2);   
    rhs_jac->matrix0[M2Ac(dU2,hv2,nb_vars)] = -1.0*energEvapdh(z_[M2],z_[x2],z_[hv2]) ;  
 
                    /* MIXTURE 1 -- VAPOR AND DROPPLETS */

    /* Vapor Mixture Volume Equation
         volMix(z_[rl1],z_[rv1],z_[M1],z_[x1]);
    */
    rhs_jac->matrix0[M2Ac(V1eq,M1,nb_vars)] =  volMixdM(z_[rl1],z_[rv1],z_[M1],z_[x1]);
    rhs_jac->matrix0[M2Ac(V1eq,x1,nb_vars)] =  volMixdx(z_[rl1],z_[rv1],z_[M1],z_[x1]);
    rhs_jac->matrix0[M2Ac(V1eq,rl1,nb_vars)] =  volMixdrl(z_[rl1],z_[rv1],z_[M1],z_[x1]);
    rhs_jac->matrix0[M2Ac(V1eq,rv1,nb_vars)] =  volMixdrv(z_[rl1],z_[rv1],z_[M1],z_[x1]);


    /* Vapor Mixture Specific Internal Energy Equation
         uMix(z_[hl1],z_[hv1],z_[rl1],z_[rv1],z_[P],z_[x1]);
    */
    rhs_jac->matrix0[M2Ac(u1eq,hl1,nb_vars)] = uMixdhl(z_[hl1],z_[hv1],z_[rl1],z_[rv1],z_[P],z_[x1]);
    rhs_jac->matrix0[M2Ac(u1eq,hv1,nb_vars)] = uMixdhv(z_[hl1],z_[hv1],z_[rl1],z_[rv1],z_[P],z_[x1]);
    rhs_jac->matrix0[M2Ac(u1eq,rl1,nb_vars)] = uMixdrl(z_[hl1],z_[hv1],z_[rl1],z_[rv1],z_[P],z_[x1]);
    rhs_jac->matrix0[M2Ac(u1eq,rv1,nb_vars)] = uMixdrv(z_[hl1],z_[hv1],z_[rl1],z_[rv1],z_[P],z_[x1]);
    rhs_jac->matrix0[M2Ac(u1eq,P,nb_vars)] = uMixdP(z_[hl1],z_[hv1],z_[rl1],z_[rv1],z_[P],z_[x1]);
    rhs_jac->matrix0[M2Ac(u1eq,x1,nb_vars)] = uMixdx(z_[hl1],z_[hv1],z_[rl1],z_[rv1],z_[P],z_[x1]);

    /* Vapor Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
     one equation for each constituant of the vapor */     
    /*  rhov(z_[P],z_[T1]);  */   
    rhs_jac->matrix0[M2Ac(rv1eq,P,nb_vars)] = laws->rhovdP(z_[P],z_[T1]);
    rhs_jac->matrix0[M2Ac(rv1eq,T1,nb_vars)] = laws->rhovdT(z_[P],z_[T1]);

    /*  rhol(z_[P],z_[T1]);  */   
    rhs_jac->matrix0[M2Ac(rl1eq,P,nb_vars)] = laws->rholdP(z_[P],z_[T1]);
    rhs_jac->matrix0[M2Ac(rl1eq,T1,nb_vars)] = laws->rholdT(z_[P],z_[T1]);

    /*  h_v(z_[P],z_[T1]);  */ 
    rhs_jac->matrix0[M2Ac(hv1eq,P,nb_vars)] = laws->h_vdP(z_[P],z_[T1]); 
    rhs_jac->matrix0[M2Ac(hv1eq,T1,nb_vars)] = laws->h_vdT(z_[P],z_[T1]); 

    /*  h_l(z_[P],z_[T1]);  */
    rhs_jac->matrix0[M2Ac(hl1eq,P,nb_vars)] = laws->h_ldP(z_[P],z_[T1]); 
    rhs_jac->matrix0[M2Ac(hl1eq,T1,nb_vars)] = laws->h_ldT(z_[P],z_[T1]);  

                /* MIXTURE 2 -- LIQUID AND BUBBLES */
    
        /* Liquid Mixture Volume Equation
        volMix(z_[rl2],z_[rv2],z_[M2],z_[x2]);
    */  
                    
    rhs_jac->matrix0[M2Ac(V2eq,rl2,nb_vars)] = volMixdrl(z_[rl2],z_[rv2],z_[M2],z_[x2]);
    rhs_jac->matrix0[M2Ac(V2eq,rv2,nb_vars)] = volMixdrv(z_[rl2],z_[rv2],z_[M2],z_[x2]);
    rhs_jac->matrix0[M2Ac(V2eq,M2,nb_vars)] = volMixdM(z_[rl2],z_[rv2],z_[M2],z_[x2]);
    rhs_jac->matrix0[M2Ac(V2eq,x2,nb_vars)] = volMixdx(z_[rl2],z_[rv2],z_[M2],z_[x2]);


    /* Liquid Mixture Specific Internal Energy Equation
        uMix(z_[hl2],z_[hv2],z_[rl2],z_[rv2],z_[P],z_[x2]); 
    */
    rhs_jac->matrix0[M2Ac(u2eq,hl2,nb_vars)] = uMixdhl(z_[hl2],z_[hv2],z_[rl2],z_[rv2],z_[P],z_[x2]); 
    rhs_jac->matrix0[M2Ac(u2eq,hv2,nb_vars)] = uMixdhv(z_[hl2],z_[hv2],z_[rl2],z_[rv2],z_[P],z_[x2]); 
    rhs_jac->matrix0[M2Ac(u2eq,rl2,nb_vars)] = uMixdrl(z_[hl2],z_[hv2],z_[rl2],z_[rv2],z_[P],z_[x2]); 
    rhs_jac->matrix0[M2Ac(u2eq,rv2,nb_vars)] = uMixdrv(z_[hl2],z_[hv2],z_[rl2],z_[rv2],z_[P],z_[x2]); 
    rhs_jac->matrix0[M2Ac(u2eq,P,nb_vars)] = uMixdP(z_[hl2],z_[hv2],z_[rl2],z_[rv2],z_[P],z_[x2]); 
    rhs_jac->matrix0[M2Ac(u2eq,x2,nb_vars)] = uMixdx(z_[hl2],z_[hv2],z_[rl2],z_[rv2],z_[P],z_[x2]);
    
    
    /* Liquid Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
     one equation for each constituant of the vapor     */
    /*  rhov(z_[P],z_[T2]); */
    rhs_jac->matrix0[M2Ac(rv2eq,P,nb_vars)] = laws->rhovdP(z_[P],z_[T2]);
    rhs_jac->matrix0[M2Ac(rv2eq,T2,nb_vars)] = laws->rhovdT(z_[P],z_[T2]);

    /*  rhol(z_[P],z_[T2]); */
    rhs_jac->matrix0[M2Ac(rl2eq,P,nb_vars)] = laws->rholdP(z_[P],z_[T2]);
    rhs_jac->matrix0[M2Ac(rl2eq,T2,nb_vars)] = laws->rholdT(z_[P],z_[T2]);

    /*  h_v(z_[P],z_[T2]); */ 
    rhs_jac->matrix0[M2Ac(hv2eq,P,nb_vars)] = laws->h_vdP(z_[P],z_[T2]);
    rhs_jac->matrix0[M2Ac(hv2eq,T2,nb_vars)] = laws->h_vdT(z_[P],z_[T2]);

    /*  h_l(z_[P],z_[T2]);  */
    rhs_jac->matrix0[M2Ac(hl2eq,P,nb_vars)] = laws->h_ldP(z_[P],z_[T2]);
    rhs_jac->matrix0[M2Ac(hl2eq,T2,nb_vars)] = laws->h_ldT(z_[P],z_[T2]);

    /////// EXTERNAL /////////
    /* Cte --> line [16] of rhs_jac full of 0    */

    /////////// COMPLEMENTARITY ///////////////// 
    /* (Vapor) mixture 1 complementarity */
    /* z_[P] - psat(z_[T1]) + z_[DP1]; */
    rhs_jac->matrix0[M2Ac(x1_nseq,P,nb_vars)] = 1.0;
    rhs_jac->matrix0[M2Ac(x1_nseq,T1,nb_vars)] = -1.0*laws->psatdT(z_[T1]);
    rhs_jac->matrix0[M2Ac(x1_nseq,DP1,nb_vars)] = 1.0;

    /* 1.0-z_[x1];*/
    rhs_jac->matrix0[M2Ac(DP1_nseq,x1,nb_vars)] = -1.0;

    /* (Liquid) mixture 2 complementarity */
    /* z_[P] - psat(z_[T2]) + z_[DP2];*/
    rhs_jac->matrix0[M2Ac(x2_ns,P,nb_vars)] = 1.0;
    rhs_jac->matrix0[M2Ac(x2_ns,T2,nb_vars)] = -1.0*laws->psatdT(z_[T2]);
    rhs_jac->matrix0[M2Ac(x2_ns,DP2,nb_vars)] = 1.0;

    /* 1.0-z_[x2];*/
    rhs_jac->matrix0[M2Ac(DP2_nseq,x2,nb_vars)] = -1.0;
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
void ModelNoWallIsochore::doubleMixtureLHS(void* env, FiniteDifference* diffMethod, double* z_, double* lhs){
    
                      /* Differential Equations */
   
    // Vapor Mixture Mass Balance Equation
    lhs[dM1] = diffMethod->finiteDifference(env,z_[M1],z[M1]);  
    // Vapor Mixture Energy Balance Equation
    lhs[dU1] = z_[M1]*diffMethod->finiteDifference(env,z_[U1],z[U1]) + z_[U1]*diffMethod->finiteDifference(env,z_[M1],z[M1]); 
    // Liquid Mixture Mass Balance Equation
    lhs[dM2] = diffMethod->finiteDifference(env,z_[M2],z[M2]);
    // Liquid Mixture Energy Balance Equation
    lhs[dU2] = z_[M2]*diffMethod->finiteDifference(env,z_[U2],z[U2]) + z_[U2]*diffMethod->finiteDifference(env,z_[M2],z[M2]);       

                        /* MIXTURE 1 -- VAPOR AND DROPPLETS */

    // Vapor Mixture Volume Equation
    lhs[V1eq] = z_[V1];
    // Vapor Mixture Specific Internal Energy Equation
    lhs[u1eq] = z_[U1]; 
    // Vapor Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor         
    lhs[rv1eq] = z_[rv1];
    lhs[rl1eq] = z_[rl1];
    lhs[hv1eq] = z_[hv1]; 
    lhs[hl1eq] = z_[hl1];   

                            /* MIXTURE 2 -- LIQUID AND BUBBLES */
     
    // Liquid Mixture Volume Equation
    lhs[V2eq] = z_[V2];

    // Liquid Mixture Specific Internal Energy Equation
    lhs[u2eq] = z_[U2]; 
    // Liquid Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor     
    lhs[rv2eq] = z_[rv2];
    lhs[rl2eq] = z_[rl2];
    lhs[hv2eq] = z_[hv2]; 
    lhs[hl2eq] = z_[hl2];  

    /////// EXTERNAL /////////
    lhs[isoV_eq] = z_[V1] + z_[V2];  // Fixed total Volume 

    /////////// COMPLEMENTARITY ///////////////// 
    // Vapor mixture complementarity
    lhs[x1_nseq] = 0;
    lhs[DP1_nseq] = 0;
    // Liquid mixture complementarity
    lhs[x2_ns] = 0;
    lhs[DP2_nseq] = 0;
}

/* Jacobian of the MCS Left Hand side
    ptrDifferentiationFunctionDz is a pointer to a function computing the partial derivative of a differentiated variable
*/
void ModelNoWallIsochore::doubleMixtureLHS_Jacobian(void* env, FiniteDifference* diffMethod, double* z_, NumericsMatrix* lhs_jac){

                            /* Differential Equations */

    // Vapor Mixture Mass Balance Equation
    lhs_jac->matrix0[M2Ac(dM1,M1,nb_vars)] = diffMethod->finiteDifferenceDz(env,z_[M1],z[M1]);

    // Vapor Mixture Energy Balance Equation
    lhs_jac->matrix0[M2Ac(dU1,M1,nb_vars)] = diffMethod->finiteDifference(env,z_[U1],z[U1]) + z_[U1]*diffMethod->finiteDifferenceDz(env,z_[M1],z[M1]); 
    lhs_jac->matrix0[M2Ac(dU1,U1,nb_vars)] = z_[M1]*diffMethod->finiteDifferenceDz(env,z_[U1],z[U1]) + diffMethod->finiteDifference(env,z_[M1],z[M1]);

    // Liquid Mixture Mass Balance Equation
    lhs_jac->matrix0[M2Ac(dM2,M2,nb_vars)] = diffMethod->finiteDifferenceDz(env,z_[M2],z[M2]);

    // Liquid Mixture Energy Balance Equation
    lhs_jac->matrix0[M2Ac(dU2,M2,nb_vars)] = diffMethod->finiteDifference(env,z_[U2],z[U2]) + z_[U2]*diffMethod->finiteDifferenceDz(env,z_[M2],z[M2]); 
    lhs_jac->matrix0[M2Ac(dU2,U2,nb_vars)] = z_[M2]*diffMethod->finiteDifferenceDz(env,z_[U2],z[U2]) + diffMethod->finiteDifference(env,z_[M2],z[M2]); 

                        /* MIXTURE 1 -- VAPOR AND DROPPLETS */

    // Vapor Mixture Volume Equation
    lhs_jac->matrix0[M2Ac(V1eq,V1,nb_vars)] = 1.0;

    // Vapor Mixture Specific Internal Energy Equation
    lhs_jac->matrix0[M2Ac(u1eq,U1,nb_vars)] = 1.0; 

    // Vapor Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor         
    lhs_jac->matrix0[M2Ac(rv1eq,rv1,nb_vars)] = 1.0;

    lhs_jac->matrix0[M2Ac(rl1eq,rl1,nb_vars)] = 1.0;

    lhs_jac->matrix0[M2Ac(hv1eq,hv1,nb_vars)] = 1.0; 

    lhs_jac->matrix0[M2Ac(hl1eq,hl1,nb_vars)] = 1.0;   

                            /* MIXTURE 2 -- LIQUID AND BUBBLES */

    // Liquid Mixture Volume Equation
    lhs_jac->matrix0[M2Ac(V2eq,V2,nb_vars)] = 1.0;

    // Liquid Mixture Specific Internal Energy Equation
    lhs_jac->matrix0[M2Ac(u2eq,U2,nb_vars)] = 1.0; 

    // Liquid Mixture State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor     
    lhs_jac->matrix0[M2Ac(rv2eq,rv2,nb_vars)] = 1.0;

    lhs_jac->matrix0[M2Ac(rl2eq,rl2,nb_vars)] = 1.0;

    lhs_jac->matrix0[M2Ac(hv2eq,hv2,nb_vars)] = 1.0; 

    lhs_jac->matrix0[M2Ac(hl2eq,hl2,nb_vars)] = 1.0;  

    /////// EXTERNAL /////////

    lhs_jac->matrix0[M2Ac(isoV_eq,V1,nb_vars)] = 1.0;  // Fixed Volume
    lhs_jac->matrix0[M2Ac(isoV_eq,V2,nb_vars)] = 1.0;  

    /////////// COMPLEMENTARITY ///////////////// 
    // LHS is set to 0 for complementarity equations    
}