#ifndef MODEL_HPP
#define MODEL_HPP

#include "model_common.hpp"
#include "laws.hpp"

// phase 1: Vapor + dropplets
// phase 2: Liquid + bubbles
// TODO need to be replaced by a map<string,int>
enum Variables {P,M1,V1,U1,rv1,rl1,hv1,hl1,T1,M2,V2,U2,rv2,rl2,hv2,hl2,T2,x1,DP1,x2,DP2};


/* Differentiation function Type
    - env : contains any information needed for the numerical differentiation (previous step, etc ...)
    - z_i : contains the latest value variable differentiated 
    - i : is the index of this variable
*/
typedef double (*ptrDifferentiationFunction)(void* env, double z_i, int i);

/*partial derivative wrt z[j] of numerical differention of the variable z[i]
  TODO  typedef double (*ptrDifferentiationFunctionDz)(void* env, double z_i, int i, int j);
  TODO  in current version i=j (finite difference)
*/
typedef double (*ptrDifferentiationFunctionDz)(void* env, double z_i, int i);

/*Constants of the model*/ 
typedef struct
{
    double Ccond; // condensation rate constant
    double Cvap; 
    double Kvl; 
    double Qext;
    double Pext;

}ModelConstants;

/* Model of condensation mass transfert */
double massCond(const double &M_, const double &x_); 
double massEvapdx(const double &M_, const double &x_);
double massEvapdM(const double &M_, const double &x_);


/* Model of evaportation mass transfert */
double massEvap(const double &M_, const double &x_); 
double massEvapdx(const double &M_, const double &x_);
double massEvapdM(const double &M_, const double &x_);


/* Model of condensation energy transfert */
double energCond(const double &M_, const double &x_, const double &h);
double energConddM(const double &M_, const double &x_, const double &h);
double energConddx(const double &M_, const double &x_, const double &h);
double energConddh(const double &M_, const double &x_, const double &h);

/* Model of evaportation energy transfert */
double energEvap(const double &M_, const double &x_, const double &h);
double energEvapdM(const double &M_, const double &x_, const double &h);
double energEvapdx(const double &M_, const double &x_, const double &h);
double energEvapdh(const double &M_, const double &x_, const double &h);


/* Model of external energy transfert  - Constant */
double energInCst(const double &Qext_);

/* Model of energy conduction at mixture interface 
FROM phase 1 TO phase 2 */
double tempConduc_1_2(const double &T1_, const double &T2_);
double tempConduc_1_2dT1(const double &T1_, const double &T2_);
double tempConduc_1_2dT2(const double &T1_, const double &T2_);

/* Volume equation of mixture */
double volMix(const double &rl_, const double &rv_, const double &M_, const double &x_);
double volMixdrl(const double &rl_, const double &rv_, const double &M_, const double &x_);
double volMixdrv(const double &rl_, const double &rv_, const double &M_, const double &x_);
double volMixdM(const double &rl_, const double &rv_, const double &M_, const double &x_);
double volMixdx(const double &rl_, const double &rv_, const double &M_, const double &x_);



/* specific internal energy of mixture */ 
double uMix(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_);
double uMixdhl(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_);
double uMixdhv(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_);
double uMixdrl(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_);
double uMixdrv(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_);
double uMixdP(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_);
double uMixdx(const double &hl_, const double &hv_, const double &rl_, const double &rv_, const double &P_, const double &x_);

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
void doubleMixtureRHS(double* z, double* rhs);

// Jacobian of the MCS right hand side
void doubleMixtureRHS_Jacobian(double* z, int sys_size, NumericsMatrix* rhs_jac);

/*Complete double mixture mixed complementarity system (MCS) Left-Hand Side 
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
void doubleMixtureLHS(void* env, ptrDifferentiationFunction diff, double* z, double* lhs);

void doubleMixtureLHS_Jacobian(void* env, int sys_size, ptrDifferentiationFunction diff, ptrDifferentiationFunctionDz diff_dz, double* z, NumericsMatrix* lhs_jac);

void display_model_state(double* z, int n, double h, int iteration);

#endif // MODEL_HPP