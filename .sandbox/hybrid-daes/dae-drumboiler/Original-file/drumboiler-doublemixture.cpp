/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <boost/timer/timer.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "SiconosKernel.hpp"
#include "MCP_cst.h"
#include "NonSmoothDrivers.h"
#include "MixedComplementarityProblem.h"
#include "SolverOptions.h"
#include "MCP_Solvers.h"
#include "NumericsVerbose.h"
#include "NumericsMatrix.h"

#include "cond.h"

using namespace std;

/* example of enum renaming usage:
M(t) = z[M], H(t) = z[H], P(t) = z[P]
names: M,H,rv,rl,hv,hl,P,T,l1,l2 are taken */
enum Variables {P,M1,V1,U1,rv1,rl1,hv1,hl1,T1,M2,V2,U2,rv2,rl2,hv2,hl2,T2,x1,DP1,x2,DP2};

/* parameters value for n-butane
    specific volume and enthalpy are molar, not in mass
*/
static double R = 461.527; // perfect gas constant J/kg/K
// static double Vtot = 0.00100306;  // Total volume (m^3) // TODO not used in this version
static double Wext = 100;   //  External energy input J*s⁻1  ; // TODO Used in this version only in Liquid water term
static double Pext = 106074.0; // (bar) Constant pressure = 1atm 
static double Cvap = 0.09;
static double Ccond = 0.1;
// static double Text1 = 600; // External temperature // TODO not used in this version
static double Kvl = 10;

/* Simulations variables */
static double t0 = 0.0; // starting time (seconds)
static double Tf = 200.00;  // Final time (seconds);
static double h = 0.01;  



/*scaling vector pointer*/
static double* scaling;

/* 2d matrix conversion in 1d array 
(Column first, indexes starting at 0) */
int M2Ac(int i, int j, int nb_line)
{
    int index = nb_line*j + i;
    return index;
}


static double * z_prev;
static double * z_curr;


//* Thermodynamics state functions and their derivatives *//


double rhov(const double &P_, const double &T_)
{
    return P_/(R*T_); 
}

double rhovdP(const double &P_, const double &T_)
{
    return 1/(R*T_);
}

double rhovdT(const double &P_, const double &T_)
{
    return -P_/(R*pow(T_,2.0));
}

/* Volumic mass of water --  */
double rhol(const double &P_, const double &T_)
{
    return (-0.0025*(T_-273.156) - 0.1992)*(T_-273.156) + 1004.4;
}

double rholdP(const double &P_, const double &T_)
{
    return 0;
}

double rholdT(const double &P_, const double &T_)
{
    return (-2*0.0025*(T_-273.156) - 0.1992);
}

/* Specific enthalpy of vapor (J/kg) in fonction of Pressure (Pa) and Temperature (K) */
double h_v(const double &P_, const double &T_)
{

    double c00 = 1.778741e+06;
    double c10 = -6.997339e-02;
    double c01 = 2.423675e+03;
    double c20 = 1.958603e-10;
    double c11 = 8.100784e-05;
    double c02 = -3.747139e-01;
    double c30 = -1.016123e-19;
    double c21 = -1.234548e-13;
    double c12 = -2.324528e-08;
    double c03 = 1.891004e-04;    
    return c00  + c10*P_ + c01*T_ 
                + c20*pow(P_,2.0) + c11*P_*T_ 
                + c02*pow(T_,2.0) + c30*pow(P_,3.0) 
                + c21*pow(P_,2.0)*T_ + c12*P_*pow(T_,2.0) 
                + c03*pow(T_,3.0); 
}

double h_vdP(const double &P_, const double &T_)
{
    double c00 = 1.778741e+06;
    double c10 = -6.997339e-02;
    double c01 = 2.423675e+03;
    double c20 = 1.958603e-10;
    double c11 = 8.100784e-05;
    double c02 = -3.747139e-01;
    double c30 = -1.016123e-19;
    double c21 = -1.234548e-13;
    double c12 = -2.324528e-08;
    double c03 = 1.891004e-04;    

    return c10 + c20*2.0*P_ + c11*T_ 
               + c30*3.0*pow(P_,2.0) + c21*2.0*P_*T_ 
               + c12*pow(T_,2.0); 
}

double h_vdT(const double &P_, const double &T_)
{
    double c00 = 1.778741e+06;
    double c10 = -6.997339e-02;
    double c01 = 2.423675e+03;
    double c20 = 1.958603e-10;
    double c11 = 8.100784e-05;
    double c02 = -3.747139e-01;
    double c30 = -1.016123e-19;
    double c21 = -1.234548e-13;
    double c12 = -2.324528e-08;
    double c03 = 1.891004e-04;    

    return c01 + c11*P_ + c02*2.0*T_
               + c21*pow(P_,2.0) + c12*P_*2.0*T_ 
               + c03*3.0*pow(T_,2.0);  
}

/* Specific enthalpy of liquid (J/kg) in fonction of Pressure (Pa) and Temperature (K) */
double h_l(const double &P_, const double &T_)
{
    double c00 = -1.210342e+07;
    double c10 = 3.932901e-02;
    double c01 = 1.310565e+05;
    double c20 = -3.425284e-10;
    double c11 = -2.572281e-04;
    double c02 = -5.801243e+02;
    double c30 = 1.974339e-19;
    double c21 = 2.427381e-12;
    double c12 = 4.966543e-07;
    double c03 = 1.314839e+00;
    double c40 = -4.256626e-27;
    double c31 = 1.512868e-21;
    double c22 = -6.054694e-15;
    double c13 = -8.389491e-11;
    double c04 = -1.484055e-03;
    double c50 = 1.597043e-35;
    double c41 = 1.356624e-31;
    double c32 = -2.492294e-24;
    double c23 = 5.082575e-18;
    double c14 = -3.822957e-13;
    double c05 = 6.712484e-07;

    return c00 + c10*P_ + c01*T_ + c20*pow(P_,2.0) 
               + c11*P_*T_ + c02*pow(T_,2.0) + c30*pow(P_,3.0) 
               + c21*pow(P_,2.0)*T_ + c12*P_*pow(T_,2.0) 
               + c03*pow(T_,3.0) + c40*pow(P_,4.0) + c31*pow(P_,3.0)*T_ 
               + c22*pow(P_,2.0)*pow(T_,2.0) + c13*P_*pow(T_,3.0) 
               + c04*pow(T_,4.0) + c50*pow(P_,5.0) + c41*pow(P_,4.0)*T_ 
               + c32*pow(P_,3.0)*pow(T_,2.0) + c23*pow(P_,2.0)*pow(T_,3.0) 
               + c14*P_*pow(T_,4.0) + c05*pow(T_,5.0);
 
}


double h_ldP(const double &P_, const double &T_)
{
    double c00 = -1.210342e+07;
    double c10 = 3.932901e-02;
    double c01 = 1.310565e+05;
    double c20 = -3.425284e-10;
    double c11 = -2.572281e-04;
    double c02 = -5.801243e+02;
    double c30 = 1.974339e-19;
    double c21 = 2.427381e-12;
    double c12 = 4.966543e-07;
    double c03 = 1.314839e+00;
    double c40 = -4.256626e-27;
    double c31 = 1.512868e-21;
    double c22 = -6.054694e-15;
    double c13 = -8.389491e-11;
    double c04 = -1.484055e-03;
    double c50 = 1.597043e-35;
    double c41 = 1.356624e-31;
    double c32 = -2.492294e-24;
    double c23 = 5.082575e-18;
    double c14 = -3.822957e-13;
    double c05 = 6.712484e-07;

    return c10 + c20*2.0*P_ + c11*T_ 
               + c30*3.0*pow(P_,2.0) 
               + c21*2.0*P_*T_ + c12*pow(T_,2.0) 
               + c40*4.0*pow(P_,3.0) + c31*3.0*pow(P_,2.0)*T_ 
               + c22*2.0*P_*pow(T_,2.0) + c13*pow(T_,3.0) 
               + c50*5.0*pow(P_,4.0) + c41*4.0*pow(P_,3.0)*T_ 
               + c32*3.0*pow(P_,2.0)*pow(T_,2.0) + c23*2.0*P_*pow(T_,3.0) 
               + c14*pow(T_,4.0); 
}

double h_ldT(const double &P_, const double &T_)
{
    double c00 = -1.210342e+07;
    double c10 = 3.932901e-02;
    double c01 = 1.310565e+05;
    double c20 = -3.425284e-10;
    double c11 = -2.572281e-04;
    double c02 = -5.801243e+02;
    double c30 = 1.974339e-19;
    double c21 = 2.427381e-12;
    double c12 = 4.966543e-07;
    double c03 = 1.314839e+00;
    double c40 = -4.256626e-27;
    double c31 = 1.512868e-21;
    double c22 = -6.054694e-15;
    double c13 = -8.389491e-11;
    double c04 = -1.484055e-03;
    double c50 = 1.597043e-35;
    double c41 = 1.356624e-31;
    double c32 = -2.492294e-24;
    double c23 = 5.082575e-18;
    double c14 = -3.822957e-13;
    double c05 = 6.712484e-07;

    return c01 + c11*P_ + c02*2.0*T_ + c21*pow(P_,2.0) 
               + c12*P_*2.0*T_ + c03*3.0*pow(T_,2.0) 
               + c31*pow(P_,3.0) + c22*pow(P_,2.0)*2.0*T_ 
               + c13*P_*3.0*pow(T_,2.0) + c04*4.0*pow(T_,3.0) 
               + c41*pow(P_,4.0) + c32*pow(P_,3.0)*2.0*T_ 
               + c23*pow(P_,2.0)*3.0*pow(T_,2.0) + c14*P_*4.0*pow(T_,3.0) 
               + c05*5.0*pow(T_,4.0); 
}

/* Saturation pressure (Pa) in fonction of Temperature (K)*/
double psat(const double &T_)
{
    double c7 = 7.677448367806697e-11;
    double c6 = -2.327974895639335e-07;
    double c5 = 2.984399245658974e-04;
    double c4 = -2.081210501212062e-01;
    double c3 = 8.527291155079814e+01;
    double c2 = -2.055993982138471e+04;
    double c1 = 2.704822454027627e+06;
    double c0 = -1.499284498173245e+08;

    return c0 + c1*T_ + c2*pow(T_,2.0) + c3*pow(T_,3.0) 
              + c4*pow(T_,4.0) + c5*pow(T_,5.0) 
              + c6*pow(T_,6.0) + c7*pow(T_,7.0);

}

double psatdT(const double &T_)
{ 
    double c7=7.677448367806697e-11;
    double c6=-2.327974895639335e-07;
    double c5=2.984399245658974e-04;
    double c4=-2.081210501212062e-01;
    double c3=8.527291155079814e+01;
    double c2=-2.055993982138471e+04;
    double c1=2.704822454027627e+06;
    double c0=-1.499284498173245e+08;

    return  c1 + 2.0*c2*T_ + 3.0*c3*pow(T_,2.0) 
            + 4.0*c4*pow(T_,3.0) + 5.0*c5*pow(T_,4.0) 
            + 6.0*c6*pow(T_,5.0) + 7.0*c7*pow(T_,6.0);
}

void F(void* env, int size, double *z, double * F);
/*
void F(void* env, int size, double *z, double * F){

    // Mixture 1 (vapor)
    F[0] = z[M1]-z_prev[M1] + h*Ccond*z[M1]*(1-z[x1]);
    F[1] = z[M1]*(z[U1]-z_prev[U1]) + z[U1]*(z[M1]-z_prev[M1]) 
            - h*10*(z[T2]-z[T1]) + h*z[hv1]*Ccond*z[M1]*(1-z[x1]);// -h*Wext;                         // "Implicit" Euler discr.
    F[2] = -z[V1] + z[M1]/z[rl1] + z[M1]*z[x1]*(1/z[rv1]-1/z[rl1]) ;
    F[3] = -z[U1] + 1.0*( z[hl1] - z[P]/z[rl1] ) 
                 + z[x1]*( (z[hv1]-z[hl1]) + (z[P]/z[rl1] - z[P]/z[rv1]) ); // Rq: T is not directly in this eq.
    F[4] = -z[rv1] + rhov(z[P],z[T1]);
    F[5] = -z[rl1] + rhol(z[P],z[T1]);
    F[6] = -z[hv1] + h_v(z[P],z[T1]); 
    F[7] = -z[hl1] + h_l(z[P],z[T1]);   

    // Mixture 2 (liquide)
    F[8] = z[M2]-z_prev[M2] + h*Cvap*z[M2]*z[x2];
    F[9] = z[M2]*(z[U2]-z_prev[U2]) + z[U2]*(z[M2]-z_prev[M2]) 
            - h*10*(z[T1]-z[T2]) + h*z[hv2]*Cvap*z[M2]*z[x2];// -h*Wext;                         // "Implicit" Euler discr.
    F[10] = -z[V2] + z[M2]/z[rl2] + z[M2]*z[x2]*(1/z[rv2]-1/z[rl2]) ;
    F[11] = -z[U2] + 1.0*( z[hl2] - z[P]/z[rl2] ) 
                 + z[x2]*( (z[hv2]-z[hl2]) + (z[P]/z[rl2] - z[P]/z[rv2]) ); // Rq: T is not directly in this eq.
    F[12] = -z[rv2] + rhov(z[P],z[T2]);
    F[13] = -z[rl2] + rhol(z[P],z[T2]);
    F[14] = -z[hv2] + h_v(z[P],z[T2]); 
    F[15] = -z[hl2] + h_l(z[P],z[T2]);  

    /////// EXTERNAL /////////
    F[16] = z[P] - Pext;                  // Fixed pressure

    /////////// COMPLEMENTARITY ///////////////// 
    F[17] = z[P] - psat(z[T1]) + z[DP1];
    F[18] = 1.0-z[x1];
    F[19] = z[P] - psat(z[T2]) + z[DP2];
    F[20] = 1.0-z[x2];
}
*/
void F(void* env, int size, double *z, double * F){
    // Numerical discretization by Backward Euler
    // Mixture 1 (vapor)
    // Vapor Mass Balance Equation
    F[0] = z[M1]-z_prev[M1] + h*Ccond*z[M1]*(1-z[x1]) - h*Cvap*z[M2]*z[x2];  
    // Vapor Energy Balance Equation
    F[1] = z[M1]*(z[U1]-z_prev[U1]) + z[U1]*(z[M1]-z_prev[M1])            // -h*Wext; removed --> no external energy 
            - h*Kvl*(z[T2]-z[T1]) + h*z[hl1]*Ccond*z[M1]*(1-z[x1])
            - h*z[hv2]*Cvap*z[M2]*z[x2] ;             
    // Vapor Volume Equation
    F[2] = -z[V1] + z[M1]/z[rl1] + z[M1]*z[x1]*(1/z[rv1]-1/z[rl1]) ;
    // Vapor Specific Internal Energy Equation
    F[3] = -z[U1] + 1.0*( z[hl1] - z[P]/z[rl1] ) 
                 + z[x1]*(   (z[hv1]-z[hl1])                // Rq: T is not directly in this eq.
                           + (z[P]/z[rl1] - z[P]/z[rv1]) 
                            ); 
    // Vapor State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor         
    F[4] = -z[rv1] + rhov(z[P],z[T1]);
    F[5] = -z[rl1] + rhol(z[P],z[T1]);
    F[6] = -z[hv1] + h_v(z[P],z[T1]); 
    F[7] = -z[hl1] + h_l(z[P],z[T1]);   

    // Mixture 2 (liquide)
    // Liquid Mass Balance Equation
    F[8] = z[M2]-z_prev[M2] + h*Cvap*z[M2]*z[x2] - h*Ccond*z[M1]*(1-z[x1]); 
    // Liquid Energy Balance Equation
    F[9] = z[M2]*(z[U2]-z_prev[U2]) + z[U2]*(z[M2]-z_prev[M2])         // -h*Wext;
            - h*Kvl*(z[T1]-z[T2]) + h*z[hv2]*Cvap*z[M2]*z[x2]
            - h*z[hl1]*Ccond*z[M1]*(1-z[x1]) - h*Wext;                 
    // Liquid Volume Equation
    F[10] = -z[V2] + z[M2]/z[rl2] + z[M2]*z[x2]*(1/z[rv2]-1/z[rl2]) ;

    // Liquid Specific Internal Energy Equation
    F[11] = -z[U2] + 1.0*( z[hl2] - z[P]/z[rl2] ) 
                 + z[x2]*( (z[hv2]-z[hl2])          // Rq: T is not directly in this eq.
                 + (z[P]/z[rl2] - z[P]/z[rv2]) ); 
    // Liquid State variables equations: volumique mass (rho) and specific enthalpie (h)
    // one equation for each constituant of the vapor     
    F[12] = -z[rv2] + rhov(z[P],z[T2]);
    F[13] = -z[rl2] + rhol(z[P],z[T2]);
    F[14] = -z[hv2] + h_v(z[P],z[T2]); 
    F[15] = -z[hl2] + h_l(z[P],z[T2]);  

    /////// EXTERNAL /////////
    F[16] = z[P] - Pext;                  // Fixed pressure

    /////////// COMPLEMENTARITY ///////////////// 
    // Vapor mixture complementarity
    F[17] = z[P] - psat(z[T1]) + z[DP1];
    F[18] = 1.0-z[x1];
    // Liquid mixture complementarity
    F[19] = z[P] - psat(z[T2]) + z[DP2];
    F[20] = 1.0-z[x2];
}

void NablaF(void * env, int size, double *z, NumericsMatrix* nablaF);
/*
void NablaF(void * env, int size, double *z, NumericsMatrix* nablaF){

   /////////////////////   MIXTURE 1  ///////////////////////
   // F[0]: z[M]-z_prev[M] 
   nablaF->matrix0[M2Ac(0,M1,size)] = 1 + h*Ccond*(1-z[x1]);
   nablaF->matrix0[M2Ac(0,x1,size)] = -h*Ccond*z[M1] ;

   // F[1]: z[H]-z_prev[H] - h*Wext; 
   nablaF->matrix0[M2Ac(1,M1,size)] = (z[U1]-z_prev[U1]) + z[U1] + h*z[hv1]*Ccond*(1-z[x1]) ;
   nablaF->matrix0[M2Ac(1,U1,size)] = (z[M1]-z_prev[M1]) + z[M1];
   nablaF->matrix0[M2Ac(1,hv1,size)] = h*Ccond*z[M1]*(1-z[x1]);
   nablaF->matrix0[M2Ac(1,x1,size)] = -h*z[hv1]*Ccond*z[M1];
   nablaF->matrix0[M2Ac(1,T1,size)] = h*10;
   nablaF->matrix0[M2Ac(1,T2,size)] = -h*10;

   // F[2] : -z[V] + z[M]/z[rl] + z[l1]*(1/z[rv]-1/z[rl]) ;
   nablaF->matrix0[M2Ac(2,V1,size)] =  -1.0;
   nablaF->matrix0[M2Ac(2,M1,size)] =  1.0/z[rl1] + z[x1]*(1.0/z[rv1]-1.0/z[rl1]);
   nablaF->matrix0[M2Ac(2,rv1,size)] = - z[M1]*z[x1]/pow(z[rv1],2.0);
   nablaF->matrix0[M2Ac(2,rl1,size)] = ( z[M1]*z[x1]-z[M1])/pow(z[rl1],2.0);
   nablaF->matrix0[M2Ac(2,x1,size)] = z[M1]/z[rv1]-z[M1]/z[rl1]; 
    
   // F[3]: -z[H] + z[M]*( z[hl] - z[P]/z[rl] ) 
   //             + z[l1]*( (z[hv]-z[hl])+(z[P]/z[rl] - z[P]/z[rv]) )
   nablaF->matrix0[M2Ac(3,U1,size)] = -1; 
   nablaF->matrix0[M2Ac(3,rv1,size)] = z[x1]*z[P]/pow(z[rv1],2.0) ;
   nablaF->matrix0[M2Ac(3,rl1,size)] = z[P]/pow(z[rl1],2.0)*(1.0-z[x1]); 
   nablaF->matrix0[M2Ac(3,hv1,size)] = z[x1];
   nablaF->matrix0[M2Ac(3,hl1,size)] = 1.0-z[x1];
   nablaF->matrix0[M2Ac(3,P,size)] = -1.0/z[rl1] + z[x1]*(1/z[rl1] - 1/z[rv1]);
   nablaF->matrix0[M2Ac(3,x1,size)] = (z[hv1]-z[hl1])+(z[P]/z[rl1] - z[P]/z[rv1]);
    
   // F[4]:  -z[rv] + rhov(z)
   nablaF->matrix0[M2Ac(4,rv1,size)] =  -1;
   nablaF->matrix0[M2Ac(4,P,size)] = rhovdP(z[P],z[T1]);
   nablaF->matrix0[M2Ac(4,T1,size)] = rhovdT(z[P],z[T1]);   

   // F[5]:-z[rl] + rhol(z);
   nablaF->matrix0[M2Ac(5,rl1,size)] = -1; 
   nablaF->matrix0[M2Ac(5,P,size)] = rholdP(z[P],z[T1]);                           
   nablaF->matrix0[M2Ac(5,T1,size)] = rholdT(z[P],z[T1]); 
  
   // F[6]: -z[hv] + h_v(z); 
   nablaF->matrix0[M2Ac(6,hv1,size)] = -1; 
   nablaF->matrix0[M2Ac(6,P,size)] = h_vdP(z[P],z[T1]);                           
   nablaF->matrix0[M2Ac(6,T1,size)] = h_vdT(z[P],z[T1]); 

   // F[7]: -z[hl] + h_l(z);
   nablaF->matrix0[M2Ac(7,hl1,size)] = -1; 
   nablaF->matrix0[M2Ac(7,P,size)] = h_ldP(z[P],z[T1]);                           
   nablaF->matrix0[M2Ac(7,T1,size)] = h_ldT(z[P],z[T1]); 
  
    //////////// MIXTURE 2 ///////////////////////
   // F[8]: z[M]-z_prev[M] 
   nablaF->matrix0[M2Ac(8,M2,size)] = 1 + h*Cvap*z[x2];
   nablaF->matrix0[M2Ac(8,x2,size)] = h*Cvap*z[M2] ;

    
   // F[9]: z[H]-z_prev[H] - h*Wext; 
   nablaF->matrix0[M2Ac(9,M2,size)] = (z[U2]-z_prev[U2]) + z[U2] + h*z[hv2]*Cvap*z[x2] ;
   nablaF->matrix0[M2Ac(9,U2,size)] = (z[M2]-z_prev[M2]) + z[M2];
   nablaF->matrix0[M2Ac(9,hv2,size)] = h*Cvap*z[M2]*z[x2];
   nablaF->matrix0[M2Ac(9,x2,size)] = h*z[hv2]*Cvap*z[M2];
   nablaF->matrix0[M2Ac(9,T1,size)] = -h*10;
   nablaF->matrix0[M2Ac(9,T2,size)] = h*10;

   // F[10] : -z[V] + z[M]/z[rl] + z[l1]*(1/z[rv]-1/z[rl]) ;
   nablaF->matrix0[M2Ac(10,V2,size)] =  -1.0;
   nablaF->matrix0[M2Ac(10,M2,size)] =  1.0/z[rl2] + z[x2]*(1.0/z[rv2]-1.0/z[rl2]);
   nablaF->matrix0[M2Ac(10,rv2,size)] = - z[M2]*z[x2]/pow(z[rv2],2.0);
   nablaF->matrix0[M2Ac(10,rl2,size)] = ( z[M2]*z[x2]-z[M2])/pow(z[rl2],2.0);
   nablaF->matrix0[M2Ac(10,x2,size)] = z[M2]/z[rv2]-z[M2]/z[rl2]; 
    
   // F[11]: -z[H] + z[M]*( z[hl] - z[P]/z[rl] ) 
   //             + z[l1]*( (z[hv]-z[hl])+(z[P]/z[rl] - z[P]/z[rv]) )
   nablaF->matrix0[M2Ac(11,U2,size)] = -1; 
   nablaF->matrix0[M2Ac(11,rv2,size)] = z[x2]*z[P]/pow(z[rv2],2.0) ;
   nablaF->matrix0[M2Ac(11,rl2,size)] = z[P]/pow(z[rl2],2.0)*(1.0-z[x2]); 
   nablaF->matrix0[M2Ac(11,hv2,size)] = z[x2];
   nablaF->matrix0[M2Ac(11,hl2,size)] = 1.0-z[x2];
   nablaF->matrix0[M2Ac(11,P,size)] = -1.0/z[rl2] + z[x2]*(1/z[rl2] - 1/z[rv2]);
   nablaF->matrix0[M2Ac(11,x2,size)] = (z[hv2]-z[hl2])+(z[P]/z[rl2] - z[P]/z[rv2]);
    
   // F[12]:  -z[rv] + rhov(z)
   nablaF->matrix0[M2Ac(12,rv2,size)] =  -1;
   nablaF->matrix0[M2Ac(12,P,size)] = rhovdP(z[P],z[T2]);
   nablaF->matrix0[M2Ac(12,T2,size)] = rhovdT(z[P],z[T2]);   

   // F[13]:-z[rl] + rhol(z);
   nablaF->matrix0[M2Ac(13,rl2,size)] = -1; 
   nablaF->matrix0[M2Ac(13,P,size)] = rholdP(z[P],z[T2]);                           
   nablaF->matrix0[M2Ac(13,T2,size)] = rholdT(z[P],z[T2]); 
  
   // F[14]: -z[hv] + h_v(z); 
   nablaF->matrix0[M2Ac(14,hv2,size)] = -1; 
   nablaF->matrix0[M2Ac(14,P,size)] = h_vdP(z[P],z[T2]);                           
   nablaF->matrix0[M2Ac(14,T2,size)] = h_vdT(z[P],z[T2]); 

   // F[15]: -z[hl] + h_l(z);
   nablaF->matrix0[M2Ac(15,hl2,size)] = -1; 
   nablaF->matrix0[M2Ac(15,P,size)] = h_ldP(z[P],z[T2]);                           
   nablaF->matrix0[M2Ac(15,T2,size)] = h_ldT(z[P],z[T2]); 

   //////////////////  EXTERNAL //////////////////////////
    
   // F[16]: P(t) - Pext;
   nablaF->matrix0[M2Ac(16,P,size)] = 1; 
   

   ////////////////// COMPLEMENTARITY ////////////////////
     
   // F[17]:
   nablaF->matrix0[M2Ac(17,P,size)] = 1; 
   nablaF->matrix0[M2Ac(17,T1,size)] = -psatdT(z[T1]) ; 
   nablaF->matrix0[M2Ac(17,DP1,size)] = 1; 
    
   // F[18]:
   nablaF->matrix0[M2Ac(18,x1,size)] = -1;  

   // F[19]:
   nablaF->matrix0[M2Ac(19,P,size)] = 1; 
   nablaF->matrix0[M2Ac(19,T2,size)] = -psatdT(z[T2]) ; 
   nablaF->matrix0[M2Ac(19,DP2,size)] = 1; 
    
   // F[20]:
   nablaF->matrix0[M2Ac(20,x2,size)] = -1;  

}
*/

void NablaF(void * env, int size, double *z, NumericsMatrix* nablaF){

   /////////////////////   MIXTURE 1  ///////////////////////
   // F[0]: z[M1]-z_prev[M1] + h*Ccond*z[M1]*(1-z[x1]) - h*Cvap*z[M2]*z[x2];  
   nablaF->matrix0[M2Ac(0,M1,size)] = 1 + h*Ccond*(1-z[x1]);  
   nablaF->matrix0[M2Ac(0,M2,size)] = -h*Cvap*z[x2];
   nablaF->matrix0[M2Ac(0,x1,size)] = -h*Ccond*z[M1] ;
   nablaF->matrix0[M2Ac(0,x2,size)] = - h*Cvap*z[M2] ;

   // F[1]: z[M1]*(z[U1]-z_prev[U1]) + z[U1]*(z[M1]-z_prev[M1])         
   //         - h*Kvl*(z[T2]-z[T1]) + h*z[hl1]*Ccond*z[M1]*(1-z[x1])
   //         - h*z[hv2]*Cvap*z[M2]*z[x2];    
   nablaF->matrix0[M2Ac(1,M1,size)] = (z[U1]-z_prev[U1]) + z[U1] + h*z[hl1]*Ccond*(1-z[x1]) ;
   nablaF->matrix0[M2Ac(1,M2,size)] = -h*z[hv2]*Cvap*z[x2] ;
   nablaF->matrix0[M2Ac(1,U1,size)] = (z[M1]-z_prev[M1]) + z[M1];
   nablaF->matrix0[M2Ac(1,hl1,size)] = h*Ccond*z[M1]*(1-z[x1]);
   nablaF->matrix0[M2Ac(1,hv2,size)] = -h*Cvap*z[M2]*z[x2];
   nablaF->matrix0[M2Ac(1,x1,size)] = -h*z[hl1]*Ccond*z[M1];
   nablaF->matrix0[M2Ac(1,x2,size)] = -h*z[hv2]*Cvap*z[M2];
   nablaF->matrix0[M2Ac(1,T1,size)] = h*Kvl;
   nablaF->matrix0[M2Ac(1,T2,size)] = -h*Kvl;

   // F[2] : -z[V1] + z[M1]/z[rl1] + z[M1]*z[x1]*(1/z[rv1]-1/z[rl1]) ;
   nablaF->matrix0[M2Ac(2,V1,size)] =  -1.0;
   nablaF->matrix0[M2Ac(2,M1,size)] =  1.0/z[rl1] + z[x1]*(1.0/z[rv1]-1.0/z[rl1]);
   nablaF->matrix0[M2Ac(2,rv1,size)] = - z[M1]*z[x1]/pow(z[rv1],2.0);
   nablaF->matrix0[M2Ac(2,rl1,size)] = ( z[M1]*z[x1]-z[M1])/pow(z[rl1],2.0);
   nablaF->matrix0[M2Ac(2,x1,size)] = z[M1]/z[rv1]-z[M1]/z[rl1]; 
    
   // F[3]: -z[U1] + 1.0*( z[hl1] - z[P]/z[rl1] ) 
   //              + z[x1]*(    (z[hv1]-z[hl1])          
   //                         + (z[P]/z[rl1] - z[P]/z[rv1]) 
   //                       ); 
   nablaF->matrix0[M2Ac(3,U1,size)] = -1; 
   nablaF->matrix0[M2Ac(3,rv1,size)] = z[x1]*z[P]/pow(z[rv1],2.0) ;
   nablaF->matrix0[M2Ac(3,rl1,size)] = z[P]/pow(z[rl1],2.0)*(1.0-z[x1]); 
   nablaF->matrix0[M2Ac(3,hv1,size)] = z[x1];
   nablaF->matrix0[M2Ac(3,hl1,size)] = 1.0-z[x1];
   nablaF->matrix0[M2Ac(3,P,size)] = -1.0/z[rl1] + z[x1]*(1/z[rl1] - 1/z[rv1]);
   nablaF->matrix0[M2Ac(3,x1,size)] = (z[hv1]-z[hl1])+(z[P]/z[rl1] - z[P]/z[rv1]);
    
   // F[4]:  -z[rv1] + rhov(z[P],z[T1]);
   nablaF->matrix0[M2Ac(4,rv1,size)] =  -1;
   nablaF->matrix0[M2Ac(4,P,size)] = rhovdP(z[P],z[T1]);
   nablaF->matrix0[M2Ac(4,T1,size)] = rhovdT(z[P],z[T1]);   

   // F[5]: -z[rl1] + rhol(z[P],z[T1]);
   nablaF->matrix0[M2Ac(5,rl1,size)] = -1; 
   nablaF->matrix0[M2Ac(5,P,size)] = rholdP(z[P],z[T1]);                           
   nablaF->matrix0[M2Ac(5,T1,size)] = rholdT(z[P],z[T1]); 
  
   // F[6]: -z[hv1] + h_v(z[P],z[T1]); 
   nablaF->matrix0[M2Ac(6,hv1,size)] = -1; 
   nablaF->matrix0[M2Ac(6,P,size)] = h_vdP(z[P],z[T1]);                           
   nablaF->matrix0[M2Ac(6,T1,size)] = h_vdT(z[P],z[T1]); 

   // F[7]: -z[hl1] + h_l(z[P],z[T1]);   
   nablaF->matrix0[M2Ac(7,hl1,size)] = -1; 
   nablaF->matrix0[M2Ac(7,P,size)] = h_ldP(z[P],z[T1]);                           
   nablaF->matrix0[M2Ac(7,T1,size)] = h_ldT(z[P],z[T1]); 
  
    //////////// MIXTURE 2 /////////////////////// 
   // F[8] = z[M2]-z_prev[M2] + h*Cvap*z[M2]*z[x2] - h*Ccond*z[M1]*(1-z[x1]);
   nablaF->matrix0[M2Ac(8,M2,size)] = 1 + h*Cvap*z[x2];
   nablaF->matrix0[M2Ac(8,x2,size)] = h*Cvap*z[M2] ;
   nablaF->matrix0[M2Ac(8,M1,size)] = - h*Ccond*(1-z[x1]) ;
   nablaF->matrix0[M2Ac(8,x1,size)] = h*Ccond*z[M1];
    
   //  F[9] = z[M2]*(z[U2]-z_prev[U2]) + z[U2]*(z[M2]-z_prev[M2])         // -h*Wext;
   //                                  - h*Kvl*(z[T1]-z[T2]) + h*z[hv2]*Cvap*z[M2]*z[x2]
   //                                  - h*z[hl1]*Ccond*z[M1]*(1-z[x1]);  
   nablaF->matrix0[M2Ac(9,M2,size)] = (z[U2]-z_prev[U2]) + z[U2] + h*z[hv2]*Cvap*z[x2] ;
   nablaF->matrix0[M2Ac(9,U2,size)] = (z[M2]-z_prev[M2]) + z[M2];
   nablaF->matrix0[M2Ac(9,hv2,size)] = h*Cvap*z[M2]*z[x2];
   nablaF->matrix0[M2Ac(9,T1,size)] = - h*Kvl;
   nablaF->matrix0[M2Ac(9,T2,size)] = h*Kvl;
   nablaF->matrix0[M2Ac(9,x2,size)] = h*z[hv2]*Cvap*z[M2];
   nablaF->matrix0[M2Ac(9,hl1,size)] = - h*Ccond*z[M1]*(1-z[x1]);
   nablaF->matrix0[M2Ac(9,M1,size)] = - h*z[hl1]*Ccond*(1-z[x1]);
   nablaF->matrix0[M2Ac(9,x1,size)] = h*z[hl1]*Ccond*z[M1];

   //  F[10] = -z[V2] + z[M2]/z[rl2] + z[M2]*z[x2]*(1/z[rv2]-1/z[rl2]) ;
   nablaF->matrix0[M2Ac(10,V2,size)] =  -1.0;
   nablaF->matrix0[M2Ac(10,M2,size)] =  1.0/z[rl2] + z[x2]*(1.0/z[rv2]-1.0/z[rl2]);
   nablaF->matrix0[M2Ac(10,rv2,size)] = - z[M2]*z[x2]/pow(z[rv2],2.0);
   nablaF->matrix0[M2Ac(10,rl2,size)] = ( z[M2]*z[x2]-z[M2])/pow(z[rl2],2.0);
   nablaF->matrix0[M2Ac(10,x2,size)] = z[M2]/z[rv2]-z[M2]/z[rl2]; 
    
   //F[11] = -z[U2] + 1.0*( z[hl2] - z[P]/z[rl2] ) 
   //               + z[x2]*( (z[hv2]-z[hl2])          // Rq: T is not directly in this eq.
   //                         + (z[P]/z[rl2] - z[P]/z[rv2]) ); 
   nablaF->matrix0[M2Ac(11,U2,size)] = -1; 
   nablaF->matrix0[M2Ac(11,rv2,size)] = z[x2]*z[P]/pow(z[rv2],2.0) ;
   nablaF->matrix0[M2Ac(11,rl2,size)] = z[P]/pow(z[rl2],2.0)*(1.0-z[x2]); 
   nablaF->matrix0[M2Ac(11,hv2,size)] = z[x2];
   nablaF->matrix0[M2Ac(11,hl2,size)] = 1.0-z[x2];
   nablaF->matrix0[M2Ac(11,P,size)] = -1.0/z[rl2] + z[x2]*(1/z[rl2] - 1/z[rv2]);
   nablaF->matrix0[M2Ac(11,x2,size)] = (z[hv2]-z[hl2])+(z[P]/z[rl2] - z[P]/z[rv2]);
    
   // F[12] = -z[rv2] + rhov(z[P],z[T2]);
   nablaF->matrix0[M2Ac(12,rv2,size)] =  -1;
   nablaF->matrix0[M2Ac(12,P,size)] = rhovdP(z[P],z[T2]);
   nablaF->matrix0[M2Ac(12,T2,size)] = rhovdT(z[P],z[T2]);   

   //  F[13] = -z[rl2] + rhol(z[P],z[T2]);
   nablaF->matrix0[M2Ac(13,rl2,size)] = -1; 
   nablaF->matrix0[M2Ac(13,P,size)] = rholdP(z[P],z[T2]);                           
   nablaF->matrix0[M2Ac(13,T2,size)] = rholdT(z[P],z[T2]); 
  
   // F[14] = -z[hv2] + h_v(z[P],z[T2]); 
   nablaF->matrix0[M2Ac(14,hv2,size)] = -1; 
   nablaF->matrix0[M2Ac(14,P,size)] = h_vdP(z[P],z[T2]);                           
   nablaF->matrix0[M2Ac(14,T2,size)] = h_vdT(z[P],z[T2]); 

   // F[15] = -z[hl2] + h_l(z[P],z[T2]); 
   nablaF->matrix0[M2Ac(15,hl2,size)] = -1; 
   nablaF->matrix0[M2Ac(15,P,size)] = h_ldP(z[P],z[T2]);                           
   nablaF->matrix0[M2Ac(15,T2,size)] = h_ldT(z[P],z[T2]); 

   //////////////////  EXTERNAL //////////////////////////
    
   // F[16]: P(t) - Pext;
   nablaF->matrix0[M2Ac(16,P,size)] = 1; 
   

   ////////////////// COMPLEMENTARITY ////////////////////
     
   // F[17] = z[P] - psat(z[T1]) + z[DP1];
   nablaF->matrix0[M2Ac(17,P,size)] = 1; 
   nablaF->matrix0[M2Ac(17,T1,size)] = -psatdT(z[T1]) ; 
   nablaF->matrix0[M2Ac(17,DP1,size)] = 1; 
    
   // F[18] = 1.0-z[x1];
   nablaF->matrix0[M2Ac(18,x1,size)] = -1;  

   // F[19] = z[P] - psat(z[T2]) + z[DP2];
   nablaF->matrix0[M2Ac(19,P,size)] = 1; 
   nablaF->matrix0[M2Ac(19,T2,size)] = -psatdT(z[T2]) ; 
   nablaF->matrix0[M2Ac(19,DP2,size)] = 1; 
    
   // F[20] = 1.0-z[x2];
   nablaF->matrix0[M2Ac(20,x2,size)] = -1;  

   NM_display(nablaF);
}

/* Initialization of the implicit euler discrete problem as an MCP */
static MixedComplementarityProblem * create_mcp_drumboiler(void)
{
    /* Create a MixedComplementarityProblem */
    MixedComplementarityProblem* problem = (MixedComplementarityProblem *)malloc(sizeof(MixedComplementarityProblem));

    int n = 21;
    
    problem->n1 = 17; //nb equality constraint
    problem->n2 = 4; // nb complementarity constraint
    problem->compute_Fmcp = &F ;
    problem->compute_nabla_Fmcp = &NablaF ;
    problem->nabla_Fmcp =  NM_create(NM_DENSE, n, n); 
    problem->env = NULL;
 
     
    return problem;
}

/* Compute Scaling weights */
double* compute_scaling(double* z, int size){

  double *scalingVector = (double *)calloc(size,sizeof(double));

  // double* F_curr = (double *)calloc(size,sizeof(double));
  // F(NULL,size,z,F_curr) ;
  // F_i norm
  // for(unsigned int i=0;j<size;i++)
  // {
  //   scalingVector[i] = fabs(F_curr);
  // }

  // NablaF_i max norm
  NumericsMatrix * nablaF = NM_create(NM_DENSE, size, size);

  NablaF(NULL, size, z, nablaF);

  NM_max_abs_by_rows(nablaF, scalingVector);

  return scalingVector;
}

/* Balanced MCP problem */
void Fbalanced(void* env, int size, double *z, double * Fout);
void Fbalanced(void* env, int size, double *z, double * Fout){
    F(env, size, z, Fout);
    for(unsigned int i=0;i<size;i++)
    {  
      Fout[i] = Fout[i]/scaling[i];
    }

}


/* Balanced Jacobian of the MCP problem */
void NablaFbalanced(void * env, int size, double *z, NumericsMatrix* nablaF);
void NablaFbalanced(void * env, int size, double *z, NumericsMatrix* nablaF){   
    NablaF(env, size, z, nablaF);
    for(unsigned int i=0;i<size;i++)
    {  
      for(unsigned int j=0;j<size;j++)
      {    
          NM_zentry(nablaF, i, j, NM_get_value(nablaF, i, j)/scaling[i] );
      }
    }
}

/* Generation of the balanced MCP problem (caled at each time step)  */
MixedComplementarityProblem * balanced_mcp_drumboiler(void)
{
    /* Create a MixedComplementarityProblem */
    MixedComplementarityProblem* problem = (MixedComplementarityProblem *)malloc(sizeof(MixedComplementarityProblem));

    int n = 21;
    
    problem->n1 = 17; //nb equality constraint
    problem->n2 = 4; // nb complementarity constraint
    problem->compute_Fmcp = &Fbalanced ;
    problem->compute_nabla_Fmcp = &NablaFbalanced ;
    problem->nabla_Fmcp =  NM_create(NM_DENSE, n, n); 
    problem->env = NULL;
 
     
    return problem;
}


/* Balanced Newton solving */
static int solve_balanced_mcp_newton(int solverId)
{
  printf("solve_balanced_mcp_newton() starts for solver %s.\n", solver_options_id_to_name(solverId));

  int info = 1 ;

  
  /*Create original (not balanced) MCP to solve -- not needed (n1 and n2 could be obtained differently)*/
  MixedComplementarityProblem* problem_orig = create_mcp_drumboiler();
  
  printf("Original MCP created \n");  

  /* Set solver options */
  //SolverOptions options; // OLD CODE methods
  // options.solverId = solverId;
  // mcp_setDefaultSolverOptions(problem, &options);

  /* FB solver */
  SolverOptions * options = solver_options_create(solverId);
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 50;
  options->dparam[SICONOS_DPARAM_TOL] = 7e-7;
  
  /* VERBOSE MODE */
  numerics_set_verbose(2);

  int size = problem_orig->n1 + problem_orig->n2 ;

  // printf("Compute scaling \n");  

  scaling = compute_scaling(z_curr, size);
  cout << "scaling = "  << endl;
  for(unsigned int i=0; i<size; i++)
  {
    cout << scaling[i] << " " ;
  }
  cout << "\n "  << endl;

  /* Create the Banlanced MCP problem */
  MixedComplementarityProblem* problem_balanced = balanced_mcp_drumboiler();

  double * z = (double *)calloc(size,sizeof(double));
  double * w = (double *)calloc(size,sizeof(double));

  Fbalanced(NULL,size,z_curr,w ); // initial guess of w = Fbalanced(z_curr)

  // display constraints values at current (soon previous) point
  cout << "Fbalanced(z_curr) = "  << endl;
  for(unsigned int i=0; i<size; i++)
  {
    cout << w[i] << " " ;
  }
  cout << endl;

  // saving current -> prev
  memcpy(z_prev, z_curr, size * sizeof(double)); // prev = curr   
  // initialization of next step newton z_prev -> z 
  memcpy(z, z_prev, size * sizeof(double)); // initial guess of z

  cout << "[solve_mcp_newton] -- z_prev = "  << endl;
  for(unsigned int i=0;i<size;i++)
  {
    cout << z_prev[i]<< " " ;
  }
  cout << endl;

  // computing next step -> stored in z 
  // info = mcp_driver(problem, z , w,  &options); // OLD CODE
  info = mcp_driver(problem_balanced, z , w,  options); 

  // storing z -> z_curr
  memcpy(z_curr, z, size * sizeof(double));   // curr = next
  free(z); // z memory
  free(w); // free w memory
  free(scaling); // free scaling vector

  double * nablaFbalanced_z = problem_balanced->nabla_Fmcp->matrix0;
  double cond_num = cond(nablaFbalanced_z, size, size);
  cout << "########## nablaFbalanced(z) conditionning: " << cond_num  << endl;
  // NM_display(problem_balanced->nabla_Fmcp);

  /* Freeing the MCP*/
  delete(problem_balanced->nabla_Fmcp);
  free(problem_balanced);
}



/* Solving one step problem*/ 
static int solve_mcp_newton(int solverId)
{
  printf("test_mcp_newton() starts for solver %s.\n", solver_options_id_to_name(solverId));

  int info = 1 ;

  MixedComplementarityProblem* problem = create_mcp_drumboiler();
  
  /* Set solver options */
  //SolverOptions options; // OLD CODE methods
  // options.solverId = solverId;
  // mcp_setDefaultSolverOptions(problem, &options);
  /* FB solver */
  SolverOptions * options = solver_options_create(solverId);
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 2;
  options->dparam[SICONOS_DPARAM_TOL] = 7e-7;
  
  /* VERBOSE MODE */
  numerics_set_verbose(3);

  int size = problem->n1 + problem->n2 ;
  double * z = (double *)calloc(size,sizeof(double));
  double * w = (double *)calloc(size,sizeof(double));

  F(NULL,size,z_curr,w ); // initial guess of w = F(z_curr)


  // display constraints values at current (soon previous) point
  cout << "F(z_curr) = "  << endl;
  for(unsigned int i=0; i<size; i++)
  {
    cout << w[i] << " " ;
  }
  cout << endl;

   cout << "err on liquid mass balance eq. , F[8] = "<< w[8] << endl;
    cout << "err on liquid energy balance  eq. , F[9] = "<< w[9] << endl;
  // saving current -> prev
  memcpy(z_prev, z_curr, size * sizeof(double)); // prev = curr   
  // initialization of next step newton z_prev -> z 
  memcpy(z, z_prev, size * sizeof(double)); // initial guess of z

  cout << "[solve_mcp_newton] -- z_prev = "  << endl;
  for(unsigned int i=0;i<size;i++)
  {
    cout << z_prev[i]<< " " ;
  }
  cout << endl;

  // computing next step -> stored in z 
  // info = mcp_driver(problem, z , w,  &options); // OLD CODE
  info = mcp_driver(problem, z , w,  options); 

  // storing z -> z_curr
  memcpy(z_curr, z, size * sizeof(double));   // curr = next
  free(z); // z memory
  free(w); // free w memory

  //  double * nablaFz = problem->nabla_Fmcp->matrix0;
  //  double cond_num = cond(nablaFz, size, size);
  //  cout << "########## nablaF(z) conditionning: " << cond_num  << endl;
//   NM_display(problem->nabla_Fmcp);

  /* Freeing the MCP*/
  delete(problem->nabla_Fmcp);
  free(problem);
}

int main(void)
{
    int n = 21;
    z_prev = (double *) calloc(n,sizeof(double));
    z_curr = (double *) calloc(n,sizeof(double));
   
    h = 10.0;//0.05; //s // TODO smaller h leads to non-convergence of the Newton
    Tf = 10.0;//2500.0; // s
    //Vtot =  0.001; //  m^3 // TODO UNUSED
    Wext = -1000.0;//1000000.0;//100000.0;  // J/s (Watt) // TODO UNUSED
    Pext = 106074.0; //1060740; // 106074.0;
    Cvap = 0.09; // evaporation mixture transition rate
    Ccond = 0.01; // condensation mixture transition rate
    // Text1 = 600; // TODO UNUSED
    // /* Initial Conditions */  
    // ######## MIXTURE 1 (vapor)
    double M10 = 0.5;//1.0; // kg
    double x10 = 1.0; //0.000589*M0; // Kg
    double T10 = 600; //600.0;  //K 
    z_prev[T1] = T10;
    z_prev[M1] = M10; 

    double P0 = psat(z_prev[T1]);
    P0 = Pext;
    z_prev[P] = P0;
    // double rv0 = 10.0; // kg/m^3
    double rv10 = rhov(z_prev[P],z_prev[T1]);   // kg/m^3
    double rl10 = rhol(z_prev[P],z_prev[T1]);  // kg/m^3
    double hv10 = h_v(z_prev[P],z_prev[T1]);   // J/kg
    double hl10 = h_l(z_prev[P],z_prev[T1]); // J/kg

    double V10 = M10*x10/rv10 + M10*(1-x10)/rl10;

    double DP10;
    if (P0-psat(z_prev[T1]) < 0.0)
      DP10 = -1.0*(P0-psat(z_prev[T1]));  // [P-Psat]⁻ in Pa
    else   
      DP10 = 0.0;  // [P-Psat]⁻ in Pa

    double U10 = (1.0-x10)*(h_l(z_prev[P],z_prev[T1]) - P0/rhol(z_prev[P],z_prev[T1])) + x10*(h_v(z_prev[P],z_prev[T1]) - P0/rhov(z_prev[P],z_prev[T1])); // J
    
    z_prev[P] = P0;  
    z_prev[V1] = V10;
    z_prev[M1] = M10;  z_prev[U1] = U10;  z_prev[rv1] = rv10;
    z_prev[rl1] = rl10;  z_prev[hv1] = hv10;  z_prev[hl1] = hl10;
    z_prev[T1] = T10; z_prev[x1] = x10;  z_prev[DP1] = DP10; 

    // ######## MIXTURE 2 (liquid)

    // double M20 = 0.5;//1.0; // kg 
    double M20 = 0.0;//1.0; // kg
    double x20 =  1.0; //0.000589*M0; // Kg
    // double T20 = 373.156;  //K 
    double T20 = 600.0;  //K 
    z_prev[T2] = T20 ;
    z_prev[M2] = M20; 

    // double rv0 = 10.0; // kg/m^3
    double rv20 = rhov(z_prev[P],z_prev[T2]);   // kg/m^3
    double rl20 = rhol(z_prev[P],z_prev[T2]);  // kg/m^3
    double hv20 = h_v(z_prev[P],z_prev[T2]);   // J/kg
    double hl20 = h_l(z_prev[P],z_prev[T2]); // J/kg

    double V20 = M20*x20/rv20 + M20*(1-x20)/rl20;

    double DP20;
    if (P0-psat(z_prev[T2]) < 0.0)
      DP20 = -1.0*(P0-psat(z_prev[T2]));  // [P-Psat]⁻ in Pa
    else   
      DP20 = 0.0;  // [P-Psat]⁻ in Pa

    double U20 = (1.0-x20)*(h_l(z_prev[P],z_prev[T2]) - P0/rhol(z_prev[P],z_prev[T2])) + x20*(h_v(z_prev[P],z_prev[T2]) - P0/rhov(z_prev[P],z_prev[T2])); // J

    z_prev[V2] = V20;
    z_prev[M2] = M20;  z_prev[U2] = U20;  z_prev[rv2] = rv20;
    z_prev[rl2] = rl20;  z_prev[hv2] = hv10;  z_prev[hl2] = hl20;
    z_prev[T2] = T20; z_prev[x2] = x20;  z_prev[DP2] = DP20; 


    memcpy(z_curr, z_prev, n*sizeof(double));

    cout << endl;
    cout << "### general varaibles" << endl;
    cout << " V  = ??" << " |  Tw = ??" <<   " |  P = " <<  z_curr[P]  << endl; 
    cout << "### vapour mixture  varaibles" << endl;
    cout << " Mvd  = " <<  z_curr[M1] << " |  uvd = " <<  z_curr[U1]  <<  " |  Vvd = "<<  z_curr[V1]  << endl;
    cout << " TV  = " <<  z_curr[T1] << " |  zV = ???"   <<  " |  xv = " <<  z_curr[x1] << endl;
    cout << " rv  = " <<  z_curr[rv1] << " |  rd = " <<  z_curr[rl1] << endl;
    cout << " hv  = " <<  z_curr[hv1] << " |  hd = " <<  z_curr[hl1] << endl;
    cout << " DPV = " << z_curr[DP1] << " |  P-PsatV = " << z_curr[P]-psat(z_curr[T1]) << endl;
    cout << "###  liquid mixture  varaibles" << endl;
    cout << " Mlb  = " <<  z_curr[M2] << " |  ulb = " <<  z_curr[U2]   << " |  V2 = " <<  z_curr[V2] <<  endl;
    cout << " TL  = " <<  z_curr[T2] << " |  zL = ???"  <<  " |  xb = " <<  z_curr[x2] << endl;
    cout << " rl  = " <<  z_curr[rl2] << " |  rb = " <<  z_curr[rv2] << endl;
    cout << " hl  = " <<  z_curr[hl2] << " |  hb = " <<  z_curr[hv2] << endl;
    cout << " DPL = " << z_curr[DP2] << " |  P-PsatL = " << z_curr[P]-psat(z_curr[T2]) << endl;
    cout << endl;

    cout << "Initial conditions: "<< endl;
    for(unsigned int i=0;i<n;i++)
    {
      cout << z_curr[i] << " " ;
    } 
    cout << endl;


    int N = ceil(Tf/h)+1;
    //int N = 1;

    unsigned int outputSize = n+1;
    SimpleMatrix dataPlot(N + 1, outputSize);
    dataPlot(0, 0) = 0.0;
    for(unsigned int i=0;i<n;i++)
    {
      dataPlot(0, i+1) = z_curr[i] ;
    }
    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    // boost::timer::auto_cpu_timer time;
    
    for(unsigned int j=0;j<N;j++)
    {
    //   int info = solve_mcp_newton(SICONOS_MCP_NEWTON_FB_FBLSA);   // FB solver option test
     int info = solve_mcp_newton(SICONOS_MCP_NEWTON_MIN_FBLSA);              // Problem without matrix balancing
    //int info = solve_balanced_mcp_newton(SICONOS_MCP_NEWTON_MIN_FBLSA); // Calling balancing method 

      if(info==0)
      {
        cout << "[in main()] next-step z_curr ="<< endl;
        for(unsigned int i=0;i<n;i++)
        {
          cout << z_curr[i] << " ";
        }
        cout << endl;
        cout << "### general varaibles" << endl;
        cout << " V  = ??" << " |  Tw = ??" <<   " |  P = " <<  z_curr[P]  << endl; 
        cout << "### vapour mixture  varaibles" << endl;
        cout << " Mvd  = " <<  z_curr[M1] << " |  uvd = " <<  z_curr[U1]  <<  " |  Vvd = "<<  z_curr[V1]  << endl;
        cout << " TV  = " <<  z_curr[T1] << " |  zV = ???"   <<  " |  xv = " <<  z_curr[x1] << endl;
        cout << " rv  = " <<  z_curr[rv1] << " |  rd = " <<  z_curr[rl1] << endl;
        cout << " hv  = " <<  z_curr[hv1] << " |  hd = " <<  z_curr[hl1] << endl;
        cout << " DPV = " << z_curr[DP1] << " |  P-PsatV = " << z_curr[P]-psat(z_curr[T1]) << endl;
        cout << "###  liquid mixture  varaibles" << endl;
        cout << " Mlb  = " <<  z_curr[M2] << " |  ulb = " <<  z_curr[U2]   << " |  V2 = " <<  z_curr[V2] <<  endl;
        cout << " TL  = " <<  z_curr[T2] << " |  zL = ???"  <<  " |  xb = " <<  z_curr[x2] << endl;
        cout << " rl  = " <<  z_curr[rl2] << " |  rb = " <<  z_curr[rv2] << endl;
        cout << " hl  = " <<  z_curr[hl2] << " |  hb = " <<  z_curr[hv2] << endl;
        cout << " DPL = " << z_curr[DP2] << " |  P-PsatL = " << z_curr[P]-psat(z_curr[T2]) << endl;
        cout << endl;
        cout << "info = " << info << endl;
        cout << "t =  " << ((j+1)*h) << " seconds " << endl;
        cout << "####\n" << endl;
        // if(z_curr[T]>374.881)
            // return 0;
      }

      dataPlot(k, 0) = (k*h);
      for(unsigned int i=0;i<n;i++)
      {
        dataPlot(k, i+1) = z_curr[i] ;
      }
      
      k++;
    }
  cout << "End of computation - Number of iterations done: " << k - 1 << endl;
  cout << "Computation Time :" << endl;
  // time.report();
  // --- Output files ---
  cout << "====> Output file writing ..." << endl;
  dataPlot.resize(k, outputSize);
  ioMatrix::write("result-evap-process-mcp-6eq-isobare-testemptyphase.dat", "ascii", dataPlot, "noDim");
      
}
