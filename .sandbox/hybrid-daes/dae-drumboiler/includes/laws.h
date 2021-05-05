#ifndef LAWS_H
#define LAWS_H

#include "model_common.h"

/* Volumic mass of vapor (kg/m^3) in function of Pressure (Pa) and Temperature (K)*/
double rhov(const double &P_, const double &T_);
double rhovdP(const double &P_, const double &T_);
double rhovdT(const double &P_, const double &T_);

/* Volumic mass of liquid water (kg/m^3) in function of Temperature (K) only */
double rhol(const double &P_, const double &T_);
double rholdP(const double &P_, const double &T_);
double rholdT(const double &P_, const double &T_);

/* Specific enthalpy of vapor (J/kg) in fonction of Pressure (Pa) and Temperature (K) */
double h_v(const double &P_, const double &T_);
double h_vdP(const double &P_, const double &T_);
double h_vdT(const double &P_, const double &T_);

/* Specific enthalpy of liquid (J/kg) in fonction of Pressure (Pa) and Temperature (K) */
double h_l(const double &P_, const double &T_);
double h_ldP(const double &P_, const double &T_);
double h_ldT(const double &P_, const double &T_);

/* Saturation pressure (Pa) in fonction of Temperature (K)*/
double psat(const double &T_);
double psatdT(const double &T_);


/* specific internal energy equation*/ 
double u(const double &P_, const double &h_, const double &r_);
double udP(const double &P_, const double &h_, const double &r_);
double udh(const double &P_, const double &h_, const double &r_);
double udr(const double &P_, const double &h_, const double &r_);

/* volume equation */
double vol(const double &r_, const double &M_);
double voldr(const double &r_, const double &M_);
double voldM(const double &r_, const double &M_);



#endif // LAWS_H