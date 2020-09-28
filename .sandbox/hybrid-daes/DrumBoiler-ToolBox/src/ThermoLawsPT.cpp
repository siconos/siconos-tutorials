
#include "ThermoLawsPT.hpp"

double ThLPT_perfGaz_incompLiquid::rhov(const double &P_, const double &T_)
{
    return P_/(R*T_); 
}

double ThLPT_perfGaz_incompLiquid::rhovdP(const double &P_, const double &T_)
{
    return 1/(R*T_);
}

double ThLPT_perfGaz_incompLiquid::rhovdT(const double &P_, const double &T_)
{
    return -P_/(R*pow(T_,2.0));
}

/* Volumic mass of water --  */
double ThLPT_perfGaz_incompLiquid::rhol(const double &P_, const double &T_)
{
    return (-0.0025*(T_-273.156) - 0.1992)*(T_-273.156) + 1004.4;
 
}

double ThLPT_perfGaz_incompLiquid::rholdP(const double &P_, const double &T_)
{
    return 0;

}

double ThLPT_perfGaz_incompLiquid::rholdT(const double &P_, const double &T_)
{
    return (-2*0.0025*(T_-273.156) - 0.1992);
}


/* Specific enthalpy of vapor (J/kg) in fonction of Pressure (Pa) and Temperature (K) */
double ThLPT_perfGaz_incompLiquid::h_v(const double &P_, const double &T_)
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

double ThLPT_perfGaz_incompLiquid::h_vdP(const double &P_, const double &T_)
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

double ThLPT_perfGaz_incompLiquid::h_vdT(const double &P_, const double &T_)
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
double ThLPT_perfGaz_incompLiquid::h_l(const double &P_, const double &T_)
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


double ThLPT_perfGaz_incompLiquid::h_ldP(const double &P_, const double &T_)
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

double ThLPT_perfGaz_incompLiquid::h_ldT(const double &P_, const double &T_)
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
double ThLPT_perfGaz_incompLiquid::psat(const double &T_)
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

double ThLPT_perfGaz_incompLiquid::psatdT(const double &T_)
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

/* specific internal energy equation*/ 
double ThLPT_perfGaz_incompLiquid::u(const double &P_, const double &h_, const double &r_)
{
    return h_- P_/r_;
}
double ThLPT_perfGaz_incompLiquid::udP(const double &P_, const double &h_, const double &r_){
    return -1.0/r_;
}
double ThLPT_perfGaz_incompLiquid::udh(const double &P_, const double &h_, const double &r_){
    return 1.0;
}
double ThLPT_perfGaz_incompLiquid::udr(const double &P_, const double &h_, const double &r_){
    return P_/pow(r_,2.0);
}

/* volume equation */
double ThLPT_perfGaz_incompLiquid::vol(const double &r_, const double &M_){
    return M_/r_;
}
double ThLPT_perfGaz_incompLiquid::voldr(const double &r_, const double &M_){
    return -1.0*M_/pow(r_,2.0);
}
double ThLPT_perfGaz_incompLiquid::voldM(const double &r_, const double &M_){
    return 1.0/r_;
}


/* Volumic mass of water --  */
double ThLPT_CompressLiquid::rhol(const double &P_, const double &T_)
{
    
    double c00=-1.829480e+03;
    double c10=3.054319e-05;
    double c01=2.781477e+01;
    double c20=3.828328e-14;
    double c11=-3.043374e-07;
    double c02=-9.811186e-02;
    double c30=8.084907e-22;
    double c21=-6.519321e-16;
    double c12=1.176861e-09;
    double c03=1.443719e-04;
    double c40=-6.330092e-32;
    double c31=-4.041731e-24;
    double c22=2.507897e-18;
    double c13=-2.063206e-12;
    double c04=-6.651181e-08;
    double c50=-3.239489e-39;
    double c41=1.530767e-33;
    double c32=4.801971e-27;
    double c23=-2.757629e-21;
    double c14=1.380801e-15;
    double c05=-1.948287e-11;

    return c00 + c10*P_ + c01*T_ + c20*pow(P_,2.0) 
               + c11*P_*T_ + c02*pow(T_,2.0) + c30*pow(P_,3.0) 
               + c21*pow(P_,2.0)*T_ + c12*P_*pow(T_,2.0) 
               + c03*pow(T_,3.0) + c40*pow(P_,4.0) + c31*pow(P_,3.0)*T_ 
               + c22*pow(P_,2.0)*pow(T_,2.0) + c13*P_*pow(T_,3.0) 
               + c04*pow(T_,4.0) + c50*pow(P_,5.0) + c41*pow(P_,4.0)*T_ 
               + c32*pow(P_,3.0)*pow(T_,2.0) + c23*pow(P_,2.0)*pow(T_,3.0) 
               + c14*P_*pow(T_,4.0) + c05*pow(T_,5.0);
}

double ThLPT_CompressLiquid::rholdP(const double &P_, const double &T_)
{
    double c00=-1.829480e+03;
    double c10=3.054319e-05;
    double c01=2.781477e+01;
    double c20=3.828328e-14;
    double c11=-3.043374e-07;
    double c02=-9.811186e-02;
    double c30=8.084907e-22;
    double c21=-6.519321e-16;
    double c12=1.176861e-09;
    double c03=1.443719e-04;
    double c40=-6.330092e-32;
    double c31=-4.041731e-24;
    double c22=2.507897e-18;
    double c13=-2.063206e-12;
    double c04=-6.651181e-08;
    double c50=-3.239489e-39;
    double c41=1.530767e-33;
    double c32=4.801971e-27;
    double c23=-2.757629e-21;
    double c14=1.380801e-15;
    double c05=-1.948287e-11;

    return c10 + c20*2.0*P_ + c11*T_ 
               + c30*3.0*pow(P_,2.0) 
               + c21*2.0*P_*T_ + c12*pow(T_,2.0) 
               + c40*4.0*pow(P_,3.0) + c31*3.0*pow(P_,2.0)*T_ 
               + c22*2.0*P_*pow(T_,2.0) + c13*pow(T_,3.0) 
               + c50*5.0*pow(P_,4.0) + c41*4.0*pow(P_,3.0)*T_ 
               + c32*3.0*pow(P_,2.0)*pow(T_,2.0) + c23*2.0*P_*pow(T_,3.0) 
               + c14*pow(T_,4.0);  


}

double ThLPT_CompressLiquid::rholdT(const double &P_, const double &T_)
{

    double c00=-1.829480e+03;
    double c10=3.054319e-05;
    double c01=2.781477e+01;
    double c20=3.828328e-14;
    double c11=-3.043374e-07;
    double c02=-9.811186e-02;
    double c30=8.084907e-22;
    double c21=-6.519321e-16;
    double c12=1.176861e-09;
    double c03=1.443719e-04;
    double c40=-6.330092e-32;
    double c31=-4.041731e-24;
    double c22=2.507897e-18;
    double c13=-2.063206e-12;
    double c04=-6.651181e-08;
    double c50=-3.239489e-39;
    double c41=1.530767e-33;
    double c32=4.801971e-27;
    double c23=-2.757629e-21;
    double c14=1.380801e-15;
    double c05=-1.948287e-11;

    return c01 + c11*P_ + c02*2.0*T_ + c21*pow(P_,2.0) 
               + c12*P_*2.0*T_ + c03*3.0*pow(T_,2.0) 
               + c31*pow(P_,3.0) + c22*pow(P_,2.0)*2.0*T_ 
               + c13*P_*3.0*pow(T_,2.0) + c04*4.0*pow(T_,3.0) 
               + c41*pow(P_,4.0) + c32*pow(P_,3.0)*2.0*T_ 
               + c23*pow(P_,2.0)*3.0*pow(T_,2.0) + c14*P_*4.0*pow(T_,3.0) 
               + c05*5.0*pow(T_,4.0); 


}


double ThLPT_fittedGaz::rhov(const double &P_, const double &T_)
{

    double c00=8.573340e+01;
    double c10=2.032089e-05;
    double c01=-3.131904e-01;
    double c20=9.742006e-15;
    double c11=-2.602836e-08;
    double c02=3.165324e-04;
    double c30=-1.477236e-22;
    double c21=7.411890e-18;
    double c12=8.461200e-12;
    double c03=-9.503033e-08;

    return c00 + c10*P_ + c01*T_ + c20*pow(P_,2.0) 
               + c11*P_*T_ + c02*pow(T_,2.0) 
               + c30*pow(P_,3.0) + c21*pow(P_,2.0)*T_ 
               + c12*P_*pow(T_,2.0) + c03*pow(T_,3.0);

  
}

double ThLPT_fittedGaz::rhovdP(const double &P_, const double &T_)
{

    double c00=8.573340e+01;
    double c10=2.032089e-05;
    double c01=-3.131904e-01;
    double c20=9.742006e-15;
    double c11=-2.602836e-08;
    double c02=3.165324e-04;
    double c30=-1.477236e-22;
    double c21=7.411890e-18;
    double c12=8.461200e-12;
    double c03=-9.503033e-08;
  
    return c10 + c20*2.0*P_ + c11*T_ 
               + c30*3.0*pow(P_,2.0) + c21*2.0*P_*T_ 
               + c12*pow(T_,2.0);
 
}

double ThLPT_fittedGaz::rhovdT(const double &P_, const double &T_)
{


    double c00=8.573340e+01;
    double c10=2.032089e-05;
    double c01=-3.131904e-01;
    double c20=9.742006e-15;
    double c11=-2.602836e-08;
    double c02=3.165324e-04;
    double c30=-1.477236e-22;
    double c21=7.411890e-18;
    double c12=8.461200e-12;
    double c03=-9.503033e-08;

    return c01 + c11*P_ + c02*2.0*T_ 
               + c21*pow(P_,2.0) + c12*P_*2.0*T_ 
               + c03*3.0*pow(T_,2.0); 
}
