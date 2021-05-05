#ifndef THERMOLAWSPT_HPP
#define THERMOLAWSPT_HPP

#include "model_assets.hpp"

class ThermoLawsPT{
    public:
        ThermoLawsPT(){};
        ~ThermoLawsPT(){};

        /* Volumic mass of vapor (kg/m^3) in function of Pressure (Pa) and Temperature (K)*/
        virtual double rhov(const double &P_, const double &T_) = 0;
        virtual double rhovdP(const double &P_, const double &T_) = 0;
        virtual double rhovdT(const double &P_, const double &T_) = 0;

        /* Volumic mass of liquid water (kg/m^3) in function of Temperature (K) only */
        virtual double rhol(const double &P_, const double &T_) = 0;
        virtual double rholdP(const double &P_, const double &T_) = 0;
        virtual double rholdT(const double &P_, const double &T_) = 0;

        /* Specific enthalpy of vapor (J/kg) in fonction of Pressure (Pa) and Temperature (K) */
        virtual double h_v(const double &P_, const double &T_) = 0;
        virtual double h_vdP(const double &P_, const double &T_) = 0;
        virtual double h_vdT(const double &P_, const double &T_) = 0;

        /* Specific enthalpy of liquid (J/kg) in fonction of Pressure (Pa) and Temperature (K) */
        virtual double h_l(const double &P_, const double &T_) = 0;
        virtual double h_ldP(const double &P_, const double &T_) = 0;
        virtual double h_ldT(const double &P_, const double &T_) = 0;

        /* Saturation pressure (Pa) in fonction of Temperature (K)*/
        virtual double psat(const double &T_) = 0;
        virtual double psatdT(const double &T_) = 0;


        /* specific internal energy equation*/ 
        virtual double u(const double &P_, const double &h_, const double &r_) = 0;
        virtual double udP(const double &P_, const double &h_, const double &r_) = 0;
        virtual double udh(const double &P_, const double &h_, const double &r_) = 0;
        virtual double udr(const double &P_, const double &h_, const double &r_) = 0;

        /* volume equation */
        virtual double vol(const double &r_, const double &M_) = 0;
        virtual double voldr(const double &r_, const double &M_) = 0;
        virtual double voldM(const double &r_, const double &M_) = 0;  
};


class ThLPT_perfGaz_incompLiquid : public ThermoLawsPT{
    public:
        ThLPT_perfGaz_incompLiquid(){};
        ~ThLPT_perfGaz_incompLiquid(){};

        /* Volumic mass of vapor (kg/m^3) in function of Pressure (Pa) and Temperature (K)*/
        virtual double rhov(const double &P_, const double &T_);
        virtual double rhovdP(const double &P_, const double &T_);
        virtual double rhovdT(const double &P_, const double &T_);

        /* Volumic mass of liquid water (kg/m^3) in function of Temperature (K) only */
        virtual double rhol(const double &P_, const double &T_);
        virtual double rholdP(const double &P_, const double &T_);
        virtual double rholdT(const double &P_, const double &T_);

        /* Specific enthalpy of vapor (J/kg) in fonction of Pressure (Pa) and Temperature (K) */
        virtual double h_v(const double &P_, const double &T_);
        virtual double h_vdP(const double &P_, const double &T_);
        virtual double h_vdT(const double &P_, const double &T_);

        /* Specific enthalpy of liquid (J/kg) in fonction of Pressure (Pa) and Temperature (K) */
        virtual double h_l(const double &P_, const double &T_);
        virtual double h_ldP(const double &P_, const double &T_);
        virtual double h_ldT(const double &P_, const double &T_);

        /* Saturation pressure (Pa) in fonction of Temperature (K)*/
        virtual double psat(const double &T_);
        virtual double psatdT(const double &T_);


        /* specific internal energy equation*/ 
        virtual double u(const double &P_, const double &h_, const double &r_);
        virtual double udP(const double &P_, const double &h_, const double &r_);
        virtual double udh(const double &P_, const double &h_, const double &r_);
        virtual double udr(const double &P_, const double &h_, const double &r_);

        /* volume equation */
        virtual double vol(const double &r_, const double &M_);
        virtual double voldr(const double &r_, const double &M_);
        virtual double voldM(const double &r_, const double &M_);


    protected:
        const double R = 461.527;  

};

class ThLPT_fittedGaz : virtual public ThLPT_perfGaz_incompLiquid{
    public:
        ThLPT_fittedGaz(){};
        ~ThLPT_fittedGaz(){};

        /* Volumic mass of vapor (kg/m^3) in function of Pressure (Pa) and Temperature (K)*/
        virtual double rhov(const double &P_, const double &T_);
        virtual double rhovdP(const double &P_, const double &T_);
        virtual double rhovdT(const double &P_, const double &T_);

};


class ThLPT_CompressLiquid : virtual public ThLPT_perfGaz_incompLiquid{
    public:
        ThLPT_CompressLiquid(){};
        ~ThLPT_CompressLiquid(){};

         /* Volumic mass of liquid water (kg/m^3) in function of Pressure (Pa) and Temperature (K) */
        virtual double rhol(const double &P_, const double &T_);
        virtual double rholdP(const double &P_, const double &T_);
        virtual double rholdT(const double &P_, const double &T_);

};

class ThLPT_Fitted : public ThLPT_fittedGaz, public ThLPT_CompressLiquid{
    public:
        ThLPT_Fitted(){};
        ~ThLPT_Fitted(){};
};


#endif // THERMOLAWSPT_HPP