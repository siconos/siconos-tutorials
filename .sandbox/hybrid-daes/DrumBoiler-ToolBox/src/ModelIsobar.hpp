#ifndef MODELISOBAR_HPP
#define MODELISOBAR_HPP

#include "Model.hpp"
#include "ThermoLawsPT.hpp"

class ModelNoWallIsobar: public Model
{
    enum Variables {P,M1,V1,U1,rv1,rl1,hv1,hl1,T1,M2,V2,U2,rv2,rl2,hv2,hl2,T2,x1,DP1,x2,DP2};
    enum Equations {dM1,dU1,dM2,dU2,V1eq,u1eq,rv1eq,rl1eq,hv1eq,hl1eq,V2eq,u2eq,rv2eq,rl2eq,
                    hv2eq,hl2eq,isoP_eq,x1_nseq,DP1_nseq,x2_ns,DP2_nseq};

    typedef struct
    {
        double Ccond; // condensation rate constant
        double Cvap; 
        double Kvl; 
        double xqext1; // specific (mass-ratio) heat energy input
        double xqext2;
        double qext1; // specific (mass) heat energy input
        double qext2;
        double Qext1; // constant heat energy input
        double Qext2;
        double Pext; // Fixed Pressure (= P(0) in isobar case)

    }ModelConstants;
    public:
        ModelNoWallIsobar();
        ~ModelNoWallIsobar(){};

        void display_model_state(double h, int iteration);

        int initialize_model(double P, double M1, double x1, double T1,
                             double M2, double x2, double T2);

        void set_ExternalHeatExchange(double xqext1, double xqext2, double qext1,
                                      double qext2, double Qext1, double Qext2);

        // TODO is this a good way to architecture the thermolaws/model links ?
        void set_lawsPT(ThermoLawsPT* laws_){laws = laws_; };

        // void get_constants(const ModelConstants & constants);


        /* Model of condensation mass transfert */
        double massCond(const double &M_, const double &x_);
        double massConddx(const double &M_, const double &x_);
        double massConddM(const double &M_, const double &x_);

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



        double energIn(const double & M_,const double & Q_);
        double energIndM(const double & M_,const double & Q_);

        /* Model of external energy transfert  - Constant = Qext 
        TODO remove if not used.
        TODO Add other external interaction 
        */
        double energIn1(const double & x_,const double & Q_);
        double energIn1dx(const double & x_,const double & Q_);

        double energIn2(const double & x_,const double & Q_);
        double energIn2dx(const double & x_,const double & Q_);


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
        virtual void doubleMixtureRHS(double* z, double* rhs);

        // Jacobian of the MCS right hand side
        virtual void doubleMixtureRHS_Jacobian(double* z, NumericsMatrix* rhs_jac);

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
        virtual void doubleMixtureLHS(void* env, FiniteDifference* diffMethod, double* z, double* lhs);

        virtual void doubleMixtureLHS_Jacobian( void* env, FiniteDifference* diffMethod,
                                                double* z, NumericsMatrix* lhs_jac);


    protected:

        ModelConstants constants;
        ThermoLawsPT* laws;

};


#endif //MODELISOBAR_HPP
