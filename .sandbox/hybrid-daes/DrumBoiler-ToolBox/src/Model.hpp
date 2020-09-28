#ifndef MODEL_HPP
#define MODEL_HPP

#include "model_assets.hpp"

class FiniteDifference
{
    public:

        FiniteDifference(){}
        ~FiniteDifference(){}

        double finiteDifference(void* env, double z_i, double zprev_i);  
        double finiteDifferenceDz(void* env, double z_i, double zprev_i);  

}; 


class Model
{
    public:

        

        Model(){};
        ~Model(){delete(z);};
        enum Variables { };
        enum Equations { };

        double* get_z();
        int get_nb_eqs();
        int get_nb_vars();
        int get_nb_diffeqs();
        int get_nb_equaeqs();
        int get_nb_compleqs();
        virtual std::vector<std::string> get_vars_name();
        virtual std::string get_var_name(const int& index);
        virtual std::vector<std::string> get_eqs_name();
        virtual std::string get_eq_name(const int& index);

        virtual void display_model_state(double h, int iteration);

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
        virtual void doubleMixtureRHS(double* z, double* rhs){};

        // Jacobian of the MCS right hand side
        virtual void doubleMixtureRHS_Jacobian(double* z, NumericsMatrix* rhs_jac){};

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
        virtual void doubleMixtureLHS(void* env, FiniteDifference* diffMethod, double* z, double* lhs){};

        virtual void doubleMixtureLHS_Jacobian(void* env, FiniteDifference* diffMethod,
                                               double* z, NumericsMatrix* lhs_jac){};


    protected:

        // std::map<std::string,int> variables;
        int nb_vars; // TODO nb_vars enforced at model definition -- not dynamic model
        int nb_tot_eqs;
        std::vector<std::string> vars_name;
        std::vector<std::string> eqs_name;
        int nb_diff_eqs;
        int nb_equa_eqs;
        int nb_compl_eqs;

        double* z; 

};

#endif //MODEL_HPP