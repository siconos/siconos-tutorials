#include "Model.hpp"



// finite difference discretisation of the differentiation operator (scaled by 'h')
double FiniteDifference::finiteDifference(void* env, double z_i, double zprev_i){

    return z_i - zprev_i;

}

// partial derivartive of finite difference 
double FiniteDifference::finiteDifferenceDz(void* env, double z_i, double zprev_i){
    return 1.0;
}

// void Model::set_vars(const std::vector<std::string> & var_list){

//     int n = var_list.size();
//     for(int i=0;i<n;i++){
//         variables[var_list[n]]=n;
//     }
// }

std::vector<std::string> Model::get_vars_name(){
    return vars_name;
}

std::string Model::get_var_name(const int& index){
    // // Iterate through all elements in std::map and search for the passed element
    // std::map<std::string, int>::iterator it = variables.begin();
    // while(it != variables.end())
    // {
    //     if(it->second == index)
    //     return it->first;
    //     it++;
    // }
    return vars_name[index];
}

std::vector<std::string> Model::get_eqs_name(){
    return eqs_name;
}

std::string Model::get_eq_name(const int& index){
    // // Iterate through all elements in std::map and search for the passed element
    // std::map<std::string, int>::iterator it = variables.begin();
    // while(it != variables.end())
    // {
    //     if(it->second == index)
    //     return it->first;
    //     it++;
    // }
    return eqs_name[index];
}



// int Model::id(const std::string& name){
//     std::map<std::string, int>::iterator it =  variables.find(name);
//     if(it != variables.end())
//         return it->second;
//     else
//     {
//         std::cout << "Variable" << name << "Not Found !\n";
//     }
           
// }

double* Model::get_z(){
    return z;
}


int Model::get_nb_eqs(){
    return nb_tot_eqs;
}
int Model::get_nb_vars(){
    return nb_vars;
}

int Model::get_nb_diffeqs(){
    return nb_diff_eqs;
}


int Model::get_nb_equaeqs(){
    return nb_equa_eqs;
}
int Model::get_nb_compleqs(){
    return nb_compl_eqs;
}

void Model::display_model_state(double h, int iteration){

    std::cout << "#### Iteration: "<< iteration << std::endl;
    for(int i=0;i<nb_vars;i++){
        if(i%2==0)
            std::cout << vars_name[i] << "  = "<<  z[i]<< " | "; 
        else
            std::cout << vars_name[i] << "  = " <<  z[i]  << std::endl;    
    }
    std::cout << std::endl;
    std::cout << "t =  " << ((iteration)*h) << " seconds " << std::endl;
    std::cout << "####\n" << std::endl;
}