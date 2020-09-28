#ifndef MODEL_ASSETS_HPP
#define MODEL_ASSETS_HPP

#include <boost/timer/timer.hpp>
#include <map>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <chrono>
#include "cond.h"

#include "NumericsMatrix.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

inline void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}


// typedef std::pair<std::string, int> varType;

	
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


/* 2d matrix conversion in 1d array 
(Column first, indexes starting at 0) */
inline int M2Ac(int i, int j, int nb_line)
{
    int index = nb_line*j + i;
    return index;
}

// inline std::ostream& operator<<(std::ostream& out, std::vector<double> v) {
// 	std::cout.precision(8);
//     out << " [ ";
// 	for (unsigned int i = 0; i < v.size(); i++)
// 		out << v[i] << " ";
//     out << " ] ";
// 	out << std::endl;
// 	return out;
// }

#endif // MODEL_ASSETS_HPP