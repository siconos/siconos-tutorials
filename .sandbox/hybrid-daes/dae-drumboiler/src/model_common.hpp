#ifndef MODEL_COMMON_HPP
#define MODEL_COMMON_HPP

#include <boost/timer/timer.hpp>
#include <map>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cond.h"

#include "NumericsMatrix.h"

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

#endif // MODEL_COMMON_HPP