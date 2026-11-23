#ifndef MODEL_COMMON_HPP
#define MODEL_COMMON_HPP

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <boost/timer/timer.hpp>
#include <iostream>
#include <map>
#include <vector>

#include "NumericsMatrix.h"
#include "cond.h"

/* 2d matrix conversion in 1d array
(Column first, indexes starting at 0) */
inline int M2Ac(int i, int j, int nb_line) {
  int index = nb_line * j + i;
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

#endif  // MODEL_COMMON_HPP
