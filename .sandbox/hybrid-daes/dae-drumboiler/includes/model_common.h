#ifndef MODEL_COMMON_H
#define MODEL_COMMON_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <boost/timer/timer.hpp>
#include <iostream>
#include <map>
#include <vector>

#include "NumericsMatrix.h"

/* 2d matrix conversion in 1d array
(Column first, indexes starting at 0) */
int M2Ac(int i, int j, int nb_line) {
  int index = nb_line * j + i;
  return index;
}

#endif  // MODEL_COMMON_H
