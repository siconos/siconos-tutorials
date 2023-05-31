/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifdef _WIN32
#define SICONOS_EXPORT extern "C" __declspec(dllexport)
#else
#define SICONOS_EXPORT extern "C"
#endif
#include <stdio.h>
#include <math.h>



SICONOS_EXPORT void TwoDofsOscillatorB(double time, unsigned int sizeOfx, double *fExt, unsigned int sizeZ, double* z)
{
  printf("In oscillatorB : sizeOfx = %i \n", sizeOfx);
  double omega = 0.3;
  double f1 =1.0;
  double f2 =0.0;
  
  fExt[0] = 0.0;
  fExt[1] = f1 * cos(omega*time);
  fExt[2] = 0.0;
  fExt[3] = f2 * cos(omega*time);

  printf("In oscillatorB : fext = %e, %e \n", fExt[1], fExt[3] );
  
}

