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

/*
!\file absolute-value and sign-mlcp.cpp
 
  Example of MLCP formulation of the absolute value and sign example.
  Can be generalized for Linear switching DAE .

  Current version simulate only 1 time step.

  TODO:

  - iterations

  - mlcp solvers have to send all solutions (currently only 1)

  Examples of mlcp in cpp are in siconos/numerics/src/MLCP/test/main_mlcp.cpp
  Examples of mlcp in cpp are in siconos/numerics/swig/test

*/

#include "SiconosKernel.hpp"
#include "MixedLinearComplementarityProblem.h"
#include "MLCP_Solvers.h"
#include "NumericsVerbose.h"

using namespace std;

/*
 ******************************************************************************
 */
void printSolution(const char *name, int n, int m, int NbLines, double *z, double *w)
{
  int i;
  printf(" *** ************************************** ***\n");
  printf(" ****** z = ********************************\n");
  for (i = 0 ; i < n + m ; i++)
    printf("%s: %14.7e  \n", name, z[i]);
  printf(" ****** w = ********************************\n");
  for (i = 0 ; i < m ; i++)
    printf("%s: %14.7e  \n", name, w[i + (NbLines - m)]);
}


/*      |A  B|
    M = |C  D|

    z = [x,lambda]

    q = |a|
        |b| 

    MLCP Type 1: Mz + q = w = [w1,
                                w2]
                0 =< w2 _|_ lambda >= 0            
    */
void generate_Type1_mlcp(MixedLinearComplementarityProblem* mlcp)
{
    mlcp->isStorageType2 = 0;
    mlcp->isStorageType1 = 1;
    int n = mlcp->n; // nbEq
    int m = mlcp->m; // nbIneq
    double * A = mlcp->A;
    double * C = mlcp->C;
    double * B = mlcp->B;
    double * D = mlcp->D;
    double * a = mlcp->a;
    double * b = mlcp->b;

    NumericsMatrix* M = NM_create(NM_DENSE, n+m, n+m);
    double * _M = M->matrix0;
    
    /* Remainder Column first */
    // setting A in M (A is a n*n mat)
    // m[i + (i/nbline[A])*nbline[C]] = A[i]
    for(unsigned int i=0; i<n*n; i++)
    {
        _M[i + (i/n)*m] = A[i];
    }

    // setting C in M (C is a mxn mat)
    for(unsigned int i=0; i<m*n; i++)
    {
        _M[i + ((i/m) + 1)*n] = C[i];
    }

    // setting B in M (B is a nxm mat)
    for(unsigned int i=0; i<n*m; i++)
    {
        _M[i + ((n+m)*n) + (i/n)*m] = B[i];
    }

    // setting D in M (D is a mxm mat)
    for(unsigned int i=0; i<m*m; i++)
    {
        _M[i + ((n+m)*n) + ((i/m) + 1)*n] = D[i];
    }

    double * q = (double *) calloc(n+m,sizeof(double));

    for(unsigned int i=0; i<n; i++)
    {
        q[i] = a[i];
    }
    for(unsigned int i=0; i<m; i++)
    {
        q[i+n] = b[i];
    }

    mlcp->M = M;
    mlcp->q = q;
}

int main(int argc, char* argv[])
{
    try
    {
        numerics_set_verbose(1);
        /* Input parametrization */
        double B1 = -1.0;
        double B2 = 0.5; 

        /* Initial Conditions */
        double x1k = -1; // previous position of x1 
        double x2k = 0;  // previous position of x2
        double h = 0.05; // time step

        /* problem structure */      
        int nbEq = 3; 
        int nbIneq = 3; 
        // REMINDER: Column first convention

        /* A = [-1 0 B1; 0 -1 B2; 1 -1 0] */ 
        double * A = (double *) calloc(nbEq*nbEq , sizeof(double)); // smarter way to init a double* in c++ ?
        A[0] = -1; A[1] = 0; A[2] = 1; // 1st column
        A[3] = 0; A[4] = -1; A[5] = -1; // 2nd column
        A[6] = B1; A[7] = B2; A[8] = 0; // 3rd column

        /* B = [0 0 0; 0 0 0; 1 -1 0]*/
        double * B = (double *) calloc(nbEq*nbIneq, sizeof(double));
        B[2] = 1.0; 
        B[5] = -1.0; 


        double * C = (double *) calloc(nbIneq*nbEq, sizeof(double));
        C[0] = 2; C[1] = 1;

        double * D = (double *) calloc(nbIneq*nbIneq, sizeof(double));
        D[0] = 1; D[5] = -1; D[7] = 1;

        /* b = [0; 0; 2]*/
        double * b = (double *) calloc(nbIneq, sizeof(double));
        b[2] = 2.0; 

        double * a = NULL; // initialization to NULL


        bool hasNextStep = true; // temp


        int info = 0;
        double *w = (double *) calloc(nbIneq+nbEq, sizeof(double));
        double *z = (double *) calloc(nbIneq+nbEq, sizeof(double));


        while(hasNextStep)
        {
            hasNextStep = false;
 
            /* definition of mlcp at step k corresponding to backward euler discritization */            
            MixedLinearComplementarityProblem* mlcp = new MixedLinearComplementarityProblem();
            
            /* Problem defined in format:
                w1 = Az1 + Bz2 + a
                w2 = Cz1 + Dz2 + b
                w1 = 0
                0 =< w2 _|_ z2 >= 0*/
            mlcp->isStorageType2 = 1;
            mlcp->m = nbIneq;
            mlcp->n = nbEq;
 
            mlcp->A = A;
            mlcp->B = B;
            mlcp->C = C;
            mlcp->D = D;
            mlcp->b = b;   
            
            /* a = [x1k + h; x2k; 1] depends of k*/
            if(a ==NULL)
                a = (double *) calloc(nbEq, sizeof(double));
            else
            {
                free(a);
                a = (double *) calloc(nbEq, sizeof(double));
            }
            a[0] = x1k + h; a[1] = x2k; a[2] = 1; 
            mlcp->a = a;         

            generate_Type1_mlcp(mlcp);
            mixedLinearComplementarity_display(mlcp);


            SolverOptions* options = new SolverOptions(); 

            // mixedLinearComplementarity_enum_setDefaultSolverOptions(mlcp, options); // seems to fails
            options->solverId = SICONOS_MLCP_DIRECT_ENUM;
            mixedLinearComplementarity_setDefaultSolverOptions(mlcp, options);
            cout << "Solver default options DONE" << endl; 
            options->dparam[0] = 1e-12;
            mlcp_driver_init(mlcp, options);
            cout << "Solver driver init DONE" << endl;      
            // cout << "dWork = " << options->dWork <<  endl;
 

            info = mlcp_driver(mlcp, z, w, options);
            // cout << "info = " <<  info << endl;
            // printSolution("ENUM", nbEq, nbIneq, nbEq+nbIneq, z, w);
            // mixedLinearComplementarity_deleteDefaultSolverOptions(mlcp, options);
            // mlcp_driver_reset(mlcp, options);
            
            // mlcp_enum(mlcp, z, w, &info, options);
            cout << "info = " <<  info << endl;
            cout << "z = [ ";
            for(unsigned int i=0; i<nbIneq+nbEq; i++)
            {
                cout << z[i] << " , ";
            }

            cout << " ]" << endl;

            cout << "w = [ ";
            for(unsigned int i=0; i<nbIneq+nbEq; i++)
            {
                cout << w[i] << " , ";
            }

            cout << " ]" << endl;


            cout << "Computation DONE" << endl;
        }
    }

    catch (SiconosException& e)    
    {
    cerr << e.report() << endl;
    return 1;
    }

    catch (...) 
    {
    cerr << "Exception caught in absolute-value-and-sign-lcp.cpp" << endl;
    return 1;
    }


}