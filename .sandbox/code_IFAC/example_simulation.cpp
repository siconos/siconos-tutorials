/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
  Example of Linear switching DAE with a constraint defined by a multivalued operator
  Here we considere the parametrization of the sliding_repulsive case
*/

#include "example_class.hpp"

int main(int argc, char* argv[])
{
    try
    {

        unsigned int dimX       = 3;    // Dimension of the system state variables
        unsigned int dimLambda  = 3;    // Dimension of the system lambda variables

        double t0       = 0.0;          // initial computation time
        double T        = 10.0;         // final computation time 

        double x1_0     = 0.0;         // initial condition in state variable x1
        double x2_0     = 0.0;          // initial condition in state variable x2
        double z_0      = 0;            // initial condition in algebraic variable z

        SP::SiconosVector init(new SiconosVector({x1_0, x2_0, z_0}));

        SP::SiconosMatrix A( new SimpleMatrix(dimX,dimX) );         
	
	// This vector B is specific to the sliding repulsive case
        double B0 = -1.0;
        double B1 = 0.5;
        A->setRow(0,SiconosVector({0.0, 0.0, B0}));
        A->setRow(1,SiconosVector({0.0, 0.0, B1}));
        A->setRow(2,SiconosVector({1.0, -1.0, 0.0}));

        SP::SimpleMatrix M(new SimpleMatrix(dimX,dimX));
        (*M)(0,0) = 1.0;
        (*M)(1,1) = 1.0;

        SP::SiconosVector b(new SiconosVector({1.0, 0.0, 1.0}));

        SP::SimpleMatrix C( new SimpleMatrix(dimLambda,dimX) );
        (*C)(0,0) = 2.0;
        (*C)(1,0) = 1.0;

        SP::SimpleMatrix D( new SimpleMatrix(dimLambda,dimLambda) );
        (*D)(0,0) = 1.0;
        (*D)(1,2) = 1.0;
        (*D)(2,1) = -1.0;

        SP::SimpleMatrix R( new SimpleMatrix(dimX,dimLambda) );
        (*R)(2,0) = 1.0;
        (*R)(2,1) = -1.0;

        SP::SiconosVector e(new SiconosVector({0.0, 0.0, 2.0}));
        ProblemType type = SLIDING_REPULSIVE; // Specific to sliding repulsive case (not critical for simulation)
        Problem* problem = new Problem( A, R, b, C, D, e, M, init, t0, T, type);
        vector<double> time_steps({ 0.009, 0.09, 0.9}); // Time-steps only needed for constructor
        ConvergenceTest test(problem, time_steps);
        int k;
        SP::SimpleMatrix results = test.simulate(problem,0.2,&k); // Time step used is 0.2 as specifed here
        cout << (*results) << endl;

        unsigned int outputSize = 7;
        results->resize(k, outputSize);
        ioMatrix::write("simulations_results.dat", "ascii", (*results), "noDim");

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
