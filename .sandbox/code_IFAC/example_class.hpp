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

#include "SiconosKernel.hpp"
#include "lcp_cst.h"

using namespace std;

enum ProblemType { SLIDING_CROSSING, SLIDING_REPULSIVE, NON_NULL_A };

class Problem
{
public:

    SP::SiconosMatrix A;
    SP::SimpleMatrix R;
    SP::SiconosVector b;
    SP::SimpleMatrix C;
    SP::SimpleMatrix D;
    SP::SiconosVector e;
    SP::SimpleMatrix M;
    SP::SiconosVector init;
    double t0;
    double T;
    ProblemType type;
    unsigned int dimX = 3; // forced for example
    unsigned int dimLambda = 3; // forced for example

    Problem(SP::SiconosMatrix lA, SP::SimpleMatrix lR, SP::SiconosVector lb,
            SP::SimpleMatrix lC, SP::SimpleMatrix lD, SP::SiconosVector le,
            SP::SimpleMatrix lM, SP::SiconosVector linit,
            double lt0, double lT, ProblemType ltype ){

    A       = lA;
    R       = lR;
    b       = lb;
    C       = lC;
    D       = lD;
    e       = le;
    M       = lM;
    init    = linit;
    t0      = lt0;
    T       = lT;
    type    = ltype;
    }

    ~Problem(){

    }

};

class ConvergenceTest
{
public:


    ConvergenceTest(Problem* lproblem, vector<double> ltime_steps){
        problem = lproblem;
        time_steps = ltime_steps;  
    }

    /*Simulate a nonsmooth DAE for a given time step
      The size of the problem is fixed  
    */   
    SP::SimpleMatrix simulate(Problem* p, double h, int *out){
        unsigned int dimX       = p->dimX;    // Dimension of the system state variables
        unsigned int dimLambda  = p->dimLambda;    // Dimension of the system lambda variables
        double t0       = p->t0;          // initial computation time
        double T        = p->T;         // final computation tim

        // -------------------------
        // --- Dynamical systems ---
        // -------------------------
        SP::FirstOrderLinearTIDS dyn(new FirstOrderLinearTIDS(p->init,p->A));
        dyn->setbPtr(p->b);
        dyn->setMPtr(p->M);
        // -------------------------
        // --- LCP Relation ---
        // -------------------------
        // Relation LCP lhs
        SP::FirstOrderLinearTIR relation(new FirstOrderLinearTIR(p->C, p->R) );
        relation->setDPtr(p->D);
        relation->setePtr(p->e);
        // NonSmooth law: LCP
        SP::NonSmoothLaw nslaw(new ComplementarityConditionNSL(dimLambda));
        // interaction 
        SP::Interaction inter(new Interaction(nslaw, relation));
        // -----------------------------
        // --- Siconos Model Entity ---
        // ----------------------------
        SP::NonSmoothDynamicalSystem switch_dae(new NonSmoothDynamicalSystem(t0, T));
        // add the dynamical system in the non smooth dynamical system
        switch_dae->insertDynamicalSystem(dyn);
        // link the interaction and the dynamical system
        switch_dae->link(inter, dyn);
        // -----------------------------
        // --- Simulation Definition ---
        // -----------------------------
        // -- (1) OneStepIntegrators --
        double theta = 1.0;
        double gamma = 1.0;
        SP::EulerMoreauOSI osi(new EulerMoreauOSI(theta,gamma));
        // -- (2) Time discretisation --
        SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));
        // -- (3) one step non smooth problem
        SP::LCP osnspb(new LCP(SICONOS_LCP_ENUM));
        osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS] = 1; 
        osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_SKIP_TRIVIAL] = SICONOS_LCP_SKIP_TRIVIAL_NO;
        osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_SEED] = 0; // SEED
        // osnspb->setNumericsVerboseMode(true);
        // -- (4) Simulation setup with (1) (2) (3)
        SP::TimeStepping s(new TimeStepping(switch_dae, td, osi, osnspb));

        // =========================== End of model definition ===========================
        //   // ================================= Computation =================================
        int N = ceil((T - t0) / h)+1; // Number of time steps
        
        // --- Get the values to be plotted ---
        // -> saved in a matrix dataPlot
        unsigned int outputSize = 7;
        SP::SimpleMatrix dataPlot( new SimpleMatrix(N + 1, outputSize));

        SP::SiconosVector x = dyn->x();
        SP::SiconosVector lambda = inter->lambda(0);

        (*dataPlot)(0, 0) = switch_dae->t0();
        (*dataPlot)(0, 1) = (*x)(0);   // x1
        (*dataPlot)(0, 2) = (*x)(1);   // x2
        (*dataPlot)(0, 3) = (*x)(2);   // z
        (*dataPlot)(0, 4) = (*lambda)(0);
        (*dataPlot)(0, 5) = (*lambda)(1);       // 
        (*dataPlot)(0, 6) = (*lambda)(2);       // 
        // --- Time loop ---
        // ==== Simulation loop - Writing without explicit event handling =====
        int k = 1;
        boost::progress_display show_progress(N);
        boost::timer time;
        time.restart();
        SiconosVector xk = SiconosVector(*(p->init));
        SiconosVector xsol = SiconosVector(*(p->init));
        SiconosVector lsol = SiconosVector({0, 0, 0});
        SiconosVector diff = SiconosVector(*(p->init));
        double norm = INFINITY;
        int nb_modes = 8;
        double norm_diff;
        while (s->hasNextEvent())
        {
            
            for(int i=0;i<nb_modes;i++)
            { 
               osnspb->numericsSolverOptions()->iparam[SICONOS_LCP_IPARAM_ENUM_SEED] = i; // SEED 
               s->computeOneStep();
               diff = *x-xk;
               norm_diff = diff.norm2();
               if(norm_diff<norm)
               {
                    xsol = *x;
                    lsol = *lambda; 
                    norm = norm_diff;
               } 
            //    osnspb->display();
            }
            (*x) = xsol;
            (*lambda) = lsol;
            // --- Get values to be plotted ---
            (*dataPlot)(k, 0) =  s->nextTime();
            (*dataPlot)(k, 1) = (*x)(0);
            (*dataPlot)(k, 2) = (*x)(1);
            (*dataPlot)(k, 3) = (*x)(2);
            (*dataPlot)(k, 4) = (*lambda)(0);
            (*dataPlot)(k, 5) = (*lambda)(1);
            (*dataPlot)(k, 6) = (*lambda)(2);
            xk = xsol;
            norm =INFINITY;
            s->nextStep();
            ++show_progress;
            k++;

        }
        // cout  << "End of computation - Number of iterations done: " << k - 1 << endl;
        cout << "Computation Time " << time.elapsed()  << endl;
        *out = k;
        return dataPlot;

    }

    /* run the simulation and compare the results with the corresponding analytic solution for a given set of time steps*/    
    vector<double> run(){
        
        vector<double> errors(time_steps.size()); 
        for(unsigned int i=0;i<time_steps.size();i++)
        {
            double h = time_steps[i]; 
            SP::SimpleMatrix current_results;
            int k;
            current_results = simulate(problem,h,&k);
            // int N = ceil((problem->T - problem->t0) / h)+1; // not satisfying
            double x1f = (*current_results)(k-1,1);
            double x2f = (*current_results)(k-1,2);
            // temporary as SiconosVector({x1f, x2f}) is not accepted by IDE
            vector<double> tmp({x1f, x2f}); 
            SiconosVector xf = SiconosVector(tmp);
            // diff.norm2()
            switch(problem->type) {
                case SLIDING_CROSSING :
                    tmp = vector<double>({6.5*2./3.,6.5*2./3. + 1});
                    break;    
                case SLIDING_REPULSIVE : 
                    tmp = vector<double>({8/3.,1.+8/3.});
                    break;    
                case NON_NULL_A :
                    tmp = vector<double>({1.8867421308653363, 2.8867421308653363});
                    break;    
                default:
                    cout << "Invalid Problem Type" << endl;
                    break;    
            }
            SiconosVector diff = SiconosVector(tmp);
            // cout << xf << endl;
            // cout << diff << endl;
            // cout << k-1 << endl;
            diff = diff - xf;
            double current_error = diff.norm2();
            errors[i] = current_error;
            cout << current_error << endl;
        }
        return errors;
    }

    ~ConvergenceTest(){
        delete(problem);
    }

private:

    vector<double> time_steps;
    Problem* problem;
};
