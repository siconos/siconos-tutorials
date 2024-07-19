#include "SiconosKernel.hpp"
#include "FirstOrderNonLinearDS.hpp"
#include "SolverOptions.h"
#include "NumericsVerbose.h"

#include <math.h>

static double t0 = 0.0;     // Start time
static double T  = 10.0;    // End time
// Systems parameters
static double A1 = -2.12e3;
static double B1 = -3.0e-3;
static double C1 = -0.707;
static double C3 = -1;
static double B3 = 1e-4;
static double B3X3B3 = 1e-4;

/* Update the offset qk of the LCP*/
int update(double *q, std::vector<double> y, double h, int k){
    try{
      /* 
        qk = (h*Dhat - B3X3B3)*Lk + hC1*x1k + h*(C3*u2dotk +  C2u2k) 
      */ 
      double tk = t0 + k*h;
      double u2dotk = -1e-4*10*sin(10*tk);
      q[0] = -B3X3B3*y[1] + C1*h*y[0] + C3*h*u2dotk; 
      return 0;
    }
    catch (...)
    {
      std::cout << "unsafe memory access in update qk function !\n";
      return 1;
    }
}

SimpleMatrix* compute_reference(double h, SolverOptions * options, NumericsMatrix* M)
{
  try
    {
        double lambda10 = 0.0;
        double x10 = 0.01;

        size_t size = 1;    
        std::vector<double> init = {x10,lambda10};
        double * q = (double *)calloc(size,sizeof(double));

        update(q, init, h, 0);

        LinearComplementarityProblem* lcp = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));
        lcp->M = M;
        lcp->q = q;
        lcp->size = size;

        size_t k = 0;
        size_t N = int((T - t0) / h) + 2;

        size_t outputSize = 1 + 4;
        SimpleMatrix* dataPlot = new SimpleMatrix(N, outputSize); 
        
        double u1k;
        double u2k;
        double u2dotk;

        (*dataPlot)(0, 0) = t0;
        (*dataPlot)(0, 1) = x10; // x1
        u2k = 1e-4*cos(10*(t0+h*k));
        (*dataPlot)(0, 2) = -B3*lambda10 - u2k; // x2 = - B3*lambda - u2(tk)
        (*dataPlot)(0, 3) = lambda10; // lambda1 
        // data_plot[k, 4] = z0 computed afterward as it is depending on \dot{lambda} 

        // ==== Simulation loop - Writing without explicit event handling =====
        k+=1;
        std::vector<double> yk = {x10,lambda10};
        double * lambdak = (double *)calloc(size,sizeof(double));
        double * wk = (double *)calloc(size,sizeof(double));
        int info = 0;
        double dlambda1 = 0.0;
        while (k<N)
        { 

          u1k = 0.301e4*sin(t0+h*(k-1)); // explicit input
          u2k = 1e-4*cos(10*(t0+h*k));
          u2dotk = -10*1e-4*sin(10*(t0+h*(k-1)));
          // update lcp
          update(q,yk,h,k-1);    // explicit input
          info = linearComplementarity_driver(lcp, lambdak , wk,  options);
          // --- new solutions ---
          yk[0] = (1+A1*h)*yk[0] + B1*h*yk[1] + h*u1k; // Fully explicit compute of x1[k+1]
          yk[1] = lambdak[0]; // Lambda1[k+1]

          (*dataPlot)(k, 0) = (*dataPlot)(k-1, 0) + h;
          (*dataPlot)(k, 1) = yk[0]; // x1
          (*dataPlot)(k, 2) = -B3*yk[1] - u2k; // x2
          (*dataPlot)(k, 3) = yk[1]; // lambda1 

          dlambda1 = ((*dataPlot)(k, 3)-(*dataPlot)(k-1, 3))/h; // \dot{lambda1} at 
          (*dataPlot)(k-1, 4) = -B3*dlambda1 - u2dotk; // z at K-1 !

          k++;
        }   
        dataPlot->resize(k-1, outputSize);
        std::cout<<"Reference: " << k << " steps !\n";        
        free(lambdak);
        free(wk); 
        free(lcp);
        free(q);
        return dataPlot;  
    }
    catch(siconos::exception& e)
      {
        siconos::exception::process(e);
        return NULL;
      }
    catch (...)
      {
        siconos::exception::process();
        std::cerr << "Exception caught in Computeref()" << std::endl;
        return NULL;
      }
}

double compute_error2ref(double h, SimpleMatrix* reference, SolverOptions * options, NumericsMatrix* M)
{
     try
    {
        double lambda10 = 0.0;
        double x10 = 0.01;

        size_t size = 1;    
        std::vector<double> init = {x10,lambda10};
        double * q = (double *)calloc(size,sizeof(double));

        update(q,init,h,0);

        LinearComplementarityProblem* lcp = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));
        lcp->M = M;
        lcp->q = q;
        lcp->size = size;

        size_t k = 0;
        size_t N = int((T - t0) / h) + 1;

        // ==== Simulation loop - Writing without explicit event handling =====
        k+=1;
        std::vector<double> yk = {x10,lambda10}; // previous LCP solution vector
        double * lambdak = (double *)calloc(size,sizeof(double)); // LCP solution vector
        double * wk = (double *)calloc(size,sizeof(double)); // LCP out solution vector
        int info = 0; // lcp solving flag
        double dlambda1 = 0.0; // d(lambda1)/dt at previous step  
                
        double u1k;
        double u2k;

        u2k = 1e-4*cos(10*(t0+h*0));

        // Tested variables for convergence: (x1,x2)
        SP::SiconosVector x_prev(new SiconosVector(2)); 
        x_prev->setValue(0, x10); //x10
        x_prev->setValue(1, B3*lambda10 - u2k);      //x20  
        SP::SiconosVector x(new SiconosVector(2)); 
        x->setValue(0, 0.0);
        x->setValue(1, 0.0);

        // difference in between current pwl at t approx and reference
        SP::SiconosVector x_diff(new SiconosVector(2));
        x_diff->setValue(0,0.0);
        x_diff->setValue(1,0.0);
        
        // direction of pwl approx
        SP::SiconosVector a(new SiconosVector(2));
        a->setValue(0,0.0);
        a->setValue(1,0.0);
        
        // offset of pwl approx
        SP::SiconosVector b(new SiconosVector(2));
        b->setValue(0,0.0);
        b->setValue(1,0.0);

        int j = 0;
        double error = -1.0; 
        double current_error;
        double comp_next_pwl;
        double tc = t0;
        // std::cout << "Starting Simulation \n";
        while (k<N)
        {   
            comp_next_pwl = true; // compute new pwl approx
            
            u1k = 0.301e4*sin(t0+h*(k-1));
            u2k = 1e-4*cos(10*(t0+h*k));

            // compute new discrete step    
            update(q,yk,h,(k-1));    
            info = linearComplementarity_driver(lcp, lambdak , wk,  options);
            if(info != 0) 
              {std::cout << "ERROR IN LCP SOLVING -- BREAK\n"; return 1e9;}
            // --- new solutions ---
            yk[0] = (1+A1*h)*yk[0] + B1*h*yk[1] + h*u1k;  // Fully explicit
            yk[1] = lambdak[0];

            // current point    
            (*x)(0) = yk[0]; //x1
            (*x)(1) = -B3*yk[1] - u2k;; // x2

            // compute piecewise linear solution
            (*a)(0) = ( (*x)(0)-(*x_prev)(0) ) / h;
            (*a)(1) = ( (*x)(1)-(*x_prev)(1) ) / h;
            (*b)(0) = (*x_prev)(0)-(*a)(0)*tc;
            (*b)(1) = (*x_prev)(1)-(*a)(1)*tc;

            tc = tc+h;
            k++;
          
            (*x_prev)(0) = (*x)(0);
            (*x_prev)(1) = (*x)(1);

            while(comp_next_pwl)
            {  
                (*x_diff)(0) = (*reference)(j,1) - (*a)(0)*(*reference)(j,0) - (*b)(0);
                (*x_diff)(1) = (*reference)(j,2) - (*a)(1)*(*reference)(j,0) - (*b)(1);   

                current_error = x_diff->norm2();

                if(error<current_error)
                {
                  error=current_error;   
                }
                j++;    

                if( (*reference)(j,0)>tc)
                    comp_next_pwl = false;
                if(j >= reference->size(0)) 
                    comp_next_pwl = false;       
            }
        }   
      
        std::cout << "\n (h,error) = (" << h << " , " << error << ")" << std::endl;  
        free(lambdak);
        free(wk); 
        free(lcp);
        free(q);
        return error;  
    }   
    catch(siconos::exception& e)
      {
        siconos::exception::process(e);
        return 1e9;
      }
    catch (...)
      {
        siconos::exception::process();
        std::cerr << "Exception caught in Computeref()" << std::endl;
        return 1e9;
      }
}

int main(int argc, char* argv[])
{

    double href =  5e-6;

    SolverOptions * options = solver_options_create(SICONOS_LCP_LEMKE);

    size_t size = 1; 
    NumericsMatrix* M = NM_create(NM_DENSE, size, size); // column first matrix storage
    M->matrix0[0] = B3X3B3;

    SimpleMatrix* dataSolref = compute_reference(href,options,M);
    ioMatrix::write("result_example1_ref.dat", "ascii", *dataSolref, "noDim");
    std::cout<<"Done !\n";

    std::vector<double> htable({ 1.00000000e-05, 1.08033152e-05, 1.16711619e-05, 1.26087241e-05,
                                    1.36216020e-05, 1.47158460e-05, 1.58979923e-05, 1.71751022e-05,
                                    1.85548042e-05, 2.00453398e-05, 2.16556124e-05, 2.33952406e-05,
                                    2.52746159e-05, 2.73049642e-05, 2.94984134e-05, 3.18680658e-05,
                                    3.44280759e-05, 3.71937355e-05, 4.01815648e-05, 4.34094109e-05,
                                    4.68965549e-05, 5.06638264e-05, 5.47337285e-05, 5.91305720e-05,
                                    6.38806207e-05, 6.90122480e-05, 7.45561067e-05, 8.05453121e-05,
                                    8.70156393e-05, 9.40057378e-05, 1.01557362e-04, 1.09715619e-04,
                                    1.18529241e-04, 1.28050875e-04, 1.38337396e-04, 1.49450249e-04,
                                    1.61455815e-04, 1.74425806e-04, 1.88437696e-04, 2.03575182e-04,
                                    2.19928686e-04, 2.37595891e-04, 2.56682330e-04, 2.77302012e-04,
                                    2.99578103e-04, 3.23643668e-04, 3.49642455e-04, 3.77729765e-04,
                                    4.08073370e-04, 4.40854524e-04, 4.76269038e-04, 5.14528453e-04,
                                    5.55861305e-04, 6.00514488e-04, 6.48754729e-04, 7.00870182e-04,
                                    7.57172149e-04, 8.17996938e-04, 8.83707874e-04});
                                     
    std::vector<double> errs; 
    errs.resize(htable.size());
    for(int i=0;i<htable.size();i++)
      {
        errs[i] = compute_error2ref(htable[i], dataSolref, options, M);
      } 

    SimpleMatrix res(1,htable.size());
    SiconosVector errVec(errs);
    res.setRow(0, errVec);
    ioMatrix::write("result_ex1_convergence.dat", "ascii", res, "noDim"); 

    delete(dataSolref);

    delete(options);
    delete(M);
}