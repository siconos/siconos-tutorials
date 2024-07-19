#include "SiconosKernel.hpp"
#include "FirstOrderNonLinearDS.hpp"
#include "SolverOptions.h"
#include "NumericsVerbose.h"

#include <math.h>

static double t0 = 0.0;
static double T  = 10.0;    // end time

double R = 1.0;    
double C =  0.1; 
double L  =  4.7e-4;


SimpleMatrix* compute_reference(double h)
{
  try
    {       
        double lambda10 = 0.0;
        double x10 = 0.01;

        size_t k = 0;
        size_t N = int((T - t0) / h) + 2;

        size_t outputSize = 1 + 4;
        SimpleMatrix* dataPlot = new SimpleMatrix(N, outputSize); 
     
        double uck;
        double u1k;
        double tk = t0;
        double ucdotk;
        uck = -(1/R)*(sin(10*tk) - 1);
        ucdotk = -(1/R)*(10*cos(10*tk));

        (*dataPlot)(0, 0) = tk;
        (*dataPlot)(0, 1) = x10; // x1
        (*dataPlot)(0, 2) = - uck; // x2 = - B3*lambda - u2(tk)
        (*dataPlot)(0, 3) = lambda10; // lambda1 
        (*dataPlot)(0, 4) = 0;  // computed afterward as it is depending on \dot{lambda} 
        k+=1; // (k in loop = k+1 in algo)

        // ==== Simulation loop - Writing without explicit event handling =====
        while (k<N)
        { 
            tk = (*dataPlot)(k-1, 0) + h;
            (*dataPlot)(k, 0) = tk;

            uck = -(1/R)*(sin(10*(tk-h)) - 1.0); // explicite (at "k")
            u1k = sin(tk-h); //explicite
            ucdotk = -(1/R)*10*cos(10*(tk-h)); //explicite
            (*dataPlot)(k, 1) = (*dataPlot)(k-1, 1) + h*((u1k/L) - uck); // x1
            (*dataPlot)(k, 3) = 0.0; // lambda1 only solution to 0 <= 1-sin(t) _|_ lambda >= 0 
            (*dataPlot)(k-1, 4) = -ucdotk - uck; // z at K-1 = -B3Ldot - B2L -v3dot -v2
            uck = -(1/R)*(sin(10*tk) - 1.0); // at "k+1"
            (*dataPlot)(k, 2) = -uck; // x2
            k++;
        }   
        dataPlot->resize(k-1, outputSize);
        std::cout<<"Reference: " << k << " steps !\n";    
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


int main(int argc, char* argv[])
{

    double href =  5e-6;

    SimpleMatrix* dataSolref = compute_reference(href);
    ioMatrix::write("result_example2_ref.dat", "ascii", *dataSolref, "noDim");
    std::cout<<"Done !\n";

}