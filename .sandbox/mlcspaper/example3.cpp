#include "SiconosKernel.hpp"
#include "FirstOrderNonLinearDS.hpp"
#include "SolverOptions.h"
#include "NumericsVerbose.h"

using namespace std;

static double B31 = 1.0;
static double B32 = -2.0;

int update(double *q, vector<double> y, double h){
    try{
    q[0] = (h-1)*y[1] + (h+3)*y[2] + 2*h*y[0];  // For this example only ! NOT SAFE !
    q[1] = (2-3*h)*y[1] + (7*h-4)*y[2] - h*y[0];
    return 0;
    }
    catch (...)
    {
        cout << "unsafe memory access !\n";
        return 1;
    }
}

SimpleMatrix* compute_reference(double h, SolverOptions * options, NumericsMatrix* M)
{
  try
    {
        double T = 5.0;
        double t0 = 0.0;

        double lambda10 = 5.0;
        double lambda20 = 0.0;
        double x10 = 5.0;


        size_t size = 2;    
        vector<double> init = {x10,lambda10,lambda20};
        double * q = (double *)calloc(size,sizeof(double));

        update(q,init,h);

        LinearComplementarityProblem* lcp = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));
        lcp->M = M;
        lcp->q = q;
        lcp->size = size;

        size_t k = 0;
        size_t N = int((T - t0) / h) + 2;

        size_t outputSize = 1+5;
        SimpleMatrix* dataPlot = new SimpleMatrix(N, outputSize);  
        (*dataPlot)(0, 0) = t0;
        (*dataPlot)(0, 1) = x10; // x1
        (*dataPlot)(0, 2) = -B31*lambda10 -B32*lambda20; // x2
        (*dataPlot)(0, 3) = lambda10; // lambda1 
        (*dataPlot)(0, 4) = lambda20; // lambda2 
         // data_plot[k, 5] = z0 computed afterward as it is depending on \dot{lambda} 

        // ==== Simulation loop - Writing without explicit event handling =====
        k+=1;
        vector<double> yk = {x10,lambda10,lambda20};
        double * lambdak = (double *)calloc(size,sizeof(double));
        double * wk = (double *)calloc(size,sizeof(double));
        int info = 0;
        double dlambda1 = 0.0;
        double dlambda2 = 0.0;
        while (k<N)
        { 
            // update lcp
            update(q,yk,h);    
            info = linearComplementarity_driver(lcp, lambdak , wk,  options);
            // --- new solutions ---
            yk[0] = (1-2.0*h)*yk[0] + h*(2*yk[1]-yk[2]); // Fully explicit
            yk[1] = lambdak[0];
            yk[2] = lambdak[1];

            (*dataPlot)(k, 0) = (*dataPlot)(k-1, 0)+h;
            (*dataPlot)(k, 1) = yk[0]; // x1
            (*dataPlot)(k, 2) = -B31*yk[1] -B32*yk[2]; // x2
            (*dataPlot)(k, 3) = yk[1]; // lambda1 
            (*dataPlot)(k, 4) = yk[2]; // lambda2 
            dlambda1 = ((*dataPlot)(k, 3)-(*dataPlot)(k-1, 3))/h; // \dot{lambda1} at 
            dlambda2 = ((*dataPlot)(k, 4)-(*dataPlot)(k-1, 4))/h; // \dot{lambda2}
            (*dataPlot)(k-1, 5) = -B31*dlambda1 -B32*dlambda2; // z at K-1 !
            k++;
        }   
        dataPlot->resize(k-1, outputSize);
        cout<<"Reference: " << k << " steps !\n";        
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
        cerr << "Exception caught in Computeref()" << endl;
        return NULL;
      }
}

int compute_WeiestrassDiscretization(double h)
{
  cout << "In compute_WeiestrassDiscretization\n";
  try
    {
        // numerics_set_verbose(3);
        
        int dimX = 2;
        int dimZ = 1;
        int dimLambda = 2;

        /* Initial Conditions */
        double x10   = 5; //
        double x20   = -5; // Free choice
        double z0    = 15.0; 

        // Lambda0 are determined by the first step of Euler 
        double t0    = 0.0;     //start time
        double T     = 5;    // end time 

        // -------------------------
        // --- Dynamical systems ---
        // -------------------------

        vector<double> x0 = {x10, x20};
        SP::SiconosVector init(new SiconosVector(x0));

        /* LinearTIDS with no inputs */
        SP::SiconosMatrix A( new SimpleMatrix(dimX,dimX) ); 
        (*A)(0,0) = -2.0;
        SP::FirstOrderLinearTIDS dyn(new FirstOrderLinearTIDS(init,A));
        SP::SimpleMatrix M(new SimpleMatrix(dimX,dimX));
        (*M)(0,0) = 1.0;
        (*M)(1,1) = 1.0;

        // -----------------------------
        // --- Siconos Model Entity ---
        // ----------------------------
        SP::NonSmoothDynamicalSystem switch_dae(new NonSmoothDynamicalSystem(t0, T));

        // add the dynamical system in the non smooth dynamical system
        switch_dae->insertDynamicalSystem(dyn);

        // -------------------------
        // --- MLCP Relation ---
        // -------------------------
        SP::SimpleMatrix C( new SimpleMatrix(dimLambda+dimZ,dimX) );
        (*C)(0,1) = 0.0;  (*C)(0,1) = 1.0; 
        (*C)(1,0) = 2.0;  (*C)(1,1) = 1.0;
        (*C)(2,0) = -1.0; (*C)(2,1) = 3.0;

        SP::SimpleMatrix D( new SimpleMatrix(dimLambda+dimZ,dimLambda+dimZ) );
        (*D)(0,0) = 0.0;  (*D)(0,1) = 1.0;  (*D)(0,2) = -2.0;
        (*D)(1,0) = -1.0; (*D)(1,1) = 2.0;  (*D)(1,2) = -1.0;
        (*D)(2,0) = 2.0;  (*D)(2,1) = 0.0;  (*D)(2,2) = 1.0;

        SP::SimpleMatrix B( new SimpleMatrix(dimX,dimLambda+dimZ) );
        (*B)(0,0) = 0.0; (*B)(0,1) = 2.0; (*B)(0,2) = -1.0;
        (*B)(1,0) = 1.0; (*B)(1,1) = 0.0; (*B)(1,2) = 0.0;

        SP::FirstOrderLinearTIR relationMLCP(new FirstOrderLinearTIR(C, B) );
        relationMLCP->setDPtr(D);

        // NonSmooth law: MCP
        SP::NonSmoothLaw nslawMLCP(new MixedComplementarityConditionNSL(dimLambda,dimZ));
        // interaction 
        SP::Interaction interMLCP(new Interaction(nslawMLCP, relationMLCP));
        // link the interaction and the dynamical system
        switch_dae->link(interMLCP, dyn);

        // -----------------------------
        // --- Simulation Definition ---
        // -----------------------------

        // -- (1) OneStepIntegrators --
        double thetaOSI = 0.0;
        double gammaOSI = 1.0;
        SP::EulerMoreauOSI osi(new EulerMoreauOSI(thetaOSI,gammaOSI));


        // -- (2) Time discretisation --
        SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));

        // -- (3) one step non smooth problem
        SP::MLCP osnspb(new MLCP());
        // osnspb->setNumericsVerboseMode(true);
        
        // -- (4) Simulation setup with (1) (2) (3)
        SP::TimeStepping s(new TimeStepping(switch_dae, td, osi, osnspb));

        int N = ceil((T - t0) / h)+1; // Number of time steps

        // --- Get the values to be plotted ---
        // -> saved in a matrix dataPlot
        unsigned int outputSize = dimX+dimZ+dimLambda+1;
        SimpleMatrix dataPlot(N, outputSize);

        double lambda10 = 5.0;
        double lambda20 = 0.0;
        vector<double> lambda0 = {z0, lambda10, lambda20};
        SiconosVector initLambda = SiconosVector(lambda0); // Numerically computed
        interMLCP->setLambda(0, initLambda);
        SP::SiconosVector x = dyn->x();
        SP::SiconosVector lambdaLCP = interMLCP->lambda(0);

        dataPlot(0, 0) = switch_dae->t0();
        dataPlot(0, 1) = (*x)(0);   // x1(t)
        dataPlot(0, 2) = (*x)(1);   // x2(t)
        dataPlot(0, 3) = (*lambdaLCP)(0);   // z(t)
        dataPlot(0, 4) = (*lambdaLCP)(1);  // lambda1(t) (NS multiplier)
        dataPlot(0, 5) = (*lambdaLCP)(2);  // lambda2(t) (NS multiplier)

        // --- Time loop ---
        
        // ==== Simulation loop - Writing without explicit event handling =====
        int k = 1;

        // s->setResetAllLambda(false);
        while (s->hasNextEvent())
        {   
            s->computeOneStep(); 

            // --- Get values to be plotted ---
            dataPlot(k, 0) =  s->nextTime();
            dataPlot(k, 1) = (*x)(0);   // x1(t)
            dataPlot(k, 2) = (*x)(1);   // x2(t)
            dataPlot(k, 3) = (*lambdaLCP)(0);   // z(t)
            dataPlot(k, 4) = (*lambdaLCP)(1);  // lambda1(t) (NS multiplier)
            dataPlot(k, 5) = (*lambdaLCP)(2);  // lambda2(t) (NS multiplier)
            s->nextStep();
            k++;

        }
        cout  << "End of computation - Number of iterations done: " << k - 1 << endl;

        // --- Output files ---
        cout << "====> Output file writing ..." << endl;
        dataPlot.resize(k, outputSize);
        ioMatrix::write("result_example3_weiestrass_ref.dat", "ascii", dataPlot, "noDim");
        return 0;
    }

    catch(siconos::exception& e)
      {
         siconos::exception::process(e);
         return -1;
      }
    catch(...)
      {
        siconos::exception::process();
        cerr << "Exception caught in example3.cpp / compute_WeiestrassDiscretization" << endl;
        return -1;
      }
   
}


double compute_error2ref(double h, SimpleMatrix* reference, SolverOptions * options, NumericsMatrix* M)
{
     try
    {
        double T = 5.0;
        double t0 = 0.0;

        double lambda10 = 5.0;
        double lambda20 = 0.0;
        double x10 = 5.0;


        size_t size = 2;    
        vector<double> init = {x10,lambda10,lambda20};
        double * q = (double *)calloc(size,sizeof(double));

        update(q,init,h);

        LinearComplementarityProblem* lcp = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));
        lcp->M = M;
        lcp->q = q;
        lcp->size = size;

        size_t k = 0;
        size_t N = int((T - t0) / h) + 1;

        // ==== Simulation loop - Writing without explicit event handling =====
        k+=1;
        vector<double> yk = {x10,lambda10,lambda20}; // previous LCP solution vector
        double * lambdak = (double *)calloc(size,sizeof(double)); // LCP solution vector
        double * wk = (double *)calloc(size,sizeof(double)); // LCP out solution vector
        int info = 0; // lcp solving flag
        double dlambda1 = 0.0; // d(lambda1)/dt at previous step
        double dlambda2 = 0.0; // d(lambda1)/dt at previous step   

        // Tested variable for convergence: (x1,x2)
        SP::SiconosVector x_prev(new SiconosVector(2)); 
        x_prev->setValue(0, x10);
        x_prev->setValue(1, -B31*lambda10 -B32*lambda20);        
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

            // compute new discrete step    
            update(q,yk,h);    
            info = linearComplementarity_driver(lcp, lambdak , wk,  options);
            // --- new solutions ---
            yk[0] = (1-2.0*h)*yk[0] + h*(2*yk[1]-yk[2]); // Fully explicit
            yk[1] = lambdak[0];
            yk[2] = lambdak[1];

            // current point    
            (*x)(0) = yk[0];
            (*x)(1) = -B31*yk[1] -B32*yk[2];

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
      
        std::cout << "\n error = " << error << std::endl;  
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
        cerr << "Exception caught in Computeref()" << endl;
        return 1e9;
      }
}


double compute_error_weiestrass2ref(double h, SimpleMatrix* reference)
{
     try
    {
        
        // numerics_set_verbose(3);
        
        int dimX = 2;
        int dimZ = 1;
        int dimLambda = 2;

        /* Initial Conditions */
        double x10   = 5; //
        double x20   = -5; // Free choice
        double z0    = 15.0; 
        double t0    = 0.0;  //start time
        double T     = 5;    // end time 

        // -------------------------
        // --- Dynamical systems ---
        // ------------------------- 
        vector<double> x0 = {x10, x20};
        SP::SiconosVector init(new SiconosVector(x0));

        /* LinearTIDS with no inputs */
        SP::SiconosMatrix A( new SimpleMatrix(dimX,dimX) ); 
        (*A)(0,0) = -2.0;
        SP::FirstOrderLinearTIDS dyn(new FirstOrderLinearTIDS(init,A));
        SP::SimpleMatrix M(new SimpleMatrix(dimX,dimX));
        (*M)(0,0) = 1.0;
        (*M)(1,1) = 1.0;
        dyn->setMPtr(M);

        // -----------------------------
        // --- Siconos Model Entity ---
        // ----------------------------
        SP::NonSmoothDynamicalSystem switch_dae(new NonSmoothDynamicalSystem(t0, T));

        // add the dynamical system in the non smooth dynamical system
        switch_dae->insertDynamicalSystem(dyn);

        // -------------------------
        // --- MLCP Relation ---
        // -------------------------
        SP::SimpleMatrix C( new SimpleMatrix(dimLambda+dimZ,dimX) );
        (*C)(0,1) = 0.0;  (*C)(0,1) = 1.0; 
        (*C)(1,0) = 2.0;  (*C)(1,1) = 1.0;
        (*C)(2,0) = -1.0; (*C)(2,1) = 3.0;

        SP::SimpleMatrix D( new SimpleMatrix(dimLambda+dimZ,dimLambda+dimZ) );
        (*D)(0,0) = 0.0;  (*D)(0,1) = 1.0;  (*D)(0,2) = -2.0;
        (*D)(1,0) = -1.0; (*D)(1,1) = 2.0;  (*D)(1,2) = -1.0;
        (*D)(2,0) = 2.0;  (*D)(2,1) = 0.0;  (*D)(2,2) = 1.0;

        SP::SimpleMatrix B( new SimpleMatrix(dimX,dimLambda+dimZ) );
        (*B)(0,0) = 0.0; (*B)(0,1) = 2.0; (*B)(0,2) = -1.0;
        (*B)(1,0) = 1.0; (*B)(1,1) = 0.0; (*B)(1,2) = 0.0;

        SP::FirstOrderLinearTIR relationMLCP(new FirstOrderLinearTIR(C, B) );
        relationMLCP->setDPtr(D);

        // NonSmooth law: MCP
        SP::NonSmoothLaw nslawMLCP(new MixedComplementarityConditionNSL(dimLambda,dimZ));
        // interaction 
        SP::Interaction interMLCP(new Interaction(nslawMLCP, relationMLCP));
        // link the interaction and the dynamical system
        switch_dae->link(interMLCP, dyn);

        // -----------------------------
        // --- Simulation Definition ---
        // -----------------------------

        // -- (1) OneStepIntegrators --
        double thetaOSI = 0.0;
        double gammaOSI = 1.0;
        SP::EulerMoreauOSI osi(new EulerMoreauOSI(thetaOSI,gammaOSI));


        // -- (2) Time discretisation --
        SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));

        // -- (3) one step non smooth problem
        SP::MLCP osnspb(new MLCP());
        // osnspb->setNumericsVerboseMode(true);
        
        // -- (4) Simulation setup with (1) (2) (3)
        SP::TimeStepping s(new TimeStepping(switch_dae, td, osi, osnspb));
        size_t N = ceil((T - t0) / h)+1; // Number of time steps
        SP::SiconosVector x = dyn->x();
        SP::SiconosVector lambda = interMLCP->lambda(0);

        // ==== Simulation loop - Writing without explicit event handling =====

        SP::SiconosVector x_prev(new SiconosVector(2));
        x_prev->setValue(0,x10);
        x_prev->setValue(1,x20);

        SP::SiconosVector x_diff(new SiconosVector(2));
        x_diff->setValue(0,0.0);
        x_diff->setValue(1,0.0);
        
        SP::SiconosVector a(new SiconosVector(2));
        a->setValue(0,0.0);
        a->setValue(1,0.0);
        
        SP::SiconosVector b(new SiconosVector(2));
        b->setValue(0,0.0);
        b->setValue(1,0.0);

        int j = 0;
        double error = -1.0; 
        double current_error;
        double comp_next_pwl;
        double tc = t0;

        while (s->hasNextEvent())
        {   
            comp_next_pwl = true;
            s->computeOneStep(); 
            (*a)(0) = ( (*x)(0)-(*x_prev)(0) ) / h;
            (*a)(1) = ( (*x)(1)-(*x_prev)(1) ) / h;
            (*b)(0) = (*x_prev)(0)-(*a)(0)*tc;
            (*b)(1) = (*x_prev)(1)-(*a)(1)*tc;

            tc = s->nextTime();
            s->nextStep();
           
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
      
        std::cout << "\n error = " << error << std::endl;  
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
        cerr << "Exception caught in Computeref()" << endl;
        return 1e9;
      }
}


int main(int argc, char* argv[])
{

    double href = 0.5e-5;

    SolverOptions * options = solver_options_create(SICONOS_LCP_LEMKE);

    size_t size = 2; 
    NumericsMatrix* M = NM_create(NM_DENSE, size, size); // column first matrix storage
    M->matrix0[0] = 1; //M11
    M->matrix0[1] = -2; // M21
    M->matrix0[2] = -2; // M12
    M->matrix0[3] = 4; // M22

    SimpleMatrix* dataSolref = compute_reference(href,options,M);
    ioMatrix::write("result_example3_ref.dat", "ascii", *dataSolref, "noDim");
    cout<<"Done !\n";

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
                                    7.57172149e-04, 8.17996938e-04, 8.83707874e-04, 9.54697470e-04,
                                    1.03138977e-03, 1.11424288e-03, 1.20375170e-03, 1.30045090e-03,
                                    1.40491810e-03, 1.51777730e-03, 1.63970266e-03, 1.77142246e-03,
                                    1.91372352e-03, 2.06745584e-03, 2.23353771e-03, 2.41296118e-03,
                                    2.60679802e-03, 2.81620607e-03, 3.04243618e-03, 3.28683970e-03,
                                    3.55087652e-03, 3.83612383e-03, 4.14428548e-03, 4.47720223e-03,
                                    4.83686269e-03, 5.22541521e-03, 5.64518076e-03, 6.09866670e-03,
                                    6.58858186e-03, 7.11785265e-03, 7.68964057e-03, 8.30736107e-03,
                                    8.97470401e-03, 9.69565562e-03, 1.04745224e-02, 1.13159567e-02,
                                    1.22249846e-02, 1.32070362e-02, 1.42679775e-02, 1.54141458e-02,
                                    1.66523876e-02, 1.79900992e-02, 1.94352711e-02, 2.09965360e-02,
                                    2.26832196e-02, 2.45053971e-02, 2.64739529e-02, 2.86006458e-02,
                                    3.08981791e-02, 3.33802767e-02, 3.60617651e-02, 3.89586614e-02,
                                    4.20882699e-02, 4.54692846e-02, 4.91219013e-02, 5.30679382e-02,
                                    5.73309663e-02, 6.19364499e-02, 6.69118990e-02, 7.22870335e-02,
                                    7.80939607e-02, 8.43673672e-02, 9.11447260e-02, 9.84665203e-02,
                                    1.06376485e-01, 1.14921870e-01, 1.24153719e-01, 1.34127175e-01,
                                    1.44901815e-01, 1.56541998e-01, 1.69117254e-01, 1.82702700e-01,
                                    1.97379486e-01, 2.13235280e-01, 2.30364794e-01, 2.48870348e-01,
                                    2.68862481e-01, 2.90460612e-01, 3.13793754e-01, 3.39001283e-01,
                                    3.66233771e-01, 3.95653887e-01, 4.27437364e-01, 4.61774057e-01,
                                    4.98869069e-01, 5.38943979e-01, 5.82238167e-01, 6.29010244e-01,
                                    6.79539592e-01, 7.34128040e-01, 7.93101660e-01, 8.56812721e-01,
                                    9.25641789e-01, 1.00000000e+00});
     
    std::vector<double> errs;
    errs.resize(htable.size());
    for(int i=0;i<htable.size();i++)
    {
        errs[i] = compute_error2ref(htable[i], dataSolref, options, M);
    }

    SimpleMatrix res(1,htable.size());
    SiconosVector errVec(errs);
    res.setRow(0, errVec);
    ioMatrix::write("result_ex3_convergence.dat", "ascii", res, "noDim"); 

    int info2 = compute_WeiestrassDiscretization(href);
    cout << "terminated with code : " << info2 << endl;

    std::vector<double> errs2;
    errs2.resize(htable.size());
    for(int i=0;i<htable.size();i++)
    {
        errs2[i] = compute_error_weiestrass2ref(htable[i], dataSolref);
    }

    SimpleMatrix res2(1,htable.size());
    SiconosVector errVec2(errs2);
    res2.setRow(0, errVec2);
    ioMatrix::write("result_ex3_convergence_weiestrass.dat", "ascii", res2, "noDim");  

    delete(dataSolref);

    delete(options);
    delete(M);
}