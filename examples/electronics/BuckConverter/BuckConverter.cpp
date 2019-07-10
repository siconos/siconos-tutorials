#include "SiconosKernel.hpp"
#include "NumericsMatrix.h"
#include "Accelerate.h"

#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <exception>
#include <sstream>

#include "const.h"

using namespace std;
// NBHYP : at least 1
#define NBHYP 2
#define SIZEX 5
#define NSLSIZE_BUCK ((4*NBHYP) + 4)
#define SIZEZ_PAR 11
#define SIZEZ_INP 2

int main(int argc, char* argv[])
{

	double t0 = 0.0;
	double h_step = 0.05E-9;		// Time step
  double theta = 1.0;

  double T = 200E-6;          // Total simulation time

  T = 300E-6;
  bool test = false ;
  if (test)
  {
    T = 200E-8;
  }
  double VrefSettlingTime = 100e-6;
  double Rload = 10.0;

  int incx=1,incy=1,rowsize;

  int info;
  clock_t LCP_CPUtime,LCP_CPUtime_std,LCP_CPUtime_bck;

	string Modeltitle = "BuckConverter";
	double tinst;

	double VI = 3.0;	// Power supply
  double Vref = 1.8;  // output voltage setpoint

  double C  = 22E-6;
  double L = 10E-6;
//     double C  = 10E-6;
//     double L = 4E-6;

//  Ramp generator parameters
//    double RampFreq = 600E3;
//    double RampPER = 1.0/RampFreq;
  double VlowRamp = 0.0;
  double VhighRamp = 0.75*VI;
  double RampTD = 0;
  double RampTR = 1.655e-6;       // not 0 !!!
  double RampTF = 10e-9;           // not 0 !!!
  double RampPW = 1e-9;
  double RampPER = 1.6667e-6;

//  compensator parameters
  double AmpliGain = 10E3;
  double AmpliBW = 30E6;
  double R11 = 15.58E3;
  double R21 = 5.613E6;
  double C11 = 20E-12;
//     double R11 = 10E3;
//     double R21 = 8E6;
//     double C11 = 10E-12;

  double R12 = 227.8E3;

  double C21 = 1.9E-12;

  double alpha = (1.0/R11) + (1.0/R12) + (1.0/R21);
  double beta = (1.0/R11) + (1.0/R12);
  double tau11 = R11*C11;
  double tau21 = R21*C21;
  double tauAmpli = 1.0/(2.0*3.14159*AmpliBW);

//  comparator parameters
  double X1Comp = -0.1;
  double X2Comp = 0.1;
  double VsatComp = VI;
  double SlopeComp = VsatComp/(X2Comp - X1Comp);

	// Buck MOS N parameters
	double Vt0N = 2.0*VI/3.0;
	double KvalN = 10.0;
	double HalfKN = KvalN/2.0;
//    double VthDN = 0.2;     // intrinsic diode threshold
  double VthDN = 0.8;     // intrinsic diode threshold

	// Buck MOS P parameters
	double Vt0P = -2.0*VI/3.0;
	double KvalP = 10.0;
	double HalfKP = KvalP/2.0;
//    double VthDP = 0.2;     // intrinsic diode threshold
  double VthDP = 0.8;     // intrinsic diode threshold

  // // open the results file
  // ofstream outFile("BuckConverter.dat");          // checks that it's opened
  // if ( !outFile.is_open() )
  // {
  //   cout << "function write error : Fail to open \"BuckConverter.dat\"" << endl;
  //   exit(-1);
  // }

/*    ostringstream oss;
      string buffeross;*/

//    unsigned int nbPlot = 14;
  unsigned int nbPlot = 11;

	// PWL approximation of f(x) = 0 if x < 0 else f(x) = x^2
//     if (NBHYP < 2)
//     {
//         cout << " Number of hyperplanes = " << NBHYP << " , should be greater than 2 !" << endl;
//         exit(-1);
//     }
	double limhyp[NBHYP];
	double varp[NBHYP];
  double ractolpwl,widthhyp;

  if (NBHYP > 1)
  {
    ractolpwl = (VI-Vt0N)/(2.0 * (1.0 + (sqrt(2.0)*(NBHYP - 1))));
    widthhyp = 2.0*sqrt(2.0)*ractolpwl;

    cout << " Error max on fPWL for x < " << VI-Vt0N << " is : " << ractolpwl*ractolpwl << endl;;
    limhyp[0] = 0;
    varp[0] = 2.0*ractolpwl;
    limhyp[1] = (1.0 + sqrt(2.0))*ractolpwl;
    varp[1] = 2.0*widthhyp;
    for (unsigned int i = 2;i < NBHYP;i++)
    {
      limhyp[i] = limhyp[i-1] + widthhyp;
      varp[i] = 2.0*widthhyp;
    }

  } else {
    limhyp[0] = 0;
    varp[0] = 2.0*(VI-Vt0N)/(1.0 + sqrt(2.0));
  }

  double pente = varp[0];
  cout << " fmodel2_1(x) = " << pente << " * x " << endl;
  for (unsigned int i = 1;i < NBHYP;i++)
  {
    pente += varp[i];
    cout << " fmodel2_" << i+1 << "(x) = " << pente << " * x + (" << (ractolpwl*ractolpwl) - (pente*pente/4) << ")" << endl;
  }
  cout << endl;
  cout << "fmodel2(x) = x < 0      ? 0     :";
  for (unsigned int i = 1;i < NBHYP;i++) cout << "\\" << endl << "             x < " << limhyp[i] << " ? fmodel2_" << i << "(x) :";

  cout << "fmodel2_" << NBHYP << "(x)" << endl;

  cout << "Approximated conductance/resistance of power PMOS : " << HalfKP*pente << " ; " << 1.0/(HalfKP*pente) << endl;
  cout << "Approximated conductance/resistance of power NMOS : " << HalfKN*pente << " ; " << 1.0/(HalfKN*pente) << endl;
  cout << "---------------------------------------------------------------------------------------" << endl;

	// Definition of PWL solving useful vectors and matrix
	SiconosVector vec1(2*NBHYP),vec2(2*NBHYP);
	SiconosVector vecHyp(2*NBHYP);
  SimpleMatrix fPWLmat(1,2*NBHYP);

	for (unsigned int i = 0;i < NBHYP;i++)
	{
		vec1.setValue(i        ,1.0);
		vec2.setValue(i+NBHYP  ,1.0);

		vecHyp.setValue(i,        limhyp[i]);
		vecHyp.setValue(i+NBHYP,  limhyp[i]);

		fPWLmat.setValue(0,i,         varp[i]);
		fPWLmat.setValue(0,i+NBHYP,  -varp[i]);

	}


  SP::SiconosVector init_stateLS(new SiconosVector(SIZEX));

  SP::SiconosMatrix LS_A (new SimpleMatrix(SIZEX,SIZEX));
  (*LS_A)(0,0) = -( (1.0/(Rload*C)) + (beta/(C*alpha*R21)));
  (*LS_A)(0,1) = 1.0/C;
  (*LS_A)(0,2) = 1.0/(C*alpha*R21*R11);
  (*LS_A)(0,3) = -beta/(C*alpha*R21);
  (*LS_A)(0,4) =  beta/(C*alpha*R21);

  (*LS_A)(1,0) = -1.0/L;

  (*LS_A)(2,0) = (1.0 - (((1.0/R11) + (1.0/R12))/alpha)) / tau11;
  (*LS_A)(2,2) = ((1.0/(alpha*R11)) - 1.0) / tau11;
  (*LS_A)(2,3) = 1.0/(alpha*R21*tau11);
  (*LS_A)(2,4) = -((*LS_A)(2,3));

  (*LS_A)(3,0) = -((1.0/R11) + (1.0/R12))/(alpha*tau21);
  (*LS_A)(3,2) = 1.0/(alpha*R11*tau21);
  (*LS_A)(3,3) = ((1.0/(alpha*R21)) - 1.0) / tau21;
  (*LS_A)(3,4) = -((*LS_A)(3,3));

  (*LS_A)(4,0) = -((1.0/R11) + (1.0/R12))*AmpliGain / (alpha*tauAmpli);
  (*LS_A)(4,2) = AmpliGain/(alpha*R11*tauAmpli);
  (*LS_A)(4,3) = AmpliGain/(alpha*R21*tauAmpli);
  (*LS_A)(4,4) = -((AmpliGain/(alpha*R21)) + 1.0) / tauAmpli;

//     SiconosVector LS_b(SIZEX);
//     LS_b(1) = -VthDN/L;

  // --- Dynamical system creation ---
  SP::FirstOrderLinearDS LSBuckConverter (new FirstOrderLinearDS(init_stateLS,LS_A));

  SP::SiconosVector paramVin (new SiconosVector(SIZEZ_PAR + SIZEZ_INP));
	paramVin->setValue(0,VlowRamp);
	paramVin->setValue(1,VhighRamp);
	paramVin->setValue(2,RampTD);
	paramVin->setValue(3,RampTR);
	paramVin->setValue(4,RampTF);
	paramVin->setValue(5,RampPW);
	paramVin->setValue(6,RampPER);
	paramVin->setValue(7,Vref);
	paramVin->setValue(8,VrefSettlingTime);
  paramVin->setValue(9,AmpliGain/tauAmpli);
  paramVin->setValue(10,-VthDN/L);
//    SimpleMatrix* LS_T = new SimpleMatrix(SIZEX,3);
//    LS_T->setValue(4,1,AmpliGain/tauAmpli);
//    LSBuckConverter->setTPtr(LS_T);

  LSBuckConverter->setzPtr(paramVin);

  LSBuckConverter->setComputebFunction("./plugins.so","Rampb");



  SP::SimpleMatrix Coltemp  (new SimpleMatrix(2*NBHYP,1));

  SP::SimpleMatrix Int_D_buck (new SimpleMatrix(NSLSIZE_BUCK,NSLSIZE_BUCK));
  Int_D_buck->eye();

  Coltemp->setCol(0,vec1+vec2);
  Int_D_buck->setBlock(2,0, SlopeComp * *Coltemp);
  Int_D_buck->setBlock(2,1,(-SlopeComp) * *Coltemp);
  Int_D_buck->setBlock(2+(2*NBHYP),0,(-SlopeComp) * *Coltemp);
  Int_D_buck->setBlock(2+(2*NBHYP),1,SlopeComp * *Coltemp);

  Coltemp->setCol(0,vec2);
  Int_D_buck->setBlock(2, NSLSIZE_BUCK-1, (-1.0) * *Coltemp);
  Int_D_buck->setBlock(2+(2*NBHYP),NSLSIZE_BUCK-1, *Coltemp);

  Int_D_buck->setValue(NSLSIZE_BUCK-2,NSLSIZE_BUCK-2,0.);
  Int_D_buck->setValue(NSLSIZE_BUCK-2,NSLSIZE_BUCK-1,-1.0);

  Int_D_buck->setBlock(NSLSIZE_BUCK-1, 2, (-HalfKP)*fPWLmat);
  Int_D_buck->setBlock(NSLSIZE_BUCK-1, 2+(2*NBHYP),   HalfKN*fPWLmat);
  Int_D_buck->setValue(NSLSIZE_BUCK-1, NSLSIZE_BUCK-2, 1.0);
  Int_D_buck->setValue(NSLSIZE_BUCK-1, NSLSIZE_BUCK-1, 0.);

  Int_D_buck->display();
  //getchar();
  SP::SimpleMatrix Int_F0_buck(new SimpleMatrix(NSLSIZE_BUCK,SIZEZ_PAR + SIZEZ_INP));
  Int_F0_buck->setValue(0,SIZEZ_PAR,-1.0);
  Int_F0_buck->setValue(1,SIZEZ_PAR,-1.0);

  SP::SiconosVector Int_e_buck ( new SiconosVector(NSLSIZE_BUCK));
  Int_e_buck->setValue(0,X1Comp);
  Int_e_buck->setValue(1,X2Comp);
  Int_e_buck->setBlock(2,(-Vt0P) * (vec1+vec2) - VI * vec1 + VthDN * vec2 + vecHyp);
  Int_e_buck->setBlock(2+(2*NBHYP),Vt0N * (vec1+vec2) - VthDN * vec2 + vecHyp);
  Int_e_buck->setValue(NSLSIZE_BUCK-2,VI+VthDP+VthDN);

  SP::SimpleMatrix Int_C_buck ( new SimpleMatrix(NSLSIZE_BUCK,SIZEX));
  Int_C_buck->setValue(0,4,1.0);
  Int_C_buck->setValue(1,4,1.0);
  Int_C_buck->setValue(NSLSIZE_BUCK-1,1,1.0);

  SP::SimpleMatrix Int_B_buck (new SimpleMatrix(SIZEX,NSLSIZE_BUCK));
  Int_B_buck->setValue(1,NSLSIZE_BUCK-1,1.0/L);

  SP::ComplementarityConditionNSL nslaw_buck (new ComplementarityConditionNSL(NSLSIZE_BUCK));

  SP::FirstOrderLinearTIR LTIRBuckConverter_buck  (new FirstOrderLinearTIR(Int_C_buck,Int_B_buck));
  LTIRBuckConverter_buck->setDPtr(Int_D_buck);
  LTIRBuckConverter_buck->setePtr(Int_e_buck);
  LTIRBuckConverter_buck->setFPtr(Int_F0_buck);

  SP::Interaction InterBuckConverter_buck ( new Interaction(nslaw_buck,LTIRBuckConverter_buck));

  InterBuckConverter_buck->display();
  //getchar();


  // --- Model creation ---
  SP::NonSmoothDynamicalSystem NSDSBuckConverter ( new NonSmoothDynamicalSystem(t0,T));
  NSDSBuckConverter->setTitle(Modeltitle);

  // add the dynamical system in the non smooth dynamical system
  NSDSBuckConverter->insertDynamicalSystem(LSBuckConverter);
  // link the interaction and the dynamical system
  NSDSBuckConverter->link(InterBuckConverter_buck, LSBuckConverter);

  // --- End of model specification ---

// --- Simulation specification---

  SP::TimeDiscretisation TiDisc (new TimeDiscretisation(t0, h_step));

  SP::EulerMoreauOSI OSI_LSBuckConverter (new EulerMoreauOSI(theta));

  //SP::LCP LCP_BuckConverter ( new LCP("LCP","RPGS", 16, 1e-4,0, "toto",0.0,2.0 ));
  SP::LCP LCP_BuckConverter ( new LCP(SICONOS_LCP_LEMKE));
  SP::TimeStepping StratBuckConverter (new TimeStepping(NSDSBuckConverter,
                                                        TiDisc,
                                                        OSI_LSBuckConverter,LCP_BuckConverter));


//    StratBuckConverter->getEventsManagerPtr()->setTick(h_step/100000.);

//     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","NSQP", 100, 1e-7 );
//    LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","Lemke", 100, 0.0001 );
//     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","LexicoLemke", 100, 1e-4,0, "toto",0.0,1.0 );
//    LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","CPG", 1000, 1e-2 );
//     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","NLGS", 5000, 1e-5 );
//    LCP_BuckConverter->setisMSparseBlock(true);

// method* solvingMethod = (LCP_BuckConverter->getSolverPtr()->getSolvingMethodPtr());
// method* solvingMethodBackup = (LCP_BuckConverter->getSolverBackupPtr()->getSolvingMethodPtr());

// cout << "solvingMethod name = " << solvingMethod->lcp.name << endl;
// cout << "solvingMethod itermax = " << solvingMethod->lcp.itermax << endl;
// cout << "solvingMethod tol = " << solvingMethod->lcp.tol << endl;
// cout << "solvingMethod rho = " << solvingMethod->lcp.rho << endl;
// solvingMethod->lcp.chat = 0;
//	cout << "solvingMethod  chat = " << solvingMethod->lcp.chat << endl;
//cout << "solvingMethodBackup name = " << solvingMethodBackup->lcp.name << endl;

//**************************************************************************************************
//  straight implementation of integration & one step solving algorithm
//**************************************************************************************************
  double *b_straight;
  double *zpar_straight;

  b_straight = LSBuckConverter->b()->getArray();
  zpar_straight = LSBuckConverter->z()->getArray();

  double z_straight[NSLSIZE_BUCK];
  double w_straight[NSLSIZE_BUCK];

  double *fPWLmat_straight;
  fPWLmat_straight = fPWLmat.getArray();

  SiconosVector LambdaP(2*NBHYP),LambdaN(2*NBHYP);

  SP::SiconosVector xpt = LSBuckConverter->x();
  SP::SiconosVector lambdapt_buck = InterBuckConverter_buck->lambda(0);

  unsigned int k = 0;
  unsigned int N = ceil((T-t0)/h_step);

//TiDisc->NSteps(); // Number of time steps
  tinst = k*h_step;

// --- Get the values to be plotted ---
// -> saved in a matrix dataPlot
  SimpleMatrix dataPlot(N+1, nbPlot);

  SP::SiconosVector x  = LSBuckConverter->x();
  SP::SiconosVector y  = InterBuckConverter_buck->y(0);
  SP::SiconosVector lambda  = InterBuckConverter_buck->lambda(0);
  double * lambda_array = lambda->getArray();

// For the initial time step:

// time
  dataPlot(k,0) = k*h_step;

// ramp voltage
  dataPlot(k,1) = (*LSBuckConverter->z())(SIZEZ_PAR);

// MOS P drain potential
  dataPlot(k,2)  =  - VthDN;

// L current
  dataPlot(k,3)  = (*x)(1);

// output voltage
  dataPlot(k,4)  = (*x)(0);

// gate voltage Vcomp = V_G
  dataPlot(k,5)  = SlopeComp * (z_straight[0] - z_straight[1]);

// error voltage
  dataPlot(k,6)  = (*x)(4);

// DPMOS current
  dataPlot(k,7)  = z_straight[NSLSIZE_BUCK-2];

// DNMOS current
  dataPlot(k,8)  = w_straight[NSLSIZE_BUCK-1];


// // PMOS Isd current
  rowsize = 2*NBHYP;
  dataPlot(k,9)  =  HalfKP * cblas_ddot( rowsize , fPWLmat_straight , incx , &(z_straight[2]) , incy );
// *(++dataPlot) = HalfKP *cblas_ddot( &rowsize , fPWLmat_straight , &incx , &(z_straight[2]) , &incy );

// // NMOS Ids current
  dataPlot(k,10) =  HalfKN *cblas_ddot( rowsize , fPWLmat_straight , incx , &(z_straight[2+(2*NBHYP)]) , incy );
// *(++dataPlot) = HalfKN *cblas_ddot( &rowsize , fPWLmat_straight , &incx , &(z_straight[2+(2*NBHYP)]) , &incy );


// nb iterations lcp solver
//    *(++dataPlot) = nbiterlcp;

// solver backup
//    *(++dataPlot) = 0;

// nb iterations lcp solver 1eq
//    *(++dataPlot) = nbiter1eq;

// --- Compute elapsed time ---
  double t1,t2,elapsed;
  struct timeval tp;
  int rtn;
  clock_t start, endloop;
  double elapsed2;
  double elapsedCPU_LCP;


  start = clock();
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;

  cout << " Start computing loop ..." << endl;
//     LCP_BuckConverter->resetCPUtime();
//     LCP_BuckConverter->resetStat();

  try{

    LCP_CPUtime = 0;
    LCP_CPUtime_std = 0;
    LCP_CPUtime_bck = 0;
    info = 0;
    // --- Time loop  ---
    while ((k < N) && (info == 0))
    {
      k++;
      // solve ...
      StratBuckConverter->computeOneStep();
      // LSBuckConverter->display();
      // x->display();
      // getchar();
      // LCP_BuckConverter->display();

      // InterBuckConverter_buck->display();
      // getchar();
      // time
      dataPlot(k,0) = k*h_step;
      // ramp voltage
      dataPlot(k,1) = (*LSBuckConverter->z())(SIZEZ_PAR);
      // MOS P drain potential
      dataPlot(k,2)  = (*lambda)(NSLSIZE_BUCK-1) - VthDN;
      // L current
      dataPlot(k,3)  = (*x)(1);
      // output voltage
      dataPlot(k,4)  = (*x)(0);
      // gate voltage
      dataPlot(k,5)  = SlopeComp * ( (*lambda)(0) - (*lambda)(1));
      // error voltage
      dataPlot(k,6)  = (*x)(4);
      // DPMOS current
      dataPlot(k,7)  = (*lambda)(NSLSIZE_BUCK-2);
      // DNMOS current
      dataPlot(k,8)  = (*y)(NSLSIZE_BUCK-1);
      // PMOS Isd current
      rowsize = 2*NBHYP;
      dataPlot(k,9)  =  HalfKP * cblas_ddot( rowsize , fPWLmat_straight , incx , &(lambda_array[2]) , incy );
      // *(++dataPlot) = HalfKP *cblas_ddot( &rowsize , fPWLmat_straight , &incx , &(z_straight[2]) , &incy );
      //  NMOS Ids current
      dataPlot(k,10) =  HalfKN *cblas_ddot( rowsize , fPWLmat_straight , incx , &(lambda_array[2+(2*NBHYP)]) , incy );
      // *(++dataPlot) = HalfKN *cblas_ddot( &rowsize , fPWLmat_straight , &incx , &(z_straight[2+(2*NBHYP)]) , &incy );

    
      StratBuckConverter->nextStep();

    
      if ((k % (N/100)) == 0)
      {
        cerr << "-------- " << (100.0*k)/N << " % achieved... ( "<< k << " steps)" << endl;
//             cerr << "nb CPU cycles LCP = " << LCP_CPUtime << endl;
      }
    }

  } // end of "try" section
  // --- Exceptions handling ---
  catch(SiconosException e)
  {
    cout << "Time step n " << k << " : " << tinst << endl;
    cout << "SiconosException" << endl;
    cout << e.report() << endl;
  }
  catch (std::exception& e)
  {
    cout << "Exception: " << e.what() << endl;
    exit(-1);
  }
  catch(...)
  {
    cout << "Exception caught " << endl;
  }

// --- elapsed time computing ---
  endloop = clock();
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  elapsed=t2-t1;
  elapsed2 = (endloop-start) / (double)CLOCKS_PER_SEC;

  elapsedCPU_LCP = LCP_CPUtime / (double)CLOCKS_PER_SEC;

  cout << "time = " << elapsed << " --- cpu time " << elapsed2 << "--- cpu time in lcp solving : " << elapsedCPU_LCP << endl;
  cerr << "time = " << elapsed << " --- cpu time " << elapsed2 << "--- cpu time in lcp solving : " << elapsedCPU_LCP << endl;

// dataPlot (ascii) output
  ioMatrix::write("BuckConverter.dat", "ascii", dataPlot, "noDim");


  cout << "End of program" << endl;
}
