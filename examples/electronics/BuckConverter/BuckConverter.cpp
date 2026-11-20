#include <NumericsMatrix.h>
#include <sys/time.h>

#include <SiconosKernel.hpp>
// #include "const.h"

using namespace std;
// NBHYP : at least 1
#define NBHYP 2
#define SIZEX 5
#define NSLSIZE_BUCK ((4 * NBHYP) + 4)
#define SIZEZ_PAR 11
#define SIZEZ_INP 2

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  double t0 = 0.0;
  double h_step = 0.05E-9;  // Time step
  double theta = 1.0;

  double T = 200E-6;  // Total simulation time

  T = 300E-6;
  bool test = true;
  if (test) {
    T = 200E-8;
  }
  double VrefSettlingTime = 100e-6;
  double Rload = 10.0;

  int incx = 1, incy = 1, rowsize;

  int info;
  clock_t LCP_CPUtime, LCP_CPUtime_std, LCP_CPUtime_bck;

  std::string Modeltitle = "BuckConverter";
  double tinst;

  double VI = 3.0;    // Power supply
  double Vref = 1.8;  // output voltage setpoint

  double C = 22E-6;
  double L = 10E-6;
  //     double C  = 10E-6;
  //     double L = 4E-6;

  //  Ramp generator parameters
  //    double RampFreq = 600E3;
  //    double RampPER = 1.0/RampFreq;
  double VlowRamp = 0.0;
  double VhighRamp = 0.75 * VI;
  double RampTD = 0;
  double RampTR = 1.655e-6;  // not 0 !!!
  double RampTF = 10e-9;     // not 0 !!!
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

  double alpha = (1.0 / R11) + (1.0 / R12) + (1.0 / R21);
  double beta = (1.0 / R11) + (1.0 / R12);
  double tau11 = R11 * C11;
  double tau21 = R21 * C21;
  double tauAmpli = 1.0 / (2.0 * 3.14159 * AmpliBW);

  //  comparator parameters
  double X1Comp = -0.1;
  double X2Comp = 0.1;
  double VsatComp = VI;
  double SlopeComp = VsatComp / (X2Comp - X1Comp);

  // Buck MOS N parameters
  double Vt0N = 2.0 * VI / 3.0;
  double KvalN = 10.0;
  double HalfKN = KvalN / 2.0;
  //    double VthDN = 0.2;     // intrinsic diode threshold
  double VthDN = 0.8;  // intrinsic diode threshold

  // Buck MOS P parameters
  double Vt0P = -2.0 * VI / 3.0;
  double KvalP = 10.0;
  double HalfKP = KvalP / 2.0;
  //    double VthDP = 0.2;     // intrinsic diode threshold
  double VthDP = 0.8;  // intrinsic diode threshold

  // // open the results file
  // ofstream outFile("BuckConverter.dat");          // checks that it's opened
  // if ( !outFile.is_open() )
  // {
  //   cout << "function write error : Fail to open \"BuckConverter.dat\"\n";
  //   exit(-1);
  // }

  /*    ostringstream oss;
        string buffeross;*/

  //    unsigned int nbPlot = 14;
  unsigned int nbPlot = 11;

  // PWL approximation of f(x) = 0 if x < 0 else f(x) = x^2
  //     if (NBHYP < 2)
  //     {
  //         cout << " Number of hyperplanes = " << NBHYP << " , should be greater than 2 !\n";
  //         exit(-1);
  //     }
  double limhyp[NBHYP];
  double varp[NBHYP];
  double ractolpwl, widthhyp;

  if (NBHYP > 1) {
    ractolpwl = (VI - Vt0N) / (2.0 * (1.0 + (sqrt(2.0) * (NBHYP - 1))));
    widthhyp = 2.0 * sqrt(2.0) * ractolpwl;

    std::cout << " Error max on fPWL for x < " << VI - Vt0N
              << " is : " << ractolpwl * ractolpwl << "\n";
    ;
    limhyp[0] = 0;
    varp[0] = 2.0 * ractolpwl;
    limhyp[1] = (1.0 + sqrt(2.0)) * ractolpwl;
    varp[1] = 2.0 * widthhyp;
    for (unsigned int i = 2; i < NBHYP; i++) {
      limhyp[i] = limhyp[i - 1] + widthhyp;
      varp[i] = 2.0 * widthhyp;
    }

  } else {
    limhyp[0] = 0;
    varp[0] = 2.0 * (VI - Vt0N) / (1.0 + sqrt(2.0));
  }

  double pente = varp[0];
  std::cout << " fmodel2_1(x) = " << pente << " * x \n";
  for (unsigned int i = 1; i < NBHYP; i++) {
    pente += varp[i];
    std::cout << " fmodel2_" << i + 1 << "(x) = " << pente << " * x + ("
              << (ractolpwl * ractolpwl) - (pente * pente / 4) << ")\n";
  }
  std::cout << endl;
  std::cout << "fmodel2(x) = x < 0      ? 0     :";
  for (unsigned int i = 1; i < NBHYP; i++)
    std::cout << "\\"
              << "\n"
              << "             x < " << limhyp[i] << " ? fmodel2_" << i << "(x) :";

  std::cout << "fmodel2_" << NBHYP << "(x)\n";

  std::cout << "Approximated conductance/resistance of power PMOS : " << HalfKP * pente
            << " ; " << 1.0 / (HalfKP * pente) << endl;
  std::cout << "Approximated conductance/resistance of power NMOS : " << HalfKN * pente
            << " ; " << 1.0 / (HalfKN * pente) << endl;
  std::cout
      << "-----------------------------------------------------------------------------------"
         "----\n";

  // Definition of PWL solving useful vectors and matrix
  Vector vec1(2 * NBHYP), vec2(2 * NBHYP);
  Vector vecHyp(2 * NBHYP);
  Matrix fPWLmat(1, 2 * NBHYP);

  for (unsigned int i = 0; i < NBHYP; i++) {
    vec1(i) = 1.0;
    vec2(i + NBHYP) = 1.0;

    vecHyp(i) = limhyp[i];
    vecHyp(i + NBHYP) = limhyp[i];

    fPWLmat.setValue(0, i, varp[i]);
    fPWLmat.setValue(0, i + NBHYP, -varp[i]);
  }

  // --- Dynamical system creation ---
  auto init_stateLS = std::make_shared<Vector>(SIZEX);
  auto LSBuckConverter = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
      *init_stateLS, siconos::algebra::alias_t);

  Matrix LS_A{SIZEX, SIZEX};
  LS_A.setZero();
  LS_A(0, 0) = -((1.0 / (Rload * C)) + (beta / (C * alpha * R21)));
  LS_A(0, 1) = 1.0 / C;
  LS_A(0, 2) = 1.0 / (C * alpha * R21 * R11);
  LS_A(0, 3) = -beta / (C * alpha * R21);
  LS_A(0, 4) = beta / (C * alpha * R21);

  LS_A(1, 0) = -1.0 / L;

  LS_A(2, 0) = (1.0 - (((1.0 / R11) + (1.0 / R12)) / alpha)) / tau11;
  LS_A(2, 2) = ((1.0 / (alpha * R11)) - 1.0) / tau11;
  LS_A(2, 3) = 1.0 / (alpha * R21 * tau11);
  LS_A(2, 4) = -(LS_A(2, 3));

  LS_A(3, 0) = -((1.0 / R11) + (1.0 / R12)) / (alpha * tau21);
  LS_A(3, 2) = 1.0 / (alpha * R11 * tau21);
  LS_A(3, 3) = ((1.0 / (alpha * R21)) - 1.0) / tau21;
  LS_A(3, 4) = -(LS_A(3, 3));

  LS_A(4, 0) = -((1.0 / R11) + (1.0 / R12)) * AmpliGain / (alpha * tauAmpli);
  LS_A(4, 2) = AmpliGain / (alpha * R11 * tauAmpli);
  LS_A(4, 3) = AmpliGain / (alpha * R21 * tauAmpli);
  LS_A(4, 4) = -((AmpliGain / (alpha * R21)) + 1.0) / tauAmpli;

  LSBuckConverter->setConstantA(LS_A, siconos::algebra::alias_t);

  //     SiconosVector LS_b(SIZEX);
  //     LS_b(1) = -VthDN/L;

  Vector paramVin{SIZEZ_PAR + SIZEZ_INP};

  paramVin(0) = VlowRamp;
  paramVin(1) = VhighRamp;
  paramVin(2) = RampTD;
  paramVin(3) = RampTR;
  paramVin(4) = RampTF;
  paramVin(5) = RampPW;
  paramVin(6) = RampPER;
  paramVin(7) = Vref;
  paramVin(8) = VrefSettlingTime;
  paramVin(9) = AmpliGain / tauAmpli;
  paramVin(10) = -VthDN / L;
  //    SiconosMatrix* LS_T = new SiconosMatrix(SIZEX,3);
  //    LS_T->setValue(4,1,AmpliGain/tauAmpli);
  //    LSBuckConverter->setTPtr(LS_T);

  LSBuckConverter->setComputebVectorFunction(
      [&paramVin](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        // Warning: paramVin is modified at each call of b(t)!

        // double epsitime = 1e-15;
        double VlowRamp = paramVin[0];
        double VhighRamp = paramVin[1];
        double RampTD = paramVin[2];
        double RampTR = paramVin[3];
        double RampTF = paramVin[4];
        double RampPW = paramVin[5];
        double RampPER = paramVin[6];
        double SlopeRiseRamp = (VhighRamp - VlowRamp) / RampTR;
        double SlopeFallRamp = (VhighRamp - VlowRamp) / RampTF;
        double phaseRamp;

        double Vref = paramVin[7];
        double VrefSettlingTime = paramVin[8];

        double UPtr[2];

        phaseRamp = std::fmod(time, RampPER);

        if (phaseRamp < RampTD)
          UPtr[0] = VlowRamp;
        else if (phaseRamp < (RampTD + RampTR))
          UPtr[0] = VlowRamp + (SlopeRiseRamp * (phaseRamp - RampTD));
        else if (phaseRamp < (RampTD + RampTR + RampPW))
          UPtr[0] = VhighRamp;
        else if (phaseRamp < (RampTD + RampTR + RampPW + RampTF))
          UPtr[0] = VhighRamp - (SlopeFallRamp * (phaseRamp - (RampTD + RampTR + RampPW)));
        else
          UPtr[0] = VlowRamp;

        /*if (time < VrefSettlingTime) UPtr[1] = Vref*time/VrefSettlingTime;
          else UPtr[1] = Vref; */
        if (time < VrefSettlingTime)
          UPtr[1] = Vref * time / VrefSettlingTime;
        else if (time > 2. * VrefSettlingTime)
          UPtr[1] = -Vref * time / VrefSettlingTime + 3 * Vref;
        else
          UPtr[1] = Vref;

        paramVin[paramVin.size() - 2] = UPtr[0];
        paramVin[paramVin.size() - 1] = UPtr[1];
        result.setZero();
        result[1] = paramVin[10];
        result[4] = paramVin[9] * UPtr[1];
      });

  Matrix int_D_buck{NSLSIZE_BUCK, NSLSIZE_BUCK};
  int_D_buck.setIdentity();

  auto Coltemp = vec1 + vec2;
  auto colsize = vec2.size();
  int_D_buck.block(2, 0, colsize, 1) = SlopeComp * Coltemp;
  int_D_buck.block(2, 1, colsize, 1) = (-SlopeComp) * Coltemp;
  int_D_buck.block(2 + (2 * NBHYP), 0, colsize, 1) = (-SlopeComp) * Coltemp;
  int_D_buck.block(2 + (2 * NBHYP), 1, colsize, 1) = SlopeComp * Coltemp;

  int_D_buck.block(2, NSLSIZE_BUCK - 1, colsize, 1) = (-1.0) * vec2;
  int_D_buck.block(2 + (2 * NBHYP), NSLSIZE_BUCK - 1, colsize, 1) = vec2;

  int_D_buck.setValue(NSLSIZE_BUCK - 2, NSLSIZE_BUCK - 2, 0.);
  int_D_buck.setValue(NSLSIZE_BUCK - 2, NSLSIZE_BUCK - 1, -1.0);

  int_D_buck.block(NSLSIZE_BUCK - 1, 2, fPWLmat.rows(), fPWLmat.cols()) = (-HalfKP) * fPWLmat;
  int_D_buck.block(NSLSIZE_BUCK - 1, 2 + (2 * NBHYP), fPWLmat.rows(), fPWLmat.cols()) =
      HalfKN * fPWLmat;
  int_D_buck.setValue(NSLSIZE_BUCK - 1, NSLSIZE_BUCK - 2, 1.);
  int_D_buck.setValue(NSLSIZE_BUCK - 1, NSLSIZE_BUCK - 1, 0.);

  //  siconos::algebra::print(int_D_buck);
  // getchar();

  Matrix int_F0_buck{NSLSIZE_BUCK, SIZEZ_PAR + SIZEZ_INP};
  int_F0_buck.setZero();
  int_F0_buck.setValue(0, SIZEZ_PAR, -1.0);
  int_F0_buck.setValue(1, SIZEZ_PAR, -1.0);

  Vector int_e_buck{NSLSIZE_BUCK};
  int_e_buck.setZero();
  int_e_buck(0) = X1Comp;
  int_e_buck(1) = X2Comp;
  int_e_buck.segment(2, vecHyp.size()) =
      (-Vt0P) * (vec1 + vec2) - VI * vec1 + VthDN * vec2 + vecHyp;
  int_e_buck.segment(2 + (2 * NBHYP), vecHyp.size()) =
      Vt0N * (vec1 + vec2) - VthDN * vec2 + vecHyp;
  int_e_buck(NSLSIZE_BUCK - 2) = VI + VthDP + VthDN;

  Matrix int_C_buck{NSLSIZE_BUCK, SIZEX};
  int_C_buck.setZero();
  int_C_buck.setValue(0, 4, 1.0);
  int_C_buck.setValue(1, 4, 1.0);
  int_C_buck.setValue(NSLSIZE_BUCK - 1, 1, 1.0);

  Matrix int_B_buck{SIZEX, NSLSIZE_BUCK};
  int_B_buck.setZero();
  int_B_buck.setValue(1, NSLSIZE_BUCK - 1, 1.0 / L);

  auto nslaw_buck =
      std::make_shared<siconos::modeling::ComplementarityConditionNSL>(NSLSIZE_BUCK);

  auto LTIRBuckConverter_buck = std::make_shared<siconos::modeling::FirstOrderLinearR>();
  LTIRBuckConverter_buck->setConstantB(int_B_buck);
  LTIRBuckConverter_buck->setConstantC(int_C_buck);
  LTIRBuckConverter_buck->setConstantD(int_D_buck);
  LTIRBuckConverter_buck->setComputeeVectorFunction(
      [&int_e_buck, &int_F0_buck, &paramVin](
          double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        // We must take into account paramVin values updated by the DS (call b(t))
        result = int_e_buck + int_F0_buck * paramVin;
      });

  auto InterBuckConverter_buck =
      std::make_shared<siconos::modeling::Interaction>(nslaw_buck, LTIRBuckConverter_buck);

  //   siconos::algebra::print(*InterBuckConverter_buck);
  // getchar();

  // --- Model creation ---
  auto NSDSBuckConverter =
      std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
  NSDSBuckConverter->setTitle(Modeltitle);

  // add the dynamical system in the non smooth dynamical system
  NSDSBuckConverter->insertDynamicalSystem(LSBuckConverter);
  // link the interaction and the dynamical system
  NSDSBuckConverter->link(InterBuckConverter_buck, LSBuckConverter);

  // --- End of model specification ---

  // --- Simulation specification---

  auto TiDisc = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h_step);

  auto OSI_LSBuckConverter = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);

  // auto LCP_BuckConverter ( new LCP("LCP","RPGS", 16, 1e-4,0, "toto",0.0,2.0 ));
  auto LCP_BuckConverter =
      std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_LEMKE);
  auto StratBuckConverter = std::make_shared<siconos::simulation::TimeStepping>(
      NSDSBuckConverter, TiDisc, OSI_LSBuckConverter, LCP_BuckConverter);

  //    StratBuckConverter->getEventsManagerPtr()->setTick(h_step/100000.);

  //     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","NSQP", 100, 1e-7 );
  //    LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","Lemke", 100, 0.0001 );
  //     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","LexicoLemke", 100, 1e-4,0,
  //     "toto",0.0,1.0 );
  //    LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","CPG", 1000, 1e-2 );
  //     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","NLGS", 5000, 1e-5 );
  //    LCP_BuckConverter->setisMSparseBlock(true);

  // method* solvingMethod = (LCP_BuckConverter->getSolverPtr()->getSolvingMethodPtr());
  // method* solvingMethodBackup =
  // (LCP_BuckConverter->getSolverBackupPtr()->getSolvingMethodPtr());

  // std::cout << "solvingMethod name = " << solvingMethod->lcp.name << endl;
  // std::cout << "solvingMethod itermax = " << solvingMethod->lcp.itermax << endl;
  // std::cout << "solvingMethod tol = " << solvingMethod->lcp.tol << endl;
  // std::cout << "solvingMethod rho = " << solvingMethod->lcp.rho << endl;
  // solvingMethod->lcp.chat = 0;
  //	std::cout << "solvingMethod  chat = " << solvingMethod->lcp.chat << endl;
  // std::cout << "solvingMethodBackup name = " << solvingMethodBackup->lcp.name << endl;

  //**************************************************************************************************
  //  straight implementation of integration & one step solving algorithm
  //**************************************************************************************************

  Vector w_straight{NSLSIZE_BUCK};
  w_straight.setZero();
  Vector z_straight{NSLSIZE_BUCK};
  z_straight.setZero();

  double* fPWLmat_straight;
  fPWLmat_straight = fPWLmat.data();

  Vector LambdaP(2 * NBHYP), LambdaN(2 * NBHYP);

  auto xpt = LSBuckConverter->x();
  auto lambdapt_buck = InterBuckConverter_buck->lambda(0);

  unsigned int k = 0;
  unsigned int N = ceil((T - t0) / h_step);

  // TiDisc->NSteps(); // Number of time steps
  tinst = k * h_step;

  // --- Get the values to be plotted ---
  // -> saved in a matrix dataPlot
  Matrix dataPlot(N + 1, nbPlot);

  auto x = LSBuckConverter->x();
  auto y = InterBuckConverter_buck->y(0);
  auto lambda = InterBuckConverter_buck->lambda(0);

  // For the initial time step:

  // time
  dataPlot(k, 0) = k * h_step;

  // ramp voltage
  dataPlot(k, 1) = paramVin(SIZEZ_PAR);

  // MOS P drain potential
  dataPlot(k, 2) = -VthDN;

  // L current
  dataPlot(k, 3) = (*x)(1);

  // output voltage
  dataPlot(k, 4) = (*x)(0);

  // gate voltage Vcomp = V_G
  dataPlot(k, 5) = SlopeComp * (z_straight(0) - z_straight(1));

  // error voltage
  dataPlot(k, 6) = (*x)(4);

  // DPMOS current
  dataPlot(k, 7) = z_straight(NSLSIZE_BUCK - 2);

  // DNMOS current
  dataPlot(k, 8) = w_straight(NSLSIZE_BUCK - 1);

  // // PMOS Isd current
  rowsize = 2 * NBHYP;
  // dataPlot(k, 9) =
  //     HalfKP * cblas_ddot(rowsize, fPWLmat_straight, incx, &(z_straight[2]), incy);
  // *(++dataPlot) = HalfKP *cblas_ddot( &rowsize , fPWLmat_straight , &incx , &(z_straight[2])
  // , &incy );
  dataPlot(k, 9) = HalfKP * fPWLmat.col(0).head(rowsize).dot(z_straight.segment(2, rowsize));

  // // NMOS Ids current
  // dataPlot(k, 10) = HalfKN * cblas_ddot(rowsize, fPWLmat_straight, incx,
  //                                       &(z_straight[2 + (2 * NBHYP)]), incy);
  dataPlot(k, 10) =
      HalfKN * fPWLmat.col(0).head(rowsize).dot(z_straight.segment(2 + (2 * NBHYP), rowsize));

  // *(++dataPlot) = HalfKN *cblas_ddot( &rowsize , fPWLmat_straight , &incx ,
  // &(z_straight[2+(2*NBHYP)]) , &incy );

  // nb iterations lcp solver
  //    *(++dataPlot) = nbiterlcp;

  // solver backup
  //    *(++dataPlot) = 0;

  // nb iterations lcp solver 1eq
  //    *(++dataPlot) = nbiter1eq;

  // --- Compute elapsed time ---
  double t1, t2, elapsed;
  struct timeval tp;
  int rtn;
  clock_t start, endloop;
  double elapsed2;
  double elapsedCPU_LCP;

  start = clock();
  rtn = gettimeofday(&tp, NULL);
  t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

  std::cout << " Start computing loop ...\n";
  //     LCP_BuckConverter->resetCPUtime();
  //     LCP_BuckConverter->resetStat();

  try {
    LCP_CPUtime = 0;
    LCP_CPUtime_std = 0;
    LCP_CPUtime_bck = 0;
    info = 0;
    // --- Time loop  ---
    while ((k < N) && (info == 0)) {
      k++;
      // solve ...
      StratBuckConverter->computeOneStep();
      // siconos::algebra::print(*LSBuckConverter);
      // siconos::algebra::print(*x);
      // getchar();
      // siconos::algebra::print(*LCP_BuckConverter);

      // siconos::algebra::print(*InterBuckConverter_buck);
      // getchar();
      // time
      dataPlot(k, 0) = k * h_step;
      // ramp voltage
      dataPlot(k, 1) = (paramVin)(SIZEZ_PAR);
      // MOS P drain potential
      dataPlot(k, 2) = (*lambda)(NSLSIZE_BUCK - 1) - VthDN;
      // L current
      dataPlot(k, 3) = (*x)(1);
      // output voltage
      dataPlot(k, 4) = (*x)(0);
      // gate voltage
      dataPlot(k, 5) = SlopeComp * ((*lambda)(0) - (*lambda)(1));
      // error voltage
      dataPlot(k, 6) = (*x)(4);
      // DPMOS current
      dataPlot(k, 7) = (*lambda)(NSLSIZE_BUCK - 2);
      // DNMOS current
      dataPlot(k, 8) = (*y)(NSLSIZE_BUCK - 1);
      // PMOS Isd current
      rowsize = 2 * NBHYP;
      // dataPlot(k, 9) =
      //     HalfKP * cblas_ddot(rowsize, fPWLmat_straight, incx, &(lambda_array[2]), incy);
      // *(++dataPlot) = HalfKP *cblas_ddot( &rowsize , fPWLmat_straight , &incx ,
      // &(z_straight[2]) , &incy );
      dataPlot(k, 9) = HalfKP * fPWLmat.col(0).head(rowsize).dot(lambda->segment(2, rowsize));

      //  NMOS Ids current
      // dataPlot(k, 10) = HalfKN * cblas_ddot(rowsize, fPWLmat_straight, incx,
      //                                       &(lambda_array[2 + (2 * NBHYP)]), incy);
      dataPlot(k, 10) =
          HalfKN * fPWLmat.col(0).head(rowsize).dot(lambda->segment(2 + (2 * NBHYP), rowsize));

      // *(++dataPlot) = HalfKN *cblas_ddot( &rowsize , fPWLmat_straight , &incx ,
      // &(z_straight[2+(2*NBHYP)]) , &incy );

      StratBuckConverter->nextStep();

      if ((k % (N / 100)) == 0) {
        std::cerr << "-------- " << (100.0 * k) / N << " % achieved... ( " << k << " steps)\n";
        //             cerr << "nb CPU cycles LCP = " << LCP_CPUtime << endl;
      }
    }

  }  // end of "try" section
  // --- Exceptions handling ---
  catch (...) {
    siconos::exception::process();
    return 1;
  }

  // --- elapsed time computing ---
  endloop = clock();
  rtn = gettimeofday(&tp, NULL);
  t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
  elapsed = t2 - t1;
  elapsed2 = (endloop - start) / (double)CLOCKS_PER_SEC;

  elapsedCPU_LCP = LCP_CPUtime / (double)CLOCKS_PER_SEC;

  std::cout << "time = " << elapsed << " --- cpu time " << elapsed2
            << "--- cpu time in lcp solving : " << elapsedCPU_LCP << endl;
  // dataPlot (ascii) output
  siconos::algebra::io::write("BuckConverter.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);

  // relative tolerance: the event-driven integration (Lsodar, event location) differs
  // with the platform (libm, blas) by up to 1e-9 in absolute value
  double error = 0.0, eps = 1e-8 * dataPlot.cwiseAbs().maxCoeff();
  if ((error = siconos::algebra::io::compareRefFile(dataPlot, "BuckConverter.ref", eps)) > eps)
    return 1;

  std::cout << "End of program\n";
  return 0;
}
