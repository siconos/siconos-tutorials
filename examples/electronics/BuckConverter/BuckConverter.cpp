#include "SiconosKernel.h"

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
    double VrefSettlingTime = 100e-6;
    double Rload = 10.0;

    char NOTRANS = 'N';
    int incx,incy,rowsize,colsize;

    double a1,b1;

    int info;
    int iter,titer;
    double err;
    clock_t startLCPsolve,startLCPuni,timeLCPuni;
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

   // open the results file
    ofstream outFile("BuckConverter.dat");          // checks that it's opened
    if ( !outFile.is_open() )
    {
        cout << "function write error : Fail to open \"BuckConverter.dat\"" << endl;
        exit(-1);
    }

/*    ostringstream oss;
    string buffeross;*/

//    unsigned int nbPlot = 14;
    unsigned int nbPlot = 11;
    double *dataPlot,*dataPlotini;
    char buffer[21];
    int nblignebuf = 40;
    char bufferligne[nblignebuf*((nbPlot*20)+2)];

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
	SimpleVector vec1(2*NBHYP),vec2(2*NBHYP);
	SimpleVector vecHyp(2*NBHYP);
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


    SimpleVector init_stateLS(SIZEX);

    SimpleMatrix LS_A(SIZEX,SIZEX);
    LS_A(0,0) = -( (1.0/(Rload*C)) + (beta/(C*alpha*R21)));
    LS_A(0,1) = 1.0/C;
    LS_A(0,2) = 1.0/(C*alpha*R21*R11);
    LS_A(0,3) = -beta/(C*alpha*R21);
    LS_A(0,4) =  beta/(C*alpha*R21);

    LS_A(1,0) = -1.0/L;

    LS_A(2,0) = (1.0 - (((1.0/R11) + (1.0/R12))/alpha)) / tau11;
    LS_A(2,2) = ((1.0/(alpha*R11)) - 1.0) / tau11;
    LS_A(2,3) = 1.0/(alpha*R21*tau11);
    LS_A(2,4) = -(LS_A(2,3));

    LS_A(3,0) = -((1.0/R11) + (1.0/R12))/(alpha*tau21);
    LS_A(3,2) = 1.0/(alpha*R11*tau21);
    LS_A(3,3) = ((1.0/(alpha*R21)) - 1.0) / tau21;
    LS_A(3,4) = -(LS_A(3,3));

    LS_A(4,0) = -((1.0/R11) + (1.0/R12))*AmpliGain / (alpha*tauAmpli);
    LS_A(4,2) = AmpliGain/(alpha*R11*tauAmpli);
    LS_A(4,3) = AmpliGain/(alpha*R21*tauAmpli);
    LS_A(4,4) = -((AmpliGain/(alpha*R21)) + 1.0) / tauAmpli;

//     SimpleVector LS_b(SIZEX);
//     LS_b(1) = -VthDN/L;

    // --- Dynamical system creation ---
    FirstOrderLinearDS* LSBuckConverter = new FirstOrderLinearDS(1,init_stateLS,LS_A);

	SimpleVector *paramVin = new SimpleVector(SIZEZ_PAR + SIZEZ_INP);
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

    LSBuckConverter->setZPtr(paramVin);

    LSBuckConverter->setComputeBFunction("./RampbPlugin.so","Rampb");


    DynamicalSystemsSet allDS;
    allDS.insert(LSBuckConverter);

    SimpleMatrix* Coltemp = new SimpleMatrix(2*NBHYP,1);

    SiconosMatrix* Int_D_buck = new SimpleMatrix(NSLSIZE_BUCK,NSLSIZE_BUCK);
    Int_D_buck->eye();

    Coltemp->setCol(0,vec1+vec2);
    Int_D_buck->matrixCopy(SlopeComp * *Coltemp,2,0);
    Int_D_buck->matrixCopy((-SlopeComp) * *Coltemp,2,1);
    Int_D_buck->matrixCopy((-SlopeComp) * *Coltemp,2+(2*NBHYP),0);
    Int_D_buck->matrixCopy(SlopeComp * *Coltemp,2+(2*NBHYP),1);

    Coltemp->setCol(0,vec2);
    Int_D_buck->matrixCopy((-1.0) * *Coltemp, 2,           NSLSIZE_BUCK-1);
    Int_D_buck->matrixCopy(*Coltemp,          2+(2*NBHYP), NSLSIZE_BUCK-1);

    Int_D_buck->setValue(NSLSIZE_BUCK-2,NSLSIZE_BUCK-2,0.);
    Int_D_buck->setValue(NSLSIZE_BUCK-2,NSLSIZE_BUCK-1,-1.0);

    Int_D_buck->matrixCopy((-HalfKP)*fPWLmat, NSLSIZE_BUCK-1, 2);
    Int_D_buck->matrixCopy(   HalfKN*fPWLmat, NSLSIZE_BUCK-1, 2+(2*NBHYP));
    Int_D_buck->setValue(NSLSIZE_BUCK-1, NSLSIZE_BUCK-2, 1.0);
    Int_D_buck->setValue(NSLSIZE_BUCK-1, NSLSIZE_BUCK-1, 0.);

    SiconosMatrix* Int_F0_buck = new SimpleMatrix(NSLSIZE_BUCK,SIZEZ_PAR + SIZEZ_INP);
    Int_F0_buck->setValue(0,SIZEZ_PAR,-1.0);
    Int_F0_buck->setValue(1,SIZEZ_PAR,-1.0);

    SimpleVector* Int_e_buck = new SimpleVector(NSLSIZE_BUCK);
    Int_e_buck->setValue(0,X1Comp);
    Int_e_buck->setValue(1,X2Comp);
    Int_e_buck->setBlock(2,(-Vt0P) * (vec1+vec2) - VI * vec1 + VthDN * vec2 + vecHyp);
    Int_e_buck->setBlock(2+(2*NBHYP),Vt0N * (vec1+vec2) - VthDN * vec2 + vecHyp);
    Int_e_buck->setValue(NSLSIZE_BUCK-2,VI+VthDP+VthDN);

    SiconosMatrix* Int_C_buck = new SimpleMatrix(NSLSIZE_BUCK,SIZEX);
    Int_C_buck->setValue(0,4,1.0);
    Int_C_buck->setValue(1,4,1.0);
    Int_C_buck->setValue(NSLSIZE_BUCK-1,1,1.0);

    SiconosMatrix* Int_B_buck = new SimpleMatrix(SIZEX,NSLSIZE_BUCK);
    Int_B_buck->setValue(1,NSLSIZE_BUCK-1,1.0/L);

    ComplementarityConditionNSL *nslaw_buck = new ComplementarityConditionNSL(NSLSIZE_BUCK);

    FirstOrderLinearTIR* LTIRBuckConverter_buck = new FirstOrderLinearTIR(Int_C_buck,Int_B_buck);
    LTIRBuckConverter_buck->setDPtr(Int_D_buck);
    LTIRBuckConverter_buck->setEPtr(Int_e_buck);
    LTIRBuckConverter_buck->setFPtr(Int_F0_buck);

    Interaction* InterBuckConverter_buck = new Interaction("InterBuckConverter_buck",allDS,0,NSLSIZE_BUCK,nslaw_buck,LTIRBuckConverter_buck);

    InteractionsSet allInter;
    allInter.insert(InterBuckConverter_buck);

    NonSmoothDynamicalSystem* NSDSBuckConverter = new NonSmoothDynamicalSystem(allDS,allInter,false);

    // --- Model creation ---
    Model BuckConverter(t0,T,Modeltitle);
    BuckConverter.setNonSmoothDynamicalSystemPtr(NSDSBuckConverter);

//	cout << " Matrice LCP : " << endl;
//	SimpleMatrix MatLCP(intersize,intersize);
//	MatLCP = *Int_D + (h_step * theta * (*Int_C * *Int_B));

/*
    ofstream outFileMat("matriceD.sce");          // checks that it's opened
    if ( !outFileMat.is_open() )
         SiconosMatrixException::selfThrow("function write error : Fail to open \"matriceD.sce\"");
//
    outFileMat << "matriceD = [.." << endl;

	for (unsigned int i = 0;i < NSLSIZE;i++)
	{
        unsigned int j;
		for (j = 0;j < NSLSIZE-1;j++)
        {
			outFileMat << setprecision(5) << (*Int_D)(i,j) << " , ";
        }
		outFileMat << (*Int_D)(i,j);
        if (i != NSLSIZE-1) outFileMat << " ;" << endl;
        else outFileMat << " ]" << endl;
	}
    outFileMat.close();
*/
	// --- End of model specification ---

    // --- Simulation specification---

    TimeDiscretisation* TiDisc = new TimeDiscretisation(h_step,&BuckConverter);

    TimeStepping* StratBuckConverter = new TimeStepping(TiDisc);
//    StratBuckConverter->getEventsManagerPtr()->setTick(h_step/100000.);

    Moreau* OSI_LSBuckConverter = new Moreau(LSBuckConverter,theta,StratBuckConverter);

//     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","NSQP", 100, 1e-7 );
//    LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","Lemke", 100, 0.0001 );
     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","RPGS", 16, 1e-4,0, "toto",0.0,2.0 );
//     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","LexicoLemke", 100, 1e-4,0, "toto",0.0,1.0 );
//    LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","CPG", 1000, 1e-2 );
//     LCP* LCP_BuckConverter = new LCP(StratBuckConverter,"LCP","NLGS", 5000, 1e-5 );
//    LCP_BuckConverter->setisMSparseBlock(true);

    method* solvingMethod = (LCP_BuckConverter->getSolverPtr()->getSolvingMethodPtr());
    method* solvingMethodBackup = (LCP_BuckConverter->getSolverBackupPtr()->getSolvingMethodPtr());

	cout << "solvingMethod name = " << solvingMethod->lcp.name << endl;
	cout << "solvingMethod itermax = " << solvingMethod->lcp.itermax << endl;
	cout << "solvingMethod tol = " << solvingMethod->lcp.tol << endl;
	cout << "solvingMethod rho = " << solvingMethod->lcp.rho << endl;
	solvingMethod->lcp.chat = 0;
//	cout << "solvingMethod  chat = " << solvingMethod->lcp.chat << endl;
    cout << "solvingMethodBackup name = " << solvingMethodBackup->lcp.name << endl;

    cout << " -----  Model description ------" << endl;
    BuckConverter.display();

    StratBuckConverter->initialize();

	cout << "End strategy initialize " << endl;

//**************************************************************************************************
//  straight implementation of integration & one step solving algorithm
//**************************************************************************************************
    double x_straight[SIZEX];
    double xtemp_straight[SIZEX];
    double xfree_straight[SIZEX];
    double *b_straight;
    double *zpar_straight;
    double r_straight[SIZEX];

    b_straight = LSBuckConverter->getBPtr()->getArray();
    zpar_straight = LSBuckConverter->getZPtr()->getArray();

    double q_buck_straight[NSLSIZE_BUCK];

    double q_straight[NSLSIZE_BUCK];
    double z_straight[NSLSIZE_BUCK];
    double w_straight[NSLSIZE_BUCK];
    double zprev_straight[NSLSIZE_BUCK];
    double wprev_straight[NSLSIZE_BUCK];

    double *M_straight;
    M_straight = LCP_BuckConverter->getMPtr()->getArray();

    OSI_LSBuckConverter->getWPtr(LSBuckConverter)->PLUInverseInPlace();
    double *W_straight;
    W_straight = OSI_LSBuckConverter->getWPtr(LSBuckConverter)->getArray();

    double *B_buck_straight;
    B_buck_straight = Int_B_buck->getArray();
    double *C_buck_straight;
    C_buck_straight = Int_C_buck->getArray();
    double *F0_buck_straight;
    F0_buck_straight = Int_F0_buck->getArray();
    double *e_buck_straight;
    e_buck_straight = Int_e_buck->getArray();

    for (unsigned int i = 0;i < SIZEX;i++) x_straight[i] = 0;
    for (unsigned int i = 0;i < NSLSIZE_BUCK;i++) z_straight[i] = 0;

    double *fPWLmat_straight;
    fPWLmat_straight = fPWLmat.getArray();

    incx = 1;
    incy = 1;

/*    int indic[7] = {0 , 1 , 6 , 7 , 8 , 9 , 11};
    int indicop[5] = {2 , 3 , 4 , 5 , 10};
    int ipiv[NSLSIZE_BUCK];
    double submatlcp[NSLSIZE_BUCK*NSLSIZE_BUCK];
    double subq[NSLSIZE_BUCK];
    double qexa[NSLSIZE_BUCK] = {-0.203955 , -0.00395453 , -1 , -0.5 , 2.8 , 3.3 , 2 , 2.5 , 1.2 , 1.7 , 4.6 , 0.278073};
    double zexa[NSLSIZE_BUCK] = {0.203956 , 0.0039556 , 0 , 0 , 0 , 0 , 1 , 0.5 , 1.03508 , 0.53508 , 0 , 0.76492};
    double wexa[NSLSIZE_BUCK];
    double bufzexa[NSLSIZE_BUCK];

    clock_t timeLU;
    int sizesublcp = 7;
    int sizelcp = NSLSIZE_BUCK;
    int infoLU;
    int nboperlu = 200000;

    timeLU = clock();
    for(int k=0;k< nboperlu;k++)
    {
        for(int j = 0;j < sizesublcp;j++)
        {
            for(int i = 0;i < sizesublcp;i++)
            {
                submatlcp[(j*sizesublcp)+i] = M_straight[(indic[j]*sizelcp)+indic[i]];
            }
        }
        dgetrfmine(&sizesublcp, &sizesublcp, submatlcp, &sizesublcp, ipiv,&infoLU);
        if (infoLU != 0) {cerr << "Pb dgetrf !" << endl;exit(-1);}
    }
    timeLU = clock() - timeLU;
    double timefactLU = (double)(timeLU)/(double)(nboperlu);
    cout << "temps " << nboperlu << " fact LU = " << timeLU << " temps elem = " << timefactLU << endl;
    cerr << "temps " << nboperlu << " fact LU = " << timeLU << " temps elem = " << timefactLU << endl;

    timeLU = clock();
    for(int k=0;k< nboperlu;k++)
    {
        dcopymine( &sizelcp , zexa , &incx , bufzexa , &incx );
        for(int i = 0;i < sizelcp;i++) zexa[i] = 0.;
        for(int i = 0;i < sizesublcp;i++)
        {
            subq[i] = -qexa[indic[i]];
        }
        dgetrsmine(&NOTRANS, &sizesublcp, &incx, submatlcp, &sizesublcp,
            ipiv, subq, &sizesublcp, &infoLU);
        if (infoLU != 0) {cerr << "Pb dgetrs !" << endl;exit(-1);}

        for(int i = 0;i < sizesublcp;i++) zexa[indic[i]] = subq[i];
//        infoLU = filter_result_LCP(sizelcp,M_straight,qexa,zexa,solvingMethod->lcp.tol,solvingMethod->lcp.chat,wexa);
    }
    timeLU = clock() - timeLU;
    double timesolveLU = (double)(timeLU)/(double)(nboperlu);
    cout << "temps " << nboperlu << " solve LU = " << timeLU << " temps elem = " << timesolveLU << endl;
    cerr << "temps " << nboperlu << " solve LU = " << timeLU << " temps elem = " << timesolveLU << endl;

    timeLU = clock();
    for(int k=0;k< nboperlu;k++)
    {
        infoLU = filter_result_LCP(sizelcp,M_straight,qexa,zexa,solvingMethod->lcp.tol,solvingMethod->lcp.chat,wexa);
    }
    timeLU = clock() - timeLU;
    double timefilterresult = (double)(timeLU)/(double)(nboperlu);
    cout << "temps " << nboperlu << " filter result = " << timeLU << " temps elem = " << timefilterresult << endl;
    cerr << "temps " << nboperlu << " filter result = " << timeLU << " temps elem = " << timefilterresult << endl;
    if (infoLU != 0) {cerr << "mauvaise sol filter !" << endl;exit(-1);}

    infoLU = lcp_solver( M_straight , qexa , &sizelcp , solvingMethod , zexa , wexa );
    if (infoLU != 0) {cerr << "mauvaise sol lcp !" << endl;exit(-1);}
    infoLU = lcp_solver( M_straight , qexa , &sizelcp , solvingMethod , zexa , wexa );
    int nbiterlcp;
    int nbitersign;
    int nbiter1eq;
    int nbiterlcpsolverrpgs = 0;
    nbiterlcp  = solvingMethod->lcp.iter;
    nbitersign = solvingMethod->lcp.nbitersign;
    nbiter1eq  = solvingMethod->lcp.nbiter1eq;
    if ( !((nbiterlcp == 1) && (nbitersign == 0) && (nbiter1eq == 0)) ) {cerr << "pas de sol pred !" << endl;exit(-1);}

    timeLU = clock();
    for(int k=0;k< nboperlu;k++)
    {
        for(int i = 0;i < sizelcp;i++) zexa[i] = 0.;
        infoLU = lcp_solver( M_straight , qexa , &sizelcp , solvingMethod , zexa , wexa );
        nbiterlcpsolverrpgs += solvingMethod->lcp.iter;
    }
    timeLU = clock() - timeLU;
//    double timelcp_solver = (double)(timeLU)/(double)(nboperlu);
    double timelcp_solver = (double)(timeLU)/(double)(nbiterlcpsolverrpgs);
    cout << "temps " << nboperlu << " lcp solver = " << timeLU << " temps elem = " << timelcp_solver << endl;
    cerr << "temps " << nboperlu << " lcp solver = " << timeLU << " temps elem = " << timelcp_solver << endl;
    cerr << "nbiterlcpsolverrpgs = " << nbiterlcpsolverrpgs << endl;

    timeLU = clock();
    for(int k=0;k< nboperlu;k++)
    {
        infoLU = lcp_solver( M_straight , qexa , &sizelcp , solvingMethodBackup , zexa , wexa );
    }
    timeLU = clock() - timeLU;
    double timesolverbackup = (double)(timeLU)/(double)(nboperlu);
    cout << "temps " << nboperlu << " lcp solver bck = " << timeLU << " temps elem = " << timesolverbackup << endl;
    cerr << "temps " << nboperlu << " lcp solver bck = " << timeLU << " temps elem = " << timesolverbackup << endl;
    if (infoLU != 0) {cerr << "mauvaise sol lcp bck !" << endl;exit(-1);}

    int nbiterrpgs = 0;
    int     iparamLCP[5];
    double  dparamLCP[5];
    timeLU = clock();
    for(int k=0;k< nboperlu;k++)
    {
        for(int i = 0;i < sizelcp;i++) zexa[i] = 0.;
        iparamLCP[0] = solvingMethod->lcp.itermax;
        iparamLCP[1] = solvingMethod->lcp.chat;
        dparamLCP[0] = solvingMethod->lcp.tol;
        dparamLCP[1] = solvingMethod->lcp.rho;

        lcp_rpgs( &sizelcp , M_straight , qexa , zexa , wexa , &infoLU , iparamLCP , dparamLCP );
        if (infoLU != 0) {cerr << "Pb rpgs !" << endl;exit(-1);}

        nbiterrpgs += iparamLCP[2];
    }
    timeLU = clock() - timeLU;
    double timerpgs = (double)(timeLU)/(double)(nbiterrpgs);
    cout << "temps " << nboperlu << " lcp solver rpgs = " << timeLU << " nb iters = " << nbiterrpgs << " temps elem = " << timerpgs << endl;
    cerr << "temps " << nboperlu << " lcp solver rpgs = " << timeLU << " nb iters = " << nbiterrpgs << " temps elem = " << timerpgs << endl;

    clock_t timeLU0;
    timeLU = clock();
    for(int k=0;k< nboperlu;k++)
    {
        timeLU0 = clock();
    }
    timeLU = clock() - timeLU;
    double timeclock = (double)(timeLU)/(double)(nboperlu);
    cout << "temps " << nboperlu << " clock() = " << timeLU << " temps elem = " << timeclock << endl;
    cerr << "temps " << nboperlu << " clock() = " << timeLU << " temps elem = " << timeclock << endl;
*/
//**************************************************************************************************
//     SparseBlockStructuredMatrix *Msparse;
//     Msparse = LCP_BuckConverter->getMspblPtr();
//     int sizeblcol;
//     double maxabsval;
//
//     ofstream outFileMat2("calcmat.sce");          // checks that it's opened
//     if ( !outFileMat2.is_open() )
//     {
//         cout << "function write error : Fail to open \"calcmat.sce\"" << endl;
//         exit(-1);
//     }
//
//     outFileMat2.precision(10);
//     cout << endl << "   LCP sparse matrix :  " << endl;
//     for (unsigned int i = 0 ; i < Msparse->nbblocks ; i++)
//     {
//         outFileMat2 << "block" << Msparse->RowIndex[i] << "_" << Msparse->ColumnIndex[i] << " = [.." << endl;
//         cout << endl << " block n° " << i << " at (" << Msparse->RowIndex[i] << "," << Msparse->ColumnIndex[i] << ")" << endl << endl;
//         sizeblcol = Msparse->blocksize[Msparse->RowIndex[i]];
//         maxabsval = 0.;
//         for (int j = 0 ; j < sizeblcol ; j++)
//         {
//             for (int k = 0 ; k < Msparse->blocksize[Msparse->ColumnIndex[i]] ; k++) {
//               cout << Msparse->block[i][(k*sizeblcol) + j] << " ";
//               outFileMat2 << Msparse->block[i][(k*sizeblcol) + j] << " ";
//               if (abs(Msparse->block[i][(k*sizeblcol) + j]) > maxabsval) maxabsval = Msparse->block[i][(k*sizeblcol) + j];
//             }
//             if (j != sizeblcol-1) outFileMat2 << " ;" << endl;
//             else outFileMat2 << " ]" << endl;
//             cout << endl;
//         }
//         cout << " max value in block = " << maxabsval << endl;
//     }
//     outFileMat2.close();

//**************************************************************************************************

    SiconosVector *xpt,*lambdapt_buck;
	SimpleVector LambdaP(2*NBHYP),LambdaN(2*NBHYP);

    xpt = LSBuckConverter->getXPtr();
    lambdapt_buck = InterBuckConverter_buck->getLambdaPtr(0);

    unsigned int k = 0;
    unsigned int N = TiDisc->getNSteps(); // Number of time steps
    tinst = k*h_step;

    int nbiterlcp = 0;
    int nbitersign = 0;
    int nbiter1eq = 0;

    int nbiterlcptotal = 0;
    int diffnbiter1eq = 0;
    int nbbackup = 0;
    char str_tinst[20];
    int nbitersupp = 0;
    int nb0iter = 0;
    int nb1iternochg = 0;
    int nbiternochg = 0;
    int longserienochg = 0;
    int nbsolvelu = 0;
    int solvebck;
    int nbiterbackup = 0;
    int sizesublcp;
    int nhisto;
    int histo[10] = {0,0,0,0,0,0,0,0,0,0};
    int nbsolvestd = 0;

    dataPlotini = (double*) malloc((1 + ((N+1)*nbPlot))*sizeof(double));
    if (dataPlotini == NULL)
    {
        cout << " Pb allocation memoire dataPlot !!" << endl;
        exit(-1);
    }
    dataPlot = dataPlotini;

    // For the initial time step:

    // time
    *(++dataPlot) = k*h_step;

    // ramp voltage
    *(++dataPlot) = (LSBuckConverter->getZ())(SIZEZ_PAR);

    // MOS P drain potential
    *(++dataPlot) =  - VthDN;

	// L current
	*(++dataPlot) = x_straight[1];

	// output voltage
	*(++dataPlot) = x_straight[0];

    // gate voltage
    *(++dataPlot) = SlopeComp * (z_straight[0] - z_straight[1]);

	// error voltage
	*(++dataPlot) = x_straight[4];

    // DPMOS current
    *(++dataPlot) = z_straight[NSLSIZE_BUCK-2];

    // DNMOS current
    *(++dataPlot) = w_straight[NSLSIZE_BUCK-1];

    // PMOS Isd current
    rowsize = 2*NBHYP;
    *(++dataPlot) = HalfKP * ddotmine( &rowsize , fPWLmat_straight , &incx , &(z_straight[2]) , &incy );

    // NMOS Ids current
    *(++dataPlot) = HalfKN * ddotmine( &rowsize , fPWLmat_straight , &incx , &(z_straight[2+(2*NBHYP)]) , &incy );

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

    double maxdiffopa = 0;

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
        tinst = k*h_step;
//         cout << "  k = " << k << endl;
        sprintf(str_tinst,"%1.10e ",tinst);

        LSBuckConverter->computeB(tinst);

        for(unsigned int i = 0;i < SIZEX;i++) xtemp_straight[i] = x_straight[i] + (h_step * b_straight[i]);

        a1 = 1.;
        b1 = 0.;
        rowsize = SIZEX;
        colsize = SIZEX;
        dgemvmine( &NOTRANS , &rowsize , &colsize , &a1 , W_straight , &rowsize , xtemp_straight ,
                &incx , &b1 , xfree_straight , &incy );

        b1 = 1.;
        rowsize = NSLSIZE_BUCK;
        colsize = SIZEZ_PAR + SIZEZ_INP;
        dcopymine( &rowsize ,e_buck_straight,&incx,q_straight,&incy );
        dgemvmine( &NOTRANS , &rowsize , &colsize , &a1 , F0_buck_straight , &rowsize , zpar_straight ,
                &incx , &b1 , q_straight , &incy );

        colsize = SIZEX;
        dgemvmine( &NOTRANS , &rowsize , &colsize , &a1 , C_buck_straight , &rowsize , xfree_straight ,
                &incx , &b1 , q_straight , &incy );

        dcopymine( &rowsize ,z_straight,&incx,zprev_straight,&incy );
        dcopymine( &rowsize ,w_straight,&incx,wprev_straight,&incy );

        solvebck = 0;
//        startLCPsolve = clock();
//        startLCPuni = clock();
        info = lcp_solver( M_straight , q_straight , &rowsize , solvingMethod , z_straight , w_straight );
//        timeLCPuni = clock() - startLCPuni;

//         nbiterlcp = solvingMethod->lcp.iter;
//         nbitersign = solvingMethod->lcp.nbitersign;
//         nbiter1eq = solvingMethod->lcp.nbiter1eq;

//         if ( !((nbiterlcp == 1) && (nbitersign == 0) && (nbiter1eq == 0)) )
//         {
//             nbiterlcptotal += nbiterlcp;
//         } else {
//             sizesublcp = solvingMethod->lcp.sizesublcp;
//             nhisto = 1;
//             while(sizesublcp >= 10*nhisto) nhisto++;
//             histo[nhisto-1]++;
//         }

        if (info != 0)
        {
            dcopymine( &rowsize ,zprev_straight,&incx,z_straight,&incy );
            dcopymine( &rowsize ,wprev_straight,&incx,w_straight,&incy );

//            startLCPuni = clock();
            info = lcp_solver( M_straight , q_straight , &rowsize , solvingMethodBackup , z_straight , w_straight );
//            timeLCPuni = clock() - startLCPuni;
//             nbiterbackup = solvingMethodBackup->lcp.iter;
//             nbbackup++;
            solvebck = 1;
//            LCP_CPUtime_bck += timeLCPuni;
        } else {
//             nbiterbackup = 0;
//          LCP_CPUtime_std += timeLCPuni;
        }

//        LCP_CPUtime += (clock() - startLCPsolve);
        if (info != 0)
        {
            cout << "!!! simulation breakdown at time " << tinst << " !!!" << endl;
            cerr << "!!! simulation breakdown at time " << tinst << " !!!" << endl;
        }

//         if (nbiterlcp > nbitersign) nbitersupp += (nbiterlcp - (nbitersign+1));
//         if ((solvebck == 0) && (nbiterlcp == 0)) nb0iter++;
//         if (nbitersign == 0) longserienochg++;
//         else {
//             if (longserienochg != 0) cout << "fin de serie no chg " << longserienochg << " à " << str_tinst << endl;
//             longserienochg = 0;
//         }

//         if ((solvebck == 0) && (nbiterlcp >= 1))
//         {
//             if ((nbiterlcp == 1) && (nbitersign == 0) && (nbiter1eq == 0)) nbsolvelu++;
//             if ((nbiterlcp == 1) && (nbitersign == 0) && (nbiter1eq == 1)) nb1iternochg++;
//             if ((nbitersign == 0) && (nbiter1eq != 0)) nbiternochg++;
//             if ((nbitersign != 0) || (nbiter1eq != 0)) nbsolvestd++;
//         }
//        if ( ((nbitersign+1) == nbiter1eq) || (nbiter1eq == 0) )
//         if (nbiter1eq == 0)
//         {
//             diffnbiter1eq = 0;
//         } else {
//             if (nbitersign > nbiter1eq)
//             {
//                 diffnbiter1eq = nbitersign - nbiter1eq;
// //                 cout << "chgt mode1 à " << str_tinst << " RPGS iters = " << solvingMethod->lcp.iter;
// //                 cout << " nbiter1eq = " << nbiter1eq << " nbitersign = " << nbitersign << endl;
//             } else {
// //                 cout << "Erreur sur lcp iter sign à " << str_tinst << endl;
// //                 info = 1;
//             }
//         }
//         if ((solvingMethod->lcp.iter != 0) && (nbiter1eq == 0))
//         {
//             cout << "chgt signe continu à " << str_tinst << " RPGS iters = " << solvingMethod->lcp.iter << endl;
//         }

        rowsize = SIZEX;
        dcopymine( &rowsize ,xfree_straight,&incx,x_straight,&incy );

        b1 = 0.;
        colsize = NSLSIZE_BUCK;
        dgemvmine( &NOTRANS , &rowsize , &colsize , &a1 , B_buck_straight , &rowsize , z_straight ,
                &incx , &b1 , r_straight , &incy );

        b1 = 1.;
        a1 = h_step;
        colsize = SIZEX;
        dgemvmine( &NOTRANS , &rowsize , &colsize , &a1 , W_straight , &rowsize , r_straight ,
                &incx , &b1 , x_straight , &incy );

        if (abs(x_straight[4]) > maxdiffopa) maxdiffopa = abs(x_straight[4]);
        // time
        *(++dataPlot) = k*h_step;

        // ramp voltage
        *(++dataPlot) = (LSBuckConverter->getZ())(SIZEZ_PAR);

        // MOS P drain potential
        *(++dataPlot) = z_straight[NSLSIZE_BUCK-1] - VthDN;

        // L current
        *(++dataPlot) = x_straight[1];

        // output voltage
        *(++dataPlot) = x_straight[0];

        // gate voltage
        *(++dataPlot) = SlopeComp * (z_straight[0] - z_straight[1]);

        // error voltage
        *(++dataPlot) = x_straight[4];

        // DPMOS current
        *(++dataPlot) = z_straight[NSLSIZE_BUCK-2];

        // DNMOS current
        *(++dataPlot) = w_straight[NSLSIZE_BUCK-1];

        // PMOS Isd current
        rowsize = 2*NBHYP;
        *(++dataPlot) = HalfKP * ddotmine( &rowsize , fPWLmat_straight , &incx , &(z_straight[2]) , &incy );

        // NMOS Ids current
        *(++dataPlot) = HalfKN * ddotmine( &rowsize , fPWLmat_straight , &incx , &(z_straight[2+(2*NBHYP)]) , &incy );

        // nb iterations lcp solver
//        *(++dataPlot) = nbiterlcp + nbiterbackup;

        // solver backup
//        *(++dataPlot) = solvebck;

        // nb iterations lcp solver 1eq
//        *(++dataPlot) = nbiter1eq;

//         cout << endl << "time = " << dataPlot(k,0) << " Vramp = " << dataPlot(k,1) << " Vd MOSP = " << dataPlot(k,2);
//         cout << " IL = " << dataPlot(k,3) << " Vout = " << dataPlot(k,4) << endl << endl;

        if ((k % (N/100)) == 0)
        {
            cerr << "-------- " << (100.0*k)/N << " % achieved..." << endl;
//             cerr << "nb CPU cycles LCP = " << LCP_CPUtime << endl;
        }
    }

  } // end of "try" section
  // --- Exceptions handling ---
  catch(SiconosException e)
  {
    cout << "Time step n° " << k << " : " << tinst << endl;
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
//     cout << "nb CPU cycles LCP = " << LCP_CPUtime;
//     cout << " , nb CPU cycles loop = " << endloop-start << endl;
//     // Number of time iterations
//     cout << "Number of time steps : " << k << " with " << nbbackup << " backup solving and " << nb0iter << " direct solution." << endl;
//     cout << "Nb sol 1 iter no chg = " << nb1iternochg << "  Nb sol no chg = " << nbiternochg << endl;
//     cout << "Nb solve lu = " << nbsolvelu << " Nb solve std = " << nbsolvestd << " check = " << nbsolvestd + nbsolvelu + nbbackup + nb0iter << endl;
//     cout << "Total number of lcp iterations = " << nbiterlcptotal << " nb iter supp = " << nbitersupp << endl;
// //    cout << "lcp backup iterations = " << nbiterbackup << " average = " << (double)(nbiterbackup)/(double)(nbbackup) << endl;
//     cout << "histogramme sizesublcp = ";
//     int sumhisto = 0;
//     for (int i=0;i<10;i++) {cout << histo[i] << "  ";sumhisto += histo[i];}
//     cout << endl;
//     cout << "sumhisto = " << sumhisto << endl;


//     double timelcpestim = 0;
//     timelcpestim += (k-nbsolvelu)*timefactLU;
//     timelcpestim += k*timelcp_solver;
//    timelcpestim += k*2*timeclock;
//     timelcpestim += nbiterlcptotal*timerpgs;
//     timelcpestim += nbbackup*timesolverbackup;
//     cout << "LCP CPU time estimation = " << timelcpestim/(double)CLOCKS_PER_SEC << " s with CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << endl;
//     cerr << "LCP CPU time estimation = " << timelcpestim/(double)CLOCKS_PER_SEC << " s with CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << endl;

//    cout << "Average number of lcp iterations = " << nbiterlcptotal/(double)k << " iterations per time step." << endl;

//    LCP_BuckConverter->printStat();
//     cout << "Average number of tried LCP matrices: ";
//     cout << (double)(LCP_BuckConverter->getStatnbessai())/(double)(LCP_BuckConverter->getStatnbsolve()) << endl;
//      cout << "Percent of lcp solved with initial algo : ";
//      cout << 100.0*(double)(LCP_BuckConverter->getStatsolverstd())/(double)(LCP_BuckConverter->getStatnbsolve()) << endl;
//    cout << "Max absolute opa output voltage : " << maxdiffopa << endl;

    dataPlot = dataPlotini;
    for (unsigned int i=0;i < k;i++)
    {
      if (!(i % nblignebuf)) sprintf(bufferligne,"");
      sprintf(buffer,"%1.10e ",*(++dataPlot));
      strcat(bufferligne,buffer);
      for (unsigned int j=1; j < nbPlot; j++)
      {
	      sprintf(buffer,"%1.4e ",*(++dataPlot));
          strcat(bufferligne,buffer);
      }
      sprintf(buffer,"\n");
      strcat(bufferligne,buffer);
      if (!((i+1) % nblignebuf)) outFile << bufferligne;
    }
    if ((k) % nblignebuf) outFile << bufferligne;
    outFile.close();

    free(dataPlotini);

//  Saving the current state for simulation restart
//     ofstream outFileState("BuckConverterstate.dat");          // checks that it's opened
//     if ( !outFileState.is_open() )
//     {
//         cout << "function write error : Fail to open \"BuckConverterstate.dat\"" << endl;
//     } else {
//         outFileState.precision(12);
//         for (unsigned int i=0;i < SIZEX;i++) outFileState << xpt->getValue(i) << endl;
//         outFileState.close();
//     }

    delete LCP_BuckConverter;
    delete OSI_LSBuckConverter;
    delete TiDisc;
    delete StratBuckConverter;
    delete NSDSBuckConverter;
    delete InterBuckConverter_buck;
    delete LTIRBuckConverter_buck;
    delete Int_B_buck;
    delete Int_C_buck;
    delete Int_F0_buck;
    delete Int_e_buck;
    delete Int_D_buck;
    delete nslaw_buck;
    delete paramVin;
    delete Coltemp;
//	delete LS_T;
    delete LSBuckConverter;
    cout << "End of program" << endl;
}
