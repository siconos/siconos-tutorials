#include "SiconosKernel.hpp"
#include <math.h>
#include "SolverOptions.h"
using namespace std;

// main program

#include <chrono>


int main(int argc, char* argv[])
{
  // Exception handling
  try
  {
    // == User-defined parameters ==
    unsigned int ndof = 2;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 100;        // Total simulation times
    double h = 1.0e-2;      // Time step
    double xinit = 0.0;
    double vinit = 0.0;
    char  filename[50] = "simu.";
    if(argc==1)
    {
      xinit = 12.0;
      vinit =6.0;
      strncpy(&filename[5],"1.0.1.0.log",7);
    }
    else if(argc==3)
    {
      //printf("argv[0] %s\n", argv[0]);
      printf("xinit is set to %f\n", atof(argv[1]));
      printf("vinit is set to %f\n", atof(argv[2]));

      xinit = atof(argv[1]);
      vinit = atof(argv[2]);
      int sizeofargv1 = strlen(argv[1]);
      // printf("sizeofargv1 %i\n",sizeofargv1);
      strncpy(&filename[5],argv[1],sizeofargv1);
      int sizeofargv2 = strlen(argv[2]);
      //printf("sizeofargv2 %i\n",sizeofargv2);
      strncpy(&filename[5+sizeofargv1],".",1);

      strncpy(&filename[5+sizeofargv1+1],argv[2],sizeofargv2);
      strncpy(&filename[5+sizeofargv1+sizeofargv2+1],".log",4);


      // printf("Output is written in filename %s\n",  filename);
    }
    else
    {
      cout << "wrong  number of arguments = " << argc << endl;
    }

    double m=1, stiffness=1;
    double alpha =1.0;

    // ================= Creation of the model =======================
    // Steps:
    // - create some Dynamical Systems
    // - create some Interactions between those Dynamical Systems
    //   Interaction = some relations (constraints) and a NonSmoothLaw
    // - create a NonSmoothDynamicalSystem with the DynamicalSystems and the Interactions
    // - add this NonSmoothDynamicalSystem into a Model
    // - add a Simulation to the model
    //  Simulation = TimeDiscretisation + OneStepIntegrator and OneStepNSProblem

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // First System:
    // dx/dt = Ax + u(t) + r
    // x(0) = x0
    // Note: r = Blambda, B defines in relation below.

    SP::SiconosMatrix A(new SimpleMatrix(ndof, ndof));
    (*A)(0, 0) = 0.0;
    (*A)(0, 1) = 1.0;
    (*A)(1, 0) = -stiffness/m;
    (*A)(1, 1) = 0.0;
    SP::SiconosVector x0(new SiconosVector(ndof));
    (*x0)(0) = xinit;
    (*x0)(1) = vinit;

    SP::FirstOrderLinearDS process(new FirstOrderLinearDS(x0, A));
    //    process->setComputebFunction("ObserverLCSPlugin","uProcess");

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 1; // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    SP::SimpleMatrix B(new SimpleMatrix(ndof, ninter));
    (*B)(0, 0) = 0.0;
    (*B)(1, 0) = alpha;

    SP::SimpleMatrix C(new SimpleMatrix(ninter, ndof));
    (*C)(0, 0) = 0.0;
    (*C)(0, 1) = 1.0;

    SP::FirstOrderLinearR myProcessRelation(new FirstOrderLinearR(C, B));

    // NonSmoothLaw
    unsigned int nslawSize = 1;
    SP::NonSmoothLaw myNslaw(new RelayNSL(nslawSize, -1., 1.));

    myNslaw->display();

    // The Interaction which involves the first DS (the process)
    SP::Interaction myProcessInteraction(new Interaction(myNslaw, myProcessRelation));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem relayOscillator(new NonSmoothDynamicalSystem(t0, T));
    relayOscillator->insertDynamicalSystem(process);
    relayOscillator->link(myProcessInteraction, process);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));
    // == Creation of the Simulation ==
    SP::TimeStepping s(new TimeStepping(relayOscillator, td));
    // -- OneStepIntegrators --
    double theta = 0.5;
    SP::EulerMoreauOSI myIntegrator(new EulerMoreauOSI(theta));
    s->insertIntegrator(myIntegrator);

    // -- OneStepNsProblem --

    SP::Relay osnspb(new Relay());


    osnspb->setSolverId(SICONOS_RELAY_LEMKE);
    osnspb->numericsSolverOptions()->dparam[0] = 1e-08;
    s->insertNonSmoothProblem(osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    unsigned int outputSize =10; // number of required data
    unsigned int N = ceil((T - t0) / h); // Number of time steps

    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector xProc = process->x();
    SP::SiconosVector lambdaProc = myProcessInteraction->lambda(0);
    SP::SiconosVector yProc = myProcessInteraction->y(0);
    SP::SiconosVector  vectorfield = process->rhs();
    unsigned int k = 0; // Current step


    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = relayOscillator->t0(); // Initial time of the model
    dataPlot(k, 1) = (*xProc)(0);
    dataPlot(k, 2) = (*xProc)(1);
    dataPlot(k, 3) = (*lambdaProc)(0);
    dataPlot(k, 4) = (*yProc)(0);
    dataPlot(k, 7) = vectorfield->getValue(0);
    dataPlot(k, 8) = vectorfield->getValue(1);


    // ==== Simulation loop =====
    cout << "====> Start computation ... " << endl << endl;

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    // Simulation loop
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    while(k < N - 1)
    {
      k++;

      //  osnspb->setNumericsVerboseMode(1);

      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      dataPlot(k, 3) = (*lambdaProc)(0);
      dataPlot(k, 4) = (*yProc)(0);
      process->computeRhs(s->nextTime());
      if(k==1)  // tricks just for display to avoid the computation of the initial Rhs
      {
        dataPlot(k-1, 7) = vectorfield->getValue(0);
        dataPlot(k-1, 8) = vectorfield->getValue(1);
      }

      dataPlot(k, 7) = vectorfield->getValue(0);
      dataPlot(k, 8) = vectorfield->getValue(1);
      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << endl;;
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << "Computation time : " << elapsed << " ms" << endl;
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("FrictionOscillator.dat", "ascii", dataPlot, "noDim");
    ioMatrix::write(filename, "ascii", dataPlot, "noDim");
    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("FrictionOscillator.ref", "ascii", dataPlotRef);
    std::cout << (dataPlot-dataPlotRef).normInf() <<std::endl;
    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      return 1;
    }


  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }
}
