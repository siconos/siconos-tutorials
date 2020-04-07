#include <chrono>
#include "SiconosKernel.hpp"

using namespace std;

// main program
int main(int argc, char* argv[])
{
  // Exception handling
  try
  {
    // == User-defined parameters ==
    unsigned int ndof = 2;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 2.0;        // Total simulation times
    double h = 1.0e-3;      // Time step
    double Vinit = 3.0;


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
    (*A)(0, 1) = 0.0;
    (*A)(1, 0) = 0.0;
    (*A)(1, 1) = 0.0;
    SP::SiconosVector x0(new SiconosVector(ndof));
    (*x0)(0) = Vinit;
    (*x0)(1) = Vinit;
    SP::FirstOrderLinearDS process(new FirstOrderLinearDS(x0, A));
    //    process->setComputebFunction("ObserverLCSPlugin","uProcess");


    double c =25.0;

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 2; // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    SP::SimpleMatrix B(new SimpleMatrix(ndof, ninter));
    (*B)(0, 0) = 1.0;
    (*B)(1, 0) = 1.0+c;
    (*B)(0, 1) = -(1.0+c);
    (*B)(1, 1) = 1.0;
    *B = 2.0 * (*B);
    SP::SimpleMatrix C(new SimpleMatrix(ninter, ndof));
    (*C)(0, 0) = 1.0;
    (*C)(1, 0) = 0.0;
    (*C)(0, 1) = 0.0;
    (*C)(1, 1) = 1.0;
    SP::FirstOrderLinearR myProcessRelation(new FirstOrderLinearR(C, B));
    SP::SimpleMatrix D(new SimpleMatrix(ninter, ninter));
    (*D)(0, 0) = 0.0;
    (*D)(0, 1) = 0.0;
    (*D)(1, 0) = 0.0;
    (*D)(1, 1) = 0.0;

    myProcessRelation->setDPtr(D);
    //myProcessRelation->setComputeEFunction("ObserverLCSPlugin","computeE");

    // Second relation, related to the observer
    // haty = C hatX + D hatLambda + E
    // hatR = B hatLambda

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    SP::NonSmoothLaw myNslaw(new RelayNSL(nslawSize));

    myNslaw->display();

    // The Interaction which involves the first DS (the process)
    string nameInter = "processInteraction"; // Name

    SP::Interaction myProcessInteraction(new Interaction(myNslaw, myProcessRelation));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem filippov(new NonSmoothDynamicalSystem(t0, T));
    filippov->insertDynamicalSystem(process);
    filippov->link(myProcessInteraction,process);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));
    // == Creation of the Simulation ==
    SP::TimeStepping s(new TimeStepping(filippov, td));
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
    unsigned int outputSize = 7; // number of required data
    unsigned int N = ceil((T - t0) / h); // Number of time steps

    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector xProc = process->x();
    SP::SiconosVector lambdaProc = myProcessInteraction->lambda(0);
    SP::SiconosVector yProc = myProcessInteraction->y(0);

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = filippov->t0(); // Initial time of the model
    dataPlot(0, 1) = (*xProc)(0);
    dataPlot(0, 2) = (*xProc)(1);
    dataPlot(0, 3) = -1.0;
    dataPlot(0, 4) = -1.0;
    dataPlot(0, 5) = (*xProc)(0);
    dataPlot(0, 6) = (*xProc)(1);



    // ==== Simulation loop =====
    cout << "====> Start computation ... " << endl << endl;

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    unsigned int k = 0; // Current step

    // Simulation loop
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
        while (k < N - 1)
    {
      k++;
      //  osnspb->setNumericsVerboseMode(1);

      //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
      s->computeOneStep();
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      dataPlot(k, 3) = (*lambdaProc)(0);
      dataPlot(k, 4) = (*lambdaProc)(1);
      dataPlot(k, 5) = (*yProc)(0);
      dataPlot(k, 6) = (*yProc)(1);
      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
cout << "Computation Time " << endl;;
    time.report();    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("Filippov.dat", "ascii", dataPlot, "noDim");
    double error=0.0, eps=1e-05;
    if ((error=ioMatrix::compareRefFile(dataPlot, "Filippov.ref", eps)) >= 0.0
        && error > eps)
      return 1;
  }

  catch (SiconosException e)
  {
    cerr << e.report() << endl;
    return 1;
  }
  catch (...)
  {
    cerr << "Exception caught in Fillipov.cpp" << endl;
    return 1;
  }
}
