/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
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
#include <SiconosKernel.hpp>
#include <chrono>

using namespace std;
using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

// main program
int main(int argc, char* argv[]) {
  // Exception handling
  try {
    // == User-defined parameters ==
    unsigned int ndof = 2;  // number of degrees of freedom of your system
    double t0 = 0.0;
    double T = 1;                 // Total simulation time
    double h = 1.0e-4;            // Time step
    double hcontroller = 1.0e-2;  // Time step
    double Vinit = 1.0;

    if (h > hcontroller) {
      THROW_EXCEPTION("hcontroller must be larger than h");
    }

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
    // Note: r = Blambda, B defined in relation below.

    auto x0 = std::make_shared<Vector>(ndof);
    (*x0)(0) = Vinit;
    (*x0)(1) = -Vinit;

    auto processDS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        *x0, siconos::algebra::SiconosDenseMatrix::Zero(x0->size(), x0->size()),
        siconos::algebra::SiconosVector::Zero(x0->size()), siconos::algebra::copy_t);

    Vector sampledControl{2};
    sampledControl.setZero();  // Will be updated in the time loop.

    processDS->setComputebVectorFunction(
        [&sampledControl](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
          auto input = sin(50 * time);
          result(0) = input + sampledControl(0);
          result(1) = -input + sampledControl(1);
        });

    auto controllerDS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(
        *x0, siconos::algebra::SiconosDenseMatrix::Zero(x0->size(), x0->size()),
        siconos::algebra::SiconosVector::Zero(x0->size()), siconos::algebra::copy_t);

    // --------------------
    // --- Interactions ---
    // --------------------
    unsigned int ninter = 2;  // dimension of your Interaction = size of y and lambda vectors

    // First relation, related to the process
    // y = Cx + Dlambda
    // r = Blambda
    auto B = std::make_shared<Matrix>(ndof, ninter);
    B->setZero();
    (*B)(0, 0) = 2.0;
    (*B)(1, 1) = 2.0;
    auto C = std::make_shared<Matrix>(ninter, ndof);
    C->setZero();
    (*C)(0, 0) = 1.0;
    (*C)(1, 1) = 1.0;

    auto myControllerRelation = std::make_shared<siconos::modeling::FirstOrderLinearR>();
    myControllerRelation->setConstantB(*B);
    myControllerRelation->setConstantC(*C);

    // NonSmoothLaw
    unsigned int nslawSize = 2;
    auto myNslaw = std::make_shared<siconos::modeling::RelayNSL>(nslawSize);

    myNslaw->display();

    // The Interaction which involves the first DS (the process)
    auto myControllerInteraction =
        std::make_shared<siconos::modeling::Interaction>(myNslaw, myControllerRelation);

    // -------------
    // --- Model process ---
    // -------------
    auto process = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    process->insertDynamicalSystem(processDS);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    auto processTD = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // == Creation of the Simulation ==
    auto processSimulation =
        std::make_shared<siconos::simulation::TimeStepping>(process, processTD, 0);
    processSimulation->setName("Simulation of the process");
    // -- OneStepIntegrators --
    double theta = 0.5;
    auto processIntegrator = std::make_shared<siconos::integrators::EulerMoreauOSI>(theta);
    processSimulation->insertIntegrator(processIntegrator);

    // -------------
    // --- Model controller ---
    // -------------
    auto controller = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    controller->insertDynamicalSystem(controllerDS);
    controller->link(myControllerInteraction, controllerDS);

    // ------------------
    // --- Simulation ---
    // ------------------
    // TimeDiscretisation
    auto controllerTD =
        std::make_shared<siconos::simulation::TimeDiscretisation>(t0, hcontroller);
    // == Creation of the Simulation ==
    auto controllerSimulation =
        std::make_shared<siconos::simulation::TimeStepping>(controller, controllerTD);
    controllerSimulation->setName("Simulation of the controller");
    // -- OneStepIntegrators --
    double controllertheta = 0.5;
    auto controllerIntegrator =
        std::make_shared<siconos::integrators::EulerMoreauOSI>(controllertheta);
    controllerSimulation->insertIntegrator(controllerIntegrator);

    // -- OneStepNsProblem --
    auto controllerLCP = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    auto controllerOSNSPB =
        std::make_shared<siconos::nonsmooth_formulations::Relay>(SICONOS_RELAY_PGS);
    controllerSimulation->insertNonSmoothProblem(controllerOSNSPB);

    // coupling the simulation

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Get the values to be plotted ---
    unsigned int outputSize = 10;         // number of required data
    unsigned int N = ceil((T - t0) / h);  // Number of time steps
    Matrix dataPlot(N, outputSize);

    auto xProc = processDS->x();

    // -> saved in a matrix dataPlot
    dataPlot(0, 0) = process->t0();  // Initial time of the model
    dataPlot(0, 1) = (*xProc)(0);

    unsigned int Ncontroller = ceil((T - t0) / hcontroller) + 1;  // Number of time steps
    Matrix dataPlotController(Ncontroller, outputSize);

    auto xController = controllerDS->x();
    auto lambda = myControllerInteraction->lambda(0);
    auto y = myControllerInteraction->y(0);

    // -> saved in a matrix dataPlot
    dataPlotController(0, 0) = controller->t0();  // Initial time of the model
    dataPlotController(0, 1) = (*xController)(0);
    dataPlotController(0, 2) = (*xController)(1);
    dataPlotController(0, 5) = (*lambda)(0);
    dataPlotController(0, 6) = (*lambda)(1);
    dataPlotController(0, 7) = (*y)(0);
    dataPlotController(0, 8) = (*y)(1);

    // ==== Simulation loop =====
    std::cout << "====> Start computation ... \n\n";

    // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    int k = 0;  // Current step
    int kcontroller = 0;
    // Simulation loop
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    while (controllerSimulation->hasNextEvent()) {
      kcontroller++;
      //      cout << "step controller--> " << kcontroller << " at time t =" <<
      //      controllerSimulation->nextTime() << endl;

      // Computation of the controller over the sampling time
      controllerSimulation->computeOneStep();

      //  input of the controller in the process thanks to z and sampledControl
      sampledControl = *B * *lambda;

      while (processSimulation->hasNextEvent() &&
             processSimulation->nextTime() < controllerSimulation->nextTime()) {
        k++;
        //        cout << "         step --> " << k  << " at time t =" <<
        //        processSimulation->nextTime() << endl;

        processSimulation->computeOneStep();
        dataPlot(k, 0) = processSimulation->nextTime();
        dataPlot(k, 1) = (*xProc)(0);
        dataPlot(k, 2) = (*xProc)(1);
        processSimulation->nextStep();
      }
      dataPlotController(kcontroller, 0) =
          controllerSimulation->nextTime();  // Initial time of the model
      dataPlotController(kcontroller, 1) = (*xController)(0);
      dataPlotController(kcontroller, 2) = (*xController)(1);
      dataPlotController(kcontroller, 5) = (*lambda)(0);
      dataPlotController(kcontroller, 6) = (*lambda)(1);
      dataPlotController(kcontroller, 7) = (*y)(0);
      dataPlotController(kcontroller, 8) = (*y)(1);

      // feedback output of the measures for the controller
      *(xController) = *(xProc);

      controllerSimulation->nextStep();
    }
    std::cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    std::cout << "Computation Time \n";
    ;
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    std::cout << "====> Output file writing ...\n";

    siconos::algebra::io::write("RBS-Controller.dat", dataPlotController,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    siconos::algebra::io::write("RBS.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlotController, "RBS-Controller.ref",
                                                      eps)) > eps)
      return 1;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "RBS.ref", eps)) > eps)
      return 1;
    return 0;
  }

  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
