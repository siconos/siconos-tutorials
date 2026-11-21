/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

// This is the program to simulate the dynamic of a rocking block by using the Siconos platform
//==================================================================================================================
#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <algorithm>
#include <chrono>
#include <numbers>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

constexpr double GGearth = 9.8100;

//---------------------------------------------------
double LengthBlock = 1.0;                       // Length of the rocking block
double HeightBlock = 0.5;                       // Height of the rocking block
unsigned int Nfreedom = 3;                      // Number of degrees of freedom
unsigned int Ncontact = 2;                      // Number of contacts
double MassBlock = 1.0;                         // Mass of the rocking block
double PosXiniPointA = 0.0;                     // Initial coordinate X of the point A
double PosYiniPointA = 0.0;                     // Initial coordinate Y of the point A
double AngleThetaIni = std::numbers::pi / 3.0;  // Initial angle theta of the block
double VelXiniPointA = 0.0;                     // Initial relative velocity Vx of the point A
double VelYiniPointA = 0.0;                     // Initial relative velocity Vy of the point A
double RotVelBlockIni = 0.0;                    // Initial angular velocity of the block
double e = 0.9;                                 // Restitution coefficient
double TimeInitial = 0.0;                       // Initial time of the simulation
double TimeFinal = 0.58;                        // Final time of the simulation
double _rho = 0.99;       // used to computer parameters for NewMark Scheme
double StepSize = 0.001;  // Time step size
unsigned int maxIter = 20000;
bool IsTreatFirstSteps = false;
bool IsHandleVelConstraint = false;
//==========================================================================================================
//                                             Main function
//==========================================================================================================
int main(int argc, char* argv[]) {
  try {
    //===========================================================================================================
    //                  I: Declare the dynamical systems
    //===========================================================================================================
    // 1. Set the mass matrix
    auto mass = std::make_shared<Matrix>(Nfreedom, Nfreedom);
    double InertiaBlock;
    InertiaBlock = (MassBlock / 12.0) * ((HeightBlock * HeightBlock) +
                                         (LengthBlock * LengthBlock));  // moment of inertia
    (*mass)(0, 0) = MassBlock;
    (*mass)(1, 1) = MassBlock;
    (*mass)(2, 2) = InertiaBlock;
    // 2. Set the initial position of the block in function of the initial position of the
    // contact point A (left-hand contact)
    auto PosIniBlock = std::make_shared<Vector>(Nfreedom);
    (*PosIniBlock)(0) = PosXiniPointA + 0.5 * LengthBlock * cos(AngleThetaIni) -
                        0.5 * HeightBlock * sin(AngleThetaIni);
    (*PosIniBlock)(1) = PosYiniPointA + 0.5 * LengthBlock * sin(AngleThetaIni) +
                        0.5 * HeightBlock * cos(AngleThetaIni);
    (*PosIniBlock)(2) = AngleThetaIni;
    std::cout.precision(15);
    std::cout << "x0: " << (*PosIniBlock)(0) << "\n";
    std::cout << "y0: " << (*PosIniBlock)(1) << "\n";
    std::cout << "theta0: " << (*PosIniBlock)(2) << "\n";
    std::cout.precision(15);
    // (*PosIniBlock)(0) = 0.5;
    // (*PosIniBlock)(1) = 0.5;
    // (*PosIniBlock)(2) = 0.0;

    // 3. Set the initial velocity of the block in function of the initial relative velocity of
    // the contact point A
    auto VelIniBlock = std::make_shared<Vector>(Nfreedom);
    (*VelIniBlock)(0) = VelXiniPointA - (0.5 * LengthBlock * sin(AngleThetaIni) +
                                         0.5 * HeightBlock * cos(AngleThetaIni)) *
                                            RotVelBlockIni;
    (*VelIniBlock)(1) = VelYiniPointA + (0.5 * LengthBlock * cos(AngleThetaIni) -
                                         0.5 * HeightBlock * sin(AngleThetaIni)) *
                                            RotVelBlockIni;
    (*VelIniBlock)(2) = RotVelBlockIni;

    // (*VelIniBlock)(0) = 0.0;
    // (*VelIniBlock)(1) = 0.0;
    // (*VelIniBlock)(2) = 0.0;

    // 4. Instantiate the object of "LagrangianTIDS"
    auto RockingBlock = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        *PosIniBlock, *VelIniBlock, *mass, siconos::algebra::alias_t);
    // 5. Set the external force
    Vector ExternalForces{Nfreedom};
    ExternalForces.setZero();
    ExternalForces(1) = -MassBlock * GGearth;
    RockingBlock->setConstantFext(ExternalForces, siconos::algebra::alias_t);  //
    std::cout << "Initial position of the rocking block:\n";
    siconos::algebra::print(*PosIniBlock);
    std::cout << "Initial velocity of the rocking block:\n";
    siconos::algebra::print(*VelIniBlock);
    std::cout << "Mass matrix of the rocking block:\n";
    siconos::algebra::print(*mass);
    std::cout << "External force applied on the rocking block:\n";
    siconos::algebra::print(ExternalForces);
    //==================================================================================================================
    //              II: Declare the relation et interaction between dynamical systems
    //==================================================================================================================
    //
    /*
    auto H= std::make_shared<Matrix>(1,Nfreedom);
    (*H)(0,1) = 1.0;
    auto E= std::make_shared<Vector>(1);
    (*E)(0) = -0.5*HeightBlock;
    */
    // Impact law
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    // Interaction at contact point 1
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation1->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = q(1) - 0.5 * LengthBlock * sin(q(2)) - 0.5 * HeightBlock * cos(q(2));
    });

    relation1->setComputeContactMatrixFunction(
        [](const siconos::algebra::BlockVector& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = 0.0;
          result(0, 1) = 1.0;
          result(0, 2) = -0.5 * LengthBlock * cos(q(2)) + 0.5 * HeightBlock * sin(q(2));
        });

    relation1->setComputeContactMatrixDotFunction(
        [](const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(0, 2) =
              (0.5 * LengthBlock * sin(q(2)) + 0.5 * HeightBlock * cos(q(2))) * qdot(2);
        });

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation1);
    // Interaction at contact point 2
    auto relation2 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation2->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = q(1) + 0.5 * LengthBlock * sin(q(2)) - 0.5 * HeightBlock * cos(q(2));
    });

    relation2->setComputeContactMatrixFunction(
        [](const siconos::algebra::BlockVector& q,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result(0, 0) = 0.0;
          result(0, 1) = 1.0;
          result(0, 2) = 0.5 * LengthBlock * cos(q(2)) + 0.5 * HeightBlock * sin(q(2));
        });
    relation2->setComputeContactMatrixDotFunction(
        [](const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot,
           Eigen::Ref<siconos::algebra::MapType> result) {
          result.setZero();
          result(0, 2) =
              (-0.5 * LengthBlock * sin(q(2)) + 0.5 * HeightBlock * cos(q(2))) * qdot(2);
        });

    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation2);
    // Interactions for the whole dynamical system
    //================================================================================================================
    //            III. Create the "model" object
    //================================================================================================================
    auto RoBlockModel =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(TimeInitial, TimeFinal);
    RoBlockModel->insertDynamicalSystem(RockingBlock);
    RoBlockModel->link(inter1, RockingBlock);
    RoBlockModel->link(inter2, RockingBlock);

    //================================================================================================================
    //            IV. Create the simulation
    //================================================================================================================
    // 1. Time discretization
    auto TimeDiscret =
        std::make_shared<siconos::simulation::TimeDiscretisation>(TimeInitial, StepSize);
    // 2. Integration solver for one step
    auto NewMarkAlpha =
        std::make_shared<siconos::integrators::NewMarkAlphaOSI>(_rho, IsHandleVelConstraint);
    // 3. Nonsmooth problem
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto position = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto acceleration = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    position->numericsSolverOptions()->dparam[0] = 1e-12;
    impact->numericsSolverOptions()->dparam[0] = 1e-12;
    acceleration->numericsSolverOptions()->dparam[0] = 1e-12;

    // 4. Simulation with (1), (2), (3)
    auto EDscheme =
        std::make_shared<siconos::simulation::EventDriven>(RoBlockModel, TimeDiscret);
    EDscheme->insertIntegrator(NewMarkAlpha);
    EDscheme->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_ED_IMPACT);
    EDscheme->insertNonSmoothProblem(acceleration,
                                     siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC);
    EDscheme->insertNonSmoothProblem(position,
                                     siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS);
    // bool check1 = EDscheme->hasOneStepNSProblem(impact);
    // bool check2 = EDscheme->hasOneStepNSProblem(acceleration);
    // std::cout << "Impact law included in the simulation: " << check1 << "\n";
    // std::cout << "LCP at acceleration level included in the simulation: " << check2 << "\n";
    //==================================================================================================================
    //                    V. Process the simulation
    //==================================================================================================================
    // -------------------------------- Simulation initialization
    // ------------------------------------------------------
    EDscheme->setPrintStat(true);
    auto eventsManager =
        EDscheme->eventsManager();  // ponters point to the "eventsManager" object
    auto PosBlock =
        RockingBlock->q_read();  // read-only view to the position vector of the rocking block
    auto VelBlock =
        RockingBlock->velocity_read();  // read-only view  to the velocity of the rocking block
    RockingBlock->initRhs(RoBlockModel->t0());
    auto AcceBlock =
        RockingBlock
            ->acceleration_read();  // read-only view to the velocity of the rocking block

    auto indexSet0 = RoBlockModel->topology()->indexSet(0);
    std::cout << "Size of IndexSet0: " << indexSet0->size() << "\n";

    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    //-------------------- Save the output during simulation
    //---------------------------------------------------------
    unsigned int NpointSave = 583;  //
    unsigned int SizeOutput = 13;   //
    Matrix DataPlot(NpointSave, SizeOutput);
    //------------- At the initial time
    //-----------------------------------------------------------------------------
    DataPlot(0, 0) = RoBlockModel->t0();
    DataPlot(0, 1) = PosBlock(0);  // Position X
    DataPlot(0, 2) = PosBlock(1);  // Position Y
    DataPlot(0, 3) = PosBlock(2);  // Angle theta
    DataPlot(0, 4) = VelBlock(0);  // Velocity Vx
    DataPlot(0, 5) = VelBlock(1);  // Velocity Vy
    DataPlot(0, 6) = VelBlock(2);  // Angular velocity
    DataPlot(0, 7) = 0.0;          // Gap at first contact
    DataPlot(0, 8) = 0.0;          // Gap at second contact
    DataPlot(0, 9) = 0.0;          // Relative velocity at first contact
    DataPlot(0, 10) = 0.0;         // Relative velocity at second contact
    DataPlot(0, 11) = 0.0;         // Force at first contact
    DataPlot(0, 12) = 0.0;         // Force at second contact
    //----------------------------------- Simulation starts
    //----------------------------------------------------------
    std::cout << "====> Start computation ... \n\n";
    bool NSEvent = false;
    unsigned int NumberNSEvent = 0;
    double alpha_m, alpha_f, beta, gamma;
    unsigned int k = 1;
    auto start = std::chrono::system_clock::now();
    while (EDscheme->hasNextEvent() && (k < NpointSave + 1)) {
      if (IsTreatFirstSteps) {
        if (k == 1)  // first step
        {
          alpha_m = 0.0;
          alpha_f = 0.0;
          gamma = 1.0 / 2.0 + 1.0 / 10.0;
          beta = 1.0 / 4.0 * (gamma + 1.0 / 2.0) * (gamma + 1.0 / 2.0);
          NewMarkAlpha->setAlpha_m(alpha_m);
          NewMarkAlpha->setAlpha_f(alpha_f);
          NewMarkAlpha->setBeta(beta);
          NewMarkAlpha->setGamma(gamma);
        } else if (k == 2) {
          alpha_m = 0.0;
          alpha_f = -1.0 / 3.0;  // -1/3 <= alpha_f <= 0
          gamma = 1.0 / 2.0 - alpha_f;
          beta = 1.0 / 4.0 * (1.0 - alpha_f) * (1.0 - alpha_f);
          NewMarkAlpha->setAlpha_m(alpha_m);
          NewMarkAlpha->setAlpha_f(alpha_f);
          NewMarkAlpha->setBeta(beta);
          NewMarkAlpha->setGamma(gamma);
        } else {
          NewMarkAlpha->setParametersFromRho_infty(_rho);
        }
      }
      EDscheme->advanceToEvent();  // lead the simulation run from one event to the next
      auto GapCon1 = inter1->y(0);
      auto GapCon2 = inter2->y(0);
      auto VelCon1 = inter1->y(1);
      auto VelCon2 = inter2->y(1);

      auto LambdaCon1 = inter1->lambda(2);
      auto LambdaCon2 = inter2->lambda(2);
      //---------- detect the statue of the current event ------------------------------------
      if (eventsManager->nextEvent()->getType() ==
          siconos::simulation::EventType::NS)  // the current event is non-smooth
      {
        NSEvent = true;
      };
      EDscheme->processEvents();  // process the current event
      //------------------- get data at the beginning of non-smooth events
      //---------------------------
      if (NSEvent) {
        DataPlot(k, 0) = EDscheme->startingTime();  // instant at non-smooth event
        DataPlot(k, 1) = RockingBlock->qMemory().getSiconosVector(0)(0);         // Position X
        DataPlot(k, 2) = RockingBlock->qMemory().getSiconosVector(0)(1);         // Position Y
        DataPlot(k, 3) = RockingBlock->qMemory().getSiconosVector(0)(2);         // Angle theta
        DataPlot(k, 4) = RockingBlock->velocityMemory().getSiconosVector(0)(0);  // Velocity Vx
        DataPlot(k, 5) = RockingBlock->velocityMemory().getSiconosVector(0)(1);  // Velocity Vy
        DataPlot(k, 6) =
            RockingBlock->velocityMemory().getSiconosVector(0)(2);  // Angular velocity
        // EDscheme->update(1);
        k++;
        ++NumberNSEvent;

        NSEvent = false;  // The next event is maybe smooth
      };
      //-------------------- get data at smooth events or at the end of non-smooth events
      DataPlot(k, 0) = EDscheme->startingTime();
      DataPlot(k, 1) = PosBlock(0);        // Position X
      DataPlot(k, 2) = PosBlock(1);        // Position Y
      DataPlot(k, 3) = PosBlock(2);        // Position theta
      DataPlot(k, 4) = VelBlock(0);        // Velocity Vx
      DataPlot(k, 5) = VelBlock(1);        // Velocity Vy
      DataPlot(k, 6) = VelBlock(2);        // Velocity Vtheta
      DataPlot(k, 7) = (*GapCon1)(0);      // Gap at first contact
      DataPlot(k, 8) = (*GapCon2)(0);      // Gap at second contact
      DataPlot(k, 9) = (*VelCon1)(0);      // Relative velocity at first contact
      DataPlot(k, 10) = (*VelCon2)(0);     // Relative velocity at second contact
      DataPlot(k, 11) = (*LambdaCon1)(0);  // Force at first contact
      DataPlot(k, 12) = (*LambdaCon2)(0);  // Force at second contact
      // go to the next time step
      k++;

      // // Display information
      // std::cout << "********At the end of integation step***************"<< (k - 1) << "\n";
      // std::cout << "Information on Dynamical System\n";
      // std::cout << "Position: ";
      // siconos::algebra::print(*PosBlock);
      // std::cout << "Velocity: ";
      // siconos::algebra::print(*VelBlock);
      // std::cout << "Acceleration: ";
      // siconos::algebra::print(*AcceBlock);
      // std::cout << "Information on contacts\n";
      // for(std::tie(ui,uiend) = indexSet0->vertices(); ui!=uiend; ++ui)
      //   {
      //     auto inter = indexSet0->bundle(*ui);
      //     std::cout << "Contact number: " << inter->number() << "\n";
      //     std::cout << "Contact gap: ";
      //     siconos::algebra::print(*inter->y(0));
      //     std::cout << "Contact relative velocity: ";
      //     siconos::algebra::print(*inter->y(1));
      //     std::cout << "Contact Force: \n";
      //     siconos::algebra::print(*inter->lambda(2));
      //   }
    };
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write("RockingBlockED_NewMarkAlpha.dat", DataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    // Comparison with the reference, in two groups of columns.
    //
    // The generalized-alpha scheme (NewMarkAlphaOSI) recovers the acceleration, and the
    // contact forces (solution of the LCP at the acceleration level), from the positions,
    // through a = (q - q_pred) / (beta h^2): a round-off difference of one ulp on q gives
    // a difference of ~1e-16 / (beta h^2) on a and on the forces, which then lives in the
    // barely damped high-frequency mode of the scheme. The trajectory is not affected:
    // positions, velocities, gaps and relative velocities stay within 1e-11 of the
    // reference. Between macOS (arm64, clang, fused multiply-add) and Linux (x86_64, gcc),
    // the contact force differs by up to 1.6e-6.
    //  - kinematic columns (t, q, v, gaps, relative velocities): absolute tolerance 1e-9;
    //  - contact forces: tolerance 1e-5 * max |force|.
    double eps = 1e-9;
    if (siconos::algebra::io::compareRefFile(DataPlot, "RockingBlockED_NewMarkAlpha.ref", eps,
                                             {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}) > eps)
      return 1;
    double eps_forces = 1e-5 * std::max(DataPlot.col(11).cwiseAbs().maxCoeff(),
                                        DataPlot.col(12).cwiseAbs().maxCoeff());
    if (siconos::algebra::io::compareRefFile(DataPlot, "RockingBlockED_NewMarkAlpha.ref",
                                             eps_forces, {11, 12}) > eps_forces)
      return 1;
    return 0;
  } catch (...) {
    siconos::exception::process();
    return 1;
  }
}
