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
#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

constexpr double GGearth = 9.8100;

//---------------------------------------------------
double LengthBlock = 1.0;                       // Length of the rocking block
double HeightBlock = 0.5 * LengthBlock;         // Height of the rocking block
unsigned int Nfreedom = 3;                      // Number of degrees of freedom
unsigned int Ncontact = 2;                      // Number of contacts
double MassBlock = 1.0;                         // Mass of the rocking block
double PosXiniPointA = 0.0;                     // Initial coordinate X of the point A
double PosYiniPointA = 0.0;                     // Initial coordinate Y of the point A
double AngleThetaIni = std::numbers::pi / 3.0;  // Initial angle theta of the block
double VelXiniPointA = 0.0;                     // Initial relative velocity Vx of the point A
double VelYiniPointA = 0.0;                     // Initial relative velocity Vy of the point A
double RotVelBlockIni = 0.0;                    // Initial angular velocity of the block
// Parameters for impact model
double e1 = 0.8;                             // Restitution coefficient at the contact A
double e2 = 0.8;                             // Restitution coefficient at the contact B
double K1 = 1.0e10;                          // Stiffness at the contact A
double K2 = 1.0e10;                          // Stiffness at the contact B
double eta = 1.0;                            // Elasticity coefficient
double DelP = 5.0e-5;                        // Step size used for impact resolution
std::string TypeContactLaw = "BiStiffness";  // Type of compliance contact law
unsigned int Nstep_save_impact =
    100;                              // Number of steps every which we save data during impact
unsigned int Npoint_save_impact = 1;  // Number of points saved during impact
unsigned int step_start = 0;
unsigned int step_end = step_start + Npoint_save_impact * Nstep_save_impact;
unsigned int Nstep_max_impact =
    1000000;  // Number maximal of steps allowed for impact computation
std::string impact_data_name = "data_impact.dat";
bool _IsSaveDataImpact = false;
// Parameter for Event-Driven simulation
double TimeInitial = 0.0;       // Initial time of the simulation
double TimeFinal = 2.0;         // Final time of the simulation
double StepSize = 0.01;         // Time step size
unsigned int NpointSave = 252;  //
unsigned int SizeOutput = 13;   //
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
    InertiaBlock =
        (MassBlock / 12.0) * (pow(HeightBlock, 2) + pow(LengthBlock, 2));  // moment of inertia
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
    auto nslaw1 = std::make_shared<siconos::modeling::MultipleImpactNSL>(e1, K1, eta);
    auto nslaw2 = std::make_shared<siconos::modeling::MultipleImpactNSL>(e2, K2, eta);
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

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw1, relation1);
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

    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw2, relation2);
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
    auto OSI = std::make_shared<siconos::integrators::LsodarOSI>();
    // 3. Nonsmooth problem
    auto impact = std::make_shared<siconos::nonsmooth_formulations::MultipleImpact>(
        TypeContactLaw, DelP);
    auto multiple_impact =
        std::dynamic_pointer_cast<siconos::nonsmooth_formulations::MultipleImpact>(impact);
    multiple_impact->SetSaveData(_IsSaveDataImpact);
    multiple_impact->SetNameOutput(impact_data_name.c_str());
    multiple_impact->SetNstepSave(Nstep_save_impact);
    multiple_impact->SetNstepMax(Nstep_max_impact);
    multiple_impact->SetSizeDataSave(Npoint_save_impact);
    multiple_impact->SetStepMinMaxSave(step_start, step_end);
    auto acceleration = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    // 4. Simulation with (1), (2), (3)
    auto EDscheme =
        std::make_shared<siconos::simulation::EventDriven>(RoBlockModel, TimeDiscret);
    EDscheme->insertIntegrator(OSI);
    EDscheme->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_ED_IMPACT);
    EDscheme->insertNonSmoothProblem(acceleration,
                                     siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC);

    // bool check1 = EDscheme->hasOneStepNSProblem(impact);
    // bool check2 = EDscheme->hasOneStepNSProblem(acceleration);
    // std::cout << "Impact law included in the simulation: " << check1 << std::endl;
    // std::cout << "LCP at acceleration level included in the simulation: " << check2 <<
    // std::endl;
    //==================================================================================================================
    //                    V. Process the simulation
    //==================================================================================================================
    // -------------------------------- Simulation initialization
    // ------------------------------------------------------
    // auto lsodar = std::static_pointer_cast<LsodarOSI>(OSI);
    // lsodar->setMinMaxStepSizes(1.0e-3,1.0e-3);
    // lsodar->setTol(1,1.0e-3,1.0e-6);
    // lsodar->setMaxOrder(2, 2);

    auto eventsManager =
        EDscheme->eventsManager();  // ponters point to the "eventsManager" object
    auto PosBlock =
        RockingBlock->q_read();  // read-only view to the position vector of the rocking block
    auto VelBlock =
        RockingBlock->velocity_read();  // read-only view  to the velocity of the rocking block
    RockingBlock->initRhs(RoBlockModel->t0());
    auto AcceBlock =
        RockingBlock
            ->acceleration_read();  // read-only view  to the velocity of the rocking block
    auto GapCon1 = inter1->y(0);
    auto GapCon2 = inter2->y(0);
    auto VelCon1 = inter1->y(1);
    auto VelCon2 = inter2->y(1);
    auto indexSet0 = RoBlockModel->topology()->indexSet(0);
    std::cout << "Size of IndexSet0: " << indexSet0->size() << std::endl;

    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    //-------------------- Save the output during simulation
    //---------------------------------------------------------
    unsigned int NpointSave = 252;  //
    unsigned int SizeOutput = 13;   //
    Matrix DataPlot(NpointSave, SizeOutput);
    //------------- At the initial time
    //-----------------------------------------------------------------------------
    DataPlot(0, 0) = RoBlockModel->t0();
    DataPlot(0, 1) = PosBlock(0);     // Position X
    DataPlot(0, 2) = PosBlock(1);     // Position Y
    DataPlot(0, 3) = PosBlock(2);     // Angle theta
    DataPlot(0, 4) = VelBlock(0);     // Velocity Vx
    DataPlot(0, 5) = VelBlock(1);     // Velocity Vy
    DataPlot(0, 6) = VelBlock(2);     // Angular velocity
    DataPlot(0, 7) = (*GapCon1)(0);   // Gap at first contact
    DataPlot(0, 8) = (*GapCon2)(0);   // Gap at second contact
    DataPlot(0, 9) = (*VelCon1)(0);   // Relative velocity at first contact
    DataPlot(0, 10) = (*VelCon2)(0);  // Relative velocity at second contact
    DataPlot(0, 11) = 0.0;            //(*LambdaCon1)(0); // Force at first contact
    DataPlot(0, 12) = 0.0;            //(*LambdaCon2)(0); // Force at second contact
    //----------------------------------- Simulation starts
    //----------------------------------------------------------
    std::cout << "====> Start computation ... \n\n";
    bool NSEvent = false;
    unsigned int NumberNSEvent = 0;
    unsigned int k = 1;
    auto start = std::chrono::system_clock::now();
    while (EDscheme->hasNextEvent() && (k < NpointSave)) {
      EDscheme->advanceToEvent();  // lead the simulation run from one event to the next
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
        DataPlot(k, 1) = PosBlock(0);               // Position X
        DataPlot(k, 2) = PosBlock(1);               // Position Y
        DataPlot(k, 3) = PosBlock(2);               // Position theta
        DataPlot(k, 4) = RockingBlock->velocityMemory().getSiconosVector(1)(0);  // Velocity Vx
        DataPlot(k, 5) = RockingBlock->velocityMemory().getSiconosVector(1)(1);  // Velocity Vy
        DataPlot(k, 6) =
            RockingBlock->velocityMemory().getSiconosVector(1)(2);  // Angular velocity
        DataPlot(k, 7) = (*GapCon1)(0);                             // Gap at first contact
        DataPlot(k, 8) = (*GapCon2)(0);                             // Gap at second contact
        DataPlot(k, 9) = (*VelCon1)(0);   // Relative velocity at first contact
        DataPlot(k, 10) = (*VelCon2)(0);  // Relative velocity at second contact
        // EDscheme->update(1);
        k++;
        ++NumberNSEvent;

        NSEvent = false;  // The next event is maybe smooth
      };
      //-------------------- get data at smooth events or at the end of non-smooth events
      //---------------
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
    siconos::algebra::io::write("RockingBlockED_Lsodar_LZBimpact.dat", DataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(
             DataPlot, "RockingBlockED_Lsodar_LZBimpact.ref", eps)) > eps)
      return 1;
    return 0;
  } catch (...) {
    siconos::exception::process();
    return 1;
  }
}
