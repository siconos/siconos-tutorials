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

// This is the program to simulate the dynamic of a rocking block by using the Siconos platform
//==================================================================================================================
#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>
#include <string>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

constexpr double GGearth = 9.8100;

//---------------------------------------------------
double LengthBlock = 1.;   // Length of the rocking block
double HeightBlock = 1.5;  // Height of the rocking block
double LBlock = 0.2;
double HBlock = 0.1;

unsigned int Nfreedom = 3;                       // Number of degrees of freedom
unsigned int Ncontact = 2;                       // Number of contacts
double MassBlock = 1.0;                          // Mass of the rocking block
double PosXiniPointA = 0.0;                      // Initial coordinate X of the point A
double PosYiniPointA = 0.5;                      // Initial coordinate Y of the point A
double AngleThetaIni = std::numbers::pi / 10.0;  // Initial angle theta of the block
double VelXiniPointA = 0.0;                      // Initial relative velocity Vx of the point A
double VelYiniPointA = 0.0;                      // Initial relative velocity Vy of the point A
double RotVelBlockIni = 0.0;                     // Initial angular velocity of the block
double e = 0.5;                                  // Restitution coefficient
double TimeInitial = 0.0;                        // Initial time of the simulation
double TimeFinal = 2.0;                          // Final time of the simulation
double StepSize = 0.01;                          // Time step size
double criterion = 0.05;
unsigned int maxIter = 20000;
//==========================================================================================================
//                                             Main function
//==========================================================================================================
// Simulation with a given update rule of the positions in the projection loop:
// - ProjectionUpdate::incremental (default): q_{j+1} = q_j + G(q_j)^T lambda_j;
// - ProjectionUpdate::from_origin: q_{j+1} = a + G(q_j)^T mu_j, exact projection of the
//   positions a given by the velocity step (as in Acary, CMAME 2013 and the Scilab code Hots).
// See docs/devel/PROJECTION_SCHEMES.md in siconos. Results in name.dat, compared with
// name.ref.
int simulate(siconos::simulation::ProjectionUpdate update, const std::string& name) {
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
    // (*PosIniBlock)(0) = PosXiniPointA + 0.5 * LengthBlock * cos(AngleThetaIni) -
    //                     0.5 * HeightBlock * sin(AngleThetaIni);
    // (*PosIniBlock)(1) = PosYiniPointA + 0.5 * LengthBlock * sin(AngleThetaIni) +
    //                     0.5 * HeightBlock * cos(AngleThetaIni);
    // (*PosIniBlock)(2) = AngleThetaIni;

    (*PosIniBlock)(0) = 0.0;
    (*PosIniBlock)(1) = 1.;
    (*PosIniBlock)(2) = 0.2;

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
    /*
    (*VelIniBlock)(0) = 0.0;
    (*VelIniBlock)(1) = 0.0;
    (*VelIniBlock)(2) = 0.0;
    */
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
    // Impact law
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    // Interaction at contact point 1
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation1->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = q(1) - 0.5 * LBlock * sin(q(2)) - 0.5 * HBlock * cos(q(2));
    });

    relation1->setComputeHFunction([](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::MapType> result) {
      result(0, 0) = 0.0;
      result(0, 1) = 1.0;
      result(0, 2) = -0.5 * LBlock * cos(q(2)) + 0.5 * HBlock * sin(q(2));
    });

    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation1);
    // Interaction at contact point 2
    auto relation2 = std::make_shared<siconos::modeling::LagrangianScleronomousR>();
    relation2->setComputehFunction([](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::SiconosVector> y) {
      y(0) = q(1) + 0.5 * LBlock * sin(q(2)) - 0.5 * HBlock * cos(q(2));
    });

    relation2->setComputeHFunction([](const siconos::algebra::BlockVector& q,
                                      Eigen::Ref<siconos::algebra::MapType> result) {
      result(0, 0) = 0.0;
      result(0, 1) = 1.0;
      result(0, 2) = 0.5 * LBlock * cos(q(2)) + 0.5 * HBlock * sin(q(2));
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
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanDirectProjectionOSI>(0.5001);
    // 3. Nonsmooth problem
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto impact_pos =
        std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>(
            SICONOS_MLCP_ENUM);
    // 4. Simulation with (1), (2), (3)
    auto TSscheme = std::make_shared<siconos::simulation::TimeSteppingDirectProjection>(
        RoBlockModel, TimeDiscret, OSI, impact, impact_pos, 0);
    TSscheme->setProjectionUpdate(update);

    //==================================================================================================================
    //                    V. Process the simulation
    //==================================================================================================================
    // -------------------------------- Simulation initialization
    // ------------------------------------------------------
    auto PosBlock =
        RockingBlock->q_read();  // read-only view to the position vector of the rocking block
    auto VelBlock =
        RockingBlock->velocity_read();  // read-only view  to the velocity of the rocking block
    //-------------------- Save the output during simulation
    //---------------------------------------------------------
    unsigned int NpointSave = 200;  //
    unsigned int SizeOutput = 9;    //
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

    Vector tmp{Nfreedom};
    tmp = *mass * VelBlock;
    double kineticEnergy = 0.5 * VelBlock.dot(tmp);
    DataPlot(0, 7) = kineticEnergy;

    Vector PosRef{Nfreedom};
    PosRef.setZero();
    PosRef(1) = HeightBlock / 2.0;
    double potentialEnergy = -(PosBlock - PosRef).dot(ExternalForces);
    DataPlot(0, 8) = potentialEnergy;

    //----------------------------------- Simulation starts
    //----------------------------------------------------------
    std::cout << "====> Start computation ... \n\n";
    unsigned int k = 1;
    auto start = std::chrono::system_clock::now();
    while (k < NpointSave) {
      TSscheme->computeOneStep();
      DataPlot(k, 0) = TSscheme->nextTime();
      DataPlot(k, 1) = PosBlock(0);  // Position X
      DataPlot(k, 2) = PosBlock(1);  // Position Y
      DataPlot(k, 3) = PosBlock(2);  // Position theta
      DataPlot(k, 4) = VelBlock(0);  // Velocity Vx
      DataPlot(k, 5) = VelBlock(1);  // Velocity Vy
      DataPlot(k, 6) = VelBlock(2);  // Velocity Vtheta

      tmp = *mass * VelBlock;
      kineticEnergy = 0.5 * VelBlock.dot(tmp);
      DataPlot(k, 7) = kineticEnergy;

      potentialEnergy = -(PosBlock - PosRef).dot(ExternalForces);
      DataPlot(k, 8) = potentialEnergy;
      // go to the next time step
      k++;

      TSscheme->nextStep();
    };
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    siconos::algebra::io::write(name + ".dat", DataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(DataPlot, name + ".ref", eps)) > eps)
      return 1;
    return 0;
  } catch (...) {
    siconos::exception::process();
    return 1;
  }
}

int main() {
  int err = simulate(siconos::simulation::ProjectionUpdate::incremental, "RockingBlockTSProj");
  err += simulate(siconos::simulation::ProjectionUpdate::from_origin,
                  "RockingBlockTSProj-FromOrigin");
  return err;
}
