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

/*!\file BouncingBallED.cpp
  \brief \ref EMBouncingBall - C++ input file, Event-Driven version - V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with an Event-Driven scheme.
*/

#include <SiconosKernel.hpp>
#include <chrono>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

constexpr double PI = 3.14159265;
constexpr double g = 0.0;  // Gravity

using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDofBall = 1;     // degrees of freedom of ball 1
    double Height = 0.2;  // Distance between impactor balls and monodisperse balls
    double V_impact = 1.0;
    // Balls in the tapered chain
    unsigned int NumberBalls = 10;  // Number
    double R_ball = 0.01;           // Ball radius
    // Material properties of balls
    double mass_density = 7780;  // mass density
    double CoefRes = 1.0;        // Restitution coefficient
    double YoungBall = 203.0e9;  // Young modulus of the balls
    double PoissonBall = 0.3;    // Poison coefficient of the balls
    double PowCompLaw = 1.5;     // Power of the compliance law: 1.0 for linear contact and 3/2
                                 // for the Hertzian contact
    std::string TypeContactLaw = "BiStiffness";  // Type of compliance contact law
    // Parameters for the global simulation
    double t0 = 0;                  // initial computation time
    double T = 0.6;                 // final computation time
    double h = 0.001;               // time step
    unsigned int Npointsave = 610;  // Number of data points to be saved
    // For impact computation
    double DelPest = 1.0e-6;  // Step size estimated for multiple impacts computation
    unsigned int Nstep_save_impact =
        100;  // Number of steps every which we save data during impact
    unsigned int Step_begin = 1;
    unsigned int Npoint_save_impact = 100;  // Number of points saved during impac
    unsigned int Step_end = Step_begin + Npoint_save_impact * Nstep_save_impact;
    unsigned int Nstep_max_impact =
        10000000;  // Number maximal of steps allowed for impact computation
    std::string impact_data_name = "data_impact.dat";
    bool _IsSaveDataImpact = false;
    //---------------------------------------
    // ---- Configuration of chaines
    //--------------------------------------
    //************* Balls ******************
    unsigned int NumberContacts = NumberBalls - 1;  // Number of contacts
    //(1) Radius of balls
    auto RadiusBalls = std::make_shared<Vector>(NumberBalls);
    for (unsigned int k = 0; k < NumberBalls; ++k) {
      (*RadiusBalls)(k) = R_ball;
    }
    // (2) Mass of balls
    auto MassBalls = std::make_shared<Vector>(NumberBalls);
    for (unsigned int id = 0; id < NumberBalls; ++id) {
      (*MassBalls)(id) = (4.0 / 3.0) * PI * pow((*RadiusBalls)(id), 3) * mass_density;
    }
    // (3) Initial position of balls
    // For the impactor balls
    auto InitPosBalls = std::make_shared<Vector>(NumberBalls);
    (*InitPosBalls)(0) = 0.0;
    (*InitPosBalls)(1) = (*RadiusBalls)(0) + Height + (*RadiusBalls)(1);
    for (unsigned int j = 2; j < NumberBalls; ++j) {
      (*InitPosBalls)(j) = (*InitPosBalls)(j - 1) + (*RadiusBalls)(j - 1) + (*RadiusBalls)(j);
    }
    // (4) Initial velocity of balls
    auto InitVelBalls = std::make_shared<Vector>(NumberBalls);
    (*InitVelBalls)(0) = V_impact;
    for (unsigned int i = 1; i < NumberBalls; ++i) {
      (*InitVelBalls)(i) = 0.0;
    }
    //****************** Contacts ******************
    // (1) Restitution coefficient at contacts
    auto ResCofContacts = std::make_shared<Vector>(NumberContacts);
    auto ElasCofContacts = std::make_shared<Vector>(NumberContacts);
    for (unsigned int id = 0; id < NumberContacts; ++id) {
      (*ResCofContacts)(id) = CoefRes;
      (*ElasCofContacts)(id) = PowCompLaw;
    }
    // (2) Stiffness at contacts
    auto StiffContacts = std::make_shared<Vector>(NumberContacts);
    double Rmoy, Emoy;
    for (unsigned int id = 0; id < NumberContacts; ++id) {
      Emoy = (2.0 / 3.0) * (YoungBall / (1.0 - pow(PoissonBall, 2)));
      Rmoy = ((*RadiusBalls)(id) * (*RadiusBalls)(id + 1)) /
             ((*RadiusBalls)(id) + (*RadiusBalls)(id + 1));
      (*StiffContacts)(id) = pow(Rmoy, 0.5) * Emoy;
    }
    // // Display and save the configuration of the chain simulated
    // cout << "Configuation of ball chains\n";
    // cout.precision(15);
    // cout << "Radius of balls: \n";
    // siconos::algebra::print(*RadiusBalls);
    // cout << "Mass of balls: \n";
    // siconos::algebra::print(*MassBalls);
    // cout << "Initial position of balls: \n";
    // siconos::algebra::print(*InitPosBalls);
    // cout << "Initial velocity of balls: \n";
    // siconos::algebra::print(*InitVelBalls);
    // cout << "Restitution coefficient at contacts:\n";
    // siconos::algebra::print(*ResCofContacts);
    // cout<< "Stiffness at contacts: \n";
    // siconos::algebra::print(*StiffContacts);
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ...\n\n";
    // -------------
    // --- Model ---
    // -------------
    auto BallChain = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::LsodarOSI>();

    std::vector<std::shared_ptr<siconos::modeling::LagrangianLinearTIDS>> VecOfallDS;
    double _Rball, _massBall, _Pos0Ball, _Vel0Ball;
    for (unsigned int i = 0; i < NumberBalls; ++i) {
      _Rball = (*RadiusBalls)(i);      // radius of the ball
      _massBall = (*MassBalls)(i);     // mass of the ball
      _Pos0Ball = (*InitPosBalls)(i);  // initial position of the ball
      _Vel0Ball = (*InitVelBalls)(i);  // initial velocity of the ball
      // Declaration of the DS in Siconos
      auto MassBall = std::make_shared<Matrix>(nDofBall, nDofBall);
      (*MassBall)(0, 0) = _massBall;
      // -- Initial positions and velocities --
      auto q0Ball = std::make_shared<Vector>(nDofBall);
      auto v0Ball = std::make_shared<Vector>(nDofBall);
      (*q0Ball)(0) = _Pos0Ball;
      (*v0Ball)(0) = _Vel0Ball;
      // -- The dynamical system --
      auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
          *q0Ball, *v0Ball, *MassBall, siconos::algebra::copy_t);
      // -- Set external forces (weight1) --

      Vector FextBall{nDofBall};
      FextBall.setZero();
      FextBall(0) = -_massBall * g;
      ball->setConstantFext(FextBall, siconos::algebra::alias_t);

      //
      VecOfallDS.push_back(ball);
      BallChain->insertDynamicalSystem(ball);
    }
    // --------------------
    // --- Interactions ---
    // --------------------
    auto H = std::make_shared<Matrix>(1, (nDofBall + nDofBall));
    double ResCoef, Stiff, ElasPow;
    (*H)(0, 0) = -1.0;
    (*H)(0, 1) = 1.0;
    auto E = std::make_shared<Vector>(1);

    for (unsigned int j = 0; j < NumberContacts; ++j) {
      ResCoef = (*ResCofContacts)(j);
      Stiff = (*StiffContacts)(j);
      ElasPow = (*ElasCofContacts)(j);
      (*E)(0) = -1.0 * ((*RadiusBalls)(j) + (*RadiusBalls)(j + 1));
      auto nslaw =
          std::make_shared<siconos::modeling::MultipleImpactNSL>(ResCoef, Stiff, ElasPow);
      auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H, *E);
      auto interaction = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);
      BallChain->link(interaction, VecOfallDS[j], VecOfallDS[j + 1]);
    }

    // ----------------
    // --- Simulation ---
    // ----------------
    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);
    // -- (3) Non smooth problem --
    auto impact = std::make_shared<siconos::nonsmooth_formulations::MultipleImpact>(
        TypeContactLaw, DelPest);
    // auto impact=
    // std::make_shared<siconos::nonsmooth_formulations::MultipleImpact>(TypeContactLaw,NestImpact);
    auto multiple_impact =
        std::dynamic_pointer_cast<siconos::nonsmooth_formulations::MultipleImpact>(impact);
    multiple_impact->SetSaveData(_IsSaveDataImpact);
    multiple_impact->SetNameOutput(impact_data_name.c_str());
    multiple_impact->SetNstepSave(Nstep_save_impact);
    multiple_impact->SetNstepMax(Nstep_max_impact);
    multiple_impact->SetStepMinMaxSave(Step_begin, Step_end);
    multiple_impact->SetSizeDataSave(Npoint_save_impact);
    auto acceleration = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::EventDriven>(BallChain, t);
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_ED_IMPACT);
    s->insertNonSmoothProblem(acceleration, siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC);

    // =========================== End of model definition ===========================
    //----------------------------------- Initialization-------------------------------
    s->setPrintStat(true);
    auto DSG0 = BallChain->topology()->dSG(0);
    auto IndexSet0 = BallChain->topology()->indexSet(0);
    // // Display topology of the system
    // cout << "Number of vectices of IndexSet0: " << IndexSet0->size() << endl;
    // cout << "Number of vectices of DSG0: " << DSG0->size() << endl;
    //
    auto eventsManager = s->eventsManager();
    // ================================= Computation =================================
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 2 * NumberBalls + 1;
    Matrix dataPlot(Npointsave, outputSize);

    // --- Time loop ---
    cout << "====> Start computation ... \n\n";
    // ==== Simulation loop - Writing without explicit event handling =====
    bool nonSmooth = false;
    unsigned int NumberOfEvents = 0;
    unsigned int NumberOfNSEvents = 0;
    unsigned int k = 0;
    siconos::graphs::DynamicalSystemsGraph::VIterator ui, uiend;
    auto start = std::chrono::system_clock::now();

    //====================================================================
    while ((k < Npointsave) & (s->hasNextEvent())) {
      dataPlot(k, 0) = s->startingTime();
      // Save state of the balls
      unsigned int col_pos = 1;
      unsigned int col_vel = NumberBalls + 1;
      for (boost::tie(ui, uiend) = DSG0->vertices(); ui != uiend; ++ui) {
        auto ds = DSG0->bundle(*ui);
        auto lag_ds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds);
        auto q = lag_ds->q();
        auto v = lag_ds->velocity();
        dataPlot(k, col_pos) = (*q)(0);
        dataPlot(k, col_vel) = (*v)(0);
        col_pos++;
        col_vel++;
      }
      ++k;
      s->advanceToEvent();  // run simulation from one event to the next
      if (eventsManager->nextEvent()->getType() == siconos::simulation::EventType::NS) {
        nonSmooth = true;
      };
      //
      s->processEvents();  // process events
      if (nonSmooth) {
        // siconos::algebra::print(*multiple_impact);
        dataPlot(k, 0) = s->startingTime();
        // Save state of the balls
        unsigned int col_pos = 1;
        unsigned int col_vel = NumberBalls + 1;
        for (boost::tie(ui, uiend) = DSG0->vertices(); ui != uiend; ++ui) {
          auto ds = DSG0->bundle(*ui);
          auto lag_ds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds);
          const auto& q = lag_ds->qMemory().getSiconosVector(1);
          const auto& v = lag_ds->velocityMemory().getSiconosVector(1);
          dataPlot(k, col_pos) = q(0);
          dataPlot(k, col_vel) = v(0);
          col_pos++;
          col_vel++;
        }
        nonSmooth = false;
        ++NumberOfNSEvents;
        ++NumberOfEvents;

        ++k;
      }
      // --- Get values to be plotted ---
      ++NumberOfEvents;
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("MonodisperseChainOfBalls-LZBModel.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "MonodisperseChainOfBalls-LZBModel.ref", eps)) > eps)
      return 1;
    return 0;

  } catch (...) {
    siconos::exception::process();
    return 1;
  }
}
