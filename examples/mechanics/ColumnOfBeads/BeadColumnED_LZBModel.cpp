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
constexpr double g = 9.81;  // Gravity

using namespace std;

int main(int argc, char* argv[]) {
  try {
    // ================= Creation of the model =======================

    // User-defined main parameters
    int nDofBall = 1;      // degrees of freedom of ball 1
    double Height = 0.05;  // Distance between impactor balls and monodisperse balls
    // Balls in the tapered chain
    size_t NumberBalls = 10;      // Number
    double R_base_taper = 0.005;  // Base radii of the tapered chain
    double q_taper = 0.0;         // Tapering factor
    // Material properties of balls
    double mass_density = 7780;  // mass density
    double Res_BallBall = 1.0;   // Restitution coefficient at contacts ball-ball
    double Res_BallWall = 1.0;   // Restitution coefficient at contacts ball-wall
    double YoungBall = 203.0e9;  // Young modulus of the balls
    double PoissonBall = 0.3;    // Poison coefficient of the balls
    double YoungWall = 203.0e9;  // Young modulus of the wall
    double PoissonWall = 0.3;    // Poison coefficient of the wall
    double PowCompLaw = 1.5;     // Power of the compliance law: 1.0 for linear contact and 3/2
                                 // for the Hertzian contact
    std::string TypeContactLaw = "BiStiffness";  // Type of compliance contact law
    // Parameters for the global simulation
    double t0 = 0;     // initial computation time
    double T = 0.5;    // final computation time
    double h = 0.001;  // time step
    // For impact computation
    double DelPest = 1.0e-6;  // Step size estimated for multiple impacts computation
    unsigned int Nstep_save_impact =
        100;  // Number of steps every which we save data during impact
    unsigned int Npoint_save_impact = 2000;  // Number of points saved during impact
    unsigned int step_start = 0;
    unsigned int step_end = step_start + Npoint_save_impact * Nstep_save_impact;
    unsigned int Nstep_max_impact =
        1000000;  // Number maximal of steps allowed for impact computation
    std::string impact_data_name = "data_impact.dat";
    bool _IsSaveDataImpact = false;
    //---------------------------------------
    // ---- Configuration of chaines
    //--------------------------------------
    //************* Balls ******************
    auto NumberContacts = NumberBalls;  // Number of contacts
    //(1) Radius of balls
    Vector RadiusBalls{NumberBalls};
    RadiusBalls =
        Eigen::VectorXd::NullaryExpr(RadiusBalls.size(), [q_taper, R_base_taper](int i) {
          return (pow((1.0 - q_taper), i + 1)) * R_base_taper;
        });

    // (2) Mass of balls
    Vector MassBalls{NumberBalls};
    MassBalls =
        Eigen::VectorXd::NullaryExpr(MassBalls.size(), [RadiusBalls, mass_density](int i) {
          return (4.0 / 3.0) * PI * pow(RadiusBalls(i), 3) * mass_density;
        });

    //****************** Contacts ******************
    // (1) Restitution coefficient at contacts
    Vector ResCofContacts{NumberContacts};
    Vector ElasCofContacts{NumberContacts};
    ResCofContacts.setConstant(Res_BallWall);
    ElasCofContacts.setConstant(PowCompLaw);

    // (2) Stiffness at contacts
    Vector StiffContacts{NumberContacts};
    double Rmoy, Emoy;
    for (unsigned int id = 0; id < NumberContacts; ++id) {
      // for ball-wall contact
      if (id == 0) {
        Emoy = (4.0 / 3.0) *
               ((YoungBall * YoungWall) / ((1.0 - pow(PoissonBall, 2)) * YoungWall +
                                           (1.0 - pow(PoissonWall, 2)) * YoungBall));
        Rmoy = RadiusBalls(0);
      }
      // Ball-ball contact
      else {
        Emoy = (2.0 / 3.0) * (YoungBall / (1.0 - pow(PoissonBall, 2)));
        Rmoy =
            (RadiusBalls(id - 1) * RadiusBalls(id)) / (RadiusBalls(id - 1) + RadiusBalls(id));
      }
      StiffContacts(id) = pow(Rmoy, 0.5) * Emoy;
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
    // ----------------
    // --- Simulation ---
    // ----------------
    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::LsodarOSI>();

    std::vector<std::shared_ptr<siconos::modeling::LagrangianLinearTIDS>> VecOfallDS(
        NumberBalls);

    Vector InitPosBalls{NumberBalls};
    InitPosBalls(0) = Height + RadiusBalls(0);
    for (size_t j = 1; j < NumberBalls; ++j) {
      InitPosBalls(j) = InitPosBalls(j - 1) + RadiusBalls(j - 1) + RadiusBalls(j);
    }
    std::vector<Vector> FextBall(NumberBalls, Vector::Zero(nDofBall));
    std::vector<Matrix> MassBall(NumberBalls, Matrix::Zero(nDofBall, nDofBall));
    std::vector<Vector> q0Ball(NumberBalls, Vector::Zero(nDofBall));
    std::vector<Vector> vel0Ball(NumberBalls, Vector::Zero(nDofBall));

    for (auto id = 0; id < NumberBalls; ++id) {
      // -- The dynamical system --
      FextBall[id](0) = -MassBalls(id) * g;
      MassBall[id](0, 0) = MassBalls(id);
      q0Ball[id](0) = InitPosBalls(id);

      VecOfallDS[id] = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
          q0Ball[id], vel0Ball[id], MassBall[id], siconos::algebra::alias_t);
      VecOfallDS[id]->setConstantFext(FextBall[id], siconos::algebra::alias_t);
      BallChain->insertDynamicalSystem(VecOfallDS[id]);
    }

    // --------------------
    // --- Interactions ---
    // --------------------

    std::vector<std::shared_ptr<siconos::modeling::Interaction>> interactions(NumberContacts);
    std::vector<std::shared_ptr<siconos::modeling::MultipleImpactNSL>> nslaws(NumberContacts);
    std::vector<std::shared_ptr<siconos::modeling::LagrangianLinearTIR>> relations(
        NumberContacts);
    std::vector<Vector> E(NumberContacts, Vector::Zero(1));
    // contact wall-ball - id = 0
    Matrix Hwall_ball{1, nDofBall};
    Hwall_ball.setZero();
    Hwall_ball(0, 0) = 1.;
    E[0] << -RadiusBalls(0);
    nslaws[0] = std::make_shared<siconos::modeling::MultipleImpactNSL>(
        ResCofContacts(0), StiffContacts(0), ElasCofContacts(0));
    relations[0] = std::make_shared<siconos::modeling::LagrangianLinearTIR>(Hwall_ball, E[0]);
    interactions[0] =
        std::make_shared<siconos::modeling::Interaction>(nslaws[0], relations[0]);
    BallChain->link(interactions[0], VecOfallDS[0]);

    Matrix H{1, 2 * nDofBall};
    H.setZero();
    H(0, 0) = -1.0;
    H(0, 1) = 1.0;

    for (auto id = 1; id < NumberContacts; ++id) {
      auto ResCoef = ResCofContacts(id);
      auto Stiff = StiffContacts(id);
      auto ElasPow = ElasCofContacts(id);
      E[id](0) = -(RadiusBalls(id - 1) + RadiusBalls(id));

      nslaws[id] =
          std::make_shared<siconos::modeling::MultipleImpactNSL>(ResCoef, Stiff, ElasPow);
      relations[id] = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H, E[id]);
      interactions[id] =
          std::make_shared<siconos::modeling::Interaction>(nslaws[id], relations[id]);
      // For ball-ball contact
      BallChain->link(interactions[id], VecOfallDS[id - 1], VecOfallDS[id]);
    }
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
    multiple_impact->SetSizeDataSave(Npoint_save_impact);
    multiple_impact->SetStepMinMaxSave(step_start, step_end);
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
    unsigned int Npointsave = 778;  // Number of data points to be saved
    Matrix dataPlot(Npointsave, outputSize);

    // --- Time loop ---
    cout << "====> Start computation ...\n\n ";
    // ==== Simulation loop - Writing without explicit event handling =====
    bool nonSmooth = false;
    unsigned int NumberOfEvents = 0;
    unsigned int NumberOfNSEvents = 0;
    siconos::algebra::Index k = 0;
    siconos::graphs::DynamicalSystemsGraph::VIterator ui, uiend;
    auto start = std::chrono::system_clock::now();

    size_t col_pos = 1;
    size_t col_vel = NumberBalls + 1;
    dataPlot(k, 0) = s->startingTime();

    for (auto id = 1; id < NumberContacts; ++id) {
      auto ds = VecOfallDS[id];
      auto lag_ds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds);
      auto q = lag_ds->q_read();
      auto v = lag_ds->velocity_read();
      dataPlot(k, col_pos) = q(0);
      dataPlot(k, col_vel) = v(0);
      col_pos++;
      col_vel++;
    }

    //====================================================================
    while ((k < Npointsave) & (s->hasNextEvent())) {
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
        col_pos = 1;
        col_vel = NumberBalls + 1;
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

      col_pos = 1;
      col_vel = NumberBalls + 1;
      dataPlot(k, 0) = s->startingTime();
      for (boost::tie(ui, uiend) = DSG0->vertices(); ui != uiend; ++ui) {
        auto ds = DSG0->bundle(*ui);
        auto lag_ds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds);
        auto q = lag_ds->q_read();
        auto v = lag_ds->velocity_read();
        dataPlot(k, col_pos) = q(0);
        dataPlot(k, col_vel) = v(0);
        col_pos++;
        col_vel++;
      }
      ++k;
      // --- Get values to be plotted ---
      ++NumberOfEvents;
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms\n";
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    siconos::algebra::io::write("BeadColumnED_LZBModel.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);

    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(dataPlot, "BeadColumnED_LZBModel.ref",
                                                      eps)) > eps)
      return 1;
    return 0;
  } catch (...) {
    siconos::exception::process();
    return 1;
  }
}
