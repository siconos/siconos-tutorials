#ifndef USER_PARAMS
#define USER_PARAMS

namespace user_defined_ref {

// User-defined main parameters
inline constexpr int nDof = 100;     // degrees of freedom for the beam
inline constexpr double t0 = 1e-8;   // initial computation time
inline constexpr double T = 0.0015;  // final computation time
inline constexpr double h = 2e-6;    // time step
inline constexpr double position_init =
    0.0000503;  // initial position (contact at 251.5 h: not on a time step)
inline constexpr double velocity_init = -.1;        // initial velocity
inline constexpr double epsilon = 0.0;              // 1e-1;
inline constexpr double theta = 1 / 2.0 + epsilon;  // theta for MoreauJeanOSI integrator
// theta = 1.0;
inline constexpr double E = 210e9;     // young Modulus
inline constexpr double S = 0.000314;  //  Beam Section 1 cm  for the diameter
// S=0.1;
inline constexpr double L = 1.0;       // length of the  beam
inline constexpr double rho = 7800.0;  // specific mass
// rho=1.0;
// double g = 9.81;  // Gravity
inline constexpr double g = 0.0;
}  // namespace user_defined_ref

namespace user_defined {

// User-defined main parameters
inline constexpr int nDof = 10;      // degrees of freedom for the beam
inline constexpr double t0 = 1e-8;   // initial computation time
inline constexpr double T = 0.0015;  // final computation time
inline constexpr double h = 1e-7;    // time step
inline constexpr double position_init =
    0.0000503;  // initial position (contact at 251.5 h: not on a time step)
inline constexpr double velocity_init = -.1;        // initial velocity
inline constexpr double epsilon = 0.0;              // 1e-1;
inline constexpr double theta = 1 / 2.0 + epsilon;  // theta for MoreauJeanOSI integrator
// theta = 1.0;
inline constexpr double E = 210e9;     // young Modulus
inline constexpr double S = 0.000314;  //  Beam Section 1 cm  for the diameter
// S=0.1;
inline constexpr double L = 1.0;       // length of the  beam
inline constexpr double rho = 7800.0;  // specific mass
// rho=1.0;
// inline constexpr double g = 9.81; // Gravity
inline constexpr double g = 0.0;
}  // namespace user_defined
#endif
