// User-defined main parameters
unsigned int nElement = 100;  // degrees of freedom for the beam
double t0 = 1e-8;               // initial computation time
double T = .1;                  // final computation time
double h = 1e-3;                // time step
double position_init = 0.000001;   // initial position of the ball
double velocity_init =  -2.;    // initial velocity of the ball
double epsilon = 0.0;           //1e-1;
double theta = 1/2.0 + epsilon; // theta for MoreauJeanOSI integrator
//theta = 0.55;
double E = 210e9;               // young Modulus
double R = 0.01;                // 1 cm  for the radius
double I = 3.14/4.*R*R*R*R ;    // Moment of inertia
double S = 3.14* R*R;           // Beam Section 
double L = 1.0;                 // length of the  beam
double rho = 7800.0 ;           // specific mass
//rho=1.0;
double g = 9.81;                // Gravity
//g=0.0;

