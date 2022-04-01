#ifndef SCCONST_H
#define SCCONST_H

namespace parameters {
// parameters according to Table 1
// geometrical characteristics
double l1 = 0.1530;
double l2 = 0.3060;
double a = 0.05;
double b = 0.025;
double c = 0.001;
double d = 2. * (b + c);

// inertial properties
double m1 = 0.038;
double m2 = 0.038;
double m3 = 0.0760;
double J1 = 7.4e-5;
double J2 = 5.9e-4;
double J3 = 2.7e-6;

// force elements
double gravity = 9.81;

}; // namespace parameters

#endif
