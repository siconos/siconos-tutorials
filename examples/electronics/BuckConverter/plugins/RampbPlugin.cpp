#include <stdio.h>
#include <math.h>
#include <iostream>
//#include <RuntimeException.h>

using namespace std;

// function to compute b in x' = A x + b(t,z) + r
extern "C" void Rampb(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double *z)
{

    double epsitime = 1e-15;
    double VlowRamp = z[0];
    double VhighRamp = z[1];
    double RampTD = z[2];
    double RampTR = z[3];
    double RampTF = z[4];
    double RampPW = z[5];
    double RampPER = z[6];
    double SlopeRiseRamp = (VhighRamp - VlowRamp)/RampTR;
    double SlopeFallRamp = (VhighRamp - VlowRamp)/RampTF;
    double phaseRamp;

    double Vref = z[7];
    double VrefSettlingTime = z[8];

    double phasePulse;
    double UPtr[2];

    phaseRamp = fmod(time,RampPER);

    if (phaseRamp < RampTD) UPtr[0] = VlowRamp;
    else if (phaseRamp < (RampTD + RampTR)) UPtr[0] = VlowRamp + (SlopeRiseRamp*(phaseRamp - RampTD));
    else if (phaseRamp < (RampTD + RampTR + RampPW)) UPtr[0] = VhighRamp;
    else if (phaseRamp < (RampTD + RampTR + RampPW + RampTF)) UPtr[0] = VhighRamp - (SlopeFallRamp*(phaseRamp - (RampTD + RampTR + RampPW)));
    else UPtr[0] = VlowRamp;

    if (time < VrefSettlingTime) UPtr[1] = Vref*time/VrefSettlingTime;
    else UPtr[1] = Vref;

    z[sizeOfZ-2] = UPtr[0];
    z[sizeOfZ-1] = UPtr[1];

    b[1] = z[10];
    b[4] = z[9]*UPtr[1];

  }
