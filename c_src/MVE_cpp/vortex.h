#ifndef VORTEX_H
#define VORTEX_H

#include <common_functions.h>

#define VORTEX_IN_FLOW 1
#define VORTEX_OUT_OF_RANGE 2
#define VORTEX_IN_BODY 3


class Vortex
{
public:
    Vortex(double X, double Y, double Vorticity);
    double x, y;
    double vorticity;

    double deployIdentifier;

    bool operator==(const Vortex& a) const;
    bool operator!=(const Vortex& a) const;

    double inducedXVelocity(double x, double y);
    double inducedYVelocity(double x, double y);
};


#endif // VORTEX_H
