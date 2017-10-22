#ifndef VORTEXMOTION_H
#define VORTEXMOTION_H

#include <vortex.h>

class VortexMotion
{
    double ux, uy, tau;
public:
    Vortex* MovingVortex;

    VortexMotion(Vortex& Object, double ux, double uy, double tau);
    VortexMotion(Vortex& Object, double tau);

    void UpdatedVortex();
    void addVelocityByX(double vx);
    void addVelocityByY(double vy);
};

#endif // VORTEXMOTION_H
