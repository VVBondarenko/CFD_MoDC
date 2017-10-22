#include "vortexmotion.h"

VortexMotion::VortexMotion(Vortex &Object, double ux, double uy, double tau)
{
    this->MovingVortex = &Object;
    this->ux = ux;
    this->uy = uy;
    this->tau = tau;
}

VortexMotion::VortexMotion(Vortex &Object, double tau)
{
    this->MovingVortex = &Object;
    this->ux = 0.;
    this->uy = 0.;
    this->tau = tau;
}

void VortexMotion::UpdatedVortex()
{
    (*MovingVortex).x += ux*tau;
    (*MovingVortex).y += uy*tau;
}

void VortexMotion::addVelocityByX(double vx)
{
    ux += vx;
}

void VortexMotion::addVelocityByY(double vy)
{
    uy += vy;
}
