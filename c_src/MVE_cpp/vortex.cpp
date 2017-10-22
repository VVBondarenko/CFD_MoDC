#include "vortex.h"

Vortex::Vortex(double X, double Y, double Vorticity)
{
    this->x = X;
    this->y = Y;
    this->vorticity = Vorticity;
}

bool Vortex::operator==(const Vortex &a) const
{
    return (x==a.x) && (y==a.y) && (vorticity==a.vorticity);
}

bool Vortex::operator!=(const Vortex &a) const
{
    return !(*this == a);
}

double Vortex::inducedXVelocity(double x, double y)
{
    return vorticity*Qfield_x(x-this->x, y-this->y);
}

double Vortex::inducedYVelocity(double x, double y)
{
    return vorticity*Qfield_y(x-this->x, y-this->y);
}
