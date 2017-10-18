#include "streamline.h"


ControlPointOfFlow::ControlPointOfFlow(double x, double y)
{
    this->x = x;
    this->y = y;
}

void ControlPointOfFlow::move(double dx, double dy)
{
    this->x += dx;
    this->y += dy;
}
