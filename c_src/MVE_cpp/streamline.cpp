#include "streamline.h"


StreamLinePoint::StreamLinePoint(double x, double y)
{
    this->x = x;
    this->y = y;
}

void StreamLinePoint::move(double dx, double dy)
{
    this->x += dx;
    this->y += dy;
}
