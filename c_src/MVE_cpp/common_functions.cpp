#include "common_functions.h"

double Qfield_x(double x, double y)
{
    double r = x*x+y*y;
    return -y/(2.*M_PI*(r+0.0001*0.0001));
}

double Qfield_y(double x, double y)
{
    double r = x*x+y*y;
    return  x/(2.*M_PI*(r+0.0001*0.0001));
}

double Pfield_x(double x, double y)
{
    double r = hypot(x,y);
    if(r == 0)
        return 1.;
    return x/r*exp(-r);
}

double Pfield_y(double x, double y)
{
    double r = hypot(x,y);
    if(r == 0)
        return 1.;
    return y/r*exp(-r);
}

double Mfield(double x, double y)
{
    double r2= x*x+y*y;
    double r = sqrt(r2);
    if(x==0 && y==0)
        return 0.;
    return (r+1.)*exp(-r)/r2;
}
