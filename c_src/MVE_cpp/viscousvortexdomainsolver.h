#ifndef VISCOUSVORTEXDOMAINSOLVER_H
#define VISCOUSVORTEXDOMAINSOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_linalg.h>

typedef struct Vortex
{
    double x, y;
    double vorticity;

    int active;
} Vortex;


class ViscousVortexDomainSolver
{
    int NumberOfPanels;
    int MaxVortexes;
    int Niterations;
    double nu;
    double tau;
    double rho;
    double vx_inf;
    double vy_inf;

    double Qfield_x(double x, double y);
    double Qfield_y(double x, double y);

    double Pfield_x(double x, double y);
    double Pfield_y(double x, double y);

    double Mfield(double x, double y);

    int VortexInBody(double x, double y);
public:
    ViscousVortexDomainSolver();

    void Solve();


};

#endif // VISCOUSVORTEXDOMAINSOLVER_H
