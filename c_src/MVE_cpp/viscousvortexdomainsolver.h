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

    double  *PanelLength;
    double  **PanelNodes,
            **PanelMids,
            **PanelNorms,
            **PanelTaus,
            **DeployPoints;

    gsl_matrix *OriginalVortexGenerationMatrix;// = gsl_matrix_alloc(NumberOfPanels+1,NumberOfPanels+1);
    gsl_matrix *VortexGenerationMatrix;// = gsl_matrix_alloc(NumberOfPanels+1,NumberOfPanels+1);

    double Qfield_x(double x, double y);
    double Qfield_y(double x, double y);

    double Pfield_x(double x, double y);
    double Pfield_y(double x, double y);

    double Mfield(double x, double y);

    int VortexInBody(double x, double y);
public:
    ViscousVortexDomainSolver();
    ~ViscousVortexDomainSolver();

    void Solve();
    void DivideProfileToPanels();
    void CompletingGeneratingMatrix();
};

#endif // VISCOUSVORTEXDOMAINSOLVER_H
