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
    int ActiveVortexesInFLow;

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

    gsl_matrix *OriginalVortexGenerationMatrix;
    gsl_matrix *VortexGenerationMatrix;

    gsl_vector *VelocityProjections;
    gsl_vector *NewVorticities;

    gsl_permutation *Permutation;

    //ToDo: replace static array with vector container
    Vortex *InFlow;
    Vortex *NextInFlow;

    double Qfield_x(double x, double y);
    double Qfield_y(double x, double y);

    double Pfield_x(double x, double y);
    double Pfield_y(double x, double y);

    double Mfield(double x, double y);

    int VortexInBody(double x, double y); //not used (why?)
public:

    ViscousVortexDomainSolver();
    ~ViscousVortexDomainSolver();

    void Solve();
    void DivideProfileToPanels();
    void CompletingGeneratingMatrix();
    void UpdateVotexPositions();
    void UpdateAssociatedVortexes();

    void Output_ParaView_Field(const char *FileName);
    void Output_ParaView_Line(const char *FileName);
};

#endif // VISCOUSVORTEXDOMAINSOLVER_H
