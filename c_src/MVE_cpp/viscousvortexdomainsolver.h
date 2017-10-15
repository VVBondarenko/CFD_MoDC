#ifndef VISCOUSVORTEXDOMAINSOLVER_H
#define VISCOUSVORTEXDOMAINSOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <vector>
#include <algorithm>
//#include <parallel/algorithm>
#include "pstl/execution"
#include "pstl/algorithm"
#include "tbb/parallel_for_each.h"
#include <mutex>

#include <omp.h>

#include <streamline.h>

#define VORTEX_IN_FLOW 1
#define VORTEX_OUT_OF_RANGE 2
#define VORTEX_IN_BODY 3


struct Vortex
{
    double x, y;
    double vorticity;

    int status;
    int generationSide;

    bool operator==(const Vortex& a) const
    {
        return (x==a.x) && (y==a.y) && (vorticity==a.vorticity);
    }
    bool operator!=(const Vortex& a) const
    {
        return !(*this == a);
    }
};

class VortexInMotion
{
    double ux, uy, tau;
public:
    Vortex Object;
    VortexInMotion(Vortex Object, double ux, double uy, double tau)
    {
        this->Object = Object;
        this->ux = ux;
        this->uy = uy;
        this->tau = tau;
    }

    VortexInMotion(Vortex Object, double tau)
    {
        this->Object = Object;
        this->ux = 0.;
        this->uy = 0.;
        this->tau = tau;
    }

    Vortex UpdatedVortex()
    {
        Object.x += ux*tau;
        Object.y += uy*tau;
        return Object;
    }
    void addVelocityByX(double vx)
    {
        ux += vx;
    }
    void addVelocityByY(double vy)
    {
        uy += vy;
    }
};



//ToDo: add method for resumption of calculations

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

    double maxEps;

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

    std::vector <Vortex> InFlow;
    std::vector <VortexInMotion> InFlowEvolution;
    std::mutex SyncFlowEvolution;


    std::vector <Vortex> InBodyVortexes;

    std::vector <StreamLinePoint> StreamLine1;
    std::vector <StreamLinePoint> StreamLine2;

    void UpdateStreamLines();

    double Qfield_x(double x, double y);
    double Qfield_y(double x, double y);

    double Pfield_x(double x, double y);
    double Pfield_y(double x, double y);

    double Mfield(double x, double y);

    int isVortexInBody(double x, double y); //not used (why?)
public:

    ViscousVortexDomainSolver();
    ~ViscousVortexDomainSolver();

    void Solve();
    void DivideProfileToPanels();
    void CompletingGeneratingMatrix();

    void Output_ParaView_Field(const char *FileName);
    void Output_ParaView_Field_with_exclusions(const char *FileName);
    void Output_ParaView_AllVortexes(const char *FileName);
    void Output_StreamLines(const char *FileName);

private:
    void UpdateVotexPositions();
    void ComputeVortexMotion(Vortex Considered);
    void ComputeAssociatedVortexes();
    double UpdateEpsilon(Vortex Considered, double InitialEpsilon);
    int LeakageControl(Vortex Checking);
    void AddNewVortexesToFlow();
    double xVelocityAt(double x, double y);
    double yVelocityAt(double x, double y);
};

#endif // VISCOUSVORTEXDOMAINSOLVER_H
