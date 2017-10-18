#ifndef VISCOUSVORTEXDOMAINSOLVER_H
#define VISCOUSVORTEXDOMAINSOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <common_functions.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
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


class Vortex
{
public:
    Vortex(double X, double Y, double Vorticity)
    {
        this->x = X;
        this->y = Y;
        this->vorticity = Vorticity;
    }
    double x, y;
    double vorticity;

    double deployIdentifier;

    bool operator==(const Vortex& a) const
    {
        return (x==a.x) && (y==a.y) && (vorticity==a.vorticity);
    }
    bool operator!=(const Vortex& a) const
    {
        return !(*this == a);
    }

    double inducedXVelocity(double x, double y)
    {
        return vorticity*Qfield_x(x-this->x, y-this->y);
    }

    double inducedYVelocity(double x, double y)
    {
        return vorticity*Qfield_y(x-this->x, y-this->y);
    }

};

class VortexMotion
{
    double ux, uy, tau;
public:
    Vortex* MovingVortex;
    VortexMotion(Vortex& Object, double ux, double uy, double tau)
    {
        this->MovingVortex = &Object;
        this->ux = ux;
        this->uy = uy;
        this->tau = tau;
    }

    VortexMotion(Vortex& Object, double tau)
    {
        this->MovingVortex = &Object;
        this->ux = 0.;
        this->uy = 0.;
        this->tau = tau;
    }

    void UpdatedVortex()
    {
        (*MovingVortex).x += ux*tau;
        (*MovingVortex).y += uy*tau;
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



//ToDo: add resumption of calculations mechanism
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

    gsl_matrix *OriginalVortexGenerationMatrix;
    gsl_matrix *VortexGenerationMatrix;
    gsl_matrix *VortexGenerationMatrix_Inversed;
    gsl_vector *VelocityProjections;
    gsl_vector *NewVorticities;
    gsl_permutation *Permutation;


    std::vector <Vortex> Flow;
    std::vector <VortexMotion> FlowEvolution;
    std::mutex SyncFlowEvolution;
    std::vector <Vortex> InBodyVortexes;


    std::vector <ControlPointOfFlow> StreamLine1;
    std::vector <ControlPointOfFlow> StreamLine2;
    void StepInTimeForStreamLines();
    void UpdateControlPointsOfStreamLine(std::vector<ControlPointOfFlow> &StreamLine);
    void MoveSingleControlPoint(ControlPointOfFlow &UpdatingPoint);

    int isVortexInBody(double x, double y); //ToDo: update to be virtual

    void DivideProfileToPanels();
    void CompletingGeneratingMatrix();
    void ComputeAssociatedVortexes();
    void AppendAssociatedVortexesToFlow();
    void UpdateVotexPositions();
    void ComputeVortexMotion(Vortex &Considered);
    double UpdateEpsilon(Vortex Considered, double InitialEpsilon);
    int LeakageControl(Vortex Checking);

public:

    ViscousVortexDomainSolver();
    ~ViscousVortexDomainSolver();

    void Solve();

    double xVelocityAt(double x, double y);
    double yVelocityAt(double x, double y);

    void Output_ParaView_Field(const char *FileName);
    void Output_ParaView_Field_with_exclusions(const char *FileName);
    void Output_ParaView_AllVortexes(const char *FileName);
    void Output_StreamLines(const char *FileName);
};

#endif // VISCOUSVORTEXDOMAINSOLVER_H
