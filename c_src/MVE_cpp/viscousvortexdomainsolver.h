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
#include <vortex.h>
#include <vortexmotion.h>

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
