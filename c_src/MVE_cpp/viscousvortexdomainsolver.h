#ifndef VISCOUSVORTEXDOMAINSOLVER_H
#define VISCOUSVORTEXDOMAINSOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_linalg.h>

#include <vector>
#include <thread>
#include <mutex>

#define VORTEX_IN_FLOW 1
#define VORTEX_OUT_OF_RANGE 2
#define VORTEX_IN_BODY 3


typedef struct Vortex
{
    double x, y;
    double vorticity;

    int status;
    int generationSide;
} Vortex;

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

    std::vector<Vortex> InFlow;
    std::vector<Vortex> NextInFlow;
    std::vector<Vortex> InBodyVortexes;

    std::mutex NextInFlow_WriteMutex;

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

    void Output_ParaView_Field(const char *FileName);
    void Output_ParaView_Line(const char *FileName);

    void Thread_UpdateVortexPosition(int ID, int ThreadNum);
    static void ThreadCrutch_UpdateVotexPositions(ViscousVortexDomainSolver *Task, int ID, int ThreadNum);

private:
    void UpdateVotexPositions();
    void UpdateAssociatedVortexes();
    double UpdateEpsilon(int i, double InitialEpsilon);
    int LeakageControl(double u_x, int i, double NextY, double u_y, double NextX);
    void AddNewVortexesToFlow();
};

#endif // VISCOUSVORTEXDOMAINSOLVER_H
