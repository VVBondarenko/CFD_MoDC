#ifndef STREAMLINE_H
#define STREAMLINE_H

#include <vector>
//#include <viscousvortexdomainsolver.h>


class ControlPointOfFlow
{
public:
    ControlPointOfFlow(double x, double y);

    double x,y;
    void move(double dx, double dy);
};
/*
class StreamLine
{
    ViscousVortexDomainSolver *Problem;
    double tau;

    std::vector <ControlPointOfFlow> LineElements;
    std::vector <ControlPointOfFlow> DeployPoints;

    void UpdateSingleElement(ControlPointOfFlow &UpdatingPoint);

public:
    StreamLine(ViscousVortexDomainSolver &Problem, double tau);
    void StepInTime();
    void UpdateElements();

    void addDeploys(double x, double y);
};
*/
#endif // STREAMLINE_H
