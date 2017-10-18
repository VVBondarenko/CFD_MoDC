#include <viscousvortexdomainsolver.h>

int main()
{

    ViscousVortexDomainSolver *TestSolver = new ViscousVortexDomainSolver();
    TestSolver->Solve();
    delete TestSolver;

    return 0;
}
