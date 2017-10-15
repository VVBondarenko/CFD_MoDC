#include <viscousvortexdomainsolver.h>

int main()
{
//    omp_set_dynamic(false);
//    omp_set_num_threads(2);

    ViscousVortexDomainSolver *TestSolver = new ViscousVortexDomainSolver();
    TestSolver->Solve();
    delete TestSolver;

    return 0;
}
