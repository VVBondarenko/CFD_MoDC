#include "viscousvortexdomainsolver.h"

ViscousVortexDomainSolver::ViscousVortexDomainSolver()
{
    NumberOfPanels = 200;
    MaxVortexes = 16500;
    Niterations = 2000;

    nu = 1./53.;
    tau = 0.01;
    rho = 1.;
    vx_inf = 3.;
    vy_inf = 0.;

    PanelLength = new double  [NumberOfPanels];
    PanelMids   = new double* [NumberOfPanels];
    PanelNodes  = new double* [NumberOfPanels];
    PanelNorms  = new double* [NumberOfPanels];
    PanelTaus   = new double* [NumberOfPanels];
    DeployPoints= new double* [NumberOfPanels];
    for(int i=0;i<NumberOfPanels;i++)
    {
        PanelMids   [i] = new double [2];
        PanelNodes  [i] = new double [2];
        PanelNorms  [i] = new double [2];
        PanelTaus   [i] = new double [2];
        DeployPoints[i] = new double [2];
    }

    OriginalVortexGenerationMatrix = gsl_matrix_alloc(NumberOfPanels+1,NumberOfPanels+1);
    VortexGenerationMatrix         = gsl_matrix_alloc(NumberOfPanels+1,NumberOfPanels+1);
    VortexGenerationMatrix_Inversed= gsl_matrix_alloc(NumberOfPanels+1,NumberOfPanels+1);

    VelocityProjections = gsl_vector_alloc(NumberOfPanels+1);
    NewVorticities      = gsl_vector_alloc(NumberOfPanels+1);

    Permutation = gsl_permutation_alloc (NumberOfPanels+1);

}

ViscousVortexDomainSolver::~ViscousVortexDomainSolver()
{
    delete [] PanelLength;
    for(int i = 0; i < NumberOfPanels; i++)
    {
        delete [] PanelMids[i];
        delete [] PanelNodes[i];
        delete [] PanelNorms[i];
        delete [] PanelTaus[i];
        delete [] DeployPoints[i];
    }
    delete [] PanelMids;
    delete [] PanelNodes;
    delete [] PanelNorms;
    delete [] PanelTaus;
    delete [] DeployPoints;

    gsl_matrix_free(OriginalVortexGenerationMatrix);
    gsl_matrix_free(VortexGenerationMatrix);
    gsl_matrix_free(VortexGenerationMatrix_Inversed);

    gsl_vector_free(VelocityProjections);
    gsl_vector_free(NewVorticities);

    gsl_permutation_free(Permutation);

}

void ViscousVortexDomainSolver::Solve()
{
    DivideProfileToPanels();
    CompletingGeneratingMatrix();

    ComputeAssociatedVortexes();

    FILE *forces;
    forces = fopen("forces.dat","w");

    for(int iterator=0; iterator<Niterations; iterator++)
    {
        double f_x = 0., f_y = 0./*, M_z = 0.*/;

        AppendAssociatedVortexesToFlow();

        printf("active vortexes \t%d,\t iteration\t%d\n",(int)Flow.size(),iterator);

        StepInTimeForStreamLines();
        UpdateVotexPositions();

        for(int m = 0; m < NumberOfPanels; m++)
        {
            f_x += -rho*(gsl_vector_get(NewVorticities,m)/tau* PanelMids[m][1]);
            f_y +=  rho*(gsl_vector_get(NewVorticities,m)/tau* PanelMids[m][0]);
        }
        fprintf(forces,"%f %f %f\n",iterator*tau, f_x, f_y);

        ComputeAssociatedVortexes();

        if(iterator%3==0 || 1)
        {
            char FileName[32];
            sprintf(FileName,"ControlPoints%6.6d.csv",iterator);
            Output_StreamLines(FileName);
        }
    }
}

void ViscousVortexDomainSolver::DivideProfileToPanels()
{
    int i;
    for(i=0;i<NumberOfPanels;i++)
    {
        double PhiStep = 2.*M_PI/((double)NumberOfPanels);
        PanelNodes[i][0] = cos(PhiStep*i);
        PanelNodes[i][1] = sin(PhiStep*i);

        DeployPoints[i][0] = (1.+PhiStep*1.e-8)*PanelNodes[i][0];
        DeployPoints[i][1] = (1.+PhiStep*1.e-8)*PanelNodes[i][1];

        PanelNorms[i][0] = cos(PhiStep*(i+0.5));
        PanelNorms[i][1] = sin(PhiStep*(i+0.5));

    }

    for(i=0;i<NumberOfPanels-1;i++)
    {
        PanelMids[i][0] = (PanelNodes[i][0]+PanelNodes[i+1][0])/2;
        PanelMids[i][1] = (PanelNodes[i][1]+PanelNodes[i+1][1])/2;

        PanelTaus[i][0] = (PanelNodes[i+1][0]-PanelNodes[i][0]);
        PanelTaus[i][1] = (PanelNodes[i+1][1]-PanelNodes[i][1]);
        PanelLength[i]  = sqrt(PanelTaus[i][0]*PanelTaus[i][0]+
                                PanelTaus[i][1]*PanelTaus[i][1]);
        PanelTaus[i][0] /= PanelLength[i];
        PanelTaus[i][1] /= PanelLength[i];
        //            PanelNorms[i][0] = -PanelTaus[i][1];
        //            PanelNorms[i][1] = PanelTaus[i][0];
    }
    PanelMids[NumberOfPanels-1][0] = (PanelNodes[NumberOfPanels-1][0]+PanelNodes[0][0])/2;
    PanelMids[NumberOfPanels-1][1] = (PanelNodes[NumberOfPanels-1][1]+PanelNodes[0][1])/2;
    PanelTaus[NumberOfPanels-1][0] = (PanelNodes[0][0]-PanelNodes[NumberOfPanels-1][0]);
    PanelTaus[NumberOfPanels-1][1] = (PanelNodes[0][1]-PanelNodes[NumberOfPanels-1][1]);
    PanelLength[NumberOfPanels-1]  = hypot(PanelTaus[NumberOfPanels-1][0],PanelTaus[NumberOfPanels-1][1]);
    PanelTaus[NumberOfPanels-1][0] /= PanelLength[NumberOfPanels-1];
    PanelTaus[NumberOfPanels-1][1] /= PanelLength[NumberOfPanels-1];
    //        PanelNorms[Np-1][0] = -PanelTaus[Np-1][1];
    //        PanelNorms[Np-1][1] = PanelTaus[Np-1][0];

}

void ViscousVortexDomainSolver::CompletingGeneratingMatrix()
{
    int i,j,k;
    for(i=0; i<NumberOfPanels; i++)
    {
        for(j=0; j<NumberOfPanels; j++)
        {
            gsl_matrix_set(VortexGenerationMatrix,i,j,
            PanelNorms[i][0]*Qfield_x(PanelMids[i][0]-DeployPoints[j][0],
                                      PanelMids[i][1]-DeployPoints[j][1])
           +PanelNorms[i][1]*Qfield_y(PanelMids[i][0]-DeployPoints[j][0],
                                      PanelMids[i][1]-DeployPoints[j][1]));
        }
    }

    for(i=0;i<NumberOfPanels;i++)
    {
        gsl_matrix_set(VortexGenerationMatrix,NumberOfPanels,i,1);
        gsl_matrix_set(VortexGenerationMatrix,i,NumberOfPanels,1);
    }
    gsl_matrix_set(VortexGenerationMatrix,NumberOfPanels,NumberOfPanels,0);
    gsl_matrix_memcpy(OriginalVortexGenerationMatrix,VortexGenerationMatrix);

    gsl_linalg_LU_decomp(VortexGenerationMatrix,Permutation,&k);
    gsl_linalg_LU_invert(VortexGenerationMatrix,Permutation,VortexGenerationMatrix_Inversed);
    gsl_matrix_memcpy(VortexGenerationMatrix, OriginalVortexGenerationMatrix);

}

void ViscousVortexDomainSolver::AppendAssociatedVortexesToFlow()
{
    for(int i=0;i<NumberOfPanels;i++)
    {
        Vortex NewVortexFromPanel = Vortex(DeployPoints[i][0],DeployPoints[i][1],gsl_vector_get(NewVorticities,i));

        if((2*i)/(NumberOfPanels))
            NewVortexFromPanel.deployIdentifier = 0.;
        else
            NewVortexFromPanel.deployIdentifier = 1.;

        Flow.push_back(NewVortexFromPanel);
    }
}


void ViscousVortexDomainSolver::ComputeAssociatedVortexes()
{
    int j,k;

    double VelocityProjInControlPoint[NumberOfPanels];
    for(j=0; j<NumberOfPanels; j++)
        VelocityProjInControlPoint[j] = PanelNorms[j][0]*vx_inf+
                                        PanelNorms[j][1]*vy_inf;

    //ToDo: think about ways of parallelizing of this part of code (OpenMP, possibly)
    //This part is obviously parallelizable over j variable
#pragma omp parallel for private(j)
    for(j=0; j<NumberOfPanels; j++)
    {

        for(Vortex Considered : Flow)
        {
            VelocityProjInControlPoint[j] +=
                    (PanelNorms[j][0]*Considered.inducedXVelocity(PanelMids[j][0],PanelMids[j][1])
                    +PanelNorms[j][1]*Considered.inducedYVelocity(PanelMids[j][0],PanelMids[j][1]));
        }
        gsl_vector_set(VelocityProjections,j,-VelocityProjInControlPoint[j]);
    }
    gsl_vector_set(VelocityProjections,NumberOfPanels,0);

    gsl_blas_dgemv (CblasNoTrans, 1., VortexGenerationMatrix_Inversed, VelocityProjections, 0., NewVorticities);

}



void ViscousVortexDomainSolver::UpdateVotexPositions()
{
    FlowEvolution.clear();
    FlowEvolution.reserve(Flow.size());

    tbb::parallel_for_each
            (Flow.begin(),Flow.end(),[this](Vortex& Considered)
    {
        ComputeVortexMotion(Considered);
    });


    tbb::parallel_for_each
            (FlowEvolution.begin(), FlowEvolution.end(),
             [this](VortexMotion i)
    {
        i.UpdatedVortex();
    }
    );

    //ToDo: modify to consider in body vortexes
    Flow.erase(
                std::remove_if(Flow.begin(),
                               Flow.end(),
                               [this](Vortex Considered)
                {
                     return LeakageControl(Considered)!=VORTEX_IN_FLOW;
                }),
                Flow.end());
}

void ViscousVortexDomainSolver::ComputeVortexMotion(Vortex& Considered)
{
    double eps = sqrt(tau*nu/2.);
    double  I0 = 0.,
            I1 = 0.,
            I2_x = 0., I2_y = 0.,
            I3_x = 0., I3_y = 0.;

    double u_x = vx_inf;
    double u_y = vy_inf;


    eps = fmax(UpdateEpsilon(Considered, eps),eps);
    I0 = 2.*M_PI*eps*eps;


    //Vortex-vortex interaction
    //  (already implemented somewhere else... possibly, should be replaced. Но это не точно)
    for(Vortex Processing : Flow)
    {
        if( Considered!=Processing )
        {
            u_x += Processing.inducedXVelocity(Considered.x,Considered.y);
            u_y += Processing.inducedYVelocity(Considered.x,Considered.y);
        }
    }

    for(Vortex Processing : Flow)
    {
        double rx = (Considered.x-Processing.x)/eps;
        double ry = (Considered.y-Processing.y)/eps;
        if( Processing != Considered)
        {
            I2_x -= Processing.vorticity*Pfield_x(rx, ry)/eps;
            I2_y -= Processing.vorticity*Pfield_y(rx, ry)/eps;

            I1   += Processing.vorticity*exp(-hypot(rx, ry));
        }
    }

    //viscous surface-vortex interaction speed
    for(int j=0;j<NumberOfPanels;j++)
    {
        double Rr = hypot((Considered.x-PanelMids[j][0]),(Considered.y-PanelMids[j][1]));

        if(Rr>3.*PanelLength[j] && Rr<20.*PanelLength[j])
        {
            I3_x -= PanelNorms[j][0]*exp(-Rr/eps)*PanelLength[j];
            I3_y -= PanelNorms[j][1]*exp(-Rr/eps)*PanelLength[j];

            I0   -= ((Considered.x-PanelMids[j][0])*PanelNorms[j][0]
                    +(Considered.y-PanelMids[j][1])*PanelNorms[j][1])
                    *PanelLength[j]
                    *Mfield((Considered.x-PanelMids[j][0])/eps,
                            (Considered.y-PanelMids[j][1])/eps);
        }
        else if(Rr>0.1*PanelLength[j])
        {
            int Ndiv = 4;
            for(int k = 0; k < Ndiv; k++)
            {
                double positionCoeff = (1.+2.*k)/(2*Ndiv);
                double xOffset = PanelTaus[j][0]*PanelLength[j]*positionCoeff;
                double yOffset = PanelTaus[j][1]*PanelLength[j]*positionCoeff;

                Rr = hypot((Considered.x-xOffset),(Considered.y-yOffset));

                I3_x -= PanelNorms[j][0]*exp(-Rr/eps)*PanelLength[j]*(1./Ndiv);
                I3_y -= PanelNorms[j][1]*exp(-Rr/eps)*PanelLength[j]*(1./Ndiv);

                I0   -= ((Considered.x-xOffset)*PanelNorms[j][0]
                        +(Considered.y-yOffset)*PanelNorms[j][1])
                        *PanelLength[j]*(1./Ndiv)
                        *Mfield((Considered.x-xOffset)/eps,
                                (Considered.y-yOffset)/eps);
            }
        }
        else
        {
            I3_x = 2.*PanelNorms[j][0]*eps*(1.-exp(-PanelLength[j]*0.5/eps));
            I3_y = 2.*PanelNorms[j][1]*eps*(1.-exp(-PanelLength[j]*0.5/eps));
            I0   = M_PI*eps*eps;
            break;
        }
    }
    u_x += nu*(-I2_x/I1+I3_x/I0);
    u_y += nu*(-I2_y/I1+I3_y/I0);

    SyncFlowEvolution.lock();
    FlowEvolution.push_back(VortexMotion(Considered,u_x,u_y,tau));
    SyncFlowEvolution.unlock();
}


double ViscousVortexDomainSolver::UpdateEpsilon(Vortex Considered, double InitialEpsilon)
{
    if((int)Flow.size()>NumberOfPanels)
    {
        double sq1=InitialEpsilon*InitialEpsilon;
        double sq2=sq1, sq3=sq1;
        for(Vortex Processing : Flow)
        {
            if(Considered!=Processing)
            {
                double temp = (Considered.x-Processing.x)*(Considered.x-Processing.x)
                             +(Considered.y-Processing.y)*(Considered.y-Processing.y);
                if(
                        temp < sq1 //&&
//                        Considered.vorticity*Processing.vorticity>0
                        )
                {
                    sq1 = temp;
                    sq2 = sq1;
                    sq3 = sq2;
                }
            }
        }
        InitialEpsilon = sqrt( (sq1+sq2+sq3)/3 );
    }
    return InitialEpsilon;
}


int ViscousVortexDomainSolver::LeakageControl(Vortex Checking)
{
    if(Checking.x < 1000 && Checking.y < 1000 )
    {
        if(isVortexInBody(Checking.x, Checking.y))
            return VORTEX_IN_FLOW;
        else
            return VORTEX_IN_BODY;
    }
    else
    {
        return VORTEX_OUT_OF_RANGE;
    }
}


int ViscousVortexDomainSolver::isVortexInBody(double x, double y)
{
    if(x*x+y*y<1)
        return 0;
    else
        return 1;
}


double ViscousVortexDomainSolver::xVelocityAt(double x, double y)
{
    double vx=vx_inf;
    for(Vortex Considered : Flow)
        vx += Considered.inducedXVelocity(x,y);

    return vx;
}

double ViscousVortexDomainSolver::yVelocityAt(double x, double y)
{
    double vy=vy_inf;
    for(Vortex Considered : Flow)
        vy += Considered.inducedYVelocity(x,y);

    return vy;
}

void ViscousVortexDomainSolver::Output_ParaView_Field_with_exclusions(const char *FileName)
{
    int i,j;
    FILE *VelocityField;
    VelocityField = fopen(FileName,"w");
    fprintf(VelocityField,"TITLE = \"Flow Model\"\n");
    fprintf(VelocityField,"VARIABLES = \"x\", \"y\", \"vx\", \"vy\"\n");
    fprintf(VelocityField,"ZONE T=\"Frame 0\", I=%d, J=%d\n", 200, 100);
    for(i=0;i<200;i++)
    {
        for(j=0;j<100;j++)
        {
            double x = 0.08*i-4,
                   y = 0.08*j-4;
            double vx = xVelocityAt(x, y),
                   vy = yVelocityAt(x, y);

            fprintf(VelocityField,"%f %f %f %f\n",x, y, vx, vy);
        }
    }
    fclose(VelocityField);
}

void ViscousVortexDomainSolver::Output_ParaView_Field(const char *FileName)
{
    int i,j,n;
    FILE *VelocityField;
    VelocityField = fopen(FileName,"w");
    fprintf(VelocityField,"TITLE = \"Flow Model\"\n");
    fprintf(VelocityField,"VARIABLES = \"x\", \"y\", \"vx\", \"vy\"\n");
    fprintf(VelocityField,"ZONE T=\"Frame 0\", I=%d, J=%d\n", 200, 100);
    for(i=0;i<200;i++)
    {
        for(j=0;j<100;j++)
        {
            double vx=vx_inf, vy=vy_inf;
            for(n=0;n<(int)Flow.size();n++)
            {
                    vx += Flow[n].vorticity*Qfield_x((0.08*i-4)-Flow[n].x,
                                                       (0.08*j-4)-Flow[n].y);
                    vy += Flow[n].vorticity*Qfield_y((0.08*i-4)-Flow[n].x,
                                                       (0.08*j-4)-Flow[n].y);
            }
            fprintf(VelocityField,"%f %f %f %f\n",0.08*i-4, 0.08*j-4, vx, vy);
        }
    }
    fclose(VelocityField);
}

void ViscousVortexDomainSolver::Output_ParaView_AllVortexes(const char *FileName)
{
    int n;
    FILE *VelocityField;
    VelocityField = fopen(FileName,"w");
    fprintf(VelocityField,"X,Y,Z,vorticity,side\n");

    double average = 0.;
    for(n=0;n<(int)Flow.size();n++)
        average += Flow[n].vorticity*Flow[n].vorticity;
    average = sqrt(average)/((double)((int)Flow.size()));

    for(n=0;n<(int)Flow.size();n++)
    {
        if(fabs(fabs(Flow[n].vorticity)-average)/average > 0.9 )
            fprintf(VelocityField,"%f,%f,%f,%f,%f\n",Flow[n].x,Flow[n].y,0.,
                                                     Flow[n].vorticity, (double)Flow[n].deployIdentifier);

    }
    fclose(VelocityField);
}



void ViscousVortexDomainSolver::UpdateControlPointsOfStreamLine(std::vector <ControlPointOfFlow> &StreamLine)
{
    tbb::parallel_for_each
            (StreamLine.begin(),StreamLine.end(),[this](ControlPointOfFlow &UpdatingPoint)
    {
        MoveSingleControlPoint(UpdatingPoint);
    });
}

void ViscousVortexDomainSolver::MoveSingleControlPoint(ControlPointOfFlow &UpdatingPoint)
{
    ControlPointOfFlow PredictionPoint = UpdatingPoint;
    PredictionPoint.move(tau*xVelocityAt(UpdatingPoint.x,UpdatingPoint.y),
                         tau*yVelocityAt(UpdatingPoint.x,UpdatingPoint.y));
    UpdatingPoint.move(0.5*tau*(xVelocityAt(UpdatingPoint.x,UpdatingPoint.y)+
                                xVelocityAt(PredictionPoint.x,PredictionPoint.y)),
                       0.5*tau*(yVelocityAt(UpdatingPoint.x,UpdatingPoint.y)+
                                yVelocityAt(PredictionPoint.x,PredictionPoint.y)));
}

void ViscousVortexDomainSolver::StepInTimeForStreamLines()
{
    //ToDo: move to exernal program
    ControlPointOfFlow NewPoint1 = ControlPointOfFlow(0.,-1.1);
    ControlPointOfFlow NewPoint2 = ControlPointOfFlow(0., 1.1);


    UpdateControlPointsOfStreamLine(StreamLine1);
    UpdateControlPointsOfStreamLine(StreamLine2);


    StreamLine1.push_back(NewPoint1);
    StreamLine2.push_back(NewPoint2);
}

void ViscousVortexDomainSolver::Output_StreamLines(const char *FileName)
{
    FILE *ControlPoints;
    ControlPoints = fopen(FileName,"w");
    fprintf(ControlPoints,"X,Y,Z,side\n");

    for(ControlPointOfFlow i : StreamLine1)
        fprintf(ControlPoints, "%f,%f,%f,%f\n", i.x, i.y, 0., 0.);

    for(ControlPointOfFlow i : StreamLine2)
        fprintf(ControlPoints, "%f,%f,%f,%f\n", i.x, i.y, 0., 1.);
    fclose(ControlPoints);
}
