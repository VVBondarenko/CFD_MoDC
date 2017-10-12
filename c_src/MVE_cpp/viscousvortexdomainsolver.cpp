#include "viscousvortexdomainsolver.h"

ViscousVortexDomainSolver::ViscousVortexDomainSolver()
{
    maxEps = 0.;
    NumberOfPanels = 200;
    MaxVortexes = 16500;
    Niterations = 500;
    ActiveVortexesInFLow=0;

    nu = 1./53.;
    tau = 0.03;
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

    gsl_vector_free(VelocityProjections);
    gsl_vector_free(NewVorticities);

    gsl_permutation_free(Permutation);

}

void ViscousVortexDomainSolver::AddNewVortexesToFlow()
{
    for(int i=0;i<NumberOfPanels;i++)
    {
        Vortex NewVortexFromPanel;
        NewVortexFromPanel.x = DeployPoints[i][0];
        NewVortexFromPanel.y = DeployPoints[i][1];
        NewVortexFromPanel.status = VORTEX_IN_FLOW;
        NewVortexFromPanel.vorticity = gsl_vector_get(NewVorticities,i);
        if((2*i)/(NumberOfPanels))
            NewVortexFromPanel.generationSide = 0;
        else
            NewVortexFromPanel.generationSide = 1;



        InFlow.push_back(NewVortexFromPanel);
    }
    ActiveVortexesInFLow = InFlow.size();
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

        AddNewVortexesToFlow();

        printf("active vortexes \t%d,\t iteration\t%d\t%e\n",ActiveVortexesInFLow,iterator,maxEps);

        UpdateStreamLines();
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
//            sprintf(FileName,"VorticityField%6.6d.csv",iterator);
//            Output_ParaView_AllVortexes(FileName);
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
    int i,j;
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

}

double ViscousVortexDomainSolver::UpdateEpsilon(int i, double InitialEpsilon)
{
    int j;
    if((int)InFlow.size()>NumberOfPanels)
    {
        double sq1=InitialEpsilon*InitialEpsilon;
        double sq2=sq1, sq3=sq1;
        for(j=0;j<(int)InFlow.size();j++)
        {
            if(i!=j)
            {
                double temp = (InFlow[i].x-InFlow[j].x)*(InFlow[i].x-InFlow[j].x)
                             +(InFlow[i].y-InFlow[j].y)*(InFlow[i].y-InFlow[j].y);
                if(
                        temp < sq1 //&&
//                        InFlow[i].vorticity*InFlow[j].vorticity>0
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

double ViscousVortexDomainSolver::UpdateEpsilon(Vortex Considered, double InitialEpsilon)
{
    if((int)InFlow.size()>NumberOfPanels)
    {
        double sq1=InitialEpsilon*InitialEpsilon;
        double sq2=sq1, sq3=sq1;
        for(Vortex Processing : InFlow)
        {
            if(Considered!=Processing)
            {
                double temp = (Considered.x-Processing.x)*(Considered.x-Processing.x)
                             +(Considered.y-Processing.y)*(Considered.y-Processing.y);
                if(
                        temp < sq1 //&&
//                        InFlow[i].vorticity*InFlow[j].vorticity>0
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



void ViscousVortexDomainSolver::ComputeAssociatedVortexes()
{
    int i,j,k;

    double VelocityProjInControlPoint[NumberOfPanels];
    for(j=0; j<NumberOfPanels; j++)
        VelocityProjInControlPoint[j] = PanelNorms[j][0]*vx_inf+
                                        PanelNorms[j][1]*vy_inf;
    if(!InFlow.empty())
    {
        for(i=0; i<(int)InFlow.size(); i++)
        {
            for(j=0; j<NumberOfPanels; j++)
            {
                VelocityProjInControlPoint[j] +=
                        (PanelNorms[j][0]*Qfield_x((PanelMids[j][0]-InFlow[i].x),
                        (PanelMids[j][1]-InFlow[i].y))
                        +PanelNorms[j][1]*Qfield_y((PanelMids[j][0]-InFlow[i].x),
                        (PanelMids[j][1]-InFlow[i].y)))
                        *InFlow[i].vorticity;
            }
        }
    }

    for(j=0; j<NumberOfPanels; j++)
        gsl_vector_set(VelocityProjections,j,-VelocityProjInControlPoint[j]);
    gsl_vector_set(VelocityProjections,NumberOfPanels,0);

    gsl_matrix_memcpy(VortexGenerationMatrix, OriginalVortexGenerationMatrix);

    gsl_linalg_LU_decomp(VortexGenerationMatrix,Permutation,&k);
    gsl_linalg_LU_solve (VortexGenerationMatrix,Permutation,VelocityProjections,NewVorticities);

}

double ViscousVortexDomainSolver::xVelocityAt(double x, double y)
{
    double vx=vx_inf;
    for(int n = 0; n < (int)InFlow.size(); n++)
    {
        if(hypot(x-InFlow[n].x, y-InFlow[n].y) > sqrt(tau*nu/2.))
            vx += InFlow[n].vorticity*Qfield_x(x-InFlow[n].x,
                                               y-InFlow[n].y);
    }

    return vx;
}

double ViscousVortexDomainSolver::yVelocityAt(double x, double y)
{
    double vy=vy_inf;
    for(int n = 0; n < (int)InFlow.size(); n++)
    {
        if(hypot(x-InFlow[n].x, y-InFlow[n].y) > sqrt(tau*nu/2.))
            vy += InFlow[n].vorticity*Qfield_y(x-InFlow[n].x,
                                               y-InFlow[n].y);
    }
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
            for(n=0;n<(int)InFlow.size();n++)
            {
                    vx += InFlow[n].vorticity*Qfield_x((0.08*i-4)-InFlow[n].x,
                                                       (0.08*j-4)-InFlow[n].y);
                    vy += InFlow[n].vorticity*Qfield_y((0.08*i-4)-InFlow[n].x,
                                                       (0.08*j-4)-InFlow[n].y);
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
    for(n=0;n<(int)InFlow.size();n++)
        average += InFlow[n].vorticity*InFlow[n].vorticity;
    average = sqrt(average)/((double)((int)InFlow.size()));

    for(n=0;n<(int)InFlow.size();n++)
    {
        if(fabs(fabs(InFlow[n].vorticity)-average)/average > 0.9 )
            fprintf(VelocityField,"%f,%f,%f,%f,%f\n",InFlow[n].x,InFlow[n].y,0.,
                                                     InFlow[n].vorticity, (double)InFlow[n].generationSide);

    }
    fclose(VelocityField);
}


void ViscousVortexDomainSolver::UpdateVotexPositions()
{
    InFlowEvolution.clear();

//    for(Vortex Considered : InFlow)
//        ComputeVortexMotion(Considered);
    std::for_each(InFlow.begin(),InFlow.end(),[this](Vortex Considered)
    {
        ComputeVortexMotion(Considered);
    });

    InFlow.clear();

    for(VortexInMotion i : InFlowEvolution)
    {
        Vortex UpdatingVortex = i.UpdatedVortex();
        int status = LeakageControl(UpdatingVortex);
        if(status == VORTEX_IN_FLOW)
            InFlow.push_back(UpdatingVortex);

        if(status == VORTEX_IN_BODY)
            InBodyVortexes.push_back(UpdatingVortex);
    }
}

void ViscousVortexDomainSolver::ComputeVortexMotion(Vortex Considered)
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
    for(Vortex Processing : InFlow)
    {
        double rx = Considered.x-Processing.x;
        double ry = Considered.y-Processing.y;
        if( Considered!=Processing )
        {
            u_x += Processing.vorticity*Qfield_x(rx,ry);
            u_y += Processing.vorticity*Qfield_y(rx,ry);

            I2_x -= Processing.vorticity*Pfield_x((rx)/eps,
                                                 (ry)/eps)/eps;

            I2_y -= Processing.vorticity*Pfield_y((rx)/eps,
                                                 (ry)/eps)/eps;

            I1   += Processing.vorticity*exp(-hypot((rx),
                                                   (ry))/eps);
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

    InFlowEvolution.push_back(VortexInMotion(Considered,u_x,u_y,tau));
}


int ViscousVortexDomainSolver::LeakageControl(Vortex Checking)
{
    double NextX = Checking.x;
    double NextY = Checking.y;

    int m;
    if(NextX < 1000 && NextY < 1000 )
    {
        //контроль протекания
        double Rmid = hypot((PanelMids[0][0]-NextX),
                            (PanelMids[0][1]-NextY));
        double projectionDistOnMidVector = (PanelMids[0][0]*NextX+PanelMids[0][1]*NextY)/
                                           (PanelMids[0][0]*PanelMids[0][0]+PanelMids[0][1]*PanelMids[0][1]);

        for(m=1;m<NumberOfPanels;m++)
        {
            double temp = sqrt((PanelMids[m][0]-NextX)*(PanelMids[m][0]-NextX)+
                               (PanelMids[m][1]-NextY)*(PanelMids[m][1]-NextY));
            if(fmin(Rmid,temp)<Rmid)
            {
                Rmid = temp;
                projectionDistOnMidVector = (PanelMids[m][0]*NextX+PanelMids[m][1]*NextY)/
                                            (PanelMids[m][0]*PanelMids[m][0]+PanelMids[m][1]*PanelMids[m][1]);
            }
        }


        if(projectionDistOnMidVector > 1. && isVortexInBody(NextX, NextY))
        {
            return VORTEX_IN_FLOW/*InFlow[i].status*/;
        }
        else
        {
            return VORTEX_IN_BODY;
        }
    }
    else
    {
        return VORTEX_OUT_OF_RANGE;
    }
}


void ViscousVortexDomainSolver::UpdateStreamLines()
{
    //ToDo: move to exernal program

    StreamLinePoint NewPoint1 = StreamLinePoint(0., 1.1);
    StreamLinePoint NewPoint2 = StreamLinePoint(0.,-1.1);

    std::vector <StreamLinePoint>  FirstStreamLineTemp = StreamLine1;
    std::vector <StreamLinePoint> SecondStreamLineTemp = StreamLine2;

    for(unsigned i = 0; i < StreamLine1.size(); i++)
    {
        FirstStreamLineTemp[i].move(tau*xVelocityAt(StreamLine1[i].x,StreamLine1[i].y),
                                    tau*yVelocityAt(StreamLine1[i].x,StreamLine1[i].y));
        StreamLine1[i].move(0.5*tau*(xVelocityAt(StreamLine1[i].x,StreamLine1[i].y)
                                    +xVelocityAt(FirstStreamLineTemp[i].x,FirstStreamLineTemp[i].y)),
                            0.5*tau*(yVelocityAt(StreamLine1[i].x,StreamLine1[i].y)
                                    +yVelocityAt(FirstStreamLineTemp[i].x,FirstStreamLineTemp[i].y)));
    }

    for(unsigned i = 0; i < StreamLine2.size(); i++)
    {
        SecondStreamLineTemp[i].move(tau*xVelocityAt(StreamLine2[i].x,StreamLine2[i].y),
                                     tau*yVelocityAt(StreamLine2[i].x,StreamLine2[i].y));
        StreamLine2[i].move(0.5*tau*(xVelocityAt(StreamLine2[i].x,StreamLine2[i].y)+
                                     xVelocityAt(SecondStreamLineTemp[i].x,SecondStreamLineTemp[i].y)),
                            0.5*tau*(yVelocityAt(StreamLine2[i].x,StreamLine2[i].y)+
                                     yVelocityAt(SecondStreamLineTemp[i].x,SecondStreamLineTemp[i].y)));
    }

    StreamLine1.push_back(NewPoint1);
    StreamLine2.push_back(NewPoint2);
}

void ViscousVortexDomainSolver::Output_StreamLines(const char *FileName)
{
    FILE *ControlPoints;
    ControlPoints = fopen(FileName,"w");
    fprintf(ControlPoints,"X,Y,Z,side\n");

    for(StreamLinePoint i : StreamLine1)
        fprintf(ControlPoints, "%f,%f,%f,%f\n", i.x, i.y, 0., 0.);

    for(StreamLinePoint i : StreamLine2)
        fprintf(ControlPoints, "%f,%f,%f,%f\n", i.x, i.y, 0., 1.);
    fclose(ControlPoints);
}

double ViscousVortexDomainSolver::Qfield_x(double x, double y)
{
    double r = x*x+y*y;
//    return -y/(2*M_PI*fmax(r,1/1024/1024));
    return -y/(2.*M_PI*(r+0.001*0.001));
}

double ViscousVortexDomainSolver::Qfield_y(double x, double y)
{
    double r = x*x+y*y;
    return  x/(2.*M_PI*(r+0.001*0.001));
}

double ViscousVortexDomainSolver::Pfield_x(double x, double y)
{
    double r = hypot(x,y);
    if(r == 0)
        return 1.;
    return x/r*exp(-r);
}

double ViscousVortexDomainSolver::Pfield_y(double x, double y)
{
    double r = hypot(x,y);
    if(r == 0)
        return 1.;
    return y/r*exp(-r);
}

double ViscousVortexDomainSolver::Mfield(double x, double y)
{
    double r2= x*x+y*y;
    double r = sqrt(r2);
    if(x==0 && y==0)
        return 0.;
    return (r+1.)*exp(-r)/r2;
}

int ViscousVortexDomainSolver::isVortexInBody(double x, double y)
{
    if(x*x+y*y<1)
        return 0;
    else
        return 1;
}
