#include "viscousvortexdomainsolver.h"

ViscousVortexDomainSolver::ViscousVortexDomainSolver()
{
    NumberOfPanels = 100;
    MaxVortexes = 16500;
    Niterations = 250;
    ActiveVortexesInFLow=0;

    nu = 1/13;
    tau = 0.001;
    rho = 1.;
    vx_inf = 1.;
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
//    delete [] OriginalVortexGenerationMatrix;
//    delete [] VortexGenerationMatrix;

    gsl_vector_free(VelocityProjections);
    gsl_vector_free(NewVorticities);
//    delete [] VelocityProjections;
//    delete [] NewVorticities;

    gsl_permutation_free(Permutation);

//    delete [] InFlow;
//    delete [] NextInFlow;
}

void ViscousVortexDomainSolver::Solve()
{
    int i,k,m, iterator;

    DivideProfileToPanels();
    CompletingGeneratingMatrix();

    // recieve initial aproximation for vortexes on surface
    for(i=0;i<NumberOfPanels;i++)
    {
        gsl_vector_set(VelocityProjections,i,-(PanelNorms[i][0]*vx_inf+PanelNorms[i][1]*vy_inf));
    }
    gsl_vector_set(VelocityProjections,NumberOfPanels,0.);

    gsl_linalg_LU_decomp(VortexGenerationMatrix,Permutation,&k);
    gsl_linalg_LU_solve (VortexGenerationMatrix,Permutation,VelocityProjections,NewVorticities);


    FILE *forces;
    forces = fopen("forces.dat","w");

    for(iterator=0; iterator<Niterations; iterator++)
    {
        double f_x = 0., f_y = 0./*, M_z = 0.*/;

        //add generated vortexes to flow
        for(i=0;i<NumberOfPanels;i++)
        {
            Vortex NewVortexFromPanel;
            NewVortexFromPanel.x = DeployPoints[i][0];
            NewVortexFromPanel.y = DeployPoints[i][1];
            NewVortexFromPanel.status = VORTEX_IN_FLOW;
            NewVortexFromPanel.vorticity = gsl_vector_get(NewVorticities,i);

            InFlow.push_back(NewVortexFromPanel);
        }
        ActiveVortexesInFLow = InFlow.size();


        // time iteration
        printf("active vortexes \t%d,\t iteration\t%d\n",ActiveVortexesInFLow,iterator);

        UpdateVotexPositions();

        for(m = 0; m < NumberOfPanels; m++)
        {
            f_x += -rho*(gsl_vector_get(NewVorticities,m)/tau* PanelMids[m][1]);
            f_y +=  rho*(gsl_vector_get(NewVorticities,m)/tau* PanelMids[m][0]);
        }

        fprintf(forces,"%f %f %f\n",iterator*tau, f_x, f_y);

        //updating coords
        InFlow = NextInFlow;

        UpdateAssociatedVortexes();
        if(iterator%3==0 || 1)
        {
            char FileName[32];
            sprintf(FileName,"VorticityField%6.6d.csv",iterator);
            Output_ParaView_Line(FileName);
        }
    }
}

void ViscousVortexDomainSolver::DivideProfileToPanels()
{
    int i;
    for(i=0;i<NumberOfPanels;i++)
    {
        PanelNodes[i][0] = cos(2.*M_PI/NumberOfPanels*i);
        PanelNodes[i][1] = sin(2.*M_PI/NumberOfPanels*i);

        DeployPoints[i][0] = 1.05*PanelNodes[i][0];
        DeployPoints[i][1] = 1.05*PanelNodes[i][1];

        PanelNorms[i][0] = cos(2.*M_PI/NumberOfPanels*(i+0.5));
        PanelNorms[i][1] = sin(2.*M_PI/NumberOfPanels*(i+0.5));

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
    //(откуда ноги растут? POLARA?)
    int j;
    if(ActiveVortexesInFLow>NumberOfPanels)
    {
        double sq1=InitialEpsilon*InitialEpsilon;
        double sq2=sq1, sq3=sq1;
        for(j=0;j<ActiveVortexesInFLow;j++)
        {
            double temp = (InFlow[i].x-InFlow[j].x)*(InFlow[i].x-InFlow[j].x)
                         +(InFlow[i].y-InFlow[j].y)*(InFlow[i].y-InFlow[j].y);
            if(
                    i!=j && temp < sq1 &&
                    InFlow[j].status==1 &&
                    InFlow[i].vorticity*InFlow[j].vorticity>0
               )
            {
                sq1 = temp;
                sq2 = sq1;
                sq3 = sq2;
            }
        }
        InitialEpsilon = sqrt( (sq1+sq2+sq3)/3 );
    }
    return InitialEpsilon;
}

int ViscousVortexDomainSolver::LeakageControl(double u_x, int i, double NextY, double u_y, double NextX)
{
    int m;
    if(tau*u_x < 1000 && tau*u_y < 1000 )
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


        if(projectionDistOnMidVector > 1.)
        {
            return InFlow[i].status;
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

void ViscousVortexDomainSolver::UpdateVotexPositions()
{
    int i;
    NextInFlow.clear();
    /*
    for(i=0;i<ActiveVortexesInFLow;i++)
    {
        double eps = sqrt(4/3)*PanelLength[0];
        double  I0 = 2*M_PI*eps*eps,
                I1 = 0.,
                I2_x = 0., I2_y = 0.,
                I3_x = 0., I3_y = 0.;

        double u_x = vx_inf;
        double u_y = vy_inf;


        eps = UpdateEpsilon(i, eps);

        //Vortex-caused part of velocity
        for(j=0;j<ActiveVortexesInFLow;j++)
        {
            double rx = InFlow[i].x-InFlow[j].x;
            double ry = InFlow[i].y-InFlow[j].y;
            if( InFlow[j].status==1 &&  i!=j )
            {
                u_x += InFlow[j].vorticity*Qfield_x(rx,ry);
                u_y += InFlow[j].vorticity*Qfield_y(rx,ry);
            }
        }


        //Viscous component of velocity
        for(j=0;j<ActiveVortexesInFLow;j++)
        {
            if( InFlow[j].status==1 )
            {
                I2_x += InFlow[j].vorticity*Pfield_x((InFlow[i].x-InFlow[j].x)/eps,(InFlow[i].y-InFlow[j].y)/eps)/eps;
                I2_y += InFlow[j].vorticity*Pfield_y((InFlow[i].x-InFlow[j].x)/eps,(InFlow[i].y-InFlow[j].y)/eps)/eps;
                I1  += InFlow[j].vorticity*exp(-hypot((InFlow[i].x-InFlow[j].x),(InFlow[i].y-InFlow[j].y))/eps);
            }
        }
        for(j=0;j<NumberOfPanels;j++)
        {
            double Rr = hypot((InFlow[i].x-PanelMids[j][0]),(InFlow[i].y-PanelMids[j][1]));
            I3_x -= PanelNorms[j][0]*exp(-Rr/eps)*PanelLength[j];
            I3_y -= PanelNorms[j][1]*exp(-Rr/eps)*PanelLength[j];
            I0   -= ((InFlow[i].x-PanelMids[j][0])*PanelNorms[j][0]
                    +(InFlow[i].y-PanelMids[j][1])*PanelNorms[j][1])
                    *PanelLength[j]
                    *Mfield((InFlow[i].x-PanelMids[j][0])/eps,
                    (InFlow[i].y-PanelMids[j][1])/eps);

        }

        u_x += nu*(-I2_x/I1+I3_x/I0);
        u_y += nu*(-I2_y/I1+I3_y/I0);

        double NextX = InFlow[i].x + tau*u_x;
        double NextY = InFlow[i].y + tau*u_y;
        double NextVorticity = InFlow[i].vorticity;

        int NextStatus = LeakageControl(u_x, i, NextY, u_y, NextX);
        Vortex NextVortex;
        NextVortex.status = NextStatus;
        NextVortex.vorticity = NextVorticity;
        NextVortex.x = NextX;
        NextVortex.y = NextY;
        if(NextStatus == VORTEX_IN_FLOW)
            NextInFlow.push_back(NextVortex);
        if(NextStatus == VORTEX_IN_BODY)
            InBodyVortexes.push_back(NextVortex);
    }
    */

    int ThreadNum = 8;

    std::vector<std::thread> ThrPool;

    for(i=0;i<ThreadNum;i++)
    {
        ThrPool.push_back(std::thread(ThreadCrutch_UpdateVotexPositions,this,i,ThreadNum));
    }

    for(i=0;i<ThreadNum;i++)
    {
        ThrPool[i].join();
    }


}

void ViscousVortexDomainSolver::UpdateAssociatedVortexes()
{
    int i,j,k;

    double VelocityProjInControlPoint[NumberOfPanels];
    for(j=0; j<NumberOfPanels; j++)
        VelocityProjInControlPoint[j] = PanelNorms[j][0]*vx_inf+
                                        PanelNorms[j][1]*vy_inf;

    for(i=0; i<ActiveVortexesInFLow; i++)
    {
        if(InFlow[i].status==1)
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
            for(n=0;n<ActiveVortexesInFLow;n++)
            {
                if(InFlow[n].status==1)
                {
                    vx += InFlow[n].vorticity*Qfield_x((0.04*i-2)-InFlow[n].x,
                                                       (0.04*j-2)-InFlow[n].y);
                    vy += InFlow[n].vorticity*Qfield_y((0.04*i-2)-InFlow[n].x,
                                                       (0.04*j-2)-InFlow[n].y);
                }
            }
            fprintf(VelocityField,"%f %f %f %f\n",0.04*i-2, 0.04*j-2, vx, vy);
        }
    }
    fclose(VelocityField);
}

void ViscousVortexDomainSolver::Output_ParaView_Line(const char *FileName)
{
    int n;
    FILE *VelocityField;
    VelocityField = fopen(FileName,"w");
    fprintf(VelocityField,"X,Y,Z,vorticity\n");
    for(n=0;n<ActiveVortexesInFLow;n++)
    {
        if(InFlow[n].status==1)
            fprintf(VelocityField,"%f,%f,%f,%f\n",InFlow[n].x,InFlow[n].y,0.,InFlow[n].vorticity);

    }
    fclose(VelocityField);
}

void ViscousVortexDomainSolver::Thread_UpdateVortexPosition(int ID, int ThreadNum)
{

    int i,j;

    int H = ActiveVortexesInFLow/ThreadNum;
    int IStart  =   ID*H;
    int IEnd    =   IStart + H;
    if(ID==(ThreadNum-1))
        IEnd = ActiveVortexesInFLow;
    for(i=IStart;i<IEnd;i++)
    {
        double eps = sqrt(4/3)*PanelLength[0];
        double  I0 = 2*M_PI*eps*eps,
                I1 = 0.,
                I2_x = 0., I2_y = 0.,
                I3_x = 0., I3_y = 0.;

        double u_x = vx_inf;
        double u_y = vy_inf;


        eps = UpdateEpsilon(i, eps);

        //Vortex-caused part of velocity
        for(j=0;j<ActiveVortexesInFLow;j++)
        {
            double rx = InFlow[i].x-InFlow[j].x;
            double ry = InFlow[i].y-InFlow[j].y;
            if( InFlow[j].status==1 &&  i!=j )
            {
                u_x += InFlow[j].vorticity*Qfield_x(rx,ry);
                u_y += InFlow[j].vorticity*Qfield_y(rx,ry);
            }
        }


        //Viscous component of velocity
        for(j=0;j<ActiveVortexesInFLow;j++)
        {
            if( InFlow[j].status==1 )
            {
                I2_x += InFlow[j].vorticity*Pfield_x((InFlow[i].x-InFlow[j].x)/eps,(InFlow[i].y-InFlow[j].y)/eps)/eps;
                I2_y += InFlow[j].vorticity*Pfield_y((InFlow[i].x-InFlow[j].x)/eps,(InFlow[i].y-InFlow[j].y)/eps)/eps;
                I1  += InFlow[j].vorticity*exp(-hypot((InFlow[i].x-InFlow[j].x),(InFlow[i].y-InFlow[j].y))/eps);
            }
        }
        for(j=0;j<NumberOfPanels;j++)
        {
            double Rr = hypot((InFlow[i].x-PanelMids[j][0]),(InFlow[i].y-PanelMids[j][1]));
            I3_x -= PanelNorms[j][0]*exp(-Rr/eps)*PanelLength[j];
            I3_y -= PanelNorms[j][1]*exp(-Rr/eps)*PanelLength[j];
            I0   -= ((InFlow[i].x-PanelMids[j][0])*PanelNorms[j][0]
                    +(InFlow[i].y-PanelMids[j][1])*PanelNorms[j][1])
                    *PanelLength[j]
                    *Mfield((InFlow[i].x-PanelMids[j][0])/eps,
                    (InFlow[i].y-PanelMids[j][1])/eps);

        }

        u_x += nu*(-I2_x/I1+I3_x/I0);
        u_y += nu*(-I2_y/I1+I3_y/I0);

        double NextX = InFlow[i].x + tau*u_x;
        double NextY = InFlow[i].y + tau*u_y;
        double NextVorticity = InFlow[i].vorticity;

        int NextStatus = LeakageControl(u_x, i, NextY, u_y, NextX);
        Vortex NextVortex;
        NextVortex.status = NextStatus;
        NextVortex.vorticity = NextVorticity;
        NextVortex.x = NextX;
        NextVortex.y = NextY;

        NextInFlow_WriteMutex.lock();
        if(NextStatus == VORTEX_IN_FLOW)
            NextInFlow.push_back(NextVortex);
        if(NextStatus == VORTEX_IN_BODY)
            InBodyVortexes.push_back(NextVortex);
        NextInFlow_WriteMutex.unlock();
    }


}

void ViscousVortexDomainSolver::ThreadCrutch_UpdateVotexPositions(ViscousVortexDomainSolver *Task, int ID, int ThreadNum)
{
    Task->Thread_UpdateVortexPosition(ID, ThreadNum);
}

double ViscousVortexDomainSolver::Qfield_x(double x, double y)
{
    double r = x*x+y*y;
    return -y/(2*M_PI*fmax(r,1/1024/1024));
}

double ViscousVortexDomainSolver::Qfield_y(double x, double y)
{
    double r = x*x+y*y;
    return  x/(2*M_PI*fmax(r,1/1024/1024));
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

int ViscousVortexDomainSolver::VortexInBody(double x, double y)
{
    if(x*x+y*y<1)
        return 1;
    else
        return 0;
}
