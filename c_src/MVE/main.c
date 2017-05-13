#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <gsl/gsl_linalg.h>

/*  Code from PGRM as samples for GSL
 *
 *
 *
    gsl_matrix_memcpy	(general_system, stream_function.sys);
    Task->sys = gsl_matrix_alloc (N*N,N*N);
    Task->rightpart	= gsl_vector_alloc(N*N);

void form_system_t (task *Task)
{
    int i, j;
    basis_args args;
    args.x = 0.;
    args.y = 0.;
    args.m = 0;
    args.n = 0;

    #pragma omp parallel for shared(Task,N) private(i,j) firstprivate(args)
    for(i = 0; i < N*N; i++)
    {
        args.m = i;
        gsl_vector_set(Task->rightpart, i, gauss_integral3(right_under_int_t,Task->area,args,2,Task));
        for(j = 0; j < N*N; j++)
        {
            args.n = j;
            gsl_matrix_set(Task->sys, i,j, gauss_integral3(left_under_int_t,Task->area,args,2,Task));
        }
    }
}

void form_right_part_t (task *Task)
{
    int i;//, j;
    basis_args args;
    args.x = 0.;
    args.y = 0.;
    args.m = 0;
    args.n = 0;

    #pragma omp parallel for shared(Task,N) private(i) firstprivate(args)
    for(i = 0; i < N*N; i++)
    {
        args.m = i;
        gsl_vector_set(Task->rightpart, i, gauss_integral3(right_under_int_t,Task->area,args,2,Task));
    }
}

void solve_task(task *Task)
{
    int i;
    gsl_permutation * p = gsl_permutation_alloc (N*N);
    gsl_linalg_LU_decomp (Task->sys,p,&i);
    gsl_linalg_LU_solve (Task->sys,p,Task->rightpart,Task->solution);
}
 *
 *
 *
 *
 *
 */

typedef struct Vortex
{
    double x, y;
    double vorticity;

    int active;
} Vortex;


/* possibly should be replaced with one function? */
double Qfield_x(double x, double y)
{
    double r = x*x+y*y;
    return -y/(2*M_PI*fmax(r,0.001));
}

double Qfield_y(double x, double y)
{
    double r = x*x+y*y;
    return  x/(2*M_PI*fmax(r,0.001));
}

double Pfield_x(double x, double y)
{
    double r = sqrt(x*x+y*y);
    if(r == 0)
        return 0;
    return x/r*exp(-r);
}

double Pfield_y(double x, double y)
{
    double r = sqrt(x*x+y*y);
    if(r == 0)
        return 0;
    return y/r*exp(-r);
}


double Mfield(double x, double y)
{
    double r2= x*x+y*y;
    double r = sqrt(r2);
    if(x==0 && y==0)
        return 0.;
    return (r+1.)*exp(-r)/r2;
}

int VortexInBody(double x, double y)
{
    if(x*x+y*y<1)
        return 1;
    else
        return 0;
}

int main(int argc, char *argv[])
{
    int Np = 40, MaxVortexes = 6000, Niterations = 120,      i,j,k, iterator;
    double nu = 0.0001, tau = 0.001, eps = 0.0001, rho = 1.;
    double vx_inf = 1., vy_inf = 0.;

    double  PanelNodes [Np][2],
            PanelMids  [Np][2],
            PanelNorms [Np][2],
            PanelTaus  [Np][2],
            PanelLength[Np];

    double DeployPoints[Np][2];

    Vortex InFlow[MaxVortexes], NextInFlow[MaxVortexes];

    int ActiveVortexesInFLow=0;



    if(1)/* cylinder profile forming & deploys */
    {

        for(i=0;i<Np;i++)
        {
            PanelNodes[i][0] = cos(2.*M_PI/Np*i);
            PanelNodes[i][1] = sin(2.*M_PI/Np*i);

            DeployPoints[i][0] = 1.02*PanelNodes[i][0];
            DeployPoints[i][1] = 1.02*PanelNodes[i][1];
        }

        for(i=0;i<Np-1;i++)
        {
            PanelMids[i][0] = (PanelNodes[i][0]+PanelNodes[i+1][0])/2;
            PanelMids[i][1] = (PanelNodes[i][1]+PanelNodes[i+1][1])/2;

            PanelTaus[i][0] = (PanelNodes[i+1][0]-PanelNodes[i][0]);
            PanelTaus[i][1] = (PanelNodes[i+1][1]-PanelNodes[i][1]);
            PanelLength[i]  = sqrt(PanelTaus[i][0]*PanelTaus[i][0]+
                                   PanelTaus[i][1]*PanelTaus[i][1]);
            PanelTaus[i][0] /= PanelLength[i];
            PanelTaus[i][1] /= PanelLength[i];
            PanelNorms[i][0] = -PanelTaus[i][1];
            PanelNorms[i][1] =  PanelTaus[i][0];
        }
        PanelMids[Np-1][0] = (PanelNodes[Np-1][0]+PanelNodes[0][0])/2;
        PanelMids[Np-1][1] = (PanelNodes[Np-1][1]+PanelNodes[0][1])/2;
        PanelTaus[Np-1][0] = (PanelNodes[0][0]-PanelNodes[Np-1][0]);
        PanelTaus[Np-1][1] = (PanelNodes[0][1]-PanelNodes[Np-1][1]);
        PanelLength[Np-1]  = sqrt(PanelTaus[Np-1][0]*PanelTaus[Np-1][0]+
                               PanelTaus[Np-1][1]*PanelTaus[Np-1][1]);
        PanelTaus[Np-1][0] /= PanelLength[Np-1];
        PanelTaus[Np-1][1] /= PanelLength[Np-1];
        PanelNorms[Np-1][0] = -PanelTaus[Np-1][1];
        PanelNorms[Np-1][1] =  PanelTaus[Np-1][0];
    }

    // defining and completing Vortex Birth Matrix
    gsl_matrix *VBM = gsl_matrix_alloc(Np+1,Np+1);
    for(i=0; i<Np; i++)
    {
        for(j=0; j<Np; j++)
        {
            gsl_matrix_set(VBM,i,j,
            PanelNorms[i][0]*Qfield_x(PanelMids[i][0]-DeployPoints[j][0],
                                      PanelMids[i][1]-DeployPoints[j][1])
           +PanelNorms[i][1]*Qfield_y(PanelMids[i][0]-DeployPoints[j][0],
                                      PanelMids[i][1]-DeployPoints[j][1]));
        }
    }
    for(i=0;i<Np;i++)
    {
        gsl_matrix_set(VBM,Np,i,1);
        gsl_matrix_set(VBM,i,Np,1);
    }
    gsl_matrix_set(VBM,Np,Np,0);



    // recieve initial aproximation for vortexes on surface
    gsl_vector *VOI = gsl_vector_alloc(Np+1); //Velocity On Infinity
    gsl_vector *VFV = gsl_vector_alloc(Np+1); //Velocity from active vortexes

    gsl_vector *velocities_on_surface = gsl_vector_alloc(Np+1);

    gsl_vector *NewVorticities = gsl_vector_alloc(Np+1);

//    gsl_vector_add()
    for(i=0;i<Np;i++)
    {
        gsl_vector_set(VOI,i,PanelNorms[i][0]*vx_inf+PanelNorms[i][1]*vy_inf);
    }
    gsl_vector_set(VOI,Np,0.);


    gsl_permutation * p = gsl_permutation_alloc (Np+1);
    gsl_linalg_LU_decomp (VBM,p,&k);
    gsl_linalg_LU_solve (VBM,p,VOI,NewVorticities);


    for(iterator=0; iterator<Niterations; iterator++)
    {

        for(i=0;i<Np;i++)
        {
            InFlow[ActiveVortexesInFLow+i].active = 1;
            InFlow[ActiveVortexesInFLow+i].vorticity = gsl_vector_get(NewVorticities,i);
            InFlow[ActiveVortexesInFLow+i].x = DeployPoints[i][0];
            InFlow[ActiveVortexesInFLow+i].y = DeployPoints[i][1];
        }
        ActiveVortexesInFLow += Np;

        // time iteration
        double u_x, u_y;
        printf("active vortexes %d\n",ActiveVortexesInFLow);
        for(i=0;i<ActiveVortexesInFLow;i++)
        {
            double  I0 = 0., I1 = 0.,
                    I2_x = 0., I2_y = 0.,
                    I3_x = 0., I3_y = 0.;

            u_x = vx_inf;
            u_y = vy_inf;
            //Curl component of velocity
            for(j=0;j<ActiveVortexesInFLow;j++)
            {
                if( InFlow[j].active==1 && (i-j)!=0 )
                {
                    u_x += InFlow[j].vorticity*Qfield_x(InFlow[i].x-InFlow[j].x,InFlow[i].y-InFlow[j].y);
                    u_y += InFlow[j].vorticity*Qfield_y(InFlow[i].x-InFlow[j].x,InFlow[i].y-InFlow[j].y);
                }
            }
            //Viscous component of velocity
            for(j=0;j<ActiveVortexesInFLow;j++)
            {
                if( InFlow[j].active==1 )
                {
                    I2_x -= InFlow[j].vorticity*Pfield_x((InFlow[i].x-InFlow[j].x)/eps,(InFlow[i].y-InFlow[j].y)/eps)/eps;
                    I2_y -= InFlow[j].vorticity*Pfield_y((InFlow[i].x-InFlow[j].x)/eps,(InFlow[i].y-InFlow[j].y)/eps)/eps;
                    I1 += InFlow[j].vorticity*exp(-sqrt((InFlow[i].x-InFlow[j].x)*(InFlow[i].x-InFlow[j].x)
                                                       +(InFlow[i].y-InFlow[j].y)*(InFlow[i].y-InFlow[j].y))/eps);
                }
            }
            I0 = 2*M_PI*eps*eps;
            for(j=0;j<Np;j++)
            {
                if( InFlow[j].active==1 )
                {
                    I3_x += PanelNorms[j][0]*exp(-sqrt((InFlow[i].x-PanelMids[j][0])*(InFlow[i].x-PanelMids[j][0])
                            +(InFlow[i].y-PanelMids[j][1])*(InFlow[i].y-PanelMids[j][1]))/eps)*PanelLength[j];
                    I3_y += PanelNorms[j][1]*exp(-sqrt((InFlow[i].x-PanelMids[j][0])*(InFlow[i].x-PanelMids[j][0])
                            +(InFlow[i].y-PanelMids[j][1])*(InFlow[i].y-PanelMids[j][1]))/eps)*PanelLength[j];
                    I0 += ((InFlow[i].x-PanelMids[j][0])*PanelNorms[j][0]+(InFlow[i].y-PanelMids[j][1])*PanelNorms[j][1])
                            *PanelLength[j]*Mfield((InFlow[i].x-PanelMids[j][0])/eps,(InFlow[i].y-PanelMids[j][1])/eps);
                }
            }

            //        printf("%f %f\n",u_x,u_y);
            u_x += nu*(-I2_x/I1+I3_x/I0);
            u_y += nu*(-I2_y/I1+I3_y/I0);
            //        printf("%f %f\n",u_x,u_y);

            if(tau*u_x < 100 && tau*u_y < 100)
            {
                NextInFlow[i].x = InFlow[i].x + tau*u_x;
                NextInFlow[i].y = InFlow[i].y + tau*u_y;
                if(!VortexInBody(NextInFlow[i].x,NextInFlow[i].y))
                {
                    NextInFlow[i].active = InFlow[i].active;
                    NextInFlow[i].vorticity = InFlow[i].vorticity;
                }
                else
                {
                    NextInFlow[i].active = 3;
                    NextInFlow[i].vorticity = InFlow[i].vorticity;
                }
            }
            else
            {
                NextInFlow[i].x = InFlow[i].x;
                NextInFlow[i].y = InFlow[i].y;
                NextInFlow[i].active = 2;
                NextInFlow[i].vorticity = InFlow[i].vorticity;
            }
        }
        //updating coords
        for(i=0;i<ActiveVortexesInFLow;i++)
        {
            InFlow[i] = NextInFlow[i];
        }

        //definition of speed from vortexes after iteration
        double ExtSpeedInMiddles[Np];
        for(j=0; j<Np; j++)
            ExtSpeedInMiddles[j] = 0;
        gsl_vector_memcpy(velocities_on_surface,VOI);

        for(i=0; i<ActiveVortexesInFLow; i++)
        {
            if(InFlow[i].active==1)
            {
                for(j=0; j<Np; j++)
                {
                    ExtSpeedInMiddles[j] +=
                            (PanelNorms[j][0]*Qfield_x(PanelMids[j][0]-InFlow[i].x,
                            PanelMids[j][1]-InFlow[i].y)
                            +PanelNorms[j][1]*Qfield_y(PanelMids[j][0]-InFlow[i].x,
                            PanelMids[j][1]-InFlow[i].y))
                            *InFlow[i].vorticity;
                }
            }
        }

        for(j=0; j<Np; j++)
            gsl_vector_set(VFV,j,ExtSpeedInMiddles[j]);
        gsl_vector_add(velocities_on_surface,VFV);
        gsl_vector_set(velocities_on_surface,Np,0);

        gsl_linalg_LU_decomp (VBM,p,&k);
        gsl_linalg_LU_solve (VBM,p,velocities_on_surface,NewVorticities);

    }

    //the most interesting part: output of results
    FILE *VortexesOut, *BodyOut;
    VortexesOut = fopen("MVE_flow.dat","w");
    BodyOut =     fopen("MVE_body.dat","w");

    for(i=0;i<ActiveVortexesInFLow;i++)
    {
        if(InFlow[i].active==1)
            fprintf(VortexesOut,"%f %f %f\n",InFlow[i].x,InFlow[i].y,InFlow[i].vorticity);
    }

    for(i=0; i<Np; i++)
        fprintf(BodyOut,"%f %f\n", DeployPoints[i][0], DeployPoints[i][1]);

    fclose(VortexesOut);
    fclose(BodyOut);
    return 0;
}
