#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"


#ifdef UNBINDING

/* Use AHF method doing Unbinding*/
particle *unbinding_process( particle *HHP, int N, int *RN ) // N=haloNp, RN = N - N_unbound
{
int i,j,k;
double sqa;
double H_of_a;
float ss[3],vs[3];
float center[3],vcenter[3];
double phi0;
int sumN=0;
double totsum=0.0;
double deltar=0.0;
float escapeV=0.0;
float VV;   // particle velocity
int NN=0;
particle *TP;


sqa = sqrt(Time);
H_of_a = Hubble * sqrt(Omega0 / (Time * Time * Time) + (1 - Omega0 - OmegaLambda) / (Time * Time) + OmegaLambda);
for(j = 0; j < 3; j++)
{
     ss[j] = 0.0;
     vs[j] = 0.0;
}
/*Calc center of mass and center of velocity*/
for (i = 0; i < N; i++)
for(j = 0; j < 3; j++)
{
ss[j] += HHP[i].pos[j];
vs[j] += HHP[i].vel[j];
}
for(j = 0; j < 3; j++)
{
    //center[j] = fof_periodic_wrap(ss[j] / N);
    center[j] =  ss[j] / N;
    vcenter[j] = vs[j] / N;
}
#ifdef DEBUG
fprintf(Logfile,"In unbinding_process, cm_pos=%g, %g, %g\n",center[0], center[1], center[2]);
fprintf(Logfile,"In unbinding_process, cm_vel=%g, %g, %g\n",vcenter[0], vcenter[1], vcenter[2]);
fflush(Logfile);
#endif
/*convert physical unit*/
for (i = 0; i < N; i++) 
{
    for (j = 0; j < 3; j++)
    {
	      HHP[i].cpos[j] = fof_periodic( HHP[i].pos[j] - center[j]);  
          HHP[i].cvel[j] = HHP[i].vel[j];         //  Comoving vel in box
     	  HHP[i].pos[j]  = fof_periodic( HHP[i].pos[j] - center[j]);  
    	  HHP[i].pos[j] *= Time;  // convert physical pos
    	  HHP[i].vel[j]  = sqa * (HHP[i].cvel[j] - vcenter[j]) + H_of_a * HHP[i].pos[j];  //add Hubble flow
    }
    HHP[i].rad = sqrt(pow(HHP[i].pos[0],2.0)+pow(HHP[i].pos[1],2.0)+pow(HHP[i].pos[2],2.0));  // need physical r
    HHP[i].Phi = 0.0;
}
qsort(HHP, N, sizeof(particle), rad_sort_particle);  // halo particles sort center to outside
/*Calc potential constant*/
for(i=0;i<N;i++)   
{
sumN ++;
if(i==0)
deltar = fabs( HHP[i].rad );
else
deltar = fabs( HHP[i].rad - HHP[i-1].rad ); 
totsum += sumN / pow(HHP[i].rad, 2.0) * deltar;
}
phi0 = G * PartMass * totsum + G * N * PartMass / HHP[N-1].rad ;  //   omit negative   
#ifdef DEBUG
fprintf(Logfile,"In unbinding_process, phi0=%g\n",phi0);
fflush(Logfile);
#endif
/*Calc every particle potential and escape Velocity*/
sumN=0;
totsum=0.0;
deltar=0.0;
for(k=0,i=0;i<N;i++)   
{
sumN ++;
if(i==0)
deltar = fabs( HHP[i].rad );
else
deltar = fabs( HHP[i].rad - HHP[i-1].rad ); 
totsum += sumN / pow(HHP[i].rad,2.0) * deltar;
HHP[i].Phi = G * PartMass * totsum - phi0;
VV = sqrt(pow(HHP[i].vel[0],2.0)+pow(HHP[i].vel[1],2.0)+pow(HHP[i].vel[2],2.0));    
escapeV = sqrt( 2.0 * fabs( HHP[i].Phi ) ); // escape velocity
if ( VV > escapeV  )
k++;
}

*RN = N - k;  //  remove unbound
 NN = N - k;
#ifdef DEBUG
fprintf(Logfile,"In unbinding_process, N=%d, Nunbound=%d, N200=%d\n",N,k,NN);
fflush(Logfile);
#endif
TP = (particle *)malloc(sizeof(particle)*NN );

for(k=0,i=0;i<N;i++)   
{
VV = sqrt(pow(HHP[i].vel[0],2.0)+pow(HHP[i].vel[1],2.0)+pow(HHP[i].vel[2],2.0));
escapeV = sqrt( 2.0 * fabs( HHP[i].Phi ) );
if ( VV <=  escapeV  )
    {
        for(j=0;j<3;j++)
        {
        TP[k].pos[j] = HHP[i].cpos[j];
        TP[k].vel[j] = HHP[i].cvel[j]; 
        }
#ifdef HALOID
        TP[k].id = HHP[i].id;
#endif
    k++;
    }
}
return TP;
}



particle *Unbinding( particle *HHP, int N, int *Nnew )  // N = halo_Np, Nnew = n200 after Unbinding
{
int i,j,k;
int NNN=0;
int iter = 0;
int sump;
particle *NHP,*HP;
particle *TP;

HP = malloc(N * sizeof(particle));

for(k=0,i=0;i<N;i++)  
{
        for(j=0;j<3;j++)
        {
        HP[k].pos[j] = HHP[i].pos[j];
        HP[k].vel[j] = HHP[i].vel[j];       
       // HP[k].cpos[j] = HHP[i].cpos[j];
       // HP[k].cvel[j] = HHP[i].cvel[j];       
        }
#ifdef HALOID
        HP[k].id = HHP[i].id;
#endif
    k++;
}
sump = N;

while(iter < 10)  //  at least loop 10
{
#ifdef DEBUG
fprintf(Logfile,"----iter = %d----\n",iter);
fflush(Logfile);
#endif
iter++;
NHP=malloc(sump*sizeof(particle));
NHP = unbinding_process( HP, sump, &NNN); // N=haloNp, NNN = N - N_unbound
#ifdef DEBUG
fprintf(Logfile,"In unbinding, N200=%d\n",NNN);
fflush(Logfile);
#endif
free(HP);
HP = malloc(NNN * sizeof(particle));
for(k=0,i=0;i<NNN;i++)  
{
        for(j=0;j<3;j++)
        {
        HP[k].pos[j] = NHP[i].pos[j];
        HP[k].vel[j] = NHP[i].vel[j];       
        }
#ifdef HALOID
        HP[k].id = NHP[i].id;
#endif
    k++;
}

free(NHP);
if(sump-NNN <= 3)   // if remove Np<=3 break loop
    break;
sump = NNN;
}
*Nnew = NNN;
TP = (particle *)malloc(sizeof(particle)*NNN );

for(k=0,i=0;i<NNN;i++)  
{
        for(j=0;j<3;j++)
        {
        TP[k].cpos[j] = HP[i].pos[j];
        TP[k].cvel[j] = HP[i].vel[j];       
        }
#ifdef HALOID
        TP[k].id = HP[i].id;
#endif
    k++;
}
return TP;
}
#endif
