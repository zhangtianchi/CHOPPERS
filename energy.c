#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"


/* Calculate halo kinetic energy */
float kinetic( void ) 
{
int i,j,k;
float sqa, ke=0.0;
float H_of_a;
float ss[3],vs[3],r200;
int N200=0;
float center[3],vcenter[3];
float dx[3],dv[3];
float *vd;
float m200=0.0;

sqa = sqrt(Time);
N200 = np;
vd = malloc( sizeof( double ) * N200 );    
memset(vd, 0.0, sizeof(double) * N200);  // initialization
m200 = N200 * PartMass;
H_of_a = Hubble * sqrt(Omega0 / (Time * Time * Time) + (1 - Omega0 - OmegaLambda) / (Time * Time) + OmegaLambda);

for(j = 0; j < 3; j++)
{
     ss[j] = 0.0;
     vs[j] = 0.0;
}
for (i = 0; i < N200; i++)
for(j = 0; j < 3; j++)
{
ss[j] += data[i].cpos[j];
vs[j] += data[i].cvel[j];
}
for(j = 0; j < 3; j++)
{
    center[j]  = fof_periodic_wrap(ss[j] / N200) ;
    vcenter[j] =  vs[j] / N200;
}

    for (i = 0; i < N200; i++) 
    for (j = 0; j < 3; j++)
	{
	  dx[j] = Time * fof_periodic( data[i].cpos[j] - center[j] );
	  dv[j] = sqa * ( data[i].cvel[j] - vcenter[j] );
	  dv[j] += H_of_a * dx[j];   // no use because dx X dv = dx X (H*dx+dv) , dx X dx =0 ，so this term =0
      vd[i] += 0.5 * PartMass * dv[j] * dv[j];
	}
    for (i = 0; i < N200; i++) 
    ke += vd[i]; 

free(vd);

return ke / m200;
}

#ifdef OUTPE 
/* Use PP method compute halo potential energy, When Np > 10^6 can cause more time, change to Octree method in future... */
float potential( void )
{
int i, j, k, N;
double r12=0.0, pe=0.0, *ppe;
double softening = 0.0;
float m200=0.0;
N = np;
m200 = N * PartMass;
ppe = malloc( sizeof( double ) * N );    
memset(ppe, 0.0, sizeof(double) * N);  // initialization

for(i=0;i<N;i++)
for(j=0;j<N;j++) 
{
if(i==j) continue;
r12 = sqrt(pow(fof_periodic(data[i].cpos[0]-data[j].cpos[0]),2.0)+pow(fof_periodic(data[i].cpos[1]-data[j].cpos[1]),2.0)+pow(fof_periodic(data[i].cpos[2]-data[j].cpos[2]),2.0));
ppe[i] += -1.0* G * PartMass * PartMass / sqrt( pow(r12,2.0)+pow(softening,2.0) );
}

for(i=0;i<N;i++)
    pe += ppe[i];

free(ppe);

return pe / m200;
}
#endif


