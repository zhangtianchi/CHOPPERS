#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"


/*Use Bullock 2001 paper calculate spin parameter */
float bspin( int gr, float sam[], float vel[], float cm[], float *vdisp ) //  return spin and angular momentum, halo velocity and center of mass, velocity dispersion
{
int i,j,k;
float sqa;
float H_of_a;
float ss[3],vs[3],r200;
int N200=0;
float center[3],vcenter[3];
float lx,ly,lz;
float hspin;
float dx[3],dv[3];
float xcenter,ycenter,zcenter;
float vd = 0.0;

lx = ly = lz = 0;
sqa = sqrt(Time);
r200 = Halo_R_Crit200[gr];
N200 = np;
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
    cm[j] = center[j];
    vcenter[j] =  vs[j] / N200;
    vel[j] = vcenter[j] * sqa ;
}

    for (i = 0; i < N200; i++) 
    {
      for(j = 0; j < 3; j++)
	{
	  dx[j] = Time * fof_periodic( data[i].cpos[j] - center[j] );
	  dv[j] = sqa * ( data[i].cvel[j] - vcenter[j] );
	  dv[j] += H_of_a * dx[j];   // no use because dx X dv = dx X (H*dx+dv) , dx X dx =0 ，so this term =0
      vd += dv[j] *dv[j];
	}
      lx += dx[1] * dv[2] - dx[2] * dv[1];
      ly += dx[2] * dv[0] - dx[0] * dv[2];
      lz += dx[0] * dv[1] - dx[1] * dv[0];
    }
  *vdisp = sqrt( vd / (3.0 * N200 ) );
//special AM
  xcenter = lx / N200;
  ycenter=  ly / N200;
  zcenter = lz / N200;
  sam[0] = xcenter;
  sam[1] = ycenter;
  sam[2] = zcenter;
// spin
hspin = sqrt(xcenter*xcenter+ycenter*ycenter+zcenter*zcenter) / ( sqrt(2.0) * sqrt(G*PartMass*N200 / r200) * r200);
return hspin;
}

#ifdef OUTPE
/*Use Peebles 1969 paper calculate spin parameter */
float pspin( int gr, float ke, float pe ) 
{
int i,j,k;
float sqa;
float H_of_a;
float ss[3],vs[3],r200;
int N200=0;
float center[3],vcenter[3];
float lx,ly,lz;
float hspin;
float dx[3],dv[3];
float xcenter,ycenter,zcenter;
float m200 = 0.0;

lx = ly = lz = 0;
sqa = sqrt(Time);
r200 = Halo_R_Crit200[gr];
N200 = np;
m200 = N200 * PartMass;
H_of_a = Hubble * sqrt(Omega0 / (Time * Time * Time) + (1 - Omega0 - OmegaLambda) / (Time * Time) + OmegaLambda);

for(j = 0; j < 3; j++)
{
     ss[j] = 0.0;
     vs[j] = 0.0;
}
for (i = 0; i < N200; i++)
for (j = 0; j < 3; j++)
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
    {
      for(j = 0; j < 3; j++)
	{
	  dx[j] = Time * fof_periodic( data[i].cpos[j] - center[j] );
	  dv[j] = sqa * ( data[i].cvel[j] - vcenter[j] );
	  dv[j] += H_of_a * dx[j];   // no use because dx X dv = dx X (H*dx+dv) , dx X dx =0 ，so this term =0
	}
      lx += PartMass * ( dx[1] * dv[2] - dx[2] * dv[1] );
      ly += PartMass * ( dx[2] * dv[0] - dx[0] * dv[2] );
      lz += PartMass * ( dx[0] * dv[1] - dx[1] * dv[0] );
    }
//special AM
  xcenter = lx;
  ycenter=  ly;
  zcenter = lz;
 // xcenter = lx / N200;
 // ycenter=  ly / N200;
 // zcenter = lz / N200;
// spin
hspin = sqrt(xcenter*xcenter+ycenter*ycenter+zcenter*zcenter) * sqrt( fabs( ke * m200 + pe * m200 ) ) / ( G * pow( m200, 5.0 / 2.0  ) )  ;
return hspin;
}
#endif


/*
#ifdef DEBUG
fprintf(Logfile,"In spin.c, N200=%d, r200=%g, datapos=%g,%g,%g\n",N200, r200,data[0].cpos[0],data[0].cpos[1],data[0].cpos[2]);
fprintf(Logfile,"In spin.c, cmpos=%g,%g,%g\n",center[0],center[1],center[2]);
fprintf(Logfile,"In spin.c, cmvel=%g,%g,%g\n",vcenter[0],vcenter[1],vcenter[2]);
fprintf(Logfile,"In spin.c, SAM=%g,%g,%g,am=%g,spin=%g\n",xcenter,ycenter,zcenter,sqrt(xcenter*xcenter+ycenter*ycenter+zcenter*zcenter),hspin);
fflush(Logfile);
#endif
*/
