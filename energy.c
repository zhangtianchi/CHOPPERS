#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include<sys/types.h>
#include<sys/timeb.h>
#include "allvars.h"
#include "proto.h"


#define Nr 1000        //  Calc Potential use 1000 random particle inside r200, Then weighted


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
memset(vd, 0.0, sizeof(double) * N200); 
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
	  dv[j] += H_of_a * dx[j];   
      vd[i] += 0.5 * dv[j] * dv[j];

	}
    for (i = 0; i < N200; i++) 
    ke += vd[i]; 

free(vd);

return ke;
}

#ifdef OCTREE  

#define TreeAllocFactor 2.0
float Potential( void )
{
  int i,j,Nnode,allbytes;
  double pott,sum=0.0, m200;
  m200 = NumPart * PartMass;

  P = malloc(sizeof(struct potdata)*NumPart);
  memset(P, 0, sizeof(struct potdata) * NumPart);
    for(i=0;i<NumPart;i++)
    for(j=0;j<3;j++)
    P[i].Pos[j] = data[i].cpos[j];

  allbytes = tree_treeallocate(TreeAllocFactor * NumPart, NumPart);


  Nnode = tree_treebuild();

#ifdef DEBUG
fprintf(Logfile,"Allocate memory for %d tree-(%d) nodes (%g MB).\n", (int)(TreeAllocFactor*NumPart), Nnode, allbytes / (1024.0 * 1024.0));
fflush(Logfile);
#endif

  for(i = 0; i < NumPart; i++)
    {
  pott = tree_treeevaluate_potential(i);
    sum += pott;
    }

    free(P);
tree_treefree();



return sum / 2.0 ;  //  twice potential need /2.0
}
#endif

/* Use PP method and SPH kernel compute halo potential energy, calc Potential use 1000 random particle inside r200, Then weighted */
float potential( void )
{
int *a,w,t;
int i, j, k, N;
double r=0.0, pot=0.0;
double picka;
float m200=0.0;
float h, mass, u, wp;
int Nrand;
int *ip;
h = 2.8 * Softening;  
m200 = np * PartMass;
mass = PartMass;

Nrand = Nr;

particle *pdata;


if(np > Nrand)
{
	N = Nrand;
	pdata = (particle *) malloc(sizeof(particle) * N);
	ip= (int *) malloc(sizeof(int) * N);
	memset(ip, 0, sizeof(int) * N); 
	memset(pdata, 0, sizeof(particle) * N); 

	srand( (unsigned)time( NULL ) ); 
	for(i=0;i<N;i++)
        ip[i] = (int)(rand() / (double)RAND_MAX * np);

	for(k=0,i=0;i<N;i++)
	for(j=0;j<3;j++)
	pdata[i].cpos[j] = data[ip[i]].cpos[j];
	free(ip);
}
else
{
N = np;

pdata = (particle *) malloc(sizeof(particle) * N);
memset(pdata, 0, sizeof(particle) * Nrand);

for(i=0;i<N;i++)
for(j=0;j<3;j++)
pdata[i].cpos[j] = data[i].cpos[j];
}




for(i=0;i<N;i++)
for(j=0;j<N;j++) 
{
r = sqrt(pow(fof_periodic(pdata[i].cpos[0]-pdata[j].cpos[0]),2.0)+pow(fof_periodic(pdata[i].cpos[1]-pdata[j].cpos[1]),2.0)+pow(fof_periodic(pdata[i].cpos[2]-pdata[j].cpos[2]),2.0));
//   use SPH kernel
      if(r >= h)
	pot -=  mass / r;  // potential energy
      else
	{
	  u = r / h;

	  if(u < 0.5)
	    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	  else
	    wp = -3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

	pot +=  mass  / h * wp;
	}

}
free(pdata);

return G / Time * pot * ( np * ( np - 1 )/( N * ( N - 1 ) ) / 2.0 ) ;       //phy unit, weighted, twice potential, need to divide 2
}




/*
#ifdef DEBUG
//fprintf(Logfile,"In energy.c, pot=%g\n",pe);
fprintf(Logfile,"p[%d].pos=%g, %g, %g; p[%d].pos=%g, %g, %g; dr= %g, pot=%g\n",i,pdata[i].cpos[0],pdata[i].cpos[1],pdata[i].cpos[2],j,pdata[j].cpos[0],pdata[j].cpos[1],pdata[j].cpos[2],r12,ppe[i]);
fflush(Logfile);
#endif
#ifdef DEBUG
fprintf(Logfile,"In energy.c, random particle list is ");
for(i=0;i<10;i++)
fprintf(Logfile,"%d\t",ip[i]);
fprintf(Logfile,"\n");
fflush(Logfile);
#endif

*/
/*
unsigned int seedVal;
struct timeb timeBuf;
ftime(&timeBuf);
seedVal=((((unsigned int)timeBuf.time&0xFFFF)+(unsigned int)timeBuf.millitm)^(unsigned int)timeBuf.millitm);  // enhance random
srand((unsigned int)seedVal);
*/
/*
a = malloc(sizeof(int)*np);
for (i=0;i<np;i++)
  a[i]=i+1;
for (i=0;i<Nrand;i++)
{
 w=rand()%(np-i)+i;
 t=a[i];
 a[i]=a[w];
 a[w]=t;
}
for(i=0;i<Nrand;i++)
ip[i] = a[i];
free(a);
*/
/*
N=np;
particle *pdata;
pdata = (particle *) malloc(sizeof(particle) * N);
memset(pdata, 0, sizeof(particle) * N); 

for(i=0;i<N;i++)
for(j=0;j<3;j++)
pdata[i].cpos[j] = data[i].cpos[j];
*/
/*
	//picka=rand()/(RAND_MAX+1.0);
		//if(picka < ( (double)N / (double)np ))
*/

