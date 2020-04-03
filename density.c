#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "allvars.h"



void density_profile(int gr)
{
    int i;			
    int j,n;			
    double r;
    int k;		
    int count;
    double rho,mvir,rvir;

    i = 1;	
    do 
      {
	i++;
	rconv = sqrt(data[i].rad);
	rho = (float) i *PartMass / (4.0 / 3.0 * M_PI * rconv * rconv * rconv);
#ifdef DeltaVir
      } while (((float) i / (8.0 * log((float) i)) * sqrt(DELTA * RHO_CRIT / rho / pow(Time,3.0)  )) < KAPPA);    // rho_c/rho_phy
#else
      } while (((float) i / (8.0 * log((float) i)) * sqrt(200.0 * RHO_CRIT / rho / pow(Time,3.0)  )) < KAPPA);    // rho_c/rho_phy
#endif

    mvir = Halo_M_Crit200[gr];
    rvir = Halo_R_Crit200[gr];
    
#ifdef MYWORK
    //double logr0 = log10(BoxSize / pow(TotNumPart,1.0/3.0) / 50.0 / rvir )+log10(rvir); //  2.8softening	
    //double logr0 = log10(0.0039); //  softening	
    //double dlogr = (log10(1.0)+log10(rvir)-logr0)/DENSNBIN;    //  [softening, r200]
    double logr0 = log10(Rpromin)+log10(rvir);	
    double dlogr = (log10(Rpromax)+log10(rvir)-logr0)/DENSNBIN;   //  output profile range [Rpromin*R200,Rpromax*R200]
#else
    double logr0 = log10(Rpromin)+log10(rvir);	
    double dlogr = (log10(Rpromax)+log10(rvir)-logr0)/DENSNBIN;   //  output profile range [Rpromin*R200,Rpromax*R200]
#endif
/*    
    k = 0; 

#ifdef DEBUG
fprintf(Logfile,"In density.c halo %d, np=%d,data[%d].rad=%g,sqrt data[%d].rad=%g\n", gr, np,np-1,data[np-1].rad,np-1,sqrt(data[np-1].rad));
fflush(Logfile);
#endif
*/

for(j=0;j<DENSNBIN;j++)  
{
n=0;
for(i=0;i<np;i++)   
{
if((log10(sqrt(data[i].rad))>=logr0+j*dlogr)&&(log10(sqrt(data[i].rad))<logr0+(j+1)*dlogr))
n++;
}
npbin[j]=n;
r=((logr0+j*dlogr)+(logr0+(j+1)*dlogr))/2.0;
denr[j]=pow(10.0,r)/Halo_R_Crit200[gr];
if(n==0)
denrho[j]=0.0;
denrho[j]=(PartMass*n)/(4.0/3.0*M_PI*(pow(pow(10,logr0+(j+1)*dlogr),3.0)-pow(pow(10,logr0+j*dlogr),3.0)))/RHO_CRIT;
//#ifdef DEBUG
//fprintf(Logfile,"In density.c npbin[%d]=%d,denr[%d]=%g,denrho[%d]=%g\n", j, n,j,denr[j],j,denrho[j]);
//fflush(Logfile);
//#endif
}







/*
    while (k < DENSNBIN) 
      {
    	i = 0;
	count = 0;
	while ( (0.5*log10(data[i].rad) < (logr0+dlogr*k)) && (i<np) ) i++;  // count spherical region particle and find boundry
	while ( (0.5*log10(data[i].rad) < (logr0+dlogr*(k+1))) && (i<np) ) 
	  {
	    i++;
	    count++;
	  }
	npbin[k] = count;

	denr[k] = pow(10., logr0 + dlogr * (k + 0.5)) / Halo_R_Crit200[gr];  //  r / r200 
	if (count>0)
          {
	  logrho[k] = log10(count * PartMass / ((4.0 / 3.0) * M_PI * ( pow(pow(10.0, (logr0 + dlogr * (k+1))), 3.0) - pow(pow(10.0, (logr0 + dlogr * k )), 3.0))));
      denrho[k] = pow(10., logrho[k])  / RHO_CRIT ;     // rho / rhoc
          }
	else
          {
	  denrho[k]=0.0;
          }
#ifdef DEBUG
fprintf(Logfile,"In density.c npbin[%d]=%d,denr[%d]=%g,denrho[%d]=%g\n", k, count,k,denr[k],k,denrho[k]);
fflush(Logfile);
#endif
	k++;
    }
*/
}



void fit_density_profile(int gr)
{
    int i;			
    int k;		
    int j,n;
    double r;
    int count;
    double rho,mvir,rvir;

    mvir = Halo_M_Crit200[gr];
    rvir = Halo_R_Crit200[gr];
    
#ifdef MYWORK
    //double logr0 = log10(BoxSize / pow(TotNumPart,1.0/3.0) / 50.0 / rvir )+log10(rvir); //  2.8softening	
    //double logr0 = log10(0.0039); //  softening	
    //double dlogr = (log10(1.0)+log10(rvir)-logr0)/DENSNBIN;    //  [softening, r200]
   double logr0 = log10(Rpromin)+log10(rvir);	
   double dlogr = (log10(Rpromax)+log10(rvir)-logr0)/DENSNBIN;   //  output profile range [Rpromin*R200,Rpromax*R200]
#else
    double logr0 = log10(Rpromin)+log10(rvir);	
    double dlogr = (log10(Rpromax)+log10(rvir)-logr0)/DENSNBIN;   //  output profile range [Rpromin*R200,Rpromax*R200]
#endif


for(j=0;j<DENSNBIN;j++)  
{
n=0;
for(i=0;i<np;i++)   
{
if((log10(sqrt(data[i].rad))>=logr0+j*dlogr)&&(log10(sqrt(data[i].rad))<logr0+(j+1)*dlogr))
n++;
}
r=((logr0+j*dlogr)+(logr0+(j+1)*dlogr))/2.0;
logr[j]=r;
if(n==0)
logrho[j]=0.0;
logrho[j]=log10((PartMass*n)/(4.0/3.0*M_PI*(pow(pow(10,logr0+(j+1)*dlogr),3.0)-pow(pow(10,logr0+j*dlogr),3.0))));
//#ifdef DEBUG
//fprintf(Logfile,"In fitnfw.c npbin[%d]=%d,denr[%d]=%g,denrho[%d]=%g\n", j, n,j,denr[j],j,denrho[j]);
//fflush(Logfile);
//#endif
}

/*
    k = 0;  


    while (k < DENSNBIN) 
      {
    	i = 0;
	    count = 0;
	while ( (0.5*log10(data[i].rad) < (logr0+dlogr*k)) && (i<np) ) i++;  // count spherical region particle and find boundry
	while ( (0.5*log10(data[i].rad) < (logr0+dlogr*(k+1))) && (i<np) ) 
	  {
	    i++;
	    count++;
	  }

	logr[k] = logr0 + dlogr * (k + 0.5);  // log r for fit NFW

	if (count>0)
          {
	  logrho[k] = log10(count * PartMass / ((4.0 / 3.0) * M_PI * ( pow(pow(10.0, (logr0 + dlogr * (k+1))), 3.0) - pow(pow(10.0, (logr0 + dlogr * k )), 3.0))));    //  log rho for fit NFW 
          }
	else
          {
	  logrho[k]=0.0;
          }
	k++;
    }
*/
}
/*
#ifdef DEBUG
fprintf(Logfile,"In density.c halo %d, in %d bin, i=%d, count=%d\n", gr, k, i, count);
fflush(Logfile);
#endif
*/
