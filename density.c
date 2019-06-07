#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "allvars.h"



void density_profile(int gr)
{
    int i;			
    int k;		
    int count;
    double rho,mvir,rvir;

    i = 1;	
    do 
      {
	i++;
	rconv = sqrt(data[i].rad);
	rho = (float) i *PartMass / (4.0 / 3.0 * M_PI * rconv * rconv * rconv);
      } while ((float) i / (8.0 * log((float) i)) * sqrt(200.0 * RHO_CRIT / (rho / pow( Time, 3.0 ) ) ) < KAPPA);  // rho_c/rho_phy


    mvir = Halo_M_Crit200[gr];
    rvir = Halo_R_Crit200[gr];
    
    double logr0 = log10(Rpromin)+log10(rvir);	
    double dlogr = (log10(Rpromax)+log10(rvir)-logr0)/DENSNBIN;   //  output profile range [Rpromin*R200,Rpromax*R200]
    
    k = 0;  /* step over bins */


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
	logr[k] = (logr0 + dlogr * (k + 0.5));  //  r for fit nfw
	denr[k] = pow(10., logr0 + dlogr * (k + 0.5)) / Halo_R_Crit200[gr];  //  r / r200 
	if (count>0)
          {
	  logrho[k] = log10(count * PartMass / ((4.0 / 3.0) * M_PI * ( pow(pow(10.0, (logr0 + dlogr * (k+1))), 3.0) - pow(pow(10.0, (logr0 + dlogr * k )), 3.0))));
      denrho[k] = pow(10., logrho[k]) / pow(Time, 3.0) / RHO_CRIT ;     // rho / rhoc
          }
	else
          {
	  logrho[k]=0.0;
	  denrho[k]=0.0;
          }
	k++;
    }
}
