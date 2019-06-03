#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "allvars.h"



void vc_profile(int gr)
{
    int i;			
    int k;		
    int count;
    double rho,mvir,rvir;

    rvir = Halo_R_Crit200[gr];
    
    double logr0 =  log10(Rpromin)+log10(rvir);	
    double dlogr = (log10(Rpromax)+log10(rvir)-logr0)/VELNBIN;   //  output profile range [Rpromin*R200,Rpromax*R200]
    
    k = 0;  /* step over bins */


    while (k < VELNBIN) 
      {
    	i = 0;
	    count = 0;
	while ( (0.5*log10(data[i].rad) < (logr0+dlogr*k)) && (i<np) ) 
	  {
	    i++;
	    count++;
	  }
	velnp[k] = count;
	velr[k] = pow(10., logr0 + dlogr * k) / Halo_R_Crit200[gr];  //  r / r200 
	if (count>0)
      velv[k] = sqrt( G * PartMass * count / pow( 10., logr0 + dlogr * k ) );     // Vc(r)
	else
	  velv[k]=0.0;
	k++;
    }
}
