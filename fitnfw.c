#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_eigen.h>
#include "allvars.h"

#define NUM_PARAM 2
#define TOL_FIT 1e-6
#define MAX_ITER 500
/*****************************************GSL LIB*****************************************************/
 int lnNFW_f (const gsl_vector * par, void *data, gsl_vector * f)
 { 
   size_t nbin = ((struct Hdata *)data)->nbins;
   double *r=((struct Hdata *)data)->r;
   double *Lrho = ((struct Hdata *)data)->rho;//log(rho)
   double *Lsigma = ((struct Hdata *) data)->sigma;//D(log(rho))
 
   double Lrhos = gsl_vector_get (par, 0);
   double rs = gsl_vector_get (par, 1);
 
   size_t i;
 
   for (i = 0; i < nbin; i++)
	 {
	   /* NFW  Yi = Rhos /rt/(1+rt)^2 */
	   double rt = r[i]/rs;
	   double Yi = Lrhos-log(rt)-2.*log(1+rt);   
	   gsl_vector_set (f, i, (Yi - Lrho[i])/Lsigma[i]);
	 }
 
   return GSL_SUCCESS;
 }
 

 int lnNFW_df (const gsl_vector * par, void *data, gsl_matrix * J)
 {
   size_t nbin = ((struct Hdata *)data)->nbins;
   double *r=((struct Hdata *)data)->r;
   double *Lsigma = ((struct Hdata *) data)->sigma;
 
   double rs = gsl_vector_get (par, 1);
 
   size_t i;
 
   for (i = 0; i < nbin; i++)
	 {
	   /* Jacobian matrix J(i,j) = dfi / dxj, */
	   /* where fi = (Yi - yi)/sigma[i],      */
	   /*   Yi=log(Rhoi),  Rhoi= Rhos /rt/(1+rt)^2  */
	   /* and the xj are the parameters (Lrhos,rs) */
	   double rt = r[i]/rs;
	   double s = Lsigma[i];
	   gsl_matrix_set (J, i, 0, 1./s);                      //  dlog(p)/dlog(ps)
	   gsl_matrix_set (J, i, 1, (1.+2.*rt/(1+rt))/rs/s);   //  dlog(p)/d(rs)
	 }
   return GSL_SUCCESS;
 }
 

 int lnNFW_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
 {
   lnNFW_f (x, data, f);
   lnNFW_df (x, data, J);
 
   return GSL_SUCCESS;
 }

int set_solver_data(HALOPROFILE *haloprof,gsl_multifit_function_fdf *f)
{
double *r, *y, *sigma, *vol;
int nbin,i;
struct Hdata *d;

nbin=haloprof->nbins;

/* This is the data to be fitted */
d=malloc(sizeof(struct Hdata));
r=malloc(sizeof(double)*nbin);
y=malloc(sizeof(double)*nbin);
sigma=malloc(sizeof(double)*nbin);

for (i = 0; i < nbin; i++)
 {
		   sigma[i]=1;
   r[i]=haloprof->rr[i];
   y[i]=haloprof->lgrho[i];

	   
 }
d->nbins=nbin;
d->r=r;
d->rho=y;
d->sigma=sigma; 
f->f = &lnNFW_f;
f->df = &lnNFW_df;
f->fdf = &lnNFW_fdf;
f->n = nbin;
f->p = NUM_PARAM;
f->params = d;
return 1;
}

void free_solver_data(gsl_multifit_function_fdf *f)
{
struct Hdata *d;
d=f->params;
free(d->r);
free(d->rho);
free(d->sigma);
free(d);
}
double fit_halo_prof(HALOPROFILE *haloprof,double par[2],double err[2])
{	
gsl_multifit_function_fdf f;
const gsl_multifit_fdfsolver_type *T;
gsl_multifit_fdfsolver *s;
gsl_matrix *covar;
int status;
unsigned int i,iter = 0;
double x_init[2];
x_init[0]=1.;        //init rho
x_init[1]=1.;		 //init rs
gsl_vector_view x = gsl_vector_view_array (x_init, NUM_PARAM);

status=set_solver_data(haloprof,&f);
if(status==0)
{
 par[0]=par[1]=0.;
 err[0]=err[1]=0.;
 return 1;
}
covar = gsl_matrix_alloc (NUM_PARAM, NUM_PARAM);
T = gsl_multifit_fdfsolver_lmsder;
s = gsl_multifit_fdfsolver_alloc (T, f.n, f.p);
gsl_multifit_fdfsolver_set (s, &f, &x.vector);
do
 {
   iter++;
   status = gsl_multifit_fdfsolver_iterate (s);

   if (status)
	 break;
   status = gsl_multifit_test_delta (s->dx, s->x, TOL_FIT, TOL_FIT);
 }
while (status == GSL_CONTINUE && iter < MAX_ITER);
gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
{ 
 double chi = gsl_blas_dnrm2(s->f);//2-norm of the residual function
 double dof = f.n-f.p;
 double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
 par[0]=exp(FIT(0));
 err[0]=par[0]*c*ERR(0);
 par[1]=FIT(1);
 err[1]=c*ERR(1);
//printf ("rhos = %.5f +/- %.5f\n", par[0], err[0]);
//printf ("rs   = %.5f +/- %.5f\n", par[1], err[1]);

rrhos=par[0];
rrs=par[1];
}
//printf ("iter=%u, status = %s\n\n", iter, gsl_strerror (status));
//printf("\nconcentration= %g\n",haloprof->rvir/par[1]);
free_solver_data(&f);
gsl_multifit_fdfsolver_free (s);			
gsl_matrix_free (covar);
return haloprof->rvir/par[1];
/*
if(status==GSL_SUCCESS)
return 0;
else
return 1;
*/
}

double cnfw(int gr)  
{
    int i, j, k;
    int nbin,nvir;
    int NN = 0;
    float *temp1,*temp2;
    double rvir,con;
    nbin = DENSNBIN;
    rvir = Halo_R_Crit200[gr] * Time;   // phy unit 
    nvir=(int)(Halo_M_Crit200[gr] / PartMass);

for(j = 0; j < nbin; j++)
    if(logrho[j] != 0.0)
    NN++;


temp1 = malloc( sizeof( float) * NN);
temp2 = malloc( sizeof( float) * NN);

for(k = 0, j = 0; j < nbin; j++)
{
    if(logrho[j] != 0.0)
    {    
    temp1[k] = pow(10.0, logr[j]) * Time;  // phy unit
    temp2[k] = log( pow( 10.0, logrho[j] ) / pow( Time, 3.0 ) );  //rho_phy=(1+z)^3*rho_comoving
    k++;
    }    
}

HALOPROFILE *pro;
pro = malloc( sizeof( HALOPROFILE ) );
pro->nbins = NN;
pro->rvir = rvir;
pro->rr = temp1;
pro->lgrho = temp2;

double pa[2], er[2];
if( nvir < Ncut_fitnfw )
    con=0.0;
else
con = fit_halo_prof(pro, pa, er);
if(con < 1.0 || con > 50.0 )
con = 0.0;
free(pro);
free(temp1);
free(temp2);
return con;
}

/*==============================================================================
 * small routine calculating cNFW according to Eq.(9) in Prada et al. (2012)
 *==============================================================================*/


/*use Newton-Raphson method calc C_vmax*/
int False(double a,double b,double eps,double(*f) ( ),double *x,double xx)
{
int m;
double fa,fb,y;
m=0;
fa=(*f)(a,xx);
fb=(*f)(b,xx);
if( fa*fb > 0 )
return(-1);
do
{
m=m+1;
*x=(a*fb-b*fa)/(fb-fa);
y=(*f)(*x,xx);
	if(y*fa < 0)
	{
	b=*x;
	fb=y;
	}
	else
	{
	a=*x;
	fa=y;
	}
}while (fabs(y)>=eps);
return m;
}

double pradafunc(double x,double xx)
{
double y;
y=0.216*x/(log(1+x)-x/(1+x))-xx;
return y;
}

float cvmax( int gr, float vmax, float v200 )
{
float cv;
int kconv;
double result;
double(*p) (double,double);
int nvir;
nvir=(int)(Halo_M_Crit200[gr] / PartMass);
if( nvir < Ncut_fitnfw )
    cv=0.0;
else
{
p = &pradafunc;
kconv=False( 2.1 , 50.0, 0.0000001, p, &result, pow((vmax/v200), 2.0) );   
cv=result;
}
return cv;
}



