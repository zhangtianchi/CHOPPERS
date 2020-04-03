
#include <stdio.h>
#include <malloc.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>


void endrun(int ierr);
peanokey peano_hilbert_key(int x, int y, int z, int bits);
int get_next_file(int begsnap, int begfilenr, int endsnap, int endfilenr, int *snapshot, int *filenr);
int load_hash_table(void);
void set_units(void);
int  get_spherical_region_coordinates(int grp);
double fof_periodic(double x);
double fof_periodic_wrap(double x);
int rad_sort_particle(const void *a, const void *b);
float get_fsub(int nhalo,int *nsub);
float Dvir(float om_0,float z);
void density_profile(int gr);
void fit_density_profile(int gr);
void vc_profile(int gr);
int lnNFW_f (const gsl_vector * par, void *data, gsl_vector * f);
int lnNFW_df (const gsl_vector * par, void *data, gsl_matrix * J);
int lnNFW_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);
int set_solver_data(HALOPROFILE *haloprof,gsl_multifit_function_fdf *f);
void free_solver_data(gsl_multifit_function_fdf *f);
double fit_halo_prof(HALOPROFILE *haloprof,double par[2],double err[2]);
double cnfw(int gr); 
void free_particle_data(void);
void load_sub_catalogue(void);

#ifdef UNBINDING
particle *unbinding_process( particle *HHP, int N, int *RN );
particle *Unbinding( particle *HHP, int N, int *Nnew );
#endif
float bspin( int gr, float sam[], float vel[], float cm[], float *vdisp ); 
void shape( float E[], float tensor[] ); 
int False(double a,double b,double eps,double(*f) ( ),double *x,double xx);
double pradafunc(double x,double xx);
float cvmax( int gr, float vmax, float v200 );
float kinetic( void ); 
float pspin( int gr, float ke, float pe ); 

#ifdef PROJECTION
int Index(int x, int y, int z);
double *Griding_CIC( particle *P, int N, int SIZE, double L );
float get_halo_cic(int grp);
#endif

float funNFW(int gr, float c, float r);
float bestfit(int gr, float c);

size_t tree_treeallocate(int maxnodes, int maxpart);
int tree_treebuild(void);
void tree_update_node_recursive(int no, int sib, int father);
double tree_treeevaluate_potential(int target);
void tree_treefree(void);
#ifdef OCTREE
float Potential( void );
#else
float potential( void );
#endif


















