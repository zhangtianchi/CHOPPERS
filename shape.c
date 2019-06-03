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

#define SQR(a) (a*a)

/* use moment of inertia calc a(sqrt(lambda1,max)), b, c */
void shape( float E[], float tensor[] )   //   return a, b, c and inertia tensor  
{

  double inertia[6];
  int i,j;
  double eval_i[3];
  double shape[9];
  double a,b,c;

    for(j = 0; j < 6; j++)
      inertia[j] = 0.0;

for(i = 0; i < np; i++)
{
inertia[0]+=SQR(data[i].pos[0]);  //Ixx
inertia[3]+=SQR(data[i].pos[1]);  //Iyy
inertia[5]+=SQR(data[i].pos[2]);  //Izz
inertia[1]+=data[i].pos[0]*data[i].pos[1]; //Ixy
inertia[2]+=data[i].pos[0]*data[i].pos[2]; //Ixz
inertia[4]+=data[i].pos[1]*data[i].pos[2]; //Iyz
}
for(i = 0; i < 6; i++)
  tensor[i] = inertia[i];

//We store only 9 components since inertia tensor is symetric 
  shape[0] = inertia[0];
  shape[1] = inertia[1];
  shape[2] = inertia[2];
  shape[3] = inertia[1];
  shape[4] = inertia[3];
  shape[5] = inertia[4];
  shape[6] = inertia[2];
  shape[7] = inertia[4];
  shape[8] = inertia[5];

    gsl_matrix_view m  = gsl_matrix_view_array (shape, 3, 3);
    gsl_vector *eval = gsl_vector_alloc (3); 
    gsl_matrix *evec = gsl_matrix_alloc (3, 3);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
    gsl_eigen_symmv (&m.matrix, eval, evec, w); 
    gsl_eigen_symmv_free (w); 
    gsl_eigen_symmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_DESC); // ABS_ASC is little->big
          for (i = 0; i < 3; i++)
        {
            eval_i[i] = gsl_vector_get (eval, i);
            eval_i[i] = sqrt(fabs(eval_i[i]));
           // gsl_vector_view evec_i = gsl_matrix_column (evec, i);
           // printf ("eigenvalue = %g\n", eval_i); 
           // printf ("eigenvector = \n"); 
          //  gsl_vector_fprintf (stdout,&evec_i.vector, "%g");
                                                                                                
        }
            gsl_vector_free (eval); 
            gsl_matrix_free (evec);

for(i=0;i<3;i++)
E[i] = eval_i[i];


}

