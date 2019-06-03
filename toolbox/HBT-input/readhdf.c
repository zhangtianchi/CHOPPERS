/*****************Read HBT2_MPI HDF5 data**********************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>

#include "allvars.h"
#include "proto.h"


int scanline(char * filename)
{
        FILE * f = 0; char line[256]=""; int lines = 0;
        f = fopen(filename, "r");
        if(!f) return 0;
        while(!feof(f)) { fgets(line, 99999999, f); lines++; }
        fclose(f);
        return lines-1;

}

double fof_periodic_wrap(double x)
{
    while (x >= BBoxSize)
	x -= BBoxSize;

    while (x < 0)
	x += BBoxSize;

    return x;
}


double fof_periodic(double x)
{
    if (x >= 0.5 * BBoxSize)
	x -= BBoxSize;

    if (x < -0.5 * BBoxSize)
	x += BBoxSize;

    return x;
}


int loadsubcat(char *Dir , int snap)
{
  int i,j;
  int *tmp_int;
  float *tmp_float;
  float *tmp3_float;
  hid_t file_id, dataset_id;
  hid_t dataspace;
  hsize_t size[1];
  herr_t err;
  char buf1[1024];
  int num;
sprintf(buf1,"%s/HBT2/SubSnap_%03d.hdf5",Dir,snap);
////////////////////////////////////////////////read HDF5////////////////////////////////////////
file_id = H5Fopen(buf1, H5F_ACC_RDONLY, H5P_DEFAULT);
///////////////////////////////////get subnum///////////////////////////////////////////////////////////////
/*  get the size of dataset */
  dataset_id = H5Dopen(file_id, "/Subhalos/Nbound", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset_id);
  H5Sget_simple_extent_dims(dataspace, size, NULL);
  err = H5Sclose(dataspace);
  err = H5Dclose(dataset_id);
  num = size[0];

//printf("Nsub(include_Nbound=1):%d\n",num);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
sub= (SUBCAT *)malloc(sizeof(SUBCAT) * num);
tmp_int=malloc(sizeof(int)* num);
tmp_float=malloc(sizeof(int)* num);
tmp3_float=malloc(sizeof(float)*num*3);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/ComovingMostBoundPosition", H5P_DEFAULT);
 err = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp3_float);
  for (i = 0; i < num; i++)
    for (j = 0; j < 3; j++)
{
      sub[i].ComovingMostBoundPosition[j] = tmp3_float[i*3+j];
}
  err = H5Dclose(dataset_id);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/SpinBullock", H5P_DEFAULT);
 err = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp3_float);
  for (i = 0; i < num; i++)
    for (j = 0; j < 3; j++)
{
      sub[i].SpinBullock[j] = tmp3_float[i*3+j];
}
  err = H5Dclose(dataset_id);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/Nbound", H5P_DEFAULT);
err = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int);
  for (i = 0; i < num; i++)
{
    sub[i].Nbound = tmp_int[i];
}
  err = H5Dclose(dataset_id);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/SnapshotIndexOfBirth", H5P_DEFAULT);
err = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int);
  for (i = 0; i < num; i++)
{
    sub[i].SnapshotIndexOfBirth = tmp_int[i];
}
  err = H5Dclose(dataset_id);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/SnapshotIndexOfDeath", H5P_DEFAULT);
err = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int);
  for (i = 0; i < num; i++)
{
    sub[i].SnapshotIndexOfDeath = tmp_int[i];
}
  err = H5Dclose(dataset_id);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/Rank", H5P_DEFAULT);
err = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int);
  for (i = 0; i < num; i++)
{
    sub[i].Rank = tmp_int[i];

}
  err = H5Dclose(dataset_id);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/TrackId", H5P_DEFAULT);
err = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int);
  for (i = 0; i < num; i++)
{
    sub[i].TrackId = tmp_int[i];

}
  err = H5Dclose(dataset_id);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/HostHaloId", H5P_DEFAULT);
err = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_int);
  for (i = 0; i < num; i++)
{
    sub[i].HostHaloId = tmp_int[i];

}
  err = H5Dclose(dataset_id);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/Mbound", H5P_DEFAULT);
err = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_float);
  for (i = 0; i < num; i++)
{
    sub[i].Mbound = tmp_float[i];

}
  err = H5Dclose(dataset_id);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/M200Crit", H5P_DEFAULT);
 err = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_float);
  for (i = 0; i < num; i++)
{
      sub[i].M200Crit = tmp_float[i];
}
  err = H5Dclose(dataset_id);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataset_id = H5Dopen(file_id, "/Subhalos/R200CritComoving", H5P_DEFAULT);
 err = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp_float);
  for (i = 0; i < num; i++)
{
      sub[i].R200CritComoving = tmp_float[i];
}
  err = H5Dclose(dataset_id);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//printf("%d\n",sub[0].Nbound);


free(tmp_int);
free(tmp_float);
free(tmp3_float);

err = H5Fclose(file_id);

return num;
}






