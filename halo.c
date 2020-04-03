#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);


void free_particle_data(void)
{

  /* free stuff from load_sub_catalogue() */

  free(NsubPerHalo);
  free(FirstSubOfHalo);
  free(SubLen);
  free(SubOffset);
  free(SubParentHalo);
  free(Halo_M_Mean200);
  free(Halo_R_Mean200);
  free(Halo_M_Crit200);
  free(Halo_R_Crit200);
  free(Halo_M_TopHat200);
  free(Halo_R_TopHat200);
  free(SubPos);
  free(SubVel);
  free(SubVelDisp);
  free(SubVmax);
  free(SubSpin);
  free(SubMostBoundID);
  free(SubHalfMass);

}

int load_hash_table(void)
{
  int i, dummy, filenr, firstcell, lastcell;
  char buf[512];
  FILE *fd;

  sprintf(buf, "%s/snapdir_%03d/%s_%03d.0", OutputDir, SnapshotNum, SnapshotFileBase, SnapshotNum);

  if(!(fd = fopen(buf, "r")))
    {
      printf("can't open file `%s'\n", buf);
      endrun(1);
    }

  my_fread(&dummy, sizeof(int), 1, fd);
  my_fread(&header, sizeof(header), 1, fd);
  my_fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);

  NFiles = header.num_files;
  HashTabSize = header.hashtabsize;
  BoxSize = header.BoxSize;
  Omega0 = header.Omega0;
  OmegaLambda = header.OmegaLambda;
  HubbleParam = header.HubbleParam;
  Time = header.time;

  TotNumPart = header.npartTotal[1] + (((long long) header.npartTotal[2]) << 32);
  PartMass = header.mass[1];
#ifdef DEBUG
  fprintf(Logfile, "TotNumPart=%d%09d\n", (int) (TotNumPart / 1000000000), (int) (TotNumPart % 1000000000));
  fprintf(Logfile, "hashtabsize=%d\n", header.hashtabsize );
  fflush(Logfile);
#endif
  HashTable = malloc(HashTabSize * sizeof(int));
  FileTable = malloc(HashTabSize * sizeof(int));
  LastHashCell = malloc(NFiles * sizeof(int));
  NInFiles = malloc(NFiles * sizeof(int));



  for(filenr = 0; filenr < NFiles; filenr++)
    {
      sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", OutputDir, SnapshotNum, SnapshotFileBase, SnapshotNum,
	      filenr);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  endrun(1);
	}
      my_fread(&dummy, sizeof(int), 1, fd);
      my_fread(&header, sizeof(header), 1, fd);
      my_fread(&dummy, sizeof(int), 1, fd);

      /* skip positions */
      my_fread(&dummy, sizeof(int), 1, fd);
      fseek(fd, dummy, SEEK_CUR);
      my_fread(&dummy, sizeof(int), 1, fd);

      /* skip velocities */
      my_fread(&dummy, sizeof(int), 1, fd);
      fseek(fd, dummy, SEEK_CUR);
      my_fread(&dummy, sizeof(int), 1, fd);

      /* skip IDs */
      my_fread(&dummy, sizeof(int), 1, fd);
      fseek(fd, dummy, SEEK_CUR);
      my_fread(&dummy, sizeof(int), 1, fd);

      my_fread(&dummy, sizeof(int), 1, fd);
      my_fread(&firstcell, sizeof(int), 1, fd);
      my_fread(&lastcell, sizeof(int), 1, fd);
      my_fread(&dummy, sizeof(int), 1, fd);

      LastHashCell[filenr] = lastcell;
      NInFiles[filenr] = header.npart[1];

      if(firstcell < 0 || firstcell >= header.hashtabsize || lastcell < 0 || lastcell >= header.hashtabsize)
	{
	  printf("Error in %d!\n",ThisTask);
	  endrun(1);
	}

      my_fread(&dummy, sizeof(int), 1, fd);
      my_fread(&HashTable[firstcell], sizeof(int), lastcell - firstcell + 1, fd);
      my_fread(&dummy, sizeof(int), 1, fd);

      for(i = 0; i < lastcell - firstcell + 1; i++)
	FileTable[firstcell + i] = filenr;

      fclose(fd);
    }

  return 0;
}

void load_sub_catalogue(void)
{

    int i,tmp;
    int filenr = 0;
    int TotNsubs = 0;
    int Nfiles = 512;
    char buf[1024];
    FILE *fd;

    /* read the total number of groups and subhalos */

    sprintf(buf, "%s/postproc_%03d/sub_tab_%03d.%d", OutputDir, SnapshotNum, SnapshotNum, FileNr);

    if (!(fd = fopen(buf, "r"))) {
      printf("can't open file %s\n", buf);
      endrun(1);
    }
    
    my_fread(&Ngroups, sizeof(int), 1, fd);
    my_fread(&Nids, sizeof(int), 1, fd);
    my_fread(&TotNgroups, sizeof(int), 1, fd);
    my_fread(&Nfiles, sizeof(int), 1, fd);
    my_fread(&Nsubhalos, sizeof(int), 1, fd);
    printf("In Task %d, Snap=%d, rhoc=%g,rhom=%g,omz=%g, Dvir=%g, Dvir=%g rho_mean, z=%g, TotNgroups=%d, Nhalo=%d\n", ThisTask, SnapshotNum, RHO_CRIT, RHO_MEAN, OMZ, DELTA, DELTA/OMZ, 1.0/Time-1.0, TotNgroups, Ngroups);
    fflush(stdout);

    NsubPerHalo = malloc(sizeof(int)*Ngroups);   // HBT2 rank  rank=0 halo, rank!=0 subhalo
    my_fread(NsubPerHalo, sizeof(int), Ngroups, fd);
    FirstSubOfHalo = malloc(sizeof(int)*Ngroups);
    my_fread(FirstSubOfHalo, sizeof(int), Ngroups, fd);

    SubLen=malloc(sizeof(int)*Nsubhalos);
    my_fread(SubLen, sizeof(int), Nsubhalos, fd);
    SubOffset=malloc(sizeof(int)*Nsubhalos);
    my_fread(SubOffset, sizeof(int), Nsubhalos, fd);

    SubParentHalo=malloc(sizeof(int)*Nsubhalos);    
    my_fread(SubParentHalo, sizeof(int), Nsubhalos, fd);

    Halo_M_Mean200=malloc(sizeof(float)*Ngroups);    //  HBT2 replaced by Dfifth
    my_fread(Halo_M_Mean200, sizeof(float), Ngroups, fd);


    Halo_R_Mean200=malloc(sizeof(float)*Ngroups);
    my_fread(Halo_R_Mean200, sizeof(float), Ngroups, fd);


    Halo_M_Crit200=malloc(sizeof(float)*Ngroups);    
    my_fread(Halo_M_Crit200, sizeof(float), Ngroups, fd);

    Halo_R_Crit200=malloc(sizeof(float)*Ngroups);
    my_fread(Halo_R_Crit200, sizeof(float), Ngroups, fd);

    Halo_M_TopHat200=malloc(sizeof(float)*Ngroups);
    my_fread(Halo_M_TopHat200, sizeof(float), Ngroups, fd);

    Halo_R_TopHat200=malloc(sizeof(float)*Ngroups);
    my_fread(Halo_R_TopHat200, sizeof(float), Ngroups, fd);

    SubPos = malloc(Nsubhalos * 3 * sizeof(float));
    my_fread(SubPos, 3 * sizeof(float), Nsubhalos, fd);

    SubVel = malloc(Nsubhalos * 3 * sizeof(float));
    my_fread(SubVel, 3 * sizeof(float), Nsubhalos, fd);

    SubVelDisp = malloc(Nsubhalos * sizeof(float));
    my_fread(SubVelDisp, sizeof(float), Nsubhalos, fd);

    SubVmax = malloc(Nsubhalos * sizeof(float));
    my_fread(SubVmax, sizeof(float), Nsubhalos, fd);

    SubSpin = malloc(Nsubhalos * 3 * sizeof(float));
    my_fread(SubSpin, 3 * sizeof(float), Nsubhalos, fd);

    SubMostBoundID = malloc(Nsubhalos * sizeof(long long));
    my_fread(SubMostBoundID, sizeof(long long), Nsubhalos, fd);

    SubHalfMass = malloc(Nsubhalos * sizeof(float));
    my_fread(SubHalfMass, sizeof(float), Nsubhalos, fd);

    fclose(fd);
}


/*! This catches I/O errors occuring for my_fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fread) on task=%d has occured.\n", ThisTask);
      fflush(stdout);
      endrun(778);
    }
  return nread;
}


int get_spherical_region_coordinates(int grp)
{
    int dummy, i, j, k, count,nvir;
    char buf[1000];
    FILE *fd;
    float rvir,mvir,vv;
    struct io_header header;

    int len,isub;
    int rank=0;
    float pos[3];

    float *ppos; /* particle positions */
    float *vvel; /* particle velocities */
#ifdef HALOID
    long long *iid;
#endif
    double sx,sy,sz;

    sx=sy=sz=0.0;

    particle *tmp;

    double cellsize;

    int cells, cx, cy, cz, ix, iy, iz, iix, iiy, iiz, base, hashbits, hashkey, hashtabsize;
    
    float refx, refy, refz;
    

    sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", OutputDir, SnapshotNum, SnapshotFileBase,SnapshotNum, 0);


    rvir=Halo_R_Crit200[grp];
    mvir=Halo_M_Crit200[grp];
    isub=FirstSubOfHalo[grp];
    rank=NsubPerHalo[grp];
    pos[0]=SubPos[3*isub];
    pos[1]=SubPos[3*isub+1];
    pos[2]=SubPos[3*isub+2];
    
    
    if (!(fd = fopen(buf, "r"))) {
	printf("can't open file %s\n", buf);
	endrun(1);
    }
    
    fread(&dummy, sizeof(int), 1, fd);
    fread(&header, sizeof(struct io_header), 1, fd);
    fread(&dummy, sizeof(int), 1, fd);
    fclose(fd);

    BoxSize = header.BoxSize;

    hashtabsize = header.hashtabsize;



    base = 1;
    while (base * base * base < hashtabsize)
	base *= 2;

    hashbits = 1;
    while ((1 << hashbits) < base)
	hashbits++;

    cellsize = BoxSize / (double) base;

#ifdef MYWORK
    cells = (int) ceil( 2.0 * rvir / cellsize + 4.0);
#else
        cells = (int) ( 3.0 * rvir / cellsize + 1.0);    //    halo range cell
      //cells = (int) ( 6.0*rvir / cellsize + 2.0);    //    halo range cell
      //cells = (int) ( 6.0*rvir / cellsize + 1.0);    //    halo range cell
#endif

#ifdef DEBUG
fprintf(Logfile,"In halo.c, pos = %g, %g, %g\n", pos[0], pos[1], pos[2]);
fprintf(Logfile,"In halo.c, cellsize=%g, cells=%d\n",cellsize, cells);
fprintf(Logfile,"In halo.c, Base=%d, hashbits=%d\n",base,hashbits);
fflush(Logfile);
#endif

    cx = (int) (pos[0] / cellsize);
    cy = (int) (pos[1] / cellsize);
    cz = (int) (pos[2] / cellsize);


    count = 0;

    for (ix = cx - cells; ix <= cx + cells; ix++)
	for (iy = cy - cells; iy <= cy + cells; iy++)
	    for (iz = cz - cells; iz <= cz + cells; iz++) {
		iix = ix;
		iiy = iy;
		iiz = iz;

		if (iix < 0)
		    iix += base;
		if (iix >= base)
		    iix -= base;

		if (iiy < 0)
		    iiy += base;
		if (iiy >= base)
		    iiy -= base;

		if (iiz < 0)
		    iiz += base;
		if (iiz >= base)
		    iiz -= base;

		hashkey = peano_hilbert_key(iix, iiy, iiz, hashbits);


		if (hashkey != LastHashCell[FileTable[hashkey]])
		    count += HashTable[hashkey + 1] - HashTable[hashkey];
		else
		    count += NInFiles[FileTable[hashkey]] - HashTable[hashkey];
	    }

#ifdef DEBUG
fprintf(Logfile,"In halo.c, allcount=%d\n",count);
fflush(Logfile);
#endif

    /* allocate space for positions and velocities */

    ppos = (float *) malloc(3 * sizeof(float) * count);
    vvel = (float *) malloc(3 * sizeof(float) * count);
#ifdef HALOID
    iid = (long long *) malloc( sizeof(long long) * count);
#endif

/*positions*/ 
    count = 0;
    for (ix = cx - cells; ix <= cx + cells; ix++)
	for (iy = cy - cells; iy <= cy + cells; iy++)
	    for (iz = cz - cells; iz <= cz + cells; iz++) {
		iix = ix;
		iiy = iy;
		iiz = iz;

		if (iix < 0)
		    iix += base;
		if (iix >= base)
		    iix -= base;

		if (iiy < 0)
		    iiy += base;
		if (iiy >= base)
		    iiy -= base;

		if (iiz < 0)
		    iiz += base;
		if (iiz >= base)
		    iiz -= base;

		hashkey = peano_hilbert_key(iix, iiy, iiz, hashbits);

		sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", OutputDir, SnapshotNum, SnapshotFileBase,SnapshotNum, FileTable[hashkey]);

		if (!(fd = fopen(buf, "r"))) {
		    printf("can't open file %s\n", buf);
		    endrun(1);
		}
		fread(&dummy, sizeof(int), 1, fd);
		fread(&header, sizeof(struct io_header), 1, fd);
		fread(&dummy, sizeof(int), 1, fd);

		BoxSize = header.BoxSize;

		fread(&dummy, sizeof(int), 1, fd);
		fseek(fd, 3 * HashTable[hashkey] * sizeof(float), SEEK_CUR);


		if (hashkey != LastHashCell[FileTable[hashkey]])
		    len = HashTable[hashkey + 1] - HashTable[hashkey];
		else
		    len = NInFiles[FileTable[hashkey]] - HashTable[hashkey];

		fread(&ppos[3 * count], sizeof(float), 3 * len, fd);

		count += len;
		fclose(fd);
	    }

/*Velocities*/
	count = 0;
    for (ix = cx - cells; ix <= cx + cells; ix++)
	for (iy = cy - cells; iy <= cy + cells; iy++)
	    for (iz = cz - cells; iz <= cz + cells; iz++) {
		iix = ix;
		iiy = iy;
		iiz = iz;

		if (iix < 0)
		    iix += base;
		if (iix >= base)
		    iix -= base;

		if (iiy < 0)
		    iiy += base;
		if (iiy >= base)
		    iiy -= base;

		if (iiz < 0)
		    iiz += base;
		if (iiz >= base)
		    iiz -= base;

		hashkey = peano_hilbert_key(iix, iiy, iiz, hashbits);

		sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", OutputDir, SnapshotNum,SnapshotFileBase,SnapshotNum, FileTable[hashkey]);

		if (!(fd = fopen(buf, "r"))) {
		    printf("can't open file `%s'\n", buf);
		    exit(0);
		}
		fread(&dummy, sizeof(int), 1, fd);
		fread(&header, sizeof(struct io_header), 1, fd);
		fread(&dummy, sizeof(int), 1, fd);


      fread(&dummy, sizeof(int), 1, fd);
      fseek(fd, dummy, SEEK_CUR);
      fread(&dummy, sizeof(int), 1, fd);
		fread(&dummy, sizeof(int), 1, fd);
		fseek(fd, 3 * HashTable[hashkey] * sizeof(float), SEEK_CUR);


		if (hashkey != LastHashCell[FileTable[hashkey]])
		    len = HashTable[hashkey + 1] - HashTable[hashkey];
		else
		    len = NInFiles[FileTable[hashkey]] - HashTable[hashkey];

		fread(&vvel[3 * count], sizeof(float), 3 * len, fd);
		count += len;
		fclose(fd);
	    }
#ifdef HALOID
/*id*/
	count = 0;
    for (ix = cx - cells; ix <= cx + cells; ix++)
	for (iy = cy - cells; iy <= cy + cells; iy++)
	    for (iz = cz - cells; iz <= cz + cells; iz++) {
		iix = ix;
		iiy = iy;
		iiz = iz;

		if (iix < 0)
		    iix += base;
		if (iix >= base)
		    iix -= base;

		if (iiy < 0)
		    iiy += base;
		if (iiy >= base)
		    iiy -= base;

		if (iiz < 0)
		    iiz += base;
		if (iiz >= base)
		    iiz -= base;

		hashkey = peano_hilbert_key(iix, iiy, iiz, hashbits);

		sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", OutputDir, SnapshotNum,SnapshotFileBase,SnapshotNum, FileTable[hashkey]);


		if (!(fd = fopen(buf, "r"))) {
		    printf("can't open file `%s'\n", buf);
		    exit(0);
		}
		fread(&dummy, sizeof(int), 1, fd);
		fread(&header, sizeof(struct io_header), 1, fd);
		fread(&dummy, sizeof(int), 1, fd);


      /* skip positions and velocities */
      fread(&dummy, sizeof(int), 1, fd);
      fseek(fd, dummy, SEEK_CUR);
      fread(&dummy, sizeof(int), 1, fd);
      fread(&dummy, sizeof(int), 1, fd);
      fseek(fd, dummy, SEEK_CUR);
      fread(&dummy, sizeof(int), 1, fd);
      /* read id */
	  fread(&dummy, sizeof(int), 1, fd);
	  fseek(fd, HashTable[hashkey] * sizeof(long long), SEEK_CUR);


		if (hashkey != LastHashCell[FileTable[hashkey]])
		    len = HashTable[hashkey + 1] - HashTable[hashkey];
		else
		    len = NInFiles[FileTable[hashkey]] - HashTable[hashkey];

		fread(&iid[count], sizeof(long long), len, fd);
		count += len;
		fclose(fd);
	    }
//////////////////////////////////////////////////////
#endif

    data = (particle *) malloc(sizeof(particle) * count);
    memset(data, 0, sizeof(particle) * count);  // initialization

    for (i = 0; i < count; i++) 
    {
        for (j = 0; j < 3; j++) 
        {
	        data[i].cpos[j] = ppos[i * 3 + j];
	        data[i].cvel[j] = vvel[i * 3 + j];
	        data[i].vel[j]  = vvel[i * 3 + j];
	        data[i].pos[j]  = fof_periodic(ppos[i * 3 + j] - pos[j]);
	        //data[i].cpos[j]  = fof_periodic(ppos[i * 3 + j] - pos[j]);
        }
#ifdef HALOID
            data[i].id = iid[i];
#endif
	data[i].rad = data[i].pos[0] * data[i].pos[0] + data[i].pos[1] * data[i].pos[1] + data[i].pos[2] * data[i].pos[2];
    }
    free(ppos);
    free(vvel);
#ifdef HALOID
    free(iid);
#endif
    /* sort by distance from the center */

    qsort(data, count, sizeof(particle), rad_sort_particle);

#ifdef DEBUG
fprintf(Logfile,"In halo %d, Orgin trackid=%d, M200=%g, r200=%g, N200=%d, rank=%d\n",grp, SubLen[isub], Halo_M_Crit200[grp], Halo_R_Crit200[grp], (int)(Halo_M_Crit200[grp]/PartMass),NsubPerHalo[grp]);
fflush(Logfile);
#endif
#ifdef USEHBTR200
if(rank==0)
{

    double tmm=0.0;    
    for( i=0; i<count; i++ )
    {
#ifdef DeltaVir
        tmm=pow( PartMass * (i+1) / ( DELTA * RHO_CRIT * 4.0 * 3.1415926 /3.0   ) , 1.0 / 3.0  ); // Comoving rvir   vir*rho_crit = virr * rho_mean  at z=0,vir=101,virr=340
        //tmm=pow( PartMass * (i+1) / ( DELTA / OMZ * RHO_MEAN * 4.0 * 3.1415926 /3.0   ) , 1.0 / 3.0  ); // Comoving rvir   vir*rho_crit = virr * rho_mean  at z=0, om=0.3  => vir=101,virr=340
#else
        tmm=pow( PartMass * (i+1) / ( 200.0 * RHO_CRIT * 4.0 * 3.1415926 /3.0   ) , 1.0 / 3.0  ); // Comoving r200
#endif
        if(sqrt(data[i].rad) - tmm > 0)  // tmm * Time is physical unit and Time * sqrt(data[i].rad) is also physical unit, so it is no problem!
        {
        Halo_R_Crit200[grp] = sqrt(data[i].rad);   // Comoving r200
        Halo_M_Crit200[grp] = (i + 1) * PartMass;
        break;
        }
    }
}
else
{
        Halo_R_Crit200[grp] = rvir;   // Comoving r200
#ifdef DeltaVir
        Halo_M_Crit200[grp] = DELTA*RHO_CRIT*4.0/3.0*3.1415926*pow(rvir,3.0);
#else
        Halo_M_Crit200[grp] = 200.0*RHO_CRIT*4.0/3.0*3.1415926*pow(rvir,3.0);
#endif
}
#else
    double tmm=0.0;    
    for( i=0; i<count; i++ )
    {
#ifdef DeltaVir
        tmm=pow( PartMass * (i+1) / ( DELTA * RHO_CRIT * 4.0 * 3.1415926 /3.0   ) , 1.0 / 3.0  ); // Comoving rvir   vir*rho_crit = virr * rho_mean  at z=0,vir=101,virr=340
        //tmm=pow( PartMass * (i+1) / ( DELTA / OMZ * RHO_MEAN * 4.0 * 3.1415926 /3.0   ) , 1.0 / 3.0  ); // Comoving rvir   vir*rho_crit = virr * rho_mean  at z=0, om=0.3  => vir=101,virr=340
#else
        tmm=pow( PartMass * (i+1) / ( 200.0 * RHO_CRIT * 4.0 * 3.1415926 /3.0   ) , 1.0 / 3.0  ); // Comoving r200
#endif
        if(sqrt(data[i].rad) - tmm > 0)  // tmm * Time is physical unit and Time * sqrt(data[i].rad) is also physical unit, so it is no problem!
        {
        Halo_R_Crit200[grp] = sqrt(data[i].rad);   // Comoving r200
        Halo_M_Crit200[grp] = (i + 1) * PartMass;
        break;
        }
    }
#endif
#ifdef DEBUG
fprintf(Logfile,"In halo %d, after calcr200, M200=%g, r200=%g, r200phy=%g, N200=%d\n",grp ,Halo_M_Crit200[grp], Halo_R_Crit200[grp], Halo_R_Crit200[grp]*Time, (int)(Halo_M_Crit200[grp]/PartMass));
fflush(Logfile);
#endif

    mvir = Halo_M_Crit200[grp];
    rvir = Halo_R_Crit200[grp];
    nvir=(int)(mvir/PartMass);

    count = 0;

#ifdef MYWORK
    //while (data[count].rad <= 4.0 * rvir * rvir) 
    while (data[count].rad <= rvir * rvir) 
	count++;
#else
    while (data[count].rad <= rvir * rvir) 
	count++;
#endif

#ifdef DEBUG
fprintf(Logfile,"In halo %d, after calcr200, count N200=%d\n",grp, count);
fflush(Logfile);
#endif

    tmp = (particle *) malloc(sizeof(particle) * count);
    memset(tmp, 0, sizeof(particle) * count);  // initialization
    for (i = 0; i < count; i++) 
      {
      for (j = 0; j < 3; j++) 
          {
	      tmp[i].pos[j] = data[i].pos[j];
	      tmp[i].cpos[j] = data[i].cpos[j];
	      tmp[i].vel[j] = data[i].vel[j];
	      tmp[i].cvel[j] = data[i].cvel[j];
          }
#ifdef HALOID
          tmp[i].id = data[i].id;
#endif
	      tmp[i].rad = data[i].rad;
      }

    free(data);
    data = (particle *) malloc(sizeof(particle) * count);
    memset(data, 0, sizeof(particle) * count);  // initialization

    for (i = 0; i < count; i++) 
      {
      for (j = 0; j < 3; j++) 
          {
	      data[i].pos[j]  = tmp[i].pos[j];
	      data[i].cpos[j] = tmp[i].cpos[j];
	      data[i].vel[j]  = tmp[i].vel[j];
	      data[i].cvel[j] = tmp[i].cvel[j];
          }
#ifdef HALOID
          data[i].id = tmp[i].id;
#endif
	      data[i].rad = tmp[i].rad;
      }
    qsort(data, count, sizeof(particle), rad_sort_particle);

    free(tmp);
#ifdef UNBINDING


int newn200=0;
tmp = (particle *) malloc(sizeof(particle) * count);
memset(tmp, 0, sizeof(particle) * count);  // initialization
tmp = Unbinding( data, count, &newn200 );  // N = halo_Np, newn200 = n200 after Unbinding
    Halo_M_Crit200[grp] = newn200 * PartMass;
    free(data);
    data = (particle *) malloc(sizeof(particle) * newn200);
    memset(data, 0, sizeof(particle) * newn200);  // initialization
    for (i = 0; i < newn200; i++) 
      {
      for (j = 0; j < 3; j++) 
          {
	      data[i].pos[j]  = fof_periodic( fof_periodic_wrap( tmp[i].cpos[j] + pos[j] ) - pos[j]);
	      data[i].cpos[j] = fof_periodic_wrap( tmp[i].cpos[j] + pos[j] );
	      data[i].vel[j]  = tmp[i].cvel[j];
	      data[i].cvel[j] = tmp[i].cvel[j];
          }
#ifdef HALOID
          data[i].id = tmp[i].id;
#endif
	      data[i].rad = data[i].pos[0] * data[i].pos[0] + data[i].pos[1] * data[i].pos[1] + data[i].pos[2] * data[i].pos[2];
      }
    free(tmp);
    qsort(data, newn200, sizeof(particle), rad_sort_particle);

#ifdef DEBUG
fprintf(Logfile,"In halo %d, after unbinding M200=%g, r200=%g, N200=%d\n",grp ,newn200*PartMass, Halo_R_Crit200[grp],newn200);
fflush(Logfile);
#endif
#endif

    mvir = Halo_M_Crit200[grp];
    nvir=(int)(mvir/PartMass);

    for( i=0; i<nvir; i++ )
        if( PartMass * i > 0.5 * Halo_M_Crit200[grp])
        {
            Rhalf = sqrt( data[i].rad );
            break;
        }

#ifdef DEBUG
fprintf(Logfile,"In halo %d, after calcr200, rhalf=%g\n",grp, Rhalf);
fflush(Logfile);
#endif

    Vmax=0;
    Rmax=0;
#ifdef UNBINDING
    for (i = 0; i < newn200; i++) 
#else
    for (i = 0; i < count; i++) 
#endif
	if(i > 20) 
	  {
	    vv=(i+1)/sqrt(data[i].rad);
          if(vv>Vmax)
          {
          Vmax = vv;
          Rmax = sqrt( data[i].rad );    
          }
	  }
    Vmax=sqrt(G*PartMass*Vmax/Time); // phy unit

#ifdef DEBUG
fprintf(Logfile,"In halo %d, after calcr200, Vmax=%g, Rmax=%g\n",grp, Vmax, Rmax);
fflush(Logfile);
#endif

#ifdef UNBINDING
    return newn200;
#else
    return count;
//return (int)(Halo_M_Crit200[grp]/PartMass);
#endif
}

int rad_sort_particle(const void *a, const void *b)
{
    if (((particle *) a)->rad > ((particle *) b)->rad)
	return +1;

    if (((particle *) a)->rad < ((particle *) b)->rad)
	return -1;

    return 0;
}

double fof_periodic_wrap(double x)
{
    while (x >= BoxSize)
	x -= BoxSize;

    while (x < 0)
	x += BoxSize;

    return x;
}


double fof_periodic(double x)
{
    if (x >= 0.5 * BoxSize)
	x -= BoxSize;

    if (x < -0.5 * BoxSize)
	x += BoxSize;

    return x;
}
/*
#ifdef DEBUG
fprintf(Logfile,"In halo %d, Np=%d, r=%g, Vr=%g\n", grp, i, sqrt(data[i].rad), sqrt(G*PartMass*vv/Time) );
fflush(Logfile);
#endif
*/
