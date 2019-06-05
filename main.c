#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "allvars.h"
#include "proto.h"


int main(int argc, char **argv)
{
  int gr, begSnapshotNum, begFileNr, endSnapshotNum, endFileNr, flag = 0, rep, i, j, realn;
  float fsub; int nsub,isub;
  float lowmass,omegaz,z;
  halostruct halo;
  char buf[500];
  FILE *fd;
#ifdef HALOID
  char hbuf[500];
  FILE *hfd;
#endif
#ifdef PROJECTION
  char pbuf[500];
  FILE *pfd;
#endif
  float *temp, *temp1, *temp2, vd;
  int nbin;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  if(argc != 5)
    {
      if(ThisTask == 0)
	{
	  fprintf(stderr,"\nusage: CHOPPERS OutputPath  Snapnumber Filenum Filename \n");
	  fprintf(stderr, "  OutputPath      the snapshot and subcatalogue path\n");
	  fprintf(stderr, "  Snapnumber      begin with this snapshot number\n");
	  fprintf(stderr, "  Filenumber      in one snapshot the subcatalogue file number\n");
	  fprintf(stderr, "  Filename        save halo name,default is halo_XXX\n");
	}
      MPI_Finalize();
      exit(1);
    }

  strcpy(OutputDir, argv[1]);
  begSnapshotNum = atoi(argv[2]);
  endSnapshotNum = begSnapshotNum;
  begFileNr = 0;
  endFileNr = atoi(argv[3]) - 1;
  strcpy(haloname, argv[4]);
  FilesPerSnapshot = endFileNr + 1;

  sprintf(SnapshotFileBase, "snap"); // snap basename

if(ThisTask == 0)
{
  printf("\nThis is CHOPPERS: Calculate gadget format spHerical halO ProPERtieS \n");
  printf("\nRunning on %d processors\n",NTask);
  printf("\nCode was compiled with setting:\n%s\n\n",COMPILETIMESETTINGS);
  fflush(stdout);
}

  set_units();

#ifdef DEBUG
  sprintf(buf, "Task_%03d", begSnapshotNum);
  mkdir(buf, 02755);
  sprintf(buf, "Task_%03d/Task_%03d.%d", begSnapshotNum,begSnapshotNum, ThisTask);
  Logfile = fopen(buf, "w");
#endif

  for(rep = 0; rep < NTask; rep++)
    {
      if(ThisTask == rep)
	flag = get_next_file(begSnapshotNum, begFileNr, endSnapshotNum, endFileNr, &SnapshotNum, &FileNr);
      MPI_Barrier(MPI_COMM_WORLD);
    }



  load_hash_table();

  lowmass=PartMass*NMIN;
  DELTA=overdensity(Omega0,1.0/Time-1.0);
  z=1.0/Time-1.0;
  omegaz = Omega0 * pow(1 + z,3) / (Omega0 * pow(1 + z, 3) + (1 - Omega0 - OmegaLambda) * pow(1 + z, 2) + OmegaLambda);
  RHO_CRIT=(3 * Omega0 * Hubble * Hubble / (8 * M_PI * G))/omegaz;

  while(flag == 0)
    {
        
#ifdef DEBUG
      fprintf(Logfile, "doing on Task=%d, SnapshotNum=%d, z=%g, rhoc=%g, FileNr=%d\n\n", ThisTask, SnapshotNum, z, RHO_CRIT,  FileNr);
      fflush(Logfile);
#endif
      realn=0;

      load_sub_catalogue();

      sprintf(buf, "%s/%s_%03d/halo_%03d.%d", OutputDir, haloname, SnapshotNum, SnapshotNum, FileNr);
      if(!(fd = fopen(buf, "w")))
	{
	  printf("can't open file `%s'\n", buf);
	  endrun(1);
	}
#ifdef HALOID
      sprintf(hbuf, "%s/%s_%03d/halo_id_%03d.%d", OutputDir, haloname, SnapshotNum, SnapshotNum, FileNr);
      if(!(hfd = fopen(hbuf, "w")))
	{
	  printf("can't open file `%s'\n", hbuf);
	  endrun(1);
	}
#endif
#ifdef PROJECTION
      sprintf(pbuf, "%s/%s_%03d/halo_pic_%03d.%d", OutputDir, haloname, SnapshotNum, SnapshotNum, FileNr);
      if(!(pfd = fopen(pbuf, "w")))
	{
	  printf("can't open file `%s'\n", pbuf);
	  endrun(1);
	}
#endif
      fwrite(&Ngroups,sizeof(int),1,fd);
      nbin = DENSNBIN;
      fwrite(&nbin,sizeof(int),1,fd);
      nbin = VELNBIN;
      fwrite(&nbin,sizeof(int),1,fd);

#ifdef HALOID
    int Ngr, idcount=0, offcount=0, Nid=0;    
    int *plen;
    int *poffset;
    long long *tmpid;
      for(Ngr=0, gr = 0; gr < Ngroups; gr++)
	  if(Halo_M_Crit200[gr] >= lowmass)
	  Ngr++;
      plen = malloc( sizeof(int) * Ngr );
      poffset = malloc( sizeof(int) * Ngr );
      tmpid = malloc( sizeof(long long) * Ngr * 2000000 );   //  assume halo have 2 million particle to allocate memory.
#ifdef DEBUG
fprintf(Logfile,"Save haloid, Ngr=%d\n", Ngr);
fflush(Logfile);
#endif
#endif

#ifdef PROJECTION
      int Pg;    
      for(Pg=0, gr = 0; gr < Ngroups; gr++)
	  if(Halo_M_Crit200[gr] >= lowmass)
	  Pg++;
#endif

      for(gr = 0; gr < Ngroups; gr++)
	{
	  if(Halo_M_Crit200[gr] < lowmass) continue;
	  np = get_spherical_region_coordinates(gr);
#ifdef HALOID
      Nid += np;
      plen[gr] = np; 
      for(j = 0; j < np; j++)
      tmpid[idcount++] = data[j].id;
      poffset[gr] = offcount; 
      offcount += np;
#endif
	  isub = FirstSubOfHalo[gr];

      for(j = 0; j < 3; j++)
	  halo.pos[j] = SubPos[3*isub+j];

	  halo.m200 = Halo_M_Crit200[gr];
	  halo.n200 = (int)(Halo_M_Crit200[gr] / PartMass );
	  halo.r200 = Halo_R_Crit200[gr];
	  halo.rhalf = Rhalf;
	  halo.v200 = sqrt( G * halo.m200 / halo.r200 );
	  halo.vmax = Vmax;
	  halo.rmax = Rmax;
      temp = malloc(sizeof(float)*3);
      temp1 = malloc(sizeof(float)*3);
      temp2 = malloc(sizeof(float)*3);
      halo.spinB = bspin(gr, temp, temp1, temp2, &vd); 
	  for(i = 0;i < 3;i++)
      {
      halo.sam[i] = temp[i];
      halo.vel[i] = temp1[i];
      halo.cm[i]  = temp2[i];
      }
      halo.vdisp = vd;
      free(temp);
      free(temp1);
      free(temp2);
      temp = malloc(sizeof(float)*3);
      temp1 = malloc(sizeof(float)*6);
      shape(  temp, temp1 );
      for(j = 0; j < 6; j++)
	  halo.InertialTensor[j] = temp1[j];
      halo.ea = temp[0];
      halo.eb = temp[1];
      halo.ec = temp[2];
      free(temp);
      free(temp1);
	  density_profile(gr);
	  halo.rconv = rconv;
      halo.c = cnfw(gr);
      halo.cvmax = cvmax( gr, halo.vmax, halo.v200 );
      halo.PE = 0.0;
      halo.KE = kinetic( );
#ifdef HBT 
      halo.fsub = SubSpin[3*isub];
#else
	  halo.fsub = get_fsub(gr,&nsub); 
#endif
      halo.soff = 0.0;
      for (j = 0; j < 3; j++) 
	  halo.soff  += pow( fof_periodic( halo.pos[j] - halo.cm[j] ), 2.0 );
      halo.soff = sqrt( halo.soff ) / halo.r200;
      halo.virial = 0.0;
#ifdef OUTPE 
      halo.PE = potential( );
      halo.virial = 2.0 * halo.KE / fabs( halo.PE );
      halo.spinP = pspin( gr, halo.KE, halo.PE );   // need calc Potential energy and Kinetic energy
#endif
	  for(i=0;i<DENSNBIN;i++)
	    {
	  halo.denr[i]=denr[i];
	  halo.denrho[i]=denrho[i];
	  halo.npbin[i]=npbin[i];
	    }
	  vc_profile(gr);
	  for(i=0;i<VELNBIN;i++)
	    {
	  halo.velr[i] = velr[i];
	  halo.velv[i] = velv[i];
	  halo.velnp[i] = velnp[i]; 
	    }
#ifdef HBT 
	  halo.trackid = SubLen[isub];
	  halo.birth = SubOffset[isub];
#else
	  halo.trackid = 0;
	  halo.birth = 0;
#endif

#ifdef DEBUG
fprintf(Logfile,"Halo %d, trackid=%d, birthsnap=%d, New M200=%g, r200=%g, N200=%d, rhalf=%g, cnfw=%g, cvmax=%g, spin=%g\n, pos=%g, %g, %g, vel=%g, %g, %g, cmpos=%g, %g, %g, sam=%g, %g, %g\nInertialTensor=",gr, halo.trackid, halo.birth, halo.m200, halo.r200, (int)(halo.m200/PartMass), halo.rhalf, halo.c, halo.cvmax, halo.spinB, halo.pos[0], halo.pos[1], halo.pos[2], halo.vel[0], halo.vel[1], halo.vel[2], halo.cm[0], halo.cm[1], halo.cm[2], halo.sam[0], halo.sam[1], halo.sam[2]);
for(j = 0; j < 6; j++)
fprintf(Logfile,"%g, ", halo.InertialTensor[j] );
fprintf(Logfile,"\nfirst axis=%g, second axis=%g, thrid axis=%g,V200=%g, Vmax=%g, rmax=%g, Vdispersion=%g, soff=%g, kinetic=%g\n\n", halo.ea, halo.eb, halo.ec, halo.v200, halo.vmax, halo.rmax, halo.vdisp, halo.soff, halo.KE);
#ifdef OUTPE
fprintf(Logfile,"potential=%g, Virial ratio=%g, Peebles spin=%g\n", halo.PE, halo.virial, halo.spinP);
#endif
fflush(Logfile);
#endif

	  realn++;
	  fwrite(&halo,sizeof(halo),1,fd);
	}
      rewind(fd);
      fwrite(&realn,sizeof(int),1,fd);
      fclose(fd);

#ifdef HALOID
        fwrite(&Ngr,sizeof(int),1,hfd);               //  Nhalo
        fwrite(&Nid,sizeof(int),1,hfd);               //  Nid
	    for(i=0;i<Ngr;i++)
        fwrite(&plen[i],sizeof(int),1,hfd);           //  Len
	    for(i=0;i<Ngr;i++)
        fwrite(&poffset[i],sizeof(int),1,hfd);        //  Offset
	    for(i=0;i<Nid;i++)
        fwrite(&tmpid[i],sizeof(long long),1,hfd);    //  Id
        free(plen);
        free(poffset);
        free(tmpid);
        fclose(hfd);
#endif
        free(data); 

#ifndef UNBINDING
#ifdef PROJECTION
        float lbox;
        lbox=2.0*HaloRadii;
        ngrid =NGRID;
        fwrite(&Pg,sizeof(int),1,pfd);        //  nhalo
        fwrite(&ngrid,sizeof(int),1,pfd);     //  ngrid
        fwrite(&lbox,sizeof(float),1,pfd);    //  lbox
        int tpic=0;    
        for(gr = 0; gr < Pg; gr++)
	    {
	    if(Halo_M_Crit200[gr] < lowmass) continue;
	    tpic=get_halo_cic(gr);
        for(i= 0; i< ngrid*ngrid; i++)       
        fwrite(&CICxy[i],sizeof(double),1,pfd);     // pic_xy
        for(i= 0; i< ngrid*ngrid; i++)       
        fwrite(&CICyz[i],sizeof(double),1,pfd);     // pic_yz
        for(i= 0; i< ngrid*ngrid; i++)       
        fwrite(&CICxz[i],sizeof(double),1,pfd);     // pic_xz
        }
        fclose(pfd);
        free(CICxy);
        free(CICxz);
        free(CICyz);
#endif
#endif

      if(realn !=0) free_particle_data();
      flag = get_next_file(begSnapshotNum, begFileNr, endSnapshotNum, endFileNr, &SnapshotNum, &FileNr);
    }
#ifdef DEBUG
  fclose(Logfile);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}


int get_next_file(int begsnap, int begfilenr, int endsnap, int endfilenr, int *snapshot, int *filenr)
{
  int snap, fnr, ret;
  char buf[1000];
  struct stat filestatus;
  FILE *fd;
  double age;
#ifdef HALOID
  char hbuf[500];
  FILE *hfd;
#endif
#ifdef PROJECTION
  char pbuf[500];
  FILE *pfd;
#endif

  snap = begsnap;
  fnr = begfilenr;
  do
    {
      sprintf(buf, "%s/%s_%03d/", OutputDir, haloname, snap);
      mkdir(buf, 02755);

      sprintf(buf, "%s/%s_%03d/halo_%03d.%d", OutputDir, haloname, snap, snap, fnr);
#ifdef HALOID
      sprintf(hbuf, "%s/%s_%03d/halo_id_%03d.%d", OutputDir, haloname, snap, snap, fnr);
#endif
#ifdef PROJECTION
      sprintf(pbuf, "%s/%s_%03d/halo_pic_%03d.%d", OutputDir, haloname, snap, snap, fnr);
#endif

      ret = stat(buf, &filestatus);

      if(ret != 0)		/* seems not to exist */
	{
	  if(!(fd = fopen(buf, "w")))
	    {
	      printf("can't create file '%s'\n", buf);
	      endrun(787);
	    }
	  fclose(fd);
#ifdef HALOID
	  if(!(hfd = fopen(hbuf, "w")))
	    {
	      printf("can't create file '%s'\n", hbuf);
	      endrun(787);
	    }
	  fclose(hfd);
#endif
#ifdef PROJECTION
	  if(!(pfd = fopen(pbuf, "w")))
	    {
	      printf("can't create file '%s'\n", pbuf);
	      endrun(787);
	    }
	  fclose(pfd);
#endif

	  *snapshot = snap;
	  *filenr = fnr;
	  return 0;
	}

      /* ok, the file existed */

      if(filestatus.st_size == 0)	/* it is only a lock file. Check how old it is */
	{
	  age = difftime(time(NULL), filestatus.st_ctime);

	  if(age > 24 * 3600.0 * 100.0 )	/* older than 100 days? Probably wasn't finished. Let's try again */
	    {
	      if(!(fd = fopen(buf, "w")))
		{
		  printf("can't create file '%s'\n", buf);
		  endrun(787);
		}
	      fclose(fd);

	      *snapshot = snap;
	      *filenr = fnr;
	      return 0;
	    }
	}

      fnr = fnr + 1;
      if(fnr >= FilesPerSnapshot)
	{
	  fnr = 0;
	  snap++;
	}
    }
  while(snap <= endsnap && fnr <= endfilenr);
  return 1;
}


void endrun(int ierr)
{
  if(ierr)
    {
      printf("task %d: endrun called with an error level of %d\n\n\n", ThisTask, ierr);
      MPI_Abort(MPI_COMM_WORLD, ierr);
      exit(0);
    }

  MPI_Finalize();
  exit(0);
}


float get_fsub(int nhalo,int *nsub)
{
  int i,subid,nn=0;
  int msub = 0;
  float dx, dy, dz;
  float xc,yc,zc,rvir,rvir2;
  float mvir;

  subid=FirstSubOfHalo[nhalo];
  xc=SubPos[3*subid];
  yc=SubPos[3*subid+1];
  zc=SubPos[3*subid+2];

  rvir=Halo_R_Crit200[nhalo];
  mvir=Halo_M_Crit200[nhalo];
  rvir2=rvir*rvir;
  
  for (i = FirstSubOfHalo[nhalo] + 1; i < FirstSubOfHalo[nhalo] + NsubPerHalo[nhalo]; i++) 
    {
      if (SubParentHalo[i] != nhalo) 
	{
	  printf("subhalo id not matched\n");
	  endrun(898);
	}
    
      dx = SubPos[3 * i + 0] - xc;
      dy = SubPos[3 * i + 1] - yc;
      dz = SubPos[3 * i + 2] - zc;
	
      if ((dx * dx + dy * dy + dz * dz) < rvir2)
	{
	  msub += SubLen[i];
	  nn++;
	}
    }
  *nsub=nn;
  return msub*PartMass/mvir;
}

float overdensity(float om_0,float z)
{
  float omz,d;
  omz=om_0*pow(1.0+z,3);
  omz/=om_0*pow(1.0+z,3)+(1-om_0);

  d=omz-1.0;
  return 18.0*M_PI*M_PI+82.0*d-39.0*d*d;
}

void set_units(void)
{
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);
  Hubble = HUBBLE * UnitTime_in_s;
  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
}
