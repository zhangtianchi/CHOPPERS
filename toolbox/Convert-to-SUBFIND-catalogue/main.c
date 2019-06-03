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

typedef struct subhalo
{
  int TrackId;
/*SUBFIND halo property*/   
  int NsubPerHalo;
  int FirstSubOfHalo;
  int SubLen;           //   replaced by TrackId
  int SubOffset;        //   replaced by Birthsnap;
  int SubParentHalo;   //    replaced by Rank
  float Halo_M_Mean200;
  float Halo_R_Mean200;
  float Halo_M_TopHat200;
  float Halo_R_TopHat200;
  float SubVel[3];
  float SubVelDisp[3];
  float SubVmax;
  float SubSpin[3];    //    replaced by fsub
  long long SubMostBoundID;
  float SubHalfMass;
  float R200;
  float M200;
  float pos[3];//default pos of sub
  int Birthsnap;
  int Rank;
  float Spin;
}  SUBCAT;

int scanline(char * filename)
{
        FILE * f = 0; char line[256]=""; int lines = 0;
        f = fopen(filename, "r");
        if(!f) return 0;
        while(!feof(f)) { fgets(line, 99999999, f); lines++; }
        fclose(f);
        return lines-1;

}

int main(int argc, char **argv)
{
  
  char buf[500];
  FILE *fd;
  SUBCAT *halo;
  int FILENR; 
  int i,j,k;
    i=j=k=0;
  int ThisTask,NTask;
  int Ngroups,Nids,TotNgroups,NFiles,Nsubhalos;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

if(argc != 5)
    {
          fprintf(stderr, "\n  usage: HBT2LSubfind+inputdir+outputdir+snap+filename\n\n");
      exit(1);
    }

/////////////////////////////////////////////////////////////////////////////////
    Nids = 1;
    NFiles = NTask;
int snap;
char output[256];
char input[256];
char ffile[256];
strcpy(input, argv[1]);
strcpy(output, argv[2]);
snap=atoi(argv[3]);
strcpy(ffile, argv[4]);
  sprintf(buf, "%s/%s_%03d",input,ffile,snap);
if((fd=fopen(buf,"r"))==NULL)
{	
		printf("Cann't open  file %s\n",buf);
		return 1;

}
  FILENR = scanline(buf);
  halo= (SUBCAT *)malloc(sizeof(SUBCAT) * FILENR);
    TotNgroups = FILENR;
    if(ThisTask == 0)
    printf("NFiles=%d,TotNgroups=%d\n",NFiles,FILENR);

    
    for(i=0;i<FILENR;i++)
    {
             fscanf(fd,"%e %e %e %e %e %e %d %d %d",&halo[i].pos[0],&halo[i].pos[1],&halo[i].pos[2],&halo[i].R200,&halo[i].M200,&halo[i].Spin,&halo[i].TrackId,&halo[i].Birthsnap,&halo[i].Rank);  
    }  
  fclose(fd);
    for(k=0,i=0;i<FILENR;i++)
    {
       halo[i].NsubPerHalo = 1;
       halo[i].SubLen = halo[i].TrackId;
       halo[i].SubOffset = halo[i].Birthsnap;
       halo[i].SubParentHalo = halo[i].Rank;
       halo[i].Halo_M_Mean200 = 0.0;
       halo[i].Halo_R_Mean200 = 0.0;
       halo[i].Halo_M_TopHat200 = 0.0;
       halo[i].Halo_R_TopHat200 = 0.0;
      for(j=0;j<3;j++)
        {
            halo[i].SubVel[j]=0.0;
            halo[i].SubVelDisp[j]=0.0;
            halo[i].SubSpin[j]=0.0;
        }
       halo[i].SubSpin[0] = halo[i].Spin;
       halo[i].SubVmax = 0.0;
       halo[i].SubMostBoundID = 0L;
       halo[i].SubHalfMass = 0.0;
    if(i % NTask == 0 )
        {
        for(j=i;j<i+NTask && j<FILENR;j++)   //  ensure j < FILENR , avoid segement error
       halo[j].FirstSubOfHalo = k ;
            k++;
        } 
    }

/********************************************every core save data******************************************************/
    for(k=0,i=ThisTask;i<FILENR;i+=NTask)
    {
        k++;
    }
    Ngroups = Nsubhalos = k;
   printf("Task %d: Ngroups = %d\n",ThisTask,Nsubhalos);
    if(ThisTask == 0)
    {
      sprintf(buf, "%s/postproc_%03d/", output, snap);
      mkdir(buf, 02755);
    }
    MPI_Barrier(MPI_COMM_WORLD);	// wait to make sure that directory has been created 

      sprintf(buf, "%s/postproc_%03d/sub_tab_%03d.%d", output, snap, snap, ThisTask);
    if((fd=fopen(buf,"w"))==NULL)
    {	
		printf("Cann't open  file %s\n",buf);
		return 1;

    }

  fwrite(&Ngroups, sizeof(int), 1, fd);
  fwrite(&Nids, sizeof(int), 1, fd);
  fwrite(&TotNgroups, sizeof(int), 1, fd);
  fwrite(&NFiles, sizeof(int), 1, fd);
  fwrite(&Nsubhalos, sizeof(int), 1, fd);
    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].NsubPerHalo, sizeof(int), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].FirstSubOfHalo, sizeof(int), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].SubLen, sizeof(int), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].SubOffset, sizeof(int), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].SubParentHalo, sizeof(int), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].Halo_M_Mean200, sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].Halo_R_Mean200, sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].M200, sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].R200, sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].Halo_M_TopHat200, sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].Halo_R_TopHat200, sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
    for(j=0;j<3;j++)
  fwrite(&halo[i].pos[j], sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
    for(j=0;j<3;j++)
  fwrite(&halo[i].SubVel[j], sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].SubVelDisp, sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].SubVmax, sizeof(float), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
    for(j=0;j<3;j++)
  fwrite(&halo[i].SubSpin[j], sizeof(float), 1, fd);
    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].SubMostBoundID, sizeof(long long), 1, fd);

    for(i=ThisTask;i<FILENR;i+=NTask)
  fwrite(&halo[i].SubHalfMass, sizeof(float), 1, fd);
   fclose(fd);
  free(halo);
  MPI_Finalize();
  return 0;
}
