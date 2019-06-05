#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <hdf5.h>
#include "allvars.h"


int main(int argc, char **argv)
{

if(argc != 4)
    {
          fprintf(stderr, "\n  usage: HBT2HaloTrace+outputdir+snap+filename\n\n");
      exit(1);
    }

/////////////////////////////////////////////////////////////////////////////////
int snap,i,j,k,FILENR;
int TotNbound=0;
char filename[256],output[256],buf[500],buf1[500];
strcpy(output, argv[1]);
snap=atoi(argv[2]);
strcpy(filename, argv[3]);
nsub=loadsubcat(output,snap);
printf("snap_%03d Nsub=%d\n",snap,nsub);
FILE *fp;
float Mp;
/*******************Load Trackid****************************/
sprintf(buf, "%s/%s",output,filename);
//printf("%s\n",buf);
if((fp=fopen(buf,"r"))==NULL)
{	
		printf("Cann't open  file %s\n",buf);
		return 1;

}
  FILENR = scanline(buf);
  Halo= (HA *)malloc(sizeof(HA) * FILENR);
    for(i=0;i<FILENR;i++)
    {
             fscanf(fp,"%d",&Halo[i].TrackId);
    }  
  fclose(fp);
sprintf(buf1,"%s/trackid_Subfind_Used_%03d",output,snap);
if((fp=fopen(buf1,"w"))==NULL)
{	
		printf("Cann't open  file\n");
		return 1;

}
////////get Mp//////////
for (i = 0; i <nsub; i++)
if( sub[i].Nbound > 3  )
    {
Mp = sub[i].Mbound / sub[i].Nbound;
        break;
    }
///////////////////////////////////Calc fsub///////////////////////////////////////////
for (j = 0; j <FILENR; j++)
{
TotNbound=0;
for (i = 0; i <nsub; i++)
{
if( sub[i].HostHaloId == sub[Halo[j].TrackId].HostHaloId  )
if((sqrt( pow( fof_periodic( sub[i].ComovingMostBoundPosition[0] - sub[Halo[j].TrackId].ComovingMostBoundPosition[0]), 2.0)  +pow( fof_periodic( sub[i].ComovingMostBoundPosition[1] - sub[Halo[j].TrackId].ComovingMostBoundPosition[1]), 2.0) + pow( fof_periodic( sub[i].ComovingMostBoundPosition[2]-sub[Halo[j].TrackId].ComovingMostBoundPosition[2]),2.0) ) <= sub[Halo[j].TrackId].R200CritComoving )        &&        (sub[i].Rank != 0)    )
TotNbound+=sub[i].Nbound;
}
Halo[j].Msub = (float)TotNbound * Mp;
}
//////////////////////////////////////////////////////////////////////////////////////
/* pick  halo  */
for (j = 0; j <FILENR; j++)
    {
i=Halo[j].TrackId;
    fprintf(fp,"%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\n",sub[i].ComovingMostBoundPosition[0],sub[i].ComovingMostBoundPosition[1],sub[i].ComovingMostBoundPosition[2],sub[i].R200CritComoving,sub[i].M200Crit, Halo[j].Msub / (sub[i].M200Crit + Halo[j].Msub), sub[i].TrackId,sub[i].SnapshotIndexOfBirth,sub[i].Rank);
    } 
fclose(fp);

return 0;
}
































