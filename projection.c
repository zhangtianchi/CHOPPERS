#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"

#ifdef PROJECTION

/***************Get Halo region CICgrid*************************/
int Index(int x, int y, int z)
{
	return x * NGRID * NGRID + y * NGRID + z;
}
double *Griding_CIC( particle *P, int N, int SIZE, double L )  // SIZE is Halo_Npart   N is grid number  L is boxsize
{
	int n;
	int ibin, jbin, kbin, ibinm, ibinp, jbinm, jbinp, kbinm, kbinp;
	double hx, hy, hz, hx0, hy0, hz0, hxp, hxm, hyp, hym, hzp, hzm;
	double value = 1.0;
	double *gridP = (double *)malloc(N * N * N * sizeof(double));
	memset(gridP, 0, sizeof(double) * N * N * N);  // initialization
	double H = (double)L / (double)N;

	for(n = 0; n < SIZE; n++)
	{
		ibin = P[n].pos[0] / H;
		jbin = P[n].pos[1] / H;
		kbin = P[n].pos[2] / H;
  
		if(ibin > (N - 1) || ibin < 0) 
        {
        printf("error: ibin exceeds boundary: halo particle id: %d, %d,%d,%d\n", n, ibin, jbin, kbin);
        printf("error: pos: %g, %g, %g\n", P[n].pos[0], P[n].pos[1], P[n].pos[2]);
        endrun(1);
        }
		if(jbin > (N - 1) || jbin < 0) 
        {
        printf("error: jbin exceeds boundary: halo particle id: %d, %d,%d,%d\n", n, ibin, jbin, kbin);
        printf("error: pos: %g, %g, %g\n", P[n].pos[0], P[n].pos[1], P[n].pos[2]);
        endrun(1);
        }
		if(kbin > (N - 1) || kbin < 0) 
        {
        printf("error: kbin exceeds boundary: halo particle id: %d, %d,%d,%d\n", n, ibin, jbin, kbin);
        printf("error: pos: %g, %g, %g\n", P[n].pos[0], P[n].pos[1], P[n].pos[2]);
        endrun(1);
        }
		hx = P[n].pos[0] / H - ibin - 0.5;
		hy = P[n].pos[1] / H - jbin - 0.5;
		hz = P[n].pos[2] / H - kbin - 0.5;

		hx0 = 1. - fabs(hx);
		hy0 = 1. - fabs(hy);
		hz0 = 1. - fabs(hz);

		if(hx > 0)
		{
			hxp = hx;
			hxm = 0.;
		}
		else
		{
			hxp = 0.;
			hxm = -hx;
		}

		if(hy > 0)
		{
			hyp = hy;
			hym = 0.;
		}
		else
		{
			hyp = 0.;
			hym = -hy;
		}

		if(hz > 0)
		{
			hzp = hz;
			hzm = 0.;
		}
		else
		{
			hzp = 0.;
			hzm = -hz;
		}

		ibinm = (ibin - 1 + N) % N;
		ibinp = (ibin + 1 + N) % N;
		jbinm = (jbin - 1 + N) % N;
		jbinp = (jbin + 1 + N) % N;
		kbinm = (kbin - 1 + N) % N;
		kbinp = (kbin + 1 + N) % N;

		gridP[Index(ibinm, jbinm, kbinm)] += hxm * hym * hzm * value;
		gridP[Index(ibinm, jbinm, kbin)] += hxm * hym * hz0 * value;
		gridP[Index(ibinm, jbinm, kbinp)] += hxm * hym * hzp * value;
		gridP[Index(ibinm, jbin , kbinm)] += hxm * hy0 * hzm * value;
		gridP[Index(ibinm, jbin , kbin)] += hxm * hy0 * hz0 * value;
		gridP[Index(ibinm, jbin , kbinp)] += hxm * hy0 * hzp * value;
		gridP[Index(ibinm, jbinp, kbinm)] += hxm * hyp * hzm * value;
		gridP[Index(ibinm, jbinp, kbin)] += hxm * hyp * hz0 * value;
		gridP[Index(ibinm, jbinp, kbinp)] += hxm * hyp * hzp * value;

		gridP[Index(ibin , jbinm, kbinm)] += hx0 * hym * hzm * value;
		gridP[Index(ibin , jbinm, kbin)] += hx0 * hym * hz0 * value;
		gridP[Index(ibin , jbinm, kbinp)] += hx0 * hym * hzp * value;
		gridP[Index(ibin , jbin , kbinm)] += hx0 * hy0 * hzm * value;
		gridP[Index(ibin , jbin , kbin)] += hx0 * hy0 * hz0 * value;
		gridP[Index(ibin , jbin , kbinp)] += hx0 * hy0 * hzp * value;
		gridP[Index(ibin , jbinp, kbinm)] += hx0 * hyp * hzm * value;
		gridP[Index(ibin , jbinp, kbin)] += hx0 * hyp * hz0 * value;
		gridP[Index(ibin , jbinp, kbinp)] += hx0 * hyp * hzp * value;

		gridP[Index(ibinp, jbinm, kbinm)] += hxp * hym * hzm * value;
		gridP[Index(ibinp, jbinm, kbin)] += hxp * hym * hz0 * value;
		gridP[Index(ibinp, jbinm, kbinp)] += hxp * hym * hzp * value;
		gridP[Index(ibinp, jbin , kbinm)] += hxp * hy0 * hzm * value;
		gridP[Index(ibinp, jbin , kbin)] += hxp * hy0 * hz0 * value;
		gridP[Index(ibinp, jbin , kbinp)] += hxp * hy0 * hzp * value;
		gridP[Index(ibinp, jbinp, kbinm)] += hxp * hyp * hzm * value;
		gridP[Index(ibinp, jbinp, kbin)] += hxp * hyp * hz0 * value;
		gridP[Index(ibinp, jbinp, kbinp)] += hxp * hyp * hzp * value;
	}

	return gridP;
}

int get_halo_cic(int grp)
{

#ifdef DEBUG
fprintf(Logfile,"In get_halo_cic function, begin Projection\n");
fflush(Logfile);
#endif

    int dummy, i,j,k, count,nvir;
    char buf[1000];
    FILE *fd;
    float rvir,mvir,vv;
    struct io_header header;
    int n200;
    int len,isub;
    float pos[3];
    particle *tm, *cicp;
    float *ppos;
    double cellsize;
    int cells, cx, cy, cz, ix, iy, iz, iix, iiy, iiz, base, hashbits, hashkey, hashtabsize;

    sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", OutputDir, SnapshotNum, SnapshotFileBase,SnapshotNum, 0);
    rvir=Halo_R_Crit200[grp];
    mvir=Halo_M_Crit200[grp];
    isub=FirstSubOfHalo[grp];
    pos[0]=SubPos[3*isub];
    pos[1]=SubPos[3*isub+1];
    pos[2]=SubPos[3*isub+2];

    if (!(fd = fopen(buf, "r"))) {
	printf("can't open file `%s'\n", buf);
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


    cells = (int) (  HaloRadii / cellsize + 2 );

#ifdef DEBUG
fprintf(Logfile,"In get_halo_cic function, cellsize=%g, HaloRadii/cellsize=%g, cells=%d\n",cellsize, HaloRadii/cellsize, cells);
fprintf(Logfile,"In get_halo_cic function, Base=%d, hashbits=%d\n",base,hashbits);
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

    ppos = (float *) malloc(3 * sizeof(float) * count);


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
		    printf("can't open file `%s'\n", buf);
		    exit(0);
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

#ifdef DEBUG
fprintf(Logfile,"In get_halo_cic function, allcount=%d\n",count);
fflush(Logfile);
#endif

    tm = (particle *) malloc(sizeof(particle) * count);
    memset(tm, 0, sizeof(particle) * count);  // initialization

    for (i = 0; i < count; i++)
    {
	tm[i].pos[0] = fof_periodic(ppos[i * 3 + 0] - pos[0]);
	tm[i].pos[1] = fof_periodic(ppos[i * 3 + 1] - pos[1]);
	tm[i].pos[2] = fof_periodic(ppos[i * 3 + 2] - pos[2]);
	tm[i].rad = tm[i].pos[0] * tm[i].pos[0] + tm[i].pos[1] * tm[i].pos[1] + tm[i].pos[2] * tm[i].pos[2];
    }
    free(ppos);
    qsort(tm, count, sizeof(particle), rad_sort_particle);

    count = 0;
    while (sqrt(tm[count].rad) < HaloRadii) 
 	count++;

#ifdef DEBUG
fprintf(Logfile,"In get_halo_cic function, less than %g Mpc/h Np=%d\n",HaloRadii,count);
fflush(Logfile);
#endif

    cicp = (particle *) malloc(sizeof(particle) * count);

    for (i = 0; i < count; i++) 
      {
	cicp[i].pos[0] = tm[i].pos[0] + HaloRadii;
	cicp[i].pos[1] = tm[i].pos[1] + HaloRadii;
	cicp[i].pos[2] = tm[i].pos[2] + HaloRadii;
      }
      free(tm);

double *CIC;
CIC = malloc( sizeof( double)*( NGRID* NGRID* NGRID) );
ngrid =NGRID; 
#ifdef DEBUG
fprintf(Logfile,"In get_halo_cic function, ngrid=%d, begin CIC\n",ngrid);
fflush(Logfile);
#endif

CIC=Griding_CIC(cicp,ngrid,count,2.0*HaloRadii);  

#ifdef DEBUG
fprintf(Logfile,"In get_halo_cic function, end CIC\n");
fflush(Logfile);
#endif
free(cicp);

CICxy=malloc(sizeof(double)*(NGRID*NGRID));
CICxz=malloc(sizeof(double)*(NGRID*NGRID));
CICyz=malloc(sizeof(double)*(NGRID*NGRID));
for(i= 0; i< NGRID; i++)       
for(j= 0; j< NGRID; j++)      
{
CICxy[i * NGRID  + j]=0.0;
CICxz[i * NGRID  + j]=0.0;
CICyz[i * NGRID  + j]=0.0;
}

for(i= 0; i< NGRID; i++)       
for(j= 0; j< NGRID; j++)      
for(k= 0; k< NGRID; k++)     
{
CICxy[i * NGRID  + j]+=CIC[i * NGRID * NGRID + j * NGRID + k];  //xy
CICxz[i * NGRID  + k]+=CIC[i * NGRID * NGRID + j * NGRID + k];  //xz
CICyz[j * NGRID  + k]+=CIC[i * NGRID * NGRID + j * NGRID + k];  //yz
}

free(CIC);
#ifdef DEBUG
fprintf(Logfile,"In get_halo_cic function, end Projection\n");
fflush(Logfile);
#endif

    return count;
}
#endif
