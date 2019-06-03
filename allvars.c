#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"

int ThisTask;
int NTask;

FILE *Logfile;

int SnapshotNum;
int FileNr;
int FilesPerSnapshot;
int HashTabSize;
int NFiles;
int *HashTable, *FileTable, *LastHashCell, *NInFiles;

int Nids;
int TotNgroups;

char OutputDir[512];
char haloname[500];
char SnapshotFileBase[512];
double Hubble;
double G;
double Omega0;
long long TotNumPart;

double OmegaLambda;
double HubbleParam;
double UnitTime_in_s;
double UnitDensity_in_cgs;
double UnitEnergy_in_cgs;

double BoxSize;

struct io_header header;

double Time;
int NumPart;
double PartMass;

int Ngroups;

int Nsubhalos, *SubParentHalo, *SubOffset, *SubLen, *FirstSubOfHalo, *NsubPerHalo;
float *SubPos, *SubVel, *SubVelDisp, *SubVmax, *SubSpin, *SubHalfMass;
float *Halo_M_Mean200, *Halo_R_Mean200, *Halo_M_Crit200, *Halo_R_Crit200, *Halo_M_TopHat200, *Halo_R_TopHat200;
long long *SubMostBoundID;


particle *data;

float Vmax, Rmax, Rhalf;

float DELTA,rconv;

int np, npbin[DENSNBIN], velnp[VELNBIN];
float logr[DENSNBIN],logrho[DENSNBIN];
float denr[DENSNBIN],denrho[DENSNBIN];
float velr[VELNBIN], velv[VELNBIN];

float RHO_CRIT;


double rrhos,rrs,connfw;

#ifdef PROJECTION
int ngrid;
double *CICxy, *CICyz, *CICxz;
#endif
