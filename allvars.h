#include <stdio.h>
//////////////////////////////////////////////////////////////////////////////////////////
#define  GRAVITY     6.672e-8
#define  HUBBLE      3.2407789e-18	/* in h/sec */
////////////////////////////Unit in Mpc/h or kpc/h/////////////////////////////////////////
#define UnitLength_in_cm                   3.08568e+24
#define UnitMass_in_g                      1.989e+43
#define UnitVelocity_in_cm_per_s           100000
///////////////////////////////////////Parameter///////////////////////////////////////////
#define NMIN             1              //   halo resolution
#define DENSNBIN         20             //   density profile bin number, fit NFW also use it.
#define KAPPA            1.0            //   P03 rconv kappa
#define Rpromin          0.05           //   output density profile (Vc) min = Rpromin*R200, fit NFW also use it.
#define Rpromax          1.0            //   output density profile (Vc) max = Rpromax*R200
#define Ncut_fitnfw      1000           //   if Np < 1000 halo  concentration=0, else use fitnfw.c or vmaxfit
#define VELNBIN          20             //   density profile (Vc) bin number, fit NFW also use it.
#define NGRID            1024           //   CIC grid number
#define HaloRadii        4.0            //   unit in Mpc/h    2.0*HaloRadii is CIC BoxSize 


extern int  ThisTask;
extern int  NTask;
extern int  FilesPerSnapshot;
extern int  SnapshotNum;
extern int  FileNr;
extern int  HashTabSize;
extern int  NFiles;
extern int  *HashTable, *FileTable, *LastHashCell, *NInFiles;

extern int  Ngroups;
extern int  Nids;
extern int  TotNgroups;

extern long long TotNumPart;

extern char OutputDir[512];
extern char haloname[500];
extern char SnapshotFileBase[512];

extern double Hubble;
extern double G;
extern double Omega0;

extern double OmegaLambda;
extern double HubbleParam;

extern double UnitTime_in_s;
extern double UnitDensity_in_cgs;
extern double UnitEnergy_in_cgs;

extern FILE *Logfile;

extern double BoxSize;

typedef long long peanokey;

extern struct io_header
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  int hashtabsize;
  char fill[84];		/* fills to 256 Bytes */
}
header;

extern double  Time;
extern int     NumPart;
extern double  PartMass;

/* subfind   */
extern int Nsubhalos, *SubParentHalo, *SubOffset, *SubLen, *FirstSubOfHalo, *NsubPerHalo;
extern float *SubPos, *SubVel, *SubVelDisp, *SubVmax, *SubSpin, *SubHalfMass;
extern float *Halo_M_Mean200, *Halo_R_Mean200, *Halo_M_Crit200, *Halo_R_Crit200, *Halo_M_TopHat200, *Halo_R_TopHat200;
extern long long *SubMostBoundID;

typedef struct 
{
  float pos[3];            /*position*/
  float vel[3];            /*velocity*/
  float cpos[3];            /*comoving position*/
  float cvel[3];            /*comoving velocity*/
  float rad;               /* (r-r_center)  */
  float Phi;               /*particle  potential for unbinding*/
#ifdef HALOID
  long long id;            /* particle id default is long long type, if use int type, change it in all code */
#endif
} particle;

typedef struct 
{
  float pos[3];	             /* halo position */
  float vel[3];              /* halo velocity in physical unit */
  float sam[3];              /* halo specific angular momentum */
  float cm[3];               /* halo center of mass */
  float InertialTensor[6];   /* halo Inertial Tensor {Ixx, Ixy, Ixz, Iyy, Iyz, Izz} */
  float m200;                /* default halo mass */
  int   n200;                /* default halo particle number*/
  float r200;                /* default halo radius  */
  float v200;                /* v200 = sqrt( G * M200 / r200 )   */
  float vmax;                /* maximum of rotation curve */
  float rmax;                /* radius of rotation curve maximum */
  float rhalf;               /* radius of half m200 */
  float spinB;               /* Bullock’ 2001 spin parameter */
  float spinP;               /* Peebles’ 1969 spin parameter*/
  float rconv;               /* Power 2003 convergence radius */
  float c;	                 /* halo concentration, use GSL fit NFW */
  float cvmax;               /* halo concentration, Prada et al. 2012 definition */
  float vdisp;               /* halo velocity dispersion */
  float ea;                  /* first largest axis of moment of inertia tensor */
  float eb;                  /* second largest axis of moment of inertia tensor */
  float ec;                  /* third largest axis of moment of inertia tensor  */
  float PE;                  /* halo potential energy */
  float KE;                  /* halo kinetic energy  */
  float fsub;	             /* Msub / M200, NOTE: use HBT2 need input this parameter  */
  float soff;                /* Thomas et al. 2001 center of mass displacement, soff= | r_minpot - r_cm | / r200 */
  float virial;              /* Neto et al. 2007 Virial ratio: virial = 2KE / | PE | */
  float denr[DENSNBIN];      /* halo density profile r = r / r200 */
  float denrho[DENSNBIN];    /* halo density profile rho = rho / rhoc in physical unit  */
  int   npbin[DENSNBIN];     /* particle number in current spherical shell*/
  float velr[VELNBIN];       /* halo circular velocity  r = r / r200 */ 
  float velv[VELNBIN];       /* halo circular velocity Vc*/
  int   velnp[VELNBIN];      /* number of particles inside sphere of radius r*/
  int   trackid;             /* NOTE: HBT2 halo trackid, subfind no use  */
  int   birth;               /* NOTE: HBT2 halo snapshot of birth, when the halo first becomes resolved, subfind no use */
} halostruct;

extern particle *data;

extern float Vmax, Rmax, Rhalf;

extern float DELTA;

extern float rconv;

extern int np, npbin[DENSNBIN], velnp[VELNBIN];
extern float logr[DENSNBIN], logrho[DENSNBIN]; // for fit nfw
extern float denr[DENSNBIN], denrho[DENSNBIN];
extern float velr[VELNBIN], velv[VELNBIN];
extern float RHO_CRIT;   

extern double rrhos,rrs;
typedef struct 
{
	int nbins;	
	float *rr;
	float *lgrho;
	float rvir;
} HALOPROFILE;

struct Hdata {
   size_t nbins;
   double * r;
   double * rho;//log(rho)
   double * sigma;//error for log(rho)
 };

#ifdef PROJECTION
extern int ngrid;
extern double *CICxy, *CICyz, *CICxz;
#endif
