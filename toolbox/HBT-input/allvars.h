
#include "proto.h"

#define BBoxSize 200

typedef struct subhalo
{
  int TrackId;
  float Mbound;
  int Nbound;
  int HostHaloId;
  int Rank;
  float R200CritComoving;
  float M200Crit;
  float ComovingMostBoundPosition[3];//default pos of sub
  int SnapshotIndexOfBirth;
  int SnapshotIndexOfDeath;
  float SpinBullock[3];
  float fsub;       //  totMsub(<r200) / M200 for halo not subhalo
}  SUBCAT;

extern SUBCAT *sub;

extern int nsub;

typedef struct alo
{
  int TrackId;
  float Msub;
}  HA;

extern HA *Halo;
