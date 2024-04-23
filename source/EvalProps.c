#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"global.h"
#include<time.h>
#include<string.h>
#define IMUL    314159269
#define IADD    453806245
#define MASK    2147483647
#define SCALE   0.4656612873e-9

#define MASTER 0

void EvalProps(){
int n;
 real v, vv;
  vSum = 0.;
  vvSum = 0.;
  //int n;
  for(n = 1 ; n <= nAtom ; n ++){
    vv = 0.;
    v = vx[n] - 0.5 * deltaT * ax[n];
    vSum += v;
    vv += Sqr(v);
    v = vy[n] - 0.5 * deltaT * ay[n];
    vSum += v;
    vv += Sqr(v);
    //v = vz[n] - 0.5 * deltaT * az[n];
    //vSum += v;
    //vv += Sqr(v);
    vvSum += vv;
  }
  kinEnergy = 0.5 * vvSum / nAtom;
  potEnergy = uSum / nAtom;
  totEnergy = kinEnergy + potEnergy;
  //pressure = density * (vvSum + virSum)/(nAtom * NDIM);



}
