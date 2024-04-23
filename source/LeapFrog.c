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

void LeapFrog(){

if(stepCount <= stepEquil){
    double gSum, varS, massS;
    temperature = 1./GAMMA;
    
   if(stepCount == 1) varS = 0.;
    double A, S1, S2, T;
    int n;
    S1 = 0.; S2 = 0.; gSum = 0.; massS = 0.1; 

    vvSum = 0.;
    double halfdt = 0.5*deltaT;
    for (n = 1; n <= nAtom; n++){
      T = vx[n] + halfdt * ax[n];
      S1 += T * ax[n];
      S2 += Sqr(T);
      
      T = vy[n] + halfdt * ay[n];
      S1 += T * ay[n];
      S2 += Sqr(T);

      //T = vz[n] + halfdt * az[n];
      //S1 += T * az[n];
      //S2 += Sqr(T);
     vvSum += (Sqr(vx[n]) +  Sqr(vy[n]));
    }

    A = -S1 / S2;
    S2 = vvSum;
   
    double C = 1 + A*deltaT ;
    double D = deltaT * (1 + 0.5 * A * deltaT);
    
    int i,j;
    real dr[NDIM+1], r, rr, ri, rrCut;
    double vv;

    double uVal, uSum, fcVal, f, AA, AASum;
    double TVal, T_Config;

    double deno, VVSum;	 
    double PE;
    PE=0.;
    deno = 0.;
    VVSum = 0.;
    keConfig = 0.;
    AASum = 0.;
   
  for(n=1;n<=nAtom; n++)
     TValSum[n] = 0.;

  rrCut = Sqr(rCut);

/*****Calculating Configarational temperature*****/
if(thermo == 'C'){
for(i = 1 ; i <= nAtom; i ++){
    for(j = i+1 ; j <= nAtom ; j ++){
      dr[1] = rx[i] - rx[j];
      if(fabs(dr[1]) > regionH[1])
	dr[1] -= SignR(region[1], dr[1]);

     dr[2] = ry[i] - ry[j];
      if(fabs(dr[2]) > regionH[2])
	dr[2] -= SignR(region[2], dr[2]);
  
     rr = Sqr(dr[1]) + Sqr(dr[2]);
     if(rr < rrCut ){
     r = sqrt(rr);
     ri = 1/r;
     uVal = ri*exp(-kappa*r);

     TVal = (1./rr + Sqr(kappa) + kappa/r)*uVal;
     TValSum[i] += TVal;
     TValSum[j] += TVal;
   } }
     AA = Sqr(ax[i]) + Sqr(ay[i]);
     AASum += AA;
     vv = Sqr(vx[i]) + Sqr(vy[i]);
     VVSum += vv; 
     deno += TValSum[i];    
} 
     PE = uSum/nAtom; 
     keConfig = AASum/deno;

     double gSumconfig, varSconfig, massSconfig;
     if(stepCount == 1) varSconfig = 0.;
     gSumconfig = 0.; massSconfig = 2.0; 
   
     gSumconfig = (AASum/deno -temperature)/massSconfig;
     varSconfig += deltaT*gSumconfig; 

      /*****Configarational Nose-Hoover thermostat*****/
   for (n = 1; n <= nAtom; n++){
      vx[n] += deltaT * ax[n];
      rx[n] += deltaT * (vx[n]);
      vy[n] += deltaT * ay[n];
      ry[n] += deltaT * (vy[n]);
    }
      /*****Kinetic Nose-Hoover thermostat*****/
  }else if(thermo == 'N'){  
    gSum = (0.5*S2 - (nAtom + 1)*temperature)/massS;
    varS += deltaT*gSum; 
   for (n = 1; n <= nAtom; n++){
      vx[n] += deltaT * (ax[n] - varS *vx[n]);
      rx[n] += deltaT * vx[n];
      vy[n] += deltaT * (ay[n] - varS *vy[n]);
      ry[n] += deltaT * vy[n];
   }
      /*****for Gaussian thermostat*****/
 }else if(thermo == 'G'){              
      for (n = 1; n <= nAtom; n++){
      vx[n] = C * vx[n] + D * ax[n];
      rx[n] += deltaT * vx[n];
      vy[n] = C * vy[n] + D * ay[n];
      ry[n] += deltaT * vy[n];
      //vz[n] = C * vz[n] + D * az[n];
      //rz[n] += deltaT * vz[n];
     }
   }
 }else{
    int n;
//#pragma acc parallel loop copyin(ax[1:nAtom], ay[1:nAtom], az[1:nAtom], deltaT) copy(vx[1:nAtom], vy[1:nAtom], vz[1:nAtom], rx[1:nAtom], ry[1:nAtom], rz[1:nAtom])
    for(n = 1 ; n <= nAtom ; n ++){
      vx[n] += deltaT * ax[n];
      rx[n] += deltaT * vx[n];
      vy[n] += deltaT * ay[n];
      ry[n] += deltaT * vy[n];
      //vz[n] += deltaT * az[n];
      //rz[n] += deltaT * vz[n];
    }
  }

}
