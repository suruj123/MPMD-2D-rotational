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


void ComputeForce(char BC){
int i,j;
int n;
double rrCut;

double  rrRange, deltaR, normFac;
  rrRange = Sqr(rangeRdf);
  deltaR = rangeRdf / sizeHistRdf;


if(taskid == MASTER){
  if(stepCount%stepRdf == 0){
   countRdf += 1;
 }
}


if(taskid < remnant){
//printf("%d\n",taskid);
for (n = 1; n <= nAtom; n++){
    ax[n] = 0.;
    ay[n] = 0.;
   // az[n] = 0.;
  }

     
rrCut = Sqr(rCut);
uSum = 0.0 ;
  
  //double drX, drY;

//#pragma acc parallel loop gang reduction(+:uSum) present(ax[1:nAtom], ay[1:nAtom], rx[1:nAtom], ry[1:nAtom], region[1:], regionH[1:]) copy(uSum) 
  for ( i = first; i < first+npoint; i++){

  // #pragma acc loop vector
    for ( j = 1 ;  j <= nAtom; j++){
     if(j>i){
     
      //Application of periodic wraparound
double      drX = rx[i] - rx[j];
      //if (drX >= regionH[1]) drX = drX - region[1];
      //if (drX <= -regionH[1]) drX = drX + region[1];
double      drY = ry[i] - ry[j];
      //if (drY >= regionH[2]) drY = drY - region[2];
      //if (drY <= -regionH[2]) drY = drY + region[2];

/*double      drZ = rz[i] - rz[j];
      if (drZ >= regionH[3]) drZ = drZ - region[3];
      if (drZ <= -regionH[3]) drZ = drZ + region[3];*/


double      rr = Sqr(drX) + Sqr(drY);

      //Interaction is not strong enough beyond "rCut"
      if (rr < rrCut){
double        r = sqrt(rr);
double        ri = 1.0 / r;
double        uVal =ri * exp(-kappa * r);
double        fcVal = ri * uVal * (ri + kappa);
double        fx = fcVal * drX;
double        fy = fcVal * drY;
//double        fz = fcVal * drZ;
//#pragma acc atomic 
{
ax[i] += fx;
}
//#pragma acc atomic
{
ax[j] -= fx;
}
//#pragma acc atomic 
{
ay[i] += fy;
}
//#pragma acc atomic
{
ay[j] -= fy;
}
/*#pragma acc atomic 
{
az[i] += fz;
}
#pragma acc atomic
{
az[j] -= fz;
}*/

        
       uSum += uVal;
        //virSum += 0.5 * fcVal * rr;
        //rfAtom += 0.5 * drX * fcVal * drY;
      }

    }}

 

}

//#pragma acc update self(ax[1:nAtom], ay[1:nAtom])
 
MPI_Reduce(ax, fax, (nAtom+1), MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
MPI_Reduce(ay, fay, (nAtom+1), MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
//MPI_Reduce(az, faz, (nAtom+1), MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
MPI_Reduce(&uSum, &fuSum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

}

else{
//printf("else = %d\n",taskid);
for (n = 1; n <= nAtom; n++){
    ax[n] = 0.;
    ay[n] = 0.;
    //az[n] = 0.;
  }

     
rrCut = Sqr(rCut);
uSum = 0.0 ;
  
  //double drX, drY;

//#pragma acc parallel loop gang reduction(+:uSum) present(ax[1:nAtom], ay[1:nAtom], rx[1:nAtom], ry[1:nAtom], region[1:], regionH[1:]) copy(uSum) 
  for ( i = first; i < first+npoints; i++){
 
 //#pragma acc loop vector
   
    for ( j = 1 ;  j <= nAtom; j++){
     if(j>i){
     
      //Application of periodic wraparound
double      drX = rx[i] - rx[j];
      //if (drX >= regionH[1]) drX = drX - region[1];
      //if (drX <= -regionH[1]) drX = drX + region[1];
double      drY = ry[i] - ry[j];
      //if (drY >= regionH[2]) drY = drY - region[2];
      //if (drY <= -regionH[2]) drY = drY + region[2];

/*double      drZ = rz[i] - rz[j];
      if (drZ >= regionH[3]) drZ = drZ - region[3];
      if (drZ <= -regionH[3]) drZ = drZ + region[3];*/


double      rr = Sqr(drX) + Sqr(drY);

      //Interaction is not strong enough beyond "rCut"
      if (rr < rrCut){
double        r = sqrt(rr);
double        ri = 1.0 / r;
double        uVal =ri * exp(-kappa * r);
double        fcVal = ri * uVal * (ri + kappa);
double        fx = fcVal * drX;
double        fy = fcVal * drY;
//double        fz = fcVal * drZ;
//#pragma acc atomic 
{
ax[i] += fx;
}
//#pragma acc atomic
{
ax[j] -= fx;
}
//#pragma acc atomic 
{
ay[i] += fy;
}
//#pragma acc atomic
{
ay[j] -= fy;
}
/*#pragma acc atomic 
{
az[i] += fz;
}
#pragma acc atomic
{
az[j] -= fz;
}*/

        
       uSum += uVal;
        //virSum += 0.5 * fcVal * rr;
        //rfAtom += 0.5 * drX * fcVal * drY;
      }

    }}
}

//#pragma acc update self(ax[1:nAtom], ay[1:nAtom])

 
MPI_Reduce(ax, fax, (nAtom+1), MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
MPI_Reduce(ay, fay, (nAtom+1), MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
//MPI_Reduce(az, faz, (nAtom+1), MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
MPI_Reduce(&uSum, &fuSum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

  }



}
