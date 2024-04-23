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
void Close(){

//fclose(fpresult);
 //fclose(fpxyz);
 

  //#pragma acc exit data copyout(rx[1:nAtom], ry[1:nAtom]) delete(ax[1:nAtom], ay[1:nAtom])
  //#pragma acc exit data copyout(region[1:], regionH[1:])
 
  //#pragma acc exit data copyout(rx[0:3], ry[0:3], ax[0:3], ay[0:3])
  //#pragma acc exit data copyout(region[0:3], regionH[0:3])
  free(rx);
  free(ry);
  //free(rz);
  free(vx);
  free(vy);
  //free(vz);
  free(ax);
  //free(az);
  free(ay);
  free(fax);
  free(fay);
  free(p);
  free(cellList);
  free(dx);
  //free(atomType);
 free(HistRdf);
 free(fHistRdf);

//Configurational temperature
  
  free(A);
  free(Ax);
  free(Ay);
  free(TValSum);


}
