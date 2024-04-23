#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"global.h"
#include<string.h>
void PrintSummary(){

fprintf( fpresult,"(timeNow) %lf (vSum) %lf (TE) %lf (KE) %lf (PE) %lf\n", 
    timeNow, vSum, sTotEnergy, sKinEnergy,sPotEnergy);
  fflush(fpresult);

}
