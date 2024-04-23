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
void AccumProps0(){

    sTotEnergy = ssTotEnergy = 0.;
    sKinEnergy = ssKinEnergy = 0.;
    //sPressure = ssPressure = 0.;
    sPotEnergy = ssPotEnergy = 0.;
    //sKEConfig = ssKEConfig = 0.;
    //svirSum = 0.;

}
