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

void AccumProps2(){
    sTotEnergy /= stepAvg;
    ssTotEnergy = sqrt(ssTotEnergy/stepAvg - Sqr(sTotEnergy));
    sKinEnergy /= stepAvg;
    ssKinEnergy = sqrt(ssKinEnergy/stepAvg - Sqr(sKinEnergy));
    //sPressure /= stepAvg;
    //ssPressure = sqrt(ssPressure/stepAvg - Sqr(sPressure));
    sPotEnergy /= stepAvg;
    ssPotEnergy = sqrt(ssPotEnergy/stepAvg - Sqr(ssPotEnergy));
    //sKEConfig /= stepAvg;
    //ssKEConfig = sqrt(ssKEConfig/stepAvg - Sqr(sKEConfig));
    //svirSum /= stepAvg;

}
