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

void AccumProps1(){
    sTotEnergy += totEnergy;
    ssTotEnergy += Sqr(totEnergy);
    sKinEnergy += kinEnergy;
    ssKinEnergy += Sqr(kinEnergy);
    //sPressure += pressure;
    //ssPressure += Sqr(pressure);
    sPotEnergy += potEnergy;
    ssPotEnergy += Sqr(potEnergy);
    //sKEConfig += keConfig;
    //ssKEConfig += Sqr(keConfig);
    //svirSum += virSum;
}
