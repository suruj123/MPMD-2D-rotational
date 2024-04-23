#include "mpi.h"
//#include <openacc.h>
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

  void Init(){

 char dummy[128];
  FILE *fp;
  fp = fopen("input-data","r");
  fscanf(fp, "%s %lf", dummy, &GAMMA);
  fscanf(fp, "%s %lf", dummy, &rCut);
  fscanf(fp, "%s %lf", dummy, &kappa);
  fscanf(fp, "%s %lf", dummy, &deltaT);
  fscanf(fp, "%s %d", dummy, &stepAvg);
  fscanf(fp, "%s %d", dummy, &stepEquil);
  fscanf(fp, "%s %d", dummy, &stepLimit);
  fscanf(fp, "%s %d", dummy, &stepDump);
  fscanf(fp, "%s %d", dummy, &stepTrajectory);
  fscanf(fp, "%s %c", dummy, &thermo);
  fscanf(fp, "%s %c", dummy, &BC);
  fscanf(fp, "%s %d", dummy, &limitCorrAv);
  fscanf(fp, "%s %d", dummy, &nBuffCorr);
  fscanf(fp, "%s %d", dummy, &nFunCorr);
  fscanf(fp, "%s %d", dummy, &nFunCorr2);
  fscanf(fp, "%s %d", dummy, &nValCorr);
  fscanf(fp, "%s %d", dummy, &stepCorr);
  fscanf(fp, "%s %d", dummy, &limitAcfAv);
  fscanf(fp, "%s %d", dummy, &nBuffAcf);
  fscanf(fp, "%s %d", dummy, &nValAcf);
  fscanf(fp, "%s %d", dummy, &stepAcf);
  fscanf(fp, "%s %lf", dummy, &rangeRdf);
  fscanf(fp, "%s %d", dummy, &limitRdf);
  fscanf(fp, "%s %d", dummy, &sizeHistRdf);
  fscanf(fp, "%s %d", dummy, &stepRdf);
  fscanf(fp, "%s %d", dummy, &limitVel);
  fscanf(fp, "%s %lf", dummy, &rangeVel);
  fscanf(fp, "%s %d", dummy, &sizeHistVel);
  fscanf(fp, "%s %d", dummy, &stepVel);
  fscanf(fp, "%s %c", dummy, &cellAlgo);
  

  fclose(fp);

  FILE *fpSTATE;
  if((fpSTATE = fopen("../STATE","r"))==NULL)
    fprintf(fpresult,"Could not open ../STATE file\n");
  fscanf(fpSTATE, "%s %lf", dummy, &timeNow);
  fscanf(fpSTATE, "%s %d", dummy, &nAtom);
  fscanf(fpSTATE, "%s %lf", dummy, &region[1]);
  fscanf(fpSTATE, "%s %lf", dummy, &region[2]);
  fscanf(fpSTATE, "%s %lf", dummy, &radout);
  fscanf(fpSTATE, "%s %lf", dummy, &radin);
  //fscanf(fpSTATE, "%s %lf", dummy, &region[3]);
  density = nAtom/(region[1]*region[2]);
  cells[1] = region[1] / rCut;
  cells[2] = region[2] / rCut;
  //cells[3] = region[3] / rCut;
  cellList = (int *)malloc((nAtom + cells[1] * cells[2] + 1) * sizeof(int));
  regionH[1] = 0.5*region[1];
  regionH[2] = 0.5*region[2];
  //regionH[3] = 0.5*region[3];
  
  	
  rx = (double*)malloc( (nAtom + 1) * sizeof(double));
  ry = (double*)malloc( (nAtom + 1) * sizeof(double));
  //rz = (double*)malloc( (nAtom + 1) * sizeof(double));
  vx = (double*)malloc( (nAtom + 1) * sizeof(double));
  vy = (double*)malloc( (nAtom + 1) * sizeof(double));
  //vz = (double*)malloc( (nAtom + 1) * sizeof(double));
  ax = (double*)malloc( (nAtom + 1) * sizeof(double));
  ay = (double*)malloc( (nAtom + 1) * sizeof(double));
  //az = (double*)malloc( (nAtom + 1) * sizeof(double));
  fax = (double*)malloc( (nAtom + 1) * sizeof(double));
  fay = (double*)malloc( (nAtom + 1) * sizeof(double));
  //faz = (double*)malloc( (nAtom + 1) * sizeof(double));
  dx = (double*)malloc( (nAtom + 1) * sizeof(double));
  p = (double*)malloc( (nAtom + 1) * sizeof(double));
  //atomType = (double*)malloc( (nAtom + 1) * sizeof(double));

  
 // #pragma acc enter data copyin(region[1:], regionH[1:])
	
  int n, idx;
   
    
    for(n = 1; n <= nAtom; n ++){
    fscanf(fpSTATE, "%lf %lf %lf %lf %lf", &dx[n], &rx[n], &ry[n], &vx[n], &vy[n]);

	//printf("error=%lf\n",rx[n]);
    }
    fclose(fpSTATE);

 // #pragma acc enter data copyin(rx[1:nAtom], ry[1:nAtom]) create(ax[1:nAtom], ay[1:nAtom]) 
  
  for(n = 1; n <= nAtom; n ++){
  //printf("%lf %lf %lf %lf %lf taskid=%d\n", dx[n], rx[n], ry[n], vx[n], vy[n],taskid);
  }
}
