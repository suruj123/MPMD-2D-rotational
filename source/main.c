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

/*const gsl_rng_type * T;
gsl_rng * rnd;*/

void Init();
void AllocArrays();
void ComputeForce();
void LeapFrog();
void ApplyBoundaryCond();
void EvalProps();
void AccumProps1();
void AccumProps2();
void PrintSummary();
void AccumProps0();
void Close();




int main(int argc, char **argv){

//mpi environment start
time_t t1, t2;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
MPI_Status status;


   // int ngpus=acc_get_num_devices(acc_device_nvidia);
    //int devicenum=taskid%ngpus;
    //acc_set_device_num(devicenum,acc_device_nvidia);

    // Call acc_init after acc_set_device_num to avoid multiple contexts on device 0 in multi GPU systems
    //acc_init(acc_device_nvidia);


if(taskid == MASTER){
sprintf(dirprefix,"../output/");
prefix = strcat(dirprefix, argv[1]);
sprintf(result, "%s.result", prefix);
fpresult = fopen (result, "w");
sprintf (xyz, "%s.xyz", prefix);
fpxyz = fopen (xyz, "w");
sprintf (rdf, "%s.rdf", prefix);
fprdf = fopen (rdf, "w");
}

  /*gsl_rng_env_setup();
  T = gsl_rng_default;
  rnd = gsl_rng_alloc (T);
  gsl_rng_set (rnd, time (0));*/

  Init();
  AllocArrays();

//Starting the initialisation
//npoint or npoints no of particles are being sent to various cores.
//This structure will be accessed by each and every core.
//*********************************************
npoints = nAtom/numtasks;
actual = npoints*numtasks;
remnant = nAtom - actual;
npoint = npoints+1;
k = 0;
int i;
for (i = 0; i < numtasks; i++) {

if(taskid < remnant){

if (taskid == i) {
first = k + 1;
         
printf ("task=%3d  first point=%5d  npoints=%4d\n", taskid, 
                 first, npoints);
}
else k += npoint;
}
else{
if (taskid == i) {
first = k + 1 + remnant;
         
printf ("task=%3d  first point=%5d  npoints=%4d\n", taskid, 
                 first, npoints);
}
else k += npoints;
}
}
//**********************************************



if(taskid == MASTER){
//fp1 = fopen("mpi_result1.txt","w");
//fpxyz = fopen("mpi_xyz1.txt","w");
t1 = time(NULL);
}

int n;
  for (n = 1; n <= nAtom; n++){
      ax[n] = 0.;
          ay[n] = 0.;
              //az[n] = 0.;
                }

//#pragma acc update device(ax[1:nAtom], ay[1:nAtom])  

for(stepCount=1;stepCount<=stepLimit;stepCount++){


/*int n;
for (n = 1; n <= nAtom; n++){
    ax[n] = 0.;
    ay[n] = 0.;
    az[n] = 0.;
  }*/
//#pragma acc data copyin(rx[1:nAtom], ry[1:nAtom], rz[1:nAtom], ax[1:nAtom], ay[1:nAtom], az[1:nAtom], region[1:], regionH[1:]) 
{  
 ComputeForce();
}

 if(taskid == MASTER){
 timeNow=stepCount*deltaT;
 memcpy(&uSum, &fuSum, sizeof(double));
 memcpy(ax, fax, (nAtom+1)*sizeof(double));
 memcpy(ay, fay, (nAtom+1)*sizeof(double));
 //memcpy(az, faz, (nAtom+1)*sizeof(double));
 //double RandR(int *seed){
  //*seed = (*seed * IMUL + IADD) & MASK;
  //return (*seed * SCALE);
// }

 LeapFrog();
 ApplyBoundaryCond();
 EvalProps();
 AccumProps1();

 if((stepCount % stepAvg)==0){
 AccumProps2();
 

 PrintSummary();
 AccumProps0();
}
int n;
if(stepCount % stepTrajectory == 0){
fprintf(fpxyz, "%d\n", nAtom);
fprintf(fpxyz, "timeNow %lf region[1] %lf region[2] %lf\n", timeNow, region[1], region[2]);

for(n = 1 ; n <= nAtom ; n ++){
fprintf(fpxyz, "%d\t %0.16lf\t %0.16lf\t %0.16lf\t %0.16lf\n", n, rx[n], ry[n], vx[n], vy[n]);
fprintf(fpxyz, "\n");
}}

if(stepCount % stepDump == 0){
char DUMP[256];
  FILE *fpDUMP;
  sprintf (DUMP, "%s.STATE", prefix);
  fpDUMP = fopen (DUMP, "w");
 
  fprintf(fpDUMP, "timeNow %lf\n", timeNow);
  fprintf(fpDUMP, "nAtom %d\n", nAtom);
  fprintf(fpDUMP, "region[1] %lf\n", region[1]);
  fprintf(fpDUMP, "region[2] %lf\n", region[2]);
  //fprintf(fpDUMP, "region[3] %lf\n", region[3]);
  
  for(n = 1 ; n <= nAtom ; n ++)
  fprintf(fpDUMP, "%d\t %0.16lf\t %0.16lf\t %0.16lf\t %0.16lf\n", n, rx[n], ry[n], vx[n], vy[n]);
  fclose(fpDUMP);
}

first = 1; 
}

MPI_Bcast(rx, (nAtom+1), MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
MPI_Bcast(ry, (nAtom+1), MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
//MPI_Bcast(rz, (nAtom+1), MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

/*#pragma acc parallel loop present(ax[1:nAtom], ay[1:nAtom])
for(n = 1 ; n <= nAtom ; n ++){
  ax[n] = 0.0;
  ay[n] = 0.0;
  //az[n] = 0.0;
}*/

//#pragma acc update device(rx[1:nAtom], ry[1:nAtom])

} //end of time loop

if(taskid ==MASTER){
    t2 = time(NULL);
    fprintf(fpresult,"Execution time %lf secs\n", difftime(t2,t1));
  
  fclose(fpxyz);
  fclose(fpresult);
  fclose(fprdf);
  }
 Close();

MPI_Finalize();
return 0;
}
