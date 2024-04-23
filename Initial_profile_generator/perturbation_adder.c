#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define Sqr(x) ((x)*(x))

#define IMUL    314159269
#define IADD    453806245
#define MASK    2147483647
#define SCALE   0.4656612873e-9
double RandR(int *seed);
double RandR(int *seed){
  *seed = (*seed * IMUL + IADD) & MASK;
  return (*seed * SCALE);
}

int main(int argc, char **argv){

  int NX;
  int NY;
  double density;
  double GAMMA;
  double vtheta1 = 1.0, vtheta2 = 2.0, r1 = 145, r2 = 155;
  double pi = 4.0*atan(1.0);

  double timeNow, radout, radin;
  int nAtom;
  double c[3], gap[3], region[3];
  int n,nx,ny;

  char dummy[25];
  FILE *fp;
  fp = fopen("STATE1","r");
//  fscanf(fp, "%s %d", dummy, &NX);
 // fscanf(fp, "%s %d", dummy, &NY);
  //fscanf(fp, "%s %lf", dummy, &density);
  //fscanf(fp, "%s %lf", dummy, &GAMMA);

  fscanf(fp, "%s %lf", dummy, &timeNow);
  fscanf(fp, "%s %d", dummy, &nAtom);
  fscanf(fp, "%s %lf", dummy, &region[1]);
  fscanf(fp, "%s %lf", dummy, &region[2]);
  fscanf(fp, "%s %lf", dummy, &radout);
  fscanf(fp, "%s %lf", dummy, &radin);

  // SET THE COORDINATES ON A SQUARE LATTICE
  double *rx, *ry;
  rx = (double *)malloc((nAtom+1)*sizeof(double));
  ry = (double *)malloc((nAtom+1)*sizeof(double));



  radin = 50.0 , radout = 250.0;
  
  /*region[1] = NX/sqrt(density);
  region[2] = NY/sqrt(density);
  
  gap[1] = region[1]/NX;
  gap[2] = region[2]/NY;
  n = 0;
  for(ny = 1 ; ny <= NY ; ny++){
    c[2] = (ny - 0.5) * gap[2] - 0.5*region[2];
    for(nx = 1 ; nx <= NX ; nx++){
      c[1] = (nx - 0.5) * gap[1] - 0.5*region[1];
      
      radius = sqrt(c[1]*c[1] + c[2]*c[2]);
      if(radius >= radin){
      if(radius <= radout){ 
      n ++;
      rx[n] = c[1];
      ry[n] = c[2];
     }}
    }
  }*/

  //int nAtom = n;
  // SET THE VELOCITIES FOR ALL PARTICLES
  double *vx, *vy, *atomType;		
  vx = (double *)malloc((nAtom+1)*sizeof(double));
  vy = (double *)malloc((nAtom+1)*sizeof(double));
  atomType = (double *)malloc((nAtom+1)*sizeof(double));
  /*double vSum[3], ang;
  int randSeed = 21;
  vSum[1] = vSum[2] = 0.0;
  double vMag = sqrt(2*(1.-1./nAtom)/GAMMA);
  
  for(n = 1 ; n <= nAtom ; n++){
    ang = 2 * M_PI * RandR(&randSeed);
    vx[n] = vMag * cos(ang);
    vy[n] = vMag * sin(ang);
    vSum[1] += vx[n];
    vSum[2] += vy[n];
  }
  vSum[1] /= nAtom;
  vSum[2] /= nAtom;
  for(n = 1 ; n <= nAtom ; n++){
    vx[n] -= vSum[1];
    vy[n] -= vSum[2];
  }*/


  for(n=1;n<=nAtom;n++){
  fscanf(fp, "%d %lf %lf %lf %lf", &atomType[n], &rx[n], &ry[n], &vx[n], &vy[n]);
  }

  fclose(fp);

  // DUMP THE STATE ON THE 'STATE' FILE
  FILE *fpSTATE;
  fpSTATE = fopen("STATE","w");
  fprintf(fpSTATE,"timeNow 0\n");
  fprintf(fpSTATE,"nAtom %d\n", nAtom);
  fprintf(fpSTATE,"region[1] %E\n", region[1]);
  fprintf(fpSTATE,"region[2] %E\n", region[2]);
  fprintf(fpSTATE,"radout %E\n", radout);
  fprintf(fpSTATE,"radin %E\n", radin);

  double vr, vtheta, vthetaoriginal, r, vthetamax=0.0, theta;
  
  double a, b;
  
  a = (vtheta2 - vtheta1)/(r2 - r1);
  
  b = vtheta2 - a*r2;

  double V_0 = 1.0, epsilon = 0.5, del = 0.5, mode_no = 6.0;

  for(n = 1 ; n <= nAtom ; n ++){

  theta = atan(ry[n]/rx[n]);

     if (rx[n] == 0. && ry[n] == 0.) theta = 0.;
     if( rx[n] == 0.&& ry[n] > 0.) theta = 0.5*pi;
     else if( rx[n] == 0.&& ry[n] < 0.) theta = 1.5*pi;
     else if( rx[n] < 0. && ry[n] == 0.) theta = pi;
     else if( rx[n] < 0. && ry[n] > 0.) theta += pi;
     else if( rx[n] < 0. && ry[n] < 0.) theta += pi;
     else if( rx[n] > 0. && ry[n] < 0.) theta += 2.0*pi;
     else if( rx[n] > 0. && ry[n] == 0.) theta = 0.0;


  r = sqrt(Sqr(rx[n]) + Sqr(ry[n]));
  vr = (rx[n]*vx[n] + ry[n]*vy[n])/r;
  vthetaoriginal =  (rx[n]*vy[n] - ry[n]*vx[n])/r;

  /*if(r>=50 && r<r1){
  vtheta = 0.0;
  vx[n] = vr*cos(theta) - (vtheta + vthetaoriginal)*sin(theta);// - 0.1;
  vy[n] = vr*sin(theta) + (vtheta + vthetaoriginal)*cos(theta);// + 0.1;
  }
  
  if(r>=r1 && r<=r2){
  //if(r>=r1 && r<=r1+10.0){
  //vtheta = a*r + b;
  vtheta = 1.0;
  vx[n] = vr*cos(theta) - (vtheta + vthetaoriginal)*sin(theta);// - 0.1;
  vy[n] = vr*sin(theta) + (vtheta + vthetaoriginal)*cos(theta);// + 0.1;
  //fprintf(fpSTATE,"%d %E %E %E %E\n", n, rx[n], ry[n], vx[n], vy[n]);
  }

  if(r>r2){
  //vtheta = a*r + b;
  vtheta = 1.0;
  vx[n] = vr*cos(theta) - (vtheta + vthetaoriginal)*sin(theta);// - 0.1;
  vy[n] = vr*sin(theta) + (vtheta + vthetaoriginal)*cos(theta);// + 0.1;
  //fprintf(fpSTATE,"%d %E %E %E %E\n", n, rx[n], ry[n], vx[n], vy[n]);
  }*/


  vtheta = (V_0)*tanh((r-150)/epsilon) + del*cos(mode_no*theta);
  vx[n] = vr*cos(theta) - (vtheta + vthetaoriginal)*sin(theta);// - 0.1;
  vy[n] = vr*sin(theta) + (vtheta + vthetaoriginal)*cos(theta);// + 0.1;  

  fprintf(fpSTATE,"%d %E %E %E %E\n", n, rx[n], ry[n], vx[n], vy[n]);
  }
  fclose(fpSTATE);

  free(rx);
  free(ry);
  free(vx);
  free(vy);
  free(atomType);
  return 0;
}
 
