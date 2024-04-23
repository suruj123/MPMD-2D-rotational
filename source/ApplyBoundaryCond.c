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
double RandR(int *seed){
  *seed = (*seed * IMUL + IADD) & MASK;
  return (*seed * SCALE);
}




void ApplyBoundaryCond(){


int n,  vSign, plus = +1, minus = -1;

  
  int randSeed=21;
double v_o,vth,pi, p1, p2, vtheta, vthetawallout, vthetawallin;
  vth=sqrt(2.0/GAMMA);
  pi=4.0*atan(1.0);
  v_o=2.0;
  p1=0.0;
  p2=0.0;
  vthetawallout = 1.5;
  vthetawallin = -2.0;
   
    for(n = 1; n <= nAtom; n++) {
       
       p[n]=RandR(&randSeed);
	//p[n]=1.0;
    }
//#pragma acc parallel loop
  for(n = 1 ; n <= nAtom ; n ++){
    // P.B.C along x and y axis
  if(BC == 'P'){
    rx[n] -= region[1]*rint(rx[n]/region[1]);
    ry[n] -= region[2]*rint(ry[n]/region[2]);
    //rz[n] -= region[3]*rint(rz[n]/region[3]);
    }
   // R.B.C along y axis
 else if (BC == 'R'){

double theta, R;
double r = sqrt(Sqr(rx[n]) + Sqr(ry[n])), r1, vr;

double vMag = sqrt(2*(1.-1./nAtom)/GAMMA);  

//outer cylinder
  
if(r>=radout){

theta = atan(ry[n]/rx[n]);

     if (rx[n] == 0. && ry[n] == 0.) theta = 0.;
     if( rx[n] == 0.&& ry[n] > 0.) theta = 0.5*pi;
     else if( rx[n] == 0.&& ry[n] < 0.) theta = 1.5*pi;
     else if( rx[n] < 0. && ry[n] == 0.) theta = pi;
     else if( rx[n] < 0. && ry[n] > 0.) theta += pi;
     else if( rx[n] < 0. && ry[n] < 0.) theta += pi;
     else if( rx[n] > 0. && ry[n] < 0.) theta += 2.0*pi;
     else if( rx[n] > 0. && ry[n] == 0.) theta = 0.0;

R = 0.999999*radout;

rx[n] = R*cos(theta);
ry[n] = R*sin(theta);

//rx[n] = rx[n] - deltaT*vx[n];
//ry[n] = ry[n] - deltaT*vy[n];

r1 = sqrt(Sqr(rx[n]) + Sqr(ry[n]));

vr = (rx[n]*vx[n] + ry[n]*vy[n])/r;

vtheta =  (rx[n]*vy[n] - ry[n]*vx[n])/r;

vr = -vr;
//vtheta = 0.0;

vx[n] = vr*cos(theta) - vtheta*sin(theta);// - 0.1;
vy[n] = vr*sin(theta) + vtheta*cos(theta);// + 0.1;

//vx[n] = vr*cos(theta) - vthetawallout*sin(theta);// - 0.1;
//vy[n] = vr*sin(theta) + vthetawallout*cos(theta);// + 0.1;

//vx[n] = -vtheta*sin(theta);// - 0.1;
//vy[n] = vtheta*cos(theta);// + 0.1;

//vx[n] = -vx[n];
//vy[n] = -vy[n];

}


//inner cylinder

if(r<=radin){


theta = atan(ry[n]/rx[n]);

     if (rx[n] == 0. && ry[n] == 0.) theta = 0.;
     if( rx[n] == 0.&& ry[n] > 0.) theta = 0.5*pi;
     else if( rx[n] == 0.&& ry[n] < 0.) theta = 1.5*pi;
     else if( rx[n] < 0. && ry[n] == 0.) theta = pi;
     else if( rx[n] < 0. && ry[n] > 0.) theta += pi;
     else if( rx[n] < 0. && ry[n] < 0.) theta += pi;
     else if( rx[n] > 0. && ry[n] < 0.) theta += 2.0*pi;
     else if( rx[n] > 0. && ry[n] == 0.) theta = 0.0;

R = 1.000001*radin;

rx[n] = R*cos(theta);
ry[n] = R*sin(theta);

//rx[n] = rx[n] - deltaT*vx[n];
//ry[n] = ry[n] - deltaT*vy[n];

//r1 = sqrt(Sqr(rx[n]) + Sqr(ry[n]));

//vx[n] = -vx[n];
//vy[n] = -vy[n];

//vr = (x*vx[n] + y*vy[n])/r;

//vx[n] = vr*cos(theta);// - vtheta*sin(theta);
//vy[n] = vr*sin(theta);// + vtheta*cos(theta);


//vx[n] = -vthetain*sin(theta);// - 0.1;
//vy[n] = vthetain*cos(theta);// + 0.1;

vr = (rx[n]*vx[n] + ry[n]*vy[n])/r;

vtheta =  (rx[n]*vy[n] - ry[n]*vx[n])/r;

vr = -vr;
vtheta = 0.0;

//vx[n] = vr*cos(theta) - vthetawallin*sin(theta);// - 0.1;
//vy[n] = vr*sin(theta) + vthetawallin*cos(theta);// + 0.1;

vx[n] = vr*cos(theta) - vtheta*sin(theta);// - 0.1;
vy[n] = vr*sin(theta) + vtheta*cos(theta);// + 0.1;

}
     
   
} 
  
  
  else if (BC == 'M'){
  for( n = 1 ; n <= nAtom ; n ++){
 rx[n] -= region[1]*rint(rx[n]/region[1]);
 //rz[n] -= region[3]*rint(rz[n]/region[3]);
    
     if (ry[n] >= regionH[2]){
     //ry[n] = ry[n]-deltaT*vy[n];
     ry[n] = regionH[2] * 0.999999;
     vy[n]=-vth*cos(p[n]*0.50*pi);
     //if (vx[n] < 0){
     vx[n]= v_o + vth*sin(p[n]*0.50*pi)*(fabs(vx[n])/vx[n]);
      
     }else if (ry[n] <= -regionH[2]){
     //ry[n] =ry[n]-deltaT*vy[n];
     ry[n] = -0.999999*regionH[2];
     vy[n]=vth*cos(p[n]*0.50*pi);
     //if (vx[n] < 0){
     vx[n]= -v_o + vth*sin(p[n]*0.50*pi)*(fabs(vx[n])/vx[n]);
      
    }
     /*if (fabs(ry[n]) >= regionH[2]){
     ry[n] = regionH[2] * 0.999999*(fabs(ry[n])/ry[n]);
     vy[n]=-vth*cos(p[n]*0.50*pi)*(fabs(ry[n])/ry[n]); 
     vx[n] = v_o*(fabs(ry[n])/ry[n]) + vth*sin(p[n]*0.50*pi)*(fabs(vx[n])/vx[n]);
    }*/ 
   }	
  }
 }
}
