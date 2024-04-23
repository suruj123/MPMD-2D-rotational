/*-----------------------------------------------------------------------------

    MPMD-v2.0 : MULTI POTENTIAL MOLECULAR DYNAMICS-version 2.0 
    A parallel classical molecular dynamics code
    Copyright (C) 2018  Harish Charan, charan.harish@gmail.com

    This program is free software but a proper permission must be taken: you can 
    redistribute it and/or modify it under the terms of the GNU General Public License as 
    published by the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
    more details.

    You should have received a copy of the GNU General Public License along 
    with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    See the README file in the top-level MPMD-v2.0 directory.

-----------------------------------------------------------------------------*/



#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"global.h"
#include<string.h>
#include <omp.h>
#define IMUL    314159269
#define IADD    453806245
#define MASK    2147483647
#define SCALE   0.4656612873e-9
double RandR(int *seed){
  *seed = (*seed * IMUL + IADD) & MASK;
  return (*seed * SCALE);
}


int main(int argc, char **argv){



    sprintf(dirprefix,"../output/");
    prefix = strcat(dirprefix, argv[1]);
    sprintf(result, "%s.result", prefix);
    fpresult = fopen (result, "w");
    sprintf (xyz, "%s.xyz", prefix);
    fpxyz = fopen (xyz, "w");

    
    //FILE*fpresult;
    //fpresult = fopen ("result", "w");
  //time_t t1, t2;
  //time_t t1, t2;
  //t1 = time(NULL);

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
  density = nAtom/(region[1]*region[2]);
  cells[1] = region[1] / rCut;
  cells[2] = region[2] / rCut;
  cellList = (int *)malloc((nAtom + cells[1] * cells[2] + 1) * sizeof(int));
  regionH[1] = 0.5*region[1];
  regionH[2] = 0.5*region[2];
  
  	/*GAMMA=50;
	kappa=0.5;
	stepEquil=50000;
	stepAvg	=1000;
	//thermo=C;	
	//BC=R;
	rCut=6;
  	stepLimit=100000;
        deltaT=0.01;
        //nAtom=100;
        
        
        //regionH[1]=region[1]*0.5;
	//regionH[2]=region[2]*0.5;*/
  //density = nAtom/(region[1]*region[2]);
  /*row=0.318;
  region[1]=sqrt((nAtom/row));
  region[2]=sqrt((nAtom/row));
  cells[1] = region[1] /rCut;
  cells[2] = region[2] /rCut;
  cellList = (int *)malloc((nAtom + cells[1] * cells[2] + 1) * sizeof(int));
  regionH[1] = 0.5*region[1];
  regionH[2] = 0.5*region[2];*/
	//printf("%lf %lf %lf\n",region[1],rCut,deltaT);
   //FILE *fpSTATE;
    //fpSTATE = fopen("STATE","r");
  rx = (real*)malloc( (nAtom + 1) * sizeof(real));
  ry = (real*)malloc( (nAtom + 1) * sizeof(real));
  vx = (real*)malloc( (nAtom + 1) * sizeof(real));
  vy = (real*)malloc( (nAtom + 1) * sizeof(real));
  ax = (real*)malloc( (nAtom + 1) * sizeof(real));
  ay = (real*)malloc( (nAtom + 1) * sizeof(real));
  dx = (real*)malloc( (nAtom + 1) * sizeof(real));
  p = (real*)malloc( (nAtom + 1) * sizeof(real));

  //#pragma acc enter data copyin(rx[0:3], ry[0:3],ax[0:3], ay[0:3])
  //#pragma acc enter data copyin(region[0:3], regionH[0:3])
	
  int n, idx;
   
    
    for(n = 1; n <= nAtom; n ++){
    fscanf(fpSTATE, "%lf %lf %lf %lf %lf", &dx[n], &rx[n], &ry[n], &vx[n], &vy[n]);

	//printf("error=%lf\n",rx[n]);
    }
    fclose(fpSTATE);

  
  for(n = 1; n <= nAtom; n ++){
  //printf("%lf %lf %lf %lf %lf\n", dx[n], rx[n], ry[n], vx[n], vy[n]);
  }
  for(n=1;n<=nAtom;n++){
	//printf("%d\n",n);


	}
 //allocarrays
// Configurational Temperature 
     A = (double *)malloc((nAtom+1)*sizeof(double));
     Ax = (double *)malloc((nAtom+1)*sizeof(double));
     Ay = (double *)malloc((nAtom+1)*sizeof(double));
     TValSum = (double *)malloc((nAtom+1)*sizeof(double)); 




//accumprops(0)
sTotEnergy = ssTotEnergy = 0.;
    sKinEnergy = ssKinEnergy = 0.;
    sPressure = ssPressure = 0.;
    sPotEnergy = ssPotEnergy = 0.;
    sKEConfig = ssKEConfig = 0.;
    svirSum = 0.;

//start time

//FILE *fpresult;
//fpresult=fopen("result","w");
double st = omp_get_wtime();

for(stepCount=1;stepCount<=stepLimit;stepCount++){
timeNow=stepCount*deltaT;

//printf("iteration%d\n",stepCount);




//computeforcecells
 if(cellAlgo == 'C'){
 double dr[NDIM+1], invWid[NDIM+1], shift[NDIM+1], f, fcVal, rr, rrCut, ri, r, uVal;
  int c, I, J, m1, m1X, m1Y, m2, m2X, m2Y, offset;
  int iofX[] = {0, 0, 1, 1, 0, -1, -1, -1, 0, 1},
      iofY[] = {0, 0, 0, 1 ,1, 1, 0, -1, -1, -1};
  
  
  rrCut = Sqr(rCut);
  invWid[1] = cells[1]/region[1];
  invWid[2] = cells[2]/region[2];
	
 
 
  for(n = nAtom+1; n <= nAtom+cells[1]*cells[2] ; n++)
    cellList[n] = 0;
  
  for(n = 1 ; n <= nAtom ; n ++){
    c = ((int)((ry[n] + regionH[2])*invWid[2]))*cells[1] + (int)((rx[n]+regionH[1])*invWid[1]) + nAtom+ 1;
    cellList[n] = cellList[c];
    cellList[c] = n;
  }
  
  for(n = 1 ; n <= nAtom ; n ++){
    ax[n] = 0.;
    ay[n] = 0.;
  }
  
  uSum = 0.0 ;
  virSum = 0.0;
  rfAtom = 0.0;

  //int start = 1 + rank*(cells[2]/size);
  //int end = (rank+1)*(cells[2]/size);

  for(m1Y = 1 ; m1Y <= cells[2] ; m1Y ++){
    for(m1X = 1 ; m1X <= cells[1] ; m1X ++){
      m1 = (m1Y-1) * cells[1] + m1X + nAtom;
      for(offset = 1 ; offset <= 9 ; offset ++){
	m2X = m1X + iofX[offset]; shift[1] = 0.;
	if(m2X > cells[1]){
	  m2X = 1; shift[1] = region[1];
	}else if(m2X == 0){
	  m2X = cells[1]; shift[1] = -region[1];
	}
	m2Y = m1Y + iofY[offset]; shift[2] = 0.;
	if(m2Y > cells[2]){
	  m2Y = 1; shift[2] = region[2];
	}else if(m2Y == 0){
	  m2Y = cells[2]; shift[2] = -region[2];
	}
	m2 = (m2Y-1)*cells[1] + m2X + nAtom;
	I = cellList[m1];
	while(I > 0){
	  J = cellList[m2];
	  while(J > 0){
	    if(m1 == m2 && J != I){
	      dr[1] = rx[I] - rx[J] - shift[1];
	      dr[2] = ry[I] - ry[J] - shift[2];
	      rr = Sqr(dr[1]) + Sqr(dr[2]);
	      if(rr < rrCut){
		r = sqrt(rr);
		ri = 1/r;
		uVal = ri*exp(-kappa*r);
		fcVal = ri*uVal*(ri+kappa);
		f = fcVal * dr[1];
		ax[I] += f;
		f = fcVal * dr[2];
		ay[I] += f;
		uSum +=  0.5 * uVal;
		virSum += 0.5 * fcVal * rr;
		rfAtom += 0.5 * dr[1] * fcVal * dr[2];
	      }
	    }else if(m1 != m2){
	      dr[1] = rx[I] - rx[J] - shift[1];
	      dr[2] = ry[I] - ry[J] - shift[2];
	      rr = Sqr(dr[1]) + Sqr(dr[2]);
	      if(rr < rrCut){
		r = sqrt(rr);
		ri = 1/r;
		uVal = ri*exp(-kappa*r);
		fcVal = ri*uVal*(ri+kappa);
		f = fcVal * dr[1];
		ax[I] += f;
		f = fcVal * dr[2];
		ay[I] += f;
		uSum +=  0.5 * uVal;
		virSum += 0.5 * fcVal * rr;
		rfAtom += 0.5 * dr[1] * fcVal * dr[2];
	      }
	    }
       	    J = cellList[J];
	  }
	  I = cellList[I];
	}
      }
    }
  }
}else if(cellAlgo == 'A'){
double dr[NDIM+1], invWid[NDIM+1], shift[NDIM+1],rrCut,f;
double sumfx,sumfy,sum3;
double virSum1,virSum2;
double rfAtom1,rfAtom2;
//int i,j;

for (n = 1; n <= nAtom; n++){
    ax[n] = 0.;
    ay[n] = 0.;
  }

     
rrCut = Sqr(rCut);
uSum = 0.0 ;
  
  //double drX, drY;
//1
//#pragma acc parallel loop gang reduction(+:uSum)

//2
//#pragma acc parallel loop collapse(2) reduction(+:uSum) //copyin(rx[0:nAtom],ry[:nAtom],region[0:2],regionH[0:2], dr[0:2]) copyout(ax[0:nAtom], ay[0:nAtom]) reduction(+:uSum) reduction(+:sumfx) reduction(+:sumfy) //reduction(+:sum3) 
  

  for (int i = 1; i <= nAtom; i++){
//1 
//#pragma acc loop independent

//or 1
//#pragma acc loop vector
    //i = j is not allowed in the interaction computation
    for (int j =  1; j <= nAtom; j++){
     if(j >i){
     //continue
      //Application of periodic wraparound
double      drX = rx[i] - rx[j];
      if (drX > regionH[1]) drX = drX - region[1];
      if (drX < -regionH[1]) drX = drX + region[1];
double      drY = ry[i] - ry[j];
      if (drY > regionH[2]) drY = drY - region[2];
      if (drY < -regionH[2]) drY = drY + region[2];

double      rr = Sqr(drX) + Sqr(drY);

      //Interaction is not strong enough beyond "rCut"
      if (rr < rrCut){
double        r = sqrt(rr);
double        ri = 1. / r;
double        uVal =ri * exp(-kappa * r);
double        fcVal = ri * uVal * (ri + kappa);
double        fx = fcVal * drX;
double        fy = fcVal * drY;
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

        
        uSum += uVal;
        //virSum += 0.5 * fcVal * rr;
        //rfAtom += 0.5 * drX * fcVal * drY;
      }
    }}
  } 
}






//leapfrogstep
if(stepCount <= stepEquil){
    double gSum, varS, massS;
    temperature = 1./GAMMA;
    
   if(stepCount == 1) varS = 0.;
    double A, S1, S2, T;
    int n;
    S1 = 0.; S2 = 0.; gSum = 0.; massS = 0.1; 

    vvSum = 0.;
    double halfdt = 0.5*deltaT;
    for (n = 1; n <= nAtom; n++){
      T = vx[n] + halfdt * ax[n];
      S1 += T * ax[n];
      S2 += Sqr(T);
      
      T = vy[n] + halfdt * ay[n];
      S1 += T * ay[n];
      S2 += Sqr(T);
     vvSum += (Sqr(vx[n]) +  Sqr(vy[n]));
    }

    A = -S1 / S2;
    S2 = vvSum;
   
    double C = 1 + A*deltaT ;
    double D = deltaT * (1 + 0.5 * A * deltaT);
    
    int i,j;
    real dr[NDIM+1], r, rr, ri, rrCut;
    double vv;

    double uVal, uSum, fcVal, f, AA, AASum;
    double TVal, T_Config;

    double deno, VVSum;	 
    double PE;
    PE=0.;
    deno = 0.;
    VVSum = 0.;
    keConfig = 0.;
    AASum = 0.;
   
  for(n=1;n<=nAtom; n++)
     TValSum[n] = 0.;

  rrCut = Sqr(rCut);

/*****Calculating Configarational temperature*****/
if(thermo == 'C'){
for(i = 1 ; i <= nAtom; i ++){
    for(j = i+1 ; j <= nAtom ; j ++){
      dr[1] = rx[i] - rx[j];
      if(fabs(dr[1]) > regionH[1])
	dr[1] -= SignR(region[1], dr[1]);

     dr[2] = ry[i] - ry[j];
      if(fabs(dr[2]) > regionH[2])
	dr[2] -= SignR(region[2], dr[2]);
  
     rr = Sqr(dr[1]) + Sqr(dr[2]);
     if(rr < rrCut ){
     r = sqrt(rr);
     ri = 1/r;
     uVal = ri*exp(-kappa*r);

     TVal = (1./rr + Sqr(kappa) + kappa/r)*uVal;
     TValSum[i] += TVal;
     TValSum[j] += TVal;
   } }
     AA = Sqr(ax[i]) + Sqr(ay[i]);
     AASum += AA;
     vv = Sqr(vx[i]) + Sqr(vy[i]);
     VVSum += vv; 
     deno += TValSum[i];    
} 
     PE = uSum/nAtom; 
     keConfig = AASum/deno;

     double gSumconfig, varSconfig, massSconfig;
     if(stepCount == 1) varSconfig = 0.;
     gSumconfig = 0.; massSconfig = 2.0; 
   
     gSumconfig = (AASum/deno -temperature)/massSconfig;
     varSconfig += deltaT*gSumconfig; 

      /*****Configarational Nose-Hoover thermostat*****/
   for (n = 1; n <= nAtom; n++){
      vx[n] += deltaT * ax[n];
      rx[n] += deltaT * (vx[n]);
      vy[n] += deltaT * ay[n];
      ry[n] += deltaT * (vy[n]);
    }
      /*****Kinetic Nose-Hoover thermostat*****/
  }else if(thermo == 'N'){  
    gSum = (0.5*S2 - (nAtom + 1)*temperature)/massS;
    varS += deltaT*gSum; 
   for (n = 1; n <= nAtom; n++){
      vx[n] += deltaT * (ax[n] - varS *vx[n]);
      rx[n] += deltaT * vx[n];
      vy[n] += deltaT * (ay[n] - varS *vy[n]);
      ry[n] += deltaT * vy[n];
   }
      /*****for Gaussian thermostat*****/
 }else if(thermo == 'G'){              
      for (n = 1; n <= nAtom; n++){
      vx[n] = C * vx[n] + D * ax[n];
      rx[n] += deltaT * vx[n];
      vy[n] = C * vy[n] + D * ay[n];
      ry[n] += deltaT * vy[n];
     }
   }
 }else{
    int n;
//#pragma acc parallel loop
    for(n = 1 ; n <= nAtom ; n ++){
      vx[n] += deltaT * ax[n];
      rx[n] += deltaT * vx[n];
      vy[n] += deltaT * ay[n];
      ry[n] += deltaT * vy[n];
    }
  }




//applyboundarycondition
int  vSign;

  
    int randSeed=21;
  double v_o,vth,pi,p1,p2;
  vth=sqrt(2.0/GAMMA);
  pi=4.0*atan(1.0);
  v_o=10.0*vth;
  p1=0.0;
  p2=0.0;
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
    }
   // R.B.C along y axis
 else if (BC == 'R'){
     rx[n] -= region[1]*rint(rx[n]/region[1]); 
     vSign = 0;
     if (ry[n] >= regionH[2]){
     ry[n] = regionH[2] * 0.999999;  vSign = -1;
     } else if (ry[n] <= -regionH[2]){
     ry[n] = -0.999999*regionH[2];       vSign = 1;  
     }if(vSign) {
     if (vy[n] * vSign < 0.)
     vy[n] = -vy[n];                                             
   }
  } 
  
  
  else if (BC == 'M'){

  
 
     rx[n] -= region[1]*rint(rx[n]/region[1]); 
     vSign = 0;
     if (ry[n] >= regionH[2]){
     ry[n] = ry[n]-deltaT*vy[n]; // vSign = -1;
     vy[n]=-vth*cos((p[n]+p1)*0.50*pi+p2);
     if (vx[n] < 0){
     vx[n]=-vth*sin((p[n]+p1)*0.50*pi+p2)+v_o;
     }
     else{ 
     vx[n]=vth*sin((p[n]+p1)*0.50*pi+p2)+v_o;
     } 
     } else if (ry[n] <= -regionH[2]){
     ry[n] =ry[n]-deltaT*vy[n];      // vSign = 1;  
     vy[n]=vth*cos((p[n]+p1)*0.50*pi+p2);
     if (vx[n] < 0){
     vx[n]=-vth*sin((p[n]+p1)*0.50*pi+p2);
     }
     else{ 
     vx[n]=vth*sin((p[n]+p1)*0.50*pi+p2);
     } 
    }
    
    

  } 
  }	
  
  





//evalprops
 real v, vv;
  vSum = 0.;
  vvSum = 0.;
  //int n;
  for(n = 1 ; n <= nAtom ; n ++){
    vv = 0.;
    v = vx[n] - 0.5 * deltaT * ax[n];
    vSum += v;
    vv += Sqr(v);
    v = vy[n] - 0.5 * deltaT * ay[n];
    vSum += v;
    vv += Sqr(v);
    vvSum += vv;
  }
  kinEnergy = 0.5 * vvSum / nAtom;
  potEnergy = uSum / nAtom;
  totEnergy = kinEnergy + potEnergy;
  pressure = density * (vvSum + virSum)/(nAtom * NDIM);





//accumprops(1)
    sTotEnergy += totEnergy;
    ssTotEnergy += Sqr(totEnergy);
    sKinEnergy += kinEnergy;
    ssKinEnergy += Sqr(kinEnergy);
    sPressure += pressure;
    ssPressure += Sqr(pressure);
    sPotEnergy += potEnergy;
    ssPotEnergy += Sqr(potEnergy);
    sKEConfig += keConfig;
    ssKEConfig += Sqr(keConfig);
    svirSum += virSum;





if((stepCount % stepAvg)==0){
//accumprops(2)
    sTotEnergy /= stepAvg;
    ssTotEnergy = sqrt(ssTotEnergy/stepAvg - Sqr(sTotEnergy));
    sKinEnergy /= stepAvg;
    ssKinEnergy = sqrt(ssKinEnergy/stepAvg - Sqr(sKinEnergy));
    sPressure /= stepAvg;
    ssPressure = sqrt(ssPressure/stepAvg - Sqr(sPressure));
    sPotEnergy /= stepAvg;
    ssPotEnergy = sqrt(ssPotEnergy/stepAvg - Sqr(ssPotEnergy));
    sKEConfig /= stepAvg;
    ssKEConfig = sqrt(ssKEConfig/stepAvg - Sqr(sKEConfig));
    svirSum /= stepAvg;







//printsummary
fprintf( fpresult,"(timeNow) %lf (vSum) %lf (TE) %lf (KE) %lf (KEConfig) %lf (PE) %lf (PR) %lf (VR) %lf\n", 
    timeNow, vSum, sTotEnergy, sKinEnergy, sKEConfig, sPotEnergy, sPressure, svirSum);
  fflush(fpresult);




//accumprops(0)
    sTotEnergy = ssTotEnergy = 0.;
    sKinEnergy = ssKinEnergy = 0.;
    sPressure = ssPressure = 0.;
    sPotEnergy = ssPotEnergy = 0.;
    sKEConfig = ssKEConfig = 0.;
    svirSum = 0.;


}
 if(stepCount % stepTrajectory == 0){
fprintf(fpxyz, "%d\n", nAtom);
  fprintf(fpxyz, "timeNow %lf region[1] %lf region[2] %lf\n", timeNow, region[1], region[2]);
  for(n = 1 ; n <= nAtom ; n ++){
    fprintf(fpxyz, "%d\t %lf\t %lf\t %lf\t %lf\n", n, rx[n], ry[n], vx[n], vy[n]);
  fprintf(fpxyz, "\n");
}
}	
}
//t2 = time(NULL);
//fprintf(fpresult, "Execution time %lf secs\n", difftime(t2,t1));

double runtime = omp_get_wtime() - st;
fprintf(fpresult," total: %f s\n", runtime);
 
 fclose(fpresult);
 fclose(fpxyz);
  
  //#pragma acc exit data copyout(rx[0:3], ry[0:3], ax[0:3], ay[0:3])
  //#pragma acc exit data copyout(region[0:3], regionH[0:3])
  free(rx);
  free(ry);
  free(vx);
  free(vy);
  free(ax);
  free(ay);
  free(p);
  free(cellList);

//Configurational temperature
  
  free(A);
  free(Ax);
  free(Ay);
  free(TValSum);
return 0; 
}
