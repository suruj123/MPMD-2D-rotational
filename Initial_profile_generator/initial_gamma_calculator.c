/* Generates a coordinates and velocities grid from trajectory data */
/* This grid is useful for plotting quiver plots */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define	NDIM 2
#define Sqr(x) ((x)*(x))

int main(int argc, char **argv){
  printf("%d\n",argc);
  if(argc != 5){
    printf("Usage : ./grid [*.xyz] [# of frames] [ngrids_X] [ngrids_Y]\n");
    exit(1);
  }

  FILE *fp;
  fp = fopen(argv[1], "r");
  int Frames = atoi(argv[2]);
  char dummy[25];
  int frn;
  int nAtom, n, atomType, c;
  double timeNow, region[NDIM+1], regionH[NDIM+1]; 
  double *rx, *ry, *vx, *vy;
  double pi = 4.0*atan(1.0);

  int ngrids_X = atoi(argv[3]);
  int ngrids_Y = atoi(argv[4]);
 
  int NHIST = 4;
  double **Grid;
  Grid = (double **)malloc((NHIST+1)*sizeof(double*));
  for(n = 0 ; n <= NHIST ; n ++)
    Grid[n] = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  int j;
  for(j = 1 ; j <= NHIST ; j ++)
    for(n = 1; n <= ngrids_Y*ngrids_X; n ++)
      Grid[j][n] = 0.;

  double *streamFun;
  streamFun = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));
  
  double *omega;
  omega = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  double *Q;
  Q = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  double gridWidth[NDIM+1];
  FILE *fptime;
  fptime = fopen ("kinenergy1","w");
 
  double angle_factor = 180.0/pi;
  double radius, angle = 2.0*pi*angle_factor, theta, r, radout, radin; 
  int maxc = 1, minc = 1;
  double maxvort = 0.0, minvort = 0.0;
  
  printf("%lf\n", angle);

  // Scan the frames 
  for(frn = 1; frn <= Frames; frn ++){
    
  fscanf(fp, "%s %lf", dummy, &timeNow);
  fscanf(fp, "%s %d", dummy, &nAtom);
  fscanf(fp, "%s %lf", dummy, &region[1]);
  fscanf(fp, "%s %lf", dummy, &region[2]);
  fscanf(fp, "%s %lf", dummy, &radout);
  fscanf(fp, "%s %lf", dummy, &radin);




    /*fscanf(fp, "%d", &nAtom);
    fscanf(fp, "%s %lf %s %lf %s %lf", dummy, &timeNow, dummy, &region[1], dummy,
	&region[2]);*/

    if(frn == 1){
      radius = radout - radin;
      gridWidth[1] = radius/ngrids_X;
      gridWidth[2] = angle/ngrids_Y;
      // gridWidth[1] = region[1]/ngrids_X;
      //gridWidth[2] = region[2]/ngrids_Y;
      regionH[1] = 0.5*region[1];
      regionH[2] = 0.5*region[2];
      rx = (double *)malloc((nAtom+1)*sizeof(double));
      ry = (double *)malloc((nAtom+1)*sizeof(double));
      vx = (double *)malloc((nAtom+1)*sizeof(double));
      vy = (double *)malloc((nAtom+1)*sizeof(double));
    }
    
    for(n = 1; n <= nAtom; n ++){

      fscanf(fp, "%d %lf %lf %lf %lf", &atomType, &rx[n], &ry[n], &vx[n], &vy[n]);

      r = sqrt( Sqr(rx[n]) + Sqr(ry[n]) );

      theta = atan(ry[n]/rx[n]);

     if (rx[n] == 0. && ry[n] == 0.) theta = 0.;
     if( rx[n] == 0.&& ry[n] > 0.) theta = 0.5*pi;
     else if( rx[n] == 0.&& ry[n] < 0.) theta = 1.5*pi;
     else if( rx[n] < 0. && ry[n] == 0.) theta = pi;
     else if( rx[n] < 0. && ry[n] > 0.) theta += pi;
     else if( rx[n] < 0. && ry[n] < 0.) theta += pi;
     else if( rx[n] > 0. && ry[n] < 0.) theta += 2.0*pi;
     else if( rx[n] > 0. && ry[n] == 0.) theta = 0.0;

     theta = theta*angle_factor;
     //printf("%lf %d\n", theta, n);

//      c = (int)((ry[n]+regionH[2])/gridWidth[2]) * ngrids_X +
//	(int)((rx[n]+regionH[1])/gridWidth[1]) + 1;
	c = (int)((theta)/gridWidth[2]) * ngrids_X +
	(int)((r - radin)/gridWidth[1]) + 1;

	/*if(n==1){
	maxc = c;
	minc = c;}

	if(c>maxc)
	maxc = c;

 	if(c<minc)
	minc = c; 

	if( c==5)
	printf("c = %d %d %lf\n", c, n, r);*/

      Grid[1][c] ++;				
      Grid[3][c] += (rx[n]*vx[n] + ry[n]*vy[n])/(r);
      Grid[4][c] += (rx[n]*vy[n] - ry[n]*vx[n])/(r);
    }

//printf("%d %d\n",maxc, minc);

    // Get the average center of mass velocity
    int nx, ny;
    for(ny = 1 ; ny <= ngrids_Y ; ny ++){
      for(nx = 1 ; nx <= ngrids_X ; nx ++){
	n = (ny - 1 )*ngrids_X + nx;
	if(Grid[1][n] > 0){
  	  Grid[3][n] /= Grid[1][n];
  	  Grid[4][n] /= Grid[1][n];
	}
      }
    }
      
    // For thermal kinetic energy

    double vr, vtheta;

    for(n = 1; n <= nAtom; n ++){

      r = sqrt( Sqr(rx[n]) + Sqr(ry[n]) );

      theta = atan(ry[n]/rx[n]);

     if (rx[n] == 0. && ry[n] == 0.) theta = 0.;
     if( rx[n] == 0.&& ry[n] > 0.) theta = 0.5*pi;
     else if( rx[n] == 0.&& ry[n] < 0.) theta = 1.5*pi;
     else if( rx[n] < 0. && ry[n] == 0.) theta = pi;
     else if( rx[n] < 0. && ry[n] > 0.) theta += pi;
     else if( rx[n] < 0. && ry[n] < 0.) theta += pi;
     else if( rx[n] > 0. && ry[n] < 0.) theta += 2.0*pi;
     else if( rx[n] > 0. && ry[n] == 0.) theta = 0.0;

     theta = theta*angle_factor;

     c = (int)((theta)/gridWidth[2]) * ngrids_X +
	(int)((r - 50)/gridWidth[1]) + 1;

      vr = (rx[n]*vx[n] + ry[n]*vy[n])/(r);
      vtheta = (rx[n]*vy[n] - ry[n]*vx[n])/(r);


      Grid[2][c] += 0.5 * ( Sqr(vr - Grid[3][c]) + Sqr(vtheta - Grid[4][c]) );	
    }

    // Vorticity 
    double Urp, Urn, Utp, Utn;
    for(ny = 1 ; ny <= ngrids_Y ; ny ++){
      for(nx = 1 ; nx <= ngrids_X ; nx ++){
	if(nx == 1){
	  Utn = Grid[4][(ny-1)*ngrids_X + 2];
	  //Utp = Grid[4][(ny-1)*ngrids_X + ngrids_X];
	  Utp = 0.0;
	}else if(nx == ngrids_X){
	  //Utn = Grid[4][(ny-1)*ngrids_X + 1];
	  Utn = 0.0;
	  Utp = Grid[4][(ny-1)*ngrids_X + ngrids_X-1];
	}else{
	  Utn = Grid[4][(ny-1)*ngrids_X + nx+1];
	  Utp = Grid[4][(ny-1)*ngrids_X + nx-1];
	}

	if(ny == 1){
	  Urn = Grid[3][ngrids_X + nx];
	  Urp = Grid[3][(ngrids_X-1)*ngrids_X + nx];
	}else if(ny == ngrids_Y){
	  Urn = Grid[3][nx];
	  Urp = Grid[3][(ngrids_Y-2)*ngrids_X + nx];
	}else{
	  Urn = Grid[3][ny*ngrids_X + nx];
	  Urp = Grid[3][(ny-2)*ngrids_X + nx];
	}

	n = (ny-1)*ngrids_X + nx;
	//omega[n] = (Uyn - Uyp - Uxn + Uxp)/(2*gridWidth[1]);
	omega[n] = (Grid[4][n])/r + (Utn - Utp)/(2*gridWidth[1]) - (Urn - Urp)/(2*gridWidth[2]*r);
      }
    }

    // Okubo-Weiss parameter
    double dxUx, dxUy, dyUy, dyUx;
    for(ny = 1 ; ny <= ngrids_Y ; ny ++){
      for(nx = 1 ; nx <= ngrids_X ; nx ++){
	n = (ny-1)*ngrids_X + nx;
	if(nx == 1){
	  dxUx = (Grid[3][(ny-1)*ngrids_X + 2] - Grid[3][(ny-1)*ngrids_X + ngrids_X])/(2*gridWidth[1]);
	  dxUy = (Grid[4][(ny-1)*ngrids_X + 2] - Grid[4][(ny-1)*ngrids_X + ngrids_X])/(2*gridWidth[1]);
	}else if(nx == ngrids_X){
	  dxUx = (Grid[3][(ny-1)*ngrids_X + 1] - Grid[3][(ny-1)*ngrids_X + ngrids_X-1])/(2*gridWidth[1]);
	  dxUy = (Grid[4][(ny-1)*ngrids_X + 1] - Grid[4][(ny-1)*ngrids_X + ngrids_X-1])/(2*gridWidth[1]);
	}else{
	  dxUx = (Grid[3][(ny-1)*ngrids_X + nx+1] - Grid[3][(ny-1)*ngrids_X + nx-1])/(2*gridWidth[1]);
	  dxUy = (Grid[4][(ny-1)*ngrids_X + nx+1] - Grid[4][(ny-1)*ngrids_X + nx-1])/(2*gridWidth[1]);
	}
	if(ny == 1){
	  dyUy = (Grid[4][ngrids_X + nx] - Grid[4][(ngrids_X-1)*ngrids_X + nx])/(2*gridWidth[2]);
	  dyUx = (Grid[3][ngrids_X + nx] - Grid[3][(ngrids_X-1)*ngrids_X + nx])/(2*gridWidth[2]);
	}else if(ny == ngrids_Y){
	  dyUy = (Grid[4][nx] - Grid[4][(ngrids_Y-2)*ngrids_X + nx])/(2*gridWidth[2]);
	  dyUx = (Grid[3][nx] - Grid[3][(ngrids_Y-2)*ngrids_X + nx])/(2*gridWidth[2]);
	}else{
	  dyUy = (Grid[4][ny*ngrids_X + nx] - Grid[4][(ny-2)*ngrids_X + nx])/(2*gridWidth[2]);
	  dyUx = (Grid[3][ny*ngrids_X + nx] - Grid[3][(ny-2)*ngrids_X + nx])/(2*gridWidth[2]);
	}
	Q[n] = Sqr(dxUx - dyUy) + Sqr(dyUx + dxUy) - Sqr(omega[n]);
      }
    }

    /*int x0=1, y0=1, xp, yp;
    for(n = 1; n <= ngrids_Y*ngrids_X; n ++)
      streamFun[n] = 0.0;
    for(ny = x0 ; ny <= ngrids_Y ; ny ++){
      for(nx = y0 ; nx <= ngrids_X ; nx ++){
	n = (ny-1)*ngrids_X + nx; 
	for(yp = y0; yp <= ny; yp ++){ 
	  streamFun[n] += Grid[3][(yp-1)*ngrids_X + nx]*gridWidth[2];
	}
	for(xp = x0; xp <= nx; xp ++){ 
	  streamFun[n] -= Grid[4][(y0-1)*ngrids_X + xp]*gridWidth[1];
	}
      }
    }*/

// Get the Gamma profile
    	double ke ,ke2; 
	int countgrid;
	ke = 0.0;
	countgrid = 0;

    for(ny = 1 ; ny <= ngrids_Y ; ny ++){
      for(nx = 1 ; nx <= ngrids_X ; nx ++){
	n = (ny - 1 )*ngrids_X + nx;
	// Kinetic energy per particle
	if (Grid[1][n]>0)
	Grid[2][n] /= Grid[1][n];
	ke += Grid[2][n];
	countgrid += 1 ;
      }
    }

  ke2 = ke/countgrid;
  fprintf(fptime,"%lf  %lf  %d %d\n", ke2, 1.0/ke2 ,frn, countgrid);
  //printf("%d\n",frn);
  //printf("%lf  %lf  %d %d\n",ke2, 1.0/ke2 ,frn, countgrid);




    // Get the average kinetic energy
    	/*double ke ,ke2; 
	int countgrid;
	ke = 0.0;
	countgrid = 0;
    	
    for(ny = 1 ; ny <= ngrids_Y ; ny ++){
      for(nx = 1 ; nx <= ngrids_X ; nx ++){
	n = (ny - 1 )*ngrids_X + nx;
	// Kinetic energy per particle
	Grid[2][n] /= Grid[1][n];
	ke += Grid[2][n];
	countgrid += 1 ;
      }
    }
  ke2 = ke/countgrid;
  fprintf(fptime,"%lf  %lf  %d %d\n",ke2, 1.0/ke2 ,frn, countgrid);*/

    double binArea = region[1]*region[2]/(ngrids_X*ngrids_Y);

    // Get the stream function
    double sFirst = 0.;
    for(ny = 1 ; ny <= ngrids_Y ; ny ++){
      for(nx = 1 ; nx <= ngrids_X ; nx ++){
	n = (ny - 1 )*ngrids_X + nx;
	if(nx == 1){
	  sFirst += (Grid[1][n]/binArea) * Grid[3][n] * gridWidth[2];
	  streamFun[n] = sFirst;
	}else
	  streamFun[n] = streamFun[n-1] - (Grid[1][n]/binArea) * Grid[4][n] *gridWidth[1];

        if(Grid[1][n] == 0){
	  streamFun[n] = 0.;
	}

      }
    }

    // Now dump the data to 'stdout' for this frame
    FILE *outfp;
    char outfile[25];
    sprintf(outfile,"QUIVER_time_1_%.0lf", timeNow);
    outfp = fopen(outfile, "w");
    //double x, y;
    int nnx , nny;
    if(ngrids_X % 2 == 0){
    nnx = (int)(0.5*ngrids_X) + 1;
    nny = (int)(0.5*ngrids_Y) + 1;
    }else{
    nnx = (int)(0.5*(ngrids_X - 1)) + 1;
    nny = (int)(0.5*(ngrids_Y - 1)) + 1;
    } 
	if(frn == 1){
	printf("nnx = %d nny = %d\n",nnx,nny);
	}

    double avvtheta = 0.0, x = radin, y = 0.0, cx, cy, y1;

for(ny = 1 ; ny <= ngrids_Y ; ny ++){
  
  if(ny>1)
  y += gridWidth[2];

  y1 = y*(pi/180.0);
// printf("%lf %lf\n", y, y1);

  for(nx = 1 ; nx <= ngrids_X ; nx ++){

  if(nx>1)
  x += gridWidth[1];

  cx = x*cos(y1);
  cy = x*sin(y1);

    //x = (nx-0.5)*gridWidth[1] + 50;
 
	n = (ny - 1 )*ngrids_X + nx;
//	y = (ny-0.5)*gridWidth[2] - regionH[2];
	//if(nx == nnx){

	//vr = sqrt(Sqr(Grid[3][n]) + Sqr(Grid[4][n]));
	//avvr += vr;
	avvtheta += Grid[4][n]; 


	if(n == 1 && frn == 1){
	maxvort = omega[n];	
	minvort  = omega[n];
	printf("%d %d\n", n, frn);
	}

	if(omega[n] > maxvort)
	maxvort = omega[n];

        if(omega[n] < minvort)
	minvort = omega[n];


	fprintf(outfp, "%lf\t %lf\t %lf\t %lf\t %lf\t%lf\n" , cx, cy,
	    Grid[3][n], Grid[4][n], Grid[1][n]/(binArea),omega[n]);
       //}
	Grid[1][n] = 0;
	Grid[2][n] = 0;
	Grid[3][n] = 0;
	Grid[4][n] = 0;
      }

	x = radin;
      //avvtheta = avvtheta/ngrids_Y;
      //fprintf(outfp, "%lf\t %lf\n" , x, avvtheta);


    }
    printf("%s\n", outfile);
    fclose(outfp);



  }

  printf("maximum vorticity = %lf minimum vorticity = %lf\n", maxvort, minvort);

  fclose(fptime);
  fclose(fp);
  free(Q);
  free(omega);
  free(streamFun);
  free(rx);
  free(ry);
  free(vx);
  free(vy);
  
  
  return 0;
}
