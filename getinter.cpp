#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#include "sphere.h"
#include "zernike.h"
#include "initial.h"
#include "calculate_Fb.h"
#include "calculate_Fe.h"

#include "EigInfoTwo.h"
#include "in_and_out.h"
#include "point.h"
#include "rotation.h"

/*********************************************************************/

int main(int argc,char *argv[]) {
  int i,n,read_ok;
  int dis=200;
  int NP=1192;
  double x[NP],y[NP],z[NP],Q_inter[5*NP],r[NP],t[NP],p[NP];
  double tx,ty,tz,ra,th,ph;
  char fname[200];
  FILE *fp;
  
  //eta = 50;
  //landau_t = -8;
  //R1=4;
  //  Ra=1.95;
  // R2=R1*Ra;

  zernike_init_fijk(32,32,16,dis-1,dis+1,dis);
  Anlm = new double[5*Basis]();
  Qijk = new double[5*Point]();
  iput(Anlm,landau_t,R1,Ra,L21,eta,N,L,M,argv[1],1,1);
  for (i=0;i<5;i++)
    calc_fijk(Anlm + i*Basis,Qijk + i*Point);

  for (i=0;i<5*Point;i++)
    Qijk[i] = Qijk[i]/Qscale;

  if ((fp = fopen("POINT/R1_1192.txt","r")) == NULL) printf("Open Point Error!\n");
  for (i=0;i<NP;i++){
    read_ok = fscanf(fp,"%lf %lf %lf\n",&x[i],&y[i],&z[i]);
    //printf("r=%f\n",*(radius+i));
    }
  fclose(fp);
  /*
   for (i=0;i<NP;i++){
    printf("x=%lf y=%lf z=%lf\n",x[i],y[i],z[i]);
    //printf("r=%f\n",*(radius+i));
    }
  */
  for (i=0;i<NP;i++){
    tx=x[i];
    ty=y[i];
    tz=z[i];
    //printf("NP=%d x=%f y=%f z=%f\n",NP,x[i],y[i],z[i]);
    ra = sqrt(tx*tx + ty*ty + tz*tz);
    th = acos(tz/ra);
    if (fabs(tx) < 1e-14 && ty >= 0)
      ph = PI/2;
    else if (fabs(tx) < 1e-14 && ty < 0)
      ph = PI*3/2;
    else if (tx > 0 && ty >= 0)
      ph = atan(ty/tx);
    else if (tx > 0 && ty < 0)
      ph = 2*PI + atan(ty/tx);
    else if (tx < 0 && ty >= 0)
      ph = PI + atan(ty/tx);
    else if (tx < 0 && ty < 0)
      ph = PI + atan(ty/tx);
    for (n=0;n<5;n++)
      Q_inter[i + n*NP] = interpolation(Qijk + n*(dis+1)*(dis+1)*dis,dis,dis+1,dis,ra,th,ph,dis);
    r[i]=ra;
    t[i]=th;
    p[i]=ph;
    //printf("i=%d\n",i);
  }

  sprintf(fname,"%s/Drawface/%s%s_%s_t_%.2f_R_%.2f_dis_%d_NP_%d_eta_%.1f_Ra_%.3f_%s_inter.txt",DIR,func_type,boundary_in,boundary_out,landau_t,R1,dis,NP,eta,Ra,argv[1]);
	  
  fp = fopen(fname,"w");

  for (i=0;i<NP;i++) {
    fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",(2*(R2-R1)*r[i]+(2*R1-R2))*sin(t[i])*cos(p[i]),(2*(R2-R1)*r[i]+(2*R1-R2))*sin(t[i])*sin(p[i]),(2*(R2-R1)*r[i]+(2*R1-R2))*cos(t[i]),1.0/3+Q_inter[i],Q_inter[i + NP],Q_inter[i + 2*NP],1.0/3+Q_inter[i + 3*NP],Q_inter[i + 4*NP],1.0/3-Q_inter[i]-Q_inter[i + 3*NP]);
  }

  fclose(fp);
  zer_destroy();
  delete [] Anlm;
  delete [] Qijk; 
  return 0;
}
