/**
 * @file   getpoint.cpp
 * @author mofified by onedimension <onedimension@onedimension-PC>
 * @date   Thu Nov 13 16:16:24 2014
 * 
 * @brief  根据展开数系数插值得到Q(r_i, theta, phi)，输出到/Drawpoint
 * 
 * 
 */
using namespace std;

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#include "global.h"
#include "Basis.h"
#include "initial.h"
#include "calculate_Fb.h"
#include "calculate_Fel.h"

#include "EigInfoTwo.h"
#include "in_and_out.h"
#include "point.h"

/*********************************************************************/
int main(int argc,char *argv[]) 
{
  int i1, i2, i3;
  double st_Rad, st_t, st_eta, ed_Rad, ed_t, ed_eta, x, y, z;
  double q[5];
  double eg[3];
  double vec[3][3];
  double beta;
  char fname[200];
  FILE *fp;
  
  //设定物理参数
  double A = - 0.172 * 1e6, B = -2.12 * 1e6, C = 1.73*1e6, L1 = 4 * 1e-11, R0 = 5 * 1e-7, W = 10;
  double D0 = 2.2 * R0;  ///球形粒子球心的距离
  
  xi0 = acosh(D0/(2*R0));
  cout << "xi0 = " << xi0 << endl; 
  D = D0/R0;
  a = sqrt(pow(D/2.0, 2) - 1);         ///Bispherical坐标的系数
  landau_t = 27 * A * C / (B * B);
  landau_t = 2 * sqrt(6) - 16/3.0;

  cout<<"t = " << landau_t << endl;
  eta = 0.5 * 27 * C * W / (B*B*R0);  //表面能的常数
  cout <<"eta = " << eta << endl;
  ksi = 27 * C * L1/(B*B*R0*R0);
  cout << "L1 = " << ksi << endl;
  Rad = 1/sqrt(ksi);
  cout << "Rad = " << Rad << endl;

  Qscale = - 3.0 * sqrt(6)/2 * C / B;
  beig = Qscale/3.0;

  double S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;

  q_far = new double[5];
  
  q_far[0] = - S0/3.0;
  q_far[1] = 0;
  q_far[2] = 0;
  q_far[3] = - S0/3.0;
  q_far[4] = 0;

  Basis_init(64,64,3,201,201,2);

  Anlm = new double[5*Basis]();
  Qijk = new double[5*Point]();
  fbulk = new double[innerPoint]();
 
  iput(Anlm,landau_t,Rad,eta,N,L,M,argv[1],1,0); 
	
  for(int i = 0; i < 5; i++)
  {
    calc_fijk(Anlm + i * Basis, Qijk + i * Point);  //由展开式系数计算Q_{ijk}
  }

  // calc_bulkenergy_landau_de(Qijk,fbulk);

  /**
  char outputQ_filename[30];
  sprintf(outputQ_filename,"Q_I_%d_J_%d_K_%d.txt",I, J, K);
  ofstream outputQ;

  outputQ.open(outputQ_filename);
  outputQ << setprecision(10);
  for (int i = 0; i < I; i++) 
  {
    for (int j = 0; j < J; j++) 
    {
      for (int k = 0; k < K; k++) 
      {
	for (int n = 0; n < 5; n++)
	{
	  outputQ << Qijk[n * Point + i * J * K + j * K + k] << endl;
	}
      }
    }
  }
  outputQ.close();
  */
	  

  /**
  for (int i = 0; i < 5 * Point; i++)
  {
    Qijk[i] = Qijk[i]/Qscale;
  }
  */
	
  sprintf(fname,"%s/Drawpoint/%s%s_t_%.2f_R_%.2f_D_%.2f_I_%d_J_%d_K_%d_eta_%.1f_L21_%+.1f_s0_%.2f_%s_tensor.txt",DIR,func_type,boundary,landau_t,Rad, D, I,J,K,eta,L21,beig,argv[1]);
	  
  // sprintf(fname,"%s/Drawpoint/%s%s_t_%.2f_R_%.2f_I_%d_J_%d_K_%d_eta_%.1f_%s_point.txt",DIR,func_type,boundary,landau_t,Rad,I,J,K,eta,argv[1]);
  //计算特征值以及特征向量，双轴性指标
  fp = fopen(fname,"w");
  // ix = 0;

  for (int k = 0; k < K; k++) 
  {
    for (int i = 0; i < I; i += 4) 
    {
      for (int j = 0; j < J; j+=4) 
      {
	for (int n = 0; n < 5; n++)
	{
	  q[n] = (Qijk[n * Point + i * J * K + j * K + k] + q_far[n])/Qscale;
	}

	// QRforEig(q,eg,vec);
	// sort(eg,i1,i2,i3);

	//	cout << mu[i] << endl;
	      
	// beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
	x = a * sin(mu[i])/(cosh(xi0 * p[j]) - cos(mu[i])) * cos(theta[k]);
	y = a * sin(mu[i])/(cosh(xi0 * p[j]) - cos(mu[i])) * sin(theta[k]);
	z = a * sinh(xi0 * p[j])/(cosh(xi0 * p[j]) - cos(mu[i]));
	fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e \n",x,y,z,q[0],q[1],q[2],q[3],q[4]);
	// ix++;
      }
    }
  }

  /**  
  for (int i = 0; i < I; i++) 
  {
    for (int j = 0; j < J; j++) 
    {
      for (int n = 0; n < 5; n++)
      {
	q[n] = Qijk[n * Point + i * J * K + j * K];
      }
	      
      QRforEig(q,eg,vec);
      sort(eg,i1,i2,i3);
      
      beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);
      x = a * sin(mu[i])/(cosh(xi0 * p[j]) - cos(mu[i])) * cos(theta[0]);
      y = a * sin(mu[i])/(cosh(xi0 * p[j]) - cos(mu[i])) * sin(theta[0]);
      z = a * sinh(xi0 * p[j])/(cosh(xi0 * p[j]) - cos(mu[i]));
      fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",x,y,z,eg[i1],eg[i2],eg[i3],vec[0][i3],vec[1][i3],vec[2][i3],beta,fbulk[i * J * K + j * K], eg[i3] - eg[i2]);
    }
  }
  */

  fclose(fp);

  Basis_destroy();
  delete [] Anlm;
  delete [] Qijk;
  delete [] fbulk;
  
  return 0;
}
