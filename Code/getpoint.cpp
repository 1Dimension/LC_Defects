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
  double A = - 0.172 * 1e6, B = -2.12 * 1e6, C = 1.73*1e6, L1 = 4 * 1e-11, R0 = 1 * 1e-5, W = 1000;
  // double A, B = - 0.816 * 3 * 1e6, C = 0.45 * 4 * 1e6, L1 = 6 * 1e-12, R0 = 1 * 1e-7, W = 10;                                                         
  //A = 0.9 * B * B / C / 27.0;  


  R1 = 0.9;
  double ep = 0.05;
  
  L21 = 0;

  // a = sqrt(pow(D/2.0, 2) - 1);         ///Bispherical     
  a = sqrt((pow(1 - R1*R1 - ep * ep, 2) - 4 * ep * ep * R1 * R1))/(2 * ep);
  // cout << "D =" << D << endl;
  cout << "a =" << a << endl;
  xi0 = asinh(a);
  xi1 = asinh (a/R1); //// 
  // xi1 = 0;
  
  cout << "Xi1 = " << xi1 << endl;
  // cout << "D =" << a * cosh(xi1)/sinh(xi1) + D/2 - 1 - R1<< endl;
  cout << "xi0 = " << xi0 << endl; 

  cout << a * cosh(xi1)/sinh(xi1) << endl;
  cout << a * cosh(xi0)/sinh(xi0) << endl;

  landau_t = 27 * A * C / (B * B);

  // landau_t = 2 * sqrt(6) - 16/3.0;
  cout<<"t = " << landau_t << endl;


  double zmax = a * cosh(xi0)/sinh(xi0) + 1.0;
  double zmin = a * cosh(xi0)/sinh(xi0) - 1.0;

  cout << zmin <<" " << zmax << endl;

  
  /**
  eta = 0.5 * 27 * C * W / (B*B*R0);  //表面能的常数
  // eta  = 1000;
  cout <<"eta = " << ea << endl;
  ksi = 27 * C * L1/(B*B*R0*R0);
  cout << "L1 = " << ksi << endl;
  Rad = 1/sqrt(ksi);
  cout << "Rad = " << Rad << endl;
  **/
  
  W = 1;
  
  Rad = 50;
  ksi = 1.0/Rad/Rad;
  // eta = 1000;
  cout << "Rad = " << Rad << endl;
  R0 = sqrt(27*C*L1/ksi/B/B);
  cout << "R0 = " << R0 << "m" << endl;
  eta = 0.5 * 27 * C * W / (B*B*R0);


  eta = 1;
 

  Qscale = - 3.0 * sqrt(6)/2.0 * C / B;
  beig = Qscale/3.0;

  double S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;

  S0 = S0/Qscale;  //Q = S0(nn - 1/3 I)


  Basis_init(128, 16, 64, 256, 33, 128);
  malloc_variables_init();
  
  // Anlm = new double[5*Basis]();
  // Qijk = new double[5*Point]();
  // fbulk = new double[innerPoint]();

  //  beig = Qscale/3;	  
  iput(Anlm,landau_t,Rad,eta,N,L,M,argv[1],1, 0); 
	
  for(int i = 0; i < 7; i++)
  {
    calc_fijk(Anlm + i * Basis, Qijk + i * Point);  //由展开式系数计算Q_{ijk}
    calc_dr_fijk( Anlm + i * Basis, dr_qijk + i * innerPoint);
    calc_dt_fijk(Anlm  + i * Basis, dt_qijk + i * innerPoint);
    calc_dp_fijk(Anlm + i * Basis, dp_qijk + i * innerPoint);
  }

  cout << "Done" << endl;
	  

  // double Fb = calc_energy_Fb(Anlm, Qijk, fbulk);
  //   double Fel = calc_F_el(Anlm);
  //cout << Fb <<' ' << Fel << endl;

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
    Qijk[i] = Qijk[i] /Qscale;  //无量纲化前的Q
  }
  */
  
 for (int i = 0; i < I; i ++) 
  {
    for (int j = 0; j < J; j ++) 
    {
      for (int k = 0; k < K; k++) 
      {
	for (int n = 0; n < 5; n++)
	{
	  Qijk[n * Point + i * J * K + j * K + k] = (Qijk[n * Point + i * J * K + j * K + k])/Qscale;
	}
      }
    }
  }
	
  sprintf(fname,"%s/Drawpoint/%s%s_t_%.2f_R_%.2f_D_%.2f_I_%d_J_%d_K_%d_eta_%.1f_L21_%+.1f_s0_%.2f_%s_point.txt",DIR,func_type,boundary,landau_t,Rad, ep, I,J,K,eta,L21,beig,argv[1]);


  // sprintf(fname,"%s/Drawpoint/%s%s_t_%.2f_R_%.2f_I_%d_J_%d_K_%d_eta_%.1f_%s_point.txt",DIR,func_type,boundary,landau_t,Rad,I,J,K,eta,argv[1]);
  //计算特征值以及特征向量，双轴性指标
  fp = fopen(fname,"w");
  // ix = 0;

  double phix, phiy, phiz;

  for (int k = 0; k < K; k++) 
  {
    for (int i = 0; i < I; i += 1) 
    {
      for (int j = 0; j < J; j += 1) 
      {
	for (int n = 0; n < 5; n++)
	{
	  q[n] = Qijk[n * Point + i * J * K + j * K + k];
	}

	double Qzz = -q[0]-q[3];	
	QRforEig(q,eg,vec);
	sort(eg,i1,i2,i3);
	      
	beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);

	int ix = i * J * K + j * K + k;

	phix = dr_qijk[ix + 5 * innerPoint] * drdx[ix] + dt_qijk[ix + 5 * innerPoint]*dtdx[ix] + dp_qijk[ix + 5 * innerPoint]*dpdx[ix];	
        phiy = dr_qijk[ix + 5 * innerPoint] * drdy[ix] + dt_qijk[ix + 5 * innerPoint]*dtdy[ix] + dp_qijk[ix + 5 * innerPoint]*dpdy[ix];
	phiz = dr_qijk[ix + 5 * innerPoint] * drdz[ix] + dt_qijk[ix + 5 * innerPoint]*dtdz[ix] + dp_qijk[ix + 5 * innerPoint]*dpdz[ix];
	
	/**
	if(mu[i] == 0)
	{
	  x = 0;
	  y = 0;
	  z = a * sinh(xi0 * p[j])/(cosh(xi0 * p[j]) - cos(mu[i])); 
	}
	else if(p[j] == 0)
	{
	  x = a * sin(mu[i])/(cosh(xi0 * p[j]) - cos(mu[i])) * cos(theta[k]);
	  y = a * sin(mu[i])/(cosh(xi0 * p[j]) - cos(mu[i])) * sin(theta[k]);
	  z = 0;
	}
	else 
	{
	}
	*/
	x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
	y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
	z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i]));
	
	fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",x,y,z,phix,phiy,phiz,vec[0][i3],vec[1][i3],vec[2][i3],beta,Qijk[6 * Point + i * J * K + j * K + k], Qijk[5 * Point + i * J * K + j * K + k], Qzz*Qzz,q[0]);
	// fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",x,y,z,eg[i1],eg[i2],eg[i3],vec[0][i3],vec[1][i3],vec[2][i3],beta,fbulk[i * J * K + j * K + k], Qijk[5 * Point + i * J * K + j * K + k], Qzz*Qzz,q[0]);
	// ix++;
      }
    }
  }
 
  for (int i = 0; i < I; i++) 
  {
    for (int j = 0; j < J; j++) 
    {
      for (int n = 0; n < 5; n++)
      {
	q[n] = Qijk[n * Point + i * J * K + j * K];
      }
	   
      double Qzz = -q[0]-q[3];   
      QRforEig(q,eg,vec);
      sort(eg,i1,i2,i3);
      
      beta = 1.0 - 6.0*pow(eg[i1]*eg[i1]*eg[i1] + eg[i2]*eg[i2]*eg[i2] + eg[i3]*eg[i3]*eg[i3],2)/pow(eg[i1]*eg[i1] + eg[i2]*eg[i2] + eg[i3]*eg[i3],3);

      int ix = i * J * K + j * K;

      phix = dr_qijk[ix + 5 * innerPoint] * drdx[ix] + dt_qijk[ix + 5 * innerPoint]*dtdx[ix] + dp_qijk[ix + 5 * innerPoint]*dpdx[ix];
      phiy = dr_qijk[ix + 5 * innerPoint] * drdy[ix] + dt_qijk[ix + 5 * innerPoint]*dtdy[ix] + dp_qijk[ix + 5 * innerPoint]*dpdy[ix];
      phiz = dr_qijk[ix + 5 * innerPoint] * drdz[ix] + dt_qijk[ix + 5 * innerPoint]*dtdz[ix] + dp_qijk[ix + 5 * innerPoint]*dpdz[ix];
      
      x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[0]);
      y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[0]);
      z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i]));
      //  fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",x,y,z,eg[i1],eg[i2],eg[i3],vec[0][i3],vec[1][i3],vec[2][i3],beta,fbulk[i * J * K + j * K] * Jacobi[i * J + j], Qijk[5 * Point + i * J * K + j * K], Qzz*Qzz, q[0]);
      	fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",x,y,z,phix,phiy,phiz,vec[0][i3],vec[1][i3],vec[2][i3],beta,Qijk[6 * Point + i * J * K + j * K], Qijk[5 * Point + i * J * K + j * K], Qzz*Qzz,q[0]);
    }
  }

  fclose(fp);

  Basis_destroy();
  var_destroy();
  
  return 0;
}
