/*
 * @file   initial.cpp
 * @author onedimension <onedimension@onedimension-PC>
 * @date   Tue Sep 30 21:33:57 2014
 * 
 * @brief 初始化，申请内存，定义一些东西
 * 
 * 
 */
using namespace std;

#include <fftw3.h>
#include <cmath>
#include "initial.h"
#include "global.h"
#include "Basis.h"

void malloc_variables_init() 
{
  Anlm = new double[7 * Basis];   //q_i对应的Zenike多项式的系数，
  //  Vnlm = new double[5 * (2 * M - 1) * L * N];  //？
  grad_Fb = new double[7 * Basis];
  grad_Energy = new double[7 * Basis];
  Qijk = new double[7 * Point];   //Q的展开式系数
  fbulk = new double[innerPoint];  
  grad_fbulk = new double[7 * innerPoint];
 
  grad_fs1 = new double[7 * I * K]();
  grad_fs2 = new double[7 * I * K]();
  grad_Fs = new double[7 * Basis]();

  grad_fs_dr = new double[I * K];
  grad_fs_dt = new double[I * K];

  drdx = new double[innerPoint]();
  drdy = new double[innerPoint]();
  drdz = new double[innerPoint]();
  dtdx = new double[innerPoint]();
  dtdy = new double[innerPoint]();
  dtdz = new double[innerPoint]();
  dpdx = new double[innerPoint]();
  dpdy = new double[innerPoint]();
  dpdz = new double[innerPoint]();

  /**
   * 定义复合求导的运算
   * 
   */
  for (int i = 0; i < I; i++) 
  {
    for (int j = 0; j < J; j++) 
    {
      for (int k = 0; k < K; k++) 
      {
	drdx[i * J * K + j * K + k] = (cosh(realxi(p[j])) * cos(mu[i]) - 1) * cos(theta[k]) / a;
	drdy[i * J * K + j * K + k] = (cosh(realxi(p[j])) * cos(mu[i]) - 1) * sin(theta[k]) / a;
	drdz[i * J * K + j * K + k] =  - sinh(realxi(p[j])) * sin(mu[i])/ a;
	dpdx[i * J * K + j * K + k] = - sinh(realxi(p[j])) * sin(mu[i]) * cos(theta[k]) / a / (xi1 - xi0) * 2;
	dpdy[i * J * K + j * K + k] = - sinh(realxi(p[j])) * sin(mu[i]) * sin(theta[k]) / a / (xi1 - xi0) * 2;
	dpdz[i * J * K + j * K + k] = (1 - cosh(realxi(p[j])) * cos(mu[i])) / a / (xi1 - xi0) * 2;

	dtdx[i * J * K + j * K + k] = - (cosh(realxi(p[j])) - cos(mu[i]))/(a*sin(mu[i])) * sin(theta[k]);
	dtdy[i * J * K + j * K + k] = (cosh(realxi(p[j])) - cos(mu[i]))/(a*sin(mu[i])) * cos(theta[k]);
	dtdz[i * J * K + j * K + k] = 0;
      }
    }
  } 
    
  dr_qijk = new double[7 * innerPoint]();
  dt_qijk = new double[7 * innerPoint]();
  dp_qijk = new double[7 * innerPoint]();

  Q1x = new double[innerPoint]();
  Q1y = new double[innerPoint]();
  Q1z = new double[innerPoint]();
  Q2x = new double[innerPoint]();
  Q2y = new double[innerPoint]();
  Q2z = new double[innerPoint]();
  Q3x = new double[innerPoint]();
  Q3y = new double[innerPoint]();
  Q3z = new double[innerPoint]();
  Q4x = new double[innerPoint]();
  Q4y = new double[innerPoint]();
  Q4z = new double[innerPoint]();
  Q5x = new double[innerPoint]();
  Q5y = new double[innerPoint]();
  Q5z = new double[innerPoint]();

  phix = new double[innerPoint]();
  phiy = new double[innerPoint]();
  phiz = new double[innerPoint]();

  rr_x = new double[innerPoint]();
  rr_y = new double[innerPoint]();
  rr_z = new double[innerPoint]();

  QL2 = new double[innerPoint]();

  grad_QL2_dr = new double[7 * innerPoint]();
  grad_QL2_dt = new double[7 * innerPoint]();
  grad_QL2_dp = new double[7 * innerPoint]();

  grad_QL1_dr = new double[7 * innerPoint]();
  grad_QL1_dt = new double[7 * innerPoint]();
  grad_QL1_dp = new double[7 * innerPoint]();

  grad_FeL1 = new double[7 * Basis]();
  grad_FeL1_r = new double[7 * Basis]();
  grad_FeL1_t = new double[7 * Basis]();
  grad_FeL1_p = new double[7 * Basis]();

  grad_FeL2 = new double[7 * Basis]();
  grad_FeL2_r = new double[7 * Basis]();
  grad_FeL2_t = new double[7 * Basis]();
  grad_FeL2_p = new double[7 * Basis]();


}

void var_destroy() 
{
  delete [] Anlm;
  delete [] Vnlm;
  delete [] grad_Fb;
  delete [] grad_Energy;
  delete [] Qijk;
  delete [] fbulk;
  delete [] grad_fbulk;

  delete [] drdx;
  delete [] drdy;
  delete [] drdz;
  delete [] dtdx;
  delete [] dtdy;
  delete [] dtdz;
  delete [] dpdx;
  delete [] dpdy;
  delete [] dpdz;

  delete [] dr_qijk;
  delete [] dt_qijk;
  delete [] dp_qijk;

  delete [] Q1x;
  delete [] Q1z;
  delete [] Q2x;
  delete [] Q2y;
  delete [] Q3x;
  delete [] Q3z;
  delete [] Q4y;
  delete [] Q4z;
  delete [] Q5y;
  delete [] Q5z;

  delete [] QL2;
  delete [] grad_QL2_dr;
  delete [] grad_QL2_dt;
  delete [] grad_QL2_dp;


  delete [] grad_QL1_dr;
  delete [] grad_QL1_dt;
  delete [] grad_QL1_dp;

  delete [] grad_FeL2;
  delete [] grad_FeL2_r;
  delete [] grad_FeL2_t;
  delete [] grad_FeL2_p;

  delete [] grad_FeL1;
  delete [] grad_FeL1_r;
  delete [] grad_FeL1_t;
  delete [] grad_FeL1_p;

  delete[] grad_Fs;
  delete[] grad_fs1;
  delete[] grad_fs2;
}

/**
void tranA2V(double *A,double *V) {
  int i,n,l,m;
  for (i = 0; i < 5; i++)
    for (m = 1 - M; m <= M - 1; m++)
      for (l = abs(m); l < L; l++)
	for (n = l; n < N; n += 2)
	  V[inx(i,n,l,m)] = (*A++);
}
*/


