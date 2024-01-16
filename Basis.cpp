/**
 * @file   zernike.cpp
 * @author onedimension <onedimension@onedimension-PC>
 * @date   Sun Dec 14 18:07:01 2014
 * 
 * @brief  
 * 
 * 
 */
using namespace std;

#include<iostream>
#include<cstdlib>
#include<fstream>
#include <iomanip>
#include<cmath>
#include "global.h"
#include "Basis.h"

double cutoff(double xi, double mu, double R_in, double R_outer)
{
  double R2;
  R2 = a * a * (sin(mu) * sin(mu) + sinh(xi) * sinh(xi))/((cosh(xi) - cos(mu)) * (cosh(xi) - cos(mu)));
  double x1 = R_outer * R_outer - R2;
  double x2 = R2 - R_in * R_in;
  if(x1 <= 0)
  {
    //    cout <<  pow(a,3)/(pow(cosh(xi) - cos(mu), 3)) << endl;
    return 0;
  }
  else
  {
    return 1;

    /**
    double y1 =  exp(-1.0/x1);
    double y2;
    if(x2 <= 0)
    {
      y2 = 0;
    }
    else
    {
      y2 = exp(-1.0/x2);
      // cout << xi << " " << mu << endl;
    }

    return y1/(y1+y2);
    */
  }
}

void sphere_param_init(int maxn,int n,int l,int m,int i,int j,int k)
{
  MaxN = maxn;
  N = n;  ///r方向上的最大多项式次数
  L = l;  ///z方向上的多项式次数
  M = m;  ///\phi方向上的多项式次数
  I = i;  ///r方向上的节点数
  J = j;  ///z方向上的节点数
  K = k;  ///\phi方向上的节点数
 
  Basis = L * ((2 * M - 1) * N - M * (M - 1));  ///基函数的个数
  Point = I * (J + 2) * K;  /// J, J + 1为边界点，分别对应 \pm 1
  innerPoint = I * J * K;   ///内部节点数
}

void Basis_init(int pr_n,int pr_l,int pr_m,int pr_i,int pr_j,int pr_k) 
{

  cout<<setprecision(20);

  int i, j, k, n, l, m, s;
  int read_ok;
  int ix;
  double x, y;
  char path[200];
  ifstream infile;
  
  sphere_param_init(128, pr_n, pr_l, pr_m, pr_i, pr_j, pr_k);

  mu = new double[I]; /// r \in [0.1] 
  lambda = new double[I];
  p = new double[J + 2];  ///xi方向上的Guass节点, J, J+1为边界
  theta = new double[K];

  Jacobi = new double[I * J];
  surface_jacobi = new double[I * 2];

  Cutoff = new double[innerPoint];
  
 //高斯积分系数，r和theta方向
  coe_mu = new double[I];
  coe_p = new double[J];

  //Znlm=P*R*X
  Rnmr = new double[I * (M * N - M * (M - 1)/2)]; /// R_n^{(m)}(r_i)
  Plp = new double[(J + 2) * L];
  Xm = new double[K * (2 * M - 1)];

  //导数
  dRnmr = new double[I * (M * N - M * (M - 1)/2)];  ///R_nm关于r的导数
  dPlp = new double[(J + 2) * L];
  dXm = new double[K * (2 * M - 1)];
  
  Kijm = new double[I * J * (2 * M - 1)];
		    
  /// FFT in \theta
  FDI_p = new double[K/2+1]();  ///input  
  FDO_p = new double[K/2+1]();  ///output
  ///Even Function Y_{k} = X_0 + (-1)^k X_{n-1} + \sum_{j = 1}^{n-2} X_{j}cos(pi * jk/(n-1))
  FFTp_p = fftw_plan_r2r_1d(K/2+1, FDI_p, FDO_p, FFTW_REDFT00,FFTW_MEASURE);

  ///Odd Function Y_{k} = 2 \sum_{j = 0}^{n-1} X_{j}sin(pi * (j+1)(k+1)/(n+1))
  FFTq_p = fftw_plan_r2r_1d(K/2-1,FDI_p, FDO_p, FFTW_RODFT00, FFTW_MEASURE);   

  /// 
  FDI_c = new double[K+1]();
  FDO_c = new double[K+1]();
  FFTp_c = fftw_plan_r2r_1d(K+1,FDI_c,FDO_c,FFTW_REDFT00,FFTW_MEASURE);  ///Even Function
  FFTq_c = fftw_plan_r2r_1d(K-1,FDI_c,FDO_c,FFTW_RODFT00,FFTW_MEASURE);  ///Odd Function

  /*********************************************************************/

  /// 读取[-1,1]上的高斯节点 cos(\mu)
  sprintf(path,"%s/Parameter/roots_I_%d.txt",DIR,I);  ///文件名
  infile.open(path);
  if ( ! infile ) printf("Open roots_I Error!\n");
  for (i = 0; i < I; i++)
  {
    infile >> lambda[i];
  }
  infile.close();

  for(i = 0; i < I; i++)
  {
    mu[i] = acos(lambda[i]);  ///  0  < mu[i] < \pi
    // cout << mu[i] << endl;
  }
 
  ///(-1,1)上的Gauss节点
  sprintf(path,"%s/Parameter/roots_J_%d.txt",DIR,J);
  infile.open(path);
  if ( ! infile ) cout << "Open roots_J Error!" << endl;
  for (j = 0; j < J + 2; j++)
  {
    infile >> p[j];
  }
  infile.close();
  
  ///读入积分系数
  ///(0, 1)上的Gauss积分系数
  sprintf(path,"%s/Parameter/coe_mu_%d.txt",DIR,I);
  infile.open(path);
  if ( ! infile ) cout << "Open coe_mu Error!" << endl;
  for(i = 0; i < I; i++)
  {
    infile >> coe_mu[i];  
    // cout << coe_mu[i] << endl;
  }
  infile.close();
  
  sprintf(path,"%s/Parameter/coe_mu_%d.txt",DIR,J);
  infile.open(path);
  if ( ! infile ) printf("Open coe_theta Error!\n");
  for (j = 0; j < J; j++)
  {
    infile >> coe_p[j];
  }
  infile.close();

  ///  
  for (k = 0; k < K; k++)
  {
    theta[k] = 2.0 * PI * k/K;
  }

  ///从文件中读去R_n^{m}(r_i)的值存到DataRnmr中  MaxN
  /**
   * Note: 文件中的存储顺序是对每一个n，m，存所有的r_i , (I+1)*maxN*maxN
   * 
   */
  double *DataRnmr = new double[I * MaxN * MaxN];
  sprintf(path,"%s/Parameter/Rnm_N_%d_I_%d.txt", DIR, MaxN, I );
  infile.open(path);
  ix = 0;
  if( !infile ) cout << "Open Rnm Error!" << endl;
  for(n = 0; n < MaxN; n++)
  {
    for(m = 0; m < MaxN; m++)
    {
      for(i = 0; i < I; i++)
      {
	infile >> DataRnmr[ix];
	ix++;
      }
    }
  }

  infile.close();


 /**
   * 用递推关系式计算dRnmr。。。
   */
  double *DatadRnmr = new double[I * MaxN * MaxN]();

  for(i = 0; i < I; i++) 
  {
    n = 0;
    DatadRnmr[n * MaxN * I + n * I + i] = 0; ///dP_0^{0}

    n = 1;
    DatadRnmr[n * MaxN * I + n * I + i] = - sqrt(3.0/4/PI) * lambda[i];
   
    for (n = 2; n < MaxN; n++)  /// P_n^{n}
    {
      DatadRnmr[n * MaxN * I + n * I + i] = - sqrt((2.0 * n + 1)/(2.0 * n)) * (sqrt(1 - lambda[i] * lambda[i]) * DatadRnmr[(n - 1) * MaxN * I + (n - 1) * I + i]  + lambda[i] * DataRnmr[(n - 1) * MaxN * I + (n - 1) * I + i]);
    }

    for (n = 1; n < MaxN; n++) /// P_n^{n - 1}
    {
      DatadRnmr[n * MaxN * I + (n - 1) * I + i] = sqrt(2 * n + 1) * (-sqrt(1 - lambda[i] * lambda[i]) * DataRnmr[(n - 1) * MaxN * I + (n - 1) * I + i] + lambda[i] * DatadRnmr[(n - 1) * MaxN * I + (n - 1) * I + i]);
  }

    for (n = 2; n < MaxN; n++) 
    {
      for(m = 0; m <= n - 2; m++) 
      {
	DatadRnmr[n * MaxN * I + m * I + i] = sqrt((2.0 * n + 1)/(n*n - m * m * 1.0)) * (sqrt(2 * n - 1) * (- sqrt(1 - lambda[i] * lambda[i])  * DataRnmr[(n - 1) * MaxN * I + m * I + i] + lambda[i] * DatadRnmr[(n - 1) * MaxN * I + m * I + i]) - sqrt(1.0 * (n - 1 + m) * (n - 1 - m)/(2 * n - 3)) * DatadRnmr[(n - 2) * MaxN * I + m * I + i]);
      }
    }
  }

  /**
   * Rnmr的存储顺序 [N-m+1] * M * I
   * 这样写不容易索引？？ 顺序与求和一致， Be careful here！
   */
  ix = 0;
  for(i = 0; i < I; i++)
  {
    for (m = 0; m < M; m++)
    {
      for(n = m; n < N; n++)
      {
	Rnmr[ix] = DataRnmr[n * MaxN * I + m * I + i];
	ix++;
      }
    }
  }

  
  /**
   * dRnmr的存储顺序 [N-m+1] * M * I
   * 这样写不容易索引？？ 顺序与求和一致， Be careful here！
   */
  ix = 0;
  for(i = 0; i < I; i++)
  {
    for (m = 0; m < M; m++)
    {
      for(n = m; n < N; n++)
      {
	dRnmr[ix] = DatadRnmr[n * MaxN * I + m * I + i];
	/**
	if(n == 0)
	{
	  cout << dRnmr[ix] << endl;
	}
	*/

	ix++;
      }
    }
  }

  delete[] DataRnmr;
  delete[] DatadRnmr;

  double *DataPl = new double[(J + 2) * MaxN]();
  sprintf(path,"%s/Parameter/P_L_%d_J_%d.txt", DIR, MaxN, J);
  infile.open(path);
  ix = 0;
  if( ! infile ) cout << "Open Pl Error!" << endl;
  for(l = 0; l < MaxN; l++)
  {
    for(j = 0; j < J + 2; j++)
    {
      infile >> DataPl[ix];
      ix++;
    }
  }
  
  infile.close();

  /**
   * Pl的存储顺序 L * J
   * 这样写不容易索引？？ 顺序与求和一致， Be careful here！
   */
  ix = 0;
  for(j = 0; j < J + 2; j++)
  {
    for(l = 0; l < L; l ++)
    {
      Plp[ix] = DataPl[l * (J + 2) + j];


      /**
      if(l == 0)
      {
	cout << Plp[ix] << endl;
      }
      */
      ix++;
    }
  }

  delete[] DataPl;

  double *DatadPl = new double[(J + 2) * MaxN]();
  sprintf(path,"%s/Parameter/dP_L_%d_J_%d.txt", DIR, MaxN, J);
  infile.open(path);
  ix = 0;
  if( ! infile ) cout << "Open dPl Error!" << endl;

  for(l = 0; l < MaxN; l++)
  {
    for(j = 0; j < J + 2; j++)
    {
      infile >> DatadPl[ix];
      ix++;
    }
  }
  
  infile.close();

  /**
   * Pl的存储顺序 L * J
   * 这样写不容易索引？？ 顺序与求和一致， Be careful here！
   */
  ix = 0;
  for(j = 0; j < J + 2; j++)
  {
    for(l = 0; l < L; l++)
    {
      dPlp[ix] = DatadPl[l * (J + 2) + j];
      ix++;
    }
  }

  delete[] DatadPl;

  ///给X_m(\theta_k) 赋值, 存储顺序，给定\theta_k,所有的m
  ix = 0;
  for(k = 0; k < K; k++) 
  {
    for(m = 1 - M; m < 0; m++)
    {
      Xm[ix] = sin(abs(m)*theta[k]);
      ix++;
    }
    for(m = 0; m < M; m++)
    {
      Xm[ix] = cos(m*theta[k]);
      ix++; 
    }
  }

  ix = 0;
  for (k = 0; k < K; k++) 
  {
    for (m = 1 - M; m < 0; m++)
    {
      dXm[ix] = abs(m) * cos(abs(m)*theta[k]);
      ix++;
    }
    for(m = 0; m < M; m++)
    {
      dXm[ix] = - m * sin(m * theta[k]);  /// m = 0, dX = 0
      ix++;
    }
  }

  for(i = 0; i < I; i++)
  {
    for(j = 0; j < J; j++)
    {

      Jacobi[i * J + j] =  pow(a,3)/pow(cosh(realxi(p[j])) - cos(mu[i]), 3);
    }
  }

  for(i = 0; i < I; i++)
  {
    //  double tmp1, tmp2;
    //  tmp1 = - 1.0 * sin(mu[i]) / sqrt(1 - pow((D/2.0 * cos(mu[i]) - 1)/(D/2.0 - cos(mu[i])),2.0));
    //  tmp2 = - a * a / pow( D/2.0 - cos(mu[i]), 2.0);
    // surface_jacobi[i] = a / (D/2.0 - cos(mu[i])) * tmp1 * tmp2;
    // surface_jacobi[i] = - tmp2;
    surface_jacobi[i] = a * a / pow( cosh(xi0) - cos(mu[i]), 2.0);
    surface_jacobi[i + I] = sinh(xi1) * sinh(xi1) / pow( cosh(xi1) - cos(mu[i]), 2.0);
  
    // surface_jacobi[i + I] = a * a / pow( cosh(xi1) - cos(mu[i]), 2.0);
  }
  // cout << endl;
  
  //边界上的坐标
  /// 表面能的平移不变性
  xb = new double[I * K * 2];
  yb = new double[I * K * 2];
  zb = new double[I * K * 2];
 
  /// 球心的坐标(0, 0, \pm D/2)

  for(i = 0; i < I; i++)
  {
    for(k = 0; k < K; k++)
    {
      ///球心(0, 0, 0)
      xb[i * K + k] = a * sin(mu[i])/(cosh(xi0) - cos(mu[i])) * cos(theta[k]);
      yb[i * K + k] = a * sin(mu[i])/(cosh(xi0) - cos(mu[i])) * sin(theta[k]);
      zb[i * K + k] = a * a/(cosh(xi0) - cos(mu[i])) - a * cosh(xi0)/sinh(xi0);
      // cout << xb[i * K + k] * xb[i * K + k] + yb[i * K + k] * yb[i * K + k] + zb[i * K + k] * zb[i * K + k] << endl;

      ///球心(0,0,D)
      xb[I * K + i * K + k] = a * sin(mu[i])/(cosh(xi1) - cos(mu[i])) * cos(theta[k]) / R1;
      yb[I * K + i * K + k] = a * sin(mu[i])/(cosh(xi1) - cos(mu[i])) * sin(theta[k]) / R1;
      zb[I * K + i * K + k] = (a * sinh(xi1) /(cosh(xi1) - cos(mu[i])) - a * cosh(xi1)/sinh(xi1)) / R1;
      // cout << xb[I * K + i * K + k] * xb[I * K + i * K + k] + yb[I * K + i * K + k] * yb[I * K + i * K + k] + zb[I * K + i * K + k] * zb[I * K + i * K + k] << endl;
    }
  }


  double test = 0;
  for(j = 0; j < J; j++)
  {
    test += Plp[j * L + 31] * Plp[j * L + 31] * coe_p[j];
  }
  cout << test << endl; 

  m = 2;
  n = 127;
  int n1 = 127;
  int m1 = 0;
  test = 0;
  int absm = abs(m);


  // test = 0; 
  double *tmp;
  tmp = new double[I]();
  for(i = 0; i < I; i++) 
  {
    for(k = 0; k < K; k++)
    {
      tmp[i] += Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n1] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n] * Xm[k * (2 * M - 1) + m + M - 1] * Xm[k * (2 * M - 1) + m + M - 1] * 2 * PI / K;
    }
    
    // tmp[i] = Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n1] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n] * 2 * PI;
    test += tmp[i] * coe_mu[i];
  }
  cout << test << endl;

  delete[] tmp;
  

}

void Basis_destroy() 
{
  delete [] mu;
  delete [] p;
  delete [] theta;
  delete [] lambda;

  delete [] Jacobi;
 
  delete [] xb;
  delete [] yb;
  delete [] zb;
  delete [] coe_mu;
  delete [] coe_p;

  delete [] Rnmr;
  delete [] Plp;
  delete [] Xm;

  delete [] dRnmr;
  delete [] dPlp;
  delete [] dXm;

  delete [] Kijm;

  fftw_destroy_plan(FFTp_p);
  fftw_destroy_plan(FFTq_p);
  delete [] FDI_p;
  delete [] FDO_p;

  fftw_destroy_plan(FFTp_c);
  fftw_destroy_plan(FFTq_c);
  delete [] FDI_c;
  delete [] FDO_c;

 delete[] surface_jacobi;
}

///已知展开系数求函数值
/// fnlm = double[Basis] L * ((2 * M - 1) * N - M * (M - 1))
void calc_fijk(double *fnlm, double *fijk) 
{
  double *value, *fnm;
  int ix, ix1;

  value = new double[(2 * M - 1) * I * (J + 2)];  
  fnm = new double[((2 * M - 1) * N - M * (M - 1)) * (J + 2)];

  ix = 0;
  for(int j = 0; j < J + 2; j++)
  {
    ix1 = 0;
    for(int m = 1 - M; m < M; m++)
    {
      for(int n = abs(m); n < N; n++)
      {
	fnm[ix] = 0;
	for(int l = 0; l < L; l++)
	{	 
	  fnm[ix] += fnlm[l * ((2 * M - 1) * N - M * (M - 1)) + ix1] *  Plp[j * L + l];	
	}
	ix++;
	ix1++;
      }
    }
  }

  ix1 = 0;
  for(int i = 0; i < I; i++)
  {
    ix = 0;
    for(int j = 0; j < J + 2; j++)
    {
      for(int m = 1 - M; m < M; m++)
      {
	value[ix1] = 0;
	int absm = abs(m);
	for(int n = absm; n < N; n++)
	{
	  value[ix1] += fnm[ix] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];  ///R_n^{(m)}(r_i)
	  ix++;
	}
	ix1++;
      }
    }
  }

  /*
  for(int i = 0; i < I; i++) 
  {
    for(int j = 0; j < J; j++) 
    {
      for (int k = 0; k < K; k++) 
      {
	fijk[i * J * K + j * K + k] = 0;
	for (int m = 0; m < 2 * M - 1; m++) 
	{
	  fijk[i * J * K + j * K + k] += value[i * (J + 2) * (2 * M - 1) + j * (2 * M -1) + m] * Xm[k * (2 * M - 1) + m];
	}
      }
    }
  }
  */
 
  if (K >= 2 * M) 
  {
    for (int i = 0; i < I; i++) 
    {
      for (int j = 0; j < J; j++) 
      {
	for (int m = 0; m < M; m++)
	  FDI_p[m] = value[i * (J + 2) * (2 * M - 1) + j * (2 * M -1) + m + M - 1];
	for (int m = M; m < K/2 + 1; m++)
	  FDI_p[m] = 0;
	fftw_execute(FFTp_p);

	for(int k = 0; k < K/2; k++)
	  fijk[i * J * K + j * K + k] = 0.5 * (FDO_p[k] + FDI_p[0]);
	
	/// pi to 2 pi
	for (int m = 1; m < M; m+=2)
	  FDI_p[m] = -FDI_p[m];
	fftw_execute(FFTp_p);
	for (int k = 0; k < K/2; k++)
	  fijk[i * J * K + j * K + k + K/2] = 0.5 * (FDO_p[k] + FDI_p[0]);
       

	for(int m = 0; m < M - 1; m++)
	{
	  FDI_p[m] = value[i * (J + 2) * (2 * M - 1) + j * (2 * M -1) + M - 2 - m];
	}
	for(int m = M - 1; m < K/2 + 1; m++)
	{
	  FDI_p[m] = 0;
	}
	fftw_execute(FFTq_p);
	for(int k = 1; k < K/2; k++)
	  fijk[i * J * K + j * K + k] += 0.5 * FDO_p[k-1];
      
	for(int m = 0; m < M - 1; m+=2)
	{
	  FDI_p[m] = -FDI_p[m];
	}
	fftw_execute(FFTq_p);
	for(int k = 1; k < K/2; k++)
	  fijk[i * J * K + j * K + k + K/2] += 0.5 * FDO_p[k-1];
      }
    }
  }
  else 
  {
    for(int i = 0; i < I; i++) 
    {
      for(int j = 0; j < J; j++) 
      {
	for(int k = 0; k < K; k++) 
	{
	  fijk[i * J * K + j * K + k] = 0;
	  for(int m = 0; m < 2 * M - 1; m++) 
	  {
	    fijk[i * J * K + j * K + k] += value[i * (J + 2) * (2 * M - 1) + j * (2 * M -1) + m] * Xm[k * (2 * M - 1) + m];
	  }
	}
      }
    }
  }

  int Boundary_index = I * J * K;
  for(int j = 0; j < 2; j++)
  {
    for(int i = 0; i < I; i++)
    {
      for(int k = 0; k < K; k++)
      {
	fijk[Boundary_index + j * I * K + i * K + k] = 0;
	for (int m = 0; m < 2 * M - 1; m++) 
	{
	  fijk[Boundary_index + j * I * K + i * K + k] += value[i * (J + 2) * (2 * M - 1) + (J + j) * (2 * M -1) + m] * Xm[k * (2 * M - 1) + m];
	}
      }
    }
  }

  delete[] value;
  delete[] fnm;
}

void calc_dr_fijk(double *fnlm, double *dr_fijk) 
{
  double *value, *fnm;
  int ix, ix1;

  value = new double[(2 * M - 1) * I * J];  ///存储R_n^{|m|}(r_i)
  fnm = new double[((2 * M - 1) * N - M * (M - 1)) * J];

  ix = 0;
  for(int j = 0; j < J; j++)
  {
    ix1 = 0;
    for(int m = 1 - M; m < M; m++)
    {
      for(int n = abs(m); n < N; n++)
      {
	fnm[ix] = 0;
	for(int l = 0; l < L; l++)
	{
	  fnm[ix] += fnlm[l * ((2 * M - 1) * N - M * (M - 1)) + ix1] * Plp[j * L + l];	
	}
	ix++;
	ix1++;
      }
    }
  }

  ix1 = 0;
  for(int i = 0; i < I; i++)
  {
    ix = 0;
    for(int j = 0; j < J; j++)
    {
      for(int m = 1 - M; m < M; m++)
      {
	value[ix1] = 0;
	for(int n = abs(m); n < N; n++)
	{
	  int absm = abs(m);
	  value[ix1] += fnm[ix] * dRnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];  ///R_n^{(m)}(r_i)
	  ix++;
	}
	ix1++;
      }
    }
  }

  /*
  for(int i = 0; i < I; i++) 
  {
    for(int j = 0; j < J; j++) 
    {
      for (int k = 0; k < K; k++) 
      {
	dr_fijk[i * J * K + j * K + k] = 0;
	for (int m = 0; m < 2 * M - 1; m++) 
	{
	  dr_fijk[i * J * K + j  * K + k] += value[i * J * (2 * M - 1) + j * (2 * M - 1) + m] * Xm[k * (2 * M - 1) + m];
	}
      }
    }
  }
  */

  if (K >= 2 * M) 
  {
    for (int i = 0; i < I; i++) 
    {
      for (int j = 0; j < J; j++) 
      {
	for (int m = 0; m < M; m++)
	  FDI_p[m] = value[i * J * (2 * M - 1) + j * (2 * M -1) + m + M - 1];
	for (int m = M; m < K/2 + 1; m++)
	  FDI_p[m] = 0;
	fftw_execute(FFTp_p);

	for(int k = 0; k < K/2; k++)
	  dr_fijk[i * J * K + j * K + k] = 0.5 * (FDO_p[k] + FDI_p[0]);
	
	/// pi to 2 pi
	for (int m = 1; m < M; m+=2)
	  FDI_p[m] = -FDI_p[m];
	fftw_execute(FFTp_p);
	for (int k = 0; k < K/2; k++)
	  dr_fijk[i * J * K + j * K + k + K/2] = 0.5 * (FDO_p[k] + FDI_p[0]);
       
	for(int m = 0; m < M - 1; m++)
	{
	  FDI_p[m] = value[i * J * (2 * M - 1) + j * (2 * M -1) + M - 2 - m];
	}
	for(int m = M - 1; m < K/2 + 1; m++)
	{
	  FDI_p[m] = 0;
	}
	fftw_execute(FFTq_p);
	for(int k = 1; k < K/2; k++)
	  dr_fijk[i * J * K + j * K + k] += 0.5 * FDO_p[k-1];
      
	for(int m = 0; m < M - 1; m+=2)
	{
	  FDI_p[m] = -FDI_p[m];
	}
	fftw_execute(FFTq_p);
	for(int k = 1; k < K/2; k++)
	  dr_fijk[i * J * K + j * K + k + K/2] += 0.5 * FDO_p[k-1];
      }
    }
  }
  else 
  {
    for(int i = 0; i < I; i++) 
    {
      for(int j = 0; j < J; j++) 
      {
	for(int k = 0; k < K; k++) 
	{
	  dr_fijk[i * J * K + j * K + k] = 0;
	  for(int m = 0; m < 2 * M - 1; m++) 
	  {
	    dr_fijk[i * J * K + j * K + k] += value[i * J * (2 * M - 1) + j * (2 * M -1) + m] * Xm[k * (2 * M - 1) + m];
	  }
	}
      }
    }
  }

  delete[] value;
  delete[] fnm;
}

void calc_dp_fijk(double *fnlm,double *dp_fijk) 
{
  double *value, *fnm;
  int ix, ix1;

  value = new double[(2 * M - 1) * I * J];  ///存储R_n^{|m|}(r_i)
  fnm = new double[((2 * M - 1) * N - M * (M - 1)) * J];

  ix = 0;
  for(int j = 0; j < J; j++)
  {
    ix1 = 0;
    for(int m = 1 - M; m < M; m++)
    {
      for(int n = abs(m); n < N; n++)
      {
	fnm[ix] = 0;
	for(int l = 0; l < L; l++)
	{
	  fnm[ix] += fnlm[l * ((2 * M - 1) * N - M * (M - 1)) + ix1] * dPlp[j * L + l];	
	}
	ix++;
	ix1++;
      }
    }
  }

  ix1 = 0;
  for(int i = 0; i < I; i++)
  {
    ix = 0;
    for(int j = 0; j < J; j++)
    {
      for(int m = 1 - M; m < M; m++)
      {
	value[ix1] = 0;
	for(int n = abs(m); n < N; n++)
	{
	  int absm = abs(m);
	  value[ix1] += fnm[ix] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];  ///R_n^{(m)}(r_i)
	  ix++;
	}
	ix1++;
      }
    }
  }

  /*
  for(int i = 0; i < I; i++) 
  {
    for(int j = 0; j < J; j++) 
    {
      for (int k = 0; k < K; k++) 
      {
	dp_fijk[i * J * K + j * K + k] = 0;
	for (int m = 0; m < 2 * M - 1; m++) 
	{
	  dp_fijk[i * J * K + j * K + k] += value[i * J * (2 * M - 1) + j * (2 * M -1) + m] * Xm[k * (2 * M - 1) + m];
	}
      }
    }
  }
  */

  if (K >= 2 * M) 
  {
    for (int i = 0; i < I; i++) 
    {
      for (int j = 0; j < J; j++) 
      {
	for (int m = 0; m < M; m++)
	  FDI_p[m] = value[i * J * (2 * M - 1) + j * (2 * M -1) + m + M - 1];
	for (int m = M; m < K/2 + 1; m++)
	  FDI_p[m] = 0;
	fftw_execute(FFTp_p);

	for(int k = 0; k < K/2; k++)
	  dp_fijk[i * J * K + j * K + k] = 0.5 * (FDO_p[k] + FDI_p[0]);
	
	/// pi to 2 pi
	for (int m = 1; m < M; m+=2)
	  FDI_p[m] = -FDI_p[m];
	fftw_execute(FFTp_p);
	for (int k = 0; k < K/2; k++)
	  dp_fijk[i * J * K + j * K + k + K/2] = 0.5 * (FDO_p[k] + FDI_p[0]);
       

	for(int m = 0; m < M - 1; m++)
	{
	  FDI_p[m] = value[i * J * (2 * M - 1) + j * (2 * M -1) + M - 2 - m];
	}
	for(int m = M - 1; m < K/2 + 1; m++)
	{
	  FDI_p[m] = 0;
	}
	fftw_execute(FFTq_p);
	for(int k = 1; k < K/2; k++)
	  dp_fijk[i * J * K + j * K + k] += 0.5 * FDO_p[k-1];
      
	for(int m = 0; m < M - 1; m+=2)
	{
	  FDI_p[m] = -FDI_p[m];
	}
	fftw_execute(FFTq_p);
	for(int k = 1; k < K/2; k++)
	  dp_fijk[i * J * K + j * K + k + K/2] += 0.5 * FDO_p[k-1];
      }
    }
  }
  else 
  {
    for(int i = 0; i < I; i++) 
    {
      for(int j = 0; j < J; j++) 
      {
	for(int k = 0; k < K; k++) 
	{
	  dp_fijk[i * J * K + j * K + k] = 0;
	  for(int m = 0; m < 2 * M - 1; m++) 
	  {
	    dp_fijk[i * J * K + j * K + k] += value[i * J * (2 * M - 1) + j * (2 * M -1) + m] * Xm[k * (2 * M - 1) + m];
	  }
	}
      }
    }
  }
 

  delete[] value;
  delete[] fnm;
}

void calc_dt_fijk(double *fnlm,double *dt_fijk) 
{
  double *value, *fnm;
  int ix, ix1;

  value = new double[(2 * M - 1) * I * J];  ///存储R_n^{|m|}(r_i)
  fnm = new double[((2 * M - 1) * N - M * (M - 1)) * J];

  ix = 0;
  for(int j = 0; j < J; j++)
  {
    ix1 = 0;
    for(int m = 1 - M; m < M; m++)
    {
      for(int n = abs(m); n < N; n++)
      {
	fnm[ix] = 0;
	for(int l = 0; l < L; l++)
	{
	  fnm[ix] += fnlm[l * ((2 * M - 1) * N - M * (M - 1)) + ix1] * Plp[j * L + l];	
	}
	ix++;
	ix1++;
      }
    }
  }

  ix1 = 0;
  for(int i = 0; i < I; i++)
  {
    ix = 0;
    for(int j = 0; j < J; j++)
    {
      for(int m = 1 - M; m < M; m++)
      {
	value[ix1] = 0;
	for(int n = abs(m); n < N; n++)
	{
	  int absm = abs(m);
	  value[ix1] += fnm[ix] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];  ///R_n^{(m)}(r_i)
	  ix++;
	}
	ix1++;
      }
    }
  }

  for(int i = 0; i < I; i++) 
  {
    for(int j = 0; j < J; j++) 
    {
      for (int k = 0; k < K; k++) 
      {
	dt_fijk[i * J * K + j * K + k] = 0;
	for (int m = 0; m < 2 * M - 1; m++) 
	{
	  dt_fijk[i * J * K + j * K + k] += value[i * J * (2 * M - 1) + j * (2 * M -1) + m] * dXm[k*(2*M - 1) + m];
	}
      }
    }
  }

  delete[] value;
  delete[] fnm;
}

double IntCylind(double *f) 
{
  int i,j,k;
  double intf = 0;
  double *G1 = new double[J*K]();
  
  intf = 0;
  for(j = 0; j < J; j++) 
  {
    for(k = 0; k < K; k++) 
    {
      G1[j * K + k] = 0;
      for(i = 0; i < I; i++)
      {
	G1[j * K + k] += coe_mu[i] * f[i * J * K + j * K + k];
      }
    }
  }

  for (k = 0; k < K; k++) 
  {
    for (j = 0; j < J; j++)
    {
      intf += coe_p[j] * G1[j * K + k];
    }
  }

  intf = intf * 2.0 * PI/K;
  
  delete[] G1;
  
  return intf;
}


void calc_fnlm(double *fijk, double *fnlm) 
{
  int i,j,k,n,l,m;
  double *value;
  int ix,ix1;

  value = new double[((2 * M - 1) * N - M * (M - 1)) * J];

  /*
  for(m = 1 - M; m < M; m++)
  {
    for(i = 0; i < I; i++) 
    {
      for(j = 0; j < J; j++) 
      {
	Kijm[(m + M - 1) * I * J + i * J + j] = 0;
	for(k = 0; k < K; k++)
	{
	  Kijm[(m + M - 1) * I * J + i * J + j] += fijk[i * J * K + j * K + k] * Xm[k * (2 * M - 1) + m + M - 1] * 2 * PI / K;
	}
      }
    }
  }
  */

  for(i = 0; i < I; i++) 
  {
    for(j = 0; j < J; j++) 
    {
      for(k = 0; k < K; k++)
  	FDI_c[k] = fijk[i * J * K + j * K + k];
      fftw_execute(FFTp_c);
      for (m = 0; m < M; m++)
      	Kijm[(m + M - 1) * I * J + i * J + j] = 0.5 * (FDO_c[2*m] + FDI_c[0]);
      
      for (k = 0; k < K - 1; k++)
      	FDI_c[k] = fijk[i * J * K + j * K + k + 1];
      fftw_execute(FFTq_c);
      for (m = 1; m < M; m++)
      	Kijm[(M - 1 - m) * I * J + i * J + j] = 0.5 * FDO_c[2 * m-1];
    }
  }
    
  ix = 0;
  for(m = 1 - M; m < M; m++)
  {
    for(n = abs(m); n < N; n++)
    {
      for(j = 0; j < J; j++)
      {
	value[ix] = 0;
	for(i = 0; i < I; i++)
	{
	  int absm = abs(m);
	  value[ix] += Kijm[(m + M - 1) * I * J + i * J + j] * coe_mu[i] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];
	}
	ix++;
      }
    }
  }

  ix = 0;
  ix1 = 0;
  for(l = 0; l < L; l++)
  {
    ix1 = 0;
    for(m = 1 - M; m < M;  m++)
    {
      for(n = abs(m); n < N; n++)
      {
	fnlm[ix] = 0;
	for(j = 0; j < J; j++)
	{ 
	  fnlm[ix] += value[ix1] * Plp[j * L + l] * coe_p[j]; 
	  ix1++;
	}
	fnlm[ix] = fnlm[ix] * 2 * PI / K;

	ix++;
      }
    }
  }

  delete[] value;
}


void calc_gnlm(double *fijk, double *fnlm) 
{
  int i, j, k, n, l, m;
  double *value;
  int ix, ix1;

  value = new double[((2 * M - 1) * N - M * (M - 1)) * J];

  /**
  for(m = 1 - M; m < M; m++)
  {
    for(i = 0; i < I; i++) 
    {
      for(j = 0; j < J; j++) 
      {
	Kijm[(m + M - 1) * I * J + i * J + j] = 0;
	for(k = 0; k < K; k++)
	{
	  Kijm[(m + M - 1) * I * J + i * J + j] += fijk[i * J * K + j * K + k] * Xm[k * (2 * M - 1) + m + M - 1] * 2 * PI / K;
	}
      }
    }
  }
  */

  for(i = 0; i < I; i++) 
  {
    for(j = 0; j < J; j++) 
    {
      for(k = 0; k < K; k++)
      {
	FDI_c[k] = fijk[i * J * K + j * K + k];

	/**
	if(fijk[i * J * K + j * K + k] > 1e-9)
	{
	  FDI_c[k] = fijk[i * J * K + j * K + k];
	}
	else
	{
	  FDI_c[k] = 0;
	}
	*/
      }
      fftw_execute(FFTp_c);
      for (m = 0; m < M; m++)
      	Kijm[(m + M - 1) * I * J + i * J + j] = 0.5 * (FDO_c[2*m] + FDI_c[0]);
      
      for (k = 0; k < K - 1; k++)
      	FDI_c[k] = fijk[i * J * K + j * K + k + 1];
      fftw_execute(FFTq_c);
      for (m = 1; m < M; m++)
      	Kijm[(M - 1 - m) * I * J + i * J + j] = 0.5 * FDO_c[2 * m-1];
    }
  }

  ix = 0;
  for(m = 1 - M; m < M; m++)
  {
    for(n = abs(m); n < N; n++)
    {
      for(j = 0; j < J; j++)
      {
	value[ix] = 0;
	for(i = 0; i < I; i++)
	{
	  int absm = abs(m);
	  value[ix] += Kijm[(m + M - 1) * I * J + i * J + j] * Jacobi[i * J + j] * coe_mu[i] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];
	}
	ix++;
      }
    }
  }
 
  ix = 0;
  ix1 = 0;
  for(l = 0; l < L; l++)
  {
    ix1 = 0;
    for(m = 1 - M; m < M;  m++)
    {
      for(n = abs(m); n < N; n++)
      {
	fnlm[ix] = 0;
	for(j = 0; j < J; j++)
	{ 
	  fnlm[ix] += value[ix1] * Plp[j * L + l] * coe_p[j]; 
	  ix1++;
	}
	fnlm[ix] = fnlm[ix] * 2 * PI / K;

       
	/*
        if(fnlm[ix] < 1e-9)
        {
          fnlm[ix] = 0;
        }
	*/

	ix++;
      }
    }
  }

  //  fnlm[const_idx] = 0;

  delete[] value;
}

void calc_dr_gnlm(double *fijk, double *dr_fnlm) 
{
  int i,j,k,n,l,m;
  double *value;
  int ix,ix1;
  int idx;

  value = new double[((2 * M - 1) * N - M * (M - 1)) * J];

  /*
  for(m = 1 - M; m < M; m++)
  {
    for(i = 0; i < I; i++) 
    {
      for(j = 0; j < J; j++) 
      {
	Kijm[(m + M - 1) * I * J + i * J + j] = 0;
	for(k = 0; k < K; k++)
	{
	  Kijm[(m + M - 1) * I * J + i * J + j] += fijk[i * J * K + j * K + k] * Xm[k * (2 * M - 1) + m + M - 1] * 2 * PI / K; 
	}
      }
    }
  }
  */

  for(i = 0; i < I; i++) 
  {
    for(j = 0; j < J; j++) 
    {
      for(k = 0; k < K; k++)
      {
	FDI_c[k] = fijk[i * J * K + j * K + k];

	/**
	if(fijk[i * J * K + j * K + k] > 1e-9)
	{
	  FDI_c[k] = fijk[i * J * K + j * K + k];
	}
	else
	{
	  FDI_c[k] = 0;
	}
	*/
      }
      fftw_execute(FFTp_c);
      for (m = 0; m < M; m++)
      	Kijm[(m + M - 1) * I * J + i * J + j] = 0.5 * (FDO_c[2*m] + FDI_c[0]);
      
      for (k = 0; k < K - 1; k++)
      	FDI_c[k] = fijk[i * J * K + j * K + k + 1];
      fftw_execute(FFTq_c);
      for (m = 1; m < M; m++)
      	Kijm[(M - 1 - m) * I * J + i * J + j] = 0.5 * FDO_c[2 * m-1];
    }
  }
  
  ix = 0;
  for(m = 1 - M; m < M; m++)
  {
    for(n = abs(m); n < N; n++)
    {
      for(j = 0; j < J; j++)
      {
	value[ix] = 0;
	for(i = 0; i < I; i++)
	{
	  int absm = abs(m);
	  value[ix] += Kijm[(m + M - 1) * I * J + i * J + j] * Jacobi[i * J + j] * coe_mu[i] * dRnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];
	}
	ix++;
      }
    }
  }
 
  ix = 0;
  for(l = 0; l < L; l++)
  {
    ix1 = 0;
    for(m = 1 - M; m < M;  m++)
    {
      for(n = abs(m); n < N; n++)
      {
	dr_fnlm[ix] = 0;
	for(j = 0; j < J; j++)
	{ 
	  dr_fnlm[ix] += value[ix1] * Plp[j * L + l] * coe_p[j]; 
	  ix1++;
	}
	dr_fnlm[ix] = dr_fnlm[ix] * 2 * PI / K;

	/*
        if(dr_fnlm[ix] < 1e-9)
        {
          dr_fnlm[ix] = 0;
        }
	*/

	ix++;
      }
    }
  }

  //  dr_fnlm[const_idx] = 0;

  delete[] value;
}

void calc_dp_gnlm(double *fijk, double *dp_fnlm) 
{
  int i,j,k,n,l,m;
  double *value;
  int ix,ix1;

  value = new double[((2 * M - 1) * N - M * (M - 1)) * J];

  /*
  for(m = 1 - M; m < M; m++)
  {
    for(i = 0; i < I; i++) 
    {
      for(j = 0; j < J; j++) 
      {
	Kijm[(m + M - 1) * I * J + i * J + j] = 0;
	for(k = 0; k < K; k++)
	{
	  Kijm[(m + M - 1) * I * J + i * J + j] += fijk[i * J * K + j * K + k] * Xm[k * (2 * M - 1) + m + M - 1] * 2 * PI / K;
	}
      }
    }
  }
  */

  for(i = 0; i < I; i++) 
  {
    for(j = 0; j < J; j++) 
    {
      for(k = 0; k < K; k++)
      {
	FDI_c[k] = fijk[i * J * K + j * K + k];

	/**
	if(fijk[i * J * K + j * K + k] > 1e-9)
	{
	  FDI_c[k] = fijk[i * J * K + j * K + k];
	}
	else
	{
	  FDI_c[k] = 0;
	}
	*/
      }
      fftw_execute(FFTp_c);
      for (m = 0; m < M; m++)
      	Kijm[(m + M - 1) * I * J + i * J + j] = 0.5 * (FDO_c[2*m] + FDI_c[0]);
      
      for (k = 0; k < K - 1; k++)
      	FDI_c[k] = fijk[i * J * K + j * K + k + 1];
      fftw_execute(FFTq_c);
      for (m = 1; m < M; m++)
      	Kijm[(M - 1 - m) * I * J + i * J + j] = 0.5 * FDO_c[2 * m-1];
    }
  }

  ix = 0;
  for(m = 1 - M; m < M; m++)
  {
    for(n = abs(m); n < N; n++)
    {
      for(j = 0; j < J; j++)
      {
	value[ix] = 0;
	for(i = 0; i < I; i++)
	{
	  int absm = abs(m);
	  /// 	  value[ix] += Kijm[(m + M - 1) * I * J + i * J + j] * pow(a,3) * Jacobi[i * J + j] * coe_mu[i] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];  
	  value[ix] += Kijm[(m + M - 1) * I * J + i * J + j] * Jacobi[i * J + j] * coe_mu[i] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];
	}
	ix++;
      }
    }
  }
 
  ix = 0;
  ix1 = 0;
  for(l = 0; l < L; l++)
  {
    ix1 = 0;
    for(m = 1 - M; m < M;  m++)
    {
      for(n = abs(m); n < N; n++)
      {
	dp_fnlm[ix] = 0;
	for(j = 0; j < J; j++)
	{ 
	  dp_fnlm[ix] += value[ix1] * dPlp[j * L + l] * coe_p[j]; 
	  ix1++;
	}
	dp_fnlm[ix] = dp_fnlm[ix] * 2 * PI / K;

	/*
        if(dp_fnlm[ix] < 1e-9)
        {
          dp_fnlm[ix] = 0;
        }
	*/

	ix++;
      }
    }
  }

  // dp_fnlm[const_idx] = 0;

  delete[] value;
}

void calc_dt_gnlm(double *fijk, double *dt_fnlm) 
{
  int i,j,k,n,l,m;
  double *value;
  int ix,ix1;
  int idx;

  value = new double[((2 * M - 1) * N - M * (M - 1)) * J]();
  
  for(m = 1 - M; m < M; m++)
  {
    for(i = 0; i < I; i++) 
    {
      for(j = 0; j < J; j++) 
      {
	Kijm[(m + M - 1) * I * J + i * J + j] = 0;
	for(k = 0; k < K; k++)
	{
	  Kijm[(m + M - 1) * I * J + i * J + j] += fijk[i * J  * K + j * K + k] * dXm[k * (2 * M - 1) + m + M - 1] * 2 * PI / K;
	  /**
	  if(fijk[i * J * K + j * K + k] > 1e-9)
	  {
	    Kijm[(m + M - 1) * I * J + i * J + j] += fijk[i * J  * K + j * K + k] * dXm[k * (2 * M - 1) + m + M - 1] * 2 * PI / K;
	  }
	  */
	}
      }
    }
  }

  ix = 0;
  for(m = 1 - M; m < M; m++)
  {
    for(n = abs(m); n < N; n++)
    {
      for(j = 0; j < J; j++)
      {
	value[ix] = 0;
	for(i = 0; i < I; i++)
	{
	  int absm = abs(m);
	  value[ix] += Kijm[(m + M - 1) * I * J + i * J + j] * Jacobi[i * J + j] * coe_mu[i] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n];
	}
	ix++;
      }
    }
  }
  
  ix = 0;
  ix1 = 0;
  for(l = 0; l < L; l++)
  {
    ix1 = 0;
    for(m = 1 - M; m < M;  m++)
    {
      for(n = abs(m); n < N; n++)
      {
	dt_fnlm[ix] = 0;
	for(j = 0; j < J; j++)
	{ 
	  dt_fnlm[ix] += value[ix1] * Plp[j * L + l] * coe_p[j]; 
	  ix1++;
	}

	/*
        if(dt_fnlm[ix] < 1e-9)
        {
          dt_fnlm[ix] = 0;
        }
	*/

	ix++;
      }
    }
  }

  // dt_fnlm[const_idx] = 0;

  delete[] value;
}

double IntBisphereJacobi(double *f) 
{
  int i,j,k;
  double intf = 0;
  double *G1 = new double[J * K];
  
  intf = 0;
  for(j = 0; j < J; j++) 
  {
    for(k = 0; k < K; k++) 
    {
      G1[j * K + k] = 0;
      for(i = 0; i < I; i++)
      {
	G1[j * K + k] += coe_mu[i] * f[i * J * K + j * K + k] *  Jacobi[i * J + j];
      }
    }
  }

  for (k = 0; k < K; k++) 
  {
    for (j = 0; j < J; j++)
    {
      intf += coe_p[j] * G1[j * K + k];
    }
  }

  intf = intf * 2 * PI/K;
  
  delete[] G1;
  
  return intf;
}

void calc_gnm(double *fik, double *f_nm, int s) 
{
  /// 球面上的积分
  // int i, k, n, m;
  int ix;
  double *Pim = new double[(2 * M - 1) * I];
  
  for(int m = 0; m < 2 * M - 1; m++) 
  {
    for (int i = 0; i < I; i++) 
    {
      Pim[m * I + i] = 0;
      for(int k = 0; k < K; k++)
      {
      	Pim[m * I + i] += fik[i * K + k] * Xm[k * (2 * M - 1) + m];
      }
      Pim[m * I + i] = 2 * PI / K * Pim[m * I + i];
    }
  }

  ix = 0;
  for(int m = 0; m < 2 * M - 1; m++) 
  {
    int absm = abs(m - M + 1);
    for(int n = absm; n < N; n++) 
    {
      f_nm[ix] = 0;
      for(int i = 0; i < I; i++)
      {
  	f_nm[ix] += Pim[m * I + i] * coe_mu[i] * surface_jacobi[i + s * I] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n]; 
      }
      ix++;
    }
  }
  
  delete [] Pim;
}

void calc_dt_gnm(double *fik, double *f_nm, int s) 
{
  /// 球面上的积分
  // int i, k, n, m;
  int ix;
  double *Pim = new double[(2 * M - 1) * I];
  
  for(int m = 0; m < 2 * M - 1; m++) 
  {
    for (int i = 0; i < I; i++) 
    {
      Pim[m * I + i] = 0;
      for(int k = 0; k < K; k++)
      {
      	Pim[m * I + i] += fik[i * K + k] * dXm[k * (2 * M - 1) + m];
      }
      Pim[m * I + i] = 2 * PI / K * Pim[m * I + i];
    }
  }

  ix = 0;
  for(int m = 0; m < 2 * M - 1; m++) 
  {
    int absm = abs(m - M + 1);
    for(int n = absm; n < N; n++) 
    {
      f_nm[ix] = 0;
      for(int i = 0; i < I; i++)
      {
  	f_nm[ix] += Pim[m * I + i] * coe_mu[i] * surface_jacobi[i + s * I] * Rnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n]; 
      }
      ix++;
    }
  }
  
  delete [] Pim;
}

void calc_dr_gnm(double *fik, double *f_nm, int s) 
{
  /// 球面上的积分
  // int i, k, n, m;
  int ix;
  double *Pim = new double[(2 * M - 1) * I];
  
  for(int m = 0; m < 2 * M - 1; m++) 
  {
    for (int i = 0; i < I; i++) 
    {
      Pim[m * I + i] = 0;
      for(int k = 0; k < K; k++)
      {
      	Pim[m * I + i] += fik[i * K + k] * Xm[k * (2 * M - 1) + m];
      }
      Pim[m * I + i] = 2 * PI / K * Pim[m * I + i];
    }
  }

  ix = 0;
  for(int m = 0; m < 2 * M - 1; m++) 
  {
    int absm = abs(m - M + 1);
    for(int n = absm; n < N; n++) 
    {
      f_nm[ix] = 0;
      for(int i = 0; i < I; i++)
      {
  	f_nm[ix] += Pim[m * I + i] * coe_mu[i] * surface_jacobi[i + s * I] * dRnmr[i * (N*M - M*(M-1)/2) + N * absm - absm * (absm + 1)/2 + n]; 
      }
      ix++;
    }
  }
  
  delete [] Pim;
}





