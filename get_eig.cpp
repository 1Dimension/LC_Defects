/**
 * @file   main.cpp
 1;3201;0c* @author onedimension <onedimension@onedimension-PC>
 * @date   Tue Jun 24 21:18:04 2014 
 * @brief 主程序  
 * 
 * 
 */


/*********************************************************************/
/** 
 *  
 * @param A  Zernike多项式的系数 
 * @param V 
 * @param Q 
 * @param Fbulk Bulk_Energy的表达式
 * @param Felas 弹性能的表达式
 * @param Fpena 罚函数(Surface Energy)
 * @param Fbeta 
 * 
 *
 */
using namespace std;

#include <iomanip>

#include <cstdio>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "lbfgs.h"

#include "global.h"
#include "Basis.h"
#include "initial.h"
#include "calculate_Fb.h"
#include "calculate_Fel.h"
#include "calculate_Fs.h"
// #include "calculate_f.h"

#include "EigInfoTwo.h"
#include "in_and_out.h"
#include "calcvector.h"
#include "point.h"
// #include "rotation.h"
// #include "check.h"


double cal_F(double *A, double *Q, double& Fbulk, double& Felas) 
{
  
  Fbulk = calc_energy_Fb(A, Q, fbulk);  //根据展开系数计算bulk energy
  Felas = calc_F_el(A); //One-constant approximation
  

  //   Fbulk =  calc_f(A, Q);

  Fs = calc_energy_Fs(Q);

  //   cout << Fbulk << " "  <<  Felas << " " <<  eta * Fs / xi0 << endl;
  // return xi0 * Fbulk + eta * Fs;
  return (xi1 - xi0)/2 * (Fbulk + Felas) + eta * Fs;
  //   return Felas;
  //   cout << Fs << endl;
  // return Fs;
  //   return Fbulk;
  //   return Fbulk + Felas;
}

/*********************************************************************/

void cal_dF(double *A, double *Q,double *grad_Energy) 
{
  // double eps = 0.1;
  //  calc_grad_bulkenergy_landau_de(Q, grad_fbulk);
  
  for (int i = 0; i < 5; i++)
    calc_gnlm(grad_fbulk + i * innerPoint, grad_Fb + i * Basis);  ///计算BULK ENERGY对展开式系数的导数 
  
  for(int i = 0; i < 5 * Basis; i++)
  {
    grad_Energy[i] = (xi1 - xi0)/2 * (grad_Fb[i] + grad_FeL1[i]) + eta * grad_Fs[i];
    //   grad_Energy[i] = (xi1 - xi0)/2 * (grad_Fb[i]) + eta * grad_Fs[i];
    //   grad_Energy[i] = grad_FeL1[i];
    // grad_Energy[i] = grad_Fs[i];
    //   grad_Energy[i] = grad_Fb[i];
    //   grad_Energy[i] = grad_Fb[i] + 0.5 * ksi * grad_FeL1[i];
  }

    //   grad_Energy[i] = (xi1 - xi0)/2 * (g
}


static lbfgsfloatval_t evaluate(void *instance,
				const double *Bnlm,
				double *grad_Energy,
				const int n,
				const lbfgsfloatval_t step) 
{
  double Energy;
  Energy = cal_F(Anlm, Qijk, Fbulk, Felas);
  cal_dF(Anlm, Qijk, grad_Energy);
  return Energy;
}

static int progress(void *instance,
		    double *Anlm,
 		    double *grad_Energy,
		    double Energy,
		    const lbfgsfloatval_t xnorm,
		    const lbfgsfloatval_t gnorm,
		    const lbfgsfloatval_t step,
		    int n,
		    int k,
		    int ls) 
{
 
 
   
  /**
  if(k % 200 == 0)
  {
    oput(Anlm, landau_t, Rad, eta, L, M, N, filename);
  }
  */

  return 0;
}

/********************************************************************/
int main(int argc,char *argv[]) 
{

  // FILE *fp = fopen("Log.txt","a+");

  char fname[200];
  //  char filename[200];
  FILE *fp;

  // filename = argv[1];


  // cout << "Filename: " <<  filename << endl;
  
  //设定物理参数
  double A = - 0.172 * 1e6, B = -2.12 * 1e6, C = 1.73 * 1e6, L1 = 4 * 1e-11, R0 = 5 * 1e-5, W = 1e-2;
  //  double A, B = - 0.816 * 3 * 1e6, C = 0.45 * 4 * 1e6, L1 = 6 * 1e-12, R0 =  * 1e-7, W = 10;
  // A = 1.01 *  * B / C / 27.0; 
  // double D0 = 0.5 * R0;  ///球心间的距离

  L21 = 5;

  R1 = 0.7;
  double ep = 0.1;

  // theta_fd = PI/4.0;

  cout << "R = " <<  R1 << endl;
  // a = sqrt(pow(D/2.0, 2) - 1);         ///Bispherical坐标的系数
  a = sqrt((pow(1.0 - R1*R1 - ep * ep, 2.0) - 4 * ep * ep * R1 * R1))/(2.0 * ep);
  // cout << "D =" << D << endl;
  cout << "a =" << a << endl;
  xi0 = asinh(a);
  xi1 = asinh(a/R1); //// 
  // xi1 = 0;

  // xi0 = -log(sqrt(a * a + 1.0) - a );
  // xi1 = -log(sqrt(a*a/R1/R1 + 1) - a / R1);
 

  cout << "xi0 = " << xi0 << endl;  
  cout << "Xi1 = " << xi1 << endl;
  // cout << "D =" << a * cosh(xi1)/sinh(xi1) + D/2 - 1 - R1<< endl;
 

  cout << a * cosh(xi0)/sinh(xi0) - a * cosh(xi1)/sinh(xi1) << endl;
  // cout <<  << endl;

  landau_t = 27 * A * C / (B * B);

  //  landau_t = sqrt(6) - 4/3.0;
  //  landau_t = 2 * sqrt(6) - 16/3.0; 
  cout<<"t = " << landau_t << endl;

  /**
  eta = 0.5 * 27 * C * W / (B*B*R0);  //表面能的常数 
  // eta = 1000;
  cout <<"eta = " << eta << endl;
  ksi = 27 * C * L1/(B*B*R0*R0);
  cout << "L1 = " << ksi << endl;
  Rad = 1/sqrt(ksi);
  cout << "Rad = " << Rad << endl; 
  **/

  
W = 1e-2;
  Rad = 20;
  // eta = 1;
  ksi = 1.0/Rad/Rad;
 //  eta = 1000;
  cout << "Rad = " << Rad << endl;
  R0 = sqrt(27*C*L1/ksi/B/B);
  cout << "R0 = " << R0 << "m" << endl;
 eta = 0.5 * 27 * C * W / (B*B*R0);  //表面能的常数 

  Qscale = - 3.0*sqrt(6)/2 * C / B;
  beig = sqrt(6)/2 * 1.2;
  beig = Qscale/3; 
  // cout << beig << endl;


  /// I = J = K,  L = M = N
  Basis_init(64, 16, 32, 64, 33, 64);
  

  malloc_variables_init();  //申请内存，定义变量

  //  theta_fd = PI/4.0;
  iput(Anlm, landau_t, Rad, eta, 64, 16, 32, argv[1], 1, 2);  ///展开系数等

  cout << "Done!" << endl;


  sprintf(fname,"%s/Eigvalue/%s%s_t_%.2f_R_%.2f_I_%d_J_%d_K_%d_eta_%.1f_L21_%+.1f_s0_%.2f_%s_Eig.txt",DIR,func_type,boundary,landau_t,Rad, I,J,K,eta,L21,beig,argv[1]);

  double Energy;

  Energy = cal_F(Anlm, Qijk, Fbulk, Felas);
  cal_dF(Anlm, Qijk, grad_Energy);
  
  double df = Norm(grad_Energy, 5*Basis);


  double* V0 = new double[5 * Basis](); // g_k_1;
  double* V1 = new double[5 * Basis](); // g_k_1;

  // double* grad_Energy_0 = new double[5 * Basis](); // g_k_1;
  double* Anlm_0 = new double[5 * Basis](); // x_k_1;
  double* Anlm_1 = new double[5 * Basis](); // x_k_1;

  double l = 1e-7;
  double tau = 1e-6;

  srand((unsigned)time(NULL));


  if ((fp = fopen(fname,"r")) == NULL)
  {
    printf("Open Error!\n");
    for(int ix = 0; ix < 5 * Basis; ix++)
    {
      V0[ix] = 2.0*rand()/RAND_MAX - 1;
    }
  }
  else
  {

  int err;

  for (int ix = 0; ix < 5 * Basis; ix++)
  {
    /// fprintf(fp,"%16.15e\n",V0[ix]);
    err = fscanf(fp,"%lf\n",&V0[ix]);
  }
  
  fclose(fp);
  }

  /**

  */

  df = Norm(V0, 5 * Basis);
  for(int ix = 0; ix < 5 * Basis; ix++)
  {
    V0[ix] = V0[ix]/df; // x_1
  }


      
  cout << "Norm(V0) = " << Norm(V0, 5*Basis) << endl;

  double Tmp;

  double* grad_Energy_0 = new double[5 * Basis]();
  double* Gv = new double[5 * Basis](); // x_k_1;
  double* RHS = new double[5 * Basis](); // x_k_1;
  double* RHS_0 = new double[5 * Basis](); // x_k_1;


  double* Sk = new double[5 * Basis](); // x_k_1;
  double* Gk = new double[5 * Basis](); // x_k_1;


  
  double s2;
  double sy;


  

  double Gvv, vv;

  for(int it = 0; it < 10000; it++)  /// it is k
  {
    /// Compute Gv
    /**
    for(int ix = 0; ix < 5 * Basis; ix++)
    {
      Anlm_0[ix] = Anlm[ix] + l * V0[ix]; // x_1
      Anlm_1[ix] = Anlm[ix] - l * V0[ix]; // x_1
    }


    Tmp = cal_F(Anlm_0, Qijk, Fbulk, Felas);
    cal_dF(Anlm_0, Qijk, grad_Energy_0);  /// grad_Energy = k-1 g_1

    Tmp = cal_F(Anlm_1, Qijk, Fbulk, Felas);
    cal_dF(Anlm_1, Qijk, grad_Energy);  /// grad_Energy = k-1 g_1

    for(int ix = 0; ix < 5 * Basis; ix++)
    {
      Gv[ix] =   (grad_Energy_0[ix] - (grad_Energy[ix]))/2.0/l;  // Anlm[ix] + l * V0[ix]; // x_1
    }

   
    vv = Inner(V0, V0, 5 * Basis);
    Gvv = Inner(Gv, V0, 5 * Basis);

    // cout << Gvv << endl;

    for(int ix = 0; ix < 5 * Basis; ix++)
    {
      RHS[ix] = - 2.0/vv * (Gv[ix] - Gvv * V0[ix] / vv);
      //  RHS[ix] = - 2.0/vv * Gv[ix] + 2.0 * Gvv * V0[ix] / vv/vv;
      V1[ix] = RHS[ix] * tau + V0[ix];
    }

    cout << Inner(RHS, V0, 5 * Basis) << endl;

    df =  Norm(RHS, 5*Basis);
    if(abs(df) < 1e-4)
	break;

    cout << "lambda = " << Gvv/vv <<  " ,vv =  " << vv  << ", df = " << df << endl; 

    
    for(int ix = 0; ix < 5 * Basis; ix++)
    {
      V0[ix] = V1[ix]; // x_1
    }
    */
    
    //Barzilai-Borwein GD
    if(it == 0)
    {
      for(int ix = 0; ix < 5 * Basis; ix++)
      {
	Anlm_0[ix] = Anlm[ix] + l * V0[ix]; // x_1
	Anlm_1[ix] = Anlm[ix] - l * V0[ix]; // x_1
      }

      Tmp = cal_F(Anlm_0, Qijk, Fbulk, Felas);
      cal_dF(Anlm_0, Qijk, grad_Energy_0);  /// grad_Energy = k-1 g_1

      Tmp = cal_F(Anlm_1, Qijk, Fbulk, Felas);
      cal_dF(Anlm_1, Qijk, grad_Energy);  /// grad_Energy = k-1 g_1

      for(int ix = 0; ix < 5 * Basis; ix++)
      {
	Gv[ix] =   (grad_Energy_0[ix] - (grad_Energy[ix]))/2.0/l;  // Anlm[ix] + l * V0[ix]; // x_1
      }
      
      vv = Inner(V0, V0, 5 * Basis);
      Gvv = Inner(Gv, V0, 5 * Basis);

      // cout << Gvv << endl;
      for(int ix = 0; ix < 5 * Basis; ix++)
      {
	RHS[ix] = 2.0/vv * (Gv[ix] - Gvv/vv * V0[ix]);
	V1[ix] = - RHS[ix] * tau + V0[ix];
	RHS_0[ix] = RHS[ix];
      }   
    }
    else
    {
      // df = norm(V0, 5 * Basis);
      
      for(int ix = 0; ix < 5 * Basis; ix++)
      {
	Anlm_0[ix] = Anlm[ix] + l * V0[ix]; // x_1
	Anlm_1[ix] = Anlm[ix] - l * V0[ix]; // x_1
      }

      Tmp = cal_F(Anlm_0, Qijk, Fbulk, Felas);
      cal_dF(Anlm_0, Qijk, grad_Energy_0);  /// grad_Energy = k-1 g_1

      Tmp = cal_F(Anlm_1, Qijk, Fbulk, Felas);
      cal_dF(Anlm_1, Qijk, grad_Energy);  /// grad_Energy = k-1 g_1

      for(int ix = 0; ix < 5 * Basis; ix++)
      {
	Gv[ix] = (grad_Energy_0[ix] - (grad_Energy[ix]))/2.0/l;  // Anlm[ix] + l * V0[ix]; // x_1
      }

      // cout <<  Inner(Gv, Gv, 5 * Basis) << endl;
	
      vv = Inner(V0, V0, 5 * Basis);
      Gvv = Inner(Gv, V0, 5 * Basis);

      // cout << Gvv << endl;
      for(int ix = 0; ix < 5 * Basis; ix++)
      {
	RHS[ix] = 2.0/vv * (Gv[ix] - Gvv/vv * V0[ix]);
      }

  
      df =  Norm(RHS, 5*Basis);
      if(abs(df) < 1e-3)
	break;
     
      if(it % 1 == 0)
      {	
	cout << "lambda = " << Gvv/vv << " ,vv =  " << vv  << ", df = " << df << endl;
      }
      
      /// compute g_k
      s2 = 0;
      sy = 0;
      for(int ix = 0; ix < 5 * Basis; ix++)
      {
   
	//	Gk[ix] =  RHS[ix] - RHS_0[ix]; // s_k_1
	// Anlm_0[ix] = Anlm[ix] - Anlm_0[ix]; // y_k_1

	s2 += Sk[ix] * Sk[ix];
	sy += Sk[ix] * (RHS[ix] - RHS_0[ix]);
       
      }


      double a_k = s2/sy;

      if(a_k * df >= 1e-2)
      {
	a_k = 1e-2/df;
	cout << a_k << endl;
      }
 
      for(int ix = 0; ix < 5 * Basis; ix++)
      {
	// 	V1[ix] = - a_k * RHS[ix] + V0[ix]; // x_k
	V1[ix] =  - a_k * RHS[ix] + V0[ix]; // x_k
      }
    }

    df = Norm(V1, 5 * Basis);
    
    for(int ix = 0; ix < 5 * Basis; ix++)
    {
      Sk[ix] = V1[ix] - V0[ix];
      // Gk[ix]
      //  V0[ix] = V1[ix]; // / df; // x_1
      V0[ix] = V1[ix] / df;
      RHS_0[ix] = RHS[ix];
    }


  }
   

  cout << "Lambda = " << Gvv/vv << " ,vv =  " << vv  << endl;

	  
  fp = fopen(fname,"w");

  for (int i = 0; i < 5 * Basis; i++)
  {
    fprintf(fp,"%16.15e\n", V0[i]);
  }
  

  /// L-BFGS
  /**
  lbfgs_parameter_t param;
  int ret = 0;
  
  Energy = cal_F(Anlm, Qijk, Fbulk, Felas);
  printf("Energy=%16.15f,  Fs = %16.15f, Fbulk = %16.15f, Fel = %16.15f, Fs_z = %16.15f \n",Energy, Fs, Fbulk, Felas, Fs_z);  

  lbfgs_parameter_init(&param); 
  param.m = 20;   ///修正Hessen矩阵的步数
  param.epsilon = 1e-4; 
  //  param.ftol = 0.1;
  // sqrt(5*Basis);  // ||a||_2 <= C||a||_{\infty}
  param.max_iterations = 1000;
  //  param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
  param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
  ret = lbfgs(5 * Basis, Anlm, &Energy, evaluate, progress, NULL, &param);   ///BFGS迭代找极值
  // oput(Anlm, landau_t, Rad, eta, N, L, M, argv[1]);
  
  printf("Energy=%16.15f,  Fs = %16.15f, Fbulk = %16.15f, Fel = %16.15f, Fs_z = %16.15f \n",Energy, Fs, Fbulk * Lambda, Felas, Fs_z);  

  oput(Anlm, landau_t, Lambda, eta, L, M, N, filename);
  **/

  

  //  fprintf(fp,"%.2f %.6f, %.10f, %s\n",landau_t, Rad, Energy, filename);

  

  return 0;
}
