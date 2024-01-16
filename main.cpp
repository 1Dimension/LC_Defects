/**
 * @file   main.cpp
 1;3201;0c* @author onedimension <onedimension@onedimension-PC>
 * @date   Tue Jun 24 21:18:04 2014
 * @brief 主程序
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
#include "calculate_Fs.h"
#include "calculate_Fel.h"

#include "EigInfoTwo.h"
#include "in_and_out.h"
#include "calc_vector.h"
#include "point.h"
// #include "rotation.h"
// #include "check.h"

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
 * @return
 */
double cal_F(double *A, double *Q, double& Fbulk, double& Felas)
{

  // Fbulk = calc_energy_Fb(A, Q, fbulk);  //根据展开系数计算bulk energy
  Felas = calc_F_el(A, Q); //One-constant approximation
  Fs = calc_energy_Fs(Q);

  //  cout <<  (xi1 - xi0)/2 * (Fbulk + Felas) << endl;

  return (xi1 - xi0)/2 * (Fbulk + Felas) + eta * Fs;

  // return eta * Fs;

}

/*********************************************************************/

void cal_dF(double *A, double *Q,double *grad_Energy)
{
  for(int i = 0; i < 5 * Basis; i++)
  {
    grad_Energy[i] = (xi1 - xi0)/2 * (grad_Fb[i] + grad_FeL1[i]) + eta * grad_Fs[i];
  }





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

  if(k % 10 == 0)
  {
    printf("Iteration %d:  ",k);
    printf("Energy = %16.15f  ",Energy);
  //     printf("normdF = %16.15f  step = %16.15f \n Fs = %16.15f, Fbulk = %16.15f, Fel = %16.15f\n",gnorm,step, Fs, Fbulk, Felas);
    printf("normdF = %16.15f  step = %16.15f Fs = %16.15f\n",gnorm,step, Fs);
  }



  if(k % 200 == 0)
  {

    // sprintf(fname_loc, "%s_T_%d", fname_arv, k);

    // cout << fname_loc << endl;


    oput(Anlm, landau_t, Rad, eta, N, L, M, fname_arv);
    cout << "Update done!" << endl;

    // delete[] filename;
  }


  return 0;
}

/*********************************************************************/
int main(int argc,char *argv[])
{
  ///  cout<<setprecision(20);
  int i, j, k, n, l, s, m, ret = 0, con = 0;
  double Energy,normdF;
  double fr_Rad,fr_t,fr_eta,st_Rad,st_t,st_eta,ed_Rad,ed_t,ed_eta;
  double S0;
  lbfgs_parameter_t param;
  FILE *fp = fopen("Log.txt","a+");


  fname_arv = argv[1];
  fname_loc = new char[30];

  //设定物理参数
  double A = - 0.172 * 1e6, B = -2.12 * 1e6, C = 1.73 * 1e6, L1 = 4 * 1e-11, R0 = 5 * 1e-5, W = 1e-2;
  //  double A, B = - 0.816 * 3 * 1e6, C = 0.45 * 4 * 1e6, L1 = 6 * 1e-12, R0 =  * 1e-7, W = 10;
  // A = 0.9 * B * B / C / 27.0;
  // double D0 = 0.5 * R0;  ///球心间的距离

  L21 = 0;

  R1 = 0.9;
  double ep = 0.05;

  // theta_fd = PI/4.0;

  cout << "R = " <<  R1 << endl;
  // a = sqrt(pow(D/2.0, 2) - 1);         ///Bispherical坐标的系数
  a = sqrt((pow(1.0 - R1*R1 - ep * ep, 2.0) - 4 * ep * ep * R1 * R1))/(2.0 * ep);
  // cout << "D =" << D << endl;
  cout << "a =" << a << endl;
  xi0 = asinh(a);
  xi1 = asinh(a/R1); ////

  /**
  cout << "xi0 = " << xi0 << endl;
  cout << "Xi1 = " << xi1 << endl;
  cout << a * cosh(xi0)/sinh(xi0) - a * cosh(xi1)/sinh(xi1) << endl;
  **/
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


  Rad = 50;
  // eta = 1;
  // W = 1000;
  ksi = 1.0/Rad/Rad;
  cout << ksi << endl;
  // eta = 1000;
  cout << "Rad = " << Rad << endl;
  R0 = sqrt(27*C*L1/ksi/B/B);
  cout << "R0 = " << R0 << "m" << endl;

  cout << "ksi = " <<  ksi << endl;

  eta = 0.5 * 27 * C * W / (B*B*R0);  //表面能的常数

  eta = 1; // 100;


  Qscale = - 3.0*sqrt(6)/2 * C / B;
  beig = sqrt(6)/2 * 1.2;
  beig = Qscale/3;
  // cout << beig << endl;

  S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;


  Basis_init(64, 16, 64, 128, 33, 128);
  malloc_variables_init();  //申请内存，定义变量

  cout << "Done!" << endl;

if(! thisAnderson.initialized)
  {
    thisAnderson.m_max = 30;
    thisAnderson.N = 6 * Basis;
    thisAnderson.Initialize();
  }



  //  theta_fd = PI/4.0;



  iput(Anlm, landau_t, Rad, eta, 64, 16, 64, argv[1], 0, 2);  ///展开系数等


  lbfgs_parameter_init(&param);
  param.m = 50;   ///修正Hessen矩阵的步数
  param.epsilon = 1e-10;
  param.delta = 1e-5;
  // param.ftol = 0.01;
  param.max_iterations = 100000;
  param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
  // param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
  // param.linesearch = LBFGS_LINESEARCH_MORETHUENTE;
  ret = lbfgs(5 * Basis, Anlm, &Energy, evaluate, progress, NULL, &param);   ///BFGS迭代找极值


  oput(Anlm, landau_t, Rad, eta, N, L, M, fname_arv);

  Energy = cal_F(Anlm, Qijk, Fbulk, Felas);  ///计算总能量
  // Energy = Energy * xi0; /// 坐标变换产生的常数
  cal_dF(Anlm, Qijk, grad_Energy);
  printf("R = %.2f t = %.2f Energy = %16.15f normdF = %16.15f ", Rad, landau_t, Energy, Norm(grad_Energy, 6*Basis));
  fprintf(fp,"%.2f %.2f %.6f, %.2f, %.10f, %.10f, %s\n",landau_t, ep, Rad, L21, Energy,  Energy * B * B * R0 * R0 / 27/ C / L1, argv[1]);

  fclose(fp);
  return 0;
}
