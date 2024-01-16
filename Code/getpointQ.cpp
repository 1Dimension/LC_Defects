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
  
  //      
  double A = - 0.172 * 1e6, B = -2.12 * 1e6, C = 1.73*1e6, L1 = 4 * 1e-11, R0 = 5 * 1e-7, W = 1e-2;
  // double A, B = - 0.816 * 3 * 1e6, C = 0.45 * 4 * 1e6, L1 = 6 * 1e-12, R0 = 1 * 1e-7, W = 10;                                                         
  //A = 1.2 * B * B / C / 27.0;  


  R1 = 0.9;
  double ep = 0.05;


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
  eta = 0.5 * 27 * C * W / (B*B*R0);  //      
  eta = 1000;
  cout <<"eta = " << eta << endl;

  /**
  ksi = 27 * C * L1/(B*B*R0*R0);
  cout << "L1 = " << ksi << endl;
  Rad = 1/sqrt(ksi);
  cout << "Rad = " << Rad << endl;
  */

  Rad = 10;
  ksi = 1.0/Rad/Rad;
  // eta = 1000;
  cout << "Rad = " << Rad << endl;
  R0 = sqrt(27*C*L1/ksi/B/B);
  cout << "R0 = " << R0 << "m" << endl;
  
  Qscale = - 3.0 * sqrt(6)/2.0 * C / B;
  beig = Qscale/3.0;

  double S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;

  S0 = S0/Qscale;  //Q = S0(nn - 1/3 I)

  Basis_init(64, 8, 32, 101, 16, 64);
  malloc_variables_init();
  
  // Anlm = new double[5*Basis]();
  // Qijk = new double[5*Point]();
  // fbulk = new double[innerPoint]();

  //  beig = Qscale/3;	  
  iput(Anlm,landau_t,Rad,eta,N,L,M,argv[1],1,0); 
	
  for(int i = 0; i < 5; i++)
  {
    calc_fijk(Anlm + i * Basis, Qijk + i * Point);  //        Q_{ijk}
  }





  // calc_bulkenergy_landau_de(Qijk,fbulk);
  //   double Fb = calc_energy_Fb(Anlm, Qijk, fbulk);
  //   double Fel = calc_F_el(Anlm);
  // cout << Fb <<' ' << Fel << endl;

  // cout << Fel << endl;
  
  char outputQ_filename[30];
  sprintf(outputQ_filename,"Q_I_%d_J_%d_K_%d.txt",I, J, K);
  ofstream outputQ;

  outputQ.open(outputQ_filename);
  outputQ << setprecision(10);
  
  int ix = 0;
  for (int i = 0; i < I; i++) 
  {
    for (int j = 0; j < J; j++) 
    {
      for (int k = 0; k < K; k++) 
      {
	/**
	double rho = a * sin(mu[i])/(cosh(xi0 * p[j]) - cos(mu[i]));
	double z = a * sinh (xi0 * p[j])/(cosh(xi0 * p[j]) - cos(mu[i]));

	double nx, ny, nz; 

	double r1 = pow(rho - D/2.0, 2) + z * z;

	S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;

	if(r1 <= 0.04)
	{
	  double sin_phi =  z/sqrt(r1);
	  double cos_phi = (rho - D/2.0)/sqrt(r1); 

	  nx = sin_phi * cos(theta[k]);
	  ny = sin_phi * sin(theta[k]);
	  nz = cos_phi; 

	  if(r1 <= 1e-4)
	  {
	    S0 = 0;
	  }


	  Qijk[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
	  outputQ << Qijk[i * J * K + j * K + k] << endl;

	  Qijk[1 * innerPoint + i * J * K + j * K + k] = S0 * nx * ny;
	  outputQ << Qijk[1 * innerPoint + i * J * K + j * K + k] << endl;
	  
	  Qijk[2 * innerPoint + i * J * K + j * K + k] = S0 * nx * nz;
	  outputQ << Qijk[2 * innerPoint + i * J * K + j * K + k] << endl;

	  Qijk[3 * innerPoint + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
	  outputQ << Qijk[3 * innerPoint + i * J * K + j * K + k] << endl;

	  Qijk[4 * innerPoint + i * J * K + j * K + k] = S0 * ny * nz;
	  outputQ << Qijk[4 * innerPoint + i * J * K + j * K + k] << endl;
	}

	else
	{
	  for (int n = 0; n < 5; n++)
	  {
	    outputQ << Qijk[n * Point + i * J * K + j * K + k] + q_far[n] << endl;
	  }
	}
	*/

	x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[0]);
	y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[0]);
	z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i]));

	if(z > 0)
	{
    	  for (int n = 0; n < 5; n++)
	  {
	    outputQ << Qijk[n * Point + i * J * K + j * K + k] + q_far[n] << endl;
	  }
	  
	  ix++;
	}
	else
	{
	}





      }
    }
  }

  outputQ.close();
  

  Basis_destroy();
  var_destroy();

  // delete [] Anlm;
  // delete [] Qijk;
  // delete [] fbulk;
  
  return 0;
}
