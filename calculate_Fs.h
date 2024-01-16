/**
 * @file   calculate_Fs.h
 * @author onedimension <onedimension@onedimension-PC>
 * @date   Mon Nov 17 19:11:08 2014
 * 
 * @brief  计算表面能(surface Energy)
 * 
 * 
 */
#include "global.h"
#include "Basis.h"
#include "initial.h"
#include <cmath>

void calc_fs(double *Qb, double *x_b, double *y_b, double *z_b, double *fs, double *grad_fs, int Ind) /*计算surf energy density*/
{
  // int n, j, k;
  double s0 = sqrt(1.5) * (3.0 + sqrt(9.0 - 8 * landau_t))/4;

  if(landau_t > 9.0/8)
  {
    s0 = 0;
  }

  double lambda = 0.01;

  double q1,q2,q3,q4,q5;
  double qq1,qq2,qq3,qq4,qq5;
  double qp1, qp2, qp3;

  double qqq1, qqq4;

  double phi;

  double dphi_x, dphi_y, dphi_z;
  double Qn;

  double p1, p2, p3, p4, p5, trQb2, trQb3;

  double temp_x, temp_y, temp_z;

  double Wp = 100;
  double Wh = 100;

  double eta2 = 1e-2;
  double et  = sqrt(eta2);
  double gamma = 0; // 5000 / eta;
  double w2 = 100;

  double g1 = 100;

  double b1, b2, b3;

  if(strcmp(boundary, "PH") == 0)  ///average tangency
  {
    if(Ind == 1)  // Out surface
    {
     	
      for(int i = 0; i < I; i++)
      {
		for(int k = 0; k < K; k++)
		{
			int ix = i * K + k;

	  		q1 = Qb[ix * 6 + 0] - s0 * (x_b[ix] * x_b[ix] - 1.0/3);
	  		q2 = Qb[ix * 6 + 1] - s0 * x_b[ix] * y_b[ix];
	  		q3 = Qb[ix * 6 + 2] - s0 * x_b[ix] * z_b[ix];
	  		q4 = Qb[ix * 6 + 3] - s0 * (y_b[ix] * y_b[ix] - 1.0/3);
	  		q5 = Qb[ix * 6 + 4] - s0 * y_b[ix] * z_b[ix];

	  
	 		qq1 = Qb[ix * 6 + 0] + 1.0/3 * s0;
	  		qq2 = Qb[ix * 6 + 1];
	  		qq3 = Qb[ix * 6 + 2];
	  		qq4 = Qb[ix * 6 + 3]  + 1.0/3 * s0;
	  		qq5 = Qb[ix * 6 + 4];

	  		qqq1 = Qb[ix * 6 + 0];
	  		qqq4 = Qb[ix * 6 + 3];

	  		trQb2 = 2*(qqq1*qqq1 + qq2*qq2 + qq3*qq3 + qqq4*qqq4 + qq5*qq5 + qqq1*qqq4);
	
	
	  		phi = Qb[ix * 6 + 5];  // phi = +1 homo, phi = -1 tang

	  		dphi_x = phix[i * J * K + (J - 1) * K + k];
	  		dphi_y = phiy[i * J * K + (J - 1) * K + k];
	  		dphi_z = phiz[i * J * K + (J - 1) * K + k];
	  		// i * J * K + j * K + k
	  

	  		qp1 = (Qb[ix * 6 + 0] + s0/3.0) * x_b[ix] + Qb[ix * 6 + 1] * y_b[ix] + Qb[ix * 6 + 2] * z_b[ix];
	  		qp2 = Qb[ix * 6 + 1] * x_b[ix] + (Qb[ix * 6 + 3] + s0/3.0) * y_b[ix] + Qb[ix * 6 + 4] * z_b[ix];
	  		qp3 = Qb[ix * 6 + 2] * x_b[ix] + Qb[ix * 6 + 4] * y_b[ix] - (Qb[ix * 6 + 0] + Qb[ix * 6 + 3] - s0/3.0) * z_b[ix];

		  // Qn = qq1 * dphi_x * dphi_x + qq4 * dphi_y * dphi_y + (s0 -qq1 - qq4) * dphi_z * dphi_z + 2 * qq2 * dphi_x * dphi_y + 2 * qq3 * dphi_x * dphi_z + 2 * qq5 * dphi_z * dphi_y;

		  b1 = (Qb[ix * 6 + 0] + s0/3.0) * dphi_x + Qb[ix * 6 + 1] * dphi_y + Qb[ix * 6 + 2] * dphi_z;
	  	  b2 = Qb[ix * 6 + 1] * dphi_x + (Qb[ix * 6 + 3] + s0/3.0) * dphi_y + Qb[ix * 6 + 4] * dphi_z;
	  	  b3 = Qb[ix * 6 + 2] * dphi_x + Qb[ix * 6 + 4] * dphi_y - (Qb[ix * 6 + 0] + Qb[ix * 6 + 3] - s0/3.0) * dphi_z;

	  	  Qn = b1 * b1 + b2 * b2 + b3 * b3;



	  
	     double dphi2 = 1/2.0 * (dphi_x * dphi_x + dphi_y * dphi_y + dphi_z * dphi_z);


	  /**
	  fs[ix] =  Wp * (qp1 * qp1 + qp2 * qp2 + qp3 * qp3)  + gamma * Qn + lambda * (dphi2 * et + (phi * phi - 1) * (phi *  phi - 1)/4.0/et) + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * w2;
	  grad_fs[ix] =  Wp * ( 2 * qp1 * x_b[ix] - 2 * qp3 * z_b[ix]) + gamma * ( 2 * b1 * dphi_x  - 2 * b3 * dphi_z) + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq1 + 2 * qqq4) * w2;
	  grad_fs[I * K + ix] =  Wp * (2 * qp1 * y_b[ix] + 2 * qp2 * x_b[ix]) + gamma * (2 * b1 * dphi_y + 2 * b2 * dphi_x)  + 2 * (trQb2 - s0*s0*2/3) * 4 * qq2 * w2;
	  grad_fs[2 * I * K + ix] =  Wp * (2 * qp1 * z_b[ix]  + 2 * qp3 * x_b[ix])  + gamma * (2 * b1 * dphi_z  + 2 * b3 * dphi_x) + 2 * (trQb2 - s0*s0*2/3) * 4 * qq3 * w2;
	  grad_fs[3 * I * K + ix] = Wp * (2 * qp2 * y_b[ix] - 2 * qp3 * z_b[ix]) + gamma * (2 * b2 * dphi_y - 2 * b3 * dphi_z) + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq4 + 2 * qqq1) * w2;
	  grad_fs[4 * I * K + ix] = Wp * (2 * qp2 * z_b[ix] + 2 * qp3 * y_b[ix])  + gamma * (2 * b2 * dphi_z + 2 * b3 * dphi_y)  + 2 * (trQb2 - s0*s0*2/3) * 4 * qq5 * w2;
	  grad_fs[5 * I * K + ix] = (phi*phi - 1) * phi / et * lambda;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);

	  //  grad_fs[5 * I * K + ix] = (phi*phi - 1) * phi / eta2 / eta;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);


	  int index = i * J * K + (J - 1) * K + k;

	  
	  temp_x = (2 * b1 * qq1 + 2 * b2 * qq2 + 2 * b3 * qq3) * (gamma) + dphi_x * (lambda * et );
	  temp_y = (2 * b1 * qq2 + 2 * b2 * qq4 + 2 * b3 * qq5) * (gamma) + dphi_y *  (lambda * et );
	  temp_z = (2 * b1 * qq3 + 2 * b2 * qq5 + 2 * b3 * (s0 - qq1 - qq4) )*(gamma) + dphi_z * (lambda * et);
	  */



	  
	 fs[ix] = Wh * 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4) * (1 + phi)/2.0 + Wp * (qp1 * qp1 + qp2 * qp2 + qp3 * qp3) * ((1 - phi)/2.0 + dphi2)  + gamma * Qn + lambda * (dphi2 * et + (phi * phi - 1) * (phi *  phi - 1)/4.0/et) + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * w2;
	  grad_fs[ix] = Wh * 2 * (2 * q1 + q4) * (1 + phi)/2.0 + Wp * ( 2 * qp1 * x_b[ix] - 2 * qp3 * z_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * ( 2 * b1 * dphi_x  - 2 * b3 * dphi_z) + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq1 + 2 * qqq4) * w2;
	  grad_fs[I * K + ix] = Wh * 2 * (2 * q2) * (1 + phi)/2.0 + Wp * (2 * qp1 * y_b[ix] + 2 * qp2 * x_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * (2 * b1 * dphi_y + 2 * b2 * dphi_x)  + 2 * (trQb2 - s0*s0*2/3) * 4 * qq2 * w2;
	  grad_fs[2 * I * K + ix] = Wh * 2*(2 * q3) * (1 + phi)/2.0 + Wp * (2 * qp1 * z_b[ix]  + 2 * qp3 * x_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * (2 * b1 * dphi_z  + 2 * b3 * dphi_x) + 2 * (trQb2 - s0*s0*2/3) * 4 * qq3 * w2;
	  grad_fs[3 * I * K + ix] = Wh * 2*(2 * q4 + q1) * (1 + phi)/2.0 + Wp * (2 * qp2 * y_b[ix] - 2 * qp3 * z_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * (2 * b2 * dphi_y - 2 * b3 * dphi_z) + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq4 + 2 * qqq1) * w2;
	  grad_fs[4 * I * K + ix] = Wh * 2*(2 * q5) *  (1 + phi)/2.0 + Wp * (2 * qp2 * z_b[ix] + 2 * qp3 * y_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * (2 * b2 * dphi_z + 2 * b3 * dphi_y)  + 2 * (trQb2 - s0*s0*2/3) * 4 * qq5 * w2;
	  grad_fs[5 * I * K + ix] = Wh * 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4)/2.0 - Wp * (qp1 * qp1 + qp2 * qp2 + qp3 * qp3)/2.0 + (phi*phi - 1) * phi / et * lambda;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);

	  //  grad_fs[5 * I * K + ix] = (phi*phi - 1) * phi / eta2 / eta;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);


	  int index = i * J * K + (J - 1) * K + k;

	  
	  temp_x = (2 * b1 * qq1 + 2 * b2 * qq2 + 2 * b3 * qq3) * (gamma) + dphi_x * (lambda * et + Wp * (qp1 * qp1 + qp2 * qp2 + qp3 * qp3));
	  temp_y = (2 * b1 * qq2 + 2 * b2 * qq4 + 2 * b3 * qq5) * (gamma) + dphi_y *  (lambda * et + Wp * (qp1 * qp1 + qp2 * qp2 + qp3 * qp3));
	  temp_z = (2 * b1 * qq3 + 2 * b2 * qq5 + 2 * b3 * (s0 - qq1 - qq4) )*(gamma) + dphi_z * (lambda * et + Wp * (qp1 * qp1 + qp2 * qp2 + qp3 * qp3));


	  /**
	  fs[ix] = Wh * 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4) * (1 + phi)/2.0 + (qp1 * qp1 + qp2 * qp2 + qp3 * qp3) * ((1 - phi)/2.0 + dphi2)  + gamma *Qn * sqrt(eta2) + lambda * (dphi2/eta + (phi * phi - 1) * (phi *  phi - 1)/4.0/eta2/eta) + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * w2;
	  grad_fs[ix] = Wh * 2 * (2 * q1 + q4) * (1 + phi)/2.0 + ( 2 * qp1 * x_b[ix] - 2 * qp3 * z_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * (dphi_x * dphi_x - dphi_z * dphi_z) * sqrt(eta2) + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq1 + 2 * qqq4) * w2;
	  grad_fs[I * K + ix] = Wh * 2 * (2 * q2) * (1 + phi)/2.0 + (2 * qp1 * y_b[ix] + 2 * qp2 * x_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * (2 * dphi_x * dphi_y) * sqrt(eta2)  + 2 * (trQb2 - s0*s0*2/3) * 4 * qq2 * w2;
	  grad_fs[2 * I * K + ix] = Wh * 2*(2 * q3) * (1 + phi)/2.0 + (2 * qp1 * z_b[ix]  + 2 * qp3 * x_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * (2 * dphi_x * dphi_z) * sqrt(eta2) + 2 * (trQb2 - s0*s0*2/3) * 4 * qq3 * w2;
	  grad_fs[3 * I * K + ix] = Wh * 2*(2 * q4 + q1) * (1 + phi)/2.0 + (2 * qp2 * y_b[ix] - 2 * qp3 * z_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * (dphi_y * dphi_y - dphi_z * dphi_z) * sqrt(eta2)  + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq4 + 2 * qqq1) * w2;
	  grad_fs[4 * I * K + ix] = Wh * 2*(2 * q5) *  (1 + phi)/2.0 + (2 * qp2 * z_b[ix] + 2 * qp3 * y_b[ix]) * ((1 - phi)/2.0 + dphi2) + gamma * (2 * dphi_z * dphi_y) * sqrt(eta2)  + 2 * (trQb2 - s0*s0*2/3) * 4 * qq5 * w2;
	  grad_fs[5 * I * K + ix] = Wh * 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4)/2.0 - (qp1 * qp1 + qp2 * qp2 + qp3 * qp3)/2.0 + (phi*phi - 1) * phi / eta2 / eta * lambda;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);

	  //  grad_fs[5 * I * K + ix] = (phi*phi - 1) * phi / eta2 / eta;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);


	  int index = i * J * K + (J - 1) * K + k;

	  
	  temp_x = (2 * qq1 * dphi_x + 2 * qq2 * dphi_y + 2 * qq3 * dphi_z) * (gamma) * sqrt(eta2) + dphi_x * (lambda /  eta + qp1 * qp1 + qp2 * qp2 + qp3 * qp3 );
	  temp_y = (2 * qq4 * dphi_y + 2 * qq2 * dphi_x + 2 * qq5 * dphi_z) * (gamma) * sqrt(eta2) + dphi_y *  (lambda /  eta + qp1 * qp1 + qp2 * qp2 + qp3 * qp3 );
	  temp_z = (2 * (s0 -qq1 - qq4) * dphi_z + 2 * qq5 * dphi_y + 2 * qq3 * dphi_x)*(gamma) * sqrt(eta2) + dphi_z * (lambda /  eta + qp1 * qp1 + qp2 * qp2 + qp3 * qp3 );
	  */
	  
	  
	  
	  // double dphi2 = 1/2.0 * (dphi_x * dphi_x + dphi_y * dphi_y + dphi_z * dphi_z);
	  /**
	  fs[ix] = Wh * 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4) * (1 + phi)/2.0 + (qp1 * qp1 + qp2 * qp2 + qp3 * qp3) * ((1 - phi)/2.0)  + gamma *Qn + lambda * (dphi2/eta + (phi * phi - 1) * (phi *  phi - 1)/4.0/eta2/eta) + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * (w2 + dphi2 * g1);
	  grad_fs[ix] = Wh * 2 * (2 * q1 + q4) * (1 + phi)/2.0 + ( 2 * qp1 * x_b[ix] - 2 * qp3 * z_b[ix]) * ((1 - phi)/2.0) + gamma * (dphi_x * dphi_x - dphi_z * dphi_z)  + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq1 + 2 * qqq4) * (w2 + dphi2 * g1);
	  grad_fs[I * K + ix] = Wh * 2 * (2 * q2) * (1 + phi)/2.0 + (2 * qp1 * y_b[ix] + 2 * qp2 * x_b[ix]) * ((1 - phi)/2.0) + gamma * (2 * dphi_x * dphi_y)   + 2 * (trQb2 - s0*s0*2/3) * 4 * qq2 * (w2 + dphi2 * g1);
	  grad_fs[2 * I * K + ix] = Wh * 2*(2 * q3) * (1 + phi)/2.0 + (2 * qp1 * z_b[ix]  + 2 * qp3 * x_b[ix]) * ((1 - phi)/2.0) + gamma * (2 * dphi_x * dphi_z)  + 2 * (trQb2 - s0*s0*2/3) * 4 * qq3 * (w2 + dphi2 * g1);
	  grad_fs[3 * I * K + ix] = Wh * 2*(2 * q4 + q1) * (1 + phi)/2.0 + (2 * qp2 * y_b[ix] - 2 * qp3 * z_b[ix]) * ((1 - phi)/2.0) + gamma * (dphi_y * dphi_y - dphi_z * dphi_z)  + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq4 + 2 * qqq1) * (w2 + dphi2 * g1);
	  grad_fs[4 * I * K + ix] = Wh * 2*(2 * q5) *  (1 + phi)/2.0 + (2 * qp2 * z_b[ix] + 2 * qp3 * y_b[ix]) * ((1 - phi)/2.0) + gamma * (2 * dphi_z * dphi_y)   + 2 * (trQb2 - s0*s0*2/3) * 4 * qq5 * (w2 + dphi2 * g1);
	  grad_fs[5 * I * K + ix] = Wh * 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4)/2.0 - (qp1 * qp1 + qp2 * qp2 + qp3 * qp3)/2.0 + (phi*phi - 1) * phi / eta2 / eta * lambda;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);

	  //  grad_fs[5 * I * K + ix] = (phi*phi - 1) * phi / eta2 / eta;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);

	  int index = i * J * K + (J - 1) * K + k;

	  
	  temp_x = (2 * qq1 * dphi_x + 2 * qq2 * dphi_y + 2 * qq3 * dphi_z) * (gamma)  + dphi_x * (lambda /  eta + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * g1);
	  temp_y = (2 * qq4 * dphi_y + 2 * qq2 * dphi_x + 2 * qq5 * dphi_z) * (gamma) + dphi_y *  (lambda /  eta + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * g1);
	  temp_z = (2 * (s0 -qq1 - qq4) * dphi_z + 2 * qq5 * dphi_y + 2 * qq3 * dphi_x)*(gamma)  + dphi_z * (lambda / eta + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * g1);
	  */


	  /**
	  fs[ix] = Wp * (qp1 * qp1 + qp2 * qp2 + qp3 * qp3) + (phi * phi - 1) * (phi *  phi - 1)/4.0/eta2/eta + gamma *Qn * sqrt(eta2) + dphi2/eta + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * w2;
	  grad_fs[ix] =  Wp * ( 2 * qp1 * x_b[ix] - 2 * qp3 * z_b[ix]) + gamma * (dphi_x * dphi_x - dphi_z * dphi_z) * sqrt(eta2)   + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq1 + 2 * qqq4) * w2;
	  grad_fs[I * K + ix] = Wp * (2 * qp1 * y_b[ix] + 2 * qp2 * x_b[ix]) + gamma * (2 * dphi_x * dphi_y) * sqrt(eta2) + 2 * (trQb2 - s0*s0*2/3) * 4 * qq2 * w2;
	  grad_fs[2 * I * K + ix] = Wp * (2 * qp1 * z_b[ix]  + 2 * qp3 * x_b[ix]) + gamma * (2 * dphi_x * dphi_z) * sqrt(eta2) + 2 * (trQb2 - s0*s0*2/3) * 4 * qq3 * w2;
	  grad_fs[3 * I * K + ix] = Wp * (2 * qp2 * y_b[ix] - 2 * qp3 * z_b[ix]) + gamma * (dphi_y * dphi_y - dphi_z * dphi_z) * sqrt(eta2) + 2 * (trQb2 - s0*s0*2/3) * (4 * qqq4 + 2 * qqq1) * w2;
	  grad_fs[4 * I * K + ix] = Wp * (2 * qp2 * z_b[ix] + 2 * qp3 * y_b[ix])  + gamma * (2 * dphi_z * dphi_y) * sqrt(eta2) + 2 * (trQb2 - s0*s0*2/3) * 4 * qq5 * w2;
	  grad_fs[5 * I * K + ix] = (phi*phi - 1) * phi / eta2 / eta;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);
	  

	  //  grad_fs[5 * I * K + ix] = (phi*phi - 1) * phi / eta2 / eta;    // + (log(phi) + 1 - log(1 - phi) - 1)/eta2 + gamma * (1 - 2 * phi);

	  int index = i * J * K + (J - 1) * K + k;

	  
	  temp_x = (2 * qq1 * dphi_x + 2 * qq2 * dphi_y + 2 * qq3 * dphi_z) * gamma * sqrt(eta2) + dphi_x / eta;
	  temp_y = (2 * qq4 * dphi_y + 2 * qq2 * dphi_x + 2 * qq5 * dphi_z) * gamma * sqrt(eta2) + dphi_y / eta;
	  temp_z = (2 * (s0 -qq1 - qq4) * dphi_z + 2 * qq5 * dphi_y + 2 * qq3 * dphi_x)* gamma * sqrt(eta2) + dphi_z / eta;
	  */


	  /**
	  temp_x = dphi_x / eta;
	  temp_y = dphi_y / eta;
	  temp_z = dphi_z / eta;
	  **/

	  grad_fs_dr[ix] = (temp_x*drdx[index] + temp_y*drdy[index] + temp_z*drdz[index]);
	  grad_fs_dt[ix] = (temp_x*dtdx[index] + temp_y*dtdy[index] + temp_z*dtdz[index]);
	
	 }
      }
    }
    else
    {
      for(int ix = 0; ix < I * K; ix++)
      {
	
	q1 = (Qb[ix * 6 + 0] + s0/3) * x_b[ix] + Qb[ix * 6 + 1] * y_b[ix] + Qb[ix * 6 + 2] * z_b[ix];
	q2 = Qb[ix * 6 + 1] * x_b[ix] + (Qb[ix * 6 + 3] + s0/3) * y_b[ix] + Qb[ix * 6 + 4] * z_b[ix];
	q3 = Qb[ix * 6 + 2] * x_b[ix] + Qb[ix * 6 + 4] * y_b[ix] - (Qb[ix * 6 + 0] + Qb[ix * 6 + 3] - s0/3) * z_b[ix];

   
	p1 = Qb[ix * 6 + 0];
	p2 = Qb[ix * 6 + 1];
	p3 = Qb[ix * 6 + 2];
	p4 = Qb[ix * 6 + 3];
	p5 = Qb[ix * 6 + 4];
	

	trQb2 = 2*(p1*p1 + p2*p2 + p3*p3 + p4*p4 + p5*p5 + p1*p4);

	double eta_2 = w2;

	fs[ix] = Wp * (q1 * q1 + q2 * q2 + q3 * q3 + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * eta_2);
      
	grad_fs[ix] = Wp * ( 2 * q1 * x_b[ix] - 2 * q3 * z_b[ix] + 2 * (trQb2 - s0*s0*2/3) * (4 * p1 + 2 * p4) * eta_2 ) ;
	grad_fs[I * K + ix] = Wp * (2 * q1 * y_b[ix] + 2 * q2 * x_b[ix] + 2 * (trQb2 - s0*s0*2/3) * 4 * p2 * eta_2 ) ;
	grad_fs[2 * I * K + ix] = Wp * (2 * q1 * z_b[ix]  + 2 * q3 * x_b[ix] + 2 * (trQb2 - s0*s0*2/3) * 4 * p3 * eta_2);
	grad_fs[3 * I * K + ix] = Wp * (2 * q2 * y_b[ix] - 2 * q3 * z_b[ix] + 2 * (trQb2 - s0*s0*2/3) * (4 * p4 + 2 * p1) * eta_2 );
	grad_fs[4 * I * K + ix] = Wp * (2 * q2 * z_b[ix] + 2 * q3 * y_b[ix] + 2 * (trQb2 - s0*s0*2/3) * 4 * p5 * eta_2 );
	grad_fs[5 * I * K + ix] = 0;
	
	/**
	fs[ix] = 0;
	 
	grad_fs[ix] = 0;
	grad_fs[I * K + ix] = 0;
	grad_fs[2 * I * K + ix] = 0;
	grad_fs[3 * I * K + ix] = 0;
	grad_fs[4 * I * K + ix] = 0;
	grad_fs[5 * I * K + ix] = 0;
	*/
       
      }
    }
  }


  else if(strcmp(boundary, "arvtan") == 0)  ///average tangency
    {
  	
      for(int i = 0; i < I; i++)
	{
	  for(int k = 0; k < K; k++)
	    {
	    
	    int ix = i * K + k;

	    q1 = (Qb[ix * 6 + 0] + s0/3) * x_b[ix] + Qb[ix * 6 + 1] * y_b[ix] + Qb[ix * 6 + 2] * z_b[ix];
	    q2 = Qb[ix * 6 + 1] * x_b[ix] + (Qb[ix * 6 + 3] + s0/3) * y_b[ix] + Qb[ix * 6 + 4] * z_b[ix];
	    q3 = Qb[ix * 6 + 2] * x_b[ix] + Qb[ix * 6 + 4] * y_b[ix] - (Qb[ix * 6 + 0] + Qb[ix * 6 + 3] - s0/3) * z_b[ix];
      
      /**
      fs[ix] = 2 * (q1 * q1 + q2 * q2 + q3 * q3);
	
      grad_fs[ix] = 2 * (2 * q1 * x_b[ix] - 2 * q3 * z_b[ix]);
      grad_fs[I * K + ix] = 2 * (2 * q1 * y_b[ix] + 2 * q2 * x_b[ix]);
      grad_fs[2 * I * K + ix] = 2 * (2 * q1 * z_b[ix]  + 2 * q3 * x_b[ix]);
      grad_fs[3 * I * K + ix] = 2 * (2 * q2 * y_b[ix] - 2 * q3 * z_b[ix]);
      grad_fs[4 * I * K + ix] = 2 * (2 * q2 * z_b[ix] + 2 * q3 * y_b[ix]);
      */

	    p1 = Qb[ix * 6 + 0];
	    p2 = Qb[ix * 6 + 1];
	    p3 = Qb[ix * 6 + 2];
	    p4 = Qb[ix * 6 + 3];
	    p5 = Qb[ix * 6 + 4];




      
      //  phi = Qb[ix * 6 + 5];  // phi = +1 homo, phi = -1 tang

	    dphi_x = phix[i * J * K + (J - 1) * K + k];
	    dphi_y = phiy[i * J * K + (J - 1) * K + k];
	    dphi_z = phiz[i * J * K + (J - 1) * K + k];
	    // i * J * K + j * K + k
	    
	    Qn = dphi_x * x_b[ix] + dphi_y * y_b[ix] + dphi_z * z_b[ix];
	    double eta_p = 1;
      
      
	 

	    trQb2 = 2*(p1*p1 + p2*p2 + p3*p3 + p4*p4 + p5*p5 + p1*p4);
	    trQb3 = 3*(2*p2*p3*p5 - (p1*p1-p2*p2+p3*p3)*p4 + p1*(p2*p2-p4*p4-p5*p5));
	    
	    double eta_2 = 0;
      
	    fs[ix] = q1 * q1 + q2 * q2 + q3 * q3  + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * eta_2; // + eta_p * Qn * Qn;
      
	    grad_fs[ix] = ( 2 * q1 * x_b[ix] - 2 * q3 * z_b[ix] + 2 * (trQb2 - s0*s0*2/3) * (4 * p1 + 2 * p4) * eta_2 );
	    grad_fs[I * K + ix] = (2 * q1 * y_b[ix] + 2 * q2 * x_b[ix] + 2 * (trQb2 - s0*s0*2/3) * 4 * p2 * eta_2 );
	    grad_fs[2 * I * K + ix] = (2 * q1 * z_b[ix]  + 2 * q3 * x_b[ix] + 2 * (trQb2 - s0*s0*2/3) * 4 * p3 * eta_2);
	    grad_fs[3 * I * K + ix] = (2 * q2 * y_b[ix] - 2 * q3 * z_b[ix] + 2 * (trQb2 - s0*s0*2/3) * (4 * p4 + 2 * p1) * eta_2 );
	    grad_fs[4 * I * K + ix] = (2 * q2 * z_b[ix] + 2 * q3 * y_b[ix] + 2 * (trQb2 - s0*s0*2/3) * 4 * p5 * eta_2 );



	    
	    double eta_U = 0;

	    // trQ2 = 2*(q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
	    // trQ3 = 3*(2*q2*q3*q5 - (q1*q1-q2*q2+q3*q3)*q4 + q1*(q2*q2-q4*q4-q5*q5));

	    fs[ix] = fs[ix] + eta_U * (trQb2*trQb2*trQb2 - 6*trQb3*trQb3);
	    
	    grad_fs[ix] = grad_fs[ix] + eta_U * (6*trQb2*trQb2*(2*p1+p4) - 36*trQb3*(p2*p2-p4*p4-p5*p5-2*p1*p4));
	    grad_fs[I * K + ix] = grad_fs[I * K + ix] + eta_U * (12*trQb2*trQb2*p2 - 72*trQb3*(p3*p5 + p2*p4 + p1*p2));
	    grad_fs[2 * I * K + ix] = grad_fs[2 * I * K + ix] + eta_U * (12*trQb2*trQb2*p3 - 72*trQb3*(p2*p5 - p3*p4));
	    grad_fs[3 * I * K + ix] = grad_fs[3 * I * K + ix] + eta_U * (6*trQb2*trQb2*(2*p4+p1) - 36*trQb3*(p2*p2-p1*p1-p3*p3-2*p1*p4));
	    grad_fs[4 * I * K + ix] = grad_fs[4 * I * K + ix] + eta_U * (12*trQb2*trQb2*p5 - 72*trQb3*(p2*p3 - p1*p5));
	    


		
	    // grad_fs[5 * I * K + ix] = 2 * Qn ;

	    /**
	    int index = i * J * K + (J - 1) * K + k;
	    
	    temp_x = 2 * eta_p * Qn * x_b[ix];
	    temp_y = 2 * eta_p * Qn * y_b[ix];
	    temp_z = 2 * eta_p * Qn * z_b[ix];

	    grad_fs_dr[ix] = (temp_x*drdx[index] + temp_y*drdy[index] + temp_z*drdz[index]);
	    grad_fs_dt[ix] = (temp_x*dtdx[index] + temp_y*dtdy[index] + temp_z*dtdz[index]);
	    // grad_fs_dp[ix] = (temp_x*dpdx[index] + temp_y*dpdy[index] + temp_z*dpdz[index]);
	    */
	    }

	} 
    }

  

 else if(strcmp(boundary, "mix2") == 0)  ///average tangency
  {
    if(Ind == 2)
    {
      for(int ix = 0; ix < I * K; ix++)
      {
    
	q1 = Qb[ix * 5 + 0] - s0 * (x_b[ix] * x_b[ix] - 1.0/3);
	q2 = Qb[ix * 5 + 1] - s0 * x_b[ix] * y_b[ix];
	q3 = Qb[ix * 5 + 2] - s0 * x_b[ix] * z_b[ix];
	q4 = Qb[ix * 5 + 3] - s0 * (y_b[ix] * y_b[ix] - 1.0/3);
	q5 = Qb[ix * 5 + 4] - s0 * y_b[ix] * z_b[ix];

	fs[ix] = 2 * (q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
    
	grad_fs[ix] = 2 * (2 * q1 + q4);
	grad_fs[I * K + ix] = 2*(2 * q2);
	grad_fs[2 * I * K + ix] = 2*(2 * q3);
	grad_fs[3 * I * K + ix] = 2*(2 * q4 + q1);
	grad_fs[4 * I * K + ix] = 2*(2 * q5);

	/**
	fs[ix] = 0;
    
	grad_fs[ix] = 0;
	grad_fs[I * K + ix] = 0;
	grad_fs[2 * I * K + ix] = 0;
	grad_fs[3 * I * K + ix] = 0;
	grad_fs[4 * I * K + ix] = 0;
	grad_fs[5 * I * K + ix] = 0;
	*/
       

      }
    }
    else
    {
      for(int ix = 0; ix < I * K; ix++)
      {
	q1 = (Qb[ix * 5 + 0] + s0/3) * x_b[ix] + Qb[ix * 5 + 1] * y_b[ix] + Qb[ix * 5 + 2] * z_b[ix];
	q2 = Qb[ix * 5 + 1] * x_b[ix] + (Qb[ix * 5 + 3] + s0/3) * y_b[ix] + Qb[ix * 5 + 4] * z_b[ix];
	q3 = Qb[ix * 5 + 2] * x_b[ix] + Qb[ix * 5 + 4] * y_b[ix] - (Qb[ix * 5 + 0] + Qb[ix * 5 + 3] - s0/3) * z_b[ix];

	fs[ix] = (q1 * q1 + q2 * q2 + q3 * q3); // / (R1 * R1);

	grad_fs[ix] = (2 * q1 * x_b[ix] - 2 * q3 * z_b[ix]); // / (R1 * R1);
	grad_fs[I * K + ix] = (2 * q1 * y_b[ix] + 2 * q2 * x_b[ix]); // / (R1 * R1);
	grad_fs[2 * I * K + ix] = (2 * q1 * z_b[ix]  + 2 * q3 * x_b[ix]); // / (R1 * R1);
	grad_fs[3 * I * K + ix] = (2 * q2 * y_b[ix] - 2 * q3 * z_b[ix]); // / (R1 * R1);
	grad_fs[4 * I * K + ix] = (2 * q2 * z_b[ix] + 2 * q3 * y_b[ix]); // / (R1 * R1);


	double eta_2 = 0;

	fs[ix] = (q1 * q1 + q2 * q2 + q3 * q3 + (trQb2 - s0*s0*2.0/3) * (trQb2 - s0*s0*2.0/3) * eta_2);
      
	grad_fs[ix] = ( 2 * q1 * x_b[ix] - 2 * q3 * z_b[ix] + 2 * (trQb2 - s0*s0*2/3) * (4 * p1 + 2 * p4) * eta_2 ) ;
	grad_fs[I * K + ix] = (2 * q1 * y_b[ix] + 2 * q2 * x_b[ix] + 2 * (trQb2 - s0*s0*2/3) * 4 * p2 * eta_2 ) ;
	grad_fs[2 * I * K + ix] = (2 * q1 * z_b[ix]  + 2 * q3 * x_b[ix] + 2 * (trQb2 - s0*s0*2/3) * 4 * p3 * eta_2);
	grad_fs[3 * I * K + ix] = (2 * q2 * y_b[ix] - 2 * q3 * z_b[ix] + 2 * (trQb2 - s0*s0*2/3) * (4 * p4 + 2 * p1) * eta_2 );
	grad_fs[4 * I * K + ix] = (2 * q2 * z_b[ix] + 2 * q3 * y_b[ix] + 2 * (trQb2 - s0*s0*2/3) * 4 * p5 * eta_2 );
      }
    }
  }

}

double calc_energy_Fs(double *Q)
{
  double Fs = 0;
  double *fs1 = new double[I * K];
  double *fs2 = new double[I * K];
  double *Qb = new double[2 * I * K * 6];
  int ix, ix1;
  
  for(int bi = 0; bi < 2; bi++)
  {
    for(int i = 0; i < I; i++)
    {
      for(int k = 0; k < K; k++)
      {
	for (int n = 0; n < 6; n++)
	{
 	  Qb[bi * 6 * I * K + i * K * 6 + k * 6 + n] = Q[I * J * K + bi * I * K + i * K + k + n * Point];
	}

	//	Qb[bi * 5 * I * K + i * K * 5 + k * 5] += q0_f;
	// Qb[bi * 5 * I * K + i * K * 5 + k * 5 + 3] += q3_f;
      }
    }
  }

  int Ind = 1;
  calc_fs(Qb, xb, yb, zb, fs1, grad_fs1, Ind);
  // calc_fs(Qb, xb, yb, zb, fs2, grad_fs2);
  Ind  = 2;
  calc_fs(Qb + 6 * I * K, xb + I * K, yb + I * K, zb + I * K, fs2, grad_fs2, Ind);

  ix = 0;
  double fs1_tmp = 0, fs2_tmp = 0; 
  for (int i = 0; i < I; i++)
  {
    fs1_tmp = 0;
    fs2_tmp = 0;
    for (int k = 0; k < K; k++)
    {
      fs1_tmp += fs1[i * K + k];
      fs2_tmp += fs2[i * K + k];
    }     
    Fs += fs1_tmp * surface_jacobi[i] * coe_mu[i] + fs2_tmp * surface_jacobi[i + I] * coe_mu[i];
    // Fs += fs1_tmp * surface_jacobi[i] * coe_r[i];
  }

  Fs = Fs * 2 * PI / K;

  // cout << Fs << endl;

  double *fnm1, *fnm2;
  fnm1 = new double[(2 * M - 1) * N - M * (M - 1)]();
  fnm2 = new double[(2 * M - 1) * N - M * (M - 1)]();

  int sig1 = 0, sig2 = 1;

  for(int s = 0; s < 5; s++)
  {
    calc_gnm(grad_fs1 + s * I * K, fnm1, sig1);
    calc_gnm(grad_fs2 + s * I * K, fnm2, sig2);



    // cout << fnm1[47] << " " << fnm2[47] << endl;

    ix = 0;
    ix1 = 0;
    if(s < 5)
    {
      for(int l = 0; l < L; l++)
	{
	  ix1 = 0;
	  for(int m = 1 - M; m < M; m++) 
	    {
	      for(int n = abs(m); n < N; n++)
		{
		  grad_Fs[s * Basis + ix] = fnm1[ix1] * Plp[J * L + l] + fnm2[ix1] * Plp[(J + 1) * L + l];
		  // grad_Fs[s * Basis + ix] = fnm1[ix1] * Plp[J * L + l];
		  ix++;
		  ix1++;
		}
	    }
	}
    }
    else
    {
  
      for(int l = 0; l < 1; l++)
      {
	ix1 = 0;
	for(int m = 1 - M; m < M; m++) 
	  {
	    for(int n = abs(m); n < N; n++)
	      {
		grad_Fs[s * Basis + ix] = 0; // fnm1[ix1] * Plp[J * L + l] + fnm2[ix1] * Plp[(J + 1) * L + l];
		// grad_Fs[s * Basis + ix] = fnm1[ix1] * Plp[J * L + l];
		ix++;
		ix1++;
	      }
	  }
      }
    }
  }

 
  double *dr_fnm1, *dt_fnm1;
  dr_fnm1 = new double[(2 * M - 1) * N - M * (M - 1)]();
  dt_fnm1 = new double[(2 * M - 1) * N - M * (M - 1)]();

  calc_dr_gnm(grad_fs_dr, dr_fnm1, sig1);
  calc_dt_gnm(grad_fs_dt, dt_fnm1, sig1);

  ix = 0;
  ix1 = 0;

  /**
  for(int l = 0; l < L; l++)
  {
    ix1 = 0;
    for(int m = 1 - M; m < M; m++) 
    {
      for(int n = abs(m); n < N; n++)
      {
	// 	cout << dt_fnm1[ix1] << endl;
	grad_Fs[5 * Basis + ix] = grad_Fs[5 * Basis + ix] + dr_fnm1[ix1] * Plp[(J-1) * L + l] + dt_fnm1[ix1] * Plp[(J-1) * L + l];
	ix++;
	ix1++;
      }
    }
  }
  */

  for(int l = 0; l < 1; l++)
  {
    ix1 = 0;
    for(int m = 1 - M; m < M; m++) 
    {
      for(int n = abs(m); n < N; n++)
      {
	// 	cout << dt_fnm1[ix1] << endl;
	grad_Fs[5 * Basis + ix] = 0; // grad_Fs[5 * Basis + ix] + dr_fnm1[ix1] * Plp[(J-1) * L + l] + dt_fnm1[ix1] * Plp[(J-1) * L + l];
	grad_Fs[6 * Basis + ix] = 0; //
	ix++;
	ix1++;
      }
    }
  }

  delete[] dr_fnm1;
  delete[] dt_fnm1;
  

  //  exit(0);

 

  // cout << Fs << endl;

  delete[] Qb;
  delete[] fs1;
  delete[] fs2;
  delete[] fnm1;
  delete[] fnm2;


  
  return Fs;
}

