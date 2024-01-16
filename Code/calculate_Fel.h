
/**
 * @file   calculate_Fel.h
 * @author onedimension <onedimension@onedimension-PC>
 * @date   Mon Nov  3 21:28:52 2014
 * 
 * @brief  
 * 
 * 
 */
#ifndef _CALCULATE_FEL_H
#define _CALCULATE_FEL_H

using namespace std;

#include <iostream>
#include "global.h"
#include "Basis.h"

double calc_F_el(double *qnlm, double *qijk)  
{

  double a,b,c;
  double trQ2,trQ3,trQ4;


  //Dimensionless的LdG泛函的参数
  a = 0.5 * landau_t;
  b = - sqrt(6);
  c = 0.5;

  double lambda1 = 0.0001; // 0.0001; // 0.0001; //  0.0001; // 0.00001; // .0001; // ; // 0.00001;// 0.0001; // .0001; // 0.0001; // 0.0001; // 0.0001; // .0001; // 0.0001; // 0.0001;
  double lambda2 = 0.0001; // .0001; // ; // 0; //.0001.0001; // parameters for smectic phase

  double kappa = 1000;

   
  int i, j, k;
  double intG;
  double temp_x,temp_y,temp_z;

  double b1, b2, b3;

  double c1, c2, c3;

  double q1, q2, q3, q4, q5;

  double qq1, qq4;

  double phi, rr;

  double L23 = 0;
  
  double pQp;
  double layer_angle = 0; //  30 / 180.0 * PI; // 30; // 30 / 180.0 * PI;

  double s0 = sqrt(1.5) * (3.0 + sqrt(9.0 - 8 * landau_t))/4;
  double f_bulkmin = s0 * s0 / 54*(9 * landau_t - 3*sqrt(6)*s0);

  double q = 25;
  
  double gamma_dr = 1;

  if(landau_t > 9.0/8)
  {
    s0 = 0;
  }

  for(i = 0; i < 7; i++)
  {
    calc_fijk(qnlm + i * Basis, qijk + i * Point);
    calc_dr_fijk(qnlm + i * Basis, dr_qijk + i * innerPoint);
    calc_dt_fijk(qnlm + i * Basis, dt_qijk + i * innerPoint);
    calc_dp_fijk(qnlm + i * Basis, dp_qijk + i * innerPoint);
  }

  for(int ix = 0; ix < innerPoint; ix++)
  {
    q1 = qijk[ix];
    q2 = qijk[Point + ix];
    q3 = qijk[2 * Point + ix];
    q4 = qijk[3 * Point + ix];
    q5 = qijk[4 * Point + ix];
                                                                    
    qq1 = q1 + 1/3.0*s0;
    qq4 = q4 + 1/3.0*s0;

    Q1x[ix] = dr_qijk[ix] * drdx[ix] + dt_qijk[ix] * dtdx[ix] + dp_qijk[ix] * dpdx[ix];
    Q1y[ix] = dr_qijk[ix] * drdy[ix] + dt_qijk[ix] * dtdy[ix] + dp_qijk[ix] * dpdy[ix];
    Q1z[ix] = dr_qijk[ix] * drdz[ix] + dt_qijk[ix] * dtdz[ix] + dp_qijk[ix] * dpdz[ix];

    Q2x[ix] = dr_qijk[ix + innerPoint] * drdx[ix] + dt_qijk[ix + innerPoint] * dtdx[ix] + dp_qijk[ix + innerPoint] * dpdx[ix];
    Q2y[ix] = dr_qijk[ix + innerPoint] * drdy[ix] + dt_qijk[ix + innerPoint] * dtdy[ix] + dp_qijk[ix + innerPoint] * dpdy[ix];
    Q2z[ix] = dr_qijk[ix + innerPoint] * drdz[ix] + dt_qijk[ix + innerPoint] * dtdz[ix] + dp_qijk[ix + innerPoint] * dpdz[ix];

    Q3x[ix] = dr_qijk[ix + 2 * innerPoint] * drdx[ix] + dt_qijk[ix + 2 * innerPoint]*dtdx[ix] + dp_qijk[ix + 2 * innerPoint]*dpdx[ix];
    Q3y[ix] = dr_qijk[ix + 2 * innerPoint] * drdy[ix] + dt_qijk[ix + 2 * innerPoint]*dtdy[ix] + dp_qijk[ix + 2 * innerPoint]*dpdy[ix];
    Q3z[ix] = dr_qijk[ix + 2 * innerPoint] * drdz[ix] + dt_qijk[ix + 2 * innerPoint]*dtdz[ix] + dp_qijk[ix + 2 * innerPoint]*dpdz[ix];

    Q4x[ix] = dr_qijk[ix + 3 * innerPoint] * drdx[ix] + dt_qijk[ix + 3 * innerPoint]*dtdx[ix] + dp_qijk[ix + 3 * innerPoint]*dpdx[ix];
    Q4y[ix] = dr_qijk[ix + 3 * innerPoint] * drdy[ix] + dt_qijk[ix + 3 * innerPoint]*dtdy[ix] + dp_qijk[ix + 3 * innerPoint]*dpdy[ix];
    Q4z[ix] = dr_qijk[ix + 3 * innerPoint] * drdz[ix] + dt_qijk[ix + 3 * innerPoint]*dtdz[ix] + dp_qijk[ix + 3 * innerPoint]*dpdz[ix];

    Q5x[ix] = dr_qijk[ix + 4 * innerPoint] * drdx[ix] + dt_qijk[ix + 4 * innerPoint]*dtdx[ix] + dp_qijk[ix + 4 * innerPoint]*dpdx[ix];	
    Q5y[ix] = dr_qijk[ix + 4 * innerPoint] * drdy[ix] + dt_qijk[ix + 4 * innerPoint]*dtdy[ix] + dp_qijk[ix + 4 * innerPoint]*dpdy[ix];
    Q5z[ix] = dr_qijk[ix + 4 * innerPoint] * drdz[ix] + dt_qijk[ix + 4 * innerPoint]*dtdz[ix] + dp_qijk[ix + 4 * innerPoint]*dpdz[ix];

    trQ2 = 2*(q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5 + q1*q4);
    trQ3 = 3*(2*q2*q3*q5 - (q1*q1-q2*q2+q3*q3)*q4 + q1*(q2*q2-q4*q4-q5*q5));
    trQ4 = trQ2*trQ2;

    fbulk[ix] = a*trQ2 + b*trQ3 + c*trQ4 - f_bulkmin;
    
    
    b1 = Q1x[ix] + Q2y[ix] + Q3z[ix];
    b2 = Q2x[ix] + Q4y[ix] + Q5z[ix];
    b3 = Q3x[ix] + Q5y[ix] + (-Q1z[ix]-Q4z[ix]);

    c1 = (q1 + s0/3.0) * b1 + q2 * b2 + q3 * b3;
    c2 = q2 * b1 + (q4 + s0/3.0) * b2 + q5 * b3;
    c3 = q3 * b1 + q5 * b2 + (-q1 - q4 + s0/3.0) * b3;



    // cout << phix[ix]* phix[ix] +  phiy[ix]*phiy[ix] + phiz[ix]* phiz[ix] << endl;

    grad_fbulk[ix] = a * (4*q1+2*q4) + b*3*(q2*q2-q4*q4-q5*q5-2*q1*q4) + c*2*trQ2*(4*q1+2*q4) + 2 * L23 * (c1 * b1 - c3 * b3)  * 0.5 * ksi + 2 * lambda1 * (pQp - s0 * q * q * cos(layer_angle) * cos(layer_angle)) * (phix[ix] * phix[ix] - phiz[ix] * phiz[ix]) * rr * rr;
    grad_fbulk[innerPoint + ix] = a * 4*q2 + b*6*(q3*q5 + q2*q4 + q1*q2) + c*2*trQ2*4*q2  + 2 * L23 * (c1 * b2 + c2 * b1) * 0.5 * ksi + 2 * lambda1 * (pQp - s0 * q * q * cos(layer_angle) * cos(layer_angle)) * (phix[ix] * phiy[ix]*2) * rr * rr;
    grad_fbulk[2 * innerPoint + ix] = a * 4*q3 + b*6*(q2*q5 - q3*q4) + c*2*trQ2*4*q3 + 2 * L23 * (c1 * b3 + c3 * b1)  * 0.5 * ksi + 2 * lambda1 * (pQp - s0 * q * q * cos(layer_angle) * cos(layer_angle)) * (phix[ix] * phiz[ix]*2) * rr * rr;
    grad_fbulk[3 * innerPoint + ix] = a * (4*q4+2*q1) + b*3*(q2*q2-q1*q1-q3*q3-2*q1*q4) + c*2*trQ2*(4*q4+2*q1) + 2 * L23 * (c2 * b2 - c3 * b3)  * 0.5 * ksi + 2 * lambda1 * (pQp - s0 * q * q * cos(layer_angle) * cos(layer_angle)) * (phiy[ix] * phiy[ix] - phiz[ix] * phiz[ix]) * rr * rr;
    grad_fbulk[4 * innerPoint + ix] = a * 4*q5 + b*6*(q2*q3 - q1*q5) + c*2*trQ2*4*q5 + 2 * L23 * (c2 * b3 + c3 * b2)  * 0.5 * ksi + 2 * lambda1 * (pQp - s0 * q * q * cos(layer_angle) * cos(layer_angle)) * (phiz[ix] * phiy[ix]*2) * rr * rr;
    
    
  

    QL2[ix] = 0.5 * ksi * (2 * (Q1x[ix]*Q1x[ix]+Q1y[ix]*Q1y[ix]+Q1z[ix]*Q1z[ix]+Q2x[ix]*Q2x[ix]+Q2y[ix]*Q2y[ix]+Q2z[ix]*Q2z[ix] + Q3x[ix]*Q3x[ix] + Q3y[ix]*Q3y[ix] + Q3z[ix]*Q3z[ix]+Q4x[ix]*Q4x[ix]+Q4y[ix]*Q4y[ix]+Q4z[ix]*Q4z[ix] + Q5x[ix]*Q5x[ix] + Q5y[ix]*Q5y[ix] + Q5z[ix]*Q5z[ix] + Q1x[ix]*Q4x[ix]+Q1y[ix]*Q4y[ix]+Q1z[ix]*Q4z[ix]) + L21 * (b1 * b1 + b2 * b2 + b3 * b3) + L23 * (c1 * c1 + c2 * c2 + c3 * c3)) + fbulk[ix] + lambda1 * rr * rr * (pQp - s0 * q * q * cos(layer_angle) * cos(layer_angle)) * (pQp - s0 * q * q * cos(layer_angle) * cos(layer_angle)) + lambda2 * rr * rr * ( phix[ix]* phix[ix] +  phiy[ix]*phiy[ix] + phiz[ix]* phiz[ix] - q * q * cos(layer_angle) * cos(layer_angle)) * ( phix[ix]* phix[ix] +  phiy[ix]*phiy[ix] + phiz[ix]* phiz[ix] - q * q * cos(layer_angle) * cos(layer_angle)) + lambda1 * gamma_dr * (rr_x[ix] * rr_x[ix] + rr_y[ix] * rr_y[ix] + rr_z[ix] * rr_z[ix]);




   
    temp_x = 2 * (2 * Q1x[ix] + Q4x[ix]) + L21 * (2 * b1) + L23 * 2 * (c1 * (q1 + s0/3.0) + c2 * q2 + c3 * q3); // \pp Q1x appear in b1
    temp_y = 2 * (2 * Q1y[ix] + Q4y[ix]);
    temp_z = 2 * (2 * Q1z[ix] + Q4z[ix]) + L21 * 2 * b3 * (-1.0) + (-2.0) * L23 * (c1 * (q3) + c2 * (q5) + c3 * (-q1 - q4 + s0/3.0)); // \pp Q1z appears in b3

    grad_QL1_dr[ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];
    
    temp_x = 4 * Q2x[ix] + L21 * 2 * b2 + L23 * 2 * (c1 * q2 + c2 * (q4 + s0/3.0) + c3 * (q5));  // Q2x in b2
    temp_y = 4 * Q2y[ix] + L21 * 2 * b1 + L23 * 2 * (c1 * (q1 + s0/3.0) + c2 * q2 + c3 * q3);     // Q2y in b1
    temp_z = 4 * Q2z[ix];

    grad_QL1_dr[innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];
	
    temp_x = 4 * Q3x[ix] + L21 * 2 * b3 + L23 * 2 * (c1 * (q3) + c2 * (q5) + c3 * (-q1 - q4 + s0/3.0));   // Q3x in b3
    temp_y = 4 * Q3y[ix];
    temp_z = 4 * Q3z[ix] + L21 * 2 * b1 + L23 * 2 * (c1 * (q1 + s0/3.0) + c2 * q2 + c3 * q3);  //  Q3z in b1

    grad_QL1_dr[2 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[2 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[2 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];

    temp_x = 2 * (2 * Q4x[ix] + Q1x[ix]);
    temp_y = 2 * (2 * Q4y[ix] + Q1y[ix]) + L21 * 2 * b2  + L23 * 2 * (c1 * q2 + c2 * (q4 + s0/3.0) + c3 * (q5));   // Q4y in b2
    temp_z = 2 * (2 * Q4z[ix] + Q1z[ix]) + L21 * 2 * b3 * (-1.0) + (-2.0) * L23 * (c1 * (q3) + c2 * (q5) + c3 * (-q1 - q4 + s0/3.0));  // Q4z in b3

    grad_QL1_dr[3 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[3 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[3 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];

    temp_x = 4 * Q5x[ix]; 
    temp_y = 4 * Q5y[ix] + L21 * 2 * b3 + L23 * 2 * (c1 * (q3) + c2 * (q5) + c3 * (-q1 - q4 + s0/3.0));    // Q5y in b3
    temp_z = 4 * Q5z[ix] + L21 * 2 * b2 + L23 * 2 * (c1 * q2 + c2 * (q4 + s0/3.0) + c3 * (q5));  // Q5z in b2

    grad_QL1_dr[4 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[4 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[4 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];

    temp_x = 2 * lambda1 * rr * rr * (pQp - s0 * q * q * cos(layer_angle) * cos(layer_angle)) * (2 * qq1 * phix[ix] + 2 * q2 * phiy[ix] + 2 * q3 * phiz[ix]) +  2 * lambda2 * rr * rr * (phix[ix]* phix[ix] +  phiy[ix]*phiy[ix] + phiz[ix]* phiz[ix] - q * q * cos(layer_angle) * cos(layer_angle)) * 2 * phix[ix];
    temp_y = 2 * lambda1 * rr * rr * (pQp - q * q * s0 * cos(layer_angle) * cos(layer_angle)) * (2 * qq4 * phiy[ix] + 2 * q2 * phix[ix] + 2 * q5 * phiz[ix]) +  2 * lambda2 * rr * rr * (phix[ix]* phix[ix] +  phiy[ix]*phiy[ix] + phiz[ix]* phiz[ix] - q * q * cos(layer_angle) * cos(layer_angle)) * 2 * phiy[ix]; 
    temp_z = 2 * lambda1 * rr * rr * (pQp - s0 * q * q * cos(layer_angle) * cos(layer_angle)) * (2 * (-q1 - q4 + 1/3.0 * s0) * phiz[ix] + 2 * q3 * phix[ix] + 2 * q5 * phiy[ix]) + 2 * lambda2 * rr * rr * ( phix[ix]* phix[ix] +  phiy[ix]*phiy[ix] + phiz[ix]* phiz[ix] - q * q * cos(layer_angle) * cos(layer_angle)) * 2 * phiz[ix];

    grad_QL1_dr[5 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[5 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[5 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];

    temp_x = rr_x[ix] * lambda1 * 2 * gamma_dr; 
    temp_y = rr_y[ix] * lambda1 * 2 * gamma_dr;
    temp_z = rr_z[ix] * lambda1 * 2 * gamma_dr;

    grad_QL1_dr[6 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[6 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[6 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];


    for(int s = 0; s < 5; s++)
    {
      grad_QL1_dr[s * innerPoint + ix] = grad_QL1_dr[s * innerPoint + ix] * 0.5 * ksi;
      grad_QL1_dt[s * innerPoint + ix] = grad_QL1_dt[s * innerPoint + ix] * 0.5 * ksi;
      grad_QL1_dp[s * innerPoint + ix] = grad_QL1_dp[s * innerPoint + ix] * 0.5 * ksi;
    }


    /**
    grad_QL1_dr[5 * innerPoint + ix] =  grad_QL1_dr[5 * innerPoint + ix];
    grad_QL1_dt[5 * innerPoint + ix] =  grad_QL1_dt[5 * innerPoint + ix];
    grad_QL1_dp[5 * innerPoint + ix] =  grad_QL1_dp[5 * innerPoint + ix];
    */

    /**
    QL2[ix] = 2 * 0.5 * (Q1x[ix]*Q1x[ix]+Q1y[ix]*Q1y[ix]+Q1z[ix]*Q1z[ix]+Q2x[ix]*Q2x[ix]+Q2y[ix]*Q2y[ix]+Q2z[ix]*Q2z[ix]
            +Q3x[ix]*Q3x[ix]+Q3y[ix]*Q3y[ix]+Q3z[ix]*Q3z[ix]+Q4x[ix]*Q4x[ix]+Q4y[ix]*Q4y[ix]+Q4z[ix]*Q4z[ix]
            +Q5x[ix]*Q5x[ix]+Q5y[ix]*Q5y[ix]+Q5z[ix]*Q5z[ix]
		 +Q1x[ix]*Q4x[ix]+Q1y[ix]*Q4y[ix]+Q1z[ix]*Q4z[ix]);
    */
    
    // cout << QL2[ix] << " " ;
  }

  // cout << endl;  exit(0);

  intG = IntBisphereJacobi(QL2);

  /**
  for (int ix = 0; ix < innerPoint; ix++) 
  {

    temp_x=2*(2*Q1x[ix]+Q4x[ix]);
    temp_y=2*(2*Q1y[ix]+Q4y[ix]);
    temp_z=2*(2*Q1z[ix]+Q4z[ix]);

    grad_QL1_dr[ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];
    
    temp_x=4*Q2x[ix];
    temp_y=4*Q2y[ix];
    temp_z=4*Q2z[ix];

    grad_QL1_dr[innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];
	
    temp_x=4*Q3x[ix];
    temp_y=4*Q3y[ix];
    temp_z=4*Q3z[ix];

    grad_QL1_dr[2 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[2 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[2 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];

    temp_x=2*(2*Q4x[ix]+Q1x[ix]);
    temp_y=2*(2*Q4y[ix]+Q1y[ix]);
    temp_z=2*(2*Q4z[ix]+Q1z[ix]);

    grad_QL1_dr[3 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[3 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[3 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];

    temp_x=4*Q5x[ix];
    temp_y=4*Q5y[ix];
    temp_z=4*Q5z[ix];

    grad_QL1_dr[4 * innerPoint + ix] = temp_x*drdx[ix] + temp_y*drdy[ix] + temp_z*drdz[ix];
    grad_QL1_dt[4 * innerPoint + ix] = temp_x*dtdx[ix] + temp_y*dtdy[ix] + temp_z*dtdz[ix];
    grad_QL1_dp[4 * innerPoint + ix] = temp_x*dpdx[ix] + temp_y*dpdy[ix] + temp_z*dpdz[ix];


    for(int s = 0; s < 5; s++)
    {
      grad_QL1_dr[s * innerPoint + ix] = grad_QL1_dr[s * innerPoint + ix] * 0.5;
      grad_QL1_dt[s * innerPoint + ix] = grad_QL1_dt[s * innerPoint + ix] * 0.5;
      grad_QL1_dp[s * innerPoint + ix] = grad_QL1_dp[s * innerPoint + ix] * 0.5;
    }

  }
  */

  for (int i = 0; i < 7; i++) 
  {
    calc_gnlm(grad_fbulk + i * innerPoint, grad_Fb + i * Basis);
    calc_dr_gnlm(grad_QL1_dr + i * innerPoint, grad_FeL1_r + i * Basis);
    calc_dt_gnlm(grad_QL1_dt + i * innerPoint, grad_FeL1_t + i * Basis);
    calc_dp_gnlm(grad_QL1_dp + i * innerPoint, grad_FeL1_p + i * Basis);
    //  calc_gnlm(grad_fbulk + i * innerPoint, grad_Fb + i * Basis);  ///计算BULK ENERGY对展开式系数的导数 
  }

  for (int i = 0; i < 7 * Basis; i++) 
  {
    grad_FeL1[i] = grad_FeL1_r[i] + grad_FeL1_t[i] + grad_FeL1_p[i];
  }

  // exit();

  return intG;

}


#endif
