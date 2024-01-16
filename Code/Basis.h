/**
 * @file   zernike.h
 * @author onedimension <onedimension@onedimension-PC>
 * @date   Tue Dec 16 18:58:29 2014
 * 
 * @brief  
 * 
 * 
 */

#ifndef _BASIS_H
#define _BASIS_H

using namespace std;

#include "global.h"
#include<iostream>
#include<fstream>

void sphere_param_init(int maxn,int n,int l,int m,int i,int j,int k);

double cutoff(double xi, double mu, double R_outer, double R_in);

void Basis_init(int pr_n,int pr_l,int pr_m,int pr_i,int pr_j,int pr_k);

void Basis_destroy(); 

///已知展开系数求函数值
/// fnlm = double[Basis] L * ((2 * M - 1) * N - M * (M - 1))
void calc_fijk(double *fnlm,double *fijk);
void calc_dr_fijk(double *fnlm, double *dr_fijk); 
void calc_dp_fijk(double *fnlm,double *dp_fijk); 
void calc_dt_fijk(double *fnlm,double *dt_fijk); 

double calc_fijk(double *fnlm, int i, int j, int k);

double IntCylind(double *f); 

void calc_fnlm(double *fijk, double *fnlm);
void calc_gnlm(double *fijk, double *fnlm); 
void calc_dr_gnlm(double *fijk, double *dr_fnlm); 
void calc_dp_gnlm(double *fijk, double *dp_fnlm); 
void calc_dt_gnlm(double *fijk, double *dt_fnlm);

void calc_dt_gnm(double *fik, double *fnm, int s);
void calc_dr_gnm(double *fik, double *fnm, int s); 

void calc_gnm(double *fik, double *fnm, int s);

double IntBisphereJacobi(double *f); 

inline double realxi(double xi)
{
  return (xi + 1)/2.0 * (xi1 - xi0) + xi0;
}


#endif
