#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <cmath>
#include <fftw3.h>

#include "Anderson.h"

#define DIR "./"
// #define PI 3.141592653589793
#define PI 3.1415926535897932384626433832795028841971693993751
#define INF 1e7
#define func_type "poly"

/**
const int RANDOM_MODE = 0;
const int RADIAL_MODE = 1;
const int SPLIT_MODE = 2;
*/

extern const char boundary[];
extern const int uniaxial;  
extern const int anisotropy;

extern double Qscale;
extern double beig;

extern double ksi;   ///无量纲化后的弹性能前的常数L1
extern double L21;
extern double eta,eta2;
extern double landau_t;    // need 9*(C/B)*(1-2*C/B) <= t <= 9/8
extern double Rad;
extern double R1;
extern double a;
extern double xi0, xi1;

extern char* filename;

extern char *fname_arv, *fname_loc;

extern double theta_fd;

extern int MaxN,N,L,M,I,J,K,Basis,Point,innerPoint;

extern int const_idx;
extern double *mu, *p, *theta, *lambda;
extern double *Jacobi;
extern double *Cutoff;
extern double *xb,*yb,*zb;
extern double *coe_mu,*coe_p;
extern double *Rnmr, *Plp, *Xm;
extern double *dRnmr, *dPlp, *dXm;

extern double *Kijm;
/// double *Knl;

extern double *FDI_p,*FDO_p,*FDI_c,*FDO_c;
extern fftw_plan FFTp_p,FFTq_p,FFTp_c,FFTq_c;

extern double *Anlm;
extern double *Vnlm;

extern double *Anlm_old;

extern double *grad_Fb;
extern double *grad_Fe;
extern double *grad_Energy;
extern double *Qijk;
extern double *fbulk;
// double *fbeta;
extern double *grad_fbulk;
// double *grad_fbeta;

extern double Fs;
extern double Fbulk, Felas;

extern double *grad_fs1;
extern double *grad_fs2;

extern double *dr_qijk,*dt_qijk,*dp_qijk;
extern double *Q1x, *Q1y, *Q1z, *Q2x, *Q2y, *Q2z, *Q3x, *Q3y, *Q3z, *Q4x, *Q4y,*Q4z, *Q5x, *Q5y,*Q5z;
extern double *QL2,*grad_QL2_dr,*grad_QL2_dt,*grad_QL2_dp;
extern double *QL1,*grad_QL1_dr,*grad_QL1_dt,*grad_QL1_dp;
extern double *grad_FeL1,*grad_FeL1_r,*grad_FeL1_t,*grad_FeL1_p;
extern double *grad_FeL2,*grad_FeL2_r,*grad_FeL2_t,*grad_FeL2_p;

extern double *phix, *phiy, *phiz;
extern double *rr_x, *rr_y, *rr_z;

extern double *grad_Fs;

extern double *grad_fs_dr;
extern double *grad_fs_dt;

extern double *drdx, *drdy, *drdz;
extern double *dtdx, *dtdy, *dtdz;
extern double *dpdx, *dpdy, *dpdz;
extern double *surface_jacobi;
// double *Jacobi; //坐标变换后作积分的Jacobi行列式的绝对值

// double *dxdr,*dxdp,*dxdt, *dydr, *dydp, *dydt, *dzdr, *dzdt, *dzdp;

//int inx(int i,int n,int l,int m);

extern Anderson thisAnderson;

#endif
