#include "global.h"
//  #include "Anderson.h"

/**
const int RANDOM_MODE = 0;
const int RADIAL_MODE = 1;
const int SPLIT_MODE = 2;
*/
// 

const char boundary[] = "arvtan";
// const char boundary[] = "vertical";
// const char boundary[] = "mix";
// const char boundary[] = "mix2";
// const char boundary[] = "PH";  // BC with PH
const int uniaxial = 0;
const int anisotropy = 0;

double Qscale;
double beig;

double ksi;   ///无量纲化后的弹性能前的常数L1
double L21;
double eta, eta2 = 10;
double landau_t;    // need 9*(C/B)*(1-2*C/B) <= t <= 9/8
double Rad;
double R1;
double a;
double xi0, xi1;

double theta_fd;

double Fs;
double Fbulk, Felas;

char *fname_arv, *fname_loc;

char *filename;

int MaxN,N,L,M,I,J,K,Basis,Point,innerPoint;

int const_idx;
double *mu, *p, *theta, *lambda;
double *Jacobi;
double *Cutoff;
double *xb,*yb,*zb;
double *coe_mu,*coe_p;
double *Rnmr, *Plp, *Xm;
double *dRnmr, *dPlp, *dXm;

double *Kijm;
/// double *Knl;

double *FDI_p,*FDO_p,*FDI_c,*FDO_c;
fftw_plan FFTp_p,FFTq_p,FFTp_c,FFTq_c;

double *Anlm;
double *Vnlm;

double *Anlm_old;

double *grad_Fb;
double *grad_Fe;
double *grad_Energy;
double *Qijk;
double *fbulk;
double *grad_fbulk;

// double Fs;
// double *fs;
double *grad_fs1;
double *grad_fs2;

double *dr_qijk,*dt_qijk,*dp_qijk;
double *Q1x, *Q1y, *Q1z, *Q2x, *Q2y, *Q2z, *Q3x, *Q3y, *Q3z, *Q4x, *Q4y,*Q4z, *Q5x, *Q5y,*Q5z;
double *QL2,*grad_QL2_dr,*grad_QL2_dt,*grad_QL2_dp;
double *QL1,*grad_QL1_dr,*grad_QL1_dt,*grad_QL1_dp;
double *grad_FeL1,*grad_FeL1_r,*grad_FeL1_t,*grad_FeL1_p;
double *grad_FeL2,*grad_FeL2_r,*grad_FeL2_t,*grad_FeL2_p;

double *phix, *phiy, *phiz;

double *rr_x, *rr_y, *rr_z;

double *grad_Fs;

double *grad_fs_dr;
double *grad_fs_dt;

double *drdx,*drdy,*drdz;
double *dtdx,*dtdy,*dtdz;
double *dpdx,*dpdy,*dpdz;
double *surface_jacobi; //坐标变换后作积分的Jacobi行列式的绝对值

Anderson thisAnderson;

// double *dxdr,*dxdp,*dxdt, *dydr, *dydp, *dydt, *dzdr, *dzdt, *dzdp;
