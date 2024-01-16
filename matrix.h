//#ifndef _ASDFD_INCLUDED_
//#define _ASDFD_INCLUDED_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <time.h>
#include "fftw3.h"
#include <stdlib.h>
using namespace std;

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif
#define PI_HALF 1.5707963267948966192
#define PI2 6.2831853071795864769
#define PI4 12.566370614359172954
#define PIPI 9.8696044010893586188

void vecdisp(double* a, int dim); // 
double vecmax(double*x, int dim);
double vecmin(double*x, int dim);
double vecmean(double*x, int dim);

void vecsum(double* res, int dim, double* x); // res = sum(x);
void vecplus(double* res, int dim, double* x, double* y); // res= x+y
void vecmulplus(double* res, int dim, double* x, double* y , double a); // res= x+a*y
void vecminus(double* res, int dim, double* x, double* y); // res=x-y
void vecdot(double* res, int dim, double* x, double* y); // res= x'*y

void vecmul(double* res, int dim, double* x, double* a); // res = a[0]*x
void vecdotmul(double* res, int dim, double* x, double* a); // res = a.*x
void vecdiv(double* res, int dim, double* x, double* a); // res = x/a[0]
void matvecmul(double* res, int m, int n, double* M, double* x); // res = M*x
void negmatvecmul(double* res, int m, int n, double* M, double* x); // res = -M*x
void initeyemat(double* res, int n); // res = eye(n);
void veccopy(double* res, int dim, double* x); // res = x
void vec2mat(double* res, int dim,double* x,double* y); // res = x*y'

double normvec(double *f, int dim); // ||f||_2
void BFGSrevise(double* H, double* s,double* y, int dim);

// index change
int PNindex2number(int *index, int *N, int D);
void number2PNindex(int *index, int n, int *N, int D);
int PNindex2number_fftw(int *index, int *N, int D);
void number2PNindex_fftw(int *index, int num, int *N, int D);
int PNindex2number_fftw_half(int *index, int *N, int D);
void number2PNindex_fftw_half(int *index, int num, int *N, int D);

// number of elements
int prod(int *N, int D);
int prod_fft(int *N, int D);

//void GradInFourierSpace(fftw_complex *grad_f[], fftw_complex *f, int *N, double *T, int dim);
//void FourierCutter(fftw_complex *complex, int *N, int *n, int dim);	// a simple tool for set zeros

void SpinMatrixCreate(double T[][3], double theta,double gamma,double beta); // spin mattrix T
double dawson(double x);

//#endif

void GaussElimination(double *u, double **A_matrix, double *rhs, int n);
void NumberSwap(double *a, double *b);
void CheckIsHaveSolution(double **A_matrix, int n, int k);
void Pivot(double **A_matrix, double *Rhs, int n, int k);