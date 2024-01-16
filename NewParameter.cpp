/**
 * @file   parameter.cpp
 * @author onedimension <onedimension@onedimension-PC>
 * @date   Thu Dec 11 11:42:56 2014
 * 
 * @brief  计算程序需要的一些系数
 *         I need courage to Code!!! We need step by step
 * 
 * 
 */

using namespace std;

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision;

int main(int argc, char* argv[])
{
  cout<<setprecision(20);

  int MaxN = 64;
  int N, M, L;  ///N >= M >= L
  int I, J, K;

  I = 32;
  N = MaxN;
  long double* r;
  r = new long double[I + 1]();

  long double *Nnms;
  Nnms = new long double[N * N * N]();
  
  char Nnmsfile[20], in_filename[20], out_filename1[20], out_filename2[20];

  sprintf(Nnmsfile, "Nnmk_N_%d.txt", N);
  sprintf(in_filename, "roots_I_%d.txt", I); //[0,1]区间的Gauss节点
  sprintf(out_filename1, "Rnm_N_%d_I_%d.txt", N, I); ///存储计算的R_n^(m)(r_i)
  sprintf(out_filename2, "dRnm_N_%d_I_%d.txt", N, I); ///存储计算的R_n^(m)(r_i)

  ifstream input;
  ifstream inNnms;
  ofstream output1;
  ofstream output2;

  inNnms.open(Nnmsfile);
  input.open(in_filename);
  output1.open(out_filename1);
  output2.open(out_filename2);
  output1 << setprecision(20);
  output2 << setprecision(20);

  if(! input)
  {
    cout << "Unable to open roots_I !";
    exit(1);
  }
  
  for(int i = 0; i <= I; i++)
  {
    input >> r[i];
  }

  input.close();

  if(! inNnms)
  {
    cout << "Unable to open Nnms !";
    exit(1);
  }
  
  for(int i = 0; i < N * N * N; i++)
  {
    inNnms >> Nnms[i];
  }

  inNnms.close();

  ///R(n, m)
  long double *R, *dR;
  R = new long double[I + 1]();
  dR = new long double[I + 1]();
  
  int ix = 0;
  for(int n = 0; n < N; n++)
  {
    for(int m = 0; m < N; m++)
    {
      for(int s = 0; s < N; s++)
      {
	for(int i = 0; i <= I; i++)
	{
	  if(Nnms[ix] != 0)
	  {
	    if(r[i] != 0)
	    {
	      long double log_R;
	      if(Nnms[ix] > 0)
	      {
		log_R =  log(Nnms[ix]) + log(r[i]) * (long double)(n - 2.0 * s);
		R[i] += exp(log_R); // * (long double)pow(-1,s)
	      }
	     
	      else
	      {
		log_R =  log(-Nnms[ix]) + log(r[i]) * (long double)(n - 2.0 * s);
		R[i] -= exp(log_R); // * (long double)pow(-1,s)
	      }
	      

	      /**
	      if(fabs(Nnms[ix]) < 1) 
	      {
		
		log_R =  log(Nnms[ix]) + log(r[i]) * (long double)(n - 2.0 * s);
		R[i] += exp(log_R); // * (long double)pow(-1,s)
		
		R[i] += Nnms[ix] * pow(r[i], (long double)(n - 2.0 * s));
	      */

	     
	     
	      
	     
	      //R[i] += Nnms[ix] * pow(r[i], (long double)(n - 2.0 * s));
	      /*
	      if(Nnms[ix] > 0)
	      {
		log_R =  log(Nnms[ix]) + log(r[i]) * (long double)(n - 2.0 * s);
		R[i] += exp(log_R); // * (long double)pow(-1,s)
	      }
	      */
	      /**
	      else
	      {
		log_R =  log(- Nnms[ix]) + log(r[i]) * (long double)(n - 2.0 * s);
		R[i] -= exp(log_R);
	      }
	      */
	
	    }

	    if(n - 2 * s == 0)
	    {
	      dR[i] += 0;
	    }
	    else
	    {
	      if(r[i] != 0)
	      {
		long double log_dR;
		if(Nnms[ix] < 0)
		{
		  log_dR = log(-Nnms[ix]) + log(r[i]) * (long double)(n - 2*s - 1.0) + log((long double)(n - 2*s));
		  dR[i] += exp(log_dR)*-1; 
		}
		else
		{
		  log_dR = log(Nnms[ix]) + log(r[i]) * (long double)(n - 2*s - 1.0) + log((long double)(n - 2*s));
		  dR[i] += exp(log_dR);
		}
	      }
	      else
	      {
		dR[i] += 0;
	      }
	    }
	  }
	}
	ix++;
      }
      for(int i = 0; i <= I; i++)
      {

	output1 << R[i] << endl;
	output2 << dR[i] << endl;
	R[i] = 0;
	dR[i] = 0;
      }
    }
  }

  delete[] R;
  delete[] dR;
  delete[] Nnms;

  output1.close();
  output2.close();

  J = 32;
  L = MaxN;

  long double *Nls;
  Nls = new long double[L * L]();

  char inNlsfile[20];

  sprintf(inNlsfile, "Nls_L_%d.txt", L);
  inNnms.open(inNlsfile);

  sprintf(in_filename, "roots_J_%d.txt", J); //[-1,1]区间的Gauss节点
  input.open(in_filename);
  sprintf(out_filename1, "P_L_%d_J_%d.txt", L, J); ///存储计算的P_l(p_j) 
  output1.open(out_filename1);
  sprintf(out_filename2, "dP_L_%d_J_%d.txt", L, J); ///存储计算的P_l(p_j) 
  output2.open(out_filename2);

  output1 << setprecision(20);
  output2 << setprecision(20);

  long double* p;
  p = new long double[J + 2]();

  if(! input)  ///打开失败
  {
    cout << "Unable to open roots_J !";
    exit(1);
  }
  
  for(int j = 0; j < J + 2; j++)
  {
    input >> p[j];
    cout << p[j] << endl;
  }
    
  input.close();

  if(!inNnms)  ///打开失败
  {
    cout << "Unable to open Nls !";
    exit(1);
  }
  
  for(int i = 0; i < L * L; i++)
  {
    inNnms >> Nls[i];
  }
    
  inNnms.close();


  ///P_l  
  ///如果用迭代？？？/
  long double* P;
  long double* dP;
  P = new long double[J + 2]();
  dP = new long double[J + 2]();

  ix = 0;
  for(int l = 0; l < L; l++)
  {
    for(int s = 0; s < L; s++)
    { 
      if(Nls[ix] != 0)
      {
	for(int j = 0; j < J + 2; j++)
	{
	  if(p[j] != 0)
	  {
	    long double log_P = log(Nls[ix]) + log(fabs(p[j]))*(long double)(l - 2*s);
	    if(p[j] < 0)
	    {
	      P[j] += exp(log_P) * pow(-1,l-2*s);
	    }
	    else
	      P[j] += exp(log_P);
	  }
	  else
	  {
	    P[j] += 0;
	  }
	  if(l - 2 * s != 0)
	  {
	    if(p[j] != 0)
	    {
	      long double log_dP = log(Nls[ix]) + log(fabs(p[j])) * (long double)(l - 2*s - 1) + log((long double)(l - 2*s));
	      if(p[j] < 0)
	      {
		dP[j] += exp(log_dP) * pow(-1,l-2*s-1) ;
	      }
	      else
		dP[j] += exp(log_dP);
	      //  dP[j] += Nls[ix] * powl(p[j], (long double)(l - 2*s - 1)) * (long double)(l - 2*s);
	    }
	    else
	    {
	      dP[j] += 0;
	    }
	  }
	  else
	  {
	    dP[j] += 0;
	  }
	}
      }
      ix++;
    }

    for(int j = 0; j < J + 2; j++)
    {
      output1 << P[j] << endl;
      output2 << dP[j] << endl;
      P[j] = 0;
      dP[j] = 0;
    }
  }

  output1.close();
  output2.close();

  delete[] P;
  delete[] dP;
  delete[] Nls;

  return 0;
}
