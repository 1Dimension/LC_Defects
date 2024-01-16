
using namespace std;

#include <iostream>
#include <fstream>
#include "Basis.h"
#include "global.h"

/** 
 * 给定迭代初值
 * 
 * @param A 展开式系数
 * @param mode 
 */
void initial(double *A,int mode) 
{
  int i,j,k,n,ix;

  double Theta = PI/2.0; // PI/6.0;

  double theta_b, phi_b;

  double x, y, z;
  double xn, yn, zn;

  int Boundary_index = I * J * K;
  // double *s = new double[I+1]();
  double *Q = new double[7 * Point]();
  double S0 = sqrt(1.5)*(3.0 + sqrt(9.0 - 8*landau_t))/4;




  if(landau_t > 9.0/8)
  {  
    S0 = 0;
  }     

  //  cout << "S0 = "<< S0/Qscale * 2.0/3 << endl;
  cout << "S0 = "<< S0 << endl;
  
  switch (mode) 
  {
  case 0:

     
    srand((unsigned)time(NULL));
    for (i = 0; i < 6 * Basis; i++)
      A[i] = 0.2*(2.0*rand()/RAND_MAX - 1);
   

    break;

  case 1:  ///指向全为(0,0,1)

       
    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
	for (k = 0; k < K; k++) 
	{

	  Q[i * J * K + j * K + k] = S0 * (2.0/3.0) ;
	  Q[Point + i * J * K + j * K + k] = 0;
	  Q[2 * Point + i * J * K + j * K + k] = 0;
	  Q[3 * Point + i * J * K + j * K + k] =  S0 * (- 1.0/3.0);
	  Q[4 * Point + i * J * K + j * K + k] = 0;

	  x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
	  y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
	  z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i]))  - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));

	  Q[5 * Point + i * J * K + j * K + k] = 1; // (2.0*rand()/RAND_MAX - 1);
	  Q[6 * Point + i * J * K + j * K + k] = 2.0*rand()/RAND_MAX - 1;   // 1; // (2.0*rand()/RAND_MAX - 1);
	  


	  /**
	  if(z > 0) // && x * x + y * y < 0.1)
	  {
	    //  Q[5 * Point + i * J * K + j * K + k] = 1;
	    Q[5 * Point + i * J * K + j * K + k] = (2.0*rand()/RAND_MAX - 1);
	  }
	  else
	  {
	    Q[5 * Point + i * J * K + j * K + k] = -1;
	  }
	  */
	  

	  
	}
      }
    }



   
    for (n = 0; n < 7; n++)
    {
      calc_fnlm(Q + n * Point, A + n * Basis);
    }

    break;

  case 2:

  
    for (j = 0; j < J; j++)
    { 
      ix = 0;

      for (i = 0; i < I; i++)
      {
	for (k = 0; k < K; k++)
	{

	  x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
	  y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
	  z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i]))  - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));

	  
	  /**
	  xn = cos(Theta) * x - sin(Theta) * z;
	  yn = y;
	  zn = sin(Theta) * x + cos(Theta) * z; 
	  */

	  
	  xn = x;
	  yn = y;
	  zn = z;

	  if(z > 0)
	  {
	    Q[i * J * K + j * K + k] = S0 * (1 * 1 - 1.0/3);
	    Q[Point + i * J * K + j * K + k] = 0;
	    Q[2 * Point + i * J * K + j * K + k] = 0;
	    Q[3 * Point + i * J * K + j * K + k] = S0 * (- 1.0/3);
	    Q[4 * Point + i * J * K + j * K + k] = 0;
	    
	    
	  }
	  else
	  {
	    Q[i * J * K + j * K + k] = S0 * (0 - 1.0/3);
	    Q[Point + i * J * K + j * K + k] = 0;
	    Q[2 * Point + i * J * K + j * K + k] = 0;
	    Q[3 * Point + i * J * K + j * K + k] = S0 * (1 - 1.0/3);
	    Q[4 * Point + i * J * K + j * K + k] = 0;
	   }
	 

	  /**
	  if(sqrt(xn*xn + yn*yn) > 1e-6)
	  {
	    // phi_b = atan2(zn, sqrt(xn*xn + yn*yn));
	    
	    phi_b = atan2(sqrt(xn*xn + yn*yn), zn);

	    theta_b = atan2(yn, xn);
	    

	    // cout << phi_b << " ";

	    
	    double nx = sin(phi_b + PI/2.0)*cos(theta_b), ny = sin(phi_b + PI/2.0) * sin(theta_b), nz = cos(phi_b + PI/2.0);
	  
	    Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
	    Q[Point + i * J * K + j * K + k] = S0 * nx * ny;
	    Q[2 * Point + i * J * K + j * K + k] = S0* nx * nz;
	    Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
	    Q[4 * Point + i * J * K + j * K + k] = S0 * ny*nz;
	  }
	  else
	  {
	    for(int s = 0; s < 5; s++)
	    {
	      Q[s * Point + i * J * K + j * K + k] = 0;
	    }
	  }
	  */

	  
	  ix++; 
	}
      }


    }    
    
    
    
    for (n = 0; n < 6; n++)
    {
      calc_fnlm(Q + n * Point, A + n * Basis);
    }

    break;

  case 3:

    // theta_fd = PI/2;
    /**
    for (i = 0; i < I; i++)
    {
       for (j = 0; j < J; j++)
       {
	 for (k = 0; k < K; k++)
	   {
	     double nx = sin(theta_fd), ny = 0, nz = cos(theta_fd);
	     
	     Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
	     Q[Point + i * J * K + j * K + k] = 0;
	     Q[2 * Point + i * J * K + j * K + k] = nx * nz;
	     Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
	     Q[4 * Point + i * J * K + j * K + k] = 0;
	   }
       }
    }
    */

   for (j = 0; j < J; j++)
    { 
      ix = 0;

      for (i = 0; i < I; i++)
      {
	for (k = 0; k < K; k++)
	{

	  x = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * cos(theta[k]);
	  y = a * sin(mu[i])/(cosh(realxi(p[j])) - cos(mu[i])) * sin(theta[k]);
	  z = a * sinh(realxi(p[j]))/(cosh(realxi(p[j])) - cos(mu[i]))  - a * cosh(realxi(p[j]))/sinh(realxi(p[j]));

	  /**	 
	  xn = cos(Theta) * x - sin(Theta) * z;
	  yn = y;
	  zn = sin(Theta) * x + cos(Theta) * z; 
	  */

	  xn = x;
	  yn = y;
	  zn = z;

	  if(sqrt(xn*xn + zn*zn) > 1e-6)
	  {
	     phi_b = atan2(yn, sqrt(xn*xn + zn*zn));
	    
	   // phi_b = acos(yn);

	    theta_b = atan2(zn, xn);
	    

	    // cout << phi_b << " ";

	    
	    double nx = sin(phi_b + PI/2)*cos(theta_b), ny = sin(phi_b + PI/2) * sin(theta_b), nz = cos(phi_b + PI/2);

	    //  double nx = sin(phi_b + PI/2.0)*cos(theta_b), ny = sin(phi_b + PI/2.0) * sin(theta_b), nz = cos(phi_b + PI/2.0);
	  
	    Q[i * J * K + j * K + k] = S0 * (nx * nx - 1.0/3);
	    Q[Point + i * J * K + j * K + k] = S0 * nx * ny;
	    Q[2 * Point + i * J * K + j * K + k] = S0* nx * nz;
	    Q[3 * Point + i * J * K + j * K + k] = S0 * (ny * ny - 1.0/3);
	    Q[4 * Point + i * J * K + j * K + k] = S0 * ny*nz;
	  }
	  else
	  {
	    for(int s = 0; s < 5; s++)
	    {
	      Q[s * Point + i * J * K + j * K + k] = 0;
	    }
	  }
	  
		  ix++; 
}
 	}	

    } 

    for (n = 0; n < 6; n++)
    {
      calc_fnlm(Q + n * Point, A + n * Basis);
    }

    break;

  case 4:  ///从文件中读取Q_{ijk}

    char filename[40];
    sprintf(filename,"Q_I_%d_J_%d_K_%d.txt", I, J, K);

    ifstream input;
    input.open(filename);

    for (i = 0; i < I; i++) 
    {
      for (j = 0; j < J; j++) 
      {
	for (k = 0; k < K; k++) 
	{
	  for (int s = 0; s < 5; s++)
	  {
	    input >> Qijk[s * Point + i * J * K + j * K + k];
	    //  Qijk[s * Point + i * J * K + j * K + k] -= q_far[s];
	    //  Qijk[s * Point + i * J * K + j * K + k] = Qijk[s * Point + i * J * K + j * K + k]; 
	  }

	  Qijk[5 * Point + i * J * K + j * K + k] = -1;
	}
      }
    }

    input.close();

    for (n = 0; n < 6; n++)
    {
      calc_fnlm(Qijk + n * Point, A + n * Basis);
    }

    break;
  }
  //   delete [] s;
  delete [] Q;
}

/** 
 * 
 * 
 * @param fname 
 * @param t 
 * @param R 
 * @param eta 
 * @param N 
 * @param L 
 * @param M 
 * @param num 
 */
void getfname(char fname[],double t,double R,double eta,int N,int L,int M,char num[]) 
{
  char strR[20];
  char strN[20];
  char strL[20];
  char strM[20];
  
  if (R < 10)
    sprintf(strR,"0%.2f",R);
  else
    sprintf(strR,"%.2f",R);

  if (N < 10)
    sprintf(strN,"0%d",N);
  else
    sprintf(strN,"%d",N);

  if (L < 10)
    sprintf(strL,"0%d",L);
  else
    sprintf(strL,"%d",L);
  
  if (M < 10)
    sprintf(strM,"0%d",M);
  else
    sprintf(strM,"%d",M);
  
  sprintf(fname,"%s/Result/%s%s_t_%+.2f_R_%s_N_%s_L_%s_M_%s_eta_%.2f_L21_%+.1f_s0_%.2f_%s.txt",DIR,func_type,boundary,t,strR,strN,strL,strM,eta,L21,beig,num);

  //  sprintf(fname,"%s/Result/%s%s_t_%+.2f_R_%s_N_%s_L_%s_M_%s_eta_%.0e_%s.txt",DIR,func_type,boundary,t,strR,strN,strL,strM,eta,num);

}

/***************************************************************************/
/** 
 * 
 * 
 * @param A 展开式系数
 * @param t 约化的温度
 * @param R 半径
 * @param eta 参数
 * @param FN 
 * @param FL 
 * @param FM 
 * @param num 
 * @param firstvalue 迭代初值 
 * @param mode 
 */
void iput(double *A, double t, double R,double eta,int FN,int FL,int FM,char num[],int firstvalue,int mode) 
{
  int i,n,l,m,ix,err;
  int Max;
  FILE *fp;
  char path[200];

  if (firstvalue == 0)
  {
    initial(A, mode);
  }

  else  //读取展开式系数 
  {
    
    if (M >= FM)
      Max = M;
    else
      Max = FM;
    
    Max = 128;

    double *DataAnlm = new double[7 * Max * Max * (2 * Max - 1)]();
    getfname(path,t,R,eta,FN,FL,FM,num);

    cout << path << endl;

    if ((fp = fopen(path,"r")) == NULL) printf("Open Anlm Error!\n");
    /**
    ix = 0;
    for(i = 0; i < 5 * Basis; i++)
    {
      err = fscanf(fp,"%lf\n", &A[ix]);
      ix++;
    }
    fclose(fp);
    */

    
    for(i = 0; i < 7; i++)
    {
      for(l = 0; l < FL; l++)
      {
	for(m = 1 - FM; m <= FM - 1; m++)
	{
	  for(n = abs(m); n < FN; n++)
	  {
	    err = fscanf(fp,"%lf\n",&DataAnlm[i * Max * Max * (2 * Max - 1) + l * (2 * Max -1) * Max + (m + Max - 1) * Max + n]);
	  }
	}
      }
    }
	 
   fclose(fp);

    ix = 0;
    for (i = 0; i < 7; i++)
    {
      for(l = 0; l < L; l++)
      {
	for(m = 1 - M; m <= M - 1; m++)
	{
	  for(n = abs(m); n < N; n++)
	  {
	    A[ix++] = DataAnlm[i * Max * Max * (2 * Max - 1) + l * (2 * Max -1) * Max + (m + Max - 1) * Max + n];
	  }
	}
      }
    }

    delete [] DataAnlm;
  }

}

/***************************************************************************/


void oput(double *A,double t,double R,double eta,int FN,int FL,int FM,char num[]) 
{
  int i,n,l,m,ix,err;
  FILE* fp;
  char path[200];
  
  getfname(path,t,R,eta,FN,FL,FM,num);
  fp = fopen(path,"w");  //path = ./Result/......，要建立Result目录

  for (int i = 0; i < 7 * Basis; i++)
  {
    fprintf(fp,"%16.15e\n",A[i]);
  }

  fclose(fp);
}



