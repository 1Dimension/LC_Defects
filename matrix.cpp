#include "matrix.h"

void vecdisp(double* a, int dim) // 
{
	int i=0;
	for (;i<dim;i++)
		cout<<a[i]<<endl;
	cout<<endl;
}

double vecmax(double*x, int dim)
{
	double m = x[0];
	for(int i=0;i<dim;i++)
		m = m>x[i] ? m:x[i];
	return m;
}

double vecmin(double*x, int dim)
{
	double m = x[0];
	for(int i=0;i<dim;i++)
		m = m<x[i] ? m:x[i];
	return m;
}

double vecmean(double*x, int dim)
{
	double m = 0;
	for(int i=0;i<dim;i++)
		m += x[i];
	return m/dim;
}

void vecsum(double* res, int dim, double* x) // res = sum(x)
{
	int i=0;
	res[0]=0;
	for (;i<dim;i++)
		res[0]=res[0]+x[i];	
}


void vecplus(double* res, int dim, double* x, double* y) // res= x+y
{	
	int i=0;
	for (;i<dim;i++)
		res[i] = x[i]+y[i];
		
}

void vecmulplus(double* res, int dim, double* x, double* y , double a) // res= x+a*y
{
    int i=0;
	for (;i<dim;i++)
		res[i] = x[i]+a*y[i];
}


void vecminus(double* res, int dim, double* x, double* y) // res=x-y
{	
	int i=0;
	for (;i<dim;i++)
		res[i] = x[i]-y[i];
		
}

void vecdot(double* res, int dim, double* x, double* y) // res= x'*y
{	
	int i=0; res[0]=0;
	for (;i<dim;i++)
		res[0] += x[i]*y[i];
		
}


void vecmul(double* res, int dim, double* x, double* a) // res = a[0]*x
{	
	int i=0; 
	for (;i<dim;i++)
		res[i] = x[i]*a[0];
		
}

void vecdotmul(double* res, int dim, double* x, double* a) // res = a.*x
{	
	int i=0; 
	for (;i<dim;i++)
		res[i] = x[i]*a[i];
		
}

void vecdiv(double* res, int dim, double* x, double* a) // res = x/a[0]
{	
	int i=0; 
	for (;i<dim;i++)
		res[i] = x[i]/a[0];
		
}

void matvecmul(double* res, int m, int n, double* M, double* x) // res = M*x
{
	int i,j;
	double* p = M;
	for (i=0;i<m;i++)
	{
		res[i]=0;
		for (j=0;j<n;j++)
		{
			res[i] += (p[0]*x[j]);
			p++;
		}
	}
}

void negmatvecmul(double* res, int m, int n, double* M, double* x) // res = -M*x
{
	int i,j;
	double* p = M;
	for (i=0;i<m;i++)
	{
		res[i]=0;
		for (j=0;j<n;j++)
		{
			res[i] -= (p[0]*x[j]);
			p++;
		}
	}
}

void initeyemat(double* res, int n) // res = eye(n)
{
	int i = 0;
	for(;i<n*n;i++)
		res[i] = (i%(n+1) == 0) ? 1 : 0;	// ????n+1???	
}

void veccopy(double* res, int dim, double* x) // res = x
{
	int i = 0;
	for(;i<dim;i++)
		res[i] = x[i];
}

void vec2mat(double* res, int dim,double* x,double* y) // res = x*y'
{
	int i,j;
	for(i=0;i<dim;i++)
	{
		for(j=0;j<dim;j++)
		{
			res[i*dim+j]=x[i]*y[j];
		}
	}
}

double normvec(double *f, int dim){ // ||f||_2
    double norm=0;
    for(int i=0;i<dim;++i){
        norm+=f[i]*f[i];
    }
    return sqrt(norm);
}

void BFGSrevise(double* H, double* s,double* y, int dim){
    double ys=0;
    for(int i=0;i<dim;++i){
        ys += y[i]*s[i];
    }
    
    double yHy=0;
    for(int i=0;i<dim;++i){
        for(int j=0;j<dim;++j){
            yHy += y[i]*H[i*dim+j]*y[j];
        }
    }
    
    double *Hy = (double*)malloc(sizeof(double)*dim);
    for(int i=0; i<dim; ++i){
        Hy[i]=0;
        for(int j=0;j<dim;++j){
            Hy[i] += H[i*dim+j]*y[j];
        }
    }
    
    for(int i=0;i<dim;++i){
        for(int j=0; j<dim; ++j){
            H[i*dim+j] += s[i]*s[j]/ys*(1+yHy/ys) - (s[i]*Hy[j]+Hy[i]*s[j])/ys;
        }
    }
	free(Hy);
}

// used in the FFT-like arrangement --- matlab function
int PNindex2number(int *index, int *N, int D)
{
	for(int i=0; i<D; i++)
		if(index[i]<0)
			index[i] += N[i];
	int num = index[D-1];
	for(int i=D-2; i>=0; i--)
		num = num*N[i]+index[i];
	return num;
}
void number2PNindex(int *index, int num, int *N, int D)
{
	for(int i=0; i<D; i++)
	{
		index[i] = num % N[i];
		num = (num-index[i])/N[i];
		if(index[i] > N[i]/2)
			index[i] -= N[i];
	}
}
// used in the FFT-like arrangement --- FFTW function
int PNindex2number_fftw(int *index, int *N, int D)
{
	for(int i=0; i<D; i++)
		if(index[i]<0)
			index[i] += N[i];  // 0 <= index[i] < N[i]

	int num = index[0];
	for(int i=1; i<D; i++)
		num = num*N[i]+index[i];
	return num;
}
void number2PNindex_fftw(int *index, int num, int *N, int D)
{
	for(int i=D-1; i>=0; i--)
	{
		index[i] = num % N[i];
		num = (num-index[i])/N[i];
		if(index[i] > N[i]/2)
			index[i] -= N[i];
	}
}
// used in the FFT-like arrangement --- FFTW c2r/r2c function
int PNindex2number_fftw_half(int *index, int *N, int D)
{
	if(index[D-1]<0){
		printf("[PNindex2number_fftw_half] error index.\n");
		return -1;
	}

	for(int i=0; i<D; i++)
		if(index[i]<0)
			index[i] += N[i];

	int n = N[D-1];
	N[D-1] = n/2+1;	// that is the difference 

	int num = index[0];
	for(int i=1; i<D; i++)
		num = num*N[i]+index[i];

	N[D-1] = n;
	return num;
}
void number2PNindex_fftw_half(int *index, int num, int *N, int D)
{
	int n = N[D-1];
	N[D-1] = n/2+1;	// that is the difference
	for(int i=D-1; i>=0; i--)
	{
		index[i] = num % N[i];
		num = (num-index[i])/N[i];
	}
	N[D-1] = n;
	for(int i=D-1; i>=0; i--)
	{
		if(index[i] > N[i]/2)
			index[i] -= N[i];
	}
}
// simple function 
int prod(int *N, int D){
	int num = 1;
	for(int i=0; i<D; i++)
		num *= N[i];
	return num;
}
int prod_fft(int *N, int D){
	int num = 1;
	for(int i=0; i<D-1; i++)
		num *= N[i];
	num *= N[D-1]/2+1;
	return num;
}

void SpinMatrixCreate(double T[][3], double theta,double gamma,double beta)
{
	double ct = cos(theta), st = sin(theta);
	double cg = cos(gamma), sg = sin(gamma);
	double cb = cos(beta), sb = sin(beta);
	T[0][0] = ct; 
	T[0][1] = -st*cg;
	T[0][2] = st*sg;

	T[1][0] = st*cb;
	T[1][1] = ct*cb*cg-sb*sg;
	T[1][2] = -ct*cb*sg-sb*cg;

	T[2][0] = st*sb;
	T[2][1] = ct*sb*cg+cb*sg;
	T[2][2] = -ct*sb*sg+cb*cg;

	/*for(int i=0; i<3; i++)
		for(int j=0; j<3; j++){
			cout<<i<<" "<<j<<" : ";
			double a=0;
			for(int m=0;m<3;m++)
				a += T[m][i]*T[m][j];
			cout<<a<<endl;
		}
	getchar();	// check T*/ 


} // spin mattrix T

/*void GradInFourierSpace(fftw_complex *grad_f[], fftw_complex *f, int *N, double *T, int dim) // 梯度
{
	int ND = prod(N,dim);
	int *index = (int*)malloc(sizeof(int)*dim);
	for(int i=0; i<ND; i++)	{
		number2PNindex_fftw_half(index, i, N,dim);	//number2PNindex, number2PNindex_fftw,number2PNindex_fftw_half
		for(int j=0; j<dim; j++){
			double k = index[j];
			double t = T[j];
			double a = 2*PI*k/t;
			grad_f[j][i][0] = -a*f[i][1];
			grad_f[j][i][1] = a*f[i][0];
		}
	}
	free(index);
}*/
/*
void FourierCutter(fftw_complex *complex, int *N, int *n, int dim)
{
	int ND = prod(N,dim);
	int *index = (int*)malloc(sizeof(int)*dim);
	double maxCut = 0;
	for(int i=0; i<ND; i++){  // (warning)
		number2PNindex_fftw_half(index, i, N, dim);	//number2PNindex, number2PNindex_fftw,number2PNindex_fftw_half
		for(int d=0; d<dim; d++)
			if( index[d] > n[d]/2 || index[d]<-n[d]/2){ // check if need
				maxCut = max( maxCut, complex[i][0]*complex[i][0] + complex[i][1]*complex[i][1] );
				complex[i][0] = 0;
				complex[i][1] = 0;
				break;
			}
	}
	if(maxCut>1e-3)
		printf("\n[FourierCutter] warning : maxCut = %g\n", maxCut);
	free(index);
}//*/
/*
void Fourier2Real(fftw_complex *complex, double *real, int *N, int dim)
{
//	cout<<"in Fourier2Real"<<endl;
	fftw_plan p = fftw_plan_dft_c2r(dim, N, complex, real, FFTW_ESTIMATE);
	fftw_execute(p); 
	fftw_destroy_plan(p);

	int ND = prod(N,dim);
	for(int i=0; i<ND; i++)
		;//real[i] /= N;  // c2r 不要除
//	cout<<"out Fourier2Real"<<endl;
}
void Real2Fourier(double *real, fftw_complex *complex, int *N, int *n, int dim) 
{
	fftw_plan p = fftw_plan_dft_r2c(dim, N, real, complex, FFTW_ESTIMATE);
	fftw_execute(p); 
	fftw_destroy_plan(p);

	int ND = prod(N,dim);
	for(int i=0; i<ND; i++){
		complex[i][0] /= ND;	// r2c 要除N
		complex[i][1] /= ND;
	}

	if(n[0]>0 || n[0]<N[0]){  // k=0: do not cut-off
		FourierCutter(complex,N,n,dim);	// de-aliasing
		getchar();
	}
}*/

//integration \int_0^x e^(t^2-x^2) dt 
double dawson(double x) 
{
	const int NMAX = 6;
	const double H = 0.4;
	const double A1 = 2.0/3.0;
	const double A2 = 0.4;
	const double A3 = 2.0/7.0;
	
	int i,n0;
	double d1,d2,e1,e2,sum,x2,xp,xx,ans;
	static double c[NMAX+1];
	static int init = 0; 
	if (init == 0) {
		init=1;
		for (i=1;i<=NMAX;i++) 
			c[i]=exp(-((2.0*i-1.0)*H)*((2.0*i-1.0)*H));
	}
	if (fabs(x) < 0.2) { 
		x2=x*x;
		ans=x*(1.0-A1*x2*(1.0-A2*x2*(1.0-A3*x2)));
	} 
	else { 
		xx=fabs(x);
		n0=2*(int)(0.5*xx/H+0.5);
		xp=xx-n0*H;
		e1=exp(2.0*xp*H);
		e2=e1*e1;
		d1=n0+1;
		d2=d1-2.0;
		sum=0.0;
		for (i=1;i<=NMAX;i++,d1+=2.0,d2-=2.0,e1*=e2)
			sum += c[i]*(e1/d1+1.0/(d2*e1));
		ans=1./sqrt(PI)*(x>0.0? fabs(exp(-xp*xp)):-fabs(exp(-xp*xp)))*sum; 
	}
	return ans;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GaussElimination(double *u, double **A_matrix, double *rhs, int n)
{
	for(int s=0; s<n-1; s++)
	{
		CheckIsHaveSolution(A_matrix, n, s);
		Pivot(A_matrix, rhs, n, s);
		for(int i=s+1; i<n; i++)
		{
			double c=A_matrix[i][s]/A_matrix[s][s];
			rhs[i]-=rhs[s]*c;
			for(int j=s+1; j<n; j++)
			{
				A_matrix[i][j]-=A_matrix[s][j]*c;
			}
		}
	}
//    for(int i = 0; i < n; i++)
//    {
//        for(int j = 0; j < n; j ++)
//        {
//            printf("A_MATRIX[%d][%d] = %e\t ", i, j , A_matrix[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
	////////////////////////////////back replace;
	int Vsize=n;
	u[Vsize-1]=rhs[Vsize-1]/A_matrix[Vsize-1][Vsize-1];
	for(int i=Vsize-2; i>=0; i--)
	{
		double S=0.0;
		for(int j=i+1; j<Vsize; j++)
		{
			S+=A_matrix[i][j]*u[j];
		}
		u[i]=(rhs[i]-S)/A_matrix[i][i];
	}
}

void NumberSwap(double *a, double *b)
{
	double tmp = *a;
	*a = *b;
	*b = tmp;
}

void CheckIsHaveSolution(double **A_matrix, int n, int k)
{
	int i=0;
	for(i=k; i<n; i++)
	{
		if(A_matrix[i][k] != 0.0)
			break;
	}
	if(i==n)
	{
		std::cout<<"No solution or no unique solution!"<<std::endl;
		std::cout<<"Please input CONTRAL+C to stop!"<<std::endl;
		getchar();
	} 
}

void Pivot(double **A_matrix, double *Rhs, int n, int k)
{
	int t=k;
	for(int j=k; j<n; j++)
	{
		if(fabs(A_matrix[j][k])>fabs(A_matrix[t][k]))
		{
			t=j;
		}
	}
	if(t!=k)
	{
		NumberSwap(Rhs+k, Rhs+t);
		for(int j=k; j<n; j++)
		{
			NumberSwap(*(A_matrix+k)+j, *(A_matrix+t)+j);
		}
	}
}

