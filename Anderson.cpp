

#include "matrix.h"
#include "Anderson.h"
#include "stdio.h"
#include "stdlib.h"


Anderson::Anderson(void)
{
	k = 0;
	m_max = 3; // 默认
	N = 10; // 默认 
	initialized = false;
	printf("\n[Anderson] a empty Anderson is created.\n");

/*
	const int n = 4;
	double x[n], xt[n], b[n];

	double **A;
	A = (double**)malloc(sizeof(double*)*n);
	A[0] = (double*)malloc(sizeof(double)*n*n);
	for(int i=1; i<n; i++){
		A[i] = A[i-1]+n;
	}

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			A[i][j] = (i+j)%n;
		}
		x[i] = i;
	}
	
	A[0][0] = 2; A[0][1] = 1;
	A[1][0] = 1; A[1][1] = 0;

	int nn = 2;
	for(int i=0; i<nn; i++){
		for(int j=0; j<nn; j++){
			printf("%g\t", A[i][j]);
		}
		printf("\n");
	}
	for(int i=0; i<nn; i++){
		b[i] = 0;
		for(int j=0; j<nn; j++){
			b[i] += A[i][j]*x[j];
		}
	}
	GaussElimination(xt, A, b, nn);
	for(int i=0; i<n; i++){
		printf("x = %g\txt = %g\tdiff = %g\n", x[i], xt[i], x[i]-xt[i]);
	}


getchar();
*/

}
Anderson::~Anderson(void)
{
	if( initialized )
		Finalize();
	printf("\n[Anderson] Anderson is distroyed.\n");
}
void Anderson::Initialize()
{
	initialized = true;
	k = 0;

	this->x_history = (double**)malloc(sizeof(double)*m_max);
	this->g_history = (double**)malloc(sizeof(double)*m_max);
	this->f_history = (double**)malloc(sizeof(double)*m_max);
	x_history[0] = (double*)malloc(sizeof(double)*m_max*N);
	g_history[0] = (double*)malloc(sizeof(double)*m_max*N);
	f_history[0] = (double*)malloc(sizeof(double)*m_max*N);
	for(int i=1; i<m_max; i++){
		x_history[i] = x_history[i-1]+N;
		g_history[i] = g_history[i-1]+N;
		f_history[i] = f_history[i-1]+N;
	}

	this->x_newadd = (double*)malloc(sizeof(double)*N);
	this->g_newadd = (double*)malloc(sizeof(double)*N);
	this->x_new = (double*)malloc(sizeof(double)*N);

	this->innarMatrix = (double**)malloc(sizeof(double*)*m_max);
	this->innarMatrix[0] = (double*)malloc(sizeof(double)*m_max*m_max);
	for(int i=1; i<m_max; i++){
		innarMatrix[i] = innarMatrix[i-1]+m_max;
	}
	this->leftVector = (double*)malloc(sizeof(double)*m_max);
	this->alpha = (double*)malloc(sizeof(double)*m_max);
}
void Anderson::Finalize()
{
	initialized = false;
	k = 0;

	free(x_history[0]);
	free(g_history[0]);
	free(f_history[0]);
	free(x_history);
	free(g_history);
	free(f_history);

	free(this->x_newadd);
	free(this->g_newadd);
	free(this->x_new);

	free(this->innarMatrix);
	free(this->leftVector);
	free(this->alpha);
}

void Anderson::AddHistory()
{
	// AndersonMixing 的step 1: k,m, update x,g,f history
	// 这里只添加history，并不做迭代
	k++; // start of k is zero
	m = k<m_max ? k : m_max;  // m = min( k, m_max)
	int id_newadd = (k-1)%m_max;
	for(int ir=0; ir<N; ir++){
		x_history[id_newadd][ir] = x_newadd[ir];
		g_history[id_newadd][ir] = g_newadd[ir];
		f_history[id_newadd][ir] = g_newadd[ir] - x_newadd[ir];
	}
}
void Anderson::AndersonMixing(double updateRatio)
{
	// 使用时，新的数据copy到g_newadd中
	// 使用Solve()之后从x_new中copy迭代结果

	/*// step 1: k,m, update x,g,f history
	k++; // start of k is zero
	m = k<m_max ? k : m_max;  // m = min( k, m_max)
	int id_newadd = (k-1)%m_max;
	for(int ir=0; ir<N; ir++){
		x_history[id_newadd][ir] = x_newadd[ir];
		g_history[id_newadd][ir] = g_newadd[ir];
		f_history[id_newadd][ir] = g_newadd[ir] - x_newadd[ir];
	}//*/

	if( k < 2 )
		goto update; // 开始时没有方程要解

	// step 2: update innarMatrix, (m-1)-by-(m-1)
	for(int i=0; i<m-1; i++){
		double temp;
		for(int j=0; j<=i; j++){
			temp = 0;
			for(int ir=0; ir<N; ir++)
				temp += (f_history[i][ir] - f_history[m-1][ir]) *(f_history[j][ir] - f_history[m-1][ir]); // innar product

			this->innarMatrix[i][j] = temp;
			this->innarMatrix[j][i] = temp; // symmetry
		}

		temp = 0;
		for(int ir=0; ir<N; ir++)
			temp += (f_history[i][ir] - f_history[m-1][ir])* f_history[m-1][ir];
		this->leftVector[i] = -temp;
	}

       
	/**
	for(int i=0; i<m-1; i++){
		double Ax = 0;
		for(int j=0; j<m-1; j++){
			printf("%g\t", innarMatrix[i][j]);
		}
		printf("| %g", leftVector[i]);
		printf("\n");
      
	}
	// */

	// step 3: solve alpha
	GaussElimination(alpha, innarMatrix, leftVector, m-1);
	//GaussElimination(double *u, double **A_matrix, double *rhs, int n)

update:
	// step 4: get x_new
	alpha[m-1] = 1;
	for(int i=0; i<m-1; i++)
		alpha[m-1] -= alpha[i]; // sum(alpha) = 1

	
	/**
	for(int i=0; i<m; i++){
		printf("alpha[%d] = %g\n", i, alpha[i]);
		}//*/

	printf("[AndersonMixing] k = %d, m = %d, updateRatio = %g: ", k, m, updateRatio);
	for(int ir=0; ir<N; ir++){
		x_new[ir] = 0;
		for(int i=0; i<m; i++)
			x_new[ir] += alpha[i] * ( x_history[i][ir] + updateRatio*f_history[i][ir] );
	}
}
