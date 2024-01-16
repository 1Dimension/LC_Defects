#pragma once

// Anderson mixing

class Anderson
{
public:
	int m, m_max; // m = min(m_max, k)
	int k; // 
	int N;
	bool initialized;
	
	double *x_newadd;
	double *g_newadd;
	double *x_new; // as output

private:
	double **x_history; 
	double **g_history;
	double **f_history; // f = g(x_i) - x_i

	double **innarMatrix; // innar product
	double *leftVector;
	double *alpha; // update factor

public:
	Anderson(void);
	virtual ~Anderson(void);

	void Initialize();
	void Finalize();
	void AddHistory();
	void AndersonMixing(double updateRatio);
};

