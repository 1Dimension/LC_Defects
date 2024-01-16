

double interpolation(double *Q,int I,int J,int K,double r,double t,double p,int n) {
  int i,j,k;
  double dr,dt,dp;
  double x,y,z;

  dr = 1.0/n;
  dt = PI/n;
  dp = 2.0*PI/n;
  i = r/dr;
  j = t/dt;
  k = p/dp;
  if (i == I)
    i = I-1;
  if (j == J-1)
    j = J-2;
  x = r/dr - i;
  y = t/dt - j;
  z = p/dp - k;

  if (k < K-1)
    return x*y*z*Q[(i+1)*J*K + (j+1)*K + (k+1)] + x*y*(1-z)*Q[(i+1)*J*K + (j+1)*K + k] + x*(1-y)*z*Q[(i+1)*J*K + j*K + (k+1)] + x*(1-y)*(1-z)*Q[(i+1)*J*K + j*K + k] + (1-x)*y*z*Q[i*J*K + (j+1)*K + (k+1)] + (1-x)*y*(1-z)*Q[i*J*K + (j+1)*K + k] + (1-x)*(1-y)*z*Q[i*J*K + j*K + (k+1)] + (1-x)*(1-y)*(1-z)*Q[i*J*K + j*K + k];
  else
    return x*y*z*Q[(i+1)*J*K + (j+1)*K + 0] + x*y*(1-z)*Q[(i+1)*J*K + (j+1)*K + k] + x*(1-y)*z*Q[(i+1)*J*K + j*K + 0] + x*(1-y)*(1-z)*Q[(i+1)*J*K + j*K + k] + (1-x)*y*z*Q[i*J*K + (j+1)*K + 0] + (1-x)*y*(1-z)*Q[i*J*K + (j+1)*K + k] + (1-x)*(1-y)*z*Q[i*J*K + j*K + 0] + (1-x)*(1-y)*(1-z)*Q[i*J*K + j*K + k];
}

void calc_matrixH(double H[3][3],double *A,double *V) {
  int n,i,j,k,i1,i2,i3;
  double Q0[5];
  double eg[3];
  double vec[3][3];
  double *Rn00 = new double[N/2]();
  Rn00[0] = sqrt(3.0);
  for (k=1;k<N/2;++k)
    Rn00[k] = Rn00[k-1]*(-1)*sqrt(4*k+3)/sqrt(4*k-1)*(2*k+1)/(2*k);
  tranA2V(A,V);
  for (int i=0;i<5;i++) {
    Q0[i] = 0;
    for (int n=0;n<N;n+=2)
      Q0[i] += V[inx(i,n,0,0)]*Rn00[n/2];
    Q0[i] = Q0[i]/sqrt(4.0*PI);
  }
  delete [] Rn00;

  QRforEig(Q0,eg,vec);
  sort(eg,i1,i2,i3);
  
  if (eg[i2] <= 0) {
    H[0][0] = vec[0][i1];
    H[0][1] = vec[1][i1];
    H[0][2] = vec[2][i1];
    H[1][0] = vec[0][i2];
    H[1][1] = vec[1][i2];
    H[1][2] = vec[2][i2];
    H[2][0] = vec[0][i3];
    H[2][1] = vec[1][i3];
    H[2][2] = vec[2][i3];
  }
  else {
    H[0][0] = vec[0][i3];
    H[0][1] = vec[1][i3];
    H[0][2] = vec[2][i3];
    H[1][0] = vec[0][i2];
    H[1][1] = vec[1][i2];
    H[1][2] = vec[2][i2];
    H[2][0] = vec[0][i1];
    H[2][1] = vec[1][i1];
    H[2][2] = vec[2][i1];
  }
}

/* dis only equals 100 or 200 */

void rotation(double *A,double *V,int dis,int RI,int RJ,int RK) {
  int n,i,j,k,l,m,s,t;
  double H[3][3],Q[3][3],tildeQ[3][3];
  
  zer_destroy();
  zernike_init(N,L,M,dis,dis+1,dis);
  double *Q1 = new double[5*Point]();
  for (i=0;i<5;i++)
    calc_fijk(A + i*Basis,Q1 + i*Point);
  
  calc_matrixH(H,A,V);

  zernike_init(N,L,M,RI,RJ,RK);
  double *Q2 = new double[5*Point]();
  double x,y,z,tx,ty,tz,ra,th,ph;
  int ix = 0;

  for (i=0;i<I;i++) {
    for (j=0;j<J;j++) {
      for (k=0;k<K;k++) {
	x = radius[i]*sin(theta[j])*cos(phi[k]);
	y = radius[i]*sin(theta[j])*sin(phi[k]);
	z = radius[i]*cos(theta[j]);
	tx = H[0][0]*x + H[1][0]*y + H[2][0]*z;
	ty = H[0][1]*x + H[1][1]*y + H[2][1]*z;
	tz = H[0][2]*x + H[1][2]*y + H[2][2]*z;
	ra = sqrt(tx*tx + ty*ty + tz*tz);
	th = acos(tz/ra);
	if (fabs(tx) < 1e-14 && ty >= 0)
	  ph = PI/2;
	else if (fabs(tx) < 1e-14 && ty < 0)
	  ph = PI*3/2;
	else if (tx > 0 && ty >= 0)
	  ph = atan(ty/tx);
	else if (tx > 0 && ty < 0)
	  ph = 2*PI + atan(ty/tx);
	else if (tx < 0 && ty >= 0)
	  ph = PI + atan(ty/tx);
	else if (tx < 0 && ty < 0)
	  ph = PI + atan(ty/tx);
	
	/* if (th < 0 || th > PI) */
	/*   printf("ErrOR\n"); */
	/* if (ph < 0 || ph > 2*PI) */
	/*   printf("ErrOR\n"); */
	/* if (fabs(ra*sin(th)*cos(ph) - tx) > 1e-14) */
	/*   printf("ErrOR\n"); */
	/* if (fabs(ra*sin(th)*sin(ph) - ty) > 1e-14) */
	/*   printf("ErrOR\n"); */
	/* if (fabs(ra*cos(th) - tz) > 1e-14) */
	/*   printf("ErrOR\n"); */
	
	for (n=0;n<5;n++)
	  Q2[ix + n*Point] = interpolation(Q1 + n*(dis+1)*(dis+1)*dis,dis,dis+1,dis,ra,th,ph,dis);
  	Q[0][0] = Q2[ix];
  	Q[0][1] = Q2[Point + ix];
  	Q[0][2] = Q2[2*Point + ix];
  	Q[1][1] = Q2[3*Point + ix];
  	Q[1][2] = Q2[4*Point + ix];
  	Q[1][0] = Q[0][1];
  	Q[2][0] = Q[0][2];
  	Q[2][1] = Q[1][2];
  	Q[2][2] = -Q[0][0] - Q[1][1];
  	for (s=0;s<2;s++) {
  	  for (t=s;t<3;t++) {
  	    tildeQ[s][t] = 0;
  	    for (l=0;l<3;l++)
  	      for (m=0;m<3;m++)
  	    	tildeQ[s][t] += H[s][m]*Q[m][l]*H[t][l];
  	  }
  	}
  	Q2[ix] = tildeQ[0][0];
  	Q2[Point + ix] = tildeQ[0][1];
  	Q2[2*Point + ix] = tildeQ[0][2];
  	Q2[3*Point + ix] = tildeQ[1][1];
  	Q2[4*Point + ix] = tildeQ[1][2];
  	++ix;
      }
    }
  }
  
  for (i=0;i<5;i++)
    calc_fnlm(Q2 + i*Point,A + i*Basis);

  delete [] Q1;
  delete [] Q2;
}
