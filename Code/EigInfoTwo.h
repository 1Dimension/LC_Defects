
bool qr(int n, double *b, double *c, double *q, double eps, int l){
  int i,j,k,m,it,u,v;
  double d,f,h,g,p,r,e,s;

  c[n-1]=0.0; d=0.0; f=0.0;

  for (j=0; j<=n-1; j++){
    it=0;
    h=eps*(fabs(b[j])+fabs(c[j]));
    if (h>d) d=h;
    m=j;
    while((m<=n-1)&&(fabs(c[m])>d)) m++;
    if(m!=j){
      do{
	if (it==l) return 1;
	it=it+1;
	g=b[j];
	p=(b[j+1]-g)/(2.0*c[j]);
	r=sqrt(p*p+1.0);
	if(p>=0.0)
	  b[j]=c[j]/(p+r);
	else
	  b[j]=c[j]/(p-r);
	h=g-b[j];
	for(i=j+1; i<=n-1; i++) b[i]=b[i]-h;
	f=f+h; p=b[m]; e=1.0; s=0.0;
	for (i=m-1; i>=j; i--){
	  g=e*c[i]; h=e*p;
	  if(fabs(p)>=fabs(c[i])){
	    e=c[i]/p;
	    r=sqrt(e*e+1.0);
	    c[i+1]=s*p*r; s=e/r; e=1.0/r;
	  }
	  else{
	    e=p/c[i]; r=sqrt(e*e+1.0);
	    c[i+1]=s*c[i]*r;
	    s=1.0/r; e=e/r;
	  }
	  p=e*b[i]-s*g;
	  b[i+1]=h+s*(e*g+s*b[i]);
	  for(k=0; k<=n-1; k++){
	    u=k*n+i+1; v=u-1;
	    h=q[u]; q[u]=s*q[v]+e*h;
	    q[v]=e*q[v]-s*h;
	  }
	}
	c[j]=s*p;
	b[j]=e*p;
      }while (fabs(c[j])>d);
    }
    b[j]=b[j]+f;
  }
  for(i=0; i<=n-1; i++){
    k=i; p=b[i];
    if(i+1<=n-1){
      j=i+1;
      while((j<=n-1)&&(b[j]<=p)){
	k=j;
	p=b[j];
	j=j+1;
      }
    }
    if(k!=i){
      b[k]=b[i];
      b[i]=p;
      for(j=0; j<=n-1; j++){
	u=j*n+i;
	v=j*n+k;
	p=q[u]; q[u]=q[v]; q[v]=p;
      }
    }
  }
  return 0;
}

void householder(double *a, int n, double *q, double *b, double *c) {
  int i,j,k,u;
  double h,f,g,h2;

  for (i=0; i<=n-1; i++){
    for (j=0; j<=n-1; j++){
      u=i*n+j;
      q[u]=a[u];
    }
  }
  for (i=n-1; i>=1; i--){
    h=0.0;
    if(i>1){
      for (k=0; k<=i-1; k++){
	u=i*n+k;
	h=h+q[u]*q[u];
      }
    }
    if (h+1.0==1.0){
      c[i]=0.0;
      if (i==1) c[i]=q[i*n+i-1];
      b[i]=0.0;
    }
    else{
      c[i]=sqrt(h);
      u=i*n+i-1;
      if (q[u]>0.0) c[i]=-c[i];
      h=h-q[u]*c[i];
      q[u]=q[u]-c[i];
      f=0.0;
      for(j=0; j<=i-1; j++){
	q[j*n+i]=q[i*n+j]/h;
	g=0.0;
	for(k=0; k<=j; k++) g=g+q[j*n+k]*q[i*n+k];
	if(j+1<=i-1){
	  for (k=j+1; k<=i-1; k++) g=g+q[k*n+j]*q[i*n+k];
	}
	c[j]=g/h;
	f=f+g*q[j*n+i];
      }
      h2=f/(h+h);
      for(j=0; j<=i-1; j++){
	f=q[i*n+j];
	g=c[j]-h2*f;
	c[j]=g;
	for(k=0; k<=j; k++){
	  u=j*n+k;
	  q[u]=q[u]-f*c[k]-g*q[i*n+k];
	}
      }
      b[i]=h;
    }
  }
  for (i=0; i<=n-2; i++) c[i]=c[i+1];
  c[n-1]=0.0;
  b[0]=0.0;
  for (i=0; i<=n-1; i++){
    if((b[i]!=0.0)&&(i-1>=0)){
      for (j=0; j<=i-1; j++){
	g=0.0;
	for (k=0; k<=i-1; k++) g=g+q[i*n+k]*q[k*n+j];
	for (k=0; k<=i-1; k++){
	  u=k*n+j;
	  q[u]=q[u]-g*q[k*n+i];
	}
      }
    }
    u=i*n+i;
    b[i]=q[u]; q[u]=1.0;
    if (i-1>=0){
      for (j=0; j<=i-1; j++){
	q[i*n+j]=0.0;
	q[j*n+i]=0.0;
      }
    }
  }
}



/***************************************
 * Q5: five numbers;
 * EigValue[3]: eigenvalues;
 * EigVector[3][3]: eigenvectors. 
 ***************************************/


void QRforEig( double *Q5, double *EigValue, double (*EigVector)[3] )
{
  double a[9],b[3],c[3],q[9];
  double eps = 8.0e-16;
  int l = 1000000;
  int n =3;
  bool QRflag=0;
	
  a[0] = Q5[0];
  a[1] = Q5[1];
  a[2] = Q5[2];
  a[3] = Q5[1];
  a[4] = Q5[3];
  a[5] = Q5[4];
  a[6] = Q5[2];
  a[7] = Q5[4];
  a[8] = - Q5[0] - Q5[3];

  householder(a, n, q, b, c);
  QRflag = qr(n, b, c, q, eps, l);

  if (QRflag)
    {
      printf("Error in QR method for eigenvalues: step is large than %d!\n",l);
      exit(1);

    }

  for (int i=0;i<3;i++)
    {
      EigValue[i] = b[i];
    }


  for (int i=0;i<9;i++)
    {
      EigVector[i/3][i%3] = q[i];
    }
}
