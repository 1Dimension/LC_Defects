#ifndef _CALCVECTOR_H
#define _CALCVECTOR_H

double Norm(double *fnlm,int DIMM) {
  double norm;
  norm = 0;
  for (int i=0;i<DIMM;i++)
    norm = norm + fnlm[i]*fnlm[i];
  
  norm = sqrt(norm);
  return norm;
}


double Inner(double *A, double *B, int DIMM) {
  double Inner;
  Inner = 0;
  for (int i = 0; i < DIMM; i++)
    Inner = Inner + A[i] * B[i];


  // Inner = Inner / DIMM;
  
  return Inner;
}



#endif
