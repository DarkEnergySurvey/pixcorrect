int foo(double *ar1, double *ar2, double coeff, int n) {
  int i;
  for (i=0; i<n; i++) {
    ar1[i] = coeff*ar1[i] + ar2[i];
  }
  return(0);
}

