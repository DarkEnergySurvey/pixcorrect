#include <math.h>

int harmonic_mean(float *ar1, float *ar2, int n) {
  int i;
  for (i=0; i<n; i++) {
    if ( (ar1[i]==0.0) || (ar2[i]==0.0) ) {
      ar1[i]=0.0;
    } else {
      ar1[i]=1.0/( (1.0/ar1[i]) + (1.0/ar2[i]) );
    }
  }
  return(0);
}

