#include "desimage.h"

/* Dummy to demonstrate application of C library to full focal plane of images*/

int fpnumber(desimage *im) {
  int hdu, i;
  for (hdu=0; hdu<CCDNUM2; hdu++) {
    for (i=0; i<im[hdu].npixels; i++) {
      im[hdu].image[i] = hdu;
    }
  }
  return(0);
}
