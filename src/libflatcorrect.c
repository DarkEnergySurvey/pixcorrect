#include "despyfits/desimage.h"
#include "despyfits/mask_bits.h"
#include "despyfits/pixsupport.h"

/* based on imcorrect.c svn-30871 lines 1769-1781 */
/*  RAG: Feb 23, 2015 */
/*  Assumes a weight plane exists and has good values   */
/*  Also assumes that the "wgt" plane of a bias has a value given as  */
/*  an uncertainty in the bias image (rather than an inverse variance */

int flat_c(desimage output, desimage flat) {
  int i;
  for (i=0; i<output.npixels; i++) {
      if (flat.image[i]>0.){
          output.image[i]/=flat.image[i];
          output.varim[i]/=Squ((double)flat.image[i]);
          if (flat.varim[i]>0.0){
              output.varim[i]+=1.0/((double)flat.varim[i]);
          }
      }else{
          output.image[i]=0.0;
          output.varim[i]=0.0;
          output.mask[i] |= BADPIX_BPM;
      }
  }
  return(0);
}
