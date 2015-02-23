#include "despyfits/desimage.h"
#include "despyfits/mask_bits.h"

/* based on imcorrect.c svn-30871 lines 1466-1475 */
int bias_c(desimage output, desimage bias) {
  int i;
  for (i=0; i<output.npixels; i++) {
    output.image[i]-=bias.image[i];
/*  RAG: Feb 23, 2015 */
/*  Assumes a weight plane exists and has good values   */
/*  Also assumes that the "wgt" plane of a bias has a value given as  */
/*  an uncertainty in the bias image (rather than an inverse variance */
    if (bias.varim[i]>0.0){
       output.varim[i]+=1.0/(double)bias.varim[i];
    }
  }
  return(0);
}
