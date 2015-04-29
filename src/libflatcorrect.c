#include "desimage.h"
#include "mask_bits.h"
// #include "pixsupport.h"

/* based on imcorrect.c svn-30871 lines 1769-1781 */
/*  RAG: Feb 23, 2015 */
/*  Assumes a weight plane exists and has good values   */
/*  Also assumes that the "wgt" plane of a bias has a value given as  */
/*  an uncertainty in the bias image (rather than an inverse variance */
/* EHN: Update with new struct, correct propogration of errors */
int flat_c(desimage output, desimage flat) {
  int i;
  for (i=0; i<output.npixels; i++) {
      if (flat.image[i]>0.){
          output.image[i]/=flat.image[i];
	  output.weight[i] = (flat.image[i]*flat.image[i])/( (1.0/output.weight[i]) + output.image[i]/flat.weight[i] );
      }else{
          output.image[i]=0.0;
          output.weight[i]=0.0;
          output.mask[i] |= BADPIX_BPM;
      }
  }
  return(0);
}
