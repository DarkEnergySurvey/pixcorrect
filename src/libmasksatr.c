#include "mask_bits.h"
#include "desimage.h"

/* based on
 * https://dessvn.cosmology.illinois.edu/svn/desdm/devel/imsupport/trunk/src/imarithsubs.c
 * Revision: 33120 
 * lines 1486-1498 */
int column_in_section(int col,int *sec) {
  int lx = sec[0];
  int ux = sec[1];
  int npix = ux - lx;
  if(npix < 0){
    npix = ux;
    ux = lx;
    lx = npix;
  }
  return((col >= lx) && (col <= ux)); 
}

/* Based on
 * https://dessvn.cosmology.illinois.edu/svn/desdm/devel/imdetrend/trunk/src/imcorrect.c
 * REVISION: 33120 
 * lines 1673-1696
 */
int mask_saturation(desimage output, int *saturatepixels) {
  int xpos;
  int i;
  float maxsaturate;
  float image_val;

  if (output.saturateA > output.saturateB){
    maxsaturate = output.saturateA;
  } else {
    maxsaturate = output.saturateB;
  }
  
  *saturatepixels=0;
  for (i=0; i<output.npixels; i++) {

    image_val = output.image[i];

    xpos=(i%output.axes[0])+1;
    if (maxsaturate > 0.){
      if (column_in_section(xpos,output.ampsecan)) {
	/* amplifier A */
	if (image_val >= output.saturateA) {
	  output.mask[i] |= BADPIX_SATURATE;
	  (*saturatepixels)++;
	}
      } else {
	/* amplifier B */
	if (output.image[i] >= output.saturateB) {
	  output.mask[i] |= BADPIX_SATURATE;
	  (*saturatepixels)++;
	}
      }
    }
  }
  return(0);
}
