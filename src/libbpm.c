#include "mask_bits.h"
#include "desimage.h"

/* based on imcorrect.c svn-30871 lines 1466-1475 */
int obpm(desimage output, desimage bpm) {
  int i;
  for (i=0; i<output.npixels; i++) {
    if (bpm.mask[i]) {
      output.mask[i] = BADPIX_BPM;
      if (bpm.mask[i] & BPMDEF_EDGE) {
	output.mask[i] |= BADPIX_EDGE;
      }
    } else {
      output.mask[i] = 0;
    }
  }
  
  return(0);
}

/* based on imcorrect.c svn-30871 lines 1480-1485 */
int bpm(desimage output, desimage bpm) {
  int i;
  for (i=0; i<output.npixels; i++) {
    if (bpm.mask[i]) {
      output.mask[i] |= BADPIX_BPM;
      if (bpm.mask[i] & BPMDEF_EDGE) {
	output.mask[i] |= BADPIX_EDGE;
      }
    }
  }
  return(0);
}

