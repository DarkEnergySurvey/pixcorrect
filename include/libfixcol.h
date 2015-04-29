#include "desimage.h"
#include "mask_bits.h"

int colAve(int icol, short int bitson,
	   desimage bpm, desimage output,
	   double *pixavg, double *pixrms, double *pavgerr);

int fixCol(desimage bpm, desimage output);

void johnQsort(int gindex[],double value[], int last);
