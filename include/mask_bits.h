/* define BADPIX bit mappings (for MASK HDU) */
#define BADPIX_BPM 1          /* set in bpm (hot/dead pixel/column)        */
#define BADPIX_SATURATE 2     /* saturated pixel                           */
#define BADPIX_INTERP 4	      /* interpolated pixel                        */
#define BADPIX_THRESHOLD 0.10 /* pixels less than this fraction of sky     */
                              /* are filtered -- helps remove failed reads */
#define BADPIX_LOW     8      /* too little signal- i.e. poor read         */
#define BADPIX_CRAY   16      /* cosmic ray pixel                          */
#define BADPIX_STAR   32      /* bright star pixel                         */
#define BADPIX_TRAIL  64      /* bleed trail pixel                         */
#define BADPIX_EDGEBLEED 128  /* edge bleed pixel                          */
#define BADPIX_SSXTALK 256    /* pixel potentially effected by xtalk
                                               from super-saturated source */
#define BADPIX_EDGE   512     /* pixel flagged to exclude CCD glowing 
                                                                     edges */
#define BADPIX_STREAK 1024    /* pixel associated with satellite 
                                                   (meteor/ufo) streak     */
#define BADPIX_FIX    2048    /* a bad pixel that was fixed                */

/* define BPMDEF bit mappings (for BPM definition) */
#define BPMDEF_FLAT_MIN 1     /* Pixels that are dull in the flats.        */
#define BPMDEF_FLAT_MAX 2     /* Pixels that are hot in the flats.         */
#define BPMDEF_FLAT_MASK 4    /* Pixels that are in the BPM for the flats. */
#define BPMDEF_BIAS_HOT 8     /* Pixels that are hot in the biases.        */
#define BPMDEF_BIAS_WARM 16   /* Pixels that are warm in the biases.       */
#define BPMDEF_BIAS_MASK 32   /* Pixels that are in the BPM for the biases */
#define BPMDEF_BIAS_COL 64    /* Pixels that are downstream of a hot pixel 
                                                              in the bias. */
#define BPMDEF_EDGE 128       /* Pixels on the glowing edges of the CCD.   */
#define BPMDEF_CORR 256       /* Correctable pixels 
                                       (usually downstream of hot pixels). */
#define BPMDEF_TAPE_BUMP 512  /* Pixels that reside in tape bumps.         */
#define BPMDEF_FUNKY_COL 1024 /* Columns with charge redistribution in sky 
                                                                exposures. */
#define BPMDEF_WACKY_PIX 2048 /* Outliers in stacked sky exposures.        */
#define BPMDEF_GENERIC   4096 /* Generic (undifferentiated) bad pixel flag */
