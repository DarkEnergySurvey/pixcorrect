#include "libfixcol.h"
#include <stdio.h>
#include <math.h>

/* based on 
 * https://dessvn.cosmology.illinois.edu/svn/desdm/devel/imsupport/trunk/src/mask_utils.c
 * revision 33120 */

int colAve(int icol,short int bitson,desimage bpm,desimage output,double *pixavg,double *pixrms,double *pavgerr)
{
  /*
  Routine to compute the sky level for a column.

  Routine returns 0 for success, non-zero is failure
  icol = compute average for this column (first column=0)
  bitson = require bpm mask bit to be set for this computation
           no other bits may be set
  bitson = 0 means no bits required to be set, no bits may  be set
  bpm = bad pixel mask image
  output = output image
  pixavg = column average
  pixrms = column rms
  pavgerr = estimated error on column average
  The routine computes the median and then clips the distribution
  and recomputes the median.  
  It seems to be fairly robust in practice, but could fail if there are
  no pixels that are at the sky level because of bleed trails or very
  large objects
  */

  //Minimum number of samples required to compute the average.
  const int MINSAMP=100;
  //Arrays are sized for DES 4096 pixel columns
  int gindex[4096];
  double pixel[4096];
  int j, n;
  int ie, in, io;
  int il, iu, lower, upper;
  int nl, nu;
  double imval, ave, sigma1, sigma2, err;
  int nsample;
  short int bitsoff;
  bitsoff = ~bitson;
  //Code assumes bpm and output have same image dimensions
  io = 0;
  //Loop over rows and find pixels with the requested bpm bits
  for (j=0;j<output.axes[1];++j)
    {
      in = output.axes[0]*j + icol;
       imval = output.image[in];
      //first check image for NAN or INF values
     if (0.0*imval!=0.0) continue;
     if ((bitson!=0) && ((bpm.mask[in]&bitson)==0)) continue;
     if (bitsoff&bpm.mask[in]) continue;
      pixel[io] = imval;
      gindex[io] = io;
      ++io;
    }
  nsample = io;
  //Check for required minimum sample size
  if (nsample<MINSAMP)  return(1);
  //Sort pixels by value
  johnQsort(gindex,pixel,nsample-1);
  //Find median
  n = 0.5*nsample;
  ie = gindex[n];
  ave = pixel[ie];
  // Find lower 1 sigma point
  n = 0.1587*nsample;
  ie = gindex[n];
  sigma1 = ave - pixel[ie];
  //Refine estimate by clipping at 4 sigma
  lower = ave - 4.0*sigma1;
  for (n=0;n<nsample;++n)
    {
      ie = gindex[n];
      if (pixel[ie]>=lower) break;
    }
  nl = n;
  upper = ave + 4.0*sigma1;
  for (n=nsample-1;n>0;--n)
    {
      ie = gindex[n];
      if (pixel[ie]<=upper) break;
    }
  nu = n;
  //Get median of clipped distribution 
  n = 0.5*(nl+nu);
  ie = gindex[n];
  ave = pixel[ie];
  //Compute lower 1 sigma point
  n = 0.1587*(nu-nl) + nl;
  il = gindex[n];
  sigma1 = ave - pixel[il];
  //and upper 1 sigma point
  n = 0.8413*(nu-nl) + nl;
  iu = gindex[n];
  sigma2 = pixel[iu] - ave;
  //Estimated error on median
  err = sigma1/sqrt(nu-nl);
  //Median
  *pixavg = ave;
  //Rms (actually +/-34% width of distribution)
  *pixrms = 0.5*(sigma1+sigma2);
  //Error on median (approximate)
  *pavgerr = err;
  return(0);
}


int fixCol(desimage bpm,desimage output)
{
  /*
  Routine to fix columns flagged in bpm
  bpm = bad pixel array (input)
  output = output image to be modified
  */
  const int NUMCOL=10;
  const double SIGTOL=10.0;
  int i, j, jo, js, jl;
  int icol, ilow, ihigh;
  /* short int bitson, bitsoff; */
  /* int bad, nhot; */
  double pixavg, pixrms, pavgerr;
  int ncol, sign;
  double wgt;
  double corr, corr_err;
  double col_ave, col_rms, col_err;
  double norm_ave, norm_rms, norm_err;
  int nsat;
  printf("Fixcol: Looking for correctable columns.\n");
  //Loop over columns looking for "fixable bad columns"
  for (i=0;i<output.axes[0];++i)
    {
      js = i;
      jl = output.axes[0]*(output.axes[1]-1) + js;
      icol = i;
      //Check for columns with correctable pixels
      if (!(bpm.mask[js]&BPMDEF_CORR) && !(bpm.mask[jl]&BPMDEF_CORR)) continue;
      printf(" Fixcol: Column=%i is correctable.\n",i);
      if (colAve(icol,BPMDEF_CORR,bpm,output,&pixavg,&pixrms,&pavgerr)) 
	{	
	  
	  continue;
	}
      //Count saturated pixels
      nsat = 0;
      for (j=0;j<output.axes[1];++j)
	{
	  jo = j*output.axes[0] + icol;
	  if (output.mask[jo]&BADPIX_SATURATE) ++nsat;
	}
      //If saturated pixels are found, skip the correction.
      //If there are saturated pixels there could be a long bleed trail
      //that would throw off the column average. It may be OK to make 
      //the correction in those cases but it seems safer to just skip the
      //small correction in order to avoid the possibility of a large
      //non-sense correction
      if (nsat>0) 
	{
	  printf(" Fixcol: Column=%i has a saturated pixel.  Skip correction.\n",i);
	  continue;
	}
      //column average sky level
      col_ave = pixavg;
      //estimated error on sky level
      col_err = pavgerr;
      //rms of the pixels in the column
      col_rms = pixrms;
      //find NUMCOL normal columns (half above and half below the column to 
      //be corrected
      sign = 1;
      ilow = icol;
      ihigh = icol;
      ncol = 0;
      norm_ave = 0.0;
      norm_err = 0.0;
      norm_rms = 0.0;
      while (ncol<NUMCOL && (ilow>=0 || ihigh<output.axes[0]))
	{
	  icol = -1;
	  if (ilow<0) sign = +1;
	  if (ihigh>=output.axes[0]) sign = -1; 
	  //Find next lower "normal" column
	  if (sign<0)
	    {
	      --ilow;
	      while (ilow>=0)
		{
		  js = ilow;
		  jl = output.axes[0]*(output.axes[1]-1) + js;
		  if ((bpm.mask[js]&BPMDEF_CORR) || (bpm.mask[jl]&BPMDEF_CORR)) --ilow;
		  else
		    {
		      icol = ilow;
		      break;
		    }
		}
	    }
	  //Find next higher "normal" column
	  if (sign>0)
	    {
	      ++ihigh;
	      while (ihigh<output.axes[0])
		{
		  js = ihigh;
		  jl = output.axes[0]*(output.axes[1]-1) + js;
		  if ((bpm.mask[js]&BPMDEF_CORR) || (bpm.mask[jl]&BPMDEF_CORR)) ++ihigh;
		  else
		    {
		      icol = ihigh;
		      break;
		    }
		}
	    }
	  if (icol<0) continue;
	  if (!colAve(icol,0,bpm,output,&pixavg,&pixrms,&pavgerr)) 
	    {
	      wgt = 1.0/(pavgerr*pavgerr);
	      norm_ave += (wgt*pixavg);
	      norm_err += wgt;
	      norm_rms += wgt*pixrms;
	      ++ncol;
	      sign *= -1;
	    }
	}
      //Just skip correction if we found fewer than NUMCOL reference columns.
      //This is pathological and should never happen
      if (ncol!=NUMCOL) 
	{
	  printf(" Fixcol: Column=%i could not be corrected because %i normal columns were not found\n",icol,NUMCOL);
	  continue;
	}
      //Compute weighted average of normal columns
      norm_ave /= norm_err;
      norm_rms /= norm_err;
      norm_err = 1.0/sqrt(norm_err);
      //Correction is the difference between the column and the nominal level
      corr = col_ave - norm_ave;
      //Errors added in quadrature
      corr_err = sqrt(col_err*col_err+norm_err*norm_err);
      //Skip correction if column rms is different than the nominal rms.
      //This is a sign of trouble in computing the column average
      if (fabs(col_rms-norm_rms)>SIGTOL*corr_err) 
	{
	  printf(" Fixcol: Column=%i skipped because col_rms=%.1f disagrees with nominal=%.1f "
		 "outside expected tolerance=%.1f \n",
		 i,col_rms,norm_rms,SIGTOL*corr_err);
	  continue;
	}
      //Skip correction is correction is less than 5 sigma
      if (fabs(corr/corr_err)<5.0){
	  	  printf(" Fixcol: Column=%i skipped because correction=%.1f +/- %.1f "
			 "is not significant\n",i,corr,corr_err);
                  printf(" Fixcol: Unsetting BADPIX_BPM for column=%i if no other flag present in BPM\n",i);
                  for (j=0;j<output.axes[1];++j){
	             jo = j*output.axes[0] + i;
                     if (bpm.mask[jo]&BPMDEF_CORR){
                        output.mask[jo]&=~BADPIX_BPM;
                     }
                  }
		  continue;
      }else{
	  printf(" Fixcol: Column=%i Correction=%.1f +/- %.1f \n",i,corr,corr_err);
	  printf(" Fixcol: Column sky=%.1f Nominal sky=%.1f Column rms=%.1f "
		 "Nominal rms=%.1f \n",
		 col_ave,norm_ave,col_rms,norm_rms);
	  for (j=0;j<output.axes[1];++j){
	     jo = j*output.axes[0] + i;
	     if (bpm.mask[jo]&BPMDEF_CORR){
                output.image[jo] -= corr;
                output.mask[jo] &= ~BADPIX_BPM;
                output.mask[jo] |= BADPIX_FIXED;
             }
         }
      }
    }
  return(0);
}

//This is John Marriner's version of quick sort 
//-- probably originally from Numerical Recipes
#define INSERTION_SORT_BOUND 16 /* boundary point to use insertion sort */
void johnQsort(int gindex[],double value[], int last)
{
  //Quick sort routine
  //Sort according to value.  Last element in array is value[last]
  //On entry gindex[n] = n.  On exit value[gindex[n]] will be ordered for n=0,1,2,...

  int stack_pointer = 0;
  int first_stack[32];
  int last_stack[32];
  int ifirst, ilast, imed, idown, iup;
  int first=0;
  for (;;)
  {
    if (last - first <= INSERTION_SORT_BOUND)
    {
      /* for small sort, use insertion sort */
      int indx;
      int prev_val = gindex[first];
      int cur_val;

      for (indx = first + 1; indx <= last; ++indx)
      {
        cur_val = gindex[indx];
        if (value[prev_val]>value[cur_val])
        {
          /* out of order: array[indx-1] > array[indx] */
          int indx2;
          gindex[indx] = prev_val; /* move up the larger item first */

          /* find the insertion point for the smaller item */
          for (indx2 = indx - 1; indx2 > first; )
          {
            int temp_val = gindex[indx2 - 1];
            if (value[temp_val]>value[cur_val])
            {
              gindex[indx2--] = temp_val;
              /* still out of order, move up 1 slot to make room */
            }
            else
              break;
          }
          gindex[indx2] = cur_val; /* insert the smaller item right here */
        }
        else
        {
          /* in order, advance to next element */
          prev_val = cur_val;
        }
      }
    }
    else
    {
      int pivot;
 
      /* try quick sort */
      {
        int temp;
        int med = (first + last) >> 1;
        /* Choose pivot from first, last, and median position. */
        /* Sort the three elements. */
        temp = gindex[first];
	ilast = gindex[last];
        if (value[temp]>value[ilast])
	  {
	    gindex[first] = gindex[last]; gindex[last] = temp;
	  }
        temp = gindex[med];
	ifirst = gindex[first];
        if (value[ifirst]>value[temp])
        {
          gindex[med] = gindex[first]; gindex[first] = temp;
        }
        temp = gindex[last];
	    imed = gindex[med];
        if (value[imed]>value[temp])
        {
          gindex[last] = gindex[med]; gindex[med] = temp;
        }
        pivot = gindex[med];
      }
      {
        int up;
        {
	  int down;
          /* First and last element will be loop stopper. */
	  /* Split array into two partitions. */
	  down = first;
	  up = last;
	  for (;;)
	  {
	    do
	    {
	      ++down;
	      idown = gindex[down];
	    } while (value[pivot]>value[idown]); 
	    do
	    {
	      --up;
	      iup = gindex[up];
	    } while (value[iup]>value[pivot]);
 
	    if (up > down)
	    {
	      int temp;
	      /* interchange L[down] and L[up] */
	      temp = gindex[down]; gindex[down]= gindex[up]; gindex[up] = temp;
	    }
	    else
	      break;
	  }
	}
	{
	  int len1; /* length of first segment */
	  int len2; /* length of second segment */
	  len1 = up - first + 1;
	  len2 = last - up;
	  /* stack the partition that is larger */
	  if (len1 >= len2)
	  {
	    first_stack[stack_pointer] = first;
	    last_stack[stack_pointer++] = up;
 
	    first = up + 1;
	    /*  tail recursion elimination of
	     *  johnQsort(gindex,fun_ptr,up + 1,last)
	     */
	  }
	  else
	  {
	    first_stack[stack_pointer] = up + 1;
	    last_stack[stack_pointer++] = last;

	    last = up;
	    /* tail recursion elimination of
	     * johnQsort(gindex,fun_ptr,first,up)
	     */
	  }
	}
        continue;
      }
      /* end of quick sort */
    }
    if (stack_pointer > 0)
    {
      /* Sort segment from stack. */
      first = first_stack[--stack_pointer];
      last = last_stack[stack_pointer];
    }
    else
      break;
  } /* end for */
  return;
}
