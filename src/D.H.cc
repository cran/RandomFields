/* 
 Authors
 Martin Schlather, martin.schlather@cu.lu

 library for unconditional simulation of stationary and isotropic random fields

 Copyright (C) 2002 - 2004 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <math.h>  
#include <stdio.h>  
#include <assert.h>
#include "RFsimu.h"

void boxcounting(double *z,  // data
		 int *lx,    // basic length of data 
		 int *repet, // "independent" repetitions of data 
		 double *factor, // rescalin of data values by factor
		 int *eps,   // box lengths
		 int *leps,  // length of eps
		 double *sum // leps * repet
    ) {
  // in sum, the boxcounting result for the first lx data values is given 
  // according to the box lengths given by eps; then the boxes are counted
  // for the next lx data values, etc. Hence, sum must have dimension
  int i,  k, lastbox, j, e, r, s, truelx, total;
  double zz, min, max, f; 
  // boxcounting expects that the first and the last values are repeated
  // in the first dimension, see D.H.R; accordingly the total length in
  // the first direction is increased by 2:
  truelx = 2 + *lx;
  total = *repet * truelx; 
  s = 0; // 0 to leps * repet -1
  for (r=0; r<total; r+=truelx) {
    for (k=0; k<*leps; k++, s++) {
      sum[s] = 0;
      e = eps[k]; // box size
      assert(e <= *lx); 
      lastbox = r + (e * (*lx / e)) + 1 - e;
      for (i = 1 + r; i <= lastbox; ) {
	// # boxen um  Daten von 0.5 * (i + i-1) bis 0.5 * (i + i+1)
	// zu ueberdecken
	min = max = 0.5 * (z[i] + z[i-1]);
	for (j=0; j<e; j++, i++) 
	  if (z[i] < min) min = z[i]; else if (z[i] > max) max = z[i];
	assert(i <= r + *lx + 1);
	if ((zz = 0.5 * (z[i-1] + z[i])) < min) min = zz; 
	else if (zz > max) max = zz;
	f = *factor / (double) e;
	sum[s] += floor(max * f) - floor(min * f) + 1.0;
      }
    }
  }
  return;
}

void periodogram(double *dat, // data
		 int *len,    // basic length of data
		 int *repet,  // # of "independent" repetitions
		 int *fftm,  // fftm[0] >=1 not >=1 (R indexing assumed)
		 //             interval of the vector given by fft of part
		 int *part,  // not the basic length but only a part is
                 //             transformed
		 int *shift, // part is then shifted by shift until end of
		 //             vector is reached; fft values are averaged
                 //             over all these parts
		 double *lambda // (fftm[1] - fftm[0] + 1) * repet
    ){ 
  int seg_dat, Xerror, k, j, total_seg, segment_r, start_k, end_k,
      segm_l, lenMpart, delta_l, r, end_l;
  double factor, *compl_number, n_inv, *taper, taper_fact, cos_factor;
  FFT_storage FFT;
  
  assert(*len >= *part);
  assert((fftm[0]>=1) && (fftm[1] <=*part) && (fftm[1] >= fftm[0]));
  assert(*shift>=1);
  
  factor = log(2 * PI * *len);
  lenMpart = *len - *part;
  start_k = 2 * (fftm[0] - 1); // of complex numbers
  end_k = 2 * (fftm[1] - 1);   // of complex numbers
  delta_l = fftm[1] - fftm[0] + 1;
  end_l = delta_l * *repet;
  n_inv = 1.0 / (double) ((int) (1.0 + (*len - *part) / *shift));
  FFT_NULL(&FFT);
  compl_number = NULL;
  taper = NULL;
  if ((compl_number = (double*) malloc(sizeof(double) * 2 * *part))==NULL){
    goto ErrorHandling;
  }
  if ((taper = (double*) malloc(sizeof(double) * *part))==NULL){
    goto ErrorHandling;
  }
  
  for (j=0; j<end_l; j++) lambda[j]=0.0;
  taper_fact = sqrt( 2.0 / (3.0 * ((double) *part + 1.0)) );
  cos_factor = 2.0 * PI / ((double) *part + 1.0);
  for (j=0; j<*part; j++) {
    // taper is used to reduce bias for high frequencies by smoothing start 
    // and end of series towards zero
    // However a high bias for the smallest one or two frequencies is introduced 
    taper[j] = taper_fact * ( 1.0 - cos (cos_factor * (double) (j + 1)));
  }

  for (segment_r=segm_l=r=0; r<*repet; r++, segment_r+=*len, segm_l+=delta_l) {
    // repet
    for (seg_dat=0; seg_dat<=lenMpart; seg_dat+=*shift) {
      // shifts
      total_seg = seg_dat + segment_r;
      for (k=0, j=0; j<*part ; j++) {
	compl_number[k++] = dat[total_seg + j] * taper[j];
	compl_number[k++] = 0.0;
      }
      if ((Xerror=fastfourier(compl_number,
			      part,   // length
			      1,      // dim 
			      (total_seg==0), // start? 
			      false,  // inverse?  
			      &FFT))!=0)     // aux 
	goto ErrorHandling;
      for (j=segm_l, k=start_k; k<end_k; j++, k+=2) {
	lambda[j] += log(compl_number[k] * compl_number[k] + 
			 compl_number[k+1] * compl_number[k+1]) - factor;
      }
    }
  }
  for (j=0; j<end_l; j++) lambda[j] *= n_inv; // averaging over shifts
  free(compl_number);
  free(taper);
  FFT_destruct(&FFT);
  return;
 ErrorHandling:
  if (compl_number!=NULL) free(compl_number);
  if (taper!=NULL) free(taper);
  for (j=0; j<end_l; j++) lambda[j] = RF_NAN;
  FFT_destruct(&FFT);
  return;
}

// both detrended fluctuation and aggreated variation  
void detrendedfluc(double *dat, // data
		   int *lx,     // basic length of data
		   int *repet,  // # of "independent" repetitions
		   int* blocks, // block (box) lengths
		   int *ldfa,   // # of blocks
		   double *dfavar,   // repet * ldfa
		   double *varmethvar// repet * ldfa
  ){
  int i, j, m, nbox, ex, endfor, r, total, last, e, k;
  double var, a, b, residual, Yt, realm, sumt, meanY, t, 
    VarMeth_old, VarMeth_mean, VarMeth_var, delta, realnbox;

  assert(*lx > 1);
  total = *lx * *repet;
  for (e = r = 0; r < total; r += *lx, e+=*ldfa) {
    last = r + *lx;
    for (j=r+1; j<last; j++) { dat[j] += dat[j-1]; }
    for (i=0; i<*ldfa; i++){
      var = 0.0;
      m = blocks[i];
      realm = (double) m;
      nbox = *lx / m;      
      realnbox = (double) nbox;
      ex = r + m * nbox;
      sumt = 0.5 * realm * (realm + 1.0); // sum_i^m i
      //sumtsq = realm * (realm + 1.0) * (2.0 * realm + 1.0) / 6.0;// sum_i^m i^2

      // aggregated variation
      if (nbox > 1) {
	VarMeth_mean = dat[ex - 1] / (realnbox);
	VarMeth_old = 0;
	VarMeth_var = 0;
	for (j=r+m-1; j<ex; j+=m) {
	  delta = (dat[j] - VarMeth_old) - VarMeth_mean;
	  VarMeth_var += delta * delta;
	  VarMeth_old = dat[j];
	}
        varmethvar[e + i] = log(VarMeth_var / ((realnbox - 1.0)));
        // varmethvar[e + i] = log(VarMeth_var / ((realnbox)));
      } else {
	varmethvar[e + i] = RF_NAN;
      }

      // detrended fluctuation
      for (j=r; j<ex; j+=m) {
	endfor = j + m;
	// mean      
	for (Yt = meanY = 0.0, k=j, t=1.0; k<endfor; k++, t+=1.0) {
	  meanY += dat[k];
	  Yt += dat[k] * t;
	}
	meanY /= realm;
	b = (Yt - meanY * sumt) * 12 / (realm * (realm + 1.0) * (realm - 1.0));
	a = meanY - b * sumt / realm;
	for (k=j, t=1.0; k<endfor; k++, t+=1.0) {
	  residual = dat[k] - (a + b * t);
	  var += residual * residual;
	}
      }
      dfavar[e + i] = log(var / ((realm - 1.0) * (double) nbox));
    }
  }
}


// range method
#define FRACT_MAXDIM 10
void minmax(double *dat,  // data
	    int *lx,      // basic length
	    int *repet,   // repetitions
	    int *boxes,   // box lengths
	    int *lb,      // # of box lengths
	    double *count // repet * lb
  ) {
  int b, tb, cb, r, totalblocks, epsilon, idx, start, end;
  double min, max;
  for (start=cb=r=0; r<*repet; r++, start += *lx) {
    for (b=0; b<*lb; b++, cb++) {
      epsilon = boxes[b];
      // boxes always overlap by one point, thus *lx - 1
      // i.e. epsilon gives the length of the interval and
      // both marginal points are taken into account
      totalblocks = (*lx-1) / epsilon;
      idx = start;
      end = start + epsilon;
      count[cb] = 0.0;
      for (tb=0; tb<totalblocks; tb++, end+=epsilon) { 
	min = max = dat[idx];
	for (; idx<end; ) {// first value already considered!!
	  idx++;
	  if (dat[idx] < min) min = dat[idx]; else 
	    if (dat[idx] > max) max = dat[idx];
	}
	count[cb] += max - min;
      }
      count[cb] = log(count[cb] / (double) epsilon);
    } // b
  } // r
};

