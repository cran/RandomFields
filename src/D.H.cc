/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 calculation of the Hurst coefficient and the fractal dimension of
 the graph of a random field

 Copyright (C) 2002 - 2015 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
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
#include <Rdefines.h>

#include "RF.h"

SEXP boxcounting(SEXP Z,  // data
		 SEXP LX,    // basic length of data 
		 SEXP Repet, // "independent" repetitions of data 
		 SEXP Factor, // rescalin of data values by factor
		 SEXP Eps   // box lengths
    ) {
  // in sum, the boxcounting result for the first lx data values is given 
  // according to the box lengths given by eps; then the boxes are counted
  // for the next lx data values, etc. Hence, sum must have dimension
  int 
    *eps = INTEGER(Eps),
    leps = length(Eps),
    repet = INTEGER(Repet)[0],
    lx = INTEGER(LX)[0];
  long i,  k,  j, e,   r, lastbox, s, // 0 to leps * repet -1
    truelx = 2 + lx, 
    total = repet * truelx;
  double zz, min, max, f, *sum,
    factor = REAL(Factor)[0],
    *z = REAL(Z); 

  SEXP Sum;
  PROTECT(Sum = allocVector(REALSXP, leps * repet));
  sum = REAL(Sum);

  // boxcounting expects that the first and the last values are repeated
  // in the first dimension, see D.H.R; accordingly the total length in
  // the first direction is increased by 2:
  for (s=r=0; r<total; r+=truelx) {
    for (k=0; k<leps; k++, s++) {
      sum[s] = 0;
      e = eps[k]; // box size
      assert(e <= lx); 
      lastbox = r + (e * (lx / e)) + 1 - e;
      for (i = 1 + r; i <= lastbox; ) {
	// # boxen um  Daten von 0.5 * (i + i-1) bis 0.5 * (i + i+1)
	// zu ueberdecken
	min = max = 0.5 * (z[i] + z[i-1]);
	for (j=0; j<e; j++, i++) 
	  if (z[i] < min) min = z[i]; else if (z[i] > max) max = z[i];
	assert(i <= r + lx + 1);
	if ((zz = 0.5 * (z[i-1] + z[i])) < min) min = zz; 
	else if (zz > max) max = zz;
	f = factor / (double) e;
	sum[s] += floor(max * f) - floor(min * f) + 1.0;
      }
    }
  }
  UNPROTECT(1);
  return Sum;
}

SEXP periodogram(SEXP Dat, // data
		 SEXP Len,    // basic length of data
		 SEXP Repet,  // # of "independent" repetitions
		 SEXP FFTM,  // fftm[0] >=1 not >=1 (R indexing assumed)
		 //             interval of the vector given by fft of part
		 SEXP Part,  // not the basic length but only a part is
                 //             transformed
		 SEXP Shift // part is then shifted by shift until end of
		 //            vector is reached; fft values are averaged
                 //             over all these parts
    ){ 
  int *fftm = INTEGER(FFTM),
    part = INTEGER(Part)[0];
  assert((fftm[0]>=1) && (fftm[1] <=part) && (fftm[1] >= fftm[0]));
  long seg_dat,  k, j, total_seg, segment_r, 
    repet = INTEGER(Repet)[0],
    start_k = 2 * (fftm[0] - 1), // of complex numbers
    end_k = 2 * (fftm[1] - 1),  // of complex numbers
    segm_l,   r, 
    delta_l = fftm[1] - fftm[0] + 1,
    end_l = delta_l * repet,
    err = ERRORFAILED,
    len = INTEGER(Len)[0],
    lenMpart = len - part,
    shift = INTEGER(Shift)[0];
  assert(lenMpart >= 0);
  assert(shift>=1);
 
  double 
    taper_fact = sqrt( 2.0 / (3.0 * ((double) part + 1.0)) ),
    cos_factor = 2.0 * PI / ((double) part + 1.0),
    *compl_number = NULL, 
    *taper = NULL,
    *lambda = NULL,
     n_inv = 1.0 / (double) ((int) (1.0 + lenMpart / (double) shift)),
    *dat = REAL(Dat),
    factor = log(2.0 * PI * len)
    ; 

  FFT_storage FFT;
  SEXP Lambda;
  
  PROTECT(Lambda = allocVector(REALSXP, end_l)); 
  lambda = REAL(Lambda);
  for (j=0; j<end_l; j++) lambda[j]=0.0;
  
  FFT_NULL(&FFT);
  if ((compl_number = (double*) MALLOC(sizeof(double) * 2 * part))==NULL || 
      (taper = (double*) MALLOC(sizeof(double) * part))==NULL){
    goto ErrorHandling;
  }
  
  for (j=0; j<part; j++) {
    // taper is used to reduce bias for high frequencies by smoothing start 
    // and end of series towards zero
    // However a high bias for the smallest one or two frequencies is introduced 
    taper[j] = taper_fact * ( 1.0 - cos (cos_factor * (double) (j + 1)));
  }

  for (segment_r=segm_l=r=0; r<repet; r++, segment_r+=len, segm_l+=delta_l) {
    // repet
    for (seg_dat=0; seg_dat<=lenMpart; seg_dat+=shift) {
      // shifts
      total_seg = seg_dat + segment_r;
      for (k=0, j=0; j<part ; j++) {
	compl_number[k++] = dat[total_seg + j] * taper[j];
	compl_number[k++] = 0.0;
      }
      if ((err=fastfourier(compl_number,
			   &part,   // length
			      1,      // dim 
			      (total_seg==0), // start? 
			      false,  // inverse?  
			   &FFT))!=0) {    // aux 
	goto ErrorHandling;
      }
      for (j=segm_l, k=start_k; k<end_k; j++, k+=2) {
	lambda[j] += log(compl_number[k] * compl_number[k] + 
			 compl_number[k+1] * compl_number[k+1]) - factor;
      }
    }
  }
 

 ErrorHandling:
  FREE(compl_number);
  FREE(taper);
  FFT_destruct(&FFT);
  UNPROTECT(1);
  if (err != NOERROR) ERR("error occured when calculating the periodogram");
  for (j=0; j<end_l; j++) lambda[j] *= n_inv;
  return Lambda;
}

// both detrended fluctuation and aggreated variation  
SEXP detrendedfluc(SEXP Dat, // data
		   SEXP LX,     // basic length of data
		   SEXP Repet,  // # of "independent" repetitions
		   SEXP Blocks, // block (box) lengths
		   SEXP Ldfa   // # of blocks
		   ){
  int *blocks = INTEGER(Blocks);
  long i, j, m, nbox, ex, endfor, r,  last, e, k,
    lx = INTEGER(LX)[0],
    repet = INTEGER(Repet)[0],
    ldfa = INTEGER(Ldfa)[0],
    total = lx * repet
    ;
   assert(lx > 1);
 double var, a, b, residual, Yt, realm, sumt, meanY, t, 
    VarMeth_old, VarMeth_mean, VarMeth_var, delta, realnbox, *lvar,
    *dat = REAL(Dat);
  SEXP Lvar; // ## rbind(l.varmeth.var, l.dfa.var) 2 * (repet * ldfa)
  PROTECT(Lvar = allocMatrix(REALSXP, 2, ldfa * repet));
  lvar = REAL(Lvar);
  
  for (e = r = 0; r < total; r += lx, e+=ldfa) {
  last = r + lx;
  for (j=r+1; j<last; j++) { dat[j] += dat[j-1]; }
  for (i=0; i<ldfa; i++){
  int idx = 2 *(e + i);
    var = 0.0;
     m = blocks[i];
      realm = (double) m;
      nbox = lx / m;      
      realnbox = (double) nbox;
      ex = r + m * nbox;
      sumt = 0.5 * realm * (realm + 1.0); // sum_i^m i

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
        lvar[idx] = log(VarMeth_var / ((realnbox - 1.0)));
      } else {
	lvar[idx] = RF_NA;
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
      lvar[idx + 1] = log(var / ((realm - 1.0) * (double) nbox));
    }  
  }
  UNPROTECT(1);
  return Lvar;
}


// range method
#define FRACT_MAXDIM 10
SEXP minmax(SEXP Dat,  // data
    SEXP LX,      // basic length
    SEXP Repet,   // repetitions
    SEXP Boxes,   // box lengths
    SEXP LB      // # of box lengths
    ) {
  int b, tb, cb, r, totalblocks, epsilon, idx, start, end,
    lx = INTEGER(LX)[0],
    repet = INTEGER(Repet)[0],
    *boxes = INTEGER(Boxes),
    lb = INTEGER(LB)[0];
  double min, max, *count,
    *dat = REAL(Dat);
  SEXP Count;
  PROTECT(Count = allocVector(REALSXP, repet * lb));
  count = REAL(Count);

  for (start=cb=r=0; r<repet; r++, start += lx) {
    for (b=0; b<lb; b++, cb++) {
      epsilon = boxes[b];
      // boxes always overlap by one point, thus lx - 1
      // i.e. epsilon gives the length of the interval and
      // both marginal points are taken into account
      totalblocks = (lx-1) / epsilon;
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

  UNPROTECT(1);
  return Count;
}
