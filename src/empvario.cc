/*
  Authors  Martin Schlather, schlather@math.uni-mannheim.de 

  calculation of the empirical variogram

  Copyright (C) 2002 - 2017 Martin Schlather, 

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
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111 - 1307, USA.
*/


#include <Rmath.h>  
#include <stdio.h>  
//#include <stdlib.h>
#include "RF.h"
 

// z coordinate run the fastest in values, x the slowest   

#define TOOLS_DIM 1
#define TOOLS_MEMORYERROR 2
#define TOOLS_XERROR 3
#define TOOLS_BIN_ERROR 4
#define TOOLS_UNKNOWN_CHAR 5
#define TOOLS_METHOD 6
#define METHOD_CROSS 0
#define METHOD_PSEUDO 1
#define METHOD_COVARIANCE 2
#define METHOD_PSEUDOMADOGRAM 3
#define METHOD_CROSSMADOGRAM 4



double EV_TIMES_LESS = 2.0; // 0:always new method; inf:always usual method
int EV_SYSTEM_LARGER = 100; // 0:always new method; inf:always usual method

int LOW = -1;

SEXP empiricalvariogram(SEXP  X, SEXP Dim, SEXP Lx, 
			SEXP Values, SEXP Repet, SEXP Grid, 
		        SEXP Bin, SEXP Nbin, 
			// note that within the subsequent algorithm
			// the sd of the Efct ist correctly calculated:
			//   instead for summing up Efct^2 in sd one should
			//   sum up a^2 + b^2 !
			//   so finally only NAs are returned in this case
			//   use emp
			SEXP Vdim, SEXP Method)
/*     x      : matrix of coordinates, ### rows define points
 *    dim    : dimension, 
 *    lx     : length of x, y, and z
 *    values : (univariate) data for the points
 *    repet  : number of repetitions (the calculated emprical variogram will 
 *             be an overall average)
 *    grid   : (boolean) if true lx must be 3 [ x==c(start, end, step) ]
 *    bin    : specification of the bins 
 *             sequence is checked whether it is strictly isotone
 *    nbin   : number of bins. i.e. length(bin) - 1
 *    vdim   : dimension of data
 *    method : defined above
 */
{  
  int dim = INTEGER(Dim)[0],
    lx = INTEGER(Lx)[0],
    repet = INTEGER(Repet)[0],
    grid = INTEGER(Grid)[0],
    nbin = INTEGER(Nbin)[0],
    vdim = INTEGER(Vdim)[0],
    method = INTEGER(Method)[0];
  double *x = REAL(X), *values = REAL(Values),
    *bin = REAL(Bin), *res = NULL, *sumhead = NULL
    , *sumtail = NULL, n;
  int totaln = nbin * vdim * vdim;
  int d, halfnbin, gridpoints[MAXVARIODIM], dimM1, err, 
    low, cur, up, curbin;
  long segment, totalpointsrepetvdim, totalpointsvdim,
    totalpoints = NA_INTEGER;  
  double  * xx[MAXVARIODIM], // maxdist[MAXVARIODIM],// dd, 
    * BinSq = NULL;
  SEXP Res;
  // res contains the variogram, sd and n.bin (res gets returned to function)
  // first column is of res is the variogram
  // second column is the sd and the third n.bin (number of data per bin)

  PROTECT(Res = allocMatrix(REALSXP, totaln, 3));
  res = REAL(Res);
  for(int i=0; i < totaln * 3; res[i++] = 0);

  if ( dim > MAXVARIODIM) {err = TOOLS_DIM; goto ErrorHandling; }
  if( vdim == 1) {
    // set cross to pseudo, if dimension is 1
    if(method == METHOD_CROSS) method = METHOD_PSEUDO;
    if(method == METHOD_CROSSMADOGRAM) method = METHOD_PSEUDOMADOGRAM;
  }
  
  for (int i = segment = 0; i < dim; i++, segment += lx)
    xx[i] = &(x[segment]); 
  if (xx[0]==NULL) {err=TOOLS_XERROR; goto ErrorHandling; }
  for (int i=0; i< nbin; i++) {
    if (bin[i] >= bin[i + 1])  {err = TOOLS_BIN_ERROR; goto ErrorHandling; }
  }
  
  dimM1 =  dim - 1; 
  halfnbin =  nbin / 2; 
   
  if ((BinSq = (double  * ) MALLOC(sizeof(double) * (nbin + 1)))==NULL) {
    err = TOOLS_MEMORYERROR; goto ErrorHandling; 
  }
  if(method == METHOD_COVARIANCE) {
    if( (sumhead = (double *) CALLOC(totaln, sizeof(double))) == NULL) {
      err = TOOLS_MEMORYERROR; goto ErrorHandling;
    }
    if( (sumtail = (double *) CALLOC(totaln, sizeof(double))) == NULL) {
      err = TOOLS_MEMORYERROR; goto ErrorHandling;
    }
  } 
  
  for (int i = 0; i <= nbin; i++){ 
    BinSq[i] = bin[i] > 0 ?  bin[i] * bin[i] : bin[i]; 
  }


  //////////////////////////////////// GRID ////////////////////////////////////
  if (grid) {    
    int d1, d2, head, tail; 
    double p1[MAXVARIODIM], p2[MAXVARIODIM], distSq, dx[MAXVARIODIM]; 
    long indextail[MAXVARIODIM], indexhead[MAXVARIODIM],
      segmentbase[MAXVARIODIM]; //SegmentbaseTwo[MAXVARIODIM],
      //SegmentbaseFour[MAXVARIODIM];

    // does not work!! :
    //GetGridSize(x, y, z, dim, &(gridpoints[0]), &(gridpoints[1]), &(gridpoints[2])); 
    // totalpoints = gridpoints[0] * gridpoints[1] * gridpoints[2]; 
    //
    // instead :
      //    dd = 0.0;
    totalpoints = 1;
    for (int i = 0; i <= dimM1; i++) {
      gridpoints[i] = (int) (xx[i][XLENGTH]); 
      totalpoints  *= gridpoints[i]; 
      // maxdist[i] = (gridpoints[i] - 1) * xx[i][XSTEP]; 
      //      dd += maxdist[i] * maxdist[i]; 
    }
    totalpointsvdim = totalpoints * vdim;
    
    for (int i = 0; i <= dimM1; i++) { dx[i] = xx[i][XSTEP] * 0.99999999; }
    totalpointsrepetvdim = totalpointsvdim * repet; 
    segmentbase[0] = 1; // here x runs the fastest;  
    //    SegmentbaseTwo[0] = segmentbase[0] << 1; 
    //SegmentbaseFour[0] = segmentbase[0] << 2
      ; 
    for (int i = 1; i <= dimM1; i++) { 
      segmentbase[i] = segmentbase[i - 1] * gridpoints[i - 1]; 
      //      SegmentbaseTwo[i] = segmentbase[i] << 1; 
      // SegmentbaseFour[i] = segmentbase[i] << 2; 
    }
    
    for (int i = 0; i <= dimM1; i++) {indexhead[i] = 0; p1[i] = 0.0; }
       
    // loop through all pair of points, except (head, tail) for head=tail, 
    // which is treated separately at the end, if necessary (not relevant for 
    // variogram itself, but for function E and V)
#define FOR(DO)								\
    for (head = 0; head<totalpoints; ) {				\
      for (int i = 0; i <= dimM1; i++) {				\
	indextail[i] = 0;						\
	p2[i] = 0.0;							\
      }									\
      for (tail = 0; tail<head; ) {					\
	distSq = 0.0;							\
	for (int i = 0; i <= dimM1; i++) {				\
	  double dx2;							\
	  dx2 = p1[i] - p2[i];						\
	  distSq  += dx2 * dx2;						\
	}								\
	if ((distSq>BinSq[0]) && (distSq <= BinSq[nbin])) {		\
	  /*  search which bin distSq in  */				\
	  low = 0; up = nbin; /*  21.2.01,  * nbin - 1  */		\
	  cur =  halfnbin;						\
	  while (low!=up) {						\
	    if (distSq> BinSq[cur]) {low = cur;} else {up = cur-1;}	\
	    cur = (up + low + 1) / 2;					\
	  }								\
	  for(int row = 0; row < vdim; row++) {				\
	    for(int col = row; col < vdim; col++) {			\
	      for (segment = 0; segment<totalpointsrepetvdim; segment += totalpointsvdim){ \
		curbin = low + nbin * (vdim * row + col);		\
		DO;							\
		res[curbin] += x2;					\
		res[curbin + totaln] += x2 * x2;			\
		res[curbin + 2*totaln]++;				\
	      }								\
	    }								\
	  }								\
	  assert(low < nbin);						\
	}								\
	tail++;								\
	d2 = dimM1;							\
	indextail[d2]++;						\
	p2[d2] +=dx[d2];						\
	while (indextail[d2] >= gridpoints[d2]) {			\
	  indextail[d2] = 0; p2[d2] = 0;				\
	  d2--;								\
	  assert(d2 >= 0);						\
	  indextail[d2]++;						\
	  p2[d2] +=dx[d2];						\
	}								\
      }									\
      head++;								\
      if (head<totalpoints) {						\
	d1 = dimM1;							\
	indexhead[d1]++; p1[d1] +=dx[d1];				\
	while (indexhead[d1] >= gridpoints[d1]) {			\
	  indexhead[d1] = 0;						\
	  p1[d1] = 0;							\
	  d1--;								\
	  assert(d1 >= 0);						\
	  indexhead[d1]++;						\
	  p1[d1] +=dx[d1];						\
	}								\
      }									\
    }

    switch(method){
    case(METHOD_PSEUDO):
      // pseudo variogram
      FOR(double x2 = values[head + totalpoints * row + segment]
	  - values[tail + totalpoints * col + segment]; x2 *= x2);
      break;
    case(METHOD_CROSS):
      // cross variogram
      FOR(double x2 = (values[head + totalpoints * row + segment]
		       -values[tail + totalpoints * row + segment]) *
	  (values[head + totalpoints * col + segment]-values[tail + totalpoints * col + segment]));	  
      break;
    case(METHOD_COVARIANCE):
      FOR(double x2 = (values[head + totalpoints * row + segment] *
		       values[tail + totalpoints * col + segment]);
	  sumhead[curbin] += values[head + totalpoints * row + segment];
	  sumtail[curbin] += values[tail + totalpoints * col + segment]);
      break;
    case(METHOD_PSEUDOMADOGRAM):
      // pseudo madogram
      FOR(double x2 = FABS(values[head + totalpoints * row + segment]
			   - values[tail + totalpoints * col + segment]));
      break;
    case(METHOD_CROSSMADOGRAM):
      // cross madogram
      FOR(double x2 = FABS(values[head + totalpoints * row + segment]
			   -values[tail + totalpoints * row + segment]) *
	  FABS(values[head + totalpoints * col + segment]
	       -values[tail + totalpoints * col + segment]); x2 = SQRT(x2));	  
      break;
    default:
      err = TOOLS_METHOD; goto ErrorHandling;
    }   

  } else {
    ////////////////////////////////////  ARBITRARY /////////////////////////////
    totalpoints =  lx;
    totalpointsvdim = totalpoints * vdim;
    totalpointsrepetvdim = totalpointsvdim * repet;
#define FORARB(DO)							\
    for (int i = 0; i<totalpoints; i++) { /* to have a better performance for large */ \
      /*                                 data sets, group the data first into blocks*/ \
      for (int j = 0; j<i; j++) {					\
        double distSq;							\
      	for (distSq = 0.0, d = 0; d<=dimM1; d++) {			\
	  double dx;							\
	  dx = xx[d][i] - xx[d][j];					\
	  distSq  += dx * dx;						\
	}								\
	/* see also above */						\
	/*distSq = SQRT(distSq); 26.2. */				\
	if (distSq>BinSq[0] && distSq<=BinSq[nbin]) {			\
	  low = 0;							\
	  up = nbin;							\
	  cur =  halfnbin;						\
	  while (low!=up) {						\
	    if (distSq> BinSq[cur]) low = cur; else up = cur - 1; /* ( * ; * ]*/ \
	    cur = (up + low + 1) / 2;					\
	  }								\
	  for(int row = 0; row < vdim; row++) {				\
	    for(int col = row; col < vdim; col++) {			\
	      for (segment = 0; segment<totalpointsrepetvdim; segment  += totalpointsvdim) { \
		curbin = low + nbin * (vdim * row + col);		\
		DO;							\
		res[curbin] += x2;					\
		res[curbin + totaln] += x2 * x2;			\
		res[curbin + 2*totaln]++;				\
	      }								\
	    }								\
	  }								\
	  assert(low < nbin);						\
	}								\
      }									\
    }
    
    switch(method){
    case(METHOD_PSEUDO):
      // pseudo variogram
      FORARB(double x2 = values[i + totalpoints * row + segment]
	     -values[j + totalpoints * col + segment]; x2 *= x2);
      break;
    case(METHOD_CROSS):
      // cross variogram
      FORARB(double x2 = (values[i + totalpoints * row + segment]
			  - values[j + totalpoints * row + segment])
	     * (values[i + totalpoints * col + segment]
		- values[j + totalpoints * col + segment]));	  
      break;
    case(METHOD_COVARIANCE):
      FORARB(double x2 = (values[i + totalpoints * row + segment]
			  * values[j + totalpoints * col + segment]);
	  sumhead[curbin] += values[i + totalpoints * row + segment];
	     sumtail[curbin] += values[j + totalpoints * col + segment]);
      break;
    case(METHOD_PSEUDOMADOGRAM):
      // pseudo madogram
      FORARB(double x2 = FABS(values[i + totalpoints * row + segment]
			      - values[j + totalpoints * col + segment]));
      break;
    case(METHOD_CROSSMADOGRAM):
      // cross madogram
      FORARB(double x2 = FABS(values[i + totalpoints * row + segment]
			      - values[j + totalpoints * row + segment]) *
	     FABS(values[i + totalpoints * col + segment]
		  - values[j + totalpoints * row + segment]); x2 = SQRT(x2));	  
      break;
    default:
      err = TOOLS_METHOD; goto ErrorHandling;
    }

  } // end else


  if(method == METHOD_COVARIANCE) {
    for(int row = 0; row < vdim; row++) {
      for(int col = row; col < vdim; col++) {
	curbin = nbin * (vdim * row + col);
	res[curbin] = RF_NAN;
	for(int i = 1; i < nbin; i++) {
	  n = res[i + curbin + 2*totaln];
	  double m0 = sumhead[i + curbin]/n;
	  double mh = sumtail[i + curbin]/n;
	  res[i + curbin] = res[i + curbin]/ (n) - m0 * mh;
	  res[i + curbin + totaln] = SQRT(res[i + curbin + totaln]
					  / ((n - 1)) - (n) / (n - 1) * res[i + curbin] * res[i + curbin]);
	}
      }
    }
  } else if(method == METHOD_CROSS || method == METHOD_PSEUDO) {
    for(int row = 0; row < vdim; row++) {
      for(int col = row; col < vdim; col++) {
	curbin = nbin * (vdim * row + col);
	for(int i = 1; i < nbin; i++) {
	  n = res[i + curbin + 2*totaln];
	  res[i + curbin] /= (2 * n);
	  res[i + curbin + totaln] = SQRT(res[i + curbin + totaln] / (4 * (n - 1)) - (n) /  (n - 1) * res[i + curbin] * res[i + curbin]);
	}
      }
    }
  } else {
    for(int row = 0; row < vdim; row++) {
      for(int col = row; col < vdim; col++) {
	curbin = nbin * (vdim * row + col);
	for(int i = 1; i < nbin; i++) {
	  n = res[i + curbin + 2*totaln];
	  res[i + curbin] /= (n);
	  res[i + curbin + totaln] = 0;
	}
      }
    }
  }

  for(int row = 1; row < vdim; row++) {
    for(int col = 0; col < row; col++) {
      int noncalcpos = nbin * (vdim * row + col);
      int calcpos = nbin * (row + vdim * col);
      if (method == METHOD_COVARIANCE) res[noncalcpos] = RF_NAN;
      for(int i = 1; i < nbin; i++) {
	res[i + noncalcpos] = res[i + calcpos];
	res[i + noncalcpos + totaln] = res[i + calcpos + totaln];
	res[i + noncalcpos + 2*totaln] = res[i + calcpos + 2*totaln];
      }
    }
  }
  
  UNPROTECT(1);
  FREE(BinSq);
  if(method == METHOD_COVARIANCE) {
    FREE(sumhead);
    FREE(sumtail);
  }
  return Res;
  
 ErrorHandling:
  FREE(BinSq);
  if(method == METHOD_COVARIANCE) {
    FREE(sumhead);
    FREE(sumtail);
  }
  UNPROTECT(1);
  switch (err) {
  case TOOLS_DIM :
    ERR("dimension exceed max dimension of empirical variogram estimation"); 
  case TOOLS_MEMORYERROR :  
    ERR("Memory alloc failed in empiricalvariogram.\n"); 
  case TOOLS_XERROR :  
    ERR("The x coordinate may not be NULL.\n"); 
  case TOOLS_BIN_ERROR :
    ERR("Bin components not an increasing sequence.\n"); 
  case TOOLS_UNKNOWN_CHAR:
    ERR("unknown type of second order characteristic");
  case TOOLS_METHOD:
    BUG;
  default : BUG;
  }
}
