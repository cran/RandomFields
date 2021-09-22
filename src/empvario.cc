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


#include <stdio.h>  
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "RF.h"
 

// z coordinate run the fastest in values, x the slowest   


#define LOOP(DOIT)						\
  for(int row = 0; row < vdim; row++) {				\
    for(int col = row; col < vdim; col++) {			\
      int curbin = nbin * (vdim * row + col),			\
	other = (row - col) * nbinvdim,				\
	endfor = nbin + curbin;					\
      for(int i = curbin; i < endfor; i++) {			\
	int i_other = i + other;				\
	double n = res[i_other + nidx] = res[i + nidx];		\
	DOIT;							\
	res[i_other + evidx] = res[i + evidx];			\
	res[i_other + sdSqidx] = res[i + sdSqidx];		\
      }								\
    }								\
  }


void calculate_means(int method, int vdim, int nbin, int totaln,
		     double *sumhead, double *sumtail, double *res) {
  int nbinvdim = nbin * (1.0 - vdim),
    nidx = totaln * EV_N,
    evidx = totaln * EV_EV,
    sdSqidx = totaln * EV_SDSQ;
  switch(method){
  case(METHOD_PSEUDO): case(METHOD_VARIOGRAM):
    LOOP(res[i + sdSqidx] = 0.25 * ( res[i + sdSqidx] / (n - 1.0) -
   		      res[i + evidx] * res[i + evidx] / (n * (n - 1.0)));
	 res[i + evidx] /= (2.0 * n)
	 );
    break;
  case(METHOD_COVARIANCE):
    LOOP(double m0 = sumhead[i] / n;
	 double mh = sumtail[i] / n;
	 // Reihenfolge!!
	 res[i + sdSqidx] = (res[i + sdSqidx] / (n - 1) -
	 		    res[i + evidx] * res[i + evidx] / (n * (n - 1)));
	 res[i + evidx] = res[i + evidx] / n - m0 * mh
	 );
    break;
  case(METHOD_PSEUDOMADOGRAM): case(METHOD_MADOGRAM):
    LOOP(res[i + sdSqidx] = (res[i + sdSqidx] / (n - 1) -
			    res[i + evidx] * res[i + evidx] / (n * (n - 1)));
	 res[i + evidx] /= n;
	 );
      break;
  default:
    PRINTF("calculate_means:\n");
    XERR(TOOLS_METHOD);
  }
}

 

SEXP empirical(SEXP  X, SEXP Dim, SEXP Lx, 
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
  double 
    *x = REAL(X),
    *values = REAL(Values),
    *bin = REAL(Bin),
    *res = NULL,
    *sumhead = NULL,
    *sumtail = NULL;
  int totaln = nbin * vdim * vdim,
    nidx = totaln * EV_N,
    evidx = totaln * EV_EV,
    sdSqidx = totaln * EV_SDSQ;
  int d, halfnbin, gridpoints[MAXVARIODIM], dimM1, 
    low, cur, up, 
    err = NOERROR;
  long segment, totalpointsrepetvdim, totalpointsvdim,
    totalpoints = NA_INTEGER;  
  double  * xx[MAXVARIODIM], // maxdist[MAXVARIODIM],// dd, 
    * BinSq = NULL;
  SEXP Res;
  // res contains the variogram (etc), variance of the variogram, and n.bin (res gets returned to function)
  // first column is of res is the variogram (etc)
  // second column is the variance of the variogram (etc)
  // third column n.bin (number of data per bin)

  PROTECT(Res = allocMatrix(REALSXP, totaln, 3));
  res = REAL(Res);
  for(int i=0; i < totaln * 3; res[i++] = 0);

  if ( dim > MAXVARIODIM || dim <= 0) {err = TOOLS_DIM; goto ErrorHandling; }
  if( vdim == 1) {
    // set cross to pseudo, if dimension is 1
    if(method == METHOD_VARIOGRAM) method = METHOD_PSEUDO;
    if(method == METHOD_MADOGRAM) method = METHOD_PSEUDOMADOGRAM;
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
    err = ERRORMEMORYALLOCATION; goto ErrorHandling; 
  }
  if(method == METHOD_COVARIANCE) {
    if( (sumhead = (double *) CALLOC(totaln, sizeof(double))) == NULL) {
      err = ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    if( (sumtail = (double *) CALLOC(totaln, sizeof(double))) == NULL) {
      err = ERRORMEMORYALLOCATION; goto ErrorHandling;
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
		int curbin = low + nbin * (vdim * row + col);		\
		DO;							\
		res[curbin + evidx] += x2;				\
		res[curbin + sdSqidx] += x2 * x2;			\
		res[curbin + nidx]++;				\
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
    case(METHOD_VARIOGRAM):
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
    case(METHOD_MADOGRAM):
      // cross madogram
      FOR(double x2 = FABS(values[head + totalpoints * row + segment]
			   -values[tail + totalpoints * row + segment]) *
	  FABS(values[head + totalpoints * col + segment]
	       -values[tail + totalpoints * col + segment]); x2 = SQRT(x2));	  
      break;
    default:
     PRINTF("empirical:\n");
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
	    if (distSq> BinSq[cur]) low = cur;				\
	    else up = cur - 1; /* ( * ; * ]*/				\
	    cur = (up + low + 1) / 2;					\
	  }								\
	  for(int row = 0; row < vdim; row++) {				\
	    for(int col = row; col < vdim; col++) {			\
	      for (segment = 0; segment<totalpointsrepetvdim;		\
		   segment  += totalpointsvdim) {			\
		int curbin = low + nbin * (vdim * row + col);		\
		DO;							\
		res[curbin + evidx] += x2;				\
		res[curbin + sdSqidx] += x2 * x2;			\
		res[curbin + nidx]++;					\
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
    case(METHOD_VARIOGRAM):
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
    case(METHOD_MADOGRAM):
      // cross madogram
      FORARB(double x2 = FABS(values[i + totalpoints * row + segment]
			      - values[j + totalpoints * row + segment]) *
	     FABS(values[i + totalpoints * col + segment]
		  - values[j + totalpoints * row + segment]); x2 = SQRT(x2));	  
      break;
    default:
      PRINTF("final emp. vario: %d\n", method);
      err = TOOLS_METHOD; goto ErrorHandling;
    }

  } // end else

  calculate_means(method, vdim, nbin, totaln, sumhead, sumtail, res);

  
 ErrorHandling:
  FREE(BinSq);
  if(method == METHOD_COVARIANCE) {
    FREE(sumhead);
    FREE(sumtail);
  }
  UNPROTECT(1);
  if (err == NOERROR)
    return(Res);

  PRINTF("In empirical variogram an error happened:");
  XERR(err);
}


SEXP empvarioXT(SEXP Xsexp, SEXP Tsexp, 
		SEXP Lx, 
	        SEXP Values, SEXP Repet, SEXP Grid,
		SEXP Bin, SEXP Nbin, 
		SEXP Phi,    // vector of a real and an integer
		SEXP Theta,  // vector of a real and an integer
		SEXP DT,   // in grid units, smallest positive value, for  
		//            the temporal lag of the emp. variogram
	        //            stepT * nstepT is the greatest time span
		SEXP Vdim, SEXP Method)
/*    X      : matrix of coordinates, ### rows define points
 *    lx     : length of x,y, and z
 *    values : (univariate) data for the points
 *    repet  : number of repetitions (the calculated emprical variogram will 
 *             be an overall average, see also EmpiricalVariogram in empvario.R)
 *    grid   : (boolean) if true lx must be 3 [ x==c(start,end,step) ]
 *    bin    : specification of the bins 
 *             sequence is checked whether it is strictly isotone
 *    nbin   : number of bins. i.e. length(bin)-1
 *    sum    : empirical variogram
 *    n      : number of pairs for each bin
 *    vdim   : dimension of data
 *    Method : 0 cross-variogram, 1 pseudo-variogram, 2 covariance
 */
{ 
 // turning SEXP to c++ variable
  int lx = INTEGER(Lx)[0], repet = INTEGER(Repet)[0],
    grid = INTEGER(Grid)[0],
    nbin = INTEGER(Nbin)[0],
    *dT = INTEGER(DT), vdim = INTEGER(Vdim)[0],
    method = INTEGER(Method)[0];
  double *X = REAL(Xsexp), *T = REAL(Tsexp),
    *values = REAL(Values), *bin = REAL(Bin),
    *phi = REAL(Phi), *theta = REAL(Theta), *res = NULL;


  long rep;
  int i ,gridpoints[4], totalbins, totalspatialbins,
    twoNphi, twoNphiNbin, Ntheta, nbinNphiNtheta, maxi[4], mini[4], ix, iy, iz,
    it, endX, startX, endY, startY, endZ, startZ, endT, low, cur, up, totaln,
    nidx, evidx, sdSqidx,
    ktheta, kphi, nstepTstepT, vec[4], deltaT, x, y, z, t,  j, stepT, nstepT,
    halfnbin = nbin / 2,
    err  = NOERROR; 
  double *xx[4], *BinSq, maxbinsquare, step[4], delta[4], psq[4], dx[4],
    startphi, invsegphi, starttheta, invsegtheta, thetadata, phidata,
    zylinderradius,
    *sumhead = NULL,
    *sumtail = NULL;
  SEXP Res; // res contains the variogram, sd and n.bin 
  // first column is of res is the variogram
  // second column is the sd^2 and the third n.bin (number of data per bin)

  
  stepT = dT[0];
  nstepT = dT[1];
  nstepTstepT = stepT * nstepT;
  twoNphi =  phi[1]==0 ?  1 : (2 * (int) phi[1]);
  twoNphiNbin = twoNphi * nbin;
  Ntheta = theta[1]==0 ? 1 : (int) theta[1];
  startphi = phi[0] - PI / (double) twoNphi; // [0, 2 pi]
  invsegphi = phi[1] / PI; // note that phi[1] can be zero!
  //starttheta = theta[0] - PIHALF / (double) Ntheta;  // [0, pi]
  starttheta = theta[0] - PIHALF;  // [0, pi]
  invsegtheta = theta[1] / PI; // note that theta[1] can be zero
  nbinNphiNtheta = twoNphiNbin * Ntheta;
  totalspatialbins =  twoNphiNbin * Ntheta;
  totalbins = totalspatialbins * (nstepT + 1);
  totaln = totalbins * vdim * vdim;
  nidx = totaln * EV_N;
  evidx = totaln * EV_EV;
  sdSqidx = totaln * EV_SDSQ;
   
  PROTECT(Res = allocMatrix(REALSXP, totaln, 3));
  res = REAL(Res);

  BinSq = NULL;
  for (rep=i=0; i<3; i++, rep += lx) xx[i]=&(X[rep]);
  xx[3] = T;

  if (xx[0]==NULL) {err=TOOLS_XERROR; goto ErrorHandling;}
  for (i=0; i < nbin; i++) {
    if (bin[i]>=bin[i+1])  {err=TOOLS_BIN_ERROR; goto ErrorHandling;}
  }
  if(vdim == 1) {
    if(method == METHOD_VARIOGRAM) method = METHOD_PSEUDO;
    if(method == METHOD_MADOGRAM) method = METHOD_PSEUDOMADOGRAM;
  }
   
  if ((BinSq = (double *) MALLOC(sizeof(double)* (nbin + 1)))==NULL) {
    err=ERRORMEMORYALLOCATION; goto ErrorHandling; 
  } 

  if(method == METHOD_COVARIANCE) {
    if( (sumhead = (double *) CALLOC(totaln, sizeof(double))) == NULL) {
      err = ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    if( (sumtail = (double *) CALLOC(totaln, sizeof(double))) == NULL) {
      err = ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
  }

  // res contains the variogram, sd and n.bin (res gets returned to function)
  // first column is of res is the variogram
  // second column is the sd^2 and the third n.bin (number of data per bin)
  //printf("total %d %d %d %d \n", nstepT,  twoNphi, *nbin,  Ntheta);

  for(i = 0; i < totaln * 3; res[i++] = 0); 
  for (i=0; i<= nbin; i++){if (bin[i]>0) BinSq[i]=bin[i] * bin[i]; 
    else BinSq[i]=bin[i];
  }

  assert(NEARBY(atan2(-1.0, 0.0) + PIHALF) == 0.0);
  assert(atan2(0.0, 0.0) == 0.0);
  maxbinsquare = BinSq[nbin];
 

  //////////////////////////////////// GRID ////////////////////////////////////
  if (grid) {
    int segmentbase[6];
    long totalpoints = 1;
    segmentbase[0]=1; // here x runs the fastest;
    
    for (i=0; i<=3; i++) {
      step[i] = xx[i][XSTEP];
      gridpoints[i] = (int) xx[i][XLENGTH];
      totalpoints *= (gridpoints[i] > 0) ? gridpoints[i] : 1;
      maxi[i] = (int) (SQRT(maxbinsquare) / step[i] + 0.1);
      if (maxi[i] >= gridpoints[i]) maxi[i] = gridpoints[i] - 1;
      mini[i] = -maxi[i];
      //      if (maxi[i]==0)  maxi[i] = 1; wrong 22.11.03
      segmentbase[i+1] =  segmentbase[i] * gridpoints[i];
    }
    segmentbase[4] *= vdim;
    segmentbase[5] = segmentbase[4] * repet;

    //printf("\n");
    //for(int i=totalpoints;i<2*totalpoints;i+=100){
    //  printf("%10g ", values[i]);
    //}
    //printf("\n");
    // sementbase: [0] : 1, [1] : length(x), [2] : length(x) * l(y),
    //  [3] : l(x) * l(y) * l(z), [4] : l(x) * l(y) * l(z) * l(T)
    //  [5] : l(x) * l(y) * l(z) * l(T) * repet

    // define needed so there is no if needed before each calculation (grid)
#define FORXT(DO)								\
    for (ix=mini[0]; ix <= maxi[0]; ix++) {				\
      delta[0] = ix * step[0];						\
      if ((psq[0] = delta[0] * delta[0]) > maxbinsquare) continue;	\
									\
      vec[0] = ix;							\
      if (ix>=0) {							\
	endX = gridpoints[0] - ix;					\
	startX = 0;							\
      } else {								\
	endX = gridpoints[0];						\
	startX = -ix;							\
      }									\
      for (iy= mini[1]; iy <= maxi[1];  iy++) {				\
	delta[1] = iy * step[1];					\
	if ((psq[1] = delta[1] * delta[1] + psq[0]) > maxbinsquare) continue; \
	vec[1] = vec[0] + iy * segmentbase[1];				\
									\
        /* angles */							\
	phidata = NEARBY(atan2(delta[1], delta[0]) - startphi);		\
	phidata = Mod(phidata, TWOPI);					\
	kphi = nbin * (int) (phidata * invsegphi);			\
									\
									\
	for (iz=mini[2]; iz <= maxi[2]; iz++) {				\
	  delta[2] = iz * step[2];					\
	  if ((psq[2] = delta[2] * delta[2] + psq[1]) > maxbinsquare) continue;\
	  vec[2] = vec[1] + iz * segmentbase[2];			\
									\
	  {								\
	    low=0; up= nbin; /* */ cur= halfnbin;			\
	    while(low!=up){						\
	      if (psq[2] > BinSq[cur]) {low=cur;} else {up=cur-1;}/*( . ; . ]*/\
	      cur=(up+low+1)/2;						\
	    }								\
	  }	/* low*/						\
									\
	  /* angles*/							\
	  thetadata = NEARBY(atan2(delta[2],SQRT(psq[1]))-starttheta);	\
	  thetadata = Mod(thetadata, PI);				\
	  ktheta = (int) (thetadata * invsegtheta);			\
									\
	  low += kphi + ktheta * twoNphiNbin;				\
									\
	  for (deltaT=0, it=low; deltaT<=nstepTstepT;			\
	       deltaT+=stepT, it+=nbinNphiNtheta) {			\
	    vec[3] = vec[2] + deltaT * segmentbase[3];			\
									\
	    assert(startX>=0);						\
	    for (x=startX; x<endX; x++) {				\
	      if (iy>=0) {						\
		endY = gridpoints[1]  * segmentbase[1] - vec[1];	\
		startY = x;						\
	      } else {							\
		endY = gridpoints[1] * segmentbase[1];			\
		startY = x - iy * segmentbase[1];			\
	      }								\
	      assert(startY>=0);					\
	      for (y=startY; y<endY; y+=segmentbase[1]) {		\
		if (iz>=0) {						\
		  endZ = gridpoints[2] * segmentbase[2] - vec[2];	\
		  startZ = y;						\
		} else {						\
		  endZ = gridpoints[2] * segmentbase[2];		\
		  startZ = y - iz * segmentbase[2];			\
		}							\
                assert(startZ>=0);					\
		for (z=startZ; z<endZ; z+=segmentbase[2]) {		\
		  endT = gridpoints[3] * segmentbase[3] - vec[3];	\
		  for (t=z; t<endT; t+=segmentbase[3]) {		\
		    for(int row = 0; row < vdim; row++) {		\
		      for(int col = row; col < vdim; col++) {		\
			int curbin = it + totalbins * (vdim * row + col);\
			for (rep = t; rep < segmentbase[5];		\
			     rep += segmentbase[4]) {			\
			  int rv;					\
			  rv = rep + vec[3];				\
			  assert(rv < segmentbase[5] && rv >=0);	\
			  /* inserts calculation from switch below */	\
			  DO;						\
			  res[curbin + evidx] += x2;			\
			  res[curbin + sdSqidx] += x2 * x2;		\
			  res[curbin + nidx]++;				\
			}						\
		      }							\
		    } /* repeat*/					\
		  } /* t */						\
		} /* z */						\
	      } /* y */							\
	    } /* x */							\
	  } /* deltaT */						\
	} /* iz */							\
      } /* iy */							\
    } /* ix */

    // checks which method shoud be used (gets inserted above)
   
    switch(method) {
    case(METHOD_PSEUDO):
      // pseudo variogram
      FORXT(double x2 = values[rep + totalpoints * row]
	  - values[rv + totalpoints * col];
	  x2 *= x2);
      break;
    case(METHOD_VARIOGRAM):
      // cross variogram
      FORXT(double x2 = (values[rep + totalpoints * row]
		       - values[rv + totalpoints * row])
	  * (values[rep + totalpoints * col]
	     - values[rv + totalpoints * col]));	  
      break;
    case(METHOD_COVARIANCE):
      FORXT(double x2 = (values[rep + totalpoints * row]
		       * values[rv + totalpoints * col]);
	  sumhead[curbin] += values[rep + totalpoints * row];
	  sumtail[curbin] += values[rv + totalpoints * col]);
      break;
    case(METHOD_PSEUDOMADOGRAM):
      // pseudo madogram
      FORXT(double x2 = FABS(values[rep + totalpoints * row]
			   - values[rv + totalpoints * col]));
      break;
    case(METHOD_MADOGRAM):
      // cross madogram
      FORXT(double x2 = FABS(values[rep + totalpoints * row]
			   - values[rv + totalpoints * row])
	  * FABS(values[rep + totalpoints * col]
		 - values[rv + totalpoints * col]); x2 = SQRT(x2));	  
      break;
    default:
     PRINTF("emvarioXT:\n");
     err = TOOLS_METHOD; goto ErrorHandling;
    }
    
  } else {
    ////////////////////////////////////  ARBITRARY /////////////////////////////
    // rows : x, columns : T 
    long totalpoints, totalpointsrepetvdim, totalpointsvdim, spatial,  jj;
    spatial = lx;
    i = 3;
    step[i] = xx[i][XSTEP];
    gridpoints[i] = (int) xx[i][XLENGTH];
    totalpoints = spatial * gridpoints[i];
    totalpointsvdim = totalpoints * vdim;
    totalpointsrepetvdim = totalpoints * repet;
    
    // define needed so there is no if needed before each calculation (arbitrary)
#define FORARBXT(DO)							\
    for (i=0;i<spatial;i++) { /* to have a better performance for large*/ \
      /*                  data sets, group the data first into blocks*/	\
      for (j=0; j<spatial; j++) {					\
	double distSq;							\
	dx[0] = xx[0][j] - xx[0][i];					\
	dx[1] = xx[1][j] - xx[1][i];					\
	dx[2] = xx[2][j] - xx[2][i];					\
	distSq = dx[0] * dx[0] + dx[1] * dx[1];				\
	zylinderradius = SQRT(distSq);					\
	distSq += dx[2] * dx[2];					\
									\
	if ((distSq>BinSq[0]) && (distSq<=BinSq[nbin])) {		\
	  low=0; up=nbin; cur=halfnbin;					\
	  while (low!=up) {						\
	    if (distSq> BinSq[cur]) {low=cur;} else {up=cur-1;} /* ( * ; * ]*/ \
	    cur=(up+low+1)/2;						\
	  }								\
									\
	  /* angles */							\
	  phidata = NEARBY(atan2(dx[1], dx[0]) - startphi);		\
	  phidata = Mod(phidata, TWOPI);				\
	  kphi = nbin * (int) (phidata * invsegphi);			\
									\
	  thetadata = NEARBY(atan2(dx[2], zylinderradius) - starttheta); \
	  thetadata = Mod(thetadata, PI);				\
	  ktheta = (int) (thetadata * invsegtheta);			\
									\
	  low += kphi + twoNphiNbin * ktheta;				\
	  assert(low>=0 && low<totalspatialbins);			\
									\
	  for (deltaT=0; deltaT<= nstepTstepT; deltaT+=stepT,		\
		 low+=totalspatialbins) {				\
	    jj = j + deltaT * spatial;					\
	    endT= totalpoints - deltaT * spatial;			\
	    for (t=0; t<endT; t+=spatial) {				\
	      for(int row = 0; row < vdim; row++) {			\
		for(int col = row; col < vdim; col++) {			\
		  for (rep = t; rep < totalpointsrepetvdim;		\
		       rep += totalpointsvdim){				\
		    int curbin = 0;					\
		    /* inserts calculation from switch below */		\
		    curbin = low + totalbins * (vdim * row + col);	\
		    DO;							\
		    res[curbin + evidx] += x2;				\
		    res[curbin + sdSqidx] += x2 * x2;			\
		    res[curbin + nidx]++;				\
		  }							\
		}							\
	      }								\
	    } /* t*/							\
	  } /* deltat*/							\
	} /* if (distSq)*/						\
      } /* j */								\
    } /* i */
   
    // checks which method shoud be used (gets inserted above) 
    switch(method) {
    case(METHOD_PSEUDO):
      // pseudo variogram
      FORARBXT(double x2 = values[rep + jj + totalpoints * row]
	     - values[i + rep + totalpoints * col];
	     x2 *= x2);
      break;
    case(METHOD_VARIOGRAM):
      // cross variogram
      FORARBXT(double x2 = (values[rep + jj + totalpoints * row]
			  - values[i + rep + totalpoints * row])
	     * (values[rep + jj + totalpoints * col]
		- values[i + rep + totalpoints * col]));	  
      break;
    case(METHOD_COVARIANCE):
      FORARBXT(double x2 = (values[rep + jj + totalpoints * row]
			  * values[i + rep + totalpoints * col]);
	     sumhead[curbin] += values[rep + jj + totalpoints * row];
	     sumtail[curbin] += values[i + rep  + totalpoints * col]);
      break;
    case(METHOD_PSEUDOMADOGRAM):
      // pseudo madogram
      FORARBXT(double x2 = FABS(values[rep + jj + totalpoints * row]
			      - values[i + rep + totalpoints * col]));
      break;
    case(METHOD_MADOGRAM):
      // cross madogram
      FORARBXT(double x2 = FABS(values[rep + jj + totalpoints * row]
				- values[i + rep + totalpoints * row])
	       * FABS(values[rep + jj + totalpoints * col]
		      - values[i + rep + totalpoints * col]); x2 = SQRT(x2));	  
      break;
    default:
     PRINTF("final emvarioXT:\n");
      err = TOOLS_METHOD; goto ErrorHandling;
    }

    
  } /* arbitrary */

  calculate_means(method, vdim, totalbins, totaln, sumhead, sumtail, res);

 ErrorHandling:
  FREE(BinSq);
  if(method == METHOD_COVARIANCE) {
    FREE(sumhead);
    FREE(sumtail);
  }
  UNPROTECT(1);
  
  if (err == NOERROR)
    return(Res);

  PRINTF("When calculating the  empirical variogram an error happened:");
  XERR(err);
} 
	 



