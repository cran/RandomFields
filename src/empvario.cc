/*
 Authors  Martin Schlather, schlather@math.uni-mannheim.de 

 calculation of the empirical variogram

 Copyright (C) 2002 - 2014 Martin Schlather, 

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


#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include "RF.h"
 

// z coordinate run the fastest in values, x the slowest   

#define TOOLS_DIM 1
#define TOOLS_MEMORYERROR 2
#define TOOLS_XERROR 3
#define TOOLS_BIN_ERROR 4
#define TOOLS_UNKNOWN_CHAR 5


double EV_TIMES_LESS = 2.0; // 0:always new method; inf:always usual method
int EV_SYSTEM_LARGER = 100; // 0:always new method; inf:always usual method

double variogram(double a, double b) {
  double dummy = a - b; 
  return dummy * dummy; 
}

double Efunction(double a, double b) {return a + b; } 
/*   see Schlather (2001), Bernoulli, 7(1), 99 - 117
    Schlather (????), Math. Nachr, conditionally accepted
    Schlather, Ribeiro, Diggle (in prep.) 
 */

double kmm(double a, double b) {return a * b; }  /*  
					     Stoyan's k_mm function
					     see Stoyan, Kendall, & Mecke, 
					     or Stoyan & Stoyan (both Wiley) 
					   */

int LOW = -1; 

void empiricalvariogram(double  * x, int  * dim, int  * lx, 
			double  * values, int  * repet, int  * grid, 
			double  * bin, int  * nbin, int  * charact, 
			double  * vario, 
			double  * sd, // 
			// note that within the subsequent algorithm
			// the sd of the Efct ist correctly calculated:
			//   instead for summing up Efct^2 in sd one should
			//   sum up a^2 + b^2 !
			//   so finally only NAs are returned in this case
			//   use emp
			int  * n)
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
 *    character: 1 : function E
 *               2 : Stoyan's kmm function
 *               default : semivariogram
 *    vario    : empirical variogram
 *    n      : number of pairs for each bin
  */
{ 
  int i, ii, j, d, halfnbin, gridpoints[MAXVARIODIM], dimM1, EV_METHOD, err, 
    low, cur, up;    
  long totalpoints, totalpointsrepet, segment, 
    factor =  -1; 
  double ( * characteristic)(double, double) = NULL; 
  double  * xx[MAXVARIODIM], maxdist[MAXVARIODIM], dd, 
     * BinSq = NULL; 

  if ( * dim > MAXVARIODIM) {err = TOOLS_DIM; goto ErrorHandling; }
    
  for (segment = i = 0; i< * dim; i++, segment += * lx) xx[i] = &(x[segment]); 
  if (xx[0]==NULL) {err=TOOLS_XERROR; goto ErrorHandling; }
  for (i=0; i< * nbin; i++) {
    if (bin[i] >= bin[i + 1])  {err = TOOLS_BIN_ERROR; goto ErrorHandling; }
  }


  dimM1 =  *dim - 1; 
  halfnbin =  * nbin / 2; 

  switch ( *charact) {
  case 0 :
    characteristic = variogram; 
    factor = 2; 
    break; 
  case 1 : 
    characteristic = Efunction; 
    factor = 2; 
    break; 
  case 2  : 
    characteristic = kmm; 
    factor = 1; 
    break; 
  default : 
    err = TOOLS_UNKNOWN_CHAR; goto ErrorHandling; 
   }
   
  if ((BinSq = (double  * ) MALLOC(sizeof(double) *  ( * nbin + 1)))==NULL) {
    err = TOOLS_MEMORYERROR; goto ErrorHandling; 
  }
  for (i = 0; i< * nbin; i++){
    sd[i] = vario[i] = 0.0; 
    n[i] = 0; 
  }
  for (i = 0; i <=  * nbin; i++){ 
    BinSq[i] = bin[i] > 0 ?  bin[i] * bin[i] : bin[i]; 
  }


  //////////////////////////////////// GRID ////////////////////////////////////
  if ( * grid) {
    
    error("use option 'fft' for data on a grid."); 

    int d1, d2, head, tail, swapstart; 
    double p1[MAXVARIODIM], p2[MAXVARIODIM], distSq, dx[MAXVARIODIM]; 
    int delta[MAXVARIODIM], deltaTwo[MAXVARIODIM], maxtail;    
    long indextail[MAXVARIODIM], indexhead[MAXVARIODIM],
      valuePtail, valuePointhead, 
      segmentbase[MAXVARIODIM], DeltaTwoSegmentbase[MAXVARIODIM], 
      DeltaFourSegmentbase[MAXVARIODIM], 
      SegmentbaseTwo[MAXVARIODIM], SegmentbaseFour[MAXVARIODIM]; 
    double psq[MAXVARIODIM], p[MAXVARIODIM], maxbinsquare; 
    bool forbid;    

    // does not work!! :
    //GetGridSize(x, y, z, dim, &(gridpoints[0]), &(gridpoints[1]), &(gridpoints[2])); 
    // totalpoints = gridpoints[0] * gridpoints[1] * gridpoints[2]; 
    //
    // instead :
    for (dd = 0.0, totalpoints = 1, i = 0; i <= dimM1; i++) {
      gridpoints[i] = (int) (xx[i][XLENGTH]); 
      totalpoints  *= gridpoints[i]; 
      maxdist[i] = (gridpoints[i] - 1) * xx[i][XSTEP]; 
      dd += maxdist[i] * maxdist[i]; 
    }
    EV_METHOD = (int) ( ((EV_TIMES_LESS * BinSq[ * nbin]) < dd) ||
			  (EV_SYSTEM_LARGER<totalpoints)); // ? large system ?
    for (i = 0; i <= dimM1; i++) { dx[i] = xx[i][XSTEP] * 0.99999999; }
    totalpointsrepet = totalpoints *  *repet; 
    segmentbase[0] = 1; // here x runs the fastest;  
    SegmentbaseTwo[0] = segmentbase[0] << 1; 
    SegmentbaseFour[0] = segmentbase[0] << 2; 
    for (i = 1; i <= dimM1; i++) { 
      segmentbase[i] = segmentbase[i - 1] * gridpoints[i - 1]; 
      SegmentbaseTwo[i] = segmentbase[i] << 1; 
      SegmentbaseFour[i] = segmentbase[i] << 2; 
    } 




    EV_METHOD = 0; /// !!! since EV_METHOD = 1 has invalid read !!!!!!
    
    switch (EV_METHOD) {
    case 0 :
      for (i = 0; i <= dimM1; i++) {indexhead[i] = 0; p1[i] = 0.0; }
       
      // loop through all pair of points, except (head, tail) for head=tail, 
      // which is treated separately at the end, if necessary (not relevant for 
      // variogram itself, but for function E and V)
      for (head = 0; head<totalpoints; ) {
	for (i = 0; i <= dimM1; i++) {
	  indextail[i] = 0; 
	  p2[i] = 0.0; 
	}
	for (tail = 0; tail<head; ) {
	  distSq = 0.0; 
	  for (i = 0; i <= dimM1; i++) {
	    double dx2; 
	    dx2 = p1[i] - p2[i]; 
	    distSq  += dx2 * dx2; 
	  }
	  
	  if ((distSq>BinSq[0]) && (distSq <= BinSq[ * nbin])) { 
	    /*  search which bin distSq in  */
	    low = 0; up =   * nbin; /*  21.2.01,  * nbin - 1  */  
	    cur =  halfnbin; 
	    while (low!=up) { 
	      if (distSq> BinSq[cur]) {low = cur;} else {up = cur-1;} // (* ; *]
	      cur = (up + low + 1) / 2; 
	    }
	    for (segment = 0; segment<totalpointsrepet; segment += totalpoints){
	      double x2; 
	      x2 = characteristic(values[head + segment],
				  values[tail + segment]); 
	      vario[low] += x2; 
	      sd[low] += x2 * x2; 
	    }
	    assert(low < *nbin); 
	    n[low]++; 
	  }	
	  tail++; 
	  d2 = dimM1; 
	  indextail[d2]++;
	  p2[d2] +=dx[d2]; 
	  while (indextail[d2] >= gridpoints[d2]) { 
	    indextail[d2] = 0; p2[d2] = 0; 
	    d2--; 
	    assert(d2 >= 0); 	
	    indextail[d2]++; 
	    p2[d2] +=dx[d2]; 
	  }
	}   
	head++; 
	if (head<totalpoints) {
	  d1 = dimM1;  
	  indexhead[d1]++; p1[d1] +=dx[d1];      
	  while (indexhead[d1] >= gridpoints[d1]) { 
	    indexhead[d1] = 0; 
	    p1[d1] = 0; 
	    d1--;  
	    assert(d1 >= 0); 
	    indexhead[d1]++; 
	    p1[d1] +=dx[d1]; 
	  }
	}
      }
      if(PL >= PL_SUBDETAILS) 
	{
	for(i = 0; i< * nbin; i++){
	  PRINTF(" %d:%f(%f, %d)",
		 i, vario[i]/(double)(factor * n[i]), vario[i], n[i]); 
	}
	PRINTF("  XXXY\n"); 
	}      
	break; 

	/* ******************************************************** */

    case 1 :
      /*  idea: 
	 1) fix a vector `delta' of positive components
	 2) search for any pair of points whose distSqance vector dv equals
	    |dv[i]|=delta
	    and some the f(value1, value2) up.
	 
	realisation :
	1) one preferred direction, namely dimM1 - direction (of 0..dimM1)
	2) for the minor directions, swap the values of delta from positive 
	   to negative (all 2^{dimM1} posibilities)   - - > indexhead
	   EXCEPTION :  delta[dimM1]==0. Then only the minor directions up to 
                        dimM1 - 2 are swapped!
                        (otherwise values will be counted twice as indextail is 
                         running through all the positions with
                         coordinate[dimM1]==0)
           if delta[dimM1 - 1] is also 0 then only minor direction up to dimM1 - 3 
           are swapped, etc.
        3) indextail running through all positions with coordinate[dimM1]==0. 
	   This is not very effective if the grid is small in comparison to the
	   right boundary of the bins!    
       */
          
      maxbinsquare = BinSq[ * nbin]; 
      for (i = 0; i <= dimM1; i++) {
	delta[i] = 0; 
	psq[i] = p[i] = 0.0; 
      }
      for (; ; ) { // delta
	do {
	  d = 0;  
	  delta[d]++;  // this excludes distSq==0 !	
	  p[d] +=dx[d]; psq[d] = p[d] * p[d]; 
	  for (distSq = 0.0, i = 0; i<=dimM1; i++) {distSq += psq[i]; }
	  while ((distSq>maxbinsquare)  || (p[d]>maxdist[d])) {/*  ( * )  */ 
	    delta[d] = 0; psq[d] = p[d] = 0.0; 
	    d++;   
	    if (d<=dimM1) { 
	      delta[d]++; 
	      p[d] +=dx[d];
	      psq[d] = p[d] * p[d]; 
	      for (distSq = 0.0, i = 0; i<=dimM1; i++) distSq += psq[i];
	    } else break; 	
	  }
	} while ((distSq<=BinSq[0]) && (d<=dimM1)); 
	if (d>dimM1) break; 
	if (PL >= PL_SUBDETAILS) 
	    PRINTF("\n {%d %d %d}", delta[0], delta[1], delta[2]); 
	assert((distSq<=BinSq[ * nbin]) && (distSq>BinSq[0])); 
	{
	  low = 0; up = *nbin; /*   */ cur =  halfnbin; 
	  while(low!=up){
	    if (distSq> BinSq[cur]) low = cur; else up = cur-1;// ( * ; * ] 
	    cur = (up + low + 1) / 2; 
	  }
	}	

	
	i = 0; 
	while ((i<dimM1) && (delta[i]==0)) i++; 
	swapstart = i + 1; 
	
	for(i = 0; i<=dimM1; i++) {
	  deltaTwo[i]=delta[i] << 1; 
	  DeltaTwoSegmentbase[i]  =  delta[i] * SegmentbaseTwo[i]; 
	  DeltaFourSegmentbase[i] =  delta[i] * SegmentbaseFour[i]; 
	}

	for (i = 0; i<=dimM1; i++) indextail[i] = 0;
	valuePtail = valuePointhead = 0; 
	
	maxtail = gridpoints[0] - delta[0];  
	for(; ; ) { // indextail
	  for (i = 0, valuePointhead = 0; i<=dimM1; i++) {
	    indexhead[i] = indextail[i] + delta[i]; /*  Do not use  - delta[i], as last 
						   dimension must always be 
						    + delta[]  */
	    valuePointhead  += indexhead[i] * segmentbase[i]; 
	  }  // one point is given by indextail, the other by indexhead	
	  for (; ; ) {// indexhead
	    forbid = false; // in case dim=1 ! (???, 21.2.01)
	    for (i = 0; i<=dimM1; i++) { // for (i = 1; ... ) should be enough...
	      if ((forbid = indexhead[i] < 0) || indexhead[i] >= gridpoints[i]) 
		break; 
	    }
	    if (!forbid) { 	  
	      long segtail, seghead; 
	      for (tail = 0; tail<maxtail; tail++) {	
		segtail = valuePtail + tail; 
		seghead = valuePointhead + tail; 
		for (segment=0; segment<totalpointsrepet; segment+=totalpoints){
		  double x2; 
		  x2 = characteristic(values[segtail + segment], 
				      values[seghead + segment]); 
		  vario[low] += x2; 
		  sd[low] += x2 * x2; 
		} 	
		assert(low< *nbin); 
		n[low]++; 
	      }	
	    }
	    
	    d = swapstart;  
	    if (d<=dimM1) {
	      indexhead[d] -=deltaTwo[d];  
	      valuePointhead -= DeltaTwoSegmentbase[d]; 
	      while (indexhead[d] < indextail[d] - delta[d] || delta[d]==0) { 
		if (delta[d] != 0) {
		  // so, indexhead[d]<indextail[d] - delta[d]; 
		  // i.e.indexhead[d] = indextail[d] - 3 * delta[d]
		  valuePointhead += DeltaFourSegmentbase[d]; 
		  indexhead[d] = indextail[d] + delta[d]; 
		}
		d++; 
		if (d<=dimM1) {
		  indexhead[d] -= deltaTwo[d]; 
		  valuePointhead -= DeltaTwoSegmentbase[d]; 
		} else break; 
	      }
	    } else break; 
	    if (d>dimM1) break; 
	  } //for(; ; ), indexhead
	  
	  d = 1;  
	  if (d<=dimM1) {    // not satisfied for one dimensional data! ?
	    indextail[d]++; 
	    valuePtail += segmentbase[d]; 
	    while (indextail[d] >= gridpoints[d]) { 
	      valuePtail -= indextail[d] * segmentbase[d]; 
	      indextail[d] = 0; 
	      d++; 
	      if (d<=dimM1) {
		indextail[d]++; 
		valuePtail += segmentbase[d]; 
	      } else break; 
	    }
	  } else break; 
	  if (d>dimM1) break; 
	} // for(; ; ), indextail
      }  //for(; ; ), delta
      if(PL >= PL_SUBDETAILS) {
	PRINTF("\n"); 
	for(i = 0; i< * nbin; i++){
	  PRINTF(" %d:%f(%f, %d)  ", 
		 i, vario[i]/(double)(factor * n[i]), vario[i], n[i]); 
	}
	PRINTF(" YYY\n"); 
      }
      break; 
    default : BUG;
    }
  } else {
    ////////////////////////////////////  ARBITRARY /////////////////////////////
    totalpoints =  * lx; 
    totalpointsrepet = totalpoints *  * repet; 
    for (i = 0; i<totalpoints; i++){ // to have a better performance for large 
      //                            data sets, group the data first into blocks
      for (j = 0; j<i; j++){
        double distSq; 
      	for (distSq = 0.0, d = 0; d<=dimM1; d++) {
	  double dx; 
	  dx = xx[d][i] - xx[d][j]; 
	  distSq  += dx * dx; 
	}		
	// see also above
	//distSq = sqrt(distSq); 26.2.
	if (distSq>BinSq[0] && distSq<=BinSq[*nbin]) { 
	  low = 0; 
	  up = *nbin; 
	  cur =  halfnbin; 
	  while (low!=up) {  
	    if (distSq> BinSq[cur]) low = cur; else up = cur - 1; // ( * ; * ]  
	    cur = (up + low + 1) / 2; 
	  } 	
	  for (segment = 0; segment<totalpointsrepet; segment  += totalpoints) { 
	    double x2; 
	    x2 = characteristic(values[i + segment], values[j + segment]); 
	    if (R_FINITE(x2)) {
	      vario[low] += x2; 
	      sd[low] += x2 * x2; 
	      n[low]++; 
	    }
	  }	
	  assert(low< * nbin); 
	}
	/*   */
      }
    }
  }
  
  // ii is used in switch by this value: !!!!
  ii = 0; while(ii <= *nbin && bin[ii] < 0) ii++;
  ii--; 

  // 11.10.03: factor ohne  * repet to be clearer
  // case ii : separat behandelt, neu geschrieben
  //  if ( * repet>1) {
  //   for (j=0; j< * nbin; j++) n[j]  *=  * repet; 
  //}

  for (i = 0; i<ii; i++) sd[i] = vario[i] =  RF_NA;
  for (i++; i< * nbin; i++){
    if (n[i]>0) { 
      vario[i] /= (double) (factor * n[i]); 
      if (n[i]>1) 
	sd[i] = sqrt(sd[i] / (factor * factor * (double) (n[i] - 1))
		     - (double) (n[i]) / (double) (n[i] - 1) * 
		     vario[i] * vario[i]); 
      else sd[i] = RF_INF; 
    } else { 
      vario[i] =  RF_NA; 
      sd[i] = RF_INF; 
    }
  }

  i = ii; // sicherheitshalber - -  eigentlich sollte unten kein i mehr auftauchen
  long jj; 
  switch ( * charact) {
  case 0 :
    if ((ii< * nbin)&&(ii >= 0)){
      n[ii]  += totalpointsrepet; 
      vario[ii] /=  (double) (factor * n[ii]); 
      sd[ii] = sqrt(sd[ii] / (factor * factor * (double) (n[ii] - 1))
		    - (double) (n[ii]) / (double) (n[ii] - 1) * vario[ii] * vario[ii]); 
    }    
    break; 
  case 1 : // E, calculating E(0) correctly, taking into account that the
    //        bin including 0 may include further points
    if ((ii< * nbin)&&(ii >= 0)){
      n[ii]  *= factor; // * 2
      for (jj = 0; jj<totalpointsrepet; jj++) {
	vario[ii]  += values[jj]; 
      }
      n[ii]  += totalpointsrepet; 
      vario[ii] = vario[ii] / (double) n[ii]; 
      for (i = 0; i< * nbin; i++) sd[i] = RF_NA; // beste Loesung, da sonst
      // nur Chaos gibt - -  da sd falsch berechnet ist.
      // siehe MPP package fuer die richtige Berechnung.
    }
    break; 
  case 2 : // kmm, calculating kmm(0} and dividing all the values by mean^2
    double mean, square, qu, x2; 
    for (jj = 0, mean = square = qu = 0.0; jj<totalpointsrepet; jj++) { 
      mean  += values[jj]; 
      square  += (x2 = values[jj] * values[jj]); 
      qu  += x2 * x2; 
    }
    mean /= (double) totalpointsrepet; 
    square /= (double) totalpointsrepet; 
    if ((ii< * nbin)&&(ii >= 0)){
      n[ii]  += totalpointsrepet;    
      vario[ii] =  (vario[ii] + square) / (double) n[ii]; // factor is 1
      sd[ii] = sqrt(((sd[ii] + qu) - (double) (n[ii]) * vario[ii] * vario[ii])
		    / ((double) (n[ii] - 1))); 
    }
    mean *= mean; 
    for (jj = 0; jj< * nbin; jj++) { vario[jj] /= mean; }
    break; 
  default : BUG;    
  }

  if (BinSq!=NULL) free(BinSq); 
  return; 
  
 ErrorHandling:
  if (BinSq!=NULL) free(BinSq); 
  for (i = 0; i< * nbin; i++){sd[i] = vario[i] = RF_NA; } 
  switch (err) {
  case TOOLS_DIM :
    error("dimension exceed max dimension of empirical variogram estimation"); 
  case TOOLS_MEMORYERROR :  
    error("Memory alloc failed in empiricalvariogram.\n"); 
  case TOOLS_XERROR :  
    error("The x coordinate may not be NULL.\n"); 
  case TOOLS_BIN_ERROR :
    error("Bin components not an increasing sequence.\n"); 
  case TOOLS_UNKNOWN_CHAR:
    error("unknown type of second order characteristic"); 
  default : BUG;
  }
}
