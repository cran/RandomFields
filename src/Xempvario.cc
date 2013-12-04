/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de 

 a second library for calculating the empirical variogram

 Copyright (C) 2002 - 2013 Martin Schlather, 

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
#include <stdlib.h>
#include "RF.h"
 

// z coordinate run the fastest in values, x the slowest   

//#define debug_tools 1
#define TOOLS_MEMORYERROR 1
#define TOOLS_XERROR 2
#define TOOLS_BIN_ERROR 3
#define NEARBY 1e15

// naechste Zeile nur notwendig, weil atan2 in Windows nicht
// ordentlich programmiert ist
#define NEARBYINT(x)  floor((x) * NEARBY + 0.5) / NEARBY
// #define NEARBYINT(x)  x

void empvarioXT(double *X, double *T, 
		int *lx, 
		double *values, int *repet, int *grid,
		double *bin, int *nbin, 
		double *phi,    // vector of a real and an integer
		double *theta,  // vector of a real and an integer
		int *dT,   // in grid units, smallest positive value, for  
		//            the temporal lag of the emp. variogram
	        //            stepT * nstepT is the greatest time span
		double *sum,   // \sum (a-b)^2 / 2
		double *sq,   // \sum (a-b)^4 / 4
		int *n)
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
 */
{ 
  long rep;
  int i, halfnbin,gridpoints[4], err, totalbins, totalspatialbins,
    twoNphi, twoNphiNbin, Ntheta, nbinNphiNtheta, maxi[4], mini[4], ix, iy, iz,
    it, endX, startX, endY, startY, endZ, startZ, endT, low, cur, up,
    ktheta, kphi, nstepTstepT, vec[4], deltaT, x, y, z, t,  j, stepT, nstepT;
  double *xx[4], *BinSq, maxbinsquare, step[4], delta[4], psq[4], dx[4],
    startphi, invsegphi, starttheta, invsegtheta, thetadata, phidata,
    zylinderradius;

  stepT = dT[0];
  nstepT = dT[1];
  nstepTstepT = stepT * nstepT;
  twoNphi =  phi[1]==0 ?  1 : (2 * (int) phi[1]);
  twoNphiNbin = twoNphi * *nbin;
  Ntheta = theta[1]==0 ? 1 : (int) theta[1];
  startphi = phi[0] - PI / (double) twoNphi; // [0, 2 pi]
  invsegphi = phi[1] / PI; // note that phi[1] can be zero!
  //starttheta = theta[0] - PIHALF / (double) Ntheta;  // [0, pi]
  starttheta = theta[0] - PIHALF;  // [0, pi]
  invsegtheta = theta[1] / PI; // note that theta[1] can be zero
  nbinNphiNtheta = twoNphiNbin * Ntheta;
 
  BinSq = NULL;
  for (rep=i=0; i<3; i++, rep+=*lx) xx[i]=&(X[rep]);
  xx[3] = T;

  if (xx[0]==NULL) {err=TOOLS_XERROR; goto ErrorHandling;}
  for (i=0; i<*nbin;i++) {
    if (bin[i]>=bin[i+1])  {err=TOOLS_BIN_ERROR; goto ErrorHandling;}
  }

  halfnbin = *nbin / 2;
   
  if ((BinSq = (double *) MALLOC(sizeof(double)* (*nbin + 1)))==NULL) {
    err=TOOLS_MEMORYERROR; goto ErrorHandling; 
  }
  totalspatialbins =  twoNphiNbin * Ntheta;
  totalbins = totalspatialbins * (nstepT + 1);
  for (i=0; i<totalbins; i++){sq[i]=sum[i]=0.0; n[i]=0;}
  for (i=0; i<=*nbin; i++){if (bin[i]>0) BinSq[i]=bin[i] * bin[i]; 
  else BinSq[i]=bin[i];
  }

  assert(NEARBYINT(atan2(-1.0, 0.0) + PIHALF) == 0.0);
  assert(atan2(0.0, 0.0) == 0.0);
  maxbinsquare = BinSq[*nbin];


  //////////////////////////////////// GRID ////////////////////////////////////
  if (*grid) {
    error("use option 'fft' for space-time data on a grid");

    int segmentbase[6];

    segmentbase[0]=1; // here x runs the fastest; 
    for (i=0; i<=3; i++) {
      step[i] = xx[i][XSTEP];

      gridpoints[i] = (int) xx[i][XLENGTH]; 
      maxi[i] = (int) (sqrt(maxbinsquare) / step[i] + 0.1);
      if (maxi[i] >= gridpoints[i]) maxi[i] = gridpoints[i] - 1;
      mini[i] = -maxi[i];
      //      if (maxi[i]==0)  maxi[i] = 1; wrong 22.11.03
      segmentbase[i+1] =  segmentbase[i] * gridpoints[i];
    }
    

    segmentbase[5] = segmentbase[4] * *repet;
    // sementbase: [0] : 1, [1] : length(x), [2] : length(x) * l(y),
    //  [3] : l(x) * l(y) * l(z), [4] : l(x) * l(y) * l(z) * l(T)
    //  [5] : l(x) * l(y) * l(z) * l(T) * repet


    for (ix=mini[0]; ix <= maxi[0]; ix++) {
      delta[0] = ix * step[0];
      if ((psq[0] = delta[0] * delta[0]) > maxbinsquare) continue;

      vec[0] = ix;
      if (ix>=0) {
	endX = gridpoints[0] - ix;
	startX = 0;
      } else {
	endX = gridpoints[0]; 
	startX = -ix;
      }
      for (iy= mini[1]; iy <= maxi[1];  iy++) {
	delta[1] = iy * step[1];
	if ((psq[1] = delta[1] * delta[1] + psq[0]) > maxbinsquare) continue; 
	vec[1] = vec[0] + iy * segmentbase[1];
	
        // angles
	phidata = NEARBYINT(atan2(delta[1], delta[0]) - startphi);
	while (phidata < 0) phidata += TWOPI;
	while (phidata >= TWOPI) phidata -= TWOPI;
	kphi = *nbin * (int) (phidata * invsegphi);


	for (iz=mini[2]; iz <= maxi[2]; iz++) {
	  delta[2] = iz * step[2];
	  if ((psq[2] = delta[2] * delta[2] + psq[1]) > maxbinsquare) continue;
	  vec[2] = vec[1] + iz * segmentbase[2];
	 
	  {
	    low=0; up= *nbin; /* */ cur= halfnbin;
	    while(low!=up){
	      if (psq[2] > BinSq[cur]) {low=cur;} else {up=cur-1;}/*( . ; . ]*/ 
	      cur=(up+low+1)/2;
	    }
	  }	// low

	  // angles
	  thetadata = NEARBYINT(atan2(delta[2], sqrt(psq[1])) - starttheta);
	  while (thetadata < 0) thetadata += PI; 
	  while (thetadata >= PI) thetadata -= PI;
	  ktheta = (int) (thetadata * invsegtheta);
	  
	  low += kphi + ktheta * twoNphiNbin;

	  for (deltaT=0, it=low; deltaT<=nstepTstepT; 
	       deltaT+=stepT, it+=nbinNphiNtheta) {
	    vec[3] = vec[2] + deltaT * segmentbase[3];
	    
	    assert(startX>=0);
	    for (x=startX; x<endX; x++) {
	      if (iy>=0) {
		endY = gridpoints[1]  * segmentbase[1] - vec[1];
		startY = x;
	      } else {
		endY = gridpoints[1] * segmentbase[1]; 
		startY = x - iy * segmentbase[1];
	      }    
	      assert(startY>=0);  
	      for (y=startY; y<endY; y+=segmentbase[1]) {
		if (iz>=0) {
		  endZ = gridpoints[2] * segmentbase[2] - vec[2];
		  startZ = y;
		} else {
		  endZ = gridpoints[2] * segmentbase[2];
		  startZ = y - iz * segmentbase[2];
		}      
                assert(startZ>=0);
		for (z=startZ; z<endZ; z+=segmentbase[2]) {		  
		  endT = gridpoints[3] * segmentbase[3] - vec[3];
		  for (t=z; t<endT; t+=segmentbase[3]) {
		    for (rep = t; rep < segmentbase[5]; rep += segmentbase[4]) {
		      double x2;
		      int rv;
		      rv = rep + vec[3];
		      assert(rv < segmentbase[5] && rv >=0);
		      x2 = values[rep] - values[rv];
		      x2 *= x2;
		      sum[it] += x2;
		      sq[it] += x2 * x2;
		    } // repeat	
		    n[it]++;
		  } // t
		} // z
	      } // y
	    } // x		 
	  } // deltaT
	} // iz
      } // iy
    } // ix
  } else {
    ////////////////////////////////////  ARBITRARY /////////////////////////////
    // rows : x, columns : T 
    long totalpoints, totalpointsrepet, spatial,  jj;
    spatial = *lx;
    i = 3;
    step[i] = xx[i][XSTEP];
    gridpoints[i] = (int) xx[i][XLENGTH];
    totalpoints = spatial * gridpoints[i];
    totalpointsrepet = totalpoints * *repet;
    for (i=0;i<spatial;i++) { // to have a better performance for large 
      //                   data sets, group the data first into blocks
      for (j=0; j<spatial; j++) {
	double distSq;
	dx[0] = xx[0][j] - xx[0][i];
	dx[1] = xx[1][j] - xx[1][i];
	dx[2] = xx[2][j] - xx[2][i];
	distSq = dx[0] * dx[0] + dx[1] * dx[1];
	zylinderradius = sqrt(distSq);
	distSq += dx[2] * dx[2];
			
	if ((distSq>BinSq[0]) && (distSq<=BinSq[*nbin])) {
	  low=0; up=*nbin; cur=halfnbin;
	  while (low!=up) {
	    if (distSq> BinSq[cur]) {low=cur;} else {up=cur-1;} // ( * ; * ]  
	    cur=(up+low+1)/2; 
	  }

	  // angles
	  phidata = NEARBYINT(atan2(dx[1], dx[0]) - startphi);
	  while (phidata < 0) phidata += TWOPI;
	  while (phidata >= TWOPI) phidata -= TWOPI;
	  kphi = *nbin * (int) (phidata * invsegphi);

	  thetadata = NEARBYINT(atan2(dx[2], zylinderradius) - starttheta);
	  while (thetadata < 0) thetadata += PI; 
	  while (thetadata >= PI) thetadata -= PI;
	  ktheta = (int) (thetadata * invsegtheta);

	  low += kphi + twoNphiNbin * ktheta;
	  assert(low>=0 && low<totalspatialbins);
	  
	  for (deltaT=0; deltaT<= nstepTstepT; deltaT+=stepT,
		 low+=totalspatialbins) {
	    jj = j + deltaT * spatial;
	    endT= totalpoints - deltaT * spatial;
	    for (t=0; t<endT; t+=spatial) {
	      for (rep = t; rep < totalpointsrepet; rep += totalpoints){ 
		double x2;
		x2 = values[rep + jj] - values[i + rep];
		if (R_FINITE(x2)) {
		  x2 *= x2;
		  sum[low] += x2;
		  sq[low] += x2 * x2;
		  n[low]++;
		}
	      }
	    } // t
	  } // deltat
	} // if (distSq)
      } // j
    } // i
  } // arbitrary
    

  for (i=0; i<totalbins; i++) {
    sum[i] *= 0.5;
    sq[i] *= 0.25;
  }
  if (BinSq!=NULL) free(BinSq);

  return;

 ErrorHandling:
  PRINTF("Error: ");
  switch (err) {
  case TOOLS_MEMORYERROR :  
    PRINTF("Memory alloc failed in empiricalvariogram.\n"); break;
  case TOOLS_XERROR :  
    PRINTF("The x coordinate may not be NULL.\n"); break;
  case TOOLS_BIN_ERROR :
    PRINTF("Bin components not an increasing sequence.\n"); break;
  default : assert(false);
  }
  if (BinSq!=NULL) free(BinSq);
  for (i=0;i<*nbin;i++){sq[i]=sum[i]=RF_NAN;} 
} 
	 



