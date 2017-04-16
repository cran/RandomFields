/* 
 Authors
 Sebastian Engelke
 Johannes Martini

 library for simulation of random fields -- init part and error messages

 Copyright (C) 2011 -- 2013 Sebastian Engelke, Johannes Martini
 Copyright (C) 2014 Sebastian Engelke, Johannes Martini, Martin Schlather
 Copyright (C) 2015 -- 2017 Martin Schlather

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

#include <Rmath.h>
#include <stdio.h>
//#include <stdlib.h>
#include "RF.h"
 
// z coordinate run the fastest in values, x the slowest   

#define TOOLS_MEMORYERROR 501
#define TOOLS_XERROR 502
#define TOOLS_BIN_ERROR 503
#define NEARBY 1e15

// naechste Zeile nur notwendig, weil atan2 in Windows nicht
// ordentlich programmiert ist
#define NEARBYINT(x)  (FLOOR((x) * (NEARBY) + 0.5) / (NEARBY))
// #define NEARBYINT(x)  x



int GetAngleBin(double angle, double startangle, double SegPerPI, double maxAngle){
  double phi;  
  int kphi;
  phi = angle - startangle;
  while (phi < 0) phi += maxAngle;
  while (phi >= maxAngle) phi -= maxAngle;
  kphi = (int) (phi * SegPerPI / PI);
  return(kphi);
}

static int maxjj = 0;
double get(double *sumvals, int jj) {
  if (jj > maxjj) {
    maxjj = jj;
    //printf("jj=%d\n", jj);
  }
  return sumvals[jj];
}


SEXP fftVario3D(SEXP Coord,
		SEXP Sumvals,
		SEXP Nbvals,
		SEXP Bin,
		SEXP Nbin,
		SEXP LenT,
		SEXP StepT,
		SEXP NstepT,
		SEXP Phi,
		SEXP Theta,
		SEXP Repet,
		SEXP Vdim,
		SEXP SegmentEmpVario,
		SEXP Pseudo
		) {

/*    coord   : 3x3 matrix of coordinates of the 3 space dimension
 *              in (new) gridtriple notation, (start,by,length)
 *    sumvals : matrix of non-averaged values of the empirical variogram
 *              as a function of the distance vector
 *    nbvals  : number of values contributing to the respective direction
 *              sumvals / nbvals = meanvals   
 *    bin     : vector of bin boundaries for the radial space part  
 *    nbin    : number of bins of the radial space part
 *    stepT   : stepsize in temporal direction (must be multiple of the time stepsize)        
 *    nstepT  : number of steps in temporal direction
 *    phi     : 2-dim vector, first component gives the starting angle and 
 *              second component the number of angles / PI (all in the xy-plane)
 *    theta   : 2-dim vector, first component gives the starting angle and 
 *    repet   : number of repetitions (the calculated emprical variogram will 
 *              be an overall average, see also EmpiricalVariogram in empvario.R)
 *              second component the number of angles / PI (all orthogonal to xy-plane)
 *    empvario: vector of bins filled with averaged empirical variogram
 *    n       : number of pairs in each bin
 *    pseudo  : if pseudo == 1, then the pseudo (cross)-variogram is computed
 */

  long rep;
  int i, binidx, halfnbin, totalbinsvdimSq,
    Nphi, NphiNbin, Ntheta, nbinNphiNtheta, Nphihalf, maxi[3],
    ix, iy, iz, startX, startY, startZ, startT, startrep, low, cur,
    ktheta, kthetaz, kphi, kphixy, kphix, kphiy, kT, nstepTstepT, deltaT, 
    x, y, z, v1, v2,
    timecomp, k,
    err = NOERROR, 
    nbin = INTEGER(Nbin)[0],
    vdim = INTEGER(Vdim)[0],
    repet = INTEGER(Repet)[0],
    segmentEmpVario = INTEGER(SegmentEmpVario)[0],
    totalbins = segmentEmpVario * vdim * vdim,
    nstepT = INTEGER(NstepT)[0],
    stepT = INTEGER(StepT)[0],
    lenT = INTEGER(LenT)[0]
   ;
  bool pseudo = LOGICAL(Pseudo)[0];
  double tolerance, *xx[3], maxbinsquare, step[3], 
    d0, d1, d2,  // delta0
    psq0, psq1, psq2,
    startphi, starttheta, thetadata, phidata, phixdata,
    phiydata, phixydata, binshift,
    *phi = REAL(Phi),
    *theta = REAL(Theta),
    *BinSq = NULL,  
    *sumvals = REAL(Sumvals),
    *nbvals=  REAL(Nbvals),
    *bin = REAL(Bin),
    *coord = REAL(Coord);
  
  SEXP back;
  back = PROTECT(allocMatrix(REALSXP, totalbins, 2) );
  
  double *emp_vario = REAL(back),
    *n = emp_vario + totalbins
    ;

#define SUMVALS_JJ sumvals[jj]
  //#define SUMVALS_JJ get(sumvals, jj)

  
  kphix = 0;
  kphixy = 0;
  kphiy = 0;
  Nphihalf = 0;

  nstepTstepT = stepT * nstepT;
  timecomp = nstepTstepT==0 ? 1 : 2;

  if((!pseudo) && (timecomp==1)){
    Nphi =  phi[1]==0 ?  1 : (int) phi[1];    
  }
  else{
    Nphihalf =  phi[1]==0 ?  0 : (int) phi[1];
    Nphi =  phi[1]==0 ?  1 : 2 * ((int) phi[1]);
  }
  Ntheta = theta[1]==0 ? 1 : (int) theta[1];
  //print("1,..");  
  NphiNbin = Nphi * nbin;
  startphi = ((!pseudo) && (timecomp==1)) ? (phi[0] - PI / (double) (2*Nphi))
    : (phi[0] - PI / (double) (Nphi) ); // [0, 2 pi]    
  //invsegphi = phi[1] / PI; // note that phi[1] can be zero!
  //starttheta = theta[0] - PI / (double) (2*theta[1]);  // [0, pi]
  starttheta = 0;
  // invsegtheta = theta[1] / PI; // note that theta[1] can be zero
  nbinNphiNtheta = NphiNbin * Ntheta;
  //print("2,..");
  for (rep=i=0; i<3; i++, rep+=3) xx[i] = coord + rep;
    if (xx[0]==NULL) {err=TOOLS_XERROR; goto ErrorHandling;}
  for (i=0; i<nbin;i++) {
    if (bin[i]>=bin[i+1])  {err=TOOLS_BIN_ERROR; goto ErrorHandling;}
  }

  halfnbin = nbin / 2;
  
  if ((BinSq = (double *) MALLOC(sizeof(double)* (nbin + 1))) ==NULL) {
    err=TOOLS_MEMORYERROR; goto ErrorHandling; 
  }
  //print("3,..");
  totalbins = nbinNphiNtheta * (nstepT + 1);
  totalbinsvdimSq = totalbins * vdim * vdim;
  for (i=0; i<totalbinsvdimSq; i++){
    //printf("%d ", i);
    emp_vario[i]= n[i]=0.0;
  }

  //print("sizeof double: %d ", (int) sizeof(double));
  if (sizeof(double) == 8 )
    binshift = bin[nbin] * bin[nbin] * (sizeof(double) == 8 ? 1e-12 : 1e-4);
  // shift all bins to the right to overcome rounding errors
  for (i=0; i<=nbin; i++) {
    BinSq[i] = bin[i]>0 ? bin[i] * bin[i] + binshift : bin[i];
  }
  
  assert(NEARBYINT(atan2(-1.0, 0.0) + PIHALF) == 0.0);
  assert(atan2(0.0, 0.0) == 0.0);
  maxbinsquare = BinSq[nbin];
  
  int segmentbase[10];
  
  segmentbase[0]=1; // here x runs the fastest; 
  for (i=0; i<=2; i++) {
    step[i] = xx[i][1];    
    maxi[i] = xx[i][2];
    segmentbase[i+1] =  segmentbase[i] * maxi[i];
  }
  
  segmentbase[4] = segmentbase[3] * lenT; 
  if (!pseudo) {
    segmentbase[5] = 4 * segmentbase[4];
    segmentbase[6] = segmentbase[5] * timecomp;
  } else {
    segmentbase[5] = 8 * segmentbase[4];
    segmentbase[6] = segmentbase[5];
  }    
  segmentbase[7] = segmentbase[6] * repet;
  segmentbase[8] = segmentbase[7] * vdim * vdim;
  // sementbase: [0] : 1, [1] : length(x), [2] : length(x) * l(y),
  //  [3] : l(x) * l(y) * l(z), [4] : l(x) * l(y) * l(z) * l(T)
  //  [5] : 4 (or 8) * l(x) * l(y) * l(z) * l(T) 
  //  [6] : 4 (or 8) * l(x) * l(y) * l(z) * l(T) * timecomp
  //  [7] : 4 (or 8) * l(x) * l(y) * l(z) * l(T) * timecomp * repet    
  //  for pseudo variogram we need the whole sphere (not only half-sphere)
  

  for (v1=0; v1<vdim; v1++) {
    for (v2=0; v2<vdim; v2++,
	   sumvals += segmentbase[7], 
	   nbvals += segmentbase[7],
	   emp_vario += segmentEmpVario,
	   n += segmentEmpVario) {
  
      //
      //for (ix=0; ix<=8; ix++) printf("%d:%d ", ix, segmentbase[ix]);
      //printf(">> ... %d %d %d\n", v1, v2, segmentEmpVario);
      startX = 0;  
      binidx = 0;  
      for (x=startX, ix=0; x < segmentbase[1]; x++, ix++) {
	
	d0 = ix * step[0];   
	startY = x;
	if ((psq0 = d0 * d0) > maxbinsquare) continue;      
	
	for (y = startY, iy=0; y < segmentbase[2];  y+=segmentbase[1], iy++){
	  
	  d1 = iy * step[1];
	  startZ = y;
	  if ((psq1 = d1 * d1 + psq0) >maxbinsquare) continue; 
	  
	  // angles
	  if(!pseudo){
	    phidata = NEARBYINT(atan2(d1, d0) - startphi);
	    kphi = GetAngleBin(phidata, 0, phi[1], PI);
	    phixdata = NEARBYINT(atan2(d1, -d0) - startphi);
	    kphix = GetAngleBin(phixdata, 0, phi[1], PI);
	    if((d0 == 0) && (d1 == 0)) kphi = 0;// points on the z-axes 
	  }
	  else{
	    phidata = NEARBYINT(atan2(d1, d0) - startphi);
	    kphi = GetAngleBin(phidata,0, phi[1], TWOPI);
	    phixdata = NEARBYINT(atan2(d1, -d0) - startphi);
	    kphix = GetAngleBin(phixdata, 0, phi[1], TWOPI);
	    phiydata = NEARBYINT(atan2(-d1, d0) - startphi);
	    kphiy = GetAngleBin(phiydata, 0, phi[1], TWOPI);
	    phixydata = NEARBYINT(atan2(-d1, -d0) - startphi);
	    kphixy = GetAngleBin(phixydata, 0, phi[1], TWOPI);
	  }
	  
	  //if(d1 == d0)
	  //print("kphi: %d and kphix: %d, phidata: %lf, startphi: %lf", kphi, kphix, phidata, startphi);
	  // kphi is index of angle bin in x,y > 0 plane
	  // kphix is index of angle bin if phi 
	  // is mirrored along the x-axis (in the xy-plane)
	  // kphiy and kphixy similar
	  
	  for (d2=0.0, z=startZ, iz=0; 
	       z < segmentbase[3]; 
	       z+=segmentbase[2], iz++, d2 += step[2]){
	    startT = z;
	    if ((psq2 = d2 * d2 + psq1) > maxbinsquare) continue;
	    
	    //print("5,..");
	    { /* finds the bins for the radial space part */
	      int up = nbin;
	      low=0;
	      cur= halfnbin;
	      while(low!=up){
		if (psq2 > BinSq[cur]) {low=cur;} else {up=cur-1;}/*( . ; . ]*/ 
		cur=(up+low+1) / 2;
	      }
	    }	// low
	    
	    // angles
	    thetadata = NEARBYINT(PIHALF - atan2(d2, SQRT(psq1)));
	    ktheta = GetAngleBin(thetadata, starttheta, theta[1], PI);
	    kthetaz = GetAngleBin(PI - thetadata, starttheta, theta[1], PI);
	    
	    for (deltaT=0, kT=0; deltaT<=nstepTstepT; 
		 deltaT += stepT, kT += nbinNphiNtheta) {
	      //print("kT %d, stepT %d, nsteptT %d, deltaT %d", kT, *stepT, *nstepT, deltaT); 	 
	      //print("kT %d --", kT);
	      //timecomp = x==0 ? 1 : 2;
	      
	      startrep = startT + deltaT * segmentbase[3];
	      for (k=1; k <= timecomp; k++, startrep += segmentbase[5]){ 
		//if a time component is given, then, for statistical efficiency, the data is reflected also along the time component
		binidx = cur + kT;
		assert( binidx < totalbins && binidx >= 0 );
		//print("ix = %d -- k = %d", ix, k);
		//for (i=0, T=(startT+(deltaT+i)*segmentbase[3]); i < endT; T+=segmentbase[3], i++) {
		
		for (rep = startrep; rep < segmentbase[7]; rep += segmentbase[6]) {
		  //print("6,..");
		  //print("ix: %d, iy: %d, iz: %d, kphi: %d, ktheta: %d, kphix: %d, kthetaz: %d, \n", ix, iy, iz, kphi, ktheta, kphix, kthetaz);
		  //if(binidx == 0) print("rep: %d , sumvals[rep]: %lf \n", rep, sumvals[rep]);	      
		  if(!pseudo) {
		    int basX = NA_INTEGER, 
		      jj = rep,
		      kM1 = k - 1,	
		      TkM2 = 2 - k,
		      base = binidx + (kphi + Nphihalf * kM1) * nbin,
		      Idx = ((Ntheta - 1 - ktheta) * kM1 + 
			     ktheta * TkM2) * NphiNbin;

		    //printf("rep =%d [%d %d; %d]  v1=%d %d [%d] x=%d %d %d deltaT=%d %d\n", 
		    //rep, startrep, segmentbase[7], segmentbase[6],
		    //v1, v2, vdim, x, y, z, deltaT, k);
		    
		    emp_vario[base + Idx] += SUMVALS_JJ;
		    n[base + Idx] += nbvals[jj];	
		    jj += segmentbase[4];
		    if (ix > 0 && (iy > 0 || iz > 0)) {
		      basX = binidx + (kphix + Nphihalf * kM1) * nbin;
		      emp_vario[basX + Idx] += SUMVALS_JJ;
		      n[basX + Idx] += nbvals[jj];
		    }
		    if (iy > 0 && iz > 0){
		      Idx = ((Ntheta -1 - kthetaz) * kM1 + kthetaz * TkM2) * 
			NphiNbin;
		      jj += segmentbase[4];
		      emp_vario[base + Idx] += SUMVALS_JJ;
		      n[base + Idx] += nbvals[jj];
		      
		      if (ix > 0){
			jj += segmentbase[4];
			emp_vario[basX + Idx] += SUMVALS_JJ;
			n[basX + Idx] += nbvals[jj];
		      }
		    }
		  } else { // not pseudo
		    int 
		      jj = rep,
		      knbin = binidx + kphi * nbin,
		      kxnbin = binidx + kphix * nbin,
		      kN = ktheta * NphiNbin,
		      kzN = kthetaz * NphiNbin,
		      Idx = knbin + kN;
		    emp_vario[Idx] += SUMVALS_JJ;
		    n[Idx] += nbvals[jj];
		    jj += segmentbase[4];
		    if (ix > 0) {
		      Idx = kxnbin + kN;
		      emp_vario[Idx] += SUMVALS_JJ;
		      n[Idx] += nbvals[jj];
		    }
		    if (iz > 0) {
		      Idx = knbin + kzN;
		      jj += segmentbase[4];
		      emp_vario[Idx] += SUMVALS_JJ;
		      n[Idx] += nbvals[jj]; 
		      
		      if (ix > 0) {
			Idx = kxnbin + kzN;
			jj += segmentbase[4];
			emp_vario[Idx] += SUMVALS_JJ;
			n[Idx] += nbvals[jj];
		      }
		    }
		    if(iy > 0) {
		      int
			kxynbin = NA_INTEGER,
			kynbin = binidx + kphiy * nbin;
		      Idx = kynbin + kN;
		      jj = rep + 4 * segmentbase[4];
		      emp_vario[Idx] += SUMVALS_JJ;
		      n[Idx] += nbvals[jj];	
		      jj += segmentbase[4];
		      if (ix > 0){ // vor ix>0 && iz > 0 !
			kxynbin = binidx +  kphixy * nbin;
			Idx = kxynbin + kN;
			emp_vario[Idx] += SUMVALS_JJ;
			n[Idx] += nbvals[jj];
		      }
		      if (iz > 0){
			Idx = kynbin + kzN;
			jj += segmentbase[4];
			emp_vario[Idx] += SUMVALS_JJ;
			n[Idx] += nbvals[jj];  
			
			if (ix > 0) {
			  Idx = kxynbin + kzN;
			  jj += segmentbase[4];
			  emp_vario[Idx] += SUMVALS_JJ;
			  n[Idx] += nbvals[jj];
			}
		      } // iz > 0
		    } // iy > 0
		  } // else (not pseudo)	  
		} // repet
	      } // k timecomponent
	    } // deltaT ; time bin
	  } // z
	} // y
      } //x	
  
      tolerance = GLOBAL.empvario.tol * segmentbase[3];  // to do. Warum?
      for(i=0; i < totalbins ;i++){
	if (FABS(emp_vario[i]) < tolerance) emp_vario[i] = 0.0; 
      }
    } // vdim1
  } // vdim
  
    
  // print("C SCRIPT END !!!!!!!!! \n");
  
 ErrorHandling:
  if (err != NOERROR) {
    switch (err) {
    case TOOLS_MEMORYERROR :  
      ERR("Memory alloc failed in empiricalvariogram.\n");
    case TOOLS_XERROR :  
      ERR("The x coordinate may not be NULL.\n"); break;
    case TOOLS_BIN_ERROR :
      ERR("Bin components not an increasing sequence.\n"); break;
    default : BUG;
    }
  }

  FREE(BinSq);
  UNPROTECT(1);
  return back;
}




/* coord   : 3xd matrix of coordinates of d dimensional space in gridtriple notation, 
 * last row is number of gridpoints
 * S : array of non-averaged values of the empirical variogram
 * as a function of the distance vector
 * N  : number of values contributing to the respective direction
 * sumvals / nbvals = meanvals   
 * bounds     : vector of bin boundaries for the radial space part  
 * nBounds    : length of the bounds vector
 * stepT   : stepsize in temporal direction (must be multiple of the time stepsize)        
 * nstepT  : number of steps in temporal direction
 * repeat   : number of repetitions
 */
// void fftVarioIsotrop(
// 	double*	spaceStepLength,
// 	int*		nSpaceSteps,
// 	int*		dim,
// 	int*		nTimeSteps,
// 	double*	S,
// 	int*		W,
// 	double*	spaceBounds,
// 	int*		nSpaceBounds,
// 	int*		timeBounds,
// 	int*		nTimeBounds,
// 	int*		repeat,
// 	double*	matrix
// )
// {
// 	/*
// 	 * ******************************************************************************************** 
// 	 * DEKLARATION
// 	 * ******************************************************************************************** 
// 	 */
	 
// 	// die bin-Matrix fuer die aufsummierten werte
// 	double sumBins[*nSpaceBounds][*nTimeBounds];
// 	// die bin-Matrix fuer die gewichte
// 	int weightBins[*nSpaceBounds][*nTimeBounds];
	
// 	// der cache fuer die quadrierten werte der bin grenzen
// 	double sqrBoundCache[*nSpaceBounds];

// 	// der cache fuer die quadrierten werte der koordinatenpunkte
// 	double sqrCache[*dim];

// 	// der koordinatenvektor
// 	int coord[*dim];
// 	// der dazugehoerige index
// 	int index;
// 	// die dazugehoerige norm quadriert (so spart man sich das staendige wurzelziehen, denn man will ja nur distanzen vergleichen)
// 	double sqrNorm;
	
// 	// der zeit index
// 	int time;

// 	// der wiederholungsindex
// 	int rep;
	
// 	// die bin nummer im raum
// 	int spaceBin;
// 	// die bin nummer in der zeit
// 	int timeBin;
	
// 	// die anzahl der koordinatenpunkte
// 	int nPoints;
// 	// die anzahl der zeit + koordinatenpunkte
// 	int nSpaceTimePoints;
	
// 	// die dimension, in die gerade gegangen wird
// 	int current;
		
// 	// flag, das anzeigt ob in der "wandernden" dimension ein wraparound stattfand
// 	int wrap;
	
// 	/*
// 	 * ******************************************************************************************** 
// 	 * INITIATION
// 	 * ******************************************************************************************** 
// 	 */
	
// 	for (int i= 0; i < *nSpaceBounds; i++)
// 	{
// 		sqrBoundCache[i]= spaceBounds[i] * spaceBounds[i];
		
// 		for (j= 0; j < *nTimeBounds; i++)
// 		{
// 			sumBins[i][j]= 0;
// 			weightBins[i][j]= 0;
// 		}
// 	}

// 	nPoints= 1;	
// 	for (int i= 0; i < *dim; i++)
// 	{
// 		nPoints*= nSpaceSteps[i];
// 		sqrCache[i]= 0;
// 		coord[i]= 0;
// 	}
// 	nSpaceTimePoints= nPoints* *nTimeSteps;

// 	sqrNorm= 0;
	
// 	// wir starten mit der x richtung
// 	current= 0;
	
// 	wrap= FALSE;
	
// 	/*
// 	 * ******************************************************************************************** 
// 	 * MAIN LOOP
// 	 * ******************************************************************************************** 
// 	 */
	
// 	for (index= 0; index < nPoints; index++)
// 	{
// 		sqrNorm= sqrCache[1]+ coord[0]* spaceStepLength[0]* coord[0]* spaceStepLength[0];
		
// 		spaceBin= getSpaceBin(sqrBoundCache, nSpaceBounds, sqrNorm);

// 		// wenn der abstand > maxbin ist, koennen wir mit dem naechsten punkt weitermachen		
// 		if (spaceBin == -1)
// 		{		
// 			for (time= 0; time < *nTimeSteps; time++)
// 			{
// 				int timeOffset= time* nPoints;
				
// 				timeBin= getTimeBin(timeBounds, nTimeBounds, time);
				
// 				// selbes argument wie oben, wir koennen hier auch mit dem naechsten raumpunkt weitermachen
// 				if (timeBin == -1)
// 					break;
				
// 				for (rep= 0, int repOffset= 0; rep < repeat; rep++; repOffset+= nSpaceTimePoints)
// 				{
// 					sumBins[spaceBin][timeBin]= S[index+ timeOffset+ repOffset]; 
// 					weightBins[spaceBin][timeBin]= W[index+ timeOffset+ repOffset]; 
// 				}
// 			}
// 		}
		
// 		//TODO
// 		//update coord and caches
// 	}
// }

// int getSpaceBin(double *sqrBoundCache, int* nSpaceBounds, double sqrNorm)
// {
// 	int bin= -1;

// 	if (sqrNorm > sqrBoundCache[*nSpaceBounds- 1])
// 		return bin;
	
// 	for (bin= 0; (bin < *nSpaceBounds) && (sqrNorm > sqrBoundCache[bin]); bin++)
// 	{
// 	}
	
// 	return bin;
// }

// int getTimeBin(*timeBounds, *nTimeBounds, time)
// {
// 	int bin= -1;
	
// 	if (*timeBounds * *nTimeBounds < time)
// 		return bin;
	
// 	bin= (time- 1)/ *timeBounds;
	
// 	return bin;
// }

