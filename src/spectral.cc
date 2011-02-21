/*
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Simulation of a random field by spectral turning bands

 Copyright (C) 2000 -- 2011 Martin Schlather, 

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
#include <stdlib.h>
#include "RF.h"
#include "Covariance.h"
#include <assert.h>
#include "MPPstandard.h"


void SetParamSpectral(int *action, int *nLines, int *grid,
		      int *ergodic,
		      int *metropolis,
		      int *nmetro, double *sigma) {
  spectral_param *gp = &(GLOBAL.spectral);
  int d;
  if (*action) {
    for (d=0; d<MAXTBMSPDIM; d++) gp->lines[d] = nLines[d];
    *grid= (int) gp->grid;
    gp->grid=(bool) *grid;
    gp->ergodic=(bool) *ergodic;
    gp->metro = *metropolis;
    gp->nmetro = *nmetro;
    gp->sigma = *sigma;
  } else {
    for (d=0; d<MAXTBMSPDIM; d++) nLines[d]=gp->lines[d];
    *grid= (int) gp->grid;
    *ergodic = (int) gp->ergodic;
    *metropolis = gp->metro;
    *nmetro = gp->nmetro;
    *sigma = gp->sigma;
  }
}

void spectral_destruct(void **S) 
{ 
  if (*S!=NULL) {
    // do NOT delete cov --- only pointer
      // spectral_storage *x; x =  *((spectral_storage**)S);
    free(*S);   
    *S = NULL;
  }
}



int init_spectral(method_type *meth){
  int err=NOERROR;
  cov_model *cov = meth->cov;
  cov_fct *C = CovList + cov->nr;
  globalparam *gp = meth->gp;
  spectral_param *lp = &(gp->spectral);
  spec_covstorage *cs = &(cov->spec);

  assert(meth->S == NULL);
  SET_DESTRUCT(spectral_destruct);
//  if ((meth->S=malloc(sizeof(spectral_storage)))==0){
//    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
//  }
//  s = (spectral_storage*)meth->S;
  if (cov->tsdim > MAXTBMSPDIM) {
    err=ERRORMAXDIMMETH; goto ErrorHandling;
  } 

  cs->sigma = lp->sigma;
  cs->nmetro = lp->nmetro;
  cs->density = NULL;

// 13.10. the following does not seem to be necessary,
// since subsequently it checked for reduceddim<=2;
//  if (key->Time) {err=ERRORTIMENOTALLOWED; goto ErrorHandling;}

  if (cov->tsdim > 4) {
      err = ERRORWRONGDIM;
      goto ErrorHandling;
  }
 
  if (cov->vdim > 1) {
    err = ERRORNOMULTIVARIATE;
    goto ErrorHandling;
  }

//  PrintModelInfo(cov);
//assert(false);
//  CovList[100000].cov((double*) &err, cov, (double*) &err);
//  // 
//assert(false);
 
  if (C->initspectral == NULL) {
      err = ERRORCOVFAILED;
  } else {
    err = C->initspectral(cov);
  }

//  printf("sp err %d %s\n", err, C->name);
  // assert(false);

  
//  printf("spectral no error\n");

 
 ErrorHandling:

  // printf("spectral error %d\n", err);
//  
//PrintModelInfo(cov);
//assert(false);

  return err;
}

void E1(spectral_storage *s, double A, double *e) {
  // ignore grid
  if (s->grid) warning("in 1d no spectral grid implemented yet");
  e[0] = A * ((double) (UNIFORM_RANDOM < 0.5) * 2.0 - 1.0);
}

void E2(spectral_storage *s, double A, double *e) {
  double phi;
  if (s->grid) {
    phi = s->phi2d += s->phistep2d;
  } else {
    phi = TWOPI*UNIFORM_RANDOM;
  }

  // printf("A=%f\n", A);
  // printf("%f\n",phi);
  //printf("%f\n", e[1]);

  e[0] = A * cos(phi);
  e[1] = A * sin(phi);

//  printf("phi=%f A=%f grid=%d step=%f phi2d=%f %f %f\n", 
//	 phi, A, (int) s->grid, 
//	 s->phistep2d, s->phi2d,
//	 e[0], e[1]);
}

void E3(spectral_storage *s, double A, double *e) {
  // ignore grid
  if (s->grid) warning("in 3d no spectral grid implemented yet");
  double phi, psi, Asinpsi;
  phi = TWOPI*UNIFORM_RANDOM;
  psi = PI*UNIFORM_RANDOM;
  Asinpsi = A * sin(psi);
  e[0] = A * cos(psi);
  e[1] = Asinpsi * cos(phi);
  e[2] = Asinpsi * sin(phi);
}

void E(int dim, spectral_storage *s, double A, double *e) {
  switch (dim) {
      case 1 : E1(s, A, e); break;
      case 2 : E2(s, A, e); break;
      case 3 : E3(s, A, e); break;
      default : assert(false);
  }
}




// ******** STANDARD *********************
void do_spectral(method_type *meth, res_type *res) 
    // in two dimensions only!
{  
  globalparam *gp = meth->gpdo;
  spectral_param *lp = &(gp->spectral);
  decision_param *dp = &(gp->decision);
  cov_model *cov = meth->cov;
  cov_fct *C = CovList + cov->nr;
  location_type *loc = meth->loc;
//  spec_covstorage *cs = &(cov->spec);
  spectral_storage s;
  int nt,
    ntot = lp->lines[cov->tsdim],
    origdim = loc->timespacedim,
    spatialdim = loc->spatialdim;
  double inct, 
    *x = loc->x;

  int d, n, nx, gridlenx, gridleny, gridlenz, gridlent;
  double E[MAXTBMSPDIM], // must always be of full dimension, even 
    // if lower dimension is simulated -- simulation algorithm for e
    // might use higher dimensional components
    inc[MAXTBMSPDIM], VV;
  long total = loc->totalpoints,
    spatialpoints = loc->spatialtotalpoints;

  int zaehler,
    every = gp->general.every,
    nthreshold = (every>0) ? every : MAXINT,
    deltathresh = nthreshold;
			
  s.grid = lp->grid;
  s.phistep2d = TWOPI/ (double) ntot; 
  s.phi2d = s.phistep2d * UNIFORM_RANDOM;
  s.ergodic = lp->ergodic;
//  printf("%d %f %f\n", s->grid, s->phistep2d, s->phi2d); assert(false);
  
  for (d=0; d<MAXTBMSPDIM; d++) E[d] = inc[d] = 0.0;
  for (zaehler=0; zaehler<total; zaehler++) res[zaehler]=0.0;
  //the very procedure:

  gridlenx = gridleny = gridlenz = gridlent = 1;
  switch (origdim) {
      case 4 : 
	  gridlent = loc->length[3];
	  // no break;
      case 3 : 
	  gridlenz = loc->length[2];
	  // no break;
     case 2 : 
	  gridleny = loc->length[1];
	  // no break;
      case 1 : 
	  gridlenx = loc->length[0];
	  break;
      default : assert(false);
  }
  for (n=0; n<ntot; n++) {
//      assert(cov->spec.nmetro > 0);

//      printf("%f\n", e[1]);assert(false);

//  PrintModelInfo(cov);
// assert(false);


    C->spectral(cov, &s, E);
    

    //   PrintModelInfo(cov);
    //    printf("spect %d %f %f \n", n, e[0], e[1]);
//     assert(false);

    VV = TWOPI * UNIFORM_RANDOM;

    if (loc->grid) {
      for (d = 0; d < origdim; d++) {
	// VV += meth->grani[d][XSTART] * E[d];
	  inc[d] = loc->xgr[d][XSTEP] * E[d];
//	  printf("E[%d]=%f, inc=%f\n", d, E[d], inc[d]);
      }
      
      if (dp->exactness != DECISION_FALSE) {
	double incx, incy, incz, segz, segy, segx;
	int ny, nz;
	long zaehler;

	incx = inc[0];            
	incy = inc[1];
	incz = inc[2];            
	inct = inc[3];
	
	zaehler = 0; 
	for (nt=0; nt < gridlent; nt++) {	
          for (segz = VV, nz = 0; nz < gridlenz; nz++) {	
	    for (segy = segz, ny = 0; ny < gridleny; ny++) {	
	      for (segx = segy, nx = 0; nx < gridlenx; nx++) {
		// printf("zaehler=%d %f\n", zaehler, segx);
		  res[zaehler++] += (res_type) cos(segx);	  
		segx += incx;
	      }
	      segy += incy;
	    }
	    segz += incz;
	  }
	  VV += inct;
	}
      } else { // ! DECISION_PARAM.exactness
 	double cix, ciy, six, siy, ciz, siz, cit, sit,
	    ct, st, cz, sz, cy, sy, cx, sx, dummy;
	int ny, nz, nt;
	long zaehler;
	
	cix = cos(inc[0]);  six = sin(inc[0]);
	ciy = cos(inc[1]);  siy = sin(inc[1]);
	ciz = cos(inc[2]);  siz = sin(inc[2]);
	cit = cos(inc[3]);  sit = sin(inc[3]);
	ct = cos(VV);
	st = sin(VV);
	
	zaehler = 0; 
	for (nt=0; nt < gridlent; nt++) {	
          for (cz=ct, sz=st, nz = 0; nz < gridlenz; nz++) {	
	    for (cy=cz, sy=sz, ny = 0; ny < gridleny; ny++) {	
	      for (cx=cy, sx=sy, nx = 0; nx < gridlenx; nx++) {
		// printf("zaehler=%d %f\n", zaehler, segx);
		  res[zaehler++] += (res_type) cx;
		dummy = cx * cix - sx * six;
		sx = cx * six + sx * cix;
		cx = dummy;
	      }
	      dummy = cy * ciy - sy * siy;
	      sy = cy * siy + sy * ciy;
	      cy = dummy;
	    }
	    dummy = cz * ciz - sz * siz;
	    sz = cz * siz + sz * ciz;
	    cz = dummy;
	  }
	  dummy = ct * cit - st * sit;
	  st = ct * sit + st * cit;  // 1.9.07: sz
	  ct = dummy;                // 1.9.07: cz
	}
	// cos(a + b) = cos(a) cos(b) - sin(a) sin(b)
        // sin(a + b) = cos(a) sin(b) + cos(b) sin(a)
      }
    } else { // no grid
      int j;
      double psi, cit, sit;
      if (loc->Time) {
	inct = loc->T[XSTEP] * E[spatialdim];
	cit = cos(inct); 
	sit = sin(inct);
	for (j=nx=0; nx<spatialpoints; nx++) { 
	  psi = VV;
	  for (d=0; d<spatialdim; d++) psi += E[d] * x[j++];
	  if (dp->exactness == DECISION_FALSE) {
	    double sx, cx, dummy;
	    cx = cos(psi);
	    sx = sin(psi);
	    for (nt = nx; nt < total; nt += spatialpoints) {
		res[nt] += (res_type) cx;
	      dummy = cx * cit - sx * sit;
	      sx = cx * sit + sx * cit;
	      cx = dummy;
	    }
	  } else {
	    for (nt = nx; nt < total; nt += spatialpoints, psi += inct){
	      res[nt] += (res_type) cos(psi);
	    }	  
	  } // ! exact
	}
      } else {
	double cp = E[0], sp = E[1];
	switch (origdim) {
	    case 1 : 
	      for (nx=0; nx<total; nx++) {	
		res[nx] += (res_type) cos(cp * x[nx] + VV);
	      }
	      break;
	    case 2 :
	      int twonx;
	      for (nx=0; nx<total; nx++) {	
		twonx = nx << 1;
		res[nx] += (res_type) cos(cp * x[twonx] + sp * x[twonx+1] + VV);
	      }
	      break;
	    default :
	      long j;
	      for (j=nx=0; nx<total; nx++) {	
		psi = VV;
		for (d=0; d<origdim; d++) psi += E[d] * x[j++];
		res[nx] += (res_type) cos(psi);
	      }
	}
      }
    } // no grid 
    
    STANDARDUSER;
    
  } // for k<NTOT

  double sqrttwodivbyn;
  C->cov(ZERO, cov, &sqrttwodivbyn);
  sqrttwodivbyn = sqrt(sqrttwodivbyn * 2.0 // 2.0 from spectral simul.
		       / (double) ntot);
  for (nx=0; nx<total; nx++) { res[nx]  *= (res_type) sqrttwodivbyn; }
}



