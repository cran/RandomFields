/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by spectral turning bands

 Copyright (C) 2000 -- 2014 Martin Schlather, 

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
#include "Covariance.h"


#define SPECTRAL_LINES (COMMON_GAUSS + 1)
#define SPECTRAL_GRID (COMMON_GAUSS + 2)
#define SPECTRAL_METRO_FACTOR (COMMON_GAUSS + 3)
#define SPECTRAL_SIGMA (COMMON_GAUSS + 4)

int check_spectral(cov_model *cov) {
#define nsel 4
  cov_model *next=cov->sub[0],
    *key =cov->key,
    *sub = key==NULL ? next : key;
 int err,
    dim = cov->tsdim; // taken[MAX DIM],
  spectral_param *gp  = &(GLOBAL.spectral);

 
  ROLE_ASSERT(ROLE_GAUSS);
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
    return ERRORDIM;

  kdefault(cov, SPECTRAL_LINES, gp->lines[dim]); // ok
  kdefault(cov, SPECTRAL_GRID, (int) gp->grid); //ok
  kdefault(cov, SPECTRAL_METRO_FACTOR, gp->prop_factor);
  kdefault(cov, SPECTRAL_SIGMA, gp->sigma); // ok
  if ((err = checkkappas(cov)) != NOERROR) return err;

  // APMI(cov);
  // printf("\n\ncheck_spectral %d %s %d\n", key==NULL, NICK(sub), cov->role);
 
  if (key == NULL) {
#define loops 3
    int i, err2[loops],
      iso[loops] = {ISOTROPIC, SPACEISOTROPIC, ZEROSPACEISO};
    //
    for (i=0; i<loops; i++) {
      if ((err2[i] = CHECK(next, dim,  dim, PosDefType, XONLY, iso[i],
			   SUBMODEL_DEP, cov->role)) == NOERROR) break;
      //PMI(cov);    
      // printf("fhler hier %d\n", cov->role);
    }
    if (i >= loops) return err2[0];
    
  
    if (cov->role != ROLE_BASE && sub->pref[SpectralTBM] == PREF_NONE) 
      return ERRORPREFNONE;    
  } else { //  key != NULL
    // falls hier gelandet, so ruft SPECTRAL SPECTRALINTERN auf!!
    // eventuell mit dazwischenliegenden RMS's

    //PMI(cov, "here");
 
    if ((err = CHECK(sub, dim, dim, ProcessType, XONLY, CARTESIAN_COORD, 
		     SUBMODEL_DEP, ROLE_GAUSS)) != NOERROR) {
      return err;
    }
  }

  // APMI(cov);
  //printf("spectral role=%d %d\n", cov->role, ROLE_BASE);

  setbackward(cov, sub);

  return NOERROR;
}


void range_spectral(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[SPECTRAL_LINES] = 1;
  range->max[SPECTRAL_LINES] = RF_INF;
  range->pmin[SPECTRAL_LINES] = 1;
  range->pmax[SPECTRAL_LINES] = 1000;
  range->openmin[SPECTRAL_LINES] = false;
  range->openmax[SPECTRAL_LINES] = true; 

  range->min[SPECTRAL_GRID] = 0;
  range->max[SPECTRAL_GRID] = 1;
  range->pmin[SPECTRAL_GRID] = 0;
  range->pmax[SPECTRAL_GRID] = 1;
  range->openmin[SPECTRAL_GRID] = false;
  range->openmax[SPECTRAL_GRID] = false; 

  range->min[SPECTRAL_METRO_FACTOR] = 0;
  range->max[SPECTRAL_METRO_FACTOR] = RF_INF;
  range->pmin[SPECTRAL_METRO_FACTOR] = 0;
  range->pmax[SPECTRAL_METRO_FACTOR] = 1e-5;
  range->openmin[SPECTRAL_METRO_FACTOR] = false;
  range->openmax[SPECTRAL_METRO_FACTOR] = true; 

  range->min[SPECTRAL_SIGMA] = 0;
  range->max[SPECTRAL_SIGMA] = RF_INF;
  range->pmin[SPECTRAL_SIGMA] = 0;
  range->pmax[SPECTRAL_SIGMA] = 1e-5;
  range->openmin[SPECTRAL_SIGMA] = false;
  range->openmax[SPECTRAL_SIGMA] = true; 
}

int struct_spectral(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel) {
  if (cov->sub[0]->pref[SpectralTBM] == PREF_NONE) {
    return ERRORPREFNONE;
  }

  ROLE_ASSERT_GAUSS;
 
  return NOERROR;
}


int init_spectral(cov_model *cov, gen_storage *S){
  cov_model *next = cov->sub[0],
    *key =cov->key,
    *sub = key == NULL ? next : key;
  int err=NOERROR;
  spec_properties *s = &(S->spec);
  location_type *loc = Loc(cov);
  //  if ((meth->S=MALLOC(sizeof(spectral_storage)))==0){
  //    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  //  }
  //  s = (spectral_storage*)meth->S;

  //  assert(false); printf("\n\nsdrfsjfd\n");

  if (cov->role == ROLE_COV) {
    return NOERROR;
  }

  ROLE_ASSERT_GAUSS;

  cov->method = SpectralTBM;
  
  if (loc->distances) return ERRORFAILED;
  if (cov->tsdim > MAXTBMSPDIM) {
    err=ERRORMAXDIMMETH; goto ErrorHandling;
  } 
  S->Sspectral.prop_factor = P0(SPECTRAL_METRO_FACTOR);
  s->sigma = P0(SPECTRAL_SIGMA);
  s->nmetro = 0;
  s->density = NULL;
  // 13.10. the following does not seem to be necessary,
  // since subsequently it checked for reduceddim<=2;
  //  if (key->Time) {err=ERRORTIME NOTALLOWED; goto ErrorHandling;}
  if (cov->tsdim >= 4) {
    err = ERRORWRONGDIM;
    goto ErrorHandling;
  }
 
  if (cov->vdim[0] > 1) {
    err = ERRORNOMULTIVARIATE;
    goto ErrorHandling;
  }
 
  //if (key == NULL) APMI(cov->calling);
 
  if ((err = INIT(sub, 0, S)) != NOERROR) 
    goto ErrorHandling;


  err = FieldReturn(cov);

 ErrorHandling:
  cov->simu.active = err == NOERROR;
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
  // print("A=%f\n", A);
  // print("%f\n",phi);
  //print("%f\n", e[1]);
  e[0] = A * cos(phi);
  e[1] = A * sin(phi);
  //  print("phi=%f A=%f grid=%d step=%f phi2d=%f %f %f\n", 
  //	 phi, A, (int) s->grid, 
  //	 s->phistep2d, s->phi2d,
  //	 e[0], e[1]);
}
void E12(spectral_storage *s, int dim, double A, double *e) {
  if (dim==2) E2(s, A, e); 
  else {
    double e1[2];
    E2(s, A, e1);
    e[0]=e1[0];
  }
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
  default: BUG;
  }
}


void do_spectral(cov_model *cov, gen_storage *S) 
// in two dimensions only!
{  
  //  gauss_param *dp = &(gp->gauss);
  double exact = GLOBAL.general.exactness;
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr; // nicht gatternr
  location_type *loc = Loc(cov);

  //PMI(cov->calling->calling);
  assert(S != NULL);
  spec_properties *cs = &(S->spec);
  spectral_storage *s = &(S->Sspectral);
  int d, n, gridlenx, gridleny, gridlenz, gridlent,
    ntot = P0INT(SPECTRAL_LINES),
    origdim = loc->timespacedim,
    spatialdim = loc->spatialdim,
    every = GLOBAL.general.every,
    nthreshold = (every>0) ? every : MAXINT,
    deltathresh = nthreshold;
  double inct, 
    *x = loc->x,
    E[MAXTBMSPDIM], // must always be of full dimension, even 
    // if lower dimension is simulated -- simulation algorithm for e
    // might use higher dimensional components
    inc[MAXTBMSPDIM], VV,
    *res = cov->rf;
  long nt, nx, total = loc->totalpoints,
    spatialpoints = loc->spatialtotalpoints;
  bool loggauss = GLOBAL.gauss.loggauss;

			
  s->grid = P0INT(SPECTRAL_GRID);
  s->phistep2d = TWOPI/ (double) ntot; 
  s->phi2d = s->phistep2d * UNIFORM_RANDOM;
//  print("%d %f %f\n", s->grid, s->phistep2d, s->phi2d); assert(false);
  
  for (d=0; d<MAXTBMSPDIM; d++) E[d] = inc[d] = 0.0;
  for (n=0; n<total; n++) res[n]=0.0;
  //the very procedure:
  gridlenx = gridleny = gridlenz = gridlent = 1;
  if (loc->grid) {
    switch (origdim) {
    case 4 : 
      gridlent = (int) loc->xgr[3][XLENGTH];
      // no break;
    case 3 : 
      gridlenz = (int) loc->xgr[2][XLENGTH];
      // no break;
    case 2 : 
      gridleny = (int) loc->xgr[1][XLENGTH];
      // no break;
    case 1 : 
      gridlenx = (int) loc->xgr[0][XLENGTH];
      break;
    default: BUG;
    }
  }



  for (n=0; n<ntot; n++) {
    C->spectral(next, S, E);
    if (PL > 6)
      PRINTF("spect: %d %f %f %d %d sigma=%f nmetro=%d\n",
	     n, E[0], E[1], loc->grid, loc->caniso  != NULL,
	     cs->sigma, cs->nmetro);
       //assert(n < 3);
   
    //   if (meth->cscale != 1.0) {
    //    double invscale = 1.0 / meth->cscale;       
    //   for (d=0; d<cov->tsdim; d++)  E[d] *= invscale;
    //}
    if (loc->caniso  != NULL) {
      double oldE[MAXTBMSPDIM],
	*A = loc->caniso;
      int m, k, j,
	nrow = origdim,
	nrowdim = nrow * cov->tsdim;
      
      for (d=0; d<cov->tsdim; d++) {
	oldE[d] = E[d];
	E[d] = 0.0;
      }
      
      for (d=0, k=0; d<nrow; d++, k++) {
	for (m=0, j=k; j<nrowdim; j+=nrow) {
	  E[d] += oldE[m++] * A[j];
	}
      }
    }
 
    //  print("spect %d %f %f %d %d %f\n",
    //	   n, E[0], E[1], loc->grid, meth->caniso  != NULL, meth->cscale);
    VV = TWOPI * UNIFORM_RANDOM;
    if (loc->grid) {
      for (d = 0; d < origdim; d++) {
	// VV += meth->grani[d][XSTART] * E[d];
	  inc[d] = loc->xgr[d][XSTEP] * E[d];
//	  print("E[%d]=%f, inc=%f\n", d, E[d], inc[d]);
      }
      
      if (!ISNA(exact) && !exact) {
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
		// print("zaehler=%d %f res=%f %f\n", zaehler, segx, res[zaehler], cos(segx) );
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
	int ny, nz;
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
		//print("zaehler=%d res=%f %f\n", zaehler, res[zaehler], cx );
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
      double psi, cit, sit;
      if (loc->Time) {
	int j;
	inct = loc->T[XSTEP] * E[spatialdim];
	cit = cos(inct); 
	sit = sin(inct);
	for (j=nx=0; nx<spatialpoints; nx++) { 
	  psi = VV;
	  for (d=0; d<spatialdim; d++) psi += E[d] * x[j++];
	  if (!ISNA(exact) && !exact) {
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
    
    if (n >= nthreshold) {				 
      /*   PRINTF("\b\b\b\b%d%%", (int) (n * 100 / ntot)); */ 
      PRINTF("%7d %3d%%\n", n, (int) (n * 100 / ntot)); 
      nthreshold += deltathresh;	       
    } 
    R_CheckUserInterrupt();
    
  } // for k<NTOT
  double sqrttwodivbyn;
  COV(ZERO, next, &sqrttwodivbyn);
  
  sqrttwodivbyn = sqrt(sqrttwodivbyn * 2.0 // 2.0 from spectral simul.
		       / (double) ntot);
  for (nx=0; nx<total; nx++) { 
    //printf("%f %f\n", res[nx], sqrttwodivbyn);
    res[nx]  *= (res_type) sqrttwodivbyn; 
  }

  if (loggauss) {
    for (nx=0; nx<total; nx++) res[nx] = exp(res[nx]);
  }

}
