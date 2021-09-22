/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by spectral turning bands

 Copyright (C) 2000 -- 2017 Martin Schlather, 

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
#include <stdio.h>  
#include "def.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "operator.h"
#include "Processes.h"


#define SPECTRAL_LINES (COMMON_GAUSS + 1)
#define SPECTRAL_GRID (COMMON_GAUSS + 2)
#define SPECTRAL_METRO_FACTOR (COMMON_GAUSS + 3)
#define SPECTRAL_SIGMA (COMMON_GAUSS + 4)

int check_spectral(model *cov) {
#define nsel 4
  model *next=cov->sub[0],
    *key =cov->key,
    *sub = key==NULL ? next : key;
  int err;
  spectral_param *gp  = &(GLOBAL.spectral);

  ASSERT_CARTESIAN;
  FRAME_ASSERT_GAUSS_INTERFACE;
  // FRAME_ASSERT(GaussMethod);

  kdefault(cov, SPECTRAL_LINES, gp->lines[OWNXDIM(0)]); // ok
  kdefault(cov, SPECTRAL_GRID, (int) gp->grid); //ok
  kdefault(cov, SPECTRAL_METRO_FACTOR, gp->prop_factor);
  kdefault(cov, SPECTRAL_SIGMA, gp->sigma); // ok
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);

  // APMI(cov);
  // printf("\n\ncheck_spectral %d %s %d\n", key==NULL, NICK(sub), cov->frame);
 
  if (key == NULL) {
#define loops 2
    int i, err2[loops];
    isotropy_type
      iso[loops] = {ISOTROPIC, DOUBLEISOTROPIC};
    //
    for (i=0; i<loops; i++) {
      if ((err2[i] = CHECK(next, OWNLOGDIM(0), OWNXDIM(0), 
			   PosDefType, XONLY, iso[i],
			   SUBMODEL_DEP, GaussMethodType)) == NOERROR) break;
      //PMI(cov);    
      // printf("fhler hier %d\n", cov->frame);
    }
    if (i >= loops) RETURN_ERR(err2[0]);
    
    if (//!hasGaussMethodFrame(cov) ||
	sub->pref[SpectralTBM] == PREF_NONE) {
      //APMI(cov);
      RETURN_ERR(ERRORPREFNONE);
    }
  } else { //  key != NULL
    // falls hier gelandet, so ruft SPECTRAL SPECTRALINTERN auf!!
    // eventuell mit dazwischenliegenden RMS's

    //PMI(cov, "here");
 
    if ((err = CHECK_PASSFRAME(sub, GaussMethodType)) != NOERROR) {
      // if ((err = CHECK(sub, dim, dim, ProcessType, XONLY, CARTESIAN_COORD, 
      //		     SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
      RETURN_ERR(err);
    }
  }

  // APMI(cov);
  //printf("spectral frame=%d %d\n", cov->frame, AnyType);

  setbackward(cov, sub);
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;
}


void range_spectral(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;

  range->min[SPECTRAL_LINES] = 1;
  range->max[SPECTRAL_LINES] = RF_INF;
  range->pmin[SPECTRAL_LINES] = 1;
  range->pmax[SPECTRAL_LINES] = 1000;
  range->openmin[SPECTRAL_LINES] = false;
  range->openmax[SPECTRAL_LINES] = true; 

  booleanRange(SPECTRAL_GRID);

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

int struct_spectral(model *cov, model VARIABLE_IS_NOT_USED **newmodel) {
  if (cov->sub[0]->pref[SpectralTBM] == PREF_NONE) {
    RETURN_ERR(ERRORPREFNONE);
  }

  RETURN_NOERROR;
}


int init_spectral(model *cov, gen_storage *S){
  model *next = cov->sub[0],
    *key =cov->key,
    *sub = key == NULL ? next : key;
  int err=NOERROR,
    dim = ANYDIM;
  spec_properties *s = &(S->spec);
  location_type *loc = Loc(cov);
  //  if ((meth->S=MALLOC(sizeof(spectral_storage)))==0){
  //    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  //  }
  //  s = (spectral_storage*)meth->S;

  //  assert(false); printf("\n\nsdrfsjfd\n");

  if (hasEvaluationFrame(cov)) {
    RETURN_NOERROR;
  }

  cov->method = SpectralTBM;
  
  if (loc->distances) RETURN_ERR(ERRORFAILED);
  if (dim > MAXTBMSPDIM) {
    err=ERRORMAXDIMMETH; goto ErrorHandling;
  } 
  S->Sspectral.prop_factor = P0(SPECTRAL_METRO_FACTOR);
  s->sigma = P0(SPECTRAL_SIGMA);
  s->nmetro = 0;
  s->density = NULL;
  // 13.10. the following does not seem to be necessary,
  // since subsequently it checked for reduceddim<=2;
  //  if (key->Time) {err=ERRORTIME NOTALLOWED; goto ErrorHandling;}
  if (dim >= 4) {
    err = ERRORWRONGDIM;
    goto ErrorHandling;
  }
 
  if (VDIM0 > 1) {
    err = ERRORNOMULTIVARIATE;
    goto ErrorHandling;
  }
 
  if ((err = INIT(sub, 0, S)) != NOERROR) 
    goto ErrorHandling;


  err = ReturnOwnField(cov);

 ErrorHandling:
  cov->simu.active = err == NOERROR;
  RETURN_ERR(err);
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
  // print("A=%10g\n", A);
  // print("%10g\n",phi);
  //print("%10g\n", e[1]);
  e[0] = A * COS(phi);
  e[1] = A * SIN(phi);
  //  print("phi=%10g A=%10g grid=%d step=%10g phi2d=%10g %10g %10g\n", 
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
  Asinpsi = A * SIN(psi);
  e[0] = A * COS(psi);
  e[1] = Asinpsi * COS(phi);
  e[2] = Asinpsi * SIN(phi);
}
void E(int dim, spectral_storage *s, double A, double *e) {
  switch (dim) {
  case 1 : E1(s, A, e); break;
  case 2 : E2(s, A, e); break;
  case 3 : E3(s, A, e); break;
  default: BUG;
  }
}


void do_spectral(model *cov, gen_storage *S) 
// in two dimensions only!
{  
  //  gauss_param *dp = &(gp->gauss);
  usr_bool exact = GLOBAL.general.exactness;
  model *next = cov->sub[0];
  defn *C = DefList + NEXTNR; // nicht COVNR!
  location_type *loc = Loc(cov);

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
  double inct, VV,
     E[MAXTBMSPDIM], // must always be of full dimension, even 
    // if lower dimension is simulated -- simulation algorithm for e
    // might use higher dimensional components
    *x = loc->x,
    oldE[MAXTBMSPDIM] = { 0, 0, 0, 0 },
    inc[MAXTBMSPDIM] = { 0, 0, 0, 0}, 
    *res = cov->rf;
  long nt, nx, total = loc->totalpoints,
    spatialpoints = loc->spatialtotalpoints;
  SAVE_GAUSS_TRAFO;
			
  s->grid = P0INT(SPECTRAL_GRID);
  s->phistep2d = TWOPI/ (double) ntot; 
  s->phi2d = s->phistep2d * UNIFORM_RANDOM;
//  print("%d %10g %10g\n", s->grid, s->phistep2d, s->phi2d); assert(false);
  
  for (d=0; d<MAXTBMSPDIM; d++) E[d] = inc[d] = 0.0;
  for (n=0; n<total; n++) res[n]=0.0;
  //the very procedure:
  gridlenx = gridleny = gridlenz = gridlent = 1;
  if (loc->grid) {
    switch (origdim) {
    case 4 : 
      gridlent = (int) loc->xgr[3][XLENGTH];
      FALLTHROUGH_OK; 		
    case 3 : 
      gridlenz = (int) loc->xgr[2][XLENGTH];
      FALLTHROUGH_OK; 		
    case 2 : 
      gridleny = (int) loc->xgr[1][XLENGTH];
      FALLTHROUGH_OK; 		
    case 1 : 
      gridlenx = (int) loc->xgr[0][XLENGTH];
      break;
    default: BUG;
    }
  }


  for (n=0; n<ntot; n++) {
    C->spectral(next, S, E);
    if (PL > 6) {
      PRINTF("spect: %d %10g %10g %d %d sigma=%10g nmetro=%d\n",
	     n, E[0], E[1], loc->grid, loc->caniso  != NULL,
	     cs->sigma, cs->nmetro);
    }
       //assert(n < 3);
   
    //   if (meth->cscale != 1.0) {
    //    double invscale = 1.0 / meth->cscale;       
    //   for (d=0; d<cov->ts dim; d++)  E[d] *= invscale;
    //}
    if (loc->caniso  != NULL) {
      double 
	*A = loc->caniso;
      int m, k, j,
	dim = ANYDIM,
	nrow = origdim;
      
      for (d=0; d<dim; d++) {
	oldE[d] = E[d];
	E[d] = 0.0;
      }      
      for (d=0, k=0; d<nrow; d++, k++) {
	for (m=0, j=k; m<dim; m++, j+=nrow) {
	  E[d] += oldE[m] * A[j];
	}
      }
    }

    // PMI(cov);   d = 1;
    // printf("E[%d]=%10g %d n=%d\n\n\n\n\n", d, E[d], loc->caniso  != NULL, n);

    //  print("spect %d %10g %10g %d %d %10g\n",
    //	   n, E[0], E[1], loc->grid, meth->caniso  != NULL, meth->cscale);
    VV = TWOPI * UNIFORM_RANDOM;
    if (loc->grid) {
      //      printf("OK\n\n\n\n"); // E!!
      for (d = 0; d < origdim; d++) {
	// VV += meth->grani[d][XSTART] * E[d];
	  inc[d] = loc->xgr[d][XSTEP] * E[d];
//;
//	  printf("E[%d]=%10g\n\n\n\n\n", d, E[d]);
///	  printf("x[%d]=%10g\n\n\n\n", d, loc->xgr[d][XSTEP]);
//	    printf("inc[%d]=%10g\n\n\n\n\n", d, inc[d]);
      }
      

      //for (; d < 4; d++) printf(">>> E[%d]=%10g, inc=%10g\n", d, E[d], inc[d]);
     
       
      if (exact == False) {
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
		// print("zaehler=%d %10g res=%10g %10g\n", zaehler, segx, res[zaehler], COS(segx) );
		  res[zaehler++] += (double) COS(segx);	  
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
	
	cix = COS(inc[0]);  six = SIN(inc[0]);
	ciy = COS(inc[1]);  siy = SIN(inc[1]);
	ciz = COS(inc[2]);  siz = SIN(inc[2]);
	cit = COS(inc[3]);  sit = SIN(inc[3]);
	ct = COS(VV);
	st = SIN(VV);
	
	zaehler = 0; 
	for (nt=0; nt < gridlent; nt++) {	
          for (cz=ct, sz=st, nz = 0; nz < gridlenz; nz++) {	
	    for (cy=cz, sy=sz, ny = 0; ny < gridleny; ny++) {	
	      for (cx=cy, sx=sy, nx = 0; nx < gridlenx; nx++) {
		//print("zaehler=%d res=%10g %10g\n", zaehler, res[zaehler], cx );
		  res[zaehler++] += (double) cx;
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
	// COS(a + b) = COS(a) COS(b) - SIN(a) SIN(b)
        // SIN(a + b) = COS(a) SIN(b) + COS(b) SIN(a)
      }
    } else { // no grid
      double psi, cit, sit;
      if (loc->Time) {
	int j;
	inct = loc->T[XSTEP] * E[spatialdim];
	cit = COS(inct); 
	sit = SIN(inct);
	for (j=nx=0; nx<spatialpoints; nx++) { 
	  psi = VV;
	  for (d=0; d<spatialdim; d++) psi += E[d] * x[j++];
	  if (exact == False) {
	    double sx, cx, dummy;
	    cx = COS(psi);
	    sx = SIN(psi);
	    for (nt = nx; nt < total; nt += spatialpoints) {
		res[nt] += (double) cx;
	      dummy = cx * cit - sx * sit;
	      sx = cx * sit + sx * cit;
	      cx = dummy;
	    }
	  } else {
	    for (nt = nx; nt < total; nt += spatialpoints, psi += inct){
	      res[nt] += (double) COS(psi);
	    }	  
	  } // ! exact
	}
      } else {
	double cp = E[0], sp = E[1];
	switch (origdim) {
	    case 1 : 
	      for (nx=0; nx<total; nx++) {	
		res[nx] += (double) COS(cp * x[nx] + VV);
	      }
	      break;
	    case 2 :
	      int twonx;
	      for (nx=0; nx<total; nx++) {	
		twonx = nx << 1;
		res[nx] += (double) COS(cp * x[twonx] + sp * x[twonx+1] + VV);
	      }
	      break;
	    default :
	      long j;
	      for (j=nx=0; nx<total; nx++) {	
		psi = VV;
		for (d=0; d<origdim; d++) psi += E[d] * x[j++];
		res[nx] += (double) COS(psi);
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
  COV(ZERO(next), next, &sqrttwodivbyn);
  
  sqrttwodivbyn = SQRT(sqrttwodivbyn * 2.0 // 2.0 from spectral simul.
		       / (double) ntot);
  for (nx=0; nx<total; nx++) { 
    //printf("%10g %10g\n", res[nx], sqrttwodivbyn);
    res[nx]  *= (double) sqrttwodivbyn; 
  }

  BOXCOX_INVERSE;
} 

