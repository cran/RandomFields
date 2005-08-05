/*
 Authors
 Martin Schlather, schlath@hsu-hh.de

 Simulation of a random field by spectral turning bands

 Copyright (C) 2000 -- 2005 Martin Schlather, 

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
#include "RFsimu.h"
#include "RFCovFcts.h"
#include <assert.h>


int SPECTRAL_LINES=500; 
bool SPECTRAL_GRID=true;

void SetParamSpectral(int *action,int *nLines, int *grid) {
  if (*action) {
    SPECTRAL_LINES=*nLines; 
    SPECTRAL_GRID=(bool) *grid;
  } else {
    *nLines=SPECTRAL_LINES; 
    *grid= (int) SPECTRAL_GRID;
  }
}

void spectral_destruct(void **S) 
{ 
  if (*S!=NULL) {
      // spectral_storage *x; x =  *((spectral_storage**)S);
    free(*S);   
    *S = NULL;
  }
}

int init_simulatespectral(key_type *key, int m) {
  methodvalue_type *meth;  
  covinfo_type *keycov;
  spectral_storage *s;
  int Xerror, v, actcov;
  double store_variance=-1.0;

  meth = &(key->meth[m]);
  SET_DESTRUCT(spectral_destruct, m);
  if ((meth->S=malloc(sizeof(spectral_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (spectral_storage*)meth->S;
  if (key->Time) {Xerror=ERRORTIMENOTALLOWED; goto ErrorHandling;}

  // noerror_repeat = key->anisotropy;
  actcov=0;
  Xerror = NOERROR;
  for (v=0; v<key->ncov; v++) {
    keycov = &(key->cov[v]);
    if ((keycov->method==SpectralTBM) && (keycov->left)) {
      cov_fct *cov;
      int timespacedim;
      timespacedim = keycov->truetimespacedim;
      if (timespacedim == 0) {Xerror=ERRORTRIVIAL; goto ErrorHandling;}
      if (timespacedim > 2) {Xerror=ERRORWRONGDIM; goto ErrorHandling;}
      assert((keycov->nr >= 0) && (keycov->nr < currentNrCov));
      assert(keycov->param[VARIANCE] >= 0.0);
      meth->covlist[actcov] = v;
      cov = &(CovList[keycov->nr]);
      if (cov->type==ISOHYPERMODEL) {
	  Xerror=ERRORHYPERNOTALLOWED; 
	  goto ErrorHandling;
      }
      if (cov->implemented[SpectralTBM] != IMPLEMENTED) {
        Xerror=ERRORNOTDEFINED; goto ErrorHandling;}
      else s->randomAmplitude[actcov] = cov->spectral;
      if ((Xerror=cov->check(keycov->param, timespacedim, SpectralTBM))!=NOERROR)
	  goto ErrorHandling;
      if (actcov>0) {
        if (key->cov[v-1].op) {
          if (key->cov[v-1].method != SpectralTBM){
            if (GENERAL_PRINTLEVEL>0) 
	      PRINTF("severe error - contact author. %d %d %d %d(%s) %d(%s)\n",
		     v, key->cov[v-1].op, key->ncov, key->cov[v-1].method,
		     METHODNAMES[key->cov[v-1].method],
		     SpectralTBM, METHODNAMES[SpectralTBM]);
	    Xerror=ERRORMETHODMIX; goto ErrorHandling;
	    }
	  Xerror=ERRORNOMULTIPLICATION; goto ErrorHandling;
	}
	if (store_variance == keycov->param[VARIANCE]) actcov++;
	else {
	  keycov->left=true;
	  continue;
	}
      } else {
	store_variance = keycov->param[VARIANCE];
	// noerror_repeat= true; // only for !key->anisotropy ??
	actcov++;
      }
      if ((Xerror=Transform2NoGrid(key, v)) != NOERROR) 
	goto ErrorHandling;
      keycov->left = false;
    }
  } // for v
  meth->actcov=actcov;

  if (actcov==0) { /* no covariance for the considered method found */
    if (Xerror==NOERROR) Xerror=NOERROR_ENDOFLIST;
    goto ErrorHandling;
  }

  // if (noerror_repeat) return NOERROR_REPEAT;
  // else return Xerror;
  return NOERROR_REPEAT;

 ErrorHandling:
  return Xerror;
}

void do_simulatespectral(key_type *key, int m, double *res ) 
  // in two dimensions only!
{  
  methodvalue_type *meth; 
  int k, nx, ny, v;
  double phi=RF_INF, phistep=RF_INF, sqrttwodivbyn; //initialisation not use
  spectral_storage *s;
  covinfo_type *kc;

  meth = &(key->meth[m]);
  assert(key->active);
  s = (spectral_storage*) meth->S;

  for (v=0; v<meth->actcov; v++) {
    kc = &(key->cov[meth->covlist[v]]);
    assert((kc->truetimespacedim>=1) && (kc->truetimespacedim<=2));
  }

  kc = &(key->cov[meth->covlist[0]]);
  sqrttwodivbyn = sqrt(kc->param[VARIANCE] * 2.0 / (double) SPECTRAL_LINES);
  
  /* are the angles of the lines situated on a grid or random?
     if grid, then the step of the angle is 2\pi / SPECTRAL_LINES.
  */
  if (SPECTRAL_GRID) {
    phistep=TWOPI/ (double) SPECTRAL_LINES; 
    phi=phistep*UNIFORM_RANDOM;
  }
  v = (int) (UNIFORM_RANDOM * (double) meth->actcov);
  for (k=0; k<key->totalpoints; k++) { res[k]=0.0; }
  //the very procedure:
  for (k=0; k<SPECTRAL_LINES; k++) {
    double cp, Amp, VV;
    if (SPECTRAL_GRID) {phi+=phistep;} else {phi=TWOPI*UNIFORM_RANDOM;}  
    v = (v + 1) % meth->actcov;
    kc = &(key->cov[meth->covlist[v]]);
    VV = TWOPI * UNIFORM_RANDOM;
    Amp= (s->randomAmplitude[v])(kc->param); 
    cp = Amp * cos(phi);
    if (kc->simugrid) {
      double incx, segx;
      int zaehler;            
      zaehler = 0; 
      incx = kc->x[XSTEPDIM1] * cp;
      if (key->timespacedim==1) {
	segx = VV + kc->x[XSTARTDIM1] * cp;
	for (nx=0; nx<key->length[0]; nx++) {
	  res[zaehler++] += cos(segx);	  
	  segx += incx;
	}
      } else { // key->timespacedim==2
	double incy, segy, sp;            
	sp = Amp * sin(phi); 
	segy = VV + kc->x[XSTARTDIM1] * cp + kc->x[XSTARTDIM2] * sp;
	incy = kc->x[XSTEPDIM2] * sp; 
	for (ny=0; ny<key->length[1]; ny++) {	
	  segx = segy;
	  for (nx=0; nx<key->length[0]; nx++) {
	    res[zaehler++] += cos(segx);	  
	    segx += incx;
	  }
	  segy += incy;
	}
      }
    } else { // no grid
      if (kc->truetimespacedim==1) {
	for (nx=0; nx<key->totalpoints; nx++) {	
	  res[nx] += cos(cp * kc->x[nx] + VV);
	}
      } else { // kc->truetimespacedim==2 
	int twonx;
	double sp;
	sp = Amp * sin(phi); 
	for (nx=0; nx<key->totalpoints; nx++) {	
	  twonx = nx * 2;
	  res[nx] += cos(cp * kc->x[twonx] + sp * kc->x[twonx+1] + VV);
	}
      }	      
    }
  } // for k<SPECTRAL_LINES
  for (nx=0; nx<key->totalpoints; nx++) { res[nx] = res[nx]*sqrttwodivbyn; }  
}
  








