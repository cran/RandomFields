/*
 Authors
 Martin Schlather, martin.schlather@cu.lu

 Simulation of a random field by spectral turning bands

 Copyright (C) 2000 -- 2004 Martin Schlather, 

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
  switch (*action) {
  case 0 :
    SPECTRAL_LINES=*nLines; 
    SPECTRAL_GRID=(bool) *grid;
    break;
  case 1 :
    *nLines=SPECTRAL_LINES; 
    *grid= (int) SPECTRAL_GRID ;
    if (GetNotPrint)   break;
  case 2 : 
    PRINTF("\nSpectral\n========\nnLines=%d\ngrid=%d\n",
	   SPECTRAL_LINES,SPECTRAL_GRID);
    break;
  default : PRINTF(" unknown action\n"); 
  } 
}

typedef struct spectral_storage {
  param_type param;
  randommeasure randomAmplitude[MAXCOV];
  Real *x;
  int actcov;
  bool grid;
  int TrueDim;
} spectral_storage;

void spectral_destruct(void **S) 
{ 
  if (*S!=NULL) {
    spectral_storage *x;
    x =  *((spectral_storage**)S);
    if (x->x!=NULL) free(x->x);
    free(*S);   
    *S = NULL;
  }
}

int init_simulatespectral(key_type *key, int m)
{
  spectral_storage *s;
  int Xerror,  v, start_param[MAXDIM], index_dim[MAXDIM];
  bool no_last_comp;
  cov_fct *cov;
  
  SET_DESTRUCT(spectral_destruct);
  if ((key->S[m]=malloc(sizeof(spectral_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (spectral_storage*)key->S[m];
  s->x = NULL;

  if (key->Time) {Xerror=ERRORTIMENOTALLOWED; goto ErrorHandling;}

  unsigned short int actcov;
  int covnr[MAXCOV];
  int multiply[MAXCOV];
  Real store_param[TOTAL_PARAM];
  long bytes;
  bool noerror_repeat;

  noerror_repeat = key->anisotropy;
  bytes = key->timespacedim * key->timespacedim;


  //FC1(SpectralTBM, spectral, s->param);
  actcov=0;
  { 
    int v;
    for (v=0; v<key->ncov; v++) {
      if ((key->method[v]==SpectralTBM) && (key->left[v])) {
	// Variance==0 is not eliminated anymore !! -- maybe this could be improved
	key->left[v] = false;
	assert((key->covnr[v] >= 0) && (key->covnr[v] < currentNrCov));
	assert(key->param[v][VARIANCE] >= 0.0);
	covnr[actcov] = key->covnr[v];
	if (CovList[covnr[actcov]].spectral==NULL)
        {Xerror=ERRORNOTDEFINED; goto ErrorHandling;}
	memcpy(s->param[actcov], key->param[v], sizeof(Real) * key->totalparam);
	if (actcov>0) {
	  if ((multiply[actcov-1] = key->op[v-1]) && 
	      key->method[v-1] != SpectralTBM){
	    if (GENERAL_PRINTLEVEL>0) 
	      PRINTF("severe error - contact author. %d %d %d %d (%s) %d (%s)\n",
		     v, key->op[v-1], key->ncov, key->method[v-1],
		     METHODNAMES[key->method[v-1]],
		     SpectralTBM, METHODNAMES[SpectralTBM]);
	    Xerror=ERRORMETHODMIX; goto ErrorHandling;
	  }
	}
  // end FC1
	if (key->anisotropy) {
	  if (actcov>0) { 
	    if (memcmp(store_param, &(s->param[actcov][ANISO]), bytes) ||
		(store_param[VARIANCE] != s->param[actcov][VARIANCE])) {
	      key->left[v]=true; actcov--;}
	  } else {
	    memcpy(store_param, &(s->param[actcov-1][ANISO]), bytes);
	    store_param[VARIANCE] = s->param[actcov][VARIANCE];
	  }
	} else {
	  assert(fabs(key->param[v][SCALE] * key->param[v][INVSCALE]-1.0)
		 < EPSILON);
	  if (actcov>0) { 
	    if (store_param[VARIANCE] != s->param[actcov][VARIANCE]){
	      key->left[v]=true; 
	      noerror_repeat= true;
	      actcov--;
	    }
	  } else {
	    store_param[VARIANCE] = s->param[actcov][VARIANCE];
	  }
	}
  // FC2;
	actcov++;
      }
    }
  }
  if (actcov==0) { /* no covariance for the considered method found */\
    if (key->traditional) Xerror=ERRORNOTINITIALIZED;\
    else Xerror=NOERROR_ENDOFLIST;\
    goto ErrorHandling;\
  }
 //end FC2 

  s->grid = key->grid && !key->anisotropy;
  s->actcov = actcov;
  for (v=0; v<actcov; v++) s->randomAmplitude[v] = CovList[covnr[v]].spectral;
  GetTrueDim(key->anisotropy, key->timespacedim, s->param[0],
	     &(s->TrueDim), &no_last_comp, start_param, index_dim); 
  if ((Xerror=Transform2NoGrid(key, s->param[0], s->TrueDim, 
			       start_param, &(s->x))) != NOERROR)
    goto ErrorHandling; 
  if (s->TrueDim > 2) {Xerror=ERRORWRONGDIM; goto ErrorHandling;}
  if (s->TrueDim == 0) {Xerror=ERRORTRIVIAL; goto ErrorHandling;}
  for (v=0; v<actcov; v++) {
    cov = &(CovList[covnr[v]]);
    if ((cov->check!=NULL) && 
	(Xerror=cov->check(s->param[v], s->TrueDim, SpectralTBM)))
      goto ErrorHandling;
  }
  if (noerror_repeat) return NOERROR_REPEAT;
  else return 0;
 
 ErrorHandling:
  return Xerror;
}

void do_simulatespectral(key_type *key, int m, Real *res ) 
  // in two dimensions only!
{  
  int k, nx, ny, v;
  Real phi=RF_INF, phistep=RF_INF, sqrttwodivbyn, VV;  //initialisation not used
  spectral_storage *s;
  
  assert(key->active);
  s = (spectral_storage*)key->S[m];
  assert((s->TrueDim==1) || (s->TrueDim==2));

  sqrttwodivbyn = sqrt(s->param[0][VARIANCE]*2.0/ (double) SPECTRAL_LINES);  
  
  /* are the angles of the lines situated on a grid or random?
     if grid, then the step of the angle is 2\pi / SPECTRAL_LINES.
  */
  if (SPECTRAL_GRID) {
    phistep=TWOPI/ (double) SPECTRAL_LINES; 
    phi=phistep*UNIFORM_RANDOM;
  }
  {
    register long i,endfor;
    endfor = key->totalpoints;
    for (i=0; i<endfor; i++) { res[i]=0.0; }
  }
  v = (int) (UNIFORM_RANDOM * (double) s->actcov);
  

  //the very procedure:
  if (s->TrueDim==1) {
   if (s->grid) {
      long zaehler;
      Real incx, segx;
      for(k=0; k<SPECTRAL_LINES; k++){
	Real cp,Amp;            
	v = (v+1) % s->actcov;
	Amp=  s->param[v][INVSCALE] * 
	  (s->randomAmplitude[v])(s->param[v] ); 
	if (SPECTRAL_GRID) {phi+=phistep;} else {phi=TWOPI*UNIFORM_RANDOM;}  
	VV = TWOPI*UNIFORM_RANDOM;
	cp = Amp * cos(phi);  
	zaehler=0;      
	incx=key->x[0][XSTEP]*cp; 
	segx = VV;
	for (nx=0; nx<key->length[0]; nx++) {
	  res[zaehler++] += cos(segx);	  
	  segx += incx;
	}
      }
    } else { // no grid
      for(k=0; k<SPECTRAL_LINES; k++){
	Real cp,Amp;
	v = (v+1) % s->actcov;
	Amp= (s->randomAmplitude[v])(s->param[v] ); 
	if (SPECTRAL_GRID) {phi+=phistep;} else {phi=TWOPI*UNIFORM_RANDOM;}  
	VV = TWOPI*UNIFORM_RANDOM;
	cp = Amp * cos(phi);    
	for (nx=0; nx<key->totalpoints; nx++) {	
	  res[nx] += cos(cp * s->x[nx] + VV);
	}
      }		      
    }
  } else {

    if (s->grid) {
      long zaehler;
      Real incx, incy, segx, segy;
      for(k=0; k<SPECTRAL_LINES; k++){
	Real cp,sp,Amp;            
	v = (v+1) % s->actcov;
	Amp=  s->param[v][INVSCALE] * 
	  (s->randomAmplitude[v])(s->param[v] ); 
	if (SPECTRAL_GRID) {phi+=phistep;} else {phi=TWOPI*UNIFORM_RANDOM;}  
	VV = TWOPI*UNIFORM_RANDOM;
	cp=Amp * cos(phi);   sp=Amp * sin(phi); 

	zaehler=0;      
	segy=VV;
	incx=key->x[0][XSTEP]*cp; incy=key->x[1][XSTEP]*sp; 
	for (ny=0; ny<key->length[1]; ny++) {	
	  segx = segy;
	  for (nx=0; nx<key->length[0]; nx++) {
	    res[zaehler++] += cos(segx);	  
	    segx+=incx;
	  }
	  segy += incy;
	}
      }
    } else { // no grid
      int twonx;
      for(k=0; k<SPECTRAL_LINES; k++){
	Real cp,sp,Amp;
	v = (v+1) % s->actcov;
	Amp= (s->randomAmplitude[v])(s->param[v] ); 
	if (SPECTRAL_GRID) {phi+=phistep;} else {phi=TWOPI*UNIFORM_RANDOM;}  
	VV = TWOPI*UNIFORM_RANDOM;
	cp=Amp * cos(phi);   sp=Amp * sin(phi); 
	for (nx=0; nx<key->totalpoints; nx++) {	
	  twonx = nx << 1;
	  res[nx] += cos(cp * s->x[twonx] + sp * s->x[twonx+1] + VV);
	}
      }		      
    }
  }
  for (nx=0; nx<key->totalpoints; nx++) { res[nx] = res[nx]*sqrttwodivbyn; }  
}
  








