/*
 Authors
 Martin Schlather, Martin.Schlather@uni-bayreuth.de

 Simulation of a random field by spectral turning bands

 Copyright (C) 2000 Martin Schlather, 

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


int init_simulatespectral(key_type *key)
{
  assert(key->cov->spectral!=NULL);
  assert(key->param[VARIANCE]>=0.0);
  assert(key->param[SILL]>=key->param[VARIANCE]); 
  assert(fabs(key->param[SCALE]*key->param[INVSCALE]-1.0) < EPSILON); 
  if (key->dim != 2)  return ERRORNOTPROGRAMMED;
  key->S = NULL;
  key->destruct = NULL;
  if(key->cov->check!=NULL) return key->cov->check(key); else return 0; 
}

void do_simulatespectral(key_type *key, Real *res END_WITH_GSLRNG) 
  // in two dimensions only!
{  
  int k,nx,ny;
  Real phi,phistep,sqrttwodivbyn,VV;  
  randommeasure randomAmplitude;
   
  assert(key->active);
#ifdef RF_GSL
  assert(RANDOM!=NULL);
#endif

  randomAmplitude = key->cov->spectral;
  sqrttwodivbyn = sqrt(key->param[VARIANCE]*2.0/ (double) SPECTRAL_LINES);  
  
  /* are the angles of the lines situated on a grid or random?
     if grid, then the step of the angle is 2\pi / SPECTRAL_LINES.
  */
  if (SPECTRAL_GRID) {
    phistep=TWOPI/ (double) SPECTRAL_LINES; 
    phi=phistep*UNIFORM_RANDOM;
  }
  
  //the very procedure:
  if (key->grid) {
    long zaehler;
    Real incx, incy, segx, segy;
    for(k=0; k<SPECTRAL_LINES; k++){
      Real cp,sp,Amp;            
      Amp=randomAmplitude(key->param END_WITH_RANDOM); 
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
  } else {
    for(k=0; k<SPECTRAL_LINES; k++){
      Real cp,sp,Amp;
      Amp=randomAmplitude(key->param END_WITH_RANDOM); 
      if (SPECTRAL_GRID) {phi+=phistep;} else {phi=TWOPI*UNIFORM_RANDOM;}  
      VV = TWOPI*UNIFORM_RANDOM;
      cp=Amp * cos(phi);   sp=Amp * sin(phi); 
      for (nx=0; nx<key->totallength; nx++) {
	res[nx] += cos(cp * key->x[0][nx] + sp * key->x[1][nx] + VV);
      }
    }		      
 }
  for (nx=0; nx<key->totallength; nx++) { res[nx] = res[nx]*sqrttwodivbyn; }  
}



