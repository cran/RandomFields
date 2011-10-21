/*
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Metropolis-Hasting for drawing from the spectral density

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
 

void metropolis(cov_model *cov, double *x) {
  spec_covstorage *s = &(cov->spec);
  spectral_density density = s->density;
  int i,d, dim = cov->tsdim,
    n = s->nmetro;
  double p, proposal[MAXMPPDIM], dx,
    sigma = s->sigma;
 
  for (i=0; i<n; i++) {
    dx = density(x, cov);
    for (d=0; d<dim; d++) {
      proposal[d] = x[d] + rnorm(0.0, sigma);
    }
    p = density(proposal, cov) / dx;
    if (p>=1 || UNIFORM_RANDOM < p) {
      for (d=0; d<dim; d++) {
	x[d] = proposal[d];
      }
    } 
  }   
}

#define nBase1 30000
#define maxSearch 100
#define relation 0.001
#define bestfactor1 0.8
#define Factor 2.0
#define maxNiter 100
#define nBase2 8
#define nBase3 500
#define bestfactor2 0.7
int search_metropolis(cov_model *cov) {
  spec_covstorage *s = &(cov->spec);
  double dummy, E, V, Var[maxSearch], Sigma[maxSearch],  
    x[MAXMPPDIM], 
    factor = Factor, 
    maxVar = 0.0,
    MaxVButOne = 0.0, 
    MaxVButTwo = 0.0;
  int n, i, j, k, d, err=NOERROR, 
    dim = cov->tsdim;
    
    i = 0;
    GetRNGstate();

  if (s->sigma <= 0.0) {
      int orig_n = s->nmetro;
      s->nmetro = 1;
      s->sigma = 1.0;
      while (i < maxSearch) {
	for (d=0; d<dim; d++) x[d] = 0.0;
	E = V = 0.0;
	for (j=0; j<nBase1; j++) {
	  metropolis(cov, x);
	  dummy = 0.0;
	  for (d=0; d<dim; d++) dummy += x[d] * x[d];
	  V += dummy;
	  E += sqrt(dummy);
	}
	E /= (double) nBase1;
	V = V / (double) nBase1 - E * E;
	if (V < relation * maxVar) {
	  if (factor > 1.0) factor = 1.0 / factor; else break;
	  s->sigma = factor;
	} else {
	  if (V > maxVar) {
	    MaxVButTwo = MaxVButOne;
	    MaxVButOne = maxVar;
	    maxVar = V;
	  }
	  Sigma[i] = s->sigma;
	  Var[i++] = V;
	  s->sigma *= factor;
	} 
      }
      
      
      if (i >= maxSearch) {
	err = ERRORMETROPOLIS;
	goto ErrorHandling;
      }
      n = 0;
      E = 0.0;
      for (j =0; j<i; j++) {
	if (Var[j] >= bestfactor1 * MaxVButTwo) {
	  if (PL > 3) 
	    PRINTF("%d. sigma=%f var=%f\n", j, Sigma[j], Var[j]); 
	  n++;
	  E += log(Sigma[j]);
	}
      }
      s->sigma = exp(E / (double) n) * 4;
      if (PL > 3) PRINTF("optimal sigma=%f \n", s->sigma); 
    
      
      s->nmetro = orig_n;
    }
    // ****** searching for optimal n
    
    if (s->nmetro <= 0) {
	n = 1;
	maxVar = 0.0;
	for (j=0; n<=maxNiter; n *= 2, j++) { // factor 2 nie aendern!
	    for (d=0; d<dim; d++) x[d] = 0.0;
	    Var[j] = 0.0;
	    for (k=0; k<nBase3; k++) {
		E = V = 0.0;
		for (i=0; i<nBase2; i++) {
		    metropolis(cov, x);
		    dummy = 0.0;
		    for (d=0; d<dim; d++) dummy += x[d] * x[d];
		    V += dummy;
		    E += sqrt(dummy);
		}
		E /= (double) nBase2;
		Var[j] += V / (double) nBase2 - E * E;
	    }
	    Var[j] /= (double) nBase3;
	    if (Var[j] > maxVar) maxVar = Var[j]; 
	}
	maxVar *= bestfactor2;
	for (i=0; i<j; i++) {
	    if (Var[i] >= maxVar) break;
	}
	s->nmetro = 1 << (i+2); // Faktor 4 mehr als ermittelt -- zur Sicherheit
	if (PL > 3) PRINTF("optimal n=%d \n", s->nmetro); 
    }

     
 
 ErrorHandling:
    PutRNGstate();
    return err;
}

		  
