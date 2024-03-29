/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Metropolis-Hasting for drawing from the spectral density

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
#include "RF.h"
 

void metropolis(model *cov, gen_storage *S, double *x) {
  spec_properties *s = &(S->spec);
  spectral_density density = s->density;
  int i,d,
    dim = total_logicaldim(OWN),
    n = s->nmetro;
  double p, proposal[MAXTBMSPDIM], dx,
    *E = s->E,
    sigma = s->sigma;
  if (dim >= MAXTBMSPDIM) BUG;
  //ERR("too high dimension for metropolis.");
 
  for (i=0; i<n; i++) {
    dx = density(E, cov);
    for (d=0; d<dim; d++) {
      proposal[d] = E[d] + rnorm(0.0, sigma);      
    }
    p = density(proposal, cov) / dx;
    if (p>=1 || UNIFORM_RANDOM < p) {
      for (d=0; d<dim; d++) {
	E[d] = proposal[d];
      }
    }
  }
  for (d=0; d<dim; d++) x[d] = E[d];
}



#define nBase1 30000
#define maxSearch 100
#define bestfactor 1.2
#define percent 0.3
#define nBase1percent 0.1
#define Factor 1.5
#define nBase2 150000


// to do ?! in rectangle umschreiben ?!

int search_metropolis(model *cov, gen_storage *S) {
  spec_properties *s = &(S->spec);
  double Sigma[maxSearch], x[MAXTBMSPDIM], oldx[MAXTBMSPDIM], log_s, p,
    prop_factor = S->Sspectral.prop_factor,
    factor = Factor;
  int n, i, j, d, ungleich, zaehler, //Z[maxSearch],  
    D[maxSearch], min, mintol,
    minzaehler, optzaehler, maxzaehler,
    err=NOERROR, 
    dim = total_logicaldim(OWN);
 
  s->nmetro = 1;
  if (s->sigma <= 0.0) {
    s->sigma = 1.0;
    optzaehler = (int) (nBase1 * percent);
    minzaehler = (int) (optzaehler * nBase1percent);
    if (minzaehler < 10) GERR("internal variables insufficient. Please contact author.");
    // GERR("please increase 'spectral.percent' and/or nBase1");
    maxzaehler = nBase1 - minzaehler;
    if (maxzaehler < minzaehler) GERR("Internal variables insufficient. Please contact author.");
    //GERR("please decrease 'spectral.percent'");
    min = nBase1 + 1;
    for (i=0; i < maxSearch; i++) {
      Sigma[i] = s->sigma;
      for (d=0; d<dim; d++) s->E[d] = oldx[d] = 0.0;
      zaehler = 0;
      for (j=0; j<nBase1; j++) {
	metropolis(cov, S, x);
	ungleich = 0;
	for (d=0; d<dim; d++) {
	  ungleich += x[d] != oldx[d];
	  oldx[d] = x[d];
	}
	if (ungleich) zaehler++;
      }
      
      //Z[i] = zaehler;
      D[i] = std::abs(zaehler - optzaehler); // integers
      if (D[i] < min) min = D[i];
      
      if (PL >= PL_DETAILS) {
	PRINTF("s=%10g: z=%d < %d [%d, %d] fact=%10g D=%d %d\n",
	     Sigma[i], zaehler, minzaehler, maxzaehler, optzaehler, factor,
		   D[i], min);
      }
	
      if (zaehler < minzaehler || zaehler > maxzaehler) {
	if (factor > 1.0) factor = 1.0 / factor; else break;
	s->sigma = factor;
      } else {      
	s->sigma *= factor; 
      }
    }
    
    if (i >= maxSearch) GERR("Metropolis search algorithm for optimal sd failed\n -- check whether the scale of the problem has been chosen appropriately");
    
    n = 0;
    log_s = 0.0;
    mintol = (int) (bestfactor * min);
    for (j =0; j<i; j++) {
      if (D[j] <= mintol) {
	if (PL >= PL_DETAILS) {	  
	  PRINTF("%d. sigma=%10g D=%d %d\n", j, Sigma[j], D[j], mintol);
	}
	n++;
	log_s += LOG(Sigma[j]);
      }
    }
      
    s->sigma = EXP(log_s / (double) n);
    if (PL >= PL_DETAILS) { PRINTF("optimal sigma=%10g \n", s->sigma); }
  }
    
  // ****** searching for optimal n
   
  for (d=0; d<dim; d++) s->E[d] = oldx[d] = 0.0;
  zaehler = 0;
  for (j=0; j<nBase2; j++) {
    metropolis(cov, S, x);
    ungleich = 0;
    for (d=0; d<dim; d++) {
      ungleich += x[d] != oldx[d];
      oldx[d] = x[d];
    }
    if (ungleich) zaehler++;
  }
  p = (double) zaehler / (double) nBase2;
  s->nmetro = 1 + (int) FABS(prop_factor / LOG(p));
  if (PL >= PL_DETAILS) {
    for (d=0; d<dim; d++)  PRINTF("d=%d E=%10g\n", d, s->E[d]);
  }
  
  if (PL >= PL_DETAILS) {
    PRINTF("opt.sigma=%10g opt.n=%d (p=%10g, id=%10e, zaehler=%d, dim=%d)\n", 
	   s->sigma, s->nmetro, p, prop_factor, zaehler, OWNLOGDIM(0));
  }
 
 ErrorHandling:
   
    RETURN_ERR(err);
}
