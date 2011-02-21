/*
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Simulation of a random field by spatial averaging method

 Copyright (C) 2006 -- 2011 Martin Schlather, 

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


// max betrag der Ableitung von exp(-x^2 + y^2) ist sqrt(2) / 2

// 3.2.1

//double ave_const2(double *p) {
// return pow(M_2_PI / (1.0 - exp(-lp->r2)) *lp->dist * lp->dist,
//	     0.5 * p[EFFECTIVEDIM]);
//}


int init_ave(method_type *meth){
  location_type *loc = meth->loc;
  cov_model *cov = meth->cov;
  bool diag, separatelast, semiseparatelast, quasidiag;
  cov_fct *C = CovList + covOhneGatter(cov)->nr;
  int idx[MAXMPPDIM], err,
    dim = cov->tsdim;

  if (meth->caniso != NULL) {
      analyse_matrix(meth->caniso, loc->timespacedim, dim, &diag, &quasidiag,
		     idx, &semiseparatelast, &separatelast);
      if (!separatelast)
	  return ERRORFAILED;
  } 

  err = init_mppave(meth, dim - 1, true);
  if (err != NOERROR) return ERRORFAILED;
  mpp_storage *s = (mpp_storage *) meth->S;

  s->aniso = getAnisoMatrix(meth);

//  PrintModelInfo(cov);
  //  printf("%d %d %d %s\n", C->aveg, cov, s, C->name); 
  approx_integral(C->avelogg, cov, s); 
  //  if DECISION_PARAM.exactness != DECISION_FALSE return ERRORNOTPOSSIBLE !!
  return NOERROR;
}


void do_ave(method_type *meth, res_type *res) {
    STANDARDDEFINITIONEN(covOhneGatter(meth->cov)); // weder schoen noch sicher!!
  mpp_param *lp = &(gp->mpp);
  // pointsampler location = NULL;

//  static spectral_storage spec = {0.0, 0.0, false, false}; 
  int locations,
      dim = s->dim, //  == spatialdim;
    time = loc->timespacedim * cov->tsdim - 1,
    spatial = loc->spatialtotalpoints;
  double timescale,
      dux[MAXMPPDIM]; // e[MAXMPPDIM], 
  ave_fct f = C->avef;
  ave_logfct logg = C->avelogg;
//  cov_model *next = cov->sub[0];
//  cov_fct *CN =  CovList + next->nr;
  // spectral_do spectral = CN->spectral; 
   
//  C->cov(ZERO, cov, &(s->var));
  
  assert(f != NULL && logg!=NULL);
  Radius =  s->effectiveRadius;
  timescale = loc->xgr[dim][XSTEP] * s->aniso[time];
  factor = sqrt(//s->var *
	       2.0 // aus spectral dichte simulation!
		      * s->integral) ;

  STANDARDINIT;

  
  if (lp->locations > MPP_AUTO) locations = lp->locations;
  else {
      cov_model *dummy=cov;
      locations = C->mpplocations;
//	    printf("%d %s %s %df %d\n",  
//		   dummy->nr, CovList[dummy->nr].name, C->name,
//		   locations, MPP_AUTO);
      while (locations == MPP_AUTO) {
	  if (dummy->nr >=GATTER && dummy->nr<=LASTGATTER) 
	      dummy = dummy->sub[0];
	  locations = CovList[dummy->nr].mpplocations;
//		printf("%s\n",  CovList[dummy->nr].name);
      }
  }
  
  // print("loc  %d %s %ld %d %d\n", lp->locations, C->name, C->mpplocations, 
//	MPP_AUTO, MPP_MISMATCH);

  switch (locations) {
      case MPP_AUTO: case MPP_MISMATCH :
	error("impossible");
	break;
      case MPP_GRID :
	grid_initu(s);
	 s->location = unif_u;
	break;
      case MPP_POISS : 
	poiss_initu(s);
	 s->location = unif_u;
	break;
      case MPP_UNIF :	unif_initu(s);
	 s->location = unif_u;
	break;
      case MPP_GAUSS :
	gauss_initu(s);
	 s->location = gauss_u; 
	break;
      default:
	assert(false);
  }


//  printf("tot=%d\n", s->ntot);
 
  int ntot = s->ntot;
  for (n=0; n<ntot; n++) {

      //     printf("%ld %ld %d \n",coin_ave, randomcoin, n);

    randomcoin(s, cov);
    double 
      spectral = s->c[AVE_SPECTRAL],
      VV =  s->c[AVE_VV],
      inct = spectral * timescale,
      cit = cos(inct), 
      sit = sin(inct);
 
    if (loc->grid && meth->type <= TypeDiag) {
      STANDARDGRID;

// print("insode %d\n", inside);
      if (inside) while (true) {
	double gu, sign;

	for (d=0; d < dim; d++) dux[d] = x[d] - u[d];
///	print("insode XX %d  %d %d \n", inside, dim, cov->tsdim);

	gu = sign = 1.0;
//	warning("logg");
	logg(dux, cov, cov->tsdim - 1, &gu, &sign); 
//print("insdd ode %d\n", inside);
	gu = sign * exp(gu + s->logInvSqrtDens);
//    print("hier\n");
	if (gu > 0.0) {
	    double segt, ct, st;

//	    print(">0.0 %d %s \n", cov->nr, CovList[cov->nr].name);
            segt = VV + spectral * f(dux, cov, cov->tsdim - 1);
// war f(u, cov, cov->tsdim - 1) ?! 20.10.08
	    //   print("hier D\n");
	    ct = gu * cos(segt);
	    st = gu * sin(segt);
	  int zt;
//   print("hier x\n");
	  for (zt = zaehler; zt < total; zt += spatial) {
	      res[zt] += (res_type) ct;	  
	    double dummy = ct * cit - st * sit;
	    st = ct * sit + st * cit;
	    ct = dummy;
	  }	  
	}
//   print("hierf\n");
	STANDARDINKREMENT;
//     print("hier 5\n");
    } // while (true) points

    } else { // no grid
       double *x, gu, sign;
      int j;
      x = meth->space; 
      for (j = zaehler = 0; zaehler<spatial; zaehler++, zaehler += dim) {
	for (d=0; d < dim; d++) 
	  dux[d] = x[j++] - u[d];
	logg(dux, cov, cov->tsdim - 1, &gu, &sign);
	gu = sign * exp(0.25 * dim * s->loginvscale + gu + s->logInvSqrtDens);

	if (gu > 0.0) {
	  double segt = VV + spectral * f(dux, cov, cov->tsdim - 1),
	    ct = gu * cos(segt),
	    st = gu * sin(segt); 
	  int zt;
	  for (zt = zaehler; zt < total; zt += spatial) {
	      res[zt] += (res_type) ct;	  
	    double dummy = ct * cit - st * sit;
	    st = ct * sit + st * cit;
	    ct = dummy;
	  }
	}
      } // for zaehler
    } // no grid

    STANDARDUSER;
  } // for n

  if (factor != 1.0) 
    for (zaehler=0; zaehler<total; zaehler++) 
	res[zaehler] *= (res_type) factor; 

}
 

 
