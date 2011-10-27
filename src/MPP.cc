

/*
 Authors 
 Martin Schlather, martin.schlather@math.uni-goettingen.de 

 Copyright (C) 2001 -- 2011 Martin Schlather, 

 Simulation of random fields by the random coin technique
 

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

/* TO DO
  * noch fast alles auf key->x programmiert statt s->x
  * time component: aenderungen siehe TBM !
*/

#include <math.h>
 
#include "RF.h"
#include "Covariance.h"
#include "MPPstandard.h"

// ----------------------------------------------------------------------
//                  3.2.2, 5.2
// ----------------------------------------------------------------------
// max betrag der Ableitung von exp(-x^2 + y^2) ist sqrt(2) / 2


void SetParamMPP(int *action,
		 int *locations,
		 double *intens, 
		 double *plus,   // ehem MPP_RADIUS ; s->addradius
		 double *relplus, 
//		 double *scale,
		 double *approxzero,
		 double *samplingdist,
		 double *samplingr,
		 double *p,
		 double *beta
		 ) {
  int d;
  mpp_param *gp = &(GLOBAL.mpp);
  if (*action) {
    gp->locations = *locations;
    for (d=0; d<MAXMPPDIM; d++) { 
	gp->intens[d] = intens[d];
	gp->plus[d] = plus[d];
	gp->relplus[d] = relplus[d];
//	gp->scale[d] = scale[d];
    }
    if (*approxzero<=0.0) {
      if (PL>0) 
	PRINTF("\nWARNING! `approx.zero' not a positive number. Ignored.\n");
    } else  gp->approxzero = *approxzero;
    gp->samplingr = *samplingr;
    gp->p = *p;
    gp->beta = *beta;
  } else {
    *locations= gp->locations; // -1:auto; 0: Grid; 1:uniform; 2:Gauss; 3:Poiss
    if (*locations < -1) *locations = -1; 
    else if (*locations>3) error("mpplocations must be in {-1,...,3}");
    for (d=0; d<MAXMPPDIM; d++) {
	intens[d] = gp->intens[d];
	plus[d] = gp->plus[d];
	relplus[d] = gp->relplus[d];
//	scale[d] = gp->scale[d];
    }
    *approxzero = gp->approxzero;
    *samplingdist = gp->samplingdist;
    *samplingr = gp->samplingr;
    *p = gp->p;
    *beta = gp->beta;
  }
}

void mpp_destruct(void **S) {
  // mpp_storage in RFsimu.h !!
  if (*S!=NULL) {
    mpp_storage *x; 
    x = *((mpp_storage**) S);
    if (x->aniso != NULL) free(x->aniso);
    // mpp_storage *s = *((mpp_storage**)S);
    free(*S);
    *S = NULL;
  }
}

void MPP_NULL(mpp_storage* s) {
  int i;
  s->integral = s->integralsq =
    s->effectiveRadius = s->effectivearea = s->plus = s->relplus =
    s->logapproxzero = s->samplingdist = s->average = s->factor = 
    RF_NAN;
  for (i=0; i<MAXMPPDIM; i++) {
    s->min[i] = s->max[i] = s->minsimu[i] = s->maxsimu[i] = s->mean[i] =
      s->sdgauss[i] = RF_NAN;
  }
  for (i=0; i<=LASTC; i++) s->c[i] = RF_NAN;
  s->dim = s->ntot = -1;  // otherwise error possible in InternalGetKeyInfo  
  s->r = s->rsq = s->spherical_adddist =  
    s->maxheight = s->integralpos = RF_NAN;
  s->invscale = 1.0;
  s->loginvscale = 0.0;
  s->logInvSqrtDens = 0.0;
  s->location = NULL;
  s->aniso = NULL;
}


////////////////////////////////////
/// only  locations U + sqrt inverse of density of spatial random point

void gauss_initu(mpp_storage *s) { 
  s->ntot = (int) (s->intensity  * s->effectivearea);
  s->trueintensity = (double) s->ntot;
}

//assert(false);
//  s->logInvSqrtDens = -0.5 * (-0.5 * V * sumsq - (double) s->dim * M_LN_SQRT_2PI 
//		     - log(sqrtdet)
//		     - 0.5 * (double) s->dim * log(V))
// + 0.25 * (double) s->dim * V


void gauss_u(mpp_storage *s, cov_model *cov) { 
  int d,
    dim = s->dim;
  double dummy, 
    sumsq=0.0,
//    V = s->invscale,
//    sqrtV = sqrt(V),
    logsqrtdet=0.0;

  s->getsd(s, cov); // get s->sd[d]
  for (d=0; d<dim; d++) {
    s->u[d] = rnorm(s->mean[d], s->sdgauss[d]);
    dummy = (s->u[d] - s->mean[d]) / s->sdgauss[d];
    sumsq += dummy * dummy;
    logsqrtdet += log(s->sdgauss[d]);
  }
  
  s->logInvSqrtDens = 
    0.5 * (0.5 * (sumsq + (double) s->dim * M_LN_SQRT_2PI) + logsqrtdet);


  if (s->logInvSqrtDens > 8.0 && PL > 5) {
    PRINTF("gaussu: %f %f %f %f\n", 
	   -0.5 * sumsq,
	   - (double) s->dim * M_LN_SQRT_2PI,
	   - logsqrtdet,
	   s->logInvSqrtDens);
  }
}

void unif_initu(mpp_storage *s) { 
  int d,
    dim = s->dim;
  double prod = 1.0;

  s->ntot = (int) (s->intensity * s->effectivearea);
  for (d=0; d<dim; d++) {
    prod *= s->lensimu[d];
  }
  s->logInvSqrtDens = -0.5 * (-log(prod)); 
  s->trueintensity = s->effectivearea / (double) s->ntot;
}

void unif_u(mpp_storage *s, cov_model *cov) {
  // double *min, double *length, double V, int dim, double *u
  int d,
    dim = s->dim;
  double
    V = s->invscale;

  for (d=0; d<dim; d++) {
    s->u[d] = (s->minsimu[d] + UNIFORM_RANDOM * s->lensimu[d]) * V; 
  }
}

void poiss_initu(mpp_storage *s) { 
  double lambda = s->effectivearea * s->intensity;
  s->ntot = (int) rpois(lambda);
  s->trueintensity = s->intensity;
}

//void poiss_u(mpp_storage *s) {
//  // double *min, double *length, double V, int dim, double *u
//}

void grid_initu(mpp_storage *s) {
  int d,
    dim = s->dim;
  double griddist,
    intensity = s->intensity,
    doubledim = (double) dim;
  s->ntot = 1;
  if (intensity > 0.0) {
    s->griddist = griddist = pow(intensity, -1.0 / doubledim);	  
    for (d=0; d<dim; d++) {
      s->u[d] = s->minsimu[d] + UNIFORM_RANDOM * griddist;	 
      s->nmax[d] = (int) ((s->maxsimu[d] - s->u[d]) / griddist) + 1;
      s->n[d] = 0;
      s->ntot *= s->nmax[d];
    }
    s->trueintensity = intensity;    
  } else { // intens == 0, i.e. only a single realisation
    s->griddist = griddist = pow(s->effectivearea, +1.0 / doubledim);
    for (d=0; d<MAXMPPDIM; d++) {
      s->u[d] = s->mean[d];
      s->nmax[d] = 1;
      s->n[d] = 0;
    }
    s->trueintensity = 1.0 / s->effectivearea;
  }
  s->logInvSqrtDens = 0.5 * doubledim * log(griddist);
}

void grid_u(mpp_storage *s, cov_model *cov) {
  int d = 0,
    dim = s->dim;
  double dist = s->griddist;
  s->u[d] += dist;
  (s->n[d])++;
  while ( s->n[d] >= s->nmax[d]){
    s->u[d] = s->minsimu[d];
    s->n[d] = 0;
    if (++d >= dim) break;
    s->u[d] += dist;
    (s->n[d])++;
  }
}

void approx_integral(ave_logfct logg, cov_model *cov, mpp_storage *s) {
// return M_2_PI / (1.0 - exp(-lp->r2) ) * lp->dist * lp->dist;
  double
    maxheight = R_NegInf,
    integral = 0.0,
    integralsq = 0.0;
  double value, sign, xstart[MAXMPPDIM], xend[MAXMPPDIM], x[MAXMPPDIM],
    inc = s->samplingdist,
    R = (s->samplingr <= 0.0) ? s->effectiveRadius : s->samplingr;
  int d,
    dim = s->dim; 
  
  for (d=0; d<dim; d++) {
    xstart[d] = -R + UNIFORM_RANDOM * inc;
    xend[d] = R;
    x[d] = xstart[d];
  }
  
  while (true) {
    logg(x, cov, dim, &value, &sign);
    value = sign * exp(value);
    integral += value; 
    integralsq += value * value;
    if (value > maxheight) maxheight = value;
    d = 0;
    x[d] += inc;
    while (x[d] >= xend[d]) {
      x[d] = xstart[d];
      if (++d >= dim) break;
      x[d] += inc;
    }
    if (d >= dim) break;
  } // while (true) points
  s->integralpos = s->integral = integral;
  s->integralsq = integralsq;
  s->maxheight = maxheight;
}


int init_mppave(method_type *meth, int dim, bool timesep) {
  location_type *loc = meth->loc;
  mpp_storage *s;
  globalparam *gp = meth->gp;
  mpp_param* lp = &(gp->mpp);
  cov_model *cov = meth->cov;
  cov_fct *C = CovList + cov->nr;
  int  err, d,
//    PL = gp->general.printlevel,
    dimM1 = dim - 1;
  long
    total = loc->totalpoints;  

  err = NOERROR;
  SET_DESTRUCT(mpp_destruct);
  assert(meth->S==NULL);

  if ((meth->S=malloc(sizeof(mpp_storage)))==0){
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  if (dim > MAXMPPDIM) {
    err=ERRORMAXDIMMETH; goto ErrorHandling;
  } 

//  printf("initmpp %ld %ld %d\n", meth->S, meth, meth->nr);

  s = (mpp_storage*) meth->S;
  MPP_NULL(s);
  s->dim = dim;              //
  s->intensity = lp->intens[dimM1]; //
  s->samplingdist = lp->samplingdist;
  s->samplingr = lp->samplingr;
  s->c[DRAWMIX_P] = lp->p; 
  s->c[DRAWMIX_BETA] = lp->beta; 
    
//  s->type = (C->average == NULL) ? Spatial 
//    : (lp->random) ? TempRandom : TempGrid;

  Transform2NoGrid(meth, timesep); // s->type != Spatial //
 
  if (meth->type <= TypeDiag &&  loc->grid) {
    for (d=0; d<dim; d++) {
      s->min[d] = meth->grani[d][XSTART];
      s->max[d] = s->min[d] +
	meth->grani[d][XSTEP] * (double) (loc->length[d] - 1);
    }
  } else {
    int ix;
      for (d=0; d<dim; d++) {
      s->min[d]=RF_INF; 
      s->max[d]=RF_NEGINF;
      }
    int endfor = total * dim;
    double *x;
    x = meth->sptime;
    for (ix=0; ix < endfor; ) {
      for (d=0; d<dim; d++, ix++) {
	if (x[ix] < s->min[d]) s->min[d] = x[ix];
	if (x[ix] > s->max[d]) s->max[d] = x[ix];
      }
    }
  }

  // the following parametes must be set BEFORE the C->mppinit call !!
  s->effectiveRadius = 
      s->plus = lp->plus[dimM1];  
  s->relplus = lp->relplus[dimM1];  
  s->logapproxzero = log(lp->approxzero);
  s->effectivearea = 1.0;
  s->integralpos = RF_NAN;


  C->mppinit(s, cov);

// PrintModelInfo(cov); vassert(false);

  s->getsd = CovList[cov->nr].sd;

  for (d=0; d<dim; d++) {
    s->minsimu[d] = s->min[d] - s->effectiveRadius;
    s->maxsimu[d] = s->max[d] + s->effectiveRadius;
    s->lensimu[d] = s->maxsimu[d] - s->minsimu[d];
    s->mean[d] = 0.5 * (s->maxsimu[d] + s->minsimu[d]);
    // s->sd[d] = 0.5 * s->lensimu[d] / lp->scale[dim];

    //   if (s->type = Spatial && s->plus[d] > 0.0 && PL>1)
    //   PRINTF("Note: window has been enlarged by fixed value plus[%d] (%f)\n plus>0 may lead to incorrect simulation results that are difficult to detect.\n", d, s->plus[d]);
    s->effectivearea *= s->lensimu[d];

//    printf("effectivearea %f %f %f %f max= %f %f %f\n", 
//	   s->effectivearea, s->lensimu[d],
//	   s->maxsimu[d], s->minsimu[d], s->max[d], s->min[d],
//	   s->effectiveRadius);
    
    //   assert(!ISNA(s->effectivearea));
  }

  return NOERROR;

 ErrorHandling:
  return err; 
}

int init_mpp(method_type *meth) {
  return init_mppave(meth, meth->cov->tsdim, false);
}


void do_addmpp(method_type *meth, res_type *res ) { 
    STANDARDDEFINITIONEN(meth->cov);
  
  mpp_param *lp = &(meth->gpdo->mpp);
  mpp_get mppget = C->mppget;
  simu_type *simu = meth->simu;
  int  dim = cov->tsdim,
    spatial = loc->totalpoints;


//  C->cov(ZERO, cov, &(s->var));
//  assert(s->var == 1.0);
//  s->var *= meth->cvar;
  // logsqrtvar = log(sqrt(var));

  STANDARDINIT;

  switch (simu->distribution) {
      case DISTR_GAUSS : 
	int locations;

//	printf("%d %d %s\n", locations, C->mpplocations, C->name);

	if (lp->locations > MPP_AUTO) locations = lp->locations;
	else {
	    cov_model *dummy=cov;
	    locations = C->mpplocations;
//	    printf("%d %s %s %df %d\n",  
//		   dummy->nr, CovList[dummy->nr].name, C->name,
//		   locations, MPP_AUTO);
	    while (locations == MPP_AUTO) {
//		printf("mpp %s %d\n", CovList[dummy->nr].name,
//		       CovList[dummy->nr].mpplocations);
		if (dummy->nr >=GATTER && dummy->nr<=LASTGATTER) 
		    dummy = dummy->sub[0];
		locations = CovList[dummy->nr].mpplocations;
//		printf("%s\n",  CovList[dummy->nr].name);
	    }
	}

//	printf("%d %d %s\n", locations, C->mpplocations, C->name);

	switch (locations) {
	    case MPP_AUTO : case MPP_MISMATCH :
	      PrintModelInfo(cov); //
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
	    case MPP_UNIF :
	      unif_initu(s);
	      s->location = unif_u;
	      break;
	    case MPP_GAUSS :
	      gauss_initu(s);
	      s->location = gauss_u; 
	      break;
	    default:
	      assert(false);
	}

	s->factor = sqrt(1.0/ (s->trueintensity * s->integralsq));
	s->average = s->trueintensity * s->integral;

//	printf("%f %f %f\n", s->factor, s->trueintensity, s->integralsq);
//	assert(false);

	break;
      case DISTR_POISSON :
	// not done yet -- it is only a fragment
	// on the R level, it is excluded that the program runs 
	// into this part
	//
	poiss_initu(s);
	s->location = unif_u;
	s->factor = 1.0;
	s->average = 0.0; 
	break;
      default : assert(false);
  }
  int ntot = s->ntot; // ja nicht nach vorne schieben, da s->ntot in switch
  //                     geaendert wird

  //3483 0.023033 1884.955592 3600.000000 0.523599 0.523599
//  printf("%d %f %f %f %f %f\n", ntot, s->factor, s->average, 
//	 s->trueintensity, s->integralsq, s->integral);
//  assert(false);

// printf("%f %f %f %d\n", s->factor, s->trueintensity, s->integralsq, ntot);
// assert(false);
  
  
  for (n=0; n<ntot; n++) 
  {
      if (PL>5) PRINTF("n=%d ntot=%d\n", n, ntot);
//    s->c[MPP_INVSCALE] = mixscale(param); // INVSCALE, da SCALE=1/r waere in
//    //                                    exp(-rU^2) f(r) \D r

//    assert(s->c[MPP_INVSCALE] > 0.0);
//    s->c[MPP_LOGINVSCALE] = log(s->c[MPP_INVSCALE]);
//    randparam(param, r);  // random parameter

    randomcoin(s, cov);
//    assert(false);

    Radius = s->effectiveRadius;
    //   printf("ok %d %d %d \n", n, loc->grid, meth->type, meth->type==TypeAny);
//    printf("%f \n", Radius);
   
    if (loc->grid && meth->type <= TypeDiag) {
      STANDARDGRID;
//      printf("%f %f %f\n", s->u[0], s->u[1], s->u[2]);
      int nn = 0;
      if (inside) while (true) {
	if (false && mppget(x, s->u, cov, s) > 0.0)
	  print("%d %d %f u=%f %f x=%f %f; %f %ld\n", //
		 n, (int)zaehler, (double) res[0], s->u[0], s->u[1], x[0],x[1],
		 mppget(x, s->u, cov, s), total);
	

	res[zaehler] += mppget(x, s->u, cov, s);
	nn++;
	if (!(res[zaehler] < 100000)) {
//	  printf("%d %f\n",  (int) zaehler, mppget(x,s->u,cov,s));
	  assert(false);
	} 
	STANDARDINKREMENT;

      } // while (true) points
//      printf("%d\n", nn);
    } else { // no grid
      double *x;
       x = meth->sptime; 
      for (zaehler=0; zaehler<spatial; zaehler++, x += dim) {
	  res[zaehler] += mppget(x, s->u, cov, s);
      }
    } // no grid
    
    STANDARDUSER;
  } // for n

  double average = s->average;
  factor = s->factor;

//  printf("averager %f %f %f\n", s->average, s->factor, res[0]); assert(false);

  if (s->factor != 1.0 || s->average != 0.0) {  
    for (zaehler=0; zaehler<total; zaehler++) {
//      printf("%d %f ", zaehler, res[zaehler]); 
      res[zaehler] = factor * (res[zaehler] - average);
//      printf(" %f\n", res[zaehler]);
    }
  }
}
  

