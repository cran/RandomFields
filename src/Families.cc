/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of multivariate distribution families

 Copyright (C) 2013 -- 2014 Martin Schlather

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
Foundation, Inc., 59 Temple Place - Suix2te 330, Boston, MA  02111-1307, USA.  
*/




#include "RF.h"
#include "families.h"
#include "Covariance.h"


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/*

  Achtung: falls in *R x gegeben wird, handelt es sich um
  eine bedingte Verteilung. Da wo NA steht wird simuliert.

  todo: Fuer *R x koennte man dies auch programmieren, ist aber noch 
  gemacht

  Falls in *P2sided x[i]==y[i] so wird fuer diese Stellen
    \frac{\partial}{\prod_{diese i's} \partial y_{i_k}} 
        F(y[j_1], y_{i_k}, ....) - F(x[j_1], y_{i_k}, ....)
  betrachtet 

  toDo : erweiterung auf asymmetrische Kurven; dann wird nicht 
  mehr in der Mitte aufgeschnitten, sondern wo linke und rechte
  dichtefunktion sich schneiden

*/




void arcsqrtP(double *x, cov_model  VARIABLE_IS_NOT_USED *cov, double *v) {
  if (*x <= PIHALF) *v = 0.0;
  else {
    *v = atan(sqrt( *x / PIHALF - 1.0)) / PIHALF; // without (5) !!
    // siehe tcf*KS-1.pdf
  }  
}

// --- not corrected yet:


void arcsqrtD(double *x, cov_model  VARIABLE_IS_NOT_USED *cov, double *v) {
  if (*x <= M_PI_2 ) *v = 0.0;
  else {
    *v = SQRT2 / (M_PI * *x * sqrt(*x / M_PI_2 - 2.0)); // wrong
  }
}

void arcsqrtDlog(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v, double *Sign) {
  arcsqrtD(x, cov, v);
  *v = log(*v);
  *Sign = 1.0;
 }

void arcsqrtDinverse(double VARIABLE_IS_NOT_USED *v, cov_model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *left, double VARIABLE_IS_NOT_USED *right) {
  error("Dinverse of arcsqrt unknown");
}

void arcsqrtQ(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  double z = tan(0.5 * PI * *x);
  *v = 2.0 * (z * z + 1.0);
}		   

void arcsqrtR(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v) {
  if (x != NULL) *v = *x;
  else {
    double u = UNIFORM_RANDOM;
    arcsqrtQ(&u, cov, v);
  }
}



//////////////////////////////////////////////////////////////////////

void addVariable(char *name, double *x, int nrow, int ncol, SEXP env) {  
  SEXP  Y;
  int j, size;
    if (ncol == 1) 
      PROTECT(Y = allocVector(REALSXP, nrow));
    else 
      PROTECT(Y = allocMatrix(REALSXP, nrow, ncol));
    size = nrow * ncol;
    //printf("name = %s %f %d %d %d\n",name,x[0],nrow,ncol,isEnvironment(env));
    for (j=0; j<size; j++) REAL(Y)[j] = x[j];
  
    defineVar(install(name), Y, env);
    UNPROTECT(1);
}




#define DISTR_DX 0
#define DISTR_PX 1
#define DISTR_QX 2
#define DISTR_RX 3
#define DISTR_NCOL 4
#define DISTR_NROW 5
#define DISTR_ENV 6
#define DISTR_NAME 7
#define DISTR_LAST DISTR_NAME
void evaluateDistr(cov_model *cov, int which, double *Res) {
  SEXP  res,
    env = PENV(DISTR_ENV)->sexp;
  int size,
    nkappas = CovList[cov->nr].kappas,
    i = nkappas;

  // PMI(cov);
  
  if (cov->ownkappanames != NULL) {
    while(cov->ownkappanames[--i] != NULL) {
      addVariable(cov->ownkappanames[i], P(i), cov->nrow[i], cov->ncol[i], 
		  env);		
    }  
  }

  //  PMI(cov);
  //  print("which = %d\n", which);
  //  eval( ((sexp_type*) P(1])->sexp , env);

  res = eval(PSEXP(which)->sexp, env);
  size = P0INT(DISTR_NROW) * P0INT(DISTR_NCOL);
  for (i=0; i<size; i++) Res[i] = REAL(res)[i]; 
}

void distrD(double *x, cov_model *cov, double *v) {
  addVariable((char*) "x", x, 1, 1, PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_DX, v);
}
void distrDlog(double *x, cov_model *cov, double *v, double *Sign) {
  distrD(x, cov, v);
  *v = log(*v);
  *Sign = 1.0;
}

void distrDinverse(double VARIABLE_IS_NOT_USED  *v, cov_model  VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED  *x, double  VARIABLE_IS_NOT_USED *y) {
  // v is here INPUT variable
  error("'distrDinverse' not programmed yet");  
}

void distrP(double *x, cov_model *cov, double *v) {
  addVariable((char*) "q", x, 1, 1,  PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_PX, v);
}
void distrP2sided(double *x, double *y, cov_model *cov, double *v) {
  int i,
    size = P0INT(DISTR_NROW) * P0INT(DISTR_NCOL),
    dim = cov->xdimown;
  double w;
  assert(dim == cov->tsdim);
  assert(dim == size);
  
  if (dim == 1) {
    double z = x != NULL ? *x : -*y;
    addVariable((char*) "q", &z, 1, 1, PENV(DISTR_ENV)->sexp);
    evaluateDistr(cov, DISTR_PX, &w);
    addVariable((char*) "q", y, 1, 1, PENV(DISTR_ENV)->sexp);
    evaluateDistr(cov, DISTR_PX, v);
    *v -= w;
  } else {
    error("multivariate families of distribution functions not programmed yet");

    double 
      Sign = 1.0,
      *z = cov->Sdollar->z;
    if (z == NULL) z = cov->Sdollar->z =
		     (double*) MALLOC(size * sizeof(double));
 
    for (i=0; i<size; i++) v[i] = 0.0;

    while (true) {
      // Siebformel !! ueber 2^size Moeglichkeiten
      for (i=0; i<size; i++) v[i] = Sign * z[i];
    }
  }
}

void distrQ(double *x, cov_model *cov, double *v) {
  addVariable((char*) "p", x, 1, 1, PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_QX, v);
}
void distrR(double *x, cov_model *cov, double *v) {
  if (x != NULL) error("Conditional distribution not allowed yet");
  addVariable((char*) "n", &ONE, 1, 1, PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_RX, v);
}
void distrR2sided(double *x, double *y, cov_model *cov, double *v) {
  if (x != NULL || y != NULL) 
    error("conditional distribution not allowed yet");
  addVariable((char*) "n", &ONE, 1, 1, PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_RX, v);
}
void kappa_distr(int i, cov_model *cov, int *nr, int *nc){
  //  int dim = cov->tsdim;
  *nc = *nr = i < CovList[cov->nr].kappas ? SIZE_NOT_DETERMINED : -1;
}
int check_distr(cov_model *cov) {
  ROLE_ASSERT(ROLE_DISTR);

  kdefault(cov, DISTR_NROW, 1);
  kdefault(cov, DISTR_NCOL, 1);
  cov->vdim2[0] = P0INT(DISTR_NROW);
  cov->vdim2[1] = P0INT(DISTR_NCOL);
 
  if (cov->Sdollar != NULL && cov->Sdollar->z != NULL)
    DOLLAR_DELETE(&(cov->Sdollar));
  if (cov->Sdollar == NULL) {
    cov->Sdollar = (dollar_storage*) MALLOC(sizeof(dollar_storage));
    DOLLAR_NULL(cov->Sdollar);
  }
  assert(cov->Sdollar->z == NULL);


  return NOERROR;
}
int init_distr(cov_model VARIABLE_IS_NOT_USED *cov, 
	       storage  VARIABLE_IS_NOT_USED *s){
  int err = NOERROR;
  cov->mpp.maxheight = RF_NAN;
  cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0;
  return err;
}
void do_distr_do(cov_model *cov, double *v){ 
  DO_PARAM_MODELS;
  distrR(NULL, cov, v);
}
void range_distr(cov_model *cov, range_type *range){
  range->min[DISTR_NROW] = 1;
  range->max[DISTR_NROW] = 10;
  range->pmin[DISTR_NROW] = 1;
  range->pmax[DISTR_NROW] = 10;
  range->openmin[DISTR_NROW] = false;
  range->openmax[DISTR_NROW] = true;
  
  range->min[DISTR_NCOL] = 1;
  range->max[DISTR_NCOL] = 10;
  range->pmin[DISTR_NCOL] = 1;
  range->pmax[DISTR_NCOL] = 10;
  range->openmin[DISTR_NCOL] = false;
  range->openmax[DISTR_NCOL] = false; 
  
  int i,
    kappas = CovList[cov->nr].kappas;
  for (i=DISTR_LAST + 1; i<kappas; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = 1e10;
    range->pmax[i] = -1e10;
    range->openmin[i] = true;
    range->openmax[i] = true; 
  }

}


#define GAUSS_PARAMETER_BASICS					\
  double VARIABLE_IS_NOT_USED *m = P(GAUSS_DISTR_MEAN),		\
    *sd = P(GAUSS_DISTR_SD)					\

#define GAUSS_PARAMETERS						\
  GAUSS_PARAMETER_BASICS;						\
  int VARIABLE_IS_NOT_USED i, VARIABLE_IS_NOT_USED mi, VARIABLE_IS_NOT_USED  si, \
    len_mean = cov->nrow[GAUSS_DISTR_MEAN],				\
    len_sd = cov->nrow[GAUSS_DISTR_SD],					\
    dim = cov->xdimown;							\
  assert(dim == cov->tsdim);	

#define FOR for(mi=si=i=0; i<dim; i++, mi=(mi + 1) % len_mean, si=(si + 1) % len_sd)

void gaussDlog(double *x, cov_model *cov, double *v, double *Sign) {
  GAUSS_PARAMETERS;
  double returnlog = P0(GAUSS_DISTR_LOG);
  *v = 0.0;
  FOR *v += dnorm(x[i], m[mi], sd[si], returnlog);
  *Sign = 1.0;
}

void gaussD(double *x, cov_model *cov, double *v) {
  GAUSS_PARAMETERS;
  double returnlog = P0(GAUSS_DISTR_LOG);
  if (returnlog) { 
    double Sign;
    gaussDlog(x, cov, v, &Sign);
  } else {        
    *v = 1.0;
    FOR *v *= dnorm(x[i], m[mi], sd[si], returnlog);
  }
}

void gaussDinverse(double *v, cov_model *cov, double *left, double *right) {
  // v is here INPUT variable !!
  GAUSS_PARAMETERS;
  FOR {
    double dummy = - 2.0 * log(*v * SQRTTWOPI * sd[si]);
    if (dummy < 0.0) {
      left[i] = right[i] = m[mi];
    } else {
      dummy = sd[mi] * sqrt(dummy);
      left[i] = m[mi] - dummy;
      right[i] = m[mi] + dummy;
    }
  }
}

void gaussP(double *x, cov_model *cov, double *v) {
  GAUSS_PARAMETERS; 
  double returnlog = P0(GAUSS_DISTR_LOG);
  if (returnlog) {
    *v = 0.0;
    FOR *v += pnorm(x[i], m[mi], sd[si], 1.0, returnlog);
  } else {
    *v = 1.0;
    FOR *v *= pnorm(x[i], m[mi], sd[si], 1.0, returnlog);
  }
}
void gaussP2sided(double *x, double *y, cov_model *cov, double *v) {
  GAUSS_PARAMETERS; 
  double returnlog = P0(GAUSS_DISTR_LOG);
  assert(y != NULL);
  
  if (x == NULL) { // symmetric 2-sided with a=-b
    if (returnlog) {
      FOR *v += y[i] != 0.0 
	? log(2.0 * pnorm(y[i], m[mi], sd[si], 1.0, 0.0) - 1.0)
	: dnorm(y[i], m[mi], sd[si], returnlog);
    } else {
      FOR *v *= y[i] != 0.0 
	? 2.0 * pnorm(y[i], m[mi], sd[si], 1.0, 0.0) - 1.0
 	: dnorm(y[i], m[mi], sd[si], returnlog);
    }
  } else {
    if (returnlog) {
      FOR *v +=  x[i] != y[i]
	? log(pnorm(y[i], m[mi], sd[si], 1.0, 0.0) -
	      pnorm(x[i], m[mi], sd[si], 1.0, 0.0))
	: dnorm(y[i], m[mi], sd[si], returnlog);
    } else {
      FOR *v *=  x[i] != y[i]
	? (pnorm(y[i], m[mi], sd[si], 1.0, 0.0) -
	   pnorm(x[i], m[mi], sd[si], 1.0, 0.0))
	: dnorm(y[i], m[mi], sd[si], returnlog);
    }
  }
}	   
void gaussQ(double *x, cov_model *cov, double *v) {
  GAUSS_PARAMETER_BASICS; 
  double returnlog = P0(GAUSS_DISTR_LOG);
  *v *= qnorm(x[0], m[0], sd[0], 1.0, returnlog);
}		   

void gaussR(double *x, cov_model *cov, double *v) {
  GAUSS_PARAMETERS; 
  if (x == NULL) {
    FOR v[i] = rnorm(m[mi], sd[si]);
  } else {
    FOR v[i] = R_FINITE(x[i]) ? x[i] : rnorm(m[mi], sd[si]);
  }
}


void gaussR2sided(double *x, double *y, cov_model *cov, double *v) {
  GAUSS_PARAMETERS; 
  if (x == NULL) {
    FOR {
      while (fabs(v[i] = rnorm(m[mi], sd[si])) > y[i]); 
    }
  } else {
    FOR { 
      while ((v[i] = rnorm(m[mi], sd[si])) < x[i] || v[i] > y[i]);
    }
  }
}

void kappa_gauss_distr(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i == GAUSS_DISTR_SD || i== GAUSS_DISTR_MEAN) {
    *nc = 1;
    *nr = SIZE_NOT_DETERMINED;
  } else
    *nc = *nr = i == GAUSS_DISTR_LOG ? 1 : -1;
}

int check_gauss_distr(cov_model *cov) {
  ROLE_ASSERT(ROLE_DISTR);
  GAUSS_PARAMETERS;

  if (cov->xdimprev != dim || dim != cov->tsdim) return ERRORDIM;

  if (m == NULL) kdefault(cov, GAUSS_DISTR_MEAN, 0.0);
  if (sd == NULL) kdefault(cov, GAUSS_DISTR_SD, 1.0);
  kdefault(cov, GAUSS_DISTR_LOG, 0.0);

  cov->vdim2[0] = cov->xdimprev;
  cov->vdim2[1] = 1;

  return NOERROR;
}

int init_gauss_distr(cov_model *cov, storage VARIABLE_IS_NOT_USED *s){
  GAUSS_PARAMETERS; 
  
  if (cov->mpp.moments >= 0) {
    cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0; 
    if (cov->mpp.moments >= 1) {
      if (dim > 1) SERR("multivariate moment cannot be calculated");
      cov->mpp.M[1] = cov->mpp.Mplus[1] = m[0];
      if (cov->mpp.moments >= 2) {
	double SD = sd == NULL ? 1.0 : sd[0];      
	cov->mpp.M[2] = cov->mpp.Mplus[2] = SD * SD + m[0] * m[0] ;
      }
    }
  }

  // maxheight based on a Gaussian shape function with the same
  // sd's but normed to 1 at the origin
  cov->mpp.maxheight = intpow(SQRTTWOPI, dim);
  FOR cov->mpp.maxheight *= sd[si];
  cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0;
  return NOERROR;
}

void do_gauss_distr(cov_model *cov, double *v){
  GAUSS_PARAMETERS;
  DO_PARAM_MODELS;
  cov->mpp.maxheight = intpow(SQRTTWOPI, -dim);
  FOR cov->mpp.maxheight /= sd[si];
  gaussR(NULL, cov, v);
}

void range_gauss_distr(cov_model VARIABLE_IS_NOT_USED  *cov, range_type *range){
  range->min[GAUSS_DISTR_MEAN] = RF_NEGINF;
  range->max[GAUSS_DISTR_MEAN] = RF_INF;
  range->pmin[GAUSS_DISTR_MEAN] = -1e8;
  range->pmax[GAUSS_DISTR_MEAN] = 1e8;
  range->openmin[GAUSS_DISTR_MEAN] = true;
  range->openmax[GAUSS_DISTR_MEAN] = true;
  
  range->min[GAUSS_DISTR_SD] = 0;
  range->max[GAUSS_DISTR_SD] = RF_INF;
  range->pmin[GAUSS_DISTR_SD] = 1e-10;
  range->pmax[GAUSS_DISTR_SD] = 1e8;
  range->openmin[GAUSS_DISTR_SD] = true;
  range->openmax[GAUSS_DISTR_SD] = true;
  
  range->min[GAUSS_DISTR_LOG] = 0;
  range->max[GAUSS_DISTR_LOG] = 1;
  range->pmin[GAUSS_DISTR_LOG] = 0;
  range->pmax[GAUSS_DISTR_LOG] = 1;
  range->openmin[GAUSS_DISTR_LOG] = false;
  range->openmax[GAUSS_DISTR_LOG] = false; 
}




#define SPHERIC_SPACEDIM 0
#define SPHERIC_BALLDIM 1
double random_spheric(int tsdim, int balldim) {
  int d;
  double r2 = -1.0;    
  while (r2 < 0.0) {
    r2 = 1.0;
    for (d=tsdim; d<balldim; d++) {
      double dummy = UNIFORM_RANDOM;
      r2 -= dummy * dummy;
    }
  }

  //if (r2 == 1.0) crash();
  // 
  //printf("\n\nrandom %f %d %d\n", r2, tsdim, balldim); //assert(false);
  return 0.5 * sqrt(r2);
}

void sphericD(double VARIABLE_IS_NOT_USED *x, cov_model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED  *v) {
  error("density of 'RRspheric' cannot be calculated yet");
}
void sphericDinverse(double VARIABLE_IS_NOT_USED  *v, cov_model VARIABLE_IS_NOT_USED  *cov, double VARIABLE_IS_NOT_USED  *left, double VARIABLE_IS_NOT_USED  *right) {
  // v is here INPUT variable !!
  error("density of 'RRspheric' cannot be calculated yet");
}
void sphericDlog(double VARIABLE_IS_NOT_USED  *x, cov_model VARIABLE_IS_NOT_USED  *cov, double VARIABLE_IS_NOT_USED  *v, double  VARIABLE_IS_NOT_USED *Sign) {
  error("density of 'RRspheric' cannot be calculated yet");
}
void sphericP(double VARIABLE_IS_NOT_USED  *x, cov_model VARIABLE_IS_NOT_USED   *cov, double VARIABLE_IS_NOT_USED  *v) {
  error("density of 'RRspheric' cannot be calculated yet");
}
void sphericQ(double VARIABLE_IS_NOT_USED *x, cov_model VARIABLE_IS_NOT_USED  *cov, double VARIABLE_IS_NOT_USED  *v) {
  error("density of 'RRspheric' cannot be calculated yet");
}
void sphericR(double *x, cov_model *cov, double *v) {
  int 
    dim = P0INT(SPHERIC_SPACEDIM),
    balldim = P0INT(SPHERIC_BALLDIM);  
  if (x==NULL) *v = random_spheric(dim, balldim); // q[>0] ist E
  else error("conditional distribution cannot be calculated for sphericP.");
}

int check_RRspheric(cov_model *cov) {
  int err;
  ROLE_ASSERT(ROLE_DISTR);

  kdefault(cov, SPHERIC_SPACEDIM, 1);
  kdefault(cov, SPHERIC_BALLDIM, P0INT(SPHERIC_SPACEDIM));
  if ((err = checkkappas(cov)) != NOERROR) return err;

  if (cov->tsdim != 1) SERR("only dimension 1 allowed");

  cov->vdim2[0] = cov->xdimprev;
  cov->vdim2[1] = 1;

  return NOERROR;
}

void range_RRspheric(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SPHERIC_SPACEDIM] = 1; 
  range->max[SPHERIC_SPACEDIM] = RF_INF;
  range->pmin[SPHERIC_SPACEDIM] = 0;
  range->pmax[SPHERIC_SPACEDIM] = 10;
  range->openmin[SPHERIC_SPACEDIM] = false;
  range->openmax[SPHERIC_SPACEDIM] = true;

  range->min[SPHERIC_BALLDIM] = 1; 
  range->max[SPHERIC_BALLDIM] = RF_INF;
  range->pmin[SPHERIC_BALLDIM] = 0;
  range->pmax[SPHERIC_BALLDIM] = 10;
  range->openmin[SPHERIC_BALLDIM] = false;
  range->openmax[SPHERIC_BALLDIM] = true;
}


int init_RRspheric(cov_model *cov, storage VARIABLE_IS_NOT_USED *s) {
  int i, m,
    nm = cov->mpp.moments,
    dim = P0INT(SPHERIC_SPACEDIM),
    balldim = P0INT(SPHERIC_BALLDIM),
    testn = GLOBAL.mpp.n_estim_E;

 double scale, 
    *M = cov->mpp.M,
    *Mplus = cov->mpp.Mplus;
 
  // printf("ball mppM2 %f\n", cov->mpp.M2);
  //printf("%d %d %d\n", testn, GLOBAL.mpp.n_estim_E,
  //	 P0INT(SPHERIC_BALLDIM));// assert(false);
  
  for (M[0] = 1.0, m=1; m<=nm; M[m++] = 0.0);
  for (i=0; i<testn; i++) { // take random sample of scales
    // printf("testn %d %f %d\n", i, scale, testn);
    scale = random_spheric(dim, balldim);
    double dummy = 1.0;
    for (m=1; m<=nm; m++) {
      dummy *= scale;
      M[m] += dummy;
    }
  }
  for (m=1; m<=nm; m++)  {
    M[m] /= (double) testn;
    Mplus[m] = M[m];
    //
    //printf("M: %d %f\n", d, M[d]);
  }

  // exakte Formel; check ob das gleiche wie simu
  if (PL > 1) 
    PRINTF("init_spheric %f %f %f\n", M[nm], 
	   exp( (balldim - dim) * M_LN_SQRT_PI // log(sqrt(pi)) 
		+ lgammafn(0.5 * cov->tsdim + 1) - lgammafn(0.5 * balldim + 1)),
	   exp( - dim * M_LN_SQRT_PI // log(sqrt(pi)) 
		+ lgammafn(0.5 * cov->tsdim + 1))
	   );
	   
  //cov->mpp.refradius = RF_NAN;
  cov->mpp.maxheight = RF_NAN;
  cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0;

  return NOERROR;
}

void do_RRspheric(cov_model *cov, double *v) {
  DO_PARAM_MODELS;
  sphericR(NULL, cov, v);
}



// ***** deterministic ****
#define DETERM_MEAN 0
#define DETERM_PARAMETERS						\
  double *mean = P(DETERM_MEAN);					\
  int i, mi, dim = cov->xdimown,					\
    len_mean = cov->nrow[DETERM_MEAN]				      

#define DETERMFOR for(mi=i=0; i<dim; i++, mi=(mi + 1) % len_mean)


void determD(double *x, cov_model *cov, double *v) { 
  DETERM_PARAMETERS;
  DETERMFOR if (x[i] != mean[mi]) {
    *v = 0.0; 
    return;
  }
  *v = RF_INF;
}

void determDlog(double *x, cov_model *cov, double *v, double *Sign) { 
  determD(x, cov, v);
  *Sign = 1.0;
  *v = log(*v);
}

void determDinverse(double VARIABLE_IS_NOT_USED *v, cov_model *cov, double *left, double *right) {
  DETERM_PARAMETERS;
  DETERMFOR left[i] = right[i] = mean[mi];
}

void determP(double *x, cov_model *cov, double *v) { 
  DETERM_PARAMETERS;
  DETERMFOR {
    if (x[i] < mean[mi]) {
      *v = 0.0; 
      return;
    }
  }
  *v = 1.0; 
}

void determP2sided(double *x, double *y, cov_model *cov, double *v) {   
  DETERM_PARAMETERS;
  assert(y != NULL);  
  if (x == NULL) {
    DETERMFOR {
      if (-y[i] > mean[mi] || y[i] < mean[mi]) {
	*v = 0.0; 
	return;
      }
    }
  } else {
    DETERMFOR {
      if (x[i] > mean[mi] || y[i] < mean[mi]) {
	*v = 0.0; 
	return;
      }
    }
  }
  *v = 1.0; 
}

void determQ(double *x, cov_model *cov, double *v) { 
  determR(x, cov, v);
}

void determR(double *x, cov_model *cov, double *v) { 
  DETERM_PARAMETERS;
  if (x==NULL) DETERMFOR v[i] = mean[i]; 
  else 
    DETERMFOR v[i] = !R_FINITE(x[i]) || x[i] == mean[mi] ? mean[mi] : RF_NAN;
}

void kappa_determ(int i, cov_model *cov, int *nr, int *nc){
  int dim = cov->tsdim;
   *nc = 1;
  *nr = i == 0 ? dim : i == 1 ? 1 : -1;
}

void determR2sided(double *x, double *y, cov_model *cov, double *v) { 
  DETERM_PARAMETERS;
  if (x == NULL) {
    DETERMFOR v[i] = fabs(y[i]) > mean[mi] ? mean[mi] : RF_NAN;    
  } else {
    DETERMFOR v[i] = x[i] < mean[mi] && y[i] > mean[mi] ? mean[mi] : RF_NAN;
  }
}

int check_determ(cov_model *cov) {
  int 
    dim = cov->xdimown;  

  if (cov->xdimprev != dim || dim != cov->tsdim) return ERRORDIM;
  
  if (PisNULL(DETERM_MEAN)) kdefault(cov, DETERM_MEAN, 0.0);
  
  cov->vdim2[0] = dim;
  cov->vdim2[1] = 1;

  return NOERROR;
}


void range_determ(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[DETERM_MEAN] = RF_NEGINF;
  range->max[DETERM_MEAN] = RF_INF;
  range->pmin[DETERM_MEAN] = -1e10;
  range->pmax[DETERM_MEAN] = 1e10;
  range->openmin[DETERM_MEAN] = true;
  range->openmax[DETERM_MEAN] = true;
}

int init_determ(cov_model VARIABLE_IS_NOT_USED *cov, storage VARIABLE_IS_NOT_USED  *s) {
  int err = NOERROR;

  // normed against an extremal gaussian process with shape that 
  // is identically constant 1 in space
  cov->mpp.maxheight = 1.0;
  cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0;
  return err;
}

void do_determ(cov_model *cov, double *v) {
  DO_PARAM_MODELS;
  determR(NULL, cov, v);
}




// ***** Set ****
void setParamD(double *x, cov_model *cov, double *v) { 
VTLG_D(x, cov->sub[SETPARAM_TO], v);
}

void setParamDlog(double *x, cov_model *cov, double *v, double *Sign) { 
  VTLG_DLOG(x, cov->sub[SETPARAM_TO], v, Sign);
}

void setParamDinverse(double *v, cov_model *cov, double *left, double *right) {
  NONSTATINVERSE_D(v, cov->sub[SETPARAM_TO], left, right);
}

void setParamP(double *x, cov_model *cov, double *v) {   
  VTLG_P(x, cov->sub[SETPARAM_TO], v);
}

void setParamP2sided(double *x, double *y, cov_model *cov, double *v) {   
  VTLG_P2SIDED(x, y, cov->sub[SETPARAM_TO], v);
}

void setParamQ(double *x, cov_model *cov, double *v) {   
  VTLG_Q(x, cov->sub[SETPARAM_TO], v);
}

void setParamR(double *x, cov_model *cov, double *v) {   
  VTLG_R(x, cov->sub[SETPARAM_TO], v);
}

void setParamR2sided(double *x, double *y, cov_model *cov, double *v) {   
  VTLG_R2SIDED(x, y, cov->sub[SETPARAM_TO], v);
}

int check_setParam(cov_model *cov) {
  cov_model *to= cov->sub[SETPARAM_TO];
  int err,
    dim = cov->xdimown,
    role = cov->role;
   domain_type dom = cov->domown;
  isotropy_type iso = cov->isoown;

  kdefault(cov, SET_PERFORMDO, true);

  if (cov->xdimprev != dim || dim != cov->tsdim) return ERRORDIM;

  if ((err = CHECK(to, dim, dim, RandomType, dom, iso, 
		     SUBMODEL_DEP, role)) != NOERROR) return err;

  setbackward(cov, to);
  cov->vdim = to->vdim;
  cov->vdim2[0] = to->vdim2[0];
  cov->vdim2[1] = to->vdim2[1];
  cov->deterministic = false;

  TaylorCopy(cov, to);

  return NOERROR;
}


int init_setParam(cov_model VARIABLE_IS_NOT_USED *cov, storage VARIABLE_IS_NOT_USED *s) {
  cov_model *to= cov->sub[SETPARAM_TO];
  int err;
  if ((err = INIT(to, cov->mpp.moments, s)) != NOERROR)
    return err;
  cov->mpp.maxheight = to->mpp.maxheight;
  return NOERROR;
}

void do_setParam(cov_model *cov, double *v) {
  cov_model *next = cov->sub[SETPARAM_TO];
  set_storage *X = cov->Sset;
  cov_model *to = cov->sub[SETPARAM_TO];
  int performDo = P0INT(SET_PERFORMDO);
  DORANDOM(next, v);

  if (performDo > 0) DORANDOM(next, v); 
  if (X->remote != NULL) {
    assert(X->set != NULL);
    X->set(cov->sub[0], X->remote, X->variant);
  }
  if (performDo < 0) DORANDOM(next, v); 


  if (performDo == 0) setParamR(NULL, cov, v);
  cov->mpp.maxheight = to->mpp.maxheight;
}


void range_setParam(cov_model *cov, range_type *range){
  range_setparam(cov, range);
  /*
    range->min[SETPARAM_VARIANT] = 0;
    range->max[SETPARAM_VARIANT] = 0;
    range->pmin[SETPARAM_VARIANT] = 0;
    range->pmax[SETPARAM_VARIANT] = 0;
    range->openmin[SETPARAM_VARIANT] = false;
    range->openmax[SETPARAM_VARIANT] = false;
  */
}




#define LOC_PARAMETER_BASICS			\
  cov_model *next = cov->sub[0];		\
  int dim=cov->xdimown;				\
  double *loc = P(LOC_LOC),			\
    *scale= P(LOC_SCALE)			        
  

#define LOC_PARAMETERS				\
  LOC_PARAMETER_BASICS;				\
  int i, mi, si,				\
    len_loc = cov->nrow[LOC_LOC],		\
    len_scale = cov->nrow[LOC_SCALE];		\
  assert(dim == cov->tsdim)

#define LOCFOR for(mi=si=i=0; i<dim; i++, mi=(mi + 1) % len_loc, si=(si + 1) % len_scale)
 
void locD(double *x, cov_model *cov, double *v) {
  LOC_PARAMETERS;
  double prod = 1.0,
    *z = cov->Sdollar->z;
  if (z == NULL) z = cov->Sdollar->z = (double*) MALLOC(dim * sizeof(double));
    
  LOCFOR { 
    z[i] = (x[i] - loc[mi]) / scale[si];
    prod *= scale[si];
  }
  VTLG_D(z, next, v);
  *v *= prod;
}

void locDlog(double *x, cov_model *cov, double *v, double *Sign) {
  locD(x, cov, v);
  *v = log(*v);
  *Sign = 1.0;
}

void locDinverse(double *v, cov_model *cov, double *left, double *right) {
  LOC_PARAMETERS;
  NONSTATINVERSE_D(v, next, left, right);
  LOCFOR {
    left[i] = left[i] * scale[si] + loc[mi];
    right[i] = right[i] * scale[si] + loc[mi];
  }
}


void locP(double *x, cov_model *cov, double *v) {
  LOC_PARAMETERS;
  double *z = cov->Sdollar->z;
  if (z == NULL) z = cov->Sdollar->z = (double*) MALLOC(dim * sizeof(double));
  LOCFOR z[i] = (x[i] - loc[mi]) / scale[si];
  VTLG_P(z, next, v);
}

void locP2sided(double *x, double *y, cov_model *cov, double *v) {
  LOC_PARAMETERS;
  double *z = cov->Sdollar->z;
  if (z == NULL) z = cov->Sdollar->z = (double*) MALLOC(dim * sizeof(double)); 
  if (x == NULL) {
    LOCFOR z[i] = (y[i] - loc[mi]) / scale[si];
    VTLG_P2SIDED(NULL, z, next, v);
  } else {    
    double *zy = cov->Sdollar->y;
    if (zy == NULL) zy = cov->Sdollar->y=(double*) MALLOC(dim * sizeof(double));
    LOCFOR { 
      z[i] = (x[i] - loc[mi]) / scale[si];
      zy[i] = (y[i] - loc[mi]) / scale[si];
    }
    VTLG_P2SIDED(z, zy, next, v);
  }
}

void locQ(double *x, cov_model *cov, double *v) {
  LOC_PARAMETER_BASICS;  
  if (dim != 1) BUG;
  VTLG_Q(x, next, v);
  v[0] = v[0] * scale[0] + loc[0]; 
}		   

void locR(double *x, cov_model *cov, double *v) {
  LOC_PARAMETERS;  
  VTLG_R(x, next, v);
  if (x == NULL) LOCFOR v[i] = v[i] * scale[si] + loc[mi]; 
  else LOCFOR v[i] = R_FINITE(x[i]) ? x[i] : v[i] * scale[si] + loc[mi]; 
}

void locR2sided(double *x, double *y, cov_model *cov, double *v) {
  LOC_PARAMETERS;  
  VTLG_R2SIDED(x, y, next, v);
  //   printf("locR dim=%d scale=%f %f %f %f -> %f %f %f\n", dim, *scale, y[0], y[1], y[2], v[0], v[1], v[2]);
  LOCFOR {
    v[i] = v[i] * scale[si] + loc[mi]; 
  }
}

void kappa_loc(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i == LOC_SCALE || i== LOC_LOC) {
    *nc = 1;
    *nr = SIZE_NOT_DETERMINED;
  } else
    *nc = *nr = -1;
}

int check_loc(cov_model *cov) {
  ROLE_ASSERT(ROLE_DISTR);
  cov_model *next = cov->sub[0];
  int
    // len_loc = cov->nrow[LOC_LOC],
    // len_scale = cov->nrow[LOC_SCALE],
    dim = cov->xdimown;
  if (cov->xdimprev != dim || dim != cov->tsdim) return ERRORDIM;
  double 
    *loc = P(LOC_LOC),
    *scale = P(LOC_SCALE);
  int err;

  if ((err = CHECK(next, dim, dim, RandomType, cov->domprev, cov->isoprev, 
		     SUBMODEL_DEP, cov->role)) != NOERROR) return err;

  // assert(false); /// not checked yet
 
  if (loc == NULL) kdefault(cov, LOC_LOC, 0.0);
  if (scale == NULL) kdefault(cov, LOC_SCALE, 1.0);
 
  cov->vdim = next->vdim;
  cov->vdim2[0] = next->vdim2[0];
  cov->vdim2[1] = next->vdim2[1];

  if (cov->Sdollar != NULL && cov->Sdollar->z != NULL)
    DOLLAR_DELETE(&(cov->Sdollar));
  if (cov->Sdollar == NULL) {
    cov->Sdollar = (dollar_storage*) MALLOC(sizeof(dollar_storage));
    DOLLAR_NULL(cov->Sdollar);
  }
  assert(cov->Sdollar->z == NULL);

  return NOERROR;
}

int init_loc(cov_model *cov, storage *s){
  // cov_fct *C = CovList + cov->nr;
  LOC_PARAMETERS;
  int err;
  
  if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR) return err;
 
  if (cov->mpp.moments >= 0) {
    cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0; 
    if (cov->mpp.moments >= 1) {
      if (dim > 1) {
	LOCFOR if (scale[si] != 1.0 || loc[mi]!=0.0)
	  SERR("multivariate moment cannot be calculated");
      }
      cov->mpp.M[1] = cov->mpp.M[1] * scale[0] + loc[0];
      cov->mpp.Mplus[1] = loc[0] == 0.0 ? cov->mpp.Mplus[1] * scale[0] : RF_NAN;
      if (cov->mpp.moments >= 2) {
	double ssq = scale[0] * scale[0];
	cov->mpp.M[2] =
	  cov->mpp.M[2] * ssq + loc[0] * (2.0 * cov->mpp.M[1] - loc[0]);
	cov->mpp.Mplus[1] = loc[0] == 0.0 ? cov->mpp.Mplus[1] * ssq : RF_NAN;
      }
    }
  }

  // normed against a shape function with the same scale modification
  cov->mpp.maxheight = next->mpp.maxheight * scale[0];
  cov->mpp.M[0] = next->mpp.M[0];
  cov->mpp.Mplus[0] = next->mpp.Mplus[0];
  return NOERROR;
}

void do_loc(cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  double *scale= P(LOC_SCALE);
  DO_PARAM_MODELS;
  DORANDOM(cov->sub[0], v);
  locR(NULL, cov, v);
  cov->mpp.maxheight = next->mpp.maxheight * scale[0];
}

void range_loc(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[LOC_LOC] = RF_NEGINF;
  range->max[LOC_LOC] = RF_INF;
  range->pmin[LOC_LOC] = -1e8;
  range->pmax[LOC_LOC] = 1e8;
  range->openmin[LOC_LOC] = true;
  range->openmax[LOC_LOC] = true;
  
  range->min[LOC_SCALE] = 0;
  range->max[LOC_SCALE] = RF_INF;
  range->pmin[LOC_SCALE] = 1e-10;
  range->pmax[LOC_SCALE] = 1e8;
  range->openmin[LOC_SCALE] = true;
  range->openmax[LOC_SCALE] = true;

}



#define UNIF_PARAMETER_BASICS			\
  double					\
  *min = P(UNIF_MIN),				\
    *max = P(UNIF_MAX)			

#define UNIF_PARAMETERS				\
  UNIF_PARAMETER_BASICS;			\
  int i, mini, maxi,				\
   len_min = cov->nrow[UNIF_MIN],		\
   len_max = cov->nrow[UNIF_MAX],		\
   dim = cov->xdimown;				\
 assert(dim == cov->tsdim)				

#define UNIFOR  for(mini=maxi=i=0; i<dim; i++, mini=(mini + 1) % len_min, \
		      maxi=(maxi + 1) % len_max)
 
// Abgeschlossene Intervalle!! wichtig !!
void unifD(double *x, cov_model *cov, double *v) {
  UNIF_PARAMETERS;
   
  double area = 1.0;
  UNIFOR {
    if (x[i] < min[mini] || x[i] > max[maxi]) { *v = 0.0; return; }
    else area *= max[maxi] - min[mini];  
  }
  *v = 1.0 / area;
}


void unifDlog(double *x, cov_model *cov, double *v, double *Sign) {
  unifD(x, cov, v);
  *v = log(*v);
  *Sign = 1.0;
}

void unifDinverse(double *v, cov_model *cov, double *left, double *right) {
   UNIF_PARAMETERS;
   double area = 1.0;
   UNIFOR area *= max[maxi] - min[mini];  
   if (*v * area > 1){
     UNIFOR left[i] = right[i] = 0.5 * (max[maxi] + min[mini]);
   } else {
     UNIFOR {
       left[i] = min[mini];
       right[i] = max[maxi];
     //   printf("unif=%d %f %f\n", i, min[mini], max[maxi]);
     }
   }
}

void unifP(double *x, cov_model *cov, double *v) {
  // distribution functions
  UNIF_PARAMETERS;
 
  double area = 1.0;
  UNIFOR {
    if (x[i] <= min[mini]) { *v = 0.0; return; }
    if (x[i] < max[maxi]) 
      area *= (x[i] - min[mini]) / (max[maxi] - min[mini]);  
    // ekse area *= 1.0;
  }
  *v = area;
}


void unifP2sided(double *x, double *y, cov_model *cov, double *v) {
  // distribution functions
  UNIF_PARAMETERS;
  double a, b;

  double area = 1.0;
  UNIFOR {
    a = x != NULL ? x[i] : -y[i];
    if (a != y[i]) {
      if (a < min[mini]) a =  min[mini];
      b = y[i] > max[maxi] ? max[maxi] : y[i];
      if (a >= b) { *v = 0.0; return; }
      if (a < b) area *= (b - a) / (max[maxi] - min[mini]);  
    } else {
      if (a < min[mini] || a > max[maxi]) { *v = 0.0; return; }
      area /= max[maxi] - min[mini];  
    }
  }
  *v = area;
}

void unifQ(double *x, cov_model *cov, double *v) {
  UNIF_PARAMETER_BASICS;
  *v = min[0] + *x * (max[0] - min[0]);
}


void unifR(double *x, cov_model *cov, double *v) { 
  UNIF_PARAMETERS;
  if (x == NULL) {
    UNIFOR v[i] = min[mini] + UNIFORM_RANDOM * (max[maxi] - min[mini]);
  } else {
    UNIFOR {
      v[i] = !R_FINITE(x[i]) 
	? min[mini] + UNIFORM_RANDOM * (max[maxi] - min[mini]) 
	: x[i] >= min[mini] && x[i] <= max[maxi] ? x[i] : RF_NAN;
    }
  }
}


void unifR2sided(double *x, double *y, cov_model *cov, double *v) { 
  UNIF_PARAMETERS;
  double a, b;
 
  UNIFOR {
    a = x != NULL ? (x[i] < min[mini] ? min[mini] : x[i])
      : (-y[i] < min[mini] ? min[mini] : -y[i]);
    b = y[i] > max[maxi] ? max[maxi] : y[i];
    if (a > b) error("parameters of 2-sided unifR out of range");      
    v[i] = a + UNIFORM_RANDOM * (b-a);
  }  
}


void kappa_unif(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i == UNIF_MIN || i == UNIF_MAX) {
    *nc = 1;
    *nr = SIZE_NOT_DETERMINED;
  } else
    *nc = *nr = -1;

  //  printf("i=%d %d %d\n", i, *nr, *nc);
}

int check_unif(cov_model *cov) {

  //PMI(cov->calling);

  ROLE_ASSERT(ROLE_DISTR);

  int 
    //   len_min = cov->nrow[UNIF_MIN],
    // len_max = cov->nrow[UNIF_MAX],
    dim = cov->xdimown;

  if (cov->xdimprev != dim || dim != cov->tsdim) return ERRORDIM;
    
  if (PisNULL(UNIF_MIN)) kdefault(cov, UNIF_MIN, 0.0);
  if (PisNULL(UNIF_MAX)) kdefault(cov, UNIF_MAX, 1.0);

  cov->vdim2[0] = cov->tsdim;
  cov->vdim2[1] = 1;

  //PMI(cov);

  return NOERROR;
}



int init_unif(cov_model *cov, storage VARIABLE_IS_NOT_USED *s){
  UNIF_PARAMETERS;
 
  if (cov->mpp.moments >= 0) {
    cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0; 
    if (cov->mpp.moments >= 1) {
      if (dim > 1) SERR("multivariate moment cannot be calculated");
      cov->mpp.M[1] = 0.5 * (min[0] + max[0]);
      cov->mpp.Mplus[1] = max[0] <= 0.0 ? 0.0 : 0.5 * max[0];
      if (cov->mpp.moments >= 2) {
	cov->mpp.M[2] = max[0] - min[0];
	cov->mpp.M[2] *= cov->mpp.M[2] / 12.0;
      }
    }
  }
   
  // normed against a shape function that is the indicator function 
  // of the respective window
  cov->mpp.maxheight = 1.0;
  UNIFOR  cov->mpp.maxheight *= max[maxi] - min[mini];
  cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0;
    
  return NOERROR;
}


void do_unif(cov_model *cov, double *v){
  UNIF_PARAMETERS;
  DO_PARAM_MODELS;
  unifR(NULL, cov, v);
  cov->mpp.maxheight = 1.0;
  UNIFOR cov->mpp.maxheight /= max[maxi] - min[mini];
}

void range_unif(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[UNIF_MIN] = RF_NEGINF;
  range->max[UNIF_MIN] = RF_INF;
  range->pmin[UNIF_MIN] = -1e8;
  range->pmax[UNIF_MIN] = 1e8;
  range->openmin[UNIF_MIN] = true;
  range->openmax[UNIF_MIN] = true;
  
  range->min[UNIF_MAX] = RF_NEGINF;
  range->max[UNIF_MAX] = RF_INF;
  range->pmin[UNIF_MAX] = -1e8;
  range->pmax[UNIF_MAX] = 1e8;
  range->openmin[UNIF_MAX] = true;
  range->openmax[UNIF_MAX] = true;
}


#define SURFACE(dim) (dim * intpow(2.0, dim)) //of a d-dim cube with edge length 2, i.e. a (dim-1) dimensional surface !


// Ring constanter intensitaet
#define TwoDSimu					     \
  X = UNIFORM_RANDOM * (start + end);			     \
  Y = (2.0 * UNIFORM_RANDOM - 1.0) * (end - start);	     \
  i0 = UNIFORM_RANDOM < 0.5;				     \
  							     \
  x[1 - i0] = Y >= 0 ? Y + start : Y - start;		     \
  x[i0] = ((Y >= 0) xor (i0==1)) ? X -  start : start - X;		     



double VolumeOfCubeRing(double *xsort, double start, double end, int dim, 
			int squeezed_parts) {
  int d,
    red_dim = dim - squeezed_parts;
  assert(red_dim > 0);
  double res = intpow(2.0, dim);
  //  printf("res = %f\n", res);
  for (d=1; d<=squeezed_parts; d++) res *= xsort[d]; 

  // printf("res = %f %f %d %f\n", res, xsort[1], red_dim,
  //	 intpow(end, red_dim) - intpow(start, red_dim));

  return res * (intpow(end, red_dim) - intpow(start, red_dim));
}

double PoweredVolOfCube(double *xsort, double start, double end, double p, 
			int dim, int squeezed_parts) {
  double red_dim = dim - squeezed_parts;
  assert(red_dim > 0);
  double res = red_dim * intpow(2.0, dim), // surface * 2^squeezed
    pPd = p + red_dim;
  int d;

  for (d=1; d<=squeezed_parts; d++) res *= xsort[d];

  //printf("res=%f %f %f %f\n", res,  pow(end, pPd), pow(start, pPd), pPd);

  //return  * intpow(2 * start, squeezed_parts) 
  return res * (pow(end, pPd) - pow(start, pPd)) / pPd;
}

double ExpVolOfCube(double start, double end, double p, double s, int dim, 
		    int squeezed_parts) {

  assert(start >= 0 && end > start);

  // aus evaluate_rectangular
  // die Dichte gleich c * s * p * x^{p-1-(d-1)} * e^{-s x^p} / b_d
  // Fuer nicht squeezed Dimensionen d' = d - squeezed_parts gilt (mit c=1)
  // int_start^end  c * s * p * x^{p-d} * e^{-s x^p} / b_d  d x 
  // =
  // b_{d'} * int_start^end c * s * p * x^{p-d} * x^{d'-1} * e^{-s x^p}/b_d d x 
  // =_{u = s x^p; du= s p x^{p-1} dx}
  // b_{d'}/b_d * int_{s start^p}^{s end^p} * x^{d'-d} * e^{-u}   d u
  // =
  // b_{d'}/b_d * s^{(d - d')/p} * 
  //                  int_{s start^p}^{s end^p} * u^{(d'-d) / p} * e^{-u}  d u

  double 
    red_dim = dim - squeezed_parts,
    a = - squeezed_parts / p + 1.0;
  assert(red_dim > 0);

  //  printf("s=%f start=%f p=%f sq=%d a=%f %f\n", s, start, p, squeezed_parts, a, s * pow(start, p));
  
 
  //  printf("expvol %f %f %f a=%f ret=%f\n", 
  //       SURFACE(red_dim) / SURFACE(dim),  pow(s, squeezed_parts / p),
  //   incomplete_gamma(s * pow(start, p), s * pow(end, p), a), a,
  //   SURFACE(red_dim) / SURFACE(dim) * pow(s, squeezed_parts / p) * 
  // incomplete_gamma(s * pow(start, p), s * pow(end, p), a));

  return 
    SURFACE(red_dim) / SURFACE(dim) * pow(s, squeezed_parts / p) * 
    incomplete_gamma(s * pow(start, p), s * pow(end, p), a);
}

void RandomPointOnCubeRing(double start, double end, int dim, double *x) {
  double twostart = 2.0 * start, X, Y;
  int i0;

  //  printf("start = %f %f\n", start, end);
  assert(start >= 0 && end > start);

  // todo linear on the ring instead of constant ?!

  // todo general dimension

  switch(dim) {
  case 1: {
    x[0] = (2.0 * UNIFORM_RANDOM - 1.0) * (end - start);
    if (x[0] < 0.0) x[0] -= start; else x[0] += start;
  }
    break;
  case 2: TwoDSimu
    break;
  case 3: {
    // Aufsplitten in Roehre mit obigen 2-d Schnitt und Laenge 2*start
    // und den beiden noch fehlenden Kappen der Groesse 4 end^2
    double massTube, massEnds, Z;

    massTube = 4.0 * (start + end) * (end - start) * twostart;
    massEnds = 2.0 * end;
    massEnds = 2.0 * massEnds * massEnds;
    assert(massTube > 0);
    assert(massEnds > 0);
    
    X = UNIFORM_RANDOM * (massTube + massEnds);
    if (X < massTube) { // Roehre
      TwoDSimu 
	x[2] = (2.0 * UNIFORM_RANDOM - 1.0) * start;
    } else { // kappen
      x[0] = (2.0 * UNIFORM_RANDOM - 1.0) * end;
      x[1] = (2.0 * UNIFORM_RANDOM - 1.0) * end;
      Z = (2.0 * UNIFORM_RANDOM - 1.0) * (end - start);
      x[2] = Z > 0 ? Z + start : Z - start;
    }
    //printf("x= %f %f %f  start=%f %f Z=%f\n", x[0], x[1],x[2], start, end, Z);
  }
    break;
  default : BUG;
  }
  
  //printf("x0= %f %f %f dim=%d\n", x[0], x[1],x[2],dim);
}
 

#define TwoDSurface		  \
  X = UNIFORM_RANDOM * 2 * dist;  \
  if (X > 4 * dist) {		  \
    if (X > 6 * dist) {		  \
      x[0] = -dist;		  \
      x[1] = X - 7 * dist;	  \
    } else {			  \
      x[1] = dist;		  \
      x[0] = X - 5 * dist;	  \
    }				  \
  } else {			  \
    if (X > 2 * dist) {		  \
      x[0] = dist;		  \
      x[1] = X - 3 * dist;	  \
    } else {			  \
      x[1] = -dist;		  \
      x[0] = X - 1 * dist;	  \
    }				  \
  }				  
 
void RandomPointOnCubeSurface(double dist, int dim, double *x) {
  double X;

   assert(dist >= 0);
 // todo general dimension

  // printf("RandomPointOnCubeSurface dim=%d\n", dim);

  switch(dim) {
  case 1:
    x[0] = UNIFORM_RANDOM < 0.5 ? dist : -dist;
    break;
  case 2:  TwoDSurface
      break;
  case 3: {
    if ( (X = UNIFORM_RANDOM * 6) > 2 ) { // Roehre
      TwoDSurface
	x[2] = (2.0 * UNIFORM_RANDOM - 1.0) * dist;
    } else { // Kappen
      x[0] = (2.0 * UNIFORM_RANDOM - 1.0) * dist;
      x[1] = (2.0 * UNIFORM_RANDOM - 1.0) * dist;
      x[2] = X > 1 ? -dist : dist;
    }
  }
    break;
  default : BUG;
  }
}
 

// rectangular distribution
// Marco's Idee um Kirstins Sachen zu implementieren

#define IDX_INNER 0
#define IDX_STEPS 1
#define IDX_OUTER NSTEP + IDX_STEPS 
#define ASSIGN_IDX_INNER (-IDX_STEPS)
#define ASSIGN_IDX_OUTER (-IDX_STEPS - 1)
#define ASSIGN_IDX_MISMATCH (-IDX_STEPS - 2)
#define INNER rect->inner
#define INNER rect->inner
#define INNER_CONST rect->inner_const
#define INNER_POW rect->inner_pow
#define OUTER rect->outer
#define OUTER_CONST rect->outer_const
#define OUTER_POW  rect->outer_pow
#define OUTER_POW_CONST rect->outer_pow_const
#define STEP   rect->step
#define NSTEP  rect->nstep
#define VALUE(IDX)  rect->value[IDX_STEPS + IDX]
#define INNER_CUM  rect->weight[IDX_INNER]
#define OUTER_CUM  rect->weight[IDX_OUTER]
#define WEIGHT(IDX)  rect->weight[IDX_STEPS + IDX]
#define TMP (rect->tmp_n)
#define TMP_WEIGHT rect->tmp_weight
#define RIGHT_END(IDX) rect->right_endpoint[IDX]
#define SQUEEZED(IDX) rect->squeezed_dim[IDX]
#define ASSIGNMENT(IDX) rect->assign[IDX] // for control only

#define RECTANGULAR_PARAMETER_BASICS		\
  rect_storage *rect = cov->Srect;		\
  int VARIABLE_IS_NOT_USED  dim = cov->xdimown;	\
  assert(dim == cov->tsdim)				

#define RECTANGULAR_PARAMETERS				\
  RECTANGULAR_PARAMETER_BASICS;				\
  cov_model VARIABLE_IS_NOT_USED *next = cov->sub[0];	\
  if (rect == NULL) BUG
 
#define RECTFOR  for(i=0; i<dim; i++)

#define DEBUG_CUMSUM false
void CumSum(double *y, bool kurz, cov_model *cov, double *cumsum) {
  // kurz: kann eine Abkurzung verwendet werden, da rect->weight
  // bereits bekannt ist ?
  RECTANGULAR_PARAMETERS;
  double min, 
    *ysort = rect->ysort,    
    *zneu = rect->z;    
  int d, kstep, dM1,
    start_cumul = 1,
    dimP1 = dim + 1,
    One = 1;
  bool use_weights = cumsum != rect->weight;
  
  //  printf("y=%f \n", y[0]);
  ysort[0] = 0.0;
  for (d=0; d<dim; d++) {
    //
    //   printf("%d y=%f\n", d, y[d]);
    if (y[d] <= 0.0) {  BUG;}
    ysort[d+1] = fabs(y[d]);
    //printf("y=%f\n", y[d]);
  }

  {
    int *i = rect->i;
    if (dim > 1) {
      Ordering(ysort, &dimP1, &One, i);
      assert(dim < 3 || 
	     (ysort[i[0]] <= ysort[i[1]] && ysort[i[1]] <= ysort[i[2]]));
    } else {
      i[0] = 0;
      i[1] = 1;      
    }
    assert(i[0] == 0 && i[1] > 0);
    //for (d=0; d<=dim; d++) printf("%d,%f ", i[d], z[d]); printf("\n");
    for(d=0; d<dimP1; d++) zneu[d] = ysort[i[d]];
    for(d=0; d<dimP1; d++) ysort[d] = zneu[d];
  }

  TMP = 0;
  d = 1;
  dM1 = d - 1;
  if (zneu[d] >= INNER && use_weights) {
    kstep = zneu[d] >= OUTER ? NSTEP : (zneu[d] - INNER) / STEP;  

    //
    if (DEBUG_CUMSUM) printf("kstep %d %f\n", kstep, (zneu[d] - INNER)/STEP);//

    if (kurz) {
      if (zneu[d] == RF_INF) {
	cumsum[TMP] = OUTER_CUM;
	TMP++;
	//PMI(cov, -1); print("%d %d %f\n",IDX_OUTER, NSTEP, OUTER_CUM);
	assert(cumsum[TMP-1] >= 0);
	return;
      }	
      RIGHT_END(TMP) = INNER + kstep * STEP;
      SQUEEZED(TMP)  = - MAXINT / 2;
      ASSIGNMENT(TMP)= ASSIGN_IDX_MISMATCH;
      //  
      if (DEBUG_CUMSUM)
	printf("0. TMP=%d kstep=%d ass=%d\n", TMP, kstep,  ASSIGNMENT(TMP));//
      cumsum[TMP++] = WEIGHT(kstep - 1);
	assert(cumsum[TMP-1] >= 0);
    } else {
      int k;
      for (k=-IDX_STEPS; k<kstep; k++) {
	RIGHT_END(TMP) = INNER + (k + IDX_STEPS) * STEP;
	SQUEEZED(TMP)  = dM1;
	ASSIGNMENT(TMP) = k == -IDX_STEPS ? ASSIGN_IDX_INNER : k;
	//
	if (DEBUG_CUMSUM)
	  printf("1. TMP=%d k=%d ass=%d\n", TMP, k,  ASSIGNMENT(TMP));//
	cumsum[TMP++]  = WEIGHT(k);
	assert(cumsum[TMP-1] >= 0);
      }
    }
    zneu[dM1] = INNER + STEP * kstep; 
    start_cumul = TMP;
  } else {
    kstep = 0;
    //  printf("AA d=%d %d %f %f\n", d, TMP, zneu[dM1], INNER);
  }
  
  for ( ; d<=dim; d++) {
    dM1 = d - 1;
    if (zneu[d] == zneu[dM1]) continue;
    //printf("%f == %f %d\n",zneu[d],  zneu[dM1], zneu[d] == zneu[dM1] );
   // printf("A d=%d %d\n", d, TMP);
 
    if (zneu[dM1] < INNER) {
      min = RIGHT_END(TMP) = zneu[d] <= INNER ? zneu[d] : INNER;
      SQUEEZED(TMP) = dM1;
      ASSIGNMENT(TMP) = ASSIGN_IDX_INNER;
      cumsum[TMP++] = INNER_CONST * 
	PoweredVolOfCube(ysort, zneu[dM1], min, INNER_POW, dim, dM1);
       // 
      if (DEBUG_CUMSUM)
	printf("2. TMP=%d kstep=%d ass=%d cum=%f inner_c=%f start=%f ende=%f pow=%f dim=%d squ=%d\n", TMP, kstep,  ASSIGNMENT(TMP), cumsum[TMP-1], //
	       INNER_CONST, zneu[dM1],
	       min,PoweredVolOfCube(ysort, zneu[dM1], min, INNER_POW, dim, dM1),
	       dim, dM1);// assert(false);
	assert(cumsum[TMP-1] >= 0);
      if (zneu[d] <= INNER) continue;
      zneu[dM1] = INNER;
    }

    if (zneu[dM1] < OUTER) {
      min = zneu[d] <= OUTER ? zneu[d] : OUTER;
      int k,
	steps = (min - zneu[dM1]) / STEP;
      double a = zneu[dM1];
      assert(steps <= NSTEP);
      
      //printf("d=%d min=%f %f %f\n", d, min, zneu[dM1], zneu[d]);
      
      //  
      assert(TMP <= NSTEP + 2 + dim);
      
      for (k=0; k<steps; k++, a += STEP, kstep++) {
	RIGHT_END(TMP) = a + STEP;
	SQUEEZED(TMP) = dM1;
	ASSIGNMENT(TMP) = kstep; // NICHT k !
	//   	double z; FCTN(&a, cov->sub[0], &z);
	//

	cumsum[TMP++] = 
	  VolumeOfCubeRing(ysort, a, a + STEP, dim, dM1) * VALUE(kstep);
	assert(cumsum[TMP-1] >= 0);	
  	//
	if (DEBUG_CUMSUM)
	  printf("3. TMP=%d kstep=%d ass=%d a=%4.2f st=%4.2f end=%4.2f val=%4.2e, tot=%4.4f\n", TMP, kstep,  ASSIGNMENT(TMP), a, STEP, RIGHT_END(TMP), //
		 VALUE(kstep), cumsum[TMP-1]);
	//	assert(dM1 != 1);
    }
      // 
      // printf("B %d %d\n", d, TMP);
      assert(TMP <= NSTEP + 2 + dim);
      
      if (kstep < NSTEP) {
	RIGHT_END(TMP) = min;
	SQUEEZED(TMP) = dM1;
	ASSIGNMENT(TMP) = kstep;
	assert(TMP > 0);
	// 
	if (DEBUG_CUMSUM)
	  printf("4. TMP=%d k=%d ass=%d\n", TMP, kstep, ASSIGNMENT(TMP));//
	cumsum[TMP++] = VolumeOfCubeRing(ysort, a, min, dim, dM1) * VALUE(kstep);
	assert(cumsum[TMP-1] >= 0);
	if (zneu[d] <= OUTER) continue;
      }
      
      // printf("outer=%f %f %f %e %d %d %d\n", OUTER, a,  cumsum[TMP-1], 
      //        OUTER - a, OUTER == min, dM1, kstep < NSTEP);
      //assert(false);
      // printf("C %d %d\n", d, TMP);
    }
      

    RIGHT_END(TMP) = zneu[d];
    SQUEEZED(TMP) = dM1;
    ASSIGNMENT(TMP) = ASSIGN_IDX_OUTER;
    //  
    if (DEBUG_CUMSUM) 
      printf("5. TMP=%d k=%d ass=%d %f p=%f right=%f\n", //
	     TMP, 9999, ASSIGNMENT(TMP), OUTER, OUTER_POW, zneu[d]);
    
    if (next->finiterange == true) {
      cumsum[TMP++] = 0.0;
      continue;
    } else if (OUTER_POW > 0) {
      //      printf("here %d C=%f Z=%f PC=%f e=%f %f\n", 
      ///	     IDX_OUTER, OUTER_CONST,zneu[d],
      //    OUTER_POW_CONST,
      //     ExpVolOfCube(OUTER, zneu[d], OUTER_POW, 
      //		  OUTER_POW_CONST, dim, dM1),
      //     OUTER_POW_CONST * ExpVolOfCube(OUTER, zneu[d], OUTER_POW, 
      //				    OUTER_POW_CONST, dim, dM1)	     
      //     );
      cumsum[TMP++] = OUTER_CONST * ExpVolOfCube(OUTER, zneu[d], OUTER_POW, 
						 OUTER_POW_CONST, dim, dM1);

     // 
      if (DEBUG_CUMSUM)
	printf("6. TMP=%d k=%d ass=%d oc=%f c=%f op=%f opc=%f dim=%d %d right=%f\n", TMP, 10001, ASSIGNMENT(TMP),OUTER_CONST, OUTER, OUTER_POW, //
	       OUTER_POW_CONST, dim, dM1, zneu[d]);
      //  printf("TMP=%d %f\n", TMP, cumsum[TMP-1]);
      // if (cumsum[TMP-1] < 0) {
      //t i;
      //fori=1; i<TMP; i++) {
      //	  if (cumsum[i] < 0)
      //	    printf("i=%d %f %f\n", i, cumsum[i], cumsum[i - 1]);
      //	}
      //	PMI(cov);
      //      }

      assert(cumsum[TMP-1] >= 0);
      assert(R_FINITE(cumsum[TMP-1]));

    } else {
      //printf("here XXX\n");
      cumsum[TMP++]= OUTER_CONST * 
	PoweredVolOfCube(ysort, OUTER, zneu[d], OUTER_POW, dim, dM1);
      assert(cumsum[TMP-1] >= 0);
    }

  }


  //PMI(cov);

  // make cummulative weights
  // 
  //d=0; printf("cumsum %d %e;   %d %e\n", d, cumsum[d], start_cumul-1, cumsum[start_cumul-1]);
  for (d=start_cumul; d<TMP; d++) {
    // printf("cumsum %d %e\n", d, cumsum[d]);
    cumsum[d] += cumsum[d-1];
    //     printf("cumsum %d c=%e %e %d %d end=%f inner=%f\n", d, cumsum[d], cumsum[d]-cumsum[d-1], ASSIGNMENT(d), ASSIGNMENT(d-1), RIGHT_END(d), INNER);
    assert(ASSIGNMENT(d)==ASSIGN_IDX_OUTER || ASSIGNMENT(d) >= ASSIGNMENT(d-1));
  }

  //APMI(cov);
}


// (unnormierte) VERTEILUNGSFUNKTION F(x,x,...,x) ist 
// c * (exp(-O^p - exp(-a x^p)) fuer grosse x, also integriert ueber 
// den Raum, gesehen als Funktion einer Variable x; O=OUTER
// Somit ist die Dichte gleich c * s * p * x^{p-d} * e^{-s x^p} / b
// wobei b die Oberflaeche des Wuerfels mit Kantenlaenge 2 bezeichnet
//
// weiter in  ExpVolOfCube
#define OUTER_DENSITY(v, x)					\
  if (OUTER_POW > 0) {						\
    double pow_x = OUTER_POW_CONST * pow(x, OUTER_POW);		\
    v = OUTER_CONST * OUTER_POW * pow_x * intpow(x, -dim) *	\
      exp(-pow_x) / SURFACE(dim);				\
  } else {							\
    /* DICHTE ist c*x^p fuer grosse x */			\
    v = OUTER_CONST * pow(x, OUTER_POW);			\
  }


void evaluate_rectangular(double *x, cov_model *cov, double *v) {  
  // *x wegen searchInverse !!
  RECTANGULAR_PARAMETERS;

  if (*x < 0.0) BUG; //crash();
  assert(*x >= 0);

  if (*x <= INNER) { // for technical reasons <= not <; see Dinverse
    // DICHTE ist c*x^p nahe Ursprung
    //  printf("inner ");
    *v = INNER_CONST * pow(*x, INNER_POW);
    //  printf("x=%f %f %f\n", *x, INNER, *v);
    return;
  } else if (*x >= OUTER) {
    //printf("outer ");
    if (next->finiterange == true) {
      *v = 0.0;
      return;
    }
    OUTER_DENSITY(*v, *x);
    return;					
  } else {
    //   printf("k=%d %f %f ", (int) ((*x - INNER) / STEP),
    //	   VALUE(0),
    //	   VALUE( (int) ((*x - INNER) / STEP)));
   *v = VALUE( (int) ((*x - INNER) / STEP));
    return;
  }
}

void rectangularD(double *x, cov_model *cov, double *v) {
  bool onesided = P0INT(RECT_ONESIDED); 
  if (onesided && *x <= 0) {
    *v = 0;
    return;
  }
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
  RECTANGULAR_PARAMETERS;
  int i;
   
  double max = RF_NEGINF;

  RECTFOR {double y=fabs(x[i]); if (y > max) max = y;}

  evaluate_rectangular(&max, cov, v);
  
  // printf("max %f %f\n", max, *v);
  // printf("rectD %f %f %d %f\n", max, *v, P0INT(RECT_NORMED), OUTER_CUM);

  if (P0INT(RECT_NORMED)) *v /= OUTER_CUM; 
  if (onesided) *v *= 2.0; // R^1 !!
}


void rectangularDlog(double *x, cov_model *cov, double *v, double *Sign) {
  rectangularD(x, cov, v);
  *v = log(*v);
  *Sign = 1.0;
}


void rectangularDinverse(double *V, cov_model *cov, double *left,
			 double *right) {
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
  RECTANGULAR_PARAMETERS;
  double  er, outer,
    x = RF_NAN,
    v = *V; 
  int i;
  bool onesided = P0INT(RECT_ONESIDED); 

  if (P0INT(RECT_NORMED)) v *= OUTER_CUM; // unnormed v
  if (onesided) v /= 2.0;
  
  //  printf("V=%e dim=%d\n", *V, dim);
    // assert(*V >= 0.0 && *V <= 1.0);
  
  if (*V <= 0.0) {
    RECTFOR {
      left[i] = RF_NEGINF;
      right[i] = RF_INF;
      //printf("i=%d %f %f\n", i, left[i], right[i]);
     }
     return;
  }

  
  if (!next->finiterange == true && OUTER_POW > 1) {
    // local maximum of density function exists for outer_pow > 1
    outer = pow( (OUTER_POW -1) / (OUTER_POW * OUTER_POW_CONST), 1 / OUTER_POW);
    if (outer < OUTER) outer = OUTER;    
  } else outer = OUTER;

  evaluate_rectangular(&outer, cov, &er); 

  //   printf("outer =%f %f %f %f\n", outer, OUTER, er, v);

  if (er > v) { // inverse is in the tail
    // first guess:
    double inverse,
      eps = 0.01;
     
    if (OUTER_POW > 0) {

      //      printf("here\n");
      // rough evaluation
      inverse = pow(- log(v / (OUTER_CONST * OUTER_POW)) / 
		    OUTER_POW_CONST, 1.0 / OUTER_POW);
      if (inverse <= outer) inverse = 2 * outer;
      x = searchInverse(evaluate_rectangular, cov, inverse, outer, v, eps);
    } else 
      //      printf("hereB\n");
      x = pow(OUTER_CONST / v, 1 / OUTER_POW);

  } else {
    //printf("hereC\n");
      int nsteps = (OUTER - INNER) / STEP;
    for (i=nsteps - 1; i>=0; i--) 
      if (VALUE(i) >= v) {
	x = INNER + (i+1) * STEP;
	break;
      }
    if (i<0) { // all values had been smaller, so the 
      // inverse is in the inner part
      //     printf("hereD\n");
     evaluate_rectangular(&(INNER), cov, &er);
      if (er >= v) x = INNER;
      else if (INNER_POW == 0) x = 0.0; 
      else if (INNER_POW < 0.0) {
	x = pow(v / INNER_CONST, 1.0 / INNER_POW);	
      } else BUG;
    }
  }
  
  RECTFOR {
    //printf("rect i=%d x=%f\n", x);
    left[i] = onesided ? 0.0 : -x;
    right[i] = x;
  }
}


void rectangularP(double VARIABLE_IS_NOT_USED *x, 
		  cov_model VARIABLE_IS_NOT_USED *cov, 
		  double VARIABLE_IS_NOT_USED *v) {
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
 // RECTANGULAR_PARAMETERS;

  // not programmed yet
  BUG;

}


void rectangularP2sided(double *x, double *y, cov_model *cov, double *v) {
  // distribution functions
  bool onesided = P0INT(RECT_ONESIDED); 
  int d;
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
  RECTANGULAR_PARAMETERS;
  if (x != NULL) BUG;
  if (onesided && *y <= 0) {
    *v = 0;
    return;
  } 
  for (d=0; d<dim; d++) {// To Do ?? Stimmt das? Oder Masse der dimension d-k ??
    assert(y[d] >= 0.0);
    if (y[d] == 0.0) {
      *v = 0.0;
      return;
    }
  }
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
     
  //  CumSum(y, false, cov, rect->tmp_weight);
  //  double tw = TMP_WEIGHT[TMP-1];

  CumSum(y, true, cov, rect->tmp_weight);

  
 
  *v = TMP_WEIGHT[TMP-1];
  // print("tw=%f %f y=%f %f %f\n", tw, *v, y[0], y[1], y[2]); 
  //  if(*v != tw) APMI(cov);
  //  assert(*v == tw);

  if (P0INT(RECT_NORMED)) *v /= OUTER_CUM;

  //  printf("p2sided %f %d outercum=%f, %d %f\n",
  //	 *y, P0INT(RECT_NORMED), OUTER_CUM, 
  //	 TMP, TMP_WEIGHT[TMP-1]);
  //  int i; for (i=0; i<TMP; i++) printf("ii=%d %f\n", i, TMP_WEIGHT[i]);

  assert(R_FINITE(*v));

}


void rectangularQ(double  VARIABLE_IS_NOT_USED *x, 
		  cov_model  VARIABLE_IS_NOT_USED *cov, 
		  double  VARIABLE_IS_NOT_USED *v) {
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
   error("rectangularQ not programmed yet");
}


#define ACCEPT_REJECT							\
  if (!P0INT(RECT_APPROX)) {				\
    int ni = dim,							\
      quot = dim + 1,							\
      bytes =  sizeof(double) * dim;					\
    double newquot, approx, truevalue,					\
      max = RF_NEGINF;							\
    RECTFOR {double z=fabs(v[i]); if (z > max) max = z;}		\
    evaluate_rectangular(&max, cov, &approx);				\
    FCTN(v, next, &truevalue);						\
    newquot = truevalue / approx;					\
    assert(quot < cov->qlen + 2);					\
    if (isMonotone(next->monotone)) { /* rejection sampling */		\
      assert(approx >= truevalue);					\
      cov->q[ni] = 0;							\
       if (UNIFORM_RANDOM >= newquot) {					\
	/* printf("."); */						\
	RECT_R;								\
      } /* else printf(":"); */						\
    } else { /* MCMC */							\
      if (R_FINITE(cov->q[ni])) {					\
	cov->q[ni] = cov->q[ni] - 1.0;					\
	if (UNIFORM_RANDOM * cov->q[quot] >= newquot) { /* keep old one */ \
	  memcpy(v, cov->q, bytes);					\
	} else { /* accept new one*/					\
	  cov->q[quot] = newquot; /* pi / Q */				\
	  memcpy(cov->q, v, bytes);					\
	}								\
      } else { /* !R_FINITE(cov->q[ni]), i.e. starting value */		\
	cov->q[ni] = P0INT(RECT_MCMC_N) - 1.0;		\
	cov->q[quot] = newquot; /* pi / Q */				\
	memcpy(cov->q, v, bytes);					\
      }									\
    }									\
    if (cov->q[ni] > 0) {						\
      RECT_R;								\
    } else {								\
      cov->q[ni] = P0INT(RECT_MCMC_N);			\
      if (final) {							\
	cov->total_n++;							\
	cov->total_sum += truevalue;					\
      }									\
      return;								\
    }									\
  } else {								\
    if (final) {							\
      double approx,							\
	max = RF_NEGINF;						\
      RECTFOR {double z=fabs(v[i]); if (z > max) max = z;}		\
      evaluate_rectangular(&max, cov, &approx);				\
      cov->total_n++;							\
      cov->total_sum += approx;						\
    }									\
  }



void rectangularR(double *x, cov_model *cov, double *v) { 
  if (x != NULL) error("put 'flat = false'");
  RECTANGULAR_PARAMETERS;
  bool final = true;
  
  int i = CeilIndex(UNIFORM_RANDOM * OUTER_CUM, rect->weight, NSTEP + 2);
  //int i = searchFirstGreater(rect->weight, NSTEP + 2, UNIFORM_RANDOM  * OUTER_CUM);

  //printf("i=%d dim=%d %d %d %s\n", i, dim, IDX_INNER, IDX_OUTER, NICK(cov->sub[0]));
 
  if (i == IDX_INNER) {
    RandomPointOnCubeSurface(pow(UNIFORM_RANDOM, 1.0 / (INNER_POW+dim)) * INNER,
			     dim, v);
  } else if (i == IDX_OUTER) {
     double u;
    assert(next->finiterange == false);
    if (OUTER_POW > 0) // sampling by inversion formula
      u = pow(pow(OUTER, OUTER_POW) - log(UNIFORM_RANDOM) / OUTER_POW_CONST,
	      1 / OUTER_POW);
    else u = pow(UNIFORM_RANDOM, 1.0 / (OUTER_POW + dim)) * OUTER; // p+d < 0 !
    RandomPointOnCubeSurface(u, dim, v);
  } else { // i >= 0
   i -= IDX_STEPS;
    double start = INNER + i * STEP;
    RandomPointOnCubeRing(start, start + STEP, dim, v); 
  }
  if (P0INT(RECT_ONESIDED)) *v = fabs(*v);
#define RECT_R rectangularR(x, cov, v)
 
  ACCEPT_REJECT;
  
}

static int zi=0, zs=0, zo=0;

void rectangularR2sided(double *x, double *y, cov_model *cov, double *v) { 
  if (x != NULL) 
    error("2-sided distribution function for rectangular not programmed yet");
  RECTANGULAR_PARAMETERS;
  int d, sel, red_dim, i,
    *I = rect->i;
  double *u, start, end,
    *ysort = rect->ysort;    
 
  CumSum(y, false, cov, rect->tmp_weight); // ACHTUNG! reck->z gesetzt
  //assert(rect->tmp_weight[TMP-1] >= 2.1 && rect->tmp_weight[TMP-1] <= 2.2);
  double random = UNIFORM_RANDOM * rect->tmp_weight[TMP-1];
  bool final = SQUEEZED(TMP-1) == 0 && 
    (!P0INT(RECT_APPROX) || !next->deterministic); // ohne Einschraenkung
  assert(!final); // nur zur Kontrolle

  sel = CeilIndex(random, rect->tmp_weight, TMP);
  //printf("r=%f inner=%f o.tmp=%f lst.step.tmp=%f\n", random,  rect->tmp_weight[IDX_INNER],  rect->tmp_weight[TMP-1],rect->tmp_weight[TMP-1]);

  red_dim = dim - SQUEEZED(sel);
  if (red_dim <= 0) BUG;
  start = sel > 0 ? RIGHT_END(sel - 1) : 0;
  end = RIGHT_END(sel);
  assert(start < end);

  // PMI(cov, -1); 
 
  // printf("%f dim=%d squ=%d sel=%d, %f %f %d\n",
  //	 y[0], dim, SQUEEZED(sel), sel, random, rect->tmp_weight[TMP-1], TMP);
 
 
  u = rect->tmp_weight; // nutzen des Platzes...
  //printf("assess %d %d %d\n", ASSIGNMENT(sel), sel, TMP);
  if (ASSIGNMENT(sel) == ASSIGN_IDX_INNER) {
    //   RandomPointOnCubeSurface(pow(UNIFORM_RANDOM, 1.0 / (INNER_POW+dim))
    //                        * INNER, dim, v);
    zi ++;
    double // a + b x^p probab distr on [s, e], i.e. b = 1/(e^p-s^p),
      p = INNER_POW + red_dim,
      sp = pow(start, p),
      ep = pow(end, p),
      binv = ep - sp,
      a = - sp / binv,
      r = pow((UNIFORM_RANDOM - a) * binv, 1 / p);
    //
    //   printf("inner %d %d ap=%f r=%f inner=%f w=%f tmp=%f, %f\n", dim, SQUEEZED(sel), ap, r, INNER, INNER_CUM, TMP_WEIGHT[IDX_INNER], TMP_WEIGHT[IDX_INNER+1]);

    RandomPointOnCubeSurface(r, red_dim, u);
  } else if (ASSIGNMENT(sel) == ASSIGN_IDX_OUTER) {
    zo++;
    // printf("outer\n");
    double w;
    assert(next->finiterange == false);
    if (OUTER_POW > 0) {// sampling by inversion formula
      double op = pow(OUTER, OUTER_POW),
	factor = 1.0 - exp(- OUTER_POW_CONST * (pow(end, OUTER_POW) - op));
      w = pow(op - log(1 - UNIFORM_RANDOM * factor) / OUTER_POW_CONST,
	      1 / OUTER_POW);
    } else w =pow(1 - UNIFORM_RANDOM *(1 - pow(end/OUTER, OUTER_POW + red_dim)),
		  1.0 / (OUTER_POW + red_dim)) * OUTER; 
		 // p+d < 0 ! 
    //printf("hXXere\n");
    RandomPointOnCubeSurface(w, red_dim, u);
  } else { // i \in steps
    zs++;
    //  printf("steps \n");
    assert(ASSIGNMENT(sel) >= 0 && ASSIGNMENT(sel) < NSTEP);
    // printf("hAere\n");
   
    // APMI(cov->calling);

    RandomPointOnCubeRing(start, end, red_dim, u); 
  }
  // printf("zi=%d s=%d o=%d tot=%d\n", zi, zs, zo, zi+zs+zo);
 
  // up to now simulation on cubes on the non-squeezed dimensions
  // Now the coordinates have to put on the right place and
  // the remaining coordinates have to be filled with uniformly scattered
  // random variables on [-ysort[d], ysort[d]
  // I[.] and ysort[.] have been set by CumSum
  assert(I[0] == 0.0);
  for (d=1; d<=SQUEEZED(sel); d++) {
    v[I[d] - 1] = (2.0 * UNIFORM_RANDOM - 1.0) * ysort[d];//I[d] indiziert ab 1
  //                                                    da ysort[0]=0 per def.
  }
  for (d=SQUEEZED(sel); d<dim; d++) {
    v[I[d+1] - 1] = u[d - SQUEEZED(sel)];
  }

  if (P0INT(RECT_ONESIDED)) *v = fabs(*v);

#undef RECT_R
#define RECT_R rectangularR2sided(x, y, cov, v)	
  ACCEPT_REJECT;

  //printf("u=%f %f %f\n", u[0], u[1], u[2]);
}


int GetMajorant(cov_model *cov) {
  RECTANGULAR_PARAMETERS;
  double v, m, delta, old_ratio,
    starting_point = 1.0,
    min_steps = P0(RECT_MINSTEPLENGTH), // middle
    safety = P0(RECT_SAFETY),
    safetyP1 = safety + 1.0,
    EMsafety = 1.0 - safety,
    dimsafety = dim * safety,
    inner_min = P0(RECT_INNERMIN),
    outer_max = P0(RECT_OUTERMAX),
    inner_factor = 0.5,
    outer_factor = 1.3
    ;
  
  assert(!PisNULL(RECT_MAXSTEPS));

  int i, j, 
    max_steps = P0INT(RECT_MAXSTEPS),
    parts = P0INT(RECT_PARTS),
    max_iterations =  P0INT(RECT_MAXIT);

  assert(next->taylor[0][TaylorConst] > 0.0);
  INNER_CONST = next->taylor[0][TaylorConst] * safetyP1;
  INNER_POW = next->taylor[0][TaylorPow] * EMsafety - dimsafety;
  

  // INNER
  i = 0;
  INNER = starting_point;
  FCTN(&(INNER), next, &v);
  m = INNER_CONST * pow(INNER, INNER_POW);
  //   PMI(next, -1);
  // printf("inner=%e m=%e v=%e inner_c=%4.3f (%4.3f) inner_p=%4.3f (%4.3f)\n", 
  //	 INNER, m, v, INNER_CONST, next->taylor[0][TaylorConst],
  //	 INNER_POW, next->taylor[0][TaylorPow]);
  while (m < v && INNER >= inner_min) {
    INNER *= inner_factor;
    INNER_CONST *= safetyP1;
    INNER_POW = INNER_POW * EMsafety - dimsafety;
    FCTN(&(INNER), next, &v);
    m = INNER_CONST * pow(INNER, INNER_POW);
    // printf("inner=%e m=%e v=%e i=%d, maxit=%d\n", 
    //	INNER, m, v, i, max_iterations);
  }
  if (INNER < inner_min && m < v) {
    //PMI(cov);
    SERR("no majorant found close to the origin");
  }

 
  i = 0;

  while (true) {
    old_ratio = m / v;
    delta = INNER / parts;
    double x = INNER - delta;

    //printf("delta %f x=%f i=%d\n", delta, x, i);

    for (j=1; j<parts; j++, x-=delta) {
      //      printf("j=%d parts=%d max_it=%d\n", j, parts, max_iterations);

      FCTN(&x, next, &v);
      m = INNER_CONST * pow(x, INNER_POW);
      double ratio = m / v;

      //printf("inner j=%d (%d) x=%4.3f m=%4.3f ratio=%4.3f old=%4.3f\n", 
      // j, max_iterations, x, m, ratio, old_ratio);

      if (ratio < old_ratio) break;
      old_ratio = ratio; 
    }

    if (j == parts) break;
    INNER *= inner_factor;
    if (INNER > x) INNER = x;
    if (++i > max_iterations) 
      SERR1("max. iterations of %d performed without success. Increase the value of 'maxit'", max_iterations);
  }

  //printf("inner j=%d (%d) inner=%4.3f m=%4.3f\n", 
      //  j, max_iterations, INNER, m);


  // OUTER
  if (next->tail[0][TaylorExpPow] <= 0.0) {
    // polynomial tail
    OUTER_POW = next->tail[0][TaylorPow] * safetyP1 + dimsafety;
    OUTER_POW_CONST = RF_NAN;
  } else { // exponential tail
    assert(next->tail[0][TaylorExpConst] > 0.0);
    OUTER_POW = next->tail[0][TaylorExpPow] / safetyP1;
    OUTER_POW_CONST = next->tail[0][TaylorExpConst] / safetyP1;
  }
  OUTER_CONST = next->tail[0][TaylorConst] * safetyP1;
 
  if (next->finiterange == true) {
    // printf("xx %s %d\n", CovList[next->sub[0]->nr].nick, CovList[next->sub[0]->nr].finiterange == 1);
    
    // APMI(cov);
    INVERSE(ZERO, next, &(OUTER));
    if (INNER >= OUTER) OUTER = INNER * (1.0 + 1e-5);
  } else { // ! finiterange
    assert(next->tail[0][TaylorConst] > 0.0);
    i = 0;
    OUTER = starting_point * safetyP1;
    FCTN(&(OUTER), next, &v);
    
    //    m = OUTER_CONST * (OUTER_POW <= 0 ? pow(OUTER, OUTER_POW) 
    //		       : exp(- OUTER_POW_CONST *  pow(OUTER, OUTER_POW)));
    OUTER_DENSITY(m, OUTER);
    while (m < v && OUTER < outer_max) {
      if (++i > max_iterations) 
	SERR("no majorant found. Function does allow for a majorant or increase 'maxit'");
      OUTER_CONST *= safetyP1;
      OUTER_POW /= safetyP1;
      OUTER *= outer_factor;
      FCTN(&(OUTER), next, &v);
      OUTER_DENSITY(m, OUTER);
      // printf("i =%d  outer %f m=%e v=%e\n", i, OUTER, m, v); assert(false);
  }
    if (OUTER > outer_max && m < v)
      SERR("No majorant found close for large arguments");
    
    //  printf("i =%d  outer %f m=%e v=%e\n", i, OUTER, m, v); //assert(false);

    //printf("i =%d  outer %f %f m=%f v=%f\n", 
    //i, OUTER, 
    //	  OUTER_CONST * (OUTER_POW <= 0 ? pow(OUTER, OUTER_POW) 
    //			 : exp(- OUTER_POW_CONST *  pow(OUTER, OUTER_POW))),
    //	  m , v); //assert(false);

    i=0;
    while (true) {
      old_ratio = m / v;
      delta = OUTER / parts;
      double x = OUTER + delta;
      for (j=1; j<parts; j++, x+=delta) {
	FCTN(&x, next, &v);
	OUTER_DENSITY(m, x);
	double ratio = m / v;
	if (ratio < old_ratio) break;
	old_ratio = ratio; 
      }
      if (j == parts) break;
      OUTER *= outer_factor;
      if (OUTER < x) OUTER = x;
      if (++i > max_iterations) BUG;
    }
  }

  //  printf("inner outer %f %f\n", INNER, OUTER);// assert(false);
  // STEPS  
  
  // todo variable Schrittweiten waeren wegen des steilen
  // Anstiegs am Pol besser
  STEP = (OUTER - INNER) / max_steps;
  if (STEP < min_steps) {
    NSTEP = (int) ((OUTER - INNER) / min_steps);
    STEP = (OUTER - INNER) / NSTEP;
  } else NSTEP = max_steps;

  // PMI(cov, "getm");

  return NOERROR; 
}

int check_rectangular(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err,
    dim = cov->xdimown;

  ROLE_ASSERT(ROLE_DISTR);

  kdefault(cov, RECT_SAFETY, GLOBAL.distr.safety);
  kdefault(cov, RECT_MINSTEPLENGTH, GLOBAL.distr.minsteplen);
  kdefault(cov, RECT_MAXSTEPS, GLOBAL.distr.maxsteps);
  kdefault(cov, RECT_PARTS, GLOBAL.distr.parts);
  kdefault(cov, RECT_MAXIT, GLOBAL.distr.maxit);
  kdefault(cov, RECT_INNERMIN, GLOBAL.distr.innermin);
  kdefault(cov, RECT_OUTERMAX, GLOBAL.distr.outermax);
  kdefault(cov, RECT_MCMC_N, GLOBAL.distr.mcmc_n);
  kdefault(cov, RECT_NORMED, true);
  kdefault(cov, RECT_APPROX, true);
  kdefault(cov, RECT_ONESIDED, false);
  //  int onesided = P0INT(RECT_ONESIDED);
  
  cov->q = (double *) CALLOC(dim + 2, sizeof(double));
  cov->qlen = 1;
  cov->q[dim] = RF_NAN;

  if ((err = CHECK(next, dim, dim, ShapeType, XONLY, 
		   dim == 1 && P0INT(RECT_ONESIDED) ? CARTESIAN_COORD 
		   : ISOTROPIC, 
		   SCALAR, ROLE_BASE)) != NOERROR) {
    return err;
  }

  if (!next->deterministic) {
    //APMI(cov);
    SERR("currently, only deterministic submodels are allowed"); // to do
  }
  
  // if (!isMonotone(next->monotone)) SERR("only monotone submodels are allowed");

  if (next->taylorN <= 0 || next->tailN <= 0) {
    //PMI(next);
    SERR1("'%s' does not have integrability information", NICK(next));
  }
  assert(R_FINITE(next->taylor[0][TaylorPow]));
  if (next->taylor[0][TaylorPow] <= -dim)
    SERR1("pole of '%s' not integrable", NICK(next));
  assert(R_FINITE(next->tail[0][TaylorPow]));

  //  PMI(next);

  if (next->tail[0][TaylorPow] >= -dim && next->tail[0][TaylorExpPow] == 0.0 &&
      next->tail[0][TaylorConst] != 0.0) 
    SERR1("tail of '%s' not integrable", NICK(next));
  assert(R_FINITE(next->tail[0][TaylorConst]));

  //PMI(next);
  if (next->taylor[0][TaylorConst] == 0.0) 
    SERR1("'%s' seems to be a trivial shape function", NICK(next));

  //  APMI(cov);
  
  if (cov->xdimprev != dim || dim != cov->tsdim) return ERRORDIM;

  cov->vdim2[0] = cov->tsdim;
  cov->vdim2[1] = 1;

  //PMI(cov);

  return NOERROR;
}



int init_rectangular(cov_model *cov, storage VARIABLE_IS_NOT_USED *s){
  int  err;

  // printf("init rectangular\n");
 
  if (cov->Srect != NULL) RECT_DELETE(&(cov->Srect)); // init_rectangular 
  // kann ein 2. Mal aufgerufen werden
  cov->Srect = (rect_storage *) MALLOC(sizeof(rect_storage));
  RECT_NULL(cov->Srect);

  RECTANGULAR_PARAMETERS; // muss nach obigen stehen und vor allem anderen

  if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR) {
    //PMI(cov);
    return err;
  }

  if ((err = GetMajorant(cov)) != NOERROR) return err; // muss genau nach MALLOC oben und MALLOC unten stehen
  if (INNER >= OUTER) BUG;

  //   APMI(cov);
   
  int d,
    abschnitte = NSTEP + 2, 
    tmp_n = NSTEP + 2 + dim;
  double x;

  if ((rect->value = (double*) MALLOC(sizeof(double) * abschnitte)) == NULL ||
      (rect->weight = (double*) MALLOC(sizeof(double) * abschnitte)) == NULL||
      (rect->tmp_weight = (double*) CALLOC(tmp_n, sizeof(double))) == NULL ||
      (rect->right_endpoint = (double*) MALLOC(sizeof(double) * tmp_n)) ==NULL||
      (rect->ysort = (double*) MALLOC(sizeof(double) * (dim + 1))) == NULL ||
      (rect->z = (double*) MALLOC(sizeof(double) * (dim + 1))) == NULL ||//dummy
      (rect->squeezed_dim = (int*) MALLOC(sizeof(int) * tmp_n)) == NULL ||
      (rect->assign = (int*) MALLOC(sizeof(int) * tmp_n)) == NULL ||
      (rect->i = (int*) MALLOC(sizeof(int) * (dim + 1))) == NULL)
   return ERRORMEMORYALLOCATION;

 
  x = INNER;
  for (d=0; d<NSTEP; d++, x+=STEP) {
    FCTN(&x, next, &(VALUE(d)));
    //printf("x=%f %f\n", x, VALUE(d));
    assert(VALUE(d) >= 0);
  }
  rect->value[IDX_INNER] = rect->value[IDX_OUTER] = RF_NAN;

  for(d=0; d<dim; d++) rect->tmp_weight[d] = RF_INF;
  
  CumSum(rect->tmp_weight, false, cov, rect->weight);

  cov->mpp.M[0] = cov->mpp.Mplus[0] = P0INT(RECT_NORMED) ? 1.0 : OUTER_CUM;
  assert(R_FINITE(OUTER_CUM));
  if (cov->mpp.moments > 0) {    
    cov->mpp.M[1] = next->mpp.M[1];
    cov->mpp.Mplus[1] = next->mpp.Mplus[1]; 
    if (!R_FINITE(cov->mpp.M[1])){
      //APMI(cov);
      BUG;
    }
  }

  // normed against cov->sub[0]
  cov->mpp.maxheight = OUTER_CUM;

  //PMI(cov);

  if (false) {
    double mini=1e-6, u,w,t, tot_u=0.0, tot_t=0.0, e,
      xx = mini;
    // 
    mini = STEP; xx=INNER-mini;
    double end;
    int k,
      n = NSTEP
      ;
    n = 1 * NSTEP;
    assert(mini >= 0.0);
    assert(xx >= 0.0);
    for (k=0, end=INNER; k<=n+1; k++, end+=STEP) {
      if (k==n+1) end *= 5;
      u = t = 0.0;
      while (xx < end) {
	FCTN(&xx, next, &w);
	double y = xx + 1e-10;
	assert(y >= 0.0);
	evaluate_rectangular(&y, cov, &e);
	//printf("%f %e %e\n", xx, w, e);
	u += w * 4 * M_PI / 3 * (powl(xx + mini, 3) - powl(xx, 3));
	t += w * (powl(2 * (xx + mini), 3) - powl(2 * xx, 3));
	xx += mini;
      }
      tot_u += u;
      tot_t += t;
      //      printf("%d %e %e tot=%f %f\n", k-1, u, t, tot_u, tot_t);
    }
  }

  //  printf("totalweight=%f\n", OUTER_CUM);
  assert(OUTER_CUM >= 1.0 || next->nr != STROKORB_MONO); 
  //
  //  printf("init rect %d %d\n", NSTEP, TMP);
  assert(NSTEP == TMP - 2);
  
  //  A  PMI(cov, -1);
  //assert(false);
  cov->total_n = 0;
  cov->total_sum = 0.0;

   return NOERROR;
}


void do_rectangular(cov_model *cov, double *v){
  RECTANGULAR_PARAMETERS;
  //
  storage s;
  STORAGE_NULL(&s);  
  DO_PARAM_MODELS;
  if (!next->deterministic) {
    DO(next, &s);
    cov->initialised = false;
    INIT(cov, cov->mpp.moments, &s);
  }
  rectangularR(NULL, cov, v);
}

void range_rectangular(cov_model *cov, range_type *range){
  range->min[RECT_SAFETY] = 0;
  range->max[RECT_SAFETY] = RF_INF;
  range->pmin[RECT_SAFETY] = 0.001;
  range->pmax[RECT_SAFETY] = 0.1;
  range->openmin[RECT_SAFETY] = true;
  range->openmax[RECT_SAFETY] = true;
  
  range->min[RECT_MINSTEPLENGTH] = 0;
  range->max[RECT_MINSTEPLENGTH] = RF_INF;
  range->pmin[RECT_MINSTEPLENGTH] = 0.01;
  range->pmax[RECT_MINSTEPLENGTH] = 10;
  range->openmin[RECT_MINSTEPLENGTH] = false;
  range->openmax[RECT_MINSTEPLENGTH] = true;
  
  range->min[RECT_MAXSTEPS] = 1;
  range->max[RECT_MAXSTEPS] = RF_INF;
  range->pmin[RECT_MAXSTEPS] = 2;
  range->pmax[RECT_MAXSTEPS] = 10000;
  range->openmin[RECT_MAXSTEPS] = false;
  range->openmax[RECT_MAXSTEPS] = true;
  
  range->min[RECT_PARTS] = 2;
  range->max[RECT_PARTS] = RF_INF;
  range->pmin[RECT_PARTS] = 3;
  range->pmax[RECT_PARTS] = 10;
  range->openmin[RECT_PARTS] = false;
  range->openmax[RECT_PARTS] = true;
  
  range->min[RECT_MAXIT] = 1;
  range->max[RECT_MAXIT] = RF_INF;
  range->pmin[RECT_MAXIT] = 2;
  range->pmax[RECT_MAXIT] = 10;
  range->openmin[RECT_MAXIT] = false;
  range->openmax[RECT_MAXIT] = true;
  
  range->min[RECT_INNERMIN] = 0;
  range->max[RECT_INNERMIN] = RF_INF;
  range->pmin[RECT_INNERMIN] = 1e-100;
  range->pmax[RECT_INNERMIN] = 1;
  range->openmin[RECT_INNERMIN] = true;
  range->openmax[RECT_INNERMIN] = true;
  
  range->min[RECT_OUTERMAX] = 0;
  range->max[RECT_OUTERMAX] = RF_INF;
  range->pmin[RECT_OUTERMAX] = 1;
  range->pmax[RECT_OUTERMAX] = 1e10;
  range->openmin[RECT_OUTERMAX] = true;
  range->openmax[RECT_OUTERMAX] = true;
  
  range->min[RECT_MCMC_N] = 1;
  range->max[RECT_MCMC_N] = MAXINT;
  range->pmin[RECT_MCMC_N] = 0;
  range->pmax[RECT_MCMC_N] = 1000;
  range->openmin[RECT_MCMC_N] = false;
  range->openmax[RECT_MCMC_N] = true;

  range->min[RECT_NORMED] = 0;
  range->max[RECT_NORMED] = 1;
  range->pmin[RECT_NORMED] = 0;
  range->pmax[RECT_NORMED] = 1;
  range->openmin[RECT_NORMED] = false;
  range->openmax[RECT_NORMED] = false;

  range->min[RECT_APPROX] = 0;
  range->max[RECT_APPROX] = 1;
  range->pmin[RECT_APPROX] = 0;
  range->pmax[RECT_APPROX] = 1;
  range->openmin[RECT_APPROX] = false;
  range->openmax[RECT_APPROX] = false;

  range->pmin[RECT_ONESIDED] = range->min[RECT_ONESIDED] = 0;
  range->pmax[RECT_ONESIDED] = range->max[RECT_ONESIDED] = cov->tsdim == 1;
  range->openmin[RECT_ONESIDED] = false;
  range->openmax[RECT_ONESIDED] = false;
}




