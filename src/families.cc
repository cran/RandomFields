/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of multivariate distribution families

 Copyright (C) 2013 -- 2017 Martin Schlather

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




#include "def.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "families.h"
#include "operator.h"


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


void arcsqrtP(double *x, model  VARIABLE_IS_NOT_USED *cov, double *v) {
  double 
    scale = 4.0 * P0(ARCSQRT_SCALE),
      y = *x / scale;
 
  if (y <= PIHALF) *v = 0.0;
  else {
    *v = ATAN(SQRT( y / PIHALF  - 1.0)) / PIHALF; // without (5) !!
    // siehe tcf*KS-1.pdf
  }  
}

// --- not corrected yet:


void arcsqrtD(double *x, model  VARIABLE_IS_NOT_USED *cov, double *v) {
  double 
    scale = 4.0 * P0(ARCSQRT_SCALE),
    y = *x / scale;
 if (y <= PIHALF ) *v = 0.0;
 else *v = 1.0 / (M_PI * scale * y * SQRT(y / PIHALF - 1.0)); 
}



void arcsqrtDlog(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  arcsqrtD(x, cov, v);
  *v = LOG(*v);
}

void arcsqrtDinverse(double VARIABLE_IS_NOT_USED *v, model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *left, double VARIABLE_IS_NOT_USED *right) {
  if (v == NULL || *v <=0 ) {
    left[0] = 0;
    right[0] = RF_INF;
  } else ERR("Dinverse of arcsqrt unknown");
}

void arcsqrtQ(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  if (*x < 0 || *x > 1) {*v = RF_NA; return;}
 double 
   scale = 4.0 * P0(ARCSQRT_SCALE);
  double z = TAN(PIHALF * *x);
  *v = PIHALF * (z * z + 1.0) * scale;
  }		   

void arcsqrtR(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  if (x != NULL) *v = *x;
  else {
    double 
       u = UNIFORM_RANDOM;
    arcsqrtQ(&u, cov, v); 
    
    //  printf("scale %10g %10g \n", scale, *v);
  
    assert(*v >= PIHALF);

    //      *v = 0.1;
  
   }
}


int init_arcsqrt(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  if (cov->mpp.moments >= 0) {
    cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0; 
  }

  // maxheight based on a Gaussian shape function with the same
  // sd's but normed to 1 at the origin
  cov->mpp.maxheights[0] = RF_NA;
  cov->mpp.unnormedmass = RF_NA;
  RETURN_NOERROR;
}

void do_arcsqrt(model *cov, double *v){ 
  arcsqrtR(NULL, cov, v);
}


void range_arcsqrt(model  VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[ARCSQRT_SCALE] = 0.0;
  range->max[ARCSQRT_SCALE] = RF_INF;
  range->pmin[ARCSQRT_SCALE] = 0.0;
  range->pmax[ARCSQRT_SCALE] = 100000;
  range->openmin[ARCSQRT_SCALE] = false;
  range->openmax[ARCSQRT_SCALE] = true;
}

//////////////////////////////////////////////////////////////////////

void addVariable(char *name, double *x, int nrow, int ncol, SEXP env) {  
  SEXP Y;
  int j,
    size= nrow * ncol;
  if (ncol == 1) PROTECT(Y = allocVector(REALSXP, nrow));
  else PROTECT(Y = allocMatrix(REALSXP, nrow, ncol));
  double *y = REAL(Y);
  for (j=0; j<size; j++) y[j] = x[j];  
  defineVar(install(name), Y, env);
  UNPROTECT(1);
}


void addIntVariable(char *name, int *x, int nrow, int ncol, SEXP env) {  
  SEXP  Y;
  int j,
    size= nrow * ncol;
  if (ncol == 1) PROTECT(Y = allocVector(INTSXP, nrow));
  else PROTECT(Y = allocMatrix(INTSXP, nrow, ncol));
  int *y = INTEGER(Y);
  for (j=0; j<size; j++) y[j] = x[j];  
  defineVar(install(name), Y, env);
  UNPROTECT(1);
}



 
void evaluateDistr(model *cov, int which, double *Res) {
  SEXP  res,
    env = PENV(DISTR_ENV)->sexp;
  int size,
    nkappas = DefList[COVNR].kappas,
   i = DISTR_LAST + 1;
  //    i = nkappas;

  // PMI(cov);
  
  //  if (cov->ownkappanames != NULL) {
  //    while(cov->ownkappanames[--i] != NULL) {
  //     addVariable(cov->ownkappanames[i], P(i), cov->nrow[i], cov->ncol[i], 
  //		  env);		
  //  }  
  // }

  
 if (cov->ownkappanames != NULL) {
    while(i < nkappas && cov->ownkappanames[i] != NULL) {
      addVariable(cov->ownkappanames[i], P(i), cov->nrow[i], cov->ncol[i], 
		  env);	
      i++;
    }  
  }

  //  PMI(cov);
  //  print("which = %d\n", which);
  //  eval( ((sexp_type*) P(1])->sexp , env);

  res = eval(PLANG(which)->sexp, env);
  size = P0INT(DISTR_NROW) * P0INT(DISTR_NCOL);
  for (i=0; i<size; i++) {
    // printf("i=%d %d\n", i, size);
    Res[i] = REAL(res)[i]; 
  }
}

void distrD(double *x, model *cov, double *v) {
  addVariable((char*) "x", x, 1, 1, PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_DX, v);
}
void distrDlog(double *x, model *cov, double *v) {
  distrD(x, cov, v);
  *v = LOG(*v);
}

void distrDinverse(double VARIABLE_IS_NOT_USED  *v, model  VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED  *x, double  VARIABLE_IS_NOT_USED *y) {
  // v is here INPUT variable
  NotProgrammedYet("distrDinverse");  
}

void distrP(double *x, model *cov, double *v) {
  addVariable((char*) "q", x, 1, 1,  PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_PX, v);
}
void distrP2sided(double *x, double *y, model *cov, double *v) {
  int i,
    size = P0INT(DISTR_NROW) * P0INT(DISTR_NCOL);
  double w;
  assert(OWNXDIM(0) == size);
 
  if (ANYOWNDIM == 1) {
   double z = x != NULL ? *x : -*y;
    addVariable((char*) "q", &z, 1, 1, PENV(DISTR_ENV)->sexp);
    evaluateDistr(cov, DISTR_PX, &w);
    addVariable((char*) "q", y, 1, 1, PENV(DISTR_ENV)->sexp);
    evaluateDistr(cov, DISTR_PX, v);
    *v -= w;
  } else {
    NotProgrammedYet("multivariate families of distribution functions");

    double 
      Sign = 1.0;
    TALLOC_XX1(z, size);
    
    for (i=0; i<size; i++) v[i] = 0.0;
    while (true) {
      // Siebformel !! ueber 2^size Moeglichkeiten
      for (i=0; i<size; i++) v[i] = Sign * z[i];
    }
    END_TALLOC_XX1;
  }
}

void distrQ(double *x, model *cov, double *v) {
  if (*x < 0 || *x > 1) {*v = RF_NA; return;}
  addVariable((char*) "p", x, 1, 1, PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_QX, v);
}
void distrR(double *x, model *cov, double *v) {
  if (x != NULL) ERR("Conditional distribution not allowed yet");
  addVariable((char*) "n", &ONE, 1, 1, PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_RX, v);
}
void distrR2sided(double *x, double *y, model *cov, double *v) {
  if (x != NULL || y != NULL) 
    ERR("conditional distribution not allowed yet");
  addVariable((char*) "n", &ONE, 1, 1, PENV(DISTR_ENV)->sexp);
  evaluateDistr(cov, DISTR_RX, v);
}
void kappa_distr(int i, model *cov, int *nr, int *nc){
  *nc = *nr = i < DefList[COVNR].kappas ? SIZE_NOT_DETERMINED : -1;
}
int check_distr(model *cov) {
  ASSERT_ONESYSTEM;
  ASSERT_CARTESIAN;

  kdefault(cov, DISTR_NROW, 1);
  kdefault(cov, DISTR_NCOL, 1);
  VDIM0 = P0INT(DISTR_NROW);
  VDIM1 = P0INT(DISTR_NCOL);
 
  EXTRA_STORAGE;

  RETURN_NOERROR;
}
int init_distr(model VARIABLE_IS_NOT_USED *cov, 
	       gen_storage VARIABLE_IS_NOT_USED *s){
  int err = NOERROR;
  cov->mpp.maxheights[0] = RF_NA;
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
  cov->mpp.unnormedmass = RF_NAN;
  RETURN_ERR(err);
}
void do_distr_do(model *cov, double *v){ 
  distrR(NULL, cov, v);
}
void range_distr(model *cov, range_type *range){
#define distr_n 5
  int idx[distr_n] = {DISTR_DX, DISTR_PX, DISTR_QX, DISTR_RX, DISTR_ENV};
  for (int i=0; i<distr_n; i++) {
    int j = idx[i];
    range->min[j] = RF_NAN;
    range->max[j] = RF_NAN;
    range->pmin[j] = RF_NAN;
    range->pmax[j] = RF_NAN;
    range->openmin[j] = false;
    range->openmax[j] = false; 
  }

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
    kappas = DefList[COVNR].kappas;
  for (i=DISTR_LAST + 1; i<kappas; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = 1e10;
    range->pmax[i] = - 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true; 
  }

}


#define GAUSS_PARAMETER_BASICS					\
  double *m = P(GAUSS_DISTR_MEAN),		\
    *sd = P(GAUSS_DISTR_SD)					\

#define GAUSS_PARAMETERS						\
  GAUSS_PARAMETER_BASICS;						\
  int  i,  mi,  si,\
     len_mean = cov->nrow[GAUSS_DISTR_MEAN],	\
     len_sd = cov->nrow[GAUSS_DISTR_SD],		\
    dim = ANYDIM	

#define FOR for(mi=si=i=0; i<dim; i++, mi=(mi + 1) % len_mean, si=(si + 1) % len_sd)

void gaussDlog(double *x, model *cov, double *v) {
  GAUSS_PARAMETERS;
  int returnlog = true;
  *v = 0.0;
  FOR {
    *v += dnorm(x[i], m[mi], sd[si], returnlog);
    // printf("i=%d x=%10g m=%10g sd=%10g  v=%10g\n", i, x[i], m[mi], sd[si],  *v);
  }

}

void gaussD(double *x, model *cov, double *v) {
  GAUSS_PARAMETERS;
  int returnlog = P0INT(GAUSS_DISTR_LOG);
  if (returnlog) { 
    gaussDlog(x, cov, v);
  } else {        
    *v = 1.0;
    FOR *v *= dnorm(x[i], m[mi], sd[si], returnlog);
  }
}

void gaussDinverse(double *v, model *cov, double *left, double *right) {
  // v is here INPUT variable !!
  GAUSS_PARAMETERS;
  FOR {
    double dummy = - 2.0 * LOG(*v * SQRTTWOPI * sd[si]);
    if (dummy < 0.0) {
      left[i] = right[i] = m[mi];
    } else {
      dummy = sd[mi] * SQRT(dummy);
      left[i] = m[mi] - dummy;
      right[i] = m[mi] + dummy;
    }
  }
}

void gaussP(double *x, model *cov, double *v) {
  GAUSS_PARAMETERS; 
  int returnlog = P0INT(GAUSS_DISTR_LOG);
  if (returnlog) {
    *v = 0.0;
    FOR *v += pnorm(x[i], m[mi], sd[si], true, returnlog);
  } else {
    *v = 1.0;
    FOR *v *= pnorm(x[i], m[mi], sd[si], true, returnlog);
  }
}
void gaussP2sided(double *x, double *y, model *cov, double *v) {
  GAUSS_PARAMETERS; 
  int returnlog = P0INT(GAUSS_DISTR_LOG);
  assert(y != NULL);
  
  if (x == NULL) { // symmetric 2-sided with a=-b
    if (returnlog) {
      *v = 0.0;
      FOR *v += y[i] != 0.0 
	? LOG(2.0 * pnorm(y[i], m[mi], sd[si], true, false) - 1.0)
	: dnorm(y[i], m[mi], sd[si], returnlog);
    } else {
      *v = 1.0;
      FOR {
	if (y[i] != 0.0)
	  *v *= 2.0 * pnorm(y[i], m[mi], sd[si], true, false) - 1.0;
	else *v *= dnorm(y[i], m[mi], sd[si], returnlog);
      }
    }
  } else {
    if (returnlog) {
      *v = 0.0;
      FOR {
	if (x[i] != y[i])
	  *v += LOG(pnorm(y[i], m[mi], sd[si], true, false) -
		    pnorm(x[i], m[mi], sd[si], true, false));
	else *v += dnorm(y[i], m[mi], sd[si], returnlog);
      }
    } else {
      *v = 1.0;
      FOR {
	if (x[i] != y[i]) *v *= (pnorm(y[i], m[mi], sd[si], true, false) -
				 pnorm(x[i], m[mi], sd[si], true, false));
	else *v *= dnorm(y[i], m[mi], sd[si], returnlog);
      }
    }
  }
//   //  printf("p2 %10g\n", *v);
}	   
void gaussQ(double *x, model *cov, double *v) {
  if (*x < 0 || *x > 1) {*v = RF_NA; return;}
  GAUSS_PARAMETER_BASICS; 
  int returnlog = P0INT(GAUSS_DISTR_LOG);
  *v = qnorm(x[0], m[0], sd[0], true, returnlog);
}		   

void gaussR(double *x, model *cov, double *v) {
  GAUSS_PARAMETERS; 
  if (x == NULL) {
    FOR v[i] = rnorm(m[mi], sd[si]);
  } else {
    FOR {
      if (R_FINITE(x[i])) v[i] = x[i];
      else v[i] = rnorm(m[mi], sd[si]);
    }
  }
}


void gaussR2sided(double *x, double *y, model *cov, double *v) {
  GAUSS_PARAMETERS; 
  if (x == NULL) {
    FOR {
      while (FABS( (v[i] = rnorm(m[mi], sd[si])) ) > y[i]); 
    }
  } else {
    FOR { 
      while ((v[i] = rnorm(m[mi], sd[si])) < x[i] || v[i] > y[i]);
    }
  }
}

void kappa_gauss_distr(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i == GAUSS_DISTR_SD || i== GAUSS_DISTR_MEAN) {
    *nc = 1;
    *nr = SIZE_NOT_DETERMINED;
  } else
    *nc = *nr = i == GAUSS_DISTR_LOG ? 1 : -1;
}

int check_gauss_distr(model *cov) {
  ASSERT_UNREDUCED;
  ASSERT_CARTESIAN;
  GAUSS_PARAMETER_BASICS;						
 
  if (m == NULL) kdefault(cov, GAUSS_DISTR_MEAN, 0.0);
  if (sd == NULL) kdefault(cov, GAUSS_DISTR_SD, 1.0);
  kdefault(cov, GAUSS_DISTR_LOG, false);

   VDIM0 = PREVXDIM(0);
  VDIM1 = 1;

  RETURN_NOERROR;
}

int init_gauss_distr(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  GAUSS_PARAMETERS; 
  if (cov->mpp.moments >= 0) {
    cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0; 
    if (cov->mpp.moments >= 1) {
      if (dim > 1) SERR("multivariate moment cannot be calculated");
      cov->mpp.mM[1] = cov->mpp.mMplus[1] = m[0];
      if (cov->mpp.moments >= 2) {
	double SD = sd == NULL ? 1.0 : sd[0];      
	cov->mpp.mM[2] = cov->mpp.mMplus[2] = SD * SD + m[0] * m[0] ;
      }
    }
  }

  // maxheight based on a Gaussian shape function with the same
  // sd's but normed to 1 at the origin
  cov->mpp.maxheights[0] = intpow(INVSQRTTWOPI, dim);
  FOR cov->mpp.maxheights[0] /= sd[si];
  //  cov->mpp.unnormedmass = RF_NA; // / cov->mpp.maxheights[0];
  cov->mpp.unnormedmass = 1.0 / cov->mpp.maxheights[0];

					   
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
  RETURN_NOERROR;
}

void do_gauss_distr(model *cov, double *v){
  double 
    *sd = P(GAUSS_DISTR_SD);				
  int  i,  mi,  si,
    len_mean = cov->nrow[GAUSS_DISTR_MEAN],	
    len_sd = cov->nrow[GAUSS_DISTR_SD],		
    dim = ANYOWNDIM;		
 
  cov->mpp.maxheights[0] = intpow(SQRTTWOPI, -dim);
  FOR cov->mpp.maxheights[0] /= sd[si];
  gaussR(NULL, cov, v);
}

void range_gauss_distr(model VARIABLE_IS_NOT_USED  *cov, range_type *range){
  range->min[GAUSS_DISTR_MEAN] = RF_NEGINF;
  range->max[GAUSS_DISTR_MEAN] = RF_INF;
  range->pmin[GAUSS_DISTR_MEAN] = - 1e8;
  range->pmax[GAUSS_DISTR_MEAN] = 1e8;
  range->openmin[GAUSS_DISTR_MEAN] = true;
  range->openmax[GAUSS_DISTR_MEAN] = true;
  
  range->min[GAUSS_DISTR_SD] = 0;
  range->max[GAUSS_DISTR_SD] = RF_INF;
  range->pmin[GAUSS_DISTR_SD] = 1e-10;
  range->pmax[GAUSS_DISTR_SD] = 1e8;
  range->openmin[GAUSS_DISTR_SD] = true;
  range->openmax[GAUSS_DISTR_SD] = true;
  
  booleanRange(GAUSS_DISTR_LOG);
}




#define SPHERIC_SPACEDIM 0
#define SPHERIC_BALLDIM 1
#define SPHERIC_RADIUS 2
#define SPHERIC_RADIUS 2
double random_spheric(int logicaldim, int balldim) {
  double r2 = - 1.0;    
  while (r2 < 0.0) {
    r2 = 1.0;
    for (int d=logicaldim; d<balldim; d++) {
      double dummy = UNIFORM_RANDOM;
      r2 -= dummy * dummy;
    }
  }
  return 0.5 * SQRT(r2);
}


void sphericD(double VARIABLE_IS_NOT_USED *x, model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED  *v) {
  ERR("density of 'RRspheric' cannot be calculated yet");
}
void sphericDinverse(double *v, model  *cov, double *left, double  *right) {
  if (v==NULL || *v <= 0.0) {
    *left = 0.0;
    *right = (0.5 * P0(SPHERIC_RADIUS));
  } else {
    //v is here INPUT variable !!
    ERR("density of 'RRspheric' cannot be calculated yet");
  }
}
void sphericDlog(double VARIABLE_IS_NOT_USED  *x, model VARIABLE_IS_NOT_USED  *cov, double VARIABLE_IS_NOT_USED  *v) {
  ERR("density of 'RRspheric' cannot be calculated yet");
}
void sphericP(double VARIABLE_IS_NOT_USED  *x, model VARIABLE_IS_NOT_USED   *cov, double VARIABLE_IS_NOT_USED  *v) {
  ERR("density of 'RRspheric' cannot be calculated yet");
}
void sphericQ(double VARIABLE_IS_NOT_USED *x, model VARIABLE_IS_NOT_USED  *cov, double VARIABLE_IS_NOT_USED  *v) {
  if (*x < 0 || *x > 1) {*v = RF_NA; return;}
  ERR("density of 'RRspheric' cannot be calculated yet");
}
void sphericR(double *x, model *cov, double *v) {
  int 
    dim = P0INT(SPHERIC_SPACEDIM),
    balldim = P0INT(SPHERIC_BALLDIM);  
  if (x==NULL) *v = random_spheric(dim, balldim) * P0(SPHERIC_RADIUS); // q[>0] ist E
  else ERR("conditional distribution cannot be calculated for sphericP.");
}

int check_RRspheric(model *cov) {
  ASSERT_UNREDUCED;
  ASSERT_CARTESIAN;
  int err;

  kdefault(cov, SPHERIC_SPACEDIM, 1);
  kdefault(cov, SPHERIC_BALLDIM, P0INT(SPHERIC_SPACEDIM));
  kdefault(cov, SPHERIC_RADIUS, 1.0);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);

  if (OWNLOGDIM(0) != 1) SERR("only dimension 1 allowed");


  VDIM0 = PREVXDIM(0);
  VDIM1 = 1;

  RETURN_NOERROR;
}

void range_RRspheric(model VARIABLE_IS_NOT_USED *cov, range_type *range){
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

  range->min[SPHERIC_RADIUS] = 0; 
  range->max[SPHERIC_RADIUS] = RF_INF;
  range->pmin[SPHERIC_RADIUS] = 0.001;
  range->pmax[SPHERIC_RADIUS] = 1000;
  range->openmin[SPHERIC_RADIUS] = true;
  range->openmax[SPHERIC_RADIUS] = true;
}


int init_RRspheric(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int i, m,
    nm = cov->mpp.moments,
    nmvdim = nm + 1, // vdim == 1 
    dim = P0INT(SPHERIC_SPACEDIM),
    balldim = P0INT(SPHERIC_BALLDIM),
     testn = GLOBAL.mpp.n_estim_E;

  double scale, dummy,
    *M = cov->mpp.mM,  
    radius = P0(SPHERIC_RADIUS),
    *Mplus = cov->mpp.mMplus;
 
  // printf("ball mppM2 %10g\n", cov->mpp.mM2);
  //printf("%d %d %d\n", testn, GLOBAL.mpp.n_estim_E,
//   //	 P0INT(SPHERIC_BALLDIM));// assert(false);
  
  for (M[0] = 1.0, m=1; m < nmvdim; M[m++] = 0.0);
  for (i=0; i<testn; i++) { // take random sample of scales
    // printf("testn %d %10g %d\n", i, scale, testn);
    scale = random_spheric(dim, balldim);
    dummy = 1.0;
    for (m=1; m < nmvdim; m++) {
      dummy *= scale;
      M[m] += dummy;
    }
  }
  for (dummy=radius, m=1; m < nmvdim; m++, dummy *= radius)  {
    M[m] = dummy * (double) testn;
    Mplus[m] = M[m];
    //
    //printf("M: %d %10g\n", d, M[d]);
  }

//   // exakte Formel; check ob das gleiche wie simu
  if (PL > 1) {
    dummy = 0.5 * OWNLOGDIM(0) + 1.0;
    PRINTF("init_spheric %10g %10g %10g\n", M[nm], 
	   EXP( (balldim - dim) * M_LN_SQRT_PI // LOG(SQRT(pi)) 
		+ lgammafn(dummy) - lgammafn(0.5 * balldim + 1)),
	   EXP( - dim * M_LN_SQRT_PI // LOG(SQRT(pi)) 
		+ lgammafn(dummy))
	   );
  }
  //cov->mpp.refradius = RF_NA;
  cov->mpp.maxheights[0] = RF_NA;
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;

  RETURN_NOERROR;
}

void do_RRspheric(model *cov, double *v) {
  sphericR(NULL, cov, v);
}



// ***** deterministic ****
#define DETERM_MEAN 0
#define DETERM_PARAMETERS						\
  double *mean = P(DETERM_MEAN);					\
  int i,  mi,			\
     dim = OWNTOTALXDIM,				\
    len_mean = cov->nrow[DETERM_MEAN]

#define DETERMFOR for(mi=i=0; i<dim; i++, mi=(mi + 1) % len_mean)


void determD(double *x, model *cov, double *v) { 
  DETERM_PARAMETERS;
  DETERMFOR if (x[i] != mean[mi]) {
    *v = 0.0; 
    return;
  }
  *v = RF_INF;
}

void determDlog(double *x, model *cov, double *v) { 
  determD(x, cov, v);
  *v = LOG(*v);
}

void determDinverse(double VARIABLE_IS_NOT_USED *v, model *cov, double *left, double *right) {
  DETERM_PARAMETERS;
  DETERMFOR left[i] = right[i] = mean[mi];
}

void determP(double *x, model *cov, double *v) { 
  DETERM_PARAMETERS;
  DETERMFOR {
    if (x[i] < mean[mi]) {
      *v = 0.0; 
      return;
    }
  }
  *v = 1.0; 
}

void determP2sided(double *x, double *y, model *cov, double *v) {   
  DETERM_PARAMETERS;
  assert(y != NULL);  
  *v = 1.0; 
  if (x == NULL) {
    DETERMFOR {
      if (y[i] == 0.0 && mean[mi] == 0.0) *v = RF_INF;
      else if (-y[i] > mean[mi] || y[i] < mean[mi]) {
	*v = 0.0; 
	return;
      }
    }
  } else {
    DETERMFOR {
      if (x[i] == y[i] && mean[mi] == x[i]) *v = RF_INF;
      else if (x[i] > mean[mi] || y[i] < mean[mi]) {
	*v = 0.0; 
	return;
      }
    }
  }
}

void determQ(double *x, model *cov, double *v) { 
  if (*x < 0 || *x > 1) {*v = RF_NA; return;}
  v[0] = P0(DETERM_MEAN);
}

void determR(double *x, model *cov, double *v) { 
  DETERM_PARAMETERS;
  if (x==NULL) DETERMFOR v[i] = mean[i]; 
  else 
    DETERMFOR v[i] = !R_FINITE(x[i]) || x[i] == mean[mi] ? mean[mi] : RF_NA;
}

void kappa_determ(int i, model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = i == 0 ? OWNTOTALXDIM  : i == 1 ? 1 : -1;
}

void determR2sided(double *x, double *y, model *cov, double *v) { 
  DETERM_PARAMETERS;
  if (x == NULL) {
    DETERMFOR v[i] = FABS(y[i]) > mean[mi] ? mean[mi] : RF_NA;    
  } else {
    DETERMFOR v[i] = x[i] < mean[mi] && y[i] > mean[mi] ? mean[mi] : RF_NA;
  }
}

int check_determ(model *cov) {
  ASSERT_UNREDUCED;
  ASSERT_CARTESIAN;
  int 
    dim = OWNTOTALXDIM;  
  
  if (PisNULL(DETERM_MEAN)) kdefault(cov, DETERM_MEAN, 0.0);
  
  VDIM0 = dim;
  VDIM1 = 1;

  RETURN_NOERROR;
}


void range_determ(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[DETERM_MEAN] = RF_NEGINF;
  range->max[DETERM_MEAN] = RF_INF;
  range->pmin[DETERM_MEAN] = - 1e10;
  range->pmax[DETERM_MEAN] = 1e10;
  range->openmin[DETERM_MEAN] = true;
  range->openmax[DETERM_MEAN] = true;
}

int init_determ(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED  *s) {
  int err = NOERROR;

  // normed against an extremal gaussian process with shape that 
  // is identically constant 1 in space
  cov->mpp.maxheights[0] = 1.0;
  cov->mpp.unnormedmass = 1.0;
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
  RETURN_ERR(err);
}

void do_determ(model *cov, double *v) {
  determR(NULL, cov, v);
}




// ***** Set ****


void setParamD(double *x, model *cov, double *v) { 
VTLG_D(x, cov->sub[SETPARAM_LOCAL], v);
}

void setParamDlog(double *x, model *cov, double *v) { 
  VTLG_DLOG(x, cov->sub[SETPARAM_LOCAL], v);
}

void setParamDinverse(double *v, model *cov, double *left, double *right) {
  NONSTATINVERSE_D(v, cov->sub[SETPARAM_LOCAL], left, right);
}

void setParamP(double *x, model *cov, double *v) {   
  VTLG_P(x, cov->sub[SETPARAM_LOCAL], v);
}

void setParamP2sided(double *x, double *y, model *cov, double *v) {   
  VTLG_P2SIDED(x, y, cov->sub[SETPARAM_LOCAL], v);
}

void setParamQ(double *x, model *cov, double *v) {   
  VTLG_Q(x, cov->sub[SETPARAM_LOCAL], v);
}

void setParamR(double *x, model *cov, double *v) {   
  VTLG_R(x, cov->sub[SETPARAM_LOCAL], v);
}

void setParamR2sided(double *x, double *y, model *cov, double *v) {   
  VTLG_R2SIDED(x, y, cov->sub[SETPARAM_LOCAL], v);
}

int check_setParam(model *cov) {
  ASSERT_UNREDUCED;
  ASSERT_CARTESIAN;
  model *next= cov->sub[SETPARAM_LOCAL];
  int err,
    dim = OWNTOTALXDIM;

  kdefault(cov, SET_PERFORMDO, true);

  if ((err = CHECK_R(next, dim)) != NOERROR) RETURN_ERR(err);

  setbackward(cov, next);
  VDIM0 = next->vdim[0];
  VDIM1 = next->vdim[1];

  TaylorCopy(cov, next);
  cov->mpp.maxheights[0] = next->mpp.maxheights[0];
  cov->mpp.unnormedmass = next->mpp.unnormedmass;
  if (cov->mpp.moments >= 1) {
    cov->mpp.mM[0] = cov->mpp.mMplus[0];     
    cov->mpp.mMplus[0] = next->mpp.mMplus[0];
  }
  
  RETURN_NOERROR;
}


int init_setParam(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  model *next= cov->sub[SETPARAM_LOCAL];
  set_storage *X = cov->Sset;
  assert(X != NULL);
  int err;
  if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR)
    RETURN_ERR(err);
  if (X->remote != NULL) {
    assert(X->set != NULL);
    X->set(cov->sub[0], X->remote, X->variant);
  }
  TaylorCopy(cov, next);
  cov->mpp.maxheights[0] = next->mpp.maxheights[0];
  cov->mpp.unnormedmass = next->mpp.unnormedmass;
   if (cov->mpp.moments >= 1) {
    cov->mpp.mM[0] = next->mpp.mM[0];   
    cov->mpp.mMplus[0] = next->mpp.mMplus[0];
  }
 RETURN_NOERROR;
}

void do_setParam(model *cov, double *v) {
  model *next = cov->sub[SETPARAM_LOCAL];
  bool performDo = P0INT(SET_PERFORMDO);
  if (performDo) { DORANDOM(next, v); } 
  //cov->mpp.maxheights[0] = next->mpp.maxheights[0];
}


void range_setParam(model *cov, range_type *range){
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
  model *next = cov->sub[0];		\
  int dim=OWNTOTALXDIM;		\
  double *loc = P(LOC_LOC),			\
    *scale= P(LOC_SCALE)			        
  

#define LOC_PARAMETERS				\
  LOC_PARAMETER_BASICS;				\
  int i, mi, si,				\
    len_loc = cov->nrow[LOC_LOC],		\
    len_scale = cov->nrow[LOC_SCALE];		\
  assert(OWNTOTALXDIM == OWNLOGDIM(0))

#define LOCFOR for(mi=si=i=0; i<dim; i++, mi=(mi + 1) % len_loc, si=(si + 1) % len_scale)
 
void locD(double *x, model *cov, double *v) {
  LOC_PARAMETERS;
  double prod = 1.0;
  TALLOC_X1(z, dim);
    
  LOCFOR { 
    z[i] = (x[i] - loc[mi]) / scale[si];
    prod *= scale[si];
  }
  VTLG_D(z, next, v);
  *v /= prod;
  END_TALLOC_X1;
}

void locDlog(double *x, model *cov, double *v) {
  locD(x, cov, v);
  *v = LOG(*v);
}

void locDinverse(double *v, model *cov, double *left, double *right) {
  LOC_PARAMETERS;
  NONSTATINVERSE_D(v, next, left, right);
  LOCFOR {
    left[i] = left[i] * scale[si] + loc[mi];
    right[i] = right[i] * scale[si] + loc[mi];
  }
}

void locP(double *x, model *cov, double *v) {
  LOC_PARAMETERS;
  TALLOC_X1(z, dim);
  LOCFOR z[i] = (x[i] - loc[mi]) / scale[si];
  VTLG_P(z, next, v);
  END_TALLOC_X1;
}

void locP2sided(double *x, double *y, model *cov, double *v) {
  LOC_PARAMETERS;
  TALLOC_X1(z, dim);
  if (x == NULL) {
    LOCFOR {
      z[i] = (y[i] - loc[mi]) / scale[si];
    }
    VTLG_P2SIDED(NULL, z, next, v);

    if (*v > 0 && *v < RF_INF) {
      LOCFOR {
	if (z[i] == 0) *v = *v / scale[si];
      }
    }
    
    //  if (!true || *v > 1.0 / 100.0 && *v != 1.0) {
    //  LOCFOR {
    //	printf("loc %d %10g %10g s=%10g %10g\n", i, y[i], z[i],  scale[si], loc[mi]);
    // }
    //  printf("v=%10e\n", *v);   
    //}
  } else {    
    TALLOC_X2(zy, dim);
    LOCFOR { 
      z[i] = (x[i] - loc[mi]) / scale[si];
      zy[i] = (y[i] - loc[mi]) / scale[si];
    }
    VTLG_P2SIDED(z, zy, next, v);

    if (*v > 0 && *v < RF_INF) {
      LOCFOR {
	if (z[i] == zy[i]) *v = *v / scale[si];
      }
    }
    END_TALLOC_X2;
  }
  END_TALLOC_X1;
}

void locQ(double *x, model *cov, double *v) {
  LOC_PARAMETER_BASICS;  
  if (dim != 1) BUG;
  VTLG_Q(x, next, v);
  v[0] = v[0] * scale[0] + loc[0]; 
}		   

void locR(double *x, model *cov, double *v) {
  LOC_PARAMETERS;
  TALLOC_DOUBLE(z1);
  if (x != NULL) {
    TALLOC_GLOBAL_X1(z1, dim);
    LOCFOR z1[i] = (x[i] - loc[mi]) / scale[si];
  }
  VTLG_R(z1, next, v);
  if (x == NULL) LOCFOR v[i] = v[i] * scale[si] + loc[mi]; 
  else LOCFOR v[i] = R_FINITE(x[i]) ? x[i] : v[i] * scale[si] + loc[mi]; 
  FREE_TALLOC(z1);
}

void locR2sided(double *x, double *y, model *cov, double *v) {
  LOC_PARAMETERS;  
  TALLOC_DOUBLE(z1);
  if (x != NULL) {
    TALLOC_GLOBAL_X1(z1, dim);
    LOCFOR z1[i] = (x[i] - loc[mi]) / scale[si];
  } // !! Z1 can be != NULL as used differently among the loc-fcts
  TALLOC_X2(z2, dim);
  assert(y != NULL);
  LOCFOR z2[i] = (y[i] - loc[mi]) / scale[si];
  
  VTLG_R2SIDED(z1, z2, next, v);
  //     printf("locR dim=%d scale=%10g %10g %10g %10g -> %10g %10g %10g %d -> %10g %10g %10g\n", dim, *scale, y[0], y[1], y[2], z2[0], z2[1], z2[2], z1==NULL, v[0], v[1], v[2]);
 
  LOCFOR {
    v[i] = v[i] * scale[si] + loc[mi]; 
  }
  FREE_TALLOC(z1);
  END_TALLOC_X2;
}

void kappa_loc(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i == LOC_SCALE || i== LOC_LOC) {
    *nc = 1;
    *nr = SIZE_NOT_DETERMINED;
  } else if (i == LOC_POWER) *nc = * nr = 1;
  else *nc = *nr = OUT_OF_RANGE;
}

int check_loc(model *cov) {
  ASSERT_UNREDUCED;
  ASSERT_CARTESIAN; 
  model *next = cov->sub[0];
  int
    dim = OWNTOTALXDIM;
  double 
    *loc = P(LOC_LOC),
    *scale = P(LOC_SCALE);
  int err;

  kdefault(cov, POWPOWER, 0.0);

  //PMI(cov);

  if ((err = CHECK_R(next, dim)) != NOERROR) RETURN_ERR(err);

  // assert(false); /// not checked yet
 
  if (loc == NULL) kdefault(cov, LOC_LOC, 0.0);
  if (scale == NULL) kdefault(cov, LOC_SCALE, 1.0);
 
  VDIM0 = next->vdim[0];
  VDIM1 = next->vdim[1];

  EXTRA_STORAGE;
  RETURN_NOERROR;
}

int init_loc(model *cov, gen_storage *s){
  // defn *C = DefList + COVNR;
  LOC_PARAMETERS;
  int err;
  double p = P0(LOC_POWER); // same frame as POWPOWER of '$power'
  
  if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR) RETURN_ERR(err);
 
  if (cov->mpp.moments >= 0) {
    cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0; 
    if (cov->mpp.moments >= 1) {
      if (dim > 1) {
	LOCFOR if (scale[si] != 1.0 || loc[mi]!=0.0)
	  SERR("multivariate moment cannot be calculated");
      }
      cov->mpp.mM[1] = cov->mpp.mM[1] * scale[0] + loc[0];
      cov->mpp.mMplus[1] = loc[0] == 0.0 ? cov->mpp.mMplus[1] * scale[0] : RF_NA;
      if (cov->mpp.moments >= 2) {
	double ssq = scale[0] * scale[0];
	cov->mpp.mM[2] =
	  cov->mpp.mM[2] * ssq + loc[0] * (2.0 * cov->mpp.mM[1] - loc[0]);
	cov->mpp.mMplus[1] = loc[0] == 0.0 ? cov->mpp.mMplus[1] * ssq : RF_NA;
      }
    }
  }

  // normed against a shape function with the same scale modification
  // if (R_FINITE(next->mpp.unnormedmass))
  cov->mpp.unnormedmass = next->mpp.unnormedmass * POW(scale[0], dim + p); 
  cov->mpp.maxheights[0] = next->mpp.maxheights[0] / POW(scale[0], dim);

  cov->mpp.mM[0] = next->mpp.mM[0];
  cov->mpp.mMplus[0] = next->mpp.mMplus[0];

  RETURN_NOERROR;
}

void do_loc(model *cov, double *v){
  //  model *next = cov->sub[0];
  //double *scale= P(LOC_SCALE);
  DORANDOM(cov->sub[0], v);
  locR(NULL, cov, v);
  // cov->mpp.maxheights[0] = next->mpp.maxheights[0] * scale[0];
}

void range_loc(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[LOC_LOC] = RF_NEGINF;
  range->max[LOC_LOC] = RF_INF;
  range->pmin[LOC_LOC] = - 1e8;
  range->pmax[LOC_LOC] = 1e8;
  range->openmin[LOC_LOC] = true;
  range->openmax[LOC_LOC] = true;
  
  range->min[LOC_SCALE] = 0;
  range->max[LOC_SCALE] = RF_INF;
  range->pmin[LOC_SCALE] = 1e-10;
  range->pmax[LOC_SCALE] = 1e8;
  range->openmin[LOC_SCALE] = true;
  range->openmax[LOC_SCALE] = true;

  range->min[LOC_POWER] = RF_NEGINF;
  range->max[LOC_POWER] = RF_INF;
  range->pmin[LOC_POWER] = -OWNTOTALXDIM;
  range->pmax[LOC_POWER] = +OWNTOTALXDIM;
  range->openmin[LOC_POWER] = true;
  range->openmax[LOC_POWER] = true;
}



#define UNIF_PARAMETER_BASICS				\
  double						\
   *min = (double*) P(UNIF_MIN),	\
     *max= (double*) P(UNIF_MAX)		
	

#define UNIF_PARAMETERS							\
  UNIF_PARAMETER_BASICS;						\
  int  i,						\
     mini,						\
     maxi,						\
     len_min = cov->nrow[UNIF_MIN],			\
     len_max = cov->nrow[UNIF_MAX],		\
    dim = ANYOWNDIM

#define UNIFOR  for(mini=maxi=i=0; i<dim; i++, mini=(mini + 1) % len_min, \
		      maxi=(maxi + 1) % len_max)
 
// Abgeschlossene Intervalle!! wichtig !!
void unifD(double *x, model *cov, double *v) {
  UNIF_PARAMETERS;
  double area = 1.0;
  bool normed = P0INT(UNIF_NORMED);
  UNIFOR {
    if (x[i] < min[mini] || x[i] > max[maxi]) { *v = 0.0; return; }
    else if (normed) area *= max[maxi] - min[mini];  
  }
  *v = 1.0 / area;
}


void unifDlog(double *x, model *cov, double *v) {
  unifD(x, cov, v);
  *v = LOG(*v);
}

void unifDinverse(double *v, model *cov, double *left, double *right) {
  UNIF_PARAMETERS;
   double area = 1.0;
   if (P0INT(UNIF_NORMED)) UNIFOR area *= max[maxi] - min[mini];  
   if (*v * area > 1){
     UNIFOR left[i] = right[i] = 0.5 * (max[maxi] + min[mini]);
   } else {
     UNIFOR {
       left[i] = min[mini];
       right[i] = max[maxi];
       //  printf("unif=%d %10g %10g\n", i, min[mini], max[maxi]);
     }
   }
}

void unifP(double *x, model *cov, double *v) {
  // distribution functions
  UNIF_PARAMETERS;
 
  double area = 1.0;
  bool normed = P0INT(UNIF_NORMED);
  UNIFOR {
    if (x[i] <= min[mini]) { *v = 0.0; return; }
    if (x[i] < max[maxi]) area *= x[i] - min[mini];
    if (normed) area /= max[maxi] - min[mini];  
    // ekse area *= 1.0;
  }
  *v = area;
}


void unifP2sided(double *x, double *y, model *cov, double *v) {
  // distribution functions
  UNIF_PARAMETERS;
  double a, b;
  bool normed = P0INT(UNIF_NORMED);

  double area = 1.0;
  UNIFOR {
    a = x != NULL ? x[i] : -y[i];
    if (a != y[i]) {
      if (a < min[mini]) a =  min[mini];
      b = y[i] > max[maxi] ? max[maxi] : y[i];
      if (a >= b) { *v = 0.0; return; }
      area *= b - a;  
    } else {
      if (a < min[mini] || a > max[maxi]) { *v = 0.0; return; }
    }
    if (normed) area /= max[maxi] - min[mini];  
  }
  *v = area;
}


void unifQ(double *x, model *cov, double *v) {
  if (*x < 0 || *x > 1) {*v = RF_NA; return;}
  UNIF_PARAMETER_BASICS;
  if (P0INT(UNIF_NORMED)) {
    *v = min[0] + *x * (max[0] - min[0]);
  } else {
    *v = min[0] + *x;
  }
}


void unifR(double *x, model *cov, double *v) { 
  UNIF_PARAMETERS;
  if (x == NULL) {
    UNIFOR {
      v[i] = min[mini] + UNIFORM_RANDOM * (max[maxi] - min[mini]);
     }
  } else {
    UNIFOR {
      if (!R_FINITE(x[i]))
	v[i] = min[mini] + UNIFORM_RANDOM * (max[maxi] - min[mini]);
      else v[i] = x[i] >= min[mini] && x[i] <= max[maxi] ? x[i] : RF_NA;
    }
  }
}


void unifR2sided(double *x, double *y, model *cov, double *v) { 
  UNIF_PARAMETERS;
  double a, b;
 
  UNIFOR {
    if (x != NULL) a =  x[i] < min[mini] ? min[mini] : x[i];
    else a = -y[i] < min[mini] ? min[mini] : -y[i];
    b = y[i] > max[maxi] ? max[maxi] : y[i];
    if (a > b) ERR("parameters of 2-sided unifR out of range");      
    v[i] = a + UNIFORM_RANDOM * (b-a);
  }  
}


void kappa_unif(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i == UNIF_MIN || i == UNIF_MAX) {
    *nc = 1;
    *nr = SIZE_NOT_DETERMINED;
  } else if (i == UNIF_NORMED) {
    *nc = *nr = 1;
  } else *nc = *nr = OUT_OF_RANGE;
}


int check_unif(model *cov) {
  ASSERT_UNREDUCED;
  ASSERT_CARTESIAN;
   
  if (PisNULL(UNIF_MIN)) kdefault(cov, UNIF_MIN, 0.0);
  if (PisNULL(UNIF_MAX)) kdefault(cov, UNIF_MAX, 1.0);
  kdefault(cov, UNIF_NORMED, true); // see call of unif by special
  // functions where value false is needed.
  // maybe to write another function would be clearer.

  VDIM0 = OWNLOGDIM(0);
  VDIM1 = 1;

  RETURN_NOERROR;
}



int init_unif(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  UNIF_PARAMETERS;

  //  printf("entering init_unif\n");
 
   
  // normed against a shape function that is the indicator function 
  // of the respective window
  cov->mpp.unnormedmass = 1.0;
  UNIFOR { cov->mpp.unnormedmass *= max[maxi] - min[mini]; }
  if (P0INT(UNIF_NORMED)) {
    cov->mpp.maxheights[0] =  1.0 / cov->mpp.unnormedmass;
    if (cov->mpp.moments >= 0) {
      cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0; 
      if (cov->mpp.moments >= 1) {
	if (dim > 1) SERR("multivariate moment cannot be calculated");
	cov->mpp.mM[1] = 0.5 * (min[0] + max[0]);
	cov->mpp.mMplus[1] = max[0] <= 0.0 ? 0.0 : 0.5 * max[0];
	if (cov->mpp.moments >= 2) {
	  cov->mpp.mM[2] = max[0] - min[0];
	cov->mpp.mM[2] *= cov->mpp.mM[2] / 12.0;
	}
      }
    }    
  } else { // not normed
    cov->mpp.maxheights[0] = 1.0;
    cov->mpp.mM[0] = cov->mpp.mMplus[0] = cov->mpp.unnormedmass;
    if (cov->mpp.moments > 0) 
      SERR("unnormed unif does not allow for higher moments");
  }
  // PMI(cov);
  //    printf("leaving init_unif\n");
   
  RETURN_NOERROR;
}


void do_unif(model *cov, double *v){
  unifR(NULL, cov, v);
  // printf("v=%10g\n", *v);
  // cov->mpp.maxheights[0] = 1.0;
  // UNIFOR cov->mpp.maxheights[0] /= max[maxi] - min[mini];
}


void range_unif(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[UNIF_MIN] = RF_NEGINF;
  range->max[UNIF_MIN] = RF_INF;
  range->pmin[UNIF_MIN] = - 1e8;
  range->pmax[UNIF_MIN] = 1e8;
  range->openmin[UNIF_MIN] = true;
  range->openmax[UNIF_MIN] = true;
  
  range->min[UNIF_MAX] = RF_NEGINF;
  range->max[UNIF_MAX] = RF_INF;
  range->pmin[UNIF_MAX] = - 1e8;
  range->pmax[UNIF_MAX] = 1e8;
  range->openmin[UNIF_MAX] = true;
  range->openmax[UNIF_MAX] = true;

  booleanRange(UNIF_NORMED);

}


#define SURFACE(dim) ((dim) * intpow(2.0, dim)) //of a d-dim cube with edge length 2, i.e. a (dim-1) dimensional surface !


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
  //  printf("res = %10g\n", res);
  for (d=1; d<=squeezed_parts; d++) res *= xsort[d]; 

  // printf("res = %10g %10g %d %10g\n", res, xsort[1], red_dim,
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

  //printf("res=%10g %10g %10g %10g\n", res,  POW(end, pPd), POW(start, pPd), pPd);

  //return  * intpow(2 * start, squeezed_parts) 
  return res * (POW(end, pPd) - POW(start, pPd)) / pPd;
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

  //  printf("s=%10g start=%10g p=%10g sq=%d a=%10g %10g\n", s, start, p, squeezed_parts, a, s * POW(start, p));
  
 
  //  printf("expvol %10g %10g %10g a=%10g ret=%10g\n", 
  //       SURFACE(red_dim) / SURFACE(dim),  POW(s, squeezed_parts / p),
  //   incomplete_gamma(s * POW(start, p), s * POW(end, p), a), a,
  //   SURFACE(red_dim) / SURFACE(dim) * POW(s, squeezed_parts / p) * 
  // incomplete_gamma(s * POW(start, p), s * POW(end, p), a));

  return 
    SURFACE(red_dim) / SURFACE(dim) * POW(s, squeezed_parts / p) * 
    incomplete_gamma(s * POW(start, p), s * POW(end, p), a);
}

void RandomPointOnCubeRing(double start, double end, int dim, double *x) {
  double twostart = 2.0 * start, X, Y;
  int i0;

  //  printf("start = %10g %10g\n", start, end);
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
    //printf("x= %10g %10g %10g  start=%10g %10g Z=%10g\n", x[0], x[1],x[2], start, end, Z);
  }
    break;
  default : BUG;
  }
  
  //printf("x0= %10g %10g %10g dim=%d\n", x[0], x[1],x[2],dim);
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
#define IDX_OUTER (NSTEP + IDX_STEPS) // change also KeyInfo if changed
#define ASSIGN_IDX_INNER (-IDX_STEPS)
#define ASSIGN_IDX_OUTER (-IDX_STEPS - 1)
#define ASSIGN_IDX_MISMATCH (-IDX_STEPS - 2)
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
#define ASSIGNMENT(IDX) rect->asSign[IDX] // for control only


#define RECTANGULAR_PARAMETER_BASICS		\
  rect_storage *rect = cov->Srect;		\
  int VARIABLE_IS_NOT_USED dim = ANYOWNDIM;				\
  if (rect == NULL) BUG

#define RECTANGULAR_PARAMETERS				\
  RECTANGULAR_PARAMETER_BASICS;				\
  model VARIABLE_IS_NOT_USED *next = cov->sub[0];	\
 
#define RECTFOR  for(i=0; i<dim; i++)

#define DEBUG_CUMSUM false
void CumSum(double *y, bool kurz, model *cov, double *cumsum) {
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
  
  //  printf("y=%10g \n", y[0]);
  ysort[0] = 0.0;
  for (d=0; d<dim; d++) {
    //
    //    printf("%d y=%10g\n", d, y[d]);
    if (y[d] < 0.0) {  BUG; }
    ysort[d+1] = FABS(y[d]);
    //printf("y=%10g\n", y[d]);
  }

  int *idx = rect->idx;
  if (dim > 1) {
    Ext_ordering(ysort, dimP1, One, idx);
    assert(dim < 3 ||  (ysort[idx[0]] <= ysort[idx[1]] &&
			ysort[idx[1]] <= ysort[idx[2]]));
  } else {
    idx[0] = 0;
    idx[1] = 1;      
  }
  assert(idx[0] == 0 && idx[1] > 0);
  //for (d=0; d<=dim; d++) printf("%d,%10g ", idx[d], z[d]); printf("\n");
  for(d=0; d<dimP1; d++) zneu[d] = ysort[idx[d]];
  for(d=0; d<dimP1; d++) ysort[d] = zneu[d];

  TMP = 0;
  d = 1;
  dM1 = d - 1;
  if (zneu[d] >= INNER && use_weights) {
    kstep = zneu[d] >= OUTER ? NSTEP : (int) ((zneu[d] - INNER) / STEP);  

    //
    if (DEBUG_CUMSUM) {
      PRINTF("kstep %d %10g\n", kstep, (zneu[d] - INNER)/STEP);
    }

    if (kurz) {
      if (zneu[d] == RF_INF) {
	cumsum[TMP] = OUTER_CUM;
	TMP++;
	assert(cumsum[TMP-1] >= 0);
	return;
      }	
      RIGHT_END(TMP) = INNER + kstep * STEP;
      SQUEEZED(TMP)  = - MAXINT / 2;
      ASSIGNMENT(TMP)= ASSIGN_IDX_MISMATCH;
      //  
      if (DEBUG_CUMSUM) {
	PRINTF("0. TMP=%d kstep=%d ass=%d\n", TMP, kstep,  ASSIGNMENT(TMP));
      }
      
      cumsum[TMP++] = WEIGHT(kstep - 1);
      assert(cumsum[TMP-1] >= 0);
    } else {
      int k;
      //printf("nicht kurz %d %d\n", IDX_STEPS, kstep);
      for (k=-IDX_STEPS; k<kstep; k++) {
	RIGHT_END(TMP) = INNER + (k + IDX_STEPS) * STEP;
	SQUEEZED(TMP)  = dM1;
	ASSIGNMENT(TMP) = k == -IDX_STEPS ? ASSIGN_IDX_INNER : k;
	//
	if (DEBUG_CUMSUM) {
	  PRINTF("1. TMP=%d k=%d ass=%d\n", TMP, k,  ASSIGNMENT(TMP));
	}
	cumsum[TMP++]  = WEIGHT(k);
	assert(cumsum[TMP-1] >= 0);
      }
    }
    zneu[dM1] = INNER + STEP * kstep; 
    start_cumul = TMP;
  } else {
    kstep = 0;
    //    printf("AA d=%d %d %10g %10g\n", d, TMP, zneu[dM1], INNER);
  }
  
  for ( ; d<=dim; d++) {
    dM1 = d - 1;
    //  printf("X d=%d %d dim=%d\n", d, TMP, dim);
     if (zneu[d] == zneu[dM1]) continue;
    //printf("%10g == %10g %d\n",zneu[d],  zneu[dM1], zneu[d] == zneu[dM1] );
   //   printf("A d=%d %d dim=%d\n", d, TMP, dim);
 
    if (zneu[dM1] < INNER) {
      min = RIGHT_END(TMP) = zneu[d] <= INNER ? zneu[d] : INNER;
      SQUEEZED(TMP) = dM1;
      ASSIGNMENT(TMP) = ASSIGN_IDX_INNER;
      cumsum[TMP++] = INNER_CONST * 
	PoweredVolOfCube(ysort, zneu[dM1], min, INNER_POW, dim, dM1);
       // 
      if (DEBUG_CUMSUM) {
	PRINTF("2. TMP=%d kstep=%d asS=%d cum=%10g inner_c=%10g start=%10g ende=%10g pow=%10g dim=%d squ=%d\n", TMP-1, kstep,  ASSIGNMENT(TMP-1), cumsum[TMP-1], //
	       INNER_CONST, zneu[dM1],
	       min,PoweredVolOfCube(ysort, zneu[dM1], min, INNER_POW, dim, dM1),
	       dim, dM1);
      }
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
      
      //printf("d=%d min=%10g %10g %10g\n", d, min, zneu[dM1], zneu[d]);
      
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
	if (DEBUG_CUMSUM) {
	  PRINTF("3. TMP=%d kstep=%d asS=%d a=%4.2f st=%4.2f end=%4.2f val=%4.2e, tot=%4.4f %d\n", TMP-1, kstep,  ASSIGNMENT(TMP-1), a, STEP, RIGHT_END(TMP-1), //
		 VALUE(kstep), cumsum[TMP-1], ASSIGNMENT(0));
	}
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
	if (DEBUG_CUMSUM) {
	  PRINTF("4. TMP=%d k=%d ass=%d\n", TMP, kstep, ASSIGNMENT(TMP));
	}
	cumsum[TMP++] = VolumeOfCubeRing(ysort, a, min, dim, dM1) * VALUE(kstep);
	assert(cumsum[TMP-1] >= 0);
	if (zneu[d] <= OUTER) continue;
      }
      
      // printf("outer=%10g %10g %10g %10e %d %d %d\n", OUTER, a,  cumsum[TMP-1], 
      //        OUTER - a, OUTER == min, dM1, kstep < NSTEP);
      //assert(false);
      // printf("C %d %d\n", d, TMP);
    }
       

    RIGHT_END(TMP) = zneu[d];
    SQUEEZED(TMP) = dM1;
    ASSIGNMENT(TMP) = ASSIGN_IDX_OUTER;
    //  
    if (DEBUG_CUMSUM) {
      PRINTF("5. TMP=%d k=%d ass=%d %10g p=%10g right=%10g\n",
	     TMP, 9999, ASSIGNMENT(TMP), OUTER, OUTER_POW, zneu[d]);
    }
    
    if (next->finiterange == wahr) {
      cumsum[TMP++] = 0.0;
      continue;
    } else if (OUTER_POW > 0) {
      //      printf("here %d C=%10g Z=%10g PC=%10g e=%10g %10g\n", 
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
      if (DEBUG_CUMSUM) {
	PRINTF("6. TMP=%d k=%d ass=%d oc=%10g c=%10g op=%10g opc=%10g dim=%d %d right=%10g\n", TMP, 10001, ASSIGNMENT(TMP-1),OUTER_CONST, OUTER, OUTER_POW,
	       OUTER_POW_CONST, dim, dM1, zneu[d]);
      }
      //  printf("TMP=%d %10g\n", TMP, cumsum[TMP-1]);
      // if (cumsum[TMP-1] < 0) {
      //t i;
      //fori=1; i<TMP; i++) {
      //	  if (cumsum[i] < 0)
      //	    printf("i=%d %10g %10g\n", i, cumsum[i], cumsum[i - 1]);
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
  //d=0; printf("cumsum %d %10e;   %d %10e\n", d, cumsum[d], start_cumul-1, cumsum[start_cumul-1]);
  for (d=start_cumul; d<TMP; d++) {
    // printf("cumsum %d %10e\n", d, cumsum[d]);
    cumsum[d] += cumsum[d-1];
     //  printf("cumsum %d [%d,%d] c=%10e %10e asSign=%d (-1)=%d end=%10g inner=%10g ASS_I_OUT=%d\n", d, start_cumul, TMP, cumsum[d], cumsum[d]-cumsum[d-1], ASSIGNMENT(d), ASSIGNMENT(d-1), RIGHT_END(d), INNER, ASSIGN_IDX_OUTER);
     //
    assert(ASSIGNMENT(d)==ASSIGN_IDX_OUTER || ASSIGNMENT(d) >= ASSIGNMENT(d-1));
  }

  //APMI(cov);  
} // CumSum


// (unnormierte) VERTEILUNGSFUNKTION F(x,x,...,x) ist 
// c * (EXP(-O^p - EXP(-a x^p)) fuer grosse x, also integriert ueber 
// den Raum, gesehen als Funktion einer Variable x; O=OUTER
// Somit ist die Dichte gleich c * s * p * x^{p-d} * e^{-s x^p} / b
// wobei b die Oberflaeche des Wuerfels mit Kantenlaenge 2 bezeichnet
//
// weiter in  ExpVolOfCube
#define OUTER_DENSITY(v, x)					\
  if (OUTER_POW > 0) {						\
    double pow_x = OUTER_POW_CONST * POW(x, OUTER_POW);		\
    v = OUTER_CONST * OUTER_POW * pow_x * intpow(x, -dim) *	\
      EXP(-pow_x) / SURFACE(dim);				\
  } else {							\
    /* DICHTE ist c*x^p fuer grosse x */			\
    v = OUTER_CONST * POW(x, OUTER_POW);			\
  }


void evaluate_rectangular(double *x, model *cov, double *v) {  
  // *x wegen searchInverse !!
  RECTANGULAR_PARAMETERS;

  if (*x < 0.0) BUG; //crash();
  assert(*x >= 0);

  if (*x <= INNER) { // for technical reasons <= not <; see Dinverse
    // DICHTE ist c*x^p nahe Ursprung
    //  printf("inner ");
    *v = INNER_CONST * POW(*x, INNER_POW);
    //    printf("inner x=%10g %10g %10g\n", *x, INNER, *v);
    return;
  } else if (*x >= OUTER) {
    //    printf("outer ");
    if (next->finiterange == wahr) {
      *v = 0.0;
      return;
    }
    OUTER_DENSITY(*v, *x);
    return;					
  } else {
    //    printf("middle k=%d %10g %10g ", (int) ((*x - INNER) / STEP), VALUE(0),   VALUE( (int) ((*x - INNER) / STEP)));
   *v = VALUE( (int) ((*x - INNER) / STEP));
    return;
  }
}

void rectangularD(double *x, model *cov, double *v) {
  bool onesided = P0INT(RECT_ONESIDED); 
  if (onesided && *x <= 0) {
    *v = 0;
    return;
  }
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
  RECTANGULAR_PARAMETER_BASICS;			       
  int i;
   
  double max = RF_NEGINF;

  RECTFOR {double y = FABS(x[i]); if (y > max) max = y;}

  evaluate_rectangular(&max, cov, v);
  
  // printf("max %10g %10g\n", max, *v);
  // printf("rectD %10g %10g %d %10g\n", max, *v, P0INT(RECT_NORMED), OUTER_CUM);

  if (P0INT(RECT_NORMED)) *v /= OUTER_CUM; 
  if (onesided) *v *= 2.0; // R^1 !!
}


void rectangularDlog(double *x, model *cov, double *v) {
  rectangularD(x, cov, v);
  // printf("log rect v=%10g %10g\n", *v,  LOG(*v));
  *v = LOG(*v);
}


void rectangularDinverse(double *V, model *cov, double *left,
			 double *right) {
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
  RECTANGULAR_PARAMETERS;
  double  er, outer,
    x = RF_NA,
    v = *V; 
  int i;
  bool onesided = P0INT(RECT_ONESIDED); 

  if (P0INT(RECT_NORMED)) v *= OUTER_CUM; // unnormed v
  if (onesided) v /= 2.0;
  
  //  printf("V=%10e dim=%d\n", *V, dim);
    // assert(*V >= 0.0 && *V <= 1.0);
  
  if (*V <= 0.0) {
    RECTFOR {
      left[i] = RF_NEGINF;
      right[i] = RF_INF;
      //printf("i=%d %10g %10g\n", i, left[i], right[i]);
     }
     return;
  }

  
  ///  if (!next->finiterange == true && OUTER_POW > 1) {
  assert((next->finiterange == falsch) == (!next->finiterange == true));
  if (next->finiterange == falsch && OUTER_POW > 1) {
    // local maximum of density function exists for outer_pow > 1
    outer = POW( (OUTER_POW -1) / (OUTER_POW * OUTER_POW_CONST), 1 / OUTER_POW);
    if (outer < OUTER) outer = OUTER;    
  } else outer = OUTER;

  evaluate_rectangular(&outer, cov, &er); 

  //   printf("outer =%10g %10g %10g %10g\n", outer, OUTER, er, v);

  if (er > v) { // inverse is in the tail
    // first guess:
    double inverse,
      eps = 0.01;
     
    if (OUTER_POW > 0) {

      //      printf("here\n");
      // rough evaluation
      inverse = POW(- LOG(v / (OUTER_CONST * OUTER_POW)) / 
		    OUTER_POW_CONST, 1.0 / OUTER_POW);
      if (inverse <= outer) inverse = 2 * outer;
      x = searchInverse(evaluate_rectangular, cov, inverse, outer, v, eps);
    } else 
      //      printf("hereB\n");
      x = POW(OUTER_CONST / v, 1 / OUTER_POW);

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
	x = POW(v / INNER_CONST, 1.0 / INNER_POW);	
      } else BUG;
    }
  }
  
  RECTFOR {
    //printf("rect i=%d x=%10g\n", x);
    left[i] = onesided ? 0.0 : -x;
    right[i] = x;
  }
}


void rectangularP(double VARIABLE_IS_NOT_USED *x, 
		  model VARIABLE_IS_NOT_USED *cov, 
		  double VARIABLE_IS_NOT_USED *v) {
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
 // RECTANGULAR_PARAMETERS;

  NotProgrammedYet("");

}


void rectangularP2sided(double *x, double *y, model *cov, double *v) {
  // distribution functions
  bool onesided = P0INT(RECT_ONESIDED); 
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
  RECTANGULAR_PARAMETER_BASICS;	

  if (x != NULL) BUG;
  if (onesided && *y <= 0) {
    *v = 0.0;
    return;
  }

  //for (d=0; d<dim; d++) {// To Do ?? Stimmt das? Oder Masse der dimension d-k ?? SCJEINT NICHT ZU STIMMEN!!
  //    assert(y[d] >= 0.0);
  //    if (y[d] == 0.0) {
  //      *v = 0.0;
  //      return;
  //    }
  //  }

  
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
     
  //  CumSum(y, false, cov, TMP_WEIGHT);
  //  double tw = TMP_WEIGHT[TMP-1];

  CumSum(y, true, cov, TMP_WEIGHT);
 
  *v = TMP_WEIGHT[TMP-1];
  if (P0INT(RECT_NORMED)) *v /= OUTER_CUM;
  assert(R_FINITE(*v));
}


void rectangularQ(double *x, model VARIABLE_IS_NOT_USED *cov, 
		  double  VARIABLE_IS_NOT_USED *v) {
  if (*x < 0 || *x > 1) {*v = RF_NA; return;}
  if (!P0INT(RECT_APPROX)) ERR("approx=FALSE only for simulation");
  NotProgrammedYet("rectangularQ");
}


#define ACCEPT_REJECT							\
  if (!P0INT(RECT_APPROX)) {						\
  int ni = dim,								\
      quot = dim + 1,							\
      bytes =  sizeof(double) * dim;					\
    double newquot, approx, truevalue,					\
      max = RF_NEGINF;							\
    RECTFOR {double z = FABS(v[i]); if (z > max) max = z;}		\
    evaluate_rectangular(&max, cov, &approx);				\
    ABSFCTN(v, next, &truevalue);					\
    newquot = truevalue / approx;					\
    assert(quot < cov->qlen + 2);					\
    if (isMonotone(next->monotone)) { /* rejection sampling */		\
      assert(approx >= truevalue);					\
      cov->q[ni] = 0.0;							\
      if (UNIFORM_RANDOM >= newquot) {					\
	/*printf(".//");*/						\
	RESTART;							\
      } /* else printf("://"); */					\
    } else { /* MCMC */							\
      if (R_FINITE(cov->q[ni])) {					\
	cov->q[ni] -= - 1.0;						\
	if (UNIFORM_RANDOM * cov->q[quot] >= newquot) { /* keep old one */ \
	  MEMCOPYX(v, cov->q, bytes);					\
	} else { /* accept new one*/					\
	  cov->q[quot] = newquot; /* pi / Q */				\
	  MEMCOPYX(cov->q, v, bytes);					\
	}								\
      } else { /* !R_FINITE(cov->q[ni]), i.e. starting value */		\
	cov->q[ni] = P0INT(RECT_MCMC_N) - 1.0;				\
	cov->q[quot] = newquot; /* pi / Q */				\
	MEMCOPYX(cov->q, v, bytes);					\
      }									\
    }									\
    if (cov->q[ni] > 0) {						\
      RESTART;								\
    } else {								\
      cov->q[ni] = P0INT(RECT_MCMC_N);					\
      /* if (final) { */						\
      /* cov->total_n++;*/						\
      /* cov->total_sum += truevalue;*/					\
      /*} */								\
      return;								\
    }									\
  } else {								\
    if (final) {							\
      double approx,							\
	max = RF_NEGINF;						\
      RECTFOR {double z = FABS(v[i]); if (z > max) max = z;}		\
      evaluate_rectangular(&max, cov, &approx);				\
      /*      cov->total_n++;	*/					\
      /* cov->total_sum += approx;	*/				\
    }									\
  }



void rectangularR(double *x, model *cov, double *v) { 
  //printf("rectR %lu\n", x);

  
  if (x != NULL) {
   //crash();
    ERR("put 'flat = false'");
  }

  RECTANGULAR_PARAMETERS;
  bool final = true;
  
 EntryPoint_R:

  int i = CeilIndex(UNIFORM_RANDOM * OUTER_CUM, rect->weight, NSTEP + 2);
  //int i = searchFirstGreater(rect->weight, NSTEP + 2, UNIFORM_RANDOM  * OUTER_CUM);

  //  printf("i=%d dim=%d %d %d %.50s\n", i, dim, IDX_INNER, IDX_OUTER, NICK(cov->sub[0])); BUG;
  
 
  if (i == IDX_INNER) {
    RandomPointOnCubeSurface(POW(UNIFORM_RANDOM, 1.0 / (INNER_POW+dim)) * INNER,
			     dim, v);
    // printf("v=%10g\n", *v);
    // BUG;
  } else if (i == IDX_OUTER) {
     double u;
     assert(next->finiterange == falsch); //
    if (OUTER_POW > 0) // sampling by inversion formula
      u = POW(POW(OUTER, OUTER_POW) - LOG(UNIFORM_RANDOM) / OUTER_POW_CONST,
	      1 / OUTER_POW);
    else u = POW(UNIFORM_RANDOM, 1.0 / (OUTER_POW + dim)) * OUTER; // p+d < 0 !
    RandomPointOnCubeSurface(u, dim, v);
  
  } else { // i >= 0
   i -= IDX_STEPS;
    double start = INNER + i * STEP;
    RandomPointOnCubeRing(start, start + STEP, dim, v); 
  }
  if (P0INT(RECT_ONESIDED)) *v = FABS(*v);
#define RESTART goto EntryPoint_R
 
  ACCEPT_REJECT;
}

// sta tic int zi=0, zs=0, zo=0;

void rectangularR2sided(double *x, double *y, model *cov, double *v) { 
  if (x != NULL) 
    NotProgrammedYet("2-sided distribution function for rectangular");
  RECTANGULAR_PARAMETERS;
  int d, sel, red_dim, i,
    *idx = rect->idx;
  double *u, start, end,
    *ysort = rect->ysort; 

 EntrancePoint_R2:

  /*
=23772==    by 0x11DC440E: rectangularR2sided(double*, double*, model*, double*) (families.cc:2160)
==23772==    by 0x11CE9CB2: do_pgs_maxstable(model*, gen_storage*) (Huetchen.cc:973)
==23772==    by 0x11CEB402: do_pts_given_shape(model*, gen_storage*) (Huetchen.cc:1051)
==23772==    by 0x11D9618F: dompp(model*, gen_storage*, double*) (extremes.cc:318)
==23772==    by 0x11EC3693: simulate(double*, model*, double*) (rf_interfaces.cc:380)
==23772==    by 0x11F0CA7A: EvaluateModel (userinterfaces.cc:661)
  */
  
  CumSum(y, false, cov, TMP_WEIGHT); // ACHTUNG! reck->z gesetzt


  //assert(TMP_WEIGHT[TMP-1] >= 2.1 && TMP_WEIGHT[TMP-1] <= 2.2);
  assert(TMP > 0);
  double random = UNIFORM_RANDOM * TMP_WEIGHT[TMP-1];
  bool final = SQUEEZED(TMP-1) == 0 && 
    (!P0INT(RECT_APPROX) || next->randomkappa); // ohne Einschraenkung
  assert(!final); // nur zur Kontrolle

  //printf("random = %10g\n", random);
  sel = CeilIndex(random, TMP_WEIGHT, TMP);
  //printf("r=%10g inner=%10g o.tmp=%10g lst.step.tmp=%10g\n", random,   TMP_WEIGHT[IDX_INNER],   TMP_WEIGHT[TMP-1], TMP_WEIGHT[TMP-1]);

  red_dim = dim - SQUEEZED(sel);
  if (red_dim <= 0) BUG;
  start = sel > 0 ? RIGHT_END(sel - 1) : 0;
  end = RIGHT_END(sel);
  assert(start < end);

  // printf("%10g dim=%d squ=%d sel=%d, %10g %10g %d\n",
  //	 y[0], dim, SQUEEZED(sel), sel, random,  TMP_WEIGHT[TMP-1], TMP);
 
 
  u = TMP_WEIGHT; // nutzen des Platzes...
  //printf("assess %d %d %d\n", ASSIGNMENT(sel), sel, TMP);
  if (ASSIGNMENT(sel) == ASSIGN_IDX_INNER) {
    //   RandomPointOnCubeSurface(POW(UNIFORM_RANDOM, 1.0 / (INNER_POW+dim))
    //                        * INNER, dim, v);
    //zi ++;
    double // a + b x^p probab distr on [s, e], i.e. b = 1/(e^p-s^p),
      p = INNER_POW + red_dim,
      sp = POW(start, p),
      ep = POW(end, p),
      binv = ep - sp,
      a = - sp / binv,
      r = POW((UNIFORM_RANDOM - a) * binv, 1 / p);
    //
    //   printf("inner %d %d ap=%10g r=%10g inner=%10g w=%10g tmp=%10g, %10g\n", dim, SQUEEZED(sel), ap, r, INNER, INNER_CUM, TMP_WEIGHT[IDX_INNER], TMP_WEIGHT[IDX_INNER+1]);

    RandomPointOnCubeSurface(r, red_dim, u);
  } else if (ASSIGNMENT(sel) == ASSIGN_IDX_OUTER) {
    //zo++;
    // printf("outer\n");
    double w;
    assert(next->finiterange == falsch);//
    if (OUTER_POW > 0) {// sampling by inversion formula
      double op = POW(OUTER, OUTER_POW),
	factor = 1.0 - EXP(- OUTER_POW_CONST * (POW(end, OUTER_POW) - op));
      w = POW(op - LOG(1.0 - UNIFORM_RANDOM * factor) / OUTER_POW_CONST,
	      1 / OUTER_POW);
    } else w = POW(1.0 - 
		  UNIFORM_RANDOM *(1.0 - POW(end/OUTER, OUTER_POW + red_dim)),
		  1.0 / (OUTER_POW + red_dim)) * OUTER; 
		 // p+d < 0 ! 
    //printf("hXXere\n");
    RandomPointOnCubeSurface(w, red_dim, u);
  } else { // i \in steps
    //    zs++;
    //  printf("steps \n");
    //printf("sel=%d asS=%d %d\n",sel, ASSIGNMENT(sel), NSTEP);
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
  // idx[.] and ysort[.] have been set by CumSum
  assert(idx[0] == 0.0);
  for (d=1; d<=SQUEEZED(sel); d++) {
    v[idx[d] - 1] = (2.0 * UNIFORM_RANDOM - 1.0) * ysort[d];//idx[d] indiziert ab 1
  //                                                    da ysort[0]=0 per def.

 
 //  printf("%10g\n", ysort[d]);
    // v[idx[d] - 1] = 2 *  UNIFORM_RANDOM - 1; /* print */
  }
  for (d=SQUEEZED(sel); d<dim; d++) {
    v[idx[d+1] - 1] = u[d - SQUEEZED(sel)];
  }

  if (P0INT(RECT_ONESIDED)) *v = FABS(*v);

#undef RESTART
#define RESTART goto EntrancePoint_R2
   ACCEPT_REJECT;

  //printf("u=%10g %10g %10g\n", u[0], u[1], u[2]);
}


int GetMajorant(model *cov) {
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
    outer_factor = 1.3,
    in_pow = next->taylor[0][TaylorPow];

  bool
    inpownot0 = in_pow != 0.0;
  
  assert(!PisNULL(RECT_MAXSTEPS));

  int i, j, 
    max_steps = P0INT(RECT_MAXSTEPS),
    parts = P0INT(RECT_PARTS),
    max_iterations =  P0INT(RECT_MAXIT);


  assert(next->taylor[0][TaylorConst] > 0.0);
  INNER_CONST = next->taylor[0][TaylorConst] * safetyP1;
  INNER_POW = inpownot0 ? in_pow * EMsafety - dimsafety : 0.0;
  

  // INNER
  i = 0;
  INNER = starting_point;
  ABSFCTN(&(INNER), next, &v);
  m = INNER_CONST * POW(INNER, INNER_POW);
  //  printf("inner=%10e m=%10e v=%10e inner_c=%4.3f (%4.3f) inner_p=%4.3f (%4.3f)\n",  INNER, m, v, INNER_CONST, next->taylor[0][TaylorConst], INNER_POW, in_pow);

  while (m < v && INNER >= inner_min) {
    INNER *= inner_factor;
    INNER_CONST *= safetyP1;
    if (inpownot0) INNER_POW = INNER_POW * EMsafety - dimsafety;
    ABSFCTN(&(INNER), next, &v);
    m = INNER_CONST * POW(INNER, INNER_POW);
    // printf("inner=%10e m=%10e v=%10e i=%d, maxit=%d\n", 
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

    //printf("delta %10g x=%10g i=%d\n", delta, x, i);

    for (j=1; j<parts; j++, x-=delta) {
      //      printf("j=%d parts=%d max_it=%d\n", j, parts, max_iterations);

      ABSFCTN(&x, next, &v);
      m = INNER_CONST * POW(x, INNER_POW);
      double ratio = m / v;

      //printf("inner j=%d (%d) x=%4.3f m=%4.3f ratio=%4.3f old=%4.3f v=%10g\n", 
      //    j, max_iterations, x, m, ratio, old_ratio, v);

      if (!(  ratio >= old_ratio || (!inpownot0 && v <= INNER_CONST)  ))
	old_ratio = ratio; 
    }

    if (j == parts) break;
    INNER *= inner_factor;
    if (INNER > x) INNER = x;
    if (++i > max_iterations) 
      SERR2("%d iterations performed without success. Increase the value of '%.50s'", max_iterations, distr[RECT_MAXIT]);
  }

  //printf("inner j=%d (%d) inner=%4.3f m=%4.3f\n", 
      //  j, max_iterations, INNER, m);


  // OUTER
  if (next->tail[0][TaylorExpPow] <= 0.0) {
    // polynomial tail
    OUTER_POW = next->tail[0][TaylorPow] * safetyP1 + dimsafety;
    OUTER_POW_CONST = RF_NA;
  } else { // exponential tail
    assert(next->tail[0][TaylorExpConst] > 0.0);
    OUTER_POW = next->tail[0][TaylorExpPow] / safetyP1;
    OUTER_POW_CONST = next->tail[0][TaylorExpConst] / safetyP1;
  }
  OUTER_CONST = next->tail[0][TaylorConst] * safetyP1;
 
  if (next->finiterange == wahr) {
    // printf("xx %.50s %d\n", DefList[next->sub[0]->nr].nick, DefList[next->sub[0]->nr].finiterange == 1);
    
    // APMI(cov);
    INVERSE(ZERO(next), next, &(OUTER));
    if (INNER >= OUTER) OUTER = INNER * (1.0 + 1e-5);
  } else { // ! finiterange
    assert(next->tail[0][TaylorConst] > 0.0);
    i = 0;
    OUTER = starting_point * safetyP1;
    ABSFCTN(&(OUTER), next, &v);
    
    //    m = OUTER_CONST * (OUTER_POW <= 0 ? POW(OUTER, OUTER_POW) 
    //		       : EXP(- OUTER_POW_CONST *  POW(OUTER, OUTER_POW)));
    OUTER_DENSITY(m, OUTER);
    while (m < v && OUTER < outer_max) {
      if (++i > max_iterations) 
	SERR1("No majorant found. Function does not allow for a majorant or increase '%.50s'", distr[RECT_MAXIT]);
      OUTER_CONST *= safetyP1;
      OUTER_POW /= safetyP1;
      OUTER *= outer_factor;
      ABSFCTN(&(OUTER), next, &v);
      OUTER_DENSITY(m, OUTER);
      // printf("i =%d  outer %10g m=%10e v=%10e\n", i, OUTER, m, v); assert(false);
  }
    if (OUTER > outer_max && m < v)
      SERR("No majorant found close for large arguments");
    
    //  printf("i =%d  outer %10g m=%10e v=%10e\n", i, OUTER, m, v); //assert(false);

    //printf("i =%d  outer %10g %10g m=%10g v=%10g\n", 
    //i, OUTER, 
    //	  OUTER_CONST * (OUTER_POW <= 0 ? POW(OUTER, OUTER_POW) 
    //			 : EXP(- OUTER_POW_CONST *  POW(OUTER, OUTER_POW))),
    //	  m , v); //assert(false);

    i=0;
    while (true) {
      old_ratio = m / v;
      delta = OUTER / parts;
      double x = OUTER + delta;
      for (j=1; j<parts; j++, x+=delta) {
	ABSFCTN(&x, next, &v);
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

  //  printf("inner outer %10g %10g\n", INNER, OUTER);// assert(false);
  // STEPS  
  
  // todo variable Schrittweiten waeren wegen des steilen
  // Anstiegs am Pol besser
  STEP = (OUTER - INNER) / max_steps;
  if (STEP < min_steps) {
    NSTEP = (int) ((OUTER - INNER) / min_steps);
    STEP = (OUTER - INNER) / NSTEP;
  } else NSTEP = max_steps;

  // PMI(cov, "getm");

  RETURN_NOERROR; 
}

int check_rectangular(model *cov) {
  model *next = cov->sub[0];
  int err,
    dim = OWNXDIM(0);

  ASSERT_CARTESIAN;

  kdefault(cov, RECT_SAFETY, GLOBAL.distr.safety);
  kdefault(cov, RECT_MINSTEPLENGTH, GLOBAL.distr.minsteplen);
  kdefault(cov, RECT_MAXSTEPS, GLOBAL.distr.maxsteps);
  kdefault(cov, RECT_PARTS, GLOBAL.distr.parts);
  kdefault(cov, RECT_MAXIT, GLOBAL.distr.maxit);
  kdefault(cov, RECT_INNERMIN, GLOBAL.distr.innermin);
  kdefault(cov, RECT_OUTERMAX, GLOBAL.distr.outermax);
  kdefault(cov, RECT_MCMC_N, GLOBAL.distr.mcmc_n);
  kdefault(cov, RECT_NORMED, true); // currently only the normed version 
  // is needed by RandomFields. Delete the parameter?? 2.2.19
  kdefault(cov, RECT_APPROX, true);
  kdefault(cov, RECT_ONESIDED, false);
  //  int onesided = P0INT(RECT_ONESIDED);
  
  if (cov->q == NULL) QALLOC(dim + 2);
  cov->q[dim] = RF_NA;

  if ((err = CHECK(next, dim, dim, ShapeType, XONLY, 
		   dim == 1 && P0INT(RECT_ONESIDED) ? CARTESIAN_COORD 
		   : ISOTROPIC, 
		   SCALAR, cov->frame)) != NOERROR) {
    RETURN_ERR(err);
  }

  if (next->randomkappa) RETURN_ERR(ERRORRANDOMKAPPA);
  
  // if (!isMonotone(next->monotone)) SERR("only monotone submodels are allowed");
  

  if (next->taylorN <= 0 || next->tailN <= 0) {
    //    APMI(next);
    SERR1("'%.50s' does not have integrability information", NICK(next));
  }
  assert(R_FINITE(next->taylor[0][TaylorPow]));
  if (next->taylor[0][TaylorPow] <= -dim)
    SERR1("pole of '%.50s' not integrable", NICK(next));
  assert(R_FINITE(next->tail[0][TaylorPow]));

  //  PMI(next);

  if (next->tail[0][TaylorPow] >= -dim && next->tail[0][TaylorExpPow] == 0.0 &&
      next->tail[0][TaylorConst] != 0.0) 
    SERR1("tail of '%.50s' not integrable", NICK(next));
  assert(R_FINITE(next->tail[0][TaylorConst]));

  //PMI(next);
  if (next->taylor[0][TaylorConst] == 0.0) {
    //PMI(next);
    SERR1("'%.50s' seems to be a trivial shape function", NICK(next));
  }

  VDIM0 = OWNLOGDIM(0);
  VDIM1 = 1;

  //PMI(cov);

  RETURN_NOERROR;

}



int init_rectangular(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  int  err;

  NEW_STORAGE(rect);
  RECTANGULAR_PARAMETERS; // muss nach obigen stehen und vor allem anderen

  if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR) {
    //PMI(cov);
    RETURN_ERR(err);
  }
  
  if ((err = GetMajorant(cov)) != NOERROR) RETURN_ERR(err); // muss genau nach MALLOC oben und MALLOC unten stehen
  if (INNER >= OUTER) BUG;

  //   APMI(cov);
   
  int d,
    abschnitte = NSTEP + 2, 
    tmp_n = NSTEP + 2 + dim;
  double x;

  if ((rect->value = (double*) MALLOC(sizeof(double) * abschnitte)) == NULL ||
      (rect->weight = (double*) MALLOC(sizeof(double) * abschnitte)) == NULL||
      (TMP_WEIGHT = (double*) CALLOC(tmp_n, sizeof(double))) == NULL ||
      (rect->right_endpoint = (double*) MALLOC(sizeof(double) * tmp_n)) ==NULL||
      (rect->ysort = (double*) MALLOC(sizeof(double) * (dim + 1))) == NULL ||
      (rect->z = (double*) MALLOC(sizeof(double) * (dim + 1))) == NULL ||//dummy
      (rect->squeezed_dim = (int*) MALLOC(sizeof(int) * tmp_n)) == NULL ||
      (rect->asSign = (int*) MALLOC(sizeof(int) * tmp_n)) == NULL ||
      (rect->idx = (int*) MALLOC(sizeof(int) * (dim + 1))) == NULL)
    RETURN_ERR(ERRORMEMORYALLOCATION);

 
  x = INNER;
  for (d=0; d<NSTEP; d++, x+=STEP) {
    ABSFCTN(&x, next, &(VALUE(d)));
    //printf("x=%10g %10g\n", x, VALUE(d));
    assert(VALUE(d) >= 0);
  }
  rect->value[IDX_INNER] = rect->value[IDX_OUTER] = RF_NA;
  EXTRA_STORAGE;
  for(d=0; d<dim; d++)  TMP_WEIGHT[d] = RF_INF;
  CumSum(TMP_WEIGHT, false, cov, rect->weight);

  cov->mpp.mM[0] = cov->mpp.mMplus[0] = P0INT(RECT_NORMED) ? 1.0 : OUTER_CUM;
  assert(R_FINITE(OUTER_CUM));
  if (cov->mpp.moments > 0) {
    double f = cov->mpp.mM[0] / next->mpp.mM[0];
    cov->mpp.mM[1] = next->mpp.mM[1] * f;
    cov->mpp.mMplus[1] = next->mpp.mMplus[1] * f; 
    if (!R_FINITE(cov->mpp.mM[1])){
      //      crash();      APMI(cov);
      BUG;
    }
  }

  // normed against cov->sub[0]
  cov->mpp.maxheights[0] = RF_NA;
  if (isMonotone(next->monotone))
    cov->mpp.maxheights[0] = INNER_POW >= 0 ? INNER_CONST : RF_INF;
  cov->mpp.unnormedmass = OUTER_CUM;
  

  // if (OUTER_CUM < 1.0) { APMI(cov); }

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
	ABSFCTN(&xx, next, &w);
	double y = xx + 1e-10;
	assert(y >= 0.0);
	evaluate_rectangular(&y, cov, &e);
	//printf("%10g %10e %10e\n", xx, w, e);
	u += w * 4 * M_PI / 3 * (powl(xx + mini, 3) - powl(xx, 3));
	t += w * (powl(2 * (xx + mini), 3) - powl(2 * xx, 3));
	xx += mini;
      }
      tot_u += u;
      tot_t += t;
      //      printf("%d %10e %10e tot=%10g %10g\n", k-1, u, t, tot_u, tot_t);
    }
  }

  //  printf("totalweight=%10g\n", OUTER_CUM);
  assert(OUTER_CUM >= 1.0 || NEXTNR != STROKORB_MONO); 
  //
  //  printf("init rect %d %d\n", NSTEP, TMP);
  assert(NSTEP == TMP - 2);
  

  //cov->total_n = 0;
  // cov->total_sum = 0.0;

   RETURN_NOERROR;
}


void do_rectangular(model *cov, double *v){
  model  *next = cov->sub[0];
  //int err;
  gen_storage s;
  gen_NULL(&s);  
  DO(next, &s);  
  rectangularR(NULL, cov, v);
}

void range_rectangular(model *cov, range_type *range){
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

  booleanRange(RECT_NORMED);
  booleanRange(RECT_APPROX);

  range->pmin[RECT_ONESIDED] = range->min[RECT_ONESIDED] = 0;
  range->pmax[RECT_ONESIDED] = range->max[RECT_ONESIDED] = OWNLOGDIM(0)==1;
  range->openmin[RECT_ONESIDED] = false;
  range->openmax[RECT_ONESIDED] = false;
}




#define MCMC_FCTN 0

#define MCMC_MCMC_N 0
#define MCMC_MCMC_SIGMA 1
#define MCMC_NORMED 2
#define MCMC_MAXDENSITY 3
#define MCMC_RANDLOC 4
#define MCMC_GIBBS 5
#define MCMC_LAST_PARAM MCMC_GIBBS

void kappa_mcmc(int i, model  VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  *nc = 1;
  *nr = (i == MCMC_MCMC_SIGMA) ? 0 : (i <= MCMC_LAST_PARAM) ? 1 : -1;
}

void mcmcIntegral(model  VARIABLE_IS_NOT_USED *cov) {
  NotProgrammedYet("mcmcIntegral");
  /*
    int n = ;
  bool normed = P0INT(MCMC_NORMED);
  double value;
  PINT(MCMC_NORMED)[0] = false;
  while (cov->q[MCMC_CURR_N] >= 1.0) {
    cov->q[MCMC_CURR_SUM] -= 1.0;
    mcmcR(NULL, unklar welche Funktion, &value);    
  }
  */
}


void mcmcD(double *x, model *cov, double *v) {
  model *next = cov->sub[MCMC_FCTN];
  //  int mcmc_n =  P0INT(MCMC_MCMC_N);
  FCTN(x, next, v);
  *v = FABS(*v);
  if (P0INT(MCMC_NORMED)) {
    BUG; // not programmed yet
  }
}

void mcmcDlog(double *x, model *cov, double *v) {
  mcmcD(x, cov, v);
  *v = LOG(*v);
}

void mcmcDinverse(double  VARIABLE_IS_NOT_USED *V, model  VARIABLE_IS_NOT_USED *cov, double  VARIABLE_IS_NOT_USED *left,
			 double  VARIABLE_IS_NOT_USED *right) {
  NotProgrammedYet("mcmcDinverse");
}


void mcmcP(double VARIABLE_IS_NOT_USED *x, 
		  model VARIABLE_IS_NOT_USED *cov, 
		  double VARIABLE_IS_NOT_USED *v) {
 NotProgrammedYet("mcmcP");
}

void mcmcP2sided(double  VARIABLE_IS_NOT_USED *x, double  VARIABLE_IS_NOT_USED *y, model VARIABLE_IS_NOT_USED  *cov, double VARIABLE_IS_NOT_USED  *v) {
  NotProgrammedYet("mcmcP2sided");
}


void mcmcQ(double *x, model VARIABLE_IS_NOT_USED *cov, 
		  double  VARIABLE_IS_NOT_USED *v) {
  if (*x < 0 || *x > 1) {*v = RF_NA; return;}
  NotProgrammedYet("mcmcQ");
}

void mcmcR(double *x, model *cov, double *v) { 
  if (x != NULL) {
    ERR("put 'flat = false'");
  }

  model *next = cov->sub[MCMC_FCTN];
  location_type *loc=Loc(cov);
  int i, d,
    vdim = total_logicaldim(OWN),
    n = P0INT(MCMC_MCMC_N);
  double proposedvalue, 
    maxdens = P0(MCMC_MAXDENSITY),
    *sigma = P(MCMC_MCMC_SIGMA),
    posvalue = cov->Smcmc->posvalue,
    *delta = cov->Smcmc->deltapos,
    *pos = cov->Smcmc->pos;
  assert(delta != NULL && pos != NULL);
  bool
    gibbs = (bool) P0INT(MCMC_GIBBS),
    randloc = (bool) P0INT(MCMC_RANDLOC);
  TALLOC_X1(proposed, vdim);
  TALLOC_X2(propdelta, vdim);

 
  for (i=0; i<n; i++) {
    for (d = 0; d < vdim; d++) propdelta[d] = delta[d];
    if (gibbs) {
      int idx = (int) (UNIFORM_RANDOM * (double) vdim);
      proposed[idx] = (propdelta[idx] += 
			rnorm(0.0, sigma[idx % cov->nrow[MCMC_MCMC_SIGMA]]));
    } else {
      for (d=0; d < vdim; d++) {
	//	printf("%d %d %10g %10g\n", d, d % cov->nrow[MCMC_MCMC_SIGMA], 
	//	       sigma[d % cov->nrow[MCMC_MCMC_SIGMA]], proposed[d]);	
	proposed[d] = (propdelta[d] 
		       += rnorm(0.0, sigma[d % cov->nrow[MCMC_MCMC_SIGMA]]));
      }
    }
    
    if (loc != NULL && randloc) { // to do
      if (loc->grid) {  
	for (d = 0; d < vdim; d++) {   
	  proposed[d] += loc->xgr[d][XSTART] + loc->xgr[d][XSTEP] *
	    (double) ((int) UNIFORM_RANDOM * (loc->xgr[d][XLENGTH] - 1.0));
	}
      } else {
	double
	  *xx = loc->x + vdim * (int) (loc->spatialtotalpoints*UNIFORM_RANDOM);
	if (loc->Time) {
	  int vM1 = vdim - 1;
	  for (d = 0; d < vM1; d++) proposed[d] += xx[d];
	  proposed[d] += loc->T[XSTART] + loc->T[XSTEP] *
	    (double) ((int) UNIFORM_RANDOM * (loc->T[XLENGTH] - 1.0));
	} else {
	  for (d = 0; d < vdim; d++) proposed[d] += xx[d];
	}
      }         
    }

 
    FCTN(proposed, next, &proposedvalue);
    
    if (proposedvalue > maxdens) { proposedvalue = maxdens; }
    if (proposedvalue > posvalue || proposedvalue > posvalue*UNIFORM_RANDOM) {
      posvalue = proposedvalue;
      for (d = 0; d < vdim; d++) {
	pos[d] = proposed[d];
	delta[d] = propdelta[d];
      }
    }
  } // for i<n
  END_TALLOC_X1;
  END_TALLOC_X2;
  
  cov->Smcmc->posvalue = posvalue;
  for (d = 0; d < vdim; d++) v[d] = pos[d];
}

void mcmcR2sided(double VARIABLE_IS_NOT_USED  *x, double  VARIABLE_IS_NOT_USED *y, model  VARIABLE_IS_NOT_USED *cov, double  VARIABLE_IS_NOT_USED *v) { 
  
}

int check_mcmc(model *cov) {
  ASSERT_UNREDUCED;
  ASSERT_CARTESIAN;
  model *next = cov->sub[MCMC_FCTN];
  int err, vdim; 

  ASSERT_CARTESIAN;

  kdefault(cov, MCMC_NORMED, false);
  if (P0INT(MCMC_NORMED)) NotProgrammedYet("mcmc (normed=TRUE)"); // OK

  vdim = total_logicaldim(OWN);
  if (vdim != OWNTOTALXDIM) SERR("inconsistent dimensions given.");
  if ((err = CHECK(next, vdim, vdim, ShapeType, XONLY, CARTESIAN_COORD,
		   SCALAR, RandomType)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  VDIM0 = vdim;
  VDIM1 = 1;
  
  if (PisNULL(MCMC_MCMC_SIGMA)) {
    location_type  *loc = Loc(next);
    if (loc != NULL && loc->grid) {
      PALLOC(MCMC_MCMC_SIGMA, vdim, 1);
      for (int d=0; d<vdim; d++)
	P(MCMC_MCMC_SIGMA)[d] = loc->xgr[d][XSTEP] * 0.1;
    } else {
      SERR1("'%.50s' must be given.", KNAME(MCMC_MCMC_SIGMA));
    }
  }

  kdefault(cov, MCMC_MCMC_N, GLOBAL.distr.mcmc_n);
  kdefault(cov, MCMC_MAXDENSITY, 1000);
  kdefault(cov, MCMC_RANDLOC, false); // scattering among 
  //                                the locations if on a grid?
  kdefault(cov, MCMC_GIBBS, false);

  NEW_STORAGE(mcmc);
  EXTRA_STORAGE;

  RETURN_NOERROR;
}



int init_mcmc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  model *next = cov->sub[MCMC_FCTN];
  location_type *loc=Loc(cov);
  int  d, err,
    vdim = total_logicaldim(OWN);
   double 
    maxdens = P0(MCMC_MAXDENSITY);
   //#ifdef DO_PARALLEL
     // if (GLOBAL_UTILS->basic.cores != 1)
   //   SERR("'RRmcmc' can be used only with 1 core (RFoptions(cores=1). Use multicores from within R.")       
   //#endif   
  if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR) RETURN_ERR(err);
  ALLC_NEW(Smcmc, pos, vdim, pos);
  ALLC_NEW(Smcmc, delta, vdim, deltapos);
  
  //PMI(cov);
  for (d=0; d<vdim; d++) pos[d] = delta[d] = 0.0;
  if (loc != NULL && loc->lx > 0) { 
    if (loc->grid) {
      assert(loc->xgr[0]!= NULL);
      for (d=0; d<vdim; d++) pos[d] = loc->xgr[d][XSTART];
    } else if (loc->Time) {
      int vM1 = vdim - 1;
      assert(loc->x != NULL);
      for (d=0; d<vM1; d++) pos[d] = loc->x[d];
      pos[d] = loc->T[XSTART];
    } else {
      assert(loc->x != NULL);
      for (d=0; d<vdim; d++) pos[d] = loc->x[d];
    }
  }

  FCTN(pos, next, &(cov->Smcmc->posvalue));

  if (cov->Smcmc->posvalue > maxdens) cov->Smcmc->posvalue = maxdens;
  RETURN_NOERROR;
}


void do_mcmc(model *cov, double *v){
  model *next = cov->sub[MCMC_FCTN];//int err;
  gen_storage s;
  gen_NULL(&s);  
  DO(next, &s);  
  mcmcR(NULL, cov, v);
}


void range_mcmc(model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[MCMC_MCMC_N] = 1;
  range->max[MCMC_MCMC_N] = RF_INF;
  range->pmin[MCMC_MCMC_N] = 20;
  range->pmax[MCMC_MCMC_N] = 50;
  range->openmin[MCMC_MCMC_N] = false;
  range->openmax[MCMC_MCMC_N] = true;
  
  range->min[MCMC_MCMC_SIGMA] = 0;
  range->max[MCMC_MCMC_SIGMA] = RF_INF;
  range->pmin[MCMC_MCMC_SIGMA] = 1e-7;
  range->pmax[MCMC_MCMC_SIGMA] = 1e7;
  range->openmin[MCMC_MCMC_SIGMA] = true;
  range->openmax[MCMC_MCMC_SIGMA] = true;

  booleanRange(MCMC_NORMED);

  range->min[MCMC_MAXDENSITY] = 0;
  range->max[MCMC_MAXDENSITY] = RF_INF;
  range->pmin[MCMC_MAXDENSITY] = 1e-7;
  range->pmax[MCMC_MAXDENSITY] = 1e7;
  range->openmin[MCMC_MAXDENSITY] = true;
  range->openmax[MCMC_MAXDENSITY] = true;

  booleanRange(MCMC_RANDLOC);
  booleanRange(MCMC_GIBBS);

 }




// #define UNIF_LOGREFAREA dim
#define SIGN_P 0
void randomSign(double *x, model *cov, double *v) { 
  model *next = cov->sub[0];
  assert(next != NULL);
  COV(x, next, v);
  (*v) *= cov->q[0];
  // printf("%10g ", cov->q[0]);
}

void lograndomSign(double *x, model *cov, double *v, double *Sign) { 
  model *next = cov->sub[0];
  assert(next != NULL);
  LOGCOV(x, next, v, Sign);
  (*Sign) *= cov->q[0];
}

void randomSignInverse(double *v, model *cov, double *x){
  model *next = cov->sub[0];
  INVERSE(v, next, x);
}
void randomSignNonstatInverse(double *v, model *cov, double *x, double *y){
  model *next = cov->sub[0];
  NONSTATINVERSE(v, next, x, y);
}

int check_randomSign(model *cov) {
  model *next = cov->sub[0];
  int err, 
    size = 1;
  if (cov->q == NULL) QALLOC(size);
  
  kdefault(cov, 0, 0.5);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  if ((err = CHECK_NOPASS(next)) != NOERROR)  RETURN_ERR(err);
  setbackward(cov, next);
  RETURN_NOERROR;
}

void range_randomSign(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SIGN_P] = 0; 
  range->max[SIGN_P] = 1;
  range->pmin[SIGN_P] = 0.0001;
  range->pmax[SIGN_P] = 0.9999;
  range->openmin[SIGN_P] = false;
  range->openmax[SIGN_P] = false;
}

int init_randomSign(model *cov, gen_storage *s) {
  model *next = cov->sub[0];
  int err;
  double Eminus;
  assert(next != NULL);
  if ((err=INIT(next, cov->mpp.moments, s))!= NOERROR) RETURN_ERR(err);
  if (next->fieldreturn == wahr && next->loggiven) 
    SERR("log return is incompatible with random Sign");
  
  if (cov->mpp.moments >= 1) {
    cov->mpp.mM[0] = next->mpp.mM[0];
    cov->mpp.mMplus[0] = next->mpp.mMplus[0];
    Eminus = cov->mpp.mMplus[1] - cov->mpp.mM[1];
    cov->mpp.mMplus[1] = 
      P(SIGN_P)[0] * (cov->mpp.mMplus[1] - Eminus) + Eminus; 
    cov->mpp.mM[1] = 0.0;
  }
  cov->mpp.maxheights[0] = next->mpp.maxheights[0];
  cov->mpp.unnormedmass = next->mpp.unnormedmass;
  ReturnOtherField(cov, next);

  // assert(cov->mpp.maxheights[0] == 1.00);

  RETURN_ERR(err);
}

void do_randomSign(model *cov, gen_storage *s) {
  model *next = cov->sub[0];
  assert(next != NULL);
  DO(next, s); // nicht gatternr
  cov->q[0] = 2.0 * (UNIFORM_RANDOM <= P0(SIGN_P)) - 1.0;
  if (cov->q[0] != 1.0 && next->fieldreturn == wahr) { 
    assert(cov->q[0] == - 1.0);
    if (next->loggiven) ERR("log return is incompatible with random Sign");
    int i,
      endfor = Loc(next)->totalpoints;
    double *rf = cov->rf;
    for (i=0; i<endfor; i++) rf[i] = -rf[i];
  }
}


int struct_randomSign(model *cov, model **newmodel) {  
  if (hasGaussMethodFrame(cov) || hasPoissonFrame(cov)) {
    int err = STRUCT(cov->sub[0], newmodel);
    //  assert(cov->sub[0]->mpp.maxheights[0] == 1.0);
    RETURN_ERR(err);
  }
  SERR1("'%.50s' not allowed in this context.", NICK(cov));
}
