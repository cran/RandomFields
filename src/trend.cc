
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Handling of the different possibilities to pass the trend

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nondomain models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc

 Copyright (C) 2011 -- 2015 Marco Oesting & Martin Schlather

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

/*Note:
the parameter 'polycoeff' is hidden in R; it is a vector consisting of the
coefficients of the corresponding trend polynomials:
the first choose(polydeg[1]+d,d) components belong to the first polynomial,
the next choose(polydeg[2]+d,d) to the second one,...
the corresponding monomial functions are of the following order:
1, x, x^2, ..., x^k, y, x y, ... x^(k-1) y, y^2, x y^2..., y^k,
z, x z, x^2 z, ...., x^(k-1) z, y z, x y z, x^2 y z, ..., z^k
*/


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <math.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

#include "RF.h"
#include "primitive.h"
#include "Operator.h"
#include "variogramAndCo.h"



//////////////////////////////////////////////////////////////////////
//    mixed
//////////////////////////////////////////////////////////////////////

#define MIXED_CONSTANT 0

void mixed(double *x, cov_model *cov, double *v) {
  // bis auf den letzten Schluss alles schon auf bel vdim geschrieben.
  cov_model  *next = cov->sub[0];
  location_type *loc = Loc(cov);
  mixed_storage *s = cov->Smixed;
  double //distTF, 
    var = 1.0, 
    *covmatrix=NULL;
  int i, err, // Xnrow, 
    //    err=NOERROR,
    vdim = cov->vdim[0],
    vdimsq = vdim * vdim;
  if (cov->vdim[0] != cov->vdim[1]) BUG;

  if (isNegDef(cov) && cov->nsub == 0) {
    for (i=0; i<vdimsq; i++) v[i] = 0.0;
    return;
  }


  if (!isShape(cov) && !isTrend(cov)) BUG;
  *v = LP0(MIXED_X) * P0(MIXED_BETA);
  return; 

  // TO DO

  BUG;

  if (PisNULL(MIXED_X)) { // no X given, but a submodel	
    COV(x, next, v);
    BUG;
    return;
  }
    
 

  if (cov->q[MIXED_CONSTANT]) {
    // Xu , u ~ N(0, C)
    // classical random effect
    cov_model *sub = next;
    while (isDollar(sub)) {
      var *= PARAM0(sub, DVAR);
      sub=sub->sub[0];
    }
    covmatrix = LPARAM(next, 99999); BUG; // FIX_M);
    if (LNCOL(MIXED_X) != LPARAMNROW(next, 99999)) { // FIX_M)) {
      GERR("mixed model: X and covariance matrix M do not match");
    }
  } else {    
    // Xu , u ~ random field
    // random effect, geostatistical covariance model
    CovarianceMatrix(cov->key, s->mixedcov);
    covmatrix = s->mixedcov;
  } // sub->nr != FIX

  // ab hier vdim =1  
  *v = XkCXtl(LP(MIXED_X), covmatrix, LNROW(MIXED_X), LNCOL(MIXED_X),
	      loc->i_row, loc->i_col);
  *v *= var;
  return;

 ErrorHandling:
  XERR(err);
  
}

void mixed_nonstat(double *x, double *y, cov_model *cov, double *v){

  NotProgrammedYet("");
  

  if (cov->nsub != 0 && PisNULL(MIXED_X)) { // no X given, but a submodel	
    cov_model *next = cov->sub[0];
    NONSTATCOV(x, y, next, v);
  } else mixed(x, cov, v);
}

void covmatrix_mixed(cov_model *cov, double *v) {
  cov_model *sub = cov->sub[0];
  
  if (cov->ncol[MIXED_X] == 0) { // || element < 0) {
    CovList[sub->nr].covmatrix(sub, v);
    return;
  }

  double *C=NULL;
  int nz=0,
    nrow = LNROW(MIXED_X), 
    ncol= LNCOL(MIXED_X),
    ncol2 = ncol * ncol; 
  C = (double*) MALLOC(sizeof(double) * ncol2);
  if (C==NULL) {
    StandardCovMatrix(cov, v);
    return;
  }
  
  CovList[sub->nr].covmatrix(sub, C);
  XCXt(LP(MIXED_X), C, v, nrow, ncol); BUG;//todo:?reprogramm XCXt with alloc here ?
  Loc(cov)->totalpoints = nrow;
  for (int j=0; j<ncol2; j++) nz += C[j] != 0.0;

  // PMI(cov);
  
  /*
    long i,j,k, tot=Loc(cov)->totalpoints;
  printf("\nStart mixed C %d  ---- check within comment\n", element); //
  for (k=i=0; i<ncol * ncol; i+=tot) {
     for (j=0; j<ncol; j++) printf("%f ", C[k++]); //
     printf("\n");//
  }

  printf("\nStart mixed t(X) %d\n", element); //
  for (k=i=0; i<tot*tot; i+=tot) {
  for (j=0; j<tot; j++) printf("%f ", X->p[element][k++]);//
   printf("\n");//
  }

  printf("\nStart mixed t(v) %d\n", element); //
    for (k=i=0; i<tot*tot; i+=tot) {
    for (j=0; j<tot; j++) printf("%f ", v[k++]);//
    printf("--- end check within comment \n");//
  }
  */
  
  FREE(C);
}

int set_mixed_constant(cov_model *cov) {
  cov_model 
    *next = cov->sub[0],
    *sub = next;   
  bool simple = true;
  location_type *loc=Loc(cov);
  long i,
    totalpoints = loc->totalpoints;
  int
    *ncol = cov->ncol,
    *nrow = cov->nrow; 
  listoftype *X = PLIST(MIXED_X);

  if (cov->q == NULL) QALLOC(1);

  while (sub != NULL && isDollar(sub) &&
	 ((simple = PARAMisNULL(sub, DPROJ)
	   && (PARAMisNULL(sub, DSCALE) || PARAM0(sub, DSCALE) == 1.0)
	   && PARAMisNULL(sub, DANISO)))) sub=sub->sub[0];
  
  BUG;
  if ((cov->q[MIXED_CONSTANT] = sub != NULL && 9999)) {// sub->nr == FIX)) {
    //next->delflag = DEL_COV - 6;
    if (isDollar(next) && next->nrow[DVAR]==0) {
      //next->delflag = DEL_COV - 6;
      if (!simple) 
	SERR1("'%s' not allowed together with an anisotropic structrue",
	     NICK(cov));
    }
    
    int constant_size;
    
    for (i=0; i<nrow[MIXED_X]; i++) {
      BUG;
      constant_size = PARAMLIST(sub, 9999) ///FIX_M)
	->nrow[i];
      
      if (ncol[MIXED_X] > 0 && X->ncol[i] != constant_size) {
	SERR7("%ldth matrix '%s' (%d x %d) and (%d x %d) constant matrix '%s' do not match", i, KNAME(MIXED_X), X->nrow[i], X->ncol[i], constant_size, constant_size, NICK(sub));
      }
    }
  } else {
    for (i=0; i<nrow[MIXED_X]; i++) {
      if (X->nrow[i] != X->ncol[i])
	SERR3("%ldth  matrix is not symmetric (%d x %d)",
	      i+1, X->nrow[i], X->ncol[i]);
    }
  }
  
  if (false) // test ist wichtig, kollidiert z.Zt. aber mit SetAndGetModelInfo
    for (i=0; i<nrow[MIXED_X]; i++) {
      if (X->nrow[i] != totalpoints) {
	//PMI(cov);
	SERR4("number of rows of entry %ld of '%s' (%d) are different from the number of locations (%ld)", i+1, KNAME(MIXED_X),  X->nrow[i], totalpoints);
      }
    }
  return NOERROR;
}

char iscovmatrix_mixed(cov_model *cov) {
  int err;
  // printf("iscov %d %ld\n", cov->qlen, cov->q);
  if (cov->q == NULL && (err = set_mixed_constant(cov)) != NOERROR) XERR(err);
  // printf("iscov %d %ld\n", cov->qlen, cov->q);
  return 2 * (int) (cov->nsub > 0) * (int) (cov->q[MIXED_CONSTANT] || 
					    cov->ncol[MIXED_X] > 0);
}

void kappamixed(int i, cov_model  VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  if (i==MIXED_DIM) *nc = *nr = 1;
  else if (i==MIXED_BETA || i==MIXED_DIST) {
    *nc = 1;
    *nr = SIZE_NOT_DETERMINED; 
  } else if (i==MIXED_X || i==MIXED_COORD) *nc = *nr = SIZE_NOT_DETERMINED; 
  else *nc = *nr = -1;
}
int checkmixed(cov_model *cov) {

 
    // cov_model *sub;
  //location_type *loc=Loc(cov);
  int i, err, 
    nkappa = CovList[cov->nr].kappas,
     nsub = cov->nsub,
    *ncol = cov->ncol,
    *nrow = cov->nrow; // taken[MAX DIM],
  char msg[300];

  //printf("%d\n", MIXED_X);
  listoftype *X = PLIST(MIXED_X);

  if (!isCartesian(cov->isoown)) return ERRORNOTCARTESIAN;

  if (cov->nrow[0] * cov->ncol[0] != 1 || cov->nrow[1] * cov->ncol[1] != 1)
    SERR("currently only constant mean possible");

  ROLE_ASSERT(ROLE_GAUSS || cov->role == ROLE_COV);
  for (i=0; i<nkappa; i++) {// NEVER; cf kriging.R, for instance
    if (cov->kappasub[i] != NULL) SERR("parameters may not be random");
    if (cov->ncol[i] * cov->nrow[i] > (i < 2))
      SERR("current version does not allow for mxied effects as this part of the package is being rewritten."); // SELECT ebenso loeschen?!
  }

  cov->vdim[0] = cov->vdim[1] = 1; //falls kein submodel vorhanden (Marco)
  cov->maxdim=INFDIM;
  cov->matrix_indep_of_x = true;

  if (ncol[MIXED_BETA] > 0) { // b is given  
    if (nsub != 0) 
      SERR("b and a covariance model may not be given at the same time");
    if (ncol[MIXED_X] == 0) SERR2("if '%s' is given '%s' must be given",
				  KNAME(MIXED_BETA), KNAME(MIXED_X));
    for (i=0; i<nrow[MIXED_X]; i++) {
      if (X->ncol[i] != nrow[MIXED_BETA]) {
	sprintf(msg,
		"%dth set: (%d x %d) matrix '%s' and (%d x %d) vector '%s' do not match",
		i, X->nrow[0], X->ncol[i], KNAME(MIXED_X),
		nrow[MIXED_BETA], ncol[MIXED_BETA], KNAME(MIXED_BETA));
	SERR(msg);
      }
    }
  } else if (nsub == 0) { // only X is given -- then only a deterministic 
	//                                 constant is given
    if (ncol[MIXED_BETA] == 0) 
      SERR1("if no covariance model is given then '%s' must be given",
	    KNAME(MIXED_BETA));
    if (ncol[MIXED_X] != 1) // deterministic effect
      SERR("X must have one column");
    kdefault(cov, MIXED_BETA, 1);
  } else { // submodel is given
    cov_model 
      *next = cov->sub[0];   
    //    double var = 1.0;

    if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
      return ERRORDIM;
    if ((err = CHECK(next, cov->tsdim, cov->xdimown, PosDefType,
		       cov->domown, 
                       cov->isoown, SUBMODEL_DEP, ROLE_COV)) != NOERROR) {
        // print("error\n");
      return err;
    }

    if (cov->q == NULL && (err = set_mixed_constant(cov)) != NOERROR) {
      return err;
    }
    
    // warning("some checks in model 'mixed' are missing");
    // ob X mit C zusammengeht.
    
    setbackward(cov, next);
  }

  if (cov->vdim[0] > 1) 
    SERR("multivariate version of mixed not programmed yet");

  if (PisNULL(MIXED_DIST) xor PisNULL(MIXED_DIM))
    SERR2("if '%s' and '%s' must be given at the same time",
	  KNAME(MIXED_DIM), KNAME(MIXED_DIST));
  if (PisNULL(MIXED_DIST) xor PisNULL(MIXED_COORD))
    SERR2("'%s' and '%s' may not be given together", 
	  KNAME(MIXED_DIST), KNAME(MIXED_COORD));

  // incorrect. but save !!
  //  cov->semiseparatelast = false; // taken[tsxdim - 1] <= 1;
  // cov->separatelast = false;     // taken[tsxdim - 1] <= 1; ?? 
  return NOERROR;
}


void rangemixed(cov_model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  int i;

  for (i=MIXED_X; i<=MIXED_COORD; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }

  i=MIXED_DIST;
  range->min[i] = 0;
  range->max[i] = RF_INF;
  range->pmin[i] = 1e-10;
  range->pmax[i] = 1e10;
  range->openmin[i] = false;
  range->openmax[i] = true;

}


int initmixed(cov_model *cov, gen_storage  VARIABLE_IS_NOT_USED *S) {
  int store = GLOBAL.general.set;
  location_type *loc = Loc(cov);
  mixed_storage *s;
  char errorloc_save[nErrorLoc];
  double 
    *b = P(MIXED_BETA),
    *coord = P(MIXED_COORD),
    *dist = P(MIXED_DIST);
  int 
    Cn = -1, 
    Cdim = -1, 
     err = NOERROR,
    vdim = cov->vdim[0],
    dim = coord!=NULL ? cov->ncol[MIXED_COORD] : P0INT(MIXED_DIM);
  long
    totalpoints = loc->totalpoints,
    total = vdim * totalpoints;

  //  listoftype *X = PLIST(MIXED_X);
  bool distTF;

  assert(cov->vdim[0] == cov->vdim[1]);

  return ERRORFAILED; // muss in zerlegt werden in init und struct
  // und unten struct aufrufen richtig praeparieren.


  ROLE_ASSERT_GAUSS;
  
  // cholesky zerlegung + bereitstellung von b

  strcpy(errorloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "%s%s: ", errorloc_save, "init mixed model");
  
  assert(cov->nr == MIXEDEFFECT);

  NEW_STORAGE(mixed);
  s = cov->Smixed;

  NotProgrammedYet("");

  FREE(s->Xb);
  if (cov->ncol[MIXED_BETA] > 0) { // b is given
    // X is given, but no covariance model
    if (cov->nrow[MIXED_X] > 1) {
      warning("using first element of list X only");
    }
    
    if ((s->Xb = (double*) MALLOC(sizeof(double) * LNROW(MIXED_X))) == 
	NULL) { err=ERRORMEMORYALLOCATION; goto ErrorHandling; }  
   
    Ax(LP(MIXED_X), b, LNROW(MIXED_X), LNCOL(MIXED_X), s->Xb);

  } else { // submodel is given

    NotProgrammedYet(""); //if (cov->q[MIXED_CONSTANT]) ERR("not 'constant'").

    cov_model *sub, *next = cov->sub[0];
    
    if ((s->Xb = (double*) MALLOC(sizeof(double) * total)) == 0) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }  
 
    if ((distTF = dist!=NULL)) {
      // Xu , u ~ geostatmodel, u i.a. viel laenger als X (tierzucht)
      // stehende Vektoren, aneinandergereiht
      Cn = cov->ncol[MIXED_DIST];
      Cn = (int) (0.5 + 0.5 * (1.0 + sqrt(1 + 8.0 * Cn)));
      Cdim = cov->nrow[MIXED_DIST];
      assert(Cdim==1);
      if ((err = covCpy(&(cov->key), next, dist, NULL, dim, dim, Cn, 
			false /*Time */,
			false, distTF)) != NOERROR) goto ErrorHandling;
    } else if (coord != NULL) {
      // Xu , u ~ geostatmodel, u i.a. viel laenger als X (tierzucht)
      // stehende Vektoren, aneinandergereiht
      Cn = cov->ncol[MIXED_COORD];
      Cdim = cov->nrow[MIXED_COORD];
      
      if (Cdim > Cn)
	warning("The dimension of the coordinates is higher than the number of points");
      if ((err = covCpy(&(cov->key), next)) != NOERROR) goto ErrorHandling;
    } else {
      // Xu , u ~ geostatmodel, X quadratisch	
      Cn = loc->totalpoints;

      BUG; // C HECK whether cov_list = false is ok
      if ((err = covCpy(&(cov->key), next, coord, NULL, dim, dim, Cn, false, 
			false, distTF)) != NOERROR) goto ErrorHandling;
    }
    if (cov->key->nr != GAUSSPROC) addModel(&(cov->key), GAUSSPROC);   
    cov->key->calling = cov->key;
    
    NEW_STORAGE(gen);
    if ((err = INIT(cov->key, 0, cov->Sgen)) != NOERROR) goto ErrorHandling;
   
    int Xnrow, Xncol;
    Xnrow = Xncol = Cn * sub->vdim[0];
    if (loc->i_row==0 && loc->i_col==0) {
      FREE(s->mixedcov);   
      s->mixedcov = (double*) MALLOC(sizeof(double) * Xnrow * Xnrow);
      if (s->mixedcov == NULL) {
	err = ERRORMEMORYALLOCATION;
	goto ErrorHandling;
      }
    }

  } // end of submodel
 
  FieldReturn(cov);

 ErrorHandling: 
  GLOBAL.general.set = store;
  cov->initialised = err == NOERROR;
  return err;
}


static int keepRandomEffect = false;
void domixed(cov_model *cov, gen_storage  VARIABLE_IS_NOT_USED *S){
  location_type *loc = Loc(cov);
  mixed_storage *s = cov->Smixed;
  double *res  = cov->rf;
  int 
    vdim = cov->vdim[0];
  long i,
    totalpoints = loc->totalpoints,
    total = vdim * totalpoints;
   assert(cov->vdim[0] == cov->vdim[1]);

  if (cov->ncol[MIXED_BETA] > 0) { // b is given
    // X is given, but no covariance model
    if (total == LNROW(MIXED_X)) {
      for (i=0; i<total; i++) res[i] = s->Xb[i];
    } else {
      assert(LNROW(MIXED_X) == 1);
      for (i=0; i<total; i++) res[i] = s->Xb[0];
    }    
  } else { // submodel is given
    cov_model *key = cov->key;
    if (!keepRandomEffect || !s->initialized) {
      do_gaussprocess(key, cov->Sgen);
    }     
    if (PisNULL(MIXED_X)) {
      double *rf = cov->key->rf;
      for (i=0; i<total; i++) res[i] = rf[i];
    } else {
      AxResType(LP(MIXED_X), cov->key->rf, LNROW(MIXED_X), LNCOL(MIXED_X), res);
    } 
  }
}
    




//////////////////////////////////////////////////////////////////////
//    trend
//////////////////////////////////////////////////////////////////////

void trend(double *x, cov_model *cov, double *v){
  cov_model *musub = cov->kappasub[TREND_MEAN];
  int i,
    vdim = cov->vdim[0],
    vSq = vdim * vdim;
  double *mu = P(TREND_MEAN);

  if (cov->role == ROLE_COV) {
    BUG;
    for (i=0; i<vSq; i++) v[i]=0.0;
  }
  if (!isShape(cov->typus) && !isTrend(cov->typus)) BUG;

  if (musub != NULL) {
    FCTN(x, musub, v);
  } else for (i=0; i<vdim; i++) v[i]=ISNA(mu[i]) || ISNAN(mu[i]) ? 1.0 :  mu[i];
  // 1.0 notwendig fuer likelihood berechnung;
}

void kappatrend(int i, cov_model *cov, int *nr, int *nc){
  // i nummer des parameters
  int k;
  
  switch(i) {
    case TREND_MEAN: //mu
      *nr = SIZE_NOT_DETERMINED; 
      *nc = 1;
    break;
  
    case TREND_LINEAR: //plane
      *nr = cov->tsdim;
      *nc = SIZE_NOT_DETERMINED;
    break;
    
    case TREND_POLY: //polydeg
      *nr = SIZE_NOT_DETERMINED;
      *nc = 1;
     break;
    
    case TREND_PARAM_POLY: //polycoeff
      if(PisNULL(TREND_POLY)) {
        *nr = -1; 
      } else {
        *nr = 0;
	
	//APMI(cov);
        for(k=0; k < cov->nrow[TREND_POLY]; k++) {
	  //printf("%d\n", *nr);
	  *nr += binomialcoeff(PINT(TREND_POLY)[k] + 
			      cov->tsdim, cov->tsdim);
	}
      }
      *nc = 1;
    break;
    
    case TREND_FCT: //arbitraryfct
      *nr = 1;
      *nc = 1;
    break;
    
    case TREND_PARAM_FCT: //fctcoeff
      *nr = 1;
      *nc = 1;
    break;

    default:
      *nr = -1; 
      *nc = -1;
  }
}


int checktrend(cov_model *cov){
  // if (cov->ncol[TREND_LINEAR] > 0 || cov->ncol[TREND_FCT]>0
  //    || cov->ncol[TREND_PARAM_FCT] > 0)
  //  return(ERRORNOTPROGRAMMED);
  cov_model 
    *prev = cov->calling,
    *musub = cov->kappasub[TREND_MEAN];
  double *mu = P(TREND_MEAN),
    *plane = P(TREND_LINEAR);
  int i,  err,
    nkappa = CovList[cov->nr].kappas,
    *polydeg = PINT(TREND_POLY),
    vdim = 0, 
    tsdim= cov->tsdim,
    basislen = 0;
  //SEXP Rx, fctbody, envir;

  if (! ((isCoordinateSystem(cov->isoown) && musub != NULL)
	 ||
	 (isIsotropic(cov->isoown) && musub == NULL)))
    return ERRORFAILED;

 
  bool ok = (musub != NULL && isCoordinateSystem(cov->isoown) && prev != NULL && 
	     prev->nr == MULT && isTrend(prev->typus)) 
    ||   // passt das so alles??
    (musub == NULL && isIsotropic(cov->isoown) &&
     (prev == NULL || prev->nr == TRENDEVAL ||
      (prev->nr == PLUS && (prev->calling == NULL || isProcess(prev->calling) ||
			    isInterface(prev->calling) || isTrend(prev->typus)))))
     ;
  if (!ok) SERR("trend model is misused.");

 
  ROLE_ASSERT(ROLE_GAUSS || cov->role == ROLE_COV);
  for (i=0; i<nkappa; i++) {// NEVER; cf kriging.R, for instance
    if (cov->kappasub[i] != NULL && (i>0 || isRandom(cov->kappasub[i])))
      SERR("parameters may not be random or sub");
  }

  if (musub != NULL) {
    assert(cov->tsdim == cov->xdimown);
    if ((err = CHECK(musub, cov->tsdim, cov->xdimown, ShapeType,
		     XONLY, cov->isoown,
		     SUBMODEL_DEP, ROLE_BASE)) != NOERROR) return err;
    if (isRandom(musub->typus)) NotProgrammedYet("mixed effects");
    vdim = musub->vdim[0];
    cov->matrix_indep_of_x = false;
  } else if (mu != NULL) {
    if (prev != NULL && prev->nr != PLUS && isProcess(prev) && 
	!isInterface(prev) && prev->nr != TRENDEVAL) {
      cov->typus = ShapeType;
     }
    vdim = cov->nrow[TREND_MEAN];
    cov->matrix_indep_of_x = true;
  } else SERR("no trend argument given");

  
  
  if (plane != NULL) {
    SERR("'plane' currently not working. Please contact author if needed");

    if(vdim>0 && (vdim != cov->ncol[TREND_LINEAR])) {
      SERR("trend parameters have different multivariate dimensions");
    } else vdim = cov->ncol[TREND_LINEAR];
    cov->matrix_indep_of_x = false;
  }
  
  if (!PisNULL(TREND_POLY)) {
    SERR("'poly' currently not working. Please contact author if needed");

     if (vdim>0) {
      if (vdim != cov->nrow[TREND_POLY])
	SERR("trend parameters have different multivariate dimensions");
      SERR("polynomials and free functions in trend may not be mixed with other trend definitions. Please use a sum of trends.");
    }
    vdim = cov->nrow[TREND_POLY];
    if (PisNULL(TREND_PARAM_POLY)) {
      for (i=0; i<vdim; i++) 
	basislen += binomialcoeff(polydeg[i] + tsdim, tsdim);      
      PALLOC(TREND_PARAM_POLY, basislen, 1);
      for(i=0; i<basislen; i++) P(TREND_PARAM_POLY)[i] = RF_NA;
    }
    cov->matrix_indep_of_x = false;
  }
  
  if (!PisNULL(TREND_FCT)) { //brauche hier: construct.fct$vdim
    cov->matrix_indep_of_x = false;

    NotProgrammedYet("arbitrary function");

    
//     kdefault(cov, TREND_PARAM_FCT, 1.0);
//     PROTECT(envir = allocSExp(ENVSXP));
//     SET_ENCLOS(envir, R_GlobalEnv);
//     
//     PROTECT(fctbody = allocSExp(CLOSXP));
//     fctbody = BODY(*((SEXP *) arbitraryfct));
//     
//     PROTECT(Rx = allocVector(REALSXP,xlen));
//     for(j=0; j<xlen; j++) REAL(Rx)[j] = 0;
//     defineVar(install("x"), Rx, envir);
//     defineVar(install("y"), ScalarReal(0), envir);
//     defineVar(install("z"), ScalarReal(0), envir);
//     defineVar(install("T"), ScalarReal(0), envir);
//     
//     int vdimproposal = length(eval(fctbody, envir));       
//     UNPROTECT(3);
//     
//     if (vdim > 0) {
//       if (vdim != vdimproposal)
// 	GERR("trend parameters have different multivariate dimensions");
//       GERR("polynomials and free functions in trend may not be mixed with other trend definitions. Please use a sum of trends.");
//     }
//     vdim = vdimproposal;
  }
  
  if (vdim <= 0) {
    vdim = prev->vdim[0];
    if (vdim <= 0) 
      SERR("multivariate dimension for trend cannot be determined.");
    PALLOC(TREND_MEAN, vdim, 1);
    for(i=0; i<vdim; i++) P(TREND_MEAN)[i] = 0.0;
  }
  cov->vdim[0] = cov->vdim[1] = vdim;
  cov->isoown = cov->matrix_indep_of_x ? IsotropicOf(cov->isoown) 
    : CoordinateSystemOf(cov->isoown);

 return NOERROR;

}



int checktrendproc(cov_model *cov){
  return checktrend(cov);
}



void rangetrend(cov_model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  //P(TREND_MEAN]: mu / mean
  
  range->min[TREND_MEAN] = RF_NEGINF;
  range->max[TREND_MEAN] = RF_INF;
  range->pmin[TREND_MEAN] = -10^10;
  range->pmax[TREND_MEAN] = 10^10;
  range->openmin[TREND_MEAN] = true;
  range->openmax[TREND_MEAN] = true;
  
  //P(TREND_LINEAR]: plane
  range->min[TREND_LINEAR] = RF_NEGINF;
  range->max[TREND_LINEAR] = RF_INF;
  range->pmin[TREND_LINEAR] = -10^10;
  range->pmax[TREND_LINEAR] = 10^10;
  range->openmin[TREND_LINEAR] = true;
  range->openmax[TREND_LINEAR] = true;
  
  //P(TREND_POLY]: polydeg / polynomial degree
  range->min[TREND_POLY] = 0;
  range->max[TREND_POLY] = MAXINT;
  range->pmin[TREND_POLY] = 0;
  range->pmax[TREND_POLY] = 10;
  range->openmin[TREND_POLY] = false;
  range->openmax[TREND_POLY] = false;
  
  //P(TREND_PARAM_POLY]: polycoeff / coefficients of polynomial
  range->min[TREND_PARAM_POLY] = RF_NEGINF;
  range->max[TREND_PARAM_POLY] = RF_INF;
  range->pmin[TREND_PARAM_POLY] = -10^10;
  range->pmax[TREND_PARAM_POLY] = 10^10;
  range->openmin[TREND_PARAM_POLY] = true;
  range->openmax[TREND_PARAM_POLY] = true;
 
  //P(TREND_FCT]: arbitraryfct / arbitrary function
  range->min[TREND_FCT] = RF_NEGINF;
  range->max[TREND_FCT] = RF_INF;
  range->pmin[TREND_FCT] = -10^10;
  range->pmax[TREND_FCT] = 10^10;
  range->openmin[TREND_FCT] = true;
  range->openmax[TREND_FCT] = true;
  
  //P(TREND_PARAM_FCT]: fctcoeff / coefficient of arbitrary function
  range->min[TREND_PARAM_FCT] = RF_NEGINF;
  range->max[TREND_PARAM_FCT] = RF_INF;
  range->pmin[TREND_PARAM_FCT] = -10^10;
  range->pmax[TREND_PARAM_FCT] = 10^10;
  range->openmin[TREND_PARAM_FCT] = true;
  range->openmax[TREND_PARAM_FCT] = true;
}

				
void GetInternalMeanI(cov_model *cov, int vdim, double *mean){
  // assuming that mean[i]=0.0 originally for all i
  int i;
  if (cov->nr == TREND) {
    if (cov->ncol[TREND_MEAN]==1) {
      if (cov->nrow[TREND_MEAN] != vdim || cov->kappasub[TREND_MEAN] != NULL) {
	for (i=0; i<vdim; i++) mean[i] = RF_NA;
	return; // only scalar allowed !
      }
      for (i=0; i<vdim; i++) mean[i] += P(TREND_MEAN)[i];
    } 
  } else if (cov->nr == CONST && isTrend(cov)) {
    for (i=0; i<vdim; i++) mean[i] += P(CONST_C)[i];
  } else if (isTrend(cov)) {
    if (cov->xdimown < MAXSIMUDIM) {
      FCTN(ZERO, cov, mean);//to be improved! to do
    } else for (i=0; i<vdim; i++) mean[i] = RF_NA; 
  }
  if (cov->nr == PLUS || cov->nr == TREND) {
    for (i=0; i<cov->nsub; i++) GetInternalMeanI(cov->sub[i], vdim, mean);
  }
}

void GetInternalMean(cov_model *cov, int vdim, double *mean){
  int i;
  for (i=0; i<vdim; i++) mean[i]=0.0;
  GetInternalMeanI(cov, vdim, mean);
}



int binomialcoeff(int n, int k) {
  //programmed as in wikipedia
  int i, res=1;
  
  if((k < 0) || (k > n)) return 0;
  if(k > n-k) k = n-k; //symmetry
  for(i=0; i<k; i++) {
     res *= n-i;
     res /= i+1; //two steps because of integer division
  }
  return res;
}


/*

int init_trend(cov_model *cov, gen_storage *S) {
  
  long err = NOERROR;
  int i, 
    *polydeg = PINT(TREND_POLY),
    basislen=0,
    tsdim = cov->tsdim,
    vdim = cov->vdim[0];
  //SEXP fctformals, argnames;
  if (cov->vdim[0] != cov->vdim[1]) BUG;

  //assert(false);

  ROLE_ASSERT_GAUSS;
   
  NEW_STORAGE(trend);
  trend_storage *s = cov->Strend;

  if(polydeg != NULL) {
    for(i=0; i<vdim; i++) 
      basislen += binomialcoeff(polydeg[i]+tsdim,tsdim);
  }
  if ((s->x = (double *) MALLOC(sizeof(double) * tsdim))==NULL ||
      (s->xi = (int *) MALLOC(sizeof(int) * tsdim))==NULL ||
      (s->evalplane = (double *) MALLOC(sizeof(double) * vdim))==NULL || 
      (basislen > 0 && 
       (s->powmatrix = (int *) MALLOC(sizeof(int) * basislen * tsdim))==NULL)){
    err=ERRORMEMORYALLOCATION;   
    goto ErrorHandling;
  }
 
  if (basislen > 0) poly_basis(cov, S); //generates basis of monomials
  //each row consists of one basis element
  //the j-th column consists of the power of the j-th space-time-dimension

 
  
  if (!PisNULL(TREND_FCT)) { //hier werden Argumente von arbitraryfct ueberprueft
    NotProgrammedYet("");
//       fctformals = getAttrib(FORMALS(*((SEXP *) arbitraryfct)), R_NamesSymbol);
//       nargs = length(fctformals);
//       PROTECT(argnames = allocVector(STRSXP,1));
//       
//       SET_STRING_ELT(argnames,0,mkChar("y"));
//       matchind = 0;
//       for(j=0; j<nargs; j++) 
// 	 matchind += (STRING_ELT(argnames,0) == STRING_ELT(fctformals,j));
//       if (matchind>0) {
//         if((meth->loc->spatialdim < 2) || (meth->loc->xvectorvalued))
// 	  GERR("The variable y does not match to the locations.\n");
//       }
//       
//       SET_STRING_ELT(argnames,0,mkChar("z"));
//       matchind = 0;
//       for(j=0; j<nargs; j++) 
// 	 matchind += (STRING_ELT(argnames,0) == STRING_ELT(fctformals,j));
//       if (matchind>0) {
// 	if((meth->loc->spatialdim < 3) || (meth->loc->xvectorvalued))
// 	  GERR("The variable z does not match to the locations.\n")	
//       }
//       
//       SET_STRING_ELT(argnames,0,mkChar("T"));
//       matchind = 0;
//       for(j=0; j<nargs; j++) 
// 	 matchind += (STRING_ELT(argnames,0) == STRING_ELT(fctformals,j));
//       if (matchind>0) {
//         if (meth->loc->Time == false) 
// 	  GERR("The variable T may be used for time only.\n")		  
//       }
//       UNPROTECT(1);
  }

  err = FieldReturn(cov);

  return NOERROR;
  
  ErrorHandling:
   return err;
}


void do_trend(cov_model *cov, gen_storage  VARIABLE_IS_NOT_USED *s){
  location_type *loc = Loc(cov);
  char errorloc_save[nErrorLoc];
  trend_storage *S = cov->Strend;
  //PMI(cov->calling->calling);
  assert(S != NULL);
  double t,
    *mu = P(TREND_MEAN),
    //  *plane   = P(TREND_LINEAR),
    //  *polycoeff = P(TREND_PARAM_POLY),
    //    *fctcoeff = P(TREND_PARAM_FCT),
    **xgr = loc->xgr,
    // *x = S->x,
    // *evalplane = S->evalplane
    ;

  int  v, w, 
    basislen, startindex,
    *polydeg = PINT(TREND_POLY),
    vdim = cov->vdim[0],
    tsdim = cov->tsdim,
    spatialdim = loc->spatialdim,
    *xi = S->xi,
    *powmatrix = S->powmatrix;
  //SEXP fctbody, tempres, envir, Rx;

  long i, j, k,
    totalpoints = loc->totalpoints,
    total = totalpoints * vdim;
 
  double *res = cov->rf;
  
  assert(cov->vdim[0] == cov->vdim[1]);

  strcpy(errorloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "%s%s: ", errorloc_save, "add trend model");
 
  // print("%s\n", ERROR_LOC);

  if (cov->kappasub[TREND_MEAN] != NULL) {
    
  } else if (mu != NULL) {  
    for (k=0; k<total; ) {
      for (v=0; v<vdim; res[k++] = mu[v++]); 
    } 
  } else {
    for (k=0; k<total; res[k++]=0.0);
  }

  if (plane != NULL) {  
    if (loc->grid) {
      for (w=0; w<tsdim; w++)  {
	x[w]=xgr[w][XSTART];
	xi[w]=0;
      }
      long endfor = totalpoints * vdim;
      for (k=0; k<endfor; ) {
	xA(x, plane, cov->nrow[TREND_LINEAR], cov->ncol[TREND_LINEAR],
	   evalplane);
	for(v=0; v<vdim; v++) res[k++] += evalplane[v];

	i = 0;
	(xi[i])++;
	x[i] += xgr[i][XSTEP];
	while(xi[i]>=loc->xgr[i][XLENGTH]) {
	  xi[i] = 0;
	  x[i] = xgr[i][XSTART];
	  if (i<tsdim-1) {
	    i++;
	    (xi[i])++;
	    x[i] += xgr[i][XSTEP];
	  } else {
	    assert(k==endfor);
	  }
	}
      }	
    } else if (loc->Time) {
      long m, 
	endfor= (long int) loc->xgr[loc->timespacedim-1][XLENGTH],
	endfor2 = loc->spatialtotalpoints;
      for(k=m=0, t=loc->T[XSTART]; m<endfor; m++, t+=loc->T[XSTEP]) {
	// for(t=loc->T[XSTART], i=0; t < loc->T[XEND]; t += loc->T[XSTEP]) {
	for(m=0, j=0; m < endfor2; m++) {
	  for(w=0; w<spatialdim; w++) x[w] = (loc->x)[j++];
	  x[spatialdim] = t;
	  xA(x, plane, cov->nrow[TREND_LINEAR], cov->ncol[TREND_LINEAR],
	     evalplane);
	  for(v=0;v<vdim;v++) res[k++] += evalplane[v];
	}
      }
    } else {
      int m,
	endfor = totalpoints * tsdim;
      for (k=m=0; m<endfor; m+=tsdim) {
	xA(loc->x + m, plane, cov->nrow[TREND_LINEAR], 
	   cov->ncol[TREND_LINEAR], evalplane);
	for(v=0;v<vdim;v++) res[k++] += evalplane[v];
      }
    }
  }
  
  if(polydeg != NULL) {
    startindex=0;
    long end_k = totalpoints * vdim;
    for(v=0; v<vdim; v++) { 
      basislen = binomialcoeff(polydeg[v] + tsdim, tsdim);
      if(isnan(polycoeff[0])) {
	ERR("Error: cannot evaluate polynomial without coefficients.\n");
      }
      if (loc->grid) {
	for (w=0; w<tsdim; w++)  {
	  x[w]=xgr[w][XSTART];
	  xi[w]=0;
	}
	for(k=v; k<end_k; k+=vdim) {
	  //evaluation of trend polynomial
	  res[k] += evalpoly(x, powmatrix + startindex*tsdim,
			     polycoeff + startindex, basislen, tsdim);
	  i = 0;
	  (xi[i])++;
	  x[i] += xgr[i][XSTEP];
	  while(xi[i]>=xgr[i][XLENGTH]) {
	    xi[i] = 0;
	    x[i] = xgr[i][XSTART];
	    if (i<tsdim-1) {
	      i++;
	      (xi[i])++;
	      x[i] += xgr[i][XSTEP];
	    } else {
	      assert(k==v + end_k - vdim);
	    }
	  }
	}	
      } else if(loc->Time) {
	 for(t=loc->T[XSTART], k=v; k<end_k; t+=loc->T[XSTEP]) {
	   long endfor = loc->spatialtotalpoints;
	   for(i=0, j=0; i<endfor; i++, k+=vdim) {
             for(w=0; w<spatialdim; w++)  x[w] = (loc->x)[j++];
	     x[spatialdim] = t;
	     //evaluation of trend polynomial
             res[k] += evalpoly(x, powmatrix + startindex*tsdim,
				    polycoeff + startindex, basislen, tsdim);
           }
	 }
      } else {
	for(k=v; k<end_k; k+=vdim) {
	  //evaluation of trend polynomial
	  res[k] += evalpoly(x, powmatrix + startindex*tsdim,
			     polycoeff + startindex, basislen, tsdim);
	}
      }
      startindex += basislen;
    }
  }
  
  if (!PisNULL(TREND_FCT)) { //muss hier arbitraryfct auswerten
    NotProgrammedYet("");
//      if (isnan(fctcoeff[0])) {
// 	ERR("Error: cannot evaluate function without coefficient.\n");
//      }
//      
//      PROTECT(envir = allocSExp(ENVSXP));
//      SET_ENCLOS(envir, R_GlobalEnv);
//      
//      PROTECT(fctbody = allocSExp(CLOSXP));
//      fctbody = BODY(*((SEXP *) arbitraryfct));
//      
//      PROTECT(tempres = allocVector(REALSXP, vdim));
//      PROTECT(Rx = allocVector(REALSXP,xlen));     
//      
//      if (loc->grid) {
//          for (w=0; w<tsdim; w++)  {
// 	    x[w]=xgr[w][XSTART];
// 	    xi[w]=0;
//          }
//          for(k=0; k<totalpoints; k++) {
// 	   //evaluation of trend polynomial
//            for(w=0; w<xlen; w++) REAL(Rx)[w] = x[w];
//            defineVar(install("x"), Rx, envir);
//            if (spatialdim>1)
//               defineVar(install("y"), ScalarReal(x[1]), envir);
//            if (spatialdim>2)
//               defineVar(install("z"), ScalarReal(x[2]), envir);
//            defineVar(install("T"), ScalarReal(x[tsdim-1]), envir);
// 	   tempres = eval(fctbody, envir);
// 	   for (v=0; v<vdim; v++)
// 	     res[k*vdim+v] += fctcoeff[0]*REAL(tempres)[v];
//            i = 0;
//            (xi[i])++;
//            x[i] += xgr[i][XSTEP];
//            while(xi[i]>=len[i]) {
//              xi[i] = 0;
//              x[i] = xgr[i][XSTART];
//              if (i<tsdim-1) {
// 	       i++;
// 	       (xi[i])++;
// 	       x[i] += xgr[i][XSTEP];
//              } else {
// 	       assert(k==totalpoints-1);
//              }
//            }
//          }	
//      } else if(loc->Time) {
// 	 int k, 
// 	   endfor= (int) loc->T[XLENGTH];
// 	 for(k=0, t=loc->T[XSTART], i=0; k<endfor; k++, t += loc->T[XSTEP]) {
// 	   //	 for(t=loc->T[XSTART], i=0; t < loc->T[XEND]; t+=loc->T[XSTEP]) {
// 	   for(k=0, j=0; k < loc->spatialtotalpoints; k++, i++) {
//              for(w=0; w<spatialdim; w++) x[w] = (loc->x)[j++];
// 	     x[spatialdim] = t;
// 	     //evaluation of trend polynomial
//              for(w=0; w<xlen; w++) REAL(Rx)[w] = x[w];
//              defineVar(install("x"), Rx, envir);
//              if (spatialdim > 1)
//                defineVar(install("y"), ScalarReal(x[1]), envir);
//              if (spatialdim > 2)
//                defineVar(install("z"), ScalarReal(x[2]), envir);
//              if (loc->Time)
//                defineVar(install("T"), ScalarReal(x[tsdim-1]), envir);
//              tempres = eval(fctbody, envir);
// 	     for (v=0; v<vdim; v++)
// 	       res[i*vdim+v] += fctcoeff[0]*REAL(tempres)[v];
//            }
// 	 }
//        } else {
// 	 for(k=0; k<totalpoints; k++) {
// 	   //evaluation of trend polynomial
// 	   for(w=0; w<spatialdim; w++) x[w] = (loc->x)[k*spatialdim+w];
//            for(w=0; w<xlen; w++) REAL(Rx)[w] = x[w];
//            defineVar(install("x"), Rx, envir);
//            if (spatialdim>1)
//               defineVar(install("y"), ScalarReal(x[1]), envir);
//            if (spatialdim>2)
//               defineVar(install("z"), ScalarReal(x[2]), envir);
//            defineVar(install("T"), ScalarReal(x[tsdim-1]), envir);
// 	   tempres = eval(fctbody, envir);
// 	   for (v=0; v<vdim; v++)
// 	     res[k*vdim+v] += fctcoeff[0]*REAL(tempres)[v];
// 	   }
//        }
//        UNPROTECT(4);
  }  
  
  
  return;

}

void poly_basis(cov_model *cov, gen_storage  VARIABLE_IS_NOT_USED *s) {
  
  trend_storage *S = cov->Strend;
  int basislen=0, powsum, d, i, j, k, v,
    dim = cov->tsdim,
    vdim = cov->vdim[0],
    *powmatrix = S->powmatrix,
      *dimi=NULL,
      err=NOERROR;
  int *polydeg = PINT(TREND_POLY);
  assert(cov->vdim[0] == cov->vdim[1]);
  
  dimi = (int *) MALLOC(dim * sizeof(int));
  if (dimi == NULL) {
     err = ERRORMEMORYALLOCATION;
     goto ErrorHandling;
  }
  
  for(v=0, j=0; v<vdim; v++) {
    basislen = binomialcoeff(polydeg[v] + dim, dim);
    //   print("v: %d\n", v);
    for(d=0; d<dim; d++) dimi[d] = 0;
    for(k=0; k<basislen; k++) {
      for(d=0; d<dim; d++) powmatrix[j++] = dimi[d];
      // for(d=0;d<dim;d++) print("%d ", powmatrix[j*dim+d]);
      // print("\n");
      i=0;
      (dimi[i])++;
      powsum = 0;
      for(d=0; d<dim; d++) powsum += dimi[d];
      while(powsum > polydeg[v]) {
        dimi[i] = 0;
        if (i < dim-1) {
	  i++;
	  (dimi[i])++;
        }
        powsum = 0;
        for(d=0; d<dim; d++) powsum += dimi[d];
      }
    }
    // print("\n");
  }
  //print("\n");
  
  ErrorHandling:
  FREE(dimi);   
  if (err != NOERROR) XERR(err);
   
  return;
  
}

double evalpoly(double *x, int *powmatrix, double *polycoeff, int basislen,
	      int dim) {
  int i, d, j;
  double res = 0, tempres;
  for(j=i=0; i<basislen; i++) {
     tempres=1;
     for(d=0; d<dim; d++) tempres *= pow(x[d], powmatrix[j++]);
     res += polycoeff[i] * tempres;
  }
  return(res);
}

void likelihood_trend(double VARIABLE_IS_NOT_USED *x, double VARIABLE_IS_NOT_USED invmatrix, cov_model *cov, 
		      double VARIABLE_IS_NOT_USED *v){ 
  if (cov->role == ROLE_GAUSS) {
    
    NotProgrammedYet("");
      //C = covarianzmatrix
      // b = (X^T C^{-1} X)^{-1} X^top C^{-1} y
      // Baysiean mit b ~ N(b_0, C_B) :
      //  b = (X^T C^{-1} X + C_b^{-1})^{-1} [X^top C^{-1} y + C_b^{-1} b_0
  } else {
    // deterministic is OK
    NotProgrammedYet("");
  }

}
*/




int checkTrendEval(cov_model *cov) { // auch fuer TrendEval
   cov_model *next = cov->sub[0];
   int err,
     dim =  Gettimespacedim(cov);

  if ((err = CHECK(next, cov->tsdim, cov->xdimown, TrendType,
		     XONLY, cov->isoown,
		     SUBMODEL_DEP, ROLE_BASE)) != NOERROR) return err;

  assert(isTrend(cov->sub[0]->typus));
  setbackward(cov, next);
  cov->vdim[0] = next->vdim[0]; 
  cov->vdim[1] = next->vdim[1];
  if (cov->vdim[0] != 1) NotProgrammedYet("");
  KAPPA_BOXCOX;

  // to do alloc_cov nicht optimal hinsichtlich speicher etc. ,-  vdim[1] waere 
  // eigentlich das richtige
  if ((err = alloc_cov(cov, dim, cov->vdim[0], cov->vdim[0])) != NOERROR) 
    return err;
  
  // EXTRA_STORAGE;

  return NOERROR;
}
 

int init_TrendEval(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s){// auch fuer TrendEval
  int err;

  if (cov->vdim[0] != 1) NotProgrammedYet("");

  if ((err = check_fctn(cov)) != NOERROR) return err;
  
  ROLE_ASSERT_GAUSS;
  err = FieldReturn(cov);
  cov->simu.active = err == NOERROR;
  if (PL>= PL_STRUCTURE) PRINTF("\n'%s' is now initialized.\n", NAME(cov));
  return err;
}

void do_TrendEval(cov_model *cov, gen_storage  VARIABLE_IS_NOT_USED *s){
  double
    *res = cov->rf;
  assert(res != NULL);
  char errorloc_save[nErrorLoc];
  strcpy(errorloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "%s%s: ", errorloc_save, "add trend model");
  Fctn(NULL, cov, res);
  strcpy(ERROR_LOC, errorloc_save);

  BOXCOX_INVERSE;
  return; 
}


void range_TrendEval(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){  
  GAUSS_COMMON_RANGE;
}
