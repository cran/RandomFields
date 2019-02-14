
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives of hypermodels

 Copyright (C) 2015 -- 2017 Martin Schlather

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

#include <Rmath.h>
#include "questions.h"
#include "primitive.h"
#include "primitive.others.h"
#include "operator.h"
#include "Processes.h"
#include "variogramAndCo.h"
#include "rf_interfaces.h"

// *sub = cov->key == NULL ? cov->sub[0] : cov->key;
#define SCALAR_P(A,B,C) Ext_scalarX(A,B,C, SCALAR_AVX)
#define START_GAUSS_INTERFACE					\
  int cR = INTEGER(model_reg)[0];				\
  if (cR < 0 || cR > MODEL_MAX) BUG;				\
  set_currentRegister(cR);			\
  model *cov = KEY()[cR];			\
  assert (cov != NULL && equalsnowInterface(cov)) 			\
  cov = cov->key != NULL ? cov->key : cov->sub[0];		\
  assert(isnowProcess(cov)) 					\
  if (COVNR != GAUSSPROC)					\
    ERR("register not initialised as Gaussian likelihood");	\
  model *calling = cov->calling;					\
  if (calling == NULL ||						\
      (CALLINGNR != LIKELIHOOD_CALL && CALLINGNR != LINEARPART_CALL)	\
      ) BUG;								\
  likelihood_storage *L = cov->Slikelihood;			\
  if (L == NULL) ERR("register not initialised as likelihood method");\
  int store = GLOBAL.general.set;			       


int matrixcopyNA(double *dest, double *src, double *cond,
		  int rows, int cols, int append_cond) {
  int k = 0;
  double *ptsrc = src;
  for (int j=0; j<cols; j++) {
    for (int i=0; i<rows; i++, ptsrc++) {
      if (!ISNAN(cond[i])) dest[k++] = *ptsrc;
    }
  }
  
  double *C = cond;
  int m = 0;
  for (int j=0; j<append_cond; j++) {
    for (int i=0; i<rows; i++, m++) {
      if (!ISNAN(C[m])) dest[k++] = C[m];
    }
  }

  if (k == 0) RFERROR("one of the data set seems to consist of NAs only");
  int notnas = k == 0 ? 0 : k / (append_cond + cols); 
  
  assert(notnas * (append_cond + cols) == k);
  return notnas;
}


void SqMatrixcopyNA(double *dest, double *src, double *cond, int rows) {
  int k = 0;
  for (int j=0; j<rows; j++) {
    if (!ISNAN(cond[j])) {
      double *ptsrc = src + j * rows;     
      for (int i=0; i<rows; i++, ptsrc++) {
	if (!ISNAN(cond[i]))  dest[k++] = *ptsrc;
      }
    }
  }
}


SEXP set_boxcox(SEXP boxcox) {
  double *bc = REAL(boxcox);
  int len = length(boxcox);
  for (int i=0; i<len; i++) GLOBAL.gauss.boxcox[i] = bc[i];
  GLOBAL.gauss.loggauss = false;
  return R_NilValue;
}

SEXP get_boxcox() {
  SEXP ans;
  int len = 2 * MAXGAUSSVDIM;
  PROTECT(ans =allocVector(REALSXP, len));
  if (GLOBAL.gauss.loggauss) {
    for (int i=0; i<len; REAL(ans)[i++] = 0.0);
  } else  for (int i=0; i<len; i++) REAL(ans)[i] = GLOBAL.gauss.boxcox[i];
  UNPROTECT(1);
  return ans;
}

 
void boxcox_inverse(double boxcox[], int vdim, double *res, int pts, 
		    int repet) {	
  // (y * L + 1)^{1/L} - mu

  for (int r=0; r<repet; r++) {
    for (int v=0; v<vdim; v++) {
      double lambda = boxcox[2 * v],
	invlambda = 1.0 / lambda,
	mu = boxcox[2 * v + 1];
      if ((!ISNA(lambda) && FABS(lambda) < 1e-20)) {	
	for (long i=0; i<pts; i++) {
	  res[i] = EXP(res[i]) - mu; 
	}
      } else if (ISNA(lambda) || lambda != RF_INF)
	for (int i=0; i<pts; i++) {
	  double dummy = lambda * res[i] + 1.0;
	  if ((dummy < 0 && lambda != CEIL(lambda)) || 
	      (dummy == 0 && invlambda <= 0) ) {
	    RFERROR("value(s) in the inverse Box-Cox transformation not positive");
	  }
	  res[i] = POW(dummy, invlambda) - mu;
	}
    }
  }
}

void boxcox_trafo(double boxcox[], int vdim, double *res, long pts, int repet){	
  // box cox: (x + mu)^L - 1) / L
  for (int r=0; r<repet; r++) {
    for (int v=0; v<vdim; v++) {
      double lambda = boxcox[2 * v],
	invlambda = 1.0 / lambda,
	mu = boxcox[2 * v + 1];
      if ((!ISNA(lambda) && FABS(lambda) < 1e-20)) {	
	// to do taylor
	for (long i=0; i<pts; i++) {
	  double dummy = res[i] + mu;
	  if (dummy < 0 || (dummy == 0 && lambda <= 0) ) {
	    RFERROR("value(s) in the Box-Cox transformation not positive");
	  }
	  res[i] = LOG(dummy);			
	}
      } else if (ISNA(lambda) || lambda != RF_INF) {
	for (long i=0; i<pts; i++) {
	  double dummy = res[i] + mu;
	  if ((dummy < 0 && lambda != CEIL(lambda)) 
	      || (dummy == 0 && lambda <= 0) ) {
	    RFERROR("value(s) in the Box-Cox transformation not positive");
	  }
	  res[i] = (POW(dummy, lambda) - 1.0) * invlambda; 
	}
      }
    }
  }
}


SEXP BoxCox_trafo(SEXP boxcox, SEXP res, SEXP Vdim, SEXP inverse){	
  int vdim = INTEGER(Vdim)[0],
    repet = isVector(res) ? 1 : ncols(res),
    pts = isVector(res) ? length(res) / vdim : nrows(res);
  if (vdim > MAXGAUSSVDIM)
    RFERROR2("multi-dimensionality, %d, exceeds maximum, %d, in Box-Cox transformation", 
	     vdim, MAXGAUSSVDIM);
  if (length(res) != pts * vdim * repet)
    RFERROR("multi-dimensionality incorrect in Box-Cox transformation");
  if (length(boxcox) < 2 * vdim) RFERROR("too few entries in 'boxcox'");
  if (LOGICAL(inverse)[0])
    boxcox_inverse(REAL(boxcox), vdim, REAL(res), pts, repet);
  else 
    boxcox_trafo(REAL(boxcox), vdim, REAL(res), pts, repet);
  return R_NilValue;
}



SEXP gauss_linearpart(SEXP model_reg, SEXP Set){
  START_GAUSS_INTERFACE;
#define ll 3
  const char *names[ll] = {"Y", "X", "vdim"};
  SEXP ans, namevec, X, Y;
  int 
    unp = 4, // ans, namevec, X, Y
    sets = L->sets, 
    element  = INTEGER(Set)[0],
    vdim = VDIM0,
    betatot = L->betas[L->fixedtrends];

  if (element > 0  && element > sets) ERR("set number out of range");

  PROTECT(ans = allocVector(VECSXP, ll));
  PROTECT(namevec = allocVector(STRSXP, ll));
  for (int k=0; k<ll; k++) SET_STRING_ELT(namevec, k, mkChar(names[k]));
  
  if (element <= 0) {
    PROTECT(Y = allocVector(VECSXP, sets)); 
    PROTECT(X = allocVector(VECSXP, sets)); 
    for (GLOBAL.general.set = 0; GLOBAL.general.set < sets;
	 GLOBAL.general.set++) {
      int
	totptsvdim = Gettotalpoints(cov) * vdim;
      SEXP partY , partX;
      PROTECT(partY = allocVector(REALSXP, totptsvdim));
      MEMCOPY(REAL(partY), L->YhatWithoutNA[GLOBAL.general.set], 
	      totptsvdim * sizeof(double));
      SET_VECTOR_ELT(Y, GLOBAL.general.set, partY);
      UNPROTECT(1);
      
      if (L->fixedtrends) {
	PROTECT(partX = allocMatrix(REALSXP, totptsvdim, betatot));
	MEMCOPY(REAL(partX), L->X[GLOBAL.general.set], 
		totptsvdim * betatot * sizeof(double));
	SET_VECTOR_ELT(X, GLOBAL.general.set, partX);
	UNPROTECT(1);
      } else {
	SET_VECTOR_ELT(ans, GLOBAL.general.set, R_NilValue);    
      }
    }
  } else {
    GLOBAL.general.set = element - 1;
    int totptsvdim = Gettotalpoints(cov) * vdim;
   
    PROTECT(Y = allocVector(REALSXP, totptsvdim));
    MEMCOPY(REAL(Y), 
	    L->YhatWithoutNA[GLOBAL.general.set], totptsvdim * sizeof(double));
  
    if (L->fixedtrends) {
      PROTECT(X = allocMatrix(REALSXP, totptsvdim, betatot));
      MEMCOPY(REAL(X), L->X[GLOBAL.general.set],
	       totptsvdim * betatot * sizeof(double));
    } else {
      X = R_NilValue; 
      unp--;
    }
  }
  
  int k=0;
  SET_VECTOR_ELT(ans, k++, Y);    
  SET_VECTOR_ELT(ans, k++, X);    
  SET_VECTOR_ELT(ans, k++, ScalarInteger(vdim));
  assert(k == ll);

  setAttrib(ans, R_NamesSymbol, namevec);
  UNPROTECT(unp);

  GLOBAL.general.set = store;
  return ans;
}

 

#define LOGLI_RESIDUALS 0
#define LOGLI_WHOLETREND 1
void get_logli_residuals(model *cov, double *work, double *ans,
			 int modus) {
  assert(isnowProcess(cov));
  likelihood_storage *L = cov->Slikelihood;
  assert(L != NULL);
  listoftype *datasets = L->datasets;
  assert(datasets != NULL);
  assert(ans != NULL);

  int
    vdim = VDIM0,
    betas = L->betas[L->fixedtrends],
    ncol = NCOL_OUT_OF(datasets),
    repet = ncol / vdim,
    nrow = NROW_OUT_OF(datasets),
    ndata = nrow * ncol,
    totptsvdim = nrow * vdim;
  double *X = L->X[GLOBAL.general.set],
    *pres = ans,
    *data = SET_OUT_OF(datasets);

  if (modus == LOGLI_RESIDUALS) {
    MEMCOPY(pres, data, ndata * sizeof(double));
    if (R_FINITE(P(GAUSS_BOXCOX)[0]) && R_FINITE(P(GAUSS_BOXCOX)[1]))
      boxcox_trafo(P(GAUSS_BOXCOX), vdim, pres, nrow, repet);
  } else {
    // BUG; printf("what's that good for?");
    for (int i=0; i<ndata; pres[i++] = 0.0);
  }

  //   for (int i=0; i<ndata; i++) printf("%10g ", pres[i]);  BUG;
  //  printf("\n%d %d %d vdim=%d\n", L->ignore_trend, L->dettrends, L->fixedtrends, vdim);

  if (L->ignore_trend) {
     return;
  }

  bool delete_work;
  if ((delete_work = work == NULL)) {
    work = (double *) MALLOC(nrow * vdim * sizeof(double));
  }

  double *betavec = L->betavec;
  int r, z;
  
  if (L->dettrends) {
    for (int i=0; i<L->dettrends; i++) {      
      if (L->nas_det[i]) {// non-linear parts without linear part
	FctnIntern(cov, L->cov_det[i], L->cov_det[i], work, true);
	for (z=r=0; r<repet; r++)
	  for (int k=0; k<totptsvdim; k++) pres[z++] -= work[k];
      }
    }
    for (z=r=0; r<repet; r++)
      for (int k=0; k<totptsvdim; k++)
	pres[z++] -= L->YhatWithoutNA[GLOBAL.general.set][k];
  }
  
  if (L->fixedtrends) {
    for (r=0; r<repet; r++, betavec += betas) {
      if (r==0 || L->betas_separate) {
	for (int i=0; i<totptsvdim; work[i++] = 0.0);
	//  for (int j=0; j<totptsvdim;j ++) printf("j=%d p %10g d %10g w %10g\n ", j, pres[j], data[j], work[j]);
	for (int i=0; i<betas; i++) {
	  double beta = betavec[i];
	  //	    printf("i=%d beta = %10g %d\n", i, beta, nrow);
	  for (int j=0; j<nrow; j++, X++) {	
	    //	      printf("i=%d %d %10g %10g\n", i, betas, beta, *X);
	    work[j] += *X * beta;
	  }
	}
      }
      
      for (int j=0; j<nrow; j++, data++, pres++) {
	  *pres -= work[j];
	  //    printf("j=%d p %10g d %10g w %10g\n ", j, *pres, *data, work[j]);
      }
    }
  }
  
  // for (int j=0; j<nrow; j++, data++, pres++) printf("%10g %10g\n", *pres, *data);
  if (modus != LOGLI_RESIDUALS) {
    //BUG; printf("whats that good for?");
    for (int i=0; i<ndata; i++) pres[i] = -pres[i];
  }
  
 
  if (delete_work) FREE(work);
}


SEXP get_logli_residuals(model *cov, int modus) {
  likelihood_storage *L = cov->Slikelihood;		       
  listoftype *datasets = L->datasets;
  int 
    vdim = VDIM0,
    sets = GET_LOC_SETS(cov);
  bool matrix = false;
  SEXP all_res, res;

  assert(!L->betas_separate || sets == 1);

  int max = 0;
  for (GLOBAL.general.set = 0; GLOBAL.general.set < sets; 
       GLOBAL.general.set++){
    int
      nrow = NROW_OUT_OF(datasets),
      totptsvdim = nrow * vdim;
    if (max < totptsvdim) max = totptsvdim;
  }
  if (L->work == NULL)  L->work = (double*) MALLOC(max * sizeof(double));
    
  PROTECT(all_res =  allocVector(VECSXP, sets));
  for (GLOBAL.general.set = 0; GLOBAL.general.set < sets; GLOBAL.general.set++)
    if ((matrix = NCOL_OUT_OF(datasets) > 1)) break;

  for (GLOBAL.general.set = 0; GLOBAL.general.set < sets; GLOBAL.general.set++){
    int 
      ncol = NCOL_OUT_OF(datasets),
      nrow = NROW_OUT_OF(datasets);

    if (matrix) {
      PROTECT(res = allocMatrix(REALSXP, nrow, ncol));
    } else {
      PROTECT(res =allocVector(REALSXP, nrow));
    }

    get_logli_residuals(cov, L->work, REAL(res), modus);
    SET_VECTOR_ELT(all_res, GLOBAL.general.set, res);
    UNPROTECT(1);
  }
  
  UNPROTECT(1);
  return all_res;
}

SEXP get_logli_residuals(SEXP model_reg) {
  START_GAUSS_INTERFACE;
  SEXP ans = get_logli_residuals(cov, LOGLI_RESIDUALS);
  GLOBAL.general.set = store;
  return ans;
}

SEXP get_logli_wholetrend(SEXP model_reg) {
  //BUG;  printf("whats that good for?");
  START_GAUSS_INTERFACE;
  SEXP ans = get_logli_residuals(cov, LOGLI_WHOLETREND);
  GLOBAL.general.set = store;
  return ans;
}




void get_fx(model *predict, model *cov, double *v, int set) {
  // muss angenommen werden, dass alles gesetzt ist.
  likelihood_storage *L = cov->Slikelihood;			
  int store = GLOBAL.general.set;
  assert (set >= 0 && set < GET_LOC_SETS(predict));
  GLOBAL.general.set = set;
  listoftype *datasets = L->datasets;
  int 
    err = NOERROR,
    betatot = L->betas[L->fixedtrends],    
    vdim = VDIM0,
    ncol = NCOL_OUT_OF(datasets),
    repet = L->betas_separate ? ncol / vdim : 1,
    pred_tot = Gettotalpoints(predict),
    predtotvdim = pred_tot * vdim,
    predtotvdimrepet = pred_tot * ncol;
  double *X = NULL;
  for (int k=0; k<predtotvdimrepet; v[k++] = 0.0);
  if (L->ignore_trend) goto ErrorHandling; // exist
  
  
  if (!L->betas_separate && repet > 1) GERR("BUG");
  
  if ((X = (double*) MALLOC(predtotvdim * sizeof(double))) == NULL) {
    // printf("%d %d \n", pred_tot, vdim); APMI(cov);
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  
  int r,m;
  for (int i=0; i<L->dettrends; i++) {      
    FctnIntern(predict, L->cov_det[i],  L->cov_det[i], X, true);
    for (r = m = 0; r < repet; r++) {
      for (int k=0; k<predtotvdim; k++) v[m++] += X[k];
    }
  }
  /*if(L->dettrends){
    for (z=r=0; r<repet; r++)
    for (int k=0; k<predtotvdim; k++)
        v[z++] += L->YhatWithoutNA[GLOBAL.general.set][k];
  }*/
  for (int i=0; i<L->fixedtrends; i++) {      
    FctnIntern(predict, L->cov_fixed[i],  L->cov_fixed[i], X, true);
    if (L->betas[i+1] - L->betas[i] != 1) BUG;
    double *betavec = L->betavec + L->betas[i];
    for (r = m = 0; r < repet; r++) {
      double beta = *betavec;
      	//printf("beta = %10g %10g\n", beta, X[0]); BUG;
      for (int k=0; k<predtotvdim; k++) v[m++] += X[k] * beta;
      if (L->betas_separate) betavec += betatot;
    }
  }
  
  ErrorHandling :
    GLOBAL.general.set = store;
  FREE(X);
  OnErrorStop(err, predict);
}


void get_F(model *cov, double *work, double *ans) {
  assert(isnowProcess(cov));
  likelihood_storage *L = cov->Slikelihood;
  assert(L != NULL);
  listoftype *datasets = L->datasets;
  assert(datasets != NULL);
  assert(ans != NULL);
  
  int
    vdim = VDIM0,
    betas = L->betas[L->fixedtrends],
    ncol = NCOL_OUT_OF(datasets),
    repet = ncol / vdim,
    nrow = NROW_OUT_OF(datasets),
    //    ndata = nrow * ncol,
    totptsvdim = nrow * vdim;
  double *X = L->X[GLOBAL.general.set],
    *pres = ans,
    *data = SET_OUT_OF(datasets);
  
  //MEMCOPY(pres, data, ndata * sizeof(double));
  if (R_FINITE(P(GAUSS_BOXCOX)[0]) && R_FINITE(P(GAUSS_BOXCOX)[1]))
    boxcox_trafo(P(GAUSS_BOXCOX), vdim, pres, nrow, repet);
  
  //   for (int i=0; i<ndata; i++) printf("%10g ", pres[i]);  BUG;
  //  printf("\n%d %d %d vdim=%d\n", L->ignore_trend, L->dettrends, L->fixedtrends, vdim);
                  
  if (L->ignore_trend) {
    return;
  }
                  
  bool delete_work;
  if ((delete_work = work == NULL)) {
    work = (double *) MALLOC(nrow * vdim * sizeof(double));
  }
  
  double *betavec = L->betavec;
  int r, z;
  
  if (L->dettrends) {
    
    for (int i=0; i<L->dettrends; i++) {      
      if (L->nas_det[i]) {
	FctnIntern(cov, L->cov_det[i], L->cov_det[i], work, true);
	for (z=r=0; r<repet; r++)
	  for (int k=0; k<totptsvdim; k++) pres[z++] += work[k];
      }
    }
    for (z=r=0; r<repet; r++)
      for (int k=0; k<totptsvdim; k++)
	pres[z++] += L->YhatWithoutNA[GLOBAL.general.set][k];
  }
  
  if (L->fixedtrends) {
    
    for (r=0; r<repet; r++, betavec += betas) {
      if (r==0 || L->betas_separate) {
	for (int i=0; i<totptsvdim; work[i++] = 0.0);
	//  for (int j=0; j<totptsvdim;j ++) printf("j=%d p %10g d %10g w %10g\n ", j, pres[j], data[j], work[j]);
	for (int i=0; i<betas; i++) {
	  double beta = betavec[i];
	  //	    printf("i=%d beta = %10g %d\n", i, beta, nrow);
	  for (int j=0; j<nrow; j++, X++) {	
	    //	      printf("i=%d %d %10g %10g\n", i, betas, beta, *X);
	    work[j] += *X * beta;
	  }
	}
      }
      
      for (int j=0; j<nrow; j++, data++, pres++) {
	*pres += work[j];
	//    printf("j=%d p %10g d %10g w %10g\n ", j, *pres, *data, work[j]);
      }
    }
  }
  
  // for (int j=0; j<nrow; j++, data++, pres++) printf("%10g %10g\n", *pres, *data);
  if (delete_work) FREE(work);
}


void gauss_predict(model *predict, model *Cov, double *v) {

  model
    *cov = Cov->key != NULL ? Cov->key : Cov->sub[0],
    //    *sub = cov->key == NULL ? cov->sub[0] : cov->key,
    *pred_cov = predict->key != NULL ? predict->key : predict->sub[0];
  
  //PMI(predict);
  assert(isnowProcess(cov) && 
	 isnowVariogram(cov->key == NULL ? cov->sub[0] : cov->key) &&
	 equalsnowInterface(predict));
 likelihood_storage *L = cov->Slikelihood;	
  assert(L != NULL);
  listoftype *datasets = L->datasets;
  int 
    err = NOERROR,
    Exterr = NOERROR,
    store = GLOBAL.general.set,
    sets = GET_LOC_SETS(cov);
  if (sets != 1) ERR("only one data set allowed when predicting");
  assert(cov->Ssolve != NULL);

  //  APMI(predict);
  assert(predict->Sextra != NULL);
 
  location_type 
    *loc = Loc(cov),
    *pred_loc = Loc(predict);   
  assert(pred_loc != NULL && loc != NULL);
  GLOBAL.general.set = 0;
  int 
    spatialdim = pred_loc->spatialdim,
    timespacedim = GetLoctsdim(predict),
    vdim = VDIM0,
    ncol =  NCOL_OUT_OF(datasets),    
    repet = ncol / vdim,
    totpts = Gettotalpoints(Cov),
    totptsvdim  = totpts * vdim,
    //
    totvdimrepet = totpts * vdim * repet,
    pred_tot = Gettotalpoints(predict),
    pred_totvdim = pred_tot * vdim;
  double 
    *residuals = NULL;

  if ((residuals = (double *) MALLOC(totvdimrepet * sizeof(double))) == NULL) {
    GLOBAL.general.set = store;
    XERR(ERRORMEMORYALLOCATION);
  }
  get_logli_residuals(cov, NULL, residuals, LOGLI_RESIDUALS);

 
  CovarianceMatrix(cov, L->C);  // tot x tot x vdim x vdim //sub -> cov 30.12.15
  
  if (GLOBAL.krige.ret_variance) {
    //
    GERR("kriging variance cannot be returned, currently");
     // end !data_nas
  } else {  
    ONCE_NEW_COV_STORAGE(predict, extra); // sicherheitshalber
    TALLOCCOV_NEW(predict, Sextra, Cx0, vdim*vdim*pred_tot, a1, 1);//da eh egal
    TALLOCCOV_NEW(predict, Sextra, y0, timespacedim, a2, XSIZE);
    gauss_trend(predict, cov, v, 0); // not sub
    assert(cov->Ssolve != NULL);

    double *ResiWithoutNA, *Ccur;
    int atonce, endfor, 
      notnas = NA_INTEGER;
    if (L->data_nas[GLOBAL.general.set]) {
      atonce = 1;
      Ccur = L->Cwork;
      ResiWithoutNA = L->Xwork;
    } else {
      notnas = totptsvdim;
      atonce = repet;
      Ccur = L->C;
      ResiWithoutNA = residuals;
    }
    endfor = repet / atonce;
 
    for (int r=0; r<endfor; r++) {
      if (L->data_nas[GLOBAL.general.set]) {
	double *resi = residuals + r * totptsvdim;
	notnas = matrixcopyNA(ResiWithoutNA, resi, resi, totptsvdim, 1, 0);
	SqMatrixcopyNA(Ccur, L->C, resi, totptsvdim);
      }

      //printf("********************\n");
  
      Exterr = Ext_solvePosDef(Ccur, notnas, 
				true, // except for intrinsic Kriging
        	//       as if ordinary Kriging (p 167; 265 in Chiles
				ResiWithoutNA, atonce,
			     NULL, cov->Ssolve);

      if (Exterr != NOERROR) goto ErrorHandling;    
      
      double inc[MAXLILIGRIDDIM+1], x[MAXLILIGRIDDIM+1],
	xstart[MAXLILIGRIDDIM+1],
	*pt_locx = loc->x,
	*pty = y0;
      assert(pred_loc != NULL && loc != NULL);
      int d, gridlen[MAXLILIGRIDDIM + 1], end[MAXLILIGRIDDIM + 1],
	start[MAXLILIGRIDDIM + 1], nx[MAXLILIGRIDDIM + 1],
	Tidx = 0,
	i_row = 0,
	Tlen = loc->Time ? loc->T[XLENGTH] : 0,
	pred_timespacedim = GetLoctsdim(predict),
	tsxdim = timespacedim,
	cumgridlen[MAXLILIGRIDDIM + 1];
      if (timespacedim > MAXLILIGRIDDIM) GERR("dimension too large");
      if (loc->grid) {
	STANDARDSTART_X;
	pty = x;
      }
      for (int i=0; i<totpts; i++) {
	if (!loc->grid) {	    
	  if (loc->Time) {
	    MEMCOPY(pty, pt_locx, spatialdim * sizeof(double));
	    pty[spatialdim] = loc->T[XSTART] + loc->T[XSTEP] * Tidx;
	    Tidx = (Tidx + 1) % Tlen;
	    if (Tidx == 0) pt_locx += spatialdim;
	  } else {
	    pty = pt_locx + spatialdim * i;
	  }
	}
	
	if (pred_loc->grid) {
	  for(d=0; d<pred_timespacedim; d++) 
	    pred_loc->ygr[d][XSTART] = pty[d];
	} else {
	  MEMCOPY(pred_loc->y, pty, pred_timespacedim * sizeof(double));
	  if (pred_loc->Time) pred_loc->T[0] = pty[spatialdim];
	}	 
	assert(pred_cov -> calling != NULL);
	DefList[MODELNR(pred_cov)].covariance(pred_cov, Cx0);//predtot x vdim^2
	double *ptCx0 = Cx0;
	for (int j=0; j<vdim; j++, ptCx0 += pred_totvdim) { // influence of the different variates of position i
	  int m=0;
	  for (int s=0; s<atonce; s++) {	
	    double factor = residuals[totptsvdim * s + totpts * j + i];
	    if (ISNAN(factor)) {
	      m += pred_totvdim;
	      continue;
	    }
	    for (int k=0; k<pred_totvdim; k++) { // alle Y_i (ohne repet)
	      v[m++] += factor * ptCx0[k];	      
	    }
	  }
	}
	if (loc->grid) { // !! geklammert!!!
	  STANDARDINKREMENT_X;
	}
      }
    } // end !data_nas
    END_TALLOC_NEW(a1);
    END_TALLOC_NEW(a2);
  } // ! krige var
 
  //  for (int i=0; i<10; i++) printf("%10g\n", v[i]);
    


 ErrorHandling :
  GLOBAL.general.set = store;
  FREE(residuals);
    
  if (Exterr != NOERROR) {
    if (Exterr == ERRORM) {
      Ext_getErrorString(predict->err_msg);
      ERR(predict->err_msg);
    } else err = Exterr;
  }

  if (err != NOERROR) XERR(err);
  
}


//////////////////////////////////////////////////////////////////////
// Gauss process

SEXP simple_residuals(SEXP model_reg){
  START_GAUSS_INTERFACE;
  likelihood_info *info = &(L->info);

  assert(!L->ignore_trend);
  
  listoftype *datasets = L->datasets;
  if (info->globalvariance && info->pt_variance != NULL) 
    *(info->pt_variance) = 1.0;
 
  int i, j, 
    fx_notnas = 0,
    sets = GET_LOC_SETS(calling),  
    err = NOERROR,
    Exterr = NOERROR,
    vdim = VDIM0,
    betatot = L->betas[L->fixedtrends],
    betaSq = betatot * betatot;
  //if (VDIM1 != 1) BUG;

  for (i=0; i<L->fixedtrends; i++) 
    fx_notnas += (!L->nas_fixed[i]) * (L->betas[i+1] - L->betas[i]);
  
  for (i=0; i<betaSq; L->XtX[i++] = 0.0);
  for (i=0; i<betatot; L->XCY[i++] = 0.0);
  for (GLOBAL.general.set = 0; GLOBAL.general.set < sets;
       GLOBAL.general.set++) {
    //  Gettotalpoints, Loc etc haengen alle von GLOBAL.general.set ab !
    int k, m, 
      ncol =  NCOL_OUT_OF(datasets),
      repet = ncol / vdim,
      nrow = NROW_OUT_OF(datasets),
      ndata = nrow * ncol,
      totpts = Gettotalpoints(calling),
      totptsvdim = totpts * vdim;
    assert(nrow == totpts); 
    double 
      *Xdata = L->X[GLOBAL.general.set] + betatot * totptsvdim;
  
    if (L->dettrend_has_nas || L->nas_boxcox) {
      MEMCOPY(Xdata, SET_OUT_OF(datasets), ndata * sizeof(double));
      if (L->nas_boxcox) 
	boxcox_trafo(P(GAUSS_BOXCOX), vdim, Xdata, nrow, repet);
    }
   
    for (i=0; i<L->fixedtrends; i++) {
      if (L->nas_fixed[i]) {
	double *Lx = L->X[GLOBAL.general.set] + L->betas[i] * totptsvdim,
	  end = totptsvdim * (L->betas[i+1] - L->betas[i]);
	for (m=0; m<end; Lx[m++] = 0.0);
      }
    }

    if (L->random > 0) BUG; // to do
 
    double *Xcur,
      *dataWithoutNA = NULL;
    int atonce, endfor,
      notnas = NA_INTEGER;
    if (L->data_nas[GLOBAL.general.set]) {
      atonce = 1;
      Xcur = L->Xwork;
    } else {
      atonce = repet;
      Xcur = L->X[GLOBAL.general.set];
      dataWithoutNA = Xdata;
      notnas = totptsvdim;
    }
    endfor = repet / atonce;
    
 
    for (int r=0; r<endfor; r++) {
      if (L->data_nas[GLOBAL.general.set]) {
	double *datacur = Xdata + r * totptsvdim;
	notnas = matrixcopyNA(Xcur, L->X[GLOBAL.general.set], 
			      datacur, totptsvdim, betatot, atonce);
	dataWithoutNA = Xcur + notnas * betatot;
      }

      if (L->fixedtrends) {
	matmulttransposed(Xcur, Xcur, L->XitXi, notnas, betatot, betatot);
	// bei simple residuals muessen alle spalten und Zeilen
	// mit nas_fixed raus:
	for (int i1=k=0; i1<L->fixedtrends; i1++) { 
	  if (L->nas_fixed[i1] == 0) {
	    int end_i = L->betas[i1 + 1];
	    for (int i2=L->betas[i1]; i2<end_i; i2++) {
	      for (int j1=0; j1<L->fixedtrends; j1++) { 
		if (L->nas_fixed[j1] == 0) {
		  int end_j = L->betas[j1 + 1];
		  for (int j2=L->betas[j1]; j2<end_j; j2++) {
		    int idx = i2 * betatot + j2;
		    L->XtX[k++] += atonce * L->XitXi[idx];
		  } // j2
		}
	      } // j1
	    } // i2
	  }
	} // i1
	
	MEMCOPY(L->sumY, dataWithoutNA, sizeof(double) * notnas);
	k = notnas;
	for (int rr=1; rr<atonce; rr++) {
	  for(int d=0; d<notnas; d++) {
	    L->sumY[d] += dataWithoutNA[k++];
	  }
	}
	
	double *dummy = L->betavec;
	matmult(L->sumY, Xcur, dummy, 1, notnas, betatot);
	for (i=0; i<betatot; i++) L->XCY[i] += dummy[i];
      }
    } // end !data_nas
  } // Global.general.set

  if (L->fixedtrends) {
    // XCY, etc NUR FUER FIXED TREND
    double *beta = L->betavec; // !!!
    for (i = 0; i<betatot;  beta[i++] = 0.0);

    if (fx_notnas > 0) {
      int k=0;
      for (int j1=0; j1 < L->fixedtrends; j1++) { 
	if (L->nas_fixed[j1] == 0) {
	  int end=L->betas[j1 + 1];
	  for (int j2=L->betas[j1]; j2<end ; j2++) {
	    beta[k++] = L->XCY[j2];
	  } // j2
	}
      }
      assert(cov->Ssolve != NULL);
      Exterr = Ext_solvePosDef(L->XtX, fx_notnas, true, beta, 1, NULL, 
				cov->Ssolve);
      if (Exterr != NOERROR) goto ErrorHandling;    
     
      if (fx_notnas < betatot) {
	int b = L->fixedtrends - 1;
	i = betatot - 1; 
	j = fx_notnas - 1;
	while (b>=0 && i >= 0) {
	  int bi = L->betas[b + 1] - L->betas[b];
	  if (L->nas_fixed[b]) {
	    for (k=0; k<bi; k++) beta[i--] = 0.0;
	  } else for (k=0; k<bi; k++) {
	      assert(i >= 0 && j >= 0);
	      beta[i--] = beta[j--];
	    }
	  assert(i>=j);
	}
        assert (i == 0 && j == 0 && b == 0);
      } else assert(fx_notnas == betatot);
    }    
  }
  
 ErrorHandling:
  GLOBAL.general.set = store;
  if (Exterr != NOERROR) {
    if (Exterr == ERRORM) {
      Ext_getErrorString(cov->err_msg);
      ERR(cov->err_msg);
    } else err = Exterr;
  }
  if (err != NOERROR) XERR(err);

  return get_logli_residuals(model_reg);
}


//////////////////////////////////////////////////////////////////////
// Gauss process

// likelihood / density
void gaussprocessDlog(double VARIABLE_IS_NOT_USED *x, model *cov, 
		      double *v){  
  assert(isnowProcess(cov));
  model *calling = cov->calling;
  //*sub = cov->key == NULL ? cov->sub[0] : cov->key;  
  assert(calling != NULL && DefList[CALLINGNR].cov == likelihood);
  likelihood_storage *L = cov->Slikelihood;
  assert(L != NULL);
  listoftype *datasets = L->datasets;
  assert(!L->ignore_trend);
  likelihood_info *info = &(L->info);
  if (info->globalvariance && info->pt_variance != NULL) 
    *(info->pt_variance) = 1.0;
  double 
    proj = 0.0,
    YCinvY = 0.0,
    logdettot = 0.0;  
  int i, j, repet,
    sets = GET_LOC_SETS(calling),  
    err = NOERROR,
    Exterr = NOERROR,
    variables = 0,
    vdim = VDIM0,
    store = GLOBAL.general.set,
    betatot = L->betas[L->fixedtrends],    
    betaSq = betatot * betatot;
  // PMI(cov);
  //  if (VDIM1 != 1) BUG;

  *v = 0.0;
  for (i=0; i<betaSq; L->XtX[i++] = 0.0);
  for (i=0; i<betatot; L->XCY[i++] = 0.0);
  for (GLOBAL.general.set = 0; GLOBAL.general.set < sets;
       GLOBAL.general.set++) {
    //  Gettotalpoints, Loc etc haengen alle von GLOBAL.general.set ab !    
    int k, m, p,
      ncol =  NCOL_OUT_OF(datasets);    
    assert(!L->betas_separate || sets == 1);
    repet = ncol / vdim;
    int nrow = NROW_OUT_OF(datasets), 
      ndata = nrow * ncol,
      totpts = Gettotalpoints(calling),
      totptsvdim = totpts * vdim;
    double 
      logdet,
      *val = L->C,
      *YhatWithoutNA = L->YhatWithoutNA[GLOBAL.general.set],
      *Xdata = L->X[GLOBAL.general.set] + betatot * totptsvdim;
  
    if (L->dettrend_has_nas || L->nas_boxcox) {
      MEMCOPY(Xdata, SET_OUT_OF(datasets), ndata * sizeof(double));
 
      if (L->nas_boxcox) {
	for (p = m = 0; m < repet; m++) {
	  for (k = 0; k < totptsvdim; k++) {
	    Xdata[p++] -= YhatWithoutNA[k];
	  }
	}
      }
    
      for (i=0; i<L->dettrends; i++) {
	if (L->nas_det[i]) {
	  FctnIntern(cov, L->cov_det[i], L->cov_det[i], val, false);	    
	  for (p = m = 0; m < repet; m++) {
	    for (k = 0; k < totptsvdim; k++) Xdata[p++] -= val[k];
	  }
	}
      }

     if (L->nas_boxcox) 
       boxcox_trafo(P(GAUSS_BOXCOX), vdim, Xdata, nrow, repet);

    } // (L->dettrend_has_nas || L->nas_boxcox)

    for (i=0; i<L->fixedtrends; i++) {
      //      printf("i=%d %d\n", i, L->nas_fixed[i]);
      if (L->nas_fixed[i]) {
	FctnIntern(cov, L->cov_fixed[i], L->cov_fixed[i],
		   L->X[GLOBAL.general.set] + L->betas[i] * totptsvdim,
		   false); // "Lueckentext" auffuelen !!
      }
    }

    if (L->random > 0) BUG; // to do
    if (info->globalvariance && info->pt_variance != NULL)
      *(info->pt_variance) = 1.0;
    CovarianceMatrix(cov, L->C); //sub -> cov 30.12.15
    // printf("\nC =  ");for (i=0; i<totptsvdim * totptsvdim; i++)printf("%10g ", L->C[i]); printf("\n"); 
      //      BUG;   APMI(sub);

    //for (int i=0; i<9; i++) printf("%10g ", L->C[i]); printf("\n");
   //   APMI0(cov);
   // BUG;
   //    printf("%10g %10g %10g %10g %10g\n", Xdata[0], Xdata[1], Xdata[2], Xdata[3], Xdata[4]); 		      

    double *Xcur, *Ccur, 
      *dataWithoutNA = NULL;
    int atonce, endfor, XYcols, 
      notnas = NA_INTEGER;
    if (L->data_nas[GLOBAL.general.set]) {
       atonce = 1;
      Ccur = L->Cwork;
      Xcur = L->Xwork;
    } else {
      notnas = totptsvdim;
      atonce = repet;
      Ccur = L->C;
      Xcur = L->X[GLOBAL.general.set];
      dataWithoutNA = Xdata;
    }
    endfor = repet / atonce;
    XYcols = betatot + atonce;

    // printf("%10g\n", L->X[GLOBAL.general.set][0]); BUG; 

    for (int r=0; r<endfor; r++) {
      if (L->data_nas[GLOBAL.general.set]) {
	double *datacur = Xdata + r * totptsvdim;
	notnas = matrixcopyNA(Xcur, L->X[GLOBAL.general.set], 
			      datacur, totptsvdim, betatot, atonce);
	SqMatrixcopyNA(Ccur, L->C, datacur, totptsvdim);
	dataWithoutNA = Xcur + notnas * betatot;
      }
 
      variables += notnas * atonce;
      MEMCOPY(L->CinvXY, Xcur, notnas * XYcols * sizeof(double));
      assert(cov->Ssolve != NULL);

      //BUG;
      //printf("XXX= %d\n", 4);

      Exterr = Ext_solvePosDef(Ccur, notnas, 
			      true,  // except for negative definite function
			      L->CinvXY, XYcols,
			      &logdet, cov->Ssolve);
      if (Exterr != NOERROR) goto ErrorHandling;  
      logdettot += logdet * atonce;
      double *CinvY = L->CinvXY + betatot * notnas;
      int end_i = notnas * atonce;
      
      // for (i=0; i<end_i; i++) YCinvY += CinvY[i] * dataWithoutNA[i];
      YCinvY += SCALAR_P(CinvY, dataWithoutNA, end_i);

      if (L->fixedtrends) {
	if (L->betas_separate) {
	  assert(sets == 1);
	  matmulttransposed(Xcur, L->CinvXY, L->XtX, notnas, betatot, betatot);
	  matmulttransposed(L->CinvXY, dataWithoutNA, L->XCY, notnas,
			    betatot,atonce);
	} else {
	  //BUG;
	  matmulttransposed(Xcur, L->CinvXY, L->XitXi, notnas, betatot,betatot);
	  for (i=0; i<betaSq; i++) L->XtX[i] += atonce * L->XitXi[i];

	  // printf("atonce=%d %10g xcur=%10g %10g\n", atonce, L->XitXi[0], Xcur[0], L->CinvXY[0]);
	  
	  MEMCOPY(L->sumY, dataWithoutNA, sizeof(double) * notnas);

	  int kk = notnas;
	  for (int rr=1; rr<atonce; rr++) {
	    for(int d=0; d<notnas; d++) {
	      L->sumY[d] += dataWithoutNA[kk++];
	    }
	  }
	  
	  double *dummy = L->betavec;
	  matmult(L->sumY, L->CinvXY, dummy, 1, notnas, betatot);	  
	  for (i=0; i<betatot; i++) L->XCY[i] += dummy[i];
	}
      } // fixedtrends	
    } // r, repet	  
  } // Global.general.set
  
  if (L->fixedtrends) {
    // XCY, etc NUR FUER FIXED TREND
    assert(!L->betas_separate || sets == 1);
    GLOBAL.general.set = 0;  // betas_separate only
    double *beta = L->betavec; // !!!
    int 
      all_betatot = betatot,
      ncol =  NCOL_OUT_OF(datasets);
    repet = ncol / vdim;
  
    if (L->betas_separate) {
      all_betatot *= repet;
    }
  
   
    MEMCOPY(beta, L->XCY, sizeof(double) * all_betatot);
    assert(!L->betas_separate || sets == 1);
    assert(cov->Ssolve != NULL);

    //  for (int h=0; h<betatot * betatot; h++)
      //printf("%10g ", L->XtX[h]);
      //    printf("\n>> %d %d %10g\n",  betatot, L->betas_separate, beta[0]);

    Exterr = Ext_solvePosDef(L->XtX, betatot, true, beta, 
		     L->betas_separate ? repet : 1, NULL, 
		     cov->Ssolve);
    if (Exterr != NOERROR) goto ErrorHandling;  
       //   printf("%.50s AC %10g %10g %10g %10g %d done\n", NAME(cov), beta[0], beta[1], beta[2], beta[3], all_betatot); BUG;
    for (j=0; j<all_betatot; j++) {
      proj += L->XCY[j] * beta[j];
    }      
    MEMCOPY(v + 1 + info->globalvariance, beta, all_betatot * sizeof(double));

    if (L->betas_separate) for(i=0; i<betatot; i++) *(L->where_fixed[i]) =RF_NA;
    else for (i=0; i<betatot; i++) *(L->where_fixed[i]) = beta[i];
  }

  *v = -0.5 * (logdettot + variables * LOG(TWOPI));

  //  printf(">> v =%10g %10g #=%d 2pi^=%10g sq=%10g %10g %d\n", *v, 0.5 *logdettot, variables , -0.5 * variables * LOG(TWOPI), YCinvY, proj, info->globalvariance); 
 

  double delta;
  delta = YCinvY - proj;
  if (info->globalvariance) {
    //   printf("%10g delta=%10g %d %10g \n", *v, delta, variables, 0.5 * variables * LOG(delta)); BUG;
    
    v[1] = delta / variables;
    *v -= 0.5 * variables * (1 + LOG(v[1]));
     if (info->pt_variance != NULL) *(info->pt_variance) = v[1];
  } else {
    //    printf("delta = %10g, %10g YCY=%10g proj=%10g\n", delta, -0.5 * delta, YCinvY,proj);
    *v -= 0.5 * delta;
  }
  //  printf("*v=%10g %10g %d\n", *v, P(GAUSS_BOXCOX)[0], info->globalvariance);
  
  if (R_FINITE(P(GAUSS_BOXCOX)[0])) {
    //printf("*v=%10g %10g %d\n", *v, P(GAUSS_BOXCOX)[0], info->globalvariance);BUG;

    for (j=0; j<vdim; j++) {
      double lambda  = P(GAUSS_BOXCOX)[2 * j];
      if (lambda != 1.0) {
	double 
	  lambdaM1 = lambda - 1.0,
	  mu = P(GAUSS_BOXCOX)[2 * j + 1],
	  sum = 0.0;
	for (GLOBAL.general.set = 0; GLOBAL.general.set < sets;
	     GLOBAL.general.set++) {
	  int 
	    ncol =  NCOL_OUT_OF(datasets),   
	    nrow = NROW_OUT_OF(datasets), 
	    totptsvdim = nrow * vdim;
	  repet = ncol / vdim;
	  double 
	    *data = SET_OUT_OF(datasets) + nrow * j;
	  for (int r=0; r < repet; r++, data += totptsvdim) {
	    for (int ii=0; ii<nrow; ii++) {
	      double dummy = data[ii];
	      if (!ISNAN(dummy))
		sum += LOG(dummy + mu);
	    }
	  }
	}
	*v += lambdaM1 * sum;
      }
    }
  }
  //  printf("      *v=%10g %10g\n", *v, P(GAUSS_BOXCOX)[0]);

  
  //  APMI(cov);  BUG;

 ErrorHandling:
  GLOBAL.general.set = store;

  if (Exterr != NOERROR) {
    if (Exterr == ERRORM) {
      Ext_getErrorString(cov->err_msg);
      ERR(cov->err_msg);
    } else err = Exterr;
  }

  if (err != NOERROR) XERR(err);
   
}

// v ok
void PutGlblVar(int *reg, double *var) {
  assert (*reg >= 0 && *reg <= MODEL_MAX);
  model *cov = KEY()[*reg];
  assert(cov != NULL && equalsnowInterface(cov));
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key; 
  likelihood_storage *L;
  if (sub == NULL || !isnowProcess(sub) || (L = sub->Slikelihood) == NULL) BUG;
  likelihood_info *info = &(L->info);
  if (info->pt_variance != NULL) {
    // == NULL may happen if globalvariance is true by user, but not
    // reflected in the model
    *(info->pt_variance) = *var;
  }
}

SEXP get_likeliinfo(SEXP model_reg) {
  START_GAUSS_INTERFACE;
  likelihood_info *info = &(L->info);
#define nn 5
  const char *names[nn] = 
    {"betas", "betanames", 
     "estimate_variance", "sum_not_isna_data", 
    "betas_separate"};
  listoftype *datasets = L->datasets;
  int k, 
    sets = GET_LOC_SETS(cov),  
    sum_not_isna_data = 0,
    betas = L->betas[L->fixedtrends];
  SEXP namevec, ans, betanames;

  for (GLOBAL.general.set = 0; GLOBAL.general.set < sets; GLOBAL.general.set++){
    int
      ncol = NCOL_OUT_OF(datasets),
      nrow = NROW_OUT_OF(datasets),
      ndata = nrow * ncol;
    sum_not_isna_data += ndata - L->data_nas[GLOBAL.general.set];
  }
  
  PROTECT(ans = allocVector(VECSXP, nn));
  PROTECT(namevec = allocVector(STRSXP, nn));
  for (k=0; k<nn; k++) SET_STRING_ELT(namevec, k, mkChar(names[k]));

  //  printf("betas %d\n", betas);
  PROTECT(betanames = allocVector(STRSXP, betas));
  for (k=0; k<betas; k++) {
    //printf("%d %.50s<\n", k, L->betanames[k]);
    SET_STRING_ELT(betanames, k, mkChar(L->betanames[k]));
  }

  k = 0;  
  SET_VECTOR_ELT(ans, k++, ScalarReal(betas));
  SET_VECTOR_ELT(ans, k++, betanames);
  SET_VECTOR_ELT(ans, k++, ScalarLogical(info->globalvariance));
  SET_VECTOR_ELT(ans, k++, ScalarInteger(sum_not_isna_data));
  SET_VECTOR_ELT(ans, k++, ScalarLogical(L->betas_separate));
  assert(k == nn);


  setAttrib(ans, R_NamesSymbol, namevec);
  UNPROTECT(3);
 
  GLOBAL.general.set = store;
  return ans;
}


int countbetas(model *cov, double ***where) {
  defn *C = DefList + COVNR;
  int i,
    sum = 0,
    kappas = C->kappas;

  for(i=0; i<kappas; i++) {
    if (cov->kappasub[i] != NULL) continue;
    if (isnowTrendParam(cov, i)) {
      assert(C->kappatype[i] == REALSXP);
      if (!PisNULL(i)) {
	int total = cov->ncol[i] * cov->nrow[i];
	double *p = P(i);
	if (ISNAN(p[0])) {
	  sum += total;
	  for (int j=0; j<total; j++) { // 0 for 'where'
	    if (!ISNAN(p[j])) ERR("trend parameters must be all NA or none");
	    if (where != NULL) {
	      **where = p + j;
	      (*where)++;	     
	    }	    
	  }
	} else {
	  for (int j=1; j<total; j++) {
	    if (ISNAN(p[j])) ERR("trend parameters must be all NA or none");
	  } // j
	} // not (ISNAN(p[0])) 
      } // (!PisNULL(i)) 
    } // is Trend
  } // for kappa
  //  printf("sum = %d\n", sum);
  return sum;
}


void GetBeta(model *cov, likelihood_storage *L, int *neffect)  {
  if (isnowProcess(cov)) {
    int nas = ((bool) ISNA((P(GAUSS_BOXCOX)[0]))) +
      ((bool) ISNA((P(GAUSS_BOXCOX)[1])));
    assert(cov->key == NULL);
    assert(!PisNULL(GAUSS_BOXCOX));
    if (nas > 0) (*neffect)++;
    GetBeta(cov->sub[0], L, neffect);
    return;
  }

  int i,
    n = COVNR == PLUS ? cov->nsub : 1;
  if (*neffect >= MAX_LIN_COMP) ERR("too many linear components");
  bool 
    plus = COVNR == PLUS;
  likelihood_info *info = &(L->info);

  //  printf("entering getbeta\n");  pmi(cov, 0);

  for (i=0; i<n; i++) {
    model *component = plus ? cov->sub[i] : cov;
    //   pmi(component, 0);
    if (MODELNR(component) == PLUS) {
      //   printf("recursive\n");
      GetBeta(component, L, neffect);
      //  printf("end recursive\n");
       continue;
    }

    // printf("%.50s %d %ld\n", NAME(cov), info->effect[*neffect], (long int) L);
    
    //  if (L->effect[*neffect] > LastMixedEffect) { } else 

    //printf("neffect %d %d\n", *neffect, info->effect[*neffect]);
    if (info->effect[*neffect]  == DetTrendEffect) { 
      L->cov_det[L->dettrends++] = component;
    } else if (info->effect[*neffect] == FixedTrendEffect) {
      int b = UNSET;
      L->betas[L->fixedtrends + 1] = L->betas[L->fixedtrends];
      L->cov_fixed[L->fixedtrends++] = component;
      if (MODELNR(component) == MULT) {
	for (int j=0; j<component->nsub; j++) {
	  if ((b = countbetas(component->sub[j], NULL))  > 0) {
	    break;
	  }
	}
      } else b = countbetas(component, NULL);	 

      if (b <= 0) {
	(*neffect)++;
	continue;
      }
      //ERR("fixed effect expected, but corresponding NA not found");
      int base =  L->betas[L->fixedtrends];
      L->betas[L->fixedtrends] += b;  
      if (L->maxbeta < b) L->maxbeta = b;

      model *comp = component;
      if (MODELNR(comp) == MULT) {
	for (int ii=0; ii<comp->nsub; ii++)
	  if (MODELNR(comp->sub[0]) == CONST && 
	      ISNA(PARAM0(comp->sub[0], CONST_C))){
	    comp = comp->sub[ii==0 && comp->nsub >= 2];
	    break;
	  }
      }
      if (isDollar(comp)) comp = comp->sub[0];
      char abbr[LENMSG];
      int bytes = (1 + GLOBAL.fit.lengthshortname) * sizeof(char);
      //      printf("%d %.50s %.50s\n", *neffect, NAME(component), NAME(comp));
      Abbreviate(NAME(comp), abbr);
      //      printf("%.50s<\n", abbr);
      //  if (L->betanames[*neffect] != NULL) {printf("%.50s %d %.50s\n", NAME(cov),  *neffect, L->betanames[*neffect]); crash();} else  printf("%.50s %d NULL\n", NAME(cov),  *neffect); 
      
      assert(L->betanames[base] == NULL);
      if (b == 1) {
	L->betanames[base] = (char*) MALLOC(bytes);
	SPRINTF(L->betanames[base], "%.50s", abbr);
      } else {
	for (int ii=0; ii<b; ii++) {
	  L->betanames[base + ii] = (char*) MALLOC(bytes);
	  SPRINTF(L->betanames[base + ii], "%.50s.%d", abbr, ii);
	}
      }

      //printf("----> %d %.50s\n", base, L->betanames[base]);

    } 
    (*neffect)++;
  }
}




void GetBeta(model *cov, likelihood_storage *L, int *neffect,
		 double ***where)  {
  if (isnowProcess(cov)) {
    int nas = ISNA((P(GAUSS_BOXCOX)[0])) + ISNA((P(GAUSS_BOXCOX)[1]));
    assert(cov->key == NULL);
    assert(!PisNULL(GAUSS_BOXCOX));
    if (nas > 0) (*neffect)++;
    GetBeta(cov->sub[0], L, neffect, where);
    return;
  }
  int i,
    n = COVNR == PLUS ? cov->nsub : 1;
  bool 
    plus = COVNR == PLUS;
  likelihood_info *info = &(L->info);
  for (i=0; i<n; i++) {
    model *component = plus ? cov->sub[i] : cov;
    if (MODELNR(component) == PLUS) {
      GetBeta(component, L, neffect, where);
      continue;
    }
    if (info->effect[*neffect] == FixedTrendEffect) {
      if (MODELNR(component) == MULT) {
	for (int j=0; j<component->nsub; j++) {
	  if ((countbetas(component->sub[j], where)) > 0) {
	    break;
	  }
	}
      } else countbetas(component, where);	  
    }
    (*neffect)++;
  }
}


#define MAX_TOTALPTS 10000
int struct_gauss_logli(model *cov) {
  assert(isnowProcess(cov));
  model 
    *calling = cov->calling,
    *sub = cov->key != NULL ? cov->key : cov->sub[0];
  location_type *loc = Loc(cov);
  int dim = ANYDIM;
  // int dim = ANYDIM;

  if (false) { // Beschleunigung, insb. im Genetik-Bereich
    // gleich ganz allgemein implementieren??
    if (isCartesian(SUB)) {
      if (isXonly(SUB)) {
	if (isIsotropic(SUB)) {
	  if (Gettotalpoints(cov) < MAX_TOTALPTS) {
	    // nur distances abspeichern
	    // in ownloc
	  }
	} else {
	  if (Gettotalpoints(cov) <
	      MAX_TOTALPTS / SQRT((double) dim)) {
	    // distance vectors abspeochen
	  }
	}
      }
    } else {
      NotProgrammedYet("non-cartesian systems")
      // to do 
    }
  }

  int i, k,  // info->neffect, 
    max_total_data_Sq, betatot,
    err = NOERROR,
    ne = 0,
    //n = SUBNR == PLUS ? sub->nsub : 1,
    store = GLOBAL.general.set,
    sets = GET_LOC_SETS(cov),  
    maxrepet = 0,
    maxndata = 0,
    maxtotpts = 0,
    vdim = VDIM0,
    two_vdim = 2 * vdim,  
    n = vdim * vdim,
    dummy_fixed=0, 
    dummy_det = 0;
  usr_bool 
    global_var = CALLINGNR == LIKELIHOOD_CALL 
    ? (usr_bool) PARAM0INT(calling, LIKELIHOOD_NA_VAR)
    : False;
  
  
  assert(calling != NULL);
  if (CALLINGNR != LIKELIHOOD_CALL && CALLINGNR != LINEARPART_CALL) BUG;
  NEW_STORAGE(likelihood);
  likelihood_storage *L = cov->Slikelihood;
  likelihood_NULL(L);
  likelihood_info *info = &(L->info);
  likelihood_info_NULL(info);
  L->ignore_trend = 
    CALLINGNR == LIKELIHOOD_CALL && PARAM0INT(calling, LIKELIHOOD_IGNORETREND); 
  // to do : betas_separate bereinigen && ignore_trend nicht alles
  // allocieren!

  L->betas_separate = 
    CALLINGNR == LIKELIHOOD_CALL && PARAM0INT(calling, LIKELIHOOD_BETASSEPARATE);
  if (L->betas_separate && sets > 1) SERR("separate estimation of the betas only for one data set possible.");

  L->nas_boxcox =  L->nas_boxcox_mu = 0;
  for (i=0; i<two_vdim; i++) {
    L->nas_boxcox += ISNA(P(GAUSS_BOXCOX)[i]);
    if (i % 2 == 1) L->nas_boxcox_mu += ISNA(P(GAUSS_BOXCOX)[i]);
  }

 
  // printf("%10g %10g\n ",  P(GAUSS_BOXCOX)[0], P(GAUSS_BOXCOX)[1]); BUG;

 
  //PMI(cov);
 
  // information about free variables: 
  // first: get effect and NAs   
  
  //  printf("globar %d %d\n", global_var, PARAM0INT(calling, LIKELIHOOD_NA_VAR));
  // assert(false);

  assert(cov->calling != NULL && cov->calling->calling == NULL);
  ONCE_NEW_COV_STORAGE(cov->calling, mle);
  err = SetAndGetModelInfo(cov, GLOBAL.fit.lengthshortname, // 3
			   LIKELI_NA_INTEGER, LIKELI_EXCLUDE_TREND, loc->xdimOZ,
			   global_var, info, original);
  if (err != NOERROR) goto ErrorHandling;

  // second: existence and "size" of known trend and fixed trend

  //  printf("\n\n\n\n .>>>>>>>> get Beta:\n");
  GetBeta(cov, L, &ne);
  // printf("\n\n\n\n .>>>>>>>> get Beta fodnex\n");

  if (L->fixedtrends + 1 < MAX_LIN_COMP)
    L->betas[L->fixedtrends + 1] = L->betas[L->fixedtrends]; 
  betatot = L->betas[L->fixedtrends];
  
  if (betatot > 0) {
    ne = 0;
    double **pt = L->where_fixed = 
      (double **) MALLOC( L->betas[L->fixedtrends] * sizeof(double**) );
    GetBeta(cov, L, &ne, &pt);  // pt wird hochgezaehlt
    //    printf("%d %d\n", betatot, pt - L->where_fixed );
    assert(pt == L->where_fixed + betatot);
  }
 
  // alloc_cov is need for both CovarianceMatrix and FctnInternal
  if (n < L->maxbeta) n = L->maxbeta;

  //printf("maxbeta %d %d %d\n", L->maxbeta, vdim, n); 

  if ((err = alloc_cov(cov, dim, n, 1)) != NOERROR) goto ErrorHandling;
  
  //  printf("betas = %d/%d %d;  %d %d %d\n", L->betas[ne], betatot, info->neffect, info->effect[0], info->effect[1], info->effect[2]);
  for(k=0; k<info->neffect; k++) {
    //    printf("k=%d %d \n", k, info->effect[k]);
    if (info->effect[k] == DetTrendEffect) {
      L->nas_det[dummy_det++] = info->nas[k];
    } else if (info->effect[k] == FixedTrendEffect) {
      L->nas_fixed[dummy_fixed++] = info->nas[k];
    } 
  }

  assert(dummy_det == L->dettrends && dummy_fixed == L->fixedtrends && 
	 0 == L->random);

  for (GLOBAL.general.set = 0; GLOBAL.general.set < sets; GLOBAL.general.set++){
    int totptsvdim = Gettotalpoints(cov) * vdim;
    if (totptsvdim > L->max_total_data) 
      L->max_total_data = totptsvdim;
  }
  L->X = (double **) CALLOC(sets, sizeof(double *));
  L->YhatWithoutNA = (double **) CALLOC(sets, sizeof(double *));
  for (i=0; i<sets; i++) L->X[i] = L->YhatWithoutNA[i] = NULL;
  L->sumY = (double *) MALLOC(L->max_total_data * sizeof(double)); 
  L->sets = sets;
  
  for (GLOBAL.general.set = 0; GLOBAL.general.set < sets; GLOBAL.general.set++){
    double 
       *val = L->sumY;	  // dummy  
    int 
      totptsvdim = Gettotalpoints(cov) * vdim;
    
    double *Yhat = L->YhatWithoutNA[GLOBAL.general.set] =
      (double*) CALLOC(totptsvdim, sizeof(double));
    
    for (i=0; i<L->dettrends; i++) {
      if (L->nas_det[i] == 0) {
  	FctnIntern(cov, L->cov_det[i], L->cov_det[i], val, false);   
	for (k = 0; k < totptsvdim; k++) Yhat[k] += val[k];
      } else L->dettrend_has_nas = true; 
    }

    //    for (k = 0; k < totptsvdim; k++) printf("%10g ", Yhat[k]); printf("\n"); BUG;
    
    if (L->random) 
      GERR("Currently the linear part of random effects cannot be determined");
  }
 
  if (CALLINGNR == LINEARPART_CALL) {
    if (L->fixedtrends) {
      for (GLOBAL.general.set = 0; GLOBAL.general.set < sets;
	   GLOBAL.general.set++){
	int tot = betatot * Gettotalpoints(cov) * vdim;
	L->X[GLOBAL.general.set] = (double*) CALLOC(tot, sizeof(double)); 
      }  
    }
  } else {
    listcpy(&(L->datasets),  PARAMLIST(calling, LIKELIHOOD_DATA), false);
    L->data_nas = (int *) MALLOC(sets * sizeof(int));
    listoftype *datasets = L->datasets;
       
    for (GLOBAL.general.set = 0; GLOBAL.general.set < sets; 
	 GLOBAL.general.set++){
      double 
	*data = SET_OUT_OF(datasets);	    
      int
	ncol = NCOL_OUT_OF(datasets),
	repet = ncol / vdim,
	nrow = NROW_OUT_OF(datasets),
	data_nas = 0,
	ndata = nrow * ncol;
 
      if (nrow != Gettotalpoints(cov)) BUG;
      if (repet > maxrepet) maxrepet = repet;
      if (nrow  > maxtotpts) maxtotpts = nrow;
      if (ndata > maxndata) maxndata = ndata;
          
      // any_data_nas = false;  
      for (k = 0; k<ndata; data_nas += ISNA(data[k++]));
      L->data_nas[GLOBAL.general.set] = data_nas;
      L->data_has_nas |= data_nas > 0;

      if (L->nas_boxcox_mu) {
	bool warn_boxcox = true;
	int z = 0;
	for (int r=0; r<repet; r++) {
	  for(i=0; i<vdim; i++) {
	    double negmu = - GLOBAL.fit.BC_lambdaLB[2 * i + 1];
	    for(k=0; k<nrow; k++, z++) {
	      if (!ISNA(data[z]) && data[z] <= negmu) {
		if (warn_boxcox) {
		  warning("data values do not ensure positive values with given lower bounds for the second parameter in the Box-Cox transformation, hence the lower bound has been increased.");
		  warn_boxcox = false;
		}
		negmu = data[k] - 1e-10;
		GLOBAL.fit.BC_lambdaLB[2 * i + 1] = -negmu;	  
	      }
	    } // k, nrow
	  } // i, vdim
	} // r, repet
      } else if (R_FINITE(P(GAUSS_BOXCOX)[0])) {
	int z = 0;
	bool warn_boxcox = false;
	for (int r=0; r<repet; r++) {
	  for(i=0; i<vdim; i++) {
	    double negmu = - GLOBAL.fit.BC_lambdaLB[2 * i + 1];
	    for(k=0; k<nrow; k++, z++) {
	      if ((warn_boxcox = !ISNA(data[z]) && data[z] <= negmu)) break;
	    }
	  }
	}
	if (warn_boxcox) warning("data values do not ensure positive values in Box-Cox transformation, hence likely errors will occur.");
      }
    }
 
    // provide necessary space and account for known information
    max_total_data_Sq = L->max_total_data * L->max_total_data;
    long totdata_bytes = max_total_data_Sq * sizeof(double);

  
    L->C = (double *) MALLOC(totdata_bytes);    
    if (L->data_has_nas) {
      L->Cwork = (double *) MALLOC(totdata_bytes);
      int tot = L->max_total_data + betatot * Gettotalpoints(cov);
      L->Xwork = (double*) MALLOC(tot * sizeof(double));
    }
    if (L->fixedtrends) {
      L->XitXi =  (double*) MALLOC(sizeof(double) * betatot * 
				   (betatot > maxrepet ? betatot : maxrepet)); 
      L->XtX = (double*) CALLOC(betatot * betatot, sizeof(double)); 
      int all_betatot = betatot;
      if (L->betas_separate) all_betatot *= maxrepet;
      L->XCY = (double *) MALLOC(all_betatot * sizeof(double)); 
      L->betavec = (double *) MALLOC(all_betatot * sizeof(double)); 
    }

    L->CinvXY = 
      (double*) MALLOC(sizeof(double) * (betatot* maxtotpts* vdim + maxndata)); 
    
    for (GLOBAL.general.set = 0; GLOBAL.general.set<sets; GLOBAL.general.set++){
      double 
	*data = SET_OUT_OF(datasets),
	*YhatWithoutNA = L->YhatWithoutNA[GLOBAL.general.set];	    
      int  m, p,
	ncol = NCOL_OUT_OF(datasets),
	repet = ncol / vdim,
	nrow = NROW_OUT_OF(datasets),
	ndata = nrow * ncol,
	totptsvdim = nrow * vdim;
      
      L->X[GLOBAL.general.set] = 
	(double*) CALLOC( (betatot + repet * vdim) * Gettotalpoints(cov), 
                           sizeof(double)); // auch die Y-Werte; vdim included!
      
      
      if (!L->dettrend_has_nas && L->nas_boxcox == 0) {
	double *Xdata = L->X[GLOBAL.general.set] + betatot * totptsvdim;
	MEMCOPY(Xdata, data, ndata * sizeof(double));

	//  printf("data=%10g %d %d %d\n", data[0], ndata, betatot, totptsvdim);	
	//  BUG;
 
	for (p = m = 0; m < repet; m++) {
	  for (k = 0; k < totptsvdim; k++) {
	    Xdata[p++] -= YhatWithoutNA[k]; // XX lesefehler
	  }
	}
      
	if (L->nas_boxcox == 0 && R_FINITE(P0(GAUSS_BOXCOX)))
	  boxcox_trafo(P(GAUSS_BOXCOX), vdim, Xdata, nrow, repet);    
	
      }

      for (i=0; i<L->random; i++) {
	BUG; // not programmed yet
	if (L->nas_random[i] == 0) {
	  BUG;
	} else L->random_has_nas = true;
      }
    } // GLOBAL.set
  } // not rftrend

  for (GLOBAL.general.set = 0; GLOBAL.general.set < sets; GLOBAL.general.set++){
    for (i=0; i<L->fixedtrends; i++) {
      if (L->nas_fixed[i] == 0) {
	FctnIntern(cov, L->cov_fixed[i], L->cov_fixed[i], 
		   L->X[GLOBAL.general.set] + L->betas[i] * Gettotalpoints(cov),
		   false);// "Lueckentext"!!
	//	printf("hier %d %d incr=%d %10g\n", L->fixedtrends, L->nas_fixed[i],
	//	       L->betas[i] * Gettotalpoints(cov),
	//	       L->X[GLOBAL.general.set] [L->betas[i] * Gettotalpoints(cov)]);
	//	APMI(cov);
	
      }
      else L->fixedtrend_has_nas = true;
    }

    //    int kk = 0;
    //    for (int h=0; h<betatot; h++) { printf("\n");
    //      for (int hh=0; hh<Gettotalpoints(cov); hh++) {
    //	printf("%d ", (int)L->X[GLOBAL.general.set][kk++]);
    //      }
    //    } BUG;

  } // GLOBAL.set

  EXT_NEW_STORAGE(solve);
  

  //  assert(cov->Slikelihood != NULL);  printf("%d\n", addressbits(cov->Slikelihood));  APMI(cov);

 ErrorHandling:
  GLOBAL.general.set = store;

  RETURN_ERR(err);
}
  
    
void gauss_trend(model *predict, model *cov, double *v, int set) {
  // muss angenommen werden, dass alles gesetzt ist.
   likelihood_storage *L = cov->Slikelihood;			
  int store = GLOBAL.general.set;
  assert (set >= 0 && set <GET_LOC_SETS(predict));
  GLOBAL.general.set = set;
  listoftype *datasets = L->datasets;
  int 
    err = NOERROR,
    betatot = L->betas[L->fixedtrends],    
    vdim = VDIM0,
    ncol = NCOL_OUT_OF(datasets),
    repet = L->betas_separate ? ncol / vdim : 1,
    pred_tot = Gettotalpoints(predict),
    predtotvdim = pred_tot * vdim,
    predtotvdimrepet = pred_tot * ncol;
  double *X = NULL;
  
  for (int k=0; k<predtotvdimrepet; v[k++] = 0.0);
  if (L->ignore_trend) goto ErrorHandling; // exist
  

  if (!L->betas_separate && repet > 1) GERR("BUG");

  if ((X = (double*) MALLOC(predtotvdim * sizeof(double))) == NULL) {
    // printf("%d %d \n", pred_tot, vdim); APMI(cov);
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  
  int r,m;
  for (int i=0; i<L->dettrends; i++) {      
    FctnIntern(predict, L->cov_det[i],  L->cov_det[i], X, true);
    for (r = m = 0; r < repet; r++) {
      for (int k=0; k<predtotvdim; k++) v[m++] += X[k];
    }
  }
  for (int i=0; i<L->fixedtrends; i++) {      
    FctnIntern(predict, L->cov_fixed[i],  L->cov_fixed[i], X, true);
    if (L->betas[i+1] - L->betas[i] != 1) BUG;
    double *betavec = L->betavec + L->betas[i];
    for (r = m = 0; r < repet; r++) {
      double beta = *betavec;
      //	printf("beta = %10g %10g\n", beta, X[0]); BUG;
      for (int k=0; k<predtotvdim; k++) v[m++] += X[k] * beta;
      if (L->betas_separate) betavec += betatot;
    }
  }

 ErrorHandling :
  GLOBAL.general.set = store;
  FREE(X);
  if (err != NOERROR) XERR(err);

}
 
