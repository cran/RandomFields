/* 
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 library for simulation of random fields -- get key strukture

 Copyright (C) 2001 -- 2011 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
RO
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

//#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
//#include <sys/timeb.h>
 
#include <string.h>
#include "RF.h"
#include <Rdefines.h>
//#include "CovFcts.h"
//#include <unistd.h>
//#include <R_ext/Utils.h>     


SEXP Logi(bool* V, int n, int max, long *mem) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(VECSXP, 0);
  (*mem) += n * sizeof(bool);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(LGLSXP, n));
  for (i=0; i<n; i++) LOGICAL(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP Num(double* V, int n, int max, long *mem) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(REALSXP, 0);
  (*mem) += n * sizeof(double);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(REALSXP, n));
  for (i=0; i<n; i++) REAL(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP Result(res_type* V, int n, int max, long *mem) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(REALSXP, 0);
  (*mem) += n * sizeof(double);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(REALSXP, n));
  for (i=0; i<n; i++) REAL(dummy)[i] = (double) V[i];
  UNPROTECT(1);
  return dummy;
}


SEXP Int(int *V, int n, int max, long *mem) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(INTSXP, 0);
  (*mem) += n * sizeof(int);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(INTSXP, n));
  for (i=0; i<n; i++) INTEGER(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP Char(char *V, int n, int max, long *mem) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(INTSXP, 0);
  (*mem) += n * sizeof(char);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(INTSXP, n));
  for (i=0; i<n; i++) INTEGER(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

//#define MAX_INT 2147483647
#define MAX_INT   2000000000
SEXP Mat(double* V, int row, int col, int max, long *mem) {
  int i, n;
  if (V==NULL) return allocMatrix(REALSXP, 0, 0);
  n = row * col;
  (*mem) += n * sizeof(double);
  if (n>max) {
      int nn[2];
      nn[0] = row;
      nn[1] = col;
      return TooLarge(nn, 2);
  }
  SEXP dummy;
  PROTECT(dummy=allocMatrix(REALSXP, row, col));
  for (i=0; i<n; i++) REAL(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP ResultMat(res_type* V, int row, int col, int max, long *mem) {
  int i, n;
  if (V==NULL) return allocMatrix(REALSXP, 0, 0);
  n = row * col;
  (*mem) += n * sizeof(double);
  if (n>max) {
      int nn[2];
      nn[0] = row;
      nn[1] = col;
      return TooLarge(nn, 2);
  }
  SEXP dummy;
  PROTECT(dummy=allocMatrix(REALSXP, row, col));
  for (i=0; i<n; i++) REAL(dummy)[i] = (double) V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP MatInt(int* V, int row, int col, int max, long *mem) {
  int i, n;
  if (V==NULL) return allocMatrix(INTSXP, 0, 0);
  n = row * col;
  (*mem) += n * sizeof(int);
  if (n>max) {
      int nn[2];
      nn[0] = row;
      nn[1] = col;
      return TooLarge(nn, 2);
  }
  SEXP dummy;
  PROTECT(dummy=allocMatrix(INTSXP, row, col));
  for (i=0; i<n; i++) INTEGER(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP Array3D(double** V, int depth, int row, int col, int max, long *mem) {
  int i, j, m, n;
  if (V==NULL) return alloc3DArray(REALSXP, 0, 0, 0);
  m = row * col;
  n = row * col * depth;
  (*mem) += n * sizeof(double);
  if (n>max) {
      int nn[3];
      nn[0] = row;
      nn[1] = col;
      nn[2] = depth;
      return TooLarge(nn, 3);
  }
  SEXP dummy;
  PROTECT(dummy=alloc3DArray(REALSXP, depth, row, col));
  for(j=0; j<depth; j++) {
    for (i=0; i<m; i++) {
      REAL(dummy)[j*m+i] = V[j][i];
    }
  }
  UNPROTECT(1);
  return dummy;
}

SEXP TooLarge(int *n, int l){
#define nTooLarge 2 // mit op
  const char *tooLarge[nTooLarge] = {"size", "msg"};
  int i;
  long mem=0;
  SEXP namevec, info;
  PROTECT(info=allocVector(VECSXP, nTooLarge));
  PROTECT(namevec = allocVector(STRSXP, nTooLarge));
  for (i=0; i<nTooLarge; i++) SET_STRING_ELT(namevec, i, mkChar(tooLarge[i]));
  setAttrib(info, R_NamesSymbol, namevec);
  i=0;
  SET_VECTOR_ELT(info, i++, Int(n, l, l, &mem));
  SET_VECTOR_ELT(info, i,
		 mkString("too many elements - increase max.elements"));
  UNPROTECT(2);
  return info;
}


SEXP Param(void* p, int nrow, int ncol, SEXPTYPE type, bool drop, long* mem) {
  SEXP dummy = NULL;	
  int i;
  if (p == NULL) {
    dummy =  allocVector(REALSXP, 0); 
  } else if (type == REALSXP) {
      dummy = (ncol==1 && drop) ? Num((double*) p, nrow, MAX_INT, mem) 
	  : Mat((double*) p, nrow, ncol, MAX_INT, mem);
  } else if (type == INTSXP) {
      dummy = (ncol==1 && drop) ? Int((int*) p, nrow, MAX_INT, mem)
	  : MatInt((int*) p, nrow, ncol, MAX_INT, mem);
  } else if (type >= LISTOF){
      dummy = allocVector(VECSXP, nrow);
      listoftype *q = (listoftype *) p;
      for (i=0; i<nrow; i++) {
	  SET_VECTOR_ELT(dummy, i, 
			 Param(q->p[i], q->nrow[i], q->ncol[i], REALSXP, false,
			       mem));
      }
  } else {
    assert(false);
  }
  return dummy;
}


SEXP GetModelInfo(cov_model *cov, int level, bool gatter, long *mem)
{
#define ninfo0 3
#define ninfo1 4
#define ninfo2 4
#define ninfo3 11
     
//  printf("1st gmi %s\n", CovList[cov->nr].name);
  if (cov == NULL) return allocVector(VECSXP, 0);
  SEXP model, submodels, nameMvec, param, pnames;
  int i, j, nmodelinfo, k = 0;
  cov_fct *C = CovList + cov->nr;
  
//    printf("2nd gmi\n");
  if (!gatter && cov->nr >= GATTER && cov->nr <= LASTGATTER) 
      return GetModelInfo(cov->sub[0], level, gatter, mem);
//   printf("3rd gmi\n");
  
  nmodelinfo = ninfo0;
  switch (level) {
      case 3 : nmodelinfo += ninfo3;
      case 2 : nmodelinfo += ninfo2;
      case 1 : nmodelinfo += ninfo1;
   }
  if (cov->nsub==0) nmodelinfo--;

  PROTECT(model = allocVector(VECSXP, nmodelinfo));
  PROTECT(nameMvec = allocVector(STRSXP, nmodelinfo));

  SET_STRING_ELT(nameMvec, k, mkChar("name"));
  cov_fct *CC = CovList + cov->nr;
  while(strncmp(CC->name, InternalName, strlen(InternalName)) ==0) CC--;
  SET_VECTOR_ELT(model, k++, mkString(CC->name));

  SET_STRING_ELT(nameMvec, k, mkChar("param"));  
  int notnull = 0;
  for (i=0; i<C->kappas; i++) {
      if (cov->nrow[i]>0 && cov->ncol[i]>0) notnull++;
  }
  
  PROTECT(param = allocVector(VECSXP, notnull));
  PROTECT(pnames = allocVector(STRSXP, notnull));
 
  for (j=i=0; i<C->kappas; i++) {
    if (cov->nrow[i]>0 && cov->ncol[i]>0) { 
      SET_STRING_ELT(pnames, j, mkChar(C->kappanames[i]));
      SET_VECTOR_ELT(param, j,
		     Param((void*) cov->p[i], cov->nrow[i], cov->ncol[i], 
			   C->kappatype[i], true, mem));    
    j++;
    }
  }
  setAttrib(param, R_NamesSymbol, pnames);
//  printf("gmi !\n");

  SET_VECTOR_ELT(model, k++, param);
  UNPROTECT(2);

//  goto END;

//  printf("start GMI 1 %d\n", level);
  
  if (level>=1) {      
      SET_STRING_ELT(nameMvec, k, mkChar("covnr"));
      SET_VECTOR_ELT(model, k++, ScalarInteger(cov->nr));

      SET_STRING_ELT(nameMvec, k, mkChar("vdim"));
      SET_VECTOR_ELT(model, k++, ScalarInteger(cov->vdim));

//  SET_STRING_ELT(nameMvec, k, mkChar("naturalscaling"));  
//  SET_VECTOR_ELT(model, k++, ScalarInteger(cov->naturalscaling));

      SET_STRING_ELT(nameMvec, k, mkChar("tsdim"));  
      SET_VECTOR_ELT(model, k++, ScalarInteger(cov->tsdim));
      
      SET_STRING_ELT(nameMvec, k, mkChar("xdim"));  
      SET_VECTOR_ELT(model, k++, ScalarInteger(cov->xdim));

  } 
   
//  printf("GMI 2\n");

  if (level>=2) {
      SET_STRING_ELT(nameMvec, k, mkChar("statIn"));  
      SET_VECTOR_ELT(model, k++, ScalarInteger(cov->statIn));
      
      SET_STRING_ELT(nameMvec, k, mkChar("isoIn"));
      SET_VECTOR_ELT(model, k++, ScalarInteger(cov->isoIn));

      SET_STRING_ELT(nameMvec, k, mkChar("internalq"));  
      SET_VECTOR_ELT(model, k++, Num(cov->q, cov->qlen, MAX_INT, mem));

      SET_STRING_ELT(nameMvec, k, mkChar("pref"));  
      SET_VECTOR_ELT(model, k++, Int(cov->pref, Nothing + 1, MAXINT, mem));
  } 

//  printf("GMI 3\n");
  if (level>=3) {
      SET_STRING_ELT(nameMvec, k, mkChar("maxdim"));  
      SET_VECTOR_ELT(model, k++, ScalarInteger(cov->maxdim));
      
      SET_STRING_ELT(nameMvec, k, mkChar("derivatives"));  
      SET_VECTOR_ELT(model, k++, ScalarInteger(cov->derivatives));
      
      SET_STRING_ELT(nameMvec, k, mkChar("normalmix"));  
      SET_VECTOR_ELT(model, k++, ScalarLogical(cov->normalmix));
      
      SET_STRING_ELT(nameMvec, k, mkChar("anyNAdown"));  
      SET_VECTOR_ELT(model, k++, ScalarInteger((int) cov->anyNAdown));
      
      SET_STRING_ELT(nameMvec, k, mkChar("anyNAscaleup"));  
      SET_VECTOR_ELT(model, k++, ScalarInteger((int) cov->anyNAscaleup));
      
      SET_STRING_ELT(nameMvec, k, mkChar("MLE"));  
      SET_VECTOR_ELT(model, k++, ScalarLogical(cov->MLE != NULL));

      SET_STRING_ELT(nameMvec, k, mkChar("finiterange"));  
      SET_VECTOR_ELT(model, k++, ScalarLogical(cov->finiterange));
      
      SET_STRING_ELT(nameMvec, k, mkChar("diag"));  
      SET_VECTOR_ELT(model, k++, ScalarLogical(cov->diag));
  
//  SET_STRING_ELT(nameMvec, k, mkChar("quasidiag"));  
//  SET_VECTOR_ELT(model, k++, ScalarLogical(cov->quasidiag));

      SET_STRING_ELT(nameMvec, k, mkChar("semisep.last"));  
      SET_VECTOR_ELT(model, k++, ScalarLogical(cov->semiseparatelast));
      
      SET_STRING_ELT(nameMvec, k, mkChar("sep.last"));  
      SET_VECTOR_ELT(model, k++, ScalarLogical(cov->separatelast));
      
//  SET_STRING_ELT(nameMvec, k, mkChar("idx"));  
//  SET_VECTOR_ELT(model, k++, Int(cov->idx, cov->tsdim, MAXINT, mem));
      SET_STRING_ELT(nameMvec, k, mkChar("user"));  
      SET_VECTOR_ELT(model, k++, Int(cov->user, Nothing + 1, MAXINT, mem));
  } 
  		 
   if (cov->nsub > 0) {
    SET_STRING_ELT(nameMvec, k, mkChar("submodels"));  
    //   printf("gmi 1\n");
    PROTECT(submodels = allocVector(VECSXP, cov->nsub));
    int zaehler = 0;
//    printf("gmi 2\n");
    for (i=0; i<MAXSUB; i++) {
//	printf("gmi x %d %d\n", i, MAXSUB);
      // printf("\n\n %d %i %s\n", zaehler, i, CovList[cov->nr].name);
      if (cov->sub[i] != NULL) {
//	printf("gmi a\n");
	if (level >= 2 || cov->sub[i]->nr < GATTER || 
	    cov->sub[i]->nr > LASTGATTER) {
//	printf("gmi b\n");
//	PrintModelInfo(cov);
//	printf("gmi XX %d %d %d %d\n", i, MAXSUB, level, gatter);
	    SET_VECTOR_ELT(submodels, zaehler, 
			   GetModelInfo(cov->sub[i], level, gatter, mem));
	}
	else {
//	printf("gmi c\n");
	  SET_VECTOR_ELT(submodels, zaehler, 
			 GetModelInfo(cov->sub[i]->sub[0], level, gatter,mem));
	}
	if (++zaehler >= cov->nsub) break;
      }
    }
    SET_VECTOR_ELT(model, k++, submodels);
    UNPROTECT(1);
 
//    printf("return GMI\n");
  }

//  SET_STRING_ELT(nameMvec, k, mkChar("method"));  
//  SET_VECTOR_ELT(model, k++, mkChar(METHODNAMES[(int) cov->usermethod]));
  
  setAttrib(model, R_NamesSymbol, nameMvec);

// END:

  assert(k==nmodelinfo); 
  UNPROTECT(2); // model + namemodelvec
  return model;
}

SEXP GetExtModelInfo(SEXP keynr, SEXP Level, SEXP gatter) {
  int knr, level;
  long mem;
  knr = INTEGER(keynr)[0];
  level = INTEGER(Level)[0];
  if (knr>=0) {
    if (knr < MAXKEYS) {
      key_type *key;
      key = &(KEY[knr]);
      return key->cov == NULL ? allocVector(VECSXP, 0) 
	  : GetModelInfo(key->cov, level, (bool) INTEGER(gatter)[0], &mem);
    }
  } else {
    knr = -knr-1;
    if (knr < MODEL_MAX && STORED_MODEL[knr] != NULL)
	return GetModelInfo(STORED_MODEL[knr], level,(bool) INTEGER(gatter)[0],
			    &mem); // 0
  } 
  return allocVector(VECSXP, 0);
}


#define nglobalinfo 0
SEXP GetGlobalInfo(globalparam global, long *mem) {
  SEXP namevec, l; 

  PROTECT(l = allocVector(VECSXP, nglobalinfo));
  PROTECT(namevec = allocVector(STRSXP, nglobalinfo));
  setAttrib(l, R_NamesSymbol, namevec);
  assert(0 == nglobalinfo);
  UNPROTECT(2);
  return l;
}

#define nlocinfo 10
SEXP GetLocationInfo(location_type *loc, long *mem) {
  if (loc == NULL) return allocVector(VECSXP, 0);
  const char *info[nlocinfo] = 
    {"timespacedim", "length", "spatialdim", "spatialtotpts", "totpts", 
     "grid", "Time", "xgr", "x", "T"};
  SEXP namevec, l;
  int k,
    tsdim = loc->timespacedim,
    spdim = loc->spatialdim;

//  printf("%d %d\n", tsdim, spdim);

  PROTECT(l = allocVector(VECSXP, nlocinfo));
  PROTECT(namevec = allocVector(STRSXP, nlocinfo));
  for (k=0; k<nlocinfo; k++)
       SET_STRING_ELT(namevec, k, mkChar(info[k]));

  k = 0;
  SET_VECTOR_ELT(l, k++, ScalarInteger(tsdim));
  SET_VECTOR_ELT(l, k++, Int(loc->length, tsdim, MAX_INT, mem));
  SET_VECTOR_ELT(l, k++, ScalarInteger(loc->spatialdim));
  SET_VECTOR_ELT(l, k++, ScalarInteger((int) loc->spatialtotalpoints));
  SET_VECTOR_ELT(l, k++, ScalarInteger((int) loc->totalpoints));
  SET_VECTOR_ELT(l, k++, ScalarLogical(loc->grid));
  SET_VECTOR_ELT(l, k++, ScalarLogical(loc->Time));
  SET_VECTOR_ELT(l, k++, 
		 Mat(loc->xgr[0], loc->grid ? 3 : 0, 
		     spdim, MAX_INT, mem));
  SET_VECTOR_ELT(l, k++, 
		 Mat(loc->x, loc->spatialdim, 
		     loc->grid ? 0 : loc->spatialtotalpoints, 
		     MAX_INT, mem));
  SET_VECTOR_ELT(l, k++, Num(loc->T, loc->Time ? 3 : 0, MAX_INT, mem));
  
  setAttrib(l, R_NamesSymbol, namevec);
  UNPROTECT(2); // l + namelvec
  assert(k == nlocinfo);
  return l;
}



#define ntrendinfo 5
SEXP GetTrendInfo(trend_type *trend, long *mem) {
  if (trend == NULL) return allocVector(VECSXP, 0);
  const char *info[ntrendinfo] = 
    {"lTrendFct", "TrendModus", "TrendFct", "mean", "LinTrend"};
  SEXP namevec, l;
  int k;


  PROTECT(l = allocVector(VECSXP, ntrendinfo));
  PROTECT(namevec = allocVector(STRSXP, ntrendinfo));
  for (k=0; k<ntrendinfo; k++) {
    SET_STRING_ELT(namevec, k, mkChar(info[k]));
  }

  k = 0;
  SET_VECTOR_ELT(l, k++, ScalarInteger(trend->lTrendFct));
  SET_VECTOR_ELT(l, k++, ScalarInteger(trend->TrendModus));
  SET_VECTOR_ELT(l, k++, mkString(trend->TrendFunction == NULL ? "" : 
				  trend->TrendFunction));
  SET_VECTOR_ELT(l, k++, ScalarReal(trend->mean));
  SET_VECTOR_ELT(l, k++, Num(trend->LinearTrend, trend->lLinTrend, MAX_INT,
			     mem));
  setAttrib(l, R_NamesSymbol, namevec);
  UNPROTECT(2); // l + namelvec
  assert(k == ntrendinfo);

  return l;
}
  

#define nsimuinfo 4
SEXP GetSimuInfo(simu_type *simu, long *mem) {
  if (simu == NULL) return allocVector(VECSXP, 0);
  const char *info[nsimuinfo] = 
    {"active", "stop", "distr", "expect.simu"};
  SEXP namevec, l;
  int k;

  PROTECT(l = allocVector(VECSXP, nsimuinfo));
  PROTECT(namevec = allocVector(STRSXP, nsimuinfo));
  for (k=0; k<nsimuinfo; k++)
       SET_STRING_ELT(namevec, k, mkChar(info[k]));

  k = 0;
  SET_VECTOR_ELT(l, k++, ScalarLogical(simu->active));
  SET_VECTOR_ELT(l, k++, ScalarLogical(simu->stop));
  SET_VECTOR_ELT(l, k++, mkString(DISTRNAMES[simu->distribution]));
  SET_VECTOR_ELT(l, k++, ScalarInteger(simu->expected_number_simu));
  assert(k==nsimuinfo); 
  
  setAttrib(l, R_NamesSymbol, namevec);
  UNPROTECT(2); // l + namelvec
  assert(k == nsimuinfo);

  return l;
}


#define nmethodinfo 14
SEXP GetMethodInfo(method_type *meth, 
		   bool ignore_active, int depth, int max, long *mem) {
  if (meth==NULL) {
    return allocVector(VECSXP, 0);
  }

  const char *methodinfo[nmethodinfo] = 
    {"name", "compatible", "cov", 
     "caniso","cproj",  "cscale", "cvar", "matrixtype", "hanging",
      "space", "sptime", 
     "S", "sub", "mem"};

 SEXP namemethodvec, method, nameSvec, S, submethods;
  location_type *loc = meth->loc;
  int i, k,l, dim,
      level = 3,
      NR = -meth->nr -1;
  long totalpoints, dummymem, timespacedim, Time, xdim; 

  if (meth->loc==NULL) {
    timespacedim = totalpoints = Time = 0;  
  } else {
    totalpoints = loc->totalpoints;
    timespacedim = loc->timespacedim;
    Time = loc->Time;
  }
  if (meth->cov==NULL) {
    xdim = dim = 0;
  } else {
    dim = meth->cov->tsdim;
    xdim = meth->cov->xdim;
  }


  const char *invm[NoFurtherInversionMethod + 1] =
    {"Cholesky", "SVD", "None"};

  PROTECT(method = allocVector(VECSXP, nmethodinfo));
  PROTECT(namemethodvec = allocVector(STRSXP, nmethodinfo));
  for (k=0; k<nmethodinfo; k++)
       SET_STRING_ELT(namemethodvec, k, mkChar(methodinfo[k]));

  l = 0;
  if (NR >= 0) {
    SET_VECTOR_ELT(method, l++, mkString(METHODNAMES[NR]));
  } else {
    cov_fct *CC = CovList + meth->nr;
    while(strcmp(CC->name, InternalName) ==0) {
	// printf("%s %s\n", CC->name, InternalName);
      CC--;
    }
    SET_VECTOR_ELT(method, l++, mkString(CC->name));
  }

  //printf("%d %d %d %d\n", MaxMpp, ExtremalGauss, meth->nr, NR); assert(false);


  SET_VECTOR_ELT(method, l++, 
		 ScalarLogical(meth->compatible));
  SET_VECTOR_ELT(method, l++, GetModelInfo(meth->cov, level, true, &dummymem));
  SET_VECTOR_ELT(method, l++, Mat(meth->caniso, timespacedim,
				  xdim, MAX_INT, mem));
  SET_VECTOR_ELT(method, l++, Int(meth->cproj, xdim, MAX_INT, mem));
  SET_VECTOR_ELT(method, l++, ScalarReal(meth->cscale));
  SET_VECTOR_ELT(method, l++, ScalarReal(meth->cvar));
  SET_VECTOR_ELT(method, l++, ScalarInteger(meth->type));
  SET_VECTOR_ELT(method, l++, ScalarLogical(meth->hanging != NULL));
  SET_VECTOR_ELT(method, l++, Num(meth->space,
				  (dim - Time) * totalpoints, 
				  MAX_INT, mem));
  SET_VECTOR_ELT(method, l++, Num(meth->sptime, dim * totalpoints, 
				  MAX_INT, mem));

 //  SET_VECTOR_ELT(method, l++, Mat(, dim, MAX_INT, mem));
  
  int nS, nSlist[Forbidden + 1] =
    {12 /* CE */, 2 /*Cutoff*/, 2 /* Intr */, 
     4 /* TBM2 */, 4 /* TBM3 */,  0 /*Spectral */,
     3 /* dir */, 12 /* sequ */,  0 /* Markov */, 
     0 /* average */,
     3 /* nug */, 24 /* coin */,  4 /* hyp */, 
     1 /* noth */, 11 /* maxmpp */,
     4 /*extremalGauss */, 1 /* Forbidden */};
  assert(Forbidden == 16);

  k = 0;
  
//  printf("nS : %d %d %d\n", meth->S!=NULL, depth, NR);

  nS = (meth->S!=NULL && depth<=4 && NR >= 0) ? nSlist[NR] : 0;
 
  PROTECT(S = allocVector(VECSXP, nS));
  PROTECT(nameSvec = allocVector(STRSXP, nS));
  if (nS != 0) switch (NR) {
      case CircEmbed : {
	CE_storage* s;
	s = (CE_storage*) meth->S;
	int vdim = s->vdim, 
	  vdimSQ = vdim * vdim;
	SET_STRING_ELT(nameSvec, k, mkChar("size"));
	SET_VECTOR_ELT(S, k++, Int(s->m, dim, MAX_INT, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("trials"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->trials));
	SET_STRING_ELT(nameSvec, k, mkChar("next_new"));
	SET_VECTOR_ELT(S, k++, ScalarLogical(s->new_simulation_next));
	SET_STRING_ELT(nameSvec, k, mkChar("curSimuPosition"));
	SET_VECTOR_ELT(S, k++, Int(s->cur_square, dim, MAX_INT, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("simupositions"));
	SET_VECTOR_ELT(S, k++, Int(s->max_squares, dim, MAX_INT, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("segmentLength"));
	{
	  int dummy[MAXCEDIM], i;
	  for (i=0; i< dim; i++) 
	    dummy[i] = s-> square_seg[i] / s->cumm[i];
	  SET_VECTOR_ELT(S, k++, Int(dummy, dim, MAX_INT, mem));
	}
	SET_STRING_ELT(nameSvec, k, mkChar("totalsize"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->mtot));
	SET_STRING_ELT(nameSvec, k, mkChar("smallest.neg.real"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->smallestRe));
	SET_STRING_ELT(nameSvec, k, mkChar("largest.abs.imag"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->largestAbsIm));
	SET_STRING_ELT(nameSvec, k, mkChar("fft"));
	SET_VECTOR_ELT(S, k++,
		       vdim == 1 
		       ? Mat(s->c[0], 2, s->mtot, max, mem)
	               : Array3D(s->c, vdimSQ, 2, s->mtot, max, mem));	    
	SET_STRING_ELT(nameSvec, k, mkChar("invfft"));
	SET_VECTOR_ELT(S, k++, 
		       (vdim == 1 || true) ? ScalarInteger(s->mtot) : true
		       ? Mat(s->d[0], 2, s->mtot, max, mem)
	               : Array3D(s->d, vdimSQ, 2, s->mtot, max, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("positivedefinite"));
	SET_VECTOR_ELT(S, k++, ScalarLogical(s->positivedefinite));
      } break;
      case CircEmbedCutoff : case CircEmbedIntrinsic : {
	localCE_storage* s;
	s = (localCE_storage*) meth->S;
	SET_STRING_ELT(nameSvec, k, mkChar("correctionTerms"));
	SET_VECTOR_ELT(S, k++, Mat((double*) s->correction, dim, dim,
				   MAX_INT, mem));

//	printf("calling form meth\n");	

	SET_STRING_ELT(nameSvec, k, mkChar("new"));
	SET_VECTOR_ELT(S, k++, InternalGetKeyInfo(&(s->key), ignore_active,
						      depth + 1, max));
//	printf("back to meth\n");	
      } break;
      case TBM2: case TBM3 : {
	TBM_storage* s;
	s = (TBM_storage*) meth->S;
	int* length = s->key.loc.length;
//	    SET_STRING_ELT(nameSvec, k, mkChar("aniso"));
//	    SET_VECTOR_ELT(S, k++, 
//			   Mat(s->aniso, dim, s->reduceddim, MAX_INT, mem));
//	    SET_STRING_ELT(nameSvec, k, mkChar("simugrid"));
//	    SET_VECTOR_ELT(S, k++, ScalarLogical(s->simugrid));	
//	    SET_STRING_ELT(nameSvec, k, mkChar("simuspatialdim"));
//	    SET_VECTOR_ELT(S, k++, ScalarInteger(s->simuspatialdim));	
	SET_STRING_ELT(nameSvec, k, mkChar("ce_dim"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->ce_dim));
//	    SET_STRING_ELT(nameSvec, k, mkChar("reduceddim"));
//	    SET_VECTOR_ELT(S, k++, ScalarInteger(s->reduceddim));		
	SET_STRING_ELT(nameSvec, k, mkChar("center"));
	SET_VECTOR_ELT(S, k++, Num(s->center, dim - 
				   (s->ce_dim==2), MAX_INT, mem));
//	SET_STRING_ELT(nameSvec, k, mkChar("x"));
//	SET_VECTOR_ELT(S, k++, Mat(s->x, dim, 
//				   loc->grid ? 3 : totalpoints, max, mem));
//	SET_STRING_ELT(nameSvec, k, mkChar("xsimgr"));
//	SET_VECTOR_ELT(S, k++, Mat(s->xsimugr[0], loc->grid ? 3 : 0,
//				   dim, MAX_INT, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("line")); // never change this
	//                                 name; otherwise the tricky
	//                                 TBM examples does not work anymore
	SET_VECTOR_ELT(S, k++, (s->ce_dim <= 1)
		       ? Result(s->simuline, length[0], max, mem)
		       : ResultMat(s->simuline, length[0], length[1], max, mem)
	  );
	SET_STRING_ELT(nameSvec, k, mkChar("new"));
	SET_VECTOR_ELT(S, k++, InternalGetKeyInfo(&(s->key), ignore_active,
						  depth + 1, max));
      } break;
      case SpectralTBM: {
	/*
	spectral_storage* s;
	s = (spectral_storage*) meth->S;
	SET_STRING_ELT(nameSvec, k, mkChar("density"));
	SET_VECTOR_ELT(S, k++, ScalarLogical(s->density != NULL));    
	SET_STRING_ELT(nameSvec, k, mkChar("sigmametro"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->sigma));    
	SET_STRING_ELT(nameSvec, k, mkChar("nmetro"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->nmetro));    
	SET_STRING_ELT(nameSvec, k, mkChar("phistep2d"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->phistep2d));    
	SET_STRING_ELT(nameSvec, k, mkChar("phi2d"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->phi2d));    
	SET_STRING_ELT(nameSvec, k, mkChar("grid"));
	SET_VECTOR_ELT(S, k++, ScalarLogical(s->grid));    
	*/	
      } break;
      case Direct : {
	direct_storage* s;
	s = (direct_storage*) meth->S;
	SET_STRING_ELT(nameSvec, k, mkChar("invmethod"));
	SET_VECTOR_ELT(S, k++, (s->method < NoFurtherInversionMethod) ?
		       mkString(invm[s->method]) :  ScalarLogical(false));
	SET_STRING_ELT(nameSvec, k, mkChar("sqrtCov"));
	SET_VECTOR_ELT(S, k++, Mat(s->U, totalpoints, totalpoints, max, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("dummy"));
	SET_VECTOR_ELT(S, k++, Num(s->G, totalpoints + 1, max, mem));	    
      }  break;
      case Sequential : {
	sequential_storage *s;
	s = (sequential_storage*) meth->S;
	SET_STRING_ELT(nameSvec, k, mkChar("back"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->back));
	SET_STRING_ELT(nameSvec, k, mkChar("totpnts"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->totpnts));
	SET_STRING_ELT(nameSvec, k, mkChar("spatialpnts"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->spatialpnts));
	SET_STRING_ELT(nameSvec, k, mkChar("ntime"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->ntime));
	SET_STRING_ELT(nameSvec, k, mkChar("initial"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->initial));
	SET_STRING_ELT(nameSvec, k, mkChar("sqrtCov22"));
	SET_VECTOR_ELT(S, k++, Mat(s->U22, s->totpnts, s->totpnts, 
				   max, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("sqrtCov11"));
	SET_VECTOR_ELT(S, k++, Mat(s->U11, s->spatialpnts, s->spatialpnts, 
				   max, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("Cov21.debug"));
	SET_VECTOR_ELT(S, k++, Mat(s->Cov21, s->totpnts, s->spatialpnts,
				   max, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("Inv22.debug"));
	SET_VECTOR_ELT(S, k++, Mat(s->Inv22, s->totpnts, s->totpnts,
				   max, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("mutransposed"));
	SET_VECTOR_ELT(S, k++, Mat(s->MuT, s->totpnts, s->spatialpnts, 
				   max, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("dummy"));
	SET_VECTOR_ELT(S, k++, Num(s->G, s->totpnts, max, mem)); 
	SET_STRING_ELT(nameSvec, k, mkChar("res0"));
	SET_VECTOR_ELT(S, k++, Result(s->res0, s->totpnts, max, mem)); 
      }
      case Markov : {
	/*
	  markov_storage* s;
	  s = (markov_storage*) meth->S;
	*/
      } break;
      case Average : {
      } break;
      case Nugget : {
	nugget_storage* s;
	s = (nugget_storage*) meth->S;
	SET_STRING_ELT(nameSvec, k, mkChar("simple"));
	SET_VECTOR_ELT(S, k++, ScalarLogical(s->simple));	
	SET_STRING_ELT(nameSvec, k, mkChar("simugrid"));
	SET_VECTOR_ELT(S, k++, ScalarLogical(s->simugrid));	
//	SET_STRING_ELT(nameSvec, k, mkChar("srqtnugget"));
//	SET_VECTOR_ELT(S, k++, ScalarReal(s->sqrtnugget));	
	SET_STRING_ELT(nameSvec, k, mkChar("internalsort"));
	SET_VECTOR_ELT(S, k++, Int(s->pos, totalpoints, max, mem));	
      }  break;
      case RandomCoin : {
	mpp_storage* s;
	s = (mpp_storage*) meth->S;

	SET_STRING_ELT(nameSvec, k, mkChar("integral"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->integral));
	SET_STRING_ELT(nameSvec, k, mkChar("integralsq"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->integralsq));
	SET_STRING_ELT(nameSvec, k, mkChar("effectiveRadius"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->effectiveRadius));
	SET_STRING_ELT(nameSvec, k, mkChar("effectivearea"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->effectivearea));
	SET_STRING_ELT(nameSvec, k, mkChar("plus"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->plus));
	SET_STRING_ELT(nameSvec, k, mkChar("relplus"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->relplus));	    
	SET_STRING_ELT(nameSvec, k, mkChar("lensimu"));
	SET_VECTOR_ELT(S, k++, Num(s->lensimu, dim, MAX_INT, mem)); 
	SET_STRING_ELT(nameSvec, k, mkChar("min"));
	SET_VECTOR_ELT(S, k++, Num(s->min, dim, MAX_INT, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("maxgrid"));
	SET_VECTOR_ELT(S, k++, Num(s->max, dim, MAX_INT, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("minsimu"));
	SET_VECTOR_ELT(S, k++, Num(s->minsimu, dim, MAX_INT, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("maxsimu"));
	SET_VECTOR_ELT(S, k++, Num(s->maxsimu, dim, MAX_INT, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("mean"));
	SET_VECTOR_ELT(S, k++, Num(s->mean, dim, MAX_INT, mem));
	SET_STRING_ELT(nameSvec, k, mkChar("sdgauss"));
	SET_VECTOR_ELT(S, k++, Num(s->sdgauss, dim, MAX_INT, mem));

	SET_STRING_ELT(nameSvec, k, mkChar("internal.constants"));
	SET_VECTOR_ELT(S, k++, Num(s->c, 6, MAX_INT, mem)); 
	SET_STRING_ELT(nameSvec, k, mkChar("u"));
	SET_VECTOR_ELT(S, k++, Num(s->u, dim, MAX_INT, mem));
//	SET_STRING_ELT(nameSvec, k, mkChar("var"));
//	SET_VECTOR_ELT(S, k++, ScalarReal(s->var));
	SET_STRING_ELT(nameSvec, k, mkChar("logapproxzero"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->logapproxzero));
	SET_STRING_ELT(nameSvec, k, mkChar("samplingdist"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->samplingdist));
	SET_STRING_ELT(nameSvec, k, mkChar("samplingr"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->samplingr));
	SET_STRING_ELT(nameSvec, k, mkChar("average"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->average));
	SET_STRING_ELT(nameSvec, k, mkChar("factor"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->factor));

	SET_STRING_ELT(nameSvec, k, mkChar("dim"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->dim));
	SET_STRING_ELT(nameSvec, k, mkChar("ntot"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->ntot));
	SET_STRING_ELT(nameSvec, k, mkChar("invscale"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->factor));
	SET_STRING_ELT(nameSvec, k, mkChar("logInvSqrtDens"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->logInvSqrtDens));
      }  break;
      case Hyperplane : {
	hyper_storage* s;
	s = (hyper_storage*) meth->S;
	SET_STRING_ELT(nameSvec, k, mkChar("rx"));
	SET_VECTOR_ELT(S, k++, Num(s->rx, dim, MAX_INT, mem)); 
	SET_STRING_ELT(nameSvec, k, mkChar("center"));
	SET_VECTOR_ELT(S, k++, Num(s->center, dim, MAX_INT, mem)); 
	SET_STRING_ELT(nameSvec, k, mkChar("radius"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->radius));	
	SET_STRING_ELT(nameSvec, k, mkChar("HyperplaneFctOK"));
	SET_VECTOR_ELT(S, k++, ScalarLogical(s->hyperplane != NULL));
      } break;
      case MaxMpp : {
	mpp_storage* s;
	s = (mpp_storage*) meth->S;
	SET_STRING_ELT(nameSvec, k, mkChar("integralpos")); 
	SET_VECTOR_ELT(S, k++, ScalarReal(s->integralpos));
	SET_STRING_ELT(nameSvec, k, mkChar("factor"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->factor));
	SET_STRING_ELT(nameSvec, k, mkChar("maxheight"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->maxheight));	    
	SET_STRING_ELT(nameSvec, k, mkChar("effectiveRadius"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->effectiveRadius));
	SET_STRING_ELT(nameSvec, k, mkChar("effectivearea"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->effectivearea));
	SET_STRING_ELT(nameSvec, k, mkChar("plus"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->plus));
	SET_STRING_ELT(nameSvec, k, mkChar("min"));
	SET_VECTOR_ELT(S, k++, Num(s->min, dim, MAX_INT, mem));	    
	SET_STRING_ELT(nameSvec, k, mkChar("lensimu"));
	SET_VECTOR_ELT(S, k++, Num(s->lensimu, dim, MAX_INT, mem)); 
	SET_STRING_ELT(nameSvec, k, mkChar("internal.constants"));
	SET_VECTOR_ELT(S, k++, Num(s->c, 6, MAX_INT, mem)); 
//	SET_STRING_ELT(nameSvec, k, mkChar("primitiveFctOK"));
//	SET_VECTOR_ELT(S, k++, ScalarLogical(s->MppFct != NULL));
	SET_STRING_ELT(nameSvec, k, mkChar("dim"));
	SET_VECTOR_ELT(S, k++, ScalarInteger(s->dim));
      } break;
      case ExtremalGauss : {
	extremes_storage* s;
	s = (extremes_storage*) meth->S;
	SET_STRING_ELT(nameSvec, k, mkChar("inv_mean_pos"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->inv_mean_pos));
	SET_STRING_ELT(nameSvec, k, mkChar("assumedmax"));
	SET_VECTOR_ELT(S, k++, ScalarReal(s->assumedmax));
	SET_STRING_ELT(nameSvec, k, mkChar("field")); 
	
	SET_VECTOR_ELT(S, k++, Result(s->rf, 0, max, mem));	
	
	/*
	  locdim = loc->grid ? dim : 1;
	  if (locdim==1) {
	    SET_VECTOR_ELT(S, k++, Num(s->rf, totalpoints, max, mem));	
	  } else {
	    PROTECT(dummy=allocVector(INTSXP, locdim));
	    for (i=0; i<locdim; i++) INTEGER(dummy)[i] = loc->length[i];
	    PROTECT(array=allocArray(REALSXP, dummy));
	    for (i=0; i<totalpoints; i++) REAL(array)[i] = s->rf[i];
	    SET_VECTOR_ELT(S, k++, array);
	    UNPROTECT(2);
	  }
	*/

	SET_STRING_ELT(nameSvec, k, mkChar("new"));
	SET_VECTOR_ELT(S, k++, InternalGetKeyInfo(&(s->key), ignore_active,
						  depth + 1, max));
      } break;
      default :  
	assert(meth->nr <0);
	SET_STRING_ELT(nameSvec, k, mkChar("calling.method from "));
	SET_VECTOR_ELT(S, k++, mkString(CovList[- meth->nr].name));   
  }

  assert(k==nS);
  setAttrib(S, R_NamesSymbol, nameSvec);
  SET_VECTOR_ELT(method, l++, S);
 
  // SET_STRING_ELT(nameSvec, , mkChar("submethods"));  
  PROTECT(submethods = allocVector(VECSXP, meth->nsub));
  for (i=0; i<meth->nsub; i++) {
    SET_VECTOR_ELT(submethods, i, 
		   GetMethodInfo(meth->sub[i], ignore_active, depth + 1, max,
				 mem));
  }
  SET_VECTOR_ELT(method, l++, submethods);
  UNPROTECT(1); 

  //ERROR: cast from ‘long int*’ to ‘int’ loses precision
  SET_VECTOR_ELT(method, l++, ScalarInteger((int) *mem)); //original: (int) mem

  setAttrib(method, R_NamesSymbol, namemethodvec);
  UNPROTECT(2); // method + namemethodvec

  assert(l == nmethodinfo);
  
  UNPROTECT(2); // method + namemethodvec
  return method;
}



#define ninfo 8
SEXP InternalGetKeyInfo(key_type *key, bool ignore_active, int depth,
		        int max){
  const char *infonames[ninfo] = 
    {"gp", "gpdo", "simu", "loc", "trend", "cov", "meth", "mem"
    };
  int ni=0, actninfo, i,
      level = 3;
  long mem=0;
  SEXP info, namevec; 

  actninfo = (key->simu.active || ignore_active) ? ninfo : 3;
  PROTECT(info=allocVector(VECSXP, actninfo));
  PROTECT(namevec = allocVector(STRSXP, actninfo));
  for (i=0; i<actninfo; i++) SET_STRING_ELT(namevec, i, mkChar(infonames[i]));
  setAttrib(info, R_NamesSymbol, namevec);
  
  //printf("gp\n");
  SET_VECTOR_ELT(info, ni++, GetGlobalInfo(key->gp, &mem));
  //printf("gpdo\n");
  SET_VECTOR_ELT(info, ni++, GetGlobalInfo(key->gpdo, &mem));
  //printf("simu\n");
  SET_VECTOR_ELT(info, ni++, GetSimuInfo(&(key->simu), &mem));

  if (actninfo > 3) {
//  printf("loc\n");
    SET_VECTOR_ELT(info, ni++, GetLocationInfo(&(key->loc), &mem));
//    printf("trnde\n");
    SET_VECTOR_ELT(info, ni++, GetTrendInfo(&(key->trend), &mem));
//    printf("cov\n");
    SET_VECTOR_ELT(info, ni++, GetModelInfo(key->cov, level, true, &mem));
//   printf("meth\n");
    
     //   printf("%d\n", key->meth);
    //   printf("%d\n",ignore_active );
    //  printf("%d\n", depth);
    //  printf("%d\n", max);
    //  printf("%d\n", mem);

    SET_VECTOR_ELT(info, ni++, GetMethodInfo(key->meth, ignore_active, 
					     depth, max, &mem));
    
    //printf("done\n");
    SET_VECTOR_ELT(info, ni++, ScalarInteger(mem));
  }
  assert(ni ==actninfo);
  UNPROTECT(2); // info + name

  //assert(false);
  // printf("end\n");
  return info;
}

SEXP GetRegisterInfo(SEXP keynr, SEXP Ignoreactive, SEXP max_elements) 
{ // extended
    int knr;   
  knr = INTEGER(keynr)[0];
  if ((knr<0) || (knr>=MAXKEYS)) {
    return allocVector(VECSXP, 0);
  }
  // PrintMethodInfo(KEY[knr].meth);
  return InternalGetKeyInfo(&(KEY[knr]), LOGICAL(Ignoreactive)[0], 0, 
			    INTEGER(max_elements)[0]);
}


void leer(int level){
  char format[255];
//  printf("level=%d\n", level);
  sprintf(format,"%%%ds", -level * 3);
  PRINTF(format, "");
}

int MAX_PMI;
void PMI(cov_model *cov, int level)
{
     
  int i, j, endfor;
  cov_fct *C = CovList + cov->nr;
#define MNlength 5
  char MN[Forbidden + 1][MNlength],
      TriNames[TriBeyond + 1][12] = {"false", "true", "false (max)", "NaN"};
  //int n = 2;

  

  for (i=0; i<=Forbidden; i++) {
    strcopyN(MN[i], METHODNAMES[i], MNlength);
  }

  leer(level);
  cov_fct *CC = C;
  while(strcmp(CC->name, InternalName) ==0) CC--;
  PRINTF(">>> %s <<< [%d]", CC->name, cov->nr);
  PRINTF("\n");

  leer(level); PRINTF("param\n");

  // printf(">> %d %s %d\n", cov->nr, C->name, C->kappas);

  for(i=0; i<C->kappas; i++) {
    leer(level + 1); PRINTF("%-10s ", C->kappanames[i]);
    if (cov->p[i] == NULL) {
      PRINTF(" NULL\n");
      continue;
    }
    if (C->kappatype[i] == REALSXP) {
      if (cov->ncol[i]==1) {
	if (cov->nrow[i]==1) {
	  PRINTF("%f", cov->p[i][0]); 
	} else {
	  PRINTF("[%d] ", cov->nrow[i]);
	  endfor = cov->nrow[i]; if (endfor > MAX_PMI) endfor = MAX_PMI;
	  for (j=0; j<endfor; j++)
	    PRINTF(" %f", cov->p[i][j]); 
	}
      } else {
	PRINTF("[%d, %d] ", cov->nrow[i],  cov->ncol[i]);
	endfor = cov->nrow[i] * cov->ncol[i]; 
	if (endfor > MAX_PMI) endfor = MAX_PMI;
	for (j=0; j<endfor; j++)
	  PRINTF(" %f", cov->p[i][j]); 
      }
    } else if (C->kappatype[i] == INTSXP) {
      if (cov->ncol[i]==1) {
	if (cov->nrow[i]==1) {
	  PRINTF("%d", ((int*) cov->p[i])[0]); 
	} else {
	  PRINTF("[%d] ", cov->nrow[i]);
	  endfor = cov->nrow[i]; 
	  if (endfor > MAX_PMI) endfor = MAX_PMI;
	  for (j=0; j<endfor; j++)
	    PRINTF(" %d", ((int*) cov->p[i])[j]); 
	}
      } else {
	PRINTF("[%d, %d] ", cov->nrow[i],  cov->ncol[i]);
	endfor = cov->nrow[i] * cov->ncol[i];
	if (endfor > MAX_PMI) endfor = MAX_PMI;
	for (j=0; j<endfor; j++)
	  PRINTF(" %d", ((int*) cov->p[i])[j]); 
      }
    } else if (C->kappatype[i] == LISTOF + REALSXP) {
      int k;
      listoftype *p=(listoftype*) (cov->p[i]);	
      PRINTF("list [%d]\n", cov->nrow[i]);
      leer(level + 2); 
      for (k=0; k<cov->nrow[i]; k++) {
	if (p->ncol[k]==1) {
          if (p->nrow[k]==1) {
	     PRINTF("%f", p->p[k][0]); 
	  } else {
	      PRINTF("[%d] ", p->nrow[k]);
	    endfor = p->nrow[k]; if (endfor > MAX_PMI) endfor = MAX_PMI;
	    for (j=0; j<endfor; j++)
	      PRINTF(" %f", p->p[k][j]); 
	  }
        } else {
          PRINTF("[%d, %d] ", p->nrow[k],  p->ncol[k]);
	  endfor = p->nrow[k] * p->ncol[k]; 
	  if (endfor > MAX_PMI) endfor = MAX_PMI;
	  for (j=0; j<endfor; j++)
	    PRINTF(" %f", p->p[k][j]); 
	}
      }
      PRINTF("\n"); 
    } else {
      assert(false);
    }
    PRINTF("\n");
  }

  leer(level); PRINTF("%-10s [%d]","internal-q", cov->qlen);  
  endfor = cov->qlen; if (endfor > MAX_PMI) endfor = MAX_PMI;
  for (i=0; i<endfor; i++) PRINTF(" %f", cov->q[i]); 
  PRINTF("\n");

  
  leer(level); PRINTF("%-10s %s\n","calling", cov->calling==NULL ?
    "NULL" : CovList[cov->calling->nr].name);  
  leer(level); PRINTF("%-10s %d (%s)\n","statIn", 
		      cov->statIn, STATNAMES[(int) cov->statIn]);  
  leer(level); PRINTF("%-10s %d (%s)\n","isoIn", 
		      cov->isoIn, ISONAMES[(int) cov->isoIn]);
//  leer(level); PRINTF("%-10s %d\n","naturalscaling", cov->naturalscaling);  
  leer(level); PRINTF("%-10s %d\n","tsdim", cov->tsdim);  
  leer(level); PRINTF("%-10s %d\n","xdim", cov->xdim);  
  leer(level); PRINTF("%-10s %d\n","vdim", cov->vdim);  
  leer(level); PRINTF("%-10s %d\n","maxdim", cov->maxdim);  
  leer(level); PRINTF("%-10s %d\n","derivatives", cov->derivatives);  
  leer(level); PRINTF("%-10s %d\n","normalmix", (int) cov->normalmix);  
  leer(level); PRINTF("%-10s %d\n","finiterng", (int) cov->finiterange);  
  leer(level); PRINTF("%-10s %d\n","diag", (int) cov->diag);  
  leer(level); PRINTF("%-10s %d\n","ssep.last", (int) cov->semiseparatelast);  
  leer(level); PRINTF("%-10s %d\n","sep.last", (int) cov->separatelast);  
  leer(level); PRINTF("%-10s %d\n","tbm2num", (int) cov->tbm2num);  
  leer(level); PRINTF("%-10s %s\n","anyNAdown", TriNames[(int) cov->anyNAdown]);  
  leer(level); PRINTF("%-10s %s\n","anyNAscaleup", 
		      TriNames[(int) cov->anyNAscaleup]);
  leer(level); PRINTF("%-10s %d\n","spec:nmetro", cov->spec.nmetro);
  leer(level); PRINTF("%-10s %f\n","spec:sigma", cov->spec.sigma);


  leer(level); PRINTF("%-10s %d\n", "MLE", cov->MLE==NULL);
 
  leer(level); PRINTF("%-10s ", "pref");  
  for (i=0; i<=Direct; i++) PRINTF("%s:%d ", MN[i], (int) cov->pref[i]);
  PRINTF("\n"); leer(level); PRINTF("%-10s ", "");  
  for (; i<=Nothing; i++) PRINTF("%s:%d ", MN[i], (int) cov->pref[i]);
  PRINTF("\n");
		 
  leer(level); PRINTF("%-10s ", "user");  
  for (i=0; i<=Direct; i++) PRINTF("%s:%d ", MN[i], (int) cov->user[i]);
  PRINTF("\n"); leer(level); PRINTF("%-10s ", "");  
  for (; i<=Nothing; i++) PRINTF("%s:%d ", MN[i], (int) cov->user[i]);
  PRINTF("\n");

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] == NULL) {
//      PRINTF(" NULL\n");
      continue;
    }
//    leer(level); PRINTF("%s%-10s %d %d:\n",C->name, "submodel", i, level);  
    leer(level); PRINTF("%-10s %d:", "submodel", i);  
    
//    printf("%d \n", cov->sub[i]->nr);

    PMI(cov->sub[i], level + 1);
  }
}

void PrintModelInfo(cov_model *cov) { // OK
  PRINTF("\n");
  if (cov == NULL) {
    PRINTF("Covariance model is empty.\n");
  }
  else PMI(cov, 0);
}

/*
  leer(level); PRINTF("%-10s %d\n","stop", (int) meth->stop);  
  leer(level); PRINTF("%-10s %s\n","distrib", DISTRNAMES[loc->distribution]); 
  leer(level); PRINTF("%-10s %d\n","T", (int) loc->T);  
  leer(level); PRINTF("%-10s %d\n","x", (int) loc->x);  
  leer(level); PRINTF("%-10s %d\n","totalpts", (int) loc->totalpoints);  

  leer(level); PRINTF("%-10s ", "length"); 
  for (j=0; j<loc->timespacedim; j++) {
    PRINTF("%d ", loc->length[j]);  
  }
  PRINTF("\n");

  leer(level); PRINTF("%-10s %d\n","exp.#.simu",
		      (int) meth->expected_number_simu);  
*/

 

void PMeI(method_type *meth, int level){
//  printf("XXX %d\n", meth);  printf("XXX %d\n", meth->cov);
//  printf("XXX %d\n", meth->cov->tsdim);
    cov_model *cov = meth->cov;
  int i,j, dim;
  PRINTF("\n");
  //leer(level); PRINTF("%-10s %d\n","gp", (int) meth->gp);  
  //leer(level); PRINTF("%-10s %d\n","gpdo", (int)  meth->gpdo);  
  //leer(level); PRINTF("%-10s %d\n","destruct", (int) meth->destruct);  
  //leer(level); PRINTF("%-10s %d\n","S", (int) meth->S);  
  //leer(level); PRINTF("%-10s %d\n","compatible", (int) meth->compatible);  
  if (cov == NULL) {
    dim = -1;
    leer(level);  PRINTF("cov model is empty\n");
  } else {
    dim = cov->tsdim;
    leer(level); PRINTF("cov [%d]", cov);
    PMI(cov, level+1);
  }
  if (meth->caniso == NULL) {
      leer(level); PRINTF("%-10s %s\n",  "caniso", "NULL"); 
  } else {
      for (i=0; i<dim; i++) {
	  leer(level); PRINTF("%-10s ", i==0 ? "caniso" : ""); 
	  for (j=0; j<dim; j++) {
	      PRINTF("%2.2f ", meth->caniso[i + j * dim]);  
	  }
	  PRINTF("\n");
      }
  }
  leer(level); PRINTF("%-10s %f\n","cvar", meth->cvar);  
  leer(level); PRINTF("%-10s %d\n","matr.type", (int) meth->type);  
  //leer(level); PRINTF("%-10s %d\n","grid", (int) ,loc->grid);  
  //leer(level); PRINTF("%-10s %d\n","xdim", dim);  
//  leer(level); PRINTF("%-10s %ld\n","S", 
//		      (POINTER) meth->S);  
//  leer(level); PRINTF("%-10s %ld\n","space",
//		      (POINTER) meth->space);  
//  leer(level); PRINTF("%-10s %ld\n","sptime", 
//		      (POINTER) meth->sptime);  
  //leer(level); PRINTF("%-10s %d\n","grani", -1);  
//  leer(level); PRINTF("%-10s %ld\n","domethod", 
//		      (POINTER) meth->domethod);
  leer(level); PRINTF("%-10s %d (%s) \n","nr", meth->nr,
		      (meth->nr < 0) ? METHODNAMES[-meth->nr-1] : 
		      (meth->nr > 1000) ? "NOT SET" : CovList[meth->nr].name);  
//  leer(level); PRINTF("%-10s %ld\n","hanging", 
//		      (POINTER) meth->hanging);  
  leer(level); PRINTF("%-10s %d\n","xdimout", meth->xdimout);  

  if (level < 5) {
    if (meth->nsub == 0) {
	leer(level); PRINTF("%-10s (none)\n","submethods");  
    } else {
      for (j=0; j<meth->nsub; j++) {
	assert(meth->sub[j] != NULL);
	leer(level); PRINTF("%-10s %d", "submeth", j); 
	PMeI(meth->sub[j], level +1);
      }
    }
  }

  if (false && level == 0) {
    for (i=0; i<Forbidden; i++)
      PRINTF("%-15s init:%d do:%d\n", METHODNAMES[i], init_method[i],
	     do_method[i]);
  }

}

void PrintMethodInfo(method_type *meth) {
  PRINTF("\n");
  if (meth == NULL) 
    PRINTF("no method given\n");
  else PMeI(meth, 0);
}


