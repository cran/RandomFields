/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- get key strukture

 Copyright (C) 2001 -- 2014 Martin Schlather, 

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

name_type FT = {"false", "true"},
  TriNames = {"mismatch", 
	      "dep. next model",
	      "dep. prev. model",
	      "param. dependent",
	      "false", "true", 
	      "normal mixture",
	      "NaN"};  

//#define MAX_INT 2147483647
#define MAX_INT   2000000000


SEXP String(char *V) {
  //int i;
  SEXP str;
  PROTECT(str = allocVector(STRSXP, 1)); 
  SET_STRING_ELT(str, 1, mkChar(V));
  UNPROTECT(1);
  return str;
}


SEXP Logi(bool* V, int n, int max) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(VECSXP, 0);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(LGLSXP, n));
  for (i=0; i<n; i++) LOGICAL(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}
SEXP Logi(bool* V, int n) {
  return Logi(V, n, MAX_INT);
}

SEXP Num(double* V, int n, int max) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(REALSXP, 0);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(REALSXP, n));
  for (i=0; i<n; i++) REAL(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}
SEXP Num(double* V, int n) {
  return Num(V, n, MAX_INT);
}

SEXP Result(res_type* V, int n, int max) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(REALSXP, 0);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(REALSXP, n));
  for (i=0; i<n; i++) REAL(dummy)[i] = (double) V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP Result(res_type* V, int n) {
  return Result(V, n, MAX_INT);
}

SEXP Int(int *V, int n, int max) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(INTSXP, 0);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(INTSXP, n));
  for (i=0; i<n; i++) INTEGER(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP Int(int* V, int n) {
  return Int(V, n, MAX_INT);
}

SEXP Char(const char **V, int n, int max) {
  int i;
  SEXP dummy;
  if (V==NULL) return allocVector(STRSXP, 0);
  if (n>max) return TooLarge(&n, 1);
  PROTECT(dummy=allocVector(STRSXP, n));
  for (i=0; i<n; i++) SET_STRING_ELT(dummy, i, mkChar(V[i]));  
  UNPROTECT(1);
  return dummy;
}

SEXP Char(const char **V, int n) {
  return Char(V, n, MAX_INT);
}

SEXP Mat(double* V, int row, int col, int max) {
  int i, n;
  if (V==NULL) return allocMatrix(REALSXP, 0, 0);
  n = row * col;
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

SEXP Mat(double* V, int row, int col) {
  return Mat(V, row, col, MAX_INT);
}

SEXP ResultMat(res_type* V, int row, int col, int max) {
  int i, n;
  if (V==NULL) return allocMatrix(REALSXP, 0, 0);
  n = row * col;
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

SEXP ResultMat(res_type* V, int row, int col) {
  return ResultMat(V, row, col, MAX_INT);
}

SEXP MatInt(int* V, int row, int col, int max) {
  int i, n;
  if (V==NULL) return allocMatrix(INTSXP, 0, 0);
  n = row * col;
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

SEXP MatInt(int* V, int row, int col) {
  return MatInt(V, row, col, MAX_INT);
}

SEXP Array3D(double** V, int depth, int row, int col, int max) {
  int i, j, m, n;
  if (V==NULL) return alloc3DArray(REALSXP, 0, 0, 0);
  m = row * col;
  n = row * col * depth;
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

SEXP Array3D(double** V, int depth, int row, int col) {
  return Array3D(V, depth, row, col, MAX_INT);
}

SEXP TooLarge(int *n, int l){
#define nTooLarge 2 // mit op
  const char *tooLarge[nTooLarge] = {"size", "msg"};
  int i;
  SEXP namevec, info;
  PROTECT(info=allocVector(VECSXP, nTooLarge));
  PROTECT(namevec = allocVector(STRSXP, nTooLarge));
  for (i=0; i<nTooLarge; i++) SET_STRING_ELT(namevec, i, mkChar(tooLarge[i]));
  setAttrib(info, R_NamesSymbol, namevec);
  i=0;
  SET_VECTOR_ELT(info, i++, Int(n, l, l));
  SET_VECTOR_ELT(info, i,
		 mkString("too many elements - increase max.elements"));
  UNPROTECT(2);
  return info;
}


SEXP Param(void* p, int nrow, int ncol, SEXPTYPE type, bool drop) {
  int i;
  if (p == NULL) {
    return allocVector(REALSXP, 0); 
  } else {
    switch(type) {
    case REALSXP :
      return (ncol==1 && drop) ? Num((double*) p, nrow) 
	: Mat((double*) p, nrow, ncol);
    case INTSXP :
      return (ncol==1 && drop) ? Int((int*) p, nrow)
	: MatInt((int*) p, nrow, ncol);
    case CLOSXP :
      BUG;
    case STRSXP :
      return String((char*) p);
    case LANGSXP : 
      {
	const char *msg[1] = {"expression given by the user"};
	return Char(msg, 1);
      }
      break; // 
    default:
      if (type >= LISTOF){
	SEXP dummy = NULL;	
	PROTECT(dummy = allocVector(VECSXP, nrow));
	listoftype *q = (listoftype *) p;
	for (i=0; i<nrow; i++) {
	  SET_VECTOR_ELT(dummy, i, 
			 Param(q->p[i], q->nrow[i], q->ncol[i], REALSXP,false));
	}
	UNPROTECT(1);
	return dummy;
      } else {
	BUG;
      }
    } // switch
  }
  BUG;
  return NULL;
}

  
#define nsimuinfo 3
SEXP GetSimuInfo(simu_type *simu) {
  if (simu == NULL) return allocVector(VECSXP, 0);
  const char *info[nsimuinfo] = 
    {"active", "pair", "expect.simu"};
  SEXP namevec, l;
  int k;

  PROTECT(l = allocVector(VECSXP, nsimuinfo));
  PROTECT(namevec = allocVector(STRSXP, nsimuinfo));
  for (k=0; k<nsimuinfo; k++)
    SET_STRING_ELT(namevec, k, mkChar(info[k]));

  k = 0;  
  SET_VECTOR_ELT(l, k++, ScalarLogical(simu->active));
  SET_VECTOR_ELT(l, k++, ScalarLogical(simu->pair));
  SET_VECTOR_ELT(l, k++, ScalarInteger(simu->expected_number_simu));
  assert(k==nsimuinfo); 
  
  setAttrib(l, R_NamesSymbol, namevec);
  UNPROTECT(2); // l + namelvec
  assert(k == nsimuinfo);

  return l;
}


#define nlocinfo 14
SEXP GetLocationInfo(location_type *loc) {
  if (loc == NULL) return allocVector(VECSXP, 0);
  const char *info[nlocinfo] = 
    {"timespacedim", "xdimOZ", "length", "spatialdim", "spatialtotpts", 
     "totpts", "distances", "grid", "Time", "xgr", "x", "T", "ygr", "y"};
  SEXP namevec, l;
  int k,
    nloc = nlocinfo,
    tsdim = loc->timespacedim,
    spdim = loc->spatialdim;

  if (loc->ly <= 0) nloc -= 2;
  PROTECT(l = allocVector(VECSXP, nloc));
  PROTECT(namevec = allocVector(STRSXP, nloc));
  for (k=0; k<nloc; k++)
    SET_STRING_ELT(namevec, k, mkChar(info[k]));

  k = 0;
  SET_VECTOR_ELT(l, k++, ScalarInteger(tsdim));
  SET_VECTOR_ELT(l, k++, ScalarInteger(loc->xdimOZ));
  SET_VECTOR_ELT(l, k++, Int(loc->length, tsdim));
  SET_VECTOR_ELT(l, k++, ScalarInteger(loc->spatialdim));
  SET_VECTOR_ELT(l, k++, ScalarInteger((int) loc->spatialtotalpoints));
  SET_VECTOR_ELT(l, k++, ScalarInteger((int) loc->totalpoints));
  SET_VECTOR_ELT(l, k++, ScalarLogical(loc->distances));
  SET_VECTOR_ELT(l, k++, ScalarLogical(loc->grid));
  SET_VECTOR_ELT(l, k++, ScalarLogical(loc->Time));
  SET_VECTOR_ELT(l, k++, 
		 Mat(loc->xgr[0], loc->grid ? 3 : 0, 
		     spdim)); // "xgr"

  // printf("spatial %d %d\n", loc->spatialtotalpoints, loc->lx);
  SET_VECTOR_ELT(l, k++, 
		 Mat(loc->x, loc->xdimOZ, 
		     loc->grid ? 0 : 
   		        loc->distances ? loc->lx * (loc->lx - 1) / 2 : loc->lx, 
		     MAX_INT)); //"x"
  SET_VECTOR_ELT(l, k++, Num(loc->T, loc->Time ? 3 : 0));// "T"
  if (loc->ly > 0) {   
    
    //printf("loc->ly %f\n", loc->ly);
 
    if (loc->distances) BUG;
    SET_VECTOR_ELT(l, k++, 
		   Mat(loc->ygr[0], loc->grid ? 3 : 0, spdim));
    SET_VECTOR_ELT(l, k++, 
		   Mat(loc->y, loc->xdimOZ,
		       loc->grid ? 0 : loc->ly));
  } else {
    if (loc->ygr[0] != NULL || loc->y != NULL) BUG;
  }
  
  setAttrib(l, R_NamesSymbol, namevec);
  UNPROTECT(2); // l + namelvec
  assert(k == nloc);
  return l;
}


SEXP GetModelInfo(cov_model *cov, int prlevel, int spConform, 
		  int whichSub, int Level) {
  // whichSub:  0=submodels, 1=keys, 2=both

#define ninfo0 4
#define ninfo1 6
#define ninfo2 8
#define ninfo3 10
#define ninfo4 3

  /*   !!!!!     ACHTUNG     !!!!!
       Wenn diese Funktion geaendert wird,
       muss auch GetExtModelInfo geaendert werden
       !!!!!                 !!!!!
  */

  //  print("1st gmi %s\n", NICK(cov));
  if (cov == NULL) return allocVector(VECSXP, 0);
  SEXP model, submodels, nameMvec, param, pnames;
  int i, j, nmodelinfo, k = 0;
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  location_type
    *loc = cov->calling == NULL ? cov->prevloc : Loc(cov),
    *callingloc = cov->calling == NULL || Level == 0 ? NULL : Loc(cov->calling);
  bool
    given_key = cov->Splus != NULL || cov->key != NULL,
    return_key = given_key && whichSub != 0,
    return_sub = cov->nsub > 0 && (whichSub != 1 || !given_key);
 
  nmodelinfo = ninfo0;
  switch(prlevel > 4 ? 4 : prlevel) {
  case 4 : nmodelinfo += ninfo4;
  case 3 : nmodelinfo += ninfo3;
  case 2 : nmodelinfo += ninfo2;
  case 1 : nmodelinfo += ninfo1;
  default: {}
  }
  if (!return_sub) nmodelinfo--;
  if (!return_key) nmodelinfo--; 

  PROTECT(model = allocVector(VECSXP, nmodelinfo));
  PROTECT(nameMvec = allocVector(STRSXP, nmodelinfo));

  SET_STRING_ELT(nameMvec, k, mkChar("name"));
  cov_fct *CC = CovList + cov->nr; // nicht gatternr
  while(strncmp(CC->name, InternalName, strlen(InternalName)) ==0) CC--;

  if (spConform) {
    char name[MAXCHAR+2];
    sprintf(name, "%s", CC->nick); 
    SET_VECTOR_ELT(model, k++, mkString(name));
  } else {
    SET_VECTOR_ELT(model, k++, mkString(CC->name));
  }
  SET_STRING_ELT(nameMvec, k, mkChar("param"));  
  int notnull = 0;
  for (i=0; i<C->kappas; i++) {
    if (cov->nrow[i]>0 && cov->ncol[i]>0) notnull++;
  }
  
  PROTECT(param = allocVector(VECSXP, notnull));
  PROTECT(pnames = allocVector(STRSXP, notnull));
 
  for (j=i=0; i<C->kappas; i++) {
    if (cov->nrow[i]>0 && cov->ncol[i]>0) {
      if (isAnyDollar(cov) && i==DANISO) {
	SET_STRING_ELT(pnames, j, mkChar("Aniso"));
	double
	  *Aniso = (double*) MALLOC(cov->nrow[i]*cov->ncol[i] * sizeof(double));
	int t,l,m;
	for (t=l=0; l<cov->nrow[i]; l++) {
	  for (m=0; m<cov->ncol[i]; m++) {
	    Aniso[t++] = P(DANISO)[m * cov->nrow[i] + l];
	  }
	}
	SET_VECTOR_ELT(param, j,
		       Param((void*) Aniso, cov->nrow[i], cov->ncol[i], 
			     C->kappatype[i], true));    
	free (Aniso);
      }
      SET_STRING_ELT(pnames, j, mkChar(!strcmp(C->kappanames[i], FREEVARIABLE)
				       && cov->ownkappanames[i] != NULL
				       ? cov->ownkappanames[i] 
				       : C->kappanames[i]));
      SET_VECTOR_ELT(param, j,
		     Param((void*) cov->px[i], cov->nrow[i], cov->ncol[i], 
			   C->kappatype[i], true));    
      j++;
    }
  }
  setAttrib(param, R_NamesSymbol, pnames);
  //  print("gmi !\n");

  SET_VECTOR_ELT(model, k++, param);
  UNPROTECT(2);

  

  //  goto END;

  //  print("start GMI 1 %d\n", prlevel);
  
  if (prlevel>=1) {      
    SET_STRING_ELT(nameMvec, k, mkChar("covnr"));
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->nr));

    SET_STRING_ELT(nameMvec, k, mkChar("vdim"));
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->vdim));

    //  SET_STRING_ELT(nameMvec, k, mkChar("naturalscaling"));  
    //  SET_VECTOR_ELT(model, k++, ScalarInteger(cov->naturalscaling));

    SET_STRING_ELT(nameMvec, k, mkChar("tsdim"));  
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->tsdim));
      
    SET_STRING_ELT(nameMvec, k, mkChar("xdimprev"));  
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->xdimprev));

    SET_STRING_ELT(nameMvec, k, mkChar("xdimown"));  
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->xdimown));

    SET_STRING_ELT(nameMvec, k, mkChar("indep.of.x"));  
    SET_VECTOR_ELT(model, k++, ScalarLogical(cov->matrix_indep_of_x));
  } 
   
  //  print("GMI 2\n");

  if (prlevel>=2) {
    SET_STRING_ELT(nameMvec, k, mkChar("type"));  
    SET_VECTOR_ELT(model, k++, Char(TYPENAMES + cov->typus, 1));

    SET_STRING_ELT(nameMvec, k, mkChar("role"));  
    SET_VECTOR_ELT(model, k++, Char(ROLENAMES + cov->role, 1));

    SET_STRING_ELT(nameMvec, k, mkChar("domown"));  
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->domown));
      
    SET_STRING_ELT(nameMvec, k, mkChar("isoown"));
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->isoown));
 
    SET_STRING_ELT(nameMvec, k, mkChar("internalq"));  
    SET_VECTOR_ELT(model, k++, Num(cov->q, cov->qlen));

    SET_STRING_ELT(nameMvec, k, mkChar("pref"));  
    SET_VECTOR_ELT(model, k++, Int(cov->pref, Nothing + 1));

    SET_STRING_ELT(nameMvec, k, mkChar("simu"));  
    SET_VECTOR_ELT(model, k++, GetSimuInfo(&(cov->simu)));
	
    SET_STRING_ELT(nameMvec, k, mkChar("loc"));  
    SET_VECTOR_ELT(model, k++, 
		   loc == callingloc ? mkString("see calling model, or there are no location information given at all") :
		   GetLocationInfo(loc));

  } 

  //  print("GMI 3\n");
  if (prlevel>=3) {
    SET_STRING_ELT(nameMvec, k, mkChar("logspeed"));
    SET_VECTOR_ELT(model, k++, ScalarReal(cov->logspeed));

    SET_STRING_ELT(nameMvec, k, mkChar("maxdim"));  
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->maxdim));
      
    SET_STRING_ELT(nameMvec, k, mkChar("full_derivs"));  
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->full_derivs));
      
    SET_STRING_ELT(nameMvec, k, mkChar("loggiven"));  
    SET_VECTOR_ELT(model, k++, ScalarLogical(cov->loggiven));

    SET_STRING_ELT(nameMvec, k, mkChar("monotone"));  
    SET_VECTOR_ELT(model, k++,
		   Char(MONOTONE_NAMES + cov->monotone - MISMATCH, 1));
            
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
    //  SET_VECTOR_ELT(model, k++, Int(cov->idx, cov->tsdim));
    //      SET_STRING_ELT(nameMvec, k, mkChar("user"));  
    // SET_VECTOR_ELT(model, k++, Int(cov->user, Nothing + 1));
  }

  if (prlevel>=4) {

    mpp_properties *mpp = &(cov->mpp);
 
    //   SET_STRING_ELT(nameMvec, k, mkChar("mpp.refradius"));  
    //  SET_VECTOR_ELT(model, k++, ScalarReal(mpp->refradius));
 
    SET_STRING_ELT(nameMvec, k, mkChar("mpp.maxheight"));  
    SET_VECTOR_ELT(model, k++, ScalarReal(mpp->maxheight));
  
    SET_STRING_ELT(nameMvec, k, mkChar("mpp.M"));
    SET_VECTOR_ELT(model, k++, cov->mpp.moments == 0 ? R_NilValue :
		   Num(mpp->M, cov->mpp.moments,  MAX_INT));

    SET_STRING_ELT(nameMvec, k, mkChar("mpp.Mplus"));
    SET_VECTOR_ELT(model, k++, cov->mpp.moments == 0 ? R_NilValue :
		   Num(mpp->Mplus, cov->mpp.moments,  MAX_INT));
  } 
   		 

  if (return_key) {
    if (cov->key != NULL) {
      SET_STRING_ELT(nameMvec, k, mkChar("key"));  
      SET_VECTOR_ELT(model, k++, GetModelInfo(cov->key, prlevel, spConform, 
					      whichSub, Level + 1));
    } else {  /// cov->nr == PLUS && cov->Splus != NULL
      int ii, n,
	subs = C->maxsub;
      SEXP keys;
      SET_STRING_ELT(nameMvec, k, mkChar("key"));  
      for (ii=n=0; ii<subs; ii++) if (cov->Splus->keys[ii] != NULL) n++;
      PROTECT(keys = allocVector(VECSXP, n));
      for (ii=n=0; ii<subs; ii++) 
	if (cov->Splus->keys[ii] != NULL)
	  SET_VECTOR_ELT(keys, n++, 
			 GetModelInfo(cov->Splus->keys[ii], prlevel, spConform,
				      whichSub, Level + 1));
      SET_VECTOR_ELT(model, k++, keys);
      UNPROTECT(1);
    }
  }
 
  if (return_sub) {
    SET_STRING_ELT(nameMvec, k, mkChar("submodels"));  
    //   print("gmi 1\n");
    //   PROTECT(namesubmodels = allocVector(STRSXP, cov->nsub));
    PROTECT(submodels = allocVector(VECSXP, cov->nsub));
    int zaehler = 0;
    //    print("gmi 2\n");
    for (i=0; i<MAXSUB; i++) {
      //	print("gmi x %d %d\n", i, MAXSUB);
      // print("\n\n %d %i %s\n", zaehler, i, NICK(cov));
      if (cov->sub[i] != NULL) {
	//	print("gmi a\n");
	//	SET_STRING_ELT(namesubmodels, zaehler, mkChar(C->subnames[i]));
	SET_VECTOR_ELT(submodels, zaehler, 
		       GetModelInfo(cov->sub[i], prlevel, spConform, 
				    whichSub, Level + 1));
	if (++zaehler >= cov->nsub) break;
      }
    }
    // setAttrib(submodels, R_NamesSymbol, namesubmodels);
    SET_VECTOR_ELT(model, k++, submodels);
    UNPROTECT(1);
  }

  //  SET_STRING_ELT(nameMvec, k, mkChar("method"));  
  //  SET_VECTOR_ELT(model, k++, mkChar(METHODNAMES[(int) cov->usermethod]));
  
  setAttrib(model, R_NamesSymbol, nameMvec);

  // END:
  //   printf("k %d %d prlevel=%d\n", k, nmodelinfo, prlevel);

  assert(k==nmodelinfo); 
  UNPROTECT(2); // model + namemodelvec
  return model;
}

SEXP GetExtModelInfo(SEXP keynr, SEXP Prlevel, SEXP spConform, SEXP whichSub) {
  int knr = INTEGER(keynr)[0],
    
    prlevel = INTEGER(Prlevel)[0] % 10;
//    # prlevel: >= 10 : wie 0-4(-9), jedoch ohne CALL_FCT zu loeschen
  bool 
    delete_call = INTEGER(Prlevel)[0] < 10;
  cov_model *cov, *orig;
  SEXP res, names;

  if (knr>=0 && knr <= MODEL_MAX && KEY[knr] != NULL) {
    orig = cov = KEY[knr];
    if (delete_call && isInterface(cov)) {
      cov = cov->key == NULL ? cov->sub[0] : cov->key;
    }
    res = GetModelInfo(cov, prlevel, (bool) INTEGER(spConform)[0], 
		       INTEGER(whichSub)[0], 0);
    if (prlevel>=1 && delete_call) {
      names = getAttrib(res, R_NamesSymbol);
      int i, len = length(names);
      for (i=0; i<len; i++) {
	const char *name = CHAR(STRING_ELT(names, i));
	if (strcmp("xdimprev", name) == 0) {
	  //  printf("xdim %d %d\n", INTEGER(VECTOR_ELT(res, i))[0],orig->xdimprev);
	  //assert(false);
	  INTEGER(VECTOR_ELT(res, i))[0] = orig->xdimprev;
	  break;
	}
      }
    }
    return res;
  }
  return allocVector(VECSXP, 0);
}


//#define nglobalinfo 0
//SEXP GetGlobalInfo(globalparam global) {
//  SEXP namevec, l; 
//
//  PROTECT(l = allocVector(VECSXP, nglobalinfo));
///  PROTECT(namevec = allocVector(STRSXP, nglobalinfo));
///  setAttrib(l, R_NamesSymbol, namevec);
//  assert(0 == nglobalinfo);
//  UNPROTECT(2);
//  return l;
//}

/*
  #define ntrendinfo 5
  SEXP GetTrendInfo(trend_type *trend) {
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
  SET_VECTOR_ELT(l, k++, Num(trend->LinearTrend, trend->lLinTrend));
  setAttrib(l, R_NamesSymbol, namevec);
  UNPROTECT(2); // l + namelvec
  assert(k == ntrendinfo);

  return l;
  }
*/




void leer(int level){
  char format[255];
  //  print("level=%d\n", level);
  sprintf(format,"%%%ds", -level * 3);
  PRINTF(format, "");
}

int MAX_PMI = 5;

void PrintPoints(location_type *loc, char *name, 
		 double *x, coord_type xgr, int lx) {
#define maxpts 100
  int i;
  if (loc->grid) {
    PRINTF("loc:%sgr    ");
    for (i=0; i<loc->timespacedim; i++) 
      PRINTF("(%3.3f, %3.3f, %2.0f) ", xgr[i][XSTART], xgr[i][XSTEP], 
	     xgr[i][XLENGTH]);
  } else {
    PRINTF("loc:%s      ", name);

    //  printf("%ld %d\n", loc->x, loc->lx);

    if (loc->lx == 0) {
      PRINTF("not given! (%d)", addressbits(loc->x));
    } else {
      int total = loc->distances ? lx * (lx-1) / 2 : lx * loc->xdimOZ,
	endfor = total;

      if (endfor > maxpts) endfor = maxpts;
      for (i=0; i<endfor; i++) {
	PRINTF("%4.3f", x[i]);
	if ((i+1) % loc->xdimOZ == 0) PRINTF(";");
	PRINTF(" ");
      }
      if (endfor < total) 
      PRINTF("... [%d not shown]", total - endfor);
    }
  }
  PRINTF("\n");
}

void PrintLoc(int level, location_type *loc, bool own) {
  int i;
  if (loc == NULL) {
    leer(level); PRINTF("%-10s %s\n", "loc:", "not given");
    return;
  }
  if (own) {
    leer(level); PRINTF("%-10s %d\n", "own is set:", addressbits(loc));
  }

  leer(level); PRINTF("%-10s %d %d %d\n","loc:ts,sp,xdimOZ", 
		      loc->timespacedim, loc->spatialdim, loc->xdimOZ);
  leer(level); PRINTF("%-10s ","loc:length");
  for (i=0; i<loc->timespacedim; i++) PRINTF("%d ", loc->length[i]);
  PRINTF("\n");
  leer(level); PRINTF("%-10s %d\n","loc:lx", loc->lx);
  leer(level); PRINTF("%-10s %d\n","loc:totpts", loc->totalpoints);
  leer(level); PRINTF("%-10s %s\n","loc:grid", FT[loc->grid]);
  leer(level); PRINTF("%-10s %s\n","loc:dist", FT[loc->distances]);
  leer(level); PRINTF("%-10s %s\n","loc:Time", FT[loc->Time]);
  leer(level); PrintPoints(loc, (char *) "x", loc->x, loc->xgr, loc->lx);
  if (loc->y!=NULL || loc->ygr[0]!=NULL) {  
    // printf("not null\n");
    leer(level); PrintPoints(loc, (char*) "y", loc->y, loc->ygr, loc->ly);
  }
  if (loc->Time) { 
    leer(level); PRINTF("%-10s (%f %f %f)\n", "loc:T", 
			loc->T[0], loc->T[1], loc->T[2]);
  }
  leer(level); PRINTF("%-10s ","loc:cansio");
  if (loc->caniso==NULL) PRINTF("null\n"); 
  else {
    int endfor = loc->cani_nrow * loc->cani_ncol;
    PRINTF(" [%d, %d] ",  loc->cani_nrow, loc->cani_ncol);
    if (endfor > MAX_PMI) endfor = MAX_PMI;
    for (i=0; i<endfor; i++) PRINTF(" %f", loc->caniso[i]); 
    PRINTF("\n");
  }
  //leer(level); PRINTF("%-10s %d\n","loc:stat", loc->domain);

}
  

static bool PMI_print_dollar = !true,
  PMI_print_mpp = true,
  PMI_print_pgs = !true, 
  PMI_print_details = !true,
  PMI_print_loc = true;
static int  PMI_print_rect = 0;  // 0, 1, 2

void pmi(cov_model *cov, char all, int level) {     
  int i, j, endfor;
  cov_fct *C = CovList + cov->nr; // nicht gatternr
#define MNlength 4
  char MN[Forbidden + 1][MNlength], name[100];
 
  //int n = 2;

  for (i=0; i<=Forbidden; i++) {
    strcopyN(MN[i], METHODNAMES[i], MNlength);
  }

  cov_fct *CC = C;
  while(strcmp(CC->name, InternalName) ==0) CC--;
  if (level == 0)  PRINTF("******   %s   ****** [%d,%d]", CC->nick, cov->nr, cov->zaehler);
  else PRINTF("    **** %s **** [%d,%d]", CC->nick, cov->nr, cov->zaehler);
  PRINTF("\n");

  leer(level); PRINTF("%-10s %s\n", "param", C->kappas == 0 ? "none" : ""); 
  // print(">> %d %s %d\n", cov->nr, C->name, C->kappas);

  for(i=0; i<C->kappas; i++) {
       //SEXP xs, ys, zs, ts;   
    //    printf("i=%d %d %d %ld %s %s\n", i, !strcmp(C->kappanames[i], FREEVARIABLE),
    //	cov->ownkappanames[i] != NULL, cov->p[i], cov->ownkappanames[i], C->kappanames[i]);

    strcpy(name, 
	   !strcmp(C->kappanames[i], FREEVARIABLE) && 
	           cov->ownkappanames != NULL && cov->ownkappanames[i] != NULL
	   ? cov->ownkappanames[i]
	   : C->kappanames[i]
	   );
    name[9] = '\0';
    leer(level + 1); PRINTF("%-10s", 
			    name); 
    if (PisNULL(i)) {
      PRINTF(" NULL");
    } else if (C->kappatype[i] == REALSXP) {
      if (cov->ncol[i]==1) {
	if (cov->nrow[i]==1) {
	  PRINTF("%f", P0(i)); 
	} else {
	  PRINTF("[%d] ", cov->nrow[i]);
	  endfor = cov->nrow[i]; if (endfor > MAX_PMI) endfor = MAX_PMI;
	  for (j=0; j<endfor; j++)
	    PRINTF(" %f", P(i)[j]); 
	}
      } else {
	PRINTF("[%d, %d] ", cov->nrow[i],  cov->ncol[i]);
	endfor = cov->nrow[i] * cov->ncol[i]; 
	if (endfor > MAX_PMI) endfor = MAX_PMI;
	for (j=0; j<endfor; j++)
	  PRINTF(" %f", P(i)[j]); 
      }
    } else if (C->kappatype[i] == INTSXP) {
      if (cov->ncol[i]==1) {
	if (cov->nrow[i]==1) {
	  PRINTF("%d", P0INT(i)); 
	} else {
	  PRINTF("[%d] ", cov->nrow[i]);
	  endfor = cov->nrow[i]; 
	  if (endfor > MAX_PMI) endfor = MAX_PMI;
	  for (j=0; j<endfor; j++)
	    PRINTF(" %d", PINT(i)[j]); 
	}
      } else {
	PRINTF("[%d, %d] ", cov->nrow[i],  cov->ncol[i]);
	endfor = cov->nrow[i] * cov->ncol[i];
	if (endfor > MAX_PMI) endfor = MAX_PMI;
	for (j=0; j<endfor; j++)
	  PRINTF(" %d", PINT(i)[j]); 
      }
    } else if (C->kappatype[i] == CLOSXP) {
      PRINTF("arbitrary function...");      
    } else if (C->kappatype[i] == LANGSXP) {
      PRINTF("arbitrary language expression ...");      
    } else if (C->kappatype[i] == LISTOF + REALSXP) {
      int k, ende=0;
      listoftype *p= PLIST(i);	
      PRINTF("list [%d]\n", cov->nrow[i]);
      leer(level + 2); 
      for (k=0; k<cov->nrow[i]; k++) {
	if (p->ncol[k]==1) {
          if (p->nrow[k]==1) {
	    PRINTF("%f", p->p[k][0]); 
	  } else {
	    PRINTF("[%d] ", p->nrow[k]);
	    ende = endfor = p->nrow[k]; 
	    if (endfor > MAX_PMI) endfor = MAX_PMI;
	    for (j=0; j<endfor; j++)
	      PRINTF(" %f", p->p[k][j]); 
	  }
        } else {
          PRINTF("[%d, %d] ", p->nrow[k],  p->ncol[k]);
	  ende = endfor = p->nrow[k] * p->ncol[k]; 
	  if (endfor > MAX_PMI) endfor = MAX_PMI;
	  for (j=0; j<endfor; j++) PRINTF(" %f", p->p[k][j]); 
	}
	if (ende > MAX_PMI) PRINTF ("...");
	if (k<cov->nrow[i]-1)  PRINTF("\n"); 
      }
    } else {
      assert(false);
    }
    if (cov->kappasub[i] != NULL) {
      PRINTF("  <=");
      pmi(cov->kappasub[i], all, level + 3);
    } else  PRINTF("\n");
  }
  if (cov->Sset != NULL) {
    cov_model *from = cov->Sset->remote;
    leer(level + 1); 
    PRINTF("%-10s%s [%d]\n", "<remote>", NICK(from), from->zaehler);     
  }

  //  for (; i<MAXPARAM; i++) if (cov->kappasub[i] != NULL) 
  //			    PRINTF(" %d != NULL !!\n", i);

  leer(level); PRINTF("%-10s [%d]","internal-q", cov->qlen);  

  endfor = cov->qlen; if (endfor > MAX_PMI) endfor = MAX_PMI;
  //  printf("qlen %d %d %d\n", cov->qlen, MAX_PMI, endfor);
  for (i=0; i<endfor; i++) PRINTF(" %f", cov->q[i]); 
  PRINTF("\n");

  if (cov->calling == NULL && level != 0) BUG;
  leer(level); PRINTF("%-10s %s\n","calling", cov->calling==NULL ?
		      "NULL" : NICK(cov->calling));  
  leer(level); PRINTF("%-10s %s (%d)\n","gatter",
		      cov->gatternr >=0 ? CovList[cov->gatternr].name : "none",
		      cov->gatternr);  
  if (cov->secondarygatternr >= 0) {
    leer(level);
    PRINTF("%-10s %s (%d)\n","secondary",
	   CovList[cov->secondarygatternr].name, cov->secondarygatternr);  
  }
  leer(level); PRINTF("%-10s %d/%d (%s/%s)\n","domprev/own", 
		      cov->domprev, cov->domown,
		      STATNAMES[(int) cov->domprev],
		      STATNAMES[(int) cov->domown]);
  leer(level); PRINTF("%-10s %d/%d (%s/%s)\n","isoprev/own", 
		      cov->isoprev, cov->isoown,
		      ISONAMES[(int) cov->isoprev],
		      ISONAMES[(int) cov->isoown]);
  leer(level); PRINTF("%-10s %d, %d:%d, %d/%d/%d\n","ts-x-v-dim",
		      cov->tsdim, cov->xdimprev, cov->xdimown, cov->vdim,
		      cov->vdim2[0], cov->vdim2[1]);  
   leer(level); PRINTF("%-10s %d\n","maxdim", cov->maxdim);  
  leer(level); PRINTF("%-10s %s (%d)\n", "type", TYPENAMES[cov->typus],
		      (int) cov->typus);
  leer(level); PRINTF("%-10s %s (%d)\n","role", ROLENAMES[cov->role],
		      (int) cov->role);
  leer(level); PRINTF("%-10s %s\n","method", METHODNAMES[cov->method]);
  leer(level); PRINTF("%-10s %s\n","initialised", FT[cov->initialised]);
  leer(level); PRINTF("%-10s %s\n","fieldret", FT[cov->fieldreturn]);
  leer(level); PRINTF("%-10s %s\n","origrf", FT[cov->origrf]);
  leer(level); PRINTF("%-10s %d\n","rf", addressbits(cov->rf));
  leer(level); PRINTF("%-10s ", "pref");  
  for (i=0; i<=Sequential; i++) PRINTF("%s:%d ", MN[i], (int) cov->pref[i]);
  PRINTF("\n"); leer(level); PRINTF("%-10s ", "");  
  for (; i<=Nothing; i++) PRINTF("%s:%d ", MN[i], (int) cov->pref[i]);
  PRINTF("\n");
  if (PMI_print_details) {
    leer(level); PRINTF("%-10s %s\n","determinst", FT[cov->deterministic]);  
    leer(level); PRINTF("%-10s %s\n","user_given", 
		      cov->user_given == ug_internal ? "internal" :
		      cov->user_given == ug_explicit ? "explicit" : "implicit");
    leer(level); PRINTF("%-10s %f\n","logspeed", cov->logspeed);  
    leer(level); PRINTF("%-10s %d, %d\n","full/rese deriv's", 
			cov->full_derivs, cov->rese_derivs);  
    leer(level); PRINTF("%-10s %s\n","loggiven", FT[cov->loggiven]);  
    leer(level); PRINTF("%-10s %s (%d)\n","monotone", 
			MONOTONE_NAMES[cov->monotone - MISMATCH], 
			cov->monotone);  
    leer(level); PRINTF("%-10s %s\n","finiterng", 
			TriNames[cov->finiterange - MISMATCH]);  
    leer(level); PRINTF("%-10s %d\n","storage", addressbits(cov->stor));
    leer(level); PRINTF("%-10s %s\n","simu:activ", FT[cov->simu.active]);
    leer(level); PRINTF("%-10s %s\n","simu:pair", FT[cov->simu.pair]);
    leer(level); PRINTF("%-10s %d\n","simu:expect", cov->simu.expected_number_simu);    
    leer(level); PRINTF("%-10s %s\n", "MLE", FT[cov->MLE==NULL]); 
  }
  //leer(level); PRINTF("%-10s %d\n","naturalscaling", cov->naturalscaling);  
  //leer(level); PRINTF("%-10s %d\n","init",CovList[cov->nr].init!=init_failed);
  //leer(level); PRINTF("%-10s %d\n","diag", (int) cov->diag);  
  //leer(level); PRINTF("%-10s %d\n","ssep.last", (int) cov->semiseparatelast);
  //leer(level); PRINTF("%-10s %d\n","sep.last", (int) cov->separatelast);  
  //leer(level); PRINTF("%-10s %d\n","tbm2num", (int) cov->tbm2num);  
  //leer(level); PRINTF("%-10s %d\n","spec:nmetro", cov->spec.nmetro);
  //leer(level); PRINTF("%-10s %f\n","spec:sigma", cov->spec.sigma);
  if (PMI_print_dollar && cov->Sdollar != NULL) {
    leer(level); PRINTF("%-10s %d\n","$:z", addressbits(cov->Sdollar->z));
    leer(level); PRINTF("%-10s %d\n","$:z2", addressbits(cov->Sdollar->z2));
    leer(level); PRINTF("%-10s %d\n","$:y", addressbits(cov->Sdollar->y));
  }


  if (PMI_print_mpp) {
    if (R_FINITE(cov->mpp.log_zhou_c)) {
      leer(level); PRINTF("%-10s %f\n","mpp:zhou_c", cov->mpp.log_zhou_c);
    } 
    leer(level); PRINTF("%-10s %f\n","mpp:maxhgt", cov->mpp.maxheight);    
    leer(level); PRINTF("%-10s ","mpp:M+");
    if (cov->mpp.Mplus == NULL) PRINTF("not initialized yet (size=%d)\n",
				       cov->mpp.moments);
    else {
      for (i=0; i<=cov->mpp.moments; i++) PRINTF("%f, ", cov->mpp.Mplus[i]);
      PRINTF("\n");    
    }
    leer(level); PRINTF("%-10s ","mpp:M");
    if (cov->mpp.Mplus == NULL) PRINTF("not initialized yet (size=%d)\n",
				       cov->mpp.moments);
    else {
      for (i=0; i<=cov->mpp.moments; i++) PRINTF("%f, ", cov->mpp.M[i]);
      PRINTF("\n");    
    }
    leer(level); PRINTF("%-10s %s\n","mpp:log",
			FT[CovList[cov->nr].log != ErrLogCov]);
    leer(level); PRINTF("%-10s %s\n","mpplgnonst", 
			FT[CovList[cov->nr].nonstatlog != ErrLogCovNonstat]);
    leer(level); PRINTF("%-10s %s\n","mpp:do",
			FT[CovList[cov->nr].Do != do_failed]);
  }

  
  if (PMI_print_pgs && cov->Spgs != NULL) {
    pgs_storage *pgs = cov->Spgs;
    int d,
      size = pgs->size,
      dim = cov->xdimown;

    leer(level); PRINTF("%-10s %f\n","pgs:mass", pgs->totalmass);
    leer(level); PRINTF("%-10s %f\n","pgs:logdens", pgs->log_density);
    
    leer(level); PRINTF("%-10s %d\n","pgs:size", size);
    leer(level); PRINTF("%-10s %d\n","pgs:dim", dim);

#define SHOWDEFAULT(Z, X, Y) if (pgs->X != NULL) {			\
      leer(level);  PRINTF("%-10s ",Z);					\
      for (d=0; d<dim; d++) PRINTF(Y, pgs->X[d]); PRINTF("\n");}
#define SHOW(Z, X) SHOWDEFAULT(Z, X, "%f ")
#define SHOWINT(Z, X) SHOWDEFAULT(Z, X, "%d ")


    SHOW("pgs:v", v);    
    SHOW("pgs:x", x);
    SHOW("pgs:xstart", xstart);
    SHOW("pgs:inc", inc);
    SHOW("pgs:suppmin", supportmin);
    SHOW("pgs:suppmax", supportmax);
    SHOWINT("pgs:gridlen", gridlen);
    SHOWINT("pgs:start", start);
    SHOWINT("pgs:end", end);
    SHOW("pgs:delta", delta);
    SHOWINT("pgs:nx", nx);
    if (pgs->pos != NULL) { // gauss
      location_type *loc = Loc(cov);
      leer(level); PrintPoints(loc, (char *) "pgs.x", loc->x, pgs->xgr,
			       loc->lx);	
      SHOW("pgs:y", y);
      SHOWINT("pgs:pos", pos);
      SHOWINT("pgs:min", min);
      SHOWINT("pgs:max", max);
    }      
    if (pgs->halfstepvector != NULL) { // max-stable 
      leer(level); PRINTF("%-10s %s\n","pgs:flat", FT[pgs->flat]);
      leer(level); PRINTF("%-10s %f\n","pgs:orig", pgs->value_orig);
      leer(level); PRINTF("%-10s %f\n","pgs:globmin", pgs->globalmin);
      leer(level); PRINTF("%-10s %f\n","pgs:cur.thres",pgs->currentthreshold);
      SHOW("pgs:half", halfstepvector);
      if (pgs->single != NULL) {
	leer(level); { PRINTF("%-10s ","pgs:single"); 
	  for (d=0; d<size; d++) PRINTF("%f ", pgs->single[d]); PRINTF("\n"); }
      }
      if (pgs->total != NULL) {
	leer(level); { PRINTF("%-10s ","pgs:total"); 
	  for (d=0; d<size; d++) PRINTF("%f ", pgs->total[d]); PRINTF("\n"); }
      }
    } else { // gauss oder poisson
      leer(level); PRINTF("%-10s %f\n","pgs:intens", pgs->intensity);
    }
  }
  
  if (PMI_print_rect) {
    int d;
    if (cov->Srect != NULL) {
      rect_storage *p = cov->Srect;
      int
	nstepP2 = p->nstep + 2,
	dim = cov->xdimown,
	dimP1 = dim + 1;
      
      leer(level); PRINTF("%-10s %f\n","rct:inner", p->inner);
      leer(level); PRINTF("%-10s %f\n","rct:in.cnst", p->inner_const);
      leer(level); PRINTF("%-10s %f\n","rct:in.pow", p->inner_pow);
      leer(level); PRINTF("%-10s %f\n","rct:outer", p->outer);
      leer(level); PRINTF("%-10s %f\n","rct:o.cnst", p->outer_const);
      leer(level); PRINTF("%-10s %f\n","rct:o.pow", p->outer_pow);
      leer(level); PRINTF("%-10s %f\n","rct:o,pow.c", p->outer_pow_const);
      leer(level); PRINTF("%-10s %f\n","rct:step", p->step);
      leer(level); PRINTF("%-10s %d\n","rct:nstep", p->nstep);
      leer(level); PRINTF("%-10s %d\n","rct:ntmp", p->tmp_n);
      
      if (PMI_print_rect > 1 && p->value != NULL) {
	leer(level); { PRINTF("%-10s ","rct:val"); 
	  for (d=0; d<nstepP2; d++) PRINTF("%4.3f ", p->value[d]); PRINTF("\n"); }
	leer(level); { PRINTF("%-10s ","rct:wght"); 
	  for (d=0; d<nstepP2; d++) PRINTF("%4.4f ", p->weight[d]);
	  PRINTF("\n"); }
	leer(level); { PRINTF("%-10s ","rct:tmp"); 
	  for (d=0; d<p->tmp_n; d++) PRINTF("%4.3f ", p->tmp_weight[d]);
	  PRINTF("\n"); }
	leer(level); { PRINTF("%-10s ","rct:right"); 
	  for (d=0; d<p->tmp_n; d++) PRINTF("%4.4f ", p->right_endpoint[d]); 
	  PRINTF("\n"); }
	leer(level); { PRINTF("%-10s ","rct:z"); 
	  assert(p->z != NULL);
	  for (d=0; d<dimP1; d++) PRINTF("%4.3f ", p->z[d]); PRINTF("\n"); }
	leer(level); { PRINTF("%-10s ","rct:squeezd"); 
	  for (d=0; d<p->tmp_n; d++) PRINTF("%d ", p->squeezed_dim[d]);
	  PRINTF("\n"); }
	leer(level); { PRINTF("%-10s ","rct:assign"); 
	  for (d=0; d<p->tmp_n; d++) PRINTF("%d ", p->assign[d]);
	  PRINTF("\n"); }
	leer(level); { PRINTF("%-10s ","rct:idx"); 
	  for (d=0; d<dimP1; d++) PRINTF("%d ", p->i[d]); PRINTF("\n"); }
      }
    }
    if (cov->taylorN > 0) {
      for (d=0; d<cov->taylorN; d++) {
	leer(level); 
	PRINTF("%-10s c%d=%f p%d=%f\n", d==0 ? "taylor" : "", 
	       d, cov->taylor[d][TaylorConst], d, cov->taylor[d][TaylorPow]);
      }
    } else {leer(level); PRINTF("%-10s %s\n", "taylor", "undetermined");}
    if (cov->tailN > 0) {
	for (d=0; d<cov->tailN; d++) {
	  leer(level);
	  PRINTF("%-10s c%d=%4.3f p%d=%4.3f ce%d=%4.3f pe%d=%4.3f\n", 
		 d==0 ? "tailtlr" : "",
		 d, cov->tail[d][TaylorConst], 
		 d, cov->tail[d][TaylorPow],
		 d, cov->tail[d][TaylorExpConst], 
		 d, cov->tail[d][TaylorExpPow]);
	}
    } else {leer(level); PRINTF("%-10s %s\n", "tailtlr", "undetermined");}
  }
    
  if (PMI_print_loc) {
    if (cov->ownloc != NULL) {
      PrintLoc(level, cov->ownloc, true);
    } else {
      leer(level);    
      if (cov->calling == NULL) {
	PRINTF("%-10s (%d)\n", "loc:extern", addressbits(cov->prevloc));
	PrintLoc(level, cov->prevloc, false);
      } else if (cov->prevloc == cov->calling->prevloc) {
	PRINTF("%-10s (%d)\n", "loc:calling->prev", addressbits(cov->prevloc));
	PrintLoc(level, cov->prevloc, false);
      } else if (cov->prevloc == cov->calling->ownloc) {
	PRINTF("%-10s (%d)\n", "loc:calling->own", addressbits(cov->prevloc));
	PrintLoc(level, cov->prevloc, false);
      } else  PRINTF("%-10s (%d)\n", "loc:MISMATCH", addressbits(cov->prevloc));
    }
  }
 
  bool givenkey = false;
  if (cov->key != NULL) {
    leer(level);
    givenkey = true;
    PRINTF("%-10s :", "key");  
    if (all >= 0) pmi(cov->key, all, level + 1);
  }

  if (cov->Splus != NULL) {
    givenkey = true;
    int ii;
    if (all >= 0) {
      for (ii=0; ii < cov->nsub; ii++) {
	cov_model *key = cov->Splus->keys[ii];
	leer(level);      
	if (key != NULL) {
	    PRINTF("%-10s ++ %d ++:", "plus.key", ii); pmi(key, all, level + 1);
	}  else PRINTF("%-10s ++ %d ++: %s\n","plus.key", ii, "empty");	
      }
    }
  }

  if ((!givenkey && all==0) || all>0) {
    for (i=0; i<C->maxsub; i++) {
      if (cov->sub[i] == NULL) {
	//      PRINTF(" NULL\n");
	continue;
      }
      leer(level); 
      PRINTF("%s %d (%s) of '%s':", "submodel", i, C->subnames[i], C->nick);  
      pmi(cov->sub[i], all, level + 1);
    }
  }
}



void iexplDollar(cov_model *cov, bool MLEnatsc_only) {
  /*    
	get the naturalscaling values and devide the preceeding scale model     
	by this value
  */
  double *p, invscale;
  cov_model *dollar = cov->calling;

  bool solve = (cov->nr == NATSC_INTERN ||
		(cov->nr == NATSC_USER && !MLEnatsc_only))
    && dollar != NULL && isDollar(dollar);

  if (solve) {
    cov_model 
      *next = cov->sub[0];
    assert(dollar!=NULL && isDollar(dollar));

    INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
    if (ISNA(invscale)) error("inverse function of in 'iexplDollar' unknown");
    
    p = PARAM(dollar, DSCALE);
    if (p != NULL) {
      p[0] /= invscale;
    } else {
      p = PARAM(dollar, DANISO);      
      if (p != NULL) { 	
	int i,
	  n = dollar->nrow[DANISO] * dollar->ncol[DANISO];
	for (i=0; i<n; i++) p[i] *= invscale;
      } else {
	assert(cov->nr == NATSC_USER);
      }
    }
  } else {
    int i;
    for (i=0; i<MAXSUB; i++) { // cov->sub[i]: luecken erlaubt bei PSgen !
      if (cov->sub[i] != NULL) iexplDollar(cov->sub[i], MLEnatsc_only);
    }
  }
}


  
SEXP IGetModel(cov_model *cov, int modus, bool spConform, bool do_notreturnparam) {
  // modus:
  //  AS_SAVED : Modell wie gespeichert
  //  DEL_NATSC : Modell unter Annahme PracticalRange>0 (natsc werden geloescht)
  //  SOLVE_NATSC : natscale soweit wie moeglich zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)
  //  DEL_MLE : nur natscale_MLE werden geloescht
  //  SOLVE_MLE : nur natscale_MLE  zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)
  //
  // modus: 10-12 : wie 0-2, jedoch ohne CALL_FCT zu loeschen 


  SEXP model, nameMvec;
  int i, nmodelinfo,
    k = 0; 
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  bool plus_mixed_models;
  
  if ((cov->nr == NATSC_INTERN && cov->nr != GETMODEL_AS_SAVED) ||
      (cov->nr == NATSC_USER && modus == GETMODEL_DEL_NATSC)) { 
    return IGetModel(cov->sub[0], modus, spConform, do_notreturnparam);
  }
  
  nmodelinfo = C->kappas + 1;
  for (i=0; i<MAXSUB; i++) if (cov->sub[i] != NULL) nmodelinfo++;
  for (i=0; i<C->kappas; i++) {
    if (PisNULL(i) ||
	(do_notreturnparam && C->paramtype(i, 0, 0) == DONOTRETURNPARAM))
      nmodelinfo--;
  }

  PROTECT(model = allocVector(VECSXP, nmodelinfo));
  PROTECT(nameMvec = allocVector(STRSXP, nmodelinfo));

  SET_STRING_ELT(nameMvec, k, mkChar("")); // name
  cov_fct *CC = CovList + cov->nr; // nicht gatternr
  while(strncmp(CC->name, InternalName, strlen(InternalName)) ==0) CC--;
  //  print("conform %d\n", GLOBAL.general.sp_conform);
  if ((plus_mixed_models = (cov->nr == PLUS))) {
    plus_mixed_models = cov->calling == NULL; // to do: stimmt nicht mehr
    if (plus_mixed_models) {
      for (i=0; i<MAXSUB; i++) {
	if ((plus_mixed_models = 
	     (cov->sub[i] != NULL && cov->sub[i]->nr == MIXEDEFFECT))) break;
      }
    }
  }

  if (spConform && !plus_mixed_models) {
    SET_VECTOR_ELT(model, k++, mkString(CC->nick));
  } else {
    SET_VECTOR_ELT(model, k++, mkString(CC->name));
  }

  for(i=0; i<C->kappas; i++) {
    // naechste Zeile geloescht 30.4.2013
    if (PisNULL(i) ||
	 (do_notreturnparam && C->paramtype(i, 0, 0) == DONOTRETURNPARAM)) {	
       // k++; // 9.1.09, needed for "$" (proj can be NULL) -- 19.1.09 k++ 
	// deleted since otherwise outcome does not fit 
	// model definition anymore (arbitrary NULL in the definition)
	// instead: nmodelinfo--; above wihtin loop
      continue;
    }     
    SET_STRING_ELT(nameMvec, k, mkChar(C->kappanames[i]));
    SET_VECTOR_ELT(model, k++, 
		   Param((void*) cov->px[i], cov->nrow[i], cov->ncol[i], 
			 C->kappatype[i], true));
  }

  // vielleicht mal spaeter !
//  SET_STRING_ELT(nameMvec, k, mkChar("user"));  
//  SET_VECTOR_ELT(model, k++, Int(cov->user, Nothing + 1));

//  print("a %d %d %d %s\n",cov->nsub, k, nmodelinfo, NICK(cov));

  int zaehler = 0;
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
// print("i=%d %d %d %s\n", i,cov->nsub,zaehler,CovList[cov->sub[i]->nr].name);
      SET_STRING_ELT(nameMvec, k, mkChar(C->subnames[i]));
     SET_VECTOR_ELT(model, k++,
		    IGetModel(cov->sub[i], modus, spConform, do_notreturnparam));
      if (++zaehler >= cov->nsub) break;
    }
  }

  assert(k == nmodelinfo);

  setAttrib(model, R_NamesSymbol, nameMvec);
  UNPROTECT(2); // model + namemodelvec

  return model;
}

SEXP GetModel(SEXP keynr, SEXP Modus, SEXP SpConform, SEXP Do_notreturnparam) {
  // modus:
  //  AS_SAVED : Modell wie gespeichert
  //  DEL_NATSC : Modell unter Annahme PracticalRange>0 (natsc werden geloescht)
  //  SOLVE_NATSC : natscale soweit wie moeglich zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)
  //  DEL_MLE : nur natscale_MLE werden geloescht
  //  SOLVE_MLE : nur natscale_MLE  zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)
  //
  // modus: 10+ : wie oben, jedoch ohne CALL_FCT zu loeschen 


// Nutzer kann 3 Modifikationen des Models in MLE laufen lassen:
//      * keinen Praktikal range oder individuell angeben
//      * extern practical range definieren
//      * intern practical range verwenden lassen (use.naturalscaling)

  int knr = INTEGER(keynr)[0],
    modus = INTEGER(Modus)[0] % 10;
  bool
    delete_call = INTEGER(Modus)[0] < 10,
    spConform = (bool) INTEGER(SpConform)[0],
    do_notreturnparam = (bool) INTEGER(Do_notreturnparam)[0];

  //printf("modus=%d %d %d %d\n", modus, delete_call, spConform, do_notreturnparam);
  cov_model *cov;
  if (knr < 0 || knr  > MODEL_MAX || KEY[knr] == NULL) 
    return allocVector(VECSXP, 0);

  cov = KEY[knr];
  if (delete_call && isInterface(cov)) {    
    cov = cov->key == NULL ? cov->sub[0] : cov->key;
  }
  //PMI(cov);

  if (modus == GETMODEL_DEL_NATSC || modus == GETMODEL_DEL_MLE) {
     return IGetModel(cov, modus, spConform, do_notreturnparam);
  } else {
    cov_model *dummy = NULL; //ACHTUNG: "=NULL" hinzugefuegt
    SEXP value;
    if (covcpy(&dummy, cov) != NOERROR) return R_NilValue;
    iexplDollar(dummy, modus == GETMODEL_SOLVE_MLE);
    if (modus == GETMODEL_SOLVE_NATSC) {
      modus = GETMODEL_DEL_NATSC;
    } else if (modus == GETMODEL_SOLVE_MLE) {
      modus = GETMODEL_DEL_MLE;
    }
    value = IGetModel(dummy, modus, spConform, do_notreturnparam);
    COV_DELETE(&dummy);
    return(value);
  }
}




void Path(cov_model *cov, cov_model *sub) {
  //  printf("%ld\n", cov);
  //  printf("covnr =%d %s\n", cov->nr, NICK(cov));
  cov_fct *C = CovList + cov->nr;
  if (cov->calling == NULL) {
    PRINTF(" *** "); 
  } else {
    //    printf("calling\n");
    Path(cov->calling, cov);
  }

  if (sub == NULL) return; 
  if (cov->key == sub) { PRINTF("%s.key.%d->", C->nick, cov->zaehler); return; }

  int i;
  for (i=0; i<C->maxsub; i++) {
    if (cov->sub[i] == sub) {
      PRINTF("%s.sub[%d].%d->", C->nick, i, cov->zaehler);
      return;
    }
  }
  
  if (cov->Splus != NULL) {
    for (i=0; i<C->maxsub; i++) {
      if (cov->Splus->keys[i] == sub) {
	PRINTF("%s.S[%d].zaehler->", C->nick, i, cov->zaehler);
	return;
      }
    }
  }

  for (i=0; i<C->kappas; i++) {
    if (cov->kappasub[i] == sub) {
      PRINTF("%s.%s.%d->", C->nick, C->kappanames[i], cov->zaehler);
      return;
    }
  }


  PRINTF("%s (UNKNOWN,%d)->", C->nick, cov->zaehler);
}

void pmi(cov_model *cov) { // OK
  PRINTF("\n");

  if (cov == NULL) {
    PRINTF("\nCovariance model is empty.\n\n");
  } else {
    Path(cov, NULL);
    pmi(cov, false, 0);
  }
  //printf("pmi done\n");
}


void pmi(cov_model *cov, const char *msg) { 
  PRINTF("\n%s", msg); pmi(cov);
}



void pmi(const char *msg) { 
  cov_model *cov;
  int i;
  PRINTF("\n\n%s\n", msg);
  for (i=0; i<=MODEL_MAX; i++) {
    cov = KEY[i];
    if (cov != NULL) {
      PRINTF("Register '%s'\n", REGNAMES[i]);
      pmi(cov);
    }
  }
}


void pmi(cov_model *cov, char all) { 
  PRINTF("\n");
  if (cov == NULL) {
    PRINTF("\nCovariance model is empty.\n\n");
  } else {
    Path(cov, NULL);
    pmi(cov, all, 0);
  }
  //printf("pmi done\n"); 
}


/*
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

 


void PSTOR(cov_model *cov, storage *x) {  
  
  assert(cov != NULL);
  
  int d,
    dim = cov->tsdim;

  if (x==NULL) { PRINTF("no storage information\n"); return; }

 
  for (d=0; d<dim; d++) {
    //   PRINTF("%d. win:[%3.3f, %3.3f] c=%3.3f info:[%3.3f, %3.3f] E=%3.3f cum=%3.3f\n",
   PRINTF("%d. c=%3.3f info:[%3.3f, %3.3f] E=%3.3f cum=%3.3f\n",
	   d, //x->window.min[d], x->window.max[d], x->window.centre[d],
	   RF_NAN, RF_NAN, // pgs->mppinfo.min[d], x->mppinfo.max[d],
	   x->spec.E[d], x->spec.sub_sd_cum[d]);  
  }

  PRINTF("spec:step=%3.3f phi=%3.3f id=%3.3f grid=%s ergo=%s sig=%3.3f dens=%3.3f nmetr=%d\n",       
	 x->Sspectral.phistep2d, x->Sspectral.phi2d, x->Sspectral.prop_factor,
	 FT[x->Sspectral.grid], FT[x->Sspectral.ergodic], 
	 x->spec.sigma, x->spec.density, x->spec.nmetro);
}
