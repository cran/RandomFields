/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- get key strukture

 Copyright (C) 2001 -- 2017 Martin Schlather, 

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

//#include <Rmath.h>  
#include <stdio.h>  
//#include <stdlib.h>
//#include <sys/timeb.h>
//#include <unistd.h>
//#include <R_ext/Utils.h>     
#include <string.h>
#include <Rdefines.h>
#include "RF.h"
#include "Coordinate_systems.h"



#define WHICH_USER 0
#define WHICH_INTERNAL 1
#define WHICH_CALL_USER 2
#define WHICH_CALL_INTERNAL 3
#define WHICH_USER_BUT_ONCE 4
#define WHICH_INTERNAL_BUT_ONCE 5
#define WHICH_USER_BUT_ONCE_JUMP 6
#define WHICH_INTERNAL_BUT_ONCE_JUMP 7
#define WHICH_ALL_USER 8
#define WHICH_ALL 3

cov_model *WhichSub(cov_model *cov, int which){    
  // printf("%s %d %d\n", NAME(cov), which, isInterface(cov));
  if (!isInterface(cov)) return cov;
  cov_model *ans = cov;

  bool internfst = which == WHICH_INTERNAL || which == WHICH_USER_BUT_ONCE || 
    which == WHICH_USER_BUT_ONCE_JUMP,
    userfst = which == WHICH_USER || which == WHICH_INTERNAL_BUT_ONCE || 
    which==WHICH_INTERNAL_BUT_ONCE_JUMP;
  if (userfst || internfst){
    if (cov->Splus != NULL) warning("for '+', it is unclear which path to take");
    ans = internfst && ans->key != NULL ? ans->key : ans->sub[0];
    if (ans == NULL) BUG;
    if (which == WHICH_USER_BUT_ONCE_JUMP) ans = ans->sub[0];
    else if (which == WHICH_INTERNAL_BUT_ONCE_JUMP) ans = ans->key;
    if (ans == NULL) BUG;
    return ans;
  }
// which == WHICH_ALL_USER
  return cov;
}



name_type FT = {"false", "true"},
  TriNames = {"mismatch", 
	      "dep. next model",
	      "dep. prev. model",
	      "param. dependent",
	      "false", "true", 
	      "normal mixture",
	      "NaN"};  



SEXP GetLocationUserInfo(location_type **loc) {  
  int len;
  if (loc == NULL || (len = (*loc)->len) <= 0) return allocVector(VECSXP, 0);
   
  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, len));
  for (int i=0; i<len; i++) {
    SEXP namevec, l;
    location_type *L = loc[i];
    int k=0,
      elmts = L->Time;
    if (L->distances) {
      int laenge = L->lx * (L->lx - 1) / 2;
      elmts += 2;
      PROTECT(namevec = allocVector(STRSXP, elmts));
      PROTECT(l = allocVector(VECSXP, elmts));

      SET_STRING_ELT(namevec, k, mkChar("distances"));
      SET_VECTOR_ELT(l, k++, L->xdimOZ == 1 ? Num(L->x, laenge) 
		     : Mat(L->x, L->xdimOZ, laenge));

      SET_STRING_ELT(namevec, k, mkChar("dim"));
      SET_VECTOR_ELT(l, k++, ScalarInteger(L->timespacedim));  
    } else {
      if (L->ly > 0) elmts++;
      elmts += 2;
      PROTECT(namevec = allocVector(STRSXP,  elmts));
      PROTECT(l = allocVector(VECSXP, elmts));
      
      SET_STRING_ELT(namevec, k, mkChar("x"));
      SET_VECTOR_ELT(l, k++, L->grid 
		     ? Mat(L->xgr[0], 3, L->spatialdim)
		     : Mat_t(L->x, L->lx, L->xdimOZ));

      if (L->ly > 0) {
	SET_STRING_ELT(namevec, k, mkChar("y"));
	SET_VECTOR_ELT(l, k++, L->grid 
		     ? Mat(L->ygr[0], 3, L->spatialdim)
		     : Mat_t(L->y, L->ly, L->xdimOZ));
      }
  
      SET_STRING_ELT(namevec, k, mkChar("grid"));
      SET_VECTOR_ELT(l, k++, ScalarLogical(L->grid));  
    }    

    if (L->Time) {
      SET_STRING_ELT(namevec, k, mkChar("T"));
      SET_VECTOR_ELT(l, k++, Num(L->T, 3));       
    }
    assert(k == elmts); 

    setAttrib(l, R_NamesSymbol, namevec);
    SET_VECTOR_ELT(ans, i, l);
    UNPROTECT(2);
  }
  
  UNPROTECT(1); // l + namelvec

  return ans;
  
}






#define nlocinfo 13
SEXP GetLocationInfo(location_type *loc) {
  if (loc == NULL) return allocVector(VECSXP, 0);
  const char *info[nlocinfo] = 
    {"timespacedim", "xdimOZ", "spatialdim", "spatialtotpts", 
     "totpts", "distances", "grid", "Zeit", // Time geht nicht, da 
     //                                Abkuerzungskollision mit "T"
     "xgr", "x", "T", "ygr", "y"};
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
  SET_VECTOR_ELT(l, k++, ScalarInteger(loc->spatialdim));
  SET_VECTOR_ELT(l, k++, ScalarInteger((int) loc->spatialtotalpoints));
  SET_VECTOR_ELT(l, k++, ScalarInteger((int) loc->totalpoints));
  SET_VECTOR_ELT(l, k++, ScalarLogical(loc->distances));
  SET_VECTOR_ELT(l, k++, ScalarLogical(loc->grid));
  SET_VECTOR_ELT(l, k++, ScalarLogical(loc->Time));
  SET_VECTOR_ELT(l, k++, Mat(loc->xgr[0], loc->grid ? 3 : 0, spdim)); // "xgr"

  SET_VECTOR_ELT(l, k++, 
    Mat(loc->x, loc->xdimOZ, 
        loc->grid ? 0 : loc->distances ? loc->lx * (loc->lx - 1) / 2 : loc->lx
    )); //"x"
  SET_VECTOR_ELT(l, k++, Num(loc->T, loc->Time ? 3 : 0));// "T"
  if (loc->ly > 0) {    
    if (loc->distances) BUG;
    SET_VECTOR_ELT(l, k++, Mat(loc->ygr[0], loc->grid ? 3 : 0, spdim));
    SET_VECTOR_ELT(l, k++, Mat(loc->y, loc->xdimOZ, loc->grid ? 0 : loc->ly));
  } else {
    if (loc->ygr[0] != NULL || loc->y != NULL) BUG;
  }
  
  setAttrib(l, R_NamesSymbol, namevec);
  UNPROTECT(2); // l + namelvec
  assert(k == nloc);
  return l;
}


SEXP Param(cov_model *cov, void* p, int nrow, int ncol, SEXPTYPE type,
	     bool drop) {
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
    case STRSXP :
      return String((char*) p);
    case CLOSXP :
      BUG;
    case VECSXP : 
      if (cov->nr == COVARIATE) {      
       return(GetLocationUserInfo(cov->Scovariate->loc));
      } else {
	const char *msg[1] = {"R list"};
	return Char(msg, 1);
      }
      break;       
    case LANGSXP : 
      {
	const char *msg[1] = {"R object"};
	return Char(msg, 1);
      }
      break; 
    case ENVSXP : 
      {
	const char *msg[1] = {"R environment"};
	return Char(msg, 1);
      }
      break; 
    default:
      if (type >= LISTOF){
	SEXP dummy = NULL;	
	PROTECT(dummy = allocVector(VECSXP, nrow));
	listoftype *q = (listoftype *) p;
	for (i=0; i<nrow; i++) {
	  SET_VECTOR_ELT(dummy, i, 
			 Param(cov, q->lpx[i], q->nrow[i], q->ncol[i], REALSXP,
			       false));
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


SEXP GetModelInfo(cov_model *cov, int prlevel, bool both, int spConform, 
		  int whichSub, int Level) {
  assert(whichSub >=WHICH_USER && whichSub <= WHICH_ALL);
  // whichSub:  0=submodels, 1=keys, 2=both

  #include "primitive.h"
 
#define ninfobase 3
#define ninfo0 0
#define ninfo1 6
#define ninfo2 8
#define ninfo3 7
#define ninfo4 4

  /*   !!!!!     ACHTUNG     !!!!!
       Wenn diese Funktion geaendert wird,
       muss auch GetExtModelInfo geaendert werden
       !!!!!                 !!!!!
  */

  if (cov == NULL) return allocVector(VECSXP, 0);
  SEXP model, submodels, nameMvec, param, pnames;
  int i, j, nmodelinfo, 
    notnull = 0,
    k = 0;
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  location_type
    *loc = cov->calling == NULL ? PrevLoc(cov) : Loc(cov);
  bool return_param,    
    given_key = cov->Splus != NULL || cov->key != NULL,
    return_key = given_key && whichSub != WHICH_USER,
    return_sub = cov->nsub > 0 && (whichSub != WHICH_INTERNAL || !given_key),
    show_param = prlevel >= 1;
 
  //  printf("%s %d %d\n", NAME(cov), return_key, return_sub);


  for (i=0; i<C->kappas; i++) {
    if (both) {
      notnull += cov->nrow[i]>0 && cov->ncol[i]>0 && show_param;
      notnull += cov->kappasub[i] != NULL;
    } else {
      notnull += (cov->nrow[i]>0 && cov->ncol[i]>0 && show_param) ||
      cov->kappasub[i] != NULL;
    }
  }
  return_param = notnull > 0;

  nmodelinfo = ninfobase + (int) return_param;
  switch(prlevel > 5 ? 5 : prlevel) {
  case 5 : nmodelinfo += ninfo4;
  case 4 : nmodelinfo += ninfo3;
  case 3 : nmodelinfo += ninfo2;
  case 2 : nmodelinfo += ninfo1;
  case 1 : nmodelinfo += ninfo0;
  default: {}
  }
  if (!return_sub) nmodelinfo--;
  if (!return_key) nmodelinfo--; 

  PROTECT(model = allocVector(VECSXP, nmodelinfo));
  PROTECT(nameMvec = allocVector(STRSXP, nmodelinfo));

  SET_STRING_ELT(nameMvec, k, mkChar("name"));
  cov_fct *CC = CovList + cov->nr; // nicht gatternr579
  while(STRNCMP(CC->name, InternalName, STRLEN(InternalName)) ==0) CC--;

  if (spConform) {
    char name[MAXCHAR+2];
    SPRINTF(name, "%s", CC->nick); 
    SET_VECTOR_ELT(model, k++, mkString(name));
  } else {
    SET_VECTOR_ELT(model, k++, mkString(CC->name));
  }

  if (return_param) {
    PROTECT(param = allocVector(VECSXP, notnull));
    PROTECT(pnames = allocVector(STRSXP, notnull));
    for (j=i=0; i<C->kappas; i++) {

      if (cov->kappasub[i] != NULL) {	  
	  SET_STRING_ELT(pnames, j, 
			 mkChar(!STRCMP(C->kappanames[i], FREEVARIABLE)
				&& cov->ownkappanames[i] != NULL
				? cov->ownkappanames[i] 
				: C->kappanames[i])); 
	  SET_VECTOR_ELT(param, j,
			 GetModelInfo(cov->kappasub[i], prlevel, both, 
				      spConform, whichSub, Level + 1)); 
	  j++;
	  if (!both) continue;
      }
 
      if (show_param && cov->nrow[i]>0 && cov->ncol[i]>0) {
	if (isAnyDollar(cov) && i==DANISO) {
	  SET_STRING_ELT(pnames, j, mkChar("Aniso"));
	  double
	    *Aniso = (double*) MALLOC(cov->nrow[i] * cov->ncol[i] * 
				      sizeof(double));
	  int t,l,m;
	  
	  for (t=l=0; l<cov->nrow[i]; l++) {
	    for (m=0; m<cov->ncol[i]; m++) {
	      Aniso[t++] = P(DANISO)[m * cov->nrow[i] + l];
	    }
	  }
	  SET_VECTOR_ELT(param, j,
			 Param(cov, (void*) Aniso, cov->nrow[i],
			       cov->ncol[i], C->kappatype[i], true));    
	  FREE(Aniso);
	  j++;
	  continue;
	}
	
	SET_STRING_ELT(pnames, j, 
		       mkChar(!STRCMP(C->kappanames[i], FREEVARIABLE)
			      && cov->ownkappanames[i] != NULL
			      ? cov->ownkappanames[i] 
			      : C->kappanames[i]));
	SET_VECTOR_ELT(param, j,
		       Param(cov, (void*) cov->px[i], cov->nrow[i], 
			     cov->ncol[i], C->kappatype[i], true));    
	j++;
      }
    }
    setAttrib(param, R_NamesSymbol, pnames);
    
    SET_STRING_ELT(nameMvec, k, mkChar("param"));  
    SET_VECTOR_ELT(model, k++, param);
    UNPROTECT(2);

  }
  
  if (prlevel>=2) {      
    SET_STRING_ELT(nameMvec, k, mkChar("covnr"));
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->nr));

    SET_STRING_ELT(nameMvec, k, mkChar("vdim"));
    SET_VECTOR_ELT(model, k++,  ScalarInteger(cov->vdim[0])
		   // cov->vdim[0] == cov->vdim[1] 
		   // ? ScalarInteger(cov->vdim[0]) 
		   // : Int(cov->vdim, 2)
		   );
 
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
   
  if (prlevel>=3) {
    SET_STRING_ELT(nameMvec, k, mkChar("type"));  
    SET_VECTOR_ELT(model, k++, Char(TYPENAMES + cov->typus, 1));

    SET_STRING_ELT(nameMvec, k, mkChar("role"));  
    SET_VECTOR_ELT(model, k++, Char(ROLENAMES + cov->role, 1));
 
    SET_STRING_ELT(nameMvec, k, mkChar("domown"));  
    SET_VECTOR_ELT(model, k++, mkString(DOMAIN_NAMES[cov->domown]));
  
    SET_STRING_ELT(nameMvec, k, mkChar("isoown"));
    SET_VECTOR_ELT(model, k++, mkString(ISONAMES[cov->isoown]));
  
    if (prlevel >= 5) {
    SET_STRING_ELT(nameMvec, k, mkChar("isoprev"));
      SET_VECTOR_ELT(model, k++, mkString(ISONAMES[cov->isoprev]));
    }
   
 
    SET_STRING_ELT(nameMvec, k, mkChar("internalq"));  
    SET_VECTOR_ELT(model, k++, Num(cov->q, cov->qlen));
  
    SET_STRING_ELT(nameMvec, k, mkChar("pref"));  
    SET_VECTOR_ELT(model, k++, Int(cov->pref, Nothing + 1));
   
    SET_STRING_ELT(nameMvec, k, mkChar("simu"));  
    SET_VECTOR_ELT(model, k++, GetSimuInfo(&(cov->simu)));
  	
    SET_STRING_ELT(nameMvec, k, mkChar("loc"));  
    SET_VECTOR_ELT(model, k++, 
		   Level == 0 || cov->calling == NULL ||
		   Loc(cov->calling) != Loc(cov) 
		   ? GetLocationInfo(loc) 
		   : Loc(cov) == NULL 
		   ? mkString("no locations given")
		   : mkString("same as calling model"));
  }

  if (prlevel>=4) {
    SET_STRING_ELT(nameMvec, k, mkChar("logspeed"));
    SET_VECTOR_ELT(model, k++, ScalarReal(cov->logspeed));

    SET_STRING_ELT(nameMvec, k, mkChar("maxdim"));  
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->maxdim));
      
    SET_STRING_ELT(nameMvec, k, mkChar("full_derivs"));  
    SET_VECTOR_ELT(model, k++, ScalarInteger(cov->full_derivs));
      
    SET_STRING_ELT(nameMvec, k, mkChar("loggiven"));  
    SET_VECTOR_ELT(model, k++, ScalarLogical(cov->loggiven));

    SET_STRING_ELT(nameMvec, k, mkChar("monotone"));  

    int idx = cov->monotone - (int) MISMATCH;
    if (idx < 0 || idx > BERNSTEIN - (int) MISMATCH) {
      PRINTF("monotone %d %d %d\n",  cov->monotone,  MISMATCH,
	     cov->monotone - (int) MISMATCH);
      BUG;
    }
    SET_VECTOR_ELT(model, k++, Char(MONOTONE_NAMES + idx, 1));
            
    SET_STRING_ELT(nameMvec, k, mkChar("MLE"));  
    SET_VECTOR_ELT(model, k++, ScalarLogical(cov->MLE != NULL));

    SET_STRING_ELT(nameMvec, k, mkChar("finiterange"));  
    SET_VECTOR_ELT(model, k++, ScalarLogical(cov->finiterange));
      
    //  SET_STRING_ELT(nameMvec, k, mkChar("idx"));  
    //  SET_VECTOR_ELT(model, k++, Int(cov->idx, cov->tsdim));
    //      SET_STRING_ELT(nameMvec, k, mkChar("user"));  
    // SET_VECTOR_ELT(model, k++, Int(cov->user, Nothing + 1));
  }

  if (prlevel>=5) {

    mpp_properties *Mpp = &(cov->mpp);
    int
      vdim = cov->vdim[0],
      nm = cov->mpp.moments,
      nmvdim = (nm + 1) * vdim;
 
    //   SET_STRING_ELT(nameMvec, k, mkChar("mpp.refradius"));  
    //  SET_VECTOR_ELT(model, k++, ScalarReal(mpp->refradius));
 
    SET_STRING_ELT(nameMvec, k, mkChar("mpp.maxheight"));  
    SET_VECTOR_ELT(model, k++, Num(Mpp->maxheights, vdim, MAXINT));
  
    SET_STRING_ELT(nameMvec, k, mkChar("mpp.M"));
    SET_VECTOR_ELT(model, k++, cov->mpp.moments == 0 ? R_NilValue :
		   Num(Mpp->mM, nmvdim, MAXINT)); 

    SET_STRING_ELT(nameMvec, k, mkChar("mpp.Mplus"));
    SET_VECTOR_ELT(model, k++, cov->mpp.moments == 0 ? R_NilValue :
		   Num(Mpp->mMplus, nmvdim, MAXINT));
  } 
   		 
  if (return_key) {
    if (cov->key != NULL) {
      SET_STRING_ELT(nameMvec, k, mkChar("internal"));  
      SET_VECTOR_ELT(model, k++, GetModelInfo(cov->key, prlevel, both, 
					      spConform, whichSub, Level + 1));
    } else {  /// cov->nr == PLUS && cov->Splus != NULL
      int ii, n,
	subs = C->maxsub;
      SEXP keys;
      SET_STRING_ELT(nameMvec, k, mkChar("internal"));  
      for (ii=n=0; ii<subs; ii++) if (cov->Splus->keys[ii] != NULL) n++;
      PROTECT(keys = allocVector(VECSXP, n));
      for (ii=n=0; ii<subs; ii++) 
	if (cov->Splus->keys[ii] != NULL)
	  SET_VECTOR_ELT(keys, n++, 
			 GetModelInfo(cov->Splus->keys[ii], prlevel, both,
				      spConform, whichSub, Level + 1));
      SET_VECTOR_ELT(model, k++, keys);
      UNPROTECT(1);
    }
  }

  if (return_sub) {
    SET_STRING_ELT(nameMvec, k, mkChar("submodels"));  
    PROTECT(submodels = allocVector(VECSXP, cov->nsub));
    int zaehler = 0;
    for (i=0; i<MAXSUB; i++) {
      if (cov->sub[i] != NULL) {
	SET_VECTOR_ELT(submodels, zaehler, 
		       GetModelInfo(cov->sub[i], prlevel, both, spConform, 
				    whichSub, Level + 1));
	if (++zaehler >= cov->nsub) break;
      }
    }
    SET_VECTOR_ELT(model, k++, submodels);
    UNPROTECT(1);
  }
  
  setAttrib(model, R_NamesSymbol, nameMvec);
  assert(k==nmodelinfo); 
  UNPROTECT(2); // model + namemodelvec
  return model;
}

SEXP GetExtModelInfo(SEXP keynr, SEXP Prlevel, SEXP spConform, SEXP whichSub) {
  int knr = INTEGER(keynr)[0],
    which = INTEGER(whichSub)[0] == WHICH_ALL_USER 
    ? WHICH_ALL : INTEGER(whichSub)[0] % 2,
    prlevel = std::abs(INTEGER(Prlevel)[0]) % 10;
  bool 
    both = INTEGER(Prlevel)[0] < 0,
    delete_call = std::abs(INTEGER(Prlevel)[0]) < 10;
  cov_model *cov, *orig;
  SEXP res, names;

  if (knr>=0 && knr <= MODEL_MAX && KEY[knr] != NULL) {
    orig = KEY[knr];
    cov = WhichSub(orig, INTEGER(whichSub)[0]);
    res = GetModelInfo(cov, prlevel, both, 
		       (bool) INTEGER(spConform)[0], which, 0);
    if (prlevel>=1 && delete_call) {
      names = getAttrib(res, R_NamesSymbol);
      int i, len = length(names);
      for (i=0; i<len; i++) {
	const char *name = CHAR(STRING_ELT(names, i));
	if (STRCMP("xdimprev", name) == 0) {
	  INTEGER(VECTOR_ELT(res, i))[0] = orig->xdimprev;
	  break;
	}
      }
    }
    return res;
  }
  return allocVector(VECSXP, 0);
}


void Path(cov_model *cov, cov_model *sub) {
  cov_fct *C = CovList + cov->nr;
  const char *sep="-> ";
  if (cov->calling == NULL) {
    PRINTF(" *** "); 
  } else {
    Path(cov->calling, cov);
  }
 
  if (sub == NULL) return; 
  if (cov->key == sub) { PRINTF("%s.key.%d%s", C->nick, cov->zaehler, sep); return; }

  int i;
  for (i=0; i<C->maxsub; i++) {
    if (cov->sub[i] == sub) {
      PRINTF("%s[%s,%d].%d%s", C->nick, C->subnames[i], i, cov->zaehler, sep);
      return;
    }
  }
  
  if (cov->Splus != NULL) {
    for (i=0; i<C->maxsub; i++) {
      if (cov->Splus->keys[i] == sub) {
	PRINTF("%s.S[%d].%d%s", C->nick, i, cov->zaehler, sep);
	return;
      }
    }
  }

  for (i=0; i<C->kappas; i++) {
    if (cov->kappasub[i] == sub) {
      PRINTF("%s.%s.%d%s", C->nick, C->kappanames[i], cov->zaehler, sep);
       return;
    }
  }


  PRINTF("%s (UNKNOWN,%d)%s", C->nick, cov->zaehler, sep);
}


void leer(int level){
  char format[255];
  if (level == 0) return;
  SPRINTF(format,"%%%ds", -level * 3);
  PRINTF(format, "");
}

int MAX_PMI = 5;

void PrintPoints(location_type *loc, char *name, 
		 double *x, coord_type xgr, long lx) {
#ifndef SHOW_ADDRESSES  
  return; // unclear error below on CRAN
#endif 

#define maxpts 100
  int i;
  if (loc->grid) { 
    PRINTF("loc:%sgr    ", name);
    for (i=0; i<loc->timespacedim - loc->Time; i++) {
      PRINTF("(%3.3f, %3.3f, %2.0f) ", xgr[i][XSTART], xgr[i][XSTEP], 
	     xgr[i][XLENGTH]);
    }
  } else {
    PRINTF("loc:%s      ", name);
    if (loc->lx <= 0) {
      PRINTF("not given! (%d)", addressbits(loc->x));
    } else {
      long total = loc->distances ? lx * (lx-1) / 2 : lx * loc->xdimOZ,
	endfor = total;
      if (endfor > maxpts) endfor = maxpts;
      for (i=0; i<endfor; i++) {
	PRINTF("%4.3f", x[i]);
	if ((i+1) % loc->xdimOZ == 0) PRINTF(";");
	PRINTF(" ");
      }
      if (endfor < total) 
	PRINTF("... [%ld not shown]", total - endfor);
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
  leer(level); PRINTF("%-10s %ld\n", "loc:lx", loc->lx);
  leer(level); PRINTF("%-10s %ld %ld\n","loc:totpts",
		      loc->spatialtotalpoints, loc->totalpoints);
  leer(level); PRINTF("%-10s %ld\n","loc:len", loc->len);
  leer(level); PRINTF("%-10s %s\n","loc:grid", FT[loc->grid]);
  leer(level); PRINTF("%-10s %s\n","loc:dist", FT[loc->distances]);
  leer(level); PRINTF("%-10s %s\n","loc:Time", FT[loc->Time]);
#ifdef SHOW_ADDRESSES    
  leer(level); PrintPoints(loc, (char *) "x", loc->x, loc->xgr, loc->lx);
  if (loc->y!=NULL || loc->ygr[0]!=NULL) {  
    assert(loc->ly > 0);
    leer(level); PrintPoints(loc, (char*) "y", loc->y, loc->ygr, loc->ly);
  } 
#else
  leer(level); PRINTF("loc:x,y\t addresses not shown\n");
#endif

  if (loc->Time) { 
    leer(level); PRINTF("%-10s (%f %f %f)\n", "loc:T", 
			loc->T[0], loc->T[1], loc->T[2]);
  }

  leer(level); PRINTF("%-10s ","loc:cansio");
  if (loc->caniso == NULL) {
    PRINTF("null\n"); 
  } else {
    int endfor = loc->cani_nrow * loc->cani_ncol;
    PRINTF(" [%d, %d] ",  loc->cani_nrow, loc->cani_ncol);
    if (endfor > MAX_PMI) endfor = MAX_PMI;
    for (i=0; i<endfor; i++) PRINTF(" %f", loc->caniso[i]); 
    PRINTF("\n");
  }
}
  

static bool PMI_print_dollar = !true,
  PMI_print_mpp = !true,
  PMI_print_pgs = !true, 
  PMI_print_details = !true,
  PMI_print_loc = true,
  PMI_print_fields = !true;
static int  PMI_print_rect = 1;  // 0, 1, 2

void pmi(cov_model *cov, char alles, int level, int maxlevel) {    
  int i, j, endfor;
  cov_fct *C = CovList + cov->nr; // nicht gatternr
#define MNlength 4
  char MN[Forbidden + 1][MNlength], name[100];
 
  //int n = 2;
  if (level > maxlevel) {
    leer(level); PRINTF("'%s' left out\n", NAME(cov)); return;
  }

  for (i=0; i<=Forbidden; i++) {
    strcopyN(MN[i], METHODNAMES[i], MNlength);
  }

  cov_fct *CC = C;
  while(STRCMP(CC->name, InternalName) ==0) CC--;
  if (level == 0) 
    PRINTF("******   %s, %s   ****** [%d,%d]", CC->nick, CC->name, cov->nr, cov->zaehler);
  else PRINTF("    **** %s, %s **** [%d,%d]", CC->nick, CC->name, cov->nr, cov->zaehler);  PRINTF("\n");

  leer(level); PRINTF("%-10s %s\n", "param", C->kappas == 0 ? "none" : ""); 

  for(i=0; i<C->kappas; i++) {
    STRCPY(name, 
	   !STRCMP(C->kappanames[i], FREEVARIABLE) && 
	           cov->ownkappanames != NULL && cov->ownkappanames[i] != NULL
	   ? cov->ownkappanames[i]
	   : C->kappanames[i]
	   );
    name[9] = '\0';
    leer(level + 1); PRINTF("%-10s", 
			    name); 
    if (PisNULL(i)) {
      PRINTF(" NULL");
    } else { // ! PisNULL
      switch(C->kappatype[i]) {
      case REALSXP : {
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
      }
	break;
      case INTSXP : {
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
      }
	break;
      
      case CLOSXP : PRINTF("< arbitrary function >"); break;
      case VECSXP : PRINTF("< some R list >"); break;
      case LANGSXP: PRINTF("< language expression >"); break;
      case ENVSXP : PRINTF("< R environment >"); break;
      case LISTOF + REALSXP : {
	int k, ende=0,
	  k_max = 20,
	  k_end = cov->nrow[i];
	listoftype *p= PLIST(i);	
	PRINTF("list [%d]\n", cov->nrow[i]);
	if (k_end > k_max) k_end = k_max;
	for (k=0; k<k_end; k++) {
	  leer(level + 2); 
   	  if (p->ncol[k]==1) {
            if (p->nrow[k]==1) {
   	      PRINTF("%f", p->lpx[k][0]); 
	    } else {
	      PRINTF("[%d] ", p->nrow[k]);
	      ende = endfor = p->nrow[k]; 
	      if (endfor > MAX_PMI) endfor = MAX_PMI;
	      for (j=0; j<endfor; j++) PRINTF(" %f", p->lpx[k][j]); 
            }
          } else {
            PRINTF("[%d, %d] ", p->nrow[k],  p->ncol[k]);
            ende = endfor = p->nrow[k] * p->ncol[k]; 
	    if (endfor > MAX_PMI) endfor = MAX_PMI;
	    for (j=0; j<endfor; j++) PRINTF(" %f", p->lpx[k][j]); 
          }
	  if (ende > MAX_PMI) PRINTF ("...");
	  if (k<cov->nrow[i]-1)  PRINTF("\n"); 
        }
	if (k_end < cov->nrow[i]) { leer(level + 2); PRINTF ("..."); }
      }
      break;
      default:  BUG;
      } // switch
    }
    if (cov->kappasub[i] != NULL) {
      PRINTF("  <=");
      pmi(cov->kappasub[i], alles, level + 3, maxlevel);
    } else  PRINTF("\n");
  } // for
  if (cov->Sset != NULL) {
    cov_model *from = cov->Sset->remote;
    leer(level + 1); 
    PRINTF("%-10s%s [%d]\n", "<remote>", Nick(from), from->zaehler);     
  }

  leer(level); PRINTF("%-10s [%d]","internal-q", cov->qlen);  

  endfor = cov->qlen; if (endfor > MAX_PMI) endfor = MAX_PMI;
  for (i=0; i<endfor; i++) PRINTF(" %f", cov->q[i]); 
  PRINTF("\n");

  if (cov->calling == NULL && level != 0) {
    PRINTF("current model is %s\n", NICK(cov)); 
    BUG;
  }
  leer(level); 
  if (cov->calling==NULL)  PRINTF("%-10s %s\n","calling", "NULL");
  else PRINTF("%-10s %s [%d]\n","calling", 
	      Nick(cov->calling), cov->calling->zaehler);  
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
		      DOMAIN_NAMES[(int) cov->domprev],
		      DOMAIN_NAMES[(int) cov->domown]);
  leer(level); PRINTF("%-10s %d/%d (%s/%s)\n","isoprev/own", 
		      cov->isoprev, cov->isoown,
		      ISONAMES[(int) cov->isoprev],
		      ISONAMES[(int) cov->isoown]);
  leer(level); PRINTF("%-10s %d, %d:%d:%d, %d/%d\n","ts-x-v-dim",
		      cov->tsdim, cov->xdimprev, cov->xdimgatter, cov->xdimown,
		      cov->vdim[0], cov->vdim[1]);  
   leer(level); PRINTF("%-10s %d\n","maxdim", cov->maxdim);  
  leer(level); PRINTF("%-10s %s (%d)\n", "type", TYPENAMES[cov->typus],
		      (int) cov->typus);
  leer(level); PRINTF("%-10s %s (%d)\n","role", ROLENAMES[cov->role],
		      (int) cov->role);
  leer(level); PRINTF("%-10s %s\n","determinst", FT[cov->deterministic]);  
  if (PMI_print_fields) {
    leer(level); PRINTF("%-10s %s\n","method", METHODNAMES[cov->method]);
    leer(level); PRINTF("%-10s %s\n","initialised", FT[cov->initialised]);
    leer(level); PRINTF("%-10s %s\n","fieldret", FT[cov->fieldreturn]);
    leer(level); PRINTF("%-10s %s\n","origrf", FT[cov->origrf]);
    leer(level); PRINTF("%-10s %d\n","rf", addressbits(cov->rf));
    leer(level); PRINTF("%-10s %s\n","checked", FT[cov->checked]);
    leer(level); PRINTF("%-10s ", "pref");  
    for (i=0; i<=Sequential; i++) PRINTF("%s:%d ", MN[i], (int) cov->pref[i]);
    PRINTF("\n"); leer(level); PRINTF("%-10s ", "");  
    for (; i<=Nothing; i++) PRINTF("%s:%d ", MN[i], (int) cov->pref[i]);
    PRINTF("\n");
  }
  if (PMI_print_details) {
    leer(level); PRINTF("%-10s %s\n","user_given", 
		      cov->user_given == ug_internal ? "internal" :
		      cov->user_given == ug_explicit ? "explicit" : "implicit");
    leer(level); PRINTF("%-10s %f\n","logspeed", cov->logspeed);  
    leer(level); PRINTF("%-10s %d, %d\n","full/rese deriv's", 
			cov->full_derivs, cov->rese_derivs);  
    leer(level); PRINTF("%-10s %s\n","loggiven", FT[cov->loggiven]);  
    int idx = cov->monotone - (int) MISMATCH;
    if (idx < 0 || idx > BERNSTEIN - (int) MISMATCH) {
      PRINTF("monotone %d %d %d\n",  cov->monotone,  MISMATCH,
	     cov->monotone - (int) MISMATCH);
      BUG;
    }
    leer(level); PRINTF("%-10s %s (%d)\n","monotone", 
			MONOTONE_NAMES[idx], cov->monotone);  
    leer(level); PRINTF("%-10s %s (%d)\n","pointwise", 
			POSITIVITY_NAMES[cov->ptwise_definite],
			cov->ptwise_definite);  
    leer(level); PRINTF("%-10s %s\n","finiterng", 
			TriNames[cov->finiterange - MISMATCH]);  
    leer(level); PRINTF("%-10s %d %d %d %d\n","stor/2/extra/like", 
			addressbits(cov->Sgen),
			addressbits(cov->Sgatter),
			addressbits(cov->Sextra),
			addressbits(cov->Slikelihood));
    leer(level); PRINTF("%-10s %s\n","simu:activ", FT[cov->simu.active]);
    leer(level); PRINTF("%-10s %s\n","simu:pair", FT[cov->simu.pair]);
    leer(level); PRINTF("%-10s %d\n","simu:expect", cov->simu.expected_number_simu);    
    leer(level); PRINTF("%-10s %s\n", "MLE", FT[cov->MLE==NULL]); 
    int d;
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
  //leer(level); PRINTF("%-10s %d\n","naturalscaling", cov->naturalscaling);  
  //leer(level); PRINTF("%-10s %d\n","init",CovList[cov->nr].init!=init_failed);
   //leer(level); PRINTF("%-10s %d\n","tbm2num", (int) cov->tbm2num);  
  //leer(level); PRINTF("%-10s %d\n","spec:nmetro", cov->spec.nmetro);
  //leer(level); PRINTF("%-10s %f\n","spec:sigma", cov->spec.sigma);
  if (PMI_print_dollar && cov->Sdollar != NULL) {
    leer(level); PRINTF("%-10s %d\n","$:z", addressbits(cov->Sdollar->z));
    leer(level); PRINTF("%-10s %d\n","$:z2", addressbits(cov->Sdollar->z2));
    leer(level); PRINTF("%-10s %d\n","$:y", addressbits(cov->Sdollar->y));
  }

  if (PMI_print_mpp) {
    int 
      nm = cov->mpp.moments,
      vdim = cov->vdim[0],
      nmvdim = (nm + 1) * vdim;
    if (R_FINITE(cov->mpp.maxheights[0])) {
      leer(level); PRINTF("%-10s ","mpp:maxhgt"); 
      for (i=0; i<=vdim; i++) PRINTF("%f, ",cov->mpp.maxheights[i]);
      PRINTF("\n");    
    }
    if (R_FINITE(cov->mpp.unnormedmass)) {
      leer(level); PRINTF("%-10s %f\n","mpp:u-mass", cov->mpp.unnormedmass);    
    }
    leer(level); PRINTF("%-10s ","mpp:M+");
    if (cov->mpp.mMplus == NULL) 
      PRINTF("not initialized yet (size=%d)\n", cov->mpp.moments);
    else {
      for (i=0; i<nmvdim; i++) PRINTF("%f, ", cov->mpp.mMplus[i]);
      PRINTF("\n");    
    }
    leer(level); PRINTF("%-10s ","mpp:M");
    if (cov->mpp.mM == NULL)
      PRINTF("not initialized yet (size=%d)\n", cov->mpp.moments);
    else {
      for (i=0; i<nmvdim; i++)  PRINTF("%f, ", cov->mpp.mM[i]);
      PRINTF("\n");    
    }
    leer(level); PRINTF("%-10s %s\n","mpp:log",
			FT[CovList[cov->nr].log != ErrLogCov]);
    leer(level); PRINTF("%-10s %s\n","mpplgnonst", 
			FT[CovList[cov->nr].nonstatlog != ErrLogCovNonstat]);
    leer(level); PRINTF("%-10s %s\n","mpp:do",
			FT[CovList[cov->nr].Do != do_failed]);
  }

  
  if (PMI_print_pgs) {
    if (cov->Spgs == NULL) {
    } else {
      pgs_storage *pgs = cov->Spgs;
      int d,
	size = pgs->size,
	dim = cov->xdimown;
      
      leer(level); PRINTF("%-10s %f\n","pgs:mass", pgs->totalmass);
      leer(level); PRINTF("%-10s %f\n","pgs:logdens", pgs->log_density);
      
      leer(level); PRINTF("%-10s %d\n","pgs:size", size);
      leer(level); PRINTF("%-10s %d\n","pgs:dim", dim);
      
#define SHOWDEFAULT(Z, X, Y) if (pgs->X != NULL) {			\
	leer(level);  PRINTF("%-10s ",Z);				\
	for (d=0; d<dim; d++) PRINTF(Y, pgs->X[d]); PRINTF("\n");}
#define SHOW(Z, X) SHOWDEFAULT(Z, X, "%f ")
#define SHOWINT(Z, X) SHOWDEFAULT(Z, X, "%d ")
      
      if (R_FINITE(pgs->zhou_c)) {
	leer(level); PRINTF("%-10s %f\n","mpp:zhou_c", pgs->zhou_c);
      } 
      
      SHOW("pgs:v", v);    
      SHOW("pgs:x", x);
      SHOW("pgs:own_start", own_grid_start);
      SHOW("pgs:own_step", own_grid_step);
      SHOW("pgs:xstart", xstart);
      SHOW("pgs:inc", inc);
      SHOW("pgs:suppmin", supportmin);
      SHOW("pgs:suppmax", supportmax);
      SHOWINT("pgs:gridlen", gridlen);
      SHOWINT("pgs:start", start);
      SHOWINT("pgs:end", end);
      SHOWINT("pgs:nx", nx);
      if (pgs->pos != NULL) { // gauss
#ifdef SHOW_ADDRESSES
	location_type *loc = Loc(cov);
	leer(level); PrintPoints(loc, (char *) "pgs.x", loc->x, pgs->xgr,loc->lx);
#endif
	SHOW("pgs:y", y);
	SHOWINT("pgs:pos", pos);
	SHOWINT("pgs:min", min);
	SHOWINT("pgs:max", max);
      }
      if (pgs->halfstepvector != NULL) { // max-stable 
	leer(level); PRINTF("%-10s %s\n","pgs:flat", FT[pgs->flat]);
	leer(level); PRINTF("%-10s %f\n","pgs:globmin", pgs->globalmin);
	leer(level); PRINTF("%-10s %f\n","pgs:cur.thres",pgs->currentthreshold);
	SHOW("pgs:half", halfstepvector);
	if (pgs->single != NULL) {
	  leer(level); { PRINTF("%-10s ","pgs:single"); 
	    for (d=0; d<size; d++) PRINTF("%f ", pgs->single[d]);PRINTF("\n"); }
	}
	if (pgs->total != NULL) {
	  leer(level); { PRINTF("%-10s ","pgs:total"); 
	    for (d=0; d<size; d++) PRINTF("%f ", pgs->total[d]); PRINTF("\n"); }
	}
      } else { // gauss oder poisson
	leer(level); PRINTF("%-10s %f\n","pgs:intens", pgs->intensity);
      }
      if (pgs->cov != NULL) {
	leer(level); 
	PRINTF("%-10s %s [%d]\n","pgs:cov", Nick(pgs->cov), pgs->cov->zaehler);
      }
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
      leer(level); PRINTF("%-10s %f\n","rct:total", p->weight[1 + p->nstep]);
      
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
	leer(level); { PRINTF("%-10s ","rct:asSign"); 
	  for (d=0; d<p->tmp_n; d++) PRINTF("%d ", p->asSign[d]);
	  PRINTF("\n"); }
	leer(level); { PRINTF("%-10s ","rct:idx"); 
	  for (d=0; d<dimP1; d++) PRINTF("%d ", p->i[d]); PRINTF("\n"); }
      }
    }
   }
    
  if (PMI_print_loc) {
    leer(level);  
    //    printf(" **** %ld %ld %ld\n", PLoc(cov), cov->ownloc, cov->prevloc);
   if (PLoc(cov) == NULL) {
      PRINTF("%-10s %s\n", "loc:", "<null>");
    } else {
     //printf(" **** %ld %ld %ld %d \n", PLoc(cov), cov->ownloc, cov->prevloc,  cov->prevloc[0]->len);
      PRINTF("%-10s %d\n", "loc:sets", Loc(cov)->len);
      if (cov->ownloc != NULL) {
	PrintLoc(level, OwnLoc(cov), true);
      } else {
	location_type *prevloc = PrevLoc(cov);
	leer(level);    
	if (cov->calling == NULL) { 
	   PRINTF("%-10s (%d)\n", "loc:extern", addressbits(prevloc));
	   PrintLoc(level, prevloc, false);
	} else if (prevloc == PrevLoc(cov->calling)) {
	  PRINTF("%-10s (%d)\n", "loc:calling->prev", addressbits(prevloc));
	  // PrintLoc(level, prevloc, false);
	} else if (prevloc == OwnLoc(cov->calling)) {
	  PRINTF("%-10s (%d)\n", "loc:calling->own", addressbits(prevloc));
	  // PrintLoc(level, prevloc, false);
	} else {
	  PRINTF("%-10s (%d)\n", "loc:MISMATCH", addressbits(prevloc));
	}
      }
    }
  }

  bool givenkey = false;
  if (cov->key != NULL) {
    leer(level);
    givenkey = true;
    PRINTF("%-10s :", "key");  
    if (alles >= 0) pmi(cov->key, alles, level + 1, maxlevel);
  }

  if (cov->Splus != NULL) {
    givenkey = true;
    int ii;
    if (alles >= 0) {
      for (ii=0; ii < cov->nsub; ii++) {
	cov_model *key = cov->Splus->keys[ii];
	leer(level);      
	if (key != NULL) {
	    PRINTF("%-10s ++ %d ++:", "plus.key", ii); 
	    pmi(key, alles, level + 1, maxlevel);
	}  else PRINTF("%-10s ++ %d ++: %s\n","plus.key", ii, "empty");	
      }
    }
  }

  if ((!givenkey && alles==0) || alles>0) {
    for (i=0; i<C->maxsub; i++) {
      if (cov->sub[i] == NULL) {
	continue;
      }
      leer(level); 
      PRINTF("%s %d (%s) of '%s':", "submodel", i, C->subnames[i], C->nick);  
      pmi(cov->sub[i], alles, level + 1, maxlevel);
    }
  }
}



void pmi(cov_model *cov, int maxlevel) { // OK
  PRINTF("\n");

  if (cov == NULL) {
    PRINTF("\nCovariance model is empty.\n\n");
  } else {
    Path(cov, NULL);
    pmi(cov, false, 0, maxlevel);
  }
}


void pmi(cov_model *cov) { // OK
  pmi(cov, 999999);
}


void iexplDollar(cov_model *cov, bool MLEnatsc_only) {
  /*    
	get the naturalscaling values and devide the preceeding scale model     
	by this value
  */
  double *p, invscale;
  cov_model *dollar = cov->calling;

 
  bool solving = (cov->nr == NATSC_INTERN ||
		(cov->nr == NATSC_USER && !MLEnatsc_only))
    && dollar != NULL && isDollar(dollar);

  if (solving) {
    cov_model 
      *next = cov->sub[0];
    assert(dollar!=NULL && isDollar(dollar));

    INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
    //invscale = 1.0;
    if (ISNAN(invscale)) ERR("inverse function of in 'iexplDollar' unknown");
    
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


  
SEXP IGetModel(cov_model *cov, int modus, int spConform, bool solveRandom,
	      bool do_notreturnparam) {
  // modus:
  //  AS_SAVED : Modell wie gespeichert
  //  DEL_NATSC : Modell unter Annahme PracticalRange>0 (natsc werden geloescht)
  //  SOLVE_NATSC : natscale soweit wie moeglich zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)
  //  DEL_MLE : nur natscale_MLE werden geloescht
  //  SOLVE_MLE : nur natscale_MLE  zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)
  //
  SEXP model, nameMvec;
  int i, nmodelinfo,
    k = 0; 
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  bool plus_mixed_models;


  if ((cov->nr == NATSC_INTERN && modus != GETMODEL_AS_SAVED) ||
      (cov->nr == NATSC_USER && modus == GETMODEL_DEL_NATSC)) { 
    return IGetModel(cov->sub[0], modus, spConform, solveRandom, 
		     do_notreturnparam);
  }
  
  nmodelinfo = C->kappas + 1;
   
   for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL && cov->sub[i]->nr != IDCOORD) // insb. !IDCOORD
      nmodelinfo++;
  }
  for (i=0; i<C->kappas; i++) {
    if ((PisNULL(i) || !STRCMP(C->kappanames[i], INTERNAL_PARAM) ||
	 (do_notreturnparam && SortOf(cov, i, 0, 0) == DONOTRETURNPARAM))
	&& cov->kappasub[i] == NULL)
      nmodelinfo--;
  }

  PROTECT(model = allocVector(VECSXP, nmodelinfo));
  PROTECT(nameMvec = allocVector(STRSXP, nmodelinfo));

  SET_STRING_ELT(nameMvec, k, mkChar("")); // name
  cov_fct *CC = CovList + cov->nr; // nicht gatternr
  while(STRNCMP(CC->name, InternalName, STRLEN(InternalName)) ==0) CC--;
  if ((plus_mixed_models = (cov->nr == PLUS))) {
    plus_mixed_models = cov->calling == NULL; // to do: stimmt nicht mehr
    if (plus_mixed_models) {
      for (i=0; i<MAXSUB; i++) {
	if ((plus_mixed_models = 
	     (cov->sub[i] != NULL && cov->sub[i]->nr == MIXEDEFFECT))) break;
      }
    }
  }

  if (spConform > 1 || (spConform && !plus_mixed_models)) {
    SET_VECTOR_ELT(model, k++, mkString(CC->nick));
  } else {
    SET_VECTOR_ELT(model, k++, mkString(CC->name));
  }

  for(i=0; i<C->kappas; i++) {
    //printf("donotreturnparam = %d\n", do_notreturnparam);

    if ( (!STRCMP(C->kappanames[i], INTERNAL_PARAM) ||
	  PisNULL(i) ||
	  (do_notreturnparam && SortOf(cov, i, 0, 0) == DONOTRETURNPARAM))
	 && cov->kappasub[i] == NULL
      ) continue;
    if (cov->kappasub[i] != NULL && (!solveRandom || PisNULL(i))) {
      SET_STRING_ELT(nameMvec, k, mkChar(C->kappanames[i]));
      SET_VECTOR_ELT(model, k++,
		     IGetModel(cov->kappasub[i], modus, spConform, solveRandom,
			       do_notreturnparam));
      continue;
    }
    SET_STRING_ELT(nameMvec, k, mkChar(C->kappanames[i]));
    SET_VECTOR_ELT(model, k++, 
		   Param(cov, (void*) cov->px[i], cov->nrow[i],
			 cov->ncol[i], C->kappatype[i], true));
  }

  int zaehler = 0;
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL && cov->sub[i]->nr != IDCOORD) {
      SET_STRING_ELT(nameMvec, k, mkChar(C->subnames[i]));
      SET_VECTOR_ELT(model, k++, IGetModel(cov->sub[i], modus, spConform, 
					   solveRandom, do_notreturnparam));
      if (++zaehler >= cov->nsub) break;
    }
  }

  //  printf("%s %d %d %d\n", NAME(cov), k, nmodelinfo, C->kappas);
  //   PMI(cov);
  assert(k == nmodelinfo);

  setAttrib(model, R_NamesSymbol, nameMvec);
  UNPROTECT(2); // model + namemodelvec

  return model;
}

SEXP GetModel(SEXP keynr, SEXP Modus, SEXP SpConform, SEXP whichSub, 
	      SEXP SolveRandom,
	      SEXP Do_notreturnparam) {
  // modus:
  //  AS_SAVED : Modell wie gespeichert
  //  DEL_NATSC : Modell unter Annahme PracticalRange>0 (natsc werden geloescht)
  //  SOLVE_NATSC : natscale soweit wie moeglich zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)
  //  DEL_MLE : nur natscale_MLE werden geloescht
  //  SOLVE_MLE : nur natscale_MLE  zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)

  // whichSub: only "user" and "call+user" are/should be used!

// Nutzer kann 3 Modifikationen des Models in MLE laufen lassen:
//      * keinen Praktikal range oder individuell angeben
//      * extern practical range definieren
//      * intern practical range verwenden lassen (use.naturalscaling)

  int err = NOERROR,
    knr = INTEGER(keynr)[0],  
    spConform = INTEGER(SpConform)[0],
    modus = INTEGER(Modus)[0];
  bool
    solveRandom = LOGICAL(SolveRandom)[0],
    do_notreturnparam = (bool) INTEGER(Do_notreturnparam)[0],
    na_ok = NAOK_RANGE;

   SEXP value = R_NilValue;

   cov_model *cov,
     *dummy = NULL; //ACHTUNG: "=NULL" hinzugefuegt

  if (knr < 0 || knr  > MODEL_MAX || KEY[knr] == NULL) {
    err = ERRORREGISTER;
    goto ErrorHandling;
  }

  //printf("her\n");

  if ((cov = KEY[knr]) == NULL) goto ErrorHandling; // return noerror, nil
  cov = WhichSub(cov, INTEGER(whichSub)[0]);
  if (cov == NULL) BUG;
  //printf("herff\n");
  // printf("cwhich ov = %s %d\n", NAME(cov), INTEGER(whichSub)[0]);

  if (modus == GETMODEL_DEL_NATSC || modus == GETMODEL_DEL_MLE) {
    value = IGetModel(cov, modus, spConform, solveRandom, do_notreturnparam);
  } else {
    if (isInterface(cov)) { 
      if ((err = covCpy(&dummy, true, cov, cov->prevloc, NULL, false, true,
			true)) != NOERROR) goto ErrorHandling;
      dummy->calling = NULL;
    } else {
      if ((err = covCpy(&dummy, cov)) != NOERROR) goto ErrorHandling;
    }
    NAOK_RANGE = true;
    int skipchecks = GLOBAL_UTILS->basic.skipchecks;
    GLOBAL_UTILS->basic.skipchecks = true;
    err = CHECK(dummy, cov->tsdim, cov->xdimprev, cov->typus,
		cov->domprev, cov->isoprev, cov->vdim, cov->role);
    GLOBAL_UTILS->basic.skipchecks = skipchecks;


    //PMI(cov); MERR(err);

    if (err != NOERROR) goto ErrorHandling;    
    iexplDollar(dummy, modus == GETMODEL_SOLVE_MLE);
    if (modus == GETMODEL_SOLVE_NATSC) {
      modus = GETMODEL_DEL_NATSC;
    } else if (modus == GETMODEL_SOLVE_MLE) {
      modus = GETMODEL_DEL_MLE;
    }
    value = IGetModel(dummy, modus, spConform, solveRandom, do_notreturnparam);
   }

 ErrorHandling:
  NAOK_RANGE = na_ok;
  if (dummy != NULL) COV_DELETE_WITHOUT_LOC(&dummy);
  if (err != NOERROR) XERR(err);

  return(value);

}



void ple_intern(cov_fct *C){
  int i;
  PRINTF("pref: ");
  for (i=0; i<Nothing; i++) PRINTF("%d ", C->pref[i]);
  PRINTF("\n");
}

void ple_(cov_model *cov) {
  PRINTF("  %s\n", NAME(cov));
  ple_intern(CovList + cov->nr);
}
 
void ple_(char *name) {
  PRINTF("PLE %s\n", name);
  ple_intern(CovList + getmodelnr(name));
}
 

void PSTOR(cov_model *cov, gen_storage *x) {  
  
  assert(cov != NULL);
  
  int d,
    dim = cov->tsdim;

  if (x==NULL) { PRINTF("no storage information\n"); return; }

 
  for (d=0; d<dim; d++) {
   PRINTF("%d. info:[%3.3f, %3.3f] E=%3.3f cum=%3.3f\n",
	   d, //x->window.min[d], x->window.max[d], x->window.centre[d],
	   RF_NA, RF_NA, // pgs->mppinfo.min[d], x->mppinfo.max[d],
	  x->spec.E[d], x->spec.sub_sd_cum[d]);  
  }

  PRINTF("spec:step=%3.3f phi=%3.3f id=%3.3f grid=%s sig=%3.3f nmetr=%d\n",
	 x->Sspectral.phistep2d, x->Sspectral.phi2d, x->Sspectral.prop_factor,
	 FT[x->Sspectral.grid], x->spec.sigma,x->spec.nmetro);
}


void NoCurrentRegister() { currentRegister=-1; }
void GetCurrentRegister(int *reg) {  *reg = currentRegister; }


bool isSameCoordSystem(isotropy_type iso, coord_sys_enum os) {
  switch(os) {
  case cartesian : case gnomonic: case orthographic : return isCartesian(iso);
  case earth : return isEarth(iso);
  case sphere : return isSpherical(iso);
  case coord_mix : return true;
  default: ERR("unknown coordinate system");
  }
}

coord_sys_enum GetCoordSystem(isotropy_type iso) {
  if (isCartesian(iso)) return cartesian;
  if (isEarth(iso)) return earth;
  if (isSpherical(iso)) return sphere;
  return coord_mix;
}

coord_sys_enum SearchCoordSystem(cov_model *cov, coord_sys_enum os, 
				 coord_sys_enum n_s) {
  if (n_s == coord_keep) {
    if (!isSameCoordSystem(cov->isoown, os)) {
      n_s = GetCoordSystem(cov->isoown);
    }
  } else {
    if (n_s == coord_mix || !isSameCoordSystem(cov->isoown, n_s)) {
      return coord_mix;
    }
  }
  int i;
  coord_sys_enum nn_s;
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      nn_s = SearchCoordSystem(cov->sub[i], os, n_s);
      if (nn_s != n_s) {
	if (n_s == coord_keep) n_s = nn_s;
	else {
	  return coord_mix;
	}
      }
    }
  }
  return n_s;
}


SEXP GetCoordSystem(SEXP keynr, SEXP oldsystem, SEXP newsystem) {
  int knr = INTEGER(keynr)[0];
  cov_model *cov;
  SEXP res;
  char CS[2][30] = {"coordinate system", "new coordinate system"};

  if (knr>=0 && knr <= MODEL_MAX && KEY[knr] != NULL) {
    cov = KEY[knr];
    coord_sys_enum
      os = (coord_sys_enum) GetName(oldsystem, CS[0], 
				    COORD_SYS_NAMES, nr_coord_sys, coord_auto);
    coord_sys_enum
      n_s = (coord_sys_enum) GetName(newsystem, CS[1], 
				    COORD_SYS_NAMES, nr_coord_sys, coord_keep);
   if (os == coord_auto) {
      os = GetCoordSystem(cov->isoprev);
   }
   if (n_s == coord_keep) {
      n_s = SearchCoordSystem(cov, os, n_s);
    }

    if (n_s == coord_mix && GLOBAL.internal.warn_coord_change) {    
      char msg[300];
      SPRINTF(msg, "the covariance model relies on at least two different coordinate systems. Check that this is not due to misspecification of the covariance model. To avoid this warning set 'RFoptions(%s=FALSE)'", 
	      internals[INTERNALS_COORD_CHANGE]);
      warning(msg);
      GLOBAL.internal.warn_coord_change = false;
    }

    bool warn = ((os != coord_auto && os != cartesian) ||
		 (n_s != coord_keep && n_s != os));

    switch (GLOBAL.general.reportcoord) {
    case reportcoord_none :  return R_NilValue;
    case reportcoord_warnings_orally : 
      if (warn) {
	char msg[200];
	SPRINTF(msg, "internal change of coordinate system from '%s' to '%s'. To avoid this message change ",
		COORD_SYS_NAMES[os], COORD_SYS_NAMES[n_s]);
	warning(msg);
      } 
      return R_NilValue;
    case reportcoord_warnings : if (!warn) return R_NilValue; // no break!
    case reportcoord_always :     
      PROTECT(res = allocVector(STRSXP, 2));
      SET_STRING_ELT(res, 0, mkChar(COORD_SYS_NAMES[os]));    
      SET_STRING_ELT(res, 1, mkChar(COORD_SYS_NAMES[n_s]));    
      UNPROTECT(1);
      break;
    default : BUG;
    }
    return res;
  }
  return R_NilValue;
}


void pci(int nr) {
  cov_fct *C = CovList + nr;
  int i;
  PRINTF("%s (%s)\n", C->name, C->nick); 

  PRINTF("  pref:");
  for (i=0; i<=Nothing; i++) PRINTF("%s:%d ", METHODNAMES[i], C->pref[i]);
  PRINTF("\n");    
}

void pci() {
  int nr;
  for (nr = 0; nr<currentNrCov; nr++) PCI(nr);
}

