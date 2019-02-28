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
#include <string.h>
#include "questions.h"
#include <Rdefines.h>
#include <stdint.h>

#include "Coordinate_systems.h"
#include "operator.h"
#include "Processes.h"
#include "startGetNset.h"



// PMI/PMI0
bool PMI_print_dollar = !true,
  PMI_print_mpp = !true,
  PMI_print_pgs = !true, 
  PMI_print_details = !true,
  PMI_print_structure = !true,
  PMI_print_loc = true,
  PMI_print_fields = !true;
int  PMI_print_rect = 1;  // 0, 1, 2
int MAX_PMI = 5;

#define isInternalKappa(i) (!STRCMP(C->kappanames[i], INTERNAL_PARAM))

const char *FTshort[4] = {"-", "-", "F", "T"};
name_type FT = {"false", "true"},
  TriNames = {"mismatch", 
	      "dep. next model",
	      "dep. prev. model",
	      "param. dependent",
	      "false", "true", 
	      "normal mixture",
	      "NaN"};  




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

model *WhichSub(model *cov, int which){    
  if (!equalsnowInterface(cov)) return cov;
  model *ans = cov;

  bool internfst = which == WHICH_INTERNAL || which == WHICH_USER_BUT_ONCE || 
    which == WHICH_USER_BUT_ONCE_JUMP,
    userfst = which == WHICH_USER || which == WHICH_INTERNAL_BUT_ONCE || 
    which==WHICH_INTERNAL_BUT_ONCE_JUMP;
  if (userfst || internfst){
    if (cov->Splus != NULL && cov->Splus->keys_given)
      warning("for '+', it is unclear which path to take");
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


SEXP ScalarExtendedLogical(ext_bool x) {
  if (x == 0 || x == 1) return ScalarLogical(x);
  return ScalarInteger(x);
}
 
#define ExtBoolNames(x)					\
  ((x) == wahr ? "true" :				\
   (x) == falsch ? "false" :				\
   (x) == Paramdep ? "parameter dependent" :		\
   (x) == Submodeldep ? "depdending on submodel(s)" : 	\
   (x) == Huetchenownsize ? "own gridsize" :		\
   "<error: unknown value>")

SEXP RedMat(double* V, int row, int col, bool drop) {
  if (drop) return Num(V, row * col);
  return Mat(V, row, col);
}

SEXP RedMatInt(int* V, int row, int col, bool drop) {
  if (drop) return Int(V, row * col);
  return MatInt(V, row, col);
}

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
     
      SET_VECTOR_ELT(l, k++, RedMat(L->x, L->xdimOZ, laenge, L->xdimOZ == 1));

      SET_STRING_ELT(namevec, k, mkChar("dim"));
      SET_VECTOR_ELT(l, k++, ScalarInteger(L->timespacedim));  
    } else {
      if (L->ly > 0) elmts++;
      elmts += 2;
      PROTECT(namevec = allocVector(STRSXP,  elmts));
      PROTECT(l = allocVector(VECSXP, elmts));
      
      SET_STRING_ELT(namevec, k, mkChar("x"));
      if (L->grid) SET_VECTOR_ELT(l, k++, Mat(L->xgr[0], 3, L->spatialdim));
      else SET_VECTOR_ELT(l, k++, Mat_t(L->x, L->lx, L->xdimOZ));

      if (L->ly > 0) {
	SET_STRING_ELT(namevec, k, mkChar("y"));
	if (L->grid) SET_VECTOR_ELT(l, k++, Mat(L->ygr[0], 3, L->spatialdim));
 	else SET_VECTOR_ELT(l, k++, Mat_t(L->y, L->ly, L->xdimOZ));
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
     "totpts", "distances", "grid", "has.time.comp", // Time geht nicht, da 
      //                                Abkuerzungskollision mit "T"
     "xgr", "x", "T", "ygr", "y"};
  SEXP namevec, l;
  int k,
    nloc = nlocinfo,
    timespacedim = loc->timespacedim,
    spdim = loc->spatialdim;

  if (loc->ly <= 0) nloc -= 2;
  PROTECT(l = allocVector(VECSXP, nloc));
  PROTECT(namevec = allocVector(STRSXP, nloc));
  for (k=0; k<nloc; k++)
    SET_STRING_ELT(namevec, k, mkChar(info[k]));

  k = 0;
  SET_VECTOR_ELT(l, k++, ScalarInteger(timespacedim));
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


SEXP Param(model *cov, void* p, int nrow, int ncol, SEXPTYPE type,
	     bool drop) {
  int i;
  if (p == NULL) {
    return allocVector(REALSXP, 0); 
  } else {
    switch(type) {
    case REALSXP :
      return RedMat((double*) p, nrow, ncol, ncol==1 && drop);
    case INTSXP :
      return RedMatInt((int*) p, nrow, ncol, ncol==1 && drop);
    case STRSXP :
      return MatString((char**) p, nrow, ncol);
    case CLOSXP :
      BUG;
    case VECSXP : 
      if (COVNR == COVARIATE) {      
       return(GetLocationUserInfo(cov->Scovariate->loc));
      } else {
	const char *info[1] = {"R list"};
	return Char(info, 1);
      }
      break;       
    case LANGSXP : 
      {
	return Rf_duplicate(((sexp_type*) p)->sexp);
	//	const char *info[1] = {"R object"}; return Char(info, 1);
      }
      break; 
    case ENVSXP : 
      {
	return Rf_duplicate(((sexp_type*) p)->sexp);
	//	const char *info[1] = {"R environment"};
	//	return Char(info, 1);
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

SEXP Param(model *cov, int i, SEXPTYPE type,  bool drop) {
  return Param(cov, (void*) cov->px[i], cov->nrow[i], cov->ncol[i], type, drop);
}

#define nsimuinfo 3
SEXP GetSimuInfo(simu_storage *simu) {
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



SEXP IGetModelInfo(model *cov, int prlevel, bool both, int spConform, 
		  int whichSub, int Level, sort_origin origin) {
  // level 6 : only isoown

  
  assert(whichSub >=WHICH_USER && whichSub <= WHICH_ALL);
  // whichSub:  0=submodels, 1=keys, 2=both
  
  #include "primitive.h"
 
#define ninfobase 3
#define ninfo0 0
#define ninfo1 6
#define ninfo2 9
#define ninfo3 6
#define ninfo4 4

  /*   !!!!!     ACHTUNG     !!!!!
       Wenn diese Funktion geaendert wird,
       muss auch GetModelInfo geaendert werden
       !!!!!                 !!!!!
  */

  if (cov == NULL) return allocVector(VECSXP, 0);
  SEXP Model, submodels, nameMvec, param, pnames;
  int i, j, nmodelinfo,
    prL = prlevel > 5 ? 0 : prlevel,
    notnull = 0,
    k = 0;
  defn *C = DefList + COVNR; // nicht gatternr
  location_type
    *loc = cov->calling == NULL ? PrevLoc(cov) : Loc(cov);
  bool return_param,
    param_ok[MAXPARAM + 1],
    given_key = cov->Splus != NULL && cov->Splus->keys_given,
    return_key = given_key && whichSub != WHICH_USER,
    return_sub = cov->nsub > 0 && (whichSub != WHICH_INTERNAL || !given_key),
    show_param = prL >= 1;
 
  for (i=0; i<C->kappas; i++) {
    int old_nn = notnull;
    sortsofparam sort = SortOf(cov, i, 1, 1, origin);
    if ((param_ok[i] = sort <= LASTRETURNED)) {    
      if (both) {
	notnull += cov->nrow[i]>0 && cov->ncol[i]>0 && show_param;
	notnull += cov->kappasub[i] != NULL;
      } else {
	notnull += (cov->nrow[i]>0 && cov->ncol[i]>0 && show_param) ||
	  cov->kappasub[i] != NULL;
      }
      param_ok[i] = notnull > old_nn;
    }
  }
  return_param = notnull > 0;

  nmodelinfo = ninfobase + (int) return_param;
  switch(prL) {
  case 5 : nmodelinfo += ninfo4; FALLTHROUGH_OK; 
  case 4 : nmodelinfo += ninfo3; FALLTHROUGH_OK; 
  case 3 : nmodelinfo += ninfo2; FALLTHROUGH_OK; 
  case 2 : nmodelinfo += ninfo1; FALLTHROUGH_OK; 
  case 1 : nmodelinfo += ninfo0; break;
  default: {
    switch(prlevel) {
    case 0 : break;
    case 6 : nmodelinfo += 1; break;
    default : ERR("Possible 'level' values range from 0 to 6.");
    }
  }
  }
  if (!return_sub) nmodelinfo--;
  if (!return_key) nmodelinfo--; 

  PROTECT(Model = allocVector(VECSXP, nmodelinfo));
  PROTECT(nameMvec = allocVector(STRSXP, nmodelinfo));

  SET_STRING_ELT(nameMvec, k, mkChar("name"));
  defn *CC = DefList + COVNR; // nicht gatternr579
  while(STRNCMP(CC->name, InternalName, STRLEN(InternalName)) ==0) CC--;


  if (true || spConform) {
    SET_VECTOR_ELT(Model, k++, mkString(CC->nick));
  } else {
    SET_VECTOR_ELT(Model, k++, mkString(CC->name));
  }

  if (prlevel == 6) {
    SET_STRING_ELT(nameMvec, k, mkChar("coordinates")); 
    SET_VECTOR_ELT(Model, k++,
		   mkString(ISO_NAMES[CoordinateSystemOf(OWNISO(0))]));
  }

 
  if (return_param) {
    PROTECT(param = allocVector(VECSXP, notnull));
    PROTECT(pnames = allocVector(STRSXP, notnull));
    for (j=i=0; i<C->kappas; i++) {
      if (!param_ok[i]) continue;
      if (cov->kappasub[i] != NULL) {	  
	  SET_STRING_ELT(pnames, j, mkChar(OWNKAPPA(C, i)));
	  SET_VECTOR_ELT(param, j,
			 IGetModelInfo(cov->kappasub[i], prlevel, both, 
				      spConform, whichSub, Level + 1, origin)); 
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
		
	SET_STRING_ELT(pnames, j, mkChar(OWNKAPPA(C, i)));
	SET_VECTOR_ELT(param, j, Param(cov, i, C->kappatype[i], true));    
	j++;
      }
    }
    setAttrib(param, R_NamesSymbol, pnames);
    
    SET_STRING_ELT(nameMvec, k, mkChar(prlevel > 5 ? "submodels" : "param"));  
    SET_VECTOR_ELT(Model, k++, param);
    UNPROTECT(2);

  }
  
  if (prL>=2) {      
    SET_STRING_ELT(nameMvec, k, mkChar("covnr"));
    SET_VECTOR_ELT(Model, k++, ScalarInteger(COVNR));

    SET_STRING_ELT(nameMvec, k, mkChar("vdim"));
    SET_VECTOR_ELT(Model, k++,  ScalarInteger(VDIM0)
		   // VDIM0 == VDIM1 
		   // ? ScalarInteger(VDIM0) 
		   // : Int(cov->vdim, 2)
		   );
 
    //  SET_STRING_ELT(nameMvec, k, mkChar("naturalscaling"));  
    //  SET_VECTOR_ELT(Model, k++, ScalarInteger(cov->naturalscaling));

    SET_STRING_ELT(nameMvec, k, mkChar("logicaldim"));  
    SET_VECTOR_ELT(Model, k++, ScalarInteger(OWNLOGDIM(0)));
      
    SET_STRING_ELT(nameMvec, k, mkChar("prev.xdim"));  
    SET_VECTOR_ELT(Model, k++, ScalarInteger(PREVXDIM(0)));

    SET_STRING_ELT(nameMvec, k, mkChar("own.xdim"));  
    SET_VECTOR_ELT(Model, k++, ScalarInteger(OWNXDIM(0)));

    SET_STRING_ELT(nameMvec, k, mkChar("indep.of.x"));  
    SET_VECTOR_ELT(Model, k++, ScalarLogical(cov->matrix_indep_of_x));
  }

  if (prL>=3) {
    SET_STRING_ELT(nameMvec, k, mkChar("type"));  
    SET_VECTOR_ELT(Model, k++, Char(TYPE_NAMES + OWNTYPE(0), 1));

    SET_STRING_ELT(nameMvec, k, mkChar("frame"));  
    SET_VECTOR_ELT(Model, k++, Char(TYPE_NAMES + cov->frame, 1));
 
    SET_STRING_ELT(nameMvec, k, mkChar("domown"));  
    SET_VECTOR_ELT(Model, k++, mkString(DOMAIN_NAMES[OWNDOM(0)]));

    SET_STRING_ELT(nameMvec, k, mkChar("isoown"));
    SET_VECTOR_ELT(Model, k++, mkString(ISO_NAMES[OWNISO(0)]));
 
     if (prL >= 5) {
      SET_STRING_ELT(nameMvec, k, mkChar("isoprev"));
      SET_VECTOR_ELT(Model, k++, mkString(ISO_NAMES[PREVISO(0)]));
    }

    SET_STRING_ELT(nameMvec, k, mkChar("internalq"));  
    SET_VECTOR_ELT(Model, k++, Num(cov->q, cov->qlen));

    SET_STRING_ELT(nameMvec, k, mkChar("storage"));
    if (cov->SlocalCE != NULL) {
      int which = (VDIM0 > 1) * 2 + (CC->check != check_co);
      // 0:co; 1:intr; 2:multvariate co
      SEXP storage, namestorage;
      int L = 0,
	storage_size = 2 + (which <= 1 ? 4 : 8);
      localCE_storage *S0 = cov->SlocalCE;
      localvariab *S = S0->q + 0;
      PROTECT(storage = allocVector(VECSXP, storage_size));
      PROTECT(namestorage = allocVector(STRSXP, storage_size));

      SET_STRING_ELT(namestorage, L, mkChar("R"));  
      SET_VECTOR_ELT(storage, L++, ScalarReal(S->R));

      SET_STRING_ELT(namestorage, L, mkChar("a"));  
      SET_VECTOR_ELT(storage, L++, ScalarInteger(S->a));
      
      switch(which) {
      case 0 :
	SET_STRING_ELT(namestorage, L, mkChar("constant"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cutoff.constant));

	SET_STRING_ELT(namestorage, L, mkChar("b"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cutoff.b));
	
	SET_STRING_ELT(namestorage, L, mkChar("asqrtr"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cutoff.asqrtr));

	SET_STRING_ELT(namestorage, L, mkChar("R_theor"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cutoff.theor));
	break;
      case 1:
	SET_STRING_ELT(namestorage, L, mkChar("A0"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->intrinsic.A0));

	SET_STRING_ELT(namestorage, L, mkChar("A2"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->intrinsic.A2));

	SET_STRING_ELT(namestorage, L, mkChar("B"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->intrinsic.B));

	SET_STRING_ELT(namestorage, L, mkChar("MAX"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->intrinsic.MAX));
	break;
      case 2:
	SET_STRING_ELT(namestorage, L, mkChar("constant"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cube.constant));

	SET_STRING_ELT(namestorage, L, mkChar("R"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cube.R));

	SET_STRING_ELT(namestorage, L, mkChar("A"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cube.A));

	SET_STRING_ELT(namestorage, L, mkChar("B"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cube.B));

	SET_STRING_ELT(namestorage, L, mkChar("C"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cube.C));

	SET_STRING_ELT(namestorage, L, mkChar("N"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cube.N));

	SET_STRING_ELT(namestorage, L, mkChar("M"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cube.M));

	SET_STRING_ELT(namestorage, L, mkChar("L"));  
	SET_VECTOR_ELT(storage, L++, ScalarReal(S->cube.L));
	break;
      default : BUG;
      }
      
      assert(L == storage_size);
      setAttrib(storage, R_NamesSymbol, namestorage);
      SET_VECTOR_ELT(Model, k++, storage);
      UNPROTECT(2);
    } else SET_VECTOR_ELT(Model, k++, allocVector(VECSXP, 0));
   
  
    SET_STRING_ELT(nameMvec, k, mkChar("pref"));  
    SET_VECTOR_ELT(Model, k++, Int(cov->pref, Nothing + 1));
   
    SET_STRING_ELT(nameMvec, k, mkChar("simu"));  
    SET_VECTOR_ELT(Model, k++, GetSimuInfo(&(cov->simu)));
  	
    SET_STRING_ELT(nameMvec, k, mkChar("loc"));
    if (Level == 0 || cov->calling == NULL || Loc(cov->calling) != Loc(cov))
      SET_VECTOR_ELT(Model, k++, GetLocationInfo(loc));
    else
      SET_VECTOR_ELT(Model, k++,
		     mkString(Loc(cov) == NULL ? "no locations given"
			      : "same as calling model"));
  }

  if (prL >= 4) {
   //  printf("k=%d\n", k);
    SET_STRING_ELT(nameMvec, k, mkChar("logspeed"));
    SET_VECTOR_ELT(Model, k++, ScalarReal(cov->logspeed));

    SET_STRING_ELT(nameMvec, k, mkChar("maxdim"));  
    SET_VECTOR_ELT(Model, k++, ScalarInteger(MAXDIM(OWN, 0)));
      
    SET_STRING_ELT(nameMvec, k, mkChar("full_derivs"));  
    SET_VECTOR_ELT(Model, k++, ScalarInteger(cov->full_derivs));
      
    SET_STRING_ELT(nameMvec, k, mkChar("loggiven"));  
    SET_VECTOR_ELT(Model, k++, ScalarLogical(cov->loggiven));

    SET_STRING_ELT(nameMvec, k, mkChar("monotone"));  
    int idx = cov->monotone - (int) MON_UNSET;
    if (idx < 0 || idx > BERNSTEIN - (int) MON_UNSET) {
      PRINTF("monotone %d %d %d\n",  cov->monotone,  MON_UNSET,
	     cov->monotone - (int) MON_UNSET);
      BUG;
    }
    SET_VECTOR_ELT(Model, k++, Char(MONOTONE_NAMES + idx, 1));
             
    SET_STRING_ELT(nameMvec, k, mkChar("finiterange"));  
    SET_VECTOR_ELT(Model, k++, ScalarExtendedLogical(cov->finiterange));
      
    //  printf("ek=%d\n", k);
  }
  
  if (prL>=5) {
    mpp_properties *Mpp = &(cov->mpp);
    int
      vdim = VDIM0,
      nm = cov->mpp.moments,
      nmvdim = (nm + 1) * vdim;
 
    //   SET_STRING_ELT(nameMvec, k, mkChar("mpp.refradius"));  
    //  SET_VECTOR_ELT(Model, k++, ScalarReal(mpp->refradius));
 
    SET_STRING_ELT(nameMvec, k, mkChar("mpp.maxheight"));  
    SET_VECTOR_ELT(Model, k++, Num(Mpp->maxheights, vdim, MAXINT));
  
    SET_STRING_ELT(nameMvec, k, mkChar("mpp.M"));
    SET_STRING_ELT(nameMvec, k+1, mkChar("mpp.Mplus"));
    if (nm == 0) {
      SET_VECTOR_ELT(Model, k++, R_NilValue); 
      SET_VECTOR_ELT(Model, k++, R_NilValue);
    } else {
      SET_VECTOR_ELT(Model, k++, Num(Mpp->mM, nmvdim, MAXINT)); 
      SET_VECTOR_ELT(Model, k++, Num(Mpp->mMplus, nmvdim, MAXINT));
    }
  } 
   		 
  if (return_key) {
    if (cov->key != NULL) {
      SET_STRING_ELT(nameMvec, k, mkChar("internal"));  
      SET_VECTOR_ELT(Model, k++,
		     IGetModelInfo(cov->key, prlevel, both, spConform,
				  whichSub, Level + 1, origin));
    } else {  /// COVNR == PLUS/M && cov->Splus != NULL
      int ii, n,
	subs = C->maxsub;
      SEXP keys;
      SET_STRING_ELT(nameMvec, k, mkChar("internal"));  
      for (ii=n=0; ii<subs; ii++) if (cov->Splus->keys[ii] != NULL) n++;
      PROTECT(keys = allocVector(VECSXP, n));
      for (ii=n=0; ii<subs; ii++) 
	if (cov->Splus->keys[ii] != NULL)
	  SET_VECTOR_ELT(keys, n++, 
			 IGetModelInfo(cov->Splus->keys[ii], prlevel, both,
				      spConform, whichSub, Level + 1, origin));
      SET_VECTOR_ELT(Model, k++, keys);
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
		       IGetModelInfo(cov->sub[i], prlevel, both, spConform, 
				    whichSub, Level + 1, origin));
	if (++zaehler >= cov->nsub) break;
      }
    }
    SET_VECTOR_ELT(Model, k++, submodels);
    UNPROTECT(1);
  }
  
  setAttrib(Model, R_NamesSymbol, nameMvec);

  assert(k==nmodelinfo); 
  UNPROTECT(2); // model + nameModelvec
  return Model;
}

SEXP GetModelInfo(SEXP keynr, SEXP Prlevel, SEXP spConform, SEXP whichSub,
		     SEXP origin) {
  int 
    knr = INTEGER(keynr)[0],
    which = INTEGER(whichSub)[0] == WHICH_ALL_USER 
    ? WHICH_ALL : INTEGER(whichSub)[0] % 2,
    prlevel = std::abs(INTEGER(Prlevel)[0]) % 10;
  bool 
    both = INTEGER(Prlevel)[0] < 0,
    delete_call = std::abs(INTEGER(Prlevel)[0]) < 10;
  model *cov, *orig;
  SEXP res, names;
  model **key = KEY();

  if (knr>=0 && knr <= MODEL_MAX && key[knr] != NULL) {
    orig = key[knr];
    cov = WhichSub(orig, INTEGER(whichSub)[0]);

    res = IGetModelInfo(cov, prlevel, both,
		       (bool) INTEGER(spConform)[0], which, 0,
		       (sort_origin) INTEGER(origin)[0]);

    if (prlevel>=1 && delete_call) {
      PROTECT(names = getAttrib(res, R_NamesSymbol));
      int i, len = length(names);
      for (i=0; i<len; i++) {
	const char *name = CHAR(STRING_ELT(names, i));
	if (STRCMP("prev.xdim", name) == 0) {
	  INTEGER(VECTOR_ELT(res, i))[0] = XDIM(PREVSYSOF(orig), 0);
	  break;
	}
      }
      UNPROTECT(1);
    }
    return res;
  }
  return allocVector(VECSXP, 0);
}


void Path(model *cov, model *sub) {
  defn *C = DefList + COVNR;
  const char *sep="-> ";
  
  //printf("%d %ld %ld\n", cov->calling == NULL, cov, cov->calling);

  if (cov->calling == NULL) {
    PRINTF(" *** "); 
  } else {
    // if (cov->calling == cov) crash();
    assert(cov->calling != cov);
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

  //assert(cov->Splus == NULL || cov->Splus->keys_given);
	 
  
  if (cov->Splus != NULL) { //&& cov->Splus->keys_given) {
    for (i=0; i<C->maxsub; i++) {
      if (cov->Splus->keys[i] == sub) {
	PRINTF("%s.S[%d].%d%s", C->nick, i, cov->zaehler, sep);
	return;
      }
    }
  }

  for (i=0; i<C->kappas; i++) {
    if (cov->kappasub[i] == sub) {
      PRINTF("%s.%s.%d%s", C->nick, OWNKAPPA(C, i), cov->zaehler, sep);
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

void PrintPoints(location_type *loc, char *name, 
		 double *x, coord_type xgr, long lx) {
#ifndef SHOW_ADDRESSES  
  return; // unclear Error below on CRAN
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
    PRINTF("loc:%s [%d] ", name, lx);
    if (lx <= 0) {
      PRINTF("not given! (%d)", addressbits(loc->x));
    } else {
      long total = loc->distances ? lx * (lx-1) / 2 : lx * loc->xdimOZ,
	endfor = total;
       if (endfor > maxpts) endfor = maxpts;
      for (i=0; i<endfor; i++) {
	PRINTF("%4.3f", x[i]);
	if ((i+1) % loc->xdimOZ == 0) { PRINTF(";"); }
	PRINTF(" ");
      }
      if (endfor < total) {
	PRINTF("... [%ld not shown]", total - endfor);
      }
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
  leer(level); PRINTF("%-10s %d\n", "loc:lx", loc->lx);
  leer(level); PRINTF("%-10s %d %d\n","loc:totpts",
		      loc->spatialtotalpoints, loc->totalpoints);
  leer(level); PRINTF("%-10s %d\n","loc:len", loc->len);
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
    leer(level); PRINTF("%-10s (%g %g %g)\n", "loc:T", 
			loc->T[0], loc->T[1], loc->T[2]);
  }

  leer(level); PRINTF("%-10s ","loc:cansio");
  if (loc->caniso == NULL) {
    PRINTF("null\n"); 
  } else {
    int endfor = loc->cani_nrow * loc->cani_ncol;
    PRINTF(" [%d, %d] ",  loc->cani_nrow, loc->cani_ncol);
    if (endfor > MAX_PMI) endfor = MAX_PMI;
    for (i=0; i<endfor; i++) PRINTF(" %g", loc->caniso[i]); 
    PRINTF("\n");
  }
}

void PRINTMAX(double *a, int n, int max) {
  if (n <= max + 2) for (int i=0; i<n; i++) { PRINTF("%g ", a[i]); }
  else {
    for (int i=0; i<max; i++) PRINTF("%g ", a[i]);
    PRINTF("(%d not printed)", max - n);
  }
}
void PRINTMAX(int *a, int n, int max) {
  if (n <= max + 2) for (int i=0; i<n; i++) { PRINTF("%d ", a[i]); }
  else {
    for (int i=0; i<max; i++) PRINTF("%d ", a[i]);
    PRINTF("(%d not printed)", max - n);
  }
}

#define PRINTREAL(X) \
  if (R_FINITE(X)) PRINTF(" %g", X);					\
  else PRINTF(ISNAN(X) ? " NaN" : ISNA(X) ? " NA" : (X) == RF_INF ? " Inf" : \
	      (X) == RF_NEGINF ? " -Inf" : " <unknown>")

#define PRINTINT(X)  \
  if ((X) == NA_INTEGER) { PRINTF(" NA(int)");} else PRINTF(" %d", X)



#define TREE_MAXSTORAGE 20
int getroot(model *root, model *storage[TREE_MAXSTORAGE]) {
  if (root->calling != NULL) {
    int n = getroot(root->calling, storage);
    if (n >= TREE_MAXSTORAGE) BUG;
    storage[n] = root;
    return n+1;
  } 
  storage[0] = root;
  return 1;
}

#define P10(X,Y) PRINTF(FABS(Y) < 1e-6 || FABS(Y) > 1e8 ? "%-10s %10e\n" : "%-10s %g\n", X, Y);

void pmi(model *cov, char all_subs, int level, int maxlevel,
	 int stor_level, model *storage[TREE_MAXSTORAGE]) {    
  int  endfor;
  defn *C = DefList + COVNR; // nicht gatternr
#define MNlength 4
  char MN[Forbidden + 1][MNlength], name[100];

  //if (storage != NULL)
  //printf("storage=%d %s %s \n", stor_level,NAME(storage[stor_level]),NAME(cov));
  if (storage != NULL && storage[stor_level] != cov) return;
 
  //int n = 2;
  if (level > maxlevel) {
    leer(level); PRINTF("'%s' [%d] left out\n", NAME(cov), cov->zaehler);
    // printI(cov);
    return;
  }

  for (int i=0; i<=Forbidden; i++) {
    strcopyN(MN[i], METHOD_NAMES[i], MNlength);
  }

  defn *CC = C;
  while(STRCMP(CC->name, InternalName) ==0) CC--;
  if (level == 0) {
     PRINTF("******   %s (%d; V=%d of %d)   ****** [%d]", CC->nick, COVNR,
	   cov->variant, CC->variants, cov->zaehler);
  } else PRINTF("    **** %s (%d; %d of %d) **** [%d]", CC->nick, COVNR, 
	      cov->variant, CC->variants, cov->zaehler);
  PRINTF("\n");

  //  printI(cov);printD(cov);
   leer(level); PRINTF("%-10s %s\n", "param", C->kappas == 0 ? "none" : ""); 

  for(int i=0; i<C->kappas; i++) {
    //    printf("[[ %s %s ]]  ", C->kappanames[i],
    //	   cov->ownkappanames != NULL ? cov->ownkappanames[i] : "---");
    
    
    STRCPY(name, OWNKAPPA(C, i));
    name[9] = '\0';
    if (PisNULL(i)) {
      if (i>0 && C->kappas > 6 && cov->kappasub[i] == NULL) {
	if (!PisNULL(i-1) || cov->kappasub[i-1] != NULL) leer(level + 1);
	else PRINTF(", ");
	PRINTF("%s", name);
	if (i+1 == C->kappas ||
	    (i+1 < C->kappas && (!PisNULL(i+1) || cov->kappasub[i+1] != NULL))){
	  PRINTF(": NULL");
	} else continue;
      } else {
	leer(level + 1); PRINTF("%-10s", name); 
	PRINTF(" NULL");
      }
    } else { // ! PisNULL
      leer(level + 1); PRINTF("%-10s", name); 
      int size = endfor = NROW(i) * NCOL(i);
      if (C->kappatype[i] < LISTOF) {
	if (endfor > MAX_PMI) endfor = MAX_PMI;
	if (NCOL(i)==1 && NROW(i) > 1) { PRINTF("[%d] ", NROW(i)); }
	else PRINTF("[%d, %d] ", NROW(i), NCOL(i));
      } 
      switch(C->kappatype[i]) {
      case REALSXP : for (int j=0; j<endfor; j++) PRINTREAL(P(i)[j]); break;
      case INTSXP : for (int j=0; j<endfor; j++) PRINTINT(PINT(i)[j]); break;
      case STRSXP : for (int j=0; j<endfor; j++) PRINTF(" %s", PCHAR(i)[j]); break;
      case CLOSXP : PRINTF("< arbitrary function >"); break;
      case VECSXP : PRINTF("< some R list >"); break;
      case LANGSXP: PRINTF("< language expression >"); break;
      case ENVSXP : PRINTF("< R environment >"); break;
      case LISTOF + REALSXP : {
	int ende=0,
	  k_max = 100,
	  k_end = NROW(i);
	listoftype *p= PLIST(i);	
	PRINTF("list [%d]\n", cov->nrow[i]);
	if (k_end > k_max) k_end = k_max;
	for (int k=0; k<k_end; k++) {
	  leer(level + 2); 
   	  if (p->ncol[k]==1) {
            if (p->nrow[k]==1) {
   	      PRINTREAL(p->lpx[k][0]); 
	    } else {
	      PRINTF("[%d] ", p->nrow[k]);
	      ende = endfor = p->nrow[k]; 
	      if (endfor > MAX_PMI) endfor = MAX_PMI;
	      for (int j=0; j<endfor; j++) PRINTREAL(p->lpx[k][j]); 
            }
          } else {
            PRINTF("[%d, %d] ", p->nrow[k],  p->ncol[k]);
            ende = endfor = p->nrow[k] * p->ncol[k]; 
	    if (endfor > MAX_PMI) endfor = MAX_PMI;
	    for (int j=0; j<endfor; j++) PRINTREAL(p->lpx[k][j]); 
          }
	  if (ende > MAX_PMI) { PRINTF("..."); }
	  if (k<cov->nrow[i]-1) { PRINTF("\n"); }
        }
	if (k_end < cov->nrow[i]) { leer(level + 2); PRINTF("...\n"); }
      }
      break;
      default:
	PRINTF("%s (%d)'s kappatype(%d) = %d\n", NAME(cov), COVNR, i,
	       C->kappatype[i]);
	BUG;
      } // switch
      if (C->kappatype[i] < LISTOF && size > MAX_PMI) { PRINTF(" ..."); }
    }
    if (cov->kappasub[i] != NULL) {
      PRINTF("  <=");
      pmi(cov->kappasub[i], all_subs, level + 3, maxlevel, stor_level, storage);
    } 
    PRINTF("\n");
  } // for
  if (cov->Sset != NULL) {
    model *from = cov->Sset->remote;
    leer(level + 1); 
    PRINTF("%-10s%s [s%d]\n", "<remote>", Nick(from), from->zaehler);     
  }

  if (cov->SlocalCE != NULL) {
    int which = (VDIM0 > 1) * 2 + (CC->check != check_co);
    localCE_storage *S0 = cov->SlocalCE;
    localvariab *q = S0->q2;
    leer(level); PRINTF("%-10s\n", "SlocalCE->q2");
    switch(which) {
    case 0 :
      leer(level + 1);PRINTF("%-10s %3.3f %3.3f %3.3f %g\n", "constant",
			     q[0].cutoff.constant, q[1].cutoff.constant,
			     q[2].cutoff.constant, q[3].cutoff.constant);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "b",
			     q[0].cutoff.b, q[1].cutoff.b, q[2].cutoff.b, q[3].cutoff.b);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "asqrtr",
			     q[0].cutoff.asqrtr, q[1].cutoff.asqrtr,
			     q[2].cutoff.asqrtr, q[3].cutoff.asqrtr);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "R_theor",
			     q[0].cutoff.theor, q[1].cutoff.theor,
			     q[2].cutoff.theor, q[3].cutoff.theor);
      break;
    case 1:
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "A0",
			     q[0].intrinsic.A0, q[1].intrinsic.A0,
			     q[2].intrinsic.A0, q[3].intrinsic.A0);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "A2",
			     q[0].intrinsic.A2, q[1].intrinsic.A2,
			     q[2].intrinsic.A2, q[3].intrinsic.A2);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "B",
			     q[0].intrinsic.B, q[1].intrinsic.B,
			     q[2].intrinsic.B, q[3].intrinsic.B);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "MAX",
			     q[0].intrinsic.MAX, q[1].intrinsic.MAX,
			     q[2].intrinsic.MAX, q[3].intrinsic.MAX);
      break;
    case 2:
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "constant",
			     q[0].cube.constant, q[1].cube.constant,
			     q[2].cube.constant, q[3].cube.constant);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "R",
			     q[0].cube.R, q[1].cube.R,
			     q[2].cube.R, q[3].cube.R);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "A",
			     q[0].cube.A, q[1].cube.A, q[2].cube.A,q[3].cube.A);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "B",
			     q[0].cube.B, q[1].cube.B, q[2].cube.B,q[3].cube.B);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "C",
			     q[0].cube.C, q[1].cube.C, q[2].cube.C,q[3].cube.C);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "N",
			     q[0].cube.N, q[1].cube.N, q[2].cube.N,q[3].cube.N);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "M",
			     q[0].cube.M, q[1].cube.M, q[2].cube.M,q[3].cube.M);
      leer(level + 1);PRINTF("%-10s %g %g %g %g\n", "L",
			     q[0].cube.L, q[1].cube.L, q[2].cube.L,q[3].cube.L);
      break;
    default : BUG;
    }
  } else if (cov->qlen > 0) {
    leer(level); PRINTF("%-10s", "internal-q");  
    endfor = cov->qlen; if (endfor > MAX_PMI) endfor = MAX_PMI;
    for (int i=0; i<endfor; i++) PRINTF(" %g", cov->q[i]); 
    PRINTF("\n");
  }

  if (cov->calling == NULL && level != 0) {
    PRINTF("current model is %s\n", NICK(cov)); 
    BUG;
  }
  leer(level); 
  if (cov->calling==NULL) { PRINTF("%-10s %s\n","calling", "NULL");}
  else PRINTF("%-10s %s [%d]\n","calling", 
	      Nick(cov->calling), cov->calling->zaehler);
  //if (cov->root == NULL) crash();
  assert(cov->root != NULL);
  assert(cov->base != NULL);
  if (level == 0 || PMI_print_structure) {
     leer(level);
     PRINTF("%-10s %s [%d]; KEY=%lu (%d)\n", "root",
	    NAME(cov->root), cov->root->zaehler,
	    (uintptr_t) KEYtypeOf(cov), KEYtypeOf(cov)->currentRegister);
  }
   
  int  o,
    g = NRi(GATTER[0]), // OK
    p = NRi(PREV[0]); // OK
  if (PMI_print_structure){
    leer(level); PRINTF("%-10s %s (%d) / %s (%d)\n","prev/#",
		      p==UNSET ? "unset" : p >=0 ? DefList[p].name : "none", p,
		      g==UNSET ? "unset" : g >=0 ? DefList[g].name : "none", g);
  }
  // ACHTUNG: NICHT DIE ABKUERZUNGEN PREVDOM ETC VERWENDEN, DA HIER
  // KEIN PLAUSIBILITAETS C H E C K   GEMACHT WERDEN SOLL
  p = (int) DOMi(PREV[0]); // OK
  g = (int) DOMi(GATTER[0]); // OK
  o = (int) DOMi(OWN[0]); // OK
  leer(level); PRINTF("%-10s %d/%d/%d (%s/%s/%s)\n","dom", p, g, o,
		      DOMAIN_NAMES[(int) p], DOMAIN_NAMES[(int) g], 
		      DOMAIN_NAMES[(int) o]);
  p = (int) ISOi(PREV[0]); // OK
  g = (int) ISOi(GATTER[0]); // OK
  o = (int) ISOi(OWN[0]); // OK
  leer(level); PRINTF("%-10s %d/%d/%d (%s/%s/%s)\n","iso", p, g, o,
		      ISO_NAMES[p], ISO_NAMES[g], ISO_NAMES[o]);
  
  p = (int) XDIMi(PREV[0]); // OK
  g = (int) XDIMi(GATTER[0]); // OK
  o = (int) XDIMi(OWN[0]); // OK
  int lp = (int) LOGDIMi(PREV[0]), // OK
    lg = (int) LOGDIMi(GATTER[0]), // OK
    lo = (int) LOGDIMi(OWN[0]); // OK
  leer(level); PRINTF("%-10s %d:%d:%d, %d:%d:%d, %d,%d\n","l-x-v-dim",
		      lp, lg, lo, p, g, o, VDIM0, VDIM1);  // OK

  if (PMI_print_structure) {
    leer(level); PRINTF("%-10s %s\n","rand.kappa", FT[cov->randomkappa]);
    
    leer(level); PRINTF("%-10s %d\n","variant", (int) cov->variant);
    o = (int) MAXDIMi(OWN[0]); // OK
    leer(level); PRINTF("%-10s %d\n","maxdim", o);
    leer(level); PRINTF("%-10s %d %d %d likeli==%d locCE=%d \n","stor/2/extra",
			addressbits(cov->Sgen),
			addressbits(cov->Sgatter),
			addressbits(cov->Sextra),
			addressbits(cov->Slikelihood),
			addressbits(cov->SlocalCE)
			);
    leer(level); PRINTF("%-10s ","S");
#define S(X) if (cov->S##X == NULL) {} else{ PRINTF("%s, ", #X); }
    S(ce); S(localCE); S(approxCE); S(direct); S(hyper); S(nugget);
    S(plus); S(sequ); S(tbm); S(trend); S(br); S(get); S(pgs); S(set);
    S(polygon); S(rect); S(dollar); S(earth); S(extra); S(solve);
    S(biwm); S(bistable); S(scatter); S(mcmc); S(gen); S(likelihood);
    S(covariate); S(mle);     
    if (cov->Sgatter == NULL) { PRINTF("gatter=NULL!, "); }
    PRINTF("\n");
    if (cov->Splus != NULL) {
      leer(level); PRINTF("%-10s ","conform");
      bool *conform = cov->Splus->conform;
      for (int i=0; i<cov->nsub; i++) PRINTF("%s", conform[i] ? "T" : "F");
      PRINTF("\n");
    }
    leer(level); PRINTF("%-10s %s\n","simu:activ", FT[cov->simu.active]);
    leer(level); PRINTF("%-10s %s\n","     pair", FT[cov->simu.pair]);
    leer(level); PRINTF("%-10s %d\n","     expect", cov->simu.expected_number_simu);    
  }
  p = (int) TYPEi(PREV[0]); // OK
  g = (int) TYPEi(GATTER[0]); // OK
  o = (int) TYPEi(OWN[0]); // OK
  leer(level); PRINTF("%-10s %d/%d/%d %s/%s/%s\n", "type", p, g, o,
		      TYPE_NAMES[p], TYPE_NAMES[g], TYPE_NAMES[o]);
  leer(level); PRINTF("%-10s %s\n","frame", TYPE_NAMES[cov->frame]);
  if (cov->err != NOERROR) {
    leer(level); PRINTF(" ! %-s ! level=%d err=%d (%s)\n","error", cov->err_level, cov->err, cov->err_msg);
  }


  if (PMI_print_details) {
    leer(level); PRINTF("%-10s %s\n","user_given", 
		      cov->user_given == ug_internal ? "internal" :
		      cov->user_given == ug_explicit ? "explicit" : "implicit");
    leer(level); P10("logspeed", cov->logspeed);  
    leer(level); PRINTF("%-10s %d, %d\n","full/rese deriv's", 
			cov->full_derivs, cov->rese_derivs);  
    leer(level); PRINTF("%-10s %s\n","loggiven", FT[cov->loggiven]);  
    int idx = cov->monotone - (int) MON_UNSET;
    if (idx < 0 || idx > BERNSTEIN - (int) MON_UNSET) {
      PRINTF("monotone %d %d %d\n",  cov->monotone,  MON_UNSET,
	     cov->monotone - (int) MON_UNSET);
      BUG;
    }
    leer(level); PRINTF("%-10s %s (%d)\n","monotone", 
			MONOTONE_NAMES[idx], cov->monotone);  
    leer(level); PRINTF("%-10s %s (%d)\n","pointwise", 
			POSITIVITY_NAMES[cov->ptwise_definite],
			cov->ptwise_definite);  
    leer(level); PRINTF("%-10s %s\n","finiterng", 
			TriNames[cov->finiterange - MISMATCH]);
   

    if (cov->taylorN > 0) {
      for (int d=0; d<cov->taylorN; d++) {
	leer(level); 
	PRINTF("%-10s c%d=%g p%d=%g\n", d==0 ? "taylor" : "", 
	       d, cov->taylor[d][TaylorConst], d, cov->taylor[d][TaylorPow]);
      }
    } else {leer(level); PRINTF("%-10s %s\n", "taylor", "undetermined");}
    if (cov->tailN > 0) {
	for (int d=0; d<cov->tailN; d++) {
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
  if (PMI_print_fields) {
    leer(level); PRINTF("%-10s %s\n","method", METHOD_NAMES[cov->method]);
    leer(level); PRINTF("%-10s %s\n","initialised", FT[cov->initialised]);
    leer(level); PRINTF("%-10s %s\n","fieldret",ExtBoolNames(cov->fieldreturn));
    leer(level); PRINTF("%-10s %s\n","origrf", FT[cov->origrf]);
    leer(level); PRINTF("%-10s %d\n","rf", addressbits(cov->rf));
    leer(level); PRINTF("%-10s %s\n","checked", FT[cov->checked]);
    leer(level); PRINTF("%-10s ", "pref");
    int i;
    for (i=0; i<=Sequential; i++) PRINTF("%s:%d ", MN[i], (int) cov->pref[i]);
    PRINTF("\n"); leer(level); PRINTF("%-10s ", "");  
    for (; i<=Nothing; i++) PRINTF("%s:%d ", MN[i], (int) cov->pref[i]);
    PRINTF("\n");
  }
  
  //leer(level); PRINTF("%-10s %d\n","naturalscaling", cov->naturalscaling);  
  //leer(level); PRINTF("%-10s %d\n","init",DefList[COVNR].init!=init_failed);
  //leer(level); PRINTF("%-10s %d\n","spec:nmetro", cov->spec.nmetro);
  //leer(level); P10("spec:sigma", cov->spec.sigma);
  if (PMI_print_dollar && cov->Sdollar != NULL) {
    leer(level); PRINTF("%-10s ", "separabel", cov->Sdollar->separable);
    // currently none
  } 
  if (PMI_print_details) {
    if (cov->Sextra != NULL) {
      extra_storage *s = cov->Sextra;
#define PE(a1)								\
      if (s->a1 == NULL) {\
	if (s->n_##a1!=0) { ERR1("%.50s=NULL; n<>0 !! ERROR !!\n", #a1); } \
      } else if (s->n_##a1 == 0) { ERR1("%.50s<>NULL; n=0 !! ERROR !!\n", #a1); } \
      else PRINTMAX(s->a1, s->n_##a1, 10)
      PE(a1);PE(a2);PE(a3);PE(b1);PE(b2);PE(c1);PE(c2);
      PE(i1);//PE(i2);
      PE(j1);//PE(j2);
      PE(k1);PE(k2);
    }
  }


  if (PMI_print_mpp) {
    int 
      nm = cov->mpp.moments,
      vdim = VDIM0,
      nmvdim = (nm + 1) * vdim;
    if (R_FINITE(cov->mpp.maxheights[0])) {
      leer(level); PRINTF("%-10s ","mpp:maxhgt"); 
      for (int i=0; i<=vdim; i++) PRINTF("%g, ",cov->mpp.maxheights[i]);
      PRINTF("\n");    
    }
    if (R_FINITE(cov->mpp.unnormedmass)) {
      leer(level); P10("mpp:u-mass", cov->mpp.unnormedmass);    
    }
    leer(level); PRINTF("%-10s ","mpp:M+");
    if (cov->mpp.mMplus == NULL) 
      PRINTF("not initialized yet (size=%d)\n", nm);
    else {
      for (int i=0; i<nmvdim; i++) PRINTF("%g, ", cov->mpp.mMplus[i]);
      PRINTF("\n");    
    }
    leer(level); PRINTF("%-10s ","mpp:M");
    if (cov->mpp.mM == NULL) {
      PRINTF("not initialized yet (size=%d)\n", nm);
    } else {
      for (int i=0; i<nmvdim; i++)  PRINTF("%g, ", cov->mpp.mM[i]);
      PRINTF("\n");    
    }
    leer(level); PRINTF("%-10s %s\n","mpp:log",
			FT[DefList[COVNR].log != ErrLogCov]);
    leer(level); PRINTF("%-10s %s\n","mpplgnonst", 
			FT[DefList[COVNR].nonstatlog != ErrLogCovNonstat]);
    leer(level); PRINTF("%-10s %s\n","mpp:do",
			FT[DefList[COVNR].Do != do_failed]);
  }

  
  if (PMI_print_pgs) {
    if (cov->Spgs == NULL) {
    } else {
      pgs_storage *pgs = cov->Spgs;
      int 
	size = pgs->size,
	dim = OWNTOTALXDIM;
      
      leer(level); P10("pgs:mass", pgs->totalmass);
      leer(level); P10("pgs:logdens", pgs->log_density);
      
      leer(level); PRINTF("%-10s %d\n","pgs:size", size);
      leer(level); PRINTF("%-10s %d\n","pgs:dim", dim);
      leer(level); P10("pgs:threshold", pgs->currentthreshold);
      
#define SHOWDEFAULT(Z, X, Y) if (pgs->X != NULL) {			\
	leer(level);  PRINTF("%-10s ",Z);				\
	for (int d=0; d<dim; d++) PRINTF(Y, pgs->X[d]); PRINTF("\n");}
#define SHOW(Z, X) SHOWDEFAULT(Z, X, "%g ")
#define SHOWINT(Z, X) SHOWDEFAULT(Z, X, "%d ")
      
      if (R_FINITE(pgs->zhou_c)) {
	leer(level); P10("mpp:zhou_c", pgs->zhou_c);
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
	leer(level); PRINTF("%-10s %s\n","pgs:flathull", FT[pgs->flathull]);
	leer(level); P10("pgs:globmin", pgs->globalmin);
	leer(level); P10("pgs:cur.thres",pgs->currentthreshold);
	SHOW("pgs:half", halfstepvector);
	if (pgs->single != NULL) {
	  leer(level); { PRINTF("%-10s ","pgs:single"); 
	    for (int d=0; d<size; d++) PRINTF("%g ", pgs->single[d]);PRINTF("\n"); }
	}
	if (pgs->total != NULL) {
	  leer(level); { PRINTF("%-10s ","pgs:total"); 
	    for (int d=0; d<size; d++) PRINTF("%g ", pgs->total[d]); PRINTF("\n"); }
	}
      } else { // gauss oder poisson
	leer(level); P10("pgs:intens", pgs->intensity);
      }
      if (pgs->cov != NULL) {
	leer(level); 
	PRINTF("%-10s %s [%d]\n","pgs:cov", Nick(pgs->cov), pgs->cov->zaehler);
      }
    }
  }
  
  if (PMI_print_rect) {
    if (cov->Srect != NULL) {
      rect_storage *r = cov->Srect;
      int
	nstepP2 = r->nstep + 2;
      
      leer(level); PRINTF("%-10s %g c=%g p=%g\n","rct:inner",
			  r->inner, r->inner_const, r->inner_pow);
      leer(level); PRINTF("%-10s %g c=%g p=%g pc=%g\n","rct:outer",
			  r->outer, r->outer_const, r->outer_pow,
			  r->outer_pow_const);
      leer(level); P10("rct:step", r->step);
      leer(level); PRINTF("%-10s %d\n","rct:nstep", r->nstep);
      leer(level); PRINTF("%-10s %d\n","rct:ntmp", r->tmp_n);
      // leer(level); P10("rct:total", r->weight[1 + r->nstep]);
      
      if (PMI_print_rect > 1 && r->value != NULL) {
	leer(level); { PRINTF("%-10s ","rct:val"); 
	  for (int d=0; d<nstepP2; d++) PRINTF("%4.3f ", r->value[d]);
	  PRINTF("\n"); }
	leer(level); { PRINTF("%-10s ","rct:wght"); 
	  for (int d=0; d<nstepP2; d++) PRINTF("%4.4f ", r->weight[d]);
	  PRINTF("\n"); }
      }
    }
  }
    
  if (PMI_print_loc) {
    leer(level);  
    //    printf(" **** %ld %ld %ld\n", PLoc(cov), cov->ownloc, cov->prevloc);
   if (PLoc(cov) == NULL) {
      PRINTF("%-10s %s\n", "loc:", "<null>");
    } else {
     //     printf(" **** %ld %ld %ld %d \n", PLoc(cov), cov->ownloc, cov->prevloc,  cov->prevloc[0]->len);
     assert(PLoc(cov)[0] != NULL);
     PRINTF("%-10s %d\n", "loc:sets", (PLoc(cov)[0])->len);
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
	if (level == 0) PrintLoc(4, prevloc, false);
      }
    }
  }

  bool givenkey;
  if ((givenkey = cov->key != NULL)) {
    leer(level);
    PRINTF("%-10s :", "key");  
    if (all_subs >= 0) pmi(cov->key, all_subs, level + 1, maxlevel,
			   stor_level, storage);
  }

  bool Splus_given = cov->Splus != NULL && cov->Splus->keys_given;
  if (Splus_given) {
    givenkey = true;
    if (all_subs >= 0) {
      for (int i=0; i < cov->nsub; i++) {
	model *key = cov->Splus->keys[i];
	leer(level);      
	if (key != NULL) {
	    PRINTF("%-10s ++ %d ++:", "plus.keys", i); 
	    pmi(key, all_subs, level + 1, maxlevel, stor_level, storage);
	}  else PRINTF("%-10s ++ %d ++: %s\n","plus.key", i, "empty");	
      }
    }
  }

  if ((!givenkey && all_subs==0) || all_subs>0) {
    for (int i=0; i<C->maxsub; i++) {
      if (cov->sub[i] == NULL) {
	continue;
      }
      leer(level); 
      PRINTF("%s %d (%s) of '%s':", "submodel", i, C->subnames[i], C->nick);  
      pmi(cov->sub[i], all_subs, level + 1, maxlevel, stor_level, storage);
    }
  }
}


void pmi(model *cov, int maxlevel) { // OK
  PRINTF("\n");

  if (cov == NULL) {
    PRINTF("\nCovariance model is empty.\n\n");
  } else {
    Path(cov, NULL);
    pmi(cov, !false, 0, maxlevel, 0, NULL);
  }
}



void pmiroot(model *cov, int maxlevel) { // OK
  model *storage[TREE_MAXSTORAGE];
  getroot(cov, storage);
  pmi(storage[0], false, 0, maxlevel, 0, NULL); //storage);
}

void iexplDollar(model *cov, bool MLEnatsc_only) {
  /*    
	get the naturalscaling values and devide the preceeding scale model     
	by this value
  */
  double *p, invscale;
  model *dollar = cov->calling;
  bool solving = dollar != NULL && isDollar(dollar) &&
    (COVNR == NATSC_INTERN || (COVNR == NATSC_USER && !MLEnatsc_only));
  
  if (solving) {
    model 
      *next = cov->sub[0];
    assert(dollar!=NULL && isDollar(dollar));

    INVERSE(&GLOBAL.gauss.approx_zero, next, &invscale);
    //invscale = 1.0;
    if (ISNAN(invscale))
      ERR("inverse function unknown when calculating explicite scale");
    
    p = PARAM(dollar, DSCALE);
    if (p != NULL) {
       p[0] /= invscale;
    } else {
       p = PARAM(dollar, DANISO);      
      if (p != NULL) { 	
	int n = dollar->nrow[DANISO] * dollar->ncol[DANISO];
	for (int i=0; i<n; i++) p[i] *= invscale;
      } else {
	assert(COVNR == NATSC_USER);
      }
    }
  } else {
    for (int i=0; i<MAXSUB; i++) { // cov->sub[i]: luecken erlaubt bei PSgen !
      if (cov->sub[i] != NULL) iexplDollar(cov->sub[i], MLEnatsc_only);
    }
  }

}


  
SEXP IGetModel(model *cov, int modus, int spConform, bool solveRandom,
	       sortsofparam whichparam, sort_origin origin) {
  // modus:
  //  AS_SAVED : Modell wie gespeichert
  //  DEL_NATSC : Modell unter Annahme PracticalRange>0 (natsc werden geloescht)
  //  SOLVE_NATSC : natscale soweit wie moeglich zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)
  //  DEL_MLE : nur natscale_MLE werden geloescht
  //  SOLVE_MLE : nur natscale_MLE  zusammengezogen (natsc werden
  //               drauf multipliziert; Rest wie gespeichert)
  //
  SEXP Model, nameMvec;
  int i, nmodelinfo,
    k = 0; 
  defn *C = DefList + COVNR; // nicht gatternr
  bool ok[MAXPARAM];

  if ((COVNR == NATSC_INTERN && modus != GETMODEL_AS_SAVED) ||
      (COVNR == NATSC_USER && modus == GETMODEL_DEL_NATSC)) { 
    return IGetModel(cov->sub[0], modus, spConform, solveRandom, whichparam,
		     origin);
  }
  
  nmodelinfo = C->kappas + 1;
   
   for (i=0; i<MAXSUB; i++) {
     if (cov->sub[i] != NULL //&& MODELNR(cov->sub[i]) != IDCOORD
	) // insb. !IDCOORD
      nmodelinfo++;
  }
  for (i=0; i<C->kappas; i++) {
    if (cov->kappasub[i] != NULL && (!solveRandom || PisNULL(i))) continue;
    ok[i] = !PisNULL(i) && whichparam != NOPARAMETERS;
    if (!ok[i]) { nmodelinfo--; continue; }
    sortsofparam sort = SortOf(cov, i, 0, 0, origin);
    ok[i] = whichparam == ALLPARAMETERS ||
      (whichparam == STANDARD && sort <= LASTUSERSORTOF) ||
      (whichparam == INCLUDENOTRETURN &&
       (sort <= LASTUSERSORTOF || sort == IGNOREPARAM)) ||
      (whichparam == INTERNALPARAMETERS && isInternalKappa(i));
    if (ok[i]) continue;
    ok[i] = SortOf(cov, i, 0, 0, origin) == whichparam;
    if (!ok[i]) nmodelinfo--;
  }

  PROTECT(Model = allocVector(VECSXP, nmodelinfo));
  PROTECT(nameMvec = allocVector(STRSXP, nmodelinfo));

  SET_STRING_ELT(nameMvec, k, mkChar("")); // name
  defn *CC = DefList + COVNR; // nicht gatternr
   while(STRNCMP(CC->name, InternalName, STRLEN(InternalName)) ==0) CC--;

   SET_VECTOR_ELT(Model, k++, mkString(CC->nick));

   //   if (spConform >= 1) { SET_VECTOR_ELT(Model, k++, mkString(CC->nick));
   // } else {  SET_VECTOR_ELT(Model, k++, mkString(CC->name)); }

  for(i=0; i<C->kappas; i++) {
    if (cov->kappasub[i] != NULL && (!solveRandom || PisNULL(i))) {
      SET_STRING_ELT(nameMvec, k, mkChar(OWNKAPPA(C, i)));
      SET_VECTOR_ELT(Model, k++,
		     IGetModel(cov->kappasub[i], modus, spConform, solveRandom,
			       whichparam, origin));
      continue;
    }
    if (!ok[i]) continue;
    
    SET_STRING_ELT(nameMvec, k, mkChar(OWNKAPPA(C, i)));    
    if (C->kappaParamType[i] == CharInputType ) {
      int value = MISMATCH;
      if (C->kappatype[i] == INTSXP) value = P0INT(i); else
	if (C->kappatype[i] == REALSXP) value = P0(i); else BUG;
      assert(NROW(i) == 1 && NCOL(i) == 1);
      SET_VECTOR_ELT(Model, k++, 
		     Param(cov, C->kappaParamTypeNames[i] + value,
			   NROW(i), NCOL(i), STRSXP, true));
    } else if (C->kappaParamType[i] == MixedInputType &&
	       ( (C->kappatype[i] == INTSXP && P0INT(i) < 0) ||
		 (C->kappatype[i] == REALSXP && P0(i) < 0.0) )) {// pos numbers and string
      int value = MISMATCH;
      if (C->kappatype[i] == INTSXP) value -= P0INT(i); else
	if (C->kappatype[i] == REALSXP) value -= P0(i); else BUG;
      assert(NROW(i) == 1 && NCOL(i) == 1);
      SET_VECTOR_ELT(Model, k++, 
		     Param(cov, C->kappaParamTypeNames[i] + value,
			   NROW(i), NCOL(i), STRSXP, true));
    } else if (C->kappaParamType[i] >= NN1) {
      assert(NROW(i) == 1);
      SET_VECTOR_ELT(Model, k++, 
		     Param(cov,
			   LIST_OF_NAMES[C->kappaParamType[i] - NN1] + P0INT(i),
			   NROW(i), NCOL(i), STRSXP, true));
    } else {
      SET_VECTOR_ELT(Model, k++, Param(cov, i, C->kappatype[i], true));
    }
  }

  int zaehler = 0;
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL // && MODELNR(cov->sub[i]) != IDCOORD
	) {
      SET_STRING_ELT(nameMvec, k, mkChar(C->subnames[i]));
      SET_VECTOR_ELT(Model, k++, IGetModel(cov->sub[i], modus, spConform, 
					   solveRandom, whichparam, original));
      if (++zaehler >= cov->nsub) break;
    }
  }

  assert(k == nmodelinfo);

  setAttrib(Model, R_NamesSymbol, nameMvec);
  UNPROTECT(2); // model + namemodelvec

  return Model;
}


SEXP GetModel(SEXP keynr, SEXP Modus, SEXP SpConform, SEXP whichSub, 
	      SEXP SolveRandom, SEXP Whichparam, SEXP origin) {
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

  int
    knr = INTEGER(keynr)[0],  
    spConform = INTEGER(SpConform)[0],
    modus = INTEGER(Modus)[0],
    whichparam = INTEGER(Whichparam)[0];
  bool
    solveRandom = LOGICAL(SolveRandom)[0];
 
  model *cov = NULL,
    **key = KEY();
 
  if (whichparam < 0 || whichparam > LASTSORTOF) XERR(ERRORSORTOF);
  if (knr < 0 || knr  > MODEL_MAX || key[knr] == NULL) XERR(ERRORREGISTER);
  if ((cov = key[knr]) == NULL) return R_NilValue; // return noError, nil
  
  cov = WhichSub(cov, INTEGER(whichSub)[0]);
  if (cov == NULL) BUG;
  bool na_ok =NAOK_RANGE;//erst nach cov definierbar; KEYtypeOf(cov)->naok_range

  // dummy erst abjetzt 
  if (modus == GETMODEL_DEL_NATSC || modus == GETMODEL_DEL_MLE)
    return IGetModel(cov, modus, spConform, solveRandom,
		     (sortsofparam) whichparam,
		     (sort_origin) INTEGER(origin)[0]);

  //    TREE(cov);
  SEXP value = R_NilValue;
  int  err = NOERROR,
    skipchecks = GLOBAL_UTILS->basic.skipchecks;
  model *dummy = NULL;
  if (equalsnowInterface(cov))
    err = covcpy(&dummy, true, cov, cov->prevloc, NULL, false, true, true);
  else err = covcpy(&dummy, cov);
  if (err != NOERROR) goto ErrorHandling;
  set_NAOK_RANGE(true);
  GLOBAL_UTILS->basic.skipchecks = true;
  SET_CALLING_NULL(dummy, cov);
  err = CHECK_ONLY(dummy);
  GLOBAL_UTILS->basic.skipchecks = skipchecks;
  // TREE(dummy);
  if (err != NOERROR) goto ErrorHandling;
  iexplDollar(dummy, modus == GETMODEL_SOLVE_MLE);
  if (modus == GETMODEL_SOLVE_NATSC) {
    modus = GETMODEL_DEL_NATSC;
  } else if (modus == GETMODEL_SOLVE_MLE) {
    modus = GETMODEL_DEL_MLE;
  }
  PROTECT(value = IGetModel(dummy, modus, spConform, solveRandom,
			    (sortsofparam) whichparam,
			    (sort_origin) INTEGER(origin)[0]));
  if (dummy != NULL) COV_DELETE_WITHOUT_LOC(&dummy, cov);
  UNPROTECT(1);
  set_NAOK_RANGE(na_ok);
  return(value);

 ErrorHandling:
  set_NAOK_RANGE(na_ok);
  if (dummy != NULL) COV_DELETE_WITHOUT_LOC(&dummy, cov);
  XERR(err);
}


void ple_intern(defn *C){
  int i;
  PRINTF("pref: ");
  for (i=0; i<Nothing; i++) PRINTF("%d ", C->pref[i]);
  PRINTF("\n");
}


void ple_(model *cov) {
  PRINTF("  %s\n", NAME(cov));
  ple_intern(DefList + COVNR);
}

 
void ple_(char *name) {
  PRINTF("PLE %s\n", name);
  ple_intern(DefList + getmodelnr(name));
}
 

void PSTOR(model *cov, gen_storage *x) {   
  if (x==NULL) { PRINTF("no storage information\n"); return; }
  assert(cov != NULL);
  
  int dim = OWNLOGDIM(0);

 
  for (int d=0; d<dim; d++) {
   PRINTF("%d. info:[%3.3f, %3.3f] E=%3.3f cum=%3.3f\n",
	   d, //x->window.min[d], x->window.max[d], x->window.centre[d],
	   RF_NA, RF_NA, // pgs->mppinfo.min[d], x->mppinfo.max[d],
	  x->spec.E[d], x->spec.sub_sd_cum[d]);  
  }

  PRINTF("spec:step=%3.3f phi=%3.3f id=%3.3f grid=%s sig=%3.3f nmetr=%d\n",
	 x->Sspectral.phistep2d, x->Sspectral.phi2d, x->Sspectral.prop_factor,
	 FT[x->Sspectral.grid], x->spec.sigma,x->spec.nmetro);
}


void NoCurrentRegister() { set_currentRegister(UNSET); }
void GetCurrentRegister(int *reg) {  *reg = currentRegister(); }


void psys(system_type *sys, bool earlybreak) {
  int last = LASTi(sys[0]);
  const char *nl1[2] = {" ", "\n"};
  const char *nl2[2] = {"\n       ", " "};
  if (last < 0) last = 0;
  for (int s=0; s<=last; s++) {
    PRINTF(" s=%d(%d): nr=%d log=%d x=%d%smax=%d cum=%d%s %d:'%s' %d:'%s' %d:'%s'\n",
	   s, LASTi(sys[s]), NRi(sys[s]), LOGDIMi(sys[s]), XDIMi(sys[s]), // OK
	   nl1[earlybreak], MAXDIMi(sys[s]),  // OK
	   CUMXMITi(sys[s]), nl2[earlybreak], // OK
	   TYPEi(sys[s]), TYPE_NAMES[TYPEi(sys[s])],  // OK
	   DOMi(sys[s]), DOMAIN_NAMES[DOMi(sys[s])],  // OK
	   ISOi(sys[s]), ISO_NAMES[ISOi(sys[s])]); // OK
  }
}

void psys(system_type *sys) { psys(sys, true); }

void psys(model *cov) {
  PRINTF("\nprev:  "); psys(PREV, false);
  PRINTF("gatter:"); psys(GATTER, false); 
  PRINTF("own:   "); psys(OWN, false);
}





void pcl(int nr) {
  defn *C = DefList + nr;
  int i;
  PRINTF("%s (%s)\n", C->name, C->nick); 

  PRINTF("  pref:");
  for (i=0; i<=Nothing; i++) {
    PRINTF("%s:%d ", METHOD_NAMES[i], C->pref[i]);
    if (i == Sequential) { PRINTF("\n\t"); }
  }
  PRINTF("\n");
  int v = C->variants;
  PRINTF("derivs full=%d rese=%d var=%d\n", C->F_derivs, C->RS_derivs, v);
  for (i=0; i<v; i++) {
    system_type *s = C->systems[v];
    PRINTF("type=%d (%s) ", s->type, TYPE_NAMES[s->type]);
    PRINTF("dom=%d (%s) ", s->dom, DOMAIN_NAMES[s->dom]);
    PRINTF("\n");
  }
}

void pcl(model *cov) { pcl(COVNR); }

void pcl() {
  int nr;
  for (nr = 0; nr<currentNrCov; nr++) pcl(nr);
}



bool tree(model *cov, int current, char all_subs, int level,    // TREE
	  model *storage[TREE_MAXSTORAGE], int n, bool alle) {
  defn *C = DefList + COVNR, // nicht gatternr
    *CC = C;
  bool found,
    Splus_given = cov->Splus != NULL; // && cov->Splus->keys_given;
  while(STRCMP(CC->name, InternalName) ==0) CC--;
  PRINTF("%s (V%d of %d) [%d", CC->name, //COVNR, 
	 cov->variant, CC->variants, cov->zaehler);
  if (true) {
    PRINTF(";%s%s", FTshort[cov->DallowedDone + 2 * (C->Dallowed != NULL)],
	   FTshort[cov->IallowedDone + 2 * (C->Iallowed != NULL)]);
  }
  if (true) {
    PRINTF(";%d%d,%d%d", CONDPREVDOM(0), CONDPREVISO(0),
	   CONDOWNDOM(0), CONDOWNISO(0));
  }
  PRINTF("] %d %d", cov->zaehler, current);
  if ((found = cov->zaehler == current)) { PRINTF("\t***"); }
  else if (level < n && storage[level] == cov) {
    if (!false) { 
      for (int j=level; j<=n; j++) PRINTF("   ");
      PRINTF("\t::");
    }
  } else if (!alle) {
    bool further = cov->nsub > 0;
    if (!further)
      for (int i=0; i<C->kappas; i++)
	if ((further = cov->kappasub[i] != NULL)) break;
    if (!further && Splus_given)
      for (int i=0; i < cov->nsub; i++)
  	if ((further = cov->Splus->keys[i] != NULL)) break;
    if (further) { PRINTF(" ..."); }
    PRINTF("\n");
    return false;
  }
  PRINTF("\n");

  for (int i=0; i<C->kappas; i++) {
    if (cov->kappasub[i] != NULL) {
      leer(level); PRINTF("%s (%d): ", C->kappanames[i], i);
      assert(cov->kappasub[i]->calling == cov);
      found |= tree(cov->kappasub[i], current, all_subs, level + 1,// TREE
		    storage, n, alle);
    } 
  }
 
  bool givenkey;
  if ((givenkey = cov->key != NULL) && all_subs >= 0) {
    leer(level); PRINTF("key: ");
    assert(cov->key->calling == cov);
    found |= tree(cov->key, current, all_subs, level+1, storage, n, alle);//TREE
  }

  //printf("%s %d %d\n", NAME(cov), Splus_given, all_subs);
  if (Splus_given && all_subs >= 0) 
    for (int i=0; i < cov->nsub; i++) {
      model *key = cov->Splus->keys[i];
      if (key != NULL) {
	givenkey = true;
	leer(level); PRINTF("array (%d): ", i); 
	assert(key->calling == cov);
	found |= tree(key, current, all_subs, level+1, storage, n, alle);//TREE
      }
    }

  if ((!givenkey && all_subs==0) || all_subs>0 || !found) 
    for (int i=0; i<C->maxsub; i++) {
      if (cov->sub[i] != NULL) {
	leer(level);  PRINTF("%s (%d): ", C->subnames[i], i);
	if (cov->sub[i]->calling != cov) {
	  PRINTF("%s(%d) -> %s(%d)\n", NAME(cov), cov->zaehler,
		 NAME(cov->sub[i]), (cov->sub[i]->zaehler)); // crash();
	}
	assert(cov->sub[i]->calling == cov);
	found |= tree(cov->sub[i], current, all_subs, level + 1,// TREE
		      storage, n, alle);
      }
    }
  return found;
}

void tree(model *cov, bool alle) {// TREE
  int current = cov->zaehler,
    all_subs = 0;
  // all_subs = 1; printf("all_subs=%d\n", all_subs);
  model *storage[TREE_MAXSTORAGE];
  int n = getroot(cov, storage);
  if (!tree(storage[0], current, all_subs, 0, storage, n, alle) || // TREE
      all_subs < 0) {
    //  printf("found = %d all_subs=%d\n",
    // 	   tree(storage[0], current, all_subs, 0, storage, n, alle), all_subs);
    // PMI(cov); crash();
    BUG;
  } //
}


void printD(allowedD_type allowedD) {//
  bool ok = false;
  for (int i = (int) FIRST_DOMAIN; i<=(int) LAST_DOMAINUSER; i++) {
    if (allowedD[i]) {
      ok = true;
      PRINTF("%s, ", DOMAIN_NAMES[i]);
    }
  }
  if (!ok) { PRINTF("no domains or all!"); }
  PRINTF("\n");
}

void printD(model *cov) {//
  PRINTF("'%s' allows ", NAME(cov));
  printD(cov->allowedD); //
}

void printI(allowedI_type allowedI) { //
  bool ok = false;
  for (int i = (int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; i++)
    if (allowedI[i]){
      ok = true;
      PRINTF("%s, ", ISO_NAMES[i]);
    }
  if (!ok) { PRINTF("no isotropies or all!"); }
  PRINTF("\n");
}

void printI(model *cov) { //
  PRINTF("'%s' allows ", NAME(cov));
  printI(cov->allowedI); //
}
