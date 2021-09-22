/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields 

 Copyright (C) 2001 -- 2017 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
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

#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "primitive.h"
#include "operator.h"
#include "rf_interfaces.h"
#include "startGetNset.h"


#define PERR(X) {LOCAL_MSG; SPRINTF(MSG, "'%.100s': %.800s", param_name, X); RFERROR(MSG);}
#define PERR1(X,Y) {LOCAL_MSG; LOCAL_ERRMSG2; SPRINTF(MSG, "'%.100s': %.800s", param_name, X); SPRINTF(MSG2, MSG, Y); RFERROR(MSG2);}


SEXP getListElement(SEXP list, char *str) {
  SEXP names, elmt = R_NilValue;
  PROTECT(names = getAttrib(list, R_NamesSymbol));
  int i;
  if (names == R_NilValue) {
    UNPROTECT(1);
    return R_NilValue;
  }
  for (i = 0; i < length(names); i++) {
    if(STRCMP(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  UNPROTECT(1);
  return elmt;
}

#define RET(X) {UNPROTECT(1); return X;}
int getListEltNr(SEXP list,const char *str) {
  // -1 NOMATCHING if no matching name is found
  // -2 MULTIPLEMATCHING if multiple matching fctns are found,
  //    without one matching exactly
  SEXP names;
  PROTECT(names = getAttrib(list, R_NamesSymbol));
  if (names == R_NilValue) RET(NOMATCHING);
  unsigned int ln;
  int Nr=0, i,
    n=length(names);
  
  ln=STRLEN(str);
  while ( Nr < n  && STRNCMP(str, CHAR(STRING_ELT(names, Nr)), ln)) Nr++;
  if (Nr < n) { 
    if (ln==STRLEN(CHAR(STRING_ELT(names, Nr)))) {
      for (i=Nr+1; i < n; i++) {
	if (STRNCMP(str, CHAR(STRING_ELT(names, i)), ln) == 0)
	  RET(MULTIPLEMATCHING);
      }
      RET(Nr);
    }
    // a matching function is found. Are there other functions that match?
    int j; 
    bool multiplematching=false;
    j=Nr+1; // if two or more covariance functions have the same name 
    //            the last one is taken 
    while (j<n) {
      while ( (j<n) && STRNCMP(str, CHAR(STRING_ELT(names, j)), ln)) {j++;}
      if (j<n) {
	if (ln==STRLEN(CHAR(STRING_ELT(names, j)))) { 
	  for (; j < n; j++)
	    if (STRNCMP(str, CHAR(STRING_ELT(names, j)), ln) == 0) {
	      RET(MULTIPLEMATCHING);
	    }
	  RET(j);
	}
	else {multiplematching=true;}
      }
      j++;
    }
    if (multiplematching) RET(MULTIPLEMATCHING);
  } else RET(NOMATCHING);
  RET(Nr);
}



void fetchParam(model *cov, model *next, int i, char *name) {
  if (!PARAMisNULL(next, i)) {
    if (next->ncol[i] != 1 || next->nrow[i] != 1) {
      ERR1("%.50s is not a scalar", name);
    }
    if (PisNULL(i)) kdefault(cov, i, PARAM0(next, i));
    else P(i)[0] = P0(i) * PARAM0(next, i);
  } 
}

void includeparam(void **qq,       // px
		  SEXPTYPE type,   // internal type
		  int len,         // future internal length
		  SEXP p,          // user's values
		  int base,     
		  char *param_name // name used in error messages
		  ) {
  int j;
  switch(type) {
  case REALSXP : 
    {
      *qq = MALLOC(sizeof(double) * len);      
      double *q = (double *) *qq;   
      for (j=0; j<len; j++) {
	q[j] = Real(p, param_name, base + j); // p[base + j]
	if (!GLOBAL_UTILS->basic.skipchecks && R_finite(q[j]) && 
	    FABS(q[j]) > MAXACCEPTED) 
	  RFERROR2("'%.50s' has an absolute value larger than %10e, what is believed to be a misspecification.", param_name, MAXACCEPTED);
	// insures that there is not conflict with PrepareModel2
      } 
    }
    break;
  case INTSXP :
    {      
      *qq = MALLOC(sizeof(int) * len);
      int * q = (int *) *qq;
      for (j=0; j<len; j++) q[j] = Integer(p, param_name, base + j);
    }
    break;
 
  case STRSXP :
    {
      int nch;
      *qq = MALLOC(sizeof(char*) * len);
      char** q = (char **) *qq;
      for (j=0; j<len; j++) {
	nch = STRLEN((char*) CHAR(STRING_ELT(p, j)));
	q[j] = (char *) MALLOC(sizeof(char) * (nch+1));
	STRCPY(q[j], (char*) CHAR(STRING_ELT(p, j)));
      }
    }
    break;
    /*
  case LANGSXP : case ENVSXP :   
     if (STRCMP("setseed", param_name) != 0 && STRCMP("env", param_name)!=0){ 
       if (GLOBAL.general.storing) {
	 RFERROR1("If models with R commands in the parameters (such as '%.50s') are used then 'storing' must be FALSE.", DefList[USER].nick);
       }
       if (!GLOBAL.internal.warn_oldstyle) {
	 RFERROR1("Models with R commands in the parameters (such as '%.50s') may not be called by obsolete functions.\nSee the notes in '?RMmodelsAdvanced' and set 'RFoldstyle(FALSE)'.", DefList[USER].nick);
       }
     }
     // no break! 
     
     //case CLOSXP :  
  case VECSXP : // wird kopiert somit warnung nicht notwendig 
    if ((int) TYPEOF(p) != (int) type) 
	RFERROR2("argument has type #%d whilst #%d was expected", TYPEOF(p), type);

    *qq = MALLOC(sizeof(sexp_type));
    sexp_type *q;
    q = (sexp_type*) *qq;
    q->Delete = true;
    
    q->sexp = p;  
    R_PreserveObject(q->sexp);
  
     break;
     */
  case LANGSXP : {
   if (STRCMP("setseed", param_name) != 0 && STRCMP("env", param_name)!=0){ 
       if (GLOBAL.general.storing) {
	 RFERROR1("If models with R commands in the parameters (such as '%.50s') are used then 'storing' must be FALSE.", DefList[USER].nick);
       }
       if (!GLOBAL.internal.warn_oldstyle) {
	 RFERROR1("Models with R commands in the parameters (such as '%.50s') may not be called by obsolete functions.\nSee the notes in '?RMmodelsAdvanced' and set 'RFoldstyle(FALSE)'.", DefList[USER].nick);
       }
     }
     // no break! 
   if ((int) TYPEOF(p) != (int) type) 
     RFERROR3("argument '%.50s' has type #%d whilst #%d was expected.",
	  param_name, TYPEOF(p), type);

    *qq = MALLOC(sizeof(sexp_type));
    sexp_type *q;
    q = (sexp_type*) *qq;
    q->Delete = true;
    
    q->sexp = p;  
    R_PreserveObject(q->sexp);
  }
    break;

  case ENVSXP :   {
     if (STRCMP("setseed", param_name) != 0 && STRCMP("env", param_name)!=0){ 
       if (GLOBAL.general.storing) {
	 RFERROR1("If models with R commands in the parameters (such as '%.50s') are used then 'storing' must be FALSE.", DefList[USER].nick);
       }
       if (!GLOBAL.internal.warn_oldstyle) {
	 RFERROR1("Models with R commands in the parameters (such as '%.50s') may not be called by obsolete functions.\nSee the notes in '?RMmodelsAdvanced' and set 'RFoldstyle(FALSE)'.", DefList[USER].nick);
       }
     }
     // no break! 
   if ((int) TYPEOF(p) != (int) type) 
     RFERROR3("argument '%.50s' has type #%d whilst #%d was expected.",
	  param_name, TYPEOF(p), type);

    *qq = MALLOC(sizeof(sexp_type));
    sexp_type *q;
    q = (sexp_type*) *qq;
    q->Delete = true;
    
    q->sexp = p;  
    R_PreserveObject(q->sexp);
  }
    break;

  case VECSXP : {// wird kopiert somit warnung nicht notwendig 
    if ((int) TYPEOF(p) != (int) type) 
      RFERROR3("argument '%.50s' has type #%d whilst #%d was expected.",
	  param_name, TYPEOF(p), type);

    *qq = MALLOC(sizeof(sexp_type));
    sexp_type *q;
    q = (sexp_type*) *qq;
    q->Delete = true;
    
    q->sexp = p;  
    R_PreserveObject(q->sexp);
  }
     break;


  case LISTOF + REALSXP : // list	 
    //  vector and matrix are turned into a list of 1 matrix	  
    {
      assert(base == 0);
      int i, locallen;
      SEXP pi;
      listoftype *L;
      bool list = TYPEOF(p) == VECSXP;
      if (!list && TYPEOF(p) != REALSXP && TYPEOF(p) != INTSXP &&
	  TYPEOF(p) != LGLSXP) {
	PRINTF("SXP type %d != %d\n", TYPEOF(p), REALSXP);
	BUG;
      }
      
      locallen = list ? len : 1;
     
      L = (listoftype*) (*qq = LIST_CREATE(locallen, type));

      for (i=0; i<locallen; i++) {
	pi = list ? VECTOR_ELT(p, i) : p;

	includeparam((void**) (L->lpx + i), REALSXP, length(pi), 
		     pi, base, param_name); 
       	
	if (isMatrix(pi)) {
	  // Achtung isVector ist true falls isMatrix true
	  L->ncol[i] = ncols(pi);
	  L->nrow[i] = nrows(pi);	
	} else if (isVector(pi)) {
	  L->ncol[i] = 1;
	  L->nrow[i] = length(pi);
	} else  {
	  PERR("list element(s) neither vector nor matrix");
	}
      }
    }
    break;


    /*
  case LISTOF + LISTOF + INTSXP : // list	 
    {
      BUG;

      assert(base == 0);
      int i,  k,
	locallen = len; // "data sets"
      SEXP pi;
      listoftype *q;
      if (TYPEOF(p) != VECSXP) PERR("not a list of list");
      if (locallen > MAX ELEMENTS) PERR("too many list elements");
 
      *qq = MA LLO C(sizeof(listoftype));
      q=(listoftype*) (*qq);
      q->deletelist = true;
      for (i=0; i<MAX ELEMENTS; i++) {
	q->px[i] = NULL;
	q->ncol[i] = 0;
	q->nrow[i] = 0;
      }

      for (k=0; k<locallen; k++) {  // "k-th data set"
	SEXP pk = VECTOR_ELT(p, k);
	freelistoftype *ll = 
	  (freelistoftype*) (q->px[k]=(double*) MALLOC(sizeof(freelistoftype)));
	int 
	  lpk = length(pk); // kinds of selection sets
	ll->fpx = (int **) CAL LOC(sizeof(int*), lpk);
	ll->ncol = (int *) MALLOC( sizeof(int*) * lpk); // sel sets of 
	ll->nrow = (int *) MALLOC( sizeof(int*) * lpk); // same kind
	q->nrow[k] = lpk; // 
	
	for (i=0; i<lpk; i++) {
	  pi = VECTOR_ELT(pk, i);
	  includeparam((void**) (ll->fpx + i), INTSXP, length(pi), 
		       pi, base, param_name); 
	
	  if (isMatrix(pi)) {
	    // Achtung isVector ist true falls isMatrix true
	    ll->ncol[i] = ncols(pi);
	    q->nrow[i] = nrows(pi);	
	  } else if (isVector(pi)) {
	    q->ncol[i] = 1;
	    q->nrow[i] = length(pi);
	  } else  {
	    PERR("list element(s) neither vector nor matrix");
	  }
	}
      }
    }
    break;
    */

  default : PERR("unmatched internal type of parameter");
  } // switch
}




model * CMbuild(SEXP Model, KEY_type *KT, int cR);
  

void CheckModel(SEXP Model, double *x, double *Y, double *T, 
		int Spatialdim, /* spatial dim only ! */
		int XdimOZ,
		int Lx, int Ly,
		bool Grid,
		bool Distances,
		bool Time, 
		SEXP xlist,
		KEY_type *KT,
		int cR) {
  bool dist_ok = false,
    
    distances = Distances,
    no_xlist = length(xlist) == 0; 
 
  int err,
    lx = Lx,
    zaehler = 0,
    ly = Ly;
  double *y = Y;
  model *cov = NULL;
  char EM2[LENERRMSG] = "";
    
  GetRNGstate(); // some parts use monte carlo to determine internal values
  while (true) {
    STRCPY(KT->error_loc, "Building the model");  
    cov = CMbuild(Model, KT, cR);
    //    PMI(cov);
    //printf("cov-errmd '%s'\n", cov->err_msg);

    STRCPY(KT->error_loc, "Having built the model");
    
    char *PREF_FAILURE = KEYtypeOf(cov)->PREF_FAILURE;
    STRCPY(PREF_FAILURE, "");

     //    PMI(cov);    TREE(cov);

    if (cov == NULL) RFERROR("no model is given");
    assert(cov->calling == NULL);
    
    dist_ok = COVNR == COVMATRIX || DefList[COVNR].check==check_simulate;
    dist_ok = true; // to do -- what has been the restriction??

    if (no_xlist) {
      assert(lx > 0 && ly >=0);

      if (distances) {
	if (dist_ok) {
	  lx = (int) (1e-9 + 0.5 * (1 + SQRT(1. + 8 * lx)));
	  if (Lx != lx * (lx - 1) / 2)
	    RFERROR("distance length not of form 'n * (n - 1) / 2'");
	} else {
	  NotProgrammedYet("currently: distances in simulations");
	  distances = false;
	}
      }

      cov->prevloc = LOCLIST_CREATE(1, XdimOZ + (int) Time);
      assert(cov->prevloc != NULL && cov->prevloc[0] != NULL);
      if ((err = loc_set(x, y, T, Spatialdim, XdimOZ, lx, ly, Time, Grid, 
			 distances, cov->prevloc + 0)
	   ) != NOERROR) goto ErrorHandling;
      //printf("loc cov-errmd '%s'\n", cov->err_msg);

    } else { // xlist
      assert(x == NULL && Y==NULL && T==NULL && lx == 0 && ly == 0);
      cov->prevloc = loc_set(xlist, dist_ok);
      assert(cov->prevloc != NULL);
    }
    
    location_type *loc;
    loc = Loc(cov);
    int spatialdim;
    spatialdim = loc->spatialdim;
    //  PREVDOM(0) = ly == 0 ? XONLY : KERNEL;
 
    switch (GLOBAL.coords.coord_system) {
    case coord_auto: 
    case cartesian: 
      set_system(PREV, 0, spatialdim + (int) loc->Time,
		 UNSET, loc->xdimOZ + (int) loc->Time, InterfaceType,
		 XONLY, // formal sind alle Interface Modelle nur von 
		 //                       (dummy) Variablen abhaengig
		 CARTESIAN_COORD);
      break; 
    case earth:  case sphere: {
      if (spatialdim < 2) // && spatialdim !=3) 
	GERR1("earth coordinate ans spherical coordinates consist of two angles [degree]. Got %d component.", spatialdim);
      
      set_system(PREV, 0, 2, UNSET,
		 /* 2, */ loc->spatialdim
#if MAXSYSTEMS == 1
		 + (int) loc->Time
#endif		 
		 , InterfaceType,
		 XONLY, // formal sind alle Interface Modelle nur von 
		 //                       (dummy) Variablen abhaengig
		 GLOBAL.coords.coord_system==earth ? EARTH_COORD
		 : SPHERICAL_COORD);
      
#if MAXSYSTEMS > 1
      //PMI0(cov);
      int remaining_dim = loc->xdimOZ - 2; // Hoehe auch gegeben??
      /*
	if (remaining_dim >= 1) { // surely height
	set_system(PREV, 0, 1,
	UNSET, 1, InterfaceType,
	XONLY, // formal sind alle Interface Modelle nur von 
	//                       (dummy) Variablen abhaengig
	LOGCART_COORDS); 
	remaining_dim--;
	}
      */
      
      remaining_dim += (int) loc->Time;
      if (remaining_dim >= 1) { 
	set_system(PREV, 1, remaining_dim,
		   UNSET, remaining_dim, InterfaceType,
		   XONLY, // formal sind alle Interface Modelle nur von 
		   //                       (dummy) Variablen abhaengig
		   CARTESIAN_COORD); 
      }
      //APMI0(cov);
#endif
      
      }
    
      break;
       //    case coord_auto: 
      //    case cartesian: cov->isoprev = CARTESIAN_COORD; break; 
      //    case earth: 
      //      if (spatialdim < 2) // && spatialdim !=3) 
      //	GERR1("an earth coordinate consists of two angles [degree]. Got %d component.",
      //	      spatialdim);
      //      cov->isoprev =  EARTH_COORD; 
      //      break;
      //    case sphere: 
      //      if (spatialdim < 2) // && spatialdim !=3)
      //	GERR1("a spherical coordinate consists of two angles [radians]. Got %d component.",
      //	      spatialdim);
      //      cov->isoprev =  SPHERICAL_COORD; 
      //      break;
    default: GERR("Unknown coordinate system.");
    }
    cov->calling=NULL;


    //   
 
    STRCPY(KT->error_loc, "Checking the model");
    if (PL >= PL_DETAILS) {
      //PMI(cov);//OK
    }

    assert(isInterface(cov));
    if ((err = CHECK_GEN(
			 cov, SUBMODEL_DEP, SUBMODEL_DEP, InterfaceType, true))
	!= NOERROR) {
     //     if ((err = CHECK(c o v, cov->tsdim, cov->xdimprev, InterfaceType,
     //		       cov->domprev, cov->isoprev, 
     //		       SUBMODEL_DEP, Any Type))
     //	!= NOERROR) {

      // printf("check cov-errmd '%s'\n", cov->err_msg);
    
      if (PL >= PL_ERRORS) {
	PMI(cov); // OK
	PRINTF("err =%d\n", err);
      }
      KT->error_causing_cov = cov;
      goto ErrorHandling;
    }

    if (PL >= PL_DETAILS) {
      PMI(cov); // OK
    }
 
    SPRINTF(KT->error_loc, "%.50s process", TYPE_NAMES[cov->frame]);
    //    PrInL=-1;
    
    // if (PL >= PL_DETAILS) { PRINTF("CheckModel Internal A\n"); }

    assert(CallingSet(cov));

    //    printf("\n\n\n userinterface.cc *************************************** \n\n\n\n");

    //   PMI(cov);
   
    err = STRUCT(cov, NULL);
    //    PMI0(cov);    PMI0(cov->key);
    //PMI(cov);
    //    printf("struct cov-errmd '%s'\n", cov->err_msg);

    assert((COVNR >= COVFCTN && COVNR <= VARIOGRAM_CALL) || cov->initialised);
    // printf("rrr\n");

    
    if (false) {
      PMI(cov); printf("COVNR= %d %d %d; %d  %d %d\n", //
		       COVNR, COVFCTN, VARIOGRAM_CALL,
		       cov->key==NULL,
		       cov->key==NULL ? cov->sub[0]->initialised :
		       cov->key->initialised,
		       err
		     // cov->key->initialised
		       );
    }
    assert(err != NOERROR ||
	   (COVNR >= COVFCTN && COVNR <= VARIOGRAM_CALL) ||
	   (cov->key==NULL && cov->sub[0]->initialised) ||
	   (cov->key!=NULL && cov->key->initialised));

    //assert(err == NOERROR);
 
   ErrorHandling:
    // printf("zaehler=%d %d '%s'\n", zaehler, err, cov->err_msg);
    
    if (err == NOERROR || zaehler>=1) {
      // printf("%d %d\n", err, zaehler);
      break;    
    } 

    char EM[LENERRMSG];
    model *causing = KT->error_causing_cov;
    // 94848191771984 94848191771984
    //    printf("err=%d %s %ld z=%d\n", err, NAME(causing), causing, causing->zaehler);
    ///printf("msg=%s\n", causing==NULL ?  (char*) "<strange failure>" : causing->err_msg);
    
    errorMSG(err, causing==NULL ?  (char*) "<strange failure>" : causing->err_msg,
	     KT, EM);
    SPRINTF(EM2, "%.50s: %.50s%.500s", KT->error_loc, PREF_FAILURE, EM);
    //    printf("location:%.50s  %.50s  %.50rs\n", KT->error_loc,  cov->err_msg, EM);
    if (lx == 0 || distances) break;
    y = x;
    ly = lx;
    zaehler++;

  }// Ende while

  // !!!! ACHTUG !!!!
  // DO NOT WRITE BACK SINCE, IF EVER, ONLY USED FOR TECHNICAL REASONS
  //
  PutRNGstate();// PutRNGstate muss vor UNPROTECT stehen!! -- hier irrelevant

  if (err != NOERROR) { 
    RFERROR(EM2); 
  }

  if (GLOBAL.internal.warn_mathdef == True) {
    GLOBAL.internal.warn_mathdef = False;
    PRINTF("To guarantee definiteness use 'RMmodel(var=const)', and not 'const * RMmodel()'.\n");
  }

  // muss ganz zum Schluss stehen, da Fehler zu einer unvollstaendigen
  // Struktur fuehren koennen
  assert(CallingSet(cov));
  //printf("rrr AAA\n");
}


model *InitIntern(int cR, SEXP Model, SEXP x, bool NA_OK) {
  set_currentRegister(cR);
  bool listoflists = (TYPEOF(x) == VECSXP && TYPEOF(VECTOR_ELT(x, 0)) ==VECSXP);
  SEXP
    set = listoflists ? VECTOR_ELT(x, 0) : x,
    xx = VECTOR_ELT(set, XLIST_X);
  bool
    grid = LOGICAL(VECTOR_ELT(set, XLIST_GRID))[0],
    Time = LOGICAL(VECTOR_ELT(set, XLIST_TIME))[0],
    distances = LOGICAL(VECTOR_ELT(set, XLIST_DIST))[0];
  int
    xdimOZ = grid ? ncols(xx) : nrows(xx),
    spatialdim = INTEGER(VECTOR_ELT(set, XLIST_SPATIALDIM))[0];
  KEY_type *KT = KEYT();
  
  setKT_NAOK_RANGE(KT, NA_OK);
  CheckModel(Model, NULL, NULL, NULL, spatialdim, xdimOZ, 0, 0,
	     false, // any value
	     distances, Time, x, KT, cR);
  return KT->KEY[cR];
}

SEXP Init(SEXP model_reg, SEXP Model, SEXP x, SEXP NA_OK) {

  // printf("initialising RndomFields model\n");
  
  model *cov = InitIntern(INTEGER(model_reg)[0], Model, x,
			  (bool) LOGICAL(NA_OK)[0]);

  if (PL >= PL_COV_STRUCTURE) {
    PMI(cov);// OK
  }

  //  PMI(cov);

  SEXP ans;
  PROTECT (ans = allocVector(INTSXP, 2));
  INTEGER(ans)[0] = (cov)->vdim[0];
  INTEGER(ans)[1] = (cov)->vdim[1];
  UNPROTECT(1);
  //  APMI(key[INTEGER(model_reg)[0]]);
 
  return ans;
}


/*
// for testing only
SEXP EvaluateModelXX(){
  PRINTF("dummy simulation\n");
  SEXP result;
  GetRNGstate();
  PROTECT(result = allocVector(REALSXP, 10));
  for (int i=0; i<10; i++) REAL(result)[i]=UNIFORM_RANDOM;
  PutRNGstate();// PutRNGstate muss vor UNPROTECT stehen!!
  UNPROTECT(1);
  return result;
}
*/

  
SEXP EvaluateModel(SEXP X, SEXP Covnr){
  SEXP result=R_NilValue, 
    dummy=R_NilValue;
  assert(currentNrCov != UNSET);  
  int d, mem, len;
  model *cov = KEY()[INTEGER(Covnr)[0]];
  KEY_type *KT = cov->base;
   
  //   PMI(cov);
  // printf("evaluating RndomFields model: %f x %s\n", REAL(X)[0], NAME(cov));
 
  STRCPY(KT->error_loc, "");
  if (cov == NULL) RFERROR("register not initialised");
  if ( (len = cov->qlen) == 0) {
    BUG;
    if (cov->fieldreturn == wahr) RFERROR("model cannot be evaluated");
    mem = VDIM0 * VDIM1;
  } else {   
     DefList[COVNR].cov(REAL(X), cov, NULL); // nicht FCTN, COV o.ae.
     // da Trafo des ersten Wertes moeglich !
    if (len > 1 && cov->q[len-1] == 1) len--; // last is for n=#{simulations}
    for (mem=1, d=0; d<len; d++) {
      mem *= (int) cov->q[d];
    }
  }
  assert(hasInterfaceFrame(cov));

  if (len <= 0) BUG;
  if (len <= 2) {
    if (len == 1) PROTECT(result = allocVector(REALSXP, mem)); 
    else PROTECT(result = allocMatrix(REALSXP, cov->q[0], cov->q[1]));

    GetRNGstate();
    DefList[COVNR].cov(REAL(X), cov, REAL(result));
    PutRNGstate(); // PutRNGstate muss vor UNPROTECT stehen!!
    UNPROTECT(1);
  } else {
    PROTECT(dummy=allocVector(INTSXP, len));
    for (d=0; d<len; d++) INTEGER(dummy)[d] = (int) cov->q[d];
    PROTECT(result=allocArray(REALSXP, dummy));
    GetRNGstate();
    DefList[COVNR].cov(REAL(X), cov, REAL(result));
    PutRNGstate();// PutRNGstate muss vor UNPROTECT stehen!!
    UNPROTECT(2);
  }

  //  for(int i=0; i<cov->q[0] * cov->q[1]; i++) REAL(result)[i] = i;
  //  if (result != NULL) UNPROTECT(1 + (int) (len > 2));return result;
  

  //  printf("length(X)=%d dim=%d %s (%d)\n", length(X), GATTERXDIM(0), NAME(cov), COVNR);
  // PMI0(cov);
  // assert(COVNR == SIMULATE || (length(X) == GATTERXDIM(0)));

  // PMI(cov);

 return result;
}



void GetCurrentNrOfModels(int *nr) {
  *nr = currentNrCov;
}


SEXP GetAllModelNames(SEXP Newnames){
   assert(currentNrCov != UNSET);
  int i, n;
  SEXP names;
  bool newnames = LOGICAL(Newnames)[0];
  for (n=i=0; i<currentNrCov; i++) 
    if (DefList[i].name[0] != '-') n++;
  PROTECT(names = allocVector(STRSXP, n)); 
  for (n=i=0; i<currentNrCov; i++) {
    if (DefList[i].name[0] != '-') {
      SET_STRING_ELT(names, n++,
		     mkChar(newnames ? DefList[i].nick : DefList[i].name));
    }
  }
  UNPROTECT(1);
  return names;
}
  
  

void GetModelName(int *nr,char **name, char **nick){
   assert(currentNrCov != UNSET);
  if ((*nr<0) ||(*nr>=currentNrCov)) {
    strcopyN(*name,"", MAXCHAR); 
    strcopyN(*nick,"", MAXCHAR); 
    return;
  }
  strcopyN(*name, DefList[*nr].name, MAXCHAR);
  strcopyN(*nick, DefList[*nr].nick, MAXCHAR);
}

void GetNrParameters(int *covnr, int* kappas) {
   assert(currentNrCov != UNSET);
  if (*covnr<0 || *covnr>=currentNrCov) {*kappas=-999;}
  else *kappas = DefList[*covnr].kappas;
  /*  else {
    defn *C = DefList + *covnr;
    int k = C->kappas,
      kk =0;
    for (int i=0; i<k; i++) {
      p rintf("%50s %50s\n", C->kappanames[i], INTERNAL_PARAM);
      kk += (str cmp(C->kappanames[i], INTERNAL_PARAM) != 0);
    }
    *kappas = kk;
    p rintf("%d \n", kk);
}
  */

}

void GetModelNr(char **name, int *nr) {
  *nr = getmodelnr(*name);
}


SEXP GetParameterNames(SEXP nr) {
   assert(currentNrCov != UNSET);
  defn *C = DefList + INTEGER(nr)[0]; // nicht gatternr
  SEXP pnames;
  int i;
  // print("hello %20s\n", C->name);

  PROTECT(pnames = allocVector(STRSXP, C->kappas)); 
  for (i=0; i<C->kappas; i++) {
       SET_STRING_ELT(pnames, i, mkChar(C->kappanames[i]));
  }
  UNPROTECT(1);
  return(pnames);
}

SEXP GetCathegoryNames() {
  SEXP pnames;
  int i;
  PROTECT(pnames = allocVector(STRSXP, (int) OtherType + 1)); 
  for (i=0; i<=OtherType; i++) {
    SET_STRING_ELT(pnames, i, mkChar(CAT_TYPE_NAMES[i]));
  }
  UNPROTECT(1);
  return(pnames);
}

SEXP GetSubNames(SEXP nr) {
  defn *C = DefList + INTEGER(nr)[0]; // nicht gatternr
  SEXP subnames, list, subintern;
  int i, j,
    nsub =  C->maxsub;

  // parameter and submodels may have identical names
  // this means instead of a numerical parameter a submodel can be
  // given. This happens, for instance, for nonstWM
  
  
  PROTECT(list = allocVector(VECSXP, 2)); 
  PROTECT(subnames = allocVector(STRSXP, nsub)); 
  PROTECT(subintern = allocVector(INTSXP, nsub)); 
  for (j=i=0; i<C->maxsub; i++) {
    if (C->subintern[i]) {
      PRINTF("%s subintern[%d]=true\n", C->nick, i);
    }
    
    // since 17 May 2014:
    INTEGER(subintern)[i] = C->subintern[i];
    SET_STRING_ELT(subnames, j++, mkChar(C->subnames[i]));
  
 
    // formely:
    //if (!C->subintern[i]) SET_STRING_ELT(subnames, j++,
    //                                     mkChar(C->subnames[i]));
    //for (nsub=i=0; i<C->maxsub; i++) nsub += !C->subintern[i];
 }
  
  SET_VECTOR_ELT(list, 0, subnames);
  SET_VECTOR_ELT(list, 1, subintern);
  UNPROTECT(3);
  return(list);
}


SEXP GetRange(){
  assert(false); // muss neu geschrieben werden, als SEXP
  return R_NilValue;

/*
  // never change double without crosschecking with fcts in RFCovFcts.cc!
  // index is increased by one except index is the largest value possible
  defn *C;
  model cov;
  range_type *r;
  int i,j, kappas;

   assert(currentNrCov != UNSET);
  if ((*nr<0) || (*nr>=currentNrCov)) goto ErrorHandling;
  C = DefList + *nr; // nicht gatternr
  kappas = DefList[*nr].kappas;
  if (*lparam > kappas || *lrange != kappas * 4) goto ErrorHandling;
  if (*index < 0) {
    getrange->n = 1;
    for (i=0; i<lparam; i++) {
      cov.p[i] = (double *) CALL OC(sizeof(double), 1);
      cov.p[i][0] = param[i];
    }
    getrange->n = 1;
    C->range(&cov, &getrange);
  }
  if (*index >= getrange.n) goto ErrorHandling;
  r = getrange.ranges + *index;
  for (j=i=0; i<kappas; i++) {
    range[j++] = r->min[i];
    range[j++] = r->max[i];
    range[j++] = r->pmin[i];
    range[j++] = r->pmax[i];
  }
  return;

ErrorHandling : 
  for (i=0; i<*lrange; i++) range[i]=RF_NA;
 *index = - 100;
 return;
*/
}


void PMLheader(char* firstcolumn, int nick) {
  int i;
  const char header1[]=" #    cir cut int TBM spe dir seq tre ave coi hyp spe\n";
  const char header2[]=" p    cul off rin     ctr ect uen nd  rag ns  erp cif\n";
  for (i=0; i<=nick; i++) PRINTF(firstcolumn, ""); 
  PRINTF("%4s", ""); PRINTF(header1);  
  for (i=0; i<=nick; i++) PRINTF(firstcolumn, ""); 
  PRINTF("%4s", ""); PRINTF(header2);
}

void PrintModelList(int *intern, int *operat, int* Nick)
{ // used in Roger's book
    int i, k, last_method, m, OP;
//char header[]="circ cut intr TBM spec dir seq Mak ave add hyp part\n";
  char coded[6][2]={"-", "X", "+", "N", "H", "S"};
//  char typenames[4][2]={"i", "f", "s", "n"};
  char specialnames[4][2]={".", "n", "f", "?"};
  char firstcolumn[20], name[MAXCHAR];
  int maxchar=10; // <= MAXCHAR=14
  assert(MAXCHAR >= maxchar);
 int type[MAXNRCOVFCTS], op[MAXNRCOVFCTS], monotone[MAXNRCOVFCTS], 
    finite[MAXNRCOVFCTS], simpleArguments[MAXNRCOVFCTS], 
    internal[MAXNRCOVFCTS], dom[MAXNRCOVFCTS], 
    iso[MAXNRCOVFCTS], vdim[MAXNRCOVFCTS], maxdim[MAXNRCOVFCTS],
    paramtype[MAXNRCOVFCTS * MAXPARAM],
    nick = *Nick;
  
  last_method = (int) Nothing; // even not special method   
  assert(currentNrCov!=UNSET);
  if (DefList==NULL) {PRINTF("There are no functions available!\n");} 
  else {
    defn *C;
    int includevariants = false, nr;

    GetAttr(NULL, type, op, monotone, finite, simpleArguments,
	    internal, dom, iso, maxdim, vdim, &includevariants,
	    paramtype,
	    &nr);
 
    SPRINTF(firstcolumn,"%%%ds", -maxchar);
    PRINTF("\n\n");
    PRINTF("%20s      List of models\n", "");
    PRINTF("%20s      ==============\n", "");
    PRINTF("%10s[See also PrintMethodList for the names of the columns();\n", 
	   "");
    PRINTF("%10s use 'operator=TRUE' to see all available models        ]\n",// OK 
	   ""); 

    for (OP = 0; OP <= *operat; OP++) {
      C = DefList;
      PRINTF("\n\n");
      if (OP) {
	PRINTF("%4s Operators\n", "");
 	PRINTF("%4s =========\n\n", "");  
     } else {
	PRINTF("%4s Simple models\n", "");
	PRINTF("%4s =============\n\n", "");  
      }
      PMLheader(firstcolumn, nick);

      for (k=1, i=0; i<currentNrCov; i++, C++) { 
	if ( (!isPosDef((Types)(type[i])) && !isManifold((Types) (type[i]))) ||
	     op[i] != OP ||
	     (!*intern && internal[i]) )
	  continue;
	strcopyN(name, C->name, maxchar);
	if (STRNCMP(C->name, InternalName, STRLEN(InternalName)) ==0) {
	  //	  printf("%s %d\n", C->name, *intern);
	  if (*intern < 2)  continue;
	}
	PRINTF("%2d. ", k++);
	PRINTF(firstcolumn, name);
	if (nick) {
	  strcopyN(name, C->nick, maxchar);
	  PRINTF(firstcolumn, name);
	}
	//	if (C->kappas > 9) PRINTF(
	PRINTF("%2d ", C->kappas);
//	PRINTF("%s", internal[i] ? specialnames[4] :
//		 (C->maxsub==0) ? specialnames[5] : specialnames[op[i]]);
	PRINTF("%s", 
	       specialnames[isNormalMixture((monotone_type) monotone[i]) ? 1
			    : finite[i] == 1 ? 2 
			    : isManifold((Types)(type[i])) || 
			    monotone[i]<0 || finite[i] < 0 ? 3
			    : 0]	                          
	       );	  
	PRINTF(" ");
	//             above works since normal mixture cannot have finite.range 
	assert(internal[i]==0 || internal[i]==1);
	
	for (m=(int) CircEmbed; m<last_method; m++)
	  if (m != Nugget) {
	    PRINTF("%3s%s", coded[(int) C->implemented[m]], " ");
	  }
	PRINTF("\n");
      }
    }
    
    PMLheader(firstcolumn, nick);
    PRINTF("\n%4sLegend:","");
    PRINTF("\n%4s=======\n","");
    PRINTF("First row after number of parameters:\n");
    PRINTF("'%s': normal mixture model\n",
	   specialnames[1]);
    PRINTF("'%s': finite range\n", 
	   specialnames[2]);
    PRINTF("'%s': neither a normal mixture nor a finite range\n", 
	   specialnames[0]);
   PRINTF("'%s': could be a normal mixture or have a finite range\n", 
	   specialnames[3]);

    PRINTF("\nAll other rows:\n");
    PRINTF("'%s': method not available\n", 
	   coded[0]);
    PRINTF("'%s': method available for at least some parameter values\n",
	   coded[1]);
    PRINTF("'%s': integral for the covariance is evaluated only numerically\n", 	   coded[2]);
    PRINTF("\n");
  }
}

void PrintModelList() {
  int wahr = 1, Zero = 0;
    PrintModelList(&wahr, &wahr, &Zero);
}
/*
void GetModelList(int* idx, int*internal) {
  int i, j, m;
   assert(currentNrCov != UNSET); 

  if (DefList==NULL) return;
  for (j=i=0; i<currentNrCov; i++) {
    if (!*internal && DefList[i].internal) continue;
    for (m=(int) CircEmbed; m<(int) Nothing; m++) {
      idx[j++] = DefList[i].implemented[m];
     }
  }
  return;
}
*/



void GetAttr(int *Nr, int *type, int *op, int *monotone, int *finiterange, 
	     int *simpleArguments,
	     int *internal, int *dom, int *iso, int *maxdim, int *vdim,
	     int *includevariants, int *paramtype, int *n) {
#define MAXPN 10 /* only used for testing purposes */
  int nr, p, k, j;
  defn *C = DefList;
  assert((*includevariants) xor (Nr == NULL));

  for (j = nr = 0; nr<currentNrCov; nr++, C++){   
    int v, variants = *includevariants ? C->variants : 1;
    for (v=0; v<variants; v++, j++) {
      type[j] = SYSTYPE(C->systems[v], 0); 
      dom[j] = DOM(C->systems[v], 0);
      iso[j] = ISO(C->systems[v], 0);
      if (*includevariants) Nr[j] = nr;
      vdim[j] = C->vdim;
      op[j] = (int) C->maxsub > 0;   
      maxdim[j] = MAXDIM(C->systems[v], 0);
      finiterange[j] = (int) C->finiterange;
      simpleArguments[j] = true;
      for (k=0; k<C->kappas; k++) 
	if (C->kappatype[k] != INTSXP && C->kappatype[k] != REALSXP) {
	  simpleArguments[j] = false;
	  break;
	}
      monotone[j] = C->Monotone;
      internal[j] = C->internal;
      for (p=0; p<C->kappas; p++) {
	paramtype[j * MAXPARAM + p] = C->kappaParamType[p];
      }
    }
  }
  *n = j;
}




int InternalGetProcessType(model *cov) {
  int nr = COVNR;
  if (isInterface(cov)) return InternalGetProcessType(cov->sub[0]);

  switch(SYSTYPE(DefList[nr].systems[0], 0)) {
  case TcfType: case PosDefType: case VariogramType : case TrendType :
  case GaussMethodType : return GAUSSPROC;
  case BrMethodType : return BROWNRESNICKPROC;
  case ProcessType :
    if (nr == DOLLAR_PROC) return InternalGetProcessType(cov->sub[0]);
    else if (nr == PLUS_PROC || nr == MULT_PROC) //|| nr == TREND_PROC
      return GAUSSPROC;
    return COVNR;
  case ManifoldType:
    if (nr == PLUS || nr == MULT || nr ==DOLLAR || nr == POWER_DOLLAR || 
	nr == USER) return GAUSSPROC;
    else BUG;
  default :
    BUG;
  }
 BUG;
}
 
model *Build_cov(SEXP model_reg, SEXP Model) {
  return CMbuild(Model, KEYT(), INTEGER(model_reg)[0]);
}

SEXP GetProcessType(SEXP model_reg, SEXP Model) {
  model *cov = Build_cov(model_reg, Model);
  int covnr = InternalGetProcessType(cov);
  SEXP ans;
  PROTECT (ans = allocVector(STRSXP, 1));
  SET_STRING_ELT(ans, 0, mkChar(DefList[covnr].nick));
  UNPROTECT(1);
  return ans;
}





void GetModelRegister(char **name, int* nr) {
  *nr = Match(*name, REG_NAMES, MODEL_MAX+1);

  //  print("%d\n", *nr);
  
  if (*nr<0 || *nr > MODEL_MAX) RFERROR("name for model register unknown");
}
  

SEXP allintparam() {
  defn *C;
  int n, i, np;
  for (np=n=0; n<currentNrCov; n++) {
    C = DefList + n;
    //printf("\n%s ", C->nick);
    for (i=0; i<C->kappas; i++) {
      if (C->kappatype[i] == INTSXP) {
	//printf("%s ", C->kappanames[i]);
	np++;
      }
    }
  }
  SEXP x;
  PROTECT (x = allocVector(STRSXP, np));
  for (np=n=0; n<currentNrCov; n++) {
    C = DefList + n;
    for (i=0; i<C->kappas; i++) {
      if (C->kappatype[i] == INTSXP) 
	SET_STRING_ELT(x, np++, mkChar(C->kappanames[i]));
    }
  }
  UNPROTECT(1);
  return x;
}




void CMbuild(SEXP Model, int level, model **Cov,
	     model *calling, KEY_type *KT, model *root) {
  int covnr, i, j, nkappas, //, meth, // methNr,  // old
    elt=0,
    len = length(Model);
#define NLEER 80
  char leer[NLEER], name[MAXCHAR], param_name[PARAMMAXCHAR];
  errorloc_type ERR_LOC;
  //  methname[METHODMAXCHAR], 
  SEXP m, p, names;
  defn *C; 
  model *cov;
  bool subleft[MAXSUB]; 

  if (TYPEOF(Model) != VECSXP) { // not list8
    RFERROR("list expected, which defines the (covariance) model.");
  }
  PROTECT(names = getAttrib(Model, R_NamesSymbol));
  PROTECT(m = VECTOR_ELT(Model, elt++));
  if (!isString(m)) {
    RFERROR("first element of each list must be the name of a (covariance) model"); 
  }
  if (length(m) != 1) RFERROR("model name must be a single word");
  strcopyN(name, (char*) CHAR(STRING_ELT(m, 0)), MAXCHAR);

  covnr = getmodelnr(name);
  //  if (covnr == NATSC && NS) ERR("natsc model and RFparameters(PracticalRange=T RUE) may not be given at the same time");
  if (covnr == MULTIPLEMATCHING)
    RFERROR1("multiple matching of (covariance) model name '%.50s'", name)
  else if (covnr == NOMATCHING) 
    RFERROR("'%.50s' is an unknown model name.\n'RFgetModelNames()' yields the list of models", name);

  *Cov = NULL;
  addModel(Cov, covnr, calling, true);
 
  cov = *Cov;
  if (root == NULL) root = cov;
  cov->root = root;
  cov->base = KT;
  cov->user_given = ug_explicit;
  strcopyN(leer, "                                                           ", 
	   (2 * level < NLEER) ? 2 * level : NLEER);
  SPRINTF(ERR_LOC, "%.50s\n%.50s%.50s... ", "", leer, name);//
  //  SPRINTF(ERR_LOC, "%.50s\n%.50s%.50s... ", KT->error_loc, leer, name);//
  if (!isDollar(cov)) {
    model *call = cov->calling;
    if (call != NULL && !isDollar(call))
      call->user_given = ug_implicit;
  }
  
  cov->nsub = 0;
  set_nr(GATTER, UNSET);
  C=DefList + COVNR; // nicht gatter nr

  nkappas = C->kappas;
  for (i=0; i<C->maxsub; i++) subleft[i] = true;
  
  for ( ; elt < len; elt++) {
  //    if (elt == methNr) continue;  // former definition, 26.8.2011
    p = VECTOR_ELT(Model, elt);
    i = UNSET;
    strcopyN(param_name, names == R_NilValue ? "" 
	     : CHAR(STRING_ELT(names, elt)), PARAMMAXCHAR);     

   bool empty = STRCMP(param_name, "") == 0;
   if  ((!empty || C->maxsub == 0) && nkappas != 0) {
     if (!empty) {
       if ((i = Match(param_name, C->kappanames, nkappas)) < 0) {
	 i = Match(param_name, STANDARDPARAM, nkappas);
       }
     }

      if (i<0 && !STRCMP(C->kappanames[nkappas - 1], FREEVARIABLE)) {
	if ((j = Match(param_name, C->subnames, C->maxsub)) < 0) {
	  j = Match(param_name, STANDARDSUB, MAXSUB);
	}
	if (j < 0) {
	  if (cov->ownkappanames == NULL) {
	    cov->ownkappanames = (char**) CALLOC(nkappas, sizeof(char*));
	  }
	  
	  // ein FREEVARIABLE-Zeichen wird nie verwendet -- Puffer fuer einfachen Algorithmus; deshalb "i-1"
	  //	  i = nkappas - 1; 	 	  
	  //	  while(i>0 && !STRCMP(C->kappanames[i - 1], FREEVARIABLE) && 
	  //		cov->ownkappanames[i] != NULL) i--;
	  //	  if (i <= 0 || cov->ownkappanames[i-1] != NULL) 
	  //  ERR("too many free parameters");
	  i = 0; 
	  while(STRCMP(C->kappanames[i], FREEVARIABLE) && i<nkappas) i++;
	  while(i<nkappas && (!PisNULL(i) || cov->kappasub[i] != NULL)) i++;
	  if (i >= nkappas || !PisNULL(i) || cov->kappasub[i] != NULL) 
	    ERR("too many free parameters");

	  cov->ownkappanames[i] = 
	    (char*) MALLOC(sizeof(char) * (1 + STRLEN(param_name)));
	STRCPY(cov->ownkappanames[i], param_name);
	}
      }
    } else {
      i = NOMATCHING;
    }

    if (i<0 // so, no parameter name recognised
	&& isVectorList(p) && isString(VECTOR_ELT(p, 0))) {

      // submodels
      if (!empty) { // named submodel ?
	if ((i = Match(param_name, C->subnames, C->maxsub)) < 0) {
	  i = Match(param_name, STANDARDSUB, MAXSUB);
	}
	if (i < 0) {
	  if (PL >= PL_IMPORTANT) {
	    int s;
            PRINTF("Got argument name '%s' in '%s'.", 
		   param_name, NAME(cov));
	    if (C->maxsub == 0) {
	      PRINTF("The models does not have genuine submodels.");
	    } else {
	      PRINTF("Allowed submodel names are: ");
              for (s=0; s<C->maxsub; s++) {
		if (s>0) { PRINTF(", "); }
		PRINTF("%s", STANDARDSUB[s]);
		if (STRCMP(STANDARDSUB[s], C->subnames[s]) != 0) {
		  PRINTF(", %s", C->subnames[s]);
		}
	      }
	    }
	    if (C->kappas > 0) {
	      PRINTF("Some of the genuine arguments might also suit: ");
	      for (s=0; s<C->kappas; s++) {
		if (s>0) { PRINTF(", "); }
		PRINTF("%s", C->kappanames[s]);
	      }
	    }
	    PRINTF("\n");
	  }
	  if (i == MULTIPLEMATCHING) {

#define MMERR(X) {ERR3("%.50s\n%.50s: %.50s", KT->error_loc, param_name, X);}

            MMERR("multiple matching of submodel names") 
	      } else {
            MMERR("unmatched covariance name for submodel")
	      }
	}
      } else { // unnamed submodel
	for (j=0; j < C->maxsub; j++) {
	  if (subleft[j]) {
	    i = j;
	    break;
	  }
	}
	if (j == C->maxsub) {
	  ERR("too many submodels");
	}
	// internal 
	if (C->subintern[i]) {
	  char info[200];
	  SPRINTF(info,
		  "submodel could not be identified. Please give the names of the parameters and the submodels explicitely. Did you mean '%.50s' ?", 
		  C->subnames[i]);
	  MMERR(info);
	}
      }

      if (!subleft[i]) MMERR("submodel given twice"); 
      
      subleft[i] = false;
      CMbuild(p, level + 1, cov->sub + i, cov, KT, root);
      if (isInterface(cov->sub[i]))
	ERR("models of interface type allowed only as mail model");
      //cov->sub[i]->calling = cov; // ja nicht SET_CALLING
      assert(cov->calling == NULL || cov->calling->root == cov->root);

      STRCPY(KT->error_loc, ERR_LOC);
      (cov->nsub)++;
    } else { // parameter (not list)
      // parameter identification
      int len_p = length(p); // originally "l"
      //      if (param_name[0]==ONEARGUMENT_NAME && STRLEN(param_name) == 1) { 
      //     if (STRCMP(param_name, "k") == 0) {
      if (param_name[0] == ONEARGUMENT_NAME && param_name[1] ==  '\0') {
	if (TYPEOF(p) != REALSXP && TYPEOF(p) != INTSXP && TYPEOF(p) != LGLSXP) 
	  PERR1("if '%c' is used as parameter name then only a numerical vector is allowed as value", ONEARGUMENT_NAME);
	if (len_p != C->kappas) 
	  PERR1("length of vector does not match number of parameters. Do not use short name '%c' in complex situtations either.", ONEARGUMENT_NAME);
	// if (C->maxsub!=0)
	// ERR("short form k of parameter name not allowed for sophistacted models");
	for (j=0; j<len_p; j++) {
	  if (!PisNULL(j)) PERR("parameter given twice"); // p[0] OK
	  if (true ||  // neu seit 6.8.14, wegen RFgui, modelParam = .Call(C_SetAndGetModelInfo",
	      C->kappatype[j] != INTSXP ||
	      Integer(p, param_name, j) != NA_INTEGER) {
	    cov->ncol[j] = cov->nrow[j] = 1;	  

	    if ( C->kappatype[j] != INTSXP && C->kappatype[j] != REALSXP)
	      PERR("model has too complicated arguments to allow abbreviations");
	    
	    includeparam((void**) (cov->px + j), C->kappatype[j], 1, 
			 p, j, param_name);
	  }
	}
	continue;
      }

      if (STRCMP(param_name, "") == 0 && i < 0) {
	if (nkappas == 1) i = 0; else ERR("Parameters must be named");
      } else {
	// i set before the submodels
	if (i < 0) {
	  if (PL >= PL_IMPORTANT) {
	    int s;
            PRINTF("allowed parameter names are: ");
	    if (C->kappas == 0) { PRINTF("<none>"); }
	    else {
	      for (s=0; s<C->kappas; s++) {
		if (STRCMP("", C->kappanames[s]) != 0) {		  
		  PRINTF("'%s', ", C->kappanames[s]);
		}
	      }
	      PRINTF("alternatively ");
	      for (s=0; s<C->kappas; s++) PRINTF("'%s', ", STANDARDPARAM[s]);
	    }
            PRINTF("allowed sub names are: ");
	    if (C->maxsub == 0) { PRINTF("<none>"); }
	    else for (s=0; s<C->maxsub; s++) {
		PRINTF("%s, ", C->subnames[s]);
	      }
	    PRINTF("\n");
	  }

	  if (i == MULTIPLEMATCHING) PERR("multiple matching of parameter");
	  
	  if (i == NOMATCHING){
	    PERR("recognized as parameter, but unknown");      
	  } else {
	    BUG;
	  }
	}
      }
  
      if (isVectorList(p) && isString(VECTOR_ELT(p, 0))) {
	PFREE(i);
	CMbuild(p, level + 1, cov->kappasub + i, cov, KT, root);
	if (isInterface(cov->kappasub[i]))
	  ERR("models of interface type allowed only as mail model");
	// cov->kappasub[i]->calling =  cov; // ja nicht SET_CALLING
 	assert(cov->calling == NULL || cov->calling->root == cov->root);
	STRCPY(KT->error_loc, ERR_LOC);
      } else {
	// first the fixed one for the dimensions !!
	if (!PisNULL(i) || cov->kappasub[i] != NULL)
	  ERR("parameter given twice");

	includeparam((void**) (cov->px + i), C->kappatype[i], len_p, 
		     p, 0, param_name);
	
	if (C->kappatype[i] > LISTOF && TYPEOF(p) != VECSXP) {
	  if (isMatrix(p) || isVector(p)) {
	    cov->ncol[i] = cov->nrow[i] = 1;
	  } else {
	    PERR("not a list of vectors or matrices");
	  }
	  
	} else {	  
	  if (isMatrix(p)) {
	    // Achtung isVector ist true falls isMatrix true
	    cov->ncol[i] = ncols(p);
	    cov->nrow[i] = nrows(p);
	  } else if (isVector(p)) {
	    cov->ncol[i] = 1;
	    cov->nrow[i] = len_p;
	  } else if (isFunction(p) || isLanguage(p) || isEnvironment(p)) {
	    cov->ncol[i] = 1;
	    cov->nrow[i] = 1;	  
	  } else  {
	    PERR("neither vector nor matrix");
	  }
	}
      }
    } // else
  } // for elt

  for (j=0; j < C->maxsub; j++) if (subleft[j]) break;
  if (j < C->maxsub) {
    bool shape_process = // covnr == SCHLATHER || covnr == SMITH ||
      covnr == RANDOMCOIN_USER || COVNR == AVERAGE_USER;
    if (j < C->minsub) { 
      if (!shape_process || subleft[C->maxsub - 1]) {
	model *call = cov->calling;
	while (call != NULL && MODELNR(call) != RFGET) //not CALLNR
	  call = call->calling;
	if (call == NULL) ERR("not enough submodels");      
      }
      if (PL >= PL_IMPORTANT && !shape_process)  {
	// bei shape_process darf von 2 untermodellen nur 1 angegeben werden
	warning("not all submodels are given");
      }
    }
  }

  if (GLOBAL.general.naturalscaling > 0) {
    for (j=0; j < C->maxsub; j++) {
      model *sub = cov->sub[0];
      if (sub != NULL) {
	defn *S = DefList + SUBNR;
	if (S->maxsub == 0 && S->inverse != ErrInverse && S->cov != nugget &&
	    S->cov != constant) {
	  if (GLOBAL.general.naturalscaling == NATSCALE_MLE && isDollar(cov)) {
	    // natsc wird nur dann eingeschoben, wenn davor scale=NA war
	    bool natsc= !PisNULL(DSCALE) && ISNAN(P0(DSCALE)) && 
	      cov->kappasub[DSCALE] == NULL;
	
	    int idx = DANISO;
	    while (!natsc) {
	      int ii,
		total = cov->ncol[idx] * cov->nrow[idx];
	      double *pp = P(idx);
	      for (ii=0; ii<total; ii++) {
		if (ISNAN(pp[ii])) {
		  natsc = true;
		  break;
		}
	      }
	      if (idx == DAUSER || natsc) break;
	      idx = DAUSER;
	    }
	    if (natsc) addModel(cov, j, NATSC_INTERN); 
	  } else {
	    addModel(cov, j, NATSC_USER);
	  }
	}
      }
    }
  }
  UNPROTECT(2);
} // CMbuild

model *CMbuild(SEXP Model, KEY_type *KT, int cR) {
  assert(currentNrCov != UNSET);
  assert(KT != NULL && KT->KEY != NULL)  
  if (cR < 0 || cR >= MODEL_MAX) BUG;
  model **Cov = KT->KEY + cR;
  if (*Cov != NULL) {
    COV_DELETE(Cov, NULL);
  }
  assert(*Cov == NULL);
 
  CMbuild(Model, 0, Cov, NULL, KT, NULL);
  model *cov = *Cov;
  assert(cov->root == cov);
  assert(cov->base == KT);
  if (!isInterface(cov)) {
    //PMI0(cov); crash();
    BUG;
  }

  //printf("cmbuild %d\n", cR);
  // if (cR == 13)
  // APMI(*Cov);
  
  return cov;
}
