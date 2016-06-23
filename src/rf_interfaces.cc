/*
Authors
Martin Schlather, schlather@math.uni-mannheim.de

main library for unconditional simulation of random fields

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2015 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Fre Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/* 

  PRINTING LEVELS
  ===============

  error messages: 1
  forcation : 2
  minor tracing information : 3--5
  large debugging information: >10

*/

/*
  calculate the transformation of the points only once and store the result in
  a register (indexed by ncov) of key, not of s->. Advantages: points have to
  calculated only once for several tries; if zonal anisotropy then specific 
  methods may be called;  
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include "RF.h"
//#include <unistd.h>
#include "Operator.h"
#include "Coordinate_systems.h"
#include "variogramAndCo.h"


/* 
 in CheckCovariance and other the following dimensions are used:
 xdimOZ       : dimension of the points given by the user.
              value is <= spacedim since in case of isotropy only
	      the distances might be given by the user
 timespacedim : (= kc->dim = key->timespacedim)
              the true dimension for the location (if necessary by 
	      explicite parameter, e.g. in CovarianceFct.)
*/



SEXP getListElement(SEXP list, char *str) {
  SEXP elmt = R_NilValue,
    names = getAttrib(list, R_NamesSymbol);
  int i;
  if (names == R_NilValue) {
    return R_NilValue;
  }
  for (i = 0; i < length(names); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  return elmt;
}

int getListEltNr(SEXP list,const char *str) {
  // == -1 if no matching name is found
  // == -2 if multiple matching fctns are found, without one matching exactly
  SEXP names = getAttrib(list, R_NamesSymbol);
  if (names == R_NilValue) return NOMATCHING;
  unsigned int ln;
  int Nr=0, i,
    n=length(names);
  
  ln=strlen(str);
  while ( Nr < n  && strncmp(str, CHAR(STRING_ELT(names, Nr)), ln)) Nr++;
  if (Nr < n) { 
    if (ln==strlen(CHAR(STRING_ELT(names, Nr)))) {
      for (i=Nr+1; i < n; i++) {
	if (strncmp(str, CHAR(STRING_ELT(names, i)), ln) == 0) {
	  return MULTIPLEMATCHING;
	}
      }
      return Nr;
    }
    // a matching function is found. Are there other functions that match?
    int j; 
    bool multiplematching=false;
    j=Nr+1; // if two or more covariance functions have the same name 
    //            the last one is taken 
    while (j<n) {
      while ( (j<n) && strncmp(str, CHAR(STRING_ELT(names, j)), ln)) {j++;}
      if (j<n) {
	if (ln==strlen(CHAR(STRING_ELT(names, j)))) { 
	  for (; j < n; j++)
	    if (strncmp(str, CHAR(STRING_ELT(names, j)), ln) == 0) {
	      return MULTIPLEMATCHING;
	    }
	  return j;
	}
	else {multiplematching=true;}
      }
      j++;
    }
    if (multiplematching) {
      return MULTIPLEMATCHING;
    }
  } else return NOMATCHING;
  return Nr;
}



void fetchParam(cov_model *cov, cov_model *next, int i, char *name) {
  if (!PARAMisNULL(next, i)) {
    if (next->ncol[i] != 1 || next->nrow[i] != 1) {
      char msg[255];
      sprintf(msg, "%s is not a scalar", name);
      ERR(msg);
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
	nch = strlen((char*) CHAR(STRING_ELT(p, j)));
	q[j] = (char *) MALLOC(sizeof(char) * (nch+1));
	strcpy(q[j], (char*) CHAR(STRING_ELT(p, j)));
      }
    }
    break;

    /*

  case LANGSXP : case ENVSXP :   
     if (strcmp("setseed", param_name) != 0 && strcmp("env", param_name)!=0){ 
       if (GLOBAL.general.storing) {
	 ERR1("If models with R commands in the parameters (such as '%s') are used then 'storing' must be FALSE.", CovList[USER].nick);
       }
       if (!GLOBAL.internal.warn_oldstyle) {
	 ERR1("Models with R commands in the parameters (such as '%s') may not be called by obsolete functions.\nSee the notes in '?RMmodelsAdvanced' and set 'RFoldstyle(FALSE)'.", CovList[USER].nick);
       }
     }
     // no break! 
     
     //case CLOSXP :  
  case VECSXP : // wird kopiert somit warnung nicht notwendig 
    if ((int) TYPEOF(p) != (int) type) 
	ERR2("argument has type #%d whilst #%d was expected", TYPEOF(p), type);

    *qq = MALLOC(sizeof(sexp_type));
    sexp_type *q;
    q = (sexp_type*) *qq;
    q->Delete = true;
    
    q->sexp = p;  
    R_PreserveObject(q->sexp);

     break;
     */
  case LANGSXP : {
   if (strcmp("setseed", param_name) != 0 && strcmp("env", param_name)!=0){ 
       if (GLOBAL.general.storing) {
	 ERR1("If models with R commands in the parameters (such as '%s') are used then 'storing' must be FALSE.", CovList[USER].nick);
       }
       if (!GLOBAL.internal.warn_oldstyle) {
	 ERR1("Models with R commands in the parameters (such as '%s') may not be called by obsolete functions.\nSee the notes in '?RMmodelsAdvanced' and set 'RFoldstyle(FALSE)'.", CovList[USER].nick);
       }
     }
     // no break! 
   if ((int) TYPEOF(p) != (int) type) 
	ERR2("argument has type #%d whilst #%d was expected", TYPEOF(p), type);

    *qq = MALLOC(sizeof(sexp_type));
    sexp_type *q;
    q = (sexp_type*) *qq;
    q->Delete = true;
    
    q->sexp = p;  
    R_PreserveObject(q->sexp);
  }
    break;

  case ENVSXP :   {
     if (strcmp("setseed", param_name) != 0 && strcmp("env", param_name)!=0){ 
       if (GLOBAL.general.storing) {
	 ERR1("If models with R commands in the parameters (such as '%s') are used then 'storing' must be FALSE.", CovList[USER].nick);
       }
       if (!GLOBAL.internal.warn_oldstyle) {
	 ERR1("Models with R commands in the parameters (such as '%s') may not be called by obsolete functions.\nSee the notes in '?RMmodelsAdvanced' and set 'RFoldstyle(FALSE)'.", CovList[USER].nick);
       }
     }
     // no break! 
   if ((int) TYPEOF(p) != (int) type) 
	ERR2("argument has type #%d whilst #%d was expected", TYPEOF(p), type);

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
	ERR2("argument has type #%d whilst #%d was expected", TYPEOF(p), type);

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
	PRINTF("type %d != %d", TYPEOF(p), REALSXP);
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
	ll->fpx = (int **) CALLOC(sizeof(int*), lpk);
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



void CMbuild(SEXP model, int level, cov_model **Cov) {
  int covnr, i, j, nkappas, //, meth, // methNr,  // old
    elt=0,
    len = length(model);
#define NLEER 80
  char leer[NLEER], name[MAXCHAR], param_name[PARAMMAXCHAR], ERR_LOC[nErrorLoc];
  //  methname[METHODMAXCHAR], 
  // msg[200];
  SEXP m, p,
    names = getAttrib(model, R_NamesSymbol);
  cov_fct *C; 
  cov_model *cov;
  bool subleft[MAXSUB]; 


  if (TYPEOF(model) != VECSXP) { // not list8
    ERR("list expected")
      }
  PROTECT(m = VECTOR_ELT(model, elt++));
  if (!isString(m)) {
    ERR("first element must be the name of a covariance model"); 
  }
  if (length(m) != 1) ERR("model name must be a single word");
  strcopyN(name, (char*) CHAR(STRING_ELT(m, 0)), MAXCHAR);

  strcopyN(leer, "                                                           ", 
	   (2 * level < NLEER) ? 2 * level : NLEER);
  sprintf(ERR_LOC, "%s\n%s%s... ", ERROR_LOC, leer, name);
  strcpy(ERROR_LOC, ERR_LOC);
  covnr = getmodelnr(name);

  //  if (covnr == NATSC && NS) ERR("natsc model and RFparameters(PracticalRange=TRUE) may not be given at the same time");
  if (covnr == -2) { ERR("multiple matching of covariance name") }
  else if (covnr == -1) { 
    ERR("unknown model.  'RFgetModelNames()' yields the list of models");
  }

  *Cov = NULL;
  addModel(Cov, covnr, NULL, true);
  
  cov = *Cov;
  cov->user_given = ug_explicit;
  if (!isDollar(cov)) {
    cov_model *calling = cov->calling;
    if (calling != NULL && !isDollar(calling))
      calling->user_given = ug_implicit;
  }
  cov->nsub = 0;
  cov->gatternr = S2S;
  C=CovList + cov->nr; // nicht gatternr

  nkappas = C->kappas;
  for (i=0; i<C->maxsub; i++) subleft[i] = true;
  
  for ( ; elt < len; elt++) {
  //    if (elt == methNr) continue;  // former definition, 26.8.2011
    p = VECTOR_ELT(model, elt);
    strcopyN(param_name, names == R_NilValue ? "" 
	     : CHAR(STRING_ELT(names, elt)), PARAMMAXCHAR);     

    if  (strcmp(param_name, "") && nkappas != 0) {
      if ((i = Match(param_name, C->kappanames, nkappas)) < 0) {
	i = Match(param_name, STANDARDPARAM, nkappas);
      }

      if (i<0 && !strcmp(C->kappanames[nkappas - 1], FREEVARIABLE)) {
	if ((j = Match(param_name, C->subnames, C->maxsub)) < 0) {
	  j = Match(param_name, STANDARDSUB, MAXSUB);
	}
	if (j < 0) {
	  if (cov->ownkappanames == NULL) {
	    cov->ownkappanames = (char**) CALLOC(nkappas, sizeof(char*));
	  }
	  
	  // ein FREEVARIABLE-Zeichen wird nie verwendet -- Puffer fuer einfachen Algorithmus; deshalb "i-1"
	  i = nkappas - 1; 
	  	 	  
	  while(i>0 && !strcmp(C->kappanames[i - 1], FREEVARIABLE) && 
		cov->ownkappanames[i] != NULL) {
	    i--;
	  }
	  if (i <= 0 || cov->ownkappanames[i-1] != NULL) 
	    ERR("too many free parameters");
	  cov->ownkappanames[i] = 
	    (char*) MALLOC(sizeof(char) * (1 + strlen(param_name)));
	strcpy(cov->ownkappanames[i], param_name);
	}
      }
    } else {
      i = NOMATCHING;
    }

    if (i<0 // so, no parameter name recognised
	&& isVectorList(p) && isString(VECTOR_ELT(p, 0))) {

      // submodels
      if (strcmp(param_name, "") != 0) { // named submodel
	if ((i = Match(param_name, C->subnames, C->maxsub)) < 0) {
	  i = Match(param_name, STANDARDSUB, MAXSUB);
	}
	if (i < 0) {
	  if (PL >= PL_IMPORTANT) {
	    int s;
            PRINTF("allowed submodel names are: ");
	    if (C->maxsub == 0) PRINTF("(none)"); 
	    else   
              for (s=0; s<C->maxsub; s++) {
		if (s>0)  PRINTF(", ");
		PRINTF("%s", STANDARDSUB[s]);
		if (strcmp(STANDARDSUB[s], C->subnames[s]) != 0) {
		  PRINTF("%s", C->subnames[s]);
		}
	      }
	    PRINTF("\n");
	  }
	  if (i == MULTIPLEMATCHING) {

#define MMERR(X) {ERR3("%s\n%s: %s", ERROR_LOC, param_name, X);}

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
	  char msg[200];
	  sprintf(msg,
		  "submodel could not be identified. Please give the names of the parameters and the submodels explicitely. Did you mean '%s' ?", 
		  C->subnames[i]);
	  MMERR(msg);
	}
      }

      if (!subleft[i]) MMERR("submodel given twice"); 
      
      subleft[i] = false;
      CMbuild(p, level + 1, cov->sub + i);
      strcpy(ERROR_LOC, ERR_LOC);
      cov->sub[i]->calling = cov;
      (cov->nsub)++;
    } else { // parameter (not list)
      // parameter identification
      int len_p = length(p); // originally "l"
      //      if (param_name[0]==ONEARGUMENT_NAME && strlen(param_name) == 1) { 
      //     if (strcmp(param_name, "k") == 0) {
      if (param_name[0] == ONEARGUMENT_NAME && param_name[1] ==  '\0') {
	if (TYPEOF(p) != REALSXP && TYPEOF(p) != INTSXP && TYPEOF(p) != LGLSXP) 
	  PERR1("if '%c' is used as parameter name then only a numerical vector is allowed as value", ONEARGUMENT_NAME);
	if (len_p != C->kappas) 
	  PERR1("length of vector does not match number of parameters. Do not use short name '%c' in complex situtations either.", ONEARGUMENT_NAME);
	// if (C->maxsub!=0)
	// ERR("short form k of parameter name not allowed for sophistacted models");
	for (j=0; j<len_p; j++) {
	  if (!PisNULL(j)) PERR("parameter given twice"); // p[0] OK
	  if (true ||  // neu seit 6.8.14, wegen RFgui, modelParam = .Call("SetAndGetModelInfo",
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
      
      if (strcmp(param_name, "") == 0) {
	if (nkappas == 1) i = 0; else ERR("parameters must be named");
      } else {
	// i set before the submodels
	if (i < 0) {
	  if (PL >= PL_IMPORTANT) {
	    int s;
            PRINTF("allowed parameter names are: ");
	    if (C->kappas == 0) PRINTF("<none>"); 
	    else for (s=0; s<C->kappas; s++) {
		if (strcmp("", C->kappanames[s]) != 0) 
		  PRINTF("'%s', ", C->kappanames[s]);
		PRINTF("'%s', ", STANDARDPARAM[s]);
	      }
            PRINTF("allowed sub names are: ");
	    if (C->maxsub == 0) PRINTF("<none>"); 
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
	CMbuild(p, level + 1, cov->kappasub + i);
	strcpy(ERROR_LOC, ERR_LOC);
	cov->kappasub[i]->calling = cov;
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
  }

  for (j=0; j < C->maxsub; j++) if (subleft[j]) break;
  if (j < C->maxsub) {
    bool shape_process = isMaxStable(cov) || 
      covnr == RANDOMCOIN_USER || cov->nr == AVERAGE_USER;
    if (j < C->minsub) { 
      if (!shape_process || subleft[C->maxsub - 1]) {
	cov_model *prev = cov->calling;
	while (prev != NULL && prev->nr != RFGET) prev = prev->calling;
	if (prev == NULL) ERR("not enough submodels");      
      }
      if (PL >= PL_IMPORTANT && !shape_process)  {
	// bei shape_process darf von 2 untermodellen nur 1 angegeben werden
	warning("not all submodels are given");
      }
    }
  }
  
  cov_model *next = cov->sub[0];
  if (next != NULL) {
    cov_fct *N = CovList + next->nr;
    if (NS > 0 && N->maxsub == 0 && N->inverse != ErrInverse &&
	N->cov != nugget &&  N->cov != constant) {
      if (NS == NATSCALE_MLE && isDollar(cov)) {
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
	if (natsc) addModel(cov, 0, NATSC_INTERN); 
      } else {
	addModel(Cov, NATSC_USER); 
      }
    }
  }
  UNPROTECT(1);
}


bool CallingSet(cov_model *cov) {
  int i;
  for (i=0; i<cov->nsub; i++) {
    cov_model *sub = cov->sub[i];
    if (sub == NULL) {
      if (CovList[cov->nr].range != range_randomcoin) {
	PMI(cov);//
	return false;
      }
    } else {
      if (sub->calling != cov) {
	PRINTF("%dth submodel\n", i);
	PMI(cov);//
	return false;
      }
      if (!CallingSet(sub)) {  
	return false;
      }
    }
  }
    for (; i < MAXSUB; i++) {
      if (cov->sub[i] != NULL) {
	char msg[200];
	sprintf(msg, "%s: %dth submodel not NULL although nsub=%d",
		NAME(cov), i, cov->nsub);
	warning(msg);
	BUG;
      }
    }

 if (cov->key != NULL && !CallingSet(cov->key)) return false;
  if (cov->Splus != NULL) {
    for (i=0; i<cov->nsub; i++) {
      cov_model *sub = cov->Splus->keys[i];
      if (sub == NULL) {
	PMI(cov); //
	return false;
      } else {
	if (sub->calling != cov) {
	  PMI(cov); //
	  return false;
	}
	if (!CallingSet(sub)) return false;
      }
    }
  }
  return true;
}

void CheckModelInternal(SEXP model, double *x, double *Y, double *T, 
			int Spatialdim, /* spatial dim only ! */
			int XdimOZ,
			int Lx, int Ly,
			bool Grid,
			bool Distances,
			bool Time, 
			SEXP xlist,
			cov_model **Cov) {
  bool dist_ok = false,
    distances = Distances,
    no_xlist = length(xlist) == 0; 
 
  int err,
    lx = Lx,
    zaehler = 0,
    ly = Ly;
  double *y = Y;
  cov_model *cov;
  char EM2[LENERRMSG] = "";
   
  GetRNGstate(); // some parts use monte carlo to determine internal values
  if (currentNrCov==-1) InitModelList();  

   while (true) {
    strcpy(ERROR_LOC, "Building the model:");
    cov = NULL;
   
    if (*Cov != NULL) {
      //APMI(*Cov);
      COV_DELETE(Cov);
    }
    assert(*Cov == NULL);
    CMbuild(model, 0, Cov);
 
    strcpy(ERROR_LOC, "Having built the model:");
    cov = *Cov;
    if (cov == NULL) ERR("no model is given");
    assert(cov->calling == NULL);

    dist_ok = cov->nr == COVMATRIX || CovList[cov->nr].check==check_simulate;
    dist_ok = true; // to do -- what has been the restriction??

    if (no_xlist) {
      cov->prevloc = LOCLIST_CREATE(1);
      assert(cov->prevloc != NULL && cov->prevloc[0] != NULL);

      assert(lx > 0 && ly >=0);

      if (distances) {
	if (dist_ok) {
	  lx = (int) (1e-9 + 0.5 * (1 + sqrt(1. + 8 * lx)));
	  if (Lx != lx * (lx - 1) / 2)
	    ERR("distance length not of form 'n * (n - 1) / 2'");
	} else {
	  NotProgrammedYet("currently: distances in simulations");
	  distances = false;
	}
      }
      if ((err = loc_set(x, y, T, Spatialdim, XdimOZ, lx, ly, Time, Grid, 
			 distances, cov->prevloc + 0)
	   ) != NOERROR) goto ErrorHandling;

    } else { // xlist
      // double AA=0; for (int ii=1; ii<10000000; ii++) for (int ii0=1; ii0<10000000; ii0++) AA+=exp(ii0) / log(ii);  printf("C\n");
      //PMI(cov); 

      // printf("%d %d %d %d %d\n", x==NULL, Y==NULL, T==NULL, lx==0, ly==0);

      assert(x == NULL && Y==NULL && T==NULL && lx == 0 && ly == 0);
      cov->prevloc = loc_set(xlist, dist_ok);
      assert(cov->prevloc != NULL);
    }
    
    location_type *loc;
    int spatialdim;
    loc = Loc(cov);
    cov->tsdim = loc->spatialdim + (int) loc->Time;
    cov->xdimprev = loc->xdimOZ + (int) loc->Time;
    spatialdim = loc->spatialdim;
    //  cov->domprev = ly == 0 ? XONLY : KERNEL;
    cov->domprev = XONLY; // formal sind alle Interface Modelle nur von 
    //                              einer (dummy) Variablen abhaengig

    //printf("co %d\n", GLOBAL.coords.coord_system);
    //   printf("%d %d\n", loc->xdimOZ,  cov->xdimprev);

    switch (GLOBAL.coords.coord_system) {
    case coord_auto: 
    case cartesian: cov->isoprev = CARTESIAN_COORD; break; 
    case earth: 
      if (spatialdim < 2) // && spatialdim !=3) 
	GERR1("an earth coordinate consists of two angles [degree]. Got %d component.",
	      spatialdim);
      cov->isoprev =  EARTH_COORDS; 
      break;
    case sphere: 
      if (spatialdim < 2) // && spatialdim !=3)
	GERR1("a spherical coordinate consists of two angles [radians]. Got %d component.",
	      spatialdim);
      cov->isoprev =  SPHERICAL_COORDS; 
      break;
    default: GERR("Unknown coordinate system.");
    }
    cov->calling=NULL;

 
   strcpy(ERROR_LOC, "Checking the model:");
    if (PL >= PL_DETAILS) {
      //PMI(cov);//OK
    }

    if ((err = CHECK(cov, cov->tsdim, cov->xdimprev, InterfaceType,
		       cov->domprev, cov->isoprev, 
		       SUBMODEL_DEP, ROLE_BASE))
	!= NOERROR) {
     ////
      if (PL >= PL_ERRORS) {
	PMI(cov); // OK
	PRINTF("err =%d\n", err);
      }
      goto ErrorHandling;
    }

    if (PL >= PL_DETAILS) {
      PMI(cov); // OK
    }
 
    sprintf(ERROR_LOC, "%s process: ", ROLENAMES[cov->role]);
    strcpy(PREF_FAILURE, "");
    PrInL=-1;
    
    if (PL >= PL_DETAILS) PRINTF("CheckModel Internal A\n");

    assert(CallingSet(cov));
  
    err = STRUCT(cov, NULL);

    //assert(err == NOERROR);
 
   ErrorHandling:

    if (err == NOERROR || zaehler>=1) {
      break;    
    } 

    char EM[LENERRMSG];   
    FinalErrorMSG(err, EM);
    sprintf(EM2, "%s%s", PREF_FAILURE, EM);
    if (lx == 0 || distances) break;
    y = x;
    ly = lx;
    zaehler++;

   }// Ende while
  
  PutRNGstate();


  if (err != NOERROR) { 
    if (PL >= PL_ERRORS) {ERRLINE; }
    RFERROR(EM2); 
  }
  
  // muss ganz zum Schluss stehen, da Fehler zu einer unvollstaendigen
  // Struktur fuehren koennen
  assert(CallingSet(cov));
}


SEXP Init(SEXP model_reg, SEXP model, SEXP x, SEXP NA_OK) {
  currentRegister = INTEGER(model_reg)[0];
  NAOK_RANGE = (bool) LOGICAL(NA_OK)[0];
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

  CheckModelInternal(model, NULL, NULL, NULL, spatialdim,
		     xdimOZ, 0, 0,
		     false, // any value
		     distances, Time,
		     x,
		     KEY + currentRegister
		     );

  cov_model *cov = KEY[currentRegister];
  NAOK_RANGE = false;
  if (PL >= PL_COV_STRUCTURE) {
    PMI(cov);// OK
  }

  SEXP ans;
  PROTECT (ans = allocVector(INTSXP, 2));
  INTEGER(ans)[0] = (cov)->vdim[0];
  INTEGER(ans)[1] = (cov)->vdim[1];
  UNPROTECT(1);
  //  APMI(KEY[INTEGER(model_reg)[0]]);
 
  return ans;
}



SEXP EvaluateModel(SEXP X, SEXP Covnr){
  
  // printf("AX0\n\n"); if (PL == 100) BUG;

  if (currentNrCov==-1) InitModelList();  
  int d, mem, len,
    err = NOERROR;
  cov_model *cov = KEY[INTEGER(Covnr)[0]];
  SEXP result=NULL, 
    dummy=NULL;
  
  //   PMI(cov);

  strcpy(ERROR_LOC, "");
  if (cov == NULL) GERR("register not initialised");
  if ( (len = cov->qlen) == 0) {
    BUG;
    if (cov->fieldreturn) GERR("model cannot be evaluated");
    mem = cov->vdim[0] * cov->vdim[1];
  } else {   
     CovList[cov->nr].cov(REAL(X), cov, NULL); // nicht FCTN, COV o.ae.
     // da Trafo des ersten Wertes moeglich !
    if (len > 1 && cov->q[len-1] == 1) len--; // last is for n=#{simulations}
    for (mem=1, d=0; d<len; d++) {
      mem *= (int) cov->q[d];
    }
  }

  if (len == 1) PROTECT(result = allocVector(REALSXP, mem)); 
  else if (len == 2)  
    PROTECT(result = allocMatrix(REALSXP, cov->q[0], cov->q[1])); 
  else {
    PROTECT(dummy=allocVector(INTSXP, len));
    for (d=0; d<len; d++) {
      INTEGER(dummy)[d] = (int) cov->q[d];
    }
    PROTECT(result=allocArray(REALSXP, dummy));
  }

  //  for(int i=0; i<cov->q[0] * cov->q[1]; i++) REAL(result)[i] = i;
  //  if (result != NULL) UNPROTECT(1 + (int) (len > 2));return result;
  

  GetRNGstate();
 
  FCTN(REAL(X), cov, REAL(result)); 
 
  PutRNGstate();

 ErrorHandling:
  if (result != NULL) UNPROTECT(1 + (int) (len > 2));
  if (err != NOERROR) XERR(err);

 return result;
}


#define DENS_LOG 0
#define DENS_SEED 1
#define DENS_ENV 2
/*
void density(double VARIABLE_IS_NOT_USED *value, cov_model *cov, double *v) {

  assert(!P0INT(DENS_LOG));
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  char errorloc_save[nErrorLoc];
  int ni = 0,
    err = NOERROR;
  double *res;
  location_type *loc = PrevLoc(cov);
  long vdimtot; vdimtot = loc->totalpoints * cov->vdim[0];
  assert(cov->vdim[0] == cov->vdim[1]);

  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array

  strcpy(errorloc_save, ERROR_LOC);

  PutRNGstate();
  ERR("stop : ni nae Zei falsch");
  double simu_seed = GLOBAL_UTILS->basic.seed + (ni - 1);
  addVariable((char*) "seed", &simu_seed, 1, 1, PENV(DENS_ENV)->sexp);
  eval(PLANG(DENS_SEED)->sexp, PENV(DENS_ENV)->sexp);
  GetRNGstate();

  sprintf(ERROR_LOC, "%s %d: ", errorloc_save, ni);
 
  assert(cov->Sgen != NULL);
  
  NotProgrammedYet("density");

  //  if (cov->nr == DENSITY) LOGDENSITY(sub, cov->Sgen); 
  //  else if (cov->nr == PROBAB) PROBAB(sub, cov->Sgen);

  if (sizeof(double) == sizeof(double) && false) {
    MEMCOPY(res, cov->rf, sizeof(double) * vdimtot);
  } else {
    int i; for (i=0; i<vdimtot; i++) {
      res[i] = cov->rf[i];
    }
  }
  
  if (!sub->simu.active) 
    GERR1("could not perform multiple simulations. Is '%s == FALSE'?",
	  general[GENERAL_STORING]);

 ErrorHandling: 
 
  PutRNGstate();
  
  if (err > NOERROR) {
    XERR(err);
  }
}


int check_density(cov_model *cov) {
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = PrevLoc(cov);
  int j, err, role, iso;
  Types type;
  //  bool vdim_close_together = GLOBAL.general.vdim_close_together;

  ASSERT_LOC_GIVEN;

 
  cov->simu.expected_number_simu = GLOBAL.general.expected_number_simu;
  if (cov->simu.expected_number_simu > 1 && 
      GLOBAL.general.expected_number_simu <= 1)
    SERR("expected number of simulations inconsistent");

  GLOBAL.internal.stored_init = GLOBAL.general.storing || 
    GLOBAL.general.expected_number_simu > 1;  
  
  if (cov->key == NULL) {
    domain_type dom = KERNEL;
    //char msg[200];
    if (isProcess(sub)) {
      role = ROLE_GAUSS;
      type = ProcessType;
      iso = UNREDUCED;
    } else {
      role = ROLE_COV;
      type = PosDefType;
      iso = SYMMETRIC;
    }
    if (cov->role == ROLE_BASE) role = ROLE_BASE;

    err = ERRORTYPECONSISTENCY;

    for (j=0; j<=2; j++) {
      if ((TypeConsistency(type, sub, 0) && 
	   (err = CHECK(sub, loc->timespacedim, cov->xdimown, type, 
			dom, iso, cov->vdim, role)) == NOERROR) 
	  || isProcess(sub)) break;

      if (j==0) type = VariogramType;
      else {
	type = TrendType;
	dom = XONLY;
	iso = CARTESIAN_COORD;
      }
    }

    if (err != NOERROR) return err;
  } else {
    BUG;
    role = role_of_process(sub->nr);
    if (role == ROLE_FAILED) BUG;

    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
		     XONLY, 
		     isCartesian(cov->isoprev) ? CARTESIAN_COORD : cov->isoprev,
		     cov->vdim, role)) != NOERROR) {
      return err;
    }
  }

  setbackward(cov, sub);  
  int subvdim = sub->vdim[0];
  cov->vdim[0]=subvdim; 
  cov->vdim[1]=sub->vdim[1]; 

  if (cov->q == NULL) {
    QALLOC(1);
    cov->q[0] = 1.0; 
  }

  return NOERROR;
}


int struct_density(cov_model *cov, cov_model VARIABLE_IS_NOT_USED  **newmodel){
  cov_model *next = cov->sub[0],
    *sub = next;
  location_type *loc = PrevLoc(cov);
  int err,
    subrole = ROLE_FAILED,
    nr = next->nr;

  //APMI(next);

  if (isVariogram(next) || isTrend(next)) {
    if ((err = covCpy(&(cov->key), next)) != NOERROR) return err;
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;
    
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
		     XONLY, 
		     isCartesian(cov->isoprev) ? CARTESIAN_COORD : cov->isoprev, 
		     cov->vdim, ROLE_GAUSS)) != NOERROR) {
      return err;
    }
    subrole = ROLE_GAUSS;    
  } else if (isBernoulliProcess(next)) subrole = ROLE_BERNOULLI;
  else if (isGaussBasedProcess(next)) subrole = ROLE_GAUSS;   
  else if (isBrownResnickProcess(next)) subrole = ROLE_BROWNRESNICK;    
  else if (nr == POISSONPROC) subrole = ROLE_POISSON;
  else if (nr == SCHLATHERPROC) subrole = ROLE_SCHLATHER;
  else if (nr == EXTREMALTPROC) subrole = ROLE_SCHLATHER;
  else if (nr == SMITHPROC) subrole = ROLE_SMITH;
  else {
    ILLEGAL_ROLE;
  }
  sub->role = subrole; // ansonsten muesste hier C-HECK aufgerufen werden
  // hoffentlich gut gehende Abkuerzung, dass S-TRUCT aufgerufen wird,
  // und danach C-HECK (was auf jeden Fall gemacht werden muss)

  cov->simu.active = next->simu.active = false; 
  sub->simu.expected_number_simu = cov->simu.expected_number_simu;

  if (PL >= PL_DETAILS) PRINTF("Struct Density\n");

  if ((err = STRUCT(sub, NULL)) != NOERROR) { return err; }
  if (PL >= PL_DETAILS) PRINTF("Checking Density\n");


  assert(cov->Sgen == NULL);
  NEW_STORAGE(gen);

  if (!sub->initialised) {
    if (PL >= PL_DETAILS) PRINTF("Struct Density C\n");
   
    if (//cov->key != NULL && // bei sub==next waere der falsche role gesetzt
	// irgenwie komisch, cov->key abzufragen und check(sub aufzurufen
	// aufgrund von Beispiel in RTpoisson geloescht
	(err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
		     cov->domprev, cov->isoprev, cov->vdim,
		     subrole)) != NOERROR) {
      //
      return err;
    }
    
    if (PL >= PL_DETAILS) {
      PRINTF("\n\nStruct Density (%s, #=%d), after 2nd check:",
	     NICK(sub), sub->gatternr); // OK
      PMI(sub); // OK
    }
 
    
    if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) {
      //APMI(cov); // !!! ?? hier weitermachen
      return err; 
    }
  }
 
  cov->initialised = true;
  cov->fieldreturn = true;
  cov->origrf = false;
  cov->rf = sub->rf;

  cov->simu.active = sub->simu.active = true; 
  return NOERROR;
}

void range_density(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[DENS_LOG] = 0;
  range->max[DENS_LOG] = 1;
  range->pmin[DENS_LOG] = 0;
  range->pmax[DENS_LOG] = 1;
  range->openmin[DENS_LOG] = false;
  range->openmax[DENS_LOG] = false;
}
*/


#define SIMU_CHECKONLY 0
#define SIMU_SEED 1
#define SIMU_ENV 2

void simulate(double *N, cov_model *cov, double *v){
  assert(!P0INT(SIMU_CHECKONLY));
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  char errorloc_save[nErrorLoc],
    format[20],
    back[]="\b\b\b\b\b\b\b\b\b\b\b", 
    prozent[]="%",
    pch = GLOBAL.general.pch;
  int ni, digits, //n, 
    nn,
    err = NOERROR,
    each = 0;
  double *res,
    realeach=0.0;
  simu_type *simu = NULL;  
  location_type *loc = PrevLoc(cov);
  long vdimtot; vdimtot = loc->totalpoints * cov->vdim[0];
  assert(cov->vdim[0] == cov->vdim[1]);

  // die folgenden Zeilen duefen nicht abgeaendert werden!!
  // Grund; *N wird eventuell durch Koordinatentrafo veraendert
  // so dass bei v=NULL ohne Umweg ueber gatter aufgerufen wird
  // und der korrekte Wert in cov->q gespeichert wird, der dann
  // spaeter wieder ausgelesen wird und nn zugeordnet wird.
  if (v == NULL) {
    cov->q[cov->qlen - 1] =  (int) *N;
    return; // EvaluateModel needs information about size
  //                      of result array
  } 
  nn = (int) cov->q[cov->qlen - 1];

 strcpy(errorloc_save, ERROR_LOC);
  simu = &(cov->simu);
  if (!simu->active) {
    err=ERRORNOTINITIALIZED; goto ErrorHandling;
  }

  if (nn>1 && pch != '\0') {
    if (pch == '!') {
      digits = (nn<900000000) 
	? 1 + (int) trunc(log((double) nn) / log(10.0)) : 9;
      back[digits] = '\0';
      each = (nn < 100) ? 1 :  nn / 100;
      sprintf(format, "%ss%s%dd", prozent, prozent, digits);
    } else if (pch == '%') {
      back[4] = '\0';
      realeach = (double) nn / 100.0;
      each = (nn < 100) ? 1 : (int) realeach;
      sprintf(format, "%ss%s%dd%ss", prozent, prozent, 3, prozent);
    } else each = 1;
  } else each = nn + 1;
  // size = cov->vdim * nn * loc->totalpoints;
  res = v;

  sub->simu.pair = false;
  for (ni=1; ni<=nn; ni++, res += vdimtot) { 
    
   if (GLOBAL_UTILS->basic.seed != NA_INTEGER && nn > 1) {
      if (PisNULL(SIMU_SEED) || PisNULL(SIMU_ENV)) {
	BUG;
      } else {
 	PutRNGstate();
	double simu_seed = GLOBAL_UTILS->basic.seed + (ni - 1);
	addVariable((char*) "seed", &simu_seed, 1, 1, PENV(SIMU_ENV)->sexp);
	eval(PLANG(SIMU_SEED)->sexp, PENV(SIMU_ENV)->sexp);
	GetRNGstate(); 
      }
    }

  sprintf(ERROR_LOC, "%s %d: ", errorloc_save, ni);
 
  R_CheckUserInterrupt();

    if (PL > 0) {
     if (ni % each == 0) {
	if (pch == '!')  
	PRINTF(format, back, ni / each);
	else if (pch == '%')
	  PRINTF(format, back, (int) (ni / realeach), prozent);
	else PRINTF("%c", pch);
      }
    }

    assert(cov->Sgen != NULL);
    DO(sub, cov->Sgen);

    if (sizeof(double) == sizeof(double) && false) {
      MEMCOPY(res, cov->rf, sizeof(double) * vdimtot);
    } else {
      int i; for (i=0; i<vdimtot; i++) {
	res[i] = cov->rf[i];
      }
    }
     
   if (!sub->simu.active) 
      GERR1("could not perform multiple simulations. Is '%s == FALSE'?",
	   general[GENERAL_STORING]);
     
  }

  if (nn>1 && pch != '\0') {
    if (pch == '!' || pch == '%') PRINTF("%s", back);
    else PRINTF("\n");
  }

  assert(simu != NULL);
  sub->simu.active = simu->active = sub->simu.active && GLOBAL.general.storing;

 ErrorHandling: 
 
  PutRNGstate();
  
  if (err > NOERROR) {
    if (simu != NULL) sub->simu.active = simu->active = false;
    XERR(err);
  }
 
}


int check_simulate(cov_model *cov) {
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = PrevLoc(cov);
  int j, d, role, iso,
    err = NOERROR;
  Types type;
  bool vdim_close_together = GLOBAL.general.vdim_close_together;
  
  ASSERT_LOC_GIVEN;
  kdefault(cov, SIMU_CHECKONLY, false);

  cov->simu.expected_number_simu = GLOBAL.general.expected_number_simu;
  if (cov->simu.expected_number_simu > 1 && 
      GLOBAL.general.expected_number_simu <= 1)
    SERR("expected number of simulations inconsistent");

  GLOBAL.internal.stored_init = GLOBAL.general.storing || 
    GLOBAL.general.expected_number_simu > 1;  
  
  if (cov->key == NULL) {
    domain_type dom; 
    char errmsg[LENERRMSG];
    //char msg[200];
    if (isProcess(sub)) {
      dom = XONLY; // 5/2015; voher beides Kerne
      role = ROLE_GAUSS;
      type = ProcessType;
      iso = UNREDUCED;
    } else {
      dom = KERNEL;
      role = ROLE_COV;
      type = PosDefType;
      iso = SymmetricOf(cov->isoprev);
    }
    assert(isCoordinateSystem(cov->isoprev));
    
    if (cov->role == ROLE_BASE) role = ROLE_BASE;

    err = ERRORTYPECONSISTENCY;
    errorMSG(err, errmsg);
 
    for (j=0; j<=2; j++) {

      if ( (TypeConsistency(type, sub, 0) && 
	   (err = CHECK(sub, loc->timespacedim, cov->xdimown, type, 
			dom, iso, cov->vdim, role)) == NOERROR) 
	   || isProcess(sub)) {// muss als zweites stehen !! 	
	break;
      }
 
      if (j==0) {
	type = VariogramType;
	errorMSG(err, errmsg);
      } else {
	type = TrendType;
	dom = XONLY;
	iso = cov->isoprev;
      }
    }
    if (err != NOERROR) SERR(errmsg);
  } else {
    role = role_of_process(sub->nr);
    if (role == ROLE_FAILED) BUG;

    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
		     XONLY, 
		     isCartesian(cov->isoprev) ? CARTESIAN_COORD : cov->isoprev,
		     cov->vdim, role)) != NOERROR) {
      return err;
    }
  }
  

   setbackward(cov, sub);  
  int subvdim = sub->vdim[0];
  cov->vdim[0]=subvdim; 
  cov->vdim[1]=sub->vdim[1]; 

  if (cov->q == NULL) {
    int len=1;
    if (loc->grid) len += loc->timespacedim; else len++;
    if (subvdim > 1) len++;
    QALLOC(len);
    cov->q[--len] = 1.0; // number of simulations
    if (subvdim > 1 && !vdim_close_together) cov->q[--len]=(double) subvdim;
    if (loc->grid) {
      for (d=loc->timespacedim-1; d>=0; d--) 
	cov->q[--len] = loc->xgr[d][XLENGTH];
    } else {
      cov->q[--len] = loc->totalpoints;
    }
    if (subvdim > 1 && vdim_close_together) {
      assert(len == 1);
      cov->q[--len]=(double) subvdim;
    }
    assert(len==0);
  }

  return NOERROR;
}


int struct_simulate(cov_model *cov, cov_model VARIABLE_IS_NOT_USED  **newmodel){
  cov_model *next = cov->sub[0],
    *sub = next;
  location_type *loc = PrevLoc(cov);
  int err,
    subrole = ROLE_FAILED,
    nr = next->nr;

  if (isVariogram(next) || isTrend(next)) {
    if ((err = covCpy(&(cov->key), next)) != NOERROR) return err;
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;

    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
		     XONLY, 
		     isCartesian(cov->isoprev) ? CARTESIAN_COORD :cov->isoprev, 
		     cov->vdim, ROLE_GAUSS)) != NOERROR) {
      return err;
    }
    subrole = ROLE_GAUSS;    
  } else if (isBernoulliProcess(next)) subrole = ROLE_BERNOULLI;
  else if (isGaussBasedProcess(next)) subrole = ROLE_GAUSS;   
  else if (isBrownResnickProcess(next)) subrole = ROLE_BROWNRESNICK;    
  else if (nr == POISSONPROC) subrole = ROLE_POISSON;
  else if (nr == SCHLATHERPROC) subrole = ROLE_SCHLATHER;
  else if (nr == EXTREMALTPROC) subrole = ROLE_SCHLATHER;
  else if (nr == SMITHPROC) subrole = ROLE_SMITH;
  // else if (isTrend(next)) subrole = ROLE_GAUSS;
  else {
    ILLEGAL_ROLE;
  }
  sub->role = subrole; // ansonsten muesste hier C-HECK aufgerufen werden
  // hoffentlich gut gehende Abkuerzung, dass S-TRUCT aufgerufen wird,
  // und danach C-HECK (was auf jeden Fall gemacht werden muss)

  cov->simu.active = next->simu.active = false; 
  sub->simu.expected_number_simu = cov->simu.expected_number_simu;
  if (P0INT(SIMU_CHECKONLY)) return NOERROR;

  if (PL >= PL_DETAILS) PRINTF("Struct Simulate\n");

  if ((err = STRUCT(sub, NULL)) != NOERROR) { return err; }
  if (PL >= PL_DETAILS) PRINTF("Checking Simulate\n");


  assert(cov->Sgen == NULL);
  NEW_STORAGE(gen);

  if (!sub->initialised) {
    if (PL >= PL_DETAILS) PRINTF("Struct Simulate C\n");
   
    if (//cov->key != NULL && // bei sub==next waere der falsche role gesetzt
	// irgenwie komisch, cov->key abzufragen und check(sub aufzurufen
	// aufgrund von Beispiel in RTpoisson geloescht
	(err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
		     cov->domprev, cov->isoprev, cov->vdim,
		     subrole)) != NOERROR) {
      //
      return err;
    }
    
    if (PL >= PL_DETAILS) {
      PRINTF("\n\nStruct Simulate (%s, #=%d), after 2nd check:",
	     NICK(sub), sub->gatternr); // OK
      PMI(sub); // OK
    }
 
    if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) {
      //APMI(cov); // !!! ?? hier weitermachen
      return err; 
    }
  }
 
  cov->initialised = true;
  cov->fieldreturn = true;
  cov->origrf = false;
  cov->rf = sub->rf;

  cov->simu.active = sub->simu.active = true; 

  return NOERROR;
}

void range_simulate(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[SIMU_CHECKONLY] = 0;
  range->max[SIMU_CHECKONLY] = 1;
  range->pmin[SIMU_CHECKONLY] = 0;
  range->pmax[SIMU_CHECKONLY] = 1;
  range->openmin[SIMU_CHECKONLY] = false;
  range->openmax[SIMU_CHECKONLY] = false;

#define simu_n 2
  int idx[simu_n] = {SIMU_SEED, SIMU_ENV};
  for (int i=0; i<simu_n; i++) {
    int j = idx[i];
    range->min[j] = RF_NAN;
    range->max[j] = RF_NAN;
    range->pmin[j] = RF_NAN;
    range->pmax[j] = RF_NAN;
    range->openmin[j] = false;
    range->openmax[j] = false; 
  }


}

//void do_simulate(cov_model *cov, gen_storage *s){
//  assert(false);
//}

/////////////////////////////////////////




double GetPriors(cov_model *cov) {
  cov_fct *C = CovList + cov->nr;
  cov_model *kap;
  int i,
    kappas = C->kappas,
    nsub = cov->nsub;
  double v,
    logli = 0.0;
  for (i=0; i<kappas; i++) {
    if ((kap = cov->kappasub[i]) != NULL) {
      if (isRandom(kap->typus)) {
	if (C->kappatype[i] < LISTOF) {
	  assert(C->kappatype[i] == REALSXP);
	  VTLG_DLOG(P(i), kap, &v);
	} else if (C->kappatype[i] < LISTOF + LISTOF) { // not tested yet
	  NotProgrammedYet("hierachical models for multiple data sets");
	  assert(C->kappatype[i] == LISTOF + REALSXP);
	  int
	    store = GLOBAL.general.set,
	    end = cov->nrow[i];
	  for (GLOBAL.general.set = 0; GLOBAL.general.set<end; 
	       GLOBAL.general.set++) {
	    VTLG_DLOG(LP(i), kap, &v);
	    logli += v;
	  }
	  GLOBAL.general.set = store;
	} else BUG;
	logli += v;
      }
      logli += GetPriors(kap);
    }
  }
  
  for (i=0; i<nsub; i++) {
    logli += GetPriors(cov->sub[i]);
  }

  return logli; 
}





/* *************** */
/* LOG LIKELIHOOD */
/* *************** */

 void kappalikelihood(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc) {
  *nr = *nc = i == LIKELIHOOD_DATA ? 0 
    : i <= LIKELIHOOD_LAST ? 1 : -1;
}

void likelihood(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v) { 
  cov_model *process = cov->key == NULL ? cov->sub[0] : cov->key;
  if (v == NULL) {
    likelihood_storage *L = process->Slikelihood;
    assert(L != NULL);
    int //betas = L->betas[L->fixedtrends],
      vdim = process->vdim[0];
    likelihood_info *info = &(L->info);
    listoftype *datasets = L->datasets;
    int
      store = GLOBAL.general.set,
      betatot = L->betas[L->fixedtrends],
      all_betatot = betatot;
    GLOBAL.general.set = 0;
    if (L->betas_separate)
      all_betatot *= NCOL_OUT_OF(datasets) / vdim;
    cov->q[0] = 1 + info->globalvariance + all_betatot;
     GLOBAL.general.set = store;
    return; 
  }

  assert(isProcess(process));
  VTLG_DLOG(NULL, process, v); 

  assert(process->key == NULL);
  *v += GetPriors(process->sub[0]); 
 }

 

int check_likelihood(cov_model *cov) {
  int err,
    sets = GET_LOC_SETS(cov);
  int store = GLOBAL.general.set; 

  if ((err = check_linearpart(cov))) return err;
 
  kdefault(cov, LIKELIHOOD_NA_VAR, GLOBAL.fit.estimate_variance);
  kdefault(cov, LIKELIHOOD_BETASSEPARATE, false);
  if (P0INT(LIKELIHOOD_BETASSEPARATE)) BUG; // to estimate beta 
  // for each independent repetition within a data set is a feature
  // that has been implemented, but currently it is not used
  kdefault(cov, LIKELIHOOD_IGNORETREND, false);
  if (PisNULL(LIKELIHOOD_DATA)) BUG;
  for (GLOBAL.general.set = 0; GLOBAL.general.set<sets; GLOBAL.general.set++) {
    long 
      datatot = LNROW(LIKELIHOOD_DATA) * LNCOL(LIKELIHOOD_DATA),
      totpts = Gettotalpoints(cov),
      vdimtot = totpts * cov->vdim[0];
    int
      repet = datatot / vdimtot;

    //printf("%d %d %d %d\n", vdimtot, totpts, repet, datatot);

    if (repet * vdimtot != datatot || repet == 0)  {
      GERR("data and coordinates do not match");
    }
    LNRC_(LIKELIHOOD_DATA, nrow) = totpts;
    LNRC_(LIKELIHOOD_DATA, ncol) = datatot / totpts;
  }

 ErrorHandling:
  GLOBAL.general.set = store;

  return err;
}

int struct_likelihood(cov_model *cov, 
		      cov_model VARIABLE_IS_NOT_USED  **newmodel){
  //return struct_rftrend(cov, newmodel);
  assert(cov->key == NULL);
  cov_model *sub = cov->sub[0];
  int err = NOERROR;
  location_type *loc = Loc(cov);

  if (isVariogram(sub)) {
    if ((err = covCpy(&(cov->key), sub)) != NOERROR) return err;
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;
    
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType, XONLY, 
		     isCartesian(cov->isoprev) ? CARTESIAN_COORD : cov->isoprev,
		     cov->vdim, ROLE_GAUSS)) != NOERROR) {
      return err;
    }
  }
  if (!isProcess(sub)) 
    SERR1("'%s' can be calculated only for processes.", NICK(cov));
  sub->role = ROLE_LIKELIHOOD;

  if ((err = STRUCT(sub, NULL)) != NOERROR) return err;
  return err;
}


void range_likelihood(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[LIKELIHOOD_DATA] = RF_NEGINF;
  range->max[LIKELIHOOD_DATA] = RF_INF;
  range->pmin[LIKELIHOOD_DATA] = -1e8;
  range->pmax[LIKELIHOOD_DATA] = 1e8;
  range->openmin[LIKELIHOOD_DATA] = true;
  range->openmax[LIKELIHOOD_DATA] = true;

  range->min[LIKELIHOOD_NA_VAR] = 0;
  range->max[LIKELIHOOD_NA_VAR] = 1;
  range->pmin[LIKELIHOOD_NA_VAR] = 0;
  range->pmax[LIKELIHOOD_NA_VAR] = 1;
  range->openmin[LIKELIHOOD_NA_VAR] = false;
  range->openmax[LIKELIHOOD_NA_VAR] = false;

  range->min[LIKELIHOOD_BETASSEPARATE] = 0;
  range->max[LIKELIHOOD_BETASSEPARATE] = 1;
  range->pmin[LIKELIHOOD_BETASSEPARATE] = 0;
  range->pmax[LIKELIHOOD_BETASSEPARATE] = 1;
  range->openmin[LIKELIHOOD_BETASSEPARATE] = false;
  range->openmax[LIKELIHOOD_BETASSEPARATE] = false;

  range->min[LIKELIHOOD_IGNORETREND] = 0;
  range->max[LIKELIHOOD_IGNORETREND] = 1;
  range->pmin[LIKELIHOOD_IGNORETREND] = 0;
  range->pmax[LIKELIHOOD_IGNORETREND] = 1;
  range->openmin[LIKELIHOOD_IGNORETREND] = false;
  range->openmax[LIKELIHOOD_IGNORETREND] = false;
}



void linearpart(double VARIABLE_IS_NOT_USED *x, cov_model  VARIABLE_IS_NOT_USED *cov, double  VARIABLE_IS_NOT_USED *v) { 
  BUG; // darf nicht aufgerufen werden
}

SEXP get_linearpart(SEXP model_reg, SEXP Set){
  currentRegister = INTEGER(model_reg)[0];			
  if (currentRegister < 0 || currentRegister > MODEL_MAX) BUG;	
  cov_model *cov = KEY[currentRegister];			
  cov = cov->key != NULL ? cov->key : cov->sub[0];		
  if (cov->nr == GAUSSPROC) 
    return gauss_linearpart(model_reg, Set);
  else BUG;
}


int check_linearpart(cov_model *cov) {
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  int j, err, role, iso;
  Types type;
 
  assert(cov->prevloc != NULL);
  
  domain_type dom = KERNEL;
  //char msg[200];
  if (isProcess(sub)) {
    role = ROLE_GAUSS; 
    type = ProcessType;
    iso = UNREDUCED;
  } else {
    role = ROLE_COV;
    type = PosDefType;
    iso = SymmetricOf(cov->isoprev);
  }
  if (cov->role == ROLE_BASE) role = ROLE_BASE;
  
  err = ERRORTYPECONSISTENCY;
  
  for (j=0; j<=1; j++) {
    if ((TypeConsistency(type, sub, 0) && 
	 (err = CHECK(sub,  Gettimespacedim(cov), cov->xdimown, type, 
		      dom, iso, cov->vdim, role)) == NOERROR)
	|| isProcess(sub)) break;
    
    if (j==0) type = VariogramType;
  }
  if (err != NOERROR) goto ErrorHandling;

  setbackward(cov, sub);  
  cov->vdim[0]=sub->vdim[0]; 
  cov->vdim[1]=sub->vdim[1]; 

  if (cov->q == NULL) QALLOC(2);
  cov->q[0] = Gettotalpoints(cov);
  cov->q[1] = cov->vdim[0];
  
 ErrorHandling:
  return err;
}



int struct_linearpart(cov_model *cov, 
    cov_model VARIABLE_IS_NOT_USED  **newmodel){
  assert(cov->key == NULL);
  cov_model *sub = cov->sub[0];
  int err = NOERROR;
  location_type *loc = Loc(cov);

  if (isVariogram(sub)) {
    if ((err = covCpy(&(cov->key), sub)) != NOERROR) return err;
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;
    
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType, XONLY, 
		     isCartesian(cov->isoprev) ? CARTESIAN_COORD : cov->isoprev,
		     cov->vdim, ROLE_GAUSS)) != NOERROR) {
      return err;
    }
  }
  if (!isProcess(sub)) 
    SERR1("'%s' can be calculated only for processes.", NICK(cov));

  sub->role = ROLE_LIKELIHOOD;

  if ((err = STRUCT(sub, NULL)) != NOERROR) return err;
  likelihood_storage *L = sub->Slikelihood;
  if (L == NULL) return ERRORFAILED;
  if (L->dettrend_has_nas || L->fixedtrend_has_nas) {
    warning("NAs detected in '%s'; hence zero's introduced", NICK(cov));
  }
  return err;
}




#define EVALDISTR_X 0 // d
#define EVALDISTR_Q 1 // p
#define EVALDISTR_P 2 // q
#define EVALDISTR_N 3 // r
#define EVALDISTR_DIM 4 // r

void kappa_EvalDistr(int i, cov_model VARIABLE_IS_NOT_USED *cov, 
		     int *nr, int *nc){
  *nc = *nr = i <= EVALDISTR_N ? 0 : i == EVALDISTR_DIM ? 1 : -1;
}

void EvalDistr(double VARIABLE_IS_NOT_USED *N, cov_model *cov, double *v){
  char errorloc_save[nErrorLoc];
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  double  *xqp;
  int i, j,
    dim = cov->tsdim,
    n = (int) (cov->q[cov->qlen - 1]);

  if (v==NULL) return; // EvaluateModel needs information about size
  strcpy(errorloc_save, ERROR_LOC);

  if (!PisNULL(EVALDISTR_X)) { // d
    xqp = P(EVALDISTR_X);
    for (j=i=0; i<n; i++, j+=dim) VTLG_D(xqp + j, sub, v + i);
  } else if (!PisNULL(EVALDISTR_Q)) { // p
    xqp = P(EVALDISTR_Q);
    for (j=i=0; i<n; i++, j+=dim) VTLG_P(xqp + i, sub, v + j);
  } else if (!PisNULL(EVALDISTR_P)) {// q
    xqp = P(EVALDISTR_P);
    for (j=i=0; i<n; i++, j+=dim) VTLG_Q(xqp + j, sub, v + i);
  } else if (!PisNULL(EVALDISTR_N)) {// r
    xqp = P(EVALDISTR_N);
    for (j=i=0; i<n; i++, j+=dim) {
      VTLG_R(NULL, sub, v + j);    
    }
  } else BUG;
}


int check_EvalDistr(cov_model *cov) {
  cov_fct *C = CovList + cov->nr;
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  int size_idx,  err, //type,  nr = sub->nr,
    *nrow = cov->nrow,
    *ncol = cov->ncol,
    dim = cov->tsdim, 
    zaehler=0;

  ROLE_ASSERT(ROLE_DISTR);
 
  if (cov->q == NULL) {
    int nn = 1; // !! 1 obwohl 2 reserviert werden !! Wg EvaluateModel
    if (dim > 1 && 
	((!PisNULL(EVALDISTR_N) && P0(EVALDISTR_N) > 1) ||
	 (!PisNULL(EVALDISTR_Q) && P0(EVALDISTR_Q) > 1)))
      nn++;
    QALLOC(nn + 1); 
    cov->qlen = nn;
   /*  cov->qlen = 1; // !! 1 obwohl 2 reserviert werden !! Wg EvaluateModel
    if (dim > 1 && 
	((!PisNULL(EVALDISTR_N) && P0(EVALDISTR_N) > 1) ||
	 (!PisNULL(EVALDISTR_Q) && P0(EVALDISTR_Q) > 1)))
      cov->qlen++;
      cov->q = (double *) MALLOC(sizeof(double) * (cov->qlen + 1)); // QALLOC NOT APPROPRIATE
   */

    cov->q[0] = dim; // is overwritten if not a matrix is returned
    size_idx = cov->qlen - 1;

    if (PisNULL(EVALDISTR_N)) {
      if (!PisNULL(EVALDISTR_X)) { // d
	if (dim > 1 && nrow[EVALDISTR_X] != dim) 
	  SERR2("dimenson of '%s' does not match '%s' ", 
		C->kappanames[EVALDISTR_X], C->kappanames[EVALDISTR_DIM]);
	cov->q[size_idx] = nrow[EVALDISTR_X] * ncol[EVALDISTR_X] / dim;
	zaehler++;
      }
      if (!PisNULL(EVALDISTR_Q)) { // p
	if (dim > 1 && nrow[EVALDISTR_Q] != dim) 
	  SERR2("dimension of '%s' does not match '%s' ", 
		C->kappanames[EVALDISTR_Q], C->kappanames[EVALDISTR_DIM]);
	cov->q[size_idx] = nrow[EVALDISTR_Q] * ncol[EVALDISTR_Q] / dim;
	zaehler++;
      } 
      if (!PisNULL(EVALDISTR_P)) { // q
	if (ncol[EVALDISTR_P] != 1) 
	  SERR1("'%s' must be a vector", C->kappanames[EVALDISTR_P]);
	cov->q[size_idx] = nrow[EVALDISTR_P]  * dim;
 	zaehler++;
      } 
    }// kein else wegen zaehler !!

    if (!PisNULL(EVALDISTR_N)) { // r      
      cov->q[size_idx] = P0(EVALDISTR_N) * dim;
      zaehler++;
    } 
    if (zaehler != 1) SERR("exactly one of the parameters must be given");

  }
 
  if (!isRandom(sub)) {
    SERR1("'%s' is not a distribution", NICK(sub));
  }

  // PMI(sub);
  
  if ((err = CHECK_R(sub, dim)) != NOERROR) return err;
 
  setbackward(cov, sub);  

  return NOERROR;
}

// S TRUCT(!isCovariance(cov) ? NULL : &(cov->key));



void range_EvalDistr(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[EVALDISTR_X] = RF_NEGINF;
  range->max[EVALDISTR_X] = RF_INF;
  range->pmin[EVALDISTR_X] = -1e8;
  range->pmax[EVALDISTR_X] = 1e8;
  range->openmin[EVALDISTR_X] = true;
  range->openmax[EVALDISTR_X] = true;

  range->min[EVALDISTR_Q] = RF_NEGINF;
  range->max[EVALDISTR_Q] = RF_INF;
  range->pmin[EVALDISTR_Q] = -1e8;
  range->pmax[EVALDISTR_Q] = 1e8;
  range->openmin[EVALDISTR_Q] = true;
  range->openmax[EVALDISTR_Q] = true;

  range->min[EVALDISTR_P] = 0;
  range->max[EVALDISTR_P] = 1;
  range->pmin[EVALDISTR_P] = 0;
  range->pmax[EVALDISTR_P] = 1;
  range->openmin[EVALDISTR_P] = false;
  range->openmax[EVALDISTR_P] = false;

  range->min[EVALDISTR_N] = 1;
  range->max[EVALDISTR_N] = RF_INF;
  range->pmin[EVALDISTR_N] = 1;
  range->pmax[EVALDISTR_N] = 1e8;
  range->openmin[EVALDISTR_N] = false;
  range->openmax[EVALDISTR_N] = false;

  range->min[EVALDISTR_DIM] = 1;
  range->max[EVALDISTR_DIM] = RF_INF;
  range->pmin[EVALDISTR_DIM] = 1;
  range->pmax[EVALDISTR_DIM] = 10;
  range->openmin[EVALDISTR_DIM] = false;
  range->openmax[EVALDISTR_DIM] = true;
}


int struct_EvalDistr(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){
  cov_model *next = cov->sub[0],
    *sub = next;
  int  err,
    dim = cov->tsdim;
   
  cov->simu.active = next->simu.active = false; 
  if (PL >= PL_DETAILS) PRINTF("Struct EvalDistr\n");

  if ((err = STRUCT(sub, NULL)) != NOERROR) { return err; }
  if (PL >= PL_DETAILS) PRINTF("Checking EvalDistr\n");
  
  if ((err = CHECK_R(sub, dim)) != NOERROR) return err;
    
  if (PL >= PL_DETAILS) {
    PRINTF("\n\nStruct EvalDistr (%s, #=%d), after 2nd check:",
	   NICK(sub), sub->gatternr); // OK
  }
  
  assert(cov->Sgen == NULL);
  NEW_STORAGE(gen);

  if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) {
    //APMI(cov); // !!! ?? hier weitermachen
    return err; 
  }

  if (cov->rf == NULL) {
    int size = cov->q[0];
    if (cov->qlen > 1) size *= cov->q[1];

    if ((cov->rf = (double*) MALLOC(sizeof(double) * size)) 
	== NULL) return ERRORMEMORYALLOCATION; 
    cov->fieldreturn = cov->origrf = true;
  }
  cov->simu.active = sub->simu.active = true; 
  return NOERROR;
}


//void do_EvalDistr(cov_model *cov, gen_storage *s){
//  assert(false);
//}


void Dummy(double VARIABLE_IS_NOT_USED *x, cov_model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *value) {
  BUG; //ERR("unexpected call of dummy function");
}


int check_dummy(cov_model *cov) {
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = PrevLoc(cov);
  int i, 
    err = NOERROR;

  ASSERT_LOC_GIVEN;

  for (i=0; i<=1; i++) {
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, VariogramType,
		     i==0 ? XONLY : KERNEL, 
		     SymmetricOf(cov->isoprev), SUBMODEL_DEP,
		     ROLE_COV)) == NOERROR) break;
  }
  if (err != NOERROR) return err;
  setbackward(cov, sub);  
  cov->vdim[0] = sub->vdim[0]; 
  cov->vdim[1] = sub->vdim[1]; 

  return NOERROR;
}


int struct_dummy(cov_model VARIABLE_IS_NOT_USED *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){    
  return NOERROR;
}


int alloc_pgs(cov_model *cov) {
  if (cov->xdimown != cov->tsdim) BUG;
  return alloc_pgs(cov, cov->xdimown);
}

int alloc_pgs(cov_model *cov, int dim) { // all what is necessary for dompp
  pgs_storage *pgs = NULL;

  assert(cov->Spgs == NULL); // NIE pgs_DELETE(&(cov->Spgs)) in Huetchen, da sonst dompp durcheinander kommt;
  NEW_STORAGE(pgs);
  pgs = cov->Spgs;

  if ((pgs->supportmin = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->supportmax = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->supportcentre = (double*) CALLOC(dim, sizeof(double))) == NULL || 
      (pgs->own_grid_start = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->own_grid_step = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->own_grid_len = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      
      (pgs->gridlen = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->end = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->start = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->delta = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->nx = (int*) CALLOC(dim, sizeof(int))) == NULL ||
   
      (pgs->xstart = (double*) CALLOC(dim, sizeof(double))) == NULL || 
      (pgs->x = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->inc = (double*) CALLOC(dim, sizeof(double))) == NULL)
    return ERRORMEMORYALLOCATION;

  return NOERROR;
}


int alloc_cov(cov_model *cov, int dim, int rows, int cols) { 
  // all what is necessary for dompp
  // but also for cov evaluation and fctn evaluation!

  int err;

   if (cov->Spgs != NULL) pgs_DELETE(&(cov->Spgs));
  if ((err = alloc_pgs(cov, dim)) != NOERROR) return err;

  // erst danach!!!
  pgs_storage *pgs = cov->Spgs;
  int max,
    rowscols = rows * cols;
  max = rows;
  if (cols > max) max = cols;

  assert(pgs != NULL && pgs->x!=NULL && pgs->endy==NULL && pgs->ptrrow==NULL);

  if ( 
      (pgs->endy = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->startny = (int*) CALLOC(dim, sizeof(int))) == NULL ||

      (pgs->ptrcol = (int*) CALLOC(max, sizeof(int))) == NULL ||
      (pgs->ptrrow = (int*) CALLOC(max, sizeof(int))) == NULL)
    return ERRORMEMORYALLOCATION;

  if (
      (pgs->C0x  = (double*) CALLOC(rowscols, sizeof(double))) == NULL || 
      (pgs->C0y  = (double*) CALLOC(rowscols, sizeof(double))) == NULL || 
      (pgs->cross= (double*) CALLOC(rowscols, sizeof(double))) == NULL || 
      (pgs->z    = (double*) CALLOC(rowscols, sizeof(double))) == NULL ||       
      (pgs->Val= (double**) CALLOC(rowscols, sizeof(double*))) == NULL 
      ) return ERRORMEMORYALLOCATION;

  pgs->rowscols = rowscols;

  return NOERROR;
}


#define RFGET_UP 0
#define RFGET_REGISTER 1
#define RFGET_SUB 0

void RFget(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v){
  get_storage *s = cov->Sget;
  cov_model *get_cov = s->get_cov;
  int i,
    nr = s->get_cov->nr,
    param_nr = s->param_nr,
    *idx = s->idx,
    size = s->size;

  if (CovList[nr].kappatype[param_nr] == REALSXP) {
    double *p = PARAM(get_cov, param_nr);    
    if (s->all) {
      for (i=0; i<size; i++) v[i] = p[i];
    } else {
      for (i=0; i<size; i++) v[i] = p[(int) idx[i]];
    }
  } else if (CovList[nr].kappatype[param_nr] == INTSXP) {    
    int *p = PARAMINT(get_cov, param_nr);
    if (s->all) {
      for (i=0; i<size; i++) v[i] = (double) p[i];
    } else {
      for (i=0; i<size; i++) v[i] = (double) p[idx[i]];
    }
  } else BUG;
}



int SearchParam(cov_model *cov, get_storage *s) {
  cov_model *orig;
  int i, subcov,
    up = P0INT(RFGET_UP);
  if (PisNULL(RFGET_REGISTER)) {
    orig = cov->calling;
    if (orig == NULL) SERR("register not given");
    for (i=0; i<up; i++) {
      orig = orig->calling;
      while (orig != NULL && orig->user_given == ug_internal)
	orig = orig->calling;
      if (orig == NULL) SERR1("value of '%s' too high", KNAME(RFGET_UP));
    }
  } else {
    int nr = P0INT(RFGET_REGISTER);
    if (nr<0 || nr>MODEL_MAX) SERR("invalid register number");
    if (up != 0) SERR1("'%s' may not be given.", KNAME(RFGET_UP));
    orig = KEY[nr];
  }
  s->orig = orig;
  
  while (true) {
    while (CovList[orig->nr].maxsub > 0 && orig != NULL &&
	     orig->user_given == ug_internal) 
      orig = (orig->nr == PLUS || orig->nr == MULT || orig->nr==MPPPLUS) 
	&& orig->Splus != NULL ? orig->Splus->keys[0]
	: orig->key != NULL ? orig->key
	: orig->sub[0];
    if (orig == NULL || CovList[orig->nr].maxsub == 0) 
      SERR("model does not match the request or is too complicated");
    if (cov->nr != orig->nr) 
      SERR2("(sub)models (%s, %s) do not match",
	    NICK(cov), NICK(orig));
    for (subcov=0; subcov < orig->nsub; subcov++) 
      if (cov->sub[subcov] != NULL) break;   
    if (subcov < orig->nsub) { // submodel found
      if (orig->sub[subcov] == NULL) 
	SERR2("(sub)models (%s; %s) do not match",
	      NICK(cov), NICK(orig));
      cov  = cov->sub[subcov];
      orig = orig->sub[subcov];
    } else {
      int
	kappas = CovList[orig->nr].kappas;
      for (i=0; i < kappas; i++) 
	if (cov->kappasub[i] != NULL) break;         
      if (i < kappas) { // param submodel found
	if (orig->kappasub[i] == NULL) 
	  SERR2("parameter models of %s and %s do not match",
		NICK(cov), NICK(orig));
        cov  = cov->kappasub[i];
	orig = orig->kappasub[i];   
      } else {
	for (i=0; i < kappas; i++) 
	  if (cov->kappasub[i] != NULL) break;         
	if (i >= kappas) SERR("no parameter given");
	cov_fct *C = CovList + cov->nr;
	if (C->kappatype[i] == REALSXP) s->all = P(i)[0] == 0;
	else if (C->kappatype[i] == INTSXP) s->all = PINT(i)[0] == 0;
	else SERR("only numeric parameters are allowed");
	if (s->all) {
	  s->vdim[0] = orig->nrow[i];
	  s->vdim[1] = orig->ncol[i];
	  s->size = s->vdim[0] * s->vdim[1];
	} else {
	  int k;
	  s->size = s->vdim[0] = cov->nrow[i];
	  s->vdim[1] = cov->ncol[i];
	  if (cov->ncol[i] != 1) SERR("only vectors of indices are allowed");
	  assert(s->idx == NULL);
	  s->idx = (int *) MALLOC(sizeof(int) * s->size);
	  if (C->kappatype[i] == REALSXP)
	    for (k=0; k<s->size; k++) s->idx[k] = ((int) P(i)[k]) - 1;
	  else 
	    for (k=0; k<s->size; k++) s->idx[k] = PINT(i)[k] - 1;
	}
	s->get_cov = orig;
	s->param_nr = i;
	return NOERROR;
      }
    }
  } // while true
  BUG;
  return ERRORFAILED; // nur fuer compiler
}
  
int check_RFget(cov_model *cov) {

  BUG; /// todo:  Code ist noch nicht ausgegoren !!

  //cov_model *orig, *sub;
  int i, err;
    //    len = ((idx != NULL) ? cov->nrow[RFGET_IDX]
    //	  : get->get_cov->ncol[param_nr] * get->get_cov->nrow[param_nr]);
  if (cov->Sget != NULL) return NOERROR;

  kdefault(cov, RFGET_UP, 0);
  NEW_STORAGE(get);
  get_storage *s = cov->Sget;

  if ((err = SearchParam(cov, s)) != NOERROR) return err;
  for (i=0; i<2; i++) cov->vdim[i] = s->vdim[i];
  return NOERROR;
}
  

void range_RFget(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[RFGET_UP] = 0;
  range->max[RFGET_UP] = RF_INF;
  range->pmin[RFGET_UP] = 0;
  range->pmax[RFGET_UP] = 100;
  range->openmin[RFGET_UP] = false;
  range->openmax[RFGET_UP] = true;

  range->min[RFGET_REGISTER] = 0;
  range->max[RFGET_REGISTER] = MODEL_MAX;
  range->pmin[RFGET_REGISTER] = 0;
  range->pmax[RFGET_REGISTER] = MODEL_USER;
  range->openmin[RFGET_REGISTER] = false;
  range->openmax[RFGET_REGISTER] = false;
}


int struct_RFget(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){
  int i, err;
  
  NEW_STORAGE(get);
  get_storage *s = cov->Sget;
  if ((err = SearchParam(cov, s)) != NOERROR) return err;
  for (i=0; i<2; i++) {
    if (cov->vdim[i] != s->vdim[i])
      SERR("mismatch of dimensions when constructing the model");
  }
  cov->fieldreturn = cov->origrf = false;  
  return NOERROR;
}




//void do_Rfget(cov_model *cov, gen_storage *s){
//  assert(false);
//}






//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

cov_model *get_around_gauss(cov_model *cov) {
  cov_model *nxt= cov;
  if (nxt->nr == SCHLATHERPROC) nxt = nxt->sub[0]; // kein else 
  if (nxt->nr == GAUSSPROC) nxt = nxt->sub[0];
  if (isGaussMethod(nxt) || isBernoulliProcess(nxt)) {
    if (nxt->nr == AVERAGE_USER){
      if (nxt->sub[0] == NULL) ERR("covariance cannot be calculated (yet) for arbitrary shape functions.");
      nxt = nxt->sub[nxt->sub[0] == NULL];
      if (nxt->nr == AVERAGE_INTERN) nxt = nxt->sub[nxt->sub[0] == NULL];
    } 
    else if (nxt->nr == CE_CUTOFFPROC_USER) {
      nxt = nxt->sub[0];      
      if (nxt->nr == CE_CUTOFFPROC_INTERN) nxt = nxt->sub[0];
    }
    else if (nxt->nr == CE_INTRINPROC_USER) {
      nxt = nxt->sub[0];   
      if (nxt->nr == CE_INTRINPROC_INTERN) nxt = nxt->sub[0];   
    }
    else if (nxt->nr == CE_INTRINPROC_USER) {
      nxt = nxt->sub[0];   
      if (nxt->nr == CE_INTRINPROC_INTERN) nxt = nxt->sub[0];   
    }
    else if (nxt->nr == HYPERPLANE_USER) {
      nxt = nxt->sub[0];   
      if (nxt->nr == HYPERPLANE_INTERN) nxt = nxt->sub[0];   
    }
    else if (nxt->nr == RANDOMCOIN_USER) {
      if (nxt->sub[0] == NULL) ERR("covariance cannot be calculated (yet) for arbitrary shape functions.");
      nxt = nxt->sub[nxt->sub[0] == NULL];   
      if (nxt->nr == AVERAGE_INTERN) nxt = nxt->sub[nxt->sub[0] == NULL];   
    } else {

    }
  }
  return nxt;
}



cov_model *get_around_max_stable(cov_model *cov) {
  cov_model *nxt = cov;

  if (isBrownResnickProcess(nxt)) {
    nxt = nxt->sub[0];
    if (nxt->calling->nr == BROWNRESNICKPROC && 
        isBRuserProcess(nxt)) {
      nxt = nxt->sub[0];
    }
  } 
  return nxt;
}



void Cov(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  CovList[sub->nr].covariance(sub, value);
}

int check_fct_intern(cov_model *cov, Types type, bool close, bool kernel,
		     int rows, int cols) {
  cov_model
    *next = cov->sub[0],
    *sub = cov->key == NULL ? next : cov->key;
  location_type *loc = Loc(cov);
  ASSERT_LOC_GIVEN;

  int err,  i, k, iso[4], // +  (int) nr_coord_sys_proj], 
    endfor = 0,
    lastdomain = kernel ? KERNEL : XONLY,
    dim = Gettimespacedim(cov);

  assert(dim == cov->tsdim);

  iso[endfor++] = type == ShapeType 
    ? CoordinateSystemOf(cov->isoown) 
    : SymmetricOf(cov->isoown);
  if (iso[endfor-1] == ISO_MISMATCH) BUG;

  for (i=0; i < endfor; i++) {
    for (k=XONLY; k<=lastdomain; k++) {
      if ((err = CHECK(sub, dim, cov->xdimown, type, k, iso[i], SUBMODEL_DEP,
		       sub!=next || isVariogram(sub) ? ROLE_COV : ROLE_BASE))
	  == NOERROR) break;
    }
    if (err == NOERROR) break;
  }
  if (err != NOERROR) return err;
  setbackward(cov, sub); 

  // memory set according to the submodel as this model will
  // be evaluated (cov clear; fctn as function; predict: cov-matrix!)
  if ((err = alloc_cov(cov, dim, cov->vdim[0], cov->vdim[1])) != NOERROR)
    return err;

  // this is how cov will forward the result
  if (rows > 0) cov->vdim[0] = rows;
  if (cols > 0) cov->vdim[1] = cols;

  if (sub->pref[Nothing] == PREF_NONE) SERR("given model cannot be evaluated")
  
  if (cov->q == NULL) {
    int d,
      len=1; // # of "simulations" (repetitions!)
    if (loc->grid) len += dim; else len ++;      
    for (i=0; i<2; i++) len += (int) (cov->vdim[i] > 1);
    QALLOC(len);

#define VDIMS  for (i=0; i<2; i++) if (cov->vdim[i] > 1) cov->q[d++] = cov->vdim[i]
#define LOCS if (loc->grid) {						\
      for (i=0; i<dim; i++) cov->q[d++] = loc->xgr[i][XLENGTH];		\
    } else {								\
      cov->q[d++] = loc->totalpoints;					\
    }      
    
    d = 0;
    if (close) {
      VDIMS;
      LOCS;	
   } else {
      LOCS;
      VDIMS;
    }
    cov->q[d] = 1;
    assert(d == len-1);
  }

  return NOERROR;
}

int check_cov_intern(cov_model *cov, Types type, bool close, bool kernel) {
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  if (isProcess(CovList[sub->nr].Typi[0])) {
    int err = CHECK(sub, cov->tsdim, cov->xdimown, ProcessType, XONLY, 
		    cov->isoown, SUBMODEL_DEP, 
		    cov->role == ROLE_BASE ? ROLE_BASE : ROLE_GAUSS);
    if (err != NOERROR) return err;
    setbackward(cov, sub);
    cov->vdim[0] = sub->vdim[0];
    cov->vdim[1] = sub->vdim[1];
    return NOERROR;
  } else return check_fct_intern(cov, type, close, kernel, 0, 0);
}

int check_cov(cov_model *cov) {  
  return check_cov_intern(cov, PosDefType, GLOBAL.general.vdim_close_together,
			  true);
}

int struct_cov(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){
  int err;
  cov_model *sub, 
    *next = cov->sub[0];
  //location_type *loc = PrevLoc(cov);
 
  sub = get_around_gauss(next);
  // fehlt : Poisson


  //printf("%s\n", NAME(cov)); assert(false);

  if (sub != next) {
    cov->key = sub; // illegaler bypass legen -- (gauss.)process darf nie 
    //                 aktiv aufgerufen werden
    err = cov->nr == COVMATRIX ? check_covmatrix(cov) : check_cov(cov);
    assert(cov->Spgs != NULL);
    cov->key = NULL; // bypass entfernen
    return err;

    /*
    if ((err = covCpy(&(cov->key), sub)) != NOERROR) return err;   
    sub = cov->key;
    sub->calling = cov;
    if (!isPosDef(sub->typus)) SERR("covariance model cannot be determined");
    
    //  printf("%ld %ld %ld %ld\n", next, sub, cov->sub[0], cov->key); 
    //  PMI(next);
    //   APMI(sub);
    return 
    */

    /*
      if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, 
		     PosDefType, 
		     loc->y == NULL && loc->ygr[0] == NULL ? XONLY : KERNEL, 
		     SymmetricOf(cov->isoprev), cov->vdim,
		     ROLE_COV)) != NOERROR) {
	return err;
      }
    */
  }
  return NOERROR;
}


void CovMatrix(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array

  //  APMI(cov);

  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  //printf("CovMatrix %s=%d %s=%d \n\n\n\n\n\n", NAME(cov), cov->Spgs != NULL, NAME(cov->sub[0]), cov->sub[0]->Spgs != NULL); 

   CovList[sub->nr].covmatrix(sub, value);
}



int check_covmatrix(cov_model *cov) {
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = PrevLoc(cov);
  int err, rows, cols,
    iso = SymmetricOf(cov->isoprev);
 
  ASSERT_LOC_GIVEN;
  if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, 
		     PosDefType, KERNEL, iso, SUBMODEL_DEP,
		   ROLE_GAUSS)) != NOERROR) { 
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, 
		     VariogramType, XONLY, iso, SUBMODEL_DEP,
		     ROLE_GAUSS)) != NOERROR) {
      return err;
    }
  }
  setbackward(cov, sub);  
  rows = cov->vdim[0] = sub->vdim[0]; 
  cols = cov->vdim[1] = sub->vdim[1]; 

  if (cov->q == NULL) {
    QALLOC(2);
    cov->q[0] = loc->totalpoints * rows;
      cov->q[1] = loc->totalpoints * cols;
  }

  int dim = loc->timespacedim;
  assert(dim == cov->tsdim);
  if ((err = alloc_cov(cov, dim, rows, cols)) != NOERROR) return err;

  //  printf("%s=%d %s=%d \n\n\n\n\n\n", NAME(cov), cov->Spgs != NULL, NAME(cov->sub[0]), cov->sub[0]->Spgs != NULL); 
 
  return NOERROR;
}



void Pseudovariogram(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  CovList[sub->nr].pseudovariogram(sub, value);
}


void Variogram(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;

  CovList[sub->nr].variogram(sub, value);
}

int check_vario(cov_model *cov) {
  return check_cov_intern(cov, VariogramType, GLOBAL.general.vdim_close_together,
			  true);
}


int struct_variogram(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){
  int err;
  cov_model *sub,
    *next = cov->sub[0];
  location_type *loc = Loc(cov);
 
  if ((sub = get_around_max_stable(next)) == next) sub = get_around_gauss(next);
  if (sub != next) {
    if ((err = covCpy(&(cov->key), sub)) != NOERROR) return err;       
    sub = cov->key;
    sub->calling = cov;
    if (!isVariogram(sub->typus)) SERR("variogram model cannot be determined");
  } else {
    if (!isVariogram(sub->typus)) SERR("not a variogram model");
  }

  if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, VariogramType,
		     loc->y == NULL && loc->ygr[0] == NULL ? XONLY : KERNEL,
		     SYMMETRIC, cov->vdim,
		     ROLE_COV)) != NOERROR) {
    return err;
  }   
  return NOERROR;
} 
  
   
// bool close = GLOBAL.general.vdim_close_together;

void Fctn(double VARIABLE_IS_NOT_USED *X, cov_model *cov, double *value) {
  cov_model *sub = cov->sub[0];
  if (sub == NULL) BUG;
  FctnIntern(cov, cov, sub, value, false);
}


void FctnIntern(cov_model *cov, cov_model *covVdim, cov_model *sub, double *value, bool ignore_y){ 
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array

  FINISH_START(cov, covVdim, ignore_y, 1);  
  double *y = ygiven ? pgs->supportmin : ZERO;
  int role = sub->role;
  if (sub->role == ROLE_COV) sub->role = ROLE_BASE;

  assert(sub != NULL); 
  assert(pgs != NULL);
  assert(cov->key == NULL);

  bool kernel = sub->domprev != XONLY;

  //INCLUDE_VAL;			

  long i, ii, v, j;							
  double **Val= pgs->Val,						
    *cross=pgs->cross;	
	
  //printf("vdimSq %d %d %d\n", vdimSq, cov->zaehler, pgs->rowscols);
  if (vdimSq > pgs->rowscols) {PMI(cov); BUG;} //
  if (vdimSq > 1) {							
    for (v = i = 0; i<vdimSqtot ; i+=vdim0tot) {			
      for (ii=i, j=0; j<vdim0; ii+=tot, v++, j++) {	
	Val[v] = value + ii;
      }									
    }									
  }
  
  if (kernel && !ygiven && PL > 0) {					
    char wrn[200]; 				
    sprintf(wrn, "'%s' is called with a single variable only, although it is used as a kernel. So, the second variable is set to zero, here.n", NICK(cov)); 
    warning(wrn);							
  }
 
#define FUNCTION FCTN(x, sub, value)
#define FUNCTION_Y NONSTATCOV(x, y, sub, value)
#define FUNCTION2 FCTN(x, sub, cross);		\
  for (v = 0; v<vdimSq; v++) { Val[v][loc->i_row] = cross[v];}
#define FUNCTION2_Y						\
  NONSTATCOV(x, y, sub, cross);				\
  for (v = 0; v<vdimSq; v++) Val[v][loc->i_row] = cross[v];	
 
  PERFORM(FUNCTION, FUNCTION2, FUNCTION_Y, FUNCTION2_Y);

  sub->role = role;
} 




int check_fctn(cov_model *cov) {
  EXTRA_STORAGE;
#define nTypes 2  
  Types T[nTypes] = {ShapeType, TrendType};
  int i, err;
  for (i=0; i<nTypes; i++) {
    err = check_fct_intern(cov, T[i], GLOBAL.general.vdim_close_together, 
			   true, 0, 0);
    if (err == NOERROR) return err;
  }
  return err;
}

 



/* ****************************** */
/*             PREDICT            */
/* ****************************** */
#define PREDICT_REGISTER 0

void predict(double VARIABLE_IS_NOT_USED *x, cov_model *predict, double *v) { 
  assert(predict != NULL && !PARAMisNULL(predict, PREDICT_REGISTER));
  cov_model 
    *cov = KEY[PARAM0INT(predict, PREDICT_REGISTER)],
    *sub = cov->key == NULL ? cov->sub[0] : cov->key;  
  assert(cov != NULL);
  if (v==NULL) {
    likelihood_storage *L = sub->Slikelihood;
    int store =  GLOBAL.general.set;
    GLOBAL.general.set = 0;
    listoftype *datasets = L->datasets;
    int
      vdim = cov->vdim[0],
      ncol =  NCOL_OUT_OF(datasets),
      repet = ncol / vdim;
    GLOBAL.general.set = store;
    assert(predict->qlen > 0 && cov->q != NULL);
    predict->q[predict->qlen - 1] = repet;
    return; // EvaluateModel needs information about size
   //                      of result array, given in cov->q
  }
   
   if (sub->nr == GAUSSPROC) {
     gauss_predict(predict, cov, v);
    return;
  }

  BUG;
}


int check_predict(cov_model *predict) {
  assert(predict != NULL);
  //PMI(predict);
  if (PARAMisNULL(predict, PREDICT_REGISTER)) SERR("'register; must be given.");
  cov_model *cov = KEY[PARAM0INT(predict, PREDICT_REGISTER)];
  location_type 
    *pred_loc = Loc(predict);
  cov_model 
    *sub = cov->key == NULL ? cov->sub[0] : cov->key;  
  int err;
  assert(sub->nr == GAUSSPROC);
  assert(pred_loc->delete_x);
  assert(pred_loc->timespacedim == Loc(cov)->timespacedim);
  assert(pred_loc->Time == Loc(cov)->Time);
  likelihood_storage *L = sub->Slikelihood;

  if (L == NULL || L->X == NULL)
    SERR1("'%s' not fully initialized", NICK(cov));

  if (cov == NULL || cov->nr != LIKELIHOOD_CALL || !cov->checked) 
    SERR1("'%s' not initialized", NICK(cov));

  if (pred_loc->y != NULL || pred_loc->ygr[0] != NULL) {
    if (predict->Sextra == NULL) // i.e. first call, so user's input
      SERR("set of y-values (kernal definition) not allowed");
  } else {  
    CONDCOV_NEW_STORAGE(predict, extra, a);
    assert(pred_loc->delete_y); // = true;
    if (pred_loc->grid) {
      int i,
	spatialdim = pred_loc->spatialdim,
	nr = spatialdim * 3;
      double *y = (double*) MALLOC(nr * sizeof(double));
      for (i=0; i<nr; i++) y[i] = 1.0;
      assert(pred_loc->ygr[0] == NULL);
      pred_loc->ly = 3;
      if ((err = setgrid(pred_loc->ygr, y, pred_loc->ly, spatialdim))!=NOERROR) 
	 return err;
      FREE(y);
      // assert(!pred_loc->Time); 
      if (pred_loc->Time) pred_loc->ygr[spatialdim] = pred_loc->T;
   } else {       
      pred_loc->ly = 1;
      // wichtig im folgenden tsdim nicht spatialdim
      pred_loc->y = (double *) MALLOC(pred_loc->timespacedim * sizeof(double));  
      pred_loc->T[XSTART] = pred_loc->T[XSTEP] = 0.0;
      pred_loc->T[XLENGTH] = 1;
    }
    assert(cov->Sextra == NULL);
  }

  err = check_fct_intern(predict, PosDefType, 
			 GLOBAL.general.vdim_close_together, true,
			 cov->vdim[0], 1);

  return err;
}

int struct_predict(cov_model *cov, cov_model VARIABLE_IS_NOT_USED  **newmodel){
  return struct_cov(cov, newmodel);
}


void range_predict(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[PREDICT_REGISTER] = 0;
  range->max[PREDICT_REGISTER] = MODEL_MAX;
  range->pmin[PREDICT_REGISTER] = 0;
  range->pmax[PREDICT_REGISTER] = MODEL_MAX;
  range->openmin[PREDICT_REGISTER] = false;
  range->openmax[PREDICT_REGISTER] = false;
}
