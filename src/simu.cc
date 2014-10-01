/*
Authors
Martin Schlather, schlather@math.uni-mannheim.de

main library for unconditional simulation of random fields

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2014 Martin Schlather

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
//#include <sys/timeb.h>
 
#include <string.h>
#include "RF.h"
#include <unistd.h>
#include <R_ext/Utils.h>     
#include "Covariance.h"


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
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
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
      //printf("%s, len = %d\n", param_name, len);
      *qq = MALLOC(sizeof(int) * len);
      int * q = (int *) *qq;
      for (j=0; j<len; j++) q[j] = Integer(p, param_name, base + j);
    }
    break;
  case CLOSXP : 
    {
      error("Not programmed yet.\n");
      assert(false); // Marco: fctn das? JA! Das sind Funktionen.
      *qq = MALLOC(sizeof(SEXP));
      SEXP *q=(SEXP*) *qq;
      *q = p;
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
   case LANGSXP : 
     {
       if (strcmp("setseed", param_name) != 0 && strcmp("env", param_name)!=0){ 
	 //      printf("here\n");
	 if (GLOBAL.general.storing) {
	   char msg[200];
	   sprintf(msg, "If models with R commands in the parameters (such as '%s') are used then 'storing' must be FALSE.", CovList[USER].nick);
	   error(msg);
	 }
	 if (!GLOBAL.internal.warn_oldstyle) {
	   char msg[300];
	   sprintf(msg, "Models with R commands in the parameters (such as '%s') may not be called by obsolete functions.\nSee the notes in '?RMmodelsAdvanced' and set 'RFoldstyle(FALSE)'.", CovList[USER].nick);
	   error(msg);
	 }
       }
       
       *qq = MALLOC(sizeof(sexp_type));
       sexp_type *q = (sexp_type*) *qq;
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
      listoftype *q;
      bool list = TYPEOF(p) == VECSXP;
      if (!list && TYPEOF(p) != REALSXP) {
	PRINTF("type %d != %d", TYPEOF(p), REALSXP);
	BUG;
      }
      

      locallen = list ? len : 1;
      //     print("type=%d,%d %d %d %d\n", 
      //	    TYPEOF(p), VECSXP, locallen, MAXELEMENTS, type,length(p));
     
      if (locallen > MAXELEMENTS) PERR("too many list elements");
      *qq = MALLOC(sizeof(listoftype));
      q=(listoftype*) (*qq);
      q->deletelist = true;
      for (i=0; i<MAXELEMENTS; i++) {
	q->p[i] = NULL;
	q->ncol[i] = 0;
	q->nrow[i] = 0;
      }
      for (i=0; i<locallen; i++) {
	pi = list ? VECTOR_ELT(p, i) : p;


	//	print("%s %d type=%d,%d ncol/nrow %d %d %d %d\n ", 
	//     param_name, i, TYPEOF(p), VECSXP,
	//	ncols(pi), nrows(pi), isMatrix(pi),  pi == p);

	includeparam((void**) (q->p + i), REALSXP, length(pi), 
		     pi, base, param_name); 
       
	
	if (isMatrix(pi)) {
	  // Achtung isVector ist true falls isMatrix true
	  q->ncol[i] = ncols(pi);
	  q->nrow[i] = nrows(pi);	
	} else if (isVector(pi)) {
	  q->ncol[i] = 1;
	  q->nrow[i] = length(pi);
	} else  {
	  PERR("list element(s) neither vector nor matrix");
	}
      }
    }
    break;
  default : PERR("unmatched internal type of parameter");
  }
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
  bool subleft[MAXSUB]; //, True=true, False=false;


  // print("hier\n");

  if (TYPEOF(model) != VECSXP) { // not list8
    ERR("list expected")
      }
  PROTECT(m = VECTOR_ELT(model, elt++));
  if (!isString(m)) {
    ERR("first element must be the name of a covariance model"); 
  }
  if (length(m) != 1) ERR("model name must be a single word");
  strcopyN(name, (char*) CHAR(STRING_ELT(m, 0)), MAXCHAR);

  // printf("%d\n", (2 * level < NLEER) ? 2 * level : NLEER); assert(false);

  strcopyN(leer, "                                                           ", 
	   (2 * level < NLEER) ? 2 * level : NLEER);
  sprintf(ERR_LOC, "%s\n%s%s... ", ERROR_LOC, leer, name);
  strcpy(ERROR_LOC, ERR_LOC);
  covnr = getmodelnr(name);

  //  if (covnr == NATSC && NS) ERR("natsc model and RFparameters(PracticalRange=TRUE) may not be given at the same time");
  if (covnr == -2) { ERR("multiple matching of covariance name") }
  else if (covnr == -1) { 
    ERR("unknown covariance model.  'RFgetModelNames()' yields the list of models");
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

      //      printf("i=%d %s %s kappas=%d %s\n", i, name, param_name, nkappas,  C->kappanames[nkappas - 1]);
      
 
      if (i<0 && !strcmp(C->kappanames[nkappas - 1], FREEVARIABLE)) {
	//print("free variable\n");
	if ((j = Match(param_name, C->subnames, C->maxsub)) < 0) {
	  j = Match(param_name, STANDARDSUB, MAXSUB);
	}
	if (j < 0) {
	  if (cov->ownkappanames == NULL) {
	    cov->ownkappanames = (char**) CALLOC(nkappas, sizeof(char*));
	  }
	  
	  // ein FREEVARIABLE-Zeichen wird nie verwendet -- Puffer fuer einfachen Algorithmus; deshalb "i-1"
	  i = nkappas - 1; 
	  
	  //PMI(cov);
	  
	  //printf("i=%d\n", i);
	  //	  printf("%s %d\n", C->kappanames[i - 1]
	  //		 , !strcmp(C->kappanames[i - 1], FREEVARIABLE));
	  //	  printf("i=%d %s\n", 0, cov->ownkappanames[i] );
	  
	  while(i>0 && !strcmp(C->kappanames[i - 1], FREEVARIABLE) && 
		cov->ownkappanames[i] != NULL) {
	    //printf("i=%d %s %ld\n", i, C->kappanames[i - 1],  
	    //	   cov->ownkappanames[i] );
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
      //print("nothing\n");
      i = NOMATCHING;
    }
    //  printf("%s %s i=%d %d %d \n", name, param_name, i, isVectorList(p), isVectorList(p) && isString(VECTOR_ELT(p, 0)));

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

#define MMERR(X) {sprintf(MSG,"%s\n%s: %s",ERROR_LOC, param_name, X);	\
	      error(MSG);}

            MMERR("multiple matching of submodel names")
	      } else {
            MMERR("unmatched covariance name for submodel")
	      }
	}
      } else {
	for (j=0; j < C->maxsub; j++) {
	  if (subleft[j]) {
	    i = j;
	    break;
	  }
	}
	if (j == C->maxsub) {
	  // PRINTF("%d %d\n", j, C->maxsub);
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
      int l = length(p);
      //      if (param_name[0]==ONEARGUMENT_NAME && strlen(param_name) == 1) { 
      //     if (strcmp(param_name, "k") == 0) {
      if (param_name[0] == ONEARGUMENT_NAME && param_name[1] ==  '\0') {
	if (TYPEOF(p) != REALSXP && TYPEOF(p) != INTSXP && TYPEOF(p) != LGLSXP) 
	  PERR1("if '%c' is used as parameter name then only a numerical vector is allowed as value", ONEARGUMENT_NAME);
	if (l != C->kappas) 
	  PERR1("length of vector does not match number of parameters. Do not use short name '%c' in complex situtations either.", ONEARGUMENT_NAME);
	// if (C->maxsub!=0)
	// ERR("short form k of parameter name not allowed for sophistacted models");
	for (j=0; j<l; j++) {
	  if (!PisNULL(j)) PERR("parameter given twice"); // p[0] OK
	  if (true ||  // neu seit 6.8.14, wegen RFgui, modelParam = .Call("SetAndGetModelInfo",
	      C->kappatype[j] != INTSXP ||
	      Integer(p, param_name, j) != NA_INTEGER) {
	    cov->ncol[j] = cov->nrow[j] = 1;	  

	    if ( C->kappatype[j] != INTSXP && C->kappatype[j] != REALSXP)
	      PERR("model has too complicated arguments to allow abbreviations");
	    
	    // printf("include %s %d %s\n", Nick(cov), j, KNAME(j));
	    includeparam((void**) (cov->px + j), C->kappatype[j], 1, 
			 p, j, param_name);
	  }
	}
	continue;
      }
      
      if (strcmp(param_name, "") == 0) {
	//printf("XX name %s\n", NICK(cov));
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
  
      //     printf("model %d %s %d %s %d %d = %f\n", elt, C->nick, i, C->kappanames[i], cov->p[i] != NULL, cov->kappasub[i] != NULL, Real(p, C->kappanames[i], 0));

       
      if (isVectorList(p) && isString(VECTOR_ELT(p, 0))) {
	PFREE(i);
	CMbuild(p, level + 1, cov->kappasub + i);
	strcpy(ERROR_LOC, ERR_LOC);
	cov->kappasub[i]->calling = cov;
      } else {
	// first the fixed one for the dimensions !!
	if (!PisNULL(i) || cov->kappasub[i] != NULL)
	  ERR("parameter given twice");

	//printf("Include %s %d %s\n", Nick(cov), i, KNAME(i));
	includeparam((void**) (cov->px + i), C->kappatype[i], l, 
		     p, 0, param_name);
	
	//    print("danach %d %d %s %ld %ld\n", i, C->kappatype[i], param_name, cov->p[i],    cov->p[i+1]);
		
	//    print(">>>>>    %s %d\n", param_name, isEnvironment(p));        	
	
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
	    cov->nrow[i] = length(p);
	  } else if (isFunction(p) || isLanguage(p) || isEnvironment(p)) {
	    cov->ncol[i] = 1;
	    cov->nrow[i] = 1;	  
	  } else  {
	    //print("type=%d\n", TYP
	    PERR("neither vector nor matrix");
	  }
	}
      }
    }
  }

  for (j=0; j < C->maxsub; j++) if (subleft[j]) break;
  if (j < C->maxsub) {
    bool shape_process = isMaxStable(cov) || 
      covnr == RANDOMCOIN_USER || cov->nr == AVERAGE_USER;
    if (j < C->minsub) { 
      if (!shape_process || subleft[C->maxsub - 1]) {
	//APMI(cov);
	cov_model *prev = cov->calling;
	while (prev != NULL && prev->nr != RFGET) prev = prev->calling;
	if (prev == NULL) ERR("not enough submodels");      
      }
      if (PL >= PL_SUBIMPORTANT && !shape_process)  {
	// da darf von 2 untermodellen nur 1 angegeben werden
	warning("not all submodels are given");
      }
    }
  }
  
  //   print("ns=%d %d %d %s\n", NS, NATSCALE_MLE, C->inverse != NULL, C->name);
  cov_model *next = cov->sub[0];
  if (next != NULL) {
    cov_fct *N = CovList + next->nr;
    if (NS > 0 && N->maxsub == 0 && N->inverse != ErrInverse &&
	N->cov != nugget &&  N->cov != constant) {
      // printf("here\n");
      if (NS == NATSCALE_MLE && isDollar(cov)) {
	// natsc wird nur dann eingeschoben, wenn davor scale=NA war
	bool natsc= !PisNULL(DSCALE) && ISNAN(P0(DSCALE));
	
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
	//	printf("aa %s %s\n", NICK(cov), cov->calling == NULL ? "none" : NICK(cov->calling));
	addModel(Cov, NATSC_USER); 
      }
    }
  }
  UNPROTECT(1);
}


/*

void CheckModelInternal(SEXP model, SEXP x, SEXP y, SEXP T, 
			int spatialdim, // spatial dim only !
			int xdimOZ,
			int *lx, 
			//			int domprev,
			bool grid,
			bool distances,
			bool Time, 
			cov_model **Cov) {
  //  int err = NOERROR,
  //    //  len = length(x),
  //tsdim = spatialdim + (int) Time;

  error("currently 'x' may not be a list");
   
}
*/


bool CallingSet(cov_model *cov) {
  int i;
  for (i=0; i<cov->nsub; i++) {
    cov_model *sub = cov->sub[i];
    if (sub == NULL) {
      if (CovList[cov->nr].range != range_randomcoin) {
	PMI(cov, "Calling point A error"); //
	return false;
      }
    } else {
      if (sub->calling != cov) {
	PMI(cov, "Calling point B error"); //
	return false;
      }
      if (!CallingSet(sub)) {  
	return false;
      }
    }
  }
  if (cov->key != NULL && !CallingSet(cov->key)) return false;
  if (cov->Splus != NULL) {
    for (i=0; i<cov->nsub; i++) {
      cov_model *sub = cov->Splus->keys[i];
      if (sub == NULL) {
	PMI(cov, "plus A error");//
	return false;
      } else {
	if (sub->calling != cov) {
	  PMI(cov, "plus B error"); //
	  return false;
	}
	if (!CallingSet(sub)) return false;
      }
    }
  }
  return true;
}

void CheckModelInternal(SEXP model, double *x, double *Y, double *T, 
			int spatialdim, /* spatial dim only ! */
			int xdimOZ,
			int Lx, int Ly,
			bool grid,
			bool Distances,
			bool Time, 
			cov_model **Cov) {
  int err,
    zaehler = 0,
    tsdim = spatialdim + (int) Time,
    ly = Ly;
  double *y = Y;
  cov_model *cov;
  char EM2[500 + 100 * Nothing] = "";

 
  GetRNGstate(); // some parts use monte carlo to determine internal values
  if (currentNrCov==-1) InitModelList();  

  // printf("start checkmodelinternal\n");
 
  while (true) {
    strcpy(ERROR_LOC, "Building the model:");
    cov = NULL;

    if (*Cov != NULL) COV_DELETE(Cov);
    assert(*Cov == NULL);
    CMbuild(model, 0, Cov);

    strcpy(ERROR_LOC, "Having built the model:");
    cov = *Cov;
    if (cov == NULL) ERR("no model is given");
    assert(cov->calling == NULL);

    
    int lx = Lx;
    bool distances = Distances;
    if (distances) {
      if (cov->nr == COVMATRIX || CovList[cov->nr].check == check_simulate) {
	lx = (int) (1e-9 + 0.5 * (1 + sqrt(1. + 8 * lx)));
	if (Lx != lx * (lx - 1) / 2)
	  ERR("distance length not of form 'n * (n - 1) / 2'");
      } else {
	distances = false;
      }
    } else if (xdimOZ < spatialdim) {
      PRINTF("xdimOZ=%d, spatialdim=%d Time=%d\n", xdimOZ, spatialdim, Time);
      GERR("dimension of 'x' does not fit the space-time dimension");
    }

    cov->tsdim = spatialdim + (int) Time;
    cov->xdimprev = xdimOZ;
    //  cov->domprev = ly == 0 ? XONLY : KERNEL;
    cov->domprev = XONLY; // formal sind alle Interface Modelle nur von 
    //                              einer (dummy) Variablen abhaengig
   
    switch (GLOBAL.coords.coord_system) {
    case coord_auto: 
    case cartesian: cov->isoprev = CARTESIAN_COORD; break;
    case earth: 
      if (spatialdim != 2) 
	GERR1("an earth coordinate consists of two angles. Got %d components.",
	      spatialdim);
      cov->isoprev =  EARTH_COORD; 
      break;
    default: GERR("unknwon coordinate system");
    }



    if (xdimOZ < spatialdim) {
      if (xdimOZ != 1) GERR("reduced dimension is not one");
      //if (isCov(cov)) cov->isoprev = ISOTROPIC;
      //assert(false);
    } else if (xdimOZ > spatialdim) {
      PRINTF("xdimOZ=%d, spatialdim=%d Time=%d\n", xdimOZ, spatialdim, Time);
      GERR("mismatch of dimensions");
    }
    cov->calling=NULL;
    cov->prevloc = NULL;
    
    if ((err = loc_set(x, y, T, spatialdim, xdimOZ, lx, ly, Time, grid, 
		       distances, &(cov->prevloc))) != NOERROR)
      goto ErrorHandling;

   strcpy(ERROR_LOC, "Checking the model:");
    if (PL >= PL_DETAILS) {
      //   assert(cov->isoprev != 6);
      PMI(cov, "\nCheckModelInternal, before checking:");//OK
    }

    if ((err = CHECK(cov, tsdim, xdimOZ + Time, InterfaceType,	 
		       cov->domprev, cov->isoprev, 
		       SUBMODEL_DEP, ROLE_BASE))
	!= NOERROR) {
     ////
      if (PL >= PL_ERRORS) {
	PMI(cov, "CheckModelInternal has failed:"); //
	PRINTF("err =%d\n", err);
      }
      goto ErrorHandling;
    }

   if (PL >= PL_DETAILS) {
      PMI(cov, "CheckModelInternal, after checking:"); // OK
    }
 
    sprintf(ERROR_LOC, "%s process: ", ROLENAMES[cov->role]);
    strcpy(PREF_FAILURE, "");
    PrInL=-1;
    
    if (PL >= PL_DETAILS) PRINTF("CheckModelInternal A\n");
    
    // APMI(cov);

    assert(CallingSet(cov));
  
   err = STRUCT(cov, NULL);
 
    //   APMI(cov);

  ErrorHandling:

    //    printf("\n>> %d ly=%d zaehler=%d lx=%d\n\n", err, ly, zaehler, lx);

   
    if (err == NOERROR || zaehler>=1) {
      break;    
    } 

   //  printf("err=%dn",err);
    // APMI(cov);
    
    char EM[500];   
    errorMSG(err, EM);
 
    sprintf(EM2, "%s%s", PREF_FAILURE, EM);
    if (lx == 0 || distances) break;

    //   printf(">> %d '%s'\n", err, EM2);

    y = x;
    ly = lx;
    zaehler++;
  }

 
  //printf("end  checkmodelinternal\n");
  PutRNGstate();

  // AERR(err);
  //PMI(cov);
  // Fehlermeldung vom ersten Durchgang wird angezeigt.
  if (err != NOERROR) { 
    //     printf("%d '%s'\n", err, EM2);
    // APMI(cov);
    error(EM2); 
  }
  
  // muss ganz zum Schluss stehen, da Fehler zu einer unvollstaendigen
  // Struktur fuehren koennen
  assert(CallingSet(cov));
  //APMI(cov);
}


SEXP Init(SEXP model_reg, SEXP model,
	  SEXP x, SEXP y, SEXP T, 
	  SEXP spatialdim,
	  SEXP Grid, SEXP Distances, SEXP Time,
	  SEXP NA_OK) {
  // MAXSIMUDIM
  bool naok = (bool) LOGICAL(NA_OK)[0],
    grid = LOGICAL(Grid)[0],
    distances = LOGICAL(Distances)[0],
    time =  LOGICAL(Time)[0];
  int
    xdimOZ = grid ? ncols(x) : nrows(x),
    lx = grid ? 3 : ncols(x),
    ly = length(y) == 0 ? 0 : grid ? 3 : ncols(y);
  SEXP ans;
 
  currentRegister = INTEGER(model_reg)[0];

  NAOK_RANGE = naok; 
  PROTECT (ans = allocVector(INTSXP, 2));


  //  print("real x %d grid=%d %d time=%d lx=%d %d ncols nc=%d nr=%d  xdimOZ%d\n", isReal(x), grid, distances, time, lx, ly, ncols(x), nrows(x), xdimOZ);

  //  printf("key = %d\n", currentRegister);
 
  if (isReal(x))
    CheckModelInternal(model, REAL(x), REAL(y), REAL(T), INTEGER(spatialdim)[0],
		       xdimOZ, lx,  ly,
		       grid, distances, time,
		       KEY + currentRegister
		       );
  else {
    BUG;
    // todo not programmed yet
    /*
      CheckModelInternal(model, x, y, T, INTEGER(spatialdim)[0],
      xdimOZ, &lx,
      grid, distances, time,
      KEY + INTEGER(model_reg)[0]
      );
    */
  }
 

  cov_model *cov = KEY[currentRegister];
  NAOK_RANGE = false;
  if (PL >= PL_COV_STRUCTURE) {
    PMI(cov);// OK
  }

  INTEGER(ans)[0] = (cov)->vdim2[0];
  INTEGER(ans)[1] = (cov)->vdim2[1];
  UNPROTECT(1);
  //  APMI(KEY[INTEGER(model_reg)[0]]);

  return ans;
}



SEXP EvaluateModel(SEXP X, SEXP Covnr){
  if (currentNrCov==-1) InitModelList();  
  int d, mem, len,
    err = NOERROR;
  cov_model *cov = KEY[INTEGER(Covnr)[0]];
  SEXP result=NULL, 
    dummy=NULL;
  

  //  {double x[5]={1,2,3,6,5}; 
  //   return evaluate(x, 5);
  // }

  strcpy(ERROR_LOC, "");
  if (cov == NULL) GERR("register not initialised");
  if ( (len = cov->qlen) == 0) {
    BUG;
    if (cov->fieldreturn) GERR("model cannot be evaluated");
    mem = cov->vdim2[0] * cov->vdim2[1];
  } else {   
    //   
    //    printf("ende %f\n", REAL(X)[0]);

    //PMI(cov);

    FCTN(REAL(X), cov, NULL); 
    if (len > 1 && cov->q[len-1] == 1) len--; // last is for n=#{simulations}
    for (mem=1, d=0; d<len; d++) mem *= (int) cov->q[d];
  }

  //     printf("len =%d mem=%d %f %f %f=%f\n", len, mem, cov->q[0], cov->q[1], cov->q[2], *REAL(X));//assert(false);

  if (len == 1) PROTECT(result = allocVector(REALSXP, mem)); 
  else if (len == 2)  
    PROTECT(result = allocMatrix(REALSXP, cov->q[0], cov->q[1])); 
  else {
    PROTECT(dummy=allocVector(INTSXP, len));
    for (d=0; d<len; d++) {
      INTEGER(dummy)[d] = (int) cov->q[d];
      //    printf("q[%d]=%d %d\n", d, (int) cov->q[d], (int) cov->q[0] * (int) cov->q[1] * (int) cov->q[2]);
    }
    PROTECT(result=allocArray(REALSXP, dummy));
  }

  //   printf(">> %d %ld %d\n",  len, REAL(result), length(result)); 
  //   printf("len =%d  %f %f %f=%f\n", len, cov->q[0], cov->q[1], cov->q[2], *REAL(X));
  //APMI(cov);
  //PMI(cov, "last");
 
  GetRNGstate();
  FCTN(REAL(X), cov, REAL(result)); 
  PutRNGstate();

  //printf("done\n");
  // APMI(cov);

 ErrorHandling:
  if (result != NULL) UNPROTECT(1 + (int) (len > 2));
  if (err != NOERROR) XERR(err);

  return result;
}

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
  location_type *loc = cov->prevloc;
  long vdimtot; vdimtot = loc->totalpoints * cov->vdim2[0];
  assert(cov->vdim2[0] == cov->vdim2[1]);

  cov->q[cov->qlen - 1] = nn = (int) *N;
  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array

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
    
    if (GLOBAL.general.seed != NA_INTEGER && nn > 1) {
      if (PisNULL(SIMU_SEED) || PisNULL(SIMU_ENV)) {
	BUG;
      } else {
 	PutRNGstate();
	double simu_seed = GLOBAL.general.seed + (ni - 1);
	addVariable((char*) "seed", &simu_seed, 1, 1, PENV(SIMU_ENV)->sexp);
	eval(PSEXP(SIMU_SEED)->sexp, PENV(SIMU_ENV)->sexp);
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
    //

    //printf("here %d %s %d\n", ni, NICK(cov), sub->stor != NULL);
    // A    PMI(cov);
    // PSTOR(cov, cov->stor); // assert(false);
    // crash(cov);
    
    assert(cov->stor != NULL);
    DO(sub, cov->stor);

    if (sizeof(double) == sizeof(res_type) && false) {
      // printf("\n\n\n\n\n **********  %ld %ld !!\n", vdimtot, res);
      // APMI(cov);
      MEMCOPY(res, cov->rf, sizeof(res_type) * vdimtot);
      //printf("!!!  ********** !!\n");
    } else {
      //  printf("heredd \n");

      // printf("++ %d %ld %ld\n", vdimtot, res, cov->rf);

      int i; for (i=0; i<vdimtot; i++) {
	//printf("i=%d %f\n", i, cov->rf[i]);
	res[i] = cov->rf[i];
      }
      //  assert(false);
      //  	 APMI(cov);
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
  location_type *loc = cov->prevloc;
  int j, d, err, role, iso;
  Types type;
  bool vdim_close_together = GLOBAL.general.vdim_close_together;

  if (loc == NULL) SERR("locations not initialised.");
  kdefault(cov, SIMU_CHECKONLY, false);

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

    //printf("\n\nrole= %s\n", ROLENAMES[role]);

    for (j=0; j<=2; j++) {
      if ((TypeConsistency(type, sub) && 
	   (err = CHECK(sub, loc->timespacedim, cov->xdimown, type, 
			dom, iso, cov->vdim2, role)) == NOERROR) 
	  || isProcess(sub)) break;

      //printf("inconsistent '%s' %s '%s'\n", TYPENAMES[type], NICK(sub),
      //     TYPENAMES[CovList[sub->nr].Type]);

      if (j==0) type = NegDefType;
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
		     cov->vdim2, role)) != NOERROR) {
      return err;
    }
  }

  setbackward(cov, sub);  
  int subvdim = sub->vdim2[0];
  cov->vdim2[0]=subvdim; 
  cov->vdim2[1]=sub->vdim2[1]; 

  if (cov->q == NULL) {
    int len=1;
    if (loc->grid) len += loc->timespacedim; else len++;
    if (subvdim > 1) len++;
    
    cov->qlen = len;
    cov->q = (double *) MALLOC(sizeof(double) * len);
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

  //  APMI(cov);
  return NOERROR;
}


int struct_simulate(cov_model *cov, cov_model VARIABLE_IS_NOT_USED  **newmodel){
  cov_model *next = cov->sub[0],
    *sub = next;
  location_type *loc = cov->prevloc;
  int err,
    subrole = ROLE_FAILED,
    nr = next->nr;

  //APMI(next);

  if (isNegDef(next) || isTrend(next)) {
    if ((err = covcpy(&(cov->key), next)) != NOERROR) return err;
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;
    
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
		     XONLY, 
		     isCartesian(cov->isoprev) ? CARTESIAN_COORD : cov->isoprev, 
		     cov->vdim2, ROLE_GAUSS)) != NOERROR) {
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
    //PMI(cov, -1);
    //    PMI(sub, -1);
    //    printf("next %d %s %d %d %d\n", nr, NICK(next), isNegDef(next), isPosDef(next), CovList[nr].Type);
    ILLEGAL_ROLE;
  }
  sub->role = subrole; // ansonsten muesste hier C-HECK aufgerufen werden
  // hoffentlich gut gehende Abkuerzung, dass S-TRUCT aufgerufen wird,
  // und danach C-HECK (was auf jeden Fall gemacht werden muss)

  cov->simu.active = next->simu.active = false; 
  sub->simu.expected_number_simu = cov->simu.expected_number_simu;
  if (P0INT(SIMU_CHECKONLY)) return NOERROR;

  if (PL >= PL_DETAILS) PRINTF("Struct Simulate\n");
  //  print("%s %ld %ld %ld\n", NICK(sub), (long int) CovList[sub->nr].Struct, (long int)struct_gaussprocess, (long int)struct_spectral);
  //	A
  //  PMI(sub, "Internal"); assert(false);
  if ((err = STRUCT(sub, NULL)) != NOERROR) { return err; }
  if (PL >= PL_DETAILS) PRINTF("Checking Simulate\n");


  //  APMI(sub);
  assert(cov->stor == NULL);
  NEW_STORAGE(stor, STORAGE, gen_storage);

  if (!sub->initialised) {
    if (PL >= PL_DETAILS) PRINTF("Struct Simulate C\n");
   
    //   APMI(cov); //, "struct simu"); print("key=%ld\n", cov->key);
    
    if (//cov->key != NULL && // bei sub==next waere der falsche role gesetzt
	// irgenwie komisch, cov->key abzufragen und check(sub aufzurufen
	// aufgrund von Beispiel in RTpoisson geloescht
	(err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
		     cov->domprev, cov->isoprev, cov->vdim2,
		     subrole)) != NOERROR) {
      //
      // XERR(err);
      // APMI(sub);
      // assert(false);
      return err;
    }
    
    if (PL >= PL_DETAILS) {
      PRINTF("\n\nStruct Simulate (%s, #=%d), after 2nd check:",
	     NICK(sub), sub->gatternr); // OK
      PMI(sub); // OK
    }
 
    
    if ((err = INIT(sub, 0, cov->stor)) != NOERROR) {
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
  range->min[SIMU_CHECKONLY] = RF_NEGINF;
  range->max[SIMU_CHECKONLY] = RF_INF;
  range->pmin[SIMU_CHECKONLY] = -1e8;
  range->pmax[SIMU_CHECKONLY] = 1e8;
  range->openmin[SIMU_CHECKONLY] = true;
  range->openmax[SIMU_CHECKONLY] = true;
}

//void do_simulate(cov_model *cov, gen_storage *s){
//  assert(false);
//}

/////////////////////////////////////////


#define LIKELIHOOD_DATA 0
#define LIKELIHOOD_LEN 1
void likelihood(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v){ 
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  
  int j;
  double value,
    *data = P(LIKELIHOOD_DATA);
  location_type *loc = cov->prevloc;
  long vdimtot = loc->totalpoints * cov->vdim2[0],
    len = P0INT(LIKELIHOOD_LEN),
    repet = len / vdimtot;
  if (repet * vdimtot != len || repet == 0) 
    ERR("data and coordinates do not match");
   
  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  for (*v=0.0, j = 0; j < repet; j++, data += vdimtot) {
    VTLG_DLOG(data, sub, &value);
    *v += value;
  }
}


int check_likelihood(cov_model *cov) {
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = cov->prevloc;
  int j, err, role, iso;
  Types type;

  
  if (loc == NULL) SERR("locations not initialised.");

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
    if ((TypeConsistency(type, sub) && 
	 (err = CHECK(sub, loc->timespacedim, cov->xdimown, type, 
		      dom, iso, cov->vdim2, role)) == NOERROR) 
	|| isProcess(sub)) break;
    
    //printf("inconsistent '%s' %s '%s'\n", TYPENAMES[type], NICK(sub),
    //     TYPENAMES[CovList[sub->nr].Type]);
    
    if (j==0) type = NegDefType;
  }
  if (err != NOERROR) return err;

  setbackward(cov, sub);  
  cov->vdim2[0]=sub->vdim2[0]; 
  cov->vdim2[1]=sub->vdim2[1]; 

  if (cov->q == NULL) {
    cov->qlen = 1;
    cov->q = (double *) MALLOC(sizeof(double));
    cov->q[0] = 1;
  }

  return NOERROR;
}


int struct_likelihood(cov_model *cov, 
		      cov_model VARIABLE_IS_NOT_USED  **newmodel){
  cov_model *next = cov->sub[0],
    *sub = next;
  location_type *loc = cov->prevloc;
  int err,
    subrole = ROLE_FAILED,
    nr = next->nr;

  //APMI(next);
  //printf("PL %d \n", PL);
  //ssert(PL == 20);
  
  if (isNegDef(next)) {
    if ((err = covcpy(&(cov->key), next)) != NOERROR) return err;
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;
    
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
		     XONLY, 
		     isCartesian(cov->isoprev) ? CARTESIAN_COORD :cov->isoprev,
		     cov->vdim2, ROLE_GAUSS)) != NOERROR) {
      return err;
    }
    subrole = ROLE_GAUSS;    
  } else if (isBernoulliProcess(next)) subrole = ROLE_BERNOULLI;
  else if (isGaussBasedProcess(next)) subrole = ROLE_GAUSS;   
  else if (isBrownResnickProcess(next)) subrole = ROLE_BROWNRESNICK;    
  else if (nr == POISSONPROC) subrole = ROLE_POISSON;
  else if (nr == SCHLATHERPROC) subrole = ROLE_SCHLATHER;
  else if (nr == SMITHPROC) subrole = ROLE_SMITH;
  else {
    //PMI(cov, -1);
    //    PMI(sub, -1);
    //    printf("next %d %s %d %d %d\n", nr, NICK(next), isNegDef(next), isPosDef(next), CovList[nr].Type);
    ILLEGAL_ROLE;
  }

  if (CovList[sub->nr].logD == ErrCov)
    SERR1("log likelihood undefined for '%s'", NICK(sub));

  sub->role = subrole; // ansonsten muesste hier C-HECK aufgerufen werden
  // hoffentlich gut gehende Abkuerzung, dass S-TRUCT aufgerufen wird,
  // und danach C HECK (was auf jeden Fall gemacht werden muss)
 
   return NOERROR;
}



void range_likelihood(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[LIKELIHOOD_DATA] = RF_NEGINF;
  range->max[LIKELIHOOD_DATA] = RF_INF;
  range->pmin[LIKELIHOOD_DATA] = -1e8;
  range->pmax[LIKELIHOOD_DATA] = 1e8;
  range->openmin[LIKELIHOOD_DATA] = true;
  range->openmax[LIKELIHOOD_DATA] = true;

  range->min[LIKELIHOOD_LEN] = 1;
  range->max[LIKELIHOOD_LEN] = INFDIM;
  range->pmin[LIKELIHOOD_LEN] = 1;
  range->pmax[LIKELIHOOD_LEN] = INFDIM;
  range->openmin[LIKELIHOOD_LEN] = false;
  range->openmax[LIKELIHOOD_LEN] = false;
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
  } else if (!PisNULL(EVALDISTR_N)) {// r+
    xqp = P(EVALDISTR_N);
    for (j=i=0; i<n; i++, j+=dim) {
      //printf("j=%d %d \n", j, i); 
      VTLG_R(NULL, sub, v + j);    
    }
  } else BUG;
}


int check_EvalDistr(cov_model *cov) {
  cov_fct *C = CovList + cov->nr;
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  //  location_type *loc = cov->prevloc;
  int size_idx,  err, //type,  nr = sub->nr,
    *nrow = cov->nrow,
    *ncol = cov->ncol,
    dim = cov->tsdim, 
    zaehler=0;

  ROLE_ASSERT(ROLE_DISTR);
 
  if (cov->q == NULL) {
    cov->qlen = 1; // !! 1 obwohl 2 reserviert werden !! Wg EvaluateModel
    if (dim > 1 && 
	((!PisNULL(EVALDISTR_N) && P0(EVALDISTR_N) > 1) ||
	 (!PisNULL(EVALDISTR_Q) && P0(EVALDISTR_Q) > 1)))
      cov->qlen++;
    cov->q = (double *) MALLOC(sizeof(double) * (cov->qlen + 1));
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
	cov->q[size_idx] = nrow[EVALDISTR_P];
 	zaehler++;
      } 
    }// kein else wegen zaehler !!

    if (!PisNULL(EVALDISTR_N)) { // r      
      cov->q[size_idx] = P0(EVALDISTR_N);
      zaehler++;
    } 
    if (zaehler != 1) SERR("exactly one of the parameters must be given");
    //printf("+ %d, %d -> %s\n", (int) cov->q[EVALDISTR_Q_TYPE], sub->nr,
    //	   NICK(sub));  
    //    printf("= %d, %d -> %s\n", (int) cov->q[EVALDISTR_Q_TYPE], sub->nr,
    //	   NICK(sub));
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
  //  location_type *loc = cov->prevloc;
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
  
  assert(cov->stor == NULL);
  NEW_STORAGE(stor, STORAGE, gen_storage);

  if ((err = INIT(sub, 0, cov->stor)) != NOERROR) {
    //APMI(cov); // !!! ?? hier weitermachen
    return err; 
  }
  // assert(false);

  // APMI(cov);

  if (cov->rf == NULL) {
    int size = cov->q[0];
    if (cov->qlen > 1) size *= cov->q[1];

    if ((cov->rf = (res_type*) MALLOC(sizeof(res_type) * size)) 
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
  BUG; //error("unexpected call of dummy function");
}


int check_dummy(cov_model *cov) {
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = cov->prevloc;
  int i, err;

  if (loc == NULL) SERR("locations not initialised .");

  for (i=0; i<=1; i++) {
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, NegDefType,
		     i==0 ? XONLY : KERNEL, SYMMETRIC, SUBMODEL_DEP,
		     ROLE_COV)) == NOERROR) break;
  }
  if (err != NOERROR) return err;
  setbackward(cov, sub);  
  cov->vdim2[0] = sub->vdim2[0]; 
  cov->vdim2[1] = sub->vdim2[1]; 

  return NOERROR;
}


int struct_dummy(cov_model VARIABLE_IS_NOT_USED *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){    
  return NOERROR;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

cov_model *get_around_gauss(cov_model *cov) {
  cov_model *nxt= cov;
  if (nxt->nr == SCHLATHERPROC) nxt = nxt->sub[0]; // kein else 
  if (nxt->nr == GAUSSPROC) nxt = nxt->sub[0];
  if (isGaussMethod(nxt) || isBernoulliProcess(nxt)) {
    if (nxt->nr == AVERAGE_USER){
      if (nxt->sub[0] == NULL) error("covariance cannot be calculated (yet) for arbitrary shape functions.");
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
      if (nxt->sub[0] == NULL) error("covariance cannot be calculated (yet) for arbitrary shape functions.");
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


int alloc_pgs(cov_model *cov) {
  if (cov->xdimown != cov->tsdim) BUG;
  return alloc_pgs(cov, cov->xdimown);
}

int alloc_pgs(cov_model *cov, int dim) { // all what is necessary for dompp
  pgs_storage *pgs = NULL;

  assert(cov->Spgs == NULL); // NIE PGS_DELETE(&(cov->Spgs)) in Huetchen, da sonst dompp durcheinander kommt;
  NEW_STORAGE(Spgs, PGS, pgs_storage);
  pgs = cov->Spgs;

  //    printf("here dim=%d\n", dim);
  //assert(dim == 2); 
  //assert(false);
  
  if ((pgs->supportmin = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->supportmax = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->own_grid_start = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->own_grid_step = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->own_grid_length = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->supportcentre = (double*) CALLOC(dim, sizeof(double))) == NULL || 
      
      (pgs->gridlen = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->end = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->start = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->delta = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->nx = (int*) CALLOC(dim, sizeof(int))) == NULL ||
   
      (pgs->xstart = (double*) CALLOC(dim, sizeof(double))) == NULL || 
      (pgs->x = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->inc = (double*) CALLOC(dim, sizeof(double))) == NULL)
    return ERRORMEMORYALLOCATION;

  // printf("here end\n");
  return NOERROR;
}


int alloc_cov(cov_model *cov, int dim, int rows, int cols) { // all what is necessary for dompp

  int err;

  if (cov->Spgs != NULL) PGS_DELETE(&(cov->Spgs));
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

  // printf("here end\n");
  return NOERROR;
}


void Cov(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  CovList[sub->nr].covariance(sub, value);
}

int check_cov_intern(cov_model *cov, Types type, bool close) {

  cov_model 
    *next = cov->sub[0],
    *sub = cov->key == NULL ? next : cov->key;
  location_type *loc = Loc(cov);
  if (loc == NULL) SERR("locations not initialised!");  

  int err,  i, iso, 
    dim =  loc->timespacedim;
  assert(dim == cov->tsdim);
  bool  
    ygiven = loc->ly > 0;

  

  if (!ygiven && isNegDef(type) && !isCartesian(cov->isoown)) {
    add_y_zero(loc);
    ygiven = true;
    //PMI(cov);
  }
  //  PMI(cov);

  //   printf("ygiven %d %d %d\n", !ygiven, isNegDef(type), !isCartesian(cov->isoown));
    
  iso = type == ShapeType ? CARTESIAN_COORD : SYMMETRIC;
  
  //PMI(cov);
  //printf("type = %s\n", TYPENAMES[type]);

  for (i=0; i<=1; i++) {
    if ((err = CHECK(sub, dim, cov->xdimown, type, 
		   i==0 ? XONLY : KERNEL, iso, SUBMODEL_DEP,
		   sub!=next || isNegDef(sub) ? ROLE_COV : ROLE_BASE))
	== NOERROR) break;
  }
  if (err != NOERROR) return err;
  setbackward(cov, sub);  

  //PMI(cov);
  //  printf("%ld %ld\n", CovList[sub->nr].check, checkcox);
  //APMI(cov);
 

  if (sub->pref[Nothing] == PREF_NONE) SERR("given model cannot be evaluated")
  
  if (cov->q == NULL) {
    int d,
      len=1; // # of "simulations"
    if (loc->grid) len += dim; else len ++;      
    for (i=0; i<2; i++) len += (int) (cov->vdim2[i] > 1); 
    cov->qlen = len;  
    cov->q = (double *) MALLOC(sizeof(double) * len);

#define VDIMS  for (i=0; i<2; i++) if (cov->vdim2[i] > 1) cov->q[d++] = cov->vdim2[i]
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

  if ((err = alloc_cov(cov, dim, cov->vdim2[0], cov->vdim2[1])) != NOERROR) return err;

  return NOERROR;
}



int check_cov(cov_model *cov) {
  return check_cov_intern(cov, PosDefType, GLOBAL.general.vdim_close_together);
}


int struct_cov(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){
  int err;
  cov_model *sub, 
    *next = cov->sub[0];
  location_type *loc = cov->prevloc;
 
  sub = get_around_gauss(next);
  // fehlt : Poisson

  if (sub != next) {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) return err;   
    sub = cov->key;
    if (!isPosDef(sub->typus)) SERR("covariance model cannot be determined");
    
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, 
		     PosDefType, 
		     loc->y == NULL && loc->ygr[0] == NULL ? XONLY : KERNEL, 
		     SYMMETRIC, cov->vdim2,
		     ROLE_COV)) != NOERROR) {
    return err;
    }
  }
  return NOERROR;
}


void CovMatrix(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  CovList[sub->nr].covmatrix(sub, value);
}



int check_covmatrix(cov_model *cov) {
  cov_model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = cov->prevloc;
  int err, rows, cols;
  
  if (loc == NULL) SERR("locations not initialised");
  if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, 
		     PosDefType, KERNEL, SYMMETRIC, SUBMODEL_DEP,
		     ROLE_GAUSS)) != NOERROR) { 
    if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, 
		     NegDefType, XONLY, SYMMETRIC, SUBMODEL_DEP,
		     ROLE_GAUSS)) != NOERROR) {
      return err;
    }
  }
  setbackward(cov, sub);  
  rows = cov->vdim2[0] = sub->vdim2[0]; 
  cols = cov->vdim2[1] = sub->vdim2[1]; 

  if (cov->q == NULL) {
    int len=2;
    //    if (loc->grid) len += 2;      
    // if (vdim > 1) len += 2;
    cov->qlen = len;  
    cov->q = (double *) MALLOC(sizeof(double) * len);
    // printf("len %d\n", len);
    cov->q[0] = loc->totalpoints * rows;
      cov->q[1] = loc->totalpoints * cols;
  }

  int dim = loc->timespacedim;
  assert(dim == cov->tsdim);
  if ((err = alloc_cov(cov, dim, rows, cols)) != NOERROR) return err;
  
 
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
  //  printf("Variogram %s %ld %ld\n", NICK(sub),
  //	 CovList[sub->nr].variogram, CovList[sub->nr].covariance);
  // assert(false);
  CovList[sub->nr].variogram(sub, value);
}

int check_vario(cov_model *cov) {
  return check_cov_intern(cov, NegDefType, GLOBAL.general.vdim_close_together);
}


int struct_variogram(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel){
  int err;
  cov_model *sub,
    *next = cov->sub[0];
  location_type *loc = Loc(cov);
  // APMI(cov);
  // assert(loc != NULL);
 
  if ((sub = get_around_max_stable(next)) == next) sub = get_around_gauss(next);
  if (sub != next) {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) return err;       
    sub = cov->key;
    sub->calling = cov;
    if (!isNegDef(sub->typus)) SERR("variogram model cannot be determined");
  } else {
    if (!isNegDef(sub->typus)) SERR("not a variogram model");
  }

  // APMI(sub);

  if ((err = CHECK(sub, loc->timespacedim, cov->xdimown, NegDefType,
		     loc->y == NULL && loc->ygr[0] == NULL ? XONLY : KERNEL,
		     SYMMETRIC, cov->vdim2,
		     ROLE_COV)) != NOERROR) {
    return err;
  }
  return NOERROR;
}




int check_fctn(cov_model *cov) {
  return check_cov_intern(cov, ShapeType, true);
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
	  s->vdim2[0] = orig->nrow[i];
	  s->vdim2[1] = orig->ncol[i];
	  s->size = s->vdim2[0] * s->vdim2[1];
	} else {
	  int k;
	  s->size = s->vdim2[0] = cov->nrow[i];
	  s->vdim2[1] = cov->ncol[i];
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
  NEW_STORAGE(Sget, GET_STORAGE, get_storage);
  get_storage *s = cov->Sget;

  if ((err = SearchParam(cov, s)) != NOERROR) return err;
  for (i=0; i<2; i++) cov->vdim2[i] = s->vdim2[i];
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
  
  NEW_STORAGE(Sget, GET_STORAGE, get_storage);
  get_storage *s = cov->Sget;
  if ((err = SearchParam(cov, s)) != NOERROR) return err;
  for (i=0; i<2; i++) {
    if (cov->vdim2[i] != s->vdim2[i])
      SERR("mismatch of dimensions when constructing the model");
  }
  cov->fieldreturn = cov->origrf = false;  
  return NOERROR;
}




//void do_Rfget(cov_model *cov, gen_storage *s){
//  assert(false);
//}


