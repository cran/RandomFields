/*
Authors
Martin Schlather, martin.schlather@math.uni-goettingen.de

main library for unconditional simulation of random fields

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2011 Martin Schlather

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
 xdim       : dimension of the points given by the user.
              value is <= timespacedim since in case of isotropy only
	      the distances might be given by the user
 timespacedim : (= kc->dim = key->timespacedim)
              the true dimension for the location (if necessary by 
	      explicite parameter, e.g. in CovarianceFct.)
 anisodim   : vector of dimensions given by the user in the definition 
              of the covariance function (by the dimension of the anisotropy 
	      matrix); could be less than timespacedim in case of isotropic 
	      covariance functions;
	      the value is always less than or equal to timespacedim;
	      it can be less than timespacedim only in submodels of hypermodels
*/


SEXP getListElement(SEXP list, char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  if (names == R_NilValue) {
    return R_NilValue;
  }
  for (i = 0; i < LENGTH(names); i++) {
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
    n=LENGTH(names);
  
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
	      //assert(false);
	      return MULTIPLEMATCHING;
	    }
	  return j;
	}
	else {multiplematching=true;}
      }
      j++;
    }
    if (multiplematching) {
      // assert(false);
      return MULTIPLEMATCHING;
    }
  } else return NOMATCHING;
  return Nr;
}



void fetchParam(cov_model *cov, cov_model *next, int i, char *name) {
  if (next->p[i] != NULL) {
    if (next->ncol[i] != 1 || next->nrow[i] != 1) {
      char MSG[255];
      sprintf(MSG, "%s is not a scalar", name);
      ERR(MSG);
    }
    if (cov->p[i] == NULL) {
      cov->p[i] = (double*) malloc(sizeof(double));
      cov->ncol[i] = cov->nrow[i] = 1;
      cov->p[i][0] = next->p[i][0];
      cov->ncol[i] = cov->nrow[i] = 1;
    } else cov->p[i][0] *= next->p[i][0];
  } 
}

void includeparam(void **qq,  SEXPTYPE type, int len, SEXP p, int base,
		  char *param_name) {
  int j;
  switch(type) {
      case REALSXP : 
      {
	double *q;
	*qq = malloc(sizeof(double) * len); 

	q = (double *) *qq;
	switch(TYPEOF(p)) {
	    case REALSXP : 
	      for (j=0; j<len; j++) q[j] = REAL(p)[base+j];
	      break;
	    case INTSXP :
	      for (j=0; j<len; j++) 
		  q[j] = INTEGER(p)[j+base]==NA_INTEGER  
		      ? NA_REAL : (double) INTEGER(p)[base+j];
	      break;
	    case LGLSXP :
	      for (j=0; j<len; j++) {
		q[j] = LOGICAL(p)[j+base]==NA_LOGICAL ? NA_REAL 
		  : (double) LOGICAL(p)[base+j];
	      }
	      break;		    
	    default :
	      PERR("unmatched type of parameter");
	}
      }
      break;
      case INTSXP :
      {
	int *q;
	*qq = malloc(sizeof(int) * len);
	q = (int *) *qq;
	switch(TYPEOF(p)) {
	    case INTSXP : 
	      for (j=0; j<len; j++) q[j] = INTEGER(p)[j+base]; 
	      break;
	    case REALSXP : 
	      for (j=0; j<len; j++) {
	      double value;
	      value = REAL(p)[j+base];
	      if (ISNAN(value)) 
		PERR("NAs not allowed for integer valued parameters");
	      if (value == trunc(value)) q[j] = (int) value; 
	      else PERR("integer value expected");
	      }
	      break;
	    case LGLSXP :
	      for (j=0; j<len; j++) {
		q[j] = LOGICAL(p)[j+base]==NA_LOGICAL ? NA_INTEGER
		  : (int) LOGICAL(p)[j+base];
	      }
	      break;		    
	    default : PERR("unmatched type of parameter");
	}
      }
      break;
      case LISTOF + REALSXP : // list	 
      //  vector and matrix are turned into a list of 1 matrix	  
      {
	  assert(base == 0);
	int i, locallen;
	SEXP pi;
	listoftype *q;

	locallen = (TYPEOF(p) == VECSXP) ? len : 1;
//	printf("%d %d %d\n", locallen, MAXELEMENTS, type);

        if (locallen > MAXELEMENTS) PERR("to223o many list elements");
        *qq = malloc(sizeof(listoftype));
        q=(listoftype*) (*qq);
	for (i=0; i<MAXELEMENTS; i++) {
	  q->p[i] = NULL;
	  q->ncol[i] = 0;
	  q->nrow[i] = 0;
	}
        for (i=0; i<locallen; i++) {
	  pi = (TYPEOF(p) == VECSXP) ? VECTOR_ELT(p, i) : p;
	  includeparam((void**) (q->p + i), REALSXP, LENGTH(pi), pi, base,
		       param_name); 
//	  printf("ncol/nrow %d %d %d\n ", ncols(pi), nrows(pi), isMatrix(pi));
	  
	  if (isMatrix(pi)) {
	    // Achtung isVector ist true falls isMatrix true
	    q->ncol[i] = ncols(pi);
	    q->nrow[i] = nrows(pi);	
	  } else if (isVector(pi)) {
	    q->ncol[i] = 1;
	    q->nrow[i] = LENGTH(pi);
	  } else  {
            PERR("list element(s) neither vector nor matrix");
	  }
	}
      }
      break;
      default : PERR("unmatched internal type of parameter");
  }
}


void includeparam(void **qq,  SEXPTYPE type, int len, SEXP p, 
		  char *param_name) {
    includeparam(qq,  type, len, p, 0, param_name);
}


void CMbuild(SEXP model, int level, cov_model **Cov, cov_model *previous) {
    int covnr, i, j, nkappas,  methNr, 
    elt=0,
    len = length(model);
  char format[30], name[MAXCHAR], param_name[PARAMMAXCHAR], ERR_LOC[nErrorLoc],
      methname[METHODMAXCHAR], msg[200];
  SEXP m, p,
    names = getAttrib(model, R_NamesSymbol);
  cov_fct *C; 
  cov_model *cov;
  bool subleft[MAXSUB]; //, True=true, False=false;


  // printf("hier\n");

  if (TYPEOF(model) != VECSXP) { // not list
    ERR("list expected")
  }
  PROTECT(m = VECTOR_ELT(model, elt++));
  if (!isString(m)) 
    ERR("first element must be the name of a covariance model");
  if (LENGTH(m) != 1) ERR("model name must be a single word");
  strcopyN(name, (char*) CHAR(STRING_ELT(m, 0)), MAXCHAR);

  // printf("hier %s\n", name);

  sprintf(format, "%%s\n%%%ds%%s... ", -level);
  sprintf(ERR_LOC, format, ERROR_LOC, "", name);
  strcpy(ERROR_LOC, ERR_LOC);
  covnr = getmodelnr(name);

  if (covnr == NATSC && NS) ERR("natsc model and RFparameters(PracticalRange=TRUE) may not be given at the same time");
  if (covnr == -2) { ERR("multiple matching of covariance name") }
  else if (covnr == -1) { 
    ERR("unknown covariance model. type `PrintModelList(TRUE, FALSE)' to get the list of models");
  }

  *Cov = NULL;
  addModel(Cov, covnr);
  cov = *Cov;
  cov->nsub = 0;
  C=CovList + cov->nr;

  nkappas = C->kappas;
  for (i=0; i<C->maxsub; i++) subleft[i] = true;
  
  // user defined method
#define NUGGET_PREFERENCE 1

//  printf(">>> %s\n", C->name);

  methNr = getListEltNr(model, "me");  
  if (methNr >= 0) {
      //     printf("methnr=%d\n", methNr);
    p = VECTOR_ELT(model, methNr);
    for (i=0; i<Forbidden; i++)
      cov->user[i] = PREF_NONE;
    if (covnr == NUGGET) cov->user[Nugget] = NUGGET_PREFERENCE;
    cov->user[Nothing] = PREF_BEST;
    strcopyN(methname, CHAR(STRING_ELT(p, 0)), METHODMAXCHAR);
    int meth = Match(methname, METHODNAMES, Forbidden + 1);
    if (meth >= 0) {
      cov->user[meth] = PREF_BEST;
      if (meth == MaxMpp) {
	cov->user[meth] = cov->user[RandomCoin] = PREF_BEST;
      } else if (meth == ExtremalGauss) {
	for (i=0; i<Nothing; i++)
	  cov->user[i] = PREF_BEST;
      }
    } else {
     sprintf(msg, "%s: unknown method\n", methname);
     ERR(msg)
    }
  } else if (methNr == MULTIPLEMATCHING) {
    ERR("method given twice");
  } else { // no method given copy method from calling model
    if (previous != NULL) {
      for (i=0; i<Forbidden; i++) {
	cov->user[i] = previous->user[i];
      }
      if (covnr == NUGGET) cov->user[Nugget] = NUGGET_PREFERENCE;
    }
  }
  

  for ( ; elt < len; elt++) {
    if (elt == methNr) continue;  
    p = VECTOR_ELT(model, elt);
    strcopyN(param_name, names == R_NilValue ? "" 
	     : CHAR(STRING_ELT(names, elt)), PARAMMAXCHAR);     


    if (isVectorList(p) && isString(VECTOR_ELT(p, 0))) {
      // submodels
      if (strcmp(param_name, "") != 0) { // named submodel
	if ((i = Match(param_name, C->subnames, C->maxsub)) < 0) {
	  i = Match(param_name, STANDARDSUB, MAXSUB);
	}
	if (i < 0) {
	  if (PL > 0) {
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

#define MERR(X) {UERR;sprintf(MSG,"%s\n%s: %s",ERROR_LOC, param_name, X);\
error(MSG);}

            MERR("multiple matching of submodel names")
          } else {
            MERR("unmatched covariance name for submodel")
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
	  PRINTF("j=%d maxsub=%d\n", j, C->maxsub);
	  ERR("too many submodels");
	}
      }

      if (!subleft[i]) {
	  MERR("submodel given twice");
      }
      else subleft[i] = false;
      CMbuild(p, level + 1, cov->sub + i, cov);
      cov->sub[i]->calling = cov;
      (cov->nsub)++;
    } else { // parameter (not list)
 
      // parameter identification
      int l = LENGTH(p);

      if (param_name[0]=='k' && strlen(param_name) == 1) { 
	if (l != C->kappas) ERR("length of vector does not match number of parameters. Do not use short name 'k' in complex situtations either.");
	if (C->maxsub!=0)
	  ERR("short form 'k' of parameter name not allowed for sophistacted models");
	if (cov->p[0] != NULL) ERR("parameter given twice");
	for (j=0; j<l; j++) {
	  cov->ncol[j] = cov->nrow[j] = 1;
 	  includeparam((void**) (cov->p + j), C->kappatype[j], 1, p, j, 
		       param_name);
	}
	break;
      }
      
      if (strcmp(param_name, "") == 0) {
	if (nkappas == 1) i = 0; else ERR("parameters must be named");
      } else {
	if ((i = Match(param_name, C->kappanames, nkappas)) < 0) {
	  i = Match(param_name, STANDARDPARAM, MAXPARAM);
	}     
	if (i < 0) {
	  if (PL > 0) {
	    int s;
            PRINTF("allowed parmeter names are: ");
	    if (C->kappas == 0) PRINTF("(none)"); 
	    else for (s=0; s<C->maxsub; s++) {
	      if (strcmp("", C->kappanames[s]) != 0) 
		  PRINTF("%s, ", C->kappanames[s]);
	       PRINTF("%s, ", STANDARDPARAM[s]);
	    }
	    PRINTF("\n");
	  }

	  if (i == MULTIPLEMATCHING) PERR("multiple matching of parameter");
	  
	  if (i == NOMATCHING){
	      PERR("unknown parameter");      
	  }
	}
      }
      if (cov->p[i] != NULL) ERR("parameter given twice");

       //     printf("include %d %d %s\n", i, C->kappatype[i], param_name);

      includeparam((void**) (cov->p + i), C->kappatype[i], l, p, param_name);


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
	  cov->nrow[i] = LENGTH(p);
	} else  {
          PERR("neither vector nor matrix");
	}
      }
    }
  }

  for (j=0; j < C->maxsub; j++) if (subleft[j]) break;
  if (j < C->maxsub) {
    if (j < C->minsub) { 
      // if (PL > 3)  PrintModelInfo(cov);
      ERR("not enough submodels");
    }
    if (PL > 1) {
      for ( ; j < C->maxsub; j++) if (!subleft[j]) break;
      if (j < C->maxsub) {warning("not all submodels are given");}
    }
  }

 

 
//  PrintModelInfo(cov);
  if (NS > 0 && NS != NATSCALE_MLE && C->naturalscale != NULL) {
      // if NS but not from MLE, then add natscale whenever primitive visible
       addModel(Cov, GATTER); 
       addModel(Cov, NATSC); 
  } else {
      if (NS == NATSCALE_MLE && 
	  cov->nr >= DOLLAR && cov->nr <= LASTDOLLAR &&
	  CovList[cov->sub[0]->sub[0]->nr].naturalscale != NULL) {
	  bool natsc;
	  double *pp;

//      printf("---> %s\n", CovList[cov->sub[0]->sub[0]->nr].name);
	  //    PrintModelInfo(cov);

	  pp = cov->p[DSCALE];
	  natsc = pp != NULL &&  (ISNA(pp[0]) || ISNAN(pp[0])); 

	  if (!natsc && pp == NULL) {
	      int i,
		  total = cov->ncol[DANISO] * cov->nrow[DANISO];	      
	      natsc = true;
	      pp = cov->p[DANISO];
	      for (i=0; i<total; i++) {
		  if (!ISNA(pp[i]) && !ISNAN(pp[i])) {
		      natsc = false;
		      break;
		  }
	      }
	  }
	  
	  if (natsc) {
	      cov_model **nextaftergatter = cov->sub[0]->sub;
	      addModel(nextaftergatter, GATTER); 
	      addModel(nextaftergatter, NATSC); 
	  }
//	  PrintModelInfo(cov); assert(false);
     }
  }
	  
 
  addModel(Cov, GATTER); 

  UNPROTECT(1);
}

void DeleteGatter(cov_model **Cov) {
  cov_model *cov = *Cov;

  if (cov->nr >= GATTER && cov->nr <= LASTGATTER && cov->tsdim <= DEL_COV) {
      removeOnly(Cov);
      DeleteGatter(Cov);
  } else {
    int i;
    for (i=0; i<MAXSUB; i++) 
      if (cov->sub[i] != NULL) DeleteGatter(cov->sub + i);
  }
}
 
void checkerror(int err){
    char EM[300+2*MAXERRORSTRING], EM2[300+2*MAXERRORSTRING];
   errorMSG(err, EM); 
   sprintf(EM2, "%s\n    ...%s", ERROR_LOC, EM);
   error(EM2);
}



void CheckModelInternal(SEXP model, int tsdim, int xdim,
			bool stationary, cov_model **Cov,
			int MaxDim) {
  int err,
    PL = GLOBAL.general.printlevel;
  
  if (*Cov != NULL) {
    COV_DELETE(Cov);
  }

  if (currentNrCov==-1) InitModelList();
  if (tsdim < 1 || tsdim > MaxDim || xdim > tsdim) {
    checkerror(ERRORDIM);
  }
  
  strcpy(ERROR_LOC, "Checking model at ");
  CMbuild(model, 0, Cov, NULL);
 
  if (PL>6) PRINTF("detected model structure is:\n%s\n", ERROR_LOC);
  strcpy(ERROR_LOC, "");
  if (PL > 6) {
    PRINTF("\nCheckModelInternal, before checking:");
    PrintModelInfo(*Cov);  // OK
  }

  cov_model *cov= *Cov;
  if (cov != NULL) {
    assert(cov->nr >= GATTER && cov->nr<=LASTGATTER);
    cov->calling=NULL;
    cov->tsdim = tsdim;
    cov->xdim = xdim;
    cov->statIn = stationary ? STATIONARY : PREVMODELS;
    if ((err = CovList[cov->nr].check(cov)) != NOERROR) {
      PRINTF("err =%d\n", err);
      XERR(err);
    }
 
    DeleteGatter(Cov);

    cov = *Cov;
  } // else 
 
  if (PL > 6) {
    PRINTF("CheckModelInternal, after checking:");
    PrintModelInfo(*Cov); // OK
  }
}

void CheckModel(SEXP model, SEXP tsdim, SEXP xdim, SEXP stationary, 
		cov_model **Cov, int MaxDim) {
 
//  assert(INTEGER(tsdim)[0] == 3);
  CheckModelInternal(model, INTEGER(tsdim)[0], INTEGER(xdim)[0], 
		     LOGICAL(stationary)[0], Cov, MaxDim);
}

void strround(double x, char *s){
  if (x==RF_INF)  sprintf(s, "Inf"); else
    if (x==RF_NEGINF)  sprintf(s, "-Inf"); else
      if (x == floor(x + 0.5)) sprintf(s, "%d", (int) x);
      else sprintf(s, "%f", x);
}

void addmsg(double value, const char *sign, double y, char *msg) {
  char str1[30], str2[30];
  strround(value, str1);
  strround(y, str2);
  sprintf(msg, "%s %s %s", str1, sign, str2);
}


void addone(char *str, double x) {
  char str2[30];
  strround(x, str2);
  sprintf(str, "%s, %s", str, str2);
}
void addpair(char *str, double x, double y) {
    if (x == floor(x + 0.5) && y==floor(y + 0.5))
	sprintf(str, "%s, (%d,%d)", str, (int) x, int(y));
    else 
	sprintf(str, "%s, (%f,%f)", str, x, y);
}

void range_default(range_arraytype *ra) {
  ra->n = 1;
  range_type *range = ra->ranges;
  range->maxdim = INFDIM;
  range->finiterange = false;
}


int check_within_range(cov_model *cov, bool NAOK) {
  // sets also maxdim and finiterange in cov !

  cov_fct *C = CovList + cov->nr;
  int n, len,
      err = NOERROR, k = 0, i = 0,   // compiler dummy values
      kappas = C->kappas;
  range_arraytype ra;
  SEXPTYPE *type = C->kappatype;
  rangefct getrange = C->range;
  char msg[255]="";
  double min, max,
    value= RF_NAN;

  if (GLOBAL.general.skipchecks) return NOERROR;

  range_default(&ra);  // default : 1 range; maxdim = inf

  //  printf("%s %ld %ld\n", C->name, C->cov, C->range);

  getrange(cov, &ra); 
  cov->finiterange = ra.ranges[0].finiterange;

  if (kappas == 0) {
    if ((cov->maxdim = ra.ranges[0].maxdim) < cov->tsdim) {
      if (PL>4) PRINTF("max dim is %d %d \n", ra.ranges[0].maxdim, cov->tsdim);
      sprintf(msg, "Max. dimension is %d. Got %d.\n",
	      ra.ranges[0].maxdim, cov->tsdim);
      err = ERRORCOVNOTALLOWED;
      goto ErrorHandling;
    }
  } else {
    char Msg[255];
    range_type *range = ra.ranges;
    for(n=0; n<ra.n; n++, range++) {
	
//	printf("ra %d %f %f\n", ra.n,  range->min[0],  range->max[0]);

       if ((cov->maxdim = range->maxdim) < cov->tsdim) {
	if (ra.n == 1) {
	  sprintf(msg, "max dim is %d. Got %d.", range->maxdim, cov->tsdim);
	  err = ERRORFAILED;
	  break;
	} else {
	  sprintf(msg, "unspecified");
	  err = ERRORFAILED;
	  continue;
	}
      }

      
      err = NOERROR;
      for (i=0; i<kappas; i++) {
	if (type[i] >= LISTOF) {
	  // to simplify things -- otherwise in simu.cc
	  // code must also be changed.
	  assert(range->min[i] == RF_NEGINF && range->max[i]==RF_INF);
	  continue;  
	}
	// full range !!
	
	len=cov->ncol[i] * cov->nrow[i];
	min = range->min[i];
	max = range->max[i];

//	printf("%s ra %d i=%d n=%d %f %f \n", 
//	       C->name, ra.n, i, n,  range->min[i],  range->max[i]);
//	min = -RF_INF;
	
	for (k=0; k<len; k++) {
	  if (type[i] == REALSXP) value = cov->p[i][k];
	  else if (type[i] == INTSXP)
	      value = ((int *)cov->p[i])[k] == NA_INTEGER 
		  ? NA_REAL : (double) ((int *)cov->p[i])[k];
	  else error("parameter value of unknown SXP type");
	  if (ISNA(value) || ISNAN(value)) {
	    if (NAOK) {
	      continue;
	    } else {
	      sprintf(Msg, "%f : finiteness", value);
	      err = ERRORUNSPECIFIED;
	      break;
	    }
	  }

//	  printf("range %d %d %d op=%d %f %f %f\n", 
//		 kappas, i, k, range->openmin[i], value, min, max);

	  if (range->openmin[i] && value <= min) { 
	    addmsg(value, ">", min, Msg);
	    break;
	  } else if (value < min) {
	    addmsg(value, ">=", min, Msg); 
	    break;
	  }
	  if (range->openmax[i] && value >= max) { 
	    addmsg(value, "<", max, Msg); 
	    break;
	  } else if (value > max) { 
	    addmsg(value, "<=", max, Msg);
	    break;
	  }
	}
	
	if (k < len) {
	  err = ERRORUNSPECIFIED;
	  break;
	}
      }
      if (err == NOERROR) break;
    } // ra.n
    
    if (err != NOERROR) {
      PRINTF("--- %d %s %d %d %d %d\n", n, C->name, i, n, ra.n, err);
      if (err != ERRORFAILED)
	sprintf(msg, 
		"%s[%d] = %s does not hold (dim=%d).",
		C->kappanames[i], k+1, Msg, cov->tsdim);
      err = ERRORFAILED;
      goto ErrorHandling;
    }
  }
  cov->normalmix = cov->maxdim == INFDIM;
  
 ErrorHandling:
  if (err!=NOERROR) ERR(msg);
  return err;
}


void get_internal_ranges(cov_model *cov, cov_model *min, cov_model *max, 
	       cov_model *openmin, cov_model *openmax,
	       bool practical) {
  cov_fct *C = CovList + cov->nr;
  int n,
      err = NOERROR, k = 0, i = 0,  // compiler dummy values
      kappas = C->kappas ;
  range_arraytype ra;
  rangefct getrange = C->range;
  SEXPTYPE *type = C->kappatype;

  if (kappas > 0) {
   
    range_default(&ra);
    getrange(cov, &ra); 
    cov->finiterange = ra.ranges[0].finiterange;
    
    range_type *range = ra.ranges;
    for(n=0; n<ra.n; n++, range++) {
      if (range->maxdim < cov->tsdim) {
	assert(ra.n > 1);
	continue;
      }
      
      err = NOERROR;
      for (i=0; i<kappas; i++) {
	int len=cov->ncol[i] * cov->nrow[i];
	double dmin, dmax, dopenmin, dopenmax,
	  value = RF_NAN;
	  if (practical) {
	    dmin = range->pmin[i];
	    dmax = range->pmax[i];
	    dopenmin = dopenmax = 0.0; 

//	    printf("range %s %f %f\n", CovList[cov->nr].name, dmin, dmax);


	  } else {
	    dmin = range->min[i];
	    dmax = range->max[i];
	    dopenmin =  1.0 * (double) range->openmin[i];
	    dopenmax =  1.0 * (double) range->openmax[i];
	  }
	  
	for (k=0; k<len; k++) {
	  if (type[i] == REALSXP) {
	    value = cov->p[i][k];
	    min->p[i][k] = dmin;
	    max->p[i][k] = dmax;
	    openmin->p[i][k] = dopenmin;
	    openmax->p[i][k] = dopenmax;
	  }
	  else if (type[i] == INTSXP) {
	    value = ((int *)cov->p[i])[k] == NA_INTEGER 
	      ? NA_REAL : (double) ((int *)cov->p[i])[k];
	    ((int *) min->p[i])[k] = (int) dmin;
	    ((int *) max->p[i])[k] = (int) dmax;
	    ((int *) openmin->p[i])[k] = (int) dopenmin;
	    ((int *) openmax->p[i])[k] = (int) dopenmax;
	  } else if (type[i] == LISTOF + REALSXP) {
	    listoftype 
		*p = (listoftype*) min->p[i];
	    int j, 
		lenj = p->nrow[k] * p->ncol[k];
	    double
		*pmin = ((listoftype*) max->p[i])->p[k],
		*pmax = ((listoftype*) max->p[i])->p[k],
		*omin = ((listoftype*) openmin->p[i])->p[k],
		*omax = ((listoftype*) openmax->p[i])->p[k];
	
	    for (j=0; j<lenj; j++) {
	      pmin[j] = dmin;
	      pmax[j] = dmax;
	      omin[j] = dopenmin;
	      omax[j] = dopenmax;
	    }
	    
	    value = RF_NAN;
            // error cannot appear as range is (-infty, +infty)
	  } else  error("parameter value of unknown SXP type");
	  if (ISNA(value) || ISNAN(value)) {
	      continue;
	  }
	  dmin = range->min[i];
	  dmax = range->max[i];
	  if (value < dmin
	      || value > dmax
	      || (range->openmin[i] && value == dmin)
	      || (range->openmax[i] && value == dmax)
	    ) break;
	}
	if (k < len) {
	  err = ERRORDUMMY;
	  break;
	}
      }
      if (err == NOERROR) break;
    } // ra.n
    assert(err == NOERROR);
  } // if (kappa > 0)

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      get_internal_ranges(cov->sub[i], min->sub[i], max->sub[i],
		openmin->sub[i], openmax->sub[i], practical);
    }
  }
}

void get_ranges(cov_model *cov, cov_model **min, cov_model **max, 
	       cov_model **openmin, cov_model **openmax,
	       bool practical) {
  // returns a reliable value only in the very first
  // entry of each parameter vector/matrix

  covcpy(min, cov, false, true);
  covcpy(max, cov, false, true);
  covcpy(openmin, cov, false, true);
  covcpy(openmax, cov, false, true);
  get_internal_ranges(cov, *min, *max, *openmin, *openmax, practical);
}


void Covariance(double *x, double *y, bool y_not_given, int lx, 
		cov_model *cov, int idx, double *value) {
  int  vsq = cov->vdim * cov->vdim,
    xdim = cov->xdim;

  if (cov->pref[Nothing] == 0  || 
      (cov->statIn != STATIONARY && (cov->statIn != COVARIANCE))) {
	ERR("Covariance cannot be calculated");
  }
  
//  PrintModelInfo(cov);

  // CovMatrixIndexActive = idx >= 0;
  CovMatrixRow = idx;
  if (y_not_given) {
    for (CovMatrixCol=0; CovMatrixCol<lx; CovMatrixCol++, x += xdim) { 
      CovList[cov->nr].cov(x, cov, value + CovMatrixCol * vsq);
    }
  } else {
    for (CovMatrixCol=0; CovMatrixCol<lx; CovMatrixCol++, x += xdim, y+=xdim)
      CovList[cov->nr].nonstat_cov(x, y, cov, value + CovMatrixCol * vsq);
  }
  //CovMatrixIndexActive = false;
  //PrintModelInfo(cov);
}



void CovMulti(double *x, double *y, int lx, cov_model *cov, double *value) {

  int i, ii, v,j,
      vdim = cov->vdim,
      vdimsq = vdim * vdim,
      vdimlx = vdim * lx,
      vdim2lx = vdimsq * lx,
//    vdlxP1 = vdimlx + 1,
      xdim = cov->xdim,
      err = NOERROR;

  double **Val=NULL, *cross=NULL, *z=NULL,
      *d=x;
  
  if ((Val = (double**) malloc(sizeof(double*) * vdimsq)) == NULL){
	err = ERRORMEMORYALLOCATION;
	goto ErrorHandling;
  }
  if ((cross= (double*) malloc(sizeof(double) * vdimsq)) == NULL){
      err = ERRORMEMORYALLOCATION;
      goto ErrorHandling;
  }
  if ((z= (double*) malloc(sizeof(double) * xdim)) == NULL){
      err = ERRORMEMORYALLOCATION;
      goto ErrorHandling;
  }
 

  // for (i=0;i<1000; i++) {
//      printf("%d \n", i);
//      value[i] = 1.0;
//  }
  
  for (v = i = 0; i<vdim2lx ; i+=vdimlx) {
    for (ii=i, j=0; j<vdim; ii+=lx, v++, j++) {
//	printf("ii = %d (%d %d)\n", ii, lx, vdim);
      Val[v] = value + ii;
    }
  }
  
  CovMatrixTotal = 0;
  if (cov->statIn == STATIONARY) {
    covfct cf = CovList[cov->nr].cov;

    CovMatrixCol = - 1;
    for (CovMatrixRow=0; CovMatrixRow<lx; CovMatrixRow++, d+=xdim) {
      cf(d, cov, cross);
      for (v = 0; v<vdimsq; v++) {
//	  printf("%d %d %d %f\n", 
//		 v, CovMatrixRow, Val[v]+CovMatrixRow-value, cross[v]);
	Val[v][CovMatrixRow] = cross[v];
      }
    }
  } else if (cov->statIn == COVARIANCE){
    nonstat_covfct cf = CovList[cov->nr].nonstat_cov;
    double *dy=y;
    for (CovMatrixRow=0; CovMatrixRow<lx; CovMatrixRow++, d+=xdim, dy+=xdim) {
      cf(d, dy, cov, cross);
      for (v = 0; v<vdimsq; v++) {
	Val[v][CovMatrixRow] = cross[v];
      }
    }
  } else {
      assert(false);
  }

ErrorHandling:
  CovMatrixCol = CovMatrixRow = RF_NAN;
  if (Val!=NULL) free(Val);
  if (cross!=NULL) free(cross);
  if (z!=NULL) free(z);
  if (err != NOERROR) XERR(err); 
}


void CovIntern(double *x, double *y, bool y_not_given, int lx, double *value) {
    cov_model *cov = STORED_MODEL[MODEL_INTERN];
    if (cov->vdim > 1) {
 	CovMulti(x, y, lx, cov, value);
    } else {
	Covariance(x, y, y_not_given, lx, cov, -1, value);
    }
}

void Cov(double *x, double *y, int *y_not_given, int *lx, double *value) {
    cov_model *cov = STORED_MODEL[MODEL_USER];
    if (cov->vdim > 1) {
 	CovMulti(x, y, *lx, cov, value);
    } else {
      Covariance(x, y,  *y_not_given, *lx, cov, -1, value);
    }   
}


void CovMatrixMulti(double *x, bool dist, int size, 
		    cov_model *cov, double *value) {
    // geblockt nach variablen -- somit koennen die univariaten
    // covariance matrizen problemlos ausgelesen werden.
//  cov->vdim = 4;

    // hier kann noch etwas beschleunigt werden, dadurch, dasss
    // die Werte gespiegelt und nciht neu berechnet werden.

    int i, ii, v, k,j,
    vdim = cov->vdim,
    vdimsq = vdim * vdim,
    vdimsize = vdim * size,
    sizesqvdim = vdimsize * size,
    vdimsize2 = vdimsize * vdimsize,
//    vdsizeP1 = vdimsize + 1,
	err = NOERROR,
    xdim = cov->xdim,
    sizeM1 = size-1;

    double **Val=NULL, *cross=NULL, *z=NULL, *d;
 
    if ((Val = (double**) malloc(sizeof(double*) * vdimsq)) == NULL) {
	err = ERRORMEMORYALLOCATION;
	goto ErrorHandling;
    }
    if ((cross= (double*) malloc(sizeof(double) * vdimsq)) == NULL){
	err = ERRORMEMORYALLOCATION;
	goto ErrorHandling;
    }
    if ((z= (double*) malloc(sizeof(double) * xdim)) == NULL){
	err = ERRORMEMORYALLOCATION;
        goto ErrorHandling;
    }
   

// printf("%d %d %d dist=%d\n", vdim, vdimsq,xdim, (int) dist);

  // PrintModelInfo(cov);
  // assert(false);
//  print("\n");

//  for (i=0; i<vdimsq; i++) { cross[i]=3.1415; Val[i] = NULL; }
//  for (i=0; i<xdim; i++) z[i] = 3.1415;

  for (v = i = 0; i<vdimsize2 ; i+=sizesqvdim) {
    for (ii=i, j=0; j<vdim; ii+=size, v++, j++) {
	Val[v] = value + ii;
    }
  }
  
  CovMatrixTotal = 0;
  if (dist) {
    covfct cf = CovList[cov->nr].cov;

    // lower triangle
    for (CovMatrixCol=0; CovMatrixCol<size; CovMatrixCol++) {
      for (CovMatrixRow=CovMatrixCol+1; CovMatrixRow<size; CovMatrixRow++) {

	  d = x + 
	      (CovMatrixCol * sizeM1 - (CovMatrixCol * (CovMatrixCol + 1)) / 2
		     + CovMatrixRow -1) * xdim;

//	  printf("A: %d %d size=%d %d %d %d %ld %f\n", CovMatrixCol, CovMatrixRow, size, xdim, vdim, CovMatrixCol * sizeM1 - (CovMatrixCol * (CovMatrixCol + 1)) / 2
//		 + CovMatrixRow, d, *d);
	  
//	  PrintModelInfo(cov);


	  for (k=0; k<xdim; k++) { z[k] = - d[k]; }
	  cf(z, cov, cross);
	  //    assert(false);

	for (v = 0; v<vdimsq; v++) {
	  Val[v][CovMatrixCol * vdimsize + CovMatrixRow] = cross[v];

//	  print("%f %f %f   %d %e %d %d %d %d %d %ld %d\n", 
//		z[0], z[1], z[2],
//		v,  cross[v], CovMatrixCol * vdimsize + CovMatrixRow,
//		CovMatrixCol, vdimsize, CovMatrixRow,
//		&(Val[v][CovMatrixCol * vdimsize + CovMatrixRow]), 
//		&(Val[v][CovMatrixCol * vdimsize + CovMatrixRow]) -value);
//	  assert(fabs(cross[v]) > 0.000001);
	}
      }
    }

    //   assert(false);

    //   assert(false);
//A: 0 1 size=376 1 2 374 1085485096
//A: 374 375 size=376 1 2 70499 1086049088
//
   // upper triangle
    for (CovMatrixCol=0; CovMatrixCol<size; CovMatrixCol++) {
      for (CovMatrixRow=0; CovMatrixRow<CovMatrixCol; CovMatrixRow++) {
//	  printf("%d %d size=%d %d %d %d %ld\n", CovMatrixCol, CovMatrixRow, size, xdim, vdim, CovMatrixRow * sizeM1 - (CovMatrixRow * (CovMatrixRow + 1)) / 2
//		 + CovMatrixCol, d);
	  d = x +(CovMatrixRow * sizeM1 - (CovMatrixRow * (CovMatrixRow + 1)) / 2
		     + CovMatrixCol - 1) * xdim;
	cf(d, cov, cross);
	for (v = 0; v<vdimsq; v++) {
	  Val[v][CovMatrixCol * vdimsize + CovMatrixRow] = cross[v];
	}
      }
    }
//assert(false);
   // diag
    for (CovMatrixCol=0; CovMatrixCol<size; CovMatrixCol++) {
      CovMatrixRow = CovMatrixCol;
      cf(ZERO, cov, cross); // can be different in each step by 
      //                   covariance models depending on CovMatrixCol, -Row
      for (v = 0; v<vdimsq; v++) {
	  Val[v][CovMatrixCol * vdimsize + CovMatrixRow] = cross[v];
      }
    }
  } else {

///      assert(false);

    nonstat_covfct cf = CovList[cov->nr].nonstat_cov;
    double *X, *Y;
    for (Y=x, CovMatrixCol=0; CovMatrixCol<size; CovMatrixCol++, Y+=xdim) {
      for (X=x, CovMatrixRow=0; CovMatrixRow<size; CovMatrixRow++, X+=xdim) {
        cf(Y, X, cov, cross);
	for (v = 0; v<vdimsq; v++) {
	    Val[v][CovMatrixCol * vdimsize + CovMatrixRow] = cross[v];
//	    printf("%f %f %f ; %f %f %f ; %d %d %f",
//		   X[0], X[1], X[2], Y[0], Y[1], Y[2], v,
//		   CovMatrixCol * vdimsize + CovMatrixRow, cross[v]);
	}
      }
    }
  }

  if (cov->statIn==VARIOGRAM) {
    double first=value[0];
    for (i=0; i<vdimsize2; value[i++] -= first);
  }

ErrorHandling:
  CovMatrixCol = CovMatrixRow = RF_NAN;
  if (Val!=NULL) free(Val);
  if (z!=NULL) free(z);
  if (cross!=NULL) free(cross);
  if (err != NOERROR) XERR(err); 
}



void CovarianceMatrix(double *x, bool dist, int size, cov_model *cov, 
		      double *value) {
   int i, endfor,  ve,
	sizeP1 = size + 1,
	xdim = cov->xdim; 
  long ho;

  // double *origx = x;

 
  if (cov->pref[Nothing] == 0) {
    ERR("Covariance cannot be calculated (forbidden).");
  }


  if (cov->vdim > 1) {
      // FEHLT: dort die index-Menge mitzunehmen!!
      CovMatrixMulti(x, dist, size, cov, value);
    return;
  }
  
  CovMatrixTotal = 0;
  if (dist) {    
    covfct cf = CovList[cov->nr].cov;
    for (CovMatrixCol=0, i=0; CovMatrixCol<size; i+=sizeP1, CovMatrixCol++) {	
      CovMatrixRow = CovMatrixCol;
      cf(ZERO, cov, value + i);
      endfor = i + (size - CovMatrixCol);
      for (ve = i + 1, ho = i + size, CovMatrixRow++; ve < endfor; 
	   ve++, ho += size, x += xdim, CovMatrixRow++, CovMatrixTotal++) {

//	    print("%d  %d %d xdim=%d %ld size=%d %ld %ld %f\n", i, ve, ho,
//		 xdim, origx, size,
//		 origx + size * (size-1) * xdim / 2,
//		 x, x[0]);
	
	cf(x, cov, value + ve);

//	printf("(%d, %d) ho=%d %d v=%f %f %d %s \n",  CovMatrixRow,CovMatrixCol,
//	       ho, ve, value[ve], x[0], cov->sub[0]->nr, 
//	       CovList[cov->sub[0]->nr].name);

	value[ho] = value[ve];
      }
    }
  } else {
    nonstat_covfct cf = CovList[cov->nr].nonstat_cov;
    double *X, *Y;
    Y = x;
    for (CovMatrixCol=0, i=0; CovMatrixCol<size;
	 i+=sizeP1, CovMatrixCol++, Y += xdim) {
      CovMatrixRow = CovMatrixCol;
      cf(Y, Y, cov, value + i);
      endfor = i + (size - CovMatrixCol);
      for (X = Y + xdim, ve = i + 1, ho = i + size, CovMatrixRow++; ve < endfor; 
	   ve++, ho += size, X += xdim, CovMatrixRow++, CovMatrixTotal++){
	cf(X, Y, cov, value + ve);
	value[ho] = value[ve];
      }
    }
  }

  if (cov->statIn==VARIOGRAM) {
    int l, sizesq = size * size;
    double first=value[0];
    for (l=0; l<sizesq; value[l++] -= first);
  }
 

  CovMatrixCol = CovMatrixRow = RF_NAN;
}

int check_within_range(cov_model *cov) {
  return check_within_range(cov, false);
}

void CovMatrix(double *x, int *dist, int *size, double *value) {
    CovarianceMatrix(x, (bool) *dist, *size, STORED_MODEL[MODEL_USER], value);
}

void CovMatrixIntern(double *x, int *dist, int *size, double *value) {
      CovarianceMatrix(x, (bool) *dist, *size, STORED_MODEL[MODEL_INTERN], value);
}



void CalculateVariogram(double *x, int lx, cov_model *cov, double *value) {
  int i, m, k,
    vsq = cov->vdim * cov->vdim,
    xdim = cov->xdim;
  double *dummy = (double *) malloc(sizeof(double) *vsq),
    *C0 = (double *) malloc(sizeof(double) *vsq);
  covfct cf = CovList[cov->nr].cov;
  

  cf(ZERO, cov, C0);

  if (cov->statIn == STATIONARY || cov->statIn == VARIOGRAM) {
    for (k=i=0; i<lx; i++, x += xdim) { 
      cf(x, cov, dummy);
      for (m=0; m<vsq; m++)
	value[k++] = C0[m] - dummy[m];
    }
  } else if (cov->statIn == COVARIANCE) {
    error("cov->stationary == COVARIANCE not programmed yet");
  } else {
    error("unvalid stationarity assumption");
  }
  free(C0);
  free(dummy);
}

void Variogram(double *x, double *y, int *y_not_given, int *lx, double *value) {
  if (!y_not_given) error("y may not be given for a variogram.");
  CalculateVariogram(x, *lx, STORED_MODEL[MODEL_USER], value);
}

void VariogramIntern(double *x, double *y, int *lx, double *value) {
  assert((SEXPREC*) y==R_NilValue);
// PrintModelInfo(STORED_MODEL[MODEL_INTERN]);
//  assert(false);

  CalculateVariogram(x, *lx, STORED_MODEL[MODEL_INTERN], value);
}

void VariogramMLE(double *x, int *lx, double *value) {
  // assert((SEXPREC*) y==R_NilValue);
// PrintModelInfo(STORED_MODEL[MODEL_INTERN]);
//  assert(false);

  CalculateVariogram(x, *lx, STORED_MODEL[MODEL_MLE], value);
}


void InitSimulateRF(double *x, double *T, 
		    int *spatialdim, /* spatial dim only ! */
		    int *lx, 
		    int *grid,
		    int *Time, 
		    int *distr, /* still unused */
		    int *keyNr, 
		    int *expected_number_simu,
		    int *err) {
  // keyNr : label where intermediate results are to be stored;
  //         can be chosen freely between 0 and MAXKEYS-1
  key_type *key = KEY + *keyNr;
  strcpy(ERROR_LOC, "");
  strcpy(PREF_FAILURE, "");
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *err=ERRORREGISTER; 
    goto ErrorHandling;
  }

  if ((SEXPREC*) x==R_NilValue) {
    *err=ERRORCOORDINATES; 
    goto ErrorHandling;
  }
  *err = internal_InitSimulateRF(x, T, *spatialdim, 
				   *lx,  (bool) *grid, (bool) *Time,  
				   *distr, key, &GLOBAL, 
				   *expected_number_simu,
				   STORED_MODEL + MODEL_SIMU);

  if (*err == NOERROR) return;
 
 ErrorHandling:
  char EM[500], EM2[500 + 100 * Nothing];
  key->simu.active=false;
  errorMSG(*err, EM);

  sprintf(EM2, "%s%s\n", PREF_FAILURE, EM);  

  ERR(EM2);
  // does not return !
} 

int internal_InitSimulateRF(double *x, double *T, 
			    int spatialdim, /* spatial dim only ! */
			    int lx,   bool grid, bool Time,
			    int distr, /* still unused */
			    key_type *key, globalparam *gp,
			    int expected_number_simu,
			    cov_model **COV)
{ 
  // NOTE: if grid and Time then the Time component is treated as if
  //       it was an additional space component
  //       On the other hand, SPACEISOTROPIC covariance function will
  //       always assume that the last calling component is the time component!


  char errloc_save[nErrorLoc];
  int err, dim;
  location_type *loc;
  
  strcpy(errloc_save, ERROR_LOC);
  
  KEY_DELETE_NOTTREND(key);
  KEY_NULLNOTTREND(key);
  key->cov = *COV;
  *COV = NULL;
  key->simu.distribution=distr; 
  err = loc_set(x, T, spatialdim, lx, Time, grid, &(key->loc), gp); 
  if (err) goto ErrorHandling;
  memcpy(&(key->gp), gp, sizeof(globalparam));
  loc = &(key->loc);
  dim = loc->timespacedim;
  if (err) goto ErrorHandling; // muss als 3. stehen

  // setmethod
  method_type *meth;
  meth = key->meth = (method_type*) malloc(sizeof(method_type));
  METHOD_NULL(meth);
  meth->gp = &(key->gp);
  meth->gpdo = &(key->gpdo);
  meth->loc = &(key->loc);
  meth->simu = &(key->simu);
  meth->cov = key->cov; // oder besser kopieren ??
  meth->xdimout = dim;
 
  meth->cvar = 1.0;
  meth->type = TypeDiag;
  meth->space = meth->sptime = NULL;


  err = initstandard(meth);
  if (err) goto ErrorHandling;
 
  key->simu.active = true;    
  strcpy(ERROR_LOC, errloc_save);
  return NOERROR;
  
  ErrorHandling:
  return err;
}



void AddTrend(int *keyNr, int *n, double *res, int *err) {
  int i;
  key_type *key;
  trend_type *trend;
  location_type *loc;

  *err = NOERROR;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *err=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  key = &(KEY[*keyNr]);
  if (!key->simu.active) {*err=ERRORNOTINITIALIZED; goto ErrorHandling;}
  trend = &(key->trend);
  loc = &(key->loc);

  switch (trend->TrendModus) {
      case TREND_MEAN :
	if (trend->mean!=0.0) {
	  long int endfor= *n * loc->totalpoints;
	  for(i=0; i<endfor; i++) res[i] += trend->mean;
	}
	break;
      case TREND_LINEAR :
      case TREND_FCT :
      case TREND_PARAM_FCT :
      default: assert(false); // not programmed yet
  }
  return;

 ErrorHandling: 
  if (PL>0) ErrorMessage(Nothing,*err);
 return;
}
 
void DoPoissonRF(int *keyNr, int *pairs, int *n, double *res, int *err)
{
  // not programmed yet
  assert(false);
}

int internal_DoSimulateRF(key_type *key, int nn, res_type *orig_res) {
  // does not assume that orig_res[...] = 0.0, but it is set
  simu_type *simu = &(key->simu);
  long i,
    vdimtot = key->loc.totalpoints * key->cov->vdim;
  res_type 
      *res = orig_res;
  double
      realeach=0.0;
  char format[20],
    back[]="\b\b\b\b\b\b\b\b\b\b\b", 
    prozent[]="%",
    pch = key->gpdo.general.pch;
  int ni, digits, err,
    each=0;
  method_type *meth=key->meth;
  
  memcpy(&(key->gpdo), &GLOBAL, sizeof(globalparam));

  if (!simu->active) {err=ERRORNOTINITIALIZED; goto ErrorHandling;}
  if (nn>1 && pch != '\0') {
    if (pch == '!') {
      digits = (nn<900000000) ? 1 + (int) trunc(log((double) nn) / log(10.0)) : 9;
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
  GetRNGstate();

  // PrintMethodInfo(meth); assert(false);
  
  for (ni=1; ni<=nn; ni++, res += vdimtot) {
    R_CheckUserInterrupt();
    if (simu->stop) {
      err=ERRORNOTINITIALIZED; goto ErrorHandling;
    }
    if (ni % each == 0) {
      if (pch == '!')  
	  PRINTF(format, back, ni / each);
      else if (pch == '%')
	  PRINTF(format, back, (int) (ni / realeach), prozent);
      else PRINTF("%c", pch);
    }
    if (meth->compatible) {
      for (i=0; i<vdimtot; i++) {
	res[i] = 0.0;
      }
    }
       
  
    assert(meth->domethod != NULL);

    meth->domethod(meth, res);
  } // for n
  PutRNGstate();
  if (nn>1 && pch != '\0') {
    if (pch == '!' || pch == '%') PRINTF("%s", back);
    else PRINTF("\n");
  }
  return NOERROR;

 ErrorHandling: 
  PutRNGstate();
  simu->active = false;
  return err;
}


void DoSimulateRF(int *keyNr, int *n, int *pairs, res_type *res, int *err) {
  // does not assume that orig_res[...] = 0.0, but it is set
  int i, internal_n;
  key_type *key=NULL;
  simu_type *simu;


  *err=NOERROR; 
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *err=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  key = &(KEY[*keyNr]);
  internal_n = *n / (1 + (*pairs!=0));
  simu = &(key->simu);

  if (key->gp.general.printlevel > 7) {
//    PrintMethodInfo(key->meth);
  }

  if (!simu->active) {
    *err=ERRORNOTINITIALIZED; goto ErrorHandling;
  }
  *err = internal_DoSimulateRF(key, internal_n, res);
  if (*err != NOERROR) goto ErrorHandling;

  if (*pairs) {
    res_type * res_pair;
    long endfor=key->loc.totalpoints * *n / 2 * key->cov->vdim;
    res_pair = res + endfor;
    for (i=0; i<endfor; i++) res_pair[i] = -res[i];
  }

  //AddTrend(keyNr, n, res, err);
  if (*err) goto ErrorHandling;
 
  if (!(simu->active = key->gp.general.storing)) {
    KEY_DELETE(key);
  }
  return; 
  
 ErrorHandling: 
  if (key!=NULL) {
    simu->active = false;
    KEY_DELETE(key); 
  }
  return;
}

// nicht in userinterfaces, da zugriff auf user_... und unchecked_...

