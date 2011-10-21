
/* 
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 library for simulation of random fields 

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

extern "C" {
#include <R.h>
#include <Rdefines.h>
}
#include <R_ext/Linpack.h>
#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
 
#include <string.h>
#include "RF.h"
#include "primitive.h"
// #include "Covariance.h"


/* 
   operator "#" that has been modified such that explanatory
   variables can be used instead of the spatial coordinates
*/
void iso2iso_MLE(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);  
  else {
    CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
}
void spiso2spiso_MLE(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);  
  else {
    CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
}
void spacetime2iso_MLE(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);  
  else {
    CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
}

void Stat2iso_MLE(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);  
  else {
    CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
}
void Nonstat2iso_MLE(double *x, double *y, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);  
  else {
    CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
}
void Stat2spacetime_MLE(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);  
  else {
    CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
}
void Nonstat2spacetime_MLE(double *x, double *y, cov_model *cov, double *v) {
  //assert(false);
  cov_model *next = cov->sub[0];
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);  
  else {
      CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
}
void Stat2Stat_MLE(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);  
  else {
    CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
}
void Nonstat2Stat_MLE(double *x, double *y, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);  
  else {
    CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
}

/*
*/

  // simple transformations as variance, scale, anisotropy matrix, etc.  
void Siso_MLE(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
    vdimSq = cov->vdim * cov->vdim;
  double **p = cov->p,
    var = p[DVAR][0];

  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);//x wird nicht verwendet
  else {
      CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
  for (i=0; i<vdimSq; i++) v[i] *= var; 
}

void Sstat_MLE(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  double 
      **p = cov->p, 
      var = p[DVAR][0];
  int i,
      vdimSq = cov->vdim * cov->vdim;
 
  if (cov->MLE == NULL) CovList[next->nr].cov(x, next, v);//x wird nicht verwendet
  else {
      CovList[next->nr].cov(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }
 
  for (i=0; i<vdimSq; i++) v[i] *= var; 
}

void Snonstat_MLE(double *x, double *y, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  double 
      **p = cov->p, 
      var = p[DVAR][0];
  int i,
      vdimSq = cov->vdim * cov->vdim;
  
  if (cov->MLE == NULL) CovList[next->nr].nonstat_cov(x, y, next, v);  
  else {
      double *xx;
      xx = cov->MLE + CovMatrixTotal * next->xdim * 2;
      CovList[next->nr].nonstat_cov(xx, xx + next->xdim, next, v);
  }

  for (i=0; i<vdimSq; i++) v[i] *= var; 
}


void DS_MLE(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
      vdimSq = cov->vdim * cov->vdim;
  double  
    **p = cov->p,
    spinvscale = (p[DANISO] == NULL) ? 1.0 / p[DSCALE][0] : p[DANISO][0], 
    varSc = p[DVAR][0] * spinvscale;
  
 
  if (cov->MLE == NULL) CovList[next->nr].D(x, next, v);//x wird nicht verwendet
  else {
      CovList[next->nr].D(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }

  for (i=0; i<vdimSq; i++) v[i] *= varSc; 
}

void DDS_MLE(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
      vdimSq = cov->vdim * cov->vdim;
   double 
    **p = cov->p,
    spinvscale = (p[DANISO] == NULL) ? 1.0 / p[DSCALE][0] : p[DANISO][0], 
    varScSq = p[DVAR][0] * spinvscale * spinvscale;
  
  if (cov->MLE == NULL) CovList[next->nr].D2(x, next, v);//x wird nicht verwendet
  else {
      CovList[next->nr].D2(cov->MLE + CovMatrixTotal * next->xdim, next, v);
  }

  for (i=0; i<vdimSq; i++) v[i] *= varScSq; 
}


// nachfolgende Funktion sowie anyNAscaleup, anyNAdow, manipulate_x
// nur fuer speed in MLE verwendet, wenn die x-Werte direkt den
// knoten vorliegen sollen 
void SetPrevToTriMaxFalse(cov_model *calling) {
  // bislang war TriFalse; jetzt nicht mehr;
  // setze letzes DOLLAR auf  TriMaxFalse;
  if (calling != NULL && calling->anyNAscaleup == TriFalse) {
    cov_model *prev = calling;
    while (prev != NULL &&
	   (prev->nr < DOLLAR || prev->nr > LASTDOLLAR)) {
	//  prev->anyNAscaleup = TriBeyond;
      prev = prev->calling;
    }
    if (prev != NULL) prev->anyNAscaleup = TriMaxFalse;
  }
}


void GetNAPosition(cov_model *cov, int *usern, int *internn, 
		   naptr_type mem, naptr_type address,
		   NAname_type names, sortsofparam *sorts,
		   sortsofparam *internsorts, bool *isnan,
		   int *covzaehler,
		   int printing, int depth) {
  /*
    printing <= 0 nothing
             1 only NA, NaN
	     2 others as "*"
	     3 others as value
s
   mem[0..internn] : addresse wo (neue) NA position im Speicher	   
   address[0..internn] 
       == NULL         : address in mem the value should be put;
       == mem[internn] :
          * internsort=ANISOFROMSCALE : point to scale (using 1/value)
   sorts[0..usern]     : see sortsofparam in RF.h
   internsorts[0..internn] : see sortsofparam in RF.h
   isnan[0..usern] : NA or NaN
  */
  int i, c, r, 
      namenr = 1,     
      *nrow =  cov->nrow,
      *ncol = cov->ncol;
  //  double *lastmem = NULL;
#define SHORTlen 4
  char shortname[SHORTlen+2], shortD[SHORTlen+2];
  bool user;
  cov_fct *C = CovList + cov->nr,
    *CC = C;
  SEXPTYPE *type = C->kappatype;

  // C is pointer to true model; CC points to the first model passed
  // in initNerror.cc by IncludeModel or IncludePrim
  while(strncmp(CC->name, InternalName, strlen(InternalName)) ==0) CC--;
  covzaehler[cov->nr]++;
  strcopyN(shortname, CC->name, SHORTlen);

   if (covzaehler[cov->nr] >= 2) {
      char dummy[SHORTlen];
      strcopyN(dummy, shortname, SHORTlen-1);
      sprintf(shortname, "%s%d", dummy, covzaehler[cov->nr]);
  }
  if (printing>0) PRINTF("%s\n", CC->name); 
  // CC needed below for the kappa.names which are given

//  printf("%d:%s %d %s\n", cov->nr, CC->name, covzaehler[cov->nr], shortname);

  if (cov->manipulating_x > 0) {
      // any scale==NA in the calling function
    cov->anyNAscaleup = cov->anyNAdown = TriTrue;
    SetPrevToTriMaxFalse(cov->calling);
  } else {
    cov->anyNAscaleup = cov->anyNAdown = TriFalse;
    if (cov->calling != NULL) cov->anyNAscaleup = cov->calling->anyNAscaleup;
  }

  depth++;
  for (i=0; i<C->kappas; i++) {
    if (nrow[i] == 0 || ncol[i] == 0) continue;
    if (printing > 0) {
      leer(depth); PRINTF("%s\n", C->kappanames[i]);
    }
    for (r=0; r<nrow[i]; r++) {
       int nv = 0; // anzahl NA in aktuellem parameter
       for (c=0; c<ncol[i]; c++) {
	double v = RF_NAN; // value in aktuellem parameter
	int idx = c * nrow[i] + r;
	if (*internn >= MAX_NA) error("maximum number of NA reached");

	if (type[i] == REALSXP) {
	  v = cov->p[i][idx];
	  mem[*internn] = &(cov->p[i][idx]);
	} else if (type[i] == INTSXP) {
	  v = ((int *) cov->p[i])[idx] == NA_INTEGER 
	      ? NA_REAL : (double) ((int *) cov->p[i])[idx];
	  if (ISNA(v) || ISNAN(v))
	    error("integer variables currently do not allowed for NA"); // !!!
	  // mem[*internn] = (*double) &(((int *) cov->p[i])[idx]);
	} else if (type[i] == LISTOF + REALSXP) {
	    listoftype *q;
	    int j, end;
	    double *p;
	    q=(listoftype*) cov->p[i];
	    p = q->p[r];
	    end = q->nrow[r] * q->ncol[r];
	    for (j=0; j<end; j++)
	      if (ISNA(p[j]) || ISNAN(p[j])) 
		error("no NAs allowed in regression ");
	    v = 0.0; // dummy
	} else {
	  error("unknown SXP type");
	}
	user  = true;

 	isnan[*usern] = ISNAN(v) && !ISNA(v);
	if (ISNA(v) || ISNAN(v)) { // entgegen Arith.h gibt ISNA nur NA an !!
	  cov->anyNAdown = TriTrue;
	  if (printing > 0) {
	    if (nv>1 || (c>0 && printing > 1)) PRINTF(", "); else leer(depth+1);
	  }
	  if (cov->nr >= DOLLAR && cov->nr <= LASTDOLLAR) {
	    // shortD partial name for R level
	    cov_model *next = cov->sub[0];
	    while((next->nr >= GATTER && next->nr <= LASTGATTER)  ||
		  (next->nr >= DOLLAR && next->nr <= LASTDOLLAR) || 
		  (next->nr == NATSC))
	      next = next->sub[0];
	    if (covzaehler[next->nr] == 0) { // next wurde noch nicht
		// untersucht, somit ist covzaehler[next->nr] um 1 niedriger
		// als covzaehler[cov->nr] !
	      strcopyN(shortD, CovList[next->nr].name, SHORTlen);
	    } else {
	      char dummy[SHORTlen];
	      strcopyN(dummy, CovList[next->nr].name, SHORTlen-1);
	      sprintf(shortD, "%s%d", dummy, covzaehler[next->nr]+1);
//	      printf("$$ > %d:%s %d %s\n", next->nr, CovList[next->nr].name, 
//		     covzaehler[next->nr], shortD);
	    }
	    
	    //  assert(DANISO == DSCALE + 1);

	    address[*internn] = NULL;
	    if (i==DVAR) {
	      internsorts[*internn] = sorts[*usern] = 
		  // zur sicherheit beide:
		  (next->nr == MIXEDEFFECT || next->nr == MLEMIXEDEFFECT) ? 
		  MIXEDVAR : next->nr == NUGGET ? NUGGETVAR : VARPARAM; 
	      sprintf(names[*usern], "%s.var", shortD); // for R level only
	    } else  {
	      cov->anyNAscaleup = TriTrue; 
	      SetPrevToTriMaxFalse(cov->calling);
	      if (i==DSCALE) {
		// lastmem = mem[*internn]; // used for all diagonal elements of
		//                          aniso
		sorts[*usern] = internsorts[*internn] = SCALEPARAM;	
		sprintf(names[*usern], "%s.s", shortD);// for R level only
	      } else {
		assert(i == DANISO);
		// no lastmem here !
		
		if (cov->ncol[DSCALE] != 0 &&
		    (ISNA(cov->p[DSCALE][0]) || ISNAN(cov->p[DSCALE][0]))) {
		    assert(false);
		    /*
		  // internal
		  address[*internn] = lastmem; // points to scale which 
		  //                           was just be before and had NA
		  internsorts[*internn] = ANISOFROMSCALE; 
		  sprintf(names[*usern], "%s.aniso", shortD);// R level only
		  // no setting of sorts[users]
		  if (printing > 1) {
		    if (printing > 2) PRINTF("(given by scale)");
		    else PRINTF("*");  
		  }
		  user = false;
		    */
		} else {
		  if (cov->q != NULL && cov->q[0] != 0.0) {
//		    printf("sdfasdf\n");
		    ERR("naturalscaling only allowed for isotropic setting");
		  }
		  if (r==c) { 
		    internsorts[*internn] = sorts[*usern] = DIAGPARAM; 
		    sprintf(names[*usern], "%s.d.%d", shortD, r);// R level only
		  } else {
		    internsorts[*internn] = sorts[*usern] = ANISOPARAM;
		    sprintf(names[*usern], "%s.d.%d.%d", shortD, r, c);// R level
		  }	
		}	
	      }
	    }
	  } else {
	    // standard setting !
	    address[*internn] = NULL; 
	    if (cov->nr==MIXEDEFFECT || cov->nr==MLEMIXEDEFFECT)
		continue;
	    else 
		sorts[*usern] = internsorts[*internn] = ANYPARAM; 
	    char kappashort[SHORTlen];
	    strcopyN(kappashort, CC->kappanames[i], SHORTlen);
	    sprintf(names[*usern], "%s.%s", shortname,  kappashort);
	    if (*usern > 0 && 0==
	        strncmp(names[*usern], names[*usern-1], strlen(names[*usern]))){
		if (namenr == 1) {
		    sprintf(names[*usern-1], "%s.%s.%d",
			    shortname, kappashort, namenr);
		} 
		namenr++; 
          	sprintf(names[*usern], "%s.%s.%d", shortname, kappashort,namenr); 
	    } else namenr=1;			
	  }

	  if (PL > 7) {
//	    PRINTF(" %d mem %ld; addr %ld; %s, %d %d\n", *usern,
//		   (POINTER) mem[*internn], (POINTER)address[*internn], 
//		   names[*usern], sorts[*usern], internsorts[*internn]);
	  } 
 
	  if (user) {
	    if (printing > 0) {
	      if (printing <= 2) PRINTF("%d",*usern + 1); 
	      else PRINTF("!%d", *usern + 1);
	    }
	    (*usern)++;
	  }
	  (*internn)++;
	  nv++;
	} else {
	  // kein NA an der Position; somit nur Anzeige
	  if (printing > 1) {
	    if (c>0) PRINTF(", "); else leer(depth+1);
	    if (printing <= 2)  PRINTF("*");
	    else {
	      if (type[i]== REALSXP) PRINTF("%6.2e", v);
	      else  PRINTF("%5d", (int) v);
	    }
	  }
	}  // else kein NA 
      } //c
      if ((printing > 0 && nv > 0) || (printing > 1 && ncol[i] > 0)) 
	  PRINTF("\n"); 
    } // r
  } // i

 
  cov_model *sub;
  for (i=0; i<MAXSUB; i++) {
    sub = cov->sub[i];
    if (sub != NULL) {
      if (printing > 0) {
	leer(depth);
	PRINTF("%s = ", C->subnames[i]);
      }
      GetNAPosition(sub, usern, internn, mem, address, 
		    names, sorts, internsorts, isnan,
		    covzaehler, printing, depth);
      if (cov->anyNAdown==TriFalse && sub->anyNAdown==TriTrue)
	 cov->anyNAdown = TriTrue;
    }
  }

  if (cov->anyNAdown == TriTrue && cov->anyNAscaleup != TriTrue) {
     for (i=0; i<MAXSUB; i++) {
      sub = cov->sub[i];
      if (sub != NULL) {
	  if (sub->anyNAdown==TriFalse && sub->anyNAscaleup != TriTrue)
	      sub->anyNAdown=TriMaxFalse;
      }
    }
  }
}



SEXP GetNAPositions(SEXP model, SEXP tsdim, SEXP xdim, SEXP stationary,
 		    SEXP Print) {
    int i, usern, internn, covzaehler[MAXNRCOVFCTS];
  naptr_type mem, address;
  sortsofparam sorts[MAX_NA], internsorts[MAX_NA];
  bool isnan[MAX_NA];
  NAname_type names;

  bool skipchecks = GLOBAL.general.skipchecks;
  GLOBAL.general.skipchecks = true;
  CheckModel(model, tsdim, xdim, stationary, STORED_MODEL + MODEL_USER,
      MAXMLEDIM);  
  
  GLOBAL.general.skipchecks = skipchecks;
  usern = internn = 0;
  for (i=0; i<MAXNRCOVFCTS;  covzaehler[i++]=0);
  GetNAPosition(STORED_MODEL[MODEL_USER], &usern, &internn, mem, address,
		names,
		sorts, internsorts, isnan,  
		covzaehler, INTEGER(Print)[0], 0);  
  SEXP ans;
  PROTECT (ans =  allocVector(INTSXP, 1));
  INTEGER(ans)[0] = usern;
  UNPROTECT(1);
  return ans;
}






int countnas(cov_model *cov) {
    int i, r,  count, endfor,
      *nrow =  cov->nrow,
      *ncol = cov->ncol;
    cov_fct *C = CovList + cov->nr;
  SEXPTYPE *type = C->kappatype;

  if (( cov->nr == MIXEDEFFECT || cov->nr == MLEMIXEDEFFECT )
      && cov->nrow[1] > 0) // b is given
      return 0;

  count= 0;
  for (i=0; i<C->kappas; i++) {
    if (nrow[i] == 0 || ncol[i] == 0) continue;
    endfor = nrow[i] * ncol[i];
    if (type[i] == REALSXP) { 
      double *p = cov->p[i];
      for (r=0; r<endfor; r++) if (ISNAN(p[r]) || ISNA(p[r])) count++;
    } else if (type[i] == INTSXP) {
	int *p = (int *) cov->p[i];
	for (r=0; r<endfor; r++) if (p[r] == NA_INTEGER) count++;
    } else if (type[i] == LISTOF + REALSXP) {
	continue; // no NAs allowed
    } else assert(false);
  }
  
  cov_model *sub;
  for (i=0; i<MAXSUB; i++) {
    sub = cov->sub[i];
    if (sub == NULL) continue; 
    count += countnas(sub);
  }

  return count;
}



void Take21internal(cov_model *cov, cov_model *cov_bound,
		       double **bounds_pointer, int *NBOUNDS) {
  /* determine the bounds for the parameters that 
     are to be estimated
   */


//    PrintModelInfo(cov_bound);
//    assert(false);

  int i, c, r, nv=0,
    *nrow = cov->nrow,
    *ncol = cov->ncol,
    *nrow2 = cov_bound->nrow,
    *ncol2 = cov_bound->ncol;
  cov_fct *C = CovList + cov->nr;
  SEXPTYPE *type = C->kappatype;
  
//  printf("\n\n\n ======================= \n");
//  PrintModelInfo(cov);
//  PrintModelInfo(cov_bound);

  for (i=0; i<C->kappas; i++) {
      //     printf("%s %s \n", C->name, C->kappanames[i]);
  if (C->kappatype[i] >= LISTOF) continue;
    
    if (nrow[i] != nrow2[i] || ncol[i] != ncol2[i]) {
        PRINTF("%s i: %d, nrow1=%d, nrow2=%d, ncol1=%d, ncol2=%d\n", 
	       C->name, i, nrow[i], nrow2[i], ncol[i], ncol2[i]);
  	error("lower/upper/user does not fit to model (size of matrix)");
    }
    
    for (r=0; r<nrow[i]; r++) {
      for (c=0; c<ncol[i]; c++) {
	  double v=RF_NAN,
	      w=RF_NAN; // value in aktuellem parameter
	int idx = c * nrow[i] + r;
	
//	printf("%d %d idx=%d %d %d %d\n", r,c, idx, type[i], REALSXP, INTSXP);

	if (type[i] == REALSXP) {
	  v = cov->p[i][idx];
	  w = cov_bound->p[i][idx];
	}
	else if (type[i] == INTSXP) {
	  v = ((int *) cov->p[i])[idx] == NA_INTEGER 
	      ? NA_REAL : (double) ((int *) cov->p[i])[idx];;
	  w = ((int *) cov_bound->p[i])[idx] == NA_INTEGER 
	      ? NA_REAL : (double) ((int *) cov_bound->p[i])[idx];	  
	}
	
//	if (cov->nr >= DOLLAR && cov->nr <= LASTDOLLAR)
//	    printf("%s ", CovList[cov->sub[0]->nr].name);
	//             printf("%s %s, r=%d, c=%d: %d %d %f\n",
//	       C->name, C->kappanames[i], r, c, nrow[i], ncol[i], v);
       

	if (ISNA(v) || ISNAN(v)) { // entgegen Arith.h gibt ISNA nur NA an !!
	    if ((cov->nr < DOLLAR || cov->nr > LASTDOLLAR || 
		 i == DVAR ||
		 (i== DSCALE && cov->q == NULL) || // ! natscaling 
		 i == DANISO) &&
		((cov->nr != MIXEDEFFECT && cov->nr != MLEMIXEDEFFECT) || 
		 i != BETAMIXED)
                )// aniso ?? ABfrage OK ??
	   {
//		 printf("%s %s, r=%d, c=%d: %d <? %d\n",
//			C->name, C->kappanames[i], r, c, nv, *NBOUNDS);
 	     if (nv >= *NBOUNDS) {
//		 PRINTF("%s %s, r=%d, c=%d: %d >= %d\n",
//			C->name, C->kappanames[i], r, c, nv, *NBOUNDS);
  	       error("lower/upper/user does not fit to model (number parameters)");
	     }
	     (*bounds_pointer)[nv] = w;
	     nv++;
           }
	}  //ISNA
      } //c
    } // r
  } // i
  *NBOUNDS = *NBOUNDS - nv;
  (*bounds_pointer) += nv;

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      Take21internal(cov->sub[i], cov_bound->sub[i], bounds_pointer, NBOUNDS);
    }
  }
}


SEXP Take2ndAtNaOf1st(SEXP model, SEXP model_bound, SEXP tsdim, SEXP xdim, 
		      SEXP stationary, SEXP nbounds, SEXP skipchecks) {
    // MAXMLEDIM
  double *bounds_pointer;
  int NBOUNDS_INT =INTEGER(nbounds)[0],
      *NBOUNDS= &NBOUNDS_INT;
  bool oldskipchecks = GLOBAL.general.skipchecks;

  NAOK_RANGE = true;
  if (LOGICAL(skipchecks)[0]) GLOBAL.general.skipchecks = true;
  CheckModel(model_bound, tsdim, xdim, stationary, STORED_MODEL + MODEL_BOUNDS,
      MAXMLEDIM);
  GLOBAL.general.skipchecks = oldskipchecks;

//  PrintModelInfo(STORED_MODEL[MODEL_BOUNDS]); assert(false);

  NAOK_RANGE = true;
  CheckModel(model, tsdim, xdim, stationary, STORED_MODEL + MODEL_INTERN,
      MAXMLEDIM);

//  print("take\n");
//  PrintModelInfo(STORED_MODEL[MODEL_INTERN]);


  NAOK_RANGE = false;
  SEXP bounds;
  PROTECT(bounds = allocVector(REALSXP, *NBOUNDS));
  bounds_pointer = NUMERIC_POINTER(bounds);
//  printf("%f %f %f\n", bounds_pointer[0], bounds_pointer[1], 
//	 bounds_pointer[2]);
//  assert(false);

//  printf("s %ld\n", bounds_pointer);
  Take21internal(STORED_MODEL[MODEL_INTERN], STORED_MODEL[MODEL_BOUNDS], 
		    &bounds_pointer, NBOUNDS);
//  printf("e %ld\n", bounds_pointer);
// assert(false);

  if (*NBOUNDS != 0) error("lower/upper does not fit to model");
  UNPROTECT(1);
  return bounds;
}




void GetNARanges(cov_model *cov, cov_model *min, cov_model *max, 
		 double *minpile, double *maxpile, int *usern) {
    /* 
       determine the ranges of the parameters to be estimated 
     */
  int i,  r,
    *nrow =  cov->nrow,
    *ncol = cov->ncol;
  cov_fct *C = CovList + cov->nr;
  SEXPTYPE *type = C->kappatype;
  double 
    v = RF_NAN,
    dmin = RF_NAN,
    dmax = RF_NAN;

  for (i=0; i<C->kappas; i++) {
    int end = nrow[i] * ncol[i];
    
    if (end == 0) continue;

    if (type[i] == REALSXP) {
	dmin = min->p[i][0];
	dmax = max->p[i][0];
    } else if (type[i] == INTSXP) {
	dmin = ((int *) min->p[i])[0] == NA_INTEGER 
	    ? NA_REAL : (double) ((int *) min->p[i])[0];
	dmax = ((int *) max->p[i])[0]  == NA_INTEGER 
	    ? NA_REAL : (double) ((int *) max->p[i])[0];
    } else if (type[i] == LISTOF + REALSXP) {
	dmin = min->p[i][0];
	dmax = max->p[i][0];
    } else {
      error("unknown SXP type in GetRanges.");
    }
    
    for (r=0; r<end; r++) {
      if (type[i] == REALSXP) {
	v = cov->p[i][r];
      }
      else if (type[i] == INTSXP) {
        v = ((int *) cov->p[i])[r] == NA_INTEGER 
	  ? NA_REAL : (double) ((int *) cov->p[i])[r];
      } else if (type[i] == LISTOF + REALSXP) {
	  continue;  // !!!!!!!!!!!
      } else assert(false);

      if ( (ISNA(v) || ISNAN(v)) 
	   && cov->nr!=MIXEDEFFECT && cov->nr!=MLEMIXEDEFFECT) {
	if (cov->nr < DOLLAR || cov->nr > LASTDOLLAR || 
	    i == DVAR ||
	    (i== DSCALE && cov->q == NULL) || // ! natscaling 
	    i == DANISO // aniso ?? ABfrage OK ??
	    ) {
	    minpile[*usern] = dmin;
	    maxpile[*usern] = dmax;
//	    printf("%f %f\n", dmin, dmax);
	    (*usern)++;
	} else {
	    // maybe internn and internpile should be set here
	  continue;
	}
      } // isna
    } // r
  } // kappas

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      GetNARanges(cov->sub[i], min->sub[i], max->sub[i], minpile, maxpile, 
		  usern);
    }
  }
}



void CovMatrixMLE(double *x, int *dist, int *lx, int *set, int *submodel, // k
		  double *value) {
    LIST_ELEMENT = *set - 1;
    ELEMENTNR_PLUS = *submodel - 1;
    CovarianceMatrix(x, (bool) *dist, *lx, STORED_MODEL[MODEL_MLE], value);
    ELEMENTNR_PLUS = LIST_ELEMENT = - 1;
}


int check_recursive_range(cov_model *cov, bool NAOK) {
    /*
      called by MLEGetModelInfo
     */
  int i, err;
  if ((err = check_within_range(cov, NAOK)) != NOERROR) return err;
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      if ((err = check_recursive_range(cov->sub[i], NAOK)) != NOERROR)
	return err;
    }
  }
  return NOERROR;
}


int CheckEffect(cov_model *cov) {
  int i,j,end;
  bool na_var = false;
  double *p;
  cov_model *next;
  if (cov->nr >= GATTER && cov->nr <= LASTGATTER) return CheckEffect(cov->sub[0]);
  else if (cov->nr == MIXEDEFFECT) {    
    cov->nr = MLEMIXEDEFFECT;  
    if (cov->nrow[1] == 0 && cov->nsub == 0) return deteffect;
    if (cov->nsub == 0) return fixedeffect;
    next = cov->sub[0];
    if (next->nr >= GATTER && next->nr <= LASTGATTER) next = next->sub[0];
    if (next->nr >= DOLLAR && next->nr <= LASTDOLLAR) { 
      if (next->ncol[DVAR] == 1 && next->nrow[DVAR] == 1) {
        na_var = (ISNA(next->p[DVAR][0]) || ISNAN(next->p[DVAR][0]));
      }
      for (i=0; i<=DMAX; i++) { 
        if (i!=DVAR) {
	  end = next->ncol[i] * next->nrow[i];
	  p = next->p[i];
	  for (j=0; j<end; j++) {
	    if (ISNA(p[j]) || ISNAN(p[j])) {
	      return (next->nr == CONSTANT) ? eff_error : spaceeffect;
		  // NA in e.g. scale for constant matrix -- does not make sense
	    }
	  }
	}
      }
      next = next->sub[0];
      if (next->nr >= GATTER && next->nr <= LASTGATTER) next = next->sub[0];
    }
    if (next->nr == CONSTANT) {
       return (cov->nrow[0] > cov->ncol[0])
	   ? (na_var ? rvareffect : randomeffect)
	   : (na_var ? lvareffect : largeeffect);
    } else return spaceeffect;
  } else return remaining;
}


static naptr_type MLE_MEM, MLE_ADDRESS;
static int MLE_USERN=-1, MLE_INTERNN = -1;
static sortsofparam MLE_SORTS[MAX_NA], MLE_INTERNSORTS[MAX_NA];
static bool MLE_ISNAN[MAX_NA];


SEXP MLEGetModelInfo(SEXP model, SEXP tsdim, SEXP xdim) { // ex MLEGetNAPos
    /*
      basic set up of covariance model, determination of the NA's
      and the corresponding ranges.
     */
    cov_model *cov, *min=NULL,
    *max=NULL, 
    *openmin=NULL, 
    *openmax=NULL;
//  bool skipchecks = GLOBAL.general.skipchecks;
  NAname_type names;
  double mle_min[MAX_NA], mle_max[MAX_NA];
  int i, covzaehler[MAXNRCOVFCTS],
      usern = 0, effect[MAXSUB], nas[MAXSUB];
  bool anyeffect;
  SEXP stationary, isotropy, Effect, Nas;
#define total 5

  PROTECT(stationary=allocVector(LGLSXP, 1));
  LOGICAL(stationary)[0] = false;
  PROTECT(isotropy=allocVector(LGLSXP, 1));
  LOGICAL(isotropy)[0] = false;

  NAOK_RANGE = true;

//  printf("here\n");
  CheckModel(model, tsdim, xdim, stationary, STORED_MODEL + MODEL_MLE,
      MAXMLEDIM);
  NAOK_RANGE = false;

  // printf("here\n");
 

  cov = STORED_MODEL[MODEL_MLE];
  if (cov->nr >= GATTER && cov->nr <= LASTGATTER) {
    cov_model *next = cov->sub[0];
    if (next->statIn < COVARIANCE) {
      LOGICAL(stationary)[0] = true;
      if (next->isoIn==ISOTROPIC) {
	LOGICAL(isotropy)[0] = true;
	INTEGER(xdim)[0] = 1;
//	removeOnly(STORED_MODEL + MODEL_MLE); nicht wegnehmen
// da derzeit negative distanczen in CovMatrixMult auftreten
      }
    }
  }


  // GLOBAL.general.skipchecks = skipchecks;
  check_recursive_range(STORED_MODEL[MODEL_MLE], true);

  get_ranges(STORED_MODEL[MODEL_MLE], &min, &max, &openmin, &openmax, 
	     true);

//  printf("min/max\n");
//  PrintModelInfo(min);
//  PrintModelInfo(max);
  
//x  PrintModelInfo(cov); assert(false);
 

  MLE_USERN = MLE_INTERNN = 0;
  for (i=0; i<MAXNRCOVFCTS;  covzaehler[i++]=0);
  GetNAPosition(STORED_MODEL[MODEL_MLE], &MLE_USERN, &MLE_INTERNN,
		MLE_MEM, MLE_ADDRESS, names, 
		MLE_SORTS, MLE_INTERNSORTS,  MLE_ISNAN, 		
		covzaehler,
		(PL > 2) + (PL > 4) + (PL > 6), 0); 
  GetNARanges(STORED_MODEL[MODEL_MLE], min, max, mle_min, mle_max, 
	      &usern);
 
  cov = STORED_MODEL[MODEL_MLE];
  if (cov->nr >= GATTER && cov->nr <= LASTGATTER) cov = cov->sub[0];

  anyeffect=false;
  if (cov->nr == PLUS) {  

     for (i=0; i<cov->nsub; i++) {
	 effect[i] = CheckEffect(cov->sub[i]);
	 nas[i] = countnas(cov->sub[i]);

      //    printf("%d %d %d\n", i, effect[i], remaining);

      if (effect[i] == eff_error) 
        ERR("scaling parameter appears with constant matrix in the mixed effect part");
      if (effect[i] != remaining) anyeffect = true;
     }
  }


  PROTECT(Effect = allocVector(INTSXP, anyeffect  ? cov->nsub : 0));
  PROTECT(Nas = allocVector(INTSXP, anyeffect  ? cov->nsub : 0));
  if (anyeffect) {
    for (i=0; i<cov->nsub;i++) INTEGER(Effect)[i] = effect[i];
    for (i=0; i<cov->nsub;i++) INTEGER(Nas)[i] = nas[i];
    assert(cov->nr == PLUS); // zur sicherheit
    cov->nr = MLEPLUS;
  }
  // printf("xx %d %d %d %d\n", effect[0], effect[1], anyeffect, LENGTH(Effect));
     

  COV_DELETE(&min);
  COV_DELETE(&max);
  COV_DELETE(&openmin);
  COV_DELETE(&openmax);

  SEXP ans, matrix, nameAns, nameMatrix;
  PROTECT(matrix =  allocMatrix(REALSXP, MLE_USERN, 4));
  PROTECT(nameMatrix = allocVector(STRSXP, MLE_USERN));
  for (i=0; i<MLE_USERN; i++) {
    REAL(matrix)[i] = mle_min[i];
    REAL(matrix)[i + MLE_USERN] = mle_max[i];
    REAL(matrix)[i + 2 * MLE_USERN] = MLE_SORTS[i];
    REAL(matrix)[i + 3 * MLE_USERN] = (int) MLE_ISNAN[i];
    SET_STRING_ELT(nameMatrix, i, mkChar(names[i]));
  }
  setAttrib(matrix, R_RowNamesSymbol, nameMatrix);

  PROTECT(ans = allocVector(VECSXP, total));
  PROTECT(nameAns = allocVector(STRSXP, total));
  i = 0;
  SET_STRING_ELT(nameAns, i, mkChar("minmax"));
  SET_VECTOR_ELT(ans, i++, matrix);
  SET_STRING_ELT(nameAns, i, mkChar("stationary"));
  SET_VECTOR_ELT(ans, i++, stationary);  
  SET_STRING_ELT(nameAns, i, mkChar("isotropy"));
  SET_VECTOR_ELT(ans, i++, isotropy);  
  SET_STRING_ELT(nameAns, i, mkChar("effect"));
  SET_VECTOR_ELT(ans, i++, Effect);  
  SET_STRING_ELT(nameAns, i, mkChar("NAs"));
  SET_VECTOR_ELT(ans, i++, Nas);  
  setAttrib(ans, R_NamesSymbol, nameAns);
  assert(i==total);

  UNPROTECT(8);

  return ans;
}
 

bool anymixed(cov_model *cov) {
  int i;
  if (cov->nr == MIXEDEFFECT || cov->nr == MLEMIXEDEFFECT) return true;
  for (i=0; i<cov->nsub; i++) {
    if (anymixed(cov->sub[i])) return true;
  }
  return false;
}

 
void MLEanymixed(int *anymix) {
    *anymix = (int) anymixed(STORED_MODEL[MODEL_MLE]);
}

void iexplDollar(cov_model *cov) {
    /*    
      get the naturalscaling values and devide the preceeding scale model     
      by this value
    */
  double *p, scale;
  if (cov->nr == NATSC) {
      cov_model 
	  *next = cov->sub[0],
	  *calling = cov->calling;
      if (calling->nr >= GATTER && calling->nr <= LASTGATTER)
	  calling = calling->calling;
      assert(calling!=NULL && calling->nr >= DOLLAR && calling->nr <= LASTDOLLAR);
      scale = CovList[next->nr].naturalscale(next);

      p = calling->p[DSCALE];
      if (p != NULL) p[0] /= scale;
      else {
	  int i,
	      n = calling->nrow[DANISO] * calling->ncol[DANISO];
	  p = calling->p[DANISO];
	  for (i=0; i<n; i++) p[i] /= scale;
      }
  } else {
    int i;
    for (i=0; i<MAXSUB; i++) { // cov->sub[i]: luecken erlaubt bei PSgen !
      if (cov->sub[i] != NULL) iexplDollar(cov->sub[i]);
    }
  }
}
  

void expliciteDollarMLE(double *values) {
    // userinterfaces.cc 
  int i, un;
 	
  // first get the naturalscaling values and devide the preceeding 
  // scale model by this value 
  if (NS==NATSCALE_MLE) {
      iexplDollar(STORED_MODEL[MODEL_MLE]);
  }

  // Then get the values out of the model
  for (un=i=0; un<MLE_USERN; i++) {
      values[un++] = MLE_MEM[i][0];
      MLE_MEM[i][0] = RF_NAN;
  }
}



void PutValuesAtNA(double *values){
  int i, un;
  // set ordinary parameters for all (sub)models
  for (un=i=0; i<MLE_INTERNN; i++) MLE_MEM[i][0] = values[un++];
  if (un != MLE_USERN) ERR("severe error in PutValuesAtNA.");
 
  // finally set all aniso from scale
  //for (i=0; i<MLE_INTERNN; i++) {
  //  if (MLE_INTERNSORTS[i] == ANISOFROMSCALE) {
  //    MLE_MEM[i][0] = 1.0 / MLE_ADDRESS[i][0];
  //   }
  // }

//  PrintModelInfo(STORED_MODEL[MODEL_INTERN]);

}




#define MAX_INT   2000000
SEXP IGetModel(cov_model *cov, int modus) {

//      # modus: 1 : Modell wie gespeichert
//      #        0 : Modell unter Annahme PracticalRange>0
//      #        [ 2 : natscale soweit wie moeglich zusammengezogen ]

  SEXP model, nameMvec, dummy;
  int i, nmodelinfo,
    k = 0; 
  long int mem=0;
  cov_fct *C = CovList + cov->nr;

  if ((modus == 0 && cov->nr == NATSC) 
      || (cov->nr>=GATTER && cov->nr<=LASTGATTER)) { 
    return IGetModel(cov->sub[0], modus);
  }
  
  nmodelinfo = C->kappas + 1;
  for (i=0; i<MAXSUB; i++) if (cov->sub[i] != NULL) nmodelinfo++;
  for (i=0; i<C->kappas; i++) if (cov->p[i] == NULL) nmodelinfo--;

  PROTECT(model = allocVector(VECSXP, nmodelinfo));
  PROTECT(nameMvec = allocVector(STRSXP, nmodelinfo));

  SET_STRING_ELT(nameMvec, k, mkChar("")); // name
  cov_fct *CC = CovList + cov->nr;
  while(strncmp(CC->name, InternalName, strlen(InternalName)) ==0) CC--;
  SET_VECTOR_ELT(model, k++, mkString(CC->name));

  for(i=0; i<C->kappas; i++) {
    if (cov->p[i] == NULL) {	
	// k++; // 9.1.09, needed for "$" (proj can be NULL) -- 19.1.09 k++ 
	// deleted since otherwise outcome does not fit 
	// model definition anymore (arbitrary NULL in the definition)
	// instead: nmodelinfo--; above wihtin loop
      continue;
    } else  dummy = Param((void*) cov->p[i], cov->nrow[i], cov->ncol[i], 
			  C->kappatype[i], true, &mem);    
    SET_STRING_ELT(nameMvec, k, mkChar(C->kappanames[i]));
    SET_VECTOR_ELT(model, k++, dummy);
  }

  // vielleicht mal spaeter !
//  SET_STRING_ELT(nameMvec, k, mkChar("user"));  
//  SET_VECTOR_ELT(model, k++, Int(cov->user, Nothing + 1, MAXINT, mem));

//  printf("a %d %d %d %s\n",cov->nsub, k, nmodelinfo, CovList[cov->nr].name);
//  PrintModelInfo(cov);

  int zaehler = 0;
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
// printf("i=%d %d %d %s\n", i,cov->nsub,zaehler,CovList[cov->sub[i]->nr].name);
      SET_VECTOR_ELT(model, k++, IGetModel(cov->sub[i], modus));
      if (++zaehler >= cov->nsub) break;
    }
  }

//  PrintModelInfo(cov);
//  printf("b %d %d %s\n", k, nmodelinfo, CovList[cov->nr].name);

  assert(k == nmodelinfo);

  setAttrib(model, R_NamesSymbol, nameMvec);
  UNPROTECT(2); // model + namemodelvec

  return model;
}

SEXP GetModel(SEXP keynr, SEXP Modus) {
//    # modus: 1 : Modell wie gespeichert
//    #        0 : Modell unter Annahme PracticalRange>0 (natsc werden geloescht)
//    #        2 : natscale soweit wie moeglich zusammengezogen (natsc werden
//                 drauf multipliziert)

// Nutzer kann 3 Modifikationen des Models in MLE laufen lassen:
//      * keinen Praktikal range oder individuell angeben
//      * extern practical range definieren
//      * intern practical range verwenden lassen (use.naturalscaling)

  int knr = INTEGER(keynr)[0],
    modus = INTEGER(Modus)[0];
  cov_model *cov;
  if (knr>=0) {
    if (knr < MAXKEYS) {
      key_type *key;
      key = &(KEY[knr]);
      if (key->cov == NULL) return allocVector(VECSXP, 0);
      cov = key->cov;
    } else return allocVector(VECSXP, 0);
  } else {
    knr = -knr-1;
    if (knr < MODEL_MAX && STORED_MODEL[knr] != NULL) cov = STORED_MODEL[knr];
    else return allocVector(VECSXP, 0);
  }
  if (modus == 0 || modus == 1) return IGetModel(cov, modus);
  else {
    cov_model *dummy = NULL; //ACHTUNG: "=NULL" hinzugefuegt
    SEXP value;
    covcpy(&dummy, cov, false, true);
    iexplDollar(dummy);
    value = IGetModel(dummy, 0);
    COV_DELETE(&dummy);
    return(value);
  }
}

