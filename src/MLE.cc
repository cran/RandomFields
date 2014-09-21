/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields 

 Copyright (C) 2001 -- 2014 Martin Schlather, 

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

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Linpack.h>
#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
 
#include <string.h>
#include "RF.h"
#include "primitive.h"
// #include "Covariance.h"


bool is_top(cov_model *cov) {
  if (cov == NULL) BUG;
  return isInterface(cov);
}

void GetNAPosition(cov_model *cov, 
		   int *NAs, naptr_type mem, int *elmnts, elptr_type mem_elmnts,
		   NAname_type names, sortsofparam *sorts, bool *isnan,
		   covptr_type covModels, 
		   int *covzaehler, int allowforintegerNA,
		   int SHORTlen, int printing, int depth, bool no_variance,
		   bool excludetrends
		   ) {
  /*
    printing <= 0 nothing
    1 only NA, NaN
    2 others as "*"
    3 others as value

    mem[0..NAs] : addresse wo (neue) NA position im Speicher	   
    sorts[0..NAs]     : see sortsofparam in RF.h
    isnan[0..NAs] : NA or NaN
  */

  // PMI(cov);

  // if (SHORTlen < 2) BUG;
 
  int i, c, r, 
    namenr = 1,     
    *nrow =  cov->nrow,
    *ncol = cov->ncol;
  // double *lastmem = NULL;
  char shortname[255], shortD[255];
  cov_fct *C = CovList + cov->nr, // nicht gatternr
    *CC = C;
  SEXPTYPE *type = C->kappatype;
	  
  if (isRandom(cov)) {
    error("Random functions not allowed in models of RFfit, yet");
  }

  //printf("depth =%d %s %d\n", depth, C->name, *elmnts);


  // C is pointer to true model; CC points to the first model passed
  // in initNerror.cc by IncludeModel or IncludePrim
  while(strncmp(CC->name, InternalName, strlen(InternalName)) ==0) CC--;
  covzaehler[cov->nr]++;
  //  
  
  //  print("SHORT %d %s\n", SHORTlen, CC->name);

  strcopyN(shortname, CC->nick + 2, SHORTlen);

  if (covzaehler[cov->nr] >= 2) {
    char dummy[255];
    strcopyN(dummy, shortname, SHORTlen-1);
    sprintf(shortname, "%s%d", dummy, covzaehler[cov->nr]);
  }
  if (printing>0) PRINTF("%s\n", CC->name); 
  // CC needed below for the kappa.names which are given


  //  print("%d:%s %d %s\n", cov->nr, CC->name, covzaehler[cov->nr], shortname);

  depth++;
  for (i=0; i<C->kappas; i++) {
    if (nrow[i] == 0 || ncol[i] == 0) continue;
    if (printing > 0) {
      leer(depth); PRINTF("%s\n", C->kappanames[i]);
    }
    if (i==0 && type[i] == INTSXP && strcmp(C->kappanames[i], ELEMENT) == 0){
      if (*elmnts >= MAX_MLE_ELMNTS ) 
	ERR("maximum number of models with variable elements reached");
      if (P0INT(i) == NA_INTEGER) ERR1("'%s' may not be NA", ELEMENT);
      mem_elmnts[(*elmnts)++] = PINT(i);
      // to do!! covModels[*NAs] hierher uebertragen !! s.u.

      //printf("elmnts %s %d %s %d\n", CC->name, i, C->kappanames[i], *elmnts);
      
      continue;
    }

    for (r=0; r<nrow[i]; r++) {
      int nv = 0; // anzahl NA in aktuellem parameter
      for (c=0; c<ncol[i]; c++) {
	double v; // value in aktuellem parameter
	int idx = c * nrow[i] + r;
	if (*NAs >= MAX_NA) ERR("maximum number of NA reached");

	// print("NAs=%d, %d %d %s\n", *NAs, i, idx, C->kappanames[i]);

	if (type[i] == REALSXP) {
	  v = P(i)[idx];
	  mem[*NAs] = P(i) + idx;
	  covModels[*NAs] = cov;
	} else if (type[i] == INTSXP) {
	  v = PINT(i)[idx] == NA_INTEGER 
	    ? RF_NA : (double) PINT(i)[idx];
	  if (ISNAN(v) && !allowforintegerNA) {
	    // crash(cov);
	    ERR("integer variables currently do not allow for NA"); // !!!
	  }
	  // mem[*NAs] = (*double) &(P0INT(i])[idx]);
	} else if (type[i] == LISTOF + REALSXP) {
	  listoftype *q;
	  int j, end;
	  double *p;
	  q= PLIST(i); 
	  p = q->p[r];
	  end = q->nrow[r] * q->ncol[r];
	  for (j=0; j<end; j++)
	    if (ISNAN(p[j])) ERR("no NAs allowed in regression ");
	  v = 0.0; // dummy
	} else if (type[i] == CLOSXP) {
	  v = 0.0; //dummy
	} else if (type[i] == LANGSXP) {
	  v = 0.0; //dummy
	} else {
	  v = RF_NA;
	  BUG;
	}

 	isnan[*NAs] = ISNAN(v) && !ISNA(v);
	if (ISNAN(v)) { // entgegen Arith.h gibt ISNA nur NA an !!
	  if (isRandom(cov)){
	    NotProgrammedYet("NA in random functions");
	  }
	  if (printing > 0) {
	    if (nv>1 || (c>0 && printing > 1)) PRINTF(", "); else leer(depth+1);
	  }
	  if (isDollar(cov)) {
	    // shortD partial name for R level
	    cov_model *next = cov->sub[0];
	    while(isDollar(next) || isNatsc(next)) next = next->sub[0];// ok
	    if (covzaehler[next->nr] == 0) { // next wurde noch nicht
	      // untersucht, somit ist covzaehler[next->nr] um 1 niedriger
	      // als covzaehler[cov->nr] !
	      strcopyN(shortD, NICK(next) + 2, SHORTlen);
	    } else {
	      char dummy[255];
	      strcopyN(dummy, NICK(next) + 2, SHORTlen-1);
	      sprintf(shortD, "%s%d", dummy, covzaehler[next->nr]+1);
	      //	      print("$$ > %d:%s %d %s\n", next->nr, NICK(to), 
	      //		     covzaehler[next->nr], shortD);
	    }
	    
	    //  assert(DANISO == DSCALE + 1);

	    if (i==DVAR) {
	      cov_model *call = cov->calling;
	      if (call == NULL) BUG;
	      if (!is_top(call) && call->nr == MIXEDEFFECT) {
		sorts[*NAs] = !is_top(call->calling) && 
		  (call->calling->nr == PLUS || call->calling->nr == SELECT)
		  && is_top(call->calling->calling) ? MIXEDVAR : ANYPARAM;
	      } else {
		while (!is_top(call) && (call->nr == PLUS || call->nr==SELECT))
		  call = call->calling;	
		if (call == NULL) BUG;
		sorts[*NAs] = is_top(call)
		  ? (next->nr == NUGGET ? NUGGETVAR : VARPARAM)
		  : ANYPARAM; 
	      }
	      if (no_variance && sorts[*NAs] <= SIGNEDSDPARAM)
		sorts[*NAs] = ANYPARAM;
	      sprintf(names[*NAs], 
		      sorts[*NAs] <= SIGNEDSDPARAM || sorts[*NAs] == NUGGETVAR 
		      || sorts[*NAs]==MIXEDVAR ? "%s.var" : "%s.vx",
		      shortD); // for R level only
	    } else  {
	      if (i==DSCALE) {
		//	        lastmem = mem[*NAs]; // used for all diagonal elements of
		//                          aniso
		sorts[*NAs] = SCALEPARAM;	
		sprintf(names[*NAs], "%s.s", shortD);// for R level only
	      } else {
		assert(i == DANISO);
		// no lastmem here !
		
		if (cov->q != NULL && cov->q[0] != 0.0) {
		  //		    print("sdfasdf\n");
		  ERR("naturalscaling only allowed for isotropic setting");
		}
		if (r==c) { 
		  sorts[*NAs] = DIAGPARAM; 
		  sprintf(names[*NAs], "%s.d.%d", shortD, r);// R level only
		} else {
		  sorts[*NAs] = ANISOPARAM;
		  sprintf(names[*NAs], "%s.d.%d.%d", shortD, r, c);// R level	
		}	
	      }
	    }
	  } else {
	    // standard setting !
	    if (cov->nr==MIXEDEFFECT || cov->nr==TREND) {
	      // || cov->nr==MLEMIXEDEFFECT)
	      if (excludetrends) continue;
	      assert(type[i] == REALSXP);	      
	      sorts[*NAs] = cov->nr==TREND ? TRENDPARAM : MIXEDPARAM;
	    } else {
	      if (type[i] == INTSXP)
		sorts[*NAs] = INTEGERPARAM; 
	      else {
		// r)ow, c)olumn; i:kappa, cov->nr, C
		sorts[*NAs] = C->paramtype(i, r, c);  // ANYPARAM; 
		if (sorts[*NAs] == IGNOREPARAM || 
		    sorts[*NAs] == DONOTRETURNPARAM) continue;

		//		print("%s k=%d no_var=%d, %d \n", C->name, i,
		//		       no_variance, sorts[*NAs]);

		if (no_variance && sorts[*NAs] <= SIGNEDSDPARAM)
		  sorts[*NAs] = ANYPARAM;
	      }
	    }
	    char kappashort[255];
	    strcopyN(kappashort, CC->kappanames[i], SHORTlen);
	    sprintf(names[*NAs], "%s.%s", shortname,  kappashort);
	    if (*NAs > 0 && 0==
	        strncmp(names[*NAs], names[*NAs-1], strlen(names[*NAs]))){
	      if (namenr == 1) {
		sprintf(names[*NAs-1], "%s.%s.%d",
			shortname, kappashort, namenr);
	      } 
	      namenr++; 
	      sprintf(names[*NAs], "%s.%s.%d", shortname, kappashort,namenr); 
	    } else namenr=1;			
	  }
 
	  if (printing > 0) {
	    if (printing <= 2) PRINTF("%d",*NAs + 1); 
	    else PRINTF("!%d", *NAs + 1);
	  }
	  //	  print("%d %ld\n", *NAs, mem[*NAs]);
	  (*NAs)++;
	  nv++;
	} else { // not is.na(v)
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
      GetNAPosition(sub, NAs, mem, elmnts, mem_elmnts,
		    names, sorts, isnan,
		    covModels,
		    covzaehler, allowforintegerNA, 
		    SHORTlen, printing, depth, 
		    C->paramtype(-i-1, 0, 0) != VARPARAM,
		    excludetrends );
    }
  }

}


/* wird (derzeit) nicht verwenden !! */
SEXP GetNAPositions(SEXP model_reg, SEXP model, SEXP spatialdim, SEXP Time,
		    SEXP xdimOZ, SEXP integerNA, SEXP Print) {
  int i, NAs, elmnts,  covzaehler[MAXNRCOVFCTS];
  naptr_type mem;
  elptr_type mem_elmnts;
  covptr_type covModels;
  sortsofparam sorts[MAX_NA];
  bool isnan[MAX_NA];
  NAname_type names;

  bool skipchecks = GLOBAL.general.skipchecks;
  GLOBAL.general.skipchecks = true;
  
  CheckModelInternal(model, ZERO, ZERO, ZERO, INTEGER(spatialdim)[0], 
		     INTEGER(xdimOZ)[0], 1, 1, false, false, 
		     LOGICAL(Time)[0], 
		     KEY + INTEGER(model_reg)[0]);  
  sprintf(ERROR_LOC, "getting positions with NA: ");

  // print("global=%d\n", GLOBAL.fit.lengthshortname);
  
  GLOBAL.general.skipchecks = skipchecks;
  NAs = 0;
  for (i=0; i<MAXNRCOVFCTS;  covzaehler[i++]=0);
  GetNAPosition(KEY[MODEL_USER], &NAs, mem, &elmnts, mem_elmnts,
		names, sorts, isnan,  
		covModels,
		covzaehler,
		INTEGER(integerNA)[0],
		GLOBAL.fit.lengthshortname, 
		INTEGER(Print)[0], 
		0, false, true);  
  SEXP ans;
  PROTECT (ans =  allocVector(INTSXP, 1));
  INTEGER(ans)[0] = NAs;
  UNPROTECT(1);
  return ans;
}



int countnas(cov_model *cov, int level) {
    int i, r,  count, endfor,
      *nrow =  cov->nrow,
      *ncol = cov->ncol;
    cov_fct *C = CovList + cov->nr; // nicht gatternr
  SEXPTYPE *type = C->kappatype;

  if (cov->nr == MIXEDEFFECT && level==0) { // || cov->nr == MLEMIXEDEFFECT )
    if (cov->nrow[MIXED_BETA] > 0) return 0; // b is given
  }
  if (cov->nr == TREND && level==0) return 0;

  count= 0;
  for (i=0; i<C->kappas; i++) {
    if (nrow[i] == 0 || ncol[i] == 0 || C->paramtype(i, 0, 0) == IGNOREPARAM
	|| C->paramtype(i, 0, 0) == DONOTRETURNPARAM )
      continue;
    endfor = nrow[i] * ncol[i];
    if (type[i] == REALSXP) { 
      double *p = P(i);
      for (r=0; r<endfor; r++) if (ISNAN(p[r])) count++;
    } else if (type[i] == INTSXP) {
      int *p = PINT(i);
      for (r=0; r<endfor; r++) if (p[r] == NA_INTEGER) count++;
    } else if (type[i] == LISTOF + REALSXP) {
	continue; // no NAs allowed
    } else BUG;
  }
  
  cov_model *sub;
  for (i=0; i<MAXSUB; i++) {
    sub = cov->sub[i];
    if (sub == NULL) continue; 
    count += countnas(sub, level + 1);
  }

  return count;
}



void Take21internal(cov_model *cov, cov_model *cov_bound,
		    double **bounds_pointer, int *NBOUNDS) {
  /* determine the bounds for the parameters that 
     are to be estimated
   */

  int i, c, r,
    nv=0, // zaehlt NAs
    *nrow = cov->nrow,
    *ncol = cov->ncol,
    *nrow2 = cov_bound->nrow,
    *ncol2 = cov_bound->ncol;
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  SEXPTYPE *type = C->kappatype;
  
  if (strcmp(Nick(cov), Nick(cov_bound)) != 0) {
    // print("%s %s\n", CovList[cov->nr].nick, CovList[cov_bound->nr].nick);  
    ERR("models do not match.");
  }
//  print("\n\n\n ======================= \n");
//   print("%s %d %s\n", C->name, C->kappas, CovList[cov_bound->nr].name);
  

  for (i=0; i<C->kappas; i++) {
      //       print("i=%d %s %s alt=%s\n", i, C->name, C->kappanames[i], 
    //CovList[cov_bound->nr].name);
    if (C->kappatype[i] >= LISTOF ||
	C->paramtype(i, 0, 0) == IGNOREPARAM ||
	C->paramtype(i, 0, 0) == DONOTRETURNPARAM) continue;
    
    if (nrow[i] != nrow2[i] || ncol[i] != ncol2[i]) {
        PRINTF("%s i: %d, nrow1=%d, nrow2=%d, ncol1=%d, ncol2=%d\n", 
	       C->name, i, nrow[i], nrow2[i], ncol[i], ncol2[i]);
  	ERR("lower/upper/user does not fit the model (size of matrix)");
    }
    
    for (r=0; r<nrow[i]; r++) {
      for (c=0; c<ncol[i]; c++) {
	  double v=RF_NA,
	      w=RF_NA; // value in aktuellem parameter
	int idx = c * nrow[i] + r;
	
// print("%d %d idx=%d %d %d %d nv=%d %d\n", r,c, idx, type[i], REALSXP, INTSXP,
//	      nv, *NBOUNDS);

	if (type[i] == REALSXP) {
	  v = P(i)[idx];
	  w = PARAM(cov_bound, i)[idx];
	}
	else if (type[i] == INTSXP) {
	  v = PINT(i)[idx] == NA_INTEGER ? RF_NA : (double) PINT(i)[idx];
	  w = PARAMINT(cov_bound, i)[idx] == NA_INTEGER 
	    ? RF_NA : (double) PARAMINT(cov_bound, i)[idx];	  
	}
	     

	if (ISNAN(v)) { // entgegen Arith.h gibt ISNA nur NA an !!
	  if ((!isDollar(cov) || 
	       i == DVAR ||
	       (i== DSCALE && cov->q == NULL) || // ! natscaling 
		 i == DANISO) &&
		(cov->nr != MIXEDEFFECT && cov->nr != TREND)
                ) // aniso ?? ABfrage OK ??
	   {
	     //		 print("%s %s, r=%d, c=%d: %d <? %d\n",
	     //		C->name, C->kappanames[i], r, c, nv, *NBOUNDS);
	     if (nv >= *NBOUNDS) {
	       PRINTF("%s %s, r=%d, c=%d: %d >= %d\n",
	       		C->name, C->kappanames[i], r, c, nv, *NBOUNDS);
  	       ERR("lower/upper/user does not fit the model (number parameters)");
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


SEXP Take2ndAtNaOf1st(SEXP model_reg, SEXP model, SEXP model_bound,
		      SEXP spatialdim, SEXP Time, SEXP xdimOZ, 
		      SEXP nbounds, SEXP skipchecks) {
    // MAXMLEDIM
  double *bounds_pointer;
  int m,
    NBOUNDS_INT =INTEGER(nbounds)[0],
    *NBOUNDS= &NBOUNDS_INT,
    nr[2] = {INTEGER(model_reg)[0], MODEL_BOUNDS};
  SEXP bounds,
    models[2] = {model, model_bound};
  bool oldskipchecks = GLOBAL.general.skipchecks;

  if (nr[0] == nr[1]) error("do not use register 'model bounds'");
  
  NAOK_RANGE = true;
  if (LOGICAL(skipchecks)[0]) GLOBAL.general.skipchecks = true;
  for (m=1; m>=0; m--) { // 2er Schleife !!
    // assert(m==1);
    CheckModelInternal(models[m], ZERO, ZERO, ZERO, INTEGER(spatialdim)[0], 
		       INTEGER(xdimOZ)[0], 1, 1, false, false, 
		       LOGICAL(Time)[0], 
		       KEY + nr[m]);  
    GLOBAL.general.skipchecks = oldskipchecks;
  }
  NAOK_RANGE = false;

  PROTECT(bounds = allocVector(REALSXP, *NBOUNDS));
  bounds_pointer = NUMERIC_POINTER(bounds);

  //  PMI(KEY[nr[0]]);
  //  PMI(KEY[nr[1]]);

  Take21internal(KEY[nr[0]], KEY[nr[1]], &bounds_pointer, NBOUNDS);

  if (*NBOUNDS != 0) ERR("lower/upper does not fit to model");
  UNPROTECT(1);
  return bounds;
}




void GetNARanges(cov_model *cov, cov_model *min, cov_model *max, 
		 double *minpile, double *maxpile, int *NAs) {
    /* 
       determine the ranges of the parameters to be estimated 
     */
  int i,  r,
    *nrow = cov->nrow,
    *ncol = cov->ncol;
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  SEXPTYPE *type = C->kappatype;
  double v, dmin, dmax;

  //  print("\ngetnaranges\n");

  for (i=0; i<C->kappas; i++) {
    int end = nrow[i] * ncol[i];
    
    if (end == 0) continue;

    if (type[i] == REALSXP) {
      dmin = PARAM0(min, i);
      dmax = PARAM0(max, i);
    } else if (type[i] == INTSXP) {
      dmin = PARAM0INT(min, i) == NA_INTEGER 
	? RF_NA : (double) PARAM0INT(min, i);
      dmax = PARAM0INT(max, i) == NA_INTEGER 
	? RF_NA : (double) PARAM0INT(max, i);
    } else if (type[i] == LISTOF + REALSXP) {
      dmin = PARAM0(min, i);
      dmax = PARAM0(max, i);
    } else if (type[i] == CLOSXP) {
	dmin = 0.0;
	dmax = 0.0;
    } else if (type[i] == LANGSXP) {
	dmin = 0.0;
	dmax = 0.0;
    } else {
      BUG;
      dmin = dmax = RF_NA;
     }
    
    for (r=0; r<end; r++) {
      v = RF_NA;
      if (type[i] == REALSXP) {
	v = P(i)[r];
      }
      else if (type[i] == INTSXP) {
        v = PINT(i)[r] == NA_INTEGER ? RF_NA : (double) PINT(i)[r];
      } else if (type[i] == LISTOF + REALSXP) {
	  continue;  // !!!!!!!!!!!
      } else if (type[i] == CLOSXP) {
          continue;  // !!!!!!!!!!!
      } else if (type[i] == LANGSXP) {
          continue;  // !!!!!!!!!!!
      } else BUG;

      if (ISNAN(v) && C->paramtype(i, 0, 0) != IGNOREPARAM
	  && C->paramtype(i, 0, 0) != DONOTRETURNPARAM
	  && cov->nr!=MIXEDEFFECT && cov->nr!=TREND) {// cov->nr!=MLEMIXEDEFFECT
	if (!isDollar(cov) || (i!=DALEFT && i!=DPROJ)) {
	    minpile[*NAs] = dmin;
	    maxpile[*NAs] = dmax;
//	    print("%s %d %f %f\n", C->kappanames[i],r, dmin, dmax);
	    (*NAs)++;
	} else {
	  // maybe NAs and internpile should be set here
	  continue;
	}
      } // isna
    } // r
  } // kappas

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      GetNARanges(cov->sub[i], 
		  min->sub[i], max->sub[i], 
		  minpile, maxpile, 
		  NAs);
    }
  }
}



int check_recursive_range(cov_model *cov, bool NAOK) {
    /*
      called by SetAndGetModelInfo
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
  int i,j,end, nr;    
  bool na_var = false;
  double *p;
  cov_model *next;
  if (cov->nr == MIXEDEFFECT) {    
    // cov->nr = MLEMIXEDEFFECT;  
    if (cov->nsub == 0) {
      return (cov->nrow[MIXED_BETA] > 0 && ISNAN(P0(MIXED_BETA)))
	? fixedeffect : deteffect;
    }
    next = cov->sub[0];
    if (isDollar(next)) { 
      if (next->ncol[DVAR] == 1 && next->nrow[DVAR] == 1) {
        na_var = ISNAN(PARAM0(next, DVAR));
      }
      for (i=0; i<=DMAX; i++) { 
        if (i!=DVAR) {
	  end = next->ncol[i] * next->nrow[i];
	  p = PARAM(next, i);
	  for (j=0; j<end; j++) {
	    if (ISNAN(p[j])) {
	      return next->nr == CONSTANT ? eff_error :
		na_var ? spvareffect : spaceeffect;
		  // NA in e.g. scale for constant matrix -- does not make sense
	    }
	  }
	}
      }
      next = next->sub[0]; // ok
    }
  
    int xx =  next->nr == CONSTANT  
      ? (cov->nrow[CONSTANT_M] > cov->ncol[CONSTANT_M]
	 ? (na_var ? rvareffect : randomeffect)
	 : (na_var ? lvareffect : largeeffect)
	 )
      : (na_var ? spvareffect : spaceeffect);
    if (xx == spvareffect || xx == spaceeffect) {
      BUG; //APMI(next);
    }


    return next->nr == CONSTANT  
      ? (cov->nrow[CONSTANT_M] > cov->ncol[CONSTANT_M]
	 ? (na_var ? rvareffect : randomeffect)
	 : (na_var ? lvareffect : largeeffect)
	 )
      : (na_var ? spvareffect : spaceeffect);
  }

  if (cov->nr == TREND) {
    int trend,
      effect = eff_error;
    bool isna ;

    for (j=0; j<2; j++) {
      trend = (j==0) ? TREND_MEAN : TREND_LINEAR;  
      //  print("trend %d %d\n", j, trend);

      if ( (nr = cov->nrow[trend] * cov->ncol[trend]) > 0) {
	p = P(trend);
	isna = ISNAN(p[0]);
	if ((effect != eff_error) && ((effect == fixedtrend) xor isna))
	  SERR1("do not mix deterministic effect with fixed effects in '%s'", 
		NICK(cov));
	for (i = 1; i<nr; i++) {
	  if (ISNAN(p[i]) xor isna) 
	    SERR("mu and linear trend:  all coefficient must be deterministic or all must be estimated");
	}
	effect = isna ? fixedtrend : dettrend;
      }
    }

    for (j=0; j<2; j++) {
      int param;
      trend = (j==0) ? TREND_POLY : TREND_FCT;
      //print("trend %d\n", trend);
      param = (j==0) ? TREND_PARAM_POLY : TREND_PARAM_FCT;
      if (cov->nrow[trend] > 0) {
	if (effect != eff_error)
	  SERR("polynomials and free functions in trend may not be mixed with other trend definitions. Please use a sum of trends.");

	nr = cov->nrow[param] * cov->ncol[param];
	if (nr > 0) {
	  p = P(param);
	  isna = ISNAN(p[0]);
	  for (i = 1; i<nr; i++) {
	    if (ISNAN(p[i]) xor isna) 
	      SERR("the coefficients in trend must be all deterministic or all coefficient are estimated");
	  }
	  effect = isna ? fixedtrend : dettrend;
	} else effect = fixedtrend;
      }
    }

    return effect;
  }

  cov_model *sub = cov;
  bool Simple = true;
  if (isDollar(sub)) {
    Simple = PARAMisNULL(sub, DPROJ) && PARAMisNULL(sub, DANISO) && 
      PARAMisNULL(sub, DALEFT);
    sub = sub->sub[0];
  }
  
  if (isNatsc(sub)) sub = sub->sub[0];
  cov_fct *C = CovList + sub->nr; // nicht gatternr
  
  if (C->maxsub == 0) {
    if (isPosDef(C->Type) && 
	C->domain == XONLY && C->isotropy == ISOTROPIC &&
	(C->vdim == 1 || C->cov==nugget) && Simple) {
      return simple; 
    }
    return primitive;
  }
  return remaining;
}


static naptr_type MEMORY[MODEL_MAX+1];
static covptr_type MEM_COVMODELS[MODEL_MAX+1];
static elptr_type MEMORY_ELMNTS[MODEL_MAX+1];
int MEM_NAS[MODEL_MAX+1], MEM_ELMNTS[MODEL_MAX+1];

SEXP SetAndGetModelInfo(SEXP model_reg, SEXP model, SEXP spatialdim, 
			SEXP distances, SEXP ygiven,
			SEXP Time,
			SEXP xdimOZ, SEXP shortlen, SEXP allowforintegerNA, 
			SEXP excludetrend) { 
  // ex MLEGetNAPos
  sortsofparam MEM_SORTS[MAX_NA];
  bool anyeffect, MEM_ISNAN[MAX_NA];
    /*
      basic set up of covariance model, determination of the NA's
      and the corresponding ranges.
     */
  cov_model *cov, *key,
    *min=NULL,  *max=NULL, 
    *pmin=NULL,  *pmax=NULL, 
    *openmin=NULL,  *openmax=NULL;
//  bool skipchecks = GLOBAL.general.skipchecks;
  NAname_type names;
  double mle_min[MAX_NA], mle_max[MAX_NA], mle_pmin[MAX_NA], mle_pmax[MAX_NA],
    mle_openmin[MAX_NA], mle_openmax[MAX_NA];
  int i,  covzaehler[MAXNRCOVFCTS], effect[MAXSUB], nas[MAXSUB],
    newxdim, neffect, NAs,
    err = NOERROR, 
    modelnr = INTEGER(model_reg)[0];
  const char *colnames[8] =
    {"pmin", "pmax", "type", "is.nan", "min", "max", "openmin", "openmax"};
 SEXP trans_inv, isotropic, Sxdim,
   matrix, nameAns, nameMatrix, RownameMatrix, 
   ans=NULL;

#define total 7

  PROTECT(trans_inv=allocVector(LGLSXP, 1));
  LOGICAL(trans_inv)[0] = false;
  PROTECT(isotropic=allocVector(LGLSXP, 1));
  LOGICAL(isotropic)[0] = false;

  NAOK_RANGE = true;
 
  CheckModelInternal(model, ZERO, ZERO, ZERO, 
		     INTEGER(spatialdim)[0], INTEGER(xdimOZ)[0],
		     1, LOGICAL(ygiven)[0] ? 1 : 0, // lx, ly
		     false, //grid
		     LOGICAL(distances)[0], //distances
		     LOGICAL(Time)[0], 
		     KEY + modelnr);  

  //  APMI(KEY[modelnr]); //, "setandget");

  cov = key = KEY[modelnr];
  if (isInterface(cov)) {
    cov = cov->key == NULL ? cov->sub[0] : cov->key;
  }

  if (key->prevloc->timespacedim > MAXMLEDIM) 
    ERR("dimension too high in MLE");

  NAOK_RANGE = false;

  if (cov->pref[Nothing] == PREF_NONE) {
    err = ERRORINVALIDMODEL;
    goto ErrorHandling;
  }
  // print("here returnx tsdim %d %d %d\n", 
  // INTEGER(tsdim)[0], LOGICAL(domain)[0],  LOGICAL(isotropic)[0]);


  // print("here\n");
 
  //  PMI(cov);

  newxdim = INTEGER(xdimOZ)[0];
  if (cov->domprev == XONLY && isPosDef(cov->typus)) {
    LOGICAL(trans_inv)[0] = true;
    if (cov->isoown==ISOTROPIC) {
      LOGICAL(isotropic)[0] = true;
      newxdim = 1;
//	removeOnly(KEY + modelnr); nicht wegnehmen
// da derzeit negative distanczen in GLOBAL.CovMatrixMult auftreten -- obsolete ?!
    }
  }
  //  print("SD here returnx tsdim %d\n", INTEGER(tsdim)[0]);


  // GLOBAL.general.skipchecks = skipchecks;
  check_recursive_range(cov, true);

  // print("AS here returnx tsdim %d\n", INTEGER(tsdim)[0]);
  if ((err = get_ranges(cov, &min, &max, &pmin, &pmax, &openmin, &openmax))
      != NOERROR) goto ErrorHandling;

  //  print("hier %d\n", INTEGER(shortlen)[0]);

  MEM_ELMNTS[modelnr] = MEM_NAS[modelnr] = 0;
  for (i=0; i<MAXNRCOVFCTS;  covzaehler[i++]=0);

  // !!
  GetNAPosition(cov, MEM_NAS + modelnr, MEMORY[modelnr], 
		MEM_ELMNTS + modelnr, MEMORY_ELMNTS[modelnr],
		names, MEM_SORTS, MEM_ISNAN, 		
		MEM_COVMODELS[modelnr],
		covzaehler,
		INTEGER(allowforintegerNA)[0],
		INTEGER(shortlen)[0],
		(PL >= PL_COV_STRUCTURE) + (PL >=PL_DETAILS) + 
		           (PL >= PL_SUBDETAILS), 
		0, false, LOGICAL(excludetrend)[0]); 

// print("XX modelnr=%d\n", modelnr);
  NAs = 0; GetNARanges(cov, min, max, mle_min, mle_max, &NAs);
  //print("XXz modelnr=%d\n", modelnr);
  NAs = 0; GetNARanges(cov, pmin, pmax, mle_pmin, mle_pmax, &NAs);
  //print("XX-------- amodelnr=%d\n", modelnr);
  NAs = 0; GetNARanges(cov, openmin, openmax, mle_openmin, mle_openmax, &NAs);
 
  anyeffect=false;
  if (cov->nr == PLUS || cov->nr == SELECT) {  
     for (i=0; i<cov->nsub; i++) {
       //  print("checkeffect: %d %d\n", i, cov->nsub);
       
       effect[i] = CheckEffect(cov->sub[i]);
       nas[i] = countnas(cov->sub[i], 0);
       
       // print("i=%d nsub=%d eff=%d remain=%d nas=%d\n", i, cov->nsub, effect[i], remaining, nas[i]);
       
       if (effect[i] == eff_error) GERR("scaling parameter appears with constant matrix in the mixed effect part");

       anyeffect |= (effect[i] != remaining && effect[i] != primitive);
     }
  } else {
    effect[0] = CheckEffect(cov);
    nas[0] = countnas(cov, 0);
  }
  anyeffect = true; // always bring back effects!!
  
  if (anyeffect) {
    if (cov->nr == PLUS || cov->nr == SELECT)
      neffect = cov->nsub;
    else neffect = 1;
  } else neffect = 0;

  // print("There returnx tsdim %d\n", INTEGER(tsdim)[0]);


 
  PROTECT(Sxdim=allocVector(INTSXP, 1));
  INTEGER(Sxdim)[0] = newxdim;

  // print("hYere returnx tsdim %d\n", INTEGER(tsdim)[0]);
  PROTECT(matrix =  allocMatrix(REALSXP, MEM_NAS[modelnr], 8));
  PROTECT(RownameMatrix = allocVector(STRSXP, MEM_NAS[modelnr]));
  PROTECT(nameMatrix = allocVector(VECSXP, 2));
  for (i=0; i<MEM_NAS[modelnr]; i++) {
    REAL(matrix)[i] = mle_pmin[i];
    REAL(matrix)[i + MEM_NAS[modelnr]] = mle_pmax[i];
    REAL(matrix)[i + 2 * MEM_NAS[modelnr]] = MEM_SORTS[i];
    REAL(matrix)[i + 3 * MEM_NAS[modelnr]] = (int) MEM_ISNAN[i];
    REAL(matrix)[i + 4 * MEM_NAS[modelnr]] = mle_min[i];
    REAL(matrix)[i + 5 * MEM_NAS[modelnr]] = mle_max[i];
    REAL(matrix)[i + 6 * MEM_NAS[modelnr]] = mle_openmin[i];
    REAL(matrix)[i + 7 * MEM_NAS[modelnr]] = mle_openmax[i];
    SET_STRING_ELT(RownameMatrix, i, mkChar(names[i]));
    // print("i=%d %s %e %e\n", i, names[i], mle_pmin[i], mle_pmax[i]);
   }
      
  SET_VECTOR_ELT(nameMatrix, 0, RownameMatrix);
  SET_VECTOR_ELT(nameMatrix, 1, Char(colnames, 8));
  setAttrib(matrix, R_DimNamesSymbol, nameMatrix);
   
  PROTECT(ans = allocVector(VECSXP, total));
  PROTECT(nameAns = allocVector(STRSXP, total));
  i = 0;
  SET_STRING_ELT(nameAns, i, mkChar("minmax"));
  SET_VECTOR_ELT(ans, i++, matrix);
  SET_STRING_ELT(nameAns, i, mkChar("trans.inv")); // !!
  SET_VECTOR_ELT(ans, i++, trans_inv);  
  SET_STRING_ELT(nameAns, i, mkChar("isotropic"));
  SET_VECTOR_ELT(ans, i++, isotropic);  
  SET_STRING_ELT(nameAns, i, mkChar("effect"));
  SET_VECTOR_ELT(ans, i++, Int(effect, neffect, MAXINT));  
  SET_STRING_ELT(nameAns, i, mkChar("NAs"));
  SET_VECTOR_ELT(ans, i++, Int(nas, neffect, MAXINT));  
  SET_STRING_ELT(nameAns, i, mkChar("xdimOZ"));
  SET_VECTOR_ELT(ans, i++, Sxdim);  
  SET_STRING_ELT(nameAns, i, mkChar("matrix.indep.of.x"));
  SET_VECTOR_ELT(ans, i++, ScalarLogical(cov->matrix_indep_of_x));  
  setAttrib(ans, R_NamesSymbol, nameAns);
  assert(i==total);

  UNPROTECT(8);
  // print(" J here returnx tsdim %d\n", INTEGER(tsdim)[0]);

 ErrorHandling:
  if (min != NULL) COV_DELETE(&min);
  if (max != NULL) COV_DELETE(&max);
  if (pmin != NULL) COV_DELETE(&pmin);
  if (pmax != NULL) COV_DELETE(&pmax);
  if (openmin != NULL) COV_DELETE(&openmin);
  if (openmax != NULL) COV_DELETE(&openmax);
  if (err != NOERROR) XERR(err);
  
  return ans;
}
 


void expliciteDollarMLE(int* ModelNr, double *values) { // 
    // userinterfaces.cc 
  int i, un,
    modelnr = *ModelNr,
    NAs = MEM_NAS[modelnr];
 	
  // first get the naturalscaling values and devide the preceeding 
  // scale model by this value 
  if (NS==NATSCALE_MLE) {
    iexplDollar(KEY[modelnr], true);
  }

  // Then get the values out of the model
  for (un=i=0; un<NAs; i++) {
    values[un++] = MEMORY[modelnr][i][0];
    MEMORY[modelnr][i][0] = RF_NA;
  }
}



void PutValuesAtNAintern(int *reg, double *values, bool init){
  int i, un,
    NAs = MEM_NAS[*reg];
  cov_fct *C = NULL;
  cov_model *cov = NULL;
  gen_storage s;
  STORAGE_NULL(&s);
  s.check = false;
  // set ordinary parameters for all (sub)models
 
  // cov_model *cov = KEY[*reg];
  //cov = cov->sub[0];
  //  print("%s %d\n", NICK(cov), *reg);
  ///assert(isDollar(cov));
  //print("%ld %f %ld %f\n", cov->p[DVAR], cov->p[DVAR][0],
  //	 cov->p[DSCALE], cov->p[DSCALE][0]);

  // printf("putvaluesatna\n");
  for (un=i=0; i<NAs; i++) {
 //   print("reg=%d i=%d %d %ld %f\n", *reg, i, NAs, MEMORY[*reg][i], values[un]);
    //     print("mem=%ld %f\n", (long int) MEMORY[*reg][i], values[un]) ;  
    //if (*reg < 0 || *reg > MODEL_MAX) XERR(ERRORREGISTER);  

    MEMORY[*reg][i][0] = values[un++];
  }

  if (init)
    for (i=0; i<NAs; i++) {
      cov = MEM_COVMODELS[*reg][i];
      C = CovList + cov->nr;       
      if (i==0 || cov != MEM_COVMODELS[*reg][i-1]) {
	if (!isDummyInit(C->Init)) {
	  //	print("i=%d %s initalising (%f)\n", i, NICK(MEM_COVMODELS[*reg][i]),   MEMORY[*reg][i][0]);
	  
	  C->Init(cov, &s);
	  //print("i=%d %s initalised (%f)\n", i, NICK(MEM_COVMODELS[*reg][i]),  
	  //   MEMORY[*reg][i][0]);
	}
      }
      //  
    }
  
  //printf("done\n");

  int one = 1;
  setListElements(reg, &one, &one, &one);
}

void PutValuesAtNA(int *reg, double *values){
  PutValuesAtNAintern(reg, values, true);
}

void PutValuesAtNAnoInit(int *reg, double *values){
  PutValuesAtNAintern(reg, values, false);
}


void setListElements(int *reg, int *i, int *k, int* len_k) {
  int len = *len_k;
  if (*reg < 0 || *reg > MODEL_MAX) XERR(ERRORREGISTER);
  cov_model *cov = KEY[*reg];

  if (cov==NULL) ERR("register is not initialised by 'RFfit'");
  if (isInterface(cov)) {
    cov = cov->key == NULL ? cov->sub[0] : cov->key;
  }
  if (cov->nr == SELECT) {
    int j;
    if (len != cov->nrow[SELECT_SUBNR]) {
      PFREE(SELECT_SUBNR);
      PALLOC(SELECT_SUBNR, len, 1);
    }
    for (j=0; j<len; j++) {
      PINT(SELECT_SUBNR)[j] = k[j] - 1;
    }

  }
  int j,
    iM1 = *i - 1,
    nr = MEM_ELMNTS[*reg];
  for (j=0; j<nr; j++) {
    MEMORY_ELMNTS[*reg][j][0] = iM1;
    // printf("setlistel %d %d %ld\n", j, iM1, MEMORY_ELMNTS[*reg][j]);
  }
  
}
