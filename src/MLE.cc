/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields 

 Copyright (C) 2001 -- 2015 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
RO
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
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
#include "Coordinate_systems.h"


bool is_top(cov_model *cov) {
  if (cov == NULL) BUG;
  return isInterface(cov) || isProcess(cov);
}

int GetNAPosition(cov_model *cov, 
		   int *NAs, naptr_type mem,
		   //		   int *elmnts, elptr_type mem_elmnts,
		   NAname_type names, sortsofparam *sorts, 
		   int *vdim1, int *vdim2, bool *isnan, bool *bayesian,
		   covptr_type covModels, 
		   int *covzaehler, int allowforintegerNA,
		   int SHORTlen, int printing, int depth, bool no_variance,
		   bool excludetrends
		   ) {

  // if (*NAs == 0) printf("\n\ngetnapos %s\n\n", NAME(cov));

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
 
  sprintf(ERROR_LOC, "'%s' : ", NICK(cov));

  cov_model *sub;
  int i, c, r, 
    namenr = 1,     
    *nrow =  cov->nrow,
    *ncol = cov->ncol;
  // double *lastmem = NULL;
  char shortname[255], shortD[255];
  cov_fct *C = CovList + cov->nr, // nicht gatternr
    *CC = C;
  SEXPTYPE *type = C->kappatype;
	 
  if (SHORTlen >= 255) SHORTlen = 254;


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
    /*
    if (i==0 && type[i] == INTSXP && strcmp(C->kappanames[i], ELEMENT) == 0){
      if (*elmnts >= MAX_MLE_ELMNTS ) 
	ERR("maximum number of models with variable elements reached");
      if (P0INT(i) == NA_INTEGER) ERR1("'%s' may not be NA", ELEMENT);
      mem_elmnts[(*elmnts)++] = PINT(i);
      // to do!! covModels[*NAs] hierher uebertragen !! s.u.

      //printf("elmnts %s %d %s %d\n", CC->name, i, C->kappanames[i], *elmnts);
      
      continue;
    }
    */
    
    if ((sub = cov->kappasub[i]) == NULL || isRandom(sub)) {    

      // printf("elmnts %s %d %s\n", C->name, i, C->kappanames[i]);
      
      if (nrow[i] == 0 || ncol[i] == 0) continue;
  
      if (sub != NULL) {
	int size = nrow[i] * ncol[i];
	for (r =0; r < size; P(i)[r++] = RF_NA);
	// pmi(cov, 0);
	//printf("nas %d \n", *NAs);
      }

      if (printing > 0) {
	leer(depth); PRINTF("%s\n", C->kappanames[i]);
      }
      for (r=0; r<nrow[i]; r++) {
	int nv = 0; // anzahl NA in aktuellem parameter
	if (*NAs >= MAX_NA) SERR("maximum number of NA reached");
	for (c=0; c<ncol[i]; c++) {
	  double v; // value in aktuellem parameter
	  int idx = c * nrow[i] + r;
	  if (*NAs >= MAX_NA) SERR("maximum number of NA reached");
	  vdim1[*NAs] = r + 1;
	  vdim2[*NAs] = c + 1;
	  bayesian[*NAs] = sub != NULL;
	  // print("NAs=%d, %d %d %s\n", *NAs, i, idx, C->kappanames[i]);
	
	  if (type[i] == REALSXP) {
	    v = P(i)[idx];
	    mem[*NAs] = P(i) + idx;
	    covModels[*NAs] = cov;
	  } else if (type[i] == INTSXP) {
	    v = PINT(i)[idx] == NA_INTEGER ? RF_NA : (double) PINT(i)[idx];
	    if (ISNAN(v) && !allowforintegerNA) {
	      // crash(cov);
	      SERR("integer variables currently not allowed to be 'NA'"); // !!!
	    }
	    // mem[*NAs] = (*double) &(P0INT(i])[idx]);
	  } else if (type[i] == LISTOF + REALSXP) {
	    listoftype *q = PLIST(i);
	    double *p = q->lpx[r];
	    int j, 
	      end = q->nrow[r] * q->ncol[r];
	    for (j=0; j<end; j++)
	      if (ISNAN(p[j])) 
		SERR("no NAs allowed in arguments that can be lists");
	    v = 0.0; // dummy
	  } else if (isRObject(type[i])) {
	    v = 0.0; //dummy
	  } else {
	    v = RF_NA;
	    BUG;
	  }
	  
	  isnan[*NAs] = ISNAN(v) && !ISNA(v);
	  if (ISNAN(v)) { // entgegen Arith.h gibt ISNA nur NA an !!
	    if (printing > 0) {
	      if (nv>1 || (c>0 && printing > 1)) PRINTF(", ");
	      else leer(depth+1);
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
	    
	      // printf("cov %s %s\n", NICK(cov), shortD);
	      
	      //  assert(DANISO == DSCALE + 1);
	      
	      if (i==DVAR) { 		
		vdim1[*NAs] = vdim2[*NAs] = r + 1;
		cov_model *call = cov->calling;
		if (call == NULL) BUG;
		if (!is_top(call) && call->nr == MIXEDEFFECT) {		
		  sorts[*NAs] = !is_top(call->calling) && 
		    (call->calling->nr == PLUS)// ||call->calling->nr == SELECT)
		    && is_top(call->calling->calling) ? MIXEDVAR : ANYPARAM;
		} else {
		  while (!is_top(call) && (call->nr ==PLUS)//||call->nr==SELECT)
			 )
		    call = call->calling;	
		  if (call == NULL) BUG;
		  sorts[*NAs] = is_top(call)
		    ? (next->nr == NUGGET ? NUGGETVAR : VARPARAM)
		    : ANYPARAM; 
		}
		//printf("bsorts %s %d %d %d %s\n", NAME(cov), sorts[*NAs], no_variance, is_top(call), NAME(call));
		if (no_variance && sorts[*NAs] <= SIGNEDSDPARAM)
		  sorts[*NAs] = ANYPARAM;
		//printf("sorts %s %d %d\n", NAME(cov),sorts[*NAs],no_variance);
		sprintf(names[*NAs], 
			sorts[*NAs] <= SIGNEDSDPARAM || sorts[*NAs] ==NUGGETVAR 
			|| sorts[*NAs]==MIXEDVAR ? "%s.var" : "%s.vx",
			shortD); // for R level only
	    } else {
		if (i==DSCALE) {
		  //   lastmem = mem[*NAs]; // used for all diagonal elements of
		  //                          aniso
		  sorts[*NAs] = SCALEPARAM;	
		  sprintf(names[*NAs], "%s.s", shortD);// for R level only
		} else {
		  assert(i == DANISO);
		  // no lastmem here !
		
		  if (cov->q != NULL && cov->q[0] != 0.0) {
		    //		    print("sdfasdf\n");
		    SERR("naturalscaling only allowed for isotropic setting");
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
	    } else { // not $
	      //printf("entering non $\n");
		// standard setting !
	      if (ParamIsTrend(cov, i)) {
		// || cov->nr==MLEMIXEDEFFECT)
		assert(type[i] == REALSXP);	      
		if (r == 0 && c==0) {
		  int rr, 
		    rowcol = nrow[i] * ncol[i];
		  for (rr=0; r<rowcol; r++) {
		    double vv = P(i)[rr];
		    if (!ISNAN(vv) && !ISNA(vv)) SERR("in case of trend parameters either none or all values must be 'NA'");
		  }
	      }
		
		if (excludetrends) continue;
		sorts[*NAs] = TRENDPARAM;
	      } else {
                sorts[*NAs] = SortOf(cov, i, r, c);  // ANYPARAM; 
                if (sorts[*NAs] == IGNOREPARAM || 
		    sorts[*NAs] == DONOTRETURNPARAM) continue;
		if (sorts[*NAs] == FORBIDDENPARAM) 
		  SERR1("The parameter '%s' may not be estimated", KNAME(i));
 	
		if (type[i] == INTSXP)
		  sorts[*NAs] = INTEGERPARAM; 
		else {
		  // r)ow, c)olumn; i:kappa, cov->nr, C		  
		  //		print("%s k=%d no_var=%d, %d \n", C->name, i,
		  //		       no_variance, sorts[*NAs]);  
		  if (no_variance && sorts[*NAs] <= SIGNEDSDPARAM)
		    sorts[*NAs] = ANYPARAM;
		}
	      }
	      char kappashort[255];
	      strcopyN(kappashort, CC->kappanames[i], SHORTlen);
	      sprintf(names[*NAs], "%s.%s", shortname,  kappashort);
	      
	      // printf("NAs %d %s\n", *NAs, names[*NAs]);
	      
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
	    //	  print("added %d %ld\n", *NAs, mem[*NAs]);
	    //printf("mem: %d %s %s\n", *NAs, NAME(cov), KNAME(i));
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
    } // sub == NULL || isRandom(sub); KEIN ELSE !!

    if (sub != NULL) {
      if (printing > 0) {
	leer(depth); PRINTF("%s\n", C->kappanames[i]);
      }
      if (printing > 0) {
	leer(depth);
	PRINTF("%s = ", C->subnames[i]);
      }
 
      int err = GetNAPosition(sub, NAs, mem,// elmnts, mem_elmnts,
			      names, sorts, vdim1, vdim2, isnan, bayesian,
			      covModels,
			      covzaehler, allowforintegerNA, 
			      SHORTlen, printing, depth, 
			      SortOf(cov, -i-1, 0, 0) != VARPARAM,
			      excludetrends);      
      if (err != NOERROR) return err;
      
      continue;
    }
  } // i
 
  for (i=0; i<MAXSUB; i++) {
    sub = cov->sub[i];

  
    if (sub != NULL) {
      // printf("cur %s sub=%s\n", NAME(cov), NAME(sub));
     if (printing > 0) {
	leer(depth);
	PRINTF("%s = ", C->subnames[i]);
      }
     int err = GetNAPosition(sub, NAs, mem, //elmnts, mem_elmnts,
			     names, sorts, vdim1, vdim2, isnan, bayesian,
			     covModels,
			     covzaehler, allowforintegerNA, 
			     SHORTlen, printing, depth, 
			     false, excludetrends);
        if (err != NOERROR) return err;
  
    }
  }

  return NOERROR;

}


/* wird (derzeit) nicht verwenden !! */
SEXP GetNAPositions(SEXP model_reg, SEXP model, SEXP spatialdim, SEXP Time,
		    SEXP xdimOZ, SEXP integerNA, SEXP Print) {
  int i, NAs, //elmnts, 
    vdim1[MAX_NA], vdim2[MAX_NA], covzaehler[MAXNRCOVFCTS];
  naptr_type mem;
  //  elptr_type mem_elmnts;
  covptr_type covModels;
  sortsofparam sorts[MAX_NA];
  bool isnan[MAX_NA], bayesian[MAX_NA];
  NAname_type names;
  currentRegister = INTEGER(model_reg)[0];

  bool skipchecks = GLOBAL.general.skipchecks;
  GLOBAL.general.skipchecks = true;
  
  CheckModelInternal(model, ZERO, ZERO, ZERO, INTEGER(spatialdim)[0], 
		     INTEGER(xdimOZ)[0], 1, 1, false, false, 
		     LOGICAL(Time)[0], 
		     R_NilValue,
		     KEY + currentRegister);  
  sprintf(ERROR_LOC, "getting positions with NA: ");

  // print("global=%d\n", GLOBAL.fit.lengthshortname);
  
  GLOBAL.general.skipchecks = skipchecks;
  NAs = 0;
  for (i=0; i<MAXNRCOVFCTS;  covzaehler[i++]=0);
  int err = GetNAPosition(KEY[currentRegister], &NAs, mem,//&elmnts, mem_elmnts,
			  names, sorts, vdim1, vdim2, isnan, bayesian, 
			  covModels,
			  covzaehler,
			  INTEGER(integerNA)[0],
			  GLOBAL.fit.lengthshortname, 
			  INTEGER(Print)[0], 
			  0, false, true);  
  if (err != NOERROR) XERR(err);
  sprintf(ERROR_LOC, "'%s' : ", NICK(KEY[currentRegister]));
  SEXP ans;
  PROTECT (ans =  allocVector(INTSXP, 1));
  INTEGER(ans)[0] = NAs;
  UNPROTECT(1);
  return ans;
}



int countnas(cov_model *cov, bool excludetrend, int level) {
    int i, r,  count, endfor,
      *nrow =  cov->nrow,
      *ncol = cov->ncol;
    cov_fct *C = CovList + cov->nr; // nicht gatternr
  SEXPTYPE *type = C->kappatype;

  count= 0;
  for (i=0; i<C->kappas; i++) {
    if (cov->kappasub[i] != NULL) {
      count += countnas(cov->kappasub[i], excludetrend, level + 1);
    }

    //    print(">> %s %d %d %d excl=%d %d\n", NAME(cov), level, isTrend(cov->typus),
    if (excludetrend &&
	(level == 0 || (level == 1 && cov->calling->nr == MULT)) &&
	ParamIsTrend(cov, i)
	//&& !PisNULL(i) // wichtig, da ja noch der Parameter durch ein Modell 
	)
      // gegeben sein koennte
      continue;
    int pt = SortOf(cov, i, 0, 0);
    
    if (nrow[i] == 0 || ncol[i] == 0 || pt == IGNOREPARAM
	|| pt == DONOTRETURNPARAM || pt == FORBIDDENPARAM)
      continue;
    endfor = nrow[i] * ncol[i];

    if (type[i] == REALSXP) { 

      //printf("%s %s %f\n", NAME(cov), KNAME(i), P0(i));

      double *p = P(i);
      for (r=0; r<endfor; r++) if (ISNAN(p[r])) count++;
    } else if (type[i] == INTSXP) {
      int *p = PINT(i);
      for (r=0; r<endfor; r++) if (p[r] == NA_INTEGER) count++;
    } else {
      continue; // no NAs allowed
    }
  } // for
  
  cov_model *sub;
  for (i=0; i<MAXSUB; i++) {
    sub = cov->sub[i];
    if (sub == NULL) continue; 
    count += countnas(sub, excludetrend, level + 1);
  }

  //  print("<< %d %d\n", level, count);

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

    if (cov->kappasub[i] != NULL) {
      Take21internal(cov->kappasub[i], cov_bound->kappasub[i],
		     bounds_pointer, NBOUNDS);
      continue;
    }

    if (ParamIsTrend(cov, i)) continue;
    
    int pt = SortOf(cov, i, 0, 0);
    if (C->kappatype[i] >= LISTOF || pt == IGNOREPARAM ||
	pt == DONOTRETURNPARAM || pt == FORBIDDENPARAM) continue;
    
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
	  if (!isDollar(cov) || 
	      i == DVAR ||
	      (i== DSCALE && cov->q == NULL) || // ! natscaling 
	      i == DANISO		
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

  if (nr[0] == nr[1]) ERR("do not use register 'model bounds'");
  
  NAOK_RANGE = true;
  if (LOGICAL(skipchecks)[0]) GLOBAL.general.skipchecks = true;
  for (m=1; m>=0; m--) { // 2er Schleife !!
    // assert(m==1);
    CheckModelInternal(models[m], ZERO, ZERO, ZERO, INTEGER(spatialdim)[0], 
		       INTEGER(xdimOZ)[0], 1, 1, false, false, 
		       LOGICAL(Time)[0], 
		       R_NilValue,
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
		 double *minpile, double *maxpile, int *NAs,
		 bool dosimulations) {
    /* 
       determine the ranges of the parameters to be estimated 
     */

  int i,  r,
    *nrow = cov->nrow,
    *ncol = cov->ncol;
  cov_fct 
    *C = CovList + cov->nr; // nicht gatternr
  SEXPTYPE *type = C->kappatype;
  double v, dmin, dmax;
  
  //    print("\ngetnaranges %s %d %d\n", NAME(cov), *NAs, C->kappas);
  //  PMI(cov);

  for (i=0; i<C->kappas; i++) {
    cov_model *sub = cov->kappasub[i];

    int end = nrow[i] * ncol[i];    
    //   print("%s %d %d \n", KNAME(i), end, sub == NULL);
  
    if (end > 0 && (sub == NULL || isRandom(sub))) {
      if (type[i] == REALSXP) {
	dmin = PARAM0(min, i);
	dmax = PARAM0(max, i);
      } else if (type[i] == INTSXP) {
	dmin = PARAM0INT(min, i) == NA_INTEGER 
	  ? RF_NA : (double) PARAM0INT(min, i);
	dmax = PARAM0INT(max, i) == NA_INTEGER 
	  ? RF_NA : (double) PARAM0INT(max, i);
      } else if (type[i] == LISTOF + REALSXP) {
	int store = GLOBAL.general.set;
	GLOBAL.general.set = 0;
	dmin = LPARAM0(min, i);
	dmax = LPARAM0(max, i);
	GLOBAL.general.set = store;
      } else if (isRObject(type[i])) {
	dmin = 0.0;
	dmax = 0.0;
      } else {
	BUG;
	dmin = dmax = RF_NA;
      }

      if (sub != NULL && end == 1 && dosimulations) {// i.e. isRandom
	int simulations = 1000; 
	double rr,
	  minr = RF_INF,
	  maxr = RF_NEGINF;
	for (int k=0; k<simulations; k++) {
	  DORANDOM(sub, &rr);
	  if (minr > rr) minr = rr;
	  if (maxr < rr) maxr = rr;
	}
	// "printf("r = %f %f %f %f %s\n", minr, maxr, dmin,dmax, KNAME(i));
	if (minr > dmin) dmin = minr;
	if (maxr < dmax) dmax = maxr;
      }
      
      int pt = SortOf(cov, i, 0, 0);
      if (pt == IGNOREPARAM || pt == DONOTRETURNPARAM || pt == FORBIDDENPARAM
	  || cov->nr==MIXEDEFFECT
	  || ParamIsTrend(cov, i))
	continue;
      
      for (r=0; r<end; r++) {
	v = RF_NA;
	if (type[i] == REALSXP) {
	  v = P(i)[r];
	} else if (type[i] == INTSXP) {
	  v = PINT(i)[r] == NA_INTEGER ? RF_NA : (double) PINT(i)[r];
	} else if (type[i] == LISTOF + REALSXP) {
	  break; //continue;  // !!!!!!!!!!!
	} else if (isRObject(type[i])) {
	  break; //  continue;  // !!!!!!!!!!!
	} else BUG;
	
	if (ISNAN(v)) {
	  if (isDollar(cov)) {
	    assert(i!=DAUSER && i!=DPROJ);
	  }
	  minpile[*NAs] = dmin;
	  maxpile[*NAs] = dmax;
	  //	    print("%s %d %f %f\n", C->kappanames[i],r, dmin, dmax);
	  (*NAs)++;
	} // isna
      } // r
    } // sub == NULL | isRandom
    
    // print("end %s %d\n", KNAME(i), sub == NULL);
    if (sub != NULL) {
      //  PMI(cov->kappasub[i]);
      GetNARanges(cov->kappasub[i], 
		  min->kappasub[i], max->kappasub[i], 
		  minpile, maxpile, 
		  NAs, dosimulations);
    }
  } // kappas
 

  for (i=0; i<MAXSUB; i++) {
    // printf("%s sub==NULL %d, %s\n", NAME(cov), cov->sub[i] == NULL, 
    //	   cov->sub[i] == NULL ?  "xxxx" : NAME(cov->sub[i]));
    if (cov->sub[i] != NULL) {
      GetNARanges(cov->sub[i], 
		  min->sub[i], max->sub[i], 
		  minpile, maxpile, 
		  NAs, dosimulations);
    }
  }

  //  print("end getnaranges %d\n", *NAs);

}



int check_recursive_range(cov_model *cov, bool NAOK) {
    /*
      called by "Set And GetModelInfo"
     */
  int i, err,
    kappa = CovList[cov->nr].kappas;
  sprintf(ERROR_LOC, "'%s' : ", NICK(cov));
  if ((err = check_within_range(cov, NAOK)) != NOERROR) return err;
   for (i=0; i<kappa; i++) 
    if (cov->kappasub[i] != NULL &&
	(err = check_recursive_range(cov->kappasub[i], NAOK)) != NOERROR)
      return err;

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL && 
	(err = check_recursive_range(cov->sub[i], NAOK)) != NOERROR)
	return err;
  }
  return NOERROR;
}


sortsofeffect getTrendEffect(cov_model *cov) {
  int nkappa = CovList[cov->nr].kappas;
  int j;
  for (j=0; j<nkappa; j++) {
    //    printf("(ParamIsTrend(cov, j =%d %d %d\n", ParamIsTrend(cov, j),  SortOf(cov, j, 0, 0), isTrend(cov->typus));
    if (ParamIsTrend(cov, j)) {
      if (!PisNULL(j)) {
	return ISNA(P0(j)) || ISNAN(P0(j)) ? FixedTrendEffect : DetTrendEffect;
      } else {	  
	if (cov->kappasub[j] == NULL) {	  
	  return DetTrendEffect; // 30.4.15 vorher FixedTrendEffect
	}
	else if (isRandom(cov->kappasub[j])) return RandomEffect;
	else if (cov->nr == TREND && j == TREND_MEAN) return DetTrendEffect;
	else ERR("model too complex");
      }
    }
  }
  //  printf("here A \n");
  return DetTrendEffect;
}


int CheckEffect(cov_model *cov) {
  int i,j,end, nr;    
  bool na_var = false;
  double *p;
  cov_model *next;
  if (cov->nr == MIXEDEFFECT) {    
    BUG;
    // cov->nr = MLEMIXEDEFFECT;  
    if (cov->nsub == 0) {
      return (cov->nrow[MIXED_BETA] > 0 && ISNAN(P0(MIXED_BETA)))
	? FixedEffect : DetTrendEffect;
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
	      return na_var ? SpVarEffect : SpaceEffect;
	      // NA in e.g. scale for constant matrix -- does not make sense
	    }
	  }
	}
      }
      next = next->sub[0]; // ok
    }
    
    int xx = na_var ? SpVarEffect : SpaceEffect;
    if (xx == SpVarEffect || xx == SpaceEffect) {
      BUG; //APMI(next);
    }
    
    
    return na_var ? SpVarEffect : SpaceEffect;
  }
  
  if (cov->nr == TREND) {
    int trend, param,
      effect = effect_error;
    bool isna ;

#define TREND_END 999    
    for (trend = TREND_MEAN; trend != TREND_END;
	 trend = trend == TREND_MEAN ? TREND_LINEAR : TREND_END) {
      if ((nr = cov->nrow[trend] * cov->ncol[trend]) > 0) {
	p = P(trend);
	isna = ISNAN(p[0]);
	if ((effect != effect_error) && ((effect==FixedTrendEffect) xor isna))
	  SERR1("do not mix deterministic effect with fixed effects in '%s'", 
		NICK(cov));
	for (i = 1; i<nr; i++) {
	  if (ISNAN(p[i]) xor isna) 
	    SERR("mu and linear trend:  all coefficient must be deterministic or all must be estimated");
	}
	effect = isna ? FixedTrendEffect : DetTrendEffect;
      } else if (cov->kappasub[trend] != NULL) effect = DetTrendEffect;
    }
    
    for (trend = TREND_POLY, param = TREND_PARAM_POLY; trend != TREND_END;
	 trend = trend == TREND_POLY ? TREND_FCT : TREND_END,
	   param = TREND_PARAM_FCT){
      if (cov->nrow[trend] > 0) {
	if (effect != effect_error)
	  SERR("polynomials and free functions in trend may not be mixed with other trend definitions. Please use a sum of trends.");
	
	nr = cov->nrow[param] * cov->ncol[param];
	if (nr > 0) {
	  p = P(param);
	  isna = ISNAN(p[0]);
	  for (i = 1; i<nr; i++) {
	    if (ISNAN(p[i]) xor isna) 
	      SERR("the coefficients in trend must be all deterministic or all coefficient are estimated");
	  }
	  effect = isna ? FixedTrendEffect : DetTrendEffect;
	} else effect = FixedTrendEffect;
      }
    }
    
    return effect;
  }

  // printf("here %s %d %d\n", NAME(cov), isTrend(cov->typus), getTrendEffect(cov));
    
  if (isTrend(cov->typus)) {
    if (cov->nr == MULT) {
      sortsofeffect current = getTrendEffect(cov->sub[0]);

      // printf("current %d\n", current);
      for (i=1; i<cov->nsub; i++) {
	sortsofeffect dummy = getTrendEffect(cov->sub[i]);
	if (current != DetTrendEffect && dummy != DetTrendEffect)
	  ERR("trend parameter to be estimated given twice");
	if (current == DetTrendEffect) current = dummy;
      }
      if (current == effect_error) ERR("trend mismatch");
      return current;
    } else return getTrendEffect(cov);
  }
 
  return RemainingError;
}


int GetEffect(cov_model *cov,  likelihood_info *info) {
  if (isProcess(cov)) {
    assert(cov->key == NULL);
    assert(!PisNULL(GAUSS_BOXCOX));
    int nas = 0,
      total = cov->nrow[GAUSS_BOXCOX] * cov->ncol[GAUSS_BOXCOX];
    for (int i=0; i<total; i++)
      nas += (int) (ISNA(P(GAUSS_BOXCOX)[i]) || ISNAN(P(GAUSS_BOXCOX)[i]));
    if (nas > 0) {

      //printf("GetEffect %d %d\n", nas, info->neffect); BUG;

      info->nas[info->neffect] = nas;
      info->effect[info->neffect] = DataEffect;
      (info->neffect)++;
    }
    return GetEffect(cov->sub[0], info);
  }
  
  bool excludetrend = true,
    plus = cov->nr == PLUS; // || cov->nr == SELECT;
  int i,
    n = plus ? cov->nsub : 1;
  if (info->neffect >= MAX_LIN_COMP) ERR("too many linear components");
  for (i=0; i<n; i++) {
    cov_model *component = plus ? cov->sub[i] : cov;

    // if (component->nr == TREND && component->kappasub[TREND_MEAN] != NULL) {
    //      GetEffect(component->kappasub[TREND_MEAN], info->effect, info->nas, neffect);
    //      continue;
    //    } else
    if (component->nr == PLUS) {
      GetEffect(component, info);
      continue;
    }
 
    info->effect[info->neffect] = CheckEffect(component);
  
    //      printf("i=%d %d %d %d\n", i,n, info->neffect, info->effect[0]);
    //  pmi(component, 0);
 

    info->nas[info->neffect] = countnas(component, excludetrend, 0);


    // printf("info %d %d\n", info->neffect, info->nas[info->neffect]);

    if (info->effect[info->neffect] == effect_error) 
      SERR("scaling parameter appears with constant matrix in the mixed effect part");
    if (info->effect[info->neffect] > LastMixedEffect) {
      //    printf("varmodel = %d %d %d  %d\n", info->varmodel,model_undefined, 
      //   model_morethan1, info->effect);
      info->varmodel =
	info->varmodel == model_undefined ? info->neffect : model_morethan1;
      info->Var = component;
    } 
    (info->neffect)++;
  }
 
  return NOERROR;
}
  


naptr_type MEMORY[MODEL_MAX+1];
covptr_type MEM_COVMODELS[MODEL_MAX+1];
double *MEM_PT_VARIANCE[MODEL_MAX+1];
//static elptr_type MEMORY_ELMNTS[MODEL_MAX+1];
int MEM_NAS[MODEL_MAX+1];
//int MEM_ELMNTS[MODEL_MAX+1];


int SetAndGetModelInfo(cov_model *key, int shortlen, 
		       int allowforintegerNA, bool excludetrend, // IN
		       int newxdim,  int globvar, // IN 
		       likelihood_info *info){
			// OUT
		       //  bool *trans_inv, bool *isotropic, double **Matrix,
		       //int *neffect, int effect[MAXSUB],
		       //int *NAs, int nas[MAXSUB], NAname_type names) { 
  int i, rows,
    mem_vdim1[MAX_NA],
    mem_vdim2[MAX_NA],
    covzaehler[MAXNRCOVFCTS], 
    jump = -1,
    err = NOERROR;
  cov_model 
    *cov = !isInterface(key) ? key : key->key == NULL ? key->sub[0] : key->key,
    *sub = !isProcess(cov) ? cov : cov->key == NULL ? cov->sub[0] : cov->key,
    *min=NULL,  *max=NULL, 
    *pmin=NULL,  *pmax=NULL, 
    *openmin=NULL,  *openmax=NULL;
  sortsofparam mem_sorts[MAX_NA];
  bool mem_isnan[MAX_NA],
    mem_bayesian[MAX_NA];
  double mle_min[MAX_NA], mle_max[MAX_NA], mle_pmin[MAX_NA], mle_pmax[MAX_NA],
    mle_openmin[MAX_NA], mle_openmax[MAX_NA],
    *matrix = NULL;

  sprintf(ERROR_LOC, "checking model : ");
  if (PrevLoc(key)->timespacedim > MAXMLEDIM) 
    ERR("dimension too high in MLE");

  if (sub->pref[Nothing] == PREF_NONE) {
    err = ERRORINVALIDMODEL;
    goto ErrorHandling;
  }
  
  info->newxdim = newxdim;
  info->trans_inv = sub->domprev == XONLY && isNegDef(sub->typus);
  info->isotropic = isAnyIsotropic(sub->isoown);
  if (info->trans_inv && info->isotropic) {
      info->newxdim = 1;
//	removeOnly(KEY + modelnr); nicht wegnehmen
// da derzeit negative distanczen in GLOBAL.CovMatrixMult auftreten -- obsolete ?!
  }

  
  // GLOBAL.general.skipchecks = skipchecks;
  check_recursive_range(key, true);

  if ((err = get_ranges(cov, &min, &max, &pmin, &pmax, &openmin, &openmax))
      != NOERROR) goto ErrorHandling;

  //  MEM_ELMNTS[currentRegister] = 
  MEM_NAS[currentRegister] = 0;
  for (i=0; i<MAXNRCOVFCTS;  covzaehler[i++]=0);

  // !! ZWINGEND VOR  GetEffect DA BAYESSCHE VARIABLE AUF NA GESETZT WERDEN
  if ((err = GetNAPosition(cov, MEM_NAS + currentRegister, MEMORY[currentRegister], 
		//MEM_ELMNTS + currentRegister, MEMORY_ELMNTS[currentRegister],
		info->names, mem_sorts, mem_vdim1, mem_vdim2, mem_isnan,
		mem_bayesian,
		MEM_COVMODELS[currentRegister],
		covzaehler,
		allowforintegerNA,
		shortlen,
		(PL >= PL_COV_STRUCTURE) + (PL >= PL_DETAILS) + 
		           (PL >= PL_SUBDETAILS), 
			   0, false, excludetrend))
       != NOERROR) goto ErrorHandling;

  sprintf(ERROR_LOC, "'%s' : ", NICK(key));

  //  PMI(sub);
  
  info->NAs = 0; GetNARanges(cov, min, max, mle_min, mle_max, &(info->NAs),
			     false);
  //printf("A\n");
  //   PMI(cov);
  //   print("XXz modelnr=%d %d\n",info->NAs, MEM_NAS[currentRegister]);
  if (info->NAs != MEM_NAS[currentRegister]) BUG;
  info->NAs = 0; GetNARanges(cov, pmin, pmax, mle_pmin, mle_pmax, &(info->NAs),
			     true);
  //printf("AB\n");
  //print("XX-------- amodelnr=%d\n",info->NAs);
 if (info->NAs != MEM_NAS[currentRegister]) BUG;
  info->NAs = 0; GetNARanges(cov, openmin, openmax, mle_openmin, mle_openmax,
			     &(info->NAs), false);
  if (info->NAs != MEM_NAS[currentRegister]) BUG;
 // print("XX-------- amodelnr=%d %d\n",info->NAs, MEM_NAS[currentRegister]); BUG;

  assert(info->varmodel == model_undefined);
  assert(info->neffect == 0);
  // GetEffect ersetzt NA-variance durch 1.0 
  if ((err = GetEffect(cov, info)) != NOERROR) goto ErrorHandling;
  if (info->NAs != MEM_NAS[currentRegister]) BUG;
  //  

  //  print("XX-------- amodelnr=%d\n",info->NAs); BUG;
 
  rows = MEM_NAS[currentRegister]; // bevor NA-variance durch 1.0 ersetzt wird
  info->globalvariance = false;
  info->pt_variance = NULL;
  assert(cov->nr == GAUSSPROC || 
	 cov->role == ROLE_COV);

  //printf(">>>>>>>>>> %s gloabvar=%d role=%d varmodel=%d\n",
  //	 NAME(cov), 	 cov->role,  globvar,  info->varmodel);

  if ((cov->nr == GAUSSPROC || cov->role == ROLE_COV) && 
      (globvar == NA_INTEGER || globvar) &&
      info->varmodel != model_undefined && info->varmodel != model_morethan1){
    // does a global variance exist?
     
    if (isDollar(info->Var) && !PARAMisNULL(info->Var, DVAR) &&
	info->Var->ncol[DVAR]==1 && info->Var->nrow[DVAR]==1) {
      double var = PARAM0(info->Var, DVAR);
       if ((info->globalvariance = ISNA(var) || ISNAN(var)) &&
	  info->Var->kappasub[DVAR] == NULL) {	
	info->pt_variance = PARAM(info->Var, DVAR);
	assert(info->pt_variance != NULL);
	*(info->pt_variance) = 1.0;
	(info->NAs)--;
	info->nas[info->varmodel]--;
	for (i=0; i<rows; i++) {
	  if (MEMORY[currentRegister][i] == info->pt_variance) {
	    jump = i;
	    for (int j=i+1; j<rows; j++) {
	      MEMORY[currentRegister][j-1] = MEMORY[currentRegister][j];
	      strcpy(info->names[j-1],info->names[j]);
	    }
	    break;
	  }
	}
	if (i >= rows) BUG;
	MEM_NAS[currentRegister]--;       
	rows--;
      }
    }
    info->globalvariance |= (globvar == true);    
  }
  if (info->NAs != MEM_NAS[currentRegister]) BUG;
  MEM_PT_VARIANCE[currentRegister] = info->pt_variance;

  //printf("globalvar = %d\n", info->globalvariance);
  // assert(info->globalvariance);
  //  assert(info->pt_variance != NULL);

  
  //PMI(cov);
  //print("XX-------- amodelnr=%d %d %d %d %d\n",info->NAs, MEM_NAS[currentRegister], info->globalvariance, globvar,  info->varmodel); BUG;

  //  pmi(cov, 1);
  //  printf("%d\n", info->NAs);
  

  if (rows == 0) {
    info->Matrix = NULL;
  } else {
    info->Matrix = 
      (double*) MALLOC(rows * MINMAX_ENTRIES * sizeof(double));
    matrix = info->Matrix;
    for (int k=i=0; i<rows; i++, k++) {
      if (i == jump) k++;
      int j = i - rows;
      // printf("i=%d j=%d %f %f %d\n", i,j,  mle_pmin[i], mle_pmax[i],  mem_sorts[i]);
      matrix[j + MINMAX_PMIN * rows ] = mle_pmin[k];
      matrix[j + MINMAX_PMAX * rows] = mle_pmax[k];
      matrix[j + MINMAX_TYPE * rows] = mem_sorts[k];
      matrix[j + MINMAX_NAN * rows] = (int) mem_isnan[k];
      matrix[j + MINMAX_MIN * rows] = mle_min[k];
      matrix[j + MINMAX_MAX * rows] = mle_max[k];
      matrix[j + MINMAX_OMIN * rows] = mle_openmin[k];
      matrix[j + MINMAX_OMAX * rows] = mle_openmax[k];
      matrix[j + MINMAX_ROWS * rows] = mem_vdim1[k];
      matrix[j + MINMAX_COLS * rows] = mem_vdim2[k];
      matrix[j + MINMAX_BAYES * rows] = mem_bayesian[k];
    }
  }

 ErrorHandling:
  if (min != NULL) COV_DELETE(&min);
  if (max != NULL) COV_DELETE(&max);
  if (pmin != NULL) COV_DELETE(&pmin);
  if (pmax != NULL) COV_DELETE(&pmax);
  if (openmin != NULL) COV_DELETE(&openmin);
  if (openmax != NULL) COV_DELETE(&openmax);
  
  return err;
}


SEXP SetAndGetModelInfo(SEXP model_reg, SEXP model, SEXP x, 
			int spatialdim, 
			bool distances, int lx, int ly, bool Time, int xdimOZ,
			int shortlen, int allowforintegerNA, 
			bool excludetrend) {
  // ex MLEGetNAPos
    /*
      basic set up of covariance model, determination of the NA's
      and the corresponding ranges.
     */
//  bool skipchecks = GLOBAL.general.skipchecks;
  int i, // NAs,
    //   unp = 0,
    //  effect[MAXSUB], nas[MAXSUB],
    //*pnas = nas,
    // neffect = 0,
    err = NOERROR;
  // double *Matrix = NULL;
  const char *colnames[MINMAX_ENTRIES] =
    {"pmin", "pmax", "type", "NAN", "min", "max", "omin", "omax",
     "col", "row", "bayes"};
  SEXP 
   matrix, nameAns, nameMatrix, RownameMatrix, 
   ans=R_NilValue;
  bool 
    do_not_del_info, isList_x = TYPEOF(x) == VECSXP;
  likelihood_info local,
    *info = NULL;
 //  isListofList_x = isList_x;
 
  if (isList_x) {
    if (length(x) == 0) BUG;
   //  isListofList_x = TYPEOF(VECTOR_ELT(x, 0)) == VECSXP;
  }

  //printf("init %d\n", xdimOZ);

  currentRegister = INTEGER(model_reg)[0];
  NAOK_RANGE = true;
  CheckModelInternal(model, 
		     length(x) == 0 
		     ? ZERO 
		     : isList_x ? NULL : REAL(x),
		     length(x) == 0 ? ZERO : NULL, 
		     length(x) == 0 ? ZERO : NULL, 
		     spatialdim, xdimOZ,
		     lx, ly, // lx, ly
		     false, //grid
		     distances, //distances
		     Time, 
		     isList_x ? x : R_NilValue,
		     KEY + currentRegister);  
  NAOK_RANGE = false;  

  cov_model *cov = KEY[currentRegister],
    *Likeli = cov;
  if (isInterface(Likeli)) {
    cov_model *process = cov->key != NULL ? Likeli->key : Likeli->sub[0];
    // printf("%s %d %s %d\n", NAME(Likeli), Likeli->Slikelihood == NULL, NAME(process), isProcess(process));
    if (Likeli->Slikelihood == NULL && isProcess(process))
      Likeli = process;
  }


  
  // PMI(Likeli);

  if ((do_not_del_info = Likeli->Slikelihood != NULL)) {    
    info = &(Likeli->Slikelihood->info);
    // pmi(cov,1);   BUG;
  } else {
    info = &local;
    likelihood_info_NULL(info);
    err = SetAndGetModelInfo(cov, shortlen, allowforintegerNA, excludetrend, 
			     xdimOZ, GLOBAL.fit.estimate_variance,
			     info);
    if (err != NOERROR) goto ErrorHandling;
    //    pmi(Likeli, 1);    BUG;
  }

  int rows;
  rows = MEM_NAS[currentRegister];
  PROTECT(matrix =  allocMatrix(REALSXP, rows, MINMAX_ENTRIES));
  PROTECT(RownameMatrix = allocVector(STRSXP, rows));
  PROTECT(nameMatrix = allocVector(VECSXP, 2));
#define total 7
  PROTECT(ans = allocVector(VECSXP, total));
  PROTECT(nameAns = allocVector(STRSXP, total));

  //printf("%d\n", info->Matrix != NULL);

  if (rows > 0) {
    MEMCOPY(REAL(matrix), info->Matrix, rows * MINMAX_ENTRIES * sizeof(double));
    for (i=0; i<rows; i++) {
      //printf("%d rows %d\n", i, rows);
      //printf(">>%d %s\n", i, info->names[i]);
      SET_STRING_ELT(RownameMatrix, i, mkChar(info->names[i]));
    }
  }

  SET_VECTOR_ELT(nameMatrix, 0, RownameMatrix);
  SET_VECTOR_ELT(nameMatrix, 1, Char(colnames, MINMAX_ENTRIES));
  setAttrib(matrix, R_DimNamesSymbol, nameMatrix);
   
  i = 0;
  SET_STRING_ELT(nameAns, i, mkChar("minmax"));
  SET_VECTOR_ELT(ans, i++, matrix);
  SET_STRING_ELT(nameAns, i, mkChar("trans.inv")); // !!
  SET_VECTOR_ELT(ans, i++, ScalarLogical(info->trans_inv));  
  SET_STRING_ELT(nameAns, i, mkChar("isotropic"));
  SET_VECTOR_ELT(ans, i++, ScalarLogical(info->isotropic));  
  SET_STRING_ELT(nameAns, i, mkChar("effect"));
  SET_VECTOR_ELT(ans, i++, Int(info->effect, info->neffect, MAXINT));  
  SET_STRING_ELT(nameAns, i, mkChar("NAs"));
  SET_VECTOR_ELT(ans, i++, Int(info->nas, info->neffect, MAXINT));  
  SET_STRING_ELT(nameAns, i, mkChar("xdimOZ"));
  SET_VECTOR_ELT(ans, i++, ScalarInteger(info->newxdim));  
  SET_STRING_ELT(nameAns, i, mkChar("matrix.indep.of.x"));
  SET_VECTOR_ELT(ans, i++,
		 ScalarLogical(KEY[currentRegister]->matrix_indep_of_x));  
  setAttrib(ans, R_NamesSymbol, nameAns);
  assert(i==total);

 ErrorHandling:
  UNPROTECT(5);
  if (info != NULL && !do_not_del_info) likelihood_info_DELETE(info);

  return ans;

}


 
 
SEXP SetAndGetModelLikeli(SEXP model_reg, SEXP model, SEXP x) {
  return SetAndGetModelInfo(model_reg, model, x, 			    
			    0,  false, 0, 0, false, 0,// only dummy variables
			    GLOBAL.fit.lengthshortname, LIKELI_NA_INTEGER,
			    LIKELI_EXCLUDE_TREND);
}



SEXP SetAndGetModelInfo(SEXP model_reg, SEXP model, SEXP spatialdim, 
			SEXP distances, SEXP ygiven, SEXP Time,
			SEXP xdimOZ, SEXP shortlen, SEXP allowforintegerNA, 
			SEXP excludetrend) {
  return SetAndGetModelInfo(model_reg, model, R_NilValue, 
			    INTEGER(spatialdim)[0], LOGICAL(distances)[0], 1,
			    LOGICAL(ygiven)[0] ? 1 : 0,
			    LOGICAL(Time)[0], INTEGER(xdimOZ)[0], 
			    INTEGER(shortlen)[0],
			    INTEGER(allowforintegerNA)[0],
			    LOGICAL(excludetrend)[0]);
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
  double *pt_variance = MEM_PT_VARIANCE[*reg];
  gen_NULL(&s);
  s.check = s.dosimulate = false;
  // set ordinary parameters for all (sub)models

  currentRegister = *reg;
   for (un=i=0; i<NAs; i++) {
 //      print("reg=%d i=%d %d %ld %f %ld\n", *reg, i, NAs, MEMORY[*reg][i], values[un], pt_variance);
    //     print("mem=%ld %f\n", (long int) MEMORY[*reg][i], values[un]) ;  
  if (MEMORY[*reg][i] == pt_variance) BUG; //continue;
    MEMORY[*reg][i][0] = values[un++];
  }
  
  if (init)
    for (i=0; i<NAs; i++) {
      cov = MEM_COVMODELS[*reg][i];
      C = CovList + cov->nr;       
      if (i==0 || cov != MEM_COVMODELS[*reg][i-1]) {
	if (!isDummyInit(C->Init)) {
	  //	  print("i=%d %s initalising (%f)\n", i, NAME(cov),   MEMORY[*reg][i][0]);
	  
	  C->Init(cov, &s);
	  //  print("i=%d %s initalised (%f)\n", i, NICK(MEM_COVMODELS[*reg][i]), MEMORY[*reg][i][0]);
	}
      }
      //  
    }
  
  //printf("done\n");

  //  int one = 1;  setListElements(reg, &one, &one, &one);
}

void PutValuesAtNA(int *reg, double *values){
  PutValuesAtNAintern(reg, values, true);
}

void PutValuesAtNAnoInit(int *reg, double *values){
  PutValuesAtNAintern(reg, values, false);
}


/*
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
*/
