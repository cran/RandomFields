

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
#include <assert.h>
#include <string.h>
#include "RF.h"
#include "primitive.h"
// #include "Covariance.h"

void DeleteKey(int *keyNr) {  
  if (PL>=4) { 
    PRINTF("deleting stored parameters of register %d ...\n", *keyNr);
  }
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {return;}
  KEY_DELETE(KEY + *keyNr);
}

void DeleteAllKeys() {
  int i;
  for (i=0; i<MAXKEYS; i++) {DeleteKey(&i);}
  for (i=0; i<MODEL_MAX; i++) {
      if (STORED_MODEL[i]!=NULL)
	  COV_DELETE(STORED_MODEL + i);
  }
}

SEXP CheckModelUser(SEXP model, SEXP tsdim, SEXP xdim, SEXP stationary) {
    // MAXCOVDIM
  SEXP ans;
  PROTECT (ans = allocVector(INTSXP, 1));

  CheckModel(model, tsdim, xdim, stationary,  STORED_MODEL + MODEL_USER,
      MAXCOVDIM);

  if (GLOBAL.general.printlevel > 5) {
    PrintModelInfo(STORED_MODEL[MODEL_USER]); // OK
    // assert(false);
  }
  INTEGER(ans)[0] =  STORED_MODEL[MODEL_USER]->vdim;
  UNPROTECT(1);
  return ans; // islist
}

SEXP CheckModelIntern(SEXP model, SEXP tsdim, SEXP xdim, SEXP stationary) {
    // MAXSIMUDIM
  SEXP ans;
  PROTECT (ans = allocVector(INTSXP, 1));
  CheckModel(model, tsdim, xdim, stationary, STORED_MODEL + MODEL_INTERN,
      MAXSIMUDIM);
  if (GLOBAL.general.printlevel > 5) {
   PrintModelInfo(STORED_MODEL[MODEL_INTERN]);// OK
  }
  INTEGER(ans)[0] =  STORED_MODEL[MODEL_INTERN]->vdim;
  UNPROTECT(1);
  return ans;
}

 
SEXP CheckModelSimu(SEXP model, SEXP tsdim, SEXP xdim, SEXP stationary) {
  //assert(false);
  SEXP ans;
  PROTECT (ans = allocVector(INTSXP, 1));

// printf("aaaaaa\n");
  // PrintModelInfo(STORED_MODEL[MODEL_SIMU]);
  // UNPROTECT(1); return ans;

  CheckModel(model, tsdim, xdim, stationary, STORED_MODEL + MODEL_SIMU,
	     MAXSIMUDIM);

  if (GLOBAL.general.printlevel > 5) {
    PrintModelInfo(STORED_MODEL[MODEL_SIMU]);// OK
  }
  INTEGER(ans)[0] =  STORED_MODEL[MODEL_SIMU]->vdim;
  UNPROTECT(1);
  return ans;
}


SEXP SetParamPch(SEXP act, SEXP pch) {
  int action = INTEGER(act)[0];
  general_param *gp = &(GLOBAL.general);
  if (action) {
    gp->pch = ((char*) CHAR(STRING_ELT(pch, 0)))[0];
//    if ((strlen(*pch)>1) && (PL>0)) 
//	PRINTF("\n`pch' has more than one character -- first character taken only\n");
    return pch;
  } else {
    SEXP neupch;
    char p[2];
    p[0] = gp->pch;
    p[1] = '\0'; 
    PROTECT (neupch =  allocVector(STRSXP, 1));
    SET_STRING_ELT(neupch, 0, mkChar(p));
    UNPROTECT(1);
    return neupch;
  }
}

void SetParam(int *action, int *storing, int *printlevel, int *naturalscaling,
	      int *skipchecks, int *every) {
  general_param *gp = &(GLOBAL.general);
  if (*action) {
    gp->storing= (bool) *storing;
    PL = gp->printlevel = *printlevel;
    
//    printf("PL=%d\n", PL);

    NS = gp->naturalscaling=*naturalscaling;
    gp->skipchecks=*skipchecks;
    gp->every=*every;
  } else {
    *storing =  gp->storing;
    *printlevel = gp->printlevel;
    *naturalscaling = gp->naturalscaling;
    *skipchecks = gp->skipchecks;
    *every      = gp->every;
  }
}
void SetParamDecision(int *action, int *stationary_only, int *exactness)
{
  int NA=-1;
  decision_param *gp = &(GLOBAL.decision);
  if (*action) {
    gp->stationary_only = (*stationary_only==NA) ?
       DECISION_CASESPEC : *stationary_only ? DECISION_TRUE : DECISION_FALSE;
    gp->exactness = (*exactness==NA) ?
      DECISION_CASESPEC : *exactness ? DECISION_TRUE : DECISION_FALSE;
  } else {
    *stationary_only = gp->stationary_only==DECISION_TRUE ? true :
      gp->stationary_only==DECISION_FALSE ? false : NA;
    *exactness = gp->exactness==DECISION_TRUE ? true :
      gp->exactness==DECISION_FALSE ? false : NA;
  }
}


void SetParamCE(int *action, int *force, double *tolRe, double *tolIm,
		int *trials, 
		double *mmin, // mmin must be vector of length MAXCEDIM!!
		int *useprimes, int *strategy, double *maxmem,
		int *dependent, ce_param *CE, const char *name){
  int d;
  if (*action) {
    CE->force=(bool)*force;

    CE->tol_re=*tolRe;
    if (CE->tol_re>0) {
      CE->tol_re=0;
      if (PL>0)
	PRINTF("\nWARNING! %s.tol_re which has been positive is set to 0.\n",name);
    }

    CE->tol_im=*tolIm; 
    if (CE->tol_im<0) {
      CE->tol_im=0; 
      if (PL>0) 
	PRINTF("\nWARNING! %s.tol_im which has been negative is set 0.\n",name);
    }

    CE->trials=*trials; 
    if (CE->trials<1) {
      CE->trials=1;
      if (PL>0) 
	PRINTF("\nWARNING! %s.trials had been less than 1\n",name);
    }

    for (d=0; d<MAXCEDIM; d++) {
      CE->mmin[d] = mmin[d];
      if (CE->mmin[d]<0.0) {
	if ((CE->mmin[d]>-1.0)) {
	    CE->mmin[d] = -1.0;
	    if (PL>0) PRINTF("\nWARNING! %s.mmin[%d] set to -1.0.\n", name, d);
	}
      }
    }
 
    CE->useprimes=(bool)*useprimes;
    CE->dependent=(bool)*dependent;

    if (*strategy>LASTSTRATEGY) {  
      if (PL>0) 
	PRINTF("\nWARNING! %s_STRATEGY not set\n",name);
    } else CE->strategy=(char) *strategy;     

    CE->maxmem = *maxmem;
  } else {
    *force=CE->force;
    *tolRe=CE->tol_re;
    *tolIm=CE->tol_im;
    *trials=CE->trials;
    for (d=0; d<MAXCEDIM; d++) {mmin[d]=CE->mmin[d];}   
    *useprimes= CE->useprimes; //
    *dependent= (int) CE->dependent;
    *strategy = (int) CE->strategy;
    *maxmem = CE->maxmem; //
  }
}
  

#define MaxMaxInts 9
void GetMaxDims(int *maxints) {  *maxints = MaxMaxInts; }

void GetrfParameters(int *covmaxchar, int *methodmaxchar, 
		     int *distrmaxchar,
		     int *covnr, int *methodnr, int *distrnr,
		     int *maxdim,
		     int *maxmodels) {
  // if (currentNrCov==-1) InitModelList();
  int i=0;
  *covmaxchar=MAXCHAR;
  *methodmaxchar=METHODMAXCHAR;
  *distrmaxchar=DISTRMAXCHAR;
  *covnr=currentNrCov;
  *methodnr=(int) Forbidden;
  *distrnr=DISTRNR;
  maxdim[i++]=MAXCOVDIM;
  maxdim[i++]=MAXMLEDIM;
  maxdim[i++]=MAXSIMUDIM;
  maxdim[i++]=MAXCEDIM;
  maxdim[i++]=MAXTBMSPDIM;
  maxdim[i++]=MAXMPPDIM;
  maxdim[i++]=MAXHYPERDIM;
  maxdim[i++]=MAXNUGGDIM;
  maxdim[i++]=MAXVARIODIM;
  assert(i == MaxMaxInts);
  *maxmodels=MAXKEYS;
}

void GetrfParametersI(int *covmaxchar, int *methodmaxchar, int *distrmaxchar,
		      int *covnr, int *methodnr, int *distrnr,
		      int *maxdim, int *maxmodels){
  if (currentNrCov==-1) InitModelList();
  GetrfParameters(covmaxchar, methodmaxchar, distrmaxchar,
		  covnr, methodnr, distrnr, maxdim, maxmodels);
}


void GetDistrName(int *nr,char **name){
  if ((*nr<0) ||(*nr>=DISTRNR)) {
    strcopyN(*name,"",DISTRMAXCHAR); 
    return;
  }
  strcopyN(*name, DISTRNAMES[*nr], DISTRMAXCHAR);
}

void GetDistrNr(char **name, int *n, int *nr) 
{
  unsigned int ln, v, nn;
  nn = (unsigned int) *n;
  // == -1 if no matching function is found
  // == -2 if multiple matching fnts are found, without one matching exactly
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)

  for (v=0; v<nn; v++) {
    nr[v]=0;
    ln=strlen(name[v]);
    while ( (nr[v] < DISTRNR) && strncmp(name[v], DISTRNAMES[nr[v]], ln)) {
      (nr[v])++;
    }
    if (nr[v]<DISTRNR) { 
      // a matching function is found. Are there other functions that match?
      int j; 
      bool exactmatching,multiplematching;
      exactmatching=(ln==strlen(DISTRNAMES[nr[v]]));
      multiplematching=false;
      j=nr[v]+1; // if two or more distributions have the same name 
      //            the last one is taken 
      while (j<DISTRNR) {
	while ( (j<DISTRNR) && strncmp(name[v],DISTRNAMES[j],ln)) {j++;}
	if (j<DISTRNR) {
	  if (ln==strlen(DISTRNAMES[j])) {nr[v]=j; exactmatching=true;} 
	  else {multiplematching=true;}
	}
	j++;
      } 
      if (!exactmatching && multiplematching) {nr[v]=-2;}
    } else nr[v]=-1;
  }
}


void GetModelName(int *nr,char **name){
  if (currentNrCov==-1) InitModelList();
  if ((*nr<0) ||(*nr>=currentNrCov)) {
    strcopyN(*name,"", MAXCHAR); 
    return;
  }
  strcopyN(*name, CovList[*nr].name, MAXCHAR);
}

void GetNrParameters(int *covnr, int* kappas) {
  if (currentNrCov==-1) InitModelList();
  if (*covnr<0 || *covnr>=currentNrCov) {*kappas=-999;}
  else *kappas = CovList[*covnr].kappas;
}


void GetModelNr(char **name, int *nr) {
  *nr = getmodelnr(*name);
}


void GetMethodNr(char **name, int *nr) {
  // == -1 if no matching method is found
  // == -2 if multiple matching methods are found, without one matching exactly
  *nr = Match(*name, METHODNAMES, (int) Forbidden);
}

void GetMethodName(int *nr, char **name) 
{   
  if ((*nr<0) ||(*nr>=(int) Forbidden)) {
    strcopyN(*name,"",METHODMAXCHAR); 
    return;
  }
  strcopyN(*name, METHODNAMES[*nr], METHODMAXCHAR);
}

void PrintMethods()
{ 
  int i;
  PRINTF("\n\n  Methods for generating Gaussian random fields\n  =============================================\n\n");  
  for (i=0;i<(int) Nothing;i++) { PRINTF("  * %s\n",METHODNAMES[i]); }
  PRINTF("\n\n  Methods for non-Gaussian random fields\n  ======================================\n");
  for (i=1+(int) Nothing;i<(int) Forbidden; i++) { 
    PRINTF("  * %s\n",METHODNAMES[i]); }
  PRINTF("\n");
}


SEXP GetRange(){
  assert(false); // muss neu geschrieben werden, als SEXP

/*
  // never change double without crosschecking with fcts in RFCovFcts.cc!
  // index is increased by one except index is the largest value possible
  cov_fct *C;
  cov_model cov;
  range_type *r;
  int i,j, kappas;

  if (currentNrCov==-1) InitModelList();
  if ((*nr<0) || (*nr>=currentNrCov)) goto ErrorHandling;
  C = CovList + *nr;
  kappas = CovList[*nr].kappas;
  if (*lparam > kappas || *lrange != kappas * 4) goto ErrorHandling;
  if (*index < 0) {
    range_default(&getrange);    
    for (i=0; i<lparam; i++) {
      cov.p[i] = (double *) calloc(sizeof(double), 1);
      cov.p[i][0] = param[i];
    }
    range_default(&getrange);
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
  for (i=0; i<*lrange; i++) range[i]=RF_NAN;
 *index = -100;
 return;
*/
}


void PMLheader(char* firstcolumn) {
  char header1[]="#    cir cut int TBM TBM spe dir seq Mar ave coi hyp\n";
  char header2[]="p    cul off rin  2   3  ctr ect uen kof rag ns  erp\n";
  PRINTF(firstcolumn, ""); PRINTF("%4s", ""); PRINTF(header1);  
  PRINTF(firstcolumn, ""); PRINTF("%4s", ""); PRINTF(header2);
}

void PrintModelList(int *intern, int *operat)
{
    int i, k, last_method, m, OP;
//char header[]="circ cut intr TBM2 TBM3 spec dir seq Mark ave add hyp part\n";
  char coded[6][2]={"-", "X", "+", "N", "H", "S"};
//  char typenames[4][2]={"i", "f", "s", "n"};
  char specialnames[6][2]={"-", "o", "n", "f", "i", " "};
  char firstcolumn[20], name[MAXCHAR];
  int maxchar=10; // <= MAXCHAR=14
  int op[MAXNRCOVFCTS], normalmix[MAXNRCOVFCTS], finite[MAXNRCOVFCTS], 
    internal[MAXNRCOVFCTS], stat[MAXNRCOVFCTS], iso[MAXNRCOVFCTS],
    vdim[MAXNRCOVFCTS], maxdim[MAXNRCOVFCTS];
  
  last_method = (int) Nothing; // even not special method   
  if (currentNrCov==-1) {
     InitModelList(); 
   // if (PL>5)  PRINTF("List of covariance functions initiated.\n"); 
  }
  if (CovList==NULL) {PRINTF("There are no functions available!\n");} 
  else {
    cov_fct *C;

    GetAttr(op, normalmix, finite, internal, stat, iso, maxdim, vdim);
 
    sprintf(firstcolumn,"%%%ds ", -maxchar);
    PRINTF("\n\n");
    PRINTF("%20s      List of models\n", "");
    PRINTF("%20s      ==============\n", "");
    PRINTF("%10s[See also PrintMethodList for the names of the columns();\n", 
	   "");
    PRINTF("%10s use `operator=TRUE' to see all available models        ]\n", 
	   "");

    for (OP = 0; OP <= *operat; OP++) {
      C = CovList;
      PRINTF("\n\n");
      if (OP) {
	PRINTF("%4s Operators\n", "");
 	PRINTF("%4s =========\n\n", "");  
     } else {
	PRINTF("%4s Simple models\n", "");
	PRINTF("%4s =============\n\n", "");  
      }
       PMLheader(firstcolumn);

      for (k=1, i=0; i<currentNrCov; i++, C++) { 
	  if ((op[i] != OP) || (!*intern && internal[i])) continue;
	strcopyN(name, C->name, maxchar);
	if (strncmp(C->name, InternalName, strlen(InternalName)) ==0 &&
	    *intern < 2)  continue;
	PRINTF("%2d. ", k++);
	PRINTF(firstcolumn, name);
	PRINTF("%1d ", C->kappas);
//	PRINTF("%s", internal[i] ? specialnames[4] :
//		 (C->maxsub==0) ? specialnames[5] : specialnames[op[i]]);
	PRINTF("%s", specialnames[(C->stationary == STATIONARY && normalmix[i])
				  * 2 + finite[i] * 3]);	  
	PRINTF(" ");
	//             above works since normal mixture cannot have finite.range 
	assert(internal[i]==0 || internal[i]==1);
	
	for (m=(int) CircEmbed; m<last_method; m++)
	    if (m != Nugget)
		PRINTF("%3s%s", coded[(int) C->implemented[m]], 
		       " ");
	PRINTF("\n");
      }
    }
    
    PMLheader(firstcolumn);
    PRINTF("\n%4sLegend:","");
    PRINTF("\n%4s=======\n","");
    PRINTF("First row after number of parameters:\n");
//    PRINTF("'%s': isotropic model\n", typenames[0]);
//    PRINTF("'%s': fully symmetric spatio-temporal model\n", typenames[1]);
//    PRINTF("'%s': stationary model\n", typenames[2]);
//    PRINTF("'%s\n': non-stationary model\n", typenames[3]);
//    if (*intern!=0) PRINTF("'%s': internal\n", specialnames[4]);
//    PRINTF("'%s': operator\n", specialnames[1]);
//    PRINTF("'%s': primitive\n", specialnames[5]);

//    PRINTF("\nSecond row:\n");
    PRINTF("'%s': normal mixture model or allows for it if operator\n",
	   specialnames[2]);
    PRINTF("'%s': finite range or allows for it if submodels suit\n", 
	   specialnames[3]);

    PRINTF("\nAll other rows:\n");
    PRINTF("'%s': method not available\n", 
	   coded[0]);
    PRINTF("'%s': method available for at least some parameter values\n",
	   coded[1]);
    PRINTF("'%s': integral for the covariance is evaluated only numerically \n",  
	   coded[2]);
//    PRINTF("'%s': method is converted to the method `nugget'\n",   coded[3]);
//    PRINTF("'%s': only internally available within the method \n", 
//	   coded[4]);
//    PRINTF("'%s': available as submodel within the method\n", 
//	   coded[5]);
    PRINTF("\n");
  }
}

void PrintModelList() {
    int True = 1;
    PrintModelList(&True, &True);
}

void GetModelList(int* idx, int*internal) {
  int i, j, m;
  if (currentNrCov==-1) InitModelList(); 

  if (CovList==NULL) return;
  for (j=i=0; i<currentNrCov; i++) {
    if (!*internal && CovList[i].internal) continue;
    for (m=(int) CircEmbed; m<(int) Nothing; m++) {
      idx[j++] = CovList[i].implemented[m];
      //  printf("%d\n", j);
    }
  }
  return;
}


void StoreTrend(int *keyNr, int *modus, char **Trend, 
		int *lt, double *lambda, int *ll, int *err)
{
  unsigned long len;
  trend_type *trend;
  key_type *key;
  if (currentNrCov==-1) InitModelList(); assert(CovList!=NULL); 
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*err=ERRORREGISTER;goto ErrorHandling;}

  key = &(KEY[*keyNr]);
  TREND_DELETE(&(key->trend));
  TREND_NULL(&(key->trend));
  trend = &(key->trend);
  switch(*modus){
      case TREND_MEAN : //0
	trend->mean = *lambda;
	assert(*lt==0 && *ll==1);
	break;
      case TREND_PARAM_FCT : case TREND_FCT: // 3, 1
	trend->LinearTrend = (double *) malloc(len = sizeof(double) * *ll);
	memcpy(trend->LinearTrend, lambda, len);
	if (*modus==1) break; 
	trend->lLinTrend = *ll;
      case TREND_LINEAR : //2
	trend->TrendFunction = (char *) malloc(len = sizeof(char) * *lt);
	memcpy(trend->TrendFunction, *Trend, len);
	trend->lTrendFct = *lt;
	break; 
      default: 
	trend->TrendModus = -1; 
 	*err=ERRORUNSPECIFIED; 
	goto ErrorHandling;
  }
  trend->TrendModus = *modus;
  *err = 0;
  return;

 ErrorHandling:
  return;
}


void GetTrendLengths(int *keyNr, int *modus, int *lt, int *ll, int *lx)
{
  // currently unused, since AddTrend in RFsimu.cc is being programmed
  assert(false);
  trend_type *trend; 

  key_type *key;
  if (currentNrCov==-1) {*modus=-1;goto ErrorHandling;}
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*modus=-2; goto ErrorHandling;}
  key = &(KEY[*keyNr]);
  if (!key->simu.active) {*modus=-3; goto ErrorHandling;}
  trend = &(key->trend);
  *modus = trend->TrendModus;
  *lt = trend->lTrendFct;
  *ll = trend->lLinTrend;
  *lx = key->loc.totalpoints;
  return;
 ErrorHandling:
  return;
}



void GetTrend(int *keyNr, char **Trend, double *lambda, double *x, int *err)
{
  // currently unused, since AddTrend in RFsimu.cc is being programmed
  assert(false);
  trend_type *trend; 
  
  key_type *key;
  if (currentNrCov==-1) {*err=-1;goto ErrorHandling;}
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*err=-2; goto ErrorHandling;}
  key = &(KEY[*keyNr]);
  if (!key->simu.active) {*err=-3; goto ErrorHandling;}
  trend = &(key->trend);
  switch(trend->TrendModus){
  case 0 : break;
  case 3 : case 1 : 
    memcpy(lambda, trend->LinearTrend, trend->lLinTrend);
    if (trend->TrendModus==1) break; 
  case 2:
    strcpy(*Trend, trend->TrendFunction);
    break; 
  default: *err=-4; goto ErrorHandling;
  }
  *err = 0;
  return;

 ErrorHandling:
  return;
}





void GetKeyInfo(int *keyNr, int *total, int *lengths, int *spatialdim, 
		int *timespacedim, 
		int *grid, int *distr, int *maxdim, int *vdim)
{
    // check with DoSimulateRF and subsequently with extremes.cc, l.170
    // if any changings !
  int d;
  key_type *key;
  simu_type *simu;
  location_type *loc;

  if (*maxdim<MAXSIMUDIM) { *total=-3; *maxdim=MAXSIMUDIM; return;} 
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*total=-1; return;} 
  key = &(KEY[*keyNr]);  
  simu = &(key->simu);
  loc = &(key->loc);

//  PrintModelInfo(key->cov);

  if (!simu->active) {*total=-2;} 
  else { 
    *total=loc->totalpoints;
    *spatialdim = loc->spatialdim;
    *timespacedim = loc->timespacedim;
    for (d=0; d<loc->timespacedim; d++) lengths[d]=loc->length[d];
    *grid = (int) loc->grid;
    *distr = (int) simu->distribution; 
    *maxdim = MAXSIMUDIM;
    *vdim = key->cov->vdim;
  }
}


void GetAttr(int *op, int *normalmix, int *finiterange, int *internal, 
	     int *stat, int *iso, int *maxdim, int *vdim) {
#define MAXPN 10 /* only used for testing purposes */
  int i;
  cov_model cov;
  range_arraytype ra;
//  range_type *range;
  cov_fct *C = CovList;
  for (i=0; i<MAXPARAM; i++) {
    cov.p[i] = (double *) calloc(sizeof(double), MAXPN);
    cov.ncol[i] = cov.nrow[i] = 1;
    cov.tsdim = 1;
  }
  cov.vdim = 1;
  cov.nsub = 2;
 
  for (i=0; i<currentNrCov; i++, C++){      
    op[i] = (int) C->maxsub > 0;   
    range_default(&ra);
    C->range(&cov, &ra);
    normalmix[i] = ra.ranges[0].maxdim == INFDIM;
    maxdim[i] = ra.ranges[0].maxdim;
    finiterange[i] = ra.ranges[0].finiterange;
    
    internal[i] = C->internal;
    stat[i] = C->stationary;
    iso[i] = C->isotropy;
    vdim[i] = C->vdim;
  }
  for (i=0; i<MAXPARAM; i++) {
    free(cov.p[i]);
  }
}





SEXP IsStatAndIsoUser(SEXP removeGatter) {
    /* 
       called by Covariance in rf.R to check whether 
       distances should be passed resp. whether y may not be given.
     */
  SEXP ans, stationary, isotropy, nameAns;
  cov_model *cov = STORED_MODEL[MODEL_USER];
  int i;

  PROTECT(stationary=allocVector(LGLSXP, 1));
  LOGICAL(stationary)[0] = cov->statIn==STATIONARY || cov->statIn==VARIOGRAM;
  PROTECT(isotropy=allocVector(LGLSXP, 1));
  LOGICAL(isotropy)[0] = cov->isoIn==ISOTROPIC;

  if (cov->nr >= GATTER && cov->nr <= LASTGATTER) {
    cov_model *next = cov->sub[0];
    if (next->statIn==STATIONARY || next->statIn==VARIOGRAM) {
      LOGICAL(stationary)[0] = true;
      if (next->isoIn==ISOTROPIC) {
	LOGICAL(isotropy)[0] = true;
	if (LOGICAL(removeGatter)[0]) 
	  removeOnly(STORED_MODEL + MODEL_USER);
      }
    }
  }

  PROTECT(ans = allocVector(VECSXP, 3));
  PROTECT(nameAns = allocVector(STRSXP, 3));
  i = 0;
  SET_STRING_ELT(nameAns, i, mkChar("stationary"));
  SET_VECTOR_ELT(ans, i++, stationary);  
  SET_STRING_ELT(nameAns, i, mkChar("isotropy"));
  SET_VECTOR_ELT(ans, i++, isotropy);  
  setAttrib(ans, R_NamesSymbol, nameAns);

  UNPROTECT(4);
  return ans;
}

