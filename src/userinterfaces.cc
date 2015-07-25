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
#include <unistd.h>
 
#include <string.h>
#include "RF.h"
#include "primitive.h"
// #include "Operator.h"


#define nOptimVar 4
const char * OPTIM_VAR_NAMES[nOptimVar] = 
  {"never", "respect bounds", "try", "yes"}; // keep yes last !!

void ResetWarnings(int * all) {
  internal_param *w = &(GLOBAL.internal);
  w->warn_oldstyle = w->warn_newstyle = w->warn_Aniso = 
    w->warn_normal_mode = w->warn_mode = w->warn_coordinates =
    w->warn_on_grid = w->warn_new_definitions = w->warn_aspect_ratio = 
    w->warn_coord_change = w->warn_color_palette = w->warn_zenit =
    w->warn_constant = w->warn_negvar = w->warn_onlyvar =
    true;
  if (*all) w->warn_ambiguous = true;
}
 
void GetCurrentNrOfModels(int *init, int *nr) {
  if (*init && currentNrCov==-1) InitModelList();
  *nr = currentNrCov;
}

#define startmode normal
#define NM ((int) startmode)   /* normal mode */
//       {careless, sloppy, easygoing, normal, precise, pedantic, neurotic}
bool 
  allowdistance0[nr_modes] = {true, false, false, false, false, false, false},
  skipchecks[nr_modes] =  {true, false, false, false, false, false, false},
  ce_force[nr_modes] =       {true, true, true, false, false, false, false},
  ce_dependent[nr_modes] =   {true, true, true, false, false, false, false},
  sp_grid[nr_modes] =           {true, true, true, true, true, false, false},
  fit_reoptimise[nr_modes] = {false, false, false, false, true, true, true},
  fit_ratiotest_approx[nr_modes]={true, true, true, true, false, false, false},
  fit_cross_refit[nr_modes] = {false, false, false, false, true, true, true}
  ;
char pch[nr_modes] =         {'\0',  '\0',  '\0',  '.',   '.',   '.',   '.'}
  ;
int locmaxn[nr_modes] =      {3000, 4000, 6000, 8000, 9000, 10000000, 10000000},
  ce_trials[nr_modes] =      {1,    2,    3,    3,    4,    5,   6},
  spectral_lines[nr_modes] = {300, 500, 1500,  2500, 3500, 5000, 10000},
  tbm_lines[nr_modes] =      {30,   40,   50,   60,   80,   100,   200},
  mpp_n_estim_E[nr_modes] =  {10000, 10000, 20000,50000,80000,100000,1000000},
  hyper_superpos[nr_modes] = {200, 300, 450, 700, 1000, 2000, 5000},
  maxclique[nr_modes] =      {1000, 2000, 3500, 5000, 6000, 7000, 10000000},
  fit_critical[nr_modes] =   {-1,  -1,   -1,   0,   1,   2,   3}, 
  fit_ncrit[nr_modes] =      { 2,   4,    5,   5,   5,  10,  20}, 
  //  fit_optim_var[nr_modes] =  { 2,   2,    2,   2,   1,   1,   1} ,
  fit_split[nr_modes] =      {10000,  5,  4, 4, 4, 3, 3},
  sequ_back_initial[nr_modes]={3,   3,    5,  10,  20,  30,  40}
  ;
double 
  matrixtolerance[nr_modes] ={1e-6,  1e-6,  1e-6,  1e-12,   1e-14, 0, 0},
  exactness[nr_modes] = {false, false, false, RF_NA, true, true, true},
  ce_tolRe[nr_modes] =       {-1e2,  -1e1,  -1e-3, -1e-7,  -1e-14, 0, 0},
  ce_tolIm[nr_modes] =       {1e2,  1e1,  1e-1,  1e-3,   1e-7, 0, 0},
  ce_approx_step[nr_modes] = {1.0,  1.0,   RF_NA, RF_NA, RF_NA, 0.0, 0.0},
  direct_tol[nr_modes] =     {1e-5, 1e-7, 1e-8, 1e-9, 1e-11, 1e-13, 0},
  nugget_tol[nr_modes] =     {1e-8, 1e-8, 1e-12, 0,     0,     0,     0},
  tbm_linefactor[nr_modes] = {1.5,  1.5,  1.7,   2.0,  3.0,  5.0,  6.0},
  mpp_intensity[nr_modes] =  {50, 50, 80, 100, 200, 500, 1000},
  mpp_zero[nr_modes] =       {1e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7},
  max_max_gauss[nr_modes] =  {2,    2,    3,    4,    5,   6,   6},
  fit_pgtol[nr_modes]       ={10, 1e-2, 1e-3, 1e-4, 1e-4, 1e-6, 0},
  fit_pgtol_recall[nr_modes]={10, 1e-1, 1e-2, 1e-3, 1e-3, 1e-4, 0},
  fit_factr[nr_modes]       ={1e14, 1e13, 1e12, 1e11, 1e11, 1e9, 1e7},
  fit_factr_recall[nr_modes]={1e15, 1e14, 1e13, 1e12, 1e12, 1e10, 1e8}
   ;

const char *f_opt[nr_modes] = {"optim", "optim", "optim", "optim", "optim", "optim", "optim"}; // to do optimx


int PL=C_PRINTLEVEL, 
    NS=NAT_SCALE;
globalparam GLOBAL = {
  general_START,
  gauss_START,
  krige_START, 
  ce_START,
  spectral_START,
  tbm_START,
  direct_START,
  sequ_START,
  ave_START,
  nugget_START,
  mpp_START,
  hyper_START,
  extreme_START,
  br_START,
  distr_START,
  fit_START, 
  empvario_START, 
  gui_START,  
  graphics_START,
  register_START,
  internal_START,
  coords_START, 
  special_START,
  solve_param_default // solve
};
 

void SetDefaultOutputModeValues(output_modes mode){
  general_param *gp = &(GLOBAL.general);
  //printf("%d %d\n", mode, output_sp);
  gp->output = mode;
  gp->sp_conform = mode == output_sp;
  gp->returncall = mode == output_geor;
  gp->reportcoord = 
    mode == output_geor ? reportcoord_always : reportcoord_warnings;
  /*
  switch(mode) {
  case output_sp : 
     break;
  case output_rf :
  break;
  case output_geor :
     break;
  default: BUG;
  }
  */
}


void SetDefaultModeValues(int old, int m){
  int i;
  //                      high  fast, normal, save, pedantic
  GLOBAL.general.skipchecks = skipchecks[m];
  GLOBAL.general.pch = pch[m];
  GLOBAL.general.exactness = exactness[m];
  GLOBAL.general.matrixtolerance = matrixtolerance[m];
  GLOBAL.general.allowdist0 = allowdistance0[m];
  GLOBAL.krige.locmaxn = locmaxn[m];
  GLOBAL.krige.locsplitn[0]  = MINSPLITN(locmaxn[m]);
  GLOBAL.krige.locsplitn[1]  = MEDSPLITN(locmaxn[m]);
  GLOBAL.krige.locsplitn[2]  = MAXSPLITN(locmaxn[m]);
  GLOBAL.ce.force = ce_force[m];
  GLOBAL.ce.tol_im = ce_tolIm[m];
  GLOBAL.ce.tol_re = ce_tolRe[m];
  GLOBAL.ce.trials = ce_trials[m];
  GLOBAL.ce.dependent = ce_dependent[m];
  GLOBAL.ce.approx_grid_step = ce_approx_step[m];
  GLOBAL.direct.svdtolerance = direct_tol[m];
  GLOBAL.nugget.tol = nugget_tol[m];
  for(i=0; i<4; GLOBAL.spectral.lines[i++] = spectral_lines[m]);
  GLOBAL.spectral.grid = sp_grid[m];
  GLOBAL.tbm.lines[1] = tbm_lines[m];
  GLOBAL.tbm.lines[2] = (tbm_lines[m] * 25) / 6;
  GLOBAL.tbm.linesimufactor = tbm_linefactor[m];
  GLOBAL.tbm.grid = sp_grid[m];
  GLOBAL.mpp.n_estim_E = mpp_n_estim_E[m];
  for(i=0; i<4; GLOBAL.mpp.intensity[i++] = mpp_intensity[m]);
  GLOBAL.mpp.about_zero = mpp_zero[m];
  GLOBAL.hyper.superpos = hyper_superpos[m];
  GLOBAL.extreme.standardmax = max_max_gauss[m];
 GLOBAL.fit.locmaxn = maxclique[m];
  GLOBAL.fit.locsplitn[0] = MINCLIQUE(maxclique[m]);
  GLOBAL.fit.locsplitn[1] = MEDCLIQUE(maxclique[m]);
  GLOBAL.fit.locsplitn[2] = MAXCLIQUE(maxclique[m]);
  fit_param *f = &(GLOBAL.fit);
  f->split = fit_split[m];
  f->reoptimise = fit_reoptimise[m];
  f->ratiotest_approx = fit_ratiotest_approx[m];
  f->cross_refit = fit_cross_refit[m];
  f->critical = fit_critical[m];
  f->n_crit = fit_ncrit[m];
  f->pgtol = fit_pgtol[m];
  f->pgtol_recall = fit_pgtol_recall[m];
  f->factr = fit_factr[m];
  f->factr_recall = fit_factr_recall[m];
  //  f->optim_var_estim = fit_optim_var[m];
   char dummy[10];
  strcpy(dummy, f_opt[m]);
  f->optimiser = Match(dummy, OPTIMISER_NAMES, nOptimiser);
  GLOBAL.sequ.initial = -(GLOBAL.sequ.back = sequ_back_initial[m]);

  //  printf("optimiser %d %s\n", f->optimiser,  OPTIMISER_NAMES[f->optimiser]);

  internal_param *w = &(GLOBAL.internal);
  w->stored_init = false;
  if (m < normal) 
    w->warn_Aniso = w->warn_ambiguous =
      w->warn_on_grid = w->warn_aspect_ratio = w->warn_color_palette = false;
  else if (m>old) { 
    w->warn_oldstyle = w->warn_Aniso = w->warn_new_definitions = 
      w->warn_aspect_ratio = true;
    if (m > normal) w->warn_ambiguous = w->warn_negvar = true;
  }
  if (m != normal && w->warn_mode) {
    PRINTF("Note that the option '%s' is still in an experimental stage, so that the behaviour may change (slightly) in future.", general[GENERAL_MODUS]);
    w->warn_mode = false;
  }
}
 


SEXP GetAllModelNames(){
  if (currentNrCov==-1) InitModelList();
  int i, n;
  SEXP names;
  for (n=i=0; i<currentNrCov; i++) 
    if (CovList[i].name[0] != '-') n++;
  PROTECT(names = allocVector(STRSXP, n)); 
  for (n=i=0; i<currentNrCov; i++) {
    if (CovList[i].name[0] != '-') {
      SET_STRING_ELT(names, n++, mkChar(CovList[i].name));
    }
  }
  UNPROTECT(1);
  return names;
}
  
  

void GetModelName(int *nr,char **name, char **nick){
  if (currentNrCov==-1) InitModelList();
  if ((*nr<0) ||(*nr>=currentNrCov)) {
    strcopyN(*name,"", MAXCHAR); 
    strcopyN(*nick,"", MAXCHAR); 
    return;
  }
  strcopyN(*name, CovList[*nr].name, MAXCHAR);
  strcopyN(*nick, CovList[*nr].nick, MAXCHAR);
}

void GetNrParameters(int *covnr, int* kappas) {
  if (currentNrCov==-1) InitModelList();
  if (*covnr<0 || *covnr>=currentNrCov) {*kappas=-999;}
  else *kappas = CovList[*covnr].kappas;
}


void GetModelNr(char **name, int *nr) {
  *nr = getmodelnr(*name);
}


SEXP GetParameterNames(SEXP nr) {
  if (currentNrCov==-1) InitModelList();
  cov_fct *C = CovList + INTEGER(nr)[0]; // nicht gatternr
  SEXP pnames;
  int i;
  // print("hello %s\n", C->name);

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
    SET_STRING_ELT(pnames, i, mkChar(CAT_TYPENAMES[i]));
  }
  UNPROTECT(1);
  return(pnames);
}

SEXP GetSubNames(SEXP nr) {
  cov_fct *C = CovList + INTEGER(nr)[0]; // nicht gatternr
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
    if (C->subintern[i]) 
      PRINTF("%s subintern[%d]=true\n", C->nick, i);
    
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
  cov_fct *C;
  cov_model cov;
  range_type *r;
  int i,j, kappas;

  if (currentNrCov==-1) InitModelList();
  if ((*nr<0) || (*nr>=currentNrCov)) goto ErrorHandling;
  C = CovList + *nr; // nicht gatternr
  kappas = CovList[*nr].kappas;
  if (*lparam > kappas || *lrange != kappas * 4) goto ErrorHandling;
  if (*index < 0) {
    getrange->n = 1;
    for (i=0; i<lparam; i++) {
      cov.p[i] = (double *) CALLOC(sizeof(double), 1);
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
 *index = -100;
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
{
    int i, k, last_method, m, OP;
//char header[]="circ cut intr TBM spec dir seq Mak ave add hyp part\n";
  char coded[6][2]={"-", "X", "+", "N", "H", "S"};
//  char typenames[4][2]={"i", "f", "s", "n"};
  char specialnames[4][2]={".", "n", "f", "?"};
  char firstcolumn[20], name[MAXCHAR];
  int maxchar=10; // <= MAXCHAR=14
  int type[MAXNRCOVFCTS], op[MAXNRCOVFCTS], monotone[MAXNRCOVFCTS], 
    finite[MAXNRCOVFCTS], simpleArguments[MAXNRCOVFCTS], 
    internal[MAXNRCOVFCTS], dom[MAXNRCOVFCTS], 
    iso[MAXNRCOVFCTS], vdim[MAXNRCOVFCTS], maxdim[MAXNRCOVFCTS],
    paramtype[MAXNRCOVFCTS * MAXPARAM],
    nick = *Nick;
  
  last_method = (int) Nothing; // even not special method   
  if (currentNrCov==-1) { 
     InitModelList(); 
   // if (PL>5)  PRINTF("List of covariance functions initiated.\n"); 
  }
  if (CovList==NULL) {PRINTF("There are no functions available!\n");} 
  else {
    cov_fct *C;
    int includevariants = false, nr;

    GetAttr(NULL, type, op, monotone, finite, simpleArguments,
	    internal, dom, iso, maxdim, vdim, &includevariants,
	    paramtype,
	    &nr);
 
    sprintf(firstcolumn,"%%%ds", -maxchar);
    PRINTF("\n\n");
    PRINTF("%20s      List of models\n", "");
    PRINTF("%20s      ==============\n", "");
    PRINTF("%10s[See also PrintMethodList for the names of the columns();\n", 
	   "");
    PRINTF("%10s use 'operator=TRUE' to see all available models        ]\n", 
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
      PMLheader(firstcolumn, nick);

      for (k=1, i=0; i<currentNrCov; i++, C++) { 
	if (!isPosDef((Types)(type[i])) && !isUndefined((Types)(type[i])))
	  continue;
	if (op[i] != OP) continue;
	if (!*intern && internal[i]) continue;
	strcopyN(name, C->name, maxchar);
	if (strncmp(C->name, InternalName, strlen(InternalName)) ==0) {
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
	       specialnames[isNormalMixture(monotone[i]) ? 1
			    : finite[i] == 1 ? 2 
			    : isUndefined((Types)(type[i])) || 
			    monotone[i]<0 || finite[i] < 0 ? 3
			    : 0]	                          
	       );	  
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
  if (currentNrCov==-1) InitModelList(); 

  if (CovList==NULL) return;
  for (j=i=0; i<currentNrCov; i++) {
    if (!*internal && CovList[i].internal) continue;
    for (m=(int) CircEmbed; m<(int) Nothing; m++) {
      idx[j++] = CovList[i].implemented[m];
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
  cov_fct *C = CovList;
  assert((*includevariants) xor (Nr == NULL));

  for (j = nr = 0; nr<currentNrCov; nr++, C++){   
    int v, variants = *includevariants ? C->variants : 1;
    for (v=0; v<variants; v++, j++) {
       type[j] = C->Typi[v]; 
       dom[j] = C->domain;
      iso[j] = C->Isotropy[v];
      if (*includevariants) Nr[j] = nr;
      vdim[j] = C->vdim;
      op[j] = (int) C->maxsub > 0;   
      maxdim[j] = C->maxdim;
      finiterange[j] = C->finiterange;
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



void getUnits(SEXP el, char VARIABLE_IS_NOT_USED *name, 
	      char units[MAXUNITS][MAXUNITSCHAR], 
	      char units2[MAXUNITS][MAXUNITSCHAR]) {  
  int i, j, 
    l = length(el);
  if (TYPEOF(el) != NILSXP && TYPEOF(el) == STRSXP && l >= 1) {
    for (i=j=0; i<MAXUNITS; i++, j=(j + 1) % l) {
      strcopyN(units[i], CHAR(STRING_ELT(el, j)), MAXUNITSCHAR);
      if (units2!=NULL) {
	strcopyN(units2[i], CHAR(STRING_ELT(el, j)), MAXUNITSCHAR);
      }
    }
  } else ERR("invalid units");
}

SEXP UNITS(char units[MAXUNITS][MAXUNITSCHAR]) {
  SEXP unitnames;
  int nn;
  PROTECT(unitnames = allocVector(STRSXP, MAXUNITS)); 
  for (nn=0; nn<MAXUNITS; nn++) {
    SET_STRING_ELT(unitnames, nn, mkChar(units[nn]));
  }
  UNPROTECT(1);
  return unitnames;
}



void CE_set(SEXP el, int j, char *name, ce_param *cp, bool isList) {
  char msg[200];
  switch(j) {
  case 0: cp->force = LOG; break;
  case 1:
    Real(el, name, cp->mmin, MAXCEDIM) ;
    int d;
    for (d=0; d<MAXCEDIM; d++) {
      if (cp->mmin[d]<0.0) {
	if ((cp->mmin[d]>-1.0)) {
	  cp->mmin[d] = -1.0;
	  sprintf(msg, "'%s' set to -1.0.\n", name);
	  warning(msg);
	}
      }
    }
    break;
  case 2: int strat;
    strat = INT;
    if (strat>LASTSTRATEGY) {  
      sprintf(msg, "%s <= %d not satisfied\n", name, LASTSTRATEGY);
      warning(msg);
    } else cp->strategy= (char) strat;          
    break;
  case 3: {
    cp->maxGB = POSNUM;
    if (!isList) cp->maxmem = MAXINT;	
  } break;
  case 4: {
    cp->maxmem = POSINT;  
    if (!isList) cp->maxGB = RF_INF;	
  }  break;
  case 5: cp->tol_im = POS0NUM; break;
  case 6: cp->tol_re = NEG0NUM; break;
  case 7: cp->trials = NUM;
    if (cp->trials<1) {
      cp->trials=1;
      sprintf(msg, "%s is set to 1\n", name);
      warning(msg);
    }
    break;
  case 8: cp->useprimes = LOG; break;
  case 9: cp->dependent = LOG; break;
  case 10: cp->approx_grid_step = POS0NUM; break;
  case 11: cp->maxgridsize = POS0INT; break;
  default: ERR("unknown parameter for GLOBAL.general");
  }
}
 


#define prefixN 25
const char * prefixlist[prefixN] = 
  {"",  // -1
   "general", "gauss", "krige", 
   "circulant", "direct",  "nugget",
   "sequ", "spectral", "tbm",
   "mpp", "hyper", "maxstable",  // 12
   "br", "distr", "fit",         // 15
   "empvario", "gui", "graphics",// 18 
   "registers", "internal", "coords", "special",
   "solve",
   "obsolete"// 24
};


// IMPORTANT: all names of general must be at least 3 letters long !!!
const char *general[generalN] =
  { "modus_operandi", "printlevel", "storing", 
    "skipchecks", "every", "gridtolerance",
    "pch", "practicalrange", "spConform",
    "cPrintlevel", "exactness", "matrix_inversion", 
    "matrix_tolerance",  "allowdistanceZero", "na_rm_lines",
    "vdim_close_together", "expected_number_simu", "seed", 
    "detailed_output", "asList", "Ttriple",
    "returncall", "output", "reportcoord", "set"};


const char *gauss[gaussN]= {"paired", "stationary_only", "approx_zero", 
			    "direct_bestvar", "loggauss", "boxcox"};

const char *krige[krigeN] = { "return_variance", "locmaxn",
			     "locsplitn", "locsplitfactor", "fillall"};

const char *CE[CEN] = {"force", "mmin", "strategy", "maxGB",
		       "maxmem", "tolIm","tolRe", 
		       "trials", "useprimes", "dependent", 
		       "approx_step", "approx_maxgrid"};

const char *direct[directN] = {"root_method", "svdtolerance", "max_variab"};

const char * pnugget[pnuggetN] ={"tol"};

const char * sequ[sequN] ={"max_variables", "back_steps", "initial"};

const char * spectral[spectralN] = {"sp_lines", "sp_grid",  
				    "prop_factor", "sigma"};

const char * pTBM[pTBMN] = {"reduceddim", "fulldim", "center", 
			    "points", "lines", "linessimufactor", 
			    "linesimustep", "layers", "grid"};

const char * mpp[mppN] = {"n_estim_E", // n to determine E by simulation
			  "intensity", 
			  // "refradius_factor", 
			  "about_zero", "shape_power",
			  "scatter_max", "scatter_step",
			  // "plus", 
			  // "samplingdist", "samplingr",// MPP_cc
			  //"p", // Gneiting_cc
                         };

const char * hyper[hyperN] = {"superpos", "maxlines", "mar_distr",
			      "mar_param"};

const char * extreme[extremeN] = 
  {"max_gauss", "maxpoints", "xi",
   "density_ratio", "check_every", "flat",
   "min_n_zhou", "max_n_zhou", "eps_zhou",
   "mcmc_zhou"};

const char * br[brN] = 
  {"maxtrendmem", "meshsize", 
   "vertnumber", "optim_mixed", "optim_mixed_tol", 
   "optim_mixed_maxpoints", "variobound", "deltaAM", "corr_factor"};

const char * distr[distrN] = 
  {"safety", "minsteplen", "maxsteps", 
   "parts", "maxit",  "innermin", 
   "outermax", "mcmc_n", "repetitions"};

const char * fit[fitN] = 
  {"bin_dist_factor", "upperbound_scale_factor", "lowerbound_scale_factor", //2
   "lowerbound_scale_ls_factor","upperbound_var_factor","lowerbound_var_factor",
   //"lowerbound_sill",
   "scale_max_relative_factor", "minbounddistance", "minboundreldist",  //8
   "approximate_functioncalls", "minmixedvar",  "maxmixedvar",  //11
   "boxcox_lb",  "boxcox_ub", //13
   "use_naturalscaling",  "bins",  "nphi", //16
   "ntheta", "ntime", "only_users", //19
   "shortnamelength",  "split",   "scale_ratio", //22
   "critical",  "n_crit",  "max_neighbours",  //25
   "cliquesize",  "splitfactor_neighbours",  "smalldataset", 
   "min_diag",  "reoptimise",   "optimiser", 
   "algorithm",  "likelihood",  "ratiotest_approx", 
   "split_refined",  "cross_refit",  "estimate_variance", 
   "pgtol", "pgtol_recall",  "factr",
   "factr_recall"
  };


const char * empvario[empvarioN] = 
  {"phi0", "theta0", "tol0",
   "pseudovariogram", "fft"};

const char * gui[guiN] = 
  {"alwaysSimulate", "simu_method", "size"};

const char *graphics[graphicsN]= 
  {"always_close_screen" ,"grPrintlevel", "height", 
   "increase_upto", "always_open_screen", "file", 
   "onefile", "filenumber", "resolution"};

const char *registers[registersN] = 
  {"register", "predict_register", "likelihood_register"};

const char * internals[internalN] =  {
  // Achtung ! warn parameter werden nicht pauschal zurueckgesetzt
  "warn_oldstyle", "warn_newstyle", "warn_newAniso", 
  "warn_ambiguous", "warn_normal_mode",  "warn_mode",
  "stored.init",  "warn_scale",  "warn_on_grid", 
  "warn_coordinates", "warn_new_definitions", "warn_aspect_ratio",
  "warn_coord_change", "warn_colour_palette", "warn_missing_zenit",
  "do_tests", "warn_constant", "warn_negvar", "warn_onlyvar",
  "examples_reduced"
};

const char *coords[coordsN] =
  { "xyz_notation", "coord_system", "new_coordunits",
    "coordunits", "varunits", "varnames",
    "coordnames", "new_coord_system", "zenit", 
    "polar_coord"};

const char * special[specialN] = {"multicopies"};

const char * solve[solveN] = 
  { "use_spam", "spam_tol", "spam_min_p", "svd_tol",
    "spam_methods", "spam_min_n", "spam_sample_n", "spam_factor",
    "spam_pivot"
  };

const char * obsolete[obsoleteN] = 
  { "oldstyle", "newstyle",  "newAniso", 
    "ambiguous", "normal_mode", "solvesigma","optim_var_elimination", "sill"
};



const char **all[] = {general, gauss, krige, CE, direct, 
		      pnugget, sequ, spectral, pTBM, mpp,
		      hyper, extreme, br, distr, fit, 
		      empvario, gui, graphics, registers, internals, 
		      coords, special, solve, obsolete};

int allN[] = {generalN, gaussN, krigeN, CEN, directN,
	      pnuggetN,  sequN, spectralN, pTBMN, mppN,
	      hyperN, extremeN, brN, distrN, fitN, 
	      empvarioN, guiN, graphicsN, registersN, internalN, 
	      coordsN, specialN, solveN, obsoleteN};

void RelaxUnknownRFoption(int *relax) {
  RELAX_UNKNOWN_RFOPTION = (bool) *relax;
}


void setparameter(SEXP el, char *prefix, char *mainname, bool isList) {  
  int i,j,ii;
  char msg[LENMSG], name[200];

  //  printf("setting\n");

  isList &= GLOBAL.general.asList;
  
  sprintf(name, "%s%s%s", prefix, strlen(prefix)==0 ? "" : ".", mainname);

  //  print("set param: %s.%s.%s\n",prefix, strlen(prefix)==0 ? "" : ".", mainname);
  //  print("relax=%d\n", RELAX_UNKNOWN_RFOPTION);

  if (mainname[0] >= 'A' && mainname[0] <= 'Z' && RELAX_UNKNOWN_RFOPTION) {
    if (PL > PL_IMPORTANT)
      PRINTF("'%s' is not considered as an RFoption, but will be passed to evaluate the model formula.\n", mainname);
    return;
  }


  i = Match(prefix, prefixlist, prefixN);
  if (i<0) {
    sprintf(msg, "option unknown (unknown prefix: %s)", prefix);
    ERR(msg);
  } 
 
  //  print("%d %d pref='%s' %s relax=%d\n", i, j,  prefix, mainname, RELAX_UNKNOWN_RFOPTION);


  if (i==0) { 
#define MinNameLength 3

    int trueprefixN = prefixN - 1;
    for (i=0; i < trueprefixN; i++) {	
      //  printf("%s\n", prefixlist[ii+1]);
      j = Match(mainname, all[i], allN[i]);
      if (j != NOMATCHING) break;     
    }

    //   if (j>=0) printf("%s %s %d %d=%d %d\n", prefixlist[i+1], all[i][j], j,
    //		     strlen(mainname), strlen(all[i][j]),
    //		     strcmp(mainname, all[i][j]) != 0); else 
    //     printf("not found");
    
    
    bool ok = j >= 0;
    if (j != NOMATCHING && (!ok || strcmp(mainname, all[i][j]) != 0))  {
      for (ii=i + 1; ii < trueprefixN; ii++) {	
	int jj = Match(mainname, all[ii], allN[ii]);
	if (jj >= 0) {
	  if ((ok = strcmp(mainname, all[ii][jj]) == 0)) {
	    i = ii;
	    j = jj;
	    break;
	  }
	}
      }
      if (!ok) {
 	sprintf(msg, "%s '%s' cannot be uniquely identified.", 
		RFOPTIONS, name);
	ERR(msg);      
      }
    }
  } else {   
    i--; // since general has to prefixed.
    j = Match(mainname, all[i], allN[i]);
  }

  if (j<0) {
    if (j==-1) sprintf(msg, "Unknown option '%s'.", name);
    else sprintf(msg, "Multiple matching for '%s'.", name);
    ERR(msg);
  }

  //printf("i=%d %d\n", i, j);

  switch(i) {
  case 0: {// general
    general_param *gp;
    gp = &(GLOBAL.general);
    switch(j) {
    case GENERAL_MODUS: {
      int old_mode = gp->mode;
      static bool warned = false;
      gp->mode = GetName(el, name, MODENAMES, nr_modes, normal);
      if (old_mode != gp->mode && !warned) {
	warned = true;
	warning("Note that behaviour of 'modus_operandi' has changed within 'RFfit' in version 3.1.0 of RandomFields. Roughly:\nwhat was called 'careless' is now called 'sloppy';\nwhat was called 'sloppy' is now called 'easygoing';\nwhat was called 'easygoing' is now called 'normal';\nwhat was called 'normal' is now called 'precise';\netc.");
      }
      SetDefaultModeValues(old_mode, gp->mode);
      break;
    }
    case 1: {
      int threshold = 1000; //PL_ERRORS;
      gp->Rprintlevel = INT;
      PL = gp->Cprintlevel = 
	gp->Rprintlevel <= threshold ? gp->Rprintlevel : threshold;


    }
      break;
    case GENERAL_STORING: {
      bool storing = LOG;
      //  print("before setting storing %d %d\n", storing, KEY[0].simu.active);
      if  (length(el) > 1) {
	if (!storing) {	  
	  int nr = Integer(el, (char*) "storing (register)", 1);
	  if (nr != NA_INTEGER) {
	    //	    print("deleting register %d\n", nr);
	    if (nr<0 || nr>MODEL_MAX) 
	      ERR1("%s: number for register is out of range", RFOPTIONS);
	    COV_DELETE(KEY + nr);
	    //	  print("xstoring = %d\n", storing);
	  }
	  
	  if (length(el) >=3) {
	    nr = Integer(el, (char*) "storing (model)", 2);  
	     if (nr != NA_INTEGER) {
	       if (nr>=0 && nr<=MODEL_MAX && KEY[nr] != NULL) 
		 COV_DELETE(KEY + nr);
	     }
	  }
	}
	//
	// print("hereX %d %d %d\n", storing, length(el), KEY[0].simu.active);
      } else {
	if (!storing) {	  
	  // delete all keys
	  for (ii=0; ii<=MODEL_MAX; ii++) {
	    if (KEY[ii]!=NULL) COV_DELETE(KEY + ii);
	  }
	  //
	}
	//
	//	print("hereR %d %d\n", storing, KEY[0].simu.active);
     }
      //
      // print("here %d\n", storing);
      gp->storing = storing;}
      break;
    case 3: gp->skipchecks = LOG;    break;
    case 4: gp->every = POS0INT;      break;
    case 5: gp->gridtolerance = NUM; break;
    case 6: gp->pch = CHR;           break;
    case 7: 
      int n;
      n =INT;  
      if (n!=0 && n!=1 && n!=2 && n!=3)
	ERR("PracticalRange out of range. It should be TRUE or FALSE.");
      NS = gp->naturalscaling = n;
      break;
    case 8: 
      SetDefaultOutputModeValues(LOG ? output_sp : output_rf);
      break;
    case GENERAL_CPRINT: PL = gp->Cprintlevel = INT;        break;
    case GENERAL_EXACTNESS: gp->exactness = NUM;        break; 
    case 11: {
      int old[MAXINVERSIONS],
	l = length(el);
      if (l > MAXINVERSIONS) ERR("matrix_inversion: vector too long");
      for (ii=0; ii<l; ii++) {
	old[ii] = gp->matrix_inversion[ii];
	gp->matrix_inversion[ii] = Integer(el, name, ii); 
	if (gp->matrix_inversion[ii] > 3) break;
      }
      if (ii<l) {
	for (; ii>=0; ii--) gp->matrix_inversion[ii] = old[ii];
	ERR("values of matrix_inversion out of range");
      } else {
	for (; ii<MAXINVERSIONS; ii++) gp->matrix_inversion[ii] = -1;
      }
      break;
    }
    case 12: gp->matrixtolerance = NUM; break;
    case 13: gp->allowdist0 = LOG; break;
    case 14: gp->na_rm_lines = LOG; break;
    case GENERAL_CLOSE: gp->vdim_close_together = LOG;    
      if (gp->vdim_close_together) {
	gp->vdim_close_together = false;
 	ERR1("'%s' will likely not be programmed", general[GENERAL_CLOSE]);
	// Achtung! gausslikeli.cc waere stark betroffen!!
      }
      break;
    case 16: gp->expected_number_simu = POS0INT; break;
    case 17: gp->seed = Integer(el, name, 0, true); break;
    case 18: gp->detailed_output = LOG; break;
    case 19: gp->asList = LOG; break;
    case 20: gp->Ttriple = INT; break;
    case 21: gp->returncall = LOG; break;
    case 22:
      SetDefaultOutputModeValues((output_modes) 
				 GetName(el, name, OUTPUTMODENAMES, 
					 nr_output_modes, gp->output));
      break;
    case 23:
      gp->reportcoord = GetName(el, name, REPORTCOORDNAMES, 
				nr_reportcoord_modes, gp->reportcoord);
      break;
    case 24: gp->set = POSINT - 1; break;
    default: BUG;
    }}
    break;
  case 1: { // gauss
    gauss_param *gp;
    gp = &(GLOBAL.gauss);
    switch(j) {
    case 0: gp->paired = LOG; break;
    case 1: gp->stationary_only = NUM; break;
    case 2: gp->approx_zero = POS0NUM; break;
    case 3: gp->direct_bestvariables = POS0INT;    break;
    case 4: gp->loggauss = LOG;
      if (gp->loggauss) gp->boxcox[0] = RF_INF;
      break;
    case GAUSS_BOXCOX_OPTION: {
      int len = length(el);
      double boxcox[2 * MAXGAUSSVDIM];      
      if (len < 1 || (len > 2 && len % 2 != 0)) 
	ERR1("'%s' not a scalar or vector/matrix of multiple of length 2", 
	     gauss[GAUSS_BOXCOX_OPTION]);
      boxcox[1] = 0.0;
      for (int L=0; L<len; L++) boxcox[L] = Real(el, name, L);
      if (len == 1) len = 2;
      for (int k=0; k<len; k++) {
	if (R_FINITE(boxcox[k])) {
	  if (boxcox[k] < GLOBAL.fit.BC_lambdaLB[k]) 
	    ERR6("%s=%s[%d] < %s[%d]=%f.", 
		 boxcox[k], gauss[GAUSS_BOXCOX_OPTION], k, fit[FIT_BC_LB],k,
		 GLOBAL.fit.BC_lambdaLB[k]);
	  if (boxcox[k] > GLOBAL.fit.BC_lambdaUB[k]) 
	    ERR6("%f=%s[%d] > %s[%d]=%f.", 
		 boxcox[k],gauss[GAUSS_BOXCOX_OPTION], k, fit[FIT_BC_UB], k,
		 GLOBAL.fit.BC_lambdaUB[k]); 
	}
      }
      for (int k=0; k<len; k++) gp->boxcox[k] = boxcox[k];
      if (ISNA(gp->boxcox[0]) || gp->boxcox[0] != RF_INF) gp->loggauss = false;	
      break;
    }
    default: BUG;
    }}
    break;
  case 2: { // krige
    krige_param *kp;
    kp = &(GLOBAL.krige);
    switch(j) { 
    case 0: kp->ret_variance = LOG; break;
    case 1: kp->locmaxn = POS0INT; break;
    case KRIGE_SPLITN: {
      int  newval[3],
	len = length(el);
      if (len < 3) ERR1("'%s' must have 3 components", krige[KRIGE_SPLITN]);
      for (ii=0; ii<len; ii++) {
	newval[ii] = Integer(el, name, ii); 
	// printf("%d %d %d\n", ii, kp->locsplitn[ii],	newval[ii]);
      }
      if (newval[0] > newval[1] || newval[1] > newval[2]) {
	ERR6("%s[1] <= %s[2] <= %s[3] not satisfied [ %d %d %d ]",
	     krige[KRIGE_SPLITN], krige[KRIGE_SPLITN], krige[KRIGE_SPLITN],
	     newval[0], newval[1], newval[2]
	     ); 
      } 
      for (ii=0; ii<len; ii++) kp->locsplitn[ii] = newval[ii];
      break;
    }
    case 3: kp->locsplitfactor = POS0INT; break;
    case 4: kp->fillall = LOG; break;
    default: BUG;
    }}
    break;
  case 3: // CE
    CE_set(el, j, name, &(GLOBAL.ce), isList); 
    break;
  case 4: { //direct
    direct_param *dp;
    dp = &(GLOBAL.direct);
    switch(j) {
    case 0: int method;
      method = INT; 
      if (method<0 || method>=(int) NoFurtherInversionMethod){
	warning("given inversion method out of range; ignored\n");
      } else  dp->inversionmethod= (InversionMethod) method;
      break;
    case 1: dp->svdtolerance = NUM; break;
    case 2: dp->maxvariables = POS0INT;    break;
    default: BUG;
    }}
    break;
  case 5:  {// pnugget, 
    nugget_param *np;
    np = &(GLOBAL.nugget) ;
    switch(j) {
    case 0: np->tol = NUM; 
      if (np->tol < 0) {
	if (PL>=PL_IMPORTANT) {
	  strcpy(msg, "negative tolerance for distance in nugget covariance not allowed; set to zero");
	  warning(msg);
	}
	np->tol = 0.0;
      }
      break;
    default: BUG;
    }}
    break;
  case 6: {//sequ,
    sequ_param *sp;
    sp = &(GLOBAL.sequ) ;
    switch(j) {
    case 0: sp->max = POS0INT;   break;
    case 1: sp->back = INT; if (sp->back < 1) sp->back=1;   break;
    case 2: sp->initial = INT; break;
    default: BUG;
    }}
    break;
  case 7: { // spectral,
    spectral_param *sp;
    sp = &(GLOBAL.spectral) ;
    switch(j) {
    case 0: Integer(el, name, sp->lines, MAXTBMSPDIM); break;
    case 1: sp->grid = LOG; break;
    case SPECTRAL_PROPFACTOR: sp->prop_factor = POS0NUM;
      if (sp->prop_factor <= 0.1) {
	sp->prop_factor=0.1;
	WARNING1("'%s' less than 0.1. Set to 0.1.", 
		 spectral[SPECTRAL_PROPFACTOR]);
      }
      break;
    case 3: sp->sigma = NUM; break;
    default:  BUG;
    }}
    break;
  case 8: {// TBM
    tbm_param *tp;
    tp = &(GLOBAL.tbm) ;
    switch(j) {
    case 0: tp->tbmdim = INT; break;
    case 1: tp->fulldim = POS0INT; break;
    case 2: Real(el, name, tp->center, MAXTBMSPDIM);  break;
    case 3: tp->points = POS0INT; break;
    case 4: Integer(el, name, tp->lines, MAXTBMDIM); break;
    case 5: 
      tp->linesimufactor = Real(el, prefix, 0);
      
      if (tp->linesimufactor < 0.0) {
	sprintf(msg,"Both %s.linesimufactor and %s.linesimustep must be non-negative\n", prefix, prefix); 
	warning(msg);
	tp->linesimufactor = 0.0;
      }
      if (tp->linesimufactor>0.0 && tp->linesimustep>0.0) tp->linesimustep=0.0;
      break;
    case 6: 
      tp->linesimustep = Real(el, prefix, 0);  
      if (tp->linesimustep < 0.0) {
	sprintf(msg,"Both %s.linesimufactor and %s.linesimustep must be non-negative\n", prefix, prefix); 
	warning(msg);
	tp->linesimustep = 0.0;
      }
      if (tp->linesimufactor>0.0 && tp->linesimustep>0.0) 
	tp->linesimufactor=0.0;
      break;
    case 7:
      tp->layers = Real(el, prefix, 0);
      break;      
    case 8: 
      tp->grid = LOG; break;
    default:  BUG;
    }}
    break;
  case 9: {//  mpp, 
    mpp_param *mp;
    mp = &(GLOBAL.mpp);
    switch(j) {
    case 0: mp->n_estim_E = POS0INT; break;
    case 1: Real(el, name, mp->intensity, MAXMPPDIM); break;
      // case 2: mp->refradius_factor = POS0NUM; break;
    case 2: mp->about_zero = POS0NUM; break;
    case 3: mp->shape_power = NUM; break;
    case 4: {
      Integer(el, name, mp->scatter_max, MAXMPPDIM) ;
      break;
    }
    case 5: {
      Real(el, name, mp->scatter_step, MAXMPPDIM) ;
      break;
    }
    default: BUG;
   }}
    break;
  case 10: {//hyper, 
    hyper_param *hp;
    hp = &(GLOBAL.hyper);
    switch(j) {
    case 0: hp->superpos = POS0INT;  break;
    case 1: hp->maxlines = POS0INT;  break;
    case 2: hp->mar_distr = INT; break;
    case 3: hp->mar_param = NUM; break;
    default:  BUG;
    }}
    break;
  case 11: {//		 extreme
    extremes_param *ep;
    ep = &(GLOBAL.extreme);
    switch(j) {
    case 0: ep->standardmax = POS0NUM; break;
    case 1: ep->maxpoints = POS0INT; break;
    case 2: ep->GEV_xi = NUM; break;
    case 3: ep->density_ratio = POS0NUM; break;
    case 4: ep->check_every = POS0INT; break;
    case EXTREME_FLAT: ep->flat = INT; 
      if (ep->flat < -1 || ep->flat > 1) 
	ERR1("illegal value for '%s'", extreme[EXTREME_FLAT]);
      if (ep->flat != FALSE) 
	ERR1("currently only '%s=FALSE' allowed", extreme[EXTREME_FLAT]);
      // Programmierfehler in Huetchen.c
      break;
    case 6: ep->min_n_zhou = POSINT; break;
    case 7: ep->max_n_zhou = POSINT; break;
    case 8: ep->eps_zhou = POSNUM; break;
    case 9: ep->mcmc_zhou = POSINT; break;
    default:  BUG;
    }}
    break;
  case 12 : { // br
    br_param *ep;
    ep = &(GLOBAL.br);
    switch(j) {
    case 0: ep->BRmaxmem = POSINT; break;
    case 1: ep->BRmeshsize = POSNUM; break;
    case 2: ep->BRvertnumber = POSINT; break;
    case 3: ep->BRoptim = POS0INT; break;
    case 4: ep->BRoptimtol = POS0NUM; break;
    case 5: ep->BRoptimmaxpoints = POS0INT; break;
    case 6: ep->variobound = NUM; break;
    case 7: ep->deltaAM = POSINT; break;
    case 8: ep->corr_factor = POSNUM; break; // in [0,1]
    default:  BUG;
    }}	  
    break;
  case 13 : {// distr
    distr_param *ep;
    ep = &(GLOBAL.distr);
    switch(j) { 
    case 0: ep->safety=POSNUM; break;
    case 1: ep->minsteplen=POS0NUM; break;
    case 2: ep->maxsteps=POSINT; break;
    case 3: ep->parts=POSINT; break;
    case 4: ep->maxit=POSINT; break;
    case 5: ep->innermin=POSNUM; break;
    case 6: ep->outermax=POSNUM; break;
    case 7: ep->mcmc_n=POSNUM; break;
    case 8: ep->repetitions=POSNUM; break;
    default:  BUG;
    }}
    break;
  case 14: { // fit
    fit_param *ef;
    ef = &(GLOBAL.fit);
    switch(j) {
    case 0: ef->bin_dist_factor = POS0NUM; break; //
    case 1: ef->upperbound_scale_factor = POS0NUM; break; //
    case 2: ef->lowerbound_scale_factor = POS0NUM; break; //
    case 3: ef->lowerbound_scale_LS_factor = POS0NUM; break; //
    case 4: ef->upperbound_var_factor = POS0NUM; break; //
    case 5: ef->lowerbound_var_factor = POS0NUM; break;//
      //    case 6: ef->lowerbound_sill = POS0NUM; break; 
    case 6: ef->scale_max_relative_factor = POS0NUM; break; //
    case 7: ef->minbounddistance = POS0NUM; break;//
    case 8: ef->minboundreldist = POS0NUM; break;// 
    case 9: ef->approximate_functioncalls = POS0INT; break; //
    case 10: ef->minmixedvar = POS0NUM; break; //
    case 11: ef->maxmixedvar = POS0NUM; break; //
      //    case 14: ef->solvesigma = NUM; break;
    case 12: Real(el, name, ef->BC_lambdaLB, 2 * MAXGAUSSVDIM); break;//
    case 13: Real(el, name, ef->BC_lambdaUB, 2 * MAXGAUSSVDIM); break;//
    case 14: ef->use_naturalscaling = LOG; break;
    case 15: ef->bins = POS0INT; break;//
    case 16: ef->nphi = POS0INT; break;//
    case 17: ef->ntheta = POS0INT; break; //
    case 18: ef->ntime = POS0INT; break;//
    case 19: ef->onlyuser = LOG; break;
    case 20: { // mle
      ii=POS0INT; 
      if (ii==0) { ii = 1; warning("shortnamelength set to 1"); }
      if (ii >= LENMSG) {
	ii = LENMSG - 1;  
	sprintf(msg, "shortnamelength set to %d", ii); 
	warning(msg);
      }
      ef->lengthshortname=ii;
    }
    break;
    case 21: ef->split = INT; break;
    case 22: ef->scale_ratio = NUM; break;//
    case 23: ef->critical = INT; break;//
    case 24: ef->n_crit = POS0INT; break;///
    case FIT_MAXNEIGHBOUR: ef->locmaxn = POS0INT; break;//
    case FIT_CLIQUE: {
      int  newval[3],
	len = length(el);      
      if (len < 1 || len > 3) 
	ERR1("'%s' must have between 1 and 3 components", fit[FIT_CLIQUE]);
      for (ii=0; ii<len; ii++) newval[ii] = Integer(el, name, ii); 
      if (len == 1) newval[1] = newval[2] = newval[0];
      else {
	newval[2] = newval[1];
	newval[1] = (int) sqrt((double) newval[1] * newval[0]);
      }

      if (newval[0] > newval[1] || newval[1] > newval[2]) 
   	ERR3("%s[1] <= %s[2] <= %s[3] not satisfied",
	     fit[FIT_CLIQUE], fit[FIT_CLIQUE], fit[FIT_CLIQUE]
	     );      
      for (ii=0; ii<3; ii++) ef->locsplitn[ii] = newval[ii];
      break;
    }
    case 27: ef->locsplitfactor = POS0INT; break;//
    case 28: ef->smalldataset = POS0INT; break;//
    case 29: ef->min_diag = NUM; break;//
    case 30: ef->reoptimise = LOG; break;
    case 31: ef->optimiser = GetName(el, name, OPTIMISER_NAMES, nOptimiser, 0); //
      break;
    case 32: 
      switch(ef->optimiser) {
      case 3 : ef->algorithm = GetName(el, name, NLOPTR_NAMES, nNLOPTR); //
	break;
      default: ef->algorithm = -1;
      }
      break;
    case 33: ef->likelihood =
	GetName(el, name, LIKELIHOOD_NAMES, nLikelihood); ; break; 
    case 34: ef->ratiotest_approx = LOG; break;
    case 35: ef->split_refined = LOG; break;
    case 36: ef->cross_refit = LOG; break;
    case 37: ef->estimate_variance = INT; break;
    case 38: ef->pgtol = POSNUM; break;//
    case 39: ef->pgtol_recall = POSNUM; break;//
    case 40: ef->factr = POSNUM; break;//
    case 41: ef->factr_recall = POSNUM; break;//
    default:  BUG;
    }}
    break;
  case 15: { // empvario
    empvario_param *ep;
    ep = &(GLOBAL.empvario);
    switch(j) {
    case 0: ep->phi0=NUM; break;
    case 1: ep->theta0=NUM; break;    
    case 2: ep->tol=NUM; break;    
    case 3: ep->pseudovariogram = LOG; break;
    case 4: ep->fft = LOG; break;
   default: BUG;
    }}
    break;
  case 16: { // gui
    gui_param *gp;
    gp = &(GLOBAL.gui);
    switch(j) {
    case 0: gp->alwaysSimulate = LOG; break;
    case 1: 
      gp->method = GetName(el, name, METHODNAMES, Forbidden + 1, Nothing);
      break;
    case GUI_SIZE: {
      int sizedummy[2];
      if (length(el) != 2) ERR1("length of '%s' must be 2", gui[GUI_SIZE]);
      Integer(el, name, sizedummy, 2);
      for (ii=0; ii<2; ii++) {	
	if (sizedummy[ii] <= 1) 
	  ERR("grid in RFgui must contain at least 2 points");
      }
      for (ii=0; ii<2; ii++) {	 gp->size[ii] = sizedummy[ii]; }
    }
      break;
    default: BUG;
    }}
    break;
 case 17: { // graphics
    graphics_param *gp = &(GLOBAL.graphics);
    switch(j) {
    case 0 : gp->always_close = LOG; break;
    case 1 : gp->PL = INT; break;
    case 2 : gp->height = NUM; break;
    case GRAPHICS_UPTO : {
      int uptodummy[2];
      if (length(el) != 2)
	ERR1("length of '%s' must be 2", graphics[GRAPHICS_UPTO]);
      Integer(el, name, uptodummy, 2);
      for (ii=0; ii<2; ii++) {
	if (uptodummy[ii] <= 0)  ERR("increase_upto must be positive");
      } 
      for (ii=0; ii<2; ii++) { gp->increase_upto[ii] = uptodummy[ii]; }
    }
    break;
    case 4 : gp->always_open = INT;       
      break;
    case 5 :  if (!isList) {
	char old[100];
	strcopyN(old, gp->filename, 100);
	STR(gp->filename, 100);
	if (!gp->onefile && strcmp(gp->filename, old) !=0) gp->number = 0;
      }
      break;
    case 6 :  if (!isList) {
      bool onefile = gp->onefile;
      gp->onefile = LOG;
      if (!gp->onefile && onefile) gp->number = 0;
    } break;
    case 7 :  if (!isList) {
	gp->number = INT; 
      }
      break;
    case 8 : gp->resolution = POSNUM; break;
    default: BUG;
    }}
   break;
       /*
 case 20: { // empvario
    empvario_param *ep;
    ep = &(GLOBAL.empvario);
    switch(j) {
    case 0: 
      
     default: ERR("unknown option (empvario)");
   }}
    break;
    */ 
  case 18 : {
    registers_param *rp = &(GLOBAL.registers);
    switch(j) {
    case 0: { // simu
      int keynr = INT;
      if ((keynr<0) || (keynr>MODEL_MAX)) ERR("register number out of range");
      rp->keynr=keynr; }
      break;
   case 1: { //  predict
      int predict = INT;
      if ((predict<0) || (predict>MODEL_MAX))
	ERR("register number out of range");
      rp->predict=predict; }
      break;
    case 2: { //  likelihood
      int likelihood = INT;
      if ((likelihood<0) || (likelihood>MODEL_MAX))
	ERR("register number out of range");
      rp->likelihood=likelihood; }
      break;
   default: BUG;
    }}
   break;

  case 19: {
    internal_param *wp = &(GLOBAL.internal);
    // ACHTUNG internal wird nicht pauschal zurueckgesetzt !
    if (!isList) {
      switch(j) {
      case 0: wp->warn_oldstyle = LOG;       break;
      case 1: wp->warn_newstyle = LOG;       break;
      case INTERNALS_NEWANISO: wp->warn_Aniso = LOG;       break;
      case 3: wp->warn_ambiguous = LOG;       break;
      case 4: wp->warn_normal_mode = LOG;       break;
      case 5: wp->warn_mode = LOG;       break;
      case 6: wp->stored_init = LOG;       break;
      case 7: wp->warn_scale = LOG;       break;
      case 8: wp->warn_coordinates = LOG;       break;
      case INTERNALS_ONGRID: wp->warn_on_grid = LOG;       break;
      case 10: wp->warn_new_definitions = LOG;       break;
      case 11: wp->warn_aspect_ratio = LOG;       break;
      case INTERNALS_COORD_CHANGE: wp->warn_coord_change = LOG;       break;
      case 13: wp->warn_color_palette = LOG;       break;
      case INTERNALS_ZENIT: wp->warn_zenit = LOG;       break;
      case 15: wp->do_tests = LOG;       break;
      case 16: wp->warn_constant = LOG;       break;
      case 17: wp->warn_negvar = LOG;       break;
      case 18: wp->warn_onlyvar = LOG;       break;
      case 19: wp->examples_reduced= POS0INT;  
	if (wp->examples_reduced > 0 && 
	    wp->examples_reduced < MAX_LEN_EXAMPLES)
	  wp->examples_reduced = MAX_LEN_EXAMPLES;
	break;
    default: BUG;
      }
    } else {
      switch(j) {
      case 6: wp->stored_init = LOG;       break;
      case INTERNALS_ONGRID :  wp->warn_on_grid = LOG; break;
      case INTERNALS_COORD_CHANGE: wp->warn_coord_change = LOG;       break;
      case INTERNALS_ZENIT: wp->warn_zenit = LOG;       break;
      case 15: wp->do_tests = LOG;       break;
      case 19: wp->examples_reduced= POS0INT;  
	if (wp->examples_reduced > 0 && 
	    wp->examples_reduced < MAX_LEN_EXAMPLES)
	  wp->examples_reduced = MAX_LEN_EXAMPLES;
	break;
    default : {} // none
      }      
    }
  }
    break;

  case 20: {
    coords_param *cp = &(GLOBAL.coords);
    switch(j) {
    case COORDS_XYZNOTATION: cp->xyz_notation = INT; break;
    case 1: {
      coord_sys_enum coord =
	(coord_sys_enum) GetName(el, name, COORD_SYS_NAMES, nr_coord_sys,
				 coord_auto);
      if (coord == coord_auto || coord == earth || coord == sphere || 
	  coord == cartesian)
	cp->coord_system = coord;
      else ERR2("Cannot take '%s' as a value for the initial '%s'", 
		CHAR(STRING_ELT(el, 0)), coords[j]);
    } 
     break; 
    case 2: getUnits(el, name, cp->newunits, NULL);
      break;
    case 3: getUnits(el, name, cp->curunits, isList ? NULL : cp->newunits);
      break;
    case 4: getUnits(el, name, cp->varunits, NULL);      
      break;
    case 5: {
      SEXPTYPE type = TYPEOF(el);
      //printf("type=%d %d %d %d\n", type, INTSXP, REALSXP, LGLSXP);
      if (type == INTSXP || type == REALSXP || type == LGLSXP) {
	Integer2(el, name, cp->data_idx);
	cp->data_nr_names = 0;
      } else {
	String(el, name, cp->data_names);
	cp->data_nr_names = length(el);
      }
    }
      break;
    case 6: {
      SEXPTYPE type = TYPEOF(el);
      if (type == INTSXP || type == REALSXP || type == LGLSXP) {
	Integer2(el, name, cp->x_idx);
	cp->x_nr_names = 0;
      } else {
	String(el, name, cp->x_names);
	cp->x_nr_names = length(el);
      }
    }
      break;
    case 7: 
      cp->new_coord_system =
	(coord_sys_enum) GetName(el, name, COORD_SYS_NAMES, nr_coord_sys,
				 coord_keep);
      if (cp->new_coord_system == coord_auto) cp->new_coord_system = coord_keep;
      break; 

    case ZENIT: {
      Real(el, name, cp->zenit, 2);
    }
    break; 
 
    case 9: cp->polar_coord = LOG; 
      break;
   
   default: BUG;
    }}
    break;

  case 21: {
    special_param *sp = &(GLOBAL.special);
    switch(j) {
    case 0: sp->multcopies = POSINT;
	break;      
    default: BUG;
    }}
    break;

  case 22: {
    solve_param *so = &(GLOBAL.solve);
    switch(j) {
    case 0: so->sparse = NUM; break;
    case 1: so->spam_tol = POS0NUM; break;      
    case 2: so->spam_min_p = POS0NUM; break;      
    case 3: so->svd_tol = POS0NUM; break;      
    case 4: {
      int 
	len = length(el),
	ende = len;
      if (len > 0) {
	if (ende > SOLVE_METHODS) ende = SOLVE_METHODS;
	for (ii=0; ii<ende; ii++) so->Methods[ii] = Integer(el, name, ii);
	for ( ; ii<SOLVE_METHODS; ii++) so->Methods[ii] = so->Methods[ii-1];
      }
    }
      break;
    case 5: so->spam_min_n = POSINT; break;      
    case 6: so->spam_sample_n = POSINT; break;      
    case 7: so->spam_factor = POSINT; break;      
    case 8: so->pivot = POSINT; 
      if (so->pivot > 2) so->pivot = PIVOT_NONE;
      break;    
    default: BUG;
    }}
    break;

  case 23: { // obsolete
    internal_param *wp = &(GLOBAL.internal);
    switch(j) {
      case 0: wp->warn_oldstyle = LOG;       break;
      case 1: wp->warn_newstyle = LOG;       break;
      case 2: wp->warn_Aniso = LOG;       break;
      case 3: wp->warn_ambiguous = LOG;       break;
      case 4: wp->warn_normal_mode = LOG;       break;
    default:
      sprintf(MSG, "RFoption '%s' has been removed.", obsolete[j]);
      warning(MSG);
    }}
    break;
    

  default: BUG;
  }
}


SEXP ExtendedInteger(double x) {
  return ScalarInteger(R_FINITE(x) ? x : NA_INTEGER);
}

SEXP ExtendedBoolean(double x) {
  return ScalarLogical(ISNAN(x) ? NA_LOGICAL : x != 0.0);
}


SEXP getRFoptions() {
  SEXP list, names, sublist[prefixN-1], subnames[prefixN-1];
  
  int i, k = 0;
  char x[2]=" ";
  int trueprefixN = prefixN - 1 - 1;// general has two options for prefixes
  //                              and #22 obsolete is not shown
  //  cov_fct *C = CovList + cov->nr;

  PROTECT(list = allocVector(VECSXP, trueprefixN));
  PROTECT(names = allocVector(STRSXP, trueprefixN));
  
  for (i=0; i<trueprefixN; i++) {   
    if (strcmp(prefixlist[i+1], "obsolete") == 0) continue;
    SET_STRING_ELT(names, i, mkChar(prefixlist[i+1]));
    PROTECT(sublist[i] = allocVector(VECSXP, allN[i]));
    PROTECT(subnames[i] = allocVector(STRSXP, allN[i]));
    int endfor = allN[i];
    for (k=0; k<endfor; k++) {
      // print("%d %d %s$%s\n", i, k, prefixlist[i+1], all[i][k]);
      SET_STRING_ELT(subnames[i], k, mkChar(all[i][k])); 
    }
  }

  //#define ADD(ELT) {printf(#ELT"\n");SET_VECTOR_ELT(sublist[i], k++, ELT);}
#define ADD(ELT) SET_VECTOR_ELT(sublist[i], k++, ELT);
#define ADDCHAR(ELT) x[0] = ELT; ADD(ScalarString(mkChar(x)));
  i = 0; { 
    // printf("OK %d\n", i);
    k = 0;
    general_param *p = &(GLOBAL.general);
    ADD(ScalarString(mkChar(MODENAMES[p->mode])));
    ADD(ScalarInteger(p->Rprintlevel));    
    ADD(ScalarLogical(p->storing));
    ADD(ScalarLogical(p->skipchecks));
    ADD(ScalarInteger(p->every));    
    ADD(ScalarReal(p->gridtolerance));
    ADDCHAR(p->pch);
    if (p->naturalscaling==0 || p->naturalscaling==1) {
      ADD(ScalarLogical(p->naturalscaling));
    } else {
      ADD(ScalarInteger(p->naturalscaling));
    }
    ADD(ScalarLogical(p->sp_conform));
    ADD(ScalarInteger(p->Cprintlevel));
    ADD(ExtendedBoolean(p->exactness));    

    int nn, n_inv=0;
    for (nn=0; nn<MAXINVERSIONS; nn++) n_inv += p->matrix_inversion[nn] >= 0;
    SET_VECTOR_ELT(sublist[i], k++, Int(p->matrix_inversion, n_inv, n_inv));
    ADD(ScalarReal(p->matrixtolerance));    
    ADD(ScalarLogical(p->allowdist0));   
    ADD(ScalarLogical(p->na_rm_lines));   
    ADD(ScalarLogical(p->vdim_close_together));
    ADD(ScalarInteger(p->expected_number_simu));    
    ADD(ScalarInteger(p->seed));    
    ADD(ScalarLogical(p->detailed_output));   
    ADD(ScalarLogical(p->asList));   
    ADD(ScalarLogical(p->Ttriple == NA_INTEGER ? NA_LOGICAL
		      : p->Ttriple != 0));
    ADD(ScalarLogical(p->returncall));   
    ADD(ScalarString(mkChar(OUTPUTMODENAMES[p->output])));
    ADD(ScalarString(mkChar(REPORTCOORDNAMES[p->reportcoord])));
    ADD(ScalarInteger(p->set + 1));   
 }
  
  //  printf("OK %d\n", i);

  i++; {
    k = 0;
    gauss_param *p = &(GLOBAL.gauss);
    // nachfolgend sollte immer >= 0 sein
    ADD(ScalarLogical(p->paired));   
    ADD(ExtendedBoolean(p->stationary_only));
    ADD(ScalarReal(p->approx_zero));    
    ADD(ScalarInteger(p->direct_bestvariables));    
    ADD(ScalarLogical(p->loggauss));   
    SET_VECTOR_ELT(sublist[i], k++, Num(p->boxcox, 2 * MAXGAUSSVDIM, 
					2 * MAXGAUSSVDIM)); 
  }

  i++; {
    k = 0;
    krige_param *p = &(GLOBAL.krige);
    ADD(ScalarLogical(p->ret_variance));   
    ADD(ScalarInteger(p->locmaxn));
    SET_VECTOR_ELT(sublist[i], k++, Int(p->locsplitn, 3, 3));
    ADD(ScalarInteger(p->locsplitfactor));
    ADD(ScalarLogical(p->fillall));   
  }

  i++; {
    k = 0;
    ce_param *p = &(GLOBAL.ce);
    ADD(ScalarLogical(p->force)); 
    SET_VECTOR_ELT(sublist[i], k++, Num(p->mmin, MAXCEDIM, MAXCEDIM)); 
    ADD(ScalarInteger(p->strategy)); 
    ADD(ScalarReal(p->maxGB)); 
    ADD(ScalarReal(p->maxmem)); 
    ADD(ScalarReal(p->tol_im)); 
    ADD(ScalarReal(p->tol_re));	      
    ADD(ScalarInteger(p->trials));    
    ADD(ScalarLogical(p->useprimes)); 
    ADD(ScalarLogical(p->dependent)); 
    ADD(ScalarReal(p->approx_grid_step)); 
    ADD(ScalarInteger(p->maxgridsize));    
  } 

  i++; {
    k = 0;
    direct_param *p = &(GLOBAL.direct);
    ADD(ScalarInteger(p->inversionmethod));    
    ADD(ScalarReal(p->svdtolerance));     
    ADD(ScalarInteger(p->maxvariables));    
  }

  i++; {
    k = 0;
    nugget_param *p = &(GLOBAL.nugget);
    ADD(ScalarReal(p->tol));
    // ADD(ScalarLogical(p->meth));
  }

  i++; {
    k = 0;
    sequ_param *p = &(GLOBAL.sequ);
    ADD(ScalarInteger(p->max)) /*  does not need Extended here */;
    ADD(ScalarInteger(p->back));    
    ADD(ScalarInteger(p->initial));    
  }

  i++; {
    k = 0;
    spectral_param *p = &(GLOBAL.spectral);
    SET_VECTOR_ELT(sublist[i], k++, Int(p->lines, MAXTBMSPDIM,MAXTBMSPDIM));
    ADD(ScalarLogical(p->grid));
    ADD(ScalarReal(p->prop_factor));
    ADD(ScalarReal(p->sigma));
  }

  i++; {
    k = 0;
    tbm_param *p = &(GLOBAL.tbm);
    // nachfolgend sollte immer >= 0 sein
    ADD(ScalarInteger(p->tbmdim));
    ADD(ScalarInteger(p->fulldim));
    SET_VECTOR_ELT(sublist[i], k++, Num(p->center, MAXTBMSPDIM, MAXTBMSPDIM));
    ADD(ScalarInteger(p->points));     
    // ADD(p->method>=0 ? ScalarString(mkChar(METHODNAMES[p->method])) : 
    //	R_NilValue);
     SET_VECTOR_ELT(sublist[i], k++, Int(p->lines, MAXTBMDIM,MAXTBMDIM));
    ADD(ScalarReal(p->linesimufactor));
    ADD(ScalarReal(p->linesimustep));
    ADD(ExtendedBoolean(p->layers));
    ADD(ScalarLogical(p->grid));
  }


  i++; {
    k = 0;
    mpp_param *p = &(GLOBAL.mpp);
     ADD(ScalarInteger(p->n_estim_E));
    SET_VECTOR_ELT(sublist[i], k++,
		   Num(p->intensity, MAXMPPDIM ,MAXMPPDIM));
    //ADD(ScalarReal(p->refradius_factor));
    ADD(ScalarReal(p->about_zero));
    ADD(ScalarReal(p->shape_power));
    SET_VECTOR_ELT(sublist[i], k++, Num(p->scatter_step, MAXMPPDIM,MAXMPPDIM)); 
    SET_VECTOR_ELT(sublist[i], k++, Int(p->scatter_max, MAXMPPDIM,MAXMPPDIM)); 
     //    SET_VECTOR_ELT(sublist[i], k++, Num(p->plus, MAXMPPDIM ,MAXMPPDIM));
    //    ADD(ScalarReal(p->approxzero));
    //   ADD(ScalarReal(p->samplingdist));
    //  ADD(ScalarReal(p->samplingr));
     //  ADD(ScalarReal(p->p));
  } 

  i++; {
    k = 0;
    hyper_param *p = &(GLOBAL.hyper);
    ADD(ScalarInteger(p->superpos));    
    ADD(ScalarInteger(p->maxlines));    
    ADD(ScalarInteger(p->mar_distr));
    ADD(ScalarReal(p->mar_param));    
  } 

  i++; {
    k = 0;
    extremes_param *p = &(GLOBAL.extreme);
    ADD(ScalarReal(p->standardmax));
    ADD(ScalarInteger(p->maxpoints));           
    ADD(ScalarReal(p->GEV_xi));
    ADD(ScalarReal(p->density_ratio));
    ADD(ScalarInteger(p->check_every));
    ADD(p->flat == 0 || p->flat == 1 ? ScalarLogical(p->flat) :
	ScalarInteger(p->flat));
    ADD(ScalarInteger(p->min_n_zhou));
    ADD(ScalarInteger(p->max_n_zhou));
    ADD(ScalarReal(p->eps_zhou));
    ADD(ScalarInteger(p->mcmc_zhou));
  }

  i++; {
    k = 0;
    br_param *p = &(GLOBAL.br);
    ADD(ScalarInteger(p->BRmaxmem)); 
    ADD(ScalarReal(p->BRmeshsize));
    ADD(ScalarInteger(p->BRvertnumber));
    ADD(ScalarInteger(p->BRoptim));
    ADD(ScalarReal(p->BRoptimtol));
    ADD(ScalarInteger(p->BRoptimmaxpoints));
    ADD(ScalarReal(p->variobound));
    ADD(ScalarInteger(p->deltaAM));
    ADD(ScalarReal(p->corr_factor));
  }

  i++; {
    k = 0;
    distr_param *p = &(GLOBAL.distr);
    ADD(ScalarReal(p->safety));
    ADD(ScalarReal(p->minsteplen));
    ADD(ScalarInteger(p->maxsteps));
    ADD(ScalarInteger(p->parts));
    ADD(ScalarInteger(p->maxit));
    ADD(ScalarReal(p->innermin));
    ADD(ScalarReal(p->outermax));
    ADD(ScalarInteger(p->mcmc_n));
    ADD(ScalarInteger(p->repetitions));
  }
  
  i++; {
    k = 0;
    fit_param *p = &(GLOBAL.fit);
    ADD(ScalarReal(p->bin_dist_factor));
    ADD(ScalarReal(p->upperbound_scale_factor)); 
    ADD(ScalarReal(p->lowerbound_scale_factor)); 
    ADD(ScalarReal(p->lowerbound_scale_LS_factor)); 
    ADD(ScalarReal(p->upperbound_var_factor)); 
    ADD(ScalarReal(p->lowerbound_var_factor));
    //  ADD(ExtendedBoolean(p->lowerbound_sill)); 
    ADD(ScalarReal(p->scale_max_relative_factor)); 
    ADD(ScalarReal(p->minbounddistance));
    ADD(ScalarReal(p->minboundreldist));  
    ADD(ScalarInteger(p->approximate_functioncalls)); 
    ADD(ScalarReal(p->minmixedvar)); 
    ADD(ScalarReal(p->maxmixedvar)); 
    //  ADD(ScalarReal(p->solvesigma));
    SET_VECTOR_ELT(sublist[i], k++, Num(p->BC_lambdaLB, 2 * MAXGAUSSVDIM, 
					2 * MAXGAUSSVDIM));
    SET_VECTOR_ELT(sublist[i], k++, Num(p->BC_lambdaUB, 2 * MAXGAUSSVDIM, 
					2 * MAXGAUSSVDIM));
    ADD(ScalarLogical(p->use_naturalscaling));
    ADD(ScalarInteger(p->bins));
    ADD(ScalarInteger(p->nphi));
    ADD(ScalarInteger(p->ntheta)); 
    ADD(ScalarInteger(p->ntime));
    //   ADD(ScalarReal(p->sill));
    ADD(ScalarLogical(p->onlyuser));
    //    ADD(ScalarString(mkChar(OPTIM_VAR_NAMES[p->optim_var_estim])));
    ADD(ScalarInteger(p->lengthshortname));
    ADD(ScalarInteger(p->split));
    ADD(ScalarReal(p->scale_ratio));
    ADD(ScalarInteger(p->critical));
    ADD(ScalarInteger(p->n_crit));
    ADD(ScalarInteger(p->locmaxn));
    SET_VECTOR_ELT(sublist[i], k++, Int(p->locsplitn,3, 3));
    ADD(ScalarInteger(p->locsplitfactor));
    ADD(ScalarInteger(p->smalldataset));
    ADD(ScalarReal(p->min_diag));
    ADD(ScalarLogical(p->reoptimise));
    ADD(p->optimiser>=0 ? ScalarString(mkChar(OPTIMISER_NAMES[p->optimiser])) 
	: R_NilValue);
    ADD(p->algorithm < 0 ? R_NilValue :
	ScalarString(mkChar(p->optimiser == 3 ? NLOPTR_NAMES[p->algorithm]
			    : "")
		     ));
    ADD(ScalarString(mkChar(LIKELIHOOD_NAMES[p->likelihood])));
    ADD(ScalarLogical(p->ratiotest_approx));
    ADD(ScalarLogical(p->split_refined));
    ADD(ScalarLogical(p->cross_refit));
    ADD(ScalarInteger(p->estimate_variance));
    ADD(ScalarReal(p->pgtol)); 
    ADD(ScalarReal(p->pgtol_recall)); 
    ADD(ScalarReal(p->factr)); 
    ADD(ScalarReal(p->factr_recall)); 
   } 

  i++; {
    k = 0;
    empvario_param *p = &(GLOBAL.empvario);   
    ADD(ScalarReal(p->phi0));
    ADD(ScalarReal(p->theta0));
    ADD(ScalarReal(p->tol));
    ADD(ScalarLogical(p->pseudovariogram));
    ADD(ScalarLogical(p->fft));
  }

  i++; {
    k = 0;
    gui_param *p = &(GLOBAL.gui);   
    ADD(ScalarLogical(p->alwaysSimulate));
    ADD(p->method>=0 ? ScalarString(mkChar(METHODNAMES[p->method])) 
	: R_NilValue);
    SET_VECTOR_ELT(sublist[i], k++, Int(p->size, 2, 2));
  }

  i++; {
    k = 0;
    graphics_param *p = &(GLOBAL.graphics);   
    ADD(ScalarLogical(p->always_close));
    ADD(ScalarInteger(p->PL));
    ADD(ScalarReal(p->height));
    SET_VECTOR_ELT(sublist[i], k++, Int(p->increase_upto, 2, 2));
    ADD(ScalarLogical(p->always_open == NA_INTEGER ? NA_LOGICAL
		      : p->always_open != 0));
    ADD(ScalarString(mkChar(p->filename)));
    ADD(ScalarLogical(p->onefile));
    ADD(ScalarInteger(p->number));
    ADD(ScalarReal(p->resolution));
   }

  i++; {
    k = 0;
    registers_param *p = &(GLOBAL.registers);
    ADD(ScalarInteger(p->keynr));    
    ADD(ScalarInteger(p->predict));    
    ADD(ScalarInteger(p->likelihood));    
  }
  
  i++; {
    k = 0;
    internal_param *p = &(GLOBAL.internal);
    ADD(ScalarLogical(p->warn_oldstyle));
    ADD(ScalarLogical(p->warn_newstyle));
    ADD(ScalarLogical(p->warn_Aniso));
    ADD(ScalarLogical(p->warn_ambiguous));
    ADD(ScalarLogical(p->warn_normal_mode));
    ADD(ScalarLogical(p->warn_mode));
    ADD(ScalarLogical(p->stored_init));
    ADD(ScalarLogical(p->warn_scale));
    ADD(ScalarLogical(p->warn_coordinates));
    ADD(ScalarLogical(p->warn_on_grid));
    ADD(ScalarLogical(p->warn_new_definitions));
    ADD(ScalarLogical(p->warn_aspect_ratio));
    ADD(ScalarLogical(p->warn_coord_change));
    ADD(ScalarLogical(p->warn_color_palette));
    ADD(ScalarLogical(p->warn_zenit));
    ADD(ScalarLogical(p->do_tests));
    ADD(ScalarLogical(p->warn_constant));
    ADD(ScalarLogical(p->warn_negvar));
    ADD(ScalarLogical(p->warn_onlyvar));
    ADD(ScalarInteger(p->examples_reduced));
 }


  i++; {
    k = 0;
    coords_param *p = &(GLOBAL.coords);
    ADD(ExtendedInteger(p->xyz_notation));    
    ADD(ScalarString(mkChar(COORD_SYS_NAMES[p->coord_system])));
    ADD(UNITS(p->newunits));
    ADD(UNITS(p->curunits));
    ADD(UNITS(p->varunits));
    if (p->data_nr_names == 0) {
      SET_VECTOR_ELT(sublist[i], k++, Int(p->data_idx, 2, 2)); 
    } else {
      SET_VECTOR_ELT(sublist[i], k++,
		     String(p->data_names, p->data_nr_names, 
			  p->data_nr_names));       
    }
    if (p->x_nr_names == 0) {
      SET_VECTOR_ELT(sublist[i], k++, Int(p->x_idx, 2, 2)); 
    } else {
      SET_VECTOR_ELT(sublist[i], k++,
		     String(p->x_names, p->x_nr_names, 
			  p->x_nr_names));	
    }
    ADD(ScalarString(mkChar(COORD_SYS_NAMES[p->new_coord_system])));
    SET_VECTOR_ELT(sublist[i], k++, Num(p->zenit, 2, 2));     
    ADD(ScalarLogical(p->polar_coord));
   }

  i++; {
    k = 0;
    special_param *p = &(GLOBAL.special);
    ADD(ScalarInteger(p->multcopies));    
  }

  i++; {
    k = 0;
    solve_param *p = &(GLOBAL.solve);
    ADD(ExtendedBoolean(p->sparse));
    ADD(ScalarReal(p->spam_tol));    
    ADD(ScalarReal(p->spam_min_p));    
    ADD(ScalarReal(p->svd_tol));    
    SET_VECTOR_ELT(sublist[i], k++, Int(p->Methods, SOLVE_METHODS, 
					SOLVE_METHODS)); 
    ADD(ScalarInteger(p->spam_min_n));    
    ADD(ScalarInteger(p->spam_sample_n));    
    ADD(ScalarInteger(p->spam_factor));    
    ADD(ScalarInteger(p->pivot));    
  }
 
  //   print("%d %d\n", i, prefixN -1);
  assert (i == trueprefixN - 1); 
    /*

    ADD(ScalarReal(p->));
    ADD(ScalarReal(p->));
 
    ADD(ScalarInteger(p->));    
    ADD(ScalarInteger(p->));    

    ADD(ScalarLogical(p->));
    ADD(ScalarLogical(p->));

    ADD(ScalarString(p->));
    ADD(ScalarString(p->));

      */

  for (i=0; i<trueprefixN; i++) {
    setAttrib(sublist[i], R_NamesSymbol, subnames[i]);
    SET_VECTOR_ELT(list, i, sublist[i]);
  }
  setAttrib(list, R_NamesSymbol, names);
  
  UNPROTECT(2 + 2 * trueprefixN);

  return list;

}




void splitAndSet(SEXP el, char *name, bool isList) {
  int i, len;
  char msg[200];
  char prefix[200], mainname[200];   
  //  printf("splitandset\n");
  len = strlen(name);
  for (i=0; i < len && name[i]!='.'; i++);
  sprintf(msg, "argument '%s' not valid\n", name);
  if (i==0) ERR(msg);
  if (i==len) {
    strcpy(prefix, "");
    strcopyN(mainname, name, 200);
  } else {
    strcopyN(prefix, name, i);
    prefix[i] = '\0';
    strcpy(mainname, name+i+1);
  }

  //  printf("i=%d %d\n", i, len);
  setparameter(el, prefix, mainname, isList);
  //   printf("ende\n");
}


int InternalGetProcessType(cov_model *cov) {
  int nr = cov->nr;
  if (isInterface(cov)) return InternalGetProcessType(cov->sub[0]);
  
  switch(CovList[nr].Typi[0]) {
  case TcfType: case PosDefType: case VariogramType : case TrendType :
  case GaussMethodType : return GAUSSPROC;
  case BrMethodType : return BROWNRESNICKPROC;
  case ProcessType :
    if (nr == DOLLAR_PROC) return InternalGetProcessType(cov->sub[0]);
    else if (nr == PLUS_PROC || nr == MULT_PROC //|| nr == TREND_PROC
	     ) 
      return GAUSSPROC;
    return cov->nr;
  case UndefinedType:
    if (nr == PLUS || nr == MULT || nr ==DOLLAR || nr == POWER_DOLLAR || 
	nr == USER) return GAUSSPROC;
    else BUG;
  default :
    BUG;
  }
 BUG;
}
 
cov_model *Build_cov(SEXP model_reg, SEXP model) {
  if(currentNrCov < 0) InitModelList();
  int nr = INTEGER(model_reg)[0];
  if (nr < 0 || nr >= MODEL_MAX) BUG;
  cov_model **Cov = KEY + nr;
  if (*Cov != NULL) COV_DELETE(Cov);
  assert(*Cov == NULL);
  CMbuild(model, 0, Cov);
  return *Cov;
}

SEXP GetProcessType(SEXP model_reg, SEXP model) {
  cov_model *cov = Build_cov(model_reg, model);
  int covnr = InternalGetProcessType(cov);
  SEXP ans;
  PROTECT (ans = allocVector(STRSXP, 1));
  SET_STRING_ELT(ans, 0, mkChar(CovList[covnr].nick));
  UNPROTECT(1);
  return ans;
}





SEXP RFoptions(SEXP options) {
  int i, j, lenlist, lensub;
  SEXP el, list, sublist, names, subnames;
  char *name, *pref;
  bool isList = false;
 /* 
     In case of strange values of a parameter, undelete
     the comment for PRINTF
  */

  
  // PRINTF("start %f\n", GLOBAL.gauss.exactness);
  options = CDR(options); /* skip 'name' */
  if (options == R_NilValue) {
    //PRINTF("before get %f\n", 1.);
    list  = getRFoptions(); 
    //   PRINTF("after get %f\n", 1);  
    return list;
  }

  name = (char*) (isNull(TAG(options)) ? "" : CHAR(PRINTNAME(TAG(options))));
  if ((isList = strcmp(name, "LIST")==0)) {   
    list = CAR(options);
    if (TYPEOF(list) != VECSXP)
      ERR1("'LIST' needs as argument the output of '%s'", RFOPTIONS);
    names = getAttrib(list, R_NamesSymbol);   
    lenlist = length(list);
    for (i=0; i<lenlist; i++) {
      int len;
      pref = (char*) CHAR(STRING_ELT(names, i));  

      //   print("%d %s warn.ambig=%d\n", i, pref, GLOBAL.warn.ambiguous);

      sublist = VECTOR_ELT(list, i);
      len = strlen(pref);
      for (j=0; j < len && pref[j]!='.'; j++);
      if (TYPEOF(sublist) == VECSXP && j==len) { // no "."
	// so, general parameters may not be lists,
	// others yes
	lensub = length(sublist);
	subnames = getAttrib(sublist, R_NamesSymbol); 
	for (j=0; j<lensub; j++) {
	  name = (char*) CHAR(STRING_ELT(subnames, j));

	  // print("    %d %s warn.ambig=%d\n", j, name, GLOBAL.warn.ambiguous);

	  //    print("%d %d %s : %f %f\n", i, j, name,
	  //	     GLOBAL.gauss.exactness, GLOBAL.TBM.linesimustep);	  
	  //
	  //print("  %d %d pref=%s name=%s\n", i, j, pref, name);
	  	  
	  setparameter(VECTOR_ELT(sublist, j), pref, name, isList);
	}
      } else {   
	splitAndSet(sublist, pref, isList);
      }
    }
    //    print("end1 %f\n", GLOBAL.TBM.linesimufactor);
  } else {
    for(i = 0; options != R_NilValue; i++, options = CDR(options)) {

      //printf("set opt i=%d\n", i);

      el = CAR(options);
      name = (char*) (isNull(TAG(options)) ? "" :CHAR(PRINTNAME(TAG(options))));
      splitAndSet(el, name, isList);
    }
    //       print("end2 %f\n", GLOBAL.gauss.exactness);
  }
  GLOBAL.general.asList = true;
  return(R_NilValue);
} 
 

void GetModelRegister(char **name, int* nr) {
  *nr = Match(*name, REGNAMES, MODEL_MAX+1);

  //  print("%d\n", *nr);

  if (*nr<0 || *nr > MODEL_MAX) ERR("name for model register unknown");
}
  

void MultiDimRange(int *model_nr, int *set, double *natscale) { 
  MultiDimRange(*set, KEY[*model_nr], natscale); 
}

SEXP countelements(SEXP Idx, SEXP N, SEXP Totparts) {
  int i, *boxes,
    *idx = INTEGER(Idx),
    nbox = INTEGER(Totparts)[0],
    n = INTEGER(N)[0]; 
  SEXP Boxes;

  PROTECT(Boxes = allocVector(INTSXP, nbox));
  boxes = INTEGER(Boxes);
  for (i=0; i<nbox; i++) boxes[i]= 0;

  for (i=0; i<n; i++) {
    boxes[idx[i]]++;
  }

  UNPROTECT(1);
  return Boxes;
}


SEXP countneighbours(SEXP Xdim, SEXP Parts, SEXP Squarelength, SEXP Cumgridlen,
		     SEXP Boxes) {
  int d,  totcumlen, relstart,  y,
    nb[MAXGETNATSCALE], loc[MAXGETNATSCALE], e,
    x = 0,
    sum = 0,
    sl = INTEGER(Squarelength)[0],
    dim = INTEGER(Xdim)[0],
    boundary = (sl - 1) / 2,
    //   total = cumgridlen[dim],
    maxn = GLOBAL.krige.locmaxn,
    *parts = INTEGER(Parts),
    *cumgridlen = INTEGER(Cumgridlen),
    nboxes = length(Boxes),
    *boxes = INTEGER(Boxes);
  assert(dim <= MAXGETNATSCALE); 
  SEXP Neighbours;
  PROTECT(Neighbours = allocVector(INTSXP, nboxes));
  int *neighbours = INTEGER(Neighbours);
 
  totcumlen = 0;
  for (d=0; d<dim; d++) {
    loc[d] = -boundary; 
    nb[d] = 0;
    totcumlen += cumgridlen[d];
    //print("%d %d %d\n", d, totcumlen, cumgridlen[d]);
  }

  relstart = totcumlen * boundary;

  d = 0;
  while(d < dim) {
    y = x - relstart;
    //    print("%d\n", x);
    neighbours[x] = 0; 
    sum = e = 0;
    while(e < dim) {
      bool inside = true;
      int j;
      for (j=0; j<dim; j++) {
	double abs = loc[j] + nb[j];
	if (abs < 0 || abs>=parts[j]) {inside=false; break;}
      }
      if (inside) {
	//     	print("inside %d (%d, %d)\n", y, loc[0], loc[1]);
	sum += boxes[y];
	neighbours[x]++;
      }
      e = 0;	
      loc[e]++;
      y++;
      while (loc[e] > boundary) {
	loc[e] = -boundary; 
	y -= cumgridlen[e] * sl;
	if (++e >= dim) break;
	loc[e]++;
	y += cumgridlen[e];
      }
      //print("e=%d\n", e);
    }
    //     print("sum=%d maxn=%d (%d %d) %d %d parts=%d  %d; %d %dl nei=%d\n",
    // 	    sum, maxn, nb[0], nb[1], totcumlen, relstart, parts[0], parts[1],
	      // 	   cumgridlen[0], cumgridlen[1], neighbours[x]);
     //     assert(false);
    if (sum > maxn) {
      UNPROTECT(1);
      return NILSXP;
    }

    d = 0;					
    nb[d]++;
    x++;
    while (nb[d] >= parts[d]) {
      nb[d] = 0;
      if (++d >= dim) break;
      nb[d]++;
    }
    // print("d=%d (%d %d %d) [%d %d %d]\n", d, nb[0], nb[1], nb[2], 
    //	   parts[0], parts[1], parts[2]);
    //  assert(false);
  }

  UNPROTECT(1);
  return Neighbours;
}


SEXP getelements(SEXP Idx, SEXP Xdim, SEXP N, SEXP Cumgridlen, SEXP Boxes) {
  int i, err = NOERROR, 
    *idx = INTEGER(Idx),
    dim = INTEGER(Xdim)[0],
    n = INTEGER(N)[0],
    *cumgridlen = INTEGER(Cumgridlen),
    *boxes = INTEGER(Boxes),
    total = cumgridlen[dim],
    *count = NULL,
    **elm = NULL;
  SEXP subel = R_NilValue;
   
  if ((elm = (int **) MALLOC(sizeof(int*) * total)) == NULL ||
      (count = (int*) MALLOC(sizeof(int) * total)) == NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  for (i=0; i<total; i++) elm[i] = NULL;
 
  for (i=0; i<total; i++) {
    //   print("%d %d l=%d \n", i, total, boxes[i]); 
    if ((elm[i] = (int *) MALLOC(sizeof(int) * boxes[i])) == NULL) {
      err = ERRORMEMORYALLOCATION;
      goto ErrorHandling;
    }
    count[i] = 0;
  }

  for (i=0; i<n; i++) {
    int k=idx[i];
    elm[k][(count[k])++] = i + 1;
  }

  PROTECT(subel = allocVector(VECSXP, total));
  
  for (i=0; i<total; i++) {
    //  print("%d %d %d\n", i, total, boxes[i]);
    SEXP el;
    int k,
      end = boxes[i];
   
    PROTECT(el = allocVector(INTSXP, boxes[i]));

    for (k=0; k<end; k++) INTEGER(el)[k] = elm[i][k];
    
    SET_VECTOR_ELT(subel, i, el);
    UNPROTECT(1);
  }

  UNPROTECT(1);
 
 ErrorHandling :
  if (elm != NULL) {
    for (i=0; i<total; i++) FREE(elm[i]);
    UNCONDFREE(elm);
  }
  FREE(count);
 
  if (err!=NOERROR) XERR(err);

  return subel; 
}



SEXP getneighbours(SEXP Xdim, SEXP Parts, SEXP Squarelength, 
		   SEXP Cumgridlen, SEXP Neighbours) {
  int i, d, sum, totcumlen, relstart, x, y,
    err = NOERROR,
    nb[MAXGETNATSCALE], loc[MAXGETNATSCALE], e,
    dim = INTEGER(Xdim)[0],
    *parts = INTEGER(Parts),
    sl = INTEGER(Squarelength)[0],
    *cumgridlen = INTEGER(Cumgridlen),
    *neighbours = INTEGER(Neighbours),
    boundary = (sl - 1) / 2,
    total = cumgridlen[dim],
    ** neighb = NULL;
  SEXP subnei = R_NilValue;

  if ( (neighb = (int **) MALLOC(sizeof(int*) * total)) == NULL) {
    err =  ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  // print("OK %d\n", total);

  for (i=0; i<total; i++) neighb[i] = NULL;
 
  for (i=0; i<total; i++) {
    //    print("%d %d l=%d %d\n", i, total, neighbours[i]); 
    // R_CheckUserInterrupt();
    if ((neighb[i] = (int *) MALLOC(sizeof(int) * neighbours[i])) == NULL) {
      err =  ERRORMEMORYALLOCATION;
      goto ErrorHandling;
    }
  }

  x = 0;
  totcumlen = 0;
  for (d=0; d<dim; d++) {
    loc[d] = -boundary; 
    nb[d] = 0;
    totcumlen += cumgridlen[d];
  }

  relstart = totcumlen * boundary;

  d = 0;
  while(d < dim) {
    y = x - relstart;
    sum = e = 0;
    while(e < dim) {
      bool inside = true;
      int j;
      for (j=0; j<dim; j++) {
	double abs = loc[j] + nb[j];
	if (abs < 0 || abs>=parts[j]) {inside=false; break;}
      }
      if (inside) {
	neighb[x][sum] = y + 1; 
	sum++;
      }

      e = 0;	
      loc[e]++;
      y++;
      while (loc[e] > boundary) {
	loc[e] = -boundary; 
	y -= cumgridlen[e] * sl;
	if (++e >= dim) break;
	loc[e]++;
	y += cumgridlen[e];
      }
    }

    d = 0;					
    nb[d]++;
    x++;
    while (nb[d] >= parts[d]) {
      nb[d] = 0;
      if (++d >= dim) break;
      nb[d]++;
    }
  }

  // print("heres\n");

  PROTECT(subnei = allocVector(VECSXP, total));
  
  for (i=0; i<total; i++) {
    //    print("%d %d %d %d\n", i, total, neighbours[i]);
    R_CheckUserInterrupt();
    SEXP nei;
    PROTECT(nei = allocVector(INTSXP, neighbours[i]));

    int k,
      end = neighbours[i];
    for (k=0; k<end; k++) INTEGER(nei)[k] = neighb[i][k];
    
    SET_VECTOR_ELT(subnei, i, nei);
    UNPROTECT(1);
  }

  UNPROTECT(1);
 

 ErrorHandling :
  if (neighb != NULL) {
    for (i=0; i<total; i++) FREE(neighb[i]);
    UNCONDFREE(neighb);
  }
  if (err!=NOERROR) XERR(err);

  return subnei;
}


void DeleteKey(int *reg) {
  COV_DELETE(KEY + *reg);
}


/*
void isAuthor(int *is) {
#ifdef WIN32
  *is = false;
#else
  #define NCHAR 5
  char host[5];
  gethostname(host, NCHAR);
  host[NCHAR-1] = '\0';
  *is = strcmp("viti", host) == 0;
#endif
}
*/


SEXP allintparam() {
  cov_fct *C;
  int n, i, np;
  for (np=n=0; n<currentNrCov; n++) {
    C = CovList + n;
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
    C = CovList + n;
    for (i=0; i<C->kappas; i++) {
      if (C->kappatype[i] == INTSXP) 
	SET_STRING_ELT(x, np++, mkChar(C->kappanames[i]));
    }
  }
  UNPROTECT(1);
  return x;
}


