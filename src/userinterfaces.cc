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
#include <unistd.h>
 
#include <string.h>
#include "RF.h"
#include "primitive.h"
// #include "Covariance.h"


#define nOptimVar 4
const char * OPTIM_VAR_NAMES[nOptimVar] = 
  {"never", "respect bounds", "try", "yes"}; // keep yes last !!

#define nOptimiser 8
const char * OPTIMISER_NAMES[nOptimiser] = 
  {"optim", "optimx", "soma", "nloptr", "GenSA", "minqa", "pso", "DEoptim"
  }; 


#define nNLOPTR 15
const char *NLOPTR_NAMES[nNLOPTR] = 
  //                     Zeiten fuer ein Bsp; besser als optim?
  // optim 3.6
  // optimx 3.5
  {"NLOPT_GN_DIRECT", //  4.1 ; +- ; laufen nicht an den Rand heran!!
   "NLOPT_GN_DIRECT_L", //
   "NLOPT_GN_DIRECT_L_RAND", // 
   "NLOPT_GN_DIRECT_NOSCAL", //
   "NLOPT_GN_DIRECT_L_NOSCAL", //
   "NLOPT_GN_DIRECT_L_RAND_NOSCAL", //
   "NLOPT_GN_ORIG_DIRECT", // 3.
   "NLOPT_GN_ORIG_DIRECT_L", //
   //"NLOPT_GD_STOGO", // benoetigt gradienten (b. g.)
   // "NLOPT_GD_STOGO_RAND", // b. g.
   // "NLOPT_LD_SLSQP", //   b. g.
   //"NLOPT_LD_LBFGS_NOCEDAL", // b. g.
   // "NLOPT_LD_LBFGS", //   b. g.
   "NLOPT_LN_PRAXIS", // 3.1; eher schlechter als optim
   // "NLOPT_LD_VAR1", // b. g.
   // "NLOPT_LD_VAR2", //  b. g.
   // "NLOPT_LD_TNEWTON", // b. g.
   // "NLOPT_LD_TNEWTON_RESTART", //  b. g.
   //"NLOPT_LD_TNEWTON_PRECOND", // b. g.
   //"NLOPT_LD_TNEWTON_PRECOND_RESTART", // b. g.
   "NLOPT_GN_CRS2_LM",// 3.8; gut; viel naeher am Rand als die obigen _GN_
   // "NLOPT_GN_MLSL", //  b. g.
   //"NLOPT_GD_MLSL", // b. g.
   //"NLOPT_GN_MLSL_LDS", //  b. g.
   // "NLOPT_GD_MLSL_LDS", // b. g.
   //"NLOPT_LD_MMA", //  b. g.
   "NLOPT_LN_COBYLA", // 1.9; laeuft vollstaendig ran
   //  "NLOPT_LN_NEWUOA", // funktioniert nicht. Grund nicht nachgeschaut
   //"NLOPT_LN_NEWUOA_BOUND", // funktioniert nicht. Grund nicht nachgeschaut
   "NLOPT_LN_NELDERMEAD", // 0.8 !!  laeuft vollstaendig ran
   "NLOPT_LN_SBPLX", // 2.6;   laeuft vollstaendig ran
   // "NLOPT_LN_AUGLAG", // b. g.
   // "NLOPT_LD_AUGLAG", // b. g.
   //   "NLOPT_LN_AUGLAG_EQ", //  b. g.
   //   "NLOPT_LD_AUGLAG_EQ", // b. g.
   "NLOPT_LN_BOBYQA", // 1.8 !! Und ist genauso gut wie optim im Bsp!!
   "NLOPT_GN_ISRES" // aehnlich NLOPT_GN_CRS2_LM
  }; 

void ResetWarnings() {
  warn_param *w = &(GLOBAL.warn);
  w->oldstyle = w->newstyle = w->Aniso = w->ambiguous = w->normal_mode =
    w->warn_mode = w->warn_colour = w->warn_coordinates = true;
}

#define MaxMaxInts 9
void GetMaxDims(int *maxints) {  *maxints = MaxMaxInts; }

void GetrfParameters(int *covmaxchar, int *methodmaxchar, 
		     int *typemaxchar,
		     int *covnr, int *methodnr, int *typenr,
		     int *maxdim,
		     int *maxmodels) {
  // if (currentNrCov==-1) InitModelList();
  int i=0;
  *covmaxchar=MAXCHAR;
  *methodmaxchar=METHODMAXCHAR;
  *typemaxchar=-1; // obsolete
  *covnr=currentNrCov;
  *methodnr=(int) Forbidden;
  *typenr=ROLE_LAST;
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
  *maxmodels=MAXFIELDS;
}

void GetrfParametersI(int *covmaxchar, int *methodmaxchar, int *typemaxchar,
		      int *covnr, int *methodnr, int *typenr,
		      int *maxdim, int *maxmodels){
  if (currentNrCov==-1) InitModelList();
  GetrfParameters(covmaxchar, methodmaxchar, typemaxchar,
		  covnr, methodnr, typenr, maxdim, maxmodels);
}


bool skipchecks[nr_modes] =  {true, false, false, false, false, false, false},
  allowdistance0[nr_modes] = {true, false, false, false, false, false, false},
  ce_force[nr_modes] =       {true, true, true, false, false, false, false},
  ce_dependent[nr_modes] =   {true, true, true, false, false, false, false},
  grid[nr_modes] =           {true, true, true, true, true, false, false},
  fit_split[nr_modes] =      {false,  true,  true, true, true, true, true},
  fit_refine[nr_modes] =     {false, false, false, true, true, true, true},
  fit_reoptimise[nr_modes] = {false, false, false, true, true, true, true}
  ;
char pch[nr_modes] =         {'\0',  '\0',  '\0',  '.',   '.',   '.',   '.'}
  ;
int locmaxn[nr_modes] =      {3000, 4000, 6000, 8000, 9000, 10000000, 10000000},
  ce_trials[nr_modes] =      {1,    1,    2,    3,    4,    5,   6},
  spectral_lines[nr_modes] = {300, 500, 1500,  2500, 3500, 5000, 10000},
  tbm_lines[nr_modes] =      {40,   40,   50,   60,   80,   100,   200},
  mpp_n_estim_E[nr_modes] =  {10000, 10000, 20000,50000,80000,100000,1000000},
  hyper_superpos[nr_modes] = {200, 300, 450, 700, 1000, 2000, 5000},
  fit_critical[nr_modes] =   {-1,  -1,   -1,   0,   1,   2,   3}, 
  fit_ncrit[nr_modes] =      { 2,   4,    5,   5,   10,  10,  20}, 
  fit_optim_var[nr_modes] =  { 2,   2,    2,   2,   1,   1,   1} 
  ;
double exactness[nr_modes] = {false, false, false, NA_REAL, true, true, true},
  matrixtolerance[nr_modes] ={1e-6,  1e-6,  1e-6,  1e-12,   1e-14, 0, 0},
  ce_tolIm[nr_modes] =       {1e6,  1e6,  1e-1,  1e-3,   1e-7, 0, 0},
  ce_tolRe[nr_modes] =       {-1e6,  -1e6,  -1e-3, -1e-7,  -1e-14, 0, 0},
  ce_approx_step[nr_modes] = {1.0,  1.0,   1.0, 1.0, 0.0, 0.0, 0.0},
  direct_tol[nr_modes] =     {1e-8, 1e-8, 1e-10, 1e-12, 1e-14, 0, 0},
  nugget_tol[nr_modes] =     {1e-8, 1e-8, 1e-12, 0,     0,     0,     0},
  tbm_linefactor[nr_modes] = {1.5,  1.5,  1.7,   2.0,  3.0,  5.0,  6.0},
  mpp_intensity[nr_modes] =  {50, 50, 80, 100, 200, 500, 1000},
  mpp_zero[nr_modes] =       {1e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7},
  max_max_gauss[nr_modes] =  {2,    2,    3,    4,    5,   6,   6}
  ;

const char *f_opt[nr_modes] = {"optim", "optim", "optim", "optim", "optimx", "optimx", "optimx"}; // to do optimx
 
void SetDefaultModeValues(int old, int m){
  int i;
  //                      high  fast, normal, save, pedantic
  GLOBAL.general.skipchecks = skipchecks[m];
  GLOBAL.general.pch = pch[m];
  GLOBAL.general.exactness = exactness[m];
  GLOBAL.general.matrixtolerance = matrixtolerance[m];
  GLOBAL.general.allowdist0 = allowdistance0[m];
  GLOBAL.krige.locmaxn = locmaxn[m];
  GLOBAL.krige.locsplitn[0]  = (locmaxn[m] * 5) / 8;
  GLOBAL.krige.locsplitn[1]  = locmaxn[m] / 40;
  GLOBAL.krige.locsplitn[2]  = locmaxn[m] / 8;
  GLOBAL.ce.force = ce_force[m];
  GLOBAL.ce.tol_im = ce_tolIm[m];
  GLOBAL.ce.tol_re = ce_tolRe[m];
  GLOBAL.ce.trials = ce_trials[m];
  GLOBAL.ce.dependent = ce_dependent[m];
  GLOBAL.ce.approx_grid_step = ce_approx_step[m];
  GLOBAL.direct.svdtolerance = direct_tol[m];
  GLOBAL.nugget.tol = nugget_tol[m];
  for(i=0; i<4; GLOBAL.spectral.lines[i++] = spectral_lines[m]);
  GLOBAL.spectral.grid = grid[m];
  GLOBAL.tbm.lines[1] = tbm_lines[m];
  GLOBAL.tbm.lines[2] = (tbm_lines[m] * 25) / 6;
  GLOBAL.tbm.linesimufactor = tbm_linefactor[m];
  GLOBAL.tbm.grid = grid[m];
  GLOBAL.mpp.n_estim_E = mpp_n_estim_E[m];
  for(i=0; i<4; GLOBAL.mpp.intensity[i++] = mpp_intensity[m]);
  GLOBAL.mpp.about_zero = mpp_zero[m];
  GLOBAL.hyper.superpos = hyper_superpos[m];
  GLOBAL.extreme.standardmax = max_max_gauss[m];
  fit_param *f = &(GLOBAL.fit);
  f->split = fit_split[m];
  f->refine_onborder = fit_refine[m];
  f->reoptimise = fit_reoptimise[m];
  f->critical = fit_critical[m];
  f->n_crit = fit_ncrit[m];
  f->optim_var_estim = fit_optim_var[m];
  char dummy[10];
  strcpy(dummy, f_opt[m]);
  f->optimiser = Match(dummy, OPTIMISER_NAMES, nOptimiser);

  //  printf("optimiser %d %s\n", f->optimiser,  OPTIMISER_NAMES[f->optimiser]);

  warn_param *w = &(GLOBAL.warn);
  w->stored_init = false;
  if (m < normal) w->oldstyle = w->newstyle = w->Aniso = w->ambiguous = false;
  else if (m>old) w->oldstyle = w->Aniso = w->ambiguous = true;
  if (m != normal && w->warn_mode) {
    PRINTF("Note that the option 'mode' is still in an experimental stage, so that the behaviour may change (slightly) in future.");
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


void GetMethodNr(char **name, int *nr) {
  // == -1 if no matching method is found
  // == -2 if multiple matching methods are found, without one matching exactly
  *nr = Match(*name, METHODNAMES, (int) Nothing);
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
  if (i==1+(int) Nothing) PRINTF("  * no methods implemented yet\n");
  PRINTF("\n\n  == end of method list ================\n\n");
  PRINTF("\n");
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
  for (i=0; i<*lrange; i++) range[i]=RF_NAN;
 *index = -100;
 return;
*/
}


void PMLheader(char* firstcolumn, int nick) {
  int i;
  char header1[]=" #    cir cut int TBM spe dir seq ave coi hyp spe\n";
  char header2[]=" p    cul off rin     ctr ect uen rag ns  erp cif\n";
  for (i=0; i<=nick; i++) PRINTF(firstcolumn, ""); 
  PRINTF("%4s", ""); PRINTF(header1);  
   for (i=0; i<=nick; i++) PRINTF(firstcolumn, ""); 
  PRINTF("%4s", ""); PRINTF(header2);
}

void PrintModelList(int *intern, int *operat, int* Nick)
{
    int i, k, last_method, m, OP;
//char header[]="circ cut intr TBM spec dir seq Mark ave add hyp part\n";
  char coded[6][2]={"-", "X", "+", "N", "H", "S"};
//  char typenames[4][2]={"i", "f", "s", "n"};
  char specialnames[4][2]={".", "n", "f", "?"};
  char firstcolumn[20], name[MAXCHAR];
  int maxchar=10; // <= MAXCHAR=14
  int type[MAXNRCOVFCTS], op[MAXNRCOVFCTS], monotone[MAXNRCOVFCTS], 
    finite[MAXNRCOVFCTS], internal[MAXNRCOVFCTS], dom[MAXNRCOVFCTS], 
    iso[MAXNRCOVFCTS], vdim[MAXNRCOVFCTS], maxdim[MAXNRCOVFCTS],
    nick = *Nick;
  
  last_method = (int) Nothing; // even not special method   
  if (currentNrCov==-1) { 
     InitModelList(); 
   // if (PL>5)  PRINTF("List of covariance functions initiated.\n"); 
  }
  if (CovList==NULL) {PRINTF("There are no functions available!\n");} 
  else {
    cov_fct *C;

    GetAttr(type, op, monotone, finite, internal, dom, iso, maxdim, vdim);
 
    sprintf(firstcolumn,"%%%ds", -maxchar);
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
      PMLheader(firstcolumn, nick);

      for (k=1, i=0; i<currentNrCov; i++, C++) { 
	if (!isPosDef((Types)(type[i])) && !isUndefinedType((Types)(type[i])))
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
			    : isUndefinedType((Types)(type[i])) || 
			    monotone[i]<0 || finite[i] < 0 ? 3
			    : 0]	                          
	       );	  
	PRINTF(" ");
	//             above works since normal mixture cannot have finite.range 
	assert(internal[i]==0 || internal[i]==1);
	
	for (m=(int) CircEmbed; m<last_method; m++)
	  if (m != Nugget && m != Markov)
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





/*
void GetKeyInfo(int *keyNr, int *total, int *lengths, int *spatialdim, 
		int *timespacedim, 
		int *grid, int *role, int *maxdim, int *vdim)
{
    // check with DoSimulateRF and subsequently with extremes.cc, l.170
    // if any changings !
  int d;
  cov_model *key;
  simu_type *simu;
  location_type *loc;

  if (*maxdim<MAXSIMUDIM) { *total=-3; *maxdim=MAXSIMUDIM; return;} 
  if ((*keyNr<0) || (*keyNr>MODEL_MAX)) {*total=-1; return;} 
  if ((key = KEY[*keyNr]) == NULL) {*total=-3; return;}  
  simu = &(key->simu);
  loc = key->prevloc;

  if (!simu->active) {*total=-2;} 
  else { 
    *total=loc->totalpoints;
    *spatialdim = loc->spatialdim;
    *timespacedim = loc->timespacedim;
    for (d=0; d<loc->timespacedim; d++) lengths[d]=loc->length[d];
    *grid = (int) loc->grid;
    *role = (int) key->role; 
    *maxdim = MAXSIMUDIM;
    *vdim = key->vdim;
  }
}
*/


void GetAttr(int *type, int *op, int *monotone, int *finiterange, 
	     int *internal, int *dom, int *iso, int *maxdim, int *vdim) {
#define MAXPN 10 /* only used for testing purposes */
  int nr, p;
  cov_model Cov,
    *cov = &Cov;
  range_type range;
//  range_type *range;
  cov_fct *C = CovList;


  for (p=0; p<C->kappas; p++) 
    cov->px[p] = (double*) CALLOC(MAXPN, sizeof(double));
  Cov.tsdim = 1;
  Cov.vdim = 1;
  Cov.nsub = 2;
 
  for (nr=0; nr<currentNrCov; nr++, C++){   
    Cov.nr = nr;
    type[nr] = C->Type;
    op[nr] = (int) C->maxsub > 0;   
    C->range(cov, &range);
    maxdim[nr] = C->maxdim;
    finiterange[nr] = C->finiterange;
    monotone[nr] = C->Monotone;
    internal[nr] = C->internal;
    dom[nr] = C->domain;
    iso[nr] = C->isotropy;
    vdim[nr] = C->vdim;
  }
  for (p=0; p<C->kappas; p++) free(cov->px[p]);
}

//static int ZZ = 0;
double Real(SEXP p, char *name, int idx) {
  char msg[200];
  // if(++ZZ==65724){printf("type=%d %d '%s'\n",ZZ,TYPEOF(p), CHAR(STRING_ELT(p,0)));cov_model *cov;crash(cov);}
  if (p != R_NilValue)
    switch (TYPEOF(p)) {
    case REALSXP :  return REAL(p)[idx];
    case INTSXP : return INTEGER(p)[idx]==NA_INTEGER  
	? NA_REAL : (double) INTEGER(p)[idx];
    case LGLSXP : return LOGICAL(p)[idx]==NA_LOGICAL ? NA_REAL 
	: (double) LOGICAL(p)[idx];
    }
  // MEMCOPY(msg, p, 300); print("%s\n", msg);
  sprintf(msg, "'%s' cannot be transformed to double! (type=%d)\n",
	  name, TYPEOF(p));  
  //printf("\n>>>> '%s'\n", CHAR(STRING_ELT(p, 0)));
  ERR(msg);
  return NA_REAL;  // to avoid warning from compiler
}



void Real(SEXP el,  char *name, double *vec, int maxn) {
  char msg[200];
   int i, j, n;
  if (el == R_NilValue) {
    sprintf(msg,"'%s' cannot be transformed to double.\n", name);
    ERR(msg);
  }
  n = length(el);
  for (j=i=0; i<maxn; i++) {
    vec[i] = Real(el, name, j);
    if (++j >= n) j=0;
  }
  return;
}

int Integer(SEXP p, char *name, int idx, bool nulltoNA) {
  char msg[200];
  if (p != R_NilValue) {
    switch(TYPEOF(p)) {
    case INTSXP : 
      return INTEGER(p)[idx]; 
    case REALSXP : 
      double value;
      value = REAL(p)[idx];      
      if (ISNAN(value)) {
	return NA_INTEGER;
	//sprintf(msg, "%s: NAs not allowed for integer valued parameters", name);
	//	ERR(msg);
      }
      if (value == trunc(value)) return (int) value; 
      else {
	sprintf(msg, "%s: integer value expected", name);
	ERR(msg);
      }
    case LGLSXP :
      return  LOGICAL(p)[idx]==NA_LOGICAL ? NA_INTEGER : (int) LOGICAL(p)[idx];
    }
  } else if (nulltoNA) return NA_INTEGER;
  sprintf(msg, "%s: unmatched type of parameter [type=%d]", name, TYPEOF(p));
  ERR(msg);
  return NA_INTEGER; // compiler warning vermeiden
}

int Integer(SEXP p, char *name, int idx) {
  return Integer(p, name, idx, false);
}


void Integer(SEXP el, char *name, int *vec, int maxn) {
  char msg[200];
  int i, j, n;
  if (el == R_NilValue) {
    sprintf(msg, "'%s' cannot be transformed to integer.\n",name);
    ERR(msg);
  }
  n = length(el);
  for (j=i=0; i<maxn; i++) {
    vec[i] = Integer(el, name, j);
    if (++j >= n) j=0;
  }
}



bool Logical(SEXP p, char *name, int idx) {
  char msg[200];
  if (p != R_NilValue)
    switch (TYPEOF(p)) {
    case REALSXP : return ISNA(REAL(p)[idx]) ? NA_LOGICAL : (bool) REAL(p)[idx];
    case INTSXP :
      return INTEGER(p)[idx]==NA_INTEGER ? NA_LOGICAL : (bool) INTEGER(p)[idx];
    case LGLSXP : return LOGICAL(p)[idx];
    }
  sprintf(msg, "'%s' cannot be transformed to logical.\n", name);  
  ERR(msg);
  return NA_LOGICAL;  // to avoid warning from compiler
}


char Char(SEXP el, char *name) {
  char msg[200];
  SEXPTYPE type;
  if (el == R_NilValue) goto ErrorHandling;
  type = TYPEOF(el);
  if (type == CHARSXP) return CHAR(el)[0];
  if (type == STRSXP) {
    if (length(el)==1) {
      if (strlen(CHAR(STRING_ELT(el,0))) == 1)
	return (CHAR(STRING_ELT(el,0)))[0];
      else if (strlen(CHAR(STRING_ELT(el,0))) == 0)
	return '\0';
    }
  }
 
 ErrorHandling:
  sprintf(msg, "'%s' cannot be transformed to character.\n",  name);  
  ERR(msg);
  return 0; // to avoid warning from compiler
}


#define INT Integer(el, name, 0)
#define LOG Logical(el, name, 0)
#define NUM Real(el, name, 0)
#define CHR Char(el, name)


double NonNegInteger(SEXP el, char *name) {
  int num;
  num = INT;
  if (num<0) {
    num=0; 
    char msg[200];
    sprintf(msg,"'%s' which has been negative is set 0.\n",name);
    warning(msg);
  }
  return num;
}

double NonNegReal(SEXP el, char *name) {
  double num;
  num = NUM;
  if (num<0.0) {
    num=0.0; 
    char msg[200];
    sprintf(msg,"%s which has been negative is set 0.\n",name);
    warning(msg);
   }
  return num;
}

double NonPosReal(SEXP el, char *name) {
  double num;
  num = NUM;
  if (num>0.0) {
    num=0.0; 
    char msg[200];
    sprintf(msg,"%s which has been positive is set 0.\n",name);
    warning(msg);
  }
  return num;
}

double PositiveInteger(SEXP el, char *name) {
  int num;
  num = INT;
  if (num<=0) {
    num=0; 
    char msg[200];
    sprintf(msg,"'%s' which has been negative is set 0.\n",name);
    warning(msg);
  }
  return num;
}

double PositiveReal(SEXP el, char *name) {
  double num;
  num = NUM;
  if (num<=0.0) {
    num=0.0; 
    char msg[200];
    sprintf(msg,"%s which has been negative is set 0.\n",name);
    warning(msg);
   }
  return num;
}

#define POS0INT NonNegInteger(el, name) /* better: non-negative */
#define POS0NUM NonNegReal(el, name)
#define NEG0NUM NonPosReal(el, name)
#define POSINT PositiveInteger(el, name) /* better: non-negative */
#define POSNUM PositiveReal(el, name)

void getUnits(SEXP el, char VARIABLE_IS_NOT_USED *name, 
	      char units[MAXUNITS][MAXUNITSCHAR], 
	      char units2[MAXUNITS][MAXUNITSCHAR]) {  
  int i, j, 
    l = length(el);
  if (TYPEOF(el) != NILSXP && TYPEOF(el) == STRSXP && l >= 1) {
    for (i=j=0; i<MAXUNITS; i++, j=(j + 1) % l) {
      strncpy(units[i], CHAR(STRING_ELT(el, j)), MAXUNITSCHAR);
      units[i][MAXUNITSCHAR - 1] ='\0';		    
      if (units2!=NULL) {
	strncpy(units2[i], CHAR(STRING_ELT(el, j)), MAXUNITSCHAR);
        units2[i][MAXUNITSCHAR - 1] ='\0';		    
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

int GetName(SEXP el, char *name, const char * List[], int n,
	    int defaultvalue) {
  char msg[1000], dummy[1000];
  int i,
    nM1 = n - 1;

  if (TYPEOF(el) == NILSXP) goto ErrorHandling;
 

  if (TYPEOF(el) == STRSXP) {
    int m = Match((char*) CHAR(STRING_ELT(el, 0)), List, n);

    if (m >= 0) return m; else {
      if (strcmp((char*) CHAR(STRING_ELT(el, 0)), " ") == 0 ||
	  strcmp((char*) CHAR(STRING_ELT(el, 0)), "") == 0) {
	goto ErrorHandling;
      }
    }
  }
  sprintf(dummy, "'%s': unknown value '%s'. Possible values are:", 
	  name, CHAR(STRING_ELT(el, 0)));
  for (i=0; i<nM1; i++) {
    sprintf(msg, "%s '%s',", dummy, List[i]);    
    strcpy(dummy, msg);
  }
  sprintf(msg,"%s '%s'.", dummy, List[i]);  
  ERR(msg);
 
 ErrorHandling:
  if (defaultvalue >= 0) return defaultvalue;
  
  sprintf(msg, "'%s': no value given.", name);
  ERR(msg);

  return 999;// to avoid warning from compiler
}

int GetName(SEXP el, char *name, const char * List[], int n) {
 return GetName(el, name, List, n, -1);
}




void CE_set(SEXP el, int j, char *name, ce_param *cp) {
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
  case 3: cp->maxmem = INT;   break;
  case 4: cp->tol_im = POS0NUM; break;
  case 5: cp->tol_re = NEG0NUM; break;
  case 6: cp->trials = NUM;
    if (cp->trials<1) {
      cp->trials=1;
      sprintf(msg, "%s is set to 1\n", name);
      warning(msg);
    }
    break;
  case 7: cp->useprimes = LOG; break;
  case 8: cp->dependent = LOG; break;
  case 9: cp->approx_grid_step = POS0NUM; break;
  case 10: cp->maxgridsize = POS0INT; break;
  default: ERR("unknown parameter for GLOBAL.general");
  }
}
 


#define prefixN 20
const char * prefixlist[prefixN] = 
  {"",  // -1
   "general", "gauss", "krige", 
   "circulant", "direct",  "nugget",//6, //"markov",
   "sequ", "spectral", "tbm",
   "mpp", "hyper", "maxstable",  // 12
   "br", "distr", "fit",         // 15
   "empvario", "gui", "graphics",// 18 
   "warn", // ACHTUNG warn wird nicht pauschal zurueckgesetzt
};
#define generalN 28
// IMPORTANT: all names of general must be at least 3 letters long !!!
const char *general[generalN] =
  { "modus_operandi", "printlevel", "storing", 
    "skipchecks", "every", "register", 
    "interpolregister", "condregister", "errregister", 
    "guiregister", "gridtolerance", "pch", 
    "practicalrange", "spConform", "cPrintlevel", 
    "exactness", "matrix_inversion", "matrix_tolerance",
    "allowdistanceZero", "na_rm_lines", "vdim_close_together",
    "expected_number_simu", "xyz_notation", "coordinate_system", 
    "new_coord_units", "coord_units", "variab_units", "seed"};

#define gaussN 5
const char *gauss[gaussN]= {"paired", "stationary_only", "approx_zero", 
			    "direct_bestvar", "loggauss"};

#define krigeN 7
const char *krige[krigeN] = {"method", "return_variance", "locmaxn",
			     "locsplitn", "locsplitfactor", "fillall", 
			     "cholesky_R"};
#define CEN 11
const char *CE[CEN] = {"force", "mmin", "strategy", 
		       "maxmem", "tolIm","tolRe", 
		       "trials", "useprimes", "dependent", 
		       "approx_step", "approx_maxgrid"};
#define directN 3
const char *direct[directN] = {"root_method", "svdtolerance", "max_variab"};
#define markovN 4
const char * markov[markovN] = {"neighbours", "precision", "cyclic", "maxmem"};
#define pnuggetN 1
const char * pnugget[pnuggetN] ={"tol"};
#define sequN 3
const char * sequ[sequN] ={"max_variables", "back_steps", "initial"};
#define spectralN 5
const char * spectral[spectralN] = {"sp_lines", "sp_grid", "ergodic", 
				    "prop_factor", "sigma"};
#define pTBMN 9
const char * pTBM[pTBMN] = {"tbmdim", "fulldim", "center", 
			    "points", "lines", "linessimufactor", 
			    "linesimustep", "layers", "grid"};
#define mppN 3
const char * mpp[mppN] = {"n_estim_E", // n to determine E by simulation
			  "intensity", 
			  // "refradius_factor", 
			  "about_zero"
			  // "plus", 
			  // "samplingdist", "samplingr",// MPP_cc
			  //"p", // Gneiting_cc
                         };
#define hyperN 4
const char * hyper[hyperN] = {"superpos", "maxlines", "mar_distr", "mar_param"};
#define extremeN 6
const char * extreme[extremeN] = 
  {"max_gauss", "maxpoints", "xi", "density_ratio", "check_every", "flat"};
#define brN 8
const char * br[brN] = 
  {"maxtrendmem", "meshsize", "lowerbound_optimarea", "vertnumber", 
   "optim_mixed", "optim_mixed_tol", "optim_mixed_maxpoints",
   "variobound"};

#define distrN 9
const char * distr[distrN] = 
  {"safety", "minsteplen", "maxsteps", "parts", "maxit", 
   "innermin", "outermax", "mcmc_n", "repetitions"};

#define fitN 39
const char * fit[fitN] = 
  {"bin_dist_factor", "upperbound_scale_factor", "lowerbound_scale_factor", 
   "lowerbound_scale_ls_factor","upperbound_var_factor","lowerbound_var_factor",
   "lowerbound_sill", "scale_max_relative_factor", "minbounddistance",
   "minboundreldist",  "approximate_functioncalls", "refine_onborder",
   "minmixedvar", "maxmixedvar", "solvesigma",
   "bc_lambda_lb", "bc_lambda_ub", "use_naturalscaling", 
   "bins", "nphi", "ntheta", 
   "ntime", "sill", "only_users", 
   "optim_var_estimation", "shortnamelength", "use_spam",
   "split", "scale_ratio",  "critical",
   "n_crit", "max_neighbours", "splitn_neighbours",
   "splitfactor_neighbours", "smalldataset", "min_diag",
   "reoptimise", "optimiser", "algorithm"
  };

#define empvarioN 5
const char * empvario[empvarioN] = 
  {"phi0", "theta0", "tol0",
   "pseudovariogram", "fft"};

#define guiN 3
const char * gui[guiN] = 
  {"alwaysSimulate", "simu_method", "size"};

#define graphicsN 4
const char *graphics[graphicsN]= {"always_close_screen" ,"grPrintlevel", "height", "increase_upto"};

#define warnN 9
const char * warns[warnN] =  { // Achtung ! warn parameter werden nicht
  // pauschal zurueckgesetzt
  "oldstyle", "newstyle", "newAniso", "ambiguous", "normal_mode", 
  "colour_palette", "warn_colour", "stored.init", "warn_coordinates"};

const char **all[] = {general, gauss, krige, CE, direct, // markov,
		      pnugget, sequ, spectral, pTBM, mpp,
		      hyper, extreme, br, distr,
		      fit, empvario, gui, 
		      graphics, warns};

int allN[] = {generalN, gaussN, krigeN, CEN, directN,// markovN,
	      pnuggetN,  sequN, spectralN, pTBMN, mppN,
	      hyperN, extremeN, brN, distrN, fitN, 
	      empvarioN, guiN, graphicsN, warnN};

void RelaxUnknownRFoption(int *relax) {
  RELAX_UNKNOWN_RFOPTION = (bool) *relax;
}


void setparameter(SEXP el, char *prefix, char *mainname, bool isList) {  
  int i,j,ii;
  char msg[200], name[200];
  
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
 	sprintf(msg, "option '%s' cannot be uniquely identified.", name);
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


  switch(i) {
  case 0: {// general
    general_param *gp;
    gp = &(GLOBAL.general);
    switch(j) {
    case 0: {
      int old_mode = gp->mode;
      SetDefaultModeValues(old_mode,
			   gp->mode = GetName(el, name, MODENAMES, nr_modes,
					      normal));
      break;
    }
    case 1: {
      int threshold = 1000; //PL_ERRORS;
      gp->Rprintlevel = INT;
      PL = gp->Cprintlevel = 
	gp->Rprintlevel <= threshold ? gp->Rprintlevel : threshold;
    }
      break;
    case 2: {
      bool storing = LOG;
      //  print("before setting storing %d %d\n", storing, KEY[0].simu.active);
      if  (length(el) > 1) {
	if (!storing) {	  
	  int nr = Integer(el, (char*) "storing (register)", 1);
	  if (nr != NA_INTEGER) {
	    //	    print("deleting register %d\n", nr);
	    if (nr<0 || nr>MODEL_MAX) 
	      ERR("RFoptions: number for register is out of range");
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
    case 5: { // simu
      int keynr = INT;
      if ((keynr<0) || (keynr>MODEL_MAX)) ERR("register number out of range");
      gp->keynr=keynr; }
      break;
    case 6: { // interpol
      int keynr;
      keynr = INT;
      if ((keynr<0) || (keynr>MODEL_MAX)) 
	ERR("interpolregister number out of range");
      gp->interpolregister=keynr;}
      break;
    case 7: { // cond
      int keynr;
      keynr = INT;
      if ((keynr<0) || (keynr>MODEL_MAX)) 
	ERR("condregister number out of range");
      gp->condregister=keynr;}
      break;
    case 8: { // err 
      int keynr;
      keynr = INT;
      if ((keynr<0) || (keynr>MODEL_MAX))
	ERR("errregister number out of range");
      gp->errregister=keynr;}
      break;
    case 9: { // gui
      int keynr;
      keynr = INT;
      if ((keynr<0) || (keynr>=MODEL_MAX)) 
	ERR("guiregister number out of range");
      gp->guiregister=keynr;}
      break;
    case 10: gp->gridtolerance = NUM; break;
    case 11: gp->pch = CHR;           break;
    case 12: 
      int n;
      n =INT;  
      if (n!=0 && n!=1 && n!=2 && n!=3)
	ERR("PracticalRange out of range. It should be TRUE or FALSE.");
      NS = gp->naturalscaling = n;
      break;
    case 13: gp->sp_conform = LOG;       break;
    case 14: PL = gp->Cprintlevel = INT;        break;
    case 15: gp->exactness = NUM;        break; 
    case 16: {
      int old[MAXINVERSIONS],
	l = length(el);
      if (l > MAXINVERSIONS) ERR("matrix_inversion: vector too long");
      for (ii=0; ii<l; ii++) {
	old[ii] = gp->matrix_inversion[ii];
	gp->matrix_inversion[ii] = Integer(el, name, ii); 
	if (gp->matrix_inversion[ii] < 0 || gp->matrix_inversion[ii] > 2) break;
      }
      if (ii<l) {
	for (; ii>=0; ii--) gp->matrix_inversion[ii] = old[ii];
	ERR("values of matrix_inversion out of range");
      } else {
	for (; ii<MAXINVERSIONS; ii++) gp->matrix_inversion[ii] = -1;
      }
      break;
    }
    case 17: gp->matrixtolerance = NUM; break;
    case 18: gp->allowdist0 = LOG; break;
    case 19: gp->na_rm_lines = LOG; break;
    case 20: gp->vdim_close_together = LOG;    
      if (gp->vdim_close_together) {
	gp->vdim_close_together = false;
 	ERR("'vdim_close_together' not programmed yet");
      }
      break;
    case 21: gp->expected_number_simu = POS0INT; break;
    case 22: gp->xyz_notation = INT; break;
    case 23: gp->coord_system =
	(coord_sys_enum) GetName(el, name, COORD_SYS_NAMES, nr_coord_sys,
				 coord_auto);
      break; 
    case 24: getUnits(el, name, gp->newunits, NULL);
      break;
    case 25: getUnits(el, name, gp->curunits, isList ? NULL : gp->newunits);
      break;
    case 26: getUnits(el, name, gp->varunits, NULL);      
    break;
    case 27: gp->seed = Integer(el, name, 0, true); break;
    break;
  default: ERR("unknown option  for 'general'");
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
    case 4: gp->loggauss = LOG; break;
    default: ERR("unknown option  for 'gauss'");
    }}
    break;
  case 2: { // krige
    krige_param *kp;
    char ret;
    kp = &(GLOBAL.krige);
    switch(j) {
    case 0: 
      ret = CHR;

      //      print("%d: %d %d %d %d; %d %d %d >%s<\n", 
      //	     ret , 'O', 'S', 'I', 'A',
      //	     TYPEOF(el), CHARSXP, STRSXP, CHAR(STRING_ELT(el,0))
      //     );

      if (ret != 'A' && // auto
	  ret != 'S' && // simple
	  ret != 'O' && // ordinary	  
	  ret != 'M' && // mean	  
	  ret != 'U' && // universal
	  ret != 'I'    // intrinsic
	  ) {  
	ERR("krige.method not in S)imple O)rdinary U)niversal I)ntrinsic and A)utomatic");
      }
      kp->method = ret;
      break;
    case 1: kp->ret_variance = LOG; break;
    case 2: kp->locmaxn = POS0INT; break;
    case 3: {
      int  old[3];
      if (length(el) < 3) ERR("krige.locsplitn must have 3 components");
      for (ii=0; ii<3; ii++) {
	old[ii] = kp->locsplitn[ii];
	kp->locsplitn[ii] = Integer(el, name, ii); 
      }
      if (kp->locsplitn[0] < kp->locsplitn[2] || 
	  kp->locsplitn[2] < kp->locsplitn[1]) {
	for (ii=0; ii<3; ii++) kp->locsplitn[i] = old[ii];
	error("locsplitn[1] >= locsplitn[3] >= locsplitn[2] not satisfied");
      }
      break;
    }
    case 4: kp->locsplitfactor = POS0INT; break;
    case 5: kp->fillall = LOG; break;
    case 6: kp->cholesky = LOG; break;
    default: ERR("unknown option for 'krige'");
    }}
    break;
  case 3: // CE
    CE_set(el, j, name, &(GLOBAL.ce)); 
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
   default: ERR("unknown option for 'direct'");
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
    default: ERR("unknown option for 'nugget'");
    }}
    break;
  case 6: {//sequ,
    sequ_param *sp;
    sp = &(GLOBAL.sequ) ;
    switch(j) {
    case 0: sp->max = POS0INT;   break;
    case 1: sp->back = INT; if (sp->back < 1) sp->back=1;   break;
    case 2: sp->initial = INT; break;
    default: ERR("unknown option for 'sequ'");
    }}
    break;
  case 7: { // spectral,
    spectral_param *sp;
    sp = &(GLOBAL.spectral) ;
    switch(j) {
    case 0: Integer(el, name, sp->lines, MAXTBMSPDIM); break;
    case 1: sp->grid = LOG; break;
    case 2: sp->ergodic = LOG; break;
    case 3: sp->prop_factor = POS0NUM;
      if (sp->prop_factor <= 0.1) {
	sp->prop_factor=0.1;
	warning("'spectral.prop.factor' less than 0.1. Set to 0.1.");
      }
      break;
    case 4: sp->sigma = NUM; break;
    default: ERR("unknown option for 'spectral'");
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
    default: ERR("unknown option for 'TBM'");
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
    default: ERR("unknown option for 'mpp'");      
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
    default: ERR("unknown option for 'hyper'");
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
    case 5: ep->flat = INT; 
      if (ep->flat < -1 || ep->flat > 1) ERR("illegal value for 'flat'");
      break;
    default: ERR("unknown option for 'maxstable'"); 
    }}
    break;
  case 12 : { // br
    br_param *ep;
    ep = &(GLOBAL.br);
    switch(j) {
    case 0: ep->BRmaxmem = POSINT; break;
    case 1: ep->BRmeshsize = POSNUM; break;
    case 2: ep->BR_LB_optim_area = NUM; break;
    case 3: ep->BRvertnumber = POSINT; break;
    case 4: ep->BRoptim = POS0INT; break;
    case 5: ep->BRoptimtol = POS0NUM; break;
    case 6: ep->BRoptimmaxpoints = POS0INT; break;
    case 7: ep->variobound = NUM; break;
    default: ERR("unknown option for 'maxstable'"); 
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
    default: ERR("unknown option for 'maxstable'"); 
    }}
    break;
  case 14: { // fit
    fit_param *ef;
    ef = &(GLOBAL.fit);
    switch(j) {
    case 0: ef->bin_dist_factor = POS0NUM; break;
    case 1: ef->upperbound_scale_factor = POS0NUM; break; 
    case 2: ef->lowerbound_scale_factor = POS0NUM; break; 
    case 3: ef->lowerbound_scale_LS_factor = POS0NUM; break; 
    case 4: ef->upperbound_var_factor = POS0NUM; break; 
    case 5: ef->lowerbound_var_factor = POS0NUM; break;
    case 6: ef->lowerbound_sill = POS0NUM; break; 
    case 7: ef->scale_max_relative_factor = POS0NUM; break; 
    case 8: ef->minbounddistance = POS0NUM; break;
    case 9: ef->minboundreldist = POS0NUM; break;  
    case 10: ef->approximate_functioncalls = POS0INT; break; 
    case 11: ef->refine_onborder = LOG; break;
    case 12: ef->minmixedvar = POS0NUM; break; 
    case 13: ef->maxmixedvar = POS0NUM; break; 
    case 14: ef->solvesigma = NUM; break;
    case 15: ef->BC_lambdaLB = NUM; break;
    case 16: ef->BC_lambdaUB = NUM; break; 
    case 17: ef->use_naturalscaling = LOG; break;
    case 18: ef->bins = POS0INT; break;
    case 19: ef->nphi = POS0INT; break;
    case 20: ef->ntheta = POS0INT; break; 
    case 21: ef->ntime = POS0INT; break;
    case 22: ef->sill = POS0NUM; break;
    case 23: ef->onlyuser = LOG; break;
    case 24: ef->optim_var_estim =
	GetName(el, name, OPTIM_VAR_NAMES, nOptimVar); ; break; 
    case 25: { // mle
      ii=POS0INT; 
      if (ii==0) { ii = 1; warning("shortnamelength set to 1"); }
      if (ii>255) { ii = 255; warning("shortnamelength set to 255"); }
      ef->lengthshortname=ii;
    }
    break;
    case 26: ef->usespam = NUM; break;
    case 27: ef->split = LOG; break;
    case 28: ef->scale_ratio = NUM; break;
    case 29: ef->critical = INT; break;
    case 30: ef->n_crit = POS0INT; break;
    case 31: ef->locmaxn = POS0INT; break;
    case 32: {
      if (length(el) < 3) ERR("fit.locsplitn must have 3 components");
      for (ii=0; ii<3; ii++)
	ef->locsplitn[ii] = Integer(el, name, ii); 
      break;
    }
    case 33: ef->locsplitfactor = POS0INT; break;
    case 34: ef->smalldataset = POS0INT; break;
    case 35: ef->min_diag = NUM; break;
    case 36: ef->reoptimise = LOG; break;
    case 37: ef->optimiser = GetName(el, name, OPTIMISER_NAMES, nOptimiser, 0); 
      break;
    case 38: 
      switch(ef->optimiser) {
      case 3 : ef->algorithm = GetName(el, name, NLOPTR_NAMES, nNLOPTR); 
	break;
      default: ef->algorithm = -1;
      }
      break;
    default: ERR("unknown option for 'fit'");
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
   default: ERR("unknown option for 'empvario'");
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
    case 2: {
      int sizedummy[2];
      if (length(el) != 2) ERR("length of 'size' must be 2");
      Integer(el, name, sizedummy, 2);
      for (ii=0; ii<2; ii++) {	
	if (sizedummy[ii] <= 1) 
	  ERR("grid size in RFgui must contain at least 2 points");
      }
      for (ii=0; ii<2; ii++) {	 gp->size[ii] = sizedummy[ii]; }
    }
      break;
    default: ERR("unknown option  for 'gui'");
    }}
    break;
 case 17: { // graphics
    graphics_param *gp = &(GLOBAL.graphics);
    switch(j) {
    case 0 : gp->always_close = LOG; break;
    case 1 : gp->PL = INT; break;
    case 2 : gp->height = NUM; break;
    case 3 : {
      int uptodummy[2];
      if (length(el) != 2) ERR("length of 'increase_upto' must be 2");
      Integer(el, name, uptodummy, 2);
      for (ii=0; ii<2; ii++) {
	if (uptodummy[ii] <= 0)  ERR("increase_upto must be positive");
      } 
      for (ii=0; ii<2; ii++) { gp->increase_upto[ii] = uptodummy[ii]; }
    }
    break;
    default: ERR("unknown option  for 'graphics'");
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
  case 18: if (!isList) {
    warn_param *wp = &(GLOBAL.warn);
    switch(j) {
    case 0: wp->oldstyle = LOG;       break;
    case 1: wp->newstyle = LOG;       break;
    case 2: wp->Aniso = LOG;       break;
    case 3: wp->ambiguous = LOG;       break;
    case 4: wp->normal_mode = LOG;       break;
    case 5: wp->warn_mode = LOG;       break;
    case 6: wp->warn_colour = LOG;       break;
    case 7: wp->stored_init = LOG;       break;
    case 8: wp->warn_coordinates = LOG;       break;
    default: ERR("unknown option for 'general'");
    }}
    break;

  default: ERR("unknown option.");  
  }
  
/*
 case 6: {// markov, 
    markov_param *mp;
    mp = &(GLOBAL.markov) ;
    switch(j) {
    case 0: mp->neighbours= INT;
      if (mp->neighbours < 2) {
	if (PL>=PL_ IMPORTANT) { 
	  warning("minimal neighbourhood for Markov is 2"); }
	mp->neighbours = 2;
      } else if (mp->neighbours > 3) {
	if (PL>=PL_ IMPORTANT) {
	  warning("maximal neighbourhood for Markov is 3"); }
	mp->neighbours = 3;
	  }
      break;
    case 1: mp->precision = INT;   break;
    case 2: mp->cyclic = INT; break;
    case 3: mp->maxmem = POS0INT; break;
    default: ERR("unknown option for 'markov'");
    }}
    break;
*/

}

SEXP ExtendedInteger(double x) {
  return ScalarInteger(R_FINITE(x) ? x : NA_INTEGER);
}

SEXP ExtendedBoolean(double x) {
  return ScalarLogical(ISNA(x) || ISNAN(x) ? NA_LOGICAL : x != 0.0);
}

SEXP getRFoptions() {
  SEXP list, names, sublist[prefixN-1], subnames[prefixN-1];
  
  int i, k = 0;
  char x[2]=" ";
  int trueprefixN = prefixN - 1;
  //  cov_fct *C = CovList + cov->nr;

  PROTECT(list = allocVector(VECSXP, trueprefixN));
  PROTECT(names = allocVector(STRSXP, trueprefixN));
  
  for (i=0; i<trueprefixN; i++) {    
    SET_STRING_ELT(names, i, mkChar(prefixlist[i+1]));
    PROTECT(sublist[i] = allocVector(VECSXP, allN[i]));
    PROTECT(subnames[i] = allocVector(STRSXP, allN[i]));
    int endfor = allN[i];
    for (k=0; k<endfor; k++) {
      // print("%d %d %s$%s\n", i, k, prefixlist[i+1], all[i][k]);
      SET_STRING_ELT(subnames[i], k, mkChar(all[i][k])); 
    }
  }

#define ADD(ELT) SET_VECTOR_ELT(sublist[i], k++, ELT)
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
    ADD(ScalarInteger(p->keynr));    
    ADD(ScalarInteger(p->interpolregister));    
    ADD(ScalarInteger(p->condregister));    
    ADD(ScalarInteger(p->errregister));    
    ADD(ScalarInteger(p->guiregister));    
    ADD(ScalarReal(p->gridtolerance));
    ADDCHAR(p->pch);
    if (p->naturalscaling==0 || p->naturalscaling==1) 
      ADD(ScalarLogical(p->naturalscaling));
    else
      ADD(ScalarInteger(p->naturalscaling));
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
    ADD(ExtendedInteger(p->xyz_notation));    
    ADD(ScalarString(mkChar(COORD_SYS_NAMES[p->coord_system])));
    ADD(UNITS(p->newunits));
    ADD(UNITS(p->curunits));
    ADD(UNITS(p->varunits));
    ADD(ScalarInteger(p->seed));
   
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
  }

  i++; {
    k = 0;
    krige_param *p = &(GLOBAL.krige);
    ADDCHAR(p->method);
    ADD(ScalarLogical(p->ret_variance));   
    ADD(ScalarInteger(p->locmaxn));
    SET_VECTOR_ELT(sublist[i], k++, Int(p->locsplitn, 3, 3));
    ADD(ScalarInteger(p->locsplitfactor));
    ADD(ScalarLogical(p->fillall));   
    ADD(ScalarLogical(p->cholesky));   
 }

  i++; {
    k = 0;
    ce_param *p = &(GLOBAL.ce);
    ADD(ScalarLogical(p->force)); 
    SET_VECTOR_ELT(sublist[i], k++, Num(p->mmin, MAXCEDIM, MAXCEDIM)); 
    ADD(ScalarInteger(p->strategy)); 
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

  /*
   i++; {
    k = 0;
    markov_param *p = &(GLOBAL.markov);
    ADD(ScalarInteger(p->neighbours));    
    ADD(ScalarReal(p->precision));   
    ADD(ScalarInteger(p->cyclic));    
  }
  */

 
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
    ADD(ScalarLogical(p->ergodic));
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
  }

  i++; {
    k = 0;
    br_param *p = &(GLOBAL.br);
    ADD(ScalarInteger(p->BRmaxmem)); 
    ADD(ScalarReal(p->BRmeshsize));
    ADD(ScalarReal(p->BR_LB_optim_area));
    ADD(ScalarInteger(p->BRvertnumber));
    ADD(ScalarInteger(p->BRoptim));
    ADD(ScalarReal(p->BRoptimtol));
    ADD(ScalarInteger(p->BRoptimmaxpoints));
    ADD(ScalarReal(p->variobound));
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
    ADD(ExtendedBoolean(p->lowerbound_sill)); 
    ADD(ScalarReal(p->scale_max_relative_factor)); 
    ADD(ScalarReal(p->minbounddistance));
    ADD(ScalarReal(p->minboundreldist));  
    ADD(ScalarInteger(p->approximate_functioncalls)); 
    ADD(ScalarLogical(p->refine_onborder));
    ADD(ScalarReal(p->minmixedvar)); 
    ADD(ScalarReal(p->maxmixedvar)); 
    ADD(ScalarReal(p->solvesigma));
    ADD(ScalarReal(p->BC_lambdaLB));
    ADD(ScalarReal(p->BC_lambdaUB)); 
    ADD(ScalarLogical(p->use_naturalscaling));
    ADD(ScalarInteger(p->bins));
    ADD(ScalarInteger(p->nphi));
    ADD(ScalarInteger(p->ntheta)); 
    ADD(ScalarInteger(p->ntime));
    ADD(ScalarReal(p->sill));
    ADD(ScalarLogical(p->onlyuser));
    ADD(ScalarString(mkChar(OPTIM_VAR_NAMES[p->optim_var_estim])));
    ADD(ScalarInteger(p->lengthshortname));
    ADD(ExtendedBoolean(p->usespam));
    ADD(ScalarLogical(p->split));
    ADD(ScalarReal(p->scale_ratio));
    ADD(ScalarInteger(p->critical));
    ADD(ScalarInteger(p->n_crit));
    ADD(ScalarInteger(p->locmaxn));
    SET_VECTOR_ELT(sublist[i], k++, Int(p->locsplitn, 3, 3));
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
 }


  i++; {
    k = 0;
    warn_param *p = &(GLOBAL.warn);
    ADD(ScalarLogical(p->oldstyle));
    ADD(ScalarLogical(p->newstyle));
    ADD(ScalarLogical(p->Aniso));
    ADD(ScalarLogical(p->ambiguous));
    ADD(ScalarLogical(p->normal_mode));
    ADD(ScalarLogical(p->warn_mode));
    ADD(ScalarLogical(p->warn_colour));
    ADD(ScalarLogical(p->stored_init));
    ADD(ScalarLogical(p->warn_coordinates));
  }


  //   print("%d %d\n", i, prefixN -1);
   assert (i == trueprefixN - 1); // general has two options for prefixes

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
    strncpy(mainname, name, 200);
  } else {
    strncpy(prefix, name, i);
    prefix[i] = '\0';
    strcpy(mainname, name+i+1);
  }

  //  printf("i=%d %d\n", i, len);
  setparameter(el, prefix, mainname, isList);
  // printf("ende\n");
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
      ERR("'LIST' needs as argument the output of RFoptions");
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
      el = CAR(options);
      name = (char*) (isNull(TAG(options)) ? "" :CHAR(PRINTNAME(TAG(options))));
      splitAndSet(el, name, isList);
    }
    //       print("end2 %f\n", GLOBAL.gauss.exactness);
  }
  return(R_NilValue);
}


void GetModelRegister(char **name, int* nr) {
  *nr = Match(*name, REGNAMES, MODEL_MAX+1);

  //  print("%d\n", *nr);

  if (*nr<0 || *nr > MODEL_MAX) error("name for model register unknown");
}
  

void MultiDimRange(int *model_nr, double *natscale) { 
  MultiDimRange(KEY[*model_nr], natscale); 
}

void countelements(int *idx, int *N, int *boxes) {
  int i,
    n = *N; 
  for (i=0; i<n; i++) {
    //  print("%d %d  cum=%d %d %d\n", 
    //	   i, idx[i], cumgridlen[0], cumgridlen[1], cumgridlen[2]);
    boxes[idx[i]]++;
  }
}


void countneighbours(int *Xdim, int *parts, int *Squarelength, int *cumgridlen,
		     int *boxes, int * neighbours, int *OK) {
  int d, sum, totcumlen, relstart, x, y,
    nb[MAXGETNATSCALE], loc[MAXGETNATSCALE], e,
    sl = *Squarelength,
    dim = *Xdim,
    boundary = (sl - 1) / 2,
    //   total = cumgridlen[dim],
    maxn = GLOBAL.krige.locmaxn;

  assert(dim <= MAXGETNATSCALE);
 
  sum = 0;
  *OK = true;
 
  x = 0;
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
      *OK = false;
      return;
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
    for (i=0; i<total; i++) {
      if (elm[i] != NULL) free(elm[i]);
    }
    free(elm);
  }
  if (count != NULL) free(count);
 
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
    for (i=0; i<total; i++) {
      if (neighb[i] != NULL) free(neighb[i]);
    }
    free(neighb);
  }
  if (err!=NOERROR) XERR(err);

  return subnei;
}


void DeleteKey(int *reg) {
  COV_DELETE(KEY + *reg);
}



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
