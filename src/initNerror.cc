
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part and error messages

 Copyright (C) 2001 -- 2015 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/



// to do: V3.1+: Randverteilungen der RPs in cov, D etc implementieren;
//        do -> random part; init wird innerhalb von do aufgerufen,
//        falls nicht initialisiert oder random submodels??
// to do: MLE: random parameters einsammeln

#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
//#include <sys/timeb.h>
 
#include <string.h>
#include "RF.h"
#include "Operator.h"
#include <unistd.h>
#include "SpherModelsInitNerror.h"

int gaussmethod[Forbidden+1];

cov_model *KEY[MODEL_MAX+1];
double  ZERO[MAXSIMUDIM], 
  ONE = 1,
//    *userdefinedCovMatrix[MAXDEFMATRIX][MAXMAKEEXPLICITE],
    *OutX=NULL, 
    *OutY=NULL;
int NATSC_INTERN,NATSC_USER,
  GENERALISEDCAUCHY, STABLE,  BROWNIAN, CAUCHY, 
  GAUSS, NUGGET, PLUS, TBM2NR, BALL, ECF, MULT, GENNSST_INTERN, PROD,
  DISTRIBUTION, DETERM_DISTR, GAUSS_DISTR, SETPARAM, COVFCTN,
  COVMATRIX, RFGET, STROKORB_MONO, STROKORB_BALL_INNER, RECTANGULAR, CONST, BIND,
  POLYGON, SCATTER, MCMC_PGS, MCMC,LIKELIHOOD_CALL, LINEARPART_CALL,
  PREDICT_CALL,
  IDCOORD, MULT_INVERSE, GNEITING_INTERN,
  TRUNCSUPPORT, SHAPESTP, SHAPEAVE, BROWNRESNICK, UNIF, MPPPLUS, CUTOFF, STEIN,
  BRSHIFTED_USER, BRMIXED_USER, BRORIGINAL_USER, 
  BRSHIFTED_INTERN, BRMIXED_INTERN, BRORIGINAL_INTERN,   
  EXTREMALGAUSSIAN, RANDOMSIGN,  
  ARCSQRT_DISTR, SHAPEPOW, POW,
  PTS_GIVEN_SHAPE, STATIONARY_SHAPE, STANDARD_SHAPE,
  LOC, SET_DISTR, SCALESPHERICAL, TBM_OP, USER,TREND, // TREND_PROC,
  COVARIATE, TRAFO,
  MIXEDEFFECT, // MLEMIXEDEFFECT,
//DENSITY, 
  VARIOGRAM_CALL, SIMULATE,
  MISSING_COV, NULL_MODEL,
  POWER_DOLLAR, DOLLAR_PROC, LASTDOLLAR, DOLLAR, PLUS_PROC,
  MPPPLUS_PROC, MULT_PROC, 
  BINARYPROC, BROWNRESNICKPROC,
  GAUSSPROC, POISSONPROC,  SCHLATHERPROC, SMITHPROC, CHI2PROC,
  EXTREMALTPROC, TPROC, TRENDEVAL,
  NUGGET_USER, NUGGET_INTERN,
  CIRCEMBED,  SPECTRAL_PROC_USER, SPECTRAL_PROC_INTERN,
  DIRECT, SEQUENTIAL, SPECIFIC, SELECTNR,
  AVERAGE_USER, AVERAGE_INTERN, HYPERPLANE_USER, HYPERPLANE_INTERN,
  RANDOMCOIN_USER, CE_CUTOFFPROC_USER, CE_CUTOFFPROC_INTERN, 
  CE_INTRINPROC_USER, CE_INTRINPROC_INTERN, TBM_PROC_USER, TBM_PROC_INTERN, 
  VECTOR,
  ISO2ISO, SP2SP, SP2ISO, S2ISO, S2SP, S2S, SId, E2E, E2EIso, 
  E2Sph, E2SphIso, Sph2Sph, Sph2SphIso, LASTGATTER,
  FIRST_PLANE, LAST_PLANE, EARTHKM2CART, EARTHMILES2CART,
  EARTHKM2GNOMONIC, EARTHMILES2GNOMONIC,
  EARTHKM2ORTHOGRAPHIC, EARTHMILES2ORTHOGRAPHIC,  
  FIRST_TRAFO, LAST_TRAFO;
// userdefinedCM_RC[MAXDEFMATRIX][MAXMAKEEXPLICITE], 
bool
    NAOK_RANGE=false;
char CovNames[MAXNRCOVFCTS][MAXCHAR], CovNickNames[MAXNRCOVFCTS][MAXCHAR];
char NEWMSG[LENERRMSG];

cov_fct *CovList=NULL;
int currentNrCov=-1;
int currentRegister=-1;
 
int True=1; // never change
int False=0; // never change
char *FREEVARIABLE= (char*) "...";


//globalorig GLOBALORIG = {false, {}};
int PrInL=-1;				


pref_type PREF_ALL = {PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, 
		      PREF_BEST, PREF_BEST, PREF_NONE, PREF_BEST, PREF_BEST,
		      PREF_BEST, PREF_BEST, PREF_NONE, // specific
		                                       PREF_BEST, PREF_BEST},
  PREF_NOTHING = {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE,
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_BEST, // nothing
		                                              PREF_NONE},
  PREF_TREND =  {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
		  PREF_NONE, PREF_NONE, PREF_BEST, PREF_NONE, PREF_NONE,
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_BEST, // nothing
		                                              PREF_NONE},
  PREF_AUX = {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
	      PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE,
	      PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE};
 



/*
  record containing all user specified directions for choosing
  a simulation method. Currently only one element is contained.
*/

double GENERAL_PRECISION = 5e-15;
double EIGENVALUE_EPS = 1e-15;

int ERRORMODELNUMBER = -1;

char ERRMSG[LENERRMSG], MSG[LENERRMSG], BUG_MSG[LENMSG], MSG2[LENERRMSG],
  ERRORSTRING[MAXERRORSTRING], ERROR_LOC[nErrorLoc];



char PREF_FAILURE[100 * Nothing];


const char *METHODNAMES[Forbidden+1]={"circulant", //0
				      "cutoff",
				      "intrinsic",
				      "tbm", 
				      "spectral", //4
				      "direct",
				      "sequential",
				      "trend",
				      "average",
				      "nugget", //9
				      "coins",
				      "hyperplane",
				      "specific",
				      "any method", // nothing
				      "forbidden"},
  *ROLENAMES[ROLE_LAST + 1] = {
    "<none>",                                                     // 0
    "covariance model", "Gauss", "max-stable", "BrownResnick", "Smith",  // 5
    "Schlather", "Poisson", "PoissonGauss", "Hyperplane", "Bernoulli", // 10
    "distribution", "likelihood",
    "<rotten>", "<undefined>"},
  *CAT_TYPENAMES[OtherType + 1] = {
    // TcfType, PosDefType, VariogramType, NegDefType, ProcessType,
    // GaussMethodType, BrMethodType, PointShapeType, RandomType, ShapeType, 
    // TrendType, InterfaceType, RandomOrShape, UnDefinedType, MathDefinition, 
    // OtherType
    "RM", "RM", "RM", "RM", "RP", 
    "RP", "RP", "RM", "RR", "RM", 
    "RM", "RF", "RM", "RM", "R.", 
    "RO"},
  *REGNAMES[MODEL_MAX+1] = {"reg0", "reg1", "reg2", "reg3", "reg4", 
			    "reg5", "reg6", "reg7", "reg8", "reg9",
			    "user", "unused", "intern", "split", "gui",
			    "mle", "mlesplit", "mletrend", "mlebounds",
		 	    "kriging", "conditional", "error model"},
  *POSITIVITY_NAMES[(int) pt_mismatch + 1] = 
			      {"pt-wise positive definite", 
			       "pt-wise indefinite", 
			       "pt-wise negative definite", "pt-wise zero", 
			       "pt-wise param dependent", 
			       "pt-wise submodel dependent", 
			       "pt-wise undefined", "pt-wise unknown", 
			       "pt-wise mismatch"};

//  int SYS_TO_ISO[nr_coord_sys_proj] =  {GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ};

char
  STANDARDPARAM[MAXPARAM][MAXCHAR],
  STANDARDSUB[MAXSUB][MAXCHAR];

bool RELAX_UNKNOWN_RFOPTION=false; // auf keinen Fall aendern!

void errorMSG(int err, char* M, int len, bool final) {
  char m[LENERRMSG], m2[LENERRMSG];
  if (err >= ERRORM && err <= ERRORMEND) err = ERRORM;

  switch (err) {
  case NOERROR : strcpy(m,"none"); break;
  case NOERROR_REPEAT : strcpy(m,"none; looking for further covariances applicable to the same method");break;
  case NOERROR_ENDOFLIST : strcpy(m,"none; end of list");break;
  case ERRORDUMMY : strcpy(m,"none (dummy)"); break;
  case ERRORNOTDEFINED :       
    strcpy(m,"specified method undefined for the given model or no simulation method found for the given model");break;
  case ERRORCARTESIAN :
    strcpy(m,"only cartesian system allowed (currently)"); break;
  case ERRORNOTPROGRAMMEDYET :    
    strcpy(m,"Not programmed yet in RandomFields Version 3. Sorry."); break;
  case ERRORVDIMNOTPROGRAMMEDYET :    
    strcpy(m,"multivariate version not programmed yet. Sorry."); break;
  case ERRORTYPECONSISTENCY :
     strcpy(m,"incorrect choice of submodel (type inconsistency)"); break;
  case ERRORFAILED: 
   strcpy(m,"algorithm failed (partially)");break;
  case ERRORMEMORYALLOCATION: 
    strcpy(m, "memory allocation error -- too much space demanded or non-positive number of bytes requested"); 
    break;
  case ERRORNOTINITIALIZED: 
    strcpy(m,"not initialized or storing=FALSE");break;
  case ERRORDECOMPOSITION:
    strcpy(m,"covariance function does not seem to be (strictly) positive definite");break;
  case ERRORNOMULTIVARIATE :
    strcpy(m, "multivariate models not allowed (yet)"); 
    break;
  case ERROR_MATRIX_SQUARE :
    strcpy(m, "square matrix expected"); break;
  case ERROR_MATRIX_VDIM :
    strcpy(m, "size of matrix is not a multiple of the multivariate dimension"); break;
  case ERROR_MATRIX_POSDEF :
    strcpy(m, "matrix does not seem to be strictly positive definite"); break;
    //  case ERROR_MATRIX_ :   strcpy(m, ""); break;
  case ERRORDIM: 
    //
    //    { printf("error dimension\n"); cov_model *cov; crash(cov); }
    sprintf(m,"dimension specification not in [1,%d] or dimension of coordinates larger than that the supposed spatio-temporal process",
	    MAXSIMUDIM);break;
  case ERRORWAVING :
    strcpy(m,"Rescaling not possible (waving or large nugget effect?)");break;
  case ERRORRESCALING:
    strcpy(m,"practical range not defined");
    break;
  case ERRORNOSTATMATCH : 
    strcpy(m,"no matching assumption found for the domains");
    break;
  case ERRORANISO:
    strcpy(m,"anisotropic call not allowed"); break; 
  case ERRORUNKNOWNMETHOD:
    strcpy(m,"Unknown method in this context or unallowed mixture of methods"); 
    break;
  case ERRORWRONGDIM:
    strcpy(m,"wrong dimension"); break;
   case ERRORUNKOWNSXPTYPE:
    strcpy(m, "parameter value of unknown SXP type");
    break;
  case ERROROUTOFMETHODLIST:
    char restrictive[100], info[150];
    sprintf(restrictive, "Are the %s() too restrictive?", RFOPTIONS);
    sprintf(info, "\n You get (more) internal information if you set %s(%s=%d) before running your code.",
	    RFOPTIONS, general[GENERAL_CPRINT], PL_DETAILSUSER);
    
    sprintf(m, 
	    "Running out of list of methods. %s%s",
	    GLOBAL.general.skipchecks
	    ? "Did you try an invalid parameter combination?"
	    : restrictive,
	    PL <= 2 ? info : "" );
    break;
  case ERRORREGISTER: 
    strcpy(m, "register number out of range");
    break;
  case ERRORWRONGISO: 
    strcpy(m, "impossible sequence of isotropies or coordinate systems");
    break;
  case ERRORKERNEL:
    strcpy(m, "the mapping 'earth -> cartesian' keeps definiteness only if it is used as a kernel.");
    break;
  case ERRORM: 
    strcpy(m, ERRORSTRING);
    break;
  case ERRORWRONGVDIM: 
    strcpy(m, ERRORSTRING);
    break;
  case ERRORBADVDIM: 
    strcpy(m, "m-dimensionality could not be detected");
    break;

  case ERRORNOTCARTESIAN:
       strcpy(m, "Currently only cartesian coordinate systems are allowed");
    break;
  case ERRORODDCOORDTRAFO:
    strcpy(m, "coordinate transformation not possible or not programmed yet");
    break;
    // extremes:
  case ERRORSUBMETHODFAILED:
    sprintf(m, "no good submethods exist");
  case  ERRORSTATVARIO:
    strcpy(m, 
	   "negative definite function expected depending on 1 variable only");
    break;
   case ERRORNOVARIOGRAM:
    strcpy(m, "Variogram model not allowed in this context");
    break;
  case ERRORNORMALMIXTURE:
    strcpy(m, "only normal mixtures as first submodel allowed (Gneiting, 2002)");
    break;
  case ERRORMAXDIMMETH:
    strcpy(m, "maximal dimension of variables for the method exceeded");
    break;
  case ERRORPREVDOLLAR:
    strcpy(m, "method may not be initialised by preceding initS");
    break;
  case ERRORSPECTRAL: 
    strcpy(m, "submodel does not have spectral representation");
    break;    
  case ERRORTBMCOMBI: 
    strcpy(m, "the given combination of 'fulldim' and 'reduceddim' is not possible yet.");
    break;    

  case ERRORINVALIDMODEL : // gauss distribution, no method
    strcpy(m, "Invalid covariance model: did you wrongly use an auxiliary function to construct the model?");
    break;    
  case ERRORODDMODEL : // gauss distribution, no method
    strcpy(m, "Odd covariance model: the use of auxiliary functions and/or your choice of the parameters lead to a covariance model for which no simulation methods exist.");
    break;    
  case ERRORANISO_T :
    sprintf(m, "'%s' may not be given at the same time with '%s' or '%s'", 
	    CovList[DOLLAR].kappanames[DANISO], 
	    CovList[DOLLAR].kappanames[DAUSER], 
	    CovList[DOLLAR].kappanames[DPROJ]);
    break;
  case ERRORDIAMETERNOTGIVEN:
    strcpy(m, "Diameter must always be given");
    break;
  case ERRORPREFNONE:
    strcpy(m, "the simulaton method does not allow for the given model.");
    break;
  case ERRORPREFNONECOV:
    strcpy(m, "the given model does not allow for calculation the covariance values.");
    break;
   
    //    case : strcpy(m,"");break;
    //
    // Poisson:
  case ERRORUNKNOWNMAXTYPE :
    strcpy(m, "unknown type of max-stable process");
    break;
 
  case ERRORATOMP :
    strcpy(m, "p must be given everywhere or nowhere");
    break;
   
  case ERRORKRIGETOL :
    strcpy(m,"sigma must be at most KRIGE_TOLERANCE");
    break;


  case MSGLOCAL_OK :
    strcpy(m,"fine");
    break;
  case MSGLOCAL_JUSTTRY :
    strcpy(m,
	   "unclear whether algorithm will work for specified parameters");
    break;
  case MSGLOCAL_NUMOK :
    strcpy(m,"fine. Algorithm should work for specified parameters");
    break;
  case MSGLOCAL_ENDOFLIST :
    strcpy(m,"end of list for variants of the algorithm");
    break;
  case MSGLOCAL_SIGNPHI :
    strcpy(m,"wrong sign of covariance function at break point");
    break;
  case MSGLOCAL_SIGNPHIFST :
    strcpy(m, "wrong sign of 1st derivative of the covariance function at the break point");
    break;
  case MSGLOCAL_SIGNPHISND :
    strcpy(m, "wrong sign of 2nd derivative of the covariance function at the break point");
    break;
  case MSGLOCAL_INITINTRINSIC :
    strcpy(m,"one of a2, b or a0+phi(0) has wrong sign");
    break;
  case ERRORUNSPECIFIED :
    strcpy(m,"(unspecified)");
    break;
   default : 
     PRINTF(" error=%d\n", err); 
     // crash();
     BUG;
  }

  if (final) sprintf(m2, "%s %s", ERROR_LOC, m);
  else strcpy(m2, m);
  
  if (strlen(m) > (unsigned int) len && len > 6) {    
    //  printf("%s %d %d\n", m, strlen(m), len);
    m2[len-2] = m2[len-3] = m2[len-4] = '.';
    m2[len-5] = ' ';
    m2[len-1] ='\0';
    // printf("out %s %d %d\n", m, strlen(m), len);
  }
  strcopyN(M, m2, len);
  if (PL >= PL_ERRORS) { 
    PRINTF("error code %d [%s]\n", err, m2);
  }
}

void errorMSG(int err, char* m) {
  errorMSG(err, m, LENERRMSG, false);
}

void FinalErrorMSG(int err, char* m) {
  errorMSG(err, m, LENERRMSG, true);
}

void ErrorStop(int err) {
  char m[LENERRMSG];
  errorMSG(err, m);
  ERR(m);
}



int checkOK(cov_model VARIABLE_IS_NOT_USED *cov){
   return NOERROR;
}

int checkMissing(cov_model *cov){
  if (cov->calling == NULL) ERR("missing may not be called by the user");
  char S[100];
  cov_model *prev=cov->calling;
  sprintf(S, "'%s' does have not enough submodels", NICK(prev));
  ERR(S);
  return ERRORFAILED; // damit compiler keine Warnung bringt
}

int checkNotOK(cov_model VARIABLE_IS_NOT_USED *cov){
   return ERRORFAILED;
}

void ScaleOne(double *x, cov_model VARIABLE_IS_NOT_USED *cov, double *v){ 
  *v = *x <= 0.05 ? 1.0 : RF_NA;
} 


sortsofparam paramtypeAny(int VARIABLE_IS_NOT_USED k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) { return ANYPARAM; }

double *EinheitsMatrix(int dim) {
  // Einheitsmatrizen
  double *mem;
  if ((mem = (double*) CALLOC(dim * dim, sizeof(double))) != NULL) {
    int d;
    for (d=0; d<dim; d+=dim+1) mem[d] = 1.0;
  }
  return mem;
}




bool CheckListcov(){
  for (int nr=0; nr<currentNrCov; nr++) {     
    cov_fct *C = CovList + nr; // nicht gatternr
    // printf("name = %s %d\n", C->nick, nr);
    if (!(C->Typi[0] != UndefinedType || C->TypeFct != NULL)) return false;
    
    for (int k=0; k<C->kappas; k++) {
      if (C->kappanames[k][0] == ONEARGUMENT_NAME
	  && C->kappanames[k][1] >= '0'
	  && C->kappanames[k][1] <= '9') return false;
    }
  }
  return true;
}
   
   

void InitModelList() {
  assert(currentNrCov=-1); // otherwise something went wrong with the call
  assert(MODEL_MAX == 21); // otherwise change REGNAMES
  assert(ROLE_LAST == 14); // otherwise change ROLENAMES
  assert(OtherType == 15); // otherwise change TYPENAMES, 
  //                                           CAT_TYPENAMES[OtherType + 1]
  //                                           TypeConsistency in getNset
  assert(LAST_ISO == 17); // otherwise change ISONAMES
  assert(MAXMPPDIM <= MAXSIMUDIM); // ZERO
  // assert(CUTOFF_THEOR == 4);/* do not change this value as used in RFmethods.Rd */

  includeXport();

  int i;

  for (i=0; i<MAXSIMUDIM; i++) ZERO[i] = 0.0;

  for (i=0; i<MAXPARAM; i++) sprintf(STANDARDPARAM[i], "k%d", i+1);
  for (i=0; i<MAXSUB; i++) sprintf(STANDARDSUB[i], "u%d", i+1);
  /* ja nicht setzen !! macht riesen aerger, da RF opt ions InitModel
    nicht aufruft:
    for (i=0; i<MAXUNITS; i++) {
    strcpy(GLOBAL.general.newunits[i], "");
    strcpy(GLOBAL.general.curunits[i], "");
    strcpy(GLOBAL.general.varunits[i], "");
  }
  */

  //  assert (KEY == NULL);
  //  KEY = (cov_model*) MALLOC(sizeof(cov_model **) * (MODEL_MAX+1));

   // init models
  for (i=0; i<=MODEL_MAX; i++) {
    KEY[i] = NULL; 
    MEM_NAS[i] = -1;
  }

  if (CovList!=NULL) {
    PRINTF("List of covariance functions looks already initiated.\n"); 
    return;
  }
  CovList = (cov_fct*) MALLOC(sizeof(cov_fct) * (MAXNRCOVFCTS+1));
  // + 1 is necessary because of COVINFO_NULL that uses the last + 
  currentNrCov = 0;

  // *******************
  // **** trend-models ****
  // *******************
 
  MIXEDEFFECT = 
    IncludeModel("mixed", TrendType, 0, 1, 5, kappamixed,
		 XONLY, PREVMODELI,  // todo !!
		 checkmixed, rangemixed, PREF_NOTHING,
		 false, SUBMODEL_DEP, SUBMODEL_DEP,
		 SUBMODEL_DEP, NOT_MONOTONE);
  make_internal();
  // if element is negative and SpaceEffect then only the covariance of the subsequent model is returned (as if 'X' were not given)
  kappanames(MIXED_X_NAME, LISTOF+REALSXP, MIXED_BETA_NAME, REALSXP,
	     "coord", REALSXP, "dist", REALSXP, "dim", INTSXP);
  change_sortof(MIXED_BETA, TRENDPARAM); 
  // here: keep_typeof(0, RandomType);
  subnames("cov");
  addCov(mixed, NULL, NULL);
  addCov(mixed_nonstat);
  RandomShape(0, initmixed, domixed);
  addReturns(NULL, NULL, covmatrix_mixed, iscovmatrix_mixed,
	     NULL, NULL, NULL, NULL, NULL, NULL);
  //  MLEMIXEDEFFECT = addFurtherCov(MLEmixed, ErrCov); 
  // addCov(MLEmixed_nonstat);


  pref_type ptrend = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5};
        //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  TREND = IncludeModel("trend", TrendType,  0, 0, 6,
		       kappatrend, 
		       XONLY, PREVMODELI,
		       checktrend, 
		       rangetrend,
		       ptrend,
		       false, PARAM_DEP, INFDIM, false, NOT_MONOTONE);
  //  kappanames("mean", REALSXP, "plane", REALSXP, "polydeg",		
  //             INTSXP, "polycoeff",					
  //REALSXP, "arbitrark2xyfct", CLOSXP, "fctcoeff", REALSXP);	
 kappanames("mean", REALSXP, INTERNAL_PARAM, REALSXP, INTERNAL_PARAM, INTSXP, 
	    INTERNAL_PARAM, REALSXP, INTERNAL_PARAM, CLOSXP,
	    INTERNAL_PARAM, REALSXP);
 change_sortof(TREND_MEAN, TRENDPARAM);
 change_typeof(TREND_MEAN, ShapeType);
 addCov(trend, NULL, NULL);
  //  addCov(trend_nonstat);

  // *******************
  // **** RO-models ****
  // *******************

  pref_type pGatter=  {5, 0, 0,  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  FIRST_TRAFO=  ISO2ISO = // 2
    IncludeModel("#",  OtherType, 1, 1, 0, NULL, PREVMODELD, PREVMODELI,
		 checkNotOK, NULL, pGatter, true, SUBMODEL_DEP,
		 SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  addCov(iso2iso, D_2, DD_2, inverse2, nonstatinverse2);
  addlogCov(logiso2iso, NULL, nonstat_loginverse2);
  RandomShape(INFTY, struct2, init2, do2, dorandom2, true, true, false); 
  addReturns(NULL, NULL, covmatrixS,  iscovmatrixS, 
  	     NULL, NULL, NULL, NULL, NULL, NULL);
    
  SP2SP = addFurtherCov(spiso2spiso, D_2, DD_2); // 3
  addlogCov(logspiso2spiso);
 
  SP2ISO = addFurtherCov(spacetime2iso, D_2, DD_2); // 4
  addlogCov(logspacetime2iso);

  S2ISO = addFurtherCov(Stat2iso, ErrCov); // 5
  addCov(Nonstat2iso);// 
  addlogCov(logStat2iso, logNonstat2iso, NULL);

  assert(CovList[S2ISO].Init != NULL);
  assert(S2ISO == 5);

  S2SP = addFurtherCov(Stat2spacetime, ErrCov);// 6
  addCov(Nonstat2spacetime);// 
  addlogCov(logStat2spacetime, logNonstat2spacetime, NULL);
  
  S2S = addFurtherCov(Stat2Stat, ErrCov);// 7
  addCov(Nonstat2Stat);// 
  addlogCov(logStat2Stat, logNonstat2Stat, NULL);
  // printf("# %ld %ld %ld\n", Stat2Stat, CovList[currentNrCov-1].cov, Stat2iso);


  //print("%d\n", currentNrCov);assert(currentNrCov==8);

  SId = addFurtherCov(Stat2Stat, ErrCov);// 8
  addCov(Nonstat2Nonstat);// 
  addlogCov(logStat2Stat, logNonstat2Nonstat, NULL);
  assert(SId == 8);

  E2EIso = addFurtherCov(EarthIso2EarthIso, ErrCov);// 9
  addCov(NonstatEarth2EarthIso);// 
  addlogCov(logEarthIso2EarthIso, logNonstatEarth2EarthIso, NULL);
 
  E2E = addFurtherCov(Earth2Earth, ErrCov);// 10
  addCov(NonstatEarth2Earth);// 
  addlogCov(logEarth2Earth, logNonstatEarth2Earth, NULL);

  E2SphIso = addFurtherCov(EarthIso2SphereIso, ErrCov);// 11
  addCov(NonstatEarth2SphereIso);// 
  addlogCov(logEarthIso2SphereIso, logNonstatEarth2SphereIso, NULL);

  E2Sph = addFurtherCov(Earth2Sphere, ErrCov);// 12
  addCov(NonstatEarth2Sphere);// 
  addlogCov(logEarth2Sphere, logNonstatEarth2Sphere, NULL);


  Sph2SphIso = addFurtherCov(SphereIso2SphereIso, ErrCov);// 13
  addCov(NonstatSphere2SphereIso);// 
  addlogCov(logSphereIso2SphereIso, logNonstatSphere2SphereIso, NULL);

  Sph2Sph = addFurtherCov(Sphere2Sphere, ErrCov);// 14
  addCov(NonstatSphere2Sphere);// 
  addlogCov(logSphere2Sphere, logNonstatSphere2Sphere, NULL);

  LASTGATTER = Sph2Sph;
  
  EARTHKM2CART= // 15
      IncludeModel(">",  OtherType, 1, 1, 0, NULL, PREVMODELD, PREVMODELI,
		   checkEarth, NULL, pGatter, true, SUBMODEL_DEP,
		   4, SUBMODEL_DEP, SUBMODEL_DEP);
  addCov(EarthKM2CartStat, NULL, NULL);
  addCov(EarthKM2Cart);// 
  addlogCov(logEarthKM2CartStat, logEarthKM2Cart, NULL);
  
  EARTHMILES2CART = addFurtherCov(EarthMiles2CartStat, ErrCov);// 16
  addCov(EarthMiles2Cart);// 
  addlogCov(logEarthMiles2CartStat, logEarthMiles2Cart, NULL);


  FIRST_PLANE = 
     EARTHKM2GNOMONIC = addFurtherCov(Earth2GnomonicStat, ErrCov);// 16
  addCov(Earth2Gnomonic);// 
  addlogCov(logEarth2GnomonicStat, logEarth2Gnomonic, NULL);

  EARTHMILES2GNOMONIC =  CopyModel(">", EARTHKM2GNOMONIC);

  EARTHKM2ORTHOGRAPHIC = addFurtherCov(EarthKM2OrthogStat, ErrCov);// 18
  addCov(EarthKM2Orthog);// 
  addlogCov(logEarthKM2OrthogStat, logEarthKM2Orthog, NULL);

  EARTHMILES2ORTHOGRAPHIC = addFurtherCov(EarthMiles2OrthogStat, ErrCov);// 19
  addCov(EarthMiles2Orthog);// 
  addlogCov(logEarthMiles2OrthogStat, logEarthMiles2Orthog, NULL);
 
  LAST_PLANE = EARTHMILES2ORTHOGRAPHIC;
  LAST_TRAFO =  EARTHMILES2ORTHOGRAPHIC; 

  //  addFurtherCov(iso2iso_MLE, D_2, DD_2); // 12
  //  addFurtherCov(spiso2spiso_MLE, D_2, DD_2); // 13
  //  addFurtherCov(spacetime2iso_MLE, D_2, DD_2); // 14
  //  addFurtherCov(Stat2iso_MLE, ErrCov); // 15
  //  addCov(Nonstat2iso_MLE);// 
  //  addFurtherCov(Stat2spacetime_MLE, ErrCov);// 16
  //  int mleGatter = addFurtherCov(Stat2Stat_MLE, ErrCov);// 17
  //  addCov(Nonstat2Stat_MLE);// 
  //  if (mleGatter != S2S + (S2S - ISO2ISO + 1)) {
  //    ERR("mleGatter has the wrong number");
  //  }

  
  MISSING_COV =
    IncludePrim("missing", OtherType, 0, XONLY, SYMMETRIC,
		checkMissing,  NULL, INFDIM, true, NOT_MONOTONE);
  make_internal(); 

  NULL_MODEL =
    IncludePrim("null", UndefinedType, 1, XONLY, ISOTROPIC,
		checkOK,  rangeNullModel, INFDIM, true, NOT_MONOTONE);
  kappanames("type", INTSXP);
  addCov(0, NullModel, NullModel, NullModel, NullModel, NULL);
  RandomShape(INFTY, structOK, initOK, doOK, do_random_ok, false, false, false);
  make_internal(); 
  addTypeFct(TypeNullModel);

  // *******************
  // **** definite functions  ****
  // *******************

  SELECTNR =  // to do: replace by parameter in '+', selecting the 'type' or
    // 'models'
    IncludeModel("select", TcfType, 1, MAXSUB, 1, NULL,
		 PREVMODELD, PREVMODELI,
		 checkselect, rangeselect, PREF_ALL,
		 true, PARAM_DEP, INFDIM, SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("subnr", INTSXP);
  addCov(select, NULL, NULL); 
  
  addReturns(NULL, NULL, covmatrix_select, iscovmatrix_select, 
	     NULL, NULL, NULL, NULL, NULL, NULL);
   
  pref_type pplus =  {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 5, 5};
  //                  CE CO CI TBM Sp di sq Tr av n mpp Hy spf any
  PLUS = 
    IncludeModel("+", UndefinedType, 1, MAXSUB, 0, NULL, PREVMODELD, PREVMODELI,
		 checkplus, NULL, pplus, 
		 false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  // Achtung in covmatrix_plus wird SELECT_SUBNR verwendet!
  nickname("plus");
  addCov(plusStat, Dplus, DDplus, NULL, NULL);
  addCov(plusNonStat);
  addTBM(NULL, spectralplus);
  RandomShape(0, structplus, initplus, doplus, false, false, true);
  addReturns(NULL, NULL, covmatrix_plus, iscovmatrix_plus, 
	     NULL, NULL, NULL, NULL, NULL, NULL);
  addTypeFct(Typeplus);
  setptwise(pt_submodeldep);

 pref_type pmal =  {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //                 CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  MULT = IncludeModel("*", UndefinedType,  1, MAXSUB, 0, NULL,
		      PREVMODELD, PREVMODELI,
		      checkmal, NULL, pmal, false, SUBMODEL_DEP,
		      SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  nickname("mult");
  addCov(malStat, Dmal, NULL);
  addCov(malNonStat);
  addlogCov(logmalStat, logmalNonStat, NULL);
  //  RandomShape(structplusmal, initmal, domal, NULL);
  addTypeFct(Typemal);
  setptwise(pt_submodeldep);
 

  pref_type pS=  {5, 0, 0,  5, 5, 5, 5, 0, 0, 5, 0, 0, 1, 5};
  //              CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  DOLLAR = IncludeModel("$",  UndefinedType, // to do: tcftype durch einen allgemeinen Type ersetzen, da auch Trend dem "$" folgen kann. Z.Z. nicht moeglich.
			1, 1, 5, kappaS, // kappadollar,
			PREVMODELD, PREVMODELI, checkS, rangeS, pS,
			false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP,
			SUBMODEL_DEP);
  // do not change Order!!
  nickname("S");
  kappanames("var", REALSXP, "scale", REALSXP, "anisoT", REALSXP,
	     "Aniso", REALSXP, "proj", INTSXP);
  change_typeof(DAUSER, ShapeType);
  change_typeof(DSCALE, RandomOrShapeType);
  change_typeof(DVAR, RandomOrShapeType);
  subnames("phi");
  addTypeFct(TypeS);
  addCov(Siso, DS, DDS, D3S, D4S, inverseS, nonstatinverseS); // unterscheidung nur wegen der 
  //  geschwindigkeit, also Siso ist ein sehr haeufiger Spezialfall von Sstat
  addCov(Snonstat);
  addlogCov(logSiso, NULL, nonstat_loginverseS);
  addLocal(coinitS, ieinitS);  
  addTBM(tbm2S, NULL, spectralS);
  nablahess(nablaS, hessS);
  RandomShape(INFTY, structS, initS, doS, true, true, true);
  addReturns(NULL, NULL, covmatrixS, iscovmatrixS, 
  	     NULL, NULL, NULL, NULL, NULL, NULL);
  Taylor(RF_NA, RF_NA, RF_NA, RF_NA);
  TailTaylor(RF_NA, RF_NA, RF_NA, RF_NA);
  setptwise(pt_submodeldep);

   
  LASTDOLLAR = addFurtherCov(Sstat, DS, DDS); // 20.8.14 aus ErrCov (wieder)
  //                                        D2 gemacht
  addCov(Snonstat);
  addlogCov(logSstat, logSnonstat, NULL);
  RandomShape(INFTY, structS, initS, doS, true, true, false);
  

  // printf("%d\n",  currentNrCov); BUG;

  pref_type pPowS=  {5, 0, 0,  5, 5, 5, 5, 0, 0, 5, 0, 0, 1, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  POWER_DOLLAR = 
    IncludeModel("$power",  UndefinedType, // to do: tcftype durch einen allgemeinen Type ersetzen, da auch Trend dem "$" folgen kann. Z.Z. nicht moeglich.
		 1, 1, 3, NULL, // kappadollar,
		 PREVMODELD, PREVMODELI, checkPowS, rangePowS, pPowS,
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP,
		 SUBMODEL_DEP);
  // do not change Order!!
  nickname("Spower");
  kappanames("var", REALSXP, "scale", REALSXP, "pow", REALSXP);
  subnames("phi");
  addTypeFct(TypePowS);
  addCov(PowSstat, NULL, inversePowS); // unterscheidung nur wegen der 
  //  geschwindigkeit, also Siso ist ein sehr haeufiger Spezialfall von Sstat
  addCov(PowSnonstat);
  addlogCov(logSstat, logSnonstat, NULL);
  // addLocal(coinitS, ieinitS);  
  RandomShape(INFTY, structPowS, initPowS, doPowS, true, true, true);
  Taylor(RF_NA, RF_NA, RF_NA, RF_NA);
  TailTaylor(RF_NA, RF_NA, RF_NA, RF_NA);
  


  // at the origin: prespecified distribution
  // old RandomFields ave1,ave2
  IncludeModel("ave",  PosDefType, 1, 1, 3, kappa_ave, XONLY, SYMMETRIC,
	       checkave, rangeave, PREF_ALL, 
	       false, SCALAR, AveMaxDim, false, NOT_MONOTONE);
  kappanames("A", REALSXP, "z", REALSXP, "spacetime", INTSXP);
  addCov(ave, NULL, NULL);
  RandomShape(structAve, true);


  pref_type pbcw = {2, 5, 5, 5, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("bcw", VariogramType, 2, XONLY, ISOTROPIC,
	      checkbcw, rangebcw, pbcw,
	      SCALAR, INFDIM, false, NORMAL_MIXTURE); // todo part is even
  // LAPLACE
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(bcw, Dbcw, DDbcw, Inversebcw);
  addLocal(coinitbcw, ieinitbcw);
  AddVariant(PosDefType, ISOTROPIC);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);


  pref_type
    pbessel = {2, 0, 0,  0, 5, 3, 3, 0, 5, 0, 5, 0, 0, 5};
  //            CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("bessel",  PosDefType, 1, XONLY, ISOTROPIC, 
	      checkBessel, rangeBessel,
	      pbessel, SCALAR, INFDIM, false, NOT_MONOTONE);
  kappanames("nu", REALSXP);
  addCov(Bessel, NULL, NULL);
  addTBM(initBessel, spectralBessel);	       
   
  IncludeModel("bigneiting", PosDefType, 0, 0, 8, kappa_biGneiting, XONLY,
	       ISOTROPIC, checkbiGneiting, rangebiGneiting, PREF_ALL, 
	       false, 2, PARAM_DEP, true, NOT_MONOTONE);
  addCov(biGneiting, DbiGneiting, DDbiGneiting, NULL, NULL);
  kappanames("kappa", INTSXP,
	     "mu", REALSXP,
	     "s", REALSXP, "sred12", REALSXP,
	     "gamma", REALSXP,
	     "cdiag", REALSXP, "rhored", REALSXP, "c", REALSXP);
  change_sortof(GNEITING_S, SCALEPARAM);
  change_sortof(GNEITING_CDIAG, VARPARAM);
  change_sortof(GNEITING_C, DONOTRETURNPARAM);
  change_sortof(GNEITING_MU, CRITICALPARAM);
  change_sortof(GNEITING_GAMMA, CRITICALPARAM);
  RandomShape(0, struct_failed, initbiGneiting, do_failed, false, true, false);
   

  pref_type
    pbernoulli = {5, 0, 0,  0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //             CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("bernoulli", TcfType, 1, 1, 3, NULL, PREVMODELD, PREVMODELI,
	       checkbinary, rangebinary, pbernoulli,
	       false, SCALAR, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  kappanames("threshold", REALSXP, "correlation", INTSXP, "centred", INTSXP);
  addCov(binary, NULL, NULL);
  

  IncludeModel("biWM",  PosDefType, 0, 0, 8, kappa_biWM, XONLY, ISOTROPIC,
	       checkbiWM2, rangebiWM2, PREF_ALL,
	       false, 2, INFDIM, false, NOT_MONOTONE);
  nickname("biwm");
  addCov(biWM2, biWM2D,  ErrInverse);
  kappanames("nudiag", REALSXP, "nured12", REALSXP, 
	     "nu", REALSXP, // or lower triangle
	     "s", REALSXP,  // lower triangle definition
	     "cdiag", REALSXP, "rhored", REALSXP,
	     "c", REALSXP,  // or lower triangle
	     "notinvnu", INTSXP);
  change_sortof(BInudiag, CRITICALPARAM );
  change_sortof(BInured, CRITICALPARAM );
  change_sortof(BInu, DONOTRETURNPARAM);
  change_sortof(BIc, DONOTRETURNPARAM);
  change_sortof(BIs, SCALEPARAM );
  change_sortof(BIcdiag, VARPARAM);
  RandomShape(0, struct_failed, initbiWM2, do_failed, false, true, false);
 

  pref_type
    pbrownresnick = {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  BROWNRESNICK =
    IncludeModel("brownresnick", TcfType, 1, 1, 0, NULL, XONLY, PREVMODELI,
		 checkbrownresnick, NULL , pbrownresnick, false,
		 SCALAR, SUBMODEL_DEP, false, SUBMODEL_DEP);
  addCov(brownresnick, Dbrownresnick, DDbrownresnick, D3brownresnick, 
	 NULL, NULL);
  RandomShape(0, struct_brownresnick, init_brownresnick, do_brownresnick);
  //  Taylor(0, 0, 0, 0, 0, 0, 0, 0);

  
  IncludeModel("br2bg",  PosDefType, 1, 1, 0, XONLY, PREVMODELI, 
	       check_BR2BG, NULL, PREF_ALL, SUBMODEL_DEP, false, SUBMODEL_DEP);
  addCov(BR2BG, NULL, NULL);


  IncludeModel("br2eg", PosDefType, 1, 1, 0,  XONLY, PREVMODELI, 
	       check_BR2EG, NULL, PREF_ALL, SUBMODEL_DEP, false, SUBMODEL_DEP);
  addCov(BR2EG, NULL, NULL);
 
 
  pref_type pcauchy=  {2, 5, 0,  2, 0, 4, 0, 0, 0, 0, 0, 0, 0, 5};
  //                   CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("cauchy", PosDefType, 1, XONLY, ISOTROPIC, 
	      checkCauchy, rangeCauchy, pcauchy, 
	      SCALAR, INFDIM, false, NORMAL_MIXTURE);
  kappanames("gamma", REALSXP);
  addCov(Cauchy, DCauchy, DDCauchy, InverseCauchy);
  //  addlogCov(logCauchy);
  addTBM(TBM2Cauchy);
  addLocal(coinitCauchy, NULL);
  addGaussMixture(DrawMixCauchy, LogMixDensCauchy);
	       
  //  pref_type pctbm={2, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //     //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  // IncludePrim("cauchytbm", PosDefType,  3, XONLY, ISOTROPIC, checkOK,
  //	      rangeCauchytbm, pctbm, SCALAR, INFDIM, false);
  // kappanames("alpha", REALSXP, "beta", REALSXP, "gamma", REALSXP);
  // addCov(Cauchytbm, DCauchytbm, InverseCauchy); // scale not correct, but
  // should be an approximation that is good enough
	      
  IncludePrim("circular",  TcfType, 0, XONLY, ISOTROPIC,
	      checkOK, NULL, 2, false, GNEITING_MON);
  addCov(circular, Dcircular, ScaleOne);
  RandomShape(structcircular);


  // IncludePrim("cone", PosDefType,  3, XONLY, ISOTROPIC, checkcone, rangecone);
  //  kappanames("r", REALSXP, "socle", REALSXP, "height", REALSXP);
  // RandomShape(init_cone, mppget_cone, sd_standard, MPP_POISS);

  IncludePrim("CDeWijsian",  VariogramType, 2, NULL, XONLY, ISOTROPIC, 
	      checkdewijsian,  rangeDeWijsian, PREF_NOTHING, 
	      SCALAR, INFDIM, false, MONOTONE); 
  nickname("cdewijs");
  make_internal();
  kappanames("alpha", REALSXP, "range", REALSXP);
  addCov(DeWijsian, NULL, NULL, InverseDeWijsian); 
  
 
   // same as RMcovariate, ausser dass RMcovariate interpoliert
    // und CONSTANT mehrere Saetze von covariaten erlaubt??!!
  IncludeModel("constant", NegDefType, 0, 0, 1, kappaconstant, XONLY,
	       PREVMODELI,
		 //  PREVMODELD, PREVMODELI, 
		 //wegen Variogramm berechnung in stat. Fall
		 checkconstant, rangeconstant, PREF_ALL,
		 false, PARAM_DEP, INFDIM, false, SUBMODEL_DEP);
  kappanames("M", REALSXP);  
  addCov(constant, NULL, NULL);
  addCov(nonstatconstant);
  //  AddVariant(TcfType, PREVMODELI);
  //  AddVariant(ShapeType, EARTH_ISOTROPIC);
  //  AddVariant(TcfType, EARTH_ISOTROPIC);
  //  AddVariant(ShapeType, SPHERICAL_ISOTROPIC);
  //  AddVariant(TcfType, SPHERICAL_ISOTROPIC);



  COVARIATE = // intern ok
    IncludeModel("covariate", ShapeType, 0, 1, 5, kappa_covariate,
		 XONLY, PREVMODELI,
		 checkcovariate, rangecovariate, PREF_NOTHING, 
		 INTERN_SHOW, PARAM_DEP, INFDIM-1, false, NOT_MONOTONE);
  subnames("norm");
  kappanames(COVARIATE_C_NAME, LISTOF + REALSXP, COVARIATE_X_NAME, VECSXP, 
	     "raw", INTSXP, COVARIATE_ADDNA_NAME, INTSXP,
	     INTERNAL_PARAM, REALSXP);
  change_sortof(COVARIATE_X, DONOTVERIFYPARAM);
  change_sortof(COVARIATE_FACTOR, TRENDPARAM);
  addCov(covariate, NULL, NULL);
  setptwise(pt_paramdep);
  AddVariant(TrendType, PREVMODELI);

  pref_type pfix={0, 0, 0,  0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //              CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("fixcov", PosDefType, 0, 1, 3, kappa_fix,
	       KERNEL, SYMMETRIC, checkfix, rangefix, pfix, 
	       INTERN_SHOW, PARAM_DEP, INFDIM-1, false, NOT_MONOTONE);
  subnames("norm");
  kappanames("M", LISTOF + REALSXP, "x", VECSXP,  "raw", INTSXP);
  change_sortof(COVARIATE_X, DONOTVERIFYPARAM);
  addCov(fix);
  setptwise(pt_paramdep);

  pref_type pcox={2, 0, 0,  0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //              CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("coxisham",  PosDefType, 1, 1, 3, kappa_cox, 
	       XONLY, ZEROSPACEISO, 
	       checkcox, rangecox, pcox,
	       false, SCALAR, CoxMaxDim, false, NOT_MONOTONE);
  kappanames("mu", REALSXP, "D", REALSXP, "beta", REALSXP);  
  addCov(cox, NULL, NULL);
  addTBM(initcox, spectralcox);
  nablahess(coxnabla, coxhess);

  IncludePrim("cubic",  TcfType, 0, XONLY, ISOTROPIC, 
	      checkOK, NULL, 3, false, MONOTONE);
  addCov(cubic, Dcubic, ScaleOne);
	       
  pref_type pcurl= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("curlfree",  PosDefType, 1, 1, 0, NULL, XONLY, SYMMETRIC,
	       checkdivcurl, NULL, pcurl,
	       false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  addCov(curl, NULL, NULL);
 
  pref_type plocal={5, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //            CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  CUTOFF =  
    IncludeModel("cutoff",  PosDefType, 1, 1,2, NULL, XONLY, ISOTROPIC,
		 check_co, range_co, plocal,
		 false, SCALAR, MAXCEDIM,  true, MONOTONE);
  kappanames("diameter", REALSXP, "a", REALSXP);  
  addCov(co, NULL, NULL);
  addCallLocal(alternativeparam_co);
 
  //  warning("siehe Primitive.cc/biWM: cutoff funktioniert nicht bei MLE, vereinheitlichung mit natsc und verbesserung von biWM notwendig");
 

  IncludeModel("dagum",  PosDefType, 0, 0, 3, NULL, XONLY, ISOTROPIC,
	       checkdagum, rangedagum, PREF_ALL, false, 1, INFDIM, false,
	       PARAM_DEP);
  kappanames("beta", REALSXP, "gamma", REALSXP, INTERNAL_PARAM, REALSXP);
  addCov(dagum, Ddagum, Inversedagum);
  RandomShape(0, struct_failed, initdagum, do_failed, false, true, false);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);
  setptwise(pt_posdef);


  IncludePrim("dampedcosine",  PosDefType, 1, XONLY, ISOTROPIC,
	      checkdampedcosine, rangedampedcosine, PARAM_DEP,
	      false, NOT_MONOTONE);
  nickname("dampedcos");
  kappanames("lambda", REALSXP);
  addCov(dampedcosine, Ddampedcosine, Inversedampedcosine);
  // addlogCov(logdampedcosine);

  /*
  pref_type pderivative= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("derivative", PosDefType, 1, 1, 1, kappamixed, XONLY, PREVMODELD,
	       checkderivative, rangeDerivative, pderivative, 
	       true, //false, 
	       SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("partial",INTSXP);
  addCov(derivative, NULL, NULL);
  */
   
  pref_type pdewijsian = {2, 5, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                     CE CO CI TBM Sp di sq Ma av n mpp Hy spf any

  //old version
  IncludePrim("DeWijsian", VariogramType,  1, XONLY, ISOTROPIC,
	      checkOK, rangedewijsian, pdewijsian,
	      SCALAR, INFDIM, false, MONOTONE);

  //copied from cauchy, must be checked

  //  IncludePrim("DeWijsian", VariogramType,  1, XONLY, ISOTROPIC,
  //        checkOK, rangedewijsian,   pdewijsian,   SCALAR, INFDIM, false, MONOTONE);


  //IncludePrim("DeWijsian", VariogramType,  1, XONLY,
   //           ISOTROPIC,
  //            checkOK, rangedewijsian,   pdewijsian,   SCALAR, INFDIM, false, MONOTONE);



  nickname("dewijsian");
  kappanames("alpha", REALSXP);
  addCov(dewijsian, Ddewijsian, DDdewijsian, D3dewijsian, D4dewijsian, Inversedewijsian);
  addLocal(coinitdewijsian, NULL);

 

  pref_type pdiv= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //               CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("divfree", PosDefType, 1, 1, 0, NULL, XONLY, SYMMETRIC, 
	       checkdivcurl, NULL, pdiv, 
	       false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  addCov(div, NULL, NULL);


 
  // epsC has been for internal reasons only ; essentially
  // the gencauchy model, except that 1 in the denominator 
  // is replaced by epsilon



  pref_type pepsC = {2, 0, 0, 5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("epsC",  PosDefType, 0, 0, 3, NULL, XONLY, ISOTROPIC,
	       checkepsC, rangeepsC, pepsC, 
	       false, SCALAR, INFDIM, false, NORMAL_MIXTURE);
  nickname("epscauchy");
  kappanames("alpha", REALSXP, "beta", REALSXP, "eps", REALSXP);
  addCov(epsC, DepsC, DDepsC, NULL, NULL);
  addlogCov(logepsC);

  IncludePrim("exponential", TcfType, 0, XONLY, ISOTROPIC,
	      checkexponential, NULL, INFDIM, false, COMPLETELY_MON);
  nickname("exp");
  addCov(0, exponential, Dexponential, DDexponential, Inverseexponential, NULL);
  addlogCov(logexponential, NULL, nonstatLogInvExp);
  addLocal(coinitExp, ieinitExp);
  addHyper(hyperexponential);
  // numerisches Ergebnis passt nicht !
  addGaussMixture(DrawMixExp, LogMixDensExp);
  addTBM(TBM2exponential, NULL, spectralexponential);
  RandomShape(1, initexponential, do_exp);
  Taylor(-1, 1.0, 0.5, 2.0);
  TailTaylor(1, 0, 1, 1);
  
  // operator, replacing any covariance fct C by exp(C) (componentwise)
  IncludeModel("Exp", 
	       PosDefType, 1, 1, 2, PREVMODELD, PREVMODELI, checkExp,
	       rangeExp, PREF_ALL, SUBMODEL_DEP, false, NOT_MONOTONE);
  nickname("exponential");
  kappanames("n", INTSXP, "standardised", INTSXP);
  addCov(Exp, DExp, DDExp, NULL, NULL);
   // setptwise(pt_paramdef);

  pref_type
    pextrgauss = {5, 0, 0,  0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //             CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("extremalgauss", TcfType, 1, 1, 0, NULL,
	       XONLY, PREVMODELI,
	       check_extrgauss, NULL, pextrgauss, false,
	       SCALAR, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  nickname("schlather");
  addCov(extrgauss, NULL, NULL);
  

  IncludePrim("FD", PosDefType,  1, XONLY, ISOTROPIC,
	      checkOK, rangeFD, 1, false, NOT_MONOTONE); 
  nickname("fractdiff");
  kappanames("a", REALSXP);
  addCov(FD, NULL, NULL);

  
   // same as RMcovariate, ausser dass RMcovariate interpoliert
    // und CONSTANT mehrere Saetze von covariaten erlaubt??!!

 IncludePrim("flatpower", VariogramType, 1, XONLY, ISOTROPIC, 
		checkoesting, rangeoesting, INFDIM, false,
		BERNSTEIN); // todo BERNSTEIN
  kappanames("alpha", REALSXP);
  addCov(2, oesting, Doesting, DDoesting, NULL, NULL);
  RandomShape(0, initoesting, do_statiso);
  Taylor(-1, 2, RF_NA, 4, RF_NA, 6);
  TailTaylor(-1, RF_NA, 0, 0);

  BROWNIAN = 
    IncludePrim("fractalB", VariogramType, 1, XONLY, ISOTROPIC, 
		checkfractalBrownian, rangefractalBrownian, INFDIM, false,
		BERNSTEIN); // todo BERNSTEIN
  nickname("fbm");
  kappanames("alpha", REALSXP);
  addCov(fractalBrownian, DfractalBrownian, DDfractalBrownian, 
	 D3fractalBrownian, D4fractalBrownian, 
	 InversefractalBrownian);
  addlogCov(logfractalBrownian);
  addLocal(NULL, ieinitBrownian);
  RandomShape(0, initfractalBrownian, do_statiso);
  Taylor(-1, RF_NA, 0, 0);
  TailTaylor(-1, RF_NA, 0, 0);
 

  IncludePrim("fractgauss", PosDefType, 1, XONLY, ISOTROPIC,
	      checkOK, rangefractGauss, 1, false, NOT_MONOTONE);
  kappanames("alpha", REALSXP);
  addCov(fractGauss, NULL, NULL);


  pref_type pgauss= {2, 0, 0, 5, 5, 5, 5, 5, 0, 0, 5, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  GAUSS = 
    IncludePrim("gauss",  PosDefType, 0, XONLY, ISOTROPIC,
		checkOK, NULL, pgauss,
		SCALAR, INFDIM, false, NORMAL_MIXTURE);
  addCov(Gauss, DGauss, DDGauss, D3Gauss, D4Gauss, InverseGauss);
  addlogCov(logGauss, NULL, nonstatLogInvGauss);
  addTBM(NULL, spectralGauss);
  RandomShape(INFTY, struct_Gauss, initGauss, do_Gauss, false, true, false);
  addGaussMixture(DrawMixGauss, LogMixDensGauss);
  Taylor(-1.0, 2.0);
  TailTaylor(1, 0, 1.0, 2.0);

  
  IncludePrim("genB", VariogramType, 2, XONLY, ISOTROPIC, 
	      checkOK, rangegenBrownian, INFDIM, false, MONOTONE);
  nickname("genfbm");
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(genBrownian, NULL, NULL, InversegenBrownian); 
  addlogCov(loggenBrownian);


  pref_type pgenc = {2, 0, 0, 2, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("gencauchy", PosDefType, 2, XONLY, ISOTROPIC,
	      checkgeneralisedCauchy, rangegeneralisedCauchy, pgenc,
	      SCALAR, INFDIM, false, PARAM_DEP); // todo part is even
  // LAPLACE
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(generalisedCauchy, DgeneralisedCauchy, DDgeneralisedCauchy,
	 InversegeneralisedCauchy);
  addlogCov(loggeneralisedCauchy);
  addLocal(coinitgenCauchy, ieinitgenCauchy);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);
  setptwise(pt_posdef);


  IncludePrim("gengneiting",  PosDefType, 2, XONLY, ISOTROPIC, 
	      checkgenGneiting, rangegenGneiting, INFDIM-1, true,
	      MONOTONE); // GNEITING_MON ??
  // not INFDIM, also not normalscale mixture and alpha will be void
  kappanames("kappa", INTSXP, "mu", REALSXP);
  addCov(genGneiting, DgenGneiting, DDgenGneiting, //NULL, NULL,
	 ScaleOne, NULL);



  IncludePrim("gengneiting",  PosDefType, 2, XONLY, ISOTROPIC, 
	      checkgenGneiting, rangegenGneiting, INFDIM-1, true,
	      MONOTONE); // GNEITING_MON ??
  // not INFDIM, also not normalscale mixture and alpha will be void
  kappanames("kappa", INTSXP, "mu", REALSXP);
  addCov(genGneiting, DgenGneiting, DDgenGneiting, //NULL, NULL,
	 ScaleOne, NULL);

  GNEITING_INTERN =
    IncludeModel("gengneit_intern", PosDefType, 0, 0, 2, NULL, 
		 XONLY, ISOTROPIC,
		 checkgenGneiting, rangegenGneiting, PREF_ALL, 
		 true, SCALAR, PARAM_DEP, true, MONOTONE);  
  nickname("gengneiting");
  kappanames("kappa", INTSXP, "mu", REALSXP);
  addCov(Gneiting, DGneiting, DDGneiting, ScaleOne);


  IncludeModel("gneiting", PosDefType, 0, 0, 1, XONLY, ISOTROPIC,
	      checkGneiting, rangeGneiting, PREF_ALL, 
	      PARAM_DEP, true, MONOTONE);  // GNEITING_MON ??
  kappanames("orig", INTSXP);
  addCov(Gneiting, DGneiting, DDGneiting, ScaleOne);

  pref_type pfgennsst= { 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                 CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("gennsst",  PosDefType, 2, 2, 0, KERNEL, PREVMODELI,
	       checkgennsst, NULL, pfgennsst,
	       SUBMODEL_DEP, false, NOT_MONOTONE);
  addCov(nonstatgennsst);
  
  GENNSST_INTERN = // to do: why internal??
  IncludeModel("gennsst_intern", 
	       PosDefType, // stimmt nicht, aber egal
	       2, 2, 1, kappa_gennsst_intern, KERNEL, SYMMETRIC,
	       checkgennsst_intern, range_gennsst_intern, pfgennsst,
	       true, SCALAR, SUBMODEL_DEP, false, NOT_MONOTONE);
  nickname("gennsst");
  kappanames("A", REALSXP);
  addCov(nonstatgennsst_intern);

  /*
 pref_type phelmholtz= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("helmholtz",  PosDefType, 1, 1, 2, kappamixed, XONLY, SYMMETRIC,
         checkhelmholtz,rangeHelmholtz, phelmholtz,
           true, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("component",REALSXP,"Aniso",REALSXP);
  addCov(helmholtz, NULL, NULL);
  */

  pref_type phyper= {2, 0, 0, 3, 0, 4, 5, 0, 5, 0, 5, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("hyperbolic",  PosDefType, 3, XONLY, ISOTROPIC,
	      checkhyperbolic, rangehyperbolic, phyper,
	      SCALAR, INFDIM, false, NORMAL_MIXTURE);
  kappanames("nu", REALSXP, "lambda", REALSXP, "delta", REALSXP);
  addCov(hyperbolic, Dhyperbolic, NULL); // InversehyperbolicSq);
  addlogCov(loghyperbolic);


  IncludePrim("iacocesare",  PosDefType, 3, XONLY, SPACEISOTROPIC, 
	      checkOK, rangeIacoCesare, INFDIM, false, NOT_MONOTONE);
  nickname("iaco");
  kappanames("nu", REALSXP, "lambda", REALSXP, "delta", REALSXP);
  addCov(IacoCesare, NULL, NULL);
  

  IncludeModel("identity", UndefinedType, 1, 1, 1, NULL, PREVMODELD, PREVMODELI,
	       checkId, rangeId, PREF_ALL, 
	       false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  nickname("id");
  kappanames("vdim", INTSXP);
  addCov(IdStat, DId, DDId, IdInverse);
  addCov(IdNonStat);
  addTBM(TBM2Id, initId, spectralId);
  addLocal(coinitId, ieinitId);
  addTypeFct(Typesetparam);


  IncludePrim("kolmogorov",  VariogramType, 0, XONLY, VECTORISOTROPIC,
	      checkKolmogorov, NULL, 3, 3, false, NOT_MONOTONE);
  addCov(Kolmogorov, NULL, NULL);

  pref_type plgd1= {2, 0, 0, 5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("lgd1",  PosDefType, 2, NULL, XONLY, ISOTROPIC, 
	      checklgd1, rangelgd1, plgd1, 
	      SCALAR, PARAM_DEP, false, MONOTONE);
  nickname("lgd");
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(lgd1, Dlgd1, NULL); // Inverselgd1);



  // stimmt so nicht, siehe Gneiting 1998, on a alpha-sym multiv. char. fct:
  //  IncludeModel("lp", PosDefType,  1, 1, 1, XONLY, SYMMETRIC, 
  //	       checklp, rangelp,
  //	       (pref_type) {5, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5}
  //	       //          CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  //             );
  //  kappanames("p", REALSXP);
  //  addCov(lp, NULL, NULL);
 
  IncludeModel("mastein",  PosDefType, 1, 1, 2, XONLY, SPACEISOTROPIC, 
	       check_MaStein, range_MaStein, PREF_ALL, 
	       SUBMODEL_DEP, false, NOT_MONOTONE);
  kappanames("nu", REALSXP, "delta", REALSXP);
  addCov(MaStein, NULL, NULL);

    
  IncludeModel("ma1", PosDefType,  1, 1, 2, XONLY, SYMMETRIC,
	       checkma1, rangema1, PREF_ALL, 
	       SUBMODEL_DEP, false, NOT_MONOTONE);
  nickname("ma");
  kappanames("alpha", REALSXP, "theta", REALSXP);
  addCov(ma1, NULL, NULL);


  IncludeModel("ma2",  PosDefType, 1, 1, 0, XONLY, SYMMETRIC,
	       checkma2, NULL, PREF_ALL, SUBMODEL_DEP, false, NOT_MONOTONE);
  nickname("intexp");
  addCov(ma2, NULL, NULL);


  IncludeModel("M",  PosDefType, 1, 1, 1, kappaM, PREVMODELD, PREVMODELI,
	       checkM, rangeM, PREF_ALL,
	       false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  nickname("matrix");
  kappanames("M", REALSXP);
  addCov(Mstat, NULL, NULL);
  addCov(Mnonstat);
  add_sortof(sortof_M);

	  
  IncludeModel("matern", PosDefType, 0, 0, 2, PREVMODELD, PREVMODELI, 
	       checkMatern, rangeWM, PREF_ALL, INFDIM, false, SUBMODEL_DEP);
  kappanames("nu", REALSXP, "notinvnu", INTSXP);
  change_sortof(WM_NU, CRITICALPARAM); 
  change_typeof(WM_NU, RandomOrShapeType);
  addCov(Matern, DMatern, DDMatern, D3Matern, D4Matern, InverseMatern);
  addCov(NonStMatern);
  addlogCov(logMatern, logNonStMatern, NULL);
  addTBM(initMatern, spectralMatern);
  addLocal(coinitWM, ieinitWM);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);
  setptwise(pt_posdef);
 

  //  addGaussMixture(DrawMixWM, LogMixDensWM);

  IncludeModel("mqam", PosDefType,
  	       2, 10, 1, kappamqam, XONLY, SYMMETRIC,
  	       checkmqam, rangemqam, PREF_ALL, 
	       false, PARAM_DEP, SUBMODEL_DEP, false, NOT_MONOTONE);
  subnames("phi");
  kappanames("theta", REALSXP);
  change_sortof(QAM_THETA, CRITICALPARAM); 
  addCov(mqam, NULL, NULL);
 

  NATSC_INTERN =NATSC_USER = 
    IncludeModel("natsc", TcfType,  1, 1, 0, NULL, XONLY, ISOTROPIC,
		 checknatsc, NULL, PREF_ALL,
		 false, 1, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  // nie einen Parameter !
  addCov(natsc, Dnatsc, DDnatsc, Inversenatsc);
  addLocal(coinitnatsc, ieinitnatsc);
  addTBM(tbm2natsc, initnatsc, spectralnatsc);
  setptwise(pt_submodeldep);

  // NATSC_INTERN = CopyModel("natsc_intern", NATSC_USER);
  // make_internal();



  pref_type pfnsst= { 4, 0, 0, 2, 0, 5, 4, 0, 0, 0, 0, 0, 0, 5};
  //                 CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("nsst",  PosDefType, 2, 2, 1, XONLY, SPACEISOTROPIC,
	       checknsst, rangensst, pfnsst,
	       SUBMODEL_DEP, false, NOT_MONOTONE);
  subnames("phi", "psi");
  kappanames("delta", REALSXP); 
  change_sortof(NSST_DELTA, CRITICALPARAM);
  addCov(nsst, Dnsst, NULL);
  addTBM(TBM2nsst);
  setptwise(pt_posdef);

  //  IncludePrim("nsst2", 7, checknsst2, SPACEISOTROPIC, 
  //		   rangensst2);
  //  addCov(nsst2, Dnsst2, NULL);
  //  addTBM(NULL, NULL /* TBM3nsst2 */);

  pref_type pfnugget= { 4, 0, 0, 0, 0, 4, 4, 0, 0, 5, 0, 0, 0, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  NUGGET  = 
    IncludeModel("nugget",  TcfType, 0, 0, 2, NULL, XONLY, ISOTROPIC,
		 check_nugget, range_nugget, pfnugget, 
		 false, PREVMODEL_DEP, INFDIM, true, MONOTONE);
  kappanames("tol", REALSXP, "vdim", INTSXP);
  change_sortof(NUGGET_TOL, FORBIDDENPARAM);
  change_sortof(NUGGET_VDIM, FORBIDDENPARAM);
  addCov(nugget, NULL, Inversenugget);
  addReturns(NULL, NULL, covmatrix_nugget, iscovmatrix_nugget, 
	     NULL, NULL, NULL, NULL, NULL, NULL);
  AddVariant(TcfType, EARTH_ISOTROPIC);
  AddVariant(TcfType, SPHERICAL_ISOTROPIC);
 
  IncludeModel("parsWM", PosDefType, 0, 0, 1, kappa_parsWM, 
	       XONLY, ISOTROPIC,
	       checkparsWM, rangeparsWM, PREF_ALL,
	       false, PARAM_DEP, INFDIM, false, NOT_MONOTONE);
  nickname("parswm");
  addCov(parsWM, parsWMD, NULL);
  kappanames("nudiag", REALSXP);
  change_sortof(PARSnudiag, CRITICALPARAM);
  setptwise(pt_posdef);

  IncludePrim("penta", PosDefType, 0, XONLY, ISOTROPIC,
	      checkOK, NULL, 3, true, MONOTONE);
  addCov(penta, Dpenta, ScaleOne);

  IncludePrim("power", PosDefType,  1, XONLY, ISOTROPIC, 
	      checkpower, rangepower, INFDIM-1, true, MONOTONE);
  nickname("askey");
  kappanames("alpha", REALSXP);
  addCov(power, Dpower, ScaleOne);	
  AddVariant(TcfType, ISOTROPIC);
  //  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);
 
  POW = IncludeModel("Pow", PosDefType, 1, 1, 1, PREVMODELD, PREVMODELI,
	       checkPow, rangePow, PREF_ALL, SUBMODEL_DEP, false, NOT_MONOTONE);
  nickname("power");
  addCov(Pow, DPow, DDPow, InversePow); 
  kappanames("alpha", REALSXP);
  setptwise(pt_posdef);

  SHAPEPOW = addFurtherCov(shapePow, DshapePow);
  addCov(shapePow, DshapePow, DDshapePow, InverseShapePow); 
 
 
   pref_type pfprod= { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 5, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
   PROD = 
     IncludeModel("prod",  PosDefType, 1, 1, 0, NULL,
	       KERNEL, SYMMETRIC,  checkprod, NULL, pfprod,
	       false, SUBMODEL_DEP, INFDIM-1, false, NOT_MONOTONE);
  addCov(nonstatprod);
  subnames("phi");  
  AddVariant(PosDefType, SPHERICAL_SYMMETRIC);
  AddVariant(PosDefType, EARTH_SYMMETRIC);
  setptwise(pt_posdef);


  IncludeModel("qam",  PosDefType, 2, MAXSUB, 1, kappaqam, XONLY, ISOTROPIC,
	       checkqam, rangeqam, PREF_ALL, 
	       false, SCALAR, SUBMODEL_DEP, false, NOT_MONOTONE);
  kappanames("theta", REALSXP);
  subnames("phi");
  addCov(qam, NULL, NULL);
  change_sortof(QAM_THETA, CRITICALPARAM);
  

  IncludePrim("qexponential",  PosDefType, 1, XONLY, ISOTROPIC, 
	      checkOK, rangeqexponential, INFDIM, false, NOT_MONOTONE);
  nickname("qexp");
  kappanames("alpha", REALSXP);
  addCov(qexponential, Dqexponential, Inverseqexponential);

  SCATTER = IncludeModel("scatter",  PosDefType, 1, 1, 2, NULL, 
	       PREVMODELD, PREVMODELI, 
	       checkScatter, rangeScatter, PREF_ALL,
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("size", REALSXP, "max", INTSXP);
  addCov(Scatter, NULL, NULL);
 
  IncludeModel("schur",  PosDefType, 1, 1, 3, kappaSchur, 
	      PREVMODELD, PREVMODELI, 
	       checkSchur, rangeSchur, PREF_ALL,
	       false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("M", REALSXP, "diag", REALSXP, "rhored", REALSXP);
  addCov(Schurstat, NULL, NULL);
  addCov(Schurnonstat);
  add_sortof(sortof_M); 

 
  IncludeModel("shift", PosDefType, 1, 1, 1, kappashift, XONLY, SYMMETRIC,
	       checkshift, rangeshift, PREF_ALL, 
	       false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  nickname("delay"); // delayeffect
  addCov(shift, NULL, NULL);
  kappanames("s", REALSXP);

  IncludePrim("spherical", TcfType, 0, NULL, XONLY, ISOTROPIC, 
	      checkOK, NULL, 3, true, GNEITING_MON);
  nickname("spheric");
  addCov(spherical, Dspherical, DDspherical, ScaleOne);
  addTBM(TBM2spherical);
  RandomShape(1, structspherical, initspherical, dospherical,
	      false, true, false);
  Taylor(-3.0, 1.0, 0.5, 3.0);


  IncludePrim("stable",  PosDefType, 1, XONLY, ISOTROPIC, 
	      checkstable, rangestable, INFDIM, false, PARAM_DEP);
  kappanames("alpha", REALSXP);
  addCov(stable, Dstable, DDstable, Inversestable);
  addlogCov(logstable, NULL, nonstatLogInversestable);
  addLocal(coinitstable, ieinitstable);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);
  setptwise(pt_posdef);

  // SPACEISOTROPIC variant of stable -- used for testing purposes only
  //  IncludePrim("stableX", 1, checkOK, SPACEISOTROPIC, 
  //		  rangestable);
  //  addCov(stableX, DstableX, Inversestable);
  //  addTBM(NULL, NULL)


  STEIN =  
    IncludeModel("Stein", PosDefType,  1, 1, 2, NULL, XONLY, ISOTROPIC,
		 check_Stein, range_Stein, plocal,
		 false, SCALAR, MAXCEDIM, true, NOT_MONOTONE);
  nickname(METHODNAMES[CircEmbedIntrinsic]);
  kappanames("diameter", REALSXP, "rawR", REALSXP);  
  change_sortof(pLOC_DIAM, FORBIDDENPARAM);
  change_sortof(pLOC_A, FORBIDDENPARAM);
  addCov(Stein, NULL, NULL);
  addCallLocal(alternativeparam_Stein);
  //  RandomShape(struct_ce_approx, init_ce_approx, do_ce_approx);

 
  IncludePrim("steinst1",  PosDefType, 2, kappaSteinST1, XONLY, SYMMETRIC,
	      checkSteinST1, rangeSteinST1, INFDIM, false, NOT_MONOTONE);
  nickname("stein");
  kappanames("nu", REALSXP, "z", REALSXP);
  addCov(SteinST1, NULL, NULL);
  addTBM(initSteinST1, spectralSteinST1);
 

  IncludeModel("stp", PosDefType, 1, 2, 3, kappa_stp, KERNEL, SYMMETRIC,
	       checkstp, rangestp, PREF_ALL,
	       false, SCALAR, StpMaxDim, false, NOT_MONOTONE);
  addCov(stp);
  kappanames("S", REALSXP, "z", REALSXP, "M", REALSXP);
  change_typeof(STP_S, ShapeType);
  RandomShape(structStp, true);
  subnames("xi", "phi"); // H ueberall wo U-x steht. dort U-H(x)
  //                           gedoppelte immer zum Schluss!
  
  TBM_OP = // old RandomFields tbm2, tbm3
    IncludeModel("tbm",  PosDefType, 1, 1, 3, NULL, XONLY, PREVMODELI,
		 checktbmop, rangetbmop, PREF_ALL,
		 false, SUBMODEL_DEP, PARAM_DEP, PARAM_DEP, NOT_MONOTONE);
  kappanames("fulldim", INTSXP, "reduceddim", INTSXP, "layers", REALSXP); 
  addCov(tbm, NULL, NULL); // Dtbm, NULL); todo


  //  { int i; for (i=0; i<=Nothing; i++) printf("%d %d\n ", i, CovList[TREND].pref[i]); assert(false); }

   pref_type pfsum= { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
   IncludeModel("sum",  NegDefType, 0, 1, 0, NULL,
	       KERNEL, SYMMETRIC,  checksum, NULL, pfsum,
	       false, SUBMODEL_DEP, INFDIM-1, false, NOT_MONOTONE);
  addCov(nonstatsum);
  subnames("phi");  
  AddVariant(NegDefType, EARTH_SYMMETRIC);
  // AddVariant(NegDefType, SPHERICAL_COORDS);


  USER =
    IncludeModel("U", UndefinedType, 0, 0, 16, kappaUser, 
		 PREVMODELD, PREVMODELI,
		 checkUser, rangeUser, PREF_AUX, 
		 true,// FREEVARIABLE vorhanden. Muss extra in SpecialRMmodel.R
		 // definiert und nicht ueber generatemodels.R
		 PARAM_DEP, INFDIM, false, // per default.
		 NOT_MONOTONE);
  nickname("user");
  kappanames("type", INTSXP, "domain", INTSXP,  "isotropy", INTSXP,
	     "vdim", INTSXP, "beta", REALSXP, "varnames", INTSXP,
	     "fctn", LANGSXP, "fst", LANGSXP, "snd", LANGSXP,
	     "envir", ENVSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP, 
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP
	     //, "trd", LANGSXP
	     ); 
  // H ueberall wo U-x steht. dort U-H(x)
  addCov(User, DUser, DDUser, NULL, NULL);
  addCov(UserNonStat);
  addTypeFct(TypeUser);


  VECTOR = 
    IncludeModel("vector",  PosDefType, 1, 1, 2, NULL, XONLY, SYMMETRIC,
		 checkvector, rangevector, PREF_ALL, 
		 false, PARAM_DEP, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  addCov(vector, NULL, NULL);
  kappanames("a", REALSXP, "Dspace", INTSXP);
  addFurtherCov(vectorAniso, NULL); 


  pref_type pwave = {2, 0, 0, 0, 5, 4, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("wave", PosDefType, 0, XONLY, ISOTROPIC, 
	      checkOK, NULL, pwave, SCALAR, 3, false, NOT_MONOTONE);
  addCov(wave, NULL, Inversewave);
  addTBM(initwave, spectralwave);

 
  IncludeModel("whittle", PosDefType, 0,0, 2, PREVMODELD, PREVMODELI, 
	       checkWM, rangeWM, PREF_ALL, INFDIM, false, NORMAL_MIXTURE);
  kappanames("nu", REALSXP, "notinvnu", INTSXP);
  change_typeof(WM_NU, RandomOrShapeType);
  change_sortof(WM_NU, CRITICALPARAM); 
  addCov(Whittle, DWhittle, DDWhittle, D3Whittle, D4Whittle, InverseWhittle);
  addCov(NonStWhittle);
  addlogCov(logWhittle, logNonStWhittle, NULL);
  addTBM(initWhittle, spectralWhittle);
  addLocal(coinitWM, ieinitWM);
  addGaussMixture(DrawMixWM, LogMixDensW);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);
  setptwise(pt_posdef);



  SpherModelsinit( );
 
  // *******************
  // **** shape types  ****
  // *******************
  IncludeModel("angle", ShapeType, 0, 0, 4, kappa_Angle, XONLY, CARTESIAN_COORD,
	       checkAngle, rangeAngle, PREF_NOTHING, 
	       false, PARAM_DEP, INFDIM, false, NOT_MONOTONE);
  nickname("angle");
  addCov(Angle, NULL, NULL, invAngle, NULL);
  kappanames("angle", REALSXP, "lat.angle", REALSXP, "ratio", 
	     REALSXP, "diag", REALSXP);
 

  BALL= IncludePrim("ball",  ShapeType, 0,  NULL, 
		    XONLY, ISOTROPIC, checkOK, NULL, PREF_NOTHING, 
		    SCALAR, INFDIM-1, true, MONOTONE);
  addCov(ball, NULL, NULL, Inverseball);
  RandomShape(INFTY, struct_ball, init_ball, do_ball);   
  Taylor(1.0, 0.0);


  IncludePrim("EAxxA", ShapeType,  2, kappa_EAxxA, XONLY, CARTESIAN_COORD,
	      checkEAxxA, rangeEAxxA, PREF_NOTHING, 
	      PARAM_DEP, EaxxaMaxDim, false, NOT_MONOTONE);
  nickname("eaxxa");
  addCov(EAxxA, NULL, NULL);
  kappanames("E", REALSXP, "A", REALSXP);
  addSpecial(minmaxEigenEAxxA);


  IncludePrim("EtAxxA",  ShapeType, 3, kappa_EtAxxA, XONLY, CARTESIAN_COORD,
	      checkEtAxxA, rangeEtAxxA, 3, EaxxaMaxDim, false, NOT_MONOTONE);
  nickname("etaxxa");
  addCov(EtAxxA, NULL, NULL);
  kappanames("E", REALSXP, "A", REALSXP, "alpha", REALSXP);
  addSpecial(minmaxEigenEtAxxA);


  IDCOORD = // ACHTUNG falls internal geaendert, auch in KeyInfo.cc aendern
    IncludeModel("idcoord", ShapeType, 0, 0, 0, NULL, XONLY, PREVMODELI,
	       checkidcoord, NULL, PREF_NOTHING, 
	       true, PARAM_DEP, INFDIM, false, NOT_MONOTONE);
  addCov(idcoord, NULL, NULL);
 

  TRAFO =
    //   IncludeModel("trafo", UndefinedType, 0, 1, 1, NULL, PREVMODELD, PREVMODELI,
   IncludeModel("trafo", UndefinedType, 0, 1, 1, NULL, PREVMODELD, PREVMODELI,
	       checktrafo, rangetrafo, PREF_ALL, 
	       false, PARAM_DEP, INFDIM-1, false, NOT_MONOTONE);
  kappanames("new", INTSXP);
  change_typeof(TRAFO_ISO, NN2); // ISONAMES
  addCov(trafo, NULL, NULL);
  addCov(nonstattrafo);// 
  addlogCov(logtrafo,lognonstattrafo, NULL);
  subnames("phi");
  addTypeFct(Typetrafo);

 
  MULT_INVERSE =
    IncludeModel("mult_inverse", ShapeType, 1, 1, 0, NULL,
		 PREVMODELD, PREVMODELI,
		 checkmult_inverse, NULL, PREF_NOTHING,
		 true, SCALAR, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  addCov(mult_inverse, NULL, NULL); 
  addCov(mult_inverseNonstat);

 
  POLYGON =
    IncludeModel("polygon",  ShapeType, 0, 0, 1, NULL, XONLY, CARTESIAN_COORD, 
	       check_polygon, range_polygon, PREF_NOTHING,
	       false, SCALAR, 2, true, MONOTONE);
  kappanames("lambda", REALSXP);
  addCov(Polygon, NULL, NULL, Inversepolygon, InversepolygonNonstat);
  RandomShape(INFTY, struct_polygon, init_polygon, doOK);   
 

  IncludePrim("rational", ShapeType, 2, kappa_rational,
	      XONLY, CARTESIAN_COORD,
	      checkrational, rangerational, INFDIM, false, NOT_MONOTONE);
  addCov(rational, NULL, NULL);
  kappanames("A", REALSXP, "a", REALSXP);
  addSpecial(minmaxEigenrational);

  
  IncludePrim("rotat",  ShapeType, 2, kappa_rotat, XONLY, CARTESIAN_COORD,
	      checkrotat, rangerotat, PREF_NOTHING, SCALAR, 3, false,
	      NOT_MONOTONE);
  addCov(rotat, NULL, NULL);
  kappanames("speed", REALSXP, "phi", REALSXP);
  addSpecial(minmaxEigenrotat);


  IncludePrim("Rotat",  ShapeType, 1, kappa_Rotat, XONLY, CARTESIAN_COORD,
	      checkRotat, rangeRotat, PARAM_DEP, 3, false, NOT_MONOTONE);
  nickname("rotation");
  addCov(Rotat, NULL, NULL);
  kappanames("phi", REALSXP);


  RANDOMSIGN = 
    IncludeModel("sign",  ShapeType, 1, 1, 1, NULL, XONLY, PREVMODELI,
		 check_randomSign, range_randomSign, PREF_NOTHING,
		 false, SCALAR, SUBMODEL_DEP, SUBMODEL_DEP, NOT_MONOTONE);
  //nickname("");
  kappanames("p", REALSXP);
  addCov(randomSign, NULL, NULL, randomSignInverse, randomSignNonstatInverse);
  addlogCov(lograndomSign);
  RandomShape(1, struct_randomSign, init_randomSign, do_randomSign,
	      true, true, false); 
 

  SETPARAM = 
    IncludeModel("setparam", UndefinedType, 1, 1, 1, NULL, 
		 PREVMODELD, PREVMODELI,
		 checksetparam,  range_setparam, PREF_ALL, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  nickname("setparam");
  kappanames("performDo", INTSXP);
  addCov(setparamStat, Dsetparam, DDsetparam, D3setparam, D4setparam, 
	 Inverse_setparam, NonstatInverse_setparam);
  addCov(setparamNonStat);
  addTBM(NULL, spectralsetparam);
  RandomShape(INFTY, struct_failed, initsetparam, dosetparam,
	      false, false, true);
  addReturns(NULL, NULL, covmatrix_setparam, iscovmatrix_setparam, 
	     NULL, NULL, NULL, NULL, NULL, NULL);
  addTypeFct(Typesetparam);
 

  SHAPEAVE =
    IncludeModel("shape.ave", ShapeType, 1, 2, 3, kappa_ave,
		 XONLY, CARTESIAN_COORD,  
		 check_shapeave, rangeave, PREF_NOTHING,
		 true, SCALAR, INFDIM-1, true, NOT_MONOTONE);
  kappanames("A", REALSXP, "z", REALSXP, "spacetime", INTSXP);
  subnames("phi", "gauss");
  addlogCov(logshapeave);
  RandomShape(0, init_shapeave, do_shapeave, true);


  SHAPESTP = 
    IncludeModel("shape.stp",  ShapeType, 1, 4, 3, kappa_stp, KERNEL, 
		 CARTESIAN_COORD, check_shapestp, rangestp, PREF_NOTHING,
		 true, SCALAR, StpMaxDim, true, NOT_MONOTONE);
  kappanames("S", REALSXP, "z", REALSXP, "M", REALSXP); 
  addlogCov(logshapestp);
  subnames("xi2", "phi", "SXX", "gauss"); // hier gedoppeltes S nicht am Schluss 
  //                                       da auf checkstp zugreifend
  RandomShape(0, init_shapestp, do_shapestp);


  STROKORB_MONO =
    IncludeModel("m2r", ShapeType, 1, 1, 0, NULL, XONLY, ISOTROPIC,
		 checkstrokorb, NULL, PREF_NOTHING,
		 false, SCALAR, 3, SUBMODEL_DEP, SUBMODEL_DEP);
  addCov(strokorb, NULL, NULL); 
  RandomShape(1, structOK, init_strokorb, do_strokorb,
	      false, false, false);


#define balltest !true  
  IncludeModel("m3b", ShapeType, 1, 1, 
	       balltest ? 2 : 0, NULL, XONLY, ISOTROPIC,
	       checkstrokorbBall, 
	       balltest ? rangestrokorbball : NULL,/*for testing only*/
	       PREF_NOTHING,
	       false, SCALAR, 3, true, MONOTONE);
  if (balltest) kappanames("min", REALSXP, "max", REALSXP); 
  addCov(strokorbBallInner, NULL, NULL); 
  RandomShape(1, struct_strokorbBall, init_failed, do_failed, do_random_failed,
	      false, false, false);
  
 
  STROKORB_BALL_INNER = // !! inverse scale gegenueber paper
    IncludeModel("r3binner", ShapeType, 1, 1, 1, NULL,
		 XONLY, CARTESIAN_COORD,
		 check_strokorbBallInner, range_strokorbBallInner, PREF_AUX,
		 true, 1, 1, true, NOT_MONOTONE);
  kappanames("dim", INTSXP);
  addCov(strokorbBallInner, NULL, NULL);
  RandomShape(1, init_strokorbBallInner, do_strokorbBallInner);

  /*
    da wiederum gewichtet und zwar mit b^2 falls b die intensitaet.
    kann nicht in dichte function g(b) reingezogen werden, da
    b^2 g(b) nicht integrierbar. Stattdessen darf f (Dichte im Raum)
    nicht die Gleichverteilung sein, sondern bei grossen b um
    die zu simulierenden Punkte zusammenschrumpfen.
    Dabei nimmt man an, dass die Radien ein Vielfaches des mittleren
    Radius nicht ueberschreiten. Dies ist OK, da ungefaehr 
    exponentielles Abfallen der WK.
  */

  IncludeModel("mps", ShapeType, 1, 1, 0, NULL, XONLY,CARTESIAN_COORD,
		 checkstrokorbPoly, NULL, PREF_AUX,
		 false, SCALAR, 2, true, MONOTONE);
  addCov(strokorbPoly, NULL, NULL); 
  RandomShape(1, struct_strokorbPoly, init_failed, do_failed, do_random_failed,
	      false, false, false);
  

  TRUNCSUPPORT =
    IncludeModel("truncsupport", ShapeType, 
		 1, 1, 1, NULL, XONLY, PREVMODELI, checktruncsupport,
		 rangetruncsupport, PREF_NOTHING, false, SCALAR,
		 SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  kappanames("radius", REALSXP); // neg value == auto
  addCov(truncsupport, NULL, NULL, truncsupportInverse, StandardInverseNonstat);
  RandomShape(0, struct_truncsupport, init_truncsupport, do_truncsupport, false,
	      false, false);
  

  //////////////////////////////////////////////////
  // families of multivariate distribution; used in 

  // ACHTUNG!! addCov muss ganz zum Schluss !!

  ARCSQRT_DISTR = 
    IncludeModel("arcsqrt", RandomType, 0, 0, 1, NULL, 
		 DOMAIN_MISMATCH, ISO_MISMATCH,
		 checkOK, range_arcsqrt, PREF_AUX, 
		 true, 1, 1, false, MISMATCH); // to do: nicht mismatch,
 //     sondern monotone im eindimensionalen
  kappanames("scale", REALSXP);
  RandomShape(0, structOK, init_arcsqrt, do_arcsqrt);
  addCov(arcsqrtD, arcsqrtDlog, arcsqrtDinverse, 
	 arcsqrtP, NULL, arcsqrtQ, arcsqrtR, NULL);
 
 
  DETERM_DISTR = 
    IncludeModel("determ", RandomType, 0, 0, 1, kappa_determ, 
		 DOMAIN_MISMATCH, ISO_MISMATCH,
		 check_determ, range_determ, PREF_AUX,
		 false, SUBMODEL_DEP, INFDIM-1, SUBMODEL_DEP, MISMATCH);
  kappanames("mean", REALSXP);
  RandomShape(INFTY, structOK, init_determ, do_determ); 
  addCov(determD, determDlog, determDinverse, determP, determP2sided, determQ, 
	 determR, determR2sided);

 
  DISTRIBUTION = // FREEVARIABLE vorhanden. Muss extra in SpecialRMmodel.R
		 // definiert und nicht ueber generatemodels.R 
    IncludeModel("distr", RandomType, 0, 0, 16, kappa_distr, 
		 DOMAIN_MISMATCH, ISO_MISMATCH,
		 check_distr, range_distr, PREF_AUX,
		 true, PARAM_DEP, INFDIM-1, false, MISMATCH);
  kappanames("ddistr", LANGSXP, 
	     "pdistr", LANGSXP,
	     "qdistr", LANGSXP, 
	     "rdistr", LANGSXP,
	     "nrow", INTSXP,
	     "ncol", INTSXP,
	     "envir", ENVSXP, 
	     FREEVARIABLE, REALSXP, // wird nie verwendet -- Puffer fuer 
	     // einfachen Algorithmus
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP, 
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP
	     ); // 7 free ones are remaining !	     
  RandomShape(0, structOK, init_distr, do_distr_do);
  addCov(distrD, distrDlog, distrDinverse, distrP, distrP2sided, distrQ,
	 distrR, distrR2sided);

 
  GAUSS_DISTR = 
    IncludeModel("normal", RandomType, 0, 0, 3, kappa_gauss_distr, 
		 DOMAIN_MISMATCH, ISO_MISMATCH,
		 check_gauss_distr, range_gauss_distr, PREF_AUX, 
		 false, PARAM_DEP, INFDIM-1, false, MISMATCH);
  nickname("gauss");
  kappanames("mu", REALSXP, "sd", REALSXP, "log", INTSXP);
  RandomShape(INFTY, structOK, init_gauss_distr, do_gauss_distr);
  addCov(gaussD, gaussDlog, gaussDinverse, 
	 gaussP, gaussP2sided, gaussQ, gaussR, gaussR2sided);
 
  
  SET_DISTR = 
    IncludeModel("setDistr", RandomType, 1, 1, 1, NULL, 
		 DOMAIN_MISMATCH, ISO_MISMATCH,
		 check_setParam, range_setParam, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  subnames("to");
  kappanames("performDo", INTSXP);
  RandomShape(INFTY, structOK, init_setParam, do_setParam);
  addCov(setParamD, setParamDlog, setParamDinverse, 
	 setParamP, setParamP2sided, setParamQ,
	 setParamR, setParamR2sided);
  

  LOC =
    IncludeModel("loc", RandomType, 1, 1, 3, kappa_loc, 
		 DOMAIN_MISMATCH, ISO_MISMATCH, 
		 check_loc, range_loc, PREF_AUX, 
		 false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  kappanames("mu", REALSXP, "scale", REALSXP, "pow", REALSXP);
  RandomShape(INFTY, structOK, init_loc, do_loc);
  addCov(locD, locDlog, locDinverse, locP, locP2sided, locQ, locR, locR2sided);

  MCMC =
    IncludeModel("mcmc", RandomType, 1, 1, 6, kappa_mcmc, 
		 // ACHTUNG! Model kann auch ueber cov->q uebergeben werden.
		 //          dies vereinfacht die Verwendung von zufaelligen
		 //          Huetchen, da keine Parameter kopiert werden 
		 //          muesen, sondern direkt auf das Huetchen zugegriffen
		 //         
		DOMAIN_MISMATCH, ISO_MISMATCH,
		check_mcmc, range_mcmc, PREF_AUX, 
		false, PARAM_DEP,  INFDIM-1, true, MISMATCH);
  kappanames(distr[RECT_MCMC_N], INTSXP,  "sigma", REALSXP, 
	     "normed", INTSXP, "maxdensity", REALSXP, "rand.loc", INTSXP,
	     "gibbs", INTSXP);
  RandomShape(INFTY, structOK, init_mcmc, do_mcmc); 
  addCov(mcmcD, mcmcDlog, mcmcDinverse, mcmcP, mcmcP2sided, mcmcQ,
	 mcmcR, mcmcR2sided);

  RECTANGULAR =
    IncludeModel("rectangular", RandomType, 1, 1, 11, NULL, 
		 // ACHTUNG! Model kann auch ueber cov->q uebergeben werden.
		 //          dies vereinfacht die Verwendung von zufaelligen
		 //          Huetchen, da keine Parameter kopiert werden 
		 //          muesen, sondern direkt auf das Huetchen zugegriffen
		 //         
		 DOMAIN_MISMATCH, ISO_MISMATCH,
		 check_rectangular, range_rectangular, PREF_AUX, 
		 false, PARAM_DEP,  INFDIM-1, true, MISMATCH);
  kappanames(distr[RECT_SAFETY], REALSXP, distr[RECT_MINSTEPLENGTH], REALSXP,
	     distr[RECT_MAXSTEPS], INTSXP, distr[RECT_PARTS], INTSXP, 
	     distr[RECT_MAXIT], INTSXP, distr[RECT_INNERMIN], REALSXP,
	     distr[RECT_OUTERMAX], REALSXP, distr[RECT_MCMC_N], INTSXP,
	     "normed", INTSXP, "approx", INTSXP, "onesided", INTSXP
	     );
  RandomShape(INFTY, structOK, init_rectangular, do_rectangular); 
  addCov(rectangularD, rectangularDlog, rectangularDinverse, rectangularP, 
	 rectangularP2sided, rectangularQ, rectangularR, rectangularR2sided);


  SCALESPHERICAL = 
    IncludeModel("spheric", RandomType, 0, 0, 3, NULL,
		 DOMAIN_MISMATCH, ISO_MISMATCH,
		 check_RRspheric, range_RRspheric, PREF_AUX,
		 false, 1, 1, true, MISMATCH);
  kappanames("spacedim", INTSXP, "balldim", INTSXP, "R", REALSXP);
  RandomShape(INFTY, structOK, init_RRspheric, do_RRspheric);
  addCov(sphericD, sphericDlog, sphericDinverse, sphericP, NULL, sphericQ,
	 sphericR, NULL);
 
 
  UNIF = IncludeModel("unif", RandomType, 0, 0, 3, kappa_unif, 
		      DOMAIN_MISMATCH, ISO_MISMATCH,
		      check_unif, range_unif, PREF_AUX, 
		      false, PARAM_DEP,  INFDIM-1, true, MISMATCH);
  kappanames("min", REALSXP, "max", REALSXP, "normed", INTSXP);
  RandomShape(INFTY, structOK, init_unif, do_unif); 
  addCov(unifD, unifDlog, unifDinverse, unifP, unifP2sided, unifQ, 
	 unifR, unifR2sided);

  
  // -----------------------------
  // shape + locations 
  // they *take* all very detailed roles like ROLE_SMITH and pass
  // ROLE_MAXSTABLE to the submodel, in general
  // storage always pgs_storage !!
  MCMC_PGS =
    IncludeModel("MCMC_PGS", PointShapeType, 2, 2, 2, NULL, 
		 XONLY, CARTESIAN_COORD,
		 check_mcmc_pgs, range_mcmc_pgs, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  kappanames("normed", INTSXP,
	     "mcmc_n", INTSXP);
  subnames("shape", "loc");
  addlogCov(logpts_given_shape);
  RandomShape(SUBMODEL_DEP, struct_mcmc_pgs, init_mcmc_pgs, 
	      do_mcmc_pgs, do_random_failed, true, true, false); 


  PTS_GIVEN_SHAPE = 
    IncludeModel("ptsGivenShape", PointShapeType, 2, 2, 5, NULL, 
		 XONLY, CARTESIAN_COORD,
		 check_pts_given_shape, range_pts_given_shape, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  kappanames("density.ratio", REALSXP,  // stemming from gauss
	     "flat", INTSXP,
	     "infinitely_small", INTSXP,
	     "normed", INTSXP,
	     "isotropic", INTSXP
	     );
  subnames("shape", "loc");
  addCov(pts_given_shape, NULL, NULL);
  addlogCov(logpts_given_shape);
  RandomShape(SUBMODEL_DEP, struct_pts_given_shape, init_pts_given_shape, 
	      do_pts_given_shape, do_random_failed, true, true, false); 


  STANDARD_SHAPE = 
    IncludeModel("standardShape", PointShapeType, 1, 2, 0, NULL, 
		 XONLY, CARTESIAN_COORD,
		 check_standard_shape, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  subnames("shape");
  addCov(standard_shape, NULL, NULL);
  addlogCov(logstandard_shape);
  RandomShape(SUBMODEL_DEP, struct_standard_shape, init_standard_shape, 
	      do_standard_shape, do_random_failed, true, true, false); 


  STATIONARY_SHAPE = 
    IncludeModel("statShape", PointShapeType, 1, 1, 0, NULL, 
		 XONLY, CARTESIAN_COORD,
		 check_stationary_shape, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  subnames("shape");
  addCov(stationary_shape, NULL, NULL);
  addlogCov(logstationary_shape);
  RandomShape(SUBMODEL_DEP, struct_stationary_shape, init_stationary_shape, 
	      do_stationary_shape, do_random_failed, true, true, false); 


  pref_type pmppp =  {0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 5, 0, 10, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  MPPPLUS =
    IncludeModel("++",  PointShapeType, 1, MAXSUB, 1, kappamppplus, 
		 PREVMODELD, PREVMODELI, // CARTESIAN_COORD,
		 checkmppplus, rangempplus, pmppp, 
		 false, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("mppplus");
  kappanames("p", REALSXP);
  addCov(mppplus, NULL, NULL);


  // -----------------------------
  // Interfaces
  // -----------------------------
  
  COVFCTN =
    IncludeModel("Cov", InterfaceType, 1, 1, 0, NULL, XONLY, UNREDUCED, 
		 check_cov, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("cov");
  addCov(Cov, NULL, NULL);
  RandomShape(struct_cov); 

  COVMATRIX = 
    IncludeModel("CovMatrix", InterfaceType, 1, 1, 0, NULL, XONLY,
		 UNREDUCED, //UNREDUCED,ISOTROPIC dependening on whether
		 // distances are givenXONLY, UNREDUCED,
		 check_covmatrix, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("covmatrix");
  addCov(CovMatrix, NULL, NULL);
  RandomShape(struct_cov); 

  IncludeModel("Dummy", InterfaceType, 1, 1, 0, NULL, XONLY, UNREDUCED, 
	       check_dummy, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("dummy");
  addCov(Dummy, NULL, NULL);
  RandomShape(struct_dummy); 

  RFGET = 
    IncludeModel("get", InterfaceType, 1, 1, 2, NULL, XONLY, UNREDUCED,
		 check_RFget, range_RFget, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);  
  kappanames("up", INTSXP, "register", INTSXP);
  addCov(RFget, NULL, NULL);
  RandomShape(struct_RFget); 

 
  IncludeModel("Fctn", InterfaceType, 1, 1, 0, NULL, XONLY, UNREDUCED, 
		 check_fctn, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("fctn");
  addCov(Fctn, NULL, NULL);
  RandomShape(structOK); 

  IncludeModel("Distr", InterfaceType, 1, 1, 5, kappa_EvalDistr,
	       XONLY, UNREDUCED,
	       check_EvalDistr, range_EvalDistr, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("distr");
  kappanames("x", REALSXP, "q", REALSXP, "p", REALSXP, "n", REALSXP,
	     "dim", INTSXP);
  addCov(EvalDistr, NULL, NULL);
  RandomShape(struct_EvalDistr); 


  IncludeModel("Pseudovariogram", InterfaceType, 1, 1, 0, NULL, 
	       XONLY, UNREDUCED,
	       check_cov, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("pseudovariogram");
  addCov(Pseudovariogram, NULL, NULL);
  RandomShape(struct_variogram); 
 
  SIMULATE =
    IncludeModel("Simulate", InterfaceType, 1, 1, 3, NULL, 
	       XONLY, UNREDUCED, 
	       check_simulate, range_simulate, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("simulate");
  kappanames("checkonly", INTSXP, "setseed", LANGSXP, "env", ENVSXP);
  addCov(simulate, NULL, NULL);
  RandomShape(struct_simulate); 

  /* 
  DENSITY =
    IncludeModel("Density", InterfaceType, 1, 1, 3, NULL, 
	       XONLY, UNREDUCED, 
	       check_density, range_density, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("density");
  kappanames("log", INTSXP, "setseed", LANGSXP, "env", ENVSXP);
  addCov(density, NULL, NULL);
  RandomShape(struct_density); 
  */
  

  VARIOGRAM_CALL =
    IncludeModel("Variogram", InterfaceType, 1, 1, 0, NULL, 
		 XONLY, UNREDUCED,
		 check_vario, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  nickname("variogram");
  addCov(Variogram, NULL, NULL);
  RandomShape(struct_variogram); 


  LIKELIHOOD_CALL =
  IncludeModel("loglikelihood", InterfaceType, 1, 1, 4, kappalikelihood, 
	       XONLY, UNREDUCED, check_likelihood, range_likelihood, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  kappanames("data", LISTOF +  REALSXP, "estimate_variance", INTSXP,
	     "betas_separate", INTSXP, "ignore_trend", INTSXP);
  addCov(likelihood, NULL, NULL);
  RandomShape(struct_likelihood); 

  LINEARPART_CALL =
  IncludeModel("linearpart", InterfaceType, 1, 1, 0, NULL,
	       XONLY, UNREDUCED, check_linearpart, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  addCov(linearpart, NULL, NULL);
  RandomShape(struct_linearpart); 
 
  PREDICT_CALL =
    IncludeModel("predict", InterfaceType, 1, 1, 1, NULL,
	       XONLY, UNREDUCED, check_predict, range_predict, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, MISMATCH);
  kappanames("register", INTSXP);
  addCov(predict, NULL, NULL);
  RandomShape(struct_predict); 
 

  
  // ----------------------------
  // processes  //        CE CO CI TBM Sp di sq Ma av n mpp Hy spf any

  DOLLAR_PROC 
    = IncludeModel("$proc", ProcessType,		   
		   1, 1, 5, kappaS, // kappadollar,
		   XONLY, UNREDUCED, checkS, rangeS, PREF_ALL,
		   true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP,
		   SUBMODEL_DEP);
  // do not change Order!!
  nickname("S");
  kappanames("var", REALSXP, "scale", REALSXP, "anisoT", REALSXP,
	     "Aniso", REALSXP, "proj", INTSXP);
  change_typeof(DAUSER, ShapeType);
  change_typeof(DSCALE, RandomOrShapeType);
  subnames("phi");
  RandomShape(2, structSproc, initSproc, doSproc, true, true, true);
  addSpecific(DOLLAR);

  
  PLUS_PROC = 
    IncludeModel("plusproc", ProcessType, 1, MAXSUB, 0, NULL, 
		  XONLY, UNREDUCED,
		 checkplusproc, NULL, PREF_ALL, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  RandomShape(2, structplusproc, initplusproc, doplusproc, false, false, true);
  addSpecific(PLUS);

 
  IncludeModel("prodproc", ProcessType, 1, 1, 0, NULL, 
	       XONLY, UNREDUCED,
	       checkprodproc, NULL, PREF_ALL, 
	       true, PARAM_DEP, INFDIM-1, false, NOT_MONOTONE);
  RandomShape(2, structprodproc, initprodproc, doprodproc, false, false, true);
  addSpecific(PROD);

  IncludeModel("trafoproc", ProcessType, 1, 1, 1, NULL, 
	       XONLY, UNREDUCED,
	       checktrafoproc, rangetrafo, PREF_ALL, 
	       true, PARAM_DEP, INFDIM-1, false, NOT_MONOTONE);
  kappanames("new", INTSXP);
  change_typeof(TRAFO_ISO, NN2); // ISONAMES
  RandomShape(2, structtrafoproc, inittrafoproc, dotrafoproc,
  	      false, false, true);
  addSpecific(TRAFO);
  

  MPPPLUS_PROC =
  IncludeModel("mppplusproc", ProcessType, 1, MAXSUB, 1, kappamppplus, 
		  XONLY, UNREDUCED,
		 checkmppplus, rangempplus, PREF_ALL, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  nickname("mppplus");
  kappanames("p", REALSXP);
  RandomShape(2, struct_mppplus, init_mppplus, do_mppplus, true, true, true);
  addSpecific(MPPPLUS);


  MULT_PROC = 
  IncludeModel("multproc", ProcessType, 1, MAXSUB, 1, NULL, 
	       XONLY, UNREDUCED,
	       checkmultproc, rangemultproc, PREF_ALL, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP, SUBMODEL_DEP);
  nickname("mult");
  kappanames("multicopies", INTSXP);
  RandomShape(2, structmultproc, initmultproc, domultproc, false, false, true);
  addSpecific(MULT);

 

  AVERAGE_USER = 
    IncludeModel(METHODNAMES[Average], GaussMethodType, 1, 2, 3, kappaGProc, 
		 XONLY, UNREDUCED,
		 check_randomcoin, range_randomcoin, PREF_NOTHING,  
		 false, SCALAR, MAXMPPDIM, false, MISMATCH);
  subnames("phi", "shape");
  kappanames("boxcox", REALSXP, "intensity", REALSXP, "method", INTSXP);  
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
 // addCov(coin, NULL, NULL, coinInverse);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 

  AVERAGE_INTERN = 
    CopyModel("averageIntern", AVERAGE_USER);
  nickname(METHODNAMES[Average]);
  make_internal(); 
  RandomShape(2, struct_randomcoin, init_randomcoin, dompp, true, true, false);


  CIRCEMBED = // und die anderen fehlen auch noch !!
    IncludeModel(METHODNAMES[CircEmbed], GaussMethodType, 1, 1, 13, kappa_ce,
		 XONLY, UNREDUCED,
		 check_ce, range_ce, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXCEDIM, false, MISMATCH);
  kappanames("boxcox", REALSXP,
	     CE[CE_FORCE - COMMON_GAUSS - 1], INTSXP,
	     CE[CE_MMIN - COMMON_GAUSS - 1], REALSXP, 
	     CE[CE_STRATEGY - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_MAXGB - COMMON_GAUSS - 1], REALSXP,
	     CE[CE_MAXMEM - COMMON_GAUSS - 1], INTSXP,
	     CE[CE_TOLIM - COMMON_GAUSS - 1], REALSXP, 
	     CE[CE_TOLRE - COMMON_GAUSS - 1], REALSXP,
	     CE[CE_TRIALS - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_USEPRIMES - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_DEPENDENT - COMMON_GAUSS - 1], INTSXP,
	     CE[CE_APPROXSTEP - COMMON_GAUSS - 1], REALSXP, 
	     CE[CE_APPROXMAXGRID - COMMON_GAUSS - 1], INTSXP);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_ce_approx, init_ce_approx, do_ce_approx);
 
  CE_CUTOFFPROC_USER =  
    IncludeModel(METHODNAMES[CircEmbedCutoff], GaussMethodType, 1, 1, 15, 
		 kappa_localproc, XONLY, UNREDUCED,
		 check_local_proc, range_co_proc, PREF_NOTHING,
		 false, SCALAR, MAXCEDIM,  false, MISMATCH);
  kappanames("boxcox", REALSXP,
	     CE[CE_FORCE - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_MMIN - COMMON_GAUSS - 1], REALSXP, 
	     CE[CE_STRATEGY - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_MAXGB - COMMON_GAUSS - 1], REALSXP,
	     CE[CE_MAXMEM - COMMON_GAUSS - 1], INTSXP,
	     CE[CE_TOLIM - COMMON_GAUSS - 1], REALSXP,
	     CE[CE_TOLRE - COMMON_GAUSS - 1], REALSXP,
	     CE[CE_TRIALS - COMMON_GAUSS - 1], INTSXP,
	     CE[CE_USEPRIMES - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_DEPENDENT - COMMON_GAUSS - 1], INTSXP,
	     CE[CE_APPROXSTEP - COMMON_GAUSS - 1], REALSXP, 
	     CE[CE_APPROXMAXGRID - COMMON_GAUSS - 1], INTSXP,
 	     "diameter", REALSXP, "a", REALSXP);  
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 
 
  CE_CUTOFFPROC_INTERN = CopyModel("cutoffIntern", CE_CUTOFFPROC_USER);
  make_internal();
  RandomShape(2, struct_ce_approx, init_ce_approx, do_ce_approx);
  nickname(METHODNAMES[CircEmbedCutoff]);


  CE_INTRINPROC_USER =  
    IncludeModel(METHODNAMES[CircEmbedIntrinsic], GaussMethodType, 
		 1, 1, 15, kappa_localproc, 
		 XONLY, UNREDUCED,
		 check_local_proc, range_intrinCE, PREF_NOTHING,
		 false, SCALAR, MAXCEDIM, false, MISMATCH);
 nickname(METHODNAMES[CircEmbedIntrinsic]);
  kappanames("boxcox", REALSXP,
	     CE[CE_FORCE - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_MMIN - COMMON_GAUSS - 1], REALSXP, 
	     CE[CE_STRATEGY - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_MAXGB - COMMON_GAUSS - 1], REALSXP,
	     CE[CE_MAXMEM - COMMON_GAUSS - 1], INTSXP,
	     CE[CE_TOLIM - COMMON_GAUSS - 1], REALSXP, 
	     CE[CE_TOLRE - COMMON_GAUSS - 1], REALSXP,
	     CE[CE_TRIALS - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_USEPRIMES - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_DEPENDENT - COMMON_GAUSS - 1], INTSXP, 
	     CE[CE_APPROXSTEP - COMMON_GAUSS - 1], REALSXP, 
	     CE[CE_APPROXMAXGRID - COMMON_GAUSS - 1], INTSXP,
	     "diameter",REALSXP, "rawR", REALSXP);  
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
   RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 
 
  CE_INTRINPROC_INTERN = CopyModel("intrinsIntern", CE_INTRINPROC_USER);
  make_internal();
  RandomShape(2, struct_ce_approx, init_ce_approx, do_ce_approx);

  DIRECT = 
    IncludeModel(METHODNAMES[Direct], GaussMethodType, 
		 1, 1, 2, kappaGProc, XONLY, UNREDUCED,
		 check_directGauss, range_direct, PREF_NOTHING,
		 false,  SUBMODEL_DEP, INFDIM-1, false, MISMATCH);
  kappanames("boxcox", REALSXP,
	     direct[DIRECT_MAXVAR - COMMON_GAUSS - 1], INTSXP);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
 RandomShape(2, init_directGauss, do_directGauss);


  HYPERPLANE_USER  = 
    IncludeModel(METHODNAMES[Hyperplane], GaussMethodType, 1, 1, 6, kappaGProc,
		 XONLY, UNREDUCED,
		 check_hyperplane, range_hyperplane, PREF_NOTHING,  
		 false, SCALAR, 2, false, MISMATCH);
  kappanames("boxcox", REALSXP,
	     "superpos", INTSXP, "maxlines", INTSXP, "mar_distr", INTSXP, 
	     "mar_param", REALSXP, "additive", INTSXP);
  //  addCov(IdStat, NULL, NULL, IdInverse);
  //  addCov(IdNonStat);
 change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 

  HYPERPLANE_INTERN =   
    CopyModel("hyperIntern", HYPERPLANE_USER, check_hyperplane_intern);
  make_internal();
  nickname(METHODNAMES[ Hyperplane]);
 RandomShape(2, struct_hyperplane, init_hyperplane, do_hyperplane);
  //  printf("%ld %ld\n", CovList[HYPERPLANE_INTERN].nonstat_cov,
  //	 IdNonStat); assert(false);

  NUGGET_USER  = 
    IncludeModel(METHODNAMES[Nugget], GaussMethodType, 
		 1, 1, 3, kappaGProc,XONLY, UNREDUCED,
		 check_nugget_proc, range_nugget_proc, PREF_NOTHING, false, 
		 PREVMODEL_DEP, INFDIM, true, MISMATCH);
  kappanames("boxcox", REALSXP,"tol", REALSXP, "vdim", INTSXP);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 


  NUGGET_INTERN =   
    CopyModel("nuggetIntern", NUGGET_USER);
  make_internal();
  RandomShape(2, struct_nugget, init_nugget, do_nugget);
  nickname(METHODNAMES[Nugget]);
  /* see simu.cc, CMbuild for special treatment of nugget when
     users choice is given */
  
  /* cf. convert.R, PrepareModel, near end of function */
  /* deleted 27.12.08
     { 
     int nr = currentNrCov - 1;
     for (i=0; i<MethodsLength; i++) {
     CovList[nr].implemented[i] = GIVEN_METH_IGNORED;
     }
     CovList[nr].implemented[Nugget] = 
     CovList[nr].implemented[Direct] =
     CovList[nr].implemented[Sequential] =
     CovList[nr].implemented[CircEmbed] =IMPLEMENTED;
     }
  */
 
  RANDOMCOIN_USER = CopyModel(METHODNAMES[RandomCoin], AVERAGE_USER);

 
  SEQUENTIAL = 
    IncludeModel(METHODNAMES[Sequential], GaussMethodType, 1, 1, 4, kappaGProc, 
		 XONLY, UNREDUCED,
		 check_sequential, range_sequential, PREF_NOTHING,  
		 false, SCALAR, INFDIM-1, false, MISMATCH);
  kappanames("boxcox", REALSXP,"max_variables", INTSXP, "back_steps", INTSXP, "initial", INTSXP);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, init_sequential, do_sequential);


  SPECTRAL_PROC_USER = 
    IncludeModel(METHODNAMES[SpectralTBM], GaussMethodType,  1, 1, 5,kappaGProc,
		 XONLY, UNREDUCED,
		 check_spectral, range_spectral, PREF_NOTHING,  
		 false, SCALAR, MAXTBMSPDIM, false, MISMATCH);
  kappanames("boxcox", REALSXP,"sp_lines", INTSXP, "sp_grid", INTSXP,
	     "prop_factor", REALSXP, "sigma", REALSXP );
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 
 
  SPECTRAL_PROC_INTERN = 
    CopyModel("spectralIntern", SPECTRAL_PROC_USER);
  make_internal();
  RandomShape(2, struct_spectral, init_spectral, do_spectral);
  nickname(METHODNAMES[SpectralTBM]);


  SPECIFIC = 
    IncludeModel(METHODNAMES[Specific], GaussMethodType, 1, 1, 1, kappaGProc,
		 XONLY, UNREDUCED,
		 check_specificGauss, range_specificGauss, PREF_NOTHING,  
		 false, SUBMODEL_DEP, MAXTBMSPDIM, false, MISMATCH);
  RandomShape(2, struct_specificGauss, init_specificGauss, do_specificGauss);
  kappanames("boxcox", REALSXP);
 change_sortof(GAUSS_BOXCOX, ANYPARAM);
 
  TBM_PROC_USER = 
    IncludeModel(METHODNAMES[TBM], GaussMethodType,  1, 1, 9, tbm_kappasproc, 
		 XONLY, UNREDUCED, 
		 checktbmproc, rangetbmproc, PREF_NOTHING,
		 false, PARAM_DEP, SUBMODEL_DEP, false, MISMATCH);
  kappanames("boxcox", REALSXP, 
	     "fulldim", INTSXP, "reduceddim", INTSXP, "layers", REALSXP,
	     "lines", INTSXP, "linessimufactor", REALSXP,
	     "linesimustep",  REALSXP,  // "grid", INTSXP,
	     "center", REALSXP, "points", INTSXP); 
  // addFurtherCov(tbm2num, NULL);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_extractdollar, init_gaussprocess, do_gaussprocess); 
 
  TBM_PROC_INTERN = 
    CopyModel("tbmIntern", TBM_PROC_USER);
  make_internal();
  RandomShape(2, struct_tbmproc, init_tbmproc, do_tbmproc); 
  nickname(METHODNAMES[TBM]);
 
  TRENDEVAL = 
    IncludeModel(METHODNAMES[TrendEval], GaussMethodType, 
		 1, 1, 1, kappaGProc, XONLY, UNREDUCED,
		 checkTrendEval, range_TrendEval, PREF_TREND,
		 false,  SUBMODEL_DEP, INFDIM-1, false, MISMATCH);
  RandomShape(2, init_TrendEval, do_TrendEval);
  kappanames("boxcox", REALSXP);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);


  gaussmethod[CircEmbed] = CIRCEMBED;
  gaussmethod[CircEmbedCutoff]= CE_CUTOFFPROC_INTERN;
  gaussmethod[CircEmbedIntrinsic] = CE_INTRINPROC_INTERN;
  gaussmethod[TBM] = TBM_PROC_INTERN;
  gaussmethod[SpectralTBM] = SPECTRAL_PROC_INTERN;
  gaussmethod[Direct] = DIRECT;
  gaussmethod[Sequential] = SEQUENTIAL;
  gaussmethod[TrendEval] = TRENDEVAL;
  gaussmethod[Average] = AVERAGE_INTERN;
  gaussmethod[Nugget] = NUGGET_INTERN;
  gaussmethod[RandomCoin] = AVERAGE_INTERN;
  gaussmethod[Hyperplane] = HYPERPLANE_INTERN;
  gaussmethod[Specific] = SPECIFIC;
  gaussmethod[Nothing] =  gaussmethod[Forbidden] = -1;
 

  // non sub-gaussian processe
  BRORIGINAL_USER =
    IncludeModel("brorig", BrMethodType, 1, 2, 3, NULL, 
		 XONLY, UNREDUCED,
		 checkBrownResnickProc, range_mpp, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP);
  subnames("phi", "tcf");
  RandomShape(0, structBRuser, initBRuser, dompp); 
  addlogD(loglikelihoodBR);

  BRMIXED_USER =
    IncludeModel("brmixed", BrMethodType, 1, 2, 11, kappaBRmixed, 
		 XONLY, UNREDUCED, 
		 check_BRmixed, range_BRmixed, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP, 
	     "meshsize", REALSXP, "vertnumber", INTSXP,
             "optim_mixed", INTSXP, "optim_mixed_tol", REALSXP, 
	     "optim_mixed_maxpoints", INTSXP, "lambda", REALSXP,
	     "areamat", REALSXP, "variobound", REALSXP);
  subnames("phi", "tcf");
    RandomShape(0, structBRuser, initBRuser, dompp);
  addlogD(loglikelihoodBR);
  
  BRSHIFTED_USER = 
    IncludeModel("brshifted", BrMethodType, 
		 1, 2, 3, NULL, 
		 XONLY, UNREDUCED, 
		 checkBrownResnickProc, range_mpp, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
 subnames("phi", "tcf");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP);
   RandomShape(0, structBRuser, initBRuser, dompp);
  addlogD(loglikelihoodBR);
  
  BRORIGINAL_INTERN =
    CopyModel("brorigIntern", BRORIGINAL_USER, PointShapeType);
  make_internal();
  nickname("brorig");
  RandomShape(SUBMODEL_DEP, structBRintern, init_BRorig, do_BRorig);
  
  BRMIXED_INTERN =
    CopyModel("brmixedIntern", BRMIXED_USER, PointShapeType); 
  make_internal();
  nickname("brmixed");
  RandomShape(SUBMODEL_DEP, structBRintern, init_BRmixed, do_BRmixed);
  
  BRSHIFTED_INTERN =
    CopyModel("brshiftIntern", BRSHIFTED_USER, PointShapeType); 
  make_internal();
   nickname("brshif");
 RandomShape(SUBMODEL_DEP, structBRintern, init_BRshifted, do_BRshifted);
  
  // distributions

  
 BINARYPROC = // direct an Gauss gekoppelt!!
    IncludeModel("binaryprocess", ProcessType, 1, 1, 3, kappa_binaryprocess, 
		 XONLY, UNREDUCED, 
		 checkbinaryprocess, rangebinaryprocess, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXSIMUDIM, false, MISMATCH);
  nickname( "bernoulli");
  kappanames(INTERNAL_PARAM, REALSXP, "stationary_only", 
	     INTSXP, "threshold", REALSXP);
  //"p", REALSXP);
  RandomShape(INFTY, struct_binaryprocess, init_binaryprocess,
	      do_binaryprocess);
  
  BROWNRESNICKPROC =
    IncludeModel("brownresnick",ProcessType,  1, 2, 3, NULL, 
		 XONLY, UNREDUCED, 
		 checkBrownResnickProc, range_mpp, PREF_NOTHING,  
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
  subnames("phi", "tcf");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP);
  //  addCov(BrownResnick, NULL, NULL);
  RandomShape(0, structBrownResnick, initBrownResnick, doBrownResnick); 
  addlogD(loglikelihoodBR);

  GAUSSPROC = 
    IncludeModel("gauss.process", ProcessType, 1, 1, 2, kappaGProc,
		 XONLY, UNREDUCED,
		 checkgaussprocess, rangegaussprocess, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXSIMUDIM, false, MISMATCH);
  nickname("gauss");
  kappanames("boxcox", REALSXP, "stationary_only", INTSXP);  
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_gaussprocess, init_gaussprocess, do_gaussprocess);
  addlogD(gaussprocessDlog);
  


  POISSONPROC =
    IncludeModel("poisson", ProcessType, 1, 1, 1, NULL,
		 XONLY, UNREDUCED,  
		 check_poisson, range_poisson, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, false, MISMATCH);
  kappanames("intensity", REALSXP);
   RandomShape(0, struct_poisson, init_poisson, dompp);


  SCHLATHERPROC =
    IncludeModel("extremalgauss", ProcessType, 1, 2, 3, NULL, 
		 XONLY, UNREDUCED,  
		 check_schlather, range_mpp, PREF_NOTHING, 
		 false, SCALAR, MAXSIMUDIM, false, MISMATCH);
  nickname("schlather");
  subnames("phi", "tcf");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP);
  addCov(extremalgaussian, NULL, NULL);
  RandomShape(0, struct_schlather, init_mpp, dompp);
  addlogD(loglikelihoodSchlather);

  

  EXTREMALTPROC =
    IncludeModel("extremalt", ProcessType, 1, 1, 4, NULL, 
		 XONLY, UNREDUCED,  
		 check_schlather, range_opitz, PREF_NOTHING, 
		 false, SCALAR, MAXSIMUDIM, false, MISMATCH);
  nickname("opitz");
  subnames("phi");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP, "alpha", REALSXP);
  addCov(extremalgaussian, NULL, NULL);
  RandomShape(0, struct_schlather, init_opitzprocess, dompp);
  addlogD(loglikelihoodSchlather);

 
  SMITHPROC =
    IncludeModel("smith", ProcessType, 1, 2, 3, NULL, XONLY, UNREDUCED,
		 check_smith, range_mpp, PREF_NOTHING,  
		 false, SCALAR, MAXMPPDIM, false, MISMATCH);
  subnames("shape", "tcf");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP);
  RandomShape(0, struct_smith, init_mpp, dompp);
 
  CHI2PROC =
    IncludeModel("chi2", ProcessType, 1, 1, 2, kappaGProc, XONLY, UNREDUCED,
		 checkchisqprocess, rangechisqprocess, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXSIMUDIM, false, MISMATCH);
  kappanames("boxcox", REALSXP, "f", INTSXP);  
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(0, struct_chisqprocess, init_chisqprocess, do_chisqprocess);

  TPROC =
    IncludeModel("t", ProcessType, 1, 1, 2, kappaGProc, XONLY, UNREDUCED,
		 checkchisqprocess, rangetprocess, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXSIMUDIM, false, MISMATCH);
  kappanames("boxcox", REALSXP, "nu", REALSXP); 
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(0, struct_chisqprocess, init_chisqprocess, do_tprocess);


  // assert(SimpleChecks());  // Zzt auskommentiert
  assert(CheckListcov());
  for (int nr=0; nr<currentNrCov; nr++) {     
    cov_fct *C = CovList + nr; 
    if (C->Typi[0] == RandomType) {    
       // printf("done = %s %d\n", CovList[nr].nick, nr);
      continue;
    }
     if (nr != NUGGET){ 	
      for (int ii=0; ii<Nothing; ii++) {
	C->pref[ii] *= C->implemented[ii] == IMPLEMENTED;
      }
      C->pref[Nothing] *= (C->cov != ErrCov || C->nonstat_cov != ErrCovNonstat);
     }
  }

  IncludeModel("minus", MathDefinition, 0, 0, 3, NULL, XONLY, PREVMODELI,
	       checkMath, rangeMath, PREF_TREND,
	       false,SCALAR, PREVMODEL_DEP, false, false);
  kappanames("a", REALSXP, "b", REALSXP, "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathminus, NULL, NULL);
  AddVariant(TrendType, PREVMODELI);
 
  IncludeModel("plus", MathDefinition, 0, 0, 3, NULL, XONLY, PREVMODELI,
	       checkMath, rangeMath, PREF_TREND, 
	      false,SCALAR, 1, false, false);
  kappanames("a", REALSXP, "b", REALSXP, "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathplus, NULL, NULL);
  AddVariant(TrendType, PREVMODELI);

  IncludeModel("div", MathDefinition, 0, 0, 3, NULL, XONLY, PREVMODELI,
	       checkMath, rangeMath, PREF_TREND, 
	      false,SCALAR, 1, false, false);
  kappanames("a", REALSXP, "b", REALSXP,  "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathdiv, NULL, NULL);
  AddVariant(TrendType, PREVMODELI);

  IncludeModel("mult", MathDefinition, 0, 0, 3, NULL, XONLY, PREVMODELI,
	       checkMath, rangeMath, PREF_TREND, 
	      false,SCALAR, 1, false, false);
  kappanames("a", REALSXP, "b", REALSXP,  "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathmult, NULL, NULL);
  AddVariant(TrendType, PREVMODELI);

  CONST =
  IncludePrim("const", MathDefinition, 1, NULL, XONLY, PREVMODELI,
	      check_c, rangec, PREF_TREND, 
	      SCALAR, PREVMODEL_DEP, false, false);
  kappanames(CONST_A_NAME, REALSXP);
  change_sortof(CONST_C, TRENDPARAM);
  addCov(Mathc, NULL, NULL);
  AddVariant(TrendType, PREVMODELI);
  AddVariant(TcfType, PREVMODELI);
 
  IncludeModel("p", MathDefinition, 0, 0, 3, NULL, XONLY, PREVMODELI,
	       checkproj, rangeproj, PREF_TREND, 
	      false, SCALAR, INFDIM-1, false, false);
  kappanames("proj",  INTSXP, "new", INTSXP, "factor",  REALSXP);
  change_typeof(PROJ_ISO, NN2);
  change_sortof(PROJ_FACTOR , TRENDPARAM);
  addCov(proj, NULL, NULL);
  AddVariant(TrendType, PREVMODELI);

  BIND = 
  IncludeModel("c", MathDefinition, 0, 0, 11, NULL, XONLY, PREVMODELI,
	       check_bind, rangeMath, PREF_TREND, 
	      false,SCALAR, 1, false, false);
  kappanames("a", REALSXP, "b", REALSXP, "c", REALSXP, 
	     "d", REALSXP, "e", REALSXP, "f", REALSXP,
	     "g", REALSXP, "h", REALSXP, "i", REALSXP,
	     "j", REALSXP,
	     "factor", REALSXP);
  change_sortof(CovList[BIND].kappas - 1, TRENDPARAM);
  addCov(Mathbind, NULL, NULL);
  AddVariant(TrendType, PREVMODELI);

  IncludeModel("is", MathDefinition, 0, 0, 3, NULL, XONLY, PREVMODELI,
	       checkMath, rangeMathIs, PREF_TREND, 
	      false, SCALAR, 1, false, false);
  kappanames("a", REALSXP, "is", INTSXP, "b", REALSXP);
  change_typeof(IS_IS, NN1);
  addCov(Mathis, NULL, NULL);
  AddVariant(TrendType, PREVMODELI);


  includeStandardMath();
  for (i=0; i<currentNrCov; i++) {
    if (CovList[i].Typi[0] == MathDefinition) {
      // printf("%s\n", CovList[i].name);
      //assert(CovList[i].variants <= 2)
      CovList[i].Typi[0] = ShapeType;
    }
  }



}


bool isDollar(cov_model *cov) {
  assert(cov != NULL);
  int nr=cov->nr;
  //  printf("---> %d %d %d %d\n", cov->nr, DOLLAR, LASTDOLLAR, nr >= DOLLAR && nr <= LASTDOLLAR);
  return nr >= DOLLAR && nr <= LASTDOLLAR;
}


bool isDollarProc(cov_model *cov) {
  assert(cov != NULL);
  int nr=cov->nr;
  return nr == DOLLAR_PROC;
}

bool isAnyDollar(cov_model *cov) {
  assert(cov != NULL);
  int nr=cov->nr;
  return (nr >= DOLLAR && nr <= LASTDOLLAR) || nr == DOLLAR_PROC;
}

bool isNatsc(cov_model *cov) {
  assert(cov != NULL);
  int nr=cov->nr;
  return nr == NATSC_INTERN || nr == NATSC_USER;
}


bool isDef(Types type, typusfct iso, cov_model *cov) {
  int i;
  cov_fct *C = CovList + cov->nr;

  if (isUndefined(C)) {
    assert(CovList[cov->nr].TypeFct != NULL);
    return CovList[cov->nr].TypeFct(type, cov, 0);
  }

  for (i=0; i<C->variants; i++) {
    if (iso(C->Typi[i]) && (cov->isoown == ISO_MISMATCH || 
			    atleastSpecialised(cov->isoown, C->Isotropy[i])))
      return true;
  }

  return false;

}

//if (false) printf("Qu: "#typus" '%s'  %d\n", TYPENAMES[type], type == cond);
#define Question(typus, cond)						\
  bool is##typus(Types type) { 					\
    return type == cond;						\
  }									\
  									\
  bool is##typus(cov_model *cov) {return isDef(typus##Type, is##typus, cov);} 
 
Question(Tcf, TcfType)
Question(PosDef, PosDefType || isTcf(type))
Question(Variogram, VariogramType || isPosDef(type))
Question(NegDef, NegDefType || isVariogram(type))
Question(Process, ProcessType || type == GaussMethodType || 
	 type == BrMethodType)
Question(GaussMethod, GaussMethodType)
Question(PointShape, PointShapeType)
Question(Random, RandomType)
Question(Shape, ShapeType || isNegDef(type))
Question(Trend, TrendType)
Question(Interface, InterfaceType)
Question(Undefined, UndefinedType)
Question(Other, OtherType)

bool isUndefined(cov_fct *C) {
  //printf("%s cvariants %d\n", C->name, C->variants);
  bool is = isUndefined(C->Typi[0]);
  assert(!is || C->variants == 1) ;
  // for(bool undef=true; int i=0; i<C->variants; i++) undef &= C->Typi[i] == Undefined;
  return is;
}



////////////

bool isGaussProcess(cov_model *cov) {
  int nr = cov->nr;
  return nr == GAUSSPROC || isGaussMethod(cov);
}

bool isBernoulliProcess(cov_model *cov) {
  int nr = cov->nr;
  return nr == BINARYPROC;
}

bool isGaussBasedProcess(cov_model *cov) {
  int nr = cov->nr;
  return isGaussProcess(cov) || nr == CHI2PROC || nr == TPROC;
}

bool isBrownResnickProcess(cov_model *cov) {
  //PMI(cov);
  Types type = CovList[cov->nr].Typi[0]; 
  int nr = cov->nr;
  bool is = type == BrMethodType || nr == BROWNRESNICKPROC;
  assert(!is || CovList[cov->nr].variants == 1);
  return is;
}

bool isBRuserProcess(cov_model *cov) {
  Types type = CovList[cov->nr].Typi[0]; 
  bool is = type == BrMethodType;
  assert(!is || CovList[cov->nr].variants == 1);
  return is;
 }

bool isBRuserProcess(Types type) {
  return type == BrMethodType;
}

bool isMaxStable(cov_model *cov) {
  int nr = cov->nr;
  return isBrownResnickProcess(cov) || nr == SMITHPROC 
    || nr == SCHLATHERPROC || nr == EXTREMALTPROC;
}



bool isCov(cov_model *cov) {
  int nr = cov->nr;
  return nr == COVFCTN || nr == COVMATRIX ;
}


bool hasNoRole(cov_model *cov) {
  int role = cov->role;
  return role == ROLE_BASE;
}

bool hasExactMaxStableRole(cov_model *cov) {
  int role = cov->role;
  return role == ROLE_MAXSTABLE;
}

bool hasMaxStableRole(cov_model *cov) {
  int role = cov->role;
  return role == ROLE_MAXSTABLE || role == ROLE_BROWNRESNICK ||
    role == ROLE_SMITH || role == ROLE_SCHLATHER;
}

bool hasPoissonRole(cov_model *cov) {
  int role = cov->role;
  return role == ROLE_POISSON_GAUSS || role == ROLE_POISSON;
}


bool hasAnyShapeRole(cov_model *cov) {
  return hasMaxStableRole(cov) || hasPoissonRole(cov) || hasDistrRole(cov);
}

bool hasDistrRole(cov_model *cov) {
  int role = cov->role;
  return role == ROLE_DISTR;
}

int role_of_process(int nr) {
  return 
    (nr == AVERAGE_USER || nr == AVERAGE_INTERN ||
     nr == RANDOMCOIN_USER) ? ROLE_POISSON
    : ((nr >= CIRCEMBED &&  nr <= CE_INTRINPROC_INTERN)
       || nr == DIRECT || nr == NUGGET || nr == NUGGET_INTERN
       || nr == SEQUENTIAL 
       || nr == SPECTRAL_PROC_USER || nr == SPECTRAL_PROC_INTERN
       || nr == TBM_PROC_USER || nr == TBM_PROC_INTERN 
       || nr == GAUSSPROC || nr == TRENDEVAL
       ) ? ROLE_GAUSS
    : nr == HYPERPLANE_USER || nr == HYPERPLANE_INTERN ? ROLE_GAUSS
    : nr == SPECIFIC ? ROLE_GAUSS
    : ( nr == BRSHIFTED_USER || nr == BRMIXED_USER || nr == BRORIGINAL_USER
	|| nr == BROWNRESNICKPROC) ? ROLE_BROWNRESNICK 
    : nr == BINARYPROC ? ROLE_BERNOULLI
    : nr == POISSONPROC ? ROLE_POISSON
    : nr == SCHLATHERPROC ||  nr == EXTREMALTPROC ? ROLE_SCHLATHER
    : nr == SMITHPROC ? ROLE_SMITH
    : ROLE_FAILED;
}


bool isMonotone(cov_model *cov) {
  int monotone = cov->monotone;
  return monotone >= MONOTONE && monotone <= COMPLETELY_MON ;
} 

bool isMonotone(int monotone) {
  return monotone >= MONOTONE && monotone <= COMPLETELY_MON;
}

bool isCompletelyMonotone(int monotone) {
  return monotone == COMPLETELY_MON;
}

bool isNormalMixture(cov_model *cov) {
  int monotone = cov->monotone;
  return monotone == NORMAL_MIXTURE || monotone == COMPLETELY_MON;
}


bool isNormalMixture(int monotone) {
  return monotone == NORMAL_MIXTURE || monotone == COMPLETELY_MON;
}


bool isBernstein(cov_model *cov) {
  int monotone = cov->monotone;
  return monotone == BERNSTEIN;
}

bool isBernstein(int monotone) {
  return monotone == BERNSTEIN;
}

bool isGneiting(cov_model *cov) {
  int monotone = cov->monotone;
  return monotone == GNEITING_MON || monotone == COMPLETELY_MON;
}

bool isGneiting(int monotone) {
  return monotone == GNEITING_MON || monotone == COMPLETELY_MON;
}




bool is_any(isotropy_type iso, cov_fct *C) {
  int i;
  for (i=0; i < C->variants; i++) if (iso == C->Isotropy[i]) return true;
  return false;
}

bool is_all(isotropy_type iso, cov_fct *C) {
  int i;
  for (i=0; i < C->variants; i++) if (iso != C->Isotropy[i]) return false;
  return true;
}

bool TrafoOK(cov_model *cov) {// check other things, too, besides gatter ?
  bool ok = ((cov->gatternr >= FIRST_TRAFO && cov->gatternr <= LASTGATTER &&
	      cov->secondarygatternr == MISMATCH)
	     || 
	     (cov->gatternr > LASTGATTER && cov->gatternr <= LAST_TRAFO &&
	      cov->secondarygatternr >= FIRST_TRAFO && 
	      cov->secondarygatternr <= LASTGATTER
	      ))
    && cov->checked;

  if (!ok) {
    PMI(cov->calling); //
    PRINTF("not: %d <= %d <= %d, %d == %d oder %d < %d <= %d, %d <= %d <= %d UND  %d\n", 
	   FIRST_TRAFO, cov->gatternr, LASTGATTER,
	   cov->secondarygatternr, MISMATCH, 
	   //
	   LASTGATTER, cov->gatternr, LAST_TRAFO,
	   FIRST_TRAFO, cov->secondarygatternr, LASTGATTER,
	   //
	   cov->checked);
    BUG;
  }
  return ok;
}

bool isRObject(int type) {
  if (type == CLOSXP) BUG;
  return type == LANGSXP || type == VECSXP || type == ENVSXP;
}


bool ParamAllowsRandomness(cov_model *cov, int i) {
  cov_fct *C = CovList + cov->nr;
  //printf("%s %d %d\n", NAME(cov), i, C->kappaParamType[i] );
  return C->kappaParamType[i] == RandomType || 
    C->kappaParamType[i] == RandomOrShapeType;
}


bool ParamAllowsShape(cov_model *cov, int i) {
  cov_fct *C = CovList + cov->nr;
  //printf("%s %d %d\n", NAME(cov), i, C->kappaParamType[i] );
  return C->kappaParamType[i] == ShapeType || 
    C->kappaParamType[i] == RandomOrShapeType;
}


bool ParamIsTrend(cov_model *cov, int i) {
  return (SortOf(cov, i, 0, 0) == TRENDPARAM) && isTrend(cov->typus);
}
