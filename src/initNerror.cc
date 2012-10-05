 /* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part and error messages

 Copyright (C) 2001 -- 2011 Martin Schlather

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


#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
//#include <sys/timeb.h>
 
#include <string.h>
#include "RF.h"
#include "Covariance.h"
#include <unistd.h>

init_meth init_method[SimulationTypeLength]= 
{
  init_circ_embed, init_circ_embed_cutoff, init_circ_embed_intr,
  init_turningbands2, init_turningbands3, init_spectral, init_directGauss,
  init_sequential, init_markov, init_ave, init_nugget, init_mpp,
  init_hyperplane, init_nothing, init_nothing
};

do_meth do_method[SimulationTypeLength]= {
  do_circ_embed,  do_circ_embed_cutoff, do_circ_embed_intr, 
  do_tbm2, do_tbm3, do_spectral, do_directGauss, 
  do_sequential, do_markov, do_ave, do_nugget, do_addmpp,
  do_hyperplane, do_nothing, do_nothing
};

key_type KEY[MAXKEYS]; //static key_type KEY[MAXKEYS];
double  ZERO[MAXSIMUDIM], 
//    *userdefinedCovMatrix[MAXDEFMATRIX][MAXMAKEEXPLICITE],
    *OutX=NULL, 
    *OutY=NULL;
int GENERALISEDCAUCHY, STABLE, WHITTLE, BROWNIAN, CAUCHY,
    EXPONENTIAL, MATERN, GAUSS, NUGGET, PLUS, MLEPLUS, TBM2NR,
    MIXEDEFFECT, MLEMIXEDEFFECT, MIXX, VECTOR,
    ELEMENTNR_PLUS = -1, LIST_ELEMENT = -1,
    ISO2ISO, SP2SP, SP2ISO, S2ISO, S2SP, S2S, LASTGATTER, MLE_ENDCOV,
// userdefinedCM_RC[MAXDEFMATRIX][MAXMAKEEXPLICITE], 
    CovMatrixRow, CovMatrixCol, CovMatrixTotal, SPLIT, 
    CovMatrixIdx = 0,
    CumIdxMakeExplicite[MAXMAKEEXPLICITE];
bool CovMatrixIndexActive=false,
    NAOK_RANGE=false;
char CovNames[MAXNRCOVFCTS][MAXCHAR];
cov_model *STORED_MODEL[MODEL_MAX];
char MSG[1000], NEWMSG[1000];

cov_fct *CovList=NULL;
int currentNrCov=-1,
    CONSTANT = -1,
    OUT = -1,
    DOLLAR = -1,
    TREND = -1,
    LASTDOLLAR = -1,
    GATTER = -1,
    NATSC = -1;

int True=1; // never change
int False=0; // never change

#define MAX_CE_MEM 16777216
#define PRINTLEVEL 1
#define NAT_SCALE 0
int PL=PRINTLEVEL, 
    NS=NAT_SCALE;


globalparam GLOBAL = {
    {false, '*', false, PRINTLEVEL, NAT_SCALE, 0 //, false /* aniso */
  },//general_param general;
  {DECISION_CASESPEC, DECISION_CASESPEC},    // decision_param ;
  {false, true, false, TRIVIALSTRATEGY, 3, 0, MAX_CE_MEM,
   -1e-7, 1e-3, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //ce_param ce, 13 NULLEN
  {false, true, false, TRIVIALSTRATEGY, 1, 0, MAX_CE_MEM,
   -1e-9, 1e-7, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, // localce; 13 NULLEN
  {false, true, false, TRIVIALSTRATEGY, 3, 0, 10000000,
   -1e-7, 1e-3, {0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, // tbmce, 13 NULLEN
  { 60, 1, 2.0, 0.0}, //TBM2, specific
  {500, 1, 2.0, 0.0}, //TBM3, specific
  {Nothing, {NA_REAL, NA_REAL, NA_REAL, NA_REAL}, 0},  // TBM3 general
  {true, false, false, 0, 0.0, {500, 500, 500, 500}}, //spectral_param
  {Cholesky, 800, 8192, 1e-12},// direct_param direct;
  {5000, 5, -10}, //sequ_param sequ;
  {2, false, 250000, 1.0}, //markov_param markov;
  {}, // ave
  {0.0, true}, // nugget_param nugget;
  {MPP_AUTO, // --- old: 2.5, 0.1, false, false, //r, dist, random, show
   {100.0, 100.0, 100.0, 100.0}, // intens
   {2.0, 2.0, 2.0, 2.0},  // plus
   {2.5, 2.5, 2.5, 2.5},  // relradius
//   {2.0, 2.0, 2.0, 2.0},  // scale,x
   0.001 /* approxzero */, 0.01, 5 /* sampling dist, r*/,
   0.2/* p */, 1.3333333333333/* beta */
  },  //mpp_param mpp;
  {300, 1000, HYPER_UNIFORM, RF_NAN}, // hyper_param hyper;
  {}, //special_param special;
  {3.0}//extremes_storage extremes;
};



pref_type PREF_ALL = {PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, 
		      PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST,
		      PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST},
  PREF_NOTHING = {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE,
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_BEST};

bool compatible_primitive[SimulationTypeLength] = 
{true, true, true, true, false, false, true, false, true, false, true, false,
 // CE  CO    CI      TBM23       Sp    di     sq    Ma     av    n     mpp 
 // Hy 
 true};

/*
  record containing all user specified directions for choosing
  a simulation method. Currently only one element is contained.
*/

double GENERAL_PRECISION = 5e-15;
double EIGENVALUE_EPS = 1e-15;

char ERRORSTRING[MAXERRORSTRING], ERRORSTRING_OK[MAXERRORSTRING], 
  ERRORSTRING_WRONG[MAXERRORSTRING],
    ERROR_LOC[nErrorLoc]="";
int ERRORMODELNUMBER = -1;

char METHODNAMES[][METHODMAXCHAR]={"circulant embed", //0
				   "cutoff CE",
				   "intrinsic CE",
                                   "TBM2", 
				   "TBM3",
				   "spectral TBM", //5
                                   "direct decomp.",
				   "sequential",
				   "Markov",
				   "average",
				   "nugget", //10
                                   "coins",
                                   "hyperplanes",
				   "any method", // nothing
				   "max.MPP",
				   "extremalGauss",
				   "forbidden"};
char ISONAMES[LAST_ISO + 1][15] =  { 
    "isotropic", "spaceisotropic", "zerospaceiso", "anisotropic", "unset", 
    "", "", "", "", "mismatch"};
char STATNAMES[LAST_STAT + 1][15] =  {
    "stationary", "variogram", "IRF",
    "aux. matrix", "aux. vector", "gen. cov.",
    "covariance", "gen. variogram", "unset", "mismatch"};
char STANDARDPARAM[MAXPARAM][MAXCHAR] = {"k1", "k2", "k3", "k4", "k5", "k6"};
char STANDARDSUB[MAXSUB][MAXCHAR] = {"u1", "u2", "u3", "u4", "u5", "u6", "u7", 
			      "u8", "u9", "u10"};
char DISTRNAMES[DISTRNR][DISTRMAXCHAR]={"Gauss", "Poisson", "MaxStable"};
//                                          ONLY used in Hyperfct.cc
//                                          (cutoff embedding and intrinsic)

void errorMSG(int err, char* m) {
  int PL = GLOBAL.general.printlevel;
  if (PL > 4) {
    if (PL > 6) {
      PRINTF("\n\n================================ (PL=%d)\n",  PL);
      PrintModelInfo(STORED_MODEL[MODEL_USER]);   // OK
      PrintModelInfo(STORED_MODEL[MODEL_INTERN]); // OK
      PrintModelInfo(STORED_MODEL[MODEL_SIMU]);   // OK
    }
    PRINTF("error code %d\n", err);
  }
  
  if (err >= ERRORM && err <= ERRORMEND) err = ERRORM;

  switch (err) {
      case ERRORCURRENTLY : strcpy(m,"currently functionality is not available"); break;
      case NOERROR : strcpy(m,"none"); break;
      case NOERROR_REPEAT : strcpy(m,"none; looking for further covariances applicable to the same method");break;
      case NOERROR_ENDOFLIST : strcpy(m,"none; end of list");break;
      case ERRORDUMMY : strcpy(m,"none (dummy)"); break;
       

      case ERRORNOTDEFINED :       
	strcpy(m,"method undefined for this model");break;
      case ERRORNOTPROGRAMMED :    
	strcpy(m,"not programmed yet. Sorry"); break;
      case ERRORCOVNOTALLOWED :    
	strcpy(m,"model not allowed for specified dimension");break;
      case ERRORFAILED: 
	strcpy(m,"algorithm failed (partially)");break;
      case ERRORMEMORYALLOCATION: 
	strcpy(m,"memory allocation error"); break;
      case ERRORNOTINITIALIZED: 
	strcpy(m,"not initialized or RFparameters()$Storing==FALSE");break;
      case ERRORKEYNROUTOFRANGE: 
	strcpy(m,"wrong key number");break;
      case ERRORDECOMPOSITION:
	strcpy(m,"matrix decomposition failed");break;
      case ERRORPRECISION: 
	strcpy(m,
	      "required precision not attained: probably invalid model.\nSee also RFparameters()$direct.svdtolerance.");
	break;
      case ERRORFOURIER: 
	strcpy(m,"fft factorization failed");break;
      case ERRORCOVFAILED: 
	sprintf(m, "model and method only valid for %s. Got %s",
		ERRORSTRING_OK,ERRORSTRING_WRONG);
	break;
      case ERRORNOMULTIVARIATE :
	strcpy(m, "multivariate models not allowed");
      case ERRORMULTIMISMATCH: 
	sprintf(m, "multivariate dimension of %s expected. Got %s",
		ERRORSTRING_OK, ERRORSTRING_WRONG);
	break;
      case ERRORANYFCT: 
	strcpy(m, "auxiliary vector expected");
	break;
      case ERRORCOVNUMERICAL: 
	sprintf(m,
		"covariance exactly calculable in TBM2 only for %s. Got %s. Try TBM2.num=TRUE",
		ERRORSTRING_OK,ERRORSTRING_WRONG);
	break;
      case ERRORCEMAXMEMORY:
	sprintf(m,
		"total real numbers needed=%s > allowed max real numbers=%s -- increase CE.maxmem and/or local.maxmem using RFparameters",
		ERRORSTRING_WRONG, ERRORSTRING_OK);
	break;
      case ERRORTOOMANYLINES:
	strcpy(m, "estimated number of lines exceeds hyper.maxlines");break;
      case ERRORREGISTER: 
	sprintf(m,"Register number out of [0,%d]",MAXKEYS-1);break;
      case ERRORCOORDINATES: 
	strcpy(m,"coordinates are not given or invalid grid specification");
	break;
      case ERRORNEGATIVEVAR: 
	strcpy(m,"Variance or nugget non-positive or not symmetric");break;
      case ERRORDIM: 
	sprintf(m,"dimension specification not in [1,%d] or xdim > tsdim",
		MAXSIMUDIM);break;
      case ERRORNEGATIVESCALE :  
	strcpy(m,"scale parameter must be positive");break;
      case ERRORWAVING :
	strcpy(m,"Rescaling not possible (waving)");break;
      case ERRORRESCALING:
	strcpy(m,"practical range not defined");
	break;
    case ERRORWRONGINIT: 
	strcpy(m,"Wrong initialisation");break;
      case ERRORNN:
	strcpy(m,"The number of points on the line is too large. Check your parameters and make sure that none of the locations are given twice");break;
      case ERRORANISOTROPIC:
	strcpy(m,"anisotropic function not allowed"); break;
      case ERRORNOSTATMATCH : 
	strcpy(m,"no matching assumption found for stationarity and isotropy");
	break;
      case ERRORTIMENOTANISO:
	strcpy(m,"time component only allowed for anisotropic fields"); break; 
      case ERRORANISO:
	strcpy(m,"anisotropic call not allowed"); break; 
     case ERRORCIRCNONPOS:
	strcpy(m,"embedded matrix has (repeatedly) negative eigenvalues and approximation\nis not allowed (*CE.force=FALSE)"); break; 
      case ERRORANISOMIX:
	strcpy(m,"for multiplicative covariance functions and TBM \
the anisotropy matrices must be identical"); break;
      case ERRORISOTROPICMETHOD:
	strcpy(m,"method (TBM) only allows for spatially isotropic functions in certain dimensions"); 
	break;
      case ERRORTIMESEPARATE:
	strcpy(m,"method (TBM) does not allow the time component being mixed up with spatial components"); 
	break;
      case ERRORTIMECOMPLETESEPARATE:
	strcpy(m,"method does not allow the time component being mixed up with spatial components or vice versa"); 
	break;
      case ERRORUNKNOWNMETHOD:
	strcpy(m,"Unknown method in this context or unallowed mixture of methods"); 
	break;
      case ERRORWITHOUTTIME:
	strcpy(m,"spatially isotropic covariance fcts that are not fully isotropic must be called with a genuine time component T"); 
	break;
      case ERRORCOVNROUTOFRANGE: 
	strcpy(m,"wrong covnr number");break;
      case ERRORWRONGDIM:
	strcpy(m,"wrong dimension"); break;
      case ERRORNOTANISO:
	strcpy(m,"model definition must use list-notation with anisotropy"); break;
      case ERRORDIMMISMATCH:
	strcpy(m,"logical dimension does not match physical dimension"); break;
      case ERRORTOOMANYPOINTS:
	strcpy(m, "too many points to determine smallest distance between the points; try TBM*.linesimustep with a positive value"); break;
      case ERRORTIMENOTALLOWED:
	strcpy(m, "time component not allowed for the specified method"); break;
      case ERRORANYLEFT:
	strcpy(m, "not all covariances could be simulated at the current step");
	break;
      case ERRORMETHODEXCLUDED:
	strcpy(m, "CE Intrinsic contradictsRFparmeters()$stationary.only=true");
	break;
      case ERROROUTOFMETHODLIST:
	sprintf(m, 
         "run out of list of methods --\n  %s%s",
		GLOBAL.general.skipchecks
		? "did you try an invalid parameter combination?"
		: "are RFparameters too restrictive?",
		PL < 4 ? 
	 "\n  Retry with RFparameters(PrintLevel=4) to get more information"
		: ""
	  );
	break;
      case ERRORTBMCENTER:
	  strcpy(m, "center contains non finite values");
	 break;  
      case ERRORTBMPOINTS:
	  strcpy(m,
		 "given number of points to simulate on the TBM line is too small");
	 break;  
      case ERRORHYPERNR :
	  strcpy(m, 
		 "announced number of submodels invalid (too large or non-positive)");
	  break;
      case ERRORHYPERMETHOD :
	  strcpy(m, 
		 "methods for the hypermodel and the submodels must be given and be equal");
	  break;
      case ERRORHYPERNOTALLOWED :
	  strcpy(m, "the method does not allow for hyper models (yet).");
	  break;
      case ERRORNCOVOUTOFRANGE:
	  strcpy(m, "the number of covariance functions is zero or too large");
	  break;
      case ERRORISO:
	  strcpy(m, "isotropic model not allowed");
	  break;
      case ERRORM: 
	strcpy(m, ERRORSTRING);
	break;
      case ERRORMSG: 
	sprintf(m, "%s. Got %s", ERRORSTRING_OK,ERRORSTRING_WRONG);
	break;
      case ERRORTBMPOINTSZERO: 
	sprintf(m, "if linesimufactor and linesimustep are naught then TBM_POINTS must be at least 4 (better between 100 and 10000)");
	break;

	// case ERRORISOTROPIC:
	//   strcpy(m,"isotropic function not allowed"); break;
	//
	// extremes:
      case ERRORFULLRANK: 
	sprintf(m, "anisotropy matrix must have full rank");
	break;
      case ERRORTIMECOMPONENT:
	sprintf(m, "last column of the anisotropy matrix in the hypermodel must consist of zeros only");
	break;	
      case ERRORSUBMETHODFAILED:
	  sprintf(m, "no good submethods exist");
      case ERRORLOWRANK :
	  strcpy(m,"confused about low rank anisotropy matrix");
	  break;
      case ERRORLOWRANKTBM :
	  strcpy(m,"confused about low rank anisotropy matrix.\nWhen using TBM then try putting the model with anisotropy matrix of\nhighest rank at the beginning");
	  break;
      case ERRORLAYERS :
	  strcpy(m, "value of RFparameters()$TBM*.layers does not match the situation");
	  break;
      case ERRORONLYISOTROPIC :
	  strcpy(m, "only isotropic fields are allowed");
	  break;
      case  ERRORSTATVARIO:
	   strcpy(m, "model can be called only by stationary arguments");
	   break;
      case ERRORMARKOVNOTINCLUDED :
	strcpy(m, "Simiulaton by Markov techniques not possible, since Havard Rue's library has not been included; see RandomFields/src/includeMarkov.h for installation instructions");
	  break;
      case ERRORNOVARIOGRAM:
	strcpy(m, "Variogram model allowed in this context");
	break;
      case ERRORONLYREGULARGRID :
	  strcpy(m, "only grids with identical spacings allowed");
	  break;
      case ERRORMARKOVPARAMETER :
	  sprintf(m,
		  "GMRF method not available for the given parameters (%s)",
		  ERRORSTRING_WRONG);
	  break;
      case ERRORMARKOVMAXMEMORY:
	sprintf(m,
		"total number of points (%s) greater than the allowed number (%s) -- increase markov.maxmem using RFparameters",
		ERRORSTRING_WRONG, ERRORSTRING_OK);
	break;
      case ERRORLASTGRID:
	strcpy(m,"last component must be truely given by a non-trivial grid");
	break;
      case ERRORTRIVIALTIMEDIRECTION:
	strcpy(m,"the grid in the last direction is too small; use method `direct' instead of `sequential'");
	break;
     case ERRORMETROPOLIS:
       strcpy(m,"Metropolis search algorithm for optimal sd failed\n -- check whether the scale of the problem has been chosen appropriately");
	break;
     case ERRORLOGMIX:
	strcpy(m,"log mixture (Average) not defined for given submodel");
	break;
      case ERRORNORMALMIXTURE:
	  strcpy(m, "only normal mixtures as first submodel allowed (Gneiting, 2002)");
	  break;
      case ERRORDIAG:
	strcpy(m, "subsequent anisotropy matrices must be diagonal");
	  break;
      case ERRORUNITVAR:
	  strcpy(m, "unit variance required");
	  break;
     case ERRORMAXDIMMETH:
	  strcpy(m, "maximal dimension of variables for the method exceeded");
	  break;
     case ERRORPREVDOLLAR:
	  strcpy(m, "method may not initialised by preceding initS");
	  break;
     case ERRORNOGENUINEMETHOD:
	  strcpy(m, "the genuine method could not be detected -- this is a strange problem; please contact the author");
	  break;

      case ERRORHANGING:
	strcpy(m, "method cannot deal with implicit space reduction -- try with RFparameters(aniso=TRUE)");
	  break;
      case  ERRORCOVUNKNOWN:
	strcpy(m, "covariance model could not be identified");
	  break;

      case  ERRORMATRIX_VECTOR:
	strcpy(m, "Z must be vector");
	break;

      case  ERRORMATRIX_Z:
	strcpy(m, "the values of Z must be in {1,...,ncol(m)}");
	break;

      case ERRORSILLNULL : 
	strcpy(m,"Vanishing sill not allowed");
	break;     
      case ERRORPAIRS : 
	strcpy(m,"No pair simulation (currently) allowed.");
	break;
      case ERRORTREND : 
	strcpy(m,"No trend  (currently) allowed");
	break;
	//    case : strcpy(m,"");break;
	//
	// Poisson:
      case ERRORVARMEAN :
	strcpy(m,"Variance differs from mean");
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
	strcpy(m,
	       "wrong sign of 1st derivative of the covariance function at the break point");
	break;
      case MSGLOCAL_SIGNPHISND :
	strcpy(m,
	       "wrong sign of 2nd derivative of the covariance function at the break point");
	break;
      case MSGLOCAL_INITINTRINSIC :
	strcpy(m,"one of a2, b or a0+phi(0) has wrong sign");
	break;
     case ERRORUNSPECIFIED :
	strcpy(m,"(unspecified)");
	break;
     default : 
	PRINTF(" error=%d\n", err); 
//	cov_model *cov; 	CovList[1000].check(cov);
	assert(false);
  }
}

void ErrorMessage(SimulationType meth, int err) {
  char MS[50], m[300+2*MAXERRORSTRING]; 
  switch(meth) {
      case CircEmbed : strcpy(MS,"circulant embedding"); break;
      case CircEmbedCutoff : strcpy(MS,"cutoff circulant embedding"); break;
      case CircEmbedIntrinsic : strcpy(MS,"intrinsic circulant embedding"); 
	break;
      case TBM2 : strcpy(MS,"2-dim. TBM"); break;
      case TBM3: strcpy(MS,"3-dim. TBM"); break;
      case SpectralTBM: strcpy(MS,"spectral TBM"); break;
      case Direct : 
	  strcpy(MS,"direct Gaussian (decomposition of cov. matrix)");
	break;
      case Sequential: strcpy(MS,"sequential"); break;
      case Markov: strcpy(MS,"Markov (GMRF)"); break;
      case Average: strcpy(MS, "random spatial averages"); break;
      case Nugget: strcpy(MS,"nugget effect modelling"); break;
      case RandomCoin: strcpy(MS,"additive MPP (random coins)"); break;
      case Hyperplane : strcpy(MS,"hyperplane tessellation"); break;
//      case Special: strcpy(MS,"special procedure"); break;
      case Nothing: 
	if (err>0) strcpy(MS,"Error");
	else strcpy(MS,"Message"); 
	break;
      case MaxMpp: strcpy(MS,"max. MPP (Boolean functions)"); break;
      case Forbidden: 
	PRINTF("m=%d (forbidden call); err=%d \n", m, err);
//	strcpy(MS,"forbidden function call"); 
	assert(false);
	break;
      default : PRINTF("m=%d \n",m); assert(false);
  }
  errorMSG(err, m);
  if (strcmp(ERROR_LOC,"") != 0) PRINTF("At %s: ", ERROR_LOC);
  if (meth!=Nothing) PRINTF("Method ");
  PRINTF("%s", MS);
  if (ERRORMODELNUMBER >= 0)
    PRINTF(" in list element #%d", ERRORMODELNUMBER + 1);
  PRINTF(": %s.\n", m);
  ERRORMODELNUMBER = -1;
//  assert(false);
}

int checkOK(cov_model *cov){
  return NOERROR;
}

void EinheitsMatrix(double **mem, int dim) {
     // Einheitsmatrizen
    int d;
    *mem = (double*) calloc(dim * dim, sizeof(double));
    if (*mem == NULL) error("memory allocation error in EinheitsMatrix");
    for (d=0; d<dim; d+=dim+1) *mem[d] = 1.0;
}


void InitModelList() {

  assert(currentNrCov=-1);
  assert(CUTOFF_THEOR == 4);/* do not change this value as used in RFmethods.Rd */

  int i;

  for (i=0; i<MAXSIMUDIM; i++) ZERO[i] = 0.0;

//  for (k=0; k<MAXMAKEEXPLICITE; k++) {
 //     CumIdxMakeExplicite[k] = 0;
 //     for (i=0; i<MAXDEFMATRIX; i++) userdefinedCovMatrix[i][k] = NULL;
 // }
 // CovMatrixIndexActive = false;

   // init models
  for (i=0;i<MAXKEYS;i++) KEY_NULL(KEY + i);

  if (CovList!=NULL) {
    PRINTF("List of covariance functions looks already initiated.\n"); 
    return;
  }
  CovList = (cov_fct*) malloc(sizeof(cov_fct) * (MAXNRCOVFCTS+1));
  // + 1 is necessary because of COVINFO_NULL that uses the last + 
  currentNrCov = 0;

  PLUS = IncludeModel("+", 2, MAXSUB, 0, NULL, PREVMODELS, PREVMODELI, checkplus, 
	       rangeplus, PREF_ALL, initplus, doplus, false, SUBMODELM);
  addCov(plusStat, Dplus, DDplus, NULL);
  addCov(plusNonStat);
  addTBM(initspectralplus, spectralplus);
  MLEPLUS = addFurtherCov(MLEplusStat, ErrCov); 
  addCov(MLEplusNonStat);

  pref_type pmal =  {5, 0, 0, 3, 5, 0, 5, 5, 0, 0, 0, 0, 0, 5};
  //                CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludeModel("*", 2, MAXSUB, 0, NULL, PREVMODELS, PREVMODELI, checkmal, 
	       rangemal, pmal, NULL, NULL, false, SUBMODELM);
  addCov(malStat, Dmal, NULL);
  addCov(malNonStat);

  
  pref_type pS= {5, 5, 5, 3, 5, 5, 5, 5, 0, 0, 5, 0, 0, 5};
 	//       CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  DOLLAR = IncludeModel("$", 1, 1, 5, kappaS, // kappadollar,
			PREVMODELS, PREVMODELI, checkS, rangeS, pS,
			initS, doS, false, SUBMODELM);
  // do not change Order!!
  kappanames("var", REALSXP, "scale", REALSXP, "anisoT", REALSXP,
	     "A", REALSXP, "proj", INTSXP);
  addCov(Siso, DS, DDS, NULL);
  addLocal(coinitS, ieinitS);
  addTBM(tbm2S, initspectralS, spectralS);
  nablahess(nablaS, hessS);
  LASTDOLLAR = addFurtherCov(Sstat, ErrCov); // 4
  addCov(Snonstat);
  nablahess(nablaS, hessS);
  addInv(invS);
  
  // for speed MLE only -- not programmed yet 
  int mleDollar = addFurtherCov(Siso_MLE, DS_MLE, DDS_MLE); // 5
  addCov(Snonstat_MLE);
  if (mleDollar != LASTDOLLAR + (LASTDOLLAR - DOLLAR)) 
    error("mleDollar not at the end of $");
 

  ISO2ISO = GATTER = // 6
      IncludeModel("#", 1, 1, 0, NULL, PREVMODELS, PREVMODELI, check2, range2,
		   PREF_ALL, init2, do2, true, SUBMODELM);
  addCov(iso2iso, D_2, DD_2, NULL);
  addTBM(initspectral2, spectral2);
  addMPP(mppinit2, coin2, mppgetiso2iso, sd2, MPP_AUTO);
  addInv(inv2);
 
  SP2SP = addFurtherCov(spiso2spiso, D_2, DD_2); // 7
  addTBM(initspectral2, spectral2);
  addMPP(mppinit2, coin2, mppgetspiso2spiso, sd2, MPP_AUTO);

  SP2ISO = addFurtherCov(spacetime2iso, D_2, DD_2); // 8
  addTBM(initspectral2, spectral2);
  addMPP(mppinit2, coin2, mppgetspacetime2iso, sd2, MPP_AUTO);

  S2ISO = addFurtherCov(Stat2iso, ErrCov); // 9
  addCov(Nonstat2iso);// 
  addTBM(initspectral2, spectral2);
  addMPP(mppinit2, coin2, mppgetStat2iso, sd2, MPP_AUTO);
  addMPP(mppinit2, coin2, mppgetNonstat2iso, sd2, MPP_AUTO);
  assert(CovList[S2ISO].mppinit != NULL);

  assert(S2ISO == 9);

  S2SP = addFurtherCov(Stat2spacetime, ErrCov);// 10
  addCov(Nonstat2spacetime);// 
  addTBM(initspectral2, spectral2);
  addMPP(mppinit2, coin2, mppgetStat2spacetime, sd2, MPP_AUTO);
  addMPP(mppinit2, coin2, mppgetNonstat2spacetime, sd2, MPP_AUTO);
 
  S2S = addFurtherCov(Stat2Stat, ErrCov);// 11
  addCov(Nonstat2Stat);// 
  addTBM(initspectral2, spectral2);
  addMPP(mppinit2, coin2, mppgetStat2Stat, sd2, MPP_AUTO);
  addMPP(mppinit2, coin2, mppgetNonstat2Stat, sd2, MPP_AUTO);

  assert(S2S == 11); //meth->hanging->nr  in InternalCov.cc
  LASTGATTER = S2S;

  addFurtherCov(iso2iso_MLE, D_2, DD_2); // 12
  addFurtherCov(spiso2spiso_MLE, D_2, DD_2); // 13
  addFurtherCov(spacetime2iso_MLE, D_2, DD_2); // 14
  addFurtherCov(Stat2iso_MLE, ErrCov); // 15
  addCov(Nonstat2iso_MLE);// 
  addFurtherCov(Stat2spacetime_MLE, ErrCov);// 16
  int mleGatter = addFurtherCov(Stat2Stat_MLE, ErrCov);// 17
  addCov(Nonstat2Stat_MLE);// 
  if (mleGatter != LASTGATTER + (LASTGATTER - GATTER + 1))
    error("mleGatter not at the end of #"); // see MLE !


//  MLE_ENDCOV = 
//      IncludeModel("^", 0, 0, 0, NULL, PREVMODELS, PREVMODELI, checkHat, 
//		   rangeHat,
//		   (pref_type) {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, NULL, 
//		   NULL, true, PARAMETERM);
//  addCov(Hatstat, NULL, NULL, NULL);
//  addCov(Hatnonstat);

//  OUT =
  //     IncludeModel("0", 0, 0, 0, NULL, PREVMODELS, PREVMODELI, checkOut, 
//		   rangeOut,
//	      (pref_type) {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, NULL, 
//	       NULL, true, PARAMETERM);
//  addCov(Outstat, NULL, NULL, NULL);
//  addCov(Outnonstat);


//  SPLIT =  IncludeModel("|", 1, MAXSUB, 1, NULL, PREVMODELS, PREVMODELI, 
//	       checksplit, rangesplit,
//	       (pref_type) {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, NULL, 
//	       NULL, true, PARAMETERM);
//  kappanames("part", INTSXP);
//  // subnames("fix", "fixed");
//  addCov(splitstat, NULL, NULL, NULL);
//  addCov(splitnonstat);

  IncludeModel("ave1", 1, 1, 2, kappa_ave1, STATIONARY, ANISOTROPIC, checkave1,
	       rangeave1, PREF_ALL);
  kappanames("A", REALSXP, "z", REALSXP);
  addCov(ave1, NULL, NULL);
  addMPP(mppinit_ave1, coin_ave, mppget_ave1, sd_ave_stp, 
	 MPP_GRID, ave1_f, ave1_logg);
  subnames("phi");


  IncludeModel("ave2", 1, 1, 2, kappa_ave2, STATIONARY, ANISOTROPIC, checkave2,
	       rangeave2, PREF_ALL);
  kappanames("A", REALSXP, "z", REALSXP);
  addCov(ave2, NULL, NULL);
  addMPP(mppinit_ave2, coin_ave, mppget_ave2, sd_ave_stp, MPP_GAUSS);
  subnames("phi");

  pref_type pbessel = {2, 0, 0, 0, 0, 5, 4, 5, 0, 5, 0, 5, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludePrim("bessel", 1, STATIONARY, ISOTROPIC, checkBessel, rangeBessel,
	      pbessel, UNIVARIATE);
  kappanames("nu", REALSXP);
  addCov(Bessel, NULL, NULL);
  addTBM(initspectralBessel, spectralBessel);	       

//  IncludePrim("biWM0", 3, kappa_biWM, STATIONARY, ISOTROPIC, checkbiWM,
//	      rangebiWM, 2);
//  addCov(biWM, NULL, NULL);
//  kappanames("nu", REALSXP, "a", REALSXP, "c", REALSXP);

 
  IncludePrim("biWM", 6, kappa_biWM, STATIONARY, ISOTROPIC, checkbiWM2,
	      rangebiWM2, 2);
  addCov(biWM2, NULL, NULL);
  kappanames("nu", REALSXP, "nured12", REALSXP, 
	     "s", REALSXP, "s12", REALSXP,
	     "c", REALSXP, "rhored", REALSXP);



  pref_type pchauchy=  {2, 0, 0, 0, 3, 0, 4, 5, 0, 0, 0, 0, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludePrim("cauchy", 1, STATIONARY, ISOTROPIC, checkCauchy, rangeCauchy,
	      pchauchy, UNIVARIATE);
  kappanames("beta", REALSXP);
  addCov(Cauchy, DCauchy, DDCauchy, ScaleCauchy);
  addTBM(TBM2Cauchy);
  addLocal(coinitCauchy, NULL);
  addGaussMixture(DrawLogMixCauchy, LogMixWeightCauchy);
  addInv(invCauchySq);
	       
  pref_type pctbm={2, 0, 0, 0, 5, 0, 5, 5, 0, 0, 0, 0, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludePrim("cauchytbm", 3, STATIONARY, ISOTROPIC, checkOK,
	      rangeCauchytbm, pctbm, UNIVARIATE);
  kappanames("alpha", REALSXP, "beta", REALSXP, "gamma", REALSXP);
  addCov(Cauchytbm, DCauchytbm, ScaleCauchy); // scale not correct, but
  // should be an approximation that is good enough
	      
  IncludePrim("circular", 0, STATIONARY, ISOTROPIC, checkOK,
		 rangecircular);
  addCov(circular, Dcircular, ScaleOne);
  addMPP(mppinit_circular, coin_circular, mppget_spherical, 
	 sd_standard, MPP_POISS);

  IncludePrim("cone", 3, STATIONARY, ISOTROPIC, checkcone, rangecone);
  kappanames("r", REALSXP, "socle", REALSXP, "height", REALSXP);
  addMPP(mppinit_cone, coin_cone, mppget_cone, sd_standard, MPP_POISS);

  CONSTANT = 
      IncludeModel("constant", 0, 0, 2, NULL, PREVMODELS, PREVMODELI, 
		   checkconstant, rangeconstant, PREF_ALL,
		   NULL, NULL, false, PARAMETERM);
  kappanames("M", LISTOF+REALSXP, "vdim", INTSXP);  
  addCov(constant, NULL, NULL); 
  addCov(constant_nonstat);

  pref_type pcox={2, 0, 5, 0, 0, 5, 5, 5, 0, 0, 0, 0, 0, 5};
     //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludeModel("coxisham", 1, 1, 3, kappa_cox1, STATIONARY, ZEROSPACEISO, 
	      checkcox1, rangecox1, pcox);
  kappanames("mu", REALSXP, "D", REALSXP, "beta", REALSXP);  
  addCov(cox1, NULL, NULL);
  addTBM(initspectralcox1, spectralcox1);
  nablahess(coxnabla, coxhess);

  IncludePrim("cubic", 0, STATIONARY, ISOTROPIC, checkOK, rangecubic);
  addCov(cubic, Dcubic, ScaleOne);
	       
  pref_type pcurl= {2, 0, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludeModel("curlfree", 
	      1, 1, 0, NULL, STATIONARY, ANISOTROPIC, checkdivcurl,
	       rangedivcurl, pcurl,
	      NULL, NULL, false, PARAMETERM);
  addCov(curl, NULL, NULL);
 
//  warning("siehe Primitive.cc/biWM: cutoff funktioniert nicht bei MLE, vereinheitlichung mit natsc und verbesserung von biWM notwendig");
  IncludeModel("cutoff", 1, 1, 2, STATIONARY, ISOTROPIC,
	       check_co, range_co, PREF_ALL);
  kappanames("diameter", REALSXP, "a", REALSXP);
  addCov(co, NULL, NULL);
  addCallLocal(alternativeparam_co);
  subnames("phi");


  IncludePrim("dagum", 2, STATIONARY, ISOTROPIC, checkOK,  
	      rangedagum);
  kappanames("beta", REALSXP, "gamma", REALSXP);
  addCov(dagum, Ddagum, Scaledagum);

  IncludePrim("dampedcosine", 1, STATIONARY, ISOTROPIC, checkOK,  
		  rangedampedcosine);
  kappanames("lambda", REALSXP);
  addCov(dampedcosine, Ddampedcosine, Scaledampedcosine);

  IncludeModel("delayeffect", 
	       1, 1, 1, kappashift, STATIONARY, ANISOTROPIC, checkshift,
	       rangeshift, PREF_ALL, NULL, NULL, false, 2);
  addCov(shift, NULL, NULL);
  kappanames("r", REALSXP);


  IncludePrim("DeWijsian", 1, VARIOGRAM, ISOTROPIC, checkOK,
		 rangedewijsian);
  kappanames("alpha", REALSXP);
  addCov(dewijsian, NULL, NULL); 

  IncludePrim("CDeWijsian", 2, NULL, STATIONARY, ISOTROPIC, checkOK,
	      rangeDeWijsian, PREF_NOTHING, 1); 
  make_internal();
  kappanames("alpha", REALSXP, "range", REALSXP);
  addCov(DeWijsian, NULL, ScaleDeWijsian); 


  pref_type pdiv= {2, 0, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
 IncludeModel("divfree", 
	      1, 1, 0, NULL, STATIONARY, ANISOTROPIC, checkdivcurl,
	      rangedivcurl, pdiv,
	      NULL, NULL, false, PARAMETERM);
  addCov(div, NULL, NULL);
  IncludePrim("EAxxA", 2, kappa_EAxxA, AUXMATRIX, ANISOTROPIC,
	      checkEAxxA, rangeEAxxA, PREF_NOTHING, PARAMETERM);
  addCov(EAxxA, NULL, NULL);
  kappanames("E", REALSXP, "A", REALSXP);
  addSpecial(minEigenEAxxA);

  // epsC has been for internal reasons only ; essentially
  // the gencauchy model, except that 1 in the denominator 
  // is replaced by epsilon
  pref_type pepsC = {2, 5, 5, 0, 5, 0, 5, 5, 0, 0, 0, 0, 0, 5};
	       //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludeModel("epsC", 0, 0, 3, NULL, STATIONARY, ISOTROPIC,
	       checkepsC, rangeepsC, pepsC,
	       NULL, NULL, true, UNIVARIATE);
  kappanames("alpha", REALSXP, "beta", REALSXP, "eps", REALSXP);
  addCov(epsC, DepsC, DDepsC, NULL);

  IncludePrim("EtAxxA", 3, kappa_EtAxxA, AUXMATRIX, ANISOTROPIC,
	      checkEtAxxA, rangeEtAxxA, PARAMETERM);
  addCov(EtAxxA, NULL, NULL);
  kappanames("E", REALSXP, "A", REALSXP, "alpha", REALSXP);
  addSpecial(minEigenEtAxxA);

  IncludePrim("exponential", 0, STATIONARY, ISOTROPIC, checkexponential,
		  rangeexponential);
  addCov(exponential, Dexponential, DDexponential, Scaleexponential);
  addTBM(TBM2exponential, initspectralexponential, spectralexponential);
  addLocal(coinitExp, ieinitExp);
  addMarkov(&EXPONENTIAL);
  addHyper(hyperexponential);
  // numerisches Ergebnis passt nicht !
  addGaussMixture(DrawLogMixExp, LogMixWeightExp);
  addInv(invexponentialSq);
  
  // operator, replacing any covariance fct C by exp(C) (componentwise)
  IncludeModel("Exp", 1, 1, 0, NULL, PREVMODELS, PREVMODELI, checkExp,
	       rangeExp, PREF_ALL);
  addCov(Exp, DExp, DDExp, NULL);
  
 
  IncludePrim("FD", 1, STATIONARY, ISOTROPIC, checkOK, rangeFD);
  kappanames("a", REALSXP);
  addCov(FD, NULL, NULL);

  IncludePrim("fractalB", 1, VARIOGRAM, ISOTROPIC, checkOK,
		  rangefractalBrownian);
  kappanames("alpha", REALSXP);
  addCov(fractalBrownian, DfractalBrownian, DDfractalBrownian, NULL);
  addLocal(NULL, ieinitBrownian);
  addInv(invfractalBrownianSq);
  
  IncludePrim("fractgauss", 1, STATIONARY, ISOTROPIC, checkOK,  
		  rangefractGauss);
  kappanames("alpha", REALSXP);
  addCov(fractGauss, NULL, NULL);

  pref_type pgauss= {2, 0, 0, 0, 5, 5, 5, 5, 5, 0, 0, 5, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  GAUSS = 
      IncludePrim("gauss", 0, STATIONARY, ISOTROPIC, checkOK, rangeGauss,
		  pgauss, UNIVARIATE);
  addCov(Gauss, DGauss, DDGauss, D3Gauss, D4Gauss, ScaleGauss);
  addTBM(initspectralGauss, spectralGauss);
  addMarkov(&GAUSS);
  addMPP(mppinit_Gauss, coin_Gauss, mppget_Gauss, sd_standard, MPP_POISS);
  addGaussMixture(DrawLogMixGauss, LogMixWeightGauss);
  addInv(invGaussSq);


  IncludePrim("gausstest", 0, NULL, STATIONARY, ANISOTROPIC, checkOK, rangeGauss,
	      pgauss, UNIVARIATE);
  make_internal();
  addMPP(mppinit_Gausstest, coin_Gausstest, mppget_Gausstest, sd_standard,
	  MPP_POISS);
 
  IncludePrim("genB", 2, VARIOGRAM, ISOTROPIC, checkOK,
	      rangegenBrownian);
  kappanames("alpha", REALSXP, "delta", REALSXP);
  addCov(genBrownian, NULL, NULL); 
  
  pref_type pgenc = {2, 5, 5, 0, 5, 0, 5, 5, 0, 0, 0, 0, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludePrim("gencauchy", 2, STATIONARY, ISOTROPIC, checkgeneralisedCauchy,
	      rangegeneralisedCauchy, pgenc,
	      UNIVARIATE);
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(generalisedCauchy, DgeneralisedCauchy, DDgeneralisedCauchy,
	 ScalegeneralisedCauchy);
  addLocal(coinitgenCauchy, ieinitgenCauchy);
  addInv(invgeneralisedCauchySq);

  IncludePrim("gengneiting", 2, STATIONARY, ISOTROPIC, checkOK,  
		  rangegenGneiting);
  kappanames("n", INTSXP, "alpha", REALSXP);
  addCov(genGneiting, DgenGneiting, ScaleOne);

  IncludePrim("gneiting", 0, STATIONARY, ISOTROPIC, checkOK,
		   rangeGneiting); 
  addCov(Gneiting, DGneiting, ScaleOne);
 
  pref_type phyper= {2, 0, 0, 0, 3, 0, 4, 5, 0, 5, 0, 5, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludePrim("hyperbolic", 3, STATIONARY, ISOTROPIC, checkhyperbolic,
	      rangehyperbolic, phyper,
	      UNIVARIATE);
  kappanames("nu", REALSXP, "lambda", REALSXP, "delta", REALSXP);
  addCov(hyperbolic, Dhyperbolic, NULL);
  addInv(invhyperbolicSq);

  IncludePrim("iacocesare", 3, STATIONARY, SPACEISOTROPIC, checkOK,
		  rangeIacoCesare);
  kappanames("nu", REALSXP, "lambda", REALSXP, "delta", REALSXP);
  addCov(IacoCesare, NULL, NULL);
  
//  IncludeModel("identity", 1, 1, 1, NULL, PREVMODELS, PREVMODELI, checkId,
//	       rangeId,	PREF_ALL, initstandard, dostandard, true, SUBMODELM);
//  kappanames("vdim", INTSXP);
//  addCov(IdStat, DId, DDId, NULL);
//  addTBM(TBM2Id, initspectralId, spectralId);
//  addLocal(coinitId, ieinitId);
//  //  addMarkov(&ID);
//  addCov(IdNonStat);

  VECTOR = IncludeModel("vector", 
	       1, 1, 2, NULL, STATIONARY, ANISOTROPIC, checkvector,
	       rangevector, PREF_ALL, NULL, NULL, false, PARAMETERM);
  addCov(vector, NULL, NULL);
  kappanames("a", REALSXP, "Dspace", INTSXP);
  addFurtherCov(vectorAniso, NULL); 

  pref_type plgd1= {2, 0, 0, 0, 5, 0, 5, 5, 0, 0, 0, 0, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludePrim("lgd1", 2, STATIONARY, ISOTROPIC, checkOK, rangelgd1, plgd1,
	      UNIVARIATE);
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(lgd1, Dlgd1, Scalelgd1);


// stimmt so nicht, siehe Gneiting 1998, on a alpha-sym multiv. char. fct:
//  IncludeModel("lp", 1, 1, 1, STATIONARY, ANISOTROPIC, 
//	       checklp, rangelp,
//	       (pref_type) {5, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 5}
//	       //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
//             );
//  kappanames("p", REALSXP);
//  addCov(lp, NULL, NULL);
 
  IncludeModel("mastein", 1, 1, 2, STATIONARY, SPACEISOTROPIC, 
	       check_MaStein, range_MaStein, PREF_ALL);
  kappanames("nu", REALSXP, "delta", REALSXP);
  addCov(MaStein, NULL, NULL);
  subnames("phi");
 

  IncludeModel("M", 1, 1, 1, kappaM,
	       PREVMODELS, PREVMODELI, checkM, rangeM, PREF_ALL,
	       initstandard, dostandard, false, PARAMETERM);
  kappanames("M", REALSXP);
  addCov(Mstat, NULL, NULL);
  addCov(Mnonstat);

	       
  IncludeModel("ma1", 1, 1, 2, STATIONARY, ANISOTROPIC,
	       checkma1, rangema1, PREF_ALL);
  kappanames("alpha", REALSXP, "theta", REALSXP);
  addCov(ma1, NULL, NULL);


  IncludeModel("ma2", 1, 1, 0, STATIONARY, ANISOTROPIC,
	       checkma2, rangema2, PREF_ALL);
  addCov(ma2, NULL, NULL);


  IncludeModel("matern", 0,0, 2, STATIONARY, ISOTROPIC, 
	       checkMatern, rangeMatern, PREF_ALL);
  kappanames("nu", REALSXP, "invnu", INTSXP);
  addCov(Matern, DMatern, DDMatern, D3Matern, D4Matern, NULL);
  addTBM(initspectralMatern, spectralMatern);
  addLocal(coinitMatern, ieinitMatern);
  addInv(invMaternSq);

//  addGaussMixture(DrawLogMixWM, LogMixWeightWM);

  
  pref_type pmixed ={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5};
		   //          CE CO CI TBM23 Sp di sq Ma av  n mpp Hy any
  MIXEDEFFECT = 
      IncludeModel("mixed", 0, 1, 2, PREVMODELS, PREVMODELI,
		   checkmixed, rangemixed, pmixed);
  kappanames("X", LISTOF+REALSXP, "b", REALSXP); // #define BETAMIXED 1
  subnames("covb");
  addCov(mixed, NULL, NULL);
  addCov(mixed_nonstat);
  MLEMIXEDEFFECT = addFurtherCov(MLEmixed, ErrCov); 
  addCov(MLEmixed_nonstat);

 

  IncludeModel("mqam", 2, MAXSUB, 2, kappamqam, STATIONARY, ANISOTROPIC,
	       checkmqam, rangemqam, PREF_ALL, NULL, NULL, false, PARAMETERM);
  kappanames("theta", REALSXP, "rho", REALSXP);
  subnames("phi");
  addCov(mqam, NULL, NULL);


  pref_type pnatsc= {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
			//          CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  NATSC = IncludeModel("natsc", 1, 1, 0, NULL, 
		       STATIONARY, ISOTROPIC, checknatsc, rangenatsc, pnatsc,
			NULL, NULL, false, 1);
  // do not change Order!!
  addCov(natsc, Dnatsc, DDnatsc, NULL);
  addLocal(coinitnatsc, ieinitnatsc);
  addTBM(tbm2natsc, initspectralnatsc, spectralnatsc);
  addInv(invnatsc);
  assert(CovList[NATSC].naturalscale == NULL); // needed in CMbuild


  IncludeModel("nonstWM", 0, 1, 1,  kappaNonStWM, COVARIANCE, ANISOTROPIC, 
	       checkNonStWM, rangeNonStWM, PREF_ALL);
  addCov(NonStWMQ); // anders herum gibt es fehler in addCov(aux_covfct auxcf),
  addCov(NonStWM);  // da auxiliary ia.. nicht mit cov hand in hand gehen kann
  kappanames("nu", REALSXP);
  subnames("Nu");
  addGaussMixture(DrawLogMixNonStWM, LogMixWeightNonStWM);


  IncludeModel("nsst", 2, 2, 1, STATIONARY, SPACEISOTROPIC,
	       checknsst, rangensst, PREF_ALL);
  kappanames("delta", REALSXP);
  subnames("phi", "psi");
  addCov(nsst, Dnsst, NULL);
  addTBM(TBM2nsst);
  
//  IncludePrim("nsst2", 7, checknsst2, SPACEISOTROPIC, 
//		   rangensst2);
//  addCov(nsst2, Dnsst2, NULL);
//  addTBM(NULL, NULL /* TBM3nsst2 */);


  pref_type pnugget= {4, 0, 0, 0, 0, 0, 4, 4, 0, 0, 5, 0, 0, 5};
	       //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  NUGGET = IncludeModel("nugget", 0, 0, 1, NULL,
	       STATIONARY, ISOTROPIC, checknugget, rangenugget, pnugget,
	       NULL, NULL, false, PARAMETERM);
  kappanames("vdim", INTSXP);
  addCov(nugget, NULL, Scalenugget);
  /* see simu.cc, CMbuild for special treatment of nugget when
     users choice is given */
  
  /* cf. convert.R, PrepareModel, near end of function */
  /* deleted 27.12.08
    { 
    int nr = currentNrCov - 1;
    for (i=0; i<SimulationTypeLength; i++) {
      CovList[nr].implemented[i] = GIVEN_METH_IGNORED;
    }
    CovList[nr].implemented[Nugget] = 
      CovList[nr].implemented[Direct] =
      CovList[nr].implemented[Sequential] =
      CovList[nr].implemented[CircEmbed] =IMPLEMENTED;
  }
  */

 IncludePrim("parsbiWM", 4, kappa_parsbiWM, STATIONARY, ISOTROPIC, 
	      checkparsbiWM2, rangeparsbiWM2, 2);
  addCov(parsbiWM2, NULL, NULL);
  kappanames("nu", REALSXP, 
	     "s", REALSXP,
	     "c", REALSXP, 
	     "rhored", REALSXP);

  IncludePrim("penta", 0, STATIONARY, ISOTROPIC, checkOK, 
                rangepenta);
  addCov(penta, Dpenta, ScaleOne);

	
  IncludePrim("power", 1, STATIONARY, ISOTROPIC, checkpower, 
		  rangepower);
  kappanames("a", REALSXP);
  addCov(power, Dpower, ScaleOne);
 
 IncludeModel("Pow", 1, 1, 1, NULL, PREVMODELS, PREVMODELI, checkPow,
	       rangePow, PREF_ALL);
 addCov(Pow, DPow, DDPow, NULL); 
 kappanames("alpha", REALSXP);
 addInv(invPow);

 IncludeModel("qam", 2, MAXSUB, 1, kappaqam, STATIONARY, ANISOTROPIC,
	       checkqam, rangeqam, PREF_ALL);
  kappanames("theta", REALSXP);
  subnames("phi");
  addCov(qam, NULL, NULL);


  IncludePrim("qexponential", 1, STATIONARY, ISOTROPIC, checkOK,
	      rangeqexponential);
  kappanames("alpha", REALSXP);
  addCov(qexponential, Dqexponential, Scaleqexponential);

  IncludePrim("rational", 2, kappa_rational, AUXMATRIX, ANISOTROPIC,
	      checkrational, rangerational);
  addCov(rational, NULL, NULL);
  kappanames("A", REALSXP, "a", REALSXP);
  addSpecial(minEigenrational);


  IncludePrim("rotat", 2, kappa_rotat, AUXMATRIX, ANISOTROPIC,
	      checkrotat, rangerotat, PREF_NOTHING, UNIVARIATE);
  make_internal();
  addCov(rotat, NULL, NULL);
  kappanames("speed", REALSXP, "phi", REALSXP);
  addSpecial(minEigenrotat);

  IncludePrim("Rotat", 1, kappa_Rotat, ANYFCT, ANISOTROPIC,
	      checkRotat, rangeRotat, PARAMETERM);
  addCov(Rotat, NULL, NULL);
  kappanames("phi", REALSXP);



  IncludePrim("spherical", 0, STATIONARY, ISOTROPIC, checkOK, 
		  rangespherical);
  addCov(spherical, Dspherical, ScaleOne);
  addTBM(TBM2spherical);
  addMPP(mppinit_spherical, coin_spherical, mppget_spherical,
	 sd_standard, MPP_POISS);

  IncludePrim("stable", 1, STATIONARY, ISOTROPIC, checkstable,
	      rangestable);
  kappanames("alpha", REALSXP);
  addCov(stable, Dstable, DDstable, Scalestable);
  addLocal(coinitstable, ieinitstable);
  addInv(invstableSq);

  // SPACEISOTROPIC variant of stable -- used for testing purposes only
//  IncludePrim("stableX", 1, checkOK, SPACEISOTROPIC, 
//		  rangestable);
//  addCov(stableX, DstableX, Scalestable);
//  addTBM(NULL, NULL);


  IncludeModel("Stein", 1, 1, 2, STATIONARY, ISOTROPIC, 
	       check_Stein, range_Stein, PREF_ALL);
  kappanames("Diameter", REALSXP, "rawR", REALSXP);
  addCov(Stein, NULL, NULL);
  addCallLocal(alternativeparam_Stein);

  IncludePrim("steinst1", 2, kappaSteinST1, STATIONARY, ANISOTROPIC, 
	      checkSteinST1, rangeSteinST1);
  kappanames("nu", REALSXP, "z", REALSXP);
  addCov(SteinST1, NULL, NULL);
  addTBM(initspectralSteinST1, spectralSteinST1);
 

  IncludeModel("stp", 1, 4, 3, kappa_stp, COVARIANCE, ANISOTROPIC,
	       checkstp, rangestp, PREF_ALL);
  addCov(stp);
  kappanames("S", REALSXP, "M", REALSXP, "z", REALSXP);
  addMPP(mppinit_stp, coin_stp, mppget_stp, sd_ave_stp, MPP_GAUSS);
  subnames("phi", "Sx", "xi2", "H"); // H ueberall wo U-x steht. dort U-H(x)

  TBM2NR = IncludeModel("tbm2", 1, 1, 0, NULL,
	       STATIONARY, PREVMODELI, checktbm2, rangetbm2, PREF_ALL,      
	       initstandard, dostandard, true, UNIVARIATE);
  addCov(tbm2, NULL, NULL);
  addFurtherCov(tbm2num, NULL);

  IncludeModel("tbm3", 1, 1, 1, NULL,
	       STATIONARY, PREVMODELI, checktbm3, rangetbm3, PREF_ALL,
	       initstandard, dostandard, true, UNIVARIATE);
  addCov(tbm3, NULL, NULL);
  kappanames("n", REALSXP);
 


  /*
  int UserMatrix =
      IncludeModel("userMatrix", 0, 0, 2, kappa_userMatrix, 
	       STATIONARY, ISOTROPIC,
	       checkuserMatrix, rangeuserMatrix,
	       (pref_type) {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5},
	       //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
	       NULL, NULL, false, UNIVARIATE);
  kappanames("nr", INTSXP, "Z", INTSXP);
  addCov(userMatrix, NULL, NULL);
  */


  pref_type pwave = {2, 0, 0, 0, 0, 5, 4, 5, 0, 0, 0, 0, 0, 5};
	      //           CE CO CI TBM23 Sp di sq Ma av n mpp Hy any
  IncludePrim("wave", 0, STATIONARY, ISOTROPIC, checkOK, rangewave, pwave,
	      UNIVARIATE);
  addCov(wave, NULL, Scalewave);
  addTBM(initspectralwave, spectralwave);
 
  IncludePrim("whittle", 1, STATIONARY, ISOTROPIC, checkWhittle,
		  rangeWhittle);
  kappanames("nu", REALSXP);
  addCov(Whittle, DWhittle, DDWhittle, D3Whittle, D4Whittle, ScaleWhittle);
  addTBM(initspectralWhittle, spectralWhittle);
  addLocal(coinitWhittle, ieinitWhittle);
  addMarkov(&WHITTLE);
  addGaussMixture(DrawLogMixWM, LogMixWeightW);
  addInv(invWhittleSq);

 
  pref_type pX= {0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 5};
		   //          CE CO CI TBM23 Sp di sq Ma av  n mpp Hy any
   MIXX = IncludeModel("X", 1, 1, 1,  PREVMODELS, PREVMODELI, checkX, rangeX, 
       pX);
  kappanames("X", LISTOF+REALSXP);
  subnames("covb");
  addCov(X, NULL, NULL);
  addCov(Xnonstat);

  // addMarkov(MarkovWhittle);
// Markovcircular, Markovexponential, Markovcubic, Markovdampedcosine, MarkovFD, 
// Markovgauss, Markovgneiting, Markovgengneiting, Markovhyperbolic, Markovpenta, 
// Markovpower, Markovqexponential, Markovspherical, Markovstable, 
  

  // must be the very last one!
  int oldcurrentNrCov = currentNrCov;
  assert(currentNrCov < MAXNRCOVFCTS); // otherwise we have already reached
  //                the maximum number
  //                of models; note that it is planned that further
  //                 functions might be added by users
  IncludePrim("undefined", 0, STATIONARY, ISOTROPIC, checkundefined,
		   rangenugget);
  while (currentNrCov < MAXNRCOVFCTS) addFurtherCov(NULL, NULL);
  currentNrCov = oldcurrentNrCov;

  addusersfunctions();

 
  cov_fct *C;
  cov_model *cov = NULL;
  char biWM[] = "biWM";
  int nr, *implement,
    biwm = getmodelnr(biWM);
  if (PL > 4) PRINTF("init: checking model definitions\n");
  for (nr=0; nr<currentNrCov; nr++) {
    
    // correct pref according to given functions
    C = CovList + nr; 
    implement = C->implemented;
//  
 
    // if (nr == UserMatrix) continue;
    
    if ((nr < GATTER || nr > LASTGATTER) && 

	(nr < DOLLAR || nr > LASTDOLLAR) &&
	nr != NUGGET) { 	

      // # has always all preferences; nugget is special case, not checked
      for (i=0; i<Nothing; i++) {
//	printf("%s %d ", METHODNAMES[i] , C->pref[i]);
	C->pref[i] *= implement[i] == IMPLEMENTED || implement[i] == NUM_APPROX;
//	printf("%d\n ", C->pref[i]);
      }
      C->pref[Nothing] =
	  PREF_BEST * (C->cov != ErrCov || C->nonstat_cov != ErrCovNonstat);
    }
//{ int i; for (i=0; i<Forbidden; i++) printf("%d %d\n", i, CovList[0].pref[i]);
    //     assert(false);}
 
    // general check

    if (nr == biwm || nr == CONSTANT || nr == MIXEDEFFECT || nr == MLEMIXEDEFFECT
	|| nr == MIXX) 
	continue;

   cov = (cov_model *) malloc(sizeof(cov_model));
    COV_NULL(cov);
    cov->nr = nr;
    cov->maxdim = 1;
  
    for (i=0; i<C->kappas; i++) {
      cov->p[i] = (double*) malloc(sizeof(double) * 4);
      int j;
      for (j=0; j<4; j++) cov->p[i][j] =0.1;
      cov->ncol[i] = cov->nrow[i] = 2;
    }

    if (C->primitive && C->stationary!=AUXMATRIX && C->stationary!=ANYFCT) {
      for (i=0; i<Nothing; i++) cov->user[i] = cov->pref[i] = 0;
      cov->tsdim = cov->xdim =  (C->cov == SteinST1) ? 2 : 1;
      //  printf("ok %ld %ld\n", C->check, checkBessel);
      C->check(cov); // no check on error;
      //     printf("done\n");
      for (i=0; i<Nothing; i++) {
	if(!C->implemented[i] && cov->pref[i]>0) {
	  PRINTF("%d: %s  pos=%d, value %d > 0\n", nr, C->name, i,
		 cov->pref[i]);
	  assert(false);
	}
      }
    }
    COV_DELETE_WITHOUTSUB(&cov);
  }
  
  if (PL > 4) PRINTF("init: end of checking model definitions\n");
}

