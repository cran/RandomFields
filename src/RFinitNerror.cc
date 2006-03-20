
/* 
 Authors
 Martin Schlather, schlath@hsu-hh.de

 library for simulation of random fields -- init part and error messages

 Copyright (C) 2001 -- 2006 Martin Schlather

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
#include <assert.h>
#include <string.h>
#include "RFsimu.h"
#include "RFCovFcts.h"
#include "MPPFcts.h"
#include "Hyperfcts.h"
#include <unistd.h>

initmethod init_method[SimulationTypeLength];
do_comp_meth do_compatiblemethod[SimulationTypeLength];
do_incomp_meth  do_incompatiblemethod[SimulationTypeLength];
key_type KEY[MAXKEYS]; //static key_type KEY[MAXKEYS];
double ZERO[MAXDIM], UNIT[MAXDIM];
cov_fct *CovList=NULL;
int currentNrCov=-1;

int  GENERAL_SKIPCHECKS=false;
char GENERAL_PCH[2]="*"; 
/*  character printed after each simulation
    just for entertainment of the user
    except for "!", then the numbers are shown
*/
int GENERAL_STORING=false; 
/* true: intermediate results are stored: might be rather memory consuming,
         but the simulation can (depending on the method chosen) be much faster
	 when done the second time with exactly the same parameters
	 do not forget to call DeleteAllKeys when all the simulation are done.
   false: intermediate results are not stored if SimulateRF is called or
         stored only until DoSimulateRF is called.
	 minimum memory consumption, but might be rather slow if many
	 simulations are performed with exactly the same parameters
	 to be safe call DeleteAllKeys when all the simulations are done
	 Warning! init_foo may depend on GENERAL_STORING! This may cause problems
	 if GENERAL_STORING is changed between init_foo and do_foo from false to
	 true and if do_foo is called more than once! (SimulateRF is safe in
	 this case)

	 indifferent behaviour for the simulation methods if parameters are 
	 changed after each simulation.
	 In case, only MAXKEYS different parameter sets are used, but after each 
	 step the parameter set is changed, use different keynr for each 
	 parametere set, and STORING==true, to get fast simulation 
*/
int GENERAL_PRINTLEVEL=1;
/* 0:no messages; 
   1:error messages; 
   2: hints when algorithm is at a forcation; 
   >=3: more and more debugging information
*/

int GENERAL_NATURALSCALING=0; 
/* is called PracticalRange in R (see Chiles&Delfiner?!) 
   Note: RFparameters() allows only 0 and 1 as values!

   has an effect only if cov (and not only cov.local, e.g. Brownian Motion) is 
   defined
   0 : using the covariance function as defined
   1 : rescaling of cov fctn such that cov(1)~=0.05, if rescaling function 
       does not exist then failure, e.g. for Bessel model
   2 : exact or approximate value (e.g. expPLUScirc)
   3 : MLE (special needs taken into account, long memory covariance functions
       have too long tails, which are shortened (so threshold 0.05 is changed to
       higher values in this case))
  +10: if any of the above fails : numerical evaluation!
   else : using rescaling if rescaling function exists, otherwise without 
          rescaling
*/        

decision_param DECISION_PARAM={DECISION_CASESPEC, DECISION_CASESPEC};
/*
  record containing all user specified directions for choosing
  a simulation method. Currently only one element is contained.
*/

double GENERAL_PRECISION = 5e-15;
double EIGENVALUE_EPS = 1e-15;

char ERRORSTRING_OK[MAXERRORSTRING], ERRORSTRING_WRONG[MAXERRORSTRING],
    ERROR_LOC[nErrorLoc]="";
int ERRORMODELNUMBER = -1;

char METHODNAMES[][METHODMAXCHAR]={"circulant embedding",
				   "cutoff CE",
				   "intrinsic CE",
                                   "TBM2", 
				   "TBM3",
				   "spectral TBM",
                                   "direct matrix decomposition",
				   "nugget",
                                   "add.MPP",
                                   "hyperplanes",
                                   "other methods(not implem.)",
				   "nothing", 
				   "max.MPP",
				   "forbidden"};
char DISTRNAMES[DISTRNR][DISTRMAXCHAR]={"Gauss", "Poisson", "MaxStable"};
char OP_SIGN[][3]={"+","*", "("};
int XSTARTD[MAXDIM] = {XSTARTDIM1, XSTARTDIM2, XSTARTDIM3, XSTARTDIM4},
    XENDD[MAXDIM] = {XENDDIM1, XENDDIM2, XENDDIM3, XENDDIM4},
    XSTEPD[MAXDIM]= {XSTEPDIM1, XSTEPDIM2, XSTEPDIM3, XSTEPDIM4};
int COVLISTALL[MAXCOV] = {0, 1, 2, 3, 4, 5, 6, 7, 8, MAXCOV-1};
//SEXP ERRORDUMP=NULL;
int GENERALISEDCAUCHY, CAUCHY, STABLE, WHITTLEMATERN, BROWNIAN, EXPONENTIAL;
//                                          ONLY used in Hyperfct.cc
//                                          (cutoff embedding and intrinsic)
int LocalCovList[nLocalCovList], nlocal;


void ErrorMessage(SimulationType m, int error) {
  char MS[50], EM[300+2*MAXERRORSTRING];
  switch(m) {
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
      case Nugget: strcpy(MS,"nugget effect modelling"); break;
      case AdditiveMpp: strcpy(MS,"additive MPP (random coins)"); break;
      case Hyperplane : strcpy(MS,"hyperplane tessellation"); break;
      case Special: strcpy(MS,"special procedure"); break;
      case Nothing: 
	if (error>0) strcpy(MS,"Error");
	else strcpy(MS,"Message"); 
	break;
      case MaxMpp: strcpy(MS,"max. MPP (Boolean functions)"); break;
      case Forbidden: 
	PRINTF("m=%d (forbidden call); error=%d \n", m, error);
//	strcpy(MS,"forbidden function call"); 
	assert(false);
	break;
      default : PRINTF("m=%d \n",m); assert(false);
  }
  switch (error) {
      case USEOLDINIT : strcpy(EM,"Using stored initialization"); break;
      case NOERROR : strcpy(EM,"none"); break;
      case NOERROR_REPEAT : strcpy(EM,"none; looking for further covariances applicable to the same method");break;
      case NOERROR_ENDOFLIST : strcpy(EM,"none; end of list");break;
      case ERRORDUMMY : strcpy(EM,"none"); break;
       

      case ERRORNOTDEFINED :       
	strcpy(EM,"method undefined for this model");break;
      case ERRORONLYGRIDALLOWED : 
	strcpy(EM,"only allowed for grids");break;
      case ERRORNOTPROGRAMMED :    
	strcpy(EM,"not programmed yet. Sorry"); break;
      case ERRORCOVNOTALLOWED :    
	strcpy(EM,"model not allowed for specified dimension");break;
      case ERRORFAILED: 
	strcpy(EM,"algorithm failed (partially)");break;
      case ERRORMEMORYALLOCATION: 
	strcpy(EM,"memory allocation error"); break;
      case ERRORNOTINITIALIZED: 
	strcpy(EM,"not initialized or RFparameters()$Storing==FALSE");break;
      case ERRORKEYNROUTOFRANGE: 
	strcpy(EM,"wrong key number");break;
      case ERRORDECOMPOSITION:
	strcpy(EM,"matrix decomposition failed");break;
      case ERRORPRECISION: 
	strcpy(EM,"required precision not attained");break;
      case ERRORRESCALING:
	strcpy(EM,"rescaling not defined; try PracticalRange=FALSE");break;
      case ERRORFOURIER: 
	strcpy(EM,"fft factorization failed");break;
      case ERRORCOVFAILED: 
	sprintf(EM,
		"model and method only valid for %s. Got %s",
		ERRORSTRING_OK,ERRORSTRING_WRONG);
	break;
      case ERRORCOVNUMERICAL: 
	sprintf(EM,
		"covariance exactly calculable in TBM2 only for %s. Got %s. Try TBM2.num=TRUE",
		ERRORSTRING_OK,ERRORSTRING_WRONG);
	break;
      case ERRORMAXMEMORY:
	sprintf(EM,
		"total real numbers needed=%s > allowed max real numbers=%s -- increase CE.maxmem and/or local.maxmem using RFparameters\n",
		ERRORSTRING_WRONG, ERRORSTRING_OK);
	break;
      case ERRORTOOMANYLINES:
	strcpy(EM, "estimated number of lines exceeds hyper.maxlines");break;
      case ERRORREGISTER: 
	sprintf(EM,"Register number out of [0,%d]",MAXKEYS-1);break;
      case ERRORCOORDINATES: 
	strcpy(EM,"coordinates are not given or invalid grid specification");
	break;
      case ERRORNEGATIVEVAR: 
	strcpy(EM,"Variance or nugget non-positive");break;
      case ERRORPARAMNUMBER: 
	strcpy(EM,"number of parameters not correct");break;
      case ERRORDIM: 
	sprintf(EM,"dimension specification not in [1,%d]",MAXDIM);break;
      case ERRORNEGATIVESCALE :  
	strcpy(EM,"scale parameter must be positive");break;
      case ERRORWAVING :
	strcpy(EM,"Rescaling not possible (waving)");break;
      case ERRORWRONGINIT: 
	strcpy(EM,"Wrong initialisation");break;
      case ERRORNN:
	strcpy(EM,"The number of points on the line is too large. Check your parameters and make sure that none of the locations are given twice");break;
      case ERRORANISOTROPIC:
	strcpy(EM,"anisotropic function not allowed"); break;
      case ERRORMETHODMIX:
	strcpy(EM,"for multiplicative covariance functions the method must be identical"); break;
      case ERRORTIMENOTANISO:
	strcpy(EM,"time component only allowed for anisotropic fields"); break; 
      case ERRORNOMULTIPLICATION:
	strcpy(EM,"multiplicative covariance functions/variograms not allowed");
	break;
      case ERRORANISOMIX:
	strcpy(EM,"for multiplicative covariance functions and TBM \
the anisotropy matrices must be identical"); break;
      case ERRORISOTROPICMETHOD:
	strcpy(EM,"method (TBM) only allows for spatially isotropic functions in certain dimensions"); 
	break;
      case ERRORTIMESEPARATE:
	strcpy(EM,"method (TBM) does not allow the time component being mixed up with spatial components"); 
	break;
      case ERRORUNKNOWNMETHOD:
	strcpy(EM,"Unknown method in this context or unallowed mixture of methods"); 
	break;
      case ERRORWITHOUTTIME:
	strcpy(EM,"spatially isotropic covariance fcts that are not fully isotropic must be called with a genuine time component T"); 
	break;
      case ERRORCOVNROUTOFRANGE: 
	strcpy(EM,"wrong covnr number");break;
      case ERRORWRONGDIM:
	strcpy(EM,"wrong dimension"); break;
      case ERRORNOTANISO:
	strcpy(EM,"model definition must use list-notation with anisotropy"); break;
      case ERRORDIMMISMATCH:
	strcpy(EM,"logical dimension does not match physical dimension"); break;
      case ERRORTOOMANYPOINTS:
	strcpy(EM, "too many points to determine smallest distance between the points; try TBM*.linesimustep with a positive value"); break;
      case ERRORTIMENOTALLOWED:
	strcpy(EM, "time component not allowed for the specified method"); break;
      case ERRORANYLEFT:
	strcpy(EM, "not all covariances could be simulated at the current step");
	break;
      case ERRORMETHODEXCLUDED:
	strcpy(EM, "CE Intrinsic contradictsRFparmeters()$stationary.only=true");
	break;
      case ERROROUTOFMETHODLIST:
	strcpy(EM, "run out of list of methods -- are RFparameters too\nrestrictive? Retry with RFparameters(PrintLevel=4) to get more information");
	break;
      case ERRORTBMCENTER:
	  strcpy(EM, "center contains non finite values");
	 break;  
      case ERRORTBMPOINTS:
	  strcpy(EM,
		 "given number of points to simulate on the TBM line is too small");
	 break;  
      case ERRORHYPERNR :
	  strcpy(EM, 
		 "announced number of submodels invalid (too large or non-positive)");
	  break;
      case ERRORHYPERMETHOD :
	  strcpy(EM, 
		 "methods for the hypermodel and the submodels must be given and be equal");
	  break;
      case ERRORHYPERNOTALLOWED :
	  strcpy(EM, "the method does not allow for hyper models (yet).");
	  break;
      case ERRORNCOVOUTOFRANGE:
	  strcpy(EM, "the number of covariance functions is zero or too large");
	  break;
      case ERRORPARENTHESIS:
	  strcpy(EM, "The operator '(' must be set after and only after hyper models");
	  break;
      case ERRORMSG: 
	sprintf(EM, "%s. Got %s", ERRORSTRING_OK,ERRORSTRING_WRONG);
	break;
      case ERRORTBMPOINTSZERO: 
	sprintf(EM, "if linesimufactor and linesimustep are naught then TBM_POINTS must be at least 4 (better between 100 and 10000)");
	break;

	// case ERRORISOTROPIC:
	//   strcpy(EM,"isotropic function not allowed"); break;
	//
	// extremes:
      case ERRORFULLRANK: 
	sprintf(EM, "anisotropy matrix must have full rank");
	break;
      case ERRORTIMECOMPONENT:
	sprintf(EM, "last column of the anisotropy matrix in the hypermodel must consist of zeros only");
	break;	
      case ERRORLOWRANK :
	  strcpy(EM,"confused about low rank anisotropy matrix");
	  break;
      case ERRORLOWRANKTBM :
	  strcpy(EM,"confused about low rank anisotropy matrix.\nWhen using TBM then try putting the model with anisotropy matrix of\nhighest rank at the beginning");
	  break;
      case ERRORHYPERNOTISO :
	  strcpy(EM, 
		 "only isotropic models allowed for submodels of hyper models");
	  break;
      case ERRORLAYERS :
	  strcpy(EM, "value of TBM*.layers does not match the situation");
	  break;


      case ERRORSILLNULL : 
	strcpy(EM,"Vanishing sill not allowed");break;     
      case ERRORPAIRS : 
	strcpy(EM,"No pair simulation (currently) allowed.");break;
      case ERRORTREND : 
	strcpy(EM,"No trend  (currently) allowed");break;
	//    case : strcpy(EM,"");break;
	//
	// Poisson:
      case ERRORVARMEAN :
	strcpy(EM,"Variance differs from mean");break;
	
      case MSGLOCAL_OK :
	strcpy(EM,"fine");break;
      case MSGLOCAL_JUSTTRY :
	strcpy(EM,
	       "unclear whether algorithm will work for specified parameters");
	break;
      case MSGLOCAL_NUMOK :
	strcpy(EM,"fine. Algorithm should work for specified parameters");
	break;
      case MSGLOCAL_ENDOFLIST :
	strcpy(EM,"end of list for variants of the algorithm");break;
      case MSGLOCAL_SIGNPHI :
	strcpy(EM,"wrong sign of covariance function at break point");break;
      case MSGLOCAL_SIGNPHIFST :
	strcpy(EM,
	       "wrong sign of 1st derivative of the covariance function at the break point");
	break;
      case MSGLOCAL_SIGNPHISND :
	strcpy(EM,
	       "wrong sign of 2nd derivative of the covariance function at the break point");
	break;
      case MSGLOCAL_INITINTRINSIC :
	strcpy(EM,"one of a2, b or a0+phi(0) has wrong sign");break;


     case ERRORUNSPECIFIED :
	strcpy(EM,"(unspecified)");break;
     default : 
	PRINTF(" m=%d error=%d\n",m,error); 
	assert(false);
  }
  if (strcmp(ERROR_LOC,"") != 0) PRINTF("In %s", ERROR_LOC);
  if (m!=Nothing)  PRINTF("Method ");
  PRINTF("%s", MS);
  if (ERRORMODELNUMBER >= 0) PRINTF(" in submodel #%d", ERRORMODELNUMBER + 1);
  PRINTF(": %s.\n", EM);
  ERRORMODELNUMBER = -1;
}

int init_nothing(key_type *key, int m) {
  return ERRORFAILED;
}

void InitModelList()
{
  int i, nr;

  for (i=0; i<MAXDIM; i++) { ZERO[i]=0.0; UNIT[i]=1.0; }

  // init simulation methods;
  init_method[CircEmbed] = init_circ_embed;
  init_method[CircEmbedCutoff] = init_circ_embed_local; 
  init_method[CircEmbedIntrinsic] = init_circ_embed_local;
  init_method[TBM2] = init_turningbands2;
  init_method[TBM3] = init_turningbands3;
  init_method[SpectralTBM] = init_simulatespectral;
  init_method[Direct] = init_directGauss;
  init_method[Nugget] = init_nugget;
  init_method[AdditiveMpp] = init_mpp;
  init_method[Hyperplane] = init_hyperplane;
  init_method[Special] = init_special;
  init_method[Nothing] = init_nothing;

  do_compatiblemethod[CircEmbed] = do_circ_embed;
  do_compatiblemethod[CircEmbedCutoff] = do_circ_embed_local; 
  do_compatiblemethod[CircEmbedIntrinsic] = do_circ_embed_local; 
  do_compatiblemethod[TBM2] = NULL;
  do_compatiblemethod[TBM3] = NULL;
  do_compatiblemethod[SpectralTBM] = NULL;
  do_compatiblemethod[Direct] = do_directGauss;
  do_compatiblemethod[Nugget] = do_nugget;
  do_compatiblemethod[AdditiveMpp] = NULL;
  do_compatiblemethod[Hyperplane] = NULL;
  do_compatiblemethod[Special] = NULL;

  do_incompatiblemethod[CircEmbed] = NULL;
  do_incompatiblemethod[CircEmbedCutoff] = NULL;     
  do_incompatiblemethod[CircEmbedIntrinsic] = NULL; 
  do_incompatiblemethod[TBM2] = do_turningbands;
  do_incompatiblemethod[TBM3] = do_turningbands;
  do_incompatiblemethod[SpectralTBM] = do_simulatespectral;
  do_incompatiblemethod[Direct] = NULL;
  do_incompatiblemethod[Nugget] = NULL;
  do_incompatiblemethod[AdditiveMpp] = do_addmpp;
  do_incompatiblemethod[Hyperplane] = do_hyperplane;
  do_incompatiblemethod[Special] = do_special;

   // init models
  for (i=0;i<MAXKEYS;i++) KEY_NULL(&(KEY[i]));

  if (CovList!=NULL) {
    PRINTF("List of covariance functions looks already initiated.\n"); 
    return;
  }
  assert(currentNrCov=-1);
  CovList = (cov_fct*) malloc(sizeof(cov_fct) * (MAXNRCOVFCTS+1));
  // + 1 is necessary because of COVINFO_NULL that uses the last + 
  currentNrCov=0;
  nlocal=0;
  
  nr=IncludeModel("bessel",1, checkBessel, ISOTROPIC, false,
		  infoBessel, rangeBessel);
  addCov(nr, Bessel, NULL, NULL);
  addTBM(nr, NULL, NULL, spectralBessel);
 	       
  nr=IncludeModel("cauchy", 1, checkCauchy, ISOTROPIC, false,
		  infoCauchy, rangeCauchy);
  addCov(nr, Cauchy, DCauchy, ScaleCauchy);
  addTBM(nr, TBM2Cauchy, TBM3Cauchy, NULL);
  addLocal(nr, true, DDCauchy, &CAUCHY);
	       
  nr=IncludeModel("cauchytbm", 3, checkOK, ISOTROPIC, false,
		  infoCauchytbm, rangeCauchytbm);
  addCov(nr,Cauchytbm, DCauchytbm, NULL);
  addTBM(nr,NULL, TBM3Cauchytbm, NULL);
	      
  nr=IncludeModel("circular", 0, checkOK,ISOTROPIC, false, 
		  infocircular, rangecircular);
  addCov(nr,circular, Dcircular, Scalecircular);
  addTBM(nr, NULL, NULL, NULL);
  addOther(nr, circular_init, circularMpp, NULL, NULL, NULL);

  nr=IncludeModel("cone",3,checkcone, ISOTROPIC, false,
		  infocone, rangecone);
  addOther(nr, cone_init, cone,  NULL, NULL, NULL);

  nr=IncludeModel("constant",0, checkOK, ISOTROPIC, false, 
		  infoconstant, rangeconstant);
  addCov(nr, constant, Dconstant, NULL);
  addTBM(nr,TBM2constant,TBM3constant, NULL);
   
  nr=IncludeModel("cubic",0, checkOK,ISOTROPIC, false, infocubic, 
		  rangecubic);
  addCov(nr, cubic, Dcubic, Scalecubic);
  addTBM(nr,NULL, TBM3cubic, NULL);
	       
  nr=IncludeHyperModel("cutoff", 3, checkNinit_co, ISOHYPERMODEL, 
		       false, info_co, range_co);
  addCov(nr, co, NULL, NULL);
  addInitLocal(nr, getcutoffparam_co, alternativeparam_co, CircEmbedCutoff);

  nr=IncludeModel("dampedcosine", 1, checkOK, ISOTROPIC, false,
		  infodampedcosine, rangedampedcosine);
  addCov(nr,dampedcosine,Ddampedcosine,Scaledampedcosine);
  addTBM(nr,NULL,TBM3dampedcosine, NULL);

  nr=IncludeModel("exponential",0, checkexponential,ISOTROPIC, false, 
		  infoexponential, rangeexponential);
  addCov(nr,exponential,Dexponential,Scaleexponential);
  //addTBM(nr,NULL,TBM3exponential,Dexponential,NULL);
  //addTBM(nr,NULL,TBM3exponential,Dexponential,spectralexponential);
  addTBM(nr, TBM2exponential, TBM3exponential, spectralexponential);
  addLocal(nr, true, DDexponential, &EXPONENTIAL);
  addOther(nr, NULL, NULL, hyperexponential, NULL, NULL);
  // numerisches Ergebnis passt nicht !
//  printf("%d %d\n", nr, CovList[nr].implemented[CircEmbedCutoff]); assert(false);

  nr=IncludeModel("FD", 1, checkOK, ISOTROPIC, false, infoFD, rangeFD);
  addCov(nr, FD, NULL, NULL);

  nr=IncludeModel("fractalB",1,checkOK,ISOTROPIC, true,
		 infofractalBrownian, rangefractalBrownian);
  addCov(nr, fractalBrownian, DfractalBrownian, NULL);
  addLocal(nr, false, DDfractalBrownian, &BROWNIAN);
  
  nr=IncludeModel("fractgauss", 1, checkOK, ISOTROPIC, false,
		  infofractGauss, rangefractGauss);
  addCov(nr, fractGauss, NULL, NULL);

  nr=IncludeModel("gauss", 0, checkOK, ISOTROPIC, false, infoGauss, 
		  rangeGauss);
  addCov(nr,Gauss, DGauss, ScaleGauss);
  addTBM(nr, NULL, TBM3Gauss, spectralGauss);
  addOther(nr, gaussmpp_init, gaussmpp, NULL,  NULL, NULL);

  nr=IncludeModel("gencauchy", 2, checkgeneralisedCauchy, ISOTROPIC, false,
		  infogeneralisedCauchy, rangegeneralisedCauchy);
  addCov(nr,generalisedCauchy, DgeneralisedCauchy, ScalegeneralisedCauchy);
  addTBM(nr, NULL, TBM3generalisedCauchy, NULL);
  addLocal(nr, true, DDgeneralisedCauchy, &GENERALISEDCAUCHY);

  nr=IncludeModel("gengneiting", 2, checkOK, ISOTROPIC, false,
		  infogenGneiting, rangegenGneiting);
  addCov(nr, genGneiting, DgenGneiting, NULL);
  addTBM(nr, NULL, TBM3genGneiting, NULL);

  nr=IncludeModel("gneiting", 0, checkOK, ISOTROPIC, false, 
		  infoGneiting, rangeGneiting); 
  addCov(nr,Gneiting,DGneiting,ScaleGneiting);
  addTBM(nr, NULL, TBM3Gneiting, NULL);

  nr=IncludeModel("hyperbolic",3,checkhyperbolic, ISOTROPIC, false,
		  infohyperbolic, rangehyperbolic);
  addCov(nr,hyperbolic,Dhyperbolic,NULL);
  addTBM(nr,NULL,TBM3hyperbolic, NULL);

  nr=IncludeModel("iacocesare", 3, checkIacoCesare, SPACEISOTROPIC, false,
		  infoIacoCesare, rangeIacoCesare);
  addCov(nr, IacoCesare, NULL, NULL);
  

  nr=IncludeModel("lgd1", 2, checkOK, ISOTROPIC, false,
		 infolgd1, rangelgd1);
  addCov(nr, lgd1, Dlgd1, Scalelgd1);
  addTBM(nr, NULL, NULL, NULL);

  nr=IncludeHyperModel("mastein", 3, checkNinit_MaStein, ANISOHYPERMODEL, false, 
		       info_MaStein, range_MaStein);
  addCov(nr, MaStein, NULL, NULL);

  nr=IncludeModel("nsst", 6, checkspacetime1, SPACEISOTROPIC, false, 
		  infospacetime1, rangespacetime1);
  addCov(nr,spacetime1,Dspacetime1,NULL);
  addTBM(nr, NULL, //TBM2spacetime1,
	 NULL, NULL);
  
  nr=IncludeModel("nsst2",7,checkspacetime2,SPACEISOTROPIC, false, 
		  infospacetime2, rangespacetime2);
  addCov(nr,spacetime2,Dspacetime2,NULL);
  addTBM(nr,NULL, NULL /* TBM3spacetime2 */, NULL);

  nr=IncludeModel("nugget",0, checknugget, ISOTROPIC, false, infonugget, 
		  rangenugget);
  addCov(nr, nugget, NULL, Scalenugget);
  modelexception(nr, true, false);
  
  /* cf. convert.R, PrepareModel, near end of function */		  
  for (i=0; i<SimulationTypeLength; i++) {
    CovList[nr].implemented[i] = GIVEN_METH_IGNORED;
  }
  CovList[nr].implemented[Nugget] = CovList[nr].implemented[Direct] =
    CovList[nr].implemented[CircEmbed] =IMPLEMENTED;
 
  nr=IncludeModel("penta",0, checkOK,
                  ISOTROPIC, false, infopenta, rangepenta);
  addCov(nr,penta,Dpenta,Scalepenta);
  addTBM(nr,NULL,TBM3penta,NULL);
	
  nr=IncludeModel("power",1,checkpower,ISOTROPIC, false, 
		  infopower, rangepower);
  addCov(nr,power,Dpower,Scalepower);
  addTBM(nr, NULL, TBM3power, NULL);
 
  nr=IncludeModel("qexponential", 1, checkOK, ISOTROPIC, false,
		  infoqexponential, rangeqexponential);
  addCov(nr,qexponential, Dqexponential, Scaleqexponential);
  addTBM(nr, NULL, TBM3qexponential, NULL);

  nr=IncludeModel("spherical",0, checkOK, ISOTROPIC, false, 
		  infospherical, rangespherical);
  addCov(nr,spherical,Dspherical,Scalespherical);
  addTBM(nr, TBM2spherical, TBM3spherical,NULL);
  addOther(nr, spherical_init, sphericalMpp, NULL, NULL, NULL);

  nr=IncludeModel("stable", 1, checkstable, ISOTROPIC, false, 
		  infostable, rangestable);
  addCov(nr, stable, Dstable,  Scalestable);
  addTBM(nr, NULL, TBM3stable, NULL);
  addLocal(nr, true, DDstable, &STABLE);
 

  // SPACEISOTROPIC variant of stable -- used for testing purposes only
  nr=IncludeModel("stableX", 1, checkOK, SPACEISOTROPIC, false, 
		  infostable, rangestable);
  addCov(nr, stableX, DstableX, Scalestable);
  addTBM(nr, NULL, NULL, NULL);


  nr=IncludeHyperModel("Stein", 3, checkNinit_Stein, ISOHYPERMODEL, false, 
		       info_Stein, range_Stein);
  addCov(nr, Stein, NULL, NULL);
  addInitLocal(nr, getintrinsicparam_Stein, alternativeparam_Stein,
	       CircEmbedIntrinsic);

  nr=IncludeModel("steinst1", kappasSteinST1, checkSteinST1, ANISOTROPIC,
		  false, infoSteinST1, rangeSteinST1);
  addCov(nr, SteinST1, NULL, NULL);

  nr=IncludeModel("wave", 0, checkOK, ISOTROPIC, false, infowave, 
		  rangewave);
  addCov(nr,wave,NULL,Scalewave);
  addTBM(nr,NULL,NULL,spectralwave);
 
  nr=IncludeModel("whittlematern", 1, checkWhittleMatern, ISOTROPIC, false,
		  infoWhittleMatern, rangeWhittleMatern);
  addCov(nr,WhittleMatern, DWhittleMatern, ScaleWhittleMatern);
  addTBM(nr,NULL,TBM3WhittleMatern, spectralWhittleMatern);
  addLocal(nr, true, DDWhittleMatern, &WHITTLEMATERN);
  
  /* 
     in case of anisotropic models: do not forget to set `addodd',
     i.e. the coordinate dimensions where the covariance function is not
     an even function!
  */

  // must be the very last one!
  nr = IncludeModel("undefined", 0, checkundefined, ISOTROPIC, false, 
		    infoundefined, rangenugget);
  assert(nr > 0); // otherwise we have already reached the maximum number
  //                of models; note that it is planned that further
  //                functions might be added by users
  currentNrCov--;
  for (++nr; nr <= MAXNRCOVFCTS; nr++) 
    memcpy(&(CovList[nr]), &(CovList[currentNrCov]), sizeof(cov_fct));
  
  addusersfunctions();

}

