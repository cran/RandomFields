
/* 
 Authors
 Martin Schlather, schlath@hsu-hh.de

 library for simulation of random fields -- init part and error messages

 Copyright (C) 2001 -- 2004 Martin Schlather

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

#include <unistd.h>

initfct init_method[SimulationTypeLength];
do_comp_fct do_compatiblemethod[SimulationTypeLength];
do_incomp_fct  do_incompatiblemethod[SimulationTypeLength];
key_type KEY[MAXKEYS]; //static key_type KEY[MAXKEYS];
double ZERO[MAXDIM], UNIT[MAXDIM];
cov_fct *CovList=NULL;
int currentNrCov=-1;


char GENERAL_PCH[2]="*"; 
/*  character printed after each simulation
    just for entertainment of the user
    except for "!", then the numbers are shown
*/
int GENERAL_STORING=true; 
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

char ERRORSTRING_OK[MAXERRORSTRING], ERRORSTRING_WRONG[MAXERRORSTRING];

char METHODNAMES[][METHODMAXCHAR]={"circulant embedding",
				   "cutoff CE",
				   "intrinsic CE",
                                   "TBM2 (2dim. turning bands)", 
				   "TBM3",
				   "spectral TBM (2dim)",
                                   "direct matrix decomposition",
				   "nugget",
                                   "add.MPP",
                                   "hyperplanes",
                                   "other methods(not implem.)",
				   "nothing", 
				   "max.MPP",
				   "forbidden"};
char DISTRNAMES[DISTRNR][DISTRMAXCHAR]={"Gauss", "Poisson", "MaxStable"};
char OP_SIGN[][2]={"+","*"};
//SEXP ERRORDUMP=NULL;

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
      case Direct : strcpy(MS,"direct Gaussian (decomposition of cov. matrix)");
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
      case NOERROR : strcpy(EM,"fine"); break;
      case NOERROR_REPEAT : strcpy(EM,"none; repeat");break;
      case NOERROR_ENDOFLIST : strcpy(EM,"end of list");break;
      case ERRORDUMMY : strcpy(EM,"none"); break;
       

      case ERRORNOTDEFINED :       
	strcpy(EM,"undefined for this model");break;
      case ERRORMETHODNOTALLOWED : 
	strcpy(EM,"not allowed for the specified parameters or points");break;
      case ERRORNOTPROGRAMMED :    
	strcpy(EM,"not programmed yet. Sorry"); break;
      case ERRORCOVNOTALLOWED :    
	strcpy(EM,"model not allowed for specified dimension");break;
      case ERRORFAILED: 
	strcpy(EM,"algorithm failed (partially)");break;
      case ERRORMEMORYALLOCATION: 
	strcpy(EM,"memory allocation error");break;
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
		"cov. of model & method calculable exactly only for %s. Got %s. Try TBM2.num=TRUE.",
		ERRORSTRING_OK,ERRORSTRING_WRONG);
	break;
      case ERRORMAXMEMORY:
	sprintf(EM,
		"total real numbers needed=%s > allowed max real numbers=%s -- increase CE.maxmem using RFparameters\n",
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
the anisotropies must be identical"); break;
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
	strcpy(EM,"spatially isotropic covariance fcts that are not fully isotropic must be called with a time component"); 
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
	strcpy(EM, "two many points to determine smallest distance between the points"); break;
      case ERRORTIMENOTALLOWED:
	strcpy(EM, "time component not allowed for the specified method"); break;
      case ERRORANYLEFT:
	strcpy(EM, "not all covariances could be simulated at the current step");
	break;
      case ERRORMETHODEXCLUDED:
	strcpy(EM, "CE Intrinsic contradictsRFparmeters()$stationary.only=true");
	break;
      case ERROROUTOFMETHODLIST:
	strcpy(EM, "run out of list of methods -- are RFparameters too restrictive?\nRetry with RFparameters(PrintLevel=4) to learn more about the situation");
	break;
	
	// case ERRORISOTROPIC:
	//   strcpy(EM,"isotropic function not allowed"); break;
	//
	// extremes:
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
	
      case ERRORTRIVIAL :
	strcpy(EM,"confused about trivial parameter specification");break;
      case ERRORUNSPECIFIED :
	strcpy(EM,"(unspecified)");break;
     default : 
	PRINTF(" m=%d error=%d\n",m,error); 
	assert(false);
  }
  if (m!=Nothing)  PRINTF("Method");
  PRINTF(" %s: %s.\n",MS,EM);
}

int init_nothing(key_type *key, int m) {
  return ERRORFAILED;
}

void InitModelList()
{
  int i,d,v, nr;

  for (i=0; i<MAXDIM; i++) { ZERO[i]=0.0; UNIT[i]=1.0; }

  // init simulation methods;
  init_method[CircEmbed] = init_circ_embed;
  init_method[CircEmbedCutoff] = init_circ_embed_cutoff; 
  init_method[CircEmbedIntrinsic] = init_circ_embed_intrinsic;
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
  do_compatiblemethod[CircEmbedCutoff] = do_circ_embed; 
  do_compatiblemethod[CircEmbedIntrinsic] = do_circ_embed_intrinsic; 
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
  for (i=0;i<MAXKEYS;i++) {
    KEY[i].active = false;
    KEY[i].n_unimeth = 0;
    for (v=0; v<MAXCOV; v++) {
      KEY[i].S[v]=NULL;
      KEY[i].destruct[v]=NULL;
      //KEY[i].cov[v] = NULL;
    }
    KEY[i].destructExtremes=NULL;
    KEY[i].SExtremes=NULL;
    for (d=0; d<MAXDIM; d++) KEY[i].x[d]=NULL;
    //    KEY[i].naturalscaling =-1;
    KEY[i].TrendModus = -1;
    KEY[i].TrendFunction = NULL;
    KEY[i].LinearTrend = NULL;
  }

  if (CovList!=NULL) {
    PRINTF("List of covariance functions looks already initiated.\n"); 
    return;
  }
  assert(currentNrCov=-1);
  CovList = (cov_fct*) malloc(sizeof(cov_fct) * MAXNRCOVFCTS);
  currentNrCov=0;
  
  nr=IncludeModel("bessel",1, checkBessel, FULLISOTROPIC, false,
		  infoBessel, rangeBessel);
  addCov(nr, Bessel, NULL, NULL);
  addTBM(nr, NULL, NULL, Nothing, spectralBessel);
 	       
  nr=IncludeModel("cauchy", 1, checkCauchy, FULLISOTROPIC, false,
		  infoCauchy, rangeCauchy);
  addCov(nr, Cauchy, DCauchy, ScaleCauchy);
  addTBM(nr, TBM2Cauchy, TBM3Cauchy, CircEmbed,NULL);
	       
  nr=IncludeModel("cauchytbm",3,checkCauchytbm,FULLISOTROPIC, false,
		  infoCauchytbm, rangeCauchytbm);
  addCov(nr,Cauchytbm, DCauchytbm, NULL);
  addTBM(nr,NULL, TBM3Cauchytbm, CircEmbed,NULL);
	      
  nr=IncludeModel("circular", 0, NULL,FULLISOTROPIC, false, 
		  infocircular, rangecircular);
  addCov(nr,circular, Dcircular, Scalecircular);
  addTBM(nr, NULL, NULL, CircEmbed, NULL);
  addOther(nr, circular_init, circularMpp, NULL, NULL, NULL);

  nr=IncludeModel("cone",3,checkcone,FULLISOTROPIC, false,
		  infocone, rangecone);
  addOther(nr, cone_init, cone,  NULL, NULL, NULL);

  nr=IncludeModel("cubic",0,NULL,FULLISOTROPIC, false, infocubic, 
		  rangecubic);
  addCov(nr, cubic, Dcubic, Scalecubic);
  addTBM(nr,NULL, TBM3cubic, CircEmbed, NULL);
	       
  nr=IncludeModel("dampedcosine",1,checkdampedcosine,FULLISOTROPIC, false,
		  infodampedcosine, rangedampedcosine);
  addCov(nr,dampedcosine,Ddampedcosine,Scaledampedcosine);
  addTBM(nr,NULL,TBM3dampedcosine,CircEmbed,NULL);

  nr=IncludeModel("exponential",0,NULL,FULLISOTROPIC, false, 
		  infoexponential, rangeexponential);
  addCov(nr,exponential,Dexponential,Scaleexponential);
  //addTBM(nr,NULL,TBM3exponential,Dexponential,CircEmbed,NULL);
  //addTBM(nr,NULL,TBM3exponential,Dexponential,CircEmbed,spectralexponential);
  addTBM(nr,TBM2exponential,TBM3exponential,CircEmbed, spectralexponential);
  addOther(nr, NULL, NULL, hyperexponential, NULL, NULL);

  nr=IncludeModel("FD", 1, checkFD, FULLISOTROPIC, false, infoFD, rangeFD);
  addCov(nr, FD, NULL, NULL);

  nr=IncludeModel("fractalB",1,checkfractalBrownian,FULLISOTROPIC, true,
		 infofractalBrownian, rangefractalBrownian);
  addCov(nr, fractalBrownian, DfractalBrownian, NULL);
  addLocal(nr,fractalBrownian_intrinsic_strategy, NULL, DDfractalBrownian);
 
  nr=IncludeModel("fractgauss", 1, checkfractGauss,FULLISOTROPIC, false,
		  infofractGauss, rangefractGauss);
  addCov(nr, fractGauss, NULL, NULL);

  nr=IncludeModel("gauss", 0, NULL, FULLISOTROPIC, false, infoGauss, 
		  rangeGauss);
  addCov(nr,Gauss, DGauss, ScaleGauss);
  addTBM(nr,NULL, TBM3Gauss, CircEmbed, spectralGauss);
  addOther(nr, gaussmpp_init, gaussmpp, NULL,  NULL, NULL);

  nr=IncludeModel("gencauchy", 2, checkgeneralisedCauchy,FULLISOTROPIC, false,
		  infogeneralisedCauchy, rangegeneralisedCauchy);
  addCov(nr,generalisedCauchy, DgeneralisedCauchy, ScalegeneralisedCauchy);
  addTBM(nr, NULL, TBM3generalisedCauchy, CircEmbed, NULL);
  addLocal(nr,generalisedCauchy_intrinsic_strategy,
	   generalisedCauchy_cutoff_strategy, DDgeneralisedCauchy);

  nr=IncludeModel("gengneiting",2,checkgenGneiting,FULLISOTROPIC, false,
		  infogenGneiting, rangegenGneiting);
  addCov(nr,genGneiting,DgenGneiting,NULL);
  addTBM(nr,NULL,TBM3genGneiting, CircEmbed,NULL);

  nr=IncludeModel("gneiting",0,NULL,FULLISOTROPIC, false, infoGneiting, 
		  rangeGneiting); 
  addCov(nr,Gneiting,DGneiting,ScaleGneiting);
  addTBM(nr,NULL,TBM3Gneiting,CircEmbed,NULL);

  nr=IncludeModel("hyperbolic",3,checkhyperbolic,FULLISOTROPIC, false,
		  infohyperbolic, rangehyperbolic);
  addCov(nr,hyperbolic,Dhyperbolic,NULL);
  addTBM(nr,NULL,TBM3hyperbolic,CircEmbed,NULL);

  nr=IncludeModel("lgd1", 2, checklgd1, FULLISOTROPIC, false,
		 infolgd1, rangelgd1);
  addCov(nr, lgd1, Dlgd1, Scalelgd1);
  addTBM(nr,NULL, NULL, CircEmbed, NULL);

  nr=IncludeModel("nsst", 6, checkspacetime1, SPACEISOTROPIC, false, 
		  infospacetime1, rangespacetime1);
  addCov(nr,spacetime1,Dspacetime1,NULL);
  addTBM(nr,TBM2spacetime1,TBM3spacetime1,CircEmbed,NULL);
  
  nr=IncludeModel("nsst2",7,checkspacetime2,SPACEISOTROPIC, false, 
		  infospacetime2, rangespacetime2);
  addCov(nr,spacetime2,Dspacetime2,NULL);
  addTBM(nr,NULL,TBM3spacetime2,CircEmbed,NULL);

//  nr=IncludeModel("nsst3",6,checkspacetime1,SPACEISOTROPIC, false, 
//		  infospacetime1, rangespacetime1);
//  addCov(nr,spacetime3,NULL);

  nr=IncludeModel("nugget",0,NULL,FULLISOTROPIC, false, infonugget, rangenugget);
  addCov(nr, nugget, NULL, Scalenugget);
  /* cf. convert.R, PrepareModel, near end of function */		  
  for (i=0; i<SimulationTypeLength; i++)
    CovList[nr].implemented[i] = GIVEN_METH_IGNORED;
  CovList[nr].implemented[Nugget] = CovList[nr].implemented[Direct] =
    CovList[nr].implemented[CircEmbed] =IMPLEMENTED;
 
  nr=IncludeModel("penta",0,NULL,FULLISOTROPIC, false, infopenta, rangepenta);
  addCov(nr,penta,Dpenta,Scalepenta);
  addTBM(nr,NULL,TBM3penta,CircEmbed,NULL);
	
  nr=IncludeModel("power",1,checkpower,FULLISOTROPIC, false, 
		  infopower, rangepower);
  addCov(nr,power,Dpower,Scalepower);
  addTBM(nr, NULL, TBM3power,CircEmbed,NULL);
 
  nr=IncludeModel("qexponential",1,checkqexponential,FULLISOTROPIC, false,
		  infoqexponential, rangeqexponential);
  addCov(nr,qexponential, Dqexponential, Scaleqexponential);
  addTBM(nr, NULL, TBM3Dexponential, CircEmbed,NULL);

  nr=IncludeModel("spherical",0, NULL, FULLISOTROPIC, false, 
		  infospherical, rangespherical);
  addCov(nr,spherical,Dspherical,Scalespherical);
  addTBM(nr,TBM2spherical,TBM3spherical,CircEmbed,NULL);
  addOther(nr, spherical_init, sphericalMpp, NULL, NULL, NULL);

  nr=IncludeModel("stable",1,checkstable,FULLISOTROPIC, false, 
		  infostable, rangestable);
  addCov(nr,stable,Dstable, Scalestable);
  addTBM(nr,NULL,TBM3stable,CircEmbed,NULL);
  addLocal(nr, stable_intrinsic_strategy, stable_cutoff_strategy, DDstable);
 
  nr=IncludeModel("wave", 0, NULL, FULLISOTROPIC, false, infowave, rangewave);
  addCov(nr,wave,NULL,Scalewave);
  addTBM(nr,NULL,NULL,Nothing,spectralwave);
 
  nr=IncludeModel("whittlematern",1,checkWhittleMatern,FULLISOTROPIC, false,
		  infoWhittleMatern, rangeWhittleMatern);
  addCov(nr,WhittleMatern,DWhittleMatern, ScaleWhittleMatern);
  addTBM(nr,NULL,TBM3WhittleMatern,CircEmbed,spectralWhittleMatern);
  addLocal(nr, WhittleMatern_intrinsic_strategy, WhittleMatern_cutoff_strategy,
	   DDWhittleMatern);


//  nr = IncludeModel("test", 0, NULL, FULLISOTROPIC, false, infotest, NULL);
//  addCov(nr, testCov, NULL);
 
  /* 
     in case of anisotropic models: do not forget to set `addodd',
     i.e. the coordinate dimensions where the covariance function is not
     an even function!
  */

}

