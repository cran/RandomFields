/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part 

 Copyright (C) 2017 -- 2018 Martin Schlather

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



#include <stdio.h>   
#include <string.h>
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "RF.h"
#include "startGetNset.h"
#include "Processes.h"
#include "shape.h"
#include "rf_interfaces.h"
#include "primitive.others.h"
#include "operator.h"
extern const char * distr[distrN];
extern const char *CE[CEN];



int  BALL,
  DISTRIBUTION, DETERM_DISTR, GAUSS_DISTR, SETPARAM, COVFCTN, ANGLE,
  COVMATRIX, RFGET, STROKORB_MONO, STROKORB_BALL_INNER, RECTANGULAR,
  CONST, BIND,
  POLYGON, SCATTER, MCMC_PGS, MCMC,LIKELIHOOD_CALL, LINEARPART_CALL,
  PREDICT_CALL,
  IDCOORD, MULT_INVERSE,
  TRUNCSUPPORT, SHAPESTP, SHAPEAVE, BROWNRESNICK, UNIF, MPPPLUS, 
  BRSHIFTED_USER, BRMIXED_USER, BRORIGINAL_USER, BRNORMED,
  BRSHIFTED_INTERN, BRMIXED_INTERN, BRORIGINAL_INTERN,   
  EXTREMALGAUSSIAN, RANDOMSIGN,  
  ARCSQRT_DISTR, SHAPEPOW,
  ZHOU, BALLANI, STATIONARY_SHAPE, STANDARD_SHAPE,
  LOC, SET_DISTR, SCALESPHERICAL, TREND, // TREND_PROC,
  COVARIATE, TRAFO, TRAFOPROC, PROJMODEL,
  VARIOGRAM_CALL, SIMULATE,
  MISSING_COV, NULL_MODEL,
  DOLLAR_PROC,R, PLUS_PROC, M_PROC,
  MPPPLUS_PROC, MULT_PROC, 
  BINARYPROC, BROWNRESNICKPROC,
  GAUSSPROC, POISSONPROC,  SCHLATHERPROC, SMITHPROC, CHI2PROC,
  EXTREMALTPROC, TPROC, TREND_PROC, PROD_PROC, VAR2COV_PROC,
  NUGGET_USER, NUGGET_INTERN,
  CIRCEMBED,  SPECTRAL_PROC_USER, SPECTRAL_PROC_INTERN,
  DIRECT, SEQUENTIAL, SPECIFIC, SELECTNR,
  AVERAGE_USER, AVERAGE_INTERN, HYPERPLANE_USER, HYPERPLANE_INTERN,
  RANDOMCOIN_USER, CE_CUTOFFPROC_USER, CE_CUTOFFPROC_INTERN, 
  CE_INTRINPROC_USER, CE_INTRINPROC_INTERN, TBM_PROC_USER, TBM_PROC_INTERN,
  SCALEMODEL, BUBBLE
  ;


 

void includeOtherModels() {  
  MISSING_COV =
    IncludePrim("missing", OtherType, 0, XONLY, CARTESIAN_COORD,
		checkMissing,  NULL, INFDIM, (ext_bool) true, NOT_MONOTONE);
  make_internal(); 

  NULL_MODEL =
    IncludePrim("null", ManifoldType, 1, XONLY, ISOTROPIC,
		checkOK, rangeNullModel, INFDIM, (ext_bool) true, NOT_MONOTONE);
  kappanames("type", INTSXP);
  addCov(0, NullModel, NullModel, NullModel, NullModel, NULL);
  RandomShape(INFTY, structOK, initOK, doOK, do_random_ok, false, false, false);
  make_internal(); 
  addTypeFct(TypeNullModel);

 // *******************
  // **** trend-models ****
  // ******************
  pref_type ptrend = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5};
        //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  TREND = IncludeModel("trend", TrendType,  0, 0, 1,
		       kappatrend, 
		       XONLY, PARAMDEP_I,
		       checktrend, 
		       rangetrend,
		       ptrend,
		       false, PARAM_DEP, INFDIM, (ext_bool) false,
		       NOT_MONOTONE);
  //  kappanames("mean", REALSXP, "plane", REALSXP, "polydeg",		
  //             INTSXP, "polycoeff",					
  //REALSXP, "arbitrark2xyfct", CLOSXP, "fctcoeff", REALSXP);	
  kappanames("mean", REALSXP);
  change_sortof(TREND_MEAN, TRENDPARAM);
  change_typeof(TREND_MEAN, ShapeType);
  addCov(trend, NULL, NULL);
  setDI(NULL, allowedItrend, settrend);
  addTypeFct(Typetrend);
  //  addCov(trend_nonstat); 


 // *******************
  // **** definite functions  ****
  // *******************

  SELECTNR =  // to do: replace by parameter in '+', selecting the 'type' or
    // 'models'
    IncludeModel("select", TcfType, 1, MAXSUB, 1, NULL, SUBMODEL_D, SUBMODEL_I,
		 checkselect, rangeselect, PREF_ALL,
		 true, PARAM_DEP, INFDIM, (ext_bool)SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("subnr", INTSXP);
  addCov(select, NULL, NULL);
  setDI(allowedDselect, allowedIselect, NULL);
  
  addReturns(covmatrix_select, iscovmatrix_select);


  // *******************
  // **** shape types  ****
  // *******************
 ANGLE = 
   IncludeModel("angle", ShapeType, 0, 0, 4, kappa_Angle, XONLY,
		CARTESIAN_COORD, checkAngle, rangeAngle, PREF_NOTHING, 
	       false, PARAM_DEP, INFDIM, (ext_bool) false, NOT_MONOTONE);
  nickname("angle");
  addCov(Angle, NULL, NULL, invAngle, NULL);
  kappanames("angle", REALSXP, "lat.angle", REALSXP, "ratio", 
	     REALSXP, "diag", REALSXP);
 

  BALL= IncludePrim("ball",  ShapeType, 0,  NULL, 
		    XONLY, ISOTROPIC, checkOK, NULL, PREF_NOTHING, 
		    SCALAR, INFDIM-1, (ext_bool) true, MONOTONE);
  addCov(ball, NULL, NULL, Inverseball);
  RandomShape(INFTY, struct_ball, init_ball, do_ball);   
  Taylor(1.0, 0.0);

 
  COVARIATE = // intern ok
    IncludeModel("covariate", ShapeType, 0, 1, 5, kappa_covariate,
		 XONLY, UNREDUCED, // zwingend bei RAW-Konstruktionen !!
		 checkcovariate, rangecovariate, PREF_NOTHING, 
		 INTERN_SHOW, PARAM_DEP, INFDIM-1, (ext_bool) false,
		 NOT_MONOTONE);
  subnames("norm");
  kappanames(COVARIATE_C_NAME, LISTOF + REALSXP,
	     COVARIATE_X_NAME, VECSXP, 
	     COVARIATE_RAW_NAME, INTSXP,
	     COVARIATE_ADDNA_NAME, INTSXP,
	     "factor", REALSXP);
  change_sortof(COVARIATE_X, DONOTVERIFYPARAM);
  change_sortof(COVARIATE_ADDNA, IGNOREPARAM);
  change_sortof(COVARIATE_FACTOR, TRENDPARAM); 
  addCov(covariate, NULL, NULL);
  setptwise(pt_paramdep);
  AddVariant(ShapeType, ISOTROPIC); // only if CovarianceMatrix with distances
  AddVariant(ShapeType, EARTH_ISOTROPIC);
  AddVariant(TrendType, UNREDUCED);
  setDI(NULL, allowedIfix, NULL);


 
  // epsC has been for internal reasons only ; essentially
  // the gencauchy model, except that 1 in the denominator 
  // is replaced by epsilon

  IncludeModel("declare", // never change name without checking all .cc, .R
	       // TrendType, // warum trend??
	       TcfType,
	       0, 0, 16, NULL, PREVMODEL_D, PREVMODEL_I,
	       checkdeclare, rangedeclare,  PREF_ALL, true,
	       PARAM_DEP, INFDIM, falsch, NOT_MONOTONE);
  //  kappanames("mean", REALSXP, "plane", REALSXP, "polydeg",		
  //             INTSXP, "polycoeff",					
  //REALSXP, "arbitrark2xyfct", CLOSXP, "fctcoeff", REALSXP);	
  kappanames(FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP, 
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP, 
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP);  
  addCov(0, declarefct, declarefct, declarefct, NULL, NULL);
  addCov(declarefctnonstat);
  AddVariant(TrendType, PREVMODEL_I);

 
  IncludePrim("EAxxA", ShapeType,  2, kappa_EAxxA, XONLY, CARTESIAN_COORD,
	      checkEAxxA, rangeEAxxA, PREF_NOTHING, 
	      PARAM_DEP, EaxxaMaxDim, (ext_bool) false, NOT_MONOTONE);
  nickname("eaxxa");
  addCov(EAxxA, NULL, NULL);
  kappanames("E", REALSXP, "A", REALSXP);
  addSpecial(minmaxEigenEAxxA);


  IncludePrim("EtAxxA",  ShapeType, 3, kappa_EtAxxA, XONLY, CARTESIAN_COORD,
	      checkEtAxxA, rangeEtAxxA, 3, EaxxaMaxDim, (ext_bool) false,
	      NOT_MONOTONE);
  nickname("etaxxa");
  addCov(EtAxxA, NULL, NULL);
  kappanames("E", REALSXP, "A", REALSXP, "alpha", REALSXP);
  addSpecial(minmaxEigenEtAxxA);


  IDCOORD = // ACHTUNG falls internal geaendert, auch in KeyInfo.cc aendern
    IncludeModel("id", ShapeType, 0, 0, 0, NULL, XONLY, PREVMODEL_I,
	       checkidcoord, NULL, PREF_NOTHING, 
	       false, PARAM_DEP, INFDIM, (ext_bool) false, NOT_MONOTONE);
  addCov(idcoord, NULL, NULL);


  TRAFO =
    //   IncludeModel("trafo", ManifoldType, 0, 1, 1, NULL, PREVMODEL_D, PREVMODEL_I,
   IncludeModel("trafo", ManifoldType, 0, 1, 1, kappatrafo,
		PARAMDEP_D, PARAMDEP_I,
		checktrafo, rangetrafo, PREF_ALL, 
		false, PARAM_DEP, INFDIM-1, (ext_bool) false, NOT_MONOTONE);
  kappanames("new", INTSXP);
  change_typeof(TRAFO_ISO, NN2); // ISO_NAMES
  addCov(trafo, NULL, NULL);
  addCov(nonstattrafo);// 
  addlogCov(logtrafo, lognonstattrafo, NULL);
  subnames("phi");
  setDI(allowedDtrafo, allowedItrafo, settrafo);
  addTypeFct(Typetrafo); 
 
  MULT_INVERSE =
    IncludeModel("mult_inverse", ShapeType, 1, 1, 0, NULL,
		 PREVMODEL_D, PREVMODEL_I,
		 checkmult_inverse, NULL, PREF_NOTHING,
		 true, SCALAR, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  addCov(mult_inverse, NULL, NULL); 
  addCov(mult_inverseNonstat);

 
  POLYGON =
    IncludeModel("polygon",  ShapeType, 0, 0, 1, NULL, XONLY, CARTESIAN_COORD, 
	       check_polygon, range_polygon, PREF_NOTHING,
	       false, SCALAR, 2, (ext_bool) true, MONOTONE);
  kappanames("lambda", REALSXP);
  addCov(Polygon, NULL, NULL, Inversepolygon, InversepolygonNonstat);
  RandomShape(INFTY, struct_polygon, init_polygon, doOK);   
 

  IncludePrim("rational", ShapeType, 2, kappa_rational,
	      XONLY, CARTESIAN_COORD,
	      checkrational, rangerational, INFDIM, (ext_bool) false, NOT_MONOTONE);
  addCov(rational, NULL, NULL);
  kappanames("A", REALSXP, "a", REALSXP);
  addSpecial(minmaxEigenrational);

  
  IncludePrim("rotat",  ShapeType, 2, kappa_rotat, XONLY, CARTESIAN_COORD,
	      checkrotat, rangerotat, PREF_NOTHING, SCALAR, 3, (ext_bool) false,
	      NOT_MONOTONE);
  addCov(rotat, NULL, NULL);
  kappanames("speed", REALSXP, "phi", REALSXP);
  addSpecial(minmaxEigenrotat);


  IncludePrim("Rotat",  ShapeType, 1, kappa_Rotat, XONLY, CARTESIAN_COORD,
	      checkRotat, rangeRotat, PARAM_DEP, 3, (ext_bool) false, NOT_MONOTONE);
  nickname("rotation");
  addCov(Rotat, NULL, NULL);
  kappanames("phi", REALSXP);

  
  SCATTER = IncludeModel("scatter",  PosDefType, 1, 1, 2, NULL, 
			 PREVMODEL_D, PREVMODEL_I, 
			 checkScatter, rangeScatter, PREF_ALL,
			 true, SUBMODEL_DEP, SUBMODEL_DEP,
			 (ext_bool) SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("size", REALSXP, "max", INTSXP);
  addCov(Scatter, NULL, NULL);
  RandomShape(1, struct_scatter, init_scatter, do_scatter,
	      false, false, false);
  //	      true, true, false); 

  RANDOMSIGN = 
    IncludeModel("sign",  ShapeType, 1, 1, 1, NULL, XONLY, PREVMODEL_I,
		 check_randomSign, range_randomSign, PREF_NOTHING,
		 false, SCALAR, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, NOT_MONOTONE);
  //nickname("");
  kappanames("p", REALSXP);
  addCov(randomSign, NULL, NULL, randomSignInverse, randomSignNonstatInverse);
  addlogCov(lograndomSign);
  RandomShape(1, struct_randomSign, init_randomSign, do_randomSign,
	      true, true, false); 
 

  SETPARAM = // nur Handling von Parametern von shape nach pts
    IncludeModel("setparam", ManifoldType, 1, 1, 1, NULL, 
		 PREVMODEL_D, PREVMODEL_I,
		 checksetparam,  range_setparam, PREF_ALL, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP,
		 (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  nickname("setparam");
  kappanames("performDo", INTSXP);
  addCov(setparamStat, Dsetparam, DDsetparam, D3setparam, D4setparam, 
	 Inverse_setparam, NonstatInverse_setparam);
  addCov(setparamNonStat);
  addTBM(NULL, spectralsetparam);
  RandomShape(INFTY, struct_failed, initsetparam, dosetparam,
	      false, false, true);
  addReturns(covmatrix_setparam, iscovmatrix_setparam);
  addTypeFct(Typesetparam);

  STROKORB_MONO =
    IncludeModel("m2r", ShapeType, 1, 1, 0, NULL, XONLY, ISOTROPIC,
		 checkstrokorb, NULL, PREF_NOTHING,
		 false, SCALAR, 3, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  addCov(strokorb, NULL, NULL); 
  RandomShape(1, structOK, init_strokorb, do_strokorb,
	      false, false, false);


#define balltest !true  
  IncludeModel("m3b", ShapeType, 1, 1, 
	       balltest ? 2 : 0, NULL, XONLY, ISOTROPIC,
	       checkstrokorbBall, 
	       balltest ? rangestrokorbball : NULL,/*for testing only*/
	       PREF_NOTHING,
	       false, SCALAR, 3, (ext_bool) true, MONOTONE);
  if (balltest) kappanames("min", REALSXP, "max", REALSXP); 
  addCov(strokorbBallInner, NULL, NULL); 
  RandomShape(1, struct_strokorbBall, init_failed, do_failed, do_random_failed,
	      false, false, false);
  
 
  STROKORB_BALL_INNER = // !! inverse scale gegenueber paper
    IncludeModel("r3binner", ShapeType, 1, 1, 1, NULL,
		 XONLY, CARTESIAN_COORD,
		 check_strokorbBallInner, range_strokorbBallInner, PREF_AUX,
		 true, 1, 1, (ext_bool) true, NOT_MONOTONE);
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
	       false, SCALAR, INFDIM, falsch, NOT_MONOTONE);
  addCov(strokorbPoly, NULL, NULL); 
  RandomShape(1, struct_strokorbPoly, init_failed, do_failed, do_random_failed,
	      false, false, false);
  

  
  TRUNCSUPPORT =
    IncludeModel("truncsupport", ShapeType, 
		 1, 1, 1, NULL, XONLY, PREVMODEL_I, checktruncsupport,
		 rangetruncsupport, PREF_NOTHING, false, SCALAR,
		 SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  kappanames("radius", REALSXP); // neg value == auto
  addCov(truncsupport, NULL, NULL, truncsupportInverse, StandardInverseNonstat);
  RandomShape(0, struct_truncsupport, init_truncsupport, do_truncsupport, false,
	      false, false);
  

  //////////////////////////////////////////////////
  // families of multivariate distribution; used in 

  // ACHTUNG!! addCov muss ganz zum Schluss !!

  ARCSQRT_DISTR = 
    IncludeModel("arcsqrt", RandomType, 0, 0, 1, NULL, 
		 DOMAIN_MISMATCH, ISO_MISMATCH, // set to "cart sys"
		 checkOK, range_arcsqrt, PREF_AUX, 
		 true, 1, 1, (ext_bool) false, MON_MISMATCH); // to do: nicht mismatch,
 //     sondern monotone im eindimensionalen
  kappanames("scale", REALSXP);
  RandomShape(0, structOK, init_arcsqrt, do_arcsqrt);
  addCov(arcsqrtD, arcsqrtDlog, arcsqrtDinverse, 
	 arcsqrtP, NULL, arcsqrtQ, arcsqrtR, NULL);
 
 
  DETERM_DISTR = 
    IncludeModel("determ", RandomType, 0, 0, 1, kappa_determ, 
		 DOMAIN_MISMATCH, ISO_MISMATCH, // set to "cart sys"
		 check_determ, range_determ, PREF_AUX,
		 false, SUBMODEL_DEP, INFDIM-1, (ext_bool) SUBMODEL_DEP, MON_MISMATCH);
  kappanames("mean", REALSXP);
  RandomShape(INFTY, structOK, init_determ, do_determ); 
  addCov(determD, determDlog, determDinverse, determP, determP2sided, determQ, 
	 determR, determR2sided);

 
  DISTRIBUTION = // FREEVARIABLE vorhanden. Muss extra in SpecialRMmodel.R
		 // definiert und nicht ueber generatemodels.R 
    IncludeModel("distr", RandomType, 0, 0, 16, kappa_distr, 
		 DOMAIN_MISMATCH, ISO_MISMATCH, // set to "cart sys"
		 check_distr, range_distr, PREF_AUX,
		 true, PARAM_DEP, INFDIM-1, (ext_bool) false, MON_MISMATCH);
  kappanames("name", STRSXP,
	     "nrow", INTSXP,
	     "ncol", INTSXP,
	     "ddistr", LANGSXP, 
	     "pdistr", LANGSXP,
	     "qdistr", LANGSXP, 
	     "rdistr", LANGSXP,
	     "envir", ENVSXP, 
	     FREEVARIABLE, REALSXP, // wird nie verwendet -- Puffer fuer 
	     // einfachen Algorithmus
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP, 
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP, FREEVARIABLE, REALSXP,
	     FREEVARIABLE, REALSXP//, FREEVARIABLE, REALSXP
	     ); // 7 free ones are remaining !	     
  
  change_sortof(DISTR_DX, IGNOREPARAM);
  change_sortof(DISTR_PX, IGNOREPARAM);
  change_sortof(DISTR_QX, IGNOREPARAM);
  change_sortof(DISTR_RX, IGNOREPARAM);
  
  RandomShape(0, structOK, init_distr, do_distr_do);
  addCov(distrD, distrDlog, distrDinverse, distrP, distrP2sided, distrQ,
	 distrR, distrR2sided);

 
  GAUSS_DISTR = // 139
    IncludeModel("normal", RandomType, 0, 0, 3, kappa_gauss_distr, 
		 DOMAIN_MISMATCH, ISO_MISMATCH, // set to "cart sys"
		 check_gauss_distr, range_gauss_distr, PREF_AUX, 
		 false, PARAM_DEP, INFDIM-1, (ext_bool) false, MON_MISMATCH);
  nickname("gauss");
  kappanames("mu", REALSXP, "sd", REALSXP, "log", INTSXP);
  RandomShape(INFTY, structOK, init_gauss_distr, do_gauss_distr);
  addCov(gaussD, gaussDlog, gaussDinverse, 
	 gaussP, gaussP2sided, gaussQ, gaussR, gaussR2sided);
 
  
  SET_DISTR = // nur Handling von Parametern von shape nach pts
    IncludeModel("setDistr", RandomType, 1, 1, 1, NULL, 
		 DOMAIN_MISMATCH, ISO_MISMATCH,
		 check_setParam, range_setParam, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_MISMATCH);
  subnames("to");
  kappanames("performDo", INTSXP);
  RandomShape(INFTY, structOK, init_setParam, do_setParam);
  addCov(setParamD, setParamDlog, setParamDinverse, 
	 setParamP, setParamP2sided, setParamQ,
	 setParamR, setParamR2sided);
  

  LOC =
    IncludeModel("loc", RandomType, 1, 1, 3, kappa_loc, 
		 DOMAIN_MISMATCH, ISO_MISMATCH,  // set to "cart sys"
		 check_loc, range_loc, PREF_AUX, 
		 false, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_MISMATCH);
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
		DOMAIN_MISMATCH, ISO_MISMATCH, // set to "cart sys"
		check_mcmc, range_mcmc, PREF_AUX, 
		false, PARAM_DEP,  INFDIM-1, (ext_bool) true, MON_MISMATCH);
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
		 DOMAIN_MISMATCH, ISO_MISMATCH, // set to "cart sys"
		 check_rectangular, range_rectangular, PREF_AUX, 
		 false, PARAM_DEP,  INFDIM-1, (ext_bool) true, MON_MISMATCH);
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
		 DOMAIN_MISMATCH, ISO_MISMATCH, // set to "cart sys"
		 check_RRspheric, range_RRspheric, PREF_AUX,
		 false, 1, 1, (ext_bool) true, MON_MISMATCH);
  kappanames("spacedim", INTSXP, "balldim", INTSXP, "R", REALSXP);
  RandomShape(INFTY, structOK, init_RRspheric, do_RRspheric);
  addCov(sphericD, sphericDlog, sphericDinverse, sphericP, NULL, sphericQ,
	 sphericR, NULL);
 
 
  UNIF = IncludeModel("unif", RandomType, 0, 0, 3, kappa_unif, 
		      DOMAIN_MISMATCH, ISO_MISMATCH, // set to "cart sys"
		      check_unif, range_unif, PREF_AUX, 
		      false, PARAM_DEP,  INFDIM-1, (ext_bool) true, MON_MISMATCH);
  kappanames("min", REALSXP, "max", REALSXP, "normed", INTSXP);
  RandomShape(INFTY, structOK, init_unif, do_unif); 
  addCov(unifD, unifDlog, unifDinverse, unifP, unifP2sided, unifQ, 
	 unifR, unifR2sided);

 
 
  // -----------------------------
  // shape + locations 
  // they *take* all very detailed frames like SmithType and pass
  // also their specific type to the submodel, in general
  // storage always pgs_storage !!
  MCMC_PGS =
    IncludeModel("MCMC_PGS", PointShapeType, 2, 2, 5, NULL, 
		 XONLY, CARTESIAN_COORD,
		 check_mcmc_pgs, range_mcmc_pgs, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_MISMATCH);
  // kappas must be identical to Zhou
  kappanames("density.ratio", REALSXP,  // stemming from gauss
	     "flat", INTSXP,
	     "infinitely_small", INTSXP,
	     "normed", INTSXP,
	     "mcmc_n", INTSXP);
  subnames("shape", "loc");
  addlogCov(logZhou);
  RandomShape(SUBMODEL_DEP, struct_mcmc_pgs, init_mcmc_pgs, 
	      do_mcmc_pgs, do_random_failed, true, true, false); 


  ZHOU = 
    IncludeModel("zhou", PointShapeType, 2, 2, 5, NULL, 
		 XONLY, CARTESIAN_COORD,
		 check_Zhou, range_Zhou, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_MISMATCH);
  kappanames("density.ratio", REALSXP,  // stemming from gauss
	     "flat", INTSXP,
	     "infinitely_small", INTSXP,
	     "normed", INTSXP,
	     "isotropic", INTSXP
	     );
  subnames("shape", "loc");
  addCov(Zhou, NULL, NULL);
  addlogCov(logZhou);
  RandomShape(SUBMODEL_DEP, struct_Zhou, init_Zhou, 
	      do_Zhou, do_random_failed, true, true, false);


  
  BALLANI = 
    IncludeModel("ballani", PointShapeType, 2, 2, 5, NULL, 
		 XONLY, CARTESIAN_COORD,
		 check_Ballani, range_Ballani, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_MISMATCH);
  kappanames("density.ratio", REALSXP,  // unklar welche parameter ueberhaupt gebraucht werden.
	     "flat", INTSXP,
	     "infinitely_small", INTSXP,
	     "normed", INTSXP,
	     "isotropic", INTSXP
	     );
  subnames("shape", "loc");
  addCov(Ballani, NULL, NULL);
  addlogCov(logBallani);
  RandomShape(SUBMODEL_DEP, struct_Ballani, init_Ballani, 
	      do_Ballani, do_random_failed, true, true, false);


  


  STANDARD_SHAPE = 
    IncludeModel("standardShape", PointShapeType, 1, 2, 0, NULL, 
		 XONLY, CARTESIAN_COORD,
		 check_standard_shape, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_MISMATCH);
  subnames("shape");
  addCov(standard_shape, NULL, NULL);
  addlogCov(logstandard_shape);
  RandomShape(SUBMODEL_DEP, struct_standard_shape, init_standard_shape, 
	      do_standard_shape, do_random_failed, true, true, false); 

  pref_type pmppp =  {0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 5, 0, 10, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
   MPPPLUS =
    IncludeModel("++",  PointShapeType, 1, MAXSUB, 1, kappamppplus, 
		 PREVMODEL_D, PREVMODEL_I, // CARTESIAN_COORD,
		 checkmppplus, rangempplus, pmppp, 
		 false, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, 
		 MON_MISMATCH);
  nickname("mppplus");
  kappanames("p", REALSXP);
  addCov(mppplus, NULL, NULL);

 
  STATIONARY_SHAPE = 
    IncludeModel("statShape", PointShapeType, 1, 1, 0, NULL, 
		 XONLY, CARTESIAN_COORD,
		 check_stationary_shape, NULL, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		 MON_MISMATCH);
  subnames("shape");
  addCov(stationary_shape, NULL, NULL);
  addlogCov(logstationary_shape);
  RandomShape(SUBMODEL_DEP, struct_stationary_shape, init_stationary_shape, 
	      do_stationary_shape, do_random_failed, true, true, false); 


  COVFCTN = // ALWAYS FIRST WITHIN FUNCTION WITHOUT GENUINE I NIT
    IncludeModel("Cov", InterfaceType, 1, 1, 0, NULL, XONLY, UNREDUCED, 
		 check_cov, NULL, PREF_AUX, 
		 INTERN_SHOW, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		 MON_MISMATCH);
  nickname("cov");
  addCov(Cov, NULL, NULL);
  RandomShape(2, struct_cov, init_cov, do_failed); 

  COVMATRIX = 
    IncludeModel("CovMatrix", InterfaceType, 1, 1, 0, NULL, XONLY,
		 UNREDUCED, //UNREDUCED,ISOTROPIC dependening on whether
		 // distances are givenXONLY, UNREDUCED,
		 check_covmatrix, NULL, PREF_AUX, 
		 INTERN_SHOW, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		 MON_MISMATCH);
  nickname("covmatrix");
  addCov(CovMatrix, NULL, NULL);
  RandomShape(struct_cov); 

  IncludeModel("Dummy", InterfaceType, 1, 1, 0, NULL, XONLY, UNREDUCED, 
	       check_dummy, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_MISMATCH);
  nickname("dummy");
  addCov(Dummy, NULL, NULL);
  RandomShape(struct_dummy); 

  RFGET = 
    IncludeModel("get", InterfaceType, 1, 1, 2, NULL, XONLY, UNREDUCED,
		 check_RFget, range_RFget, PREF_AUX, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		 MON_MISMATCH);  
  kappanames("up", INTSXP, "register", INTSXP);
  addCov(RFget, NULL, NULL);
  RandomShape(struct_RFget); 

 
  IncludeModel("Fctn", InterfaceType, 1, 1, 0, NULL, XONLY, UNREDUCED, 
		 check_fctn, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_MISMATCH);
  nickname("fctn");
  addCov(Fctn, NULL, NULL);
  RandomShape(structOK); 

  IncludeModel("Distr", InterfaceType, 1, 1, 5, kappa_EvalDistr,
	       XONLY, UNREDUCED,
	       check_EvalDistr, range_EvalDistr, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_MISMATCH);
  nickname("distr");
  kappanames("x", REALSXP, "q", REALSXP, "p", REALSXP, "n", REALSXP,
	     "dim", INTSXP);
  addCov(EvalDistr, NULL, NULL);
  RandomShape(struct_EvalDistr); 


  LIKELIHOOD_CALL =
  IncludeModel("loglikelihood", InterfaceType, 1, 1, 4, kappalikelihood, 
	       XONLY, UNREDUCED, check_likelihood, range_likelihood, PREF_AUX, 
	       INTERN_SHOW, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_MISMATCH);
  kappanames("data", LISTOF +  REALSXP, "estimate_variance", INTSXP,
	     "betas_separate", INTSXP, "ignore_trend", INTSXP);
  addCov(likelihood, NULL, NULL);
  RandomShape(struct_likelihood); 

  LINEARPART_CALL =
  IncludeModel("linearpart", InterfaceType, 1, 1, 0, NULL,
	       XONLY, UNREDUCED, check_linearpart, NULL, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_MISMATCH);
  addCov(linearpart, NULL, NULL);
  RandomShape(struct_linearpart); 
 
  PREDICT_CALL =
    IncludeModel("predict", InterfaceType, 0, 1, 1, NULL,
	       XONLY, UNREDUCED, check_predict, range_predict, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_MISMATCH);
  kappanames("register", INTSXP);
  addCov(predict, NULL, NULL);
  RandomShape(struct_predict); 
 

  
  IncludeModel("Pseudovariogram", InterfaceType, 1, 1, 0, NULL, 
	       XONLY, UNREDUCED,
	       check_cov, NULL, PREF_AUX, 
	       INTERN_SHOW, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_MISMATCH);
  nickname("pseudovariogram");
  addCov(Pseudovariogram, NULL, NULL);
  RandomShape(struct_variogram); 
 
  VARIOGRAM_CALL = // ALWAYS WITHIN FUNCTION WITHOUT GENUINE I NIT
    IncludeModel("Variogram", InterfaceType, 1, 1, 0, NULL, 
		 XONLY, UNREDUCED,
		 check_vario, NULL, PREF_AUX, 
		 INTERN_SHOW, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		 MON_MISMATCH);
  nickname("variogram");
  addCov(Variogram, NULL, NULL);
  RandomShape(struct_variogram);

  


  SIMULATE =
    IncludeModel("Simulate", InterfaceType, 1, 1, 3, NULL, 
		 XONLY, UNREDUCED, 
		 check_simulate, range_simulate, PREF_AUX, 
		 INTERN_SHOW, SUBMODEL_DEP, SUBMODEL_DEP,
		 (ext_bool) SUBMODEL_DEP,
		 MON_MISMATCH);
  nickname("simulate");
  kappanames("checkonly", INTSXP, "setseed", LANGSXP, "env", ENVSXP);
  addCov(simulate, NULL, NULL);
  RandomShape(struct_simulate); 

  /* 
  DENSITY =
    IncludeModel("Density", InterfaceType, 1, 1, 3, NULL, 
	       XONLY, UNREDUCED, 
	       check_density, range_density, PREF_AUX, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_MISMATCH);
  nickname("density");
  kappanames("log", INTSXP, "setseed", LANGSXP, "env", ENVSXP);
  addCov(density, NULL, NULL);
  RandomShape(struct_density); 
  */
  


  
  // ----------------------------
  // processes  //        CE CO CI TBM Sp di sq Ma av n mpp Hy spf any

  DOLLAR_PROC 
    = IncludeModel("$proc", ProcessType,		   
		   1, 1, 5, kappaS, // kappadollar,
		   XONLY, UNREDUCED, checkS, rangeS, PREF_ALL,
		   true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		   MON_SUB_DEP);
  // do not change Order!!
  addSpecific(DOLLAR);
  RandomShape(2, structSproc, initSproc, doSproc, true, true, true);
  AddVariant(GaussMethodType, UNREDUCED);

  
  PLUS_PROC = 
    IncludeModel("plusproc", ProcessType, 1, MAXSUB, 0, NULL, 
		  XONLY, UNREDUCED,
		 checkplusproc, NULL, PREF_ALL, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		 MON_SUB_DEP);
  addSpecific(PLUS);
  RandomShape(2, structplusproc, initplusproc, doplusproc, false, false, true);
  AddVariant(GaussMethodType, UNREDUCED);

  // 
  PROD_PROC =
  IncludeModel("prodproc", ProcessType, 1, 1, 0, NULL, 
	       XONLY, UNREDUCED,
	       checkprodproc, NULL, PREF_ALL, 
	       true, PARAM_DEP, INFDIM-1, (ext_bool) false, NOT_MONOTONE);
  addSpecific(PROD);
  RandomShape(2, structprodproc, initprodproc, doprodproc, false, false, true);
  AddVariant(GaussMethodType, UNREDUCED);

  IncludeModel("trafoproc", ProcessType, 1, 1, 1, kappatrafo, 
	       XONLY, UNREDUCED,
	       checktrafoproc, rangetrafo, PREF_ALL, 
	       true, PARAM_DEP, INFDIM-1, (ext_bool) false, NOT_MONOTONE);
  addSpecific(TRAFO);
  RandomShape(2, structtrafoproc, inittrafoproc, dotrafoproc,
  	      false, false, true);
  AddVariant(GaussMethodType, UNREDUCED);
   
  
  MPPPLUS_PROC =
  IncludeModel("mppplusproc", ProcessType, 1, MAXSUB, 1, kappamppplus, 
		  XONLY, UNREDUCED,
		 checkmppplus, rangempplus, PREF_ALL, 
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_SUB_DEP);
  nickname("mppplus");
  // addSpecific(MPPPLUS);
  kappanames("p", REALSXP);
  RandomShape(2, struct_mppplus, init_mppplus, do_mppplus, true, true, true);
   

  MULT_PROC = 
  IncludeModel("multproc", ProcessType, 1, MAXSUB, 1, NULL, 
	       XONLY, UNREDUCED,
	       checkmultproc, rangemultproc, PREF_ALL, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_SUB_DEP);
  addSpecific(MULT);
  kappanames("multicopies", INTSXP);
  RandomShape(2, structmultproc, initmultproc, domultproc, false, false, true);
  AddVariant(GaussMethodType, UNREDUCED);

  
  M_PROC = 
    IncludeModel("matrixproc", ProcessType, 1, MAXSUB, 2, kappaM,
	       XONLY, UNREDUCED,
	       checkM, rangeM, PREF_ALL, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_SUB_DEP);
  addSpecific(MATRIX);
  RandomShape(2, structMproc, initMproc, doMproc, false, false, true);
  AddVariant(GaussMethodType, UNREDUCED);
 
  
  VAR2COV_PROC = 
  IncludeModel("covproc", ProcessType, 1, MAXSUB, 2, NULL, 
	       XONLY, UNREDUCED,
	       checkvar2covproc, rangevariogram2cov, PREF_ALL, 
	       true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_SUB_DEP);
  addSpecific(VARIOGRAM2COV);
  RandomShape(2, structvar2covproc, initvar2covproc, dovar2covproc,
	      false, false, true);
  AddVariant(GaussMethodType, UNREDUCED);
 
 
  TREND_PROC = 
    IncludeModel(METHOD_NAMES[Trendproc], ProcessType, 
		 0, 1, 1, kappatrend, XONLY, UNREDUCED,
		 checkTrendproc, rangetrend, PREF_TREND,
		 false, SUBMODEL_DEP, INFDIM-1, (ext_bool) false,
		 MON_MISMATCH);
  addSpecific(TREND);
  RandomShape(2, init_Trendproc, do_Trendproc);
  AddVariant(GaussMethodType, UNREDUCED);


  //////////// STANDARD PROCESSES ////////////////

  
  AVERAGE_USER = 
    IncludeModel(METHOD_NAMES[Average], GaussMethodType, 1, 2, 3, kappaGProc, 
		 XONLY, UNREDUCED,
		 check_randomcoin, range_randomcoin, PREF_NOTHING,  
		 false, SCALAR, MAXMPPDIM, (ext_bool) false, MON_MISMATCH);
  subnames("phi", "shape");
  kappanames("boxcox", REALSXP, "intensity", REALSXP, "method", INTSXP);  
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  //  change_typeof(RANDOMCOIN_INTENSITY, RandomType);
 // addCov(coin, NULL, NULL, coinInverse);
  RandomShape(2, struct_gaussmethod, init_gaussprocess, do_gaussprocess); 

  
  
  RANDOMCOIN_USER = CopyModel(METHOD_NAMES[RandomCoin], AVERAGE_USER);
  
  AVERAGE_INTERN =  CopyModel("averageIntern", AVERAGE_USER);
  addSpecific(RANDOMCOIN_USER, false);
  addSpecific(AVERAGE_USER, false);  // see gauss.cc for the use of Specific!
  RandomShape(2, struct_randomcoin, init_randomcoin, dompp, true, true, false);


  CIRCEMBED = // und die anderen fehlen auch noch !!
    IncludeModel(METHOD_NAMES[CircEmbed], GaussMethodType, 1, 1, 13, kappa_ce,
		 XONLY, UNREDUCED,
		 check_ce, range_ce, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXCEDIM, (ext_bool) false, MON_MISMATCH);
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
    IncludeModel(METHOD_NAMES[CircEmbedCutoff], GaussMethodType, 1, 1, 15, 
		 kappa_localproc, XONLY, UNREDUCED,
		 check_local_proc, range_co_proc, PREF_NOTHING,
		 false, SCALAR, MAXCEDIM,  (ext_bool) false, MON_MISMATCH);
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
  RandomShape(2, struct_gaussmethod, init_gaussprocess, do_gaussprocess); 
 
  CE_CUTOFFPROC_INTERN = CopyModel("cutoffIntern", CE_CUTOFFPROC_USER);
  addSpecific(CE_CUTOFFPROC_USER, false);
  RandomShape(2, struct_ce_approx, init_ce_approx, do_ce_approx);


  CE_INTRINPROC_USER =  
    IncludeModel(METHOD_NAMES[CircEmbedIntrinsic], GaussMethodType, 
		 1, 1, 15, kappa_localproc, 
		 XONLY, UNREDUCED,
		 check_local_proc, range_intrinCE, PREF_NOTHING,
		 false, SCALAR, MAXCEDIM, (ext_bool) false, MON_MISMATCH);
  nickname(METHOD_NAMES[CircEmbedIntrinsic]);
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
  RandomShape(2, struct_gaussmethod, init_gaussprocess, do_gaussprocess); 
 
  CE_INTRINPROC_INTERN = CopyModel("intrinsIntern", CE_INTRINPROC_USER);
  addSpecific(CE_INTRINPROC_USER, false);
  RandomShape(2, struct_ce_approx, init_ce_approx, do_ce_approx);
 

  DIRECT = 
    IncludeModel(METHOD_NAMES[Direct], GaussMethodType, 
		 1, 1, 1, kappaGProc, XONLY, UNREDUCED,
		 check_directGauss, range_direct, PREF_NOTHING,
		 false,  SUBMODEL_DEP, INFDIM-1, (ext_bool) false,
		 MON_MISMATCH);
  kappanames("boxcox", REALSXP);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, init_directGauss, do_directGauss);


  HYPERPLANE_USER  = 
    IncludeModel(METHOD_NAMES[Hyperplane], GaussMethodType, 1, 1, 6, kappaGProc,
		 XONLY, UNREDUCED,
		 check_hyperplane, range_hyperplane, PREF_NOTHING,  
		 false, SCALAR, 2, (ext_bool) false, MON_MISMATCH);
  kappanames("boxcox", REALSXP,
	     "superpos", INTSXP, "maxlines", INTSXP, "mar_distr", INTSXP, 
	     "mar_param", REALSXP, "additive", INTSXP);
  //  addCov(IdStat, NULL, NULL, IdInverse);
  //  addCov(IdNonStat);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_gaussmethod, init_gaussprocess, do_gaussprocess); 


  HYPERPLANE_INTERN =   
    CopyModel("hyperIntern", HYPERPLANE_USER, check_hyperplane_intern);
  addSpecific(HYPERPLANE_USER, false);
  RandomShape(2, struct_hyperplane, init_hyperplane, do_hyperplane);
  
  
  NUGGET_USER  = 
    IncludeModel(METHOD_NAMES[Nugget], GaussMethodType, 
		 1, 1, 3, kappaGProc, XONLY, UNREDUCED,
		 check_nugget_proc, range_nugget_proc, PREF_NOTHING, false, 
		 PREVMODEL_DEP, INFDIM, (ext_bool) true, MON_MISMATCH);
  kappanames("boxcox", REALSXP,"tol", REALSXP, "vdim", INTSXP);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_gaussmethod, init_gaussprocess, do_gaussprocess); 


  NUGGET_INTERN = CopyModel("nuggetIntern", NUGGET_USER);
  addSpecific(NUGGET_USER, false);
  RandomShape(2, struct_nugget, init_nugget, do_nugget);
  /* see simu.cc, CMbuild for special treatment of nugget when
     users choice is given */  
  /* cf. convert.R, PrepareModel, near end of function */
  

 
  SEQUENTIAL = 
    IncludeModel(METHOD_NAMES[Sequential], GaussMethodType, 1, 1, 3,
		 kappaGProc, 
		 XONLY, UNREDUCED,
		 check_sequential, range_sequential, PREF_NOTHING,  
		 false, SCALAR, INFDIM-1, (ext_bool) false, MON_MISMATCH);
  kappanames("boxcox", REALSXP, "back_steps", INTSXP, "initial", INTSXP);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, init_sequential, do_sequential);


  SPECTRAL_PROC_USER = 
    IncludeModel(METHOD_NAMES[SpectralTBM], GaussMethodType,  1, 1, 5,
		 kappaGProc,
		 XONLY, UNREDUCED,
		 check_spectral, range_spectral, PREF_NOTHING,  
		 false, SCALAR, MAXTBMSPDIM, (ext_bool) false, MON_MISMATCH);
  kappanames("boxcox", REALSXP,"sp_lines", INTSXP, "sp_grid", INTSXP,
	     "prop_factor", REALSXP, "sigma", REALSXP );
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(2, struct_gaussmethod, init_gaussprocess, do_gaussprocess); 
 

  SPECTRAL_PROC_INTERN = CopyModel("spectralIntern", SPECTRAL_PROC_USER);
  addSpecific(SPECTRAL_PROC_USER, false);
  RandomShape(2, struct_spectral, init_spectral, do_spectral);


  SPECIFIC = 
    IncludeModel(METHOD_NAMES[Specific], GaussMethodType, 1, 1, 1, kappaGProc,
		 XONLY, UNREDUCED,
		 check_specificGauss, range_specificGauss, PREF_NOTHING,  
		 false, SUBMODEL_DEP, MAXTBMSPDIM, (ext_bool) false,
		 MON_MISMATCH);
  RandomShape(2, struct_specificGauss, init_specificGauss, do_specificGauss);
  kappanames("boxcox", REALSXP);
  change_sortof(GAUSS_BOXCOX, ANYPARAM);

  
  TBM_PROC_USER = 
    IncludeModel(METHOD_NAMES[TBM], GaussMethodType,  1, 1, 9, tbm_kappasproc, 
		 XONLY, UNREDUCED, 
		 checktbmproc, rangetbmproc, PREF_NOTHING,
		 false, PARAM_DEP, SUBMODEL_DEP, (ext_bool) false,
		 MON_MISMATCH);
  kappanames("boxcox", REALSXP, 
	     "fulldim", INTSXP, "reduceddim", INTSXP, "layers", INTSXP,
	     "lines", INTSXP, "linessimufactor", REALSXP,
	     "linesimustep",  REALSXP,  // "grid", INTSXP,
	     "center", REALSXP, "points", INTSXP); 
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  change_sortof(TBM_LAYERS, ONLYRETURN); // NA will not be estimated
  RandomShape(2, struct_gaussmethod, init_gaussprocess, do_gaussprocess); 
 
  TBM_PROC_INTERN = CopyModel("tbmIntern", TBM_PROC_USER);
  addSpecific(TBM_PROC_USER);
  RandomShape(2, struct_tbmproc, init_tbmproc, do_tbmproc);
  

  gaussmethod[CircEmbed] = CIRCEMBED;
  gaussmethod[CircEmbedCutoff]= CE_CUTOFFPROC_INTERN;
  gaussmethod[CircEmbedIntrinsic] = CE_INTRINPROC_INTERN;
  gaussmethod[TBM] = TBM_PROC_INTERN;
  gaussmethod[SpectralTBM] = SPECTRAL_PROC_INTERN;
  gaussmethod[Direct] = DIRECT;
  gaussmethod[Sequential] = SEQUENTIAL;
  gaussmethod[Trendproc] = TREND_PROC;
  gaussmethod[Average] = AVERAGE_INTERN;
  gaussmethod[Nugget] = NUGGET_INTERN;
  gaussmethod[RandomCoin] = AVERAGE_INTERN;
  gaussmethod[Hyperplane] = HYPERPLANE_INTERN;
  gaussmethod[Specific] = SPECIFIC;
  gaussmethod[Nothing] =  gaussmethod[Forbidden] = MISMATCH;
 

  BRNORMED =
    IncludeModel("loggaussnormed", NormedProcessType, 1, 1, 5, kappabrnormed, 
		 XONLY, UNREDUCED,
		 check_brnormed, range_brnormed, PREF_NOTHING, 
		 false, SUBMODEL_DEP, MAXMPPDIM, (ext_bool) false,
		 MON_MISMATCH);
  kappanames("prob", REALSXP,
	    "optimize_p", INTSXP,
	    "nth", INTSXP,
	    "burn.in",  INTSXP,
	    "rejection", INTSXP
	    );
  subnames("variogram");
  addlogCov(logZhou);
  //RandomShape(SUBMODEL_DEP, struct_brnormed, init_brnormed, // 2.2.19
  //	      do_brnormed, do_random_failed, true, true, false);
  RandomShape(0, struct_brnormed, init_brnormed, do_brnormed, finalmaxstable);
  
  // non sub-gaussian processe
   

  BRORIGINAL_USER =
    IncludeModel("brorig", BrMethodType, 1, 2, 3, NULL, 
		 XONLY, UNREDUCED,
		 checkBrownResnickProc, range_mpp, PREF_NOTHING,
		 false, SCALAR, MAXMPPDIM, (ext_bool) false, MON_MISMATCH);
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP);
  subnames("phi", "tcf");
  RandomShape(0, structBRuser, initBRuser, dompp, finalmaxstable); 
  addlogD(loglikelihoodBR);
 
  BRORIGINAL_INTERN =
    CopyModel("brorigIntern", BRORIGINAL_USER, PointShapeType);
  make_internal();
  nickname("brorig");
  RandomShape(SUBMODEL_DEP, structBRintern, init_BRorig, do_BRorig);

  
  BRMIXED_USER =
    IncludeModel("brmixed", BrMethodType, 1, 2, 10, kappaBRmixed, 
		 XONLY, UNREDUCED, 
		 check_BRmixed, range_BRmixed, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, (ext_bool) false,
		 MON_MISMATCH);
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP, 
	     "meshsize", REALSXP, "vertnumber", INTSXP,
             "optim_mixed", INTSXP, "optim_mixed_tol", REALSXP, 
             "lambda", REALSXP, "areamat", REALSXP, "variobound", REALSXP);
  subnames("phi", "tcf");
  RandomShape(0, structBRuser, initBRuser, dompp, finalmaxstable);
  addlogD(loglikelihoodBR);
    
  BRMIXED_INTERN =
    CopyModel("brmixedIntern", BRMIXED_USER, PointShapeType); 
  make_internal();
  nickname("brmixed");
  RandomShape(SUBMODEL_DEP, structBRintern, init_BRmixed, do_BRmixed);

  
  BRSHIFTED_USER = 
    IncludeModel("brshifted", BrMethodType, 
		 1, 2, 3, NULL, 
		 XONLY, UNREDUCED, 
		 checkBrownResnickProc, range_mpp, PREF_NOTHING,
		 false, SCALAR, MAXMPPDIM, (ext_bool) false, MON_MISMATCH);
  subnames("phi", "tcf");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP);
  RandomShape(0, structBRuser, initBRuser, dompp, finalmaxstable);
  addlogD(loglikelihoodBR);
 
  BRSHIFTED_INTERN =
    CopyModel("brshiftIntern", BRSHIFTED_USER, PointShapeType); 
  make_internal();
  nickname("brshif");
  RandomShape(SUBMODEL_DEP, structBRintern, init_BRshifted, do_BRshifted);
  
   
  BROWNRESNICKPROC =
    IncludeModel("brownresnick", BrMethodType, // ProcessType
		 1, 2, 3, NULL, 
		 XONLY, UNREDUCED, 
		 checkBrownResnickProc, range_mpp, PREF_NOTHING,  
		 false, SCALAR, MAXMPPDIM, (ext_bool) false, MON_MISMATCH);
  subnames("phi", "tcf");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP);
  //  addCov(BrownResnick, NULL, NULL);
  RandomShape(0, structBrownResnick, initBrownResnick, doBrownResnick,
	      finaldoBrownResnick); 
  addlogD(loglikelihoodBR);


  // distributions

  
  BINARYPROC = // direct an Gauss gekoppelt!!
    IncludeModel("binaryprocess", NormedProcessType, 1, 1, 3,
		 kappa_binaryprocess, 
		 XONLY, UNREDUCED, 
		 checkbinaryprocess, rangebinaryprocess, PREF_NOTHING,
		 false, SUBMODEL_DEP, INFDIM, (ext_bool) false, MON_MISMATCH);
  nickname( "bernoulli");
  kappanames(INTERNAL_PARAM, REALSXP, "stationary_only", 
	     INTSXP, "threshold", REALSXP);
  //"p", REALSXP);
  RandomShape(INFTY, struct_binaryprocess, init_binaryprocess,
	      do_binaryprocess);

  GAUSSPROC = // never change names. See fitgauss.R, for instance
    IncludeModel("gauss.process", GaussMethodType, // formerly Processtype
		 1, 1, 2, kappaGProc,
		 XONLY, UNREDUCED,
		 checkgaussprocess, rangegaussprocess, PREF_NOTHING,
		 false, SUBMODEL_DEP, INFDIM, (ext_bool) false, MON_MISMATCH);
  nickname("gauss");
  kappanames("boxcox", REALSXP, "stationary_only", INTSXP);  
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  change_sortof(GAUSSPROC_STATONLY, ONLYRETURN); // NA will not be estimated
  RandomShape(2, struct_gaussprocess, init_gaussprocess, do_gaussprocess);
  addlogD(gaussprocessDlog);
  
  
  POISSONPROC =
    IncludeModel("poisson", PoissonType, 1, 1, 1, NULL,
		 XONLY, UNREDUCED,  
		 check_poisson, range_poisson, PREF_NOTHING,
		 false, SUBMODEL_DEP, MAXMPPDIM, (ext_bool) false,
		 MON_MISMATCH);
  kappanames("intensity", REALSXP);
  RandomShape(0, struct_poisson, init_poisson, dompp);
 

  SCHLATHERPROC =
    IncludeModel("extremalgauss", SchlatherType, 1, 2, 3, NULL, 
		 XONLY, UNREDUCED,  
		 check_schlather, range_mpp, PREF_NOTHING, 
		 false, SCALAR, INFDIM, (ext_bool) false, MON_MISMATCH);
  nickname("schlather");
  subnames("phi", "tcf");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP);
  addCov(extremalgaussian, NULL, NULL);
  RandomShape(0, struct_schlather, init_mpp, dompp, finalmaxstable);
  addlogD(loglikelihoodSchlather);

  

  EXTREMALTPROC =
    IncludeModel("extremalt", SchlatherType, 1, 1, 4, NULL, 
		 XONLY, UNREDUCED,  
		 check_schlather, range_opitz, PREF_NOTHING, 
		 false, SCALAR, INFDIM, (ext_bool) false, MON_MISMATCH);
  nickname("opitz");
  subnames("phi");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP, "alpha", REALSXP);
  addCov(extremalgaussian, NULL, NULL);
  RandomShape(0, struct_schlather, init_opitzprocess, dompp, finalmaxstable);
  addlogD(loglikelihoodSchlather);

 
  SMITHPROC =
    IncludeModel("smith", SmithType, 1, 2, 3, NULL, XONLY, UNREDUCED,
		 check_smith, range_mpp, PREF_NOTHING,  
		 false, SCALAR, MAXMPPDIM, (ext_bool) false, MON_MISMATCH);
  subnames("shape", "tcf");
  kappanames("xi", REALSXP, "mu", REALSXP, "s",  REALSXP
	     //, "intensity", REALSXP
	     );
  //  change_typeof(SMITH_INTENSITY, RandomType);
  RandomShape(0, struct_smith, init_mpp, dompp, finalmaxstable);

  
  CHI2PROC =
    IncludeModel("chi2", ProcessType, 1, 1, 2, kappaGProc, XONLY, UNREDUCED,
		 checkchisqprocess, rangechisqprocess, PREF_NOTHING,
		 false, SUBMODEL_DEP, INFDIM, (ext_bool) false, MON_MISMATCH);
  kappanames("boxcox", REALSXP, "f", INTSXP);  
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(0, struct_chisqprocess, init_chisqprocess, do_chisqprocess);

  
  TPROC =
    IncludeModel("t", ProcessType, 1, 1, 2, kappaGProc, XONLY, UNREDUCED,
		 checkchisqprocess, rangetprocess, PREF_NOTHING,
		 false, SUBMODEL_DEP, INFDIM, (ext_bool) false, MON_MISMATCH);
  kappanames("boxcox", REALSXP, "nu", REALSXP); 
  change_sortof(GAUSS_BOXCOX, ANYPARAM);
  RandomShape(0, struct_chisqprocess, init_chisqprocess, do_tprocess);
  
}


