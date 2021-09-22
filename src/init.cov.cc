
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part

 Copyright (C) 2001 -- 2017 Martin Schlather

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


#include <stdio.h>  
#include <string.h>
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "RF.h"
#include "operator.h"
#include "startGetNset.h"
#include "primitive.h"
#include "shape.h" // wegen fixcov

int
BROWNIAN,
  CAUCHY,
  CUTOFF,
  GAUSS,
  GENERALISEDCAUCHY,
  GENNSST_INTERN, 
  GNEITING_INTERN,
  MATRIX,
  NATSC_INTERN, NATSC_USER,
  NUGGET,
  MULT,
  POW,
  POWER_DOLLAR,
  PROD,
  STABLE, 
  STEIN,
  TBM_OP,
  TBM2NR,
  VARIOGRAM2COV,
  VECTOR;


void includeCovModels() {
  // at the origin: prespecified distribution
  // old RandomFields ave1,ave2
  IncludeModel("ave",  PosDefType, 1, 1, 3, kappa_ave, XONLY, SYMMETRIC,
	       checkave, rangeave, PREF_ALL, 
	       false, SCALAR, AveMaxDim, (ext_bool) false, NOT_MONOTONE);
  kappanames("A", REALSXP, "z", REALSXP, "spacetime", INTSXP);
  addCov(ave, NULL, NULL);
  RandomShape(structAve, true);


  SHAPEAVE =
    IncludeModel("shape.ave", ShapeType, 1, 2, 3, kappa_ave,
		 XONLY, CARTESIAN_COORD,  
		 check_shapeave, rangeave, PREF_NOTHING,
		 true, SCALAR, INFDIM-1, (ext_bool) true, NOT_MONOTONE);
  kappanames("A", REALSXP, "z", REALSXP, "spacetime", INTSXP);
  subnames("phi", "gauss");
  addlogCov(logshapeave);
  RandomShape(0, init_shapeave, do_shapeave, true);


  pref_type pbcw = {2, 5, 5, 5, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeScalar("bcw", VariogramType, 0, 0, 3, XONLY, ISOTROPIC,
	       checkbcw, rangebcw, pbcw,
	       INFDIM, (ext_bool) false, NORMAL_MIXTURE); // todo part is even
  // LAPLACE
  kappanames("alpha", REALSXP, "beta", REALSXP, "c", REALSXP);
  addCov(bcw, Dbcw, DDbcw, Inversebcw);
  addLocal(coinitbcw, ieinitbcw);
  AddVariant(PosDefType, ISOTROPIC);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);

  IncludeScalar("locstatfbm", PosDefType, 0, 0, 2, XONLY, ISOTROPIC, 
	       checklsfbm, rangelsfbm, PREF_ALL, INFDIM, (ext_bool) false,
	       MONOTONE); 
  nickname("lsfbm");
  kappanames("alpha", REALSXP, "const", REALSXP);
  addCov(lsfbm, Dlsfbm, DDlsfbm, D3lsfbm, D4lsfbm, Inverselsfbm);
  Taylor(-1, RF_NA, 0, 0);
  TailTaylor(-1, RF_NA, 0, 0);
   RandomShape(0, struct_failed, initlsfbm, do_failed, false, true, false);



  pref_type
    pbessel = {2, 0, 0,  0, 5, 3, 3, 0, 5, 0, 5, 0, 0, 5};
  //            CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("bessel",  PosDefType, 1, XONLY, ISOTROPIC, 
	      checkBessel, rangeBessel,
	      pbessel, SCALAR, INFDIM, (ext_bool) false, NOT_MONOTONE);
  kappanames("nu", REALSXP);
  addCov(Bessel, NULL, NULL);
  addTBM(initBessel, spectralBessel);	       
  RandomShape(0, struct_failed, initBessel, do_failed, false, true, false);
   
  IncludeModel("bigneiting", PosDefType, 0, 0, 8, kappa_biGneiting, XONLY,
	       ISOTROPIC, checkbiGneiting, rangebiGneiting, PREF_ALL, 
	       false, 2, PARAM_DEP, (ext_bool) true, NOT_MONOTONE);
  addCov(biGneiting, DbiGneiting, DDbiGneiting, NULL, NULL);
  kappanames("kappa", INTSXP,
	     "mu", REALSXP,
	     "s", REALSXP, "sred12", REALSXP,
	     "gamma", REALSXP,
	     "cdiag", REALSXP, "rhored", REALSXP, "c", REALSXP);
  add_sortof(sortof_biGneiting);
  RandomShape(0, struct_failed, initbiGneiting, do_failed, false, true, false);
   

  pref_type
    pbernoulli = {5, 0, 0,  0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //             CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("bernoulli", TcfType, 1, 1, 3, NULL, SUBMODEL_D, SUBMODEL_I,
	       checkbinary, rangebinary, pbernoulli,
	       false, SCALAR, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  kappanames("threshold", REALSXP, "correlation", INTSXP, "centred", INTSXP);
  addCov(binary, NULL, NULL);

  IncludeModel("biWM",  PosDefType, 0, 0, 8, kappa_biWM, XONLY, ISOTROPIC,
	       checkbiWM2, rangebiWM2, PREF_ALL,
	       false, 2, INFDIM, (ext_bool) false, NOT_MONOTONE);
  nickname("biwm");
  addCov(biWM2, biWM2D, biWM2DD, biWM2D3, biWM2D4,  ErrInverse);
  addLocal(coinitbiWM2, NULL);
  kappanames("nudiag", REALSXP, "nured12", REALSXP, 
	     "nu", REALSXP, // or lower triangle
	     "s", REALSXP,  // lower triangle definition
	     "cdiag", REALSXP, "rhored", REALSXP,
	     "c", REALSXP,  // or lower triangle
	     "notinvnu", INTSXP);
  RandomShape(0, struct_failed, initbiWM2, do_failed, false, true, false);
 
  // bivariate stable or bivariate exponetial model
  //IncludeModel("bistable", PosDefType, 0,0, 4, kappa_biStable, XONLY, ISOTROPIC,
  IncludeModel("bistable", PosDefType, 0,0, 7, kappa_biStable, XONLY, ISOTROPIC,
	       checkbiStable, rangebiStable, PREF_ALL,
	       false, 2, 3, (ext_bool) false, MONOTONE);
  kappanames("alpha", REALSXP,
             "s", REALSXP,  // lower triangle definition
	     "cdiag", REALSXP,
             "rho", REALSXP,
	     "rhored", REALSXP,
             "betared", REALSXP,
             "alphadiag", REALSXP
             );
  addCov(biStable, DbiStable, DDbiStable, D3biStable, D4biStable, NULL);
  addLocal(coinitbiStable, NULL);
  RandomShape(0, struct_failed, initbiStable, do_failed, false, true, false); //what is this???

  pref_type pblend = {  0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                   CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("blend", PosDefType, 2, 2, 1, kappablend, KERNEL, SYMMETRIC,
		checkblend, rangeblend, pblend,
	       false, SCALAR, SUBMODEL_DEP, falsch, NOT_MONOTONE);
  subnames("multi", "blend");
  kappanames("thresholds", REALSXP); // neg value == auto
  addCov(nonstatblend); 
  change_sortof(BLEND_THRES, ONLYRETURN);
 

  pref_type
    pbrownresnick = {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  BROWNRESNICK =
    IncludeModel("brownresnick", TcfType, 1, 1, 0, NULL, XONLY, SUBMODEL_I,
		 checkbrownresnick, NULL , pbrownresnick, false,
		 SCALAR, SUBMODEL_DEP, (ext_bool) false, MON_SUB_DEP);
  addCov(brownresnick, Dbrownresnick, DDbrownresnick, D3brownresnick, 
	 NULL, NULL);
  RandomShape(0, struct_brownresnick, init_brownresnick, do_brownresnick);
  //  Taylor(0, 0, 0, 0, 0, 0, 0, 0);

  
  IncludeScalar("br2bg",  PosDefType, 1, 1, 0, XONLY, SUBMODEL_I, 
		 check_BR2BG, NULL, PREF_ALL, SUBMODEL_DEP, (ext_bool) false, MON_SUB_DEP);
  addCov(BR2BG, NULL, NULL);


  IncludeScalar("br2eg", PosDefType, 1, 1, 0,  XONLY, SUBMODEL_I, 
	       check_BR2EG, NULL, PREF_ALL, SUBMODEL_DEP, (ext_bool) false, MON_SUB_DEP);
  addCov(BR2EG, NULL, NULL);
 
  
  pref_type pbubble = {  0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                   CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  BUBBLE =
    IncludeModel("bubble", PosDefType, 1, 2, 4, kappabubble, KERNEL, UNREDUCED,
		 checkbubble, rangebubble, pbubble,
		 false, SCALAR, INFDIM - 1, falsch, NOT_MONOTONE);
  subnames("phi", "scaling");
  kappanames("z", REALSXP, "weight", REALSXP, "minscale", REALSXP, "barycentre",
	     INTSXP);
  addCov(nonstatbubble);
  change_sortof(BUBBLE_MINSCALE, ONLYRETURN);
  change_sortof(BUBBLE_WEIGHT, ONLYRETURN);
  change_sortof(BUBBLE_Z, ONLYRETURN);
  // RandomShape(0, struct_failed, initBubble, do_failed, false, true, false);


  
  pref_type pcauchy=  {2, 5, 0,  2, 0, 4, 0, 0, 0, 0, 0, 0, 0, 5};
  //                   CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("cauchy", PosDefType, 1, XONLY, ISOTROPIC, 
	      checkCauchy, rangeCauchy, pcauchy, 
	      SCALAR, INFDIM, (ext_bool) false, NORMAL_MIXTURE);
  kappanames("gamma", REALSXP);
  addCov(Cauchy, DCauchy, DDCauchy, InverseCauchy);
  //  addlogCov(logCauchy);
  addTBM(TBM2Cauchy);
  addLocal(coinitCauchy, NULL);
  addGaussMixture(DrawMixCauchy, LogMixDensCauchy);
	       
 // todo : WARUM HABE ICH DIESES MODELL CODIERT, ABER NICHT ONLINE??     
  //  pref_type pctbm={2, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //     //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  // IncludePrim("cauchytbm", PosDefType,  3, XONLY, ISOTROPIC, checkOK,
  //	      rangeCauchytbm, pctbm, SCALAR, INFDIM, (ext_bool) false);
  // kappanames("alpha", REALSXP, "beta", REALSXP, "gamma", REALSXP);
  // addCov(Cauchytbm, DCauchytbm, InverseCauchy); // scale not correct, but
  // should be an approximation that is good enough
	      
  IncludePrim("circular",  TcfType, 0, XONLY, ISOTROPIC,
	      checkOK, NULL, 2, (ext_bool) false, GNEITING_MON);
  addCov(circular, Dcircular, ScaleOne);
  RandomShape(structcircular);


  // IncludePrim("cone", PosDefType,  3, XONLY, ISOTROPIC, checkcone, rangecone);
  //  kappanames("r", REALSXP, "socle", REALSXP, "height", REALSXP);
  // RandomShape(init_cone, mppget_cone, sd_standard, MPP_POISS);

  IncludePrim("CDeWijsian",  VariogramType, 2, NULL, XONLY, ISOTROPIC, 
	      checkdewijsian,  rangeDeWijsian, PREF_NOTHING, 
	      SCALAR, INFDIM, (ext_bool) false, MONOTONE); 
  nickname("cdewijs");
  make_internal();
  kappanames("alpha", REALSXP, "range", REALSXP);
  addCov(DeWijsian, NULL, NULL, InverseDeWijsian); 
  
 
  // next model is constant in space, possibly multivariate
  IncludeModel("constant", PosDefType, 0, 0, 1, kappaconstant, XONLY,
	       PREVMODEL_I,
		 //  PREVMODEL_D, PREVMODEL_I, 
		 //wegen Variogramm berechnung in stat. Fall
		 checkconstant, rangeconstant, PREF_ALL,
		 false, PARAM_DEP, INFDIM, (ext_bool) false, MON_SUB_DEP);
  kappanames("M", REALSXP);  
  addCov(constant, NULL, NULL);
  addCov(nonstatconstant);
  AddVariant(NegDefType, ISOTROPIC);
  //  AddVariant(TcfType, PREVMODEL_I);
  //  AddVariant(ShapeType, EARTH_ISOTROPIC);
  //  AddVariant(TcfType, EARTH_ISOTROPIC);
  //  AddVariant(ShapeType, SPHERICAL_ISOTROPIC);
  //  AddVariant(TcfType, SPHERICAL_ISOTROPIC);



  pref_type pfix={0, 0, 0,  0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //              CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("fixcov", PosDefType, 0, 1, 3, kappa_fix,
	       KERNEL, UNREDUCED, // zwingend bei RAW-Konstruktionen !!
	       checkfix, rangefix, pfix, 
	       INTERN_SHOW, PARAM_DEP, INFDIM-1, (ext_bool) false,
	       NOT_MONOTONE);
  subnames("norm");
  kappanames("M", LISTOF + REALSXP,
	     COVARIATE_X_NAME, VECSXP,
	     COVARIATE_RAW_NAME, INTSXP);
  change_sortof(COVARIATE_X, DONOTVERIFYPARAM);
  addCov(fixStat, NULL, NULL);
  addCov(fix);
  setptwise(pt_paramdep);
  addTypeFct(Typefix);
  setDI(NULL, allowedIfix, NULL); 
  AddVariant(PosDefType, ISOTROPIC);
  AddVariant(PosDefType, EARTH_ISOTROPIC);

  pref_type pcox={2, 0, 0,  0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //              CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("coxisham",  PosDefType, 1, 1, 3, kappa_cox, 
	       XONLY, SYMMETRIC,
	       checkcox, rangecox, pcox,
	       false, SCALAR, CoxMaxDim, (ext_bool) false, NOT_MONOTONE);
  kappanames("mu", REALSXP, "D", REALSXP, "beta", REALSXP);  
  addCov(cox, NULL, NULL);
  addTBM(initcox, spectralcox);
  nablahess(coxnabla, coxhess);

  IncludePrim("cubic",  TcfType, 0, XONLY, ISOTROPIC, 
	      checkOK, NULL, 3, (ext_bool) false, MONOTONE);
  addCov(cubic, Dcubic, ScaleOne);
	       
  pref_type pcurl= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("curlfree",  PosDefType, 1, 1, 1, kappadivcurl, XONLY,
	       CARTESIAN_COORD,
	       checkdivcurl, rangedivcurl, pcurl,
	       false, PARAM_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("which", INTSXP);  
  addCov(curl, NULL, NULL);
 
  pref_type plocal={5, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //            CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  CUTOFF =  
    IncludeModel("cutoff",  PosDefType, 1, 1,2, NULL, XONLY, ISOTROPIC,
		 check_co, range_co, plocal,
		 false, SUBMODEL_DEP, MAXCEDIM,  (ext_bool) true, MONOTONE);
  kappanames("diameter", REALSXP, "a", REALSXP);  
  addCov(co, NULL, NULL);
  addCallLocal(alternativeparam_co);
 
  //  warning("siehe Primitive.cc/biWM: cutoff funktioniert nicht bei MLE, vereinheitlichung mit natsc und verbesserung von biWM notwendig");
 

  IncludeModel("dagum",  PosDefType, 0, 0, 3, NULL, XONLY, ISOTROPIC,
	       checkdagum, rangedagum, PREF_ALL, false, 1, INFDIM, (ext_bool) false,
	       MON_PARAMETER);
  kappanames("beta", REALSXP, "gamma", REALSXP, INTERNAL_PARAM, REALSXP);
  addCov(dagum, Ddagum, Inversedagum);
  RandomShape(0, struct_failed, initdagum, do_failed, false, true, false);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);
  setptwise(pt_posdef);


  IncludePrim("dampedcosine",  PosDefType, 1, XONLY, ISOTROPIC,
	      checkdampedcosine, rangedampedcosine, PARAM_DEP,
	      (ext_bool) false, NOT_MONOTONE);
  nickname("dampedcos");
  kappanames("lambda", REALSXP);
  addCov(dampedcosine, Ddampedcosine, Inversedampedcosine);
  // addlogCov(logdampedcosine);

 
  pref_type pderivative= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("deriv", PosDefType, 1, 1, 1, kappaderivative,
	       XONLY, CARTESIAN_COORD,
	       checkderivative, rangederivative, pderivative, 
	       false, PARAM_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       NOT_MONOTONE);
  kappanames("which",INTSXP);
  addCov(derivative, NULL, NULL);
  
   
  pref_type pdewijsian = {2, 5, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                     CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("DeWijsian", VariogramType,  1, XONLY, ISOTROPIC,
  	      checkOK, rangedewijsian, pdewijsian,
  	      SCALAR, INFDIM, (ext_bool) false, MONOTONE);
  nickname("dewijsian");
  kappanames("alpha", REALSXP);
  addCov(dewijsian, Ddewijsian, DDdewijsian, D3dewijsian, D4dewijsian,
	 Inversedewijsian);
  addLocal(coinitdewijsian, NULL);

 

  pref_type pdiv= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //               CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("divfree", PosDefType, 1, 1, 1, kappadivcurl, XONLY,
	       CARTESIAN_COORD,
	       checkdivcurl, rangedivcurl, pdiv, 
	       false, PARAM_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       NOT_MONOTONE);
  kappanames("which", INTSXP);  
  addCov(diverge, NULL, NULL);



  pref_type pepsC = {2, 0, 0, 5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("epsC",  PosDefType, 0, 0, 3, NULL, XONLY, ISOTROPIC,
	       checkepsC, rangeepsC, pepsC, 
	       false, SCALAR, INFDIM, (ext_bool) false, NORMAL_MIXTURE);
  nickname("epscauchy");
  kappanames("alpha", REALSXP, "beta", REALSXP, "eps", REALSXP);
  addCov(epsC, DepsC, DDepsC, NULL, NULL);
  addlogCov(logepsC);

  IncludePrim("exponential", TcfType, 0, XONLY, ISOTROPIC,
	      checkexponential, NULL, INFDIM, (ext_bool) false, COMPLETELY_MON);
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
  
  // operator, replacing any covariance fct C by EXP(C) (componentwise)
  IncludeScalar("Exp", 
		PosDefType, 1, 1, 2, SUBMODEL_D, SUBMODEL_I, checkExp,
		rangeExp, PREF_ALL, SUBMODEL_DEP, (ext_bool) false,
		NOT_MONOTONE);
  nickname("exponential");
  kappanames("n", INTSXP, "standardised", INTSXP);
  addCov(Exp, DExp, DDExp, NULL, NULL);
   // setptwise(pt_paramdef);

  pref_type
    pextrgauss = {5, 0, 0,  0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //             CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("extremalgauss", TcfType, 1, 1, 0, NULL,
	       XONLY, SUBMODEL_I,
	       check_extrgauss, NULL, pextrgauss, false,
	       SCALAR, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, NOT_MONOTONE);
  nickname("schlather");
  addCov(extrgauss, NULL, NULL);
  
  //X

  IncludePrim("FD", PosDefType,  1, XONLY, ISOTROPIC,
	      checkOK, rangeFD, 1, (ext_bool) false, NOT_MONOTONE); 
  nickname("fractdiff");
  kappanames("a", REALSXP);
  addCov(FD, NULL, NULL);

  
   // same as RMcovariate, ausser dass RMcovariate interpoliert
    // und CONSTANT mehrere Saetze von covariaten erlaubt??!!

 IncludePrim("flatpower", VariogramType, 1, XONLY, ISOTROPIC, 
		checkoesting, rangeoesting, INFDIM, (ext_bool) false,
		BERNSTEIN); // todo BERNSTEIN
  kappanames("alpha", REALSXP);
  addCov(2, oesting, Doesting, DDoesting, NULL, NULL);
  RandomShape(0, initoesting, do_statiso);
  Taylor(-1, 2, RF_NA, 4, RF_NA, 6);
  TailTaylor(-1, RF_NA, 0, 0);

  BROWNIAN = 
    IncludePrim("fractalB", VariogramType, 1, XONLY, ISOTROPIC, 
		checkfractalBrownian, rangefractalBrownian, INFDIM, (ext_bool) false,
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
	      checkOK, rangefractGauss, 1, (ext_bool) false, NOT_MONOTONE);
  kappanames("alpha", REALSXP);
  addCov(fractGauss, NULL, NULL);



  pref_type pgauss= {2, 0, 0, 5, 5, 5, 5, 5, 0, 0, 5, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  GAUSS = // 61
    IncludePrim("gauss",  PosDefType, 0, XONLY, ISOTROPIC,
		checkOK, NULL, pgauss,
		SCALAR, INFDIM, (ext_bool) false, NORMAL_MIXTURE);
  addCov(Gauss, DGauss, DDGauss, D3Gauss, D4Gauss, InverseGauss);
  addlogCov(logGauss, NULL, nonstatLogInvGauss);
  addTBM(NULL, spectralGauss);
  RandomShape(INFTY, struct_Gauss, initGauss, do_Gauss, false, true, false);
  addGaussMixture(DrawMixGauss, LogMixDensGauss);
  Taylor(-1.0, 2.0);
  TailTaylor(1, 0, 1.0, 2.0);

  
  IncludePrim("genB", VariogramType, 2, XONLY, ISOTROPIC, 
	      checkOK, rangegenBrownian, INFDIM, (ext_bool) false, MONOTONE);
  nickname("genfbm");
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(genBrownian, NULL, NULL, InversegenBrownian); 
  addlogCov(loggenBrownian);


  pref_type pgenc = {2, 0, 0, 2, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("gencauchy", PosDefType, 2, XONLY, ISOTROPIC,
	      checkgeneralisedCauchy, rangegeneralisedCauchy, pgenc,
	      SCALAR, INFDIM, (ext_bool) false, MON_PARAMETER); // todo part is even
  // LAPLACE
  kappanames("alpha", REALSXP, "beta", REALSXP);
  addCov(generalisedCauchy, DgeneralisedCauchy, DDgeneralisedCauchy,
	 InversegeneralisedCauchy);
  addlogCov(loggeneralisedCauchy);
  addLocal(coinitgenCauchy, ieinitgenCauchy);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);
  setptwise(pt_posdef);


  // bivariate stable or bivariate exponnetial model
  IncludePrim("bicauchy",  PosDefType, 4, kappa_biCauchy, XONLY, ISOTROPIC,
           checkbiCauchy, rangebiCauchy, PREF_ALL,
           2, 3, (ext_bool) false, MONOTONE);
   kappanames("alpha", REALSXP,
             "beta", REALSXP,
             "s", REALSXP,  // lower triangle definition
         "rho", REALSXP  // lower triangl
             );
 addCov(biCauchy, DbiCauchy, DDbiCauchy, D3biCauchy, D4biCauchy, NULL);
  addLocal(coinitbiCauchy, NULL);

  //Y
 

  IncludePrim("gengneiting",  PosDefType, 2, XONLY, ISOTROPIC, 
	      checkgenGneiting, rangegenGneiting, INFDIM-1, (ext_bool) true,
	      MONOTONE); // GNEITING_MON ??
  // not INFDIM, also not normalscale mixture and alpha will be void
  kappanames("kappa", INTSXP, "mu", REALSXP);
  addCov(genGneiting, DgenGneiting, DDgenGneiting, //NULL, NULL,
	 ScaleOne, NULL);



  IncludePrim("gengneiting",  PosDefType, 2, XONLY, ISOTROPIC, 
	      checkgenGneiting, rangegenGneiting, INFDIM-1, (ext_bool) true,
	      MONOTONE); // GNEITING_MON ??
  // not INFDIM, also not normalscale mixture and alpha will be void
  kappanames("kappa", INTSXP, "mu", REALSXP);
  addCov(genGneiting, DgenGneiting, DDgenGneiting, //NULL, NULL,
	 ScaleOne, NULL);

  GNEITING_INTERN =
    IncludeModel("gengneit_intern", PosDefType, 0, 0, 2, NULL, 
		 XONLY, ISOTROPIC,
		 checkgenGneiting, rangegenGneiting, PREF_ALL, 
		 true, SCALAR, PARAM_DEP, (ext_bool) true, MONOTONE);  
  nickname("gengneiting");
  kappanames("kappa", INTSXP, "mu", REALSXP);
  addCov(Gneiting, DGneiting, DDGneiting, ScaleOne);


  IncludeScalar("gneiting", PosDefType, 0, 0, 1, XONLY, ISOTROPIC,
	      checkGneiting, rangeGneiting, PREF_ALL, 
	      PARAM_DEP, (ext_bool) true, MONOTONE);  // GNEITING_MON ??
  kappanames("orig", INTSXP);
  addCov(Gneiting, DGneiting, DDGneiting, ScaleOne);

  pref_type pfgennsst= { 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                 CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeScalar("gennsst",  PosDefType, 2, 2, 1, SUBMODEL_D,
	       SYMMETRIC, //SUBMODEL_I,
	       checkgennsst, rangegennsst, pfgennsst,
	       SUBMODEL_DEP, (ext_bool) false, NOT_MONOTONE);
  kappanames("dim_u", INTSXP);
  addCov(gennsst, NULL, NULL);
  addCov(nonstatgennsst);
  setDI(allowedDgennsst, NULL, NULL); //allowedIgennsst);
  
  //Y
  GENNSST_INTERN = // to do: why internal??
  IncludeModel("gennsst_intern", 
	       PosDefType, // stimmt nicht, aber egal
	       2, 2, 1, kappa_gennsst_intern, XONLY, SYMMETRIC,
	       checkgennsst_intern, range_gennsst_intern, pfgennsst,
	       true, SCALAR, SUBMODEL_DEP, (ext_bool) false, NOT_MONOTONE);
  nickname("gennsst");
  kappanames("A", REALSXP);
  addCov(gennsst_intern, NULL, NULL);

  /*
 pref_type phelmholtz= {2, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeModel("helmholtz",  PosDefType, 1, 1, 2, kappamixed, XONLY, SYMMETRIC,
         checkhelmholtz,rangeHelmholtz, phelmholtz,
           true, PARAM_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("component",REALSXP,"Aniso",REALSXP);
  addCov(helmholtz, NULL, NULL);
  */

  pref_type phyper= {2, 0, 0, 3, 0, 4, 5, 0, 5, 0, 5, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("hyperbolic",  PosDefType, 3, XONLY, ISOTROPIC,
	      checkhyperbolic, rangehyperbolic, phyper,
	      SCALAR, INFDIM, (ext_bool) false, NORMAL_MIXTURE);
  kappanames("nu", REALSXP, "lambda", REALSXP, "delta", REALSXP);
  addCov(hyperbolic, Dhyperbolic, NULL); // InversehyperbolicSq);
  addlogCov(loghyperbolic);
  RandomShape(0, struct_failed, inithyperbolic, do_failed, false, true, false);


  IncludePrim("iacocesare",  PosDefType, 3, XONLY, DOUBLEISOTROPIC, 
	      checkOK, rangeIacoCesare, INFDIM, (ext_bool) false, NOT_MONOTONE);
  nickname("iaco");
  kappanames("nu", REALSXP, "lambda", REALSXP, "delta", REALSXP);
  addCov(IacoCesare, NULL, NULL);
  

  IncludeModel("identity", ManifoldType, 1,1, 1, NULL, PREVMODEL_D, PREVMODEL_I,
	       checkId, rangeId, PREF_ALL, 
	       false, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       MON_SUB_DEP);
  nickname("idmodel");
  kappanames("vdim", INTSXP);
  addCov(IdStat, DId, DDId, IdInverse);
  addCov(IdNonStat);
  addTBM(TBM2Id, initId, spectralId);
  addLocal(coinitId, ieinitId);
  addTypeFct(TypeId);


  IncludePrim("kolmogorov",  VariogramType, 0, XONLY, VECTORISOTROPIC,
	      checkKolmogorov, NULL, 3, 3, (ext_bool) false, NOT_MONOTONE);
  addCov(Kolmogorov, NULL, NULL);

  pref_type plgd1= {2, 0, 0, 5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("lgd1",  PosDefType, 2, NULL, XONLY, ISOTROPIC, 
	      checklgd1, rangelgd1, plgd1, 
	      SCALAR, PARAM_DEP, (ext_bool) false, MONOTONE);
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
 
  IncludeScalar("mastein",  PosDefType, 1, 1, 2, XONLY, DOUBLEISOTROPIC, 
	       check_MaStein, range_MaStein, PREF_ALL, 
	       SUBMODEL_DEP, (ext_bool) false, NOT_MONOTONE);
  kappanames("nu", REALSXP, "delta", REALSXP);
  addCov(MaStein, NULL, NULL);

    
  IncludeScalar("ma1", PosDefType,  1, 1, 2, XONLY, SYMMETRIC,
	       checkma1, rangema1, PREF_ALL, 
	       SUBMODEL_DEP, (ext_bool) false, NOT_MONOTONE);
  nickname("ma");
  kappanames("alpha", REALSXP, "theta", REALSXP);
  addCov(ma1, NULL, NULL);


  IncludeScalar("ma2",  PosDefType, 1, 1, 0, XONLY, SYMMETRIC, checkma2, NULL,
		PREF_ALL, SUBMODEL_DEP, (ext_bool) false, NOT_MONOTONE);
  nickname("intexp");
  addCov(ma2, NULL, NULL);

  //X
   pref_type pmatrix = {3, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 3, 5};
  //                  CE CO CI TBM Sp di sq Tr av n mpp Hy spf any
  MATRIX=
    IncludeModel("M",  ManifoldType, 1, MAXSUB, 2, kappaM, SUBMODEL_D,
		 SUBMODEL_I,
		 checkM, rangeM, pmatrix, false, PARAM_DEP, SUBMODEL_DEP,
		 (ext_bool) SUBMODEL_DEP, NOT_MONOTONE);
  nickname("matrix");
  kappanames("M", REALSXP, "vdim", INTSXP); // vdim ist das interne!
  change_typeof(M_M, ShapeType);
  addCov(Mstat, NULL, NULL);
  addCov(Mnonstat);
  add_sortof(sortof_M);
  addTypeFct(TypeM);
  RandomShape(0, struct_failed, initM, do_failed, false, true, false);
  setDI(allowedDM, allowedIM, NULL);


  IncludeScalar("matern", PosDefType, 0, 0, 2, PARAMDEP_D, PARAMDEP_I,
	       checkMatern, rangeWM, PREF_ALL, INFDIM, (ext_bool) false,
	       MON_SUB_DEP);
  kappanames("nu", REALSXP, "notinvnu", INTSXP);
  change_sortof(WM_NU, CRITICALPARAM); 
  change_typeof(WM_NU, RandomOrShapeType);
  change_sortof(WM_NOTINV, ONLYRETURN);
  addCov(Matern, DMatern, DDMatern, D3Matern, D4Matern, InverseMatern);
  addCov(NonStMatern);
  addlogCov(logMatern, logNonStMatern, NULL);
  addTBM(initMatern, spectralMatern);
  addLocal(coinitWM, ieinitWM);
  setptwise(pt_posdef);
  setDI(allowedDWM, allowedIWM, setWM);
  addTypeFct(TypeWM);
  RandomShape(0, struct_failed, initWM, do_failed, false, true, false);
   
  //  addGaussMixture(DrawMixWM, LogMixDensWM);

  IncludeModel("mqam", PosDefType,
  	       2, 10, 1, kappamqam, XONLY, SYMMETRIC,
  	       checkmqam, rangemqam, PREF_ALL, 
	       false, PARAM_DEP, SUBMODEL_DEP, (ext_bool) false, NOT_MONOTONE);
  subnames("phi");
  kappanames("theta", REALSXP);
  change_sortof(QAM_THETA, CRITICALPARAM); 
  addCov(mqam, NULL, NULL);

  
  pref_type
    pmultiquad = {0, 0, 0,  0, 0, 5, 0, 0, 0, 0,  0, 0,  0,  5};
  //        CE CO CI TBM Sp di sq Ma av  n mpp Hy spf any
  IncludePrim("multiquad", PosDefType, 2, NULL,
	      XONLY, SPHERICAL_ISOTROPIC, checkOK, rangemultiquad,
	      pmultiquad, SCALAR, 2, (ext_bool) false, MONOTONE);
  kappanames("delta", REALSXP, "tau", REALSXP);
  addCov(multiquad, NULL, NULL);


  NATSC_INTERN =NATSC_USER = 
    IncludeModel("natsc", PosDefType,  1, 1, 0, NULL, XONLY, ISOTROPIC,
		 checknatsc, NULL, PREF_ALL,
		 false, 1, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  // nie einen Parameter !
  addCov(natsc, Dnatsc, DDnatsc, Inversenatsc);
  addLocal(coinitnatsc, ieinitnatsc);
  addTBM(tbm2natsc, initnatsc, spectralnatsc);
  setptwise(pt_submodeldep);
  AddVariant(TcfType, ISOTROPIC);

  // NATSC_INTERN = CopyModel("natsc_intern", NATSC_USER);
  // make_internal();



  pref_type pfnsst= { 4, 0, 0, 2, 0, 5, 4, 0, 0, 0, 0, 0, 0, 5};
  //                 CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludeScalar("nsst",  PosDefType, 2, 2, 1, XONLY, DOUBLEISOTROPIC,
	       checknsst, rangensst, pfnsst,
	       SUBMODEL_DEP, (ext_bool) false, NOT_MONOTONE);
  subnames("phi", "psi");
  kappanames("delta", REALSXP); 
  change_sortof(NSST_DELTA, CRITICALPARAM);
  addCov(nsst, Dnsst, NULL);
  addTBM(TBM2nsst);
  setptwise(pt_posdef);

  //  IncludePrim("nsst2", 7, checknsst2, DOUBLEISOTROPIC, 
  //		   rangensst2);
  //  addCov(nsst2, Dnsst2, NULL);
  //  addTBM(NULL, NULL /* TBM3nsst2 */);

 
  IncludeModel("parsWM", PosDefType, 0, 0, 1, kappa_parsWM, 
	       XONLY, ISOTROPIC,
	       checkparsWM, rangeparsWM, PREF_ALL,
	       false, PARAM_DEP, INFDIM, (ext_bool) false, NOT_MONOTONE);
  nickname("parswm");
  addCov(parsWM, parsWMD, NULL);
  kappanames("nudiag", REALSXP);
  change_sortof(PARSnudiag, CRITICALPARAM);
  setptwise(pt_posdef);
  RandomShape(0, struct_failed, initparsWM, do_failed, false, true, false);


  IncludePrim("penta", PosDefType, 0, XONLY, ISOTROPIC,
	      checkOK, NULL, 3, (ext_bool) true, MONOTONE);
  addCov(penta, Dpenta, ScaleOne);

  IncludePrim("power", PosDefType,  1, XONLY, ISOTROPIC, 
	      checkpower, rangepower, INFDIM-1, (ext_bool) true, MONOTONE);
  nickname("askey");
  kappanames("alpha", REALSXP);
  addCov(power, Dpower, ScaleOne);	
  AddVariant(TcfType, ISOTROPIC);
  //  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);

  // FORMERLY DIFFERENR DEFs OF PARAMs
  POW = IncludeScalar("Pow", ShapeType, 1, 1, 1, XONLY, SUBMODEL_I,
		     checkPow, rangePow, PREF_ALL, SUBMODEL_DEP,
		     Submodeldep, NOT_MONOTONE);
  nickname("power");
  addCov(Pow, DPow, DDPow, InversePow); 
  kappanames("alpha", REALSXP);
  setptwise(pt_posdef);
  AddVariant(NegDefType, SUBMODEL_I);
  AddVariant(PosDefType, SUBMODEL_I);
  AddVariant(TcfType, ISOTROPIC);
  RandomShape(0, struct_failed, initPow, do_failed, false, false, false);
  
 
   pref_type pfprod= { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
   PROD = 
     IncludeModel("prod",  PosDefType, 1, 1, 0, NULL,
		  KERNEL, UNREDUCED, 
		  checkprod, NULL, pfprod,
	       false, SUBMODEL_DEP, INFDIM-1, (ext_bool) false, NOT_MONOTONE);
  addCov(nonstatprod);
  subnames("phi");  
  //  AddVariant(PosDefType, SPHERICAL_COORD);
  // AddVariant(PosDefType, EARTH_COORD);
  setptwise(pt_posdef);
  //  printf("nick=%.50s\n",  DefList [currentNrCov - 1].nick);


  IncludeModel("qam",  PosDefType, 2, MAXSUB, 1, kappaqam, XONLY, ISOTROPIC,
	       checkqam, rangeqam, PREF_ALL, 
	       false, SCALAR, SUBMODEL_DEP, (ext_bool) false, NOT_MONOTONE);
  kappanames("theta", REALSXP);
  subnames("phi");
  addCov(qam, NULL, NULL);
  change_sortof(QAM_THETA, CRITICALPARAM);
  

  IncludePrim("qexponential",  PosDefType, 1, XONLY, ISOTROPIC, 
	      checkOK, rangeqexponential, INFDIM, (ext_bool) false, NOT_MONOTONE);
  nickname("qexp");
  kappanames("alpha", REALSXP);
  addCov(qexponential, Dqexponential, Inverseqexponential);


  pref_type pscale = {  0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                   CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  SCALEMODEL =
    IncludeModel("scale", PosDefType, 2, 3, 0, NULL, KERNEL, SYMMETRIC,
		 checkscale, NULL, pscale, false, SCALAR, INFDIM - 1, falsch,
		 NOT_MONOTONE);
  subnames("phi", "scaling", "penalty");
  addCov(nonstatscale); 

  
  IncludeModel("schur",  PosDefType, 1, 1, 3, kappaSchur, 
	       PREVMODEL_D, PREVMODEL_I, 
	       checkSchur, rangeSchur, PREF_ALL,
	       false, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, NOT_MONOTONE);
  kappanames("M", REALSXP, "diag", REALSXP, "rhored", REALSXP);
  addCov(Schurstat, DSchur, D2Schur, D3Schur, D4Schur, NULL);
  addCov(Schurnonstat);
  add_sortof(sortof_M); 

 
  IncludeModel("shift", PosDefType, 1, 1, 1, kappashift, XONLY, CARTESIAN_COORD,
	       checkshift, rangeshift, PREF_ALL, 
	       false, PARAM_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
	       NOT_MONOTONE);
  nickname("delay"); // delayeffect
  addCov(shift, NULL, NULL);
  kappanames("s", REALSXP);

  pref_type
    psinepower = {0, 0, 0,  0, 0, 5, 0, 0, 0, 0,  0, 0,  0,  5};
  //             CE CO CI TBM Sp di sq Ma av  n mpp Hy spf any
  IncludePrim("sinepower", PosDefType, 1, NULL,
	      XONLY, SPHERICAL_ISOTROPIC, checkOK, rangesinepower,
	      psinepower, SCALAR, 2, (ext_bool) false, MONOTONE);
  kappanames("alpha", REALSXP);
  addCov(sinepower, NULL, NULL);


  IncludePrim("spherical", TcfType, 0, NULL, XONLY, ISOTROPIC, 
	      checkOK, NULL, 3, (ext_bool) true, GNEITING_MON);
  nickname("spheric");
  addCov(spherical, Dspherical, DDspherical, ScaleOne);
  addTBM(TBM2spherical);
  RandomShape(1, structspherical, initspherical, dospherical,
	      false, true, false);
  Taylor(-3.0, 1.0, 0.5, 3.0);


  IncludePrim("stable",  PosDefType, 1, XONLY, ISOTROPIC, 
	      checkstable, rangestable, INFDIM, (ext_bool) false,
	      MON_PARAMETER);
  kappanames("alpha", REALSXP);
  addCov(stable, Dstable, DDstable, Inversestable);
  addlogCov(logstable, NULL, nonstatLogInversestable);
  addLocal(coinitstable, ieinitstable);
  AddVariant(TcfType, ISOTROPIC);
  AddVariant(PosDefType, SPHERICAL_ISOTROPIC);
  setptwise(pt_posdef);

  // DOUBLEISOTROPIC variant of stable -- used for testing purposes only
  //  IncludePrim("stableX", 1, checkOK, DOUBLEISOTROPIC, 
  //		  rangestable);
  //  addCov(stableX, DstableX, Inversestable);
  //  addTBM(NULL, NULL)


  STEIN =  
    IncludeModel("Stein", PosDefType,  1, 1, 2, NULL, XONLY, ISOTROPIC,
		 check_Stein, range_Stein, plocal,
		 false, SCALAR, MAXCEDIM, (ext_bool) true, NOT_MONOTONE);
  nickname(METHOD_NAMES[CircEmbedIntrinsic]);
  kappanames("diameter", REALSXP, "rawR", REALSXP);  
  change_sortof(pLOC_DIAM, FORBIDDENPARAM);
  change_sortof(pLOC_A, FORBIDDENPARAM);
  addCov(Stein, NULL, NULL);
  addCallLocal(alternativeparam_Stein);
  //  RandomShape(struct_ce_approx, init_ce_approx, do_ce_approx);

 
  IncludePrim("steinst1",  PosDefType, 2, kappaSteinST1, XONLY, SYMMETRIC,
	      checkSteinST1, rangeSteinST1, INFDIM, (ext_bool) false, NOT_MONOTONE);
  nickname("stein");
  kappanames("nu", REALSXP, "z", REALSXP);
  addCov(SteinST1, NULL, NULL);
  addTBM(initSteinST1, spectralSteinST1);
  RandomShape(0, struct_failed, initSteinST1, do_failed, false, true, false);


  IncludeModel("stp", PosDefType, 1, 2, 3, kappa_stp, KERNEL, SYMMETRIC,
	       checkstp, rangestp, PREF_ALL,
	       false, SCALAR, StpMaxDim, (ext_bool) false, NOT_MONOTONE);
  addCov(stp);
  kappanames("S", REALSXP, "z", REALSXP, "M", REALSXP);
  change_typeof(STP_S, ShapeType);
  RandomShape(structStp, true);
  subnames("xi", "phi"); // H ueberall wo U-x steht. dort U-H(x)
  //                           gedoppelte immer zum Schluss!
  
 
  SHAPESTP = 
    IncludeModel("shape.stp",  ShapeType, 1, 4, 3, kappa_stp, KERNEL, 
		 CARTESIAN_COORD, check_shapestp, rangestp, PREF_NOTHING,
		 true, SCALAR, StpMaxDim, (ext_bool) true, NOT_MONOTONE);
  kappanames("S", REALSXP, "z", REALSXP, "M", REALSXP); 
  addlogCov(logshapestp);
  subnames("xi2", "phi", "SXX", "gauss"); // hier gedoppeltes S nicht am Schluss 
  //                                       da auf checkstp zugreifend
  RandomShape(0, init_shapestp, do_shapestp);
  
  TBM_OP = // old RandomFields tbm2, tbm3
    IncludeModel("tbm", ManifoldType, 1, 1, 3, NULL, XONLY, PARAMDEP_I,
		 checktbmop, rangetbmop, PREF_ALL, false,
		 SUBMODEL_DEP, PARAM_DEP, (ext_bool) PARAM_DEP, NOT_MONOTONE);
  kappanames("fulldim", INTSXP, "reduceddim", INTSXP, "layers", INTSXP); 
  change_sortof(TBMOP_LAYERS, ONLYRETURN); // NA will not be estimated
  addCov(tbm, NULL, NULL); // Dtbm, NULL); todo
  setDI(NULL, allowedItbm, settbm);
  addTypeFct(Typetbm);


  //  { int i; for (i=0; i<=Nothing; i++) printf("%d %d\n ", i, DefList[TREND].pref[i]); assert(false); }

   pref_type pfsum= { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
   IncludeModel("sum",  NegDefType, 0, 1, 0, NULL,
		KERNEL, UNREDUCED,  checksum, NULL, pfsum,
	       false, SUBMODEL_DEP, INFDIM-1, (ext_bool) false, NOT_MONOTONE);
  addCov(nonstatsum);
  subnames("phi");  


  USER =
    IncludeModel("U", ManifoldType, 0, 0, 16, kappaUser, 
		 PARAMDEP_D, PARAMDEP_I,
		 checkUser, rangeUser, PREF_AUX, 
		 true,// FREEVARIABLE vorhanden. Muss extra in SpecialRMmodel.R
		 // definiert und nicht ueber generatemodels.R
		 PARAM_DEP, INFDIM, (ext_bool) false, // per default.
		 NOT_MONOTONE);
  nickname("user");
  kappanames("type", INTSXP, "domain", INTSXP,  "isotropy", INTSXP,
	     "vdim", INTSXP, "beta", REALSXP, "coordnames", INTSXP,
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
  setDI(allowedDuser, allowedIuser, setUser);
  addTypeFct(TypeUser);


  /*
    IncludeModel("gamma2cov",  PosDefType, 1, 1, 1, 
	       kappaGammaToCov, KERNEL, SYMMETRIC,
	       checkGammaToCov, rangeGammaToCov, pg2c, 
		false, SCALAR, SUBMODEL_DEP, (ext_bool) false, NOT_MONOTONE);
  kappanames("x0", REALSXP);
  addCov(GammaToCov);
  */


  
  pref_type pg2c = {  0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
        //           CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  VARIOGRAM2COV = // intern ok  
    IncludeModel("cov", PosDefType, 1, 1, 2, kappavariogram2cov,
		 KERNEL, CARTESIAN_COORD,
		 checkvariogram2cov, rangevariogram2cov, pg2c, 
		 INTERN_SHOW, SCALAR, INFDIM-1, (ext_bool) false, NOT_MONOTONE);
  subnames("gamma");
  kappanames("x", VECSXP, "a", REALSXP);
  // change_typeof(VAR2COV_X, MixedInputType, RMCOV_X); 
  change_sortof(VAR2COV_X, DONOTVERIFYPARAM);
  addCov(variogram2cov);
  
  VECTOR = 
    IncludeModel("vector",  PosDefType, 1, 1, 2, NULL, XONLY, CARTESIAN_COORD,
		 checkvector, rangevector, PREF_ALL, 
		 false, PARAM_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		 NOT_MONOTONE);
  addCov(vector, NULL, NULL);
  kappanames("a", REALSXP, "Dspace", INTSXP);
  addFurtherCov(vectorAniso, NULL); 


  pref_type pwave = {2, 0, 0, 0, 5, 4, 5, 0, 0, 0, 0, 0, 0, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("wave", PosDefType, 0, XONLY, ISOTROPIC, 
	      checkOK, NULL, pwave, SCALAR, 3, (ext_bool) false, NOT_MONOTONE);
  addCov(wave, NULL, Inversewave);
  addTBM(initwave, spectralwave);

 
  IncludeScalar("whittle", PosDefType, 0,0, 2, PARAMDEP_D, PARAMDEP_I,  
	       checkWM, rangeWM, PREF_ALL, INFDIM, (ext_bool) false,
	       NORMAL_MIXTURE);
  kappanames("nu", REALSXP, "notinvnu", INTSXP);
  change_typeof(WM_NU, RandomOrShapeType);
  change_sortof(WM_NU, CRITICALPARAM); 
  addCov(Whittle, DWhittle, DDWhittle, D3Whittle, D4Whittle, InverseWhittle);
  addCov(NonStWhittle);
  addlogCov(logWhittle, logNonStWhittle, NULL);
  addTBM(initWhittle, spectralWhittle);
  addLocal(coinitWM, ieinitWM);
  addGaussMixture(DrawMixWM, LogMixDensW);
  setptwise(pt_posdef);
  setDI(allowedDWM, allowedIWM, setWM);
  addTypeFct(TypeWM);
  RandomShape(0, struct_failed, initWM, do_failed, false, true, false);


  //////////////////////////////////////////////////////////////////////
  
  for (int nr=0; nr<currentNrCov; nr++) {     
    defn *C = DefList + nr;    
    //    printf("\n %.50s\t", C->name);
    for (int i=0; i<Nothing; i++) {
    //      int old = C->pref[i];
      C->pref[i] *= C->implemented[i]==IMPLEMENTED;
      //      if (old != C->pref[i] && i!=Nugget) printf("%.50s:%d%d ", METHOD_NAMES[i], old, C->pref[i]); // && i!=Specific
    }
    C->pref[Nothing] *= (C->cov != ErrCov || C->nonstat_cov != ErrCovNonstat);
  }
    
  //////////////////////////////////////////////////////////////////////
  
  pref_type pfnugget= { 1, 0, 0, 0, 0, 1, 1, 0, 0, 5, 0, 0, 0, 5};
  //                  CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  NUGGET  = 
    IncludeModel("nugget", TcfType, 0, 0, 2, NULL,
		 // XONLY, ISOTROPIC,
		 PARAMDEP_D, PARAMDEP_I,
		 check_nugget, range_nugget, pfnugget, 
		 false, PREVMODEL_DEP, INFDIM, (ext_bool) true, MONOTONE);
  kappanames("tol", REALSXP, "vdim", INTSXP);
  change_sortof(NUGGET_TOL, FORBIDDENPARAM);
  change_sortof(NUGGET_VDIM, FORBIDDENPARAM);
  addCov(nugget, NULL, Inversenugget);
  addCov(nuggetnonstat);
  addReturns(covmatrix_nugget, iscovmatrix_nugget);
  //  AddVariant(TcfType, EARTH_ISOTROPIC);
  //AddVariant(TcfType, SPHERICAL_ISOTROPIC);
  setDI(allowedDnugget, allowedInugget, setnugget);
  addTypeFct(Typenugget); 
}
