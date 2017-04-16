

/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 -- 2017 Martin Schlather

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




#ifndef RFsimu_public_H
#define RFsimu_public_H 1
#include "basic.h"

// #include <f2c.h> /* otherwise false/true not define or conflicting... */


#define NATSCALE_EXACT 1 /* or approx or mle */
#define NATSCALE_ORNUMERIC 2
//#define NATSCALE_APPROX 2
#define NATSCALE_MLE 3 /* check fitvario when changing; best not to change !!! */


#ifdef __cplusplus
extern "C" {
#endif 
  void GetCurrentNrOfModels(int *init, int *nr);
  void GetModelNr(char **name, int *n); /* transforms string "covfct" 
					   into covnr of CovList */
  void GetModelName(int *nr, char **name, char **nick); /* transforms covnr into string */
  void GetAttr(int *nr, int *type, int *op, int *monotone, int *finiterange, 
	       int * simpleArguments,
	       int *internal, int *stat, int *iso, int *maxdim, int *vdim,
	       int *includevariants, int *paramtype, int*n);
  void GetNrParameters(int *covnr, int* kappas);
  void PrintModelList(int *local, int *operat, int *nick); // used in Roger's book
  //  void GetModelList(int* idx, int*internal); // used in generate.R
 
  /* Natural Scaling -- in contrast to the RFparameters(PracticalRange)
     GENERAL_NATURALSCALING may take the following values (furthermore,
     0: no rescaling, +10 numerical determination allowed
  */
  
  //void GetNaturalScaling(int *covnr, double *q,  /* KAPPAS only  !! */
  //			      int *naturalscaling, double *natscale, int *error);
  
  void ResetWarnings(int *all);

  void PutValuesAtNA(int *reg, double *values);
  void PutValuesAtNAnoInit(int *reg, double *values);
  
  void expliciteDollarMLE(int * modelnr, double *values);
  
 
  void GetModelRegister(char **name, int* nr);
  
  void MultiDimRange(int *model_nr, int *set, double *natscale);
 
  void NoCurrentRegister();
  void GetCurrentRegister(int *reg);
 void PutGlblVar(int *reg, double *var);

  void DeleteKey(int *reg);
  

  void attachRFoptionsRandomFields();
  void detachRFoptionsRandomFields();
  void RelaxUnknownRFoption(int *);


  
  //  void GetKeyInfo(int *keyNr, int *total, int *lengths, int *dim, 
  //		  int *timespacedim,int *grid,int *type,int *maxdim,int *vdim);
  
  SEXP GetParameterNames(SEXP nr);
  SEXP GetSubNames(SEXP nr);
  SEXP GetAllModelNames();
  SEXP GetCoordSystem(SEXP keynr, SEXP oldsystem, SEXP newsystem);
  SEXP GetExtModelInfo(SEXP keynr, SEXP level, SEXP spconform, SEXP whichSub);
  SEXP GetModel(SEXP keynr, SEXP modus, SEXP spconform, SEXP whichSub,
		SEXP SolveRandom, SEXP do_notreturnparam);
  
  /* 
     check with InitSimulateRF in case of any changes !!
     both InitSimulateRF and DoGauss
  */

  // PROCESSESNAMES
  SEXP Init(SEXP model_reg, SEXP model, SEXP x, SEXP NA_OK);
  
  SEXP EvaluateModel(SEXP X, SEXP Covnr);
  
  SEXP GetProcessType(SEXP model_reg, SEXP model);
  // storing : current value
  // printlevel : current value
  // naturalscaling : fixed value
  
  /* GENERAL_STORING: 
     true: intermediate results are stored: might be rather memory consuming,
     but the simulation can (depending on the method chosen) be much faster
     when done the second time with exactly the same parameters
     do not forget to call DeleteAllKeys when all the simulation are done.

     false: intermediate results are not stored if SimulateRF is called or
     stored only until DoGauss is called.
     minimum memory consumption, but might be rather slow if many
     simulations are performed with exactly the same parameters
     to be safe call DeleteAllKeys when all the simulations are done
     Warning! init_foo may depend on GENERAL_STORING! This may cause problems
     if GENERAL_STORING is changed between init_foo and do_foo from false to
     true and if do_foo is called more than once! (SimulateRF is safe in
     this case)
     
     indifferent behaviour for the simulation methods if parameters are 
     changed after each simulation.
     In case, only MAXKEYS different parameter sets are used, but after 
     each step the parameter set is changed, use different keynr for each 
     parametere set, and STORING==true, to get fast simulation 
     PL:
     0: no messages; 
     1: error messages; 
     2: hints when algorithm is at a forcation; 
     >3: more and more debugging information
     GENERAL_NATURALSCALING: (boolean)
     is called PracticalRange in R  (see Chiles&Delfiner?!) 
  */


  // for isotropic spatial data
  SEXP empiricalvariogram(SEXP  X, SEXP Dim, SEXP Lx, 
			SEXP Values, SEXP Repet, SEXP Grid, 
		        SEXP Bin, SEXP Nbin, SEXP Vdim, SEXP Method);

  SEXP empvarioXT(SEXP Xsexp, SEXP Tsexp, 
		  SEXP Lx, 
		  SEXP Values, SEXP Repet, SEXP Grid,
		  SEXP Bin, SEXP Nbin, 
		  SEXP Phi,    // vector of a real and an integer
		  SEXP Theta,  // vector of a real and an integer
		  SEXP DT,   // in grid units, smallest positive value, for  
		  //            the temporal lag of the emp. variogram
		  //            stepT * nstepT is the greatest time span
		  SEXP Vdim, SEXP Method);
  SEXP fftVario3D(SEXP coord,
		  SEXP sumvals,
		  SEXP nbvals,
		  SEXP bin,
		  SEXP nbin,
		  SEXP lenT,
		  SEXP stepT,
		  SEXP nstepT,
		  SEXP phi,
		  SEXP theta,
		  SEXP repet,
		  SEXP vdim,
		  SEXP segmentEmpVario,
		  SEXP pseudo);
  
  
  // fractal dimension, 
  SEXP boxcounting(SEXP Z, SEXP LX,  SEXP Repet, SEXP Factor, SEXP Eps);
  SEXP detrendedfluc(SEXP dat, SEXP lx, SEXP repet, SEXP boxes, SEXP ldfa);
  SEXP periodogram(SEXP dat, SEXP len, SEXP repet, SEXP fftm, 
		   SEXP part, SEXP shift);
  SEXP minmax(SEXP dat, SEXP dim, SEXP ldim, SEXP boxes, SEXP lb);
  

  //  void Cov(int *reg, double *result);
  //  void CovMatrix(int *reg, double *result);
  // void Variogram(int *reg, double *value);
  // void Pseudovariogram(int* reg, double *result);
  //  SEXP CovMatrixIntern(SEXP reg, SEXP x, SEXP dist, SEXP grid,
  //		       SEXP lx, SEXP result, SEXP nonzeros);

  SEXP VariogramIntern(SEXP reg);
  //  SEXP CovMatrixLoc(SEXP reg, SEXP x, SEXP dist, SEXP xdim, SEXP lx,
  //		    SEXP result, SEXP nonzeros);
  SEXP CovLoc(SEXP reg, SEXP x, SEXP y, SEXP xdim, SEXP lx, SEXP result);
	      

  //  SEXP CovMatrixSelectedLoc(SEXP reg, SEXP x, SEXP dist, SEXP xdim, SEXP lx,
  //			    SEXP selected, SEXP nsel, SEXP result,
  //			    SEXP nonzeros);
  //
  //  SEXP CovMatrixSelected(SEXP reg, SEXP selected, SEXP nsel, SEXP result,
  //			 SEXP nonzeros);


  
  SEXP SetAndGetModelInfo(SEXP model_reg, SEXP model, 
			  SEXP spatialdim, SEXP distances, 
			  SEXP ygiven, // TRUE is the standard variant
			  SEXP Time, SEXP xdim, SEXP shortlen,
			  SEXP allowforintegerNA, SEXP excludetrend);
  
 SEXP SetAndGetModelLikeli(SEXP model_reg, SEXP model, SEXP x);

  SEXP Take2ndAtNaOf1st(SEXP model_reg, SEXP model, SEXP model_bound,
			SEXP spatialdim, SEXP Time, SEXP xdim, 
			SEXP nbounds, SEXP skipchecks);
    
 
  SEXP countelements(SEXP idx, SEXP N, SEXP Totparts);
  SEXP countneighbours(SEXP Xdim, SEXP parts, SEXP Squarelength, 
		       SEXP cumgridlen, SEXP boxes);
  SEXP getelements(SEXP Idx, SEXP Xdim, SEXP N, SEXP Cumgridlen, SEXP Boxes);
  SEXP getneighbours(SEXP Xdim, SEXP Parts, SEXP Squarelength, SEXP Cumgridlen,
		     SEXP Neighbours);

  SEXP set_boxcox(SEXP boxcox);
  SEXP get_boxcox();
  //  SEXP BoxCox_inverse(SEXP boxcox, SEXP res);
  SEXP BoxCox_trafo(SEXP boxcox, SEXP res, SEXP Vdim, SEXP inverse);	

  SEXP get_logli_residuals(SEXP Reg);
  SEXP get_likeliinfo(SEXP model_reg);
  SEXP simple_residuals(SEXP model_reg);
  SEXP get_linearpart(SEXP model_reg, SEXP Set);
  
#ifdef __cplusplus
}
#endif



#endif /* RF_simu_PUBLIC_H*/










