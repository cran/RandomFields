#ifndef RFsimu_public_H
#define RFsimu_public_H 1
#include "basic.h"

// #include <f2c.h> /* otherwise false/true not define or conflicting... */

// 
#define DoublePrecision 1
#ifdef DoublePrecision
typedef double res_type;
#else
typedef float res_type;
#endif


#define NATSCALE_EXACT 1 /* or approx or mle */
#define NATSCALE_ORNUMERIC 2
//#define NATSCALE_APPROX 2
#define NATSCALE_MLE 3 /* check fitvario when changing; best not to change !!! */


extern "C" {

  void getCov(int *n, double *COV);
  
  void test();
  
  void niceFFTnumber(int *n);
  
  void RNGtest();
  
  void GetrfParameters(int *covmaxchar, int *methodmaxchar,
		       int *typemaxchar,
		       int *covnr, int *methodnr, int *typenr,
		       int *maxdim, int *maxmodels);
  void GetrfParametersI(int *covmaxchar, int *methodmaxchar,
			int *typemaxchar,
			int *covnr, int *methodnr, int *typenr,
			int *maxdim, int *maxmodels);
  void GetModelNr(char **name, int *n); /* transforms string "covfct" 
					   into covnr of CovList */
  void GetModelName(int *nr, char **name, char **nick); /* transforms covnr into string */
  SEXP GetParameterNames(SEXP nr);
  SEXP GetSubNames(SEXP nr);
  void GetAttr(int *type, int *op, int *monotone, int *finiterange, 
	       int *internal, int *stat, int *iso, int *maxdim, int *vdim);
  void GetRange(int *nr, int *dim, int *index, double *params, int *lparam,
		double *range, int *lrange);
  void GetNrParameters(int *covnr, int* kappas);
  void PrintModelList(int *local, int *operat, int *nick);
  void GetModelList(int* idx, int*internal);
  void PrintMethods();
  void GetMethodName(int *nr, char **name);
  SEXP GetAllModelNames();
  void GetMethodNr(char **name, int *nr);
  //void GetGridSize(double*,double *,double *,int *,int *,int *,int *); 

  /* Natural Scaling -- in contrast to the RFparameters(PracticalRange)
     GENERAL_NATURALSCALING may take the following values (furthermore,
     0: no rescaling, +10 numerical determination allowed
  */
  
  //void GetNaturalScaling(int *covnr, double *q,  /* KAPPAS only  !! */
  //			      int *naturalscaling, double *natscale, int *error);
  
  SEXP RFoptions(SEXP options);
  void RelaxUnknownRFoption(int *relax);
  void ResetWarnings();
  
  //  void GetKeyInfo(int *keyNr, int *total, int *lengths, int *dim, 
  //		  int *timespacedim,int *grid,int *type,int *maxdim,int *vdim);
  
  SEXP GetExtModelInfo(SEXP keynr, SEXP level, SEXP spconform, SEXP whichSub);
  SEXP GetModel(SEXP keynr, SEXP modus, SEXP spconform,
		SEXP do_notreturnparam);
  
  /* 
     check with InitSimulateRF in case of any changes !!
     both InitSimulateRF and DoGauss
  */

  // PROCESSESNAMES
  SEXP Init(SEXP model_reg, SEXP model, SEXP x, SEXP y, SEXP T, 
	    SEXP spatialdim, SEXP grid, SEXP distances, SEXP time,
	    SEXP NA_OK);
  
  SEXP EvaluateModel(SEXP X, SEXP Covnr);
  
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
  void empiricalvariogram(double *x, int *dim, int *lx,
			  double *values, int *repet,
			  int *grid, double *bin, int *nbin, int *charact, 
			  double *res, double *sd, int *n);

  void empvarioXT(double *X, double *T, 
		  int *lx, 
		  double *values, int *repet, int *grid,
		  double *bin, int *nbin, 
		  double *phi,    // vector of a real and an integer
		  double *theta,  // vector of a real and an integer
		  int *dT,   // 
		  double *sum,   // \sum (a-b)^2 / 2
		  double *sq,   // \sum (a-b)^4 / 4
		  int *n);
  void fftVario3D(double *coord,
		  double *sumvals,
		  double *nbvals,
		  double *bin,
		  int *nbin,
		  int *lenT,
		  int *stepT,
		  int *nstepT,
		  double *phi,
		  double *theta,
		  int *repet,
		  double *empvario,	
		  double *n,
		  bool *pseudo);
  
  // kriging methods
  SEXP simpleKriging(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Invcov, SEXP Notna,
		     SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig);
  SEXP simpleKriging2(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Data,
		      SEXP Invcov, SEXP Notna,
		      SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
		      SEXP Sigma2);
  SEXP ordinaryKriging(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Invcov, SEXP Notna,
		       SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig);
  SEXP ordinaryKriging2(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Data,
			SEXP Invcov, SEXP Notna,
			SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
			SEXP Sigma2);
  SEXP universalKriging(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Invcov, SEXP Notna,
			SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
			SEXP Nfct, SEXP trend_expr, SEXP trend_envir);
  SEXP universalKriging2(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Data,
			 SEXP Invcov, SEXP Notna,
			 SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
			 SEXP Sigma2,  SEXP Nfct, SEXP trend_expr,
			 SEXP trend_envir);
  SEXP intrinsicKriging(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Invcov, SEXP Notna,
			SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig,
			SEXP Polydeg);
  SEXP intrinsicKriging2(SEXP Reg, SEXP Tgiven, SEXP X, SEXP Data,
			 SEXP Invcov, SEXP Notna,
			 SEXP Nx, SEXP Ngiven, SEXP Dim, SEXP Rep, SEXP Krig, 
			 SEXP Sigma2, SEXP Polydeg);

  void poly_basis_extern(int *Dim, int *Deg, int *powmatrix);
  
  // fractal dimension, 
  void boxcounting(double *z, int *lx, int * repet, double *factor,
		   int *eps, int *leps, double *sum);
  void detrendedfluc(double *dat, int *lx, int *repet, int* boxes,
		     int *ldfa, double *dfavar, double *varmethvar);
  void periodogram(double *dat, int *len, int *repet, int *fftm, 
		   int *part, int *shift, double *lambda);
  void minmax(double *dat, int *dim, int *ldim, int *boxes, int *lb, 
	      double *count);
  

  //  void Cov(int *reg, double *result);
  //  void CovMatrix(int *reg, double *result);
  // void Variogram(int *reg, double *value);
  // void Pseudovariogram(int* reg, double *result);
  SEXP CovMatrixIntern(SEXP reg, SEXP x, SEXP dist, SEXP grid,
		       SEXP lx, SEXP result);

  SEXP VariogramIntern(SEXP reg, SEXP x, SEXP lx, SEXP result);
  SEXP CovMatrixLoc(SEXP reg, SEXP x, SEXP dist, SEXP xdim, SEXP lx,
		    SEXP result);
  SEXP CovLoc(SEXP reg, SEXP x, SEXP y, SEXP xdim, SEXP lx, SEXP result);
	      

  SEXP CovMatrixSelectedLoc(SEXP reg, SEXP x, SEXP dist, SEXP xdim, SEXP lx,
			    SEXP selected, SEXP nsel, SEXP result);

  SEXP CovMatrixSelected(SEXP reg, SEXP selected, SEXP nsel, SEXP result);

  SEXP Delete_y(SEXP reg);
  void DeleteKey(int *reg);
  

  SEXP GetNAPositions(SEXP model, SEXP tsdim, SEXP xdim, SEXP integerNA, 
		      SEXP Print);
  
  SEXP SetAndGetModelInfo(SEXP model_reg, SEXP model, SEXP spatialdim,
			  SEXP distances, 
			  SEXP ygiven, // TRUE is the standard variant
			  SEXP Time, SEXP xdim, SEXP shortlen,
			  SEXP allowforintegerNA, SEXP excludetrend);
  
  
  void PutValuesAtNA(int *reg, double *values);
  void setListElements(int *reg, int *i, int *k, int *len_k);
  
  void expliciteDollarMLE(int * modelnr, double *values);

  SEXP Take2ndAtNaOf1st(SEXP model_reg, SEXP model, SEXP model_bound,
			SEXP spatialdim, SEXP Time, SEXP xdim, 
			SEXP nbounds, SEXP skipchecks);
  void UserGetNatScaling(double *natscale);
  
  void  GOUE(double * aniso, int *dim, double *grid_ext);
  
  
  //SEXP GetChar(SEXP N, SEXP Choice, SEXP Shorter, SEXP Beep, SEXP Show);
  
  void defineCovariancematrix(int *nr, int * idx, int *ncol, double *m);
  void deleteCovariancematrix();
  
  // SEXP IsStatAndIsoUser(SEXP removeGatter);
  
  
  
  
  void distInt(int *x, int*N, int *Genes, double *res);
  // void MLEMakeExplicite(double *dist, int *Lx, int *idx, int *totald);


  void GetMaxDims(int *maxints);

  void MLEanymixed(int *anymixed);
  void GetModelRegister(char **name, int* nr);
  
  void MultiDimRange(int *model_nr, double *natscale);
  void countelements(int *idx, int *N, int *boxes);
  void countneighbours(int *Xdim, int *parts, int *Squarelength, int*cumgridlen,
		       int *boxes, int * neighbours, int *OK);
  SEXP getelements(SEXP Idx, SEXP Xdim, SEXP N, SEXP Cumgridlen, SEXP Boxes);
  SEXP getneighbours(SEXP Xdim, SEXP Parts, SEXP Squarelength, SEXP Cumgridlen,
		     SEXP Neighbours);
  void duenen(int *Z, int *LL, int *BB, int *HH,
	      int *states, int *change, double *lambda, int *l_lambda,
	      int *maxtime, int *continuing, int *periodic);
  
  void Xduenen(int *ZZ, double *lambdas, double *lambdaw, double *lambdawm,
	       double *lambdawp, double *lambdad, double *lambdag, 
	       double *a, double *maxtime);
  
  
  void GetLH(int *LL, int *BB, int *HH);
  
  void intEV(int* x, int* z, int* Len, int *k, int *sumsq, int *n,
	     int *pos);
  
  SEXP GetCathegoryNames();
  void isAuthor(int *is);
  SEXP allintparam();
 
  void sleepMicro(int *micro);
  void sleepMilli(int *milli);
  void hostname(char **h, int *i);
  void pid(int *i);



}

#endif /* RF_simu_PUBLIC_H*/










