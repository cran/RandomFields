#ifndef RFsimu_public_H
#define RFsimu_public_H 1
#include "GSLvsR.h"

//#include <f2c.h> /* otherwise false/true not define or conflicting... */

#ifdef RF_FLOAT
typedef float Real;
#else
typedef double Real;
#endif 


#ifdef CPP2C 

EXTERN void testfct();

EXTERN void RNGtest();

EXTERN void GetParameterIndices(int *mean, int *variance, int *nugget, 
				int *scale, int *kappa, 
				int *lastkappa, int *invscale, int *sill);
EXTERN void GetrfParameters(int *covmaxchar, int *methodmaxchar,
			    int *covnr, int *methodnr);
EXTERN void GetModelNr(char **name, int *nr) ; /* transforms string "covfct" 
						into covnr of CovList */
EXTERN void GetModelName( int *nr,char **name) ; /* transforms covnr into string */
EXTERN void GetNrParameters(int *covnr,int* kappas);

EXTERN void PrintModelList();
EXTERN void PrintModelListDetailed();
// GetGridSize(*x,...) x=(start,end,step) => lengthx~=(end-start)/step+1
// (( only if grid==true !))
EXTERN void PrintMethods();
EXTERN void GetMethodName(int *nr, char **name);
EXTERN void GetMethodNr(char **name, int *nr);
EXTERN void GetGridSize(Real*,Real *,Real *,int *,int *,int *,int *); 
EXTERN void DeleteKey(int *key); // delete intermediate result for specified key
EXTERN void DeleteAllKeys();
// Covariance returns covariance values for covariance function covnr, 
// Parametervector p, and  vector of lx values of distances x; 
EXTERN void GetNaturalScaling(int *covnr, Real *p, int *actparam, 
			      int *naturalscaling, Real *natscale, int *error);
EXTERN void Covariance(Real *x,int *lx,int *covnr,Real *p, int *np, int *dim, 
		       Real *result); 
EXTERN void UncheckedCovFct(Real *x, int *n, int *covnr, Real *p, int *np, 
			    Real *result);
EXTERN void CovarianceNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
			    Real *result, int *naturalscaling); 
EXTERN void CovarianceMatrix(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
			     Real *result); 
EXTERN void CovarianceMatrixNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
			     Real *result, int *naturalscaling); 
EXTERN void Variogram(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
		      Real *result);  
EXTERN void VariogramNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
		      Real *result, int *naturalscaling); 
//EXTERN void VariogramMatrix(Real *x,int *lx,int *covnr,Real *p, Real *result); 
                                 

// SetParam****: 
//         action: 0:set; 1:get; 2:print
EXTERN void SetParamSpectral(int *action,int *nLines, int *grid);
EXTERN void SetParamDirectGauss(int *action,int *method,int *checkprecision,
				Real *requiredprecision, int *maxvariables);

// simulation procedure
EXTERN void SetParamCircEmbed( int *action, int *force, Real *tolRe, 
			       Real *tolIm,int *trials, int *mmin);
// either linesimufactor or linesimustep is used; linesimufactor is used 
// if linesimustep <= 0.0
//
// when calling init_tbm... linesimu* are stored and the stored ones are used 
// later on in the simulations !
// in contrast, always the current nLines are used!
// For further details, also on the other parameter, see library(rf),
// ?RFparameters, in R, http://www.R-project.org/
EXTERN void SetParamTBM2(int *action,int *nLines, Real *linesimufactor, 
			 Real *linesimustep );
// for simulation of 2dim RF by 3dim TBM:
EXTERN void SetParamTBM3D2(int *action,int *nLines, Real *linesimufactor, 
			   Real *linesimustep); 
EXTERN void SetParamTBM3D3(int *action,int *nLines, Real *linesimufactor, 
			   Real *linesimustep );
// for simulation of line by circular embedding
EXTERN void SetParamTBMCE(int *action,int *force,Real *tolRe,Real *tolIm,
			  int *trials, int *mtot ); 
EXTERN void SetParamTBM(int *action,int *tbm_method);

// for  MPP
EXTERN void SetMPP(int *action, Real *approxzero, 
			  Real *realisations, Real *radius);
EXTERN void printMPPinfo(int *keyNr);
EXTERN void GetMPPInfo(int *keyNr, Real *integral, Real *integralsq,
		  Real *integralpos, Real *maxheight);

// extremes
EXTERN void SetExtremes(int *action, Real *standardGausMax);

// simulation methods, Parameter setting

EXTERN void GetKeyInfo(int *keyNr,int *total, int *lengths, int *dim, 
			      int *grid, int *distr);
EXTERN void InitSimulateRF(Real *x, Real *y, Real *z, int *dim,int *lx, 
			   int *grid, /* lx=length(x) */
			   int *covnr, Real *ParamList, int *nParam, 
			   /* nParam=length(ParamList) */
			   int *method, 
			   int *distr, /* not used yet */
			   int *keyNr,
			   int *error
			   /* error>0 if error occured;
                                   =0 if fine
				   <0 message, e.g StoredInitUsed
			   */
			   );
EXTERN void DoSimulateRF(int *keyNr, Real *res, int *error); 
/* SimulateRF==InitSimulate && DoSimulateRF;
   SimulateRF is for convenience; it is slightly slower, if STORING==true and 
   several simulations are made for the same set of parameters, since it will
   always check whether the parameters have changed ((the fast way is to
   use InitSimulate only once and then call DoSimulateRF several times)) */
EXTERN void SimulateRF(Real *x, Real *y, Real *z, int *dim, int *lx, int *grid,
		       int *covnr, Real *ParamList, int *nParam,      
		       int *method, int *distr, int *keyNr,
		       Real *res, int *error);

// storing : cuurent value
// printlevel : current value
// naturalscaling : fixed value
EXTERN void SetParam(int *action,int *storing, int *printlevel,
		     int *naturalscaling, char **pch); 
/* GENERAL_STORING: 
   true: intermediate results are stored: might be rather memory consuming,
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
	 In case, only MAXKEYS different parameter sets are used, but after 
	 each step the parameter set is changed, use different keynr for each 
	 parametere set, and STORING==true, to get fast simulation 
  GENERAL_PRINTLEVEL:
	 0: no messages; 
	 1: error messages; 
         2: hints when algorithm is at a forcation; 
	 >3: more and more debugging information
  GENERAL_NATURALSCALING: (boolean)
         is called PracticalRange in R  (see Chiles&Delfiner?!) 
*/
EXTERN void GetRandomSize(int *l); // is ==-1 if not RF_GSL
#ifdef RF_GSL
EXTERN void StoreSeed();
EXTERN void RestoreSeed();
EXTERN void XGetSeed(int *seed,int *l);
EXTERN void SetSeed(int *seed, int *l, int *error);
EXTERN void GetSeed(int *seed, int *l, int *error);
#endif

EXTERN void empiricalvariogram(Real *x,Real *y,Real *z,int *dim,int *lx,
			       Real *values,  int *repet,
			       int *grid, Real *bin, int *nbin,int *charact, 
			       Real *res);
// obsolete soon?
EXTERN void binnedvariogram(Real *x, Real *y, Real *z, int *lx, Real *step,
		            Real *res, int *n, int *lr);


// obsolete !
EXTERN void empiricalXvariogram(Real *x,Real *y,Real *z,int *dim,int *lx,
				Real *values,  int *repet,
				int *grid, Real *bin, int *nbin,int *charact, 
				Real *res); // obsolete !

#else
extern void GetParameterIndices(int *mean, int *variance, int *nugget, 
				int *scale, int *kappa, 
				int *lastkappa, int *invscale, int *sill);
extern void GetrfParameters(int *covmaxchar, int *methodmaxchar,
			    int *covnr, int *methodnr);
extern void GetModelNr(char **name, int *nr) ; 
extern void GetModelName( int *nr,char **name) ; 
extern void GetNrParameters(int *covnr,int* kappas);
extern void PrintModelList();
extern void PrintModelListDetailed();
extern void PrintMethods();
extern void GetMethodName(int *nr, char **name);
extern void GetMethodNr(char **name, int *nr);
extern void GetGridSize(Real*,Real *,Real *,int *,int *,int *,int *); 
extern void DeleteKey(int *key); 
extern void DeleteAllKeys();
extern void GetNaturalScaling(int *covnr, Real *p, int *actparam, 
			      int *naturalscaling,Real *natscale, int *error);
extern void Covariance(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
		       Real *result); 
extern void CovarianceNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
			    Real *result, int *naturalscaling); 
extern void CovarianceMatrix(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
			     Real *result, int *naturalscaling); 
extern void CovarianceMatrixNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,
			     int *dim, Real *result, int *naturalscaling); 
extern void Variogram(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
		      Real *result); 
extern void VariogramNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,int *dim,
		      Real *result, int *naturalscaling); 
//extern void VariogramMatrix(Real *x,int *lx,int *covnr,Real *p, Real *result); 
extern void SetParamSpectral(int *action,int *nLines, int *grid);
extern void SetParamDirectGauss(int *action,int *method,int *checkprecision,
				Real *requiredprecision, int *maxvariables);
extern void SetParamCircEmbed(int *action, int *force, Real *tolRe, Real *tolIm,
			      int *trials);
extern void SetParamTBM2(int *action,int *nLines, Real *linesimufactor, 
			 Real *linesimustep);
extern void SetParamTBM3D2(int *action,int *nLines, Real *linesimufactor, 
			   Real *linesimustep); 
extern void SetParamTBM3D3(int *action,int *nLines, Real *linesimufactor,
			   Real *linesimustep);
extern void SetParamTBMCE(int *action,int *force,Real *tolRe,Real *tolTm,
			  int *trials);  
extern void SetParamTBM(int *action,int *tbm_method);
extern void SetMPP(int *action, Real *approxzero, 
			       Real *realisations, Real *radius);
extern void GetMPPInfo(int *keyNr, Real *integral, Real *integralsq,
		  Real *integralpos, Real *maxheight);
extern void SetExtremes(int *action, Real *standardGausMax);
extern void InitSimulateRF(Real *x, Real *y, Real *z, int *dim, int *lx, 
			   int *grid, int *covnr, Real *ParamList, int *nParam,
			   int *method, int* distr, int *keyNr, int *error);
extern void GetKeyInfo(int *keyNr,int *total, int *lengths, int *dim, 
			      int *grid, int *distr);
extern void DoSimulateRF(int *keyNr, Real *res, int *error); 
extern void SimulateRF(Real *x, Real *y, Real *z,  int *dim,int *lx, int *grid,
		       int *covnr, Real *ParamList, int *nParam,      
		       int *method, int *distr, int *keyNr,
		       Real *res, int *error);
extern void SetParam(int *action,int *storing, int *printlevel,
		     int *naturalscaling); 
extern void GetRandomSize(int *l);
#ifdef RF_GSL
extern void StoreSeed();
extern void RestoreSeed();
extern void XGetSeed(int *seed,int *l);
extern void SetSeed(int *seed, int *l, int *error);
extern void GetSeed(int *seed, int *l, int *error);
#endif
extern void empiricalvariogram(Real *x,Real *y,Real *z, int *dim,int *lx,
			       Real *values,  int *repet,int *grid, Real *bin, 
			       int *nbin,int *charact,Real *res);

// obsolete soon?
extern void binnedvariogram(Real *x, Real *y, Real *z, int *lx, Real *step,
		            Real *res, int *n, int *lr);

//extern void ResetRandom(int *fixed);
#endif

#endif /* RF_simu_PUBLIC_H*/



