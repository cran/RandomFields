
#ifndef RFsimu_public_H
#define RFsimu_public_H 1
#include "GSLvsR.h"

//#include <f2c.h> /* otherwise false/true not define or conflicting... */

EXTERN void getCov(int *n, double *COV);

EXTERN void test();

EXTERN void niceFFTnumber(int *n);

EXTERN void RNGtest();

EXTERN void GetParameterIndices(int *variance, int *scale, int *kappa, 
				int *lastkappa, int *invscale, int aniso);
EXTERN void GetrfParameters(int *covmaxchar, int *methodmaxchar,
			    int *distrmaxchar,
			    int *covnr, int *methodnr, int *distrnr,
			    int *maxdim, int *maxmodels);
EXTERN void GetModelNr(char **name, int *n, int *nr) ; /* transforms string "covfct" 
						into covnr of CovList */
EXTERN void GetModelName(int *nr,char **name) ; /* transforms covnr into string */
EXTERN void GetRange(int *nr, int *dim, int *index, Real *range, int *lrange);
EXTERN void GetNrParameters(int *covnr, int *n, int* kappas);
EXTERN void GetDistrName(int *nr,char **name);
EXTERN void GetDistrNr(char **name, int *n, int *nr);

EXTERN void PrintModelList();
EXTERN void GetModelList(int* idx);
EXTERN void PrintModelListDetailed();
// GetGridSize(*x,...) x=(start,end,step) => lengthx~=(end-start)/step+1
// (( only if grid==true !))
EXTERN void PrintMethods();
EXTERN void GetMethodName(int *nr, char **name);
EXTERN void GetMethodNr(char **name, int *n, int *nr);
EXTERN void GetGridSize(Real*,Real *,Real *,int *,int *,int *,int *); 
EXTERN void DeleteKey(int *key); // delete intermediate result for specified key
EXTERN void DeleteAllKeys();
EXTERN void printKEY(int *keyNr);
/* Natural Scaling -- in contrast to the RFparameters(PracticalRange)
   GENERAL_NATURALSCALING may take the following values (furthermore,
   0: no rescaling, +10 numerical determination allowed
*/
#define NATSCALE_EXACT 1
#define NATSCALE_APPROX 2
#define NATSCALE_MLE 3 /* check mleRF when changing !! */
EXTERN void GetNaturalScaling(int *covnr, Real *q,  /* KAPPAS only  !! */
			      int *naturalscaling, Real *natscale, int *error);
EXTERN void Covariance(Real *x,int *lx,int *covnr,Real *p, int *np, 
		       int *logicaldim, int *xdim, 
		       int *ncov,  int *anisotropy, int *op, Real *result); 
EXTERN void CovarianceNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,  
			    int *logicaldim, int *xdim, 
			    int *ncov,  int *anisotropy, int *op,
			    Real *result, int *naturalscaling); 
EXTERN void InitUncheckedCovFct(int *covnr, Real *p, int *np, 
				int *logicaldim, int *dim, 
				int *ncov, int *anisotropy, int *op,
			        int *naturalscaling,
				int *error); 
EXTERN void UncheckedCovFct(Real *x, int *n, Real *result);

EXTERN void CovarianceMatrix(Real *x,int *lx,int *covnr,Real *p,int *np, 
			     int *logicaldim, int *xdim, 
			     int *ncov, int *anisotropy, int *op,Real *result); 
EXTERN void CovarianceMatrixNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,
				  int *logicaldim, int *xdim, 
				  int *ncov, int *anisotropy, int *op,
				  Real *result, int *naturalscaling); 
EXTERN void Variogram(Real *x,int *lx,int *covnr,Real *p,int *np,
		      int *logicaldim, int *xdim, 
		      int *ncov, int *anisotropy, int *op,Real *result);  
EXTERN void VariogramNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,
			   int *logicaldim, int *xdim, 
			   int *ncov, int *anisotropy, int *op,
			   Real *result, int *naturalscaling); 
//EXTERN void VariogramMatrix(Real *x,int *lx,int *covnr,Real *p, Real *result); 
                                 

// SetParam****: 
//         action: 0:set; 1:get; 2:print
EXTERN void SetParamSpectral(int *action,int *nLines, int *grid);
EXTERN void SetParamDirectGauss(int *action,int *method,int *checkprecision,
				Real *requiredprecision, int *maxvariables);

// simulation procedure
EXTERN void SetParamCircEmbed(int *action, int *force, Real *tolRe, Real *tolIm,
			      int *trials, 
			      int *mmin, /* vector of length MAXDIM */
			      int *userfft, int *strategy);

// either linesimufactor or linesimustep is used; linesimufactor is used 
// if linesimustep <= 0.0
//
// when calling init_tbm... linesimu* are stored and the stored ones are used 
// later on in the simulations !
// in contrast, always the current nLines are used!
// For further details, also on the other parameter, see library(rf),
// ?RFparameters, in R, http://www.R-project.org/
EXTERN void SetParamTBM2(int *action, int *nLines, Real *linesimufactor, 
			 Real *linesimustep, int *every);
// for simulation of 2dim RF by 3dim TBM:
EXTERN void SetParamTBM3D2(int *action, int *nLines, Real *linesimufactor, 
			   Real *linesimustep, int *every); 
EXTERN void SetParamTBM3D3(int *action, int *nLines, Real *linesimufactor, 
			   Real *linesimustep, int *every);
// for simulation of line by circular embedding
EXTERN void SetParamTBMCE(int *action,int *force, Real *tolRe,Real *tolIm,
			  int *trials, 
			  int *mmin, /* vector of length MAXDIM */ 
			  int *userfft, int *strategy); 
EXTERN void SetParamTBM(int *action,int *tbm_method);
EXTERN void pokeTBM(int *out, int *in, int *err);

// for  MPP
EXTERN void SetMPP(int *action, Real *approxzero, 
			  Real *realisations, Real *radius);
// extremes
EXTERN void SetExtremes(int *action, Real *standardGausMax);

// simulation methods, Parameter setting
EXTERN void StoreTrend(int *keyNr, int *modus, char **trend, 
		       int *lt, Real *lambda, int *ll, int *error);
EXTERN void GetTrendLengths(int *keyNr, int *modus, int *lt, int *ll, int *lx);
EXTERN void GetTrend(int *keyNr, char **trend, Real *lambda, Real *x,int *error);
EXTERN void GetCoordinates(int *keyNr, Real *x, int *error);

EXTERN void InitMaxStableRF(Real *x, Real *T, int *dim, int *lx, int *grid, 
		     int *Time,
		     int *covnr, Real *ParamList, int *nParam, 
		     Real *mean,
		     int *ncov, int *anisotropy, int *op,
		     int *method, 
		     int *distr,
		     int *keyNr,
		     int *error);

EXTERN void DoMaxStableRF(int *keyNr, Real *res, int *error); 

EXTERN void GetKeyInfo(int *keyNr,int *total, int *lengths, int *dim, 
		       int *timespacedim, int *grid, int *distr, int *maxdim);


/* 
   check with InitSimulateRF in case of any changes !!
   both InitSimulateRF and DoSimulateRF
*/
#define DISTR_GAUSS 0   
#define DISTR_POISSON 1
#define DISTR_MAXSTABLE 2
EXTERN void InitSimulateRF(Real *x, Real *T, int *dim,
			   int *lx,  /* lx=length(x) */
			   int *grid, int *Time,			   
			   int *covnr, Real *ParamList, 
			   int *nParam, /* nParam=length(ParamList) */
			   Real *mean, int *ncov, int *anisotropy, int *op,
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
EXTERN void SimulateRF(Real *x, Real *T, int *dim,int *lx, int *grid, 
		       int *Time, int *covnr, Real *ParamList, int *nParam,
		       Real *mean, int *ncov, int *anisotropy, int *op,
		       int *method, int *distr, int *keyNr,
		       Real *res, int *error);
// storing : current value
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
	 >3: more and more debugg
ing information
  GENERAL_NATURALSCALING: (boolean)
         is called PracticalRange in R  (see Chiles&Delfiner?!) 
*/

// for isotropic spatial data
EXTERN void empiricalvariogram(Real *x, int *dim, int *lx,
			       Real *values, int *repet,
			       int *grid, Real *bin, int *nbin, int *charact, 
			       Real *res, Real *sd, int *n);

// for space-time and anisotropies
EXTERN void empvarioXT(Real *X, Real *T, 
		int *lx, 
		Real *values, int *repet, int *grid,
		Real *bin, int *nbin, 
		Real *phi,    // vector of a real and an integer
		Real *theta,  // vector of a real and an integer
		int *dT,   // 
		Real *sum,   // \sum (a-b)^2 / 2
		Real *sq,   // \sum (a-b)^4 / 4
		int *n);


// fractal dimension, 
EXTERN void boxcounting(double *z, int *lx, int * repet, double *factor,
			int *eps, int *leps, double *sum);
EXTERN void detrendedfluc(double *dat, int *lx, int *repet, int* boxes,
			  int *ldfa, double *dfavar, double *varmethvar);
EXTERN void periodogram(double *dat, int *len, int *repet, int *fftm, 
			int *part, int *shift, double *lambda);
EXTERN void minmax(double *dat, int *dim, int *ldim, int *boxes, int *lb, 
		   double *count);

#endif /* RF_simu_PUBLIC_H*/










