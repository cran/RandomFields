
#ifndef RFsimu_public_H
#define RFsimu_public_H 1
#include "GSLvsR.h"

//#include <f2c.h> /* otherwise false/true not define or conflicting... */

EXTERN void getCov(int *n, double *COV);

EXTERN void test();

EXTERN void niceFFTnumber(int *n);

EXTERN void RNGtest();

EXTERN void GetrfParameters(int *covmaxchar, int *methodmaxchar,
			    int *distrmaxchar,
			    int *covnr, int *methodnr, int *distrnr,
			    int *maxdim, int *maxmodels);
EXTERN void GetParamterPos(int *variance, int *kappa, int* lastkappa, 
		    int *tbm2num, int *hyperinternal, int *lasthyperinternal,
		    int *scale, int *aniso, int *hypernr, int *localdiameter,
		    int *local_r, int *cutoff_theo_r, int *hyperkappa,
		    int *total);
EXTERN void GetModelNr(char **name, int *n, int *nr) ; /* transforms string "covfct" 
						into covnr of CovList */
EXTERN void GetModelName(int *nr,char **name) ; /* transforms covnr into string */
EXTERN void GetRange(int *nr, int *dim, int *index, double *range, int *lrange);
EXTERN void GetNrParameters(int *covnr, int *n, int* dim, int* kappas);
EXTERN void GetDistrName(int *nr,char **name);
EXTERN void GetDistrNr(char **name, int *n, int *nr);
EXTERN void GetModelException(int *nr, int *n, int *scale, int*aniso);

EXTERN void PrintModelList();
EXTERN void GetModelList(int* idx);
EXTERN void PrintModelListDetailed();
// GetGridSize(*x,...) x=(start,end,step) => lengthx~=(end-start)/step+1
// (( only if grid==true !))
EXTERN void PrintMethods();
EXTERN void GetMethodName(int *nr, char **name);
EXTERN void GetMethodNr(char **name, int *n, int *nr);
EXTERN void GetGridSize(double*,double *,double *,int *,int *,int *,int *); 
EXTERN void DeleteKey(int *key); // delete intermediate result for specified key
EXTERN void DeleteAllKeys();
EXTERN void printKEY(int *keyNr);
/* Natural Scaling -- in contrast to the RFparameters(PracticalRange)
   GENERAL_NATURALSCALING may take the following values (furthermore,
   0: no rescaling, +10 numerical determination allowed
*/
#define NATSCALE_EXACT 1
#define NATSCALE_APPROX 2
#define NATSCALE_MLE 3 /* check fitvario when changing !! */
EXTERN void GetNaturalScaling(int *covnr, double *q,  /* KAPPAS only  !! */
			      int *naturalscaling, double *natscale, int *error);
EXTERN void CheckAndCompleteParameters(int *covnr, double *p, int *np, 
				       int *dim, /* timespacedim ! */
				       int *ncov, int *anisotropy, int *op,
				       double *q, int* error);
EXTERN void Covariance(double *x,int *lx,int *covnr,double *p, int *np, 
		       int *logicaldim, int *xdim, 
		       int *ncov,  int *anisotropy, int *op, double *result); 
EXTERN void CovarianceNatSc(double *x,int *lx,int *covnr,double *p,int *np,  
			    int *logicaldim, int *xdim, 
			    int *ncov,  int *anisotropy, int *op,
			    double *result, int *naturalscaling); 
EXTERN void InitUncheckedCovFct(int *covnr, double *p, int *np, 
				int *logicaldim, int *dim, 
				int *ncov, int *anisotropy, int *op,
			        int *naturalscaling,
				int *error); 
EXTERN void UncheckedCovFct(double *x, int *n, double *result);
EXTERN void UncheckedCovMatrix(double *dist, int *lx, double *result);

EXTERN void CovarianceMatrix(double *x,int *lx,int *covnr,double *p,int *np, 
			     int *logicaldim, int *xdim, 
			     int *ncov, int *anisotropy, int *op,double *result); 
EXTERN void CovarianceMatrixNatSc(double *x,int *lx,int *covnr,double *p,int *np,
				  int *logicaldim, int *xdim, 
				  int *ncov, int *anisotropy, int *op,
				  double *result, int *naturalscaling); 
EXTERN void Variogram(double *x,int *lx,int *covnr,double *p,int *np,
		      int *logicaldim, int *xdim, 
		      int *ncov, int *anisotropy, int *op,double *result);  
EXTERN void VariogramNatSc(double *x,int *lx,int *covnr,double *p,int *np,
			   int *logicaldim, int *xdim, 
			   int *ncov, int *anisotropy, int *op,
			   double *result, int *naturalscaling); 
//EXTERN void VariogramMatrix(double *x,int *lx,int *covnr,double *p, double *result); 
                                 

// SetParam****: 
//         action: 0:set; 1:get; 2:print
EXTERN void SetParamDecision( int *action, int *stationary_only, 
			      int *exactness);
EXTERN void SetParamSpectral(int *action,int *nLines, int *grid);
EXTERN void SetParamDirectGauss(int *action,int *method,int *checkprecision,
				double *requiredprecision, int *bestvariables,
				int *maxvariables);

// simulation procedure
EXTERN void SetParamCircEmbed(int *action, int *force, double *tolRe, 
			      double *tolIm,
			      int *trials, 
			      double *mmin, /* vector of length MAXDIM */
			      int *userfft, int *strategy, double*maxmem,
			      int *dependent);

// for local circulant embedding
EXTERN void SetParamLocal( int *action, int *force, double *tolRe, 
			   double *tolIm,
			   int *trials, double *mmin, int *useprimes, 
			   int *strategy,
			   double *maxmem, int *dependent);

// either linesimufactor or linesimustep is used; linesimufactor is used 
// if linesimustep <= 0.0
//
// when calling init_tbm... linesimu* are stored and the stored ones are used 
// later on in the simulations !
// in contrast, always the current nLines are used!
// For further details, also on the other parameter, see library(rf),
// ?RFparameters, in R, http://www.R-project.org/
EXTERN void SetParamTBM2(int *action, int *nLines, double *linesimufactor, 
			 double *linesimustep, int *every, int *tbm2num);
EXTERN void SetParamTBM3(int *action, int *nLines, double *linesimufactor, 
			   double *linesimustep, int *every);
// for simulation of line by circular embedding
EXTERN void SetParamTBMCE(int *action,int *force, double *tolRe,double *tolIm,
			  int *trials, 
			  double *mmin, /* vector of length MAXDIM */ 
			  int *userfft, int *strategy, double*maxmem,
			  int *dependent); 
EXTERN void SetParamTBM(int *action,int *tbm_method, double *center, 
			int *points);
EXTERN void SetParamHyperplane(int *action, int *superpos, int *maxlines, 
			       int *normalise, int *mar_distr, 
			       double *mar_param);
// for  MPP
EXTERN void SetMPP(int *action, double *approxzero, 
			  double *realisations, double *radius);
// extremes
EXTERN void SetExtremes(int *action, double *standardGausMax);
EXTERN void SetParamNugget(int *action, double *nuggettol);

// simulation methods, Parameter setting
EXTERN void StoreTrend(int *keyNr, int *modus, char **trend, 
		       int *lt, double *lambda, int *ll, int *error);
EXTERN void GetTrendLengths(int *keyNr, int *modus, int *lt, int *ll, int *lx);
EXTERN void GetTrend(int *keyNr, char **trend, double *lambda, double *x,int *error);
EXTERN void GetCoordinates(int *keyNr, double *x, int *error);

EXTERN void InitMaxStableRF(double *x, double *T, int *dim, int *lx, int *grid, 
		     int *Time,
		     int *covnr, double *ParamList, int *nParam, 
		     int *ncov, int *anisotropy, int *op,
		     int *method, 
		     int *distr,
		     int *keyNr,
		     int *error);

EXTERN void DoMaxStableRF(int *keyNr, int *n, int *pairs,
			  double *res, int *error); 

EXTERN void GetKeyInfo(int *keyNr,int *total, int *lengths, int *dim, 
		       int *timespacedim, int *grid, int *distr, int *maxdim);
EXTERN SEXP GetExtKeyInfo(SEXP keynr, SEXP Ignoreactive);
EXTERN SEXP GetExtModelInfo(SEXP keynr);

/* 
   check with InitSimulateRF in case of any changes !!
   both InitSimulateRF and DoSimulateRF
*/
#define DISTR_GAUSS 0   
#define DISTR_POISSON 1
#define DISTR_MAXSTABLE 2
EXTERN void InitSimulateRF(double *x, double *T, int *dim,
			   int *lx,  /* lx=length(x) */
			   int *grid, int *Time,			   
			   int *covnr, double *ParamList, 
			   int *nParam, /* nParam=length(ParamList) */
			   int *ncov, int *anisotropy, int *op,
			   int *method, 
			   int *distr, /* not used yet */
			   int *keyNr,
			   int *error
			   /* error>0 if error occured;
                                   =0 if fine
				   <0 message, e.g StoredInitUsed
			   */
			   );
EXTERN void DoSimulateRF(int *keyNr, int *n, int *pairs,
			 double *res, int *error); 
EXTERN void AddTrend(int *keyNr, int *n, double *res, int *error);

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
EXTERN void empiricalvariogram(double *x, int *dim, int *lx,
			       double *values, int *repet,
			       int *grid, double *bin, int *nbin, int *charact, 
			       double *res, double *sd, int *n);

// for space-time and anisotropies
EXTERN void empvarioXT(double *X, double *T, 
		int *lx, 
		double *values, int *repet, int *grid,
		double *bin, int *nbin, 
		double *phi,    // vector of a real and an integer
		double *theta,  // vector of a real and an integer
		int *dT,   // 
		double *sum,   // \sum (a-b)^2 / 2
		double *sq,   // \sum (a-b)^4 / 4
		int *n);


// kriging methods
EXTERN void simpleKriging(double *tgiven, double *x, double *invcov, 
			  int *Len_x, int *NN,
			  int *Dim, int *Rep, double *Mean, double *Res);
EXTERN void simpleKriging2(double *tgiven, double *x, double *data,
			   double *invcov,
			   int *Len_x, int *NN, int *Dim, int *Rep, 
			   double *Mean, double *Res, double *sigma2);
EXTERN void ordinaryKriging(double *tgiven, double *x, double *data, 
			    int *Len_x, int *NN,
			    int *Dim, int *Rep, double *Res);
EXTERN void ordinaryKriging2(double *tgiven, double *x, double *data, 
			     double *invcov,
			     int *Len_x, int *NN, int *Dim, int *Rep, 
			     double *Res, double *sigma2);

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










