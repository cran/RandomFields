
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

extern "C" void getCov(int *n, double *COV);

extern "C" void test();

extern "C" void niceFFTnumber(int *n);

extern "C" void RNGtest();

extern "C" void GetrfParameters(int *covmaxchar, int *methodmaxchar,
			    int *distrmaxchar,
			    int *covnr, int *methodnr, int *distrnr,
			    int *maxdim, int *maxmodels);
extern "C" void GetrfParametersI(int *covmaxchar, int *methodmaxchar,
			    int *distrmaxchar,
			    int *covnr, int *methodnr, int *distrnr,
			    int *maxdim, int *maxmodels);
extern "C" void GetModelNr(char **name, int *n); /* transforms string "covfct" 
						into covnr of CovList */
extern "C" void GetModelName(int *nr,char **name); /* transforms covnr into string */
extern "C" void GetAttr(int *hyper, int *normalmix, int *finiterange, 
		    int *internal, int *stat, int *iso, int *maxdim, int *vdim);
extern "C" void GetRange(int *nr, int *dim, int *index, double *params, int *lparam,
		      double *range, int *lrange);
extern "C" void GetNrParameters(int *covnr, int* kappas);
extern "C" void GetDistrName(int *nr,char **name);
extern "C" void GetDistrNr(char **name, int *n, int *nr);

extern "C" void PrintModelList(int *local, int *operat);
extern "C" void GetModelList(int* idx, int*internal);
extern "C" void PrintMethods();
extern "C" void GetMethodName(int *nr, char **name);
extern "C" void GetMethodNr(char **name, int *nr);
//extern "C" void GetGridSize(double*,double *,double *,int *,int *,int *,int *); 
extern "C" void DeleteKey(int *key); // delete intermediate result for specified key
extern "C" void DeleteAllKeys();
extern "C" void printKEY(int *keyNr);
/* Natural Scaling -- in contrast to the RFparameters(PracticalRange)
   GENERAL_NATURALSCALING may take the following values (furthermore,
   0: no rescaling, +10 numerical determination allowed
*/
#define NATSCALE_EXACT 1 /* or approx or mle */
#define NATSCALE_ORNUMERIC 2
//#define NATSCALE_APPROX 2
#define NATSCALE_MLE 3 /* check fitvario when changing; best not to change !!! */

//extern "C" void GetNaturalScaling(int *covnr, double *q,  /* KAPPAS only  !! */
//			      int *naturalscaling, double *natscale, int *error);



// SetParam****: 
//         action: 0:set; 1:get; 

extern "C" void SetParamAve(int *action,  
			double *spatialMAr, double *spatialMAdist
    );


extern "C" void SetParamMPP(int *action, int *locations,
			double *intens, double *plus,  
			double *relplus, 
			//	double *scale,
			double *approxzero,
			double *samplingdist,
			double *samplingr,
			double *p, double *beta
  );
// simulation procedure
extern "C" void SetParamCircEmbed(int *action, int *force, double *tolRe, 
			      double *tolIm,
			      int *trials,  
			      double *mmin, /* vector of length MAXDIM */
			      int *useprimes, int *strategy, double*maxmem,
			      int *dependent, int *method);

// for local circulant embedding
extern "C" void SetParamLocal( int *action, int *force, double *tolRe, 
			   double *tolIm,  
			   double *mmin, int *useprimes, 
			   double *maxmem, int *dependent);
extern "C" void SetParamDirect(int *action,int *method, int *bestvariables,
				int *maxvariables, double *svdtolerance);

extern "C" void SetParamMarkov(int *action,int *neighbours, double *precision,
			   int *cyclic, int *maxmem);

extern "C" void SetParamNugget(int *action, double *nuggettol, int *nuggetmeth);
extern "C" void SetParamSequential(int *action, int *maxvariables, int *back, 
			       int *start);
// either linesimufactor or linesimustep is used; linesimufactor is used 
// if linesimustep <= 0.0
//
// when calling init_tbm... linesimu* are stored and the stored ones are used 
// later on in the simulations !
// in contrast, always the current nLines are used!
// For further details, also on the other parameter, see library(rf),
// ?RFparameters, in R, http://www.R-project.org/
extern "C" void SetParamSpectral(int *action,int *nLines, int *grid, 
			     int *ergodic,
			     int *metropolis, int *nmetro, double *sigma);
// for simulation of line by circular embedding
extern "C" void SetParamTBMCE(int *action,int *force, double *tolRe,double *tolIm,
			  int *trials,
			  double *mmin, /* vector of length MAXDIM */ 
			  int *userfft, int *strategy, double*maxmem,
			  int *dependent); 
extern "C" void SetParamTBM2(int *action, int *nLines, double *linesimufactor, 
			 double *linesimustep, 
			 int *layers);
extern "C" void SetParamTBM3(int *action, int *nLines, double *linesimufactor, 
			 double *linesimustep,  
			 int *layers);
extern "C" void SetParamTBM(int *action,int *tbm_method, double *center, 
			int *points);
// for  MPP
extern "C" void SetParamHyperplane(int *action, int *superpos, int *maxlines, 
			       int *mar_distr, double *mar_param);
//extern "C" void SetParamSpecial(int *action, double *spatialMAr, 
//			    double *spatialMAdist);

// extremes
extern "C" void SetExtremes(int *action, double *standardGausMax);

// simulation methods, Parameter setting
extern "C" void StoreTrend(int *keyNr, int *modus, char **trend, 
		       int *lt, double *lambda, int *ll, int *error);
extern "C" void GetTrendLengths(int *keyNr, int *modus, int *lt, int *ll, int *lx);
extern "C" void GetTrend(int *keyNr, char **trend, double *lambda, double *x,int *error);
extern "C" void GetCoordinates(int *keyNr, double *x, int *error);

extern "C" void InitMaxStableRF(double *x, double *T, int *dim, int *lx, int *grid, 
		     int *Time,
		     int *distr,
		     int *keyNr,
		     int expected_number_simu,
		     int *err);

extern "C" void DoMaxStableRF(int *keyNr, int *n, int *pairs,
			  double *res, int *error); 

extern "C" void GetKeyInfo(int *keyNr,int *total, int *lengths, int *dim, 
		       int *timespacedim,
		       int *grid, int *distr, int *maxdim, int *vdim);
extern "C" SEXP GetRegisterInfo(SEXP keynr, SEXP Ignoreactive, SEXP max_element);
extern "C" SEXP GetExtModelInfo(SEXP keynr, SEXP level, SEXP gatter);
extern "C" SEXP GetModel(SEXP keynr, SEXP natscale);

/* 
   check with InitSimulateRF in case of any changes !!
   both InitSimulateRF and DoSimulateRF
*/
#define DISTR_GAUSS 0   
#define DISTR_POISSON 1
#define DISTR_MAXSTABLE 2
extern "C" void InitSimulateRF(double *x, double *T, 
		    int *spatialdim, /* spatial dim only ! */
		    int *lx, 
		    int *grid,
		    int *Time, 
		    int *distr, /* still unused */
		    int *keyNr, 
		    int *expected_number_simu,
		    int *error
		      /* error>0 if error occured;
                                   =0 if fine
				   <0 message, e.g StoredInitUsed
		      */
  );
extern "C" void DoSimulateRF(int *keyNr, int *n, int *pairs,
			     res_type *res, int *error);  
extern "C" void AddTrend(int *keyNr, int *n, double *res, int *error);

// storing : current value
// printlevel : current value
// naturalscaling : fixed value
extern "C" SEXP  SetParamPch(SEXP act, SEXP pch);
extern "C" void SetParam(int *action,int *storing, int *printlevel,
		     int *naturalscaling, int *skipchecks,
		     int *every); 
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
  PL:
	 0: no messages; 
	 1: error messages; 
         2: hints when algorithm is at a forcation; 
	 >3: more and more debugg
ing information
  GENERAL_NATURALSCALING: (boolean)
         is called PracticalRange in R  (see Chiles&Delfiner?!) 
*/
extern "C" void SetParamDecision(int *action, int *stationary_only,
			     int *exactness);

// for isotropic spatial data
extern "C" void empiricalvariogram(double *x, int *dim, int *lx,
			       double *values, int *repet,
			       int *grid, double *bin, int *nbin, int *charact, 
			       double *res, double *sd, int *n);

// for space-time and anisotropies
extern "C" void empvarioXT(double *X, double *T, 
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
extern "C" void simpleKriging(double *tgiven, double *x, double *invcov, 
			      int *notna, int *Len_x, int *NN,
			  int *Dim, int *Rep, double *Mean, double *Res);
extern "C" void simpleKriging2(double *tgiven, double *x, double *data,
			   double *invcov,
			   int *Len_x, int *NN, int *Dim, int *Rep, 
			   double *Mean, double *Res, double *sigma2);
extern "C" void ordinaryKriging(double *tgiven, double *x, double *data, 
			    int *Len_x, int *NN,
			    int *Dim, int *Rep, double *Res);
extern "C" void ordinaryKriging2(double *tgiven, double *x, double *data, 
			     double *invcov,
			     int *Len_x, int *NN, int *Dim, int *Rep, 
			     double *Res, double *sigma2);
extern "C" void ordinaryKriging(double *tgiven, double *x, double *data, 
			    int *Len_x, int *NN,
			    int *Dim, int *Rep, double *Res);
extern "C" void ordinaryKriging2(double *tgiven, double *x, double *data, 
			     double *invcov,
			     int *Len_x, int *NN, int *Dim, int *Rep, 
			     double *Res, double *sigma2);
extern "C" void universalKriging(double *tgiven, double *x, double *invcov,
			     int *Len_x, int *NN, int *Dim, int *Rep,
			     double *Res, double *tFmatrix, int *NFCT,
			     int *iR, int *jR, SEXP f_expr, SEXP rho);
extern "C" void universalKriging2(double *tgiven, double *x, double *data,
                              double *invcov, int *Len_x, int *NN, int *Dim,
                              int *Rep, double *Res, double *sigma2,
                              double *tFmatrix, int *NFCT,
			      int *iR, int *jR, SEXP f_expr, SEXP rho);
extern "C" void universalKriging3(double *tgiven, double *x, double *invcov,
			      int *Len_x, int *NN, int *Dim, int *Rep, 
			      double *Res, double *tPowers, int *NFCT);
extern "C" void universalKriging4(double *tgiven, double *x, double *data,
			      double *invcov, int *Len_x, int *NN, int *Dim,
			      int *Rep, double *Res, double *sigma2,
			      double *tPowers, int *NFCT);
extern "C" void intrinsicKriging(double *tgiven, double *x, double *invcov,
                             int *Len_x, int *NN, int *Dim, int *Rep,
                             double *Res, double *tPowers, int *NFCT);
extern "C" void intrinsicKriging2(double *tgiven, double *x, double *data,
                              double *invcov, int *Len_x, int *NN, int *Dim,
                              int *Rep, double *Res, double *sigma2, 
                              double *tPowers, int *NFCT,
			      SEXP inv_expr, SEXP rho);

// fractal dimension, 
extern "C" void boxcounting(double *z, int *lx, int * repet, double *factor,
			int *eps, int *leps, double *sum);
extern "C" void detrendedfluc(double *dat, int *lx, int *repet, int* boxes,
			  int *ldfa, double *dfavar, double *varmethvar);
extern "C" void periodogram(double *dat, int *len, int *repet, int *fftm, 
			int *part, int *shift, double *lambda);
extern "C" void minmax(double *dat, int *dim, int *ldim, int *boxes, int *lb, 
		   double *count);


extern "C" void Cov(double *x, double *y, int *lx, double *result);
extern "C" void Variogram(double *x,double *y, int *lx, double *result);
extern "C" void VariogramIntern(double *x,double *y, int *lx, double *result);
extern "C" void VariogramMLE(double *x, double *y, int *lx, double *result);

extern "C" void CovMatrix(double *x, int * dist, int *lx, double *result);
extern "C" void CovMatrixIntern(double *x, int * dist, int *lx, double *result);
extern "C" void CovMatrixMLE(double *x, int *dist, int *lx,  int *set,
			     int *sub, double *result);
extern "C" SEXP CheckModelUser(SEXP model, SEXP tsdim, SEXP xdim, 
			       SEXP stationary);
extern "C" SEXP CheckModelIntern(SEXP model, SEXP tsdim, SEXP xdim, 
				 SEXP stationary);
extern "C" SEXP CheckModelSimu(SEXP model, SEXP tsdim, SEXP xdim, 
			       SEXP stationary);

extern "C" SEXP GetNAPositions(SEXP model, SEXP tsdim, SEXP xdim,
			       SEXP stationary, 
			       SEXP print);
extern "C" SEXP MLEGetModelInfo(SEXP model, SEXP tsdim, SEXP xdim); // ex GetNA
extern "C" void PutValuesAtNA(double *values);
extern "C" void expliciteDollarMLE(double *values);
extern "C" SEXP Take2ndAtNaOf1st(SEXP model, SEXP model_bound, SEXP tsdim, 
				 SEXP xdim, SEXP stationary, SEXP nbounds, 
				 SEXP skipchecks);
extern "C" void UserGetNatScaling(double *natscale);

extern "C" void  GOUE(double * aniso, int *dim, double *grid_ext);


extern "C" SEXP GetChar(SEXP N, SEXP Choice, SEXP Shorter, SEXP Beep, SEXP Show);

extern "C" void defineCovariancematrix(int *nr, int * idx, int *ncol, double *m);
extern "C" void deleteCovariancematrix();

extern "C" SEXP IsStatAndIsoUser(SEXP removeGatter);
 



extern "C" void distInt(int *x, int*N, int *Genes, double *res);
// extern "C" void MLEMakeExplicite(double *dist, int *Lx, int *idx, int *totald);

extern "C" void meth0(double *poisson, int *NUMBER, double *trend,
                      int *LEN, int *KMAX, int *NULLINDEX, 
                      double *result, int *OUTPUT);
extern "C" void meth1(double *poisson, int *NUMBER, double *trend, int *TDIM,
                      int *N, int *KMAX, int *NULLINDEX,
	              int *TRENDNULLINDEX, int *hindex, int *trendhindex,
	              int *STARTINDEX, int *TRENDSTARTINDEX,
                      int *reslenR, int *proclenR, int *trendlenR, 
	              double *result, int *OUTPUT);
extern "C" void meth2(double *poisson, double *t, int *TDIM,
                      int *NUMBER, double *trend,
                      double *trendt, int *proclenR, int *trendlenR,
                      int *KMAX, int *TRENDNULLINDEX,
                      int *STARTINDEX, int *TRENDSTARTINDEX,
                      double *lower, double *upper,
                      int *reslenR, double *result, int *OUTPUT);
extern "C" void meth3(double *poisson, int *NUMBER, double *trend,
                      int *nR, int *TDIM, int *KMAX, int *TRENDNULLINDEX,
                      int *STARTINDEX, int *M, int *reslenR, double *result,
                      int *OUTPUT);
extern "C" void meth4(double *LAMBDA, double *poisson, int *NUMBER, int *TDIM,
                      double *trend, double *t, int *proclenR, int *nullindexR,
                      int *STARTINDEX, int *TRENDPROCSTARTINDEX,
                      int *trendlenR, int *intstartindexR, double *lower,
                      double *upper, int *intlenR, int *reslenR, int *KMAX,
                      double *result, int *OUTPUT);

extern "C" void GetMaxDims(int *maxints);


extern "C" void analyseForst(int* keynr, 
			     double *gausspercent, int *ngaussp,
			     int *resol, int *nres, 
			     int *edges, int *nedg,
			     double *percent, int *nperc,
			     double *result);

extern "C" void analyseForstImages(int* keynr, 
			double *gausspercent, int *ngaussp,
			int *resol, int *nres, 
			int *edges, int *nedg,
			double *percent, int *nperc,
			double *image,
			int *binary,
			int * decreased,
			int *nrdecr,
			int *ncdecr,
			int *refarea, int* areathreshold,
			double *result		
    );

extern "C" void MLEanymixed(int *anymixed);

#endif /* RF_simu_PUBLIC_H*/










