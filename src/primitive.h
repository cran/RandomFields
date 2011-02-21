#ifndef Primitives_H
#define Primitives_H 1

#define UNIT_EPSILON 1E-13

double ScaleOne(cov_model *cov);

void Bessel(double *x, cov_model *cov, double *v);
int initspectralBessel(cov_model *cov);
void spectralBessel(cov_model *cov, spectral_storage *s, double *e); 
int checkBessel(cov_model *cov);
void rangeBessel(cov_model *cov, range_arraytype *ra);


/* Cauchy models */
void Cauchy(double *x, cov_model *cov, double *v);
double ScaleCauchy(cov_model *cov);
void TBM2Cauchy(double *x, cov_model *cov, double *v);
void DCauchy(double *x, cov_model *cov, double *v);
void DDCauchy(double *x, cov_model *cov, double *v);
void invCauchySq(double *x, cov_model *cov, double *v);
int checkCauchy(cov_model *cov);
void rangeCauchy(cov_model *cov, range_arraytype *ra);
void coinitCauchy(cov_model *cov, localinfotype *li);
double DrawLogMixCauchy(cov_model *cov, mpp_storage *s); 
double LogMixWeightCauchy(double *x, cov_model *cov, mpp_storage *s);


/* another Cauchy model */
void Cauchytbm(double *x, cov_model *cov, double *v);
void DCauchytbm(double *x, cov_model *cov, double *v);
void rangeCauchytbm(cov_model *cov, range_arraytype *ra);


/* circular model */
void circular(double *x, cov_model *cov, double *v);
double Scalecircular(cov_model *cov);
void Dcircular(double *x, cov_model *cov, double *v);
void rangecircular(cov_model *cov, range_arraytype* ra);
void mppinit_circular(mpp_storage *s, cov_model *cov);
void coin_circular(mpp_storage *s, cov_model *cov); 
// res_type mppget_spherical(double *y, cov_model *cov, mpp_storage *s);


/* cone */
int checkcone(cov_model *cov);
void rangecone(cov_model *cov, range_arraytype* ra);
void mppinit_cone(mpp_storage *s, cov_model *cov);
void coin_cone(mpp_storage *s, cov_model *cov);
res_type mppget_cone(double *y, cov_model *cov, mpp_storage *s);

void constant(double *x, cov_model *cov, double *v);
void constant_nonstat(double *x, double *y, cov_model *cov, double *v);
void rangeconstant(cov_model *cov, range_arraytype* ra);
int checkconstant(cov_model *cov) ;


/* coxgauss, cmp with nsst1 !! */
// see NonIsoCovFct.cc

/* cubic */
void cubic(double *x, cov_model *cov, double *v);
void Dcubic(double *x, cov_model *cov, double *v); 
double Scalecubic(cov_model *cov);
void rangecubic(cov_model *cov, range_arraytype* ra);


/* cutoff */
// see Hypermodel.cc

/* dagum */
void dagum(double *x, cov_model *cov, double *v);
double Scaledagum(cov_model *cov); 
void Ddagum(double *x, cov_model *cov, double *v);
void rangedagum(cov_model *cov, range_arraytype* ra);
 

/*  damped cosine -- derivative of exponential:*/
void dampedcosine(double *x, cov_model *cov, double *v);
double Scaledampedcosine(cov_model *cov); 
void Ddampedcosine(double *x, cov_model *cov, double *v);
void rangedampedcosine(cov_model *cov, range_arraytype* ra);


/* De Wijsian */
void dewijsian(double *x, cov_model *cov, double *v);
void rangedewijsian(cov_model *cov, range_arraytype* ra);
/* De Wijsian */
void DeWijsian(double *x, cov_model *cov, double *v);
void rangeDeWijsian(cov_model *cov, range_arraytype* ra);
double ScaleDeWijsian(cov_model *cov); 


/* exponential model */
void exponential(double *x, cov_model *cov, double *v);
double Scaleexponential(cov_model *cov); 
void TBM2exponential(double *x, cov_model *cov, double *v);
void Dexponential(double *x, cov_model *cov, double *v);
void DDexponential(double *x, cov_model *cov, double *v);
void invexponentialSq(double *x, cov_model *cov, double *v);
int initspectralexponential(cov_model *cov);
void spectralexponential(cov_model *cov, spectral_storage *s, double *e);
void rangeexponential(cov_model *cov, range_arraytype* ra);
int checkexponential(cov_model *cov);
int hyperexponential(double radius, double *center, double *rx,
		     cov_model *cov, bool simulate, 
		     double ** Hx, double ** Hy, double ** Hr);
void coinitExp(cov_model *cov, localinfotype *li);
void ieinitExp(cov_model *cov, localinfotype *li);
double DrawLogMixExp(cov_model *cov, mpp_storage *s);
double LogMixWeightExp(double *x, cov_model *cov, mpp_storage *s);


// Brownian motion 
void fractalBrownian(double *x, cov_model *cov, double *v);
void DfractalBrownian(double *x, cov_model *cov, double *v);
void DDfractalBrownian(double *x, cov_model *cov, double *v);
void rangefractalBrownian(cov_model *cov, range_arraytype* ra);
void ieinitBrownian(cov_model *cov, localinfotype *li);
void invfractalBrownianSq(double *x, cov_model *cov, double *v);



/* FD model */
void FD(double *x, cov_model *cov, double *v);
void rangeFD(cov_model *cov, range_arraytype* ra);



/* fractgauss */
void fractGauss(double *x, cov_model *cov, double *v);
void rangefractGauss(cov_model *cov, range_arraytype* ra);


/* Gausian model */
void Gauss(double *x, cov_model *cov, double *v);
double ScaleGauss(cov_model *cov);
void DGauss(double *x, cov_model *cov, double *v);
void DDGauss(double *x, cov_model *cov, double *v);
void D3Gauss(double *x, cov_model *cov, double *v);
void D4Gauss(double *x, cov_model *cov, double *v);
int initspectralGauss(cov_model *cov);
void spectralGauss(cov_model *cov, spectral_storage *s, double *e);   
void rangeGauss(cov_model *cov, range_arraytype* ra);
double DrawLogMixGauss(cov_model *cov, mpp_storage *s);
double LogMixWeightGauss(double *x, cov_model *cov, mpp_storage *s);
void invGaussSq(double *x, cov_model *cov, double *v);
void mppinit_Gauss(mpp_storage *s, cov_model *cov);
void coin_Gauss(mpp_storage *s, cov_model *cov); 
res_type mppget_Gauss(double *y, cov_model *cov, mpp_storage *s);

void mppinit_Gausstest(mpp_storage *s, cov_model *cov);
void coin_Gausstest(mpp_storage *s, cov_model *cov); 
res_type mppget_Gausstest(double *y, cov_model *cov, mpp_storage *s);


/* generalised fractal Brownian motion */
void genBrownian(double *x, cov_model *cov, double *v);
void rangegenBrownian(cov_model *cov, range_arraytype* ra);


/* epsC */
void epsC(double *x, cov_model *cov, double *v);
void DepsC(double *x, cov_model *cov, double *v);
void DDepsC(double *x, cov_model *cov, double *v);
int checkepsC(cov_model *cov);
void rangeepsC(cov_model *cov, range_arraytype* ra);
void invepsCSq(double *x, cov_model *cov, double *v);


/* gencauchy */
void generalisedCauchy(double *x, cov_model *cov, double *v);
double ScalegeneralisedCauchy(cov_model *cov);
void DgeneralisedCauchy(double *x, cov_model *cov, double *v);
void DDgeneralisedCauchy(double *x, cov_model *cov, double *v);
int checkgeneralisedCauchy(cov_model *cov);
void rangegeneralisedCauchy(cov_model *cov, range_arraytype* ra);
void coinitgenCauchy(cov_model *cov, localinfotype *li);
void ieinitgenCauchy(cov_model *cov, localinfotype *li);
void invgeneralisedCauchySq(double *x, cov_model *cov, double *v);


/* gengneiting */
void genGneiting(double *x, cov_model *cov, double *v);
void DgenGneiting(double *x, cov_model *cov, double *v);
void rangegenGneiting(cov_model *cov, range_arraytype* ra);


/* Gneiting's functions -- alternative to Gaussian */
// #define Sqrt2TenD47 0.30089650263257344820 /* approcx 0.3 ?? */
void Gneiting(double *x, cov_model *cov, double *v); 
double ScaleGneiting(cov_model *cov);
void DGneiting(double *x, cov_model *cov, double *v); 
void rangeGneiting(cov_model *cov, range_arraytype* ra);


/* hyperbolic */
void hyperbolic(double *x, cov_model *cov, double *v); 
void Dhyperbolic(double *x, cov_model *cov, double *v); 
int checkhyperbolic(cov_model *cov);
void rangehyperbolic(cov_model *cov, range_arraytype* ra);
void invhyperbolicSq(double *x, cov_model *cov, double *v); 



/* iaco cesare model */
void IacoCesare(double *x, cov_model *cov, double *v);
//double ScaleIacoCesare(cov_model *cov); return 1.0; } 
void rangeIacoCesare(cov_model *cov, range_arraytype* ra);
//int checkIacoCesare(cov_model *cov); 


/* local-global distinguisher */
void lgd1(double *x, cov_model *cov, double *v);
double Scalelgd1(cov_model *cov);
void Dlgd1(double *x, cov_model *cov, double *v);
void rangelgd1(cov_model *cov, range_arraytype* ra);


/* mastein */
// see Hypermodel.cc


/* Whittle-Matern or Whittle or Besset ---- rescaled form of Whittle-Matern,
    see also there */ 
double Scalematern(cov_model *cov);
void Matern(double *x, cov_model *cov, double *v);
void DMatern(double *x, cov_model *cov, double *v);
void DDMatern(double *x, cov_model *cov, double *v);
void D3Matern(double *x, cov_model *cov, double *v);
void D4Matern(double *x, cov_model *cov, double *v);
int checkMatern(cov_model *cov);
void rangeMatern(cov_model *cov, range_arraytype* ra);
void ieinitMatern(cov_model *cov, localinfotype *li);
void coinitMatern(cov_model *cov, localinfotype *li);
double densityMatern(double *x, cov_model *cov);
int initspectralMatern(cov_model *cov);
void spectralMatern(cov_model *cov, spectral_storage *s, double *e); 
void invMaternSq(double *x, cov_model *cov, double *v);


/* nugget effect model */
void nugget(double *x, cov_model *cov, double *v);
double Scalenugget(cov_model *cov); 
void rangenugget(cov_model *cov, range_arraytype* ra);  
int checknugget(cov_model *cov);


/* Paciore und Stein */
void kappa_paciorek(int i, cov_model *cov, int *nr, int *nc);
void paciorek(double *x,  double *y, cov_model *cov, double *v);
int checkpaciorek(cov_model *cov);
void rangepaciorek(cov_model *cov, range_arraytype* ra);
void mppinit_paciorek(mpp_storage *s, cov_model *cov);
void coin_paciorek(mpp_storage *s, cov_model *cov);
res_type mppget_paciorek(double *x, cov_model *cov, mpp_storage *s);

/* penta */
void penta(double *x, cov_model *cov, double *v);
void Dpenta(double *x, cov_model *cov, double *v); 
double Scalepenta(cov_model *cov);
void rangepenta(cov_model *cov, range_arraytype* ra);


/* power model */ 
void power(double *x, cov_model *cov, double *v);
double Scalepower(cov_model *cov); 
void TBM2power(double *x, cov_model *cov, double *v);
void Dpower(double *x, cov_model *cov, double *v);
int checkpower(cov_model *cov);
void rangepower(cov_model *cov, range_arraytype* ra);


/* qexponential -- derivative of exponential */
void qexponential(double *x, cov_model *cov, double *v);
double Scaleqexponential(cov_model *cov);
void Dqexponential(double *x, cov_model *cov, double *v);
void rangeqexponential(cov_model *cov, range_arraytype* ra);

void kappa_rational(int i, cov_model *cov, int *nr, int *nc);
void rational(double *x, cov_model *cov, double *v);
int checkrational(cov_model *cov);
void rangerational(cov_model *cov, range_arraytype* ra);

/* spherical model */ 
void spherical(double *x, cov_model *cov, double *v);
double Scalespherical(cov_model *cov);
void TBM2spherical(double *x, cov_model *cov, double *v);
void Dspherical(double *x, cov_model *cov, double *v);
void rangespherical(cov_model *cov, range_arraytype* ra);
void mppinit_spherical(mpp_storage *s, cov_model *cov);
void coin_spherical(mpp_storage *s, cov_model *cov); 
res_type mppget_spherical(double *y, cov_model *cov, mpp_storage *s);


/* stable model */
void stable(double *x, cov_model *cov, double *v);
double Scalestable(cov_model *cov);
void Dstable(double *x, cov_model *cov, double *v);
void DDstable(double *x, cov_model *cov, double *v);
int checkstable(cov_model *cov);  
void rangestable(cov_model *cov, range_arraytype* ra);
void coinitstable(cov_model *cov, localinfotype *li);
void ieinitstable(cov_model *cov, localinfotype *li);
void invstableSq(double *x, cov_model *cov, double *v);


/* SPACEISOTROPIC stable model for testing purposes only */
void stableX(double *x, cov_model *cov, double *v);
void DstableX(double *x, cov_model *cov, double *v);

/* Stein */
// see Hypermodel.cc

/* stein space-time model */
void kappaSteinST1(int i, cov_model *cov, int *nr, int *nc);
void SteinST1(double *x, cov_model *cov, double *v);
int initspectralSteinST1(cov_model *cov);
void spectralSteinST1(cov_model *cov, spectral_storage *s, double *e);
void rangeSteinST1(cov_model *cov, range_arraytype* ra);
int checkSteinST1(cov_model *cov);  

/* wave */
void wave(double *x, cov_model *cov, double *v);
double Scalewave(cov_model *cov);
int initspectralwave(cov_model *cov);
void spectralwave(cov_model *cov, spectral_storage *s, double *e); 
void rangewave(cov_model *cov, range_arraytype* ra);

extern double Besselupperbound[Nothing + 1];

/* Whittle-Matern or Whittle or Besset */ 
void Whittle(double *x, cov_model *cov, double *v);
double ScaleWhittle(cov_model *cov);
void TBM2Whittle(double *x, cov_model *cov, double *v);
void DWhittle(double *x, cov_model *cov, double *v);
void DDWhittle(double *x, cov_model *cov, double *v);
void D3Whittle(double *x, cov_model *cov, double *v);
void D4Whittle(double *x, cov_model *cov, double *v);
double densityWhittle(double *x, cov_model *cov);
int initspectralWhittle(cov_model *cov);
void spectralWhittle(cov_model *cov, spectral_storage *s, double *e); 
int checkWhittle(cov_model *cov); 
void rangeWhittle(cov_model *cov, range_arraytype* ra);
void ieinitWhittle(cov_model *cov, localinfotype *li);
void coinitWhittle(cov_model *cov, localinfotype *li);
double DrawLogMixWM(cov_model *cov, mpp_storage *s);
double LogMixWeightW(double *x, cov_model *cov, mpp_storage *s);
res_type mppget_Whittle(double *y, cov_model *cov, mpp_storage *s);
void mppinit_Whittle(mpp_storage *s, cov_model *cov);
void coin_Whittle(mpp_storage *s, cov_model *cov);
void invWhittleSq(double *x, cov_model *cov, double *v);


void kappa_biWM(int i, cov_model *cov, int *nr, int *nc);
//void biWM(double *x, cov_model *cov, double *v);
//int checkbiWM(cov_model *cov);
//void rangebiWM(cov_model *cov, range_arraytype* ra);


void kappa_parsbiWM(int i, cov_model *cov, int *nr, int *nc);
void parsbiWM2(double *x, cov_model *cov, double *v);
int checkparsbiWM2(cov_model *cov);
void rangeparsbiWM2(cov_model *cov, range_arraytype* ra);


void biWM2(double *x, cov_model *cov, double *v);
int checkbiWM2(cov_model *cov);
void rangebiWM2(cov_model *cov, range_arraytype* ra);

/*
void userMatrix(double *x, cov_model *cov, double *v);
int checkuserMatrix(cov_model *cov);
void rangeuserMatrix(cov_model *cov, range_arraytype* ra);
void kappa_userMatrix(int i, cov_model *cov, int *nr, int *nc);
*/

void Nonestat(double *x, cov_model *cov, double *v);
void Nonenonstat(double *x, double *y, cov_model *cov, double *v);
void rangeNone(cov_model *cov, range_arraytype* ra);


#endif /* Primitives_H*/

