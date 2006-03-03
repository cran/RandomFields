#ifndef CovFcts_H
#define CovFcts_H 1


// TODO: TBM3bessel
double Bessel(double *x,double *p);
//double ScaleBessel(double *p,int scaling); missing
double spectralBessel(double *p );
int checkBessel(double *param, int timespacedim, SimulationType method);
void rangeBessel(int dim, int *index, double* range);
void infoBessel(double *p, int *maxdim, int *CEbadlybehaved);

double Cauchy(double *x, double *p);
double ScaleCauchy(double *p,int scaling);
double TBM2Cauchy(double *x, double *p);
double TBM3Cauchy(double *x, double *p);
double DCauchy(double *x, double *p);
double DDCauchy(double *x, double *p);
int checkCauchy(double *param, int timespacedim, SimulationType method);
void rangeCauchy(int dim, int *index, double* range);
void infoCauchy(double *p, int *maxdim, int *CEbadlybehaved);

double Cauchytbm(double *x, double *p);
//double ScaleCauchytbm(double *p,int scaling, int dim); missing
double TBM3Cauchytbm(double *x, double *p);
double DCauchytbm(double *x, double *p);
int checkCauchytbm(double *param, int timespacedim, SimulationType method);
void rangeCauchytbm(int dim, int *index, double* range);
void infoCauchytbm(double *p, int *maxdim, int *CEbadlybehaved);

double circular(double *x, double *p);
double Scalecircular(double *p,int scaling);
double Dcircular(double *x, double *p);
void rangecircular(int dim, int *index, double* range);
void infocircular(double *p, int *maxdim, int *CEbadlybehaved);

double constant(double *x, double *p);
double TBM2constant(double *x, double *p);
double TBM3constant(double *x, double *p);
double Dconstant(double *x, double *p);
void rangeconstant(int dim, int *index, double* range);
void infoconstant(double *p, int *maxdim, int *CEbadlybehaved);

int checkcone(double *param, int timespacedim, SimulationType method);
void rangecone(int dim, int *index, double* range);
void infocone(double *p, int *maxdim, int *CEbadlybehaved);

double cubic(double *x, double *p);
double Scalecubic(double *p,int scaling);
double TBM3cubic(double *x, double *p);
double Dcubic(double *x, double *p);
void rangecubic(int dim, int *index, double* range);
void infocubic(double *p, int *maxdim, int *CEbadlybehaved);

double dampedcosine(double *x, double*p);
double Scaledampedcosine(double *p,int scaling);
double TBM3dampedcosine(double *x, double *p);
double Ddampedcosine(double *x, double *p);
void rangedampedcosine(int dim, int *index, double* range);
void infodampedcosine(double *p, int *maxdim, int *CEbadlybehaved);

double exponential(double *x,double *p);
double Scaleexponential(double *p,int scaling);
double TBM2exponential(double *x, double *p);
double TBM3exponential(double *x, double *p);
double Dexponential(double *x, double *p);
double DDexponential(double *x, double *p);
double spectralexponential(double *p );
int hyperexponential(double radius, double *center, double *rx,
		     int dim, bool simulate, 
		     double** hx, double** hy, double** hr);
void rangeexponential(int dim, int *index, double* range);
void infoexponential(double *p, int *maxdim, int *CEbadlybehaved);
int checkexponential(double *param, int timespacedim, SimulationType method);

double fractalBrownian(double*x, double *p);
double DfractalBrownian(double *x, double*p);
double DDfractalBrownian(double *x, double*p);
void rangefractalBrownian(int dim, int *index, double* range);
void infofractalBrownian(double *p, int *maxdim, int *CEbadlybehaved);

double FD(double *x,double *p);
void rangeFD(int dim, int *index, double* range);
void infoFD(double *p, int *maxdim, int *CEbadlybehaved);

double fractGauss(double *x,double *p);
void rangefractGauss(int dim, int *index, double* range);
void infofractGauss(double *p, int *maxdim, int *CEbadlybehaved);

double Gauss(double *x, double*p);
double ScaleGauss(double*p,int scaling);
double TBM3Gauss(double *x, double*p);
double DGauss(double *x, double*p);
double spectralGauss(double *p );
void rangeGauss(int dim, int *index, double* range);
void infoGauss(double *p, int *maxdim, int *CEbadlybehaved);

double generalisedCauchy(double *x, double *p);
double ScalegeneralisedCauchy(double *p,int scaling);
double TBM3generalisedCauchy(double *x, double *p);
double DgeneralisedCauchy(double *x, double *p);
double DDgeneralisedCauchy(double *x, double *p);
int checkgeneralisedCauchy(double *param, int timespacedim, 
			   SimulationType method);
void rangegeneralisedCauchy(int dim, int *index, double* range);
void infogeneralisedCauchy(double *p, int *maxdim, int *CEbadlybehaved);

double genGneiting(double *x, double *p);
double TBM3genGneiting(double *x, double *p);
double DgenGneiting(double *x, double *p);
void rangegenGneiting(int dim, int *index, double* range);
void infogenGneiting(double *p, int *maxdim, int *CEbadlybehaved);

double Gneiting(double *x, double *p);
double ScaleGneiting(double *p,int scaling);
double TBM3Gneiting(double *x, double *p);
double DGneiting(double *x, double *p);
void rangeGneiting(int dim, int *index, double* range);
void infoGneiting(double *p, int *maxdim, int *CEbadlybehaved);

double hyperbolic(double *x, double*p);
// Scalehyperbolic(double *p,int scaling); missing
int checkhyperbolic(double *param, int timespacedim, SimulationType method);
double TBM3hyperbolic(double *x, double*p);
double Dhyperbolic(double *x, double*p);
void rangehyperbolic(int dim, int *index, double* range);
void infohyperbolic(double *p, int *maxdim, int *CEbadlybehaved);

double IacoCesare(double *x, double *p);
//double ScaleIacoCesare(double *p, int scaling);
void rangeIacoCesare(int dim, int *index, double* range);
void infoIacoCesare(double *p, int *maxdim, int *CEbadlybehaved);
int checkIacoCesare(double *param, int timespacedim, SimulationType method);

double lgd1(double *x, double*p);
double Scalelgd1(double *p,int scaling);
double Dlgd1(double *x, double *p);
void rangelgd1(int dim, int *index, double* range);
void infolgd1(double *p, int *maxdim, int *CEbadlybehaved);

// nsst
double spacetime1(double *x,double *p); 
double TBM2spacetime1(double *x, double *p);
//double TBM3spacetime1(double *x, double *p);
double Dspacetime1(double *x, double *p);
int checkspacetime1(double *param, int timespacedim, SimulationType method);
void rangespacetime1(int dim, int *index, double* range);
void infospacetime1(double *p, int *maxdim, int *CEbadlybehaved);

 // nsst2
double spacetime2(double *x,double *p);
//double TBM3spacetime2(double *x, double *p);
double Dspacetime2(double *x, double *p);
int checkspacetime2(double *param, int timespacedim, SimulationType method);
void rangespacetime2(int dim, int *index, double* range);
void infospacetime2(double *p, int *maxdim, int *CEbadlybehaved);

double nugget(double *x, double *p);
double Scalenugget(double *p, int scaling);
void rangenugget(int dim, int *index, double* range);
void infonugget(double *p, int *maxdim, int *CEbadlybehaved);
int checknugget(double *param, int timespacedim, SimulationType method);

double penta(double *x, double *p);
double Scalepenta(double *p,int scaling);
double TBM3penta(double *x, double *p);
double Dpenta(double *x, double *p);
void rangepenta(int dim, int *index, double* range);
void infopenta(double *p, int *maxdim, int *CEbadlybehaved);

double power(double *x, double *p);
double Scalepower(double *p,int scaling);
double TBM2power(double *x, double *p);
double TBM3power(double *x, double *p);
double Dpower(double *x, double *p);
int checkpower(double *param, int timespacedim, SimulationType method);
void rangepower(int dim, int *index, double* range);
void infopower(double *p, int *maxdim, int *CEbadlybehaved);

// TODO: hier laesst sicherlich noch einiges machen!!
double qexponential(double *x,double *p);
double Scaleqexponential(double *p,int scaling);
double TBM3qexponential(double *x, double *p);
double Dqexponential(double *x,double *p);
void rangeqexponential(int dim, int *index, double* range);
void infoqexponential(double *p, int *maxdim, int *CEbadlybehaved);

double spherical(double *x, double *p);
double Scalespherical(double *p,int scaling);
double TBM2spherical(double *x, double *p);
double TBM3spherical(double *x, double *p);
double Dspherical(double *x, double *p);
void rangespherical(int dim, int *index, double* range);
void infospherical(double *p, int *maxdim, int *CEbadlybehaved);

double stable(double *x, double *p);
double Scalestable(double *p,int scaling);
double TBM3stable(double *x, double *p);
double Dstable(double *x, double *p);
double DDstable(double *x, double *p);
int checkstable(double *param, int timespacedim, SimulationType method);
void rangestable(int dim, int *index, double* range);
void infostable(double *p, int *maxdim, int *CEbadlybehaved);

double stableX(double *x, double *p);
double DstableX(double *x, double *p);

double SteinST1(double *x, double *p);
//double ScaleSteinST1(double *p, int scaling);
int kappasSteinST1(int dim);
void rangeSteinST1(int dim, int *index, double* range);
void infoSteinST1(double *p, int *maxdim, int *CEbadlybehaved);
int checkSteinST1(double *param, int timespacedim, SimulationType method);

void infoundefined(double *p, int *maxdim, int *CEbadlybehaved);
int checkundefined(double *param, int timespacedim, SimulationType method);
int checkOK(double *param, int timespacedim, SimulationType method);

// TODO: TBM3wave
double wave(double *x, double *p);
double Scalewave(double *p,int scaling);
double spectralwave(double *p );
void rangewave(int dim, int *index, double* range);
void infowave(double *p, int *maxdim, int *CEbadlybehaved);

double WhittleMatern(double *x, double *p);
double ScaleWhittleMatern(double *p,int scaling);
double spectralWhittleMatern(double *p );
double TBM2WhittleMatern(double *x, double *p);
double TBM3WhittleMatern(double *x, double *p);
double DWhittleMatern(double *x, double *p);
double DDWhittleMatern(double *x, double *p);
int checkWhittleMatern(double *param, int timespacedim, SimulationType method);
void rangeWhittleMatern(int dim, int *index, double* range);
void infoWhittleMatern(double *p, int *maxdim, int *CEbadlybehaved);


#endif /* CovFcts_H*/

