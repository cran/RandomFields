#ifndef CovFcts_H
#define CovFcts_H 1


double constant(double *x, double *p, int effectivedim);
double TBM2constant(double *x, double *p, int effectivedim);
double TBM3constant(double *x, double *p, int effectivedim);
double Dconstant(double *x, double *p, int effectivedim);
void rangeconstant(int dim, int *index, double* range);
int checkconstant(double *param, int timespacedim, SimulationType method);
void infoconstant(double *p, int *maxdim, int *CEbadlybehaved);

double exponential(double *x,double *p, int dim);
double Scaleexponential(double *p,int scaling);
double TBM2exponential(double *x, double *p, int effectivedim);
double TBM3exponential(double *x, double *p, int dim);
double Dexponential(double *x, double *p, int dim);
double spectralexponential(double *p );
int hyperexponential(double radius, double *center, double *rx,
		     int dim, bool simulate, 
		     double** hx, double** hy, double** hr);
void rangeexponential(int spatialdim, int *index, double* range);
void infoexponential(double *p, int *maxdim, int *CEbadlybehaved);
int checkexponential(double *param, int timespacedim, SimulationType method);

// TODO: hier laesst sicherlich noch einiges machen!!
double qexponential(double *x,double *p, int dim);
double Scaleqexponential(double *p,int scaling);
double TBM3qexponential(double *x, double *p, int effectivedim);
double Dqexponential(double *x,double *p, int effectivedim);
int checkqexponential(double *param, int timespacedim, SimulationType method);
void rangeqexponential(int spatialdim, int *index, double* range);
void infoqexponential(double *p, int *maxdim, int *CEbadlybehaved);

double dampedcosine(double *x, double*p, int dim);
double Scaledampedcosine(double *p,int scaling);
double TBM3dampedcosine(double *x, double *p, int dim);
double Ddampedcosine(double *x, double *p, int dim);
int checkdampedcosine(double *param, int timespacedim, SimulationType method);
void rangedampedcosine(int spatialdim, int *index, double* range);
void infodampedcosine(double *p, int *maxdim, int *CEbadlybehaved);

double circular(double *x, double *p, int dim);
double Scalecircular(double *p,int scaling);
double Dcircular(double *x, double *p, int dim);
int checkcircular(double *param, int timespacedim, SimulationType method);
void rangecircular(int spatialdim, int *index, double* range);
void infocircular(double *p, int *maxdim, int *CEbadlybehaved);

double spherical(double *x, double *p, int dim);
double Scalespherical(double *p,int scaling);
double TBM2spherical(double *x, double *p, int dim);
double TBM3spherical(double *x, double *p, int dim);
double Dspherical(double *x, double *p, int dim);
void rangespherical(int spatialdim, int *index, double* range);
void infospherical(double *p, int *maxdim, int *CEbadlybehaved);
int checkspherical(double *param, int timespacedim, SimulationType method);

double power(double *x, double *p, int dim);
double Scalepower(double *p,int scaling);
double TBM2power(double *x, double *p, int dim);
double TBM3power(double *x, double *p, int dim);
double Dpower(double *x, double *p, int dim);
int checkpower(double *param, int timespacedim, SimulationType method);
void rangepower(int spatialdim, int *index, double* range);
void infopower(double *p, int *maxdim, int *CEbadlybehaved);

double stable(double *x,double *p, int dim);
double Scalestable(double *p,int scaling);
double TBM3stable(double *x, double *p, int dim);
double Dstable(double *x, double *p, int dim);
double DDstable(double *x, double *p, int effectivedim);
int checkstable(double *param, int timespacedim, SimulationType method);
void rangestable(int spatialdim, int *index, double* range);
void infostable(double *p, int *maxdim, int *CEbadlybehaved);

double WhittleMatern(double *x, double *p, int dim);
double ScaleWhittleMatern(double *p,int scaling);
double spectralWhittleMatern(double *p );
double TBM2WhittleMatern(double *x, double *p, int dim);
double TBM3WhittleMatern(double *x, double *p, int dim);
double DWhittleMatern(double *x, double *p, int dim);
double DDWhittleMatern(double *x, double *p, int effectivedim);
int checkWhittleMatern(double *param, int timespacedim, SimulationType method);
void rangeWhittleMatern(int spatialdim, int *index, double* range);
void infoWhittleMatern(double *p, int *maxdim, int *CEbadlybehaved);

double hyperbolic(double *x, double*p, int dim);
// Scalehyperbolic(double *p,int scaling); missing
int checkhyperbolic(double *param, int timespacedim, SimulationType method);
double TBM3hyperbolic(double *x, double*p, int dim);
double Dhyperbolic(double *x, double*p, int dim);
void rangehyperbolic(int spatialdim, int *index, double* range);
void infohyperbolic(double *p, int *maxdim, int *CEbadlybehaved);

double Gneiting(double *x, double *p, int dim);
double ScaleGneiting(double *p,int scaling);
double TBM3Gneiting(double *x, double *p, int dim);
double DGneiting(double *x, double *p, int dim);
int checkGneiting(double *param, int timespacedim, SimulationType method);
void rangeGneiting(int spatialdim, int *index, double* range);
void infoGneiting(double *p, int *maxdim, int *CEbadlybehaved);

double genGneiting(double *x, double *p, int dim);
double TBM3genGneiting(double *x, double *p, int dim);
double DgenGneiting(double *x, double *p, int dim);
int checkgenGneiting(double *param, int timespacedim, SimulationType method);
void rangegenGneiting(int spatialdim, int *index, double* range);
void infogenGneiting(double *p, int *maxdim, int *CEbadlybehaved);

double Gauss(double *x, double*p, int dim);
double ScaleGauss(double*p,int scaling);
double TBM3Gauss(double *x, double*p, int dim);
double DGauss(double *x, double*p, int dim);
double spectralGauss(double *p );
void rangeGauss(int spatialdim, int *index, double* range);
void infoGauss(double *p, int *maxdim, int *CEbadlybehaved);
int checkGauss(double *param, int timespacedim, SimulationType method);

double Cauchy(double *x, double *p, int dim);
double ScaleCauchy(double *p,int scaling);
double TBM2Cauchy(double *x, double *p, int dim);
double TBM3Cauchy(double *x, double *p, int dim);
double DCauchy(double *x, double *p, int dim);
int checkCauchy(double *param, int timespacedim, SimulationType method);
void rangeCauchy(int spatialdim, int *index, double* range);
void infoCauchy(double *p, int *maxdim, int *CEbadlybehaved);

double generalisedCauchy(double *x, double *p, int dim);
double ScalegeneralisedCauchy(double *p,int scaling);
double TBM3generalisedCauchy(double *x, double *p, int dim);
double DgeneralisedCauchy(double *x, double *p, int dim);
double DDgeneralisedCauchy(double *x, double *p, int effectivedim);
int checkgeneralisedCauchy(double *param, int timespacedim, 
			   SimulationType method);
void rangegeneralisedCauchy(int spatialdim, int *index, double* range);
void infogeneralisedCauchy(double *p, int *maxdim, int *CEbadlybehaved);

double Cauchytbm(double *x, double *p, int dim);
//double ScaleCauchytbm(double *p,int scaling, int dim); missing
double TBM3Cauchytbm(double *x, double *p, int dim);
double DCauchytbm(double *x, double *p, int dim);
int checkCauchytbm(double *param, int timespacedim, SimulationType method);
void rangeCauchytbm(int spatialdim, int *index, double* range);
void infoCauchytbm(double *p, int *maxdim, int *CEbadlybehaved);

// TODO: TBM3bessel
double Bessel(double *x,double *p, int dim);
//double ScaleBessel(double *p,int scaling); missing
double spectralBessel(double *p );
int checkBessel(double *param, int timespacedim, SimulationType method);
void rangeBessel(int spatialdim, int *index, double* range);
void infoBessel(double *p, int *maxdim, int *CEbadlybehaved);

// TODO: TBM3wave
double wave(double *x, double *p, int dim);
double Scalewave(double *p,int scaling);
double spectralwave(double *p );
void rangewave(int spatialdim, int *index, double* range);
void infowave(double *p, int *maxdim, int *CEbadlybehaved);
int checkwave(double *param, int timespacedim, SimulationType method);

double cubic(double *x, double *p, int dim);
double Scalecubic(double *p,int scaling);
double TBM3cubic(double *x, double *p, int dim);
double Dcubic(double *x, double *p, int dim);
int checkcubic(double *param, int timespacedim, SimulationType method);
void rangecubic(int spatialdim, int *index, double* range);
void infocubic(double *p, int *maxdim, int *CEbadlybehaved);

double penta(double *x, double *p, int dim);
double Scalepenta(double *p,int scaling);
double TBM3penta(double *x, double *p, int dim);
double Dpenta(double *x, double *p, int dim);
void rangepenta(int spatialdim, int *index, double* range);
void infopenta(double *p, int *maxdim, int *CEbadlybehaved);
int checkpenta(double *param, int timespacedim, SimulationType method);

double spacetime1(double *x,double *p, int dim);
double TBM2spacetime1(double *x, double *p, int dim);
double TBM3spacetime1(double *x, double *p, int dim);
double Dspacetime1(double *x, double *p, int dim);
int checkspacetime1(double *param, int timespacedim, SimulationType method);
void rangespacetime1(int spatialdim, int *index, double* range);
void infospacetime1(double *p, int *maxdim, int *CEbadlybehaved);

double spacetime2(double *x,double *p, int dim);
double TBM3spacetime2(double *x, double *p, int dim);
double Dspacetime2(double *x, double *p, int dim);
int checkspacetime2(double *param, int timespacedim, SimulationType method);
void rangespacetime2(int spatialdim, int *index, double* range);
void infospacetime2(double *p, int *maxdim, int *CEbadlybehaved);

double spacetime3(double *x,double *p, int dim);


double fractalBrownian(double*x, double *p, int effectivdim);
double DfractalBrownian(double *x, double*p, int effectivedim);
double DDfractalBrownian(double *x, double*p, int effectivedim);
int checkfractalBrownian(double *param, int timespacedim, SimulationType method);
void rangefractalBrownian(int dim, int *index, double* range);
void infofractalBrownian(double *p, int *maxdim, int *CEbadlybehaved);

double fractGauss(double *x,double *p, int dim);
int checkfractGauss(double *param, int timespacedim, SimulationType method);
void rangefractGauss(int spatialdim, int *index, double* range);
void infofractGauss(double *p, int *maxdim, int *CEbadlybehaved);

double lgd1(double *x, double*p, int effectivedim);
double Scalelgd1(double *p,int scaling);
double Dlgd1(double *x, double *p, int dim);
int checklgd1(double *param, int timespacedim, SimulationType method);
void rangelgd1(int spatialdim, int *index, double* range);
void infolgd1(double *p, int *maxdim, int *CEbadlybehaved);

double FD(double *x,double *p, int effectivedim);
int checkFD(double *param, int timespacedim, SimulationType method);
void rangeFD(int spatialdim, int *index, double* range);
void infoFD(double *p, int *maxdim, int *CEbadlybehaved);

double IacoCesare(double *x, double *p, int dim);
//double ScaleIacoCesare(double *p, int scaling);
void rangeIacoCesare(int dim, int *index, double* range);
void infoIacoCesare(double *p, int *maxdim, int *CEbadlybehaved);
int checkIacoCesare(double *param, int timespacedim, SimulationType method);

double nugget(double *x, double *p, int dim);
double Scalenugget(double *p, int scaling);
void rangenugget(int dim, int *index, double* range);
void infonugget(double *p, int *maxdim, int *CEbadlybehaved);
int checknugget(double *param, int timespacedim, SimulationType method);

double SteinST1(double *x, double *p, int dim);
//double ScaleSteinST1(double *p, int scaling);
int kappasSteinST1(int dim);
void rangeSteinST1(int dim, int *index, double* range);
void infoSteinST1(double *p, int *maxdim, int *CEbadlybehaved);
int checkSteinST1(double *param, int timespacedim, SimulationType method);

double nugget(double *x, double *p, int dim);
double Scalenugget(double *p, int scaling);
void rangenugget(int dim, int *index, double* range);
void infonugget(double *p, int *maxdim, int *CEbadlybehaved);
int checknugget(double *param, int timespacedim, SimulationType method);

void infoundefined(double *p, int *maxdim, int *CEbadlybehaved);
int checkundefined(double *param, int timespacedim, SimulationType method);

#endif /* CovFcts_H*/

