#ifndef CovFcts_H
#define CovFcts_H 1

SimulationType methodnugget(int spacedim, bool grid);

Real testCov(Real *x,Real *p, int dim);
SimulationType methodtest(int spacedim, bool grid);

Real exponential(Real *x,Real *p, int dim);
Real Scaleexponential(Real *p,int scaling);
Real TBM2exponential(Real *x, Real *p, int effectivedim);
Real TBM3exponential(Real *x, Real *p, int dim);
Real Dexponential(Real *x, Real *p, int dim);
Real spectralexponential(Real *p );
SimulationType methodexponential(int spacedim, bool grid);
void rangeexponential(int spatialdim, int *index, Real* range);

// TODO: hier laesst sicherlich noch einiges machen!!
Real qexponential(Real *x,Real *p, int dim);
Real Scaleqexponential(Real *p,int scaling);
Real Dqexponential(Real *x,Real *p, int effectivedim);
int checkqexponential(Real *param, int timespacedim, SimulationType method);
SimulationType methodqexponential(int spacedim, bool grid);
void rangeqexponential(int spatialdim, int *index, Real* range);

Real dampedcosine(Real *x, Real*p, int dim);
Real Scaledampedcosine(Real *p,int scaling);
Real TBM3dampedcosine(Real *x, Real *p, int dim);
Real Ddampedcosine(Real *x, Real *p, int dim);
int checkdampedcosine(Real *param, int timespacedim, SimulationType method);
SimulationType methoddampedcosine(int spacedim, bool grid);
void rangedampedcosine(int spatialdim, int *index, Real* range);

Real circular(Real *x, Real *p, int dim);
Real Scalecircular(Real *p,int scaling);
Real Dcircular(Real *x, Real *p, int dim);
SimulationType methodcircular(int spacedim, bool grid);
void rangecircular(int spatialdim, int *index, Real* range);

Real spherical(Real *x, Real *p, int dim);
Real Scalespherical(Real *p,int scaling);
Real TBM2spherical(Real *x, Real *p, int dim);
Real TBM3spherical(Real *x, Real *p, int dim);
Real Dspherical(Real *x, Real *p, int dim);
SimulationType methodspherical(int spacedim, bool grid);
void rangespherical(int spatialdim, int *index, Real* range);

Real power(Real *x, Real *p, int dim);
Real Scalepower(Real *p,int scaling);
Real TBM2power(Real *x, Real *p, int dim);
Real TBM3power(Real *x, Real *p, int dim);
Real Dpower(Real *x, Real *p, int dim);
int checkpower(Real *param, int timespacedim, SimulationType method);
SimulationType methodpower(int spacedim, bool grid);
void rangepower(int spatialdim, int *index, Real* range);

Real stable(Real *x,Real *p, int dim);
Real Scalestable(Real *p,int scaling);
Real TBM3stable(Real *x, Real *p, int dim);
Real Dstable(Real *x, Real *p, int dim);
int checkstable(Real *param, int timespacedim, SimulationType method);
SimulationType methodstable(int spacedim, bool grid);
void rangestable(int spatialdim, int *index, Real* range);

Real WhittleMatern(Real *x, Real *p, int dim);
Real ScaleWhittleMatern(Real *p,int scaling);
Real spectralWhittleMatern(Real *p );
Real TBM2WhittleMatern(Real *x, Real *p, int dim);
Real TBM3WhittleMatern(Real *x, Real *p, int dim);
Real DWhittleMatern(Real *x, Real *p, int dim);
int checkWhittleMatern(Real *param, int timespacedim, SimulationType method);
SimulationType methodWhittleMatern(int spacedim, bool grid);
void rangeWhittleMatern(int spatialdim, int *index, Real* range);

Real hyperbolic(Real *x, Real*p, int dim);
// Scalehyperbolic(Real *p,int scaling); missing
int checkhyperbolic(Real *param, int timespacedim, SimulationType method);
Real TBM3hyperbolic(Real *x, Real*p, int dim);
Real Dhyperbolic(Real *x, Real*p, int dim);
SimulationType methodhyperbolic(int spacedim, bool grid);
void rangehyperbolic(int spatialdim, int *index, Real* range);

Real Gneiting(Real *x, Real *p, int dim);
Real ScaleGneiting(Real *p,int scaling);
Real TBM3Gneiting(Real *x, Real *p, int dim);
Real DGneiting(Real *x, Real *p, int dim);
SimulationType methodGneiting(int spacedim, bool grid);
void rangeGneiting(int spatialdim, int *index, Real* range);

Real genGneiting(Real *x, Real *p, int dim);
Real TBM3genGneiting(Real *x, Real *p, int dim);
Real DgenGneiting(Real *x, Real *p, int dim);
int checkgenGneiting(Real *param, int timespacedim, SimulationType method);
SimulationType methodgenGneiting(int spacedim, bool grid);
void rangegenGneiting(int spatialdim, int *index, Real* range);

Real Gauss(Real *x, Real*p, int dim);
Real ScaleGauss(Real*p,int scaling);
Real TBM3Gauss(Real *x, Real*p, int dim);
Real DGauss(Real *x, Real*p, int dim);
Real spectralGauss(Real *p );
SimulationType methodGauss(int spacedim, bool grid);
void rangeGauss(int spatialdim, int *index, Real* range);

Real Cauchy(Real *x, Real *p, int dim);
Real ScaleCauchy(Real *p,int scaling);
Real TBM2Cauchy(Real *x, Real *p, int dim);
Real TBM3Cauchy(Real *x, Real *p, int dim);
Real DCauchy(Real *x, Real *p, int dim);
int checkCauchy(Real *param, int timespacedim, SimulationType method);
SimulationType methodCauchy(int spacedim, bool grid);
void rangeCauchy(int spatialdim, int *index, Real* range);

Real generalisedCauchy(Real *x, Real *p, int dim);
Real ScalegeneralisedCauchy(Real *p,int scaling);
Real TBM3generalisedCauchy(Real *x, Real *p, int dim);
Real DgeneralisedCauchy(Real *x, Real *p, int dim);
int checkgeneralisedCauchy(Real *param, int timespacedim, SimulationType method);
SimulationType methodgeneralisedCauchy(int spacedim, bool grid);
void rangegeneralisedCauchy(int spatialdim, int *index, Real* range);

Real Cauchytbm(Real *x, Real *p, int dim);
//Real ScaleCauchytbm(Real *p,int scaling, int dim); missing
Real TBM3Cauchytbm(Real *x, Real *p, int dim);
Real DCauchytbm(Real *x, Real *p, int dim);
int checkCauchytbm(Real *param, int timespacedim, SimulationType method);
SimulationType methodCauchytbm(int spacedim, bool grid);
void rangeCauchytbm(int spatialdim, int *index, Real* range);

// TODO: TBM3bessel
Real Bessel(Real *x,Real *p, int dim);
//Real ScaleBessel(Real *p,int scaling); missing
Real spectralBessel(Real *p );
int checkBessel(Real *param, int timespacedim, SimulationType method);
SimulationType methodBessel(int spacedim, bool grid);
void rangeBessel(int spatialdim, int *index, Real* range);

// TODO: TBM3wave
Real wave(Real *x, Real *p, int dim);
Real Scalewave(Real *p,int scaling);
Real spectralwave(Real *p );
SimulationType methodwave(int spacedim, bool grid);
void rangewave(int spatialdim, int *index, Real* range);

Real cubic(Real *x, Real *p, int dim);
Real Scalecubic(Real *p,int scaling);
Real TBM3cubic(Real *x, Real *p, int dim);
Real Dcubic(Real *x, Real *p, int dim);
SimulationType methodcubic(int spacedim, bool grid);
void rangecubic(int spatialdim, int *index, Real* range);

Real penta(Real *x, Real *p, int dim);
Real Scalepenta(Real *p,int scaling);
Real TBM3penta(Real *x, Real *p, int dim);
Real Dpenta(Real *x, Real *p, int dim);
SimulationType methodpenta(int spacedim, bool grid);
void rangepenta(int spatialdim, int *index, Real* range);

Real spacetime1(Real *x,Real *p, int dim);
Real TBM2spacetime1(Real *x, Real *p, int dim);
Real TBM3spacetime1(Real *x, Real *p, int dim);
Real Dspacetime1(Real *x, Real *p, int dim);
int checkspacetime1(Real *param, int timespacedim, SimulationType method);
SimulationType methodspacetime1(int spacedim, bool grid);
void rangespacetime1(int spatialdim, int *index, Real* range);

Real spacetime2(Real *x,Real *p, int dim);
Real TBM3spacetime2(Real *x, Real *p, int dim);
Real Dspacetime2(Real *x, Real *p, int dim);
int checkspacetime2(Real *param, int timespacedim, SimulationType method);
SimulationType methodspacetime2(int spacedim, bool grid);
void rangespacetime2(int spatialdim, int *index, Real* range);

Real spacetime3(Real *x,Real *p, int dim);

Real fractalBrownian(Real *x, Real *p, int dim);
Real twodimfractalBrownianlocal(Real *x, Real *p, int dim);
Real twodimfractalBrownianS(Real *p);
int check2dfractalBrownian(Real *param, int timespacedim, SimulationType method);
SimulationType method2dfractalBrownian(int spacedim, bool grid);
void range2dfractalBrownian(int spatialdim, int *index, Real* range);

Real threedimfractalBrownianlocal(Real *x, Real *p, int dim);
Real threedimfractalBrownianS(Real *p);
int check3dfractalBrownian(Real *param, int timespacedim, SimulationType method);
SimulationType method3dfractalBrownian(int spacedim, bool grid);
void range3dfractalBrownian(int spatialdim, int *index, Real* range);

Real fractGauss(Real *x,Real *p, int dim);
int checkfractGauss(Real *param, int timespacedim, SimulationType method);
SimulationType methodfractGauss(int spacedim, bool grid);
void rangefractGauss(int spatialdim, int *index, Real* range);

Real lgd1(Real *x, Real*p, int effectivedim);
Real Scalelgd1(Real *p,int scaling);
Real Dlgd1(Real *x, Real *p, int dim);
SimulationType methodlgd1(int spacedim, bool grid);
int checklgd1(Real *param, int timespacedim, SimulationType method);
void rangelgd1(int spatialdim, int *index, Real* range);

Real FD(Real *x,Real *p, int effectivedim);
SimulationType methodFD(int spacedim, bool grid);
int checkFD(Real *param, int timespacedim, SimulationType method);
void rangeFD(int spatialdim, int *index, Real* range);

#endif /* CovFcts_H*/

