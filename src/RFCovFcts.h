#ifndef CovFcts_H
#define CovFcts_H 1

Real nugget(Real x,Real *p);
Real Scalenugget(Real *p, int scaling);
int checknugget(key_type *key);

Real stable(Real x,Real *p);
Real Scalestable(Real *p,int scaling);
Real TBM3stable(Real x, Real *p);
int checkstable(key_type *key);

Real Gauss(Real x, Real*p);
Real ScaleGauss(Real*p,int scaling);
Real TBM3Gauss(Real x, Real*p);
Real spectralGauss(Real *p END_WITH_GSLRNG);

Real exponential(Real x,Real *p);
Real Scaleexponential(Real *p,int scaling);
Real TBM3exponential(Real x, Real *p);

Real qexponential(Real x,Real *p);
Real Scaleqexponential(Real *p,int scaling);
//Real TBM3exponential(Real x, Real *p);
int checkqexponential(key_type *key);

Real holeeffect(Real x, Real*p);
Real Scaleholeeffect(Real *p,int scaling);
Real TBM3holeeffect(Real x, Real *p);
int checkholeeffect(key_type *key);
 

Real hyperbolic(Real x, Real*p);
// Scalehyperbolic(Real *p,int scaling); missing
int checkhyperbolic(key_type *key);
Real TBM3hyperbolic(Real x, Real*p);


Real trianglemodel(Real x,Real *p);
Real Scaletrianglemodel(Real *p,int scaling);

Real circular(Real x, Real *p);
Real Scalecircular(Real *p,int scaling);

Real spherical(Real x, Real *p);
Real Scalespherical(Real *p,int scaling);
Real TBM2spherical(Real x, Real *p);
Real TBM3spherical(Real x, Real *p);


Real power(Real x, Real *p);
Real Scalepower(Real *p,int scaling);
Real TBM2power(Real x, Real *p);
Real TBM3power(Real x, Real *p);
int checkpower(key_type *key);


Real WhittleMatern(Real x, Real *p);
Real ScaleWhittleMatern(Real *p,int scaling);
Real spectralWhittleMatern(Real *p END_WITH_GSLRNG);
Real TBM3WhittleMatern(Real x, Real *p);
int checkWhittleMatern(key_type *key);


Real Gneiting(Real x, Real *p);
Real ScaleGneiting(Real *p,int scaling);
Real TBM3Gneiting(Real x, Real *p);

Real Gneitingdiff(Real x, Real *p);
//Real ScaleGneitingdiff(Real *p,int scaling); missing
Real TBM3Gneitingdiff(Real x, Real *p);
int checkGneitingdiff(key_type *key);

Real genGneiting(Real x, Real *p);
Real TBM3genGneiting(Real x, Real *p);
int checkgenGneiting(key_type *key);

Real Cauchy(Real x, Real *p);
Real ScaleCauchy(Real *p,int scaling);
Real TBM2Cauchy(Real x, Real *p);
Real TBM3Cauchy(Real x, Real *p);
int checkCauchy(key_type *key);

Real generalisedCauchy(Real x, Real *p);
Real ScalegeneralisedCauchy(Real *p,int scaling);
Real TBM3generalisedCauchy(Real x, Real *p);
int checkgeneralisedCauchy(key_type *key);

Real Cauchytbm(Real x, Real *p);
//Real ScaleCauchytbm(Real *p,int scaling); missing
Real TBM3Cauchytbm(Real x, Real *p);
int checkCauchytbm(key_type *key);


Real Bessel(Real x,Real *p);
//Real ScaleBessel(Real *p,int scaling); missing
Real spectralBessel(Real *p END_WITH_GSLRNG);
int checkBessel(key_type *key);

Real wave(Real x, Real *p);
Real Scalewave(Real *p,int scaling);
Real spectralwave(Real *p END_WITH_GSLRNG);


Real expPLUScirc(Real x, Real *p);
Real ScaleexpPLUScirc(Real *p,int scaling);
int checkexpPLUScirc(key_type *key);


Real cubic(Real x, Real *p);
Real Scalecubic(Real *p,int scaling);
Real TBM3cubic(Real x, Real *p);

Real penta(Real x, Real *p);
Real Scalepenta(Real *p,int scaling);
Real TBM3penta(Real x, Real *p);


Real fractalBrownian(Real x, Real *p);
Real twodimfractalBrownianlocal(Real x, Real *p);
Real twodimfractalBrownianR(Real *p, Real x, Real *newp);
Real twodimfractalBrownianS(Real *p);
int checkfractalBrownian(key_type *key);

Real threedimfractalBrownianlocal(Real x, Real *p);
Real threedimfractalBrownianR(Real *p, Real x, Real *newp);
Real threedimfractalBrownianS(Real *p);
int checkthreedimfractalBrownian(key_type *key);

#endif /* CovFcts_H*/

/* End. */














