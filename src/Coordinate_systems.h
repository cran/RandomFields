#ifndef Coordinates_H
#define Coordinates_H 1


void iso2iso(double *x, cov_model *cov, double *v);
void spaceDiso2iso(double *x, cov_model *cov, double *v);
void spiso2spiso(double *x, cov_model *cov, double *v);
void spaceDspiso2spiso(double *x, cov_model *cov, double *v);
void spacetime2iso(double *x, cov_model *cov, double *v);
void spaceDspiso2iso(double *x, cov_model *cov, double *v);
void Stat2iso(double *x, cov_model *cov, double *v);
void Stat2spacetime(double *x, cov_model *cov, double *v);
void Stat2Stat(double *x, cov_model *cov, double *v);
void Nonstat2iso(double *x, double *y, cov_model *cov, double *v);
void Nonstat2spacetime(double *x, double *y, cov_model *cov, double *v);
void Nonstat2Stat(double *x, double *y, cov_model *cov, double *v);
void logiso2iso(double *x, cov_model *cov, double *v, double *Sign);
void logspiso2spiso(double *x, cov_model *cov, double *v, double *Sign);
void logspacetime2iso(double *x, cov_model *cov, double *v, double *Sign);
void logStat2iso(double *x, cov_model *cov, double *v, double *Sign);
void logNonstat2iso(double *x, double *y, cov_model *cov, double *v,
		    double *Sign);
void logStat2spacetime(double *x, cov_model *cov, double *v, double *Sign);
void logNonstat2spacetime(double *x, double *y, cov_model *cov, double *v,
			  double *Sign);
void logStat2Stat(double *x, cov_model *cov, double *v, double *Sign);
void logNonstat2Stat(double *x, double *y, cov_model *cov, double *v, 
		     double *Sign);
void Nonstat2Nonstat(double *x, double *y, cov_model *cov, double *v);
void logNonstat2Nonstat(double *x, double *y, cov_model *cov, double *v, 
			double *Sign);

void EarthIso2EarthIso(double *x, cov_model *cov, double *v);
void logEarthIso2EarthIso(double *x, cov_model *cov, double *v, double *Sign);
void NonstatEarth2EarthIso(double *x, double *y, cov_model *cov, double *v);
void logNonstatEarth2EarthIso(double *x, double *y, cov_model *cov, double *v, 
			      double *Sign);
void Earth2Earth(double *x, cov_model *cov, double *v);
void logEarth2Earth(double *x, cov_model *cov, double *v, double *Sign);
void NonstatEarth2Earth(double *x, double *y, cov_model *cov, double *v);
void logNonstatEarth2Earth(double *x, double *y, cov_model *cov, double *v, 
			   double *Sign);

void EarthIso2SphereIso(double *x, cov_model *cov, double *v);
void logEarthIso2SphereIso(double *x, cov_model *cov, double *v, double *Sign);
void NonstatEarth2SphereIso(double *x, double *y, cov_model *cov, double *v);
void logNonstatEarth2SphereIso(double *x, double *y, cov_model *cov, double *v, 
			       double *Sign);
void Earth2Sphere(double *x, cov_model *cov, double *v);
void logEarth2Sphere(double *x, cov_model *cov, double *v, double *Sign);
void NonstatEarth2Sphere(double *x, double *y, cov_model *cov, double *v);
void logNonstatEarth2Sphere(double *x, double *y, cov_model *cov, double *v, 
			    double *Sign);

void SphereIso2SphereIso(double *x, cov_model *cov, double *v);
void logSphereIso2SphereIso(double *x, cov_model *cov, double *v, double *Sign);
void NonstatSphere2SphereIso(double *x, double *y, cov_model *cov,double *v);
void logNonstatSphere2SphereIso(double *x, double *y, cov_model *cov, 
				double *v, double *Sign);
void Sphere2Sphere(double *x, cov_model *cov, double *v);
void logSphere2Sphere(double *x, cov_model *cov, double *v, double *Sign);
void NonstatSphere2Sphere(double *x, double *y, cov_model *cov, double *v);
void logNonstatSphere2Sphere(double *x, double *y, cov_model *cov,
			     double *v, double *Sign);


void EarthKM2CartStat(double *x, cov_model *cov, double *v);
void logEarthKM2CartStat(double *x, cov_model *cov, double *v, double *Sign);
void EarthKM2Cart(double *x, double *y, cov_model *cov, double *v);
void logEarthKM2Cart(double *x, double *y, cov_model *cov, double *v,
		     double *Sign);
void EarthMiles2CartStat(double *x, cov_model *cov, double *v);
void logEarthMiles2CartStat(double *x, cov_model *cov, double *v,double*Sign);
void EarthMiles2Cart(double *x, double *y, cov_model *cov, double *v);
void logEarthMiles2Cart(double *x, double *y, cov_model *cov, double *v,
		     double *Sign);
int checkEarth(cov_model *cov);

void EarthKM2OrthogStat(double *x, cov_model *cov, double *v);
void logEarthKM2OrthogStat(double *x, cov_model *cov, double *v, double *Sign);
void EarthKM2Orthog(double *x, double *y, cov_model *cov, double *v);
void logEarthKM2Orthog(double *x, double *y, cov_model *cov, double *v,
		     double *Sign);
void EarthMiles2OrthogStat(double *x, cov_model *cov, double *v);
void logEarthMiles2OrthogStat(double *x, cov_model *cov, double *v,double*Sign);
void EarthMiles2Orthog(double *x, double *y, cov_model *cov, double *v);
void logEarthMiles2Orthog(double *x, double *y, cov_model *cov, double *v,
		     double *Sign);


void Earth2GnomonicStat(double *x, cov_model *cov, double *v);
void logEarth2GnomonicStat(double *x, cov_model *cov, double *v, double*Sign);
void Earth2Gnomonic(double *x, double *y, cov_model *cov, double *v);
void logEarth2Gnomonic(double *x, double *y, cov_model *cov, double *v,
		     double *Sign);


#define Mod(ZZ, modulus) ((ZZ) - floor((ZZ) / (modulus)) * (modulus))
double mod(double x, double modulus);
double lonmod(double x, double modulus); 
double latmod(double x, double modulus);


#define STATMODE_BASE(X, ZZ, lon, lat)				\
  X[0]=lonmod(ZZ[0], lon);					\
  X[1]=latmod(ZZ[1], lat)					\

void statmod2(double *x, double lon, double lat, double *y); // 2dim only


bool isIsotropic(isotropy_type iso);
bool isAnyIsotropic(isotropy_type iso);
bool isSpaceIsotropic(isotropy_type iso);
bool isZeroSpaceIsotropic(isotropy_type iso);
bool isVectorIsotropic(isotropy_type iso);
bool isSymmetric(isotropy_type iso);
bool isCartesian(isotropy_type iso);
bool isSpherical(isotropy_type iso);
bool isCylinder(isotropy_type iso);
bool isEarth(isotropy_type iso);
bool isAnySpherical(isotropy_type iso);
bool isAnySphericalIso(isotropy_type iso);
bool isPrevModelI(cov_fct *C);
bool isUnreduced(cov_fct *C);
bool isCoordinateSystem(isotropy_type iso);
bool atleastSpecialised(isotropy_type iso, isotropy_type as);
isotropy_type UpgradeToCoordinateSystem( isotropy_type iso);

bool equal_coordinate_system(isotropy_type iso1, isotropy_type iso2);
bool equal_coordinate_system(isotropy_type iso1, isotropy_type iso2, 
			     bool refined);
bool is_any(isotropy_type iso, cov_fct *C);
bool is_all(isotropy_type iso, cov_fct *C);
typedef bool (*isofct)(isotropy_type iso); /* h, cov, result */ 
bool is_any(isofct iso, cov_fct *C);
bool is_all(isofct iso, cov_fct *C);

#define MAX_SYMMETRIC(iso) (iso) == CARTESIAN_COORD ? SYMMETRIC : (iso)

#endif
