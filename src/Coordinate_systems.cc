/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2014 -- 2017 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/



#include "def.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include "extern.h"

#include <zzz_RandomFieldsUtils.h>
#include "RF.h"
#include "Coordinate_systems.h"
#include "questions.h"
#include "kleinkram.h"
extern const char *internals[internalN];
extern const char *general[generalN];
extern const char *coords[coordsN];

/* simplifying functions
   turn vector of x into length ||x|| and similar
   reduces C(x,y) to C(x - y)

   needs preceeding analysis of submodels!
   (same for preferred methods; bottom-up-analysis needed)
*/

// keep always the orderung

#define MAXEARTHXDIM 6

bool isSameCoordSystem(isotropy_type iso, coord_sys_enum os) {
  switch(os) {
  case cartesian : case gnomonic: case orthographic : return isCartesian(iso);
  case earth : return isEarth(iso);
  case sphere : return isSpherical(iso);
  case coord_mix : return true;
  default: BUG;
  }
}



coord_sys_enum SearchCoordSystem(model *cov, coord_sys_enum os, 
				 coord_sys_enum n_s) {
  if (n_s == coord_keep) {
    if (!isSameCoordSystem(OWNISO(0), os)) {
      n_s = GetCoordSystem(OWNISO(0));
    }
  } else {
    if (n_s == coord_mix || !isSameCoordSystem(OWNISO(0), n_s)) {
      return coord_mix;
    }
  }
  int i;
  coord_sys_enum nn_s;
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      nn_s = SearchCoordSystem(cov->sub[i], os, n_s);
      if (nn_s != n_s) {
	if (n_s == coord_keep) n_s = nn_s;
	else {
	  return coord_mix;
	}
      }
    }
  }
  return n_s;
}


//////////////////////////////////////////////////////////////////////
// Earth coordinate systems : stay within the system
//////////////////////////////////////////////////////////////////////
bool hasEarthHeight(system_type *sys) {
  return false; // to do !!
  return LASTSYSTEM(sys) > 0 && isLogCart(sys, 1) && LOGDIM(sys, 1) == 1;
}

double mod(double x, double modulus) { return (x - FLOOR(x / modulus) * modulus); }
double isomod(double x, double modulus) {  
  double 
    twomodulus = 2.0 * modulus;
  return modulus - FABS(Mod(x, twomodulus) - modulus);
}


double latmod(double x, double modulus) {  
  double 
    halfmodulus = 0.5 * modulus,
    twomodulus = 2.0 * modulus,
    y = x - halfmodulus;
  return FABS(Mod(y, twomodulus) - modulus) - halfmodulus;
}


#define STATMOD(ZZ, X, lon, lat)				\
    STATMODE_BASE(X, ZZ, lon, lat);				\
    for (int d = 2; d < dim; d++) X[d] = ZZ[d]

void statmod2(double *x, double lon, double lat, double *y) {
  //STATMOD(y, x, lon, lat);
  STATMODE_BASE(x, y, lon, lat);
}
	

#define ISOMOD(ZZ, X, maxangle)				\
  int dim=PREVTOTALXDIM;					\
  X[0]=isomod(ZZ[0], maxangle);					\
  for (int d = 1; d < dim; d++) X[d] = ZZ[d]




#define ESnonstat2iso(x, y)						\
  int last=PREVLASTSYSTEM;						\
  double xx[2] = {x[0], x[1]},						\
    yy[2] = {y[0], y[1]},						\
    cosine = SIN(xx[1]) * SIN(yy[1]) +					\
             (COS(xx[0]) * COS(yy[0]) + SIN(xx[0]) * SIN(yy[0])) *      \
             COS(xx[1]) * COS(yy[1]);				       	\
    cosine = cosine <= 1.0 ? (cosine < -1.0 ? -1.0 : cosine) : 1.0;	\
    X[0] = ACOS(cosine);						\
    for (int s=1; s<last; s++) {					\
      int dim = PREVXDIM(s),						\
	base = PREVCUMXOHNE(s);					\
      if (isCartesian(PREVISO(s)))					\
	for (int d=0; d<dim; d++) X[base + d - 1] = x[base + d] - y[base + d]; \
      else if (isLogCart(PREVISO(s)))					\
	for (int d=0; d<dim; d++) X[base + d - 1] = x[base + d] / y[base + d]; \
    }


void EarthIso2EarthIso(double *x) { x[0]=isomod(x[0], 180); }
void Earth2Earth(double *x) { STATMODE_BASE(x, x, 360, 180); }

void NonstatEarth2EarthIso(double *x, double *y, model *cov, 
			   double *X) {
  ESnonstat2iso(piD180 * x, piD180 * y);
  X[0] *= H80Dpi;
}


void NonstatEarth2Earth(double *x, double *y, model *cov, 
			double *X, double *Y) {
  int dim=PREVTOTALXDIM;					
  STATMOD(x, X, 360, 180);					
  STATMOD(y, Y, 360, 180);					
}


void Sphere2Earth(double *x, model *cov, double *X) {
  int dim=PREVTOTALXDIM;			       
  STATMOD(H80Dpi * x, X, 360, 180);
}


void EarthIso2SphereIso(double *x, model *cov, double *X) {
  ISOMOD(piD180 * x, X, M_PI);
}

void NonstatEarth2SphereIso(double *x, double *y, model *cov, double *X) {
  ESnonstat2iso(piD180 * x, piD180 * y);
}


void Earth2Sphere(double *x, model *cov, double *X) {
  int dim=PREVTOTALXDIM;			       
  STATMODE_BASE(X, piD180 * x, M_2_PI, M_PI);
  for (int d = 2; d < dim; d++) X[d] = x[d];
}

void NonstatEarth2Sphere(double *x, double *y, model *cov, 
			 double *X, double *Y) {
  int dim=PREVTOTALXDIM;	
  STATMOD(piD180 * x, X, M_2_PI, M_PI);
  STATMOD(piD180 * y, Y, M_2_PI, M_PI);
}

void SphereIso2SphereIso(double *x) {x[0]=isomod(x[0], M_PI); }
void Sphere2Sphere(double *x) { STATMODE_BASE(x, x, M_2_PI, M_PI); }

void NonstatSphere2SphereIso(double *x, double *y, model *cov,
			     double *X) {
  ESnonstat2iso(x, y);
}

void NonstatSphere2Sphere(double *x, double *y, model *cov, 
			  double *X, double *Y) {
  int dim=PREVTOTALXDIM;			       
  STATMOD(x, X, M_2_PI, M_PI);
  STATMOD(y, Y, M_2_PI, M_PI); 
}



//////////////////////////////////////////////////////////////////////
// Earth coordinate systems : change to cartesian system
//////////////////////////////////////////////////////////////////////



#define EARTH_LONGITUDE 0
#define EARTH_LATITUDE 1
#define EARTH_HEIGHT 2
#define pi180 0.017453292519943295474
#define EARTH_TRAFO(X, ZZ, raequ, rpol)			\
  Rcos = (raequ) * COS(ZZ[EARTH_LATITUDE] * pi180);	\
  X[0] = Rcos * COS(ZZ[EARTH_LONGITUDE] * pi180);	\
  X[1] = Rcos * SIN(ZZ[EARTH_LONGITUDE] * pi180);	\
  X[2] = (rpol) * SIN(ZZ[EARTH_LATITUDE] * pi180)		

#if MAXSYSTEMS == 1
#define ASSERT_EARTH /*not necessarily equal*/	\
assert(PREVXDIM(0) == 2 || PREVXDIM(0) == 3);	\
   assert(GATTERLASTSYSTEM == 0)
#else
#define ASSERT_EARTH  assert(PREVXDIM(0) == 2 || PREVXDIM(0) == 3); \
   assert(GATTERLASTSYSTEM == 0)
#endif


// earthIso2cartIso ist direkt moeglich, aber nicht programmiert


#define earth2cartInner(RAEQU, RPOL)					\
  ASSERT_EARTH;								\
  int					\
    dim = PREVTOTALXDIM,						\
    base = 2;							\
  double Rcos;								\
  if (hasEarthHeight(PREV)) {		BUG;				\
    base++;								\
    EARTH_TRAFO(X, x, (RAEQU) + x[EARTH_HEIGHT],( RPOL) + x[EARTH_HEIGHT]); \
    EARTH_TRAFO(Y, y, (RAEQU) + y[EARTH_HEIGHT], (RPOL) + y[EARTH_HEIGHT]); \
  } else {								\
    EARTH_TRAFO(X, x, RAEQU, RPOL);					\
    EARTH_TRAFO(Y, y, RAEQU, RPOL);					\
  }									\
  for (int d=base; d<dim; d++) {					\
    X[d - base + 3] = x[d];						\
    Y[d - base + 3] = y[d];						\
  }


 // !!! if anything is changed here change also structtrafoproc !!!!


#define earth2cartInnerStat(RAEQU, RPOL)				\
  ASSERT_EARTH;								\
  int									\
    dim = PREVTOTALXDIM,						\
    base = 2;								\
  double Rcos;								\
  if (hasEarthHeight(PREV)) { BUG; /* height */				\
    base++;								\
    EARTH_TRAFO(X, x, (RAEQU) + x[EARTH_HEIGHT],( RPOL) + x[EARTH_HEIGHT]); \
  } else {								\
    EARTH_TRAFO(X, x, RAEQU, RPOL);					\
  }									\
  for (int d=base; d<dim; d++) X[d - base + 3] = x[d]


void Earth2Cart(model *cov, double RAEQU, double RPOL, double *yy){
  ASSERT_EARTH;								
  location_type *loc = Loc(cov);
  bool height = hasEarthHeight(PREV);
  double  Rcos, X[4],
    *x = loc->x;
  int spatialdim = loc->spatialdim,
    spatialpts = loc->spatialtotalpoints,
    base = 2 + (int) height,
    restdim = spatialdim - base,
    bytes = restdim * sizeof(double);
  if (height) { /* height */				
    for (int i=0; i<spatialpts; i++, x+=spatialdim) {
      EARTH_TRAFO(X, x, (RAEQU) + x[EARTH_HEIGHT],( RPOL) + x[EARTH_HEIGHT]);
      MEMCOPY(yy, X, 3 *sizeof(double));
      yy += 3;
      if (restdim > 0) {
	MEMCOPY(yy, x + base, bytes);
	yy += restdim;
      }
    }
  } else {								 
    for (int i=0; i<spatialpts; i++, x+=spatialdim) {
      EARTH_TRAFO(X, x, RAEQU, RPOL);					
      MEMCOPY(yy, X, 3 *sizeof(double));
      yy += 3;
      if (restdim > 0) {
	MEMCOPY(yy, x + base, bytes);
	yy += restdim;
      }
    }
  }								
}

void EarthKM2CartStat(double *x, model *cov, double *X) {
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
}

void EarthKM2Cart(double *x, double *y, model *cov, double *X, double *Y) {
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
}


void EarthMiles2CartStat(double *x, model *cov, double *X) {
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
}

void EarthMiles2Cart(double *x, double *y, model *cov, double *X, double*Y){
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
}

#define QP 0
#define QZenit 9
#define Qtotal 12
int checkEarth(model *cov){
  // ACHTUNG! KEIN AUFRUF VON SUB[0] !
  if (equalsXonly(PREVDOM(0)) // 20.2.14: warum war es vorher OWNDOM(0) == X ONLY?
      && isSymmetric(PREVISO(0))) {
    //darf nie auf "X ONLY-stat" transformiert sein. Es kann aber zB bei 'shape'
    // sein, dass ein X ONLY ankommt.
    RETURN_ERR(ERRORKERNEL);
  }
  
  //  COND_NEW_STORAGE(earth, X);
  ONCE_NEW_STORAGE(earth);

  if (!isEarth(PREVISO(0))) SERR("earth system expected in first component");
  if (TRAFONR >= FIRST_PLANE && TRAFONR <= LAST_PLANE) {
    assert(GATTERXDIM(0) == 2 && GATTERXDIM(0) == GATTERLOGDIM(0));
    if (!R_FINITE(GLOBAL.coords.zenit[0]) || !R_FINITE(GLOBAL.coords.zenit[1])){
      //if (GLOBAL.internal.warn_zenit) {
	//	GLOBAL.internal.warn_zenit = false;
	//	WARN1("tried to use non-finite values of '%.50s' in a coordinate transformation. Is the zenit set?\n", coords[ZENIT]);
      //}
      SERR1("Tried to use non-finite values of '%.50s' in a coordinate transformation. Is the zenit set?", coords[ZENIT]);
    }
    
    if (TRAFONR == EARTHKM2GNOMONIC || 
	TRAFONR == EARTHMILES2GNOMONIC) {
     double Rcos, X[4];
     assert(cov->Searth != NULL);
     //if (cov->Searth == NULL) NEW_STORAGE(earth);
      if (TRAFONR == EARTHKM2GNOMONIC) {
	EARTH_TRAFO(X, GLOBAL.coords.zenit, radiuskm_aequ, radiuskm_pol); 
      } else {
 	EARTH_TRAFO(X, GLOBAL.coords.zenit, radiusmiles_aequ, radiusmiles_pol);
      }
      double Rsq = 0.0,
	*zenit = cov->Searth->cart_zenit;
      for (int d=0; d<=2; d++) Rsq += X[d] * X[d];
      for (int d=0; d<=2; d++) {
	zenit[d] = X[d] / Rsq; // achtung! nicht s q r t(Rsq) !
      }
    }
    double Zenit[2] = { GLOBAL.coords.zenit[0] * pi180, 
			GLOBAL.coords.zenit[1] * pi180},
      sin0 = SIN(Zenit[0]),
	sin1 = SIN(Zenit[1]),
	cos0 = COS(Zenit[0]),
	cos1 = COS(Zenit[1]),
	*P = cov->Searth->P;
     

    P[0] = -sin0;
    P[1] = cos0;
    P[2] = 0.0;
    P[3] = -cos0 * sin1;
    P[4] = -sin0 * sin1;
    P[5] = cos1;
    P[6] = cos0 * cos1; 
    P[7] = sin0 * cos1;
    P[8] = sin1;    

  } else {
    assert(TRAFONR >= EARTHKM2CART && TRAFONR <= EARTHMILES2CART);
  }

  if (GATTERTOTALXDIM > MAXEARTHXDIM) SERR("dimension exceeded");
    
  RETURN_NOERROR;
}




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Orthogonal
#define orthDefStat						\
  double X[MAXEARTHXDIM + 1];					\
  double *P = cov->Searth->P; assert(P!=NULL);
 

#define orthTrafoStat							\
  int m = 0;								\
  U[0] = 0.0;    for(int k=0; k<3; k++, m++) U[0] += P[m] * X[k];	\
  U[1] = 0.0;    for(int k=0; k<3; k++, m++) U[1] += P[m] * X[k];	\
  double U2 = 0; for(int k=0; k<3; k++, m++) U2   += P[m] * X[k];	\
  assert(m==9 && dim <= 2);						\
  if (U2 < 0.0) ERR("location(s) not in direction of the zenit");	\
  for (int d = 2; d<dim; d++) U[d] = x[d]

#define orthDef orthDefStat;			\
  double Y[MAXEARTHXDIM + 1];

#define orthTrafo							\
  int m=0;								\
  U[0] = V[0] = 0.0;							\
  for(int k=0; k<3; k++, m++) { U[0] += P[m] * X[k]; V[0] += P[m] * Y[k]; } \
  U[1] = V[1] = 0.0;							\
  for(int k=0; k<3; k++, m++) { U[1] += P[m] * X[k]; V[1] += P[m] * Y[k]; } \
  double U2 = 0.0, V2 = 0.0;						\
  for(int k=0; k<3; k++, m++) { U2   += P[m] * X[k]; V2   += P[m] * Y[k]; } \
  if (U2 < 0.0 || V2 < 0.0) ERR("location(s) not in direction of the zenit"); \
  for (int d = 2; d<dim; d++) {U[d] = x[d]; V[d] = y[d];}


void EarthKM2OrthogStat(double *x, model *cov, double *U) {
  assert(TRAFONR == EARTHKM2ORTHOGRAPHIC);
  orthDefStat;
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol); // results in X[0], X[1]
  orthTrafoStat;
 }

  void EarthKM2Orthog(double *x, double *y, model *cov, double *U, 
		      double *V) {
  assert(TRAFONR == EARTHKM2ORTHOGRAPHIC);
  orthDef;
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthTrafo;
}
 

void EarthMiles2OrthogStat(double *x, model *cov, double *U) {
  assert(TRAFONR == EARTHMILES2ORTHOGRAPHIC);
  orthDefStat;
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);//results in X[0], X[1]
  orthTrafoStat;
}

 void EarthMiles2Orthog(double *x, double *y, model *cov, 
			double *U, double *V) {
  assert(TRAFONR == EARTHMILES2ORTHOGRAPHIC);
   orthDef;
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  orthTrafo;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Gnomonic

#define GnomonicStat						\
  double factor = 0.0,						\
    *cart_zenit = cov->Searth->cart_zenit;				\
  for (int d=0; d<3; d++) factor += cart_zenit[d] * X[d];		\
  if (factor <= 0.0)						\
    ERR1("locations not on the half-sphere given by the '%.50s'.", \
	 coords[ZENIT]);					   \
  for (int d=0; d<3; d++) X[d] /= factor

#define Gnomonic							\
  double factor = 0.0,							\
    factorY = 0.0,							\
    *cart_zenit = cov->Searth->cart_zenit;				\
  for (int d=0; d<3; d++) {						\
    factor += cart_zenit[d] * X[d];					\
    factorY += cart_zenit[d] * Y[d];					\
  }									\
  if (factor <= 0.0 || factorY <= 0.0) 					\
    ERR1("locations not on the half-sphere given by the '%.50s'.",\
	 coords[ZENIT]);						\
  for (int d=0; d<3; d++) {						\
    X[d] /= factor;							\
    Y[d] /= factorY;							\
  }



void Earth2GnomonicStat(double *x, model *cov, double *U) {
  assert(TRAFONR == EARTHKM2GNOMONIC || TRAFONR == EARTHMILES2GNOMONIC);
  orthDefStat; 
  earth2cartInnerStat(1.0, radiuskm_pol / radiuskm_aequ);//Umwandlung in CartSys
  GnomonicStat; // Projektion des Punktes auf die Ebene, die im R^3 liegt
  orthTrafoStat; // Projektion dieser Ebene auf R^2
}

 void Earth2Gnomonic(double *x, double *y, model *cov, double *U, double*V){
  assert(TRAFONR == EARTHKM2GNOMONIC || TRAFONR == EARTHMILES2GNOMONIC);
  orthDef;
  earth2cartInner(1.0, radiuskm_pol / radiuskm_aequ);
  Gnomonic;
  orthTrafo;
}




coord_sys_enum GetCoordSystem(isotropy_type iso) {
  if (isCartesian(iso)) return cartesian;
  if (isEarth(iso)) return earth;
  if (isSpherical(iso)) return sphere;
  return coord_mix;
}



SEXP GetCoordSystem(SEXP keynr, SEXP oldsystem, SEXP newsystem) {
  int knr = INTEGER(keynr)[0];
  model *cov;
  SEXP res;
  char CS[2][30] = {"coordinate system", "new coordinate system"};
  model **key = KEY();

  if (knr>=0 && knr <= MODEL_MAX && key[knr] != NULL) {
    cov = key[knr];
    coord_sys_enum
      os = (coord_sys_enum) GetName(oldsystem, CS[0], 
				    COORD_SYS_NAMES, nr_coord_sys, coord_auto),
      n_s = (coord_sys_enum) GetName(newsystem, CS[1], 
				    COORD_SYS_NAMES, nr_coord_sys, coord_keep);
   if (os == coord_auto) {
      os = GetCoordSystem(PREVISO(0));
    }
    if (n_s == coord_keep) {
      n_s = SearchCoordSystem(cov, os, n_s);
    }

    if (n_s == coord_mix && GLOBAL.internal.warn_coord_change) {    
      WARN1("the covariance model relies on at least two different coordinate systems. Use RFgetModelInfo(level=6) and check that this is not due to misspecification of the covariance model. To avoid this warning set 'RFoptions(%.50s=FALSE)'", // OK
	    internals[INTERNALS_COORD_CHANGE]);
      GLOBAL.internal.warn_coord_change = false;
    }

    bool warn = ((os != coord_auto && os != cartesian) ||
		 (n_s != coord_keep && n_s != os));

    switch (GLOBAL.general.reportcoord) {
    case reportcoord_none :  return R_NilValue;
    case reportcoord_warnings_orally : 
      if (warn) {
	WARN3("internal change of coordinate system from '%.50s' to '%.50s'. To avoid this message change the value of '%.50s' by 'RFoptions'.",
	      COORD_SYS_NAMES[os], COORD_SYS_NAMES[n_s], 
	      general[GENERAL_REPORTCOORD]);
      }
      return R_NilValue;
    case reportcoord_warnings :
      if (!warn) return R_NilValue;
      FALLTHROUGH_OK;
    case reportcoord_always :     
      PROTECT(res = allocVector(STRSXP, 2));
      SET_STRING_ELT(res, 0, mkChar(COORD_SYS_NAMES[os]));    
      SET_STRING_ELT(res, 1, mkChar(COORD_SYS_NAMES[n_s]));    
      UNPROTECT(1);
      break;
    default : BUG;
    }
    return res;
  }
  return R_NilValue;
}



isotropy_type CoordinateSystemOf(isotropy_type iso) {
  if (isCartesian(iso)) return CARTESIAN_COORD;
  if (isEarth(iso)) return EARTH_COORD;
  return isSpherical(iso) ? SPHERICAL_COORD : ISO_MISMATCH;
}

isotropy_type EssentialCoordinateSystemOf(isotropy_type iso) {
  if (isCartesian(iso)) return CARTESIAN_COORD;
  return isAnySpherical(iso) ? SPHERICAL_COORD : ISO_MISMATCH;
}

isotropy_type SymmetricOf(isotropy_type iso) {
  if (isCartesian(iso)) return SYMMETRIC;
  if (isEarth(iso)) return EARTH_SYMMETRIC;
  return isSpherical(iso) ? SPHERICAL_SYMMETRIC : ISO_MISMATCH;
}

isotropy_type IsotropicOf(isotropy_type iso) {
  if (isCartesian(iso)) return ISOTROPIC;
  if (isEarth(iso)) return EARTH_ISOTROPIC;
  return isSpherical(iso) ? SPHERICAL_ISOTROPIC : ISO_MISMATCH;
}
