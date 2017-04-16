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



#include <Rmath.h>
#include "RF.h"
#include "Operator.h"

/* simplifying functions
   turn vector of x into length ||x|| and similar
   reduces C(x,y) to C(x - y)

   needs preceeding analysis of submodels!
   (same for preferred methods; bottom-up-analysis needed)
*/

// keep always the orderung


//////////////////////////////////////////////////////////////////////
// Earth coordinate systems : stay within the system
//////////////////////////////////////////////////////////////////////

double mod(double x, double modulus) { return (x - FLOOR(x / modulus) * modulus); }
double isomod(double x, double modulus) {  
  double 
    twomodulus = 2.0 * modulus;
  return modulus - FABS(Mod(x, twomodulus) - modulus);
}
double lonmod(double x, double modulus) {  
  double 
    halfmodulus = 0.5 * modulus,
    y = x + halfmodulus;
  return Mod(y, modulus) - halfmodulus;
}

double latmod(double x, double modulus) {  
  double 
    halfmodulus = 0.5 * modulus,
    twomodulus = 2.0 * modulus,
    y = x - halfmodulus;
  return FABS(Mod(y, twomodulus) - modulus) - halfmodulus;
}


#define STATMOD(ZZ, X, lon, lat)				\
    ALLOC_NEW(Searth, X, dim+1, X);				\
    STATMODE_BASE(X, ZZ, lon, lat);				\
    for (d = 2; d < dim; d++) X[d] = ZZ[d]

void statmod2(double *x, double lon, double lat, double *y) {
  STATMODE_BASE(y, x, lon, lat);
}


#define piD180 (M_PI * 0.00555555555555555555555)
#define H80Dpi (180.0 * INVPI)
#define ESnonstat2iso(x, y)						\
  int d, dim=cov->xdimgatter;						\
  ALLOC_NEW(Searth, X, dim+1, X);					\
  double xx[2] = {x[0], x[1]},						\
    yy[2] = {y[0], y[1]},						\
    cosine = SIN(xx[1]) * SIN(yy[1]) +					\
             (COS(xx[0]) * COS(yy[0]) + SIN(xx[0]) * SIN(yy[0])) *      \
             COS(xx[1]) * COS(yy[1]);				       	\
    cosine = cosine <= 1.0 ? (cosine < -1.0 ? -1.0 : cosine) : 1.0;	\
    X[0] = ACOS(cosine);						\
    for (d = 2; d < dim; d++) X[d-1] = x[d] - y[d];	
			

#define ISOMOD(ZZ, X, maxangle)				\
  int d, dim=cov->xdimgatter;				\
  ALLOC_NEW(Searth, X, dim+1, X);				\
  X[0]=isomod(ZZ[0], maxangle);					\
  for (d = 1; d < dim; d++) X[d] = ZZ[d]


void EarthIso2EarthIso(double *x, cov_model *cov, double *v) {
  ISOMOD(x, X, 180);
  CovList[cov->nr].cov(X, cov, v); // nicht gatternr
}
void logEarthIso2EarthIso(double *x, cov_model *cov, double *v, double *Sign) {
  ISOMOD(x, X, 180);
  CovList[cov->nr].log(X, cov, v, Sign);// nicht gatternr
}
void NonstatEarth2EarthIso(double *x, double *y, cov_model *cov, double *v) {
  ESnonstat2iso(piD180 * x, piD180 * y);
  X[0] *= H80Dpi;
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logNonstatEarth2EarthIso(double *x, double *y, cov_model *cov, double *v, 
		     double *Sign) {
  ESnonstat2iso(piD180 * x, piD180 * y);
  X[0] *= H80Dpi;
  CovList[cov->nr].log(X, cov, v, Sign);// nicht gatternr
}

void Earth2Earth(double *x,  cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;				
  STATMOD(x, X, 360, 180);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logEarth2Earth(double *x, cov_model *cov, double *v, 
		    double *Sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, 360, 180);
  CovList[cov->nr].log(X, cov, v, Sign);// nicht gatternr
}
void NonstatEarth2Earth(double *x, double *y, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;					
  STATMOD(x, X, 360, 180);					
  STATMOD(y, Y, 360, 180);					
  CovList[cov->nr].nonstat_cov(X, Y, cov, v);// nicht gatternr	
}
void logNonstatEarth2Earth(double *x, double *y, cov_model *cov, double *v, 
		     double *Sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, 360, 180);
  STATMOD(y, Y, 360, 180);
  CovList[cov->nr].nonstatlog(X, Y, cov, v, Sign);// nicht gatternr
}

// BERRETH :: sphere2earth
void Sphere2Earth(double *x, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(H80Dpi * x, X, 360, 180);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}


void EarthIso2SphereIso(double *x, cov_model *cov, double *v) {
  ISOMOD(piD180 * x, X, M_PI);
  CovList[cov->nr].cov(X, cov, v); // nicht gatternr
}
void logEarthIso2SphereIso(double *x, cov_model *cov, double *v, double *Sign) {
  ISOMOD(piD180 * x, X, M_PI);
  CovList[cov->nr].log(X, cov, v, Sign);// nicht gatternr
}
void NonstatEarth2SphereIso(double *x, double *y, cov_model *cov, double *v) {
  ESnonstat2iso(piD180 * x, piD180 * y);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logNonstatEarth2SphereIso(double *x, double *y, cov_model *cov, double *v, 
		     double *Sign) {
  ESnonstat2iso(piD180 * x, piD180 * y);
  CovList[cov->nr].log(X, cov, v, Sign);// nicht gatternr
}

void Earth2Sphere(double *x, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;			       
  //  STATMOD(piD180 * x, X, M_2_PI, M_PI);
  
  //  PMI(cov->calling);

     ALLOC_NEW(Searth, X, dim+1, X);				
      STATMODE_BASE(X, piD180 * x, M_2_PI, M_PI);
      for (d = 2; d < dim; d++) X[d] = piD180 * x[d];

 
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logEarth2Sphere(double *x, cov_model *cov, double *v, double *Sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(piD180 * x, X, M_2_PI, M_PI);
  CovList[cov->nr].log(X, cov, v, Sign);// nicht gatternr
}
void NonstatEarth2Sphere(double *x, double *y, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;	
  STATMOD(piD180 * x, X, M_2_PI, M_PI);
  STATMOD(piD180 * y, Y, M_2_PI, M_PI);
  CovList[cov->nr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logNonstatEarth2Sphere(double *x, double *y, cov_model *cov, double *v, 
		     double *Sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(piD180 * x, X, M_2_PI, M_PI);
  STATMOD(piD180 * y, Y, M_2_PI, M_PI);
  CovList[cov->nr].nonstatlog(X, Y, cov, v, Sign);// nicht gatternr
}



void SphereIso2SphereIso(double *x, cov_model *cov, double *v) {
  ISOMOD(x, X, M_PI);
  CovList[cov->nr].cov(X, cov, v); // nicht gatternr
}
void logSphereIso2SphereIso(double *x, cov_model *cov, double *v, 
			    double *Sign) {
  ISOMOD(x, X, M_PI);
  CovList[cov->nr].log(X, cov, v, Sign);// nic+ht gatternr
}	
void NonstatSphere2SphereIso(double *x, double *y, cov_model *cov, double *v) {
  ESnonstat2iso(x, y);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logNonstatSphere2SphereIso(double *x, double *y, cov_model *cov,
				double *v, double *Sign) {
  ESnonstat2iso(x, y);
  CovList[cov->nr].log(X, cov, v, Sign);// nicht gatternr
}

void Sphere2Sphere(double *x, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, M_2_PI, M_PI);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logSphere2Sphere(double *x, cov_model *cov,
				double *v, double *Sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, M_2_PI, M_PI);
  CovList[cov->nr].log(X, cov, v, Sign);// nicht gatternr
}
void NonstatSphere2Sphere(double *x, double *y, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, M_2_PI, M_PI);
  STATMOD(y, Y, M_2_PI, M_PI); 
  CovList[cov->nr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logNonstatSphere2Sphere(double *x, double *y, cov_model *cov,
				double *v, double *Sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, M_2_PI, M_PI);
  STATMOD(y, Y, M_2_PI, M_PI);   
  CovList[cov->nr].nonstatlog(X, Y, cov, v, Sign);// nicht gatternr
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

#define earth2cartInner(RAEQU, RPOL)					\
  assert(cov->xdimprev >= 2 && cov->xdimgatter >=2);/*not necessarily equal*/ \
  int origdim=cov->xdimprev;						\
  bool Time =  Loc(cov)->Time,						\
    height = origdim > 2 + (int) Time;					\
  double Rcos, X[4], Y[4];						\
  EARTH_TRAFO(X, x, height ? ((RAEQU) + x[EARTH_HEIGHT]) : (RAEQU),	\
	      height ? ((RPOL) + x[EARTH_HEIGHT]) : (RPOL));		\
  EARTH_TRAFO(Y, y, height ? ((RAEQU) + y[EARTH_HEIGHT]) : (RAEQU),	\
	      height ? ((RPOL) + y[EARTH_HEIGHT]) : (RPOL));		\
  if (Time) X[3] = x[origdim - 1];


 // !!! if anything is changed here change also structtrafoproc !!!!
#define earth2cartInnerStatX(RAEQU, RPOL)				\
   assert(cov->xdimprev >= 2 && cov->xdimgatter >= 2);			\
  int origdim=cov->xdimprev;						\
  bool Time =  Loc(cov)->Time,						\
    height = origdim > 2 + (int) Time;					\
  double Rcos;							\
  EARTH_TRAFO(X, x, height ? ((RAEQU) + x[EARTH_HEIGHT]) : (RAEQU),	\
	      height ? ((RPOL) + x[EARTH_HEIGHT]) : (RPOL));		\
  if (Time) X[3] = x[origdim - 1];


#define earth2cartInnerStat(RAEQU, RPOL)				\
  double X[4];								\
  earth2cartInnerStatX(RAEQU, RPOL)

 
void Earth2Cart(double *x, cov_model *cov,double RAEQU, double RPOL, 
		double *X){
  // !!! if anything is changed here change also structtrafoproc !!!!
  earth2cartInnerStatX(RAEQU, RPOL)
}

void EarthKM2CartStat(double *x, cov_model *cov, double *v) {
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].cov(X, cov, v);// nicht gatternr
}
void logEarthKM2CartStat(double *x, cov_model *cov, double *v, double *Sign) {
   earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
   CovList[cov->secondarygatternr].log(X, cov, v, Sign);// nicht gatternr
}
void EarthKM2Cart(double *x, double *y, cov_model *cov, double *v) {
  //PMI(cov->calling);
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logEarthKM2Cart(double *x, double *y, cov_model *cov, double *v,
		     double *Sign) {
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].nonstatlog(X, Y, cov, v, Sign);// nicht gatternr
}

void EarthMiles2CartStat(double *x, cov_model *cov, double *v) {
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].cov(X, cov, v);// nicht gatternr
}
void logEarthMiles2CartStat(double *x, cov_model *cov, double *v, double *Sign) {
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].log(X, cov, v, Sign);// nicht gatternr
}
void EarthMiles2Cart(double *x, double *y, cov_model *cov, double *v) {
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logEarthMiles2Cart(double *x, double *y, cov_model *cov, double *v,
		     double *Sign) {
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].nonstatlog(X, Y, cov, v, Sign);// nicht gatternr
}

#define QP 0
#define QZenit 9
#define Qtotal 12
int checkEarth(cov_model *cov){
  // ACHTUNG! KEIN AUFRUF VON SUB[0] !
  if (cov->domprev == XONLY// 20.2.14: warum war es vorher cov->domown == XONLY?
      && isSymmetric(cov->isoprev)) {
    // darf nie auf "XONLY-stat" transformiert sein. Es kann aber zB bei 'shape'
    // sein, dass ein XONLY ankommt.
    return ERRORKERNEL;
  }
  
  COND_NEW_STORAGE(earth, X);

  if (cov->gatternr >= FIRST_PLANE && cov->gatternr <= LAST_PLANE) {
    assert(cov->tsdim>=2 && cov->xdimown == cov->tsdim); // oder 3, 4 !! regeln!
    if (!R_FINITE(GLOBAL.coords.zenit[0]) || !R_FINITE(GLOBAL.coords.zenit[1])){
      if (GLOBAL.internal.warn_zenit) {
	GLOBAL.internal.warn_zenit = false;
	char msg[255];
	SPRINTF(msg, "tried to use non-finite values of '%s' in a coordinate transformation\n", coords[ZENIT]);
	warning(msg);
      } 
      SERR1("'%s' not finite!", coords[ZENIT]);
    }
    
   if (cov->gatternr == EARTHKM2GNOMONIC || 
       cov->gatternr == EARTHMILES2GNOMONIC) {
     int d;
     double Rcos, X[4];		       
      if (cov->Searth == NULL) NEW_STORAGE(earth);
      if (cov->gatternr == EARTHKM2GNOMONIC) {
	EARTH_TRAFO(X, GLOBAL.coords.zenit, radiuskm_aequ, radiuskm_pol); 
      } else {
 	EARTH_TRAFO(X, GLOBAL.coords.zenit, radiusmiles_aequ, radiusmiles_pol);
      }
      double Rsq = 0.0,
	*zenit = cov->Searth->cart_zenit;
      for (d=0; d<=2; d++) Rsq += X[d] * X[d];
      for (d=0; d<=2; d++) {
	zenit[d] = X[d] / Rsq;
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
    assert(cov->gatternr >= EARTHKM2CART && cov->gatternr <= EARTHMILES2CART);
    
  }

  return NOERROR;
}




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Orthogonal
#define orthDefStat					\
  int d, k, m,						\
    dim=cov->xdimgatter;				\
  assert(dim >= 2);					\
  ALLOC_NEW(Searth, U, dim+1, U);			\
  double *P = cov->Searth->P
 

#define orthTrafoStat		\
  for (m=d=0; d<3; d++) {			\
    U[d] = 0.0;					\
    for(k=0; k<3; k++, m++) {			\
      U[d] += P[m] * X[k];	\
   }						\
  }						\
  if (U[2] < 0.0) ERR("location(s) not in direction of the zenit");\
  for (d = 2; d<dim; d++) U[d] = x[d]

#define orthDef orthDefStat;			\
  ALLOC_NEW(Searth, V, dim+1, V)         	

#define orthTrafo					\
  for (m=d=0; d<3; d++) {				\
    U[d] = V[d] = 0.0;					\
    for(k=0; k<3; k++, m++) {				\
      U[d] += P[m] * X[k];				\
      V[d] += P[m] * Y[k];				\
    }							\
  }							\
  if (U[2] < 0.0 || V[2] < 0.0)				\
    ERR("location(s) not in direction of the zenit");	\
  for (d = 2; d<dim; d++) {U[d]	= x[d]; V[d] = y[d];}

void EarthKM2OrthogStat(double *x, cov_model *cov, double *v) {
  assert(cov->gatternr == EARTHKM2ORTHOGRAPHIC);
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol); // results in X[0], X[1]
  orthDefStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].cov(U, cov, v);// nicht gatternr
}
void logEarthKM2OrthogStat(double *x, cov_model *cov, double *v,
			       double *Sign) {
  assert(cov->gatternr == EARTHKM2ORTHOGRAPHIC);
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  orthDefStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].log(U, cov, v, Sign);// nicht gatternr
}
void EarthKM2Orthog(double *x, double *y, cov_model *cov, double *v) {
  assert(cov->gatternr == EARTHKM2ORTHOGRAPHIC);
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthDef;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstat_cov(U, V, cov, v);// nicht gatternr
}
void logEarthKM2Orthog(double *x, double *y, cov_model *cov, double *v,
		     double *Sign) {
  assert(cov->gatternr == EARTHKM2ORTHOGRAPHIC);
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthDef;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstatlog(U, V, cov, v, Sign);// nicht gatternr
}



void EarthMiles2OrthogStat(double *x, cov_model *cov, double *v) {
  assert(cov->gatternr == EARTHMILES2ORTHOGRAPHIC);
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol); // results in X[0], X[1]
 orthDefStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].cov(U, cov, v);// nicht gatternr
}
void logEarthMiles2OrthogStat(double *x, cov_model *cov, double *v,
			       double *Sign) {
  assert(cov->gatternr == EARTHMILES2ORTHOGRAPHIC);
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  orthDefStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].log(U, cov, v, Sign);// nicht gatternr
}
void EarthMiles2Orthog(double *x, double *y, cov_model *cov, double *v) {
   assert(cov->gatternr == EARTHMILES2ORTHOGRAPHIC);
 earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  orthDef;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstat_cov(U, V, cov, v);// nicht gatternr
}
void logEarthMiles2Orthog(double *x, double *y, cov_model *cov, double *v,
		     double *Sign) {
  assert(cov->gatternr == EARTHMILES2ORTHOGRAPHIC);
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
 orthDef;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstatlog(U, V, cov, v, Sign);// nicht gatternr
}




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Gnomonic

#define GnomonicStat						\
  double factor = 0.0,						\
    *cart_zenit = cov->Searth->cart_zenit;				\
  for (d=0; d<3; d++) factor += cart_zenit[d] * X[d];		\
  if (factor <= 0.0)						\
    ERR1("locations not on the half-sphere given by the '%s'.", coords[ZENIT]);\
  for (d=0; d<3; d++) X[d] /= factor

#define Gnomonic							\
  double factor = 0.0,							\
    factorY = 0.0,							\
    *cart_zenit = cov->Searth->cart_zenit;				\
  for (d=0; d<3; d++) {							\
    factor += cart_zenit[d] * X[d];					\
    factorY += cart_zenit[d] * Y[d];					\
  }									\
  if (factor <= 0.0 || factorY <= 0.0) 					\
    ERR1("locations not on the half-sphere given by the '%s'.", coords[ZENIT]);\
  for (d=0; d<3; d++) {							\
    X[d] /= factor;							\
    Y[d] /= factorY;							\
  }

void Earth2GnomonicStat(double *x, cov_model *cov, double *v) {
  assert(cov->gatternr == EARTHKM2GNOMONIC || 
	 cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol); // results in X[0], X[1]
  orthDefStat;
  GnomonicStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].cov(U, cov, v);// nicht gatternr
}
void logEarth2GnomonicStat(double *x, cov_model *cov, double *v,
			       double *Sign) {
  assert(cov->gatternr == EARTHKM2GNOMONIC || 
	 cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  orthDefStat;
  GnomonicStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].log(U, cov, v, Sign);// nicht gatternr
}
void Earth2Gnomonic(double *x, double *y, cov_model *cov, double *v) {
  assert(cov->gatternr == EARTHKM2GNOMONIC || 
	 cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthDef;
  Gnomonic;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstat_cov(U, V, cov, v);// nicht gatternr
}
void logEarth2Gnomonic(double *x, double *y, cov_model *cov, double *v,
		     double *Sign) {
  assert(cov->gatternr == EARTHKM2GNOMONIC || 
	 cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthDef;
  Gnomonic;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstatlog(U, V, cov, v, Sign);// nicht gatternr
}





bool is_any(isofct iso, cov_fct *C) {
  int i;
  for (i=0; i < C->variants; i++) if (iso(C->Isotropy[i])) return true;
  return false;
}

bool is_all(isofct iso, cov_fct *C) {
  int i;
  for (i=0; i < C->variants; i++) if (!iso(C->Isotropy[i])) return false;
  return true;
}

bool is_any(typusfct t, cov_fct *C) {
  int i;
  for (i=0; i < C->variants; i++) if (t(C->Typi[i])) return true;
  return false;
}

bool is_all(typusfct t, cov_fct *C) {
  int i;
  for (i=0; i < C->variants; i++) if (!t(C->Typi[i])) return false;
  return true;
}


bool isIsotropic(isotropy_type iso) {
  return iso == ISOTROPIC;
}
bool isAnyIsotropic(isotropy_type iso) {
  return iso == ISOTROPIC || iso == EARTH_ISOTROPIC || 
    iso == SPHERICAL_ISOTROPIC;
}

bool isSpaceIsotropic(isotropy_type iso) {
  return iso <= SPACEISOTROPIC;
}

bool isZeroSpaceIsotropic(isotropy_type iso) {
  return iso <= ZEROSPACEISO;
}

bool isVectorIsotropic(isotropy_type iso) {
  return iso == VECTORISOTROPIC;
}

bool isSymmetric(isotropy_type iso) {
  return iso <= SYMMETRIC;
}

bool isCartesian(isotropy_type iso) {
  return iso <= LAST_CARTESIAN;
}

bool isSpherical(isotropy_type iso) {
  return iso >= SPHERICAL_ISOTROPIC && iso <= SPHERICAL_COORDS;
}

bool isEarth(isotropy_type iso) {
  return iso >= EARTH_ISOTROPIC && iso <= EARTH_COORDS ;
}
bool isAnySpherical(isotropy_type iso) {
  // printf("any %d < %d  <= %d\n",  LAST_CARTESIAN, iso, LAST_SPHERICAL);
  return iso > LAST_CARTESIAN && iso <= LAST_ANYSPHERICAL;
}
bool isAnySphericalIso(isotropy_type iso) {
  return iso == SPHERICAL_ISOTROPIC || iso == EARTH_ISOTROPIC;
}
bool isAnySphericalNotIso(isotropy_type iso) {
  return iso != SPHERICAL_ISOTROPIC && iso != EARTH_ISOTROPIC ;
}
bool isPrevModelI(cov_fct *C){
  bool is = C->Isotropy[0] == PREVMODELI;
  //  assert(!is || C->variants == 1 ||
  //	 (C->variants == 2 &&  C->Isotropy[0] == C->Isotropy[1]));
  return is;
}
bool isUnreduced(cov_fct *C){
  bool is = C->Isotropy[0] == UNREDUCED;
  assert(!is || C->variants == 1);
  return is;
}
bool isCoordinateSystem(isotropy_type iso){
  return iso == CARTESIAN_COORD || iso == SPHERICAL_COORDS || 
    iso == EARTH_COORDS;
}
bool atleastSpecialised(isotropy_type iso, isotropy_type as) {
  // printf("ATLEAST '%s' at least specialised as '%s'? \n", ISONAMES[iso], ISONAMES[as]);
  bool less = iso <= as;
  if (isCartesian(as)) return less; // also iso < 0 OK
  if (isSpherical(as)) return isSpherical(iso) && less; //
  if (isEarth(as)) {
     //    printf("YYY  %s  %d \n", ISONAMES[iso], isAnySpherical(iso));
    if (isEarth(iso)) return less;
    if (isSpherical(iso)) {
      //printf("iso=%d   %d %d %d  %d %d\n", iso, as, EARTH_COORDS, SPHERICAL_COORDS, as - EARTH_COORDS + SPHERICAL_COORDS, iso <= as - EARTH_COORDS + SPHERICAL_COORDS);
      return iso <= as - EARTH_COORDS + SPHERICAL_COORDS;
    }
    if (isCartesian(iso)) {
      return (as == EARTH_COORDS && iso == CARTESIAN_COORD) ||
	(as == EARTH_SYMMETRIC && iso == SYMMETRIC);
    }
    return false;
  }
  if (as == UNREDUCED) {
    return isCoordinateSystem(iso);
  }
  if (as == PREVMODELI) return true;

  //  PRINTF("'%s' not at least specialised as %s\n", ISONAMES[iso], ISONAMES[as]);
  // crash();
  BUG;
}

bool isCylinder(isotropy_type iso) {
  BUG;
  return iso == CYLINDER_COORD;
}

bool equal_coordinate_system(isotropy_type iso1, isotropy_type iso2) {
  return (isCartesian(iso1) && isCartesian(iso2))
    || (isAnySpherical(iso1) && isAnySpherical(iso2))
    // || (isCylinder(iso1) && isCylinder(iso2)) 
    || iso1==UNREDUCED;
}

bool equal_coordinate_system(isotropy_type iso1, isotropy_type iso2, 
			     bool refined) {
  if (!refined) return equal_coordinate_system(iso1, iso2);

  return (isCartesian(iso1) && isCartesian(iso2))
    || (isSpherical(iso1) && isSpherical(iso2))
    || (isEarth(iso1) && isEarth(iso2))
    // || (isCylinder(iso1) && isCylinder(iso2)) 
    || iso1==UNREDUCED;
}

isotropy_type UpgradeToCoordinateSystem(isotropy_type iso) {
  return iso == ZEROSPACEISO || iso == VECTORISOTROPIC|| iso == SYMMETRIC
    ? CARTESIAN_COORD 
    : iso == EARTH_SYMMETRIC ? EARTH_COORDS
    : iso == SPHERICAL_SYMMETRIC ? SPHERICAL_COORDS
    : isCoordinateSystem(iso) ? iso
    : ISO_MISMATCH;
}

isotropy_type CoordinateSystemOf(isotropy_type iso) {
  return isCartesian(iso) ? CARTESIAN_COORD 
    : isEarth(iso) ? EARTH_COORDS
    : isSpherical(iso) ? SPHERICAL_COORDS
    : ISO_MISMATCH;
}

isotropy_type IsotropicOf(isotropy_type iso) {
  return isCartesian(iso) ? ISOTROPIC
    : isEarth(iso) ? EARTH_ISOTROPIC
    : isSpherical(iso) ? SPHERICAL_ISOTROPIC
    : ISO_MISMATCH;
}


isotropy_type SymmetricOf(isotropy_type iso) {
  return isCartesian(iso) ? SYMMETRIC
    : isEarth(iso) ? EARTH_SYMMETRIC
    : isSpherical(iso) ? SPHERICAL_SYMMETRIC
    : ISO_MISMATCH;
}



int change_coord_system(isotropy_type callingisoown,
			     isotropy_type isoprev,
			     int tsdimprev, int xdimprev,
			     int *nr, isotropy_type  *newisoprev,
			     int *newtsdim, int *newxdim, bool time) {
  // isoprev is replaced by a coordinate transformation and a nwe 

  switch(callingisoown) {
  case EARTH_COORDS : case EARTH_SYMMETRIC:
    if (isCartesian(isoprev)) {
      if (xdimprev != tsdimprev) SERR("reduced coordinates not allowed");
      if (STRCMP(GLOBAL.coords.newunits[0], UNITS_NAMES[units_km]) == 0){
	*nr = isoprev == GNOMONIC_PROJ ? EARTHKM2GNOMONIC  
	  : isoprev == ORTHOGRAPHIC_PROJ ?  EARTHKM2ORTHOGRAPHIC
	  : EARTHKM2CART;
      } else if (STRCMP(GLOBAL.coords.newunits[0], 
			UNITS_NAMES[units_miles]) == 0) {
	*nr = isoprev == GNOMONIC_PROJ ? EARTHMILES2GNOMONIC  
	  : isoprev == ORTHOGRAPHIC_PROJ ?  EARTHMILES2ORTHOGRAPHIC
	  :EARTHMILES2CART;
      } else {
	SERR4("only units '%s' and '%s' are allowed. Got '%s' (user's '%s').",
	      UNITS_NAMES[units_km], UNITS_NAMES[units_miles], 
	      GLOBAL.coords.newunits[0], GLOBAL.coords.curunits[0]);
      }
      if (isoprev == GNOMONIC_PROJ || isoprev == ORTHOGRAPHIC_PROJ) {
	*newtsdim = tsdimprev;
	*newxdim = xdimprev;
	*newisoprev = isoprev;
      } else if (isCartesian(isoprev)) {
	*newisoprev = callingisoown == EARTH_COORDS ? CARTESIAN_COORD
	  :  callingisoown == EARTH_SYMMETRIC ? SYMMETRIC : ISO_MISMATCH;
	if (tsdimprev == 2) {
	  //assert(xdimprev == 2);
	  *newtsdim = *newxdim = 3;
	} else if (tsdimprev == 3 && time) {
	   *newtsdim = *newxdim = 4;
	} else { 
	  *newtsdim = *newxdim = xdimprev;	  
	}	
      } else {
	BUG;
      }
    } else {
      //if (!isSpherical(callingisoown) && callingisoown != EARTH_ISOTROPIC) {
      //	printf("%s -->  %s\n", ISONAMES[callingisoown], ISONAMES[isoprev]);
      //BUG;
      //      }
      return ERRORODDCOORDTRAFO; // NotProgrammedYet("");
    }      
    break;
  case EARTH_ISOTROPIC : case SPHERICAL_ISOTROPIC :
    //    if (!isSpherical(callingisoown) && callingisoown != EARTH_ISOTROPIC) {
    //	printf("\n\n%s -->  %s\n", ISONAMES[callingisoown], ISONAMES[isoprev]);
    //	BUG;
    //    }
    return ERRORWRONGISO;
  default:    
    // printf("\n\n\n\n\n\n%s -->  %s\n", ISONAMES[callingisoown], ISONAMES[isoprev]);
    return ERRORODDCOORDTRAFO; // NotProgrammedYet("");
  }

  //assert(*newtsdim == 4);
  return NOERROR;
}


int SetGatter(domain_type domprev, domain_type domnext, 
	      isotropy_type isoprev, isotropy_type isonext, 
	      //	      bool first, 
	      int *nr, int *delflag) {  
  if (domprev < domnext) 
    SERR2("Cannot call more complex models ('%s') from simpler ones ('%s')",
	  DOMAIN_NAMES[(int) domnext], DOMAIN_NAMES[(int) domprev]);
  if (isAnySpherical(isonext)) {
    if(!isAnySpherical(isoprev)) BUG; //to do?: back from cartesian not possible
    //                                   currently
    if (isoprev == UNREDUCED) {
      if (isoprev != isonext || domprev != domnext)
	SERR("unclear transformation of earth/spherical coordinates");
      *nr = isEarth(isoprev) ? E2E : Sph2Sph;
    }
    if (domnext == XONLY) {     
      if (isAnySphericalIso(isonext)) {
	if ( (isAnySphericalNotIso(isoprev) && domprev != KERNEL) || 
	     (isAnySphericalIso(isoprev) && domprev != XONLY))
	  SERR("impossible spherical trafo");
	*nr = isEarth(isoprev) 
	  ? (isSpherical(isonext) ? E2SphIso : E2EIso) 
	  : (isSpherical(isonext) ? Sph2SphIso : ISO_MISMATCH); 	
      } else { // next: xonly, not iso; e.g. for shape functions
	assert(isAnySphericalNotIso(isonext));
	if (domprev == KERNEL) SERR("earth trafo not possible"); // mathemati-
	// cally likely not correct
	assert(domprev == XONLY);
	if (!isAnySphericalNotIso(isoprev))  
	  SERR("mathematically not correct coordinate transformation"); // mathematically likely not correct. to do
	*nr = isEarth(isoprev) 
	  ? (isSpherical(isonext) ? E2Sph : E2E)
	  : (isSpherical(isonext) ? Sph2Sph : ISO_MISMATCH);	
      }
    } else { // domnext == KERNEL
      assert(isAnySphericalNotIso(isonext));
      assert(domprev == KERNEL);
      assert(isAnySphericalNotIso(isoprev));
      
      *nr = isEarth(isoprev) 
	? (isSpherical(isonext) ? E2Sph : E2E)
	: (isSpherical(isonext) ? Sph2Sph : ISO_MISMATCH);
     }

     return NOERROR;
  }


    
  if (isoprev == CYLINDER_COORD || isonext ==  CYLINDER_COORD)
    SERR("general spherical coordinates not programmed yet");

  //////////////
  // cartesian
  //////////////

  assert(isCartesian(isonext));
  bool isoOK = isoprev == isonext || 
    (isoprev > isonext && isonext <= CARTESIAN_COORD); // keine Trafo innerhalb
  // der projektionen -- bislang zumindest. ? to do
  if (!isoOK) 
    SERR2("cannot call more complex models ('%s') from simpler ones ('%s')",
	  ISONAMES[(int) domnext], ISONAMES[(int) domprev]);
  // if (isoprec == GNOMONIC_PROJ &&
     
  if (domnext == XONLY && 
      (domprev == KERNEL || 
       isoprev ==  CARTESIAN_COORD || 
       isoprev == SYMMETRIC ||
       isoprev == VECTORISOTROPIC || 
       isoprev == GNOMONIC_PROJ || 
       isoprev == ORTHOGRAPHIC_PROJ ||
       isoprev ==  ZEROSPACEISO)) {	
    switch (isonext) {
    case ISOTROPIC :
      *nr = S2ISO;
      break;
    case SPACEISOTROPIC :
      *nr = S2SP;
      break;
    case ZEROSPACEISO: case VECTORISOTROPIC: case SYMMETRIC: 
    case CARTESIAN_COORD: case GNOMONIC_PROJ: case ORTHOGRAPHIC_PROJ:
      *nr = S2S;
      // *delflag = DEL_COV - 5; ///
      break;
    case UNREDUCED:
      *nr = SId;
      break;
    default: BUG;
    }
  } else {    
    if (domprev == XONLY) {
      switch(isoprev) {
      case ISOTROPIC :
	*nr = ISO2ISO; 
	break;
      case SPACEISOTROPIC :
	*nr = (isonext == ISOTROPIC) ? SP2ISO : SP2SP;
	break;	 
      default: 
	PRINTF("SetGatter prev=%s; %s\n          next=%s %s\n",
	       DOMAIN_NAMES[domprev], ISONAMES[isoprev], 
	       DOMAIN_NAMES[domnext], ISONAMES[isonext]);
	BUG;
      }
    } else { // KERNEL domprev und domnext
      *nr = SId;
      *delflag = DEL_COV - 4;//
    }
  }

  return NOERROR;
}
