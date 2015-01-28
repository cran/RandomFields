/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2014-2015 Martin Schlather

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



#include <math.h>
 
#include "RF.h"
#include "Operator.h"

/* simplifying functions
   turn vector of x into length ||x|| and similar
   reduces C(x,y) to C(x - y)

   needs preceeding analysis of submodels!
   (same for preferred methods; bottom-up-analysis needed)
*/

// keep always the orderung
void iso2iso(double *x, cov_model *cov, double *v) {
  double y=fabs(*x);
  CovList[cov->nr].cov(&y, cov, v); // nicht gatternr
}
void logiso2iso(double *x, cov_model *cov, double *v, double *Sign) {
  double y=fabs(*x);
  CovList[cov->nr].log(&y, cov, v, Sign);// nicht gatternr
}
void spiso2spiso(double *x, cov_model *cov, double *v) {
  double y[2];
  y[0] = fabs(x[0]);
  y[1] = fabs(x[1]);
  CovList[cov->nr].cov(y, cov, v);// nicht gatternr
}
void logspiso2spiso(double *x, cov_model *cov, double *v, double *Sign) {
  double y[2];
  y[0] = fabs(x[0]);
  y[1] = fabs(x[1]);
  CovList[cov->nr].log(y, cov, v, Sign);// nicht gatternr
}
void spacetime2iso(double *x, cov_model *cov, double *v) {
  double y=sqrt(x[0] * x[0] + x[1] * x[1]);
  CovList[cov->nr].cov(&y, cov, v);// nicht gatternr
}
void logspacetime2iso(double *x, cov_model *cov, double *v, double *Sign) {
  double y=sqrt(x[0] * x[0] + x[1] * x[1]);
  CovList[cov->nr].log(&y, cov, v, Sign);// nicht gatternr
}

void Stat2iso(double *x, cov_model *cov, double *v) {
  double b = 0.0;
  int i,
    dim=cov->xdimgatter;  
  for (i=0; i<dim; i++) {
    //
    b += x[i] * x[i];
    //  
  }
  b = sqrt(b);
  CovList[cov->nr].cov(&b, cov, v);// nicht gatternr
}
void logStat2iso(double *x, cov_model *cov, double *v, double *Sign) {
  double b = 0.0;
  int i,
    dim=cov->xdimgatter;  
  for (i=0; i<dim; i++) {
    b += x[i] *x[i];
  }
  b = sqrt(b);
  CovList[cov->nr].log(&b, cov, v, Sign);// nicht gatternr
}
void Nonstat2iso(double *x, double *y, cov_model *cov, double *v) {
  double a, b;
  int d,
    dim=cov->xdimgatter;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d];
    b += a * a;
  }
  b = sqrt(b);
  CovList[cov->nr].cov(&b, cov, v);// nicht gatternr
}
void logNonstat2iso(double *x, double *y, cov_model *cov, double *v,
		    double *Sign) {
  double a, b;
  int d, 
    dim=cov->xdimgatter;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d];
    b += a * a;
  }
  b = sqrt(b);
  CovList[cov->nr].log(&b, cov, v, Sign);// nicht gatternr
}
void Stat2spacetime(double *x, cov_model *cov, double *v) {
  double b, z[2];
  int i,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, i=0; i<dim; i++)  b += x[i] * x[i];
  z[0] = sqrt(b);
  z[1] = fabs(x[dim]);
  CovList[cov->nr].cov(z, cov, v);// nicht gatternr
}
void logStat2spacetime(double *x, cov_model *cov, double *v, double *Sign) {
  double b, z[2];
  int i,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, i=0; i<dim; i++)  b += x[i] * x[i];
  z[0] = sqrt(b);
  z[1] = fabs(x[dim]);
  CovList[cov->nr].log(z, cov, v, Sign);// nicht gatternr
}
void Nonstat2spacetime(double *x, double *y, cov_model *cov, double *v) {
  double a, b, z[2];
  int d,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d];
    b += a * a;
  }
  z[0] = sqrt(b);
  z[1] = fabs(x[dim] - y[dim]);
  CovList[cov->nr].cov(z, cov, v);// nicht gatternr
}
void logNonstat2spacetime(double *x, double *y, cov_model *cov, double *v,
			  double *Sign) {
  double a, b, z[2];
  int d,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d]; 
    b += a * a;
  }
  z[0] = sqrt(b);
  z[1] = fabs(x[dim] - y[dim]);
  CovList[cov->nr].log(z, cov, v, Sign);// nicht gatternr
}
void Stat2Stat(double *x, cov_model *cov, double *v) {
  CovList[cov->nr].cov(x, cov, v);// nicht gatternr
}
void logStat2Stat(double *x, cov_model *cov, double *v, double *Sign) {
  CovList[cov->nr].log(x, cov, v, Sign);// nicht gatternr
}

#define nonstat2statInner						\
  int d,								\
    dim=cov->xdimgatter;						\
  ALLOC_NEW(Sgatter, z, dim, z);					\
  for (d=0; d<dim; d++) z[d] = x[d] - y[d]
 
void Nonstat2Stat(double *x, double *y, cov_model *cov, double *v) {
  nonstat2statInner;
  CovList[cov->nr].cov(z, cov, v);// nicht gatternr
}
void logNonstat2Stat(double *x, double *y, cov_model *cov, double *v, 
		     double *Sign) {
  nonstat2statInner;
  CovList[cov->nr].log(z, cov, v, Sign);// nicht gatternr
}

void Nonstat2Nonstat(double *x, double *y, cov_model *cov, double *v) {
  CovList[cov->nr].nonstat_cov(x, y, cov, v);// nicht gatternr
}

void logNonstat2Nonstat(double *x, double *y, cov_model *cov, double *v, 
			double *Sign) {
  CovList[cov->nr].nonstatlog(x, y, cov, v, Sign);// nicht gatternr
}


//////////////////////////////////////////////////////////////////////
// Earth coordinate systems : stay within the system
//////////////////////////////////////////////////////////////////////

double mod(double x, double modulus) { return (x - floor(x / modulus) * modulus); }
double isomod(double x, double modulus) {  
  double 
    twomodulus = 2.0 * modulus;
  return modulus - fabs(Mod(x, twomodulus) - modulus);
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
  return fabs(Mod(y, twomodulus) - modulus) - halfmodulus;
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
#define ESnonstat2iso(x, y)					\
  int d, dim=cov->xdimgatter;						\
  ALLOC_NEW(Searth, X, dim+1, X);					\
  double xx[2] = {x[0], x[1]},			\
    yy[2] = {y[0], y[1]};				\
    X[0] = acos(sin(xx[1]) * sin(yy[1]) +				\
	      (cos(xx[0]) * cos(yy[0]) + sin(xx[0]) * sin(yy[0])) *	\
	      cos(xx[1]) * cos(yy[1]) );				\
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

// BERRETH ::

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
  STATMOD(piD180 * x, X, M_2_PI, M_PI);
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
#define pi180 0.017453292519943295474
#define EARTH_TRAFO(X, ZZ, raequ, rpol)			\
  Rcos = raequ * cos(ZZ[EARTH_LATITUDE] * pi180);	\
  X[0] = Rcos * cos(ZZ[EARTH_LONGITUDE] * pi180);	\
  X[1] = Rcos * sin(ZZ[EARTH_LONGITUDE] * pi180);	\
  X[2] = rpol * sin(ZZ[EARTH_LATITUDE] * pi180)		

#define earth2cartInner(RAEQU, RPOL)					\
 assert(cov->xdimprev >= 2 && cov->xdimgatter >= 2);			\
  double Rcos, X[3], Y[3];			\
  EARTH_TRAFO(X, x, RAEQU, RPOL);			\
  EARTH_TRAFO(Y, y, RAEQU, RPOL)

#define earth2cartInnerStat(RAEQU, RPOL)				\
  assert(cov->xdimprev >= 2 && cov->xdimgatter >= 2);			\
  double Rcos, X[3];							\
  EARTH_TRAFO(X, x, RAEQU, RPOL)

#define radiuskm_aequ 6378.1
#define radiuskm_pol 6356.8
#define radiusmiles_aequ 3963.17
#define radiusmiles_pol 3949.93

void EarthKM2CartStat(double *x, cov_model *cov, double *v) {
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].cov(X, cov, v);// nicht gatternr
}
void logEarthKM2CartStat(double *x, cov_model *cov, double *v, double *Sign) {
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].log(X, cov, v, Sign);// nicht gatternr
}
void EarthKM2Cart(double *x, double *y, cov_model *cov, double *v) {
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
    assert(cov->tsdim>=2 && cov->xdimown == cov->tsdim); // oder 3 !! regeln!
    if (!R_FINITE(GLOBAL.coords.zenit[0]) || !R_FINITE(GLOBAL.coords.zenit[1])){
      if (GLOBAL.internal.warn_zenit) {
	GLOBAL.internal.warn_zenit = false;
	char msg[255];
	sprintf(msg, "tried to use non-finite values of '%s' in a coordinate transformation\n", coords[ZENIT]);
	warning(msg);
      } 
      SERR1("'%s' not finite!", coords[ZENIT]);
    }
    
   if (cov->gatternr == EARTHKM2GNOMONIC || 
       cov->gatternr == EARTHMILES2GNOMONIC) {
     int d;
     double Rcos, X[3];		       
      if (cov->Searth == NULL) NEW_STORAGE(earth);
      if (EARTHKM2GNOMONIC) {
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
      sin0 = sin(Zenit[0]),
	sin1 = sin(Zenit[1]),
	cos0 = cos(Zenit[0]),
	cos1 = cos(Zenit[1]),
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

void EarthKM2GnomonicStat(double *x, cov_model *cov, double *v) {
  assert(cov->gatternr == EARTHKM2GNOMONIC);
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol); // results in X[0], X[1]
  orthDefStat;
  GnomonicStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].cov(U, cov, v);// nicht gatternr
}
void logEarthKM2GnomonicStat(double *x, cov_model *cov, double *v,
			       double *Sign) {
    assert(cov->gatternr == EARTHKM2GNOMONIC);
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
   orthDefStat;
 GnomonicStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].log(U, cov, v, Sign);// nicht gatternr
}
void EarthKM2Gnomonic(double *x, double *y, cov_model *cov, double *v) {
    assert(cov->gatternr == EARTHKM2GNOMONIC);
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthDef;
  Gnomonic;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstat_cov(U, V, cov, v);// nicht gatternr
}
void logEarthKM2Gnomonic(double *x, double *y, cov_model *cov, double *v,
		     double *Sign) {
    assert(cov->gatternr == EARTHKM2GNOMONIC);
 earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthDef;
  Gnomonic;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstatlog(U, V, cov, v, Sign);// nicht gatternr
}


void EarthMiles2GnomonicStat(double *x, cov_model *cov, double *v) {
   assert(cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);//results in X[0], X[1]
  orthDefStat;
  GnomonicStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].cov(U, cov, v);// nicht gatternr
}
void logEarthMiles2GnomonicStat(double *x, cov_model *cov, double *v,
			       double *Sign) {
   assert(cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  orthDefStat;
  GnomonicStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].log(U, cov, v, Sign);// nicht gatternr .
}
void EarthMiles2Gnomonic(double *x, double *y, cov_model *cov, double *v) {
   assert(cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  orthDef;
  Gnomonic;
  orthTrafo;
 CovList[cov->secondarygatternr].nonstat_cov(U, V, cov, v);// nicht gatternr
}
void logEarthMiles2Gnomonic(double *x, double *y, cov_model *cov, double *v,
		     double *Sign) {
   assert(cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
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
  //  return iso <= CARTESIAN_COORD;
  return iso <= LAST_CARTESIAN;
}

bool isSpherical(isotropy_type iso) {
  return iso >= SPHERICAL_ISOTROPIC && iso <= SPHERICAL_COORD;
}

bool isEarth(isotropy_type iso) {
  return iso >= EARTH_ISOTROPIC && iso <= EARTH_COORD ;
}
bool isAnySpherical(isotropy_type iso) {
  // printf("any %d < %d  <= %d\n",  LAST_CARTESIAN, iso, LAST_SPHERICAL);
  return iso > LAST_CARTESIAN && iso <= LAST_SPHERICAL;
}
bool isAnySphericalIso(isotropy_type iso) {
  return iso == SPHERICAL_ISOTROPIC || iso == EARTH_ISOTROPIC;
}
bool isPrevModelI(cov_fct *C){
  bool is = C->Isotropy[0] == PREVMODELI;
  assert(!is || C->variants == 1 ||
	 (C->variants == 2 &&  C->Isotropy[0] == C->Isotropy[1]));
  return is;
}
bool isUnreduced(cov_fct *C){
  bool is = C->Isotropy[0] == UNREDUCED;
  assert(!is || C->variants == 1);
  return is;
}
bool isCoordinateSystem(isotropy_type iso){
  return iso == CARTESIAN_COORD || iso == SPHERICAL_COORD || iso == EARTH_COORD;
}
bool atleastSpecialised(isotropy_type iso, isotropy_type as) {
  // printf("ATLEAST '%s' at least specialised as '%s'? \n", ISONAMES[iso], ISONAMES[as]);
  if (isCartesian(as)) return iso <= as; // also iso < 0 OK
  if (isSpherical(as)) {
    if (as == SPHERICAL_COORD) return isSpherical(iso); // isCartesian(iso) ||
    else if (as == SPHERICAL_ISOTROPIC) 
      return iso == SPHERICAL_ISOTROPIC; // || iso == ISOTROPIC;
    else {BUG}    
  }
  if (isEarth(as)) {
    //    printf("YYY as=%d %s  %d \n", as == EARTH_COORD, ISONAMES[iso], isAnySpherical(iso));
    if (as == EARTH_COORD) return isAnySpherical(iso); //isCartesian(iso) || 
    else if (as == EARTH_ISOTROPIC) {
      return iso == SPHERICAL_ISOTROPIC || iso == EARTH_ISOTROPIC ;
    }
    else {BUG}
  }
  if (as == UNREDUCED) {
    return isCoordinateSystem(iso);
  }
  if (as == PREVMODELI) return true;

  //  PRINTF("subsequent coordinate system / isotropy:  %s\n", ISONAMES[as]);
  //  crash();
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

isotropy_type UpgradeToCoordinateSystem( isotropy_type iso) {
  return iso == SYMMETRIC
    ? CARTESIAN_COORD 
    : isCoordinateSystem(iso) ? iso
    : ISO_MISMATCH;
}
