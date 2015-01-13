/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2014 Martin Schlather

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
#include "Covariance.h"

/* simplifying functions
   turn vector of x into length ||x|| and similar
   reduces C(x,y) to C(x - y)

   needs preceeding analysis of submodels!
   (same for preferred methods; bottom-up-analysis needed)
*/

// #
// keep always the orderung
void iso2iso(double *x, cov_model *cov, double *v) {
  double y=fabs(*x);
  CovList[cov->nr].cov(&y, cov, v); // nicht gatternr
}
void logiso2iso(double *x, cov_model *cov, double *v, double *sign) {
  double y=fabs(*x);
  CovList[cov->nr].log(&y, cov, v, sign);// nicht gatternr
}
void spiso2spiso(double *x, cov_model *cov, double *v) {
  double y[2];
  y[0] = fabs(x[0]);
  y[1] = fabs(x[1]);
  CovList[cov->nr].cov(y, cov, v);// nicht gatternr
}
void logspiso2spiso(double *x, cov_model *cov, double *v, double *sign) {
  double y[2];
  y[0] = fabs(x[0]);
  y[1] = fabs(x[1]);
  CovList[cov->nr].log(y, cov, v, sign);// nicht gatternr
}
void spacetime2iso(double *x, cov_model *cov, double *v) {
  double y=sqrt(x[0] * x[0] + x[1] * x[1]);
  CovList[cov->nr].cov(&y, cov, v);// nicht gatternr
}
void logspacetime2iso(double *x, cov_model *cov, double *v, double *sign) {
  double y=sqrt(x[0] * x[0] + x[1] * x[1]);
  CovList[cov->nr].log(&y, cov, v, sign);// nicht gatternr
}

void Stat2iso(double *x, cov_model *cov, double *v) {
  double b = 0.0;
  int i,
    dim=cov->xdimgatter;  
  
  //   print("%d\n", dim); assert(false);
  //  print("dim=%d, %f %f, %d %ld\n", dim, x[0], x[1], cov->nr, cov->calling);
  //   print("%s\n", NAME(cov));

  //     PMI(cov->calling);
  // 
 
  for (i=0; i<dim; i++) {
    //
    b += x[i] * x[i];
    //  
  }
  b = sqrt(b);
  //  APMI(cov->calling->calling);  crash(cov);
  // 
  CovList[cov->nr].cov(&b, cov, v);// nicht gatternr
  //
  //printf("aaxxx\n");
  // print("%f\n", *v);
}
void logStat2iso(double *x, cov_model *cov, double *v, double *sign) {
  double b = 0.0;
  int i,
    dim=cov->xdimgatter;  
  for (i=0; i<dim; i++) {
    b += x[i] *x[i];
  }
  b = sqrt(b);
  CovList[cov->nr].log(&b, cov, v, sign);// nicht gatternr
  // print("%f\n", *v);
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
  //if (dim != 2) APMI(cov);
  //assert(dim == 2);
  CovList[cov->nr].cov(&b, cov, v);// nicht gatternr
  //printf("b=%f %f %s\n", b, *v, NAME(cov));
}
void logNonstat2iso(double *x, double *y, cov_model *cov, double *v,
		    double *sign) {
  double a, b;
  int d, 
    dim=cov->xdimgatter;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d];
    b += a * a;
  }
  b = sqrt(b);
  CovList[cov->nr].log(&b, cov, v, sign);// nicht gatternr
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
void logStat2spacetime(double *x, cov_model *cov, double *v, double *sign) {
  double b, z[2];
  int i,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, i=0; i<dim; i++)  b += x[i] * x[i];
  z[0] = sqrt(b);
  z[1] = fabs(x[dim]);
  CovList[cov->nr].log(z, cov, v, sign);// nicht gatternr
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
			  double *sign) {
  double a, b, z[2];
  int d,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d]; 
    b += a * a;
  }
  z[0] = sqrt(b);
  z[1] = fabs(x[dim] - y[dim]);
  CovList[cov->nr].log(z, cov, v, sign);// nicht gatternr
}
void Stat2Stat(double *x, cov_model *cov, double *v) {
  CovList[cov->nr].cov(x, cov, v);// nicht gatternr
}
void logStat2Stat(double *x, cov_model *cov, double *v, double *sign) {
  CovList[cov->nr].log(x, cov, v, sign);// nicht gatternr
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
		     double *sign) {
  nonstat2statInner;
  CovList[cov->nr].log(z, cov, v, sign);// nicht gatternr
}

void Nonstat2Nonstat(double *x, double *y, cov_model *cov, double *v) {
  CovList[cov->nr].nonstat_cov(x, y, cov, v);// nicht gatternr
}

void logNonstat2Nonstat(double *x, double *y, cov_model *cov, double *v, 
			double *sign) {
  CovList[cov->nr].nonstatlog(x, y, cov, v, sign);// nicht gatternr
}


//////////////////////////////////////////////////////////////////////
// Earth coordinate systems : stay within the system
//////////////////////////////////////////////////////////////////////

#define Mod(x, modulus) (x - floor(x / modulus) * modulus)
double mod(double x, double modulus) { return Mod(x, modulus); }
double negmod(double x, double modulus) {  
  double 
    halfmodulus = 0.5 * modulus,
    twomodulus = 2.0 * modulus,
    y = x - halfmodulus;
  return fabs(Mod(y, twomodulus) - modulus) - halfmodulus;
}
double posmod(double x, double modulus) {  
  double 
    twomodulus = 2.0 * modulus;
  return modulus - fabs(Mod(x, twomodulus) - modulus);
}

#define piD180 M_PI * 0.00555555555555555555555
#define H80Dpi 180.0 * INVPI
#define ESnonstat2iso(x, y)					\
  int d, dim=cov->xdimgatter;						\
  ALLOC_NEW(Searth, X, dim+1, X);					\
  double xx[2] = {x[0], x[1]},			\
    yy[2] = {y[0], y[1]};				\
    X[0] = acos(sin(xx[1]) * sin(yy[1]) +				\
	      (cos(xx[0]) * cos(yy[0]) + sin(xx[0]) * sin(yy[0])) *	\
	      cos(xx[1]) * cos(yy[1]) );				\
    for (d = 2; d < dim; d++) X[d-1] = x[d] - y[d];	
			
#define STATMOD(x, X, lon, lat)				\
  ALLOC_NEW(Searth, X, dim+1, X);				\
  X[0]=Mod(x[0], lon);					\
  X[1]=negmod(x[1], lat);				\
  for (d = 2; d < dim; d++) X[d] = x[d]

#define ISOMOD(x, X, maxangle)				\
  int d, dim=cov->xdimgatter;				\
  ALLOC_NEW(Searth, X, dim+1, X);				\
  X[0]=posmod(x[0], maxangle);					\
  for (d = 1; d < dim; d++) X[d] = x[d]


void EarthIso2EarthIso(double *x, cov_model *cov, double *v) {
  ISOMOD(x, X, 180);
  CovList[cov->nr].cov(X, cov, v); // nicht gatternr
}
void logEarthIso2EarthIso(double *x, cov_model *cov, double *v, double *sign) {
  ISOMOD(x, X, 180);
  CovList[cov->nr].log(X, cov, v, sign);// nicht gatternr
}
void NonstatEarth2EarthIso(double *x, double *y, cov_model *cov, double *v) {
  ESnonstat2iso(piD180 * x, piD180 * y);
  X[0] *= H80Dpi;
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logNonstatEarth2EarthIso(double *x, double *y, cov_model *cov, double *v, 
		     double *sign) {
  ESnonstat2iso(piD180 * x, piD180 * y);
  X[0] *= H80Dpi;
  CovList[cov->nr].log(X, cov, v, sign);// nicht gatternr
}

void Earth2Earth(double *x,  cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;				
  STATMOD(x, X, 360, 180);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logEarth2Earth(double *x, cov_model *cov, double *v, 
		    double *sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, 360, 180);
  CovList[cov->nr].log(X, cov, v, sign);// nicht gatternr
}
void NonstatEarth2Earth(double *x, double *y, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;					
  STATMOD(x, X, 360, 180);					
  STATMOD(y, Y, 360, 180);					
  CovList[cov->nr].nonstat_cov(X, Y, cov, v);// nicht gatternr	
}
void logNonstatEarth2Earth(double *x, double *y, cov_model *cov, double *v, 
		     double *sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, 360, 180);
  STATMOD(y, Y, 360, 180);
  CovList[cov->nr].nonstatlog(X, Y, cov, v, sign);// nicht gatternr
}

// BERRETH ::

void EarthIso2SphereIso(double *x, cov_model *cov, double *v) {
  //printf("iso angle x=%f\n", *x);
  ISOMOD(piD180 * x, X, M_PI);
  CovList[cov->nr].cov(X, cov, v); // nicht gatternr
}
void logEarthIso2SphereIso(double *x, cov_model *cov, double *v, double *sign) {
  ISOMOD(piD180 * x, X, M_PI);
  CovList[cov->nr].log(X, cov, v, sign);// nicht gatternr
}
void NonstatEarth2SphereIso(double *x, double *y, cov_model *cov, double *v) {
  ESnonstat2iso(piD180 * x, piD180 * y);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logNonstatEarth2SphereIso(double *x, double *y, cov_model *cov, double *v, 
		     double *sign) {
  ESnonstat2iso(piD180 * x, piD180 * y);
  CovList[cov->nr].log(X, cov, v, sign);// nicht gatternr
}

void Earth2Sphere(double *x, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(piD180 * x, X, M_2_PI, M_PI);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logEarth2Sphere(double *x, cov_model *cov, double *v, double *sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(piD180 * x, X, M_2_PI, M_PI);
  CovList[cov->nr].log(X, cov, v, sign);// nicht gatternr
}
void NonstatEarth2Sphere(double *x, double *y, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(piD180 * x, X, M_2_PI, M_PI);
  STATMOD(piD180 * y, Y, M_2_PI, M_PI);
  CovList[cov->nr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logNonstatEarth2Sphere(double *x, double *y, cov_model *cov, double *v, 
		     double *sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(piD180 * x, X, M_2_PI, M_PI);
  STATMOD(piD180 * y, Y, M_2_PI, M_PI);
  CovList[cov->nr].nonstatlog(X, Y, cov, v, sign);// nicht gatternr
}



void SphereIso2SphereIso(double *x, cov_model *cov, double *v) {
  ISOMOD(x, X, M_PI);
  CovList[cov->nr].cov(X, cov, v); // nicht gatternr
}
void logSphereIso2SphereIso(double *x, cov_model *cov, double *v, 
			    double *sign) {
  ISOMOD(x, X, M_PI);
  CovList[cov->nr].log(X, cov, v, sign);// nic+ht gatternr
}	
void NonstatSphere2SphereIso(double *x, double *y, cov_model *cov, double *v) {
  ESnonstat2iso(x, y);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logNonstatSphere2SphereIso(double *x, double *y, cov_model *cov,
				double *v, double *sign) {
  ESnonstat2iso(x, y);
  CovList[cov->nr].log(X, cov, v, sign);// nicht gatternr
}

void Sphere2Sphere(double *x, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, M_2_PI, M_PI);
  CovList[cov->nr].cov(X, cov, v);// nicht gatternr
}
void logSphere2Sphere(double *x, cov_model *cov,
				double *v, double *sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, M_2_PI, M_PI);
  CovList[cov->nr].log(X, cov, v, sign);// nicht gatternr
}
void NonstatSphere2Sphere(double *x, double *y, cov_model *cov, double *v) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, M_2_PI, M_PI);
  STATMOD(y, Y, M_2_PI, M_PI); 
  CovList[cov->nr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logNonstatSphere2Sphere(double *x, double *y, cov_model *cov,
				double *v, double *sign) {
  int d, dim=cov->xdimgatter;			       
  STATMOD(x, X, M_2_PI, M_PI);
  STATMOD(y, Y, M_2_PI, M_PI);   
  CovList[cov->nr].nonstatlog(X, Y, cov, v, sign);// nicht gatternr
}






//////////////////////////////////////////////////////////////////////
// Earth coordinate systems : change to cartesian system
//////////////////////////////////////////////////////////////////////



#define EARTH_LONGITUDE 0
#define EARTH_LATITUDE 1
#define pi180 0.017453292519943295474
#define EARTH_TRAFO(X, x, raequ, rpol)			\
  Rcos = raequ * cos(x[EARTH_LATITUDE] * pi180);	\
  X[0] = Rcos * cos(x[EARTH_LONGITUDE] * pi180);	\
  X[1] = Rcos * sin(x[EARTH_LONGITUDE] * pi180);	\
  X[2] = rpol * sin(x[EARTH_LATITUDE] * pi180)		

//   ; printf("x=%f =%f %f %f\n", x[0], x[1], raequ, rpol)

#define earth2cartInner(raequ, rpol)					\
 assert(cov->xdimprev >= 2 && cov->xdimgatter >= 2);			\
  double Rcos, X[3], Y[3];			\
  EARTH_TRAFO(X, x, raequ, rpol);			\
  EARTH_TRAFO(Y, y, raequ, rpol)

#define earth2cartInnerStat(raequ, rpol)				\
  assert(cov->xdimprev >= 2 && cov->xdimgatter >= 2);			\
  double Rcos, X[3];							\
  EARTH_TRAFO(X, x, raequ, rpol)

#define radiuskm_aequ 6378.1
#define radiuskm_pol 6356.8
#define radiusmiles_aequ 3963.17
#define radiusmiles_pol 3949.93

void EarthKM2CartStat(double *x, cov_model *cov, double *v) {
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].cov(X, cov, v);// nicht gatternr
}
void logEarthKM2CartStat(double *x, cov_model *cov, double *v, double *sign) {
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].log(X, cov, v, sign);// nicht gatternr
}
void EarthKM2Cart(double *x, double *y, cov_model *cov, double *v) {
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  //  printf("earth : %4.4f %4.4f y=%4.4f %4.4f X=%4.4f %4.4f %4.4f  Y=%4.4f %4.4f %4.4f\n", x[0], x[1], y[0], y[1], X[0], X[1], X[2], Y[0], Y[1], Y[2]);
  CovList[cov->secondarygatternr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logEarthKM2Cart(double *x, double *y, cov_model *cov, double *v,
		     double *sign) {
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].nonstatlog(X, Y, cov, v, sign);// nicht gatternr
}

void EarthMiles2CartStat(double *x, cov_model *cov, double *v) {
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  //printf("KM x=%f =%f %f %f\n", x[0], x[1], radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].cov(X, cov, v);// nicht gatternr
}
void logEarthMiles2CartStat(double *x, cov_model *cov, double *v, double *sign) {
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].log(X, cov, v, sign);// nicht gatternr
}
void EarthMiles2Cart(double *x, double *y, cov_model *cov, double *v) {
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logEarthMiles2Cart(double *x, double *y, cov_model *cov, double *v,
		     double *sign) {
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].nonstatlog(X, Y, cov, v, sign);// nicht gatternr
}

#define QP 0
#define QZenit 9
#define Qtotal 12
int checkEarth(cov_model *cov){
  // ACHTUNG! KEIN AUFRUF VON SUB[0] !

  //  print("earth:\n"); printf("%s\n", CovList[cov->gatternr].name);

  if (cov->domprev == XONLY// 20.2.14: warum war es vorher cov->domown == XONLY?
      && isSymmetric(cov->isoprev)) {
    // darf nie auf "XONLY-stat" transformiert sein. Es kann aber zB bei 'shape'
    // sein, dass ein XONLY ankommt.
    return ERRORKERNEL;
  }
  
  COND_NEW_STORAGE(earth, X);
  // 
  // PMI(cov);

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
 

#define orthTrafoStat				\
  for (m=d=0; d<3; d++) {			\
    U[d] = 0.0;					\
    for(k=0; k<3; k++, m++) {			\
       U[d] += P[m] * X[k];			\
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
			       double *sign) {
  assert(cov->gatternr == EARTHKM2ORTHOGRAPHIC);
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  orthDefStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].log(U, cov, v, sign);// nicht gatternr
}
void EarthKM2Orthog(double *x, double *y, cov_model *cov, double *v) {
  assert(cov->gatternr == EARTHKM2ORTHOGRAPHIC);
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthDef;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstat_cov(U, V, cov, v);// nicht gatternr
}
void logEarthKM2Orthog(double *x, double *y, cov_model *cov, double *v,
		     double *sign) {
  assert(cov->gatternr == EARTHKM2ORTHOGRAPHIC);
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthDef;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstatlog(U, V, cov, v, sign);// nicht gatternr
}



void EarthMiles2OrthogStat(double *x, cov_model *cov, double *v) {
  assert(cov->gatternr == EARTHMILES2ORTHOGRAPHIC);
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol); // results in X[0], X[1]
 orthDefStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].cov(U, cov, v);// nicht gatternr
}
void logEarthMiles2OrthogStat(double *x, cov_model *cov, double *v,
			       double *sign) {
  assert(cov->gatternr == EARTHMILES2ORTHOGRAPHIC);
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  orthDefStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].log(U, cov, v, sign);// nicht gatternr
}
void EarthMiles2Orthog(double *x, double *y, cov_model *cov, double *v) {
   assert(cov->gatternr == EARTHMILES2ORTHOGRAPHIC);
 earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  orthDef;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstat_cov(U, V, cov, v);// nicht gatternr
}
void logEarthMiles2Orthog(double *x, double *y, cov_model *cov, double *v,
		     double *sign) {
  assert(cov->gatternr == EARTHMILES2ORTHOGRAPHIC);
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
 orthDef;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstatlog(U, V, cov, v, sign);// nicht gatternr
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
			       double *sign) {
    assert(cov->gatternr == EARTHKM2GNOMONIC);
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
   orthDefStat;
 GnomonicStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].log(U, cov, v, sign);// nicht gatternr
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
		     double *sign) {
    assert(cov->gatternr == EARTHKM2GNOMONIC);
 earth2cartInner(radiuskm_aequ, radiuskm_pol);
  orthDef;
  Gnomonic;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstatlog(U, V, cov, v, sign);// nicht gatternr
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
			       double *sign) {
   assert(cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  orthDefStat;
  GnomonicStat;
  orthTrafoStat;
  CovList[cov->secondarygatternr].log(U, cov, v, sign);// nicht gatternr .
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
		     double *sign) {
   assert(cov->gatternr == EARTHMILES2GNOMONIC);
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  orthDef;
  Gnomonic;
  orthTrafo;
  CovList[cov->secondarygatternr].nonstatlog(U, V, cov, v, sign);// nicht gatternr
}

