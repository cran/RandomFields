


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 -- 2017  Martin Schlather

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


#ifndef Coordinates_H
#define Coordinates_H 1


#include "RF.h"



#define piD180 (M_PI * 0.00555555555555555555555)
#define H80Dpi (180.0 * INVPI)

void stat2(double *x, model *cov, double *v);
void nonstat2(double *x, double *y, model *cov, double *v);
void logstat2(double *x, model *cov, double *v, double *Sign);
void lognonstat2(double *x, double *y, model *cov, double *v, 
		     double *Sign);

void D_2(double *x, model *cov, double *v);
void DD_2(double *x, model *cov, double *v);
//int check2(model *cov);
void inverse2(double *x, model *cov, double *v);
void nonstatinverse2(double *v, model *cov, double *left, double *right);
void nonstat_loginverse2(double *v, model *cov, double *x, double *y);
//int struct2(model *cov, model **newmodel);
int struct2(model *cov, model **newmodel);
int init2(model *cov, gen_storage *s);
void do2(model *cov, gen_storage *s);
void dorandom2(model *cov, double *v);

void EarthIso2EarthIso(double *x);
void Earth2Earth(double *x);
void NonstatEarth2EarthIso(double *x, double *y, model *cov, double *v);

void NonstatEarth2Earth(double *x, double *y, model *cov, double *v, double *);

void EarthIso2SphereIso(double *x, model *cov, double *v);
void NonstatEarth2SphereIso(double *x, double *y, model *cov, double *v);

void Earth2Sphere(double *x, model *cov, double *v);
void NonstatEarth2Sphere(double *x, double *y, model *cov, double *v, double *);

void SphereIso2SphereIso(double *x);
void NonstatSphere2SphereIso(double *x, double *y, model *cov, double *v);

void Sphere2Sphere(double *x);
void NonstatSphere2Sphere(double *x, double *y, model *cov, double *v, double *);

void EarthKM2CartStat(double *x, model *cov, double *v);
void EarthKM2Cart(double *x, double *y, model *cov, double *v, double *);

void EarthMiles2CartStat(double *x, model *cov, double *v);
void EarthMiles2Cart(double *x, double *y, model *cov, double *v, double *);

int checkEarth(model *cov);

void EarthKM2OrthogStat(double *x, model *cov, double *v);
void EarthKM2Orthog(double *x, double *y, model *cov, double *v, double *);

void EarthMiles2OrthogStat(double *x, model *cov, double *v);
void EarthMiles2Orthog(double *x, double *y, model *cov, double *v, double *);


void Earth2GnomonicStat(double *x, model *cov, double *v);
void Earth2Gnomonic(double *x, double *y, model *cov, double *v, double *);



double mod(double x, double modulus);
double lonmod(double x, double modulus); 
double latmod(double x, double modulus);
void statmod2(double *x, double lon, double lat, double *y);


// ACHTUNG! Nachfolgend nicht (ZZ) !!
#define STATMODE_BASE(X, ZZ, lon, lat)				\
  (X)[0]=lonmod(ZZ[0], lon);					\
  (X)[1]=latmod(ZZ[1], lat)					\

void statmod2(double *x, double lon, double lat, double *y); // 2dim only



#define radiuskm_aequ 6378.1
#define radiuskm_pol 6356.8
#define radiusmiles_aequ 3963.17
#define radiusmiles_pol 3949.93

void Earth2Cart(model *cov, double RAEQU, double RPOL, double *y);


coord_sys_enum GetCoordSystem(isotropy_type iso);


isotropy_type CoordinateSystemOf(isotropy_type iso);
isotropy_type EssentialCoordinateSystemOf(isotropy_type iso);
isotropy_type SymmetricOf(isotropy_type iso);
isotropy_type IsotropicOf(isotropy_type iso);

coord_sys_enum SearchCoordSystem(model *cov, coord_sys_enum os, 
				 coord_sys_enum n_s);



#endif
