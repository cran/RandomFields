/*
 Authors 
 Christoph Berreth, cberreth@mail.uni-mannheim.de

 Definition of spherical correlation functions 

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 



 Copyright (C) 2014 -- 2014 Christoph Berreth

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

/* header files */
#include <stdio.h> /* Standard Ein-/Ausgabefunktionen */
#include <math.h>  /* for fct acos(), floor() */
#include <R_ext/Constants.h>  /* PI, DOUBLE_EPS, etc */
#include "RF.h"    /* for struct cov_model, macro P0,
		    enthält u.A. basic.h für vordefinierte Konstanten*/

#include "SpherModels.h"   /* muss nach RF.h stehen, da cov_model 
			      enthalten. */



////////////////////////////////////////////////////////////////////////
/* model of the sine power function of Soubeyrand, Enjalbert and Sache,
   taken out of the paper written by Tilmann Gneiting 
   'Strictly and non-strictly positive definite functions on spheres'
   (released in Bernoulli 19(4) in 2013) */
 

int i=1;
#define SINEPOWER_ALPHA 0
void SinePower(double *x, cov_model *cov, double *v){
  double alpha = P0(SINEPOWER_ALPHA);  // Auslesen des Parameters aus cov  
  double y = *x;

  // printf("i= %d   x= %.2f\n",i++,y); // => RFsimu läuft Gitter ab.
  // => Argument ist Norm(x)
  // => Übergibt User Punkte auf SPhäre oder Gitter in  R_2??

/* User übergibt Punkte auf SPhäre anhand von Winkeln, bspw. durch
 x = seq(0,pi,length=10); y=x. Hierbei handelt es sich um den oberen
 Kugelauschnitt im 1. Quadranten. Die Cov-Fkt. benötigt zur Berechnung
 den zwischen zweier solcher Gitterpunkte liegenden Innenwinkel. Dieser 
 lässt sich mittels der Formel (siehe Block) berechnen. In Beispiel muss
 dieser (logischerweise) kleiner pi bzw. 90° sein. Möglich sind aber
 auch Innenwinkel jenseits 90°. Zur Berechnung der Cov.-Werte muss der 
 Rest der Division durch pi berechnet werden (mod in C über fmod)... 
*/ 
  

//  printf("Remainder of %f / %f is %lf\n", 9.2, 3.7, fmod(9.2,3.7));
  *v = 1.0;
  if( y >=  PI ) *v = 0.0;
  else  *v = 1.0 - pow(sin(y * 0.5), alpha);
  return;
}



/* range des Parameters alpha */
void rangeSinePower(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SINEPOWER_ALPHA] = 0.0;  // true range
  range->max[SINEPOWER_ALPHA] = 2.0;   
  range->pmin[SINEPOWER_ALPHA] = 0.0001; // practical range (assumed)
  range->pmax[SINEPOWER_ALPHA] = 2.0; 
  range->openmin[SINEPOWER_ALPHA] = true;  // lower endpoint included 
  range->openmax[SINEPOWER_ALPHA] = false; // upper endpoint excluded  
}



////////////////////////////////////////////////////////////////////////
/* model according to the equation (16), 
   taken out of the paper written by Tilmann Gneiting 
   'Strictly and non-strictly positive definite functions on spheres'
   (released in Bernoulli 19(4) in 2013)
   The special cases with fixed parameter tau at 0.5 and 1.5 are known as
   the inverse multiquadric and the Poisson spline.
 */
 

#define EQ16_DELTA 0
#define EQ16_TAU 1
#define EQ16_EPS 1e-7
#define QUAD(a) (a)*(a)
void Eq16(double *x, cov_model *cov, double *v){
  
  double delta = P0(EQ16_DELTA), // Auslesen der Parameter aus cov  
         tau = P0(EQ16_TAU); 
  double y = *x;
  
  *v = 1.0;
  if( y >=  PI ) *v = 0.0;  // ?? Passt das??
  else  *v = pow(1.0-delta, 2*tau) / pow(1+QUAD(delta)-2*delta*cos(y), tau);
  return;
}


/* range der Parameter delta & tau */
void rangeEq16(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[EQ16_DELTA] = 0.0;  // true range
  range->max[EQ16_DELTA] = 1.0;   
  range->pmin[EQ16_DELTA] = EQ16_EPS; // practical range (assumed)
  range->pmax[EQ16_DELTA] = 1.0 - EQ16_EPS; 
  range->openmin[EQ16_DELTA] = true;  // lower endpoint included 
  range->openmax[EQ16_DELTA] = true; // upper endpoint excluded 

 range->min[EQ16_TAU] = 0;  // true range
  range->max[EQ16_TAU] = RF_INF;   
  range->pmin[EQ16_TAU] = 0.0001; // practical range (assumed)
  range->pmax[EQ16_TAU] = 500; 
  range->openmin[EQ16_TAU] = true;  // lower endpoint included 
  range->openmax[EQ16_TAU] = true; // upper endpoint excluded  
}










