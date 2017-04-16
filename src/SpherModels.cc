/*
 Authors 
 Christoph Berreth, cberreth@mail.uni-mannheim.de
 Martin Schlather

 Definition of spherical correlation functions 

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 



 Copyright (C) 2014 -- 2015 Christoph Berreth
               2016 -- 2017 Martin Schlather

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
#include <Rmath.h>  /* for fct ACOS(), FLOOR() */
#include <R_ext/Constants.h>  /* PI, DOUBLE_EPS, etc */
#include "RF.h"    /* for struct cov_model, macro P0,
		    enthält u.A. basic.h für vordefinierte Konstanten*/

#include "SpherModels.h"   /* muss nach RF.h stehen, da cov_model 
			      enthalten. */



////////////////////////////////////////////////////////////////////////
/* model of the sine power function of Soubeyrand, Enjalbert and Sache,
   taken out of the paper written by Tilmann Gneiting 
   'Strictly and non-strictly positive definite functions on spheres'
   (released in Bernoulli 19(4) in 2013 -- equation (17) ) */
 

int i=1;
#define SINEPOWER_ALPHA 0
void SinePower(double *x, cov_model *cov, double *v){
  double alpha = P0(SINEPOWER_ALPHA);  // Auslesen des Parameters aus cov  
  double y = *x;
  *v = y >=  PI ? 0.0 : 1.0 - POW(SIN(y * 0.5), alpha);
  return;
}



/* range des Parameters alpha */
void rangeSinePower(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SINEPOWER_ALPHA] = 0.0;  // true range
  range->max[SINEPOWER_ALPHA] = 2.0;   
  range->pmin[SINEPOWER_ALPHA] = 0.001; // practical range (assumed)
  range->pmax[SINEPOWER_ALPHA] = 2.0; 
  range->openmin[SINEPOWER_ALPHA] = true;  // lower endpoint included 
  range->openmax[SINEPOWER_ALPHA] = false; // upper endpoint excluded  
}



////////////////////////////////////////////////////////////////////////
/* model of the multiquadric familiy, 
   taken out of the paper written by Tilmann Gneiting 
   'Strictly and non-strictly positive definite functions on spheres'
   (released in Bernoulli 19(4) in 2013 -- equation (16) )
   The special cases with fixed parameter tau at 0.5 and 1.5 are known as
   the inverse multiquadric and the Poisson spline. */
 

#define MULTIQUAD_DELTA 0
#define MULTIQUAD_TAU 1
#define MULTIQUAD_EPS 1e-7
void Multiquad(double *x, cov_model *cov, double *v){
  
  double delta = P0(MULTIQUAD_DELTA), // Auslesen der Parameter aus cov  
         tau = P0(MULTIQUAD_TAU); 
  double y = *x;
  
  y = y >= PI ? -1 : COS(y);
  *v = POW(1.0-delta, 2*tau) / POW(1+ delta * delta - 2*delta*y, tau);
  return;
}


/* range der Parameter delta & tau */
void rangeMultiquad(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[MULTIQUAD_DELTA] = 0.0;  // true range
  range->max[MULTIQUAD_DELTA] = 1.0;   
  range->pmin[MULTIQUAD_DELTA] = MULTIQUAD_EPS; // practical range (assumed)
  range->pmax[MULTIQUAD_DELTA] = 1.0 - MULTIQUAD_EPS; 
  range->openmin[MULTIQUAD_DELTA] = true;  // lower endpoint included 
  range->openmax[MULTIQUAD_DELTA] = true; // upper endpoint excluded 

  range->min[MULTIQUAD_TAU] = 0;  // true range
  range->max[MULTIQUAD_TAU] = RF_INF;   
  range->pmin[MULTIQUAD_TAU] = 0.0001; // practical range (assumed)
  range->pmax[MULTIQUAD_TAU] = 500; 
  range->openmin[MULTIQUAD_TAU] = true;  // lower endpoint included 
  range->openmax[MULTIQUAD_TAU] = true; // upper endpoint excluded  
}



////////////////////////////////////////////////////////////////////////
/* Schoenberg's representation for d=2 respectively Choquet's representation, 
   taken out of the paper written by Tilmann Gneiting 
   'Strictly and non-strictly positive definite functions on spheres'
   (released in Bernoulli 19(4) in 2013 -- equation (15) )
*/ 


double Legendre(int n, double x)
{  
  /* Idea by Kai Zhang, 
     Course Math225 Scientific Computing II, 
     Yale, Spring 2011 
  */
  
  double Pk_1,Pk_2,Pk; // P_{k-1}(x), P_{k-2}(x), P_k(x)

  Pk_2 = 0.0;
  Pk_1 = 1.0;
  Pk = 1.0;
  
  for(int k=1;k<=n;k++)
    {
     Pk = (2.0*k-1.0)/k*x*Pk_1 - (k-1.0)/k*Pk_2; 
     Pk_2 = Pk_1;
     Pk_1 = Pk;
    }
 
  return Pk;
}


/* Bsp. */

//  int N=5;
//  printf("%lf\n", Legendre(N,-1));
//  printf("%lf\n", Legendre(N,0.5));
//  printf("%lf\n", Legendre(N,0));
//  printf("%lf\n", Legendre(N,0.5));
//  printf("%lf\n", Legendre(N,1));



/*
void kappa_choquet(int i, cov_model *cov, int *nr, int *nc){
    *nc = SIZE_NOT_DETERMINED;
      *nr = i < CovList[cov->nr].kappas  ? SIZE_NOT_DETERMINED : -1;
  }

#define CHOQUET_PARA 0
void Choquet(double *x, cov_model *cov, double *v){
  
  // Benutzer übergibt Gewichte in double-Vektor b, sodass  sum(b)=1
  // Übergabe mittels cov ??  
  // Überprüfung summe = 1 ??
  
  long int n; 
  double y = *x;
  double b[4];
  b[0]=b[1]=b[2]=b[3]=0.25; //???
  
  n = sizeof(*b)/ sizeof(double);
  
  *v = 0.0;
  if(y>=PI) for(int k=0; k<n;k++) *v += b[k]/(k+1)*Legendre(0,COS(PI));
  // für y>=PI stetige Fortsetzung mit konstantem Wert psi(PI)!!
  else for(int k=0; k<n;k++) *v += b[k]/(k+1)*Legendre(0,COS(y));
  return;
}





void rangeChoquet(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range)
{

  range->min[CHOQUET_PARA] = 0.0;  // true range
  range->max[CHOQUET_PARA] = 1.0;   
  range->pmin[CHOQUET_PARA] = MULTIQUAD_EPS; // practical range (assumed)
  range->pmax[CHOQUET_PARA] = 1.0 - MULTIQUAD_EPS; 
  range->openmin[CHOQUET_PARA] = true;  // lower endpoint included 
  range->openmax[CHOQUET_PARA] = true; // upper endpoint excluded 

  }
// * Erklärung:
// * range-Angaben gelten bei Vektoren für jede Komponente  
// * true range ist wahrer Definitionsbereich der Parameter
// * practical range sind Grenzen für Optimierer. (Optimierungsintervall)
// *

int checkChoquet(cov_model *cov){

  


  return NOERROR;
}

*/
