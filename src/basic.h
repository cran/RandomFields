

/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 Martin Schlather

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





#ifndef GSL_VS_R_H
#define GSL_VS_R_H 1

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <errno.h>
#include <R_ext/Complex.h>

#ifdef SCHLATHERS_MACHINE
#undef SCHLATHERS_MACHINE
#endif

#ifdef RANDOMFIELDS_DEBUGGING
#undef RANDOMFIELDS_DEBUGGING
#endif

#ifdef showfree
#undef showfree
#endif

#ifdef DOPRINT
#undef DOPRINT
#endif



#define MAXCHAR 18 // max number of characters for (covariance) names  
#define MAXNRCOVFCTS 300
#define MAXUNITSCHAR 10
#define MAXINVERSIONS 2


// Formerly in <R_ext/Applic.h>LinkedTo: 
void fft_factor_(int n, int *pmaxf, int *pmaxp);
Rboolean fft_work_(double *a, double *b, int nseg, int n, int nspn,
		  int isn, double *work, int *iwork);/* TRUE: success */

#define MIN(A,B) ((A) < (B) ? (A) : (B));
#define MAX(A,B) ((A) > (B) ? (A) : (B));

#define LENGTH length // safety, in order not to use LENGTH defined by R
#define complex Rcomplex
#define DOT "."
#define GAUSS_RANDOM(SIGMA) rnorm(0.0, SIGMA)
#define UNIFORM_RANDOM unif_rand()
#define POISSON_RANDOM(x) rpois(x)
#define SQRT2 M_SQRT2
#define SQRTPI M_SQRT_PI
#define INVPI M_1_PI
#define PIHALF M_PI_2 
#define ONETHIRD 0.333333333333333333
#define TWOTHIRD 0.66666666666666666667
#define TWOPI 6.283185307179586476925286766559
#define INVLOG2 1.442695040888963
#define INVSQRTTWO 0.70710678118654752440084436210
#define INVSQRTTWOPI 0.39894228040143270286
#define SQRTTWOPI 2.5066282746310002416
#define SQRTINVLOG005 0.5777613700268771079749
//#define LOG05 -0.69314718055994528623
#define LOG3 1.0986122886681096913952452369225257046474905578227
#define LOG2 M_LN2

#define EPSILON     0.00000000001
#define EPSILON1000 0.000000001
#define MAXINT 2147483647
#define INFDIM MAXINT
#define INFTY INFDIM



//
// 
// 
// 



#endif /* GSL_VS_R_H */


