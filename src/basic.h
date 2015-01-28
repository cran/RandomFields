#ifndef GSL_VS_R_H
#define GSL_VS_R_H 1

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <errno.h>
#include <R_ext/Complex.h>


// Formerly in <R_ext/Applic.h>
void fft_factor_(int n, int *pmaxf, int *pmaxp);
Rboolean fft_work_(double *a, double *b, int nseg, int n, int nspn,
		  int isn, double *work, int *iwork);/* TRUE: success */

#define MIN(A,B) ((A) < (B) ? (A) : (B));
#define MAX(A,B) ((A) > (B) ? (A) : (B));

#define LENGTH length // safety, in order not to use LENGTH defined by R
#define complex Rcomplex
#define DOT "."
#define print PRINTF /* // */
#define RF_NA NA_REAL 
#define RF_NAN R_NaN
#define RF_NEGINF R_NegInf
#define RF_INF R_PosInf
#define GAUSS_RANDOM(SIGMA) rnorm(0.0, SIGMA)
#define UNIFORM_RANDOM unif_rand()
#define POISSON_RANDOM(x) rpois(x)
#define SQRT2 M_SQRT2
#define SQRTPI M_SQRT_PI
#define INVPI M_1_PI
#define PIHALF M_PI_2 
#define T_PI M_2_PI
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


#define ERRLINES assert({PRINTF("(ERROR in %s, line %d)\n", __FILE__, __LINE__); true;})

//
// 
// 
// 

// 1

extern char BUG_MSG[250];
#ifdef LOCAL_MACHINE
// __extension__ unterdrueckt Fehlermeldung wegen geklammerter Argumente
#define PRINTF Rprintf
#define INTERNAL  \
  sprintf(BUG_MSG, \
	  "made to be an internal function '%s' ('%s', line %d).", /* // */ \
	  __FUNCTION__, __FILE__, __LINE__);				\
  warning(BUG_MSG)
 
#define assert(X) if (!__extension__ (X )) {				\
    sprintf(BUG_MSG,							\
	    "'assert(%s)' failed in function '%s' (file '%s', line %d).", \
	    #X,__FUNCTION__, __FILE__, __LINE__);			\
    error(BUG_MSG);							\
  }
#define VARIABLE_IS_NOT_USED __attribute__ ((unused))
#define SHOW_ADDRESSES 1
#define BUG {								\
    sprintf(BUG_MSG, "BUG in '%s' ('%s', line %d).", \
	    __FUNCTION__, __FILE__, __LINE__);				\
    error(BUG_MSG);							\
  }									
#define DO_TESTS true

#else 


#define PRINTF Rprintf
#define INTERNAL SERR("Sorry. This functionality does not exist currently. There is work in progress at the moment by the maintainer.")
#define assert(X) {}
#define VARIABLE_IS_NOT_USED
#define BUG {								\
    sprintf(BUG_MSG, "Severe error occured in function '%s' (file '%s', line %d). Please contact maintainer martin.schlather@math.uni-mannheim.de .", \
	    __FUNCTION__, __FILE__, __LINE__);				\
    error(BUG_MSG);							\
  }									
#define DO_TESTS false

#endif

#endif /* GSL_VS_R_H */


