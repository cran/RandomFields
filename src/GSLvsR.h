

#ifndef GSL_VS_R_H
#define GSL_VS_R_H 1

#ifdef RF_GSL /* RF_GSL */
NOT ALLOWED ANYMORE !
#endif



#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

//#define 
//#define 
//#define  
#define EXTERN extern "C"
#define PRINTF Rprintf
#define RF_NAN NA_REAL 
#define RF_NEGINF R_NegInf
#define RF_INF R_PosInf
#define RF_ISNA ISNAN
#define GAUSS_RANDOM(SIGMA) rnorm(0.0, SIGMA)
#define UNIFORM_RANDOM unif_rand()
#define POISSON_RANDOM(x) rpois(x)
#define SQRT2 M_SQRT2
#define SQRTPI M_SQRT_PI
#define INVPI M_1_PI
#define RF_M_SQRT_3 M_SQRT_3
//#define TWOPI M_2PI
#define TWOPI 6.283185307179586476925286766559
#define PIHALF M_PI_2 
#define T_PI M_2_PI

#endif /* GSL_VS_R_H */


