
#ifndef GSL_VS_R_H
#define GSL_VS_R_H 1

#ifdef RF_GSL /* RF_GSL */

#ifndef MATHLIB_STANDALONE
#define MATHLIB_STANDALONE 1
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#define RANDOMNUMBERGENERATOR gsl_rng_mt19937
//                            gsl_rng_taus,...
#define PI 3.14159265358979323846264338328
#define END_WITH_RANDOM ,RANDOM
#define END_WITH_GSLRNG ,gsl_rng *RANDOM
#define END_WITH_GSL ,gsl_rng *
#define UNIFORM_RANDOM gsl_rng_uniform(RANDOM)
#define PRINTF printf
#define RF_NAN 0/0
#define RF_INF 1E308
#define RF_NEGINF -1E308
#define RF_ISNA ********
extern gsl_rng *RANDOM;



#else        /* RF_GSL */
#define CPP2C 1

#include <R.h>
#include <Rinternals.h>

#define END_WITH_RANDOM 
#define END_WITH_GSLRNG
#define END_WITH_GSL
#define UNIFORM_RANDOM unif_rand()
#define PRINTF Rprintf
#define RF_NAN NA_REAL 
#define RF_INF R_PosInf
#define RF_NEGINF R_NegInf
#define RF_ISNA ISNAN
#endif      /* RF_GSL */

#endif /* GSL_VS_R_H */

