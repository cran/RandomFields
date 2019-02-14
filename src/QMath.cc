// This file has been created automatically by 'rfGenerateMaths'
#include <R.h>
#include <Rmath.h>
#include "questions.h"
#include "primitive.others.h"
#include "startGetNset.h"
void MathACos(double *x, model *cov, double *v){
MATH_DEFAULT
*v = ACOS(w[0]); 
}


void MathASin(double *x, model *cov, double *v){
MATH_DEFAULT
*v = ASIN(w[0]); 
}


void MathATan(double *x, model *cov, double *v){
MATH_DEFAULT
*v = ATAN(w[0]); 
}


void MathAtan2(double *x, model *cov, double *v){
MATH_DEFAULT
*v = ATAN2(w[0], w[1]); 
}


void MathCos(double *x, model *cov, double *v){
MATH_DEFAULT
*v = COS(w[0]); 
}

void MathSin(double *x, model *cov, double *v){
MATH_DEFAULT
*v = SIN(w[0]); 
}


void MathTan(double *x, model *cov, double *v){
MATH_DEFAULT
*v = TAN(w[0]); 
}


void MathAcosh(double *x, model *cov, double *v){
MATH_DEFAULT
#if defined _GLIBCXX_CMATH
  *v = acosh(w[0]); // ACOSH()
#else
 *v = LOG(w[0] + SQRT(w[0] * w[0] - 1.0));
#endif
}


void MathAsinh(double *x, model *cov, double *v){
MATH_DEFAULT
#if defined _GLIBCXX_CMATH
  *v = asinh(w[0]); // ASINH()
// printf("%10g %10g\n", *v,  LOG(w[0] + SQRT(w[0] * w[0] + 1.0)));
#else
 *v = LOG(w[0] + SQRT(w[0] * w[0] + 1.0));
#endif
}


void MathAtanh(double *x, model *cov, double *v){
MATH_DEFAULT
#if defined _GLIBCXX_CMATH
  *v = atanh(w[0]); // ATANH()
#else
 *v = 0.5 * LOG( (1 + w[0]) / (1 - w[0]) );
#endif

}
 

void MathCosh(double *x, model *cov, double *v){
MATH_DEFAULT
*v = COSH(w[0]); 
}


void MathSinh(double *x, model *cov, double *v){
MATH_DEFAULT
*v = SINH(w[0]); 
}


void MathTanh(double *x, model *cov, double *v){
MATH_DEFAULT
*v = TANH(w[0]); 
}


void MathExp(double *x, model *cov, double *v){
MATH_DEFAULT
*v = EXP(w[0]); 
}


void MathLog(double *x, model *cov, double *v){
MATH_DEFAULT
*v = LOG(w[0]); 
}


void MathExpm1(double *x, model *cov, double *v){
MATH_DEFAULT 
#if defined _GLIBCXX_CMATH 
  *v = expm1(w[0]); // EXPM1()
#else
 *v = EXP(w[0]) - 1.0;
#endif
}


void MathLog1p(double *x, model *cov, double *v){
MATH_DEFAULT
#if defined _GLIBCXX_CMATH 
  *v = log1p(w[0]); // LOG1P()
#else
 *v = LOG(1.0 + w[0]);
#endif
}


void MathLogb(double *x, model *cov, double *v){
MATH_DEFAULT
  *v = LOGB(w[0]);
//  printf("logb %10g\t%10g\n", *v, logb(w[0]));
}


void MathExp2(double *x, model *cov, double *v){
MATH_DEFAULT
*v = EXP2(w[0]);
// printf("e %10g\t%10g\n", *v, EXP2(w[0])); 
}


void MathLog2(double *x, model *cov, double *v){
MATH_DEFAULT
*v = LOG_BASE2(w[0]);
// printf("e %10g\t%10g\n", *v, LOG_BASE2(w[0]));  
}


void MathPow(double *x, model *cov, double *v){
MATH_DEFAULT
*v = POW(w[0], w[1]); 
}


void MathSqrt(double *x, model *cov, double *v){
MATH_DEFAULT
*v = SQRT(w[0]); 
}


void MathHypot(double *x, model *cov, double *v){
MATH_DEFAULT
#if defined _GLIBCXX_CMATH
  *v = hypot(w[0], w[1]);
#else
 *v = SQRT(w[0] * w[0] + w[1] * w[1]);
#endif
}


void MathCbrt(double *x, model *cov, double *v){
MATH_DEFAULT
*v = CBRT(w[0]);
// printf("e %10g\t%10g\n", *v, CBRT(w[0])); 
}


void MathCeil(double *x, model *cov, double *v){
MATH_DEFAULT
*v = CEIL(w[0]); 
}


void MathFABS(double *x, model *cov, double *v){
MATH_DEFAULT
*v = FABS(w[0]); 
}


void MathFloor(double *x, model *cov, double *v){
MATH_DEFAULT
*v = FLOOR(w[0]); 
}


void MathFmod(double *x, model *cov, double *v){
MATH_DEFAULT
  *v = FMOD(w[0], w[1]); 
//printf("e %10g\t%10g\n", *v, FMOD(w[0], w[1]));
}



void MathRound(double *x, model *cov, double *v){
MATH_DEFAULT
*v = ROUND(w[0]); 
}


void MathTrunc(double *x, model *cov, double *v){
MATH_DEFAULT
*v = TRUNC(w[0]); 
}


/*

void Mathnearbyint(double *x, model *cov, double *v){
MATH_DEFAULT
  *v = std::nearbyint(w[0]); 
}

void MathLrint(double *x, model *cov, double *v){
MATH_DEFAULT
  *v = std::lrint(w[0]); 
}


void MathLlrint(double *x, model *cov, double *v){
MATH_DEFAULT
*v = llrint(w[0]); 
}

void MathLRound(double *x, model *cov, double *v){
  MATH_DEFAULT
  *v = l ROUND(w[0]); 
}


void MathLLRound(double *x, model *cov, double *v){
MATH_DEFAULT
*v = ll ROUND(w[0]); 
}

void MathCopysign(double *x, model *cov, double *v){
MATH_DEFAULT
*v = copysign(w[0], w[1]); 
}

void MathRint(double *x, model *cov, double *v){
MATH_DEFAULT
*v = rint(w[0]); 
}

void MathNextafter(double *x, model *cov, double *v){
MATH_DEFAULT
*v = nextafter(w[0], w[1]); 
}


void MathNexttoward(double *x, model *cov, double *v){
MATH_DEFAULT
*v = nexttoward(w[0], w[1]); 
}
*/

void MathErf(double *x, model *cov, double *v){
MATH_DEFAULT
 *v = 1.0 - 2.0 * pnorm(w[0], 0.0, INVSQRTTWO, 0, 0);
}


void MathErfc(double *x, model *cov, double *v){
MATH_DEFAULT
  *v = 2.0 * pnorm(w[0], 0.0, INVSQRTTWO, 0, 0);
}


void MathGamma(double *x, model *cov, double *v){
MATH_DEFAULT
  *v = gammafn(w[0]); 
}


void MathLgamma(double *x, model *cov, double *v){
MATH_DEFAULT
  *v = lgammafn(w[0]);
}


void MathRemainder(double *x, model *cov, double *v){
MATH_DEFAULT
*v = REMAINDER(w[0], w[1]);
// printf("e %10g\t%10g\t%10g\n", *v, REMAINDER(w[0], w[1]), fr ound(w[1], 0));
}


void MathFdim(double *x, model *cov, double *v){
MATH_DEFAULT
*v = FDIM(w[0], w[1]);
// printf("e %10g\t%10g\n", *v, FDIM(w[0], w[1]));

}


void MathFmax(double *x, model *cov, double *v){
MATH_DEFAULT
*v = FMAX(w[0], w[1]);
 
}


void MathFmin(double *x, model *cov, double *v){
MATH_DEFAULT
*v = FMIN(w[0], w[1]); 
}


bool allowedImaths(model *cov) {
  int z = 0,
    k = DefList[COVNR].kappas;
  
  model *sub[MAXSUB];
  for (int i=0; i<k; i++)
    if (cov->kappasub[i] != NULL) sub[z++] = cov->kappasub[i];
  return allowedIsubs(cov, sub, z);
}


#define ADDCOV(X)				\
  addCov(X, NULL, NULL);			\
  AddVariant(TrendType, PREVMODEL_I);		\
  setDI(NULL, allowedImaths, NULL);


void includeStandardMath() {
  int first = currentNrCov;

IncludeModel(".asin", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("asin");
kappanames("x", REALSXP);
 ADDCOV(MathASin);

IncludeModel(".atan", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("atan");
kappanames("x", REALSXP);
ADDCOV(MathATan);

IncludeModel(".atan2", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("atan2");
kappanames("y", REALSXP, "x", REALSXP);
ADDCOV(MathAtan2);

IncludeModel(".cos", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("cos");
kappanames("x", REALSXP);
ADDCOV(MathCos);

IncludeModel(".sin", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("sin");
kappanames("x", REALSXP);
ADDCOV(MathSin);

IncludeModel(".tan", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("tan");
kappanames("x", REALSXP);
ADDCOV(MathTan);

IncludeModel(".asinh", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("asinh");
kappanames("x", REALSXP);
ADDCOV(MathAsinh);

IncludeModel(".atanh", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("atanh");
kappanames("x", REALSXP);
ADDCOV(MathAtanh);

IncludeModel(".cosh", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("cosh");
kappanames("x", REALSXP);
ADDCOV(MathCosh);

IncludeModel(".sinh", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("sinh");
kappanames("x", REALSXP);
ADDCOV(MathSinh);

IncludeModel(".tanh", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("tanh");
kappanames("x", REALSXP);
ADDCOV(MathTanh);

IncludeModel(".log", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("log");
kappanames("x", REALSXP);
ADDCOV(MathLog);

IncludeModel(".expm1", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("expm1");
kappanames("x", REALSXP);
ADDCOV(MathExpm1);

IncludeModel(".log1p", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("log1p");
kappanames("x", REALSXP);
ADDCOV(MathLog1p);

/*
IncludeModel(".logb", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	     !false, // macht ein haufen Aerger
	     SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("logb");
kappanames("x", REALSXP);
ADDCOV(MathLogb);
*/

IncludeModel(".exp2", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("exp2");
kappanames("x", REALSXP);
ADDCOV(MathExp2);

IncludeModel(".log2", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("log2");
kappanames("x", REALSXP);
ADDCOV(MathLog2);

IncludeModel(".hypot", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("hypot");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathHypot);

IncludeModel(".cbrt", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("cbrt");
kappanames("x", REALSXP);
ADDCOV(MathCbrt);

IncludeModel(".ceil", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("ceil");
kappanames("x", REALSXP);
ADDCOV(MathCeil);

IncludeModel(".floor", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("floor");
kappanames("x", REALSXP);
ADDCOV(MathFloor);

IncludeModel(".fmod", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("fmod");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathFmod);


IncludeModel(".round", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("round");
kappanames("x", REALSXP);
ADDCOV(MathRound);

IncludeModel(".trunc", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("trunc");
kappanames("x", REALSXP);
ADDCOV(MathTrunc);


/*
IncludeModel(".nearbyint", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, NOT_MONTONE); 
nickname("nearbyint");
kappanames("x", REALSXP);
ADDCOV(MathNearbyint);

IncludeModel(".lrint", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("lrint");
kappanames("x", REALSXP);
ADDCOV(Mathlrint);

IncludeModel(".llrint", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("llrint");
kappanames("x", REALSXP);
ADDCOV(MathLlrint);

IncludeModel(".lround", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("lround");
kappanames("x", REALSXP);
ADDCOV(MathLRound);

IncludeModel(".llround", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("llround");
kappanames("x", REALSXP);
ADDCOV(MathLLRound);

IncludeModel(".copysign", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, NOT_MONTONE); 
nickname("copysign");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathCopysign);

IncludeModel(".rint", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, NOT_MONTONE); 
nickname("rint");
kappanames("x", REALSXP);
ADDCOV(MathRint);

IncludeModel(".nextafter", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, NOT_MONTONE); 
nickname("nextafter");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathNextafter);

IncludeModel(".nexttoward", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, NOT_MONTONE); 
nickname("nexttoward");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathNexttoward);
*/

IncludeModel(".erfc", MathDefType, 0, 0, 1, NULL, XONLY,
	     PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	     false, SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE); 
nickname("erfc");
kappanames("x", REALSXP);
ADDCOV(MathErfc);

IncludeModel(".lgamma", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE); 
nickname("lgamma");
kappanames("x", REALSXP);
ADDCOV(MathLgamma);


IncludeModel(".remainder", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE); 
nickname("remainder");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathRemainder);

IncludeModel(".fdim", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE); 
nickname("fdim");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathFdim);

IncludeModel(".fmax", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE); 
nickname("fmax");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathFmax);

IncludeModel(".fmin", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE); 
nickname("fmin");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathFmin);

////////////////////////////////////////////////////////////////////////
 
for (int nr=first; nr<currentNrCov; nr++) {
  assert(isMathDef(DefList + nr));
   set_type(DefList[nr].systems[0], 0, ShapeType);
}

////////////////////////////////////////////////////////////////////////
 
IncludeModel(".gamma", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_MATHDEF,
	     false, SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE); 
nickname("gamma");
kappanames("x", REALSXP);
ADDCOV(MathGamma);


IncludeModel(".exp", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_MATHDEF,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("exp");
kappanames("x", REALSXP);
ADDCOV(MathExp);


IncludeModel(".erf", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_MATHDEF,
	false, SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE); 
nickname("erf");
kappanames("x", REALSXP);
ADDCOV(MathErf);

 
IncludeModel(".fabs", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_MATHDEF,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("fabs");
kappanames("x", REALSXP);
ADDCOV(MathFABS);

IncludeModel(".acos", MathDefType, 0, 0, 1, NULL, XONLY,
	     PREVMODEL_I,checkMath,rangeMath, PREF_MATHDEF,
	     false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("acos");
kappanames("x", REALSXP);
ADDCOV(MathACos);

IncludeModel(".acosh", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_MATHDEF,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("acosh");
kappanames("x", REALSXP);
ADDCOV(MathAcosh);

IncludeModel(".pow", MathDefType, 0, 0, 2, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_MATHDEF,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("pow");
kappanames("x", REALSXP, "y", REALSXP);
ADDCOV(MathPow);

IncludeModel(".sqrt", MathDefType, 0, 0, 1, NULL, XONLY,
	 PREVMODEL_I,checkMath,rangeMath, PREF_MATHDEF,
	false, SCALAR, PREVMODEL_DEP, (ext_bool) false, NOT_MONOTONE); 
nickname("sqrt");
kappanames("x", REALSXP);
ADDCOV(MathSqrt);

}
