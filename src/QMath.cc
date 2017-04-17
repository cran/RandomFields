// This file has been created automatically by 'rfGenerateMaths'
#include <Rmath.h>
#include "RF.h"
#include "primitive.h"
void MathACos(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = ACOS(w[0]); 
}


void MathASin(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = ASIN(w[0]); 
}


void MathATan(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = ATAN(w[0]); 
}


void MathAtan2(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = ATAN2(w[0], w[1]); 
}


void MathCos(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = COS(w[0]); 
}


void MathSin(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = SIN(w[0]); 
}


void MathTan(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = TAN(w[0]); 
}


void MathAcosh(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = ACOSH(w[0]); 
}


void MathAsinh(double *x, cov_model *cov, double *v){
MATH_DEFAULT
  *v = ASINH(w[0]); 
}


void MathAtanh(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = ATANH(w[0]); 
}


void MathCosh(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = COSH(w[0]); 
}


void MathSinh(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = SINH(w[0]); 
}


void MathTanh(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = TANH(w[0]); 
}


void MathExp(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = EXP(w[0]); 
}


void MathLog(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = LOG(w[0]); 
}


void MathExpm1(double *x, cov_model *cov, double *v){
MATH_DEFAULT
  *v = EXPM1(w[0]);
}


void MathLog1p(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = LOG1P(w[0]); 
}


void MathLogb(double *x, cov_model *cov, double *v){
MATH_DEFAULT
  *v = LOGB(w[0]);
//  printf("logb %f\t%f\n", *v, logb(w[0]));
}


void MathExp2(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = EXP(w[0]);
// printf("e %f\t%f\n", *v, EXP2(w[0])); 
}


void MathLog2(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = LOG_BASE2(w[0]);
// printf("e %f\t%f\n", *v, LOG_BASE2(w[0]));  
}


void MathPow(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = POW(w[0], w[1]); 
}


void MathSqrt(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = SQRT(w[0]); 
}


void MathHypot(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = HYPOT(w[0], w[1]); 
// printf("e %f\t%f\n", *v, HYPOT(w[0], w[1]));
}


void MathCbrt(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = CBRT(w[0]);
// printf("e %f\t%f\n", *v, CBRT(w[0])); 
}


void MathCeil(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = CEIL(w[0]); 
}


void MathFABS(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = FABS(w[0]); 
}


void MathFloor(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = FLOOR(w[0]); 
}


void MathFmod(double *x, cov_model *cov, double *v){
MATH_DEFAULT
  *v = FMOD(w[0], w[1]); 
//printf("e %f\t%f\n", *v, FMOD(w[0], w[1]));
}



void MathRound(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = ROUND(w[0]); 
}


void MathTrunc(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = TRUNC(w[0]); 
}


/*

void Mathnearbyint(double *x, cov_model *cov, double *v){
MATH_DEFAULT
  *v = std::nearbyint(w[0]); 
}

void MathLrint(double *x, cov_model *cov, double *v){
MATH_DEFAULT
  *v = std::lrint(w[0]); 
}


void MathLlrint(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = llrint(w[0]); 
}

void MathLRound(double *x, cov_model *cov, double *v){
  MATH_DEFAULT
  *v = l ROUND(w[0]); 
}


void MathLLRound(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = ll ROUND(w[0]); 
}

void MathCopysign(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = copysign(w[0], w[1]); 
}

void MathRint(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = rint(w[0]); 
}

void MathNextafter(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = nextafter(w[0], w[1]); 
}


void MathNexttoward(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = nexttoward(w[0], w[1]); 
}
*/

void MathErf(double *x, cov_model *cov, double *v){
MATH_DEFAULT
 *v = 1.0 - 2.0 * pnorm(w[0], 0.0, INVSQRTTWO, 0, 0);
}


void MathErfc(double *x, cov_model *cov, double *v){
MATH_DEFAULT
  *v = 2.0 * pnorm(w[0], 0.0, INVSQRTTWO, 0, 0);
}


void MathGamma(double *x, cov_model *cov, double *v){
MATH_DEFAULT
  *v = gammafn(w[0]); 
}


void MathLgamma(double *x, cov_model *cov, double *v){
MATH_DEFAULT
  *v = lgammafn(w[0]);
}


void MathRemainder(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = REMAINDER(w[0], w[1]);
// printf("e %f\t%f\t%f\n", *v, REMAINDER(w[0], w[1]), fr ound(w[1], 0));
}


void MathFdim(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = FDIM(w[0], w[1]);
// printf("e %f\t%f\n", *v, FDIM(w[0], w[1]));

}


void MathFmax(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = FMAX(w[0], w[1]);
 
}


void MathFmin(double *x, cov_model *cov, double *v){
MATH_DEFAULT
*v = FMIN(w[0], w[1]); 
}


void includeStandardMath() {
IncludeModel(".acos", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("acos");
kappanames("a", REALSXP);
addCov(MathACos, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".asin", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("asin");
kappanames("a", REALSXP);
addCov(MathASin, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".atan", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("atan");
kappanames("a", REALSXP);
addCov(MathATan, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".atan2", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("atan2");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathAtan2, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".cos", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("cos");
kappanames("a", REALSXP);
addCov(MathCos, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".sin", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("sin");
kappanames("a", REALSXP);
addCov(MathSin, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".tan", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("tan");
kappanames("a", REALSXP);
addCov(MathTan, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".acosh", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("acosh");
kappanames("a", REALSXP);
addCov(MathAcosh, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".asinh", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("asinh");
kappanames("a", REALSXP);
addCov(MathAsinh, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".atanh", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("atanh");
kappanames("a", REALSXP);
addCov(MathAtanh, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".cosh", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("cosh");
kappanames("a", REALSXP);
addCov(MathCosh, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".sinh", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("sinh");
kappanames("a", REALSXP);
addCov(MathSinh, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".tanh", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("tanh");
kappanames("a", REALSXP);
addCov(MathTanh, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".exp", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("exp");
kappanames("a", REALSXP);
addCov(MathExp, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".log", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("log");
kappanames("a", REALSXP);
addCov(MathLog, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".expm1", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("expm1");
kappanames("a", REALSXP);
addCov(MathExpm1, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".log1p", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("log1p");
kappanames("a", REALSXP);
addCov(MathLog1p, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".logb", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("logb");
kappanames("a", REALSXP);
addCov(MathLogb, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".exp2", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("exp2");
kappanames("a", REALSXP);
addCov(MathExp2, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".log2", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("log2");
kappanames("a", REALSXP);
addCov(MathLog2, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".pow", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("pow");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathPow, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".sqrt", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("sqrt");
kappanames("a", REALSXP);
addCov(MathSqrt, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".hypot", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("hypot");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathHypot, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".cbrt", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("cbrt");
kappanames("a", REALSXP);
addCov(MathCbrt, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".ceil", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("ceil");
kappanames("a", REALSXP);
addCov(MathCeil, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".fabs", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("fabs");
kappanames("a", REALSXP);
addCov(MathFABS, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".floor", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("floor");
kappanames("a", REALSXP);
addCov(MathFloor, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".fmod", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("fmod");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathFmod, NULL, NULL);
AddVariant(TrendType, PREVMODELI);


IncludeModel(".round", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("round");
kappanames("a", REALSXP);
addCov(MathRound, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".trunc", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("trunc");
kappanames("a", REALSXP);
addCov(MathTrunc, NULL, NULL);
AddVariant(TrendType, PREVMODELI);


/*
IncludeModel(".nearbyint", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("nearbyint");
kappanames("a", REALSXP);
addCov(MathNearbyint, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".lrint", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("lrint");
kappanames("a", REALSXP);
addCov(MathLrint, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".llrint", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("llrint");
kappanames("a", REALSXP);
addCov(MathLlrint, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".lround", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("lround");
kappanames("a", REALSXP);
addCov(MathLRound, NULL, NULL);
AddVariant(TrendType, PREVMODELI);


IncludeModel(".llround", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("llround");
kappanames("a", REALSXP);
addCov(MathLLRound, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".copysign", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("copysign");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathCopysign, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".rint", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("rint");
kappanames("a", REALSXP);
addCov(MathRint, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".nextafter", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("nextafter");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathNextafter, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".nexttoward", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("nexttoward");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathNexttoward, NULL, NULL);
AddVariant(TrendType, PREVMODELI);
*/

IncludeModel(".erf", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("erf");
kappanames("a", REALSXP);
addCov(MathErf, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".erfc", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("erfc");
kappanames("a", REALSXP);
addCov(MathErfc, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".gamma", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("gamma");
kappanames("a", REALSXP);
addCov(MathGamma, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".lgamma", MathDefinition, 0, 0, 1, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("lgamma");
kappanames("a", REALSXP);
addCov(MathLgamma, NULL, NULL);
AddVariant(TrendType, PREVMODELI);


IncludeModel(".remainder", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("remainder");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathRemainder, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".fdim", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("fdim");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathFdim, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".fmax", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("fmax");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathFmax, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

IncludeModel(".fmin", MathDefinition, 0, 0, 2, NULL, XONLY,
	 PREVMODELI,checkMath,rangeMath, PREF_TREND,
	false, SCALAR, PREVMODEL_DEP, false, false); 
nickname("fmin");
kappanames("a", REALSXP, "b", REALSXP);
addCov(MathFmin, NULL, NULL);
AddVariant(TrendType, PREVMODELI);

}
