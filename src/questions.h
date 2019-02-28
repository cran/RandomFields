



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


#ifndef Questions_H
#define Questions_H 1

#include "RF.h"
#include "primitive.others.h"

typedef bool (*systemfct)(system_type *sys); /* h, cov, result */ 
typedef bool (*typefct)(Types type); /* h, cov, result */ 

bool equalsXonly(domain_type dom); 
bool isXonly(system_type *s); 
bool equalsKernel(domain_type dom); 
bool isKernel(system_type *s); 
bool equalsDomMismatch(domain_type dom); 
bool equalsPrevModelD(domain_type dom);
bool isPrevModelD(defn *C);
bool equalsParamDepD(domain_type dom);
bool isParamDepD(defn *C);
bool equalsSubModelD(domain_type dom);
bool isSubModelD(defn *C);
bool isFixed(domain_type dom);

//bool isIsotropicXonly(system_type s);
bool isIsotropicXonly(system_type *s);
bool isIsotropicXonly(model *cov); 
//bool isIsoCovFct(Systems_type s);
//bool isIsoCovFct(model *cov); 
bool equalsIsotropic(isotropy_type iso); 
//bool isIsotropic(system_type sys);
bool isIsotropic(system_type *sys);
bool isSpaceIsotropic(isotropy_type iso);
bool isSpaceIsotropic(system_type *sys);
bool equalsSpaceIsotropic(isotropy_type iso);
bool equalsSpaceIsotropic(system_type* sys);
 
bool isAnyIsotropic(isotropy_type iso);
bool equalsVectorIsotropic(isotropy_type iso); 
bool equalsSymmetric(isotropy_type iso); // ok
bool isSymmetric(isotropy_type iso);  // ok
bool equalsCartCoord(isotropy_type iso);
bool isCartesian(isotropy_type iso); 
bool isCartesian(system_type *sys, int s);
bool isCartesian(system_type *sys);
bool isGenuineCartesian(isotropy_type iso);
bool hasFullCartesian(isotropy_type iso);
//bool isUnreducedCartesian(isotropy_type iso);
//bool isUnreducedSpherical(isotropy_type iso);
//bool isUnreducedEarth(isotropy_type iso);
//bool isUnreducedLogCart(isotropy_type iso);
//bool isUnreduced(isotropy_type iso);

bool equalsGnomonic(isotropy_type iso); 
bool equalsOrthographic(isotropy_type iso);
bool isEarthProjection(isotropy_type iso);

bool equalsSphericalIsotropic(isotropy_type iso);
bool equalsSphericalSymmetric(isotropy_type iso);
bool equalsSphericalCoord(isotropy_type iso);
bool isSphericalSymmetric(isotropy_type iso);
bool equalsEarthIsotropic(isotropy_type iso);
bool equalsEarthSymmetric(isotropy_type iso);
bool equalsEarthCoord(isotropy_type iso);
bool isEarthSymmetric(isotropy_type iso);
bool isSpherical(isotropy_type iso);
bool isEarth(isotropy_type iso);
bool isAnySpherical(isotropy_type iso);
bool isAnySphericalIso(isotropy_type iso);
bool isLogCart(isotropy_type iso);
//bool isLogCart(system_type sys); 
bool isLogCart(system_type *sys, int s);
bool isCylinder(isotropy_type iso);
bool equalsAnySymmetric(isotropy_type iso);


bool equalsPrevModelI(isotropy_type iso); 
bool isPrevModelI(defn *C);
bool equalsSubModelI(isotropy_type iso); 
bool isSubModelI(defn *C);
bool equalsParamDepI(isotropy_type iso);
bool isParamDepI(defn *C);
bool equalsIsoMismatch(isotropy_type iso);
bool equalsUnreduced(isotropy_type iso); 
bool equalsUnreduced(defn *C);
bool equalsCoordinateSystem(isotropy_type iso);
bool isFixed(isotropy_type iso); 


bool HaveSameSystems(system_type * fst, system_type * snd);
bool everyCoord(systemfct iso, model *cov);
bool anyVariant(systemfct sf, defn *C);
bool anyVariant(typefct tf, defn *C);



bool isDummyInit(initfct Init);
bool isMtimesep(matrix_type type);
bool isMproj(matrix_type type);
bool isMdiag(matrix_type type);
bool isMiso(matrix_type type);

bool isTcf(Types type);
bool isPosDef(Types type);
bool isNegDef(Types type);
bool isVariogram(Types type);
bool isProcess(Types type);
bool isNormedProcess(Types type);
bool isGaussMethod(Types type);
bool isBrMethod(Types type);
bool isPointShape(Types type);
bool isRandom(Types type);
bool isShape(Types type);
bool isTrend(Types type);
bool isInterface(Types type);
bool isMathDef(Types type);
bool isSameAsPrev(Types type);
bool isPoisson(Types type);
bool isSchlather(Types type);
bool isPoissonGauss(Types type);
bool isBad(Types type);
bool isManifold(Types type);
//bool isOther(Types type);
//bool isLikelihood(Types type);

bool isManifold(defn *C);
bool isMathDef(defn *C);

// dangerous to add others!!
// see init.others.cc: CheckListmodel & questions.cc:isDefCL
bool isProcess(model *cov);
bool isNormedProcess(model *cov);
bool isGaussMethod(model *cov);
bool isBrMethod(model *cov);
bool isInterface(model *cov);
bool isPointShape(model *cov);
//bool isShape(model *cov);
bool isRandom(model *cov);
bool isTrend(model *cov);

bool isnowTcf(model *cov);
bool isnowPosDef(model *cov);
bool isnowNegDef(model *cov);
bool isnowVariogram(model *cov);
bool isnowProcess(model *cov);
bool isnowPoisson(model *cov);
bool isnowMaxStable(model *cov);
//bool isnowGaussMethod(model *cov);
bool isnowPoissonGauss(model *cov);
bool isnowRandom(model *cov);
bool isnowShape(model *cov);
bool isnowTrend(model *cov);
bool isnowSchlather(model *cov);
//bool isnowInterface(model *cov);
//bool isnowManifold(model *cov);
//bool isnowSameAsPrev(model *cov);
//bool isnowMathDef(model *cov);
//bool isnowBad(model *cov);
//bool isnowOther(model *cov);

bool equalsnowPosDef(model *cov);
bool equalsnowNegDef(model *cov);
bool equalsnowVariogram(model *cov);
bool equalsnowProcess(model *cov);
bool equalsnowPoisson(model *cov);
bool equalsnowPoissonGauss(model *cov);
bool equalsnowGaussMethod(model *cov);
bool equalsnowPointShape(model *cov);
bool equalsnowRandom(model *cov);
// bool equalsnowShape(model *cov);
bool equalsnowTrend(model *cov);
bool equalsnowInterface(model *cov);
bool equalsnowSameAsPrev(model *cov);
bool equalsnowBad(model *cov);
bool equalsnowMathDef(model *cov);
bool equalsnowSchlather(model *cov);
bool equalsnowSmith(model *cov);
//bool equalsnowOther(model *cov);

bool equalsPosDef(Types type);
bool equalsNegDef(Types type);
bool equalsVariogram(Types type);
bool equalsProcess(Types type);
bool equalsNormedProcess(Types type);
bool equalsGaussMethod(Types type);
bool equalsPointShape(Types type);
bool equalsRandom(Types type);
bool equalsShape(Types type);
bool equalsTrend(Types type);
bool equalsInterface(Types type);
bool equalsSameAsPrev(Types type);
bool equalsSmith(Types type);
bool equalsSchlather(Types type);
bool equalsBrMethod(Types type);
bool equalsPoisson(Types type);
bool equalsPoissonGauss(Types type);
bool equalsBad(Types type);
//bool equalsOther(Types type);
bool equalsLikelihood(Types type);
bool equalsEvaluation(Types type);


bool hasGaussMethodFrame(model *cov);
bool hasNormedProcessFrame(model *cov);
bool hasProcessFrame(model *cov);
bool hasPoissonFrame(model *cov);
bool hasPoissonGaussFrame(model *cov);
bool hasRandomFrame(model *cov);
bool hasInterfaceFrame(model *cov);
bool hasTrendFrame(model *cov);
bool hasEvaluationFrame(model *cov);
bool hasLikelihoodFrame(model *cov);
bool hasBrMethodFrame(model *cov);
bool hasSmithFrame(model *cov);
bool hasSchlatherFrame(model *cov);
bool hasMaxStableFrame(model *cov);
bool hasAnyPoissonFrame(model *cov);
bool hasAnyProcessFrame(model *cov);
bool hasAnyEvaluationFrame(model *cov);


bool isMonotone(model *cov);
bool isMonotone(monotone_type montone);
bool isNormalMixture(monotone_type montone);
bool isNormalMixture(model *cov);
bool isBernstein(model *cov);
bool isBernstein(monotone_type monotone);
bool isCompletelyMonotone(monotone_type monotone);
bool isGneiting(model *cov);
bool isGneiting(monotone_type monotone);


bool isRObject(int type);
bool allowsRandomParam(model *cov, int i);
bool allowsShapeParam(model *cov, int i);
bool isnowTrendParam(model *cov, int i);

bool isverysimple(model *cov);

bool QuasiOneSystem(model *cov);


bool isDummyInit(initfct Init);

bool isMtimesep(matrix_type type);
bool isMproj(matrix_type type);
bool isMdiag(matrix_type type);
bool isMiso(matrix_type type);

bool isDollar(model *cov);
bool isDollarProc(model *cov);
bool isAnyDollar(model *cov);
bool isNatsc(model *cov);



//bool isCov(model *cov);
typedef bool (*typusfct)(Types type); /* h, cov, result */ 
//bool is_any(typusfct t, defn *C);
//bool is_all(typusfct t, defn *C);
bool isRObject(int type);


int maxdim_ok(model *cov);
int maxdim_notok(model *cov);


bool equal_coordinate_system(isotropy_type iso1, isotropy_type iso2);
bool equal_coordinate_system(isotropy_type iso1, isotropy_type iso2, 
			     bool refined);
bool equal_coordinate_systems(system_type * sys1, system_type * sys2);


bool atleastSpecialised(isotropy_type iso, isotropy_type as);
bool isAngle(model *Aniso);
bool hasFullXdim(isotropy_type iso);


bool isDollar(model *cov);
bool isDollarProc(model *cov);
bool isAnyDollar(model *cov);

bool equalsNugget(int nr);
bool equalsNugget(model *cov);
bool equalsNuggetProc(model *cov);
bool isAnyNugget(model *cov);
bool isAnyNugget(int nr);

bool isNatsc(model *cov);
bool equalsBernoulliProcess(model *cov);

//bool isPlusMal(model *cov);


#endif

/* 
bool isIsotropic(isotropy_type iso);
bool isAnyIsotropic(isotropy_type iso);
bool isSpaceIsotropic(isotropy_type iso);
bool isZeroSpaceIsotropic(isotropy_type iso);
bool isVectorIsotropic(isotropy_type iso);
bool isCartesian(isotropy_type iso);
bool isSpherical(isotropy_type iso);
bool isCylinder(isotropy_type iso);
bool isEarth(isotropy_type iso);
bool isAnySpherical(isotropy_type iso);
bool isAnySphericalIso(isotropy_type iso);
bool isPrevModelI(defn *C);
bool isUnreduced(defn *C);


bool is_any(isotropy_type iso, defn *C);
bool is_all(isotropy_type iso, defn *C);
typedef bool (*isofct)(isotropy_type iso); // h, cov, result
bool is_any(isofct iso, defn *C);
bool is_all(isofct iso, defn *C);

*/
