
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2001 -- 2017 Martin Schlather

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


#include "questions.h"
#include "Coordinate_systems.h"
#include "operator.h"
#include "Processes.h"

// domain
bool equalsXonly(domain_type dom) { return dom == XONLY; }
bool isXonly(system_type *s) { 
  assert({int _i_; 
      for (_i_=1; _i_<=LASTSYSTEM(s); _i_++) if (DOM(s, _i_) != DOM(s,0)) break;
      _i_ > LASTSYSTEM(s);});
  return equalsXonly(DOM(s, 0));
}
bool equalsKernel(domain_type dom) { return dom == KERNEL; }
bool isKernel(system_type *s) { 
  assert({int _i_; 
      for (_i_=1; _i_<=LASTSYSTEM(s); _i_++) if (DOM(s, _i_) != DOM(s,0)) break;
      _i_ > LASTSYSTEM(s);});
  return equalsKernel(DOM(s, 0));
} 

bool equalsDomMismatch(domain_type dom) { return dom == DOMAIN_MISMATCH; }
bool equalsPrevModelD(domain_type dom) { return dom == PREVMODEL_D;}
bool isPrevModelD(defn *C){
  assert(LASTSYSTEM(C->systems[0]) >= 0);
  bool is = equalsPrevModelD(DOM(C->systems[0], 0));
  assert(!is || LASTSYSTEM(C->systems[0]) == 0);
 return is;
}

bool equalsSubModelD(domain_type dom) { return dom == SUBMODEL_D;}
bool isSubModelD(defn *C){
  assert(LASTSYSTEM(C->systems[0]) >= 0);
  bool is = equalsSubModelD(DOM(C->systems[0], 0));
  assert(!is || LASTSYSTEM(C->systems[0]) == 0);
  return is;
}

bool equalsParamDepD(domain_type dom) { return dom == PARAMDEP_D;}
bool isParamDepD(defn *C){
  assert(LASTSYSTEM(C->systems[0]) >= 0);
  bool is = equalsParamDepD(DOM(C->systems[0], 0));
  assert(!is || LASTSYSTEM(C->systems[0]) == 0);
  return is;
}

bool equalsParamDepI(isotropy_type iso) { return iso == PARAMDEP_I;}
bool isParamDepI(defn *C){
  assert(LASTSYSTEM(C->systems[0]) >= 0);
  bool is = equalsParamDepI(ISO(C->systems[0], 0));
  assert(!is || LASTSYSTEM(C->systems[0]) == 0);
  return is;
}


// isotropy
bool equalsIsotropic(isotropy_type iso) { return iso == ISOTROPIC; }
bool isIsotropic(system_type *sys) {
  return LASTSYSTEM(sys) == 0 && equalsIsotropic(ISO(sys, 0)); }
bool isIsotropic(model *cov) { return isIsotropic(OWN); }

bool isIsotropicXonly(system_type *s) {
  return LASTSYSTEM(s) == 0 && equalsIsotropic(ISO(s, 0)) &&
    equalsXonly(DOM(s, 0)); }
bool isIsotropicXonly(model *cov) { return isIsotropicXonly(OWN);}

bool equalsSpaceIsotropic(isotropy_type iso) { return iso == DOUBLEISOTROPIC; }
bool equalsSpaceIsotropic(system_type *sys) {
#if MAXSYSTEMS == 1
  return equalsSpaceIsotropic(ISO(sys, 0));
#else   
  return LASTSYSTEM(sys) == 1 && equalsIsotropic(ISO(sys, 0)) && 
    equalsIsotropic(ISO(sys, 1)) && LOGDIM(sys, 1) == 1;
#endif  
}

bool isSpaceIsotropic(isotropy_type iso) {
  return equalsIsotropic(iso) || equalsSpaceIsotropic(iso);
}
bool isSpaceIsotropic(system_type *sys) {
  return isIsotropic(sys) || equalsSpaceIsotropic(sys);
}

bool isAnyIsotropic(isotropy_type iso) {
  return iso == ISOTROPIC || iso==EARTH_ISOTROPIC || iso==SPHERICAL_ISOTROPIC;}

bool equalsVectorIsotropic(isotropy_type iso) { return iso == VECTORISOTROPIC; }

bool isSymmetric(isotropy_type iso) { return iso <= SYMMETRIC; }
bool equalsSymmetric(isotropy_type iso) { return iso == SYMMETRIC; }

bool equalsCartCoord(isotropy_type iso) { return iso == CARTESIAN_COORD; }

bool isCartesian(isotropy_type iso) { return iso <= LAST_CARTESIAN; }
bool isCartesian(system_type *sys, int s) {
  return LASTSYSTEM(sys) >= s && isCartesian(ISO(sys, s));}
bool isCartesian(system_type *sys) {
  int L = LASTSYSTEM(sys);
  for (int s=0; s<=L; s++)
    if (!isCartesian(ISO(sys, s))) return false;
  return true;
}

bool equalsAnySymmetric(isotropy_type iso) {
  return equalsSymmetric(iso) || equalsEarthSymmetric(iso) ||
    equalsSphericalSymmetric(iso);
}


bool isGenuineCartesian(isotropy_type iso) { return iso <= CARTESIAN_COORD; }


bool equalsGnomonic(isotropy_type iso) { return iso == GNOMONIC_PROJ; }
bool equalsOrthographic(isotropy_type iso) { return iso == ORTHOGRAPHIC_PROJ; }
bool isEarthProjection(isotropy_type iso) {
  return equalsGnomonic(iso) || equalsOrthographic(iso);
}


bool equalsSphericalIsotropic(isotropy_type iso) {
  return iso == SPHERICAL_ISOTROPIC;
}
bool equalsSphericalSymmetric(isotropy_type iso) {
  return iso == SPHERICAL_SYMMETRIC;
}
bool equalsSphericalCoord(isotropy_type iso) {
  return iso == SPHERICAL_COORD;
}
bool isSphericalSymmetric(isotropy_type iso) {
  return equalsSphericalSymmetric(iso) || equalsSphericalIsotropic(iso);
}
bool isSpherical(isotropy_type iso) {
  return iso >= FIRST_SPHERICAL && iso <= LAST_TRUESPHERICAL; }


bool equalsEarthIsotropic(isotropy_type iso) {
  return iso == EARTH_ISOTROPIC;
}
bool equalsEarthSymmetric(isotropy_type iso) {
  return iso == EARTH_SYMMETRIC;
}
bool equalsEarthCoord(isotropy_type iso) {
  return iso == EARTH_COORD;
}
bool isEarthSymmetric(isotropy_type iso) {
  return equalsEarthSymmetric(iso) || equalsEarthIsotropic(iso);
}
bool isEarth(isotropy_type iso) {
  return iso >= FIRST_EARTH && iso <= LAST_EARTH; }




bool isAnySpherical(isotropy_type iso) {
  return iso >= FIRST_ANYSPHERICAL && iso <= LAST_ANYSPHERICAL; }
bool isAnySphericalIso(isotropy_type iso) {
  return iso == SPHERICAL_ISOTROPIC || iso == EARTH_ISOTROPIC; }



bool isLogCart(isotropy_type VARIABLE_IS_NOT_USED iso) {
  return false; // iso >= FIRST_LOGCART && iso <= LAST_LOGCART;
}
bool isLogCart(system_type VARIABLE_IS_NOT_USED *sys, int VARIABLE_IS_NOT_USED s) {
  return false; // LASTSYSTEM(sys) >= s && isLogCart(ISO(sys, s));
}
bool equalsLogCartSymmetric(isotropy_type VARIABLE_IS_NOT_USED iso) {
  return false; //iso == LOGCART_SYMMETRIC;
}
bool equalsLogCartCoord(isotropy_type VARIABLE_IS_NOT_USED  iso) {
  return false; //iso == LOGCART_COORDS;
}


bool isUnreducedLogCart(isotropy_type iso) {
  return // equalsLogCartSymmetric(iso) ||
    equalsLogCartCoord(iso); 
}
bool isUnreducedEarth(isotropy_type iso) {
  return // equalsEarthSymmetric(iso) ||
    equalsEarthCoord(iso); 
}
bool isUnreducedSpherical(isotropy_type iso) {
  return //equalsSphericalSymmetric(iso) ||
    equalsSphericalCoord(iso); 
}
bool isUnreducedCartesian(isotropy_type iso) {
  return isCartesian(iso) && iso > LAST_REDUCEDXDIM_CART; }
bool isUnreduced(isotropy_type iso) { // neue def, siehe hasFullxdim
  return isUnreducedCartesian(iso) ||
    isUnreducedEarth(iso) ||
    isUnreducedSpherical(iso) ||
    isUnreducedLogCart(iso); 
} 


bool equalsPrevModelI(isotropy_type iso) { return iso == PREVMODEL_I;}
bool isPrevModelI(defn *C){
  assert(LASTSYSTEM(C->systems[0]) >= 0);
  bool is = equalsPrevModelI(ISO(C->systems[0], 0));
  assert(!is || LASTSYSTEM(C->systems[0]) == 0);
 return is;
}

bool equalsSubModelI(isotropy_type iso) { return iso == SUBMODEL_I;}
bool isSubModelI(defn *C){
  assert(LASTSYSTEM(C->systems[0]) >= 0);
  bool is = equalsSubModelI(ISO(C->systems[0], 0));
  assert(!is || LASTSYSTEM(C->systems[0]) == 0);
  return is;
}


bool equalsIsoMismatch(isotropy_type iso) { return iso == ISO_MISMATCH; }
bool equalsUnreduced(isotropy_type iso) { return iso == UNREDUCED; }


bool equalsUnreduced(defn *C){
  bool is = equalsUnreduced(ISO(C->systems[0], 0));
  return is;
}

bool equalsCoordinateSystem(isotropy_type iso) {
  return iso==CARTESIAN_COORD || iso==SPHERICAL_COORD ||  iso==EARTH_COORD; }


bool isFixed(isotropy_type iso) {
  return iso <= LAST_ISOUSER;
}

bool isFixed(domain_type dom) {
  return dom <= LAST_DOMAINUSER;
}

bool isCylinder(isotropy_type iso) {
   return iso == CYLINDER_COORD;
}



// monotonicity
bool isMonotone(monotone_type  monotone) {
  return monotone >= MONOTONE && monotone <= COMPLETELY_MON; }
bool isMonotone(model *cov) { return isMonotone(cov->monotone); }

bool isCompletelyMonotone(monotone_type monotone) { return monotone == COMPLETELY_MON; }

bool isNormalMixture(monotone_type monotone) {
  return monotone == NORMAL_MIXTURE || isCompletelyMonotone(monotone);}
bool isNormalMixture(model *cov) { return isNormalMixture(cov->monotone); }

bool isBernstein(monotone_type monotone) { return monotone == BERNSTEIN; }
bool isBernstein(model *cov) { return isBernstein(cov->monotone); }

bool isGneiting(monotone_type monotone) {
  return monotone == GNEITING_MON || isCompletelyMonotone(monotone);}
bool isGneiting(model *cov) { return isGneiting(cov->monotone); }

bool isManifold(defn *C) {
  bool is = isManifold(SYSTYPE(C->systems[0], 0));
  assert(!is || C->variants == 1) ;
  return is;
}

bool isTrend(model *cov) {
  return COVNR == TREND_PROC || isMathDef(DefList + COVNR) || COVNR==COVARIATE;
}

bool isMathDef(defn *C) {
  return LASTSYSTEM(C->systems[0]) == 0 &&
    SYSTYPE(C->systems[0], 0) == MathDefType;
}



bool isDefCL(typefct isTypus, model *cov, bool OneSystem) {
  // CL ^= DefList
  defn *C = DefList + COVNR;
  //  printf("isDefCL %.50s %d %d\n", C->name, cov->variant, C->variants);
  assert(cov->variant == UNSET ||
	 (cov->variant >= 0 && cov->variant <= C->variants));
  system_type *sys = C->systems[cov->variant == UNSET //|| C->variants 
				? 0 : cov->variant];
  //  if (LASTi((sys)[0]) < 0) { printf(" %d %d\n",cov->variant, (LASTi((sys)[0])));TREE0(cov);PMI0(cov); crash();}
  int n = SYSTEMS(sys);
  if (OneSystem && n != 1) return false;
  if (C->TypeFct != NULL) {
    //    printf("call of TypeFct in isDefCL by is%.50s(%.50s)\n", TYPE_NAMES[type], NAME(cov));
    return false;
  }
  if (!isTypus(SYSTYPE(sys, 0))) return false;
  for (int s=1; s<n; s++) if (!isSameAsPrev(SYSTYPE(sys, s))) return false;
  return true;
}

bool isNow(typefct isTypus, model *cov, bool OneSystem) {
  // Now ^= cov
  int n = OWNSYSTEMS;
  if ((OneSystem && n != 1) || !isTypus(OWNTYPE(0))) return false;
  for (int s=1; s<n; s++) if (!isSameAsPrev(OWNTYPE(s))) return false;
 return true;
}


//if (false) printf("Qu: "#typus" '%.50s'  %d\n", TYPE_NAMES[type], type == cond);

#define QuestionGenNowNoIsNow(typus, cond, condFrame, OneSystem)	\
  bool is##typus(Types type) {						\
    return type == cond;						\
  }									\
									\
  bool equals##typus(Types type) {					\
    return type == typus##Type;						\
  }									\
									\
  bool equalsnow##typus(model *cov) {					\
    return isNow(equals##typus, cov, OneSystem);			\
  } 									\
									\
  bool has##typus##Frame(model *cov) {					\
    return cov->frame /* naechste Zeile weiter wegen findall */		\
      == typus##Type;							\
  } 									\
									\
  bool hasAny##typus##Frame(model *cov) {				\
    Types type = cov->frame;						\
    return type == condFrame;						\
  }


#define QuestionGenNow(typus, cond, condFrame, OneSystem)		\
  QuestionGenNowNoIsNow(typus, cond, condFrame, OneSystem)		\
  bool isnow##typus(model *cov) {					\
    return isNow(is##typus, cov, OneSystem);				\
  } 									
  

#define QuestionGen(typus, cond, condFrame, OneSystem)			\
  QuestionGenNow(typus, cond, condFrame, OneSystem)			\
  bool is##typus(model *cov) {						\
    return cov != NULL && isDefCL(is##typus, cov, OneSystem);		\
  }									
								       
#define Question(typus, OneSystem)					\
  QuestionGenNow(typus,typus##Type, typus##Type, OneSystem)		\
  bool is##typus(model *cov) {						\
    return isDefCL(is##typus, cov, OneSystem);				\
  }									
 

#define QuestionIsNow(typus, cond, OneSystem) QuestionGen(typus,cond,cond,OneSystem)
#define QuestionNow(typus, cond, OneSys) QuestionGenNow(typus,cond,cond,OneSys)


bool isMaxStable(Types type) {
  return type == BrMethodType || type == SmithType || type == SchlatherType; }
bool isnowMaxStable(model *cov) { return isNow(isMaxStable, cov, true); }
bool hasMaxStableFrame(model *cov) { return isMaxStable(cov->frame); }

QuestionIsNow(Tcf, TcfType, true)
QuestionNow(PosDef, PosDefType || isTcf(type), false)
QuestionNow(Variogram, VariogramType || isPosDef(type), false)
QuestionNow(NegDef, NegDefType || isVariogram(type), false)
QuestionIsNow(Process, ProcessType || type == GaussMethodType ||
	      type == PoissonType || 
	      isMaxStable(type) || type == NormedProcessType, true)

Question(GaussMethod, true) 
Question(NormedProcess, true)
Question(BrMethod, true)//
Question(PointShape, true)
QuestionIsNow(Random, RandomType || isProcess(type), true) 
QuestionGen(Shape,
	    ShapeType || isNegDef(type) || isMathDef(type),
	    RandomType || isMaxStable(type) || type == PoissonType,
	    true)
QuestionNow(Trend, TrendType || isMathDef(type) ||  type == ShapeType , false)
Question(Interface, true)
Question(RandomOrShape, true) // vermutlich obsolete !!
Question(Manifold, true) 
Question(MathDef, true)

Question(Other, true)
Question(Bad, true)
Question(SameAsPrev, false) // (model *cov) is senseless
// QuestionIsNow(Any, AnyType || true, true)

Question(Smith, true)//
Question(Schlather, true)
QuestionGen(Poisson, PoissonType, PoissonType || type == PoissonGaussType, true)
Question(PoissonGauss, true) 

Question(Likelihood, true)
QuestionNow(Evaluation, EvaluationType || type == LikelihoodType, true)


// R Types

bool isRObject(int type) {
  if (type == CLOSXP) BUG;
  return type == LANGSXP || type == VECSXP || type == ENVSXP;
}



// Parameter types
bool allowsRandomParam(model *cov, int i) {
  Types type = DefList[COVNR].kappaParamType[i];
  return isRandom(type) || isRandomOrShape(type);
}

bool allowsShapeParam(model *cov, int i) {
  Types type = DefList[COVNR].kappaParamType[i];
  return isShape(type) || isRandomOrShape(type);
}


bool isnowTrendParam(model *cov, int i) {
  //  PMI(cov);
  //  assert(cov->initialised);
  return isnowTrend(cov) && (SortOf(cov, i, 0, 0, original) == TRENDPARAM);
}


bool anyVariant(systemfct sf, defn *C) {
  int i;
  for (i=0; i < C->variants; i++) if (sf(C->systems[i])) return i;
  return false;
}
bool anyVariant(typefct tf, defn *C) {
  int i;
  for (i=0; i < C->variants; i++) {
    if (tf(SYSTYPE(C->systems[i], 0))) {
      int s,
	last = LASTSYSTEM(C->systems[i]);
      for (s=1; s<=last; s++) 
	if (SYSTYPE(C->systems[i], s) != SameAsPrevType) break;
      if (s > last) return true;
    }
  }
  return false;
}

/*
bool everyVariant(systemfct sf, defn *C, int s) {
  int i;
  for (i=0; i < C->variants; i++) if (!sf(C->systems[i])) return false;
  return true;
}
*/


bool everyCoord(systemfct sysfct, model *cov) {
  int 
    n = OWNLASTSYSTEM;
  for (int s=0; s<=n; s++) if (!sysfct(OWN + s)) return false;
  return true;
}

bool HaveSameSystems(system_type * fst, system_type * snd) { 
  int n = LASTSYSTEM(fst);
  if (n != LASTSYSTEM(snd)) return false;
  for (int s=0; s<n; s++) {
    if (CoordinateSystemOf(ISO(fst, s)) != CoordinateSystemOf(ISO(snd, s)))
      return false;
  }
  return true;
}


int total_logicaldim(system_type * sys) {
  assert(sys[0].last >= 0);
  int tot =LOGDIM(sys, 0),
    last = LASTSYSTEM(sys);
  for (int s=1; s<=last; s++) tot += LOGDIM(sys, s);
return tot;
}

int maxdim_notok(model *cov) {
  int last = OWNLASTSYSTEM;
  for (int s=0; s<=last; s++)
    if (OWNMAXDIM(s) >= 0 &&  OWNMAXDIM(s) < OWNLOGDIM(s)) return s;
  return -1;
}



bool QuasiOneSystem(model *cov) {
  int n = OWNSYSTEMS;
  if (n == 1) return true;
  //  isotropy_type iso = OWNISO(0);
  domain_type dom = OWNDOM(0);
  if (!equalsnowSameAsPrev(cov)) return false;
  for (int s=1; s<n; s++) {
    if (dom != OWNDOM(s)) return false;
  }
  return true;
}




/*
#define NOTMATCHING -1
int matchingVariant(systems_type sys, model *cov, int depth) {
   int i;
  for (i=0; i<C->variants; i++) {
    if (Type Consistency(requiredtype, SYSTYPE(C->systems[i], xxx))
	&& (depth <= 0 || atleast Specialised(OWNISO(0), ISO(C->systems[i], xxx)))
	) {
      //         printf("ok ? %d \n", i +1);
      return i;
    }
  }
  return NOTMATCHING;
}

int SysConsistency(systems_type sys, model *cov, int depth) {
  assert(cov != NULL);
  int last = SYSTEMS(sys);
  defn *C = DefList + COVNR; // nicht gatternr
  if (isManifold(C)) {
    assert(C->TypeFct != NULL);    
    // printf(" undefined!\n");
    int tc = C->TypeFct(sys, cov, depth);
    //     printf("tc=%d %.50s\n", tc, NAME(cov));
    return tc;
  }  
  int n = OWNSYSYEMS;
  for (int i=0; i<C->variants; i++) {
    book notok = true;
    int j=0;
    Types type = TYPE(C->systems[i], 0);
    for (int s=0; s<n; s++) {
      if (TYPE(C->systems[i], s) != SameAsPrevType) {
	Types newtype = TYPE(C->systems[i], s);
	j++;
	if (newtype != type) {
	  while (j <= last && requiredtype[j] == requiredtype[j-1]) j++;
	} 
	if (j > last) return NOTMATCHING;
	type = newtype;
      }
      notok = !Type Consistency(requiredtype[j], type)
	||  (depth > 0 && !atleast Specialised(OWNISO(s), ISO(C->systems[i], zzz), s)));
      if (notok) break;       
    }
    if (!notok) return i;
  }
  return NOTMATCHING;
}


int SysConsistency(Types requiredtype, model *cov, int depth) {
  system_type sys;
  RESETLAST(sys);
  SYSTYPE(sys, 0) = requiredtype;
  return Type Consistency(sys, 1, cov, depth);
}

*/




bool isMiso(matrix_type type) { 
  return type == TypeMiso;
}


bool isMproj(matrix_type type) { 
  assert(TypeMtimesepproj == 2);
 return type == TypeMproj || type <= TypeMtimesepproj;
}


bool isMdiag(matrix_type type) { 
  return type <= TypeMdiag;
}

bool isMtimesep(matrix_type type) { 
  assert(TypeMtimesep == 3);
  return type <= TypeMtimesep;
}


bool isverysimple(model *cov) {
  assert(cov != NULL);
  defn *C = DefList + COVNR; // kein gatternr, da isotrop
  int i, k, total,
    kappas = C->kappas;
  for (k=0; k<kappas; k++) {
    if (cov->kappasub[k] != NULL) return false;
    total = cov->ncol[k] * cov->nrow[k];
    if (C->kappatype[k] == REALSXP) {
      for (i=0; i<total; i++) 
	if (ISNAN(P(k)[i])) return false;
    } else  if (C->kappatype[k] == INTSXP) {
      for (i=0; i<total; i++) if (P(k)[i] == NA_INTEGER) return false;
    } else return false;
  }
  return true;
}


bool isDollar(model *cov) {
  assert(cov != NULL);
  int nr=COVNR;
  return nr >= DOLLAR && nr <= LASTDOLLAR;
}


bool isDollarProc(model *cov) {
  assert(cov != NULL);
  int nr=COVNR;
  return nr == DOLLAR_PROC;
}

bool isAnyDollar(model *cov) {
  assert(cov != NULL);
  int nr=COVNR;
  return (nr >= DOLLAR && nr <= LASTDOLLAR) || nr == DOLLAR_PROC;
}


bool equalsNugget(int nr) {
  return nr == NUGGET;
}
bool equalsNugget(model *cov) {
  assert(cov != NULL);
  return equalsNugget(COVNR);
}

bool equalsNuggetProc(model *cov) {
  assert(cov != NULL);
  int nr=COVNR;
  return nr ==  NUGGET_USER || nr == NUGGET_INTERN;
}


bool isAnyNugget(int nr) {
  return equalsNugget(nr) || nr ==  NUGGET_USER || nr == NUGGET_INTERN;
}

bool isAnyNugget(model *cov) {
  assert(cov != NULL);
  return equalsNugget(cov) || equalsNuggetProc(cov);
}

bool isNatsc(model *cov) {
  assert(cov != NULL);
  int nr=COVNR;
  return nr == NATSC_INTERN || nr == NATSC_USER;
}



bool equalsBernoulliProcess(model *cov) {
  int nr = COVNR;
  return nr == BINARYPROC;
}



bool isCov(model *cov) {
  int nr = COVNR;
  return nr == COVFCTN || nr == COVMATRIX ;
}


int maxdim_ok(model *cov) {
  int last = OWNLASTSYSTEM;
  for (int s=0; s<=last; s++)
    if (OWNMAXDIM(s) >= 0 &&  OWNMAXDIM(s) < OWNLOGDIM(s)) return -s;
  return 1;
}




bool equal_coordinate_system(isotropy_type iso1, isotropy_type iso2) {
  return (isCartesian(iso1) && isCartesian(iso2))
    || (isAnySpherical(iso1) && isAnySpherical(iso2))
    // || (isCylinder(iso1) && isCylinder(iso2)) 
    || hasFullXdim(iso1) || hasFullXdim(iso2);
}

bool equal_coordinate_system(isotropy_type iso1, isotropy_type iso2, 
			     bool refined) {
  if (!refined) return equal_coordinate_system(iso1, iso2);

  return (isCartesian(iso1) && isCartesian(iso2))
    || (isSpherical(iso1) && isSpherical(iso2))
    || (isEarth(iso1) && isEarth(iso2))
    // || (isCylinder(iso1) && isCylinder(iso2)) 
    || (equalsUnreduced(iso1) && equalsUnreduced(iso2));
}


bool equal_coordinate_systems(system_type * sys1, system_type * sys2, 
			     bool refined) {
  int last = sys1->last;
  if (last == UNSET) BUG;
  if (last != sys2->last) return false;
  for (int s=0; s<last; s++) {
    if (!equal_coordinate_system(ISO(sys1, s), ISO(sys2, s), refined))
      return false;
  }
  return true;
}


bool equal_coordinate_systems(system_type * sys1, system_type * sys2) {
  return equal_coordinate_systems(sys1, sys2, false);
}


bool atleastSpecialised(isotropy_type iso, isotropy_type as) {
  if (iso == as) return true; // e.g. MISMATCH; PARAM_DEP, etc.
  if (equalsPrevModelI(iso) || equalsSubModelI(iso)) return true;
  if (equalsUnreduced(as)) return true;  
  if (equalsUnreduced(iso)) return isUnreduced(as);
 
  //  printf("atleast %.50s as %.50s; isprevmodeli=%d\n", ISO_NAMES[iso], ISO_NAMES[as],	 isPrevModelI(iso));

  bool iso_more_specialised = iso < as; 
  if (isCartesian(as)) // komme von cartesian; Umwandlung zu sphere/earth
    //                    derzeit nicht moeglich
    return equalsVectorIsotropic(as) ? false : iso_more_specialised; 
  if (isSpherical(as)) { // 
    return isSpherical(iso) && iso_more_specialised; // keine umwandlung von
    //                                spherical zu irgendwas anderem derzeit 
  }
  if (isEarth(as)) {
    if (isEarth(iso)) return iso_more_specialised;
    if (isSpherical(iso)) return iso <= as - EARTH_COORD + EARTH_ISOTROPIC;
    if (isCartesian(iso)) {
      return !(equalsEarthIsotropic(as) ||  // FEHLT HIER WAS????
	       (equalsCartCoord(iso) && isEarthSymmetric(as)));
    }
    //  printf("%.50s %.50s\n", ISO_NAMES[iso], ISO_NAMES[as]);
    //  crash();
    BUG; 
  }
  BUG;
}



Types TypeConsistency(Types requiredtype, Types deliveredtype) {
  //printf("tc: %.50s %.50s\n", TYPE_NAMES[requiredtype], TYPE_NAMES[deliveredtype]);
  if (isBad(deliveredtype)) BUG;
  if (isManifold(deliveredtype)) BUG;

  switch(requiredtype) {
  case TcfType:      return isTcf(deliveredtype) ? deliveredtype : BadType;
  case PosDefType:   return isPosDef(deliveredtype) ? deliveredtype : BadType;
  case VariogramType:return isVariogram(deliveredtype) ? deliveredtype :BadType;
  case NegDefType:   return isNegDef(deliveredtype) ? deliveredtype : BadType;
    
  case PointShapeType: return isPointShape(deliveredtype) ? deliveredtype
      : BadType;
  case ShapeType : return isShape(deliveredtype) ? deliveredtype : BadType;
  case TrendType : return isTrend(deliveredtype) ? deliveredtype : BadType;
    // case RandomOrShapeType : BUG; // sollte nur bei den Argumenten auftauchen
    // case ManifoldType: BUG;

  case ProcessType : return isProcess(deliveredtype) || 
		            isTrend(deliveredtype) ? deliveredtype : BadType;
  case GaussMethodType: return isGaussMethod(deliveredtype) ? deliveredtype
                            : BadType;
  case NormedProcessType: return isNormedProcess(deliveredtype) ? deliveredtype
                            : BadType;
  case BrMethodType: return isBrMethod(deliveredtype) ? deliveredtype : BadType;
  case SmithType : return isSmith(deliveredtype) ? deliveredtype : BadType;
  case SchlatherType: return isSchlather(deliveredtype) ?deliveredtype :BadType;
  case PoissonType : return isPoisson(deliveredtype) ? deliveredtype : BadType;
  case PoissonGaussType : return isPoissonGauss(deliveredtype) ? deliveredtype
      : BadType;
    
  case RandomType :  return isRandom(deliveredtype) ? deliveredtype : BadType;
  case InterfaceType:return isInterface(deliveredtype) ? deliveredtype :BadType;
    // case MathDefType: BUG;
  case OtherType :   return isOther(deliveredtype) ? deliveredtype : BadType;
  default : // printf("%.50s\n", TYPE_NAMES[requiredtype]);
    BUG;
    // rest sollte nur bei den Argumenten auftauchen
  }
  BUG;
  return BadType;
}



// Types TypeConsistencyCL(Types requiredtype, model *cov) bis Version 3.2.0.23

//#define TypeDebug 1
Types TypeConsistency(Types requiredtype, model *cov, isotropy_type requirediso) {
  assert(!isManifold(requiredtype));
  assert(!isBad(requiredtype));
  assert(!equalsParamDepI(requirediso));
  assert(!equalsSubModelI(requirediso));
  assert(!equalsIsoMismatch(requirediso));
  assert(requirediso != KEEPCOPY_ISO);
  assert(!equalsIsoMismatch(requirediso));

  //  PMI0(cov);
#ifdef TypeDebug
  PRINTF(">>>! typeconsist: %.50s requ=('%.50s' %.50s) %.50s\n", NAME(cov), TYPE_NAMES[requiredtype], ISO_NAMES[requirediso], ISO_NAMES[OWNISO(0)]);
 #endif
 
  assert(cov != NULL);
  defn *C = DefList + COVNR; // nicht gatternr
  if (C->TypeFct != NULL) {
#ifdef TypeDebug
    PRINTF("TypeFct != NULL, %d ", atleastSpecialised(OWNISO(0), requirediso));
#endif
    if (!atleastSpecialised(OWNISO(0), requirediso)) return BadType;

      
  //   PMI(cov); 
    //
#ifdef TypeDebug
    PRINTF("MANIFOLD %.50s requ:%.50s, %.50s (own=%.50s)\n", NAME(cov), TYPE_NAMES[requiredtype], ISO_NAMES[requirediso], ISO_NAMES[OWNISO(0)]);
#endif
    Types tc = C->TypeFct(requiredtype, cov,  OWNISO(0));
#ifdef TypeDebug    
    PRINTF("tc=%.50s\n", TYPE_NAMES[tc]);
    if (isBad(tc)) {
      PMI(cov); //
      //  BUG;
    }
#endif
    
    if (!isBad(tc) && isnowManifold(cov)) set_type(OWN, 0, (Types) tc);
    //    PMI(cov);
    //    printf("tc returned =%d\n", tc);
    return tc;
  }

 
  int V = cov->variant;

  if (V == UNSET) {
    //    printf("V unset %.50s\n", NAME(cov));
    int minV = 0,
    maxV = C->variants - 1;
    for (V=minV; V<=maxV; V++) {
#ifdef TypeDebug
      PRINTF("V for=%d(%d %d) %d %.50s atleast=%d\n", V, minV, maxV, TypeConsistency(requiredtype, SYSTYPE(C->systems[V], 0)), ISO_NAMES[ISO(C->systems[V],0)], atleastSpecialised(ISO(C->systems[V], 0), requirediso));
#endif

      Types type = SYSTYPE(C->systems[V], 0);
      if (!isBad(TypeConsistency(requiredtype, type)) &&
	  atleastSpecialised(ISO(C->systems[V], 0), requirediso)) {
	return type;
      }
    }

#ifdef TypeDebug
    PRINTF("bad type\n");
#endif
    return BadType;
  }
  
  assert(V >= 0 && V < MAXVARIANTS);
#ifdef TypeDebug
  //  PMI(cov->calling);
  PRINTF("U=%d %d\n", V, cov->variant);
#endif
  
  Types type = SYSTYPE(C->systems[V], 0);
  isotropy_type iso = ISO(C->systems[V], 0);
  if (C->setDI != NULL) {
#ifdef TypeDebug
    PRINTF("setDI !\n");
 #endif
    iso = OWNISO(0);
  }
  bool consist =  !isBad(TypeConsistency(requiredtype, type)) &&
    atleastSpecialised(iso, requirediso);
  
#ifdef TypeDebug
  PRINTF("W=%d type=%.50s requ=%.50s consist=%d (%d!=%d & %d) && %.50s requ.iso=%.50s\n", V,  TYPE_NAMES[type],TYPE_NAMES[requiredtype], consist, TypeConsistency(requiredtype, type), BadType, atleastSpecialised(iso, requirediso), ISO_NAMES[iso], ISO_NAMES[requirediso]);
#endif

  // if (!consist && equalsTrend(PREVTYPE(0))) {PMI(cov); BUG;}
  
  return consist ? type : BadType;
}

bool isAngle(model *Aniso) {
  return Aniso!=NULL && DefList[MODELNR(Aniso)].check == checkAngle;
}

bool hasFullXdim(isotropy_type iso) {
  return equalsCoordinateSystem(iso) || equalsAnySymmetric(iso) ||
    equalsVectorIsotropic(iso);
}

bool hasFullCartesian(isotropy_type iso) {
  // printf("is=%d %d %d\n", isCartesian(iso), iso, LAST_REDUCED_CART);
  return isCartesian(iso) && iso > LAST_REDUCEDXDIM_CART;
}

/*
bool isPlusMal(model *cov) {
  BUG;
  defn *C = DefList + COVNR;
  return C->check == checkplus ;
}
*/
