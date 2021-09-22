
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of other simple function

Copyright (C) 2017 -- 2018 Martin Schlather
 

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


#include "def.h"
#include <Basic_utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>

#include "questions.h"
#include "RF.h"
#include "primitive.others.h"
#include "Coordinate_systems.h"
#include "operator.h"
#include "QMath.h"
//#include "AutoRandomFields.h"
//#include "shape.h"
//#include "Processes.h"


// Sigma(x) = diag>0 + A'xx'A
void kappa_Angle(int i, model *cov, int *nr, int *nc){
  int dim = OWNXDIM(0);
  *nr = i == ANGLE_DIAG ? dim : 1;
  *nc = i <= ANGLE_DIAG && (dim==3 || i != ANGLE_LATANGLE) &&  
    (dim==2 || i != ANGLE_RATIO)  ? 1 : -1;
}

void AngleMatrix(model *cov, double *A) {
  double c, s, 
     *diag = P(ANGLE_DIAG);
  if (GLOBAL.coords.anglemode == radians) {
    c = COS(P0(ANGLE_ANGLE));
    s = SIN(P0(ANGLE_ANGLE));
  } else {
    assert(GLOBAL.coords.anglemode == degree);
    c = COS(P0(ANGLE_ANGLE) * piD180);
    s = SIN(P0(ANGLE_ANGLE) * piD180); 
  }
  int 
    dim = OWNXDIM(0) ;
  assert(dim ==2 || dim == 3);
  
  if (dim == 2) {
    A[0] = c;
    A[1] = s;
    A[2] = -s;
    A[3] = c;
  } else {
    double pc, ps;
    if (GLOBAL.coords.anglemode == radians) {
      pc = COS(P0(ANGLE_LATANGLE));
      ps = SIN(P0(ANGLE_LATANGLE));
    } else {
      pc = COS(P0(ANGLE_LATANGLE) * piD180);
      ps = SIN(P0(ANGLE_LATANGLE) * piD180);
    }
    /*
    c -s 0    pc 0 -ps   c*pc  -s  -c*ps
    s  c 0  *  0 1  0  = s*pc   c  -s*ps
    0  0 1    ps 0  pc     ps   0   pc 
    */
    A[0] = c * pc;
    A[1] = s * pc;
    A[2] = ps;
    A[3] = -s;
    A[4] = c;
    A[5] = 0.0;
    A[6] = -c * ps;
    A[7] = -s * ps;
    A[8] = pc;
  }

  if (diag!= NULL) {
    int k,d,j;
    for (k=d=0; d<dim; d++) {
      for (j=0; j<dim; j++) {
	A[k++] *= diag[j];
      }
    }
  } else {
    double ratio = P0(ANGLE_RATIO);
    A[1] /= ratio;
    A[3] /= ratio;
  }  
}



void Angle(double *x, model *cov, double *v) { /// to do: extend to 3D!
  double A[9];
  int 
    dim = OWNXDIM(0) ;
  
  AngleMatrix(cov, A);
  Ax(A, x, dim, dim, v); 
}


void invAngle(double *x, model *cov, double *v) { /// to do: extend to 3D!
  double A[9],
    c = COS(P0(ANGLE_ANGLE)),
    s = SIN(P0(ANGLE_ANGLE)),
    *diag = P(ANGLE_DIAG);
  int d, 
      dim = OWNXDIM(0) ;
  
  bool inf = x[0] == RF_INF;
  for (d=1; d<dim; d++) inf &= x[d] == RF_INF;
  if (inf) {
    for (d=0; d<dim; v[d++] = RF_INF);
    return;
  }

  bool neginf = x[0] == RF_NEGINF;
  for (d=1; d<dim; d++) neginf &= x[d] == RF_NEGINF;
  if (neginf) {
    for (d=0; d<dim; v[d++] = RF_NEGINF);
    return;
  }
  
  if (dim == 2) {
    A[0] = c;
    A[1] = -s;
    A[2] = s;
    A[3] = c;
  } else {
    double pc = COS(P0(ANGLE_LATANGLE)),
      ps = SIN(P0(ANGLE_LATANGLE));

    A[0] = c * pc;
    A[1] = -s;
    A[2] = -c * ps;
    A[3] = s * pc;
    A[4] = c;
    A[5] = -s * ps;
    A[6] = ps;
    A[7] = 0.0;
    A[8] = pc;
  }

  if (diag!= NULL) {
    int  j, k;
    for (k=d=0; d<dim; d++) {
      for (j=0; j<dim; j++) {
	A[k++] /= diag[d];
      }
    }
  } else {
    double ratio = P0(ANGLE_RATIO);
    A[2] *= ratio;
    A[3] *= ratio;
  }

  Ax(A, x, dim, dim, v);
}

int checkAngle(model *cov){
  int dim = OWNXDIM(0);

  if (dim != 2 && dim != 3)
    SERR1("'%.50s' only works for 2 and 3 dimensions", NICK(cov));

  if (PisNULL(ANGLE_DIAG)) {
    if (PisNULL(ANGLE_RATIO)) {
      SERR2("either '%.50s' or '%.50s' must be given",
	    KNAME(ANGLE_RATIO), KNAME(ANGLE_DIAG))
    } else if (dim != 2) { 
      SERR1("'%.50s' may be given only if dim=2",  KNAME(ANGLE_RATIO))
    }
  } else {
    if (!PisNULL(ANGLE_RATIO)) 
      SERR2("'%.50s' and '%.50s' may not given at the same time",
	    KNAME(ANGLE_RATIO), KNAME(ANGLE_DIAG));
  }
  VDIM0 = dim;
  VDIM1 = 1;
  cov->mpp.maxheights[0] = RF_NA;
  cov->matrix_indep_of_x = true;
  RETURN_NOERROR;
}
 
void rangeAngle(model *cov, range_type* range){
  model *calling = cov->calling;
  assert(calling != NULL);
  range->min[ANGLE_ANGLE] = 0.0;
  range->max[ANGLE_ANGLE] = 
    calling->vdim[0] == SCALAR && isDollar(calling) && 
    calling->sub[0] != cov && // so a parameter of dollar
    SYSTEMS(CALLING) == 1 
    && isSymmetric(ISO(CALLING, 0)) ? PI : TWOPI;
  
  range->pmin[ANGLE_ANGLE] = 0.0;
  range->pmax[ANGLE_ANGLE] = range->max[ANGLE_ANGLE];
  range->openmin[ANGLE_ANGLE] = false;
  range->openmax[ANGLE_ANGLE] = true;

  range->min[ANGLE_LATANGLE] = 0.0;
  range->max[ANGLE_LATANGLE] = PI;
  range->pmin[ANGLE_LATANGLE] = 0.0;
  range->pmax[ANGLE_LATANGLE] = PI;
  range->openmin[ANGLE_LATANGLE] = false;
  range->openmax[ANGLE_LATANGLE] = true;

  range->min[ANGLE_RATIO] = 0;
  range->max[ANGLE_RATIO] = RF_INF;
  range->pmin[ANGLE_RATIO] = 1e-5;
  range->pmax[ANGLE_RATIO] = 1e5;
  range->openmin[ANGLE_RATIO] = false;
  range->openmax[ANGLE_RATIO] = true;

  range->min[ANGLE_DIAG] = 0;
  range->max[ANGLE_DIAG] = RF_INF;
  range->pmin[ANGLE_DIAG] = 1e-5;
  range->pmax[ANGLE_DIAG] = 1e5;
  range->openmin[ANGLE_DIAG] = false;
  range->openmax[ANGLE_DIAG] = true;
}
 


void idcoord(double *x, model *cov, double *v){
  int d,
    vdim = VDIM0;
  for (d=0; d<vdim; d++) v[d]=x[d];
}
int checkidcoord(model *cov){
  if (PREVISO(0) != OWNISO(0)) SERR("unequal iso's");
  VDIM0 = OWNTOTALXDIM;
  VDIM1 = 1;
  RETURN_NOERROR;
}



// obsolete??!!
#define NULL_TYPE 0
void NullModel(double VARIABLE_IS_NOT_USED *x, 
	       model VARIABLE_IS_NOT_USED *cov, 
	       double VARIABLE_IS_NOT_USED *v){ return; }
void logNullModel(double VARIABLE_IS_NOT_USED *x, 
		  model VARIABLE_IS_NOT_USED *cov, 
		  double VARIABLE_IS_NOT_USED *v, 
		  int VARIABLE_IS_NOT_USED *Sign){ return; }
Types TypeNullModel(Types required, model *cov,
		  isotropy_type VARIABLE_IS_NOT_USED i) {
  return TypeConsistency(isManifold(required) ? PREVTYPE(0) : required,
			 (Types) P0INT(NULL_TYPE));
}
void rangeNullModel(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[NULL_TYPE] = TcfType;
  range->max[NULL_TYPE] = OtherType;
  range->pmin[NULL_TYPE] = TcfType;
  range->pmax[NULL_TYPE] = OtherType;
  range->openmin[NULL_TYPE] = false;
  range->openmax[NULL_TYPE] = false;
}


int checkMath(model *cov){
  int i,
    variables = DefList[COVNR].kappas,
    err = NOERROR;
  if (variables > 2) kdefault(cov, variables-1, 1.0);
  if (isEarth(OWNISO(0))) {
    defn *C = DefList + COVNR;
    if ((C->cov==MathCos || C->cov==MathSin ||  C->cov==MathTan)) {
      //      printI(cov);
      //      printf("%.50s\n", TYPE_NAMES[DEFTYPE(0)]);
      //      APMI0(cov);
      SERR1("only radians as angular system allowed for '%.50s'.", NICK(cov));
    }
  }

  for (i=0; i<variables; i++) {
    model *sub = cov->kappasub[i];
    if (sub != NULL) {
      if (i >= 2) SERR("only numerics allowed");
      bool plus = DefList[SUBNR].cov == Mathplus ||
	DefList[SUBNR].check == checkplus ||
	DefList[SUBNR].cov == Mathminus
	;
      if ((err = CHECK_PASSTF(sub, plus ? OWNTYPE(0) : ShapeType, 
		       // auch falls cov = TrendType ist
		       1, cov->frame)) != NOERROR){
     //       if ((err = CHECK(sub, cov->tsdim, cov->xdimown, 
///			plus ? cov->typus : ShapeType, 
//			// auch falls cov = TrendType ist
//			OWNDOM(0), OWNISO(0),
//			1, cov->frame)) != NOERROR){
	RETURN_ERR(err);  
      }
      if  (sub->vdim[0] != 1 || sub->vdim[1] != 1)
	SERR("only scalar functions are allowed");
      setbackward(cov, sub); 
    } else if (PisNULL(i)) {
      if (i==0 || (DefList[COVNR].cov!=Mathplus && 
		   DefList[COVNR].cov!=Mathminus && 
		   DefList[COVNR].cov!=Mathbind
		   )) SERR("not enough arguments given")
      else break;
    }
  }

  cov->pref[Trendproc] = 5;
  cov->pref[Direct] = 1;
  RETURN_NOERROR;
}


void rangeMath(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  int i,
    variables = DefList[COVNR].kappas;
  
  set_maxdim(OWN, 0, OWNLOGDIM(0));
  assert(MAXDIM(OWN, 0) > 0);
  for (i=0; i<variables; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -1e5;
    range->pmax[i] = 1e5;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}

void Mathminus(double *x, model *cov, double *v){
  MATH_DEFAULT;
  double f = P0(MATH_FACTOR); 
  if (ISNAN(f)) f = 1.0;
  *v = ( (PisNULL(1) && cov->kappasub[1]==NULL) ? -w[0] : w[0]- w[1]) * f;
}

void Mathplus(double *x, model *cov, double *v){
  MATH_DEFAULT;
  double f = P0(MATH_FACTOR); 
  if (ISNAN(f)) f = 1.0;
  *v = ( (PisNULL(1) && cov->kappasub[1]==NULL)? w[0] : w[0] + w[1]) * f;
}

void Mathdiv(double *x, model *cov, double *v){
  MATH_DEFAULT;
   double f = P0(MATH_FACTOR); 
   if (ISNAN(f)) f = 1.0;
   *v =  (w[0] / w[1]) * f;
}

void Mathmult(double *x, model *cov, double *v){
  MATH_DEFAULT;
   double f = P0(MATH_FACTOR); 
   if (ISNAN(f)) f = 1.0;
   *v = (w[0] * w[1]) * f;
}
 

#define IS_IS 1
void Mathis(double *x, model *cov, double *v){
  int i,								
    variables = DefList[COVNR].kappas; 		
  double w[3];							
  for (i=0; i<variables; i++) {					
    model *sub = cov->kappasub[i];				
    if (sub != NULL) {						
      COV(x, sub, w + i);						
    } else {  
      w[i] = i == IS_IS ? P0INT(i) : P0(i);	      
    }					
  }									
  switch((int) w[IS_IS]) {
  case 0 : *v = (double) (FABS(w[0] - w[2]) <= GLOBAL.nugget.tol); break;
  case 1 : *v = (double) (FABS(w[0] - w[2]) > GLOBAL.nugget.tol); break;
  case 2 : *v = (double) (w[2] + GLOBAL.nugget.tol >= w[0]); break;
  case 3 : *v = (double) (w[2] + GLOBAL.nugget.tol > w[0]); break;
  case 4 : *v = (double) (w[0] + GLOBAL.nugget.tol >= w[2]); break;
  case 5 : *v = (double) (w[0] + GLOBAL.nugget.tol > w[2]); break;
  default : BUG;
  }
}

void rangeMathIs(model *cov, range_type *range){
  rangeMath(cov, range); 
  
  range->min[IS_IS] = 0;
  range->max[IS_IS] = 5;
  range->pmin[IS_IS] = 0;
  range->pmax[IS_IS] = 5;
  range->openmin[IS_IS] = false;
  range->openmax[IS_IS] = false;
}

void Mathbind(double *x, model *cov, double *v){
  int vdim = VDIM0;
  MATH_DEFAULT_0(vdim); 
  double f = P0(DefList[COVNR].kappas - 1); 
  f = ISNAN(f) ? 1.0 : f;
  for (i=0; i<vdim; i++) v[i] = w[i] * f; 
}

int check_bind(model *cov) {   // R.c
  int err;
  
  if ((err = checkMath(cov)) != NOERROR) RETURN_ERR(err);
  kdefault(cov, BIND_NCOL, 1);
  
  int
    ncol = P0INT(BIND_NCOL),
    variables = BIND_VARIABLES; 
 
  while (variables > 0 && cov->nrow[variables - 1] == 0 &&
	 cov->kappasub[variables - 1] == NULL) variables--;
  VDIM0 = variables / ncol;
  VDIM1 = ncol;
  if (VDIM0 * VDIM1 != variables) 
    SERR1("'%.50s' does not fit the number of components given", KNAME(BIND_NCOL));
  cov->ptwise_definite = pt_mismatch;

 RETURN_NOERROR;
}


bool allowedDbind(model *cov) {
  defn *C = DefList + COVNR;
  int i,
    kappas = C->kappas;
  for (i=0; i<kappas; i++) if (cov->kappasub[i] != NULL) break;
  if (i<kappas) {
    bool allowed = true,
      *D = cov->allowedD;
    for (int j=(int) FIRST_DOMAIN; j<=(int) LAST_DOMAINUSER; D[j++] = false);
    for (; i<kappas; i++) {
      model *k = cov->kappasub[i];
      if (k != NULL) {
	allowed &= allowedD(k);
	for (int j=(int) FIRST_DOMAIN; j<=(int) LAST_DOMAINUSER; j++)
	  D[j] |= k->allowedD[j];
      }
    }
    //printD(cov);printf("allowd.c=%d\n", allowed);
    return allowed;
  }
  //printf("Allowd.c=%d\n", allowedItrue(cov));
  return allowedItrue(cov);
}


bool allowedIbind(model *cov) {
  defn *C = DefList + COVNR;
  int i,
    kappas = C->kappas;
  for (i=0; i<kappas; i++) if (cov->kappasub[i] != NULL) break;
  if (i<kappas) {
    bool allowed = true,
      *I = cov->allowedI;
    for (int j=(int) FIRST_ISOUSER; j<=(int) LAST_ISOUSER; I[j++] = false);
    for (; i<kappas; i++) {
      model *k = cov->kappasub[i];
      if (k != NULL) {
	allowed |= allowedI(k);
	for (int j=(int) FIRST_ISOUSER; j<=(int) LAST_ISOUSER; j++)
	  I[j] &= k->allowedI[j];
      }
    }
    // printI(cov);printf("allowdI.c=%d\n", allowed);
    return allowed;
  }
  return allowedItrue(cov);
}


Types Typebind(Types required, model VARIABLE_IS_NOT_USED *cov,
		  isotropy_type VARIABLE_IS_NOT_USED iso) {
  return required;
  /*
  defn *C = DefList + COVNR;
  int i,
    kappas = C->kappas;
  for (i=0; i<kappas; i++) if (cov->kappasub[i] != NULL) break;
  if (i<kappas) {
    // if (is(required))
    // printf("bind: requ=%.50s %.50s\n",  TYPE_NAMES[required], ISO_NAMES[iso]);
    P MI0(cov);
    BUG;
    return BadType;
  } else {
    return required;
 }
  */
}

void Mathc(double VARIABLE_IS_NOT_USED *x, model *cov, double *v) {
  double f = P0(CONST_C); 
  *v =  ISNAN(f) ?  1.0 : f;
}

int check_c(model *cov){
  bool negdef = isnowNegDef(cov);
  bool tcf = isnowTcf(cov) || equalsSphericalIsotropic(OWNISO(0));

  if (negdef) {
    if (cov->calling == NULL) BUG;
    model *pp = cov->calling->calling;

    // the following defines a (positive) constant on the level of a trend 
    // to be a trend, and not a positive definite function
    // For spatially constant covariance functions, see RMconstant      
 
    if (pp == NULL || 
	(CALLINGNR == PLUS && !isnowNegDef(pp) && !isnowTrend(pp))){
      RETURN_ERR(ERRORFAILED); // by definition,
    }
  }
  if (cov->kappasub[0] != NULL) SERR("only numerics allowed");
  cov->ptwise_definite =  
    P0(CONST_C) > 0 ? pt_posdef : P0(CONST_C) < 0 ? pt_negdef : pt_zero;
  if (tcf) MEMCOPY(cov->pref, PREF_ALL, sizeof(pref_shorttype));

   GLOBAL.internal.warn_mathdef = GLOBAL.internal.warn_mathdef == False ? False
     : isNegDef(PREVTYPE(0)) ? True : Nan;
  
  RETURN_NOERROR;
}

void rangec(model *cov, range_type *range){
  rangeconstant(cov, range);
}

bool allowedIp(model *cov) {
  bool *I = cov->allowedI;
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  if (!PisNULL(PROJ_ISO)) {  
    isotropy_type iso = (isotropy_type) P0INT(PROJ_ISO);
    I[iso] = true;
    switch(iso) {
    case ISOTROPIC: case EARTH_ISOTROPIC: case SPHERICAL_ISOTROPIC :
      I[iso] = true;
      break;
    case DOUBLEISOTROPIC : 
    case VECTORISOTROPIC :
      ERR("'VECTORISOTROPY' not programmed yet");      
    case SYMMETRIC : case EARTH_SYMMETRIC: case SPHERICAL_SYMMETRIC :
      ERR2("Use '%.50s' within arbitrarty mathematical definitions (i.e. in '%.50s') or use the argument 'proj' within definite functions)",
	   CoordinateSystemOf(iso), NICK(cov));
      break;
    case CARTESIAN_COORD: case SPHERICAL_COORD:  case EARTH_COORD:
      I[iso] =  true;
      break;
    case GNOMONIC_PROJ  : case ORTHOGRAPHIC_PROJ : 
      ERR("Do not use projection in 'R,p', but use 'RMtrafo' instead.");
      break;
    case UNREDUCED : 
      I[CARTESIAN_COORD] = I[SPHERICAL_COORD] = I[EARTH_COORD] = true;
      break;
    default : ERR2("'%.50s' not allowed for '%.50s'", ISO_NAMES[iso],KNAME(PROJ_ISO));
    }
  } else {
    I[CARTESIAN_COORD] = I[SPHERICAL_COORD] = I[EARTH_COORD] = true;
  }
  return false;
}


void proj(double *x, model *cov, double *v){
  double f = P0(PROJ_FACTOR); 
  if (ISNAN(f)) f = 1.0;
  int p = P0INT(PROJ_PROJ);
  if (p >= 1) {
    *v = x[p - 1] * f;
    // PMIR(cov); BUG;
    // printf("nr=%d f=%10g p=%d x=%10g v=%10g\n", cov->zaehler, f, p, x[p-1], *v);
    return;
  }
  if (p == PROJ_TIME) *v = x[OWNTOTALXDIM - 1] * f;
  else if (p == PROJ_SPACE) {
    double dummy = 0.0;
    int spdim = OWNTOTALXDIM - 1;
    for (int d = 0; d<spdim; d++) dummy += x[d] * x[d];
    *v = SQRT(dummy) * f;
  } else BUG;
}

bool setproj(model *cov){
  isotropy_type
    iso = PisNULL(PROJ_ISO) ? PREVISO(0) : (isotropy_type) P0INT(PROJ_ISO);

  if (!isFixed(iso)) return false; // siehe unten -- nicht immer gebraucht!!
  
  domain_type dom = PREVDOM(0);
  Types type = PREVTYPE(0);
  bool fixed = isFixed(dom);
  
  if (equalsSpaceIsotropic(iso)) {
    NotProgrammedYet("");
    if (fixed) {
#if MAXSYSTEMS == 1
      set_system(OWN, 0, PREVLOGDIM(0), MAXINT, 2, type, dom, iso);
#else    
      BUG ?? return false;
      set_system(OWN, 0, totlogicaldim - 1, MAXINT, 1, type, dom, ISOTROPIC);
      set_system(OWN, 1, 1, 1, 1, OWNTYPE(1), OWNDOM(1), ISOTROPIC);
#endif    
    } else {
#if MAXSYSTEMS == 1
      set_iso(OWN, 0, iso); // rest from the definition
#else
      set_iso(OWN, 0, iso);
      not programmed yet
#endif    
    }
  } else if (isAnySpherical(iso)) {
    int s = 0;
    set_system(OWN, s, PREVLOGDIM(s), MAXINT,	 PREVLOGDIM(s), PREVTYPE(s),
	       PREVDOM(s), CoordinateSystemOf(PREVISO(s)));    
  } else if (equalsUnreduced(iso)) {
    if (!fixed) return false;
    int last = PREVLASTSYSTEM;
    for (int s=0; s<=last; s++) 
      set_system(OWN, s, PREVLOGDIM(s), MAXINT,	 PREVLOGDIM(s), PREVTYPE(s),
		 PREVDOM(s), CoordinateSystemOf(PREVISO(s)));
  } else {
     if (fixed) set_system(OWN, 0, PREVLOGDIM(0), MAXINT, 1, type, dom, iso);
    else set_iso(OWN, 0, iso);
  }
  return true;
}


int checkproj(model *cov){
   kdefault(cov, PROJ_FACTOR, 1.0);
   kdefault(cov, PROJ_PROJ, 1);

   //  if (isoown != iso && (iso != UNREDUCED || !equalsCoordinateSystem(isoown))) {
  //    SERR2("Offered system ('%.50s') does not match the required one ('%.50s')",
  //	  ISO_NAMES[isoown], ISO_NAMES[iso]);
  //  }

   if (P0INT(PROJ_PROJ) < 0 && !GetTime(cov)) 
    SERR2("'%.50s' or '%.50s' used in a context that is not spatio-temporal.",
	  PROJECTION_NAMES[0], PROJECTION_NAMES[1]);

  RETURN_NOERROR;
}


Types Typeproj(Types required, model *cov,
	       isotropy_type VARIABLE_IS_NOT_USED required_iso) {
  if (isBad(TypeConsistency(required, ShapeType)) && // ex-mathdeftype
      isBad(TypeConsistency(required, TrendType))) return BadType;
  return atleastSpecialised(OWNISO(0), required_iso) ? required : BadType;
}

void rangeproj(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[PROJ_PROJ] = -PROJECTIONS;
  range->max[PROJ_PROJ] = OWNXDIM(0);
  range->pmin[PROJ_PROJ] = 1;
  range->pmax[PROJ_PROJ] = OWNXDIM(0);
  range->openmin[PROJ_PROJ] = false;
  range->openmax[PROJ_PROJ] = false;

  range->min[PROJ_ISO] = FIRST_ISOUSER;
  range->max[PROJ_ISO] = LAST_ISOUSER;
  range->pmin[PROJ_ISO] = FIRST_ISOUSER;
  range->pmax[PROJ_ISO] = LAST_ISOUSER;
  range->openmin[PROJ_ISO] = false;
  range->openmax[PROJ_ISO] = false;

  range->min[PROJ_FACTOR] = RF_NEGINF;
  range->max[PROJ_FACTOR] = RF_INF;
  range->pmin[PROJ_FACTOR] = - 1e5;
  range->pmax[PROJ_FACTOR] = 1e5;
  range->openmin[PROJ_FACTOR] = true;
  range->openmax[PROJ_FACTOR] = true;
}


 

void declarefct(double  VARIABLE_IS_NOT_USED *x, model *cov, double *v){
  int vdimSq = VDIM0 * VDIM1;
  for (int i=0; i<vdimSq; v[i++] = 0.0);
}

void declarefctnonstat(double  VARIABLE_IS_NOT_USED *x,
		       double  VARIABLE_IS_NOT_USED *y, model *cov, double *v){
  int vdimSq = VDIM0 * VDIM1;
  for (int i=0; i<vdimSq; v[i++] = 0.0);
}


void rangedeclare(model *cov, range_type *range){
  int k = DefList[COVNR].kappas;
  for (int i=0; i<k; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = RF_NEGINF;
    range->pmax[i] = RF_INF;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}


int checkdeclare(model *cov) { 
  //  defn *C = DefList + COVNR;     
  int v = cov->calling->vdim[0];
  //  kappas = C->kappas;
  // for (int i=0; i<kappas; i++) if (!PisNULL(i)) P(i)[0]=RF_NA;
  if (v <= 0) v=1;
  VDIM0 = VDIM1 = v;
  RETURN_NOERROR;
}
