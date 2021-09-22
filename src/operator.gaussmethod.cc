/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives of hypermodels

 Copyright (C) 2005 -- 2017 Martin Schlather
               2015-2017 Olga Moreva (cutoff, modified)
 
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
#include <R_ext/Applic.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "operator.h"
#include "Processes.h"
#include "cubicsolver.h"


typedef int (*set_local_q_type) (model *next, double a, double d);



void tbm3(double *x, model *cov, double *v, double tbmdim){
  model *next = cov->sub[TBMOP_COV];
  int i,
    vdim = VDIM0,
    vdimsq = vdim * vdim;
  assert(VDIM0 == VDIM1);
  TALLOC_XX1(v1, vdimsq);
  COV(x, next, v); // x has dim 2, if turning planes
  if (x[0] != 0.0) {
    Abl1(x, next, v1);
    for (i=0; i<vdimsq; i++) v[i] += *x * v1[i] / tbmdim;
  }
  END_TALLOC_XX1;
}

typedef struct TBM2_integr{
  model *cov;
  double *x;
} TBM2_integr;

void TBM2NumIntegrFct(double *u,  int n, void *ex) {
  int i;
  double z[2];
  TBM2_integr *info;
  info = (TBM2_integr*) ex;
  model *cov=info->cov;
  double *x = info->x;
    
  for (i=0; i<n; i++) {
    z[0] = x[0] * SQRT(1.0 - u[i] * u[i]);
    tbm3(z, cov, u + i, 1.0);
  }
}

void tbm2num(double *x, model *cov, double *v){
  TBM2_integr info;
#define MAXSUBDIVISIONS 100
  double a = 0, 
    b = 1, 
    eps = 1e-10; // 1e-15
  int maxsubdivisions = MAXSUBDIVISIONS, 
      lenw = 4 * MAXSUBDIVISIONS;
  double abserr, work[4 * MAXSUBDIVISIONS];
  int subdivisions, integralevaluations, err, iwork[MAXSUBDIVISIONS];

  info.cov = cov; // a bit strange: this is tbm2, but will call tbm3...
  //                 in TBM2NumIntegrFct ...
  info.x = x;
  Rdqags(TBM2NumIntegrFct, (void *) &info, &a, &b, &eps, &eps,
	 v, &abserr, &integralevaluations, &err,
	 &maxsubdivisions, &lenw, &subdivisions, iwork, work);
}


#define DO_NUMERIC 0
void tbm(double *x, model *cov, double *v){
  model *next = cov->sub[TBMOP_COV];
  int 
    fulldim = P0INT(TBMOP_FULLDIM),
    tbmdim = P0INT(TBMOP_TBMDIM);
    //vdimsq = cov->vdim * cov->vdim;

  // if (!hasGaussMethodFrame(cov) && !hasAnyEvaluationFrame(cov)) { BUG; } else
  if (fulldim == tbmdim + 2) tbm3(x, cov, v, (double) tbmdim);
  else if (fulldim == 2 && tbmdim == 1) {
    assert(cov->q);
    if (cov->q[DO_NUMERIC]) tbm2num(x, cov, v);
    else {	 
      // APMI(next);
      TBM2CALL(x, next, v);
      //  ASSERT_GATTER(next);
      // DefList[MODELNR(next)].tbm2(x, next, v);
    }
  }

  else XERR(ERRORTBMCOMBI);

}

void Dtbm(double *x, model *cov, double *v){
  model *next = cov->sub[TBMOP_COV];
  int i,
    fulldim = P0INT(TBMOP_FULLDIM),
    vdim = VDIM0,
    vdimsq = vdim * vdim;
  assert(VDIM0 == VDIM1);
  TALLOC_XX1(v1, vdimsq);
  double 
   tbmdim = (double) P0INT(TBMOP_TBMDIM),
    f = 1.0 + 1.0 / tbmdim;

  BUG;// to do : Taylor-Entwicklung im Ursprung

  if (fulldim == tbmdim + 2) {
    Abl1(x, next, v); // x has dim 2, if turning planes
    Abl2(x, next, v1);
    for (i=0; i<vdimsq; i++) {
      v[i] = v[i] * f + *x * v1[i] / tbmdim;
    }
  }

  else XERR(ERRORTBMCOMBI);

  END_TALLOC_XX1;
  
}

bool numeric_tbm(model *cov) {
  int n = cov->nsub;
  for (int i=0; i<n; i++) if (numeric_tbm(cov->sub[i])) return true;
  return !DefList[COVNR].tbm2;
}

int checktbmop(model *cov) {
  // printf("entering checktbmop\n");
  model *next=cov->sub[TBMOP_COV];
  tbm_param *gp  = &(GLOBAL.tbm);
  int err; 
  ASSERT_ONESYSTEM;

  kdefault(cov, TBMOP_FULLDIM, PisNULL(TBMOP_TBMDIM) || gp->tbmdim >= 0
	   ? gp->fulldim : P0INT(TBMOP_TBMDIM) - gp->tbmdim);
  kdefault(cov, TBMOP_TBMDIM, gp->tbmdim > 0 
	   ? gp->tbmdim : P0INT(TBMOP_FULLDIM) + gp->tbmdim);
  kdefault(cov, TBMOP_LAYERS, (int) gp->layers);
  if ((err = checkkappas(cov)) != NOERROR)  RETURN_ERR(err);
  //PMI0(cov);
  if (!isVariogram(OWNTYPE(0))) {
    // APMI(cov);
    SERR("must be a variogram");
  }

  int 
    tbmdim = P0INT(TBMOP_TBMDIM),
    fulldim = P0INT(TBMOP_FULLDIM),
    vdim = VDIM0,
    storedlayer = P0INT(TBMOP_LAYERS);
  bool layers = storedlayer != NA_LOGICAL ? storedlayer :
    OWNXDIM(0) == tbmdim + 1 && equalsSpaceIsotropic(OWN);
  if(VDIM0 != VDIM1) BUG;
  if (tbmdim >= fulldim)
     SERR4("'%.50s' (=%d) must be less than '%.50s' (=%d)", 
	   KNAME(TBMOP_TBMDIM), tbmdim, KNAME(TBMOP_FULLDIM), fulldim);
  if (OWNLOGDIM(0) > fulldim + layers) RETURN_ERR(ERRORWRONGDIM);

  if (OWNXDIM(0) > tbmdim + layers) {
    SERR("dimension of coordinates does not match reduced dimension of tbm");
  }
 
  if ((err = CHECK_PASSFRAME(next, EvaluationType)) != NOERROR) {
    // APMI(cov);
    RETURN_ERR(err);
  }
  if (next->pref[TBM] == PREF_NONE) RETURN_ERR(ERRORPREFNONE);

  assert(isSpaceIsotropic(OWN));
  assert(isnowVariogram(cov));
  assert(isXonly(OWN));


  set_maxdim(OWN, 0, 0);
  setbackward(cov, next);
  cov->monotone=NOT_MONOTONE;
  set_maxdim(OWN, 0, fulldim + layers);  
  cov->rese_derivs = next->rese_derivs - 1;
  cov->finiterange = 
    (ext_bool) (((fulldim - tbmdim) % 2 == 0) && next->finiterange == true);

  if (vdim > MAXTBMVDIM) 
    SERR2("vdim (%d) exceeds max. value of vdim in tbm3 (%d)", vdim,MAXTBMVDIM);

  // only after being sure that the subsequent model does not cause
  // problems. So the time dimension should be fixed.
  // This is not absolutely safe programming, but should be OK.
  // But difficult to get around for MLE calls that do not allow 
  // for NAs values in integer variables.
  PINT(TBMOP_LAYERS)[0] = layers;
 
  if (fulldim == 2 && tbmdim == 1) QALLOC1(numeric_tbm(cov));

  EXTRA_STORAGE;
  RETURN_NOERROR;
}


Types Typetbm(Types required, model *cov, isotropy_type i) {
  bool layers = P0INT(TBMOP_LAYERS);
  assert(cov->sub[0] != NULL);
  //  printf("i=%d %.50s\n",i, ISO_NAMES[i]);
  if (!isCartesian(i) ||
      ((OWNXDIM(0) == 1) xor equalsIsotropic(i)) ||
       ((OWNXDIM(0) == 2) xor equalsSpaceIsotropic(i)) ||
      (OWNXDIM(0) > 2) ||
      (layers != NA_INTEGER && layers && !equalsSpaceIsotropic(i)) ||
      !equalsXonly(OWNDOM(0))) return BadType;
  return TypeConsistency(required, cov->sub[0], i);
}


bool settbm(model *cov) {
  isotropy_type iso = CONDPREVISO(0);
  if (!isFixed(iso)) return false;

  tbm_param *gp  = &(GLOBAL.tbm);
  // critical parameter TBMOP_LAYER is already treated by SUBMODEL_I
  kdefault(cov, TBMOP_LAYERS, (int) gp->layers);
  set_type(OWN, 0, PREVTYPE(0));
  set_iso(OWN, 0, P0INT(TBMOP_LAYERS) ? DOUBLEISOTROPIC : ISOTROPIC);
  return true; // da nur PREVISO kritisch ist
}

 
bool allowedItbm(model *cov) {
  tbm_param *gp  = &(GLOBAL.tbm);
  bool *I = cov->allowedI;
  kdefault(cov, TBMOP_LAYERS, (int) gp->layers);
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  I[P0INT(TBMOP_LAYERS) ? DOUBLEISOTROPIC : ISOTROPIC] = true;
  return false;
}


void rangetbm_common(model VARIABLE_IS_NOT_USED *cov, range_type *range, bool tbmop ){ 
  int 
    TBMDIM = tbmop ? TBMOP_TBMDIM : TBM_TBMDIM,
    FULLDIM = tbmop ? TBMOP_FULLDIM : TBM_FULLDIM,
    LAYERS = tbmop ? TBMOP_LAYERS : TBM_LAYERS;
 
  range->min[FULLDIM] = 1.0;
  range->max[FULLDIM] = RF_INF;
  range->pmin[FULLDIM] = 1.0;
  range->pmax[FULLDIM] = 100;
  range->openmin[FULLDIM] = false;
  range->openmax[FULLDIM] = true;

  range->min[TBMDIM] = RF_NEGINF;
  range->max[TBMDIM] = RF_INF;
  range->pmin[TBMDIM] = RF_NEGINF;
  range->pmax[TBMDIM] = 100;
  range->openmin[TBMDIM] = false;
  range->openmax[TBMDIM] = true;

  booleanRange(LAYERS);
}

void rangetbmop(model *cov, range_type *range){ 
  rangetbm_common(cov, range, true);
}





int set_stein_q(model *cov, double r, double d) {
  localCE_storage *s = cov->calling->SlocalCE;
  double C0, phi0, phi1, phi2, 
    zero = 0.0, 
    rP1 = r + 1.0,
    rM1 = r - 1.0, 
    dsq = d * d;
  localvariab *q = s->q2;

  COV(&zero, cov, &C0);
  COV(&d, cov, &phi0);
  Abl1(&d, cov, &phi1); // derivative of
  phi1 *= d;             // scaled fctn at 1
  Abl2(&d, cov, &phi2); // 2nd derivative of
  phi2 *= dsq;            // scaled fctn at 1
  q->R =  r * d;   
  double dummy = (phi2 - phi1) / (3.0 * r * rP1);
  q->intrinsic.B  = (r == 1.0) ? 0.0 : dummy / (rM1 * dsq);
  q->intrinsic.A2 = (dummy - phi1 / 3.0 - phi2 / 6.0) / dsq;
  q->intrinsic.A0 = 0.5 * rM1 / rP1 * phi2 + phi1 / rP1 - phi0;

  if ((q->intrinsic.B < 0.0) || (q->intrinsic.A2 < 0.0) ||
      (q->intrinsic.A0 + C0 < 0.0)) {
    
    return MSGLOCAL_INITINTRINSIC;
  }

  RETURN_NOERROR;
}

void getlocalconstants(double phi0, double phi1, double phi2, double phi3,
		       double radius, double m, double n, double l,
		       double *aa, double *bb, double *cc, double *constant) {
  double phi3r2 = phi3*radius*radius,
    phi2r = phi2 * radius;

	*aa = -(phi3r2 + phi2r * (m + l - 3) + phi1 *
	   (m - 1) * (l - 1)) / (n * (m - n) * (l - n) * POW(radius, n - 1));
	*bb = -(phi3r2 + phi2r * (n + l - 3) + phi1 *
	   (n - 1) * (l - 1)) / (m * (n - m)*(l - m) * POW(radius, m - 1));
	*cc = -(phi3r2 + phi2r*(m + n - 3) + phi1*
		(m - 1)*(n - 1)) / (l*(m - l)*(n - l) * POW(radius, l - 1));
	*constant = -phi0 + *aa * POW(radius, n) +
	  *bb * POW(radius, m) + *cc * POW(radius, l);
}

int set_cutoff_q2variate(model *cov, double a, double d) {
  // auf modell ebene, d.h. co->sub[0] oder stein->sub[0]
  double 
    phi0[4], phi1[4], phi2[4], 
    a2 = a * a;     
  // radius = MAXINT;
  //double roots[3][2];
  
  
  assert(cov->calling->SlocalCE != NULL);
  localCE_storage *s = cov->calling->SlocalCE;
  localvariab *q = s->q2;
  
  COV(&d, cov, phi0);
  // COV(&d,  cov->sub[0], &phi01);
  Abl1(&d, cov, phi1);
  Abl2(&d, cov, phi2);
  
  //bivariate Whittle  Cauchy Stable

  //assert(false);
   
  if (VDIM0 > LOCALCE_MAXVDIM || VDIM1 > LOCALCE_MAXVDIM) BUG;
  if (phi1[1] != phi1[2] )
    return MSGLOCAL_NOTSYMMETRICMULTIVARIATE;
  
  for (int j = 0; j < 4; j++) {
    q[j].a = a;
    if ( phi0[j] > 0 ) {
      q[j].cutoff.b = phi2[j]*phi2[j]*phi2[j]/(108.0*phi1[j]*phi1[j]);
      q[j].cutoff.theor = POW(1.0 - 3.0 * a2 * phi1[j] / (d* phi2[j]), 1.0 / a);
      q[j].R = d * q[j].cutoff.theor;
      q[j].cutoff.asqrtr = POW(q[j].R, a);
      q[j].cutoff.constant = phi0[j] -q[j].cutoff.b*POW(-3.0*phi1[j]/phi2[j], 4);
    } else {
      q[j].cutoff.b = 0;
      q[j].cutoff.theor = 0;
      q[j].R = 0;
      q[j].cutoff.asqrtr = 0;
      q[j].cutoff.constant = 0;
    }
    
  }

  
  if (q[1].R > q[0].R || q[1].R > q[3].R ) {
    // printf("%10g > %10g | %10g > %10g\n", q[1].R, q[0].R, q[1].R, q[3].R);
    //    APMI(cov->calling);
    return MSGLOCAL_WRONGRADII;
  }
  
  
  //PMI(cov);

  
  RETURN_NOERROR;
}





int set_cutoff_q(model *cov, double a, double d) {
  assert(VDIM0 > 0 && VDIM0 <= LOCALCE_MAXVDIM);
  if (VDIM0 > 1) return set_cutoff_q2variate(cov, a, d);
  localCE_storage *s = cov->calling->SlocalCE;

    // auf modell ebene, d.h. co->sub[0] oder stein->sub[0]
  double
    phi0, phi1,
    phi2 = RF_NA,
    phi3 = RF_NA,
    phi4 = RF_NA,
    a2 = a * a,
    //   d2 = d*d,
    radius = UNSET;
  localvariab *q = s->q2;
  q->a = a;

  //  if ((err = covcpy(&neu, cov)) != NOERROR) RETURN_ERR(err);
  // DefList[cov->gatternr].cov(&d, neu, &phi0);
  // DefList[cov->gatternr].D(&d, neu, &phi1);
  //FREE(neu);
  
  COV(&d, cov, &phi0);
  Abl1(&d, cov, &phi1);


   
  //If it is a variogram, then we must determine a constant such that Const - Variogram is a valied covariance
  // The constant depends on a and d
  //    cov->calling->typus
  
    //  if (cov->typus == V ariogramType) {
  if (equalsnowVariogram(cov)) {
    
    if (a == 0.5) {
      //COV(&d, cov, &q[CUTOFF_CONSTANT]);
      
      //q[CUTOFF_B] = -2*phi1*SQRT(d);
      //q[CUTOFF_THEOR] =  POW(1.0 - 0.5* (q[CUTOFF_CONSTANT] + phi0) / phi1 /d, 1/a);
      //q[LOCAL_R] = d * q[CUTOFF_THEOR];

      q->cutoff.constant = -phi0;
      q->R = d;      
      q->cutoff.asqrtr = POW(q->R, a);
      
    } else if (a == 1 ) {
      //second derivative
      Abl2(&d, cov, &phi2);
      if (phi2 <= 0 || phi1 >= 0 || phi0 >= 0 ) {
                return MSGLOCAL_SIGNPHISND;
      }
      //if second derivative is positive at d, then
      //q[CUTOFF_CONSTANT] = -phi0+phi1*phi1/(2*phi2) + EPSILON_C;
      //q[CUTOFF_B] = 0.25*phi1*phi1/(q[CUTOFF_CONSTANT] + phi0);
      //q[CUTOFF_THEOR] = POW(1.0 - 2 * ( q[CUTOFF_CONSTANT] + phi0) / phi1 /d, 1/a);
      
      q->cutoff.constant = -phi0+phi1*phi1/(2*phi2);
      q->cutoff.b =  phi2/2;
      q->R = d  - phi1/phi2;
      q->cutoff.asqrtr = POW(q->R, a);
      
    } else if (a == CUTOFF_THIRD_CONDITION) {
      
      Abl2(&d, cov, &phi2);
      Abl3(&d, cov, &phi3);
      Abl4(&d, cov, &phi4);
      
      double roots[3][2];
      
      //de Wijsian model
      int 
	n = 4,
	m = 6,
	l = 7;
      
      q->cube.N = n;
      q->cube.M = m;
      q->cube.L = l;
      
      cubicsolver( phi4, (n+m+l-6)*phi3,
		   (n*m + n*l + m*l - 3*(n+m+l) + 7)*phi2,
		   (n-1)*(m-1)*(l-1)*phi1 ,
		   roots);
      
      
      //find a positive root of the cubic polynomial
      for (int i = 0; i < 3; i++) {
	if (roots[i][1] == 0 && roots[i][0] > radius ) {
	  radius = roots[i][0];
	}
      }
      //if no positive root, cannot do anything
      if (radius <= 0.0) return MSGLOCAL_NOPOSITIVEROOT;
      
      

	  getlocalconstants(phi0, phi1, phi2, phi3, radius, m, n, l,
			&(q->cube.A), &(q->cube.B), &(q->cube.C),
			&(q->cube.constant));
      
      if (q->cube.constant <= 0.0) return MSGLOCAL_SIGNPHI;
      q->R = radius + d;
      // q->cutoff.asqrtr = radius;  Olga, line looks strange
      
    } else BUG;
    
  } else {
    
    //now we are in a positive definite function case
    if (phi0 <= 0.0) return MSGLOCAL_SIGNPHI;
    if (phi1 >= 0.0) return MSGLOCAL_SIGNPHIFST;
    
    
    
    //find  parameters of the method
    if (a != CUTOFF_THIRD_CONDITION) {
      phi1 *= d;
      
      if (phi1 >= 0.0) return MSGLOCAL_SIGNPHIFST;
      
      q->cutoff.b =//sound qconstant even if variance of submodel is not 1
	POW(-phi1 / (2.0 * a2 * phi0), 2.0 * a) * phi0 / POW(d, 2.0 * a2);
      //  assert(false);
      q->cutoff.theor = POW(1.0 - 2.0 * a2 * phi0 / phi1, 1.0 / a);
      q->R = d * q->cutoff.theor;
      q->cutoff.asqrtr = POW(q->R, a);
      
    } else if (a == CUTOFF_THIRD_CONDITION){
      
      //a = CUTOFF_THIRD_CONDITION

      Abl2(&d, cov, &phi2);
      Abl3(&d, cov, &phi3);
      Abl4(&d, cov, &phi4);
      
      
      double roots[3][2];
      
      //Whittle  model
      int
	n = 5,
	m = 6,
	l = 7;
      q->cube.N = n;
      q->cube.M = m;
      q->cube.L = l;
      
      cubicsolver( phi4, (n+m+l-6)*phi3,
		   (n*m + n*l + m*l - 3*(n+m-l) + 7)*phi2,
		   (n-1)*(m-1)*(l-1)*phi1 ,
		   roots);
      
      //find a positive root of the cubic polynomial
      for (int i = 0; i < 3; i++) {
	if (roots[i][1] == 0 && roots[i][0] > radius ) {
	  radius = roots[i][0];
	}
      }
      
      //if no positive root, cannot do anything
      if (radius <= 0.0) return MSGLOCAL_NOPOSITIVEROOT;
      
      getlocalconstants(phi0, phi1, phi2, phi3, radius, m, n, l,
			&(q[0].cube.A), &(q[0].cube.B), &(q[0].cube.C),
			&(q[0].cube.constant));
      
      if (q->cube.constant <= 0.0) return MSGLOCAL_SIGNPHI;
      
      q->R = radius + d;
      
      
    } else BUG;
  }
  RETURN_NOERROR;
}


int check_local(model *cov,  
		 Methods method, 
		 getlocalparam init, // next->coinit
		 set_local_q_type set_local) {
  // PMI(cov);
  localCE_storage *s = cov->SlocalCE;
  assert(s != NULL);
  location_type *loc = Loc(cov);
  int i, msg, 
    dim = OWNLOGDIM(0),
    err=NOERROR;
  localvariab
    *q = s->q,
    *q2 = s->q2;
  double d=RF_NA;
  model *next = cov->sub[0];
  localinfotype li;
  //ce_param *gp  = &(GLOBAL.localce); // ok

  //printf("\n \n entering check local from %d:%.50s\n", cov->CALLINGNR,
  //DefList[cov->CALLINGNR].name);

  if ((err = CHECK(next, dim,  1,
		   method == CircEmbedCutoff ? PosDefType : VariogramType,
		   OWNDOM(0), OWNISO(0), SUBMODEL_DEP, EvaluationType)) != NOERROR) {
    //  printf("\n \n \n ::: checl_local :: if :::  dim = %d --- next->vdim[0] = %d::: \n \n \n", dim, next->vdim[0]);
    if ( method != CircEmbedCutoff ||
	 (err = CHECK(next, dim,  1,
		      VariogramType,
		      OWNDOM(0), OWNISO(0), SUBMODEL_DEP, EvaluationType)) != NOERROR)
      RETURN_ERR(err);
  }
  
  if (VDIM0 > LOCALCE_MAXVDIM) SERR("vdim of submodel must be less than 3");


  // no setbackward ?!
  setbackward(cov, next); 
  if (next->pref[method] == PREF_NONE) RETURN_ERR(ERRORPREFNONE);

  if (init == NULL) RETURN_ERR(ERRORUNKNOWNMETHOD);

  if (PisNULL(pLOC_DIAM)) {  
    double 
      diameter = GetDiameter(loc);
    if (PL>=PL_DETAILS) { LPRINT("diameter %10g\n", diameter); }
    kdefault(cov, pLOC_DIAM, diameter);
  } else {
    d = P0(pLOC_DIAM);
  }

  if (PisNULL(pLOC_A)) {
    if (DefList[NEXTNR].implemented[method]) {
      assert(init != NULL);

      init(next, &li);
      if (li.instances == 0) {
	SERR("parameter values do not allow for finding second parameter");
      }
      q->R = RF_INF;

      msg = XMSGLOCAL_FAILED;
      for (i=0; i<li.instances; i++) {
	//	printf("li %d %d %d %10g %d\n", i, li.instances, li.msg[i], li.value[i],  msg);
	if (li.msg[i] <= msg) {
	  //	
	  err = set_local(next, li.value[i], d);
	  //
	  if (err == NOERROR && (li.msg[i] < msg || q2->R < q->R)){
	    MEMCOPY(q, q2, sizeof(localvariabArray));
	    msg = li.msg[i];
	    if (PisNULL(pLOC_A)) {
	      PALLOC(pLOC_A, 1, 1); 
 	    }
	    P(pLOC_A)[0] = li.value[i]; // co: a; ie: rawr
	  }
	}
      }
      q->msg = msg;

      if (msg == MSGLOCAL_FAILED) {
	// APMI(cov);
	SERR2("'%.50s' did not work out correctly. Likely, invalid parameter \
values have been chosen for '%.50s'", NICK(cov), NICK(next));
      }
    } else {
      SERR("2nd parameter is neither given nor can be found automatically");
    }

  } else {
    if (cov->ncol[pLOC_A] != 1 || cov->nrow[pLOC_A] != 1) 
      SERR1("'%.50s' must be a scale", KNAME(pLOC_A));
    err = set_local(next, P0(pLOC_A), d);
    MEMCOPY(q, q2, sizeof(localvariabArray));
  }

  cov->pref[CircEmbed] = 5;

  RETURN_ERR(err);
}


void co(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  localCE_storage *s = cov->SlocalCE;
  localvariab *q = s->q;
  double  y=*x, 
    diameter = P0(pLOC_DIAM),  
    a = P0(pLOC_A);
   
  // assert(hasGaussMethodFrame(cov) || hasAnyEvaluationFrame(cov));

  //VDIM0 or next->vdim[0]?
  //if ( cov->SlocalCE->is_multivariate_cutoff == true ) {
  
  if ( VDIM0 > 1 ) {
    if (y <= diameter) {
      COV(x, next, v);
      for (int i = 0; i < 4; i++) {
        v[i] = v[i] - q[i].cutoff.constant;
      }
      } else {
      for (int i = 0; i < 4; i++) {
	if (y >= q[i].R) v[i] = 0.0;
	else v[i] = q[i].cutoff.b*POW(q[i].cutoff.asqrtr - POW(y, a), 4.0 * a);
      }
    }
    
  } 
  else {
    if (y <= diameter) {
      //If it is a variogram (actually -variogram), add a constant. If it is positive definite, go below
      if (isnowVariogram(next)) {
	COV(x, next, v);
	*v = *v + q->cutoff.constant;
      } else
	{
	  COV(x, next, v)
            }
      
    }
    
    else {
      if (y >= q->R) *v =  0.0;
      else {
	if (a != CUTOFF_THIRD_CONDITION)
	  *v = q->cutoff.b * POW(q->cutoff.asqrtr - POW(y, a), 2.0 * a);
	else 
	  *v = q->cube.A * POW((q->R - y), q->cube.N) +
	    q->cube.B * POW((q->R - y), q->cube.M) +
	    q->cube.C * POW((q->R - y), q->cube.L);
	//BUG;
      }
    }    
  }
}
int check_co(model *cov) {
  model *next = cov->sub[0];
  int err;
  
  NEW_STORAGE(localCE);
  assert(cov->SlocalCE != NULL);
  // 1 ? 2? More than 2? or inside set_cutoff_q2variate
  //PMI(cov);assert(next->vdim[0]>0);

   err = check_local(cov, CircEmbedCutoff,
		    DefList[NEXTNR].coinit, set_cutoff_q);

  RETURN_ERR(err);
}


bool alternativeparam_co(model VARIABLE_IS_NOT_USED *cov){
  return false;
}

void range_co(model VARIABLE_IS_NOT_USED *cov, range_type *range){

  range->min[pLOC_DIAM] = 0.0; //  CUTOFF_DIAM
  range->max[pLOC_DIAM] = RF_INF;
  range->pmin[pLOC_DIAM] = 1e-10;
  range->pmax[pLOC_DIAM] = 1e10;
  range->openmin[pLOC_DIAM] = true;
  range->openmax[pLOC_DIAM] = true;

  range->min[pLOC_A] = 0.0; // cutoff_a
  range->max[pLOC_A] = RF_INF;
  range->pmin[pLOC_A] = 0.5;
  range->pmax[pLOC_A] = 2.0;
  range->openmin[pLOC_A] = true;
  range->openmax[pLOC_A] = true;
}


void Stein(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  localCE_storage *s = cov->SlocalCE;
  localvariab *q = s->q;
  double y=*x, z, 
    diameter = P0(pLOC_DIAM);

  // assert(hasAnyEvaluationFrame(cov) || hasGaussMethodFrame(cov));

  if (y <= diameter) {
    COV(x, next, v);
    *v += q->intrinsic.A0 + q->intrinsic.A2 * y * y;
  } else {
     z = q->R - y;
     if (z <= 0.0) *v = 0.0;
     else *v = q->intrinsic.B * z * z * z / y;
  }
}

int check_Stein(model *cov)
{  
  model *next = cov->sub[0]; // dito
  NEW_STORAGE(localCE);
  assert(cov->SlocalCE != NULL);
  return check_local(cov, CircEmbedIntrinsic, DefList[NEXTNR].ieinit,
		    set_stein_q);
}  

bool alternativeparam_Stein(model *cov) {
  P(pLOC_A)[0] *= 2.0;
  return true;
}

void range_Stein(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[pLOC_DIAM] = 0.0; 
  range->max[pLOC_DIAM] = RF_INF;
  range->pmin[pLOC_DIAM] = 0.01;
  range->pmax[pLOC_DIAM] = 100;
  range->openmin[pLOC_DIAM] = true;
  range->openmax[pLOC_DIAM] = true;

  range->min[pLOC_R] = 1.0; // stein_r
  range->max[pLOC_R] = RF_INF;
  range->pmin[pLOC_R] = 1.0;
  range->pmax[pLOC_R] = 20.0;
  range->openmin[pLOC_R] = false;
  range->openmax[pLOC_R] = true;
}

