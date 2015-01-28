/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Gneiting's space-time covariance models and related models

 Copyright (C) 2006 -- 2015 Martin Schlather

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


#include <math.h>
 
#include "RF.h"
#include "Operator.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


#define AVESTP_MINEIGEN 2
#define AVESTP_LOGDET 3
#define AVESTP_V 4
#define AVESTP_LOGV 5
#define AVESTP_LOGMIXDENS 6 
#define TOTALAVESTP (AVESTP_LOGMIXDENS + 1)

#define AVE_A 0
#define AVE_Z 1
#define AVE_SPACETIME 2
#define AVE_PHI 0
#define AVE_GAUSS 1


void kappa_ave(int i, cov_model *cov, int *nr, int *nc){
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  int dim = spacetime ? cov->tsdim-1 : cov->tsdim;
  *nr = (i==AVE_A || i==AVE_Z) ? dim : 1;
  *nc = (i==AVE_A) ? dim : i < CovList[cov->nr].kappas ? 1 : -1;
}


void ave(double *h, cov_model *cov, double *v) {
  // f = uAu +zu; 
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  cov_model *next = cov->sub[0];
  int i,j,k,d,
    dim = spacetime ? cov->tsdim - 1 : cov->tsdim;
  double  detEplus2B, Ah[AveMaxDim], Eplus2B[AveMaxDim], 
    dummy, 
    hh,
    *A = P(AVE_A),
    *z = P(AVE_Z),
    c = spacetime ? h[cov->tsdim-1] : 0.0; // Sockelwert fuer c

  hh = 0.0;
  for (k=d=0; d<dim; d++) {
    for (dummy = 0.0, j=0; j<dim; j++) dummy += h[j] * A[k++];
    Ah[d] = dummy;
    c += z[d] * h[d];
    hh += h[d] * h[d]; 
  }

  for (j=d=0; d<dim; d++) {
    for (i=0; i<dim; i++) {
      Eplus2B[j++] = 2.0 * Ah[d] * Ah[i];
    }
    Eplus2B[j - dim + d] += 1.0;
  }

  det_UpperInv(Eplus2B, &detEplus2B, dim); // Eplus2B is not its inverse !

  double y = sqrt(0.5 * hh  + c * c * (1.0 - 2.0 * xUx(Ah, Eplus2B, dim)));
  COV(&y, next, v);
  *v /= sqrt(detEplus2B);
}


int checkave(cov_model *cov) {
  cov_model *next = cov->sub[0];
  bool
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  int i, j, err,
    dim =  cov->tsdim,
    spdim = spacetime ? dim - 1 : dim;
  double 
    *A = P(AVE_A);
  char msg[2][4] = {"d", "d-1"};

   if (cov->xdimown < 2) SERR("The spatial dimension must be at least 2.");
 
 
  if (dim > AveMaxDim)
    SERR2("For technical reasons max. dimension for ave is %d. Got %d.", 
	  AveMaxDim, dim);

  if (cov->ncol[AVE_A] != spdim || cov->nrow[AVE_A] != spdim) 
    SERR5("A not %sx%s matrix, but %dx%d (dim=%d)", msg[spacetime], 
	  msg[spacetime], cov->ncol[AVE_A], cov->nrow[AVE_A], spdim);

  if (cov->ncol[AVE_Z] != 1 || cov->nrow[AVE_Z] != spdim) 
    SERR1("z not (%s)-dim vector", msg[spacetime]);

  for (i=0; i<spdim; i++)
    for (j=i+1; j<spdim; j++)
      if (A[i + j * spdim] != A[j + i * spdim]) {
	A[j + i * spdim] = A[i + j * spdim];
	warning("A is not symmetric -- lower part used");
      }
  // naechste Zeile stimmt nicht mit Bernoulli ueberein

  kdefault(cov, AVE_SPACETIME, TRUE);
  if ((err = checkkappas(cov)) != NOERROR) return err;

  if (cov->xdimprev != cov->tsdim || cov->xdimprev != cov->tsdim)
    return ERRORDIM;
  if ((err = CHECK(next, dim, 1, PosDefType, XONLY, ISOTROPIC,
		   SCALAR, cov->role // ROLE_COV changed 20.7.14 wg spectral
		   )) != NOERROR) return err;
  next->delflag = DEL_COV; // set gatternr=nr, since non-negativity ensured
  if (!isNormalMixture(next->monotone)) return ERRORNORMALMIXTURE;
  if (CovList[next->nr].spectral == NULL) return ERRORSPECTRAL; // nicht gatter

  //  updatepref(cov, next); ## gute idee?
  if (next->pref[SpectralTBM] == PREF_NONE) 
    cov->pref[RandomCoin] = cov->pref[Average] = PREF_NONE;

  // kein setbackward??
  
  // no setbackard
  return NOERROR;
}



void rangeave(cov_model VARIABLE_IS_NOT_USED *cov, range_type* ra){
  int i;
  
  for (i=0; i<=1; i++) {
    ra->min[i] = RF_NEGINF;
    ra->max[i] = RF_INF;
    ra->pmin[i] = -10.0;
    ra->pmax[i] = 10.0;
    ra->openmin[i] = true;
    ra->openmax[i] = true;
  }  
  ra->min[2] = 0;
  ra->max[2] = 1;
  ra->pmin[2] = 0;
  ra->pmax[2] = 1;
  ra->openmin[2] = false;
  ra->openmax[2] = false;
}



void sd_avestp(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *S, int dim, double *sd){
  /////  window_info *w = &(S->window);
  int d;
  double b, alphamin, x2, InvSqrt2a, EmA,
    *q = cov->q;
  // see article/GEOSTATS/simuspacetime/simuspacetime2008/simuspacetime.tex 
  // for the reasoning of these calculations

  BUG; 

  assert(cov->role == Average);
  q[AVESTP_LOGV] = log(q[AVESTP_V]);
  for (x2=0.0, d=0; d<dim; d++) {
    double lensimu = RF_NA; ////  w->max[d] - w->min[d];
    x2 += lensimu * lensimu;
  }
  // x2 *= 0.25;
  b = 3.0 * q[AVESTP_V] * x2 / dim;
  alphamin = (4.0 + 4.0 * b - 2.0 * sqrt(4 * b * b + 8.0 * b + 1.0)) / 3.0;
  InvSqrt2a = 1.0 / sqrt(2.0 * alphamin * 6.0 * q[AVESTP_V]);
  *sd = InvSqrt2a;
  EmA = 1.0 - alphamin;
  cov->mpp.maxheights[0] = exp(-0.5 * log(EmA) - 0.25 * log(alphamin) + b / EmA -
			   2 * x2); // proportional zum dritten Moment !

  /*
   double radius = 
    sqrt((-9 // so e^{-9} as threshold
	  - 0.25 * dim * (q[AVESTP_LOGV] - 1.14473) // log pi
	  - 0.25 * q[AVESTP_LOGDET]
	  //+ 0.5 * cov_a->logdens
	  -  q[AVESTP_LOGMIXDENS]
	  ) / ( - q[AVESTP_MINEIGEN] * q[AVESTP_V]) ); // ???
  assert(radius > 0);
  */

}


int structAve(cov_model *cov, cov_model **newmodel) { 
  cov_model *shape, *gaussmix;
  int err;
 
  ASSERT_NEWMODEL_NOT_NULL;
  
  if (cov->role != Average) ILLEGAL_ROLE;

  if ((err = covcpy(newmodel, cov)) != NOERROR) return err;
  shape = *newmodel;
  shape->nr = SHAPEAVE;
  addModel(shape, AVE_GAUSS, GAUSS);
  gaussmix = shape->sub[AVE_GAUSS];
  gaussmix->tsdim = 1;
  gaussmix->role = ROLE_GAUSS;
  gaussmix->method = SpectralTBM;
  return NOERROR;
}




void  logshapeave(double *x, cov_model *cov, double *v, double *Sign) {
    // nur stationaer
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  int d, j, k,
    dim = spacetime ? cov->tsdim - 1 : cov->tsdim;
  double f, dummy, r2,
    *A = P(AVE_A),
    *z = P(AVE_Z),
    t = spacetime ? x[cov->tsdim-1] : 0.0,
    *q = cov->q;
 
  f = r2 = 0.0;
  for (k=d=0; d<dim; d++) {
    r2 += x[d] * x[d];
    for (dummy = z[d], j=0; j<dim; j++) dummy += x[j] * A[k++];
    f += dummy * x[d];
  }

  static bool avewarning=true;
  if (avewarning) warning("is exponent of V correct?"); avewarning=false;
   
  v[0] = 0.25 * dim * q[AVESTP_LOGV] // V^{d/4 oder d/2} !!!!!!!!!!!
    - 0.5 * (LOG2 - dim * M_LN_SQRT_PId2) - r2
    // LOG2 : wegen spectral  simulation; 
    // zweiter Term fuer logg 
    //+ CovList[phi->nr].logmixdens(x, q[AVESTP_LOGV], phi); /* g */// nicht gatternr
    ;
  Sign[0] = 1.0;
  double phase = q[AVERAGE_YPHASE] + q[AVERAGE_YFREQ] * (f - t); // Y
  Sign[1] =  phase > 0.0 ? 1.0 : phase < 0.0 ? -1.0 : 0.0;
  v[1] = log(fabs(phase));
}

int check_shapeave(cov_model *cov) {
  if (cov->sub[AVE_GAUSS] == NULL)
    SERR1("both submodels must be set to '%s'", CovList[GAUSS].nick);
  cov->mpp.maxheights[0] = RF_NA;
  return checkave(cov); // !! not next
}

int init_shapeave(cov_model *cov, gen_storage *s) { 
  ASSERT_GAUSS_METHOD(Average);
  cov_model
    *gaussmix = cov->sub[AVE_GAUSS];
  double sd,
    *q = cov->q;
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  int err = NOERROR,
    dim = spacetime ? cov->tsdim - 1 : cov->tsdim;
  

  q[AVESTP_V] = 0.0;
  q[AVESTP_MINEIGEN] = 1.0; 
  q[AVESTP_LOGDET] = 0.0;
  sd_avestp(cov, s, dim, &sd); // sd->gauss

  assert(cov->vdim[0] == 1);  assert(cov->vdim[0] == cov->vdim[1]);

  if (cov->mpp.moments >= 0) {
    cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0; 
    if (cov->mpp.moments >= 1) {
      if ((err = INIT(gaussmix, cov->mpp.moments, s)) != NOERROR) return err;
      if (cov->mpp.moments >= 2) {
	cov->mpp.mM[2] = 1.0;
      }
    }
  }
  return err;
}


void do_shapeave(cov_model *cov, gen_storage *S) { 
  // Simulation of V; sopee Bernoulli Sec. 4.2
   cov_model    
    *aveGAUSS = cov->sub[AVE_GAUSS],
    *phi = cov->sub[AVE_PHI];
  double spec_ret[StpMaxDim], sd,
    *q = cov->q;
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  int 
    dim = spacetime ? cov->tsdim - 1 : cov->tsdim;
  
  CovList[phi->nr].drawmix(phi, q + AVESTP_V); // nicht gatternr
  sd_avestp(cov, S, dim, &sd); // sd->gauss
  
  BUG;

  SPECTRAL(aveGAUSS, S, spec_ret);  // nicht gatternr
  q[AVERAGE_YFREQ] = *spec_ret * q[AVESTP_V];  
  q[AVERAGE_YPHASE] = TWOPI * UNIFORM_RANDOM;

  
  BUG; // what to do with the next line?
 }





/* coxgauss, cmp with nsst1 !! */
// C = 2 (C + 4 M H M), H = h h^t
// a = t - h M h - zh
// exp(- 0.5 * (h *h + 2 a^2 - mu C mu)) // stimmen die Vorzeichen??
// mu = h - 2 a M h
/* cox, cmp with nsst1 !! */
// coxisham
#define COX_MU 0
#define COX_D 1
#define COX_BETA 2

void GetEu2Dinv(cov_model *cov, double *x, int dim, 
		double *det, double *Eu2Dinv,
		double *newxsq, double *newx, double *z) {
    double 
	y[CoxMaxDim],
      *V = P(COX_MU),
      *D= P(COX_D),
      beta = P0(COX_BETA),
      t = x[dim],
      t2 = pow(fabs(t), beta); // standard t^2
    int d,
	dimP1 = dim + 1,
	dimsq = dim * dim;
    for (d=0; d<dim; d++) {
      y[d] = x[d] - t * V[d];
  }
  
  for (d=0; d<dimsq; d++) {
      Eu2Dinv[d] = t2 * D[d];
  }
  for (d=0; d<dimsq; d+=dimP1)  Eu2Dinv[d] += 1.0; // D + E
  det_UpperInv(Eu2Dinv, det, dim);
  *newxsq = xUxz(y, Eu2Dinv, dim, z);
  *newx = sqrt(*newxsq);
}

void cpyUf(double *Eu2Dinv, double factor, int dim, int tsdim, double *v) {
    // Eu2Dinv has dimension dim^2; v dimension tsdim^2
    // Eu2Dinv is copied into the upper left part of v and 
    // multiplied by factor
  int d, i, k, j,
	tsdimsq = tsdim * tsdim;
  
  for (i=0; i<tsdimsq; v[i++] = 0.0);
  for (i=0; i<dim; i++) {
      for (d=i * tsdim, k=i * dim, j=0; j<=i; j++)
	  v[d++] = Eu2Dinv[k++] * factor; 
      for (k += dim - 1; j<dim; j++, k+=dim) { 
	  v[d++] = Eu2Dinv[k] * factor;
      }
  }
}

void addzzT(double *v, double factor, double *z, int dim, int tsdim) {
    // z has dimension dim; v dimension tsdim^2
    // zzT is copied into the upper left part of v after being 
    // multiplied by factor
   
    int i,j,k;
    for (i=0; i<dim; i++) {
	k = i * tsdim;
	for (j=0; j<dim; j++) {
	    v[k++] += z[i] * z[j] * factor;
	}
    }
}


void kappa_cox(int i, cov_model *cov, int *nr, int *nc){
    switch (i) {
	case COX_MU :
	    *nc = 1;
	    *nr = cov->tsdim - 1;
	    break;
	case COX_D  :
	    *nc = *nr = cov->tsdim - 1;
	    break;
	case COX_BETA :
	    *nc = *nr = 1;
	    break;
	default:  *nc = *nr = -1;
    }
}

void cox(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int dim = cov->tsdim - 1,
      dimsq = dim * dim;
  double det, newx,  newxsq;
  ALLOC_EXTRA(Eu2Dinv, dimsq);
  GetEu2Dinv(cov, x, dim, &det, Eu2Dinv, &newxsq, &newx, NULL);
   
  COV(&newx, next, v);
  *v /= sqrt(det);
} 

void coxhess(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int tsdim = cov->tsdim,
      dim = tsdim - 1,
      dimsq = dim * dim;
  double z[CoxMaxDim], det, newx, newxsq, phiD, phiD2;

  ALLOC_EXTRA(Eu2Dinv, dimsq);

  GetEu2Dinv(cov, x, dim, &det, Eu2Dinv, &newxsq, &newx, z);

  Abl2(&newx, next, &phiD2);  
  if (newxsq == 0.0) {
    cpyUf(Eu2Dinv, phiD2 / sqrt(det), dim, tsdim, v);
  } else {
    Abl1(&newx, next, &phiD);
    cpyUf(Eu2Dinv, phiD / (sqrt(det) * newx), dim, tsdim, v);
    addzzT(v, (phiD2 - phiD/newx) / (sqrt(det) * newxsq), z, dim, tsdim);
  }
} 
 

void coxnabla(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int d,
      tsdim = cov->tsdim,
      dim = tsdim - 1,
      dimsq=dim * dim;
  double z[CoxMaxDim], det, newx, newxsq, phiD, factor;
  
  ALLOC_EXTRA(Eu2Dinv, dimsq);
  GetEu2Dinv(cov, x, dim, &det, Eu2Dinv, &newxsq, &newx, z); 

  if (newxsq == 0.0) {
      for (d=0; d<=dim; d++)  v[d] = 0.0;
  } else {    
    newx = sqrt(newxsq);
    Abl1(&newx, next, &phiD);
    factor = phiD / (det * newx); 
    for (d=0; d<dim; d++) {
	v[d] = factor * z[d];
    }
    for (d=0; d<tsdim; v[d++]=0.0);
  }

 }



int checkcox(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err, i, 
    dim = cov->tsdim - 1,
    dimsq = dim * dim; 

  if (cov->xdimown < 2) SERR("The space-time dimension must be at least 2.");
   
  if (cov->ncol[COX_MU] != 1 || cov->nrow[COX_MU] != dim) {
    if (cov->ncol[COX_MU] == dim && cov->nrow[COX_MU] == 1) {
      cov->nrow[COX_MU] = dim;
      cov->ncol[COX_MU] = 1; 
    } else {
      SERR3("mu is not given or not a vector of dimen. %d (nrow=%d ncol=%d)",
	    dim, cov->nrow[COX_MU], cov->ncol[COX_MU]);
    }
  }

  // is matrix positive definite?
  if (PisNULL(COX_D)) {
    PALLOC(COX_D, dim, dim);
    for (i=0; i<dimsq; i++) P(COX_D)[i] = 1.0;
  } else {
    if (!is_positive_definite(P(COX_D), dim))
      SERR("D is not (strictly) positive definite");
  }

  kdefault(cov, COX_BETA, 2.0);

  if ((err = checkkappas(cov)) != NOERROR) return err;
  if ((err = CHECK(next, dim,  1, PosDefType, XONLY, ISOTROPIC,
		     SCALAR, cov->role // ROLE_COV changed 20.7.14 wg spectral
		   )) != NOERROR)
    return err;

  if (cov->tsdim != 3)  cov->pref[SpectralTBM] = PREF_NONE;;

  next->delflag = DEL_COV; // set gatternr=nr, since non-negativity ensured
  if (!isNormalMixture(next->monotone)) return ERRORNORMALMIXTURE;
  if (CovList[next->nr].spectral == NULL) return ERRORSPECTRAL; // nicht gatternr
  
  // no setbackard
  updatepref(cov, next);
  if (P0(COX_BETA) != 2.0) cov->pref[SpectralTBM] = 0;

  cov->hess = true;

  EXTRA_STORAGE;	 
  return NOERROR;
} 

void rangecox(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){

  range->min[COX_MU] = RF_NEGINF;
  range->max[COX_MU] = RF_INF;
  range->pmin[COX_MU] = -100.0;
  range->pmax[COX_MU] = +100.0;
  range->openmin[COX_MU] = true;
  range->openmax[COX_MU] = true;

  range->min[COX_D] = RF_NEGINF;
  range->max[COX_D] = RF_INF;
  range->pmin[COX_D] = -100.0;
  range->pmax[COX_D] = +100.0;
  range->openmin[COX_D] = false;
  range->openmax[COX_D] = false;  

  range->min[COX_BETA] = 0.0;
  range->max[COX_BETA] = 2.0;
  range->pmin[COX_BETA] = 0.1;
  range->pmax[COX_BETA] = 2.0;
  range->openmin[COX_BETA] = true;
  range->openmax[COX_BETA] = false;  
}

int initcox(cov_model *cov, gen_storage *s) {
  ASSERT_GAUSS_METHOD(SpectralTBM);
  cov_model *next = cov->sub[0];
  return INIT(next, 0, s);
}

void spectralcox(cov_model *cov, gen_storage *s, double *e) { 
  cov_model *next = cov->sub[0];
  int d,
    dim = cov->tsdim - 1;
  double t, v[CoxMaxDim],
    *V = P(COX_MU),
    rho= P0(COX_D);
  SPECTRAL(next, s, e); // nicht gatternr
  
  v[0] = rnorm(0.0, INVSQRTTWO);
  v[1] = rho * v[0] + sqrt(1 - rho * rho) * rnorm(0.0, INVSQRTTWO);
 
  for (t = 0.0, d=0; d<dim; d++) {
    t += (v[d] + V[d]) * e[d];
  }
  e[dim] = -t;
}




// GenPS (generalisation of paciore and stein)
// Q = (x-y) Sy (Sx + Sy)^{-1} Sx (x-y) (weicht etwas von Stein ab)

#define STP_S 0
#define STP_Z 1
#define STP_M 2

#define STP_XI2 0
#define STP_PHI 1
#define STP_GAUSS 3
void kappa_stp(int i, cov_model *cov, int *nr, int *nc){
  *nc = (i == STP_S || i == STP_M) ? cov->tsdim : 1;
  *nr = i < CovList[cov->nr].kappas ? cov->tsdim : -1;
}

void stp(double *x,  double *y, cov_model *cov, double *v) {
  int d, j, k,
      dim =cov->tsdim,
      dimsq = dim * dim;
  double h[StpMaxDim], 
    Mh[StpMaxDim], hSx[StpMaxDim],
    Syh[StpMaxDim], xi2x, xi2y, 
    detA, hMh, cxy, zh, Q, Amux[StpMaxDim], Amuy[StpMaxDim],
    *Sc = P(STP_S),
    *M = P(STP_M),
    *z = P(STP_Z);
  cov_model *phi = cov->sub[STP_PHI],
    *Sf = cov->kappasub[STP_S],
    *xi2 =cov->sub[STP_XI2];

  ALLOC_EXTRA(Sx, dimsq);
  ALLOC_EXTRA2(Sy, dimsq);
  ALLOC_EXTRA3(A, dimsq);
    
  if (Sf != NULL) {
    FCTN(x, Sf, Sx); // symmetric, pos definite !!
    FCTN(y, Sf, Sy);
    //
    //    if (false) {
    //      int ii;
    //      printf("x=%f %f y=%f %f\n", x[0], x[1], y[0], y[1]);
    //for (ii=0; ii<4; ii++) printf("%f ", Sx[ii]);
    //      printf("\n");
    //      for (ii=0; ii<4; ii++) printf("%f ", Sy[ii]);
    //      printf("\n");
    //    }

  } else {
    int bytes = sizeof(double) * dimsq;
    MEMCOPY(Sx, Sc, bytes);
    MEMCOPY(Sy, Sc, bytes);
  }

  if (xi2 != NULL) {
    FCTN(x, xi2, &xi2x);
    FCTN(y, xi2, &xi2y);
  } else {
    xi2x = xi2y = 0.0;
  }

  for (k=0, d=0; d<dim; d++) {
    h[d] = x[d] - y[d];
  } 

  zh = hMh = 0.0;
  for (k=0, d=0; d<dim; d++) {
    Mh[d] = hSx[d] = Syh[d] = 0.0;
    for (j=0; j<dim; j++, k++) {
     Mh[d] += h[j] * M[k];
     hSx[d] += h[j] * Sx[k];
     Syh[d] += h[j] * Sy[k]; // uses that S is symmetric
    }
    zh += z[d] * h[d];
    hMh += Mh[d] * h[d];
  }
  cxy = xi2x - xi2y - zh;
  
  for (k=d=0; d<dim; d++) {
    for (j=0; j<dim; j++, k++) {
      A[k] = Sx[k] + Sy[k] + 4.0 * Mh[d] * Mh[j];
      assert(R_FINITE(A[k]));
    }
    Amux[d] = hSx[d] + 2.0 * (hMh + cxy) * Mh[d]; // uses that M is symmetric
    Amuy[d] = Syh[d] + 2.0 * (hMh - cxy) * Mh[d];
  }


  det_UpperInv(A, &detA, dim); // here

  Q = cxy * cxy - hMh * hMh + xUy(Amux, A, Amuy, dim);
  if (Q < 0.0) {
    PRINTF("x=%f,%f y=%f,%f detA=%f\n", 
	   x[0], x[1], y[0], y[1], detA);
    PRINTF("cxy=%4f hMh=%f Amux=%f A[0]=%f Amuy=%f\ndim=%d h=%f,%f hSx=%f,%f, xUy=%f Q=%f\n", 
	   cxy, hMh, Amux[0], A[0], Amuy[0], 
	   dim, h[0], h[1], hSx[0], hSx[1], xUy(Amux, A, Amuy, dim), Q);
    BUG;
  }

  Q = sqrt(Q);

  aux_covfct auxcf;
  if ((auxcf = CovList[phi->gatternr].aux_cov) != NULL) 
    auxcf(x, y, Q, phi, v);
  else 
    FCTN(&Q, phi, v);

  double 
    dx = detU(Sx, dim), 
    dy = detU(Sy, dim);
  
  *v *=  pow(2.0, 0.5 * double(dim)) * pow(dx * dy / (detA * detA), 0.25);
}



int checkstp(cov_model *cov){
  cov_model 
    *phi = cov->sub[STP_PHI],
    *Sf = cov->kappasub[STP_S],
    *xi2 =cov->sub[STP_XI2];
  int err;
  
  int dim = cov->tsdim;
 if (dim > StpMaxDim)
    SERR2("For technical reasons max. dimension for ave is %d. Got %d.", 
	  StpMaxDim, cov->xdimprev);

 if (PisNULL(STP_S) && Sf==NULL) {  // Sc
   if ((cov->px[STP_S] = EinheitsMatrix(dim)) == NULL) 
      return ERRORMEMORYALLOCATION;
    cov->ncol[STP_S] = cov->nrow[STP_S] = dim;
  }
 if (PisNULL(STP_M)) { // M
   if ((cov->px[STP_M] = EinheitsMatrix(dim)) == NULL)
      return ERRORMEMORYALLOCATION;
    cov->ncol[STP_M] = cov->nrow[STP_M] = dim;
  }
 if (PisNULL(STP_Z)) { // z
   PALLOC(STP_Z, dim, 1);
 }

 if (cov->xdimprev != cov->tsdim || cov->xdimprev != cov->tsdim)
    return ERRORDIM;
   
  if ((err = CHECK(phi, dim,  1, PosDefType, XONLY, ISOTROPIC,
		     SCALAR, cov->role // ROLE_COV changed 20.7.14 wg spectral
		   )) != NOERROR)
    return err;
  if (!isNormalMixture(phi->monotone)) return ERRORNORMALMIXTURE;

  cov->pref[Average] = 5;

  if (Sf != NULL) {
    if ((err = CHECK(Sf, dim,  dim, ShapeType, XONLY, CARTESIAN_COORD,
		       dim, cov->role // ROLE_COV changed 20.7.14 wg spectral
)) != NOERROR) 
      return err;
  }
  

  if (xi2 != NULL) {
   if ((err = CHECK(xi2, dim, dim,  ShapeType, XONLY, CARTESIAN_COORD,
		    SCALAR, cov->role // ROLE_COV changed 20.7.14 wg spectral
		    )) != NOERROR)
     return err;
  }

  // kein setbackward??

  EXTRA_STORAGE;

  cov->mpp.maxheights[0] = RF_NA;
  return NOERROR;
}

void rangestp(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  int i;
  for (i=0; i<=2; i++) { /* S, M, z */
    range->min[i] = RF_NEGINF;
    range->max[i]= RF_INF;
    range->pmin[i] = -1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}



int structStp(cov_model *cov, cov_model **newmodel) { 
  cov_model *shape;
  int err;
  
  ASSERT_NEWMODEL_NOT_NULL;
  
  if (cov->role != Average) ILLEGAL_ROLE;
  
  if ((err = covcpy(newmodel, cov)) != NOERROR) return err;
  shape = *newmodel;
  shape->nr = SHAPESTP;
  addModel(shape, STP_GAUSS, GAUSS);
  shape->sub[STP_GAUSS]->tsdim = 1;
  return NOERROR;
}


int check_shapestp(cov_model *cov) {  
  if (cov->sub[AVE_GAUSS] == NULL)
    SERR1("both submodels must be set to '%s'", CovList[GAUSS].nick);
  EXTRA_STORAGE;
  return checkstp(cov); // !! not next
}


int init_shapestp(cov_model *cov, gen_storage *s) {
  ASSERT_GAUSS_METHOD(Average);

  cov_model
    *Sf = cov->kappasub[STP_S],
    *gaussmix = cov->sub[STP_GAUSS];
  double sd,
    *q = cov->q;
  int err = NOERROR;

  if (Sf != NULL) {
    double minmax[2];
    assert(CovList[Sf->nr].minmaxeigenvalue != NULL);
    CovList[Sf->nr].minmaxeigenvalue(Sf, minmax);
    if (minmax[0] <= 0.0) error("neg eigenvalue in shape function of 'stp'");
    q[AVESTP_MINEIGEN] = minmax[0];
    q[AVESTP_LOGDET] = (double) cov->xdimprev * log(minmax[1]);
  } else {
#define dummyN (5 * StpMaxDim)
    double value[StpMaxDim], ivalue[StpMaxDim], dummy[dummyN], det,
      min = RF_INF;
    int i, Ferr,       
      dim = cov->tsdim,
      ndummy = dummyN;
 
    F77_NAME(dgeev)("No", "No", &dim, P(STP_S), &dim, 
		    value, ivalue, NULL, &dim, NULL, &dim,
		    dummy, &ndummy, &Ferr);
    if (Ferr != 0) SERR("error in F77 function call");
    det =  1.0;
    for (i=0; i<dim; i++) {
      double v = fabs(value[i]);
      det *= v;
      if (min > v) min = v;
    }
    q[AVESTP_MINEIGEN] = min;
    q[AVESTP_LOGDET] = log(det);
  }


  q[AVESTP_LOGV] = 0.0; 
  q[AVESTP_LOGMIXDENS] = 0.0;
  sd_avestp(cov, s, cov->tsdim, &sd); // sd->gauss

  assert(cov->vdim[0] == 1);  assert(cov->vdim[0] == cov->vdim[1]);

  if (cov->mpp.moments >= 0) {
    cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0; //// ??? notwendig 
    if (cov->mpp.moments >= 1) {
      if ((err = INIT(gaussmix, 2, s)) != NOERROR) return err;
      if (cov->mpp.moments >= 2) cov->mpp.mM[2] = 1.0; 
    }
  }
  return err;
}

void do_shapestp(cov_model *cov, gen_storage *s) {
  // Simulation of V; see Bernoulli Sec. 4.2
  cov_model 
    *stpGAUSS = cov->sub[STP_GAUSS],
    *phi = cov->sub[STP_PHI];
  double  spec_ret[StpMaxDim], sd,
    *q = cov->q;

  CovList[phi->nr].drawmix(phi, &(q[AVESTP_V]));
  sd_avestp(cov, s, cov->tsdim, &sd); // sd->gauss
  
  BUG;

  SPECTRAL(stpGAUSS, s, spec_ret);  // nicht gatternr
  q[AVERAGE_YFREQ] = *spec_ret * sqrt(q[AVESTP_V]);  
  q[AVERAGE_YPHASE] = TWOPI * UNIFORM_RANDOM;


  BUG; /// what to do with the next line?
  // info->logdens = CovList[phi->nr].logmixdens(ZERO, q[AVESTP_LOGV], phi);


  //info->radius = RF_INF;
  // info-sd s.o.
}


void logshapestp(double *x, double *u, cov_model *cov, double *v, double *Sign){
  // kann um ca. Faktor 2 beschleunigt werden, wenn
  // Sx , logdetU, Hx fuer alle x abgespeichert werden
  // und die Werte dann aus dem Speicher gelesen werden
  // jedoch sehr Speicherintensiv. MEMCOPY braucht man auch nicht
  cov_model 
    *Sf = cov->kappasub[STP_S],
    *xi2 =cov->sub[STP_XI2];
  int j, k, d, 
    dim= cov->xdimprev,
    dimsq = dim * dim,
    bytes = sizeof(double) * dimsq;
  double h[StpMaxDim], hSxh, hSx, xi, Mhd, 
    *Sc = P(STP_S),
    *M = P(STP_M),
    *z = P(STP_Z),
    *q = cov->q;
  
  ALLOC_EXTRA(Sx, dimsq); 

  if (Sf == NULL) {
    MEMCOPY(Sx, Sc, bytes);
  } else {
    FCTN(x, Sf, Sx); // symmetric, pos definite !!
  }

  if (xi2 == NULL) {
    xi = 0.0;
  } else {
    FCTN(x, xi2, &xi);
  }

  for (k=0, d=0; d<dim; d++) {
    h[d] = u[d] - x[d];
  }

  hSxh = 0.0;
  for (k=0, d=0; d<dim; d++) {
    Mhd = hSx = 0.0;
    for (j=0; j<dim; j++) {
     Mhd += h[j] * M[k];
     hSx += h[j] * Sx[k++];
     
    }
    xi += Mhd * h[d] + z[d] * h[d];
    hSxh += hSx * h[d];
  }
  
  double exponent =
    0.25 * dim * (// M_LN2 +  ??? !!! Rechnung!!! 
		  q[AVESTP_LOGV] - 2.0 * M_LN_SQRT_PI) // (2V/pi)^{d/4}
    + 0.25 *log(detU(Sx, dim))                          // Sx ^1/4 
      - q[AVESTP_V] * hSxh             // exp(-V(U-x) S (U-x)) 
    // + CovList[phi->nr].logmixdens(x, q[AVESTP_LOGV], phi) // g //nicht gatternr
    //    - 0.5 * cov_a->logdens // f 
    ;                // 1 / sqrt(f) 
  
  if (!(exponent < 5.0) && PL >= PL_DETAILS) {
    if (!(exponent < 6.0)) // could be NA, too
     PRINTF("\n%f logDetU=%f %f expon=%f",
          0.25 * dim * (// M_LN2 +  ??? !!! Rechnung!!! 
			q[AVESTP_LOGV] - 2.0 * M_LN_SQRT_PI) // (2V/pi)^{d/4}
	    , 0.25 * log(detU(Sx, dim))                         /// Sx ^1/4 
	    , -q[AVESTP_V]* hSxh             // exp(-V(U-x) S (U-x)) 
	    // , CovList[phi->nr].logmixdens(x, q[AVESTP_LOGV],  phi)// g
	    //, - 0.5 * cov_a->logdens // f 
	    , exponent);
    else PRINTF("!");
  };
  
  assert(exp(exponent) < 10000000.0);
  
 
  double cos_value = cos(q[AVERAGE_YPHASE] + q[AVERAGE_YFREQ] * xi);
  *v = exponent + log(fabs(cos_value)) ;  // Y 
  *Sign = cos_value > 0.0 ? 1.0 : cos_value < 0.0 ? -1.0 : 0.0;
}



/* Whittle-Matern or Whittle or Besset */ 
#define RATIONAL_A 0
#define RATIONAL_a 1
void kappa_rational(int i, cov_model *cov, int *nr, int *nc){
  *nc = (i == RATIONAL_A) ? cov->tsdim : 1;
  *nr = (i == RATIONAL_A) ? cov->tsdim : (i==RATIONAL_a) ? 2 : -1;
}
void minmaxEigenrational(cov_model *cov, double *mm) {
  double *a = P(RATIONAL_a);
  if (a[0] < a[1]) {
    mm[0] = a[0];
    mm[1] = a[1];
  } else {
    mm[0] = a[1];
    mm[1] = a[0];
  }
}
double maxEigenrational(cov_model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *mm) {
  double *a = P(RATIONAL_a);
  return (a[0] > a[1]) ? a[0] : a[1];
}
void rational(double *x, cov_model *cov, double *v) {
  int i, k, j, 
    dim = cov->tsdim;
  double nu,
    *A = P(RATIONAL_A),
    *a = P(RATIONAL_a);
  nu = 0.0;
  for (k=0, i=0; i<dim; i++) {
    double xTC;
    xTC =  0.0;
    for (j=0; j<dim; j++) {
      xTC += x[j] * A[k++];
    }
    nu += xTC * xTC;
  }
  *v = (a[0] + a[1] * nu) / (1 + nu);
}
 
int checkrational(cov_model *cov){
  int err;
  if (cov->nrow[RATIONAL_a] == 1) {
    double dummy = P0(RATIONAL_a);
    PFREE(RATIONAL_a);
    PALLOC(RATIONAL_a, 2, 1);
    P(RATIONAL_a)[0] = dummy;
    P(RATIONAL_a)[1] = 0.0;
  }
  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->mpp.maxheights[0] =  P(RATIONAL_a)[0] > P(RATIONAL_a)[1] 
    ? P(RATIONAL_a)[0] : P(RATIONAL_a)[1];
  return NOERROR;
}

void rangerational(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){

  range->min[RATIONAL_A] = RF_NEGINF;
  range->max[RATIONAL_A] = RF_INF;
  range->pmin[RATIONAL_A] = -1e10;
  range->pmax[RATIONAL_A] = 1e10;
  range->openmin[RATIONAL_A] = true;
  range->openmax[RATIONAL_A] = true;

  range->min[RATIONAL_a] = 0.0;
  range->max[RATIONAL_a] = RF_INF;
  range->pmin[RATIONAL_a] = 0.0;
  range->pmax[RATIONAL_a] = 10;
  range->openmin[RATIONAL_a] = false;
  range->openmax[RATIONAL_a] = true;
}


// Sigma(x) = diag>0 + A'xx'A
#define EAXXA_E 0
#define EAXXA_A 1
#define ETAXXA_ALPHA 2
void kappa_EAxxA(int i, cov_model *cov, int *nr, int *nc){
  *nc = (EAXXA_A == i) ? cov->tsdim : 1;
  *nr = i < CovList[cov->nr].kappas ? cov->tsdim : -1;
}
void EAxxA(double *x, cov_model *cov, double *v) {
  int d, k, j, 
    dim = cov->tsdim;
  double xA[EaxxaMaxDim],
    *E = P(EAXXA_E),
    *A = P(EAXXA_A);
  for (k=0, d=0; d<dim; d++) {
    xA[d] =  0.0;
    for (j=0; j<dim; j++) {
      xA[d] += x[j] * A[k++];
    }
  }
  for (k=d=0; d<dim; d++) {
    double xAd = xA[d];
    for (j=0; j<=d; j++) {
      v[k++] = xAd * xA[j];
    }
    v[k-1] += E[d];
    for ( ; j<dim; j++) {
      v[k++] = xAd * xA[j];
    }
  }
}

void minmaxEigenEAxxA(cov_model *cov, double *mm) {
  double 
    *E = P(EAXXA_E);
  int i,
    dim = cov->tsdim;
  for (mm[0] = RF_INF, mm[1]=RF_NEGINF, i=0; i<dim; i++) {
    if (E[i] < mm[0]) mm[0] = E[i];
    if (E[i] > mm[1]) mm[1] = E[i];
  }
}
 
int checkEAxxA(cov_model *cov){
  int err;
    if (cov->xdimown > EaxxaMaxDim)
      SERR2("For technical reasons max. dimension for ave is %d. Got %d.", 
	    EaxxaMaxDim, cov->xdimown);
    
  if ((err = checkkappas(cov)) != NOERROR) return err;

  cov->vdim[0] = cov->vdim[1] = cov->tsdim;
  cov->mpp.maxheights[0] = RF_NA;
 return NOERROR;
}
 
void rangeEAxxA(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){

  range->min[EAXXA_E] = 0.0;
  range->max[EAXXA_E] = RF_INF;
  range->pmin[EAXXA_E] = 0.0001;
  range->pmax[EAXXA_E] = 10;
  range->openmin[EAXXA_E] = true;
  range->openmax[EAXXA_E] = true;

  range->min[EAXXA_A] = RF_NEGINF;
  range->max[EAXXA_A] = RF_INF;
  range->pmin[EAXXA_A] = -1e10;
  range->pmax[EAXXA_A] = 1e10;
  range->openmin[EAXXA_A] = true;
  range->openmax[EAXXA_A] = true;
}





// Sigma(x) = diag>0 + A'xx'A
void kappa_EtAxxA(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  int tsdim = 3; 
  *nc = (i == EAXXA_A) ? tsdim : 1;
  *nr = (i == EAXXA_E || i==EAXXA_A) ? tsdim : (i==ETAXXA_ALPHA) ? 1 : -1;
}
void EtAxxA(double *x, cov_model *cov, double *v) {
  int d, k, j, 
    dim = cov->tsdim,
    time = dim - 1;
  double xAR[EaxxaMaxDim], R[9],
    *E = P(EAXXA_E),
    *A = P(EAXXA_A),
    phi = P0(ETAXXA_ALPHA),
    c =  cos(phi * x[time]),
    s = sin(phi * x[time]); 
     
  R[0] = R[4] = c;
  R[1] = s;
  R[3] = -s;
  R[2] = R[5] = R[6] = R[7] = 0.0;
  R[8] = 1.0;
 
  {
    double xA[EaxxaMaxDim];
    for (k=0, d=0; d<dim; d++) {
      xA[d] =  0.0;
      for (j=0; j<dim; j++) {
	xA[d] += x[j] * A[k++];
      }
    }
    for (k=0, d=0; d<dim; d++) {
      xAR[d] =  0.0;
      for (j=0; j<dim; j++) {
	xAR[d] += xA[j] * R[k++];
      }
    }
  }


  for (k=d=0; d<dim; d++) {
    double xAd = xAR[d];
    for (j=0; j<=d; j++) {
      v[k++] = xAd * xAR[j];
    }
    v[k-1] += E[d]; // nur korrekt falls E Vielfaches der EH-Matrix
    for ( ; j<dim; j++) {
      v[k++] = xAd * xAR[j];
    }
  }


}

void minmaxEigenEtAxxA(cov_model *cov, double *mm) {
  double 
    *E = P(EAXXA_E);
  int i,
    dim = cov->tsdim;
  for (mm[0] = RF_INF, mm[1]=RF_NEGINF, i=0; i<dim; i++) {
    if (E[i] < mm[0]) mm[0] = E[i];
    if (E[i] > mm[1]) mm[1] = E[i];
  }
}
 
int checkEtAxxA(cov_model *cov){
  int err;
  if (cov->xdimown != 3) SERR("The space-time dimension must be 3.");
  cov->vdim[0] = cov->vdim[1] = cov->tsdim;
  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->mpp.maxheights[0] = RF_NA;
 return NOERROR;
}
 
void rangeEtAxxA(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  int i;

  for (i=0; i<=2; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }

  range->min[EAXXA_E] = 0.0;
  range->max[EAXXA_E] = RF_INF;
  range->pmin[EAXXA_E] = 0.0001;
  range->pmax[EAXXA_E] = 10;
  range->openmin[EAXXA_E] = true;
  range->openmax[EAXXA_E] = true;
}




// Sigma(x) = diag>0 + A'xx'A
#define ROTAT_PHI 0 // both rotat and Rotat
#define ROTAT_SPEED 1
void kappa_rotat(int i, cov_model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = i < CovList[cov->nr].kappas ? 1 : -1;
}
void rotat(double *x, cov_model *cov, double *v) {
  int
    dim = cov->tsdim,
    time = dim - 1;
  double
    speed = P0(ROTAT_SPEED),
    phi = P0(ROTAT_PHI),
    absx = sqrt(x[0] * x[0] + x[1] * x[1]);
  *v = (absx == 0.0) ? 0.0
    : speed * (cos(phi * x[time]) * x[0] + sin(phi * x[time]) * x[1]) / absx;
}

void minmaxEigenrotat(cov_model VARIABLE_IS_NOT_USED *cov, double *mm) {
  mm[0] = -1;
  mm[1] = 1;
}
 
int checkrotat(cov_model *cov){
  int err;
  if (cov->xdimown != 3) SERR("The space-time dimension must be 3.");
   if ((err = checkkappas(cov)) != NOERROR) return err;
 cov->mpp.maxheights[0] = RF_NA;
  return NOERROR;
}
 
void rangerotat(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
  int i;

  for (i=0; i<2; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}


// Sigma(x) = diag>0 + A'xx'A
void kappa_Rotat(int i, cov_model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = i < CovList[cov->nr].kappas ?  1 : -1;
}
void Rotat(double *x, cov_model *cov, double *v) {
  int d, k, j, 
    dim = cov->tsdim,
    time = dim - 1;
  double
    phi = P0(ROTAT_PHI),
      c =  cos(phi * x[time]),
      s = sin(phi * x[time]),
      R[9]; assert(dim ==3);
   
  R[0] = R[4] = c;
  R[1] = s;
  R[3] = -s;
  R[2] = R[5] = R[6] = R[7] = 0.0;
  R[8] = 1.0;
 
  for (k=0, d=0; d<dim; d++) {
    v[d] =  0.0;
    for (j=0; j<dim; j++) {
      v[d] += x[j] * R[k++];
    }
  }
} 
int checkRotat(cov_model *cov){
  int err;
  if (cov->xdimown != 3) SERR("The space-time dimension must be 3.");
  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->vdim[0] = cov->vdim[1] = cov->tsdim;
  cov->mpp.maxheights[0] = RF_NA;
  return NOERROR;
}

void rangeRotat(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){
 
  range->min[ROTAT_PHI] = RF_NEGINF;
  range->max[ROTAT_PHI] = RF_INF;
  range->pmin[ROTAT_PHI] = -1e10;
  range->pmax[ROTAT_PHI] = 1e10;
  range->openmin[ROTAT_PHI] = true;
  range->openmax[ROTAT_PHI] = true;
}



/* Whittle-Matern or Whittle or Besset */ 
// statt 0 Parameter: 2 Parameter, M und z fuer xi
void kappaNonStWM(int i, cov_model *cov, int *nr, int *nc){
  *nc =  1;
  *nr = i < CovList[cov->nr].kappas ? 1 : -1;
}

void NonStWMQ(double *x, double *y, double sqrtQ, cov_model *cov, double *v){
// check calling functions, like hyperbolic and gneiting if any changings !!
  double loggamma, nuxy, nux, nuy;
  cov_model *nu = cov->kappasub[WM_NU];

  if (nu == NULL) {
    nuxy = P0(WM_NU);
    loggamma = lgammafn(nuxy);
  } else {
    FCTN(x, nu, &nux);
    FCTN(y, nu, &nuy);
    nuxy = 0.5 * (nux + nuy);
    loggamma = 0.5 * (lgammafn(nux) + lgammafn(nuy));
  }
    
  if (sqrtQ == 0.0) {
    *v = 1.0;
    return;
  }
 
  *v = 2.0 * exp(nuxy * log(0.5 * sqrtQ) - loggamma
		 + log(bessel_k(sqrtQ, nuxy, 2.0)) - sqrtQ);
}

void NonStWM(double *x, double *y, cov_model *cov, double *v){
// check calling functions, like hyperbolic and gneiting if any changings !!
  int d,
    dim = cov->tsdim;
  double Q=0.0;
  for (d=0; d<dim; d++) {
    double delta = x[d] - y[d];
    Q += delta * delta;
  }
  NonStWMQ(x, y, sqrt(Q), cov, v);
}

int checkNonStWM(cov_model *cov) { 
  cov_model *nu = cov->kappasub[WM_NU];
  int err,
    dim = cov->tsdim;

  NotProgrammedYet("");
 

  if (PisNULL(WM_NU) && nu==NULL) SERR1("'%s' is missing", KNAME(WM_NU));
  if (isRandom((Types) CovList[cov->nr].kappaParamType[WM_NU])) 
    SERR1("only deterministic models for '%s' are allowed.\nHowever these models can have random parameters.", KNAME(WM_NU));
  
  if (nu != NULL) {
    if ((err = CHECK(nu, dim, dim, ShapeType, XONLY, CARTESIAN_COORD,
		       SCALAR, cov->role // ROLE_COV changed 20.7.14 wg spectral
		     )) != NOERROR) 
      return err;
    if (nu->tsdim != cov->tsdim) return ERRORWRONGDIM;
  }
  // no setbackard !!
  return NOERROR;
}

sortsofparam paramtype_nonstWM(int VARIABLE_IS_NOT_USED  k, int VARIABLE_IS_NOT_USED  row, int VARIABLE_IS_NOT_USED col) {
  return CRITICALPARAM;
}


void rangeNonStWM(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range){ 
  range->min[WM_NU] = 0.0;
  range->max[WM_NU] = RF_INF;
  range->pmin[WM_NU] = 1e-2;
  range->pmax[WM_NU] = 10.0;
  range->openmin[WM_NU] = true;
  range->openmax[WM_NU] = true;
}

// using nu^(-1-nu+a)/2 for g and v^-a e^{-1/4v} as density instead of frechet
// the bound 1/3 can be dropped
// static double eM025 = exp(-0.25);
//void DrawMixNonStWM(cov_model *cov, double *random) { // inv scale
//  // V ~ F in stp
//  cov_model *nu = cov->sub[WM_NU];  
//  double minnu;
//  double alpha;
//
//  if (nu == NULL) {
//    minnu = P(WM_NU][0];
//  } else {
//    double minmax[2];
//    CovList[nu->nr].minmaxeigenvalue(nu, minmax);
//    minnu = minmax[0];
//  }
//  alpha = 1.0 + 0.5 /* 0< . < 1*/ * (3.0 * minnu - 0.5 * cov->tsdim);
//  if (alpha > 2.0) alpha = 2.0; // original choice
//  if (alpha <= 1.0) ERR("minimal nu too low or dimension too high");
//
// error("logmixdensnonstwm not programmed yet");
 /* 
  double beta = GLOBAL.mpp.beta,
    p = GLOBAL.mpp.p,
    logU;
  if (UNIFORM_RANDOM < p){
    cov_a->WMalpha = beta;
    logU =  log(UNIFORM_RANDOM * eM025);
    cov_a->WMfactor = -0.5 * log(0.25 * p * (beta - 1.0)) + 0.25;
  } else {
    cov_a->WMalpha = alpha;
    logU = log(eM025 + UNIFORM_RANDOM * (1.0 - eM025));
    cov_a->WMfactor = -0.5 * log(0.25 * (1.0 - p) * (alpha - 1.0));
  } 
  
  logmix!!

  *random = log(-0.25 / logU) / (cov_a->WMalpha - 1.0); //=invscale
  */
//}




//double LogMixDensNonStWM(double *x, double logV, cov_model *cov) {
//  // g(v,x) in stp
//  double z = 0.0;
//  error("logmixdensnonstwm not programmed yet");
  // wmfactor ist kompletter unsinn; die 2 Teildichten muessen addiert werden
  /*
  cov_model *calling = cov->calling,
    *Nu = cov->sub[0];
  if (calling == NULL) BUG;
   double nu,
    alpha = cov_a->WMalpha,
    logV = cov_a->logV,
    V = cov_a->V;
  
  if (Nu == NULL) 
    nu = P(WM_NU][0];
  else 
     FCTN(x, Nu, &nu);


   z = - nu  * M_LN2 // in g0  // eine 2 kuerzt sich raus
    + 0.5 * ((1.0 - nu) // in g0
    + alpha // lambda
    - 2.0 //fre*
    ) * logV
    - 0.5 * lgammafn(nu)  // in g0
    + cov_a->WMfactor // lambda
    - 0.125 / V   // g: Frechet
    + 0.125 * pow(V, 1.0 - alpha); // lambda: frechet

  if (!(z < 7.0)) {
    static double storage = 0.0; 
    if (gen_storage != logV) {
      if (PL >= PL_DETAILS) 
	PRINTF("alpha=%f, is=%f, cnst=%f logi=%f lgam=%f loga=%f invlogs=%f pow=%f z=%f\n",
	       alpha,V,
	       (1.0 - nu) * M_LN2 
	       , + ((1.0 - nu) * 0.5 + alpha - 2.0) * logV
	       ,- 0.5 * lgammafn(nu) 
	       , -cov_a->WMfactor
	       ,- 0.25 / V 
	       , + 0.25 * pow(V, - alpha)
	       , z);
      storage = logV;
    }
    //assert(z < 10.0);
  }
*/
//  return z;
//				      
//}




/* nsst */
/* Tilmann Gneiting's space time models, part I */
#define NSST_DELTA 0
void nsst(double *x, cov_model *cov, double *v) {
  cov_model *subphi = cov->sub[0];
  cov_model *subpsi = cov->sub[1];
  double v1, v2, psi, y;
  
  COV(ZERO, subpsi, &v1);
  COV(x + 1, subpsi, &v2);
  psi = sqrt(1.0 + v1 - v2);  // C0 : C(0) oder 0 // Cx : C(x) oder -gamma(x)
  y = x[0] / psi;
  COV(&y, subphi, v);
  *v *= pow(psi, -P0(NSST_DELTA));
}

void Dnsst(double *x, cov_model *cov, double *v) {
  cov_model *subphi = cov->sub[0];
  cov_model *subpsi = cov->sub[1];
  double v1, v2, psi, y;

  COV(ZERO, subpsi, &v1);
  COV(x + 1, subpsi, &v2);
  psi = sqrt(1.0 + v1 - v2);  // C0 : C(0) oder 0 // Cx : C(x) oder -gamma(x)
  y = x[0] / psi;
  Abl1(&y, subphi, v);
  *v *= pow(psi, -P0(NSST_DELTA) - 1.0);
}

void TBM2nsst(double *x, cov_model *cov, double *v) {
  cov_model *subphi = cov->sub[0];
  cov_model *subpsi = cov->sub[1];
  double v1, v2, psi, y;

  assert(false);

  COV(ZERO, subpsi, &v1);
  COV(x + 1, subpsi, &v2);
  psi = sqrt(1.0 + v1 - v2);  // C0 : C(0) oder 0 // Cx : C(x) oder -gamma(x)
  y = x[0] / psi;
  TBM2CALL(&y, subphi, v);
  *v *= pow(psi, -P0(NSST_DELTA));
}

int checknsst(cov_model *cov) {
  cov_model *subphi = cov->sub[0];
  cov_model *subpsi = cov->sub[1];
  int err;
      
  if (cov->xdimown != 2) SERR("reduced dimension must be 2");

  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->finiterange = false;
  if ((err = CHECK(subphi, cov->tsdim, 1, PosDefType, XONLY, ISOTROPIC, 
		     SCALAR, cov->role // ROLE_COV changed 20.7.14 wg spectral
		   )) != NOERROR) 
    return err;
  
  if (!isNormalMixture(subphi->monotone)) return(ERRORNORMALMIXTURE);
  setbackward(cov, subphi);
  assert(cov->finiterange != true); //

  if ((err = CHECK(subpsi, 1, 1, VariogramType, XONLY, ISOTROPIC, 
		     SCALAR, cov->role // ROLE_COV changed 20.7.14 wg spectral
		   )) != NOERROR) 
    return err;

  subphi->delflag = subpsi->delflag = DEL_COV-20;

  // kein setbackward(cov, subpsi);
  return NOERROR;
}

sortsofparam paramtype_nsst(int k, int VARIABLE_IS_NOT_USED row, int VARIABLE_IS_NOT_USED col) {
  return k==-1 ? VARPARAM : k==0 ? CRITICALPARAM : ANYPARAM;
}

void rangensst(cov_model *cov, range_type* range){ 
  range->min[NSST_DELTA] = cov->tsdim - 1;
  range->max[NSST_DELTA] = RF_INF;
  range->pmin[NSST_DELTA] = cov->tsdim - 1;
  range->pmax[NSST_DELTA] = 10.0;
  range->openmin[NSST_DELTA] = false;
  range->openmax[NSST_DELTA] = true;
}




/* gennsst */
/* generalisation in schlather, bernoulli 2010 */

#define GENNSST_INTERN_A 0

void kappa_gennsst_intern(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  *nc = *nr = i == 0 ? 0 : -1;
}

void nonstatgennsst_intern(double *x,double *y, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double  det, z,
    *A = P(GENNSST_INTERN_A);
  int d, 
    dim = cov->tsdim;
  ALLOC_EXTRA(ds, dim);
 
  //  printf("\nx=%f %f %f \ny=%f %f %f\nz=%f %f %f\n %f %f %f\n %f %f %f\n",  	 x[0], x[1], x[2], y[0],y[1],y[2],	 A[0],  A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8]);

  //  PMI(cov->calling->sub[1]);

  det_UpperInv(A, &det, dim);
  for (d=0; d<dim; d++) ds[d] = x[d] - y[d];
  z = sqrt(xUx(ds, A, dim));
  COV(&z, next, v);
  *v /= sqrt(det);
}

int checkgennsst_intern(cov_model *cov) {
  // internes modell muss zwischengeschaltet werden
  // um Earth->Cartesian->isotropic hinzubekommen
  // und         ^- hier noch mit A zu multiplizieren
  cov_model *next = cov->sub[0];
  int err,
    dim = cov->xdimown;

  

  assert(cov->tsdim == cov->xdimown);
  if ((err = CHECK(next, cov->tsdim, 1, PosDefType, 
		   XONLY, ISOTROPIC, SCALAR, cov->role )) != NOERROR) 
    return err;
  if (!isNormalMixture(next->monotone)) return ERRORNORMALMIXTURE;

  //APMI(cov->calling);

  if (PisNULL(GENNSST_INTERN_A)) {  
    PALLOC(GENNSST_INTERN_A, dim, dim);
  } else if (cov->nrow[GENNSST_INTERN_A] != dim) {
    PFREE(GENNSST_INTERN_A);
    PALLOC(GENNSST_INTERN_A, dim, dim);
  }

  cov->finiterange = false;
  setbackward(cov, next);
  assert(!cov->finiterange); 
  cov->vdim[0] = cov->vdim[1] = 1;

  EXTRA_STORAGE;
  return NOERROR;
}


void range_gennsst_intern(cov_model VARIABLE_IS_NOT_USED *cov,
			  range_type* range){
  range->min[GENNSST_INTERN_A] = RF_NEGINF;
  range->max[GENNSST_INTERN_A] = RF_INF;
  range->pmin[GENNSST_INTERN_A] = -100.0;
  range->pmax[GENNSST_INTERN_A] = +100.0;
  range->openmin[GENNSST_INTERN_A] = true;
  range->openmax[GENNSST_INTERN_A] = true;
}


#define GENNSST_DELTA 0
void nonstatgennsst(double *x,double *y, cov_model *cov, double *v) {
  cov_model *subphi = cov->key;
  cov_model *subpsi = cov->sub[1];
  int i, 
    dim = cov->tsdim,
    dimSq = dim * dim;
    
  //  print("A %s\n", TYPENAMES[subpsi->typus]);
  NONSTATCOV(x, y, subpsi, PARAM(subphi, GENNSST_INTERN_A));
  if (isVariogram(subpsi)) {
    ALLOC_EXTRA(z, dimSq);
    NONSTATCOV(ZERO, ZERO, subpsi, z);
    for (i=0; i<dimSq; i++) 
      PARAM(subphi, GENNSST_INTERN_A)[i] = 
	z[i] - PARAM(subphi, GENNSST_INTERN_A)[i]; 
  } else if (subpsi->typus != NegDefType) BUG;


  //    print("%f %f %f %f  dim = %d %f %f %f %f %s\n", x[0], x[1], y[0], y[1], dim, v1[0], v1[1], v1[2], v1[3],NAME(subpsi));

  NONSTATCOV(x, y, subphi, v);
}
  

int checkgennsst(cov_model *cov) {
  cov_model *subphi = cov->sub[0];
  cov_model *subpsi = cov->sub[1];
#define ncs 3
  int err,
    UCS = UpgradeToCoordinateSystem(cov->isoown);
   assert(cov->tsdim == cov->xdimown);

   if (cov->q == NULL) {
     QALLOC(1);
     cov->q[0]=NOERROR;
   }

   if (isSpherical(cov->isoown)) 
     return cov->q[0]!=NOERROR ? (int) cov->q[0] : ERRORFAILED;
   if (UCS == ISO_MISMATCH) 
     return cov->q[0]!=NOERROR ? (int) cov->q[0] : ERRORODDCOORDTRAFO;
   if (cov->xdimown != cov->tsdim) SERR("logical and physical dimension differ");
 
   //STOPAFTER(cov->isoown==EARTH_COORD || true, PMI(cov));
   //   assert(cov->isoown ==EARTH_ISOTROPIC || __stopafter__ ==1);
 
   if (cov->key == NULL) {
     if ((err = covcpy(&(cov->key), subphi)) != NOERROR) return err;
     addModel(&(cov->key), GENNSST_INTERN);
   }

  if ((cov->q[0] =err = CHECK(cov->key, cov->tsdim, cov->tsdim, PosDefType, 
			      KERNEL, SYMMETRIC,  // viel zu schwach, geht aber 
			      //                     z.Zt. nicht anders
			      // hinsichtlich der Koordinatentransformation
			      // -- wird das interne modell abgefangen
			      SCALAR, cov->role )) != NOERROR) {
    return err;
  }

  //  PL = GLOBAL.general.Cprintlevel = 7;

  
  // MUSSS ZWINGEND ALS ZWEITES KOMMEN
  if ((err = CHECK(subpsi, cov->tsdim, cov->tsdim, NegDefType,
		   KERNEL, UCS, 
		   cov->key->xdimown, cov->role)) != NOERROR) {
     return err;
  }
 

  cov->finiterange = false;
  setbackward(cov, cov->key);
  assert(!cov->finiterange); //

  cov->vdim[0] = cov->vdim[1] = 1;
  
  EXTRA_STORAGE;

  COV_DELETE(cov->sub + 0);
  if ((err = covcpy(cov->sub + 0, cov->key->sub[0])) != NOERROR) BUG;
  cov->sub[0]->calling = cov;

  // PMI(cov->calling);   assert(err == NOERROR);BUG;
  // ERR(err); BUG;

 
  // kein setbackward(cov, subpsi);
  return NOERROR;
}

