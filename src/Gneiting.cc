/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Gneiting's space-time covariance models and related models

 Copyright (C) 2006 -- 2017 Martin Schlather

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
#include "operator.h"
#include "shape.h"
#include "Processes.h"
#include "startGetNset.h"

#define LOG2 M_LN2


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


void kappa_ave(int i, model *cov, int *nr, int *nc){
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
 int dim = OWNLOGDIM(0);
  if (spacetime) dim--;
  *nr = (i==AVE_A || i==AVE_Z) ? dim : 1;
  *nc = (i==AVE_A) ? dim : i < DefList[COVNR].kappas ? 1 : -1;
}


void ave(double *h, model *cov, double *v) {
  // f = uAu +zu; 
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  model *next = cov->sub[0];
  int i,j,k,d,
    dim = OWNLOGDIM(0);
  if (spacetime) dim--;

  double Ah[AveMaxDim], Eplus2B[AveMaxDim], 
    dummy, 
    hh,
    *A = P(AVE_A),
    *z = P(AVE_Z),
    c = spacetime ? h[dim] : 0.0; // Sockelwert fuer c

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

  double det, AhEinvAh;
  Ext_XCinvXdet(Eplus2B, dim, Ah, 1, &AhEinvAh, &det, false, NULL);
  double y = SQRT(0.5 * hh  + c * c * (1.0 - 2.0 * AhEinvAh));
  COV(&y, next, v);
  *v /= SQRT(det);
}


int checkave(model *cov) {
   ASSERT_UNREDUCED;
 model *next = cov->sub[0];
  bool
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  int i, j, err,
   dim =  OWNLOGDIM(0),
    spdim = spacetime ? dim - 1 : dim;
  double 
    *A = P(AVE_A);
  char dim_info[2][4] = {"d", "d-1"};

  if (OWNTOTALXDIM < 2) SERR("The spatial dimension must be at least 2.");
 
  if (dim > AveMaxDim)
    SERR2("For technical reasons max. dimension for ave is %d. Got %d.", 
	  AveMaxDim, dim);

  if (cov->ncol[AVE_A] != spdim || cov->nrow[AVE_A] != spdim) 
    SERR5("A not %.50sx%.50s matrix, but %dx%d (dim=%d)", dim_info[spacetime], 
	  dim_info[spacetime], cov->ncol[AVE_A], cov->nrow[AVE_A], spdim);

  if (cov->ncol[AVE_Z] != 1 || cov->nrow[AVE_Z] != spdim) 
    SERR1("z not (%.50s)-dim vector", dim_info[spacetime]);

  for (i=0; i<spdim; i++)
    for (j=i+1; j<spdim; j++)
      if (A[i + j * spdim] != A[j + i * spdim]) {
	A[j + i * spdim] = A[i + j * spdim];
	warning("A is not symmetric -- lower part used");
      }
  // naechste Zeile stimmt nicht mit Bernoulli ueberein

  kdefault(cov, AVE_SPACETIME, true);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);

  if ((err = CHECK(next, dim, 1, PosDefType, XONLY, ISOTROPIC,
		   SCALAR, cov->frame // VariogramType changed 20.7.14 wg spectral
		   )) != NOERROR) RETURN_ERR(err);
  if (!isNormalMixture(next->monotone)) RETURN_ERR(ERRORNORMALMIXTURE);
  if (DefList[NEXTNR].spectral == NULL) RETURN_ERR(ERRORSPECTRAL); // nicht gatter

  //  updatepref(cov, next); ## gute idee?
  if (next->pref[SpectralTBM] == PREF_NONE) 
    cov->pref[RandomCoin] = cov->pref[Average] = PREF_NONE;

  // kein setbackward??
  
  // no setbackard
  RETURN_NOERROR;
}



void rangeave(model VARIABLE_IS_NOT_USED *cov, range_type* ra){
  int i;
  
  for (i=0; i<=1; i++) {
    ra->min[i] = RF_NEGINF;
    ra->max[i] = RF_INF;
    ra->pmin[i] = - 10.0;
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



void sd_avestp(model *cov, gen_storage VARIABLE_IS_NOT_USED *S, int dim, double *sd){
  /////  window_info *w = &(S->window);
  int d;
  double b, alphamin, x2, InvSqrt2a, EmA,
    *q = cov->q;
  // see article/GEOSTATS/simuspacetime/simuspacetime2008/simuspacetime.tex 
  // for the reasoning of these calculations

  BUG; // auf keinen fall parallel tauglich

  q[AVESTP_LOGV] = LOG(q[AVESTP_V]);
  for (x2=0.0, d=0; d<dim; d++) {
    double lensimu = RF_NA; ////  w->max[d] - w->min[d];
    x2 += lensimu * lensimu;
  }
  // x2 *= 0.25;
  b = 3.0 * q[AVESTP_V] * x2 / dim;
  alphamin = (4.0 + 4.0 * b - 2.0 * SQRT(4.0 * b * b + 8.0 * b + 1.0)) / 3.0;
  InvSqrt2a = 1.0 / SQRT(2.0 * alphamin * 6.0 * q[AVESTP_V]);
  *sd = InvSqrt2a;
  EmA = 1.0 - alphamin;
  cov->mpp.maxheights[0] = EXP(-0.5 * LOG(EmA) - 0.25 * LOG(alphamin) + b / EmA -
			   2 * x2); // proportional zum dritten Moment !

  /*
   double radius = 
    SQRT((-9 // so e^{-9} as threshold
	  - 0.25 * dim * (q[AVESTP_LOGV] - 1.14473) // log pi
	  - 0.25 * q[AVESTP_LOGDET]
	  //+ 0.5 * cov_a->logdens
	  -  q[AVESTP_LOGMIXDENS]
	  ) / ( - q[AVESTP_MINEIGEN] * q[AVESTP_V]) ); // ???
  assert(radius > 0);
  */

}


int structAve(model *cov, model **newmodel) { 
  model *shape; //, *gaussmix=NULL;
  int err;
 
  ASSERT_NEWMODEL_NOT_NULL;
   ERR0("'ave' currently does not work"); // what to do with the next lines??
 

  if ((err = covcpy(newmodel, cov)) != NOERROR) RETURN_ERR(err);
  shape = *newmodel;
  SET_NR(shape, SHAPEAVE);
  addModel(shape, AVE_GAUSS, GAUSS);
  
  ERR0("'ave' currently does not work");
// what to do  ERR1("'ave' currently does not work (% ld)", (long int) gaussmix); // what to do with the next lines??
  /*
  gaussmix = shape->sub[AVE_GAUSS];
  gaussmix->ts dim = 1;
  gaussmix->frame = GaussMethodType;
  gaussmix->method = SpectralTBM;
  */
   RETURN_NOERROR;
}




void  logshapeave(double *x, model *cov, double *v, double *Sign) {
    // nur stationaer
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  int d, j, k,
    dim = OWNLOGDIM(0);
  if (spacetime) dim--;
  double f, dummy, r2,
    *A = P(AVE_A),
    *z = P(AVE_Z),
    t = spacetime ? x[dim] : 0.0,
    *q = cov->q;
 
  f = r2 = 0.0;
  for (k=d=0; d<dim; d++) {
    r2 += x[d] * x[d];
    for (dummy = z[d], j=0; j<dim; j++) dummy += x[j] * A[k++];
    f += dummy * x[d];
  }

#ifdef SCHLATHERS_MACHINE
  warning("is exponent of V correct?"); 
#endif  
   
  v[0] = 0.25 * dim * q[AVESTP_LOGV] // V^{d/4 oder d/2} !!!!!!!!!!!
    - 0.5 * (LOG2 - dim * M_LN_SQRT_PId2) - r2
    // LOG2 : wegen spectral  simulation; 
    // zweiter Term fuer logg 
    //+ DefList[phi->nr].logmixdens(x, q[AVESTP_LOGV], phi); /* g */// nicht gatternr
    ;
  Sign[0] = 1.0;
  double phase = q[AVERAGE_YPHASE] + q[AVERAGE_YFREQ] * (f - t); // Y
  Sign[1] =  phase > 0.0 ? 1.0 : phase < 0.0 ? -1.0 : 0.0;
  v[1] = LOG(FABS(phase));
}

int check_shapeave(model *cov) {
  if (cov->sub[AVE_GAUSS] == NULL)
    SERR1("both submodels must be set to '%.50s'", DefList[GAUSS].nick);
  cov->mpp.maxheights[0] = RF_NA;
  return checkave(cov); // !! not next
}

int init_shapeave(model *cov, gen_storage *s) { 
  ASSERT_GAUSS_METHOD(Average);
  model
    *gaussmix = cov->sub[AVE_GAUSS];
  double sd,
    *q = cov->q;
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  int err = NOERROR,
   dim = OWNLOGDIM(0);
  if (spacetime) dim--;
  

  q[AVESTP_V] = 0.0;
  q[AVESTP_MINEIGEN] = 1.0; 
  q[AVESTP_LOGDET] = 0.0;
  sd_avestp(cov, s, dim, &sd); // sd->gauss

  assert(VDIM0 == 1);  assert(VDIM0 == VDIM1);

  if (cov->mpp.moments >= 0) {
    cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0; 
    if (cov->mpp.moments >= 1) {
      if ((err = INIT(gaussmix, cov->mpp.moments, s)) != NOERROR) RETURN_ERR(err);
      if (cov->mpp.moments >= 2) {
	cov->mpp.mM[2] = 1.0;
      }
    }
  }
  RETURN_ERR(err);
}


void do_shapeave(model *cov, gen_storage *S) { 
  // Simulation of V; sopee Bernoulli Sec. 4.2
   model    
    *aveGAUSS = cov->sub[AVE_GAUSS],
    *phi = cov->sub[AVE_PHI];
  double spec_ret[StpMaxDim], sd,
    *q = cov->q;
  bool 
    spacetime = (bool) (PisNULL(AVE_SPACETIME) || P0INT(AVE_SPACETIME));
  int 
    dim = OWNLOGDIM(0);
  if (spacetime) dim--;
   
  DefList[MODELNR(phi)].drawmix(phi, q + AVESTP_V); // nicht gatternr
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
// EXP(- 0.5 * (h *h + 2 a^2 - mu C mu)) // stimmen die Vorzeichen??
// mu = h - 2 a M h
/* cox, cmp with nsst1 !! */
// coxisham
#define COX_MU 0
#define COX_D 1
#define COX_BETA 2

void GetEu2Dinv(model *cov, double *x, int dim, 
		double *det, double *Eu2Dinv,
		double *newxsq, double *newx, double *z) {
    double 
	y[CoxMaxDim],
      *V = P(COX_MU),
      *D= P(COX_D),
      beta = P0(COX_BETA),
      t = x[dim],
      t2 = POW(FABS(t), beta); // standard t^2
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

  solve_storage *PT =  NULL;
  if (z !=NULL) PT->result = z;
  int err = Ext_XCinvXdet(Eu2Dinv, dim, y, 1, newxsq, det, false, PT);
  if (z != NULL) PT->result = NULL;
  if (err != NOERROR) ERR0("error occuredin 'GetEu2Dinv'");
  *newx = SQRT(*newxsq);
}

void cpyUf(double *Eu2Dinv, double factor, int dim, int logicaldim, double *v) {
    // Eu2Dinv has dimension dim^2; v dimension logicaldim^2
    // Eu2Dinv is copied into the upper left part of v and 
    // multiplied by factor
  int d, i, k, j,
	logicaldimsq = logicaldim * logicaldim;
  
  for (i=0; i<logicaldimsq; v[i++] = 0.0);
  for (i=0; i<dim; i++) {
      for (d=i * logicaldim, k=i * dim, j=0; j<=i; j++)
	  v[d++] = Eu2Dinv[k++] * factor; 
      for (k += dim - 1; j<dim; j++, k+=dim) { 
	  v[d++] = Eu2Dinv[k] * factor;
      }
  }
}

void addzzT(double *v, double factor, double *z, int dim, int logicaldim) {
    // z has dimension dim; v dimension logicaldim^2
    // zzT is copied into the upper left part of v after being 
    // multiplied by factor
   
    int i,j,k;
    for (i=0; i<dim; i++) {
	k = i * logicaldim;
	for (j=0; j<dim; j++) {
	    v[k++] += z[i] * z[j] * factor;
	}
    }
}


void kappa_cox(int i, model *cov, int *nr, int *nc){
    switch (i) {
	case COX_MU :
	    *nc = 1;
	    *nr = OWNLOGDIM(0)  - 1;
	    break;
	case COX_D  :
	    *nc = *nr = OWNLOGDIM(0)  - 1;
	    break;
	case COX_BETA :
	    *nc = *nr = 1;
	    break;
	default:  *nc = *nr = OUT_OF_RANGE;
    }
}

void cox(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
 int dim = OWNLOGDIM(0) - 1,
      dimsq = dim * dim;
  double det, newx,  newxsq;
  TALLOC_XX1(Eu2Dinv, dimsq);
  GetEu2Dinv(cov, x, dim, &det, Eu2Dinv, &newxsq, &newx, NULL);
  END_TALLOC_XX1;
  // printf("cox det=%10g %10g\n", det, newx);
  
  COV(&newx, next, v);
  *v /= SQRT(det);
 } 

void coxhess(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  int logicaldim = OWNLOGDIM(0),
      dim = logicaldim - 1,
      dimsq = dim * dim;
  double z[CoxMaxDim], det, newx, newxsq, phiD, phiD2;

  TALLOC_XX1(Eu2Dinv, dimsq);

  GetEu2Dinv(cov, x, dim, &det, Eu2Dinv, &newxsq, &newx, z);

  Abl2(&newx, next, &phiD2);  
  if (newxsq == 0.0) {
    cpyUf(Eu2Dinv, phiD2 / SQRT(det), dim, logicaldim, v);
  } else {
    Abl1(&newx, next, &phiD);
    cpyUf(Eu2Dinv, phiD / (SQRT(det) * newx), dim, logicaldim, v);
    addzzT(v, (phiD2 - phiD/newx) / (SQRT(det) * newxsq), z, dim, logicaldim);
  }
  END_TALLOC_XX1;
} 
 

void coxnabla(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  int d,
    logicaldim = OWNLOGDIM(0),
    dim = logicaldim - 1,
    dimsq=dim * dim;
  double z[CoxMaxDim], det, newx, newxsq, phiD, factor;
  
  TALLOC_XX1(Eu2Dinv, dimsq);
  GetEu2Dinv(cov, x, dim, &det, Eu2Dinv, &newxsq, &newx, z); 
  END_TALLOC_XX1;

  if (newxsq == 0.0) {
      for (d=0; d<=dim; d++)  v[d] = 0.0;
  } else {    
    newx = SQRT(newxsq);
    Abl1(&newx, next, &phiD);
    factor = phiD / (det * newx); 
    for (d=0; d<dim; d++) {
	v[d] = factor * z[d];
    }
    for (d=0; d<logicaldim; v[d++]=0.0);
  }
 }



int checkcox(model *cov) {
  model *next = cov->sub[0];
  int err, i, 
    logicaldim = OWNLOGDIM(0),
    dim = logicaldim - 1,
    dimsq = dim * dim; 

  if (OWNTOTALXDIM < 2) SERR("The space-time dimension must be at least 2.");  
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
    if (!Ext_is_positive_definite(P(COX_D), dim))
      SERR("D is not (strictly) positive definite");
  }

  kdefault(cov, COX_BETA, 2.0);

  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  if ((err = CHECK(next, dim,  1, PosDefType, XONLY, ISOTROPIC,
		     SCALAR, cov->frame // VariogramType changed 20.7.14 wg spectral
		   )) != NOERROR)
    RETURN_ERR(err);

  if (logicaldim != 3)  cov->pref[SpectralTBM] = PREF_NONE;

  if (!isNormalMixture(next->monotone)) RETURN_ERR(ERRORNORMALMIXTURE);
  if (DefList[NEXTNR].spectral == NULL) RETURN_ERR(ERRORSPECTRAL); // nicht gatternr
  
  // no setbackard
  updatepref(cov, next);
  if (P0(COX_BETA) != 2.0) cov->pref[SpectralTBM] = 0;

  cov->hess = true;
  if (cov->Ssolve == NULL) {EXT_NEW_STORAGE(solve);}
  
  EXTRA_STORAGE;	 
  RETURN_NOERROR;
} 

void rangecox(model VARIABLE_IS_NOT_USED *cov, range_type* range){

  range->min[COX_MU] = RF_NEGINF;
  range->max[COX_MU] = RF_INF;
  range->pmin[COX_MU] = - 100.0;
  range->pmax[COX_MU] = +100.0;
  range->openmin[COX_MU] = true;
  range->openmax[COX_MU] = true;

  range->min[COX_D] = RF_NEGINF;
  range->max[COX_D] = RF_INF;
  range->pmin[COX_D] = - 100.0;
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

int initcox(model *cov, gen_storage *s) {
  ASSERT_GAUSS_METHOD(SpectralTBM);
  model *next = cov->sub[0];
  return INIT(next, 0, s);
}

void spectralcox(model *cov, gen_storage *s, double *e) { 
  model *next = cov->sub[0];
  int d,
    logicaldim = OWNLOGDIM(0),
    dim = logicaldim - 1;
  double t, v[CoxMaxDim],
    *V = P(COX_MU),
    rho= P0(COX_D);
  SPECTRAL(next, s, e); // nicht gatternr
  
  v[0] = rnorm(0.0, INVSQRTTWO);
  v[1] = rho * v[0] + SQRT(1.0 - rho * rho) * rnorm(0.0, INVSQRTTWO);
 
  for (t = 0.0, d=0; d<dim; d++) {
    t += (v[d] + V[d]) * e[d];
  }
  e[dim] = -t;
}




// GenPS (generalisation of paciore and stein)
// Q = (x-y) Sy (Sx + Sy)^{-1} Sx (x-y) (weicht etwas von Stein ab)


#define STP_XI2 0
#define STP_PHI 1
#define STP_GAUSS 3
void kappa_stp(int i, model *cov, int *nr, int *nc){
  *nc = (i == STP_S || i == STP_M) ? OWNLOGDIM(0) : 1;
  *nr = i < DefList[COVNR].kappas ? OWNLOGDIM(0) : -1;
}

void stp(double *x,  double *y, model *cov, double *v) {
  int d, j, k, err,
    dim =  OWNLOGDIM(0),
    dimsq = dim * dim;
  double h[StpMaxDim], 
    Mh[StpMaxDim], hSx[StpMaxDim],
    Syh[StpMaxDim], xi2x, xi2y, 
    detA, hMh, cxy, zh, Q, Amux[StpMaxDim], Amuy[StpMaxDim],
    *Sc = P(STP_S),
    *M = P(STP_M),
    *z = P(STP_Z);
  model *phi = cov->sub[STP_PHI],
    *Sf = cov->kappasub[STP_S],
    *xi2 =cov->sub[STP_XI2];

  TALLOC_X1(Sx, dimsq);
  TALLOC_X2(Sy, dimsq);
  TALLOC_X3(A, dimsq);
    
  if (Sf != NULL) {
    FCTN(x, Sf, Sx); // symmetric, pos definite !!
    FCTN(y, Sf, Sy);
    //
    //    if (false) {
    //      int ii;
    //      printf("x=%10g %10g y=%10g %10g\n", x[0], x[1], y[0], y[1]);
    //for (ii=0; ii<4; ii++) printf("%10g ", Sx[ii]);
    //      printf("\n");
    //      for (ii=0; ii<4; ii++) printf("%10g ", Sy[ii]);
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
  double 
    dx = Ext_detPosDef(Sx, dim), 
    dy = Ext_detPosDef(Sy, dim);
  END_TALLOC_X1;
  END_TALLOC_X2;

  
  // det_UpperInv(A, &detA, dim); // here
  double muxAmuy; // changed 5.4.2018, 3.2.0.28
  err = Ext_XCinvYdet(A, dim, true, Amux, Amuy, 1, &muxAmuy, &detA, false,NULL);
  END_TALLOC_X3;
  if (err != NOERROR) ERR0("error occuredin 'GetEu2Dinv'");
 
  Q = cxy * cxy - hMh * hMh + muxAmuy;
  if (Q < 0.0) {
    PRINTF("x=%10g,%10g y=%10g,%10g detA=%10g\n", 
	   x[0], x[1], y[0], y[1], detA);
    PRINTF("cxy=%4f hMh=%10g Amux=%10g Amuy=%10g\ndim=%d h=%10g,%10g hSx=%10g,%10g, Q=%10g\n", 
	   cxy, hMh, Amux[0], Amuy[0], 
	   dim, h[0], h[1], hSx[0], hSx[1], Q);
    BUG;
  }

  Q = SQRT(Q);

  aux_covfct auxcf;
  if ((auxcf = DefList[MODELNR(phi)].aux_cov) != NULL) {
    BUG;
    auxcf(x, y, Q, phi, v);
  } else 
    FCTN(&Q, phi, v);

  
  *v *=  POW(2.0, 0.5 * double(dim)) * POW(dx * dy / (detA * detA), 0.25);

}



int checkstp(model *cov){
  ASSERT_UNREDUCED;
  model 
    *phi = cov->sub[STP_PHI],
    *Sf = cov->kappasub[STP_S],
    *xi2 =cov->sub[STP_XI2];
  int err,
    dim = OWNLOGDIM(0);



  ASSERT_CARTESIAN;
  if (dim > StpMaxDim)
    SERR2("For technical reasons max. dimension for ave is %d. Requested is %d",
	  StpMaxDim, GATTERXDIM(0));

 if (PisNULL(STP_S) && Sf==NULL) {  // Sc
   if ((cov->px[STP_S] = EinheitsMatrix(dim)) == NULL) 
     RETURN_ERR(ERRORMEMORYALLOCATION);
    cov->ncol[STP_S] = cov->nrow[STP_S] = dim;
  }
 if (PisNULL(STP_M)) { // M
   if ((cov->px[STP_M] = EinheitsMatrix(dim)) == NULL)
     RETURN_ERR(ERRORMEMORYALLOCATION);
    cov->ncol[STP_M] = cov->nrow[STP_M] = dim;
  }
 if (PisNULL(STP_Z)) { // z
   PALLOC(STP_Z, dim, 1);
 }
   
  if ((err = CHECK(phi, dim,  1, PosDefType, XONLY, ISOTROPIC,
		     SCALAR, cov->frame // VariogramType changed 20.7.14 wg spectral
		   )) != NOERROR)
    RETURN_ERR(err);
  if (!isNormalMixture(phi->monotone)) RETURN_ERR(ERRORNORMALMIXTURE);

  cov->pref[Average] = 5;

  if (Sf != NULL) {
    if ((err = CHECK(Sf, dim,  dim, ShapeType, XONLY, CARTESIAN_COORD,
		       dim, cov->frame // VariogramType changed 20.7.14 wg spectral
)) != NOERROR) 
      RETURN_ERR(err);
  }
  

  if (xi2 != NULL) {
   if ((err = CHECK(xi2, dim, dim,  ShapeType, XONLY, CARTESIAN_COORD,
		    SCALAR, cov->frame // VariogramType changed 20.7.14 wg spectral
		    )) != NOERROR)
     RETURN_ERR(err);
  }

  // kein setbackward??

  EXTRA_STORAGE;

  cov->mpp.maxheights[0] = RF_NA;
  RETURN_NOERROR;
}

void rangestp(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  int i;
  for (i=0; i<=2; i++) { /* S, M, z */
    range->min[i] = RF_NEGINF;
    range->max[i]= RF_INF;
    range->pmin[i] = - 1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}



int structStp(model *cov, model **newmodel) { 
  model *shape;
  int err;
  
  ASSERT_NEWMODEL_NOT_NULL;
  if ((err = covcpy(newmodel, cov)) != NOERROR) RETURN_ERR(err);
  shape = *newmodel;
  SET_NR(shape, SHAPESTP);
  addModel(shape, STP_GAUSS, GAUSS);
   ERR0("'stp' currently does not work"); // what to do with the next line??
  /* 
     shape->sub[STP_GAUSS]->ts dim = 1;
  */
   RETURN_NOERROR;
}


void logshapestp(double *x, double *u, model *cov, double *v, double *Sign){
  // kann um ca. Faktor 2 beschleunigt werden, wenn
  // Sx , log det U, Hx fuer alle x abgespeichert werden
  // und die Werte dann aus dem Speicher gelesen werden
  // jedoch sehr Speicherintensiv. MEMCOPY braucht man auch nicht
  model 
    *Sf = cov->kappasub[STP_S],
    *xi2 =cov->sub[STP_XI2];
  int j, k, d, 
    //    dim= cov->xdimprev,
    dim= OWNXDIM(0),
    dimsq = dim * dim,
    bytes = sizeof(double) * dimsq;
  double h[StpMaxDim], hSxh, hSx, xi, Mhd, 
    *Sc = P(STP_S),
    *M = P(STP_M),
    *z = P(STP_Z),
    *q = cov->q;
  assert(q != NULL);
  
  TALLOC_X1(Sx, dimsq); 

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
    + 0.25 * LOG(Ext_detPosDef(Sx, dim))                          // Sx ^1/4 
      - q[AVESTP_V] * hSxh             // EXP(-V(U-x) S (U-x)) 
    // + DefList[phi->nr].logmixdens(x, q[AVESTP_LOGV], phi) // g //nicht gatternr
    //    - 0.5 * cov_a->logdens // f 
    ;                // 1 / SQRT(f) 
  
  if (!(exponent < 5.0) && PL >= PL_DETAILS) {
    if (!(exponent < 6.0)) // could be NA, too
     PRINTF("\n%10g log Det U=%10g %10g expon=%10g",
          0.25 * dim * (// M_LN2 +  ??? !!! Rechnung!!! 
			q[AVESTP_LOGV] - 2.0 * M_LN_SQRT_PI) // (2V/pi)^{d/4}
	    , 0.25 * LOG(Ext_detPosDef(Sx, dim))                /// Sx ^1/4 
	    , -q[AVESTP_V]* hSxh             // EXP(-V(U-x) S (U-x)) 
	    // , DefList[phi->nr].logmixdens(x, q[AVESTP_LOGV],  phi)// g
	    //, - 0.5 * cov_a->logdens // f 
	    , exponent);
    else PRINTF("!");
  };
  
  assert(EXP(exponent) < 10000000.0);
  
 
  double cos_value = COS(q[AVERAGE_YPHASE] + q[AVERAGE_YFREQ] * xi);
  *v = exponent + LOG(FABS(cos_value)) ;  // Y 
  *Sign = cos_value > 0.0 ? 1.0 : cos_value < 0.0 ? -1.0 : 0.0;
  END_TALLOC_X1;
}


int check_shapestp(model *cov) {  
  if (cov->sub[AVE_GAUSS] == NULL)
    SERR1("both submodels must be set to '%.50s'", DefList[GAUSS].nick);
  EXTRA_STORAGE;
  return checkstp(cov); // !! not next
}


int init_shapestp(model *cov, gen_storage *s) {
  ASSERT_GAUSS_METHOD(Average);

  model
    *Sf = cov->kappasub[STP_S],
    *gaussmix = cov->sub[STP_GAUSS];
  double sd,
    *q = cov->q;
  int err = NOERROR;

  if (Sf != NULL) {
    double minmax[2];
    assert(DefList[MODELNR(Sf)].minmaxeigenvalue != NULL);
    DefList[MODELNR(Sf)].minmaxeigenvalue(Sf, minmax);
    if (minmax[0] <= 0.0) ERR0("neg eigenvalue in shape function of 'stp'");
    q[AVESTP_MINEIGEN] = minmax[0];
    q[AVESTP_LOGDET] = (double) //cov->xdimprev
      OWNXDIM(0) * LOG(minmax[1]);
  } else {
#define dummyN (5 * StpMaxDim)
    double value[StpMaxDim], ivalue[StpMaxDim], dummy[dummyN], det,
      min = RF_INF;
    int i, Ferr,       
      dim = OWNLOGDIM(0),
      ndummy = dummyN;
 
    F77dgeev("No", "No", &dim, P(STP_S), // SVD
		    &dim, value, ivalue, NULL, &dim, NULL, &dim,
		    dummy, &ndummy, &Ferr FCONE FCONE);
    if (Ferr != 0) RETURN_ERR(ERRORSVD);
    det =  1.0;
    for (i=0; i<dim; i++) {
      double v = FABS(value[i]);
      det *= v;
      if (min > v) min = v;
    }
    q[AVESTP_MINEIGEN] = min;
    q[AVESTP_LOGDET] = LOG(det);
  }


  q[AVESTP_LOGV] = 0.0; 
  q[AVESTP_LOGMIXDENS] = 0.0;
  sd_avestp(cov, s, OWNLOGDIM(0), &sd); // sd->gauss

  assert(VDIM0 == 1);  assert(VDIM0 == VDIM1);

  if (cov->mpp.moments >= 0) {
    cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0; //// ??? notwendig 
    if (cov->mpp.moments >= 1) {
      if ((err = INIT(gaussmix, 2, s)) != NOERROR) RETURN_ERR(err);
      if (cov->mpp.moments >= 2) cov->mpp.mM[2] = 1.0; 
    }
  }
  RETURN_ERR(err);
}

void do_shapestp(model *cov, gen_storage *s) {
  // Simulation of V; see Bernoulli Sec. 4.2
  model 
    *stpGAUSS = cov->sub[STP_GAUSS],
    *phi = cov->sub[STP_PHI];
  double  spec_ret[StpMaxDim], sd,
    *q = cov->q;


  DefList[MODELNR(phi)].drawmix(phi, &(q[AVESTP_V]));
  sd_avestp(cov, s, OWNLOGDIM(0), &sd); // sd->gauss
    
  BUG;

  SPECTRAL(stpGAUSS, s, spec_ret);  // nicht gatternr
  q[AVERAGE_YFREQ] = *spec_ret * SQRT(q[AVESTP_V]);  
  q[AVERAGE_YPHASE] = TWOPI * UNIFORM_RANDOM;


  BUG; /// what to do with the next line?
  // info->logdens = DefList[phi->nr].logmixdens(ZERO, q[AVESTP_LOGV], phi);


  //info->radius = RF_INF;
  // info-sd s.o.
}




/* nsst */
/* Tilmann Gneiting's space time models, part I */
void nsst(double *x, model *cov, double *v) {
  model *subphi = cov->sub[0];
  model *subpsi = cov->sub[1];
  double v1, v2, psi, y;
  
  COV(ZERO(subpsi), subpsi, &v1);
  COV(x + 1, subpsi, &v2);
  psi = SQRT(1.0 + v1 - v2);  // C0 : C(0) oder 0 // Cx : C(x) oder -gamma(x)
  y = x[0] / psi;
  COV(&y, subphi, v);
  *v *= POW(psi, -P0(NSST_DELTA));
}

void Dnsst(double *x, model *cov, double *v) {
  model *subphi = cov->sub[0];
  model *subpsi = cov->sub[1];
  double v1, v2, psi, y;

  COV(ZERO(subpsi), subpsi, &v1);
  COV(x + 1, subpsi, &v2);
  psi = SQRT(1.0 + v1 - v2);  // C0 : C(0) oder 0 // Cx : C(x) oder -gamma(x)
  y = x[0] / psi;
  Abl1(&y, subphi, v);
  *v *= POW(psi, -P0(NSST_DELTA) - 1.0);
}

void TBM2nsst(double *x, model *cov, double *v) {
  model *subphi = cov->sub[0];
  model *subpsi = cov->sub[1];
  double v1, v2, psi, y;

  assert(false);

  COV(ZERO(subpsi), subpsi, &v1);
  COV(x + 1, subpsi, &v2);
  psi = SQRT(1.0 + v1 - v2);  // C0 : C(0) oder 0 // Cx : C(x) oder -gamma(x)
  y = x[0] / psi;
  TBM2CALL(&y, subphi, v);
  *v *= POW(psi, -P0(NSST_DELTA));
}

int checknsst(model *cov) {
  model *subphi = cov->sub[0];
  model *subpsi = cov->sub[1];
  int err;
      
  if (OWNXDIM(0) != 2) SERR("reduced dimension must be 2");

  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  cov->finiterange = falsch;
  if ((err = CHECK(subphi, OWNLOGDIM(0), 1, PosDefType, XONLY, ISOTROPIC, 
		     SCALAR, cov->frame // VariogramType changed 20.7.14 wg spectral
		   )) != NOERROR) 
    RETURN_ERR(err);
  
  if (!isNormalMixture(subphi->monotone)) return(ERRORNORMALMIXTURE);
  setbackward(cov, subphi);
  assert(cov->finiterange != wahr); //

  if ((err = CHECK(subpsi, 1, 1, VariogramType, XONLY, ISOTROPIC, 
		     SCALAR, cov->frame // VariogramType changed 20.7.14 wg spectral
		   )) != NOERROR) 
    RETURN_ERR(err);

  // kein setbackward(cov, subpsi);
  RETURN_NOERROR;
}


void rangensst(model *cov, range_type* range){ 
  range->min[NSST_DELTA] = OWNLOGDIM(0) - 1;
  range->max[NSST_DELTA] = RF_INF;
  range->pmin[NSST_DELTA] = range->min[NSST_DELTA];
  range->pmax[NSST_DELTA] = 10.0;
  range->openmin[NSST_DELTA] = false;
  range->openmax[NSST_DELTA] = true;
}




/* gennsst */
/* generalisation in schlather, bernoulli 2010 */

#define GENNSST_INTERN_A 0

void kappa_gennsst_intern(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  *nc = *nr = i == 0 ? 0 : -1;
}

void gennsst_intern(double *x, model *cov, double *v) {
  model *next = cov->sub[0];
  double  det, z,
    *A = P(GENNSST_INTERN_A);
  int dim = OWNLOGDIM(0);
 
  //  printf("\nx=%10g %10g %10g \ny=%10g %10g %10g\nz=%10g %10g %10g\n %10g %10g %10g\n %10g %10g %10g\n",  	 x[0], x[1], x[2], y[0],y[1],y[2],	 A[0],  A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8]);

  //  PMI(cov->calling->sub[1]);

  //  det_UpperInv(A, &det, dim);
  if (Ext_XCinvXdet(A, dim, x, 1, &z, &det, false, NULL) != NOERROR) {
    *v = RF_NAN;
    return;
  }
  
  z = SQRT(z);
  COV(&z, next, v);
  *v /= SQRT(det);
}

int checkgennsst_intern(model *cov) {
  // internes modell muss zwischengeschaltet werden
  // um Earth->Cartesian->isotropic hinzubekommen
  // und         ^- hier noch mit A zu multiplizieren
  model *next = cov->sub[0];
  int err,
    dim = OWNXDIM(0);

  assert(OWNLOGDIM(0) == OWNXDIM(0));
  if ((err = CHECK(next, OWNLOGDIM(0), 1, PosDefType, 
		   XONLY, ISOTROPIC, SCALAR, cov->frame)) != NOERROR) 
    RETURN_ERR(err);
  if (!isNormalMixture(next->monotone)) RETURN_ERR(ERRORNORMALMIXTURE);

  if (PisNULL(GENNSST_INTERN_A)) {  
    PALLOC(GENNSST_INTERN_A, dim, dim);
  } else if (cov->nrow[GENNSST_INTERN_A] != dim) {
    PFREE(GENNSST_INTERN_A);
    PALLOC(GENNSST_INTERN_A, dim, dim);
  }

  cov->finiterange = falsch;
  setbackward(cov, next);
  assert(cov->finiterange != wahr); 
  VDIM0 = VDIM1 = 1;

  EXTRA_STORAGE;
  RETURN_NOERROR;
}


void range_gennsst_intern(model VARIABLE_IS_NOT_USED *cov,
			  range_type* range){
  range->min[GENNSST_INTERN_A] = RF_NEGINF;
  range->max[GENNSST_INTERN_A] = RF_INF;
  range->pmin[GENNSST_INTERN_A] = -100.0;
  range->pmax[GENNSST_INTERN_A] = +100.0;
  range->openmin[GENNSST_INTERN_A] = true;
  range->openmax[GENNSST_INTERN_A] = true;
}


#define GENNSST_DIM_U 0
void gennsst(double *x,  model *cov, double *v) {
  model *subphi = cov->key;
  model *subpsi = cov->sub[1];
  int 
    totaldim = ANYDIM, 
    udim = P0INT(GENNSST_DIM_U),
    hdim = totaldim - udim,
    hdimSq = hdim * hdim;
    
  COV(x + hdim, subpsi, PARAM(subphi, GENNSST_INTERN_A));
  if (isnowVariogram(subpsi)) {
    TALLOC_X1(z, hdimSq);
    double *zero = ZERO(subpsi),
      *p = PARAM(subphi, GENNSST_INTERN_A);
    COV(zero, subpsi, z);
    for (int i=0; i<hdimSq; i++) p[i] = z[i] - p[i];
    END_TALLOC_X1;
  } else if (!equalsnowNegDef(subpsi)) BUG;

  COV(x, subphi, v);
  if (ISNAN(*v)) ERR0("error occuredin 'GetEu2Dinv'");
}

void nonstatgennsst(double *x,double *y, model *cov, double *v) {
  model *subphi = cov->key;
  model *subpsi = cov->sub[1];
  int 
    totaldim = ANYDIM, 
    udim = P0INT(GENNSST_DIM_U),
    hdim = totaldim - udim,
    hdimSq = hdim * hdim;
    
  NONSTATCOV(x + hdim, y + hdim, subpsi, PARAM(subphi, GENNSST_INTERN_A));
  TALLOC_X1(z, hdimSq);
  if (isnowVariogram(subpsi)) {
    double *zero = ZERO(subpsi),
      *p = PARAM(subphi, GENNSST_INTERN_A);
    NONSTATCOV(zero, zero, subpsi, z);
    for (int i=0; i<hdimSq; i++) p[i] = z[i] - p[i];
  } else if (!equalsnowNegDef(subpsi)) BUG;
  
  for (int d=0; d<hdim; d++) z[d] = x[d] - y[d];
  COV(z, subphi, v);
  END_TALLOC_X1;
  if (ISNAN(*v)) ERR0("error occuredin 'GetEu2Dinv'");
 }
  

#define GENNSST_DIM_U 0
int checkgennsst(model *cov) {
  model *subphi = cov->sub[0];
  model *subpsi = cov->sub[1];
#define ncs 3
  int err;

  kdefault(cov, GENNSST_DIM_U, 1);  
  if (OWNLOGDIM(0)!=OWNXDIM(0)) SERR("logical and physical dimension differ");
  int udim = P0INT(GENNSST_DIM_U),
    hdim = OWNLOGDIM(0) - udim;

  if (cov->q == NULL) {
    QALLOC(1);
    cov->q[0]=NOERROR;
  }
  
  if (isSpherical(OWNISO(0))) 
    return cov->q[0] != NOERROR ? (int) cov->q[0] : ERRORFAILED;
  if (cov->key == NULL) {
    if ((err = covcpy(&(cov->key), subphi)) != NOERROR) RETURN_ERR(err);
    addModel(&(cov->key), GENNSST_INTERN);
  }
  
  if ((cov->q[0] = err = CHECK(cov->key, hdim, hdim, PosDefType, 
			       XONLY, SYMMETRIC, // viel zu schwach, geht aber
			       //                     z.Zt. nicht anders
			       // hinsichtlich der Koordinatentransformation
			       // -- wird das interne modell abgefangen
			       SCALAR, cov->frame)) != NOERROR) {
    RETURN_ERR(err);
  }
  if (!isNormalMixture(cov->key->sub[0])) SERR("'phi' not a normal mixture.");
  if (TOTALXDIM(SYSOF(cov->key)) != hdim)
    SERR("given dim does not match dimension required for 'phi'");
  
  // MUSSS ZWINGEND ALS ZWEITES KOMMEN
  domain_type type = XONLY;
  while (true) {
    if ((err = CHECK(subpsi, udim, udim, NegDefType, type,
		     OWNISO(0), hdim, cov->frame)) == NOERROR) break;
    if (type == XONLY) type = KERNEL; else RETURN_ERR(err);
  }
  
  if ((!equalsSpaceIsotropic(OWNISO(0)) ||
       MODELNR(subpsi) != MATRIX || PisNULL(M_M) ||
       subpsi->kappasub[M_M] != NULL || subpsi->nsub > 1 ||
       subpsi->sub[0]->vdim[0] != 1)
      && !equalsSymmetric(OWNISO(0))) RETURN_ERR(ERRORWRONGISO);
      
  cov->finiterange = falsch;
  setbackward(cov, cov->key);
  assert(cov->finiterange != wahr); //

  VDIM0 = VDIM1 = 1;
  
  
  COV_DELETE(cov->sub + 0, cov);
  if ((err = covcpy(cov->sub + 0, cov->key->sub[0])) != NOERROR) BUG;
  SET_CALLING(cov->sub[0], cov);
 
  // kein setbackward(cov, subpsi);
  EXTRA_STORAGE;
  RETURN_NOERROR;
}



bool allowedDgennsst(model *cov) {
  allowedD(cov->sub[1]);
  MEMCOPY(cov->allowedD, cov->sub[1]->allowedD, sizeof(allowedD_type));
  return false;
}

bool allowedIgennsst(model *cov) {
  bool *I = cov->allowedI;
  model *subpsi = cov->sub[1];
  
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  I[(int) SYMMETRIC] = true;
  I[(int) DOUBLEISOTROPIC] = MODELNR(subpsi) == MATRIX && PisNULL(M_M) &&
    subpsi->kappasub[M_M] == NULL && subpsi->nsub > 1;
  return false;
}

void rangegennsst(model VARIABLE_IS_NOT_USED *cov, range_type* ra){
  ra->min[GENNSST_DIM_U] = 1;
  ra->max[GENNSST_DIM_U] = OWNLOGDIM(0) - 1;
  ra->pmin[GENNSST_DIM_U] = 1;
  ra->pmax[GENNSST_DIM_U] = ra->max[GENNSST_DIM_U];
  ra->openmin[GENNSST_DIM_U] = false;
  ra->openmax[GENNSST_DIM_U] = false;
}

