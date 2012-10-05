/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives of hypermodels

 Copyright (C) 2005 -- 2011 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
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
#include "Covariance.h"
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>


// ------- end internal functions

void co(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double  y=*x, *q=cov->q,
    diameter = cov->p[pLOC_DIAM][0],  
    a = cov->p[pLOC_A][0];
  if (y <= diameter) CovList[next->nr].cov(x, next, v);
  else {
    *v = (y >= q[LOCAL_R]) ? 0.0
      : q[CUTOFF_B] * pow(q[CUTOFF_ASQRTR] - pow(y, a), 2.0 * a);
  }
}

bool alternativeparam_co(cov_model *cov){
  return false;
}
void range_co(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0; //  CUTOFF_DIAM
  range->max[0] = RF_INF;
  range->pmin[0] = 1e-10;
  range->pmax[0] = 1e10;
  range->openmin[0] = true;
  range->openmax[0] = true;

  range->min[1] = 0.0; // cutoff_a
  range->max[1] = RF_INF;
  range->pmin[1] = 0.5;
  range->pmax[1] = 2.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
}

int check_local(cov_model *cov,  
		 SimulationType type, int maxq, 
		 getlocalparam init, // next->coinit
		 set_local_q_type set_local) {

  int i, err=NOERROR, msg; // timespacedim -- needed ?;
  double  *q, q2[LOCAL_MAX],
    d=RF_NAN, **p = cov->p;
  cov_model *next = cov->sub[0],
    *truenext = next->sub[0];
  localinfotype li;
  cov_fct *C;
  cov->qlen = maxq;

 
  //  checkkappas(cov) is not called, but checked directly !!
  if ((err = check2X(next, cov->tsdim, 1, 
		     type == CircEmbedCutoff ? STATIONARY : VARIOGRAM, 
		     ISOTROPIC, UNIVARIATE)) != NOERROR)
    return err;

  next->tsdim = DEL_COV;
//  setbackward(cov, next); 
  if (next->maxdim < cov->maxdim) cov->maxdim = next->maxdim;
  cov->manipulating_x=true;
  cov->normalmix = false;
  // cov->finiterange = true by default
  cov->diag &= next->diag;
  //  cov->quasidiag &= next->quasidiag;
  cov->separatelast &= next->separatelast; 
  cov->semiseparatelast &= next->semiseparatelast; 
  // cov->variogram = false; 
  // no changing in pref by submodel !!
  if (cov->q != NULL) {
      free(cov->q);
      // PrintModelInfo(cov);
      // ERR("q not NULL in check_local -- ask author");
  }
  q = cov->q = (double*) malloc(maxq * sizeof(double));
  for (i = 0; i < maxq; i++) q2[i] = RF_NAN; // q2 will be copied to cov->q;
//  {int i; printf("%d %f\n", cov->qlen); for(i=0;i<9;i++)printf("%f ", cov->q[i]); printf("\n");  assert(false);}
  
  
  C = CovList + truenext->nr; // next is GATTER only

  if (p[pLOC_DIAM] == NULL) {
    ERR("Diameter must always be given");
    // cov->p[pLOC_DIAM] = (double*) malloc(sizeof(double));
    // cov->ncol[0] = cov->nrow[0] = 1;
    // cov->p[pLOC_DIAM][0] = d = GetScaledDiameter(cov);
  } else {
    if (cov->ncol[pLOC_DIAM] != 1 || cov->nrow[pLOC_DIAM] != 1) 
      ERR("diameter must be a scale");
    d = cov->p[pLOC_DIAM][0];
  }

  if (p[pLOC_A] == NULL) {
    if (C->implemented[type]) {
      int i;
      assert(init != NULL);

      // PrintModelInfo(truenext);

      init(truenext, &li);
      if (li.instances == 0) {
	ERR("parameter values do not allow for finding second parameter");
      }
      cov->p[pLOC_A] = (double*) malloc(sizeof(double));
      cov->ncol[pLOC_A] = cov->nrow[pLOC_A] = 1;
      q[LOCAL_R] = R_PosInf;

      // printf("%f\n", d); assert(false);

     msg = MSGLOCAL_FAILED;
      for (i=0; i<li.instances; i++) {
	if (li.msg[i] <= msg) {
	  err = set_local(truenext, li.value[i], d, q2);
	  if (err == NOERROR && (li.msg[i] < msg || q2[LOCAL_R] < q[LOCAL_R])){
	    memcpy(q, q2, sizeof(double) * maxq);
	    msg = li.msg[i];
	    cov->p[pLOC_A][0] = li.value[i]; // co: a; ie: rawr
	  }
	}
      }
      cov->q[LOCAL_MSG] = msg;
      if (msg == MSGLOCAL_FAILED) err = ERRORFAILED;
    } else {
	// PRINTF("%s\n", C->name);
      ERR("2nd parameter is neither given nor can be found automatically");
    }
  } else {
    if (cov->ncol[pLOC_A] != 1 || cov->nrow[pLOC_A] != 1) 
      ERR("`a' must be a scale");
    err = set_local(truenext, cov->p[pLOC_A][0], d, q2);
    memcpy(q, q2, sizeof(double) * maxq);
  }

//  printf("pdq %f %f %f\n", cov->p[pLOC_A][0], d, q2);

  if (err != NOERROR) { XERR(err)
    //   char Msg[255];
    //   errorMessage(err);
    // error(Msg);
  }
  return NOERROR;
}

int set_cutoff_q(cov_model *next, double a, double d, double *q) {
  double phi0, phi1, a2; //, one=1.0;
  a2 = a * a;
  cov_fct *C = CovList + next->nr;
  
//  PrintModelInfo(next); 

  C->cov(&d, next, &phi0);
  C->D(&d, next, &phi1);
  phi1 *= d;
  if (phi0 <= 0.0) return MSGLOCAL_SIGNPHI;
  if (phi1 >= 0.0) return MSGLOCAL_SIGNPHIFST;

//  printf("check_co phi0=%f phi1=%f a=%f d=%f a2=%f\n", phi0, phi1, a, d, a2);
//  assert(false);
//  printf("p %f\n", phi1);
//  printf("a2 %f\n", a2);
//  printf("p0 %f\n", phi0);
//  printf("a %f\n", a);
//  printf("d %f\n", d);
  q[CUTOFF_B] =//sound qconstant even if variance of submodel is not 1
    pow(- phi1 / (2.0 * a2 * phi0), 2.0 * a) * phi0 / pow(d, 2.0 * a2);
//  assert(false);
  q[CUTOFF_THEOR] = pow(1.0 - 2.0 * a2 * phi0 / phi1, 1.0 / a);
  q[LOCAL_R] = d * q[CUTOFF_THEOR];
  q[CUTOFF_ASQRTR] = pow(q[LOCAL_R], a);

//  printf("phi0=%f %f %f %f\n", phi0, phi1, a, d);
//   printf("%f\n", d); assert(false);
//  printf("%f %f %f %f\n", q[CUTOFF_B], q[CUTOFF_THEOR], q[LOCAL_R], q[CUTOFF_ASQRTR]); //assert(false);

  return NOERROR;
}
int check_co(cov_model *cov) {
  cov_model *next = cov->sub[0]->sub[0];// sub[0] is GATTER ONLY
  cov_fct *C = CovList + next->nr;  
  return check_local(cov, CircEmbedCutoff, CUTOFF_MAX, C->coinit,
		     set_cutoff_q);
}


void Stein(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;
  double y=*x, z, *q=cov->q,
    diameter = cov->p[pLOC_DIAM][0];
  if (y <= diameter) {
    C->cov(x, next, v);
    *v += q[INTRINSIC_A0] + q[INTRINSIC_A2] * y * y;

    //   printf("%f %f %f %f %d %d \n ", 
//	   y, v[0], q[INTRINSIC_A0], q[INTRINSIC_A2], INTRINSIC_A0, INTRINSIC_A2);
//    PrintModelInfo(cov);
//    assert(v[0] < RF_INF); 
  } else {
     z = q[LOCAL_R] - y;
     *v = (z <= 0.0) ? 0.0 : q[INTRINSIC_B] * z * z * z / y;
  }
}

bool alternativeparam_Stein(cov_model *cov) {
  cov->p[pLOC_A][0] *= 2.0;
  return true;
}

void range_Stein(cov_model *cov, range_arraytype* ra) 
{
  range_type *range = ra->ranges;
  range->min[0] = 0.0; //  CUTOFF_DIAM
  range->max[0] = RF_INF;
  range->pmin[0] = 0.01;
  range->pmax[0] = 100;
  range->openmin[0] = true;
  range->openmax[0] = true;

  range->min[1] = 1.0; // stein_r
  range->max[1] = RF_INF;
  range->pmin[1] = 1.0;
  range->pmax[1] = 20.0;
  range->openmin[1] = false;
  range->openmax[1] = true;
}

int set_stein_q(cov_model *next, double r, double d, double *q) {
  double C0, phi0, phi1, phi2, 
    zero = 0.0, 
    rP1 = r + 1.0,
    rM1 = r - 1.0, 
    dsq = d * d;
  cov_fct *C = CovList + next->nr;
//  printf("%ld %s\n", &zero, CovList[next->nr].name);
  // PrintModelInfo(next);
  C->cov(&zero, next, &C0);

  C->cov(&d, next, &phi0);
  C->D(&d, next, &phi1); // derivative of
  phi1 *= d;             // scaled fctn at 1
  C->D2(&d, next, &phi2); // 2nd derivative of
  phi2 *= dsq;            // scaled fctn at 1
  q[LOCAL_R] =  r * d;   
  q[INTRINSIC_A2] = (phi2 - phi1) / (3.0 * r * rP1) ;
  q[INTRINSIC_B]  = (r == 1.0) ? 0.0 : q[INTRINSIC_A2] / (rM1 * dsq);
  q[INTRINSIC_A2] = (q[INTRINSIC_A2] - phi1 / 3.0 - phi2 / 6.0) / dsq;
  q[INTRINSIC_A0] = 0.5 * rM1 / rP1 * phi2 + phi1 / rP1 - phi0;


//  printf("check: %f %e r=%f intr:%f %f  %f\n", 
//	  phi2, phi1, r, q[INTRINSIC_B],
//	  q[INTRINSIC_A0],  q[INTRINSIC_A2]);


//   assert(false);

  if ((q[INTRINSIC_B]  < 0.0) || (q[INTRINSIC_A2] < 0.0) ||
      (q[INTRINSIC_A0] + C0 < 0.0)) 
    return MSGLOCAL_INITINTRINSIC;
  return NOERROR;
}

int check_Stein(cov_model *cov)
{  
  cov_model *next = cov->sub[0]->sub[0]; // dito
  cov_fct *C = CovList + next->nr;

//  printf("******************* check stein\n");
//  PrintModelInfo(cov);
//  assert(false);

  return check_local(cov, CircEmbedIntrinsic, INTRINSIC_MAX, C->ieinit,
		     set_stein_q);
}     


void MaStein(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  double  nuG, gammas, **p = cov->p, v1, v2,
    nu = p[0][0], delta = p[1][0];
  cov_fct *C = CovList + next->nr;

  C->cov(ZERO, next, &v1);
  C->cov(x + 1, next, &v2);
  nuG = nu + v1 + v2;
  if (nuG >= 80.0) {
      error("Whittle Matern function cannot be evaluated with parameter value b+g(t) greater than 80.");
  }
  gammas = lgammafn(nu + delta) - lgammafn(nu) -lgammafn(nuG + delta);
  double s = *x;
  *v = (s==0.0) ? exp(lgammafn(nuG) + gammas)
    : 2.0 * exp(nuG * log(0.5 * s) + gammas + log(bessel_k(s, nuG, 2.0)) - s);
}
int check_MaStein(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;

  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->normalmix = cov->finiterange = false;
  cov->manipulating_x=true;
  if ((err = check2X(next, cov->tsdim, cov->xdim,
		     STATIONARY, ANISOTROPIC,
		     UNIVARIATE)) != NOERROR) return err;
  if (cov->ncol[0] != 1 || cov->nrow[0] != 1) ERR("nu not scalar");
  if (cov->ncol[1] != 1 || cov->nrow[1] != 1) ERR("d not scalar");
//  if (cov->ncol[2] != 1 || cov->nrow[2] != cov->xdim - 1) {
//    sprintf(ERRORSTRING, "%s",
//       "dimension of vector V must be 1 less than the dim of the coordinates");
//    return ERRORM;
  // }
  return NOERROR;
  // no setbackward !
}
void range_MaStein(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0; // 
  range->max[0] = RF_INF;
  range->pmin[0] = 1e-2;
  range->pmax[0] = 10.0;
  range->openmin[0] = true;
  range->openmax[0] = true;

  range->min[1] = 0.5 * (double) (cov->tsdim - 1); // d
  range->max[1] = RF_INF;
  range->pmin[1] = range->min[1];
  range->pmax[1] = 10;
  range->openmin[1] = false;
  range->openmax[1] = true;
}



/* bivariate shift model */
void kappashift(int i, cov_model *cov, int *nr, int *nc){
    *nc = 1;
    *nr = (i<=0) ? cov->tsdim : -1;
}
void shift(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;
  double y[ShiftMaxDim],
    *h = cov->p[0];
  int i,
    dim = cov->tsdim;
  
  C->cov(x, next, v + 0);
  for (i=0; i<dim; i++) y[i] = x[i] + h[i];
  C->cov(y, next, v + 1);
  for (i=0; i<dim; i++) y[i] = x[i] - h[i];
  C->cov(y, next, v + 2);
  v[3] = v[0];
}


int checkshift(cov_model *cov) {
  cov_model *next = cov->sub[0];
  // cov_fct *C = CovList + next->nr;
  int err,
    dim = cov->tsdim;
      
  if (cov->xdim > ShiftMaxDim) {
      sprintf(ERRORSTRING, 
	      "For technical reasons max. dimension for ave is %d. Got %d.", 
	      StpMaxDim, cov->xdim);
    return ERRORMSG;
  }
  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->manipulating_x=true;
  cov->normalmix = cov->finiterange = false;
  if ((err = check2X(next, dim, dim, STATIONARY, 
		     (dim > 1) ? ANISOTROPIC : ISOTROPIC, 
		     UNIVARIATE)) != NOERROR) 
    return err;
  setbackward(cov, next);
  return NOERROR;
}

void rangeshift(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = -RF_INF;
  range->max[0] = RF_INF;
  range->pmin[0] = -1000;
  range->pmax[0] = 1000;
  range->openmin[0] = true;
  range->openmax[0] = true;
}



/* Vector - Delta Delta^T */
void vector(double *x, cov_model *cov, double *v) {
  /*
    r = \| x\|
    d/d x_i f(r) = f'(r) x_i / r
    d^2/d x_i^2 f(r) = f''(r) x_i^2/r^2 + f'(r) (r^2 - x_i^2) / r^3
    d^2/d x_i x_j f(r) = f''(r) x_i x_j/r^2 + f'(r) (- x_i x_j) / r^3

    r = 0:
    d/d x_i f(r) |_0 = 0
    d^2/d x_i^2 f(r) |_0 = f''(0)
    d^2/d x_i x_j f(r) = 0

    Vector (n = dimension of spatial field)
    r=0
    L = n * f''(0)

    operator:  -0.5 * (a + 1) laplace + a * hessian,  a in [-1,1]; a=1
  */
  cov_model *next = cov->sub[0];
  cov_fct *N = CovList + next->nr;
  double norm[2], normSq0, normL2, normT2, D, D2, diag,
    a = cov->p[0][0],
    b = - 0.5 * (1.0 + a);
  int i,
    td = cov->tsdim,
    dim = ((int *)cov->p[1])[0],
    dimP1 = dim + 1,
    dimsq = dim * dim;
  
  normL2 = normT2 = 0.0;
  for (i=0; i<dim; i++) normL2 += x[i] * x[i];
  for (; i<td; i++) normT2 += x[i] * x[i];
  if (next->isoIn == ISOTROPIC) {
    normSq0 = normL2 + normT2;
  } else {
    normSq0 = normL2;
    norm[1] = sqrt(normT2); 
  }
  norm[0] = sqrt(normSq0);
  N->D(norm, next, &D);
  N->D2(norm, next, &D2);

  if (normSq0 == 0.0) {
    diag = (b * dim + a) * D2;
    for (i=0; i<dimsq; i++) {
      //    printf("%f D2 %f\n", diag, D2);
      v[i] = (i % dimP1 == 0) ? diag : 0.0; 
    }
  } else {
    double 
      D2nsq = D2 / normSq0,
      D1n = D / norm[0],
      D1n3 = D / (norm[0] * normSq0),
      delta = a * (D2nsq - D1n3),
      Vector = b * ((D2nsq - D1n3) * normL2 + dim * D1n);
    int k,l;
    
    diag = Vector + a * D1n;    
    for (i=k=0; k<dim; k++) {
      for (l=0; l<dim; l++, i++) {
	v[i] = (i % dimP1 == 0) ? diag : 0.0; 
	v[i] += delta * x[k] * x[l];
      }
    }
  }
}


/* Vector - Delta Delta^T */
void vectorAniso(double *x, cov_model *cov, double *v) {
  /*
    operator:  -0.5 * (a + 1) laplace + a * hessian,  a in [-1,1]; a=1
  */

  cov_model *next = cov->sub[0];
  cov_fct *N = CovList + next->nr;
  double *D,  laplace,
    a = cov->p[0][0],
    b = - 0.5 * (1.0 + a);
  int i, j, k, endfor,
      dim = ((int *)cov->p[1])[0],    
      dimP1 = dim + 1,
      dimsq = dim * dim,
      xdim = cov->tsdim,
      xdimsq = xdim * xdim,
      xdimP1 = xdim + 1,
      dimxdim = dim * xdim;
  
  D = (double*) malloc(sizeof(double) * xdimsq); 

  //printf("%d %d\n", dim, xdim); assert(false);

  // should pass back the hessian
  N->hess(x, next, D);
  laplace = 0.0;

  // PrintModelInfo(next);

//  printf("return matrix %d\n", xdimsq);
//  for (i=0; i<xdimsq; i+=xdim) {
//      for (j=i; j <i+xdim; j++) {
//	  printf("%f;%d  ", D[j], j);
//      }
//      printf("\n");
//  }


  for (i=0; i<dimxdim; i+=xdimP1) {
      laplace += D[i];
//      printf("lplace %d %f\n", i, D[i]);
  }
  laplace *= b;

  k = 0;
  for (i=0; i<dimxdim; i+=xdim) {  
      endfor = i + dim; 
      for (j=i; j<endfor; j++) {
//	  printf(">> %d %d %f %f %f\n", k, j, D[j], a, laplace );
	  v[k++] = D[j] * a;
      }
  }
  
  for (i=0; i<dimsq; i+=dimP1) {
     v[i] += laplace;
  }

  free(D);
}

int checkvector(cov_model *cov) {
    cov_model *next = cov->sub[0],
	*genuinenext = next->sub[0];
//  cov_fct *N = CovList + next->sub[0]->nr;
  int err,
    dim = cov->tsdim;
  isotropy_type isotropy;

  cov->nr = VECTOR; // wegen nr++ unten;
  cov->normalmix = cov->finiterange = false;
  cov->manipulating_x=true;
  if ((err = check2X(next, dim, 1, STATIONARY, ISOTROPIC,
		     UNIVARIATE)) != NOERROR) {
//     printf("err=%d\n", err);
    if ((err = check2X(next, dim, dim, STATIONARY, ZEROSPACEISO,
		       UNIVARIATE)) != NOERROR) {
//	printf("err2=%d\n", err);
     if ((err = check2X(next, dim, dim, STATIONARY, ANISOTROPIC,
			 UNIVARIATE)) != NOERROR) {
//	 printf("err3=%d\n", err);
	  return err;
      }
    }
  }
  assert(cov->xdim == cov->tsdim);
  if (genuinenext->derivatives < 2 && !genuinenext->hess) {
//      PrintModelInfo(cov);
      //     printf("%s (%d): derivatives %d \n", CovList[genuinenext->nr].name, 
//	     genuinenext->nr, genuinenext->derivatives);
      // next is GATTER
      ERR("2nd derivative of submodel not defined (for the given paramters)");
  }
  isotropy = next->isoIn;

  if (isotropy != ISOTROPIC && isotropy != SPACEISOTROPIC) {
      if (!genuinenext->hess) ERR("hess matrix not defined");
      cov->nr++;
  }
  kdefault(cov, 0, 0.5); // a
  kdefault(cov, 1, 
	   (isotropy==SPACEISOTROPIC  || isotropy==ZEROSPACEISO) 
	   ? dim - 1 : dim); // Dspace
  if ( (isotropy==SPACEISOTROPIC || isotropy==ZEROSPACEISO) 
       && ((int*) (cov->p[1]))[0] != dim - 1) {
//      printf(".... %d %d\n", ((int*) (cov->p[1]))[0], dim);
    ERR("for spatiotemporal models 'vector' must be applied to spatial part");
  }
  
//  PrintModelInfo(cov);

  setbackward(cov, next);
  cov->vdim = ((int*) (cov->p[1]))[0];
  next->tsdim = DEL_COV;
  return NOERROR;
}

void rangevector(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = -1.0;
  range->max[0] = 1.0;
  range->pmin[0] = -1.0;
  range->pmax[0] = 1.0;
  range->openmin[0] = false;
  range->openmax[0] = false;

  range->min[1] = 1.0;
  range->max[1] = cov->tsdim;
  range->pmin[1] = 1.0;
  range->pmax[1] = cov->tsdim;
  range->openmin[1] = false;
  range->openmax[1] = false;
}



int zzz=0;

/* Rot - Delta Delta^T */
void curl(double *x, cov_model *cov, double *v) {
  /*
    r = \| x\|
    d/d x_i f(r) = f'(r) x_i / r
    d^2/d x_i^2 f(r) = f''(r) x_i^2/r^2 + f'(r) (r^2 - x_i^2) / r^3
    d^2/d x_i x_j f(r) = f''(r) x_i x_j/r^2 + f'(r) (- x_i x_j) / r^3

    r = 0:
    d/d x_i f(r) |_0 = 0
    d^2/d x_i^2 f(r) |_0 = f''(0)
    d^2/d x_i x_j f(r) = 0

    Rot (n = dimension of spatial field)
    r=0
    L = n * f''(0)

    operator:  -0.5 * (a + 1) rot + a * hessian,  a in [-1,1]; a=1/2
  */
 
  cov_model *next = cov->sub[0];
  cov_fct *N = CovList + next->nr;
  double norm[2], normSq0, normL2, normT2, D, D2, D3, diag,
      a = -1.0, // curl free
      b = - 0.5 * (1.0 + a);
  int i,
      td = cov->tsdim,
      dim = td, // ((int *)cov->p[1])[0],
      // currently no space-time allowed -> BA thesis
      dimP1 = dim + 1,
      dimP2 = dim + 2,
      dimP3 = dim + 3,
      dimP2sqM1 =  dimP2 * dimP2 -1; // for three dimensions much
  // more complicated
  
  normL2 = normT2 = 0.0;
  for (i=0; i<dim; i++) normL2 += x[i] * x[i];
  for (; i<td; i++) normT2 += x[i] * x[i];
  if (next->isoIn == ISOTROPIC) {
    normSq0 = normL2 + normT2;
  } else {
    normSq0 = normL2;
    norm[1] = sqrt(normT2); 
  }
  norm[0] = sqrt(normSq0);
  N->D(norm, next, &D);
  N->D2(norm, next, &D2);
  N->D3(norm, next, &D3);


  if (normSq0 == 0.0) {
    for (i=0; i<=dimP2sqM1; i++) v[i] = 0.0;

    N->cov(norm, next, v); // v[0]
    
    diag = (b * dim + a) * D2;
    for (i=dimP3; i<dimP2sqM1; i+=dimP3) v[i] = diag; // diag 2,3

    N->D2(norm, next, v + dimP1);
    v[dimP1] *= 2.0;
    v[dimP2 * dimP1] = v[dimP1];
    N->D4(norm, next, v + dimP2sqM1);
    v[dimP2sqM1] *= (8.0 / 3.0); 
  } else {
    double 
      z[2],
      D2nsq = D2 / normSq0,
      D1n = D / norm[0],
      D3n = D3 / norm[0],
      D1n3 = D / (norm[0] * normSq0),
      delta = a * (D2nsq - D1n3),
      div = b * ((D2nsq - D1n3) * normL2 + dim * D1n);
    int k,l;
    
    N->cov(norm, next, v); // v[0]
    
    z[0] = x[0];
    z[1] = x[1];
    for (i=1; i<=dim; i++) { // kante links und oben
	v[i] = -(v[i * dimP2] = z[i-1] * D1n);
    }
 
    diag = div + a * D1n;    
    for (i= dimP3, k=0; k<dim; k++, i+=2) {
      for (l=0; l<dim; l++, i++) {
	v[i] = (i % dimP3 == 0) ? diag : 0.0; 
	v[i] += delta * x[k] * x[l];
      }
    }

    v[dimP1 * dimP2] = v[dimP1] = - v[dimP3] - v[2 * dimP3]; // ecken

    for (i=1; i<=dim; i++) { // kanten unten und hinten
	v[dimP2 * dimP1 + i] =
	    - (v[dimP2 * (i+1) - 1] = z[i-1] * (D3n + D2nsq - D1n3));
    }

    N->D4(norm, next, v + dimP2sqM1); // rechts unten 
    v[dimP2sqM1] += 2.0 * D3n - D2nsq + D1n3;
 
  } 

// for (i=0; i<=15; i++) {
//      if (i!=15) v[i] = 0.0;
//  }

}


/* Div - Delta Delta^T */
void div(double *x, cov_model *cov, double *v) { // div -free !!
  /*
    r = \| x\|
    d/d x_i f(r) = f'(r) x_i / r
    d^2/d x_i^2 f(r) = f''(r) x_i^2/r^2 + f'(r) (r^2 - x_i^2) / r^3
    d^2/d x_i x_j f(r) = f''(r) x_i x_j/r^2 + f'(r) (- x_i x_j) / r^3

    r = 0:
    d/d x_i f(r) |_0 = 0
    d^2/d x_i^2 f(r) |_0 = f''(0)
    d^2/d x_i x_j f(r) = 0

    Div (n = dimension of spatial field)
    r=0
    L = n * f''(0)

    operator:  -0.5 * (a + 1) div + a * hessian,  a in [-1,1]; a=1
  */
  cov_model *next = cov->sub[0];
  cov_fct *N = CovList + next->nr;
  double norm[2], normSq0, normL2, normT2, D, D2, D3, diag,
      a = 1.0, // divergenz free
      b = - 0.5 * (1.0 + a);
  int i,
      td = cov->tsdim,
      dim = td, // ((int *)cov->p[1])[0],
      // currently no space-time allowed -> BA thesis
      dimP1 = dim + 1,
      dimP2 = dim + 2,
      dimP3 = dim + 3,
      dimP2sqM1 =  dimP2 * dimP2 -1; // for three dimensions much
  // more complicated
  
  normL2 = normT2 = 0.0;
  for (i=0; i<dim; i++) normL2 += x[i] * x[i];
  for (; i<td; i++) normT2 += x[i] * x[i];
  if (next->isoIn == ISOTROPIC) {
    normSq0 = normL2 + normT2;
  } else {
    normSq0 = normL2;
    norm[1] = sqrt(normT2); 
  }
  norm[0] = sqrt(normSq0);
  N->D(norm, next, &D);
  N->D2(norm, next, &D2);
  N->D3(norm, next, &D3);


  if (normSq0 == 0.0) {
    for (i=0; i<=dimP2sqM1; i++) v[i] = 0.0;

    N->cov(norm, next, v); // v[0]
    
    diag = (b * dim + a) * D2;
    for (i=dimP3; i<dimP2sqM1; i+=dimP3) v[i] = diag;  // diag 2,3

    N->D2(norm, next, v + dimP1);
    v[dimP1] *= 2.0;
    v[dimP2 * dimP1] = v[dimP1];
    N->D4(norm, next, v + dimP2sqM1);
    v[dimP2sqM1] *= (8.0 / 3.0);    
  } else {
    double 
      z[2],
      D2nsq = D2 / normSq0,
      D1n = D / norm[0],
      D3n = D3 / norm[0],
      D1n3 = D / (norm[0] * normSq0),
      delta = a * (D2nsq - D1n3),
      div = b * ((D2nsq - D1n3) * normL2 + dim * D1n);
    int k,l;
    
    N->cov(norm, next, v); // v[0]
    
    z[0] = -x[1];
    z[1] = x[0];
    for (i=1; i<=dim; i++) { // kante links und oben
	v[i] = -(v[i * dimP2] = z[i-1] * D1n);
    }
     
    diag = div + a * D1n;    
    for (i= dimP3, k=0; k<dim; k++, i+=2) {
      for (l=0; l<dim; l++, i++) {
	v[i] = (i % dimP3 == 0) ? diag : 0.0; 
	v[i] += delta * x[k] * x[l];
      }
    }

    v[dimP1 * dimP2] = v[dimP1] = - v[dimP3] - v[2 * dimP3]; // ecken

    for (i=1; i<=dim; i++) { // kanten unten und hinten
	v[dimP2 * dimP1 + i] =
	    -(v[dimP2 * (i+1) - 1] =  z[i-1] * (D3n + D2nsq - D1n3));
    }

    N->D4(norm, next, v + dimP2sqM1); // rechts unten 
    v[dimP2sqM1] += 2.0 * D3n - D2nsq + D1n3;
 
  }

//  for (i=0; i<=15; i++) {
//      if (i!=15) v[i] = 0.0;
//  }

}

int checkdivcurl(cov_model *cov) {
  cov_model  *next = cov->sub[0],
      *nexttrue=next->sub[0]; // next is gatter and derivatives have not
  // been programmed yet 
//  cov_fct *N = CovList + next->sub[0]->nr;
  int err,
      dim = cov->tsdim,
      spacedim;
  isotropy_type isotropy;

  cov->normalmix = cov->finiterange = false;
  cov->manipulating_x=true;
  if ((err = check2X(next, dim, 1, STATIONARY, ISOTROPIC,
		     UNIVARIATE)) != NOERROR) {
    if ((err = check2X(next, dim, 1, STATIONARY, SPACEISOTROPIC,
		       UNIVARIATE)) != NOERROR)
    return err;
  }

  
  if (nexttrue->derivatives < 4) ERR("4th derivative of submodel not defined");
  if (cov->tsdim != 2) ERR("currently coded only for dim=2");
  isotropy = next->isoIn;
  spacedim = dim - (isotropy == SPACEISOTROPIC);

  if (isotropy != ISOTROPIC && isotropy != SPACEISOTROPIC)
      ERR("submodel must be spaceisotropic");
  if (spacedim != 2)
       ERR("model currently only defined for the plane");
  
  setbackward(cov, next);
  cov->vdim = spacedim + 2; // only correct for spacedim == 2!!!
  assert(spacedim == 2);
  next->tsdim = DEL_COV;
  return NOERROR;
}

void rangedivcurl(cov_model *cov, range_arraytype* ra){ }




void lp(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];  
  double z,
    p = cov->p[0][0];
  int d,
    dim = cov->tsdim;
  for (z=0.0, d=0; d<dim; d++) {
    z += pow(fabs(x[d]), p);
  }
  z = pow(z, 1 / p);
  
  CovList[next->nr].cov(&z, next, v);
}
int checklp(cov_model *cov) {
  bool skipchecks = GLOBAL.general.skipchecks;
  cov_fct *C = CovList + cov->nr;
  cov_model *next = cov->sub[0];  
  double p;
  int err;

  kdefault(cov, 0, 1.0);
  p =cov->p[0][0];
  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->manipulating_x=true;
      
  strcpy(ERROR_LOC, C->name);
  
  if ((err = check2X(next, cov->tsdim + 1, 1,  STATIONARY, ISOTROPIC,
		     UNIVARIATE)) 
	!= NOERROR) return err;
  next->tsdim = DEL_COV;

  if (p==1) {
    cov->maxdim = next->maxdim - 1;  // true for p==1
    if (cov->maxdim == 0) cov->maxdim = 1;
  } else if (p==2) {
    cov->maxdim = next->maxdim;
  } else if (!skipchecks) ERR("p must be 1 or 2.");
  cov->normalmix = next->normalmix;
  cov->semiseparatelast = false;
  updatepref(cov, next);

  return NOERROR;
}
void rangelp(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = 2.0;
  range->pmin[0] = 0.01;
  range->pmax[0] = 2.0;
  range->openmin[0] = true;
  range->openmax[0] = false;
}



void ma1(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];  
  double z,
      alpha = cov->p[0][0],
    theta = cov->p[1][0];
  CovList[next->nr].cov(x, next, &z);
  *v = pow(theta / (1 - (1-theta) * z), alpha);  
}
int checkma1(cov_model *cov) {
//  bool skipchecks = GLOBAL.general.skipchecks;
  cov_fct *C = CovList + cov->nr;
  cov_model *next = cov->sub[0];  
//  double p;
  int err;

  kdefault(cov, 0, 1.0);
  kdefault(cov, 1, 0.5);
  if ((err = checkkappas(cov)) != NOERROR) return err;
      
  strcpy(ERROR_LOC, C->name);
  
  if ((err = check2X(next, cov->tsdim,  cov->xdim, cov->statIn,
		     cov->isoIn,  UNIVARIATE)) 
	!= NOERROR) return err;

  cov->semiseparatelast = false;
//  updatepref(cov, next);
  return NOERROR;
}
void rangema1(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->min[0] = 0.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 0.01;
  range->pmax[0] = 10;
  range->openmin[0] = true;
  range->openmax[0] = true;

  range->min[1] = 0.0;
  range->max[1] = 1.0;
  range->pmin[1] = 0.0;
  range->pmax[1] = 1.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
}


void ma2(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];  
  double z,z0;
  CovList[next->nr].cov(ZERO, next, &z0);
  CovList[next->nr].cov(x, next, &z);
  z = z0 - z;
  *v = (z == 0) ? 1.0 : (1.0 - exp(-z)) / z;  
}
int checkma2(cov_model *cov) {
//  bool skipchecks = GLOBAL.general.skipchecks;
  cov_fct *C = CovList + cov->nr;
  cov_model *next = cov->sub[0];  
//  double p;
  int err;

  if ((err = checkkappas(cov)) != NOERROR) return err;
  cov->manipulating_x = true;
    
  strcpy(ERROR_LOC, C->name);
  
  if ((err = check2X(next, cov->tsdim,  cov->xdim, VARIOGRAM,
		     cov->isoIn,  UNIVARIATE)) 
	!= NOERROR) return err;

  cov->semiseparatelast = false;
//  updatepref(cov, next);
  return NOERROR;
}
void rangema2(cov_model *cov, range_arraytype* ra){ 
}





void kappaM(int i, cov_model *cov, int *nr, int *nc){
  *nc = 0;
  *nr = i < 1 ? 0 : -1;
}

void M(cov_model *cov, double *Z, double *V) {
  // M^t C M
  cov_model *next = cov->sub[0];
  int 
    *ncol = cov->ncol,
    *nrow = cov->nrow;
  double *MtZ, 
    alpha = 1.0,
    beta = 0.0,
    *M = cov->p[0];
//      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
// C := alpha*op( A )*op( B ) + beta*C,
// M : rows of op(A)
// N : cols of op(B)
// K : cols of op(A)= rows of op(B)
// LD-X : first dimension of X

  if (next->vdim == 1) {
    MtZ = (double*) malloc(sizeof(double) * *ncol * *ncol);
    F77_CALL(dgemm)("T", "N", ncol, ncol, nrow, Z, M, nrow, M, nrow, 
		    &beta, V, ncol); 
  } else {
    MtZ = (double*) malloc(sizeof(double) * *nrow * *ncol);
    F77_CALL(dgemm)("T", "N", ncol, nrow, nrow, &alpha, M, nrow, Z, nrow, 
		    &beta, MtZ, ncol); 
    F77_CALL(dgemm)("N", "N", ncol, ncol, nrow, &alpha, MtZ, ncol, M, nrow, 
		    &beta, V, ncol);
  }
  free(MtZ);
}
void Mstat(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int    ncol = cov->ncol[0];
  double *z = (double*) malloc(sizeof(double) * ncol * ncol);
  CovList[next->nr].cov(x, next, z);  
  M(cov, z, v);
  free(z);
}
void Mnonstat(double *x, double *y, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int    ncol = cov->ncol[0];
  double *z = (double*) malloc(sizeof(double) * ncol * ncol);
  CovList[next->nr].nonstat_cov(x, y, next, z);  
  M(cov, z, v);
  free(z);
}
int checkM(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;

  if ((err = check2X(next, cov->tsdim, cov->xdim, cov->statIn, 
		     cov->isoIn, cov->nrow[0])) != NOERROR) return err;
  setbackward(cov, next);
  cov->vdim = cov->ncol[0];
  err = checkkappas(cov);
  return err;
}
void rangeM(cov_model *cov, range_arraytype* ra){ 
  range_type *range = ra->ranges;
  
  range->min[0] = RF_NEGINF;
  range->max[0] = RF_INF;
  range->pmin[0] = -1e10;
  range->pmax[0] = 1e10;
  range->openmin[0] = true;
  range->openmax[0] = true;
  
  range->finiterange = true;
}


void IdStat(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  CovList[next->nr].cov(x, next, v);
}
void IdNonStat(double *x, double *y, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  CovList[next->nr].nonstat_cov(x, y, next, v);
}

int checkId(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;

  cov->vdim = (cov->nrow[0] > 0 && cov->nrow[0] > 0) ?
      ((int*) cov->p[0])[0] : SUBMODELM;
  if ((err = check2X(next, cov->tsdim, cov->xdim, cov->statIn,
		     cov->isoIn, cov->vdim    
		     )) != NOERROR) return err;
  if (cov->vdim == SUBMODELM) cov->vdim = next->vdim;
  setbackward(cov, next);
   
  return NOERROR;
}

void DId(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  CovList[next->nr].D(x, next, v);
}

void DDId(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  CovList[next->nr].D2(x, next, v);
}
void TBM2Id(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  CovList[next->nr].tbm2(x, next, v);
}
int initspectralId(cov_model *cov) {
  cov_model *next = cov->sub[0];
  return CovList[next->nr].initspectral(next);
}
void spectralId(cov_model *cov, spectral_storage *s, double *e) { 
  cov_model *next = cov->sub[0];
  return CovList[next->nr].spectral(next, s, e);
}
void coinitId(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0];
  CovList[next->nr].coinit(next, li);
}
void ieinitId(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0];
  CovList[next->nr].ieinit(next, li);
}
void rangeId(cov_model *cov, range_arraytype* ra){ 
  range_type *range = ra->ranges;
  range->min[0] = 1.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 1.0;
  range->pmax[0] = 10;
  range->openmin[0] = false;
  range->openmax[0] = true;

  range->finiterange = true;
  range->maxdim = INFDIM;
}


void Exp(double *x, cov_model *cov, double *v){
    double v0, v1;
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;

  C->cov(ZERO, next, &v0);
  C->cov(x, next, &v1);
  *v = exp(v1-v0); // == - (v0 - v1)
}

void DExp(double *x, cov_model *cov, double *v){
    double D;
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;
  C->D(x, next, &D);
  Exp(x, cov, v);
  *v *= -D;
}

void DDExp(double *x, cov_model *cov, double *v){
    double D, D2;
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;

  C->D(x, next, &D);
  C->D2(x, next, &D2);
  Exp(x, cov, v);
  *v *=  D * D + D2; // Achtung + D2 nicht - D2, da D2 Ableitung einer Cov.
}

int checkExp(cov_model *cov) {
  cov_model *next = cov->sub[0];
  int err;
  if (cov->statIn != STATIONARY && cov->statIn != VARIOGRAM)
      return ERRORSTATVARIO;
  if ((err = check2X(next, cov->tsdim, cov->xdim, cov->statIn, 
		     cov->isoIn, UNIVARIATE)) != NOERROR) return err;  
  cov->manipulating_x=true;
  setbackward(cov, next);
  cov->normalmix = false;
  return NOERROR;
}
void rangeExp(cov_model *cov, range_arraytype* ra) { }


void Pow(double *x, cov_model *cov, double *v){
    double v0, v1,
	alpha = cov->p[0][0];
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;

  C->cov(ZERO, next, &v0);
  C->cov(x, next, &v1);
  *v = v0 - pow(v0 - v1, alpha); 
//  printf("Pow %f %f %f v=%f\n", v0, v1, alpha, *v);
}

void DPow(double *x, cov_model *cov, double *v){
    double v0, v1, gamma,
	alpha = cov->p[0][0];
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;

  C->D(x, next, v);
  if (alpha == 1.0) return;
  C->cov(ZERO, next, &v0);
  C->cov(x, next, &v1);
  gamma =  v0 - v1;
  *v *= - alpha * pow(gamma, alpha -1.0);
  // Achtung  "-" , da covarianzfunktion angenommen
}

void DDPow(double *x, cov_model *cov, double *v){
    double D, v0, v1, gamma,
	alpha = cov->p[0][0];
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;

  C->D2(x, next, v);
  if (alpha == 1.0) return;
  C->D(x, next, &D);
  C->cov(ZERO, next, &v0);
  C->cov(x, next, &v1);
  gamma =  v0 - v1;
  *v *= - alpha *  pow(gamma, alpha - 2.0) * ((alpha - 1.0) * D + gamma * (*v));
  // Achtung  "-" , da covarianzfunktion angenommen
}

void invPow(double *x, cov_model *cov, double *v) {
    // invPow only defined for submodel having variance 1 !!!!
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;
  
  C->cov(x, next, v);
  double y = 1.0 - *v;
  if (y < 0.0 || y > 1.0) {
    if (y > -1e-14 && y < 0.0)  y=0.0;
    else if (y < 1.0 + 1e-14)  y=1.0;
    else {
	//PRINTF("covariance value %e (1 - %e = %e) at %e\n", *v, *v, 1.0-*v, *x);
      ERR("invPow valid only for non-negative covariance models with variance 1");
    }
  }
  
  *v = 1.0 - pow(y, 1.0 /  cov->p[0][0]);
//  printf("%f %f %f\n", y, *v, cov->p[0][0]);
//  printf("Inv %f %f %f v=%f\n", *x,  cov->p[0][0],y, *v);
}

int checkPow(cov_model *cov) {
  cov_model *next = cov->sub[0];
//  cov_fct *C = CovList + next->nr;
  int err;
  //printf("OK %d %d\n", cov->statIn, cov->isoIn);
  cov->manipulating_x=true;
 if (cov->statIn != STATIONARY && cov->statIn != VARIOGRAM)
      return ERRORSTATVARIO;
 if ((err = check2X(next, cov->tsdim, cov->xdim, cov->statIn, 
		    cov->isoIn, UNIVARIATE)) != NOERROR){
     // printf("OK %d\n",err);
      return err;  
 }
 // printf("xOK %d\n",err);
  setbackward(cov, next);      
  cov->normalmix = false;
  if (cov->statIn == STATIONARY && next->normalmix) {
      cov->normalmix = true;
  }
  return NOERROR;
}

void rangePow(cov_model *cov, range_arraytype* ra){ 
  range_type *range = ra->ranges;
  range->finiterange = true;
  range->maxdim = INFDIM;

  range->min[0] = 0.0;
  range->max[0] = 1.0;
  range->pmin[0] = 0.01;
  range->pmax[0] = 1.0;
  range->openmin[0] = true;
  range->openmax[0] = false;
}





void tbm3(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  double v1;
  cov_fct *C = CovList + next->nr;
  C->cov(x, next, v);
  C->D(x, next, &v1);
  (*v) += *x * v1 / cov->p[0][0];
}

int checktbm3(cov_model *cov) {
  cov_model *next = cov->sub[0];
  assert(next != NULL);
  int err, newdim;

  kdefault(cov, 0, 1); 
  newdim = (int) cov->p[0][0];
  if (cov->tsdim > newdim) return ERRORWRONGDIM; /// ??????
  
  if (cov->isoIn != ISOTROPIC && cov->isoIn != SPACEISOTROPIC) 
      return ERRORANISO;
  if ((err = check2X(next, 2 + newdim, cov->xdim, STATIONARY, cov->isoIn,
		     UNIVARIATE))
      != NOERROR) return err; 
  if (cov->isoIn != SPACEISOTROPIC || next->sub[0]->isoIn == SPACEISOTROPIC)
      next->tsdim = DEL_COV;
  cov->derivatives = next->derivatives - 1;
  if (cov->derivatives < 0)  {
      PrintModelInfo(cov); // ok 
      ERR("derivative for the submodel does not exist or is unknown");
  }
  cov->manipulating_x=true;
  setbackward(cov, next);
  cov->pref[CircEmbed] = next->pref[CircEmbed] + 1;
  if (cov->pref[CircEmbed] > PREF_BEST) cov->pref[CircEmbed] = PREF_BEST;
  return NOERROR;
}
void rangetbm3(cov_model *cov, range_arraytype* ra){ 
  range_type *range = ra->ranges;
  range->finiterange = true;
  range->maxdim = INFDIM;

  range->min[0] = 1.0;
  range->max[0] = RF_INF;
  range->pmin[0] = 1.0;
  range->pmax[0] = 1000000.0;
  range->openmin[0] = false;
  range->openmax[0] = true;
}



void tbm2(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;
  if (C->tbm2 == NULL) 
      ERR("tbm2 operator is unavailable for the specified submodel");
  C->tbm2(x, next, v);
}
int checktbm2(cov_model *cov){
  cov_model *next = cov->sub[0];
  assert(next != NULL);
  bool layer = cov->xdim == 2;
  int err;

  cov->nr = TBM2NR;
  cov->manipulating_x=true;
  if (cov->tsdim > 2 + layer) return ERRORWRONGDIM;
  assert(cov->xdim == 1 + layer);

  if ((err = check2X(next, 2 + layer, cov->xdim, STATIONARY, cov->isoIn,
		     UNIVARIATE)) != NOERROR) return err; 
  setbackward(cov, next);

  if (next->pref[TBM2] <= 0) {
    return ERRORCURRENTLY;
    cov->nr ++; // numerical evaluation ...
  }
  return NOERROR;
}
void rangetbm2(cov_model *cov, range_arraytype* ra){ 
  range_type *range = ra->ranges;
  range->finiterange = false;
  range->maxdim = 1;
}


typedef struct TBM2_integr{
  cov_model *cov;
  double *x;
  bool layer;
} TBM2_integr;
void TBM2NumIntegrFct(double *u,  int n, void *ex) {
  int i;
  double z[2];
  TBM2_integr *info;
  info = (TBM2_integr*) ex;
  cov_model *cov=info->cov;
  double *x = info->x;
    
  for (i=0; i<n; i++) {
    z[0] = x[0] * sqrt(1 - u[i] * u[i]);
    tbm3(z, cov, u + i);
  }
}
void tbm2num(double *x, cov_model *cov, double *v){
  TBM2_integr info;
#define MAXSUBDIVISIONS 100
  static double a = 0, 
      b = 1, 
      eps = 1e-10; // 1e-15
  static int maxsubdivisions = MAXSUBDIVISIONS, 
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


/* undefined model -- for technical reasons only */
void rangeundefined(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  range->maxdim = 0; 
  assert(false);
}
int checkundefined(cov_model *cov) {
  assert(false);
  return NA_INTEGER;
//  info->CEbad = true;
}





/* qam */
void kappaqam(int i, cov_model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = (i==0) ? cov->nsub-1 : -1;
}

void qam(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;
  int i, 
    nsub = cov->nsub;
  double sum, s, w,
    *theta = cov->p[0];
  
  sum = 0.0;
  for (i=1; i<nsub; i++) {
    cov_model * sub = cov->sub[i];
    CovList[sub->nr].cov(x, sub, &s);
    C->invSq(&s, next, &w);
    sum += theta[i - 1] * w;
  }

  sum = sqrt(sum);
  C->cov(&sum, next, v); 
}

int checkqam(cov_model *cov) {
  cov_model *next = cov->sub[0],
    *sub;
  int i, err, 
    nsub = cov->nsub,
    nsubM1 = nsub - 1;
  double v, sum;
      
  if ((err = checkkappas(cov)) != NOERROR) return err;
  
  sum = 0.0;
  for (i=0; i<nsubM1; i++) {
      sum += cov->p[0][i];
   }
  if (fabs(sum - 1.0) > 1e-14) ERR("theta must sum up to 1");    

  cov->manipulating_x=true;
  cov->normalmix = true;
  cov->finiterange = false;
  if ((err = check2X(next, 1, 1, STATIONARY, ISOTROPIC, 
		     UNIVARIATE)) != NOERROR)
    return err;
  if (!next->normalmix) ERR("phi is not completely monotone");   
 
  for (i=1; i<nsub; i++) {
    sub = cov->sub[i];
    if ((err = check2X(sub, cov->tsdim, cov->tsdim, STATIONARY, ANISOTROPIC, 
		       UNIVARIATE)) != NOERROR)
      return err;
    CovList[sub->nr].cov(ZERO, sub, &v);
    if (v != 1.0) return ERRORUNITVAR;
    setbackward(cov, sub);
  } 
  return NOERROR;
}

void rangeqam(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  
  range->min[0] = 0.0;
  range->max[0] = 1.0;
  range->pmin[0] = 0.0;
  range->pmax[0] = 1.0;
  range->openmin[0] = false;
  range->openmax[0] = false;
}



/* mqam */
void kappamqam(int i, cov_model *cov, int *nr, int *nc) {
  int nsub = cov->nsub - 1;
  *nc = (i==0) ? 1 : nsub;
  *nr = (i <= 1) ? nsub : -1;
  
  // printf("kappa %d nsub=%d %d %d\n", i, nsub, *nc, *nr);

}
void mqam(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;
  int i, j, k, l,
    nsub = cov->nsub,
    nsubM1 = nsub - 1;
  double s0,
    theta = cov->p[0][0],
    EMtheta = 1.0 - theta,
    *rho = cov->p[1],
    s[MAXSUB];
  
  for (i=1; i<nsub; i++) {
    cov_model * sub = cov->sub[i];
    CovList[sub->nr].cov(x, sub, &s0);
    C->invSq(&s0, next, s + i - 1);
  }
    
  for (j=0; j<nsubM1; j++) {
    l = k = j * nsub; // nicht nsubM1
    for (i=j; i<nsubM1; i++, k++, l+=nsubM1) {
      s0 = theta * s[i] + EMtheta * s[j];
      s0 = sqrt(s0);
      C->cov(&s0, next, v + k);
 
      v[l] = v[k] * rho[l];
      if (k != l) v[k] *= rho[k];
      assert (rho[l] == rho[k]) ;
      
      // printf("%d %d %d %d %f %f %e %e\n", i, j, k, l, 
//	     v[k],  rho[l], v[l]-v[k], rho[k]-rho[l]);
      
    }
  }

}

int checkmqam(cov_model *cov) {
//  cov_model *next = cov->sub[0];
  int err,
      nsub = cov->nsub,
      nsubM1 = nsub - 1;
      
  if ((err = checkqam(cov)) != NOERROR) return err;
  cov->vdim = nsubM1;
  cov->manipulating_x=true;
  return NOERROR;
}

void rangemqam(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges;
  double pmax = 2.0 / (double)  (cov->nsub-1);

  range->min[0] = 0.0;
  range->max[0] = 1.0;
  range->pmin[0] = 0.0;
  range->pmax[0] = pmax;
  range->openmin[0] = false;
  range->openmax[0] = false;
  
  range->min[1] = RF_NEGINF;
  range->max[1] = RF_INF;
  range->pmin[1] = -10.0;
  range->pmax[1] = 10.0;
  range->openmin[1] = true;
  range->openmax[1] = true;
}



/////////////// NATSC

void natsc(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  double scale, y;

  if (next->nr >= GATTER && next->nr <= LASTGATTER) {
      /* diese Abfrage muss bleiben, da roher aufruf ohne
	 voriges DeleteGatter aus stein_st_check heraus
	 erfolgt. 
	 Wenn Gatter gar nicht eingefuehrt wuerde, 
	 kann natsc nicht ueber checkX2 das primitive aufrufen
      */

      next = next->sub[0];
  }
  assert(CovList[next->nr].naturalscale != NULL);
  scale = CovList[next->nr].naturalscale(next);
  y = *x / scale;
  // letzteres darf nur passieren wenn dim = 1!!
  CovList[next->nr].cov(&y, next, v);
}

void Dnatsc(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
      vdimSq = cov->vdim * cov->vdim;
  double scale, y;

  if (next->nr >= GATTER && next->nr <= LASTGATTER)
      next = next->sub[0];
  assert(CovList[next->nr].naturalscale != NULL);
  scale = CovList[next->nr].naturalscale(next);
  y = *x / scale;

  CovList[next->nr].D(&y, next, v);
  for (i=0; i<vdimSq; i++) v[i] /= scale; 
}

void DDnatsc(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i,
      vdimSq = cov->vdim * cov->vdim;
  double scale, y, ScSq;

  if (next->nr >= GATTER && next->nr <= LASTGATTER)
      next = next->sub[0];
  assert(CovList[next->nr].naturalscale != NULL);
  scale = CovList[next->nr].naturalscale(next);
  y = *x / scale;
  ScSq = scale * scale;

  CovList[next->nr].D2(&y, next, v);
  for (i=0; i<vdimSq; i++) v[i] /= ScSq; 
}

void invnatsc(double *x, cov_model *cov, double *v) {
  cov_model *next = cov->sub[0];
  int i,
      vdimSq = cov->vdim * cov->vdim;
   double 
      scale = CovList[next->nr].naturalscale(next);
   
  CovList[next->nr].D(x, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= scale; 

}

void coinitnatsc(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0]->sub[0];
  if ( CovList[next->nr].coinit == NULL)
    error("# cannot find coinit -- please inform author");
  CovList[next->nr].coinit(next, li);
}

void ieinitnatsc(cov_model *cov, localinfotype *li) {
  cov_model *next = cov->sub[0]->sub[0];
 
  if ( CovList[next->nr].ieinit == NULL)
    error("# cannot find ieinit -- please inform author");
  CovList[next->nr].ieinit(next, li);
}

void tbm2natsc(double *x, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + next->nr;
  double scale = CovList[next->nr].naturalscale(next),
      y = x[0] / scale; // purely spatial , 
      
  C->tbm2(&y, next, v);
}

int checknatsc(cov_model *cov) {
  cov_model *next = cov->sub[0];
  cov_fct *C = CovList + cov->nr;
  int err;

  strcpy(ERROR_LOC, C->name);
  assert(cov->nr == NATSC);

  cov->manipulating_x = true; 

  if ((err = check2X(next, cov->tsdim, cov->xdim, cov->statIn, cov->isoIn,
		     SUBMODELM))
      != NOERROR) {
      return err;
  }

 
  if (next->sub[0]->statIn == cov->statIn &&
      next->sub[0]->isoIn == cov->isoIn) {
      next->tsdim = DEL_COV - 12;
  }
  // PrintModelInfo(cov);
//  printf(" %d %d %d %d \n", next->sub[0]->statIn, cov->statIn,
//	 next->sub[0]->isoIn, cov->isoIn); // assert(false);

  // now, next is shifted behind gatter !!:
  if (next->nr >= GATTER && next->nr <= LASTGATTER) {
      next = next->sub[0];
  } else ERR("natsc: # expected -- contact author");

  if (CovList[next->nr].naturalscale == NULL) {
      sprintf(ERRORSTRING, "natural scaling is not defined for %s",
	  CovList[next->nr].name);
      return ERRORFAILED;
  }
 
  cov->vdim = next->vdim;
  setbackward(cov, next);
  return NOERROR;
}

    
int initspectralnatsc(cov_model *cov) {
  cov_model *next = cov->sub[0];

  if (CovList[next->nr].initspectral == NULL) {
      //     PRINTF("natsc: %ld %s\n", 
//	     (POINTER) CovList[next->nr].initspectral, CovList[next->nr].name);
//      PrintModelInfo(cov);
      assert(CovList[next->nr].initspectral != NULL);
  }

  return CovList[next->nr].initspectral(next);
}

void spectralnatsc(cov_model *cov, spectral_storage *s, double *e) {
  cov_model *next = cov->sub[0];
  int d,
      dim = cov->tsdim;
  double 
      scale = CovList[next->nr].naturalscale(next);
  
  CovList[next->nr].spectral(next, s, e);
  for (d=0; d<dim; d++)  e[d] /=  scale;
}

void rangenatsc(cov_model *cov, range_arraytype* ra){ }

int initnatsc(method_type *meth){
    cov_model *cov=meth->cov,
	*next = cov->sub[0];
//  globalparam *gp = meth->gp;
  int err, i,
      origdim = meth->loc->timespacedim,
      total = origdim * meth->xdimout;
  double 
      scale = CovList[next->nr].naturalscale(next);

  // setAniso(meth); // schon im S-Teil ist die neue c-aniso gesetzt --
  //                 S-Teil darf somit nicht mit c-aniso aufgerufen werden
  //                 siehe doS, das sofort naechsten Teil aufruft.

  meth->sub[0] = (method_type*) malloc(sizeof(method_type));
  METHOD_NULL(meth->sub[0]);
  meth->nsub = 1;
  method_type *sub = meth->sub[0];
  cpyMethod(meth, sub, false);
  // cum_dollar(meth, meth->loc->timespacedim, cov, sub); :
  if (meth->caniso != NULL) {
      sub->caniso = (double *) malloc(total * sizeof(double));
      for (i=0; i<total; i++) sub->caniso[i] = meth->caniso[i] / scale;
  } else {
      sub->cscale = meth->cscale / scale;
  }
  sub->cov = cov->sub[0];

  err = initstandard(sub);

  meth->compatible = sub->compatible;
  meth->domethod = doS;
  meth->S = NULL;
  meth->nr = DOLLAR;
  return err;
}

void donatsc(method_type *meth, res_type *res){
    do_meth dometh = meth->sub[0]->domethod;
    dometh(meth->sub[0], res);   
}

////////////////////////////////////////////////////////////////////


void X(double *x, cov_model *cov, double *v){
  cov_model *sub, 
      *next = cov->sub[0];
  listoftype *list = (listoftype*) (cov->p[0]),
      *sublist;
  double var;
  int i, 
      end=cov->vdim*cov->vdim;
  
  if (LIST_ELEMENT < 0 || LIST_ELEMENT >= cov->nrow[0])
      ERR("illegal call of 'X' (it may be used only within fitvario)");

  if (cov->nrow[0] == 0) { // no X given	
      CovList[next->nr].cov(x, cov, v);
      return;
  } else {
    if (!(bool*) cov->q) ERR("model 'const' expected"); //  i.e. CONST
         
    if (CovMatrixRow >= list->nrow[LIST_ELEMENT] ||
	CovMatrixCol >= list->ncol[LIST_ELEMENT]) {
	ERR("matrices are not sufficently large");
    }

    sub = next;
    var = 1.0;
    if (sub->nr >= GATTER && sub->nr <= LASTGATTER) sub=sub->sub[0];
    if (sub->nr >= DOLLAR && sub->nr <= LASTDOLLAR) {
	var = sub->p[DVAR][0];
	sub=sub->sub[0];
	if (sub->nr >= GATTER && sub->nr <= LASTGATTER) sub=sub->sub[0];
    }
    sublist = (listoftype*) (sub->p[0]);

    if (list->ncol[LIST_ELEMENT] != sublist->nrow[LIST_ELEMENT]) 
	ERR("mixed model: X and covariance matrix M do not match");

    *v = XkCXtl(list->p[LIST_ELEMENT], sublist->p[LIST_ELEMENT], 
		list->nrow[LIST_ELEMENT], list->ncol[LIST_ELEMENT], 
		CovMatrixRow, CovMatrixCol);
    for (i=0; i<end; v[i++] *= var);
  }
}

void Xnonstat(double *x, double *y, cov_model *cov, double *v) {
  if (cov->nsub != 0 && cov->nrow[0] == 0) { // no X given, but a submodel	
      cov_model *next = cov->sub[0];
      CovList[next->nr].nonstat_cov(x, y, next, v);
      return;
  }
  X(x, cov, v);
}

// ToDO : write simulation init, do !!!


int checkX(cov_model *cov) {
  cov_model *sub, 
      *next = cov->sub[0];
  int  err,
       *nrow = cov->nrow; // taken[MAX DIM],
//      *ncol = cov->ncol,
//  listoftype *list = (listoftype*) (cov->p[0]);
  bool *Const;

  CovMatrixRow = CovMatrixCol = INT_MAX;

  cov->vdim = 1;
  if (cov->nsub!=1) error("nsub must be 1");

  if (cov->q == NULL) cov->q = (double*) malloc(sizeof(bool));
  Const = ((bool*) cov->q);
  
  //  PrintModelInfo(next);
  sub = next;
  if (sub->nr >= GATTER && sub->nr <= LASTGATTER) sub=sub->sub[0];
  if (sub->nr >= DOLLAR && sub->nr <= LASTDOLLAR) {
      sub=sub->sub[0];
      if (sub->nr >= GATTER && sub->nr <= LASTGATTER) sub=sub->sub[0];
  }
  *Const = sub->nr == CONSTANT;
  
//    printf("call: %s %s\n", CovList[next->nr].name, CovList[next->sub[0]->nr].name);
  if ((err = check2X(next, cov->tsdim, cov->xdim, cov->statIn, cov->isoIn,
		       SUBMODELM)) != NOERROR) {
      // printf("error\n");
      return err;
  }
  if (*Const) {
      next->tsdim = DEL_COV - 6;
      if (next->nr >= DOLLAR && next->nr<=LASTDOLLAR)
	  next->sub[0]->tsdim = DEL_COV - 6;
      
	//   if (ncol[0]==0)
//	ERR("if the cov. model is `const' in trend, then X must be given");
  } else {
    if (nrow[0]>0)
      ERR("if the cov. model is not `const' in trend, then X may not be given");
  }
  cov->vdim=sub->vdim; 
  setbackward(cov, next);
  
  
  if (cov->vdim > 1) 
    error("multivariate version not programmed yet");

  // incorrect. but save !!
  cov->semiseparatelast = false; // taken[xdim - 1] <= 1;
  cov->separatelast = false;     // taken[xdim - 1] <= 1; ?? 
  return NOERROR;
}


void rangeX(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;
  range->min[0] = RF_NEGINF;
  range->max[0] = RF_INF;
  range->pmin[0] = -1e10;
  range->pmax[0] = 1e10;
  range->openmin[0] = true;
  range->openmax[0] = true;

  range->finiterange = false;
}


void mixed(double *x, cov_model *cov, double *v){
    ERR("do not use 'mixed' as a component for the error -- use model 'X' instead");
}

void mixed_nonstat(double *x, double *y, cov_model *cov, double *v){
  mixed(x, cov, v);
}


void MLEmixed(double *x,  cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i, end = cov->vdim * cov->vdim;
  if (cov->nsub == 0) {
    for (i=0; i<end; i++) v[i] = 0.0; 
  } else {
      //   printf("mm %s\n", CovList[cov->nr].name);
    CovList[next->nr].cov(x, next, v);
  }
}

void MLEmixed_nonstat(double *x, double *y, cov_model *cov, double *v){
  cov_model *next = cov->sub[0];
  int i, 
      end = cov->vdim * cov->vdim;
  if (cov->nsub == 0) {
    for (i=0; i<end; i++) v[i] = 0.0; 
  } else {
     CovList[next->nr].nonstat_cov(x, y, next, v);
  }
}


int checkmixed(cov_model *cov) {
    // cov_model *sub;
  int i, err,
      nsub = cov->nsub,
      *ncol = cov->ncol,
      *nrow = cov->nrow; // taken[MAX DIM],
  char msg[50];
  listoftype *list = (listoftype*) (cov->p[0]);

  CovMatrixRow = CovMatrixCol = INT_MAX;

  cov->vdim = 1;
  if (ncol[1] > 0) { // b is given
    if (ncol[0] == 0) ERR("if b is given X must be given");
    if (nsub != 0) ERR("either b or a covariance model must be given in trend");
    for (i=0; i<nrow[0]; i++) {
      if (list->ncol[i] != nrow[1]) {
	sprintf(msg,
	      "%dth set: (%d x %d) matrix X and (%d x %d) vector b do not match",
		i, list->nrow[0], list->ncol[i], nrow[1], ncol[1]);
	ERR(msg);
      }
    }  
  } else if (nsub == 0) { // only X is given -- then only a deterministic 
	//                                 constant is given
    if (ncol[1] == 0) 
      ERR("if no covariance model is given in mixed model then b must be given");
    if (ncol[0] != 1) // deterministic effect
      ERR("X must have one column");
  } else { // submodel is given
      cov_model *sub, *next = cov->sub[0];
    bool *Const;
    if (cov->q == NULL) cov->q = (double*) malloc(sizeof(bool));
    Const = ((bool*) cov->q);
    
    //  PrintModelInfo(next);
    sub = next;
    if (sub->nr >= GATTER && sub->nr <= LASTGATTER) sub=sub->sub[0];
    if (sub->nr >= DOLLAR && sub->nr <= LASTDOLLAR) {
	sub=sub->sub[0];
	if (sub->nr >= GATTER && sub->nr <= LASTGATTER) sub=sub->sub[0];
    }
    *Const = sub->nr == CONSTANT;
    
//    printf("call: %s %s\n", CovList[next->nr].name, CovList[next->sub[0]->nr].name);
    if ((err = check2X(next, cov->tsdim, cov->xdim, cov->statIn, cov->isoIn,
		       SUBMODELM)) != NOERROR) {
        // printf("error\n");
	return err;
    }
    if (*Const) {
	next->tsdim = DEL_COV - 6;
	if (next->nr >= DOLLAR && next->nr<=LASTDOLLAR)
	    next->sub[0]->tsdim = DEL_COV - 6;

	//   if (ncol[0]==0)
//	ERR("if the cov. model is `const' in trend, then X must be given");
    } else {
      if (nrow[0]>0)
	ERR("if the cov. model is not `const' in trend, then X may not be given");
    }
    cov->vdim=sub->vdim; 
    setbackward(cov, next);
  }

  if (cov->vdim > 1) 
      error("multivariate version not programmed yet");

  // incorrect. but save !!
  cov->semiseparatelast = false; // taken[xdim - 1] <= 1;
  cov->separatelast = false;     // taken[xdim - 1] <= 1; ?? 
  return NOERROR;
}


void rangemixed(cov_model *cov, range_arraytype* ra){
  range_type *range = ra->ranges + 0;
  int i;

  for (i=0; i<=1; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = -1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }

  range->finiterange = false;
}


typedef struct mixed_storage {
  // bool genuine_dim[MAX362DIM];
  key_type key;
    double *b;
} mixed_storage;


void mixed_destruct(void ** S) {
  if (*S!=NULL) {
   mixed_storage *s = *((mixed_storage**) S);
   if (s->b != NULL) free(s->b);
   free(*S);
   *S = NULL;
  }
}


void mixed_NULL(mixed_storage* s) {
  KEY_NULL(&(s->key));
  s->b = NULL;
}

int initmixed(method_type *meth) {
  return ERRORNOTPROGRAMMED; 
  /*
  location_type *loc = meth->loc;
  mixed_storage *s;
  globalparam *gp = meth->gp;
  cov_model *cov= meth->cov,
      *next = cov->sub[0];
  char errorloc_save[nErrorLoc];
  int err;


 

  // cholesky zerlegung + bereitstellung von b
  
  SET_DESTRUCT(mixed_destruct);
  meth->domethod = dotrend;
  meth->compatible = true;

  if (next==NULL) return NOERROR;
  
  // submodel exists:
  if ((meth->S = malloc(sizeof(mixed_storage)))==0){
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }  
  s = (mixed_storage*) meth->S;
  mixed_NULL(s);
  if ((s->b = malloc(sizeof(double) * loc->totalpoints)) == 0) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  
  strcpy(errorloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "%s%s: ", errorloc_save, METHODNAMES[method]);
  err = internal_InitSimulateRF(loc->x, loc->T, loc->spatialdim,
				loc->lx, loc->grid, loc->Time,
				meth->simu.distribution, &(s->key),
				&global, 1, next);
  strcpy(ERROR_LOC, errorloc_save);
  if (err == NOERROR) return;
  
ErrorHandling: 
  return err;
  */
}



void domixed(method_type *meth, res_type *res){
  error("domixed not programmed yet");
  /*
    int i;
    cov_model *cov = meth->cov,
	 *next = cov->sub[0];
    double *b;
  */
}
    


