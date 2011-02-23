/*
 Authors 
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Simulation of a random field by circulant embedding
 (see Wood and Chan, or Dietrich and Newsam for the theory)
 and variants by Stein and by Gneiting

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2005 Yindeng Jiang & Martin Schlather
 Copyright (C) 2006 -- 2011 Martin Schlather

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
#include <stdio.h>  
#include <stdlib.h>
#include "RF.h"
#include <assert.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include <R_ext/Utils.h>     
#include <R_ext/Lapack.h> // MULT


void FFT_NULL(FFT_storage *FFT) 
{
  FFT->work = NULL;
  FFT->iwork = NULL;
}

void FFT_destruct(FFT_storage *FFT)
{
  if (FFT->iwork!=NULL) {free(FFT->iwork);}
  if (FFT->work!=NULL) {free(FFT->work);} 
  FFT_NULL(FFT);
}

void CE_destruct(void **S) 
{
  if (*S!=NULL) {
    CE_storage *x;
    x = *((CE_storage**)S);
    int l, vdimSQ = x->vdim * x->vdim;
    if (x->c!=NULL) {
      for(l=0; l<vdimSQ; l++) if(x->c[l]!=NULL) free(x->c[l]);
      free(x->c);
    }
    if (x->d!=NULL) {
      for(l=0; l<x->vdim; l++) if(x->d[l]!=NULL) free(x->d[l]);
      free(x->d);
    }
    FFT_destruct(&(x->FFT));
    if (x->aniso != NULL) free(x->aniso);
    free(*S);
    *S = NULL;
  }
}

void LOCAL_NULL(localCE_storage* x){
  // int i;  for (i=0; i<MAXCEDIM; i++) 
    x->correction = NULL;
    x->aniso = NULL;
    KEY_NULL(&(x->key));
}


void CE_NULL(CE_storage* x){
  FFT_NULL(&(x->FFT));  
  x->positivedefinite = FALSE;
  x->trials = -1;
  x->c = x->d = NULL;
  x->smallestRe = x->largestAbsIm = NA_REAL;
  x->aniso = NULL;
//  int i;
  //for (i=0; i<MAXCEDIM; i++) x->[i] = NULL;
}


void localCE_destruct(void **S) 
{
  if (*S!=NULL) {
    localCE_storage* x;
    x = *((localCE_storage**) S);
    if (x->correction != NULL) free(x->correction);
    if (x->aniso != NULL) free(x->aniso);
    KEY_DELETE(&(x->key)); 
    free(*S);
    *S = NULL;
  }
}


/*********************************************************************/
/*           CIRCULANT EMBEDDING METHOD (1994) ALGORITHM             */
/*  (it will always be refered to the paper of Wood & Chan 1994)     */
/*********************************************************************/

void SetParamCircEmbed( int *action, int *force, double *tolRe, double *tolIm,
			int *trials, 
			double *mmin, int *useprimes, int *strategy,
		        double *maxmem, int *dependent, int *method) 
{
  SetParamCE(action, force, tolRe, tolIm, trials,
	     mmin, useprimes, strategy,
	     maxmem, dependent, &(GLOBAL.ce), "CIRCEMBED");
  GLOBAL.ce.method = *method;
}

void SetParamLocal( int *action, int *force, double *tolRe, double *tolIm,
		     double *mmin, int *useprimes,
		    double *maxmem, int *dependent) 
{
  int ONE=1;
  SetParamCE(action, force, tolRe, tolIm, &ONE, 
	     mmin, useprimes, 
	     &ONE /* anything */, maxmem, dependent, &(GLOBAL.localce),
	     "LOCAL");
  GLOBAL.ce.method = 0;
//  assert(false);
}

int fastfourier(double *data, int *m, int dim, bool first, bool inverse,
		FFT_storage *FFT)
  /* this function is taken from the fft function by Robert Gentleman 
     and Ross Ihaka, in R */
{
    long int inv, nseg, n,nspn,i;
  int maxf, maxp, err;
  if (first) {
   int maxmaxf,maxmaxp;
   nseg = maxmaxf =  maxmaxp = 1;
   /* do whole loop just for err checking and maxmax[fp] .. */

   for (i = 0; i<dim; i++) {
     if (m[i] > 1) {
       fft_factor(m[i], &maxf, &maxp);
       if (maxf == 0) {err=ERRORFOURIER; goto ErrorHandling;}	
       if (maxf > maxmaxf) maxmaxf = maxf;
       if (maxp > maxmaxp) maxmaxp = maxp;
       nseg *= m[i];
     }
   }
 
   assert(FFT->work==NULL);
   if ((FFT->work = (double*) malloc(4 * maxmaxf * sizeof(double)))==NULL) { 
     err=ERRORMEMORYALLOCATION; goto ErrorHandling; 
   }
   assert(FFT->iwork==NULL);
   if ((FFT->iwork = (int*) malloc(maxmaxp  * sizeof(int)))==NULL) { 
     err=ERRORMEMORYALLOCATION; goto ErrorHandling;
   }
   FFT->nseg = nseg;
   //  nseg = LENGTH(z); see loop above
  }
  inv = (inverse) ? 2 : -2;
  n = 1;
  nspn = 1;  
  nseg = FFT->nseg;
  for (i = 0; i < dim; i++) {
    if (m[i] > 1) {
      nspn *= n;
      n = m[i];
      nseg /= n;
      fft_factor(n, &maxf, &maxp);
      fft_work(&(data[0]), &(data[1]), nseg, n, nspn, inv, FFT->work,FFT->iwork);
    }
  }
  return NOERROR;
 ErrorHandling:
  FFT_destruct(FFT);
  return err;
}

int fastfourier(double *data, int *m, int dim, bool first, FFT_storage *FFT){
  return fastfourier(data, m, dim, first, !first, FFT);
}


int init_circ_embed(method_type *meth)
{
    int err=NOERROR, d,
	*sum = NULL,
	index[MAXCEDIM], dummy, j, l, 
	*mm=NULL,
	*cumm=NULL, 
	*halfm=NULL; // MULT added vars j,l;
    double hx[MAXCEDIM], realmtot, steps[MAXCEDIM], 
      **c = NULL, // For multivariate *c->**c
      **Lambda = NULL,
      *tmp = NULL,
      *rwork=NULL, 
      *tmpLambda=NULL,
      *smallest=NULL;
  Rcomplex optim_lwork, 
      *R = NULL,
      *work=NULL;
  long int i,k,twoi,twoi_plus1, mtot, hilfsm[MAXCEDIM],
      *idx = NULL;
  bool cur_crit;
  matrix_type type;

  CE_storage *s;
  globalparam *gp = meth->gp;
  ce_param* lp = &(gp->ce);
  cov_model *cov = meth->cov;
  covfct cf = CovList[cov->nr].cov;
  int dim = cov->tsdim,
    PL = gp->general.printlevel,
    vdim   = cov->vdim,
    vdimSQ = vdim * vdim; // PM 12/12/2008
//  bool storing = gp->general.storing || cov->vdim > 1; // does not work anymore
  // 9.1.09
  location_type *loc = meth->loc;

  assert(meth->S == NULL);
  s = NULL;

  SET_DESTRUCT(CE_destruct);
  if ((meth->S = malloc(sizeof(CE_storage)))==0){
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  } 
  s = (CE_storage*) meth->S;
  CE_NULL(s);

  if (dim > MAXCEDIM) {
    err=ERRORMAXDIMMETH ; goto ErrorHandling;
  } 
  for (d=0; d<dim; d++) {
    s->nn[d]=loc->length[d]; 
    steps[d]=loc->xgr[d][XSTEP];
  }
  s->vdim = vdim;

  // PM 12/12/2008
  // in case of multiple calls: possible to use imaginary part of
  // field generated in preceding run. If s->dependent is set, fields
  // are read off in the following order:
  //  1. real field 1st square
  //  2. imag field 1st square
  //  3. real field 2nd square
  //  4. imag field 2nd square
  //  ...
  //  2k+1. real field (2k+1)th square
  //  2k+2. imag field (2k+1)th square
  //
  //  i.e. read off imag fields in even calls. Generate new ones in
  //  odd calls resp. read off real squares if s->dependent==1

  s->cur_call_odd = 0;  // first one is set to be odd at beginning of
                        // do_circ_embed

  if( (tmp = (double *) malloc(vdimSQ * sizeof(double))) == NULL) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }

  s->aniso = getAnisoMatrix(meth);
  type = meth->type;
  mtot=0;
  mm = s->m;
  halfm= s->halfm;
  cumm = s->cumm;
  if (PL>=5) PRINTF("calculating the Fourier transform\n");

  /* cumm[i+1]=\prod_{j=0}^i m[j] 
     cumm is used for fast transforming the matrix indices into an  
     index of the vector (the way the matrix is stored) corresponding 
     to the matrix                   */                               
  /* calculate the dimensions of the matrix C, eq. (2.2) in W&C */       

  //  ("CE missing strategy for matrix entension in case of anisotropic fields %d\n)
  //  for (i=0;i<dim;i++) { // size of matrix at the beginning  
  //	printf("%d %f\n", i, lp->mmin[i]);
  //  }
  //  assert(false);

  bool multivariate;
  multivariate = cov->vdim > 1;

  // multivariate = !true; // checken !!

  for (i=0;i<dim;i++) { // size of matrix at the beginning  
    double hilfsm_d;
    if (fabs(lp->mmin[i]) > 10000.0) {       
      sprintf(ERRORSTRING_OK, "maximimal modulus of mmin is 10000");
      sprintf(ERRORSTRING_WRONG,"%f", lp->mmin[i]);
      return ERRORMSG;
    }
    hilfsm_d = (double) s->nn[i];
    if (lp->mmin[i]>0.0) {
      if (hilfsm_d > (1 + (int) ceil(lp->mmin[i])) / 2) { // plus 1 since 
	// mmin might be odd; so the next even number should be used
	sprintf(ERRORSTRING_OK, "Minimum size in direction %d is %f", 
		(int) i, hilfsm_d);
	sprintf(ERRORSTRING_WRONG,"%d", (int) ceil(lp->mmin[i]));
	return ERRORMSG;
      }
      hilfsm_d = (double) ( (1 + (int) ceil(lp->mmin[i])) / 2);
    } else if (lp->mmin[i] < 0.0) {
      assert(lp->mmin[i] <= -1.0);
      hilfsm_d = ceil((double) hilfsm_d * -lp->mmin[i]);

      //	printf("%f %f\n", hilfsm_d, lp->mmin[i]);
    }
    if (hilfsm_d >=  MAXINT) return ERRORMEMORYALLOCATION;
    if (lp->useprimes) {
      // Note! algorithm fails if hilfsm_d is not a multiple of 2
      //       2 may not be put into NiceFFTNumber since it does not
      //       guarantee that the result is even even if the input is even !
	 hilfsm_d = 2.0 * NiceFFTNumber((int) hilfsm_d);
    } else {
      hilfsm_d = multivariate
	  ? round(pow(3.0, 1.0 + ceil(log(hilfsm_d) / LOG3 - EPSILON1000)))
       : pow(2.0, 1.0 + ceil(log(hilfsm_d) * INVLOG2 -EPSILON1000));
    }
    if (hilfsm_d >=  MAXINT) return ERRORMEMORYALLOCATION;
    hilfsm[i] = (int) hilfsm_d;
  }         
   
  //   assert(false);


  s->positivedefinite = false;     
  /* Eq. (3.12) shows that only j\in I(m) [cf. (3.2)] is needed,
     so only the first two rows of (3.9) (without the taking the
     modulus of h in the first row)
     The following variable `index' corresponds to h(l) in the following
way: index[l]=h[l]        if 0<=h[l]<=mm[l]/2
index[l]=h[l]-mm[l]   if mm[l]/2+1<=h[l]<=mm[l]-1     
Then h[l]=(index[l]+mm[l]) mod mm[l] !!
   */


  /* The algorithm below is as follows:
     while (!s->positivedefinite && (s->trials<lp->trials)){
     (s->trials)++;
     calculate the covariance values "c" according to the given "m"
     fastfourier(c)
     if (!lp->force || (s->trials<lp->trials)) {
     check if positive definite
     if (!s->positivedefinite && (s->trials<lp->trials)) {
     enlarge "m"
     }
     } else 
     print "forced" //
     }
   */

  s->trials=0;
  while (!s->positivedefinite && (s->trials<lp->trials)) {
    (s->trials)++;
    cumm[0]=1; 
    realmtot = 1.0;
    for(i=0; i<dim; i++) {
      if (hilfsm[i] < MAXINT) {
	mm[i] = hilfsm[i];
      } else {
	sprintf(ERRORSTRING_OK, 
	    "Each direction allows for at most %d points", MAXINT);
	sprintf(ERRORSTRING_WRONG,"%ld in the %ld-th direction", hilfsm[i], i);
	err = ERRORMSG; goto ErrorHandling;
      }
      halfm[i]=mm[i] / 2; 
      index[i]=1-halfm[i]; 
      cumm[i+1]=cumm[i] * mm[i]; // only relevant up to i<dim-1 !!
      realmtot *= (double) mm[i];
    }
    s->mtot = mtot = cumm[dim];

    if (PL>=2) {
      for (i=0;i<dim;i++) PRINTF("mm[%d]=%d, ",i,mm[i]);
      PRINTF("mtot=%d\n ", mtot  * vdimSQ);
    }

    if (realmtot * vdimSQ > (double) lp->maxmem) {
      sprintf(ERRORSTRING_OK, "%.0f", lp->maxmem);
      sprintf(ERRORSTRING_WRONG,"%.0f", realmtot * vdimSQ);
      err= ERRORCEMAXMEMORY; goto ErrorHandling;
    }


    if (c != NULL) {
      for(l=0; l<vdimSQ; l++) if(c[l]!=NULL) free(c[l]);
      free(c);
    }
    // for the following, see the paper by Wood and Chan!
    // meaning of following variable c, see eq. (3.8)
    if ((c = (double **) malloc(vdimSQ * sizeof(double *))) == NULL) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    for(l=0; l<vdimSQ; l++) {
      if( (c[l] = (double *) malloc(2 * mtot * sizeof(double))) == NULL) {
	err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
    }

    for(i=0; i<dim; i++){index[i]=0;}


    for (i=0; i<mtot; i++) {
      cur_crit = false;
      for (k=0; k<dim; k++) {
	hx[k] = steps[k] * 
	  (double) ((index[k] <= halfm[k]) ? index[k] : index[k] - mm[k]);
	cur_crit |= (index[k]==halfm[k]);
      }	
      dummy = i << 1;

      //      printf("ci %d %d %e %e \n", k, dim, hx[0], hx[1]);

      //cf(hx, cov, c + dummy);
      cf(hx, cov, tmp);

      //   assert(i < 5);

      for(l=0; l<vdimSQ; l++) {
	  c[l][dummy] = tmp[l];
	  c[l][dummy+1] = 0.0;

//	  printf("%d %d %f\n", i, mtot,  c[l][dummy]);
      }

      //    printf("%f %f %f\n", hx[0], hx[1], c[dummy]);

      if (cur_crit) {
	for (k=0; k<dim; k++) {
	  if (index[k]==halfm[k]) hx[k] -= steps[k] * (double) mm[k];
	}
	cf(hx, cov, tmp);
	for(l=0; l<vdimSQ; l++) {
	  c[l][dummy] = 0.5 *(c[l][dummy] + tmp[l]);
	}
      }
      //        printf("%10f,%10f:%10f \n", hx[0], hx[1], c[dummy]);

      k=0; while( (k<dim) && (++(index[k]) >= mm[k])) {
	index[k]=0;
	k++;
      }
      assert( (k<dim) || (i==mtot-1));

    }
    //   assert(false);


  
    bool first;
    for(j=0; j<vdim; j++) {
      for(k=0; k<vdim; k++) {
	l= j*vdim + k;
	if (false) {
	  int ii;
	  for (ii =0; ii<mtot * 2; ii++) {
	    if (ISNAN(c[l][ii]) || ISNA(c[l][ii])) {
	      assert(false);
	    }
	  }
	}

	if(k>=j) { //// martin 29.12.
		   // just compute FFT for upper right triangle
	  first = (l==0);
	  if (PL>6) PRINTF("FFT...");
	  if ((err=fastfourier(c[l], mm, dim, first, &(s->FFT)))!=0) 
	    goto ErrorHandling;
	}

	if (false) {
	  int ii;
	  for (ii =0; ii<mtot * 2; ii++) {
	    if (ISNAN(c[l][ii]) || ISNA(c[l][ii])) {
	      assert(false);
	    }
	  }
	}
      }
    }

    if (PL>6) PRINTF("finished\n");

    // here the A(k) have been constructed; A(k) =      c[0][k], c[1][k], ...,  c[p-1][k]
    //                                                  c[p][k], ...            c[2p-1][k]
    //                                                  ...                     ...
    //                                                  c[(p-1)*p][k], ...,     c[p*p-1][k]
    // NOTE: only the upper triangular part of A(k) has actually been filled; A is Hermitian.

    // Now:
    // A(k) is placed into R (cf W&C 1999)
    // application of zheev puts eigenvectors
    // of A(k) into R

    if( (R = (Rcomplex *) malloc(vdimSQ * sizeof(Rcomplex))) == NULL) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }

    if(Lambda!=NULL) {
      for(l = 0; l<vdim; l++) if(Lambda[l]!=NULL) free(Lambda[l]);
      free(Lambda);
    }
    if( (Lambda = (double **) malloc(vdim * sizeof(double *))) == NULL) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    // Lambda[l][k], l=0,...,p-1, will contain the p eigenvalues
    // of A(k) after application of zheev
  
    int lwork, info;
    if( (rwork = (double *) malloc(sizeof(double) * (3*vdim - 2) )) == NULL) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    if( (tmpLambda = (double *) malloc(sizeof(double) * vdim )) == NULL) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    for(l=0;l<vdim;l++) {
      if( (Lambda[l] = (double *) malloc(mtot * sizeof(double))) == NULL) {
	err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
    }

    if(vdim==1) { // in case vdim==1, the eigenvalues we need are already contained in c[0]
      for(i = 0; i<mtot; i++) {
	Lambda[0][i] = c[0][2*i];
	c[0][2*i] = 1.0;
	c[0][2*i+1] = 0.0;
      }
    }
    else {
      int index1, index2, sign;
      for(i = 0; i<mtot; i++) {
	twoi=2*i;
	twoi_plus1=twoi+1;

	// construct R (vdim x vdim matrix) (which is idential to A[i] (cf W&C))
	for(k=0; k<vdim; k++) {
	  for(l=0; l<vdim; l++) {
	      
	      
	    index1 = vdim * k + l;
	    if(l>=k) { index2=index1; sign=1;}  // obtain lower triangular bit of c as well with
	    else{ index2= vdim * l + k; sign=-1;}  // (-1) times imaginary part of upper triangular bit

	    R[index1].r = c[index2][twoi];
	    R[index1].i = sign*c[index2][twoi_plus1];

///	if (i < 5) 
//	    printf("%d %f %f; vdim=%d\n", index1, R[index1].r, R[index1].i, vdim);
	  }
	}

//	if (i < 5) 
//	      printf("2 : i=%d mtot=%d k=%d l=%d vdim=%d \n   ", i, mtot, k, l, vdim);

	//for(k=1; k<2; k++)  R[k].r = R[k].i = 7.0;

	// diagonalize R via Fortran call zheev
	// initialization
	if(i==0) {
	  lwork = -1; // for initialization call to get optimal size of array 'work'
	  // martin.1.1.09
	  // if (lp->method==0) { cholesky } else { SVD }


	  F77_CALL(zheev)("V", "U", &vdim, R, &vdim, tmpLambda, &optim_lwork,
			  &lwork, rwork, &info); 



           //// martin 29.12.: nur einmal aufrufen]
	  ////                                ganz weit vorne, oder haengt
	  //// optim_lwork.r  von den Werten der Matrix ab ??

	  lwork = (int) optim_lwork.r;

	  if( (work = (Rcomplex *) calloc(lwork, sizeof(Rcomplex))) == NULL ) {
	    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
	  }
	}

	F77_CALL(zheev)("V", "U", &vdim, R, &vdim, tmpLambda, work, &lwork,
			rwork, &info);

//	  assert(false);

	for(l=0;l<vdim;l++) Lambda[l][i] = tmpLambda[l];


	// martin.1.1.09 --- it only works technically!
//	F77_CALL(zpotf2)("Upper", &vdim, R, &vdim, &info);
 


	// now R contains the eigenvectors of A[i]
	// and Lambda[l][i], l=0,...,vdim-1, contains the eigenvalues of A[i]

	// the following step writes the vdim vdim-dimensional eigenvectors just computed into c


	int j; // martin 13.12.
	for(j=k=0; k<vdim; k++) {
	  int endfor = vdimSQ + k;
	  for(l=k; l<endfor; l+=vdim, j++) {
	    c[j][twoi] = R[l].r;
	    c[j][twoi_plus1] = R[l].i;
	  }
	}

	// peter 28.12. doch wieder
	if(false)
	  for(l=0; l<vdimSQ; l++) {
	    c[l][twoi] = R[l].r;
	    c[l][twoi_plus1] = R[l].i;
	  }


	// next: check if all eigenvalues are >= 0 (non-negative definiteness)
	// and do check this within this loop... in case of some lambda <0 there is no point in going
	// on diagonalizing if there are free trial... or not? Hmm.
      }
    }

    free(R);
    R = NULL;
    free(rwork);
    rwork = NULL;
    if (work!=NULL) free(work);
    work = NULL;
    free(tmpLambda);
    tmpLambda = NULL;

    // check if positive definite. If not: enlarge and restart 
    if (!lp->force || (s->trials<lp->trials)) {

      for(i=0; i<mtot; i++) { // MULT

	// Check if all eigenvalues of A(i) are > 0:
	l=0;
	while( (l<vdim) &&
	    ( s->positivedefinite=Lambda[l][i]>=lp->tol_re ) )            // eigenvalues ARE real: drop former line
	{l++;}                                                      // s->positivedefinite = (lambda.real>=cepar->tol_re) && (lambda.imag<=cepar->tol_im)

	if ( !s->positivedefinite) {
	  if (PL>=2) {
	    if(vdim==1) PRINTF("non-positive eigenvalue: Lambda[%d]=%e.\n",
		i, Lambda[l][i]);
	    else PRINTF("non-pos. eigenvalue in dim %d: Lambda[%d][%d]=%e.\n",
		l, i, l, Lambda[l][i]);
	  }
	  if (PL>=4) { // just for printing the smallest 
	    //                            eigenvalue (min(c))
 	    if( (sum = (int *) calloc(vdim, sizeof(int))) == NULL) {
	      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
	    }
	    if( (smallest = (double *)calloc(vdim, sizeof(double))) == NULL) {
	      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
	    }
	    if( (idx = (long int *) calloc(vdim, sizeof(long int))) == NULL) {
	      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
	    }
	    char percent[]="%";
	    for(j=i; j<mtot; j++) {
	      for(k=0; k<vdim; k++) {
		if(Lambda[k][j] < 0) {
		  sum[k]++;
		  if(Lambda[k][j] < smallest[k]) {
		    idx[k] = j;
		    smallest[k] = Lambda[k][j];}
		}
	      }
	    }
	    if(vdim==1) {
	      PRINTF("   There are %d negative eigenvalues (%5.3f%s).\n",
		  sum[0], 100.0 * (double) sum[0] / (double) mtot, percent);
	      PRINTF("   The smallest eigenvalue is Lambda[%d] =  %e\n",
		  idx[0], smallest[0]);
	    }
	    else {
	      for(l=0; l<vdim; l++) {
		PRINTF("   There are %d negative eigenvalues in dimension %d (%5.3f%s).\n",
		    sum[l], l, 100.0 * (double) sum[l] / (double) mtot, percent);
		PRINTF("   The smallest eigenvalue in dimension %d is Lambda[%d][%d] =  %e\n",
		    l, idx[l], l, smallest[l]);
	      }
	    }
	    free(sum);
	    sum = NULL;
	    free(smallest);
	    smallest = NULL;
	    free(idx);
	    idx = NULL;
	  }
	  break; // break loop 'for(i=0 ...)'
	}
      }

      if (!s->positivedefinite && (s->trials<lp->trials)) { 
	FFT_destruct(&(s->FFT));
	switch (lp->strategy) {
	  case 0 :
	    for (i=0; i<dim; i++) { /* enlarge uniformly in each direction, maybe 
				       this should be modified if the grid has 
				       different sizes in different directions */
		hilfsm[i] *= multivariate ? 3 : 2;
	    }
	    break;
	  case 1 :  
	    double cc, maxcc, hx[MAXCEDIM];

	    int maxi;
	    maxi = -1;
	    maxcc = RF_NEGINF;
	    for (i=0; i<dim; i++) hx[i] = 0.0;
	    for (i=0; i<dim; i++) {

	      // MUSS HIER ETWAS GEAENDERT WERDEN FUER MULTIVARIAT?

	      hx[i] = steps[i] * mm[i]; 

	      cc=0;
	      cf(hx, cov, tmp);
	      for(l=0; l<vdimSQ; l++) cc += fabs(tmp[l]);
	      if (PL>2) PRINTF("%d cc=%e (%e)",i,cc,hx[i]);
	      if (cc>maxcc) {
		maxcc = cc;
		maxi = i; 
	      }
	      hx[i] = 0.0;
	    }
	    assert(maxi>=0);
	    hilfsm[maxi] *= multivariate ? 3 : 2;
	    break;
	  default:
	    assert(false);
	}
      }
    } else {if (PL>=2) PRINTF("forced\n");}
    R_CheckUserInterrupt();
  } // while (!s->positivedefinite && (s->trials<lp->trials)) 354-719


  assert(mtot>0);
  if (s->positivedefinite || lp->force) { 

    // martin.1.1.09
    // if (lp->method == 0) { the following, with svd !} else { failed }

    // correct theoretically impossible values, that are still within 
    // tolerance CIRCEMBED.tol_re/CIRCEMBED.tol_im 
    double r;
    r = Lambda[0][0]; // MULT

    for(i=0,twoi=0;i<mtot;i++) {
      twoi_plus1 = twoi+1;
      for(l=0;l<vdim;l++) {
	if(Lambda[l][i] > 0.0) Lambda[l][i] = sqrt(Lambda[l][i]);
	else {
	  if(Lambda[l][i] < r) r = Lambda[l][i];
	  Lambda[l][i] = 0;
	}
      }

      // now compute R[i]*Lambda^(1/2)[i]
      for(j=0;j<vdim;j++) {
	for(l=0; l<vdim; l++) {
	  c[j*vdim+l][twoi] *= Lambda[l][i];
	  c[j*vdim+l][twoi_plus1] *= Lambda[l][i];
	}
      }

      twoi+=2;

    }
    if (PL>1) {
      if (r < -GENERAL_PRECISION) {
	PRINTF("using approximating circulant embedding:\n");
	PRINTF("\tsmallest real part has been %e \n", r);
      }
    }
    s->smallestRe = (r > 0.0) ? NA_REAL : r;
    s->largestAbsIm = 0.0; // REMOVE THIS
  } else {
    err = ERRORCIRCNONPOS;
    goto ErrorHandling;
  }
  if (PL >= 10) {
    for (i=0; i < 2 * mtot; i++)
      for (l=0; l<vdimSQ; l++)
      {PRINTF("%f ",c[l][i]);} PRINTF("\n");
  }  


  s->dependent = lp->dependent;
  s->new_simulation_next = true;
  for(i=0; i<dim; i++) { // set anyway -- does not cost too much
    s->cur_square[i] = 0;
    s->max_squares[i] = hilfsm[i] / s->nn[i];
    s->square_seg[i] = cumm[i] * 
      (s->nn[i] + (hilfsm[i] - s->max_squares[i] * s->nn[i]) / 
       s->max_squares[i]);
  }

  if ( (s->d=(double **) calloc(vdim, sizeof(double *))) == NULL ) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  for(l=0; l<vdim;l++) {
      if( (s->d[l] = (double *) calloc(2 * mtot, sizeof(double))) == NULL) {
	  err=ERRORMEMORYALLOCATION;goto ErrorHandling;
      }
  }

  s->c = c;

 
ErrorHandling:
  if (Lambda != NULL) {
    for(l = 0; l<vdim; l++) if(Lambda[l] != NULL) free(Lambda[l]);
    free(Lambda);
  }

  if (R != NULL) free(R);
  if (rwork != NULL) free(rwork);
  if (work!=NULL) free(work);
  if (tmpLambda != NULL)   free(tmpLambda);
  if (sum != NULL) free(sum);
  if (smallest != NULL) free(smallest);
  if (idx != NULL) free(index);
  if (tmp != NULL) free(tmp);
  
  if (err != NOERROR) {
      if (s != NULL) s->c = c;
      else if (c!=NULL) {
	  for(l=0;l<vdimSQ;l++) if(c[l]!=NULL) free(c[l]);
	  free(c);
      }
  }

  return err;
}



void do_circ_embed(method_type *meth, res_type *res){
  cov_model *cov = meth->cov;
  int  i, j, k, l, pos,HalfMp1[MAXCEDIM], HalfMaM[2][MAXCEDIM], index[MAXCEDIM], 
    dim, *mm, *cumm, *halfm, // MULT added vars l, pos
    vdim;
  double invsqrtmtot,
	 **c=NULL,
	 **d=NULL,
	 *dd=NULL; // MULT *c->**c
  bool vfree[MAXCEDIM+1], noexception; // MULT varname free->vfree
  long mtot, start[MAXCEDIM], end[MAXCEDIM];
  CE_storage *s;
//  globalparam *gp = meth->gp;
//  ce_param *lp = &(gp->ce);
  location_type *loc = meth->loc;
  simu_type *simu = meth->simu;

  s = (CE_storage*) meth->S;
 
/* 
     implemented here only for rotationsinvariant covariance functions
     for arbitrary dimensions;
     (so it works only for even covariance functions in the sense of 
     Wood and Chan,p. 415, although they have suggested a more general 
     algorithm;) 
     
*/  
  
  dim = cov->tsdim;
  mm = s->m;
  cumm = s->cumm;
  halfm = s->halfm;
  c = s->c;
  d = s->d;
  mtot= s->mtot;
  vdim = s->vdim; // so p=vdim

  for (i=0; i<dim; i++) {
    HalfMp1[i] =  ((mm[i] % 2)==1) ? -1 : halfm[i];
    HalfMaM[0][i] = halfm[i];
    HalfMaM[1][i] = mm[i] - 1;
  }

  //if there are doubts that the algorithm below reaches all elements 
  //of the matrix, uncomment the lines containing the variable xx
  //
  //bool *xx; xx = (bool*) malloc(sizeof(bool) * mtot);
  //for (i=0; i<mtot;i++) xx[i]=true;
  
  invsqrtmtot = 1/sqrt((double) mtot);

  if (PL>=10) PRINTF("Creating Gaussian variables... \n");
  /* now the Gaussian r.v. have to defined and multiplied with sqrt(FFT(c))*/

  for (i=0; i<dim; i++) {
    index[i]=0; 
    vfree[i] = false;
  }
  vfree[MAXCEDIM] = false;

  // PM 12/12/2008
  // don't simulate newly in case this is an even call
  s->cur_call_odd = !(s->cur_call_odd);
  s->new_simulation_next = (s->new_simulation_next && s->cur_call_odd);

  if (s->new_simulation_next) {

    // MULT
    Rcomplex *gauss1, *gauss2;
    gauss1 = (Rcomplex *) malloc(sizeof(Rcomplex) * vdim);
    gauss2 = (Rcomplex *) malloc(sizeof(Rcomplex) * vdim);

    
    // martin: 13.12. check
    //for(k=0; k<vdim; k++) {
    //      dd = d[k];
    //  for (j=0; j<mtot; j++) dd[j] = RF_NAN;
    //    }

    
    for(;;) {      
      i = j = 0;
      noexception = false;
      for (k=0; k<dim; k++) {
	i += cumm[k] * index[k];
	if ((index[k]==0) || (index[k]==HalfMp1[k])) {
	  j += cumm[k] * index[k];
	} else {
	  noexception = true; // not case 2 in prop 3 of W&C
	  j += cumm[k] * (mm[k] - index[k]);
	}
      }
      
//A 0 0 0 0 0 0.000000 0.000000
//B 0 0 0 0 0 -26.340334 0.000000

//    if (index[1] ==0) printf("A %d %d %d %d %d %f %f\n", 
//	   index[0], index[1], i, j, noexception, d[i], d[i+1]);
//      if (PL>=10) PRINTF("cumm...");
      i <<= 1; // since we have to index imaginary numbers
      j <<= 1;

//      printf("ij %d %d (%d)\n", i, j, noexception);

      if (noexception) { // case 3 in prop 3 of W&C

        for(k=0; k<vdim; k++) {
            gauss1[k].r = GAUSS_RANDOM(1.0);
            gauss1[k].i = GAUSS_RANDOM(1.0);
            gauss2[k].r = GAUSS_RANDOM(1.0);
            gauss2[k].i = GAUSS_RANDOM(1.0);
        }

	
        for(k=0; k<vdim; k++) {
	  dd = d[k];
	  dd[j] = dd[i] = dd[j+1] = dd[i+1] = 0.0;// martin:13.12. neu eingefuegt
	  for(l=0; l<vdim; l++) {
	    pos = k * vdim + l;
	    dd[i] += c[pos][i] * gauss1[l].r - c[pos][i+1] * gauss1[l].i;   
      // real part
	    dd[i+1] += c[pos][i+1] * gauss1[l].r + c[pos][i] * gauss1[l].i; 
      // imaginary part 
	    
	    dd[j] += c[pos][j] * gauss2[l].r - c[pos][j+1] * gauss2[l].i;  
       // real part    // martin: 13.12. geaendert von = + zu +=
	    dd[j+1] += c[pos][j+1] * gauss2[l].r + c[pos][j] * gauss2[l].i;  
	    // imaginary part
	  }
        }

      } else { // case 2 in prop 3 of W&C

        for(k=0; k<vdim; k++) {
            gauss1[k].r = GAUSS_RANDOM(1.0);
            gauss1[k].i = GAUSS_RANDOM(1.0);
        }

        for(k=0; k<vdim; k++) {
	  dd = d[k];
	  dd[i+1] = dd[i] = 0.0; //  martin: 13.12. neu eingefuegt
	  for(l=0; l<vdim; l++) {
	    pos = k * vdim + l;
	    dd[i] += c[pos][i] * gauss1[l].r - c[pos][i+1] * gauss1[l].i;   
	    // real part
	    dd[i+1] += c[pos][i+1] * gauss1[l].r + c[pos][i] * gauss1[l].i;  
	    // imaginary part
	  }
        }	
      }

      
//    if (i==224|| i==226 || i==228) {
//printf("B %d %d %d %d %d %f %f\n", 
//	   index[0], index[1], i, j, noexception, d[i], d[i+1]);
// assert(false); 
//    }
      
//      if (PL>=10) PRINTF("k=%d ", k);
      
//    if (index[1] == 0) printf("B %d %d %d %d %d %f %f\n", 
//	   index[0], index[1], i, j, noexception, d[i], d[i+1]);
      /* 
	 this is the difficult part.
	 We have to run over roughly half the points, but we should
	 not run over variables twice (time lost)
	 Due to case 2, we must include halfm.
	 
	 idea is:
	 for (i1=0 to halfm[dim-1])
           if (i1==0) or (i1==halfm[dim-1]) then endfor2=halfm[dim-2] 
           else endfor2=mm[dim-2]
         for (i2=0 to endfor2)
	   if ((i1==0) or (i1==halfm[dim-1])) and
	      ((i2==0) or (i2==halfm[dim-2]))
	    then endfor3=halfm[dim-3] 
	    else endfor3=mm[dim-3]
	 for (i3=0 to endfor3)
	 ....
	 
	 i.e. the first one that is not 0 or halfm (regarded from dim-1 to 0)
	 runs over 0..halfm, all the others over 0..m
	 
	 this is realised in the following allowing for arbitrary value of dim
	 
	 vfree==true   <=>   endfor==mm[]	 
      */
      k=0; 
      if (++index[k]>HalfMaM[vfree[k]][k]) {
	// in case k increases the number of indices that run over 0..m increases
	vfree[k] = true;
	index[k]= 0; 
	k++;
	while((k<dim) && (++index[k]>HalfMaM[vfree[k]][k])) {
	  vfree[k] = true;
	  index[k]= 0; 
	  k++;
	}
	if (k>=dim) break;
	// except the very last (new) number is halfm and the next index is
	// restricted to 0..halfm
	// then k decreases as long as the index[k] is 0 or halfm
	if (!vfree[k] && (index[k]==halfm[k])){//index restricted to 0..halfm?
	  // first: index[k] is halfm? (test on ==0 is superfluent)
	  k--;
	  while ( (k>=0) && ((index[k]==0) || (index[k]==halfm[k]))) {
	    // second and following: index[k] is 0 or halfm?
	    vfree[k] = false;
	  k--;
	  }
	}
      }
    }
    
    // MULT
    if(gauss1!=NULL) free(gauss1);
    if(gauss2!=NULL) free(gauss2);

//  double zz=0.0;
//  for (i =0; i<2 * mtot; i++) {
//    if (true || i<=100) printf("%d %d\n", i, 2 * mtot);
//    zz += d[i];
//     assert(i!=230);
//     // if (i>100); break;
//  }
// printf("%f\n", zz);

    // MULT
    for(k=0; k<vdim; k++) {
      fastfourier(d[k], mm, dim, false, &(s->FFT));
    }
  } // if ** simulate **



  /* now we correct the result of the fastfourier transformation
     by the factor 1/sqrt(mtot) and read the relevant matrix out of 
     the large vector c */
  for(i=0; i<dim; i++) {
    index[i] = start[i] = s->cur_square[i] * s->square_seg[i];
    end[i] = start[i] + cumm[i] * s->nn[i];
//    printf("%d start=%d end=%d cur=%d max=%d seq=%d cumm=%d, nn=%d\n",
//	   i, (int) start[i], (int) end[i],
//	   s->cur_square[i], s->max_squares[i],
//	   (int) s->square_seg[i], cumm[i], s->nn[i]);
  }
  int totpts = loc->totalpoints; // MULT

/*
  for(l=0; l<vdim; l++) { // MULT
    for(k=0; k<dim; k++) index[k] = start[k];

    for (i=l*totpts; i<(l+1)*totpts; i++){
      j=0; for (k=0; k<dim; k++) {j+=index[k];}
      //res[i] += d[l][2 * j + !(s->cur_call_odd)] * invsqrtmtot;
      res[i] += d[l][2 * j] * invsqrtmtot;
      for(k=0; (k<dim) && ((index[k] += cumm[k]) >= end[k]); k++) {
        index[k]=start[k];
      }
    }
  }
*/

  for(k=0; k<dim; k++) index[k] = start[k];

  for (i=0; i<totpts; i++) {
    j=0; for (k=0; k<dim; k++) {j+=index[k];}
    for (l=0; l<vdim; l++) { // MULT
      res[vdim*i+l] += (res_type) 
	  (d[l][2 * j + !(s->cur_call_odd)] * invsqrtmtot);
      
//      if (false)
//	if ( !( res[vdim*i+l] < 10) || i==0 && l==0)
//	    printf("!! %d %d %d %f -- %d %f %f\n",
//		 i, l, vdim, res[vdim*i+l],
//		 2 * j + !(s->cur_call_odd),
//		 d[l][2 * j + !(s->cur_call_odd)],
//		 invsqrtmtot);
      
      //res[vdim*i+l] += d[l][2 * j] * invsqrtmtot;
      
      //PRINTF("2*j+..=%d\n", 2 * j + !(s->cur_call_odd));
      //PRINTF("i=%d, l=%d, assigning %f to res[%d]\n", 
      //         i, l, d[l][2 * j] * invsqrtmtot, vdim*i+l);
    }
    for(k=0; (k<dim) && ((index[k] += cumm[k]) >= end[k]); k++) {
      index[k]=start[k];
      
    }
  }

  s->new_simulation_next = true;
  if (s->dependent && !(s->cur_call_odd)) {
  //if (s->dependent) {
//  if MaxStableRF calls Cutoff, c and look uninitialised
    k=0; 
    while(k<dim && (++(s->cur_square[k]) >= s->max_squares[k])) {
      s->cur_square[k++]=0;
    }
    s->new_simulation_next = k==dim; 
  }
  //printf("%d %d\n", s->cur_square[0],s->cur_square[1]);
  simu->stop |= s->new_simulation_next && s->d == NULL;
}





int GetOrthogonalUnitExtensions(double * aniso, int dim, double *grid_ext) {
  int k,i,j,l,m, job=01, err, dimsq, ev0, jump, endfor;
  double *s=NULL, G[MAXCEDIM+1], e[MAXCEDIM], D[MAXCEDIM], *V=NULL;
  dimsq = dim * dim;

  assert(aniso != NULL);
  
  s = (double*) malloc(dimsq * sizeof(double));
  V = (double*) malloc(dimsq * sizeof(double));
  for (k=0; k<dim; k++) {
    // kte-Zeile von aniso wird auf null gesetzt/gestrichen -> aniso_k
    // s = aniso %*% aniso_k 
    // s hat somit mindestens Eigenwert der 0 ist.
    //der zugeheorige EVektor wird miit der k-ten Spalte von aniso multipliziert
    // und ergibt dann den Korrekturfaktor
    for (i=0; i<dim; i++) {
      for (l=j=0; j<dimsq; j+=dim) {
	s[j + i] = 0.0;
	jump = l + k;
	endfor = l + dim;
	for (m=i; l<endfor; l++, m += dim) {
            if (l!=jump) s[i + j] +=  aniso[m] * aniso[l];
	}
      }
    }
//    for (printf("\n\n s\n"), j=0; j<dim; j++, printf("\n")) 
//	for (i=0; i<dimsq; i+=dim) printf("%f ", s[i+j]);
//    for (printf("\n\n aniso\n"), j=0; j<dim; j++, printf("\n")) 
//	for (i=0; i<dimsq; i+=dim) printf("%f ", aniso[i+j]);	
    
    // note! s will be distroyed by dsvdc!
    F77_NAME(dsvdc)(s, &dim, &dim, &dim, D, e, NULL /* U */,
		    &dim, V, &dim, G, &job, &err);
    if (err!=NOERROR) { err=-err;  goto ErrorHandling; }
    ev0 = -1;
    for (i=0; i<dim; i++) {
      if (fabs(D[i]) <= EIGENVALUE_EPS) {
	if (ev0==-1) ev0=i;
	else {
	    err = ERRORFULLRANK; goto ErrorHandling;
	} 
      }
    }
    
//    for (printf("\n\n V\n"), j=0; j<dim; j++, printf("\n")) 
//	for (i=0; i<dimsq; i+=dim) printf("%f ", V[i+j]);

    grid_ext[k] = 0.0;
    ev0 *= dim;
    for (i=0; i<dim; i++) {
//	printf("%d : %f %f %f\n", k,
//	       V[ev0 + i], aniso[k + i * dim], V[ev0 + i] / aniso[k + i * dim]);
      grid_ext[k] += V[ev0 + i] * aniso[k + i * dim];
    }
    grid_ext[k] = fabs(grid_ext[k]);
  }
  free(V);
  free(s);
  return NOERROR;

 ErrorHandling:
  if (err<0) {
    PRINTF("F77 error in GetOrthogonalExtensions: %d\n", -err);
    err = ERRORFAILED;
  }
  free(V);
  free(s);
  return err;
}

//void GOUE(double * aniso, int *dim, double *grid_ext) {
//    int err = GetOrthogonalUnitExtensions(aniso, *dim, grid_ext);
//    // if (err) printf("ERROR!");
//}



int init_circ_embed_local(method_type *meth, SimulationType method){
  // being here, leadin $ have already been put away
  // so a new $ is included below
 
  cov_model *dummy = NULL,
    *truecov = NULL,
    *cov = meth->cov;
  location_type *loc = meth->loc;
  simu_type *simu = meth->simu;
  globalparam *gp = meth->gp;
  int instance, i, d, dimM1, bytes,
    timespacedim = loc->timespacedim,   
    newcol,
    local_nr[Nothing + 1], 
    PL = gp->general.printlevel,
//    nrow = timespacedim,
    cncol = cov->tsdim,
    cxdim = cov->xdim,
    rowcol = timespacedim * cncol,
     err = NOERROR;
  localCE_storage *s=NULL;
  double diameter, grid_ext[MAXCEDIM],
      min[MAXCEDIM], max[MAXCEDIM];
  cov_model *localcov=NULL;
//  ce_param* lp = &(gp->ce);
  globalparam newparam;

  char co_name[]="cutoff", ie_name[]="Stein";
  local_nr[CircEmbedCutoff] = getmodelnr(co_name);
  local_nr[CircEmbedIntrinsic] = getmodelnr(ie_name);

//  static pref_type pref =
//    {5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  // CE CO CI TBM23 Sp di sq Ma av n mpp Hy - any

  SET_DESTRUCT(localCE_destruct);
  assert(meth->S==NULL);
  if ((meth->S=malloc(sizeof(localCE_storage)))==0) 
    return ERRORMEMORYALLOCATION; 
  s = (localCE_storage*) meth->S;
  LOCAL_NULL(s);

//  if (meth->caniso != NULL || meth->cproj!=NULL || meth->cscale != 1.0)
//      return ERRORPREVDOLLAR;
  newcol = cov->tsdim; 


// printf("local circulant embedding currently does not work\n");
// err = ERRORFAILED; //see assert(false); !!! below
//   goto ErrorHandling;


  if (cov->nr >= GATTER && cov->nr <= LASTGATTER) {
    if (cov->nr != S2ISO && cov->nr != ISO2ISO) {
      return ERRORUNKNOWNMETHOD; // besser machen
    }
    cov = cov->sub[0];
  }

   method_type *newmeth=NULL; 
  {
    method_type *oldmeth=NULL;
    oldmeth = (method_type*) malloc(sizeof(method_type));
    METHOD_NULL(oldmeth);
    cpyMethod(meth, oldmeth, true);
    //   printf("oldmeth xdim out %d\n", meth->xdimout);
    //   printf("oldmeth xdim out %d\n", oldmeth->xdimout);
    while (cov->nr >= DOLLAR && cov->nr <= LASTDOLLAR) {
      newmeth = (method_type *) malloc(sizeof(method_type));
      METHOD_NULL(newmeth);
      cum_dollar(oldmeth, timespacedim, cov, newmeth);
      METHOD_DELETE(&oldmeth);
      oldmeth = newmeth;

      //    printf("xdim out %d\n", oldmeth->xdimout);

      cov = cov->sub[0];
      if (cov->nr >= GATTER && cov->nr <= LASTGATTER) {
	  if (cov->nr != S2ISO && cov->nr != ISO2ISO) {
	      err = ERRORUNKNOWNMETHOD; // besser machen
	      goto ErrorHandling;
	  }
	  cov = cov->sub[0];
      }
    }
    newmeth = oldmeth;  
  }
//    printf("newmeth xdim out %d\n", newmeth->xdimout);
//  assert(false);


  if (cov->vdim > 1) {
    err = ERRORNOMULTIVARIATE;
    goto ErrorHandling;
  }
  switch(method) {
      case CircEmbedCutoff :
	  if (CovList[cov->nr].coinit == NULL) {
	      err = ERRORUNKNOWNMETHOD;
	      goto ErrorHandling;
	  }
	break;
      case CircEmbedIntrinsic :
	  if (CovList[cov->nr].ieinit == NULL) {
	      err =  ERRORUNKNOWNMETHOD;
	      goto ErrorHandling;
  	  }
	break;
      default : assert(false);
  }

//  PrintModelInfo(cov);
//  assert(cov->xdim == cov->tsdim);
  dimM1 = cov->xdim + 1;
 
  memcpy(&newparam, meth->gp, sizeof(globalparam));
  memcpy(&newparam.ce, &(newparam.localce), sizeof(ce_param));

  newmeth->loc = meth->loc;
 

//  printf("%d\n", newmeth->caniso); 
//  for (i=0; i<100; i++) printf("%d %f\n", i, newmeth->caniso[i]); 
//  assert(false);

  s->aniso = getAnisoMatrix(newmeth);
 
  {
      double  origmin[MAXCEDIM], origmax[MAXCEDIM], center[MAXCEDIM], dummy[MAXCEDIM];
      GetMinMax(meth, origmin, origmax, center, MAXCEDIM);
      diameter = GetDiameter(center, origmin, origmax, s->aniso, 
			     loc->timespacedim, newcol, 
			     min, max, dummy);
  }

  if (PL>4) PRINTF("diameter %f\n", diameter);

  err =  GetOrthogonalUnitExtensions(s->aniso, timespacedim, grid_ext);
  if (err != NOERROR) {
    goto ErrorHandling;
  }

  //    assert(false);

  covcpy(&localcov,
	 (cov->nr >= GATTER && cov->nr <= LASTGATTER) ? cov->sub[0] : cov, 
	 true, false); 
  
  addModel(&localcov, local_nr[method]);    
  
  localcov->ncol[pLOC_DIAM] = localcov->nrow[pLOC_DIAM]  = 1;
  localcov->p[pLOC_DIAM] = (double*) calloc(1, sizeof(double));
  truecov = localcov;
  
  //  assert(false);
  addModel(&localcov, GATTER);
  //  assert(false);

  addModel(&localcov, DOLLAR);
  localcov->tsdim = 
    localcov->xdim = timespacedim;

  if (cncol > 1 || newcol > 1 || true) {
      localcov->nrow[DANISO] = cncol;
      localcov->ncol[DANISO] = newcol;
      bytes = sizeof(double) * localcov->nrow[DANISO] * localcov->ncol[DANISO];
      localcov->p[DANISO] = (double*) calloc(1, bytes);
      memcpy(localcov->p[DANISO], s->aniso, bytes);
  } else {
      assert(cncol == 1 && newcol == 1);
      localcov->nrow[DSCALE] = localcov->ncol[DSCALE] = 1;
      localcov->p[DSCALE] = (double*) malloc(sizeof(double));
      localcov->p[DSCALE][0] = 1.0 / s->aniso[0];
  }


  // PrintModelInfo(localcov);
  // printf("hier geht was schief\n");
  // assert(false);
    
  kdefault(localcov, DVAR, meth->cvar);

  //printf("%d %d %d %d\n", cncol, cxdim, cov->xdim, newcol); //assert(false);

  addModel(&localcov, GATTER);
  localcov->calling=NULL;
  localcov->tsdim = cncol;
  localcov->xdim = cxdim;
  localcov->statIn = STATIONARY;

//  PrintModelInfo(localcov); printf("\n\n\n\n");
  //PrintModelInfo(localcov); assert(false);
//
//  PrintModelInfo(localcov); assert(false);

  truecov->p[pLOC_DIAM][0] = diameter;
  err = CovList[localcov->nr].check(localcov);
  DeleteGatter(&localcov);

//  PrintModelInfo(localcov); assert(false);

  dummy = localcov;

  for (;dummy!=NULL;) {
    for (i=0; i<Forbidden; i++) 
      dummy ->pref[i] = dummy->user[i] = PREF_NONE;
    dummy->pref[CircEmbed] = localcov->user[CircEmbed] = PREF_BEST;
    if (! ((dummy->nr >= GATTER && dummy->nr <= LASTGATTER)
	   ||
	   (dummy->nr >= DOLLAR && dummy->nr <= LASTDOLLAR)
	  )) break;
    dummy = dummy->sub[0];
  }

//PrintModelInfo(localcov); assert(false);
  if (err < MSGLOCAL_OK && err != NOERROR) goto ErrorHandling;
//  printf("circ err %d %d\n", err, MSGLOCAL_OK);
  
  int first_instance;
  first_instance = err != NOERROR;
  for (instance = first_instance; instance < 2; instance++) {
    KEY_DELETE(&(s->key));
    KEY_NULL(&(s->key));

//    printf("xxx %f %d\n", truecov->q[LOCAL_MSG], MSGLOCAL_OK);
    //   printf("inst %d\n, ", instance);
    
    if (instance == 1) {
      if (truecov->q[LOCAL_MSG] != MSGLOCAL_OK) {
	if (!CovList[truecov->nr].alternative(truecov)
	    // || (Checking(cov) != NOERROR)
	  ) break;
      } else {
	  if (first_instance == 0 && err != NOERROR) goto ErrorHandling;
	  else assert(false);
      }
    } else assert(instance == 0);

    for (d=0; d<timespacedim; d++) {
      if (newparam.ce.mmin[d]==0.0) {
//	double diff = max[d] - min[d];
	newparam.ce.mmin[d] = - truecov->q[LOCAL_R] / 
	  (grid_ext[d] * (double) (loc->length[d] - 1) * loc->xgr[d][XSTEP]) ;

//	printf("%d %f %f \n", d, newparam.ce.mmin[d], truecov->q[LOCAL_R]);

//	  ((double) (loc->length[d] - 1) * loc->xgr[d][XSTEP]);
	if (newparam.ce.mmin[d] > -1.0) newparam.ce.mmin[d] = -1.0;
//      printf("%d mmin=%f local=%f len=%d ste[=%f, diff=%f\n", 
//	     d, newparam.ce.mmin[d], truecov->q[LOCAL_R], 
//	    loc->length[d], loc->xgr[d][XSTEP], diff
//	     ); 
      }
    }
    
    //   assert(false);

/*
    if (PL>7) {  
      int i; 
      PRINTF("cov->q\n", cov->qlen); for(i=0;i<9;i++) PRINTF("%f ", cov->q[i]); PRINTF("\n");  assert(false);}
  

      PRINTF("nc=%d inst=%d err=%d hyp.kappa=%f, #=%d, diam=%f\n r=%f curmin_r=%f\n", 
	     v, nc, instance, Xerror, sc->param[HYPERKAPPAII],
	     (int) sc->param[HYPERNR],
	     sc->param[DIAMETER],  sc->param[LOCAL_R],
	     store_param[LOCAL_R]);
      if (hyper->implemented[CircEmbedIntrinsic] == HYPERIMPLEMENTED)
	PRINTF("lcoalr=%f, msg=%d, a0=%f, a2=%f, b=%f\n",
	       cov->q[LOCAL_R],
	       (int) cov->q[LOCAL_MSG],
	       cov->q[INTRINSIC_A0], cov->q[INTRINSIC_A2],
	       cov->q[INTRINSIC_B]
	  );
*/
//    printf("var nach internel%ld\n", &newparam);//assert(false);
    //   PrintModelInfo(localcov);
//    PrintModelInfo(truecov);
    //   assert(false);

    cov_model *keycov=NULL;
    covcpy(&keycov, localcov, false, true); 
 
//   PrintModelInfo(keycov);
//   cov_model *Stein=keycov->sub[0]->sub[0]->sub[0];
//   printf("%d %f %f %f\n", err, Stein->q[INTRINSIC_A0], Stein->q[INTRINSIC_A2], Stein->q[INTRINSIC_B]);  // 0.500000 0.003086 0.000000
// assert(false);

//   PrintModelInfo(keycov); //assert(false);
//   printf("GO !!\n\n\n");

    err = internal_InitSimulateRF(loc->xgr[0], loc->T, 
				  timespacedim - (loc->Time),
				  3, true, loc->Time,
				  DISTR_GAUSS, &(s->key), 
				  &newparam, simu->expected_number_simu,
				  &keycov); // localcov is moved to s->key

    //   printf("ci err %d\n", err);
    //   warning(" assert(false);");

    if (err == NOERROR) break;
  }
  if (err != NOERROR) goto ErrorHandling;

  double sqrt2a2, *stein_aniso;
  if (method == CircEmbedIntrinsic) {
    sqrt2a2 = sqrt(2.0 * truecov->q[INTRINSIC_A2]); // see Stein (2002)
    if ((s->correction = malloc(sizeof(double) * rowcol))==NULL){
      err = ERRORMEMORYALLOCATION;
      goto ErrorHandling;
    }
    stein_aniso = (double*) s->correction;
    for (i=0; i<rowcol; i++)  stein_aniso[i] = sqrt2a2 * s->aniso[i];
  } else {
    assert(method == CircEmbedCutoff);
  }
 
  assert(err == NOERROR);

 ErrorHandling: 
  if (newmeth!=NULL) METHOD_DELETE(&newmeth);
  if (localcov!=NULL) COV_DELETE(&localcov);
  return err;
}

int init_circ_embed_cutoff(method_type *meth) {
  return init_circ_embed_local(meth, CircEmbedCutoff);
}
int init_circ_embed_intr(method_type *meth) {
  return init_circ_embed_local(meth, CircEmbedIntrinsic);
}

void do_circ_embed_cutoff(method_type *meth, res_type *res) {  
  localCE_storage *s;
  s = (localCE_storage*) meth->S;
  assert(s->correction == NULL);
  internal_DoSimulateRF(&(s->key), 1, res);
}

void do_circ_embed_intr(method_type *meth, res_type *res) {  
  cov_model *cov = meth->cov;
  location_type *loc = meth->loc;
  double x[MAXCEDIM], dx[MAXCEDIM], *stein_aniso;
  long index[MAXCEDIM], r;
  int i, l, k,
    row = loc->timespacedim,
    col = cov->tsdim,
    rowcol =  row * col;
  localCE_storage *s;

  s = (localCE_storage*) meth->S;
  internal_DoSimulateRF(&(s->key), 1, res);
  
  for (k=0; k<row; k++) {
    index[k] = 0;
    dx[k] = x[k] = 0.0;
  }

  stein_aniso = (double*) s->correction;
  for (i=0; i<rowcol; ) {
    double gauss = GAUSS_RANDOM(1.0);
    for (l=0; l<row; l++)
      dx[l] += stein_aniso[i++] * gauss;
  }
  for (k=0; k<row; k++) dx[k] *= loc->xgr[k][XSTEP];
  for(r=0; ; ) { 
    for (k=0; k<row; k++) 
	res[r] += (res_type) x[k]; 
    r++;
    k=0;
    while( (k<row) && (++index[k]>=loc->length[k])) {
      index[k]=0;
      x[k] = 0.0;
      k++;
    }
    if (k>=row) break;
    x[k] += dx[k];
  }
}
