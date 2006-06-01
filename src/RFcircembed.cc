/*
 Authors 
 Yindeng, Jiang, jiangyindeng@gmail.com
 Martin Schlather, schlath@hsu-hh.de

 Simulation of a random field by circulant embedding
 (see Wood and Chan, or Dietrich and Newsam for the theory )

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2005 Yindeng Jiang & Martin Schlather
 Copyright (C) 2006 -- 2006 Martin Schlather

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
#include "RFsimu.h"
#include <assert.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>

#define MAX_CE_MEM 16777216
ce_param CIRCEMBED={false, true, false, true, TRIVIALSTRATEGY, 3, MAX_CE_MEM, 
		    -1e-7, 1e-3, 0, 0, 0, 0};
ce_param LOCAL_CE={false, true, false, true, TRIVIALSTRATEGY, 1, MAX_CE_MEM, 
		   -1e-9, 1e-7, 0, 0, 0, 0};


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
    if (x->c!=NULL) free(x->c);
    if (x->d!=NULL) free(x->d);
    FFT_destruct(&(x->FFT));
    free(*S);
    *S = NULL;
  }
}

void LOCAL_NULL(localCE_storage* x){
  int i;
  for (i=0; i<MAXDIM; i++) x->correction[i] = NULL;
}

void CE_NULL(CE_storage* x){
  FFT_NULL(&(x->FFT));  
  x->positivedefinite = FALSE;
  x->trials = -1;
  x->c = x->d = NULL;
  x->smallestRe = x->largestAbsIm = NA_REAL;
//  int i;
  //for (i=0; i<MAXDIM; i++) x->[i] = NULL;
}
 
void localCE_destruct(void **S) 
{
  int i;
  if (*S!=NULL) {
    localCE_storage* x;
    x = *((localCE_storage**) S);
    for (i=0; i<MAXDIM; i++) 
      if (x->correction[i] != NULL) free(x->correction[i]);
    LOCAL_NULL(x);
    DeleteKeyNotTrend(&(x->key)); 
    free(*S);
    *S = NULL;
  }
}


/*********************************************************************/
/*           CIRCULANT EMBEDDING METHOD (1994) ALGORITHM             */
/*  (it will always be refered to the paper of Wood & Chan 1994)     */
/*********************************************************************/

void SetParamCircEmbed( int *action, int *force, double *tolRe, double *tolIm,
			int *trials, int *severalrealisations,
			double *mmin, int *useprimes, int *strategy,
		        double *maxmem, int *dependent) 
{
  SetParamCE(action, force, tolRe, tolIm, trials, severalrealisations,
	     mmin, useprimes, strategy,
	     maxmem, dependent, &CIRCEMBED, "CIRCEMBED");
}

void SetParamLocal( int *action, int *force, double *tolRe, double *tolIm,
		     int *severalrealisations, double *mmin, int *useprimes,
		    double *maxmem, int *dependent) 
{
  int ONE=1;
  SetParamCE(action, force, tolRe, tolIm, &ONE, severalrealisations,
	     mmin, useprimes, 
	     &ONE /* anything */, maxmem, dependent, &LOCAL_CE, "LOCAL");
}

int fastfourier(double *data, int *m, int dim, bool first, bool inverse,
		FFT_storage *FFT)
  /* this function is taken from the fft function by Robert Gentleman 
     and Ross Ihaka, in R */
{
    long int inv, nseg, n,nspn,i;
  int maxf, maxp, Xerror;
  if (first) {
   int maxmaxf,maxmaxp;
   nseg = maxmaxf =  maxmaxp = 1;
   /* do whole loop just for Xerror checking and maxmax[fp] .. */

   for (i = 0; i<dim; i++) {
     if (m[i] > 1) {
       fft_factor(m[i], &maxf, &maxp);
       if (maxf == 0) {Xerror=ERRORFOURIER; goto ErrorHandling;}	
       if (maxf > maxmaxf) maxmaxf = maxf;
       if (maxp > maxmaxp) maxmaxp = maxp;
       nseg *= m[i];
     }
   }
 
   assert(FFT->work==NULL);
   if ((FFT->work = (double*) malloc(4 * maxmaxf * sizeof(double)))==NULL) { 
     Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling; 
   }
   assert(FFT->iwork==NULL);
   if ((FFT->iwork = (int*) malloc(maxmaxp  * sizeof(int)))==NULL) { 
     Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
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
  return Xerror;
}

int fastfourier(double *data, int *m, int dim, bool first, FFT_storage *FFT){
  return fastfourier(data, m, dim, first, !first, FFT);
}


int init_circ_embed(key_type *key, int m)
{
  methodvalue_type *meth;  
  int Xerror=NOERROR, d, actcov;
  double steps[MAXDIM], *c;
  CE_storage *s;
  ce_param* cepar;

  c = NULL;
  meth = &(key->meth[m]);
  cepar = &CIRCEMBED;
  s = NULL;
  if (!key->grid) {Xerror=ERRORONLYGRIDALLOWED;goto ErrorHandling;}
  SET_DESTRUCT(CE_destruct, m);

  if ((meth->S = malloc(sizeof(CE_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  } 
  s =  (CE_storage*) meth->S;
  CE_NULL(s);
  

  if ((Xerror = FirstCheck_Cov(key, m, true)) != NOERROR)  goto ErrorHandling;
  actcov = meth->actcov;

  for (d=0; d<key->timespacedim; d++) {
    s->nn[d]=key->length[d]; 
    steps[d]=key->x[d][XSTEP];
  }

  int *mm, *cumm, *halfm, dim;
  double hx[MAXDIM];
  int index[MAXDIM], dummy;
  long mtot,i,k,twoi;
  bool cur_crit;
  
  mtot=-1;
  mm = s->m;
  halfm= s->halfm;
  cumm = s->cumm;
  dim = key->timespacedim;
  c=NULL;
  if (GENERAL_PRINTLEVEL>=5) PRINTF("calculating the Fourier transform\n");

  /* cumm[i+1]=\prod_{j=0}^i m[j] 
     cumm is used for fast transforming the matrix indices into an  
       index of the vector (the way the matrix is stored) corresponding 
       to the matrix                   */                               
  /* calculate the dimensions of the matrix C, eq. (2.2) in W&C */       
                                                                         
 // ("CE missing strategy for matrix entension in case of anisotropic fields %d\n
  for (i=0;i<dim;i++){ // size of matrix at the beginning  
    if (fabs(cepar->mmin[i]) > 10000.0) {       
	sprintf(ERRORSTRING_OK, "maximimal modulus of mmin is 10000");
	sprintf(ERRORSTRING_WRONG,"%f", cepar->mmin[i]);
	return ERRORMSG;
    }
    mm[i] = s->nn[i];
    if (cepar->mmin[i]>0.0) {
      if (mm[i] > (1 + (int) ceil(cepar->mmin[i])) / 2) { // plus 1 since 
	// mmin might be odd; so the next even number should be used
	sprintf(ERRORSTRING_OK, "Minimum size in direction %d is %d", 
		(int) i, mm[i]);
	sprintf(ERRORSTRING_WRONG,"%d", (int) ceil(cepar->mmin[i]));
	return ERRORMSG;
      }
      mm[i] = (1 + (int) ceil(cepar->mmin[i])) / 2;
    } else if (cepar->mmin[i] < 0.0) {
	assert(cepar->mmin[i] <= -1.0);
	mm[i] = (int) ceil((double) mm[i] * -cepar->mmin[i]);
    }
    if (cepar->useprimes) {
      // Note! algorithm fails if mm[i] is not a multiple of 2
      //       2 may not be put into NiceFFTNumber since it does not
      //       guarantee that the result is even even if the input is even !
      mm[i] = 2 * NiceFFTNumber((unsigned long) mm[i]);
    } else {
      mm[i] = (1 << 1 + (int) ceil(log((double) mm[i]) * INVLOG2 - EPSILON1000));
    }
  }                             


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
  while (!s->positivedefinite && (s->trials<cepar->trials)){
    (s->trials)++;
    calculate the covariance values "c" according to the given "m"
    fastfourier(c)
    if (!cepar->force || (s->trials<cepar->trials)) {
      check if positive definite
      if (!s->positivedefinite && (s->trials<cepar->trials)) {
        enlarge "m"
      }
    } else 
      print "forced"
  }
*/

  s->trials=0;
  while (!s->positivedefinite && (s->trials<cepar->trials)){ 
    (s->trials)++;
    cumm[0]=1; 
    for(i=0;i<dim;i++){
      halfm[i]=mm[i]/2; 
      index[i]=1-halfm[i]; 
      cumm[i+1]=cumm[i] * mm[i]; // only relevant up to i<dim-1 !!
    }
    s->mtot = mtot = cumm[dim-1] * mm[dim-1]; 

    if (GENERAL_PRINTLEVEL>=2) {
      for (i=0;i<dim;i++) PRINTF("mm[%d]=%d, ",i,mm[i]);
      PRINTF("mtot=%d\n ",mtot);
    }

    if (mtot > cepar->maxmem) {
	sprintf(ERRORSTRING_OK, "%f", cepar->maxmem);
	sprintf(ERRORSTRING_WRONG,"%f", (double) mtot);
	return ERRORMAXMEMORY;
    }


    if (c != NULL) free(c); 
    // for the following, see the paper by Wood and Chan!
    // meaning of following variable c, see eq. (3.8)
    if ((c = (double*) malloc(2 * mtot * sizeof(double))) == NULL) {
      Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
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

      c[dummy] =      
         key->covFct(hx, dim, key->cov, meth->covlist, actcov, key->anisotropy);
      if (cur_crit) {
	for (k=0; k<dim; k++) {
          if (index[k]==halfm[k]) hx[k] -= steps[k] * (double) mm[k];
	}
	c[dummy] = 
	    0.5 *(c[dummy] + key->covFct(hx, dim, key->cov, meth->covlist,
					 actcov, key->anisotropy));
      }
           
      c[dummy+1] = 0.0;
      k=0; while( (k<dim) && (++(index[k]) >= mm[k])) {
	index[k]=0;
	k++;
      }
      assert( (k<dim) || (i==mtot-1));
    }

    if (GENERAL_PRINTLEVEL>6) PRINTF("FFT...");
    if ((Xerror=fastfourier(c, mm, dim, true, &(s->FFT)))!=0) 
	goto ErrorHandling;  
    if (GENERAL_PRINTLEVEL>6) PRINTF("finished\n");
   
    // check if positive definite. If not: enlarge and restart 
    if (!cepar->force || (s->trials<cepar->trials)) { 
      long int mtot2;
      mtot2 = mtot * 2; 
      twoi=0;
      // 16.9. < cepar.tol.im  changed to <=
      while ((twoi<mtot2) && (s->positivedefinite=(c[twoi]>=cepar->tol_re) && 
					  (fabs(c[twoi+1])<=cepar->tol_im)))
	{twoi+=2;}      

      if ( !s->positivedefinite) {
        if (GENERAL_PRINTLEVEL>=2)
	  // 1.1.71: %f changed to %e because c[twoi+1] is usually very small
	  PRINTF("non-positive eigenvalue: c[%d])=%e + %e i.\n",
		 (int) (twoi/2), c[twoi], c[twoi+1]);
	if (GENERAL_PRINTLEVEL>=4) { // just for printing the smallest 
	  //                            eigenvalue (min(c))
	  long int index=twoi, sum=0;
	  double smallest=c[twoi];
	  char percent[]="%";
	  for (twoi=0; twoi<mtot2; twoi += 2) {
	    if (c[twoi] < 0) {
	      sum++;
	      if (c[twoi]<smallest) { index=twoi; smallest=c[index]; }
	    } 
	  }
	  PRINTF("   There are %d negative eigenvalues (%5.3f%s).\n", 
		 sum, 100.0 * (double) sum / (double) mtot, percent);
	  PRINTF("   The smallest eigenvalue is c[%d] =  %e\n",
		 index / 2, smallest);
	}
      }

      if (!s->positivedefinite && (s->trials<cepar->trials)) { 
	FFT_destruct(&(s->FFT));
	switch (cepar->strategy) {
	case 0 :
	  for (i=0; i<dim; i++) { /* enlarge uniformly in each direction, maybe 
				   this should be modified if the grid has 
				   different sizes in different directions */
	    mm[i] <<= 1;
	  }
	  break;
	case 1 :  
	  double cc, maxcc, hx[MAXDIM];
	  int maxi;
	  maxi = -1;
	  maxcc = RF_NEGINF;
	  for (i=0; i<dim; i++) hx[i] = 0.0;
	  for (i=0; i<dim; i++) {
	    hx[i] = steps[i] * mm[i]; 
	    cc = fabs(key->covFct(hx, dim, key->cov, meth->covlist, actcov,
			     key->anisotropy)); 
	    if (GENERAL_PRINTLEVEL>2) PRINTF("%d cc=%e (%e)",i,cc,hx[i]);
	    if (cc>maxcc) {
	      maxcc = cc;
	      maxi = i; 
	    }
	    hx[i] = 0.0;
	  }
	  assert(maxi>=0);
	  mm[maxi] <<= 1;
	  break;
	default:
	  assert(false);
	}
      }
    } else {if (GENERAL_PRINTLEVEL>=2) PRINTF("forced\n");}
  } // while (!s->positivedefinite && (s->trials<cepar->trials))
  assert(mtot>0);
  if (s->positivedefinite || cepar->force) { 
    // correct theoretically impossible values, that are still within 
    // tolerance CIRCEMBED.tol_re/CIRCEMBED.tol_im 
    double r, imag;
    r = c[0];
    imag = 0.0;    
    for(i=0,twoi=0;i<mtot;i++) {
      if (c[twoi] > 0.0) {
	c[twoi] = sqrt(c[twoi]);
      } else {
	if (c[twoi] < r) r = c[twoi];
	c[twoi] = 0.0;
      }
      {
	register double a;
	if ((a=fabs(c[twoi+1])) > imag) imag = a;
      }
      c[twoi+1] = 0.0;
      twoi+=2;
    }
    if (GENERAL_PRINTLEVEL>1) {
      if (r < -GENERAL_PRECISION || imag > GENERAL_PRECISION) {
	PRINTF("using approximating circulant embedding:\n");
	if (r < -GENERAL_PRECISION) 
	  PRINTF("\tsmallest real part has been %e \n", r);
	if (imag > GENERAL_PRECISION) 
	  PRINTF("\tlargest modulus of the imaginary part has been %e \n",imag);
      }
    }
    s->smallestRe = (r > 0.0) ? NA_REAL : r;
    s->largestAbsIm = imag;
  } else {
    Xerror = ERRORFAILED;
    goto ErrorHandling;
  }
  if (GENERAL_PRINTLEVEL >= 10) {
    for (i=0; i < 2 * mtot; i++) {PRINTF("%f ",c[i]);} PRINTF("\n");
  }  
  
  s->dependent = cepar->dependent;
  s->new_simulation_next = true;
  for(i=0; i<dim; i++) { // set anyway -- does not cost too much
    s->cur_square[i] = 0;
    s->max_squares[i] = mm[i] / s->nn[i];
    s-> square_seg[i] = cumm[i] * (s->nn[i] + (mm[i] - s->max_squares[i] * 
					       s->nn[i]) / s->max_squares[i]);
  }
  
  if (cepar->severalrealisations) {
    if ((s->d=(double *) calloc(2 * mtot, sizeof(double)))==0) {
      Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;} //d
  }
  s->c = c;

//  printf("xxxx %d", (int) s);

  return NOERROR;
  
 ErrorHandling:
  if (GENERAL_STORING && s != NULL) s->c = c;
  else if (c!=NULL) {free(c);}
  return Xerror;
}


void covcpy(covinfo_type *dest, covinfo_type *source){
  assert(source->x==NULL);
  memcpy(dest, source, sizeof(covinfo_type));
}

double GetScaledDiameter(key_type *key, covinfo_type *kc) {
// SCALE and ANISO is considered as space trafo and envolved here
  double diameter; 
  diameter = 0.0;
  if (key->anisotropy) { // code can be reduced to the anisotropic case
      //      with the obvious loss of time in the isotropic case...
    int i, j, k, dim, ncorner_dim;
    double distsq;
    double sx[ZWEIHOCHMAXDIM * MAXDIM];
    dim = kc->reduceddim;
    ncorner_dim = (1 << key->timespacedim) * dim;
    GetCornersOfGrid(key, dim, kc->aniso, sx);
    for (i=dim; i<ncorner_dim; i+=dim) {
      for (j=0; j<i; j+=dim) {
        distsq = 0.0;
	for (k=0; k<dim; k++) { 
	  register double dummy;
	  dummy = fabs(sx[i + k] - sx[j + k]);
//	  printf("getscale %d %d %d %f\n", i, j, k, dummy);
	  distsq += dummy * dummy;
	}
	if (distsq > diameter) diameter = distsq;
      }
    }
  } else { // see above
    int d;
    for (d=0; d<key->timespacedim; d++) {
      double dummy;
      dummy = key->x[d][XSTEP] * (double) (key->length[d] - 1) * kc->aniso[0];
      diameter += dummy * dummy; 
    }
  }
//  printf("diameter=%f", sqrt(diameter));
  return sqrt(diameter);
}


int GetOrthogonalUnitExtensions(aniso_type aniso, int dim, double *grid_ext) {
  int k,i,j,l,m, job=01, err, dimsq, ev0, jump, endfor;
  double s[MAXDIMSQ], G[MAXDIM+1], e[MAXDIM], D[MAXDIM], V[MAXDIMSQ];
  dimsq = dim * dim;
  for (k=0; k<dim; k++) {
    for (i=0; i<dim; i++) {
      for (l=j=0; j<dimsq; j+=dim) {
	s[j + i] = 0.0;
	jump = l + k;
	endfor = l + dim;
	for (m=i; l<endfor; l++, m += dim) {
	  if (l!=jump) s[i + j] += aniso[l] * aniso[m];
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
      grid_ext[k] += V[ev0 + i] * aniso[k + i * dim];
    }
    grid_ext[k] = fabs(grid_ext[k]);
  }
  return NOERROR;

 ErrorHandling:
  if (err<0) {
    PRINTF("F77 error in GetOrthogonalExtensions: %d\n", -err);
    err = ERRORFAILED;
  }
  return err;
}


int init_circ_embed_local(key_type *key, int m)
{
  methodvalue_type *meth;  
  int Xerror=NOERROR, v, store_msg[MAXDIM], msg, simuactcov, instance,
      actcov, nc, store_nr, i, dimsq,  n_aniso, err;
  localCE_storage *s;
  cov_fct *hyper;
  covinfo_type *sc;
  bool selectlocal[Forbidden + 1];
  param_type store_param;
  double rawRmax[MAXDIM], grid_ext[MAXDIM];

  if (key->covFct != CovFct) { Xerror=ERRORNOTPROGRAMMED; goto ErrorHandling;}
  SET_DESTRUCT(localCE_destruct, m);
  if (!key->grid) { Xerror=ERRORONLYGRIDALLOWED; goto ErrorHandling;}
  meth = &(key->meth[m]);
  assert(meth->S==NULL);
   if ((meth->S=malloc(sizeof(localCE_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  } 
  s = (localCE_storage*) meth->S;
  LOCAL_NULL(s);
  KEY_NULL(&(s->key));

  // prepared for more sophisticated selections if init_circ_embed_local
  // and init_circ_embed are merged
  for (v=0; v < (int) Forbidden; v++) selectlocal[v]=false;
  selectlocal[meth->unimeth] = true;

  if ((Xerror = FirstCheck_Cov(key, m, false)) != NOERROR) goto ErrorHandling;
  actcov = meth->actcov;
  if (2 * actcov > MAXCOV) {
    Xerror=ERRORNCOVOUTOFRANGE; goto ErrorHandling;
  }

  n_aniso = key->anisotropy ? key->timespacedim * key->timespacedim : 1;

  for (i=0; i<key->timespacedim; i++) rawRmax[i] = 0.0;
  for (simuactcov=v=0; v<actcov; v++) {
    covinfo_type *kc;
    store_nr = -1;
    kc = &(key->cov[meth->covlist[v]]);
//    printf("** %d\n", kc->nr);
    sc = &(s->key.cov[simuactcov]);
    sc->method = CircEmbed;
    sc->dim = key->timespacedim;
    sc->param[VARIANCE] = 1.0; 
    sc->param[HYPERNR] = 1; // as only addition between models allowed
    // SCALE is considered as space trafo and envolved here
    sc->param[DIAMETER] = GetScaledDiameter(key, kc);
    if (GENERAL_PRINTLEVEL>7) PRINTF("diameter %f\n", sc->param[DIAMETER]);
    memcpy(&(sc->param[ANISO]), &(kc->param[ANISO]), sizeof(double) * n_aniso);
    sc->op = 2;
    covcpy(sc + 1, kc);

    instance = -1;
    store_param[LOCAL_R] =  R_PosInf;
    store_msg[simuactcov] = MSGLOCAL_JUSTTRY; // lowest preference
    for (nc=0; nc < nlocal; nc++) {// in case more cutoff methods are found
      instance = 0;
      sc->nr = LocalCovList[nc];
      hyper = &(CovList[sc->nr]);
      if (GENERAL_PRINTLEVEL>7) PRINTF("%d %s:\n", nc, hyper->name);
      if (!selectlocal[hyper->localtype]) continue;
      // for future use if more hyper functions for cutoff and intrinsic
     // are available
      // if (PREFERENCE[hyper->localtype] >= 0 && 
      //	PREFERENCE[hyper->localtype]!=sc->nr) continue; 
      while ((msg = hyper->getparam(sc + 1, sc->param, instance++))
	     != MSGLOCAL_ENDOFLIST) {
//	printf("%d %s %d \n", instance, hyper->name, msg);

	if ((Xerror=
	     hyper->checkNinit(sc, &(COVLISTALL[v]), 
			       key->ncov - meth->covlist[v], 
			       CircEmbed, key->anisotropy)) == NOERROR
	    && (
		(int) msg < (int) store_msg[simuactcov] 
		||
		((int) msg == (int) store_msg[simuactcov] &&
		 sc->param[LOCAL_R] < store_param[LOCAL_R])
		)
	    )
	{	   
	    memcpy(store_param, sc->param, sizeof(param_type));
	    store_msg[simuactcov] = msg;
	    store_nr = sc->nr;
	} // if

//	printf("%d %d %f \n",  Xerror, store_msg[simuactcov], store_param[LOCAL_R] );

        if (GENERAL_PRINTLEVEL>7) {
	    PRINTF("v=%d nc=%d inst=%d err=%d hyp.kappa=%f, #=%d, diam=%f\n r=%f curmin_r=%f\n", 
		   v, nc, instance, Xerror, sc->param[HYPERKAPPAII],
		   (int) sc->param[HYPERNR],
		   sc->param[DIAMETER],  sc->param[LOCAL_R],
		   store_param[LOCAL_R]);
	    if (hyper->implemented[CircEmbedIntrinsic] == HYPERIMPLEMENTED)
		PRINTF("userr=%f, a0=%f, a2=%f, b=%f\n",
		       sc->param[INTRINSIC_RAWR],
		       sc->param[INTRINSIC_A0], sc->param[INTRINSIC_A2],
		       sc->param[INTRINSIC_B]
		    );
	    else {
//	      printf("%s %d\n", 
//		     hyper->name, hyper->implemented[CircEmbedCutoff]);
	      assert(hyper->implemented[CircEmbedCutoff] == HYPERIMPLEMENTED);
	      PRINTF("a=%f, a_sqrt.r=%f, b=%f\n",
		     sc->param[CUTOFF_A],
		     sc->param[CUTOFF_ASQRTR], sc->param[CUTOFF_B]
		);
	    }
	}
      } // while
      err =  GetOrthogonalUnitExtensions(kc->aniso, key->timespacedim, grid_ext);
      if (err != NOERROR) {
	  // err is extra since Xerror at this stage can be both NOERROR and
	  // and an error message; the latter could be overwritten by
	  // NOERROR; Xerror is only of relavance if all loops in while
	  // fail; so it cannot be changed; see also assert(Xerror!=NOERROR)
	  // below where consistency is checked.
	  Xerror=err; 
	  goto ErrorHandling;
      }

      for (i=0; i<key->timespacedim; i++) {
        register double dummy;   
	dummy = store_param[LOCAL_R] / 
	    (grid_ext[i] * (double) (key->length[i] - 1) * key->x[i][XSTEP]);
//  	printf("\nXXXX %d grid_ext=%f dummy=%f, raw=%f local_r=%f %d %f\n",
//  	       i, grid_ext[i], dummy,  rawRmax[i], store_param[LOCAL_R],
//	       key->length[i], key->x[i][XSTEP]);
	if (rawRmax[i] < dummy) rawRmax[i] = dummy;
      }
    } // nc
    if (!R_FINITE(store_param[LOCAL_R])) { 
      if (GENERAL_PRINTLEVEL>3) {
	PRINTF("v=%d nc=%d, inst=%d err=%d hyp.kappa=%f, #=%d, diam=%f\nr=%f curmin_r=%f\n", 
	       v, nc, instance, Xerror, sc->param[HYPERKAPPAII],
	       (int) sc->param[HYPERNR],
	     sc->param[DIAMETER],  sc->param[LOCAL_R],
	       store_param[LOCAL_R]);
      }
      assert(Xerror!=NOERROR);
      goto ErrorHandling;
    }
    sc->nr = store_nr;
    assert(sc->nr >= 0);
    memcpy(sc->param, store_param, sizeof(param_type));
    simuactcov += (short int) sc->param[HYPERNR] + 1; 
    sc += (int) sc->param[HYPERNR];  
  } // v

  // prepare for call of internal_InitSimulateRF
  int covnr[MAXCOV], op[MAXCOV], CEMethod[MAXCOV], cum_nParam;
  double ParamList[MAXCOV * TOTAL_PARAM];
  char errorloc_save[nErrorLoc];
  cum_nParam = 0;
  for (v=0; v<simuactcov; v++) {
    int kappas;
    sc = &(s->key.cov[v]);
    covnr[v] = sc->nr;
    op[v] = sc->op;
    CEMethod[v] = (int) CircEmbed;
    ParamList[cum_nParam++] = sc->param[VARIANCE];
    kappas = CovList[sc->nr].kappas(key->timespacedim);
    for (i=0; i<kappas; i++)
	ParamList[cum_nParam++] = sc->param[KAPPA + i];
    for (i=0; i<n_aniso; i++)
	ParamList[cum_nParam++] = sc->param[ANISO + i];
  }

  strcpy(errorloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "%s%s ", errorloc_save, "local circ. embed.: ");
  ce_param ce_save;
  memcpy(&ce_save, &CIRCEMBED, sizeof(ce_param));
  memcpy(&CIRCEMBED, &LOCAL_CE, sizeof(ce_param));
 
  for (i=0; i<key->timespacedim; i++) {
    if (CIRCEMBED.mmin[i]==0.0) CIRCEMBED.mmin[i] = - rawRmax[i];
//    printf("%d %f\n", i, CIRCEMBED.mmin[i]);
  }

  instance = 0;
  for (;;) {
    bool anychangings;
    Xerror = internal_InitSimulateRF(key->x[0], key->T, key->spatialdim,
				     3, key->grid, key->Time,
				     covnr, ParamList, cum_nParam,
				     simuactcov, key->anisotropy, op,
				     CEMethod, DISTR_GAUSS, &(s->key), 
				     0 /* natural scaling */, CovFct);
    anychangings = false;
    if (Xerror == NOERROR || Xerror == USEOLDINIT) {
      break;
    } else {
      for (v=0; v<simuactcov; v++, instance++) {
	sc = &(s->key.cov[v]);
	hyper = &(CovList[sc->nr]);
	if ((hyper->localtype == CircEmbedCutoff || 
	     hyper->localtype == CircEmbedIntrinsic)) {
	  if (store_msg[v] != MSGLOCAL_OK) {
	    bool improved;
	    memcpy(store_param, sc->param, sizeof(param_type));
	    improved = hyper->alternative(sc, instance) && 
		(hyper->checkNinit(sc, &(COVLISTALL[v]), 
				   (int) sc->param[HYPERNR], 
				   CircEmbed, key->anisotropy) != NOERROR);
	    if (improved) memcpy(sc->param, store_param, sizeof(param_type));
	    anychangings |=  improved;
	  }
	} else {
	  assert((v % 2) == 1); // if the method is generalised to
	  // multiplicative submodels this assert reminds that
	  // this part must be changed.
	}
      }
      if (!anychangings) break;
    }
  }

  memcpy(&CIRCEMBED, &ce_save, sizeof(ce_param));
  strcpy(ERROR_LOC, errorloc_save);
  if (Xerror != NOERROR) goto ErrorHandling;

  dimsq = key->timespacedim * key->timespacedim;
  for (v=0; v<simuactcov; v++) {
    double dummy, *stein_aniso;
    sc = &(s->key.cov[v]);
    hyper = &(CovList[sc->nr]);
    if (hyper->localtype == CircEmbedIntrinsic) {
      assert((v % 2) == 0); // if the method is generalised to
	  // multiplicative submodels this assert reminds that
	  // this part must be changed.
      dummy = sqrt(2.0 * sc->param[INTRINSIC_A2]); // see Stein (2002)
      sc = &(s->key.cov[v]);
      if ((s->correction[v] = malloc(sizeof(double) * dimsq))==NULL){
        Xerror=ERRORMEMORYALLOCATION;
	goto ErrorHandling;
      }
      stein_aniso = (double*) s->correction[v];
      if (key->anisotropy) { // distinction necessary,
	  // and sc->aniso may not be used, but the 
	  // orginial param[ANISO] -- it did not work otherwise
	for (i=0; i<dimsq; i++) {
	  stein_aniso[i] = dummy * sc->param[ANISO + i];
//    printf("aniso %d %f %f %f\n", i, stein_aniso[i], dummy, 
//	   sc->param[ANISO + i]);
	}
      } else {
	  int i;
	for (i=0; i<dimsq; stein_aniso[i++] = 0.0);
	for (i=0; i<dimsq; i += key->timespacedim + 1) 
          stein_aniso[i] = dummy / sc->param[SCALE];
      }
      v++;
    } else {
      assert(hyper->localtype == CircEmbedCutoff);
      v++;
    }
  }
  return NOERROR;
 
 ErrorHandling:
  return Xerror;
}


void do_circ_embed_local(key_type *key, int m, double *res )
{  
  double x[MAXDIM], dx[MAXDIM], *stein_aniso;
  long index[MAXDIM], r;
  int v, ncov, dim, dimsq, i, l, k;
  localCE_storage *s;
  bool stein_correction;
  covinfo_type *sc;
  cov_fct *hyper;

  s = (localCE_storage*) key->meth[m].S;
  internal_DoSimulateRF(&(s->key), 1, res);
  
  dim = key->timespacedim;
  dimsq = dim * dim;
  ncov  = s->key.ncov;
  for (k=0; k<key->timespacedim; k++) {
    index[k] = 0;
    dx[k] = x[k] = 0.0;
  }
  stein_correction = false;
  for (v=0; v<ncov; v++) {
    sc = &(s->key.cov[v]);
    hyper = &(CovList[sc->nr]);
    switch(hyper->localtype) {
	case CircEmbedCutoff : 	    
	    assert(s->correction[v++] == NULL);
	    break;
	case CircEmbedIntrinsic:
	    stein_correction = true;
	    stein_aniso = (double*) s->correction[v++];
	    for (i=0; i<dimsq; ) {
	      double gauss = GAUSS_RANDOM(1.0);
	      for (l=0; l<dim; l++)
	        dx[l] += stein_aniso[i++] * gauss;
	    }
	    break;
	default: assert(false);
    }
  } // v
  if (stein_correction) {
    for (k=0; k<dim; k++) dx[k] *= key->x[k][XSTEP];
    for(r=0;;) { 
      for (k=0; k<dim; k++) res[r] += x[k]; 
      r++;
      k=0;
      while( (k<dim) && (++index[k]>=key->length[k])) {
	index[k]=0;
	x[k] = 0.0;
	k++;
      }
      if (k>=dim) break;
      x[k] += dx[k];
    }
  }
}

void do_circ_embed(key_type *key, int m, double *res){
  int  i, j, k, HalfMp1[MAXDIM], HalfMaM[2][MAXDIM], index[MAXDIM], dim,
    *mm, *cumm, *halfm;
  double XX,YY,invsqrtmtot, *c, *d;
  bool first, free[MAXDIM+1], noexception;
  long mtot, start[MAXDIM], end[MAXDIM];
  CE_storage *s;

//  printf("%d %d %d\n", m, (int) key->meth[m].S, key->n_unimeth);
// printf("key->n_unimeth %d\n", key->n_unimeth);

  s = (CE_storage*)key->meth[m].S;
  if (s->d==NULL) { /* overwrite the intermediate result directly
		       (algorithm allows for that) */
    d=s->c;
  }
  else {
    d=s->d;
  }

/* 
     implemented here only for rotationsinvariant covariance functions
     for arbitrary dimensions;
     (so it works only for even covariance functions in the sense of 
     Wood and Chan,p. 415, although they have suggested a more general 
     algorithm;) 
     
*/  
  
  dim = key->timespacedim;
  mm = s->m;
  cumm = s->cumm;
  halfm = s->halfm;
  c = s->c;
  mtot= s->mtot; 
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

  if (GENERAL_PRINTLEVEL>=10) PRINTF("Creating Gaussian variables... \n");
  /* now the Gaussian r.v. have to defined and multiplied with sqrt(FFT(c))*/

  for (i=0; i<dim; i++) {
    index[i]=0; 
    free[i] = false;
  }
  free[MAXDIM] = false;

  if (s->new_simulation_next) {
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
//      if (GENERAL_PRINTLEVEL>=10) PRINTF("cumm...");
      i <<= 1; // since we have to index imaginary numbers
      j <<= 1;
      if (noexception) { // case 3 in prop 3 of W&C
	XX = GAUSS_RANDOM(INVSQRTTWO);
	YY = GAUSS_RANDOM(INVSQRTTWO);
	d[i] = d[i+1] = c[i];   d[i] *= XX;  d[i+1] *= YY;   
	d[j] = d[j+1] = c[j];   d[j] *= XX;  d[j+1] *= -YY; 
      } else { // case 2 in prop 3 of W&C
	d[i]   = c[i] * GAUSS_RANDOM(1.0);
	d[i+1] = 0.0;
      }
      
//    if (i==224|| i==226 || i==228) {
//printf("B %d %d %d %d %d %f %f\n", 
//	   index[0], index[1], i, j, noexception, d[i], d[i+1]);
// assert(false); 
//    }
      
//      if (GENERAL_PRINTLEVEL>=10) PRINTF("k=%d ", k);
      
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
	 
	 free==true   <=>   endfor==mm[]	 
      */
      k=0; 
      if (++index[k]>HalfMaM[free[k]][k]) {
	// in case k increases the number of indices that run over 0..m increases
	free[k] = true;
	index[k]= 0; 
	k++;
	while((k<dim) && (++index[k]>HalfMaM[free[k]][k])) {
	  free[k] = true;
	  index[k]= 0; 
	  k++;
	}
	if (k>=dim) break;
	// except the very last (new) number is halfm and the next index is
	// restricted to 0..halfm
	// then k decreases as long as the index[k] is 0 or halfm
	if (!free[k] && (index[k]==halfm[k])){//index restricted to 0..halfm?
	  // first: index[k] is halfm? (test on ==0 is superfluent)
	  k--;
	  while ( (k>=0) && ((index[k]==0) || (index[k]==halfm[k]))) {
	    // second and following: index[k] is 0 or halfm?
	    free[k] = false;
	  k--;
	  }
	}
      }
  }
    

//  double zz=0.0;
//  for (i =0; i<2 * mtot; i++) {
//    if (true || i<=100) printf("%d %d\n", i, 2 * mtot);
//    zz += d[i];
//     assert(i!=230);
//     // if (i>100); break;
//  }
// printf("%f\n", zz);

    fastfourier(d, mm, dim, false, &(s->FFT));   
  } // if simulate

  /* now we correct the result of the fastfourier transformation
     by the factor 1/sqrt(mtot) and read the relevant matrix out of 
     the large vector c */
  first = true;
  for(i=0; i<dim; i++) {
    index[i] = start[i] = s->cur_square[i] * s->square_seg[i];
    end[i] = start[i] + cumm[i] * s->nn[i];
//    printf("%d start=%d end=%d cur=%d max=%d seq=%d cumm=%d, nn=%d\n",
//	   i, (int) start[i], (int) end[i],
//	   s->cur_square[i], s->max_squares[i],
//	   (int) s->square_seg[i], cumm[i], s->nn[i]);
  }
  int totpts = key->totalpoints;
  for (i=0; i<totpts; i++){
    j=0; for (k=0; k<dim; k++) {j+=index[k];} 
    res[i] += d[2 * j] * invsqrtmtot; 
    if ((fabs(d[2*j+1])>CIRCEMBED.tol_im) && 
        ((GENERAL_PRINTLEVEL>=2 && first) || GENERAL_PRINTLEVEL>=6)){ 
      PRINTF("IMAGINARY PART <> 0, %e\n", d[2 * j + 1]);
      first=false; 
    } 
    for(k=0; (k<dim) && ((index[k] += cumm[k]) >= end[k]); k++) {
      index[k]=start[k]; 
    }
  }

  if (s->dependent) {
//  if MaxStableRF calls Cutoff, c and look uninitialised
    k=0; 
    while(k<dim && (++(s->cur_square[k]) >= s->max_squares[k])) {
      s->cur_square[k++]=0;
    }
    s->new_simulation_next = k==dim; 
  }
  //printf("%d %d\n", s->cur_square[0],s->cur_square[1]);
  key->stop |= s->new_simulation_next && s->d == NULL;
}
