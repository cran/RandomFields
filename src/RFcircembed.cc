/*
 Authors 
 Martin Schlather, martin.schlather@cu.lu 

 Simulation of a random field by circulant embedding
 (see Wood and Chan, or Dietrich and Newsam for the theory )

 Copyright (C) 2001 -- 2004 Martin Schlather, 

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

#ifdef doubleReal
ce_param CIRCEMBED={false, false, TRIVIALSTARTEGY, -1e-7, 1e-3, 3, 0, 0, 0, 0}
#else
ce_param CIRCEMBED={false, false, TRIVIALSTARTEGY, -1e-5, 1e-3, 3, 0, 0, 0, 0};
#endif

typedef struct CE_storage {
  int m[MAXDIM],halfm[MAXDIM],nn[MAXDIM],cumm[MAXDIM+1]; /* !!!! **** */  
  param_type param; /* only used in local CE */
  double *c,*d;
  Real *local, sqrt_var_squarefactor[MAXCOV]; /* only used in local CE */
  FFT_storage FFT;
  long totalpoints;
} CE_storage;

void FFT_destruct(FFT_storage *FFT) 
{
  if (FFT->iwork!=NULL) {free(FFT->iwork); FFT->iwork=NULL;}
  if (FFT->work!=NULL) {free(FFT->work); FFT->work=NULL;} //?
}

void FFT_NULL(FFT_storage *FFT) 
{
  FFT->work = NULL;
  FFT->iwork = NULL;
}

void CE_destruct(void **S) 
{
  if (*S!=NULL) {
    CE_storage *x;
    x = *((CE_storage**)S);
    if (x->c!=NULL) free(x->c);
    if (x->d!=NULL) free(x->d);
    if (x->local!=NULL) free(x->local);
    FFT_destruct(&(x->FFT));
    free(*S);
    *S = NULL;
  }
}

/*********************************************************************/
/*           CIRCULANT EMBEDDING METHOD (1994) ALGORITHM             */
/*  (it will always be refered to the paper of Wood & Chan 1994)     */
/*********************************************************************/

void SetParamCircEmbed( int *action, int *force, Real *tolRe, Real *tolIm,
			int *trials, int *mmin, int *userfft, int *strategy) 
{
  SetParamCE(action, force, tolRe, tolIm, trials, mmin, userfft, strategy,
	     &CIRCEMBED, "CIRCEMBED");
}




int fastfourier(double *data, int *m, int dim, bool first, bool inverse,
		FFT_storage *FFT)
  /* this function is taken from the fft function by Robert Gentleman 
     and Ross Ihaka, in R */
{
  int inv, nseg, n,nspn,i,maxf,maxp,Xerror;
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
 
   if ((FFT->work = (double*) malloc(4 * maxmaxf * sizeof(double)))==NULL) { 
     Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling; 
   }
   if ((FFT->iwork = (int*) malloc( maxmaxp  * sizeof(int)))==NULL) { 
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
  return 0;
 ErrorHandling:
  FFT_destruct(FFT);
  return Xerror;
}

int fastfourier(double *data, int *m, int dim, bool first, FFT_storage *FFT){
  fastfourier(data, m, dim, first, !first, FFT);
}


int internal_init_circ_embed(Real *steps, bool anisotropy,
			     int *covnr, int *op, param_type param,
			     int *nn, int *m, int *cumm, int *halfm, 
			     int dim, int actcov,
			     CovFctType CovFct,
			     ce_param* cepar,
			     FFT_storage *FFT,
			     long *twoRealmtot,
			     double **cc)
{
  double *c;
  Real hx[MAXDIM], totalm;
  int H[MAXDIM];
  int  Xerror,trials,index[MAXDIM],dummy;
  long maxm, MAXM=10000000;  /* the maximum grid size: 
					      MAXM[dim-1]^dim */
  long mtot=-1,i,j,k,twoi;
  bool positivedefinite, cur_crit, critical, Critical[MAXDIM];

  c=NULL;
  if (cepar->userfft) maxm=1000000000;  else maxm=MAXM;
  if (GENERAL_PRINTLEVEL>=5) PRINTF("calculating the Fourier transform\n");
  
  /* cumm[i+1]=\prod_{j=0}^i m[j] 
     cumm is used for fast transforming the matrix indices into an
       index of the vector (the way the matrix is stored) corresponding
       to the matrix                   */
  /* calculate the dimensions of the matrix C, eq. (2.2) in W&C */
  
  // ("CE missing strategy for matrix entension in case of anisotropic fields ! %d \n", cepar->userfft); Namely if in case the field is not quadratic or the covariance function does not have the same range in all directions!

  totalm = 1.0;  
  for (i=0;i<dim;i++){ // size of matrix at the beginning
    int factor;
    if (cepar->userfft) {
      factor = (cepar->mmin[i] > -2) ? 2 : -cepar->mmin[i];
      m[i] =  factor * NiceFFTNumber((unsigned long) nn[i]); 
    } else {
      factor = (cepar->mmin[i] > -1) ? 1 : -cepar->mmin[i];
      m[i]= factor *
	(1 << (1 + (int) ceil(log((Real) nn[i]) * INVLOG2 - EPSILON1000))); 
    }
    if (m[i]<cepar->mmin[i]) {m[i]=cepar->mmin[i];} 
    totalm *= (Real) m[i];
  }
  if (totalm>maxm) {
    if (GENERAL_PRINTLEVEL>=5)
      PRINTF("totalm=%f > maxm=%d userfft=%d\n",totalm,maxm,cepar->userfft);
    Xerror=ERRORFAILED;  goto ErrorHandling;
  } 
    
  positivedefinite = false;     
    /* Eq. (3.12) shows that only j\in I(m) [cf. (3.2)] is needed,
       so only the first two rows of (3.9) (without the taking the
       modulus of h in the first row)
       The following variable `index' corresponds to h(l) in the following
       way: index[l]=h[l]        if 0<=h[l]<=m[l]/2
       index[l]=h[l]-m[l]   if m[l]/2+1<=h[l]<=m[l]-1     
       Then h[l]=(index[l]+m[l]) mod m[l] !!
    */
  
  trials=0;
  while (!positivedefinite && (trials<cepar->trials)){ 
    trials++;
    cumm[0]=1; 
    for(i=0;i<dim;i++){
      halfm[i]=m[i]/2; 
      index[i]=1-halfm[i]; 
      cumm[i+1]=cumm[i] * m[i]; // only relevant up to i<dim-1 !!
    } 

    mtot=cumm[dim-1] * m[dim-1]; 
    if (GENERAL_PRINTLEVEL>=2) {
      for (i=0;i<dim;i++) PRINTF("m[%d]=%d, ",i,m[i]);
      PRINTF("mtot=%d\n ",mtot);
    }
    *twoRealmtot = 2 * sizeof(double) * mtot;

    // for the following, see the paper by Wood and Chan!
    // meaning of following variable c, see eq. (3.8)
    if ((c = (double*) malloc(*twoRealmtot)) == 0) {
      Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }

    // first determine whether there are odd covariance functions
    critical = false;
    for (k=0; k<dim; k++) Critical[k]=false;
    for (i=0; i<actcov; i++) {
      if (CovList[covnr[i]].isotropic==ANISOTROPIC) {
	if (!CovList[covnr[i]].even) {
	  for(k=0; k<dim; k++) {
	    Critical[k] |= CovList[covnr[i]].odd[k];
	  }
	}
      }
    }
    // ***  critical odd not used yet
    
    for(i=0; i<dim; i++){index[i]=0;}
    
    for (i=0; i<mtot; i++){
      cur_crit = false;
      for (k=0; k<dim; k++) {
	hx[k] = steps[k] * 
	  (Real) ((index[k] <= halfm[k]) ? index[k] : index[k] - m[k]);
	cur_crit |= (index[k]==halfm[k]);
      }	
      dummy=i << 1;

      c[dummy] = (critical && cur_crit) ?  0.0 :
	CovFct(hx, dim, covnr, op, param, actcov, anisotropy);

      c[dummy+1] = 0.0;
      k=0; while( (k<dim) && (++(index[k]) >= m[k])) {
	index[k]=0;
	k++;
      }
      assert( (k<dim) || (i==mtot-1));
    }
            
    if (GENERAL_PRINTLEVEL>6) PRINTF("FFT...");
    if ((Xerror=fastfourier(c, m, dim, true, FFT))!=0) goto ErrorHandling;  
    if (GENERAL_PRINTLEVEL>6) PRINTF("finished\n");
   
    // check if positive definite. If not: enlarge and restart 
    if (!cepar->force || (trials<cepar->trials)) {
      i=0; twoi=0;
      // 16.9. < cepar.tol.im  changed to <=
      while ((i<mtot) && (positivedefinite=(c[twoi]>=cepar->tol_re) && 
					  (fabs(c[twoi+1])<=cepar->tol_im)))
	{i++; twoi+=2;}
      if (!positivedefinite) {
	if (GENERAL_PRINTLEVEL>=2)
	  PRINTF(" nonpos %d  %f %f \n",i,c[twoi],c[twoi+1]);
	FFT_destruct(FFT);
	free(c); c=NULL;

	totalm = 1.0;
	switch (cepar->strategy) {
	case 0 :
	  for (i=0;i<dim;i++) { /* enlarge uniformly in each direction, maybe 
				   this should be modified if the grid has 
				   different sizes in different directions */
	    m[i] <<= 1;
	    totalm *= (Real) m[i];
	  }
	  break;
	case 1 :  
	  Real cc, maxcc, hx[MAXDIM];
	  int maxi;
	  maxcc = RF_NEGINF;
	  for (i=0; i<dim; i++) hx[i] = 0.0;
	  for (i=0; i<dim; i++) {
	    hx[i] = steps[i] * m[i]; 
	    cc = fabs(CovFct(hx , dim, covnr, op, param, actcov, anisotropy)); 
	    if (GENERAL_PRINTLEVEL>2) PRINTF("%d cc=%e (%e)",i,cc,hx[i]);
	    if (cc>maxcc) {
	      maxcc = cc;
	      maxi = i; 
	    }
	    hx[i] = 0.0;
	  }	
	  m[maxi] <<= 1;
	  for (i=0;i<dim; i++) totalm *= (Real) m[i];
	  break;
	default:
	  assert(false);
	}
	if (totalm>maxm) {Xerror=ERRORFAILED;goto ErrorHandling;}
	//	assert(false);
      }
    } else {if (GENERAL_PRINTLEVEL>=2) PRINTF("forced\n");}
  }
  assert(mtot>0);
 
  if (GENERAL_PRINTLEVEL>4) {  // delete this part later on
    Real cc, maxcc, hx[MAXDIM];
    int maxi;
    maxcc = RF_NEGINF;
    for (i=0; i<dim; i++) hx[i] = 0.0;
    for (i=0; i<dim; i++) {
      hx[i] = steps[i] * m[i]; 
      cc = fabs(CovFct(hx , dim, covnr, op, param, actcov, anisotropy)); 
      if (cc>maxcc) {
	maxcc = cc;
	maxi = i; 
      }
      hx[i] = 0.0;
    }	
  } 


  if (positivedefinite || cepar->force) { 
    // correct theoretically impossible values, that are still within 
    // tolerance CIRCEMBED.tol_re/CIRCEMBED.tol_im 
    Real r, imag;
    r = imag = 0.0;    
    for(i=0,twoi=0;i<mtot;i++) {
      if (c[twoi] > 0.0) {
	c[twoi] = sqrt(c[twoi]);
      } else {
	if (c[twoi] < r) r = c[twoi];
	c[twoi] = 0.0;
      }
      {
	register Real a;
	if ((a=fabs(c[twoi+1])) > imag) imag = a;
      }
      c[twoi+1] = 0.0;
      twoi+=2;
    }
    if (GENERAL_PRINTLEVEL>1) {
      if (r<0.0 || imag>0.0) {
	PRINTF("using approximating circulant embedding:\n");
	if (r<0.0) PRINTF("\tsmallest real part has been %e \n", r);
	if (imag>0.0) 
	  PRINTF("\tlargest modulus of the imaginary part has been %e \n", imag);
      }
    }
  } else {Xerror=ERRORFAILED;goto ErrorHandling;}
  if (GENERAL_PRINTLEVEL>=10) {
    for (i=0;i<2*mtot;i++) {PRINTF("%f ",c[i]);} PRINTF("\n");
  }  
  *cc = c;
  return NOERROR;
  
 ErrorHandling:
  if (c!=NULL) {free(c);}
  return Xerror;
}



int init_circ_embed(key_type * key, int m)
{
  param_type param;
  int Xerror,i,d, start_param[MAXDIM], index_dim[MAXDIM];
  long twoRealmtot;
  double *c; 
  Real steps[MAXDIM];
  CE_storage *s;

  if (!key->grid) {Xerror=ERRORMETHODNOTALLOWED;goto ErrorHandling;}
  SET_DESTRUCT(CE_destruct);

  if ((key->S[m]=malloc(sizeof(CE_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  } 
  s =  (CE_storage*)key->S[m];
  s->c =NULL;
  s->d =NULL;
  s->local=NULL;
  FFT_NULL(&(s->FFT));
  key->destruct[m] = CE_destruct;
  
  FIRSTCHECK_COV(CircEmbed,cov,param);
  { 
    int timespacedim,v;
    bool no_last_comp;
    cov_fct *cov;
    // are methods and parameters fine ?
    for (v=0; v<actcov; v++) {
      GetTrueDim(key->anisotropy, key->timespacedim, param[v], 
		 &timespacedim, &no_last_comp, start_param, index_dim);
      cov = &(CovList[covnr[v]]);
      if ((key->Time) && !no_last_comp && (cov->isotropic==SPACEISOTROPIC))
	{timespacedim--;}
      else
	if ((cov->check!=NULL) && 
	    ((Xerror=cov->check(param[v], timespacedim, CircEmbed)))!=0)
	goto ErrorHandling;
    }
  }
  for (d=0; d<key->timespacedim; d++) {
    s->nn[d]=key->length[d]; 
    steps[d]=key->x[d][XSTEP];
  }

  if ((Xerror=internal_init_circ_embed(steps, 
				       key->anisotropy,
				       covnr, multiply, param,
				       s->nn, s->m, s->cumm, s->halfm,
				       key->timespacedim, actcov,
				       CovFct, 
				       &CIRCEMBED,
				       &(s->FFT),
				       &twoRealmtot,&c))!=0)
    goto ErrorHandling;

  // here: never replace GENERAL_STORING by key->storing
  // since, in MaxStable process sampling, GENERAL_STORING is set to true
  // whatever the value of GENERAL_STORING has been!
  if (GENERAL_STORING) {
    if ((s->d=(double *)malloc(twoRealmtot))==0){
      Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;} //d
  }
  s->c=c; 
  return 0;
  
 ErrorHandling:
  return Xerror; 
}

void internal_do_circ_embed(int *nn, int *m, int *cumm, int *halfm,
			    double *c, double *d, int Ntot, int dim,
			    FFT_storage *FFT_heap, bool add,
			    Real *res )
/* 
     implemented here only for rotationsinvariant covariance functions
     for arbitrary dimensions;
     (so it works only for even covariance functions in the sense of 
     Wood and Chan,p. 415, although they have suggested a more general 
     algorithm;) 
     
     Warning! If GENERAL_STORUNG==false when calling init_circ_embed
     and GENERAL_STORUNG==true when calling do_circ_embed, the
     complete programme will fail, since the initialization depends on
     the value of GENERAL_STORUNG
  */
{  
  int  i, j, k, HalfMp1[MAXDIM], HalfMaM[2][MAXDIM], index[MAXDIM];
  Real XX,YY,invsqrtmtot;
  bool first, free[MAXDIM+1], noexception;
  long mtot;

   
  mtot=cumm[dim-1] * m[dim-1]; 
  for (i=0; i<dim; i++) {
    HalfMp1[i] =  ((m[i] % 2)==1) ? -1 : halfm[i];
    HalfMaM[0][i] = halfm[i];
    HalfMaM[1][i] = m[i] - 1;
  }

  //if there are doubts that the algorithm below reaches all elements 
  //of the matrix, uncomment the lines containing the variable xx
  //
  //bool *xx; xx = (bool*) malloc(sizeof(bool) * mtot);
  //for (i=0; i<mtot;i++) xx[i]=true;
  
  invsqrtmtot = 1/sqrt((Real)mtot);

  if (GENERAL_PRINTLEVEL>=10) PRINTF("Creating Gaussian variables... \n");
  /* now the Gaussian r.v. have to defined and multiplied with sqrt(FFT(c))*/

  for (i=0; i<dim; i++) {
    index[i]=0; 
    free[i] = false;
  }
  free[MAXDIM] = false;
  for(;;) {      
    i = j = 0;
    noexception = false;
    for (k=0; k<dim; k++) {
      i += cumm[k] * index[k];
      if ((index[k]==0) || (index[k]==HalfMp1[k])) {
	j += cumm[k] * index[k];
      } else {
	noexception = true; // not case 2 in prop 3 of W&C
	j += cumm[k] * (m[k] - index[k]);
      }
    }
    if (GENERAL_PRINTLEVEL>=10) PRINTF("cumm...");
    i <<= 1; // since we have to index imaginary numbers
    j <<= 1;
    if (noexception) { // case 3 in prop 3 of W&C
      XX = GAUSS_RANDOM(INVSQRTTWO);
      YY = GAUSS_RANDOM(INVSQRTTWO);
      d[i] = d[i+1] = c[i];   d[i] *= XX;  d[i+1] *= YY;   
      d[j] = d[j+1] = c[j];   d[j] *= XX;  d[j+1] *= -YY; 
    } else { // case 2 in prop 3 of W&C
      d[i]   = c[i] * GAUSS_RANDOM(1.0);
      d[i+1] = 0;
    }
 
    if (GENERAL_PRINTLEVEL>=10) PRINTF("k=%d ", k);
    
    /* 
       this is the difficult part.
       We have to run over roughly half the points, but we should
       not run over variables twice (time lost)
       Due to case 2, we must include halfm.
       
       idea is:
       for (i1=0 to halfm[dim-1])
         if (i1==0) or (i1==halfm[dim-1]) then endfor2=halfm[dim-2] 
         else endfor2=m[dim-2]
         for (i2=0 to endfor2)
           if ((i1==0) or (i1==halfm[dim-1])) and
              ((i2==0) or (i2==halfm[dim-2]))
           then endfor3=halfm[dim-3] 
           else endfor3=m[dim-3]
           for (i3=0 to endfor3)
           ....
       
       i.e. the first one that is not 0 or halfm (regarded from dim-1 to 0)
       runs over 0..halfm, all the others over 0..m
       
       this is realised in the following allowing for arbitrary value of dim
       
	 free==true   <=>   endfor==m[]	 
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

  fastfourier(d, m, dim, false, FFT_heap);   
  

  /* now we correct the result of the fastfourier transformation
     by the factor 1/sqrt(mtot) and read the relevant matrix out of 
     the large vector c */
  first = true;
  for(i=0;i<dim;i++){index[i]=0;}

#define RESULT(OP)\
  for (i=0;i<Ntot;i++){ \
    j=0; for (k=0; k<dim; k++){j+=cumm[k] * index[k];} \
    res[i] OP d[2*j]*invsqrtmtot; \
    if ((fabs(d[2*j+1])>CIRCEMBED.tol_im) && \
        ((GENERAL_PRINTLEVEL>=2 && first) || GENERAL_PRINTLEVEL>=6)){ \
      PRINTF("IMAGINARY PART <> 0, %e\n",d[2*j+1]); first=false; \
    } \
    k=0; while((k<dim) && (++index[k]>=nn[k])) {index[k++]=0;} \
  }  
  if (add) {RESULT(+=)} else {RESULT(=)}
}


void do_circ_embed(key_type *key, bool add, int m, Real *res ){
  double *d;
  CE_storage *s;
  s = (CE_storage*)key->S[m];
  if (s->d==NULL) {d=s->c;} /* overwrite the intermediate 
			       result directly (algorithm 
			       allows for that) */
  else {
    assert(key->storing);
    d=s->d;
  }
  assert(key->active);

  internal_do_circ_embed(s->nn, s->m, s->cumm, s->halfm,
			 s->c, d, key->totalpoints,
			 key->timespacedim, &(s->FFT), add, res);
}



// hinr?! nicht doppeltes Feld generieren, sondern 1faches?
// da sowieso nur 1/3 verwandt wird??
// nein. an sich schon doppeltes Feld; aber bei kompakten Traeger [0,1], muss
// maxabstand der Punkte nur 1/2 betragen (durch Verdoppelung der Matrix
// wird range erreicht oder so aehnlich).
// --> beruecksichtigung bei lokaler simuation

// ?? was hatten die folgenden Kommentare zu bedeuten?
// for future development : make sure that the directions are set up correctly!
// make sure that internal_init is called with a quadratic scheme
// make sure that the relevant part is cut out correctly



int init_circ_embed_local(key_type *key, int m)
{
  int Xerror,i,d, start_param[MAXDIM], index_dim[MAXDIM],timespacedim;
  Real discretediameter,diameter, steps[MAXDIM], factor;
  double *c;
  CE_storage *s;
  long twoRealmtot;
  char actcov;
  int covnr[MAXCOV];
  int multiply[MAXCOV];

  if (!key->grid) {Xerror=ERRORMETHODNOTALLOWED;goto ErrorHandling;}
  SET_DESTRUCT(CE_destruct);

  assert(key->S[m]==NULL);
  if ((key->S[m]=malloc(sizeof(CE_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  } 
  s = (CE_storage*)key->S[m];
  s->c =NULL;
  s->d =NULL;
  s->local=NULL;
  FFT_NULL(&(s->FFT));
  key->destruct[m] = CE_destruct;

  // FC1(CircEmbedLocal,cov_loc,s->param);
  /* break the for loop directly after first finding of a local circulant 
     embedding method; check whether this function is not part of 
     a multiplicative definition -- probably this can be relaxed in future --
     currently it is unclear whether this makes sense and what the 
     extended programming would look like
  */
  actcov=0;
  { 
    int v;
    for (v=0; v<key->ncov; v++) {
      if ((key->method[v]==CircEmbedLocal) && (key->left[v])) {
      // Variance==0 is not eliminated anymore !! - maybe this could be improved
	key->left[v] = false;
	assert((key->covnr[v] >= 0) && (key->covnr[v] < currentNrCov));
	assert(key->param[v][VARIANCE] >= 0.0);
	covnr[actcov] = key->covnr[v];
	if (CovList[covnr[actcov]].cov_loc==NULL) {
	  Xerror=ERRORNOTDEFINED; goto ErrorHandling;}
	memcpy(s->param[actcov], key->param[v], sizeof(Real) * key->totalparam);
	if (actcov>0) {
	  if ((multiply[actcov-1] = key->op[v-1]) && 
	      key->method[v-1] != CircEmbedLocal){
	    if (GENERAL_PRINTLEVEL>0) 
	      PRINTF("severe error - contact author. %d %d %d %d (%s) %d (%s)\n",
		     v, key->op[v-1], key->ncov, key->method[v-1],
		     METHODNAMES[key->method[v-1]],
		     CircEmbedLocal, METHODNAMES[CircEmbedLocal]);
	    Xerror=ERRORMETHODMIX; goto ErrorHandling;
	  }
	}
  // end FC1
	if (key->op[v]) { Xerror=ERRORNOMULTIPLICATION; goto ErrorHandling; }
	else {actcov++; break;} 
  // FC2;
	actcov++;
      }
    }
  }
  if (actcov==0) { /* no covariance for the considered method found */
    if (key->traditional) Xerror=ERRORNOTINITIALIZED;
    else Xerror=NOERROR_ENDOFLIST;
    goto ErrorHandling;
  }
  // end FC2;

  { 
    int v;
    bool no_last_comp;
    cov_fct *cov;
    // are methods and parameters fine ?
    GetTrueDim(key->anisotropy, key->timespacedim, s->param[0], 
	       &timespacedim, &no_last_comp, start_param, index_dim);
    for (v=0; v<actcov; v++) {
      cov = &(CovList[covnr[v]]);
      if ((key->Time) && no_last_comp && (cov->isotropic!=FULLISOTROPIC))
	{ Xerror= ERRORWITHOUTTIME;goto ErrorHandling;}
      else {
	if ((cov->check!=NULL) && 
	    ((Xerror=cov->check(s->param[v], timespacedim, CircEmbedLocal)))!=0)
	goto ErrorHandling;
      }
    }
  }
  if (key->anisotropy) {
    Real sxx[ZWEIHOCHMAXDIM * MAXDIM], dummy, maxeigenvalue, ev;
    double D[MAXDIM];
    int d;
     GetCornersOfGrid(key, timespacedim, start_param, s->param[0], sxx);
     GetRangeCornerDistances(key, sxx, timespacedim, timespacedim,
		             &dummy, &diameter);
    { 
      Real x;
      x = 1.0 / diameter;
      // .cov uses only KAPPA out of s->param[0]
      s->param[0][VARIANCE] /= CovList[covnr[0]].cov(&x, s->param[0], 1);
      for (d=ANISO; d<key->totalparam; d++) {
	s->param[0][d] *= x;
      }
    }
  } else {
    assert(fabs(s->param[0][SCALE] * s->param[0][INVSCALE]-1.0) < EPSILON);
    diameter = 0.0;
    register Real dummy;
    dummy=0.0;
    for (d=0; d<key->timespacedim; d++) {
      dummy = key->x[d][XSTEP] * (Real) (key->length[d]-1);
      diameter += dummy * dummy; 
    }
    diameter = sqrt(diameter);

    { 
      Real x;
      x = s->param[0][SCALE] / diameter;
      s->param[0][VARIANCE] /= CovList[covnr[0]].cov(&x, s->param[0], 1);
      s->param[0][INVSCALE] = 1.0 / (s->param[0][SCALE]=diameter);
    }
  }
  if (GENERAL_PRINTLEVEL>7) PRINTF("diameter  %f \n",diameter);

  // new scale is chosen such that for the largest grid distance (diameter)
  // the local covariance function still is in the correct domain for 
  // simulating the desired variogram or covariance function (next but 
  // one line)
  //
  // The variance of the local covariance has to be chosen by def such that
  //         c * gamma_0(scale) = variance
  // if we rescale by `diameter' then 
  //         c * gamma_0(scale/diameter) != variance, i.e.,

  s->sqrt_var_squarefactor[0] = 
    // square.factor uses only KAPPA
    sqrt(4.0 * CovList[covnr[0]].square_factor(s->param[0]) *
	 s->param[0][VARIANCE]);
  // 4.0 :: 2 because of stein
  //        another 2 because sqrt( - 2.0 log(UNIFORM_RANDOM)) gives standard gaussian.
  s->totalpoints = 1;
  for (d=0; d<key->timespacedim; d++) {
    s->nn[d] = key->length[d];
    (s->totalpoints) *= s->nn[d];
    steps[d]=key->x[d][XSTEP];
  }

  if ((s->local=(Real*) malloc(sizeof(Real) * s->totalpoints))==0) 
    {Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;} 

  if ((Xerror=internal_init_circ_embed(steps, key->anisotropy,
				       covnr, multiply, s->param,
				       s->nn, s->m, s->cumm, s->halfm,
				       key->timespacedim, actcov,
				       LocalCovFct,
				       &CIRCEMBED,
				       &(s->FFT),
				       &twoRealmtot,&c))!=0)
    goto ErrorHandling;
  if (GENERAL_STORING) {
    if ((s->d=(double *)malloc(twoRealmtot))==0){
      Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;} //d
  }
  s->c=c; 
  return 0;
 
 ErrorHandling:
  return Xerror;
}

void do_circ_embed_local(key_type *key, bool add, int m, Real *res )
{  
  Real x[MAXDIM], dx[MAXDIM], factor;
  long index[MAXDIM];
  double *d, sum;
  long k,r;
  int v=0; // check if this can be relaxed
  CE_storage *s;

   s =  (CE_storage*)key->S[m];
  if (s->d==NULL) {d=s->c;} /* overwrite the intermediate 
			       result directly (algorithm 
			       allows for that) */
  else{d=s->d;}
  assert(key->active);

  internal_do_circ_embed(s->nn, s->m, s->cumm, s->halfm,
			 s->c, d, s->totalpoints,
			 key->timespacedim, &(s->FFT), false,  s->local);

  factor = INVSQRTTWO * s->sqrt_var_squarefactor[v] * s->param[v][INVSCALE];
  for (k=0; k<key->timespacedim; k++) {
    index[k]=0;
    dx[k]= GAUSS_RANDOM(factor * key->x[k][XSTEP]);
    x[k]= 0.0;
  }
  
  for(r=0;;) { 
    sum = s->local[r];
    for (k=0; k<key->timespacedim; k++)  sum += x[k]; 
    if (add) res[r++] += sum;  else res[r++] = sum; 
    k=0;
    while( (k<key->timespacedim) && (++index[k]>=key->length[k])) {
      index[k]=0;
      x[k] = 0.0;
      k++;
    }
    if (k>=key->timespacedim) break;
    x[k] += dx[k];
  }
}

