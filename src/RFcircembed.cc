/*
 Authors 
 Martin Schlather, Martin.Schlather@uni-bayreuth.de 

 Simulation of a random field by circulant embedding
 (see Wood and Chan, or Dietrich and Newsam for the theory )

 Copyright (C) 2001 Martin Schlather, 

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

bool CIRCEMBED_FORCE=false; /* after CIRCEMBED_TRIALS, the result will be 
			       accepted, although not within TOL */
int  CIRCEMBED_TRIALS=3;/* how often grid is enlarged to hopefully get a 
			   positive definite matrix */
int CIRCEMBED_MMIN = 0; /* cf TBMCE_MTOT, minimum M in each direction */
#ifdef doubleReal
Real CIRCEMBED_TOL_RE=-1e-7;/* <0; real part must be greater than TOL_RE to 
			       be considered as zero */
Real CIRCEMBED_TOL_IM=1e-3;/* >0; maximum module of imaginary part that is 
			      considered as zero */
#else
Real CIRCEMBED_TOL_RE=-1e-5;
Real CIRCEMBED_TOL_IM=1e-3;
#endif


typedef struct CE_storage {
  int m[MAXDIM],halfm[MAXDIM],cumm[MAXDIM],nn[MAXDIM];
  long totallength;
  covfct cov;
  Real param[TOTAL_PARAM];
  double *c,*d,*local;
  FFT_storage FFT;
} CE_storage;

void FFT_destruct(FFT_storage *FFT) 
{
  if (FFT->work!=NULL) {free(FFT->work); FFT->work=NULL;}
  if (FFT->iwork!=NULL) {free(FFT->iwork); FFT->iwork=NULL;}
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
			int *trials, int *mmin)
{
  switch(*action) {
  case 0 :
    CIRCEMBED_TRIALS=*trials; 
    if (CIRCEMBED_TRIALS<1) {
      CIRCEMBED_TRIALS=1;
      if (GENERAL_PRINTLEVEL>0) PRINTF("\nWARNING! CIRCEMBED_TRIALS had been less than 1\n");
    }
    CIRCEMBED_FORCE=(bool)*force;
    CIRCEMBED_TOL_RE=*tolRe;
    if (CIRCEMBED_TOL_RE>0) {
      CIRCEMBED_TOL_RE=0;
      if (GENERAL_PRINTLEVEL>0) PRINTF("\nWARNING! CIRCEMBED_TOL_RE had been positive.\n");
    }
    CIRCEMBED_TOL_IM=*tolIm; 
    if (CIRCEMBED_TOL_IM<0) {
      CIRCEMBED_TOL_IM=0; 
      if (GENERAL_PRINTLEVEL>0) PRINTF("\nWARNING! CIRCEMBED_TOL_IM had been neagtive.\n");
    }
    CIRCEMBED_MMIN=*mmin;
    if (CIRCEMBED_MMIN>0) {
      CIRCEMBED_MMIN = 
	1 << (1+(int) (log(((double) CIRCEMBED_MMIN)-0.1)/log(2.0))); 
      if ((CIRCEMBED_MMIN!=*mmin) && (GENERAL_PRINTLEVEL>0))
	PRINTF("\nWARNING! CIRCEMBED_MMIN set to %d.\n",CIRCEMBED_MMIN);
    }
    break;
  case 1 :
    *force=CIRCEMBED_FORCE;
    *tolRe=CIRCEMBED_TOL_RE;
    *tolIm=CIRCEMBED_TOL_IM;
    *trials=CIRCEMBED_TRIALS;
    *mmin=CIRCEMBED_MMIN;   
    if (GetNotPrint) break;
  case 2 :
    PRINTF("\nCIRC_EMBED\n==========\nforce=%d,\ntol_Re=%e\ntol_Im=%e\ntrials=%d\nmmin=%d\n",
	   CIRCEMBED_FORCE,CIRCEMBED_TOL_RE,CIRCEMBED_TOL_IM,
	   CIRCEMBED_TRIALS,CIRCEMBED_MMIN);
    break;
  default : PRINTF(" unknown action\n");
  }
}


int fastfourier(double *data, int *m, int dim, bool first, FFT_storage *S)
#ifdef RF_GSL
  //   fftw_plan fftw_create_plan(int n, fftw_direction dir,int flags);
  //   fftwnd_plan fftwnd_create_plan(int rank, const int *n,
  //                                  fftw_direction dir, int flags);
  //   void fftwnd(fftwnd_plan plan, int howmany, fftw_complex *in,int istride, 
  //               int idist, fftw_complex *out, int ostride, int odist);
{ 
  assert(false);
}
#else
  /* this function is taken from the fft function by Robert Gentleman 
     and Ross Ihaka, in R */
{
  int inv, nseg, n,nspn,i,maxf,maxp,error;
  if (first) {
   int maxmaxf,maxmaxp;
   inv = -2;
   nseg = maxmaxf =  maxmaxp = 1;
   /* do whole loop just for error checking and maxmax[fp] .. */
   for (i = 0; i <dim; i++) {
     if (m[i] > 1) {
       fft_factor(m[i], &maxf, &maxp);
       if (maxf == 0) {error=ERRORFOURIER; goto ErrorHandling;}	
       if (maxf > maxmaxf) maxmaxf = maxf;
       if (maxp > maxmaxp) maxmaxp = maxp;
       nseg *= m[i];
     }
   }
   if ((S->work = (double*)malloc(4 * maxmaxf * sizeof(double)))==NULL) { 
     error=ERRORMEMORYALLOCATION; goto ErrorHandling; 
   }
   if ((S->iwork = (int*)malloc(maxmaxp * sizeof(int)))==NULL) { 
     error=ERRORMEMORYALLOCATION; goto ErrorHandling;
   }
   S->nseg = nseg;
   //  nseg = LENGTH(z); see loop above
  } else {
    inv = 2;
  }
  n = 1;
  nspn = 1;  
  nseg = S->nseg;
  for (i = 0; i < dim; i++) {
    if (m[i] > 1) {
      nspn *= n;
      n = m[i];
      nseg /= n;
      fft_factor(n, &maxf, &maxp);
      fft_work(&(data[0]), &(data[1]), nseg, n, nspn, inv, S->work, S->iwork);
    }
  }
  return 0;
  ErrorHandling:
  FFT_destruct(S);
  return error;
}
#endif



int internal_init_circulantembedding(key_type * key)
{
  Real *c,h,invn2[MAXDIM],*param;
  int *m, *halfm, *cumm, *nn, error,trials,index[MAXDIM],dim;
  int maxm,MAXM[]={67108864,8192,256,64};  /* the maximum grid size: 
					      MAXM[dim-1]^dim */
  long mtot,i,j,k,twoi,twoRealmtot;
  bool positivedefinite;
  covfct cov;

  c=NULL;

  if(key->cov->check!=NULL) {
    if(error=key->cov->check(key)) goto ErrorHandling;
  }
  dim=key->dim;
  maxm=MAXM[dim -1];
  nn = ((CE_storage*)key->S)->nn;
  cov = ((CE_storage*)key->S)->cov;
  param = ((CE_storage*)key->S)->param;
  m = ((CE_storage*)key->S)->m;
  cumm = ((CE_storage*)key->S)->cumm;
  halfm = ((CE_storage*)key->S)->halfm;
  
  if (GENERAL_PRINTLEVEL>=5) PRINTF("calculating the Fourier transform\n");
  
  /* cumm[i+1]=\prod_{j=0}^i m[j] 
     cumm is used for fast transforming the matrix indices into an
       index of the vector (the way the matrix is stored) corresponding
       to the matrix                   */
  /* calculate the dimensions of the matrix C, eq. (2.2) in W&C */
  
  for (i=0;i<dim;i++){
    m[i] = 1 << (1 + (int) ceil(log(((Real) nn[i])-EPSILON1000) * INVLOG2));
    if (m[i]<CIRCEMBED_MMIN) {m[i]=CIRCEMBED_MMIN;} // maybe CIRCEMBED_MMIN[i] ...
    halfm[i] = m[i]/2; 
    index[i]=1-halfm[i]; //
    if (m[i]>maxm) {error=ERRORFAILED;  goto ErrorHandling;} 
    invn2[i] = key->x[i][XSTEP] * key->x[i][XSTEP]; /* These are the nominators in (3.1).
                     But here a rectangle nn[0] * step x ... x nn[dim-1] * step
		     is used instead of the [0,1]^d cube.
		     "*step" is already squared as finally the Euclidean distance 
		     has to calculated. 
		     Here, the notation is more general than needed (for extension
		     to non-isotropic random fields  */
  }
    
  positivedefinite = false;     
    /* Eq. (3.12) shows that only j\in I(m) [cf. (3.2)] is needed,
       so only the first to rows of (3.9) (without the taking the
       modulus of h in the first row)
       The following variable `index' corresponds to h(l) in the following
       way: index[l]=h[l]        if 0<=h[l]<=m[l]/2
       index[l]=h[l]-m[l]   if m[l]/2+1<=h[l]<=m[l]-1     
       Then h[l]=(index[l]+m[l]) mod m[l] !!
    */


  /* find Fourier transform of matrix; the latter might be enlarged 
     automatically  */
  
  trials=0;
  
  while (!positivedefinite && (trials<CIRCEMBED_TRIALS)){ 
    trials++;
    cumm[0]=1; for(i=0;i<dim-1;i++){cumm[i+1]=cumm[i] * m[i];} 
    mtot=cumm[dim-1] * m[dim-1]; 
    if (GENERAL_PRINTLEVEL>=2) PRINTF("mtot=%d\n ",mtot);
    twoRealmtot = 2*sizeof(Real) * mtot;
 
    //// for the following, see the paper by Wood and Chan!
    // meaning of following variable c, see eq. (3.8)
    if ((c =(Real*) malloc(twoRealmtot)) ==0){
      error=ERRORMEMORYALLOCATION;goto ErrorHandling;
    }
    for (i=0;i<mtot;i++){ // fill in c(h) column-wise
      register int dummy;
      j=0; 
      for (k=0;k<dim;k++) {j+=cumm[k] * ((index[k]+m[k]) % m[k]);}
      // here specialisation to rotationinvariant cov.-functions: 
      h=0.0; // Euclidean distance for the right hand side of eq (3.8) 
      for (k=0;k<dim;k++) {h+=invn2[k] * (Real) (index[k] * index[k]); }
      c[dummy=2*j]=cov(sqrt(h),param); c[dummy+1]=0.0; /* cf. eq. 3.8; first Real 
							  part, then imaginary 
							  part */
      k=0;while((k<dim)&&(++index[k]>halfm[k])) {index[k]=1-halfm[k]; k++;}
    }
        
    if (GENERAL_PRINTLEVEL>=10) { // checks -- to be deleted
      int M=10000;
      for (i=0;i<mtot;i++) {
	if ((fabs(c[2*i])>M) || (fabs(c[2*i+1])>M)) {
	  PRINTF("CY(%d %f %f %f %f)\n",i,c[2*i],c[2*i+1],
		 fabs(c[2*i])>M,fabs(c[2*i+1])>M); /* exit(1); */}
      }
    }
    
    if (error=fastfourier(c,m,dim,TRUE,&(((CE_storage*)key->S)->FFT)))
      goto ErrorHandling;    

    if (GENERAL_PRINTLEVEL>=10) { // checks -- to be deleted
      int M=10000; 
      for (i=0;i<mtot;i++) {
	if ((fabs(c[2*i])>M) || (fabs(c[2*i+1])>M)) {
	  PRINTF("CX(%d %f %f %f %f)\n",i,c[2*i],c[2*i+1],
		 fabs(c[2*i])>M,fabs(c[2*i+1])>M); /* exit(1); */
	}
      }
    }    
    
    // check if positive definite. If not: enlarge and restart 
    if (!CIRCEMBED_FORCE || (trials<CIRCEMBED_TRIALS)) {
      i=0; twoi=0;
      while ((i<mtot)&&(positivedefinite=(c[twoi]>=CIRCEMBED_TOL_RE&& 
					  fabs(c[twoi+1])<CIRCEMBED_TOL_IM)))
	{i++; twoi+=2;}
      if (!positivedefinite) {
	if (GENERAL_PRINTLEVEL>=2)
	  PRINTF(" nonpos %d  %f %f \n",i,c[twoi],c[twoi+1]);
	free(c); c=NULL;
	FFT_destruct(&(((CE_storage*)key->S)->FFT));
	for (i=0;i<dim;i++) { /* enlarge uniformly in each direction, maybe this 
				 should be modified if the grid has different 
				 sizes in different directions */
	  m[i] <<= 1; 
	  halfm[i]=m[i]/2; 
	  index[i]=1-halfm[i]; 
	  if (m[i]>maxm) {error=ERRORFAILED;goto ErrorHandling;}
	}           
      }
    } else {if (GENERAL_PRINTLEVEL>=2) PRINTF("forced\n");}
  }    
  if (positivedefinite || CIRCEMBED_FORCE) { 
    // correct theoretically impossible values, that are still within 
    // tolerance CIRCEMBED_TOL_RE/CIRCEMBED_TOL_IM 
    for(i=0,twoi=0;i<mtot;i++) {
      c[twoi] = (c[twoi]>0.0 ? c[twoi] : 0.0);
      c[twoi+1] = 0.0;
      twoi+=2;
    }
  } else {error=ERRORFAILED;goto ErrorHandling;}
  if (GENERAL_STORING) {
    if ((((CE_storage*)key->S)->d=(double *)malloc(twoRealmtot))==0){
      error=ERRORMEMORYALLOCATION;goto ErrorHandling;} //d
  }
  ((CE_storage*)key->S)->c=c; 
  
  if (GENERAL_PRINTLEVEL>=10) {
    for (i=0;i<2*mtot;i++) {PRINTF("%f ",c[i]);} PRINTF("\n");
  }
  return 0;
  
 ErrorHandling:
  if (c!=NULL) {free(c);}
  assert(key->destruct!=NULL);
  key->destruct(&key->S);
  key->destruct=NULL;
  return error;
}



int init_circulantembedding(key_type * key)
{
  int error,i,d;
  assert(key->active==false);
  assert((key->covnr>=0) &&(key->covnr<currentNrCov));
  assert(key->cov->cov!=NULL); // I do not know an example yet, 
  //                              where cov is unknown, but should be simulated
  assert(key->param[VARIANCE]>=0.0);
  assert(key->param[SILL]>=key->param[VARIANCE]); 
  assert( fabs(key->param[SCALE]*key->param[INVSCALE]-1.0) < EPSILON);
  assert((key->S==NULL) && (key->destruct==NULL));

  if ((key->S=malloc(sizeof(CE_storage)))==0){
    error=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  ((CE_storage*)key->S)->c =NULL;
  ((CE_storage*)key->S)->d =NULL;
  ((CE_storage*)key->S)->local=NULL;
  FFT_NULL(&(((CE_storage*)key->S)->FFT));
  key->destruct = CE_destruct;
  
  // this circulant embedding procedure is programmed only for regular grids 
  // and isotropic covariance functions
  if (!key->grid) {error=ERRORMETHODNOTALLOWED;goto ErrorHandling;}
  for (d=1; d<key->dim; d++){
    if ((key->length[d]!=1) && // length.==1 iff grid is of lower dimension
	((key->x[d][XSTEP]!=key->x[0][XSTEP]) || 
	 (key->length[d]!=key->length[0]))) 
	// [d[2] is the step length 
	// of the grid, check which dimension
    {
      error=ERRORNOTPROGRAMMED; goto ErrorHandling;; 
    }
  }
	
  for (i=0; i<TOTAL_PARAM; i++) {((CE_storage*)key->S )->param[i]=key->param[i];} 
  ((CE_storage*)key->S )->cov = key->cov->cov;

  // this is still desastrous!--- in do_circ, an enlarged res should be used then
  // so write an internal do_circ function

  for (d=0; d<key->dim; d++) {
    ((CE_storage*)key->S)->nn[d]=key->length[d]; 
  }
  ((CE_storage*)key->S)->totallength=key->totallength;
  
  if (error=internal_init_circulantembedding(key)) goto ErrorHandling;
  return 0;
  
 ErrorHandling:
  if (key->destruct!=NULL) key->destruct(&key->S);
  key->destruct=NULL;
  key->active = false;
  return error; 
}

void do_circulantembedding(key_type *key, Real *res END_WITH_GSLRNG)
/* 
     implemented here only for rotationsinvariant covariance functions
     for arbitrary dimensions;
     (so it works only for even covariance functions in the sense of 
     Wood and Chan,p. 415, although they have suggested a more general 
     algorithm;) 
     
     Warning! If GENERAL_STORUNG==false when calling init_circulantembedding
     and GENERAL_STORUNG==true when calling do_circulantembedding, the
     complete programme will fail, since the initialization depends on
     the value of GENERAL_STORUNG
  */
{
  int   i,*halfm, Ntot, k, kk,*m, *nn,
    *cumm,index[MAXDIM],halfmtot,lastnonzeroORmiddle,modtot;
  Real UU,VV,XX,YY,invsqrtmtot;
  Real *c,*d;  
  bool zeroORmiddle,first;
  int error,dim;
  long mtot,twoRealmtot;

#ifdef RF_GSL
  assert(RANDOM!=NULL);
#endif
  assert(key->active);
  c=((CE_storage*)key->S)->c;
  m= ((CE_storage*)key->S)->m; 
  dim=key->dim;
  halfm=((CE_storage*)key->S)->halfm; 
  cumm=((CE_storage*)key->S)->cumm; //// 
  
  mtot=cumm[dim-1] * m[dim-1]; 
  twoRealmtot = 2*sizeof(Real) * mtot;
  invsqrtmtot = 1/sqrt((Real)mtot);
  halfmtot =  mtot /2;
  modtot = cumm[dim-1] * (halfm[dim-1]+1); 
  Ntot=((CE_storage*)key->S)->totallength;
  nn = ((CE_storage*)key->S)->nn;

  if (((CE_storage*)key->S)->d==NULL) {d=c;} /* overwrite the intermediate 
						result directly (algorithm 
						allows for that) */
  else{d=((CE_storage*)key->S)->d;}

  if (GENERAL_PRINTLEVEL>=5) PRINTF("Creating Gaussian variables... \n");
  /* now the Gaussian r.v. have to defined and multiplied with sqrt(FFT(c))*/
  for(i=0;i<dim;i++){index[i]=0;} // superfluent?
    
  //// for details see Wood and Chan:
  for (i=0; i<modtot; i++) { /* about half of the loop is clearly unnecessary, 
          as pairs of r.v. are created, see inequality * below; 
          therefore modtot is used instead of mtot;
          as the procedure is complicated, it is easier to let the loop run
          over a bit more values of i than using exactly mtot/2 values */ 
    zeroORmiddle = true;
    lastnonzeroORmiddle=0; /* theoretically, it is undefined, if
                              index[k] = (0 or halfm[k])  for all k;
                              but this would lead an error in the inequality *
                              below; in the case all the indices are 0 or
                              halfm[k] then inequality * below is obviously
                              satisfied with lastnonzeroORmiddle=0
                           */            
    kk=0; /* kk takes the role of k[l] in Prop. 3; but here kk is already 
             transformed from matrix indices to a vector index */  
    for (k=0; k<dim; k++){
       if (index[k]==0 || index[k]==halfm[k]) {kk += cumm[k] * index[k];}
       else {
         zeroORmiddle=false;
         kk += cumm[k] * (m[k]-index[k]);
         lastnonzeroORmiddle=k;
       }
    }
    if (index[lastnonzeroORmiddle]<=halfm[lastnonzeroORmiddle]){ /* * */
       /* By Prop. 3, pairs of r.v. exist which are (neg.) dependent;
          so, these pairs are created together what is done in
          the following. The above inequality  guaratees, that
          those pairs are created only once ---
          
          Assume m=4 for all k and dimension=3; then the above inequality 
          becomes clear considering the following pair of tripels whose 
          corresponding random variables are simulated simultaneously.
          [(0, 0, 0); (0, 0, 2)];
          [(0, 1, 0); (0, 3, 0)]; 
          [(0, 3, 0); ...] is ruled out ! 
          [(1, 1, 1); (3, 3, 3)];
          [(1, 3, 1); (3, 1, 3)];
          [(3, 1, 3); ...] is ruled out!
          [(3, 3, 3); ...] is ruled out!
    
          [(0, 0, 2); ...] is NOT ruled out by *, but by the following 
                           inequality (i<halfmtot)
       */

      if (zeroORmiddle){ /* step 2 of Prop. 3 */
        if (i<halfmtot) { 
	  int twoi;
        /* create two Gaussian random variables XX & YY  with variance 1/2 */ 
	  UU= sqrt(-log(1.0 - UNIFORM_RANDOM));
	  VV =  TWOPI * UNIFORM_RANDOM; 
          XX = UU * sin(VV);
          YY = UU * cos(VV);
         
          kk = (i + halfmtot)<<1;
	  twoi = i << 1;
          d[twoi] = sqrt(c[twoi]) * XX; d[twoi+1]=0.0; 
          d[kk] = sqrt(c[kk]) * YY; d[kk+1]=0.0;
	}
      } else { /* step 3 of Prop. 3 */
	int twoi, twok;
	UU= sqrt(-log(1.0 - UNIFORM_RANDOM));
	VV =  TWOPI * UNIFORM_RANDOM; 
	XX = UU * sin(VV);
        YY = UU * cos(VV);	
	twoi = i << 1;
	twok = kk<< 1;

        d[twoi+1]=d[2*i]=sqrt(c[twoi]); d[twoi]*=XX; d[twoi+1]*=YY;
        d[twok+1]=d[2*kk]=sqrt(c[twok]); d[twok]*=XX; d[twok+1]*= -YY;
      }
    }
    k=0; while((k<dim) && (++index[k]>=m[k])) {index[k++]=0;}
  }

  fastfourier(d,m,dim,FALSE,&(((CE_storage*)key->S)->FFT));   

  /* now we correct the result of the fastfourier transformation
     by the factor 1/sqrt(mtot) and read the relevant matrix out of 
     the large vector c */
  first = true;
  for(i=0;i<dim;i++){index[i]=0;}
  for (i=0;i<Ntot;i++){
    kk=0;for (k=0; k<dim; k++){kk+=cumm[k] * index[k];}	 
    res[i]=d[2*kk]*invsqrtmtot;
    if ((fabs(d[2*kk+1])>CIRCEMBED_TOL_IM) && 
	((GENERAL_PRINTLEVEL>=2 && first) || GENERAL_PRINTLEVEL>=3)){
      PRINTF("IMAGINARY PART <> 0, %f\n",d[2*kk+1]); first=false;
    }
    k=0; while((k<dim) && (++index[k]>=nn[k])) {index[k++]=0; }
  }

  if (!(key->active=(GENERAL_STORING && (((CE_storage*)key->S)->d!=NULL)))) { 
    assert(key->destruct!=NULL);
    key->destruct(&key->S);
    key->destruct=NULL;
  }
}

// see /home/martin/article/C/RFce_local.cc
int init_circ_embed_local(key_type * key) {assert(false); /* not public yet */}
// make sure that internal_init is called with a quadratic scheme
// make sure that the relevant part is cut out correctly

// hinr?! nicht doppeltes Feld generieren, sondern 1faches?
// da sowieso nur 1/3 verwandt wird??
// nein. an sich schon doppeltes Feld; aber bei kompakten Traeger [0,1], muss
// maxabstand der Punkte nur 1/2 betragen (durch Verdoppelung der Matrix
// wird range erreicht oder so aehnlich).
// --> beruecksichtigung bei lokaler simuation

//
// for future development : make sure that the directions are set up correctly!
void do_circ_embed_local(key_type *key, Real *res END_WITH_GSLRNG) {
  assert(false); /* not public yet */ 
}








