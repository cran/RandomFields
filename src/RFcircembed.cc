/*
 Authors 
 Yindeng, Jiang, jiangyindeng@gmail.com
 Martin Schlather, schlath@hsu-hh.de

 Simulation of a random field by circulant embedding
 (see Wood and Chan, or Dietrich and Newsam for the theory )

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2005 Yindeng Jiang & Martin Schlather

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

ce_param CIRCEMBED={false, false, TRIVIALSTARTEGY, -1e-7, 1e-3, 3, 20000000, 
		    0, 0, 0, 0};
local_user_param LOCAL_USER_PARAM={0, 0};

typedef struct CE_storage {
  int m[MAXDIM],halfm[MAXDIM],nn[MAXDIM],cumm[MAXDIM+1]; /* !!!! **** */  
  double *c,*d;
  double factor; /* only used in local CE */
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
    FFT_destruct(&(x->FFT));
    free(*S);
    *S = NULL;
  }
}

/*********************************************************************/
/*           CIRCULANT EMBEDDING METHOD (1994) ALGORITHM             */
/*  (it will always be refered to the paper of Wood & Chan 1994)     */
/*********************************************************************/

void SetParamCircEmbed( int *action, int *force, double *tolRe, double *tolIm,
			int *trials, int *mmin, int *userfft, int *strategy,
		        double *maxmem) 
{
  SetParamCE(action, force, tolRe, tolIm, trials, mmin, userfft, strategy,
	     maxmem, &CIRCEMBED, "CIRCEMBED");
}

void SetParamLocal( int *action, double *cutoff_a, double *intrinsic_r )
{
  switch(*action) {
  case 0 :
    LOCAL_USER_PARAM.cutoff_a=*cutoff_a; 
    if (LOCAL_USER_PARAM.cutoff_a<0) {
      LOCAL_USER_PARAM.cutoff_a=0; 
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("\nWARNING! cutoff_a had been negative.\n");
    }
    LOCAL_USER_PARAM.intrinsic_r=*intrinsic_r; 
    if (LOCAL_USER_PARAM.intrinsic_r<0 || 
	(0<LOCAL_USER_PARAM.intrinsic_r && LOCAL_USER_PARAM.intrinsic_r<1)) {
      LOCAL_USER_PARAM.intrinsic_r=0; 
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("\nWARNING! intrinsic_r had been negative or between 0 and 1\n");
    }
    break;
  case 1 :
    *cutoff_a=LOCAL_USER_PARAM.cutoff_a;
    *intrinsic_r=LOCAL_USER_PARAM.intrinsic_r;
    if (GetNotPrint) break;
  case 2 :
    PRINTF("\nLOCAL_USER_PARAM\n===============\ncutoff_a=%f\nintrinsic_r=%f\n",
	   LOCAL_USER_PARAM.cutoff_a,  LOCAL_USER_PARAM.intrinsic_r);
    break;
  default : PRINTF(" unknown action\n");
  }
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
  return fastfourier(data, m, dim, first, !first, FFT);
}

double local_cov(double *t, double *param, double *localparam, int covnr, int dim, EmbedType embed)
{
    double x = (*t) * param[DIAMETER];
    if ((*t) > localparam[LOCAL_R]) return 0.0;
    if ((*t) <= 1) 
      switch(embed)
    {
      case Cutoff:
        return CovList[covnr].variogram ? 
		1.0-CovList[covnr].cov(&x, param, dim) : CovList[covnr].cov(&x, param, dim);
      case Intrinsic:
	return localparam[INTRINSIC_A0]  + 
	    (CovList[covnr].variogram ? 
		1.0-CovList[covnr].cov(&x, param, dim) : CovList[covnr].cov(&x, param, dim)) + 
	    localparam[INTRINSIC_A2]*(*t)*(*t);
      default:
        assert(false);
    }
    else
      switch(embed)
    {
      case Cutoff:
        return localparam[CUTOFF_B] * pow( pow(localparam[CUTOFF_R], 
              localparam[CUTOFF_A])-pow(*t, localparam[CUTOFF_A]),
            2.0 * localparam[CUTOFF_A]);
      case Intrinsic:
	return localparam[INTRINSIC_B]*pow(localparam[INTRINSIC_R]-(*t),3.0)/(*t);
      default:
        assert(false);
    }
}

double LocalCovFct(double *x, int dim, int *covnr, int *op,
  param_type param, local_param_type localparam, int ncov, bool anisotropy,
  EmbedType embed) 
{
  assert(!anisotropy);
  int v, j;
  double  var,  zz, zw, result, d;
  zw = result = 0.0;
  if (dim==1) d=fabs(x[0]);
  else {
    for (d=0.0, j=0; j<dim; j++) {d += x[j] * x[j];}
    d = sqrt(d);
  }
  for (v=0; v<ncov; v++) {
    zw = var = 1.0;
    for(; v<ncov; v++) {
  zz = d * param[v][INVSCALE]/param[v][DIAMETER];
  var *= param[v][VARIANCE];
  zw *= local_cov(&zz,param[v],localparam[v],covnr[v],dim,embed);
  if (op[v]==0) break; /* v=ncov does not matter since stopped anyway */
    }
    result += var*zw;
  }
 return result;
}


int local_get_initial_m(int *nn, int *m, int dim, ce_param *cepar, double Rmax,
			double inter_scaled_spacing)
{
  double totalm;
  int i;
  
  if (GENERAL_PRINTLEVEL>=5) {
    PRINTF("calculating initial m...\n");
    PRINTF("Rmax: %f; Rmax/inter_scaled_spacing: %f\n", Rmax, 
	   Rmax/inter_scaled_spacing);
  }

  m[0]= 1 << (1 + (int) ceil(log(Rmax/inter_scaled_spacing) * INVLOG2 - 
			     EPSILON1000)); 
  for (i=1;i<dim;i++){
    m[i]=m[0];
  }
  totalm=pow(m[0], dim);
  if (totalm > cepar->maxmem) {
    sprintf(ERRORSTRING_OK, "%f", cepar->maxmem);
    sprintf(ERRORSTRING_WRONG,"%f", totalm);
    return ERRORMAXMEMORY;
  }
  else return 0;
}

int circ_embed_get_initial_m(int *nn, int *m, int dim, ce_param* cepar)
{
  double totalm;
  long i;

  totalm = 1.0;  
  for (i=0;i<dim;i++){
    int factor;
    if (cepar->userfft) {
      factor = (cepar->mmin[i] > -2) ? 2 : -cepar->mmin[i];
      m[i] =  factor * NiceFFTNumber((unsigned long) nn[i]); 
    } else {
      factor = (cepar->mmin[i] > -1) ? 1 : -cepar->mmin[i];
      m[i]= factor *
	(1 << (1 + (int) ceil(log((double) nn[i]) * INVLOG2 - EPSILON1000))); 
    }
    if (m[i]<cepar->mmin[i]) {m[i]=cepar->mmin[i];} 
    totalm *= (double) m[i];
  }
  if (totalm > cepar->maxmem) {
    sprintf(ERRORSTRING_OK, "%f", cepar->maxmem);
    sprintf(ERRORSTRING_WRONG,"%f", totalm);
    return ERRORMAXMEMORY;
  }
  else return 0;
}


int circ_embed_with_initial_m(double *steps, bool anisotropy,
			     int *covnr, int *op, param_type param,
			     int *nn, int *m, int *cumm, int *halfm, 
			     int dim, int actcov,
			     CovFctType CovFct,
			     ce_param* cepar,
			     FFT_storage *FFT,
			     long *twoRealmtot,
			     double **cc,
                             local_param_type localparam, // these two arguments 
			      //                             are used in 
                             EmbedType embed)             // local circ embed
{
  double *c;
  double hx[MAXDIM], totalm;
  int  Xerror,trials,index[MAXDIM],dummy;
  long mtot=-1,i,k,twoi;
  bool positivedefinite, cur_crit, critical, Critical[MAXDIM];

  c=NULL;
  if (GENERAL_PRINTLEVEL>=5) PRINTF("calculating the Fourier transform\n");

  positivedefinite = false;     
    /* Eq. (3.12) shows that only j\in I(m) [cf. (3.2)] is needed,
       so only the first two rows of (3.9) (without the taking the
       modulus of h in the first row)
       The following variable `index' corresponds to h(l) in the following
       way: index[l]=h[l]        if 0<=h[l]<=m[l]/2
       index[l]=h[l]-m[l]   if m[l]/2+1<=h[l]<=m[l]-1     
       Then h[l]=(index[l]+m[l]) mod m[l] !!
    */

 /* The algorithm below:
  while (!positivedefinite && (trials<cepar->trials)){
    trials++;
    calculate the covariance values "c" according to the given "m"
    fastfourier(c)
    if (!cepar->force || (trials<cepar->trials)) {
      check if positive definite
      if (!positivedefinite && (trials<cepar->trials)) {
        enlarge "m"
      }
    } else 
      print "forced"
  }
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
	  (double) ((index[k] <= halfm[k]) ? index[k] : index[k] - m[k]);
	cur_crit |= (index[k]==halfm[k]);
      }	
      dummy=i << 1;

      if (critical && cur_crit){
        c[dummy] = 0.0;
      }else{
        if (embed==Standard)
          c[dummy] = CovFct(hx, dim, covnr, op, param, actcov, anisotropy);
        else
          c[dummy] = LocalCovFct(hx, dim, covnr, op, param, localparam, 
				 actcov, anisotropy, embed);
      }

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
      long int mtot2 = mtot * 2; 
      i=0; twoi=0;
      // 16.9. < cepar.tol.im  changed to <=
      while ((twoi<mtot2) && (positivedefinite=(c[twoi]>=cepar->tol_re) && 
					  (fabs(c[twoi+1])<=cepar->tol_im)))
	{twoi+=2;}      

      if ( !positivedefinite) {
        if (GENERAL_PRINTLEVEL>=2)
	  // 1.1.71: %f changed to %e because c[twoi+1] is usually very small
	  PRINTF("non-positive eigenvalue: c[%d])=%f + %e i.\n",
		 (int) (twoi/2), c[twoi], (int) (twoi/2), c[twoi+1]);
	if (GENERAL_PRINTLEVEL>=4) { // just for printing the smallest 
	  //                            eigenvalue (min(c))
	  int index=twoi, sum=0;
	  double smallest=c[twoi];
	  char percent[]="%";
	  for (twoi=0; twoi<mtot2; twoi += 2) {
	    if (c[twoi] < 0) {
	      sum++;
	      if (c[twoi]<smallest) { index=twoi; smallest=c[index]; }
	    } 
	  }
	  PRINTF("   The smallest eigenvalue is c[%d] =  %e.\n",
		 index / 2, smallest);
	  PRINTF("   The percentage of negative eigenvalues is %5.3f%s.\n", 
		 100.0 * (double) sum / (double) mtot, percent);
	}
      }

      if (!positivedefinite && (trials<cepar->trials)) { 
        assert( embed != Cutoff && embed != Intrinsic );// in these two 
	//                                         embeddings we only do 1 trial
	
	FFT_destruct(FFT);
	free(c); c=NULL;

	totalm = 1.0;
	switch (cepar->strategy) {
	case 0 :
	  for (i=0;i<dim;i++) { /* enlarge uniformly in each direction, maybe 
				   this should be modified if the grid has 
				   different sizes in different directions */
	    m[i] <<= 1;
	    totalm *= (double) m[i];
	  }
	  break;
	case 1 :  
	  double cc, maxcc, hx[MAXDIM];
	  int maxi;
	  maxi = -1;
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
	  assert(maxi>=0);
	  m[maxi] <<= 1;
	  for (i=0;i<dim; i++) totalm *= (double) m[i];
	  break;
	default:
	  assert(false);
	}
	if (totalm>cepar->maxmem) {    
	  sprintf(ERRORSTRING_OK, "%f", cepar->maxmem);
	  sprintf(ERRORSTRING_WRONG,"%f", totalm);
	  Xerror=ERRORMAXMEMORY;
	  goto ErrorHandling;
	}
	//	assert(false);
      }
    } else {if (GENERAL_PRINTLEVEL>=2) PRINTF("forced\n");}
  }
  assert(mtot>0);
 

  if (positivedefinite || cepar->force) { 
    // correct theoretically impossible values, that are still within 
    // tolerance CIRCEMBED.tol_re/CIRCEMBED.tol_im 
    double r, imag;
    r = imag = 0.0;    
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

int internal_init_circ_embed(double *steps, bool anisotropy,
			     int *covnr, int *op, param_type param,
			     int *nn, int *m, int *cumm, int *halfm, 
			     int dim, int actcov,
			     CovFctType CovFct,
			     ce_param* cepar,
			     FFT_storage *FFT,
			     long *twoRealmtot,
			     double **cc)
{
  int  Xerror;
  if ( (Xerror=circ_embed_get_initial_m(nn, m, dim, cepar)) != 0  ||
       (Xerror=circ_embed_with_initial_m(steps, anisotropy,
				       covnr, op, param,
				       nn, m, cumm, halfm,
				       dim, actcov,
				       CovFct,
				       cepar,
				       FFT,
				       twoRealmtot,cc,
                                       NULL,
                                       Standard)) != 0 ) return Xerror;
  return 0;
}


int init_circ_embed(key_type * key, int m)
{
  param_type param;
  int Xerror=NOERROR, d, start_param[MAXDIM], index_dim[MAXDIM];
  long twoRealmtot;
  double *c; 
  double steps[MAXDIM];
  CE_storage *s;
  int multiply[MAXCOV], covnr[MAXCOV];
  unsigned short int actcov;

  if (!key->grid) {Xerror=ERRORMETHODNOTALLOWED;goto ErrorHandling;}
  SET_DESTRUCT(CE_destruct);

  if ((key->S[m]=malloc(sizeof(CE_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  } 
  s =  (CE_storage*)key->S[m];
  s->c =NULL;
  s->d =NULL;
  FFT_NULL(&(s->FFT));
  key->destruct[m] = CE_destruct;
  
  if ((Xerror = FirstCheck_Cov(key, CircEmbed, param, false, 
			      covnr, multiply, &actcov)) != NOERROR)
    goto ErrorHandling;
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

  return NOERROR;
  
 ErrorHandling:
  return Xerror; 
}

void internal_do_circ_embed(int *nn, int *m, int *cumm, int *halfm,
			    double *c, double *d, int Ntot, int dim,
			    FFT_storage *FFT_heap, bool add,
			    double *res )
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
  double XX,YY,invsqrtmtot;
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
  
  invsqrtmtot = 1/sqrt((double) mtot);

  if (GENERAL_PRINTLEVEL>=10) PRINTF("Creating Gaussian variables... \n");
  /* now the Gaussian r.v. have to defined and multiplied with sqrt(FFT(c))*/

  for (i=0; i<dim; i++) {
    index[i]=0; 
    free[i] = false;
  }
  free[MAXDIM] = false;

//  if MaxStableRF calls Cutofff, c and look uninitialised

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


void do_circ_embed(key_type *key, bool add, int m, double *res ){
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

//begin
// necessary parameters for cutoff local covariance function
void set_cutoff_param(double *localparam, double *param, cov_fct *cov,
		      int timespacedim)
{
  double firstderiv = 
    cov->derivative(&param[DIAMETER], param, -1) * param[DIAMETER];
  double phi = cov->variogram
    ? (1 - cov->cov(&param[DIAMETER],param,timespacedim))
    : cov->cov(&param[DIAMETER],param,timespacedim);
  
  assert(localparam[CUTOFF_A]>0);
  localparam[CUTOFF_R] = 
    pow(1 - 2.0*localparam[CUTOFF_A]*localparam[CUTOFF_A]*phi/firstderiv, 
	1.0/localparam[CUTOFF_A]);
  localparam[CUTOFF_B] =
    pow(-firstderiv/(2.0*localparam[CUTOFF_A]*localparam[CUTOFF_A]*phi), 
	2.0*localparam[CUTOFF_A]) * phi;
}
// necessary parameters for intrinsic local covariance function
void set_intrinsic_param(double *localparam, double *param, cov_fct *cov, 
			 int timespacedim)
{
  double firstderiv = 
    cov->derivative(&param[DIAMETER], param, -1) * param[DIAMETER];
  double secondderiv = cov->secondderivt(&param[DIAMETER], param, -1) * 
    param[DIAMETER] * param[DIAMETER];
  double phi = cov->variogram 
    ? (1 - cov->cov(&param[DIAMETER],param,timespacedim)) 
    : cov->cov(&param[DIAMETER],param,timespacedim);

  assert(localparam[INTRINSIC_R]>0);
  localparam[INTRINSIC_A0] = 
    (localparam[INTRINSIC_R]-1)/(2.0*(localparam[INTRINSIC_R]+1))*secondderiv + 
    1.0/(localparam[INTRINSIC_R]+1)*firstderiv - phi;
  if (localparam[INTRINSIC_R] == 1) 
    localparam[INTRINSIC_B] = 0;
  else
    localparam[INTRINSIC_B] = (secondderiv - firstderiv) / 
      (3.0*localparam[INTRINSIC_R] *
       (localparam[INTRINSIC_R]*localparam[INTRINSIC_R]-1.0));
  localparam[INTRINSIC_A2] =
    -0.5 * (localparam[INTRINSIC_B]*(localparam[INTRINSIC_R]-1.0) *
	    (localparam[INTRINSIC_R]-1.0)*(localparam[INTRINSIC_R]+2.0)+
	    firstderiv);

}
//end

int init_circ_embed_cutoff(key_type *key, int m)
{
  int Xerror=NOERROR,d,v;
  double diameter, steps[MAXDIM];
  double *c;
  CE_storage *s;
  long twoRealmtot;
  double Rmax;
  cov_fct *cov;
  bool sameSpacing=true; 
  double inter_scaled_spacing;
  param_type param; 
  local_strategy_type strategy, overall_strategy;
  local_param_type localparam; // the extra parameters in local cov fcts  
  int multiply[MAXCOV], covnr[MAXCOV];
  unsigned short int actcov;

  ce_param cepar=CIRCEMBED;

  if (GENERAL_PRINTLEVEL>=5)
    PRINTF("Initiating cutoff...\n");

  // check whether it's square grid
  if (!key->grid || key->anisotropy)
  { 
    // printf("\n\ngrid=%d aniso=%d \n", !key->grid, key->anisotropy);
    Xerror=ERRORMETHODNOTALLOWED;goto ErrorHandling;}
  for (d=1; d<key->timespacedim; d++) 
  {
    if( fabs(key->x[d][XSTEP] - key->x[0][XSTEP])>EPSILON )
    {
      sameSpacing=false;
      break;
    }
  }
  if (!sameSpacing)
  {
    if (GENERAL_PRINTLEVEL>=2)
      PRINTF("Currently only grids of same spacing are allowed for local circulant embedding. \n");
    Xerror=ERRORMETHODNOTALLOWED;goto ErrorHandling;
  }

  SET_DESTRUCT(CE_destruct);

  assert(key->S[m]==NULL);
  if ((key->S[m]=malloc(sizeof(CE_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  } 
  s = (CE_storage*)key->S[m];
  s->c =NULL;
  s->d =NULL;
  FFT_NULL(&(s->FFT));
  key->destruct[m] = CE_destruct;

  /* break the for loop directly after first finding of a local circulant 
     embedding method; check whether this function is not part of 
     a multiplicative definition -- probably this can be relaxed in future --
     currently it is unclear whether this makes sense and what the 
     extended programming would look like
  */
  if ((Xerror = FirstCheck_Cov(key, CircEmbedCutoff, param,  
		     false, covnr, multiply, &actcov)) != NOERROR)
    goto ErrorHandling;

  assert(!key->anisotropy);
  
  diameter = 0.0;
  register double dummy;
  dummy=0.0;
  s->totalpoints = 1;
  for (d=0; d<key->timespacedim; d++) {
    s->nn[d] = key->length[d];
    (s->totalpoints) *= s->nn[d];
    steps[d]=key->x[d][XSTEP];
    dummy = steps[d] * (double) (key->length[d]-1);
    diameter += dummy * dummy; 
  }
  diameter = sqrt(diameter);
  if (GENERAL_PRINTLEVEL>7) PRINTF("diameter  %f \n",diameter);
  
  overall_strategy=TheoGuaranteed;
  inter_scaled_spacing = steps[0]/diameter;
  Rmax = 1.0;
  if (GENERAL_PRINTLEVEL>=5)
    PRINTF("inter_scaled_spacing: %f\n", inter_scaled_spacing);

  for (v=0; v<actcov; v++) {
    assert(fabs(param[v][SCALE] * param[v][INVSCALE]-1.0) < EPSILON);
    cov = &(CovList[covnr[v]]);
    param[v][DIAMETER] = diameter * param[v][INVSCALE];
    if ((cov->check!=NULL) &&
      ((Xerror=cov->check(param[v], key->timespacedim, CircEmbedCutoff)))!=0)
        goto ErrorHandling;
    assert(cov->cutoff_strategy!=NULL);
    if (LOCAL_USER_PARAM.cutoff_a>0) 
    {
      strategy = UserSpecified;
      localparam[v][CUTOFF_A] = LOCAL_USER_PARAM.cutoff_a;
    } else {
      strategy =cov->cutoff_strategy(param[v], localparam[v], TellMeTheStrategy);
      if (strategy == CallMeAgain) {
	double temp_A, temp_R;
	set_cutoff_param( localparam[v], param[v], cov, key->timespacedim );
	temp_A = localparam[v][CUTOFF_A];
	temp_R = localparam[v][CUTOFF_R];

	strategy =cov->cutoff_strategy(param[v], localparam[v], IncreaseCutoffA);
	set_cutoff_param( localparam[v], param[v], cov, key->timespacedim );
	if (temp_R < localparam[v][CUTOFF_R]) 
	  localparam[v][CUTOFF_A] = temp_A;
      }
    }
    assert(strategy != NumeGuaranteed && strategy != SearchR);
    if ((int) overall_strategy<(int) strategy) overall_strategy=strategy;
     
    set_cutoff_param( localparam[v], param[v], cov, key->timespacedim );
    if (Rmax<localparam[v][CUTOFF_R]) Rmax=localparam[v][CUTOFF_R];
  }
  if (GENERAL_PRINTLEVEL>=5)
    PRINTF("OverallStrategy: %d; Rmax: %f\n", overall_strategy, Rmax);

  cepar.force=false; 
  cepar.strategy=TRIVIALSTARTEGY; 
  cepar.trials=1;

  if ( (Xerror=local_get_initial_m(s->nn, s->m, key->timespacedim,
        &cepar, Rmax, inter_scaled_spacing)) !=0 )
  {
    if (GENERAL_PRINTLEVEL>=2) {
      PRINTF("The size of the grid or r is too big!\n");
      for (d=0;d<key->timespacedim;d++) {PRINTF("n[%d]=%d, ",d,s->nn[d]);}
      PRINTF("\n");
      PRINTF("r=%f", Rmax);
      PRINTF("\n");
      for (d=0;d<key->timespacedim;d++) {PRINTF("m[%d]=%d, ",d,s->m[d]);}
      PRINTF("\n");
    }
    goto ErrorHandling;
  }
  if ( (Xerror=circ_embed_with_initial_m(steps, key->anisotropy,
                     covnr, multiply, param,
                     s->nn, s->m, s->cumm, s->halfm,
                     key->timespacedim, actcov,
                     NULL,
                     &cepar,
                     &(s->FFT),
                     &twoRealmtot,&c,
                     localparam,
                     Cutoff)) !=0 )
  {
    switch(overall_strategy){
	case TheoGuaranteed:
	  if (GENERAL_PRINTLEVEL>=2){
	    PRINTF("Theoretically impossible error! Please try again or contact author with the following info.\n\n");
	    printkey(key);
	  }
	  break;
	case JustTry:
	  if (GENERAL_PRINTLEVEL>=2){
	    PRINTF("Cutoff embedding does not work for the model specified.\n");
	  }
	  break;
	case UserSpecified:
	  if (GENERAL_PRINTLEVEL>=2){
	    PRINTF("The specified a does not work for the given model(s).\n");
	  }
	  break;
	default : assert(false); // the other cases are not considered
    }
    goto ErrorHandling;
  }

  
  if (GENERAL_STORING) {
    if ((s->d=(double *)malloc(twoRealmtot))==0){
      Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;}
  }
  s->c=c; 
  return 0;
 
 ErrorHandling:
  return Xerror;
}

int init_circ_embed_intrinsic(key_type *key, int m)
{
  int Xerror=NOERROR, d, v;
  double diameter, steps[MAXDIM];
  double *c;
  CE_storage *s;
  long twoRealmtot;
  double Rmax;
  cov_fct *cov;
  bool sameSpacing=true, first_iteration_of_R; 
  double inter_scaled_spacing;
  param_type param; 
  local_strategy_type strategy, overall_strategy;
  int SearchNumber, SearchIndex[MAXCOV]; // which covs need to search for 
  //                                        different values of r?
  local_param_type localparam; // the extra parameters in local cov fcts  
  int multiply[MAXCOV], covnr[MAXCOV];
  unsigned short int actcov;

  ce_param cepar=CIRCEMBED;

  if (GENERAL_PRINTLEVEL>=5)
    PRINTF("Initiating intrinsic...\n");

  // check whether it's square grid
  if (!key->grid || key->anisotropy)
  { //printf("grid=%d aniso=%d \n\n\n", !key->grid, key->anisotropy);
       Xerror=ERRORMETHODNOTALLOWED;goto ErrorHandling;}
  for (d=1; d<key->timespacedim; d++) 
  {
    if( fabs(key->x[d][XSTEP] - key->x[0][XSTEP])>EPSILON )
    {
      sameSpacing=false;
      break;
    }
  }
  if (!sameSpacing)
  {
    if (GENERAL_PRINTLEVEL>=2)
      PRINTF("Currently only grids of same spacing are allowed for local circulant embedding. \n");
    Xerror=ERRORMETHODNOTALLOWED;goto ErrorHandling;
  }

  SET_DESTRUCT(CE_destruct);

  assert(key->S[m]==NULL);
  if ((key->S[m]=malloc(sizeof(CE_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  } 
  s = (CE_storage*)key->S[m];
  s->c =NULL;
  s->d =NULL;
  FFT_NULL(&(s->FFT));
  key->destruct[m] = CE_destruct;

  /* break the for loop directly after first finding of a local circulant 
     embedding method; check whether this function is not part of 
     a multiplicative definition -- probably this can be relaxed in future --
     currently it is unclear whether this makes sense and what the 
     extended programming would look like
  */
  if ((Xerror = FirstCheck_Cov(key, CircEmbedIntrinsic, param, 
		     true, covnr, multiply, &actcov)) != NOERROR)
    goto ErrorHandling;

  assert(!key->anisotropy);
  
  diameter = 0.0;
  register double dummy;
  dummy=0.0;
  s->totalpoints = 1;
  for (d=0; d<key->timespacedim; d++) {
    s->nn[d] = key->length[d];
    (s->totalpoints) *= s->nn[d];
    steps[d]=key->x[d][XSTEP];
    dummy = steps[d] * (double) (key->length[d]-1);
    diameter += dummy * dummy; 
  }
  diameter = sqrt(diameter);
  if (GENERAL_PRINTLEVEL>7) PRINTF("diameter  %f \n",diameter);
    
  overall_strategy=TheoGuaranteed;
  SearchNumber=0;
  inter_scaled_spacing = steps[0]/diameter;
  Rmax = 1.0;
  if (GENERAL_PRINTLEVEL>=5)
    PRINTF("inter_scaled_spacing: %f\n", inter_scaled_spacing);

  for (v=0; v<actcov; v++) {
    assert(fabs(param[v][SCALE] * param[v][INVSCALE]-1.0) < EPSILON);
    cov = &(CovList[covnr[v]]);
    param[v][DIAMETER] = diameter * param[v][INVSCALE];
    if ((cov->check!=NULL) &&
      ((Xerror=cov->check(param[v], key->timespacedim, CircEmbedIntrinsic)))!=0)
        goto ErrorHandling;
    assert(cov->intrinsic_strategy!=NULL);
    if (LOCAL_USER_PARAM.intrinsic_r>0) 
    {
      strategy = UserSpecified;
      localparam[v][INTRINSIC_R] = LOCAL_USER_PARAM.intrinsic_r;
    }else{
      strategy = 
	cov->intrinsic_strategy(param[v], inter_scaled_spacing, 
				 key->timespacedim, localparam[v]); 
    }
    if ((int) overall_strategy<(int) strategy) overall_strategy=strategy;
    if (strategy==SearchR) 
      SearchIndex[SearchNumber++]=v; // In this case we'll change all the R's 
      //                                to be Rmax later
    else 
      set_intrinsic_param( localparam[v], param[v], cov, key->timespacedim );
    if (Rmax<localparam[v][INTRINSIC_R]) Rmax=localparam[v][INTRINSIC_R];
  }
  if (GENERAL_PRINTLEVEL>=5)
    PRINTF("OverallStrategy: %d; SearchNumber: %d\n",
	   overall_strategy, SearchNumber);

  cepar.force=false; 
  cepar.strategy=TRIVIALSTARTEGY; 
  cepar.trials=1;
  first_iteration_of_R=true;
  for (;;)
  {
    if ( (Xerror=local_get_initial_m(s->nn, s->m, key->timespacedim, 
          &cepar, Rmax, inter_scaled_spacing)) !=0 )
    {
      // trying greater r will not help in this case
      if (GENERAL_PRINTLEVEL>=2) {
        PRINTF("The size of the grid or r is too big!\n");
        for (d=0;d<key->timespacedim;d++) {PRINTF("n[%d]=%d, ",d,s->nn[d]);}
        PRINTF("\n");
        PRINTF("r=%f", Rmax);
        PRINTF("\n");
        for (d=0;d<key->timespacedim;d++) {PRINTF("m[%d]=%d, ",d,s->m[d]);}
        PRINTF("\n");
      }
      goto ErrorHandling;
    }
    if(SearchNumber>0) { // i.e.OverallStrategy==SearchR
      Rmax = s->m[0]*inter_scaled_spacing/2.0; // Hana's suggestion
      if (GENERAL_PRINTLEVEL>=5) {
        if (first_iteration_of_R)
	  PRINTF("First iteration: ");
	else
	  PRINTF("Second iteration: ");
	PRINTF("the actual Rmax used was: %f; Rmax/inter_scaled_spacing: %f\n", 
	       Rmax, Rmax/inter_scaled_spacing);
      }
      for (v=0;v<SearchNumber;v++) 
      {
        localparam[SearchIndex[v]][INTRINSIC_R]=Rmax;
        set_intrinsic_param( localparam[SearchIndex[v]], param[SearchIndex[v]], 
          &(CovList[covnr[SearchIndex[v]]]), key->timespacedim );
      }
    }
    Xerror=circ_embed_with_initial_m(steps, key->anisotropy,
				       covnr, multiply, param,
				       s->nn, s->m, s->cumm, s->halfm,
				       key->timespacedim, actcov,
				       NULL,
				       &cepar,
				       &(s->FFT),
				       &twoRealmtot,&c,
                                       localparam,
                                       Intrinsic);
    if (Xerror==0) break;
    else{
      switch(overall_strategy){
	  case TheoGuaranteed:
	    if (first_iteration_of_R){
	      PRINTF("\nWarning: A theoretically impossible error has occured. Please contact author with the following info.\n\n");
	      printkey(key);
	    }
	    // goto next_iteration;
            // note no break !!! hence !TheoGuaranteed in the next case
	  case NumeGuaranteed:
	    if (first_iteration_of_R && !TheoGuaranteed){ 
	      PRINTF("\nWarning: A numerical error has occured. Please contact author with the following info.\n\n");
	      printkey(key);
	    }
	    // goto next_iteration;
	  case SearchR: 
            // next_iteration:
	    if (first_iteration_of_R) {
	      first_iteration_of_R = false;
	      Rmax *= 2.0;
	      continue;
	    }
	    break;
	  case UserSpecified:
	    if (GENERAL_PRINTLEVEL>=2){
	      PRINTF("The specified r does not work for the given model(s).\n");
	    }
	    break;
	  default : assert(false); // remaining cases my not appear !
      }
      goto ErrorHandling;
    }
  }

  s->factor = 0.0;
  for (v=0; v<actcov; v++) {
    // these steps are necessary, see the two cases in 3dBrownian!
    s->factor += 2.0 * localparam[v][INTRINSIC_A2] *
      // 2.0 : see Stein (2002)
      param[v][VARIANCE];
  }
  s->factor = sqrt(s->factor)/diameter;  // standard deviation of the Gaussian 
  //                                        variables in do_...

  if (GENERAL_STORING) {
    if ((s->d=(double *)malloc(twoRealmtot))==0){
      Xerror=ERRORMEMORYALLOCATION;goto ErrorHandling;}
  }
  s->c=c; 
  PRINTF("Attention: Intrinsic embedding is to be applied, so the simulated random field will not be stationary.\n");
  return 0;
 
 ErrorHandling:
  return Xerror;
}

void do_circ_embed_intrinsic(key_type *key, bool add, int m, double *res )
{  
  double x[MAXDIM], dx[MAXDIM], *d;
  long index[MAXDIM], k, r;
  CE_storage *s;

  s =  (CE_storage*)key->S[m];
  if (s->d==NULL) {d=s->c;} /* overwrite the intermediate 
			       result directly (algorithm 
			       allows for that) */
  else{d=s->d;}
  assert(key->active);

  internal_do_circ_embed(s->nn, s->m, s->cumm, s->halfm,
			 s->c, d, s->totalpoints,
			 key->timespacedim, &(s->FFT), add,  res);

  for (k=0; k<key->timespacedim; k++) {
    index[k]=0;
    dx[k]= GAUSS_RANDOM(s->factor * key->x[k][XSTEP]);
    x[k]= 0.0;
  }
  
  for(r=0;;) { 
    for (k=0; k<key->timespacedim; k++)  res[r] += x[k]; 
    r++;
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
