
/*
 Authors
 Martin Schlather, Martin.Schlather@uni-bayreuth.de

 Simulation of a random field by turning bands;
 see RFspectral.cc for spectral turning bands

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
#include <sys/timeb.h>
#include <unistd.h>
#include <assert.h>
#include "RFsimu.h"
#include "RFCovFcts.h"

#define MAXNN 100000000.0 /* max number of points simulated on a line */

int TBM2_LINES=60;// number of lines simulated
int TBM3D2_LINES=500;
int TBM3D3_LINES=500;
Real TBM2_LINESIMUFACTOR=2.0; /*factor by which simulation on the line is finer 
				than the grid */
Real TBM2_LINESIMUSTEP=0.0; /* grid lag defined absolutely */
Real TBM3D2_LINESIMUFACTOR=2.0;
Real TBM3D2_LINESIMUSTEP=0.0;
Real TBM3D3_LINESIMUSTEP=0.0;
Real TBM3D3_LINESIMUFACTOR=2.0;
bool TBMCE_FORCE=false; /* if TRIALS times failed then do as if it was OK?  
			   see also CIRCEMBED_* */
int  TBMCE_TRIALS=3;/* how often size is doubled in case of failure before 
		       giving up */
int  TBMCE_MTOT=0;  /* minimal mtot */
#ifdef doubleReal 
Real TBMCE_TOL_RE=-1e-7;
Real TBMCE_TOL_IM=1e-3;
#else
Real TBMCE_TOL_RE=-1e-5;
Real TBMCE_TOL_IM=1e-3;
#endif
SimulationType TBM_METHOD=Nothing;

typedef struct TBM_storage {
  simu1dim method; /* simulation method on the line (up to now only 
		      circulantembedding) */
  int nn;
  Real linesimuscale,half[MAXDIM];
  int mtot;
  double *c,*d,*simuline;
  FFT_storage FFT;
} TBM_storage;

void TBM_destruct(void **S) 
{
  if (*S!=NULL) {
    TBM_storage *x;
    x = *((TBM_storage**)S);
    if (x->c!=NULL) free(x->c);
    if (x->d!=NULL) free(x->d);
    if (x->simuline!=NULL) free(x->simuline);
    FFT_destruct(&(x->FFT));
    free(*S);
    *S = NULL;
  }
}

// the following two procedures are the same as in RFcircembed, except simplified
// to dimension 1
void do_circulantembedding1dim(key_type *key,Real *res END_WITH_GSLRNG)
{ 
  int dim, i, kk,
    index,halfmtot,zeroORmiddle,
    modtot,mtot;
  Real *d,*c,UU,VV,XX,YY,invsqrtmtot;
  int TwoIndex,TwoKK,twoRealmtot,nn;
  //  int OneDimMaxm=67108864;
  bool first;

#ifdef RF_GSL
  assert(RANDOM!=NULL);
#endif
  assert(key->active);
  nn = ((TBM_storage*)key->S)->nn;
  mtot = ((TBM_storage*)key->S)->mtot;
  // PRINTF(" %d mtot =%d \n",nn,mtot);
  c= ((TBM_storage*)key->S)->c;
  d= ((TBM_storage*)key->S)->d;
 
  twoRealmtot = 2*sizeof(Real) * mtot;
  halfmtot = mtot/2;
  modtot = halfmtot + 1; 
  invsqrtmtot = 1/sqrt((Real) mtot);
  dim=1;
      
  /* now the Gaussian r.v. have to defined and multiplied with sqrt(FFT(c))*/

  for (index=0; index<modtot; index++) { 
    /* about half of the loop is clearly unnecessary,    
       as pairs of r.v. are created, see inequality * below; 
       therefore modtot is used instead of mtot;
       as the procedure is complicated, it is easier to let the loop run
       over a bit more values of i than using exactly mtot/2 values */ 
    zeroORmiddle = 1;
    kk=0; /* kk takes the role of k[l] in Prop. 3; but here kk is already 
             transformed from matrix indices to a vector index */  

    if (index==0 || index==halfmtot) {kk += index;}
    else {
      zeroORmiddle=0;
      kk += mtot-index;
    }
    
    if (zeroORmiddle){ /* step 2 of Prop. 3 */
      if (index<halfmtot) { 
	/* create two Gaussian random variables XX & YY  with variance 1/2 */
	UU= sqrt(-log(1.0 - UNIFORM_RANDOM));
	VV = TWOPI * UNIFORM_RANDOM;
	XX = UU * sin(VV);
	YY = UU * cos(VV);

	kk = index + halfmtot;
	TwoIndex = index <<1;
	d[TwoIndex] = c[TwoIndex] * XX; d[TwoIndex+1]=0.0;
	TwoKK=kk <<1;
	d[TwoKK] = c[TwoKK] * YY; d[TwoKK+1]=0.0;
      }
    } else { /* step 3 of Prop. 3 */
      UU= sqrt(-log(1.0 - UNIFORM_RANDOM));
      VV = TWOPI * UNIFORM_RANDOM;
      XX = UU * sin(VV);
      YY = UU * cos(VV);

      TwoIndex= index <<1;
      d[TwoIndex+1]=(d[TwoIndex]=c[TwoIndex])*YY; d[TwoIndex]*=XX;       
      TwoKK=kk <<1;
      d[TwoKK+1]=(d[TwoKK]=c[TwoKK]) * (-YY);    d[TwoKK]*=XX; 
    }
  }

  fastfourier(d,&mtot,dim,FALSE,&(((TBM_storage*)key->S)->FFT));   

  first =true;
  for (kk=i=0;i<nn;i++,kk++){
    res[i]= (Real) d[kk++] *invsqrtmtot; 
    if ((fabs(d[kk])>TBMCE_TOL_IM) && (first)){
       if (GENERAL_PRINTLEVEL>=2) 
	 PRINTF("IMAGINARY PART <> 0; %d %f \n",kk,d[kk]); first=false;
    }
  }
}

int init_circulantembedding1dim(key_type *key, Real *p)
{ 
  int dim=1, i,j,index,halfmtot, mtot,TwoI,trials,nn;
  Real *c;
  int twoRealmtot,error;
  int MAXM=67108864;
  bool positivedefinite;
  covfct cov;

  c = NULL;
  nn = ((TBM_storage*)key->S)->nn;
  mtot = 1 << (1 + (int) ceil(log(((Real) nn)-EPSILON1000) * INVLOG2)); // nn<=2^(mtot-1)
  if (TBMCE_MTOT>mtot) {mtot=TBMCE_MTOT;}
 
  if (mtot>MAXM) { error=ERRORFAILED; goto ErrorHandling;} 

  if (key->method==TBM2) { cov=key->cov->cov_tbm2;}
  else { 
    assert(key->method==TBM3); 
    cov=key->cov->cov_tbm3;
  }

  positivedefinite = false;     
  trials =0;

  while (!positivedefinite && (trials<TBMCE_TRIALS)){
    halfmtot = mtot/2;  
    index=1-halfmtot;   
    trials++; 
   
    if (GENERAL_PRINTLEVEL>=2) {
      PRINTF("trial %d of %d\n",trials,TBMCE_TRIALS);
      PRINTF("mtot=%d\n",mtot);
    }
    twoRealmtot = 2*sizeof(Real) * mtot;
    
    if ((c =(Real*) malloc(twoRealmtot)) ==0){
      error=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    
    for (i=0;i<mtot;i++){ /* fill in c(h) column-wise*/
      j= 2 * ((index+mtot) % mtot);  
      c[j]=cov((Real) index,p); c[j+1]=0; 
      if (++index>halfmtot) {index=1-halfmtot;}
    }
    
     
    if (error=(fastfourier(c,&mtot,(int) dim,TRUE,&(((TBM_storage*)key->S)->FFT))))
      goto ErrorHandling; 

    /* check if positive definite. If not enlarge and restart */
    if (!TBMCE_FORCE || (trials<TBMCE_TRIALS)) {
      i=0;      TwoI = 0;
      while ((i<mtot)&&(positivedefinite=(c[TwoI]>=TBMCE_TOL_RE && 
					  fabs(c[TwoI+1])<TBMCE_TOL_IM)))
	{// Approximation if in tolerance bounds
	  i++;  TwoI += 2;
	}
      if (!positivedefinite) {
	if (GENERAL_PRINTLEVEL>=2) 
	  PRINTF(" nonpos %d  %f %f \n",i,c[TwoI],c[TwoI+1]);
	free(c); c=NULL;
	mtot <<= 1; 
	if (mtot>MAXM) { error=ERRORFAILED; goto ErrorHandling;}
      }           
    } else {if (GENERAL_PRINTLEVEL>=2) PRINTF("forced\n");}
  }
  // note: mtot<<=1 above; so now mtot equals number of allocated units of d
  
  if (positivedefinite || TBMCE_FORCE) {
      for(i=0,TwoI=0;i<mtot;i++,TwoI+=2) {
	c[TwoI] = (c[TwoI]>0.0 ? sqrt(c[TwoI]) : 0.0);
	c[TwoI+1] = 0.0;
      }
  } else {
    error=ERRORFAILED; goto ErrorHandling;
  }

  if ((((TBM_storage*)key->S)->simuline= (Real *)malloc(nn * sizeof(Real)))==0){
    error=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  if ((((TBM_storage*)key->S)->d = (Real *)malloc(twoRealmtot)) ==0){
    error=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  ((TBM_storage*)key->S)->mtot = mtot;
  ((TBM_storage*)key->S)->c = c;  
  return 0;
  
 ErrorHandling:
  if (c!=NULL) free(c);
  assert(key->destruct!=NULL);
  key->destruct(&key->S);
  key->destruct=NULL;
  return error;
}


void SetParamTBMCE(int *action,int *force,Real *tolRe,Real *tolIm,int *trials, 
		   int *mtot) {
  switch (*action) {
  case 0 :
    TBMCE_FORCE=(bool) *force;
    TBMCE_TRIALS = *trials;
    if (TBMCE_TRIALS<1) {
      TBMCE_TRIALS=1;
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("\nWARNING! TBMCE_TRIALS had been less than 1\n");
    }
    TBMCE_TOL_RE=*tolRe;
    if (TBMCE_TOL_RE>0) {
      TBMCE_TOL_RE=0; 
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("\nWARNING! TBMCE_TOL_RE had been positive.\n");
    }
    TBMCE_TOL_IM=*tolIm; 
    if (TBMCE_TOL_IM<0) {
      TBMCE_TOL_IM=0; 
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("\nWARNING! TBMCE_TOL_IM had been neagtive.\n");
    }
    TBMCE_MTOT=*mtot;
    if (TBMCE_MTOT>0) {
      TBMCE_MTOT = 1 << (1+(int) (log(((double) TBMCE_MTOT)-0.1)/log(2.0))); 
      if ((TBMCE_MTOT!=*mtot)  && (GENERAL_PRINTLEVEL>0))
	PRINTF("\nWARNING! TBMCE_MTOT set to %d.\n",TBMCE_MTOT);
    }
   break;
  case 1 : 
    *force=TBMCE_FORCE;
    *tolRe=TBMCE_TOL_RE;
    *tolIm=TBMCE_TOL_IM;
    *trials=TBMCE_TRIALS;
    *mtot=TBMCE_MTOT;
    if (GetNotPrint)  break;
  case 2 :
    PRINTF("\nTBM_CE\n======\nforce=%d\ntol_Re=%e\ntol_Im=%e\ntrials=%d\nmtot=%d\n",
	   TBMCE_FORCE,TBMCE_TOL_RE,TBMCE_TOL_IM,TBMCE_TRIALS,TBMCE_MTOT);
     break;
  default : PRINTF(" unknown action\n");  
  }
}

void SetParamTBM2(int *action,int *nLines, Real *linesimufactor, 
		  Real *linesimustep) {
  switch (*action) {
  case 0 :
    TBM2_LINES=*nLines;
    if ((*linesimufactor>0.0) ^ (*linesimustep>0.0)) {
      TBM2_LINESIMUFACTOR=*linesimufactor;
      TBM2_LINESIMUSTEP=*linesimustep;    
    } else { 
      PRINTF(" Exactly one of `linesimufactor' and `linesimustep' should be positive."); 
    }
    break;
  case 1 :
    *nLines=TBM2_LINES;
    *linesimufactor=TBM2_LINESIMUFACTOR;
    *linesimustep=TBM2_LINESIMUSTEP;        
       if (GetNotPrint)  break;
  case 2:
   PRINTF("\nTBM2\n====\nLines=%d\nlinesimufactor=%f\nlinesimustep=%f\n",
	   TBM2_LINES,TBM2_LINESIMUFACTOR,TBM2_LINESIMUSTEP);
    break;
  default : PRINTF(" unknown action\n"); 
  }
}

void SetParamTBM3D2(int *action,int *nLines, Real *linesimufactor, 
		    Real *linesimustep){
  switch (*action) {
  case 0 :
    TBM3D2_LINES=*nLines;
    if ((*linesimufactor>0.0)^(*linesimustep>0.0)) {
      TBM3D2_LINESIMUFACTOR=*linesimufactor;
      TBM3D2_LINESIMUSTEP=*linesimustep;
    } else { 
      PRINTF(" Exactly one of `linesimufactor' and `linesimustep' should be positive"); 
    }
    break;
  case 1 :
    *nLines=TBM3D2_LINES;
    *linesimufactor=TBM3D2_LINESIMUFACTOR;
    *linesimustep=TBM3D2_LINESIMUSTEP;
     if (GetNotPrint)   break;
  case 2:
    PRINTF("\nTBM3D2\n======\nLines=%d\nlinesimufactor=%f\nlinesimustep=%f\n",
	   TBM3D2_LINES,TBM3D2_LINESIMUFACTOR,TBM3D2_LINESIMUSTEP);
    break;
  default : PRINTF(" unknown action\n"); 
  }
}

void SetParamTBM3D3(int *action,int *nLines, Real *linesimufactor, 
		    Real *linesimustep){
  switch (*action) {
  case 0 :
    TBM3D3_LINES=*nLines;
    if ((*linesimufactor>0.0)^(*linesimustep>0.0)) {
      TBM3D3_LINESIMUFACTOR=*linesimufactor;
      TBM3D3_LINESIMUSTEP=*linesimustep;
    } else { 
      PRINTF(" Exactly one of `linesimufactor' and `linesimustep' should be positive"); 
    }
    break;
  case 1 :
    *nLines=TBM3D3_LINES;
    *linesimufactor=TBM3D3_LINESIMUFACTOR;
    *linesimustep=TBM3D3_LINESIMUSTEP;
    if (GetNotPrint)  break;
  case 2:
    PRINTF("\nTBM3D3\n======\nLines=%d\nlinesimufactor=%f\nlinesimustep=%f\n",
	   TBM3D3_LINES,TBM3D3_LINESIMUFACTOR,TBM3D3_LINESIMUSTEP);
    break;
  default : PRINTF(" unknown action\n"); 
  }
}

void SetParamTBM(int *action, int *tbm_method) {
  switch (*action) {
  case 0 :
    if ((*tbm_method<=(int) Nothing) && (*tbm_method>=0)) {      
      SimulationType m;
      m=(SimulationType) *tbm_method;
      if ((m!=TBM2)&&(m!=TBM3)&&(m!=SpectralTBM)){
	TBM_METHOD = m;
	return;
      }
    }  
    if (GENERAL_PRINTLEVEL>0) {
      PRINTF("Specified method [%d] is not allowed. Method set to `Nothing'\n",
	     *tbm_method);
    }
    TBM_METHOD = Nothing;
    break;
  case 1 :
    *tbm_method = (int) TBM_METHOD;
    if (GetNotPrint)  break;
  case 2 :
      PRINTF("\nTBM:\n====\nuser defined method=%d\n",TBM_METHOD);
   break;
  default : PRINTF(" unknown action\n"); 
  }
}



int init_turningbands(key_type *key)
{
  int error,i,d;
  Real p[TOTAL_PARAM],linesimufactor,linesimustep,dist,mindelta;

  assert((key->covnr>=0) &&(key->covnr<currentNrCov));
  assert(key->active==false);
  assert(key->param[VARIANCE]>=0.0);
  assert(key->param[SILL]>=key->param[VARIANCE]); 
  assert( fabs(key->param[SCALE]*key->param[INVSCALE]-1.0) < EPSILON);
  assert((key->S==NULL) && (key->destruct==NULL));
  
  if(key->cov->check!=NULL) if(error=key->cov->check(key)) goto ErrorHandling;
  if ((key->S=malloc(sizeof(TBM_storage)))==0){
    error=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  ((TBM_storage*)key->S)->c =NULL;
  ((TBM_storage*)key->S)->d =NULL;
  ((TBM_storage*)key->S)->simuline=NULL;
  FFT_NULL(&(((TBM_storage*)key->S)->FFT));
  key->destruct = TBM_destruct;


  if (key->dim==2) {
    if (key->method==TBM2) {
      linesimufactor=TBM2_LINESIMUFACTOR;
      linesimustep = TBM2_LINESIMUSTEP;
    }
    else {
      assert(key->method==TBM3);
      linesimufactor=TBM3D2_LINESIMUFACTOR;
      linesimustep = TBM3D2_LINESIMUSTEP;
    }
  } else if (key->dim==3) {
      assert(key->method==TBM3);
      linesimufactor=TBM3D3_LINESIMUFACTOR;
      linesimustep = TBM3D3_LINESIMUSTEP;    
  } else assert(false);
  

  for (i=0;i<TOTAL_PARAM;i++) {p[i]=key->param[i];}
  if (key->grid) {
    {
      register Real dummy;
      dist = 0.0;
      for (d=0; d<key->dim; d++) {
	dummy = key->x[d][XSTEP] * (Real) (key->length[d] - 1);
	dist += dummy * dummy;
      }
    }
    if (linesimustep>0.0) { 
      ((TBM_storage*)key->S)->linesimuscale =  1/linesimustep; 
    }
    else {
      mindelta = RF_INF;
      for (d=0; d<key->dim; d++) {
	if ((key->x[d][XSTEP]<mindelta) && (key->x[d][XSTEP]>0)) 
	  {mindelta=key->x[d][XSTEP];}
      }
      ((TBM_storage*)key->S)->linesimuscale = linesimufactor/mindelta;
    }
  } else { // not a grid
    Real min[MAXDIM],max[MAXDIM]; 
    int dim,j,d;
    dim = key->dim;
    if (linesimustep>0.0) {
      ((TBM_storage*)key->S)->linesimuscale =  1/linesimustep; 
    }
    else {      
      mindelta = RF_INF; 
      for (i=0;i<key->totallength;i++) {
	for (j=0;j<i; j++) {
	  register Real diff,dist;
	  for(dist=0.0,d=0;d<dim;d++) {
	    diff = key->x[d][i]- key->x[d][j]; 
	    dist += diff * diff;
	  }
	  if (dist<mindelta) mindelta=dist; 
	}
      }
      ((TBM_storage*)key->S)->linesimuscale = linesimufactor/sqrt(mindelta);
    }
    for (d=0;d<MAXDIM;d++) {min[d]=RF_INF; max[d]=RF_NEGINF;}
    for (i=0;i<key->totallength;i++) {
      for (d=0;d<dim;d++) {
	if (key->x[d][i]<min[d]) min[d] = key->x[d][i];
	if (key->x[d][i]>max[d]) max[d] = key->x[d][i];
      }
    }    
    dist = 0.0;
    for (d=0; d<key->dim; d++) {
      dist += (max[d]-min[d])*(max[d]-min[d]);
      ((TBM_storage*)key->S)->half[d] = (max[d]+min[d])/2.0;
    }  
  }  
  { 
    Real dummy;
    dummy = (sqrt(dist)* ((TBM_storage*)key->S)->linesimuscale);
    if (dummy > MAXNN) {error=ERRORNN; goto ErrorHandling;}
    ((TBM_storage*)key->S)->nn= 5+ (long) dummy;//+5 for safety
  }
  

  p[SCALE]=key->param[SCALE] * ((TBM_storage*)key->S)->linesimuscale ;
  p[INVSCALE] = 1.0/(p[SCALE]);
  if (GENERAL_PRINTLEVEL>=4) {
    PRINTF(" nn =%d\n", ((TBM_storage*)key->S)->nn);
    PRINTF("TBM invscale= %f; %f %f %f %f \n\n",
	   p[INVSCALE],linesimustep,key->param[SCALE],
	   linesimufactor,key->x[0][XSTEP]);
  }
  // nn=number of pixels to be simulated on the lines [==(length of diagonal 
  // through grid that is to be simulated) * (linesimufactor=="precision")]
  // 5 is for safety

  // maybe also loop and run through preference lists ...
  switch (key->tbm_method) {
  case CircEmbed :  
    ((TBM_storage*)key->S)->method = do_circulantembedding1dim; 
    error=init_circulantembedding1dim(key,p); 
    break;
  default : error=ERRORNOTPROGRAMMED; 
  }
  if (!error) return 0;
 ErrorHandling:
  if (key->destruct!=NULL) key->destruct(&key->S);
  key->destruct=NULL;
  key->active = false;
  return error;
}



void do_turningbands(key_type *key,Real *res END_WITH_GSLRNG)
{ 
  Real phi,InvSqrtNlines;
  Real *simuline;
  long n,nn,nnhalf; 
  simu1dim simulate;

#ifdef RF_GSL
  assert(RANDOM!=NULL);
#endif
  assert(key->active);

  simulate = ((TBM_storage*)key->S)->method;
  nn = ((TBM_storage*)key->S)->nn;
  //  printf("nnn = %d\n",nn);

  nnhalf = nn/2; 
  simuline = ((TBM_storage*)key->S)->simuline; 
  
  {
    register long i,endfor;
    endfor = key->totallength;
    for (i=0;i<endfor;i++) { res[i]=0.0; }
  }
 
  if (key->grid) {
    Real deltax, yoffset,deltay, nxhalf, nyhalf,xoffset,
      linesimuscaleX,linesimuscaleY,linesimuscaleZ;  
    int nx,zaehler,ny;
    nxhalf=((Real)(key->length[0]-1))/2.0; 
    nyhalf=((Real)(key->length[1]-1))/2.0;
    linesimuscaleX = ((TBM_storage*)key->S)->linesimuscale * key->x[0][XSTEP];
    linesimuscaleY = ((TBM_storage*)key->S)->linesimuscale * key->x[1][XSTEP];
    if (key->method==TBM2) { 
      Real deltaphi;
      InvSqrtNlines=1.0 / sqrt((Real) TBM2_LINES);
      deltaphi = PI / (Real) TBM2_LINES;
      phi = deltaphi * UNIFORM_RANDOM;        
      for (n=0;n<TBM2_LINES;n++){
	deltax=sin(phi) * linesimuscaleX;
	deltay=cos(phi) * linesimuscaleY;
	simulate(key,simuline END_WITH_RANDOM); /* depending on deltax and deltay,
						   shorter lines might be 
						   simulated! -> improvement in 
						   speed?! */
	/* then nnhalf must be kept variable, too: */
	yoffset= (double)nnhalf - nyhalf * deltay - nxhalf * deltax;
	zaehler = 0;
	for (ny=0;ny<key->length[1]; ny++) {	  
	  xoffset = yoffset;
	  for(nx=0;nx<key->length[0];  nx++) {
	    assert((xoffset<nn) && (xoffset>=0) );
	    res[zaehler++] += simuline[(long) xoffset];
	    xoffset += deltax;
	  }
	  yoffset += deltay;
	}
	phi += deltaphi;
      }
    } else {
      assert(key->method==TBM3);
      if (key->dim==2) { // TBM3
	InvSqrtNlines=1.0 / sqrt((Real) TBM3D2_LINES);
	for (n=0;n<TBM3D2_LINES;n++){
	  if ((GENERAL_PRINTLEVEL>=3) && (n % 100 == 0)) PRINTF("%d \n",n);
	  deltax= UNIFORM_RANDOM;// see Martin's tech rep for details
	  deltay= 
	    sqrt(1.0-deltax*deltax) * sin(UNIFORM_RANDOM*TWOPI) * linesimuscaleY;
	  deltax *= linesimuscaleX;
	  simulate(key,simuline END_WITH_RANDOM); 
	  /* depending on deltax and deltay, shorter lines might be 
	     simulated! -> improvement in speed!
	  */
	  yoffset=  (double)nnhalf - nyhalf * deltay - nxhalf * deltax ;
	  zaehler = 0;
	  for (ny=0;ny<key->length[1];ny++) {
	    xoffset = yoffset;
	    for(nx=0;nx<key->length[0];nx++) {
	      assert((xoffset<nn) && (xoffset>=0) );
	      res[zaehler++] += simuline[(long) xoffset];	  
	      if (!(fabs(simuline[(long) xoffset])<10000.0)) 
		PRINTF("s:%f ",simuline[(long) xoffset]);
	      xoffset += deltax;
	    }
	    yoffset += deltay;
	  }
	}  
      } else {
	Real dummy,nzhalf,zoffset,deltaz;
	int nz;
	assert(key->dim==3);      
	linesimuscaleZ= ((TBM_storage*)key->S)->linesimuscale * key->x[2][XSTEP];
	nzhalf=((Real)(key->length[2]-1))/2.0;
	InvSqrtNlines=1.0 / sqrt((Real) TBM3D3_LINES);
	for (n=0;n<TBM3D3_LINES;n++){
	  if ((GENERAL_PRINTLEVEL>=3) && (n % 100 == 0)) PRINTF("%d \n",n);
	  dummy = UNIFORM_RANDOM * PI;
	  deltaz = sin(dummy) * linesimuscaleZ;
	  dummy = cos(dummy);
	  deltay= UNIFORM_RANDOM * TWOPI;
	  deltax= cos(deltay) * dummy * linesimuscaleX;
	  deltay= sin(deltay) * dummy * linesimuscaleY;
	  simulate(key,simuline END_WITH_RANDOM);
	  zoffset= 
	    (double)nnhalf - nzhalf * deltaz - nyhalf * deltay - nxhalf * deltax; 
          zaehler = 0;    
	  for (nz=0;nz<key->length[2];nz++) {
	    yoffset = zoffset;
	    for (ny=0;ny<key->length[1];ny++) {
	      xoffset = yoffset;
	      for(nx=0;nx<key->length[0];nx++) {
		///
		if ((xoffset>=nn) || (xoffset<0))  
		 PRINTF(" %d %d %d | %d %d %d \n %f %f %f\n %d %f %f %f\n %d %f\n",
			 nx,ny,nz,key->length[0],key->length[1],key->length[2],
			 deltax,deltay,deltaz,
			 nnhalf,nxhalf,nyhalf,nzhalf,
                         nn,xoffset);
		///
		assert( (xoffset<nn) && (xoffset>=0) );
		res[zaehler++] += simuline[(long) xoffset];	  
		if (!(fabs(simuline[(long) xoffset])<10000.0)) 
		  PRINTF("s:%f ",simuline[(long) xoffset]);
		xoffset += deltax;
	      }
	      yoffset += deltay;
	    }
	    zoffset += deltaz;
	  }
	}  
      }
    }
  } else { // not grid
    Real halfx,halfy,ex,ey,offset,linescale;
    halfx = ((TBM_storage*)key->S)->half[0];
    halfy = ((TBM_storage*)key->S)->half[1];
    linescale = ((TBM_storage*)key->S)->linesimuscale;
    int i;
    if (key->method==TBM2) { 
      Real deltaphi;
      InvSqrtNlines=1.0 / sqrt((Real) TBM2_LINES);
      deltaphi = PI / (Real) TBM2_LINES;
      phi = deltaphi * UNIFORM_RANDOM;        
      for (n=0;n<TBM2_LINES;n++){
	ex=sin(phi) * linescale;
	ey=cos(phi) * linescale;
	simulate(key,simuline END_WITH_RANDOM); 
	/* depending on deltax and deltay, shorter lines might be simulated!
	   -> improvement in speed! */
	offset= nnhalf - halfy * ey - halfx * ex;  
	for (i=0;i<key->totallength;i++) {
	  register long index;
	  index = (long) (offset + key->x[0][i] * ex + key->x[1][i] * ey);
	  assert((index<nn) &&  (index>=0)); 
	  res[i] += simuline[index];
	}
	phi += deltaphi;
      }
    } else { // TBM3 
      assert(key->method==TBM3);
      if (key->dim==2) { 
	Real linesimusq;
	linesimusq = linescale * linescale;
	InvSqrtNlines=1.0 / sqrt((Real) TBM3D2_LINES);
	for (n=0;n<TBM3D2_LINES;n++){
	  if ((GENERAL_PRINTLEVEL>=3) && (n % 100 == 0)) PRINTF("%d \n",n);
	  ex= UNIFORM_RANDOM * linescale;// see Martin's tech rep for details
	  ey= sqrt(linesimusq - ex * ex) * sin(UNIFORM_RANDOM * TWOPI);
	  simulate(key,simuline END_WITH_RANDOM); 
	  /* depending on ex and ey, shorter lines might be simulated! -> 
	     improvement in speed! */
	  offset= (double) nnhalf - halfy * ey - halfx * ex; 
	  for (i=0;i<key->totallength;i++) {
	    register long index;
	    index = (long) (offset + key->x[0][i] * ex + key->x[1][i] * ey);
	    assert((index<nn) &&  (index>=0)); 
	    res[i] += simuline[index];
	  }	  
	}  
      } else { // TBM3, dim=3
	Real dummy,halfz,ez;
	assert(key->dim==3);      
	halfz = ((TBM_storage*)key->S)->half[2];
	InvSqrtNlines=1.0 / sqrt((Real) TBM3D3_LINES);
	for (n=0;n<TBM3D3_LINES;n++){
	  if ((GENERAL_PRINTLEVEL>=3) && (n % 100 == 0)) PRINTF("%d \n",n);
	  dummy = UNIFORM_RANDOM * PI;
	  ez = sin(dummy) * linescale;
	  dummy = cos(dummy) * linescale;
	  ey= UNIFORM_RANDOM * TWOPI;
	  ex= cos(ey) * dummy;
	  ey= sin(ey) * dummy;
	  simulate(key,simuline END_WITH_RANDOM);
	  offset= nnhalf - halfz * ez  - halfy * ey - halfx * ex;
	  for (i=0;i<key->totallength;i++) {
	    register long index;
	    index = 
	      (long) (offset + key->x[0][i] * ex + key->x[1][i] * ey + key->x[2][i] * ez);
	    assert((index<nn) &&  (index>=0)); 
	    res[i] += simuline[index];
	  }  
	}
      }
    }
  }

  { 
    register long i;    
    for(i=0;i<key->totallength;i++) { res[i] = res[i] * InvSqrtNlines; }
  }

  if (!(key->active=GENERAL_STORING)) {
    assert(key->destruct!=NULL);
    key->destruct(&key->S);
    key->destruct=NULL;
  }
  return;
}









