
/*
 Authors
 Martin Schlather, martin.schlather@cu.lu

 library for unconditional simulation of stationary and isotropic random fields

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


/*

  PRINTING LEVELS
  ===============
  error messages: 1
  forcation : 2
  minor tracing information : 3--5
  large debugging information: >10

*/

/*
  calculate the transformation of the points only once and store the result in
  a register (indexed by ncov) of key, not of s->. Advantages: points have to
  calculated only once for several tries; if zonal anisotropy then specific 
  methods may be called;  
*/

#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include <sys/timeb.h>
#include <assert.h>
#include <string.h>
#include "RFsimu.h"
#include <unistd.h>

void DeleteKey(int *keyNr)
{  
  key_type *key;
 
  int d, m;
  if (GENERAL_PRINTLEVEL>=4) { 
    PRINTF("deleting stored parameters of register %d ...\n", *keyNr);
  }
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {return;}
  key = &(KEY[*keyNr]);
  if (key->x[0]!=NULL) free(key->x[0]);
  for (d=0; d<MAXDIM; d++) {key->x[d]=NULL;}
  for (m=0; m<MAXCOV; m++) {
    if (key->destruct[m]!=NULL) key->destruct[m](&(key->S[m]));
    key->destruct[m]=NULL;
  }
  key->TrendModus = -1;
  if (key->TrendFunction!=NULL) free(key->TrendFunction); 
  key->TrendFunction=NULL;
  if (key->LinearTrend!=NULL) free(key->LinearTrend); 
  key->LinearTrend=NULL;
  if (key->destructX!=NULL) key->destructX(&(key->SX));
  key->destructX=NULL;
  key->active=false;
}

void DeleteAllKeys() {
  int i;
  for (i=0;i<MAXKEYS;i++) {DeleteKey(&i);}
}


void printkey(key_type *key) {
  int i,j;  
  char op_sign[][2]={"+","*"};

  PRINTF("\ngrid=%d, active=%d, anisotropy=%d, tbm_method=%d lx=%d sp_dim=%d,\ntotalpoints=%d, distr=%s Time=%d [%f %f %f]\ntimespacedim=%d traditional=%d compatible=%d mean=%f ncov=%d\n",	
	 key->grid, key->active, key->anisotropy,
	 key->tbm_method, key->lx, key->spatialdim, key->totalpoints,
	 DISTRNAMES[key->distribution], 
	 key->Time, key->T[0], key->T[1], key->T[2], 
	 key->timespacedim, key->traditional, key->compatible,
	 key->mean, key->ncov);
  PRINTF("\n");
  for (i=0; i<key->ncov; i++) {
    PRINTF("%d: meth=%s, cov=%d (%s) ^S=%d ^destruct=%d\n left=%d unimeth=%d (%s) [%d], param=", i,
	   METHODNAMES[key->method[i]],
	   key->covnr[i], CovList[key->covnr[i]].name,
	   key->S[i],key->destruct[i], 
	   key->left[i], 
	   key->unimeth[i], METHODNAMES[key->unimeth[i]], key->unimeth_used[i]
	   );
    for (j=0; j<key->totalparam; j++) PRINTF("%f,",key->param[i][j]);
    if (i < key->ncov - 1) PRINTF("\n%s\n",op_sign[key->op[i]]); 
    else PRINTF("\n\n");
  }
  //
  for (i=0; i<MAXDIM; i++) 
  PRINTF("%d : len=%d \n", i, key->length[i]);

  PRINTF("^x=%d ^y=%d ^z=%d ^SX=%d ^destructX=%d \n", 
	 key->x[0], key->x[1], key->x[2], key->SX, key->destructX
	 );

}

void printKEY(int *keyNr) {
  key_type *key;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    PRINTF("\nkeynr out of range!\n");
    return;
  }
  key = &(KEY[*keyNr]);
  PRINTF("\nNr=%d", *keyNr);
  printkey(key);
}


void printAllKeys() {
  int i;
  for (i=0;i<MAXKEYS;i++) {printKEY(&i);PRINTF("\n");}
}

int Transform2NoGrid(key_type *key, Real* param, int TrueDim, 
		     int *start_param, Real **xx)
{ 
  /* 
     this function transforms the coordinates according to the anisotropy
     matrix given in param (usually param is the param the first anisotropy
     matrix of the covariance product -- all the other anisotropy matrices
     have to multiple of the first, and this parameter is put into the
     covariance function as scale parameter 
     
     
     input  : key, param, TrueDim (got from GetTrueDim)
              start_param (vector of places where the parameters start)
     output : 
         **xx
           isotropic   : grid : NULL
                         arbitrary : x-coordinates multiplied with INVSCALE
           anisotropic : dim-reduced, x-coordinates multiplied with aniso matrix.
                         The index for the  temporal component runs the
                         slowest, if there is any

  */
  int i,j,k,d,w;
  long endfor, total;
  int dimM1;
  Real t, *x;
  dimM1 = key->timespacedim - 1;
  assert(xx!=NULL);
  if (key->anisotropy) {
    // used: TrueDim, start_param
    total = TrueDim * key->totalpoints;
    // total number of coordinates to store

    if ((x = *xx = (Real*) malloc(sizeof(Real) * total))==0)
      return ERRORMEMORYALLOCATION;
    /* determine points explicitly */
    if (key->Time && !key->grid) { // time component, no grid
      endfor = key->length[key->spatialdim];
      k=0;
      for (j=0, t=key->T[XSTART]; j<endfor; j++, t+=key->T[XSTEP])
	for (i=0; i<key->length[0]; i++)
	  for (d=0; d<TrueDim; d++, k++) {
 	    x[k] = 0.0;
	    for(w=0; w<key->spatialdim; w++)
	      x[k] += param[start_param[d]+w] * key->x[w][i];
	    x[k] += param[start_param[d]+w] * t;
	  }
    } else { 
      if (key->grid) {/* grid; with or without time component */
	Real y[MAXDIM]; /* current point within grid, but without
			   anisotropy transformation */
	int yi[MAXDIM]; /* counter for the current position in the grid */
	for (w=0; w<key->timespacedim; w++) {y[w]=key->x[w][XSTART]; yi[w]=0;}
	for (k=0; k<total; ){
	  for (d=0; d<TrueDim; d++, k++) {
	    x[k] = 0.0;
	    for(w=0; w<key->timespacedim; w++)
	      x[k] += param[start_param[d]+w] * y[w];
	  }
	  i = 0;
	  (yi[i])++;
	  y[i] += key->x[i][XSTEP];
	  while(yi[i]>=key->length[i]) {
	    yi[i]=0;
	    y[i] = key->x[i][XSTART];
	    if (i<dimM1) {
	      i++;
	      (yi[i])++;
	      y[i] += key->x[i][XSTEP];
	    } else {
	      assert(k==total);
	    }
	  }
	}
      } else {/* no grid, no time component */
	k = 0;
	for (i=0; i<key->length[0]; i++)
	  for (d=0; d<TrueDim; d++, k++) {
	    x[k] = 0.0;
	    for(w=0; w<key->spatialdim; w++)
	      x[k] += param[start_param[d]+w] * key->x[w][i];
	  }
      }
    }
  } else { /* isotropic, no time component */
    assert(fabs(param[INVSCALE] * param[SCALE] - 1.0)<EPSILON);
    if (key->grid) {    
      *xx = NULL;
    } else{
      if ((x = *xx = (Real*) malloc(sizeof(Real) * TrueDim * key->totalpoints)
	   )==0)
	// *xx is returned !!
	return ERRORMEMORYALLOCATION;      
      for (k=i=0; i<key->length[0]; i++) {
	for (j=0; j<key->spatialdim; j++) {
	  x[k++] = key->x[j][i] * param[INVSCALE];
	}
      }
    }
  }
  return 0;
}


int InternalGetGridSize(Real *x[MAXDIM], int *dim, int *lx)
{  
  int d;
  for (d=0; d<*dim; d++) {
    if ((x[d][XEND]<=x[d][XSTART]) || (x[d][XSTEP]<=0))
      return ERRORCOORDINATES;
    lx[d] = 1 + (int)((x[d][XEND]-x[d][XSTART])/x[d][XSTEP]); 
  }
  return NOERROR;
}


// obsolete?!
void GetGridSize(Real *x,Real *y,Real *z,int *dim,int *lx,int *ly,int *lz)
{ 
  // should now be consistent with `seq' of R
  // if not, please report!!!!!!!!!!!! 
  Real *xx[MAXDIM];
  int lxx[MAXDIM];
  xx[0]=x; xx[1]=y; xx[2]=z;
  if (InternalGetGridSize(xx,dim,lxx)) {
    *lx = *ly = *lz = 0;
  } else {
    *lx = lxx[0];  *ly = lxx[1];  *lz=lxx[2];
  }
}

#define COV_VARIO(FORMULA,FCT)\
  int v, aniso, k, j, endfor;\
  Real  var,  zz, zw, result, d;\
  cov_fct *cov=NULL; /* no meaning -- just avoids warning from -Wall */ \
  zw = result = 0.0;\
  if (anisotropy) {\
    Real z[MAXDIM];\
    int cov_isotropy, ncovM1;\
    ncovM1 = ncov - 1;\
    for (v=0; v<ncov; v++) {\
      zw = var = 1.0;\
      for(; v<ncov; v++) {\
        cov = &(CovList[covnr[v]]);\
        cov_isotropy=cov->isotropic;\
	var *= param[v][VARIANCE];\
	for (k=0; k<dim; k++) z[k]=0.0;\
	for(aniso=ANISO, k=0; k<dim; k++){\
	  for (j=0; j<dim; j++, aniso++) {\
	    z[k] += param[v][aniso] * x[j];\
            }\
        }\
 	if (cov_isotropy==ANISOTROPIC) zw *= cov->FCT(z,param[v],dim);\
          /* derzeit is dim in cov->FCT bis auf assert-checks unbenutzt */\
	else {\
          /* the following definiton allows for calculating the covariance */\
          /* funcion for spatialdim=0 and SPACEISOTROPIC model, namely as a */\
          /* purely temporal model; however, in CheckAndRescale such kind of */\
          /* calculations are prevend as being TRIVIAL(generally, independent */\
          /* of being SPACEISOTROPIC or not) -- unclear whether relaxation in */\
          /* CheckAndRescale possible */\
	  zz = 0.0;\
	  endfor = dim - 1;\
	  for (j=0; j<endfor; j++) zz += z[j] * z[j];\
	  if (cov_isotropy==FULLISOTROPIC) zz += z[endfor] * z[endfor];\
	  else z[1]=z[endfor];\
	  z[0] = sqrt(zz);\
	  zw *= cov->FCT(z, param[v], 1 + (int) (cov_isotropy==SPACEISOTROPIC));\
	}\
	if ( (v<ncovM1) && (op[v]==0)) break;\
      }\
      result += FORMULA;\
    }\
  } else {\
    if (dim==1) d=x[0];\
    else {\
      for (d=0.0, j=0; j<dim; j++) {d += x[j] * x[j];}\
      d = sqrt(d);\
    }\
    for (v=0; v<ncov; v++) {\
      zw = var = 1.0;\
      for(; v<ncov; v++) {\
	zz = d * param[v][INVSCALE];\
	var *= param[v][VARIANCE];\
        cov = &(CovList[covnr[v]]); /* cov needed in FORMULA for Vario */\
	zw *= cov->FCT(&zz,param[v],1);\
	if (op[v]==0) break; /* v=ncov does not matter since stopped anyway */\
      }\
      result += FORMULA;\
    }\
  }\
 return result;


Real CovFct(Real *x, int dim, int *covnr, int *op,
  param_type param, int ncov, bool anisotropy) {
  COV_VARIO(var * zw, cov);
}

Real LocalCovFct(Real *x, int dim, int *covnr, int *op,
  param_type param, int ncov, bool anisotropy) {
  COV_VARIO(var * zw, cov_loc);
}

Real Vario(Real *x, int dim, int *covnr, int *op,
	    param_type param, int ncov, bool anisotropy){
  // Vario defines some kind of multiplication (namely as the mulitplication
  // of the respective covariance functions)
  // 
  // make sure that no genuine variogram model are used in multiplication terms
  //
  // first part says how to add models
  // cov.variogram checks whether last model in the block is a variogram
  // this may allow relaxations towards multiplications of varioagrams in
  // the sense that the corresponding covariance functions are multiplied
  // and then 1-block is used as value for variogram (only covariance functions
  // may appear in this case
  COV_VARIO(cov->variogram ? var * zw : var * (1.0 - zw), cov); // 30.1.02
}


// checks a *single* basic covariance function
// p and param must be of the same kind (Real)!!! 
int CheckAndRescale( int covnr, int naturalscaling, Real *p,
		    int timespacedim, int anisotropy, Real* param) 
{ 
  // IMPORTANT!: p is the parameter vector of the R interface !!
  // output param : equals p, except naturalscaling
  int actparam, error, endfor=0 /*no meaning, avoids -Wall warning*/, i, j;
  Real newscale;

  if ((covnr>=currentNrCov) || (covnr<0)) return ERRORNOTDEFINED;
  if (p[VARIANCE]<0.0) return ERRORNEGATIVEVAR;
  
  actparam = KAPPA + CovList[covnr].kappas;
  if (anisotropy) {
    endfor = ANISO + timespacedim * timespacedim;
    for (j=actparam, i=ANISO; i<endfor; i++, j++) param[i] = p[j];     
  } else {
    if (CovList[covnr].isotropic!=FULLISOTROPIC) return ERRORANISOTROPIC;
    param[SCALE] = p[actparam];	
  }
  
  { // is dimension OK ? 
    int TrueDim, start_param[MAXDIM], index_dim[MAXDIM];
    bool no_last_comp;
    switch (CovList[covnr].isotropic) {
    case FULLISOTROPIC : 
      GetTrueDim(anisotropy, timespacedim, param,  &TrueDim, &no_last_comp, 
		 start_param, index_dim);
      break;
    case SPACEISOTROPIC :
      if (anisotropy) {
	GetTrueDim(anisotropy, timespacedim, param,  &TrueDim, &no_last_comp, 
		 start_param, index_dim);
	if (!no_last_comp) TrueDim--;
      } else {
	return ERRORNOTANISO;
      }
      break;
    case ANISOTROPIC :
      if (!anisotropy) return ERRORNOTANISO;
      TrueDim = timespacedim;
      break;
    default : assert(false);
    } 
    if (TrueDim==0) return ERRORTRIVIAL;
    if (CovList[covnr].method(TrueDim, false)==Forbidden)
      return ERRORCOVNOTALLOWED;
  }

  // specific parameter check for the covariance function
  if ((CovList[covnr].check!=NULL) && 
      ((error=(CovList[covnr].check(p, timespacedim, Nothing))))!=NOERROR)
    return error;
 
  // getnaturalscaling
  GetNaturalScaling(&covnr,&(p[KAPPA]),&naturalscaling,&newscale,&error);
  if (error!=NOERROR) return error;

  if (anisotropy) {
    newscale = 1.0 / newscale;
    for (j=actparam, i=ANISO; i<endfor; i++, j++) param[i] = p[j] * newscale;
  } else {
    param[SCALE] *= newscale;	
    if (CovList[covnr].cov==nugget && param[SCALE]==0) param[SCALE] = 1;
    param[INVSCALE] = 1.0 / param[SCALE];
    if (param[SCALE] <= 0.0) return ERRORNEGATIVESCALE;
  }
  for (i=0; i<actparam; i++) param[i] = p[i]; /* variance, kappas */

  return NOERROR;
}


#define FIRST_CHECK_COV_VARIO(ERROR)\
  if (currentNrCov==-1) InitModelList(); assert(CovList!=NULL);\
  if ((*logicaldim!=*xdim) && (*anisotropy || (*xdim!=1))) \
    {ERROR=ERRORDIMMISMATCH; goto ErrorHandling;}\
  if ((*xdim<1) || (*xdim>MAXDIM)) {ERROR=ERRORDIM; goto ErrorHandling;}\
  if (*anisotropy) { \
    pAdd = *xdim * *xdim + 1; /* VARIANCE */\
  } else {\
    pAdd = 2;\
  }\
  for (i=pAdd * *ncov, v=0; v<*ncov; v++) {\
    if (CovList[covnr[v]].cov==NULL) {\
      ERROR=ERRORNOTPROGRAMMED; goto ErrorHandling;}\
    i += CovList[covnr[v]].kappas;\
    }\
  if (i!=*np) {ERROR=ERRORPARAMNUMBER; goto ErrorHandling;}

// note: number of parameters is not checked whether it is correct!!
// done by R function `PrepareModel'
void CovarianceNatSc(Real *x, int *lx, int *covnr, Real *p, int *np, 
		     int *logicaldim, /* timespacedim ! */
		     int *xdim,
		     int *ncov, int *anisotropy, int *op,
		     Real *result, int *naturalscaling)
{
  // note: column of x is point !!
  int error, i,pAdd, v;
  param_type param;

  error = NOERROR;
  FIRST_CHECK_COV_VARIO(error);  

  for (v=0; v<*ncov; v++) {
    if (((error=(CheckAndRescale(covnr[v], *naturalscaling, p, *logicaldim,
				*anisotropy, param[v])))!=NOERROR) ||
	CovList[covnr[v]].variogram) 
      goto ErrorHandling;    
    // increase of address !
    p += pAdd+CovList[covnr[v]].kappas; // pointer addition!
  }
  
  for (i=0; i<*lx; i++, x += *xdim) {
    result[i] = CovFct(x, *xdim, covnr, op, param, *ncov, *anisotropy);
  }
  return;
  
 ErrorHandling:
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,error);
  for (i=0;i<*lx;i++) { result[i] = RF_NAN; }
  return;	
}

void Covariance(Real *x,int *lx, int *covnr, Real *p, int *np, 
		int *logicaldim, /* timespacedim ! */
		int *xdim,
		int *ncov, int *anisotropy, int *op,
		Real *result)
{
  CovarianceNatSc(x,lx,covnr,p,np,logicaldim,xdim,
		  ncov, anisotropy, op,result, &GENERAL_NATURALSCALING);
}

static param_type Unchecked_param;
static int Unchecked_covnr[MAXCOV], Unchecked_DIM, Unchecked_op[MAXDIM],
  Unchecked_ncov;
static bool Unchecked_anisotropy;

void InitUncheckedCovFct(int *covnr,Real *p,int *np, 
			 int *logicaldim, /* timespacedim ! */
			 int *xdim,
			 int *ncov, int *anisotropy, int *op,
			 int *naturalscaling, int *error){
  int i, pAdd, v;
 
  FIRST_CHECK_COV_VARIO(*error);
  for (v=0; v<*ncov; v++) {
    if ((*error=(CheckAndRescale(covnr[v], *naturalscaling, p, *logicaldim,
				 *anisotropy, Unchecked_param[v]))) || 
	CovList[covnr[v]].variogram) 
      goto ErrorHandling;
    p += pAdd+CovList[covnr[v]].kappas;
  }
  memcpy(Unchecked_covnr, covnr, sizeof(int) * *ncov);
  memcpy(Unchecked_op, op, sizeof(int) * *ncov);
  Unchecked_ncov = *ncov;
  Unchecked_anisotropy = (bool) *anisotropy;
  Unchecked_DIM = *xdim;
  return;
  
 ErrorHandling:
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,*error);
  return;	
}

void UncheckedCovFct(Real *x, int *lx, Real *result)
{
  int i;
  /*
     WARNING:  make sure that CovList is initialised,
                              covnr is correct,
                              p is fine according to covnr (length(p), value)
                              PracticalRange is included in p[SCALE] if 
			      necessary
  */
  for (i=0; i<*lx; i++, x += Unchecked_DIM) {
    result[i] = CovFct(x, Unchecked_DIM, Unchecked_covnr, 
		       Unchecked_op, Unchecked_param, 
		       Unchecked_ncov, Unchecked_anisotropy);
  }
}

void VariogramNatSc(Real *x,int *lx,int *covnr,Real *p,int *np,
		    int *logicaldim, /* timespacedim ! */
		    int *xdim,
		    int *ncov, int *anisotropy, int *op,
		    Real *result, int *naturalscaling)
{
  // note: column of x is point !!
  int error=NOERROR, i, pAdd, v;
  param_type param;

  FIRST_CHECK_COV_VARIO(error);

  for (v=0; v<*ncov; v++) {
    if ((error=CheckAndRescale(covnr[v], *naturalscaling, p, *logicaldim, 
			      *anisotropy,
			      param[v]))!=NOERROR) goto ErrorHandling;
    if ( (v>0) && (op[v-1]) && 
	 (CovList[covnr[v]].variogram || CovList[covnr[v-1]].variogram)) {
      // only covariance functions may be multiplied!
      error=ERRORNOMULTIPLICATION; goto ErrorHandling;
    } 
    p += pAdd+CovList[covnr[v]].kappas;
  }
  for (i=0; i<*lx; i++, x += *xdim) {
    result[i] = Vario(x, *xdim, covnr, op, param, *ncov, *anisotropy);
  }
  return;
   
 ErrorHandling:
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,error);
  for (i=0;i<*lx;i++) { result[i] = RF_NAN; }
  return;	 
}

void Variogram(Real *x,int *lx,int *covnr,Real *p,int *np, 
	       int *logicaldim, /* timespacedim ! */
	       int *xdim,
	       int *ncov, int *anisotropy, int *op,
	       Real *result)
{
  VariogramNatSc(x,lx,covnr,p,np,logicaldim,xdim,ncov,anisotropy,op,
		 result,&GENERAL_NATURALSCALING);
}
  
void CovarianceMatrixNatSc(Real *dist,int *lx, int *covnr,Real *p, int *np, 
			   int *logicaldim, /* timespacedim ! */
			   int *xdim,
			   int *ncov, int *anisotropy, int *op,
			   Real *result, int *naturalscaling)
{
  // note: column of x is point !!
  int error=NOERROR, i, ii, pAdd, v, endfor, lxP1;
  long lxq, ve, ho;
  param_type param;
  Real var;
 
 
  FIRST_CHECK_COV_VARIO(error);  
  for (v=0; v<*ncov; v++) {
    if ( ((error=(CheckAndRescale(covnr[v], *naturalscaling, p, *logicaldim, 
				 *anisotropy, param[v]))) != NOERROR) || 
	 CovList[covnr[v]].variogram) 
      goto ErrorHandling;
    p += pAdd + CovList[covnr[v]].kappas;  // pointer increment!
  }
   
  var = CovFct(ZERO, *xdim, covnr, op, param, *ncov, *anisotropy);
  lxP1 = *lx + 1;
  for (ii=*lx, i=0; ii>0; i+=lxP1, ii--) {
    result[i] = var;
    endfor = i + ii;
    for (ve = i + 1, ho = i + *lx; ve < endfor; ve++, ho += *lx, dist += *xdim){
      result[ve] = result[ho] = 
	CovFct(dist, *xdim, covnr, op, param, *ncov, *anisotropy);
    }
  }
  return;
  
 ErrorHandling:
  lxq = *lx * *lx;
  for(i=0; i<lxq; i++) result[i] = RF_NAN;
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,error); 
  return;	
}

void CovarianceMatrix(Real *dist, int *lx,int *covnr, Real *p, int *np,
		      int *logicaldim, /* timespacedim ! */
		      int *xdim,
		      int *ncov, int *anisotropy, int *op,
		      Real *result)
{
  CovarianceMatrixNatSc(dist, lx, covnr, p, np, logicaldim, xdim, ncov, 
			anisotropy, op,	result, &GENERAL_NATURALSCALING);
}


void VariogramMatrix(Real *dist,int *lx,int *covnr,Real *p,int *np, Real *result)
{ 
  // not programmed yet
  assert(false);
}

void insert_method(int* last_incompatible, key_type *key,
			 SimulationType m) {
//  printf("insert method\n");
  int idx;
  assert(key->S[key->n_unimeth]==NULL);
  assert(key->destruct[key->n_unimeth]==NULL);
  if (do_incompatiblemethod[m]!=NULL) {
    // first the non-compatible (i.e. which do not allow direct 
    // adding to the random field), then those which allow for 
    // direct adding
    // --if attention wouldn't be payed to the ordering of the methods,
    // the simulation algorithm had to create an intermediate field 
    // more frequently -- sometimes the fields are pretty large,
    // and then this procedure is quite helpful
    (*last_incompatible)++;

    key->unimeth_used[key->n_unimeth]=key->unimeth_used[*last_incompatible];
    key->unimeth[key->n_unimeth]=key->unimeth[*last_incompatible];
    key->S[key->n_unimeth]=key->S[*last_incompatible]; // necessary
    // only if traditional and further methods are tried due to
    // partial failure
    key->destruct[key->n_unimeth]=key->destruct[*last_incompatible];

    key->S[*last_incompatible] = NULL;
    key->destruct[*last_incompatible] = NULL;
    idx = *last_incompatible;
  } else idx=key->n_unimeth;
  key->unimeth[idx] = m;
  key->unimeth_used[idx] = false;
  key->n_unimeth++;
//  printf("end insert method\n");
}

void delete_method(int *M, int *last_incompatible, key_type *key) {
  // note that key->n_unimeth had been decreased before this function
  // is called
  int idx;

//  printf("delete method\n");
  if (key->destruct[*M]!=NULL) key->destruct[*M](&(key->S[*M]));
  key->destruct[*M]=NULL;
  if (do_incompatiblemethod[key->unimeth[*M]]!=NULL) {
    assert(*last_incompatible>=0);
    assert(*M<=*last_incompatible);
    key->unimeth[*M]=key->unimeth[*last_incompatible];
    key->unimeth_used[*M] = key->unimeth_used[*last_incompatible];
    key->S[*M] = key->S[*last_incompatible];
    key->S[*last_incompatible]=NULL;
    key->destruct[*M]=key->destruct[*last_incompatible];
    key->destruct[*last_incompatible]=NULL;
    idx = *last_incompatible;
    (*last_incompatible)--;
  } else idx = *M;
//  printf("idx=%d %d %d %d\n", idx, *M ,*last_incompatible, key->n_unimeth);
  key->n_unimeth--;
  if (idx==key->n_unimeth){ // deleting of the very last entry
    assert(*M==key->n_unimeth); 
  } else {
    assert(idx>=0 && idx<key->n_unimeth);
    key->unimeth[idx] = key->unimeth[key->n_unimeth];
    key->unimeth_used[idx] = key->unimeth_used[key->n_unimeth];
    key->S[idx] = key->S[key->n_unimeth];
    key->S[key->n_unimeth]=NULL;
    key->destruct[idx] = key->destruct[key->n_unimeth];
    key->destruct[key->n_unimeth]=NULL;
  }
  key->unimeth_used[key->n_unimeth] = false;
  (*M)--;
//  printf("end delete method\n");
}

void InitSimulateRF(Real *x, Real *T, 
		    int *spatialdim, /* spacial dim only ! */
		    int *lx, int *grid, 
		    int *Time,
	            int *covnr, Real *ParamList, int *nParam,
		    Real *mean,
		    int *ncov, int *anisotropy, int *op,
		    int *method, 
		    int *distr, /* still unused */
		    int *keyNr , int *error)
{ // grid  : boolean
  // keyNr : label where intermediate results are to be stored;
  //         can be chosen freely between 0 and MAXKEYS-1

  // NOTE: if grid and Time then the Time component is treated as if
  //       it was an additional space component
  //       On the other hand, SPACEISOTROPIC covariance function will
  //       always assume that the second calling component is a time component!

  bool user_defined,finished;
  int i,act_number, last_incompatible, extndd_nr, d, endfor, 
    pAdd, v, m, M;
  unsigned long totalBytes;
  Real p[TOTAL_PARAM], *PL;
  bool method_used[(int) Nothing];
  // preference lists, distinguished by grid==true/false and dimension
  // lists must end with Nothing!
  SimulationType Merr=Forbidden;
  SimulationType pg1[]={CircEmbed, CircEmbedLocal, Direct, AdditiveMpp, Nothing}; 
  SimulationType pg2[]={CircEmbed, CircEmbedLocal, TBM2, SpectralTBM,  
			 AdditiveMpp, TBM3, Direct, Nothing};
  SimulationType pg3[]={CircEmbed, CircEmbedLocal, TBM3, AdditiveMpp, 
			Direct, Nothing};
  SimulationType png1[]={Direct, AdditiveMpp, Nothing};
  SimulationType png2[]={TBM2, SpectralTBM, AdditiveMpp, TBM3, Direct, 
			 Nothing};
  SimulationType png3[]={TBM3, AdditiveMpp, Direct, Nothing};
  SimulationType *preference_list,first_method[MAXCOV];


  key_type *key;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORREGISTER; goto ErrorHandling;}
  key = &(KEY[*keyNr]);
  key->storing = GENERAL_STORING;
 
  // check parameters

  key->timespacedim = *spatialdim + (int) (bool) (*Time);
  if (currentNrCov==-1) InitModelList(); assert(CovList!=NULL); 
  if ((*spatialdim<1) || (key->timespacedim>MAXDIM) || (*Time>1) || (*Time<0)){
    *error=ERRORDIM; goto ErrorHandling;}    
  if (x==NULL) {*error=ERRORCOORDINATES;goto ErrorHandling;}
  for (i=0; i<*ncov; i++) 
    if ((covnr[i]<0) | (covnr[i]>=currentNrCov)) {
      *error=ERRORCOVNROUTOFRANGE; goto ErrorHandling; }
  if (*anisotropy) {
   pAdd = 1 + key->timespacedim * key->timespacedim;
   // not literally "total"param; there might be gaps in the sequence of the 
   // parameters (if not all seven KAPPA parameters are used)
   key->totalparam = ANISO + key->timespacedim * key->timespacedim; 
  }
  else {
    pAdd = 2;
    key->totalparam = SCALE + 2; // SCALE, INVSCALE 
  }
  
  // also checked by R function `PrepareModel'
  for (i=pAdd * *ncov, v=0; v<*ncov; v++) {
    i += CovList[covnr[v]].kappas;
  }
  if (i!=*nParam) {*error=ERRORPARAMNUMBER; goto ErrorHandling;}
  // bytes to be stored for each vector of coordinates:

  totalBytes =  sizeof(Real) * *lx * *spatialdim;
  if (key->active) {   
    assert(key->x[0]!=NULL);
    key->active = 
      (key->lx==*lx) && 
      (key->spatialdim==*spatialdim) && 
      (key->grid==(bool)*grid) && 
      (key->distribution==*distr) &&
      (key->ncov==*ncov) &&
      (!memcmp(key->covnr, covnr, sizeof(int) * *ncov)) &&
      (key->anisotropy==*anisotropy) &&
      (key->mean==*mean) &&
      ( (*Time == key->Time) &&
	(!(*Time) || ! memcmp(key->T, T, sizeof(Real) * 3))) &&
      (!memcmp(key->op, op, (*ncov>0) ? *ncov-1 : 0)) &&
      (!memcmp(key->x[0], x, totalBytes)) 
      ;       
  }

  if (!key->active) {
    // if not active, all the parameters have to be 
    //                        stored again
    // save all the parameters if key has not been active   
    key->lx=*lx;
    key->spatialdim = *spatialdim;
    key->grid=(bool) *grid;
    key->distribution=*distr; 
    key->ncov=*ncov;
    memcpy(key->covnr, covnr, sizeof(int) * *ncov);
    key->anisotropy = (bool) *anisotropy;    
    key->mean = *mean;
    // operators +, *
    memcpy(key->op, op, (*ncov>0) ? *ncov-1 : 0);
    // coordinates

    if (key->x[0]!=NULL) free(key->x[0]);
    if ((key->x[0]=(Real*) malloc(totalBytes))==NULL){
      *error=ERRORMEMORYALLOCATION;goto ErrorHandling;
    }
    memcpy(key->x[0], x, totalBytes);

    for (d=1; d<*spatialdim; d++) key->x[d]= &(key->x[0][d * *lx]);
    for (; d<MAXDIM; d++)  key->x[d]=NULL;
    
    //Time
    if (key->Time = *Time) {
      if (!key->anisotropy) { *error = ERRORTIMENOTANISO; goto ErrorHandling; }
      memcpy(key->T, T, sizeof(Real) * 3);
      if (key->grid) {
	key->x[key->spatialdim] = key->T;
      }
    }
  } else { // key->active
    // CHECK:
    for (d = key->timespacedim; d<MAXDIM; d++) {assert(key->x[d]==NULL);} 
  }
 
  // check the parameter lists -- note that naturalscaling changes the values !
  for (v=0, PL=ParamList; v<key->ncov; v++){
    int j;
    if ((*error=(CheckAndRescale(covnr[v], GENERAL_NATURALSCALING, PL, 
				 key->timespacedim, *anisotropy, p))) != NOERROR)
      goto ErrorHandling; 
    for (j=KAPPA+CovList[covnr[v]].kappas; j<=LASTKAPPA; j++) p[j]=0.0;
    if (!(key->active &=
	  (key->covnr[v]==covnr[v]) &&
	  (!memcmp(key->param[v], p, sizeof(Real) * key->totalparam))
	  )) {
      key->covnr[v]=covnr[v];
      memcpy(key->param[v], p, sizeof(Real) * key->totalparam);
    }
   PL+=pAdd+CovList[covnr[v]].kappas;
  }

  // method list must be treated separately since method is integer coded
  // when passed to InitSimulate !
  if (user_defined = (method[0]>=0)) {//either all or none shld be user defined
    for (v=0; v<key->ncov; v++) {
      if (method[v] >= (int) Nothing) {
	*error=ERRORNOTDEFINED; goto ErrorHandling;}
      if ((key->method[v]!=(SimulationType) method[v])) {	
	key->active = false;
	key->method[v]=(SimulationType) method[v];
      }
    }
  } else {
    // TO DO: not consistent -- unclear which method to choose
    if (!key->active) {
      for (v=0; v<key->ncov; v++) {
	key->method[v]=CovList[covnr[v]].method(key->spatialdim, key->grid);
      }	
    }
  }   
  
  if (key->active) { // detected that all relevant parameters agree, 
    //              no initialization is necessary
    if (GENERAL_PRINTLEVEL>=2) ErrorMessage(Nothing, USEOLDINIT);
    *error=NOERROR; 
   return; /* ***** END ***** */ 
  }

  // Delete stored, intermediate result if there is any
  for (v=0; v<MAXCOV; v++) {
    if (key->destruct[v]!=NULL) {
      key->destruct[v](&(key->S[v]));
      key->destruct[v] = NULL;
    }
  }
  if (key->destructX!=NULL) {
    key->destructX(&(key->SX));
    key->destructX = NULL;
  }
  
  // Check
  if (!key->anisotropy)
    for (v=0; v<*ncov; v++) {
      assert(fabs(key->param[v][SCALE] * key->param[v][INVSCALE]-1.0) < EPSILON);
   }

  // minor settings, usually overtaken from previous run, if all the 
  // other parameters are identical
  
  // no further choice of TBM_METHOD !!!!!!!!!!!!! TO BE CHANGED IN FUTURE
  // if (TBM_METHOD==Nothing) {key->tbm_method=CovList[*covnr].tbm_method;}
  //else 
  { key->tbm_method=TBM_METHOD;}
 
  if (key->grid) { 
    if ((*lx!=3) || 
	(InternalGetGridSize(key->x,&key->timespacedim,key->length))) {
      *error=ERRORCOORDINATES; goto ErrorHandling;
    }
    for (d=key->timespacedim; d<MAXDIM; d++) key->length[d] = 1;
    for (key->spatialtotalpoints=1, d=0; d<key->spatialdim; d++) {
      key->spatialtotalpoints *= key->length[d];
    }
    key->totalpoints = key->spatialtotalpoints;
    if (key->Time) key->totalpoints *= key->length[key->spatialdim];

    // preference list
    switch (key->timespacedim) {
    case 1 : preference_list = pg1; break;
    case 2 : preference_list = pg2; break;
    case 3: case 4 : preference_list = pg3; break;
    default : assert(false);
    }
  } else { // not grid
    key->totalpoints = key->spatialtotalpoints = key->length[0]= *lx;
   if (key->Time) {
      int Tdim; Real *Tx[MAXDIM];
      Tdim = 1;
      Tx[0] = key->T;
      if (InternalGetGridSize(Tx,&Tdim,&(key->length[key->spatialdim]))) {
	*error=ERRORCOORDINATES; goto ErrorHandling;
      } 
      key->totalpoints *= key->length[key->spatialdim];
    } 
    // just to say that considering these values does not make sense
    for (d=1; d<key->spatialdim; d++) key->length[d]=0; 
    // correct value for higher dimensions (never used)
    for (d=key->timespacedim; d<MAXDIM; d++) key->length[d] = 1;
    switch (key->timespacedim) {
    case 1 : preference_list = png1; break;
    case 2 : preference_list = png2; break;
    case 3 : case 4 : preference_list = png3; break;
    default : printkey(key); assert(false);
    }
  }

  // simple (traditional in the sence of version 1.0 of RandomFields) 
  // definition of the covariance function?
  key->traditional = (!key->anisotropy) && ((key->ncov==1) || 
    ((key->ncov==2) && ( (key->method[0]!=Nugget) &&
		     key->method[1]==Nugget) ));

  // check (and set) method in case of multiplicative covariance function
  // here, the internal first prefence might be overwritten
  endfor = key->ncov - 1;
  if (user_defined) {
    for (v=0; v<endfor; v++) {  
      if (op[v]) { // *
	if (method[v]!=method[v+1]) {*error=ERRORMETHODMIX; goto ErrorHandling;}
      }
    } 
  } else {
    for (v=0; v<endfor; v++) {  
      if (op[v]) { // *
	key->method[v] = key->method[v+1] = (key->grid) ? CircEmbed : TBM3; 
      }
    }
  }

  for (v=0; v<*ncov; v++) first_method[v] = key->method[v]; 
  // a first prefered method is tried, before the preference_list
  // is consulted; first_method ist storred to avoid that 
  // it is tried a second time when consulting the preference_list

  finished = false; 
  assert(*error==NOERROR); 
  extndd_nr = 0; // used if !traditional; has same role as act_number 
  //               in the traditional case

  for (v=0; v<*ncov; v++) key->left[v]=true;
  // indicates which covariance function are left to be initialized

  key->n_unimeth = 0;
  act_number = -1;
  last_incompatible = -1;

  while (!finished) {
    int err_occurred;

// printf("\nYY %d %d %d %d %d\n", key->n_unimeth, M, last_incompatible, user_defined, key->traditional); for (v=0; v<MAXCOV; v++) printf(" %d ", key->destruct[v]);

    *error = err_occurred = NOERROR;// necessary in case ncov=0, since then 
    //                          the following loop is not entered


    // create list of methods used
    for (i=0; i<(int) Nothing; i++) { method_used[i] = false; }
    for (v=0; v<key->ncov; v++) {
      if (key->left[v]) { 
      method_used[key->method[v]] = true;  
      // if user_defined take over the user_defined one
      // if !user_defined take specifically preferred one (for 
      //                     each covariance function)
      }
    }

    for (m=CircEmbed; m<Nothing; m++) 
      if (method_used[m])
	insert_method(&last_incompatible, key, (SimulationType) m); 
        // increments key->n_unimeth

    assert((key->destructX==NULL) && (key->SX==NULL));
    for (M=0; M<key->n_unimeth; M++) if (!key->unimeth_used[M]) {
      bool left[MAXCOV];
      key->unimeth_used[M] = true; // optimistic, corrected by delete_method,
      //                              if method fails
      for (v=0; v<key->ncov; v++) left[v] = key->left[v]; // sicherung,
      // da init_method ueberschreibt und *danach* entscheidet, ob
      // die methode funktioniert -- ueber key->method ist es zu schwierig
      // da eine methode mehrmals verwendet werden kann
      if (GENERAL_PRINTLEVEL>=4)
	PRINTF("Init %d %s %d %s\n",M,
	       METHODNAMES[key->unimeth[M]], key->n_unimeth, "--");

// printf("\nXX %d %d %d %d %d\n", key->n_unimeth, M, last_incompatible, user_defined, key->traditional); for (v=0; v<MAXCOV; v++) printf(" %d ", key->destruct[v])

      assert((key->destruct[M]==NULL));
      assert((key->S[M]==NULL));
      if (key->unimeth[M]>Special) {
	err_occurred=ERRORMETHODNOTALLOWED; goto ErrorHandling;}
      err_occurred = init_method[key->unimeth[M]](key, M);
      if (GENERAL_PRINTLEVEL>=5) {
	PRINTF("return code %d :", err_occurred);
	ErrorMessage(Nothing, err_occurred);
      }
      if (err_occurred) {
	if (err_occurred==NOERROR_REPEAT) {
	  insert_method(&last_incompatible, key, key->unimeth[M]);
	  err_occurred = NOERROR;
	  // incompatible here means that partial simulation results
	  // are not directly added to the final result by the algorithm
	  // E.g. TBM is incompatible since at the final a division is performed
	} else if (err_occurred==NOERROR_ENDOFLIST) {
	  delete_method(&M, &last_incompatible, key);
	  err_occurred = NOERROR;
	} else { 
	  Merr = key->unimeth[M];
	  *error = err_occurred; // if not traditional, error may occur anywhere
	  //                     and final *error might be 0

// printf("\n%d %d %d %d %d\n",key->n_unimeth,  M, last_incompatible, user_defined, key->traditional); for (v=0; v<MAXCOV; v++) printf(" %d ", key->destruct[v]);

	  if (user_defined || key->traditional) break;
	  delete_method(&M, &last_incompatible, key);

// printf("\nZZ %d %d %d %d %d\n",key->n_unimeth,  M, last_incompatible, user_defined, key->traditional); for (v=0; v<MAXCOV; v++) printf(" %d ",key->destruct[v]);

	  for (v=0; v<key->ncov; v++) key->left[v] = left[v];
	}
      }
    }
    key->compatible = last_incompatible <= 0; // if at most one incompatible
    // method, then no intermediate results must be stored

    if (!(finished=user_defined || *error==NOERROR)) {
      if (key->traditional) {
	// tabula rasa, since the (true) covariance or the nugget may
	// have failed; just to be sure...
	for (m=0; m<key->n_unimeth; m++) {
	  if (key->destruct[m]!=NULL) key->destruct[m](&(key->S[m]));
	  key->destruct[m]=NULL;
	}
	if (key->destructX!=NULL) key->destructX(&(key->SX));
	key->destructX=NULL;
	key->n_unimeth = 0;

	for (v=0; v<*ncov; v++) key->left[v]=true;
	// indicates which covariance function are left to be initialized

	// there is a first_method given by the programmer when defining
        // a cov model; this method is tried first; afterwards the 
        // standard list of methods is checked; however the first_method
	// should not be reconsidered!
	if (GENERAL_PRINTLEVEL>1) ErrorMessage(key->method[0], *error);
	act_number++; // within the preference list, act_number is initially -1
	if (preference_list[act_number]==first_method[0]) {
	  // the new(next or first) method out of the standard list might be 
	  // the first_method
	  act_number++;
	}
	key->method[0]=preference_list[act_number];
	// end of the list?
	finished=finished || (preference_list[act_number]==Nothing);
	//*error = ERRORFAILED;
      } else { // not traditional
	assert(key->S[M]==NULL);
	act_number++;
	for (v=0; v<key->ncov; v++) {
	  if (key->left[v]) key->method[v] = preference_list[act_number];
	}
	finished=finished || (preference_list[act_number]==Nothing);
      }
    } //if !finished || err_occurred
  } // while !finished

  key->active = !*error;
  if ((*error && (GENERAL_PRINTLEVEL>0)) ||
      (!user_defined && GENERAL_PRINTLEVEL>1)) { 
    if (!user_defined // && !key->traditional
	&& *error)
      PRINTF("algorithm failed. Last error has been: ");
    if (*error && GENERAL_PRINTLEVEL>0) {
      ErrorMessage(Merr, *error);
      PRINTF("\n");
    }
  }
  if (!key->active) goto ErrorHandling;

  if (GENERAL_PRINTLEVEL>=4) PRINTF("Successfully initialised.\n");
  return;
  
  ErrorHandling:
  DeleteKey(keyNr);  
}

  
void DoPoissonRF(int *keyNr, Real *res, int *error)
{
  // not programmed yet
  assert(false);
}

void DoSimulateRF(int *keyNr, Real *res, int *error)
{
  Real  *part_result;
  long i, m;
  key_type *key;

  part_result = NULL;
  key = &(KEY[*keyNr]);
  *error=NOERROR; 

  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  if (!key->active) {*error=ERRORNOTINITIALIZED; goto ErrorHandling;}

  if (key->n_unimeth==0) {
    for(i=0; i<key->totalpoints; i++) res[i] = key->mean;
    return;
  } 
  
  GetRNGstate();
  if (!key->compatible) {
    if ((part_result=(Real *) malloc(sizeof(Real) * key->totalpoints))
	==NULL) {*error=ERRORMEMORYALLOCATION; goto ErrorHandling;}
  }
  for(m=0; m<key->n_unimeth; m++) {
    if (GENERAL_PRINTLEVEL>=7)
      PRINTF("DoSimu %d %s\n",m, METHODNAMES[key->unimeth[m]]);
    if (do_compatiblemethod[key->unimeth[m]]!=NULL) {
      if (GENERAL_PRINTLEVEL>=7) PRINTF("compatible: add=%d\n",m!=0);
      do_compatiblemethod[key->unimeth[m]] (key, m!=0, m, res ); 
    } else {
      if (m==0) {
	if (GENERAL_PRINTLEVEL>=7) PRINTF("incomp\n");
	do_incompatiblemethod[key->unimeth[m]](key, m, res );
      }
      else {
	if (GENERAL_PRINTLEVEL>=7) PRINTF("extra %d\n",key->compatible);
	do_incompatiblemethod[key->unimeth[m]](key, m, part_result
					       );
	for (i=0; i<key->totalpoints; i++) {
	  res[i] += part_result[i];
	}
      }
    }
  } // for m
  
  if (key->mean!=0.0) {
    for(i=0; i<key->totalpoints; i++) res[i] += key->mean;
  }  
  
  if (part_result!=NULL) free(part_result);
  if (! (key->active = GENERAL_STORING && key->storing)) DeleteKey(keyNr);
  PutRNGstate();
  return;
  
 ErrorHandling: 
  PutRNGstate();
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,*error);
  if (part_result!=NULL) free(part_result);
  key->active = false;
  DeleteKey(keyNr);
  return;
}

void SimulateRF(Real *x, Real *T, int *dim,int *lx, int *grid, 
		int *Time, 
		int *covnr, Real *ParamList, int *nParam,
		Real *mean,
		int *ncov, int *anisotropy, int *op,
		int *method, int *distr,
		int *keyNr,
		Real *res, int *error)
{
  // grid   : boolean; if true then x,y,z are supposed to have three values: 
  //                      start, end, step
  //                   if false then (x,y,z) is a list of points where the 
  //                      value of the random field should be simulated
  // x,y,z  : see grid; if z==NULL it is assumed that the field is at most 
  //                    2 dimensional
  //                    if y==z==NULL, the random field is 1 dimensional
  // covnr  : nr of covariance function by getCovFct,
  // ParamList: vector with 5 to 8 elements
  //          see RFsimu.h, the definition of MEAN..KAPPAII, for the meaning 
  //          of the components
  // method : if <0 then programme chooses automatically for the method
  //          if ==-1 it starts with the preferred ("first") method for the 
  //                  specific covariance fct
  //            else starts search with previous successful method (if any)
  //          otherwise see list `SimulationType' in RFsimu.h
  // keyNR  : label where intermediate results should be stored
  // res : array/grid of simulation results

  InitSimulateRF(x, T, dim, lx, grid, Time, covnr, ParamList, nParam,
                 mean, ncov, anisotropy, op, method, distr, keyNr, error);
  if (*error==NOERROR) {
    DoSimulateRF(keyNr,res,error);
  } else {    
    if (GENERAL_PRINTLEVEL>0) PRINTF(" ** All methods have failed. **\n");
  } 
}
 
