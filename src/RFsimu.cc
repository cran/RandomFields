/*
 Authors
 Yindeng, Jiang, jiangyindeng@gmail.com
 Martin Schlather, schlath@hsu-hh.de

 library for unconditional simulation of stationary and isotropic random fields

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
//#include <sys/timeb.h>
#include <assert.h>
#include <string.h>
#include "RFsimu.h"
#include "RFCovFcts.h"
#include <unistd.h>


int FirstCheck_Cov(key_type *key, SimulationType Method, param_type param, 
		   bool No_Multiply, int *covnr, int *multiply, 
		   unsigned short int *actcov) {
  int v;
  *actcov=0;
  for (v=0; v<key->ncov; v++) {
    // printf("%d %d %d %d\n", v, key->method[v], Method, (int) key->left[v]);
    if ((key->method[v]==Method) && (key->left[v])) {
      /* Variance==0 is not eliminated anymore! -- maybe this could be improved
	 && (key->param[v][VARIANCE]>0)) { 
	 do not remove the parenths around Method! */
      key->left[v] = false;
      assert((key->covnr[v] >= 0) && (key->covnr[v] < currentNrCov));
      assert(key->param[v][VARIANCE] >= 0.0);
      covnr[*actcov] = key->covnr[v];
      if (CovList[covnr[*actcov]].implemented[Method] <= NOT_IMPLEMENTED) {
	return ERRORNOTDEFINED;
      }
      memcpy(param[*actcov], key->param[v], sizeof(double) * key->totalparam);
      if (*actcov>0) { /* *actcov>0 not v>0 since next line *actcov-1 --
			 check not for all v*/
	if (No_Multiply) {
	  if ((multiply[*actcov-1] = key->op[v-1])) { 
	    if (key->method[v-1] != Method) {
	      PRINTF("severe error - contact author. %d %d %d %d (%s) %d (%s)n",
		     v, key->op[v-1], key->ncov, key->method[v-1],
		     METHODNAMES[key->method[v-1]],
		     Method, METHODNAMES[Method]);
	      assert(false);
	    }
	    return ERRORNOMULTIPLICATION;
	  }
	} else { 
	  if ((multiply[*actcov-1] = key->op[v-1]) && key->method[v-1]!=Method){
	    if (GENERAL_PRINTLEVEL>0) 
	      PRINTF("severe error - contact author. %d %d %d %d (%s) %d (%s)n",
		     v, key->op[v-1], key->ncov, key->method[v-1],
		     METHODNAMES[key->method[v-1]],
		   Method, METHODNAMES[Method]);
	    return ERRORMETHODMIX; 
	  }
	} // not (No_Multiply)
      } // *actcov > 0
      if (!key->anisotropy)
	assert(fabs(key->param[v][SCALE] * key->param[v][INVSCALE]-1.0)<EPSILON);
      (*actcov)++;
    }
  }
  return (*actcov==0) ? /* no covariance for the considered method found ? */
    NOERROR_ENDOFLIST : NOERROR;
}

void DeleteKeyNotTrend(key_type *key)
{
  int d, m;
  if (key->x[0]!=NULL) free(key->x[0]);
  for (d=0; d<MAXDIM; d++) {key->x[d]=NULL;}
  for (m=0; m<MAXCOV; m++) {
    if (key->destruct[m]!=NULL) {
      key->destruct[m](&(key->S[m]));
      key->destruct[m]=NULL;
    }
  }
  if (key->destructExtremes!=NULL) {
    key->destructExtremes(&(key->SExtremes));
    key->destructExtremes=NULL;
  }
}

void DeleteKeyTrend(key_type *key)
{ 
  key->TrendModus = -1;
  if (key->TrendFunction!=NULL) free(key->TrendFunction); 
  key->TrendFunction=NULL;
  if (key->LinearTrend!=NULL) free(key->LinearTrend); 
  key->LinearTrend=NULL;
  key->lTrendFct = 0;
  key->lLinTrend = 0;
}

void DeleteKey(int *keyNr)
{  
   if (GENERAL_PRINTLEVEL>=4) { 
    PRINTF("deleting stored parameters of register %d ...\n", *keyNr);
  }
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {return;}
  key_type *key = &(KEY[*keyNr]);
  DeleteKeyNotTrend(key);
  DeleteKeyTrend(key);
  key->active=false;
}


void DeleteAllKeys() {
  int i;
  for (i=0;i<MAXKEYS;i++) {DeleteKey(&i);}
}


void printkey(key_type *key) {
  int i,j;  

  PRINTF("\ngrid=%d, active=%d, anisotropy=%d, tbm_method=%d lx=%d sp_dim=%d,\ntotalpoints=%d, distr=%s Time=%d [%f %f %f]\ntimespacedim=%d compatible=%d ncov=%d\n",	
	 key->grid, key->active, key->anisotropy,
	 key->tbm_method, key->lx, key->spatialdim, key->totalpoints,
	 DISTRNAMES[key->distribution], 
	 key->Time, key->T[0], key->T[1], key->T[2], 
	 key->timespacedim, key->compatible,
	 key->ncov);
  PRINTF("\n");
  for (i=0; i<key->ncov; i++) {
    PRINTF("%d: meth=%s, cov=%d (%s) ^S=%d ^destruct=%d\n left=%d unimeth=%d (%s) [%d], param=", i,
	   METHODNAMES[key->method[i]],
	   key->covnr[i], CovList[key->covnr[i]].name,
	   key->S[i],key->destruct[i], 
	   key->left[i], 
	   key->unimeth[i], METHODNAMES[key->unimeth[i]],
	   key->unimeth_alreadyInUse[i]
	   );
    for (j=0; j<key->totalparam; j++) PRINTF("%f,", key->param[i][j]);
    if (i < key->ncov - 1) PRINTF("\n%s\n",OP_SIGN[key->op[i]]); 
    else PRINTF("\n\n");
  }
  //
  for (i=0; i<MAXDIM; i++) 
  PRINTF("%d : len=%d \n", i, key->length[i]);

  PRINTF("^x=%d ^y=%d ^z=%d ^SExtremes=%d ^destructExtremes=%d \n", 
	 key->x[0], key->x[1], key->x[2], key->SExtremes, key->destructExtremes
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

int Transform2NoGrid(key_type *key, double* param, int TrueDim, 
		     int *start_param, double **xx)
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
  double t, *x;
  dimM1 = key->timespacedim - 1;
  assert(xx!=NULL);
  if (key->anisotropy) {
    // used: TrueDim, start_param
    total = TrueDim * key->totalpoints;
    // total number of coordinates to store

    if ((x = *xx = (double*) malloc(sizeof(double) * total))==NULL)
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
	double y[MAXDIM]; /* current point within grid, but without
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
      if ((x = *xx = (double*) malloc(sizeof(double) * TrueDim * key->totalpoints)
	   )==NULL)
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


int InternalGetGridSize(double *x[MAXDIM], int *dim, int *lx)
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
void GetGridSize(double *x,double *y,double *z,int *dim,int *lx,int *ly,int *lz)
{ 
  // should now be consistent with `seq' of R
  // if not, please report!!!!!!!!!!!! 
  double *xx[MAXDIM];
  int lxx[MAXDIM];
  xx[0]=x; xx[1]=y; xx[2]=z;
  if (InternalGetGridSize(xx,dim,lxx)) {
    *lx = *ly = *lz = 0;
  } else {
    *lx = lxx[0];  *ly = lxx[1];  *lz=lxx[2];
  }
}

#define COV_VARIO(FORMULA,FCT)\
  int v, aniso, k, j, endfor, ncovM1;\
  double  var,  zz, zw, result, d;\
  cov_fct *cov=NULL; /* no meaning -- just avoids warning from -Wall */ \
  zw = result = 0.0;\
  ncovM1 = ncov - 1;\
  if (anisotropy) {\
    double z[MAXDIM];\
    int cov_isotropy;\
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
    if (dim==1) d=fabs(x[0]);\
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
	zw *= cov->FCT(&zz,param[v], 1);\
	if ((v<ncovM1) && (op[v]==0)) break;\
      }\
      result += FORMULA;\
    }\
  }\
 return result;


double CovFct(double *x, int dim, int *covnr, int *op,
  param_type param, int ncov, bool anisotropy) {
  COV_VARIO(var * zw, cov);
}




double Vario(double *x, int dim, int *covnr, int *op,
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
// p and param must be of the same kind (double)!!! 
int CheckAndRescale( int covnr, int naturalscaling, double *p,
		    int timespacedim, int anisotropy, double* param) 
{ 
  // IMPORTANT!: p is the parameter vector of the R interface !!
  // output param : equals p, except naturalscaling
  int actparam, error, endfor=0 /*no meaning, avoids -Wall warning*/, i, j;
  double newscale;

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
    int maxdim, dummy;
    CovList[covnr].info(p, &maxdim, &dummy);
    if (maxdim < TrueDim) return ERRORCOVNOTALLOWED;
  }

  // specific parameter check for the covariance function
  if ((CovList[covnr].check!=NULL) && 
      ((error=(CovList[covnr].check(p, timespacedim, Nothing))))!=NOERROR){
    return error;
  }
 
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
void CovarianceNatSc(double *x, int *lx, int *covnr, double *p, int *np, 
		     int *logicaldim, /* timespacedim ! */
		     int *xdim,
		     int *ncov, int *anisotropy, int *op,
		     double *result, int *naturalscaling)
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

void Covariance(double *x,int *lx, int *covnr, double *p, int *np, 
		int *logicaldim, /* timespacedim ! */
		int *xdim,
		int *ncov, int *anisotropy, int *op,
		double *result)
{
  CovarianceNatSc(x,lx,covnr,p,np,logicaldim,xdim,
		  ncov, anisotropy, op,result, &GENERAL_NATURALSCALING);
}

static param_type Unchecked_param;
static int Unchecked_covnr[MAXCOV], Unchecked_DIM, Unchecked_op[MAXDIM],
  Unchecked_ncov;
static bool Unchecked_anisotropy;

void InitUncheckedCovFct(int *covnr, double *p,int *np, 
			 int *logicaldim, /* timespacedim ! */
			 int *xdim, int *ncov, int *anisotropy, int *op,
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

void UncheckedCovFct(double *x, int *lx, double *result)
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

void UncheckedCovMatrix(double *dist, int *lx, double *result)
// corresponds to CovarianceMatrixNatSc
{
  // note: column of x is point !!
  int i, ii, endfor, lxP1;
  long ve, ho;
  double var;
  
  var = CovFct(ZERO, Unchecked_DIM, Unchecked_covnr, 
	       Unchecked_op, Unchecked_param, 
	       Unchecked_ncov, Unchecked_anisotropy);
  lxP1 = *lx + 1;
  for (ii=*lx, i=0; ii>0; i+=lxP1, ii--) {
    result[i] = var;
    endfor = i + ii;
    for (ve = i + 1, ho = i + *lx; ve < endfor; ve++, ho += *lx, 
	   dist += Unchecked_DIM){
      result[ve] = result[ho] = 
	CovFct(dist, Unchecked_DIM, Unchecked_covnr, 
	       Unchecked_op, Unchecked_param, 
	       Unchecked_ncov, Unchecked_anisotropy);
    }
  }
  return;
}

void VariogramNatSc(double *x,int *lx,int *covnr,double *p,int *np,
		    int *logicaldim, /* timespacedim ! */
		    int *xdim,
		    int *ncov, int *anisotropy, int *op,
		    double *result, int *naturalscaling)
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

void Variogram(double *x,int *lx,int *covnr,double *p,int *np, 
	       int *logicaldim, /* timespacedim ! */
	       int *xdim,
	       int *ncov, int *anisotropy, int *op,
	       double *result)
{
  VariogramNatSc(x,lx,covnr,p,np,logicaldim,xdim,ncov,anisotropy,op,
		 result,&GENERAL_NATURALSCALING);
}
  
void CovarianceMatrixNatSc(double *dist,int *lx, int *covnr,double *p, int *np, 
			   int *logicaldim, /* timespacedim ! */
			   int *xdim,
			   int *ncov, int *anisotropy, int *op,
			   double *result, int *naturalscaling)
{
  // note: column of x is point !!
  int error=NOERROR, i, ii, pAdd, v, endfor, lxP1;
  long lxq, ve, ho;
  param_type param;
  double var;
 
 
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

void CovarianceMatrix(double *dist, int *lx,int *covnr, double *p, int *np,
		      int *logicaldim, /* timespacedim ! */
		      int *xdim,
		      int *ncov, int *anisotropy, int *op,
		      double *result)
{
  CovarianceMatrixNatSc(dist, lx, covnr, p, np, logicaldim, xdim, ncov, 
			anisotropy, op,	result, &GENERAL_NATURALSCALING);
}


void VariogramMatrix(double *dist,int *lx,int *covnr,double *p,int *np, double *result)
{ 
  // not programmed yet
  assert(false);
}

int insert_method(int* last_incompatible, key_type *key,
			 SimulationType m) {
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
    //
    // da waehrend der for-schleife unten eingefuegt wird:
    // immer moeglichst weit hinten einfuegen!
    (*last_incompatible)++;

    key->unimeth_alreadyInUse[key->n_unimeth]=
      key->unimeth_alreadyInUse[*last_incompatible];
    key->unimeth[key->n_unimeth]=key->unimeth[*last_incompatible];
    key->S[key->n_unimeth]=key->S[*last_incompatible]; // necessary
    // only if further methods are tried due to partial failure
    key->destruct[key->n_unimeth]=key->destruct[*last_incompatible];

    key->S[*last_incompatible] = NULL;
    key->destruct[*last_incompatible] = NULL;
    idx = *last_incompatible;
  } else idx=key->n_unimeth;
  key->unimeth[idx] = m;
  key->unimeth_alreadyInUse[idx] = false;
  key->n_unimeth++;
  return idx;
}

void delete_method(int *M, int *last_incompatible, key_type *key) {
  // note that key->n_unimeth had been decreased before this function
  // is called
  int idx;
  if (key->destruct[*M]!=NULL) key->destruct[*M](&(key->S[*M]));
  assert(key->S[*M]==NULL);
  key->destruct[*M]=NULL;
  if (do_incompatiblemethod[key->unimeth[*M]]!=NULL) {
    assert(*last_incompatible>=0);
    assert(*M<=*last_incompatible);
    key->unimeth[*M]=key->unimeth[*last_incompatible];
    key->unimeth_alreadyInUse[*M]= key->unimeth_alreadyInUse[*last_incompatible];
    key->S[*M] = key->S[*last_incompatible];
    key->S[*last_incompatible]=NULL;
    key->destruct[*M]=key->destruct[*last_incompatible];
    key->destruct[*last_incompatible]=NULL;
    idx = *last_incompatible;
    (*last_incompatible)--;
  } else idx = *M;
  key->n_unimeth--;
  if (idx==key->n_unimeth){ // deleting of the very last entry
    assert(*M==key->n_unimeth || *last_incompatible + 1 == key->n_unimeth); 
  } else {
    assert(idx>=0 && idx<key->n_unimeth);
    key->unimeth[idx] = key->unimeth[key->n_unimeth];
    key->unimeth_alreadyInUse[idx] = key->unimeth_alreadyInUse[key->n_unimeth];
    key->S[idx] = key->S[key->n_unimeth];
    key->S[key->n_unimeth]=NULL;
    key->destruct[idx] = key->destruct[key->n_unimeth];
    key->destruct[key->n_unimeth]=NULL;
  }
  key->unimeth_alreadyInUse[key->n_unimeth] = false;
  (*M)--;
}


/*

speed, exactness, stationarity

direct (lower -- always, upper -- never)


grid: circembed
nongrid: spectralTBM, TBM2, TBM3,


Ausnahmen: dimension beschraenkt

*/

SimulationType ML[MAXCOV][SimulationTypeLength]; // Method List
int NML[MAXCOV], IML[MAXCOV];
void init_method_list(key_type *key){
  int best_dirct=DIRECTGAUSS_BESTVARIABLES, 
    max_dirct=DIRECTGAUSS_MAXVARIABLES, 
    maxmem=500000000, 
    i, v, nr, nml;
  bool anyFirstDirect = false, anyFirstCE = false, CEbadlybehaved;

#define nLastDefault 2
  SimulationType LastDefault[nLastDefault] = {AdditiveMpp, Hyperplane}; 
  SimulationType DirectFst = key->totalpoints <= best_dirct ? Direct : Forbidden;
  SimulationType DirectLst = key->totalpoints <= max_dirct ? Direct : Forbidden;
#define nStandard 6
  SimulationType Standard[nStandard] = 
    {CircEmbed, CircEmbedIntrinsic, CircEmbedCutoff, SpectralTBM, TBM2, TBM3};

  bool Allowed[SimulationTypeLength], allowed[SimulationTypeLength];
  for (i=0; i<Nothing; i++) Allowed[i] = true; 
  if (DECISION_PARAM.stationary_only==DECISION_TRUE) 
    Allowed[CircEmbedIntrinsic]=false;
  if (!key->grid) Allowed[CircEmbed] = Allowed[CircEmbedCutoff] = 
		    Allowed[CircEmbedIntrinsic] = false;
  if (key->timespacedim > 2) 
    Allowed[SpectralTBM] = Allowed[TBM2] =  Allowed[Hyperplane] = false;
  if (DECISION_PARAM.exactness==DECISION_TRUE)
    Allowed[TBM2] = Allowed[TBM3] = Allowed[SpectralTBM] = 
      Allowed[AdditiveMpp] = Allowed[Hyperplane] = false;

#define nraw 2 * SimulationTypeLength
  SimulationType raw[nraw];
  for (v=0; v<key->ncov; v++) {
    IML[v] = -1;
    // special case: Nugget (part 1)
    if (CovList[key->covnr[v]].cov==nugget) continue;
    
    for (i=0; i<nraw; i++) raw[i] = Forbidden;
    int *implemented=CovList[key->covnr[v]].implemented, maxdim, CEbadlybehaved;
    for (i=0; i<Nothing; i++) 
      allowed[i] = Allowed[i] && implemented[i] >= IMPLEMENTED;
    if (DECISION_PARAM.stationary_only==DECISION_CASESPEC && 
	!CovList[key->covnr[v]].variogram) allowed[CircEmbedIntrinsic] = false;
    if (DECISION_PARAM.exactness != DECISION_FALSE && implemented[i]==NUM_APPROX)
      allowed[TBM2] = false;

    // multiplicative models
    if (v<key->ncov-1 && key->op[v]) {
      // nur die Reihenfolge fuer den ersten Faktor wird verwendet -- 
      // der letzte Faktor wird gesetzt wie wenn nicht Faktor (spaeter ignoriert)
      // alle anderen Faktoren werden als Faktoren gesetzt aber Methodenfolge
      //    wiederum ignoriert
      // next_element_of_preference_list achtet auf die anderen Faktoren
      for (i=0; i<Nothing; i++) 
	if (i!=TBM3 && i!=CircEmbed && i!=Direct) allowed[i] = false;
    }

    CovList[key->covnr[v]].info(key->param[v], &maxdim, &CEbadlybehaved);
    assert(key->timespacedim <= maxdim);
    if (key->grid && CEbadlybehaved && DECISION_PARAM.exactness==DECISION_TRUE) {
      // preference to Direct whenever possible?
      DirectFst = DirectLst;
    }
    nr = 0;
    raw[nr++] = DirectFst;
    if ((((key->ncov==1 || DECISION_PARAM.exactness==DECISION_FALSE) && 
	  CEbadlybehaved) || CEbadlybehaved > 1) && key->grid) {
      // key->ncov > 1 laesst hoffen, dass die andere Kovarianzfunktion
      // sich freundlicher verhaelt (z.B. wenn Nugget), so dass die Matrix
      // des Gesamtmodells sich gut bzgl. CE verhaelt.
      if (key->timespacedim==2) {
	raw[nr++] = SpectralTBM;
	if (DECISION_PARAM.exactness==DECISION_FALSE) raw[nr++] = TBM2;
      }
      if (key->timespacedim==3 && CEbadlybehaved>1) raw[nr++] = TBM3;
    }
    if (key->grid && DECISION_PARAM.exactness==DECISION_FALSE && 
	key->totalpoints * (1 << key->timespacedim) * 2*sizeof(double) > maxmem){
      raw[nr++] = SpectralTBM;
      raw[nr++] = TBM2;
      raw[nr++] = TBM3;
    }
    for (i=0; i<nStandard; i++) raw[nr++] = Standard[i];
    raw[nr++] = DirectLst;
    for (i=0; i<nLastDefault; i++) raw[nr++] = LastDefault[i];
    nml = 0;
    for (i=0; i<nr; i++) {
      if (raw[i]!=Forbidden && allowed[raw[i]]) {
	ML[v][nml++] = raw[i];
	allowed[raw[i]] = false;
      }
    }
    NML[v] = nml;
    switch (raw[0]) {
	case Direct : anyFirstDirect = true; break;
	case CircEmbed : anyFirstCE = true; break;
    }
  }

  // nugget, part 2
  for (v=0; v<key->ncov; v++) if (CovList[key->covnr[v]].cov==nugget){
    nml = 0;
    if (!key->anisotropy) {
      // CE first considered, since nugget may have a positive effect on
      // the CE simuation.
      if (anyFirstCE) ML[v][nml++] = CircEmbed;
      if (anyFirstDirect) ML[v][nml++] = Direct;
    }
    ML[v][nml++] = Nugget;
    NML[v] = nml;
  }
  
  if (GENERAL_PRINTLEVEL>5) {
    PRINTF("Lists of methods for the simulation\n===================================\n");
    for (v=0; v<key->ncov; v++) {
      PRINTF("%s: ", CovList[key->covnr[v]].name);
      for (i=0; i < NML[v]; i++) PRINTF("%s, ", METHODNAMES[ML[v][i]]);
      PRINTF("\n");
    }
  }
}
  

int next_element_of_method_list(key_type *key){
  int v;
  for (v=0; v<key->ncov; v++) {
    if (key->left[v]) {
      if (v>0 && key->op[v-1]) key->method[v] = key->method[v-1];
      else {
	(IML[v])++;
	if (IML[v] >= NML[v]) return v + 1;
	key->method[v] = ML[v][IML[v]];
      }
    }
  }
  return 0;
}


void InitSimulateRF(double *x, double *T, 
		    int *spatialdim, /* spacial dim only ! */
		    int *lx, int *grid, 
		    int *Time,
	            int *covnr, double *ParamList, int *nParam,
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

  bool user_defined, simple_definition;
  int i, last_incompatible, d, endfor, pAdd, v, m, M, err_occurred;
  unsigned long totalBytes;
  double p[TOTAL_PARAM], *PL;
  int method_used[(int) Nothing + 1];
  // preference lists, distinguished by grid==true/false and dimension
  // lists must end with Nothing!
  SimulationType Merr=Nothing;

  key_type *key;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORREGISTER; 
    goto ErrorHandling;
  }
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
 
  totalBytes =  sizeof(double) * *lx * *spatialdim;
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
      ( (*Time == key->Time) &&
	(!(*Time) || ! memcmp(key->T, T, sizeof(double) * 3))) &&
      (!memcmp(key->op, op, sizeof(int) * ((*ncov>1) ? *ncov-1 : 0))) &&
      (!memcmp(key->x[0], x, totalBytes)) 
      ;
    if (GENERAL_PRINTLEVEL>5 && !key->active) {
      PRINTF("failure with previous initialisation at lx=%d spatialdim=%d grid=%d distrib=%d ncov=%d covnr=%d   aniso=%d time=%d key->T=%d op=%d x=%d\n",
	     key->lx==*lx, key->spatialdim==*spatialdim, key->grid==(bool)*grid,
	     key->distribution==*distr, key->ncov==*ncov, 
	     !memcmp(key->covnr, covnr, sizeof(int) * *ncov), 
	     key->anisotropy==*anisotropy,
	     *Time == key->Time, 
	     !(*Time) || ! memcmp(key->T, T, sizeof(double) * 3),
	     !memcmp(key->op, op, sizeof(int) * ((*ncov>1) ? *ncov-1 : 0)),
	     !memcmp(key->x[0], x, totalBytes));
    }
  }

  if (!key->active) {
    DeleteKeyNotTrend(key);
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
    // operators +, *
    memcpy(key->op, op,  sizeof(int) * ((*ncov>1) ? *ncov-1 : 0));
    // coordinates

    if (key->x[0]!=NULL) free(key->x[0]);
    if ((key->x[0]=(double*) malloc(totalBytes))==NULL){
      *error=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    memcpy(key->x[0], x, totalBytes);

    for (d=1; d<*spatialdim; d++) key->x[d]= &(key->x[0][d * *lx]);
    for (; d<MAXDIM; d++)  key->x[d]=NULL;
    
    //Time
    if (key->Time = *Time) {
      if (!key->anisotropy) { *error = ERRORTIMENOTANISO; goto ErrorHandling; }
      memcpy(key->T, T, sizeof(double) * 3);
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
				 key->timespacedim, *anisotropy, p)))!= NOERROR){
       goto ErrorHandling; 
    }
    for (j=KAPPA+CovList[covnr[v]].kappas; j<=LASTKAPPA; j++) p[j]=0.0;
    if (!(key->active &=
	  (key->covnr[v]==covnr[v]) &&
	  (!memcmp(key->param[v], p, sizeof(double) * key->totalparam))
	  )) {
      key->covnr[v]=covnr[v];
      memcpy(key->param[v], p, sizeof(double) * key->totalparam);
    }
   PL+=pAdd+CovList[covnr[v]].kappas;
  }

  // method list must be treated separately since method is integer coded
  // when passed to InitSimulate !
  if (user_defined = (method[0]>=0)) {//either all or none shld be user defined
    for (v=0; v<key->ncov; v++) {
      if (method[v] >= (int) Nothing) {
	*error=ERRORNOTDEFINED; goto ErrorHandling;}
      if (v>0 && op[v-1] && method[v]!=method[v-1]) {
	*error=ERRORMETHODMIX; goto ErrorHandling;}
      if (DECISION_PARAM.stationary_only==DECISION_TRUE && 
	  key->method[v]==CircEmbedIntrinsic) {
	*error=ERRORMETHODEXCLUDED; goto ErrorHandling;}
      NML[v] = 1;
      IML[v] = -1;
      ML[v][0] = (SimulationType) method[v];
      if (key->method[v] != ML[v][0]) key->active = false; 
      // sonst: ML info wird nicht verwendet...
    }
  }

  if (key->active) { // detected that all relevant parameters agree, 
    //              no initialization is necessary
    if (GENERAL_PRINTLEVEL>=2) ErrorMessage(Nothing, USEOLDINIT);
    *error=NOERROR; 
   return; /* ***** END ***** */ 
  }

  // Delete stored, intermediate result if there is any
  // folgende Zeilen unschoen, da weiter vorne bereits DeleteKeyNotTrend
  // aufgerufen wurde
  for (v=0; v<MAXCOV; v++) {
    if (key->destruct[v]!=NULL) {
      key->destruct[v](&(key->S[v]));
      key->destruct[v] = NULL;
    }
  }
  if (key->destructExtremes!=NULL) {
    key->destructExtremes(&(key->SExtremes));
    key->destructExtremes = NULL;
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
    for (key->spatialtotalpoints=1, d=0; d<key->spatialdim; d++) {
      key->spatialtotalpoints *= key->length[d];
    }
    key->totalpoints = key->spatialtotalpoints;
    if (key->Time) key->totalpoints *= key->length[key->spatialdim];
  } else { // not grid
    key->totalpoints = key->spatialtotalpoints = key->length[0]= *lx;
    if (key->Time) {
      int Tdim; double *Tx[MAXDIM];
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
  }
  for (d=key->timespacedim; d<MAXDIM; d++) key->length[d] = (int) RF_NAN; // 1

  // simple definition of the covariance function?
  // (the only definition used up to version 1.0 of RandomFields) 
  simple_definition = key->ncov==1 ||
    (key->ncov==2 && !key->anisotropy && (key->method[0]==Nugget ||
					  key->method[1]==Nugget));

  err_occurred = NOERROR;
  *error = ERROROUTOFMETHODLIST;
  key->n_unimeth = 0;
  last_incompatible = -1;
  for (v=0; v<*ncov; v++) key->left[v]=true;
  //           indicates which covariance function are left to be initialized
  if (!user_defined) init_method_list(key); // nicht weiter vorne,
  // da es auf Daten von key zugreift, die weiter vorne noch nicht
  // vorhanden sind
  
  while (*error != NOERROR) {
    int error_covnr;
    if (*error != ERRORANYLEFT)
      error_covnr = next_element_of_method_list(key);
    if (error_covnr) {
      if (GENERAL_PRINTLEVEL>1 && !user_defined)
	PRINTF("covariance nr %d (%s) run out of method list.\n", 
	       error_covnr, CovList[key->covnr[error_covnr - 1]].name);
      *error = user_defined ? err_occurred : ERROROUTOFMETHODLIST;
      break;
    }
      
    *error = err_occurred = NOERROR;// necessary in case ncov=0, since then 
    //                          the following loop is not entered
    // create list of methods used
    for (i=0; i<=(int) Nothing; i++) method_used[i] = 0;
    for (v=0; v<key->ncov; v++) method_used[key->method[v]] += key->left[v]; 
    for (m=CircEmbed; m<=Nothing; m++) 
      if (method_used[m]) {
	insert_method(&last_incompatible, key, (SimulationType) m); 
	//                increments key->n_unimeth
	if (method_used[m]==1) {
	  v=0; 
	  while(key->method[v]!=m || !key->left[v]) v++;
	  ML[v][IML[v]] = Nothing; 
	  // es lohnt sich nicht, die Methoden durchzuleiern,
	  // die auf die Kovarianzfunktion bereits angewandt wurden und
	  // diese Kovarianzfunktion dabei die einzige war
	}
      }
    
    for (M=0; M<key->n_unimeth; M++) if (!key->unimeth_alreadyInUse[M]) {
      bool left[MAXCOV];
      key->unimeth_alreadyInUse[M] = true;
      for (v=0; v<key->ncov; v++) left[v] = key->left[v]; // sicherung,
      // da init_method ueberschreibt und *danach* entscheidet, ob
      // die methode funktioniert -- ueber key->method ist es zu schwierig
      // da eine methode mehrmals verwendet werden kann
      if (GENERAL_PRINTLEVEL>=4)
	PRINTF("\nInit:  method %d (%s); currently %d method(s) in use\n",
	       M, METHODNAMES[key->unimeth[M]], key->n_unimeth);
      assert(key->destruct[M]==NULL);
      assert(key->S[M]==NULL);
      if (key->unimeth[M]>Nothing) {
	*error=ERRORMETHODNOTALLOWED; 
	goto ErrorHandling;
      }
      err_occurred = init_method[key->unimeth[M]](key, M);
      if (GENERAL_PRINTLEVEL>=3) {
	PRINTF("return code %d :", err_occurred);
	ErrorMessage(Nothing, err_occurred);
      }
      switch (err_occurred) {
	  case NOERROR_REPEAT:
	    insert_method(&last_incompatible, key, key->unimeth[M]); // ??!!
	    err_occurred = NOERROR; 
	    break;
	  case NOERROR:
	    break;
	  case NOERROR_ENDOFLIST:
	    delete_method(&M, &last_incompatible, key);
	    err_occurred = NOERROR;
	    break;
	  default:
	    Merr = key->unimeth[M]; // do not shift after delete_method, since 
	    // the latter changes the value of M !
	    delete_method(&M, &last_incompatible, key);
	    *error = err_occurred; // error may occur anywhere
	    //                     and final err_occured might be 0
	    for (v=0; v<key->ncov; v++) key->left[v] = left[v];	 
	    if (simple_definition && GENERAL_PRINTLEVEL>2) 
	      ErrorMessage(Merr,*error);
      }
    } // for (M...)
    
    if (*error==NOERROR)
      for (v=0; v<key->ncov; v++) if (key->left[v]) {
	*error = ERRORANYLEFT;
	Merr = Nothing;
	break;
      }  // if (*error == ERRORANYLEFT) continue;

  } // while *error!=NOERROR
   
  if (*error!=NOERROR && !user_defined) {
    int zaehler;
    if (GENERAL_PRINTLEVEL>1) {
      PRINTF("algorithm partially failed. Retrying. Last error has been:\n");
      ErrorMessage(Merr, err_occurred);   
    }
    zaehler = 0;
    for (v=0; v<key->ncov; v++) if (key->left[v]) {
      key->method[v] = Nothing;
      zaehler++;
    }
    assert(zaehler > 0);
    for (v=0; v<key->ncov; v++) if (key->left[v]) {
      if (NML[v]==0) {
	*error = ERROROUTOFMETHODLIST;
	goto ErrorHandling;
      }
      if (GENERAL_PRINTLEVEL>=5) 
	PRINTF("'%s' not initiated yet\n", CovList[key->covnr[v]].name);
      for (i=0; i<NML[v]; i++) {
	key->method[v] = ML[v][i];
	if (key->method[v] == Nothing) {
	  if (GENERAL_PRINTLEVEL>=6) PRINTF("method #%d jumped\n", i);
	  continue;
	}
	M = insert_method(&last_incompatible, key, ML[v][i]); 
	*error = init_method[key->unimeth[M]](key, M);
	if (GENERAL_PRINTLEVEL>=3) {
	  PRINTF("%s has return code %d;", METHODNAMES[ML[v][i]], *error);
	  ErrorMessage(Nothing, *error);
	}
	if (*error != NOERROR && *error != NOERROR_REPEAT) {
	  Merr = key->method[v];
	  delete_method(&M, &last_incompatible, key);
	  key->left[v] = true;
	} else {
	  *error = NOERROR;
	  break;
	}
      }
      if (*error != NOERROR) {
	if (GENERAL_PRINTLEVEL>1) {
	  PRINTF("algorithm failed. Last error has been:\n");
	  ErrorMessage(Merr, *error);   
	}
	break;
      }
    }
  }
   
  if (!(key->active = *error==NOERROR)) goto ErrorHandling;
  if (GENERAL_PRINTLEVEL>=4) PRINTF("Successfully initialised.\n");
  key->compatible = last_incompatible <= 0; // if at most one incompatible
  // method, then no intermediate results must be stored
  return;
  
  ErrorHandling:
  key->active = false;
  if (GENERAL_PRINTLEVEL>0) {
    ErrorMessage(Merr, *error);
    PRINTF("\n");
  }
}

void AddTrend(int *keyNr, int *n, double *res, int *error)
   {
  int i;
  key_type *key;

  *error = NOERROR;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  key = &(KEY[*keyNr]);
  if (!key->active) {*error=ERRORNOTINITIALIZED; goto ErrorHandling;}

  switch (key->TrendModus) {
      case TREND_MEAN :
	if (key->mean!=0.0) {
	  long int endfor= *n * key->totalpoints;
	  for(i=0; i<endfor; i++) res[i] += key->mean;
	}
	break;
      case TREND_LINEAR :
      case TREND_FCT :
      case TREND_PARAM_FCT :
      default: assert(false); // not programmed yet
  }
  return;

 ErrorHandling: 
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,*error);
 return;
}
 
void DoPoissonRF(int *keyNr, int *pairs, int *n, double *res, int *error)
{
  // not programmed yet
  assert(false);
}

void DoSimulateRF(int *keyNr, int *n, int *pairs, double *res, int *error)
{
  // uses that res[...] = 0.0
  double  *part_result, *orig_res=res;
  long i, m;
  int ni, nn, digits, each=0;
  char back[]="\b\b\b\b\b\b\b\b\b\b\b", format[20], prozent[]="%";
  key_type *key=NULL;

  part_result = NULL;
  nn = *n;
  *error=NOERROR; 
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  key = &(KEY[*keyNr]);
  if (!key->active) {*error=ERRORNOTINITIALIZED; goto ErrorHandling;}
  if (key->n_unimeth==0) return; // using that res[...] = 0.0 by initialisation
  
  if (nn>1 && GENERAL_PCH[0] == '!') {
    digits = (nn<900000000) ? 1 + (int) trunc(log(nn) / log(10)) : 9;
    back[digits] = '\0';
    each = (nn<100) ? 1 : nn / 100;
    sprintf(format, "%ss%s%dd", prozent, prozent, digits);
  }

  GetRNGstate();
  if (!key->compatible) {
    if ((part_result=(double *) malloc(sizeof(double) * key->totalpoints))
	==NULL) {*error=ERRORMEMORYALLOCATION; goto ErrorHandling;}
  }

  for (ni=0; ni<nn; ni++, res += key->totalpoints) {
    if (nn>1)
      if (GENERAL_PCH[0] != '!') PRINTF("%s", GENERAL_PCH);
      else if (ni % each ==0) PRINTF(format, back, ni);

    for(m=0; m<key->n_unimeth; m++) {
      if (GENERAL_PRINTLEVEL>=7)
	PRINTF("DoSimu %d %s\n",m, METHODNAMES[key->unimeth[m]]);
      if (do_compatiblemethod[key->unimeth[m]]!=NULL) {
	if (GENERAL_PRINTLEVEL>=7) PRINTF("compatible: add=%d\n",m!=0);
	do_compatiblemethod[key->unimeth[m]] (key, m!=0, m, res); 
      } else {
	if (m==0) {
	  if (GENERAL_PRINTLEVEL>=7) PRINTF("incomp\n");
	  do_incompatiblemethod[key->unimeth[m]](key, m, res);
	}
	else {
	  if (GENERAL_PRINTLEVEL>=7) PRINTF("extra %d\n",key->compatible);
	  do_incompatiblemethod[key->unimeth[m]](key, m, part_result);
	  for (i=0; i<key->totalpoints; i++) res[i] += part_result[i];
	}
      }
    } // for m

    if (*pairs) {
      double* res_pair;
      res_pair = res + key->totalpoints;
      for (i=0; i<key->totalpoints; i++) res_pair[i] = -res[i];
      res = res_pair;
    }
  } // for n

  if (nn>1 && GENERAL_PCH[0] == '!') PRINTF("%s", back);

  if (*pairs) nn = 2 * nn;
  AddTrend(keyNr, &nn, orig_res, error);
  if (*error) goto ErrorHandling;
 
  if (part_result!=NULL) free(part_result);
  if (!(key->active = GENERAL_STORING && key->storing)) DeleteKey(keyNr);
  PutRNGstate();
  return; 
  
 ErrorHandling: 
  PutRNGstate();
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,*error);
  if (part_result!=NULL) free(part_result);
  if (key!=NULL) {
    key->active = false;
    DeleteKey(keyNr);
  }
  return;
}

