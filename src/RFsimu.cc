/*
Authors
 Yindeng, Jiang, jiangyindeng@gmail.com
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 library for unconditional simulation of stationary and isotropic random fields

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2006 Martin Schlather

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


int FirstCheck_Cov(key_type *key, int m, bool MultiplyAndHyper)
{
  int v, error;
  SimulationType Method;
  covinfo_type *kc;
  methodvalue_type *meth;

  meth = &(key->meth[m]);
  Method = meth->unimeth;
  meth->actcov=0;
  error = NOERROR;
  for (v=0; v<key->ncov; v++) {
    ERRORMODELNUMBER = v;	
    cov_fct *cov;
    kc = &(key->cov[v]);
    cov = &(CovList[kc->nr]); 
    meth->covlist[meth->actcov] = v;
    // printf("%d %d %d %d\n", v, kc->method, Method, (int) kc->left);
    if (kc->method==Method && kc->left) {
      assert((kc->nr >= 0) && (kc->nr < currentNrCov));
      assert(kc->param[VARIANCE] >= 0.0);
      if (cov->implemented[Method] != IMPLEMENTED) return ERRORNOTDEFINED; 
      if (cov->type==ISOHYPERMODEL || cov->type==ANISOHYPERMODEL) {
	if (!MultiplyAndHyper) error = ERRORHYPERNOTALLOWED;
	else {
	  error = 
	      cov->checkNinit(kc, &(COVLISTALL[v]), key->ncov - v - 1, Method, 
			      key->anisotropy);
	}
      } else {
	error = cov->check(kc->param, kc->reduceddim, kc->dim, Method);
      }
      if (error != NOERROR) return error;
//      printf("v=%d %d %d \n", v,  key->ncov, key->cov[v].op);
      if (v < key->ncov -1 && key->cov[v].op) { 
	if (key->cov[v+1].method != Method) {
	  if (MultiplyAndHyper) {
	    if (GENERAL_PRINTLEVEL>0) 
	      PRINTF("severe error -contact author. %d %d %d %d (%s) %d (%s)n",
		     v, key->cov[v-1].op, key->ncov, key->cov[v-1].method,
		     METHODNAMES[key->cov[v-1].method],
		     Method, METHODNAMES[Method]);
	    return ERRORMETHODMIX; 
	  } else {
	    PRINTF("severe error - contact author. %d %d %d %d (%s) %d (%s)n",
		   v, key->cov[v-1].op, key->ncov, key->cov[v-1],
		   METHODNAMES[key->cov[v-1].method],
		   Method, METHODNAMES[Method]);
	    assert(false);
	  }
	}
	if (!MultiplyAndHyper) return ERRORNOMULTIPLICATION;
      } // meth->actcov > 0
      int w, nsub;
      nsub = (cov->type==ISOHYPERMODEL || cov->type==ANISOHYPERMODEL) ? 
	  (int) kc->param[HYPERNR] : 0;
      for (w=0; w<nsub; w++) {
	kc->left = false;
	kc++;
	(meth->actcov)++; 
	v++;
	meth->covlist[meth->actcov] = v;
      }
      kc->left = false;
      (meth->actcov)++; 
    } // left
  }
  ERRORMODELNUMBER = -1;	

  return (meth->actcov==0 
      ? (error==NOERROR ? NOERROR_ENDOFLIST : error)
      : (error==NOERROR ? NOERROR : NOERROR_REPEAT));
}

void COVINFO_NULL(covinfo_type *cov){
  int m, d;
  double *param;
  covinfo_type *kc;
  for (m=0; m<MAXCOV; m++) {
    kc = &(cov[m]);
    kc->x = NULL;
    kc->reduceddim = kc->op = 0;
    kc->nr = currentNrCov; // undefined model
    kc->method = Forbidden;
    kc->genuine_last_dimension = kc->simugrid = false;
    kc->left = true;
    param = kc->param;
    for (d=0; d<MAXDIM; d++) {
      kc->genuine_dim[d] = false;
      kc->length[d] = kc->idx[d] = -1;
    }
    for (d=MAXDIM * 3 - 1; d>=0; d--) kc->xsimugr[d] = RF_NAN;
    for (d=0; d<MAXDIMSQ; d++) param[d] = kc->aniso[d] = RF_NAN;
    for (; d<TOTAL_PARAM; param[d++]=RF_NAN);
  }
}

void KEY_NULL(key_type *key)
{
  int m, d;
  key->stop = key->active = false;
  key->ncov = key->n_unimeth = 0;
  key->spatialdim = key->timespacedim = key->totalparam = -1;
  key->totalpoints = key->spatialtotalpoints = 0;
  for (m=0; m<MAXCOV; m++) {
    key->meth[m].S = NULL;
    key->meth[m].destruct = NULL;
  }
  COVINFO_NULL(key->cov);
  for (d=0; d<MAXDIM; d++) {
      key->x[d]=NULL;
      key->length[d] = -1;
  }
  key->T[0] = key->T[1] = key->T[2] = 0.0; 
  //    key->naturalscaling =-1;
  key->TrendModus = -1;
  key->TrendFunction = NULL;
  key->LinearTrend = NULL;
  key->mean = 0.0;
}

void DeleteKeyNotTrend(key_type *key)
{
  int d, m;
  if (key->x[0]!=NULL) free(key->x[0]);
  for (d=0; d<MAXDIM; d++) {key->x[d]=NULL;}
  for (m=0; m<MAXCOV; m++) {
    if (key->cov[m].x != NULL) {
      free(key->cov[m].x);
      key->cov[m].x = NULL;
    }
    if (key->meth[m].destruct!=NULL) {
      key->meth[m].destruct(&(key->meth[m].S));
      key->meth[m].destruct=NULL;
    }
  }
  key->active = false;
  key->ncov = key->n_unimeth = 0;
  key->spatialdim = key->timespacedim = key->totalparam = -1;
  key->totalpoints = key->spatialtotalpoints = 0;
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
  covinfo_type *keycov;
  methodvalue_type *meth;

  PRINTF("\ngrid=%d, active=%d, anisotropy=%d, sp_dim=%d,\ntotalpoints=%d, distr=%s Time=%d [%f %f %f]\ntimespacedim=%d compatible=%d ncov=%d\n",	
	 key->grid, key->active, key->anisotropy,
	 key->spatialdim, key->totalpoints,
	 DISTRNAMES[key->distribution], 
	 key->Time, key->T[0], key->T[1], key->T[2], 
	 key->timespacedim, key->compatible,
	 key->ncov);
  PRINTF("\n");
  for (i=0; i<key->ncov; i++) {
    keycov = &(key->cov[i]);
    PRINTF("%d: meth=%s, cov=%d (%s) left=%d, param=", i,
	   METHODNAMES[keycov->method],
	   keycov->nr, CovList[keycov->nr].name,
	   keycov->left);
    for (j=0; j<key->totalparam; j++) PRINTF("%f,", keycov->param[j]);
    if (i < key->ncov - 1) PRINTF("\n%s\n",OP_SIGN[keycov->op]); 
    else PRINTF("\n\n");
  }
  for (i=0; i<key->ncov; i++) {
    meth = &(key->meth[i]);
    PRINTF("%d:  ^S=%d ^destruct=%d\n unimeth=%d (%s) [%d]\n", i,
	   meth->S, meth->destruct, 
	   meth->unimeth, METHODNAMES[meth->unimeth],
	   meth->unimeth_alreadyInUse
	   );
  }
  //
  for (i=0; i<MAXDIM; i++) 
  PRINTF("%d : len=%d \n", i, key->length[i]);

  PRINTF("^x=%d ^y=%d ^z=%d \n", key->x[0], key->x[1], key->x[2]);

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

int Transform2NoGrid(key_type *key, int v) {
  int d;
  covinfo_type *kc;
  kc = &(key->cov[v]);
  if (kc->simugrid) {
    for (d=0; d<key->timespacedim; d++) {
      kc->length[d] = kc->genuine_dim[d] ? key->length[d] : 1;
//      printf("%d %d %d %d \n", d, 
//	       kc->length[d], kc->genuine_dim[d], key->length[d]);
    }
    Getxsimugr(key->x, kc->param, key->timespacedim, key->anisotropy, 
	       kc->xsimugr);
  }
  return Transform2NoGrid(key, kc->aniso, kc->reduceddim,
			  kc->simugrid, &(kc->x));
}

int Transform2NoGrid(key_type *key, aniso_type aniso, int reduceddim,
		     bool simugrid, double **xx)
{
  /* 
     this function transforms the coordinates according to the anisotropy
     matrix given in param (usually param is the param the first anisotropy
     matrix of the covariance product -- all the other anisotropy matrices
     have to multiple of the first, and this parameter is put into the
     covariance function as scale parameter 

     developtime only considered if !key->simugrid && key->Time
     
     output : keycov->xx
           isotropic   : grid : triple devided by scale ==multip. with aniso[0]!!
                         arbitrary : x-coordinates devided by scale
           anisotropic : if grid and aniso=diag then
                            triples are multiplied with orig param[ANISO...]
                         else 
                            x-coordinates multiplied explicitely with aniso: x*a;
                            dim-reduced
  */
  int i,j,k,n,d,w,dimM1;
  long endtime, total;
  double t, *x;

//  printf("%d %d %d %d %f %f %f %d\n",
//	 key->timespacedim,key->anisotropy,key->Time,key->grid,
//	 key->T[0],key->T[1],key->T[2],key->length[0]
//    );

  dimM1 = key->timespacedim - 1;
//  assert(*xx == NULL);
  if (*xx != NULL) {
    if (GENERAL_PRINTLEVEL>-4) 
      PRINTF("Note: x values had already been transformed.");
    return NOERROR;
  }
  // total number of coordinates to store
  total = simugrid 
    ? (key->timespacedim * 3)
    : (reduceddim * key->totalpoints);
//  printf("2nogr %d %d %d %d \n", total, reduceddim, simugrid, key->totalpoints);
//  printf("total %d %d %d %d\n", total, key->totalpoints, reduceddim, key->timespacedim);

  if ((x = *xx = (double*) malloc(sizeof(double) * total)) == NULL)
    return ERRORMEMORYALLOCATION;

  if (key->anisotropy) {
    /* determine points explicitly */
    if (key->Time && !key->grid) { // time component, no grid
      endtime =  key->length[key->spatialdim]; 
      for (k=j=0, t=key->T[XSTART]; j<endtime; j++, t += key->T[XSTEP]) {
	for (i=0; i<key->length[0]; i++) {
	  for (n=d=0; d<reduceddim; d++, k++) {
 	    x[k] = 0.0;
	    for(w=0; w<key->spatialdim; w++) {
//	      printf("k=%d n=%d, w=%d, i=%d\n", k, n, w, i);
//	      printf("xk=%f aniso=%f %f\n", x[k], aniso[n], key->x[w][i]);
	      x[k] += aniso[n++] * key->x[w][i];
	    }
	    x[k] += aniso[n++] * t; //auf keinen Fall nach vorne holen
	  }
	}
      }
    } else { 
      if (key->grid) {/* grid; with or without time component */
	if (simugrid) { // necessarily stemming from diag matrix, but 
	  // compressed according to zeros in the diagonal
	  for (i=0; i<3; i++) {
	    k = i;
	    for (n=d=0; d<key->timespacedim; d++, k+=3) {
 //	      printf("2nogrid %d %f\n", k, x[k]);
 	      x[k] = 0.0;
	      for(w=0; w<key->timespacedim; w++) 
		x[k] += aniso[n++] * key->x[w][i]; 
	    }
	  }
	} else { // grid but not simugrid
	  double y[MAXDIM]; /* current point within grid, but without
			       anisotropy transformation */
	  int yi[MAXDIM]; /* counter for the current position in the grid */
	  for (w=0; w<key->timespacedim; w++) {y[w]=key->x[w][XSTART]; yi[w]=0;}
	  for (k=0; k<total; ){
//	    printf("%f %f %f %d %d %d\n", y[0], y[1], y[2], yi[0], yi[1], yi[2]);
	    for (n=d=0; d<reduceddim; d++, k++) {
	      x[k] = 0.0;
	      for(w=0; w<key->timespacedim; w++) {
		x[k] += aniso[n++] * y[w];
	      }
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
	}
      } else { /* anisoptropic, no grid, no time component */
	for (k=i=0; i<key->totalpoints; i++)
	  for (n=d=0; d<reduceddim; d++, k++) {
	    x[k] = 0.0;
	    for(w=0; w<key->timespacedim; w++) x[k] += aniso[n++] * key->x[w][i];
	  }
      }
    }
  } else { /* not key->anisotropy, no time component */   
      //     printf("transform2 %d %d %d %d\n", key->anisotropy, key->grid,
//	     reduceddim, key->timespacedim);
    assert(reduceddim == key->timespacedim);
    if (key->grid) {
      assert(simugrid);
      for (k=d=0; d<key->timespacedim; d++) {
	for (i=0; i<3; i++) {
	  x[k++] = key->x[d][i] * aniso[0]; 
	}
      }
    } else {
      for (k=i=0; i<key->totalpoints; i++) {
        for (j=0; j<key->timespacedim; j++) x[k++] = key->x[j][i] * aniso[0];
      }
    }
  }
  return NOERROR;
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

double SndDerivCovFct(double *x, int dim, covinfo_arraytype keycov,
		      covlist_type covlist, int ncov) {
    // only for isotropic models and isotropy 
// considers variance and the scale parameter
  double fct[MAXCOV], abl[MAXCOV], snd[MAXCOV], zw, z, result = 0.0;
  int v = 0, vold, w1, w2, i;
  covinfo_type *kc;
  cov_fct *cov=NULL; /* no meaning -- just avoids warning from -Wall */

  assert(dim==1);
  while (v<ncov) {
    vold = v;
    v++; while ((v<ncov) && keycov[covlist[v-1]].op) v++;
    for (w1=vold; w1<v; w1++) {
      kc = &(keycov[covlist[w1]]);
      cov = &(CovList[kc->nr]);
      z = fabs(x[0] * kc->aniso[0]);
      fct[w1] = kc->param[VARIANCE] * cov->cov(&z, kc->param); 
      abl[w1] = kc->param[VARIANCE] * cov->derivative(&z, kc->param)
	* kc->aniso[0]; 
      snd[w1] = kc->param[VARIANCE] * cov->secondderivt(&z, kc->param)
	  * kc->aniso[0] * kc->aniso[0];
//      printf("2nd: %d %f %f %f %f %s %d\n", w1,fct[w1], abl[w1], snd[w1],
//	     cov->secondderivt(&z, kc->param),  cov->name, cov->secondderivt);
    }
    for (w1=vold; w1<v; w1++) {
      for (w2=w1; w2<v; w2++) {
	zw = (w1==w2 ? snd[w1] : (2.0 * abl[w1] * abl[w2])); 
	for (i=vold; i<v; i++)
	  if (i!=w1 && i!=w2) zw *= fct[i];
	result += zw;
      }
    }
  }
  return result;
}

double DerivCovFct(double *x, int dim, covinfo_arraytype keycov,
		   covlist_type covlist, int ncov) {
    // only for isotropic models and isotropy 
// considers variance and the scale parameter
  double fct[MAXCOV], abl[MAXCOV], zw, z, result = 0.0;
  int v = 0, vold, i, w;
  covinfo_type *kc;
  cov_fct *cov=NULL; /* no meaning -- just avoids warning from -Wall */

  assert(dim==1);
  while (v<ncov) {
    vold = v;
    v++; while ((v<ncov) && keycov[covlist[v-1]].op) v++;
    for (w=vold; w<v; w++) {
      kc = &(keycov[covlist[w]]);
      cov = &(CovList[kc->nr]);
      z = fabs(x[0] * kc->aniso[0]);

//      printf("phi1-C: %f %f %f %f\n",  kc->param[VARIANCE], kc->aniso[0], x[0], z);

      fct[w] = kc->param[VARIANCE] * cov->cov(&z, kc->param); 
      abl[w] = kc->param[VARIANCE] * cov->derivative(&z, kc->param)
	  * kc->aniso[0];
    }
    for (w=vold; w<v; w++) {
      zw = abl[w];
      for (i=vold; i<v; i++) if (i!=w) zw *= fct[i];
      result += zw;
    }
  }
  return result;
}

//int zaehler = 0;
double CovFct(double *x, int dim, covinfo_arraytype keycov,
		 covlist_type covlist, int ncov, bool anisotropy) {
  int v, ani, k, j, endfor, ncovM1,  nsub, subdim;
  double  var,  zz, zw, result, d, z[MAXDIM + 2];
  static double ZERO[MAXDIM] = {0,0,0,0};
  covinfo_type *kc;
  cov_fct *cov=NULL; /* no meaning -- just avoids warning from -Wall */

  z[0] = z[1] = RF_NAN;

  zw = result = 0.0;
  ncovM1 = ncov - 1;
//  printf("ani %d \n", anisotropy);
  if (anisotropy) {
    int cov_type;
    for (v=0; v<ncov; v++) {
      zw = var = 1.0;
      for(; v<ncov; v++) {
	kc = &(keycov[covlist[v]]);
        cov = &(CovList[kc->nr]);
	var *= kc->param[VARIANCE];
	cov_type=cov->type;
	for (k=0; k<kc->reduceddim; k++) z[k]=0.0;
	for (ani=0, k=0; k<kc->reduceddim; k++) {
	  for (j=0; j<dim; j++) {
	    z[k] += kc->aniso[ani++] * x[j];
	  }
	}
//	printf("%d %f %f\n", dim, z[0], z[1]);
	if (cov_type==ANISOTROPIC) 
	    zw *= cov->cov(z, kc->param);
	/* derzeit ist dim in cov->FCT bis auf assert-checks unbenutzt */
	else {
	  /* the following definition allows for calculating the covariance */
	  /* function for spatialdim=0 and SPACEISOTROPIC model, namely as a */
          /* purely temporal model; however, in CheckAndBuildCov such kind of */
          /* calculations are prevend as being TRIVIAL(generally, independent*/
          /* of being SPACEISOTROPIC or not) -- unclear whether relaxation in*/
          /* CheckAndBuildCov possible */
	  if (cov_type != ANISOHYPERMODEL) {
            //ISOTROPIC, ISOHYPERMODEL, spaceisotropic
            zz = 0.0;
//	    printf("red=%d %f %f %f %f %f %f -- ", 
//		   kc->reduceddim, z[0], z[1], z[2], 
//		   kc->aniso[6], kc->aniso[7], kc->aniso[8]);
	    endfor = kc->reduceddim - 1;
	    for (j=0; j<endfor; j++) zz += z[j] * z[j];
	    if (cov_type==ISOTROPIC || cov_type==ISOHYPERMODEL) 
		zz += z[endfor] * z[endfor];
	    else z[1]=z[endfor]; // spaceisotropic
	    z[0] = sqrt(zz);
	  }
	  if (cov_type==ISOHYPERMODEL || cov_type==ANISOHYPERMODEL) {
	    nsub = (int) kc->param[HYPERNR]; 
	    subdim = (int) kc->param[EFFECTIVEDIM];
	    // subdim + 1 : C(0)
	    // subdim : C(x)
//	   printf("A %d %f\n", subdim, z[0]);
	    z[subdim] =CovFct(x, dim, keycov, &(covlist[v+1]), nsub, anisotropy);
//	   printf("B\n");
	    z[subdim+1] = CovFct(ZERO, dim, keycov, &(covlist[v+1]), nsub, 
				 anisotropy);
	    // printf("%f %f %f %f %d \n",zw,  z[0], z[1], z[2], subdim);
//	    printf("zw=%f ", zw);
	    zw *= cov->cov(z, kc->param);
//	    printf("%f %f %f var=%f %f %f %f subd=%d %s %f %f\n", 
//		   z[0], z[1], z[2], kc->param[VARIANCE], 
//		   kc->param[KAPPA1], kc->param[KAPPA2], kc->param[KAPPA3],
//		   subdim, cov->name, cov->cov(z, kc->param, subdim), zw);
//
//	    assert(zw < 100000000.00); // nicht loeschen
	    v += nsub;
	    kc = &(keycov[covlist[v]]);
	  } else {
	    zw *= cov->cov(z, kc->param);
	  }
	}
	if ((v<ncovM1) && (kc->op==0)) break;
      }
      result += var * zw; // 
    }
  } else {
    if (dim==1) d=fabs(x[0]);
    else {
      for (d=0.0, j=0; j<dim; j++) {d += x[j] * x[j];}
      d = sqrt(d);
    }
    for (v=0; v<ncov; v++) {
      zw = var = 1.0;
      for(; v<ncov; v++) {
        //printf("covfct %d %d\n", v, covlist[v]);

	kc = &(keycov[covlist[v]]);
//	printf("covfct %d %d %s %f\n", v, kc->nr, CovList[kc->nr].name,
//	       kc->aniso[0]); 
        cov = &(CovList[kc->nr]); /* cov needed in FORMULA for Vario */
	if (cov->type==ISOHYPERMODEL) {
	   nsub = (int) kc->param[HYPERNR];
           z[0] = d * kc->aniso[0];
	   // subdim + 1 : C(0)
	   // subdim : C(x)
	   z[1] = CovFct(x, dim, keycov, &(covlist[v+1]), nsub,
			      anisotropy);
	   z[2] = CovFct(ZERO, dim, keycov, &(covlist[v+1]), nsub,
				 anisotropy);
//	   if (zaehler < 100) printf("-----   %f %f %f %f\n", z[0], z[1],
//				     d, kc->aniso[0]);
	   zw *= cov->cov(z, kc->param);
	   v += nsub;
	   kc = &(keycov[covlist[v]]);
	} else {	
          zz = d * kc->aniso[0];
	  var *= kc->param[VARIANCE];
	  zw *= cov->cov(&zz, kc->param);
	}
	if ((v<ncovM1) && (kc->op==0)) break;
      }
      result += var * zw; //
//      if ( zaehler++ < 20)   
//	printf("-->> %f %f %f %f %f %d %d %d\n",   result, var, zw, zz,
//	       kc->param[KAPPA], v, covlist[v], kc->nr);
    }
  }
//  if (  zaehler++ < 100)
  //   printf(" %f (%f, %f) %f d=%f %f %d %d %d -- zz=%f, var=%f zw=%f\n", 
//	   *x, z[0], z[1], result, d, kc->aniso[0],
//	   ncov, keycov[0].nr, keycov[ncov==1 ? 0 : 1].nr
//	   ,zz, var, zw
//      );

//  if (dim==2) 
//    printf("X%d %f %f, %f %f: %f %f \n",dim, x[0], x[1], z[0], z[1], zw, result);
//  else 
//    printf("-%d %f, %f : %f %f \n",dim, x[0], z[0], zw, result);
// assert(result > 0.4);

  return result;
}

char STR[30];
char *strround(double x){
  if (x == floor(x + 0.5) ) sprintf(STR, "%d", (int) x);
  else sprintf(STR, "%f", x);
  return STR;
}
void addleft(char *str, char *sign, double x) {
  sprintf(str, "%s%s%s", str, strround(x), sign);
}
void addright(char *str, char *sign, double x) {
  sprintf(str, "%s%s%s", str, sign, strround(x));
}
void addone(char *str, double x) {
  sprintf(str, "%s, %s", str, strround(x));
}
void addpair(char *str, double x, double y) {
    if (x == floor(x + 0.5) && y==floor(y + 0.5))
	sprintf(str, "%s, (%d,%d)", str, (int) x, int(y));
    else 
	sprintf(str, "%s, (%f,%f)", str, x, y);
}

int check_within_range(param_type param, cov_fct *cov, int timespacedim, 
		       char *txt) {
  int error, index, i4, i, kappas, maxdim;
  double range[LASTKAPPA * 4];
  bool truesegm;
  rangefct getrange;

  if (GENERAL_SKIPCHECKS) return NOERROR;

  getrange = cov->range;
  kappas = cov->kappas(timespacedim); 
  index = 1;
  for(;;) {
    getrange(timespacedim, &index, range); 
    if (index==RANGE_INVALIDDIM) {
      maxdim = timespacedim;
      while (index==RANGE_INVALIDDIM) { // check when dim starts being valid
	maxdim--;
	assert(maxdim > 0);
	index = 1;
	getrange(maxdim, &index, range); // range not of interest anymore
      }
      sprintf(ERRORSTRING_WRONG, "genuine dimension = %d%s", 
	      timespacedim, txt == NULL ? "" : txt);
      sprintf(ERRORSTRING_OK, "dim <= %d", maxdim);
      return ERRORCOVFAILED;  
    }
    error = -1;
    truesegm = true;
    for (i=0; i<kappas; i++) {
      i4 = i * 4;
      if (range[i4] == floor(range[i4]) + OPEN) {
	  if (param[KAPPA1 + i] <= range[i4] - OPEN) error=i;
      } else {
	  if (param[KAPPA1 + i] < range[i4]) error=i;
      }
      i4++;
      if (range[i4] + OPEN == ceil(range[i4])) {
	  if (param[KAPPA1 + i] >= range[i4] + OPEN) error=i;
      } else {
	  if (param[KAPPA1 + i] > range[i4]) error=i;
      }	   
      truesegm &= range[i4] != range[i4-1] || param[KAPPA1 + i] == range[i4];
      
/*
// printf("%d %d %d %d %d %f %f %f %d\n", i,
	param[KAPPA1 + i] < range[--i4],
	param[KAPPA1 + i] == range[i4] &&range[i4] == floor(range[i4]) + OPEN ,
	param[KAPPA1 + i] > range[++i4] ,
	param[KAPPA1 + i] == range[i4] && range[i4] + OPEN == ceil(range[i4]),
	range[i4-1], range[i4], param[KAPPA1 + i], truesegm
	);
*/
    }
    if (truesegm) { // exakte indizierung darf nur 1x auftauchen 
      //               kleine Einschr\"ankung, aber bisher kein Problem
      if (error>=0) {
	i = error;
	i4 = i * 4;
	sprintf(ERRORSTRING_WRONG, "kappa%d=%f", i+1, param[KAPPA1 + i]);
	strcpy(ERRORSTRING_OK, "");
	if (range[i4] > RF_NEGINF) {
	  if (range[i4] == floor(range[i4]) + OPEN)
	    addleft(ERRORSTRING_OK, "<", floor(range[i4]));
	  else
	    addleft(ERRORSTRING_OK, "<=", range[i4]);
	}
	i4++;
	sprintf(ERRORSTRING_OK, "%skappa%d", ERRORSTRING_OK, i+1);
	if (range[i4] < RF_INF) {

	  printf("%f %f %d %e\n", range[i4] + OPEN, ceil(range[i4]),
		 range[i4] + OPEN == ceil(range[i4]),
		 range[i4] + OPEN - ceil(range[i4]));

	  if (range[i4] + OPEN == ceil(range[i4]))
	    addright(ERRORSTRING_OK, "<", ceil(range[i4]));
	  else
	    addright(ERRORSTRING_OK, "<=", range[i4]);
	}
	sprintf(ERRORSTRING_OK, "%s when %s dim=%d%s", 
		ERRORSTRING_OK,
		(cov->type == SPACEISOTROPIC) ? "genuine spatial" :
		(cov->type == ISOTROPIC || cov->type==ISOHYPERMODEL) 
		? "genuine" : "",
		timespacedim - (int) (cov->type == SPACEISOTROPIC),
		txt == NULL ? "" : txt);
//	  printf("%e %e\n", range[i4-1], range[i4]);
//	assert(false);
	return ERRORCOVFAILED;
      } else return NOERROR;
    } else {
      if (index==RANGE_LASTELEMENT) {
	int n, idx[LASTKAPPA];
	n = 0;
	for (i4=i=0; i<kappas; i++, i4+=4) 
	  if (range[i4] == range[i4+1]) {
	    idx[n] = i;
	    n++;
	  }
	index = 1;
	assert(n>0);
	if (n==1) {
	  sprintf(ERRORSTRING_OK, "kappa%d =", idx[0]+1);
	} else {
	  sprintf(ERRORSTRING_OK,"(kappa%d, kappa%d) =", idx[0]+1, idx[1]+1);
	}
	while (index > 0) {
	  cov->range(timespacedim, &index, range);
	  if (n==1) {
	    addone(ERRORSTRING_OK, range[idx[0] * 4]);
	  } else {
	    addpair(ERRORSTRING_OK, range[idx[0] * 4], range[idx[1] *4]);
	  }
	}
	
	if (n==1) {
	  sprintf(ERRORSTRING_WRONG, "kappa%d =", idx[0]+1);
	  addone(ERRORSTRING_WRONG, param[KAPPA1]);
	} else {
	  sprintf(ERRORSTRING_WRONG,"(kappa%d, kappa%d) =", 
		  idx[0]+1, idx[1]+1);
	  addpair(ERRORSTRING_WRONG,
		  param[KAPPA1+idx[0]], param[KAPPA1+idx[1]]);
	}
	return ERRORCOVFAILED;  
      }
    }
  }
}

// checks a *single* basic covariance function
// p and param must be of the same kind (double)!!! 
int CheckAndBuildCov(int *covnr, int *op, int ncov,
		     double *ParamList, int nParam,
		     int naturalscaling, bool Time, bool grid,
		     int timespacedim, bool anisotropy, 
		     covinfo_arraytype keycov, bool *equal)
{ 
  // IMPORTANT!: p is the parameter vector of the R interface !!
  // output param : equals p, except naturalscaling
  int error, v, i, totkappas, pAniso, d;
  double newscale;
  covinfo_type *kc;
  cov_fct *cov;
  param_type param;

  if (currentNrCov==-1) InitModelList(); assert(CovList!=NULL); 
  if (ncov<0 || ncov>MAXCOV) return ERRORNCOVOUTOFRANGE;
  if (timespacedim < 1 || timespacedim > MAXDIM) return ERRORDIM; 
  pAniso = anisotropy ? timespacedim * timespacedim : 1;
  if (Time && !anisotropy) return ERRORTIMENOTANISO; 

  totkappas = 0;
  for (v=0; v<ncov; v++){
    ERRORMODELNUMBER = v;	
    if (covnr[v]<0 || covnr[v]>=currentNrCov) return ERRORCOVNROUTOFRANGE;     
    *equal &= keycov[v].nr == covnr[v];
    keycov[v].nr = covnr[v];
    kc = &(keycov[v]);
    kc->dim = timespacedim;
    cov = &(CovList[kc->nr]);
    if (GENERAL_PRINTLEVEL > 7)
       PRINTF("Check&Build %d %d covnr=%d op=%d\n", v, ncov, covnr[v],
	      (v<ncov-1) ? op[v] : -1);

    // Variance and Kappa
    int kappas, kappasP1;
    kappas = cov->kappas(timespacedim);
    kappasP1 = kappas + 1;
    *equal &= !memcmp(kc->param, ParamList, sizeof(double) * kappasP1);
    memcpy(kc->param, ParamList, sizeof(double)* kappasP1);
    ParamList += kappasP1;
    totkappas += kappas;
    if (kc->param[VARIANCE]<0.0) return ERRORNEGATIVEVAR;
    switch(cov->type) {
	case ISOTROPIC : case ISOHYPERMODEL :
	    kc->param[EFFECTIVEDIM] = 1.0;
	    break;
	case SPACEISOTROPIC :
	    kc->param[EFFECTIVEDIM] = 2.0;    
	    if (!Time) return ERRORWITHOUTTIME;
	    break;
	case ANISOTROPIC : case ANISOHYPERMODEL :
	    kc->param[EFFECTIVEDIM] = (double) timespacedim;
	    break;
    }

    // op
    if (v < ncov - 1) {
      *equal &= kc->op == op[v];
      kc->op = op[v];
      //     printf(" type %d %d %d %d %d %d %d\n", v, cov->type, ISOHYPERMODEL,
//	      kc->op, false xor false, false xor true, true xor true);
      if ((cov->type == ISOHYPERMODEL || cov->type == ANISOHYPERMODEL)
	  xor (kc->op == 2))
	  return ERRORPARENTHESIS;
    }

    // not hypermodel: ?
    // natural scaling
    GetNaturalScaling(&(kc->nr), &(kc->param[KAPPA]),
		      &naturalscaling, &newscale, &error);
//  printf("&naturalscaling = %d \n", naturalscaling);
    if (error != NOERROR) return error;
    if (anisotropy) newscale = 1.0 / newscale;
    for (i=pAniso - 1; i>=0; i--)
	param[ANISO + i] = newscale * ParamList[i];
    if (!anisotropy) {
      if (cov->type==SPACEISOTROPIC || cov->type==ANISOTROPIC || 
	  cov->type==ANISOHYPERMODEL)
	  return ERRORNOTANISO;
      if (cov->cov==nugget && param[SCALE]==0.0) param[SCALE] = 1.0;
      if (param[SCALE] <= 0.0) return ERRORNEGATIVESCALE;
    }
    ParamList += pAniso;
    *equal &= !memcmp(&(kc->param[ANISO]), &(param[ANISO]), 
		      sizeof(double)* pAniso);
    memcpy(&(kc->param[ANISO]), &(param[ANISO]), sizeof(double) * pAniso);
    
    // reduced anisotropy matrix and true grid, time and dimension
///////
    // for hypermodels it is very unclear, how the dimension should be taken
    GetTrueDim(anisotropy, timespacedim, kc->param, cov->type, 
	       &(kc->genuine_last_dimension), 
	       &(kc->reduceddim), kc->aniso);

    if (kc->reduceddim==0 || // e.g. used in RFnugget and Spectral
	(cov->type==SPACEISOTROPIC && 
	 (!kc->genuine_last_dimension || kc->reduceddim<2))
	) return ERRORLOWRANK;

    // check value of parameters
    error = check_within_range(kc->param, cov, kc->reduceddim, NULL);
    if (error != NOERROR) return error;

    int maxdim, dummy;
    cov->info(kc->param, &maxdim, &dummy);
//    printf("%d %d\n", maxdim, kc->reduceddim); 
    if (maxdim < kc->reduceddim && !GENERAL_SKIPCHECKS)
	return ERRORCOVNOTALLOWED;
    kc->simugrid = // kc->param[ANISO] !!
	grid && (!anisotropy || is_diag(&(kc->param[ANISO]), timespacedim));
    if (kc->simugrid) {
      if (anisotropy) {
	int n;
	for (i=n=d=0; d<timespacedim; d++, n += timespacedim + 1) {
	    if ((kc->genuine_dim[d] = param[ANISO + n] != 0.0)) {
	      // *** funktioniert nur deshalb weil nachgeprueft 
	      // worden ist das diag matrix ***
	      kc->idx[d] = i++;
	      //  WICHTIG! Dadurch wird Ordnung
	      //  der Dimensionen eingehalten!!
	      // siehe auch Transform2nogrid
	    } else {
	      kc->idx[d] = kc->reduceddim + d - i; 
	    }
	}
      } else {
	for (d=0; d<timespacedim; d++) {
	  kc->genuine_dim[d] = true;
	  kc->idx[d] = d;
	}
      }
    }
    if (cov->type == ISOHYPERMODEL || cov->type == ANISOHYPERMODEL) {
////////////////////////
	kc->simugrid = false; // not programmed yet -- save version !!!!!!
///////////////////////
    } else {
      error = cov->check(kc->param, kc->reduceddim, kc->dim, Nothing);
      if (error != NOERROR) return error; 
    }
  }
  ERRORMODELNUMBER = -1;

  for (v=ncov-1; v>=0; v--) {
    int w;
    ERRORMODELNUMBER = v;	
    // this loop cannot be intregrated above. The subsequent
    // submodels must be initialised first, since  
    //   1) checkNinit relies on the knowledge of the submodels
    //   2) parameters are taken over from the first submodel
    // this loop must be in reverse order for the same reasons.
    kc = &(keycov[v]);
    cov = &(CovList[kc->nr]);
    if (cov->type == ISOHYPERMODEL || cov->type == ANISOHYPERMODEL) {
      // take over the anisotropy structure from the first submodel
// 28.9.05: not a good idea !!!
//      memcpy(&(kc->param[ANISO]), &(keycov[v+1].param[ANISO]),
//	     sizeof(double) * pAniso);

/*
  memcpy(&(kc->param[ANISO]), &(keycov[v+1].param[ANISO]),
	     sizeof(double) * pAniso);
     GetTrueDim(anisotropy, timespacedim, kc->param,
		 CovList[keycov[v+1].nr].type, 
		&(kc->genuine_last_dimension), //
		&(kc->reduceddim), kc->aniso);//
     kc->aniso[0] = RF_NAN;
*/

      for (w = v + (int) (kc->param[HYPERNR]); w>v; w--) {
	 if (CovList[keycov[w].nr].type != ISOTROPIC) {
	     // z.B. tbm viel complizierter
	     error = ERRORHYPERNOTISO;
	 }
      }
      error = cov->checkNinit(kc, &(COVLISTALL[v]), ncov-v-1, Nothing,
			      anisotropy);
      if (error != NOERROR) {
	  return error;
      }
    }
  }
  ERRORMODELNUMBER = -1;
	
  if (ncov * (1 + pAniso) + totkappas != nParam) return ERRORPARAMNUMBER;
  return NOERROR;
}


static int user_dim=-1, user_ncov=-1;
static bool user_anisotropy=false;
static covinfo_arraytype user_keycov;
void CheckAndCompleteParameters(int *covnr, double *p, int *np, 
				int *dim, /* timespacedim ! */
				int *ncov, int *anisotropy, int *op,
				double *q, int* error) {
  int v;
  bool dummyequal = false;
  // note: column of x is point !!
  COVINFO_NULL(user_keycov);
  user_dim = *dim;
  user_ncov=*ncov;
  user_anisotropy=(bool) *anisotropy;
  *error = CheckAndBuildCov(covnr, op, user_ncov, p, *np, GENERAL_NATURALSCALING,
			    user_anisotropy /* Time */, false /* grid*/,
			    user_dim, user_anisotropy, user_keycov, 
			    &dummyequal);
  if (*error != NOERROR) {
      if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing, *error);
      return;
  }
  for (v=0; v<user_ncov; v++, q += TOTAL_PARAM) {
      memcpy(q, user_keycov[v].param, sizeof(param_type));
  }
}


int CheckCovariance(int *covnr, double *p, int np, 
		    int logicaldim, /* timespacedim ! */
		    int xdim, int ncov, bool anisotropy, int *op,
		    int naturalscaling, covinfo_arraytype keycov)
{
  int v, error;
  bool dummyequal = false;
  cov_fct *cov;

  if ((logicaldim!=xdim) && (anisotropy || (xdim!=1))) return ERRORDIMMISMATCH;
  error = CheckAndBuildCov(covnr, op, ncov, p, np, naturalscaling, 
			   anisotropy /* Time */,  false /* grid*/,
			  logicaldim, anisotropy, keycov, &dummyequal);
  if (error != NOERROR) return error;

  for (v=0; v<ncov; v++) {
    ERRORMODELNUMBER = v;
    covinfo_type *kc;
    kc = &(keycov[v]);
    cov = &(CovList[kc->nr]);
    if (cov->cov==NULL) return ERRORNOTPROGRAMMED;
    if (cov->type==ISOHYPERMODEL || cov->type==ANISOHYPERMODEL) {
	v += (int) kc->param[HYPERNR];
    } else {
	if (cov->variogram) return ERRORNOTDEFINED; 
    }
  }
  ERRORMODELNUMBER = -1;
  return NOERROR;
}

void CalculateCovariance(double *x, int lx, covinfo_arraytype keycov, int xdim,
		    int ncov, bool anisotropy, double *result) {
  int i;
  for (i=0; i<lx; i++, x += xdim)
    result[i] = CovFct(x, xdim, keycov, COVLISTALL, ncov, anisotropy);
}


// note: number of parameters is not checked whether it is correct!!
// done by R function `PrepareModel'
void CovarianceNatSc(double *x, int *lx, int *covnr, double *p, int *np, 
		     int *logicaldim, /* timespacedim ! */
		     int *xdim,
		     int *ncov, int *anisotropy, int *op,
		     double *result, int *naturalscaling)
{
  int error, i;
  // note: column of x is point !!
  COVINFO_NULL(user_keycov);
  user_dim = *logicaldim;
  user_ncov=*ncov;
  user_anisotropy=(bool) *anisotropy;
  if ((error = CheckCovariance(covnr, p, *np, user_dim, *xdim, user_ncov,
			       user_anisotropy,
			       op, *naturalscaling, user_keycov)) != NOERROR) {
    if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing, error);
    for (i=0; i<*lx; i++) result[i] = RF_NAN;
    return;	
  }
  // die nachfolgende Verschraenkung von logicaldim und xdim
  // funktioniert, siehe code von CheckAndBuildCov und GetTrueDim
  CalculateCovariance(x, *lx, user_keycov, *xdim, user_ncov, user_anisotropy, 
		      result);
}

void Covariance(double *x,int *lx, int *covnr, double *p, int *np, 
		int *logicaldim, /* timespacedim ! */
		int *xdim, int *ncov, int *anisotropy, int *op,
		double *result)
{

  CovarianceNatSc(x, lx, covnr, p, np, logicaldim, xdim,
		  ncov, anisotropy, op, result, &GENERAL_NATURALSCALING);
}

void CalculateCovMatrix(double *dist, int lx, covinfo_arraytype keycov,
			int xdim, int ncov, bool anisotropy, double *result) {
  // note: column of x is point !!
  int i, ii, endfor, lxP1;
  long ve, ho;
  double var;

  var = CovFct(ZERO, xdim, keycov, COVLISTALL, ncov, anisotropy);
  lxP1 = lx + 1;
  for (ii=lx, i=0; ii>0; i+=lxP1, ii--) {
    result[i] = var;
    endfor = i + ii;
    for (ve = i + 1, ho = i + lx; ve < endfor; ve++, ho += lx, dist += xdim) {
      result[ve] = result[ho] = 
	  CovFct(dist, xdim, keycov, COVLISTALL, ncov, anisotropy);
    }
  }
}

void CovarianceMatrixNatSc(double *dist,int *lx, int *covnr,double *p, 
			   int *np, int *logicaldim, /* timespacedim ! */
			   int *xdim, int *ncov, int *anisotropy, int *op,
			   double *result, int *naturalscaling)
{
  // note: column of x is point !!
  int error;
  long i, lxq;
 
  COVINFO_NULL(user_keycov);
  user_dim = *logicaldim;
  user_ncov=*ncov;
  user_anisotropy=(bool) *anisotropy;
  error = CheckCovariance(covnr, p, *np, user_dim, *xdim, user_ncov, 
			  user_anisotropy, op, *naturalscaling, user_keycov);
  if (error != NOERROR) {
      if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing, error);
      lxq = *lx * *lx;
      for (i=0; i<lxq; i++) result[i] = RF_NAN;
      return;	
  }
  CalculateCovMatrix(dist, *lx, user_keycov, *xdim, user_ncov, 
		     user_anisotropy, result);
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


static int Unchecked_dim=-1, Unchecked_ncov=-1;
static bool Unchecked_anisotropy=false;
static covinfo_arraytype Unchecked_keycov;
void InitUncheckedCovFct(int *covnr, double *p, int *np, 
			 int *logicaldim, /* timespacedim ! */
			 int *xdim, int *ncov, int *anisotropy, int *op,
			 int *naturalscaling, int *error)
{
  if (Unchecked_dim==-1) COVINFO_NULL(Unchecked_keycov);
  if ((*error = CheckCovariance(covnr, p, *np, *logicaldim, *xdim, *ncov, 
				(bool) *anisotropy, op, *naturalscaling, 
				Unchecked_keycov)) != NOERROR) {
    if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing, *error);
    return;	
  }
  Unchecked_ncov = *ncov;
  Unchecked_anisotropy = (bool) *anisotropy;
  Unchecked_dim = *xdim;
  return;
}

void UncheckedCovFct(double *x, int *lx, double *result)
{
  CalculateCovariance(x, *lx, Unchecked_keycov, Unchecked_dim, Unchecked_ncov,
		      Unchecked_anisotropy, result);
}

void UncheckedCovMatrix(double *dist, int *lx, double *result)
// corresponds to CovarianceMatrixNatSc
{
  CalculateCovMatrix(dist, *lx, Unchecked_keycov, Unchecked_dim, Unchecked_ncov,
		     Unchecked_anisotropy, result);
}


int CheckVariogram(int *covnr, double *p, int np, 
		    int logicaldim, /* timespacedim ! */
		    int xdim, int ncov, bool anisotropy, int *op,
		    int naturalscaling, covinfo_arraytype keycov)
{
  int v, error;
  bool dummyequal = false;
  cov_fct *cov;

  if ((logicaldim!=xdim) && (anisotropy || (xdim!=1))) return ERRORDIMMISMATCH;

//  printf("Cv %d %d %d\n",  logicaldim, xdim, anisotropy);

  error = CheckAndBuildCov(covnr, op, ncov, p, np, naturalscaling, 
			  anisotropy, /* Time */ false, /* grid*/
			  logicaldim, anisotropy, keycov, &dummyequal);
  if (error != NOERROR) return error;

  for (v=0; v<ncov; v++) {
    ERRORMODELNUMBER = v;
    cov = &(CovList[covnr[v]]);
    if (cov->cov==NULL) return ERRORNOTPROGRAMMED;
    // below is too restrictive since there could be a hyper model
    // with a multiplicative model as submodel
    if (v>0 && op[v-1] && (cov->variogram || CovList[covnr[v-1]].variogram))
	return ERRORNOMULTIPLICATION;
  }
  ERRORMODELNUMBER = -1;
  return NOERROR;
}

void CalculateVariogram(double *x, int lx, covinfo_arraytype keycov, int xdim,
			int ncov, bool anisotropy, double *result) {
  int i;
  double C0;
  C0 = CovFct(ZERO, xdim, keycov, COVLISTALL, ncov, anisotropy);
  for (i=0; i<lx; i++, x += xdim)
    result[i] = C0 - CovFct(x, xdim, keycov, COVLISTALL, ncov, anisotropy);
}

void VariogramNatSc(double *x,int *lx,int *covnr,double *p,int *np,
		    int *logicaldim, /* timespacedim ! */
		    int *xdim,
		    int *ncov, int *anisotropy, int *op,
		    double *result, int *naturalscaling)
{
  int error, i;

// printf("Cv %d %d %d\n",  *logicaldim, *xdim, *anisotropy);
  COVINFO_NULL(user_keycov);
  user_dim = *logicaldim;
  user_ncov=*ncov;
  user_anisotropy=(bool) *anisotropy;
  if ((error = CheckVariogram(covnr, p, *np, user_dim, *xdim, user_ncov,
			      user_anisotropy,
			      op, *naturalscaling, user_keycov)) != NOERROR) {
    if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing, error);
    for (i=0; i<*lx; i++) result[i] = RF_NAN;
    return;	
  }
  CalculateVariogram(x, *lx, user_keycov, *xdim, user_ncov, user_anisotropy, 
		     result);
}
  
void Variogram(double *x,int *lx,int *covnr,double *p,int *np, 
	       int *logicaldim, /* timespacedim ! */
	       int *xdim,
	       int *ncov, int *anisotropy, int *op,
	       double *result)
{
  VariogramNatSc(x, lx, covnr, p, np, logicaldim, xdim, ncov, anisotropy, op,
		 result, &GENERAL_NATURALSCALING);
}

void VariogramMatrix(double *dist, int *lx, int *covnr, double *p, int *np, 
		     double *result)
{ 
  // not programmed yet
  assert(false);
}

/*
int insert_method(int* last_incompatible, key_type *key,
			 SimulationType m, bool incompatible) {
  int idx;

  assert(m != Nugget || !incompatible);

  assert(key->meth[key->n_unimeth].S==NULL);
  assert(key->meth[key->n_unimeth].destruct==NULL);
  if (incompatible) {
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

    key->meth[key->n_unimeth].unimeth_alreadyInUse =
      key->meth[*last_incompatible].unimeth_alreadyInUse;
    key->meth[key->n_unimeth].unimeth=key->meth[*last_incompatible].unimeth;
    key->meth[key->n_unimeth].S=key->meth[*last_incompatible].S; // necessary
    // only if further methods are tried due to partial failure
    key->meth[key->n_unimeth].destruct=key->meth[*last_incompatible].destruct;
    key->meth[key->n_unimeth].incompatible = 
	key->meth[*last_incompatible].incompatible;
    key->meth[*last_incompatible].S = NULL;
    key->meth[*last_incompatible].destruct = NULL;
    idx = *last_incompatible;
  } else idx=key->n_unimeth;
  key->meth[idx].unimeth = m;
  key->meth[idx].incompatible = incompatible;
  key->meth[idx].unimeth_alreadyInUse = false;
  key->n_unimeth++;
  return idx;
}

void delete_method(int *M, int *last_incompatible, key_type *key) {
  // note that key->n_unimeth had been decreased before this function
  // is called
  int idx;
  if (key->meth[*M].destruct != NULL) key->meth[*M].destruct(&(key->meth[*M].S));
  assert(key->meth[*M].S==NULL);
  key->meth[*M].destruct=NULL;
  if (key->meth[*M].incompatible) {
    assert(*last_incompatible>=0);
    assert(*M<=*last_incompatible);
    key->meth[*M].unimeth = key->meth[*last_incompatible].unimeth;
    key->meth[*M].unimeth_alreadyInUse =
	key->meth[*last_incompatible].unimeth_alreadyInUse;
    key->meth[*M].S = key->meth[*last_incompatible].S;
    key->meth[*last_incompatible].S=NULL;
    key->meth[*M].incompatible = key->meth[*last_incompatible].incompatible;
    key->meth[*M].destruct=key->meth[*last_incompatible].destruct;
    key->meth[*last_incompatible].destruct=NULL;
    idx = *last_incompatible;
    (*last_incompatible)--;
  } else idx = *M;
  key->n_unimeth--;
  if (idx==key->n_unimeth){ // deleting of the very last entry
    assert(*M==key->n_unimeth || *last_incompatible + 1 == key->n_unimeth); 
  } else {
    assert(idx>=0 && idx<key->n_unimeth);
    key->meth[idx].unimeth = key->meth[key->n_unimeth].unimeth;
    key->meth[idx].unimeth_alreadyInUse = 
	key->meth[key->n_unimeth].unimeth_alreadyInUse;
    key->meth[idx].S = key->meth[key->n_unimeth].S;
    key->meth[key->n_unimeth].S=NULL;
    key->meth[idx].incompatible = key->meth[key->n_unimeth].incompatible;
    key->meth[idx].destruct = key->meth[key->n_unimeth].destruct;
    key->meth[key->n_unimeth].destruct=NULL;
  }
  key->meth[key->n_unimeth].unimeth_alreadyInUse = false;
  (*M)--;
}
*/

int insert_method(int* last_incompatible, key_type *key,
			 SimulationType m, bool incompatible) {
  int idx;

  assert(m != Nugget || !incompatible);

  if (key->meth[key->n_unimeth].destruct != NULL)
    key->meth[key->n_unimeth].destruct(&(key->meth[key->n_unimeth].S));
  assert(key->meth[key->n_unimeth].S==NULL);
  key->meth[key->n_unimeth].destruct=NULL;

  if (incompatible) {
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

    key->meth[key->n_unimeth].unimeth_alreadyInUse =
      key->meth[*last_incompatible].unimeth_alreadyInUse;
    key->meth[key->n_unimeth].unimeth=key->meth[*last_incompatible].unimeth;
    key->meth[key->n_unimeth].S=key->meth[*last_incompatible].S; // necessary
    // only if further methods are tried due to partial failure
    key->meth[key->n_unimeth].destruct=key->meth[*last_incompatible].destruct;
    key->meth[key->n_unimeth].incompatible = 
	key->meth[*last_incompatible].incompatible;
    key->meth[*last_incompatible].S = NULL;
    key->meth[*last_incompatible].destruct = NULL;
    idx = *last_incompatible;
  } else idx=key->n_unimeth;
  key->meth[idx].unimeth = m;
  key->meth[idx].incompatible = incompatible;
  key->meth[idx].unimeth_alreadyInUse = false;
  key->n_unimeth++;
  return idx;
}

void delete_method(int *M, int *last_incompatible, key_type *key) {
  // note that key->n_unimeth had been decreased before this function
  // is called
  int idx, i;
  methodvalue_type dummy;

  dummy.S = key->meth[*M].S;
  dummy.destruct= key->meth[*M].destruct;
  dummy.unimeth = key->meth[*M].unimeth;
  dummy.incompatible = key->meth[*M].incompatible;
  if (key->meth[*M].incompatible) {
    assert(*last_incompatible>=0);
    assert(*M<=*last_incompatible);
    key->meth[*M].unimeth = key->meth[*last_incompatible].unimeth;
    key->meth[*M].unimeth_alreadyInUse =
	key->meth[*last_incompatible].unimeth_alreadyInUse;
    key->meth[*M].S = key->meth[*last_incompatible].S;
    key->meth[*last_incompatible].S=NULL;
    key->meth[*M].incompatible = key->meth[*last_incompatible].incompatible;
    key->meth[*M].destruct=key->meth[*last_incompatible].destruct;
    key->meth[*last_incompatible].destruct=NULL;
    idx = *last_incompatible;
    (*last_incompatible)--;
  } else idx = *M;
  key->n_unimeth--;
  if (idx==key->n_unimeth){ // deleting of the very last entry
    assert(*M==key->n_unimeth || *last_incompatible + 1 == key->n_unimeth); 
  } else {
    assert(idx>=0 && idx<key->n_unimeth);
    key->meth[idx].unimeth = key->meth[key->n_unimeth].unimeth;
    key->meth[idx].unimeth_alreadyInUse = 
	key->meth[key->n_unimeth].unimeth_alreadyInUse;
    key->meth[idx].S = key->meth[key->n_unimeth].S;
    key->meth[key->n_unimeth].S=NULL;
    key->meth[idx].incompatible = key->meth[key->n_unimeth].incompatible;
    key->meth[idx].destruct = key->meth[key->n_unimeth].destruct;
    key->meth[key->n_unimeth].destruct=NULL;
  }

  key->meth[key->n_unimeth].unimeth_alreadyInUse = false;
  key->meth[key->n_unimeth].S = dummy.S;
  key->meth[key->n_unimeth].destruct = dummy.destruct;
  key->meth[key->n_unimeth].unimeth = dummy.unimeth;
  key->meth[key->n_unimeth].incompatible = dummy.incompatible;
  key->meth[key->n_unimeth].actcov = key->meth[*M].actcov;
  for (i=key->meth[*M].actcov-1; i>=0; i--)
      key->meth[key->n_unimeth].covlist[i] = key->meth[*M].covlist[i];
  (*M)--;
}

/*

speed, exactness, stationarity

direct (lower -- always, upper -- never)


grid: circembed
nongrid: spectralTBM, TBM2, TBM3,


Ausnahmen: dimension beschraenkt

*/

void init_method_list(key_type *key)
{
  unsigned int maxmem=500000000;
  int best_dirct=DIRECTGAUSS_BESTVARIABLES, 
    max_dirct=DIRECTGAUSS_MAXVARIABLES, 
    i, v, nr, nml;
  bool anyFirstDirect = false, anyFirstCE = false;

#define nLastDefault 2
  SimulationType LastDefault[nLastDefault] = {AdditiveMpp, Hyperplane}; 
  SimulationType DirectFst = key->totalpoints <= best_dirct ? Direct : Forbidden;
  SimulationType DirectLst = key->totalpoints <= max_dirct ? Direct : Forbidden;
#define nStandard 6
#define TBM3pos 5
#define TBM2pos 4
  SimulationType Standard[nStandard] = 
    {CircEmbed, CircEmbedIntrinsic, CircEmbedCutoff, SpectralTBM, TBM2, TBM3};

  bool Allowed[SimulationTypeLength], allowed[SimulationTypeLength];
  for (i=0; i<Nothing; i++) Allowed[i] = true; 
  if (DECISION_PARAM.stationary_only==DECISION_TRUE) 
    Allowed[CircEmbedIntrinsic]=false;
  if (!key->grid) Allowed[CircEmbed] = Allowed[CircEmbedCutoff] = 
		    Allowed[CircEmbedIntrinsic] = false;
  if (DECISION_PARAM.exactness==DECISION_TRUE)
    Allowed[TBM2] = Allowed[TBM3] = Allowed[SpectralTBM] = 
      Allowed[AdditiveMpp] = Allowed[Hyperplane] = false;

#define nraw 2 * SimulationTypeLength
  SimulationType raw[nraw];
  for (v=0; v<key->ncov; v++) {
    cov_fct *cov;
    cov = &(CovList[key->cov[v].nr]);
    key->IML[v] = -1;
    // special case: Nugget (part 1)
    if (cov->cov==nugget) continue;
    
    Allowed[TBM2] = Allowed[SpectralTBM] =Allowed[Hyperplane] =
	key->cov[v].reduceddim <= 2;
    if (key->cov[v].reduceddim == 3 && key->grid) {
	Allowed[TBM2] = true;
	Standard[TBM2pos] = TBM3;
	Standard[TBM3pos] = TBM2;
    }
    for (i=0; i<nraw; i++) raw[i] = Forbidden;
    int *implemented=cov->implemented, maxdim, CEbadlybehaved;
    for (i=0; i<Nothing; i++) 
      allowed[i] = Allowed[i] && (implemented[i] == IMPLEMENTED || 
				  implemented[i] == NUM_APPROX);
    if (DECISION_PARAM.stationary_only==DECISION_CASESPEC && 
	!cov->variogram) allowed[CircEmbedIntrinsic] = false;
    if (DECISION_PARAM.exactness != DECISION_FALSE && 
	implemented[i]==NUM_APPROX)
      allowed[TBM2] = false;

    // multiplicative models
    if (v<key->ncov-1 && key->cov[v].op) {
      // nur die Reihenfolge fuer den ersten Faktor wird verwendet -- 
      // der letzte Faktor wird gesetzt wie wenn nicht Faktor (spaeter ignoriert)
      // alle anderen Faktoren werden als Faktoren gesetzt aber Methodenfolge
      //    wiederum ignoriert
      // next_element_of_preference_list achtet auf die anderen Faktoren
      for (i=0; i<Nothing; i++) 
	if (i!=TBM3 && i!=CircEmbed && i!=Direct) allowed[i] = false;
    }

    cov->info(key->cov[v].param, &maxdim, &CEbadlybehaved);
//  printf("timespacedim %d %d %s\n", key->timespacedim, key->cov[v].reduceddim, cov->name);
    assert(key->cov[v].reduceddim <= maxdim);

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
      if (key->cov[v].reduceddim==2 || CEbadlybehaved > 1) {
	raw[nr++] = SpectralTBM;
	if (DECISION_PARAM.exactness==DECISION_FALSE) raw[nr++] = TBM2;
      }
      if (key->cov[v].reduceddim==3 && CEbadlybehaved>1) raw[nr++] = TBM3;
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
	key->ML[v][nml++] = raw[i];
	allowed[raw[i]] = false;
      }
    }
    key->NML[v] = nml;
    switch (raw[0]) {
	case Direct : anyFirstDirect = true; break;
	case CircEmbed : anyFirstCE = true; break;
	default: break;
    }
  }

  // nugget, part 2
  for (v=0; v<key->ncov; v++) if (CovList[key->cov[v].nr].cov==nugget){
    nml = 0;
    if (!key->anisotropy) {
      // CE first considered, since nugget may have a positive effect on
      // the CE simuation.
      if (anyFirstCE) key->ML[v][nml++] = CircEmbed;
      if (anyFirstDirect) key->ML[v][nml++] = Direct;
    }
    key->ML[v][nml++] = Nugget;
    key->NML[v] = nml;
  }
}
  

int next_element_of_method_list(key_type *key)
{
  int v;
  for (v=0; v<key->ncov; v++) {
    if (key->cov[v].left) {
      if (v>0 && key->cov[v-1].op) key->cov[v].method = key->cov[v-1].method;
      else {
	(key->IML[v])++;
//	printf("next %d %d %d %d \n", 
//	       v, key->IML[v], key->NML[v], key->ML[v][key->IML[v]]); 
	if (key->IML[v] >= key->NML[v]) return v + 1;
	key->cov[v].method = key->ML[v][key->IML[v]];
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
{
  // keyNr : label where intermediate results are to be stored;
  //         can be chosen freely between 0 and MAXKEYS-1
  key_type *key;

  key = NULL;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORREGISTER; 
    goto ErrorHandling;
  }
  key = &(KEY[*keyNr]);
  if (x==NULL) {*error=ERRORCOORDINATES; goto ErrorHandling;}
  strcpy(ERROR_LOC, "");

  *error = 
      internal_InitSimulateRF(x, T, *spatialdim, 
			      *lx, (bool) *grid, (bool) *Time, 
			      covnr, ParamList, *nParam,
			      *ncov, (bool) *anisotropy, op, method,
			      *distr, key, 
			      GENERAL_NATURALSCALING, CovFct);
  if (*error != NOERROR) goto ErrorHandling;
  return;
 
  ErrorHandling:
  if (GENERAL_PRINTLEVEL>0) {
    ErrorMessage(Nothing, *error);
    PRINTF("\n");
  }
}

int internal_InitSimulateRF(double *x, double *T, 
			    int spatialdim, /* spatial dim only ! */
			    int lx, bool grid, bool Time,
			    int *covnr, double *ParamList, int nParam,
			    int ncov, bool anisotropy, int *op,
			    int *method, 
			    int distr, /* still unused */
			    key_type *key,
			    int naturalscaling,
			    CovFctType covFct)
{ // grid  : boolean
  // NOTE: if grid and Time then the Time component is treated as if
  //       it was an additional space component
  //       On the other hand, SPACEISOTROPIC covariance function will
  //       always assume that the second calling component is a time component!

//    printf("internal %d %d \n", covnr[0], covnr[1]);

  bool user_defined, simple_definition;
  int i, last_incompatible, d,  v, M, err_occurred, error, timespacedim;
  unsigned long totalBytes;
  int method_used[(int) Nothing + 1];
  // preference lists, distinguished by grid==true/false and dimension
  // lists must end with Nothing!
  SimulationType Merr=Nothing;
  char errorloc_save[nErrorLoc];

  strcpy(errorloc_save, ERROR_LOC);

  timespacedim = spatialdim + (int) Time; // may not set to key->timespacedim,
  //                      since key->timespacedim possibly by DeleteKeyNotTrend
//  printf("timespacedim %d %d %d\n", timespacedim , spatialdim , (int) Time);
  if (spatialdim<1 || timespacedim>MAXDIM){
    error=ERRORDIM; goto ErrorHandling;}  
  user_defined = method[0] >= 0;
  totalBytes =  sizeof(double) * lx * spatialdim;
  key->active &= !key->stop;
  // also checked by R function `PrepareModel'

  if ((error = CheckAndBuildCov(covnr, op, ncov,  ParamList, nParam,
		                naturalscaling, Time, grid,
				timespacedim, anisotropy,
				key->cov, &(key->active))) != NOERROR)
      goto ErrorHandling;


  if (user_defined) {
    for (v=0; v<ncov; v++) {
      SimulationType meth;
      covinfo_type *kc;
      meth = (SimulationType) method[v];
//      printf("method %d\n", method[v]);
      kc = &(key->cov[v]);
      if (CovList[kc->nr].cov == nugget && meth!=Nugget && ncov>1 &&
	  ((meth != CircEmbed && meth != Direct) ||
	   kc->reduceddim < timespacedim)) {
	  if (GENERAL_PRINTLEVEL > 5) 
	      PRINTF("method for nugget (nr. %d) set to nugget (orig=%s)",
		     v, METHODNAMES[meth]);
	      meth = Nugget;
      }
      key->active &= kc->method == meth;
      key->cov[v].method = meth;
      //    printf("user %d %d %d \n", v, method[v], method);
    }
  }


  if (key->active) {   
     assert(key->x[0]!=NULL);
     key->active = 
	 (covFct == key->covFct) &&
	 (key->grid == grid) && 
	 ((key->grid && lx==3) || (!key->grid && key->totalpoints==lx)) &&
	 (key->spatialdim==spatialdim) && 
	 (key->distribution==distr) &&
	 (key->ncov==ncov) &&
	 (key->anisotropy==anisotropy) &&
	 (Time == key->Time && 
	  (!Time || !memcmp(key->T, T, sizeof(double) * 3))) &&
	 (!memcmp(key->x[0], x, totalBytes));
//     if (GENERAL_PRINTLEVEL>5 && !key->active) {
//       PRINTF("failure with previous initialisation at spatialdim=%d grid=%d distrib=%d ncov=%d  aniso=%d time=%d key->T=%d op=%d x=%d\n",
//	      key->spatialdim==spatialdim, key->grid == grid,
//	      key->distribution==distr, key->ncov==ncov, 
//	      key->anisotropy==anisotropy,
//	      Time == key->Time, 
//	      !(Time) || ! memcmp(key->T, T, sizeof(double) * 3),
//	      !memcmp(key->x[0], x, totalBytes));
//     }
   }

  if (key->active) {// detected that all relevant parameters agree, 
    //              no initialization is necessary
    // CHECK:
    for (d = key->timespacedim; d<MAXDIM; d++) {assert(key->x[d]==NULL);} 
    if (GENERAL_PRINTLEVEL>=2) ErrorMessage(Nothing, USEOLDINIT);
    return NOERROR; /* ***** END ***** */ 
  } else {// ! key->active
    DeleteKeyNotTrend(key);
    // if not active, all the parameters have to be 
    //                        stored again
    // save all the parameters if key has not been active 
    key->timespacedim = timespacedim;
    key->spatialdim = spatialdim;
    key->grid= grid;
    key->anisotropy = anisotropy;    
    key->distribution=distr; 
    key->ncov=ncov;
    // coordinates
    if (key->x[0]!=NULL) free(key->x[0]);
    if ((key->x[0]=(double*) malloc(totalBytes))==NULL){
      error=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    memcpy(key->x[0], x, totalBytes);
    for (d=1; d<key->spatialdim; d++) key->x[d]= &(key->x[0][d * lx]);
    for (; d<MAXDIM; d++)  key->x[d]=NULL;
    
    //Time
    if (key->Time = Time) {
      if (!key->anisotropy) { error = ERRORTIMENOTANISO; goto ErrorHandling; }
      memcpy(key->T, T, sizeof(double) * 3);
      if (key->grid) {
	key->x[key->spatialdim] = key->T;
      }
    }
  }

  key->covFct = covFct; // CovFct, CovFctTBM2, CovFctTBM2Num, CovFctTBM3

   // not literally "total"param; there might be gaps in the sequence of the 
   // parameters (if not all seven KAPPA parameters are used)
   key->totalparam = anisotropy 
       ? ANISO + timespacedim * timespacedim
       : SCALE + 1; // 0..SCALE
 

  // Delete stored, intermediate result if there are any
  // folgende Zeilen unschoen, da weiter vorne bereits DeleteKeyNotTrend
  // aufgerufen wurde
//  for (v=0; v<MAXCOV; v++) {
//    if (key->meth[v].destruct!=NULL) {
//      key->meth[v].destruct(&(key->meth[v].S));
//      key->meth[v].destruct = NULL;
//    }
//  }
  
   // minor settings, usually overtaken from previous run, if all the 
  // other parameters are identical
  
  if (key->grid) { 
    if ((lx!=3) || 
	(InternalGetGridSize(key->x,&key->timespacedim,key->length))) {
      error=ERRORCOORDINATES; goto ErrorHandling;
    }
    for (key->spatialtotalpoints=1, d=0; d<key->spatialdim; d++) {
      key->spatialtotalpoints *= key->length[d];
    }
    key->totalpoints = key->spatialtotalpoints;
    if (key->Time) key->totalpoints *= key->length[key->spatialdim];
  } else { // not grid
    key->totalpoints = key->spatialtotalpoints = key->length[0]= lx;
    if (key->Time) {
      int Tdim; double *Tx[MAXDIM];
      Tdim = 1;
      Tx[0] = key->T;
      if (InternalGetGridSize(Tx,&Tdim,&(key->length[key->spatialdim]))) {
	error=ERRORCOORDINATES; goto ErrorHandling;
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
  simple_definition =
      key->ncov==1 ||
      (key->ncov==2 && !key->anisotropy && op[0]==0 &&
       (key->cov[0].method==Nugget || key->cov[1].method==Nugget));
  err_occurred = NOERROR;
  error = ERROROUTOFMETHODLIST;
  key->n_unimeth = 0;
  last_incompatible = -1;
  for (v=0; v<key->ncov; v++) key->cov[v].left=true;
  //           indicates which covariance function are left to be initialized


  // method list must be treated separately since method is integer coded
  // when passed to InitSimulate !
  if (user_defined) {
    //either all or none should be user defined
    for (v=0; v<key->ncov; v++) {
      ERRORMODELNUMBER = v;
      covinfo_type *kc;
      kc = &(key->cov[v]);
      if (kc->method >= (int) Nothing) {
//	printf("kc->method=%d %d %d \n", v, kc->method, (int) Nothing);
	error=ERRORNOTDEFINED; goto ErrorHandling;}
      if (CovList[kc->nr].type == ISOHYPERMODEL || 
	  CovList[kc->nr].type == ANISOHYPERMODEL) {
	  int i, nc = (int) kc->param[HYPERNR];
	  if (nc <= 0 || nc + v >= key->ncov) {
	      error=ERRORHYPERNR; goto ErrorHandling;}
	  for (i=v+1; i<nc+v; i++)
	      if (kc->method != key->cov[i].method) {
		  error=ERRORHYPERMETHOD; goto ErrorHandling;}
      }
      if (v>0 && op[v-1] && kc->method!=key->cov[v-1].method) {
	error=ERRORMETHODMIX; goto ErrorHandling;}
      if (DECISION_PARAM.stationary_only==DECISION_TRUE && 
	  kc->method==CircEmbedIntrinsic) {
	error=ERRORMETHODEXCLUDED; goto ErrorHandling;}
      key->NML[v] = 1;
      key->IML[v] = -1;
      key->ML[v][0] = kc->method; 
      // sonst: ML info wird nicht verwendet...
    }
    ERRORMODELNUMBER = -1;
  } else {
    for (v=0; v<key->ncov; v++) {
      ERRORMODELNUMBER = v;
      if (CovList[key->cov[v].nr].type == ISOHYPERMODEL ||
	  CovList[key->cov[v].nr].type == ANISOHYPERMODEL) {
	error=ERRORHYPERMETHOD; goto ErrorHandling; 
      }
    }
    ERRORMODELNUMBER = -1;
    init_method_list(key); // nicht weiter vorne,
    // da es auf Daten von key zugreift, die weiter vorne noch nicht
    // vorhanden sind 
  }
  if (GENERAL_PRINTLEVEL>2) {
    PRINTF("Lists of methods for the simulation\n===================================\n");
    for (v=0; v<key->ncov; v++) {
      PRINTF("%s: ", CovList[key->cov[v].nr].name);
      for (i=0; i < key->NML[v]; i++) PRINTF("%s, ", METHODNAMES[key->ML[v][i]]);
      PRINTF("\n");
    }
  }
  
  
  while (error != NOERROR) {
    int error_covnr, m;
    error_covnr = NOERROR;

    if (error != ERRORANYLEFT) {
        error_covnr = next_element_of_method_list(key); 
    }

    if (error_covnr) {
      if (GENERAL_PRINTLEVEL>2)
	PRINTF("covariance nr %d (%s) run out of method list.\n", 
	       error_covnr, CovList[key->cov[error_covnr - 1].nr].name);
      if (!user_defined) error = ERROROUTOFMETHODLIST;
      break;
    }
      
    error = err_occurred = NOERROR;// necessary in case ncov=0, since then 
    //                          the following loop is not entered
    // create list of methods used
    for (i=0; i<=(int) Nothing; i++) method_used[i] = 0;
    for (v=0; v<key->ncov; v++)
      method_used[key->cov[v].method] += key->cov[v].left; 
    for (m=CircEmbed; m<=Nothing; m++) {
      if (method_used[m]) {
	insert_method(&last_incompatible, key, (SimulationType) m,
		      do_incompatiblemethod[m]!=NULL); 
	if (method_used[m]==1) {
	  v=0; 
	  while(key->cov[v].method!=m || !key->cov[v].left) v++;
	  key->ML[v][key->IML[v]] = Nothing; 
	  // es lohnt sich nicht, die Methoden durchzuleiern,
	  // die auf die Kovarianzfunktion bereits angewandt wurden und
	  // diese Kovarianzfunktion dabei die einzige war
	}
      }
    }
    
    for (M=0; M<key->n_unimeth; M++) {
      if (!key->meth[M].unimeth_alreadyInUse) {
	// unimeth_alreadyInUse: to avoid that method is called again
	// when next while loop, which would leed to double use of key->meth[M]
	bool left[MAXCOV];
	methodvalue_type *meth;
	meth =  &key->meth[M];
	meth->unimeth_alreadyInUse = true;
	for (v=0; v<key->ncov; v++) left[v] = key->cov[v].left; // sicherung,
	// da init_method ueberschreibt und *danach* entscheidet, ob
	// die methode funktioniert -- ueber key->method ist es zu schwierig
	// da eine methode mehrmals verwendet werden kann
	if (GENERAL_PRINTLEVEL>=4) {
          PRINTF("\n%sInit:  method %d (%s); %d method instance(s) in use\n",
		 errorloc_save, M, METHODNAMES[meth->unimeth], key->n_unimeth);
	  if (GENERAL_PRINTLEVEL>=5)
	    for (v=0; v<key->ncov; v++)
	      PRINTF("%d (%s): left=%s, method=%s\n",
		     v, CovList[key->cov[v].nr].name,
		     key->cov[v].left ? "true" : "false", 
		     METHODNAMES[key->cov[v].method]);
	}
	assert(meth->destruct==NULL && meth->S==NULL);
	if (meth->unimeth > Nothing) {
	  error = ERRORUNKNOWNMETHOD; 
	  goto ErrorHandling;
	}

	err_occurred = init_method[meth->unimeth](key, M);
	if (GENERAL_PRINTLEVEL>=5) {
	    ErrorMessage(Nothing, err_occurred);
	}
	switch (err_occurred) {
	    case NOERROR_REPEAT:
	      insert_method(&last_incompatible, key, meth->unimeth,
			    do_incompatiblemethod[meth->unimeth]!=NULL); // ??!!
	      err_occurred = NOERROR; 
	      break;
	    case NOERROR:
	      break;
	    case NOERROR_ENDOFLIST:
	      delete_method(&M, &last_incompatible, key);
	      err_occurred = NOERROR;
	      break;
	    default:
	      Merr = meth->unimeth; // do not shift after delete_method, since 
	      // the latter changes the value of M !
	      delete_method(&M, &last_incompatible, key); 
	      for (v=0; v<key->ncov; v++) {
		covinfo_type *kc;
		kc = &(key->cov[v]);
		kc->left = left[v];
	      }
	      if (simple_definition && GENERAL_PRINTLEVEL>2) {
		ErrorMessage(Merr,err_occurred);
	      }
	      error = err_occurred; // error may occur anywhere
	      //                     and final err_occured might be 0
	}// switch
      } // if not already in Use
    } // for (M...)

    if (error==NOERROR) {
      for (v=0; v<key->ncov; v++) {
	  if (key->cov[v].left) {
	      error = ERRORANYLEFT;
	      Merr = Nothing;
	      break;
	  }
      }  // if (error == ERRORANYLEFT) continue;
    }
  } // while error!=NOERROR

  if (error!=NOERROR && !simple_definition) {
    int zaehler, vplus;
    if (GENERAL_PRINTLEVEL>1) {
      PRINTF("algorithm partially failed. Retrying. Last error has been:\n");
      ErrorMessage(Merr, error);   
    }
    zaehler = 0;
    for (v=0; v<key->ncov; v++) if (key->cov[v].left) {
      key->cov[v].method = Nothing;
      zaehler++;
    }
    assert(zaehler > 0);
    for (v=0; v<key->ncov; v+=vplus) {
      bool left[MAXCOV];
      covinfo_type *kc, *kcdummy;

      vplus = 1;
      kc = &(key->cov[v]);
      if (kc->left) {
	int w;
        if (key->NML[v]==0) {
	  if (error != ERRORANYLEFT) key->n_unimeth++;
          error = ERROROUTOFMETHODLIST;
	  goto ErrorHandling;
	}
	if (GENERAL_PRINTLEVEL>=5) 
	    PRINTF("'%s' not initiated yet\n", CovList[kc->nr].name);
	for (i=0; i<key->NML[v]; i++) {
	  if (key->ML[v][i] == Nothing) {
	    if (GENERAL_PRINTLEVEL>=6) 
		PRINTF("method '%s' (#%d of %d) jumped\n",
		       METHODNAMES[i], i, key->NML[v]);
	    continue;
	  }
	  if (GENERAL_PRINTLEVEL>=6) 
	      PRINTF("trying '%s' (#%d of %d)\n", 
		     METHODNAMES[key->ML[v][i]], i, key->NML[v]); 

	  for (w=0; w<key->ncov; w++) left[w] = key->cov[w].left; // sicherung,
	  kcdummy = kc;
	  while (v + vplus < key->ncov && kcdummy->op) {
	    if (CovList[kcdummy->nr].type==ISOHYPERMODEL || 
		CovList[kcdummy->nr].type==ANISOHYPERMODEL) {
	      // too complicated for the moment; should be improved
	      if (error != ERRORANYLEFT) key->n_unimeth++;
	      error = ERRORFAILED;
	      goto ErrorHandling;
	    }
	    // multiplication is left
	    kcdummy->method = key->ML[v][i];
	    kcdummy = &(key->cov[v + vplus]);
	    vplus++;
	  }
	  kcdummy->method = key->ML[v][i];
	    
	  M = insert_method(&last_incompatible, key, kcdummy->method,
			    do_incompatiblemethod[kcdummy->method]!=NULL);
	  error = init_method[key->meth[M].unimeth](key, M);
	  if (GENERAL_PRINTLEVEL>=3) {
	    PRINTF("%s has return code %d;", 
		   METHODNAMES[kcdummy->method], error);
	    ErrorMessage(Nothing, error);
	  }
	  if (error != NOERROR && error != NOERROR_REPEAT) {
	    Merr = kc->method;
	    delete_method(&M, &last_incompatible, key);
	    for (w=0; w<key->ncov; w++) {
	      kcdummy = &(key->cov[w]);
	      kc->left = left[w];
	    }
	  } else {
            error = NOERROR;
	    break;
	  }
	} // NML
	if (error != NOERROR) {
	  if (key->ncov>2 && GENERAL_PRINTLEVEL>0) 
	    PRINTF("The given model seems to be complex; try specifying explicitely the simulation methods, see help('GaussRF') and help('PrintMethodList') or even split up in subsequent calls of GaussRF\n");
	  if (GENERAL_PRINTLEVEL>3) {
	    PRINTF("All algorithm for %s failed. The last error has been:\n",
		   CovList[kc->nr].name);
	    ErrorMessage(Merr, error);   
	  }
	  break;
	}
      } // left
    } // for v < key->ncov
  } // error & not simple definition
   
  if (!(key->active = error==NOERROR)) {
    if (error != ERRORANYLEFT) key->n_unimeth++;
    goto ErrorHandling;
  }
  if (GENERAL_PRINTLEVEL>=4) 
      PRINTF("%sSuccessfully initialised.\n", errorloc_save);
  key->compatible = last_incompatible <= 0; // if at most one incompatible
  // method, then no intermediate results must be stored
  strcpy(ERROR_LOC, errorloc_save);
  return NOERROR;
  
  ErrorHandling:
  key->active = false;
  if (GENERAL_PRINTLEVEL>3) {
    ErrorMessage(Merr, error);
    PRINTF("\n");
  }
  return error;
}


void AddTrend(int *keyNr, int *n, double *res, int *error) {
  int i;
  key_type *key;

  *error = NOERROR;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  key = &(KEY[*keyNr]);
  if (!key->active) {*error=ERRORNOTINITIALIZED; goto ErrorHandling;}

//  printf("addtrend %f \n", key->mean);
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



int internal_DoSimulateRF(key_type *key, int nn, double *orig_res) {
  // does not assume that orig_res[...] = 0.0, but it is set
  int error;
  long i, m;
  double  *part_result, *res, realeach=0.0;
  char back[]="\b\b\b\b\b\b\b\b\b\b\b", format[20], prozent[]="%";
  int ni, digits, each=0;
 
  res = orig_res;
  part_result = NULL;
  if (!key->active) {error=ERRORNOTINITIALIZED; goto ErrorHandling;}
  if (nn>1 && GENERAL_PCH[0] != '\0') {
    if (GENERAL_PCH[0] == '!') {
	digits = (nn<900000000) ? 1 + (int) trunc(log((double) nn) / log(10.0)) : 9;
      back[digits] = '\0';
      each = (nn < 100) ? 1 :  nn / 100;
      sprintf(format, "%ss%s%dd", prozent, prozent, digits);
    } else if (GENERAL_PCH[0] == '%') {
      back[4] = '\0';
      realeach = (double) nn / 100.0;
      each = (nn < 100) ? 1 : (int) realeach;
      sprintf(format, "%ss%s%dd%ss", prozent, prozent, 3, prozent);
    } else each = 1;
  } else each = nn + 1;
  GetRNGstate();
  if (!key->compatible) {
      if ((part_result=(double *) malloc(sizeof(double) * key->totalpoints))
	  == NULL) {error=ERRORMEMORYALLOCATION; goto ErrorHandling;}
  }
  for (ni=1; ni<=nn; ni++, res += key->totalpoints) {
    if (key->stop) {error=ERRORNOTINITIALIZED; goto ErrorHandling;}
    if (ni % each == 0) {
      if (GENERAL_PCH[0] == '!')  
	  PRINTF(format, back, ni / each);
      else if (GENERAL_PCH[0] == '%')
	  PRINTF(format, back, (int) (ni / realeach), prozent);
      else PRINTF("%s", GENERAL_PCH);
    }
    for(m=0; m<key->n_unimeth; m++) {
      methodvalue_type *meth;
      meth = &(key->meth[m]);
      if (GENERAL_PRINTLEVEL>=7)
	PRINTF("DoSimu %d %s\n",m, METHODNAMES[meth->unimeth]);
      if (meth->incompatible) {
	if (m==0) {
	  if (GENERAL_PRINTLEVEL>=7) PRINTF("incomp\n");
	  do_incompatiblemethod[meth->unimeth](key, m, res);
	} else {
	  if (GENERAL_PRINTLEVEL>=7) PRINTF("extra %d\n", key->compatible);
	  do_incompatiblemethod[meth->unimeth](key, m, part_result);
	  for (i=0; i<key->totalpoints; i++) res[i] += part_result[i];
	}
      } else {
	do_compatiblemethod[meth->unimeth] (key, m, res); 
      }
    } // for m
  } // for n
  PutRNGstate();
  if (part_result!=NULL) free(part_result);
  if (nn>1 && GENERAL_PCH[0] != '\0') {
    if (GENERAL_PCH[0] == '!' || GENERAL_PCH[0] == '%') PRINTF("%s", back);
    else PRINTF("\n");
  }
  return NOERROR;

 ErrorHandling: 
  PutRNGstate();
  if (part_result!=NULL) free(part_result);
  key->active = false;
  return error;
}


void DoSimulateRF(int *keyNr, int *n, int *pairs, double *res, int *error) {
  // does not assume that orig_res[...] = 0.0, but it is set
  int i, internal_n;
  key_type *key=NULL;

  *error=NOERROR; 
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  key = &(KEY[*keyNr]);
  internal_n = *n / (1 + (*pairs!=0));

  if (!key->active) {
    *error=ERRORNOTINITIALIZED; goto ErrorHandling;
  }

  if (key->n_unimeth == 0 || !key->meth[0].incompatible) {
    long total =internal_n * key->totalpoints ;
   for (i=0; i<total; i++)res[i] = 0.0;
  }
  
  if ((*error = internal_DoSimulateRF(key, internal_n, res)) != NOERROR)
    goto ErrorHandling;

  if (*pairs) {
    double* res_pair;
    long endfor=key->totalpoints * *n / 2;
    res_pair = res + endfor;
    for (i=0; i<endfor; i++) res_pair[i] = -res[i];
  }

  AddTrend(keyNr, n, res, error);
  if (*error) goto ErrorHandling;
 
  if (!(key->active = GENERAL_STORING))
    DeleteKey(keyNr);
  return; 
  
 ErrorHandling: 
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing, *error);
  if (key!=NULL) {
    key->active = false;
    DeleteKey(keyNr);
  }
  return;
}

int InternalSimulate(key_type *key, double *orig_res) {
   if (key->n_unimeth == 0 || !key->meth[0].incompatible) {
     long i, total = key->totalpoints;
     double *sres = orig_res;
     for (i=0; i<total; i++) sres[i] = 0.0;
   }
   return internal_DoSimulateRF(key, 1, orig_res);
}



SEXP GetExtModelInfo(SEXP keynr) {
  int knr;   
  knr = INTEGER(keynr)[0];
  if (knr>=0) {
    if (knr < MAXKEYS) {
      key_type *key;
      key = &(KEY[knr]);
      return GetModelInfo(key->cov, key->ncov, key->totalparam,
			  key->totalpoints);
    }
  } else if (knr == -1) {
    if (user_ncov>0)
      return GetModelInfo(user_keycov, user_ncov, user_anisotropy ?  
			  ANISO + user_dim * user_dim : SCALE, 0);
  } else if (knr == -2) {
    if (Unchecked_ncov > -1)
      return GetModelInfo(Unchecked_keycov, Unchecked_ncov, 
			  Unchecked_anisotropy 
			  ? ANISO + Unchecked_dim * Unchecked_dim : SCALE, 0);  
  }
  return allocVector(VECSXP, 0);
}
