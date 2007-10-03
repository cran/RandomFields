/*
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Simulation of a random field by turning bands;
 see RFspectral.cc for spectral turning bands

 Copyright (C) 2001 -- 2006 Martin Schlather, 

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
//#include <sys/timeb.h>
#include <unistd.h>
#include <assert.h>
#include <R_ext/Applic.h>
#include "RFsimu.h"
#include "RFCovFcts.h"
#include "win_linux_aux.h"

#define MAXNN 100000000.0 /* max number of points simulated on a line */

tbm_param tbm2 =   { 60, 0, 1, 2.0, 0.0, false};
tbm_param tbm3 =   {500, 0, 1, 2.0, 0.0, false};
int TBM_POINTS = 0;
double TBM_CENTER[MAXDIM] = {NA_REAL, NA_REAL, NA_REAL, NA_REAL};

ce_param TBMCE = {false, true, false, true, TRIVIALSTRATEGY, 3, 10000000, 
		  -1e-7, 1e-3, {0, 0, 0, 0}};

SimulationType TBM_METHOD=CircEmbed;

void unitvector3D(int projectiondim, double *deltax, double *deltay, 
		double *deltaz) {
  switch (projectiondim) { // create random unit vector
      case 1 : 
	*deltax= UNIFORM_RANDOM;
	*deltaz = *deltay = 0.0;
	break;
      case 2 :
	*deltaz = 0.0;
	*deltax= UNIFORM_RANDOM;// see Martin's tech rep for details
	*deltay= sqrt(1.0 - *deltax * *deltax) * sin(UNIFORM_RANDOM*TWOPI);
	break;
      case 3 : 
	double dummy;
	*deltaz = 2.0 * UNIFORM_RANDOM - 1.0;
	dummy = sqrt(1.0 - *deltaz * *deltaz);
	*deltay = UNIFORM_RANDOM * TWOPI;
	*deltax = cos(*deltay) * dummy;
	*deltay = sin(*deltay) * dummy;
	break;
      default : assert(false);
  }
}


#define TBMTIME 3
double CovFctTBM3(double *x, int dim, covinfo_arraytype keycov,
		  covlist_type covlist, int ncov, bool anisotropy){
  int v, vold, w, y;
  double dummy, zw, result;
  double r, z[2]; //'' !! 2 !
  double fct[MAXCOV], abl[MAXCOV];
  covinfo_type *kc;
  cov_fct *cov=NULL; /* no meaning -- just avoids warning from -Wall */
  result = 0.0;
  v = 0;

  if (dim==2) { // => s->anisotropy! 
    while (v<ncov) {
      vold = v;
      v++; while ((v<ncov) && keycov[covlist[v-1]].op) v++;
      zw = 1.0;
      for (w=vold; w<v; w++) {
	kc = &(keycov[covlist[w]]);
        cov = &(CovList[kc->nr]);
	z[0] = kc->aniso[0] * x[0];
	z[1] = kc->aniso[TBMTIME] * x[1];
//	printf("%f %f\n", kc->aniso[0], kc->aniso[TBMTIME]); assert(false);
	if (cov->type==ISOTROPIC) {
	  r = sqrt(z[0] * z[0] + z[1] * z[1]);
	  if (r==0) abl[w] = 0;
	  else abl[w] = kc->param[VARIANCE] * kc->aniso[0] * fabs(z[0]) / 
		   r * cov->derivative(&r, kc->param);
	  zw *= (fct[w] = kc->param[VARIANCE] * cov->cov(&r, kc->param));
	} else { 
	  assert(cov->type==SPACEISOTROPIC);
	  abl[w] = kc->param[VARIANCE] * 
	    kc->aniso[0] * cov->derivative(z, kc->param);
	  zw *= (fct[w] = kc->param[VARIANCE] * cov->cov(z, kc->param));
	}
      }
      result += zw;
      if (x[0]!=0.0) {
	// TBM-OP : C + h C' -- C'(0) is often an exception for calculation
	//          but h C' is always zero, so do not consider z[0]=0.0
	dummy = 0.0;
	for (w=vold; w<v; w++) {
	  zw = abl[w];
	  for (y=vold; y<v; y++) if (y!=w) zw *= fct[y];
	  dummy += zw;
	}
	result += fabs(x[0]) * dummy;
      }
    }
  } else {
    assert(dim==1);
    kc = &(keycov[covlist[0]]);
    cov = &(CovList[kc->nr]);
    if (ncov==1 && cov->cov_tbm3!=NULL) {
      z[1] = 0.0; //in case of (CovList[covnr[v]].type==SPACESOTROPIC) 
      z[0] = x[0] * kc->aniso[0];
      return kc->param[VARIANCE] * cov->cov_tbm3(z, kc->param);
    }
    result = CovFct(x, dim, keycov, covlist, ncov, anisotropy);
    if (x[0] != 0.0) {
      result += fabs(x[0]) * DerivCovFct(x, dim, keycov, covlist, ncov);
   }
    // printf("%f %f\n", x[0], result);
  }
//if (x[0]==0.0) printf("%f %f, %f %f: %f \n", x[0], x[1], z[0], z[1], result);
  return result;
}

double CovFctTBM2(double *x, int dim, covinfo_arraytype keycov,
		     covlist_type covlist, int ncov, bool anisotropy) {
// in dieser und der naechste Funktion hat $x$ entweder
// 1 Element (Geradenpunkt der turning bands) oder 2 Elemente
// (Punkt auf Flaeche der turning layers)
  int v, vold;
  double result;
  covinfo_type *kc;
  cov_fct *cov=NULL; /* no meaning -- just avoids warning from -Wall */
  result = 0.0;
  v = 0;
  if (dim==2) { // => s->anisotropy! (time-space simulation)
    int w;
    double z[2], zw; //'' !! 2 !
    while (v<ncov) {
      vold = v;
      v++; while ((v<ncov) && keycov[covlist[v-1]].op) v++; 
      zw = 1.0;
      for (w=vold; w<v; w++) { 
	// there could be multiplication by pure time components (z[0]=0)!
	kc = &(keycov[covlist[w]]);
        cov = &(CovList[kc->nr]);
	z[0] = kc->aniso[0] * x[0];
	z[1] = kc->aniso[TBMTIME] * x[1];
	if (cov->type==ISOTROPIC) {
	  assert(z[0]==0.0);
	  zw *= cov->cov(&(z[1]), kc->param); 
	} else {
	  assert(cov->type==SPACEISOTROPIC);
	  if (z[0]==0.0) zw *= cov->cov(z, kc->param);
	  else zw *= cov->cov_tbm2(z, kc->param);
	}
	zw *= kc->param[VARIANCE];
      }
      result += zw;
    }
  } else {
    assert(dim==1);
    double z;
    for (v=0; v<ncov; v++) {
      //    printf("%d %d %d\n", v, ncov, covlist[v]);
      kc = &(keycov[covlist[v]]);
      cov = &(CovList[kc->nr]);
      // the field must be isotropic, so the x's have been transformed
      // multiplication of covariance functions is not allowed,
      // in contrast to 3dim TBM
      z = x[0] * kc->aniso[0];
      result += (z==0.0) ? kc->param[VARIANCE] 
	  : kc->param[VARIANCE] * cov->cov_tbm2(&z, kc->param);
    }
  }
  return result;
}

typedef struct TBM2_integr {
    int dim, ncov;
    covinfo_type *keycov;
    int *covlist;
    bool anisotropy;
    double x, time;
} TBM2_integr;


void TBM2NumIntegrFct(double *u, int n, void *ex) {
  int i;
  TBM2_integr *info;
  double z[2];
  info = (TBM2_integr*) ex;
  z[1] = info->time;
  for (i=0; i<n; i++) {
       z[0] = info->x * sqrt(1 - u[i] * u[i]);
       u[i] = CovFctTBM3(z, info->dim, info->keycov, info->covlist, 
			   info->ncov, info->anisotropy);
  }
}

double CovFctTBM2Num(double *x, int dim, covinfo_arraytype keycov,
		     covlist_type covlist, int ncov, bool anisotropy) {
// in dieser und der naechste Funktion hat $x$ entweder
// 1 Element (Geradenpunkt der turning bands) oder 2 Elemente
// (Punkt auf Flaeche der turning layers)
  TBM2_integr info;
#define MAXSUBDIVISIONS 100
  static double a = 0, 
      b = 1, 
      eps = 1e-10; // 1e-15
  static int maxsubdivisions = MAXSUBDIVISIONS, 
      lenw = 4 * MAXSUBDIVISIONS;
  double result, abserr, work[4 * MAXSUBDIVISIONS];
  int subdivisions, integralevaluations, Xerror, iwork[MAXSUBDIVISIONS];

//  printf("dim=%d \n", dim);
  info.dim = dim;
  info.keycov = keycov;
  info.covlist = covlist;
  info.ncov = ncov;
  info.anisotropy = anisotropy;
  info.x = x[0];
  info.time = (dim==2) ? x[1] : 0.0;

  Rdqags(TBM2NumIntegrFct, (void *) &info, &a, &b, &eps, &eps,
	 &result, &abserr, &integralevaluations, &Xerror,
	 &maxsubdivisions, &lenw, &subdivisions, iwork, work);
//   printf("tbm2num %f %f %d %d %e %f\n", x[0], x[1], integralevaluations, subdivisions, abserr, result);

  return result;
}


void TBM_destruct(void **S) 
{
  if (*S!=NULL) {
   TBM_storage *s;
    s = *((TBM_storage**)S);
    if (s->simuline != NULL) free(s->simuline);
    if (s->x != NULL) free(s->x);
    DeleteKeyNotTrend(&(s->key));
    free(*S);
    *S = NULL;
  }
}

void SetParamTBMCE(int *action, int *force, double *tolRe, double *tolIm,
		   int *trials, double *mmin, int *useprimes, int *strategy,
		   double *maxmem, int *dependent)
{
  int TRUE=1;
  SetParamCE(action, force, tolRe, tolIm, trials, &TRUE, 
	     mmin, useprimes, strategy,
	     maxmem, dependent, &TBMCE,"TBMCE");
}


void SetParamLines(int *action,int *nLines, double *linesimufactor, 
		   double *linesimustep, int *every, 
		   int *layers, int*tbm2num,
		   tbm_param *tbm, 
		   char *name, double *linesimustep_o, double *linesimufactor_o) 
{
  if (*action) {
    tbm->lines=*nLines;
    if (*linesimufactor < 0.0 || *linesimustep < 0.0) 
	PRINTF("Both %s.linesimustep and %s.linesimufactor and must be non-negative\n", name, name); 
    tbm->linesimufactor=*linesimufactor < 0.0 ? 0.0 : *linesimufactor;
    tbm->linesimustep=*linesimustep < 0.0 ? 0.0 : *linesimustep;    
    if (*linesimufactor!=0.0 && *linesimustep!=0.0 && GENERAL_PRINTLEVEL>0 &&
	(*linesimustep_o != *linesimustep || 
	 *linesimufactor_o != *linesimufactor)) {
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("%s.linesimufactor is ignored if %s.linesimustep is positive!\n", 
	       name, name); 
      *linesimustep_o = *linesimustep;
      *linesimufactor_o = *linesimufactor;
      // else, i.e. both are zero, the old values are kept!
    }
    tbm->every=*every;
    tbm->layers = *layers;
    tbm->tbm2num=(bool) *tbm2num;
  } else {
    *nLines=tbm->lines;
    *linesimufactor=tbm->linesimufactor;
    *linesimustep=tbm->linesimustep;  
    *every=tbm->every;
    *layers = tbm->layers;
    *tbm2num = (int) tbm->tbm2num;
  }
}

void SetParamTBM2(int *action, int *nLines, double *linesimufactor, 
		  double *linesimustep, int *every, int *tbm2num,
		  int *layers) {
  static double linesimustep_old=-1, linesimufactor_old=-1;
  SetParamLines(action, nLines, linesimufactor, linesimustep, every,
		layers, tbm2num,
		 &tbm2, "TBM2", &linesimustep_old, &linesimufactor_old);
}

void SetParamTBM3(int *action,int *nLines, double *linesimufactor, 
		  double *linesimustep, int *every,
		  int *layers) {
  static double linesimustep_old=-1, linesimufactor_old=-1;
  int tbm2num=0;
  SetParamLines(action, nLines, linesimufactor, linesimustep, every,
		layers, &tbm2num,
		&tbm3, "TBM3", &linesimustep_old, &linesimufactor_old);
}

void SetParamTBM(int *action, int *tbm_method, double *center, int *points) {
  int d;
  if (*action) {
   for (d=0; d<MAXDIM; d++) TBM_CENTER[d] = center[d];
   TBM_POINTS = *points;
   // method must be the last one because of the return in the middle
   if ((*tbm_method<=(int) Nothing) && (*tbm_method>=0)) {      
      SimulationType m;
      m=(SimulationType) *tbm_method;
      if ((m==CircEmbed) || (m==Direct) || (m==Nothing)){
	TBM_METHOD = m;
	return;
      }
    }
    if (GENERAL_PRINTLEVEL>0) {
      PRINTF("Specified method [%d] is not allowed. Method set to `Nothing'\n",
	     *tbm_method);
    }
    TBM_METHOD = Nothing;
  } else {
    *tbm_method = (int) TBM_METHOD;
    for (d=0; d<MAXDIM; d++) center[d] = TBM_CENTER[d];
    *points = TBM_POINTS;
  }
}



// aufpassen auf ( * * 0 )
//               ( * * 0 )
//               ( 0 0 1 )

int init_turningbands(key_type *key, SimulationType method, int m)
{
  methodvalue_type *meth; 
  covinfo_type *kc, *first;
  int Xerror=NOERROR, d, tbm_dim, v, nonzero_pos, lasttimecomponent=-1, 
      firsttimecomponent=-1, lastmatching=-1, firstmatching=-1, 
      covnr[MAXCOV], cum_nParam[MAXCOV+1], op[MAXCOV], actcov,
      totaltimespacedim, w, Aniso, loop, iloop;
  double ParamList[TOTAL_PARAM * MAXCOV], f_epsilon = 1e-15, 
      linesimuscale, quot, diameter;
  tbm_param *tbm;
  SimulationType tbm_method;
  bool ce_dim2=false, // just to avoid error messages 
      newadditive, equal, spacecomponent_found, tbm2num; 
  char errorloc_save[nErrorLoc], msg[15];
  TBM_storage *s;
  param_type simuparam;

  assert(key->timespacedim <= 4);
  meth = &(key->meth[m]);
   /* multiplication of the anisotropy matrix from the right */
  nonzero_pos = -1;
  SET_DESTRUCT(TBM_destruct, m);

#define BUFFER 3.0

//  assert(key->cov[0].aniso[0] != 0.0);
//  printf("tbm ansio %f\n", key->cov[0].aniso[0] );

  assert(key->covFct == CovFct);
  if ((meth->S = malloc(sizeof(TBM_storage)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (TBM_storage*) meth->S;
  s->simuline=NULL;
  s->x = NULL;
  KEY_NULL(&(s->key));
  first = NULL;
  s->timespacedim = key->timespacedim;
  tbm2num = false; 

  // get list of covariance functions that are of interest
  // CovFctTBM2 might be changed to CovFctTBM2Num
  if (method==TBM2) { tbm_dim = 2; tbm = &tbm2; }
  else {              tbm_dim = 3; tbm = &tbm3; }

 
  /****************************************************************/
  /*            Extraction of matching covariances                */
  /****************************************************************/
  // einfache checks
  // extraktion von covarianzfunktion mit gleichartiger anisotropie
  //    time-komponente wird nicht ge-checked da letzter parameter true
  // bereitstellung von param, quotient, multiply
  
  /* in tbm, time component must be ALWAYS alone, i.e. matrix looks like
   * * 0 
   * * 0
   * * x
   -- this is checked late on
   
   here x may be any constant. In case x is zero 
   (for all matrices (0..actcov-1), the last dimension
   vanishes (performed in the next step)
   
   the asterixes build submatrices for (v=1..actcov-1) that are multiple of the 
   submatrix for v=0, already checked by FIRST_CHECK_
  */
  
  if (GENERAL_PRINTLEVEL>4) 
      PRINTF("determining first nonzero_pos and TBM simulation dimension...\n");
  s->ce_dim = 0;
  for (v=0; v<key->ncov; v++) {
    ERRORMODELNUMBER = v;
    kc = &(key->cov[v]);
    if ((kc->method==method) && (kc->left)) {
      first = kc; 
//  printf("tbm ansio first %f\n", first->aniso[0] );
      if (CovList[first->nr].type==ISOHYPERMODEL || 
	  CovList[first->nr].type==ANISOHYPERMODEL) {
	  Xerror=ERRORHYPERNOTALLOWED; 
	  goto ErrorHandling;
      }
      if (key->anisotropy) {
	  firstmatching = ANISO;
	  lastmatching = key->totalparam - 1;     
      } else {
	 firstmatching = lastmatching = SCALE; 
      }
      for (nonzero_pos=firstmatching; 
	   nonzero_pos<=lastmatching && first->param[nonzero_pos]==0.0; 
	   nonzero_pos++);
      if (nonzero_pos > lastmatching) {
	  Xerror=ERRORLOWRANKTBM; goto ErrorHandling;
      }
      s->simuspatialdim = first->reduceddim;
      if (key->Time) {
	ce_dim2 = tbm->layers > 0;
	lasttimecomponent = lastmatching;
        firsttimecomponent = key->totalparam - key->timespacedim;
 	ce_dim2 |= CovList[first->nr].type==SPACEISOTROPIC ||
	    first->reduceddim == 1 + tbm_dim;
//	printf("%d\n", ce_dim2); assert(false);
	while (!ce_dim2 && (++v)<key->ncov && kc->op) {
	  // here: only check whether within the multiplicative model
	  // the anisotropy matrices are not multiplicatives of each other.
	  // If so, ce_dim2 is necessarily true.
	  // At this point it is not checked whether the matrix combinations
	  // are allowed ones.
	  double quot;
	  int w;
	  kc = &(key->cov[v]);
          ce_dim2 |= CovList[kc->nr].type==SPACEISOTROPIC;
	  quot = kc->param[nonzero_pos] / kc->param[nonzero_pos];
          for (w=ANISO; w<=lasttimecomponent; w++) {
	    if (fabs(kc->param[w] * quot - kc->param[w]) >= 
		(fabs(kc->param[w]) + (double)(kc->param[w]==0.0)) *
		f_epsilon) {
	      ce_dim2 = true;
	      break;
	    }
	  }
	}
	if (ce_dim2) {
	  aniso_type aniso;
	  matrixrotat(&(first->param[ANISO]), first->dim - 1, first->dim, 
		      &(s->simuspatialdim), aniso);
	  if (tbm->layers <= 0) {
	    Xerror=ERRORLAYERS; goto ErrorHandling;
	  }
	}
      } else {
	ce_dim2 = false;
	if (tbm->layers > 1) {
	    Xerror=ERRORLAYERS; goto ErrorHandling;
	}
      }

      if (s->simuspatialdim > tbm_dim) {
	  Xerror=ERRORWRONGDIM; goto ErrorHandling;
      }
      // printf("%d %d %d \n",s->simuspatialdim, tbm_dim. 1);  assert(false);

      s->ce_dim = 1 + (int) ce_dim2;
      if (ce_dim2) {
	lastmatching = firsttimecomponent - 1;
	if (nonzero_pos >= firsttimecomponent) {
	   Xerror=ERRORLOWRANKTBM; goto ErrorHandling;}
      }
      break;
    }
  }
  ERRORMODELNUMBER = -1;

  if (first == NULL) {
    Xerror=ERRORFAILED;
    goto ErrorHandling;
  }
  assert(nonzero_pos > 0);
  if (s->ce_dim == 0) { /* no covariance for the considered method found */
    Xerror=NOERROR_ENDOFLIST;
    goto ErrorHandling;
  }

  // printf("tbm ansio first %f %d\n", first->aniso[0], ce_dim2 );

  /****************************************************************/
  /*          determine length and resolution of the band         */
  /****************************************************************/
 // sort out which finer grid should be used in the spatial domain for
  // low dimensional simulation
  if (GENERAL_PRINTLEVEL>4) PRINTF("\ndetermination of the band length...");

//  printf("tbm ansio first %f\n", first->aniso[0] );

  // minimum distance between the points   
  loop = 1;
  if (tbm->linesimustep  > 0.0) linesimuscale = 1.0 / tbm->linesimustep;
  else if (tbm->linesimufactor > 0) {
    double mindelta;
    mindelta = RF_INF; 
    if (key->grid) {
      // if linesimufactor, the fine grid scale is proportional to smallest 
      // distance in the grid. 
      for (d=0; d<s->simuspatialdim; d++) {
        if ((key->x[d][XSTEP]<mindelta) && (key->x[d][XSTEP]>0)) 
          {mindelta=key->x[d][XSTEP];}
      }
    } else {
      int j, i, d;	
      if (key->totalpoints > 50000) {
	  Xerror=ERRORTOOMANYPOINTS; goto ErrorHandling;
	  /* algorithmus kann verbessert werden, so dass diese Fehlermeldung 
	     nicht mehr notwendig ist! */ 
      }
      for (i=0; i<key->totalpoints; i++) {
	 for (j=0; j<i; j++) {
	    register double diff, dist;
	    for (dist=0.0, d=0; d<key->timespacedim; d++) {
		// if (i > 2525) printf("d=%d %d %d\n", d, i, j);
	      diff = key->x[d][i] - key->x[d][j]; 
	      dist += diff * diff;
	    }
	    if (dist>0 && dist<mindelta) mindelta=dist; 
	  }
	} // if linesimustep==0.0
      mindelta = sqrt(mindelta);
    }
    linesimuscale = tbm->linesimufactor / mindelta;
  } else {
    linesimuscale=1.0;
    loop = 2;
    if (TBM_POINTS < 4) { Xerror=ERRORTBMPOINTSZERO; goto ErrorHandling; }
  }
  //  printf("tbm ansio first %f\n", first->aniso[0] );

  // multiplication der x-skalenparameter so dass letztendlich
  // auf der TBM-Geraden mit Abstand 1 der Punkte simuliert werden kann
  // mit dem Vorteil des einfachen Zugriffs auf die simulierten Werte 
  s->simugrid = first->simugrid;
  diameter = -1.0;
  for (iloop=0; iloop<loop; iloop++) {
    int i;
    double dummylx[MAXDIM];
    // the loop must be performed twice only if 
    // linesimuscale is defined by means of TBM_POINTS and diameter.
    // the latter is dertermined in the first loop
    if (s->x != NULL) {// i.e, if iloop==1
      free(s->x);
      s->x = NULL;
    }
    if (key->anisotropy)
      for (i=firstmatching; i<=lastmatching; i++) 
	simuparam[i] = first->param[i] * linesimuscale;
    else {
      simuparam[SCALE] = first->param[SCALE] / linesimuscale;
    }

    bool Timedummy;
    Timedummy = false;
    if (ce_dim2) {
      for (i=lastmatching + 1; i<lasttimecomponent; simuparam[i++]=0.0);
      simuparam[lasttimecomponent] = 1.0;
      GetTrueDim(key->anisotropy, key->timespacedim, simuparam,
		 SPACEISOTROPIC, &Timedummy, &totaltimespacedim, s->aniso);
      assert(s->simuspatialdim == totaltimespacedim - 1);
      for (i=lastmatching + 1; i<lasttimecomponent; i++)
	  simuparam[i] = first->param[i] * linesimuscale; // only used in Getxsimugr
    } else {
      GetTrueDim(key->anisotropy, key->timespacedim, simuparam,
		 ISOTROPIC, &Timedummy, &totaltimespacedim, s->aniso);
      assert(totaltimespacedim = first->reduceddim);
    }
    s->reduceddim = totaltimespacedim;

    // for (i=ANISO; i<ANISO + 4; i++) printf("%f ", simuparam[i]);  printf("-----\n");
//    printf("%d %d %d\n",first->reduceddim,s->simuspatialdim, totaltimespacedim);
    
    // ******************************
    // diameter of the simulation area 
    assert(s->x == NULL);
    if ((Xerror=Transform2NoGrid(key, s->aniso,  totaltimespacedim, 
				s->simugrid, &(s->x))) != NOERROR)
      goto ErrorHandling;
 
/*
  printf("TRANSF %d %d %d %d  %d %d\n", iloop, key->anisotropy, key->Time,
	   key->grid, totaltimespacedim, s->simugrid);
  printf("s->x %f %f %f\n\n", s->x[0], s->x[1] , s->x[2]);
  printf("aniso %f %f %f %f %f %f %f %f %f \n\n", 
	 s->aniso[0], s->aniso[1], s->aniso[2], s->aniso[3],
	 s->aniso[4], s->aniso[5], s->aniso[6], s->aniso[7],
	 s->aniso[8]);
  printf("first %f %f %f %f \n\n", 
	 first->aniso[0], first->aniso[1], first->aniso[2], first->aniso[3]);
  for (i=0; i<12*3; i++) printf("%f ", s->x[i]);
  printf("\n");
  // assert(false);
  */
    
    GetCenterAndDiameter(key, s->simugrid, s->simuspatialdim,  
			 totaltimespacedim, s->x, s->aniso,
			 s->center, dummylx, &diameter);

    //  printf("center %f %f %f\n", s->center[0], s->center[1], s->center[2]);

    if (!ISNA(TBM_CENTER[0])) {
      double center[MAXDIM], diff;
      int n,k,g;
      for (g=1; g<key->timespacedim; g++) 
	if (!R_FINITE(TBM_CENTER[g])) {
	  Xerror=ERRORTBMCENTER; goto ErrorHandling; 
	}
      for (n=k=0; k<s->simuspatialdim; k++) {
	register double dummy;
	dummy = 0.0;
	for (g=0; g<key->timespacedim; g++) {
	  dummy += s->aniso[n++] * TBM_CENTER[g];
//	printf("center %d %d %d %f %f %f\n", k, g, n,dummy, s->aniso[n-1], TBM_CENTER[g]);
	}
	center[k] = dummy;
      }
      diff = 0.0;
      for (k=0; k<s->simuspatialdim; k++) {
	diff += (center[k] - s->center[k]) * (center[k] - s->center[k]);
	s->center[k] = center[k];
      }
      diameter += 2.0 * sqrt(diff);
    } // else printf("ISNA CENTER ** \n");
// printf("tbm: %d %d %d %f %f %f %d %d\n", 
//        iloop, loop, TBM_POINTS, diameter, linesimuscale, s->aniso[0],
//        key->grid, s->simugrid);
    if (loop==2 && iloop==0) {
      linesimuscale = (TBM_POINTS - BUFFER) / diameter;
    } else diameter = trunc(BUFFER + diameter);
    
    //   printf("diameter %f \n", diameter); 
    //   assert(false);
    
  } // loop


  if (diameter<=0.0) {
    assert(diameter==0.0);
    // Xerror=ERRORTRIVIAL; goto ErrorHandling;
  }
  if (TBM_POINTS > 0) {
    if (diameter > TBM_POINTS) { 
      PRINTF("tbm points minimum=%f, but TBM_POINTS is %d.\n", 
	     diameter, TBM_POINTS);
      Xerror=ERRORTBMPOINTS; goto ErrorHandling;
    } else diameter = (double) TBM_POINTS;
  }
  if (diameter > MAXNN) {Xerror=ERRORNN; goto ErrorHandling;}


  if (s->simugrid) {
    Getxsimugr(key->x, simuparam, key->timespacedim,  key->anisotropy, s->xsimugr); // see do_tbm
    for (d=0; d<key->timespacedim; d++) {
      s->genuine_dim[d] = first->genuine_dim[d];
    }
  }

//////////////////////////////////////////////////////////////////////
  if (GENERAL_PRINTLEVEL>4) PRINTF("extracting matching covariances...");
  actcov=0;
  cum_nParam[actcov]=0;
  spacecomponent_found = newadditive = true;
  for (v=0; v<key->ncov; v++) {
    ERRORMODELNUMBER = v;
    kc = &(key->cov[v]);
    if ((kc->method==method) && (kc->left)) {
      cov_fct *cov;
      cov = &(CovList[kc->nr]);
      assert((kc->nr>=0) && (kc->nr<currentNrCov));
      assert(kc->param[VARIANCE]>=0.0);
      meth->covlist[actcov] = v;
      if (v<key->ncov-1 && kc->op) {
	if (key->cov[v+1].method!=method) {
          if (GENERAL_PRINTLEVEL>0) 
            PRINTF("severe error -- contact author, please");
	  Xerror=ERRORMETHODMIX; goto ErrorHandling;
	}
	if (method==TBM2) {
	  Xerror=ERRORNOMULTIPLICATION; goto ErrorHandling;
	}
      }
      if (cov->type==ISOHYPERMODEL || cov->type==ANISOHYPERMODEL) {
	  Xerror=ERRORHYPERNOTALLOWED; goto ErrorHandling;}
      if (cov->implemented[method] != IMPLEMENTED &&
	  cov->implemented[method] != NUM_APPROX) {
	  Xerror = ERRORNOTDEFINED; goto ErrorHandling;}
      if ((!key->anisotropy && cov->type!=ISOTROPIC) ||
	  cov->type==ANISOTROPIC) { 
	  Xerror = ERRORISOTROPICMETHOD; goto ErrorHandling;}
      if ((cov->type==SPACEISOTROPIC) && (s->ce_dim==1)) {
	  // current initialisation does not match this further (additive)
	  // component, so left to another initialisation of tbm
	  kc->left = true;
	  while (actcov>0 && key->cov[meth->covlist[actcov-1]].op) {
	      key->cov[meth->covlist[actcov-1]].left = true;
	      actcov--;
	  }
	  assert(actcov>=0);
	  while (v<key->ncov && key->cov[v].op) v++;  
	  continue;
      }
      sprintf(msg, " (%s)", METHODNAMES[method]);
      if ((Xerror=check_within_range(kc->param, cov, tbm_dim + (int) ce_dim2, 
				     msg)) != NOERROR) {
	goto ErrorHandling;
      }    
      Xerror=cov->check(kc->param, tbm_dim + (int) ce_dim2, 
			tbm_dim + (int) ce_dim2, method);
      if (method==TBM2) {
	if (Xerror && (Xerror!=ERRORCOVNUMERICAL || !tbm->tbm2num)) {
	  goto ErrorHandling;
	}
	if (!tbm->tbm2num && cov->implemented[method]==NUM_APPROX) { 
//	    printf("%d %d %d\n", Xerror, ERRORCOVNUMERICAL, tbm->tbm2num)
//	  assert(false);
	  Xerror = ERRORNOTDEFINED; goto ErrorHandling;
	}
	tbm2num |= Xerror || (cov->implemented[method]==NUM_APPROX) ||
	    (ce_dim2 && cov->type==ISOTROPIC);
      } else { // method = TBM3
	if (Xerror) goto ErrorHandling;
      }

      /* check whether the multiplicative construction are consistent,*/
      /* i.e. that the SPATIAL part of the anisotropy matrices are    */
      /* multiplicatives of each other (more precisely, the remaining */
      /* ones are multiplicatives of the first.)                      */
      /* The time part (if ce_dim2) of the matrix must be of the form */
      /* (0,..,0,*).  */
      /* A further goal of this part is the collection of additive    */
      /* blocks that have the same anisotropy matrix structure, in    */
      /* order to minimise the simulation time                        */
      
      // nonzero_pos gives the position of the first element in the 
      // first matrix of a multiplicative model that is not zero
      quot = kc->param[nonzero_pos] / first->param[nonzero_pos];
      assert(firstmatching <= lastmatching);
      for (w=firstmatching; w<=lastmatching; w++) {
        equal = fabs(first->param[w] * quot - kc->param[w])
	    < (fabs(kc->param[w]) + (double)(kc->param[w]==0.0)) *
	    f_epsilon;
	if (!equal) break; /* even not equal up to a multiplicative constant */
      } 
      if (!equal) {
	kc->left = true;
	while (actcov>0 && key->cov[meth->covlist[actcov-1]].op) {
	  key->cov[meth->covlist[actcov-1]].left = true;
	  actcov--;
	}
	if (actcov==0) { Xerror=ERRORANISOMIX; goto ErrorHandling; }
	while (v<key->ncov && key->cov[v].op) v++;  
	continue;
      }
    
      int kappas;
      kappas = cov->kappas(1);
      assert(kappas == cov->kappas(3)); // kappas should be independent of dim!
      kappas++;
      cum_nParam[actcov + 1] = 
	cum_nParam[actcov] + kappas + 1 + 3 * (int) ce_dim2;
      memcpy(&(ParamList[cum_nParam[actcov]]), kc->param, 
	     sizeof(double) * kappas);
      Aniso = cum_nParam[actcov] + kappas;

      ParamList[Aniso] = 
	  ((key->anisotropy) ? quot : 1.0 / quot) / linesimuscale;

//      printf("quot %f\n", quot);

      covnr[actcov] = kc->nr;
      op[actcov] = kc->op;

      if (ce_dim2) { 
        /* is it of the form   * * 0
	   .                   * * 0       ?
	   .                   * * x 
	*/	
	for (w=firsttimecomponent; w<lasttimecomponent; w++) {
          // hier < , sonst <= !!
	  if (kc->param[w]!=0.0) {
	      Xerror = ERRORTIMESEPARATE;
	      goto ErrorHandling;
	  }
	}
	ParamList[Aniso + 1] = ParamList[Aniso + 2] = 0.0;
	ParamList[Aniso + 3] = kc->param[lasttimecomponent];

	if (method==TBM2 && !tbm2num) {
	  // TBM2: per multiplication block at most 1 covariance function 
	  //       with non vanishing spatial components, otherwise 
	  //       numerical integration
	  if (newadditive) spacecomponent_found = false;
	  if (spacecomponent_found) {
	    for (w=firstmatching; w<=lastmatching; w++)
	      if (kc->param[w] != 0.0) {
		  tbm2num = true;
		  break;
	      }
	  } else {
	    for (w=firstmatching; w<=lastmatching; w++)
	      if (spacecomponent_found = (kc->param[w]!=0.0)) break;
	  }
	} // method==TBM2
      } // cedim
      kc->left=false;
      actcov++;
    } // kc->left
    newadditive = kc->op == 0; // geaendert 22.1.06
  }  // v
  if (tbm2num && GENERAL_PRINTLEVEL > 1)
    PRINTF("\tnumerical evaluation of the TBM operator\n");

  ERRORMODELNUMBER = -1;
  meth->actcov = actcov;

  // xline[XEND]: number of points on the TBM line
  double xline[3];
  xline[XSTART] = 1.0;
  xline[XSTEP] = 1.0;
  xline[XEND] = (long) diameter; // diameter > MAXNN must be first since
  //                                risk of overflow !
  // CovList.tbm_method is not used yet
  // if (TBM_METHOD==Nothing) {tbmmethod[i]=CovList[*covnr].tbm_method;}
   tbm_method = TBM_METHOD;
   strcpy(errorloc_save, ERROR_LOC);
   sprintf(ERROR_LOC, "%s%s: ", errorloc_save, METHODNAMES[method]);
   ce_param ce_save;
   memcpy(&ce_save, &CIRCEMBED, sizeof(ce_param));
   memcpy(&CIRCEMBED, &TBMCE, sizeof(ce_param));
   while (tbm_method != Forbidden) {
     int tbmmethod[MAXCOV], i;
     for (i=0; i<meth->actcov; i++) tbmmethod[i]=tbm_method;

     //   printf("tbm %d %d\n", covnr[0], meth->actcov);

     Xerror = 
	 internal_InitSimulateRF(xline, key->T, 1, 3, true, 
				 ce_dim2 /* Time */, covnr, ParamList, 
				 cum_nParam[actcov], meth->actcov,
				 true /* anisotropy */, op, tbmmethod,
				 DISTR_GAUSS, &(s->key),
				 0 /* natural scaling */,
				 (method == TBM3) ?  CovFctTBM3 :
				 tbm2num ? CovFctTBM2Num : CovFctTBM2
	     );
     if (Xerror==NOERROR) break;
      switch(tbm_method) {
	  case Direct : tbm_method=CircEmbed; break;
	  case CircEmbed : tbm_method=Forbidden; break;
          default: Xerror=ERRORNOTDEFINED;
      }
   }
   memcpy(&CIRCEMBED, &ce_save, sizeof(ce_param));
   strcpy(ERROR_LOC, errorloc_save);
   if (Xerror != NOERROR) goto ErrorHandling; // do not put this line before the
   // two preceeding ones!
  
  if ((s->simuline=(double *) calloc(s->key.totalpoints, sizeof(double)))==0){
    Xerror=ERRORMEMORYALLOCATION; goto ErrorHandling;}
  if (GENERAL_PRINTLEVEL>4) PRINTF("\ntbm is now initialized.\n");
 
  // im folgenden unterscheidung bzgl. anisotropy, da 
  // oben alle Kovarianzfunktionen zusammengefasst werden mit derselben
  // (An)isotropie-Struktur, hei3t nur bei echter Anisotropie wird TBM2
  // etc. mehrfach aufgerufen!
  if (key->anisotropy) return NOERROR_REPEAT;
  else return NOERROR;
  
 ErrorHandling:
  return Xerror;
}

int init_turningbands2(key_type *key, int m) {
 return(init_turningbands(key, TBM2, m));
}

int init_turningbands3(key_type *key, int m) {
  return(init_turningbands(key, TBM3, m));
}

void do_turningbands(key_type *key, int m, double *res)
{ 
  double phi, *simuline, centerx, centery, centerz, nnhalf, 
      incx, incy,  incz, inct=0.0, stepx, stepy, stepz, stept;
  long n, nn, totpoints, ntot;
  int nt, simutimespacedim, idx, gridlent;

  tbm_param *tbm;
  methodvalue_type *meth; 
  TBM_storage *s;

  meth = &(key->meth[m]);
  s = (TBM_storage*) meth->S;
  nn = s->key.length[0];
  ntot = s->key.totalpoints;
  simutimespacedim = s->reduceddim;
  nnhalf = 0.5 * (double) nn; 
  simuline = s->simuline; 
  tbm = (meth->unimeth==TBM2) ? &tbm2 : &tbm3;

  for (n=0; n<key->totalpoints; n++) res[n]=0.0; 
  if (s->simugrid) { // old form, isotropic field
    int nx, ny, nz, gridlenx, gridleny, gridlenz, ix, iy, iz, it, keyidx;
    long zaehler;
    double xoffset, yoffset, zoffset,  toffset, e[3 + 1], deltaphi;
#define ezero 3 // must match e[3 + 1] above
    e[0] = e[1] = e[2] = e[ezero] = 0.0;
    stepx = stepy = stepz = stept = centerx = centery = centerz = 0.0;
    gridlenx = gridleny = gridlenz = gridlent = 1;
    ix = iy = iz = it = ezero;
    idx = simutimespacedim;
    keyidx = key->timespacedim;
    if (s->ce_dim==2) {
      gridlent = key->length[--keyidx]; // could be one !!
      stept = s->xsimugr[XSTEPD[keyidx]];	    
      inct = (double) nn;
      idx--; // 20.1.06 eingefuegt
    }

//  printf("%d %f %f %f \n", 0, s->xsimugr[XSTARTD[0]], s->xsimugr[XSTEPD[0]], s->xsimugr[XENDD[0]]);
//    printf("%d %f %f %f \n", 1, s->xsimugr[XSTARTD[1]], s->xsimugr[XSTEPD[1]], s->xsimugr[XENDD[1]]);
//    printf("%d %f %f %f \n", 0, s->xsimugr[XSTARTD[2]], s->xsimugr[XSTEPD[2]], s->xsimugr[XENDD[2]]);
//    printf("%d %d %d\n", idx, s->timespacedim, keyidx); //assert(false);
    switch (keyidx) {
	case 4 : 
	  gridlent = key->length[--keyidx];
	  if (s->genuine_dim[keyidx]) {
	    stept = s->xsimugr[XSTEPD[keyidx]];	    
	    it = --idx;	
	  } // no break;
	case 3 : 
//	  printf("keyidxz %d \n", keyidx);
	  gridlenz = key->length[--keyidx];
	  if (s->genuine_dim[keyidx]) {
	    stepz = s->xsimugr[XSTEPD[keyidx]];	  
	    iz = --idx;
	  }  // no break;
	case 2 : 
//	  printf("keyidxy %d idx %d \n", keyidx, idx);
	  gridleny = key->length[--keyidx];
	  if (s->genuine_dim[keyidx]) {
	    stepy = s->xsimugr[XSTEPD[keyidx]];	 
//	  printf("step y %f \n", stepy);
	    iy = --idx;
	  }  // no break;
	case 1 : 
//	  printf("keyidxx %d \n", keyidx);
	  gridlenx = key->length[--keyidx];
	  if (s->genuine_dim[keyidx]) {
	    stepx = s->xsimugr[XSTEPD[keyidx]];	  
	    ix = --idx;
	  } 
	  break;
	default : assert(false);
    }
    
    switch (s->simuspatialdim) {
	case 3 : 
	    centerz = s->center[2] - s->x[XSTARTD[2]]; // no break;
	case 2 : 
	    centery = s->center[1] - s->x[XSTARTD[1]]; // no break;
	case 1 : 
	    centerx = s->center[0] - s->x[XSTARTD[0]];
	    break;
	default : assert(false);
    }
    assert(meth->unimeth==TBM2 || meth->unimeth==TBM3);
     
    deltaphi = PI / (double) tbm->lines;
    phi = deltaphi * UNIFORM_RANDOM; 
    for (n=0; n<tbm->lines; n++) {
      if (tbm->every>0  && (n % tbm->every == 0)) PRINTF("%d \n",n);
      if (meth->unimeth==TBM2) {
	phi += deltaphi;
	e[0] = sin(phi); 
	e[1] = cos(phi);
	toffset= nnhalf - centery * e[1] - centerx * e[0];
      } else {
  	unitvector3D(s->simuspatialdim, &(e[0]), &(e[1]), &(e[2]));
 	toffset= nnhalf - centerz * e[2] - centery * e[1] - centerx * e[0];
      }
      incx = e[ix] * stepx;
      incy = e[iy] * stepy;
      incz = e[iz] * stepz;
      if (s->ce_dim == 1) inct = e[it] * stept; // else (double) nn !!

      if (s->key.n_unimeth == 0 || !s->key.meth[0].incompatible) {
	int i;
	long total = s->key.totalpoints;
	double *sres = simuline;
	for (i=0; i<total; i++) sres[i] = 0.0;
      }

      internal_DoSimulateRF(&(s->key), 1, simuline);
      zaehler = 0;
      for (nt=0; nt<gridlent; nt++) {
	zoffset = toffset;
	for (nz=0; nz<gridlenz; nz++) {	  
	  yoffset = zoffset;
	  for (ny=0; ny<gridleny; ny++) {	  
	    xoffset = yoffset;
	    for (nx=0; nx<gridlenx; nx++) {
              register long longxoffset = (long) xoffset;
	      
//     printf("\n");
//     printf("grindlength %d %d %d %d\n", gridlenx, gridleny, gridlenz, gridlent);
//     printf("s->center:%f %f %f\n", s->center[0], s->center[1], s->center[2]);
//     printf("center:%f %f %f\n", centerx, centery, centerz);
//     printf("ix %d %d %d %d\n", ix, iy, iz, it);
//     printf("e: %f %f %f %f\n", e[ix], e[iy], e[iz], e[it]);
//     printf("step %f %f %f %f\n", stepx, stepy, stepz, stept);
//     printf("incr %f %f %f %f \n", incx,  incy, incz, inct);
//     printf("nn=%d ntot=%d nt=%d nnhalf=%f: xoff=%f %d zaehl=%d toffset=%d (%d %d %d)\n",
// 	  (int) nn, (int) ntot, (int) nt, nnhalf,  xoffset, (int) longxoffset,
//  	  (int)zaehler, (int)toffset, (int)(toffset + e[ix] * stepx * gridlenx),
//  	  (int) (toffset + e[iy] * stepy * gridleny),
//  	  (int) (toffset + e[ix]* stepx* gridlenx + e[iy]* stepy* gridleny));
// assert(zaehler < 100);

	      assert((longxoffset<ntot) && (longxoffset>=0) );
	      res[zaehler++] += simuline[longxoffset];
	      xoffset += incx;
	    }
	    yoffset += incy;
	  }
	  zoffset += incz;
	}
	toffset += inct;
      }
    } // n
  } else { 
    // not simugrid, could be time-model!
    // both old and new form included
    double ex=0, ey=0, ez=0, offset, inct;
    long v;
    int i, end;

#define TBMST(UNITYVECTOR, OFFSET, INDEX) {\
      for (n=0; n<tbm->lines; n++){\
        if (tbm->every>0  && (n % tbm->every == 0)) PRINTF("%d \n",n); \
        UNITYVECTOR;\
        internal_DoSimulateRF(&(s->key), 1, simuline);\
        offset= nnhalf - (OFFSET); \
        for (v = nt = i = 0, end=totpoints; nt<gridlent; nt++, end+=totpoints){\
          for (; i < end; i++) {\
            register long index;\
	    index = (long) (offset + INDEX); \
 /* if (true || !((index<ntot) && (index>=0))) { \
    printf("\n%f %f %f (%f %f %f; %1.4f %1.4f %1.4f)\n", \
       s->x[v], s->x[v+1] , s->x[v+2], ex, ey, ez, centerx, centery, centerz); \
    printf("index=%d nn=%d ntot=%d v=%d nt=%d OFF=%f IDX=%f inct=%f e=%d\n", \
       index, nn, ntot, v, nt, OFFSET, INDEX, inct, end); \
       } */ \
            assert((index<ntot) && (index>=0)); \
	    res[i] += simuline[index]; \
	    v += simutimespacedim; \
          }\
          offset += inct;\
        } assert(true);\
	} \
    }

    inct = (double) nn;
    if (s->ce_dim == 1) {
      gridlent = 1;
      totpoints = key->totalpoints;
    } else { // ce_dim==2 nur erlaubt fuer echte Zeitkomponente
      assert(s->ce_dim==2); 
      gridlent = s->key.length[1];
      totpoints = key->spatialtotalpoints;
      //     printf("gridlent=%d totpoints=%d\n", gridlent, totpoints);
//      assert(false);
    }
    assert(gridlent * nn == ntot);
    
//    printf("%d %d %d %d\n", nn, totpoints, s->simuspatialdim, simutimespacedim);    assert(false);

    centery = stepy = incy = centerz = stepz = incz = 0.0;
    switch (s->simuspatialdim) {
	case 3 : 
	    centerz = s->center[2]; 
	    stepz = s->x[XSTEPD[2]];	    // no break;
	case 2 : 
	    centery = s->center[1];
	    stepy = s->x[XSTEPD[1]];	    // no break;
	case 1 : 
	    centerx = s->center[0];
	    stepx = s->x[XSTEPD[0]]; 
	    break;
	default : assert(false);
    }

    if (meth->unimeth==TBM2) {
      double deltaphi;
      deltaphi = PI / (double) tbm->lines;
      phi = deltaphi * UNIFORM_RANDOM; 
      if (s->simuspatialdim==1) {// dim == 1, TBM 2, arbitrary locations
	TBMST(phi += deltaphi; ex=sin(phi), // UNITVECTOR
	      centerx * ex,    // OFFSET
	      s->x[v] * ex);  //  INDEX
      } else {  // dim == 2, TBM 2
	TBMST(phi += deltaphi; ex=sin(phi); ey=cos(phi),
	      centerx * ex + centery * ey,
	      s->x[v] * ex + s->x[v+1] * ey);
      }
    } else { // TBM3, not grid, dim=1 or 2
      assert(meth->unimeth == TBM3);

//      printf("xxxx %d %d %d n", totpoints, gridlent, simutimespacedim);     
//      for (v=0; v<56; v++) printf("%f ", s->x[v]);

      switch (s->simuspatialdim) {
	  case 1 : //see Martin's techrep f. details
	    TBMST(unitvector3D(1, &ex, &ey, &ez),
		  centerx * ex, 
		  s->x[v] * ex);
	    break;
	  case 2: 
	    TBMST(unitvector3D(2, &ex, &ey, &ez),
		  centerx * ex + centery * ey,
		  s->x[v] * ex + s->x[v+1] * ey);
	    break;
	  case 3:
	    TBMST(unitvector3D(3, &ex, &ey, &ez),
		  centerx * ex + centery * ey + centerz * ez,
		  s->x[v] * ex + s->x[v+1] * ey + s->x[v+2] * ez);
	    break;
	  default : assert(false);
      }
    }
  }

  register long i;
  register double InvSqrtNlines;   
  InvSqrtNlines=1.0 / sqrt((double) tbm->lines);
  for(i=0;i<key->totalpoints;i++) {
    res[i] *= InvSqrtNlines; 
    // printf("%d, %f\n", i, res[i]);
  }
}
