/*
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Simulation of a random field by turning bands;
 see RFspectral.cc for spectral turning bands

 Copyright (C) 2001 -- 2011 Martin Schlather, 

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
 
#include <R_ext/Applic.h>
#include "RF.h"
//#include "Covariance.h"
#include "win_linux_aux.h"
#include <R_ext/Utils.h>     

#define MAXNN 100000000.0 /* max number of points simulated on a line */

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


void TBM_destruct(void **S) 
{
  if (*S!=NULL) {
   TBM_storage *s;
    s = *((TBM_storage**)S);

    if (s->simuline != NULL) free(s->simuline);
    if (s->aniso != NULL) free(s->aniso);
    KEY_DELETE(&(s->key));
    free(*S);
    *S = NULL;
  }
}

void TBM_NULL(TBM_storage* s) {
  s->simuline=NULL;
  KEY_NULL(&(s->key));
  s->ce_dim =  0;
}

void SetParamTBMCE(int *action, int *force, double *tolRe, double *tolIm,
		   int *trials, double *mmin, int *useprimes, int *strategy,
		   double *maxmem, int *dependent)
{
  SetParamCE(action, force, tolRe, tolIm, trials,
	     mmin, useprimes, strategy,
	     maxmem, dependent, &(GLOBAL.tbmce), "TBMCE");
  GLOBAL.tbmce.method = 0;
}


void SetParamLines(int *action,int *nLines, double *linesimufactor, 
		   double *linesimustep, int *layers, tbm_param *tbm, 
		   const char *name, double *linesimustep_o, double *linesimufactor_o) 
{
  if (*action) {
    tbm->lines=*nLines;
    if (*linesimufactor < 0.0 || *linesimustep < 0.0) 
	PRINTF("Both %s.linesimustep and %s.linesimufactor and must be non-negative\n", name, name); 
    tbm->linesimufactor=*linesimufactor < 0.0 ? 0.0 : *linesimufactor;
    tbm->linesimustep=*linesimustep < 0.0 ? 0.0 : *linesimustep;    
    if (*linesimufactor!=0.0 && *linesimustep!=0.0 && PL>0 &&
	(*linesimustep_o != *linesimustep || 
	 *linesimufactor_o != *linesimufactor)) {
      if (PL>0) 
	PRINTF("%s.linesimufactor is ignored if %s.linesimustep is positive!\n", 
	       name, name); 
      *linesimustep_o = *linesimustep;
      *linesimufactor_o = *linesimufactor;
      // else, i.e. both are zero, the old values are kept!
    }
    tbm->layers = *layers;
  } else {
    *nLines=tbm->lines;
    *linesimufactor=tbm->linesimufactor;
    *linesimustep=tbm->linesimustep;  
    *layers = tbm->layers;
  }
}

void SetParamTBM2(int *action, int *nLines, double *linesimufactor, 
		  double *linesimustep, int *layers) {
  static double linesimustep_old=-1, linesimufactor_old=-1;
  SetParamLines(action, nLines, linesimufactor, linesimustep,  layers, 
		&(GLOBAL.TBM2), "TBM2", &linesimustep_old, &linesimufactor_old);
}

void SetParamTBM3(int *action,int *nLines, double *linesimufactor, 
		  double *linesimustep,  int *layers) {
  static double linesimustep_old=-1, linesimufactor_old=-1;
  SetParamLines(action, nLines, linesimufactor, linesimustep, layers, 
		&(GLOBAL.TBM3), "TBM3", &linesimustep_old, &linesimufactor_old);
}

void SetParamTBM(int *action, int *tbm_method, double *center, int *points) {
  int d;
  tbmgen_param *gp = &(GLOBAL.TBMgen);
  if (*action) {
   for (d=0; d<MAXTBMSPDIM; d++) gp->center[d] = center[d];
   gp->points = *points;

   // method must be the last one because of the return in the middle
   if ((*tbm_method<=(int) Nothing) && (*tbm_method>=0)) {      
      SimulationType m;
      m=(SimulationType) *tbm_method;
      if ((m==CircEmbed) || (m==Direct) || (m==Nothing)){
	gp->method = m;
	return;
      }
    }
    if (PL>0) {
      PRINTF("Specified method [%d] is not allowed. Method set to `Nothing'\n",
	     *tbm_method);
    }
    gp->method = Nothing;
  } else {
    
    *tbm_method = (int) gp->method;
    for (d=0; d<MAXTBMSPDIM; d++) center[d] = gp->center[d];
    *points = gp->points;
  }
}



// aufpassen auf ( * * 0 )
//               ( * * 0 )
//               ( 0 0 1 )

int init_turningbands(method_type *meth, SimulationType method)
{
  location_type *loc = meth->loc;
  TBM_storage *s;
  globalparam *gp = meth->gp;
  tbm_param *tbm = (method == TBM2) ? &(gp->TBM2) : &(gp->TBM3);
  tbmgen_param *tp = &(gp->TBMgen);
  cov_model *Cov= meth->cov,
    *cov = (Cov->nr == GATTER) ? Cov->sub[0] : Cov,
    *newcov=NULL, 
    *dollar=NULL;
  int tbm_points = tp->points;
  double min[MAXTBMSPDIM], max[MAXTBMSPDIM], Center[MAXTBMSPDIM],
    *tbm_center = tp->center;
  SimulationType  tbm_method = tp->method;
  int tbm_dim = (method==TBM2) ? 2 : 3, // MAXTBMDIM !
    // totalpoints = loc->totalpoints,
    spatialpoints = loc->spatialtotalpoints,
    origdim = loc->timespacedim,
    origspatialdim = loc->spatialdim,
    endaniso = origdim * origdim - 1,
    effectivedim // the isotropic part (i.e. time might be not considered)
    ;
  char tbm2char[] = "tbm2", tbm3char[] = "tbm3";


  int err=NOERROR, d, k;
  double linesimuscale,  diameter;
  bool ce_dim2=false, // just to avoid error messages 
    tbm2num = false; 
  char errorloc_save[nErrorLoc];

  assert(meth->S == NULL);
   /* multiplication of the anisotropy matrix from the right */
  SET_DESTRUCT(TBM_destruct);

#define BUFFER 3.0

  if ((meth->S = malloc(sizeof(TBM_storage)))==0){
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  s = (TBM_storage*) meth->S;
  TBM_NULL(s);
  if (cov->tsdim > MAXTBMSPDIM) {
    err=ERRORMAXDIMMETH; goto ErrorHandling;
  } 


  /****************************************************************/
  /*               Determination of the dimensions                */
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
  
  if (origdim > 4) {
      err = ERRORWRONGDIM;
      goto ErrorHandling;
  }

  if (PL>4) 
    PRINTF("determining first nonzero_pos and TBM simulation dimension...\n");
  s->ce_dim = 0;
  effectivedim = cov->tsdim;
  if (cov->statIn != STATIONARY && cov->statIn != VARIOGRAM) { 
    // in theory it works also for variograms!
    err=ERRORNOVARIOGRAM; goto ErrorHandling;
  }
  if (loc->Time) {
    ce_dim2 = tbm->layers > 0 || cov->isoIn==SPACEISOTROPIC ||
      effectivedim == 1 + tbm_dim;
    effectivedim -= ce_dim2;  
    if (ce_dim2 && tbm->layers < 0) {
      err=ERRORLAYERS; goto ErrorHandling;
    }
  } else {
    ce_dim2 = false;
    if (tbm->layers > 1) {
      err=ERRORLAYERS; goto ErrorHandling;
    }
  }
  if (effectivedim > tbm_dim) {
    err=ERRORWRONGDIM; goto ErrorHandling;
  }
  s->ce_dim = 1 + (int) ce_dim2;
  s->simuspatialdim = effectivedim;

//  { int dd;
//  for (dd=0; dd<12; dd++) printf("dd=%d %f\n", dd, meth->c aniso[d]);
//  }
//  PrintModelInfo(cov);
//  printf("%d %d\n", origdim, cov->tsdim);

  
  s->aniso = getAnisoMatrix(meth);

  /****************************************************************/
  /*          determine length and resolution of the band         */
  /****************************************************************/


 // sort out which finer grid should be used in the spatial domain for
  // low dimensional simulation
  if (PL>4) PRINTF("\ndetermination of the band length...");

  // PrintMethodInfo(meth);

  // minimum distance between the points   
  linesimuscale = -1.0;
  if (tbm->linesimustep  > 0.0) linesimuscale = 1.0 / tbm->linesimustep;
  else if (tbm->linesimufactor > 0) {
    double mindelta;
    mindelta = RF_INF; 
    if (loc->grid) {
      // if linesimufactor, the fine grid scale is proportional to smallest 
      // distance in the grid. 
      for (d=0; d<origdim; d++) {
        if ((loc->xgr[d][XSTEP]<mindelta) && (loc->xgr[d][XSTEP]>0)) {
	  mindelta=loc->xgr[d][XSTEP];
	}
      }
    } else {
      int j, i, d, k0, k1, k2;	
      if (spatialpoints > 50000) {
	  err=ERRORTOOMANYPOINTS; goto ErrorHandling;
	  /* algorithmus kann verbessert werden, so dass diese Fehlermeldung 
	     nicht mehr notwendig ist! */ 
      }
      double diff, dist;
      for (k0=i=0; i<spatialpoints; i++, k0+=origspatialdim) {
	for (k2=j=0; j<i; j++) {
	  k1 = k0;
	  for (dist=0.0, d=0; d<origspatialdim; d++) {
	    diff = loc->x[k1++] - loc->x[k2++]; 
	    dist += diff * diff;
	    }
	  if (dist>0 && dist<mindelta) mindelta=dist; 
	}
      } // if linesimustep==0.0
      if (origdim > origspatialdim) {
	mindelta += loc->T[XSTEP] * loc->T[XSTEP];
      }
      mindelta = sqrt(mindelta);
    }
    linesimuscale = tbm->linesimufactor / mindelta;
  } else {
    if (tbm_points < 4) { err=ERRORTBMPOINTSZERO; goto ErrorHandling; }
  }


  // multiplication der x-skalenparameter so dass letztendlich
  // auf der TBM-Geraden mit Abstand 1 der Punkte simuliert werden kann
  // mit dem Vorteil des einfachen Zugriffs auf die simulierten Werte 


  // untersucht die Originaldaten
  GetMinMax(meth, min, max, Center, MAXTBMSPDIM);
  if (!ISNA(tbm_center[0])) {
    for (d=0; d < origdim; d++) {
      if (!R_FINITE(tbm_center[d])) {
	err=ERRORTBMCENTER; goto ErrorHandling; 
      }
    }
    for (d=0; d<origdim; d++) { // not here no time component!
      Center[d] = tbm_center[d]; // in orig coordinates !
    }
  }
  
  diameter = GetDiameter(Center, min, max, s->aniso, loc->timespacedim,
			 effectivedim); 

  for (d=0; d<origdim; d++) {
    s->center[d] = Center[d];
    if (loc->grid || (loc->Time && d==loc->timespacedim-1)) 
      s->center[d] -= loc->xgr[d][XSTART];
  }

//  printf("diam %f %d %f %f\n", diameter, tbm_points, BUFFER, linesimuscale );

  if (linesimuscale < 0) {
      linesimuscale = (tbm_points - BUFFER) / diameter;
  }


  diameter = trunc(BUFFER + diameter * linesimuscale); // in pts

  // printf("diam %f %f\n", diameter, linesimuscale);

  s->linesimuscale = linesimuscale;

  if (tbm_points > 0) { 
    if (diameter > tbm_points) { // never happens if linesimuscale has been < 0
      PRINTF("tbm points minimum=%f, but tbm_points is %d.\n", 
	     diameter, tbm_points);
      err=ERRORTBMPOINTS; goto ErrorHandling;
    } else diameter = (double) tbm_points;
  }
  if (diameter > MAXNN) {err=ERRORNN; goto ErrorHandling;}

// printf("diam %f %f\n", diameter, linesimuscale);


  // only needed for nongrid !
  s->spatialtotalpts = loc->totalpoints;
  if (s->ce_dim == 2) {
    s->spatialtotalpts /= loc->length[loc->timespacedim - 1];
  }
 


  /****************************************************************/
  /*          set up the covariance function on the line          */
  /****************************************************************/
  if (PL>4) PRINTF("\nsetting up of the covariance funciton on the line...\n");
  addModel(&newcov, method==TBM2 ? getmodelnr(tbm2char) : getmodelnr(tbm3char));
  newcov->tsdim = newcov->xdim = s->ce_dim;
  if (method==TBM3) { // tbm3 has an additional parameter, but not tbm2
      cov->p[0] = (double*) malloc(sizeof(double)); 
      cov->p[0][0] = (double) s->ce_dim;
  }
  for (k=0; k<Nothing; k++) newcov->user[k] = PREF_NONE;
  if (tbm_method == Nothing) {
    newcov->user[CircEmbed] = newcov->user[Direct] = PREF_BEST;
  } else {
    newcov->user[tbm_method] = PREF_BEST;
  }
  covcpy(newcov->sub + 0, 
	 (cov->nr >= GATTER && cov->nr <= LASTGATTER) ? cov->sub[0] : cov, true, false);
  newcov->sub[0]->calling = newcov;

  addModel(&newcov, GATTER); // vor tbm2/3
  newcov->calling=NULL;
  newcov->tsdim = 1;
  newcov->xdim = 1;
  newcov->statIn = STATIONARY;
  newcov->isoIn = ISOTROPIC;


  newcov->tsdim = newcov->xdim = s->ce_dim;
  addModel(newcov->sub, GATTER); 
  addModel(newcov->sub, DOLLAR);
  dollar = newcov->sub[0];
  dollar->p[DVAR] = (double*) malloc(sizeof(double)); 
  dollar->p[DVAR][0] = 1.0;
  dollar->ncol[DVAR] = dollar->nrow[DVAR] = 1; 
  dollar->p[DANISO] = (double*) malloc(sizeof(double) * s->ce_dim * s->ce_dim);
  dollar->ncol[DANISO] = dollar->nrow[DANISO] = s->ce_dim; 
  dollar->p[DANISO][0] = 1.0 / linesimuscale;

  assert(newcov->xdim == s->ce_dim);
  if (s->ce_dim == 2) {
    dollar->p[DANISO][1] = dollar->p[DANISO][2] = 0.0;
    dollar->p[DANISO][3] = s->aniso[endaniso] * loc->T[XSTEP];
  } else { 
    //  assert(s->ce_dim == 1);
  }
 
  //
//  PrintModelInfo(Cov);
  // PrintModelInfo(newcov);


  if ((err =  CovList[newcov->nr].check(newcov)) != NOERROR) goto ErrorHandling;
  DeleteGatter(&newcov);


  // PrintModelInfo(newcov);  assert(false);
  

  /****************************************************************/
  /*                set up the process on the line                */
  /****************************************************************/

  if (PL>4) PRINTF("\nsetting up the process on the line...");

  if (method==TBM2) {
    tbm2num = cov->tbm2num || (ce_dim2 && cov->isoIn == ISOTROPIC);
    if (tbm2num && PL > 1)
      PRINTF("\tnumerical evaluation of the TBM operator\n");
  }
  

  // xline[XEND]: number of points on the TBM line
  double xline[3];
  xline[XSTART] = 1.0;
  xline[XSTEP] = 1.0;
  xline[XEND] = (double) diameter; // diameter > MAXNN must be first since
  //                                risk of overflow !

  if (diameter < 10) {
    char MSG[255];
    PRINTF("\n%s%d%s%d%s",
	    "too few points on the line -- modify TBM",
	    (method == TBM2) ? 2 : 3,
	    ".linesimufactor or TBM",
	    (method == TBM2) ? 2 : 3,
	    ".linesimustep in RFparameters");
    ERR("too few points on the line");
  }

  // CovList.tbm_method is not used yet
  // if (TBM_METHOD==Nothing) {tbmmethod[i]=CovList[*covnr].tbm_method;}

  
  globalparam global;
  memcpy(&global, &GLOBAL, sizeof(globalparam));
  memcpy(&(global.ce), &(GLOBAL.tbmce), sizeof(ce_param));
 
  strcpy(errorloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "%s%s: ", errorloc_save, METHODNAMES[method]);
  err = internal_InitSimulateRF(xline, loc->T, 1, 3, true, 
				ce_dim2 /* Time */,
				DISTR_GAUSS, &(s->key),
				&global, tbm->lines, &newcov);
  // note that, by internal_Init, newcov is set to NULL

  strcpy(ERROR_LOC, errorloc_save);
  if (err != NOERROR) goto ErrorHandling; 

  
  // printf("calloc %d\n", s->key.loc.totalpoints);

  if ((s->simuline=(res_type *) calloc(s->key.loc.totalpoints, 
				     sizeof(res_type)))==0){
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;}
  if (PL>4) PRINTF("\ntbm is now initialized.\n");
 
  return NOERROR;
  
 ErrorHandling:
  return err;
}

int init_turningbands2(method_type *meth) {
 return(init_turningbands(meth, TBM2));
}

int init_turningbands3(method_type *meth) {
  return(init_turningbands(meth, TBM3));
}

void GetE(SimulationType unimeth, TBM_storage *s, int dim,  bool Time, 
	  double *phi, double deltaphi,	  
	  double *offset, double *ex, double *ey, double *ez, double *et) {
  double sube[MAXTBMSPDIM], e[MAXTBMSPDIM];
  int d, j, idx, k,
      spatialdim = s->simuspatialdim;
//      total = spatialdim * dim;
  // dim is the users's dim of the field
  // this might be reduced by zonal anisotropies leading to spatialdim

  // for debugging only
  for(d=0; d<MAXTBMSPDIM; d++) sube[d] = e[d] = RF_NEGINF;

  if (unimeth==TBM2) {
    // no choice between random and non-random
    (*phi) += deltaphi;
    sube[0] = sin(*phi); 
    sube[1] = cos(*phi);
  } else {
    unitvector3D(spatialdim, sube, sube +1, sube + 2);
//    printf("spatdim=%d %d  sube=%f  %f %f\n", 
//	   spatialdim, s->simuspatialdim, sube[0], sube[1], sube[2]);
  }

  *offset = 0.5 * (double) s->key.loc.length[0];
//    printf("offset %d\n", *offset);
  for (d=0; d<dim; d++) e[d] = 0.0;
  for (k=j=0; j<spatialdim; j++) {
    for (d=0; d<dim; d++) {
      e[d] += sube[j] * s->aniso[k++];
//      printf("e,s,a, (%d, %d) %f     %f     %f\n",
//	     d, j, e[d], sube[j], s->aniso[j + d]);
    }
  }
  for (d=0; d<dim; d++) {
    e[d] *= s->linesimuscale;
    *offset -= s->center[d] * e[d];
    //  
  }

//  for (d=0; d<dim; d++) {
//      printf("e, lines %d offset=%f e=%f %d center=%f %f %d\n", 
//	     d, *offset, e[d], s->key.loc.length[0], 
//	     s->center[d],  s->linesimuscale, s->ce_dim);
// }
 
    
  idx = spatialdim;
  if (Time && s->ce_dim==1) {
    *et = e[--idx];
  }
  switch (idx) {
      case 4 : assert(false);
      case 3 : *ez = e[--idx];
      case 2 : *ey = e[--idx];
      case 1 : *ex = e[--idx];
	break;
      default: assert(false);
  }
}

void do_internal_tbm(method_type *meth, res_type *res, SimulationType unimeth,
		     tbm_param *tbm) { 
  location_type *loc = meth->loc;
  cov_model *cov = meth->cov;
  TBM_storage *s = (TBM_storage*) meth->S;
  long nn = s->key.loc.length[0],
    ntot = s->key.loc.totalpoints;
  int dim = cov->tsdim,
    spatialdim = s->simuspatialdim,
    every = meth->gp->general.every;
    
  
  res_type
    *simuline = s->simuline;
  
  double phi, deltaphi, 
      incx=0.0, incy=0.0, incz=0.0, inct=0.0, stepx, stepy, stepz, stept;
  long n, totpoints;
  int nt, idx, gridlent;
    
  for (n=0; n<loc->totalpoints; n++) res[n]=0.0; 
  deltaphi = PI / (double) tbm->lines; // only used for tbm2
  phi = deltaphi * UNIFORM_RANDOM;     // only used for tbm2

  if (loc->grid) {
    int nx, ny, nz, gridlenx, gridleny, gridlenz;
    long zaehler;
    double xoffset, yoffset, zoffset,  toffset;

    stepx = stepy = stepz = stept =  0.0;
    gridlenx = gridleny = gridlenz = gridlent = 1;
  
    idx = dim;
    if (s->ce_dim==2) {
      gridlent = loc->length[--idx]; // could be one !!
      stept = loc->xgr[idx][XSTEP];	    
      inct = (double) nn;
    } else if (loc->Time) {
      gridlent = loc->length[--idx]; // could be one !!
      stept = loc->T[XSTEP];	    
    }

    switch (idx) {
	case 4 : 
	  assert(false);
	case 3 : 
	  gridlenz = loc->length[--idx];
	  stepz = loc->xgr[idx][XSTEP];	  
	  // no break;
	case 2 : 
	  gridleny = loc->length[--idx];
	  stepy = loc->xgr[idx][XSTEP];	 
	  // no break;
	case 1 : 
	  gridlenx = loc->length[--idx];
	  stepx = loc->xgr[idx][XSTEP];	  
	  break;
	default : assert(false);
    }
         
    for (n=0; n<tbm->lines; n++) {
      R_CheckUserInterrupt();
      if (every>0  && (n % every == 0)) PRINTF("%d \n",n);
      GetE(unimeth, s, dim, loc->Time, &phi, deltaphi,
	   &toffset, &incx, &incy, &incz, &inct);
      incx *= stepx;
      incy *= stepy;
      incz *= stepz;
      if (s->ce_dim == 1) inct *= stept;
      internal_DoSimulateRF(&(s->key), 1, simuline);

//      if (n<5) {
//	  int i;
//	  for (i=0; i<8; i++) printf("%2.3f ", simuline[i]);
//	  printf("%d \n", ntot);
//      }

      zaehler = 0;
      for (nt=0; nt<gridlent; nt++) {
	zoffset = toffset;
	for (nz=0; nz<gridlenz; nz++) {	  
	  yoffset = zoffset;
	  for (ny=0; ny<gridleny; ny++) {	  
	    xoffset = yoffset;
	    for (nx=0; nx<gridlenx; nx++) {
              long longxoffset = (long) xoffset;
	      if (!((longxoffset<ntot) && (longxoffset>=0))) {
		PRINTF("xx ntot %ld %ld %ld -- offsets %f %f %f %f\n x:%d %d y:%d %d z:%d %d t:%d %d\n", 
		       ntot, longxoffset, zaehler,
		       xoffset, yoffset, zoffset, toffset,
		       nx, gridlenx, ny, gridleny, nz, gridlenz, nt, gridlent
		    );
//		continue;
		assert((longxoffset<ntot) && (longxoffset>=0) );
	      }
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
    // not grid, could be time-model!
    // both old and new form included
    double offset;
    long v;
    int i, end;

    PrintMethodInfo(s->key.meth);

#define TBMST(INDEX) { \
      for (n=0; n<tbm->lines; n++){ \
        if (every>0  && (n % every == 0)) PRINTF("%d \n",n); \
        GetE(unimeth, s, dim, loc->Time, &phi, deltaphi, \
	   &offset, &incx, &incy, &incz, &inct);  \
        internal_DoSimulateRF(&(s->key), 1, simuline); \
        for (v = nt = i = 0, end=totpoints; nt<gridlent; nt++, end+=totpoints){ \
          for (; i < end; i++) { \
            long index; \
	    index = (long) (offset + INDEX); \
  if (!((index<ntot) && (index>=0))) { \
    PRINTF("\n%f %f %f (%f %f %f))\n", \
       loc->x[v], loc->x[v+1] , loc->x[v+2], incx, incy, incz); \
    PRINTF("index=%d nn=%d ntot=%d v=%d nt=%d OFF=%f IDX=%f inct=%f e=%d\n", \
       index, nn, ntot, v, nt, offset, INDEX, inct, end); \
       }  \
            assert((index<ntot) && (index>=0)); \
	    res[i] += simuline[index]; \
	    v += dim; \
          } \
          offset += inct; \
        } assert(true); \
      } \
    }

    if (s->ce_dim == 1) {
      gridlent = 1;
      totpoints = loc->totalpoints;
      inct = RF_NAN; // should be set by GetE
    } else { // ce_dim==2 nur erlaubt fuer echte Zeitkomponente
      assert(s->ce_dim==2); 
      gridlent = s->key.loc.length[1];
      totpoints = s->spatialtotalpts;
      inct = (double) nn;
    }
    assert(gridlent * nn == ntot);
    
    switch (spatialdim) {
	case 1 : //see Martin's techrep f. details
	  TBMST(loc->x[v] * incx);
	  break;
	case 2: 
	  TBMST(loc->x[v] * incx + loc->x[v+1] * incy);
	  break;
	case 3:
	  TBMST(loc->x[v] * incx + loc->x[v+1] * incy + loc->x[v+2] * incz);
	  break;
	default : assert(false);
      }
  } // end not grid

//  {double z; int i; for(i=0, z=0.0; i<ntot; i++) z+=simuline[i]*simuline[i];
//      printf("%f %f %f\n", 
//	     simuline[0], z / ntot, res[10] / sqrt((double) tbm->lines));
//  }

  long i;
  double InvSqrtNlines;   
  InvSqrtNlines = 1.0 / sqrt((double) tbm->lines);
  for(i=0;i<loc->totalpoints;i++) {
      res[i] *= (res_type) InvSqrtNlines; 
  }

//  {double z; int i; for(i=0, z=0.0; i<loc->totalpoints; i++) z+=res[i]*res[i];
//     printf("%f %f %f\n", 
//	     simuline[0], z /loc->totalpoints , res[10] );
//  }
  
}

void do_tbm2(method_type *meth, res_type *res) {
  do_internal_tbm(meth, res, TBM2, &(meth->gpdo->TBM2));
} 

void do_tbm3(method_type *meth, res_type *res) {
  do_internal_tbm(meth, res, TBM3, &(meth->gpdo->TBM3)); 
} 
