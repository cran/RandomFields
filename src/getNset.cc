

/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 (library for simulation of random fields)

 Copyright (C) 2001 -- 2011 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
RO
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <R.h>
#include <Rdefines.h>
#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
 
#include <string.h>
#include "RF.h"
#include "primitive.h"
#include <R_ext/Linpack.h>



void KEY_DELETE(key_type *key) {
  KEY_DELETE_NOTTREND(key);
  TREND_DELETE(&(key->trend));
}

void KEY_NULL(key_type *key)
{
  KEY_NULLNOTTREND(key);
  TREND_NULL(&(key->trend));
}


void KEY_DELETE_NOTTREND(key_type *key)
{
  LOC_DELETE(&(key->loc));
  METHOD_DELETE(&(key->meth));

  // PrintModelInfo(key->cov);

  if (key->cov!=NULL) COV_DELETE(&(key->cov));
//  KEY_NULLNOTTREND(key);
}

void KEY_NULLNOTTREND(key_type *key) {
  memcpy(&(key->gp), &GLOBAL, sizeof(globalparam)); 
  memcpy(&(key->gpdo), &GLOBAL, sizeof(globalparam));
  simu_type *simu = &(key->simu);
  simu->stop = simu->active = false;
  simu->expected_number_simu = 0;
  LOC_NULL(&(key->loc));
  // trend left as is
  key->cov = NULL;
  key->meth= NULL;
 }

void TREND_DELETE(trend_type *trend)
{ 
 if (trend->TrendFunction!=NULL) free(trend->TrendFunction); 
  if (trend->LinearTrend!=NULL) free(trend->LinearTrend); 
//  TREND_NULL(trend);
}

void TREND_NULL(trend_type *trend)
{
  trend->mean = 0.0;
  trend->TrendModus = -1;
  trend->lTrendFct = 0;
  trend->lLinTrend = 0;
  trend->TrendFunction = NULL;
  trend->LinearTrend = NULL;
}

void LOC_NULL(location_type *loc) {
  int d;
  loc->spatialdim = loc->timespacedim = -1;
  loc->totalpoints = loc->spatialtotalpoints = 0;
  for (d=0; d<MAXSIMUDIM; d++) {
    loc->xgr[d]=NULL;
    loc->length[d] = -1;
  }
  loc->x = NULL;
  loc->T[0] = loc->T[1] = loc->T[2] = 0.0; 
}
     
void LOC_DELETE(location_type *loc) {
  if (loc->x!=NULL) {
    free(loc->x); 
    loc->x = NULL;
  } else 
  if (loc->xgr[0] != NULL) {
    free(loc->xgr[0]);
    loc->xgr[0] = NULL;
  }
}

void METHOD_DELETE(method_type **Meth) {
  method_type *meth = *Meth;
  int  i;

  if (meth==NULL) {
    return;
  }

  for (i=0; i<MAXSUB; i++) {
    if (meth->sub[i] != NULL) METHOD_DELETE(meth->sub + i);
  }
  if (meth->destruct!=NULL) meth->destruct(&(meth->S));
  if (meth->caniso != NULL) free(meth->caniso);
  if (meth->cproj != NULL) free(meth->cproj);
  if (meth->loc != NULL) {
    if (meth->loc->Time) {
      if (meth->space != NULL) free(meth->space);
      if (meth->sptime != NULL) free(meth->sptime);
    } else {
      assert(meth->space == meth->sptime);
      if (meth->space != NULL) {
	free(meth->space);
      }
    }
  }
  free(*Meth);
  *Meth = NULL;
}

void METHOD_NULL( method_type *meth) {
  int  i;
  meth->gp = meth->gpdo = NULL; // never leave unset -- used in METHOD_SET
  meth->simu = NULL;
  meth->loc = NULL;

  // method
  meth->nr = -Forbidden-1;
  meth->compatible = false;
  meth->nsub = 0;
  for (i=0; i<MAXSUB; i++) meth->sub[i] = NULL;
  meth->domethod = NULL;
  meth->destruct = NULL;
  meth->S = NULL;

  meth->cov = NULL;
  meth->caniso = NULL;
  meth->cscale = 1.0;
  meth->cvar = -1.0;
  meth->cproj = NULL;
  meth->xdimout = -1;
  meth->type = TypeAny;
  meth->hanging = NULL;
  meth->space = meth->sptime = NULL;
  // grani
}


// Achtung! removeonly muss zwingend mit dem originalen SUB + i aufgerufen
// werden. Der **Cov ueberschrieben wird. Ansonsten muessen die
// Kommentarzeilen geloescht werden !
void removeOnly(cov_model **Cov) {
  cov_model *cov = *Cov,
    *next = cov->sub[0];
 
  if (cov->calling == NULL) next->calling = NULL;
  else {
    cov_model *calling = cov->calling;
//    for (i=0; i<MAXSUB; i++) if (calling->sub[i] == cov) break;
//    assert (i<MAXSUB);
//    calling->sub[i] = cov->sub[0];
    next->calling = calling;
  }
  *Cov = next;
  COV_DELETE_WITHOUTSUB(&cov);
}


void COV_DELETE_WITHOUTSUB(cov_model **Cov) {
  cov_model *cov = *Cov;
  int i, j,
    last = (cov->nr < 0) ? MAXPARAM : CovList[cov->nr].kappas; 
  if (cov->q != NULL) {
    free(cov->q);
    cov->q = NULL;
    cov->qlen = 0;
  }

  if (cov->MLE != NULL) free(cov->MLE);

  // print("last = %d  %s\n", last, CovList[cov->nr].name);
  
  for (i=0; i<last; i++) {
    if (cov->p[i] != NULL) {
      if (CovList[cov->nr].kappatype[i] >= LISTOF) {
	  listoftype *list = (listoftype*) (cov->p[i]);
	for (j=0; j<cov->nrow[i]; j++) {
	    free(list->p[j]);  
	}
      }
      free(cov->p[i]); 
      cov->p[i] = NULL;
      cov->ncol[i] = cov->nrow[i] = 0;
    }	
  }
  // important check in combination with above; can be easily removed or 
  // generalised  !!!!
  for (; i<MAXPARAM; i++) {
    if (cov->p[i] != NULL) {
	// print("last=%d current=%d max=%d\n", last, i, MAXPARAM);
	// PrintModelInfo(*Cov); //
      // assert(cov->p[i] == NULL);
    }
  }
  free(*Cov);
  *Cov = NULL;
}


void COV_DELETE(cov_model **Cov) { 
  cov_model *cov = *Cov;
  int i;

  if (cov == NULL) {
    warning("*cov is NULL in COV_DELETE");
    return;
  }

  for (i=0; i<MAXSUB; i++) { // cov->sub[i] // seit 1.10.07: luecken erlaubt
    //                         bei PSgen !
    if (cov->sub[i] != NULL) COV_DELETE(cov->sub + i);
  }
//  for (; i<MAXSUB; i++) assert(cov->sub[i] == NULL);
  COV_DELETE_WITHOUTSUB(Cov);
}



void COV_NULL(cov_model *cov) {
  int i;
  for (i=0; i<MAXPARAM; i++) {
    cov->p[i] = NULL;
    cov->ncol[i] = cov->nrow[i] = 0;
  }
  cov->q = NULL;
  cov->qlen = 0;
  for (i=0; i<MAXSUB; i++) cov->sub[i] = NULL;
  cov->nsub = 0;
  cov->calling = NULL;
  cov->tsdim = cov->xdim =cov->maxdim = UNSET;
  cov->normalmix = cov->finiterange = cov->diag = cov->semiseparatelast
    = cov->separatelast = false;
  for (i=0; i<=Nothing; i++) cov->user[i] = cov->pref[i] = PREF_BEST;
  // mpp und extrG muessen auf NONE sein !
  for (; i<=Forbidden; i++) cov->user[i] = cov->pref[i] = PREF_NONE;
  cov->derivatives = -1;
  cov->tbm2num = false;
  cov->spec.density = NULL;
  cov->statIn = STAT_MISMATCH;
  cov->isoIn = ISO_MISMATCH;
  cov->vdim = M_MISMATCH;
  cov->nr = -1;
  cov->MLE = NULL;
  cov->anyNAdown = cov->anyNAscaleup = TriBeyond;
  cov->manipulating_x=false;
  cov->spec.nmetro = -1;
  cov->spec.sigma = -1.0;
  cov->hess = false;
}

void addModel(cov_model **pcov, int covnr) {
  cov_model *cov;
  int i;
  if (covnr >= GATTER && covnr <= LASTGATTER && 
      *pcov != NULL && (*pcov)->nr == GATTER)
    return; // keine 2 Gatter hintereinander
  if ((cov = (cov_model*) malloc(sizeof(cov_model)))==0)
    ERR("memory allocation error");
  COV_NULL(cov);
  cov->nr = covnr;
  cov->nsub = 1;
  if (*pcov!=NULL) {
    cov->calling = (*pcov)->calling;
    (*pcov)->calling = cov;
  } 
  cov->sub[0] = *pcov;
  *pcov = cov;
  cov->tbm2num = NOT_IMPLEMENTED;

  if (covnr == GATTER && *pcov != NULL)
      for (i=0; i<=Forbidden; i++) {
	  cov->user[i] = cov->sub[0]->user[i];
	  cov->pref[i] = cov->sub[0]->pref[i];
      }
}


int loc_set(double *x, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int lx,  bool Time, bool grid,
	    location_type *loc, globalparam *gp) {
  int d, i, j;
  unsigned long totalBytes;
  // preference lists, distinguished by grid==true/false and dimension
  // lists must end with Nothing!
  
  loc->timespacedim = spatialdim + (int) Time;
  loc->spatialdim = spatialdim;
  loc->lx = lx;
  if (spatialdim<1 ||  loc->timespacedim>MAXSIMUDIM) return ERRORDIM;

  loc->grid= grid;
  // coordinates
  assert(loc->xgr[0]==NULL);
  assert(loc->x == NULL);


  totalBytes =  sizeof(double) * lx * spatialdim;
  if (loc->grid) {
    if ((loc->xgr[0]=(double*) malloc(totalBytes))==NULL){
      return ERRORMEMORYALLOCATION; 
    }
    memcpy(loc->xgr[0], x, totalBytes);
    for (d=1; d<spatialdim; d++) loc->xgr[d]= &(loc->xgr[0][d * lx]);
    for (; d<MAXSIMUDIM; d++)  loc->xgr[d]=NULL;
  } else {
    if ((loc->x=(double*) malloc(totalBytes))==NULL){
      return ERRORMEMORYALLOCATION; 
    }
    for (d=0; d<spatialdim; d++) loc->xgr[d] = x + d * lx;
    for (j=i=0; i<lx; i++) {
      for (d=0; d<spatialdim; d++) {
	loc->x[j++] = loc->xgr[d][i];
      }
    }
    for (d=0; d<MAXSIMUDIM; d++) loc->xgr[d]=NULL;
  }
    
  //Time
  if ((loc->Time = Time)) {
    memcpy(loc->T, T, sizeof(double) * 3);
    //if (loc->grid) 
    loc->xgr[loc->spatialdim] = loc->T;
  }

  if (loc->grid) { 
    if ((lx!=3) || 
	(InternalGetGridSize(loc->xgr, loc->timespacedim, loc->length))) {
      PRINTF("%d\n", lx);
      ERR("problem with the coordiates");
      return ERRORCOORDINATES; 
    }

    for (loc->spatialtotalpoints=1, d=0; d<loc->spatialdim; d++) {
      loc->spatialtotalpoints *= loc->length[d];
    }
    loc->totalpoints = loc->spatialtotalpoints;
    if (loc->Time) loc->totalpoints *= loc->length[loc->spatialdim];
  } else { // not grid
    loc->totalpoints = loc->spatialtotalpoints = loc->length[0]= lx;
    if (loc->Time) {
      double *Tx[MAXSIMUDIM];
      Tx[0] = loc->T;
      if (InternalGetGridSize(Tx, 1 , loc->length + loc->spatialdim)) {
	return ERRORCOORDINATES; 
      }
      loc->totalpoints *= loc->length[loc->spatialdim];
    }
    // just to say that considering these values does not make sense
    for (d=1; d<loc->spatialdim; d++) loc->length[d]=0; 
    // correct value for higher dimensions (never used)
  }
  for (d=loc->timespacedim; d<MAXSIMUDIM; d++)
      loc->length[d] = (int) RF_NAN; // 1

  return NOERROR;
}




void Aniso(double *x, double *aniso, int origdim, int dim, double *newx) {
  double dummy;
  int i, j, k;

  // hier kann die Geschwindigkeit verbessert werden,
  // indem die Typen unterschieden werden
  if (aniso == NULL) {
      assert(dim == origdim);
      memcpy(newx, x, sizeof(double) * origdim);
  } else {
    for (k=i=0; i<dim; i++) {
      dummy = 0.0;
      for (j=0; j<origdim; j++) {
        dummy += x[j] * aniso[k++];
      }
      newx[i] = dummy;
    }
  }
}



void analyse_matrix(double *aniso, int row, int col,
		    bool *diag, bool *quasidiag, int *idx,
		    bool *semiseparatelast,
		    bool *separatelast) {
  // diag, falls durch einfuegen von spalten diag-Matrix erhalten
  // werden kann
  
  // see also Type -- can it be unified ????

  /* 
       -> i
   |  * 0 0
   j  0 0 *
      0 * 0
  */
    bool notquasidiag=true, *taken;
  int j, k, startidx, i, last;

  if (aniso == NULL) {
    *diag = *quasidiag = *separatelast = *semiseparatelast;
    for (i=0; i<col; i++) idx[i] = i;
    return;
  }

  
  taken = (bool *) malloc(row * sizeof(bool));

  for (j=0; j<row; j++) {
    taken[j]=false;
    idx[j] = -1;
  }
  for (k=startidx=i=0; i<col; i++) {
    for (j=0; j<row; j++, k++) if (aniso[k] != 0.0) break;
    if (j < row) {
      if ((notquasidiag = taken[j])) break;
      taken[j] = true;
      idx[j] = i;
      for (j++, k++ ; j<row; j++) {
	if ((notquasidiag = aniso[k++] != 0.0)) break;
      }
    }
    if (notquasidiag)  break;
  }
  if ((*diag = *quasidiag = !notquasidiag)) {
    if (idx[0] == -1) idx[0] = 0;
    for (j=1; j<row; j++) {
      if (idx[j] <= idx[j-1]) {
	if (idx[j] == -1) idx[j] = idx[j-1] + 1; else break; 
      }
    }
    *diag = j >= row;
  }
  if (!(*semiseparatelast = *diag)) {
    last = col * row - 1;
    for (k=last - row + 1; k<last; k++)
      if (!(*separatelast = aniso[k++] == 0.0)) break;
  }
  if (!(*separatelast = *semiseparatelast)) {
    // entered here only if last is set!
    for (k=row - 1; k<last; k+=row)
      if (!(*separatelast = aniso[k++] == 0.0)) break;
  }

  free(taken);
}

void addkappa(int i, const char *n, SEXPTYPE t) {
  cov_fct *C = CovList + currentNrCov - 1;
  assert(n[0] != '\0' && 
	 n[0] != 'k' && // reserved for standard parameter names
	 n[0] != 'u' && // reserved for standard submodel names
	 (n[0] != 'm' || n[1] != 'e') // reserved for me-thod
    );
  assert(i < MAXPARAM);
  strcopyN(C->kappanames[i], n, PARAMMAXCHAR);
  C->kappatype[i] = t;
}

void kappanames(const char* n1, SEXPTYPE t1) {
  addkappa(0, n1, t1);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2) {
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3) {
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4) {
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
  addkappa(3, n4, t4);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5) {
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
  addkappa(3, n4, t4);
  addkappa(4, n5, t5);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6) {
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
  addkappa(3, n4, t4);
  addkappa(4, n5, t5);
  addkappa(5, n6, t6);
}

/*
real    0m47.752s
user    0m47.479s
sys     0m0.204s
make: *** [vrun] Error 1

real    0m47.996s
user    0m47.559s
sys     0m0.364s


real    4m18.616s
user    4m12.044s
sys     0m2.284s
make: *** [vrun] Error 1

real    4m25.319s
user    4m17.472s
sys     0m3.008s
schlather@DD787F2J:~/R/RF> time make v
echo -e "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"

*/

void addsub(int i, const char *n) {
  cov_fct *C = CovList + currentNrCov - 1;
  assert(n[0] != 'k' && n[0] != 'u');
  assert(i < MAXSUB);
  strcopyN(C->subnames[i], n, PARAMMAXCHAR);
}

void subnames(const char* n1) {
  addsub(0, n1);
}
void subnames(const char* n1, const char* n2) {
  addsub(0, n1);
  addsub(1, n2);
}
void subnames(const char* n1, const char* n2, const char* n3) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
}
void subnames(const char* n1, const char* n2, const char* n3, const char* n4) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
  addsub(3, n4);
}
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
  addsub(3, n4);
  addsub(4, n5);
}
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5, const char* n6) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
  addsub(3, n4);
  addsub(4, n5);
  addsub(5, n6);
}


void ErrCov(double *x, cov_model *cov, double *v) {
    // int *f; PRINTF("%d\n", f[0]);
  error("unallowed or undefined call of function");
}
void ErrCovNonstat(double *x, double *y, cov_model *cov, double *v) {
  error("unallowed or undefined call of non-stationary function");
}

char InternalName[]="-";
void kappasize1(int i, cov_model *cov, int *nrow, int *ncol) {
  *nrow = *ncol = 1;
}
void createmodel(const char *name, int kappas, size_fct kappasize,	
		 stationary_type stationary, isotropy_type isotropy,
		 checkfct check, rangefct range, init_meth initmethod,
		 do_meth domethod, int vdim, pref_shorttype pref) {
  int i;

  if (currentNrCov==-1) InitModelList(); 
  assert(CovList!=NULL);
  assert(currentNrCov>=0);
  if (currentNrCov>=MAXNRCOVFCTS) {
    PRINTF("Error. Model list full.\n");
    assert(false);
  }
  
  //printf("%d\n", PL);
  if (PL > 8) PRINTF("%d %s vdim=%d\n", currentNrCov, name, vdim); 

  strcopyN(CovList[currentNrCov].name, name, MAXCHAR);
  strcopyN(CovNames[currentNrCov], name, MAXCHAR);
  assert(strcmp(InternalName, name));

  CovList[currentNrCov].name[MAXCHAR-1]='\0';
  if (strlen(name)>=MAXCHAR) {
    PRINTF("Warning! Covariance name is truncated to `%s'",
	   CovList[currentNrCov].name);
  }

  cov_fct *C = CovList + currentNrCov;
  assert(kappas >= 0);
  C->kappas = kappas;
  C->minsub = C->maxsub = 0;
  C->primitive = true;
  memcpy(C->pref, pref, sizeof(pref_shorttype));

  C->stationary = stationary; 
  C->isotropy = isotropy;
  C->vdim = vdim;
  for (i=0; i<kappas; i++) {
    sprintf(C->kappanames[i], "k%d", i); // default (repeated twice)
    C->kappatype[i] = REALSXP;
  }
  C->kappasize = (kappasize == NULL) ? kappasize1 : kappasize;
  C->range= range; 
  assert(range != NULL);
  C->check = check;
  assert(C->check != NULL);
  for (i=0; i<SimulationTypeLength; i++)
    C->implemented[i] = NOT_IMPLEMENTED;

  C->internal = false;
  C->naturalscale=NULL;
  C->invSq = C->cov = C->D = C->D2 = C->tbm2 = C->D3 = C->D4 = ErrCov;
  C->derivs = -1;
  C->nonstat_cov = ErrCovNonstat;
  C->nabla = C->hess =  ErrCov;
  C->coinit = C->ieinit = NULL;
  C->alternative = NULL;

  C->spectral=NULL;
  C-> initspectral=NULL;

  C->mppinit = NULL;
  C->randomcoin = NULL;
  C->mppget = NULL;
  C->mppgetstat = NULL;
  C->mpplocations = MPP_MISMATCH;
  C->avef = NULL;
  C->avelogg = NULL;
  C->mineigenvalue = NULL;

  C->drawlogmix = NULL;
  C->logmixweight = NULL;

  C->hyperplane=NULL;
//  C->initspecial=NULL;
//  C->special=NULL;
  
  C->initmethod=initmethod;
  C->domethod=domethod;

  currentNrCov++;
}

int IncludePrim(const char *name, int kappas, 
		 stationary_type stationary, isotropy_type isotropy,	
		 checkfct check, rangefct range) {  
    createmodel(name, kappas, NULL, stationary, isotropy, check, range,
		NULL, NULL, UNIVARIATE, PREF_ALL);
    return currentNrCov - 1;
}
int IncludePrim(const char *name, int kappas, size_fct kappasize,
		 stationary_type stationary, isotropy_type isotropy,	
		 checkfct check, rangefct range) {  
    createmodel(name, kappas, kappasize, stationary, isotropy, check, range,
		NULL, NULL, UNIVARIATE, PREF_ALL);
    return currentNrCov - 1;
}

int IncludePrim(const char *name, int kappas, 
		 stationary_type stationary, isotropy_type isotropy,	
		 checkfct check, rangefct range, int vdim) {  
    createmodel(name, kappas, NULL, stationary, isotropy, check, range,
		NULL, NULL, vdim, PREF_ALL);
    return currentNrCov - 1;
}
int IncludePrim(const char *name, int kappas, size_fct kappasize,
		 stationary_type stationary, isotropy_type isotropy,	
		 checkfct check, rangefct range, int vdim) {  
    createmodel(name, kappas, kappasize, stationary, isotropy, check, range,
		NULL, NULL, vdim, PREF_ALL);
    return currentNrCov - 1;
}

int IncludePrim(const char *name, int kappas, 
		 stationary_type stationary, isotropy_type isotropy,	
		 checkfct check, rangefct range, pref_type pref,
		 int vdim) {  
    createmodel(name, kappas, NULL, stationary, isotropy, check, range,
		NULL, NULL, vdim, pref);
    return currentNrCov - 1;
}
int IncludePrim(const char *name, int kappas, size_fct kappasize,
		 stationary_type stationary, isotropy_type isotropy,	
		 checkfct check, rangefct range,pref_type pref, int vdim) {  
    createmodel(name, kappas, kappasize, stationary, isotropy, check, range,
		NULL, NULL, vdim, pref);
    return currentNrCov - 1;
}


void make_internal() {  
    int nr = currentNrCov - 1;  
    cov_fct *C = CovList + nr;
    C->internal = true;
}


// extern ?!
int IncludeModel(const char *name, char minsub, char maxsub, int kappas,
		  size_fct kappasize,
		  stationary_type stationary, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref, 
		  init_meth initmethod, do_meth domethod, bool internal,
		  int vdim) {  
    createmodel(name, kappas, kappasize, stationary, isotropy, check, range,
		initmethod, domethod, vdim, pref);
//    assert(maxsub > 0); // check deleted 25. nov 2008 due to nugget 
    assert(maxsub >= minsub && maxsub <= MAXSUB);

    int i, 
      nr = currentNrCov - 1;  
    cov_fct *C = CovList + nr;
    C->minsub = minsub;
    C->maxsub = maxsub;  
    C->primitive = false;
    C->internal = internal;

    for (i=0; i<maxsub; i++) {
      sprintf(C->subnames[i], "u%d", i+1); // default (repeated twice)
    }
    return nr;
}
int IncludeModel(const char *name, char minsub, char maxsub, int kappas, 
		  size_fct kappasize,
		  stationary_type stationary, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref) {
  return 
      IncludeModel(name, minsub, maxsub, kappas, kappasize, stationary, isotropy,
		   check, range, pref, NULL, NULL, false, UNIVARIATE);
}
int IncludeModel(const char *name, char minsub, char maxsub, int kappas, 
		  stationary_type stationary, isotropy_type isotropy,
		  checkfct check, rangefct range, pref_type pref) {
  return
      IncludeModel(name, minsub, maxsub, kappas, NULL, stationary, isotropy,
	       check, range, pref, NULL, NULL, false, UNIVARIATE);
}


void addCov(covfct cf, covfct D, covfct D2, natscalefct naturalscale) {
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov) && cf!=NULL);
  cov_fct *C = CovList + nr;
  assert(cf != NULL);
  C->cov = cf;
  if (C->derivs < 0) C->derivs = 0;
  assert(C->nonstat_cov == ErrCovNonstat);
  assert(C->vdim == UNIVARIATE || C->vdim == SUBMODELM || D==NULL);

  if (D != NULL) {
    assert(cf != NULL);
   if (C->cov!=NULL && C->derivs < 1) C->derivs = 1;
    C->D=D;    
    assert(C->isotropy == ISOTROPIC || C->isotropy == SPACEISOTROPIC ||
	   C->isotropy == PREVMODELI);
    C->implemented[TBM2] = NUM_APPROX;
    C->implemented[TBM3] = true; 
  }
  if (D2 != NULL) {
    assert(D != NULL);
    C->D2 = D2;
   if (C->cov!=NULL && C->D != NULL && C->derivs < 2)  C->derivs = 2;
  }
  if (C->maxsub > 0 && currentNrCov>GATTER+10) 
      assert(naturalscale == NULL);
  C->naturalscale=naturalscale;
  C->implemented[Direct] = cf != NULL;
  C->implemented[CircEmbed] = cf != NULL &&
    (C->stationary == STATIONARY || C->stationary == PREVMODELS);
  C->implemented[Sequential] = C->implemented[CircEmbed] && C->vdim <= 1;
}

void addCov(covfct cf, covfct D, natscalefct naturalscale) {
  addCov(cf, D, NULL, naturalscale);
}


void addCov(covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    natscalefct naturalscale) {

    addCov(cf, D, D2, naturalscale);

 
    cov_fct *C = CovList + currentNrCov - 1;
    C->D3 = D3;
    C->D4 = D4;
    assert(C->derivs == 2 && D3!=NULL && D4!=NULL);
    C->derivs = 4;
}



void addCov(nonstat_covfct cf) {
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov) && cf!=NULL);
  cov_fct *C = CovList + nr;

  assert(cf != NULL);
  C->nonstat_cov = cf;
  C->implemented[Direct] = true;
  C->implemented[CircEmbed] = 
      C->stationary == STATIONARY || C->stationary == PREVMODELS;
  C->implemented[Sequential] = C->implemented[CircEmbed] && C->vdim <= 1;

  if (C->derivs < 0) C->derivs = 0;
  
}

void nablahess(covfct nabla, covfct hess) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;

  assert((nr>=0) && (nr<currentNrCov));
  assert(C->vdim == UNIVARIATE || C->vdim == SUBMODELM || nabla==NULL);
  assert(C->cov != NULL && nabla!=NULL && hess != NULL);
  
  C->nabla=nabla;    
  C->hess = hess;

//  if (C->derivs>0) {
//     printf("%s: %d\n", C->name, C->derivs);
//      assert(C->derivs >= 2);
//  }
}

void addCov(aux_covfct auxcf)
{

  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov) && auxcf!=NULL);
  cov_fct *C = CovList + nr;

  assert(C->cov == ErrCov && C->nonstat_cov==ErrCovNonstat && 
	 (C->stationary == AUXMATRIX || C->stationary == ANYFCT || true));
  C->aux_cov = auxcf;
  if (C->derivs < 0) C->derivs = 0;
}

int addFurtherCov(covfct cf, covfct D, covfct D2) {
  assert(currentNrCov > 0);
  cov_fct *C = CovList + currentNrCov;
  memcpy(C, C - 1, sizeof(cov_fct));
  assert(C->vdim == UNIVARIATE || C->vdim == SUBMODELM || D == NULL);

  strcopyN(CovNames[currentNrCov], InternalName, MAXCHAR);
  C->name[0] = InternalName[0];
  strcopyN(C->name + 1, CovList[currentNrCov-1].name, MAXCHAR - 1);
  if (cf != NULL) {
    C->cov = cf;
    C->derivs = 0;
  }
  if (D != NULL) {
    assert(cf != NULL);
    C->D = D;
    C->derivs = 1;

    //   printf("%s %d\n", C->name, C->isotropy);

    assert(C->isotropy == ISOTROPIC || C->isotropy == SPACEISOTROPIC ||
	   C->isotropy == ZEROSPACEISO ||
	   C->isotropy == PREVMODELI);
 
    C->implemented[TBM2] = NUM_APPROX;
    C->implemented[TBM3] = true; 
  }
  if (D2 != NULL) {
    assert(D != NULL);
    C->D2 = D2;
    C->derivs = 2;
  }
  C->internal = true; // addCov is used without previous call of IncludeModel
  currentNrCov++;
  return currentNrCov - 1;
}

int addFurtherCov(covfct cf, covfct D) {
  return addFurtherCov(cf, D, NULL);
}

int addFurtherCov(nonstat_covfct cf) {
  assert(currentNrCov > 0);
  cov_fct *C = CovList + currentNrCov;
  memcpy(C, C - 1, sizeof(cov_fct));
  strcopyN(CovNames[currentNrCov], InternalName, MAXCHAR);
  C->name[0] = InternalName[0];
  strcopyN(C->name + 1, CovList[currentNrCov-1].name, MAXCHAR - 1);
  C->derivs = -1;
  if (cf != NULL) {
    C->nonstat_cov = cf;
    C->derivs = 0;
  }
  C->D = ErrCov;
  C->internal = true; // addCov is used without previous call of IncludeModel
  C->mppget = NULL;
  C->mppgetstat = NULL;
  currentNrCov++;
  return currentNrCov - 1;
}

int addFurtherCov(nonstat_covfct cf, covfct D, covfct D2) {
  assert(currentNrCov > 0);
  cov_fct *C = CovList + currentNrCov;
  memcpy(C, C - 1, sizeof(cov_fct));
  strcopyN(CovNames[currentNrCov], InternalName, MAXCHAR);
  C->name[0] = InternalName[0];
  strcopyN(C->name + 1, CovList[currentNrCov-1].name, MAXCHAR - 1);
  C->derivs = -1;
  if (cf != NULL) {
    C->nonstat_cov = cf;
    C->derivs = 0;
  }
  C->D = ErrCov;
  C->internal = true; // addCov is used without previous call of IncludeModel
  C->mppget = NULL;
  C->mppgetstat = NULL;
  currentNrCov++;
  return currentNrCov - 1;
}


void addLocal(getlocalparam coinit, getlocalparam ieinit) {
  int nr = currentNrCov - 1;
  assert(nr>=0 && nr < currentNrCov) ;
  cov_fct *C = CovList + nr;
  assert(C->D!=ErrCov);
  if ((C->implemented[CircEmbedIntrinsic] = ieinit != NULL)) {
    assert(C->D2 != NULL);
    C->ieinit = ieinit;
  } 
  if ((C->implemented[CircEmbedCutoff] = coinit != NULL)) {
    C->coinit = coinit;
  }
}
void addCallLocal(altlocalparam alt) {
  int nr = currentNrCov - 1;
  assert(nr>=0 && nr < currentNrCov) ;
  cov_fct *C = CovList + nr;
  C->alternative = alt;
}

int addTBM(covfct tbm2) {
  // must be called always AFTER addCov !!!!
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov));
  cov_fct *C = CovList + nr;
  assert(C->stationary == STATIONARY || C->stationary == PREVMODELS);
  assert(C->vdim == UNIVARIATE || C->vdim == SUBMODELM);
  C->tbm2=tbm2;
  if (tbm2 != NULL) {
    // addTBM is called from the other addTBM's -- so tbm2 might
    // be NULL
    assert(C->isotropy==ISOTROPIC || C->isotropy==SPACEISOTROPIC || 
	   C->isotropy==PREVMODELI);
    assert(C->D != ErrCov);
    C->implemented[TBM2] = IMPLEMENTED;
  }
  // IMPLEMENTED must imply the NUM_APPROX to simplify the choice
  // between TBM2 and Tbm2Num
  return nr;
}

void addTBM(covfct tbm2, spectral_init initspectral, spectral_do spectral) {
  int nr = addTBM(tbm2);
  cov_fct *C = CovList + nr;
  assert(C->vdim == UNIVARIATE || C->vdim == SUBMODELM);

  C->initspectral=initspectral;
  C->spectral=spectral;
  C->implemented[SpectralTBM] = true;
  assert(C->stationary == STATIONARY || C->stationary == PREVMODELS);
}

void addTBM(spectral_init initspectral, spectral_do spectral) {
  addTBM((covfct) NULL, initspectral, spectral);
}
	
void addInv(covfct invSq) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;
  assert((nr>=0) && (nr<currentNrCov));
  
  // printf("%s %d  %d\n", C->name, invSq, invstableSq);

  C->invSq=invSq;
}
	
void addHyper(hyper_pp_fct hyper_pp) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;
  assert((nr>=0) && (nr<currentNrCov));
  C->hyperplane=hyper_pp;
  C->implemented[Hyperplane] = hyper_pp!=NULL;
}
		   
//void addSpecialMeth(initstandard initspecial, dometh special)  {
///  int nr = currentNrCov - 1;
//  cov_fct *C = CovList + nr;
//  C->initspecial=initspecial;
//  C->special=special;
//  if ((special!=NULL) || (initspecial!=NULL)) 
//    assert((special!=NULL) && (initspecial!=NULL));
//  C->implemented[Special] = true;
//}

void addMarkov(int *variable) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;
  C->implemented[Markov] = true; 
  *variable = nr;
}

void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_get mppget,
	    sd_fct sd, int loc, bool timesep){
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;

//  print("addmpp %s\n", C->name);

  assert(mppinit!=NULL &&  randomcoin!=NULL && mppget!=NULL);
  C->mppinit = mppinit;
  C->randomcoin = randomcoin;
  C->mppget = mppget;
  assert(C->stationary >= COVARIANCE);
  C->sd = sd;
  C->mpplocations = loc;
  C->implemented[RandomCoin] = true; 
}

void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_get mppget, 
	    sd_fct sd, int loc){
 addMPP(mppinit, randomcoin, mppget, sd, loc, false);
}

void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_get mppget, 
	    sd_fct sd, int loc,
	    ave_fct avef, ave_logfct avelogg) {
  addMPP(mppinit, randomcoin, mppget, sd, loc, true);
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;
  assert(avef != NULL && avelogg != NULL);
   assert(C->stationary == STATIONARY);
  C->avef = avef;
  C->avelogg = avelogg;
  C->implemented[Average] = true; 

}

void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_getstat mppgetstat,
	    sd_fct sd, int loc, bool timesep){
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;
  assert(mppinit!=NULL &&  randomcoin!=NULL && mppgetstat!=NULL);
  C->mppinit = mppinit;
  C->randomcoin = randomcoin;
  C->mppgetstat = mppgetstat;
  // printf("%d\n", C->stationary);
  assert(C->stationary < COVARIANCE || C->stationary == PREVMODELS );
  C->sd = sd;
  C->mpplocations = loc;
  C->implemented[RandomCoin] = true; 
}

void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_getstat mppgetstat, 
	    sd_fct sd, int loc){
 addMPP(mppinit, randomcoin, mppgetstat, sd, loc, false);
}

void addMPP(mpp_init mppinit, mpp_model randomcoin, mpp_getstat mppgetstat, 
	    sd_fct sd, int loc,
	    ave_fct avef, ave_logfct avelogg) {
  addMPP(mppinit, randomcoin, mppgetstat, sd, loc, true);
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;
  assert(avef != NULL && avelogg != NULL);
  assert(C->stationary == STATIONARY || C->stationary == PREVMODELS);
  C->avef = avef;
  C->avelogg = avelogg;
  C->implemented[Average] = true; 

}

void addSpecial(realfct mineigen){
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;
  C->mineigenvalue = mineigen;
}

void addGaussMixture(draw_logmix drawlogmix,
		     log_mixweight logmixweight) {
  int nr = currentNrCov - 1;
  cov_fct *C = CovList + nr;
  assert(drawlogmix != NULL && logmixweight != NULL);
  C->drawlogmix = drawlogmix;
  C->logmixweight = logmixweight;
}

int getmodelnr(char *name) {
  // == -1 if no matching function is found
  // == -2 if multiple matching fctns are found, without one matching exactly
  // == -3 if internal name is passed
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)

  if (currentNrCov==-1) InitModelList();
  if (!strcmp(name, InternalName)) return MATCHESINTERNAL;

  return Match(name, CovNames, currentNrCov);
}

int Match(char *name, name_type List, int n) {
  // == -1 if no matching name is found
  // == -2 if multiple matching fctns are found, without one matching exactly
  unsigned int ln;
  int Nr;
  Nr=0;
  ln=strlen(name);
  while ( Nr < n  && strncmp(name, List[Nr], ln)) Nr++;
  if (Nr < n) { 
    if (ln==strlen(List[Nr])) // exactmatching -- take first -- changed 1/7/07
      return Nr;
    // a matching function is found. Are there other functions that match?
    int j; 
    bool multiplematching=false;
    j=Nr+1; // if two or more covariance functions have the same name 
    //            the last one is taken 
    while (j<n) {
      while ( (j<n) && strncmp(name, List[j], ln)) {j++;}
      if (j<n) {
	if (ln==strlen(List[j])) { // exactmatching -- take first 
	  return j;
	}
	else {multiplematching=true;}
      }
      j++;
    }
    if (multiplematching) {return MULTIPLEMATCHING;}
  } else return NOMATCHING;
  return Nr;
}


void UserGetNatScaling(double *natscale) {
  GetNaturalScaling(STORED_MODEL[MODEL_USER], natscale);
}


void GetNaturalScaling(cov_model *cov, double *natscale)
{ // called also by R 

  // values of naturalscaling:
  //#define NATSCALE_EXACT 1   
  //#define NATSCALE_APPROX 2
  //#define NATSCALE_MLE 3 /* check fitvario when changing !! */
  // +10 if numerical is allowed
   if (cov->nr >= GATTER && cov->nr <= LASTGATTER) cov = cov->sub[0];
  cov_fct *C = CovList + cov->nr;
  *natscale = 0.0;


  // printf("%s %d\n", C->name, C->maxsub);


  if (C->maxsub!=0) XERR(ERRORFAILED); 
 
  if (C->isotropy != ISOTROPIC || C->stationary != STATIONARY ||
      C->vdim != UNIVARIATE) 
      XERR(ERRORANISOTROPIC); 
	 
  
  if (C->naturalscale!=NULL) { 
    *natscale = C->naturalscale(cov);
    if (ISNAN(*natscale) || ISNA(*natscale) || *natscale != 0.0) {
      return;
    }
  }
    
  if (NS != NATSCALE_ORNUMERIC) XERR(ERRORRESCALING); 

  double x,newx,yold,y,newy;
  covfct cf;
  int wave, i;
  if ((C->cov)==nugget)  XERR(ERRORRESCALING); 
  if ((cf=C->cov)==NULL) XERR(ERRORNOTDEFINED);
      
      // already calculated ?
//      parami=KAPPA; // do not compare mean,variance, etc.
//      if (oldcovnr==*covnr) {
//	for (; parami<=LASTKAPPA; parami++) {
//	  if (oldp[parami]!=p[parami]) {
//	    break;
//	  }
//	}
//	if (parami > LASTKAPPA) {
//	  *natscale=OldNatScale; 
//	  return;
//	}
//      }
//      for (; parami<=LASTKAPPA; parami++) {oldp[parami]=p[parami];}

      /* **************************************
	 Now, find a quick and good solution for NewInvScale --
      */
  
  wave  = 0;
  x = 1.0;
  cf(&x, cov, &yold);
  if (ISNA(yold) || ISNAN(yold)) {*natscale = RF_NAN; return;}
  if (yold > 0.05) {
    double leftx;
    x *= 2.0;
    cf(&x, cov, &y);
    while (y > 0.05) {  
      if (yold<y){ 
	  wave++;  
	  if (wave>10) XERR(ERRORWAVING); 
      }
      yold = y;
      x *= 2.0;
      if (x>1E30) XERR(ERRORRESCALING); // or use a separate ERROR
      cf(&x, cov, &y);
    } 
    leftx = x * 0.5;
    for (i=0; i<3 /* good choice?? */ ;i++) {          
	if (y==yold) XERR(ERRORWAVING); // should never appear
      newx = x + (x-leftx)/(y-yold)*(0.05-y);
      cf(&newx, cov, &newy);
      if (newy > 0.05) {
	leftx = newx;
	yold  = newy;
      } else {
	x = newx;
	y = newy;
      }
    }
    if (y==yold) XERR(ERRORWAVING); // should never appear
    *natscale = 1.0 / ( x + (x-leftx)/(y-yold)*(0.05-y) );
  } else {
    double rightx;
    x *= 0.5;
    cf(&x, cov, &y);
    while (y < 0.05) {  
	if (yold>y){ wave++;  if (wave>10) XERR(ERRORWAVING); }
      yold = y;
      x *= 0.5;
      if (x<1E-30) XERR(ERRORRESCALING); //or use a separate ERROR
      cf(&x, cov, &y);
    }    
    rightx = x * 2.0;
    for (i=0; i<3 /* good choice?? */ ;i++) {          
	if (y==yold) XERR(ERRORWAVING); // should never appear
      newx = x + (x-rightx)/(y-yold)*(0.05-y);
      cf(&x, cov, &newy);
      if (newy < 0.05) {
	rightx = newx;
	yold   = newy;
      } else {
	x = newx;
	y = newy;
      }
    }
    if (y==yold) XERR(ERRORWAVING); // should never appear
    *natscale = 1.0 / ( x + (x-rightx)/(y-yold)*(0.05-y));
  }
}





void Getxsimugr(coord_type x, double *aniso, int timespacedim, double **xsimugr) {
  // bisher nur fuer Diagonalmatrizen 
  int n, i, w;
  if (aniso == NULL) {
    for(w=0; w<timespacedim; w++) {
      for (i=0; i<3; i++) {
	  xsimugr[w][i] = x[w][i];
      }
    } 
  } else {  
      for(n=w=0; w<timespacedim; w++, n+=timespacedim+1) {
	  for (i=0; i<3; i++) {
	      xsimugr[w][i] = aniso[n] * x[w][i];
	  }
      }
  }
}

void GetCoordinates(int *keyNr, double *x, int *err)
{
  key_type *key;
  if (currentNrCov==-1) {*err=-1;goto ErrorHandling;}
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*err=-2; goto ErrorHandling;}
  key = &(KEY[*keyNr]);
  if (!key->simu.active) {*err=-3; goto ErrorHandling;}
  if (key->trend.TrendModus==0) {*err=-4; goto ErrorHandling;}
 
  *err = ERRORNOTPROGRAMMED;  // TODO !!

  return;

 ErrorHandling:
  return;
}


void covcpy(cov_model **localcov, cov_model *cov, bool insertgatter,
	    bool keepuser) {
  assert(cov != NULL);
  int i,
    covnr = cov->nr,
    n = -1;
  cov_model *current;
  cov_fct *C = CovList + cov->nr;
  if (insertgatter && covnr >= GATTER && covnr <= LASTGATTER)
    error("covcpy detects # at forbidden place -- please contact author");
    
  if (*localcov != NULL) error("local cov not NULL");
  if ((*localcov = (cov_model*) malloc(sizeof(cov_model)))==0) 
    ERR("memory allocation error");
  if (insertgatter) {  // covnr != DOLLAR && 
    COV_NULL(*localcov);
    (*localcov)->nr = GATTER;
    (*localcov)->nsub = 1;
    (*localcov)->isoIn = cov->isoIn;
    (*localcov)->statIn = cov->statIn;
    (*localcov)->vdim = cov->vdim;
    if (((*localcov)->sub[0] = (cov_model*) malloc(sizeof(cov_model)))==0)
      ERR("memory allocation error");
   current = (*localcov)->sub[0];
  } else current = *localcov;
  
  memcpy(current, cov, sizeof(cov_model)); // replaces COV_NULL(*localcov);
  if (!keepuser) {
    for (i=0; i<Forbidden; i++)
      current->user[i] = PREF_BEST;
  }

  current->calling = *localcov; 
  
  for (i=0; i<MAXPARAM; i++) {
    if (cov->p[i] == NULL) {
      continue;
    }
    if (C->kappatype[i] >= LISTOF) {
      current->p[i] = (double*) malloc(sizeof(listoftype));
      int j, len = cov->nrow[i];
      listoftype *p = (listoftype*) (cov->p[i]),
	  *q = (listoftype*) (current->p[i]);	
      for (j=0; j<len; j++) {
	if (C->kappatype[i]== LISTOF + REALSXP) {
	  n = sizeof(double);
	} else assert(false);
	n *= p->nrow[j] * p->ncol[j];
	q->nrow[j] = p->nrow[j];
	q->ncol[j] = p->ncol[j];
	q->p[j] = (double*) malloc(n);
	memcpy(q->p[i], p->p[i], n);	    
      }
    } else {
      if (C->kappatype[i]==REALSXP) {
        n = sizeof(double);
      } else if (C->kappatype[i]==INTSXP) {
        n = sizeof(int);
      } else assert(false);
      n *= cov->nrow[i] * cov->ncol[i];
      current->p[i] = (double*) malloc(n);
      memcpy(current->p[i], cov->p[i], n);
    }
  }
  if (cov->q != NULL) {
    n = sizeof(double) * current->qlen;
    current->q = (double*) malloc(n);
    memcpy(current->q, cov->q, n);
  } else assert(current->qlen==0);
  
  for (i=0; i<MAXSUB; i++) {
    current -> sub[i] = NULL;
    if (cov->sub[i] == NULL) continue;
    if (insertgatter && 
	cov->sub[i]->nr >= GATTER && cov->sub[i]->nr <= LASTGATTER) 
      covcpy(current->sub + i, cov->sub[i]->sub[0], insertgatter, keepuser);
    else 
      covcpy(current->sub + i, cov->sub[i], insertgatter, keepuser);
    current->sub[i]->calling = current; 
  }
}



//void covcpy(cov_model **localcov, cov_model *cov) {
//  covcpy(localcov, cov, true, true);
//}

double *getAnisoMatrix(method_type *meth) {
  double *ani;
  int i, bytes, 
      origdim = meth->loc->timespacedim,
      dimP1 = origdim + 1;
  if (meth->caniso != NULL) {
    bytes = sizeof(double) * origdim * meth->xdimout;
    ani = (double *) malloc(bytes);

//      PrintMethodInfo(meth);
//      printf("%ld  %ld %d %d %d\n", ani, meth->caniso, bytes,
//	     origdim,  meth->xdimout); assert(false);

    memcpy(ani, meth->caniso, bytes);
  } else {
    double a = 1.0 / meth->cscale;  
    if (meth->cproj == NULL) {
      bytes = origdim * origdim;
      ani = (double *) calloc(bytes,  sizeof(double));
      for (i=0; i<bytes; i+=dimP1) ani[i] = a; 
    } else {
      bytes = origdim * meth->xdimout;
      ani = (double *) calloc(bytes,  sizeof(double));
      for (i=0; i<meth->xdimout; i++) ani[i + meth->cproj[i] * origdim] = a; 	
    }
  }
  return ani;
}
 
double *getAnisoMatrix(cov_model *cov, int dim) {
  double *ani;
  int i, bytes, dimout = dim,
      dimP1 = dim + 1,
       *proj = (int *)cov->p[DPROJ];
  if (cov->p[DANISO] != NULL) {
      bytes = sizeof(double) * dim * dimout;
      ani = (double *) malloc(bytes);
      memcpy(ani, cov->p[DANISO], bytes);
  } else {
    double a;
    a = (cov->p[DSCALE] == NULL) ? 1.0 : 1.0 / cov->p[DSCALE][0];  
    if (proj == NULL) {
      bytes = dim * dimout;
      ani = (double *) calloc(bytes,  sizeof(double));
      for (i=0; i<bytes; i+=dimP1) ani[i] = a; 
    } else {
      dimout = cov->nrow[DPROJ];
      bytes = dim * dimout;
      ani = (double *) calloc(bytes,  sizeof(double));
      for (i=0; i<dimout; i++) ani[i + proj[i] * dim] = a; 	
    }
  }
  return ani;
}

double GetDiameter(double *origcenter, double *origmin, double *origmax, 
		   double *aniso, int origdim, int spatialdim,
		   double *min, double *max, double *center) {
  
  // calculates twice the distance between origcenter %*% aniso and
  // all 2^spatialdim corners spanned by min and max; returns maximum distance.

  // NOTE: origcenter is not alway (min + max)/2 since origcenter might be
  //       given by the user and not calculated automatically !

  double *lx, *sx, distsq, dummy,
    diameter=0.0;
  bool *j;
  int d;
 

  j = (bool*) malloc( (origdim + 1) * sizeof(double));
  lx = (double*) malloc( origdim * sizeof(double));
  sx = (double*) malloc(spatialdim * sizeof(double));

  Aniso(origcenter, aniso, origdim, spatialdim, center);
  for (d=0; d<origdim; d++) {
    j[d]=false;
    lx[d]=origmin[d];
  }
  j[origdim] = false;
  
  for (d=0; d<spatialdim; d++) {
      min[d] = R_PosInf;
      max[d] = R_NegInf;
  }

  while(true) {
    d=0; 
    while(j[d]) {
      lx[d]=origmin[d];
      j[d++]=false;
    }
    if (d==origdim) break;
    j[d]=true;
    lx[d]=origmax[d];
    Aniso(lx, aniso, origdim, spatialdim, sx);
    
    // suche maximale Distanz von center zu transformierter Ecke
    distsq = 0.0;
    for (d=0; d<spatialdim; d++) {
      if (min[d] > sx[d]) min[d] = sx[d];
      if (max[d] < sx[d]) max[d] = sx[d];
      dummy = center[d] - sx[d];
      distsq += dummy * dummy;
      
      //  printf("%d sx=%f %f %f diam=%f\n", d, sx[d], min[d], max[d], distsq);

    }
    if (distsq > diameter) diameter = distsq;
  }

  free(j);
  free(lx);
  free(sx);
  return 2.0 * sqrt(diameter);
}

double GetDiameter(double *origcenter, double *min, double *max, 
		   double *aniso, int tsdim, int spatialdim) { 
    double diam,
	*dummymin, *dummymax, *dummycenter;
  
  dummymin = (double*) malloc(spatialdim * sizeof(double));
  dummymax = (double*) malloc(spatialdim * sizeof(double));
  dummycenter = (double*) malloc(spatialdim * sizeof(double));
  diam = GetDiameter(origcenter, min, max, 
		     aniso, tsdim, // meth->loc->timespacedim,
		     spatialdim,
		     dummymin, dummymax, dummycenter);
  free(dummymin);
  free(dummymax);
  free(dummycenter);
  return diam;
}


void GetMinMax(method_type *meth, double *min, double *max, double *center,
    int MaxDim) {
  // untersucht die Original-Werte fuer Koordinaten fuer center!!
  // aber gibt diameter in den neuen Koordinaten wieder.

  location_type *loc = meth->loc;
  int d,
    origdim = loc->timespacedim,  
    spatialdim = loc->spatialdim;  

  if (loc->grid) {
     // diameter of the grid
      
      //   printf("origdim =%d\n", origdim);

    for (d=0; d<origdim; d++) {   
      if (loc->xgr[d][XEND] > loc->xgr[d][XSTART]) {
	min[d] = loc->xgr[d][XSTART];
	max[d] = loc->xgr[d][XSTART] + loc->xgr[d][XSTEP] * 
	  (double) (loc->length[d] - 1);
      } else {
	min[d] = loc->xgr[d][XSTART] + loc->xgr[d][XSTEP] * 
	  (double) (loc->length[d] - 1);
	max[d] = loc->xgr[d][XSTART];
      }
    }

    // assert(!false);

  } else { // not simugrid
    double *xx=loc->x; 
    int i,
      endfor = loc->length[0] * loc->timespacedim;

    for (d=0; d<MaxDim; d++) {min[d]=RF_INF; max[d]=RF_NEGINF;}
    // to determine the diameter of the grid determine as approxmation
    // componentwise min and max corner
    for (i=0; i<endfor; ) {
      for (d=0; d<spatialdim; d++, i++) {
        //temporal part need not be considered, but for ease included
	if (xx[i]<min[d]) min[d] = xx[i];
	if (xx[i]>max[d]) max[d] = xx[i];
      }
    }
        
    if (loc->Time) {
      assert(d == origdim - 1);
      if (loc->T[XEND] > loc->T[XSTART]) {
	min[d] = loc->T[XSTART];
	max[d] = loc->T[XSTART] + loc->T[XSTEP] * 
	  (double) (loc->length[d] - 1);
      } else {
	min[d] = loc->T[XSTART] + loc->T[XSTEP] * (double) (loc->length[d] - 1);
	max[d] = loc->T[XSTART];
      }
    }
  } // !simugrid


  for (d=0; d<origdim; d++) { // note here no time component!
      center[d] = 0.5 * (max[d] + min[d]); 
      //    printf("%d center %f\n", d, center[d]);
  }
  for ( ; d<MaxDim; d++) {
      min[d] = 0.0;
      max[d] = 0.0;
      center[d] = 0.0;
  }

  // assert(false);
}


bool ok_n(int n, int *f, int nf) // taken from fourier.c of R
{
  int i;
  for (i = 0; i < nf; i++)
    while(n % f[i] == 0) if ((n /= f[i]) == 1) return TRUE;
  return n == 1;
}
int nextn(int n, int *f, int nf) // taken from fourier.c of R
{
    while(!ok_n(n, f, nf)) { n++; }
  return n;
}
#define F_NUMBERS1 3
#define F_NUMBERS2 3
bool HOMEMADE_NICEFFT=false;
unsigned long NiceFFTNumber(unsigned long n) {
  int f[F_NUMBERS1]={2,3,5}; 
  unsigned long i,ii,j,jj,l,ll,min, m=1;
  if (HOMEMADE_NICEFFT) {
    if (n<=1) return n;
    for (i=0; i<F_NUMBERS1; i++) 
      while ( ((n % f[i])==0) && (n>10000)) { m*=f[i]; n/=f[i]; }
    if (n>10000) {
      while (n>10000) {m*=10; n/=10;}
      n++;
    }
    min = 10000000;
    for (i=0, ii=1; i<=14; i++, ii<<=1) { // 2^14>10.000
      if (ii>=n) { if (ii<min) min=ii; break; }
      for (j=0, jj=ii; j<=9; j++, jj*=3) {// 3^9>10.000
	if (jj>=n) { if (jj<min) min=jj; break; }
	
	//for (k=0, kk=jj; k<=6; k++, kk*=5) {// 5^6>10.000
	//if (kk>=n) { if (kk<min) min=kk; break; }
	//for (l=0, ll=kk; l<=5; l++, ll*=7) {// 7^5>10.000
	
	// instead of (if 7 is included)
	for (l=0, ll=jj; l<=6; l++, ll*=5) {// 5^5>10.000
	  	  
	  if (ll>=n) { 
	    if (ll<min) min=ll;
	    break;
	  }
	  //}
	}
      }
    }
    return m*min;
  } else { // not HOMEMADE_NICEFFT
    int f[F_NUMBERS2]={2,3,5}; 
    return nextn(n, f, F_NUMBERS2);
  }
}

void cpyMethod(method_type *meth, method_type *ziel, bool cp_aniso) {
  int xdimout = meth->xdimout,
      tsdim = meth->loc->timespacedim;
  
  ziel->gp   = meth->gp;
  ziel->gpdo = meth->gpdo;
  ziel->simu = meth->simu;
  ziel->loc  = meth->loc;
  // nr  (by METHOD_NULL)
  // compatible  (by METHOD_NULL)
  // nsub  (by METHOD_NULL)
  // domethod  (by METHOD_NULL)
  // destruct  (by METHOD_NULL)
  // S  (by METHOD_NULL)
  // cov  (by METHOD_NULL)
  if (cp_aniso) {
    if (meth->cproj != NULL) {
      ziel->cproj = (int*) malloc(xdimout * sizeof(int));
      assert(ziel->cproj != NULL);
      memcpy(ziel->cproj, meth->cproj, xdimout * sizeof(int));
    }
    if (meth->caniso != NULL) {
	ziel->caniso = (double*) malloc(xdimout * tsdim * sizeof(double));
	assert(ziel->caniso != NULL);
	memcpy(ziel->caniso, meth->caniso, xdimout * tsdim * sizeof(double));
    }
  }
  ziel->cscale = meth->cscale;
  ziel->cvar = meth->cvar;
  // cproj  (by METHOD_NULL)
  ziel->xdimout = meth->xdimout;
  ziel->type = meth->type;
  // hanging  (by METHOD_NULL)
  // space  (by METHOD_NULL)
  // sptime  (by METHOD_NULL)
  // grani  (not set !)
} 

void expandgrid(coord_type xgr, int *len, double **xx, int nrow){
  double *x,* y; /* current point within grid, but without
		       anisotropy transformation */
  int pts, w, k, total, i, d, dimM1 = nrow - 1,
      *yi, /* counter for the current position in the grid */
      ncol = nrow;
  for (pts=1, i=0; i<nrow; i++) {
    pts *= len[i];
  }

  y = (double*) malloc(ncol * sizeof(double));
  yi = (int*) malloc(ncol * sizeof(int));

  total = ncol * pts;

  x = *xx = (double*) malloc(sizeof(double) * total);

  for (w=0; w<ncol; w++) {y[w]=xgr[w][XSTART]; yi[w]=0;}
  for (k=0; k<total; ){
    for (d=0; d<ncol; d++, k++) {
      x[k] = y[d];
    }
    i = 0;
    (yi[i])++;
    y[i] += xgr[i][XSTEP];
    while(yi[i]>=len[i]) {
      yi[i]=0;
      y[i] = xgr[i][XSTART];
      if (i<dimM1) {
	i++;
	(yi[i])++;
	y[i] += xgr[i][XSTEP];
      } else {
	assert(k==total);
      }
    }
  }
  free(y);
  free(yi);
}


void expandgrid(coord_type xgr, int *len, double **xx, double* aniso, 
		int nrow, int ncol){
  double *x, * y; /* current point within grid, but without
		       anisotropy transformation */
  int pts, w, k, total, n, i, d, 
      *yi,   /* counter for the current position in the grid */
      dimM1 = nrow - 1; 
  for (pts=1, i=0; i<nrow; i++) {
    pts *= len[i];
  }

  total = ncol * pts;
  x = *xx = (double*) malloc(sizeof(double) * total);
  y = (double*) malloc(ncol * sizeof(double));
  yi = (int*) malloc(ncol * sizeof(int));


  for (w=0; w<nrow; w++) {y[w]=xgr[w][XSTART]; yi[w]=0;}
  for (k=0; k<total; ){
    if (aniso==NULL) {
      for (d=0; d<ncol; d++) x[d] = y[d];
    } else {
      for (n=d=0; d<ncol; d++, k++) {
        x[k] = 0.0;
	for(w=0; w<nrow; w++) {
	  x[k] += aniso[n++] * y[w];
	}
      }
    }
    i = 0;
    (yi[i])++;
    y[i] += xgr[i][XSTEP];
    while(yi[i]>=len[i]) {
      yi[i]=0;
      y[i] = xgr[i][XSTART];
      if (i<dimM1) {
	i++;
	(yi[i])++;
	y[i] += xgr[i][XSTEP];
      } else {
	assert(k==total);
      }
    }
  }
  free(y);
  free(yi);
}

void grid2grid(coord_type xgr, static_grid grani, double *aniso, int dim) {
  // diag is assumed
  double factor;
  int n, d, i;
  if (aniso == NULL) {
    for(d=0; d<dim; d++) {
      for (i=0; i<3; i++) {
        grani[d][i] = xgr[d][i];
      }
    }
  } else {
    for(n=d=0; d<dim; d++, n+=dim+1) {
      factor = aniso[n];
      // factor = d < cov->nrow ? aniso[n] : 0.0; // hier zwingend diag !!
      for (i=0; i<3; i++) {
        grani[d][i] = xgr[d][i] * factor;
      }
    }
  }
}

void xtime2x(double *x, int nx, double *T, int timelen, 
	     double **newx, int timespacedim) {
  double *y, // umhaengen zwingend notwendig, da u.U. **xx und *newx
    // der gleiche aufrufende Pointer ist, d.h. es wird "ueberschrieben"
    *z,  t;
  int j, k, i, d, 
    spatialdim = timespacedim - 1;
  z = *newx = (double*) malloc(sizeof(double) * timespacedim * nx * timelen);
  for (k=j=0, t=T[XSTART]; j<timelen; j++, t += T[XSTEP]){
   y = x;
   for (i=0; i<nx; i++) {
       for (d=0; d<spatialdim; d++, y++) {
	z[k++] = *y; 
      }
      z[k++] = t;
    }
  }
}

void xtime2x(double *x, int nx, double *T, int len, double **newx,  
	     double *aniso, int nrow, int ncol) {
  double *y, // umhaengen zwingend notwendig, da u.U. **xx und *newx
    // der gleiche aufrufende Pointer ist, d.h. es wird "ueberschrieben"
    *z, dummy, t;
  int j, k, i, d, n, endfor, w,
      spatialdim = nrow - 1,
      nxspdim = nx * spatialdim,
      ncolM1 = ncol - 1;
  z = *newx = (double*) malloc(sizeof(double) * ncol * nx); 
  if (aniso == NULL) {
    assert(ncol == nrow);
    for (k=j=0, t=T[XSTART]; j<len; j++, t += T[XSTEP]) {
      y = x;
      for (i=0; i<nx; i++) {
	  for (d=0; d<ncolM1; d++, y++) {
	  z[k++] = *y; //auf keinen Fall nach vorne holen
	}
	z[k++] = t;
      }
    }
  } else {
    for (k=j=0, t=T[XSTART]; j<len; j++, t += T[XSTEP]){
      y = x;
      for (i=0; i<nxspdim; i+=spatialdim) {
        endfor = i + spatialdim;
        for (n=d=0; d<ncol; d++) {
	  dummy = 0.0;
          for(w=i; w<endfor; w++) {
	    dummy += aniso[n++] * y[w];
	  }
//  printf("%d %d %f\n", k, n, t);
	  z[k++] = dummy + aniso[n++] * t; //auf keinen Fall nach vorne holen
	}
      }
    }
  }
}

void x2x(double *x, int nx, double **newx, 
	 double *aniso, int nrow, int ncol) {
  double *y = x, // umhaengen zwingend notwendig, da u.U. **xx und *newx
    // der gleiche aufrufende Pointer ist, d.h. es wird "ueberschrieben"
    *z, dummy;
  int k, i, d, n, endfor, w,
      nxnrow = nx * nrow;

  z = *newx = (double*) malloc(sizeof(double) * ncol * nx);  
  if (aniso == NULL) {
    assert(nrow == ncol);
    memcpy(z, y, sizeof(double) * nrow * ncol);
  } else {
    for (k=i=0; i<nxnrow; i+=nrow) {
      endfor = i + nrow;
      for (n=d=0; d<ncol; d++) {
        dummy = 0.0;
	for(w=i; w<endfor; w++) {
	  dummy += aniso[n++] * y[w];
	}
	z[k++] = dummy; 
      }
    }
  }
}


double *matrixmult(double *m1, double *m2, int dim1, int dim2, int dim3) {
    double dummy, 
	*m0 = (double*) malloc(sizeof(double) * dim1 * dim3);
  int i,j,k;
  for (i=0; i<dim1; i++) {
    for (k=0; k<dim3; k++) {
      dummy = 0.0;
      for (j=0; j<dim2; j++) {
	dummy += m1[i + j * dim1] * m2[j + k * dim2];
      }
      m0[i + dim1 * k] = dummy;
    }
  }
  return m0;
}

matrix_type Type(double *m, int nrow, int ncol) {

    // see also analyse_matrix: can it be unified ???

  matrix_type type = TypeIso; // default
  double elmt;
  int i, endfor, endfor2;

  if (m==NULL) {
      assert(ncol == nrow);
      return TypeIso;
  }
  
  if (ncol > nrow) {
    endfor =  nrow * ncol;
    for (i=nrow * ncol; i < endfor; i++) 
      if (m[i]!= 0.0) return TypeAny;
    ncol = nrow;
  }

  i = nrow * ncol -1;
  if (i==0) return TypeIso;
  if (ncol==1) {
    endfor = nrow - 1;
    for (i=1; i<endfor; i++) if (m[i] != 0.0) return TypeAny;
    if (m[endfor] == 0.0) return TypeIso;
    if (m[0] == 0.0) return TypeIso;
    return TypeAny;
  }

  elmt = m[i];
  endfor = i - nrow;

//  printf("endfor %d, i=%d\n", endfor, i);

  for (i--; i>endfor; i--) { // last row but very last element: not all zero ? 
    if (m[i] != 0.0) return TypeAny;
  }
  endfor--;

  for (i = 0; i < endfor; ) { 
    if (m[i++] != elmt) type=TypeDiag; // not iso
    endfor2 = i + nrow;
    for ( ; i<endfor2; )
      if (m[i++] != 0.0) return TypeTimesep; // time compent of x is not 
    //   mixed with spatial components of x; (other way round may be possible
  }
  if (m[i++] != elmt) type=TypeDiag; // vorletzte Spalte, vorletzte Zeile,
  //                        d.h. vorletztes Diagonalelement
  if(m[i++] != 0.0) return TypeTimesep; // vorletzte Spalte, letzte Zeile
  
//  printf("%d %d %d %d\n", i , nrow, (ncol-1), type);

  assert(i == nrow * (ncol-1));
  return (type);
}



void Transform2NoGridExt(method_type *meth, bool timesep, bool gridexpand, 
		      double **Space, double **Sptime) {
    // this function transforms the coordinates according to the anisotropy
  // matrix 
  
  location_type *loc = meth->loc;
  int origdim = loc->timespacedim,
    dim = meth->cov->tsdim,
    *length = loc->length;
  double *aniso=NULL, 
    *T = loc->T,
    *x = loc->x;

  aniso = getAnisoMatrix(meth);
  if (loc->grid && meth->type <= TypeDiag) {
    if (gridexpand) {
      if (*Sptime == NULL) {
	expandgrid(loc->xgr, length, Sptime, aniso, origdim, dim);
     }
      if (!loc->Time) *Space = *Sptime;
    } else {
      grid2grid(loc->xgr, meth->grani, aniso, origdim);
    }
  } else if (loc->grid) { // aniso not diag
    if (!loc->Time) { // grid no time
      assert(*Sptime == *Space);
      if (*Space == NULL) {
	expandgrid(loc->xgr, length, Space, aniso, origdim, dim);
	*Sptime = *Space;
      }
    } else { // grid and time 
      if (timesep) {
	if (*Space == NULL) 
	  expandgrid(loc->xgr, length, Space, aniso, origdim, dim - 1);
      } else {
	if (*Sptime == NULL) 
	  expandgrid(loc->xgr, length, Sptime, aniso, origdim, dim);
      }
    }
  } else { // nogrid
    if (!loc->Time) { // no grid no time
       assert(*Sptime == *Space);
       if (*Space == NULL) {
	 x2x(x, length[0], Space, aniso, origdim, dim); 
	 *Sptime = *Space;
       }
    } else { // no grid but time
      if (timesep) {
	 x2x(x, length[0], Space, aniso, origdim, dim-1);
      } else {
	xtime2x(x, length[0], T, length[dim-1],Sptime, aniso, origdim, dim);
      }
    }
  }
  if (aniso != NULL) free(aniso);
}

void Transform2NoGrid(method_type *meth, bool timesep) {
  Transform2NoGridExt(meth, timesep, false, &(meth->space), &(meth->sptime));
}

int InternalGetGridSize(double *x[MAXSIMUDIM], int dim, int *lx)
{  
  int d;
  for (d=0; d<dim; d++) {
    if ((x[d][XEND]<=x[d][XSTART]) || (x[d][XSTEP]<=0.0)) {
      return ERRORCOORDINATES;
    }
    lx[d] = 1 + (int)((x[d][XEND]-x[d][XSTART])/x[d][XSTEP]); 
  }
  return NOERROR;
}


double *selectlines(double *m, int *sel, int nsel, int nrow, int ncol) {
    // selects lines in matrix m according to sel
  int j;  
  double *red_matrix,
      *Red_Matrix = (double *) malloc(sizeof(double) * nsel * ncol),
      *endfor = Red_Matrix + nsel * ncol;
  
  for (red_matrix=Red_Matrix; red_matrix<endfor; m+=nrow) {
      for (j=0; j<nsel; j++, red_matrix++) {
	  *red_matrix = m[sel[j]];
      }
  }
  return Red_Matrix;
} 


int *selectlines(int *m, int *sel, int nsel, int nrow, int ncol) {
  int j;  
  int *red_matrix,
      *Red_Matrix = (int *) malloc(sizeof(int) * nsel * ncol),
      *endfor = Red_Matrix + nsel * ncol;
  
  for (red_matrix=Red_Matrix; red_matrix<endfor; m+=nrow) {
      for (j=0; j<nsel; j++, red_matrix++) {
	  *red_matrix = m[sel[j]];
      }
  }
  return Red_Matrix;
}



void cum_dollar(method_type *meth, int origdim, cov_model *cov, 
		   method_type *ziel) {
  /* 
       cpy meth to ziel with dollar...
       schiebe scale/proj/aniso in meth->cscale/cproj/caniso
       Sowie: setze meth->xdimout
       
  */
    
    int i, total,
      nproj = cov->nrow[DPROJ],
      *proj = (int *) cov->p[DPROJ];
  assert(meth != ziel); 

  if (nproj > 0) {  // DSCALE is given at the same time !
    if (meth->caniso == NULL) {
	ziel->cscale = meth->cscale * cov->p[DSCALE][0];
      if (meth->cproj == NULL) {
	ziel->cproj = (int*) malloc(sizeof(int) * nproj);
	memcpy(ziel->cproj, proj, nproj * sizeof(int));
      } else {
        ziel->cproj = selectlines(meth->cproj, proj, nproj, -9999, 1);
      }
      ziel->type = TypeDiag;
    } else { //cov->p[DPROJ] != NULL, mcaniso != NULL
      ziel->caniso = selectlines(meth->caniso, proj, 
				 nproj, origdim, meth->xdimout);
    }
    ziel->xdimout = nproj;
  } else {
    if (cov->p[DANISO] == NULL) {
      total = origdim * meth->xdimout;
      ziel->xdimout = meth->xdimout;
      if (meth->caniso != NULL) {  
        double invscale = 1.0 / cov->p[DSCALE][0];
	ziel->caniso = (double *) malloc(total * sizeof(double));
	for (i=0; i<total; i++) ziel->caniso[i] = meth->caniso[i] * invscale;
      } else {
	ziel->cscale = meth->cscale * cov->p[DSCALE][0]; 
      }
    } else { // cov->p[DANISO] != NULL
      int ncol = cov->ncol[DANISO];
      total = origdim * cov->ncol[DANISO];
      if (meth->caniso == NULL) {
        if (meth->cproj == NULL) { 
	  ziel->caniso = (double *) malloc(total * sizeof(double));
	  //  printf("%d %d %d\n", total, ncol, cov->nrow[DANISO]);
	  // assert(false);
	  memcpy(ziel->caniso, cov->p[DANISO], total * sizeof(double));
	} else {
	  double *ca, *endfor,
	      *aniso = cov->p[DANISO];
	  ziel->caniso =
	      (double *) calloc(origdim * ncol, sizeof(double)); 
	  //                       timespacedim x ncol
	  endfor = ziel->caniso + ncol * origdim;
	  for (ca = ziel->caniso; ca<endfor;  ca += origdim) {
	    for (i=0; i<nproj; i++, aniso++) {
	      ca[meth->cproj[i]] = *aniso; 
	    }
	  }
	  free(meth->cproj);
	  meth->cproj = NULL;
	}
      } else {
	  assert(meth->xdimout == cov->nrow[DANISO]);
	  ziel->caniso = matrixmult(meth->caniso, cov->p[DANISO], origdim, 
				    meth->xdimout, ncol);
	  // caniso = caniso %*%cov->p
      }
      if (origdim != ncol) ziel->type = TypeAny;
      else ziel->type = Type(ziel->caniso, origdim, ncol);
      ziel->xdimout = ncol;
    }
  }
}
