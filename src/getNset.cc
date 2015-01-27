
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 (library for simulation of random fields)
 Copyright (C) 2001 -- 2014 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

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
#include <R_ext/Linpack.h>

#include "RF.h"
#include "primitive.h"
#include "Coordinate_systems.h"


void LOC_NULL(location_type *loc) {
  int d;
  loc->spatialdim = loc->timespacedim = loc->lx = loc->xdimOZ = -1;
  for (d=0; d<MAXSIMUDIM; d++) {
    loc->xgr[d] = loc->ygr[d] = NULL;
  }
  loc->totalpoints = loc->spatialtotalpoints = 0;
  loc->grid = loc->distances = loc->Time = false;
  loc->delete_x = true;
  loc->x = loc->y = loc->caniso = NULL;
  loc->T[0] = loc->T[1] = loc->T[2] = 0.0;
  loc->i_row = loc->i_col = I_COL_NA;
  loc->cani_ncol = loc->cani_nrow = NA_INTEGER;
}
 
     
void LOC_DELETE(location_type **Loc) {
  // print("LKOC DELETE!");
  location_type *loc = *Loc;
  if (loc == NULL) return;
    
  if (loc->x != NULL) {
    if (loc->delete_x) {
      if (loc->y != NULL && loc->y != loc->x) UNCONDFREE(loc->y);
      UNCONDFREE(loc->x); 
    }     
  }
  FREE(loc->caniso);
  // it may happen that both are set ! Especially after calling 
  // partial_loc_set in variogramAndCo.cc
  if (loc->xgr[0] != NULL && loc->spatialdim>0) {
    if (loc->ygr[0] != NULL && loc->ygr[0] != loc->xgr[0]) 
      UNCONDFREE(loc->ygr[0]);
    UNCONDFREE(loc->xgr[0]); ///
  }
  
  UNCONDFREE(*Loc);
}


//void GLOBAL_NULL() {
//  ReSetGlobal();
//}
//
//void SetTemporarilyGlobal(globalparam *gp) {
//  if (!GLOBALORIG.set) {
//    MEMCOPY(&(GLOBALORIG.param), &GLOBAL, sizeof(globalparam));
//    GLOBALORIG.set = true;
//  }
//  MEMCOPY(&GLOBAL, gp, sizeof(globalparam));
//}
//
//void ReSetGlobal() {
//  if (GLOBALORIG.set) {
//    MEMCOPY(&GLOBAL, &(GLOBALORIG.param), sizeof(globalparam));
//    GLOBALORIG.set = false;
//  }
//}

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

 void MPPPROPERTIES_NULL(mpp_properties *Mpp) {
   // Mpp->refradius= 
   int i;
   for (i=0; i<MAXMPPVDIM; i++) Mpp->maxheights[i] = RF_INF;
   Mpp->unnormedmass = RF_NA;
   Mpp->mM = Mpp->mMplus = NULL;
   //   Mpp->refsd = 
   //Mpp->totalmass = RF_NA;
   //Mpp->methnr = -1;
   //Mpp->loc_done = false;
 }

 void MPPPROPERTIES_DELETE(mpp_properties *Mpp) {   
   FREE(Mpp->mM);
   FREE(Mpp->mMplus);
 }


void COV_DELETE_WITHOUTSUB(cov_model **Cov) {
  cov_model *cov = *Cov;

  //printf("deleting %s\n", NAME(cov));
  assert(cov != NULL);

  int i, j,
    last = (cov->nr < 0) ? MAXPARAM : CovList[cov->nr].kappas; 
  
  for (i=0; i<last; i++) {
    if (!PisNULL(i)) {
      if (CovList[cov->nr].kappatype[i] == LANGSXP) {
	sexp_type *S = PSEXP(i);
	if (S->Delete) R_ReleaseObject(S->sexp);	
      } else if (CovList[cov->nr].kappatype[i] >= LISTOF) {
	listoftype *list = PLIST(i);
	if (list->deletelist) {
	  for (j=0; j<cov->nrow[i]; j++) {
	    UNCONDFREE(list->p[j]);  
	  }
	}
      }
      //printf("    - %d %d\n", i, last);
      PFREE(i);
      cov->ncol[i] = cov->nrow[i] = SIZE_NOT_DETERMINED; // ==0
    }
  }

  MPPPROPERTIES_DELETE(&(cov->mpp));

  if (cov->ownkappanames != NULL) {
    int kappas = CovList[cov->nr].kappas;
    for (j=0; j<kappas; j++) FREE(cov->ownkappanames[j]);
    UNCONDFREE(cov->ownkappanames);
  }
  
  QFREE;
 
  // important check in combination with above; can be easily removed or 
  // generalised  !!!!
 
  FREE(cov->MLE);

  //MPPPROPERTIES_DELETE(&(cov->mpp));

  cov->prevloc = NULL;
  LOC_DELETE(&(cov->ownloc));
  if (cov->key != NULL) {
    // printf("deleting key %s\n", CovList[cov->key->nr].name);
    COV_DELETE(&(cov->key));
  }
  if (cov->rf != NULL && cov->origrf) UNCONDFREE(cov->rf);

  ce_DELETE(&(cov->Sce));
  localCE_DELETE(&(cov->SlocalCE));
  approxCE_DELETE(&(cov->SapproxCE));
  direct_DELETE(&(cov->Sdirect));
  hyper_DELETE(&(cov->Shyper));  
  mixed_DELETE(&(cov->Smixed));
  nugget_DELETE(&(cov->Snugget));

  plus_DELETE(&(cov->Splus));
  sequ_DELETE(&(cov->Ssequ));
  //  SPECTRAL_DELETE(&(cov->Sspectral));
  trend_DELETE(&(cov->Strend));
  tbm_DELETE(&(cov->Stbm));
  br_DELETE(&(cov->Sbr));
  pgs_DELETE(&(cov->Spgs));
  set_DELETE(&(cov->Sset));
  polygon_DELETE(&(cov->Spolygon));
  rect_DELETE(&(cov->Srect));
  dollar_DELETE(&(cov->Sdollar));
  gatter_DELETE(&(cov->Sgatter));
  earth_DELETE(&(cov->Searth));
  extra_DELETE(&(cov->Sextra));
  biwm_DELETE(&(cov->Sbiwm));
  inv_DELETE(&(cov->Sinv));
  scatter_DELETE(&(cov->Sscatter));
  mcmc_DELETE(&(cov->Smcmc));
  get_DELETE(&(cov->Sget));
  //  SELECT_DELETE(&(cov->Sselect));
  gen_DELETE(&(cov->Sgen));  
  
  simu_type *simu = &(cov->simu);
  simu->active = simu->pair = false;
  simu->expected_number_simu = 0;

  UNCONDFREE(*Cov);
}

void COV_DELETE_WITHOUT_LOC(cov_model **Cov) { 
  cov_model *cov = *Cov;
  int i,
    nsub = CovList[cov->nr].maxsub;

  if (cov == NULL) {
    warning("*cov is NULL in COV_DELETE");
    return;
  }

  // PMI(cov, "delete");
  assert(cov != NULL);

  for (i=0; i<MAXPARAM; i++) { // cov->sub[i] // seit 1.10.07: luecken erlaubt
    //                         bei PSgen !
    if (cov->kappasub[i] != NULL) {
      //   printf("%d %d %ld\n",i, MAXPARAM, cov->kappasub[i]);
      //PMI(cov);
      //assert(false);
      COV_DELETE_WITHOUT_LOC(cov->kappasub + i);
    }
  }

  for (i=0; i<nsub; i++) { // cov->sub[i] // seit 1.10.07: luecken erlaubt
    //                         bei PSgen !
    // Achtung nach nsub koennen "blinde" Untermodelle kommen, die
    // nicht geloescht werden duerfen!!
    if (cov->sub[i] != NULL) COV_DELETE_WITHOUT_LOC(cov->sub + i);
  }
  COV_DELETE_WITHOUTSUB(Cov);
}


void COV_DELETE(cov_model **Cov) { 
  cov_model *cov = *Cov;

  //  if (*Cov == NULL) crash(*Cov);

  assert(*Cov != NULL);

  //    BUG;PMI(cov);
  // printf("del %s\n", NAME(cov));

  if (cov->calling == NULL) LOC_DELETE(&(cov->prevloc));
  COV_DELETE_WITHOUT_LOC(Cov);
  *Cov = NULL;
}


void SIMU_NULL(simu_type *simu) {
  simu->pair = simu->active = false;
  simu->expected_number_simu =0;
}


static int CovZaehler = 0;
void COV_ALWAYS_NULL(cov_model *cov) {
  cov->zaehler = CovZaehler++;
  cov->calling = NULL;
  cov->prevloc = cov->ownloc = NULL;
  cov->checked = false;

  cov->MLE = NULL;
  cov->key = NULL;
  cov->rf = NULL;

  cov->mpp.mM = cov->mpp.mMplus = NULL;
  cov->mpp.moments = MISMATCH;
  
  cov->Sce = NULL;
  cov->SlocalCE = NULL;
  cov->SapproxCE = NULL;
  cov->Sdirect = NULL;
  cov->Shyper = NULL;
  cov->Smixed = NULL;
  cov->Snugget = NULL;
  cov->Splus = NULL;
  cov->Ssequ = NULL;
  cov->Stbm = NULL;
  cov->Strend = NULL;
  cov->Sbr = NULL;
  cov->Sget = NULL;
  cov->Spgs = NULL;
  cov->Sset = NULL;
  cov->Spolygon = NULL;
  cov->Srect = NULL;
  cov->Sdollar = NULL;
  cov->Sgatter = NULL;
  cov->Searth = NULL;
  cov->Sextra = NULL;
  cov->Sbiwm = NULL;
  cov->Sinv = NULL;
  cov->Sscatter = NULL;
  cov->Smcmc = NULL;
  //cov->Sselect = NULL;
  cov->Sgen = NULL;

  cov->fieldreturn = cov->origrf = false;
  cov->initialised = false;
}

void COV_NULL(cov_model *cov) {
  int i, j;
  COV_ALWAYS_NULL(cov);

  cov->nr = cov->gatternr =  cov->secondarygatternr = MISMATCH;
  for (i=0; i<MAXPARAM; i++) {
    PtoNULL(i);
    cov->ncol[i] = cov->nrow[i] = 0; 
    cov->kappasub[i] = NULL;
  }
  for (i=0; i<MAXTAYLOR; i++) {
    for (j=(int) TaylorConst; j<=(int) TaylorPow; j++) 
      cov->taylor[i][j] = cov->tail[i][j] = 0.0;
     for ( ; j<=(int) TaylorExpPow; j++) 
       cov->tail[i][j] = 0.0;
  }

  cov->q = NULL;
  cov->qlen = 0;
  for (i=0; i<MAXSUB; i++) cov->sub[i] = NULL;
  cov->nsub = 0;
  cov->user_given = ug_internal; // to know which models are given by the
  // user and which ones have been added internally

  cov->role = ROLE_UNDEFINED;
  cov->typus = UndefinedType;
  cov->method = Forbidden;
  cov->tsdim = cov->xdimprev = cov->xdimown =cov->maxdim = cov->xdimgatter = 
    UNSET;
  cov->vdim[0] = cov->vdim[1] = MISMATCH;

 
  cov->ownkappanames = NULL;

  cov->domown = cov->domprev = DOMAIN_MISMATCH;
  cov->isoown = cov->isoprev = ISO_MISMATCH;
  cov->logspeed = RF_NA;
  cov->delflag = 0;
  cov->full_derivs = cov->rese_derivs = MISMATCH;
  cov->ptwise_definite = pt_undefined;

  cov->deterministic = true;
  cov->monotone = MISMATCH; 
  //  cov->total_n = MISMATCH;
  //  cov->total_sum = 0.0;
  cov->finiterange = cov->loggiven =  
    // cov->diag = cov->semiseparatelast  = cov->separatelast = 
    cov->matrix_indep_of_x = false;
  cov->hess = cov->tbm2num = NOT_IMPLEMENTED;
  for (i=0; i<=Nothing; i++) cov->pref[i] = PREF_BEST;
  // mpp und extrG muessen auf NONE sein !
  for (; i<=Forbidden; i++) cov->pref[i] = PREF_NONE;

  cov->taylorN = cov->tailN = 0;

  MPPPROPERTIES_NULL(&(cov->mpp));
  SIMU_NULL(&(cov->simu));

  // cov->localcov=NULL;
}



void FFT_NULL(FFT_storage *FFT) 
{
  if (FFT == NULL) return;
  FFT->work = NULL;
  FFT->iwork = NULL;
}

void FFT_destruct(FFT_storage *FFT)
{
  FREE(FFT->iwork);
  FREE(FFT->work);
  FFT_NULL(FFT);
}


void ce_DELETE(ce_storage **S) {
  ce_storage *x = *S;
  if (x != NULL) {
    int l,       
      vdim = x->vdim,      
      vdimSQ = vdim * vdim;
    if (x->c!=NULL) {
      for(l=0; l<vdimSQ; l++) FREE(x->c[l]);
      UNCONDFREE(x->c);
    }
    if (x->d!=NULL) {
      for(l=0; l<vdim; l++) FREE(x->d[l]);
      UNCONDFREE(x->d);
    }
    FFT_destruct(&(x->FFT));
    FREE(x->aniso);
    FREE(x->gauss1);
    FREE(x->gauss2);
    UNCONDFREE(*S);
  }
}

void ce_NULL(ce_storage* x){
  if (x == NULL) return;
  FFT_NULL(&(x->FFT));  
  x->positivedefinite = FALSE;
  x->trials = -1;
  x->c = x->d = NULL;
  x->smallestRe = x->largestAbsIm = RF_NA;
  x->aniso = NULL;
  x->stop = false;
  x->gauss1 = x->gauss2 = NULL;
//  int i;
  //for (i=0; i<MAXCEDIM; i++) x->[i] = NULL;
}



void localCE_DELETE(localCE_storage**S) 
{
  localCE_storage* x = *S;
  if (x!=NULL) {
    FREE(x->correction);
    UNCONDFREE(*S);
  }
}

void localCE_NULL(localCE_storage* x){
  // int i;  for (i=0; i<MAXCEDIM; i++) 
  if (x == NULL) return;
  x->correction = NULL;
}


void approxCE_DELETE(approxCE_storage **S) {
  approxCE_storage* x = * S;
  if (x != NULL) {
    FREE(x->idx);
    UNCONDFREE(*S);
  }
}

void approxCE_NULL(approxCE_storage* x){
  if (x == NULL) return;
  x->idx = NULL;
}


void direct_DELETE(direct_storage  ** S) {
  direct_storage *x = *S;
  if (x!=NULL) {
    FREE(x->U); 
    FREE(x->G);
    UNCONDFREE(*S);
  }
}

void direct_NULL(direct_storage  *x) {
  if (x == NULL) return;
  x->U = NULL;
  x->G = NULL;
}




void hyper_DELETE(hyper_storage  **S) {
  hyper_storage *x = *S; 
  if (x != NULL) {
    UNCONDFREE(*S);
  }
}

void hyper_NULL(hyper_storage* x) {
  if (x == NULL) return;
}


void mixed_DELETE(mixed_storage ** S) {
  mixed_storage *x = *S;
  if (x!=NULL) {
    FREE(x->Xb);
    FREE(x->mixedcov);
    UNCONDFREE(*S);
  }
}

void mixed_NULL(mixed_storage* x) {
  x->Xb = NULL;
  x->mixedcov = NULL;
  x->initialized = false;
}


void nugget_NULL(nugget_storage *x){
  if (x == NULL) return;
  x->pos = NULL;
  x->red_field = NULL;
}

void nugget_DELETE(nugget_storage ** S)
{
  nugget_storage *x = *S;
  if (x != NULL) {
    FREE(x->pos);
    FREE(x->red_field);
    UNCONDFREE(*S);
  }
}

void plus_NULL(plus_storage *x){
  if (x != NULL) {
    int i;
    for (i=0; i<MAXSUB; i++) x->keys[i] = NULL;
  }
}

void plus_DELETE(plus_storage ** S){
  plus_storage *x = *S;
  if (x != NULL) {
    int i;
    for (i=0; i<MAXSUB; i++)
      if (x->keys[i] != NULL) COV_DELETE(x->keys + i);
    UNCONDFREE(*S);
  }
}

void sequ_NULL(sequ_storage *x){
  if (x == NULL) return;
  x->U11 = NULL;
  x->U22 = NULL;
  x->MuT = NULL;
  x->G = NULL;
  x->Cov21 = NULL;
  x->Inv22 = NULL;
  x->res0 = NULL;
}

void sequ_DELETE(sequ_storage ** S){
  sequ_storage *x = *S;
  if (x!=NULL) {
    FREE(x->U11); 
    FREE(x->U22); 
    FREE(x->MuT); 
    FREE(x->G);
    FREE(x->Cov21);
    FREE(x->Inv22);
    FREE(x->res0);
    UNCONDFREE(*S);
  }
}

void spectral_DELETE(spectral_storage **S) 
{ 
  spectral_storage *x = *S;
  if (x!=NULL) {
    // do NOT delete cov --- only pointer
      // spectral_storage *x; x =  *((spectral_storage**)S);
    UNCONDFREE(*S);   
  }
}


void spectral_NULL(spectral_storage *x) 
{ 
  if (x!=NULL) {
  }
}


void trend_DELETE(trend_storage ** S) {
  trend_storage *x = *S;
  if (x!=NULL) {
    FREE(x->x);
    FREE(x->xi);
    FREE(x->evalplane);
    FREE(x->powmatrix);    
    UNCONDFREE(*S);
  }
}

void trend_NULL(trend_storage* x) {
  x->x = NULL;
  x->xi = NULL;
  x->evalplane = NULL;
  x->powmatrix = NULL;
}

void tbm_DELETE(tbm_storage **S) 
{
  tbm_storage *x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}

void tbm_NULL(tbm_storage* x) {
  if (x == NULL) return;
  x->ce_dim =  0;
  x->err = ERRORFAILED;
}

void BRTREND_DELETE(double **BRtrend, int trendlen) {
   int j;
   if (BRtrend == NULL) return;
   for (j=0; j<trendlen; j++) FREE(BRtrend[j]);
 }

void br_DELETE(br_storage **S) {
  br_storage *brS = *S;  
  if (brS != NULL) {
    int i;
    if (brS->trend != NULL) {
      BRTREND_DELETE(brS->trend, brS->trendlen); 
      UNCONDFREE(brS->trend);
    }
    FREE(brS->shiftedloc);     
    FREE(brS->loc2mem);

    
    if (brS->countvector != NULL) {
      for (i=0; i<brS->vertnumber; i++) FREE(brS->countvector[i]);
      UNCONDFREE(brS->countvector);
    }
 
    if (brS->areamatrix != NULL) {
      for (i=0; i<brS->vertnumber; i++) FREE(brS->areamatrix[i]);
      UNCONDFREE(brS->areamatrix);
    }
    FREE(brS->logvertnumber);
    FREE(brS->locindex);
    FREE(brS->suppmin);
    FREE(brS->suppmax);
    FREE(brS->locmin);
    FREE(brS->locmax);
    FREE(brS->loccentre);
    
    FREE(brS->mem2loc);
    FREE(brS->newx);
    if (brS->vario != NULL) COV_DELETE(&(brS->vario));
    for(i=0; i<MAXSUB; i++) {
      FREE(brS->lowerbounds[i]);
      if (brS->sub[i] != NULL) COV_DELETE(brS->sub + i);
    }
    if (brS->submodel != NULL) COV_DELETE(&(brS->submodel));
    UNCONDFREE(*S);
  }
}

void br_NULL(br_storage *brS) {
  brS->trend = NULL;
  brS->trendlen = 0;
  brS->next_am_check = NA_INTEGER;
  brS->shiftedloc = NULL;
  brS->mem2loc = brS->locindex = brS->loc2mem = NULL;
  brS->memcounter = 0;
  brS->newx = brS->suppmin = brS->suppmax = 
    brS->locmin = brS->locmax = brS->loccentre = NULL;
  int i;
  for (i=0; i<MAXSUB; i++) {
    brS->lowerbounds[i] = NULL;
    brS->radii[i] = RF_NA;
    brS->thresholds[i] = RF_INF;
    brS->sub[i] = NULL;
    brS->zeropos[i] = NA_INTEGER;
  }
  brS->vario = brS->submodel = NULL;  
  brS->idx = 0;
  brS->countvector = NULL;
  brS->areamatrix = NULL;
  brS->logvertnumber = NULL;
}


void pgs_DELETE(pgs_storage **S) 
{
  pgs_storage *x = *S;
  if (x!=NULL) {

    // Huetchen
    FREE(x->v);
    FREE(x->y);
    FREE(x->xgr[0]);
    
    FREE(x->pos);
    FREE(x->min);
    FREE(x->max); 
    
    FREE(x->single);
    FREE(x->total);
    FREE(x->halfstepvector);    
    
    FREE(x->localmin);
    FREE(x->localmax);
    FREE(x->minmean);
    FREE(x->maxmean);
    
    
    // rf_interface.cc
    FREE(x->supportmin);
    FREE(x->supportmax);
    FREE(x->supportcentre);   
    FREE(x->own_grid_start);
    FREE(x->own_grid_step);
    FREE(x->own_grid_len);
    
    FREE(x->gridlen);
    FREE(x->end);
    FREE(x->start);
    FREE(x->delta);
    FREE(x->nx);

    FREE(x->xstart);    
    FREE(x->x);
    FREE(x->inc);    
 

    // variogramAndCo.cc:
    FREE(x->endy);
    FREE(x->startny);
    FREE(x->ptrcol);
    FREE(x->ptrrow);
    FREE(x->C0x);
    FREE(x->C0y);
    FREE(x->cross);
    FREE(x->z);
    FREE(x->Val);
    
    if (x->cov != NULL) {
      cov_model *dummy = x->cov;
      if (x->cov->Spgs != NULL && x->cov->Spgs->cov != NULL && 
	  x->cov->Spgs->cov->Spgs == x) {
	 x->cov->Spgs->cov = NULL;// um "Unendlich"-Schleife zu vermeiden!
      } else {
	assert(x->cov->Spgs == NULL || x->cov->Spgs->cov == NULL || 
	       x->cov->Spgs->cov->Spgs == x);
      }
      x->cov = NULL; // Sicherheitshalber
      COV_DELETE(&(dummy));
    }

    UNCONDFREE(*S);
  }
}


void pgs_NULL(pgs_storage* x) {
  if (x == NULL) return;
   x->log_density = x->intensity = x->globalmin = x->currentthreshold = 
    RF_NEGINF;
  x->flat = x->estimated_zhou_c = x->logmean = false;
  x->totalmass = RF_NA;
  x->alpha = 1.0; // !! standard

  x->size = -1;
  x->zhou_c = RF_NA;
  x->sq_zhou_c = x->sum_zhou_c = 0;
  x->n_zhou_c = 0;
 

  // Huetchen
  x->v = x->y = x->xgr[0] = NULL;
  x->pos = x->min = x->max = NULL;
  x->single = x->total = x->halfstepvector = NULL;  
 
  x->localmin = x->localmax =x->minmean = x->maxmean = NULL;


  // rf_interface.cc
  x->supportmin = x->supportmax = x->supportcentre = //
    x->own_grid_start = x->own_grid_step =  x->own_grid_len = NULL; //

  x->gridlen = x->end = x->start = x->delta = x->nx = NULL; //

  x->xstart = x->x =  x->inc = NULL;


  // variogramAndCo.cc
  x->endy = x->startny = x->ptrcol = x->ptrrow = NULL;
  x->C0x = x->C0y = x->cross = x->z = NULL;
  x->Val = NULL;

  //x->param_set = NULL;
  x->cov = NULL;

}


void set_DELETE(set_storage **S) {
  set_storage *x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}

void set_NULL(set_storage* x) {
  if (x == NULL) return;
  x->remote = NULL;
  x->set = NULL;
  x->variant = 0;
}

void polygon_DELETE(polygon_storage **S) 
{
  polygon_storage *x = *S;
  if (x != NULL) {
    if (x->vdual != NULL) {
      int i;
      for (i=0; i<x->n_vdual; i++) FREE(x->vdual[i]);
      UNCONDFREE(x->vdual);
    }
    FREE(x->vprim);
    if (x->P != NULL) {
      freePolygon(x->P);
      UNCONDFREE(x->P);
    }
  }
  UNCONDFREE(*S);
}

void polygon_NULL(polygon_storage* x) {
  if (x == NULL) return;
  x->vprim = NULL;
  x->vdual = NULL;
  x->n_vdual = x->n_vertex = x->n_v = 0;

  polygon *P = x->P;
  if (P == NULL) BUG;
  P->e = NULL;
  P->v = NULL;
  P->n = 0;
}

void rect_DELETE(rect_storage **S){
  rect_storage *x = *S;
  if (x!=NULL) {
    FREE(x->value);
    FREE(x->weight);
    FREE(x->tmp_weight);
    FREE(x->right_endpoint);
    FREE(x->ysort);
    FREE(x->z);
    
    FREE(x->squeezed_dim);
    FREE(x->assign);
    FREE(x->i);

    UNCONDFREE(*S);
  }
}

void rect_NULL(rect_storage* x) {
  if (x == NULL) return;
  x->tmp_n = x->nstep = -999;
  x->value = x->weight = x->tmp_weight = 
    x->right_endpoint = x->ysort = x->z = NULL;


  // for PMI
  //  x->inner = x->inner_const = x->inner_pow = x->outer =
  //    x->outer_const = x->outer_pow = x->outer_pow_const = x->step = RF_NA;

  x->squeezed_dim = x->assign = x->i = NULL;
}


void dollar_DELETE(dollar_storage **S) 
{
  dollar_storage *x = *S;
  if (x!=NULL) {
    FREE(x->z);
    FREE(x->z2);
    FREE(x->y);
    FREE(x->y2);
    FREE(x->save_aniso);
    FREE(x->inv_aniso);
    FREE(x->nx);
    FREE(x->len);
    FREE(x->total);
    FREE(x->cumsum);
    UNCONDFREE(*S);
  }
}

void dollar_NULL(dollar_storage* x) {
  if (x == NULL) return;
  x->z = x->z2 = x->y =  x->y2 = x->save_aniso = x->inv_aniso = NULL;
  x->cumsum = x->nx = x->total = x->len = NULL;
}



void gatter_DELETE(gatter_storage **S) 
{
  gatter_storage *x = *S;
  if (x!=NULL) {
    FREE(x->z);
    UNCONDFREE(*S);
  }
}

void gatter_NULL(gatter_storage* x) {
  if (x == NULL) return;
  x->z = NULL;
  //  x->zsys = NULL;
}


void earth_DELETE(earth_storage **S) 
{
  earth_storage *x = *S;
  if (x!=NULL) {
    FREE(x->X);
    FREE(x->Y);
    FREE(x->U);
    FREE(x->V);
    UNCONDFREE(*S);
  }
}

void earth_NULL(earth_storage* x) {
  if (x == NULL) return;
  x->X = x->Y = x->U = x->V = NULL;
}



void extra_DELETE(extra_storage **S) 
{
  extra_storage *x = *S;
  if (x!=NULL) {
    FREE(x->a);
    FREE(x->b);
    FREE(x->c);
    UNCONDFREE(*S);
  }
}

void extra_NULL(extra_storage* x) {
  if (x == NULL) return;
  x->a = x->b = x->c = NULL;
}



void biwm_DELETE(biwm_storage **S) 
{
  biwm_storage *x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}

void biwm_NULL(biwm_storage* x) {
  if (x == NULL) return;  
}


void inv_DELETE(inv_storage **S) {
  inv_storage *x = *S;
  if (x!=NULL) {
    FREE(x->v);
    FREE(x->wert);
    UNCONDFREE(*S);
  }
}

void inv_NULL(inv_storage* x) {
  if (x == NULL) return;  
  x->v = NULL;
  x->wert = NULL;
}


void scatter_DELETE(scatter_storage **S) {
  scatter_storage *x = *S;
  if (x!=NULL) {
    FREE(x->nx);
    FREE(x->min);
    FREE(x->max);
    FREE(x->step);
    FREE(x->x);
    FREE(x->xmin);
    FREE(x->value);
    UNCONDFREE(*S);
  }
}

void scatter_NULL(scatter_storage* x) {
  if (x == NULL) return; 
  x->dim = x->vdim = -1; 
  x->min = x->max = x->nx = NULL;
  x->step = x->value = x->xmin = x->x = NULL;
}

void mcmc_DELETE(mcmc_storage **S) {
  mcmc_storage *x = *S;
  if (x!=NULL) {
    FREE(x->pos);
    FREE(x->deltapos);
    FREE(x->propdelta);
    FREE(x->proposed);
    UNCONDFREE(*S);
  }
}

void mcmc_NULL(mcmc_storage* x) {
  if (x == NULL) return; 
  x->pos = x->deltapos = x->propdelta = x->proposed = NULL;
  //x->done = false;
  x->integral = 0.0;
}


void get_NULL(get_storage *s){
  s->orig = s->get_cov = NULL;
  s->idx = NULL;
  s->param_nr = -1;
}

void get_DELETE(get_storage **S){
  get_storage *x = *S;
  if (x != NULL) {
    FREE(x->idx);
    UNCONDFREE(*S);
  }
}


void gen_NULL(gen_storage *x) {
  int d;
  if (x == NULL) return;
  //  x->mpp.newx = NULL;

  // for (d=0; d<MAXMPPDIM; d++) 
  //   x->window.min[d] = x->window.max[d] = x->window.centre[d] = RF_NA;
  x->check = true;

  x->Sspectral.phistep2d = x->Sspectral.phi2d = x->Sspectral.prop_factor
    = RF_NA;
  x->Sspectral.grid = false;

  x->spec.nmetro = -1;
  x->spec.sigma = -1.0;
  x->spec.density = NULL;
  for (d=0; d<MAXTBMSPDIM; d++) {
    x->spec.E[d] = x->spec.sub_sd_cum[d] = RF_NA;
  }
}

void gen_DELETE(gen_storage **S) {
  gen_storage *x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}

void addModel(cov_model **pcov, int covnr, cov_model *calling, bool nullOK) {
  cov_model *cov;
  int i;
  
  // printf("addmodel: %s %ld\n", CovList[covnr].name, (long int) *pcov);

  cov = (cov_model*) MALLOC(sizeof(cov_model));
  COV_NULL(cov);
  //assert(cov->calling == NULL);
  cov->nr = covnr;

  if (*pcov != NULL) {
    cov->nsub = 1;
    cov->calling = (*pcov)->calling;
    (*pcov)->calling = cov;
    cov->sub[0] = *pcov;
    for (i=0; i<=Forbidden; i++) {	
      // cov->user[i] = cov->sub[0]->user[i];
      cov->pref[i] = cov->sub[0]->pref[i];
    }
  }


  //  if (calling != NULL)   
  //printf("addmodel calling %s cov %s\n", NAME(calling), NAME(cov));


  if (calling != NULL) cov->calling = calling;
  else if (!nullOK && *pcov == NULL) {
    PRINTF("Missing link for model '%s'.\n", NICK(cov));
    BUG;
  }
  
  *pcov = cov;
}

void addModel(cov_model **pcov, int covnr) {
  addModel(pcov, covnr, NULL, false);
}

void addModel(cov_model *pcov, int subnr, int covnr) {
  assert(pcov != NULL);
  bool newsub = pcov->sub[subnr] == NULL;
  addModel(pcov->sub + subnr, covnr, pcov, false);
  pcov->nsub += (int) newsub;
}

void addModelKappa(cov_model *pcov, int subnr, int covnr) {
  assert(pcov != NULL);
  addModel(pcov->kappasub + subnr, covnr, pcov, false);
}

void addModel(cov_model **pcov, int covnr, cov_model *calling) {
  addModel(pcov, covnr, calling, false);
}

int setgrid(coord_type xgr, double *x, long lx, int spatialdim) {
  if (lx!=3) {      
    //    PRINTF("%d\n", lx);
    //    crash();
    SERR("Problem with the coordinates (non-integer number of locations or non-positive step)")      
  }
  
  int d;
  unsigned long
    totalBytes = sizeof(double) * lx * spatialdim;
 
  if (xgr[0] == NULL && (xgr[0] =(double*) MALLOC(totalBytes))==NULL)
    return ERRORMEMORYALLOCATION; 

  //printf("setgrid %d\n", totalBytes);
  MEMCOPY(xgr[0], x, totalBytes);
  
  // folgende Zeile nur beim ersten Mal zu setzen, aber egal
  for (d=1; d<spatialdim; d++) {
    xgr[d] = &(xgr[0][d * lx]); 
    if (xgr[d][XLENGTH] != (int) xgr[d][XLENGTH])
      SERR2("grid length must be integer valued. Got %f in dimension %d.",
	    xgr[d][XLENGTH], d);
    if (xgr[d][XLENGTH] < 1.0) 
      SERR2("grid length must be positive. Got %f in dimension %d.",
	    xgr[d][XLENGTH], d);
  }
  for (; d<MAXSIMUDIM; d++)  xgr[d]=NULL;    
  
  return NOERROR;
}


int add_y_zero(location_type *loc) {
  int d;
  if (loc->ly > 0) BUG;
  if (loc->distances) { SERR("distances are allowed only for cartesian systems"); } 
  if (loc->grid) {
    loc->ly = 3;
    double *ygrid = (double*) MALLOC(sizeof(double) * loc->ly *loc->spatialdim);
    int i;
    for (i=d=0; d<loc->spatialdim; d++) {
      ygrid[i++] = 0.0;
      ygrid[i++] = 0.0;
      ygrid[i++] = 1.0;
    }
    setgrid(loc->ygr, ygrid, loc->ly, loc->spatialdim);
  } else {
    loc->ly = 1;
    if ((loc->y=(double*) CALLOC(loc->ly * loc->xdimOZ, sizeof(double) ))==NULL)
      return ERRORMEMORYALLOCATION; 
  } 

  if (loc->Time) {
    if (loc->grid) loc->ygr[loc->spatialdim] = loc->T;
  }

  return NOERROR;
}

int partial_loc_set(location_type *loc, double *x, double *y,
		    long lx, long ly, bool dist, int xdimOZ, double *T,
		    bool grid, bool cpy){
  int d, err;
  unsigned long totalBytes;

  if (((loc->x != NULL) && ((loc->y == NULL) xor (ly==0))) ||
      ((loc->xgr[0] != NULL) && ((loc->ygr[0] == NULL) xor (ly==0)))) {
    // 28783440 51590192 0; 0 0 0 dist=0

   //printf("%ld %ld %d; %ld %ld %d dist=%d\n",
    //  (long int) loc->x, (long int)loc->y, ly,
    //     (long int) loc->xgr[0], (long int) loc->ygr[0], ly, loc->distances);


  // assert(false);
  // 
  //crash();
  SERR("domain structure of the first and second call do not match");
  }

  assert(x != NULL);

  loc->xdimOZ = xdimOZ; // ohne Zeit !!
  loc->lx = lx;
  loc->ly = ly;
  if (ly  > 0) { 
    //crash();
    if (dist) { SERR("distances are not allowed if y is given"); } 
  }
  
  loc->grid = grid;
  loc->distances = dist;
  if (loc->delete_x) {
    if (loc->y != NULL) {
      if (loc->y != loc->x) UNCONDFREE(loc->y); 
    }
    FREE(loc->x);
  }
  loc->delete_x = cpy;
  if (lx == 0) return NOERROR;
  
  if (grid) {
    if ((err = setgrid(loc->xgr, x, lx, loc->spatialdim)) !=NOERROR) return err;
    if (ly>0) {
      if (x == y) for(d=0; d<loc->spatialdim; d++) loc->ygr[d] = loc->xgr[d];
      else {
	if ((err = setgrid(loc->ygr, y, ly, loc->spatialdim)) !=NOERROR) 
	  return err;
      }
    }
    
    for (loc->spatialtotalpoints=1, d=0; d<loc->spatialdim; d++) {
      loc->spatialtotalpoints *= (int) loc->xgr[d][XLENGTH];
    }
    loc->totalpoints = loc->spatialtotalpoints;
  } 

  else if (dist) {
    // lx : anzahl der Puntke involviert, d.h. x muss die Laenge lx ( lx -1) / 2
    // haben.
    if (lx > 0) {
      if (cpy) {
	totalBytes =  sizeof(double) * lx * (lx - 1) / 2 * xdimOZ;
	if ((loc->x=(double*) MALLOC(totalBytes))==NULL){
	  return ERRORMEMORYALLOCATION; 
	}
	MEMCOPY(loc->x, x, totalBytes);
      } else {
	loc->x = x;
      }
    }
    loc->totalpoints = loc->spatialtotalpoints = lx;
    //     (int) (1e-9 + 0.5 * (1.0 + sqrt(1.0 + 8.0 * lx)));
    //   if (0.5 * (loc->totalpoints * (loc->totalpoints - 1.0)) != lx) {   
    //printf("tot=%d %d %d\n", loc->totalpoints, (int) ( 0.5 * (loc->totalpoints * (loc->totalpoints - 1.0))), (int) lx); assert(false);
    //   SERR("distances do not have expected size");
    //}
  }
  
  else { // not grid, not distances
    if (cpy) {
      totalBytes =  sizeof(double) * lx * loc->xdimOZ;
      
      //       printf("totalbyets %d %ld %d\n", totalBytes, lx, loc->xdimOZ);
      
      if ((loc->x=(double*) MALLOC(totalBytes))==NULL){
	return ERRORMEMORYALLOCATION; 
      }

      //printf("%d %d\n", loc->x==NULL, x==NULL);

      MEMCOPY(loc->x, x, totalBytes);
      if (loc->ly>0) {
	if (x == y) {
	  loc->y = loc->x;
	} else {
	  totalBytes =  sizeof(double) * ly * loc->xdimOZ;
	  //printf("totalbytes y %d %ld %d\n", totalBytes, ly, loc->xdimOZ);
	  if ((loc->y=(double*) MALLOC(totalBytes))==NULL){
	    return ERRORMEMORYALLOCATION; 
	  }
	  MEMCOPY(loc->y, y, totalBytes);	
	}
      }
    } else {
      loc->x = x;
      loc->y = y;
    }
    
    //printf("getnset %d\n", lx);
    
    loc->totalpoints = loc->spatialtotalpoints = lx;
    // just to say that considering these values does not make sense
  }
  
  if ((loc->Time) xor (T!=NULL)) {    
    //cov_model *cov;  crash(cov);
    //  AERR(1);
    SERR("partial_loc: time mismatch");
  }
  if (loc->Time) {
    MEMCOPY(loc->T, T, sizeof(double) * 3);
    if (grid) {
      loc->xgr[loc->spatialdim] = loc->T;
      if (ly>0) loc->ygr[loc->spatialdim] = loc->T;
    }
    
    if (loc->T[XLENGTH] <= 0) {
      //crash();
      SERR1("The number of temporal points is not positive. Check the triple definition of 'T' in the man pages of '%s'.", CovList[SIMULATE].nick)
	}
    loc->totalpoints *= (int) loc->T[XLENGTH];
  }

  return NOERROR;
}

int loc_set(double *x, double *y, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdimOZ,
	    long lx, long ly, bool Time, bool grid,
	    bool distances,
	    location_type **Loc) {
  int err;
  //unsigned long totalBytes;
  // preference lists, distinguished by grid==true/false and dimension
  // lists must end with Nothing!
 
  if (*Loc != NULL) LOC_DELETE(Loc);
 
  location_type *loc = *Loc = (location_type*) MALLOC(sizeof(location_type));
  LOC_NULL(loc);

  loc->timespacedim = spatialdim + (int) Time;
  loc->spatialdim = spatialdim;
  loc->Time = Time; 
  if (spatialdim<1 || loc->timespacedim>MAXSIMUDIM) return ERRORDIM;

  //printf("loc %ld %d\n", lx, xdimOZ);
  assert(x != NULL);
  
  if ((err = partial_loc_set(*Loc, x, y, lx, ly, distances, xdimOZ,
			     Time ? T : NULL,
			     grid, true)) != NOERROR) XERR(err);
 
  return NOERROR;
}


int loc_set(cov_model *cov, long totalpoints) {
  location_type *loc;
  if (cov->ownloc == NULL) {
    loc = cov->ownloc = (location_type*) MALLOC(sizeof(location_type));
    LOC_NULL(loc);
    loc->delete_x = false;
  } else {
    loc = Loc(cov);
    if (loc->xgr[0] != NULL || loc->x != NULL) BUG;
  }
  cov->ownloc->totalpoints = totalpoints;
  return NOERROR;
}

int loc_set(double *x, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdimOZ, /* original ! */
	    long lx, bool Time, bool grid,
	    bool distances,
	    location_type **Loc) {
  return loc_set(x, NULL, T, spatialdim, xdimOZ, lx, 0, Time, grid,
		 distances, Loc);    
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
  bool notquasidiag=true, *taken=NULL;
  int j, k, startidx, i, last;

  if (aniso == NULL) {
    *diag = *quasidiag = *separatelast = *semiseparatelast = true;
    for (i=0; i<col; i++) idx[i] = i;
    return;
  }

  taken = (bool *) MALLOC(row * sizeof(bool));

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
    /*
     * * 0
     * * 0
     * * *
     */
    last = col * row - 1;
    for (k=last - row + 1; k<last; k++)
      if (!(*separatelast = aniso[k] == 0.0)) break;
  }
  if (!(*separatelast = *semiseparatelast)) {
    /*
     * * 0
     * * 0
     0 0 *
     */  
    for (k=row - 1; k<last; k+=row)
      if (!(*separatelast = aniso[k] == 0.0)) break;
  }
  UNCONDFREE(taken);
}

int getmodelnr(char *name) {
  // == -1 if no matching function is found
  // == -2 if multiple matching fctns are found, without one matching exactly
  // == -3 if internal name is passed
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)
  int match;

  if (currentNrCov==-1) { InitModelList(); }
  if (!strcmp(name, InternalName)) return MATCHESINTERNAL;
  if ((match = Match(name, CovNickNames, currentNrCov)) >= 0) return match;
  return Match(name, CovNames, currentNrCov);
}

void MultiDimRange(cov_model *cov, double *natscale) {
  int wave, i, redxdim, d, idx, 
    xdimprev = cov->xdimprev,
    vdim = cov->vdim[0],
    err = NOERROR;
  double y, yold, x[MAXGETNATSCALE], threshold, natsc, factor, sign,
    newx, xsave, newy, 
    *dummy =  NULL,
    rel_threshold = 0.05;

  
  bool domain = cov->domown ==XONLY;
  //  covfct cf=NULL;
  // nonstat_covfct ncf=NULL;
  
  assert(MAXGETNATSCALE <= MAXCOVDIM); // ZERO
  
  redxdim = cov->xdimown;
         
  if (redxdim > MAXGETNATSCALE) {
    err = -1; goto ErrorHandling;
  }
 
  if ((dummy = (double*) MALLOC(sizeof(double) * vdim * vdim)) == NULL) {
    err = -2; goto ErrorHandling;
  }
 
  if (cov->full_derivs < 0) { err=ERRORNOTDEFINED; goto ErrorHandling; }
  if (domain) {
    COV(ZERO, cov, dummy);
  } else {
    NONSTATCOV(ZERO, ZERO, cov, dummy);
  }
  threshold = rel_threshold * dummy[0];
  
  for (d=0; d<redxdim; d++) {
    wave  = 0;
    for (i=0; i<xdimprev; i++) x[i] = 0.0;
    
    idx = (redxdim == xdimprev || d==0) ? d : xdimprev-1; 
    x[idx] = 1.0;

     
    if (domain) COV(x, cov, dummy) else NONSTATCOV(ZERO, x, cov, dummy); 
    yold = dummy[0];
    if (ISNAN(yold)) { err = -3; goto ErrorHandling;}
    // print("(%f %f) thres=%f yold=%f %d\n", x[0], x[1], threshold, yold,
    //	   yold > threshold);
    if (yold > threshold) {
      factor = 2.0;
      sign = 1.0;
    } else {
      factor = 0.5;
      sign = -1.0;
    } 
    
    double otherx;
    x[idx] *= factor;
    if (domain) COV(x, cov, dummy) else NONSTATCOV(ZERO, x, cov, dummy);
    y = dummy[0];
    
    // print("%f %f : y=%f diff=%f %f\n", x[0], x[1], y, y - threshold,
    //	   sign * (y - threshold));
    while (sign * (y - threshold) > 0) {  
      if (yold<y){ 
	if (wave++>10) { err=ERRORWAVING; goto ErrorHandling; }
      }
      yold = y;
      x[idx] *= factor;
      if (x[idx]>1E30) { err=ERRORRESCALING; goto ErrorHandling; }
      if (domain) COV(x, cov, dummy) else NONSTATCOV(ZERO, x, cov, dummy);
      y = dummy[0];
    }


    otherx = x[idx] / factor;

    //   print("%f %f : y=%f diff=%f %f\n", x[0], x[1], y, y - threshold,
    //	   sign * (y - threshold));

    for (i=0; i<3 /* good choice?? */ ;i++) {       

      //   print("%d %f %f; other=%f\n", i, x[0], x[1], otherx); 
   
      if (y==yold) { err=ERRORWAVING; goto ErrorHandling; }
      newx = x[idx] + (x[idx]-otherx)/(y-yold)*(threshold-y);
      // print("newx=%f xidx=%f, oth=%f y=%f yold=%f, thr=%f\n",
      //	     newx, x[idx], otherx, y, yold, threshold);
      //  print("thres %f\n", threshold);
      xsave = x[idx];
      x[idx] = newx;
      if (domain) COV(x, cov, dummy) else NONSTATCOV(ZERO, x, cov, dummy);
      newy = dummy[0];
      x[idx] = xsave;
      
      if (sign * (newy - threshold) > 0) {
	otherx = newx;
	yold  = newy;
      } else {
	x[idx] = newx;
	y = newy;
      }
    }

    // print("%f %f\n", y, threshold);

    if (y==yold)  { err=ERRORWAVING; goto ErrorHandling; }
    natsc = 1.0 / ( x[idx] + 
		    (x[idx]-otherx)/(y-yold)*(threshold-y) );
    
    if (redxdim == xdimprev || d==0) {
      natscale[d] = natsc;
    } else {
      int beg, end;
      if (redxdim == 2) {
	if (d==0) {
	  beg = 0;
	  end = xdimprev-1;
	} else {
	  beg = end = xdimprev-1;
	}
      } else {
	beg = 0;
	end = xdimprev;
      }
      for (i=beg; i<end; natscale[i++] = natsc);
    }
  }
  
 ErrorHandling:
  FREE(dummy);
  switch(err) {
  case NOERROR : break;
  case -1 : 
    ERR("dimension of x-coordinates too high to detect natural scaling.");
  case -2 :
    ERR("not enough memory when determining natural scaling.");
  case -3 :
    ERR("NA in model evaluation detected"); 
  default: XERR(err);
  } 
}


void UserGetNatScaling(double *natscale) {
  GetNaturalScaling(KEY[MODEL_USER], natscale);
}


void GetNaturalScaling(cov_model *cov, double *natscale)
{ // called also by R 

  // values of naturalscaling:
  //#define NATSCALE_EXACT 1   
  //#define NATSCALE_APPROX 2
  //#define NATSCALE_MLE 3 /* check fitvario when changing !! */
  
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  *natscale = 0.0;

  if (C->maxsub!=0) XERR(ERRORFAILED); 
 
  if (!is_any(ISOTROPIC, C) || cov->isoown != ISOTROPIC
      || C->domain != XONLY ||
      !isPosDef(cov->typus) 
      || C->vdim != SCALAR) 
    ERR("anisotropic function not allowed");
	 
  if (C->finiterange == true) {
    *natscale = 1.0;
    return;
  }
  
  if (C->inverse!=NULL) { 
    C->inverse(&GLOBAL.gauss.approx_zero, cov, natscale);
    *natscale = 1.0 / *natscale;
    if (ISNAN(*natscale) || *natscale != 0.0) {
      return;
    }
  }
    
  if (NS != NATSCALE_ORNUMERIC) XERR(ERRORRESCALING); 

  if ((C->cov)==nugget)  XERR(ERRORRESCALING); 
      
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
  
  assert(cov->xdimown == 1); 
  MultiDimRange(cov, natscale);
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

void TaylorCopy(cov_model *to, cov_model *from) {
  int i, j;
  to->taylorN = from->taylorN;
  to->tailN = from->tailN;
  for(i=0; i<to->taylorN; i++) {
    for (j=0; j<=TaylorPow; j++) to->taylor[i][j] = from->taylor[i][j];
  }
  for(i=0; i<to->tailN; i++) {
    for (j=0; j<=TaylorExpPow; j++) to->tail[i][j] = from->tail[i][j];
  }
}


void paramcpy(cov_model *to, cov_model *from, bool freeing,
	      bool allocating, bool copy_lists, bool recursive,
	      bool copy_mpp) {
  cov_fct *C = CovList + from->nr; // nicht gatternr
  double **pto = to->px,
    **pfrom = from->px;
  int i, v,
     n = -1,
     *to_col = to->ncol,
     *to_row = to->nrow,
     *from_col = from->ncol,
     *from_row = from->nrow;

  bool same_model = abs(to->nr - from->nr) <= 1 ||
    (isDollar(to) && isDollar(from));

  if (!same_model) {
    BUG;      
  }    
  if (freeing && !allocating) BUG;

  for (i=0; i<MAXPARAM; i++) {

 
    if (pfrom[i] == NULL) {
      continue;
    }
   
    //    printf("param %s.%d %d\n", NICK(from), i, freeing);

    if (freeing) {
      // if (strcmp(to->kappanames[i], from->kappanames[i]) != 0 ||
      //	  type != from->kappatype[i]) return false;
      FREE(pto[i]);
      to_col[i] = from_col[i];
      to_row[i] = from_row[i];
    }

  
    if (C->kappatype[i] >= LISTOF) {      
      if (allocating) pto[i] = (double*) MALLOC(sizeof(listoftype));
      int j, len = from_row[i];
      listoftype *p = PARAMLIST(from, i),
	*q = (listoftype*) (pto[i]);	
      if ((q->deletelist = copy_lists)) {      
	for (j=0; j<len; j++) {
	  if (C->kappatype[i]== LISTOF + REALSXP) {
	    n = sizeof(double);
	  } else BUG;
	  n *= p->nrow[j] * p->ncol[j];
	  q->nrow[j] = p->nrow[j];
	  q->ncol[j] = p->ncol[j];
	  if (allocating) q->p[j] = (double*) MALLOC(n);
	  MEMCOPY(q->p[j], p->p[j], n);	    
	}
      } else {
 	for (j=0; j<len; j++) {
	  q->nrow[j] = p->nrow[j];
	  q->ncol[j] = p->ncol[j];
	  q->p[j] =  p->p[j];	    
	}
      }
    } else if (C->kappatype[i]==LANGSXP) {
      n = sizeof(sexp_type);
      if (allocating) pto[i] = (double*) MALLOC(n);
      MEMCOPY(pto[i], pfrom[i], n);
      ((sexp_type *) pto[i])->Delete = false;
    } else {
      if (C->kappatype[i]==CLOSXP) {
	BUG;  // Marco brauchen wir CLSXP ??
	n = sizeof(SEXP);
      } else {
	if (C->kappatype[i]==REALSXP) {
	  n = sizeof(double);
	} else if (C->kappatype[i]==INTSXP) {
	  n = sizeof(int);
	} else BUG;
	n *= from_row[i] * from_col[i];
      }
      
      //           printf("i=%d to=%s %s from=%s %s %d\n", i, 
      //NICK(to), CovList[to->nr].kappanames[i], 
      //	     NICK(from), CovList[from->nr].kappanames[i], n);
      
      if (allocating) pto[i] = (double*) MALLOC(n); 
      MEMCOPY(pto[i], pfrom[i], n);
    }
  }

  if (copy_mpp) {
    //  PMI(to); PMI(from);
    if (to->mpp.moments < 0 && alloc_mpp_M(to, from->mpp.moments)!=NOERROR) 
      error("error in allocating memory for Poisson point process data");
    if (to->mpp.moments != from->mpp.moments) BUG;
    //printf("size = %d\n", sizeof(mpp_properties));
    assert(sizeof(mpp_properties) == 96);
    mpp_properties *To = &(to->mpp), *From=&(from->mpp);
    assert(To != NULL &&  From != NULL);
    //    To->sum_zhou_c = From->sum_zhou_c;
    //    To->sq_zhou_c = From->sq_zhou_c;
    int vdim = from->vdim[0];
    for (v=0; v<vdim; v++) To->maxheights[v] = From->maxheights[v];
    To->unnormedmass = From->unnormedmass;
    //    To->zhou_c = From->zhou_c;    
    assert(To->mM != NULL && To->mMplus != NULL);
    int nmP1 = To->moments + 1;
    MEMCOPY(To->mM, From->mM, nmP1 * sizeof(double));
    MEMCOPY(To->mMplus, From->mMplus, nmP1 * sizeof(double));

    if (to->qlen != from->qlen) BUG;
    if (from->qlen > 0) {
      assert(to->q != NULL);
      MEMCOPY(to->q, from->q, (to->qlen)* sizeof(double));
    }
  }

  if (recursive) {
    for (i=0; i<MAXSUB; i++) if (from->sub[i] != NULL) {
	assert(to->sub[i] != NULL);
	paramcpy(to->sub[i], from->sub[i], freeing,
		 allocating, copy_lists, recursive, copy_mpp);
      }
  }
  //PMI(to, "to");
}


int covcpy(cov_model **localcov, bool sub, cov_model *cov, // err
	   location_type *prevloc, location_type *ownloc,
	   bool copy_lists, bool copy_randomparam) {
  return covcpy(localcov, sub, cov, prevloc, ownloc, copy_lists, // err
		 copy_randomparam, false);
}

int covcpy(cov_model **localcov, bool sub, cov_model *cov, // err
	   location_type *prevloc, location_type *ownloc,
	   bool copy_lists, bool copy_randomparam, 
	   bool allowCopyingInterface) {
  assert(cov != NULL);
  int i,
    n = -1;
  //cov_fct *C = CovList + cov->nr; // nicht gatternr

  //PMI(cov);
 
  if ((*localcov = (cov_model*) MALLOC(sizeof(cov_model)))==0) 
    return ERRORMEMORYALLOCATION;
  cov_model *current = *localcov;

  MEMCOPY(current, cov, sizeof(cov_model)); // replaces COV_NULL(*localcov);
  COV_ALWAYS_NULL(current);
  current->calling = NULL; 
  paramcpy(current, cov, false, true, copy_lists, false, false);

  if (cov->ownkappanames != NULL) {
    int nkappas = CovList[cov->nr].kappas;
    current->ownkappanames = (char**)  CALLOC(nkappas, sizeof(char*));
    for (i=0; i<nkappas; i++) {
      if (cov->ownkappanames[i] != NULL) {
	current->ownkappanames[i] =
	  (char*) MALLOC(sizeof(char) * (1 + strlen(cov->ownkappanames[i])));
	strcpy(current->ownkappanames[i], cov->ownkappanames[i]);
      }
    }
  }

  //  PMI(current, "err = covcpy");
  
  if (cov->q != NULL) {    
    n = sizeof(double) * current->qlen;
    current->q = (double*) MALLOC(n); // QALLOC NOT APPROPRIATE
    MEMCOPY(current->q, cov->q, n);
  } else assert(current->qlen==0);

  current->prevloc = ownloc != NULL ? ownloc : 
    cov->prevloc == prevloc ? prevloc : NULL;
  if (current->prevloc == cov->prevloc && cov->calling==NULL) {
    if (!isInterface(cov)) {
      BUG;
    }
    if (!allowCopyingInterface) {PRINTF("\n\n***** unallowed copying ******\n"); BUG;}
  }
  
  for (i=0; i<MAXPARAM; i++) {
    int err;
    current->kappasub[i] = NULL;
    if (cov->kappasub[i] == NULL || !copy_randomparam) continue;
    err = covcpy(current->kappasub + i, true, cov->kappasub[i], 
		 prevloc, ownloc, copy_lists, copy_randomparam);
    if (err != NOERROR) return err;
    current->kappasub[i]->calling = current;       
  }
 
  if (sub) {
    for (i=0; i<MAXSUB; i++) {
      int err;
      current->sub[i] = NULL;
      if (cov->sub[i] == NULL) continue;
      err = covcpy(current->sub + i, sub, cov->sub[i], prevloc, ownloc,
		   copy_lists, copy_randomparam);
      if (err != NOERROR) return err;
      current->sub[i]->calling = current;       
    }
  } else {
    for (i=0; i<MAXSUB; i++) current->sub[i] = NULL;
  }

  // PMI(*localcov, "end err=covcpy");
  return NOERROR;
}


int covcpy(cov_model **localcov, bool sub, cov_model *cov, // err
	   location_type *prevloc, location_type *ownloc,
	   bool copy_lists) {
  
  int 
    err =covcpy(localcov, sub, cov, prevloc, ownloc, copy_lists, true); //err
  return err;
}
 
int covcpy(cov_model **localcov, cov_model *cov) { //err
  bool cov2key = &(cov->key)==localcov;
  //  PMI(cov); 
  int 
    err = covcpy(localcov, true, cov, cov->prevloc, NULL, true, true);//err
  //  PMI(*localcov);
  if (err == NOERROR) 
    (*localcov)->calling = cov2key || cov->calling==NULL ? cov : cov->calling;
  // falls !cov2key && cov->calling == NULL; dann gibt es nur einen Verweis
  // rueckwaerts! Muss aber sein, da sonst LOC_DELETE bei cov->calling==NULL
  // versucht cov->prevloc zu loeschen. oder es muesste prevloc auf NULL
  // gesetzt zu werden. Dies scheint eher contraproduktiv zu sein.
  return err;
}

int covcpyWithoutRandomParam(cov_model **localcov, cov_model *cov) {//err
  bool cov2key = &(cov->key)==localcov;
  int 
    err = covcpy(localcov, true, cov, cov->prevloc, NULL, true, false);
  if (err == NOERROR) 
    (*localcov)->calling = cov2key || cov->calling==NULL ? cov : cov->calling;
  // falls !cov2key && cov->calling == NULL; dann gibt es nur einen Verweis
  // rueckwaerts! Muss aber sein, da sonst LOC_DELETE bei cov->calling==NULL
  // versucht cov->prevloc zu loeschen. oder es muesste prevloc auf NULL
  // gesetzt zu werden. Dies scheint eher contraproduktiv zu sein.
  return err;
}


int covcpy(cov_model **localcov, cov_model *cov, bool copy_lists) {//err
  bool cov2key = &(cov->key)==localcov;
  int 
    err = covcpy(localcov, true, cov, cov->prevloc, NULL, copy_lists, true);
  if (err == NOERROR) 
    (*localcov)->calling = cov2key || cov->calling==NULL ? cov : cov->calling;
  // falls !cov2key && cov->calling == NULL; dann gibt es nur einen Verweis
  // rueckwaerts! Muss aber sein, da sonst LOC_DELETE bei cov->calling==NULL
  // versucht cov->prevloc zu loeschen. oder es muesste prevloc auf NULL
  // gesetzt zu werden. Dies scheint eher contraproduktiv zu sein.
  return err;
}

int covcpy(cov_model **localcov, cov_model *cov, //err
	   double *x, double *T, int spatialdim, int xdimOZ, long lx, bool Time, 
	   bool grid, bool distances) {
  bool cov2key = &(cov->key)==localcov;
  int err;
  location_type *loc= NULL;

  err = loc_set(x, T, spatialdim, xdimOZ, lx, Time, grid, distances, &loc); 
  if (err != NOERROR) return err;
  if ((err = covcpy(localcov, true, cov, loc, NULL, true, true)) != NOERROR)
    return err;
  (*localcov)->prevloc = cov->prevloc;
  (*localcov)->ownloc=loc;
  (*localcov)->calling = cov2key || cov->calling==NULL ? cov : cov->calling;
  (*localcov)->ownloc->totalpoints = loc->totalpoints;
  return NOERROR;
}


cov_model *getRemote(cov_model *remotecov, cov_model *rmt, cov_model *target) {
  cov_model *found;
  int i;
  if (rmt == target) return remotecov;
  
  for (i=0; i<MAXPARAM; i++) {
    if (rmt->kappasub[i] != NULL) {
      if (remotecov->kappasub[i] == NULL) BUG;
      found = getRemote(remotecov->kappasub[i], rmt->kappasub[i], target);
      if (found != NULL) return found;
    }
  }
 
  for (i=0; i<MAXSUB; i++) {
    if (rmt->sub[i] != NULL) {
      if (remotecov->sub[i] == NULL) BUG;
      found = getRemote(remotecov->sub[i], rmt->sub[i], target);
      if (found != NULL) return found;
    }
  }
  return NULL;  
}


void Ssetcpy(cov_model *localcov, cov_model *remotecov, cov_model *cov,
	    cov_model *rmt) {
  int i;
  if (cov->Sset != NULL) {
    localcov->Sset = (set_storage*) MALLOC(sizeof(set_storage));
    MEMCOPY(localcov->Sset, cov->Sset, sizeof(set_storage));
    localcov->Sset->remote = getRemote(remotecov, rmt, cov->Sset->remote);
    if (localcov->Sset->remote == NULL) BUG;
  }
  
  for (i=0; i<MAXPARAM; i++) {
    if (cov->kappasub[i] != NULL) {
      if (localcov->kappasub[i] == NULL) BUG;      
      Ssetcpy(localcov->kappasub[i], remotecov, cov->kappasub[i], rmt);
    }
  }
 
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      if (localcov->sub[i] == NULL) BUG;      
      Ssetcpy(localcov->sub[i], remotecov, cov->sub[i], rmt);
    }
  }
}



int newmodel_covcpy(cov_model **localcov, int model, cov_model *cov, //err
		    double *x, double *y, double *T, 
	            int spatialdim, /* spatial dim only ! */
	            int xdimOZ, long lx, long ly, bool Time, bool grid,
	            bool distances) {
  int i, err;
  assert(*localcov == NULL);
  addModel(localcov, model, NULL, true);
  cov_model *neu = *localcov;
  Types type;
  
  loc_set(x, y, T, spatialdim, xdimOZ, lx, ly, Time, grid, distances,
	  &(neu->ownloc));
  if ((err = covcpy(neu->sub + 0, cov)) != NOERROR) return err;

  // PMI(neu, "hier");

  neu->sub[0]->calling = neu;

  //   PMI(neu);
  //printf("newmodel %s %s\n", NICK(cov), TYPENAMES[CovList[cov->nr].Type]);

  type = CovList[neu->nr].Typi[0];
  assert(type == InterfaceType && CovList[neu->nr].variants == 1); // anderes z.zt nicht benutzt
  for (i=0; i<2; i++) {
    //print("i=%d %s\n", i, TYPENAMES[cov->typus]);
    if ((err = CHECK(neu, cov->tsdim, cov->xdimprev,
		     // CovList[cov->nr].Type,Martin:changed 27.11.13 
		     //cov->typus,// auch nicht 20.5.14
		     type,
		     type == InterfaceType ? XONLY : cov->domprev, 
		     type == InterfaceType ? CARTESIAN_COORD : cov->isoprev, 
		     cov->vdim, ROLE_BASE)) != NOERROR) {
      return err;
    }
    if (i==0 && (err =  STRUCT(neu, NULL)) != NOERROR) return err;
  }
  return NOERROR;
}

int newmodel_covcpy(cov_model **localcov, int model, cov_model *cov) {//err
  
  int err;
  location_type *loc = Loc(cov);
  
  err = newmodel_covcpy(localcov, model, cov,
               loc->grid ? loc->xgr[0] : loc->x, 
	       loc->grid ? loc->ygr[0] : loc->y, 
               loc->grid ? loc->xgr[0] + 3 * loc->spatialdim : loc->T, 
	       loc->spatialdim, loc->xdimOZ, 
	       loc->grid ? 3 : loc->totalpoints, 
	       loc->ly == 0 ? 0 : loc->grid ? 3 : loc->totalpoints,
	       loc->Time, loc->grid, loc->distances);
  return err;
  
}


double *getAnisoMatrix(cov_model *cov, bool null_if_id, int *nrow, int *ncol) {
  int 
    origdim = cov->prevloc->timespacedim;
  if (!isAnyDollar(cov) && null_if_id) {
    *nrow = *ncol = origdim;
    return NULL;
  }

  double *ani,
    *aniso = P(DANISO),
    a = PisNULL(DSCALE) ? 1.0 : 1.0 / P0(DSCALE);   
  int i, total,
    *proj = PINT(DPROJ),
    dimP1 = origdim + 1;

  if (aniso != NULL) {
    total = origdim * cov->ncol[DANISO];
    int bytes = total * sizeof(double);
    ani = (double *) MALLOC(bytes);
    MEMCOPY(ani, aniso, bytes); 
    for (i=0; i<total; i++) {
      //     printf("anui %d %f %f\n", i, ani[i], a);
      ani[i] *= a;    
    }
    *nrow = cov->nrow[DANISO];
    *ncol = cov->ncol[DANISO];

    // printf("A %d %d\n", *nrow, *ncol);

  } else if (proj != NULL) {
    total = origdim * cov->nrow[DPROJ];
    //
    //    print("origdim %d %f %d\n", origdim, a, proj[0]);
    ani = (double *) CALLOC(total, sizeof(double));
    for (i=0; i<total; i+=dimP1) ani[i * origdim + proj[i] - 1] = a; 
    *nrow = origdim;
    *ncol = cov->nrow[DPROJ];
    //  printf("BA %d %d\n", *nrow, *ncol); assert(false);
  } else {
    if (a == 1.0 && null_if_id) {
      *nrow = *ncol = origdim;
      return NULL;
    }
    total = origdim * origdim;
    ani = (double *) CALLOC(total, sizeof(double));
    for (i=0; i<total; i+=dimP1) ani[i] = a;
    *nrow = *ncol = origdim;
    //    printf("CA %d %d\n", *nrow, *ncol);
 }
  //   printf("DA %d %d\n", *nrow, *ncol);

  return ani;
}

double *getAnisoMatrix(cov_model *cov, int *nrow, int *ncol) {
  return getAnisoMatrix(cov, false, nrow, ncol);
}


double GetDiameter(location_type *loc,
		   double *min, double *max, double *center) {
  
  // calculates twice the distance between origcenter %*% aniso and
  // all 2^spatialdim corners spanned by min and max; returns maximum distance.

  // NOTE: origcenter is not alway (min + max)/2 since origcenter might be
  //       given by the user and not calculated automatically !

   bool *j=NULL;
  int d,
    origdim = loc->timespacedim,  
    spatialdim = loc->spatialdim;
  double distsq, dummy,
    *lx=NULL, *sx=NULL, 
   diameter=0.0;
 
  if (loc->grid) {
    double 
      *origmin = (double*) MALLOC(origdim * sizeof(double)),
      *origmax = (double*) MALLOC(origdim * sizeof(double)),
      *origcenter = (double*) MALLOC(origdim * sizeof(double));

    for (d=0; d<origdim; d++) {   
      if (loc->xgr[d][XSTEP] > 0) {
	origmin[d] = loc->xgr[d][XSTART];
	origmax[d] = loc->xgr[d][XSTART] + loc->xgr[d][XSTEP] * 
	  (loc->xgr[d][XLENGTH] - 1.0);
      } else {
	origmin[d] = loc->xgr[d][XSTART] + loc->xgr[d][XSTEP] * 
	  (double) (loc->xgr[d][XLENGTH] - 1.0);
	origmax[d] = loc->xgr[d][XSTART];
      }
      origcenter[d] = 0.5 * (origmin[d] + origmax[d]);
    }

    if (loc->caniso == NULL) {
      for  (d=0; d<origdim; d++) {   
	center[d] = origcenter[d];
	min[d] = origmin[d];
	max[d] = origmax[d];
	dummy = max[d] - min[d];
	diameter += dummy * dummy;
      }
    } else { // caniso != NULL
      j = (bool*) MALLOC( (origdim + 1) * sizeof(double));
      lx = (double*) MALLOC( origdim * sizeof(double));
      sx = (double*) MALLOC(spatialdim * sizeof(double));
      
      xA(origcenter, loc->caniso, origdim, spatialdim, center);
      for (d=0; d<origdim; d++) {
	j[d]=false;
	lx[d]=origmin[d];
      }
      j[origdim] = false;
      
      for (d=0; d<spatialdim; d++) {
	min[d] = RF_INF;
	max[d] = RF_NEGINF;
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
	xA(lx, loc->caniso, origdim, spatialdim, sx);
	
	// suche maximale Distanz von center zu transformierter Ecke
	distsq = 0.0;
	for (d=0; d<spatialdim; d++) {
	  if (min[d] > sx[d]) min[d] = sx[d];
	  if (max[d] < sx[d]) max[d] = sx[d];
	  dummy = center[d] - sx[d];
	  distsq += dummy * dummy;	  
	  //  print("%d sx=%f %f %f diam=%f\n", d, sx[d], min[d], max[d], distsq);	  
	}
	if (distsq > diameter) diameter = distsq;
      }

      UNCONDFREE(j);
      UNCONDFREE(lx);
      UNCONDFREE(sx);
    }
  
    UNCONDFREE(origmin);
    UNCONDFREE(origmax);
    UNCONDFREE(origcenter);

  } else { // not loc->grid

    if (loc->caniso != NULL) BUG;

    double *xx=loc->x; 
    int i,
      endfor = loc->spatialtotalpoints * spatialdim;

    for (d=0; d<spatialdim; d++) {
      min[d]=RF_INF; 
      max[d]=RF_NEGINF; 
    }

    // to determine the diameter of the grid determine as approxmation
    // componentwise min and max corner
    for (i=0; i<endfor; ) {
      for (d=0; d<spatialdim; d++, i++) {
        //temporal part need not be considered, but for ease included#
	//	printf("i=%d %d %d\n", i, d, endfor);
	if (xx[i] < min[d]) min[d] = xx[i];
	if (xx[i] > max[d]) max[d] = xx[i];
      }
    }
       
    if (loc->Time) {
      assert(d == origdim - 1);
      if (loc->T[XSTEP] > 0) {
	min[d] = loc->T[XSTART];
	max[d] = loc->T[XSTART] + loc->T[XSTEP] * (loc->T[XLENGTH] - 1.0);
      } else {
	min[d] = loc->T[XSTART] + loc->T[XSTEP] * (loc->T[XLENGTH] - 1.0);
	max[d] = loc->T[XSTART];
      }
    }
    
    for (diameter=0.0, d=0; d<origdim; d++) {
      center[d] = 0.5 * (max[d] + min[d]); 
      dummy = max[d] - min[d];
      diameter += dummy * dummy;
      // print("%d center %f %f %f %f\n", d, center[d], min[d], max[d], diameter);
    }
  } // !simugrid
    
  return 2.0 * sqrt(diameter);
}

double GetDiameter(location_type *loc) { 
  double diam,
    *dummymin=NULL, *dummymax=NULL, *dummycenter=NULL;
  int origdim = loc->timespacedim; 
  
  dummymin = (double*) MALLOC(origdim * sizeof(double));
  dummymax = (double*) MALLOC(origdim * sizeof(double));
  dummycenter = (double*) MALLOC(origdim * sizeof(double));
  diam = GetDiameter(loc, dummymin, dummymax, dummycenter);
  UNCONDFREE(dummymin);
  UNCONDFREE(dummymax);
  UNCONDFREE(dummycenter);
  return diam;
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
  unsigned long i,ii,j,jj,l,ll,min, m=1;
  if (HOMEMADE_NICEFFT) {
    int f[F_NUMBERS1]={2,3,5}; 
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


void expandgrid(coord_type xgr, double **xx, int nrow){
  double *x=NULL, *y=NULL; /* current point within grid, but without
		       anisotropy transformation */
  int 
    dimM1 = nrow - 1,
    *yi=NULL, /* counter for the current position in the grid */
    ncol = nrow;

  long d, pts, w, k, total, i; 
  for (pts=1, i=0; i<nrow; i++) {
    pts *= (long int) xgr[i][XLENGTH];
  }

  y = (double*) MALLOC(ncol * sizeof(double));
  yi = (int*) MALLOC(ncol * sizeof(int));

  total = ncol * pts;

  x = *xx = (double*) MALLOC(sizeof(double) * total);

  for (w=0; w<ncol; w++) {y[w]=xgr[w][XSTART]; yi[w]=0;}
  for (k=0; k<total; ){
    for (d=0; d<ncol; d++, k++) {
      x[k] = y[d];
    }
    i = 0;
    (yi[i])++;
    y[i] += xgr[i][XSTEP];
    while(yi[i]>= xgr[i][XLENGTH]) {
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
  UNCONDFREE(y);
  UNCONDFREE(yi);
}



void expandgrid(coord_type xgr, double **xx, double* aniso, int nrow, int ncol){
  double *x=NULL, * y=NULL; /* current point within grid, but without
		       anisotropy transformation */
  int
    *yi=NULL,   /* counter for the current position in the grid */
    dimM1 = nrow - 1; 
  long pts, w, k, total, n, i, d;

  if (aniso == NULL) {
    expandgrid(xgr, xx, nrow);
    return;
  }

  for (pts=1, i=0; i<nrow; i++) {
    pts *= (long int) xgr[i][XLENGTH];
  }

  total = ncol * pts;
  x = *xx = (double*) MALLOC(sizeof(double) * total);
  y = (double*) MALLOC(ncol * sizeof(double));
  yi = (int*) MALLOC(ncol * sizeof(int));


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
    while(yi[i]>=xgr[i][XLENGTH]) {
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
  UNCONDFREE(y);
  UNCONDFREE(yi);
}



void grid2grid(coord_type xgr, double **grani, double *aniso, int origdim,
	       int dim) {
  // diag is assumed
  double *pgr, *A;
  int d, i,
    origdimM1 = origdim - 1;

  (*grani) = (double *) MALLOC(sizeof(double) * 3 * dim);
  pgr = *grani;

  if (aniso == NULL) {
    for(d=0; d<dim; d++) {
      for (i=0; i<3; i++, pgr++) {
        *pgr = xgr[d][i];
      }
    }
  } else {

    //printf("%d %d %f %f\n", origdim, dim, aniso[0], aniso[1]);

    for (d=0; d<dim; d++, pgr += 3) {
      A = aniso + d * origdim;
      for (i=0; i<origdimM1; i++, A++) if (*A != 0.0) break;
      //printf("%d %f %d %d \n", d, *A, i, origdimM1);
      // factor = d < cov->nrow ? aniso[n] : 0.0; // hier zwingend diag !!
      
      pgr[XSTART] = xgr[i][XSTART] * *A;
      pgr[XSTEP] = xgr[i][XSTEP] * *A;
      pgr[XLENGTH] = xgr[i][XLENGTH];
    }
  }
}


void xtime2x(double *x, int nx, double *T,
	     double **newx, int timespacedim) {
  double *y, // umhaengen zwingend notwendig, da u.U. **xx und *newx
    // der gleiche aufrufende Pointer ist, d.h. es wird "ueberschrieben"
    *z,  t;
  int j, k, i, d, 
    timelen = (int) T[XLENGTH],
    spatialdim = timespacedim - 1;
  z = *newx = (double*) MALLOC(sizeof(double) * timespacedim * nx * timelen);
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


void xtime2x(double *x, int nx, double *T,
	     double **newx, double *aniso, int nrow, int ncol) {
  double *y, // umhaengen zwingend notwendig, da u.U. **xx und *newx
    // der gleiche aufrufende Pointer ist, d.h. es wird "ueberschrieben"
    *z, dummy, t;
  int j, k, i, d, n, endfor, w,
    spatialdim = nrow - 1,
    timelen = T[XLENGTH],
    nxspdim = nx * spatialdim;

  if (aniso == NULL) {
    assert(nrow == ncol);
    xtime2x(x, nx, T, newx, nrow);
    return;
  }

  z = *newx = (double*) MALLOC(sizeof(double) * ncol * nx * timelen); 
 
  for (k=j=0, t=T[XSTART]; j<timelen; j++, t += T[XSTEP]){
    y = x;
    for (i=0; i<nxspdim; i+=spatialdim) {
      endfor = i + spatialdim;
      for (n=d=0; d<ncol; d++) {
	dummy = 0.0;
	for(w=i; w<endfor; w++) {
	  dummy += aniso[n++] * y[w];
	}
	//	print("%d %d %f\n", k, n, t);
	z[k++] = dummy + aniso[n++] * t; //auf keinen Fall nach vorne holen
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

  z = *newx = (double*) MALLOC(sizeof(double) * ncol * nx);  
  if (aniso == NULL) {
    assert(nrow == ncol);
    MEMCOPY(z, x, sizeof(double) * nx * ncol);
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
	*m0 = (double*) MALLOC(sizeof(double) * dim1 * dim3);
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




bool isMiso(matrix_type type) { 
  return type == TypeMiso;
}


bool isMproj(matrix_type type) { 
  assert(TypeMtimesepproj == 2);
 return type == TypeMproj || type <= TypeMtimesepproj;
}


bool isMdiag(matrix_type type) { 
  return type <= TypeMdiag;
}

bool isMtimesep(matrix_type type) { 
  assert(TypeMtimesep == 3);
  return type <= TypeMtimesep;
}



matrix_type Type(double *M, int nrow, int ncol) {

    // see also analyse_matrix: can it be unified ???

  matrix_type type = TypeMiso; // default
  //  double elmt;
  int i, j, k,
    endfor = nrow * ncol;

  double *m = M;

  if (m==NULL || (ncol == 1 && nrow == 1)) {
    assert(ncol == nrow);
    return type;
  }
   
  if (ncol > nrow) {
    for (i=ncol * ncol; i < endfor; i++) 
      if (m[i]!= 0.0) return TypeMany;
    ncol = nrow; // !!
  }

  //  for(i=0; i<ncol * nrow; i++) printf("%f ", M[i]);
  //  print("\n%d %d\n", ncol, nrow);

  for (k=0; k<ncol; ) {
    for (i=0; i<nrow; i++) if (!R_FINITE(m[i]) || m[i] != 0.0) break;
    for (j=i+1; j<nrow; j++) 
      if (!R_FINITE(m[j]) || m[j] != 0.0) goto TimeInvestigation;
    if (k == ncol - 1) { 
      /* // nur fuer strikte projektion auf kleinere Koordinaten
	 if (type == TypeMproj) {
	int ii, kk;
	for (ii=0; ii<nrow; ii++) {
	  double sum = 0.0;
	  for (kk=ii; kk<endfor; kk+= nrow) sum += M[kk];
	  printf("sum %d %f %d\n", ii, sum, type); //
	  if (sum != 1.0 && sum != 0.0) {
	    type = TypeMany;
	    break;
	  }
	}
      }
      */
      if (i != k && (!R_FINITE(m[i]) || m[i] != 0.0))
	return type == TypeMproj && R_FINITE(m[i]) && m[i] == 1.0 
	  ? TypeMproj : TypeMany;
      else return type == TypeMproj 
	     ? (R_FINITE(m[i]) && (m[i] == 0.0 || m[i] == 1.0)
		? TypeMtimesepproj : TypeMany)
	     : type == TypeMiso && (!R_FINITE(m[i]) || m[i] != M[0])
	     ? TypeMdiag : type;
    }

    // k < ncol - 1; type != timesep and != timesepproj garantied
    matrix_type newtype;
    if ( i!=k && m[i] != 0.0) {
      if (m[i] != 1.0) goto TimeInvestigation; 
      newtype = TypeMproj;
    } else newtype = R_FINITE(m[i]) && m[i]==M[0] ? TypeMiso : TypeMdiag;
    if (newtype > type) type = newtype;

    k++;
    m += nrow;
    continue;

  TimeInvestigation: 
    if (k == ncol - 1) return TypeMany;
    type = TypeMtimesep;
    k = ncol - 1;
    m = M + k * nrow;
  }
  return type;
}


 
void Transform2NoGridExt(cov_model *cov, bool timesep, int gridexpand, 
			 double **grani, double **SpaceTime, 
			 double **caniso, int *Nrow, int *Ncol,
			 bool *Time, bool *grid, int *newdim, bool takeX) {
  // this function transforms the coordinates according to the anisotropy
  // matrix 

  if (!isCartesian(cov->isoown)) 
    ERR("transformation only possible within Cartesian system");
  
  location_type *loc = Loc(cov);
  bool isdollar = isAnyDollar(cov);

  matrix_type type;
  int 
    nrow = -1,
    ncol = -1,
    origdim = loc->caniso == NULL ? loc->timespacedim : loc->cani_nrow,
    dim = (isdollar ? (!PisNULL(DANISO) ? cov->ncol[DANISO] : 
		       !PisNULL(DPROJ) ? cov->nrow[DPROJ] : origdim )
	   : origdim
	   );
  double *aniso,
    *T = loc->T, *x = takeX ? loc->x : loc->y;
  coord_type 
    *xgr= takeX ? &(loc->xgr) : &(loc->ygr);
  
  if (x==NULL && (*xgr)[0] ==NULL) ERR("locations are all NULL");
 
  *newdim = dim;
  *caniso = NULL;
  *Nrow = *Ncol = -1;

  // print("newdim dim %d %d\n", *newdim, dim);

  aniso = getAnisoMatrix(cov, true, &nrow, &ncol); // null if not dollar

  

  // 
  //printf("bb NCOL %d %d %ld loc->cani_nrow %d %d aniso=%f %f \n", 
  //          nrow, ncol, aniso, loc->cani_nrow, loc->cani_ncol, aniso[0], aniso[1]);
  // assert(false);

  if (loc->caniso != NULL) {    
    if (aniso == NULL) {
      unsigned long bytes = sizeof(double) * loc->cani_nrow * loc->cani_ncol;
      aniso =(double*) MALLOC(bytes);
      MEMCOPY(aniso, loc->caniso, bytes);
      nrow = loc->cani_nrow;
      ncol = loc->cani_ncol;
    } else {
      double *aniso_old = aniso;
      assert(loc->cani_ncol == nrow);
      aniso = matrixmult(loc->caniso, aniso_old, loc->cani_nrow, nrow, ncol);
      nrow = loc->cani_nrow;
      UNCONDFREE(aniso_old);
    }
  }

  //print("transfor\n");

  type = aniso == NULL ? TypeMiso : Type(aniso, origdim, dim);
  *Time = loc->Time;
  *grid = loc->grid && !gridexpand;

  //print("grid %d %d %d\n", *grid, loc->grid, gridexpand);
  
  //  printf("nrow %d origdim=%d dim=%d %d type=%d grid=%d type=%d\n",
  //	 nrow, origdim, dim, ncol, type, loc->grid, isMproj(type)); 
  //   assert(false);
  //PMI(cov); if (origdim != nrow) crash(cov);
  assert(origdim == nrow);
  assert(dim == ncol);
  if (loc->grid) {
    assert(*xgr != NULL);
    if (isMproj(type)) {
      if (gridexpand == true) {
	expandgrid(*xgr, SpaceTime, aniso, nrow, ncol);
	*Time = false;
      } else {
	// nur aniso auf grid multipliziert
	grid2grid(*xgr, grani, aniso, nrow, ncol);
	*grid = true;	
      }
    } else if (gridexpand) {
      if (!loc->Time) { // grid no time
	expandgrid(*xgr, SpaceTime, aniso, nrow, ncol);
      } else { // grid and time 
	if (timesep && isMtimesep(type)) {
	  // space
	  expandgrid(*xgr, SpaceTime, aniso, nrow, ncol - 1); 
	  // time
	  grid2grid(*xgr + loc->spatialdim, grani, 
		    aniso + nrow * nrow - 1, 1, 1);
	} else {
	  expandgrid(*xgr, SpaceTime, aniso, nrow, ncol);
	  *Time = false;
	}
      }
    } else {     
      //      assert(false);
      //      crash();
      int i,d,k;
      (*grani) = (double *) MALLOC(sizeof(double) * 3 * origdim);
      for (k=d=0; d<origdim; d++) {
	for (i=0; i<3; i++) {
	  //   printf(" !!!!!!! xgr %d %d %d %f\n", takeX, d, i,(*xgr)[d][i]); 
	  (*grani)[k++] = (*xgr)[d][i];
	}
      }
      *newdim = dim;     
      *caniso = aniso;
      *Nrow = nrow;
      *Ncol = ncol;
      return; // hier rausgehen!
    }
  } else { // nogrid
    if (!loc->Time) { // no grid no time
	 x2x(x, loc->spatialtotalpoints, SpaceTime, aniso, nrow, ncol); 
    } else { // no grid but time
      if (timesep && isMtimesep(type)) {
	x2x(x, loc->spatialtotalpoints, SpaceTime, aniso, nrow, ncol-1);
	grid2grid((*xgr) + loc->spatialdim, grani, 
		  aniso + nrow * nrow - 1, 1, 1);
      } else {
	xtime2x(x, loc->spatialtotalpoints, T, SpaceTime,
		aniso, nrow, ncol);
	*Time = false;
      }
    }
  }

  FREE(aniso);
  
  // printf("NCOL %d\n", *ncol);
}


void Transform2NoGrid(cov_model *cov, bool timesep, int gridexpand) {
  location_type *loc = cov->prevloc; // transform2nogrid
  assert(loc != NULL);
  bool Time, grid;
  int err,
    spacedim=-1, 
    nrow=-1,
    ncol=-1;
  double  *xgr=NULL, 
    *x = NULL,
    *caniso = NULL;

  if ((loc->y != NULL && loc->y != loc->x) || 
      (loc->ygr[0] != NULL && loc->ygr[0] != loc->xgr[0])) {
    //printf("ly=%d %ld %ld %ld\n", loc->ly, loc->y, loc->ygr, loc->xgr);
    //    PMI(cov->calling->calling);
    ERR("unexpected y coordinates");
  }

  // APMI(cov->calling);
  //print("spacedim %d time=%d %f\n", spacedim, Time, (double) loc->totalpoints);

  assert(cov->prevloc != NULL);
  Transform2NoGridExt(cov, timesep, gridexpand, &xgr, &x, 
		      &caniso, &nrow, &ncol, &Time, &grid, &spacedim, true);
  
  //   print("spacedim %d time=%d %ld %d\n", spacedim, Time, loc->totalpoints, grid);  assert(false);
  //APMI(cov->calling);

  //  PMI(cov, "trafo");

  if (Time) spacedim--;
  assert(cov->ownloc == NULL);
  err = loc_set(grid ? xgr : x, 
		grid ? xgr + 3 * spacedim : xgr, 
		spacedim, 
		spacedim,
		grid ? 3 : loc->totalpoints, 
		Time, grid, false, &(cov->ownloc));
  //print("AA spacedim %d time=%d\n", spacedim, Time);


  assert(cov->ownloc->caniso == NULL);
  cov->ownloc->caniso = caniso;
  cov->ownloc->cani_nrow = nrow;
  cov->ownloc->cani_ncol = ncol;

  //  printf("nrow %d %d %ld %ld\n", nrow, ncol, caniso, cov->ownloc);
  //  PMI(cov, "transform to no grid");
  //nrow 0 256501392 256504320 256504528
  // assert(ncol < 10);

  caniso = NULL;

  FREE(x);
  FREE(xgr);
  if (err != NOERROR) ERR("error when transforming to no grid");
}


void Transform2NoGrid(cov_model *cov, double **xx) {
  bool Time, grid;
  int newdim, nrow, ncol;
  double *caniso = NULL;
  Transform2NoGridExt(cov, false, true, NULL, xx, 
		      &caniso, &nrow, &ncol, &Time, &grid, &newdim, true); 
  assert(caniso == NULL);
}


void Transform2NoGrid(cov_model *cov, double **xx, double **yy) {
  location_type *loc = Loc(cov);
  bool Time, grid;
  int newdim, nrow, ncol;
  double *caniso = NULL;
  Transform2NoGridExt(cov, false, true, NULL, xx,
		      &caniso,&nrow, &ncol, &Time, &grid, &newdim, true); 
  assert(caniso == NULL);
  if (loc->y == NULL && loc->ygr[0] == NULL) *yy = NULL;
  else Transform2NoGridExt(cov, false, true, NULL, yy, 
			   &caniso, &nrow, &ncol, &Time, &grid, &newdim, false);
  assert(caniso == NULL);
}



double *selectlines(double *m, int *sel, int nsel, int nrow, int ncol) {
    // selects lines in matrix m according to sel
  int j;  
  double *red_matrix,
      *Red_Matrix = (double *) MALLOC(sizeof(double) * nsel * ncol),
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
      *Red_Matrix = (int *) MALLOC(sizeof(int) * nsel * ncol),
      *endfor = Red_Matrix + nsel * ncol;
  
  for (red_matrix=Red_Matrix; red_matrix<endfor; m+=nrow) {
      for (j=0; j<nsel; j++, red_matrix++) {
	  *red_matrix = m[sel[j]];
      }
  }
  return Red_Matrix;
}

bool TypeConsistency(Types requiredtype, Types deliveredtype) {
  if (deliveredtype == UndefinedType) BUG;
  switch(requiredtype) {
  case TcfType:        return isTcf(deliveredtype);
  case PosDefType:     return isPosDef(deliveredtype);
  case VariogramType:  return isVariogram(deliveredtype);
  case NegDefType:     return isNegDef(deliveredtype);
  case ProcessType :   return (isProcess(deliveredtype) ||
			       isTrend(deliveredtype));
  case GaussMethodType:return isGaussMethod(deliveredtype);
  case BrMethodType:   return isBRuserProcess(deliveredtype);
  case PointShapeType: return isPointShape(deliveredtype);
  case RandomType :    return isRandom(deliveredtype);
  case ShapeType :     return isShape(deliveredtype);
  case TrendType :     return isTrend(deliveredtype);
  case InterfaceType:  return isInterface(deliveredtype);
  case OtherType :     return isOther(deliveredtype);
  default : BUG;
  }
  BUG;
  return FALSE;
}

int TypeConsistency(Types requiredtype, cov_model *cov, int depth) {
  //  if (cov == NULL) crash();
  assert(cov != NULL);
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  // printf("entering TC %s requiredtype=%s undefined!\n", NAME(cov), TYPENAMES[requiredtype]);
  if (isUndefined(C)) {
    assert(C->TypeFct != NULL);

    // printf(" undefined!\n");
    int tc = C->TypeFct(requiredtype, cov, depth);
    //  printf("%d\n", tc);
    return tc;
  }  
  //   printf("not consist %s type=%d, %d;   %s %d\n", NAME(cov), C->Typi[0], cov->typus,	 CovList[150].name, CovList[150].Typi[0]);
  int i;
  for (i=0; i<C->variants; i++) {
    //    printf("TC %d %s: %d %s %s %d %d \n", i, NAME(cov), TypeConsistency(requiredtype, C->Typi[i]), TYPENAMES[requiredtype], TYPENAMES[C->Typi[i]],  depth, depth<=0 );
    if (TypeConsistency(requiredtype, C->Typi[i])
	&& (depth <= 0 || atleastSpecialised(cov->isoown, C->Isotropy[i]))
	) {
      //         printf("ok ? %d \n", i +1);
      return i + 1;
    }
 }

  //  printf("not ok\n");

  return false;
}


int change_coordinate_system(isotropy_type callingisoown,
			     isotropy_type isoprev,
			     int tsdimprev, int xdimprev,
			     int *nr, isotropy_type  *newisoprev,
			     int *newtsdim, int *newxdim) {
  // isoprev is replaced by a coordinate transformation and a nwe 

  // printf("ISO calling %s; prev %s\n", ISONAMES[callingisoown], ISONAMES[isoprev]);

  switch(callingisoown) {
  case EARTH_COORD :
    if (isCartesian(isoprev)) {
      if (xdimprev != tsdimprev) SERR("reduced coordinates not allowed");
      if (strcmp(GLOBAL.coords.newunits[0], UNITS_NAMES[units_km]) == 0){
	*nr = isoprev == GNOMONIC_PROJ ? EARTHKM2GNOMONIC  
	  : isoprev == ORTHOGRAPHIC_PROJ ?  EARTHKM2ORTHOGRAPHIC
	  : EARTHKM2CART;
      } else if (strcmp(GLOBAL.coords.newunits[0], 
			UNITS_NAMES[units_miles]) == 0) {
	*nr = isoprev == GNOMONIC_PROJ ? EARTHMILES2GNOMONIC  
	  : isoprev == ORTHOGRAPHIC_PROJ ?  EARTHMILES2ORTHOGRAPHIC
	  :EARTHMILES2CART;
      } else {
	SERR4("only units '%s' and '%s' are allowed. Got '%s' (user's '%s').",
	      UNITS_NAMES[units_km], UNITS_NAMES[units_miles], 
	      GLOBAL.coords.newunits[0], GLOBAL.coords.curunits[0]);
      }
      if (isoprev == GNOMONIC_PROJ || isoprev == ORTHOGRAPHIC_PROJ) {
	*newtsdim = tsdimprev;
	*newxdim = xdimprev;
	*newisoprev = isoprev;
      } else {
	*newisoprev = CARTESIAN_COORD;
	if (tsdimprev == 2) {
	  //assert(xdimprev == 2);
	  *newtsdim = *newxdim = 3;
	} else {
	  *newtsdim = tsdimprev;
	  *newxdim = xdimprev;	  
	}	
      }
    } else {
      return ERRORODDCOORDTRAFO; // NotProgrammedYet("");
    }      
    break;
  case EARTH_ISOTROPIC : case SPHERICAL_ISOTROPIC :
    //printf("--- wrong iso\n");
    return ERRORWRONGISO;
  default:    
     return ERRORODDCOORDTRAFO; // NotProgrammedYet("");
  }
  return NOERROR;
}


int SetGatter(domain_type domprev, domain_type domnext, 
	      isotropy_type isoprev, isotropy_type isonext, 
	      //	      bool first, 
	      int *nr, int *delflag) {  
  if (domprev < domnext) 
    SERR2("Cannot call more complex models ('%s') from simpler ones ('%s')",
	  DOMAIN_NAMES[(int) domnext], DOMAIN_NAMES[(int) domprev]);

  //   printf("gatter dom=%d %d %s %s anySph=%d\n", domprev, domnext, ISONAMES[isoprev], ISONAMES[isonext], isAnySpherical(isonext));
  if (isAnySpherical(isonext)) {
    if (isoprev == UNREDUCED) {
      *nr = isEarth(isoprev) ? E2E : Sph2Sph;
      return NOERROR;
    }
    if(!isAnySpherical(isoprev)) BUG;
    if (domnext== XONLY) {
      if ((isoprev == EARTH_COORD || isoprev == SPHERICAL_COORD) 
	  && domprev == XONLY) {
	if (isAnySphericalIso(isonext) || isoprev == SPHERICAL_ISOTROPIC)
	  SERR("Change to isotropic earth coordinates not possible");
	*nr = isEarth(isoprev) 
	  ? (isSpherical(isonext) ? E2Sph : E2E)
	  : (isSpherical(isonext) ? Sph2Sph : ISO_MISMATCH);
	return NOERROR;
      }
      *nr = isEarth(isoprev) 
       ? (isSpherical(isonext) ? E2SphIso : E2EIso) 
       : (isSpherical(isonext) ? Sph2SphIso : ISO_MISMATCH);
    } else { // domnext == KERNEL, d.h. auch == domprev
      *nr = isEarth(isoprev) 
	? (isSpherical(isonext) ? E2Sph : E2E)
	: (isSpherical(isonext) ? Sph2Sph : ISO_MISMATCH);
    }
    if (domprev == XONLY && isoprev == isonext) *delflag=DEL_COV - 4;
    if (*nr == ISO_MISMATCH)
      SERR("change from spherical coordinates to earth coord's not possible");
    return NOERROR;
  }
    
  if (isoprev == CYLINDER_COORD || isonext ==  CYLINDER_COORD)
    SERR("general sphericaUNREDUCED,UNREDUCED,UNREDUCED,l coordinates not programmed yet");

  //////////////
  // cartesian
  //////////////

  assert(isCartesian(isoprev) || isoprev == UNREDUCED);
  bool isoOK = isoprev == isonext || 
    (isoprev > isonext && isonext <= CARTESIAN_COORD); // keine Trafo innerhalb
  // der projektionen -- bislang zumindest. ? to do
  if (!isoOK) 
    SERR2("cannot call more complex models ('%s') from simpler ones ('%s')",
	  ISONAMES[(int) domnext], ISONAMES[(int) domprev]);
  // if (isoprec == GNOMONIC_PROJ &&
     
  if (domnext == XONLY && 
      (domprev == KERNEL || 
       isoprev ==  CARTESIAN_COORD || 
       isoprev == SYMMETRIC ||
       isoprev == VECTORISOTROPIC || 
       isoprev == GNOMONIC_PROJ || 
       isoprev == ORTHOGRAPHIC_PROJ ||
       isoprev ==  ZEROSPACEISO)) {	
    switch (isonext) {
    case ISOTROPIC :
      *nr = S2ISO;
      break;
    case SPACEISOTROPIC :
      *nr = S2SP;
      break;
    case ZEROSPACEISO: case VECTORISOTROPIC: case SYMMETRIC: 
    case CARTESIAN_COORD: case GNOMONIC_PROJ: case ORTHOGRAPHIC_PROJ:
      *nr = S2S;
      // *delflag = DEL_COV - 5; ///
      break;
    case UNREDUCED:
      *nr = SId;
      break;
    default: BUG;
    }
  } else {    
    if (domprev == XONLY) {
      switch(isoprev) {
      case ISOTROPIC :
	*nr = ISO2ISO; 
	break;
      case SPACEISOTROPIC :
	*nr = (isonext == ISOTROPIC) ? SP2ISO : SP2SP;
	break;	 
      default: 
	PRINTF("GetGatter prev=%s; %s\n          next=%s %s\n",
	       DOMAIN_NAMES[domprev], ISONAMES[isoprev], 
	       DOMAIN_NAMES[domnext], ISONAMES[isonext]);
	BUG;
      }
    } else { // KERNEL domprev und domnext
      *nr = SId;
      *delflag = DEL_COV - 4;//
    }
  }

  // printf("prev %d %d ; next %d %d nr = %d\n", domprev, isoprev, domnext, isonext, *nr);

  return NOERROR;
}



int get_internal_ranges(cov_model *cov, cov_model *min, cov_model *max, 
			cov_model *pmin, cov_model *pmax, 
			cov_model *openmin, cov_model *openmax) {
  cov_fct *C = CovList + cov->nr; // nicht gatternr
  int 
      err = NOERROR, k = 0, i = 0,  // compiler dummy values
      kappas = C->kappas ;
  range_type range;
  rangefct getrange = C->range;
  SEXPTYPE *type = C->kappatype;
  
  //  print("entry\n");
  //PMI(cov, "gir");

  if (kappas > 0) {
    getrange(cov, &range); 
    
    // PMI(cov);
    assert (cov->maxdim >= cov->tsdim);
      
    err = NOERROR;
    for (i=0; i<kappas; i++) {
      int len=cov->ncol[i] * cov->nrow[i];
      double dmin, dmax, dpmin, dpmax, dopenmin, dopenmax,
	value =  RF_NA;
      dpmin = range.pmin[i];
      dpmax = range.pmax[i];
      dmin = range.min[i];
      dmax = range.max[i];
      dopenmin = (double) range.openmin[i];
      dopenmax = (double) range.openmax[i];

      if (type[i] == INTSXP) {
	if (dmin < -MAXINT) {
#ifdef RANDOMFIELDS_DEBUGGING
	  PRINTF("%s <min %s\n", NICK(cov), KNAME(i));
	  BUG;
#else 
	  dmin = (double) -MAXINT;
#endif
	}
	if (dmax > MAXINT) {
#ifdef RANDOMFIELDS_DEBUGGING
	  PRINTF("%s >max %s\n", NICK(cov), KNAME(i));
	  BUG;
#else 
	  dmax = (double) MAXINT;	  
#endif
	}
	
      }
     //    print("i=%d %d %s %f %d %d\n", i, len, C->name, 
      //	   PisNULL(i) ? RF_NA : P0(i][0], type[i], REALSXP);
      
      for (k=0; k<len; k++) {
	if (type[i] == REALSXP) {
	  value = P(i)[k];
	  PARAM(min, i)[k] = dmin;
	  PARAM(max, i)[k] = dmax;
	  PARAM(pmin, i)[k] = dpmin;
	  PARAM(pmax, i)[k] = dpmax;
	  PARAM(openmin, i)[k] = dopenmin;
	  PARAM(openmax, i)[k] = dopenmax;
	} else if (type[i] == INTSXP) {
	  value = PINT(i)[k] == NA_INTEGER ? RF_NA : (double) PINT(i)[k];
	  PARAMINT(min, i)[k] = dmin;
	  PARAMINT(max, i)[k] = dmax;
	  PARAMINT(pmin, i)[k] = dpmin;
	  PARAMINT(pmax, i)[k] = dpmax;
	  PARAMINT(openmin, i)[k] = dopenmin;
	  PARAMINT(openmax, i)[k] = dopenmax;
	} else if (type[i] == LISTOF + REALSXP) {
	  listoftype 
	    *p = PARAMLIST(min, i);

	  if (p->deletelist) {
	    int j, 
	      lenj = p->nrow[k] * p->ncol[k];
	    double
	      *qmin = PARAMLIST(min, i)->p[k],
	      *qmax = PARAMLIST(max, i)->p[k],
	      *ppmin = PARAMLIST(pmin, i)->p[k],
	      *ppmax = PARAMLIST(pmax, i)->p[k],
	      *omin = PARAMLIST(openmin, i)->p[k],
	      *omax = PARAMLIST(openmax, i)->p[k];
	    
	    for (j=0; j<lenj; j++) {
	      qmin[j] = dmin;
	      qmax[j] = dmax;
	      ppmin[j] = dpmin;
	      ppmax[j] = dpmax;
	      omin[j] = dopenmin;
	      omax[j] = dopenmax;
	    }
	  } // elses list had not been copied

	  value = RF_NA;
	  // error cannot appear as range is (-infty, +infty)
	} else if (type[i] == CLOSXP) {
          continue;   
	} else if (type[i] == LANGSXP) {
	  continue;   
	} else return ERRORUNKOWNSXPTYPE;

	if (ISNAN(value)) {
	  continue;
	}

	dmin = range.min[i];
	dmax = range.max[i];
	if (value < dmin
	    || value > dmax
	    || (range.openmin[i] && value == dmin)
	    || (range.openmax[i] && value == dmax)
	    ) {
	  sprintf(ERRORSTRING,
		  "value (%f) of '%s' in '%s' not within interval %s%f,%f%s",
		  value, C->kappanames[i], NICK(cov), 
		  range.openmin[i] ? "(" : "[", dmin, dmax,
		  range.openmax[i] ? ")" : "]"
		  );
	  return ERRORM;
	}
      } // for k
    } // for i in kappas;
  } // if (kappa > 0)

  for (i=0; i<MAXPARAM; i++) {
    if (cov->kappasub[i] != NULL) {
      err = get_internal_ranges(cov->kappasub[i], 
				min->kappasub[i], max->kappasub[i],  
				pmin->kappasub[i], pmax->kappasub[i],
				openmin->kappasub[i], openmax->kappasub[i]);
      if (err != NOERROR) return err;
    }
  }
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      err = get_internal_ranges(cov->sub[i], min->sub[i], max->sub[i], 
			       pmin->sub[i], pmax->sub[i],
			       openmin->sub[i], openmax->sub[i]);
      if (err != NOERROR) return err;
    }
  }
  return NOERROR;
}


void strround(double x, char *s){
  if (x==RF_INF)  sprintf(s, "Inf"); else
    if (x==RF_NEGINF)  sprintf(s, "-Inf"); else
      if (x == floor(x + 0.5)) sprintf(s, "%d", (int) x);
      else sprintf(s, "%f", x);
}




void addmsg(double value, const char *sign, double y, char *msg) {
  char str1[30], str2[30];
  strround(value, str1);
  strround(y, str2);
  sprintf(msg, "%s %s %s", str1, sign, str2);
}


void addone(char *str, double x) {
  char str2[30];
  strround(x, str2);
  sprintf(str, "%s, %s", str, str2);
}
void addpair(char *str, double x, double y) {
    if (x == floor(x + 0.5) && y==floor(y + 0.5))
	sprintf(str, "%s, (%d,%d)", str, (int) x, int(y));
    else 
	sprintf(str, "%s, (%f,%f)", str, x, y);
}


int check_within_range(cov_model *cov, bool NAOK) {
  // sets also maxdim and finiterange in cov !

  cov_fct *C = CovList + cov->nr; //nicht gatternr
  int len,
    err = NOERROR,
    k = 0, 
    i = 0,   // compiler dummy values
    kappas = C->kappas;
  range_type range;
  SEXPTYPE *type = C->kappatype;
  rangefct getrange = C->range;
  char Msg[255] = "";
  double min, max,
    value= RF_NA;

  if (GLOBAL.general.skipchecks) return NOERROR;

  getrange(cov, &range); 

  //  printf("maxdim = %d %d %s\n", cov->maxdim, cov->tsdim, Nick(cov));

  if (cov->maxdim >=0 && cov->maxdim < cov->tsdim) {
    GERR2("Max. dimension is %d. Got %d", cov->maxdim, cov->tsdim);
  }

   for (i=0; i<kappas; i++) {
     if (!strcmp(C->kappanames[i], FREEVARIABLE)) {
       // i.e. equal
       if (PisNULL(i)) continue;
     }
     if (type[i] >= LISTOF) {
       // to simplify things -- otherwise in simu.cc
       // code must also be changed.
       assert(range.min[i] == RF_NEGINF && range.max[i]==RF_INF);
       continue;
     }
    // full range !!
    
    len = cov->ncol[i] * cov->nrow[i];
    min = range.min[i];
    max = range.max[i];
    
    //  print("%s i=%d [%f, %f] size=%d %s\n", 
    //	    C->name, i,  range.min[i],  range.max[i], len, C->kappanames[i]);
  
    for (k=0; k<len; k++) {
      if (type[i] == REALSXP) value = P(i)[k];
      else if (type[i] == INTSXP)
	value = PINT(i)[k] == NA_INTEGER ? RF_NA : (double) PINT(i)[k];
      else if (type[i] == CLOSXP) continue;
      else if (type[i] == LANGSXP) continue;
      else GERR2("%s [%d] is not finite", C->kappanames[i], k+1);
      if (ISNAN(value)) {
	if (NAOK) {
	  continue;
	} else GERR2("%s[%d] is not finite.", C->kappanames[i], k+1);
      }
      
      //  print("k=%d %f\n", k, value);
	 // print("range %d %d %d op=%d %f %f %f\n", 
         // kappas, i, k, range.openmin[i], value, min, max);
      
      err = ERRORUNSPECIFIED;
 
      //  printf("%s %d\n", NAME(cov), i);
      //  printf("%d\n", range.openmin[i]);
      // printf("%f %f\n", value, min); 

      if (range.openmin[i] && value <= min) { 
	//crash();
	addmsg(value, ">", min, Msg);
	goto ErrorHandling;
      } else if (value < min) {
	addmsg(value, ">=", min, Msg); 
        goto ErrorHandling;
      }
      if (range.openmax[i] && value >= max) { 
	addmsg(value, "<", max, Msg); 
        goto ErrorHandling;
      } else if (value > max) { 
	addmsg(value, "<=", max, Msg);
	goto ErrorHandling;
      }
      err = NOERROR; 
    }
  } // kappas
  
ErrorHandling:

   // PMI(cov, "range");
 
  if (err != NOERROR) {
    // printf("i=%d\n", i);
    if (PL>PL_ERRORS)
      PRINTF("error in check range: %s kappa%d err=%d ('%s' does not hold.)\n",
	     C->name, i+1, err, Msg);
     if (err == ERRORUNSPECIFIED)
      SERR4("%s[%d] = %s does not hold (dim=%d).",
 	     C->kappanames[i], k+1, Msg, cov->tsdim); // + return
  }

  // print("done %s\n", NICK(cov));
 
  return err;
}

int get_ranges(cov_model *cov, cov_model **min, cov_model **max, 
		cov_model **pmin, cov_model **pmax, 
		cov_model **openmin, cov_model **openmax) {
  int err;
  // returns a reliable value only in the very first
  // entry of each parameter vector/matrix int err;
  if ((err = covcpy(min, cov, false)) != NOERROR) return err; 
  if ((err = covcpy(max, cov, false)) != NOERROR) return err; 
  if ((err = covcpy(pmin, cov, false)) != NOERROR) return err ;
  if ((err = covcpy(pmax, cov, false)) != NOERROR) return err;
  if ((err = covcpy(openmin, cov, false)) != NOERROR) return err; 
  if (( err = covcpy(openmax, cov, false)) != NOERROR) return err; 
  (*min)->calling = (*max)->calling = (*pmin)->calling = (*pmax)->calling = 
    (*openmin)->calling = (*openmax)->calling = cov;
  //  (*min)->calling = NULL; ja nicht auf NULL setzen, da sonst angenommen
  //  wird, dass prevloc ein Original ist

  return get_internal_ranges(cov, *min, *max, *pmin, *pmax, *openmin, *openmax);
}



int FieldReturn(cov_model *cov) {
   location_type *loc = Loc(cov);
  if (cov->rf != NULL) {
    //PMI(cov, "fieldreturn");
    assert(cov->fieldreturn && cov->origrf);
    if (cov->origrf) {
      UNCONDFREE(cov->rf);
    }
  } 

  //   printf("FieldReturn %d %d (%d %d)  %s\n",sizeof(res_type), loc->totalpoints, cov->vdim[0], cov->vdim[1], NICK(cov));
  if ((cov->rf = 
       (res_type*) MALLOC(sizeof(res_type) * loc->totalpoints * cov->vdim[0]))
      == NULL) return ERRORMEMORYALLOCATION;
  cov->fieldreturn = cov->origrf = true;
  return NOERROR;
}


void SetLoc2NewLoc(cov_model *cov, location_type *loc) {
  int i,
    maxsub =  CovList[cov->nr].maxsub;
  // printf("s2n %s %ld\n", NICK(cov), (long int) loc);
  if (cov->ownloc != NULL) {
    // printf("ownloc not null %s\n", NICK(cov));
    return;
  }
  for (i=0; i<MAXPARAM; i++) {
    if (cov->kappasub[i] != NULL) SetLoc2NewLoc(cov->kappasub[i], loc);
  }
  cov->prevloc = loc;
  for (i=0; i<maxsub; i++) {
    if (cov->sub[i] != NULL) SetLoc2NewLoc(cov->sub[i], loc);
  }
  if (cov->key != NULL)  SetLoc2NewLoc(cov->key, loc);
}



bool verysimple(cov_model *cov) {
  assert(cov != NULL);
  cov_fct *C = CovList + cov->nr; // kein gatternr, da isotrop
  int i, k, total,
    kappas = C->kappas;
  for (k=0; k<kappas; k++) {
    if (cov->kappasub[k] != NULL) return false;
    total = cov->ncol[k] * cov->nrow[k];
    if (C->kappatype[k] == REALSXP) {
      for (i=0; i<total; i++) 
	if (ISNAN(P(k)[i]) || ISNA(P(k)[i])) return false;
    } else  if (C->kappatype[k] == INTSXP) {
      for (i=0; i<total; i++) if (P(k)[i] == NA_INTEGER) return false;
    } else return false;
  }
  return true;
}
