/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of auxiliary correlation functions 

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is gno error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nonsta     tionary models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2014 Martin Schlather

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


#include <math.h>
 
#include "RF.h"
#include "Covariance.h"
//#include <R_ext/Lapack.h>
//#include <R_ext/Applic.h>
//#include <R_ext/Utils.h>     
//#include <R_ext/BLAS.h> 



void kdefault(cov_model *cov, int i, double v) {

  //printf("kdefault %s\n", NICK(cov));

  cov_fct *C = CovList + cov->nr; // nicht gatternr
  if (PisNULL(i)) {
    if (C->kappatype[i]==REALSXP) {
      PALLOC(i, 1, 1);
      P(i)[0] = v;
    } else if (C->kappatype[i]==INTSXP) {
      PALLOC(i, 1, 1); 
      if (ISNAN(v)) { BUG }
      else if (ISNA(v)) PINT(i)[0] = NA_INTEGER;
      else if (v > MAXINT) {BUG}
      else if (v < -MAXINT) {BUG}
      else PINT(i)[0] = (int) v;
    } else if (C->kappatype[i] == LISTOF + REALSXP) {
      //char msg[100];      
      PRINTF("%s:%s (%d) unexpected list\n", NICK(cov), C->kappanames[i], i);
      BUG;
    } else {
      //char msg[100];
      PRINTF("%s:%s (%d) is not defined\n", NICK(cov), C->kappanames[i], i);
      BUG;
    }
    cov->nrow[i] = cov->ncol[i] = 1;
  } else if (!GLOBAL.general.skipchecks) {    
    if (cov->nrow[i] != 1 || cov->ncol[i] != 1) { 
      LPRINT("%d %s %d nrow=%d, ncol=%d\n", 
	     cov->nr, NAME(cov), i, cov->nrow[i], cov->ncol[i]);
      int j; for (j=0; j<cov->ncol[i] * cov->nrow[i]; j++) {
	LPRINT("%f\n", P(i)[j]);
      }
      char param_name[100]; 
      strcpy(param_name, KNAME(i)); 
      PERR("parameter not scalar -- contact author.");
    }
  }
}


void updatepref(cov_model *cov, cov_model *sub) {
  int i;
  for (i=0; i< Forbidden; i++) {
    if (i==Specific) continue;
    if (sub->pref[i] < cov->pref[i]) {
      cov->pref[i] = sub->pref[i];
    }
  }
}


void setdefault(cov_model *cov) {
  // Achilles-Ferse: setdefault wird aufgerufen bevor die Varianten
  // der Modelle und Operatoren feststehen -- die Parameter unten muessen, falls
  // notwendig von Hand in den checks gesetzt werden.

//  cov_model *prev=cov->calling;
  int i;
  cov_fct *C = CovList + cov->nr;// nicht gatternr

  //  printf("setdefault %s\n", NAME(cov));
 
  cov->logspeed = RF_NA;

  cov->maxdim = C->maxdim;
  assert(C->vdim != MISMATCH);
  
  if ((C->vdim != PREVMODEL_DEP  &&  C->vdim != SUBMODEL_DEP)) {
     cov->vdim[0] = cov->vdim[1] = C->vdim;     
  }

  if (isPosDef(cov)) {
    for (i=0; i<MAXMPPVDIM; i++) cov->mpp.maxheights[i] = 1.0;    
  }

  cov->xdimown = cov->xdimprev;
  if (is_any(ISOTROPIC, C) && cov->isoown == ISOTROPIC && 
      isPosDef(cov->typus) &&
      C->domain==XONLY)
    cov->logspeed = 0.0;
 
  cov->finiterange = C->finiterange;  
  cov->monotone=C->Monotone;

  //   cov->diag = 
  //   cov->semiseparatelast = 
  //   cov->separatelast = isIsotropic(cov->isoown); // war false vor 27.9.14

  MEMCOPY(cov->pref, C->pref, sizeof(pref_shorttype)); 

  for (i=Nothing+1; i<Forbidden; i++) cov->pref[i] = PREF_NONE;

  //cov->deterministic = !isRandom(cov); // nur die die zuf. Param haben
  cov->method = Forbidden;

  cov->taylorN = C->TaylorN;
  cov->tailN = C->TailN;

  for (i=0; i<cov->taylorN; i++) {
    //    printf("taylor %d %d\n", i, TaylorConst);
    cov->taylor[i][TaylorConst] = C->Taylor[i][TaylorConst];
    cov->taylor[i][TaylorPow] = C->Taylor[i][TaylorPow];
  }
  for (i=0; i<cov->tailN; i++) {
    cov->tail[i][TaylorConst] = C->Tail[i][TaylorConst];
    cov->tail[i][TaylorPow] = C->Tail[i][TaylorPow];
    cov->tail[i][TaylorExpConst] = C->Tail[i][TaylorExpConst];
    cov->tail[i][TaylorExpPow] = C->Tail[i][TaylorExpPow];
  }

}


void set_integer(int *val, int sub){
  if (*val == SUBMODEL_DEP) *val = sub;
  else {    
    //    if (*val == PARAM_DEP) {cov_model *cov; crash(cov);}
    assert(*val != PARAM_DEP);
    if (sub < *val) *val = sub;
  }
}

void set_extbool(ext_bool *val, ext_bool sub){
  if (*val == SUBMODEL_DEP) *val = sub;
  else {
    if (sub < *val) *val = sub;
  }
}


// no setbackard
void setbackward(cov_model *cov, cov_model *sub) {
   cov_fct *C = CovList + cov->nr;// nicht gatternr
 // see also check_co when changing

  assert(cov != NULL);
  assert(sub != NULL);
 
  //  printf("prae maxdim %s %s %s sub=%d %d\n", C->nick, NAME(sub),
  //	 NAME(sub),	 sub->maxdim, cov->maxdim);
  set_integer(&(cov->maxdim), sub->maxdim);

  //  printf("maxdim %s sub=%d %d\n", C->nick, sub->maxdim, cov->maxdim);
  // printf("%d %d\n", (cov->monotone), sub->monotone);
  //printf("%d %d\n", (cov->finiterange), sub->finiterange);

  set_extbool(&(cov->monotone), sub->monotone);
  set_extbool(&(cov->finiterange), sub->finiterange);
  
  if (sub->full_derivs < cov->full_derivs)
      cov->full_derivs = sub->full_derivs;

  //PMI(cov);
  assert (cov->full_derivs >= -1 || 
	  (cov->full_derivs == MISMATCH && isRandom(cov) && isRandom(sub)));
  if (sub->rese_derivs < cov->rese_derivs)
      cov->rese_derivs = sub->rese_derivs;
  cov->loggiven &= sub->loggiven;
  
  updatepref(cov, sub);
  cov->tbm2num |= sub->tbm2num;

  if (C->vdim == SUBMODEL_DEP && (sub==cov->sub[0] || sub==cov->key)) {
    cov->vdim[0] = sub->vdim[0];
    cov->vdim[1] = sub->vdim[1];
  }

  cov->hess = (CovList[cov->nr].hess != NULL && sub->hess);
  cov->deterministic &= sub->deterministic;
}

int checkkappas(cov_model *cov, bool errornull){
  cov_fct *C = CovList + cov->nr; // nicht gatternr

  int i,nr, nc,
    kappas= C->kappas,
    *ncol = cov->ncol,
    *nrow = cov->nrow;
  char param_name[PARAMMAXCHAR]; // used in QERR

  //   printf("checkkappas %s\n", NAME(cov));
 
  for (i=0; i<kappas; i++) {
    strcpy(param_name, 
	   cov->ownkappanames != NULL && cov->ownkappanames[i]!=NULL 
	   ? cov->ownkappanames[i]
	   : C->kappanames[i]);
    cov_model *ks = cov->kappasub[i];     
  
    //    printf("checkkappas %s %s %d %d %d %d\n", NAME(cov), param_name, ks==NULL,
    //	   PisNULL(i), cov->nrow[i], cov->ncol[i]);
   //
    if (ks != NULL) {
      //printf("%d %d\n", isRandom(C->kappaParamType[i]), isRandom(ks));
      if (isRandom(ks)) {
	//printf("isRsndom %s %d %d %s %d \n", param_name, nrow[i], ncol[i], NAME(cov), i);
	cov->deterministic = false;
	int err, len;       
	
	nr = nrow[i];
	nc = ncol[i];
	if (nr==SIZE_NOT_DETERMINED || nc==SIZE_NOT_DETERMINED) 
	  C->kappasize(i, cov, &nr, &nc); 
	if (nr==SIZE_NOT_DETERMINED || nc==SIZE_NOT_DETERMINED) {
	  int d;
	  //printf("size: %s %d %d %s\n", param_name, nrow[i], ncol[i], NAME(cov));
	  for (d=9; d>=1; d--) {
	    //printf("d=%d\n", d);
	    err = CHECK_R(ks, d);
	    if (err != NOERROR && ks->vdim[0] > 0 && ks->vdim[0] != d) {
	      err = CHECK_R(ks, ks->vdim[0]);
	      d = 0;
	      //printf("d==%d %d\n", d, err);
	    }
	    if (err == NOERROR) {
	      nr = ks->vdim[0];
	      nc = ks->vdim[1];
	      len = nr * nc;
	      if (ks->nr == DISTRIBUTION) {
		if (PARAMINT(ks, DISTR_NROW) != NULL) 
		  nr = PARAM0INT(ks, DISTR_NROW);
		if (PARAMINT(ks, DISTR_NCOL) != NULL) 
		  nc = PARAM0INT(ks, DISTR_NCOL);
	      }
	      break;
	    }
	  }
	  if (err != NOERROR) return err;
	} else {
	  //if (nr==SIZE_NOT_DETERMINED || nc==SIZE_NOT_DETERMINED)
	  // QERR("size of random parameter could not be determined -- please give the size explicitely");	
	  len = nr * nc;	
	  if ((err = CHECK_R(ks, len)) != NOERROR)
	    QERRX(err, "random parameter not well defined");
	}
	
	if ( (ks->vdim[0] != nr || ks->vdim[1] != nc) &&
	     (ks->vdim[0] != len || ks->vdim[1] != 1) )
	  QERR("required size of random parameter does not match the model");
     
	if (cov->Sgen == NULL) NEW_STORAGE(gen);
	
	if (PisNULL(i)) {
	  PALLOC(i, nr, nc);
	}

       	//printf("nr=%d %d %f %s.%s\n", nr, nc, P(i)[0], NAME(cov), NAME(ks));
	
	if (!cov->initialised && 
	    (err = INIT_RANDOM(ks, 0, cov->Sgen, P(i))) != NOERROR) { 
	  // wird spaeter
	  // gegebememfalls nochmals initialisiert mit richtigen moments
	  //X4
	  QERRX(err, "random parameter cannot be initialized");	  
	}	
      } else { // not random, e.g. Aniso
	// no automatic check possible
	// if (!ks->checked) BUG; // fails
	// ks->checked = false; // fails as well
      }
    } // end ks != NULL (function given)


    if (PisNULL(i)) {
      if (errornull) { QERR("unset"); }
      else continue;    
    }
    
 
    C->kappasize(i, cov, &nr, &nc);   
    if ( (nc < 1 && nc != SIZE_NOT_DETERMINED) || 
	 (nr < 1 && nr != SIZE_NOT_DETERMINED)) { 
      BUG;
    }
    
    if (nc == 1 && ncol[i] != nc) {
      if (nr == 1 && nrow[i] != 1) QERR("must be a scalar");
      QERR("must be a vector, not a matrix");
    }
    if (nc > 1 && ncol[i] == 1) QERR("parameter must be a (larger) matrix");
    
    if ((nc > 0 && ncol[i] != nc) || (nr > 0 && nrow[i] != nr)) {
      
      // nc==0, nr==0 is coded as SIZE_NOT_DETERMINED
      char msg[255], msg2[255];
      sprintf(msg, "not of the given size: (%d, %d) instead of (",
	      nrow[i], ncol[i]);
      if (nr!=SIZE_NOT_DETERMINED) sprintf(msg2, "%s%d, ", msg, nr);
      else sprintf(msg2, "%sundetermined, ", msg);
      if (nc!=SIZE_NOT_DETERMINED) sprintf(msg, "%s%d)", msg2, nc);
      else sprintf(msg, "%sundetermined)", msg2);
      QERR(msg);
    }
    // if nc==0 (nr==0) then this is undetermined.
  } // for i < kappas

  // printf("end checkkappas\n");

  return NOERROR;
}

int checkkappas(cov_model *cov){
  return checkkappas(cov, true);
}



int alloc_mpp_M(cov_model *cov, int moments) {
  int maxmoments = CovList[cov->nr].maxmoments;
  // insbesondere fuer cov_models die selbst vom Random-Type sind
  assert(moments >= 0);
  assert(maxmoments != MISMATCH && maxmoments != PARAM_DEP);

  //  printf("alloc mpp %d %s %d \n", moments, NAME(cov), maxmoments);
  //if (!(maxmoments != MISMATCH && maxmoments != PARAM_DEP)) BUG;
 

  if (moments > maxmoments && maxmoments != SUBMODEL_DEP) {
    SERR2("required moments (%d) exceeds the coded moments (%d)",
	  moments, maxmoments);
  }
  if (moments <= cov->mpp.moments) return NOERROR;
  if (cov->mpp.mM != NULL) free_mpp_M(cov);
  cov->mpp.moments = moments;

  int i,
    vdim = cov->vdim[0],
    nm = cov->mpp.moments,
    nmvdim = (nm + 1) * vdim,
    bytes = sizeof(double) * nmvdim;
  if (vdim <= 0 || vdim > MAXMPPVDIM) BUG;
 
  //  printf("%d %d %d\n", cov->mpp.moments, bytes, vdim);

  cov->mpp.mM = (double*) MALLOC(bytes);
  cov->mpp.mMplus = (double*) MALLOC(bytes);

  //  printf("here\n");

  //  assert(nm < 100);
  int nmP1 = cov->mpp.moments + 1;
  for (i=0; i<nmvdim; i++) cov->mpp.mMplus[i] = cov->mpp.mM[i] = RF_NA; 
  for (i=0; i<vdim; i++) {
    int idx = i * nmP1;
    cov->mpp.mMplus[idx + 0] = cov->mpp.mM[idx + 0] = RF_INF;
  }

  //  printf("done\n");

  // cov->mpp.mMplus[0] = cov->mpp.mM[0] = 1.0;
  return NOERROR;
}

void free_mpp_M(cov_model *cov) {
  free(cov->mpp.mM);
  free(cov->mpp.mMplus);
  cov->mpp.mM = cov->mpp.mMplus = NULL;
}

int UpdateMPPprev(cov_model * cov, int moments) {
  cov_model *prev = cov->calling;
  int i, nm, err,
    nmvdim,
    vdim = cov->vdim[0];

  nm = cov->mpp.moments < prev->mpp.moments ? cov->mpp.moments
	: prev->mpp.moments;
  nmvdim = (nm + 1) * vdim;
  if (moments >= 0 && prev != NULL) {
    if (prev->mpp.moments == SUBMODEL_DEP &&
	(err = alloc_mpp_M(prev, moments)) != NOERROR) return err;
    for (i=0; i<nmvdim; i++) {
      //print("nm=%d %d %d %s %d\n", nm, cov->mpp.moments, prev->mpp.moments, NAME(prev), i);
      prev->mpp.mMplus[i] = cov->mpp.mMplus[i];
      prev->mpp.mM[i]     = cov->mpp.mM[i];
    }
  }
  
  // nachfolgende Zeilen so lassen, da sonst unerwuenscht
  // maxheight etc. nach oben gegeben werden.
  // prev->mpp.maxheight = cov->mpp.maxheight; 
  // prev->mpp.unnormedmass = cov->mpp.unnormedmass;
  
  return NOERROR;
}


int INIT_intern(cov_model *cov, int moments, gen_storage *s) { // kein err  
  if (!cov->checked) BUG;
  if (cov->initialised) return NOERROR;
  assert(cov != NULL);

  //PMI(cov, -1);

  assert(cov->gatternr >= ISO2ISO && cov->gatternr <= LASTGATTER);

  cov_fct *C = CovList + cov->nr;
  int err = NOERROR;
  sprintf(ERROR_LOC, "in %s: ", NICK(cov));
  // printf("errorloc %s\n", ERROR_LOC);
  
  if (cov->mpp.moments != SUBMODEL_DEP && cov->mpp.moments != PARAM_DEP) {
    if ((err = alloc_mpp_M(cov, moments)) != NOERROR) return err;
  } else {
    BUG; //assert(false); // passt das hier??
    if (cov->mpp.moments == PARAM_DEP) cov->mpp.moments = moments;
  }
 
  if (C->maxmoments >= 0 && moments > C->maxmoments) {
      SERR3("moments known up to order %d for '%s', but order %d required",
	    C->maxmoments, NICK(cov), moments);
  }
  
  sprintf(ERROR_LOC, "%s : ", cov->calling == NULL ? "initiating the model"
	  : NICK(cov->calling));
  ASSERT_GATTER(cov);
  if ((err = CovList[cov->gatternr].Init(cov, s)) != NOERROR) {
    return err;   
  }
   
  if ((err = UpdateMPPprev(cov, moments)) != NOERROR) {
    return err;
  }

  cov->initialised = true;
  return NOERROR;
}

void set_initialised_false(cov_model *cov){
  set_initialised_false(cov, false);
}

void set_initialised_false(cov_model *cov, bool init_deterministic){
  int i;
  if (!init_deterministic && cov->deterministic) return;
  cov->initialised = false;

  for (i=0; i<MAXPARAM; i++) {
    if (cov->kappasub[i] != NULL) {
      set_initialised_false(cov->kappasub[i]);
    }
  }

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL)
      set_initialised_false(cov->sub[i]);
  }
}


int REINIT_intern(cov_model *cov, int moments, gen_storage *s) { // kein err
  int err;
  set_initialised_false(cov);
  err = INIT_intern(cov, moments, s);
  return err;
}


int INIT_RANDOM_intern(cov_model *cov, int moments, gen_storage *s, // kein err
		       double *p) {
  if (!cov->checked) BUG;
  if (!cov->initialised) {
    int err = NOERROR;
    
    sprintf(ERROR_LOC, "in %s : ", NICK(cov));
     
    assert(cov != NULL);
    if (moments < 0) SERR("moments expected to be positive");
    if (CovList[cov->nr].maxmoments >= 0 &&
	moments > CovList[cov->nr].maxmoments) SERR("Moments do not match");
    
    if (cov->mpp.moments != SUBMODEL_DEP && cov->mpp.moments != PARAM_DEP) {
      if ((err = alloc_mpp_M(cov, moments)) != NOERROR) return err;
    } else {
      BUG; // passt das hier??
      if (cov->mpp.moments == PARAM_DEP) cov->mpp.moments = moments;
    }
    
    sprintf(ERROR_LOC, "%s:", cov->calling == NULL ? "initiating the model"
	    : NICK(cov->calling));
    ASSERT_GATTER(cov);
    if ((err = CovList[cov->gatternr].Init(cov, s)) != NOERROR) return err;   
    if (ISNAN(cov->mpp.mM[moments])) {
      SERR1("%s is not a random function", NICK(cov));
    }
    
    if ((err = UpdateMPPprev(cov, moments)) != NOERROR) return err;
    
    cov->initialised = true;
  }

  //  switch (CovList[cov->nr].kappatype[param_nr]) {
  //  case REALSXP :   
  DORANDOM(cov, p);
    

    //    break;
    //  case INTSXP :
    //    int j, len;
    //    double *dummy;
    //    dummy = (double*) MALLOC(sizeof(double) * len);
    //    DORANDOM(cov, dummy);
    //    for (j=0; j<len; j++) p[j] = (int) dummy[j];
    //    free(dummy);
    //    break;
    //  default : SERR("random parameter only allowed for numerical values");
    //  }
  
  return NOERROR;
}




void D_2(double *x, cov_model *cov, double *v){
  cov_fct *C = CovList + cov->nr;// nicht gatternr


  //printf("\n here %d %s \n", cov->isoown, NAME(cov)); 

  if (cov->xdimprev == 1) {
    assert(cov->isoown == ISOTROPIC);
    double y = fabs(*x);

    //   cov_model *prev = cov;
    //while (prev->calling != NULL) prev = prev->calling;
    
    // iso2iso ueberpruefen !!!
    // dollar + scale // aniso > 0 und dim = 1
    // nach tbm eigentlich auch nicht

    C->D(&y, cov, v);// nicht gatternr
  } else {
    assert(cov->xdimprev == 2);
    if (cov->xdimown == 1) {
      assert(cov->isoown==ISOTROPIC);
      double y=sqrt(x[0] * x[0] + x[1] * x[1]);    
      C->D(&y, cov, v);
      if (y!=0.0) *v *= x[0] / y;
    } else {
      assert(cov->xdimown == 2);
      double y[2];
      y[0] = fabs(x[0]);
      y[1] = fabs(x[1]);
      C->D(y, cov, v); 
    }
  }
}
void DD_2(double *x, cov_model *cov, double *v) {
  cov_fct *C = CovList + cov->nr;// nicht gatternr
  if (cov->isoown == ISOTROPIC) {
    double y = fabs(*x);

    // iso2iso ueberpruefen !!!
    // dollar + scale // aniso > 0 und dim = 1
    // nach tbm eigentlich auch nicht
    C->D2(&y, cov, v);// nicht gatternr
  } else {
    assert(cov->isoown == SPACEISOTROPIC);
    assert(cov->xdimprev == 2);
     
    if (is_all(ISOTROPIC, C)) {
      double w,
	xSq = x[0] * x[0],
	tSq = x[1] * x[1],
	ySq = xSq + tSq,
	y   = sqrt(ySq); 
      
      // (c'(r) * x/r)' = c''(r) * x^2/r^2 + c'(r) [ 1/r - x^2 / r^3]
      C->D2(&y, cov, v);// nicht gatternr
      if (y != 0.0) {
	C->D(&y, cov, &w);
	w /= y;
	*v = (*v - w) * xSq / ySq + w;
      } else {    
	// nothing to do ?
	// *v = x[0] / y;
      }
    } else if (is_all(SPACEISOTROPIC, C)) {
      double y[2];
      y[0] = fabs(x[0]);
      y[1] = fabs(x[1]);
      C->D2(y, cov, v); // nicht gatternr
    } else BUG;
  }
}

void DD_3(double VARIABLE_IS_NOT_USED *x, cov_model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *v) {
  error("DD_3 to be programmed\n");
}


void inverse2(double *x, cov_model *cov, double *v) {
  cov_fct *C = CovList + cov->nr;// nicht gatternr
  C->inverse(x, cov, v);//  nicht gatternr
}
   

void nonstatinverse2(double *v, cov_model *cov, double *x, double *y){
  cov_fct *C = CovList + cov->nr;// nicht gatternr

  C->nonstat_inverse(v, cov, x, y);//  nicht gatternr

  //  printf("inverse2 %f %f %f\n", v[0], x[0], y[0]);
}

void nonstat_loginverse2(double *v, cov_model *cov, double *x, double *y){
  cov_fct *C = CovList + cov->nr;// nicht gatternr

  C->nonstat_loginverse(v, cov, x, y);//  nicht gatternr

  //  printf("inverse2 %f %f %f\n", v[0], x[0], y[0]);
}
   
   

int struct2(cov_model *cov, cov_model **newmodel) {
  int err;
  char errloc_save[nErrorLoc];
  if (!cov->checked) {
    BUG;
  }
  strcpy(errloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "In %s : ", NICK(cov));

  // printf("\nstart struct %s\n", CovList[cov->nr].nick);
  

  err = CovList[cov->nr].Struct(cov, newmodel);
  if (newmodel != NULL && (*newmodel) != NULL) {
    (*newmodel)->calling = cov->calling != NULL ? cov->calling : cov;
  }

  //  assert(cov->nr != || cov->key->nr != );
  //  
  //printf("\nstruct %s", CovList[cov->nr].nick);

  if (err == NOERROR) strcpy(ERROR_LOC, errloc_save);

  return err;
}

int init2(cov_model *cov, gen_storage *s){ // s wird durchgereicht!

  // printf("init2 %s %s\n", NAME(cov), ROLENAMES[cov->role]);

  cov_fct *C = CovList + cov->nr; //  nicht gatternr
  cov_model
    *prev = cov->calling == NULL ? cov : cov->calling;
  int i,
    err = NOERROR,
    kappas = CovList[cov->nr].kappas;
  char errloc_save[nErrorLoc];
  strcpy(errloc_save, ERROR_LOC);
   
  PrInL++;

  for (i=0; i<kappas; i++) {
    cov_model *param  = cov->kappasub[i];
    if (param != NULL && isRandom(param)) {
      //printf("%f ", P0(i));
      if ((err = INIT_RANDOM(param, 0, s, P(i))) != NOERROR) return err;
    }
  }

  //  printf("prev %s %d %d %ld\n", C->nick, prev->method,  Forbidden, (long int) s);

  if (cov->method == Forbidden) cov->method = prev->method;
  if (cov->role == ROLE_GAUSS) {
    //  printf("here A\n");
    if (cov->method==SpectralTBM) {
      if (cov->calling == NULL && cov->nr != SPECTRAL_PROC_USER &&
	  cov->nr != SPECTRAL_PROC_INTERN) {
	SERR("unexpected value in init2");
      }
    }
   
    if (!cov->initialised && (err=C->Init(cov, s)) != NOERROR) {
      //printf("fehelr init\n");
      goto ErrorHandling;
    }
  } 

  else if (cov->role == ROLE_BERNOULLI) {
    if (!cov->initialised && (err=C->Init(cov, s)) != NOERROR) {
      //printf("fehelr init\n");
      goto ErrorHandling;
    }
  } 

  else if (hasAnyShapeRole(cov)) {
    cov->origrf = false;
    assert((cov->mpp.moments < 0) xor (cov->mpp.mM != NULL));

    //    if (cov->mpp.moments >= 0 && isRandom(cov->typus) && cov->mpp.mM == NULL) {
    //// insbesondere fuer cov_models die selbst vom Random-Type sind
      //alloc_mpp_M(cov);
    //    }
    
    sprintf(ERROR_LOC, "In %s: ", NICK(cov));
    if (!cov->initialised && (err = C->Init(cov, s)) != NOERROR) {
 	goto ErrorHandling;
    }
 
    sprintf(ERROR_LOC, "'%s': ", NICK(prev));//  nicht gatternr   
    err = NOERROR;  
  }

  else if (hasNoRole(cov)) {
    if (!cov->initialised && (err=C->Init(cov, s)) != NOERROR) {
      goto ErrorHandling;      
    }
  }

  else {
    ILLEGAL_ROLE;
  }

  prev->fieldreturn    = cov->fieldreturn;
  
 ErrorHandling :
  PrInL--;
  if (err == NOERROR) strcpy(ERROR_LOC, errloc_save);
  cov->initialised = err == NOERROR;
  return err;
}

void do2(cov_model *cov, gen_storage *s){
  //   cov_model *prev = cov->calling == NULL ? cov : cov->calling;
  //
  //  int i,
  //    kappas = CovList[cov->nr].kappas;

  // statt nachfolgende Zeilen: siehe init2
  //  for (i=0; i<kappas; i++) {
  //    cov_model *param  = cov->kappasub[i];
  //    if (param != NULL && isRandom(param)) DORANDOM(param, P(i));
  //  }


//printf("Do prev %s %d %d %ld\n",
//       NAME(prev), prev->method, Forbidden, (long int) s);
//printf("Do current %s %d %d %ld\n",
//       NAME(cov), cov->method, Forbidden, (long int) s);
 CovList[cov->nr].Do(cov, s); // ok

  // assert(false);
}


void dorandom2(cov_model *cov, double *v){
 CovList[cov->nr].DoRandom(cov, v); // ok
}


// from pos def to neg def
int CheckPD2ND(cov_model *cov, int tsdim, int tsxdim, isotropy_type isoprev,
	       int vdim, int role) {
#define nsel 2
  int i, j, err=NOERROR, 
    statselect[nsel] = {XONLY, KERNEL};
  // statselect[nsel]={STATIONARY, VARIOGRAM, COVARIANCE, GEN_VARIOGRAM};
  Types typeselect[nsel] = {PosDefType, NegDefType};
  
  for (i=0; i<nsel; i++) {
    for (j=0; j<nsel; j++) {
      if ((err = check2X(cov, tsdim, tsxdim, typeselect[j], statselect[i],
			 isoprev, vdim, role)) == NOERROR) return err;
    }
  }
  return err;
}

// to do: Trend komplett anders behandeln -- im moment laeuft er als Anhaengsel zur Kovarianzstruktur. Sollte der Trend seperat behandelt werden ??!!






int check2X(cov_model *cov, int tsdim, int tsxdim,
	    Types type, domain_type domprev, isotropy_type isoprev,
	    int vdim, int role) {
  return check2X(cov, tsdim, tsxdim, type, domprev, isoprev, vdim, vdim, role);
}

int check2X(cov_model *cov, int tsdim, int tsxdim,
	    Types type, domain_type domprev, isotropy_type isoprev,
	    int *vdim, int role) {
  return check2X(cov, tsdim, tsxdim, type, domprev, isoprev, vdim[0], vdim[1],
		 role);
}


int check2X(cov_model *cov, int tsdim, int tsxdim,
	    Types type, domain_type domprev, isotropy_type isoprev,
	    int vdim0, int vdim1, int role) {
  //if (cov == NULL) crash();
  //printf(">>>> %s: isoprev=%s dowprev=%s \n", NAME(cov), ISONAMES[isoprev], DOMAIN_NAMES[domprev]); 
  //APMI(cov);

  assert(cov != NULL);
  assert(vdim0 != 0 && vdim1 != 0);
  int UnUsedDeleteFlag, err;
  cov_model *prev = cov->calling == NULL ? cov : cov->calling;
  cov_fct *P = CovList + prev->nr, //  nicht gatternr
    *C = CovList + cov->nr;        //  nicht gatternr
  isotropy_type iso0; // iso
  domain_type dom, first_dom, last_dom; 
  bool checkerror = false,
    skipchecks = GLOBAL.general.skipchecks;
  char msg[1000] = "";

  // erst bei check unten
  sprintf(ERROR_LOC, "'%s' : ", NICK(cov));
  if (PL >= PL_COV_STRUCTURE) { 
    if (cov->calling == NULL) PRINTF("\n");
    LPRINT("%s\n", ERROR_LOC); 
  }

  //  
 
  
  cov->domprev = domprev;
  cov->isoprev = isoprev;
  cov->tsdim = tsdim; // muss wegen checkkappas gesetzt werden
  cov->role = role;  
  cov->typus = type;
  cov->xdimprev = cov->xdimgatter = tsxdim; //if cov is isotropy or 
  //                                          spaceisotropic it is set to 1 or 2
 
 
  if (tsxdim < 1) { 
    SERR("dimension less than 1"); 
  }
  if (role == ROLE_UNDEFINED) SERR("role undefined");
  if (isUndefined(type)) {    
    SERR("type undefined");
  }

  if (cov->calling != NULL && isInterface(cov)) {
    SERR1("'%s' may be used only as top model", NICK(cov));
  }

  
  //  printf("neg=%d pos=%d iso=%d SUMM=%d %s\n", isNegDef(type),  isPosDef(type), isoprev, SYMMETRIC, NAME(cov));

  //  assert(!isNegDef(type) || isPosDef(type) || isoprev <= SYMMETRIC || ({PMI(cov); false;})); //


  if (cov->calling != NULL) {
    cov->prevloc = Loc(prev);
  }
  //  printf("cov->calling %ld %ld\n", cov->calling, cov->prevloc);
 

  if (PL >= PL_STRUCTURE) {
    LPRINT("#[%s -> %s] (%d; %d (%d)): \n", cov->calling == NULL ? "NULL" : 
	   Nick(prev), Nick(cov), domprev, isoprev, C->Isotropy[0]); 
  }
 
  isotropy_type isolist[LAST_ISO + MAXVARIANTS];
  int origlastiso, 
    idxiso=-1;

  if (isPrevModelI(C)) {    
    isotropy_type last_iso;
    // printf("A\n");
    if (isCartesian(isoprev)) {
      for (last_iso=ISOTROPIC; last_iso<=LAST_CARTESIAN; last_iso++) {
	isolist[++idxiso] = last_iso;
      }
    } else if (isSpherical(isoprev)) {
      for (last_iso=SPHERICAL_ISOTROPIC; last_iso<=SPHERICAL_COORD; last_iso++){
	isolist[++idxiso] = last_iso;
      } 
    } else if (isEarth(isoprev)) {
      // wichtig, dass so lange wie moeglich altes System beibehalten
      // wird, damit 'scale' des Benutzers in der richtigen Einheit
      // verarbeitet wird
      for (last_iso=EARTH_ISOTROPIC; last_iso<=EARTH_COORD; last_iso++){
	isolist[++idxiso] = last_iso;
      } 
      for (last_iso=SPHERICAL_ISOTROPIC; last_iso<=SPHERICAL_COORD; last_iso++){
	isolist[++idxiso] = last_iso;
      } 
 

      //assert(isoprev >= 0);
      /*
      if (C->domain != KERNEL && isAnySpherical(isoprev)) {
	if (isEarth(isoprev)) isolist[++idxiso] = EARTH_ISOTROPIC;
	isolist[++idxiso] = SPHERICAL_ISOTROPIC;
	if (isEarth(isoprev)) isolist[++idxiso] = EARTH_COORD;
	isolist[++idxiso] = SPHERICAL_COORD;
      }
      */
    } else BUG;
  } else if (isUnreduced(C)) {
    //printf("BA\n");
    if (isCartesian(isoprev)) {     
      isolist[++idxiso] = CARTESIAN_COORD;      
    } else if (isSpherical(isoprev)) {
      isolist[++idxiso] = SPHERICAL_COORD;      
    } else if (isEarth(isoprev)) {
      isolist[++idxiso] = EARTH_COORD;      
    } else {
      isolist[++idxiso] = isoprev != UNREDUCED ? isoprev :
      cov->calling == NULL ? MISMATCH : cov->calling->isoown;
      assert(isolist[0] != UNREDUCED);
    }
  } else { // including C->isotropy == ISO_MISMATCH
    int i;
    for(i=0; i<C->variants; i++) 
      if (idxiso<0 ||  C->Isotropy[i] != isolist[idxiso])
	isolist[++idxiso] = C->Isotropy[i];  
    assert(idxiso<0 || 
	   (isolist[idxiso] != ISO_MISMATCH && isolist[idxiso] != PREVMODELI));
  }

  origlastiso = idxiso; // only for messages
  //  while (idxiso >= 0 && !atleastSpecialised(isolist[idxiso], isoprev)) 
  //   idxiso--;
  // if (idxiso < 0)
  //   SERR2("required isotropy '%s' cannot be fullfilled by '%s'",
  //	  ISONAMES[isoprev], NAME(cov));
  int //t,
 i = 0;
  //print("dd %d %d\n", i, idxiso);

  while (i <= idxiso) { // die Typen ausschliessen, die sicher nicht gehen
    // printf("CC i=%d idxiso=%d; %d %d %d\n", i, idxiso, isolist[idxiso], isoprev, atleastSpecialised(isolist[i], isoprev));
    cov->isoown = isolist[i];
    if ( TypeConsistency(type, cov, 0) && 
	equal_coordinate_system(isoprev, isolist[i]) &&
	atleastSpecialised(isolist[i], isoprev)
	) i++;
    else {
      int j;
      for (j=i; j<idxiso; j++) isolist[j] = isolist[j+1];
      idxiso--;
    }
  }
   
  //    printf("%s: %d prev=%d, first_iso %d (%s) last=%d %s C=%s\n", NAME(cov), idxiso, isoprev, isolist[0], ISONAMES[isolist[0]], isolist[idxiso], C->nick, ISONAMES[C->Isotropy[0]]);

  if (idxiso < 0) {
    if (PL >= PL_COV_STRUCTURE) 
      PRINTF("error as non-isotropic model cannot be called from isotropic one (%s -> %s)\n", ISONAMES[(int) isoprev], ISONAMES[(int) cov->isoown]);

    if (cov->calling == NULL) SERR("basic isotropy assumption does not match");
    if (cov->calling->calling == NULL) {
      SERR2("model is a '%s' function, but at least a '%s' function is required for the given specification.", 
	    ISONAMES[(int) C->Isotropy[0]], ISONAMES[(int) isoprev]);
    } else {
      SERR5("model '%s' has property '%s'. It cannot be called by '%s' which requires the property '%s' (%d)",
	    NICK(cov),
	    ISONAMES[isolist[0]], Nick(prev), ISONAMES[isolist[origlastiso]],
	    isolist[origlastiso]);
    }
  }
  
  
  first_dom = last_dom = cov->domown = C->domain;
  if (first_dom == PREVMODELD) {
    first_dom = XONLY;
    last_dom = KERNEL; // 10.10.11: GENERALISEDCOVARIANCE;
  }  
  if (last_dom > domprev) last_dom = domprev;


  //   printf("last %d %d dom=%d %d orig=%d (%s) prev=%d\n", isolist[idxiso], isolist[0],last_dom,first_dom, C->domain, C->name, domprev);
 
  if (last_dom < first_dom) {
    if (cov->calling == NULL) {
       BUG;
    }

    if (PL >= PL_COV_STRUCTURE) 
      PRINTF("model called from less complex one (%s:%s;%s -> %s:%s [%s;%s])\n",
	     cov->calling == NULL ? "NULL" : 
	     Nick(prev),  DOMAIN_NAMES[(int)CovList[prev->nr].domain], 
	     DOMAIN_NAMES[(int)domprev], NAME(cov), 
	     DOMAIN_NAMES[(int) C->domain],
	     DOMAIN_NAMES[(int) first_dom],  DOMAIN_NAMES[(int) last_dom]);
    if (cov->calling->calling == NULL) {
      SERR2("model is a '%s', but at least a '%s' is required for the given specification.",  DOMAIN_NAMES[(int) C->domain], DOMAIN_NAMES[(int) domprev]);
    } else {
      SERR1("Model cannot be called from '%s'", Nick(prev));
    }
  }
  if (PL >= PL_STRUCTURE) {
    LPRINT("(dom.start=%d, end=%d, iso.start=%d, end=%d)\n",
	   first_dom, last_dom, isolist[0], isolist[idxiso]);
  }
  
  err = ERRORNOSTATMATCH;
  int *nr = NULL;
  isotropy_type  newisoprev = MISMATCH;
  for (dom = first_dom; dom <= last_dom; dom++) {
    char checkmsg[LENERRMSG];
    int t, err3;
    cov->domown = dom;
    for (i=0; i<=idxiso; i++) {
      cov->isoown = iso0 = isolist[i];
     //     printf("dom=%d i=%d (%s) of %d\n", dom, i, ISONAMES[isolist[i]], idxiso);

      cov->full_derivs = C->F_derivs;
      cov->rese_derivs = C->RS_derivs;
      cov->loggiven = C->log != ErrLogCov; 

      nr = &(cov->gatternr); 
      newisoprev = isoprev;
      cov->tsdim = tsdim;
      cov->vdim[0] = vdim0;
      cov->vdim[1] = vdim1;
      
      //print("xx %s  xxx type = %s\n", NICK(cov), TYPENAMES[type]);
      if ((t = TypeConsistency(type, cov, 1))) {
	//printf("HERE %d %s\n", t, TYPENAMES[C->Typi[t - 1]]);
	if (!isUndefined(C->Typi[t - 1])) {
	  cov->typus = C->Typi[t - 1];
	}
      } else {
  	SERR3("required type '%s' does not match the type '%s' of '%s'",
	      TYPENAMES[type], TYPENAMES[CovList[cov->nr].Typi[t]], 
	      NICK(cov));
      }
   
      setdefault(cov); // braucht cov->isoown && typus!

      if ((err3 = checkkappas(cov, C->primitive)) != NOERROR)  return err3;

      if (PL >= PL_STRUCTURE) {
	if ((dom>first_dom || i>0)) {
	  LPRINT("");
	  //	  MERR(err); //
	}
	if (first_dom==last_dom) {
	  LPRINT("[%s]; [%s] sys=%d,%d\n", DOMAIN_NAMES[(int) first_dom], 
		 ISONAMES[isolist[i]], isoprev, cov->isoown );
	} else {
	  LPRINT("[%s..%s]:%s; [%s..%s]:%s sys=%d,%d\n",
	       DOMAIN_NAMES[(int) first_dom], DOMAIN_NAMES[(int) last_dom],  
	       DOMAIN_NAMES[(int) dom], 
		 ISONAMES[isolist[i]], ISONAMES[isolist[idxiso]],
	       ISONAMES[(int) iso0], isoprev, cov->isoown);
	}
      } 

     
      //        printf(">> %s: prev=%s new.prev=%s own:%s\n", NAME(cov), ISONAMES[isoprev], ISONAMES[newisoprev], ISONAMES[iso0]); 
      //APMI(cov);

      assert(equal_coordinate_system(isoprev, cov->isoown));
      int err2,
	newtsdim = tsxdim;
      if (cov->calling != NULL && 
	  !equal_coordinate_system(cov->calling->isoown, cov->isoown)) {

	// 	printf("%s: changing system: %s (%d) to %s (%d)\n", NAME(cov), ISONAMES[cov->calling->isoown], cov->calling->isoown,  ISONAMES[isoprev] ,isoprev);
	
	if ((err2 = change_coordinate_system(cov->calling->isoown, isoprev, 
					    tsdim, cov->xdimprev,
					    nr, &newisoprev, &newtsdim, 
					    &(cov->xdimgatter)))
	    != NOERROR) {
	  if (err == NOERROR || err == ERRORNOSTATMATCH) err = err2;
	  continue;
	//	printf("done\n");
	}

	if (isEarth(cov->calling->isoown) && (err = checkEarth(cov)) !=NOERROR){
	  //	  printf("earth error %d\n", err); MERR(err);
	  continue;
	}
	//printf("done\n");
	cov->xdimown = cov->tsdim = newtsdim;
	nr = &(cov->secondarygatternr);     
      }
      
      cov->xdimown = dom == KERNEL ? newtsdim
	: iso0 == ISOTROPIC ? 1 : iso0 == SPACEISOTROPIC ? 2 	
	: iso0 == EARTH_ISOTROPIC ||  iso0 == SPHERICAL_ISOTROPIC
	? cov->xdimprev - 1	
	: newtsdim;

       //printf("iso0=%d nmewtsdim=%d dom=%d\n", iso0, newtsdim, dom);

      if (cov->xdimown > cov->xdimprev && newtsdim <= tsxdim) { // appear if spaceiso called by iso
	if (checkerror) {
	   SERR2("%s: %s", NICK(cov), checkmsg);
	 } else {
	   SERR2("dimension at least %d needed. Got %d dimension.", 
		 cov->xdimown, cov->xdimprev);
	 }
       }

      //printf("check = %d %s %ld %ld\n", err, NICK(cov), (long int) C->check, (long int) check_directGauss);
        err = C->check(cov); // CHECK !
	//   	printf("err = %d %s\n", err, NICK(cov)); MERR(err);

 
       checkerror = err != NOERROR;
       if (checkerror) {
	 //	 PMI(cov);
	 errorMSG(err, checkmsg, LENERRMSG);
       } else { 
	 if (C->maxdim>=0 && cov->maxdim > C->maxdim) {
	   cov->maxdim = C->maxdim;
	 }
	 if (cov->vdim[0] <= 0 || cov->vdim[1] <= 0) {
	   return ERRORBADVDIM;
	 } 
	
	if ((vdim0 > 0 && cov->vdim[0] != vdim0) || 
	    (vdim1 > 0 && cov->vdim[1] != vdim1)) {
	  sprintf(ERRORSTRING,
		  "multivariate dimension (of submodel '%s'), which is %d x %d, does not match the specification of '%s', which is %d x %d and is required by '%s'",
		NICK(cov), cov->vdim[0], cov->vdim[1], C->name, vdim0, vdim1, 
		cov->calling == NULL ? "-- none --" :  P->name); 
	  return ERRORWRONGVDIM; // needed as value!
	  checkerror = true;
	}
	if (cov->monotone == PARAM_DEP) BUG;

	break;
      }
  
       //      if (err == ERRORINCORRECTSTATISO) {
       //	if (strcmp("", msg) != 0) {
       //	  err = ERRORM; 
       //	  strcpy(ERRORSTRING, msg);
       //	  strcpy(ERROR_LOC, "");
       //	}
       //	continue;
       //      } else 
       if (err > NOERROR) {
	 errorMSG(err, msg); 
       }
       
    } // for iso

    if (err == NOERROR) break;
  } // dom

   //  printf("ok\n");
 
  if (PL > PL_COV_STRUCTURE && cov->calling == NULL) {
    LPRINT("%s: end look ", Nick(cov)); 
    if (err != NOERROR) PRINTF("err = %d\n", err); else MERR(err); //
  }

  if (err != NOERROR) return err;


  if (PL >= PL_COV_STRUCTURE) {
    LPRINT("Continuing '%s' (no error):\n", Nick(cov)); 
  }

  if (!skipchecks && (err = check_within_range(cov, NAOK_RANGE)) != NOERROR) { 
    return err;
  }
 
  if (isoprev == SPACEISOTROPIC) {
    //print("\n\n\n");
    cov_model *cv = cov;
    while(cv->calling != NULL) cv = cv->calling;
   
    if (cov->xdimown != 2) {
      return ERRORDIM;
    }
    if (cov->tsdim < 2)  {
      return ERRORDIM;
    }
  }

 
  //      printf("gattered %s %d %d %d %d\n", NAME(cov), domprev, cov->domown, newisoprev, cov->isoown);
  //     PMI(cov, -1);
 
  err = SetGatter(domprev, cov->domown, 
		  newisoprev, cov->isoown, 
		  //  cov->calling == NULL,
		  nr, &UnUsedDeleteFlag); //

   //  printf("setgatter  err = %d\n", err);
  //if (iso0 == GNOMONIC_PROJ) { printf("\n\n\n\n\n !!!! \n\n\n\n\n\n\n\n");
    //APMI(cov);  }
  ASSERT_GATTERONLY(cov);

  if (PL > PL_COV_STRUCTURE) {
    LPRINT("leaving '%s' for '%s'  SetGatter error=%d deriv=%d,%d \n",
	   Nick(cov),
	   cov->calling == NULL ? "none" : Nick(prev),
	   err, cov->full_derivs, cov->rese_derivs);
  }


  sprintf(ERROR_LOC, "\"%s\": ", cov->calling == NULL ? "parsing the model"
	  : Nick(prev));

  // printf("end err = %d\n", err);
  COND_NEW_STORAGE(gatter, z); 
  if (isAnySpherical(cov->isoown)) COND_NEW_STORAGE(earth, X); 

  // printf("2err = %d %s <= %s\n", err, NAME(cov), cov->calling == NULL ? "----" : NAME(cov->calling));
  cov->checked = err == NOERROR;
  assert(err == NOERROR || (cov->vdim[0] > 0 && cov->vdim[1] > 0));
  
  return(err);
   
}

