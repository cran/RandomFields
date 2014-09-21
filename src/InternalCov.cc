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

  // PMI(cov);

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
	     cov->nr, NICK(cov), i, cov->nrow[i], cov->ncol[i]);
      int j; for (j=0; j<cov->ncol[i] * cov->nrow[i]; j++) {
	LPRINT("%f\n", P(i)[j]);
      }
      char param_name[100]; 
      strcpy(param_name, CovList[cov->nr].kappanames[i]); 
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
     cov->vdim2[0] = cov->vdim2[1] = C->vdim;     
  }

  if (isPosDef(cov)) {
    for (i=0; i<MAXMPPVDIM; i++) cov->mpp.maxheights[i] = 1.0;    
  }

  cov->xdimown = cov->xdimprev;
  if (C->isotropy == ISOTROPIC && isPosDef(cov->typus) &&
      C->domain==XONLY) 
    cov->logspeed = 0.0;
 
  cov->finiterange = C->finiterange;  
  cov->monotone=C->Monotone;

  cov->diag = 
    cov->semiseparatelast = 
    cov->separatelast = true;

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
 //  PMI(cov);
 //   PMI(sub);
 //   if (sub == NULL) crash();

  assert(cov != NULL);
  assert(sub != NULL);
 
  //  printf("prae maxdim %s %s %s sub=%d %d\n", C->nick, NICK(sub),
  //	 NICK(sub),	 sub->maxdim, cov->maxdim);
  set_integer(&(cov->maxdim), sub->maxdim);

  //  printf("maxdim %s sub=%d %d\n", C->nick, sub->maxdim, cov->maxdim);
  // printf("%d %d\n", (cov->monotone), sub->monotone);
  //printf("%d %d\n", (cov->finiterange), sub->finiterange);

  //printf("%s mono %d %d\n", NICK(cov), (cov->monotone), sub->monotone);
  set_extbool(&(cov->monotone), sub->monotone);
  //printf("range %d %d\n", (cov->finiterange), sub->finiterange);
  set_extbool(&(cov->finiterange), sub->finiterange);
  
  cov->diag &= sub->diag;
//  cov->quasidiag &= sub->quasidiag;
  cov->separatelast &= sub->separatelast;
  cov->semiseparatelast &= sub->semiseparatelast;  
  if (sub->full_derivs < cov->full_derivs)
      cov->full_derivs = sub->full_derivs;
  //  PMI(cov);
  assert (cov->full_derivs >= -1 || 
	  (cov->full_derivs == MISMATCH && isRandom(cov) && isRandom(sub)));
  if (sub->rese_derivs < cov->rese_derivs)
      cov->rese_derivs = sub->rese_derivs;
  cov->loggiven &= sub->loggiven;
  
  updatepref(cov, sub);
  cov->tbm2num |= sub->tbm2num;

  if (C->vdim == SUBMODEL_DEP && (sub==cov->sub[0] || sub==cov->key)) {
    cov->vdim2[0] = sub->vdim2[0];
    cov->vdim2[1] = sub->vdim2[1];
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

  //   printf("checkkappas %s\n", NICK(cov));
 
  for (i=0; i<kappas; i++) {
    strcpy(param_name, 
	   cov->ownkappanames != NULL && cov->ownkappanames[i]!=NULL 
	   ? cov->ownkappanames[i]
	   : C->kappanames[i]);
    cov_model *ks = cov->kappasub[i];     
  
    //    printf("checkkappas %s %s %d %d %d %d\n", NICK(cov), param_name, ks==NULL,
    //	   PisNULL(i), cov->nrow[i], cov->ncol[i]);
   //
    if (ks != NULL) {
      //printf("%d %d\n", isRandom(C->kappaParamType[i]), isRandom(ks));
      if (isRandom(ks)) {
	//	PMI(cov);
	//printf("isRsndom %s %d %d %s %d \n", param_name, nrow[i], ncol[i], NICK(cov), i);
	cov->deterministic = false;
	int err, len;       
	
	nr = nrow[i];
	nc = ncol[i];
	if (nr==SIZE_NOT_DETERMINED || nc==SIZE_NOT_DETERMINED) 
	  C->kappasize(i, cov, &nr, &nc); 
	if (nr==SIZE_NOT_DETERMINED || nc==SIZE_NOT_DETERMINED) {
	  int d;
	  //printf("size: %s %d %d %s\n", param_name, nrow[i], ncol[i], NICK(cov));
	  for (d=9; d>=1; d--) {
	    //printf("d=%d\n", d);
	    err = CHECK_R(ks, d);
	    if (err != NOERROR && ks->vdim2[0] > 0 && ks->vdim2[0] != d) {
	      err = CHECK_R(ks, ks->vdim2[0]);
	      d = 0;
	      //printf("d==%d %d\n", d, err);
	    }
	    if (err == NOERROR) {
	      nr = ks->vdim2[0];
	      nc = ks->vdim2[1];
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
	
	if ( (ks->vdim2[0] != nr || ks->vdim2[1] != nc) &&
	     (ks->vdim2[0] != len || ks->vdim2[1] != 1) )
	  QERR("required size of random parameter does not match the model");
     
	if (cov->stor == NULL)
	  cov->stor = (gen_storage *) MALLOC(sizeof(gen_storage));
	
	if (PisNULL(i)) {
	  PALLOC(i, nr, nc);
	}

	
	// APMI(cov); 	
	//printf("nr=%d %d %f %s.%s\n", nr, nc, P(i)[0], NICK(cov), NICK(ks));
	
	if (!cov->initialised && 
	    (err = INIT_RANDOM(ks, 0, cov->stor, P(i))) != NOERROR) { 
	  // wird spaeter
	  // gegebememfalls nochmals initialisiert mit richtigen moments
	  //X4
	  //PMI(ks); XERR(err);
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
      // PMI(cov);
      
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

  //  printf("alloc mpp %d %s %d \n", moments, NICK(cov), maxmoments);
  //if (!(maxmoments != MISMATCH && maxmoments != PARAM_DEP)) BUG;
 

  if (moments > maxmoments && maxmoments != SUBMODEL_DEP) {
    //APMI(cov);
    SERR2("required moments (%d) exceeds the coded moments (%d)",
	  moments, maxmoments);
  }
  if (moments <= cov->mpp.moments) return NOERROR;
  if (cov->mpp.mM != NULL) free_mpp_M(cov);
  cov->mpp.moments = moments;

  int i,
    vdim = cov->vdim2[0],
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
    vdim = cov->vdim2[0];

  //PMI(prev);

  nm = cov->mpp.moments < prev->mpp.moments ? cov->mpp.moments
	: prev->mpp.moments;
  nmvdim = (nm + 1) * vdim;
  //   print("nm=%d %d %d %s\n", nm, cov->mpp.moments, prev->mpp.moments, NICK(prev)); crash(); APMI(prev);
  if (moments >= 0 && prev != NULL) {
    if (prev->mpp.moments == SUBMODEL_DEP &&
	(err = alloc_mpp_M(prev, moments)) != NOERROR) return err;
    // PMI(prev);
    for (i=0; i<nmvdim; i++) {
      //print("nm=%d %d %d %s %d\n", nm, cov->mpp.moments, prev->mpp.moments, NICK(prev), i);
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
  assert(cov->gatternr >= ISO2ISO && cov->gatternr <= SId);

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
 
  //  printf("\n >>>> moment %s %d %d\n\n", NICK(cov), moments, C->maxmoments); //PMI(cov);

  if (C->maxmoments >= 0 && moments > C->maxmoments) {
      //printf("moment %s %d %d\n", NICK(cov), moments, C->maxmoments);PMI(cov);
      SERR3("moments known up to order %d for '%s', but order %d required",
	    C->maxmoments, NICK(cov), moments);
  }
  
  sprintf(ERROR_LOC, "%s: ", cov->calling == NULL ? "initiating the model"
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
    
    sprintf(ERROR_LOC, "in %s: ", NICK(cov));
     
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
    
    sprintf(ERROR_LOC, "%s: ", cov->calling == NULL ? "initiating the model"
	    : NICK(cov->calling));
    //PMI(cov);
    ASSERT_GATTER(cov);
    if ((err = CovList[cov->gatternr].Init(cov, s)) != NOERROR) return err;   
    if (ISNAN(cov->mpp.mM[moments])) {
      // printf("moments = %d %s\n", moments, NICK(cov));
      //APMI(cov);
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



/* simplifying functions
   turn vector of x into length ||x|| and similar
   reduces C(x,y) to C(x - y)

   needs preceeding analysis of submodels!
   (same for preferred methods; bottom-up-analysis needed)
*/

// #
// keep always the orderung
void iso2iso(double *x, cov_model *cov, double *v) {
  double y=fabs(*x);
  CovList[cov->nr].cov(&y, cov, v); // nicht gatternr
}
void logiso2iso(double *x, cov_model *cov, double *v, double *sign) {
  double y=fabs(*x);
  CovList[cov->nr].log(&y, cov, v, sign);// nicht gatternr
}
void spiso2spiso(double *x, cov_model *cov, double *v) {
  double y[2];
  y[0] = fabs(x[0]);
  y[1] = fabs(x[1]);
  CovList[cov->nr].cov(y, cov, v);// nicht gatternr
}
void logspiso2spiso(double *x, cov_model *cov, double *v, double *sign) {
  double y[2];
  y[0] = fabs(x[0]);
  y[1] = fabs(x[1]);
  CovList[cov->nr].log(y, cov, v, sign);// nicht gatternr
}
void spacetime2iso(double *x, cov_model *cov, double *v) {
  double y=sqrt(x[0] * x[0] + x[1] * x[1]);
  CovList[cov->nr].cov(&y, cov, v);// nicht gatternr
}
void logspacetime2iso(double *x, cov_model *cov, double *v, double *sign) {
  double y=sqrt(x[0] * x[0] + x[1] * x[1]);
  CovList[cov->nr].log(&y, cov, v, sign);// nicht gatternr
}

void Stat2iso(double *x, cov_model *cov, double *v) {
  double b = 0.0;
  int i,
    dim=cov->xdimgatter;  
  
  //   print("%d\n", dim); assert(false);
  //  print("dim=%d, %f %f, %d %ld\n", dim, x[0], x[1], cov->nr, cov->calling);
  //   print("%s\n", NICK(cov));

  //     PMI(cov->calling);
  //  print("end PMI %d\n", dim);

  for (i=0; i<dim; i++) {
    //printf("d=%d b=%f %f\n", i, b, x[i]);
    b += x[i] *x[i];
    //  
  }
  b = sqrt(b);
  //  APMI(cov->calling->calling);  crash(cov);
  // printf("xxx\n");
  CovList[cov->nr].cov(&b, cov, v);// nicht gatternr
  //
  //printf("aaxxx\n");
  // print("%f\n", *v);
}
void logStat2iso(double *x, cov_model *cov, double *v, double *sign) {
  double b = 0.0;
  int i,
    dim=cov->xdimgatter;  
  for (i=0; i<dim; i++) {
    b += x[i] *x[i];
  }
  b = sqrt(b);
  CovList[cov->nr].log(&b, cov, v, sign);// nicht gatternr
  // print("%f\n", *v);
}
void Nonstat2iso(double *x, double *y, cov_model *cov, double *v) {
  double a, b;
  int d,
    dim=cov->xdimgatter;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d];
    b += a * a;
  }
  b = sqrt(b);
  //if (dim != 2) APMI(cov);
  //assert(dim == 2);
  CovList[cov->nr].cov(&b, cov, v);// nicht gatternr
  //printf("b=%f %f %s\n", b, *v, NICK(cov));
}
void logNonstat2iso(double *x, double *y, cov_model *cov, double *v,
		    double *sign) {
  double a, b;
  int d, 
    dim=cov->xdimgatter;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d];
    b += a * a;
  }
  b = sqrt(b);
  CovList[cov->nr].log(&b, cov, v, sign);// nicht gatternr
}
void Stat2spacetime(double *x, cov_model *cov, double *v) {
  double b, z[2];
  int i,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, i=0; i<dim; i++)  b += x[i] * x[i];
  z[0] = sqrt(b);
  z[1] = fabs(x[dim]);
  CovList[cov->nr].cov(z, cov, v);// nicht gatternr
}
void logStat2spacetime(double *x, cov_model *cov, double *v, double *sign) {
  double b, z[2];
  int i,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, i=0; i<dim; i++)  b += x[i] * x[i];
  z[0] = sqrt(b);
  z[1] = fabs(x[dim]);
  CovList[cov->nr].log(z, cov, v, sign);// nicht gatternr
}
void Nonstat2spacetime(double *x, double *y, cov_model *cov, double *v) {
  double a, b, z[2];
  int d,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d];
    b += a * a;
  }
  z[0] = sqrt(b);
  z[1] = fabs(x[dim] - y[dim]);
  CovList[cov->nr].cov(z, cov, v);// nicht gatternr
}
void logNonstat2spacetime(double *x, double *y, cov_model *cov, double *v,
			  double *sign) {
  double a, b, z[2];
  int d,
    dim=cov->xdimgatter - 1;  
  for (b=0.0, d=0; d<dim; d++) {
    a = x[d] - y[d];
    b += a * a;
  }
  z[0] = sqrt(b);
  z[1] = fabs(x[dim] - y[dim]);
  CovList[cov->nr].log(z, cov, v, sign);// nicht gatternr
}
void Stat2Stat(double *x, cov_model *cov, double *v) {
  CovList[cov->nr].cov(x, cov, v);// nicht gatternr
}
void logStat2Stat(double *x, cov_model *cov, double *v, double *sign) {
  CovList[cov->nr].log(x, cov, v, sign);// nicht gatternr
}

#define nonstat2statInner						\
  int d,								\
    dim=cov->xdimgatter;							\
  double *z = cov->S2->z;						\
  if (z == NULL) z = cov->S2->z = (double*) MALLOC(dim * sizeof(double)); \
  for (d=0; d<dim; d++) z[d] = x[d] - y[d]
 
void Nonstat2Stat(double *x, double *y, cov_model *cov, double *v) {
  nonstat2statInner;
  CovList[cov->nr].cov(z, cov, v);// nicht gatternr
}
void logNonstat2Stat(double *x, double *y, cov_model *cov, double *v, 
		     double *sign) {
  nonstat2statInner;
  CovList[cov->nr].log(z, cov, v, sign);// nicht gatternr
}

void Nonstat2Nonstat(double *x, double *y, cov_model *cov, double *v) {
  CovList[cov->nr].nonstat_cov(x, y, cov, v);// nicht gatternr
}

void logNonstat2Nonstat(double *x, double *y, cov_model *cov, double *v, 
			double *sign) {
  CovList[cov->nr].nonstatlog(x, y, cov, v, sign);// nicht gatternr
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
    
    if (C->isotropy==ISOTROPIC) {
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
    } else {
      assert(C->isotropy == SPACEISOTROPIC);
      double y[2];
      y[0] = fabs(x[0]);
      y[1] = fabs(x[1]);
      C->D2(y, cov, v); // nicht gatternr
    }
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
    //PMI(cov->calling->calling);
    //    crash();
    BUG;
  }
  strcpy(errloc_save, ERROR_LOC);
  sprintf(ERROR_LOC, "in %s: ", NICK(cov));

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

  // printf("init2 %s %s\n", NICK(cov), ROLENAMES[cov->role]);

  cov_fct *C = CovList + cov->nr; //  nicht gatternr
  cov_model
    *prev = cov->calling == NULL ? cov : cov->calling;
  int i,
    err = NOERROR,
    kappas = CovList[cov->nr].kappas;
  char errloc_save[nErrorLoc];
  strcpy(errloc_save, ERROR_LOC);
   
  PrInL++;

  //  PMI(cov);

  for (i=0; i<kappas; i++) {
    cov_model *param  = cov->kappasub[i];
    if (param != NULL && isRandom(param)) {
      //printf("%f ", P0(i));
      if ((err = INIT_RANDOM(param, 0, s, P(i))) != NOERROR) return err;
      //printf("init2 %s->%s->%s %s %d  i=%d %f \n", NICK(cov->calling), NICK(cov), NICK(param), KNAME(i), param != NULL && isRandom(param),  i, P0(i));
      
      // PMI(cov);

    }
  }

  //  printf("prev %s %d %d %ld\n", C->nick, prev->method,  Forbidden, (long int) s);

  if (cov->method == Forbidden) cov->method = prev->method;

  //    printf("role == %d %d", cov->role, ROLE_COV);
  //  PMI(cov);

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
     
    ///cov->mpp.methnr = prev->mpp.methnr;
    cov->origrf = false;
   
    //printf("here B\n");
    //PMI(cov, -1);

    assert((cov->mpp.moments < 0) xor (cov->mpp.mM != NULL));

    //    if (cov->mpp.moments >= 0 && isRandom(cov->typus) && cov->mpp.mM == NULL) {
    //// insbesondere fuer cov_models die selbst vom Random-Type sind
      //alloc_mpp_M(cov);
    //    }
    
    sprintf(ERROR_LOC, "in %s: ", NICK(cov));
    // if (C->maxsub == 0) cov->mpp.loc_done = false;

    //APMI(cov);

    if (!cov->initialised && (err = C->Init(cov, s)) != NOERROR) {
 	goto ErrorHandling;
    }
 
    sprintf(ERROR_LOC, "%s: ", NICK(prev));//  nicht gatternr   
    err = NOERROR;  
  }

  else if (hasNoRole(cov)) {
    //   printf("here C\n");
    if (!cov->initialised && (err=C->Init(cov, s)) != NOERROR) {
	//	printf("fehelr init\n");
      goto ErrorHandling;      
    }
  }

  else {
    //     printf("here D\n");
    // APMI(cov->calling);
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
//       NICK(prev), prev->method, Forbidden, (long int) s);
//printf("Do current %s %d %d %ld\n",
//       NICK(cov), cov->method, Forbidden, (long int) s);
 CovList[cov->nr].Do(cov, s); // ok

  // assert(false);
}


void dorandom2(cov_model *cov, double *v){
  // cov_model    *prev = cov->calling == NULL ? cov : cov->calling;
  //

//printf("Do prev %s %d %d %ld\n",
//       NICK(prev), prev->method, Forbidden, (long int) s);
//printf("Do current %s %d %d %ld\n",
//       NICK(cov), cov->method, Forbidden, (long int) s);
  //PMI(cov);

 CovList[cov->nr].DoRandom(cov, v); // ok

  // assert(false);
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
  assert(cov != NULL);
  assert(vdim0 != 0 && vdim1 != 0);
  int UnUsedDeleteFlag, err;
  cov_model *prev = cov->calling == NULL ? cov : cov->calling;
  cov_fct *P = CovList + prev->nr, //  nicht gatternr
    *C = CovList + cov->nr;        //  nicht gatternr
  isotropy_type first_iso, last_iso, iso0; // iso
  domain_type dom, first_dom, last_dom; 
  bool checkerror = false,
    skipchecks = GLOBAL.general.skipchecks;
  char msg[1000] = "";

  // erst bei check unten
  sprintf(ERROR_LOC, "%s: ", NICK(cov));
  if (PL >= PL_COV_STRUCTURE ) { 
    if (cov->calling == NULL) PRINTF("\n");
    LPRINT("%s\n", ERROR_LOC); 
  }

  cov->domprev = domprev;
  cov->isoprev = isoprev;
  cov->tsdim = tsdim; // muss wegen checkkappas gesetzt werden
  cov->role = role;  
  cov->typus = type;
  cov->xdimprev = cov->xdimgatter = tsxdim; //if cov is isotropy or 
  //                                          spaceisotropic it is set to 1 or 2
  
 
  if (tsxdim < 1) { 
    //APMI(cov);
    SERR("dimension less than 1"); 
  }
  if (role == ROLE_UNDEFINED) SERR("role undefined");
  if (type == UndefinedType) {    
    //APMI(cov);
    SERR("type undefined");
  }

  if (cov->calling != NULL && isInterface(cov)) {
    //APMI(cov);
    SERR1("'%s' may be used only as top model", NICK(cov));
  }
  if (!TypeConsistency(type, cov)) {
    //     printf("type =%d %s\n", type, TYPENAMES[type]);    PMI(cov->calling);
     //    crash();
    SERR3("required type '%s' does not match the type '%s' of '%s'",
	  TYPENAMES[type], TYPENAMES[CovList[cov->nr].Type], NICK(cov));
  }

  //  if (isRandom(type) && isoprev != CARTESIAN_COORD) SERR("");
  //
  //     PMI(cov);


  //  printf("neg=%d pos=%d iso=%d SUMM=%d %s\n", isNegDef(type),  isPosDef(type), isoprev, SYMMETRIC, NICK(cov));
  // if (!(!isNegDef(type) || isoprev <= SYMMETRIC))    APMI(cov);

  assert(!isNegDef(type) || isPosDef(type) || isoprev <= SYMMETRIC || ({PMI(cov); false;})); //


  if (cov->calling != NULL) {
    cov->prevloc = Loc(prev);
   }
  //  printf("cov->calling %ld %ld\n", cov->calling, cov->prevloc);
 

  if (PL >= PL_STRUCTURE) {
    LPRINT("#[%s -> %s] (%d %d): \n", cov->calling == NULL ? "NULL" : 
	   Nick(prev), Nick(cov), domprev, isoprev); 
  }
 
  //
  cov->isoown = C->isotropy;
  if (C->isotropy == PREVMODELI) {
    first_iso = ISOTROPIC;
    last_iso = CARTESIAN_COORD;
  } else if (C->isotropy == UNREDUCED) {    
    last_iso = first_iso = 
      isCartesian(isoprev) ? CARTESIAN_COORD :
      isoprev != UNREDUCED ? isoprev :
      cov->calling == NULL ? MISMATCH : cov->calling->isoown;
    assert(last_iso != UNREDUCED);
  } else first_iso = last_iso = C->isotropy;

  

  //   printf("prev=%d, first_iso %d (%s) last=%d %s C=%s\n", isoprev,first_iso,  ISONAMES[first_iso], last_iso, C->nick, ISONAMES[C->isotropy]);

  if (last_iso > isoprev) {
    if (isCartesian(isoprev)) last_iso = isoprev;
    else NotProgrammedYet("non Cartesian");
  }
  if (last_iso < first_iso) {
    //     PMI(cov->calling->calling, "iso"); assert(false);
    // printf("last %d %d prev=%d %s %d\n", last_iso, first_iso, isoprev, C->name, C->isotropy);
    //  crash();
    if (PL >= PL_COV_STRUCTURE) 
      PRINTF("error as non-isotropic model cannot be called from isotropic one (%s -> %s)\n", ISONAMES[(int) isoprev], ISONAMES[(int) cov->isoown]);

    // APMI(cov);
    if (cov->calling == NULL) SERR("basic isotropy assumption does not match");
    if (cov->calling->calling == NULL) {
      SERR2("model is a '%s' function, but at least a '%s' function is required for the given specification.", 
	    ISONAMES[(int) C->isotropy], ISONAMES[(int) isoprev]);
    } else {
      // PMI(cov, "cannot be called");
      //  crash();
      SERR5("model '%s' has property '%s'. It cannot be called by '%s' which requires the property '%s' (%d)",
	    NICK(cov),
	    ISONAMES[(int) first_iso], Nick(prev), ISONAMES[(int) last_iso],
	    (int) last_iso);
    }
  }
  
  first_dom = last_dom = cov->domown = C->domain;
  if (first_dom == PREVMODELD) {
    first_dom = XONLY;
    last_dom = KERNEL; // 10.10.11: GENERALISEDCOVARIANCE;
  }  
  if (last_dom > domprev) last_dom = domprev;


  //printf("last %d %d dom=%d %d orig=%d (%s) prev=%d\n",
    //	   last_iso, first_iso,last_dom,first_dom,
    //	    C->domain, C->name, domprev);



  if (last_dom < first_dom) {

    //PMI(cov, "here");
    //  printf("last %d %d dom=%d %d orig=%d (%s) prev=%d\n",
    //	   last_iso, first_iso,last_dom,first_dom,
    //	   C->domain, C->name, domprev);

    // crash(cov);

    if (cov->calling == NULL) {
      //  
      //APMI(cov); //     crash(cov);
      BUG;
    }

    ///  APMI(cov);
    if (PL >= PL_COV_STRUCTURE) 
      PRINTF("model called from less complex one (%s:%s;%s -> %s:%s [%s;%s])\n",
	     cov->calling == NULL ? "NULL" : 
	     Nick(prev), 
	     STATNAMES[(int)CovList[prev->nr].domain], 
	     STATNAMES[(int)domprev],
	     NICK(cov), 
	     STATNAMES[(int) C->domain],
	     STATNAMES[(int) first_dom],  STATNAMES[(int) last_dom]
	  );
    if (cov->calling->calling == NULL) {
      //APMI(cov);
      SERR2("model is a '%s', but at least a '%s' is required for the given specification.",  STATNAMES[(int) C->domain], STATNAMES[(int) domprev]);
    } else {
      SERR1("Model cannot be called from '%s'", Nick(prev));
    }
  }
  
  err = ERRORNOSTATMATCH;
  if (PL >= PL_STRUCTURE) {
    LPRINT("(dom.start=%d, end=%d, iso.start=%d, end=%d)\n",
	   first_dom, last_dom, first_iso, last_iso);
  }

 
  int *nr = NULL;
  isotropy_type  newisoprev = MISMATCH;
  for (dom = first_dom; dom <= last_dom; dom++) {
    char checkmsg[LENERRMSG];
    cov->domown = dom;
    for (iso0=first_iso; iso0 <= last_iso; iso0++) {

      // printf("iso=%d\n", iso0);

      cov->full_derivs = C->F_derivs;
      cov->rese_derivs = C->RS_derivs;
      cov->loggiven = C->log != ErrLogCov; 

      nr = &(cov->gatternr); 
      newisoprev = isoprev;
      cov->tsdim = tsdim;
      cov->isoown = iso0;
      cov->vdim2[0] = vdim0;
      cov->vdim2[1] = vdim1;
      setdefault(cov);
 
      if ((err = checkkappas(cov, C->primitive)) != NOERROR)  return err;

      if (PL >= PL_STRUCTURE) {
	if ((dom>first_dom || iso0>first_iso)) {
	  LPRINT("");
	  MERR(err);
	}
	
	LPRINT("[%s%s%s]%s%s; [%s%s%s]%s%s sys=%d,%d\n",
	       STATNAMES[(int) first_dom], 
	       first_dom==last_dom ? "" : "..",
	       first_dom==last_dom ? "" : STATNAMES[(int) last_dom],  
	       first_dom==last_dom ? "" : ":",
	       first_dom==last_dom ? "" : STATNAMES[(int) dom], 
	       ISONAMES[(int) first_iso],
	       first_iso==last_iso ? "" : "..",
	       first_iso==last_iso ? "" : ISONAMES[(int) last_iso],
	       first_iso==last_iso ? "" : ":",
	       first_iso==last_iso ? "" : ISONAMES[(int) iso0],//, ISONAMES[iso]
	       isoprev, cov->isoown
	       );
      }

      //    printf("call=%ld\n", cov->calling);
      //   printf("%s %s \n", ISONAMES[isoprev], ISONAMES[cov->isoown]);

      assert(equal_coordinate_system(isoprev, cov->isoown));
      int newtsdim = tsxdim;
      if (cov->calling != NULL && 
	  !equal_coordinate_system(cov->calling->isoown, cov->isoown)) {
       	//

	//PMI(cov->calling);

	//	printf("%d %d %s %s\n", cov->calling->isoown ,isoprev, ISONAMES[cov->calling->isoown], ISONAMES[isoprev]);

	
	if ((err = change_coordinate_system(cov->calling->isoown, isoprev, nr, 
					    &newisoprev, &newtsdim, 
					    &(cov->xdimgatter)))
	    != NOERROR) return err;
	if (isEarth(cov->calling->isoown) && (err = checkEarth(cov)) != NOERROR)
	  continue;

	cov->xdimown = cov->tsdim = newtsdim;
	nr = &(cov->secondarygatternr);     
      }
      
       cov->xdimown = dom == KERNEL ? newtsdim
	: iso0 == ISOTROPIC ? 1 : iso0 == SPACEISOTROPIC ? 2 
	: newtsdim;

       //printf("iso0=%d nmewtsdim=%d dom=%d\n", iso0, newtsdim, dom);

       if (cov->xdimown > cov->xdimprev && newtsdim <= tsxdim) { // appear if spaceiso called by iso
	 //	 PMI(cov);
	 if (checkerror) {
	   SERR2("%s: %s", NICK(cov), checkmsg);
	 } else {
	   SERR2("dimension at least %d needed. Got %d dimension.", 
		 cov->xdimown, cov->xdimprev);
	 }
       }

       //       cov_model *calling = cov->calling;
       err = C->check(cov);

       //       printf("err = %d\n", err);

       checkerror = err != NOERROR;
       if (checkerror) {
	 errorMSG(err, checkmsg, LENERRMSG);
       } else {

	 if (C->maxdim>=0 && cov->maxdim > C->maxdim) {
	   cov->maxdim = C->maxdim;
	 }

	 if (cov->vdim2[0] <= 0 || cov->vdim2[1] <= 0) {
	   // 
	   //	   PMI(cov);
	   //	   printf("%d %d\n", vdim0, vdim1);
	   return ERRORBADVDIM;
	 } 
	
	if ((vdim0 > 0 && cov->vdim2[0] != vdim0) || 
	    (vdim1 > 0 && cov->vdim2[1] != vdim1)) {
	  sprintf(ERRORSTRING,
		  "multivariate dimension (of submodel '%s'), which is %d x %d, does not match the specification of '%s', which is %d x %d and is required by '%s'",
		NICK(cov), cov->vdim2[0], cov->vdim2[1], C->name, vdim0, vdim1, 
		cov->calling == NULL ? "-- none --" :  P->name); 
	  return ERRORWRONGVDIM; // needed as value!
	  checkerror = true;
	}
	
	break;
      }
  
      if (err == ERRORINCORRECTSTATISO) {
	if (strcmp("", msg) != 0) {
	  err = ERRORM; 
	  strcpy(ERRORSTRING, msg);
	  strcpy(ERROR_LOC, "");
	}
	continue;
      } else if (err > NOERROR) {

       /// printf("err %d\n", err);
	//PMI(cov);

	errorMSG(err, msg); 
      }
      
    } // for iso

    if (err == NOERROR) break;
  } // dom

  //  printf("ok\n");
 
  if (PL > PL_COV_STRUCTURE && cov->calling == NULL) {
    LPRINT("%s: end look ", Nick(cov)); 
    if (err != NOERROR) PRINTF("err = %d\n", err); else MERR(err);
  }

  if (err != NOERROR) return err;


  if (PL >= PL_COV_STRUCTURE) {
    LPRINT("Continuing `%s' (no error):\n", Nick(cov)); 
  }

  if (!skipchecks && (err = check_within_range(cov, NAOK_RANGE)) != NOERROR) { 
    return err;
 }
 
  if (isoprev == SPACEISOTROPIC) {
    //print("\n\n\n");
    cov_model *cv = cov;
    while(cv->calling != NULL) cv = cv->calling;
    
    if (cov->xdimprev != 2) {
      return ERRORDIM;
    }
    if (cov->tsdim < 2)  {
      return ERRORDIM;
    }
    // cov->pref[TBM2] = PREF_NONE;
  }

  //  PMI(cov);
 
  //    printf("gattered %s %d %d %d %d\n", NAME(cov), domprev, cov->domown, newisoprev, cov->isoown);
  
 
  err = SetGatter(domprev, cov->domown, 
		  newisoprev, cov->isoown, 
		  nr, &UnUsedDeleteFlag); //

  //printf("setgatter  err = %d\n", err);
  ASSERT_GATTERONLY(cov);

  if (PL > PL_COV_STRUCTURE) {
    LPRINT("leaving '%s' for '%s'  SetGatter error=%d deriv=%d,%d \n",
	   Nick(cov),
	   cov->calling == NULL ? "none" : Nick(prev),
	   err, cov->full_derivs, cov->rese_derivs);
  }

  sprintf(ERROR_LOC, "%s: ", cov->calling == NULL ? "parsing the model"
	  : Nick(prev));

  // printf("end err = %d\n", err);
  COND_NEW_STORAGE(S2, GATTER, gatter_storage, z);

  // printf("2err = %d %s <= %s\n", err, NICK(cov), cov->calling == NULL ? "----" : NICK(cov->calling));
  cov->checked = err == NOERROR;
  assert(err == NOERROR || (cov->vdim2[0] > 0 && cov->vdim2[1] > 0));
  
  return(err);
   
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////




#define EARTH_LONGITUDE 0
#define EARTH_LATITUDE 1
#define pi180 0.017453292519943295474
#define EARTH_TRAFO(X, x, raequ, rpol)			\
  Rcos = raequ * cos(x[EARTH_LATITUDE] * pi180);	\
  X[0] = Rcos * cos(x[EARTH_LONGITUDE] * pi180);	\
  X[1] = Rcos * sin(x[EARTH_LONGITUDE] * pi180);	\
  X[2] = rpol * sin(x[EARTH_LATITUDE] * pi180)		

//   ; printf("x=%f =%f %f %f\n", x[0], x[1], raequ, rpol)

#define earth2cartInner(raequ, rpol)					\
  assert(cov->xdimprev == 2 && cov->xdimgatter == 3);			\
  double Rcos, X[3], Y[3];			\
  EARTH_TRAFO(X, x, raequ, rpol);			\
  EARTH_TRAFO(Y, y, raequ, rpol)

#define earth2cartInnerStat(raequ, rpol)				\
  assert(cov->xdimprev == 2 && cov->xdimgatter == 3);			\
  double Rcos, X[3];							\
  EARTH_TRAFO(X, x, raequ, rpol)

#define radiuskm_aequ 6378.1
#define radiuskm_pol 6356.8
#define radiusmiles_aequ 3963.17
#define radiusmiles_pol 3949.93

void EarthKM2CartStat(double *x, cov_model *cov, double *v) {
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].cov(X, cov, v);// nicht gatternr
}
void logEarthKM2CartStat(double *x, cov_model *cov, double *v, double *sign) {
  earth2cartInnerStat(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].log(X, cov, v, sign);// nicht gatternr
}
void EarthKM2Cart(double *x, double *y, cov_model *cov, double *v) {
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  //  printf("earth : %4.4f %4.4f y=%4.4f %4.4f X=%4.4f %4.4f %4.4f  Y=%4.4f %4.4f %4.4f\n", x[0], x[1], y[0], y[1], X[0], X[1], X[2], Y[0], Y[1], Y[2]);
  CovList[cov->secondarygatternr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logEarthKM2Cart(double *x, double *y, cov_model *cov, double *v,
		     double *sign) {
  earth2cartInner(radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].nonstatlog(X, Y, cov, v, sign);// nicht gatternr
}

void EarthMiles2CartStat(double *x, cov_model *cov, double *v) {
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  //printf("KM x=%f =%f %f %f\n", x[0], x[1], radiuskm_aequ, radiuskm_pol);
  CovList[cov->secondarygatternr].cov(X, cov, v);// nicht gatternr
}
void logEarthMiles2CartStat(double *x, cov_model *cov, double *v, double *sign) {
  earth2cartInnerStat(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].log(X, cov, v, sign);// nicht gatternr
}
void EarthMiles2Cart(double *x, double *y, cov_model *cov, double *v) {
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].nonstat_cov(X, Y, cov, v);// nicht gatternr
}
void logEarthMiles2Cart(double *x, double *y, cov_model *cov, double *v,
		     double *sign) {
  earth2cartInner(radiusmiles_aequ, radiusmiles_pol);
  CovList[cov->secondarygatternr].nonstatlog(X, Y, cov, v, sign);// nicht gatternr
}


int checkEarth(cov_model *cov){
  // ACHTUNG! KEIN AUFRUF VON SUB[0] !

  if (cov->domprev == XONLY// 20.2.14: warum war es vorher cov->domown == XONLY?
      && isSymmetric(cov->isoprev)) {
    // PMI(cov->calling);
    return ERRORKERNEL;
  }

  return NOERROR;
}
