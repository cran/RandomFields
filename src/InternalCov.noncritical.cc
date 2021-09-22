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

 Copyright (C) 2005 -- 2017 Martin Schlather

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



#include "def.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "Processes.h"
#include "families.h"
#include "Coordinate_systems.h"


void kdefault(model *cov, int i, double v) {

  defn *C = DefList + COVNR; // nicht gatternr
  if (PisNULL(i)) {
    switch(C->kappatype[i]) {
    case REALSXP :
      PALLOC(i, 1, 1);
      P(i)[0] = v;
      break;
    case INTSXP :
      PALLOC(i, 1, 1); 
      if (v == NA_INTEGER) PINT(i)[0] = NA_INTEGER;
      else if (!R_FINITE(v)) { BUG }
      else if (v > MAXINT) { BUG}
      else if (v < -MAXINT) { BUG}
      else PINT(i)[0] = (int) v;
      break;
    case STRSXP :
      ERR2("parameter '%.50s' in '%.50s' is undefined.", KNAME(i), NAME(cov));
      break;
    case LISTOF + REALSXP :
      PRINTF("%.50s:%.50s (%d) unexpected list\n", NICK(cov), C->kappanames[i], i);
      BUG;
    default : 
      PRINTF("%.50s:%.50s (%d) is not defined\n", NICK(cov), C->kappanames[i], i);
      BUG;
    }
    cov->nrow[i] = cov->ncol[i] = 1;
  } else if (!GLOBAL_UTILS->basic.skipchecks) {    
    if (cov->nrow[i] != 1 || cov->ncol[i] != 1) { 
      LPRINT("%d %.50s %d nrow=%d, ncol=%d\n", 
	     COVNR, NAME(cov), i, cov->nrow[i], cov->ncol[i]);
      int j; for (j=0; j<cov->ncol[i] * cov->nrow[i]; j++) {
	LPRINT("%10g\n", P(i)[j]);
      }
      ERR2("parameter '%.50s' in '%.50s' is not scalar -- pls contact author.",
	   KNAME(i), NAME(cov));
    }
  }
}


void updatepref(model *cov, model *sub) {
  int i;
  for (i=0; i< Forbidden; i++) {
    if (i==Specific) continue;
    if (sub->pref[i] < cov->pref[i]) {
      cov->pref[i] = sub->pref[i];
    }
  }
}


void setdefault(model *cov, int vdim0, int vdim1) {
  // Achilles-Ferse: setdefault wird aufgerufen bevor die Varianten
  // der Modelle und Operatoren feststehen -- die Parameter unten muessen, falls
  // notwendig von Hand in den checks gesetzt werden.

  defn *C = DefList + COVNR;// nicht gatternr 
  system_type *def = DEF;
  int
    n = SYSTEMS(def),
    last = LASTSYSTEM(def);

  //  printf("'%.50s' ----------------\n", NAME(cov)); 
  
  cov->full_derivs = C->F_derivs;
  cov->rese_derivs = C->RS_derivs;
  cov->loggiven = (ext_bool) (C->log != ErrLogCov); 
 
  for (int s=0; s<n; s++) {
    set_last(OWN, s, last);
    set_maxdim(OWN, s, MAXDIM(def, s)); 
  }

  cov->logspeed = RF_NA;
  
  if ((C->vdim == PREVMODEL_DEP  ||  C->vdim == SUBMODEL_DEP)) {    
    VDIM0 = vdim0;
    VDIM1 = vdim1;
  } else {
    assert(C->vdim != MISMATCH); // ???
    VDIM0 = VDIM1 = C->vdim;     
  } 

  if (isnowPosDef(cov)) { // Achtung in vielen Faellen nicht entscheidbar
    //                       da 'now' noch nicht definiert
    for (int i=0; i<MAXMPPVDIM; i++) cov->mpp.maxheights[i] = 1.0;    
  }

  if (isIsotropic(def) && isIsotropic(OWN) && 
      isPosDef(OWNTYPE(0)) && isXonly(C->systems[0]))
    cov->logspeed = 0.0;
  
  cov->finiterange = C->finiterange;  
  cov->monotone=C->Monotone;
  cov->ptwise_definite=C->ptwise_definite;

  //   cov->diag = 
  //   cov->semiseparatelast = 
  //   cov->separatelast = isIsotropic(OWNISO(0)); // war false vor 27.9.14

  MEMCOPY(cov->pref, C->pref, sizeof(pref_shorttype)); 
  
  //for (Methods i=Nothing+1; i<Forbidden; i++) cov->pref[i] = PREF_NONE;

  cov->method = Forbidden; 

  cov->taylorN = C->TaylorN;
  cov->tailN = C->TailN;

  for (int i=0; i<cov->taylorN; i++) {
    cov->taylor[i][TaylorConst] = C->Taylor[i][TaylorConst];
    cov->taylor[i][TaylorPow] = C->Taylor[i][TaylorPow];
  }
  for (int i=0; i<cov->tailN; i++) {
    cov->tail[i][TaylorConst] = C->Tail[i][TaylorConst];
    cov->tail[i][TaylorPow] = C->Tail[i][TaylorPow];
    cov->tail[i][TaylorExpConst] = C->Tail[i][TaylorExpConst];
    cov->tail[i][TaylorExpPow] = C->Tail[i][TaylorExpPow];
  }
}


int merge_integer(int val, int sub){
  if (val == SUBMODEL_DEP) return sub;
  else {    
    //    if (*val == PARAM_DEP) {model *cov; crash(cov);}
    assert(val != PARAM_DEP);
    if (sub < val) return sub;
  }
  return val;
}


ext_bool merge_extbool(ext_bool val, ext_bool sub){
  if (val == SUBMODEL_DEP) return sub;
  else {
    assert((int) Paramdep < 0 && (int) Submodeldep < 0);
    if ((int) sub < (int) val) return sub;
  }
  return val;
}


monotone_type merge_monotone(monotone_type val, monotone_type sub){
  if (val == MON_SUB_DEP) return sub;
  else {
    if ((int) sub < (int) val) return sub;
  }
  return val;
}


// no setbackard
void setbackward(model *cov, model *sub) {
   defn *C = DefList + COVNR;// nicht gatternr
 // see also check_co when changing

  assert(cov != NULL);
  assert(sub != NULL);
 
  int last = OWNLASTSYSTEM;
  if (last == SUBLASTSYSTEM) {
    for (int s=0; s<=last; s++) { 
      //    set_maxdim(OWN, s, merge_integer(MAXDIM(OWN, s), MAXDIM(SUB, s)));
    }
  }
  cov->monotone = merge_monotone(cov->monotone, sub->monotone);  
  cov->finiterange = merge_extbool(cov->finiterange, sub->finiterange);
  
  if (sub->full_derivs < cov->full_derivs || cov->full_derivs == UNSET)
      cov->full_derivs = sub->full_derivs;

  assert(cov->full_derivs >= 0 || 
	 (cov->full_derivs == MISMATCH && isRandom(cov) && isRandom(sub)));
  if (sub->rese_derivs < cov->rese_derivs || cov->rese_derivs == UNSET)
    cov->rese_derivs = sub->rese_derivs;
  if (cov->loggiven != falsch) cov->loggiven = sub->loggiven;
  
  updatepref(cov, sub);

  if (sub==cov->sub[0] || sub==cov->key) {
    if (C->vdim == SUBMODEL_DEP) {
      VDIM0 = sub->vdim[0];
      VDIM1 = sub->vdim[1];
    }
    if (C->ptwise_definite == pt_submodeldep) {
      cov->ptwise_definite = sub->ptwise_definite;
      // assert((cov->ptwise_definite != pt_paramdep &&  // to do !
      //	      cov->ptwise_definite != pt_submodeldep &&
      //	      cov->ptwise_definite != pt_undefined));
   }
  } else {
    if (cov->ptwise_definite != sub->ptwise_definite) {
       // assert((cov->ptwise_definite != pt_paramdep && // to do !
       //     cov->ptwise_definite != pt_submodeldep &&
       ///     cov->ptwise_definite != pt_undefined &&
       //     sub->ptwise_definite != pt_paramdep &&
       //     sub->ptwise_definite != pt_submodeldep &&
       //     sub->ptwise_definite != pt_undefined));
      if (cov->ptwise_definite ==pt_mismatch ||
	  sub->ptwise_definite==pt_mismatch)
	cov->ptwise_definite = pt_mismatch;
      else if (cov->ptwise_definite==pt_unknown ||
	       sub->ptwise_definite==pt_unknown)
	cov->ptwise_definite = pt_unknown;
      else 	  
	cov->ptwise_definite = 
	  cov->ptwise_definite == pt_zero ? sub->ptwise_definite
	  : sub->ptwise_definite == pt_zero ? cov->ptwise_definite      
	  : pt_indef;
    }
  }

  cov->hess = (DefList[COVNR].hess != NULL && sub->hess);
  cov->randomkappa |= sub->randomkappa;
}


int alloc_mpp_M(model *cov, int moments) {
  int maxmoments = DefList[COVNR].maxmoments;
  // insbesondere fuer models die selbst vom Random-Type sind
  assert(moments >= 0);
  assert(maxmoments != MISMATCH && maxmoments != PARAM_DEP);


  if (moments > maxmoments && maxmoments != SUBMODEL_DEP) {
    SERR2("required moments (%d) exceeds the coded moments (%d)",
	  moments, maxmoments);
  }
  if (moments <= cov->mpp.moments) RETURN_NOERROR;
  if (cov->mpp.mM != NULL) free_mpp_M(cov);
  cov->mpp.moments = moments;

  int i,
    vdim = VDIM0,
    nm = cov->mpp.moments,
    nmvdim = (nm + 1) * vdim,
    bytes = sizeof(double) * nmvdim;

  if (vdim <= 0) BUG;
  if (vdim > MAXMPPVDIM) SERR1("multivariate dimension (%d) too large", vdim);
 
  cov->mpp.mM = (double*) MALLOC(bytes);
  cov->mpp.mMplus = (double*) MALLOC(bytes);

  //  assert(nm < 100);
  int nmP1 = cov->mpp.moments + 1;
  for (i=0; i<nmvdim; i++) cov->mpp.mMplus[i] = cov->mpp.mM[i] = RF_NA; 
  for (i=0; i<vdim; i++) {
    int idx = i * nmP1;
    cov->mpp.mMplus[idx + 0] = cov->mpp.mM[idx + 0] = RF_INF;
    cov->mpp.maxheights[i] = RF_NAN;
  }

  // cov->mpp.mMplus[0] = cov->mpp.mM[0] = 1.0;
  RETURN_NOERROR;
}

void free_mpp_M(model *cov) {
  FREE(cov->mpp.mM);
  FREE(cov->mpp.mMplus);
  cov->mpp.mM = cov->mpp.mMplus = NULL;
}

int UpdateMPPprev(model * cov, int moments) {
  model *calling = cov->calling;
  int i, nm, err,
    nmvdim,
    vdim = VDIM0;

  nm = cov->mpp.moments < calling->mpp.moments ? cov->mpp.moments
	: calling->mpp.moments;
  nmvdim = (nm + 1) * vdim;
  if (moments >= 0 && calling != NULL) {
    if (calling->mpp.moments == SUBMODEL_DEP &&
	(err = alloc_mpp_M(calling, moments)) != NOERROR) RETURN_ERR(err);
    for (i=0; i<nmvdim; i++) {
      calling->mpp.mMplus[i] = cov->mpp.mMplus[i];
      calling->mpp.mM[i]     = cov->mpp.mM[i];
    }
  }
  
  // nachfolgende Zeilen so lassen, da sonst unerwuenscht
  // maxheight etc. nach oben gegeben werden.
  // calling->mpp.maxheight = cov->mpp.maxheight; 
  // calling->mpp.unnormedmass = cov->mpp.unnormedmass;
  
  RETURN_NOERROR;
}


int INIT_intern(model *cov, int moments, gen_storage *s) { // kein err  
  if (!cov->checked) BUG;
  if (cov->initialised) RETURN_NOERROR;
  assert(cov != NULL);
  ASSERT_GATTER(cov);
  

  defn *C = DefList + COVNR;
  int err = NOERROR;
  KEY_type *KT = cov->base;

  SPRINTF(KT->error_loc, "initializing %.50s", NICK(cov));

  // printf("Iintern %.50s = %d mom=%d, vdim=%d\n", NAME(cov), cov->mpp.moments, moments, VDIM0);
  
  if (cov->mpp.moments != SUBMODEL_DEP && cov->mpp.moments != PARAM_DEP) {
    if ((err = alloc_mpp_M(cov, moments)) != NOERROR) RETURN_ERR(err);
  } else {
    BUG; //assert(false); // passt das hier??
    if (cov->mpp.moments == PARAM_DEP) cov->mpp.moments = moments;
  }
 
  if (C->maxmoments >= 0 && moments > C->maxmoments) {
      SERR3("moments known up to order %d for '%.50s', but order %d required",
	    C->maxmoments, NICK(cov), moments);
  }
  
  SPRINTF(KT->error_loc, "%.50s", cov->calling == NULL ? "initiating the model"
	  : NICK(cov->calling));
  ASSERT_GATTER(cov);
  if ((err = DefList[GATTERNR].Init(cov, s)) != NOERROR) {
    RETURN_ERR(err);   
  }
   
  if ((err = UpdateMPPprev(cov, moments)) != NOERROR) {
    RETURN_ERR(err);
  }

  cov->initialised = true;
  RETURN_NOERROR;
}

void set_initialised_false(model *cov){
  int i;
  if (!cov->randomkappa) return;
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


int REINIT_intern(model *cov, int moments, gen_storage *s) { // kein err
  int err;
  set_initialised_false(cov);
  err = INIT_intern(cov, moments, s);
  RETURN_ERR(err);
}


int INIT_RANDOM_intern(model *cov, int moments, gen_storage *s, // kein err
		       double *p) {
  if (!cov->checked) BUG;
  if (!cov->initialised) {
    int err = NOERROR;
    KEY_type *KT = cov->base;

    SPRINTF(KT->error_loc, "initializing %.50s", NICK(cov));
     
    assert(cov != NULL);
    if (moments < 0) SERR("moments expected to be positive");
    if (DefList[COVNR].maxmoments >= 0 &&
	moments > DefList[COVNR].maxmoments) SERR("Moments do not match");
    
    if (cov->mpp.moments != SUBMODEL_DEP && cov->mpp.moments != PARAM_DEP) {
      if ((err = alloc_mpp_M(cov, moments)) != NOERROR) RETURN_ERR(err);
    } else {
      BUG; // passt das hier??
      if (cov->mpp.moments == PARAM_DEP) cov->mpp.moments = moments;
    }
    
    SPRINTF(KT->error_loc, "%.50s", cov->calling == NULL ? "initiating the model"
	    : NICK(cov->calling));
    ASSERT_GATTER(cov);
    if ((err = DefList[GATTERNR].Init(cov, s)) != NOERROR) RETURN_ERR(err);   
    if (ISNAN(cov->mpp.mM[moments])) {
      SERR1("%.50s is not a random function", NICK(cov));
    }
    
    if ((err = UpdateMPPprev(cov, moments)) != NOERROR) RETURN_ERR(err);
    
    cov->initialised = true;
  }

  //  switch (DefList[COVNR].kappatype[param_nr]) {
  //  case REALSXP :   
  if (s->dosimulate) DORANDOM(cov, p);
    

    //    break;
    //  case INTSXP :
    //    int j, len;
    //    double *dummy;
    //    dummy = (double*) MALLOC(sizeof(double) * len);
    //    DORANDOM(cov, dummy);
    //    for (j=0; j<len; j++) p[j] = (int) dummy[j];
    //    FREE(dummy);
    //    break;
    //  default : SERR("random parameter only allowed for numerical values");
    //  }
  
  RETURN_NOERROR;
}


void stat2_Intern(double *x, model *cov, double **Z) {
  //    PMI(cov->calling);
  double *z = *Z;
  TALLOC_DOUBLE(z1);
  bool trafo = cov->calling != NULL && TRAFONR != UNSET;
  if (trafo) {
    TALLOC_GATTER_GLOBAL(z1, GATTERTOTALXDIM);
    DefList[TRAFONR].cov(x, cov, z1);
    x = z1;
  }

  int n = GATTERLASTSYSTEM;
  //  PMI0(cov);  printf("n=%d\n", n);
  assert(n == 0); // **Z is a workaround for more complicated systems (i_nrow)
  for (int s=0; s<=n; s++) {
    int gnr = NRi(GATTER[s]), // OK
      dim = GATTERXDIM(s);
     //    printf("stat2 %.50s\n", NAME(cov));
    assert(dim > 0 && gnr >= FIRSTGATTER && gnr <= LASTGATTER);
    
    switch(gnr) {
    case ISO2ISO :
      *x=FABS(*x);
      if (trafo) *z = *x; else *Z = x;
      break; // crucuial for i_row, r_col
    case S2ISO : case SP2ISO : {
      double b = 0.0;
      for (int d=0; d<dim; d++) b += x[d] * x[d];
      *z = SQRT(b);
      break;
    }
    case S2S : case SId :
      //      printf("dim = %d %.50s %d\n", dim, NAME(cov), trafo);
      if (trafo) MEMCOPY(z, x, sizeof(double) * dim); else *Z = x;
     break;
    case SP2SP :
      x[0]=FABS(x[0]); x[1]=FABS(x[1]);
      if (trafo) MEMCOPY(z, x, sizeof(double) * dim); else *Z = x;
      break;
    case S2SP : {
      double b = 0.0;
      int dimM1 = dim - 1;
      for (int d=0; d<dimM1; d++) b += x[d] * x[d];
      z[0] = SQRT(b);
      z[1] = FABS(x[dimM1]);
      break;
    }      
    case E2EIso : EarthIso2EarthIso(x);
      if (trafo) *z = *x; else *Z = x;
      break;
    case E2E : Earth2Earth(x);
      if (trafo) MEMCOPY(z, x, sizeof(double) * dim); else *Z = x;
      break;
    case E2SphIso : EarthIso2SphereIso(x, cov, z); break;
    case E2Sph : Earth2Sphere(x, cov, z); break;
    case Sph2SphIso : SphereIso2SphereIso(x);
      if (trafo) *z = *x; else *Z = x;
      break;    
    case Sph2Sph : Sphere2Sphere(x);
      if (trafo) MEMCOPY(z, x, sizeof(double) * dim); else *Z = x;
      break;  
    default :
      //      printf("gnr=%d\n", gnr);
      //      PMI(cov);
      BUG;
    }
    z += OWNXDIM(s);
    x += dim;
  }
  FREE_TALLOC(z1);
}

void stat2(double *x, model *cov, double *v) {
  TALLOC_GATTER(z, OWNTOTALXDIM);
  double **z1 = &z;
  stat2_Intern(x, cov, z1);
  //  if (OWNTOTALXDIM == 2)
  //printf("z after =%ld %10g %10g %10g; x=%ld %d xdim=%d\n", (*z1),(*z1)[0],(*z1)[1],
  //	   (*z1)[2], x, VDIM0, OWNTOTALXDIM);
  //  printf("stat2: %.50s\n", NAME(cov));
  DefList[COVNR].cov(*z1, cov, v);// nicht gatternr
  END_TALLOC_z; // very crucial that not z is freed!
}

void logstat2(double *x, model *cov, double *v, double *Sign) {
  TALLOC_GATTER(z, OWNTOTALXDIM);
  double **z1 = &z;
  stat2_Intern(x, cov, z1);
  DefList[COVNR].log(*z1, cov, v, Sign);// nicht gatternr
  END_TALLOC_z;
}


void nonstat2stat(double *x, double *y, model *cov, double *z) {
  int n = GATTERLASTSYSTEM;
  for (int s=0; s<=n; s++) {
    int gnr = NRi(GATTER[s]), // OK
      dim = GATTERXDIM(s);

    //    printf("nonstat2 %.50s\n", NAME(cov));
    assert(dim > 0 && gnr >= FIRSTGATTER && gnr <= LASTGATTER);
    switch(gnr) {
    case S2ISO : {
      double b = 0.0;
      for (int d=0; d<dim; d++) {
	//	printf("d=%d  %10g %10g %d\n", d, x[d], y[d], dim);
	//	if (dim == 1) APMI(cov->calling);
	
	double a = x[d] - y[d];
	b += a * a;
      }
      *z = SQRT(b);
      break;
    }
    case S2S :case SId:
      for (int d=0; d<dim; d++) z[d] = x[d] - y[d];
      break;
    case S2SP :{
      double b = 0.0;
      int dimM1 = dim - 1;
      for (int d=0; d<dimM1; d++) {
	double a = x[d] - y[d];
	b += a * a;
      }
      z[0] = SQRT(b);
      z[1] = FABS(x[dimM1] - y[dimM1]);
      break;
    }    
    case E2EIso : NonstatEarth2EarthIso(x, y, cov, z); break;
      //    case E2E : NonstatEarth2Earth(x, y, cov, z); break;
    case E2SphIso : NonstatEarth2SphereIso(x, y, cov, z); break;
      //    case E2Sph : NonstatEarth2Sphere(x, y, cov, z); break;
    case Sph2SphIso : NonstatSphere2SphereIso(x, y, cov, z); break;     
      //    case Sph2Sph : NonstatSphere2Sphere(x, y, cov, z); break;  
    default : //e.g. SId, E2E, ...
      // S2FORBIDDEN :
      //      printf("gnr=%d\n", gnr);
      //      PMI(cov);
      BUG;
    }
    z += OWNXDIM(s);
    x += dim;
    y += dim;
  }

}

void nonstat2(double *x, double *y, model *cov, double *v) {
  TALLOC_DOUBLE(z1);
  TALLOC_DOUBLE(z2);
  if (cov->calling != NULL && TRAFONR != UNSET) {
    TALLOC_GATTER_GLOBAL(z1, GATTERTOTALXDIM);
    TALLOC_GATTER_GLOBAL(z2, GATTERTOTALXDIM);
    DefList[TRAFONR].cov(x, cov, z1);
    DefList[TRAFONR].cov(y, cov, z2);
    x = z1;
    y = z2;
  }
  if (equalsKernel(OWNDOM(0))) {
    assert(DefList[COVNR].nonstat_cov != nonstat2);
    DefList[COVNR].nonstat_cov(x, y, cov, v);// nicht gatternr
  } else {
    //  printf("> %.50s ", NAME(cov));
    // kernel2xonly:
    TALLOC_GATTER(z, OWNTOTALXDIM);
    nonstat2stat(x, y, cov, z);
    DefList[COVNR].cov(z, cov, v);// nicht gatternr 
    END_TALLOC_z;
  }
  FREE_TALLOC(z1);
  FREE_TALLOC(z2);
}

void lognonstat2(double *x, double *y, model *cov,
		       double *v, double *Sign) {
  TALLOC_DOUBLE(z1);
  TALLOC_DOUBLE(z2);
  if (cov->calling != NULL && TRAFONR != UNSET) {
    TALLOC_GATTER_GLOBAL(z1, GATTERTOTALXDIM);
    TALLOC_GATTER_GLOBAL(z2, GATTERTOTALXDIM);
    DefList[TRAFONR].cov(x, cov, z1);
    DefList[TRAFONR].cov(y, cov, z2);
    x = z1;
    y = z2;
  }
  if (equalsKernel(OWNDOM(0))) {
    DefList[COVNR].nonstatlog(x, y, cov, v, Sign);// nicht gatternr
    return;
  } else {
    // kernel2xonly:
    TALLOC_GATTER(z, OWNTOTALXDIM);
    nonstat2stat(x, y, cov, z);
    DefList[COVNR].log(z, cov, v, Sign);// nicht gatternr
    END_TALLOC_z;
  }
  FREE_TALLOC(z1);
  FREE_TALLOC(z2);
}


void D_2(double *x, model *cov, double *v){  
  assert(everyCoord(isSpaceIsotropic, cov));
  defn *C = DefList + COVNR;// nicht gatternr
  int dim = GATTERXDIM(0); // frueher prevxdim
  if (dim == 1) {
    assert(equalsIsotropic(OWNISO(0)));
    double y = FABS(*x);

    // iso2iso ueberpruefen !!!
    // dollar + scale // aniso > 0 und dim = 1
    // nach tbm eigentlich auch nicht

    C->D(&y, cov, v);// nicht gatternr
  } else {
    assert(dim == 2);
    if (OWNTOTALXDIM == 1) {
      assert(equalsIsotropic(OWNISO(0)));
      double y=SQRT(x[0] * x[0] + x[1] * x[1]);    
      C->D(&y, cov, v);
      if (y!=0.0) *v *= x[0] / y;
    } else {
      assert(OWNTOTALXDIM == 2);
      double y[2];
      y[0] = FABS(x[0]);
      y[1] = FABS(x[1]);
      C->D(y, cov, v); 
    }
  }
}
void DD_2(double *x, model *cov, double *v) {
  assert(everyCoord(isSpaceIsotropic, cov));
  defn *C = DefList + COVNR;// nicht gatternr
  assert(everyCoord(isIsotropic, cov));
  int dim = GATTERXDIM(0); // frueher prevxdim
  if (dim == 1) {
    double y = FABS(*x);

    // iso2iso ueberpruefen !!!
    // dollar + scale // aniso > 0 und dim = 1
    // nach tbm eigentlich auch nicht
    C->D2(&y, cov, v);// nicht gatternr
  } else {
    //   assert(PREVTOTALXDIM == 2);
    assert(dim == 2);
    system_type *def = DEF;
    if (isIsotropic(def)) {
      double w,
	xSq = x[0] * x[0],
	tSq = x[1] * x[1],
	ySq = xSq + tSq,
	y   = SQRT(ySq); 
      
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
    } else if (equalsSpaceIsotropic(def)) {
      double y[2];
      y[0] = FABS(x[0]);
      y[1] = FABS(x[1]);
      C->D2(y, cov, v); // nicht gatternr
    } else BUG;
  }
}

void DD_3(double VARIABLE_IS_NOT_USED *x, model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *v) {
  assert(everyCoord(isSpaceIsotropic, cov));
  ERR("DD_3 to be programmed\n");
}


void inverse2(double *x, model *cov, double *v) {
  defn *C = DefList + COVNR;// nicht gatternr
  C->inverse(x, cov, v);//  nicht gatternr
}
   

void nonstatinverse2(double *v, model *cov, double *x, double *y){
  defn *C = DefList + COVNR;// nicht gatternr

  C->nonstat_inverse(v, cov, x, y);//  nicht gatternr
}

void nonstat_loginverse2(double *v, model *cov, double *x, double *y){
  defn *C = DefList + COVNR;// nicht gatternr

  C->nonstat_loginverse(v, cov, x, y);//  nicht gatternr
}
   
   

int struct2(model *cov, model **newmodel) {
  int err;
  errorloc_type errloc_save;
  KEY_type *KT = cov->base;

  if (!cov->checked) {
    BUG;
  }
  STRCPY(errloc_save, KT->error_loc);
  SPRINTF(KT->error_loc, "setting up %.50s", NICK(cov));
 
  err = DefList[COVNR].Struct(cov, newmodel);
   if (newmodel != NULL && (*newmodel) != NULL) {
    SET_CALLING(*newmodel, cov->calling != NULL ? cov->calling : cov);
  }
 
  if (err == NOERROR) STRCPY(KT->error_loc, errloc_save);

  RETURN_ERR(err);
}

int init2(model *cov, gen_storage *s){ // s wird durchgereicht!

  defn *C = DefList + COVNR; //  nicht gatternr
  model
    *calling = cov->calling == NULL ? cov : cov->calling;
  int i,
    err = NOERROR,
    kappas = DefList[COVNR].kappas;
  errorloc_type errloc_save;
  KEY_type *KT = cov->base;
  STRCPY(errloc_save, KT->error_loc);
   
  //  PrInL++;

  for (i=0; i<kappas; i++) {
    model *param  = cov->kappasub[i];
    if (param != NULL) {
      //      PMI0(param);
      if (isnowRandom(param)) {
	if ((err = INIT_RANDOM(param, 0, s, P(i))) != NOERROR) RETURN_ERR(err);
      } else if (!isnowShape(param) && (err = INIT(param, 0, s)) != NOERROR)
	RETURN_ERR(err);
    }
  }

  if (cov->method == Forbidden) { cov->method = calling->method; }
  
  SPRINTF(KT->error_loc, "Initializing %.50s", NICK(cov));
  if (!equalsBernoulliProcess(cov)) {
    switch(cov->frame) {
    case GaussMethodType :
      if (cov->method==SpectralTBM) {
	if (cov->calling == NULL && COVNR != SPECTRAL_PROC_USER &&
	    COVNR != SPECTRAL_PROC_INTERN) SERR("unexpected value in init2");
      }
      break;
    case InterfaceType: break;
    case BrMethodType: case SmithType: case SchlatherType :
    case PoissonType: case PoissonGaussType: case RandomType:
      cov->origrf = false;
      assert((cov->mpp.moments < 0) xor (cov->mpp.mM != NULL));
      break;
    case TrendType: case EvaluationType: case LikelihoodType: break;
    case NormedProcessType: break;
    default :
      // PMI0(cov);      crash();
      ILLEGAL_FRAME;
    }
  }

  if (!cov->initialised && (err = C->Init(cov, s)) != NOERROR)
    goto ErrorHandling;
  
  calling->fieldreturn = cov->fieldreturn;
  
 ErrorHandling :
  //  PrInL--;
  if (err == NOERROR) STRCPY(KT->error_loc, errloc_save);
  cov->initialised = err == NOERROR;
  SPRINTF(KT->error_loc, "'%.50s'", NICK(calling));//  nicht gatternr   
  RETURN_ERR(err);
}

void do2(model *cov, gen_storage *s){
  //   model *calling = cov->calling == NULL ? cov : cov->calling;
  //
  //  int i,
  //    kappas = DefList[COVNR].kappas;

  // statt nachfolgende Zeilen: siehe init2
  //  for (i=0; i<kappas; i++) {
  //    model *param  = cov->kappasub[i];
  //    if (param != NULL && isnowRandom(param)) DORANDOM(param, P(i));
  //  }

  DefList[COVNR].Do(cov, s); // ok

  // assert(false);
}


void dorandom2(model *cov, double *v){
 DefList[COVNR].DoRandom(cov, v); // ok
}





/*
void checkfor(model *cov) {// served by CHECKFOR; nothing to do with CHECK & err
  while (cov!=NULL && !isProcess(cov)) {
    PRINTF("%.50s ", NAME(cov));
    cov=cov->calling;
  }
  if (cov == NULL) { PRINTF("not process(es)\n"); return; }
  if (OWNISO(0) != PREVISO(0) && OWNISO(0) != ISO_MISMATCH) { APMI(cov); } //
}
*/


sortsofparam SortOf(model *cov,int k, int row, int col, sort_origin origin) {
  defn *C = DefList + COVNR;
  if (C->sortof != NULL) return C->sortof(cov, k, row, col, origin);
  // for MLE
// k non-negative: k-th parameter
// k negative: (i-1)th submodel
  if (k >= C->kappas) BUG;

  return k<0 ? VARPARAM : C->sortof_tab[k];
// k<0: varparam used to indicate that variance is allowed for submodel,
// see recursive call getnaposition; E.g. not allowed for second submodel
// of nsst: the variance parameter is something between scale and variance
}


void assert_sys(system_type VARIABLE_IS_NOT_USED * sys){assert(sys != NULL);}
//void assert_cov(model VARIABLE_IS_NOT_USED *cov) {assert(cov != NULL || __extension__({crash(); false;}));}
//
void assert_cov(model VARIABLE_IS_NOT_USED *cov) {
  //if (cov==NULL) crash();
  assert(cov != NULL);
}



bool allowedDfalse(model *cov) {
  for(int i=FIRST_DOMAIN; i<= LAST_DOMAINUSER; i++) cov->allowedD[i] = true;
  return false;
}


bool allowedDtrue(model *cov) {
  for(int i=FIRST_DOMAIN; i<=LAST_DOMAINUSER; i++) cov->allowedD[i] = true;
  return true;
}


bool allowedD(model *cov) {
  // if (cov->DallowedDone) crash();
  assert(!cov->DallowedDone);
  // muss zwingend als erstes stehen
  cov->DallowedDone = cov->calling == NULL ? true : cov->calling->DallowedDone;
  //                                  da ParamDepD fehlschlagen kann
  //printf("Dsub: %.50s\n", NAME(cov));
  assert((int) FIRST_DOMAIN == 0);
  // keine Varianten bei DOMAIN
  defn *C = DefList + COVNR;
  bool *a = cov->allowedD;

  cov->variant=0;
  assert(!isSubModelD(C) || C->Dallowed != NULL);
  if ((C->Dallowed != NULL)) {
    //printf("subi\n");
    return C->Dallowed(cov);
  }

  domain_type dom = DEFDOM(0);
  if (isParamDepD(C) && C->setDI !=NULL && !isFixed(dom) &&
      !C->setDI(cov)){
    cov->DallowedDone = false;
    return  allowedDfalse(cov);
  }
  if (isFixed(dom)) {
    for (int i=(int) FIRST_DOMAIN; i<=(int) LAST_DOMAINUSER; a[i++] = false);
    a[dom] = true;    
    return false;
  }
  return allowedDfalse(cov);
}


bool allowedIfalse(model *cov) {
  //printf("allowedfalse %.50s\n", NAME(cov));
  bool *I = cov->allowedI;
  for(int i=FIRST_ISOUSER; i<= LAST_ISOUSER; i++) I[i] = true;
  return false;
}


bool allowedItrue(model *cov) {
  //printf("allowedtrue %.50s\n", NAME(cov));
  bool *I = cov->allowedI;
  for(int i=FIRST_ISOUSER; i<=LAST_ISOUSER; i++) I[i] = true;
  return true;
}


bool allowedI(model *cov) {
  if (cov->IallowedDone) return false;

  //printf("allowedI:%.50s %d\n", NAME(cov), cov->zaehler);
  cov->IallowedDone = cov->calling == NULL ? true : cov->calling->IallowedDone;
  //                                  da ParamDepI fehlschlagen kann
  //  printf("Isub: %.50s (%d)\n", NAME(cov), cov->zaehler);
  assert((int) FIRST_ISOUSER == 0);
  defn *C = DefList + COVNR;
  int variants = C->variants;
  bool *I = cov->allowedI;

  cov->variant=0;
  assert(!isSubModelI(C) || C->Iallowed != NULL);
  if (C->Iallowed != NULL) return C->Iallowed(cov);
  
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  
  isotropy_type iso = DEFISO(0);
  if (isParamDepI(C) && C->setDI!=NULL && !isFixed(iso) && !C->setDI(cov)){
    cov->IallowedDone = false;
    return allowedIfalse(cov);
  }

  if (isFixed(iso)) {
    //    printf("fixed %.50s\n", ISO_NAMES[iso]);
    I[iso] = true;
    if (equalsUnreduced(iso)) {
      //I[SYMMETRIC] = I[EARTH_SYMMETRIC] = I[SPHERICAL_SYMMETRIC] = false;
      I[CARTESIAN_COORD] =I[EARTH_COORD] = I[SPHERICAL_COORD] = true;
    }
  } else {
    return allowedIfalse(cov);
  }

  for (cov->variant++ ; cov->variant < variants; cov->variant++) {
    // hier nur noch fixe werte zulaessig!
    //   printf("ABBx %d\n", cov->variant);
    iso = DEFISO(0);
    assert(!isParamDepI(C) && !isSubModelI(C));
    assert(isFixed(iso));
    I[iso] = true;  
  }
  
  //  printf("ABB\n");
   cov->variant=0;

  return false;
}


bool allowedIsubs(model *cov, model **sub, int z){
  //  printf("allowed check for %.50s\n", NAME(cov));
  
  bool *I = cov->allowedI;
  int
    j=0,
    idx_C = FIRST_CARTESIAN,
    idx_E = FIRST_EARTH,
    idx_S = FIRST_SPHERICAL;

  if (z == 0) {
    return allowedItrue(cov);
  }

  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  //  printf("J=%d %d\n", j, sub[1]->IallowedDone);
  while (true) {
    if (!allowedI(sub[j])) break;
    if (++j >= z) return allowedItrue(cov);
  }
  MEMCOPY(I, sub[j]->allowedI, sizeof(allowedI_type));
  //for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; i++)  printf("%d \n", I[i]);
  //printf("%.50s j=%d\n", NAME(sub[j]), j);  printI(sub[j]);
      
  while (idx_C <= (int) CARTESIAN_COORD && !I[idx_C]) idx_C++;
  while (idx_E <= (int) EARTH_COORD && !I[idx_E]) idx_E++;
  while (idx_S <= (int) SPHERICAL_COORD && !I[idx_S]) idx_S++;
  //PMI(sub[j]);
  assert(idx_C <= (int) CARTESIAN_COORD || idx_E <= (int) EARTH_COORD ||
	 idx_S <= (int) SPHERICAL_COORD);
  if (idx_C > CARTESIAN_COORD) {
    if (idx_E > EARTH_COORD) {
      idx_E = idx_S-FIRST_SPHERICAL+FIRST_EARTH;
      for (int i=idx_S; i<=SPHERICAL_COORD; i++)
	I[i-FIRST_SPHERICAL+FIRST_EARTH] = I[i];
    }
  }
    
  for (j++; j<z; j++) {
    assert (sub[j] != NULL);
    // printf("J=%d %d\n", j, sub[1]->IallowedDone);
    if (allowedI(sub[j])) continue;
    bool *subI = sub[j]->allowedI;
    int sub_C = FIRST_CARTESIAN;
    while (sub_C <= (int) CARTESIAN_COORD && !subI[sub_C]) sub_C++;
    int sub_E = FIRST_EARTH;
    while (sub_E <= (int) EARTH_COORD && !subI[sub_E]) sub_E++;
    int sub_S = FIRST_SPHERICAL;
    while (sub_S <= (int) SPHERICAL_COORD && !subI[sub_S]) sub_S++;
 
    if (sub_C <= (int) CARTESIAN_COORD) {
      if (idx_C <= (int) CARTESIAN_COORD) {
	if (equalsVectorIsotropic((isotropy_type) idx_C) && sub_C < idx_C)
	  idx_C = sub_C;
	for (; idx_C < sub_C; I[idx_C++] = false);
	for (int i = idx_C; i<=(int) CARTESIAN_COORD; i++) I[i] |= subI[i];
      }
    } else {
      idx_C = sub_C;
    }
    if (sub_E > EARTH_COORD) {
      if (sub_S <= SPHERICAL_COORD) {
	sub_E = sub_S-FIRST_SPHERICAL+FIRST_EARTH;
	for (int i=sub_S; i<=SPHERICAL_COORD; i++)
	  subI[i-FIRST_SPHERICAL+FIRST_EARTH] = subI[i];
      } else idx_S = sub_S;
    }

    if (sub_S <= (int) SPHERICAL_COORD) {
      if (idx_S <= (int) SPHERICAL_COORD) {
	for (; idx_S < sub_S; I[idx_S++] = false);
	for (int i = idx_S; i<=(int) SPHERICAL_COORD; i++) I[i] |= subI[i];
      } // else {
      //for (; sub_S <=(int) SPHERICAL_COORD; sub_S++) I[sub_S] = subI[sub_S];
      //}
    }
      
    if (sub_E <= (int) EARTH_COORD) {
      if (idx_E <= (int) EARTH_COORD) {
	for (; idx_E < sub_E; I[idx_E++] = false);
	for (int i = idx_E; i<=(int) EARTH_COORD; i++) I[i] |= subI[i];
      } //else {
      //for ( ; sub_E <=(int) EARTH_COORD; sub_E++) I[sub_E] = subI[sub_E];
      // }
    }
    if (idx_C >= (int) CARTESIAN_COORD && idx_E >= (int) EARTH_COORD &&
	idx_S >= (int) SPHERICAL_COORD) break;
  }
  if (idx_C <= CARTESIAN_COORD){
    if (equalsSymmetric((isotropy_type) idx_C)) I[EARTH_SYMMETRIC] = true;
    else I[EARTH_COORD] = I[EARTH_SYMMETRIC] = true;
  }
  if (idx_S <= SPHERICAL_COORD) {
    if (equalsSphericalCoord((isotropy_type) idx_S) &&
	isSymmetric((isotropy_type) idx_C)) // !! not equalsSymmetric!!
      idx_S = SPHERICAL_SYMMETRIC;    
    for (int i=idx_S; i<=SPHERICAL_COORD; i++)
      I[i-FIRST_SPHERICAL+FIRST_EARTH] = I[i] = true;
  } else if (idx_C > CARTESIAN_COORD) {
    if (idx_E >= EARTH_COORD) idx_E = EARTH_SYMMETRIC;
    for (int i=idx_E; i<=EARTH_COORD; i++) I[i] = true;
  }

  //  printf("%.50s:", NAME(cov));
  //  for (int i=0; i<= EARTH_COORD; i++) {printf("%d.", I[i]);} printf("\n");
  
  return false;
}

