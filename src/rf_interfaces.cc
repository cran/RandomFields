/*
Authors
Martin Schlather, schlather@math.uni-mannheim.de

main library for unconditional simulation of random fields

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
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
along with this program; if not, write to the Fre Software
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

#include "def.h"
#ifdef DO_PARALLEL
#include <omp.h>
#endif
#include <Basic_utils.h>

#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "extern.h"
#include "questions.h"
#include "Coordinate_systems.h"
#include "variogramAndCo.h"
#include "Processes.h"
#include "rf_interfaces.h"
extern const char *general[generalN];


/* 
 in CheckCovariance and other the following dimensions are used:
 xdimOZ       : dimension of the points given by the user.
              value is <= spacedim since in case of isotropy only
	      the distances might be given by the user
 timespacedim : (= kc->dim = key->timespacedim)
              the true dimension for the location (if necessary by 
	      explicite parameter, e.g. in CovarianceFct.)
*/


#define DENS_LOG 0
#define DENS_SEED 1
#define DENS_ENV 2
/*
void density(double VARIABLE_IS_NOT_USED *value, model *cov, double *v) {
  KEY_type *KT = cov->base;
  assert(!P0INT(DENS_LOG));
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  char errorloc_save[nErr orLoc];
  int ni = 0,
    err = NOERROR;
  double *res;
  location_type *loc = PrevLoc(cov);
  long vdimtot; vdimtot = loc->totalpoints * VDIM0;
  assert(VDIM0 == VDIM1);

  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array

  STRCPY(errorloc_save, KT->error_loc);

  PutRN Gstate();
  ERR0("stop : ni nae Zei falsch");
  double simu_seed = GLOBAL_UTILS->basic.seed + (ni - 1);
  addVariable((char*) "seed", &simu_seed, 1, 1, PENV(DENS_ENV)->sexp);
  eval(PLANG(DENS_SEED)->sexp, PENV(DENS_ENV)->sexp);
  GetRN Gstate();

  SPRINTF(KT->error_loc, "%.50s %d", errorloc_save, ni);
 
  assert(cov->Sgen != NULL);
  
  NotProgrammedYet("density");

  //  if (COVNR == DENSITY) LOGDENSITY(sub, cov->Sgen); 
  //  else if (COVNR == PROBAB) PROBAB(sub, cov->Sgen);

  if (sizeof(double) == sizeof(double) && false) {
    MEMCOPY(res, cov->rf, sizeof(double) * vdimtot);
  } else {
    int i; for (i=0; i<vdimtot; i++) {
      res[i] = cov->rf[i];
    }
  }
  
  if (!sub->simu.ac tive) 
    GERR1("could not perform multiple simulations. Is '%.50s == FALSE'?",
	  general[GENERAL_STORING]);

 ErrorHandling: 
 
  PutRN Gstate();
  
  if (err > NOERROR) {
    XERR(err);
  }
}


int check_density(model *cov) {
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = PrevLoc(cov);
  int j, err, frame;
  isotropy_type iso;
  Types type;
  //  bool vdim_close_together = GLOBAL.general.vdim_close_together;

  ASSERT_LOC_GIVEN;

 
  cov->simu.expected_number_simu = GLOBAL.general.expected_number_simu;
  if (cov->simu.expected_number_simu > 1 && 
      GLOBAL.general.expected_number_simu <= 1)
    SERR("expected number of simulations inconsistent");

  GLOBAL.internal.stored_init = GLOBAL.general.storing || 
    GLOBAL.general.expected_number_simu > 1;  
  
  if (cov->key == NULL) {
    domain_type dom = KERNEL;
    if (isProcess(sub)) {
      frame = GaussMethodType;
      type = ProcessType;
      iso = UNREDUCED;
    } else {
      frame = V ariogramType;
      type = PosDefType;
      iso = S YMMETRIC;
    }
    if (cov->frame = = Any Type) frame = Any Type;

    err = ERRORTYPECONSISTENCY;

    for (j=0; j<=2; j++) {
      if ((Type Consistency(type, sub) && 
	   (err = C HECK(sub, loc->timespacedim, OWNXDIM(0), type, 
			dom, iso, cov->vdim, frame)) == NOERROR) 
	  || isProcess(sub)) break;

      if (j==0) type = VariogramType;
      else {
	type = TrendType;
	dom = XONLY;
	iso = CARTESIAN_COORD;
      }
    }

    if (err != NOERROR) RETURN_ERR(err);
  } else {
    BUG;
    frame = fram e_o f_p rocess(SUBNR);
    if (fram e == R OLE_FAILED) BUG;

    if ((err = C HECK(sub, loc->timespacedim, OWNXDIM(0), ProcessType,
		     XONLY, 
		     isCartesian(PREVISO(0)) ? CARTESIAN_COORD : PREVISO(0),
		     cov->vdim, frame)) != NOERROR) {
      RETURN_ERR(err);
    }
  }

  setbackward(cov, sub);  
  int subvdim = sub->vdim[0];
  VDIM0=subvdim; 
  VDIM1=sub->vdim[1]; 

  if (cov->q == NULL) {
    QALLOC(1);
    cov->q[0] = 1.0; 
  }

  RETURN_NOERROR;
}


int struct_density(model *cov, model VARIABLE_IS_NOT_USED  **newmodel){
  model *next = cov->sub[0],
    *sub = next;
  location_type *loc = PrevLoc(cov);
  int err,
    subframe = R OLE_FAILED,
    nr = NEXTNR;

  //APMI(next);

  if (isVariogram(next) || is Trend(next)) {
    if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;
    
    if ((err = C HECK(sub, loc->timespacedim, OWNXDIM(0), ProcessType,
		     XONLY, 
		     isCartesian(PREVISO(0)) ? CARTESIAN_COORD : PREVISO(0), 
		     cov->vdim, GaussMethodType)) != NOERROR) {
      RETURN_ERR(err);
    }
    subframe = GaussMethodType;    
  } else if (is BernoulliProcess(next)) subframe = FRAME_ BERNOULLI;
  else if (isGaussMethod(next)) subframe = GaussMethodType;   
  else if (isBrMethod(next)) subframe = BrMethodType;    
  else if (nr = = POISSONPROC) subframe = Poisson;
  else if (nr = = SCHLATHERPROC) subframe = SchlatherType;
  else if (nr = = EXTREMALTPROC) subframe = SchlatherType;
  else if (nr = = SMIT HPROC) subframe = SmithType;
  else {
    ILLEGAL_FRAME;
  }
  sub->frame = subframe; // ansonsten muesste hier C-HECK aufgerufen werden
  // hoffentlich gut gehende Abkuerzung, dass S-TRUCT aufgerufen wird,
  // und danach C-HECK (was auf jeden Fall gemacht werden muss)

  cov->simu.act ive = next->simu.act ive = false; 
  sub->simu.expected_number_simu = cov->simu.expected_number_simu;

  if (PL >= PL_DETAILS) PRINTF("Struct Density\n");

  if ((err = STRUCT(sub, NULL)) != NOERROR) { RETURN_ERR(err); }
  if (PL >= PL_DETAILS) PRINTF("Checking Density\n");


  assert(cov->Sgen == NULL);
  NEW_STORAGE(gen);

  if (!sub->initialised) {
    if (PL >= PL_DETAILS) PRINTF("Struct Density C\n");
   
    if (//cov->key ! = NULL && // bei sub= =next waere der falsche frame gesetzt
	// irgenwie komisch, cov->key abzufragen und check(sub aufzurufen
	// aufgrund von Beispiel in RTpoisson geloescht
	(err = C HECK(sub, loc->timespacedim, OWNXDIM(0), ProcessType,
		     PREVDOM(0), PREVISO(0), cov->vdim,
		     subframe)) != NOERROR) {
      //
      RETURN_ERR(err);
    }
    
    if (PL >= PL_DETAILS) {
      PRINTF("\n\nStruct Density (%s, #=%d), after 2nd check:",
	     NICK(sub), sub-> gatter nr); // OK
      PMI(sub); // OK
    }
 
    
    if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) {
      RETURN_ERR(err); 
    }
  }
 
  cov->initialised = true;
  ReturnOtherField(cov, sub);

  cov->simu.act ive = sub->simu.act ive = true; 
  RETURN_NOERROR;
}

void range_density(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  booleanRange(DENS_LOG);
}
*/


#define SIMU_CHECKONLY 0
#define SIMU_SEED 1
#define SIMU_ENV 2

void simulate(double *N, model *cov, double *v){
  assert(!P0INT(SIMU_CHECKONLY));
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  errorloc_type errorloc_save;
  char
    format[20],
    back[]="\b\b\b\b\b\b\b\b\b\b\b", 
    prozent[]="%",
    pch = GLOBAL.general.pch;
  int ni, digits, //n, 
    nn,
    err = NOERROR,
    each = 0;
  double *res,
    realeach=0.0;
  simu_storage *simu = NULL;  
  location_type *loc = PrevLoc(cov);
  finaldofct finalDo = DefList[SUBNR].FinalDo;
  long vdimtot = loc->totalpoints * VDIM0;
  assert(VDIM0 == VDIM1);
  KEY_type *KT = cov->base;

  // die folgenden Zeilen duefen nicht abgeaendert werden!!
  // Grund; *N wird eventuell durch Koordinatentrafo veraendert
  // so dass bei v=NULL ohne Umweg ueber gatter aufgerufen wird
  // und der korrekte Wert in cov->q gespeichert wird, der dann
  // spaeter wieder ausgelesen wird und nn zugeordnet wird.
  if (v == NULL) {
    cov->q[cov->qlen - 1] =  (int) *N;
    return; // EvaluateModel needs information about size
  //                      of result array
  } 
  nn = (int) cov->q[cov->qlen - 1];

  STRCPY(errorloc_save, KT->error_loc);
  simu = &(cov->simu);
  if (!simu->active) {
    err=ERRORNOTINITIALIZED; goto ErrorHandling;
  }

  if (nn>1 && pch != '\0') {
    if (pch == '!') {
      digits = (nn<900000000) 
	? 1 + (int) TRUNC(LOG((double) nn) / LOG(10.0)) : 9;
      back[digits] = '\0';
      each = (nn < 100) ? 1 :  nn / 100;
      SPRINTF(format, "%.2ss%.2s%dd", prozent, prozent, digits);
    } else if (pch == '%') {
      back[4] = '\0';
      realeach = (double) nn / 100.0;
      each = (nn < 100) ? 1 : (int) realeach;
      SPRINTF(format, "%.2ss%.2s%dd%.2ss", prozent, prozent, 3, prozent);
    } else each = 1;
  } else each = nn + 1;
  // size = cov->vdim * nn * loc->totalpoints;
  res = v;

  sub->simu.pair = false;
  for (ni=1; ni<=nn; ni++, res += vdimtot) {
#ifdef DO_PARALLEL
    //omp_set_num_threads(1);
#endif

    if (GLOBAL_UTILS->basic.seed != NA_INTEGER &&
	(nn > 1 || GLOBAL.general.seed_incr != 0)) {
      assert(!PisNULL(SIMU_SEED) && !PisNULL(SIMU_ENV));
      PutRNGstate(); // PutRNGstate muss vor UNPROTECT stehen!! -- hier irrelevant
      // if seed is increased by (ni -1) only, a row of simulations
      // of max-stable fields will look very much the same !!
      int simu_seed = GLOBAL_UTILS->basic.seed + (ni - 1) *
	GLOBAL.general.seed_sub_incr
	+ nn * GLOBAL.general.seed_incr;

      // printf("seed=%d %d %d\n",simu_seed, GLOBAL_UTILS->basic.seed,GLOBAL.general.seed_incr);
      
      addIntVariable((char*) "seed", &simu_seed, 1, 1, PENV(SIMU_ENV)->sexp);
      eval(PLANG(SIMU_SEED)->sexp, PENV(SIMU_ENV)->sexp);
      GetRNGstate(); // start interner modus
    }

    
    SPRINTF(KT->error_loc, "%.50s %d", errorloc_save, ni); 
    
    if (PL > 0) {
      if (ni % each == 0) {
	if (pch == '!')  
	  PRINTF(format, back, ni / each);
	else if (pch == '%')
	  PRINTF(format, back, (int) (ni / realeach), prozent);
      else PRINTF("%c", pch);
      }
    }
  
    assert(cov->Sgen != NULL);
   
    DO(sub, cov->Sgen);

    //    printf("ni=%d vdimtot=%d %ld sub=%ld\n", ni, vdimtot, res, sub->rf);
    //if(false) for (int k=0; k<vdimtot; k++)  if (sub->rf[k] <= 0.0) APMI(sub);

    
#ifdef DO_PARALLEL
    if (GLOBAL_UTILS->basic.cores>1)
      omp_set_num_threads(GLOBAL_UTILS->basic.cores);
    R_CheckUserInterrupt();
    if (GLOBAL_UTILS->basic.cores>1) omp_set_num_threads(1); // !! wichtig,
    // auch fuer nachfolgende DO()!!
#endif    
    
    MEMCOPY(res, cov->rf, sizeof(double) * vdimtot);
    // for (int i=0; i<vdimtot; i++) res[i] = cov->rf[i];
  
    if (!sub->simu.active) 
      GERR1("could not perform multiple simulations. Is '%.50s == FALSE'?",
	    general[GENERAL_STORING]);
  }
  
  if (finalDo != NULL) {
    finalDo(sub, v, nn, cov->Sgen);
  }

#ifdef DO_PARALLEL
  omp_set_num_threads(GLOBAL_UTILS->basic.cores);
#endif    

  if (nn>1 && pch != '\0') {
    if (pch == '!' || pch == '%') PRINTF("%s", back);
    else PRINTF("\n");
  }

  assert(simu != NULL);
  sub->simu.active = simu->active = sub->simu.active && GLOBAL.general.storing;

 ErrorHandling:

  //  PMI(cov);
   
  if (err > NOERROR) {
    if (simu != NULL) sub->simu.active = simu->active = false;
    XERR(err);
  }
}


int check_simulate(model *cov) {
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = PrevLoc(cov);
  int j, d, 
    err = NOERROR,
    dim = GetLoctsdim(cov);
  isotropy_type iso;
  Types type, frame;
  bool vdim_close_together = GLOBAL.general.vdim_close_together;
  
  cov->initialised = false;
  ASSERT_LOC_GIVEN;
  kdefault(cov,
	   SIMU_CHECKONLY, false);

  cov->simu.expected_number_simu = GLOBAL.general.expected_number_simu;
  if (cov->simu.expected_number_simu > 1 && 
      GLOBAL.general.expected_number_simu <= 1)
    SERR("expected number of simulations inconsistent");

  if (GLOBAL.general.seed_incr != 0 && GLOBAL_UTILS->basic.seed == NA_INTEGER) {
    SERR("'seed' must be set if 'seed_incr' is different from 0.");
  }

  GLOBAL.internal.stored_init = GLOBAL.general.storing || 
    GLOBAL.general.expected_number_simu > 1;  
  
  if (cov->key == NULL) {
    domain_type dom; 
    if (isProcess(sub)) {
      dom = XONLY; // 5/2015; voher beides Kerne
      frame = InterfaceType, //InterfaceType; // ProcessType; // Interface
      type = ProcessType;
      // iso = CoordinateSystemOf(PREVISO(0));
    } else {
      dom = KERNEL;
      frame = EvaluationType;
      type = PosDefType;
    }
    iso = PREVISO(0);
    assert(equalsCoordinateSystem(iso));
    assert(equalsCoordinateSystem(PREVISO(0)));

    if (hasAnyEvaluationFrame(cov)) BUG;
    // if (cov->frame = = Any Type) frame = Any Type;

    int errold = ERRORTYPECONSISTENCY;
    for (j=0; j<=2; j++) {
       //      PMI0(cov);     printf("j=%d ownxdim = %d\n", j, OWNXDIM(0));
      // BUG;
     if ( /// auskommenieren ? 8.8.17 ?? Type Consistency(type, sub, 0) &&
	 (err = CHECK(sub, dim, OWNXDIM(0), type, 
		       dom, iso, cov->vdim, frame)) == NOERROR) {
	break;
      }
      if (isProcess(sub)) {// muss als zweites stehen !!
	break;
      }

      if (j==0) {
	errold = err; // Fehler(meldungen) werden abstrus, falsch model falsch
	type = VariogramType;
	errold = err;
      } else {
	type = TrendType;
	dom = XONLY;
	iso = PREVISO(0);
      }
    }
    if (err != NOERROR) {
      err = j==0 ? err : errold;
      RETURN_ERR(err);
    }
  } else {
    frame = InterfaceType;
    assert(equalsCoordinateSystem(PREVISO(0)));

    if ((err = CHECK(sub, dim, OWNXDIM(0), ProcessType, XONLY, PREVISO(0),
		     cov->vdim, frame)) != NOERROR) {
      RETURN_ERR(err);
    }
  }
  
  setbackward(cov, sub);  
  int subvdim = sub->vdim[0];
  VDIM0=subvdim; 
  VDIM1=sub->vdim[1]; 

  if (cov->q == NULL) {
    int len=1;
    if (loc->grid) len += loc->timespacedim; else len++;
    if (subvdim > 1) len++;
    QALLOC(len);
    cov->q[--len] = 1.0; // number of simulations
    if (subvdim > 1 && !vdim_close_together) cov->q[--len]=(double) subvdim;
    if (loc->grid) {
      for (d=loc->timespacedim-1; d>=0; d--) 
	cov->q[--len] = loc->xgr[d][XLENGTH];
    } else {
      cov->q[--len] = loc->totalpoints;
    }
    if (subvdim > 1 && vdim_close_together) {
      assert(len == 1);
      cov->q[--len]=(double) subvdim;
    }
    assert(len==0);
  }

  cov->initialised = true;
  RETURN_NOERROR;
}


int struct_simulate(model *cov, model VARIABLE_IS_NOT_USED  **newmodel){  
  model *next = cov->sub[0],
    *sub = next;
  //location_type *loc = PrevLoc(cov);
  int err;
  Types subframe = InterfaceType;

  if (isnowVariogram(next) || isnowTrend(next)) {
    if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
    assert(cov->key->calling == next->calling); 
    addModel(&(cov->key), isnowVariogram(next) ? GAUSSPROC : TREND_PROC);
    sub = cov->key;

    if ((err = CHECK(sub, GetLoctsdim(cov), OWNXDIM(0), ProcessType, XONLY, 
		     isCartesian(PREVISO(0)) ? CARTESIAN_COORD :PREVISO(0), //??
		     cov->vdim, InterfaceType)) != NOERROR) {
       RETURN_ERR(err);
    }    
  }
  
  sub->frame = subframe; // ansonsten muesste hier C-HECK aufgerufen werden
  // hoffentlich gut gehende Abkuerzung, dass S-TRUCT aufgerufen wird,
  // und danach C-HECK (was auf jeden Fall gemacht werden muss)

  cov->simu.active = next->simu.active = false; 
  sub->simu.expected_number_simu = cov->simu.expected_number_simu;
  if (P0INT(SIMU_CHECKONLY)) RETURN_NOERROR;

  if (PL >= PL_DETAILS) { PRINTF("Struct Simulate\n"); }
 
  if ((err = STRUCT(sub, NULL)) != NOERROR) {
    // printf(">>> %s\n", cov->base->error_loc);
    RETURN_ERR(err);
  }
   
  if (PL >= PL_DETAILS) { PRINTF("Checking Simulate\n");}


  assert(cov->Sgen == NULL);
  NEW_STORAGE(gen);

    
   if (!sub->initialised) {
    if (PL >= PL_DETAILS) { PRINTF("Struct Simulate C\n");}
   
    if (//cov->key ! = NULL && // bei sub= =next waere der falsche frame gesetzt
	// irgenwie komisch, cov->key abzufragen und check(sub aufzurufen
	// aufgrund von Beispiel in RTpoisson geloescht
	(err = CHECK_PASSTF(sub, ProcessType, cov->vdim[0], subframe))
	!= NOERROR)
      //	(err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
      //	     cov->domprev, cov->isoprev, cov->vdim,
      //	     subframe)) != NOERROR)
      RETURN_ERR(err);   
    
    if (PL >= PL_DETAILS) {
      PRINTF("\n\nStruct Simulate (%s, #=%d), after 2nd check:",
	     NICK(sub), SUBNR); // OK
      PMI(sub); // OK
    }
 
    if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) {
      RETURN_ERR(err); 
    }
  }
 
  cov->initialised = true;
  ReturnOtherField(cov, sub);
  
  //  PMI(sub);
  assert( sub->simu.active );
  cov->simu.active = sub->simu.active; 

  RETURN_NOERROR;
}

void range_simulate(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  booleanRange(SIMU_CHECKONLY);

#define simu_n 2
  int idx[simu_n] = {SIMU_SEED, SIMU_ENV};
  for (int i=0; i<simu_n; i++) {
    int j = idx[i];
    range->min[j] = RF_NAN;
    range->max[j] = RF_NAN;
    range->pmin[j] = RF_NAN;
    range->pmax[j] = RF_NAN;
    range->openmin[j] = false;
    range->openmax[j] = false; 
  }
}

//void do_simulate(model *cov, gen_storage *s){
//  assert(false);
//}

/////////////////////////////////////////




double GetPriors(model *cov) {
  defn *C = DefList + COVNR;
  model *kap;
  int i,
    kappas = C->kappas,
    nsub = cov->nsub;
  double v,
    logli = 0.0;
  for (i=0; i<kappas; i++) {
    if ((kap = cov->kappasub[i]) != NULL) {
      if (isnowRandom(kap)) {
	if (C->kappatype[i] < LISTOF) {
	  assert(C->kappatype[i] == REALSXP);
	  VTLG_DLOG(P(i), kap, &v);
	} else if (C->kappatype[i] < LISTOF + LISTOF) { // not tested yet
	  NotProgrammedYet("hierachical models for multiple data sets");
	  assert(C->kappatype[i] == LISTOF + REALSXP);
	  int
	    store = GLOBAL.general.set,
	    end = cov->nrow[i];
	  for (GLOBAL.general.set = 0; GLOBAL.general.set<end; 
	       GLOBAL.general.set++) {
	    VTLG_DLOG(LP(i), kap, &v);
	    logli += v;
	  }
	  GLOBAL.general.set = store;
	} else BUG;
	logli += v;
      }
      logli += GetPriors(kap);
    }
  }
  
  for (i=0; i<nsub; i++) {
    logli += GetPriors(cov->sub[i]);
  }

  return logli; 
}





/* *************** */
/* LOG LIKELIHOOD */
/* *************** */

 void kappalikelihood(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc) {
  *nr = *nc = i == LIKELIHOOD_DATA ? 0 
    : i <= LIKELIHOOD_LAST ? 1 : -1;
}


void likelihood(double VARIABLE_IS_NOT_USED *x, model *cov, double *v) { 
  model *process = cov->key == NULL ? cov->sub[0] : cov->key;
  if (v == NULL) {
    likelihood_storage *L = process->Slikelihood;
    assert(L != NULL);
    int //betas = L->betas[L->fixedtrends],
      vdim = process->vdim[0];
    likelihood_info *info = &(L->info);
    listoftype *datasets = L->datasets;
    int
      store = GLOBAL.general.set,
      betatot = L->betas[L->fixedtrends],
      all_betatot = betatot;
    GLOBAL.general.set = 0;
    if (L->betas_separate)
      all_betatot *= NCOL_OUT_OF(datasets) / vdim;
    cov->q[0] = 1 + info->globalvariance + all_betatot;
    cov->q[1] = 1;
    GLOBAL.general.set = store;
    return; 
  }

  assert(isProcess(process));
  VTLG_DLOG(NULL, process, v); 

  assert(process->key == NULL);
  *v += GetPriors(process->sub[0]); 
 }

 

int check_likelihood(model *cov) {
  int err,
    sets = GET_LOC_SETS(cov);
  int store = GLOBAL.general.set; 

  if ((err = check_linearpart(cov))) RETURN_ERR(err);
  
  kdefault(cov, LIKELIHOOD_NA_VAR, GLOBAL.fit.estimate_variance);
  kdefault(cov, LIKELIHOOD_BETASSEPARATE, false);
  if (P0INT(LIKELIHOOD_BETASSEPARATE)) BUG; // to estimate beta 
  // for each independent repetition within a data set is a feature
  // that has been implemented, but currently it is not used
  kdefault(cov, LIKELIHOOD_IGNORETREND, false);
  if (PisNULL(LIKELIHOOD_DATA)) BUG;
  for (GLOBAL.general.set = 0; GLOBAL.general.set<sets; GLOBAL.general.set++) {
    long 
      datatot = LNROW(LIKELIHOOD_DATA) * LNCOL(LIKELIHOOD_DATA),
      totpts = Gettotalpoints(cov),
      vdimtot = totpts * VDIM0;
    int
      repet = datatot / vdimtot;

    //  printf("repet = %d %d %d %d\n", repet,  vdimtot , VDIM0,datatot);
    
    if (repet * vdimtot != datatot || repet == 0)  {
      GERR("data and coordinates do not match");
    }
    LNRC_(LIKELIHOOD_DATA, nrow) = totpts;
    LNRC_(LIKELIHOOD_DATA, ncol) = datatot / totpts;
  }

 ErrorHandling:
  GLOBAL.general.set = store;

  RETURN_ERR(err);
}

int struct_likelihood(model *cov, 
		      model VARIABLE_IS_NOT_USED  **newmodel){
  //return struct_rftrend(cov, newmodel);
  assert(cov->key == NULL);
  model *sub = cov->sub[0];
  int err = NOERROR;
  location_type *loc = Loc(cov);

  if (isnowVariogram(sub)) {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) RETURN_ERR(err);
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;
    if ((err = CHECK(sub, loc->timespacedim, OWNXDIM(0), ProcessType, XONLY, 
		     isCartesian(PREVISO(0)) ? CARTESIAN_COORD : PREVISO(0),
		     cov->vdim, LikelihoodType)) != NOERROR) {
      RETURN_ERR(err);
    }
  } else sub->frame = LikelihoodType;

  if (!isnowProcess(sub)) 
    SERR1("'%.50s' can be calculated only for processes.", NICK(cov));
  
  if ((err = STRUCT(sub, NULL)) != NOERROR) RETURN_ERR(err);
  NEW_STORAGE(gen);
  if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) RETURN_ERR(err);
  
  RETURN_NOERROR;
}


void range_likelihood(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[LIKELIHOOD_DATA] = RF_NEGINF;
  range->max[LIKELIHOOD_DATA] = RF_INF;
  range->pmin[LIKELIHOOD_DATA] = - 1e8;
  range->pmax[LIKELIHOOD_DATA] = 1e8;
  range->openmin[LIKELIHOOD_DATA] = true;
  range->openmax[LIKELIHOOD_DATA] = true;
  
  booleanRange(LIKELIHOOD_NA_VAR);
  booleanRange(LIKELIHOOD_BETASSEPARATE);
  booleanRange(LIKELIHOOD_IGNORETREND);
}



void linearpart(double VARIABLE_IS_NOT_USED *x, model  VARIABLE_IS_NOT_USED *cov, double  VARIABLE_IS_NOT_USED *v) { 
  BUG; // darf nicht aufgerufen werden
}

SEXP get_linearpart(SEXP model_reg, SEXP Set){
   int cR = INTEGER(model_reg)[0];
   set_currentRegister(cR);
  if (cR < 0 || cR > MODEL_MAX) BUG;	
  model *cov = KEY()[cR];			
  cov = cov->key != NULL ? cov->key : cov->sub[0];
  //  TREE(cov);
  if (COVNR == GAUSSPROC) 
    return gauss_linearpart(model_reg, Set);
  else BUG;
  //  TREE(cov); 
}


int check_linearpart(model *cov) {
  assert(cov->prevloc != NULL);
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  int err, 
    dim = GetLoctsdim(cov);
  domain_type dom = KERNEL;
  Types frame = LikelihoodType;
    // SUBNR == PLUS || isProcess(sub) ? LikelihoodType : EvaluationType;
  
  //if (isEvaluationType(cov)) frame = EvaluationType;
  
  err = ERRORTYPECONSISTENCY;
  if (DistancesGiven(cov)) 
    SERR1("'%.50s' may not be used when distances are given", NAME(cov));

  if (isProcess(sub)) {
    err = CHECK(sub, dim, dim, ProcessType, dom, UNREDUCED, cov->vdim, frame);
  } else {
    err = CHECK(sub, dim, dim, PosDefType, dom,
		CoordinateSystemOf(PREVISO(0)),
		cov->vdim,
		frame);
    if (err != NOERROR)
      err = CHECK(sub, dim, dim, NegDefType, dom, SymmetricOf(PREVISO(0)),
		  cov->vdim,
		  frame);
  }
  if (err != NOERROR) RETURN_ERR(err);

  setbackward(cov, sub);  
  VDIM0=sub->vdim[0]; 
  VDIM1=sub->vdim[1]; 

  if (cov->q == NULL) QALLOC(2);
  cov->q[0] = Gettotalpoints(cov);
  cov->q[1] = VDIM0;

  RETURN_ERR(err);
}



int struct_linearpart(model *cov, 
    model VARIABLE_IS_NOT_USED  **newmodel){
  assert(cov->key == NULL);
  model *sub = cov->sub[0];
  int err = NOERROR;
  location_type *loc = Loc(cov);

  if (isnowVariogram(sub)) {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) RETURN_ERR(err);
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;
    
    if ((err = CHECK(sub, loc->timespacedim, OWNXDIM(0), ProcessType, XONLY, 
		     isCartesian(PREVISO(0)) ? CARTESIAN_COORD : PREVISO(0),
		     cov->vdim, LikelihoodType)) != NOERROR) {
      RETURN_ERR(err);
    }
  } else sub->frame = LikelihoodType;
  
  if (!isnowProcess(sub)) 
    SERR1("'%.50s' can be calculated only for processes.", NICK(cov));


  if ((err = STRUCT(sub, NULL)) != NOERROR) RETURN_ERR(err);
  likelihood_storage *L = sub->Slikelihood;
  if (L == NULL) RETURN_ERR(ERRORFAILED);
  if (L->dettrend_has_nas || L->fixedtrend_has_nas) {
    warning("NAs detected in '%20s'; hence zero's introduced", NICK(cov)); 
  }
  RETURN_ERR(err);
}




#define EVALDISTR_X 0 // d
#define EVALDISTR_Q 1 // p
#define EVALDISTR_P 2 // q
#define EVALDISTR_N 3 // r
#define EVALDISTR_DIM 4 // r

void kappa_EvalDistr(int i, model VARIABLE_IS_NOT_USED *cov, 
		     int *nr, int *nc){
  *nc = *nr = i <= EVALDISTR_N ? 0 : i == EVALDISTR_DIM ? 1 : -1;
}

void EvalDistr(double VARIABLE_IS_NOT_USED *N, model *cov, double *v){
  errorloc_type errorloc_save;
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  double  *xqp;
  int i, j,
    dim = ANYDIM,
    n = (int) (cov->q[cov->qlen - 1]);
  KEY_type *KT = cov->base;

  if (v==NULL) return; // EvaluateModel needs information about size
  STRCPY(errorloc_save, KT->error_loc);

  if (!PisNULL(EVALDISTR_X)) { // d
    xqp = P(EVALDISTR_X);
    for (j=i=0; i<n; i++, j+=dim) VTLG_D(xqp + j, sub, v + i);
  } else if (!PisNULL(EVALDISTR_Q)) { // p
    xqp = P(EVALDISTR_Q);
    for (j=i=0; i<n; i++, j+=dim) VTLG_P(xqp + i, sub, v + j);
  } else if (!PisNULL(EVALDISTR_P)) {// q
    xqp = P(EVALDISTR_P);
    for (j=i=0; i<n; i++, j+=dim) VTLG_Q(xqp + j, sub, v + i);
  } else if (!PisNULL(EVALDISTR_N)) {// r
    xqp = P(EVALDISTR_N);
    for (j=i=0; i<n; i++, j+=dim) {
      VTLG_R(NULL, sub, v + j);    
    }
  } else BUG;
}


int check_EvalDistr(model *cov) {
  defn *C = DefList + COVNR;
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  int size_idx,  err, //type,  nr = SUBNR,
    *nrow = cov->nrow,
    *ncol = cov->ncol,
    dim = ANYDIM, 
    zaehler=0;

  
  if (cov->q == NULL) {
    int nn = 1; // !! 1 obwohl 2 reserviert werden !! Wg EvaluateModel
    if (dim > 1 && 
	((!PisNULL(EVALDISTR_N) && P0(EVALDISTR_N) > 1) ||
	 (!PisNULL(EVALDISTR_Q) && P0(EVALDISTR_Q) > 1)))
      nn++;
    QALLOC(nn + 1); 
    cov->qlen = nn;
   /*  cov->qlen = 1; // !! 1 obwohl 2 reserviert werden !! Wg EvaluateModel
    if (dim > 1 && 
	((!PisNULL(EVALDISTR_N) && P0(EVALDISTR_N) > 1) ||
	 (!PisNULL(EVALDISTR_Q) && P0(EVALDISTR_Q) > 1)))
      cov->qlen++;
      cov->q = (double *) MALLOC(sizeof(double) * (cov->qlen + 1)); // QALLOC NOT APPROPRIATE
   */

    cov->q[0] = dim; // is overwritten if not a matrix is returned
    size_idx = cov->qlen - 1;

    if (PisNULL(EVALDISTR_N)) {
      if (!PisNULL(EVALDISTR_X)) { // d
	if (dim > 1 && nrow[EVALDISTR_X] != dim) 
	  SERR2("dimenson of '%.50s' does not match '%.50s' ", 
		C->kappanames[EVALDISTR_X], C->kappanames[EVALDISTR_DIM]);
	cov->q[size_idx] = nrow[EVALDISTR_X] * ncol[EVALDISTR_X] / dim;
	zaehler++;
      }
      if (!PisNULL(EVALDISTR_Q)) { // p
	if (dim > 1 && nrow[EVALDISTR_Q] != dim) 
	  SERR2("dimension of '%.50s' does not match '%.50s' ", 
		C->kappanames[EVALDISTR_Q], C->kappanames[EVALDISTR_DIM]);
	cov->q[size_idx] = nrow[EVALDISTR_Q] * ncol[EVALDISTR_Q] / dim;
	zaehler++;
      } 
      if (!PisNULL(EVALDISTR_P)) { // q
	if (ncol[EVALDISTR_P] != 1) 
	  SERR1("'%.50s' must be a vector", C->kappanames[EVALDISTR_P]);
	cov->q[size_idx] = nrow[EVALDISTR_P]  * dim;
 	zaehler++;
      } 
    }// kein else wegen zaehler !!

    if (!PisNULL(EVALDISTR_N)) { // r      
      cov->q[size_idx] = P0(EVALDISTR_N) * dim;
      zaehler++;
    } 
    if (zaehler != 1) SERR("exactly one of the parameters must be given");

  }
 
  //  if (!isRandom(sub)) SERR1("'%.50s' is not a distribution", NICK(sub));

  // PMI(sub);
  
  if ((err = CHECK_R(sub, dim)) != NOERROR) RETURN_ERR(err);
  assert(isnowRandom(sub));
 

  setbackward(cov, sub);  

  RETURN_NOERROR;
}

// S TRUCT(!isCovariance(cov) ? NULL : &(cov->key));



void range_EvalDistr(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[EVALDISTR_X] = RF_NEGINF;
  range->max[EVALDISTR_X] = RF_INF;
  range->pmin[EVALDISTR_X] = - 1e8;
  range->pmax[EVALDISTR_X] = 1e8;
  range->openmin[EVALDISTR_X] = true;
  range->openmax[EVALDISTR_X] = true;

  range->min[EVALDISTR_Q] = RF_NEGINF;
  range->max[EVALDISTR_Q] = RF_INF;
  range->pmin[EVALDISTR_Q] = - 1e8;
  range->pmax[EVALDISTR_Q] = 1e8;
  range->openmin[EVALDISTR_Q] = true;
  range->openmax[EVALDISTR_Q] = true;

  range->min[EVALDISTR_P] = 0;
  range->max[EVALDISTR_P] = 1;
  range->pmin[EVALDISTR_P] = 0;
  range->pmax[EVALDISTR_P] = 1;
  range->openmin[EVALDISTR_P] = false;
  range->openmax[EVALDISTR_P] = false;

  range->min[EVALDISTR_N] = 1;
  range->max[EVALDISTR_N] = RF_INF;
  range->pmin[EVALDISTR_N] = 1;
  range->pmax[EVALDISTR_N] = 1e8;
  range->openmin[EVALDISTR_N] = false;
  range->openmax[EVALDISTR_N] = false;

  range->min[EVALDISTR_DIM] = 1;
  range->max[EVALDISTR_DIM] = RF_INF;
  range->pmin[EVALDISTR_DIM] = 1;
  range->pmax[EVALDISTR_DIM] = 10;
  range->openmin[EVALDISTR_DIM] = false;
  range->openmax[EVALDISTR_DIM] = true;
}


int struct_EvalDistr(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  model *next = cov->sub[0],
    *sub = next;
  int  err,
    dim = ANYDIM;
   
  //  cov->simu.active = next->simu.active = false; 
  if (PL >= PL_DETAILS) { PRINTF("Struct EvalDistr\n"); }

  if ((err = STRUCT(sub, NULL)) != NOERROR) { RETURN_ERR(err); }
  if (PL >= PL_DETAILS) { PRINTF("Checking EvalDistr\n"); }
  
  if ((err = CHECK_R(sub, dim)) != NOERROR) RETURN_ERR(err);
    
  if (PL >= PL_DETAILS) {
    PRINTF("\n\nStruct EvalDistr (%s, #=%d), after 2nd check:",
	   NICK(sub), SUBNR); 
  }
  
  assert(cov->Sgen == NULL);
  NEW_STORAGE(gen);

  if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) {
    RETURN_ERR(err); 
  }

  if (cov->rf == NULL) {
    int size = cov->q[0];
    if (cov->qlen > 1) size *= cov->q[1];

    if ((cov->rf = (double*) MALLOC(sizeof(double) * size)) 
	== NULL) RETURN_ERR(ERRORMEMORYALLOCATION); 
    cov->fieldreturn = wahr;
    cov->origrf = true;
  }
  
  //  cov->simu.active = sub->simu.active = true; 
  RETURN_NOERROR;
}


//void do_EvalDistr(model *cov, gen_storage *s){
//  assert(false);
//}


void Dummy(double VARIABLE_IS_NOT_USED *x, model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *value) {
  BUG; //ERR("unexpected call of dummy function");
}


int check_dummy(model *cov) {
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = PrevLoc(cov);
  int err = NOERROR;
#define nTypes 2
  Types types[nTypes] = {NegDefType, ProcessType};
#define nFrames 2
  Types frames[nTypes] = {EvaluationType, GaussMethodType};
 
  ASSERT_LOC_GIVEN;

  for (int f=0; f<nFrames; f++) {
    for (int t=0; t<nTypes; t++) {
      for (int k=FIRST_DOMAIN; k<=LAST_DOMAINUSER; k++) {
	if ((err = CHECK(sub, loc->timespacedim, OWNXDIM(0), types[t],
			 (domain_type) k, CoordinateSystemOf(PREVISO(0)),
			 SUBMODEL_DEP, frames[f]))
	    == NOERROR) break;
      }
      if (err == NOERROR) break;
    }
    if (err == NOERROR) break;
  }
  if (err != NOERROR) RETURN_ERR(err);
  setbackward(cov, sub);  
  VDIM0 = sub->vdim[0]; 
  VDIM1 = sub->vdim[1]; 

  RETURN_NOERROR;
}


int struct_dummy(model VARIABLE_IS_NOT_USED *cov, model VARIABLE_IS_NOT_USED **newmodel){    
  RETURN_NOERROR;
}


int alloc_pgs(model *cov) {
  return alloc_pgs(cov, ANYDIM);
}

int alloc_pgs(model *cov, int dim) { // all what is necessary for dompp 
  int dimP1 = dim + 1;

  assert(cov->Spgs == NULL); // NIE pgs_DELETE(&(cov->Spgs)) in Huetchen, da sonst dompp durcheinander kommt;
  NEW_STORAGE_WITH_SAVE(pgs);
  pgs_storage *pgs = cov->Spgs;
  assert(pgs->z == NULL);

  // printf("alloc dim=%d %d\n", dim , dimP1);

  // to be save: everywhere +1 for loc->i_row although most not needed
  if ((pgs->supportmin = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
      (pgs->supportmax = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
      (pgs->supportcentre = (double*) CALLOC(dimP1, sizeof(double))) == NULL || 
      (pgs->own_grid_start = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
      (pgs->own_grid_step = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
      (pgs->own_grid_len = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
      
      (pgs->gridlen = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
      (pgs->end = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
      (pgs->start = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
      (pgs->delta = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
      (pgs->nx = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
   
      (pgs->xstart = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
       (pgs->x = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
      (pgs->xgr = (double**) CALLOC(dimP1, sizeof(double*))) == NULL ||
      (pgs->inc = (double*) CALLOC(dimP1, sizeof(double))) == NULL)
    RETURN_ERR(ERRORMEMORYALLOCATION);

  RETURN_NOERROR;
}

//93948675616576 z=
int alloc_cov(model *cov, int dim, int rows, int cols) { 
  // all what is necessary for dompp
  // but also for cov evaluation and fctn evaluation!

  int err;

  if (cov->Spgs != NULL) pgs_DELETE(&(cov->Spgs), cov);
  if ((err = alloc_pgs(cov, dim)) != NOERROR) RETURN_ERR(err);

  // erst danach!!!
  pgs_storage *pgs = cov->Spgs;
  int max,
    rowscols = rows * cols;
  max = rows;
  if (cols > max) max = cols;

  assert(pgs != NULL && pgs->x!=NULL && pgs->endy==NULL && pgs->ptrrow==NULL);

  if ( 
      (pgs->endy = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->startny = (int*) CALLOC(dim, sizeof(int))) == NULL ||

      (pgs->ptrcol = (int*) CALLOC(max, sizeof(int))) == NULL ||
      (pgs->ptrrow = (int*) CALLOC(max, sizeof(int))) == NULL)
    RETURN_ERR(ERRORMEMORYALLOCATION);

  if (
      (pgs->C0x  = (double*) CALLOC(rowscols, sizeof(double))) == NULL || 
      (pgs->C0y  = (double*) CALLOC(rowscols, sizeof(double))) == NULL || 
      (pgs->cross= (double*) CALLOC(rowscols, sizeof(double))) == NULL || 
      (pgs->z    = (double*) CALLOC(rowscols, sizeof(double))) == NULL ||       
      (pgs->Val= (double**) CALLOC(rowscols, sizeof(double*))) == NULL 
      ) RETURN_ERR(ERRORMEMORYALLOCATION);

  pgs->rowscols = rowscols;

  RETURN_NOERROR;
}


#define RFGET_UP 0
#define RFGET_REGISTER 1
#define RFGET_SUB 0
void RFget(double VARIABLE_IS_NOT_USED *x, model *cov, double *v){
  get_storage *s = cov->Sget;
  model *get_cov = s->get_cov;
  int i,
    nr = MODELNR(s->get_cov),
    param_nr = s->param_nr,
    *idx = s->idx,
    size = s->size;

  if (DefList[nr].kappatype[param_nr] == REALSXP) {
    double *p = PARAM(get_cov, param_nr);    
    if (s->all) {
      for (i=0; i<size; i++) v[i] = p[i];
    } else {
      for (i=0; i<size; i++) v[i] = p[(int) idx[i]];
    }
  } else if (DefList[nr].kappatype[param_nr] == INTSXP) {    
    int *p = PARAMINT(get_cov, param_nr);
    if (s->all) {
      for (i=0; i<size; i++) v[i] = (double) p[i];
    } else {
      for (i=0; i<size; i++) v[i] = (double) p[idx[i]];
    }
  } else BUG;
}



int SearchParam(model *cov, get_storage *s) {
  model *orig;
  int i, subcov,
    up = P0INT(RFGET_UP);
  if (PisNULL(RFGET_REGISTER)) {
    orig = cov->calling;
    if (orig == NULL) SERR("register not given");
    for (i=0; i<up; i++) {
      orig = orig->calling;
      while (orig != NULL && orig->user_given == ug_internal)
	orig = orig->calling;
      if (orig == NULL) SERR1("value of '%.50s' too high", KNAME(RFGET_UP));
    }
  } else {
    int nr = P0INT(RFGET_REGISTER);
    if (nr<0 || nr>MODEL_MAX) SERR("invalid register number");
    if (up != 0) SERR1("'%.50s' may not be given.", KNAME(RFGET_UP));
    orig = KEY()[nr];
  }
  s->orig = orig;
  
  while (true) {
    while (DefList[MODELNR(orig)].maxsub > 0 && orig != NULL &&
	     orig->user_given == ug_internal) 
      orig = (MODELNR(orig) == PLUS || MODELNR(orig) == MULT || MODELNR(orig)==MPPPLUS) 
	&& orig->Splus != NULL && orig->Splus->keys_given ? orig->Splus->keys[0]
	: orig->key != NULL ? orig->key
	: orig->sub[0];
    if (orig == NULL || DefList[MODELNR(orig)].maxsub == 0) 
      SERR("model does not match the request or is too complicated");
    if (COVNR != MODELNR(orig)) 
      SERR2("(sub)models (%.50s, %.50s) do not match",
	    NICK(cov), NICK(orig));
    for (subcov=0; subcov < orig->nsub; subcov++) 
      if (cov->sub[subcov] != NULL) break;   
    if (subcov < orig->nsub) { // submodel found
      if (orig->sub[subcov] == NULL) 
	SERR2("(sub)models (%.50s; %.50s) do not match",
	      NICK(cov), NICK(orig));
      cov  = cov->sub[subcov];
      orig = orig->sub[subcov];
    } else {
      int
	kappas = DefList[MODELNR(orig)].kappas;
      for (i=0; i < kappas; i++) 
	if (cov->kappasub[i] != NULL) break;         
      if (i < kappas) { // param submodel found
	if (orig->kappasub[i] == NULL) 
	  SERR2("parameter models of %.50s and %.50s do not match",
		NICK(cov), NICK(orig));
        cov  = cov->kappasub[i];
	orig = orig->kappasub[i];   
      } else {
	for (i=0; i < kappas; i++) 
	  if (cov->kappasub[i] != NULL) break;         
	if (i >= kappas) SERR("no parameter given");
	defn *C = DefList + COVNR;
	if (C->kappatype[i] == REALSXP) s->all = P(i)[0] == 0;
	else if (C->kappatype[i] == INTSXP) s->all = PINT(i)[0] == 0;
	else SERR("only numeric parameters are allowed");
	if (s->all) {
	  s->vdim[0] = orig->nrow[i];
	  s->vdim[1] = orig->ncol[i];
	  s->size = s->vdim[0] * s->vdim[1];
	} else {
	  int k;
	  s->size = s->vdim[0] = cov->nrow[i];
	  s->vdim[1] = cov->ncol[i];
	  if (cov->ncol[i] != 1) SERR("only vectors of indices are allowed");
	  assert(s->idx == NULL);
	  s->idx = (int *) MALLOC(sizeof(int) * s->size);
	  if (C->kappatype[i] == REALSXP)
	    for (k=0; k<s->size; k++) s->idx[k] = ((int) P(i)[k]) - 1;
	  else 
	    for (k=0; k<s->size; k++) s->idx[k] = PINT(i)[k] - 1;
	}
	s->get_cov = orig;
	s->param_nr = i;
	RETURN_NOERROR;
      }
    }
  } // while true
  BUG;
  RETURN_ERR(ERRORFAILED); // nur fuer compiler
}
  
int check_RFget(model *cov) {

  BUG; /// todo:  Code ist noch nicht ausgegoren !!

  //model *orig, *sub;
  int i, err;
    //    len = ((idx != NULL) ? cov->nrow[RFGET_IDX]
    //	  : get->get_cov->ncol[param_nr] * get->get_cov->nrow[param_nr]);
  if (cov->Sget != NULL) RETURN_NOERROR;

  kdefault(cov, RFGET_UP, 0);
  NEW_STORAGE(get);
  get_storage *s = cov->Sget;

  if ((err = SearchParam(cov, s)) != NOERROR) RETURN_ERR(err);
  for (i=0; i<2; i++) cov->vdim[i] = s->vdim[i];
  RETURN_NOERROR;
}
  

void range_RFget(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[RFGET_UP] = 0;
  range->max[RFGET_UP] = RF_INF;
  range->pmin[RFGET_UP] = 0;
  range->pmax[RFGET_UP] = 100;
  range->openmin[RFGET_UP] = false;
  range->openmax[RFGET_UP] = true;

  range->min[RFGET_REGISTER] = 0;
  range->max[RFGET_REGISTER] = MODEL_MAX;
  range->pmin[RFGET_REGISTER] = 0;
  range->pmax[RFGET_REGISTER] = MODEL_USER;
  range->openmin[RFGET_REGISTER] = false;
  range->openmax[RFGET_REGISTER] = false;
}


int struct_RFget(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  int i, err;
  
  NEW_STORAGE(get);
  get_storage *s = cov->Sget;
  if ((err = SearchParam(cov, s)) != NOERROR) RETURN_ERR(err);
  for (i=0; i<2; i++) {
    if (cov->vdim[i] != s->vdim[i])
      SERR("mismatch of dimensions when constructing the model");
  }
  cov->fieldreturn = wahr; // seltsam !
  cov->origrf = false;  
  RETURN_NOERROR;
}




//void do_Rfget(model *cov, gen_storage *s){
//  assert(false);
//}






//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

model *get_around_gauss(model *cov) {
  model *next= cov;
  if (NEXTNR == SCHLATHERPROC) next = next->sub[0]; // kein else 
  if (NEXTNR == GAUSSPROC) next = next->sub[0];

  if (isGaussMethod(next) || equalsBernoulliProcess(next)) {
    if (NEXTNR == AVERAGE_USER){
      if (next->sub[0] == NULL) ERR0("covariance cannot be calculated (yet) for arbitrary shape functions.");
      next = next->sub[next->sub[0] == NULL];
      if (NEXTNR == AVERAGE_INTERN) next = next->sub[next->sub[0] == NULL];
    } 
    else if (NEXTNR == CE_CUTOFFPROC_USER) {
      next = next->sub[0];      
      if (NEXTNR == CE_CUTOFFPROC_INTERN) next = next->sub[0];
    }
    else if (NEXTNR == CE_INTRINPROC_USER) {
      next = next->sub[0];   
      if (NEXTNR == CE_INTRINPROC_INTERN) next = next->sub[0];   
    }
    else if (NEXTNR == CE_INTRINPROC_USER) {
      next = next->sub[0];   
      if (NEXTNR == CE_INTRINPROC_INTERN) next = next->sub[0];   
    }
    else if (NEXTNR == HYPERPLANE_USER) {
      next = next->sub[0];   
      if (NEXTNR == HYPERPLANE_INTERN) next = next->sub[0];   
    }
    else if (NEXTNR == RANDOMCOIN_USER) {
      if (next->sub[0] == NULL) ERR0("covariance cannot be calculated (yet) for arbitrary shape functions.");
      next = next->sub[next->sub[0] == NULL];   
      if (NEXTNR == AVERAGE_INTERN) next = next->sub[next->sub[0] == NULL];
    } else {
      BUG;
    }
  }
  return next;
}



model *get_around_max_stable(model *cov) {
  model *next = cov;

  if (isBrMethod(next)) {
    next = next->sub[0];
    if (MODELNR(next->calling) == BROWNRESNICKPROC && isBrMethod(next)) {
      next = next->sub[0];
    }
  } 
  return next;
}



int check_fct_intern(model *cov, Types type, bool close, bool kernel,
		     int rows, int cols, Types frame) {
  model
    *next = cov->sub[0],
    *sub = cov->key == NULL ? next : cov->key;
  location_type *loc = Loc(cov);
  ASSERT_LOC_GIVEN;

  int err = ERRORFAILED,   // +  (int) nr_coord_sys_proj], 
     dim =  GetLoctsdim(cov);
 assert(dim == OWNLOGDIM(0));

 domain_type
   firstdomain = isNegDef(type) && isAnySpherical(OWNISO(0)) ? KERNEL : XONLY,
   lastdomain = kernel && !(isTrend(type) || isProcess(type)) ? KERNEL : XONLY;

 // assert(firstdomain == XONLY);

   /*
  int endfor = 0,
  isotropy_type iso[4];
  iso[endfor++] = type == ShapeType 
    ? CoordinateSystemOf(OWNISO(0)) 
    : S ymmetricOf(OWNISO(0));
  
   if (isIsoMismatch(iso[endfor-1])) BUG;
  int end_frame = sub = = next;
  for (int i=0; i < endfor; i++) {
    for (k=XONLY; k<=lastdomain; k++) {
      Types frame = V ariogramType;
      for (int f=0; f<=end_frame; f++) {
	if ((err = CHECK(sub, dim, OWNXDIM(0), type, 
			 (domain_type) k, iso[i], SUBMODEL_DEP, frame))
	    //sub!=next || isVariogram(sub) ? V ariogramType : Any Type))
	    == NOERROR) break;
	if (err == NOERROR) break;
	frame = Any Type;
      }
      if (err == NOERROR) break;
    }
    if (err == NOERROR) break;
  }
  */
 isotropy_type iso = equalsVariogram(type) || equalsNegDef(type)
    ? SymmetricOf(OWNISO(0))
    : CoordinateSystemOf(OWNISO(0));  
  if (equalsIsoMismatch(iso)) BUG;
  
  //  Types frame = isProcess(type) ? EvaluationType : EvaluationType;
  //  int end_frame = sub == next && equalsEvaluation(frame);

  //  assert(frame = EvaluationType);

  
  for (int k=firstdomain; k<=lastdomain; k++) {
    // for (int f=0; f<=end_frame; f++) {
      if ((err = CHECK(sub, dim, OWNXDIM(0),
		       type, 
		       (domain_type) k, iso, SUBMODEL_DEP, frame))
	  //sub!=next || isVariogram(sub) ? V ariogramType : Any Type))
	  == NOERROR) break;

      //PMI(cov);
      //
      //BUG;
      //    TREE0(cov);
      //  if (equalsEvaluation(frame)) BUG;      
      //frame = EvaluationType;
      ///}
      // if (err == NOERROR) break;
  }
  
  if (err != NOERROR) RETURN_ERR(err);
  setbackward(cov, sub); 

  // memory set according to the submodel as this model will
  // be evaluated (cov clear; fctn as function; predict: cov-matrix!)
  if ((err = alloc_cov(cov, dim, VDIM0, VDIM1)) != NOERROR)
    RETURN_ERR(err);

  // this is how cov will forward the result
  if (rows > 0) VDIM0 = rows;
  if (cols > 0) VDIM1 = cols;

  if (sub->pref[Nothing] == PREF_NONE) SERR("given model cannot be evaluated")
  
  if (cov->q == NULL) {
    int d,
      len=1; // # of "simulations" (repetitions!)
    if (loc->grid) len += dim; else len ++;      
    for (int i=0; i<2; i++) len += (int) (cov->vdim[i] > 1);
    QALLOC(len);

#define VDIMS  for (int i=0; i<2; i++) if (cov->vdim[i] > 1) cov->q[d++] = cov->vdim[i]
#define LOCS if (loc->grid) {						\
      for (int i=0; i<dim; i++) cov->q[d++] = loc->xgr[i][XLENGTH];	\
    } else {								\
      cov->q[d++] = loc->totalpoints;					\
    }      
    
    d = 0;
    if (close) {
      VDIMS;
      LOCS;	
   } else {
      LOCS;
      VDIMS;
    }
    cov->q[d] = 1;
    assert(d == len-1);
  }

  RETURN_NOERROR;
}


void Cov(double VARIABLE_IS_NOT_USED *x, model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  DefList[SUBNR].covariance(sub, value);
}

isotropy_type which_system[nr_coord_sys] =
  { CARTESIAN_COORD, ISO_MISMATCH, CARTESIAN_COORD, EARTH_COORD,
    SPHERICAL_COORD, CARTESIAN_COORD, CARTESIAN_COORD, ISO_MISMATCH };

int check_cov_intern(model *cov, Types type, bool close, bool kernel) {
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  if (isProcess(sub)) {
    int err,
      dim =  GetLoctsdim(cov);
    
     err = CHECK_THROUGHOUT(sub, cov, ProcessType, XONLY,
			    which_system[GLOBAL.coords.coord_system],
			    SUBMODEL_DEP,
			    EvaluationType);
     //  int err = CHECK(sub, cov->tsdim, cov->xdimown, ProcessType, XONLY, 
     //		    OWNISO(0), SUBMODEL_DEP, 
     //		    cov->frame = = Any Type ? Any Type : GaussMethodType);
    if (err != NOERROR) RETURN_ERR(err);
    setbackward(cov, sub);
    VDIM0 = sub->vdim[0];
    VDIM1 = sub->vdim[1];
    if ((err = alloc_cov(cov, dim, VDIM0, VDIM1)) != NOERROR) RETURN_ERR(err);
    RETURN_NOERROR;
  } else return check_fct_intern(cov, type, close, kernel, 0,0, EvaluationType);
}

int check_cov(model *cov) {
  return check_cov_intern(cov, PosDefType, GLOBAL.general.vdim_close_together,
			  true);
}

int struct_cov(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  int err;
  model
    *next = cov->sub[0],
    *sub = get_around_gauss(next);

  // fehlt : Poisson

  if (sub != next) {
    if ((err = COVNR == COVMATRIX ? check_covmatrix(cov) : check_cov(cov))
	!= NOERROR) RETURN_ERR(err);


    //    if (cov->Spgs == NULL) { PMI0(cov); crash(); }
    
    assert(cov->Spgs != NULL);
    
    ONCE_NEW_STORAGE(gen);

   if ((err = INIT(next, 0, cov->Sgen)) != NOERROR) RETURN_ERR(err);
    RETURN_ERR(err);
  }
  RETURN_NOERROR;
}


int init_cov(model *cov, gen_storage *s) {
  // darf nur von Likelihood aus aufgerufen werden
  if (hasAnyEvaluationFrame(cov)) {
    BUG;
    assert(cov->key == NULL);
    return INIT(cov->sub[0], 0, s);
  }

  RETURN_ERR(ERRORFAILED);
}

void CovMatrix(double VARIABLE_IS_NOT_USED *x, model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  DefList[SUBNR].covmatrix(sub, value);
}



int check_covmatrix(model *cov) {
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  location_type *loc = PrevLoc(cov);
  ASSERT_LOC_GIVEN;
  int err, rows, cols,
    dim = loc->timespacedim, // !! not GettLoctsdim
    total = loc->totalpoints; // !! not Gettotalpoints
  isotropy_type iso;

  assert(dim == OWNLOGDIM(0));
  if (loc->distances) {
    iso = PREVISO(0);
      iso = isCartesian(iso) ? ISOTROPIC
      : isEarth(iso) ? EARTH_ISOTROPIC
      : isSpherical(iso) ? SPHERICAL_SYMMETRIC
      : ISO_MISMATCH;
  } else {
    if (PREVXDIM(0) != PREVLOGDIM(0)) BUG;
  }

 
  if ((err = CHECK(sub, dim, OWNXDIM(0), 
		   PosDefType, KERNEL, CARTESIAN_COORD, SUBMODEL_DEP,
		   EvaluationType)) != NOERROR) { 
    if ((err = CHECK(sub, dim, OWNXDIM(0), 
		     VariogramType, XONLY, SymmetricOf(PREVISO(0)),
		     SUBMODEL_DEP,
		     EvaluationType)) != NOERROR) {
      RETURN_ERR(err);
    }
  }
  setbackward(cov, sub);  
  rows = VDIM0 = sub->vdim[0]; 
  cols = VDIM1 = sub->vdim[1]; 

  if (cov->q == NULL) {
    QALLOC(2);
    cov->q[0] = total * rows;
    cov->q[1] = total * cols;
  }

  if ((err = alloc_cov(cov, dim, rows, cols)) != NOERROR) RETURN_ERR(err);
  
 
  RETURN_NOERROR;
}



void Pseudovariogram(double VARIABLE_IS_NOT_USED *x, model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  DefList[SUBNR].pseudovariogram(sub, value);
}


void Variogram(double VARIABLE_IS_NOT_USED *x, model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;

  DefList[SUBNR].variogram(sub, value);
}

int check_vario(model *cov) {
  return check_cov_intern(cov, VariogramType,
			  GLOBAL.general.vdim_close_together,
			  true);
}


int struct_variogram(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  int err;
  model *sub,
    *next = cov->sub[0];
  location_type *loc = Loc(cov);
 
  if ((sub = get_around_max_stable(next)) == next) sub = get_around_gauss(next);
  if (sub != next) {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) RETURN_ERR(err);       
    sub = cov->key;
    SET_CALLING(sub, cov);
  }
 
  if ((err = CHECK(sub, loc->timespacedim, OWNXDIM(0), VariogramType,
		     loc->y == NULL && loc->ygr[0] == NULL ? XONLY : KERNEL,
		     SYMMETRIC, cov->vdim,
		     EvaluationType)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  if (!isnowVariogram(sub))
    SERR(sub != next ? "variogram model cannot be determined"
	 :"not a variogram model");
  
  RETURN_NOERROR;
} 
  
   
// bool close = GLOBAL.general.vdim_close_together;

void Fctn(double VARIABLE_IS_NOT_USED *X, model *cov, double *value) {
  model *sub = cov->sub[0];
  if (sub == NULL) BUG;
  FctnIntern(cov, cov, sub, value, false);
}


void FctnIntern(model *cov, model *covVdim, model *sub,
		double *value, bool ignore_y){ 
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  //  PMI0(cov); // PMI0(sub->calling);  PMI0(sub);
  assert(XDIM(PREVSYSOF(sub), 0) == LOGDIM(PREVSYSOF(sub), 0));
  assert(hasLikelihoodFrame(cov) || hasInterfaceFrame(cov) || isProcess(cov));
  
  assert(sub != NULL); 
  assert(cov->key == NULL);
  FINISH_START(cov, covVdim, ignore_y, 1);  
  assert(pgs != NULL);
  double
    *zero = ZERO(covVdim),
    *y = ygiven ? pgs->supportmin : ZERO(cov);
  Types frame = sub->frame;
  if (hasAnyEvaluationFrame(sub)) sub->frame = ShapeType;

  bool kernel = !isXonly(PREVSYSOF(sub));
  if (vdimSq > pgs->rowscols) {PMI(cov); BUG;} //
  INCLUDE_VAL;
 
#define FUNCTION FCTN(x, sub, value)
#define FUNCTION_Y NONSTATCOV(x, y, sub, value)
#define FUNCTION2 FCTN(x, sub, cross);		\
  for (v = 0; v<vdimSq; v++) { Val[v][i_row] = cross[v];}
#define FUNCTION2_Y						\
  NONSTATCOV(x, y, sub, cross);				\
  for (v = 0; v<vdimSq; v++) Val[v][i_row] = cross[v];
  ///  PMI(cov); PMI0(sub); crash();

  assert(hasFullXdim(ISO(PREVSYSOF(sub), 0))
	 ||
	 (ISO(PREVSYSOF(sub), 0) == SUBISO(0)
	  &&  isTrend(SYSTYPE(PREVSYSOF(sub), 0) )
	  ) 
	 );
 
  PERFORM(FUNCTION, FUNCTION2, FUNCTION_Y, FUNCTION2_Y);

  sub->frame = frame; // ja nicht loeschen!!
}


void FctnExtern(model *cov, model *covVdim, model *sub,
		double *value, bool ignore_y){
  Types frame = cov->frame;
  int dim =  GetLoctsdim(cov);
  if (alloc_cov(cov, dim, VDIM0, VDIM1) != NOERROR) XERR(ERRORMEMORYALLOCATION);
  cov->frame = LikelihoodType; // dummy, just to please next function
  FctnIntern(cov, covVdim, sub, value, ignore_y);
  cov->frame = frame;
  pgs_DELETE(&(cov->Spgs), cov);
}

 
int check_fctn(model *cov) {
  EXTRA_STORAGE;
#define nTypes 2  
  Types T[nTypes] = {TrendType, ShapeType},
    F[nTypes] = {TrendType, LikelihoodType};
  int i, err;
  for (i=0; i<nTypes; i++) {
    err = check_fct_intern(cov, T[i], GLOBAL.general.vdim_close_together, 
			   true, 0, 0, F[i]);
    if (err == NOERROR) RETURN_ERR(err);
  }
  // BUG vor 3.4.2018 -- warum?
  RETURN_ERR(err);
}

 



/* ****************************** */
/*             PREDICT            */
/* ****************************** */
#define PREDICT_REGISTER 0

void predict(double VARIABLE_IS_NOT_USED *x, model *predict, double *v) { 
  assert(predict != NULL && !PARAMisNULL(predict, PREDICT_REGISTER));
  model 
    *cov = KEY()[PARAM0INT(predict, PREDICT_REGISTER)],
    *sub = cov->key == NULL ? cov->sub[0] : cov->key;  
  assert(cov != NULL);
  if (v==NULL) {
    likelihood_storage *L = sub->Slikelihood;
    int store =  GLOBAL.general.set;
    GLOBAL.general.set = 0;
    listoftype *datasets = L->datasets;
    int
      vdim = VDIM0,
      ncol =  NCOL_OUT_OF(datasets),
      repet = ncol / vdim;
    GLOBAL.general.set = store;
    assert(predict->qlen > 0 && cov->q != NULL);
    predict->q[predict->qlen - 1] = repet;
    return; // EvaluateModel needs information about size
   //                      of result array, given in cov->q
  }
   
   if (SUBNR == GAUSSPROC) {
     gauss_predict(predict, cov, v);
    return;
  }

  BUG;
}


int check_predict(model *predict) {
  assert(predict != NULL);
  //PMI(predict);
  if (PARAMisNULL(predict, PREDICT_REGISTER))
    RFERROR("'register number not given.");
  model *cov = KEY()[PARAM0INT(predict, PREDICT_REGISTER)];
  location_type 
    *pred_loc = Loc(predict);
  model 
    *sub = cov->key == NULL ? cov->sub[0] : cov->key;  
  int err;
  assert(SUBNR == GAUSSPROC);
  assert(pred_loc->delete_x);
  assert(pred_loc->timespacedim == Loc(cov)->timespacedim);
  assert(pred_loc->Time == Loc(cov)->Time);
  likelihood_storage *L = sub->Slikelihood;

  if (L == NULL || L->X == NULL)
    SERR1("'%.50s' not fully initialized", NICK(cov));

  if (cov == NULL || COVNR != LIKELIHOOD_CALL || !cov->checked) 
    SERR1("'%.50s' not initialized", NICK(cov));

  if (pred_loc->y != NULL || pred_loc->ygr[0] != NULL) {
    if (predict->Sextra == NULL) // i.e. first call, so user's input
      SERR("set of y-values (kernal definition) not allowed");
  } else {  
    CONDCOV_NEW_STORAGE(predict, extra, a1);
    assert(pred_loc->delete_y); // = true;
    if (pred_loc->grid) {
      int i,
	spatialdim = pred_loc->spatialdim,
	nr = spatialdim * 3;
      double *y = (double*) MALLOC(nr * sizeof(double));
      for (i=0; i<nr; i++) y[i] = 1.0;
      assert(pred_loc->ygr[0] == NULL);
      pred_loc->ly = 3;
      if ((err = setgrid(pred_loc->ygr, y, spatialdim))!=NOERROR) 
	 RETURN_ERR(err);
      FREE(y);
      // assert(!pred_loc->Time); 
      if (pred_loc->Time) pred_loc->ygr[spatialdim] = pred_loc->T;
   } else {       
      pred_loc->ly = 1;
      // wichtig im folgenden timespacedim nicht spatialdim
      pred_loc->y = (double *) MALLOC(pred_loc->timespacedim * sizeof(double));  
      pred_loc->T[XSTART] = pred_loc->T[XSTEP] = 0.0;
      pred_loc->T[XLENGTH] = 1;
    }
    assert(cov->Sextra == NULL);
  }

  err = check_fct_intern(predict,
			 isProcess(predict->sub[0]) ? ProcessType : PosDefType,
			 GLOBAL.general.vdim_close_together, true,
			 VDIM0, 1, LikelihoodType);

  RETURN_ERR(err);
}

int struct_predict(model *predict, model VARIABLE_IS_NOT_USED  **newmodel){
  return struct_cov(predict, newmodel);
}


void range_predict(model VARIABLE_IS_NOT_USED *predict, range_type* range){
  range->min[PREDICT_REGISTER] = 0;
  range->max[PREDICT_REGISTER] = MODEL_MAX;
  range->pmin[PREDICT_REGISTER] = 0;
  range->pmax[PREDICT_REGISTER] = MODEL_MAX;
  range->openmin[PREDICT_REGISTER] = false;
  range->openmax[PREDICT_REGISTER] = false;
}
