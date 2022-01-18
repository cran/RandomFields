/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives (spectral measures, 
 tbm operators)

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nondomain models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc


 Copyright (C) 2017 -- 2017 Martin Schlather

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
#include <R_ext/Linpack.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>

#include "extern.h"
#include "questions.h"
#include "shape.h"
#include "operator.h"
#include "Coordinate_systems.h"





/* Covariate models */

int cutidx(double Idx, int len) {
  int idx = (int) ROUND(Idx);
  if (idx < 0) idx = 0;
  if (idx >= len) idx = len - 1;
  return idx;
}


#define GET_LOC_COVARIATE \
  assert(cov->Scovariate != NULL);					\
  location_type **local = P0INT(COVARIATE_RAW) || PisNULL(COVARIATE_X)	\
    ? PLoc(cov) : cov->Scovariate->loc;					\
  assert(local != NULL);						\
  location_type *loc =  LocLoc(local);					\
  assert(loc != NULL)
 

int get_index(double *x, model *cov) {
  GET_LOC_COVARIATE;

  //  PMI(cov);
  //  printf("%d %d to=%d last:%d==%d loc=%d\n",OWNTOTALXDIM, PREVTOTALXDIM,
  //	 total_logicaldim(SYSOF(cov)), PREVLASTSYSTEM, OWNLASTSYSTEM,
  //	 GetLoctsdim(cov));
  
  double X[2];
  int d, idx, i,
    nr = 0,
    cummul = 1.0,  
    ntot = loc->totalpoints,
    dim = OWNTOTALXDIM;
  assert(OWNTOTALXDIM == PREVTOTALXDIM && PREVLASTSYSTEM  == OWNLASTSYSTEM &&
	 OWNTOTALXDIM == total_logicaldim(SYSOF(cov)) &&   
	 loc->timespacedim == OWNTOTALXDIM);		       
  
  //printf("grid = %d\n", loc->grid);

  if (loc->grid) {
    for (d = 0; d<dim; d++) {
      int len = loc->xgr[d][XLENGTH];
      double
	step = loc->xgr[d][XSTEP];                
      if (d > 1 || !isAnySpherical(OWNISO(0))) {
	idx = cutidx((x[d] - loc->xgr[d][XSTART]) / step, len);
      } else {
	if (d == 0) { // to do: report technique?
	  double full, half, y[2];
	  int idx2,
	    xdim = 2;
	  for (i=0; i<xdim; i++) y[i] = loc->xgr[i][XSTART];
	  if (isSpherical(OWNISO(0))) {
	    full = M_2_PI;
	    half = M_PI;
	    if (GLOBAL.coords.polar_coord) NotProgrammedYet("");
	  } else if (isEarth(OWNISO(0))) {
	    full = 360.0;
	    half = 180.0;
	  } else BUG;
	  statmod2(y, full, half, X);
	  
	  idx = cutidx((x[0] - X[0]) / step, len);
	  double X0 = X[0] + full * (2.0 * (double) (x[0] > X[0]) - 1);
	  idx2 = cutidx((x[0] - X0) / step, len);
	  if (FABS(x[0] - (X0 + idx2 * step)) <
	      FABS(x[0] - (X[0] + idx * step))) idx = idx2;
	} else { 
	  assert(d==1);
	  idx = cutidx((x[d] - X[1]) / step, len);
	}
      }      
      nr += cummul * idx;

      //  printf("nr = %d\n", nr);

      cummul *= len;  
    }           
  } else { // to do: effizienterer Zugriff ueber Kaestchen eines Gitters, 
    // in dem die jeweiligen Punkte gesammelt werden. Dann muessen nur
    // 3^d Kaestchen durchsucht werden.
    
    model *next = cov->sub[0];
    double distance,
      mindist = RF_INF,
      *given = loc->x;
    for (i=0; i<ntot; i++, given+=dim) {
      assert(next->checked);
      NONSTATCOV(x, given, next, &distance);
      if (distance < mindist) {
	mindist = distance;
	nr = i;
      }
    }
  }
  return nr;
}


void kappa_covariate(int i, model VARIABLE_IS_NOT_USED *cov,
		     int *nr, int *nc){
  *nc = *nr = i <= COVARIATE_X  || i == COVARIATE_FACTOR ? SIZE_NOT_DETERMINED
    : i <= COVARIATE_FACTOR ? 1 
    : -1;
}


void covariate(double *x, model *cov, double *v){
  // ACHTUNG!! FALLS NAs verdeckt da sind, i.e. COVARIATE_ADDNA = TRUE:
  // HIER WIRD GETRICKST: vdim[0]=1 fuer Kovarianz, das hier nur 
  //                      vim[0] abgefragt wird und de facto univariat ist
  // und vdim[1]=Anzahl der covariaten fuer matrix berechnung, da
  //     fctnintern das produkt betrachtet und somit die dimensionen der
  //     design matrix reflektiert.
  // Fuer COVARIATE_ADDNA = FALSE haben wir ganz normals Verhalten
   GET_LOC_COVARIATE;
  double 
    *p = LP(COVARIATE_C);
  bool addna = (bool) cov->q[1];
  assert(cov->vdim[!addna] == 1);
  int  nr,
    dim = ANYDIM,
    vdim = cov->vdim[addna],
    ntot = loc->totalpoints;

  //  PMI0(cov);
  
  if (hasAnyEvaluationFrame(cov)) {
    for (int i=0; i<vdim; v[i++] = 0.0);
    return;
  } else assert(isnowTrend(cov));

  if (P0INT(COVARIATE_RAW)) {
    // should happen only in a matrix context!
    nr = (int) x[dim];
    if (nr * vdim >= LNROW(COVARIATE_C) * LNCOL(COVARIATE_C))
      ERR0("illegal access -- 'raw' should be FALSE");
  } else { 
    // interpolate: here nearest neighbour/voronoi
    nr = get_index(x, cov);

    //   printf("x=%10g %10g: %d\n",x[0], x[1], nr);
    
  }
  //  printf("%10g %10g %10g %d nr=%d q=%10g %d vdim=%d\n", x[0],x[1],x[2], dim, nr, cov->q[0], GLOBAL.general.vdim_close_together, vdim);
  
  if (cov->q[0] == 0) {
    if (GLOBAL.general.vdim_close_together) {
      p += nr * vdim;
      for (int i=0; i<vdim; i++, p++) v[i] = *p;
    } else {
      p += nr;
      for (int i=0; i<vdim; i++, p+=ntot) v[i] = *p;
      //      printf("%10g\n", v[0]);
      //     APMI(cov->calling->calling);
    }
  } else {
    if (GLOBAL.general.vdim_close_together) {
      double dummy = 0.0;
      p += nr * vdim;      
      for (int i=0; i<vdim; i++, p++) dummy += *p * P(COVARIATE_FACTOR)[i];
      *v = dummy;
   } else {
      p += nr;
      //   PMI(cov);
      //    printf("%d %d %d\n", nr, vdim, ntot);
      for (int i=0; i<vdim; i++, p+=ntot) v[i] = *p * P(COVARIATE_FACTOR)[i];
    }
  }
  //printf("%d %d\n", (int) v[0], (int) v[1]);
}



int check_fix_covariate(model *cov,  location_type ***local){
  int err,
    store = GLOBAL.general.set;
  GLOBAL.general.set = 0;
  bool
    globalXT = PisNULL(COVARIATE_X);
  coord_sys_enum ccs = GLOBAL.coords.coord_system;
 
  if (cov->Scovariate != NULL &&
      cov->Scovariate->matrix_err != MATRIX_NOT_CHECK_YET && 
      cov->Scovariate->matrix_err != NOERROR)
    return cov->Scovariate->matrix_err; // #define RETURN

  assert(cov->calling != NULL); // globalXT && cov->calling != NULL && isnowTrend(cov)
  kdefault(cov, COVARIATE_RAW, globalXT 
	   //&& cov->calling != NULL && isnowTrend(cov)
	   ); 

  bool raw = P0INT(COVARIATE_RAW);
  if (raw) {
    model *calling = cov->calling;
    assert(calling != NULL);    
    if (!globalXT || (// calling != NULL &&
		      isDollar(calling) && !hasVarOnly(calling) &&
		      calling->kappasub[DSCALE] != cov) 
	) {
      //printf("%d %d\n", isDollar(calling), !hasVarOnly(calling));
      // APMI0(cov); 
      SERR("if 'raw' then none of {'x', 'T', 'Aniso', 'proj', 'scale'} may be given");
    }
    assert(cov->ownloc == NULL);  
    if (cov->Scovariate == NULL) NEW_STORAGE(covariate); 
    *local = PLoc(cov);
  } else {
    if ((ccs == cartesian && !equalsCartCoord(OWNISO(0))) ||
	(ccs == earth && !equalsEarthCoord(OWNISO(0))) ||
	(ccs == sphere && !equalsSphericalCoord(OWNISO(0))))
      SERR2("'%.50s' not the global coordinate system ('%.50s')",
	    ISO_NAMES[OWNISO(0)], COORD_SYS_NAMES[GLOBAL.coords.coord_system]);
    
    if (globalXT) {
      if (cov->Scovariate == NULL) NEW_STORAGE(covariate); 
      *local = PLoc(cov);
    } else { // neither raw nor globalXT
      bool doset = cov->Scovariate == NULL;
      if (!doset) {
	location_type *loc = LocLoc(cov->Scovariate->loc);
	assert(loc != NULL);
	doset = Loc(cov)->spatialdim != loc->spatialdim || 
	  Loc(cov)->xdimOZ != loc->xdimOZ;
      }
      if (doset) {
	NEW_STORAGE(covariate);
	covariate_storage *S = cov->Scovariate;
	GLOBAL.general.set = store;
	S->loc = loc_set(PVEC(COVARIATE_X), false);
	GLOBAL.general.set = 0;
 	assert(S->loc != NULL);
	
	if (S->loc[0]->timespacedim != OWNLOGDIM(0)) 
	  SERR1("number of columns of '%.50s' do not equal the dimension",
		KNAME(COVARIATE_X));
     }
    *local = cov->Scovariate->loc;
    } // neither raw nor globalXT
  }
    
  ASSERT_ONESYSTEM;
  if (!raw) {
    if (cov->sub[0] == NULL) {
      addModel(cov, 0, TRAFO);
      kdefault(cov->sub[0], TRAFO_ISO, IsotropicOf(PREVISO(0)));
    } else {
      model *next = cov->sub[0];
      PARAMINT(next, TRAFO_ISO)[0] = IsotropicOf(PREVISO(0));
    }
 
    if ((err = CHECK(cov->sub[0], OWNLOGDIM(0), OWNXDIM(0), ShapeType,
		     KERNEL, OWNISO(0), 1, cov->frame)) != NOERROR) {
      //APMI(cov);
      RETURN_ERR(err);
    }
  }
    
  RETURN_NOERROR;
}



int checkcovariate(model *cov){
 assert(cov != NULL);
  int store = GLOBAL.general.set;
  GLOBAL.general.set = 0;
  int len,
    vdim = UNSET, 
    err = NOERROR;
  location_type **local = NULL;

  
  // double value = 1.0;
  //  covariate_storage *S;
  //  bool addna = !PisNULL(COVARIATE_ADDNA) && P0INT(COVARIATE_ADDNA);

  bool addna; 
  if (cov->q == NULL) {
    addna = (!PisNULL(COVARIATE_ADDNA) && P0INT(COVARIATE_ADDNA)) || 
      (!PisNULL(COVARIATE_FACTOR) && ISNAN(P0(COVARIATE_FACTOR)));
    QALLOC(2);
    // cov->q[0] == 0 iff vdim vectors are returned
    cov->q[1] = addna;
  } else addna = (bool) cov->q[1];

 
  if ((err = check_fix_covariate(cov, &local)) != NOERROR) goto ErrorHandling;
  assert(local != NULL);
 
  len = local[0]->len;
  //S = cov->Scovariate;  assert(S != NULL);
  if (len <= 0) BUG;
  for (GLOBAL.general.set=0;  GLOBAL.general.set<len;  GLOBAL.general.set++) {
    int
      ndata = LNROW(COVARIATE_C) * LNCOL(COVARIATE_C),
      ntot = LocLoc(local)->totalpoints;
    if (GLOBAL.general.set == 0) {
      vdim = ndata/ntot;
      //      PMI0(cov);
      //printf("vdim=%d %d %d\n", vdim, ndata, ntot); BUG;
    }

    //PMI(cov->calling);
    
    if (ntot <= 0) GERR("no locations given");
    if (vdim * ntot != ndata) {
      //BUG;
      GERR3("Number of data (%d) not a multiple of\n'the number (%d) of locations times the number of components (%d)'.", ndata, ntot, vdim);      
    }
  }
  assert(vdim > 0);

  
  if (addna) {
    if (!isnowTrend(cov))
      GERR2("'%.50s' can be true only if '%.50s' is used as a trend",
	    KNAME(COVARIATE_ADDNA), NAME(cov));
    
    model *calling = cov->calling;
    if (PisNULL(COVARIATE_FACTOR)) {
      while (calling != NULL && CALLINGNR != LIKELIHOOD_CALL &&
	     CALLINGNR != LINEARPART_CALL) {
	calling = CALLINGNR == PLUS ? calling->calling : NULL;	
      }
    }
    //   if (calling == NULL)  {
    //  GERR1("%.50s is used with NAs outside a trend definition.", NAME(cov));
    // }
    if (PisNULL(COVARIATE_FACTOR)) PALLOC(COVARIATE_FACTOR, vdim, 1);
    for (int i=0; i<vdim; i++) P(COVARIATE_FACTOR)[i] = RF_NAN;
  }

  /*
   if (PisNULL(COVARIATE_FACTOR)) {
    PALLOC(COVARIATE_FACTOR, vdim, 1);
    for (int i=0; i<vdim; i++) P(COVARIATE_FACTOR)[i] = value;
  }
  
  kdefault(cov, COVARIATE_ADDNA, 
  !PisNULL(COVARIATE_FACTOR) && ISNAN(P0(COVARIATE_FACTOR)));

  Brasilien:
  kdefault(cov, COVARIATE_ADDNA, !PisNULL(COVARIATE_FACTOR) &&
	   (ISNA(P0(COVARIATE_FACTOR)) || ISNAN(P0(COVARIATE_FACTOR))));
  */
   
  cov->q[0] = ((bool) cov->q[1]) || PisNULL(COVARIATE_FACTOR) ? 0 : vdim;
  cov->vdim[!addna] = 1; 
  cov->vdim[addna] = cov->q[0] == 0.0 ? vdim : 1;

  assert( VDIM0 > 0 && VDIM1 >0);
  

  if (hasAnyEvaluationFrame(cov) && VDIM0 != 1)
    GERR1("'%.50s' used in a wrong context", NICK(cov));

  if ((err = checkkappas(cov, false)) != NOERROR) goto ErrorHandling;
    
  cov->mpp.maxheights[0] = RF_NA;
  
 ErrorHandling: 
  GLOBAL.general.set = store;

  PFREE(COVARIATE_ADDNA); 
  RETURN_ERR(err);
}


void rangecovariate(model *cov, range_type *range){
  rangefix(cov, range);

  booleanRange(COVARIATE_ADDNA);

  range->min[COVARIATE_FACTOR] = RF_NEGINF;
  range->max[COVARIATE_FACTOR] = RF_INF;
  range->pmin[COVARIATE_FACTOR] = - 1e10;
  range->pmax[COVARIATE_FACTOR] = 1e10;
  range->openmin[COVARIATE_FACTOR] = true;
  range->openmax[COVARIATE_FACTOR] = true;
}


void kappa_fix(int i, model VARIABLE_IS_NOT_USED *cov,
		     int *nr, int *nc){
    *nc = *nr = i <= FIXCOV_X ? 0 
      : i <= FIXCOV_RAW ? 1 
      : -1;
}



void fix(double *x, double *y, model *cov, double *v){
  GET_LOC_COVARIATE;
   int nrx, nry,
     dim = OWNTOTALXDIM,
     ntot = loc->totalpoints,
    vdim = VDIM0;
  double
    *p = LP(FIXCOV_M);

  assert(VDIM0 == VDIM1);
  if (P0INT(FIXCOV_RAW)) {
    // should happen only in a matrix context!
    nrx = (int) x[dim];
    nry = y == NULL ? x[dim+1] : (int) y[dim];
    if (nrx * vdim >= LNROW(FIXCOV_M) ||
	nry * vdim >= LNCOL(FIXCOV_M))
      ERR0("illegal access -- 'raw' should be FALSE");
  } else {    
    nrx = get_index(x, cov);
    nry = get_index(y, cov);
  }

  int i, j, k, 
    ntotvdim = ntot * vdim;
  if (GLOBAL.general.vdim_close_together) {
    p += (ntotvdim * nry + nrx) * vdim;
    for (k=i=0; i<vdim; i++, p += ntotvdim) {
      double *q = p;
      for (j=0; j<vdim; j++, q++) v[k++] = *q;
    }
  } else {
    int ntotSqvdim = ntotvdim * ntot;
    p += ntotvdim * nry + nrx;
    for (k=i=0; i<vdim; i++, p += ntotSqvdim) {
      double *q = p;
      for (j=0; j<vdim; j++, q+=ntot) v[k++] = *q;
    }
  }
}


void fixStat(double *x, model *cov, double *v){
  assert(P0INT(FIXCOV_RAW));
  assert(Loc(cov)->distances);
  fix(x, NULL, cov, v);
}





/*

bool allowedDfix(model *cov) {
  kdefault(cov, FIXCOV_RAW, PisNULL(FIXCOV_X));
  // && cov->calling != NULL && isnowTrend(cov)); 
  if (P0INT(FIXCOV_RAW)) return allowedDtrue;
  for(int i=FIRST_DOMAIN; i<= LAST_DOMAINUSER; i++) cov->allowedD[i] = false;
  cov->allowedD[KERNEL] = true;
  return false;
}

*/

bool allowedIfix(model *cov) {
  model *dummy = cov;
  location_type *loc = Loc(dummy);
  while (loc == NULL && dummy->calling != NULL) {
    dummy = dummy->calling;
    loc = Loc(dummy);
  }
  if (loc == NULL) BUG;
  bool dist = loc->distances,
    *I = cov->allowedI;
  //printf("name dummy = %.50s dist = %d\n", NAME(dummy), dist);
  kdefault(cov, FIXCOV_RAW, PisNULL(FIXCOV_X));
  // && cov->calling != NULL && isnowTrend(cov)); 
  for(int i=FIRST_ISOUSER; i<=LAST_ISOUSER; i++) I[i] = false;
  if (dist) I[ISOTROPIC] = I[EARTH_ISOTROPIC] = true;
    // not spherical since it involves a transformation that looses i_nrow
  else I[CARTESIAN_COORD] = I[EARTH_COORD] = true;
  return false;
}

/*
bool setfix(model *cov) {
  kdefault(cov, FIXCOV_RAW,
	   PisNULL(FIXCOV_X) && cov->calling != NULL && isnowTrend(cov)); 
  if (P0INT(FIXCOV_RAW)) {
    isotropy_type iso = CONDPREVISO(0); 
    if (!isFixed(iso)) return false;
    set_iso(OWN, 0, iso);
    set_dom(OWN, 0, PREVDOM(0));    
  } else {
    set_iso(OWN, 0, CoordinateSystemOf(PREVISO(0)));
    set_dom(OWN, 0, KERNEL);
  }
  return true;
}
*/

Types Typefix(Types required, model *cov, isotropy_type requ_iso){
  if (isBad(TypeConsistency(required, OWNTYPE(0)))) return BadType;
  if (P0INT(FIXCOV_RAW)) return required;
  return equalsCartCoord(requ_iso) ? required : BadType;
}


int checkfix(model *cov){
  assert(cov != NULL);
  int store = GLOBAL.general.set;
  GLOBAL.general.set = 0;
  int 
    vdim = UNSET, 
    vdimSq = UNSET,
    err = NOERROR;
  location_type **local = NULL;
  bool isAnyIso = isAnyIsotropic(OWNISO(0));
  assert(!isAnyIso || Loc(cov)->distances);

  if ((err = check_fix_covariate(cov, &local)) != NOERROR) goto ErrorHandling;
  if (!P0INT(FIXCOV_RAW)) {
    if (isAnyIso)
      SERR1("distances only allowed for '%.50s' = TRUE", KNAME(FIXCOV_RAW));
    if (!(equalsKernel(OWNDOM(0)) && equalsCartCoord(OWNISO(0))))
      SERR2("Model only allowed within positive definite kernels, not within positive definite functions. (Got %.50s, %.50s.)",  DOMAIN_NAMES[OWNDOM(0)], ISO_NAMES[OWNISO(0)]);
  }
  assert(local != NULL);
  covariate_storage *S;
  int len;
  len = local[0]->len;
  S = cov->Scovariate;
  assert(S != NULL);
  if (len <= 0) BUG;
  for (GLOBAL.general.set=0;  GLOBAL.general.set<len;  GLOBAL.general.set++) {
    int
      ndata = LNROW(FIXCOV_M) * LNCOL(FIXCOV_M),
      ntot = LocLoc(local)->totalpoints;
    if (GLOBAL.general.set == 0) {
      vdim = (int) SQRT((double) ndata / (double) (ntot * ntot));
      vdimSq = vdim * vdim;
    }

    if (ntot <= 0) GERR("no locations given");
    if (cov->nrow[FIXCOV_M] != cov->ncol[FIXCOV_M])
      GERR("square matrix expected");
    if (vdimSq * ntot * ntot != ndata)
      GERR3("number of data (%d) not a multiple of the number of locations (%d x %d)^2", ndata, ntot, vdim);      
  }
  assert(vdim > 0);
  VDIM0 = VDIM1 = vdim;

  if ((err = checkkappas(cov)) != NOERROR) goto ErrorHandling;
     
  if (vdim == 1 && len == 1) {
    GLOBAL.general.set = 0;
    double *c = LP(FIXCOV_M);   
    int i,
      ntot = LocLoc(local)->totalpoints;
    for (i=0; i<ntot; i++) if (c[i] != 0.0) break;
    cov->ptwise_definite = pt_zero;
    if (i<ntot) {
      if (c[i] > 0.0) {
	cov->ptwise_definite = pt_posdef;
        for (i++; i<ntot; i++)
	  if (c[i] < 0.0) {
	    cov->ptwise_definite = pt_indef;
	    break;
	  }
      } else { // < 0.0
	cov->ptwise_definite = pt_negdef;
        for (i++; i<ntot; i++)
	  if (c[i] > 0.0) {
	    cov->ptwise_definite = pt_indef;
	    break;
	  }
      }
    }
    
  } else cov->ptwise_definite = pt_unknown; // to do ?! alle Matrizen ueberpruefen...

     
  if (S->matrix_err == MATRIX_NOT_CHECK_YET) {
    S->matrix_err = NOERROR;
    for (GLOBAL.general.set=0; GLOBAL.general.set<len; GLOBAL.general.set++){
      int j, k,
	nrow = LNROW(FIXCOV_M),
	ncol = LNCOL(FIXCOV_M);
      double
	*C =  LP(FIXCOV_M);
      if (nrow != ncol || nrow == 0) {
	S->matrix_err = err = ERROR_MATRIX_SQUARE; 
	goto ErrorHandling;
      }    

      for (k=0; k<nrow; k++) {	
	for (j=k+1; j<ncol; j++) {
	  if (C[k + nrow * j] != C[j + nrow * k]) 
	    GERR("matrix not symmetric");
	}
      }
      
   if (nrow < 3000) {
	if (!Ext_is_positive_definite(C, nrow)) {
	  S->matrix_err = err = ERROR_MATRIX_POSDEF;
	  goto ErrorHandling;
	}
      } else {
	if (len==1) 
	  PRINTF("covariance matrix is large, hence not verified to be positive definite.");
	else PRINTF("covariance matrix of %d%.50s set is large, hence not verified to be positive definite.", GLOBAL.general.set + 1, TH(GLOBAL.general.set + 1));
      }
      
    }
    if (err != NOERROR) goto ErrorHandling;
  }
       
  cov->matrix_indep_of_x = true;
  cov->mpp.maxheights[0] = RF_NA;
  
 ErrorHandling: 
  GLOBAL.general.set = store;

  RETURN_ERR(err);
}


void rangefix(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[FIXCOV_M] = RF_NEGINF;
  range->max[FIXCOV_M] = RF_INF;
  range->pmin[FIXCOV_M] = - 1e10;
  range->pmax[FIXCOV_M] = 1e10;
  range->openmin[FIXCOV_M] = true;
  range->openmax[FIXCOV_M] = true;

  range->min[FIXCOV_X] = RF_NA;
  range->max[FIXCOV_X] = RF_NA;
  range->pmin[FIXCOV_X] = RF_NA;
  range->pmax[FIXCOV_X] = RF_NA;
  range->openmin[FIXCOV_X] = true;
  range->openmax[FIXCOV_X] = true;

  booleanRange(FIXCOV_RAW);

}




// ----------------------------------------------------------------------
// ball

void ball(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) {
  //APMI(cov);
  // isotropic function, expecting only distance; BALL_RADIUS=1.0
  assert(*x >= 0);
  *v = (double) (*x <= BALL_RADIUS);
}

void Inverseball(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){
  //  crash();
  *v = *x > 1.0 ? 0.0 : BALL_RADIUS;
}


int struct_ball(model *cov, model **newmodel){
  ASSERT_NEWMODEL_NOT_NULL;
  if (hasSmithFrame(cov)) {
    return addUnifModel(cov, BALL_RADIUS, newmodel);
  } else {
    ILLEGAL_FRAME;
  }
  RETURN_NOERROR;
}

int init_ball(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {

  //  TREE0(cov);  PMI0(cov);
  
  assert(s != NULL);
  if (hasAnyEvaluationFrame(cov)) RETURN_NOERROR;
  
  if (hasSmithFrame(cov) || hasAnyPoissonFrame(cov)) {
    cov->mpp.maxheights[0] = 1.0;
    
    if (cov->mpp.moments >= 1) {
      cov->mpp.mM[1] = cov->mpp.mMplus[1] =
	VolumeBall(OWNLOGDIM(0), BALL_RADIUS);
      int i;
      for (i=2; i<=cov->mpp.moments; i++)  
	cov->mpp.mM[i] = cov->mpp.mMplus[i] = cov->mpp.mM[1];
    }
  }

  else if (hasRandomFrame(cov)) { RETURN_NOERROR; }
 
  else ILLEGAL_FRAME;


  RETURN_NOERROR;
}


void do_ball(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s) { 
  assert(s != NULL);
 
}








// ----------------------------------------------------------------------
// Poisson polygons


double meanVolPolygon(int dim, double beta) {
  double kd = VolumeBall(dim, 1.0),
    kdM1 = VolumeBall(dim-1, 1.0);
  return intpow(dim * kd / (kdM1 * beta), dim) / kd;
}

void Polygon(double *x, model VARIABLE_IS_NOT_USED *cov, double *v) { 
  polygon_storage *ps = cov->Spolygon;
  assert(ps != NULL);
  *v = (double) isInside(ps->P, x);
}
 
void Inversepolygon(double VARIABLE_IS_NOT_USED *x, model *cov, double *v){
  polygon_storage *ps = cov->Spolygon;
  int d,
    dim = OWNLOGDIM(0);

  if (ps == NULL) {
    *v = RF_NA;
    return;
  }
  polygon *P = ps->P;
  if (P != NULL) {
    double max = 0.0;
    for (d=0; d<dim; d++) {
      double y = FABS(P->box0[d]);
      if (y > max) max = y;
      y = FABS(P->box1[d]);
      if (y > max) max = y;
    }
  } else {
    BUG;

    *v = POW(meanVolPolygon(dim, P0(POLYGON_BETA)) / VolumeBall(dim, 1.0), 
	     1.0/dim);
    // to do kann man noch mit factoren multiplizieren, siehe
    // strokorb/schlather, max
    // um unabhaengigkeit von der Dimension zu erlangen
  }
}

void InversepolygonNonstat(double VARIABLE_IS_NOT_USED *v, model *cov,
			   double *min, double *max){
  polygon_storage *ps = cov->Spolygon;
  int d,
    dim = OWNLOGDIM(0);
  assert(ps != NULL);
  if (ps == NULL) {
    for (d=0; d<dim; d++) min[d] = max[d] = RF_NA;
    return;
  }
  polygon *P = ps->P;
  if (P != NULL) {
     for (d=0; d<dim; d++) {
      min[d] = P->box0[d];
      max[d] = P->box1[d];   
    }
  } else { // gibt aquivalenzradius eines "mittleres" Polygon zurueck

    BUG;

    double size = POW(meanVolPolygon(dim, P0(POLYGON_BETA)) / 
		      VolumeBall(dim, 1.0), 1.0/dim);    
    // to do kann man noch mit factoren multiplizieren, siehe
    // strokorb/schlather, max-stabile Modelle mit gleichen tcf
    for (d=0; d<dim; d++) {
      min[d] = -size;
      max[d] = size;
    }
  }
}

int check_polygon(model *cov) {
  int err,
    dim = ANYOWNDIM;
  
  assert(dim <= MAXDIM_POLY);
  if (dim != 2)
    SERR("random polygons only defined for 2 dimensions"); // to do
  kdefault(cov, POLYGON_BETA, 1.0);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  cov->randomkappa = true;
  RETURN_NOERROR;
}

void range_polygon(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[POLYGON_BETA] = 0;
  range->max[POLYGON_BETA] = RF_INF;
  range->pmin[POLYGON_BETA] = 1e-3;
  range->pmax[POLYGON_BETA] = 1e3;
  range->openmin[POLYGON_BETA] = true;
  range->openmax[POLYGON_BETA] = true;
}


int struct_polygon(model VARIABLE_IS_NOT_USED *cov,
		   model VARIABLE_IS_NOT_USED **newmodel){
  /*
  ASSERT_NEWMODEL_NOT_NULL;
  double beta = P0(POLYGON_BETA);
  if (hasAny ShapeFrame(cov)) {
    double 
      dim = OWNLOGDIM(0),
      safety = P0(POLYGON_SAFETY), // to do : zhou approach !
      equiv_cube_length = POW(beta, -1.0/dim);
    return addUnifModel(cov,  // to do : zhou approach !
			0.5 * equiv_cube_length * safety,
			newmodel);
  } 

  else if (hasRandomFrame(cov)) RETURN_NOERROR;

  else {
    ILLEGAL_FRAME;
  }
  */
  BUG;

  RETURN_NOERROR;
}

polygon_storage *create_polygon() {
  polygon_storage *ps;
  if ((ps = (polygon_storage*) MALLOC(sizeof(polygon_storage))) == NULL)
    return NULL;
  if ((ps->P = (polygon*)  MALLOC(sizeof(polygon))) == NULL) {
    FREE(ps);
    return NULL;
  }
  polygon_NULL(ps);
  return ps;
}

int init_polygon(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int i, err,
    dim = OWNLOGDIM(0);
  double beta = P0(POLYGON_BETA),
    lambda = beta; // / VolumeBall(dim - 1, 1.0) * (dim * VolumeBall(dim, 1.0));
  assert(s != NULL);
  if (cov->Spolygon == NULL) {
    if ((cov->Spolygon = create_polygon()) ==NULL) RETURN_ERR(ERRORMEMORYALLOCATION);
  }
   
  // nicht nur aber auch zur Probe
  struct polygon_storage *ps = cov->Spolygon;
  assert(ps != NULL && ps->P != NULL);
  
  if (!false) {
    freePolygon(ps->P); 
    if ((err=rPoissonPolygon(ps->P, lambda, true)) != NOERROR)
      SERR1("poisson polygon cannot be simulated (error=%d)", err);
  } else {
    BUG; // valgrind memalloc loss ! + lange Laufzeiten ?!
    if ((err=rPoissonPolygon2(ps, lambda, true)) != NOERROR)
      SERR1("Poisson polygon cannot be simulated (error=%d)", err);
  }
  
  if (hasSmithFrame(cov)) {
    double c = meanVolPolygon(dim, beta);
    cov->mpp.maxheights[0] = 1.0; 
   for (i=1; i<=cov->mpp.moments; i++) cov->mpp.mM[i] = cov->mpp.mMplus[i] = c;	
  }

  else ILLEGAL_FRAME;

  RETURN_NOERROR;
}




// ----------------------------------------------------------------------
//
#define RATIONAL_A 0
#define RATIONAL_a 1
void kappa_rational(int i, model *cov, int *nr, int *nc){
  *nc = (i == RATIONAL_A) ? OWNLOGDIM(0) : 1;
  *nr = (i == RATIONAL_A) ? OWNLOGDIM(0) : (i==RATIONAL_a) ? 2 : -1;
}
void minmaxEigenrational(model *cov, double *mm) {
  double *a = P(RATIONAL_a);
  if (a[0] < a[1]) {
    mm[0] = a[0];
    mm[1] = a[1];
  } else {
    mm[0] = a[1];
    mm[1] = a[0];
  }
}
double maxEigenrational(model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *mm) {
  double *a = P(RATIONAL_a);
  return (a[0] > a[1]) ? a[0] : a[1];
}
void rational(double *x, model *cov, double *v) {
  int i, k, j, 
    dim = OWNLOGDIM(0);
  double nu,
    *A = P(RATIONAL_A),
    *a = P(RATIONAL_a);
  nu = 0.0;
  for (k=0, i=0; i<dim; i++) {
    double xTC;
    xTC =  0.0;
    for (j=0; j<dim; j++) {
      xTC += x[j] * A[k++];
    }
    nu += xTC * xTC;
  }
  *v = (a[0] + a[1] * nu) / (1 + nu);
}
 
int checkrational(model *cov){
  int err;
  if (cov->nrow[RATIONAL_a] == 1) {
    double dummy = P0(RATIONAL_a);
    PFREE(RATIONAL_a);
    PALLOC(RATIONAL_a, 2, 1);
    P(RATIONAL_a)[0] = dummy;
    P(RATIONAL_a)[1] = 0.0;
  }
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  cov->mpp.maxheights[0] =  P(RATIONAL_a)[0] > P(RATIONAL_a)[1] 
    ? P(RATIONAL_a)[0] : P(RATIONAL_a)[1];
  RETURN_NOERROR;
}

void rangerational(model VARIABLE_IS_NOT_USED *cov, range_type* range){

  range->min[RATIONAL_A] = RF_NEGINF;
  range->max[RATIONAL_A] = RF_INF;
  range->pmin[RATIONAL_A] = - 1e10;
  range->pmax[RATIONAL_A] = 1e10;
  range->openmin[RATIONAL_A] = true;
  range->openmax[RATIONAL_A] = true;

  range->min[RATIONAL_a] = 0.0;
  range->max[RATIONAL_a] = RF_INF;
  range->pmin[RATIONAL_a] = 0.0;
  range->pmax[RATIONAL_a] = 10;
  range->openmin[RATIONAL_a] = false;
  range->openmax[RATIONAL_a] = true;
}



// ----------------------------------------------------------------------
//
// Sigma(x) = diag>0 + A'xx'A
#define EAXXA_E 0
#define EAXXA_A 1
#define ETAXXA_ALPHA 2
void kappa_EAxxA(int i, model *cov, int *nr, int *nc){
  *nc = (EAXXA_A == i) ? OWNLOGDIM(0) : 1;
  *nr = i < DefList[COVNR].kappas ? OWNLOGDIM(0) : -1;
}
void EAxxA(double *x, model *cov, double *v) {
  int d, k, j, 
    dim = OWNLOGDIM(0);
  double xA[EaxxaMaxDim],
    *E = P(EAXXA_E),
    *A = P(EAXXA_A);
  for (k=0, d=0; d<dim; d++) {
    xA[d] =  0.0;
    for (j=0; j<dim; j++) {
      xA[d] += x[j] * A[k++];
    }
  }
  for (k=d=0; d<dim; d++) {
    double xAd = xA[d];
    for (j=0; j<=d; j++) {
      v[k++] = xAd * xA[j];
    }
    v[k-1] += E[d];
    for ( ; j<dim; j++) {
      v[k++] = xAd * xA[j];
    }
  }
}

void minmaxEigenEAxxA(model *cov, double *mm) {
  double 
    *E = P(EAXXA_E);
  int i,
    dim = OWNLOGDIM(0);
  for (mm[0] = RF_INF, mm[1]=RF_NEGINF, i=0; i<dim; i++) {
    if (E[i] < mm[0]) mm[0] = E[i];
    if (E[i] > mm[1]) mm[1] = E[i];
  }
}
 
int checkEAxxA(model *cov){
  int err;
  if (OWNXDIM(0) > EaxxaMaxDim)
      SERR2("For technical reasons max. dimension for ave is %d. Got %d.", 
	    EaxxaMaxDim, OWNXDIM(0));
    
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);

  VDIM0 = VDIM1 = OWNLOGDIM(0);
  cov->mpp.maxheights[0] = RF_NA;
 RETURN_NOERROR;
}
 
void rangeEAxxA(model VARIABLE_IS_NOT_USED *cov, range_type* range){

  range->min[EAXXA_E] = 0.0;
  range->max[EAXXA_E] = RF_INF;
  range->pmin[EAXXA_E] = 0.0001;
  range->pmax[EAXXA_E] = 10;
  range->openmin[EAXXA_E] = true;
  range->openmax[EAXXA_E] = true;

  range->min[EAXXA_A] = RF_NEGINF;
  range->max[EAXXA_A] = RF_INF;
  range->pmin[EAXXA_A] = - 1e10;
  range->pmax[EAXXA_A] = 1e10;
  range->openmin[EAXXA_A] = true;
  range->openmax[EAXXA_A] = true;
}






// ----------------------------------------------------------------------
//
// Sigma(x) = diag>0 + A'xx'A
void kappa_EtAxxA(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  int logicaldim = 3; 
  *nc = (i == EAXXA_A) ? logicaldim : 1;
  *nr = (i == EAXXA_E || i==EAXXA_A) ? logicaldim : (i==ETAXXA_ALPHA) ? 1 : -1;
}
void EtAxxA(double *x, model *cov, double *v) {
  int d, k, j, 
    dim = OWNLOGDIM(0),
    time = dim - 1;
  double xAR[EaxxaMaxDim], R[9],
    *E = P(EAXXA_E),
    *A = P(EAXXA_A),
    phi = P0(ETAXXA_ALPHA),
    c =  COS(phi * x[time]),
    s = SIN(phi * x[time]); 
     
  R[0] = R[4] = c;
  R[1] = s;
  R[3] = -s;
  R[2] = R[5] = R[6] = R[7] = 0.0;
  R[8] = 1.0;
 
  {
    double xA[EaxxaMaxDim];
    for (k=0, d=0; d<dim; d++) {
      xA[d] =  0.0;
      for (j=0; j<dim; j++) {
	xA[d] += x[j] * A[k++];
      }
    }
    for (k=0, d=0; d<dim; d++) {
      xAR[d] =  0.0;
      for (j=0; j<dim; j++) {
	xAR[d] += xA[j] * R[k++];
      }
    }
  }


  for (k=d=0; d<dim; d++) {
    double xAd = xAR[d];
    for (j=0; j<=d; j++) {
      v[k++] = xAd * xAR[j];
    }
    v[k-1] += E[d]; // nur korrekt falls E Vielfaches der EH-Matrix
    for ( ; j<dim; j++) {
      v[k++] = xAd * xAR[j];
    }
  }


}

void minmaxEigenEtAxxA(model *cov, double *mm) {
  double 
    *E = P(EAXXA_E);
  int i,
    dim = OWNLOGDIM(0);
  for (mm[0] = RF_INF, mm[1]=RF_NEGINF, i=0; i<dim; i++) {
    if (E[i] < mm[0]) mm[0] = E[i];
    if (E[i] > mm[1]) mm[1] = E[i];
  }
}
 
int checkEtAxxA(model *cov){
  int err;
  if (OWNXDIM(0) != 3) SERR("The space-time dimension must be 3.");
  VDIM0 = VDIM1 = OWNLOGDIM(0);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  cov->mpp.maxheights[0] = RF_NA;
 RETURN_NOERROR;
}
 
void rangeEtAxxA(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  int i;

  for (i=0; i<=2; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = - 1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }

  range->min[EAXXA_E] = 0.0;
  range->max[EAXXA_E] = RF_INF;
  range->pmin[EAXXA_E] = 0.0001;
  range->pmax[EAXXA_E] = 10;
  range->openmin[EAXXA_E] = true;
  range->openmax[EAXXA_E] = true;
}





// ----------------------------------------------------------------------
//
// Sigma(x) = diag>0 + A'xx'A
#define ROTAT_PHI 0 // both rotat and Rotat
#define ROTAT_SPEED 1
void kappa_rotat(int i, model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = i < DefList[COVNR].kappas ? 1 : -1;
}
void rotat(double *x, model *cov, double *v) {
  int
    dim = OWNLOGDIM(0),
    time = dim - 1;
  double
    speed = P0(ROTAT_SPEED),
    phi = P0(ROTAT_PHI),
    absx = SQRT(x[0] * x[0] + x[1] * x[1]);
  *v = (absx == 0.0) ? 0.0
    : speed * (COS(phi * x[time]) * x[0] + SIN(phi * x[time]) * x[1]) / absx;
}

void minmaxEigenrotat(model VARIABLE_IS_NOT_USED *cov, double *mm) {
  mm[0] = - 1;
  mm[1] = 1;
}
 
int checkrotat(model *cov){
  int err;
  if (OWNXDIM(0) != 3) SERR("The space-time dimension must be 3.");
   if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
 cov->mpp.maxheights[0] = RF_NA;
  RETURN_NOERROR;
}
 
void rangerotat(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  int i;

  for (i=0; i<2; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = - 1e10;
    range->pmax[i] = 1e10;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
}


// Sigma(x) = diag>0 + A'xx'A
#define ROTAT_PHI 0 // both rotat and Rotat
#define ROTAT_SPEED 1
void kappa_Rotat(int i, model *cov, int *nr, int *nc){
  *nc = 1;
  *nr = i < DefList[COVNR].kappas ?  1 : -1;
}
void Rotat(double *x, model *cov, double *v) {
  int d, k, j, 
    dim = OWNLOGDIM(0),
    time = dim - 1;
  double
    phi = P0(ROTAT_PHI),
      c =  COS(phi * x[time]),
      s = SIN(phi * x[time]),
      R[9]; assert(dim ==3);
   
  R[0] = R[4] = c;
  R[1] = s;
  R[3] = -s;
  R[2] = R[5] = R[6] = R[7] = 0.0;
  R[8] = 1.0;
 
  for (k=0, d=0; d<dim; d++) {
    v[d] =  0.0;
    for (j=0; j<dim; j++) {
      v[d] += x[j] * R[k++];
    }
  }
} 
int checkRotat(model *cov){
  int err;
  if (OWNXDIM(0) != 3) SERR("The space-time dimension must be 3.");
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  VDIM0 = VDIM1 = OWNLOGDIM(0);
  cov->mpp.maxheights[0] = RF_NA;
  RETURN_NOERROR;
}

void rangeRotat(model VARIABLE_IS_NOT_USED *cov, range_type* range){
 
  range->min[ROTAT_PHI] = RF_NEGINF;
  range->max[ROTAT_PHI] = RF_INF;
  range->pmin[ROTAT_PHI] = - 1e10;
  range->pmax[ROTAT_PHI] = 1e10;
  range->openmin[ROTAT_PHI] = true;
  range->openmax[ROTAT_PHI] = true;
}




// ----------------------------------------------------------------------
//

void truncsupport(double *x, model *cov, double *v){
  // truncsupport erwartet dass vorher $ kommt, sofern skalenmischung
  model *next = cov->sub[0];
  int xdimown = OWNTOTALXDIM;
  double dist,
    radius = P0(TRUNC_RADIUS); // default -1
  if (xdimown > 1) {
    int i;
    dist = 0.0;
    for (i=0; i<xdimown; i++) dist += x[i] * x[i];
    dist = SQRT(dist);
  } else dist = FABS(*x);

  if (radius>=0 && dist > radius) { *v=0.0; return; }
  FCTN(x, next, v);
}

int checktruncsupport(model *cov) {
  model *next=cov->sub[0];
  int err;
  //    dim = OWNLOGDIM(0); // taken[MAX DIM],
  ASSERT_UNREDUCED;
  ASSERT_ONESYSTEM;
  set_maxdim(OWN, 0, INFDIM);
  cov->monotone = isMonotone(next->monotone) ? MONOTONE : NOT_MONOTONE;

  if ((err = CHECK_PASSTF(next, ShapeType, SUBMODEL_DEP, cov->frame)) != NOERROR) {
    //       print("error !!\n");
    RETURN_ERR(err);
  }
  setbackward(cov, next);
  RETURN_NOERROR;
}


void truncsupportInverse(double VARIABLE_IS_NOT_USED *x,
			 model *cov, double *v){
  *v = P0(TRUNC_RADIUS);
}

void rangetruncsupport(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[TRUNC_RADIUS] = 0;      // < 0 : unset
  range->max[TRUNC_RADIUS] = RF_INF;
  range->pmin[TRUNC_RADIUS] = 0.00001;
  range->pmax[TRUNC_RADIUS] = 100.0;
  range->openmin[TRUNC_RADIUS] = true;
  range->openmax[TRUNC_RADIUS] = true;
}

int struct_truncsupport(model *cov, model **newmodel) {  
  int err;

  ASSERT_NEWMODEL_NOT_NULL;

  if (hasPoissonFrame(cov) || hasSmithFrame(cov)) {
    if ((err = addUnifModel(cov, P0(TRUNC_RADIUS), newmodel)) != NOERROR)
      RETURN_ERR(err);
  } else ILLEGAL_FRAME_STRUCT;

  RETURN_NOERROR;
}

int init_truncsupport(model *cov, gen_storage *s) {
  int i, err,
    vdim = VDIM0;
  assert(VDIM0 == VDIM1);

  if (hasSmithFrame(cov) || hasAnyPoissonFrame(cov)) {
    model *next = cov->sub[0];
    //    double
    //      radius = P0(TRUNC_RADIUS); // default -1
    
    if ((err = INIT(next, cov->mpp.moments, s)) != NOERROR) RETURN_ERR(err);
    //   if (radius>=0 && radius < cov->mpp.refradius) cov->mpp.refradius = radius;
    
    // Eplus, M2 are assumed to be still precise !!
    for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = next->mpp.maxheights[i];

    RETURN_NOERROR;
  }

  else ILLEGAL_FRAME;
}


void do_truncsupport(model *cov, gen_storage *s) {
  //  mppinfotype *info = &(s->mppinfo);
  model *next = cov->sub[0];
  int i, 
    vdim = VDIM0;
 //  double
  //    radius = P0(TRUNC_RADIUS); // default -1
  assert(VDIM0 == VDIM1);
  
  DO(next, s);  
  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = next->mpp.maxheights[i];

  // if (radius>=0 && radius < info->radius) info->radius = radius;
}


