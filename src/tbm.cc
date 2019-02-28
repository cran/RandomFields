/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by turning bands;
 see RFspectral.cc for spectral turning bands

 Copyright (C) 2001 -- 2017 Martin Schlather, 

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

#include <Rmath.h>  
#include <stdio.h>  
//#include <stdlib.h>
//#include <sys/timeb.h>
#include <unistd.h>
 
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>     

#include "questions.h"
#include "Processes.h"
#include "operator.h"

#define MAXNN 100000000.0 /* max number of points simulated on a line */

#define TBM_LINES (COMMON_GAUSS + 4)
#define TBM_LINESIMUFACTOR (COMMON_GAUSS + 5)
#define TBM_LINESIMUSTEP (COMMON_GAUSS + 6)
#define TBM_CENTER (COMMON_GAUSS + 7)
#define TBM_POINTS (COMMON_GAUSS + 8)


#define TBM_COV 0

#define BUFFER 4.0


void unitvector3D(int projectiondim, double *deltax, double *deltay, 
		double *deltaz) {
  switch (projectiondim) { // create random unit vector
  case 1 : 
    *deltax= 2.0 * UNIFORM_RANDOM - 1.0;
    *deltaz = *deltay = 0.0;
    break;
  case 2 :
    *deltaz = 0.0;
    *deltax= 2.0 *UNIFORM_RANDOM - 1.0;// see Martin's tech rep for details
    *deltay= SQRT(1.0 - *deltax * *deltax) * SIN(UNIFORM_RANDOM*TWOPI);
    break;
  case 3 : 
    double dummy;
    *deltaz = 2.0 * UNIFORM_RANDOM - 1.0; // ir case of multivariate
    dummy = SQRT(1.0 - *deltaz * *deltaz);
    *deltay = UNIFORM_RANDOM * TWOPI;
    *deltax = COS(*deltay) * dummy;
    *deltay = SIN(*deltay) * dummy;
    break;
  default : BUG;
  }
}


model* get_user_input(model *cov) {
  // loc_user sind nicht die prevloc data, sondern, die die der Nutzer
  // eingegeben hatte !! Zumindest in erster Naehrung  
  // Problem hier falls intern aufgerufen wurde -- Abhilfe: layer muss
  // zwingend gesetzt sein
  model *user=cov;
  while (user->calling != NULL) user = user->calling;
  return user;
}


int get_subdim(model *cov, bool Time, bool *ce_dim2, int *ce_dim,
	       int *effectivedim) {
  model *next = cov->sub[TBM_COV];
  int fulldim = P0INT(TBM_FULLDIM),
    layers = P0INT(TBM_LAYERS);
  *effectivedim = ANYDIM; 
  if (Time) {
    *ce_dim2 = (layers == (int) True) || equalsSpaceIsotropic(NEXT) ||
      *effectivedim == 1 + fulldim;
    *effectivedim -= *ce_dim2;  
    if (*ce_dim2 && layers == (int) False) {
      SERR1("value of '%.50s' does not match the situation", KNAME(TBM_LAYERS))
    }
  } else {
    *ce_dim2 = false;
  }
  if (*effectivedim > fulldim) RETURN_ERR(ERRORWRONGDIM);
  *ce_dim = 1 + (int) *ce_dim2;
  RETURN_NOERROR;
}


void tbm_kappasproc(int i, model *cov, int *nr, int *nc){
  kappaGProc(i, cov, nr, nc);
  int dim = ANYDIM; 
  //  printf("tbm kappas %d %d %d\n", i, TBM_CENTER, dim);
  if (i==TBM_CENTER) *nr = dim;
}


void rangetbmproc(model *cov, range_type *range){ 
  GAUSS_COMMON_RANGE;
  rangetbm_common(cov, range, false);

  range->min[TBM_LINES] = 1.0;
  range->max[TBM_LINES] = RF_INF;
  range->pmin[TBM_LINES] = 1.0;
  range->pmax[TBM_LINES] = 10000;
  range->openmin[TBM_LINES] = false;
  range->openmax[TBM_LINES] = true;

  range->min[TBM_LINESIMUFACTOR] = 0.0;
  range->max[TBM_LINESIMUFACTOR] = RF_INF;
  range->pmin[TBM_LINESIMUFACTOR] = 0.0;
  range->pmax[TBM_LINESIMUFACTOR] = 10.0;
  range->openmin[TBM_LINESIMUFACTOR] = false;
  range->openmax[TBM_LINESIMUFACTOR] = true;

  range->min[TBM_LINESIMUSTEP] = 0.0;
  range->max[TBM_LINESIMUSTEP] = RF_INF;
  range->pmin[TBM_LINESIMUSTEP] = 0.0;
  range->pmax[TBM_LINESIMUSTEP] = 10.0;
  range->openmin[TBM_LINESIMUSTEP] = false;
  range->openmax[TBM_LINESIMUSTEP] = true;

  // range->min[TBM_GRID] = 0;
  //range->max[TBM_GRID] = 1;
  // range->pmin[TBM_GRID] = 0;
  // range->pmax[TBM_GRID] = 1;
  // range->openmin[TBM_GRID] = false;
  //  range->openmax[TBM_GRID] = false;

  range->min[TBM_CENTER] = RF_NEGINF;
  range->max[TBM_CENTER] = RF_INF;
  range->pmin[TBM_CENTER] = - 1000000;
  range->pmax[TBM_CENTER] = 1000000;
  range->openmin[TBM_CENTER] = false;
  range->openmax[TBM_CENTER] = true;

  range->min[TBM_POINTS] = 0;
  range->max[TBM_POINTS] = RF_INF;
  range->pmin[TBM_POINTS] = 0;
  range->pmax[TBM_POINTS] = 1000;
  range->openmin[TBM_POINTS] = false;
  range->openmax[TBM_POINTS] = true;

}


int checktbmproc(model *cov) {

  FRAME_ASSERT_GAUSS_INTERFACE;

  // TO DO : split TBM in a RPmodel that expects
  // a process as submodel and performs there just 
  // the required lines numbers and resolution
  // And in a RPmodel that sets up from a covariance function C
  // a gaussian process for tbm operator(C)

#define NSEL 2
  model 
    *next=cov->sub[TBM_COV],
    *key =cov->key,
    *sub = key==NULL ? next : key;
  int  err = NOERROR;  // taken[MAX DIM],
  isotropy_type
    isoselect[NSEL]={ISOTROPIC, DOUBLEISOTROPIC};
  tbm_param *gp  = &(GLOBAL.tbm);
  KEY_type *KT = cov->base;

  ASSERT_CARTESIAN;
  ASSERT_ONESYSTEM;
  
  if (GLOBAL.general.vdim_close_together && VDIM0 > 1)
    SERR1("TBM only allowed if '%.50s=FALSE'", general[GENERAL_CLOSE]);

  // printf("%.50s fulldim %d %d\n", NICK(cov), fulldim, 0);
  //   printf("%.50s fulldim %d %d\n", NICK(cov),
  // 	 fulldim, gp->lines[fulldim-1]);assert(false);

  kdefault(cov, TBM_FULLDIM, gp->fulldim);
  kdefault(cov, TBM_FULLDIM, PisNULL(TBM_TBMDIM) || gp->tbmdim >= 0
	   ? gp->fulldim : P0INT(TBM_TBMDIM) - gp->tbmdim);
    kdefault(cov, TBM_TBMDIM, gp->tbmdim > 0 
	   ? gp->tbmdim : P0INT(TBM_FULLDIM) + gp->tbmdim);
    kdefault(cov, TBM_LAYERS, (int) gp->layers);
 int 
    tbmdim = P0INT(TBM_TBMDIM),
    fulldim = P0INT(TBM_FULLDIM);
  if (tbmdim >= fulldim) {
     SERR4("'%.50s' (=%d) must be less than '%.50s' (=%d)", 
	   KNAME(TBMOP_TBMDIM), tbmdim, KNAME(TBM_FULLDIM), fulldim);
  }
  kdefault(cov, TBM_LINES, gp->lines[fulldim-1]);
  kdefault(cov, TBM_LINESIMUFACTOR, gp->linesimufactor); 
  kdefault(cov, TBM_LINESIMUSTEP, gp->linesimustep);
  //  kdefault(cov, TBM_GRID, gp->grid);
  // if ( P0INT(TBM_GRID))
  //  warning("grid parameter for tbm not programmed yet");
  if (PisNULL(TBM_CENTER)) {    
    int Dim = OWNTOTALXDIM;
      PALLOC(TBM_CENTER, Dim, 1);
      for (int d=0; d<Dim; d++) {
      P(TBM_CENTER)[d] = gp->center[d];
    }
  } else {
    if (cov->nrow[TBM_CENTER] < OWNTOTALXDIM) 
      SERR1("vector for '%.50s' too short", KNAME(TBM_CENTER));
  }

  kdefault(cov, TBM_POINTS, gp->points);
  if ((err = checkkappas(cov, false)) != NOERROR)  RETURN_ERR(err);
 
  if (key == NULL && !isProcess(sub)) { // thus isVariogram(sub)
    // Abfolge Tbm $(Aniso) iso-model braucht Moeglichkeit des 
    // anisotropen Modells
    if (hasInterfaceFrame(cov)) { // || COVNR == TBM_PROC_USER) nsel++;

      // PMI(cov);
      //      if (COVNR != TBM_PROC_USER) crash();
      
      assert(COVNR == TBM_PROC_USER);
      err = CHECK(sub, OWNLOGDIM(0), OWNXDIM(0), VariogramType, 
		       KERNEL, // wegen nutzer eingabe aniso und dem
		       // allerersten check
		       SYMMETRIC, SUBMODEL_DEP, GaussMethodType);
    } else {
      assert(COVNR != TBM_PROC_USER);
      for (int i=0; i<NSEL; i++) {	
	if ((err = CHECK(sub, OWNLOGDIM(0), OWNXDIM(0), VariogramType, XONLY, 
			 isoselect[i], SUBMODEL_DEP, GaussMethodType))
	    == NOERROR) break;
      }
    }

    if (err != NOERROR) {
       SPRINTF(KT->error_loc, "%.50s", NICK(cov));
      SERR("Its submodel could not be identified as isotropic or space-isotropic and positive definite function.");
    }
    
    
  } else { // SUBNR > FIRSTGAUSSPROC || key != NULL

    bool dummy;
    int dummydim, newdim;

    if (key != NULL) {
      // falls hier gelandet, so ruft TBMINTERN TBM auf!!
      if (COVNR == TBM_PROC_USER) { 
	model *intern = sub;
	while (intern != NULL && //MODELNR(intern) != TBM_PROC_INTERN &&
	       (isAnyDollar(intern) || MODELNR(intern) == TBM_PROC_USER  // 
		))
	  intern = intern->key != NULL ? intern->key : intern->sub[0];
	if (intern == NULL) {
	  BUG;
	} else if (intern != cov) 
	  paramcpy(intern, cov, true, true, false, false, false);
      }
    }

    if (COVNR == TBM_PROC_USER) newdim = OWNXDIM(0); 
    else if ((err = get_subdim(cov, PrevLoc(get_user_input(cov))->Time, &dummy,
			       &newdim, &dummydim)) != NOERROR) RETURN_ERR(err);

    if ((err = CHECK(sub, newdim, newdim, ProcessType, XONLY, CARTESIAN_COORD,
		     SUBMODEL_DEP, 
		     GaussMethodType))
	!= NOERROR) {
      RETURN_ERR(err);
    }
  }
  setbackward(cov, sub);
  VDIM0 = next->vdim[0];
  VDIM1 = next->vdim[1];
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov)) != NOERROR)  RETURN_ERR(err);

  RETURN_NOERROR;
}




int struct_tbmproc(model *cov, model **newmodel) { 
  model
    *next = cov->sub[TBM_COV];
  ASSERT_NEWMODEL_NULL;
  if (next->pref[TBM] == PREF_NONE) {
    RETURN_ERR(ERRORPREFNONE);
  }

  assert(COVNR == TBM_PROC_INTERN);

  model *user = get_user_input(cov);
  location_type *loc_user = PrevLoc(user);

  double linesimufactor, 
    tbm_linesimustep = P0(TBM_LINESIMUSTEP),
    tbm_linesimufactor = P0(TBM_LINESIMUFACTOR);
  int d, 
    err=NOERROR,
    //  tbm_lines = P0INT(TBM_LINES),
    fulldim = P0INT(TBM_FULLDIM), 
    tbmdim = P0INT(TBM_TBMDIM),
    //  endaniso = origdim * origdim - 1,
    *points = PINT(TBM_POINTS),
    user_dim = loc_user->timespacedim;

  long
    user_spatialpoints = loc_user->spatialtotalpoints,
    user_spatialdim = loc_user->spatialdim
    // the isotropic part (i.e. time might be not considered)
    ;
  
  bool ce_dim2=false; 


  if (!hasGaussMethodFrame(cov) && !hasInterfaceFrame(cov)) {
    RETURN_ERR(ERRORFAILED);
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
   submatrix for v=0, already checked by Check
  */
 
  if (tbmdim != 1 || (fulldim != 2 && fulldim != 3))
    SERR5("'%.50s' only works for %.50s =1 and %.50s=2 or 3. Got %d %d",
	  NICK(cov), KNAME(TBM_TBMDIM), KNAME(TBM_FULLDIM), tbmdim, fulldim);
  if (ANYDIM > MAXTBMSPDIM) RETURN_ERR(ERRORMAXDIMMETH);


  // in theory it works also for variograms!
 
  if(!isVariogram(NEXTTYPE(0)) || !isXonly(NEXT)) RETURN_ERR(ERRORNOVARIOGRAM);

  if (user_dim == 1) 
    SERR("dimension must currently be at least 2 and at most 4 for TBM");

  ONCE_NEW_STORAGE(tbm);
  tbm_storage *s = cov->Stbm;


  if ((err=get_subdim(cov, loc_user->Time, &ce_dim2, &(s->ce_dim),
		      &(s->simuspatialdim)))
      != NOERROR) RETURN_ERR(err);


   /****************************************************************/
  /*          determine grid distance for the line                */
  /****************************************************************/

 // sort out which finer grid should be used in the spatial domain for
  // low dimensional simulation
  if (PL>=PL_STRUCTURE) {PRINTF("\ndetermining the length of line segment...");}

  // PrintMethodInfo(meth);

  // minimum distance between the points   
  linesimufactor = - 1.0;
  if (tbm_linesimustep  > 0.0) linesimufactor = 1.0 / tbm_linesimustep;
  else if (tbm_linesimufactor > 0) {
    double mindelta = RF_INF;
    int spDim = user_dim - (int) ce_dim2; // eigentlich muesste man
    // es das Urbild der Zeitachse abziehen... Egal. Fuer den Ausbau reicht es.
    if (loc_user->grid) {
      // if linesimufactor, the fine grid scale is proportional to smallest 
      // distance in the grid. 
      for (d=0; d<spDim; d++) {
        if ((loc_user->xgr[d][XSTEP]<mindelta) && (loc_user->xgr[d][XSTEP]>0)) {
	  mindelta=loc_user->xgr[d][XSTEP];
	}
      }
    } else {
      long j, i, k0, k1, k2;	
      if (user_spatialpoints > 10000) {
	SERR1("too many points to determine smallest distance between the points in a reasonable time; try '%.50s' with a positive value", KNAME(TBM_LINESIMUSTEP));
	  /* algorithmus kann verbessert werden, so dass diese Fehlermeldung 
	     nicht mehr notwendig ist! */ 
      }
      double diff, dist;
      for (k0=i=0; i<user_spatialpoints; i++, k0+=user_spatialdim) {
	for (k2=j=0; j<i; j++) {
	  k1 = k0;
	  for (dist=0.0, d=0; d<user_spatialdim; d++) {
	    diff = loc_user->x[k1++] - loc_user->x[k2++]; 
	    dist += diff * diff;
	    }
	  if (dist>0 && dist<mindelta) mindelta=dist; 
	}
      } // if linesimustep==0.0
      if (user_dim == spDim && loc_user->Time) {
	mindelta += loc_user->T[XSTEP] * loc_user->T[XSTEP];
      }
      mindelta = SQRT(mindelta);
    }
    linesimufactor = tbm_linesimufactor / mindelta;
  } else {
    // both, tbm_linesimustep and tbm_linesimufactor are zeror or negative
    if (*points < BUFFER) { 
      SERR("if linesimufactor and linesimustep are naught then TBM.points must be at least 4 (better between 100 and 10000)");
    }
  }


  // multiplication der x-skalenparameter so dass letztendlich
  // auf der TBM-Geraden mit Abstand 1 der Punkte simuliert werden kann
  // mit dem Vorteil des einfachen Zugriffs auf die simulierten Werte 

  
  /****************************************************************/
  /*          determine length and resolution of the band         */
  /****************************************************************/

  location_type *loc = Loc(cov);
  assert(cov->ownloc == NULL);
  if ((loc->Time && !ce_dim2)) {
    TransformLoc(cov, ce_dim2, False, false); 
    loc = Loc(cov);
  }

  int 
    dim = GetLoctsdim(cov),
    spdim = dim - ce_dim2,
    ens = GLOBAL.general.expected_number_simu;


  double diameter, min[MAXTBMSPDIM], max[MAXTBMSPDIM], Center[MAXTBMSPDIM],
    *timecomp = loc->grid ? loc->xgr[dim - 1] : loc->T,
    min_center = 0.0,
    max_center = 0.0;
  
  GetDiameter(loc, min, max, Center, true, false, NULL);
  for (d=0; d<spdim; d++) { 
    // necessary to calculate for min and max in case of tbm_center is given
    double dummy;
    // printf("%d %10g %10g\n", d, Center[d], min[d]); 
    dummy = Center[d] - min[d];
    min_center += dummy * dummy;
    dummy = Center[d] - max[d];
    max_center += dummy * dummy;    
  }
  diameter = 2.0 * SQRT(min_center > max_center ? min_center : max_center);

  //    print("diam %10g pts=%d %10g %10g\n", diameter, *points, BUFFER, linesimufactor ); //assert(false);
  // 4.231878 0 3.000000 1054.407676


  if (linesimufactor < 0) { // i.e. points given directly by user
      linesimufactor = (*points - BUFFER) / diameter;
  }
  s->linesimufactor = linesimufactor;

  diameter = TRUNC(BUFFER + diameter * linesimufactor); // in pts

  // print("diam %10g %10g\n", diameter, linesimufactor);

  if ( *points > 0) { 
    if (diameter > *points) { // never happens if linesimufactor has been < 0
      PRINTF("tbm points minimum=%10g, but tbm_points is %d.\n", 
	     diameter,  *points);
      SERR( "given number of points to simulate on the TBM line is too small");
    } else diameter = (double) *points;
  } 

// print("diam %10g %10g %d %d\n", diameter, linesimufactor, *points, GLOBAL.tbm.points);

  if (diameter > MAXNN) {
    SERR("The number of points on the line is too large. Check your parameters and make sure that none of the locations are given twice");
  } else if (diameter < 10) {
    SERR2("too few points on the line -- modify '%.50s' or '%.50s'",
	  KNAME(TBM_LINESIMUFACTOR), KNAME(TBM_LINESIMUSTEP));
  }
 

  // only needed for nongrid !
  s->spatialtotalpts = loc->totalpoints;
  if (s->ce_dim == 2) {
    s->spatialtotalpts /= timecomp[XLENGTH];
    //    assert(false);
  }
 

  /****************************************************************/
  /*          determine the center                                 */
  /****************************************************************/

  double 
     *tbm_center = P(TBM_CENTER);
  bool setable = !isDollar(cov->calling) && user_dim == dim
    && user_spatialdim == GetLoctsdim(cov);

  GetDiameter(loc, min, max, Center, false, false, NULL);

  if (!ISNAN(tbm_center[0])) {
    assert(cov->calling != NULL);
    if (!setable) 
      SERR1("'%.50s' can be set only for isotropic fields", KNAME(TBM_CENTER));
    for (d=0; d < user_dim; d++) {
      if (!R_FINITE(tbm_center[d])) SERR("center contains non finite values");
    }
    for (d=0; d<user_dim; d++) { // not here no time component!
      Center[d] = tbm_center[d]; // in user_ coordinates !
    }    
  } else {
    if (setable) for (d=0; d<user_dim; d++) tbm_center[d] = Center[d];   
  }



  for (d=0; d<dim; d++) {
    s->center[d] = Center[d];
    if (loc->grid) s->center[d] -= loc->xgr[d][XSTART]; 
    else if (loc->Time && d==dim-1 && ce_dim2) {
      s->center[d] -= loc->T[XSTART];
    }
  }



  /****************************************************************/
  /*          set up the covariance function on the line          */
  /****************************************************************/

   if (PL>=PL_STRUCTURE) {
    PRINTF("\nrestructing tbm ");
    PRINTF(ce_dim2 ? "plane...\n" : "line...\n");
  }

  model  *modelB1;

  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);
  if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
  modelB1 = cov->key;
  while (isGaussMethod(modelB1)) modelB1 = modelB1->sub[0];
  if (modelB1 == cov->key) {
    addModel(&(cov->key), GAUSSPROC);
    modelB1 = cov->key;
  } else {
    modelB1 = modelB1->calling;
  }
 
  addModel(modelB1, 0, DOLLAR);
  if (MODELNR(modelB1)==BROWNIAN){
    kdefault(cov, DVAR, 1.0 + PARAM0(next, BROWN_ALPHA) / tbmdim);
  } else {
    model *dollar = modelB1->sub[0];
    set_both_systems(dollar, s->ce_dim, PosDefType );
    kdefault(dollar, DVAR, 1.0);
    PARAMALLOC(dollar, DANISO, s->ce_dim, s->ce_dim);
    PARAM(dollar, DANISO)[0] = 1.0 / linesimufactor;

    if (s->ce_dim == 2) {
      PARAM(dollar, DANISO)[1] = PARAM(dollar, DANISO)[2] = 0.0;
      PARAM(dollar, DANISO)[3] = loc->T[XSTEP];
    } 
    addModel(modelB1, 0, TBM_OP);
    if (!PisNULL(TBMOP_TBMDIM))
      kdefault(modelB1->sub[0], TBMOP_TBMDIM, P0INT(TBM_TBMDIM));
    if (!PisNULL(TBMOP_FULLDIM))
      kdefault(modelB1->sub[0], TBMOP_FULLDIM, P0INT(TBM_FULLDIM));
    if (!PisNULL(TBMOP_LAYERS))
      kdefault(modelB1->sub[0], TBMOP_LAYERS, P0INT(TBM_LAYERS));
  }
 

  /****************************************************************/
  /*                set up the process on the line                */
  /****************************************************************/

  if (PL>=PL_STRUCTURE) {
    //   printf("PL=%d %d\n", PL, PL_STRUCTURE);
    PRINTF("\nsetting up the process on the ");
    PRINTF(ce_dim2 ? "plane..." : "line...");
  }
  
  // xline[XEND]: number of points on the TBM line
  double xline[3];
  xline[XSTART] = 1.0;
  xline[XSTEP] = 1.0;
  s->xline_length  = xline[XLENGTH] = (double) diameter;// diameter > MAXNN must be first since
  //                                risk of overflow !

  // printf("length %10g %10g\n", diameter, 1.0 / linesimufactor); assert(false);
  // 4465  0.000948


  //printf("xline %10g %10g %10g\n", xline[XSTART], xline[XSTEP], xline[XLENGTH]);
  
  err = loc_set(xline, timecomp, 1, 1, 3, ce_dim2 /* time */, 
	  true /* grid */, false, cov->key);
  model *key = cov->key;
  if (err != NOERROR) goto ErrorHandling;
 
  GLOBAL.general.expected_number_simu = 100;
  if ((err = CHECK(key, 1 + (int) ce_dim2, 1 + (int) ce_dim2, 
		   ProcessType, XONLY, CARTESIAN_COORD,
		   cov->vdim, GaussMethodType)) != NOERROR) {
    goto ErrorHandling;  
  }

  
  assert(XDIM(SYSOF(key), 0) == s->ce_dim);

  if ((err = STRUCT(key, NULL)) != NOERROR) goto ErrorHandling;  
  if (!key->initialised) {
    if (key->key != NULL &&
	(err = CHECK(key, 1 + (int) ce_dim2, 1 + (int) ce_dim2, 
		     PosDefType, XONLY, SYMMETRIC,
		     SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
      goto ErrorHandling;  
    }
    if (PL >= PL_DETAILS) {
      PRINTF("\n\nCheckModel Internal (%.50s, #=%d), after 2nd check:", 
	     NICK(key), MODELNR(key)); // ok
    }
  }

  NEW_STORAGE(gen);

  
  if ((err = INIT(key, 0, cov->Sgen)) != NOERROR) goto ErrorHandling;

  s->err = err;

 ErrorHandling :
  //cov->initialised = err==NOERROR && key->initialised;
  GLOBAL.general.expected_number_simu = ens;
  RETURN_ERR(err);
}


// aufpassen auf ( * * 0 )
//               ( * * 0 )
//               ( 0 0 1 )


int init_tbmproc(model *cov, gen_storage *S) {
  location_type *loc = Loc(cov);
  int err=NOERROR;
  errorloc_type errorloc_save;
  model *key = cov->key; // == NULL ? cov->sub[0] : cov->key;
  assert(key != NULL);  
  tbm_storage *s = cov->Stbm;
  KEY_type *KT = cov->base;
   
  assert(s != NULL);
  assert(KT != NULL);
  STRCPY(errorloc_save, KT->error_loc);
  SPRINTF(KT->error_loc, "%.500s %.50s", errorloc_save, NAME(cov));
  cov->method = TBM;

  FRAME_ASSERT_GAUSS;

  if (s->err==NOERROR) {
    err = INIT(key, 0, S); 
  } else {
    assert(s->err < NOERROR);
  }

  STRCPY(KT->error_loc, errorloc_save);
  if (err != NOERROR) RETURN_ERR(err); 

  if (loc->distances) RETURN_ERR(ERRORFAILED);
   
  err = ReturnOwnField(cov);
  cov->simu.active = err == NOERROR;
 
  if (PL>= PL_STRUCTURE) {PRINTF("\n'%.50s' is now initialized.\n", NAME(cov));}

  RETURN_ERR(err);
  
}


void GetE(int fulldim, tbm_storage *s, int origdim, bool Time, 
	  double *phi, double deltaphi,	 double *aniso,
	  double *offset, double *ex, double *ey, double *ez, double *et) {
  double sube[MAXTBMSPDIM], e[MAXTBMSPDIM];
  int d, j, idx, k,
    spatialdim = s->simuspatialdim;

//      total = spatialdim * dim;
  // dim is the users's dim of the field
  // this might be reduced by zonal anisotropies leading to spatialdim

  //  printf("gete fulldim=%d %d %d\n", fulldim, origdim, spatialdim);

  // for debugging only
  for(d=0; d<MAXTBMSPDIM; d++) sube[d] = e[d] = RF_NEGINF;
  assert(fulldim >= spatialdim); 

  if (fulldim == 2) {
    if (deltaphi!=0.0) (*phi) += deltaphi;
    else (*phi) = UNIFORM_RANDOM * M_2_PI; // not pi, when multivariate
    sube[0] = SIN(*phi); 
    sube[1] = COS(*phi);
  }
  
  else if (fulldim == 3) {
    unitvector3D(spatialdim, sube, sube +1, sube + 2);
//    print("spatdim=%d %d  sube=%10g  %10g %10g\n", 
//	   spatialdim, s->simuspatialdim, sube[0], sube[1], sube[2]);
  }

  else RFERROR("wrong full dimension in 'GetE'");

  *offset = 0.5 * s->xline_length;

  //  print("\noffset %10g %10g aniso=%d %10g %10g fact=%10g\n", *offset,  s->xline_length, aniso==NULL, s->center[0], s->center[1], s->linesimufactor);
  

  if (aniso == NULL) {
    for (d=0; d<spatialdim; d++) e[d] = sube[d];    
  } else {
    for (d=0; d<spatialdim; d++) e[d] = 0.0;
    for (k=j=0; j<spatialdim; j++) {
      for (d=0; d<origdim; d++) {
	e[d] += sube[j] * aniso[k++];
	//      print("e,s,a, (%d, %d) %10g     %10g     %10g\n",
	//	     d, j, e[d], sube[j], aniso[j + d]);
      }
    }
  }
  for (d=0; d<spatialdim; d++) {
    e[d] *= s->linesimufactor;
    *offset -= s->center[d] * e[d];    
    //       printf("d=%d %10g sube=%10g %10g %10g %10g\n", d, e[d], sube[d],s->linesimufactor, s->center[d], *offset);
    //  
  }

  //   for (d=0; d<dim; d++) {
  //     print("e, lines %d offset=%10g e=%10g center=%10g %10g %d\n", 
  // 	     d, *offset, e[d], 
  // 	     s->center[d],  s->linesimufactor, s->ce_dim);
  //  }

  idx = spatialdim;
  if (Time && s->ce_dim==1) {
    *et = e[--idx];
  }
  switch (idx) {
  case 4 : BUG;
  case 3 : *ez = e[--idx]; FALLTHROUGH_OK; 
  case 2 : *ey = e[--idx]; FALLTHROUGH_OK; 
  case 1 : *ex = e[--idx];
    break;
  default: BUG;
  }
}



void do_tbmproc(model *cov, gen_storage  VARIABLE_IS_NOT_USED *S) {
  model 
    *key = cov->key; 
  assert(key != NULL); // == NULL ? cov->sub[0] : cov->key;
  location_type 
    *loc = Loc(cov),
    *keyloc = Loc(key);
  tbm_storage *s =  cov->Stbm;

 int vdim = VDIM0;
 double
   nn = keyloc->spatialtotalpoints;
 long
   ntot = keyloc->totalpoints * vdim,
   totvdim = loc->totalpoints * vdim;

  int
    origdim = loc->caniso == NULL ? ANYDIM : loc->cani_nrow,
    simuspatialdim = s->simuspatialdim, // for non-grids only
    every = GLOBAL.general.every,
    tbm_lines = P0INT(TBM_LINES);
   
  double
    *res = cov->rf,
    *simuline = key->rf;
  
  assert(res != NULL);
  assert(simuline != NULL);



  double phi, deltaphi, stepx, stepy, stepz, stept, 
    incx=0.0,
    incy=0.0, 
    incz=0.0,
    inct=0.0;
  long n;
  int v, nt, idx, gridlent,
    fulldim = P0INT(TBM_FULLDIM);
  SAVE_GAUSS_TRAFO;

  assert(cov->Sgen != NULL);
     
  for (n=0; n<totvdim; n++) {
    res[n]=0.0; 
  }
  deltaphi = PI / (double) tbm_lines;// only used for tbm2
  phi = deltaphi * UNIFORM_RANDOM;     // only used for tbm2
  if (!GLOBAL.tbm.grid) deltaphi = 0.0;

  if (loc->grid) {
    int nx, ny, nz, gridlenx, gridleny, gridlenz;
 
    stepx = stepy = stepz = stept =  0.0;
    gridlenx = gridleny = gridlenz = gridlent = 1;
  
    idx = origdim;
    if (s->ce_dim==2) {
      gridlent = (int) loc->xgr[--idx][XLENGTH]; // could be one !!
      stept = loc->xgr[idx][XSTEP];	    
      inct = nn; // wegen assert unten
    } else if (loc->Time) {
      idx--;
      gridlent = (int) loc->T[XLENGTH]; // could be one !!
      stept = loc->T[XSTEP];	    
    }

    switch (idx) {
	case 4 : 
	  BUG;
	case 3 : 
	  gridlenz = (int) loc->xgr[--idx][XLENGTH];
	  stepz = loc->xgr[idx][XSTEP];	  
	  FALLTHROUGH_OK; 
	case 2 : 
	  gridleny = (int) loc->xgr[--idx][XLENGTH];
	  stepy = loc->xgr[idx][XSTEP];	 
	  FALLTHROUGH_OK; 
	case 1 : 
	  gridlenx = (int) loc->xgr[--idx][XLENGTH];
	  stepx = loc->xgr[idx][XSTEP];	  
	  break;
    default: BUG;
    }
         
    for (n=0; n<tbm_lines; n++) {
      double inittoffset;
      R_CheckUserInterrupt();
      if (every>0  && n > 0 && (n % every == 0)) { PRINTF("%ld \n", n); }

      GetE(fulldim, s, origdim, loc->Time, &phi, deltaphi,
	   loc->caniso, &inittoffset, &incx, &incy, &incz, &inct);
      incx *= stepx;
      incy *= stepy;
      incz *= stepz;
      if (s->ce_dim == 1) inct *= stept;
      inittoffset += UNIFORM_RANDOM - 0.5;

      DO(key, cov->Sgen);

      // for (int i=0; i<10; i++) printf("%10e ", simuline[i]); printf("\n");  BUG;
      
    
 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) collapse(4)
#endif
     for (v=0; v<vdim; v++) {
	for (nt=0; nt<gridlent; nt++) {
	  for (nz=0; nz<gridlenz; nz++) {	  
	    for (ny=0; ny<gridleny; ny++) {	  
	      double xoffset = inittoffset + s->xline_length * v +
		inct * nt + incz * nz + incy * ny,
		endxoffset = xoffset + (gridlenx - 1) * incx;
	      long zaehler = gridlenx *
		(ny + gridleny * (nz + gridlenz * (nt + v * gridlent)));
	      if (xoffset >= ntot+1 || xoffset < 0 ||
		  endxoffset >= ntot+1 || endxoffset < 0 )  BUG;
	      for (nx=0; nx<gridlenx; nx++) {
		// long longxoffset = (long) xoffset;
		// if (longxoffset >= ntot || longxoffset < 0)  BUG;
		res[zaehler++] += simuline[(long) xoffset];
		xoffset += incx;
	      }
	      
	    }
	  }
	}
      } // vdim
    } // n

  } else { // not loc->grid
    //     assert(false);
    // not grid, could be time-model!
    // both old and new form included
    int spatialdim = loc->spatialdim;
    long totpoints = s->spatialtotalpts,
      end = spatialdim * totpoints;

#define TBMSTART {							\
      for (n=0; n<tbm_lines; n++){/*printf("n=%ld tbm.lines%d\n",n,tbm_lines);//*/ \
	R_CheckUserInterrupt();						\
        if (every>0  && (n % every == 0)) { PRINTF("%ld \n",n); }	\
	double offsetinit;						\
	GetE(fulldim, s, origdim, loc->Time, &phi, deltaphi,		\
	     loc->caniso, &offsetinit, &incx, &incy, &incz, &inct);	\
        offsetinit += UNIFORM_RANDOM - 0.5;				\
	DO(key, cov->Sgen)

#define TBMEND }} break
    
#define TBMFOR								\
   	for (v=0; v<vdim; v++) {					\
	  for (nt=0; nt<gridlent; nt++) {				\
	    double offset = offsetinit + s->xline_length * v + inct * nt; \
	    long zaehler = spatialdim * (nt + gridlent * v);		\
	    for (int xi = 0; xi < end; xi+=spatialdim) {

#define TBMFOREND }}  

#define TBMFOR0							\
    for (int zaehler = 0; zaehler < totpoints; zaehler++) {	\
      int xi = zaehler * spatialdim;
   

#define TBMST(INDEX,ZAEHLER)						\
	      long index;						\
	      index = (long) (offset + INDEX);				\
	      if (index >= ntot || index < 0) {				\
		PRINTF("\n %10g %10g %10g (%10g %10g %10g))\n",			\
		       loc->x[xi], loc->x[xi+1] , loc->x[xi+2],		\
		       incx, incy, incz);				\
		PRINTF("n=%ld index=%ld nn=%10g ntot=%ld xi=%d \n", \
		       n, index, nn, ntot, xi);				\
		PRINTF("OFF=%10g IDX=%10g inct=%10glenT=%d spatialdim=%d\n",	\
		       offset, INDEX, inct, gridlent, spatialdim);	\
	       	BUG; /* to do: kann irgendwann mal rausgenommen werden, wenn's stabil laeuft */	\
	      }								\
 	      res[ZAEHLER] += simuline[index];				\
          }

    inct = nn; // number of spatial points
    gridlent = loc->Time ? loc->T[XLENGTH] : 1; 
    assert(gridlent * nn * vdim == ntot);

    switch (simuspatialdim) {
    case 1 : //see Martin's techrep f. details
      TBMSTART;
#ifdef DO_PARALLEL
      if (vdim * gridlent == 1) {
	double offset = offsetinit;
#pragma omp parallel for num_threads(CORES) 
	TBMFOR0;
 	TBMST(loc->x[xi] * incx, zaehler);
      } else {
#pragma omp parallel for num_threads(CORES) collapse(2)
#endif
	TBMFOR;
	TBMST(loc->x[xi] * incx, zaehler++);
        TBMFOREND
#ifdef DO_PARALLEL
	} 
#endif
      TBMEND;
     
    case 2: 
      TBMSTART;
#ifdef DO_PARALLEL
      if (vdim * gridlent == 1) {
	double offset = offsetinit;
#pragma omp parallel for num_threads(CORES) 
	TBMFOR0;
 	TBMST(loc->x[xi] * incx + loc->x[xi+1] * incy, zaehler);
      } else {
#pragma omp parallel for num_threads(CORES) collapse(2)
#endif
	TBMFOR;
	TBMST(loc->x[xi] * incx + loc->x[xi+1] * incy, zaehler++);
	TBMFOREND
#ifdef DO_PARALLEL
      }
#endif
      TBMEND;
      
    case 3:
#define INC3 loc->x[xi] * incx + loc->x[xi+1] * incy + loc->x[xi+2] * incz
      TBMSTART;
#ifdef DO_PARALLEL
      if (vdim * gridlent == 1) {
	double offset = offsetinit;
#pragma omp parallel for num_threads(CORES) 
	TBMFOR0;
 	TBMST(INC3, zaehler);	
      } else {
#pragma omp parallel for num_threads(CORES) collapse(2)
#endif
	TBMFOR;
	TBMST(INC3, zaehler++);
	TBMFOREND
#ifdef DO_PARALLEL
      }
#endif
      TBMEND;
     
    default : BUG;
    }
  } // end not grid
  
//  {double z; int i; for(i=0, z=0.0; i<ntot; i++) z+=simuline[i]*simuline[i];
//      print("%10g %10g %10g\n", 
//	     simuline[0], z / ntot, res[10] / SQRT((double) tbm_lines));
//  }

  long i;
  double InvSqrtNlines;   
  InvSqrtNlines = 1.0 / SQRT((double) tbm_lines);

 
  for(i=0; i<totvdim; res[i++] *= (double) InvSqrtNlines);
  BOXCOX_INVERSE;
}

