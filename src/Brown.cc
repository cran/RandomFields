/* 
 Authors
 Marco Oesting,
 Martin Schlather, schlather@math.uni-mannheim.de

 simulation of Brown-Resnick processes

 Copyright (C) 2009 -- 2010 Martin Schlather 
 Copyright (C) 2011 -- 2017 Marco Oesting & Martin Schlather 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

// #include <inttypes.h>


#include <stdio.h>
#include "def.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "Processes.h"
#include "operator.h"
#include "rf_interfaces.h"
#include "variogramAndCo.h"
#include "extern.h"

#define BR_MESHSIZE (LAST_MAXSTABLE + 1)
#define BR_VERTNUMBER (LAST_MAXSTABLE + 2)
#define BR_OPTIM (LAST_MAXSTABLE + 3)
#define BR_OPTIMTOL (LAST_MAXSTABLE + 4)
#define BR_LAMBDA (LAST_MAXSTABLE + 5)
#define BR_OPTIMAREA (LAST_MAXSTABLE + 6)
#define BR_VARIOBOUND (LAST_MAXSTABLE + 7)





int newmodel_covcpy(model **localcov, int modelnr, model *cov, //err
		    double *x, double *y, double *T, 
	            int spatialdim, /* spatial dim only ! */
	            int xdimOZ, long lx, long ly, bool Time, bool grid,
	            bool distances) {
  Types type = SYSTYPE(DefList[modelnr].systems[0], 0);
 assert(type == InterfaceType && DefList[modelnr].variants == 1); // anderes z.zt nicht benutzt

  int i, err;
  assert(*localcov == NULL);
  assert(cov->calling != NULL);
  addModel(localcov, modelnr, NULL, true);
  model *neu = *localcov;
  neu->base = cov->base;
  neu->root = neu;

  //  PMI(neu);
  
  assert(neu->ownloc==NULL && neu->prevloc==NULL);
  assert(type == InterfaceType); // braucht prevloc !! 
  neu->prevloc = LOCLIST_CREATE(1, xdimOZ + (int) Time); // locown
  assert(cov->prevloc != NULL);
  assert(PLoc(cov) == cov->prevloc);
  assert(!((y==NULL) xor (ly == 0)));
 
  loc_set(x, y, T, spatialdim, xdimOZ, lx, ly, Time, grid, distances, neu);
  if ((err = covcpy(neu->sub + 0, cov)) != NOERROR) RETURN_ERR(err);

  SET_CALLING(neu->sub[0], neu);

  for (i=0; i<2; i++) {
    if ((err = CHECK(neu, OWNLOGDIM(0), PREVXDIM(0),
		     // DefList[COVNR].t ype,Martin:changed 27.11.13 
		     //co v->ty pus,// auch nicht 20.5.14
		     type,
		     type == InterfaceType ? XONLY : PREVDOM(0), 
		     type == InterfaceType ? CARTESIAN_COORD : PREVISO(0), 
		     cov->vdim, EvaluationType)) != NOERROR) {
      RETURN_ERR(err);
    }
    if (i==0 && (err =  STRUCT(neu, NULL)) != NOERROR) RETURN_ERR(err);
  }
  RETURN_NOERROR;
}

int newmodel_covcpy(model **localcov, int modelnr, model *cov) {//err
  
  int err,
    store = GLOBAL.general.set;
  GLOBAL.general.set = 0;
  location_type *loc = Loc(cov);
  
  err = newmodel_covcpy(localcov, modelnr, cov,
               loc->grid ? loc->xgr[0] : loc->x, 
	       loc->grid ? loc->ygr[0] : loc->y, 
               loc->grid ? loc->xgr[0] + 3 * loc->spatialdim : loc->T, 
	       loc->spatialdim, loc->xdimOZ, 
	       loc->grid ? 3 : loc->totalpoints, 
	       loc->ly == 0 ? 0 : loc->grid ? 3 : loc->totalpoints,
	       loc->Time, loc->grid, loc->distances);
  
  GLOBAL.general.set = store;
  RETURN_ERR(err);
  
}


// **********************************************************************
// Brown Resnick

int checkBrownResnickProc(model *cov) {

  //NotProgrammedYet("at the moment Brown-Resnick processes");

  model  
    *key = cov->key,
    *sub = key != NULL ? key :
            cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  int err;
  isotropy_type isoprev;
  Types type, frame;
  
  ASSERT_ONESYSTEM;
  ASSERT_CARTESIAN;
  ASSERT_ONE_SUBMODEL(cov);
  
  if ((err = SetGEVetc(cov)) != NOERROR) RETURN_ERR(err);
    
  type = isProcess(sub) || isPointShape(sub)
    ? SYSTYPE(DEFSYS(sub), 0) : VariogramType;
  frame = isVariogram(type) ? EvaluationType : BrMethodType;     
  isoprev = equalsVariogram(frame) ? SYMMETRIC : CARTESIAN_COORD;
  
  if ((err = CHECK(sub, OWNLOGDIM(0), OWNXDIM(0),
		   type,//Martin: changed 27.11.13 DefList[COVNR].Type,
		   XONLY, isoprev, SCALAR, frame)) != NOERROR) {
      RETURN_ERR(err);
  }

 
  setbackward(cov, sub);
  
  RETURN_NOERROR;  
}

int general_init(model* cov, int trendlen, gen_storage *s) {
  model *key=cov->key;
  br_storage *sBR = cov->Sbr;
  int err,
    total = Gettotalpoints(key),
    dim = ANYDIM; // ownxdim(0)
  bool keygrid = Getgrid(key);
  
  key->simu.expected_number_simu = cov->simu.expected_number_simu;
  if ((err = INIT(key, COVNR != BRNORMED, s)) != NOERROR) {
    RETURN_ERR(err);
  }

  // PMI(cov);
  
  assert(key->simu.active);
  
  assert(sBR->trend == NULL);
  sBR->trendlen = trendlen;
  if ( (sBR->trend = (double**) CALLOC(trendlen, sizeof(double*))) == NULL)
    RETURN_ERR(ERRORMEMORYALLOCATION);
  for (int j=0; j < trendlen; j++) {
    if ((sBR->trend[j] =  (double*) MALLOC(total * sizeof(double)))==NULL) 
      RETURN_ERR(ERRORMEMORYALLOCATION);
  }

  //  PMI(key);

  
  double *x;
  int len;
  if (keygrid) {
    x = Getxgr(key)[0];
    len = 3;
  } else {
    x = Getx(key);
    len = Gettotalpoints(key);
  }
  if ((err = loc_set(x, NULL, NULL, dim, dim, len, 0, false, keygrid,
		     DistancesGiven(key), sBR->vario)) > NOERROR)
    RETURN_ERR(err);
  //  printf("A\n");
  
  if (sBR->vario->sub[0] != NULL) 
    SetLoc2NewLoc(sBR->vario->sub[0], PLoc(sBR->vario));
  // printf("AB\n");
   
  cov->loggiven = wahr;
  
  pgs_storage *pgs = cov->Spgs;
  for (int d=0; d<dim; d++) {
    pgs->supportmin[d] = RF_NEGINF; // 4 * for debugging...
    pgs->supportmax[d] = RF_INF;
    pgs->supportcentre[d] = RF_NA;      
  }

  //  printf("done\n");
  
  RETURN_NOERROR;
}


#define ORIG_IDX 0
int init_BRorig(model *cov, gen_storage *s){
  assert(cov->mpp.moments >= 1);
  model *key = cov->key;
  if (key == NULL) BUG;
  int err, 
    dim = OWNXDIM(0);
  br_storage *sBR = cov->Sbr;
  pgs_storage *pgs = NULL;

  //PMI(cov->calling->calling);
  
  if ((err = alloc_cov(cov, dim, 1, 1)) != NOERROR) RETURN_ERR(err);
  pgs = cov->Spgs;
  if ((err = general_init(cov, 1, s)) != NOERROR) goto ErrorHandling;
  Variogram(NULL, sBR->vario, sBR->trend[ORIG_IDX]);
   
   assert(MODELNR(key) == GAUSSPROC);
  assert(key->mpp.moments >= 1);
  
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
  cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1.0;
  cov->mpp.maxheights[0] = EXP(GLOBAL.extreme.standardmax);
  
  pgs->zhou_c = 1.0;
  
  
  if ((err = ReturnOwnField(cov)) != NOERROR)  // must be later than INIT !
   goto ErrorHandling;  

 ErrorHandling:
  if (err != NOERROR) br_DELETE(&(cov->Sbr), cov);
  cov->simu.active = cov->initialised = err == NOERROR;
 RETURN_ERR(err);
}

void do_BRorig(model *cov, gen_storage *s) {
  model *key = cov->key;
  assert(key != NULL);
  br_storage *sBR = cov->Sbr;
  assert(sBR != NULL && sBR->trend != NULL);
  double *res = cov->rf,
    *trend = sBR->trend[ORIG_IDX];
  int 
    zeropos = sBR->zeropos;
  long
    totalpoints = Gettotalpoints(cov);
  assert(totalpoints > 0);
  assert(zeropos >= 0 && zeropos <totalpoints);
  assert(cov->origrf);

  DO(key, s); // nicht gatternr,
  double *lgres = key->rf;
  double lgreszeropos = (double) lgres[zeropos]; // wird durch DO veraendert!

  for (int i=0; i<totalpoints; i++) {
    res[i] = lgres[i] - lgreszeropos - trend[i];
    //printf("res = %10g %10g\n", lgres[i], res[i]);
  }
}

int init_BRshifted(model *cov, gen_storage *s) {
  model *key=cov->key;
  if (key == NULL) RETURN_NOERROR;
  int err = NOERROR,
    dim = ANYDIM;
  bool keygrid = Getgrid(key);
  long 
    keytotal = Gettotalpoints(key),
    shiftedloclen = keygrid ? 3 : keytotal,
    trendlenmax = (int) CEIL((double) GLOBAL.br.BRmaxmem / keytotal),
    trendlenneeded = MIN(keytotal, cov->simu.expected_number_simu);
  br_storage *sBR = cov->Sbr;
  pgs_storage *pgs = NULL;
  
  if ((err = alloc_cov(cov, dim, 1, 1)) != NOERROR) goto ErrorHandling;
  pgs = cov->Spgs;
  if ((err = general_init(cov, MIN(trendlenmax, trendlenneeded), s))
      != NOERROR) goto ErrorHandling;
    
  assert(cov->mpp.moments >= 1);
  assert(MODELNR(key) == GAUSSPROC);
  assert(key->mpp.moments >= 1);
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
  cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1.0;
  cov->mpp.maxheights[0] = EXP(GLOBAL.extreme.standardmax);
  pgs->zhou_c = 1.0;
    
  if ((sBR->shift.loc = (double*)
       MALLOC(dim*shiftedloclen*sizeof(double))) == NULL ||
      (sBR->shift.locindex = (int*) MALLOC(sizeof(int) * dim))==NULL
      ) {
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  
  trendlenmax = (int) CEIL((double) GLOBAL.br.BRmaxmem / keytotal);
  trendlenneeded = MIN(keytotal, cov->simu.expected_number_simu);
  
  sBR->shift.memcounter = 0;
  if ((sBR->shift.loc2mem=(int*) MALLOC(sizeof(int)*keytotal))==NULL) {
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  for (long j=0; j<keytotal; j++) sBR->shift.loc2mem[j] = UNSET;    
  
  if ((sBR->shift.mem2loc=(int*) MALLOC(sizeof(int)*sBR->trendlen))==NULL) {
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
  
  for (long j=0; j < sBR->trendlen; j++) sBR->shift.mem2loc[j] = UNSET;
  
  
  err = ReturnOwnField(cov);
  
 ErrorHandling:
  if (err != NOERROR) br_DELETE(&(cov->Sbr), cov);
  cov->simu.active = cov->initialised = err == NOERROR;
 RETURN_ERR(err);
  
}

void indextrafo(long onedimindex, double ** xgr, int dim, int *multidimindex) {
  int d;
  for (d=0; d<dim; d++) {
    multidimindex[d] = onedimindex % (int) xgr[d][XLENGTH];
    onedimindex = onedimindex / (int) xgr[d][XLENGTH];    
  }
}

void do_BRshifted(model *cov, gen_storage *s) {
  br_storage *sBR = cov->Sbr;
  model *key = cov->key;
  assert(cov->key != NULL);
  
  location_type *keyloc = Loc(key);
  long i, k, zeropos, zeroposMdim,
    keytotal = Gettotalpoints(key);
  int d,  trendindex,
    dim = ANYDIM,
    trendlen = sBR->trendlen,
    *locindex = sBR->shift.locindex,
    *mem2loc = sBR->shift.mem2loc,
    *loc2mem = sBR->shift.loc2mem;
  bool keygrid = Getgrid(key); 
  double *shiftedloc = sBR->shift.loc,
    **xgr = keyloc->xgr,
    **trend = sBR->trend,
    *res = cov->rf,
    *lgres = cov->key->rf;
  assert(cov->origrf);

  DO(key, s);
  zeropos = (long) FLOOR(UNIFORM_RANDOM * keytotal);
  
  if (loc2mem[zeropos] != UNSET) {
    trendindex = loc2mem[zeropos];
    if (mem2loc[trendindex] != zeropos) BUG;
  } else {
    if (sBR->shift.memcounter<trendlen) {
      trendindex = sBR->shift.memcounter; 
      sBR->shift.memcounter++;
    } else {
      trendindex = trendlen - 1; 
      loc2mem[mem2loc[trendlen-1]] = UNSET;
      mem2loc[trendlen-1] = UNSET;
    }
    if (keygrid) {
      indextrafo(zeropos, keyloc->xgr, dim, locindex); // to do: ersetzen
      for (d=0; d<dim; d++) {
	shiftedloc[3*d+XSTART]  = -locindex[d] * xgr[d][XSTEP];
	shiftedloc[3*d+XLENGTH] = xgr[d][XLENGTH];
	shiftedloc[3*d+XSTEP]   = xgr[d][XSTEP];
      }
    } else {
      zeroposMdim = zeropos*dim;
      for (k=i=0; i<keytotal; i++) {
	for (d=0; d<dim; d++, k++) {
	  shiftedloc[k] = keyloc->x[k] - keyloc->x[zeroposMdim+d];
	}
      }
    }
 
    partial_loc_set(Loc(sBR->vario), shiftedloc, NULL, 
		    keygrid ? 3: keytotal, 0, keyloc->distances,
		    dim, NULL, keygrid, true);
    if (sBR->vario->sub[0] != NULL) 
        SetLoc2NewLoc(sBR->vario->sub[0], PLoc(sBR->vario));
    Variogram(NULL, sBR->vario, sBR->trend[trendindex]);
    
    mem2loc[trendindex] = zeropos; // todo ? long instead of int
    loc2mem[zeropos] = trendindex;
  }
  
  for (i=0; i<keytotal; i++) {
    res[i] = lgres[i] - lgres[zeropos] - trend[trendindex][i];  
  }
  
}

int check_BRmixed(model *cov) {
  //NotProgrammedYet("at the moment Brown-Resnick processes");
  
  int err;
  br_param *bp = &(GLOBAL.br);
   
  if (!cov->logspeed) SERR("BrownResnick requires a variogram model as submodel that tends to infinity [t rate of at least 4log(h) for being compatible with BRmixed");
  
  kdefault(cov, BR_MESHSIZE, bp->BRmeshsize);
  kdefault(cov, BR_VERTNUMBER, bp->BRvertnumber);
  kdefault(cov, BR_OPTIM, bp->BRoptim);
  kdefault(cov, BR_OPTIMTOL, bp->BRoptimtol);
  kdefault(cov, BR_VARIOBOUND, bp->variobound);
  if (COVNR == BRMIXED_USER && cov->key == NULL && P0INT(BR_OPTIM) > 0) {
    if (!PisNULL(BR_LAMBDA)) {
      if (PisNULL(BR_OPTIMAREA)) SERR1("'%.50s' not given", KNAME(BR_OPTIMAREA));
      if (PL > 0) { PRINTF("'%.50s' set to '0'", KNAME(BR_OPTIM));}
      PINT(BR_OPTIM)[0] = 0;
    } else if (P0INT(BR_OPTIM) == 2 && !PisNULL(BR_OPTIMAREA))
      if (PL > 0) {PRINTF("'%.50s' set to '1'", KNAME(BR_OPTIM));}
  }
  
  if (cov->key != NULL && P0INT(BR_OPTIM) == 2) {
    if (!isIsotropic(SYSOF(cov->key))) {
      //     PMI(cov->key);     
      SERR("area optimisation implemented for the isotropic case only"); //@MARTIN: das scheint nicht zu funktionieren, wenn ich ein Variogramm eingebe
    }  
  }    

  kdefault(cov, BR_LAMBDA, RF_NA);
  if (PisNULL(BR_OPTIMAREA)) kdefault(cov, BR_OPTIMAREA, 0.0);
  
  if ((err = checkBrownResnickProc(cov)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov, true)) != NOERROR) RETURN_ERR(err);
  if (VDIM0 != 1) SERR("BR only works in the univariate case");

  RETURN_NOERROR; 
}

void kappaBRmixed(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  // i nummer des parameters
  
  switch(i) {
    
    case GEV_XI: case GEV_MU: case GEV_S: case BR_MESHSIZE:
    case BR_VERTNUMBER: case BR_OPTIM: case BR_OPTIMTOL: 
    case BR_LAMBDA: case BR_VARIOBOUND:
      *nr = 1;
      *nc = 1;
    break;
    
    case BR_OPTIMAREA:
      *nr = 1;
      *nc = SIZE_NOT_DETERMINED;
    break;

    default:
      *nr = *nc = OUT_OF_RANGE;
  }
}

double IdxDistance(int maxind, int zeropos, double **xgr, int dim) {
  int d,
    delta = 0,
    x = maxind,
    y = zeropos;
  for (d=0; d<dim; d++) {
    delta += ( (x % (int) xgr[d][XLENGTH]) - (y % (int) xgr[d][XLENGTH]))* ( (x % (int) xgr[d][XLENGTH]) - (y % (int) xgr[d][XLENGTH]));
    x /= xgr[d][XLENGTH];
    y /= xgr[d][XLENGTH];
  }
  return SQRT(delta);
}

void OptimArea(model *cov) {
  // side effect auf P(BR_OPTIMAREA) !!
    
  br_storage *sBR = cov->Sbr;
  pgs_storage *pgs = cov->Spgs;
  model *key = sBR->m3.sub[0];
  location_type *keyloc = Loc(key);
  double **xgr = keyloc->xgr,
    *optimarea, dummy, Error, lambda,
    maxErrorbound, Errorbound, Errorboundtmp, 
    Errortol = P0(BR_OPTIMTOL),
    step = P0(BR_MESHSIZE),
    invstepdim, **am=NULL;
  int d, j, k, cellnumber, 
    *cellcounter, idxdist,
    n_zhou_c = pgs->n_zhou_c,
    vertnumber = P0INT(BR_VERTNUMBER),
    minradius = (int) (sBR->m3.minradius/step),
    **countvector = sBR->m3.countvector,
    zeropos = sBR->zeropos,
    dim = ANYDIM,
    keytotal = Gettotalpoints(key);

  if (minradius==0) return;  
    
  if ((cellcounter = (int*) MALLOC((minradius+1) * sizeof(int))) == NULL)
     goto ErrorHandling;
  for (k=0; k<=minradius; k++) cellcounter[k]=0;
  
  if ((am = (double**) CALLOC(vertnumber, sizeof(double*)))==NULL)
    goto ErrorHandling;
  for (j=0; j<vertnumber; j++) {
    if ((am[j] = (double*) MALLOC((minradius+1) * sizeof(double))) == NULL)
      goto ErrorHandling;
  }    
  
  
  for (k=0; k<keytotal; k++) {
     idxdist = (int) CEIL(IdxDistance(k, zeropos, xgr, dim));
     if (idxdist<=minradius) (cellcounter[idxdist])++;
  }
  
  assert(cellcounter[0]==1)  
  
  invstepdim = intpow(step, -dim);
  lambda = 0.0;
  for (j=0; j<vertnumber; j++) lambda += countvector[j][0]*invstepdim / n_zhou_c;
  
  for (d=0; d<=minradius; d++) { //monotonic in each column
    for (j=vertnumber-1; j>=0; j--) {
      am[j][d] = countvector[j][d] * invstepdim / (n_zhou_c * cellcounter[d]);  
      for (k=j+1; k<vertnumber; k++) {
	if (am[j][d] < am[k][d]) {
	  dummy = am[j][d];
	  am[j][d] = am[k][d];
	  am[k][d] = dummy;
	}
      }
      am[j][d] = FMIN(am[j][d], lambda/vertnumber); //cutoff
    }
  }
  
  for (j=0; j<vertnumber; j++) { //monotonic in each row
    for (d=0; d<=minradius; d++) {
      for (k=d+1; k<=minradius; k++) {
	if (am[j][d] < am[j][k]) {
	  dummy = am[j][d];
	  am[j][d] = am[j][k];
	  am[j][k] = dummy;
	}
      }
    }
  }
  
  maxErrorbound = 0.0;
  for (d=0; d<=minradius; d++) { //areamatrix => Error matrix
    for (j=0; j<vertnumber; j++) {
      am[j][d] = lambda/vertnumber - am[j][d];
      maxErrorbound = FMAX(maxErrorbound, am[j][d]);
    }
  }
  
  cellnumber = 0;
  Error = 0.0;
  Errorboundtmp = Errorbound = 0.0;  
  while (Errorboundtmp < maxErrorbound && Error < Errortol) {
    Errorbound = Errorboundtmp; 
    Errorboundtmp = maxErrorbound;
    for (d=0; d<=minradius; d++) {
      for (j=0; j<vertnumber; j++) {
	if (am[j][d] > Errorbound) {
	  Errorboundtmp = FMIN(Errorboundtmp, am[j][d]); 
	}
      }      
    }
    Error = 0.0;
    for (d=0; d<=minradius; d++) {
      for (j=0; j<vertnumber; j++) {
	if (am[j][d] <= Errorboundtmp + 1e-6) {
	  Error += (double) cellcounter[d] * am[j][d];
	  cellnumber += cellcounter[d];
	}
      }
    }
    Error = Error / (cellnumber*lambda/vertnumber);
  }

  PFREE(BR_OPTIMAREA);
  PALLOC(BR_OPTIMAREA, 1, minradius+1);
  optimarea = P(BR_OPTIMAREA);
  for (d=0; d<=minradius; d++) {
    optimarea[d] = 0.0;
    for (j=0; j<vertnumber; j++) {
      if (am[j][d] <= Errorbound + 1e-6) {
        optimarea[d] += 1.0/vertnumber;
      }
    }
    //printf("%10g ", optimarea[d]);
  }
  //printf("\n");
  
  
  
  ErrorHandling:
  if (cellcounter != NULL) FREE(cellcounter);
  if (am != NULL) {
    for (j=0; j<vertnumber; j++) {
      if (am[j] != 0) FREE(am[j]);   
    }  
    FREE(am);  
  }    
  
  return;
}


void set_lowerbounds(model *cov) {    
  br_storage *sBR = cov->Sbr;
  assert(sBR != NULL);
  assert(sBR->m3.sub[0] != NULL);
  double step = P0(BR_MESHSIZE),
    *optimarea = P(BR_OPTIMAREA);    
  int j, k,
    dim = ANYDIM,
    minradius = (int) (sBR->m3.minradius / step);
  model *key = sBR->m3.sub[0];
  location_type *keyloc = Loc(key);
  double **xgr = keyloc->xgr;
  long keytotal = Gettotalpoints(key);
 
  for (j=0; j<keytotal; j++) {
    sBR->m3.lowerbounds[j] = RF_INF;
    k = (int) CEIL(IdxDistance(j, sBR->zeropos, xgr, dim));
    if (k <= minradius && optimarea[k]>1e-5) {
      sBR->m3.lowerbounds[j] = -LOG(optimarea[k]);  
    }  
    //printf("%10g ", sBR->m3.lowerbounds[j]);
  }
  //printf("\n");
}

int prepareBRoptim(model *cov) {
  br_storage *sBR = cov->Sbr;
  model *key = sBR->m3.sub[0];
  location_type *keyloc = Loc(key);
  double  step = P0(BR_MESHSIZE),
        **xgr = keyloc->xgr;
  int i, j, d,
    vertnumber = P0INT(BR_VERTNUMBER),
    dim = ANYDIM,
    maxradius = 1,
    minradius = (int) (sBR->m3.minradius/step);

  for (d=0; d<dim; d++) {
    maxradius *= (int) FLOOR(xgr[d][XLENGTH] / 2.0) + 1;
  }
 
  switch(P0INT(BR_OPTIM)) {
  case 0:   
    if (ISNAN(P0(BR_LAMBDA))) P(BR_LAMBDA)[0] = 1.0;
    break;
  case 1: // nothing to do
    break;
  case 2:
    sBR->m3.vertnumber = vertnumber; 
    
    if (sBR->m3.countvector != NULL || sBR->m3.areamatrix != NULL) BUG;
    if ((sBR->m3.countvector = (int**) CALLOC(vertnumber, sizeof(int*)))==NULL || 
        (sBR->m3.logvertnumber = (double *) MALLOC(vertnumber * sizeof(double)))
         == NULL)
      RETURN_ERR(ERRORMEMORYALLOCATION);
    for (j=0; j<vertnumber; j++) {
      if ((sBR->m3.countvector[j] = (int*) CALLOC(minradius + 1, sizeof(int))) 
          == NULL)
        RETURN_ERR(ERRORMEMORYALLOCATION);
      for (i=0; i<=minradius; i++) sBR->m3.countvector[j][i] = 0;
    }
    for (j=0; j<vertnumber; j++)
      sBR->m3.logvertnumber[j] = - LOG((double) (j+1)/vertnumber);    
    break;

  default:
    SERR("optimization might not be used here\n");
  }
  
  if ((sBR->m3.areamatrix  = (double *) MALLOC((minradius + 1)* sizeof(double))) 
         == NULL) {
    RETURN_ERR(ERRORMEMORYALLOCATION);
  }
  sBR->m3.areamatrix[0] = 1.0;
  if (minradius > 0) {
    for (i=1; i<=minradius; i++) {
      if (i <= cov->ncol[BR_OPTIMAREA]) {
        sBR->m3.areamatrix[i] = P(BR_OPTIMAREA)[i-1];  
      }  
      else sBR->m3.areamatrix[i] = 0.0;
    }
  }
  
  PFREE(BR_OPTIMAREA);
  PALLOC(BR_OPTIMAREA, 1, minradius + 1);
  double *optimarea = P(BR_OPTIMAREA);
  for (i=0; i<=minradius; i++) {
    optimarea[i] = sBR->m3.areamatrix[i];
  }
  set_lowerbounds(cov);
  
  if (PL >= PL_STRUCTURE) {PRINTF("BR optimisation finished...\n");}

 RETURN_NOERROR;	
}



int init_BRmixed(model *cov, gen_storage *s) {
  location_type *loc = Loc(cov);
  br_storage *sBR = cov->Sbr;
  assert(sBR != NULL);
  assert(sBR->m3.sub[0] != NULL);
  model *key = sBR->m3.sub[0];
  location_type *keyloc = Loc(key);
  int  d, err = NOERROR, 
    dim = ANYDIM,
    bytes = sizeof(double) * dim,
    keytotal = Gettotalpoints(key);
  double area = 1.0,
    step = P0(BR_MESHSIZE);
  pgs_storage *pgs = NULL;
   
  assert(isPointShape(cov));
  assert(dim > 0);
  
  if ((err = alloc_cov(cov, dim, 1, 1)) != NOERROR) goto ErrorHandling;//nur pgs
  pgs = cov->Spgs; // nach alloc_cov !!
  if ((err = general_init(cov, 1, s)) != NOERROR) goto ErrorHandling;
  Variogram(NULL, sBR->vario, sBR->trend[0]);


  if ((sBR->m3.suppmin = (double*) MALLOC(bytes))==NULL ||
      (sBR->m3.suppmax = (double*) MALLOC(bytes))==NULL ||      
      (sBR->m3.loccentre = (double*) MALLOC(bytes))==NULL 
      //     (sBR->shift.locindex = (int*) MALLOC(sizeof(int) * dim))==NULL ||
     ) {
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
 
  assert(Getcaniso(cov) == NULL);
  GetDiameter(loc, sBR->m3.suppmin, sBR->m3.suppmax, sBR->m3.loccentre);
  assert(pgs->supportmin != NULL);
  for (d=0; d<dim; d++) {
    pgs->own_grid_cumsum[d] =
      d==0 ? 1 : pgs->own_grid_cumsum[d-1] * pgs->own_grid_len[d-1];      
    sBR->m3.suppmin[d] = FLOOR((sBR->m3.suppmin[d] - sBR->m3.radius -
				sBR->m3.minradius) / step) * step - step / 2;  
    sBR->m3.suppmax[d] = CEIL((sBR->m3.suppmax[d] + sBR->m3.radius +
			       sBR->m3.minradius) / step) * step + step / 2;  
    area *= sBR->m3.suppmax[d] - sBR->m3.suppmin[d];
    pgs->own_grid_start[d] = RF_NEGINF;
    pgs->own_grid_step[d] = keyloc->xgr[d][XSTEP];
    pgs->own_grid_len[d] = keyloc->xgr[d][XLENGTH];
  }
  
  for (d=0; d<=cov->mpp.moments; d++) {
    cov->mpp.mM[d] = cov->mpp.mMplus[d] = 1.0;
  }
  cov->mpp.maxheights[0] = EXP(0.0);
    
  assert(MODELNR(key) == GAUSSPROC);

  assert(cov->mpp.moments >= 1);
  assert(Getgrid(key));
  
  assert(key->mpp.moments >= 1);
  key->mpp.mM[0] = key->mpp.mMplus[0] = 1.0;
  key->mpp.mM[1] = key->mpp.mMplus[1] = 1.0;

  
  if ((sBR->m3.lowerbounds = (double*) MALLOC(keytotal*sizeof(double))) == NULL) { 
    err=ERRORMEMORYALLOCATION; 
    goto ErrorHandling; 
  }
  
  if ((err = prepareBRoptim(cov)) != NOERROR) {
     goto ErrorHandling;
  }

  assert(keyloc != NULL);
  key->simu.active = true;

  set_lowerbounds(cov);   
  cov->rf = sBR->m3.sub[0]->rf;
  cov->origrf = false;
  cov->fieldreturn = Huetchenownsize;
  pgs->estimated_zhou_c = (bool) ISNAN(P0(BR_LAMBDA));
  pgs->zhou_c =  pgs->estimated_zhou_c ? 0.0 : P0(BR_LAMBDA)*area; //@MARTIN: area kann weg, falls in logdens
  pgs->logmean = false;
  pgs->sq_zhou_c = pgs->sum_zhou_c = 0.0;
  pgs->n_zhou_c = 0;
  sBR->m3.next_am_check = GLOBAL.br.deltaAM;

 
 ErrorHandling:
  if (err != NOERROR) br_DELETE(&(cov->Sbr), cov);
  cov->simu.active = cov->initialised = err == NOERROR;

  RETURN_ERR(err);

}

void do_BRmixed(model *cov, gen_storage *s) {
  // to do: improve simulation speed by dynamic sizes
  assert(cov->key!=NULL);
  br_storage *sBR = cov->Sbr;
  model  *key = sBR->m3.sub[0];
  assert(cov->rf == key->rf);
  location_type *keyloc = Loc(key);
  assert(Getgrid(key)); 
  pgs_storage *pgs = cov->Spgs;
  assert(pgs != NULL); 

  int d, 
    i, j, maxind, idxdist,
    dim = ANYDIM,
    hatnumber=0,
    lgtotalpoints = Gettotalpoints(key),
    zeropos = sBR->zeropos,
    vertnumber = P0INT(BR_VERTNUMBER);

  double step = P0(BR_MESHSIZE),
    invstepdim = intpow(step, -dim),
    uplusmaxval , maxval, u[MAXMPPDIM], 
    ** xgr = keyloc->xgr,
    *lowerbounds = sBR->m3.lowerbounds,
    area=1.0, *lgres = key->rf, //@MARTIN: area obsolet, falls in logdens
    *trend = sBR->trend[0];

  int minradius = (int) (sBR->m3.minradius/step); 

      
  if (P0INT(BR_OPTIM) == 2 &&  pgs->n_zhou_c >= sBR->m3.next_am_check) {      
    sBR->m3.next_am_check += GLOBAL.br.deltaAM; 
    OptimArea(cov); 
    set_lowerbounds(cov);
  }
  
  for (d=0; d<dim; d++) {
    u[d] = ROUND((UNIFORM_RANDOM*(sBR->m3.suppmax[d] - sBR->m3.suppmin[d]) +
		  sBR->m3.suppmin[d]) / step) * step;
    area *= sBR->m3.suppmax[d] - sBR->m3.suppmin[d];
    pgs->supportmin[d] = u[d] - sBR->m3.minradius - sBR->m3.radius; 
    pgs->supportmax[d] = u[d] + sBR->m3.minradius + sBR->m3.radius; 
    pgs->supportcentre[d] = u[d];
    pgs->own_grid_start[d] = keyloc->xgr[d][XSTART] + u[d];
  }  
  
  while(true) {
    DO(key, s);
    hatnumber++;
    maxval = RF_NEGINF;
    maxind = 0;
    for (i=0; i<lgtotalpoints; i++) {
      lgres[i] -= trend[i];
      if (lgres[i] > maxval) {
	maxval = lgres[i];
	maxind = i;
      }
    }
    if (maxind == zeropos) {
      pgs->sq_zhou_c += area * invstepdim * area * invstepdim; // @MARTIN: s.o
      pgs->sum_zhou_c += area * invstepdim;                    //
    }
    
    uplusmaxval = maxval - lgres[zeropos] - LOG(UNIFORM_RANDOM);

    if (P0INT(BR_OPTIM) == 2) {
      for (j=0; j<vertnumber; j++) {
        if (uplusmaxval > sBR->m3.logvertnumber[j]) {
          idxdist = (int) CEIL(IdxDistance(maxind, zeropos, xgr, dim));
          if (idxdist<=minradius) (sBR->m3.countvector[j][idxdist])++;
          break;
        }
      }
    }

    if (uplusmaxval > lowerbounds[maxind]) break;    
  }
  
  pgs->n_zhou_c += hatnumber;

  if (PL >= PL_STRUCTURE && hatnumber > 300) {
    PRINTF("note: large hat number (%d) might indicate numerically suboptimal framework\n",
	   hatnumber);
  }
  
  //shifting maximum to origin is not necessary because of stationarity
  //(conditional on T=t, the "correct" shape function is shifted which yields
  //the same stationary max-stable process; OK because there is a finite number
  //of values for t!)

  for (i=0; i<lgtotalpoints; i++) {
    lgres[i] -= maxval;
  }

  return;
    
}


void range_BRmixed(model *cov, range_type *range) {
  range_mpp(cov, range);
  
  range->min[BR_MESHSIZE] = 0;
  range->max[BR_MESHSIZE] = RF_INF;
  range->pmin[BR_MESHSIZE] = 0;
  range->pmax[BR_MESHSIZE] = RF_INF;
  range->openmin[BR_MESHSIZE] = true;
  range->openmax[BR_MESHSIZE] = true; 
  
  range->min[BR_VERTNUMBER] = 1;
  range->max[BR_VERTNUMBER] = RF_INF;
  range->pmin[BR_VERTNUMBER] = 1;
  range->pmax[BR_VERTNUMBER] = 50;
  range->openmin[BR_VERTNUMBER] = false;
  range->openmax[BR_VERTNUMBER] = false;
  
  range->min[BR_OPTIM] = 0;
  range->max[BR_OPTIM] = 2;
  range->pmin[BR_OPTIM] = 0;
  range->pmax[BR_OPTIM] = 2;
  range->openmin[BR_OPTIM] = false;
  range->openmax[BR_OPTIM] = false;
  
  range->min[BR_OPTIMTOL] = 0;
  range->max[BR_OPTIMTOL] = 1;
  range->pmin[BR_OPTIMTOL] = 0;
  range->pmax[BR_OPTIMTOL] = 0.1;
  range->openmin[BR_OPTIMTOL] = true;
  range->openmax[BR_OPTIMTOL] = true;
  
  range->min[BR_LAMBDA] = 0;
  range->max[BR_LAMBDA] = RF_INF;
  range->pmin[BR_LAMBDA] = 0;
  range->pmax[BR_LAMBDA] = RF_INF;
  range->openmin[BR_LAMBDA] = true;
  range->openmax[BR_LAMBDA] = true;
  
  range->min[BR_OPTIMAREA] = 0;
  range->max[BR_OPTIMAREA] = 1;
  range->pmin[BR_OPTIMAREA] = 0;
  range->pmax[BR_OPTIMAREA] = 1;
  range->openmin[BR_OPTIMAREA] = false;
  range->openmax[BR_OPTIMAREA] = false;
  
  range->min[BR_VARIOBOUND] = 0;
  range->max[BR_VARIOBOUND] = RF_INF;
  range->pmin[BR_VARIOBOUND] = 2;
  range->pmax[BR_VARIOBOUND] = 25;
  range->openmin[BR_VARIOBOUND] = false;
  range->openmax[BR_VARIOBOUND] = true;
}

int structBRuser(model *cov, model **newmodel) {
 
  location_type *loc = Loc(cov);
  model *sub = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  int i, d, err, model_intern, 
    dim = ANYDIMOF(sub),
    newxlen;
  bool grid;
  double centreloc[MAXMPPDIM], minloc[MAXMPPDIM], maxloc[MAXMPPDIM],
    *newx= NULL, 
    **xgr = loc->xgr;

  ASSERT_NEWMODEL_NULL;
  assert(isBrMethod(cov));  
  
  if (loc->Time || (Getgrid(cov) && loc->caniso != NULL)) {
    TransformLoc(cov, false, GRIDEXPAND_AVOID, false); // changes loc !
    SetLoc2NewLoc(sub, PLoc(cov));
  }
  
  loc = Loc(cov);
  grid = Getgrid(cov);
  
  model_intern = (COVNR == BRORIGINAL_USER) ? BRORIGINAL_INTERN
               : (COVNR == BRMIXED_USER) ? BRMIXED_INTERN
               : (COVNR == BRSHIFTED_USER) ? BRSHIFTED_INTERN
               : BRORIGINAL_USER;
	       
  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);// should be ok
  ONCE_NEW_STORAGE(gen);
  
  assert(Getcaniso(cov) == NULL);
  GetDiameter(loc, minloc, maxloc, centreloc); 
  int totalpoints = loc->totalpoints;
  newxlen = grid ? 3 : totalpoints;
  if ((newx = (double*) MALLOC(dim * newxlen * sizeof(double))) == NULL) {
     GERR("Memory allocation failed.\n"); 
  } 
  if (grid) {
    for (d=0; d<dim; d++) {
      newx[3 * d + XSTART] = xgr[d][XSTART] - centreloc[d] 
           + (((int) xgr[d][XLENGTH]) % 2 == 0) * xgr[d][XSTEP]/2;
      newx[3 * d + XSTEP] = xgr[d][XSTEP];
      newx[3 * d + XLENGTH] = xgr[d][XLENGTH];
   }
  } else {
    for (i=0; i<totalpoints; i++)
      for(d=0; d<dim; d++) newx[i*dim+d] = loc->x[i*dim+d] - centreloc[d];
  }

  if ((err = loc_set(newx, NULL, dim, dim, newxlen, false, grid, // loc->grid
                     loc->distances,  cov)) > NOERROR)
    goto ErrorHandling;
  SetLoc2NewLoc(sub, PLoc(cov));

  if ((err=covcpy(&(cov->key), sub)) > NOERROR) goto ErrorHandling;
  
  if (cov->sub[MPP_TCF] != NULL) {
     if ((err = STRUCT(sub, &(cov->key))) > NOERROR) goto ErrorHandling;
     assert(cov->key->calling == cov);
  } 

  addModel(&(cov->key), model_intern);
  assert(cov->key->calling == cov);
  
  kdefault(cov->key, GEV_XI, P0(GEV_XI));
  kdefault(cov->key, GEV_MU, P0(GEV_MU));
  kdefault(cov->key, GEV_S, P0(GEV_S));

  if (COVNR == BRMIXED_USER) {
    kdefault(cov->key, BR_MESHSIZE, P0(BR_MESHSIZE));
    kdefault(cov->key, BR_VERTNUMBER, P0INT(BR_VERTNUMBER));
    kdefault(cov->key, BR_OPTIM, P0INT(BR_OPTIM));
    kdefault(cov->key, BR_OPTIMTOL, P0(BR_OPTIMTOL));
    kdefault(cov->key, BR_VARIOBOUND, P0(BR_VARIOBOUND));
    kdefault(cov->key, BR_LAMBDA, P0(BR_LAMBDA));
    if (!PisNULL(BR_OPTIMAREA)) {
      PARAMALLOC(cov->key, BR_OPTIMAREA, cov->nrow[BR_OPTIMAREA],
                 cov->ncol[BR_OPTIMAREA]);

      PCOPY(cov->key, cov, BR_OPTIMAREA);
    }
  }
  
  if ((err = CHECK(cov->key, OWNLOGDIM(0), OWNXDIM(0), PointShapeType,
		   OWNDOM(0), OWNISO(0), 1, BrMethodType)) == NOERROR) {
    if ((err = STRUCT(cov->key, NULL)) <= NOERROR) {
      err = CHECK(cov->key, OWNLOGDIM(0), OWNXDIM(0),
		  PointShapeType, OWNDOM(0), OWNISO(0), 1, BrMethodType);
    }
  }
 
   
  ErrorHandling:
  FREE(newx);
  RETURN_ERR(err);
  
}


// new
int structBRintern(model *cov, model **newmodel) {
  model *sub = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  location_type *loc = Loc(cov);
  int i, d, j, err,
    dim = ANYDIM, // sub->ts dim, 
    totaldim = Gettotalpoints(cov) * dim,
    zeropos = 0, // default, mostly overwritten
    newxlen = 3; // default (grid)
  bool grid = Getgrid(cov);
  double step, mindist, dist,
    *newx = NULL,
    **xgr = loc->xgr;
  br_storage *sBR = NULL;
  model *submodel = NULL;
  
  ASSERT_NEWMODEL_NULL;
 
  assert(isPointShape(cov));
    
  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);// should be ok
  ONCE_NEW_STORAGE(gen);
  NEW_STORAGE_WITH_SAVE(br);
  sBR = cov->Sbr;
  sBR->nr = COVNR;

  if (cov->sub[MPP_TCF] != NULL) {
    if ((err = STRUCT(sub, &(cov->key))) > NOERROR) goto ErrorHandling;
    assert(cov->key->calling == cov);
    // SET_CALLING(cov->key, cov);
 } else if ((err=covcpy(&(cov->key), sub)) > NOERROR) goto ErrorHandling;
 
  if ((err = CHECK(cov->key,OWNLOGDIM(0), OWNXDIM(0), VariogramType, OWNDOM(0),
		   SYMMETRIC, 1, EvaluationType)) != NOERROR)
    goto ErrorHandling;

  if ((err = covcpy(&submodel, cov->key)) != NOERROR)
    goto ErrorHandling;

  if ((err = CHECK(submodel, 1, 1, VariogramType, XONLY, 
		   ISOTROPIC, 1, EvaluationType)) != NOERROR)
    goto ErrorHandling;
  
  if ((err = newmodel_covcpy(&(sBR->vario), VARIOGRAM_CALL, cov->key))!=NOERROR)
      goto ErrorHandling;
  if ((err = alloc_cov(sBR->vario, dim, 1, 1)) != NOERROR) goto ErrorHandling;
	
  addModel(&(cov->key), GAUSSPROC, cov);
  assert(cov->key->ownloc == NULL);
  
  if (COVNR == BRORIGINAL_INTERN) {
    if (!grid) {
      int totalpoints = loc->totalpoints;
      zeropos = totalpoints;
      for (i=0; i<totalpoints; i++) {
	double norm = 0.0;
	for (d=0; d<dim; d++) {
	  double dummy = loc->x[i*dim+d];
	  norm += dummy * dummy;
	}
	if (norm < 1e-8) zeropos = i;
      }
      int newtotpoints = totalpoints + (zeropos == totalpoints);
      if ((newx = (double*) MALLOC(dim*newtotpoints*sizeof(double))) == NULL) {
        err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      for(d=0; d<dim; d++) {
        for (i=0; i<totalpoints; i++) newx[i*dim+d] = loc->x[i*dim+d];
	if (zeropos == totalpoints) newx[totalpoints*dim+d] = 0.0;
      }
    } 
  } else if (COVNR == BRMIXED_INTERN) {

    step = P0(BR_MESHSIZE);
    mindist = RF_INF;
    if (grid) {
      for (d=0; d<dim; d++)  
        if (xgr[d][XSTEP] < mindist)
	  mindist = xgr[d][XSTEP];
    } else {
      int maxtotal = totaldim;
      if (maxtotal > 1000) {
	warning("minimal distance is only estimated");
	maxtotal = 1000; // to do
      }
      for (i=0; i<maxtotal; i+=dim)
        for (j=i+dim; j<maxtotal; )
	  for (d=0; d<dim; d++, j++) {
          dist = FABS(loc->x[i+d] - loc->x[j]);
          if (dist > 1e-15 && dist < mindist) mindist = dist;
	  }
      grid = true;
    }

    newxlen = 3;
 
    if (mindist < step) {     
      PRINTF("Warning! meshsize is larger than resolution of simulation! meshsize is automatically decreased to %10g.\n", mindist);
      P(BR_MESHSIZE)[0] = step = mindist;
    }
    
    if (!PisNULL(BR_OPTIMAREA)) {
      sBR->m3.minradius =  cov->ncol[BR_OPTIMAREA] * step;
    } else {
      sBR->m3.minradius = 0;
    }    
    
    double yy, C0, gammamin,
      alpha,
      xx=step * 1e-6;
    COV(ZERO(submodel), submodel, &C0);
    COV(&xx, submodel, &yy);
    alpha = LOG(C0 - yy) / LOG(xx);
    if (alpha > 2.0) alpha = 2.0;
    gammamin = 4.0 - 1.5 * alpha;
    
    INVERSE(&gammamin, submodel, &xx);
    xx = CEIL(xx/step) * step;
    sBR->m3.minradius = FMAX(sBR->m3.minradius, xx); 

    yy = P0(BR_VARIOBOUND);
    INVERSE(&yy, submodel, &(sBR->m3.radius));
    
    sBR->m3.radius = CEIL(sBR->m3.radius / step) * step;

 
    if ((newx = (double*) MALLOC(newxlen*dim*sizeof(double))) == NULL) {
      err = ERRORMEMORYALLOCATION; goto ErrorHandling;
    }

    if ((err = covcpy(sBR->m3.sub, cov->key)) != NOERROR) goto ErrorHandling;
    for (d=0; d<dim; d++) {
      newx[3*d+XSTART] = -sBR->m3.radius - sBR->m3.minradius;
      newx[3*d+XLENGTH] =
	2 * ((int) ROUND((sBR->m3.radius + sBR->m3.minradius) / step)) + 1;
      newx[3*d+XSTEP] = step;
    }

    err = loc_set(newx, NULL, dim, dim, newxlen, false, grid, 
                  false, sBR->m3.sub[0]);
    double **subxgr = Loc(sBR->m3.sub[0])->xgr;
    zeropos = 0;
    for (d = dim; d > 0; d--) {
      double len =  subxgr[d-1][XLENGTH];
      zeropos = zeropos * (int) len + (int) CEIL(len * 0.5) - 1;
    }
    sBR->zeropos = zeropos;

    if ((err = CHECK(sBR->m3.sub[0], OWNLOGDIM(0), OWNXDIM(0), ProcessType,
		     OWNDOM(0), OWNISO(0), 1, GaussMethodType)) == NOERROR) {
       if ((err = STRUCT(sBR->m3.sub[0], NULL)) <= NOERROR) {
         err = CHECK(sBR->m3.sub[0], OWNLOGDIM(0), OWNXDIM(0),  ProcessType,
		     OWNDOM(0), OWNISO(0), 1, GaussMethodType); 
       }
     }
     if (err > NOERROR) goto ErrorHandling;
     assert(sBR->m3.sub[0]->calling == cov);

  } else { //  END BRMIXED;  START SHIFTED    
    assert(COVNR == BRSHIFTED_INTERN);
  }
 

  if (newx == NULL) {
    if ((err = loc_set(grid ? * loc->xgr : loc->x, NULL, dim, dim,
		       grid ? 3 : Gettotalpoints(cov), false, grid,
		       loc->distances, cov->key)) != NOERROR)
      goto ErrorHandling;
  } else {
    if ((err = loc_set(newx, NULL, dim, dim, newxlen, false, grid,
		       false, cov->key)) != NOERROR) goto ErrorHandling;
  }
  xgr = loc->xgr;

 
  if (COVNR != BRMIXED_INTERN && grid) {
    double **subxgr = Loc(cov->key)->xgr;
    for (d=dim; d>0; d--) {
      double len =  subxgr[d-1][XLENGTH];
      zeropos = zeropos * len + (int) CEIL(len * 0.5) - 1;
    }
    sBR->zeropos = zeropos;
  }
  
  if ((err = CHECK(cov->key, OWNLOGDIM(0), OWNXDIM(0), ProcessType, OWNDOM(0),
		   OWNISO(0), 1, GaussMethodType)) == NOERROR) {
    if ((err = STRUCT(cov->key, NULL)) <= NOERROR) {
      err = CHECK(cov->key, OWNLOGDIM(0), OWNXDIM(0), ProcessType, OWNDOM(0),
		  OWNISO(0), 1, GaussMethodType); 
    }
  }
  if (err > NOERROR) goto ErrorHandling;

  assert(cov->key->calling == cov);
  
  ErrorHandling:

  if (submodel != NULL) COV_DELETE(&submodel, cov);  // ok
  
  FREE(newx);

  RETURN_ERR(err);
}

int structBrownResnick(model *cov, model **newmodel) {
  
  int d, err, meth, 
    dim = ANYDIM;
  double  maxcov,
      minloc[MAXMPPDIM], maxloc[MAXMPPDIM],
      centreloc[MAXMPPDIM], maxdist[MAXMPPDIM];      
  model *next = cov->sub[0];
  location_type *loc = Loc(cov);

  if (loc->Time || (Getgrid(cov) && loc->caniso != NULL)) {
    TransformLoc(cov, false, GRIDEXPAND_AVOID, false);
    SetLoc2NewLoc(next, PLoc(cov));
  }
  loc = Loc(cov);
  
  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);// should be ok
    
  if (hasSmithFrame(cov)) {
    if (!cov->logspeed) 
      SERR2("'%.50s' requires a variogram model as submodel that tends to infinity with rate of at least 4log(h) for being compatible with '%.50s'", NICK(cov), DefList[SMITHPROC].nick);
    model *calling = cov->calling;
    double newscale, 
      factor =  INVSQRTTWO;
      
    ASSERT_NEWMODEL_NULL;
    if (next->full_derivs < 0) SERR("given submodel does not make sense");
 
    while (isDollar(next)) {
      addModel(&(cov->key), DOLLAR, cov);
      newscale = 1.0;
      if (!PARAMisNULL(next, DSCALE) ) newscale *= PARAM0(next, DSCALE);
      if (!PARAMisNULL(next, DVAR) )  newscale /= SQRT(PARAM0(next, DVAR));
      if (factor != 1.0) {
	newscale *= factor;
	factor = 1.0;
      }
      RETURN_ERR(ERRORNOTPROGRAMMEDYET);
      

      kdefault(calling, DSCALE, newscale);
      next = next->sub[0]; // ok
    }

    if (cov->sub[MPP_TCF] != NULL) {
      RETURN_ERR(ERRORNOTPROGRAMMEDYET);
    } 

    if (NEXTNR == BROWNIAN && PARAM0(next, BROWN_ALPHA) == 2.0) {
      addModel(&(cov->key), GAUSS, cov);   // ?? 
      if (factor != 1.0) {
	addModel(&(cov->key), DOLLAR, cov);
	kdefault(cov->key, DSCALE, factor);      
      }
    } else {
      SERR("Smith process with BrownResnick tcf only possible for fractal Brownian motion with alpha=2");
    }
  } else if (hasBrMethodFrame(cov) || hasInterfaceFrame(cov) ||
	     hasNormedProcessFrame(cov)) {
    if (hasBrMethodFrame(next)) {
      SERR1("submodel of '%.50s' must be a covariance model or tcf", 
	    NICK(cov));
    } else {
      ASSERT_FRAME_DEFINED(next);  
      Types frame = isnowVariogram(next) ? EvaluationType : BadType;
	
      if (((err = covcpy(&(cov->key), next)) != NOERROR)
	  || ((err = CHECK(cov->key, OWNLOGDIM(0), OWNXDIM(0),
			  VariogramType, XONLY,
			  SYMMETRIC, 1, frame)) != NOERROR)) {
	RETURN_ERR(err);
      }

      assert(Getcaniso(cov) == NULL);
      GetDiameter(loc, minloc, maxloc, centreloc);
      for (d=0; d<MAXMPPDIM; d++) maxdist[d] = 0.5*(maxloc[d] - minloc[d]);
   
      model *K = NULL;
      if ((err = newmodel_covcpy(&K, VARIOGRAM_CALL, cov->key, maxdist, NULL,
				 NULL, dim, dim, 1, 0, false, false, false))
	  != NOERROR) RETURN_ERR(err);
      if ((err = alloc_cov(K, dim, 1, 1)) != NOERROR) RETURN_ERR(err);
      if (K->sub[0] != NULL) SetLoc2NewLoc(K->sub[0], PLoc(K));
      Variogram(NULL, K, &maxcov);

      COV_DELETE(&K, cov);
      if (isnowPosDef(next) || maxcov <= 4.0) {
	meth = BRORIGINAL_USER;  
      } else if (!next->logspeed || next->logspeed <= 4.0 || maxcov <= 10.0) {
	meth = BRSHIFTED_USER;
      } else {
	meth = BRMIXED_USER;
      }

      addModel(&(cov->key), meth, cov);
      model *key = cov->key;
      key->prevloc = PLoc(cov);
      
      kdefault(key, GEV_XI, P0(GEV_XI));
      kdefault(key, GEV_MU, P0(GEV_MU));
      kdefault(key, GEV_S, P0(GEV_S));
       
      if ((err =  CHECK(key, OWNLOGDIM(0), OWNXDIM(0), BrMethodType,
			OWNDOM(0), OWNISO(0),
	                1, BrMethodType)) == NOERROR) {
      if ((err = STRUCT(key, NULL)) <= NOERROR) {
        err = CHECK(key, OWNLOGDIM(0), OWNXDIM(0), BrMethodType,
		    OWNDOM(0), OWNISO(0),
		    1, BrMethodType);
	}
      }
      if (err > NOERROR) RETURN_ERR(err);
    }
  } else {
    ILLEGAL_FRAME;
  }

  // need check to be called?

  RETURN_NOERROR;

  
}

int initBrownResnick (model *cov, gen_storage *S) {

  model *sub = cov->key;
  if (sub == NULL)
    sub = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF]; 
  int err;

  assert(COVNR == BROWNRESNICKPROC);

  if (cov->key != NULL) {
    sub->simu.active = true;
    sub->simu.expected_number_simu = cov->simu.expected_number_simu;
    if ((err = INIT(sub, 0, S)) != NOERROR) RETURN_ERR(err);
    ReturnOtherField(cov, sub);
    //      cov->fieldreturn = wahr;
    //cov->origrf = false;
    //cov->rf = sub->rf;
  }

  cov->simu.active = cov->initialised = true;
 RETURN_NOERROR;

}

void doBrownResnick(model *cov, gen_storage *s) {

  assert(!cov->origrf);
  assert(cov->key != NULL);
  model *key = cov->key;
 
  PL++;
  DO(key, s); // nicht gatternr
  PL--;

}

void finaldoBrownResnick(model *cov, double *res, int n, gen_storage *s) {
  model *key = cov->key;
  finalmaxstable(key, res, n, s);
}


int initBRuser (model *cov, gen_storage *S) {
  location_type *loc = Loc(cov);
  model *sub = cov->key;
  if (sub == NULL) 
    sub = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  int  err, maxpoints = GLOBAL.extreme.maxpoints;
  
  assert(isBrMethod(cov));

  if (loc->distances) RETURN_ERR(ERRORFAILED);
  
  if (cov->key != NULL) {
    sub->simu.active = true;
    double ens = ((double) cov->simu.expected_number_simu) * maxpoints;
    sub->simu.expected_number_simu = (int) MIN(ens, (double) MAXINT);
    
    if ((err = INIT(sub, 1, S)) != NOERROR) RETURN_ERR(err);
    ReturnOwnField(cov); 
  }
  
  cov->simu.active = cov->initialised = true;
  RETURN_NOERROR;
}

#define ABS std::abs

#define NORMED_PROB 0 
#define NORMED_OPTIMP 1 // true/false
#define NORMED_NTH 2
#define NORMED_BURNIN 3
#define NORMED_REJECTION 4

#define NORMED_MAX 0
#define NORMED_VALUES 4


void kappabrnormed(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  // i nummer des parameters
  *nc = *nr = 1;
  if (i == NORMED_PROB) *nr = SIZE_NOT_DETERMINED;
  else if (i > NORMED_REJECTION) *nr = *nc = OUT_OF_RANGE;
  //  printf("kappa%d = %d %d\n", i,*nr ,*nc);
}


void brnormed() {
  BUG;
}

int check_brnormed(model *cov) {
  model  
    *key = cov->key,
    *next = cov->sub[0],
    *sub = key != NULL ? key : next;

  ASSERT_ONESYSTEM;
  ASSERT_CARTESIAN;

  //  printf("start\n");
  kdefault(cov, NORMED_REJECTION, true);
  kdefault(cov, NORMED_OPTIMP, false);
  kdefault(cov, NORMED_NTH, NA_INTEGER); // i.e. adaptive
  kdefault(cov, NORMED_BURNIN, NA_INTEGER);// i.e. adaptive
  int
    err = NOERROR,
    total = Gettotalpoints(cov);

  if (total <= 1)
    SERR1("'%.50s' only works with at least 2 locations.", NICK(cov));
  if (!PisNULL(NORMED_PROB)) SERR1("'%.50s' must be given.", KNAME(NORMED_PROB));
  if (NROW(NORMED_PROB) !=1 && NROW(NORMED_PROB) != total)
    SERR1("length of '%.50s' must equal either 1 or the number of locations",
	  KNAME(NORMED_PROB));
 
  
  if (NROW(NORMED_PROB) == total) {
    double sum = 0.0;
    double *p = P(NORMED_PROB);
    for (int k=1; k<total; sum += p[k++]);
    sum = 1.0 - sum;
    if ((!ISNA(p[0]) && ABS(sum - p[0]) > 1e-5) || sum < 0.0)
      SERR1("Values of '%.50s' do not sum up to 1.", KNAME(NORMED_PROB));
    p[0] = sum;
  }
  
  // printf("normed null%d\n",PisNULL(NORMED_PROB));

   
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);

  Types type = isProcess(sub) ? SYSTYPE(DEFSYS(sub), 0) : VariogramType,
    frame = isVariogram(type) ? EvaluationType : BrMethodType;
 
  isotropy_type isoprev = isVariogram(type) ? SYMMETRIC : CARTESIAN_COORD;

  assert( OWNXDIM(0) *  OWNLOGDIM(0) ==1);
  cov->mpp.maxheights[0] = 1.0;
  
  if ((err = CHECK(sub, OWNLOGDIM(0), OWNXDIM(0),
		   type,//Martin: changed 27.11.13 DefList[COVNR].Type,
		   XONLY, isoprev, SCALAR, frame)) != NOERROR) {
    RETURN_ERR(err);
  }

  setbackward(cov, sub);

  //  printf("start ende\n");
  RETURN_ERR(err);
}


int struct_brnormed(model *cov, model **newmodel) {
  model *next = cov->sub[0];
  ASSERT_NEWMODEL_NULL;
  if (!hasSchlatherFrame(cov)) SERR2("'%.50s' may only called by '%.50s'",
				     NAME(cov), DefList[SCHLATHERPROC].name);

  if (GetTime(cov) || (Getgrid(cov) && Getcaniso(cov) != NULL)) {
    TransformLoc(cov, false, GRIDEXPAND_AVOID, false); // changes loc !
    SetLoc2NewLoc(next, PLoc(cov));
  }


  location_type *loc = Loc(cov);
  int i, d, err,
    centreidx[MAXMPPDIM],
    totalpoints = loc->totalpoints,
    newxlen = loc->grid ? 3 : totalpoints,
    dim = ANYDIMOF(next);
  bool grid = Getgrid(cov);
  double centreloc[MAXMPPDIM], minloc[MAXMPPDIM], maxloc[MAXMPPDIM],
    *newx= NULL, 
    **xgr = Getxgr(cov);

  if (cov->key != NULL) COV_DELETE(&(cov->key), cov); // should be ok
  NEW_STORAGE_WITH_SAVE(br);
  br_storage *sBR = cov->Sbr;

  assert(Getcaniso(cov) == NULL);
  GetDiameter(loc, minloc, maxloc, centreloc, false, true, centreidx);
  if ((newx = (double*) MALLOC(dim * newxlen * sizeof(double))) == NULL) {
     GERR("Memory allocation failed.\n"); 
  } 
  if (grid) {
    for (d=0; d<dim; d++) {
      newx[3 * d + XSTART] = xgr[d][XSTART] - centreloc[d];
      newx[3 * d + XSTEP] = xgr[d][XSTEP];
      newx[3 * d + XLENGTH] = xgr[d][XLENGTH];
    }
  } else {
    int endfor = totalpoints * dim;
    for (i=0; i<endfor; i+=dim)
      for (d=0; d<dim; d++) newx[i + d] = loc->x[i + d] - centreloc[d];
  }

  if ((err = loc_set(newx, NULL, dim, dim, newxlen, false, grid, // loc->grid
                     loc->distances,  cov)) > NOERROR)
    goto ErrorHandling;
  SetLoc2NewLoc(next, PLoc(cov));

  
  if ((err=covcpy(&(cov->key), next)) > NOERROR) goto ErrorHandling;
  if ((err = newmodel_covcpy(&(sBR->vario), VARIOGRAM_CALL, cov->key))!=NOERROR)
    goto ErrorHandling;
  if ((err = alloc_cov(sBR->vario, dim, 1, 1)) != NOERROR) goto ErrorHandling;
  if (isnowVariogram(next)) addModel(&(cov->key), GAUSSPROC, cov);
  assert(cov->key->calling == cov);

  //  PMI(cov);

  if ((err = CHECK(sBR->vario->sub[0], 1, 1, VariogramType, XONLY, 
		   ISOTROPIC, 1, EvaluationType)) != NOERROR) goto ErrorHandling;

  if ((err = CHECK_PASSTF(cov->key, ProcessType, cov->vdim[0],
			  GaussMethodType)) != NOERROR)
    goto ErrorHandling;  
  if ((err = STRUCT(cov->key, NULL)) > NOERROR) goto ErrorHandling;

  
  ErrorHandling:
  FREE(newx);
 
 RETURN_ERR(err);
  
}


double* getCi(model *cov, int i) {
  br_storage *sBR = cov->Sbr;
  if (sBR->normed.C[i] != NULL) return sBR->normed.C[i];
  double **Ci = &(sBR->normed.dummyCi);
  if (sBR->normed.nCis < sBR->normed.maxCi) {
    Ci = sBR->normed.C + i;
    sBR->normed.nCis++;
  }
  
  bool null;
  if ((null = *Ci == NULL))
    *Ci = (double*) MALLOC(sBR->normed.total * sizeof(double));

  //  printf("%d\n", null);
  
  if (null ||
      (sBR->normed.nCis >= sBR->normed.maxCi && i != sBR->normed.current_i)) {
    CovarianceMatrixCol(sBR->vario->sub[0], i, *Ci);
    sBR->normed.current_i = i;
  }
  return *Ci;
}



void NormedSimulation(model *cov, gen_storage *s) {
  br_storage *sBR = cov->Sbr;
  model *key = cov->key;
  double *field = key->rf,
    *res = cov->rf,
    *p = P(NORMED_PROB),
    *trend = sBR->trend[ORIG_IDX];
  pgs_storage *pgs = cov->Spgs;
 
  if (P0INT(NORMED_REJECTION)) {
    BUG;
  } else { // mcmc
    int
      total = sBR->normed.total,
      endfor = sBR->normed.nth,
      zeropos = sBR->zeropos;
    for (int n=0; n<endfor; n++) {
      sBR->normed.zaehler++;
      double u = UNIFORM_RANDOM;
      int koben,
	i = sBR->normed.total / 2;  // kunten
      while (p[i] >= u && i > 0) i /= 2;
      koben = i * 2 + 1;
      if (koben >= n) koben = n - 1;
      while (i <= koben) {
	int kneu = (koben + i) / 2;
	if (p[kneu] >= u) koben = kneu;
	else i = kneu + 1;
      }
      double
	sum = 0.0,
	max = RF_NEGINF,
	*Ci = getCi(cov, i);
      DO(key, s); // nicht gatternr,
      double fieldzeropos = (double) field[zeropos];// wird durch DO veraendert!
      for (int k=0; k<total; k++) {
	//	printf("%d %d %d %10g\n", k, i, n, field[k] );
	//	printf("%10g\n", Ci[k]);
	//	printf("%10g\n", fieldzeropos );
	//	printf("%10g\n", trend[k] );
	double
	  dummy = field[k] = EXP(field[k] + Ci[k] - fieldzeropos - trend[k]);
	if (dummy > max) max = dummy;
	sum += p[k] * dummy;
      }
      // eigentlich nur das Feld nach den endfor Simulationen verwenden,
      // da innerhalb abhaengig. Deshalb n_zhou nur um eins erhoeht und
      // der Mittelwert ueber alle endfor Fehlder zu sum_zhou addiert.
      pgs->sum_zhou_c += max / (double) endfor;
      double fmaxDfprop = max / sum,
	ratio = fmaxDfprop / sBR->normed.fmaxDfprop;
      if (ratio >= 1.0 || UNIFORM_RANDOM < ratio) {
	for (int k=0; k<total; k++) res[k] = field[k] / max;
	sBR->normed.fmaxDfprop = fmaxDfprop;
	sBR->normed.max = max;
	sBR->normed.accepted++;
      }
    } // end for
  }
  pgs->n_zhou_c ++;
}
 

// total, nth, zeropos, zaehler, fmaxDfprop, max, accepted, C, dummyCi, maxCi

int init_brnormed(model *cov, gen_storage *s){

  br_storage *sBR = cov->Sbr;
  pgs_storage *pgs = NULL; 
  int err,
    total = Gettotalpoints(cov),    
    dim = ANYDIM,
    burnin = P0INT(NORMED_BURNIN),
    nth = P0INT(NORMED_NTH),
    zeropos = sBR->zeropos;
  double *trend = sBR->trend[ORIG_IDX];
  if ((err = alloc_cov(cov, dim, 1, 1)) != NOERROR) RETURN_ERR(err);
  pgs = cov->Spgs;
  if ((err = general_init(cov, 1, s)) != NOERROR) RETURN_ERR(err);
  Variogram(NULL, sBR->vario, trend);

  sBR->normed.total = Gettotalpoints(cov);
  assert(sBR->normed.current_prob == NULL);
  sBR->normed.current_prob =
    (double*) MALLOC(sBR->normed.total * sizeof(double));
  assert(sBR->normed.C == NULL);
  sBR->normed.C = (double**) CALLOC(sBR->normed.total, sizeof(double*));
  //  assert(sBR->normed.field == NULL);
  //  sBR->normed.field = (double*) MALLOC(sBR->normed.total * sizeof(double));
  assert(sBR->normed.dummyCi == NULL);

  if (NROW(NORMED_PROB) == 1) {
    int start = (int) P0(NORMED_PROB);
    PFREE(NORMED_PROB);
    PALLOC(NORMED_PROB, total, 1);
     if (start == BRP_COV && total > 20000) start = BRP_UNIF;
    else if (start == BRP_KRIG && total > 1000) start = BRP_UNIF;
     if (start != P0INT(NORMED_PROB) && PL >=  PL_IMPORTANT ) {
	PRINTF("value of '%.50s' has changed from %d to %d",
	       KNAME(NORMED_PROB), (int) P0(NORMED_PROB), start);
     }
    if (start == BRP_UNIF) {
      double unif = 1.0 / (double) total;
      for (int i=0; i<total; P(NORMED_PROB)[i++] = unif);
    } else {
      assert(start == BRP_KRIG || start == BRP_COV);
      int totalSq = total * total;
      sBR->normed.do_not_delete_C = true;     
      double
	*p = P(NORMED_PROB),
	*c  =  sBR->normed.C[0] = (double*) MALLOC(totalSq * sizeof(double));
      
      double *one = (double *) MALLOC(total * sizeof(double));   
      one[0] = 1.0;
      for (int i=1; i<total; i++) {
	sBR->normed.C[i] = c + i * total;
	one[i] = 1.0;
      }

      if (start == BRP_COV) {
	CovarianceMatrix(sBR->vario->sub[0], c);
	for (int i=0; i<totalSq; i++) c[i] = EXP(c[i]);
      } else { // start == BRP_KRIG 
#define nSigma 100
	model *key = cov->key;
	if (key == NULL) BUG;
	double *field = key->rf,
	  *neu  = (double*) CALLOC(totalSq, sizeof(double));
	for (int m=0; m<nSigma; m++) {
	  DO(key, s); // nicht gatternr,
	  double fieldzeropos = (double) field[zeropos];// d. DO veraendert!
	  for (int k=0; k<total; k++) field[k] -= fieldzeropos + trend[k];
	  for(int i=0; i<total; i++) {
	    double
	      *sigma_ik = c + i,
	      *sigma_ki = c + i * total,
	      max = RF_NEGINF,
	      *Ci = getCi(cov, i);
	    for (int k=0; k<total; k++) {
	      double dummy = neu[k] = field[k] + Ci[k];
	      if (dummy > max) max = dummy;
	    }
	    // untere dreiecksmatrix	    
	    for (int k=0; k<=i; k++) sigma_ik[k * total] += EXP(neu[k] - max);
	    for (int k=i+1; k<total; k++) sigma_ki[k] += EXP(neu[k] - max);
	  }
	}
	double minv = 1.0 / (double) nSigma;
	for (int i=0; i<total; i++) {
	  for (int k=0; k<=i; k++) {	    
	    c[k + i * total] = (c[i + k * total] *= minv);
	  }
	}
	FREE(neu);
      }
      err = Ext_solvePosDef(c, total, true, one, 1, p, NULL);
      FREE(one);
      if (err != NOERROR) RETURN_ERR(err);
      double sum = 0.0;
      for (int i=0; i<total; i++) if (p[i] > 0.0) sum += p[i]; else p[i] = 0.0;
      for (int i=0; i<total; i++) p[i] /= sum;
    }
  }
 
  pgs->zhou_c = 1.0;
  pgs->n_zhou_c = 0;
  pgs->sum_zhou_c = 0.0;
  pgs->estimated_zhou_c = true;
  sBR->normed.maxCi = (int) CEIL((double) GLOBAL.br.BRmaxmem / dim);  
   
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0; // 1.0 is import dummy value 
  cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1.0;
  cov->mpp.maxheights[0] = 1.0;

  if ((err = ReturnOwnField(cov)) != NOERROR)  // must be later than INIT !
    RETURN_ERR(err);
  
#define min_nth 2
#define burninDnth 50
#define burnin_min_rejections 5
#define factor 3

  sBR->normed.zaehler = 0;
  sBR->normed.accepted = 0;
  sBR->normed.nth = 1;
  sBR->normed.current_i = UNSET;
  if (P0(NORMED_REJECTION)) {
    
  } else {
    if (burnin == NA_INTEGER) { // adaptive burnin
      burnin = min_nth * burninDnth;
      int zaehler = 0;
      while (zaehler < burnin) {
	// (1 - q - k sigma) n =  burnin_min_rejections
	// with q the estimated acceptance rate
	//sigma^2 = q(1-q) / m  variance of q based on the simulation up to now	
	// and k an arbitrary factor, at least or about 2

	NormedSimulation(cov, s);
	
	if (++zaehler == burnin) {
	  double
	    q = sBR->normed.accepted / sBR->normed.zaehler,      
	    N = burnin_min_rejections /
	    (q - factor * SQRT(q * (1 - q) / sBR->normed.zaehler));
	  if (!R_finite(N) || N < 0) burnin *= 2;
	  else if (N > burnin) burnin = N;
	  else break;
	}
      }
    } else {// fixed burnin
      for (int zaehler=0; zaehler<burnin; zaehler++) NormedSimulation(cov, s);
    }

    if ((sBR->normed.adaptive_nth = nth == NA_INTEGER)) {
      nth = ROUND((double) burnin / burninDnth);
      if (nth <= 0)
	SERR1("'%.50s' cannot be determined without burnin", KNAME(NORMED_NTH));
      if (nth > 10000)
	WARN2("'%.50s'=%d is very large", KNAME(NORMED_NTH), nth);
    }
    sBR->normed.burnin = burnin;
    sBR->normed.nth = nth;
  }

  double *p = P(NORMED_PROB);
  int totalM1 = sBR->normed.total - 1;
  sBR->normed.current_cumprob[0] = p[0];
  for (int k=1; k<totalM1; k++)
    sBR->normed.current_cumprob[k] = sBR->normed.current_cumprob[k-1] + p[k];
  sBR->normed.current_cumprob[totalM1] = 1.0;

  cov->simu.active = cov->initialised = err == NOERROR;
  RETURN_ERR(err);
}


#define NORMED_EACH 100
void do_brnormed(model *cov, gen_storage *s) {
  assert(cov->key != NULL);
  br_storage *sBR = cov->Sbr;
  assert(sBR != NULL && sBR->trend != NULL);
  //  double *res = cov->rf,
  //    *trend = sBR->trend[ORIG_IDX];
  unsigned int each = NORMED_EACH * sBR->normed.nth;
  assert(sBR->normed.total > 0);
  assert(sBR->zeropos >= 0 && sBR->zeropos < sBR->normed.total);
  assert(cov->origrf);

  NormedSimulation(cov, s);

  if (sBR->normed.zaehler % each == 0) {
    if (P0INT(NORMED_OPTIMP) && true) {
      // update the probabilities
      BUG;
      // sBR->normed.current_prob
    }
    if (sBR->normed.adaptive_nth) {      
      double q = sBR->normed.accepted / sBR->normed.zaehler,
	N = burnin_min_rejections /
  	    (q - factor * SQRT(q * (1 - q) / sBR->normed.zaehler)); 
      sBR->normed.nth = ROUND(N / burninDnth);
    }
  }
  
  
}
 

void range_brnormed(model VARIABLE_IS_NOT_USED *cov,
    range_type *range){  
  range->min[NORMED_PROB] = 0;
  range->max[NORMED_PROB] = 1;
  range->pmin[NORMED_PROB] = 0;
  range->pmax[NORMED_PROB] = 1;
  range->openmin[NORMED_PROB] = true;
  range->openmax[NORMED_PROB] = true;
  
  booleanRange(NORMED_OPTIMP);
  
  range->min[NORMED_NTH] = 1;
  range->max[NORMED_NTH] = RF_INF;
  range->pmin[NORMED_NTH] = 1;
  range->pmax[NORMED_NTH] = 1e4;
  range->openmin[NORMED_NTH] = false;
  range->openmax[NORMED_NTH] = true;
  
  range->min[NORMED_BURNIN] = 0;
  range->max[NORMED_BURNIN] = RF_INF;
  range->pmin[NORMED_BURNIN] = 1;
  range->pmax[NORMED_BURNIN] = 1e6;
  range->openmin[NORMED_BURNIN] = false;
  range->openmax[NORMED_BURNIN] = true;
    
  booleanRange(NORMED_REJECTION);
}
