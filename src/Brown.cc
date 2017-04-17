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

#include <Rmath.h>
#include <stdio.h>
#include "RF.h"
#include "Operator.h"

#define BR_MESHSIZE (LAST_MAXSTABLE + 1)
#define BR_VERTNUMBER (LAST_MAXSTABLE + 2)
#define BR_OPTIM (LAST_MAXSTABLE + 3)
#define BR_OPTIMTOL (LAST_MAXSTABLE + 4)
#define BR_LAMBDA (LAST_MAXSTABLE + 5)
#define BR_OPTIMAREA (LAST_MAXSTABLE + 6)
#define BR_VARIOBOUND (LAST_MAXSTABLE + 7)


// **********************************************************************
// Brown Resnick

int checkBrownResnickProc(cov_model *cov) {

  //NotProgrammedYet("at the moment Brown-Resnick processes");

  cov_model  
    *key = cov->key,
    *sub = key != NULL ? key :
            cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  int err, role, 
      dim = cov->tsdim;
  isotropy_type isoprev;
  Types type;
  
  ASSERT_CARTESIAN;

  ASSERT_ONE_SUBMODEL(cov);
  
  if ((err = SetGEVetc(cov, ROLE_BROWNRESNICK)) != NOERROR) return err;
  
  role = isVariogram(sub) ? ROLE_COV
    : (isGaussProcess(sub) &&  isPointShape(cov)) ? ROLE_GAUSS
    : isBrownResnickProcess(sub) ? ROLE_BROWNRESNICK
    : isPointShape(sub) ? ROLE_BROWNRESNICK
    : ROLE_UNDEFINED;    

  type = isProcess(sub) || isPointShape(sub) 
    ? CovList[sub->nr].Typi[0] : VariogramType;
    
  ASSERT_ROLE_DEFINED(sub);  
  
  isoprev = role == ROLE_COV ? SYMMETRIC : CARTESIAN_COORD;
  
  if ((err = CHECK(sub, dim, dim, 
		   type,//Martin: changed 27.11.13 CovList[cov->nr].Type,
		   XONLY, isoprev, 1, role)) != NOERROR) {
      return err;
  }
  
  setbackward(cov, sub);
  if (cov->vdim[0] != 1) SERR("BR only works in the univariate case");  
  
  return NOERROR;  
  isoprev = role == ROLE_COV ? SYMMETRIC : CARTESIAN_COORD;
  
  if ((err = CHECK(sub, dim, dim, 
		   type,//Martin: changed 27.11.13 CovList[cov->nr].Type,
		   XONLY, isoprev, 1, role)) != NOERROR) {
      return err;
  }
  
  setbackward(cov, sub);
  if (cov->vdim[0] != 1) SERR("BR only works in the univariate case");  
  
  return NOERROR;  
}

#define ORIG_IDX 0
int init_BRorig(cov_model *cov, gen_storage *s){
  cov_model *key = cov->key;
  location_type *keyloc = NULL;
  int err, d, 
    dim = cov->tsdim;
  bool keygrid;
  br_storage *sBR = NULL;

  if (cov->role != ROLE_BROWNRESNICK) ILLEGAL_ROLE;    
  if (key == NULL) BUG;
  if ((err = alloc_cov(cov, dim, 1, 1)) != NOERROR) return err;
  pgs_storage *pgs = cov->Spgs;

  for (d=0; d<dim; d++) {
    pgs->supportmin[d] = RF_NEGINF; // 4 * for debugging...
    pgs->supportmax[d] = RF_INF;
    pgs->supportcentre[d] = RF_NA;      
  }
  
  pgs->log_density = 0.0;
  
  keyloc = Loc(key);
  keygrid = keyloc->grid;

  key->simu.active = true;
  key->simu.expected_number_simu = cov->simu.expected_number_simu;
  
  assert(cov->mpp.moments >= 1);
  if ((err = INIT(key, 1, s)) != NOERROR) goto ErrorHandling;  
  
  cov->loggiven = true;
  
  
  assert(key->nr == GAUSSPROC);
  assert(key->mpp.moments >= 1);
  
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
  cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1.0;
  cov->mpp.maxheights[0] = EXP(GLOBAL.extreme.standardmax);

  pgs->zhou_c = 1.0;
  
  sBR = cov->Sbr;
  sBR->trendlen = 1;
  if ((sBR->trend = (double**) MALLOC(sizeof(double*)))==NULL) {
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
  }    
  if ((sBR->trend[ORIG_IDX] = 
       (double*) MALLOC(keyloc->totalpoints*sizeof(double)) ) == NULL) {
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  
  if ((err = loc_set(keygrid ? keyloc->xgr[0] : keyloc->x, NULL, NULL, dim,
	             dim, keygrid ? 3 : keyloc->totalpoints, 0, false, keygrid,
		     keyloc->distances, cov->Sbr->vario)) > NOERROR)
    goto ErrorHandling;
  if (sBR->vario->sub[0] != NULL) 
    SetLoc2NewLoc(sBR->vario->sub[0], PLoc(sBR->vario));
  Variogram(NULL, sBR->vario, sBR->trend[ORIG_IDX]);
  
  if ((err = FieldReturn(cov)) != NOERROR)  // must be later than INIT !
   goto ErrorHandling;  

 ErrorHandling:
  if (err != NOERROR) br_DELETE(&(cov->Sbr));
  return err;
}

void do_BRorig(cov_model *cov, gen_storage *s) {
  cov_model *key = cov->key;
  assert(key != NULL);
  br_storage *sBR = cov->Sbr;
  assert(sBR != NULL && sBR->trend != NULL);
  double *res = cov->rf;
#define ORIG_IDX 0
  int i,
    zeropos = sBR->zeropos;
  long
    totalpoints = Loc(cov)->totalpoints;
  assert(totalpoints > 0);
  assert(zeropos >= 0 && zeropos <totalpoints);
  double *trend = sBR->trend[ORIG_IDX];
  assert(cov->origrf);

  DO(key, s); // nicht gatternr,
  double *lgres = key->rf;
  double lgreszeropos = (double) lgres[zeropos]; // wird durch DO veraendert!
  for (i=0; i<totalpoints; i++) {
    res[i] = lgres[i] - lgreszeropos - trend[i];
  }
}

int init_BRshifted(cov_model *cov, gen_storage *s) {
  cov_model *key=cov->key;
  location_type *keyloc = NULL;
  int d, dim,
    err = NOERROR;
  long j, shiftedloclen, keytotal, trendlenmax, trendlenneeded;
  bool keygrid;
  br_storage *sBR = NULL;
  
  if (cov->role == ROLE_BROWNRESNICK) {
    
    if (key != NULL) {
      dim = cov->tsdim; 
      if ((err = alloc_cov(cov, dim, 1, 1)) != NOERROR) return err;
      pgs_storage *pgs = cov->Spgs;
      for (d=0; d<dim; d++) {
	pgs->supportmin[d] = RF_NEGINF; // 4 * for debugging;
	pgs->supportmax[d] = RF_INF;
	pgs->supportcentre[d] = RF_NA;
      }
      pgs->log_density = 0.0;
      
      keyloc = Loc(key);
      keygrid = keyloc->grid;
      keytotal = keyloc->totalpoints;
      
      key->simu.active = true;
      key->simu.expected_number_simu = cov->simu.expected_number_simu; 
      assert(cov->mpp.moments >= 1);
      if ((err = INIT(key, 1, s)) != NOERROR) return err;  
      
      cov->loggiven = true;
      
      assert(key->nr == GAUSSPROC);
      assert(key->mpp.moments >= 1);
      cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
      cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1.0;
      cov->mpp.maxheights[0] = EXP(GLOBAL.extreme.standardmax);
      pgs->zhou_c = 1.0;
      
      sBR = cov->Sbr;
       
      shiftedloclen = keygrid ? 3 : keytotal;
      if ((sBR->shiftedloc = (double*)
	   MALLOC(dim*shiftedloclen*sizeof(double))) == NULL ||
	  (sBR->locindex = (int*) MALLOC(sizeof(int) * dim))==NULL
	  ) {
        err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      
      trendlenmax = (int) CEIL((double) GLOBAL.br.BRmaxmem / keytotal);
      trendlenneeded = MIN(keytotal, cov->simu.expected_number_simu);
      sBR->trendlen = MIN(trendlenmax, trendlenneeded);
   
      assert(sBR->trendlen > 0);
      
      sBR->memcounter = 0;
      if ((sBR->loc2mem=(int*) MALLOC(sizeof(int)*keytotal))==NULL) {
        err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      for (j=0; j<keytotal; j++) sBR->loc2mem[j] = -1;    
    
      if ((sBR->mem2loc=(int*) MALLOC(sizeof(int)*sBR->trendlen))==NULL) {
        err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      if ((sBR->trend=(double**) MALLOC(sizeof(double*)*sBR->trendlen))
           ==NULL) {
       err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      for (j=0; j<sBR->trendlen; j++) {
        sBR->mem2loc[j] = -1;
        if ((sBR->trend[j] = 
           (double*) MALLOC(keytotal*sizeof(double)))==NULL) {
          err = ERRORMEMORYALLOCATION; goto ErrorHandling;
        }
      }
      
      if ((err = loc_set(keygrid ? keyloc->xgr[0] : keyloc->x, NULL, NULL, dim,
	                 dim, keygrid ? 3 : keytotal, 0, false, keygrid,
	                 keyloc->distances, sBR->vario)) > NOERROR)
        return err;
      if (sBR->vario->sub[0] != NULL) 
        SetLoc2NewLoc(sBR->vario->sub[0], PLoc(sBR->vario));   
      
      if ((err = FieldReturn(cov)) != NOERROR)  // must be later than INIT !
        return err;
    }
        
    goto ErrorHandling; // no error
    
  }
  
  else ILLEGAL_ROLE;
  
  ErrorHandling:
    if (err != NOERROR) br_DELETE(&(cov->Sbr));
    return err;

}

void indextrafo(long onedimindex, double ** xgr, int dim, int *multidimindex) {
  int d;
  for (d=0; d<dim; d++) {
    multidimindex[d] = onedimindex % (int) xgr[d][XLENGTH];
    onedimindex = onedimindex / (int) xgr[d][XLENGTH];    
  }
}

void do_BRshifted(cov_model *cov, gen_storage *s) {
  br_storage *sBR = cov->Sbr;
  cov_model *key = cov->key;
  assert(cov->key != NULL);
  
  location_type *keyloc = Loc(key);
  long i, k, zeropos, zeroposMdim, keytotal = keyloc->totalpoints;
  int d,  trendindex,
    dim = cov->tsdim,
    trendlen = sBR->trendlen,
    *locindex = sBR->locindex,
    *mem2loc = sBR->mem2loc,
    *loc2mem = sBR->loc2mem;
  bool keygrid = keyloc->grid;
  double *shiftedloc = sBR->shiftedloc,
    **xgr = keyloc->xgr,
    **trend = sBR->trend;
  assert(cov->origrf);
  double *res = cov->rf;
  double *lgres = cov->key->rf;

  DO(key, s);
  zeropos = (long) FLOOR(UNIFORM_RANDOM * keytotal);
  
  if (loc2mem[zeropos] > -1) {
    trendindex = loc2mem[zeropos];
    if (mem2loc[trendindex] != zeropos) BUG;
  } else {
    if (sBR->memcounter<trendlen) {
      trendindex = sBR->memcounter; 
      sBR->memcounter++;
    } else {
      trendindex = trendlen - 1; 
      loc2mem[mem2loc[trendlen-1]] = -1;
      mem2loc[trendlen-1] = -1;
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

int check_BRmixed(cov_model *cov) {
  //NotProgrammedYet("at the moment Brown-Resnick processes");
  
  int err;
  br_param *bp = &(GLOBAL.br);
   
  ASSERT_ONE_SUBMODEL(cov);
  
  if (!cov->logspeed) SERR("BrownResnick requires a variogram model as submodel that tends to infinity [t rate of at least 4log(h) for being compatible with BRmixed");
  
  kdefault(cov, BR_MESHSIZE, bp->BRmeshsize);
  kdefault(cov, BR_VERTNUMBER, bp->BRvertnumber);
  kdefault(cov, BR_OPTIM, bp->BRoptim);
  kdefault(cov, BR_OPTIMTOL, bp->BRoptimtol);
  kdefault(cov, BR_VARIOBOUND, bp->variobound);
  if (cov->nr == BRMIXED_USER && cov->key == NULL && P0INT(BR_OPTIM) > 0) {
    if (!PisNULL(BR_LAMBDA)) {
      if (PisNULL(BR_OPTIMAREA)) SERR1("'%s' not given", KNAME(BR_OPTIMAREA));
      if (PL > 0) PRINTF("'%s' set to '0'", KNAME(BR_OPTIM));
      PINT(BR_OPTIM)[0] = 0;
    } else if (P0INT(BR_OPTIM) == 2 && !PisNULL(BR_OPTIMAREA))
      if (PL > 0) PRINTF("'%s' set to '1'", KNAME(BR_OPTIM));
  }
  
  if (cov->key != NULL && P0INT(BR_OPTIM) == 2) {
    if (!isIsotropic(cov->key->isoown)) {
     // SERR("area optimisation implemented for the isotropic case only"); //@MARTIN: das scheint nicht zu funktionieren, wenn ich ein Variogramm eingebe
    }  
  }    

  kdefault(cov, BR_LAMBDA, RF_NA);
  if (PisNULL(BR_OPTIMAREA))
    kdefault(cov, BR_OPTIMAREA, 0.0);
  
  if ((err = checkBrownResnickProc(cov)) != NOERROR) return err;
  
  if ((err = checkkappas(cov, true)) != NOERROR) return err;
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
    return ERRORDIM;
  if (cov->vdim[0] != 1) SERR("BR only works in the univariate case");

  return NOERROR; 
}

void kappaBRmixed(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
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
      *nr = -1; 
      *nc = -1;
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

void OptimArea(cov_model *cov) {
  // side effect auf P(BR_OPTIMAREA) !!
    
  br_storage *sBR = cov->Sbr;
  pgs_storage *pgs = cov->Spgs;
  cov_model *key = sBR->sub[0];
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
    minradius = (int) (sBR->minradius/step),
    **countvector = sBR->countvector,
    zeropos = sBR->zeropos,
    dim = cov->tsdim,
    keytotal = keyloc->totalpoints;

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
    //printf("%f ", optimarea[d]);
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


void set_lowerbounds(cov_model *cov) {    
  br_storage *sBR = cov->Sbr;
  assert(sBR != NULL);
  assert(sBR->sub[0] != NULL);
  double step = P0(BR_MESHSIZE),
    *optimarea = P(BR_OPTIMAREA);    
  int j, k, dim = cov->tsdim,
    minradius = (int) (sBR->minradius / step);
  cov_model *key = sBR->sub[0];
  location_type *keyloc = Loc(key);
  double **xgr = keyloc->xgr;
  long keytotal = keyloc->totalpoints;
 
  for (j=0; j<keytotal; j++) {
    sBR->lowerbounds[j] = RF_INF;
    k = (int) CEIL(IdxDistance(j, sBR->zeropos, xgr, dim));
    if (k <= minradius && optimarea[k]>1e-5) {
      sBR->lowerbounds[j] = -LOG(optimarea[k]);  
    }  
    //printf("%f ", sBR->lowerbounds[j]);
  }
  //printf("\n");
}

int prepareBRoptim(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  br_storage *sBR = cov->Sbr;
  cov_model *key = sBR->sub[0];
  location_type *keyloc = Loc(key);
  double  step = P0(BR_MESHSIZE),
        **xgr = keyloc->xgr;
  int i, j, d,
    vertnumber = P0INT(BR_VERTNUMBER),
    dim = cov->tsdim,
    maxradius = 1,
    minradius = (int) (sBR->minradius/step);

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
    sBR->vertnumber = vertnumber; 
    
    if (sBR->countvector != NULL || sBR->areamatrix != NULL) BUG;
    if ((sBR->countvector = (int**) CALLOC(vertnumber, sizeof(int*)))==NULL || 
        (sBR->logvertnumber = (double *) MALLOC(vertnumber * sizeof(double)))
         == NULL)
      return ERRORMEMORYALLOCATION;
    for (j=0; j<vertnumber; j++) {
      if ((sBR->countvector[j] = (int*) CALLOC(minradius + 1, sizeof(int))) 
          == NULL)
        return ERRORMEMORYALLOCATION;
      for (i=0; i<=minradius; i++) sBR->countvector[j][i] = 0;
    }
    for (j=0; j<vertnumber; j++)
      sBR->logvertnumber[j] = - LOG((double) (j+1)/vertnumber);    
    break;

  default:
    SERR("optimization might not be used here\n");
  }
  
  if ((sBR->areamatrix  = (double *) MALLOC((minradius + 1)* sizeof(double))) 
         == NULL) {
    return ERRORMEMORYALLOCATION;
  }
  sBR->areamatrix[0] = 1.0;
  if (minradius > 0) {
    for (i=1; i<=minradius; i++) {
      if (i <= cov->ncol[BR_OPTIMAREA]) {
        sBR->areamatrix[i] = P(BR_OPTIMAREA)[i-1];  
      }  
      else sBR->areamatrix[i] = 0.0;
    }
  }
  
  PFREE(BR_OPTIMAREA);
  PALLOC(BR_OPTIMAREA, 1, minradius + 1);
  double *optimarea = P(BR_OPTIMAREA);
  for (i=0; i<=minradius; i++) {
    optimarea[i] = sBR->areamatrix[i];
  }
  set_lowerbounds(cov);
  
  if (PL >= PL_STRUCTURE) PRINTF("BR optimisation finished...\n");

 return NOERROR;	
}



int init_BRmixed(cov_model *cov, gen_storage *s) {
  location_type *loc = Loc(cov);
  br_storage *sBR = cov->Sbr;
  assert(sBR != NULL);
  assert(sBR->sub[0] != NULL);
  cov_model *key = sBR->sub[0];
  pgs_storage *pgs = NULL;
  location_type *keyloc = Loc(key);
  int  d, err = NOERROR, 
    dim = cov->tsdim,
    bytes = sizeof(double) * dim,
    keytotal = keyloc->totalpoints;
  double area = 1.0,
    step = P0(BR_MESHSIZE);
    
  sBR->trendlen = 1;
  assert(isPointShape(cov));
  if (cov->role != ROLE_BROWNRESNICK) ILLEGAL_ROLE;
  assert(dim > 0);
  
  if ((err = alloc_cov(cov, dim, 1, 1)) != NOERROR) goto ErrorHandling;; // nur pgs
  pgs = cov->Spgs; // nach alloc_cov !!

  if ((sBR->suppmin = (double*) MALLOC(bytes))==NULL ||
      (sBR->suppmax = (double*) MALLOC(bytes))==NULL ||
      (sBR->locmin = (double*) MALLOC(bytes))==NULL ||
      (sBR->locmax = (double*) MALLOC(bytes))==NULL ||
      (sBR->loccentre = (double*) MALLOC(bytes))==NULL ||
      (sBR->locindex = (int*) MALLOC(sizeof(int) * dim))==NULL ||
      (sBR->trend == NULL && 
      (sBR->trend = (double**) CALLOC(sBR->trendlen, sizeof(double*)))==NULL) 
     ) {
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
 
  GetDiameter(loc, sBR->locmin, sBR->locmax, sBR->loccentre);
  assert(pgs->supportmin != NULL);
  for (d=0; d<dim; d++) {
    pgs->own_grid_cumsum[d] = d==0 ? 1 : pgs->own_grid_cumsum[d-1]*pgs->own_grid_len[d-1];      
    sBR->suppmin[d] =  FLOOR((sBR->locmin[d] - sBR->radius - sBR->minradius)/step)*step - step/2;  
    sBR->suppmax[d] =  CEIL((sBR->locmax[d] + sBR->radius + sBR->minradius)/step)*step + step/2;  
    area *= sBR->suppmax[d] - sBR->suppmin[d];
    pgs->supportmin[d] = RF_NEGINF;
    pgs->supportmax[d] = RF_INF;
    pgs->own_grid_start[d] = RF_NEGINF;
    pgs->own_grid_step[d] = keyloc->xgr[d][XSTEP];
    pgs->own_grid_len[d] = keyloc->xgr[d][XLENGTH];
  }
  
  pgs->log_density = 0.0; //@MARTIN: besser: +/-LOG(area)? Rolle von logdens nicht ganz klar (s. extremes.cc)
  for (d=0; d<=cov->mpp.moments; d++) {
    cov->mpp.mM[d] = cov->mpp.mMplus[d] = 1.0;
  }
  cov->mpp.maxheights[0] = EXP(0.0);
    
  assert(keyloc->grid);
  assert(key->nr == GAUSSPROC);

  key->simu.expected_number_simu = cov->simu.expected_number_simu;     
  assert(cov->mpp.moments >= 1);
  if ((err = INIT(key, 1, s)) != NOERROR) goto ErrorHandling;
  key->simu.active = true;
  cov->loggiven = true;
    
  assert(key->mpp.moments >= 1);
  key->mpp.mM[0] = key->mpp.mMplus[0] = 1.0;
  key->mpp.mM[1] = key->mpp.mMplus[1] = 1.0;
    
  if (sBR->trend[0] == NULL && 
     (sBR->trend[0] =(double*) MALLOC(keytotal*sizeof(double)))==NULL){
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
    
  if ((err = loc_set(keyloc->xgr[0], NULL, NULL, dim, dim, 3, 0,
       false, true, keyloc->distances, sBR->vario)) > NOERROR) 
     goto ErrorHandling;
  
  if (sBR->vario->sub[0] != NULL)
    SetLoc2NewLoc(sBR->vario->sub[0], PLoc(sBR->vario));
  Variogram(NULL, sBR->vario, sBR->trend[0]);
  
  if ((sBR->lowerbounds = (double*) MALLOC(keytotal*sizeof(double))) == NULL) { 
    err=ERRORMEMORYALLOCATION; 
    goto ErrorHandling; 
  }
  
  if ((err = prepareBRoptim(cov, s)) != NOERROR) {
     goto ErrorHandling;
  }

  assert(keyloc != NULL);
  key->simu.active = true;

  set_lowerbounds(cov);   
  cov->rf = sBR->sub[0]->rf;
  cov->origrf = false;
  cov->fieldreturn = HUETCHEN_OWN_GRIDSIZE;
  pgs->estimated_zhou_c = ISNAN(P0(BR_LAMBDA));
  pgs->zhou_c =  pgs->estimated_zhou_c ? 0.0 : P0(BR_LAMBDA)*area; //@MARTIN: area kann weg, falls in logdens
  pgs->logmean = false;
  pgs->sq_zhou_c = pgs->sum_zhou_c = 0.0;
  pgs->n_zhou_c = 0;
  sBR->next_am_check = GLOBAL.br.deltaAM;

 
 ErrorHandling:
  if (err != NOERROR) br_DELETE(&(cov->Sbr));

  return err;

}

void do_BRmixed(cov_model *cov, gen_storage *s) {
  // to do: improve simulation speed by dynamic sizes
  assert(cov->key!=NULL);
  br_storage *sBR = cov->Sbr;
  cov_model  *key = sBR->sub[0];
  assert(cov->rf == key->rf);
  location_type *keyloc = Loc(key);
  assert(keyloc->grid); 
  pgs_storage *pgs = cov->Spgs;
  assert(pgs != NULL); 

  int d, dim = cov->tsdim, minradius,
      i, j, maxind, idxdist, hatnumber=0,
      lgtotalpoints = keyloc->totalpoints,
      zeropos = sBR->zeropos,
      vertnumber = P0INT(BR_VERTNUMBER);

  double step = P0(BR_MESHSIZE),
    invstepdim = intpow(step, -dim),
    uplusmaxval , maxval, u[MAXMPPDIM], 
    ** xgr = keyloc->xgr,
    *lowerbounds = sBR->lowerbounds,
    area=1.0, *lgres = key->rf, //@MARTIN: area obsolet, falls in logdens
    *trend = sBR->trend[0];
      
  if (P0INT(BR_OPTIM) == 2 &&  pgs->n_zhou_c >= sBR->next_am_check) {      
    sBR->next_am_check += GLOBAL.br.deltaAM; 
    OptimArea(cov); 
    set_lowerbounds(cov);
  }

  minradius = (int) (sBR->minradius/step); 
  
  for (d=0; d<dim; d++) {
    u[d] = ROUND((UNIFORM_RANDOM*(sBR->suppmax[d] - sBR->suppmin[d]) +
		  sBR->suppmin[d]) / step) * step;
    area *= sBR->suppmax[d] - sBR->suppmin[d];
    pgs->supportmin[d] = u[d] - sBR->minradius - sBR->radius; 
    pgs->supportmax[d] = u[d] + sBR->minradius + sBR->radius; 
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
        if (uplusmaxval > sBR->logvertnumber[j]) {
          idxdist = (int) CEIL(IdxDistance(maxind, zeropos, xgr, dim));
          if (idxdist<=minradius) (sBR->countvector[j][idxdist])++;
          break;
        }
      }
    }

    if (uplusmaxval > lowerbounds[maxind]) break;    
  }
  
  pgs->n_zhou_c += hatnumber;

  if (PL >= PL_STRUCTURE && hatnumber > 300)
    PRINTF("note: large hat number (%d) might indicate numerically suboptimal framework\n",
	   hatnumber);
  
  //shifting maximum to origin is not necessary because of stationarity
  //(conditional on T=t, the "correct" shape function is shifted which yields
  //the same stationary max-stable process; OK because there is a finite number
  //of values for t!

  for (i=0; i<lgtotalpoints; i++) {
    lgres[i] -= maxval;
  }

  return;
    
}


void range_BRmixed(cov_model *cov, range_type *range) {
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

int structBRuser(cov_model *cov, cov_model **newmodel) {
 
  location_type *loc = Loc(cov);
  cov_model *sub = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  int i, d, err, model_intern, 
      dim = sub->tsdim,
      newxlen;
  bool grid;
  double centreloc[MAXMPPDIM], minloc[MAXMPPDIM], maxloc[MAXMPPDIM],
    *newx= NULL, 
    **xgr = loc->xgr;

  ASSERT_NEWMODEL_NULL;
  if (cov->role != ROLE_BROWNRESNICK) BUG;
  
  assert(isBRuserProcess(cov));  
  
  if (loc->Time || (loc->grid && loc->caniso != NULL)) {
    TransformLoc(cov, false, GRIDEXPAND_AVOID, false);
    SetLoc2NewLoc(sub, PLoc(cov));
  }
  
  loc = Loc(cov);
  grid = loc->grid;
  
  model_intern = (cov->nr == BRORIGINAL_USER) ? BRORIGINAL_INTERN
               : (cov->nr == BRMIXED_USER) ? BRMIXED_INTERN
               : (cov->nr == BRSHIFTED_USER) ? BRSHIFTED_INTERN
               : BRORIGINAL_USER;
	       
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if (cov->Sgen == NULL) NEW_STORAGE(gen);
  
  GetDiameter(loc, minloc, maxloc, centreloc);
  newxlen = loc->lx;
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
    for (i=0; i<loc->lx; i++)
      for(d=0; d<dim; d++) newx[i*dim+d] = loc->x[i*dim+d] - centreloc[d];
  }

  if ((err = loc_set(newx, NULL, dim, dim, newxlen, false, loc->grid,
                     loc->distances,  cov)) > NOERROR)
    goto ErrorHandling;
  SetLoc2NewLoc(sub, PLoc(cov));

  if ((err=covCpy(&(cov->key), sub)) > NOERROR) goto ErrorHandling;
  
  if (cov->sub[MPP_TCF] != NULL) {
     if ((err = STRUCT(sub, &(cov->key))) > NOERROR) goto ErrorHandling;
     assert(cov->key->calling == cov);
  } 

  addModel(&(cov->key), model_intern);
  
  kdefault(cov->key, GEV_XI, P0(GEV_XI));
  kdefault(cov->key, GEV_MU, P0(GEV_MU));
  kdefault(cov->key, GEV_S, P0(GEV_S));

  if (cov->nr == BRMIXED_USER) {
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
  
  cov->key->calling = cov;
  if ((err = CHECK(cov->key, dim, dim, PointShapeType, cov->domown, 
                   cov->isoown, 1, ROLE_BROWNRESNICK)) == NOERROR) {
    if ((err = STRUCT(cov->key, NULL)) <= NOERROR) {
      err = CHECK(cov->key, dim, dim, PointShapeType, cov->domown, cov->isoown,
		  1, ROLE_BROWNRESNICK);
    }
  }

  assert(cov->key->calling == cov);
  
   
  ErrorHandling:
  FREE(newx);
  return err;
  
}


// new
int structBRintern(cov_model *cov, cov_model **newmodel) {
  cov_model *sub = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  location_type *loc = Loc(cov);
  int i, d, j, err,
    dim = sub->tsdim, 
    totaldim = loc->totalpoints*dim,
    zeropos = 0, // default, mostly overwritten
    newxlen = 3; // default (grid)
  bool grid = loc->grid;
  double  norm, step, mindist, dist,
    **xgr = loc->xgr;
  br_storage *sBR = NULL;

  ASSERT_NEWMODEL_NULL;
  if (cov->role != ROLE_BROWNRESNICK) BUG;

  assert(isPointShape(cov));
    
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if (cov->Sgen == NULL) NEW_STORAGE(gen);
  NEW_STORAGE(br);
  sBR = cov->Sbr;

  
  if (cov->sub[MPP_TCF] != NULL) {
    if ((err = STRUCT(sub, &(cov->key))) > NOERROR) goto ErrorHandling;
    cov->key->calling = cov;
  } else if ((err=covCpy(&(cov->key), sub)) > NOERROR) goto ErrorHandling;
 
  if ((err = CHECK(cov->key, dim, dim, VariogramType, cov->domown, SYMMETRIC, 
		   1, ROLE_COV)) != NOERROR) goto ErrorHandling;  

  if ((err = covCpy(&(sBR->submodel), cov->key)) != NOERROR) goto ErrorHandling;
  if ((err = CHECK(sBR->submodel, 1, 1, VariogramType, XONLY, 
		   ISOTROPIC, 1, ROLE_COV)) != NOERROR) goto ErrorHandling;
  
  if ((err = newmodel_covCpy(&(sBR->vario), VARIOGRAM_CALL, cov->key))!=NOERROR)
      goto ErrorHandling;
  if ((err = alloc_cov(sBR->vario, dim, 1, 1)) != NOERROR) goto ErrorHandling;
	
  addModel(&(cov->key), GAUSSPROC, cov);
  assert(cov->key->ownloc == NULL);
  
  if (cov->nr == BRORIGINAL_INTERN) {
    if (!grid) {
      zeropos = loc->lx;
      for (i=0; i<loc->lx; i++) {
	norm = 0.0;
	for (d=0; d<dim; d++) {
	  norm += loc->x[i*dim+d]*loc->x[i*dim+d];
	}
	if (norm<1e-8) zeropos = i;
      }
      newxlen = loc->lx + (zeropos == loc->lx);
      if ((sBR->newx = (double*) MALLOC(dim*newxlen*sizeof(double))) == NULL) {
        err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      for(d=0; d<dim; d++) {
        for (i=0; i<loc->lx; i++) sBR->newx[i*dim+d] = loc->x[i*dim+d];
	if (zeropos == loc->lx) sBR->newx[loc->lx*dim+d] = 0.0;
      }  
    } 
  } else if (cov->nr == BRMIXED_INTERN) {

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
    }
    if (mindist < step) {     
      PRINTF("Warning! meshsize is larger than resolution of simulation! meshsize is automatically decreased to %f.\n", mindist);
      P(BR_MESHSIZE)[0] = step = mindist;
    }
    
    if (!PisNULL(BR_OPTIMAREA)) {
      sBR->minradius =  cov->ncol[BR_OPTIMAREA] * step;
    } else {
      sBR->minradius = 0;
    }    
    
    double alpha = -1,  yy, C0, gammamin, xx=step * 1e-6;
    COV(ZERO, sBR->submodel, &C0);
    COV(&xx, sBR->submodel, &yy);
    alpha = LOG(C0 - yy) / LOG(xx);
    if (alpha > 2.0) alpha = 2.0;
    gammamin = 4.0 - 1.5 * alpha;
    
    INVERSE(&gammamin, sBR->submodel, &xx);
    xx = CEIL(xx/step) * step;
    sBR->minradius = FMAX(sBR->minradius, xx); 

    yy = P0(BR_VARIOBOUND);
    INVERSE(&yy, sBR->submodel, &(sBR->radius));
    
    sBR->radius = CEIL(sBR->radius / step) * step;

    newxlen = 3;
    grid = true;

    if ((sBR->newx = (double*) MALLOC(newxlen*dim*sizeof(double))) == NULL) {
      err = ERRORMEMORYALLOCATION; goto ErrorHandling;
    }

    if ((err = covCpy(sBR->sub, cov->key)) != NOERROR) goto ErrorHandling;
    for (d=0; d<dim; d++) {
      sBR->newx[3*d+XSTART] = -sBR->radius - sBR->minradius;
      sBR->newx[3*d+XLENGTH] =
	2 * ((int) ROUND((sBR->radius + sBR->minradius) / step)) + 1;
      sBR->newx[3*d+XSTEP] = step;
    }

    err = loc_set(sBR->newx, NULL, dim, dim, newxlen, false, grid, 
                  false, sBR->sub[0]);
    double **subxgr = Loc(sBR->sub[0])->xgr;
    zeropos = 0;
    for (d = dim; d > 0; d--) {
      double len =  subxgr[d-1][XLENGTH];
      zeropos = zeropos * (int) len + (int) CEIL(len * 0.5) - 1;
    }
    sBR->zeropos = zeropos;

    if ((err = CHECK(sBR->sub[0], dim, dim, ProcessType, cov->domown, 
 		       cov->isoown, 1, ROLE_GAUSS)) == NOERROR) {
       if ((err = STRUCT(sBR->sub[0], NULL)) <= NOERROR) {
         err = CHECK(sBR->sub[0], dim, dim, ProcessType, cov->domown,
                     cov->isoown, 1, ROLE_GAUSS); 
       }
     }
     if (err > NOERROR) goto ErrorHandling;
     assert(sBR->sub[0]->calling == cov);

  } else { //  END BRMIXED;  START SHIFTED    
    assert(cov->nr == BRSHIFTED_INTERN);
  }
 

  if (sBR->newx == NULL) {
    if ((err = loc_set(grid ? * loc->xgr : loc->x, NULL, dim, dim,
		       grid ? 3 : loc->totalpoints, false, grid,
		       loc->distances, cov->key)) != NOERROR)
      goto ErrorHandling;
  } else {
    if ((err = loc_set(sBR->newx, NULL, dim, dim, newxlen, false, grid, false,
		       cov->key)) != NOERROR) goto ErrorHandling;
  }
  xgr = loc->xgr;

 
  if (cov->nr != BRMIXED_INTERN && grid) {
    double **subxgr = Loc(cov->key)->xgr;
    for (d=dim; d>0; d--) {
      double len =  subxgr[d-1][XLENGTH];
      zeropos = zeropos * len + (int) CEIL(len * 0.5) - 1;
    }
    sBR->zeropos = zeropos;
  }
  
  if ((err = CHECK(cov->key, dim, dim, ProcessType, cov->domown, 
                   cov->isoown, 1, ROLE_GAUSS)) == NOERROR) {
    if ((err = STRUCT(cov->key, NULL)) <= NOERROR) {
      err = CHECK(cov->key, dim, dim, ProcessType, cov->domown,
		  cov->isoown, 1, ROLE_GAUSS); 
    }
  }
  if (err > NOERROR) goto ErrorHandling;

  assert(cov->key->calling == cov);
  
  ErrorHandling:
  if (err != NOERROR) br_DELETE(&(cov->Sbr));    
  return err;
}

int structBrownResnick(cov_model *cov, cov_model **newmodel) {
  
  int d, err, meth, role,
      dim = cov->tsdim;
  double  maxcov,
      minloc[MAXMPPDIM], maxloc[MAXMPPDIM],
      centreloc[MAXMPPDIM], maxdist[MAXMPPDIM];      
  cov_model *next = cov->sub[0];
  location_type *loc = Loc(cov); 
  
  if (cov->role != ROLE_BROWNRESNICK) BUG;
  if (loc->Time || (loc->grid && loc->caniso != NULL)) {
    TransformLoc(cov, false, GRIDEXPAND_AVOID, false);
    SetLoc2NewLoc(next, PLoc(cov));
  }
  loc = Loc(cov);
  
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
       return ERRORDIM;
  if (cov->key != NULL)
    COV_DELETE(&(cov->key));
    
  if (cov->role == ROLE_SMITH) {
    if (!cov->logspeed) 
      SERR2("'%s' requires a variogram model as submodel that tends to infinity with rate of at least 4log(h) for being compatible with '%s'", NICK(cov), CovList[SMITHPROC].nick);
    cov_model *calling = cov->calling;
    double newscale, 
      factor =  INVSQRTTWO;
      
    ASSERT_NEWMODEL_NULL;
    if (next->full_derivs < 0) SERR("given submodel does not make sense");
 
    while (isDollar(next)) {
      addModel(&(cov->key), DOLLAR);
      newscale = 1.0;
      if (!PARAMisNULL(next, DSCALE) ) newscale *= PARAM0(next, DSCALE);
      if (!PARAMisNULL(next, DVAR) )  newscale /= SQRT(PARAM0(next, DVAR));
      if (factor != 1.0) {
	newscale *= factor;
	factor = 1.0;
      }
      return ERRORNOTPROGRAMMEDYET;
      

      kdefault(calling, DSCALE, newscale);
      next = next->sub[0]; // ok
    }

    if (cov->sub[MPP_TCF] != NULL) {
      return ERRORNOTPROGRAMMEDYET;
    } 

    if (next->nr == BROWNIAN && PARAM0(next, BROWN_ALPHA) == 2.0) {
      addModel(&(cov->key), GAUSS);   // ?? 
      if (factor != 1.0) {
	addModel(&(cov->key), DOLLAR);
	kdefault(cov->key, DSCALE, factor);      
      }
    } else {
      SERR("Smith process with BrownResnick tcf only possible for fractal Brownian motion with alpha=2");
    }
  } else if (cov->role == ROLE_BROWNRESNICK) {
    if (next->role == ROLE_BROWNRESNICK) {
      SERR1("submodel of '%s' must be a covariance model or tcf", 
	    NICK(cov));
    } else {
      role = isVariogram(next) ? ROLE_COV : ROLE_UNDEFINED;
	
      ASSERT_ROLE_DEFINED(next);  
      if (((err = covCpy(&(cov->key), next)) != NOERROR)
         || ((err = CHECK(cov->key, dim, dim, VariogramType, XONLY,
			  SYMMETRIC, 1, role)) != NOERROR)) {
	return err;
      }

      GetDiameter(loc, minloc, maxloc, centreloc);
      for (d=0; d<MAXMPPDIM; d++) maxdist[d] = 0.5*(maxloc[d] - minloc[d]);
   
      cov_model *K = NULL;
      if ((err = newmodel_covCpy(&K, VARIOGRAM_CALL, cov->key, maxdist, NULL,
				 NULL, dim, dim, 1, 0, false, false, false))
	  != NOERROR) return err;
      if ((err = alloc_cov(K, dim, 1, 1)) != NOERROR) return err;
      if (K->sub[0] != NULL) SetLoc2NewLoc(K->sub[0], PLoc(K));
      Variogram(NULL, K, &maxcov);
      COV_DELETE(&K);
      if (isPosDef(next) || maxcov <= 4.0) {
	meth = BRORIGINAL_USER;  
      } else if (!next->logspeed || next->logspeed <= 4.0 || maxcov <= 10.0) {
	meth = BRSHIFTED_USER;
      } else {
	meth = BRMIXED_USER;
      }

      addModel(&(cov->key), meth, cov);
      cov_model *key = cov->key;
      key->prevloc = PLoc(cov);
      
      kdefault(key, GEV_XI, P0(GEV_XI));
      kdefault(key, GEV_MU, P0(GEV_MU));
      kdefault(key, GEV_S, P0(GEV_S));
       
      if ((err =  CHECK(key, dim, dim, BrMethodType, cov->domown, cov->isoown,
	                1, ROLE_BROWNRESNICK)) == NOERROR) {
      if ((err = STRUCT(key, NULL)) <= NOERROR) {
        err = CHECK(key, dim, dim, BrMethodType, cov->domown, cov->isoown,
		    1, ROLE_BROWNRESNICK);
	}
      }
      if (err > NOERROR) return(err);
    }
  } else {
    ILLEGAL_ROLE;
  }

  // need check to be called?

  return NOERROR;

  
}

int initBrownResnick (cov_model *cov, gen_storage *S) {

  cov_model *sub = cov->key != NULL ? cov->key :
                   cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF]; 
  int err;

  assert(cov->nr == BROWNRESNICKPROC);

  if (cov->role == ROLE_BROWNRESNICK) {
    if (cov->key != NULL) {
      sub->simu.active = true;
      sub->simu.expected_number_simu = cov->simu.expected_number_simu;
      if ((err = INIT(sub, 0, S)) != NOERROR) return err;
      cov->fieldreturn = true;
      cov->origrf = false;
      cov->rf = sub->rf;
    }
  } else ILLEGAL_ROLE;
   
  return NOERROR;

}

void doBrownResnick(cov_model *cov, gen_storage *s) {

  assert(!cov->origrf);
  assert(cov->key != NULL);
  cov_model *key = cov->key;
 
  PL++;
  DO(key, s); // nicht gatternr
  PL--;

}


int initBRuser (cov_model *cov, gen_storage *S) {
  location_type *loc = Loc(cov);
  cov_model *sub = cov->key != NULL ? cov->key :
                  cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];   
  int  err, maxpoints = GLOBAL.extreme.maxpoints;
  
  assert(isBRuserProcess(cov));

  if (cov->role == ROLE_BROWNRESNICK) {
    if (loc->distances) return ERRORFAILED;
    
    if (cov->key != NULL) {
      sub->simu.active = true;
      double ens = ((double) cov->simu.expected_number_simu) * maxpoints;
      sub->simu.expected_number_simu = (int) MIN(ens, (double) MAXINT);

      if ((err = INIT(sub, 1, S)) != NOERROR) return err;
      FieldReturn(cov); 
    }
    
    return NOERROR;
  }

  else ILLEGAL_ROLE;
   
}
