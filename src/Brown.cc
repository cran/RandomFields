/* 
 Authors
 Marco Oesting,
 Martin Schlather, schlather@math.uni-mannheim.de

 simulation of Brown-Resnick processes

 Copyright (C) 2009 -- 2010 Martin Schlather 
 Copyright (C) 2011 -- 2014 Marco Oesting & Martin Schlather 

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

#include <math.h>
#include <stdio.h>
#include "RF.h"
#include "Covariance.h"

#define BR_MESHSIZE LAST_MAXSTABLE + 1
#define BR_VERTNUMBER LAST_MAXSTABLE + 2
#define BR_OPTIM LAST_MAXSTABLE + 3
#define BR_OPTIMTOL LAST_MAXSTABLE + 4
#define BR_OPTIMMAX LAST_MAXSTABLE + 5
#define BR_LAMBDA LAST_MAXSTABLE + 6
#define BR_OPTIMAREA LAST_MAXSTABLE + 7
#define BR_VARIOBOUND LAST_MAXSTABLE + 8


//////////////////////////////////////////////////////////////////////
// Brown Resnick

int checkBrownResnickProc(cov_model *cov) {
  cov_model  
    *key = cov->key,
    *sub = key != NULL ? key :
            cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  int err, role, 
      dim = cov->tsdim;
  isotropy_type isoprev;
  Types type;

  
  ASSERT_ONE_SUBMODEL(cov);
  if ((err = SetGEVetc(cov, ROLE_BROWNRESNICK)) != NOERROR) return err;
  
  role = isNegDef(sub) ? ROLE_COV
    : (isGaussProcess(sub) &&  isPointShape(cov)) ? ROLE_GAUSS
    : isBrownResnickProcess(sub) ? ROLE_BROWNRESNICK
    : isPointShape(sub) ? ROLE_BROWNRESNICK
    : ROLE_UNDEFINED;    

  type = isProcess(sub) || isPointShape(sub) 
    ? CovList[sub->nr].Type : NegDefType;
    
  ASSERT_ROLE_DEFINED(sub);  
  
  isoprev = role == ROLE_COV ? SYMMETRIC : CARTESIAN_COORD;
  
  if ((err = CHECK(sub, dim, dim, 
		   type,//Martin: changed 27.11.13 CovList[cov->nr].Type,
		   XONLY, isoprev, 1, role)) != NOERROR) {
      return err;
  }
  
  setbackward(cov, sub);
  if (cov->vdim2[0] != 1) SERR("BR only works in the univariate case");

  return NOERROR;  
}

#define ORIG_IDX 0
int init_BRorig(cov_model *cov, gen_storage *s){
  cov_model *key = cov->key;
  location_type *keyloc = NULL;
  int err, d, 
    dim = cov->tsdim;
  bool keygrid;
  BR_storage *sBR = NULL;

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
  
  GLOBAL.gauss.loggauss = false; 
  key->simu.active = true;
  key->simu.expected_number_simu = cov->simu.expected_number_simu;
  
  assert(cov->mpp.moments >= 1);
  if ((err = INIT(key, 1, s)) != NOERROR) goto ErrorHandling;  
  
  cov->loggiven = true;
  
  
  assert(key->nr == GAUSSPROC);
  assert(key->mpp.moments >= 1);
  
  cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
  cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1.0;
  cov->mpp.maxheights[0] = exp(GLOBAL.extreme.standardmax);

  pgs->zhou_c = 1.0;
  
  sBR = cov->SBR;
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
		     keyloc->distances, &Loc(cov->SBR->vario))) > NOERROR)
    goto ErrorHandling;
  if (sBR->vario->sub[0] != NULL) 
    SetLoc2NewLoc(sBR->vario->sub[0], Loc(sBR->vario));
  // assert(false);
  Variogram(NULL, sBR->vario, sBR->trend[ORIG_IDX]);
  
  if ((err = FieldReturn(cov)) != NOERROR)  // must be later than INIT !
   goto ErrorHandling;  

 ErrorHandling:
  if (err != NOERROR) BR_DELETE(&(cov->SBR));
  return err;
}

void do_BRorig(cov_model *cov, gen_storage *s) {
  cov_model *key = cov->key;
  assert(key != NULL);
  BR_storage *sBR = cov->SBR;
  assert(sBR != NULL && sBR->trend != NULL);
  res_type *res = cov->rf;
#define ORIG_IDX 0
  int i,
    zeropos = sBR->zeropos[ORIG_IDX];
  long
    totalpoints = Loc(cov)->totalpoints;
  assert(totalpoints > 0);
  assert(zeropos >= 0 && zeropos <totalpoints);
  double *trend = sBR->trend[ORIG_IDX];
  assert(cov->origrf);

  
  DO(key, s); // nicht gatternr,
  res_type *lgres = key->rf;
  double lgreszeropos = (double) lgres[zeropos]; // wird durch DO veraendert!
  for (i=0; i<totalpoints; i++) {
    // printf("zeropos=%d i=%d %f %f %f\n", zeropos, i, trend[i], lgres[i], lgres[zeropos]);
    res[i] = lgres[i] - lgreszeropos - trend[i];
    //    printf("%d %e %e %e %f lgres=%lu %d %f %f\n", i, (double) res[i], (double) lgres[i], lgreszeropos , trend[i],   lgres, zeropos, (double) lgres[zeropos], (double) lgres[0]);
    //    assert(i < 2);
  }
}

int init_BRshifted(cov_model *cov, gen_storage *s) {
  cov_model *key=cov->key;
  location_type *keyloc = NULL;
  int d, dim, err;
  long j, shiftedloclen, keytotal, trendlenmax, trendlenneeded;
  bool keygrid;
  BR_storage *sBR = NULL;
  
  
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
      
      GLOBAL.gauss.loggauss = false;
      key->simu.active = true;
      key->simu.expected_number_simu = cov->simu.expected_number_simu; 
      assert(cov->mpp.moments >= 1);
      if ((err = INIT(key, 1, s)) != NOERROR) return err;  
      
      cov->loggiven = true;
      
      assert(key->nr == GAUSSPROC);
      assert(key->mpp.moments >= 1);
      cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
      cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1.0;
      cov->mpp.maxheights[0] = exp(GLOBAL.extreme.standardmax);
      pgs->zhou_c = 1.0;
      
      sBR = cov->SBR;
       
      shiftedloclen = keygrid ? 3 : keytotal;
      if ((sBR->shiftedloc = (double*)
	   MALLOC(dim*shiftedloclen*sizeof(double))) == NULL ||
	  (sBR->locindex = (int*) MALLOC(sizeof(int) * dim))==NULL
	  ) {
        err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      
      trendlenmax = (int) ceil((double) GLOBAL.br.BRmaxmem / keytotal);
      trendlenneeded = MIN(keytotal, cov->simu.expected_number_simu);
      sBR->trendlen = MIN(trendlenmax, trendlenneeded);
   
      //PMI(cov->calling);
      //      printf("len=%d %d %d mem=%d %d exp=%d %d\n", sBR->trendlen, trendlenmax, trendlenneeded, GLOBAL.br.BRmaxmem, keytotal, cov->simu.expected_number_simu, GLOBAL.general.expected_number_simu);



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
	                 keyloc->distances, &Loc(sBR->vario))) > NOERROR)
        return err;
      if (sBR->vario->sub[0] != NULL) 
        SetLoc2NewLoc(sBR->vario->sub[0], Loc(sBR->vario));   
      
      if ((err = FieldReturn(cov)) != NOERROR)  // must be later than INIT !
        return err;
    }
        
    return NOERROR;
  }
  
  else ILLEGAL_ROLE;
  
  ErrorHandling:
    BR_DELETE(&(cov->SBR));
    return err;

}

void do_BRshifted(cov_model *cov, gen_storage *s) {
  BR_storage *sBR = cov->SBR;
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
         **trend = sBR->trend;
  assert(cov->origrf);
  res_type *res = cov->rf;
  res_type *lgres = cov->key->rf;

  DO(key, s);
  zeropos = (long) floor(UNIFORM_RANDOM * keytotal);
  
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
      indextrafo(zeropos, keyloc->length, dim, locindex); // to do: ersetzen
      for (d=0; d<dim; d++) {
	shiftedloc[3*d+XSTART]  = -locindex[d]*keyloc->xgr[d][XSTEP];
	shiftedloc[3*d+XLENGTH] = keyloc->xgr[d][XLENGTH];
	shiftedloc[3*d+XSTEP]   = keyloc->xgr[d][XSTEP];
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
        SetLoc2NewLoc(sBR->vario->sub[0], Loc(sBR->vario));
    //assert(false); 
    Variogram(NULL, sBR->vario, sBR->trend[trendindex]);
    
    mem2loc[trendindex] = zeropos; // todo ? long instead of int
    loc2mem[zeropos] = trendindex;      
  }
  
  for (i=0; i<keytotal; i++) {
    res[i] = lgres[i] - lgres[zeropos] - trend[trendindex][i];  
  }
  
}

int check_BRmixed(cov_model *cov) {
  
  int err;
  br_param *bp = &(GLOBAL.br);
   
  ASSERT_ONE_SUBMODEL(cov);
  
  if (!cov->logspeed) SERR("BrownResnick requires a variogram model as submodel that tends to infinity [t rate of at least 4log(h) for being compatible with BRmixed");
  
  kdefault(cov, BR_MESHSIZE, bp->BRmeshsize);
  kdefault(cov, BR_VERTNUMBER, bp->BRvertnumber);
  kdefault(cov, BR_OPTIM, bp->BRoptim);
  kdefault(cov, BR_OPTIMTOL, bp->BRoptimtol);
  kdefault(cov, BR_OPTIMMAX, bp->BRoptimmaxpoints);
  kdefault(cov, BR_VARIOBOUND, bp->variobound);
  if (cov->nr == BRMIXED_USER && cov->key == NULL && P0INT(BR_OPTIM) > 0) {
    if (!PisNULL(BR_LAMBDA)) {
      if (PisNULL(BR_OPTIMAREA)) SERR1("'%s' not given", KNAME(BR_OPTIMAREA));
      if (PL > 0) PRINTF("'%s' set to '0'", KNAME(BR_OPTIM));
      PINT(BR_OPTIM)[0] = 0;
    } else if (P0INT(BR_OPTIM) == 2 && !PisNULL(BR_OPTIMAREA))
      if (PL > 0) PRINTF("'%s' set to '1'", KNAME(BR_OPTIM));
  }

  kdefault(cov, BR_LAMBDA, RF_NA);
  if (PisNULL(BR_OPTIMAREA)) //necessary as areamat might not be scalar
    kdefault(cov, BR_OPTIMAREA, 1.0);
  
  if ((err = checkBrownResnickProc(cov)) != NOERROR) return err;
  
  if ((err = checkkappas(cov, true)) != NOERROR) return err;
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
    return ERRORDIM;
 if (cov->vdim2[0] != 1) SERR("BR only works in the univariate case");

  return NOERROR; 
}

void kappaBRmixed(int i, cov_model *cov, int *nr, int *nc){
  // i nummer des parameters

  int dim = cov->tsdim;
  
  switch(i) {
    
    case GEV_XI: case GEV_MU: case GEV_S: case BR_MESHSIZE:
    case BR_VERTNUMBER: case BR_OPTIM: case BR_OPTIMTOL: 
    case BR_OPTIMMAX: case BR_LAMBDA: case BR_VARIOBOUND:
      *nr = 1;
      *nc = 1;
    break;
    
    case BR_OPTIMAREA:
      switch(dim) {
	case 1:
	  *nr = 1;
          *nc = SIZE_NOT_DETERMINED;
        break;
	case 2:
	  *nr = SIZE_NOT_DETERMINED;
	  *nc = SIZE_NOT_DETERMINED;
	break;
	default:
	  *nr = 1;
	  *nc = 1;
      }
    break;

    default:
      *nr = -1; 
      *nc = -1;
  }
}



void OptimArea(cov_model *cov, int idx) {
  // side effect auf P(BR_OPTIMAREA) !!
  BR_storage *sBR = cov->SBR;
  cov_model *key = sBR->sub[idx];
  location_type *keyloc = Loc(key);
  int
    d, i, j, k, l, //maxind, 
    // totalcount = 0,
    //    limit, 
    zeroxpos, zeroypos, cellcounter,
    *len = keyloc->length,
    radius = (int) floor(keyloc->length[0] / 2.0),
    vertnumber = P0INT(BR_VERTNUMBER),
    radiusP1 = radius+1,
    //err = NOERROR,
    dim = cov->tsdim,
    //totalpoints = keyloc->totalpoints,
    optimmax = P0INT(BR_OPTIMMAX)
    //zeropos = sBR->zeropos[idx]
    ;
  double 
    dummy, //maxval,
    Error, 
    maxErrorbound, Errorbound, Errorboundtmp, 
    Errortol = P0(BR_OPTIMTOL),
    step = P0(BR_MESHSIZE),
    stepdim= intpow(step, dim),
    **am = sBR->areamatrix
    //*trend = sBR->trend[idx]
    ;

  // MARCO: ok??
  // zwingend eine Kopie erstellen, da Vertauschen durchgefuehrt werden
  for (j=0; j<vertnumber; j++) { 
    //  totalcount = 0;
    double *areamatrix = am[j];
    int *countvector = sBR->countvector[j];
    for (d=0; d<radiusP1; d++) {
      areamatrix[d] = (double) countvector[d];
      //printf("d=%d %d\n", d, countvector[d]);
    }
  }
    
  double factor = optimmax * stepdim;
  for (d=0; d<radiusP1; d++) {
    double addnumber = (d==0) ? 1 : ((dim==1) ? 2 : 4*d);
    double invfactor = 1.0 / (factor * addnumber);
    for (j=0; j<vertnumber; j++) am[j][d] *= invfactor;
  }
  
  double lambda=0.0;
  for (j=0; j<vertnumber; lambda += am[j++][0]);
  for (j=0; j<vertnumber; j++) am[j][0] = lambda / vertnumber;
  
  for (d=0; d<radius; d++) { //monotonic in each column
    for (j=0; j<vertnumber; j++) {
      for (k=j+1; k<vertnumber; k++) {
	if (am[j][d] < am[k][d]) {
	  dummy = am[j][d];
	  am[j][d] = am[k][d];
	  am[k][d] = dummy;
	}
      }
      am[j][d] = fmin(am[j][d], lambda/vertnumber); //cutoff
    }
  }
  
  for (j=0; j<vertnumber; j++) { //monotonic in each row
    for (d=0; d<radiusP1; d++) {
      for (k=d+1; k<radiusP1; k++) {
	if (am[j][d] < am[j][k]) {
	  dummy = am[j][d];
	  am[j][d] = am[j][k];
	  am[j][k] = dummy;
	}
      }
    }
  }
  
  double newlambda = 0.0; //correction for bias because of cutoff
  for (d=0; d<radiusP1; d++) {
    double addnumber = (d==0) ? 1 : ((dim==1) ? 2 : 4*d);  
      for (j=0; j<vertnumber; j++) {
	newlambda += am[j][d] * addnumber;      
      }
  }
  double invstepdimlambda = 1.0 / (stepdim*newlambda);
  for (d=0; d<radiusP1; d++) {
    for (j=0; j<vertnumber; j++) {
      am[j][d] *= invstepdimlambda;
    }
  }
  lambda *= invstepdimlambda;
  
  maxErrorbound = 0.0;
  for (d=0; d<radiusP1; d++) { //areamatrix => Error matrix
    for (j=0; j<vertnumber; j++) {
      am[j][d] = lambda/vertnumber - am[j][d];
      maxErrorbound = fmax(maxErrorbound, am[j][d]);
    }
  }
  
  Error = 0.0;
  Errorboundtmp = Errorbound = 0.0;  
  while (Errorboundtmp < maxErrorbound && Error < Errortol) {
    Errorbound = Errorboundtmp; 
    Errorboundtmp = maxErrorbound;
    for (d=0; d<radiusP1; d++) {
      for (j=0; j<vertnumber; j++) {
	if (am[j][d] > Errorbound) {
	  Errorboundtmp = fmin(Errorboundtmp, am[j][d]); 
	}
      }      
    }
    Error = 0.0;
    cellcounter = 0;
    for (d=0; d<radiusP1; d++) {
      long addnumber = (d==0) ? 1 : ((dim==1) ? 2 : 4*d);
      for (j=0; j<vertnumber; j++) {
	if (am[j][d] <= Errorboundtmp + 1e-6) {
	  Error += (double) addnumber * am[j][d]; // geaendert
	  cellcounter += addnumber;
	}
      }
    }
    Error = Error / (cellcounter*lambda/vertnumber);
  }
  //Error = PP(U_i < x_0 | (U_i + M_i, T_i) \in x_0 + E)
  // = EE(#{i: U_i<0, (U_i + M_i, T_i) \in E})/EE(#{i: (U_i + M_i, T_i) \in E})
  
  zeroxpos = (int) floor(len[0] / 2.0);
  zeroypos = (dim==1) ? 0 : (int) floor(len[1] / 2.0);
  newlambda = 0.0;
  cellcounter = 0;
  //printf("optimarea: (%d, %d)\n", cov->nrow[BR_OPTIMAREA], cov->ncol[BR_OPTIMAREA]);
  
  double invvertnumber = 1.0 / vertnumber;
  PFREE(BR_OPTIMAREA);
  PALLOC(BR_OPTIMAREA, (int) (dim==1) ? 1 : len[1], len[0]);
  double 
    *optimarea = P(BR_OPTIMAREA);
  for (l=i=0; i<cov->nrow[BR_OPTIMAREA]; i++) {
    for (k=0; k<cov->ncol[BR_OPTIMAREA]; k++, l++) {
      d = abs(k-zeroxpos) + abs(i-zeroypos); // abs OK
      optimarea[l] = 0.0;
      for (j=0; j<vertnumber; j++) {
	if (d<radiusP1 && (am[j][d] <= Errorbound + 1e-6)) {
	  optimarea[l] += invvertnumber;
	  newlambda += lambda - am[j][d]*vertnumber;
	  cellcounter++;
	}
      }
      // if (P(BR_OPTIMAREA][l] > 0) printf("%4.2f ", P(BR_OPTIMAREA][l]);
    }
    //printf("\n");
  }


  
  // return newlambda / cellcounter; // not used
}



int prepareBRoptim(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  BR_storage *sBR = cov->SBR;
  int idx = 0; // first is the largest
  cov_model *key = sBR->sub[idx];
  location_type *keyloc = Loc(key);
  int i, j, d, sum_directions,
    radius = (int) floor(keyloc->length[0] / 2.0),
    radiusP1 = radius+1,
    vertnumber = P0INT(BR_VERTNUMBER),
    //    totalpoints = keyloc->totalpoints,
    *len = keyloc->length,
    dim = cov->tsdim;
 
 
  // if (PL >= PL_STRUCTURE) PRINTF("starting BR optimiation with %d simulations [optim=%d]...\n", optimmax, optim); APMI(cov); assert(false);

  switch(P0INT(BR_OPTIM)) {
  case 0:   
    if (ISNAN(P0(BR_LAMBDA))) P(BR_LAMBDA)[0] = 1.0;
    break;
  case 1:
    // nothing to do
    break;
  case 2:
    if (dim > 2) BUG;
    sBR->vertnumber = vertnumber; 
    sum_directions  = 0;
    for (d=0; d<dim; sum_directions += len[d++]);
    //    printf("dim=%d %d\n", dim, sum_directions); PMI(key);
   

    if (sBR->countvector != NULL || sBR->areamatrix != NULL) BUG;
    if ((sBR->countvector = (int**) CALLOC(vertnumber, sizeof(int*))) == NULL || 
	(sBR->areamatrix = (double**) CALLOC(vertnumber, sizeof(double*))) == NULL||
	(sBR->logvertnumber = (double *) MALLOC(vertnumber * sizeof(double))) ==NULL
	) return ERRORMEMORYALLOCATION;
    for (i=0; i< vertnumber; i++) {
      if ((sBR->countvector[i] = (int*) CALLOC(sum_directions, sizeof(int))) 
	  == NULL ||
	  (sBR->areamatrix[i] = (double*) CALLOC(radiusP1, sizeof(double)))
	  == NULL)
	return ERRORMEMORYALLOCATION;
    }
    for (j=0; j<vertnumber; j++)
      sBR->logvertnumber[j] = - log((double) (j+1)/vertnumber);    
    break;

  default:
      SERR("optimization might not be used here\n");
  }

  if (PL >= PL_STRUCTURE) PRINTF("BR optimisation finished...\n");

 return NOERROR;	
}

void indextrafo(long onedimindex, int *length, int dim, int *multidimindex) {
  int d;
  for (d=0; d<dim; d++) {
    multidimindex[d] = onedimindex % length[d];
    onedimindex = onedimindex / length[d];    
  }
}


void set_lowerbounds(cov_model *cov) {    
  BR_storage *sBR = cov->SBR;
  assert(sBR != NULL);
  assert(sBR->sub[0] != NULL);
  double 
    *optimarea = P(BR_OPTIMAREA);
  //  assert(optimarea[0] == 1.0); Marco: wann ist dies erfuellt ???
  int 
    i, l, ii,
    xbound = (int) floor(cov->ncol[BR_OPTIMAREA] * 0.5),
    ybound = (int) floor(cov->nrow[BR_OPTIMAREA] * 0.5);

  //printf("lowerbounds: xy= %d %d\n", xbound, ybound);

  for (i=0; i<=sBR->maxidx; i++) {
    cov_model *key = sBR->sub[i];
    location_type *keyloc = Loc(key);
    long j, k, keytotal = keyloc->totalpoints,  
      len = keyloc->length[0];
 
    for (j=0; j<keytotal; j++) sBR->lowerbounds[i][j] = RF_INF;    
    for (l=0, ii=-ybound; ii<=ybound; ii++) {
      k = sBR->zeropos[i] + ii*len - xbound;
      for (j=-xbound; j<=xbound; j++, l++, k++) {
	if (optimarea[l]>1e-5) {
	  sBR->lowerbounds[i][k] = -log(optimarea[l]); 
	  //	  printf("-> idx=%d ii=%d j=%d lower=%f\n", i, ii,j,  
	  //	 sBR->lowerbounds[i][k] );
	}
	//	printf("k=%d i=%d %d %f l=%d %f %f\n", k, i, sBR->zeropos[i] , optimarea[l], l, sBR->lowerbounds[i][k], P0(BR_OPTIMAREA));
      }
    }
  }
}
 


int init_BRmixed(cov_model *cov, gen_storage *s) {

  //printf("entering init_brmixed\n");

  location_type 
    *loc = Loc(cov);
  //*keyloc = Loc(key)
  int  i, d, err, 
    //keytotal = keyloc->totalpoints,
    //*lglength = keyloc->length,
     dim = cov->tsdim,
    // optim = P0INT(BR_OPTIM),
    bytes = sizeof(double) * dim;
  BR_storage *sBR = cov->SBR;
  assert(sBR != NULL);
  assert(sBR->sub[0] != NULL);
  // step = P0(BR_MESHSIZE);
  pgs_storage *pgs = NULL;

  
  assert(isPointShape(cov));
  if (cov->role != ROLE_BROWNRESNICK) ILLEGAL_ROLE;
  assert(dim > 0);
   
  sBR->idx = -1;
  sBR->trendlen = sBR->maxidx + 1;

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
    pgs->supportmin[d] = RF_NEGINF;
    pgs->supportmax[d] = RF_INF;
     
   //printf("support %f %f radius=%f %f \n", sBR->locmax[d],   sBR->locmin[d], sBR->radius, step ); //assert(false); // woer kommt sBR->radius
    // support funktioniert nicht

    
  }
  for (d=0; d<=cov->mpp.moments; d++) {
    cov->mpp.mM[d] = cov->mpp.mMplus[d] = 1.0;
  }
  


 if ((err = prepareBRoptim(cov, s)) != NOERROR) {
    //printf("err %d\n", err);
     goto ErrorHandling;
  }

  for (i=0; i<=sBR->maxidx; i++) {
    cov_model *key = sBR->sub[i];
    location_type *keyloc = Loc(key);
    int keytotal = keyloc->totalpoints;
    //      *lglength = keyloc->length;
    //bool keygrid = keyloc->grid;
    assert(keyloc->grid);
 
    assert(key->nr == GAUSSPROC);
    GLOBAL.gauss.loggauss = false;
    key->simu.expected_number_simu = cov->simu.expected_number_simu;     
    assert(cov->mpp.moments >= 1);
    if ((err = INIT(key, 1, s)) != NOERROR) goto ErrorHandling;
    key->simu.active = true;
    cov->loggiven = true;
    
    assert(key->mpp.moments >= 1);
    cov->mpp.mM[0] = cov->mpp.mMplus[0] = 1.0;
    cov->mpp.mM[1] = cov->mpp.mMplus[1] = 1.0;
    
    
    if (sBR->trend[i] == NULL && 
	(sBR->trend[i] =(double*) MALLOC(keytotal*sizeof(double)))==NULL){
      err = ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    
    //printf("i=%d\n", i);

    if ((err = loc_set(keyloc->xgr[0], NULL, NULL, dim, dim, 3, 0,
		       false, true, keyloc->distances, 
		       &Loc(sBR->vario))) > NOERROR) goto ErrorHandling;

    if (sBR->vario->sub[0] != NULL)
      SetLoc2NewLoc(sBR->vario->sub[0], Loc(sBR->vario));
    Variogram(NULL, sBR->vario, sBR->trend[i]);

   
    if ((sBR->lowerbounds[i] = (double*) MALLOC(keytotal*sizeof(double))) 
	== NULL) { err=ERRORMEMORYALLOCATION; goto ErrorHandling; }
    assert(keyloc != NULL);
    key->simu.active = true;
  } // for i<=maxid2

   set_lowerbounds(cov);
   
    // assert(false);

  //printf("keytotal=%d\n", 45);
  cov->rf = sBR->sub[0]->rf;
  cov->origrf = false;
  cov->fieldreturn = HUETCHEN_OWN_GRIDSIZE;
  pgs->estimated_zhou_c = ISNAN(P0(BR_LAMBDA));
  pgs->zhou_c =  pgs->estimated_zhou_c ? 0.0 : P0(BR_LAMBDA);
  pgs->logmean = true;
  pgs->sq_zhou_c = pgs->sum_zhou_c = 0.0;
  pgs->n_zhou_c = 0;
  sBR->next_am_check = GLOBAL.br.deltaAM;

  return NOERROR;

 ErrorHandling:
  BR_DELETE(&(cov->SBR));

  return err;

}

int IdxDistance(int maxind, int zeropos, int *length, int dim) {
  int d,
    delta = 0,
    x = maxind,
    y = zeropos;
  for (d=0; d<dim; d++) {
    delta += abs( (x % length[d]) - (y % length[d])); // abs OK
    x /= length[d];
    y /= length[d];
  }
  return delta;
}


void do_BRmixed(cov_model *cov, gen_storage *s) {
  // to do: improve simulation speed by dynamic sizes
  assert(cov->key!=NULL);
  BR_storage *sBR = cov->SBR;
  int d,
     dim = cov->tsdim;
  pgs_storage *pgs = cov->Spgs;
  assert(pgs != NULL);
  bool changedIdx;
  double area = 1.0,
    step = P0(BR_MESHSIZE),
    invstepdim = intpow(step, -dim);

  //  printf("unif %f\n", UNIFORM_RANDOM);

  if ((changedIdx = pgs->currentthreshold == RF_NEGINF && sBR->idx != 0)) {
    sBR->idx = 0; // um zu gewaehrleisten,
  // dass insbesondere bei mehreren Simulationen wieder von vorne angefangen
  // wird
  } else if ((changedIdx = sBR->idx < sBR->maxidx &&
	      //pgs->currentthreshold < 0 !
	      pgs->currentthreshold >= sBR->thresholds[sBR->idx + 1])) {
    (sBR->idx)++;
  }

  assert(sBR->idx <= sBR->maxidx && sBR->maxidx < MAXSUB);
  cov_model  *key = sBR->sub[sBR->idx]; // !
  location_type //*loc = Loc(cov),
                *keyloc = Loc(key);
  double *lowerbounds = sBR->lowerbounds[sBR->idx],
    delta = sBR->radii[sBR->idx] + step;

  if (changedIdx) {    
    if (PL > 5)
      PRINTF("current level in BRmixed is %d", sBR->idx); 
    cov_model *prev = cov;
    while (prev->fieldreturn && !prev->origrf) {
      prev->rf = key->rf;
      prev = prev->calling;
      if (prev == NULL) break;
    }
 
    int *cumsum = pgs->own_grid_cumsum,
      *lglength = keyloc->length;

    cumsum[0] = 1;
    for (d=0; d<dim; d++) {
      pgs->own_grid_length[d] = keyloc->xgr[d][XLENGTH];
      pgs->own_grid_step[d] = keyloc->xgr[d][XSTEP];
      cumsum[d + 1] = cumsum[d] * lglength[d]; 
    }

    for (d=0; d<dim; d++) {
      sBR->suppmin[d] = sBR->locmin[d] - delta;
      sBR->suppmax[d] = sBR->locmax[d] + delta;
      area *= sBR->suppmax[d] - sBR->suppmin[d];
    }
    pgs->log_density = - log(area);
    cov->mpp.maxheights[0] = area;
    assert(cov->vdim2[0] == 1 && cov->vdim2[1] == 1);

    // MARCO lambda verstehe ich nicht -- wo/wie taucht lambda 
    // in den einzelnen Huetchen auf? 
    // @MARTIN: s. Zeile 957
    
  }

  if (PL > 5)
    PRINTF("idx=%d %d  %d zhou_n=%d %d %d\n", 
  	 sBR->idx, changedIdx, P0INT(BR_OPTIM), pgs->n_zhou_c, sBR->next_am_check,
  	 GLOBAL.br.deltaAM);
  
  if (P0INT(BR_OPTIM) == 2 &&  pgs->n_zhou_c >= sBR->next_am_check) {
    sBR->next_am_check += GLOBAL.br.deltaAM; 
    //    assert(false);
    OptimArea(cov, sBR->idx); 
    set_lowerbounds(cov);
  }

  res_type  // *res = cov->rf,
    *lgres = key->rf; // == cov->rf
  assert(cov->rf == key->rf);
  assert(keyloc->grid);

  int i, j, maxind, hatnumber,
    lgtotalpoints = keyloc->totalpoints,
    zeropos = sBR->zeropos[sBR->idx],
    vertnumber = P0INT(BR_VERTNUMBER);
  double maxval,
      u[MAXMPPDIM],
      *trend = sBR->trend[sBR->idx],     
    radius = sBR->radii[sBR->idx];


  for (d=0; d<dim; d++) {
    u[d] = sBR->suppmin[d] + UNIFORM_RANDOM*(sBR->suppmax[d] - sBR->suppmin[d]);
    pgs->supportmin[d] = u[d] - radius;
    pgs->supportmax[d] = u[d] + radius;    
    pgs->own_grid_start[d] = keyloc->xgr[d][XSTART] + u[d];
  }

  hatnumber=0;
  while(true) {
    //   printf("while unif %f\n", UNIFORM_RANDOM);
   
    DO(key, s);
    //   printf("while unif 2 %f\n", UNIFORM_RANDOM);
    maxval = RF_NEGINF;
    maxind = 0;
    for (i=0; i<lgtotalpoints; i++) {
      // MARCO: hier wird trend abgezogen; 
      lgres[i] -= trend[i];
      if (lgres[i] > maxval) {
	maxval = lgres[i];
	maxind = i;
      }
    }
    if (maxind == zeropos) {
      pgs->sq_zhou_c += invstepdim * invstepdim;
      pgs->sum_zhou_c += invstepdim;

    }

    //  printf("maxind = %d\n", maxind);

  

    //printf("here\n");
    if (P0INT(BR_OPTIM) == 2) {
      int  idxdist;
      double 
	uplusmaxval = maxval - lgres[zeropos] - log(UNIFORM_RANDOM);
      //  for(j=0; j<vertnumber; j++) printf("%f ", sBR->logvertnumber[j]); printf("uplus=%f\n", uplusmaxval);
    for (j=0; j<vertnumber; j++) {
	//	printf("here %d %d\n", j, uplusmaxval > sBR->logvertnumber[j]);
	if (uplusmaxval > sBR->logvertnumber[j]) {

	  //  printf("%f %f %d\n", uplusmaxval, sBR->logvertnumber[j], j);

	  // MARCO: ok? Warum wird die Schleife nach 1x abgebrochen?
	  // ('first' braucht es somit eigentlich nicht.
	  idxdist = IdxDistance(maxind, zeropos, keyloc->length, dim);
	  (sBR->countvector[j][idxdist])++;
	  

	  if (false && PL > 5) {
	    int dd;
	    for (dd=0; dd<75; dd++) {
	      int value = sBR->countvector[j][dd];
	      if (value < 0) BUG; // memory zugriffs-fehler ueber valgrind sichtbar
	      if (value<=9) PRINTF("%d", value);
	      else if (value<=50) {
		PRINTF(value<=20 ? "A" : value<=30 ? "B" : value<=40 ? "C" : "D"
		       );
	      }
	      else PRINTF("#");
	    }
	    PRINTF("j=%d\n", j);

	    //   if (maxind == zeropos) {
	    //     PRINTF("%d %d %d %d\n", maxind, zeropos, idxdist, j); assert(false);
	    //  }
	  }
	  //	  printf("idxdist=%d %d dist=%d val=%d len=%d %d %f %f %d %d hat=%d\n", maxind, zeropos, 
	  //	 idxdist, sBR->countvector[j][idxdist], 
	  //	 keyloc->length[0], keyloc->length[1],
	  //	 maxval, lowerbounds[maxind], sBR->idx,  pgs->n_zhou_c, hatnumber );
	  break;
	} // else assert(maxval <= lowerbounds[maxind]); // nicht korrekt.
      }
    } else {
      assert(UNIFORM_RANDOM > -1.0); // nur umvergleichbarkeit von optim 1 und 
      //                              optim 2 herzustellen 
    }

    //  printf("while unif 3 %f\n", UNIFORM_RANDOM); assert(false);

    //  assert(maxind != zeropos);
    if (maxind == zeropos && PL > 5) {
      PRINTF("zeropos maxval=%f lower=%f hat=%d zhou_n=%d optim=%d\n",
	     maxval, lowerbounds[maxind], hatnumber, pgs->n_zhou_c, P0INT(BR_OPTIM));
      //assert(false);
    }

    // MARCO, dann es sein, dass die folgende Bedingung erfuellt ist,
    // aber uplusmaxval > sBR->logvertnumber[j] nicht erfuellt ist?
    if (maxval > lowerbounds[maxind]) break;    
    hatnumber++;
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

// nd thres: -6.697333 -6.000000 1.000000 63241312 63241312 RPbrmixed
//end thres: -5.326407 -4.756986 0.700000 63241312 63241312 RPbrmixed


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
  
  range->min[BR_OPTIMMAX] = 1;
  range->max[BR_OPTIMMAX] = RF_INF;
  range->pmin[BR_OPTIMMAX] = 1000;
  range->pmax[BR_OPTIMMAX] = 1000000000;
  range->openmin[BR_OPTIMMAX] = false;
  range->openmax[BR_OPTIMMAX] = true;
  
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
  double *newx= NULL,
    centreloc[MAXMPPDIM], 
    minloc[MAXMPPDIM], 
    maxloc[MAXMPPDIM];

  //  SERR("Brown Resnick not available in version 3.0.34");

  ASSERT_NEWMODEL_NULL;
  if (cov->role != ROLE_BROWNRESNICK) BUG;
  
  assert(isBRuserProcess(cov));  
  
  if (loc->Time || (loc->grid && loc->caniso != NULL)) {
    Transform2NoGrid(cov, false, GRIDEXPAND_AVOID);
    SetLoc2NewLoc(sub, Loc(cov));
  }
  
  loc = Loc(cov);
  grid = loc->grid;
  
  model_intern = (cov->nr == BRORIGINAL_USER) ? BRORIGINAL_INTERN
               : (cov->nr == BRMIXED_USER) ? BRMIXED_INTERN
               : (cov->nr == BRSHIFTED_USER) ? BRSHIFTED_INTERN
               : BRORIGINAL_USER;
	       
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if (cov->stor == NULL) NEW_STORAGE(stor, STORAGE, gen_storage);
  
  GetDiameter(loc, minloc, maxloc, centreloc);
  newxlen = loc->lx;
  if ((newx = (double*) MALLOC(dim * newxlen * sizeof(double))) == NULL) {
     GERR("Memory allocation failed.\n"); 
  } 
  if (grid) {
    for (d=0; d<dim; d++) {
      newx[3 * d + XSTART] = loc->xgr[d][XSTART] - centreloc[d] 
           + (((int) loc->xgr[d][XLENGTH]) % 2 == 0) * loc->xgr[d][XSTEP]/2;
      newx[3 * d + XSTEP] = loc->xgr[d][XSTEP];
      newx[3 * d + XLENGTH] = loc->xgr[d][XLENGTH];
   }
  } else {
    for (i=0; i<loc->lx; i++)
      for(d=0; d<dim; d++) newx[i*dim+d] = loc->x[i*dim+d] - centreloc[d];
  }

  if ((err = loc_set(newx, NULL, dim, dim, newxlen, false, loc->grid,
                     loc->distances,  &(cov->ownloc))) > NOERROR)
    goto ErrorHandling;
  SetLoc2NewLoc(sub, cov->ownloc);

  if ((err=covcpy(&(cov->key), sub)) > NOERROR) goto ErrorHandling;
  
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
    kdefault(cov->key, BR_OPTIMMAX, P0INT(BR_OPTIMMAX));
    kdefault(cov->key, BR_LAMBDA, P0(BR_LAMBDA));
    if (!PisNULL(BR_OPTIMAREA)) {
      if ( cov->nrow[BR_OPTIMAREA] % 2 == 0 ||
	   cov->ncol[BR_OPTIMAREA] % 2 == 0 ) {
	SERR("number of rows and columns of areamat need to be odd")
      }      
      PARAMALLOC(cov->key, BR_OPTIMAREA, cov->nrow[BR_OPTIMAREA],
                 cov->ncol[BR_OPTIMAREA]);

      //PMI(cov);
      
      PCOPY(cov->key, cov, BR_OPTIMAREA);
    
    }
  }
  
  cov->key->calling = cov;
  if ((err = CHECK(cov->key, dim, dim, PointShapeType, cov->domown, 
                   cov->isoown, 1, ROLE_BROWNRESNICK)) == NOERROR) {
    // printf("Kopie erstellen\n");
    if ((err = STRUCT(cov->key, NULL)) <= NOERROR) {
      err = CHECK(cov->key, dim, dim, PointShapeType, cov->domown, cov->isoown,
		  1, ROLE_BROWNRESNICK);
    }
  }

  assert(cov->key->calling == cov);
  
   
  ErrorHandling:
     if (newx != NULL) free(newx);
     newx = NULL;
     return err;
  
}



// new
int structBRintern(cov_model *cov, cov_model **newmodel) {
  cov_model *sub = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  location_type *loc = Loc(cov);
  int i, d, j, err,  *length,
    dim = sub->tsdim, 
    totaldim = loc->totalpoints*dim,
    zeropos = 0, // default, mostly overwritten
    newxlen = 3; // default (grid)
  bool grid = loc->grid;
  double 
    //\00covval,
    norm, step,
    mindist, dist;
  BR_storage *sBR = NULL;
  //printf("xxxxxxxxxxx\n\n\nxxxxxxxxxxxxxx\n");
  
  ASSERT_NEWMODEL_NULL;
  if (cov->role != ROLE_BROWNRESNICK) BUG;

  assert(isPointShape(cov));
    
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if (cov->stor == NULL) NEW_STORAGE(stor, STORAGE, gen_storage);
  NEW_STORAGE(SBR, BR, BR_storage);
  sBR = cov->SBR;

  
  if (cov->sub[MPP_TCF] != NULL) {
    if ((err = STRUCT(sub, &(cov->key))) > NOERROR) goto ErrorHandling;
    cov->key->calling = cov;
  } else if ((err=covcpy(&(cov->key), sub)) > NOERROR) goto ErrorHandling;
 
  if ((err = CHECK(cov->key, dim, dim, NegDefType, cov->domown, SYMMETRIC, 
		   1, ROLE_COV)) != NOERROR) goto ErrorHandling;  

  if ((err = covcpy(&(sBR->submodel), cov->key)) != NOERROR) goto ErrorHandling;
  if ((err = CHECK(sBR->submodel, 1, 1, NegDefType, XONLY, 
		   ISOTROPIC, 1, ROLE_COV)) != NOERROR) goto ErrorHandling;
  
  if ((err = newmodel_covcpy(&(sBR->vario), VARIOGRAM_CALL, cov->key))!=NOERROR)
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
        if (loc->xgr[d][XSTEP] < mindist)
	  mindist = loc->xgr[d][XSTEP];
    } else {
      int maxtotal = totaldim;
      if (maxtotal > 1000) {
	warning("minimal distance is only estimated");
	maxtotal = 1000; // to do
      }
      for (i=0; i<maxtotal; i+=dim)
        for (j=i+dim; j<maxtotal; )
 	  for (d=0; d<dim; d++, j++) {
	    dist = fabs(loc->x[i+d] - loc->x[j]);
	    if (dist > 1e-15 && dist < mindist) mindist = dist;
	  }
    }
    if (mindist < step) {     
      PRINTF("Warning! meshsize is larger than resolution of simulation! meshsize is automatically decreased to %f.\n", step);
      P(BR_MESHSIZE)[0] = step = mindist;
    }
    
    if (!PisNULL(BR_OPTIMAREA)) {
      sBR->minradius = 
	floor(fmax(cov->nrow[BR_OPTIMAREA], cov->ncol[BR_OPTIMAREA])*0.5) * step;
    }

    double alpha = -1,  yy, C0, gammamin, min_radius,
      //factor = 0.1,
      xx=step * 1e-6;
    COV(ZERO, sBR->submodel, &C0);
    COV(&xx, sBR->submodel, &yy);
    alpha = log(C0 - yy) / log(xx);
    // printf("alpha = %f %e %e\n", alpha, xx, C0-yy);
    if (alpha > 2.0) alpha = 2.0;
    gammamin = 4.0 - 1.5 * alpha;
    INVERSE(&gammamin, sBR->submodel, &min_radius);
 
    xx = gammamin + P0(BR_VARIOBOUND);
    INVERSE(&xx, sBR->submodel, sBR->radii);
    sBR->radii[0] = ceil(sBR->radii[0] / step) * step;
    if (sBR->minradius > sBR->radii[0]) sBR->radii[0] = sBR->minradius;
    sBR->thresholds[0] = RF_INF;
 
    newxlen = 3;
    grid = true;
    double smallest_grid_lengths = 3.0, //kleinestes Feld, dass simuliert werden 
      //                       soll, in step Einheiten in 1 koordinaten richtung
      min_radius_step = floor(min_radius / step),
      onesided  = (double) (smallest_grid_lengths - 1) * 0.5;
    //   printf("onesided %f %f %f\n", onesided, smallest_grid_lengths, min_radius_step);
    if (onesided < min_radius_step) onesided = min_radius_step;  
     double corrected_minradius = 
       sBR->minradius + GLOBAL.br.corr_factor * (sBR->radii[0] - sBR->minradius);

    if (onesided * step < corrected_minradius) 
	onesided = floor(corrected_minradius / step);  

    // MARCO: @MARTIN: um sicherzugehen, dass man den durch areamat abgedeckten
     // Bereich komplett simuliert  

    double
      min_reduction_factor = 2.0,
      max_reduction_factor = 3.0,
      rel_steps = (sBR->radii[0] / step) / onesided,
      reduction_factor = pow(rel_steps, (double) dim / ((double) MAXSUB -1.0));
    if (reduction_factor < min_reduction_factor) 
      reduction_factor = min_reduction_factor;
    else if (reduction_factor > max_reduction_factor)
      reduction_factor = max_reduction_factor;
    sBR->maxidx = (int) ((double) dim * log(rel_steps) / log(reduction_factor));
    if (sBR->maxidx < 0) sBR->maxidx = 0; // to do!!
    else if (sBR->maxidx >= MAXSUB) sBR->maxidx = MAXSUB - 1;
  
  
//  printf("maxidx %d  %f r=%f step=%f red=%f %f relste=%f red=%f %f %f\n", sBR->maxidx, smallest_grid_lengths, sBR->radii[0], step, reduction_factor, min_reduction_factor, rel_steps, pow(rel_steps, (double) dim / ((double) MAXSUB - 1.0)), (sBR->radii[0] / step) / smallest_grid_lengths, sBR->radii[0] / step);
  
 
 
    for (i=1; i<=sBR->maxidx; i++) {
      // to do gleichmaessige reduktion ist vermutlich nicht optimal
      // mitte vermutlich groessere Schritte und dafuer kleiner runter
      // bis zum minimum.
      // to do optimale groesse fuer reduction_factor z.Zt. in [2,3]?!

      //printf("i=%d %d %f\n", i, sBR->maxidx, sBR->radii[i-1]); assert(false);

       sBR->radii[i] = round(sBR->radii[i-1] * pow(2.0, -1.0/dim) / step) * step;
       
       //   printf("radii %f \n",    sBR->radii[i] /step);

       COV(sBR->radii + i, sBR->submodel, sBR->thresholds + i);
    }
    
    //xassert(false);
    

   if ((sBR->newx = (double*) MALLOC(newxlen*dim*sizeof(double))) == NULL) {
      err = ERRORMEMORYALLOCATION; goto ErrorHandling;
   }
   for (i=0; i<=sBR->maxidx; i++) {
      if ((err = covcpy(sBR->sub + i, cov->key)) != NOERROR) goto ErrorHandling;
      for (d=0; d<dim; d++) {
	sBR->newx[3*d+XSTART] = -sBR->radii[i];
	sBR->newx[3*d+XLENGTH] = 2 * ((int) round(sBR->radii[i] / step)) + 1;
	sBR->newx[3*d+XSTEP] = step;


//	printf("dim=%d %f %f %f radii=%f step=%f\n", dim, sBR->newx[3*d+XSTART],
//	       sBR->newx[3*d+XLENGTH], sBR->newx[3*d+XSTEP],
//	       sBR->radii[i], step
//	       );
	
	
      }

      // assert(false);

      err = loc_set(sBR->newx, NULL, dim, dim, newxlen, false, grid, 
                  false, &(sBR->sub[i]->ownloc));
      length = Loc(sBR->sub[i])->length;
      for (zeropos = 0, d = dim; d > 0; d--)
	zeropos = zeropos * length[d-1] + (int) ceil(length[d-1] * 0.5) - 1;
      sBR->zeropos[i] = zeropos;

      //printf("i=%d\n", i);assert(false);
      if ((err = CHECK(sBR->sub[i], dim, dim, ProcessType, cov->domown, 
		       cov->isoown, 1, ROLE_GAUSS)) == NOERROR) {
	if ((err = STRUCT(sBR->sub[i], NULL)) <= NOERROR) {
	  err = CHECK(sBR->sub[i], dim, dim, ProcessType, cov->domown,
		      cov->isoown, 1, ROLE_GAUSS); 
	}
      }
      if (err > NOERROR) goto ErrorHandling;
      assert(sBR->sub[i]->calling == cov);
    }

   goto ErrorHandling; // ueberspringe immer den shifted & original teil !

  } else { //  END BRMIXED;  START SHIFTED    
    assert(cov->nr == BRSHIFTED_INTERN);
  }
 

  if (sBR->newx == NULL) {
    if ((err = loc_set(grid ? * loc->xgr : loc->x, NULL, dim, dim,
		       grid ? 3 : loc->totalpoints, false, grid,
		       loc->distances, &(cov->key->ownloc))) != NOERROR)
      goto ErrorHandling;
  } else {
    if ((err = loc_set(sBR->newx, NULL, dim, dim, newxlen, false, grid, false,
		       &(cov->key->ownloc))) != NOERROR) goto ErrorHandling;
  }

 
  // ONLY SHIFTED AND ORIGINAL
  if (grid) {
    length = Loc(cov->key)->length;
    for (d=dim; d>0; d--)
      zeropos = zeropos*length[d-1] + (int) ceil(length[d-1] * 0.5) - 1;
  }
  sBR->zeropos[0] = zeropos;
   
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
  if (err != NOERROR) BR_DELETE(&(cov->SBR));    
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
    Transform2NoGrid(cov, false, GRIDEXPAND_AVOID);
    SetLoc2NewLoc(next, Loc(cov));
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
      if (!PARAMisNULL(next, DVAR) )  newscale /= sqrt(PARAM0(next, DVAR));
      if (factor != 1.0) {
	newscale *= factor;
	factor = 1.0;
      }
      //      print("here3\n");
      return ERRORNOTPROGRAMMEDYET;
      

      kdefault(calling, DSCALE, newscale);
      next = next->sub[0]; // ok
      //     print("here4\n");
    }

    if (cov->sub[MPP_TCF] != NULL) {
      // print("here4\n");
      return ERRORNOTPROGRAMMEDYET;
    } 

    // PMI(atom->sub[ATOMSHAPE]); PMI(*shape); assert(false);
    // print("here5\n");
    
    if (next->nr == BROWNIAN && PARAM0(next, BROWN_ALPHA) == 2.0) {
      addModel(&(cov->key), GAUSS);   // ?? 
      if (factor != 1.0) {
	addModel(&(cov->key), DOLLAR);
	kdefault(cov->key, DSCALE, factor);      
	// (newmod->calling = atom;
      }
    } else {
      //COV_DELETE(atom->sub + ATOMSHAPE);
      SERR("Smith process with BrownResnick tcf only possible for fractal Brownian motion with alpha=2");
    }
  } else if (cov->role == ROLE_BROWNRESNICK) {
    if (next->role == ROLE_BROWNRESNICK) {
      SERR1("submodel of '%s' must be a covariance model or tcf", 
	    NICK(cov));
    } else {
      role = isNegDef(next) ? ROLE_COV : ROLE_UNDEFINED;
	
      ASSERT_ROLE_DEFINED(next);  

      // PMI(next);
      
      if (((err = covcpy(&(cov->key), next)) != NOERROR)
         || ((err = CHECK(cov->key, dim, dim, NegDefType, XONLY,
			  SYMMETRIC, 1, role)) != NOERROR)) {
	//		XERR(err); printf("xx"); 
	//assert(false);
	return err;
      }

      GetDiameter(loc, minloc, maxloc, centreloc);
      for (d=0; d<MAXMPPDIM; d++) maxdist[d] = 0.5*(maxloc[d] - minloc[d]);
   
      cov_model *K = NULL;

      //     PMI(cov->key);
      //      PMI(K);

      if ((err = newmodel_covcpy(&K, VARIOGRAM_CALL, cov->key, maxdist, NULL,
				 NULL, dim, dim, 1, 0, false, false, false))
	  != NOERROR) return err;
      if ((err = alloc_cov(K, dim, 1, 1)) != NOERROR) return err;
      if (K->sub[0] != NULL) SetLoc2NewLoc(K->sub[0], Loc(K));
      // PMI()
      Variogram(NULL, K, &maxcov);
      COV_DELETE(&K);

      //printf("maxcov: %f\n", maxcov);
      
      if (isPosDef(next) || maxcov <= 4.0) {
	meth = BRORIGINAL_USER;  
      } else if (!next->logspeed || next->logspeed <= 4.0 || maxcov <= 10.0) {
	meth = BRSHIFTED_USER;
      } else {
	meth = BRMIXED_USER;
      }

      addModel(&(cov->key), meth, cov);
      cov_model *key = cov->key;
      key->prevloc = loc;
      
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
  
  // PMI(atom); assert(false);

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
