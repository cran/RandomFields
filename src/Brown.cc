/* 
 Authors
 Marco Oesting,
 Martin Schlather, schlather@math.uni-mannheim.de

 simulation of Brown-Resnick processes

 Copyright (C) 2009 -- 2010 Martin Schlather 
 Copyright (C) 2011 -- 2013 Marco Oesting & Martin Schlather 

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
#define BR_LB_OPTIMAREA LAST_MAXSTABLE + 2
#define BR_VERTNUMBER LAST_MAXSTABLE + 3
#define BR_OPTIM LAST_MAXSTABLE + 4
#define BR_OPTIMTOL LAST_MAXSTABLE + 5
#define BR_OPTIMMAX LAST_MAXSTABLE + 6
#define BR_LAMBDA LAST_MAXSTABLE + 7
#define BR_OPTIMAREA LAST_MAXSTABLE + 8
#define BR_VARIOBOUND LAST_MAXSTABLE + 9

//////////////////////////////////////////////////////////////////////
// Brown Resnick

int checkBrownResnick(cov_model *cov) {
  cov_model  
    *key = cov->key,
    *sub = key != NULL ? key :
            cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  int err, role, isoprev,
      dim = cov->tsdim;
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
  
  isoprev = role == ROLE_COV ? SYMMETRIC : NO_ROTAT_INV;
  
  if ((err = CHECK(sub, dim, dim, 
		   type,//Martin: changed 27.11.13 CovList[cov->nr].Type,
		   XONLY, isoprev, 1, role)) != NOERROR) {
      return err;
  }
  
  setbackward(cov, sub);
  
  return NOERROR;  
}

int init_BRorig(cov_model *cov, storage *s){
 
  cov_model *key = cov->key;
  location_type *keyloc = NULL;
  int err, d, dim;
  bool keygrid;
  BR_storage *sBR = NULL;
  
  if (cov->role == ROLE_BROWNRESNICK) {
    
    if (key != NULL) {
      dim = cov->tsdim;
      if ((err = alloc_cov(cov, dim, 1)) != NOERROR) return err;
      pgs_storage *pgs = cov->Spgs;
      for (d=0; d<dim; d++) {
        pgs->supportmin[d] = RF_NEGINF; // 4 * for debugging...
        pgs->supportmax[d] = RF_INF;
        pgs->supportcentre[d] = RF_NAN;
      }
      pgs->log_density = 0.0;
      
      keyloc = Loc(key);
      keygrid = keyloc->grid;
      
      ((int*) cov->key->p[LOG_GAUSS])[0] = false;  
      key->simu.active = true;
      key->simu.expected_number_simu = cov->simu.expected_number_simu;
      
      assert(cov->mpp.moments >= 1);
      if ((err = INIT(key, 1, s)) != NOERROR) return err;  
      
      cov->loggiven = true;
      
      assert(key->nr == GAUSSPROC);
      assert(key->mpp.moments >= 1);
      cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0;
      cov->mpp.M[1] = cov->mpp.Mplus[1] = 1.0;
      cov->mpp.maxheight = GLOBAL.extreme.standardmax;
      cov->mpp.log_zhou_c = 0.0;
      
      sBR = cov->SBR;
      sBR->trendlen = 1;
      if ((cov->SBR->trend = (double**) MALLOC(sizeof(double*)))==NULL) {
        err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }    
      if ((cov->SBR->trend[0] = 
           (double*) MALLOC(keyloc->totalpoints*sizeof(double)) ) == NULL) {
        err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }

      if ((err = loc_set(keygrid ? keyloc->xgr[0] : keyloc->x, NULL, NULL, dim,
	             dim, keygrid ? 3 : keyloc->totalpoints, 0, false, keygrid,
		     keyloc->distances, &Loc(cov->SBR->vario))) > NOERROR)
	return err;
      if (sBR->vario->sub[0] != NULL) 
	SetLoc2NewLoc(sBR->vario->sub[0], Loc(sBR->vario));
      Variogram(NULL, sBR->vario, sBR->trend[0]);
      
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

void do_BRorig(cov_model *cov, storage *s) {
  int i, 
     totalpoints = Loc(cov)->totalpoints;
  assert(cov->origrf);
  res_type *res = cov->rf;
  assert(cov->key != NULL);
  cov_model *key = cov->key;
  res_type *lgres = key->rf;
  BR_storage *sBR = cov->SBR;
  assert(sBR != NULL);
  assert(sBR->trend != NULL);
  double **trend = sBR->trend;
  
  DO(key, s); // nicht gatternr
  for (i=0; i<totalpoints; i++) {
    res[i] = lgres[i] - lgres[sBR->zeropos] - trend[0][i];
  }
}

int init_BRshifted(cov_model *cov, storage *s) {
  cov_model *key=cov->key;
  location_type *keyloc = NULL;
  int d, j, err, dim,
      shiftedloclen, keytotal,
      trendlenmax, trendlenneeded;
  bool keygrid;
  BR_storage *sBR = NULL;
  
  if (cov->role == ROLE_BROWNRESNICK) {
    
    if (key != NULL) {
      dim = cov->tsdim; 
      if ((err = alloc_cov(cov, dim, 1)) != NOERROR) return err;
      pgs_storage *pgs = cov->Spgs;
      for (d=0; d<dim; d++) {
	pgs->supportmin[d] = RF_NEGINF; // 4 * for debugging;
	pgs->supportmax[d] = RF_INF;
	pgs->supportcentre[d] = RF_NAN;
      }
      pgs->log_density = 0.0;
      
      keyloc = Loc(key);
      keygrid = keyloc->grid;
      keytotal = keyloc->totalpoints;
      
      ((int*) cov->key->p[LOG_GAUSS])[0] = false;
      key->simu.active = true;
      key->simu.expected_number_simu = cov->simu.expected_number_simu; 
      assert(cov->mpp.moments >= 1);
      if ((err = INIT(key, 1, s)) != NOERROR) return err;  
      
      cov->loggiven = true;
      
      assert(key->nr == GAUSSPROC);
      assert(key->mpp.moments >= 1);
      cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0;
      cov->mpp.M[1] = cov->mpp.Mplus[1] = 1.0;
      cov->mpp.maxheight = GLOBAL.extreme.standardmax;
      cov->mpp.log_zhou_c = 0.0;
      
      sBR = cov->SBR;
      
      shiftedloclen = keygrid ? 3 : keytotal;
      if ((cov->SBR->shiftedloc = (double*) 
            MALLOC(dim*shiftedloclen*sizeof(double))) == NULL) {
        err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }

      trendlenmax = (int) ceil((double) GLOBAL.br.BRmaxmem / keytotal);
      trendlenneeded = MIN(keytotal, cov->simu.expected_number_simu);
      sBR->trendlen = MIN(trendlenmax, trendlenneeded);
   
      //PMI(cov->calling);
      //      printf("len=%d %d %d mem=%d %d exp=%d %d\n", sBR->trendlen, trendlenmax, trendlenneeded, GLOBAL.br.BRmaxmem, keytotal, cov->simu.expected_number_simu, GLOBAL.general.expected_number_simu);



      assert(sBR->trendlen > 0);
      
      sBR->memcounter = 0;
      if ((cov->SBR->loc2mem=(int*) MALLOC(sizeof(int)*keytotal))==NULL) {
        err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      for (j=0; j<keytotal; j++) sBR->loc2mem[j] = -1;    
    
      if ((cov->SBR->mem2loc=(int*) MALLOC(sizeof(int)*sBR->trendlen))==NULL) {
        err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      if ((cov->SBR->trend=(double**) MALLOC(sizeof(double*)*sBR->trendlen))
           ==NULL) {
       err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      for (j=0; j<sBR->trendlen; j++) {
        sBR->mem2loc[j] = -1;
        if ((cov->SBR->trend[j] = 
           (double*) MALLOC(keytotal*sizeof(double)))==NULL) {
          err = ERRORMEMORYALLOCATION; goto ErrorHandling;
        }
      }
      
      if ((err = loc_set(keygrid ? keyloc->xgr[0] : keyloc->x, NULL, NULL, dim,
	                 dim, keygrid ? 3 : keytotal, 0, false, keygrid,
	                 keyloc->distances, &Loc(cov->SBR->vario))) > NOERROR)
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

void do_BRshifted(cov_model *cov, storage *s) {
  BR_storage *sBR = cov->SBR;
  cov_model *key = cov->key;
  assert(cov->key != NULL);
  
  location_type *keyloc = Loc(key);
  int d, i, k, trendindex,
     zeropos, zeroposMdim,
     dim = cov->tsdim,
     keytotal = keyloc->totalpoints,
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
  zeropos = (int) floor(UNIFORM_RANDOM * keytotal);
  
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
      indextrafo(zeropos, keyloc->length, dim, locindex);
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
    Variogram(NULL, sBR->vario, sBR->trend[trendindex]);
    
    mem2loc[trendindex] = zeropos;
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
  kdefault(cov, BR_LB_OPTIMAREA, bp->BR_LB_optim_area);
  kdefault(cov, BR_VERTNUMBER, bp->BRvertnumber);
  kdefault(cov, BR_OPTIM, bp->BRoptim);
  kdefault(cov, BR_OPTIMTOL, bp->BRoptimtol);
  kdefault(cov, BR_OPTIMMAX, bp->BRoptimmaxpoints);
  kdefault(cov, BR_VARIOBOUND, bp->variobound);
  if ( (cov->nr == BRMIXED_USER) && (cov->key == NULL) && 
       (((int*) cov->p[BR_OPTIM])[0] > 0) && 
       ((cov->p[BR_LAMBDA] != NULL) || (cov->p[BR_OPTIMAREA] != NULL)) ) {
    ERR("BR optimization requires lambda and areamat to be unset");
  }
  kdefault(cov, BR_LAMBDA, 1.0);
  if (cov->p[BR_OPTIMAREA] == NULL) //necessary as areamat might not be scalar
    kdefault(cov, BR_OPTIMAREA, 1.0);
  
  if ((err = checkBrownResnick(cov)) != NOERROR) return err;
  
  if ((err = checkkappas(cov, true)) != NOERROR) return err;
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
    return ERRORDIM;

  return NOERROR; 
}

void kappaBRmixed(int i, cov_model *cov, int *nr, int *nc){
  // i nummer des parameters

  int dim = cov->tsdim;
  
  switch(i) {
    
    case GEV_XI: case GEV_MU: case GEV_S: case BR_MESHSIZE:
    case BR_LB_OPTIMAREA: case BR_VERTNUMBER: case BR_OPTIM:
    case BR_OPTIMTOL: case BR_OPTIMMAX: case BR_LAMBDA: case BR_VARIOBOUND:
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

int BRoptim(cov_model *cov, storage *s) {
  cov_model *key = cov->key;
  location_type *gridloc = Loc(key);
  res_type *lgres = key->rf;
  BR_storage *sBR = cov->SBR;
  int d, i, j, k, l, maxind, 
      **countmatrix = NULL,
      totalcount = 0, limit, addnumber,
      zeroxpos, zeroypos, cellcounter,
      *len = gridloc->length,
      err = NOERROR,
      dim = cov->tsdim,
      totalpoints = gridloc->totalpoints,
      radius = (int) floor(gridloc->length[0] / 2.0),
      vertnumber = ((int*) cov->p[BR_VERTNUMBER])[0],
      radiusP1 = radius+1,
      optimmax = ((int*) cov->p[BR_OPTIMMAX])[0],
      optim = ((int*) cov->p[BR_OPTIM])[0],
      zeropos = sBR->zeropos;
  double lambda=0.0, newlambda=0.0,
         dummy, error, errortol = cov->p[BR_OPTIMTOL][0],
         maxerrorbound, errorbound, errorboundtmp, 
	 **areamatrix = NULL,
         step = cov->p[BR_MESHSIZE][0],
	 stepdim= intpow(step, dim),
	 *trend = sBR->trend[0],
	 xmin = cov->p[BR_LB_OPTIMAREA][0],
	 summand, maxval, nullval;
	 
  switch(optim) {
    case 1:
      totalcount = 0;
      for (k=0; k<optimmax; k++) {
   	DO(key, s);
	maxind = 0;
	maxval = RF_NEGINF;
	for (i=0; i<totalpoints; i++) {
	  if (lgres[i] - lgres[zeropos] - trend[i] > maxval) {
	    maxval = lgres[i] - lgres[zeropos] - trend[i];
	    maxind = i;
	  }
	}
	if (maxind == zeropos) totalcount++;
      }
      cov->p[BR_LAMBDA][0] = totalcount / (optimmax * stepdim);
      //   printf("lambda: %f\n", cov->p[BR_LAMBDA][0]); 
      return NOERROR;	
    case 2:
      assert(dim <= 2);
	 
      if ((countmatrix = (int**) MALLOC(totalpoints*sizeof(int*))) == NULL) {
        err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      for (i=0; i<totalpoints; i++) {
        if ((countmatrix[i]=(int*) MALLOC(vertnumber*sizeof(int))) == NULL) {
          err = ERRORMEMORYALLOCATION; goto ErrorHandling;
        }
        for (j=0; j<vertnumber; j++) countmatrix[i][j] = 0;
      }
  
      if ((areamatrix = (double**) MALLOC(radiusP1*sizeof(double*))) == NULL) {
        err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
  
      for (i=0; i<radiusP1; i++) {
        if ( (areamatrix[i]=(double*) MALLOC(vertnumber*sizeof(double)))
	                                                             == NULL) {
          err = ERRORMEMORYALLOCATION; goto ErrorHandling;
        } 
        for (j=0; j<vertnumber; j++) areamatrix[i][j] = 0.0;
      }
  
      if (cov->p[BR_OPTIMAREA] != NULL) free(cov->p[BR_OPTIMAREA]);
      cov->p[BR_OPTIMAREA] = NULL;
      if ((cov->p[BR_OPTIMAREA] = (double*) MALLOC(totalpoints*sizeof(double)))
                                                                    == NULL) { 
        err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      cov->nrow[BR_OPTIMAREA] = (int) (dim==1) ? 1 : len[1];
      cov->ncol[BR_OPTIMAREA] = len[0]; 
  
      for (k=0; k<optimmax; k++) {
        summand = xmin - log(UNIFORM_RANDOM);
        DO(key, s);
        maxval = RF_NEGINF;
        maxind = 0;
        nullval = lgres[zeropos];
        for (i=0; i<totalpoints; i++) {
          lgres[i] += summand - nullval - trend[i];
          if (lgres[i] > maxval) {
	    maxval = lgres[i];
	    maxind = i;
          }
        }
        for (j=0; j<vertnumber; j++) {
          if (maxval > xmin - log((double) (j+1)/vertnumber)) {
	    countmatrix[maxind][j]++;
	    break;
          }
        }
      }
  
      for (d=0; d<radiusP1; d++) { //symmetrize
	addnumber = (d==0) ? 1 : ((dim==1) ? 2 : 4*d);
	for (j=0; j<vertnumber; j++) { 
          totalcount = 0;
          limit = (dim==1) ? 0 : (int) fmin(d,len[1]);
          for (i=-limit; i<=limit; i++) {
            if (i == d || i == -d) {
	      totalcount += countmatrix[zeropos + i*len[0]][j]; 
            } else {
	      int idx = zeropos + i*len[0] + (int) (d-abs(i));
	      totalcount += countmatrix[idx][j];
	      totalcount += countmatrix[idx][j];
	    }
          }
          areamatrix[d][j] = totalcount/(optimmax*addnumber*stepdim);
        }
      }
      for (j=0; j<vertnumber; j++) lambda += areamatrix[0][j];
      for (j=0; j<vertnumber; j++) areamatrix[0][j] = lambda / vertnumber;
  
      for (d=0; d<radius; d++) { //monotonic in each column
        for (j=0; j<vertnumber; j++) {
          for (k=j+1; k<vertnumber; k++) {
	    if (areamatrix[d][j] < areamatrix[d][k]) {
	       dummy = areamatrix[d][j];
	       areamatrix[d][j] = areamatrix[d][k];
	       areamatrix[d][k] = dummy;
	    }
          }
          areamatrix[d][j] = fmin(areamatrix[d][j], lambda/vertnumber); //cutoff
        }
      }
  
      for (j=0; j<vertnumber; j++) { //monotonic in each row
        for (d=0; d<radiusP1; d++) {
          for (k=d+1; k<radiusP1; k++) {
	    if (areamatrix[d][j] < areamatrix[k][j]) {
	      dummy = areamatrix[d][j];
	      areamatrix[d][j] = areamatrix[k][j];
	      areamatrix[k][j] = dummy;
	    }
          }
        }
      }
  
      newlambda = 0.0; //correction for bias because of cutoff
      for (d=0; d<radiusP1; d++) {
        addnumber = (d==0) ? 1 : ((dim==1) ? 2 : 4*d);  
        for (j=0; j<vertnumber; j++) {
          newlambda += areamatrix[d][j] * addnumber;      
        }
      }
      for (d=0; d<radiusP1; d++) {
        for (j=0; j<vertnumber; j++) {
          areamatrix[d][j] *= 1/(stepdim*newlambda);
        }
      }
      lambda *= 1/(stepdim*newlambda);
  
      maxerrorbound = 0.0;
      for (d=0; d<radiusP1; d++) { //arematrix => error matrix
        for (j=0; j<vertnumber; j++) {
          areamatrix[d][j] = lambda/vertnumber - areamatrix[d][j];
          maxerrorbound = fmax(maxerrorbound, areamatrix[d][j]);
        }
      }

      error = 0.0;
      errorboundtmp = errorbound = 0.0;  
      while (errorboundtmp < maxerrorbound && error < errortol) {
	errorbound = errorboundtmp; 
        errorboundtmp = maxerrorbound;
        for (d=0; d<radiusP1; d++) {
          for (j=0; j<vertnumber; j++) {
	    if (areamatrix[d][j] > errorbound) {
	      errorboundtmp = fmin(errorboundtmp, areamatrix[d][j]); 
	    }
          }      
        }
        error = 0.0;
        cellcounter = 0;
        for (d=0; d<radiusP1; d++) {
          addnumber = (d==0) ? 1 : ((dim==1) ? 2 : 4*d);
          for (j=0; j<vertnumber; j++) {
	    if (areamatrix[d][j] <= errorboundtmp + 1e-6) {
	      error += addnumber*areamatrix[d][j]*stepdim*exp(-xmin);
	      cellcounter += addnumber;
	    }
          }
        }
        error = (1-exp(-error)) /
                  (1-exp(-cellcounter*lambda/vertnumber*stepdim*exp(-xmin)));
      }
  
      zeroxpos = floor(len[0] / 2.0);
      zeroypos = (dim==1) ? 0 : floor(len[1] / 2.0);
      newlambda = 0.0;
      cellcounter = 0;
      //printf("optimarea: (%d, %d)\n", cov->nrow[BR_OPTIMAREA], cov->ncol[BR_OPTIMAREA]);
      for (l=i=0; i<cov->nrow[BR_OPTIMAREA]; i++) {
        for (k=0; k<cov->ncol[BR_OPTIMAREA]; k++, l++) {
          d = (int) abs(k-zeroxpos) + (int) abs(i-zeroypos);
          cov->p[BR_OPTIMAREA][l] = 0.0;
          for (j=0; j<vertnumber; j++) {
            if (d<radiusP1 && (areamatrix[d][j] <= errorbound + 1e-6)) {
	      cov->p[BR_OPTIMAREA][l] += 1.0/vertnumber;
	      newlambda += lambda - areamatrix[d][j]*vertnumber;
	      cellcounter++;
	    }
          }
          // if (cov->p[BR_OPTIMAREA][l] > 0) printf("%4.2f ", cov->p[BR_OPTIMAREA][l]);
        }
        //printf("\n");
      }
      cov->p[BR_LAMBDA][0] = newlambda / cellcounter;
      //printf("lambda: %f\n", cov->p[BR_LAMBDA][0]);
  
      ErrorHandling:
      if (countmatrix != NULL) {
        for (i=0; i<totalpoints; i++) {
          if (countmatrix[i] != NULL) free(countmatrix[i]);
          countmatrix[i] = NULL;
	}
        free(countmatrix);
      }
      countmatrix=NULL;

      if (areamatrix != NULL) {
        for (i=0; i<radiusP1; i++) {
          if (areamatrix[i] != NULL) free(areamatrix[i]);
          areamatrix[i] = NULL;
        }
        free(areamatrix);
      }
      areamatrix=NULL;
      return err;

    case 0: default:
      SERR("optimization might not be used here\n");
  }
  
}

int init_BRmixed(cov_model *cov, storage *s) {
  cov_model *key = cov->key;
  location_type *loc = Loc(cov), *keyloc;
  int i, j, k, l, d, err, dim,
      shiftedloclen, keytotal,
      xbound, ybound,
      optim = ((int*) cov->p[BR_OPTIM])[0];
  BR_storage *sBR;
  double area = 1.0, *optimarea, *lowerbounds,
         step = cov->p[BR_MESHSIZE][0],
         xmin = cov->p[BR_LB_OPTIMAREA][0];
  bool grid = loc->grid;

  assert(isPointShape(cov));
  
  if (cov->role == ROLE_BROWNRESNICK) {
      
    if (key != NULL) {
      dim = cov->tsdim;      
      if ((err = alloc_cov(cov, dim, 1)) != NOERROR) return err;
      pgs_storage *pgs = cov->Spgs;
      for (d=0; d<dim; d++)
	pgs->supportmin[d]=pgs->supportmax[d]=pgs->supportcentre[d]=RF_NAN;
      GetDiameter(loc, pgs->supportmin, pgs->supportmax, pgs->supportcentre);
      sBR = cov->SBR;
      for (d=0; d<dim; d++) {
        pgs->supportmax[d] += sBR->radius + step;
        pgs->supportmin[d] -= sBR->radius + step;
      }
      for (d=0; d<cov->mpp.moments; d++) {
        cov->mpp.M[d] = cov->mpp.Mplus[d] = 1.0;
      }
      pgs->log_density = 0.0;

      keyloc = Loc(key);
      keytotal = keyloc->totalpoints;
      ((int*) cov->key->p[LOG_GAUSS])[0] = false;
      key->simu.active = true;
      key->simu.expected_number_simu = cov->simu.expected_number_simu;     
      assert(cov->mpp.moments >= 1);
      if ((err = INIT(key, 1, s)) != NOERROR) return err;  
      
      cov->loggiven = true;
      
      assert(key->nr == GAUSSPROC);
      assert(key->mpp.moments >= 1);
      cov->mpp.M[0] = cov->mpp.Mplus[0] = 1.0;
      cov->mpp.M[1] = cov->mpp.Mplus[1] = 1.0;
            
      sBR = cov->SBR;
      sBR->trendlen = 1;
      if ((cov->SBR->trend = (double**) MALLOC(sizeof(double*)))==NULL) {
         err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }    
      if ((cov->SBR->trend[0] = 
             (double*) MALLOC(keytotal*sizeof(double)) )==NULL) {
         err = ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      if ((err = loc_set(keyloc->xgr[0], NULL, NULL, dim, dim, 3, 0,
	                 false, true, keyloc->distances, 
			 &Loc(cov->SBR->vario))) > NOERROR) return err;     
      if (sBR->vario->sub[0] != NULL)
	SetLoc2NewLoc(sBR->vario->sub[0], Loc(sBR->vario));
      Variogram(NULL, sBR->vario, sBR->trend[0]);
    
      shiftedloclen = grid ? 3 : loc->totalpoints;
      if ( (cov->SBR->shiftedloc = (double*) 
              MALLOC(dim*shiftedloclen*sizeof(double))) == NULL) {
         err = ERRORMEMORYALLOCATION; goto ErrorHandling; 
      }
      for (d=0; d<dim; d++) {
        if (grid) {
	  sBR->shiftedloc[3*d+XSTART]  = loc->xgr[d][XSTART];
	  sBR->shiftedloc[3*d+XLENGTH] = loc->xgr[d][XLENGTH];
	  sBR->shiftedloc[3*d+XSTEP]   = loc->xgr[d][XSTEP];
        } else {
          for (i=d; i<shiftedloclen*dim; i+=dim) {
	    sBR->shiftedloc[i] = loc->x[i];  
	  }
        }
      }
      
      if ((err = FieldReturn(cov)) != NOERROR)  // must be later than INIT !
        return err;
      
      if (optim && (err = BRoptim(cov, s)) != NOERROR) return err;
      
      xbound = (int) floor(cov->ncol[BR_OPTIMAREA]/2.0);
      ybound = (int) floor(cov->nrow[BR_OPTIMAREA]/2.0);
      optimarea = cov->p[BR_OPTIMAREA];
      
      if ((cov->SBR->lowerbounds=(double*) MALLOC(keytotal*sizeof(double)))
          == NULL) {
        err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      lowerbounds = sBR->lowerbounds;
      for (i=0; i<keytotal; i++) lowerbounds[i] = RF_INF;
      
      //printf("lowerbounds:\n");
      for (l=0, i=-ybound; i<=ybound; i++) {
	k = sBR->zeropos + i*keyloc->length[0] - xbound;
        for (j=-xbound; j<=xbound; j++, l++, k++) {
	  if (optimarea[l]>1e-5) {
	    lowerbounds[k] = xmin - log(optimarea[l]);
	    //printf("%f ", lowerbounds[k]); 
	  }
        }
        //printf("\n");
      }
    
      for (d=0; d<dim; d++) area *= pgs->supportmax[d] - pgs->supportmin[d];
      cov->mpp.maxheight = log(cov->p[BR_LAMBDA][0]);
      cov->mpp.log_zhou_c = log(area);
   }
    
    return NOERROR;
  }
  
  else ILLEGAL_ROLE;
  
  ErrorHandling:
    BR_DELETE(&(cov->SBR));
    return err;
  
}

void do_BRmixed(cov_model *cov, storage *s) {
  assert(cov->key!=NULL);
  
  BR_storage *sBR = cov->SBR;
  cov_model  *key = cov->key;
  res_type   *res = cov->rf,
           *lgres = key->rf;
  assert(cov->origrf);
  location_type *loc = Loc(cov),
                *gridloc = Loc(key);
  assert(gridloc->grid);

  pgs_storage *pgs = cov->Spgs;
  assert(pgs != NULL);
  
  int i, d, maxind, index, multiindex[MAXMPPDIM],
      totalpoints = loc->totalpoints,
      lgtotalpoints = gridloc->totalpoints,
      *lglength = Loc(key)->length,
      zeropos = sBR->zeropos,
      dim = cov->tsdim;
  bool hat_simu_active = true,
      grid = loc->grid,
      withinradius;
  double nullval, maxval,
      summand, u, x, 
      xmin = cov->p[BR_LB_OPTIMAREA][0],
      *shiftedloc = sBR->shiftedloc,
      *lowerbounds = sBR->lowerbounds,
      *trend = sBR->trend[0],
      lambda = cov->p[BR_LAMBDA][0],
      radius = sBR->radius,
      step = cov->p[BR_MESHSIZE][0];

  sBR->hatnumber=0;
  
  for (d=0; d<dim; d++) {
    u = pgs->supportmin[d] + 
      UNIFORM_RANDOM*(pgs->supportmax[d]-pgs->supportmin[d]);
    if (grid) {
      shiftedloc[3*d+XSTART] = loc->xgr[d][XSTART] - u;
    } else {
      for (i=0; i<totalpoints; i++)
        shiftedloc[i*dim+d] = loc->x[i*dim+d] - u;
    }
  }
  
  while(hat_simu_active) {
    summand = xmin - log(UNIFORM_RANDOM);
    DO(key, s);
    maxval = RF_NEGINF;
    maxind = 0;
    nullval = lgres[zeropos];
    for (i=0; i<lgtotalpoints; i++) {
      lgres[i] += summand - nullval - trend[i];
      if (lgres[i] > maxval) {
	maxval = lgres[i];
	maxind = i;
      }
    }
    if (maxval > lowerbounds[maxind]) {
      hat_simu_active = false;
    } else {
      sBR->hatnumber++;
    }
  }
  
  //shifting maximum to origin is not necessary because of stationarity
  //(conditional on T=t, the "correct" shape function is shifted which yields
  //the same stationary max-stable process; OK because there is a finite number
  //of values for t!
  for (i=0; i<lgtotalpoints; i++) lgres[i] += log(lambda) - maxval;
  
  for (i=0; i<totalpoints; i++) {
    res[i] = 0.0;
    withinradius = true;
    if (grid) indextrafo(i, loc->length, dim, multiindex);
    index = 0;
    for (d=dim-1; d>=0; d--) {
      x = grid ? shiftedloc[3*d+XSTART] + multiindex[d]*shiftedloc[3*d+XSTEP] :
                 shiftedloc[i*dim+d];
      if (fabs(step*floor(x/step)) > radius) withinradius = false;
      index = index*lglength[d] + (int) (floor(x/step) + floor(lglength[d]/2.0));
    }
    res[i] = withinradius ? lgres[index] : RF_NEGINF;
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
  
  range->min[BR_LB_OPTIMAREA] = -RF_INF;
  range->max[BR_LB_OPTIMAREA] = RF_INF;
  range->pmin[BR_LB_OPTIMAREA] = -6;
  range->pmax[BR_LB_OPTIMAREA] = 0;
  range->openmin[BR_LB_OPTIMAREA] = true;
  range->openmax[BR_LB_OPTIMAREA] = true;
  
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
  range->pmax[BR_VARIOBOUND] = 6;
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
  double *newx=NULL,
         centreloc[MAXMPPDIM], 
         minloc[MAXMPPDIM], 
	 maxloc[MAXMPPDIM];
  
  if (newmodel != NULL) SERR("unexpected call of structBRuser"); /// ?????
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
  if (cov->stor == NULL) cov->stor = (storage *) MALLOC(sizeof(storage));    
  STORAGE_NULL(cov->stor);
 
  GetDiameter(loc, minloc, maxloc, centreloc);
  newxlen = loc->lx;
  if ((newx = (double*) MALLOC(dim * newxlen * sizeof(double))) == NULL) {
     GERR("Memory allocation failed.\n"); 
  } 
  if (grid) {
    for (d=0; d<dim; d++) {
      newx[3 * d + XSTART] = loc->xgr[d][XSTART] - centreloc[d] 
           + (((int) loc->xgr[d][XLENGTH]) % 2 == 0)*loc->xgr[d][XSTEP]/2;
      newx[3 * d + XSTEP] = loc->xgr[d][XSTEP];
      newx[3 * d + XLENGTH] = loc->xgr[d][XLENGTH];
   }
  } else {
    for (i=0; i<loc->lx; i++)
      for(d=0; d<dim; d++) newx[i*dim+d] = loc->x[i*dim+d] - centreloc[d];
  }

  if ((err = loc_set(newx, NULL, dim, dim, newxlen, false, loc->grid,
                     loc->distances,  &(cov->ownloc))) > NOERROR)
    return err;
  SetLoc2NewLoc(sub, cov->ownloc);
      
  if (newx!=NULL) free(newx);
  newx = NULL;

  if ((err=covcpy(&(cov->key), sub)) > NOERROR) return err;
  
  if (cov->sub[MPP_TCF] != NULL) {
     if ((err = STRUCT(sub, &(cov->key))) > NOERROR) return err;
     cov->key->calling = cov;
  } 

  addModel(&(cov->key), model_intern);
  
  kdefault(cov->key, GEV_XI, cov->p[GEV_XI][0]);
  kdefault(cov->key, GEV_MU, cov->p[GEV_MU][0]);
  kdefault(cov->key, GEV_S, cov->p[GEV_S][0]);

  if (cov->nr == BRMIXED_USER) {
    kdefault(cov->key, BR_MESHSIZE, cov->p[BR_MESHSIZE][0]);
    kdefault(cov->key, BR_LB_OPTIMAREA, cov->p[BR_LB_OPTIMAREA][0]);    
    kdefault(cov->key, BR_VERTNUMBER, ((int*) cov->p[BR_VERTNUMBER])[0]);
    kdefault(cov->key, BR_OPTIM, ((int*) cov->p[BR_OPTIM])[0]);
    kdefault(cov->key, BR_OPTIMTOL, cov->p[BR_OPTIMTOL][0]);
    kdefault(cov->key, BR_VARIOBOUND, cov->p[BR_VARIOBOUND][0]);  
    kdefault(cov->key, BR_OPTIMMAX, ((int*) cov->p[BR_OPTIMMAX])[0]);
    kdefault(cov->key, BR_LAMBDA, cov->p[BR_LAMBDA][0]);
    if (cov->p[BR_OPTIMAREA] != NULL) {
      cov->key->nrow[BR_OPTIMAREA] = cov->nrow[BR_OPTIMAREA]; 
      cov->key->ncol[BR_OPTIMAREA] = cov->ncol[BR_OPTIMAREA];
      if ( cov->nrow[BR_OPTIMAREA] % 2 == 0 ||
	   cov->ncol[BR_OPTIMAREA] % 2 == 0 ) {
	SERR("number of rows and columns of areamat need to be odd")
      }      
      cov->key->p[BR_OPTIMAREA] = (double*) 
      MALLOC(sizeof(double)*cov->nrow[BR_OPTIMAREA]*cov->ncol[BR_OPTIMAREA]);
      MEMCOPY(cov->key->p[BR_OPTIMAREA], cov->p[BR_OPTIMAREA], 
	      sizeof(double)*cov->nrow[BR_OPTIMAREA]*cov->ncol[BR_OPTIMAREA]);
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
  if (err > NOERROR) return err;
  assert(cov->key->calling == cov);
  
  return NOERROR;
  
  ErrorHandling:
     if (newx != NULL) free(newx);
     newx = NULL;
     return err;
  
}

int structBRintern(cov_model *cov, cov_model **newmodel) {
  cov_model *sub = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE: MPP_TCF];
  location_type *loc = Loc(cov);
  int i, d, j, err,  *length,
      dim = sub->tsdim, 
      totaldim = loc->totalpoints*dim,
      zeropos = 0, newxlen;
  bool grid = loc->grid;
  double *newx = NULL, 
        covval, norm, step,
	mindist, dist,
	radius[MAXMPPDIM];
  BR_storage *sBR = NULL;
  
  if (newmodel != NULL) SERR("unexpected call of structBR"); /// ?????
  if (cov->role != ROLE_BROWNRESNICK) BUG;

  assert(isPointShape(cov));
    
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if (cov->stor == NULL) cov->stor = (storage *) MALLOC(sizeof(storage));    
  STORAGE_NULL(cov->stor);
  
  if (cov->SBR != NULL) BR_DELETE(&(cov->SBR));
  if ((cov->SBR = (BR_storage*) MALLOC(sizeof(BR_storage)))==NULL) {
     err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  sBR = cov->SBR;
  BR_NULL(sBR);
  
  if ((err=covcpy(&(cov->key), sub)) > NOERROR) return err;
  
  if (cov->sub[MPP_TCF] != NULL) {
    if ((err = STRUCT(sub, &(cov->key))) > NOERROR) 
	return err;
    cov->key->calling = cov;
  }
  
  CHECK(cov->key, dim, dim, NegDefType, cov->domown, SYMMETRIC, 1, ROLE_COV);
  
  newmodel_covcpy(&(cov->SBR->vario), VARIOGRAM_CALL, cov->key);
  if ((err = alloc_cov(cov->SBR->vario, dim, 1)) != NOERROR) return err;
	
  addModel(&(cov->key), GAUSSPROC);
  cov->key->calling = cov;
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
      if ((newx = (double*) MALLOC(dim*newxlen*sizeof(double))) == NULL) {
        err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      for(d=0; d<dim; d++) {
        for (i=0; i<loc->lx; i++) newx[i*dim+d] = loc->x[i*dim+d];
	if (zeropos == loc->lx) newx[loc->lx*dim+d] = 0.0;
      }  
    }
  } else if (cov->nr == BRMIXED_INTERN) {
    step = cov->p[BR_MESHSIZE][0];
    mindist = RF_INF;
    if (grid) {
      for (d=0; d<dim; d++)  
        if (loc->xgr[d][XSTEP] < mindist)
	  mindist = loc->xgr[d][XSTEP];
    } else {
      for (i=0; i<totaldim; i+=dim)
        for (j=i+dim; j<totaldim; )
	  for (d=0; d<dim; d++, j++) {
	    dist = fabs(loc->x[i+d] - loc->x[j]);
	    if (dist > 1e-15 && dist < mindist) mindist = dist;
	  }
    } 
    if (mindist < step) {
      cov->p[BR_MESHSIZE][0] = mindist;
      step = cov->p[BR_MESHSIZE][0];
      PRINTF("Warning! meshsize is larger than resolution of simulation! meshsize is automatically decreased to %f.\n", step);
    }
    for (d=0; d<MAXMPPDIM; d++) radius[d] = 0.0;
    if (cov->p[BR_OPTIMAREA] != NULL) {
      radius[0] = step*((int) floor(fmax(cov->nrow[BR_OPTIMAREA],
                                         cov->ncol[BR_OPTIMAREA])/2.0) - 1);
    }
    while (true) {
      radius[0] += step; 
      if ((err = loc_set(radius, NULL, NULL, dim, dim, 1, 0, false, false,
	                 false, &Loc(sBR->vario))) > NOERROR) return err;
      if (sBR->vario->sub[0] != NULL) 
	 SetLoc2NewLoc(sBR->vario->sub[0], Loc(sBR->vario));
      Variogram(NULL, sBR->vario, &covval);
      if (covval >= cov->p[BR_VARIOBOUND][0]) break;
    }
    sBR->radius = radius[0];
    //printf("radius: %f\n", radius[0]);
    newxlen = 3;
    grid = true;
    
    if ((newx = (double*) MALLOC(newxlen*dim*sizeof(double))) == NULL) {
      err = ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    for (d=0; d<dim; d++) {
      newx[3*d+XSTART] = -sBR->radius;
      newx[3*d+XLENGTH] = 2*((int) round(sBR->radius/step))+1;
      newx[3*d+XSTEP] = step;
    }
  }
  
  if (newx == NULL) {
    err = loc_set(grid ? *loc->xgr : loc->x, NULL, dim, dim,
	          grid ? 3 : loc->totalpoints, false, grid,
	          loc->distances, &(cov->key->ownloc));
  } else {
    err = loc_set(newx, NULL, dim, dim, newxlen, false, grid, 
                  false, &(cov->key->ownloc));
    free(newx);
    newx=NULL;
  }
  if (err > NOERROR) return(err);
  
  if (grid) {
    length = Loc(cov->key)->length;
    for (d=dim; d>0; d--)
      zeropos = zeropos*length[d-1] + ceil(length[d-1]/2.0) - 1;
  }
  sBR->zeropos = zeropos;

  if ((err = CHECK(cov->key, dim, dim, ProcessType, cov->domown, 
                   cov->isoown, 1, ROLE_GAUSS)) == NOERROR) {
    if ((err = STRUCT(cov->key, NULL)) <= NOERROR) {
      err = CHECK(cov->key, dim, dim, ProcessType, cov->domown,
		  cov->isoown, 1, ROLE_GAUSS); 
    }
  }
  if (err > NOERROR) return err;

  assert(cov->key->calling == cov);
  
  return NOERROR;

  ErrorHandling:
     if (newx != NULL) free(newx);
     newx = NULL;
     BR_DELETE(&(cov->SBR));    
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
      
    if (newmodel != NULL) SERR("unexpected call of structBR "); /// ?????
    if (next->full_derivs < 0) SERR("given submodel does not make sense");
 
    while (isDollar(next)) {
      addModel(&(cov->key), DOLLAR);
      //     cov->key->calling = 
      //      calling->sub[ATOMSHAPE]->calling = calling;
      //      calling = (*shape)->sub[0];
      // shape = calling->sub + 0;
      // print("here2\n");
       
      //  printf("model %s %d %d \n", CovList[(*shape)->nr].name,
      //	     next->p[DSCALE] != NULL, next->p[DVAR] != NULL);

      newscale = 1.0;
      if (next->p[DSCALE] != NULL) newscale *= next->p[DSCALE][0];
      if (next->p[DVAR] != NULL)  newscale /= sqrt(next->p[DVAR][0]);
      if (factor != 1.0) {
	newscale *= factor;
	factor = 1.0;
      }
      //      print("here3\n");
      return ERRORNOTPROGRAMMED;
      kdefault(calling, DSCALE, newscale);
      next = next->sub[0]; // ok
      //     print("here4\n");
    }

    if (cov->sub[MPP_TCF] != NULL) {
      // print("here4\n");
      return ERRORNOTPROGRAMMED;
    } 

    // PMI(atom->sub[ATOMSHAPE]); PMI(*shape); assert(false);
    // print("here5\n");
    
    if (next->nr == BROWNIAN && next->p[BROWN_ALPHA][0] == 2.0) {
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
      role = isNegDef(next) ? ROLE_COV 
	     : ROLE_UNDEFINED;
	
      ASSERT_ROLE_DEFINED(next);  
      
      if (((err = covcpy(&(cov->key), next)) != NOERROR) 
         || ((err = CHECK(cov->key, dim, dim, NegDefType, XONLY,
              SYMMETRIC, 1, role)) != NOERROR)) 
	 return err;
      GetDiameter(loc, minloc, maxloc, centreloc);
      for (d=0; d<MAXMPPDIM; d++) maxdist[d] = 0.5*(maxloc[d] - minloc[d]);
   
      cov_model *K = NULL;
      newmodel_covcpy(&K, VARIOGRAM_CALL, cov->key, maxdist, NULL, NULL, 
		      dim, dim, 1, 0, false, false, false);
      if ((err = alloc_cov(K, dim, 1)) != NOERROR) return err;
      if (K->sub[0] != NULL) SetLoc2NewLoc(K->sub[0], Loc(K));
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

      addModel(&(cov->key), meth);
      cov->key->calling = cov;
      cov_model *key = cov->key;
      key->prevloc = loc;
      
      kdefault(key, GEV_XI, cov->p[GEV_XI][0]);
      kdefault(key, GEV_MU, cov->p[GEV_MU][0]);
      kdefault(key, GEV_S, cov->p[GEV_S][0]);
       
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

int initBrownResnick (cov_model *cov, storage *S) {

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

void doBrownResnick(cov_model *cov, storage *s) {

  assert(!cov->origrf);
  assert(cov->key != NULL);
  cov_model *key = cov->key;
 
  DO(key, s); // nicht gatternr

}


int initBRuser (cov_model *cov, storage *S) {
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
      sub->simu.expected_number_simu = MIN(ens, MAXINT);

      if ((err = INIT(sub, 1, S)) != NOERROR) return err;
      FieldReturn(cov); 
    }
    
    return NOERROR;
  }

  else ILLEGAL_ROLE;
   
}
