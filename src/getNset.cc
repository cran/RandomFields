
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 (library for simulation of random fields)
 Copyright (C) 2001 -- 2017 Martin Schlather, 

This program is free software; you can redist ribute it and/or
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

#include <R.h>
#include <Rdefines.h>
#include <stdio.h>  
 #include <string.h>
#include <R_ext/Linpack.h>

#include "questions.h"
#include "primitive.h"
#include "Coordinate_systems.h"
#include "Processes.h"
#include "xport_import.h"
#include "operator.h"
#include "PoissonPolygon.h"
#include "families.h"



void listpt(listoftype **To, listoftype *p, int len, SEXPTYPE Rtype,
	    bool force_allocating) {
  if (*To == NULL || force_allocating) {
    *To = (listoftype *) MALLOC(sizeof(listoftype));
  }
  listoftype *q = *To;
  q->lpx = p->lpx;
  q->ncol = p->ncol;
  q->nrow = p->nrow;
  q->deletelist = false;
  q->len = len;
  q->Rtype = (int) Rtype;
}

void listcpy(listoftype **To, listoftype *p, bool force_allocating) {
  // force_allocating in case of "c ovcpy"
  int size,
    len = p->len,
    sizeint = len * sizeof(int);
  if (p->Rtype == LISTOF + REALSXP) {
    size = sizeof(double);
  } else BUG;
   
  if (*To == NULL || force_allocating) *To = LIST_CREATE(len, p->Rtype);
  listoftype *q = *To;

  for (int j=0; j<len; j++) {
    int n = size * p->nrow[j] * p->ncol[j];
    if (q->lpx[j] == NULL) q->lpx[j] = (double*) MALLOC(n);
    MEMCOPY(q->lpx[j], p->lpx[j], n);	    
  }
  MEMCOPY(q->nrow, p->nrow, sizeint);
  MEMCOPY(q->ncol, p->ncol, sizeint);
}



// Achtung! removeonly muss zwingend mit dem originalen SUB + i aufgerufen
// werden. Der **Cov ueberschrieben wird. Ansonsten muessen die
// Kommentarzeilen geloescht werden !
void removeOnly(model **Cov) {
  model *cov = *Cov,
    *next = cov->sub[0];
 
  if (cov->calling == NULL) next->calling = NULL;  
  else {
    model *calling = cov->calling;
//    for (i=0; i<MAXSUB; i++) if (calling->sub[i] == cov) break;
//    assert (i<MAXSUB);
//    calling->sub[i] = cov->sub[0];
    SET_CALLING(next, calling);
  }
  *Cov = next;
  COV_DELETE_WITHOUTSUB(&cov, next);
}



void addModel(model **pcov, int covnr, model *calling, bool nullOK) {
  model *cov;
  int i;
  cov = (model*) MALLOC(sizeof(model));
  COV_NULL(cov, calling == NULL ? NULL : calling->base);
  //
  assert(cov->calling == NULL);
  set_nr(OWN, covnr);
  set_last(OWN, 0, 0);
 
  if (*pcov != NULL) {
    cov->nsub = 1;
    cov->sub[0] = *pcov;
    calling = calling == NULL ? (*pcov)->calling : calling;
    SET_CALLING(cov, calling);
    (*pcov)->calling = cov; // NOT SET_CALLING, sice pcov->calling already set
    for (i=0; i<=Forbidden; i++) {	
      // cov->user[i] = cov->sub[0]->user[i];
      cov->pref[i] = cov->sub[0]->pref[i];
    }
  } else {
    SET_CALLING(cov, calling);
  }
 
  if (calling == NULL && !nullOK) {
    PRINTF("Missing link for model '%.50s'. Inform author.\n", NICK(cov));
    BUG;
  }
  
 *pcov = cov;
}

void addModel(model **pcov, int covnr) {
  assert(*pcov != NULL && (*pcov)->calling != NULL);
  addModel(pcov, covnr, NULL, false);
}

void addModel(model *pcov, int subnr, int covnr) {
  assert(pcov != NULL);
  bool newsub = pcov->sub[subnr] == NULL;
  addModel(pcov->sub + subnr, covnr, pcov, false);
  pcov->nsub += (int) newsub;
}

void addModelKappa(model *pcov, int subnr, int covnr) {
  assert(pcov != NULL);
  addModel(pcov->kappasub + subnr, covnr, pcov, false);
}

void addModel(model **pcov, int covnr, model *calling) {
  addModel(pcov, covnr, calling, false);
}


int addUnifModel(model *cov, double radius, model **newmodel) {
  addModel(newmodel, UNIF, cov);
  kdefault(*newmodel, UNIF_MIN, -radius);
  kdefault(*newmodel, UNIF_MAX, radius);
  RETURN_NOERROR;
}


int setgrid(coord_type xgr, double *x, int spatialdim) {
  int lx = 3;
  
  int d;
  unsigned long
    totalBytes = sizeof(double) * lx * spatialdim; // nothing necessary for time
 
  if (xgr[0] == NULL && (xgr[0] =(double*) MALLOC(totalBytes))==NULL)
    return ERRORMEMORYALLOCATION;

  MEMCOPY(xgr[0], x, totalBytes);
  
  // folgende Zeile nur beim ersten Mal zu setzen, aber egal
  for (d=1; d<spatialdim; d++) {
    xgr[d] = &(xgr[0][d * lx]); 
    if (xgr[d][XLENGTH] != (int) xgr[d][XLENGTH]) 
      FAILED2("grid length must be integer valued. Got %10e in dimension %d.",
 	    xgr[d][XLENGTH], d);
     if (xgr[d][XLENGTH] < 1.0) 
      FAILED2("grid length must be positive. Got %10e in dimension %d.",
 	    xgr[d][XLENGTH], d);
   }
  
  /*
  if (GLOBAL.internal.examples_reduced) {    
    for (d=0; d<spatialdim; d++) {
      if (xgr[d][XLENGTH] > GLOBAL.internal.examples_reduced) {
	warning("the size of the example has been reduced");
	xgr[d][XLENGTH] = GLOBAL.internal.examples_reduced;
      }
    }
  }
  */
  
  return NOERROR;
}

int partial_loc_set(location_type *loc, double *x, double *y,
		    long lx, long ly, bool dist, int xdimOZ, double *T,
		    bool grid, bool cpy){
  int d, Err;
  unsigned long totalBytes;

  if (lx >= MAXINT ||  ly >= MAXINT) return XERRORTOOMANYLOC;

  if (((loc->x != NULL) && ((loc->y == NULL) xor (ly==0))) ||
      ((loc->xgr[0] != NULL) && ((loc->ygr[0] == NULL) xor (ly==0))))     
    FAILED("domain structure of the first and second call do not match");
 

  assert(x != NULL);

  /*
  if (GLOBAL.internal.examples_reduced) {
    if (lx > GLOBAL.internal.examples_reduced) {
      //  lx = GLOBAL.internal.examples_reduced;  
      warning("The size of the example has been reduced.");
    }

    if (ly > GLOBAL.internal.examples_reduced) {
      ly = GLOBAL.internal.examples_reduced;  
      warning("The size of y-coordinates in the example has been reduced.");
    }
  }
  */

  loc->xdimOZ = xdimOZ; // ohne Zeit !!
  loc->lx = lx;
  loc->ly = ly;
  if (ly  > 0) { 
    //crash();
    if (dist) FAILED("distances are not allowed if y is given");
  }

  loc->grid = grid;
  loc->distances = dist;
  if (loc->delete_y && loc->y != loc->x) FREE(loc->y); 
  if (loc->delete_x) FREE(loc->x);
  loc->delete_x = loc->delete_y = cpy;
  if (lx == 0) return NOERROR;
  
  if (grid) {
    loc->delete_x = true;
    assert(lx == 3);
    if ((Err = setgrid(loc->xgr, x, loc->spatialdim)) != NOERROR) return Err;
    if (ly>0) {
      if (x == y) {
	for(d=0; d<loc->spatialdim; d++) loc->ygr[d] = loc->xgr[d];
	loc->delete_y = false;
      } else {
	assert(ly == 3);
	if ((Err = setgrid(loc->ygr, y, loc->spatialdim)) != NOERROR) return Err;
      }
    }
    
    double len = 1.0;
    for (d=0; d<loc->spatialdim; d++) {
      len *= loc->xgr[d][XLENGTH];
    }
    // if (len != (int) len) crash();
    assert(len == (int) len);
    if (len < MAXINT) loc->totalpoints = loc->spatialtotalpoints = (int) len;
    else return XERRORTOOMANYLOC;
  } 

  else if (dist) {
    // lx : anzahl der Puntke involviert, d.h. x muss die Laenge lx ( lx -1) / 2
    // haben.
    if (lx > 0) {
      if (cpy) {
	totalBytes =  sizeof(double) * lx * (lx - 1) / 2 * xdimOZ;
	if ((loc->x=(double*) MALLOC(totalBytes))==NULL)
	  return ERRORMEMORYALLOCATION;
	MEMCOPY(loc->x, x, totalBytes);
      } else {
	loc->x = x;
      }
    }
    loc->totalpoints = loc->spatialtotalpoints = lx;
    //     (int) (1e-9 + 0.5 * (1.0 + SQRT(1.0 + 8.0 * lx)));
    //   if (0.5 * (loc->totalpoints * (loc->totalpoints - 1.0)) != lx) {   
    //printf("tot=%d %d %d\n", loc->totalpoints, (int) ( 0.5 * (loc->totalpoints * (loc->totalpoints - 1.0))), (int) lx); assert(false);
    //   SERR("distances do not have expected size");
    //}
  }
  
  else { // not grid, not distances
    if (cpy) {
      totalBytes =  sizeof(double) * lx * loc->xdimOZ;
      assert(x != NULL);
      assert(loc != NULL && loc->x == NULL); 
      if ((loc->x=(double*) MALLOC(totalBytes)) == NULL)
	return ERRORMEMORYALLOCATION; 
      assert(loc->x != NULL);

      // 
      MEMCOPY(loc->x, x, totalBytes);
      // for (int k=0; k<lx *loc->xdimOZ; k++) loc->x[k] = x[k];


      if (loc->ly>0) {
	if (x == y) {
	  loc->y = loc->x;
	  loc->delete_y = false;
	} else {
	  totalBytes =  sizeof(double) * ly * loc->xdimOZ;
	  if ((loc->y=(double*) MALLOC(totalBytes))==NULL)
	    return ERRORMEMORYALLOCATION; 
	  MEMCOPY(loc->y, y, totalBytes);	
	}
      }
    } else {
      loc->x = x;
      loc->y = y;
    }
    
    loc->totalpoints = loc->spatialtotalpoints = lx;
    // just to say that considering these values does not make sense
  }
  
  if ((loc->Time) xor (T!=NULL)) {    
    //model *cov;  crash(cov);
    //  AERR(1);
    FAILED("partial_loc: time mismatch");
  }
  if (loc->Time) {
    MEMCOPY(loc->T, T, sizeof(double) * 3);
    if (grid) {
      loc->xgr[loc->spatialdim] = loc->T;
      if (ly>0) loc->ygr[loc->spatialdim] = loc->T;
    }
    
    if (loc->T[XLENGTH] <= 0) {
      //crash();
      FAILED1("The number of temporal points is not positive. Check the triple definition of 'T' in the man pages of '%.50s'.", DefList[SIMULATE].nick);
    }

    /* else if (GLOBAL.internal.examples_reduced &&
	       loc->T[XLENGTH] > GLOBAL.internal.examples_reduced) {
      loc->T[XLENGTH] = GLOBAL.internal.examples_reduced;  
      warning("The length of the tiem component in the example has been reduced.");
      }
    */

    if ((double) loc->totalpoints * loc->T[XLENGTH] >= MAXINT) 
      FAILED("too many space-time locations");
    loc->totalpoints *= (int) loc->T[XLENGTH];
  }

  return NOERROR; 
}

int loc_set(double *x, double *y, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdimOZ,
	    long lx, long ly, bool Time, bool grid,
	    bool distances,
	    location_type **LLoc) {
  int Err;
  //unsigned long totalBytes;
  // preference lists, distinguished by grid==true/false and dimension
  // lists must end with Nothing!
 
  if (xdimOZ < spatialdim) {
    if (distances) {
      if (xdimOZ != 1) FAILED("reduced dimension is not one");
    } else FAILED3("dim (%d) of 'x' does not fit the spatial dim (%d); Time=%d",
		   xdimOZ,spatialdim, Time);
  } else if (xdimOZ > spatialdim)
    FAILED3("mismatch of dimensions (xdim=%d > space=%d; Time=%d)",
	   xdimOZ, spatialdim, Time);

  //int len = *LLoc == NULL ? 1 : (*LLoc)->len;
  if (*LLoc != NULL && (*LLoc)->lx > 0) {
    BUG; // OK
    LOC_SINGLE_DELETE(LLoc);
    *LLoc = (location_type*) MALLOC(sizeof(location_type));
  }

  location_type *loc = *LLoc;
  assert(loc->xgr != NULL && loc->xgr[0] == NULL);
  assert(loc->ygr != NULL && loc->ygr[0] == NULL);
 
  loc->timespacedim = spatialdim + (int) Time;
  loc->spatialdim = spatialdim;
  loc->Time = Time; 

  if (spatialdim<1) return ERRORDIM;

  assert(x != NULL);

  
  if ((Err = partial_loc_set(*LLoc, x, y, lx, ly, distances, xdimOZ,
			     Time ? T : NULL,
			     grid, true)) != NOERROR) XERR(Err);
 
  return NOERROR;
}




int loc_set(double *x, double *y, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdimOZ,
	    long lx, long ly,
	    bool Time, bool grid,  bool distances,
	    model *cov) {
  int Err,
    store = GLOBAL.general.set;
  location_type **oldloc = cov->ownloc;
  GLOBAL.general.set = 0;

  cov->ownloc = LOCLIST_CREATE(1, xdimOZ + (int) Time); // locown
  assert(cov->ownloc != NULL);
  assert(PLoc(cov) != cov->prevloc);

  Err = loc_set(x, y, T, spatialdim, xdimOZ, lx, ly, Time, grid,
		distances, cov->ownloc);    
  // Errorhandling:
  GLOBAL.general.set = store;
  LOC_DELETE(&oldloc);
  assert(grid xor !Loc(cov)->grid);
  return Err;
}



int loc_set(double *x, double *T, 
	    int spatialdim, // spatial dim only !
	    int xdimOZ, // original ! 
	    long lx, bool Time, bool grid,
	    bool distances,
	    location_type **Loc) {
  return loc_set(x, NULL, T, spatialdim, xdimOZ, lx, 0, Time, grid,
		 distances, Loc);    
} 




int loc_set(double *x, double *T, 
	    int spatialdim, /* spatial dim only ! */
	    int xdimOZ, /* original ! */
	    long lx, bool Time, bool grid,
	    bool distances,
	    model *cov) {
  return loc_set(x, NULL, T, spatialdim, xdimOZ, lx, 0, Time, grid,
		 distances, cov);    
} 


location_type ** loc_set(SEXP xlist, bool distances_ok){
  bool
    listoflists = (TYPEOF(xlist) == VECSXP &&
		   TYPEOF(VECTOR_ELT(xlist, 0)) == VECSXP);
  int lx, Err,
    spatialdim = NA_INTEGER, 
    xdimOZ = UNSET,
    sets = !listoflists ? 1 : length(xlist);
  bool
    Time = false,
    distances = false;
  location_type **loc;

  for (int i=0; i<sets; i++) {
    SEXP
      set = listoflists ? VECTOR_ELT(xlist, i) : xlist,
      xx = VECTOR_ELT(set, XLIST_X),
      yy = VECTOR_ELT(set, XLIST_Y),
      TT = VECTOR_ELT(set, XLIST_T);   
    bool
      ggrid = LOGICAL(VECTOR_ELT(set, XLIST_GRID))[0]; 
    int
      xxdimOZ = ggrid ? ncols(xx) : nrows(xx),
      llx = ggrid ? 3 : ncols(xx),
      lly = length(yy) == 0 ? 0 : ggrid ? 3 : ncols(yy);      
    
    if (i==0) {
      xdimOZ = xxdimOZ;
      spatialdim = INTEGER(VECTOR_ELT(set, XLIST_SPATIALDIM))[0];
      Time =  LOGICAL(VECTOR_ELT(set, XLIST_TIME))[0];
      distances = LOGICAL(VECTOR_ELT(set, XLIST_DIST))[0];  
      loc = LOCLIST_CREATE(sets, xdimOZ + (int) Time);   
    } else {
      if (xdimOZ != xxdimOZ ||
	  spatialdim != INTEGER(VECTOR_ELT(set, XLIST_SPATIALDIM))[0] ||
	  Time != LOGICAL(VECTOR_ELT(set, XLIST_TIME))[0] ||
	  distances != LOGICAL(VECTOR_ELT(set, XLIST_DIST))[0]
	  ) BUG;
    }
  
    if (distances) {
      if (distances_ok){
	lx = (int) (1e-9 + 0.5 * (1 + SQRT(1. + 8 * llx)));
	if (llx != lx * (lx - 1) / 2) 
	  RFERROR("distance length not of form 'n * (n - 1) / 2'");	
      } else {
	NotProgrammedYet("currently: distances in simulations");
	distances = false;
      }
    } else {
      lx = llx;
      if (xxdimOZ != xdimOZ) 
	RFERROR("dimensions of the coordinates are not the same for the different sets");
    }
 
   if ((Err = loc_set(REAL(xx), REAL(yy), REAL(TT), spatialdim, xdimOZ, 
		       lx, lly, Time, ggrid, 
		       distances, loc + i)) != NOERROR) {
      LOC_DELETE(&loc); 
      XERR(Err);
    }
  }

  return loc;
}



int getmodelnr(char *name) {
  // NOMATCHING, -1, if no matching function is found
  // MULTIPLEMATCHING, -2, if multiple matching fctns are found,
  //  without one matching exactly
  // MATCHESINTERNAL, -3 if internal name is passed
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)
  int match;

  // to do : alphabetisch geordnet namen + geordnet nach Anfangsbuchstaben bring en Beschleunigung. Falls nicht mit "R" anfaengt, billige Suche 

  assert(currentNrCov != UNSET);
  
  if (!STRCMP(name, InternalName)) return MATCHESINTERNAL;
  if ((match = Match(name, CovNickNames, currentNrCov)) >= 0) return match;
  return Match(name, CovNames, currentNrCov);
}

void MultiDimRange(int set, model *cov, double *natscale) {
  int wave, redxdim, d, idx, 
    xdimprev = PREVTOTALXDIM,
    vdim = VDIM0,
    lastsystem = PREVLASTSYSTEM,
    store = GLOBAL.general.set,
    err = NOERROR;
  double y, yold, threshold, natsc, factor, Sign,
    newx, xsave, newy, 
    *x = NULL,
    *dummy =  NULL,
    rel_threshold = 0.05;
  bool islogcart[MAXSYSTEMS],
    domain = isXonly(OWN);
  //  covfct cf=NULL;
  // nonstat_covfct ncf=NULL;
  
  assert( HaveSameSystems(PREV, OWN));
  redxdim = OWNTOTALXDIM;
  GLOBAL.general.set = set;   

  //  if (redxdim > xdimprev) {
  if (redxdim != xdimprev) 
    GERR("dimension of x-coordinates too high to detect natural scaling.");
  
  if (cov->full_derivs < 0) { err=ERRORNOTDEFINED; goto ErrorHandling; }
 
  if ((dummy = (double*) MALLOC(sizeof(double) * vdim * vdim)) == NULL ||
      (x = (double *) MALLOC(sizeof(double) * xdimprev)) == NULL) 
    GERR("not enough memory when determining natural scaling.");
 
  if (cov->full_derivs < 0) { err=ERRORNOTDEFINED; goto ErrorHandling; }
  if (domain) {
    COV(ZERO(cov), cov, dummy);
  } else {
    double *zero = ZERO(cov);
    NONSTATCOV(zero, zero, cov, dummy);
  }
  threshold = rel_threshold * dummy[0];
  //
  int cumi;
  cumi = 0;
  for (int i=0; i<redxdim; i++) islogcart[i] = false;
  assert(lastsystem >=0);
  for (int s=0; s<=lastsystem; s++) {
    bool isLog = isLogCart(PREV, s);
    int xdim = PREVXDIM(s);
    for (int i=0; i<xdim; i++) islogcart[cumi++] = isLog;
  }
  assert(cumi == xdimprev);
  //

  // for (int s=0; s<=lastsystem; s++) {
  for (d=0; d<redxdim; d++) {
    wave  = 0;
    for (int i=0; i<xdimprev; i++) x[i] = (double) islogcart[i];     
    idx = (redxdim == xdimprev || d==0) ? d : xdimprev-1; 
    x[idx] = islogcart[idx] ? M_E : 1.0; // wrong compiler warning
    
    if (domain) COV(x, cov, dummy) else NONSTATCOV(ZERO(cov), x, cov, dummy); 
    yold = dummy[0];
    if (ISNAN(yold)) GERR("NA in model evaluation detected");
    if (yold > threshold) {
      factor = 2.0;
      Sign = 1.0;
    } else {
      factor = 0.5;
      Sign = - 1.0;
    } 
    
    double otherx;
    if (islogcart[idx]) x[idx] = POW(x[idx], factor); else x[idx] *= factor;
    if (domain) COV(x, cov, dummy) else NONSTATCOV(ZERO(cov), x, cov, dummy);
    y = dummy[0];
    
    while (Sign * (y - threshold) > 0) {  
      if (yold<y){ 
	if (wave++>10) { err=ERRORWAVING; goto ErrorHandling; }
      }
      yold = y;
      if (islogcart[idx]) x[idx] = POW(x[idx], factor); else x[idx] *= factor;
      if (x[idx]>1E30) { err=ERRORRESCALING; goto ErrorHandling; }
      if (domain) COV(x, cov, dummy) else NONSTATCOV(ZERO(cov), x, cov, dummy);
      y = dummy[0];
    }
    
    if (islogcart[idx]) otherx = POW(x[idx], 1 / factor);
    else otherx = x[idx] / factor;
    
    for (int i=0; i<3 /* good choice?? */ ;i++) {       
      if (y==yold) { err=ERRORWAVING; goto ErrorHandling; }
      double f = (threshold-y) / (y-yold);
      
      if (islogcart[idx]) newx = x[idx] * POW(x[idx] / otherx, f);
      else newx = x[idx] + (x[idx] - otherx) * f;
      xsave = x[idx];
      x[idx] = newx;
      if (domain) COV(x, cov, dummy) else NONSTATCOV(ZERO(cov), x, cov, dummy);
      newy = dummy[0];
      x[idx] = xsave;
      
      if (Sign * (newy - threshold) > 0) {
	otherx = newx;
	yold  = newy;
      } else {
	x[idx] = newx;
	y = newy;
      }
    }
    
    if (y==yold)  { err=ERRORWAVING; goto ErrorHandling; }
    natsc = 1.0 / ( x[idx] + 
		    (x[idx]-otherx)/(y-yold)*(threshold-y) );
    
    if (redxdim == xdimprev || d==0) {
      natscale[d] = natsc;
    } else {
      int beg, end;
      if (redxdim == 2) {
	if (d==0) {
	  beg = 0;
	  end = xdimprev-1;
	} else {
	  beg = end = xdimprev-1;
	}
      } else {
	beg = 0;
	end = xdimprev;
      }
      for (int i=beg; i<end; natscale[i++] = natsc);
    }
  } // reddim

 ErrorHandling:
  FREE(dummy);
  FREE(x);
  GLOBAL.general.set = store;
  if (err != NOERROR) XERR(err);
}

void MultiDimRange(int *model_nr, int *set, double *natscale) { 
  MultiDimRange(*set, KEY()[*model_nr], natscale); 
}


void GetNaturalScaling(model *cov, double *natscale)
{ // called also by R 

  // values of naturalscaling:
  //#define NATSCALE_EXACT 1   
  //#define NATSCALE_APPROX 2
  //#define NATSCALE_MLE 3 /* check fitvario when changing !! */
  
  defn *C = DefList + COVNR; // nicht gatternr
  *natscale = 0.0;

  if (C->maxsub!=0) XERR(ERRORFAILED); 
 
  if (!equalsIsotropic(DEFISO(0))  || 
      !equalsIsotropic(OWNISO(0)) || 
      !equalsXonly(OWNDOM(0)) || 
      !isPosDef(OWNTYPE(0)) || 
      C->vdim != SCALAR)
    ERR("anisotropic function not allowed");
	 
  if (C->finiterange == wahr) {
    *natscale = 1.0;
    return;
  }
  
  if (C->inverse!=NULL) { 
    C->inverse(&GLOBAL.gauss.approx_zero, cov, natscale);
    *natscale = 1.0 / *natscale;
    if (ISNAN(*natscale) || *natscale != 0.0) {
      return;
    }
  }
    
  if (GLOBAL.general.naturalscaling != NATSCALE_ORNUMERIC)
    XERR(ERRORRESCALING); 

  if ((C->cov)==nugget)  XERR(ERRORRESCALING); 
  if ( ! HaveSameSystems(PREV, OWN))
    ERR("coordinate system changes not allowed");
     
  // already calculated ?
  //      parami=KAPPA; // do not compare mean,variance, etc.
  //      if (oldcovnr==*covnr) {
  //	for (; parami<=LASTKAPPA; parami++) {
  //	  if (oldp[parami]!=p[parami]) {
  //	    break;
  //	  }
  //	}
  //	if (parami > LASTKAPPA) {
  //	  *natscale=OldNatScale; 
  //	  return;
  //	}
  //      }
  //      for (; parami<=LASTKAPPA; parami++) {oldp[parami]=p[parami];}

  /* **************************************
     Now, find a quick and good solution for NewInvScale --
  */
  
  assert(OWNTOTALXDIM == 1); 
  MultiDimRange(0, cov, natscale);
}



//void UserGetNatScaling(double *natscale) {
//  GetNaturalScaling(KEY()[MODEL_USER], natscale);
//}


void Getxsimugr(coord_type x, double *aniso, int timespacedim, double **xsimugr) {
  // bisher nur fuer Diagonalmatrizen 
  int n, i, w;
  if (aniso == NULL) {
    for(w=0; w<timespacedim; w++) {
      for (i=0; i<3; i++) {
	xsimugr[w][i] = x[w][i];
      }
    } 
  } else {  
    for(n=w=0; w<timespacedim; w++, n+=timespacedim+1) {
      for (i=0; i<3; i++) {
	xsimugr[w][i] = aniso[n] * x[w][i];
      }
    }
  }
}

void TaylorCopy(model *to, model *from) {
  int i, j;
  to->taylorN = from->taylorN;
  to->tailN = from->tailN;
  for(i=0; i<to->taylorN; i++) {
    for (j=0; j<=TaylorPow; j++) to->taylor[i][j] = from->taylor[i][j];
  }
  for(i=0; i<to->tailN; i++) {
    for (j=0; j<=TaylorExpPow; j++) to->tail[i][j] = from->tail[i][j];
  }
}
   
void paramcpy(model *to, model *from, 
	      bool freeing,     // of all the parameters
	      bool force_allocating, // in "c ovcpy" notwendig 
	      bool copy_lists,  // die Unterlisten der LISTOF-Parameter
	      bool recursive, bool copy_mpp) {
  defn *C = DefList + MODELNR(from); // nicht gatternr
  double **pto = to->px,
    **pfrom = from->px;
  int i, v,
     n = UNSET,
     *to_col = to->ncol,
     *to_row = to->nrow,
     *from_col = from->ncol,
     *from_row = from->nrow;

  bool same_model = std::abs(MODELNR(to) - MODELNR(from)) <= 1 ||
    (isDollar(to) && isDollar(from));
  
  if (!same_model) {
    BUG;      
  }

  for (i=0; i<MAXPARAM; i++) {

    if (pfrom[i] == NULL) {
      continue;
    }
   
    if (freeing) {
      PARAMFREE(to, i);
      to_col[i] = from_col[i];
      to_row[i] = from_row[i];
    }
    
    SEXPTYPE type = C->kappatype[i];
    if (type >= LISTOF) {      
      int len = from_row[i];
      listoftype *p = PARAMLIST(from, i);
      if (copy_lists) {      
	listcpy((listoftype **) (pto + i), p, force_allocating);
      } else {
	listpt((listoftype **) (pto + i), p, len, type, force_allocating);
      }
    } else if (isRObject(type)) {
      n = sizeof(sexp_type);
      if (pto[i] == NULL || force_allocating) pto[i] = (double*) MALLOC(n);
      MEMCOPY(pto[i], pfrom[i], n);
      ((sexp_type *) pto[i])->Delete = false;
    } else {
      switch(type) {
      case REALSXP : n = sizeof(double); break;
      case INTSXP : n = sizeof(int); break;
      case STRSXP : n = sizeof(char**); break;
      default : BUG;
      }
      int total = from_row[i] * from_col[i];
      if (pto[i] == NULL || force_allocating) pto[i]=(double*) CALLOC(total, n);
      if (type == STRSXP) {
	for (int k=0; k<total; k++) {
	  char *ptok = PARAMCHAR(to, i)[k],
	    *pfromk = PARAMCHAR(from, i)[k];
	  int bytes = (STRLEN(pfromk) + 1) * sizeof(char);
	  if (ptok != NULL) FREE(ptok);
	  ptok = ((char**) (pto[i]))[k] = (char*) MALLOC(bytes);
	  MEMCOPY(ptok, pfromk, bytes);
	}
      } else MEMCOPY(pto[i], pfrom[i], n * total);
    }
  } // for i

  if (copy_mpp) {
    if (to->mpp.moments < 0 && alloc_mpp_M(to, from->mpp.moments)!=NOERROR) 
      RFERROR("error in allocating memory for Poisson point process data");
    if (to->mpp.moments != from->mpp.moments) BUG;
    assert(sizeof(mpp_properties) == 96);
    mpp_properties *To = &(to->mpp), *From=&(from->mpp);
    assert(To != NULL &&  From != NULL);
    //    To->sum_zhou_c = From->sum_zhou_c;
    //    To->sq_zhou_c = From->sq_zhou_c;
    int vdim = from->vdim[0];
    for (v=0; v<vdim; v++) To->maxheights[v] = From->maxheights[v];
    To->unnormedmass = From->unnormedmass;
    //    To->zhou_c = From->zhou_c;    
    assert(To->mM != NULL && To->mMplus != NULL);
    int nmP1 = To->moments + 1;
    MEMCOPY(To->mM, From->mM, nmP1 * sizeof(double));
    MEMCOPY(To->mMplus, From->mMplus, nmP1 * sizeof(double));

    if (to->qlen != from->qlen) BUG;
    if (from->qlen > 0) {
      assert(to->q != NULL);
      MEMCOPY(to->q, from->q, (to->qlen)* sizeof(double));
    }
  }

  if (recursive) {
    for (i=0; i<MAXSUB; i++) if (from->sub[i] != NULL) {
	assert(to->sub[i] != NULL);
	paramcpy(to->sub[i], from->sub[i], freeing,
		 force_allocating, copy_lists, recursive, copy_mpp);
      }
  }
}


int covcpy(model **localcov, bool sub, model *cov, // err
	   location_type **prevloc, location_type **ownloc,
	   bool copy_lists, bool copy_randomparam, 
	   bool allowCopyingInterface) {
  assert(cov != NULL);
  int i,
    n = UNSET;
  //defn *C = DefList + COVNR; // nicht gatternr

  if ((*localcov = (model*) MALLOC(sizeof(model)))==0)
    RETURN_ERR(ERRORMEMORYALLOCATION);
  model *current = *localcov;

  MEMCOPY(current, cov, sizeof(model)); // replaces COV_NULL(*localcov);
  COV_ALWAYS_NULL(current);
  SET_CALLING_NULL(current, cov);
  paramcpy(current, cov, false, true, copy_lists, false, false);

  if (cov->ownkappanames != NULL) {
    int nkappas = DefList[COVNR].kappas;
    current->ownkappanames = (char**)  CALLOC(nkappas, sizeof(char*));
    for (i=0; i<nkappas; i++) {
      if (cov->ownkappanames[i] != NULL) {
	current->ownkappanames[i] =
	  (char*) MALLOC(sizeof(char) * (1 + STRLEN(cov->ownkappanames[i])));
	STRCPY(current->ownkappanames[i], cov->ownkappanames[i]);
      }
    }
  }
  
  if (cov->q != NULL) {    
    n = sizeof(double) * current->qlen;
    current->q = (double*) MALLOC(n); // QALLOC NOT APPROPRIATE
    MEMCOPY(current->q, cov->q, n);
  } else assert(current->qlen==0);

  current->prevloc = ownloc != NULL ? ownloc : 
    cov->prevloc == prevloc ? prevloc : NULL;
  if (current->prevloc == cov->prevloc && cov->calling==NULL) {
    if (!equalsnowInterface(cov)) {
      BUG;
    }
    if (!allowCopyingInterface) {PRINTF("\n\n***** unallowed copying ******\n"); BUG;}
  }
  
  for (i=0; i<MAXPARAM; i++) {
    int err;
    current->kappasub[i] = NULL;
    if (cov->kappasub[i] == NULL || !copy_randomparam) continue;
    err = covcpy(current->kappasub + i, true, cov->kappasub[i], 
		 prevloc, ownloc, copy_lists, copy_randomparam, false);
    if (err != NOERROR) RETURN_ERR(err);
    SET_CALLING(current->kappasub[i], current);
    
  }
 
  if (sub) {
    for (i=0; i<MAXSUB; i++) {
      int err;
      current->sub[i] = NULL;
      if (cov->sub[i] == NULL) continue;
      err = covcpy(current->sub + i, sub, cov->sub[i], prevloc, ownloc,
		   copy_lists, copy_randomparam, false);
      if (err != NOERROR) RETURN_ERR(err);
      SET_CALLING(current->sub[i], current); 
    }
  } else {
    for (i=0; i<MAXSUB; i++) current->sub[i] = NULL;
  }
  return NOERROR;
}


 
int covcpy(model **localcov, model *cov) { //err
  bool cov2key = &(cov->key)==localcov;
  int 
    err = covcpy(localcov, true, cov, 
		 cov->prevloc, NULL, 
		 false, true, false);//err
  if (err == NOERROR) {
    model *calling = cov2key || cov->calling==NULL ? cov : cov->calling;
    SET_CALLING(*localcov, calling);
  }
  // falls !cov2key && cov->calling == NULL; dann gibt es nur einen Verweis
  // rueckwaerts! Muss aber sein, da sonst LOC_DELETE bei cov->calling==NULL
  // versucht cov->prevloc zu loeschen. oder es muesste prevloc auf NULL
  // gesetzt zu werden. Dies scheint eher contraproduktiv zu sein.
  RETURN_ERR(err);
}

int covcpyWithoutRandomParam(model **localcov, model *cov) {//err
  bool cov2key = &(cov->key)==localcov;
  int 
    err = covcpy(localcov, true, cov, cov->prevloc, NULL, false, false, false);
  if (err == NOERROR) {
    model *calling = cov2key || cov->calling==NULL ? cov : cov->calling;
    SET_CALLING(*localcov, calling);
  }
  // falls !cov2key && cov->calling == NULL; dann gibt es nur einen Verweis
  // rueckwaerts! Muss aber sein, da sonst LOC_DELETE bei cov->calling==NULL
  // versucht cov->prevloc zu loeschen. oder es muesste prevloc auf NULL
  // gesetzt zu werden. Dies scheint eher contraproduktiv zu sein.
  RETURN_ERR(err);
}


int covcpy(model **localcov, model *cov, bool copy_lists) {//err
  bool cov2key = &(cov->key)==localcov;
  int 
    err = covcpy(localcov, true, cov, cov->prevloc, NULL, copy_lists, true,
		 false);
  if (err == NOERROR) {
    model *calling = cov2key || cov->calling==NULL ? cov : cov->calling;
    SET_CALLING(*localcov, calling);
  }
  // falls !cov2key && cov->calling == NULL; dann gibt es nur einen Verweis
  // rueckwaerts! Muss aber sein, da sonst LOC_DELETE bei cov->calling==NULL
  // versucht cov->prevloc zu loeschen. oder es muesste prevloc auf NULL
  // gesetzt zu werden. Dies scheint eher contraproduktiv zu sein.
  RETURN_ERR(err);
}


int covcpy(model **localcov, model *cov, //err
	   double *x, double *T, int spatialdim, int xdimOZ, long lx,
	   bool Time,  bool grid, bool distances) {
  bool cov2key = &(cov->key)==localcov;
  int err;
  location_type **loc = LOCLIST_CREATE(1, xdimOZ + (int) Time);
 model *calling = cov2key || cov->calling==NULL ? cov : cov->calling;

  if ((err = loc_set(x, T, spatialdim, xdimOZ, lx, Time, grid, distances, loc))
      != NOERROR) goto ErrorHandling;
  if ((err = covcpy(localcov, true, cov, loc, NULL, false, true, false))
      != NOERROR) goto ErrorHandling;
  (*localcov)->prevloc = cov->prevloc;
  (*localcov)->ownloc = loc;
   SET_CALLING(*localcov, calling);

 ErrorHandling:
  if (err != NOERROR) LOC_DELETE(&loc);

  RETURN_ERR(err);
}


model *getRemote(model *remotecov, model *rmt, model *target) {
  model *found;
  int i;
  if (rmt == target) return remotecov;
  
  for (i=0; i<MAXPARAM; i++) {
    if (rmt->kappasub[i] != NULL) {
      if (remotecov->kappasub[i] == NULL) BUG;
      found = getRemote(remotecov->kappasub[i], rmt->kappasub[i], target);
      if (found != NULL) return found;
    }
  }
 
  for (i=0; i<MAXSUB; i++) {
    if (rmt->sub[i] != NULL) {
      if (remotecov->sub[i] == NULL) BUG;
      found = getRemote(remotecov->sub[i], rmt->sub[i], target);
      if (found != NULL) return found;
    }
  }
  return NULL;  
}


void Ssetcpy(model *localcov, model *remotecov, model *cov,
	    model *rmt) {
  int i;
  if (cov->Sset != NULL) {
    localcov->Sset = (set_storage*) MALLOC(sizeof(set_storage));
    MEMCOPY(localcov->Sset, cov->Sset, sizeof(set_storage));
    localcov->Sset->remote = getRemote(remotecov, rmt, cov->Sset->remote);
    if (localcov->Sset->remote == NULL) BUG;
  }
  
  for (i=0; i<MAXPARAM; i++) {
    if (cov->kappasub[i] != NULL) {
      if (localcov->kappasub[i] == NULL) BUG;      
      Ssetcpy(localcov->kappasub[i], remotecov, cov->kappasub[i], rmt);
    }
  }
 
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      if (localcov->sub[i] == NULL) BUG;      
      Ssetcpy(localcov->sub[i], remotecov, cov->sub[i], rmt);
    }
  }
}



double *getAnisoMatrix(model *cov, bool null_if_id, int *nrow, int *ncol) {
  int 
    origdim = PrevLoc(cov)->timespacedim;
  if (!isAnyDollar(cov) && null_if_id) { // probably not used anymore
    *nrow = *ncol = origdim;
    return NULL;
  }

  double *ani,
    *aniso = P(DANISO),
    a = PisNULL(DSCALE) ? 1.0 : 1.0 / P0(DSCALE);   
  int i, total,
    dimP1 = origdim + 1;

  if (aniso != NULL) {
    total = origdim * cov->ncol[DANISO];
    long bytes = total * sizeof(double);
    ani = (double *) MALLOC(bytes);
    MEMCOPY(ani, aniso, bytes); 
    for (i=0; i<total; i++) {
      ani[i] *= a;    
    }
    *nrow = cov->nrow[DANISO];
    *ncol = cov->ncol[DANISO];
  } else if (!PisNULL(DPROJ)) {
    int nproj = Nproj,
      *proj = PPROJ;
    total = origdim * nproj;
    ani = (double *) CALLOC(total, sizeof(double));
    for (i=0; i<nproj; i++) {
      ani[i * origdim + proj[i] - 1] = a; 
    }
    *nrow = origdim;
    *ncol = nproj;
  } else {
    if (a == 1.0 && null_if_id) {
      *nrow = *ncol = origdim;
      return NULL;
    }
    total = origdim * origdim;
    ani = (double *) CALLOC(total, sizeof(double));
    for (i=0; i<total; i+=dimP1) ani[i] = a;
    *nrow = *ncol = origdim;
  }

  return ani;
}

double *getAnisoMatrix(model *cov, int *nrow, int *ncol) {
  return getAnisoMatrix(cov, false, nrow, ncol);
}


double GetDiameter(location_type *loc, double *min, double *max,
		   double *center, bool docaniso, bool center_on_loc,
		   int *position) {
  
  // calculates twice the distance between origcenter %*% aniso and
  // all 2^spatialdim corners spanned by min and max; returns maximum distance.

  // NOTE: origcenter is not alway (min + max)/2 since origcenter might be
  //       given by the user and not calculated automatically !

   bool *j=NULL;
  int d,
    origdim = loc->timespacedim,  
    spatialdim = loc->spatialdim;
  double dummy,
    *lx=NULL, *sx=NULL, 
   radiusSq=0.0;
 
  if (loc->grid) {
    double
      *origmin = (double*) MALLOC(origdim * sizeof(double)),
      *origmax = (double*) MALLOC(origdim * sizeof(double)),
      *origcenter = (double*) MALLOC(origdim * sizeof(double));

    for (d=0; d<origdim; d++) {   
      if (loc->xgr[d][XSTEP] > 0) {
	origmin[d] = loc->xgr[d][XSTART];
	origmax[d] = loc->xgr[d][XSTART] + loc->xgr[d][XSTEP] * 
	  (loc->xgr[d][XLENGTH] - 1.0);
      } else {
	origmin[d] = loc->xgr[d][XSTART] + loc->xgr[d][XSTEP] * 
	  (loc->xgr[d][XLENGTH] - 1.0);
	origmax[d] = loc->xgr[d][XSTART];
      }
      if (center_on_loc) {
	double pos = FLOOR(0.5 * loc->xgr[d][XLENGTH]);
	origcenter[d] = origmin[d] + // do not change floor, see below
	  FABS(loc->xgr[d][XSTEP]) * FLOOR(0.5 * loc->xgr[d][XLENGTH]);
	if (position != NULL) position[d] = (int) pos;
      }
      else origcenter[d] = 0.5 * (origmin[d] + origmax[d]);
    }

    if (!docaniso || loc->caniso == NULL) {
      for  (d=0; d<origdim; d++) {   
	center[d] = origcenter[d];
	min[d] = origmin[d];
	max[d] = origmax[d];
	dummy = center[d] - min[d]; // min not max, according to floor above
	radiusSq += dummy * dummy;
      }
    } else { // caniso != NULL
      j = (bool*) MALLOC( (origdim + 1) * sizeof(double));
      lx = (double*) MALLOC(origdim * sizeof(double));
      sx = (double*) MALLOC(spatialdim * sizeof(double));
      
      xA(origcenter, loc->caniso, origdim, spatialdim, center);
      for (d=0; d<origdim; d++) {
	j[d]=false;
	lx[d]=origmin[d];
      }
      j[origdim] = false;
      
      for (d=0; d<spatialdim; d++) {
	min[d] = RF_INF;
	max[d] = RF_NEGINF;
      }
      
      while(true) {
	d=0; 
	while(j[d]) {
	  lx[d]=origmin[d];
	  j[d++]=false;
	}
	if (d==origdim) break;
	j[d]=true;
	lx[d]=origmax[d];
	xA(lx, loc->caniso, origdim, spatialdim, sx);
	
	// suche maximale Distanz von center zu transformierter Ecke
	double distsq = 0.0;
	for (d=0; d<spatialdim; d++) {
	  if (min[d] > sx[d]) min[d] = sx[d];
	  if (max[d] < sx[d]) max[d] = sx[d];
	  dummy = center[d] - sx[d];
	  distsq += dummy * dummy;	  
	}
	if (distsq > radiusSq) radiusSq = distsq;
      }

      UNCONDFREE(j);
      UNCONDFREE(lx);
      UNCONDFREE(sx);
    }
  
    UNCONDFREE(origmin);
    UNCONDFREE(origmax);
    UNCONDFREE(origcenter);

  } else { // not loc->grid
    if (loc->caniso != NULL) BUG;

    double *xx=loc->x; 
    int i,
      endfor = loc->spatialtotalpoints * spatialdim;

    for (d=0; d<spatialdim; d++) {
      min[d]=RF_INF; 
      max[d]=RF_NEGINF; 
    }

    for (i=0; i<endfor; ) {
      for (d=0; d<spatialdim; d++, i++) {
	//temporal part need not be considered, but for ease included#
	if (xx[i] < min[d]) min[d] = xx[i];
	if (xx[i] > max[d]) max[d] = xx[i];
      }
    }
    
    if (loc->Time) {
      assert(spatialdim == origdim - 1);
      if (loc->T[XSTEP] > 0) {
	min[spatialdim] = loc->T[XSTART];
	max[spatialdim] =
	  loc->T[XSTART] + loc->T[XSTEP] * (loc->T[XLENGTH] - 1.0);
      } else {
	min[spatialdim] =
	  loc->T[XSTART] + loc->T[XSTEP] * (loc->T[XLENGTH] - 1.0);
	max[spatialdim] = loc->T[XSTART];
      }
    }

    
    for (radiusSq=0.0, d=0; d<origdim; d++) {
      center[d] = 0.5 * (max[d] + min[d]); 
      dummy = max[d] - center[d];
      radiusSq += dummy * dummy;
    }

    if (center_on_loc) {
      double rSQ,
	minrSQ = RF_INF;
 
      int minx = -999999;
     if (loc->Time) {
	assert(spatialdim == origdim - 1);
	center[spatialdim] = min[spatialdim] + // do not change floor, see below
	  FABS(loc->T[XSTEP]) * FLOOR(0.5 * loc->T[XLENGTH]);
	dummy = center[spatialdim] - min[spatialdim];
      }

      for (i=0; i<endfor; ) {
	for (rSQ = 0.0, d=0; d<spatialdim; d++, i++) {
	  dummy = xx[i] - center[d];
	  rSQ += dummy * dummy;
	}
	if (rSQ < minrSQ) {
	  minrSQ = rSQ;
	  minx = i - spatialdim;
	}
      }
      for (d=0; d<spatialdim; d++) {
	center[d] = xx[minx + d];
	dummy = max[d] - center[d];
	radiusSq += dummy * dummy;
      }
      if (position != NULL) position[0] = minx / spatialdim;
    }
  }    
  return 2.0 * SQRT(radiusSq);
}

double GetDiameter(location_type *loc, double *min, double *max,
		   double *center) {
  return GetDiameter(loc, min, max, center, true, false, NULL);
}


double GetDiameter(location_type *loc) { 
  double diam,
    *dummymin=NULL, *dummymax=NULL, *dummycenter=NULL;
  int origdim = loc->timespacedim; 
  
  dummymin = (double*) MALLOC(origdim * sizeof(double));
  dummymax = (double*) MALLOC(origdim * sizeof(double));
  dummycenter = (double*) MALLOC(origdim * sizeof(double));
  diam = GetDiameter(loc, dummymin, dummymax, dummycenter, true, false, NULL);
  UNCONDFREE(dummymin);
  UNCONDFREE(dummymax);
  UNCONDFREE(dummycenter);
  return diam;
}



bool ok_n(int n, int *f, int nf) // taken from fourier.c of R
{
  int i;
  for (i = 0; i < nf; i++)
    while(n % f[i] == 0) if ((n /= f[i]) == 1) return true;
  return n == 1;
}
int nextn(int n, int *f, int nf) // taken from fourier.c of R
{
    while(!ok_n(n, f, nf)) { n++; }
  return n;
}
#define F_NUMBERS1 3
#define F_NUMBERS2 3
bool HOMEMADE_NICEFFT=false;
unsigned long NiceFFTNumber(unsigned long n) {
  unsigned long i,ii,j,jj,l,ll,min, m=1;
  if (HOMEMADE_NICEFFT) {
    int f[F_NUMBERS1]={2,3,5}; 
    if (n<=1) return n;
    for (i=0; i<F_NUMBERS1; i++) 
      while ( ((n % f[i])==0) && (n>10000)) { m*=f[i]; n/=f[i]; }
    if (n>10000) {
      while (n>10000) {m*=10; n/=10;}
      n++;
    }
    min = 10000000;
    for (i=0, ii=1; i<=14; i++, ii<<=1) { // 2^14>10.000
      if (ii>=n) { if (ii<min) min=ii; break; }
      for (j=0, jj=ii; j<=9; j++, jj*=3) {// 3^9>10.000
	if (jj>=n) { if (jj<min) min=jj; break; }
	
	//for (k=0, kk=jj; k<=6; k++, kk*=5) {// 5^6>10.000
	//if (kk>=n) { if (kk<min) min=kk; break; }
	//for (l=0, ll=kk; l<=5; l++, ll*=7) {// 7^5>10.000
	
	// instead of (if 7 is included)
	for (l=0, ll=jj; l<=6; l++, ll*=5) {// 5^5>10.000
	  	  
	  if (ll>=n) { 
	    if (ll<min) min=ll;
	    break;
	  }
	  //}
	}
      }
    }
    return m*min;
  } else { // not HOMEMADE_NICEFFT
    int f[F_NUMBERS2]={2,3,5}; 
    return nextn(n, f, F_NUMBERS2);
  }
}


void expandgrid(coord_type xgr, double **xx, double* aniso, 
		int olddim, 
		int nrow, // matrix size of old dim
		int ncol  // new dim
		){
  double *x=NULL, * y=NULL; /* current point within grid, but without
		       anisotropy transformation */
  int
    *yi=NULL,   /* counter for the current position in the grid */
    dimM1 = olddim - 1; 
  long pts, w, k, total, n, i, d;
  assert(olddim <= nrow);
  if (aniso == NULL && olddim != ncol) BUG;

  for (pts=1, i=0; i<olddim; i++) pts *= (long int) xgr[i][XLENGTH];

  total = ncol * pts;
  x = *xx = (double*) MALLOC(sizeof(double) * total);
  y = (double*) MALLOC(olddim * sizeof(double));
  yi = (int*) MALLOC(olddim * sizeof(int));


  for (w=0; w<olddim; w++) {y[w]=xgr[w][XSTART]; yi[w]=0;}
  for (k=0; k<total; ){
    if (aniso==NULL) {
      for (d=0; d<ncol; d++, k++) x[k] = y[d];
    } else {
      for (n=d=0; d<ncol; d++, k++, n+=nrow - olddim) {
        x[k] = 0.0;
	for(w=0; w<olddim; w++) {
	  x[k] += aniso[n++] * y[w];
	}
      }
    }
    i = 0;
    (yi[i])++;
    y[i] += xgr[i][XSTEP];

    while(yi[i]>=xgr[i][XLENGTH]) {
      yi[i]=0;
      y[i] = xgr[i][XSTART];
      if (i<dimM1) {
	i++;
	(yi[i])++;
	y[i] += xgr[i][XSTEP];
      } else {
	assert(k==total);
      }
    }
  }
  UNCONDFREE(y);
  UNCONDFREE(yi);
}



void grid2grid(coord_type xgr, double **grani, double *aniso, int origdim, int dim) {
  // diag (+permutations) is assumed, so each row and each col has at most one non-zero
  // component
  double *pgr, *A;
  int d, i,
    origdimM1 = origdim - 1;

  (*grani) = (double *) MALLOC(sizeof(double) * 3 * dim);
  pgr = *grani;

  if (aniso == NULL) {
    for(d=0; d<dim; d++) {
      for (i=0; i<3; i++, pgr++) {
        *pgr = xgr[d][i];
      }
    }
  } else {
    for (d=0; d<dim; d++, pgr += 3) {
      A = aniso + d * origdim;
     for (i=0; i<origdimM1; i++, A++) if (*A != 0.0) break;    
      pgr[XSTART] = xgr[i][XSTART] * *A;
      pgr[XSTEP] = xgr[i][XSTEP] * *A;
      pgr[XLENGTH] = xgr[i][XLENGTH];
    }
  }
}


void xtime2x(double *x, int nx, double *T,
	     double **newx, int timespacedim) {
  double *y, // umhaengen zwingend notwendig, da u.U. **xx und *newx
    // der gleiche aufrufende Pointer ist, d.h. es wird "ueberschrieben"
    *z,  t;
  int j, k, i, d, 
    timelen = (int) T[XLENGTH],
    spatialdim = timespacedim - 1;
  z = *newx = (double*) MALLOC(sizeof(double) * timespacedim * nx * timelen);
  for (k=j=0, t=T[XSTART]; j<timelen; j++, t += T[XSTEP]){
   y = x;
   for (i=0; i<nx; i++) {
       for (d=0; d<spatialdim; d++, y++) {
	z[k++] = *y; 
      }
      z[k++] = t;
    }
  }
}


void xtime2x(double *x, int nx, double *T,
	     double **newx, double *aniso, int nrow, int ncol) {
  double *y, // umhaengen zwingend notwendig, da u.U. **xx und *newx
    // der gleiche aufrufende Pointer ist, d.h. es wird "ueberschrieben"
    *z, dummy, t;
  int j, k, i, d, n, endfor, w,
    spatialdim = nrow - 1,
    timelen = T[XLENGTH],
    nxspdim = nx * spatialdim;

  if (aniso == NULL) {
    assert(nrow == ncol);
    xtime2x(x, nx, T, newx, nrow);
    return;
  }

  z = *newx = (double*) MALLOC(sizeof(double) * ncol * nx * timelen); 
 
  for (k=j=0, t=T[XSTART]; j<timelen; j++, t += T[XSTEP]){
    y = x;
    for (i=0; i<nxspdim; i+=spatialdim) {
      endfor = i + spatialdim;
      for (n=d=0; d<ncol; d++) {
	dummy = 0.0;
	for(w=i; w<endfor; w++) {
	  dummy += aniso[n++] * y[w];
	}
	z[k++] = dummy + aniso[n++] * t; //auf keinen Fall nach vorne holen
      }
    }
  }
}


void x2x(double *x, int nx, double **newx, 
	 double *aniso, int physical_nrow, int nrow, int ncol) {
  double *y = x, // umhaengen zwingend notwendig, da u.U. **xx und *newx
    // der gleiche aufrufende Pointer ist, d.h. es wird "ueberschrieben"
    *z = *newx = (double*) MALLOC(sizeof(double) * ncol * nx);  
  if (aniso == NULL) {
    assert(nrow == ncol);
    MEMCOPY(z, x, sizeof(double) * nx * ncol);
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif
    for (int ii=0; ii<nx; ii++) {
      int i = ii * nrow,
	k = ii * ncol,
	endfor = i + nrow;
      for (int d=0; d<ncol; d++) {
	int n = physical_nrow * d;
        double dummy = 0.0;
	for(int w=i; w<endfor; w++) {
	  dummy += aniso[n++] * y[w];
	}
	z[k++] = dummy; 
      }
    }
  }
}


matrix_type Type(double *M, int nrow, int ncol) {

    // see also analyse_matrix: can it be unified ???

  matrix_type type = TypeMiso; // default
  //  double elmt;
  int i, j, k,
    ncolM1 = ncol - 1,
    endfor = nrow * ncol;

  double *m = M;

  if (m==NULL || (ncol == 1 && nrow == 1)) {
    assert(ncol == nrow);
    return type;
  }
   
  if (ncol > nrow) {
    for (i=ncol * ncol; i < endfor; i++) 
      if (m[i]!= 0.0) return TypeMany;
    ncol = nrow; // !!
  }


  for (k=0; k<ncol; ) {
    for (i=0; i<nrow; i++) if (!R_FINITE(m[i]) || m[i] != 0.0) break;
    for (j=i+1; j<nrow; j++) // mehr als 1 Eintrag != 0
      if (!R_FINITE(m[j]) || m[j] != 0.0) goto TimeInvestigation;

    
    matrix_type newtype;
    if (i!=k && i<nrow) {
      newtype = TypeMproj;
    } else { // entweder i==k oder nur Nullen (dann auch i=k setzbar).
      // Um fall i>=nrow zu vermeiden, wechsel zu k
      newtype = !R_FINITE(M[0]) || !R_FINITE(m[k]) || m[k]!=M[0] ? TypeMdiag
	: TypeMiso;
    } 
    
    if (k < ncolM1) type = newtype > type ? newtype : type;
    else {
      if (type == TypeMtimesep) return i>=nrow-1 ? type : TypeMany;
      if (type == TypeMproj) return i>=nrow-1 ? TypeMtimesepproj : type;
      return newtype > type ? newtype : type;
    }    

    k++;
    m += nrow;
    continue;

  TimeInvestigation: 
    if (k == ncol - 1) return TypeMany;
    type = TypeMtimesep;
    k = ncol - 1;
    m = M + k * nrow;
  }
  return type;
}


 
void TransformLocExt(model *cov,  location_type *loc, bool timesep, 
		     usr_bool gridexpand, // false, true, GRIDEXPAND_AVOID
		     bool same_nr_of_points, // derzeit immer true!!
		     double **grani, double **SpaceTime, 
		     double **caniso, int *Nrow, int *Ncol,//caniso obsolete
		     bool *Time, bool *grid, int *newdim, bool takeX,
		     bool involvedollar) {
  // this fctn transforms the coordinates according to the anisotropy matrix 
  
  location_type *locCani = Loc(cov);
  loc = (loc == NULL) ? locCani : loc;
  assert(locCani->timespacedim == loc->timespacedim);
   bool isdollar = isAnyDollar(cov) && involvedollar;
  matrix_type type;
  int 
    nrow = UNSET,
    ncol = UNSET,
    origdim = locCani->caniso == NULL ? loc->timespacedim : locCani->cani_nrow,
    dim = origdim;
  if (isdollar) {
    if (!PisNULL(DANISO)) dim = cov->ncol[DANISO];
    else if (!PisNULL(DPROJ)) dim = Nproj;
  }
  double *aniso,
    *T = loc->T,
    *x = takeX ? loc->x : loc->y;
  coord_type 
    *xgr= takeX ? &(loc->xgr) : &(loc->ygr);

  
  if (x==NULL && (*xgr)[0] ==NULL) ERR("locations are all NULL");
 
  *newdim = dim;
  *caniso = NULL;
  *Nrow = *Ncol = UNSET;
  *grani = NULL;
  *SpaceTime = NULL;

  if (isdollar) {
    aniso = getAnisoMatrix(cov, true, &nrow, &ncol); 
  } else {
    aniso = NULL;
    nrow = ncol = PrevLoc(cov)->timespacedim;
  }

  if (locCani->caniso != NULL) {    
    if (aniso == NULL) {
      unsigned long
	bytes = sizeof(double) * locCani->cani_nrow * locCani->cani_ncol;
      aniso =(double*) MALLOC(bytes);
      MEMCOPY(aniso, locCani->caniso, bytes);
      nrow = locCani->cani_nrow;
      ncol = locCani->cani_ncol;
    } else {
      double *aniso_old = aniso;
      assert(locCani->cani_ncol == nrow);
      aniso = matrixmult(locCani->caniso, aniso_old, locCani->cani_nrow,
			 nrow, ncol);
      nrow = locCani->cani_nrow;
      UNCONDFREE(aniso_old);
    }
  }

  type = TypeMiso;
  if (aniso != NULL) type = Type(aniso, origdim, dim);
  *Time = loc->Time;
  *grid = false;

  assert(origdim == nrow);
  assert(dim == ncol);
  if (loc->grid) {
    assert(*xgr != NULL);
    if (gridexpand==True || (gridexpand==GRIDEXPAND_AVOID && !isMproj(type))) {
      if (timesep && isMtimesep(type) && *Time) {	
	// space
	//                             not nrow
	expandgrid(*xgr, SpaceTime, aniso, nrow - 1, nrow, ncol - 1);	  
	// time
	grid2grid(*xgr + loc->spatialdim, grani,
		  aniso == NULL ? NULL : aniso + nrow * (ncol-1), 
		  1, 1);
      } else {
	*Time = false;// time is also expanded, if given
	expandgrid(*xgr, SpaceTime, aniso, nrow, nrow, ncol); 
      }
    } else {
      *grid = true;	
      if (isMproj(type) && (!same_nr_of_points || ncol==nrow)) {
	// grid wird multipliziert und/oder umsortiert
	grid2grid(*xgr, grani, aniso, nrow, ncol);
	*Time = false; // no need to have time extra
      } else { // !gridexpand, !isMproj, but still grid
	// z.B. plusmalS.cc falls TBM_INTERN
	// nur aniso auf grid multipliziert
	int i,d,k;
	(*grani) = (double *) MALLOC(sizeof(double) * 3 * origdim);
	for (k=d=0; d<origdim; d++) {
	  for (i=0; i<3; i++) {
	    (*grani)[k++] = (*xgr)[d][i];
	  }
	}
	*caniso = aniso; 
	*Nrow = nrow;
	*Ncol = ncol;
	return; // hier rausgehen --- aniso darf nicht geloescht werden!
      }
    }
  } else { // nogrid
    if (! *Time) { // no grid no time
      x2x(x, loc->spatialtotalpoints, SpaceTime, aniso, nrow, nrow, ncol); 
    } else if (timesep && isMtimesep(type)) {  // no grid, but timesep
      if (same_nr_of_points && ncol!=nrow) { // do not reduce
	x2x(x, loc->spatialtotalpoints, SpaceTime, aniso, nrow, nrow-1, ncol-1);
	grid2grid(&T, grani, aniso==NULL ? NULL : aniso + nrow*ncol - 1, 1, 1);	
      } else if (aniso != NULL && R_finite(aniso[nrow * ncol - 1]) &&
		 aniso[nrow * ncol - 1]==0.0) {//no time comp 
	x2x(x, loc->spatialtotalpoints, SpaceTime, aniso, nrow, nrow - 1, ncol);
	*Time = false;
      } else { // both time and space are left, and time is separated
	x2x(x, loc->spatialtotalpoints, SpaceTime, aniso, nrow, nrow - 1, 
	    ncol - 1);
	grid2grid(&T, grani, aniso==NULL ? NULL : aniso + nrow*ncol - 1, 1, 1);
      }
    } else if (ncol == 1 && type == TypeMproj &&  
	       R_finite(aniso[nrow - 1]) && aniso[nrow - 1] != 0.0 && 
	       !same_nr_of_points) { //only time component left
      assert(nrow != ncol);
      *Time = false;
      *grid = true;
      grid2grid(&T, grani, aniso==NULL ? NULL : aniso + nrow*ncol - 1, 1, 1);
    } else { // no grid, but time, no timesep || not isMtimesep(type)
      xtime2x(x, loc->spatialtotalpoints, T, SpaceTime, aniso, nrow, ncol);
      *Time = false;
    }
  }

  FREE(aniso);
}


void TransformCovLoc(model *cov, bool timesep, usr_bool gridexpand, 
		     bool same_nr_of_points, // in case of projection
		     bool involvedollar) {
  location_type *loc = PrevLoc(cov); // transform2nogrid
  assert(loc != NULL);
  bool Time, grid;
  int err,
    spacedim=UNSET, 
    nrow=UNSET,
    ncol=UNSET;
  double  *xgr=NULL, 
    *x = NULL,
    *caniso = NULL;
  if ((loc->y != NULL && loc->y != loc->x) || 
      (loc->ygr[0] != NULL && loc->ygr[0] != loc->xgr[0])) {
    ERR("unexpected y coordinates");
  }

  assert(cov->prevloc != NULL);
  TransformLocExt(cov, NULL, timesep, gridexpand,  same_nr_of_points, &xgr, &x, 
		      &caniso, &nrow, &ncol, &Time, &grid, &spacedim, true,
		      involvedollar);
  
  if (Time) spacedim--;
  assert(cov->ownloc == NULL);
  if (spacedim > 0) {
    err = loc_set(grid ? xgr : x, 
		  grid ? xgr + 3 * spacedim : xgr, 
		  spacedim, 
		  spacedim,
		  grid ? 3 : loc->spatialtotalpoints, 
		  Time, grid, false, cov);
    assert(grid xor !Loc(cov)->grid);
  } else {
    assert(Time);
    err = loc_set(xgr, NULL, 1, 1, 3, false, true, false, cov);  
  }


 
  // falls not gridexpand und nicht diag bzw. proj
  location_type *ownloc = Loc(cov);
  assert(ownloc->caniso == NULL);
  ownloc->caniso = caniso;
  ownloc->cani_nrow = nrow;
  ownloc->cani_ncol = ncol;
  caniso = NULL;

  FREE(x);
  FREE(xgr);
  if (err != NOERROR) ERR("error when transforming to no grid");
}


void TransformLoc(model *cov, bool timesep, usr_bool gridexpand, 
		  bool involvedollar) {
  TransformCovLoc(cov, timesep, gridexpand, true, involvedollar);
}


//void TransformLocReduce(model *cov, bool timesep, usr_bool gridexpand, 
//			bool involvedollar) {
//  TransformCovLoc(cov, timesep, gridexpand, false, involvedollar);
//}


int TransformLoc(model *cov, location_type *loc, double **xx) {
  bool Time, grid;
  int newdim, nrow, ncol;
  double *caniso = NULL,
    *SpaceTime = NULL;
  assert(!isAnyDollar(cov));
  TransformLocExt(cov, loc, false, True, true, &SpaceTime, xx, 
		  &caniso, &nrow, &ncol, &Time, &grid, &newdim, true, false);
  assert(!grid)
  assert(caniso == NULL && SpaceTime ==NULL);
  return newdim;
}

int TransformLoc(model *cov, double **xx, 
		      bool involvedollar) {
  bool Time, grid;
  int newdim, nrow, ncol;
  double *caniso = NULL,
    *SpaceTime = NULL;
  TransformLocExt(cov, NULL, false, True, true, &SpaceTime, xx, 
		  &caniso, &nrow, &ncol, &Time, &grid, &newdim, true,
		  involvedollar);
  assert(!grid)
  assert(caniso == NULL && SpaceTime ==NULL);
  return newdim;
}


int TransformLoc(model *cov, double **xx, double **yy, 
		      bool involvedollar) {
  location_type *loc = Loc(cov);
  bool Time, grid;
  int newdim, nrow, ncol;
  double *caniso = NULL,
    *SpaceTime = NULL;
  TransformLocExt(cov, NULL, false, True, true, &SpaceTime, xx,
		      &caniso,&nrow, &ncol, &Time, &grid, &newdim, true,
		      involvedollar); 
  assert(!grid)
  assert(caniso == NULL && SpaceTime ==NULL);
  if (loc->y == NULL && loc->ygr[0] == NULL) *yy = NULL;
  else TransformLocExt(cov, NULL, false, True, true, &SpaceTime, yy, 
			   &caniso, &nrow, &ncol, &Time, &grid, &newdim, false,
			   involvedollar);
  assert(!grid)
  assert(caniso == NULL && SpaceTime ==NULL);
  return newdim;
}



double *selectlines(double *m, int *sel, int nsel, int nrow, int ncol) {
    // selects lines in matrix m according to sel
  int j;  
  double *red_matrix,
      *Red_Matrix = (double *) MALLOC(sizeof(double) * nsel * ncol),
      *endfor = Red_Matrix + nsel * ncol;
  
  for (red_matrix=Red_Matrix; red_matrix<endfor; m+=nrow) {
      for (j=0; j<nsel; j++, red_matrix++) {
	  *red_matrix = m[sel[j]];
      }
  }
  return Red_Matrix;
} 


int *selectlines(int *m, int *sel, int nsel, int nrow, int ncol) {
  int j;  
  int *red_matrix,
      *Red_Matrix = (int *) MALLOC(sizeof(int) * nsel * ncol),
      *endfor = Red_Matrix + nsel * ncol;
  
  for (red_matrix=Red_Matrix; red_matrix<endfor; m+=nrow) {
      for (j=0; j<nsel; j++, red_matrix++) {
	  *red_matrix = m[sel[j]];
      }
  }
  return Red_Matrix;
}


bool CallingSet(model *cov) { // is calling model correctly set in 
  // subsequent models
  int i;
  for (i=0; i<cov->nsub; i++) {
    model *sub = cov->sub[i];
    if (sub == NULL) {
      if (DefList[COVNR].range != range_randomcoin) {
	PMI(cov);//
	return false;
      }
    } else {
      if (sub->calling != cov) {
	PRINTF("%d%.50s submodel\n", i+1, TH(i+1));
	PMI(cov);//
	return false;
      }
      if (!CallingSet(sub)) {  
	return false;
      }
    }
  }
    for (; i < MAXSUB; i++) {
      if (cov->sub[i] != NULL) {
	WARN4("%.50s: %d%.50s submodel not NULL although nsub=%d",
		NAME(cov), i+1, TH(i+1), cov->nsub);
	BUG;
      }
    }

 if (cov->key != NULL && !CallingSet(cov->key)) return false;
  if (cov->Splus != NULL && cov->Splus->keys_given) {
    for (i=0; i<cov->nsub; i++) {
      model *sub = cov->Splus->keys[i];
      if (sub == NULL) { 
	if (i > 0 && cov->Splus->keys[i-1] != NULL) { //either all or none NULL!
	  PMI(cov); //
	  return false;
	}
      } else {
	if (sub->calling != cov) {
	  PMI(cov); //
	  return false;
	}
	if (!CallingSet(sub)) return false;
      }
    }
  }
  return true;
}



int get_internal_ranges(model *cov, model *min, model *max, 
			model *pmin, model *pmax, 
			model *openmin, model *openmax) {
  defn *C = DefList + COVNR; // nicht gatternr
  int 
      err = NOERROR, k = 0, i = 0,  // compiler dummy values
      kappas = C->kappas ;
  range_type range;
  rangefct getrange = C->range;
  SEXPTYPE *type = C->kappatype;
  
  if (kappas > 0) {
    getrange(cov, &range); 
    //    printf("***** >>>>>>>>> %d\n", INFDIM);
    assert(maxdim_ok(cov));
      
    err = NOERROR;
    for (i=0; i<kappas; i++) {
      if (isRObject(type[i]) || type[i]==STRSXP)continue;
      int len=cov->ncol[i] * cov->nrow[i];
      double dmin, dmax, dpmin, dpmax, dopenmin, dopenmax,
	value =  RF_NA;
      dpmin = range.pmin[i];
      dpmax = range.pmax[i];
      dmin = range.min[i];     
      dmax = range.max[i];

      //     printf("cov=%s %d\n", NAME(cov), i);
      dopenmin = (double) range.openmin[i];//2006
      dopenmax = (double) range.openmax[i];

      if (type[i] == INTSXP) {
	if (dmin < -MAXINT) {
	  dmin = (double) -MAXINT;
	}
	if (dmax > MAXINT) {
	  dmax = (double) MAXINT;	  
	}	
      }
      for (k=0; k<len; k++) {
	switch(type[i]) {
	case REALSXP :
	  value = P(i)[k];
	  PARAM(min, i)[k] = dmin;
	  PARAM(max, i)[k] = dmax;
	  PARAM(pmin, i)[k] = dpmin;
	  PARAM(pmax, i)[k] = dpmax;
	  PARAM(openmin, i)[k] = dopenmin;
	  PARAM(openmax, i)[k] = dopenmax;
	  break;
	case INTSXP :
	  value = PINT(i)[k] == NA_INTEGER ? RF_NA : (double) PINT(i)[k];
	  PARAMINT(min, i)[k] = dmin;
	  PARAMINT(max, i)[k] = dmax;
	  PARAMINT(pmin, i)[k] = dpmin;
	  PARAMINT(pmax, i)[k] = dpmax;
	  PARAMINT(openmin, i)[k] = dopenmin;
	  PARAMINT(openmax, i)[k] = dopenmax;
	  break;
	case LISTOF + REALSXP : {
	  listoftype *p = PARAMLIST(min, i);
	  if (p->deletelist) {
	    int j, 
	      lenj = p->nrow[k] * p->ncol[k];
	    double
	      *qmin = PARAMLIST(min, i)->lpx[k],
	      *qmax = PARAMLIST(max, i)->lpx[k],
	      *ppmin = PARAMLIST(pmin, i)->lpx[k],
	      *ppmax = PARAMLIST(pmax, i)->lpx[k],
	      *omin = PARAMLIST(openmin, i)->lpx[k],
	      *omax = PARAMLIST(openmax, i)->lpx[k];
	    
	    for (j=0; j<lenj; j++) {
	      qmin[j] = dmin;
	      qmax[j] = dmax;
	      ppmin[j] = dpmin;
	      ppmax[j] = dpmax;
	      omin[j] = dopenmin;
	      omax[j] = dopenmax;
	    }
	  } // elses list had not been copied	  
	  value = RF_NA;
	}
	  break;
	default : 
	  RETURN_ERR(ERRORUNKOWNSXPTYPE);
	}
	
	if (ISNAN(value)) continue;

	dmin = range.min[i];
	dmax = range.max[i];
	if (value < dmin
	    || value > dmax
	    || (range.openmin[i] && value == dmin)
	    || (range.openmax[i] && value == dmax)
	    ) {
	  SERR7("value (%10g) of '%.50s' in '%.50s' not within interval %.50s%10g,%10g%.50s", 
	       value, C->kappanames[i], NICK(cov), 
	       range.openmin[i] ? "(" : "[", dmin, dmax,
	       range.openmax[i] ? ")" : "]"
	       );
	  
	  RETURN_ERR(ERRORM)
	}
      } // for k
    } // for i in kappas;
  } // if (kappa > 0)

  for (i=0; i<MAXPARAM; i++) {
    if (cov->kappasub[i] != NULL) {
      err = get_internal_ranges(cov->kappasub[i], 
	 			min->kappasub[i], max->kappasub[i],  
				pmin->kappasub[i], pmax->kappasub[i],
				openmin->kappasub[i], openmax->kappasub[i]);
      if (err != NOERROR) RETURN_ERR(err);
    }
  }
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      err = get_internal_ranges(cov->sub[i], min->sub[i], max->sub[i], 
			       pmin->sub[i], pmax->sub[i],
			       openmin->sub[i], openmax->sub[i]);
      if (err != NOERROR) RETURN_ERR(err);
    }
  }
  RETURN_NOERROR;
}





void strround(double x, char *s){
  if (x==RF_INF)  SPRINTF(s, "Inf"); else
    if (x==RF_NEGINF)  SPRINTF(s, "-Inf"); else
      if (x == FLOOR(x + 0.5)) SPRINTF(s, "%d", (int) x);
      else SPRINTF(s, "%10g", x);
}
void addmsg(double value, const char *Sign, double y, char *msg) {
  char str1[30], str2[30];
  if ( FABS(value - y) <= 1e-8 *  y) {
    SPRINTF(msg, "%12.12e %.5s %12.12e", value, Sign, y);
  } else {
    strround(value, str1);
    strround(y, str2);
    SPRINTF(msg, "%.50s %.5s %.50s", str1, Sign, str2);
  }
}
/*
void addone(char *str, double x) {
  char str2[30];
  strround(x, str2);
  SPRINTF(str, "%.50s, %.50s", str, str2);
}
void addpair(char *str, double x, double y) {
    if (x == FLOOR(x + 0.5) && y==FLOOR(y + 0.5))
	SPRINTF(str, "%.50s, (%d,%d)", str, (int) x, int(y));
    else 
	SPRINTF(str, "%.50s, (%10g,%10g)", str, x, y);
}
*/


int get_ranges(model *cov, model **min, model **max, 
		model **pmin, model **pmax, 
		model **openmin, model **openmax) {
  int err;
  // returns a reliable value only in the very first
  // entry of each parameter vector/matrix int err;
  if ((err = covcpy(min, cov, true)) != NOERROR) RETURN_ERR(err); 
  if ((err = covcpy(max, cov, true)) != NOERROR) RETURN_ERR(err); 
  if ((err = covcpy(pmin, cov, true)) != NOERROR) RETURN_ERR(err);
  if ((err = covcpy(pmax, cov, true)) != NOERROR) RETURN_ERR(err);
  if ((err = covcpy(openmin, cov, true)) != NOERROR) RETURN_ERR(err); 
  if (( err = covcpy(openmax, cov, true)) != NOERROR) RETURN_ERR(err); 
  SET_CALLING(*min, cov);
  SET_CALLING(*max, cov);
  SET_CALLING(*pmin, cov);
  SET_CALLING(*pmax, cov); 
  SET_CALLING(*openmin, cov);
  SET_CALLING(*openmax, cov);
  //  (*min)->c alling = NULL; ja nicht auf NULL setzen, da sonst angenommen
  //  wird, dass prevloc ein Original ist

  return get_internal_ranges(cov, *min, *max, *pmin, *pmax, *openmin, *openmax);
}



int ReturnOwnField(model *cov) {
   location_type *loc = Loc(cov);
  if (cov->rf != NULL) {
    assert(cov->fieldreturn == wahr && cov->origrf);
    if (cov->origrf) {
      UNCONDFREE(cov->rf);
    }
  } 

  if ((cov->rf = 
       (double*) MALLOC(sizeof(double) * loc->totalpoints * VDIM0))
      == NULL) RETURN_ERR(ERRORMEMORYALLOCATION)
  cov->fieldreturn = wahr;
  cov->origrf = true;
  RETURN_NOERROR;
}


int ReturnOtherField(model *cov, model *which) {
  assert(cov->rf == NULL || cov->rf == which->rf);
  cov->fieldreturn = which->fieldreturn;
  cov->origrf = false;
  cov->rf = which->rf;
  RETURN_NOERROR;
}


void SetLoc2NewLoc(model *cov, location_type **loc) {
  int i,
    maxsub =  DefList[COVNR].maxsub;
  if (cov->ownloc != NULL) return;
  
  for (i=0; i<MAXPARAM; i++) 
    if (cov->kappasub[i] != NULL) SetLoc2NewLoc(cov->kappasub[i], loc);
  
  cov->prevloc = loc;
  for (i=0; i<maxsub; i++) 
    if (cov->sub[i] != NULL) SetLoc2NewLoc(cov->sub[i], loc);
  
  if (cov->key != NULL)  SetLoc2NewLoc(cov->key, loc);
  if (cov->Splus != NULL && cov->Splus->keys_given)
    for (i=0; i<maxsub; i++)
      if (cov->sub[i] != NULL) SetLoc2NewLoc(cov->sub[i], loc);
  if (cov->Sbr != NULL || cov->Sget != NULL || cov->Spgs != NULL ||
      cov->Sset != NULL || cov->Slikelihood != NULL) BUG;
}
  







void set_xdim_intern(system_type *sys, int s, int value) {
  int last = LASTSYSTEM(sys);
  if (s > last) {
    if (s > last + 1)
      RFERROR("improper index found when setting the dimension");
    for (int j=0; j<=s; sys[j++].last = s);
  }
  XDIMi(sys[s]) = value; // OK
  if (s==0) {
    set_cumxmit(sys, s, value);
    s++;
  }
  for (int j=(MAX(s, 1)); j<=last; j++) {
    set_cumxmit(sys, j, CUMXMIT(sys, j-1) + XDIM(sys, j));
  }
}




void set_system(system_type * sys, int idx, int logicaldim, int maxdim, 
		int xdim, Types type, domain_type dom, isotropy_type iso,
		bool check_unset) {

  //  if (idx >= MAXSYSTEMS) { printf("%d %d\n", idx, MAXSYSTEMS); crash(); }  
  assert(idx < MAXSYSTEMS);
  int last = isSetLastSystem(sys) ? LASTSYSTEM(sys) : 0;

  set_logdim(sys, idx, logicaldim);
  set_maxdim(sys, idx, maxdim);
  XDIMi(sys[idx]) = xdim; // OK
  set_type(sys, idx, type);
  set_dom(sys, idx, dom);
  set_iso(sys, idx, iso);

  if (idx >= last) {
    for (int s=0; s<=last; s++) {
      set_last(sys, s, idx);
      if (check_unset && (LOGDIM(sys, s) == UNSET || XDIM(sys, s) == UNSET)) {
	 BUG;
      // maxdim could be still unset !
      }
    }
  }
 
  if (idx==0) {
    set_cumxmit(sys, idx, xdim);
    idx++;
  }
  for (int s=idx; s<=last; s++)
    set_cumxmit(sys, s, CUMXMIT(sys, s-1) + XDIM(sys, s));
}

void set_system(system_type * sys, int idx, int logicaldim, int maxdim, 
		int xdim, Types type, domain_type dom, isotropy_type iso) {
  set_system(sys, idx, logicaldim, maxdim, xdim, type, dom, iso, true);
}


void set_both_systems(model *cov, int idx, int logicaldim, int maxdim, 
		      int xdim, Types type, domain_type dom, isotropy_type iso){
  set_system(OWN, idx, logicaldim, maxdim, xdim, type, dom, iso);
  set_system(PREV, idx, logicaldim, maxdim, xdim, type, dom, iso);
}

void set_both_systems(model *cov, int dim, Types type){
#if MAXSYSTEMS == 1
  if (dim == 1)
    set_both_systems(cov, 0, dim, dim, dim, type, XONLY, ISOTROPIC);
  else if (dim ==2)
    set_both_systems(cov, 0, dim, dim, dim, type, XONLY, DOUBLEISOTROPIC);
  else BUG;
#else
  for (int d=0; d<dim; d++)
    set_both_systems(cov, d, 1, 1, 1, type, XONLY, ISOTROPIC);
#endif
}

void set_system_type(system_type *sys, Types type) {
  int last = LASTSYSTEM(sys);
  if (last == UNSET) BUG;
  set_type(sys, 0, type);
  for (int j=1; j<=last; j++) set_type(sys, j, SameAsPrevType);
}


void set_system_domain(system_type *sys, domain_type dom) {
  int last = LASTSYSTEM(sys);
  if (last == UNSET) BUG;
  for (int j=0; j<=last; j++) set_dom(sys, j, dom);
}

