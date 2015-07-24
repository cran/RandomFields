/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 simulation of a random field by hyperplane tessellation

 Copyright (C) 2001 -- 2015 Martin Schlather, 

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
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USV.
*/

#include <math.h>  
#include <stdio.h>  
//#include <stdlib.h>
 
//#include <string.h>
#include "RF.h"
#include "shape_processes.h"
#include "avltr_modified.h"
#include "Coordinate_systems.h"


//#include <unistd.h>

// #define HYPER_UNIFORM 0   see RF.h

#define BLOCKSIZE 1000
#define BLOCKS 1000
typedef double *colour_type[BLOCKS];
typedef unsigned int *code_type[BLOCKS];
typedef double (*randomvar_type)(double p);

#define HYPER_SUPERPOS (COMMON_GAUSS + 1)
#define HYPER_MAXLINES (COMMON_GAUSS + 2)
#define HYPER_MAR_DISTR (COMMON_GAUSS + 3)
#define HYPER_MAR_PARAM (COMMON_GAUSS + 4)
#define HYPER_ADDITIVE (COMMON_GAUSS + 5)



double uniform(double VARIABLE_IS_NOT_USED p) {return UNIFORM_RANDOM;}
double frechet(double p) { 
  return exp(log(-1.0/log(UNIFORM_RANDOM)) / p);
}
double bernoulli(double p) {
  return (double) (UNIFORM_RANDOM <= p);
}


int check_hyperplane(cov_model *cov) {
 cov_model 
   *key = cov->key,
   *next= cov->sub[0],
   *sub = key != NULL ? key : next;
 int err,
   dim = cov->tsdim
    ; // taken[MAX DIM],
  hyper_param *gp  = &(GLOBAL.hyper);

  ASSERT_CARTESIAN;
  ROLE_ASSERT(ROLE_GAUSS);

  kdefault(cov, HYPER_SUPERPOS, gp->superpos);
  kdefault(cov, HYPER_MAXLINES, gp->maxlines);
  kdefault(cov, HYPER_MAR_DISTR, gp->mar_distr);
  kdefault(cov, HYPER_MAR_PARAM, gp->mar_param);
  kdefault(cov, HYPER_ADDITIVE, true);
  if ((err = checkkappas(cov, false)) != NOERROR) return err;
 

  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
    return ERRORDIM;

  int iso = (cov->calling != NULL && cov->calling->nr == HYPERPLANE_INTERN) ? ISOTROPIC
    : SYMMETRIC;

  if (key == NULL) {
    if ((err = CHECK(next, dim, dim, PosDefType, XONLY, iso, SCALAR, ROLE_COV)) != NOERROR) {
      return err;
    }
  } else { // wenn dies eintritt dann ruft HyperIntern Hyper auf
    cov_model *intern = sub;
    while (intern != NULL && //intern->nr != HYPERPLANE_INTERN && 
	   isAnyDollar(intern)) {
      intern = intern->key != NULL ? intern->key : intern->sub[0];
    }
    if (intern == NULL || intern->nr != HYPERPLANE_INTERN) {
      BUG;
    } else if (intern != cov) 
      paramcpy(intern, cov, true, true, false, false, false);
 
    if ((err = CHECK(sub, dim,  dim, ProcessType, XONLY, CARTESIAN_COORD, 
		       SCALAR, cov->role)) != NOERROR) {
      return err;
    }
  }


  setbackward(cov, sub);
  KAPPA_BOXCOX;
 if ((err = checkkappas(cov)) != NOERROR) return err;

  return NOERROR;
}


int check_hyperplane_intern(cov_model *cov) {  
  ROLE_ASSERT(ROLE_GAUSS);
  assert(cov->key == NULL);

  cov_model *next= cov->sub[0];  
  int err,
    dim = cov->tsdim;    
  hyper_param *gp  = &(GLOBAL.hyper);

  kdefault(cov, HYPER_SUPERPOS, gp->superpos);
  kdefault(cov, HYPER_MAXLINES, gp->maxlines);
  kdefault(cov, HYPER_MAR_DISTR, gp->mar_distr);
  kdefault(cov, HYPER_MAR_PARAM, gp->mar_param);
  kdefault(cov, HYPER_ADDITIVE, true);
 
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) 
    return ERRORDIM;

  if ((err = CHECK(next, dim,  dim, PosDefType, XONLY, UNREDUCED, //ISOTROPIC, UNREDUCED, 
		   SCALAR, ROLE_COV)) != NOERROR) {
    return err;
  }

  if (cov->role == ROLE_GAUSS && next->pref[Hyperplane] == PREF_NONE)
    //      ((next->nr >= DOLLAR && next->nr <= LASTDOLLAR) || 
    //       next->sub[0]->pref[Hyperplane] == PREF_NONE))
    return ERRORPREFNONE;

  setbackward(cov, next);

  return NOERROR;
}


void range_hyperplane(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;

  range->min[HYPER_SUPERPOS] = 1;
  range->max[HYPER_SUPERPOS] = RF_INF;
  range->pmin[HYPER_SUPERPOS] = 100;
  range->pmax[HYPER_SUPERPOS] = 1000;
  range->openmin[HYPER_SUPERPOS] = false;
  range->openmax[HYPER_SUPERPOS] = true; 

  range->min[HYPER_MAXLINES] = 1;
  range->max[HYPER_MAXLINES] = RF_INF;
  range->pmin[HYPER_MAXLINES] = 1;
  range->pmax[HYPER_MAXLINES] = 5000;
  range->openmin[HYPER_MAXLINES] = false;
  range->openmax[HYPER_MAXLINES] = true; 

  range->min[HYPER_MAR_DISTR] = HYPER_UNIFORM;
  range->max[HYPER_MAR_DISTR] = HYPER_BERNOULLI;
  range->pmin[HYPER_MAR_DISTR] = HYPER_UNIFORM;
  range->pmax[HYPER_MAR_DISTR] = HYPER_BERNOULLI;
  range->openmin[HYPER_MAR_DISTR] = false;
  range->openmax[HYPER_MAR_DISTR] = false; 

  range->min[HYPER_MAR_PARAM] = 0;
  range->max[HYPER_MAR_PARAM] = RF_INF;
  range->pmin[HYPER_MAR_PARAM] = 0.01;
  range->pmax[HYPER_MAR_PARAM] = 10;
  range->openmin[HYPER_MAR_PARAM] = false;
  range->openmax[HYPER_MAR_PARAM] = false; 

  range->min[HYPER_ADDITIVE] = 0;
  range->max[HYPER_ADDITIVE] = 1;
  range->pmin[HYPER_ADDITIVE] = 0;
  range->pmax[HYPER_ADDITIVE] = 1;
  range->openmin[HYPER_ADDITIVE] = false;
  range->openmax[HYPER_ADDITIVE] = false; 
}

int struct_hyperplane(cov_model *cov, cov_model VARIABLE_IS_NOT_USED **newmodel) {
  //  PMI(cov, "structhyper");
  if (cov->sub[0]->pref[Hyperplane] == PREF_NONE) {
    return ERRORPREFNONE;
  }
  
  //PMI(cov->calling);

  ROLE_ASSERT(ROLE_GAUSS);
  return NOERROR;
}

int init_hyperplane(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *S){
  cov_model *next = cov->sub[0];
  location_type *loc = Loc(cov);
  hyper_storage *s = NULL;
  long int lines;

  assert(!PisNULL(HYPER_MAXLINES));
  int d,
    maxlines = P0INT(HYPER_MAXLINES),
    dim = cov->tsdim,
    err = NOERROR,  
    optdim=2;// falls dies gelockert wird, so kc->idx[d] nicht vergessen!
  double 
    min[MAXHYPERDIM], max[MAXHYPERDIM],
    *hx = NULL, *hy = NULL, *hz = NULL;

  if (sizeof(unsigned int) != 4)
    SERR("method is written for machines with sizeof(unsigned int) == 4");

  ROLE_ASSERT_GAUSS;
  cov->method = Hyperplane;
 

  if (loc->distances) return ERRORFAILED;
  if (dim > MAXHYPERDIM) {
    err=ERRORMAXDIMMETH; goto ErrorHandling;
  } 
 
  /* n == number of fields superposed 
     lx+1== grid points on x-axis
     ly+1==                y-
     gridlength== grid distance
     lambda== intensity for the Poisson hyperplanes
     res == result (matrix of simulated values)
     normalize==1 then mean=0, var=1 of the random field
     ==0, n>1 then field is devided by sqrt(n)
  */
  
  
  /* cells are coded by a sequence of binary information.
     Each binary unit indicates if the cell is on the "left"
     hand side or the "right" hand side of a line 
  */


  NEW_STORAGE(hyper);
  s = cov->Shyper;
    
  
/****************************************************************/
  /*            Extraction of matching covariances                */
  /****************************************************************/
  if (cov->tsdim == 1) {
    strcpy(ERRORSTRING_OK,"dim=2");
    sprintf(ERRORSTRING_WRONG,
	    "genuine dim=1; this has not been programmed yet.");
    err = ERRORCOVFAILED;
    goto ErrorHandling;
  }

  if (dim > optdim || dim < 1) { 
    err = ERRORWRONGDIM;
    goto ErrorHandling;
  }

  if (!loc->grid) {
    GERR("Hyperplane currently only allows for grids and anisotropies along the axes");
  }
  
  ERRORMODELNUMBER = -1;	

 
  /****************************************************************/
  /*            determine size of surrounding rectangle           */
  /****************************************************************/

  assert(loc->caniso == NULL ||
	 (loc->timespacedim == loc->cani_nrow && 
	  cov->xdimprev == loc->cani_ncol));

  s->radius = 0.5 * GetDiameter(loc, min, max, s->center);

  for (d=0; d-dim; d++) {
    s->rx[d] = 0.5 * (max[d] - min[d]);
  }
 
  s->hyperplane =  CovList[next->nr].hyperplane;

   if (s->hyperplane == NULL)  {
    err = ERRORFAILED;
     goto ErrorHandling;
  }
  lines = s->hyperplane(s->radius, s->center, s->rx, cov /* !! */, 
			    false, &hx, &hy, &hz);

  if (lines > maxlines) {
    GERR("estimated number of lines exceeds hyper.maxlines");  
  } else if (lines < 0) { err = -lines; goto ErrorHandling; }
  
  err = FieldReturn(cov);
  // 
  
ErrorHandling:
  FREE(hx);
  FREE(hy);
  FREE(hz);

  cov->simu.active = err == NOERROR;
  return err;
}



int cmpcells(cell_type *aa, cell_type *bb, int *n) {
    return memcmp(aa->code,  bb->code, *n * sizeof(unsigned int));
}

void delcell(cell_type* aa, int VARIABLE_IS_NOT_USED *n) {
    FREE(aa->code);
    FREE(aa);
}

cell_type *determine_cell(double gx, double gy, double* hx, double* hy, 
			  double* hr, int *integers, avltr_tree **tree, 
			  randomvar_type randomvar, double p)
{ 
    // gx, gy : coordinates
    // hx, hy, hr: list of line coordinates
    // tree: tree for storing already detected cells with at least one
    //       one point in it
    // randomvar: random constant for each cell
    // resindex: numbering of the point with the coordinates (gx, gy)
    // add: maximum or addition of the values for overlaying 
    //      hyperplane tessellations?
    // res: result vector (random field)
  int tt, index, bb;
  unsigned int *cd;
  cell_type *cell;
  static cell_type *lastcell=NULL;
  
  if ((cell = (cell_type*) MALLOC(sizeof(cell_type))) == NULL){
      goto ErrorHandling;
  }
  cell->code = NULL;
  if ((cell->code = (unsigned int*) 
       MALLOC(*integers * sizeof(unsigned int))) == NULL)
      goto ErrorHandling;
 
  cd = cell->code;
  /* calculate the code; if a grid element with */
  /* this code exists, then take this colour    */
  /* (as the two points are then in the same    */
  /* cell) else create a new colour (new cell!) */
  for (tt=index=0; tt<*integers; tt++) {    /* calculate the code... */
    cd[tt] = 0; // of no value
    for (bb=0; bb<32; bb++, index++) {
      cd[tt] <<= 1;
      cd[tt] |= ((hx[index] * gx + hy[index] * gy) < hr[index]);
    }
  }
  if (*tree==NULL) { /* is it the very first point ? */    
      *tree = avltr_create((avl_comparison_func) cmpcells, integers);
      cell->colour = (double) randomvar(p);
      avltr_insert(*tree, cell);
      lastcell = cell; 
  } else { /* search if the calculated code has already appeared (makes
	      sense as it is not the very first point calculated */
      /* as the grid points are visited successively, there is good
	 chance that the previous (calculated) point belongs to the 
	 same cell. Let us check that first! */
      if (memcmp(lastcell->code, cell->code, 
		 *integers * sizeof(unsigned int))
	  && ((lastcell = (cell_type*) *avltr_probe(*tree, cell)) == cell)) {
//	  print("here %d %d \n", 
//		 (int) cell->code[0], (int) lastcell->code[0]);
	  lastcell->colour = (double) randomvar(p);
      } else {
	  delcell(cell, NULL); 
      }
  }
  return lastcell;
  
 ErrorHandling:
  if (cell != NULL) {
     FREE(cell->code);
     UNCONDFREE(cell);
  }
  return NULL;
}

void do_hyperplane(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *S) {
  location_type
    *loc = Loc(cov);
  int    integers, bits, q, endfor, err, len[MAXHYPERDIM],
    vdim = cov->vdim[0],
    dim = cov->tsdim,
    totvdim = loc->totalpoints * vdim,
    superpos = P0INT(HYPER_SUPERPOS),
    mar_distr = P0INT(HYPER_MAR_DISTR);
  double
    gx, gy, 
    variance = isDollar(cov) ? P0(DVAR) : 1.0, 
    *hx = NULL,
    *hy = NULL, 
    *hr = NULL, 
     *res = cov->rf,
    E=RF_NA,
    sd=RF_NA,
    mar_param = P0(HYPER_MAR_PARAM);
  long i, j, resindex;
  randomvar_type randomvar=NULL;
  hyper_storage *s = cov->Shyper;
  avltr_tree *tree = NULL;
  cell_type *cell;
  bool
    additive = (bool) P0INT(HYPER_ADDITIVE);
  SAVE_GAUSS_TRAFO;

  
  switch (mar_distr) {
      case HYPER_UNIFORM : randomvar=uniform; break;
      case HYPER_FRECHET : randomvar=frechet; break;
      case HYPER_BERNOULLI : randomvar=bernoulli; break;
      default : ERR("random var of unknown type");
  }


  if (additive) for (i=0; i < totvdim; res[i++]=0.0);
  else for (i=0; i < totvdim; res[i++]= (double) RF_NEGINF);// max-stable 
  /* how many Poisson Hyperplanes maximal (on circle x [0,rmax]) ?  --> p */
  

  //  warning("return");  return;
 

  switch (dim) {
      case 1 :
      	ERR("wrong dimension (1) in hyperplane\n");
      case 2 :
	int nn;
	double deltax, deltay;

	deltax = loc->xgr[0][XSTEP];
	deltay = loc->xgr[1][XSTEP];
	for (i=0; i<dim; i++) len[i] = (int) loc->xgr[i][XLENGTH];

	for(nn=0; nn<superpos; nn++){
	  q = s->hyperplane(s->radius, s->center, s->rx,
			    cov, true, &hx, &hy, &hr);
	  
	  /* as the length of the codes for the cells are naturally a multiple 
	     of number of bits of an integer variable, some lines are added to
	     the simulated ones in order to get a number of lines that is a 
	     multiple of the number of bits of an integer --- the lines are 
	     placed in 2*rmax, so outside the rectangle */
	  bits = 8 * sizeof(unsigned int);
	  integers = (int) (q / bits);
	  if (integers * bits < q) {
	    integers++;
	    endfor = integers * bits;
	    for(; q < endfor; q++) {
	      hx[q] = hy[q] = 0; 
	      hr[q] = 2.0 * s->radius;
	    } 
	  }

	  /* temporary code */
	  if (isMdiag(Type(loc->caniso, loc->cani_nrow, loc->cani_ncol))) {
	    for (gy=loc->xgr[1][XSTART], resindex=j=0; j<len[1]; 
		 j++) {
	      for (gx= loc->xgr[0][XSTART], i=0; i<len[0]; i++,
		     resindex++) {
//		  print("\n%f %f\n", gx, gy);  
		if ((cell = determine_cell(gx, gy, hx, hy, hr, &integers,
					   &tree, randomvar, mar_param))
		     == NULL) {
		      err = ERRORMEMORYALLOCATION;
		      goto ErrorHandling;
		  }
// 		print("%d %d %f %d %d %d\n",
//		       resindex, avltr_count(tree),cell->colour, add, 
//		       integers, q);
		if (additive) res[resindex] +=  cell->colour;
		  else if (res[resindex] <  cell->colour)
		      res[resindex] = cell->colour;
		gx += deltax;
	      }
	      gy += deltay;
	    }  
	  } else {
	    for (j=resindex=0; resindex < loc->totalpoints; resindex++) {
	      if ((cell=determine_cell(loc->x[j], loc->x[j+1], 
				       hx, hy, hr,
				       &integers, &tree, randomvar, 
				       mar_param))==NULL){
		  err = ERRORMEMORYALLOCATION;
		  goto ErrorHandling;
	      }
	      if (additive) res[resindex] += cell->colour;
	      else if (res[resindex] < cell->colour)
		  res[resindex] = cell->colour;
	      j += dim;
	    }
	  }
	  FREE(hx); FREE(hy); FREE(hr); 
	  avltr_destroy(tree, delcell);
	  tree = NULL;
	}/* for nn */
	break;
      default:	
	ERR("wrong dimension (>2) in hyperplane\n"); 
  } // switch  (dim)


  return;
  if (additive) {
    if (cov->role == ROLE_GAUSS) {
      switch (mar_distr) {
      case HYPER_UNIFORM : 
	E = 0.5; 
	sd = 1.0 / 12.0;
	break;
      case HYPER_FRECHET :
	assert(mar_param > 2);
	NotProgrammedYet("frechet");
	break;
      case HYPER_BERNOULLI : 
	E = mar_param;
	sd = mar_param * (1.0 - mar_param);
	break;
      default : ERR("distribution unknown in hyperplane\n");
      }
      sd = sqrt(variance / (superpos * sd));
      for(i=0; i<loc->totalpoints; i++) 
	res[i] = (double) (((double) res[i] - superpos * E) * sd);    
    
      BOXCOX_INVERSE;
    } else {
      // no standardization up to now -- later alpha-stable ?? TODO
    }
  } else {
    // no standardization up to now -- later max-stable ?? TODO
  }

  return;

 ErrorHandling: 
//  if (PL>0)
  FREE(hx);
  FREE(hy);
  FREE(hr);
  if (tree!=NULL) avltr_destroy(tree, delcell);
  XERR(err); 
  ERR("hyperplane failed\n");
}
                      
		   


       
   
