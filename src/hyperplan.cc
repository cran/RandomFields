/*
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 simulation of a random field by hyperplane tessellation

 Copyright (C) 2001 -- 2011 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <math.h>  
#include <stdio.h>  
//#include <stdlib.h>
 
//#include <string.h>
#include "RF.h"
#include "avltr_modified.h"
//#include <unistd.h>

// #define HYPER_UNIFORM 0   see RF.h
#define HYPER_FRECHET 1
#define HYPER_BERNOULLI 2

#define BLOCKSIZE 1000
#define BLOCKS 1000
typedef double *colour_type[BLOCKS];
typedef unsigned int *code_type[BLOCKS];
typedef double (*randomvar_type)(double p);


double uniform(double p) {return UNIFORM_RANDOM;}
double frechet(double p) { 
  return exp(log(-1.0/log(UNIFORM_RANDOM)) / p);
}
double bernoulli(double p) {
  return (double) (UNIFORM_RANDOM <= p);
}


void hyper_destruct(void **S)
{
  if (*S != NULL) {
    hyper_storage *x; 
    x = *((hyper_storage**)S);
    if (x->aniso != NULL) free(x->aniso); 
    free(*S);
    *S = NULL; 
  }
}

void hyper_NULL(hyper_storage* s) {
    s->aniso = NULL;
}

void SetParamHyperplane(int *action, int *superpos, int *maxlines, 
		        int *mar_distr, double *mar_param)
{
  hyper_param *gp = &(GLOBAL.hyper);
  if (*action) {
    gp->superpos = *superpos;
    gp->maxlines = *maxlines;
    gp->max_distr = *mar_distr;
    gp->mar_param = *mar_param;
  } else {
    *superpos = gp->superpos;
    *maxlines = gp->maxlines;
    *mar_distr = gp->max_distr;
    *mar_param = gp->mar_param; 
  }
}

int init_hyperplane(method_type *meth){

  location_type *loc = meth->loc;
  hyper_storage *s;
  globalparam *gp = meth->gp;
  hyper_param *lp = &(gp->hyper);
  cov_model *cov = meth->cov;
  int dim = cov->tsdim;

  int err,  
    optdim=2; // falls dies gelockert wird, so kc->idx[d] nicht vergessen!
  
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
     hand side or the "right" hand side of a line */


  SET_DESTRUCT(hyper_destruct);
  if ((meth->S=malloc(sizeof(hyper_storage)))==0) {
    err=ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  } 
  s = (hyper_storage*) meth->S;
  hyper_NULL(s);
  if (dim > MAXHYPERDIM) {
    err=ERRORMAXDIMMETH; goto ErrorHandling;
  } 
    
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

  
  Transform2NoGrid(meth, false);
  s->aniso = getAnisoMatrix(meth);

  ERRORMODELNUMBER = -1;	

 
  /****************************************************************/
  /*            determine size of surrounding rectangle           */
  /****************************************************************/
  double origcenter[MAXHYPERDIM], 
      origmin[MAXHYPERDIM], origmax[MAXHYPERDIM],
      min[MAXHYPERDIM], max[MAXHYPERDIM];
  int d;

 
  GetMinMax(meth, origmin, origmax, origcenter, MAXHYPERDIM);
  s->radius = 0.5 * GetDiameter(origcenter, origmin, origmax, 
				s->aniso, loc->timespacedim, 
				loc->timespacedim,
				min, max, s->center);
   for (d=0; d-dim; d++) {
    s->rx[d] = 0.5 * (max[d] - min[d]);
  }

  double *h;

 
  s->hyperplane =  CovList[cov->nr].hyperplane;

  h=NULL;
  if (s->hyperplane == NULL)  {
      err = ERRORFAILED;
      goto ErrorHandling;
  }
  if (s->hyperplane(s->radius, s->center, s->rx, cov, false, &h, &h, &h)
      > lp->maxlines) {
    err = ERRORTOOMANYLINES;
    goto ErrorHandling;
  }

  return NOERROR;
  
 ErrorHandling:
  return err;
}


typedef struct cell_type {
    unsigned int *code;
    res_type colour;
} cell_type;


int cmpcells(void *a, void *b, void *param) {
    cell_type *aa, *bb;
    int *n;
    aa = (cell_type*) a;
    bb = (cell_type*) b;
    n = (int*) param;
    return memcmp(aa->code,  bb->code, *n * sizeof(unsigned int));
}

void delcell(void *a, void *param) {
    cell_type *aa;
    aa = (cell_type*) a;
    free(aa->code);
    free(aa);
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
  
  if ((cell = (cell_type*) malloc(sizeof(cell_type))) == NULL){
      goto ErrorHandling;
  }
  cell->code = NULL;
  if ((cell->code = (unsigned int*) 
       malloc(*integers * sizeof(unsigned int))) == NULL) 
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
      cell->colour = (res_type) randomvar(p);
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
//	  printf("here %d %d \n", 
//		 (int) cell->code[0], (int) lastcell->code[0]);
	  lastcell->colour = (res_type) randomvar(p);
      } else {
	  delcell(cell, NULL); 
      }
  }
  return lastcell;
  
 ErrorHandling:
  if (cell != NULL) {
      if (cell->code != NULL) free(cell->code);
      free(cell);
  }
  return NULL;
}

void do_hyperplane(method_type *meth, res_type *res)
{
  location_type *loc = meth->loc;
  cov_model *cov = meth->cov;
  int dim = cov->tsdim;
  hyper_param *lp = &(meth->gpdo->hyper);

  assert(sizeof(unsigned int) == 4);
 
  double gx, gy, *hx, *hy, *hr, variance,
    E=RF_NAN, 
    sd=RF_NAN;
  int resindex, integers, bits, q, endfor, i, err, j;
  randomvar_type randomvar=NULL;
  hyper_storage *s;
  bool add=FALSE;
  avltr_tree *tree;
  cell_type *cell;

  hx = hy = hr = NULL;
  s = (hyper_storage*) meth->S;
  variance = 
    (cov->nr >= DOLLAR && cov->nr <= LASTDOLLAR) ? cov->p[DVAR][0] : 1.0;
  tree = NULL;
  
  switch (lp->max_distr) {
      case HYPER_UNIFORM : randomvar=uniform; break;
      case HYPER_FRECHET : randomvar=frechet; break;
      case HYPER_BERNOULLI : randomvar=bernoulli; break;
  default : error("random var of unknown type");
  }
  
  switch (meth->simu->distribution) {
      case DISTR_GAUSS : 
	add = true;
	break;
      case DISTR_POISSON : 
	add = true; 
	break;
      case DISTR_MAXSTABLE : 
	add = false;
	break;
      default : 
	error("unknown distribution in hyperplane algorthim\n");
  }

  if (add) for (i=0; i < loc->totalpoints; res[i++]=0.0);
  else // max-stable 
      for (i=0; i < loc->totalpoints; res[i++]= (res_type) R_NegInf);
  /* how many Poisson Hyperplanes maximal (on circle x [0,rmax]) ?  --> p */

  switch (dim) {
      case 1 :
	error("wrong dimension (1) in hyperplane\n");
      case 2 :
	int nn;
	double deltax, deltay;

	deltax = meth->grani[0][XSTEP];
	deltay = meth->grani[1][XSTEP];

	for(nn=0; nn<lp->superpos; nn++){
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
	  if (meth->type <= TypeDiag) {
	    for (gy=meth->grani[1][XSTART], resindex=j=0; j<loc->length[1]; 
		 j++) {
	      for (gx= meth->grani[0][XSTART], i=0; i<loc->length[0]; i++,
		     resindex++) {
//		  printf("\n%f %f\n", gx, gy);  
		if ((cell = determine_cell(gx, gy, hx, hy, hr, &integers,
					   &tree, randomvar, lp->mar_param))
		     == NULL) {
		      err = ERRORMEMORYALLOCATION;
		      goto ErrorHandling;
		  }
// 		printf("%d %d %f %d %d %d\n",
//		       resindex, avltr_count(tree),cell->colour, add, 
//		       integers, q);
		if (add) res[resindex] +=  cell->colour;
		  else if (res[resindex] <  cell->colour)
		      res[resindex] = cell->colour;
		gx += deltax;
	      }
	      gy += deltay;
	    }  
	  } else {
	    for (j=resindex=0; resindex < loc->totalpoints; resindex++) {
	      if ((cell=determine_cell(meth->space[j], meth->space[j+1], 
				       hx, hy, hr,
				       &integers, &tree, randomvar, 
				       lp->mar_param))==NULL){
		  err = ERRORMEMORYALLOCATION;
		  goto ErrorHandling;
	      }
	      if (add) res[resindex] += cell->colour;
	      else if (res[resindex] < cell->colour)
		  res[resindex] = cell->colour;
	      j += dim;
	    }
	  }
	  free(hx); free(hy); free(hr); 
	  hx = hy = hr = NULL;
	  avltr_destroy(tree, delcell);
	  tree = NULL;
	}/* for nn */
	break;
      default: 
	error("wrong dimension (>2) in hyperplane\n"); 
  } // switch  (dim)
  switch (meth->simu->distribution) {
    case DISTR_GAUSS :   
      switch (lp->max_distr) {
        case HYPER_UNIFORM : 
          E = 0.5; 
          sd = 1.0 / 12.0;
          break;
        case HYPER_FRECHET :
          assert(lp->mar_param > 2);
          error("hyper_frechet not programmed yet\n");
          break;
        case HYPER_BERNOULLI : 
          E = lp->mar_param;
          sd = lp->mar_param * (1.0 - lp->mar_param);
          break;
      default : error("distribution unknown in hyperplane\n");
      }
      sd = sqrt(variance / (lp->superpos * sd));
      for(i=0; i<loc->totalpoints; i++) 
	  res[i] = (res_type) (((double) res[i] - lp->superpos * E) * sd);    
      break;
    case DISTR_POISSON : 
     error("Poission not allowed in hyperplane\n");
      break;
    case DISTR_MAXSTABLE : 
      error("Maxstable not allowed in hyperplane\n");
      break;
      default : 
	error("Distribution unknown in hyperplane\n");
  }
  return;

 ErrorHandling: 
//  if (PL>0)
  if (hx != NULL) free(hx);
  if (hy != NULL) free(hy);
  if (hr != NULL) free(hr);
  if (tree!=NULL) avltr_destroy(tree, delcell);
  ErrorMessage(Hyperplane, err);
  error("hyperplane failed\n");
}
                      
		   


       
   
