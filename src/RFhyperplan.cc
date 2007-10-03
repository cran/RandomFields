/*
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 library for unconditional simulation of stationary and isotropic random fields

 Copyright (C) 2001 -- 2006 Martin Schlather, 

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
#include <assert.h>
//#include <string.h>
#include "RFsimu.h"
#include "avltr.h"
//#include <unistd.h>

#define HYPER_UNIFORM 0   // without parameter
#define HYPER_FRECHET 1
#define HYPER_BERNOULLI 2

#define BLOCKSIZE 1000
#define BLOCKS 1000
typedef double *colour_type[BLOCKS];
typedef unsigned int *code_type[BLOCKS];
typedef double (*randomvar_type)();

int HYPERPLANE_SUPERPOS = 300;
int HYPERPLANE_MAXLINES = 1000;
int HYPERPLANE_MAR_DISTR = HYPER_UNIFORM;
double HYPERPLANE_MAR_PARAM = RF_NAN;

double uniform() {return UNIFORM_RANDOM;}
double frechet() { 
  return exp(log(-1.0/log(UNIFORM_RANDOM))/HYPERPLANE_MAR_PARAM);
}
double bernoulli() {
  return (double) (UNIFORM_RANDOM <= HYPERPLANE_MAR_PARAM);
}


void hyper_destruct(void **S)
{
  if (*S != NULL) {
    // hyper_storage *s; s = *((hyper_storage**)S);
    free(*S);
    *S = NULL; 
  }
}

void SetParamHyperplane(int *action, int *superpos, int *maxlines, 
			int *normalise, int *mar_distr, double *mar_param)
{
  if (*action) {
    HYPERPLANE_SUPERPOS = *superpos;
    HYPERPLANE_MAXLINES = *maxlines;
    HYPERPLANE_MAR_DISTR = *mar_distr;
    HYPERPLANE_MAR_PARAM = *mar_param;
  } else {
    *superpos = HYPERPLANE_SUPERPOS;
    *maxlines = HYPERPLANE_MAXLINES;
    *mar_distr = HYPERPLANE_MAR_DISTR;
    *mar_param = HYPERPLANE_MAR_PARAM; 
  }
}

int init_hyperplane(key_type *key, int m)
{
  methodvalue_type *meth; 
  covinfo_type *kc=NULL;
  hyper_storage *s;
  int error, reduceddim, 
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

  
  meth = &(key->meth[m]);
  SET_DESTRUCT(hyper_destruct, m);
  if ((meth->S=malloc(sizeof(hyper_storage)))==0) {
    error=ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  } 
  s = (hyper_storage*) meth->S;
    

  /****************************************************************/
  /*            Extraction of matching covariances                */
  /****************************************************************/
  int v, d, actcov;
  actcov=0;
  for (v=0; v <key->ncov; v++) {
    ERRORMODELNUMBER = v;	
    kc = &(key->cov[v]);
    if (kc->method==Hyperplane && kc->left) {
      cov_fct *cov;
      meth->covlist[actcov] = v;
      cov = &(CovList[kc->nr]);
      assert(kc->nr>=0 && kc->nr<currentNrCov);
      assert(kc->param[VARIANCE] >= 0.0);
      if ((key->ncov>v+1 && kc->op) || 
	  (v>0 && key->cov[v-1].op)) {
	  error=ERRORNOMULTIPLICATION; goto ErrorHandling;
      }
      
      /*    investigation of the param structure and the dimension    */
      /*             check parameter of covariance function           */
      reduceddim = kc->reduceddim;
      if (cov->type==ISOHYPERMODEL || cov->type==ANISOHYPERMODEL) {
	  v += (int) kc->param[HYPERNR];
	  error=ERRORHYPERNOTALLOWED; 
	  goto ErrorHandling;
      }
      if (cov->implemented[Hyperplane] != IMPLEMENTED) { 
	error = ERRORNOTDEFINED;
	goto ErrorHandling;
      }
      if ((error = cov->check(kc->param, reduceddim, key->timespacedim,
			      Hyperplane)) 
	  != NOERROR) {
	ERRORMODELNUMBER = v;	
	goto ErrorHandling;
      }
      if (reduceddim == 1) {
	strcpy(ERRORSTRING_OK,"dim=2");
	sprintf(ERRORSTRING_WRONG,
		"genuine dim=1; this has not been programmed yet.");
	error = ERRORCOVFAILED;
	goto ErrorHandling;
      }
      if (reduceddim > optdim || reduceddim < 1) { 
	error = ERRORWRONGDIM;
	goto ErrorHandling;
      }

     if ((error=Transform2NoGrid(key, v)) != NOERROR) goto ErrorHandling;
     s->hyperplane = cov->hyperplane;
     kc->left = false;
     break;
    }
  } 
  ERRORMODELNUMBER = -1;	
  if (v==key->ncov) { /* no covariance for the considered method found */
    error=NOERROR_ENDOFLIST;
    goto ErrorHandling;
  } else assert(kc!=NULL);
  meth->actcov = 1;
  reduceddim = kc->reduceddim;

  /****************************************************************/
  /*            determine size of surrounding rectangle           */
  /****************************************************************/
  
  GetCenterAndDiameter(key, kc->simugrid, reduceddim, kc->reduceddim,
		       kc->x, kc->aniso, s->center, s->rx, &(s->radius));
  s->radius *= 0.5;
  for (d=0; d<kc->reduceddim; d++) s->rx[d] *= 0.5;

  double *h;
  h=NULL;
  if (s->hyperplane(s->radius, s->center, s->rx, reduceddim, false, &h, &h, &h)
      > HYPERPLANE_MAXLINES) {
    error = ERRORTOOMANYLINES;
    goto ErrorHandling;
  }

  if (key->anisotropy) return NOERROR_REPEAT;
  else return NOERROR;
  
 ErrorHandling:
  return error;
}


typedef struct cell_type {
    unsigned int *code;
    double colour;
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
			  randomvar_type randomvar)
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
      cell->colour = randomvar();
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
	  lastcell->colour = randomvar();
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

void do_hyperplane(key_type *key, int m, double *res)
{
  methodvalue_type *meth; 
  covinfo_type *kc;
  double gx, gy, *hx, *hy, *hr, E, sd, variance;
  int resindex, integers, bits, q, endfor, i,
    xerror, j, reduceddim;
  randomvar_type randomvar;
  hyper_storage *s;
  bool add;
  avltr_tree *tree;
  cell_type *cell;

  hx = hy = hr = NULL;
  assert(sizeof(unsigned int) == 4);
  meth = &(key->meth[m]);
  assert(meth->actcov == 1);
  kc = &(key->cov[meth->covlist[0]]);
  reduceddim = kc->reduceddim;
  s = (hyper_storage*) meth->S;
  assert(meth->actcov == 1);
  variance = key->cov[meth->covlist[0]].param[VARIANCE];
  tree = NULL;
  
  switch (HYPERPLANE_MAR_DISTR) {
      case HYPER_UNIFORM : randomvar=uniform; break;
      case HYPER_FRECHET : randomvar=frechet; break;
      case HYPER_BERNOULLI : randomvar=bernoulli; break;
      default : assert(false);
  }
  
  switch (key->distribution) {
      case DISTR_GAUSS : 
	add = true;
	break;
      case DISTR_POISSON : 
	add = true; 
	break;
      case DISTR_MAXSTABLE : 
	add = false;
	break;
      default : assert(false);
  }

  if (add) for (i=0; i<key->totalpoints; res[i++]=0.0);
  else // max-stable 
    for (i=0; i<key->totalpoints; res[i++]=R_NegInf);
  /* how many Poisson Hyperplanes maximal (on circle x [0,rmax]) ?  --> p */

  switch (reduceddim) {
      case 1 :
	assert(false);
      case 2 :
	int nn;
	double deltax, deltay;

	deltax = kc->xsimugr[XSTEPD[0]];
	deltay = kc->xsimugr[XSTEPD[1]];

	for(nn=0; nn<HYPERPLANE_SUPERPOS; nn++){
	  q = s->hyperplane(s->radius, s->center, s->rx,
			    reduceddim, true, &hx, &hy, &hr);
	  
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
	  if (kc->simugrid) {
	    for (gy=kc->xsimugr[XSTARTD[1]], resindex=j=0; j<key->length[1]; 
		 j++) {
	      for (gx= kc->xsimugr[XSTARTD[0]], i=0; i<key->length[0]; i++,
		     resindex++) {
//		  printf("\n%f %f\n", gx, gy);  
		if ((cell = determine_cell(gx, gy, hx, hy, hr, &integers,
					   &tree, randomvar)) == NULL) {
		      xerror = ERRORMEMORYALLOCATION;
		      goto ErrorHandling;
		  }
// 		printf("%d %d %f %d %d %d\n",
//		       resindex, avltr_count(tree),cell->colour, add, 
//		       integers, q);
		  if (add) res[resindex] += cell->colour;
		  else if (res[resindex] < cell->colour)
		      res[resindex] = cell->colour;
		gx += deltax;
	      }
	      gy += deltay;
	    }  
	  } else {
	    for (j=resindex=0; resindex<key->totalpoints; resindex++) {
	      if ((cell=determine_cell(kc->x[j], kc->x[j+1], hx, hy, hr,
				       &integers, &tree, randomvar))==NULL){
		  xerror = ERRORMEMORYALLOCATION;
		  goto ErrorHandling;
	      }
	      if (add) res[resindex] += cell->colour;
	      else if (res[resindex] < cell->colour)
		  res[resindex] = cell->colour;
	      j += reduceddim;
	    }
	  }
	  free(hx); free(hy); free(hr); 
	  hx = hy = hr = NULL;
	  avltr_destroy(tree, delcell);
	  tree = NULL;
	}/* for nn */
	break;
      default: assert(false);
  } // switch  (reduceddim)
  switch (key->distribution) {
    case DISTR_GAUSS :   
      switch (HYPERPLANE_MAR_DISTR) {
        case HYPER_UNIFORM : 
          E = 0.5; 
          sd = 1.0 / 12.0;
          break;
        case HYPER_FRECHET :
          assert(HYPERPLANE_MAR_PARAM > 2);
          assert(false); 
          break;
        case HYPER_BERNOULLI : 
          E = HYPERPLANE_MAR_PARAM;
          sd = HYPERPLANE_MAR_PARAM * (1.0 - HYPERPLANE_MAR_PARAM);
          break;
      default : assert(false);
      }
      sd = sqrt(variance / (HYPERPLANE_SUPERPOS * sd));
      for(i=0; i<key->totalpoints; i++) 
	res[i] = (res[i] - HYPERPLANE_SUPERPOS * E) * sd;    
      break;
    case DISTR_POISSON : 
      assert(false); 
      break;
    case DISTR_MAXSTABLE : 
      assert(false);
      break;
      default : assert(false);
  }
  return;

 ErrorHandling: 
//  if (GENERAL_PRINTLEVEL>0)
  if (hx != NULL) free(hx);
  if (hy != NULL) free(hy);
  if (hr != NULL) free(hr);
  if (tree!=NULL) avltr_destroy(tree, delcell);
  ErrorMessage(Hyperplane, xerror); assert(false);
}
                      
		   


       
   
