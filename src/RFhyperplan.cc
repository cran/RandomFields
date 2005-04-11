/*
 Authors
 Martin Schlather, schlath@hsu-hh.de

 library for unconditional simulation of stationary and isotropic random fields

 Copyright (C) 2001 -- 2004 Martin Schlather, 

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


typedef struct hyper_storage{
  int timespacedim; 
  double deltax[MAXDIM], *x, param[TOTAL_PARAM], rmax, lx[MAXDIM], mx[MAXDIM]; 
  bool grid;
  hyper_pp_fct hyperplane;
} hyper_storage;


void hyper_destruct(void **S)
{
  if (*S != NULL) {
    hyper_storage *s;
    s = *((hyper_storage**)S);
    if (s->x != NULL) free(s->x);
    free(*S);
    *S = NULL; 
  }
}

void SetParamHyperplane(int *action, int *superpos, int *maxlines, 
			int *normalise, int *mar_distr, double *mar_param)
{
  switch (*action) {
  case 0 :
    HYPERPLANE_SUPERPOS = *superpos;
    HYPERPLANE_MAXLINES = *maxlines;
    HYPERPLANE_MAR_DISTR = *mar_distr;
    HYPERPLANE_MAR_PARAM = *mar_param;
   break;
  case 1 :
    *superpos = HYPERPLANE_SUPERPOS;
    *maxlines = HYPERPLANE_MAXLINES;
    *mar_distr = HYPERPLANE_MAR_DISTR;
    *mar_param = HYPERPLANE_MAR_PARAM; 
    if (GetNotPrint) break;
  case 2 : 
    PRINTF("\nHyperplane tessellation\n=======================\nsuperpositions=%d\nmaxlines=%d\nmarginal distribution=%d\nbmarginal parameter=%d\n",
	    HYPERPLANE_SUPERPOS, HYPERPLANE_MAXLINES, HYPERPLANE_MAR_DISTR,
	    HYPERPLANE_MAR_PARAM);
     break;
  default : PRINTF(" unknown action\n"); 
  }
}


int init_hyperplane(key_type *key, int m)
{
  hyper_storage *s;
  int error, start_aniso[MAXDIM];

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
  if ((key->S[m]=malloc(sizeof(hyper_storage)))==0) {
    error=ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  } 
  s = (hyper_storage*)key->S[m];
  s->x = NULL;
    

  /****************************************************************/
  /*            Extraction of matching covariances                */
  /****************************************************************/
  int actcov;
  for (actcov=0; actcov < key->ncov; actcov++) {
    if ((key->method[actcov]==Hyperplane) && (key->left[actcov])) {
      key->left[actcov]=false;
      assert((key->covnr[actcov]>=0) && (key->covnr[actcov]<currentNrCov));
      assert(key->param[actcov][VARIANCE] >= 0.0);
      memcpy(s->param, key->param[actcov], sizeof(double) * key->totalparam);
      if ((key->ncov>actcov+1 && key->op[actcov]) || 
	  (actcov>0 && key->op[actcov-1])) {
        error=ERRORNOMULTIPLICATION; goto ErrorHandling;
      }
      break;
    }
  } 
  if (actcov==key->ncov) { /* no covariance for the considered method found */
    error=NOERROR_ENDOFLIST;
    goto ErrorHandling;
  }
  // do not change actcov !!! used again later on



  /****************************************************************/
  /*    investigation of the param structure and the dimension    */
  /*             check parameter of covariance function           */
  /****************************************************************/
  // determine the reduced dimension of the space
  if (GENERAL_PRINTLEVEL>4) PRINTF("\nchecking parameter structure...");
  bool no_last_comp;
  int index_dim[MAXDIM];
  cov_fct *cov;
  GetTrueDim(key->anisotropy, key->timespacedim, 
	     s->param, // note param is modified, having as very last component
	     //           ==1 if any matrix has very last component <>0
	     &s->timespacedim,
	     &no_last_comp, // vanishing for *all* v=0..actcov-1 ?! 
	     start_aniso, index_dim);
  if (s->timespacedim > 2 || s->timespacedim==1) {
    error = ERRORWRONGDIM;
    goto ErrorHandling;
  } 
  s->grid = key->grid && !key->anisotropy;
  if ((error=Transform2NoGrid(key, s->param, s->timespacedim,
			      start_aniso, &(s->x))) != NOERROR)
    goto ErrorHandling;
  cov = &(CovList[key->covnr[actcov]]);
  if (cov->implemented[Hyperplane] <= NOT_IMPLEMENTED) { 
    error = ERRORNOTDEFINED;
    goto ErrorHandling;
  }
  if (cov->check!=NULL && 
      ((error=cov->check(s->param, s->timespacedim, Hyperplane))) != NOERROR) 
    goto ErrorHandling;
  
  
  /****************************************************************/
  /*            determine size of surrounding rectangle           */
  /****************************************************************/
  int i;
  if (s->grid) {
    for (i=0; i < s->timespacedim; i++) {
      s->mx[i] = s->lx[i] = 0.5 * (key->x[i][XEND] - key->x[i][XSTART]);
      s->deltax[i] = key->x[i][XSTEP] * s->param[INVSCALE];
    }
  } else {
    double min[MAXDIM],max[MAXDIM]; 
    int d, ix, endfor;
    for (d=0; d<MAXDIM; d++) {min[d]=RF_INF; max[d]=RF_NEGINF;}    
    if (key->grid) { // key->grid     
      double sxx[ZWEIHOCHMAXDIM * MAXDIM];
      // unsorted, reduced for param[0...0], #=2^Dim, 
      endfor = 1 << key->timespacedim;   
      // to determine the diameter of the grid determine first 
      // componentwise min and max corner
      GetCornersOfGrid(key, s->timespacedim, start_aniso, s->param, sxx);      
      for (ix=i=0; i<endfor; i++, ix+=s->timespacedim) { 
	for(d=0; d<s->timespacedim; d++) {
	  if (sxx[ix+d] < min[d]) min[d] = sxx[ix+d];
	  if (sxx[ix+d] > max[d]) max[d] = sxx[ix+d];
	}
      }
    } else { // not key->grid
      for (ix=i=0; i<key->totalpoints; i++, ix+=s->timespacedim) {
	if ((GENERAL_PRINTLEVEL>4) && ((i % 10000)==0))
	  PRINTF(" %d [%d]\n",i,key->totalpoints);
	// determine componentwise min and max (for the diameter)
	for(d=0; d<s->timespacedim; d++){//temporal part need not be considered
	  if (s->x[ix+d] < min[d]) min[d] = s->x[ix+d];
	  if (s->x[ix+d] > max[d]) max[d] = s->x[ix+d];
	}
      }
    }
    for (i=0; i < s->timespacedim; i++) {
      s->mx[i] = 0.5 * (min[i] + max[i]);
      s->lx[i] = 0.5 * (max[i] - min[i]);
    }
  } // not s->grid

  s->hyperplane = cov->hyperplane;
  double*h;
  h=NULL;
  if (s->hyperplane(s->lx, s->mx, s->timespacedim, false, &h, &h, &h) >
      HYPERPLANE_MAXLINES) {
    error = ERRORTOOMANYLINES;
    goto ErrorHandling;
  }

  if (key->anisotropy) return NOERROR_REPEAT;
  else return NOERROR;
  
 ErrorHandling:
  return error;
}


#define UPDATE_RES \
   if (add) res[resindex] += colour[nrblock][nrelement]; \
   else if (res[resindex] < colour[nrblock][nrelement]) \
	    res[resindex] = colour[nrblock][nrelement] \

int determine_cell(double gx, double gy, double* hx, double* hy, double* hr,
               int integers, unsigned int *cd, int *Codenr, code_type code,
               colour_type colour, int*Currentblock, randomvar_type randomvar,
	       int resindex, bool add, double *res)
{ 
  int tt, bb, block, element, currentblock, codenr, uppercd, lowercd, blockend,
      lastelement, error, endfor;
  static int nr, nrblock, nrelement, index; // index for the code
  //                                             of the previous point
  bool found,  cd_before_code;
  
  error = NOERROR;
  currentblock = *Currentblock;
  codenr = *Codenr;
  blockend = (BLOCKSIZE - 1) * integers;

  /* calculate the code; if a grid element with */
  /* this code exists, then take this colour    */
  /* (as the two points are then in the same    */
  /* cell) else create a new colour (new cell!) */
  for (tt=index=0; tt<integers; tt++){    /* calculate the code... */
    for (bb=0; bb<32; bb++){
      cd[tt] <<= 1;
      cd[tt] |= ((hx[index] * gx + hy[index] * gy) < hr[index]);
      index++;
    }	
  }
  block = (int) (codenr / BLOCKSIZE);
  element = codenr % BLOCKSIZE;
  assert(block < BLOCKS);

  if (block > currentblock) {
    currentblock++;
    if ((code[currentblock] = 
	 (unsigned int*) malloc(integers * sizeof(int) * BLOCKSIZE))
	== NULL) {error=ERRORMEMORYALLOCATION; goto ErrorHandling;}   
    /* ORDERED list of the codes for all the CELLS */
    if ((colour[currentblock] = 
	 (double*) malloc(sizeof(double) * BLOCKSIZE))
	== NULL) {error=ERRORMEMORYALLOCATION; goto ErrorHandling;} 
    /* the colour for all the cells */
    
  }
  if (codenr==0) { /* is it the very first point ? */   
    for(tt =0; tt<integers; tt++) {code[0][tt]=cd[tt]; }        
    nr = nrblock = nrelement = 0;
    colour[nrblock][nrelement] = randomvar();
    UPDATE_RES;
    codenr++; 
  } else { /* search if the calculated code has already appeared (makes
	      sense as it is not the very first point calculated */
    /* as the grid points are visited successively, there is good
       chance that the previous (calculated) point belongs to the 
       same cell. Let us check that first! */
    found=true; tt=0;
    while (found & (tt<integers)){ /* directly found: point `nr' alright
				      then (still pointing on the previous
				      found cell) */
      found = (cd[tt]==code[nrblock][nrelement * integers+tt]); 
      tt++;
    }
    if (found) {
      UPDATE_RES;
    }
    else { /* search through the ordered list, by the usual 
	      "halfening the intervall [0..nr]"-algorithm*/
      nr=(codenr-1)/2; 
      nrblock = (int) (nr / BLOCKSIZE);
      nrelement = nr % BLOCKSIZE;
      uppercd = codenr - 1;
      lowercd = 0; 
      while (true) {
	for(tt=0; cd[tt]==code[nrblock][nrelement*integers+tt] && 
	      tt<integers; tt++);
	if (tt==integers) {
	  /* here we were successful, the cell already exists */
	  UPDATE_RES;
	  break;
	} 
	cd_before_code = (cd[tt] < code[nrblock][nrelement * integers + tt]);
	/* bad luck, it is definitively a new cell which has to 
	   be created */
	if (cd_before_code) uppercd=nr-1; else lowercd=nr+1;
	if (uppercd < lowercd) {
	  if (!cd_before_code) {
	    nr++; /* now nr shows the position where 
		     the new code should be inserted 
		  */
	    nrblock = (int) nr / BLOCKSIZE;
	    nrelement = nr % BLOCKSIZE;
	  }
	  /* we have to make space for the new code inserted into
	     the ordered list of codes ... */	       	       
	  lastelement = element;
	  if (nrblock < currentblock) {
	    for (bb=block; bb>nrblock; bb--) {
	      for (tt=lastelement * integers - 1; tt>=0; tt--)
		code[bb][tt + integers] = code[bb][tt];
	      for (tt=integers - 1; tt>=0; tt--)
		code[bb][tt] = code[bb-1][blockend + tt];
	      for(tt=lastelement - 1; tt>=0; tt-- ) 
		colour[bb][tt+1]=colour[bb][tt];
	      colour[bb][0] =colour[bb-1][BLOCKSIZE-1];
	      lastelement = BLOCKSIZE;
	    }             
	  }
	  endfor = nrelement * integers;
	  for (tt=lastelement * integers - 1; tt>=endfor; tt--)
	    code[nrblock][tt+integers] = code[nrblock][tt];
	  for (tt=lastelement - 1; tt>=nrelement; tt--)
	    colour[nrblock][tt+1]=colour[nrblock][tt];
	  
	  for (tt=0; tt<integers; tt++) /* ... inserting ... */
	    code[nrblock][nrelement * integers + tt] = cd[tt];
	  colour[nrblock][nrelement] = randomvar();

	  UPDATE_RES;
	  codenr++;

	  if (false)
	  { 
	    int cn, cnblock, cnelement;
	    for (cn=0; cn<codenr; cn++) {
	      cnblock = (int) cn / BLOCKSIZE;
	      cnelement = cn % BLOCKSIZE;
	      for (tt=0; tt<integers; tt++) {
		printf("%u ", code[cnblock][cnelement * integers + tt]);
	      }
	      printf("%f", colour[cnblock][cnelement]);
	      printf("\n");
	    }
	    printf("\n");
	  }

	  break;
	}
        // not uppercd==lowercd
	nr = (uppercd + lowercd) / 2; /* neither successful nor */
	/*          definitely no success, so continuing searching */
	nrblock = (int) (nr / BLOCKSIZE);
	nrelement = nr % BLOCKSIZE;
      } /* while */
    } /* else not directly found */
  } /* } else { not the very first one */
  *Currentblock = currentblock;
  *Codenr = codenr;
  return error;

 ErrorHandling:
  return error;
}

void do_hyperplane(key_type *key, int m, double *res)
{
  double gx, gy, *hx,*hy,*hr, E, sd;
  int codenr, resindex, integers, bits, q, endfor, i, currentblock,
    xerror, j;
  colour_type colour;
  unsigned int *cd;
  code_type code; 
  randomvar_type randomvar;
  hyper_storage *s;
  bool add;
  
  assert(key->active);
  s = (hyper_storage*) key->S[m];
  
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

  if (add) for (i=0; i<key->totalpoints; res[i++]=0);
  else for (i=0; i<key->totalpoints; res[i++]=R_NegInf);
  /* how many Poisson Hyperplanes maximal (on circle x [0,rmax]) ?  --> p */
  
  switch (s->timespacedim) {
      case 1 :
	assert(false);
      case 2 :
	hx = hy = hr = NULL;
	int nn;
	for(nn=0; nn<HYPERPLANE_SUPERPOS; nn++){
	  q = s->hyperplane(s->lx, s->mx, s->timespacedim, true, &hx, &hy, &hr);
	  
	  /* as the length of the codes for the cells are naturally a multiple 
	     of number of bits of an integer variable, some lines are added to
	     the simulated ones in order to get a number of lines that is a 
	     multiple of the number of bits of an integer --- the lines are 
	     placed in 2*rmax, so outside the rectangle */
	  bits = 8 * sizeof(int);
	  integers = (int) (q / bits);
	  if (integers * bits < q) {
	    integers++;
	    endfor = integers * bits;
	    for(; q < endfor; q++) {
	      hx[q] = hy[q] = 0; 
	      hr[q] = 2.0 * (s->lx[0] + s->lx[1]);
	    } 
	  }
	  if ((cd = (unsigned int*) malloc(integers * sizeof(int))) == NULL){
	    xerror=ERRORMEMORYALLOCATION;
	    goto ErrorHandling;
	  }    
	  /* temporary code */
	  
	  codenr = 0;
	  currentblock = -1;
	  if (s->grid) {
	    for (gy=0.0, resindex=j=0; j<key->length[1]; j++) {          
	      for (gx=0.0, i=0; i<key->length[0]; i++){  
		if ((xerror=determine_cell(gx, gy, hx, hy, hr, integers, cd, 
					  &codenr, code, colour, &currentblock, 
					  randomvar, resindex++, add, 
					  res)) != NOERROR)
		goto ErrorHandling;
		gx += s->deltax[0];
	      }
	      gy += s->deltax[1];
	    }  
	  } else {
	    for (j=resindex=0; resindex<key->totalpoints; resindex++) {
	      if ((xerror=determine_cell(s->x[j], s->x[j+1], hx, hy, hr, 
					integers, cd, &codenr, code, colour, 
					&currentblock, randomvar, 
					resindex, add, res)) != NOERROR)
		goto ErrorHandling;
	      j += s->timespacedim;
	    }	 
	  }
	  free(hx); free(hy); free(hr); free(cd);
	  hx = hy = hr = NULL;
	  cd = NULL;
	  for (q=0; q<=currentblock; q++) {
	    free(code[q]); code[q] = NULL; 
	    free(colour[q]); colour[q]=NULL;
	  }
	} /* for nn */
  } // switch  (s->timespacedim)
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
      sd = sqrt(s->param[VARIANCE] / (HYPERPLANE_SUPERPOS * sd));
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
  ErrorMessage(Hyperplane, xerror);
  assert(false);
}
                      
		   


       
   
