/* 
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 wrapper for Havard Rue's GMRF library

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

// was passiert wenn geo-coef nicht verfuegbar?
// wie wird *problem geloescht ??


#include <math.h>
#include <stdio.h>
 

#include "RF.h"
#ifdef INCLUDEMARKOV
#include "GMRFLib/GMRFLib.h"
#endif


#define MARKOV_DEBUG false
// #define MARKOV_DEBUG true



typedef struct markov_storage {
// typedef cannot be moved to RF.h since
// it requires GMRFLib/GMRFLib.h and the inclusion
// in RF.h leads to conflicts with blas.h definitions
// used in circulant.cc
int rowcol[2];
#ifdef INCLUDEMARKOV
    GMRFLib_problem_tp *problem;
#endif
} markov_storage;


void Markov_destruct(void **S) 
{ 
  if (*S!=NULL) {
#ifdef INCLUDEMARKOV
    markov_storage *s;
    s = *((markov_storage**)S);
    GMRFLib_free_problem(s->problem);
#endif
   free(*S);   
    *S = NULL;
  }
}


void SetParamMarkov(int *action,int *neighbours, double *precision,
		    int *cyclic, int *maxmem) {
  markov_param *lp = &(GLOBAL.markov);
  if (*action) {
    lp->neighbours = *neighbours; 
    if (lp->neighbours < 2) {
	if (PL>0) PRINTF("minimal neighbourhood is 2");
	lp->neighbours = 2;
    } 
    if (lp->neighbours > 3) {
	if (PL>0) PRINTF("maximal neighbourhood is 3");
	lp->neighbours = 3;
    }
    lp->precision = *precision;
    lp->cyclic = *cyclic;
    lp->maxmem = *maxmem;
  } else {
    *neighbours = lp->neighbours; 
    *precision = lp->precision;
    *cyclic = lp->cyclic;
    *maxmem = lp->maxmem;
  }
}



int init_markov(method_type *meth){ 
  SET_DESTRUCT(Markov_destruct);

  int err;
#ifdef INCLUDEMARKOV
  int name, endfor, v, d;
  double param, range, gridstep=0.0, scale=0.0;
  markov_storage* s;
  GMRFLib_geo_problem_tp *geo_problem = NULL;
  int reduceddim reduceddim = meth->reduceddim;
  cov_model *cov=meth->cov;
  cov_fct *cf = &(CovList[cov->nr]);


  if (MARKOV_DEBUG) {
      FILE *fp;
      if ((fp = fopen("available-geo-coefs.dat", "w")))
      {
	  PRINTF("print available geo-coefs to file\n");
	  GMRFLib_print_geo_coefs(fp);
	  fclose(fp); 
      }
  }
     
  SET_DESTRUCT(Markov_destruct);
  if ((meth->S = malloc(sizeof(markov_storage)))==0){
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }  
  s = (markov_storage*) meth->S;
  s->problem = NULL;
  s->rowcol[1] = 1; // falls nur 1-dimensionales Feld

//////////////////////////////////////////////////////////////////////
  err = ERRRORCURRENTLY;
  goto ErrorHandling;
//////////////////////////////////////////////////////////////////////
 
  if (reduceddim > 2) {err=ERRORWRONGDIM; goto ErrorHandling;}     

  if (loc->totalpoints > lp->maxmem) {
      	sprintf(ERRORSTRING_OK, "%d", lp->maxmem);
	sprintf(ERRORSTRING_WRONG,"%ld", loc->totalpoints);
	err = ERRORMARKOVMAXMEMORY;
	goto ErrorHandling;
  } 
 

  err = NOERROR; 
  if (cf->implemented[Markov] != IMPLEMENTED) err=ERRORNOTDEFINED; 
  if (!cov->diag) err=ERRORONLYISOTROPIC; 
  if (err!=NOERROR) goto ErrorHandling;
 
  // current restriction:
  gridstep = meth->xorig[0][XSTEP];
  assert(gridstep != 0.0);
  for (d=1; d<cov->tsdim; d++) {
    if (meth->xorig[d][XSTEP] != gridstep) {
      err=ERRORONLYREGULARGRID; goto ErrorHandling;
    }
  }

  double *aniso;
  aniso = getAnisoMatrix(meth);
  if (meth->type == TypeIso && aniso[0] != 0.0) { 
      scale = 1.0 / aniso[0];
      free(aniso);
  } else {
      free(aniso);
      err=ERRORONLYISOTROPIC; goto ErrorHandling; 
  }
     
  param = 1.0; // dummy 
  name = GMRFLib_CORTP_EXP; // exponential
  range = 3.0;
  if (kc->nr == GAUSS) {
      name =  GMRFLib_CORTP_GAUSS;
      range = sqrt(3.0);
  } else if (kc->nr == WHITTLE) {
#define ONEMINUSEPS 1.0 - 1e-10
    if (gp->general.naturalscaling) {
      GetNaturalScaling(cov, &range, &err);
      error("markov not programmed anymore");
      if (err != NOERROR) goto ErrorHandling;
      range = 1.0 / range;
    } else {
      double xmax, x, xmin, value;
      name = GMRFLib_CORTP, MATERN;
      xmax = 1.0;
      WM(&xmax, cov->p[0][0], &value);
      while (value > 0.05) {
	xmax *= 2.0;
	WM(&xmax, cov->p[0][0], &value);
      }
      xmin = (xmax > 1.0) ? 0.5 * xmax : 0.0;
      while (ONEMINUSEPS > xmin / xmax) {
	x = 0.5 * (xmin + xmax);
	WM(&x, cov->p[0][0], &value);
	if (value > 0.05) xmin = x; 
	else xmax = x;
      }
      range = 1.0 / x;
    }
  } else assert(kc->nr == EXPONENTIAL);
 

  range *= scale / gridstep;
  v = 0;
  for (d=0; d<cov->tsdim; d++) {
    assert(v <= 1);
    s->rowcol[v++] = loc->length[d];
  }
  
  
  GMRFLib_set_error_handler_off();
  if (GMRFLib_is_geo_coefs(name, lp->neighbours, param, range) 
      != GMRFLib_SUCCESS) { 
      sprintf(ERRORSTRING_WRONG,
	      "GMRFLib_is_geo_coefs: err=#%d; param=%f, range=%f", 
	      err, param, range);
      err = ERRORMARKOVPARAMETER;
      goto ErrorHandling;
  }

  err = GMRFLib_init_geo_problem(&geo_problem, name, 
			   lp->neighbours, param, range, 
			   s->rowcol[0], s->rowcol[1], 
			   &lp->precision, lp->cyclic); 
  if (err != GMRFLib_SUCCESS) {
      sprintf(ERRORSTRING_WRONG,
	      "GMRFLib_init_geo_problem: err=#%d; param=%f, range=%f",
	      err, param, range);
      err = ERRORMARKOVPARAMETER;
      goto ErrorHandling;
  }

  GetRNGstate();
  GMRFLib_uniform_init((unsigned int) (4294967296.0 * UNIFORM_RANDOM));
  PutRNGstate();
  err = GMRFLib_init_problem(&(s->problem), NULL, NULL, NULL, NULL, 
				geo_problem->graph,
				geo_problem->Qfunc, geo_problem->Qfunc_arg, 
				NULL, NULL, GMRFLib_NEW_PROBLEM);
  if (err != GMRFLib_SUCCESS) {
      sprintf(ERRORSTRING_WRONG,
	      "GMRFLib_init_problem: err=#%d; param=%f, range=%f", 
	      err, param, range);
      err = ERRORMARKOVPARAMETER;
      goto ErrorHandling;
  }
  err = NOERROR;
  
  if (MARKOV_DEBUG) {
      FILE *fp;
      int i, j;
     /* 
	 print the Q-matrix
      */
      fp = fopen("Q.dat", "w");
      PRINTF("write [Q.dat]\n");    
      for(i = 0; i < geo_problem->graph->n; i++)
      {
	  fprintf(fp, "%d %d %.16f\n", i, i,
		  geo_problem->Qfunc(i, i, geo_problem->Qfunc_arg));
	  for(j = 0; j < geo_problem->graph->nnbs[i];  j++)
	  {
	      fprintf(fp, "%d %d %.16f\n", i, 
		      geo_problem->graph->nbs[i][j],
		      geo_problem -> Qfunc(i,
					      geo_problem->graph->nbs[i][j],
					      geo_problem->Qfunc_arg));
	  } 
      }
      fclose(fp);
  }

  // if (noerror_repeat) return NOERROR_REPEAT;
  // else return err;
  if (geo_problem != NULL) GMRFLib_free_geo_problem(geo_problem);
  return NOERROR_REPEAT;

#else
  err = ERRORMARKOVNOTINCLUDED; goto ErrorHandling;
#endif

 ErrorHandling:
#ifdef INCLUDEMARKOV
  if (geo_problem != NULL) GMRFLib_free_geo_problem(geo_problem);
#endif
  return err;
} 

    

void do_markov(method_type *meth, res_type *res ) {
#ifdef INCLUDEMARKOV
  markov_storage *S;
  long k, i, j;
  int ncol, nrow, idx;
  GMRFLib_problem_tp *problem;

  S = (markov_storage*) meth->S;
  problem = S->problem;
  nrow = S->rowcol[0];
  ncol = S->rowcol[1];
   
  GMRFLib_sample(problem);
  
  k = 0;
  for(j=0; j<ncol; j++) {
      for(i=0; i<nrow; i++, k++) {
	  GMRFLib_lattice2node(&idx, i, j, nrow, ncol);
	  res[k] += (res_type) problem->sample[idx];
      }
  }
  
  if (MARKOV_DEBUG) {//////////////////////
      FILE *fp;
      double *marg_var;
      fp = fopen("sample.dat", "w");
      PRINTF("write [sample.dat]\n");
      for(i=0; i<nrow; i++)
      {
	  for(j=0; j<ncol; j++)
	  {
	      GMRFLib_lattice2node(&idx, i, j, nrow, ncol);
	      fprintf(fp, " %lf", problem->sample[idx]);
	  }
	  fprintf(fp, "\n");
      }
      fclose(fp);
      
      GMRFLib_Qinv(problem, GMRFLib_QINV_ALL);
      fp = fopen("variances.dat", "w");
      PRINTF("write [variances.dat]\n");
      for(i=0; i<nrow; i++)
      {
	  for(j=0; j<ncol; j++)
	  {
	      GMRFLib_lattice2node(&idx, i, j, nrow, ncol);
	      marg_var = GMRFLib_Qinv_get(problem, idx, idx);
	      fprintf(fp, " %lf", *marg_var);
	  }
	  fprintf(fp, "\n");
      }
      fclose(fp);
  } //////////////////////////////
#else 
  assert(false);
#endif
}

// ruwe ??


