/*
Authors
Martin Schlather, schlather@math.uni-mannheim.de

main library for unconditional simulation of random fields

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2014 Martin Schlather
s
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
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
//#include <sys/timeb.h>
 
#include <string.h>
#include "RF.h"
//#include <unistd.h>
#include "Covariance.h"

void genuineStatOwn(cov_model *cov, domain_type *stat, Types *type) {
  cov_model *sub = cov;
  while (isGaussMethod(sub)) sub = sub->sub[0];
  *stat = sub->domown;
  *type = sub->typus;
}

#define STANDARDSTART							\
  cov_model *prev = cov->calling;					\
  assert(prev != NULL);							\
  cov_fct *P = CovList + prev->nr;					\
  if (prev->nr != COVFCTN && prev->nr != VARIOGRAM_CALL &&		\
      P->check != check_cov && prev->nr != COVMATRIX &&			\
      P->check != check_fctn && prev->nr != DIRECT &&			\
      P->check != checkS) BUG;						\
  assert(prev->Spgs != NULL);						\
  pgs_storage *pgs= prev->Spgs;						\
  assert(pgs->x != NULL);						\
  location_type *loc = Loc(cov);					\
  assert(loc != NULL);							\
  bool grid =loc->grid && !loc->Time && loc->caniso==NULL;		\
  bool trafo = (loc->Time || loc->caniso != NULL) && !loc->distances;	\
  bool VARIABLE_IS_NOT_USED dist = loc->distances;			\
  long cumgridlen[MAXMPPDIM +1],					\
    err = NOERROR;							\
  int d,								\
    *gridlen=pgs->gridlen,						\
    *end=pgs->end,							\
    *endy VARIABLE_IS_NOT_USED =pgs->endy,				\
    *start=pgs->start,							\
    *startny VARIABLE_IS_NOT_USED =pgs->startny,			\
    *nx=pgs->nx,							\
    *ny VARIABLE_IS_NOT_USED = pgs->delta,				\
    tot = loc->totalpoints,						\
    VARIABLE_IS_NOT_USED totM1 = tot - 1,				\
    tsdim VARIABLE_IS_NOT_USED = loc->timespacedim,			\
    vdim VARIABLE_IS_NOT_USED = cov->vdim,				\
    VARIABLE_IS_NOT_USED vdimP1 = vdim + 1,				\
    VARIABLE_IS_NOT_USED vdimSq = vdim * vdim,				\
    vdimtot = vdim * tot,						\
    VARIABLE_IS_NOT_USED vdimSqtot = vdim * vdimtot,			\
    VARIABLE_IS_NOT_USED vdimtotSq = vdimtot * tot,			\
    VARIABLE_IS_NOT_USED vdimSqtotSq = vdimtot * vdimtot,		\
    VARIABLE_IS_NOT_USED xdimOZ = loc->xdimOZ, /* weicht u.U. von tsdim bei dist=true ab */ \
    tsxdim = cov->xdimprev; /* weicht u.U. von tsdim bei dist=true ab */ \
  double *x = pgs->x,/* freed only if grid */				\
    *xstart= pgs->xstart,						\
    *inc=pgs->inc,							\
    *y = pgs->supportmin,	/* freed only if grid */		\
    *ystart VARIABLE_IS_NOT_USED =pgs->supportmax,			\
    *yy = NULL,		/* set by Transform2NoGrid */			\
    *xx = NULL,		/* dito			   */			\
    *x0 VARIABLE_IS_NOT_USED = NULL,	/* never free it  */		\
    *incy VARIABLE_IS_NOT_USED =pgs->supportcentre,			\
    *z VARIABLE_IS_NOT_USED = pgs->z;					\
  bool ygiven = loc->y != NULL || loc->ygr[0] != NULL;/* might be changed !*/ \
  if (grid) {								\
    cumgridlen[0] = 1;							\
    for (d=0; d<tsxdim; d++) {						\
      inc[d] = loc->xgr[d][XSTEP];					\
      gridlen[d] = loc->length[d];					\
      cumgridlen[d+1] = gridlen[d] * cumgridlen[d];			\
									\
      start[d] = 0;							\
      end[d] = gridlen[d];						\
      nx[d] = start[d];							\
      x[d] = xstart[d] = loc->xgr[d][XSTART];				\
    }									\
  }									\
  loc->i_col = loc->i_row = 0					
									

// ny wird nur inv CovMatrix verwendet! Stattdessen dort weder
// ystart noch incy!
#define STANDARDSTART_Y_SUPPL						\
  for (d=0; d<tsxdim; d++){						\
    incy[d] = loc->ygr[d]==NULL ? loc->xgr[d][XSTEP] : loc->ygr[d][XSTEP]; \
    y[d] = ystart[d] =							\
      loc->ygr[d]==NULL ? loc->xgr[d][XSTART] : loc->ygr[d][XSTART];	\
    endy[d] = loc->ygr[d]==NULL ? loc->xgr[d][XLENGTH]:loc->ygr[d][XLENGTH]; \
    startny[d] = 0;							\
    ny[d] = startny[d];							\
  }



#define STANDARDINKREMENT	\
  d = 0;			\
  nx[d]++;			\
  x[d] += inc[d];		\
  while (nx[d] >= end[d]) {	\
    nx[d] = start[d];		\
    x[d] = xstart[d];		\
    if (++d >= tsxdim) break;	\
    nx[d]++;			\
    x[d] += inc[d];		\
   }				\
  if (d >= tsxdim) break;		
  // end define StandardInkrement

				
// Achtung hier *nicht* incy, starty
#define STANDARDINKREMENT_Y	\
  d = 0;			\
  ny[d]++;			\
  y[d] += incy[d];		\
   while (ny[d] >= endy[d]) {	\
     ny[d] = startny[d];	\
     y[d] = ystart[d];		\
     if (++d >= tsxdim) break;	\
    ny[d]++;			\
    y[d] += incy[d];		\
   }				\
  // end define StandardInkrement_Y

#define RECYCLE_Y					\
  if (d>tsxdim)						\
    for (d=0; d<tsxdim; d++) {				\
      y[d] = ystart[d] = loc->ygr[d][XSTART];		\
      ny[d] = startny[d];				\
    }
				

// loc->grid bei SelectedCovMatrix notwendig
#define STANDARD_ENDE				\
  if (grid) {					\
    if (xx != NULL) free(xx);			\
    if (yy != NULL) free(yy);			\
  }						\
  loc->i_col = loc->i_row = 0;
  // standard ende

 

void CovVario(cov_model *cov, bool is_cov, bool pseudo, double *value) { 
  STANDARDSTART;
  domain_type domown;
  Types type;
  int i, ii, v, j, m, n;
  bool stat;
  double 
    *C0x = pgs->C0x,
    *C0y = pgs->C0y,
    **Val= pgs->Val,
    *cross=pgs->cross;
  
  // if(isCartesian(prev->isoown) || ygiven);

  genuineStatOwn(cov, &domown, &type);
  stat = cov->domown == XONLY && isNegDef(type);

  if (is_cov) {
    assert(({/*PMI(cov, "Cov/vario"); */ cov->pref[Nothing] != PREF_NONE && isShape(type);})); //
  } else {
    if (cov->pref[Nothing] == PREF_NONE || !isNegDef(type)) {
      assert(({PMI(cov, "cov/Vario"); true;})); //
      ERR("given model is not a variogram");
    }
    if (stat && isCartesian(prev->isoown)) {
      COV(ZERO, cov, C0y);
      
      // PMI(cov->calling);
      //printf("C0y %f %f %f %f\n", C0y[0], C0y[1], C0y[2], C0y[3]);
      //assert(false);
      
    } else {
      if (cov->vdim > 1 && !stat) 
	ERR("multivariate variogram only possible for stationary models");
      NONSTATCOV(ZERO, ZERO, cov, C0y);
    } 
  }
  

   if (cov->vdim > 1) { 
    for (v = i = 0; i<vdimSqtot ; i+=vdimtot) {
      for (ii=i, j=0; j<vdim; ii+=tot, v++, j++) {
	Val[v] = value + ii;
      }
    }    
  }

  //    PMI(cov->calling);


#define UNIVAR COV(x, cov, value + loc->i_row)


#define UNIVAR_Y NONSTATCOV(x, y, cov, value + loc->i_row * vdimSq)

#define MULT {					\
    COV(x, cov, cross);				\
    for (v = 0; v<vdimSq; v++) {		\
      Val[v][loc->i_row] = cross[v];	        \
      /* printf("v=%d\n", v);  // */		\
    }						\
  }

#define MULT_Y	{				\
    NONSTATCOV(x, y, cov, cross);		\
     for (v = 0; v<vdimSq; v++) {		\
      Val[v][loc->i_row] = cross[v];		\
      /* printf("MULT_Y %d %f\n", v, cross[v]); // */ 	\
    }						\
   } 
 
#define VARIO_UNIVAR {							\
    COV(x, cov, cross);							\
    value[loc->i_row] = *C0y - *cross;					\
  }									\

#define VARIO_UNIVAR_Y {						\
    NONSTATCOV(x, x, cov, C0x);						\
    NONSTATCOV(y, y, cov, C0y);						\
    NONSTATCOV(x, y, cov, cross);					\
    value[loc->i_row] =	0.5 * (*C0x + *C0y) - *cross;			\
  }									\
  
  
#define VARIO_MULT {							\
    COV(x, cov, cross);							\
    /* printf("vdim =%d\n", vdim);//	*/				\
    for  (v=0, m=0; m<vdim; m++) {					\
      for (n=0; n<vdim; n++, v++) {					\
       Val[v][loc->i_row] =						\
	 /* 0.5 * (C0x[v] + C0y[v] - cross[m*vdim+n] -cross[n*vdim+m]); */ \
	 C0y[v] - 0.5 * (cross[v] + cross[n*vdim+m]);		\
       /* printf("v=%d %d: %f %f\n", v, n, cross[v], C0y[v]);//	*/	\
     }									\
    }									\
  }									\
  									\

#define PSEUDO_MULT {						\
    COV(x, cov, cross);						\
    for (v=0; v<vdimSq; v++) {						\
      Val[v][loc->i_row] = C0y[v] - cross[v];				\
    }									\
  }									\
  
#define VARIO_MULT_Y {							\
    NONSTATCOV(x, x, cov, C0x);						\
    NONSTATCOV(y, y, cov, C0y);						\
    NONSTATCOV(x, y, cov, cross);					\
   for  (v=m=0; m<vdim; m++) {						\
      for (n=0; n<vdim; n++, v++) {					\
	Val[v][loc->i_row] =						\
	  0.5 * (C0x[v] + C0y[v] - cross[v] -cross[n*vdim+m]);	\
      }									\
    }									\
  }									\

#define PSEUDO_MULT_Y {							\
    NONSTATCOV(x, x, cov, C0x);						\
    NONSTATCOV(y, y, cov, C0y);						\
    NONSTATCOV(x, y, cov, cross);					\
    for (v=0; v<vdimSq; v++) {						\
      Val[v][loc->i_row] = 0.5 * (C0x[v] + C0y[v]) - cross[v];		\
    }									\
  }									\
  
  if (grid) {
    if (ygiven) {							
      STANDARDSTART_Y_SUPPL;
      if (cov->vdim == 1) {
	assert(loc->i_row==0);
	while (true) {
	  if (is_cov) UNIVAR_Y else VARIO_UNIVAR_Y; 
	  loc->i_row++; 
	  STANDARDINKREMENT_Y;
	  RECYCLE_Y;
	  STANDARDINKREMENT; 
	}
      } else { 
	assert(Val != NULL && loc->i_row==0);
	while (true) {
	  if (is_cov) MULT_Y else if (pseudo) PSEUDO_MULT_Y else VARIO_MULT_Y;
	  loc->i_row++; 
	  STANDARDINKREMENT_Y;
	  RECYCLE_Y;
	  STANDARDINKREMENT;
	}
      }
    } else {	// grid, y not given
      if (cov->vdim == 1) { 
	//printf("Val = %d %d\n", Val, loc->i_row);

	assert(loc->i_row==0);
	while (true) {
	  //assert(!is_cov);
	  if (is_cov) UNIVAR else VARIO_UNIVAR;  
	  //printf("%d\n", loc->i_row);
	  //printf("x= %f %d %d v=%f \n", *x,loc->i_row, vdimSq,
	  // (double) value[ loc->i_row * vdimSq]	 );
	  loc->i_row++;  
	  STANDARDINKREMENT; 
	} 
      } else { 
	assert(loc->i_row==0);
	while (true) { 
	  if (is_cov) MULT else if (pseudo) PSEUDO_MULT else VARIO_MULT; 
	  loc->i_row++;
	  STANDARDINKREMENT; 
	} 
      }
    }
  } else { // not a grid

    //   PMI(cov);
    if (trafo) {
      Transform2NoGrid(cov, &xx, &yy);   
      x = xx;
      y = yy;
    } else { x=loc->x; y=loc->y;  }

    if (ygiven) {

      double *y0, *yend;
      y0 = y;
      yend = y + tsxdim * loc->ly;
      if (cov->vdim == 1) { // hier
 	for (loc->i_row=0; loc->i_row<tot; loc->i_row++, x+=tsxdim, y +=tsxdim){
	  if (y >= yend) y = y0;
	  if (is_cov) UNIVAR_Y else VARIO_UNIVAR_Y; 
	  //  printf("tsxdim=%d %f %f, %f %f\n", 
	  //	 tsxdim, x[0], x[1], y[0], y[1], value[loc->i_row]);
	}
	//assert(false);
      } else {
 	for (loc->i_row=0; loc->i_row<tot; loc->i_row++, x+=tsxdim, y +=tsxdim){
	  if (y >= yend) y = y0;
	  // printf("xdimOZ=%d %f %f\n", xdimOZ, x[0], y[0]);
	  //print("tsxdim = %d %d y=%ld y0=%ld end=%ld ly=%d\n", 
	  // tsxdim, loc->i_row, (long int) y, (long int) y0, (long int) yend,
	  //		loc->ly);
	  if (is_cov) MULT_Y else if (pseudo) PSEUDO_MULT_Y else VARIO_MULT_Y;
	}
      }
    } else {

      if (cov->vdim == 1) {
	for (loc->i_row=0; loc->i_row<tot; loc->i_row++, x+=tsxdim) {
	  if (is_cov) UNIVAR else VARIO_UNIVAR;  
	}
      } else {
	for (loc->i_row=0; loc->i_row<tot; loc->i_row++, x+=tsxdim) {
	  if (is_cov) MULT else if (pseudo) PSEUDO_MULT else VARIO_MULT;
	}
      }
    }
  }

  STANDARD_ENDE;
  if (err != NOERROR) XERR(err); 
}



void CovarianceMatrix(cov_model *cov, double *v) {
  domain_type domown;
  Types type;
  genuineStatOwn(cov, &domown, &type);
  if (cov->pref[Nothing] == PREF_NONE ||
      (!isPosDef(type) && (!isNegDef(type) || domown!=XONLY))) {
    // Variogramme sind hier definitiv erlaubt, da durch Addition einer 
    // Konstanten, das Zeug zu einer Kovarianzmatrix gemacht werden kann
    // siehe direct.cc
    assert(({PMI(cov, "cov matrix"); true;})); //
    error("covariance matrix: given model is not a covariance function");
  }
  int l,n,m, VDIM, NEND, NINCR, MINCR, ENDFORINCR;
  double *C = NULL;
  bool vdim_closetogether = GLOBAL.general.vdim_close_together;
  //assert(false);
  
  STANDARDSTART;
  if (grid) {
    STANDARDSTART_Y_SUPPL;
  }

  if (ygiven && (loc->x != loc->y || loc->xgr[0] != loc->ygr[0])) {
    //printf("%d x=%ld y=%ld xgr=%ld  ygr=%ld\n", 
    //	   ygiven, (long int) loc->x, (long int) loc->y, (long int) loc->xgr[0], (long int) loc->ygr[0]);
    // ein Paerchen ist NULL;
    GERR("for the covariance matrix, no y-value may be given");
  }
  
  if (vdim_closetogether) {
    // v-dimensions close together
    VDIM = vdim;
    NEND = vdimSqtot;
    NINCR = vdimtot;
    ENDFORINCR = vdim;
    MINCR = 1;
  } else {
    // values of any single multivariate component close together
    // default in GLOBAL.CovMatrixMulti
    VDIM = 1;
    NEND = vdimSqtotSq;
    NINCR = vdimtotSq;
    ENDFORINCR = vdimtot;
    MINCR = tot;
  }
  
#define MULTICOV							\
  C = v + VDIM * (loc->i_col + loc->i_row * vdimtot);			\
  for (l=n=0; n<NEND; n+=NINCR) {					\
    int endfor = n + ENDFORINCR;					\
    for (m=n; m<endfor; m+=MINCR) {  /* int k= VDIM * (loc->i_col + loc->i_row * vdimtot) + m ; 	printf("MULTICOV %d,%d k=%d z=%d s=%d l=%d %f\n", loc->i_row,loc->i_col, k, k % vdimtot, (int) (k / vdimtot), l, z[l]); assert(R_FINITE(C[0])); // */  \
      C[m] = z[l++];							\
    }									\
  }									\
									\
  if (loc->i_col != loc->i_row) {					\
    C = v + VDIM * (loc->i_row + loc->i_col * vdimtot);			\
    for (l=m=0; m<ENDFORINCR; m+=MINCR) {				\
      for (n=m; n<NEND; n+=NINCR) { /*int k= VDIM * (loc->i_row + loc->i_col * vdimtot) + n ; 	printf("MULTICOV %d,%d k=%d z=%d s=%d l=%d %f\n", loc->i_row,loc->i_col, k, k % vdimtot, (int) (k / vdimtot), l, z[l]); // */ \
	C[n] = z[l++];							\
      }									\
    }									\
  }
  
  
  if (grid) {
    loc->i_col = 0;
    while (true) {
      loc->i_row = loc->i_col;
      for (d=0; d<tsxdim; d++) {
	y[d] = x[d];
	ystart[d] = xstart[d];
	ny[d] = nx[d];
      }   
      while (true) {
	NONSTATCOV(x, y, cov, z);
	// printf("%f %f %f\n", x[0], y[0], z[0]); 
	//printf("xyz (%f, %f), (%f %f) %f\n", x[0], x[1], y[0], y[1], *z);
	//	printf("xyz %f %f %f\n", *x, *y, *z);
	MULTICOV; 
	loc->i_row++;
	STANDARDINKREMENT_Y;
	if (d >= tsxdim) break;
      }
      loc->i_col++;  
      STANDARDINKREMENT;
    }
  } else {
    if (trafo) {
      Transform2NoGrid(cov, &xx);    
      x0 = xx;
    } else x0 = loc->x;
    x = x0;

    for (loc->i_col=0; loc->i_col<tot; loc->i_col++, x+=tsxdim) {
      for (y=x, loc->i_row=loc->i_col; loc->i_row<tot; 
	   loc->i_row++, y+=tsxdim) {
	if (dist) {	  	  
	  // here x and y are ignored	  
	  
	  if (loc->i_row==loc->i_col) {
	    COV(ZERO, cov, z); // !! Achtung; wegen der Mixed models
	    //                    koennen hier verschiedene Werte auftreten
	    //                   aufgrund von loc->i_col und loc->i_row !!
	  } else {
	    COV(x0 + (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2
		      + loc->i_row -1) * tsxdim , cov, z);


	    //	 PMI(cov->calling->calling->calling);
	    //	  printf("loc=%d %d idx=%d %f z=%f\n", loc->i_row, loc->i_col, 
	    //	 (loc->i_col * totM1 - 
	    //	  (loc->i_col * (loc->i_col + 1)) / 2
	    //	  + loc->i_row -1) * tsxdim, 
	    //	 x0[ (loc->i_col * totM1 - 
	    //	      (loc->i_col * (loc->i_col + 1)) / 2
	    //	      + loc->i_row -1) * tsxdim], *z);
	   
	  //	 
	  assert(tsxdim == 1);

	  }

	  	 
	  if (false)
	    if (loc->i_col < 10 && loc->i_row<10) 
	      printf("cm %d %d totM1=%d %d %f %f \n", (int) loc->i_row, (int)loc->i_col, (int) totM1, (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2 + loc->i_row -1) * tsxdim, x0[ (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2 + loc->i_row -1) * tsxdim] , *z); // else assert(false);

	} else {
	  
	  //	  printf("%d %s\n", cov->gatternr, NICK(cov));
	  //assert(cov->gatternr == 5);

	  // PMI(cov->calling);

	  NONSTATCOV(x, y, cov, z);
	  //	 printf("x=%f %f y=%f %f z=%f\n", x[0], x[1], y[0], y[1], z[0]);
	  //PMI(cov);
	  assert(R_FINITE(z[0]));
	}
	MULTICOV;
      }
    }
  }
  

  if (false) {
    for (m=0; m<27; m++) {
      for (n=0; n<27; n++) {
	//printf("%+2.2f ", v[n * 27 + m]);
      }
      //printf("\n"); 
    }
  }
  //assert(false);
  
 ErrorHandling: 
  STANDARD_ENDE;
  if (err!=NOERROR) XERR(err); 
  
  //  int i,j,k;
  //  for (k=i=0; i<tot*tot; i+=tot) {
  //    for (j=0; j<tot;j++) printf("%f ", v[k++]);
  //    printf("\n");  }
  
} // CovarianzMatrix




int ptrStart(int *ptr, int *selected, int nsel, int tot, int vdim) {
  int v,i, step, start, newstart;
  v = start = step = ptr[0] = 0;
  if (selected[ptr[v]] >= step + tot) ptr[v] = -1;
  for (v=1; v<vdim; v++) {
    i = (int) ( (nsel - ptr[v-1]) / (vdim - v + 1) );
    step = tot * v;
    while (i<nsel && selected[i] < step) i++;
    i--;
    while (i>=0 && selected[i] >= step) i--;
    ptr[v] = i+1; ///
    if (ptr[v]>=nsel || selected[ptr[v]] >= step + tot) ptr[v] = -1;
    else {
      newstart = selected[ptr[v]] % tot;
      if (newstart < start) start = newstart;
    }
  }
  
  // for (v=0; v<vdim; v++) print("p%d: %d", v, ptr[v]); 
  //  print(" min=%d\n", start); 
  
  return start;
}


void ptrNext(int *ptr, int *selected, int nsel, int tot, int vdim, int* min) {
  int v,
    step = tot,
    m = *min;
  
  *min = tot; 
  for (v=0; v<vdim; v++, step+=tot) {
    // print("v=%d\n", v);
    if (ptr[v] < 0) continue;
    if (selected[ptr[v]] % tot == m) {
      ptr[v]++;
      if (ptr[v] >= nsel || selected[ptr[v]] >= step) {
	  ptr[v] = -1;
	  continue;
      }
    }
    int newmin = selected[ptr[v]] % tot;
    if (newmin < *min) *min = newmin;
  }
  // for (v=0; v<vdim; v++) print("p%d: %d", v, ptr[v]); 
  //print(" min=%d\n", *min); 
  }

void split(int i, int tsxdim, long int *cum, double *inc, double *x) {
  int j, k;
  for (k=tsxdim-1; k>=0; k--) {
    j = i / cum[k];
    i -= j * cum[k];
    x[k] = (double) j * inc[k];
  }
}

void SelectedCovMatrix(cov_model *cov,
		       int *selected /* nr! */, int nsel, double *v) {
  bool vdim_close_together = GLOBAL.general.vdim_close_together;
  domain_type domown;
  Types type;
  int l,vrow, vcol, 
    oldCMC = RF_NAN, 
    oldCMR = RF_NAN;
    //   vdimselrow, vdimselcol, selrow, selcol, 

  genuineStatOwn(cov, &domown, &type);
  if (cov->pref[Nothing] == PREF_NONE || !isPosDef(type)) {
    assert(({PMI(cov, "cov matrix"); true;})); //
    error("cov. matrix: given model is not a covariance function");
  }
  if (vdim_close_together) {
    error("vdim_closetogether not programmed yet for selected access on covariance matrices");
  }
  
   
  STANDARDSTART;
  int *ptrcol = pgs->ptrcol, 
    *ptrrow = pgs->ptrrow;

  oldCMC = loc->i_col;
  oldCMR = loc->i_row;
  
  if (ygiven && (loc->x != loc->y || loc->xgr[0] != loc->ygr[0])) {
    // ein Paerchen ist NULL;
    //printf("%d x=%ld y=%ld xgr=%ld  ygr=%ld\n", 
    //	   ygiven, (long int) loc->x, (long int) loc->y, (long int) loc->xgr[0], (long int) loc->ygr[0]);
    GERR("for the covariance matrix, no y-value may be given");
  }
  
  if (!grid) {
    if (trafo) {
      Transform2NoGrid(cov, &xx);    
      x0 = xx;
    }
    else x0 = loc->x;
  }
  
  loc->i_col = ptrStart(ptrcol, selected, nsel, tot, vdim);
  
  while (loc->i_col < tot) {
    //  print("col=%d\n", loc->i_col);
    if (grid) split(loc->i_col, tsxdim, cumgridlen, inc, x);
    else x = x0 + tsxdim * loc->i_col; 
    MEMCOPY(ptrrow, ptrcol, sizeof(int) * vdim);
    loc->i_row=loc->i_col;
    while (loc->i_row < tot){
      // print("row=%d\n", loc->i_row);
      // print("col=%d row=%d\n", loc->i_col, loc->i_row);
      if (dist) {
	if (loc->i_row==loc->i_col) {
	  COV(ZERO, cov, z);  // Achtung aufgrund loc->i_row  und loc->i_col
	  // koennen hier verschiedene Werte auftreten (mixed models) !!
	  } else {
	  COV(x0 + 
	      (loc->i_col * totM1 - 
	       (loc->i_col * (loc->i_col + 1)) / 2
	       + loc->i_row -1) * tsxdim , cov, z);
	  
	}
      } else {
	//	assert(false);
	if (grid) split(loc->i_row, tsxdim, cumgridlen, inc, y);
	else y = x0 + tsxdim * loc->i_row; 
	NONSTATCOV(x, y, cov, z);
      }
      // print("hello\n");
      
      for (vcol=l=0; vcol<vdim; vcol++) {
	//
	//	print(">> %d %d sel=%d tot=%d i=%d %d j=%d\n",
	//	      vcol, ptrcol[vcol], selected[ptrcol[vcol]],
	//	      tot, loc->i_col, selected[ptrcol[vcol]] % tot == loc->i_col,
	//     loc->i_row);
	
	if (ptrcol[vcol] >= 0 && 
	    selected[ptrcol[vcol]] % tot == loc->i_col) {
	  for (vrow=0; vrow<vdim; vrow++) {
	    //	    print("vcol=%d, vrow=%d %d %d j=%d %d z=%f, %d (%d %d)\n",
	    //		  vcol, vrow, ptrrow[vrow],
	    //	   	  selected[ptrrow[vrow]], loc->i_row, 
	    //		  selected[ptrrow[vrow]] % tot == loc->i_row,
	    //	  z[l], l,
	    //		  ptrcol[vcol] * nsel + ptrrow[vrow],
	    //	  ptrcol[vcol] + ptrrow[vrow] * nsel );
	    if (ptrrow[vrow] >= 0 && 
		selected[ptrrow[vrow]] % tot == loc->i_row) {
	      v[ptrcol[vcol] * nsel + ptrrow[vrow]] = 
		v[ptrcol[vcol] + ptrrow[vrow] * nsel] = z[l++];
	      //   print("%d %d %d %d \n", loc->i_col, loc->i_row,
	      //          ptrcol[vcol], ptrrow[vrow]);
	    }
	  } // for vrow
	}
      } // for vcol
      ptrNext(ptrrow, selected, nsel, tot, vdim, &loc->i_row);
      } // while loc->i_row < tot
    ptrNext(ptrcol, selected, nsel, tot, vdim, &loc->i_col);  
      // assert( loc->i_col < 5);
  } // while loc->i_col < tot
  
  if (cov->domown==XONLY && isNegDef(cov->typus) && 
      !isPosDef(cov->typus)) {
    double first=v[0];
    for (l=0; l<vdimSqtotSq; v[l++] -= first);
  }

  
 ErrorHandling:
  STANDARD_ENDE; // 	ErrorHandling:
  
  loc->i_col=oldCMC;
  loc->i_row=oldCMR; 
  
  if (err!=NOERROR) XERR(err); 
    
}


void partial_loc_set_matrix(cov_model *cov, double *x, int lx, bool dist,
			    bool grid){
  location_type *loc = Loc(cov);
  int err;
  bool ygiven = !dist && loc->ly != 0;
  if ((err = partial_loc_set(loc, x, ygiven ? x : NULL,
			     lx, ygiven ? lx : 0, dist,
			     loc->xdimOZ,
			     NULL, grid, false)) 
	!= NOERROR) XERR(err);
}

void partial_loc_set_matrixOZ(cov_model *cov, double *x, int lx, bool dist,
			    int *xdimOZ){ 
  // *xdimOZ to distinguish from the previous partial_loc_set definition
  location_type *loc = Loc(cov);
  int err;
  bool ygiven = !dist && loc->ly != 0;
    
  if ((err = partial_loc_set(loc, x, ygiven ? x : NULL, 
			     lx, ygiven ? lx : 0, dist, *xdimOZ, 
			     NULL, loc->grid, false)) 
      != NOERROR) XERR(err);
}



void partial_loc_set(cov_model *cov, double *x, int lx, bool dist, bool grid){
  location_type *loc = Loc(cov);
  int err;
  bool cartesian = isCartesian(cov->isoown);
  if (!cartesian && loc->ly==0) add_y_zero(loc);
  if ((err = partial_loc_set(loc, x, cartesian ? NULL : ZERO, 
			     lx, !cartesian, dist, loc->xdimOZ,
			     NULL, grid, false)) 
      != NOERROR) XERR(err);
}


void partial_loc_setOZ(cov_model *cov, double *x, double *y, 
		     int lx, int ly, bool dist, int *xdimOZ){
  // *xdimOZ to distinguish from the previous partial_loc_set definition
  location_type *loc = Loc(cov);
  int err;
  
  //  printf("partial_loc_set dist = %d %d \n", dist, loc->ly);
  
  if ((err = partial_loc_set(loc, x, y, lx, ly, dist, *xdimOZ, 
			     NULL, loc->grid, false)) 
      != NOERROR) XERR(err);
}


void partial_loc_setOZ(cov_model *cov, double *x, int lx, bool dist, int *xdimOZ){
  // *xdimOZ to distinguish from the previous partial_loc_set definition
  location_type *loc = Loc(cov);
  int err;
  bool cartesian = isCartesian(cov->isoown);
  if (!cartesian && loc->ly==0) add_y_zero(loc);
  
  //  printf("partial_loc_set dist = %d %d \n", dist, loc->ly);
   
  if ((err = partial_loc_set(loc, x, cartesian ? NULL : ZERO,  
			     lx, !cartesian, dist, *xdimOZ, 
			     NULL, loc->grid, false)) 
      != NOERROR) XERR(err);
  //PMI(cov);
}



void partial_loc_setXY(cov_model *cov, double *x, double *y, int lx, int ly) {
  location_type *loc = Loc(cov);
  int err;
  assert(y != NULL);
  if ((err = partial_loc_set(loc, x, y, lx, ly, false,
			     loc->xdimOZ, NULL, loc->grid, 
			     false)) 
      != NOERROR) XERR(err);
}


void partial_loc_setXY(cov_model *cov, double *x, double *y, int lx) {
  location_type *loc = Loc(cov);
  int err;

  //  PMI(cov);
  
  if ((err = partial_loc_set(loc, x, y, lx, y == NULL ? 0 : lx, false,
			     loc->xdimOZ, NULL, loc->grid, 
			     false)) 
      != NOERROR) XERR(err);
}


void partial_loc_null(cov_model *cov) {
  location_type *loc = Loc(cov);
  loc->lx = loc->ly = 0;
  loc->x = NULL;
  loc->y = NULL;
}


void InverseCovMatrix(cov_model *cov, double *v) {
  location_type *loc = Loc(cov);
  int vdimtot = loc->totalpoints * cov->vdim;  
  CovList[cov->nr].covariance(cov, v);
  invertMatrix(v, vdimtot);
}


//////////////////////////////////////////////////////////////////////
// Schnittstellen

#define STANDARDINTERN					\
  if (reg < 0 || reg > MODEL_MAX) XERR(ERRORREGISTER);	\
  if (currentNrCov==-1) InitModelList();		\
  cov_model *cov = KEY[reg];				\
  if (cov == NULL) { ERR("register not initialised") }	\
  cov_model *truecov = !isInterface(cov) ?		\
    cov : cov->key == NULL ? cov->sub[0] : cov->key
 
#define STANDARDINTERN_SEXP						\
  if (INTEGER(reg)[0] < 0 || INTEGER(reg)[0] > MODEL_MAX) XERR(ERRORREGISTER); \
  if (currentNrCov==-1) InitModelList();				\
  cov_model *cov = KEY[INTEGER(reg)[0]];			\
  if (cov == NULL) { ERR("register not initialised") }		\
  cov_model VARIABLE_IS_NOT_USED *truecov =  !isInterface(cov)	?	\
    cov : cov->key == NULL ? cov->sub[0] : cov->key

//  if (cov->pref[Nothing] == PREF_NONE) { PMI(cov); XERR(ERRORINVALIDMODEL) }

SEXP Delete_y(SEXP reg) {
  STANDARDINTERN_SEXP; 
  int d;
  location_type *loc = Loc(cov);
  if (loc->y != NULL) {
    if (loc->y != loc->x) free(loc->y);
    loc->y = NULL;
  }
  if (loc->ygr[0] != NULL) {
    if (loc->ygr[0] != loc->xgr[0]) free(loc->ygr[0]);
    for (d=0; d<MAXSIMUDIM; d++) loc->ygr[d] = NULL;
  }
  loc->ly = 0;
  return NULL;
}

SEXP CovLoc(SEXP reg, SEXP x, SEXP y, SEXP xdimOZ, SEXP lx,
	    SEXP result) {
  STANDARDINTERN_SEXP;

  partial_loc_setXY(cov, REAL(x), TYPEOF(y) == NILSXP ? NULL : REAL(y),
		  INTEGER(lx)[0]);
  CovList[truecov->nr].covariance(truecov, REAL(result));
  // PMI(cov, "covloc"); printf("xdimOZ %d \n", INTEGER(xdimOZ)[0]);
  partial_loc_null(cov);
  if (Loc(cov)->xdimOZ != INTEGER(xdimOZ)[0]) BUG;
  
  return NULL;
}

SEXP CovMatrixSelectedLoc(SEXP reg, SEXP x, SEXP dist, SEXP xdimOZ, SEXP lx,
			  SEXP selected, SEXP nsel, 
			  SEXP result) {
  STANDARDINTERN_SEXP;
  partial_loc_set_matrixOZ(cov, REAL(x), INTEGER(lx)[0], LOGICAL(dist)[0], 
		  INTEGER(xdimOZ));
  CovList[truecov->nr].selectedcovmatrix(truecov, INTEGER(selected), 
					 INTEGER(nsel)[0], REAL(result));
  partial_loc_null(cov);
  if (Loc(cov)->xdimOZ != INTEGER(xdimOZ)[0]) BUG;
  
  return NULL;
}


SEXP CovMatrixSelected(SEXP reg, SEXP selected, SEXP nsel,
		       SEXP result) {
  
  STANDARDINTERN_SEXP;
  CovList[truecov->nr].selectedcovmatrix(truecov, INTEGER(selected), 
					 INTEGER(nsel)[0], REAL(result));
  
  return NULL;
}

/*

         .C("setListElements",Reg, nteger(1), as.integer(1),
             PACKAGE="RandomFields");

*/      

void CovIntern(int reg, double *x, double *y, int lx, int ly, double *value) {
  STANDARDINTERN;
  partial_loc_setXY(cov, x, y, lx, ly);
  CovList[truecov->nr].covariance(truecov, value);
  partial_loc_null(cov);
}

void CovIntern(int reg, double *x, double *value) {
  STANDARDINTERN;
  location_type *loc = Loc(cov);
  bool cartesian = isCartesian(cov->isoown);
  if (!cartesian && loc->ly==0) add_y_zero(loc);

  partial_loc_setXY(cov, x, cartesian ? NULL : ZERO, 1, !cartesian);
  CovList[truecov->nr].covariance(truecov, value);
  partial_loc_null(cov);
}


SEXP CovMatrixLoc(SEXP reg, SEXP x, SEXP dist, SEXP xdimOZ,
		  SEXP lx, SEXP result) {
  STANDARDINTERN_SEXP;
  partial_loc_set_matrixOZ(cov, REAL(x), INTEGER(lx)[0], LOGICAL(dist)[0], 
			 INTEGER(xdimOZ)); //

  //  PMI(cov, "covmatrixloc");

  CovList[truecov->nr].covmatrix(truecov, REAL(result));
 
 //  printf("***************************************************\n");
  //  { int i,j,k, tot=Loc(cov)->totalpoints;
  //    printf("covmatrixloc %ld %ld %d\n", CovList[cov->nr].covmatrix, covmatrix_select, INTEGER(lx)[0]);
  //  for (k=i=0; i<tot*tot; i+=tot) {
  //   for (j=0; j<tot; j++) printf("%f ", value[k++]);
  //   printf("\n");  }}


  partial_loc_null(cov);
  
  return(NULL);
}


SEXP CovMatrixIntern(SEXP reg, SEXP x, SEXP dist, SEXP grid,
		     SEXP lx, SEXP result) {
  
  STANDARDINTERN_SEXP;
  //PMI(truecov);
  partial_loc_set_matrix(cov, REAL(x), INTEGER(lx)[0], LOGICAL(dist)[0],
			 LOGICAL(grid)[0]);
  //PMI(cov);  
  //PMI(truecov);
  //assert(cov->xdimprev == 2);
  CovList[truecov->nr].covmatrix(truecov, REAL(result));
  partial_loc_null(cov);
  return(NULL);
}


SEXP VariogramIntern(SEXP reg, SEXP x, SEXP lx, SEXP result) {
  
  STANDARDINTERN_SEXP;
  location_type *loc = Loc(cov);

  //   printf("vario intern dist=%d\n", false);
 
  partial_loc_setOZ(cov, REAL(x), INTEGER(lx)[0], false, &(loc->xdimOZ));
  CovList[truecov->nr].variogram(truecov, REAL(result));
  partial_loc_null(cov);
  
  //PMI(truecov->calling, -1); 
  //assert(false);
return(NULL);
}


void PseudovariogramIntern(int reg, double *x, double *y,
			   int lx, int ly, double *value) {
  STANDARDINTERN;
  location_type *loc = Loc(cov);
  partial_loc_setOZ(cov, x, y, lx, ly, false, &(loc->xdimOZ));
  CovList[truecov->nr].pseudovariogram(truecov, value);
  partial_loc_null(cov);
}

void PseudovariogramIntern(int reg, double *x, double *value) {
  STANDARDINTERN;
  location_type *loc = Loc(cov);
  bool cartesian = isCartesian(cov->isoown);
  if (!cartesian && loc->ly==0) add_y_zero(loc);
  partial_loc_setOZ(cov, x, cartesian ? NULL : ZERO, 
		    1, !cartesian, false, &(loc->xdimOZ));
  CovList[truecov->nr].pseudovariogram(truecov, value);
  partial_loc_null(cov);
}






// bool close = GLOBAL.general.vdim_close_together;

void Fctn(double VARIABLE_IS_NOT_USED *X, cov_model *cov, double *value){
  if (value==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  cov_model *sub =  cov->sub[0];
  assert(cov->key == NULL);
  int 
    vdim12 = cov->vdim2[0] * cov->vdim2[1];
 
  assert(cov->Spgs != NULL);						
  pgs_storage *pgs= cov->Spgs;						
  assert(pgs->x != NULL);						
  location_type *loc = Loc(cov);					
  assert(loc != NULL);							
  bool grid =loc->grid && !loc->Time && loc->caniso==NULL;		
  bool trafo = (loc->Time || loc->caniso != NULL) && !loc->distances;	
  bool ygiven = loc->y != NULL || loc->ygr[0] != NULL;			
  long cumgridlen[MAXMPPDIM +1],					
    err = NOERROR;							
  int d,								
    *gridlen=pgs->gridlen,						
    *end=pgs->end,							
    *endy =pgs->endy,				
    *start=pgs->start,							
    *startny =pgs->startny,			
    *nx=pgs->nx,							
    *ny = pgs->delta,				
    tot = loc->totalpoints,						
    tsxdim = cov->xdimown; /* weicht u.U. von tsdim bei dist=true ab */ 
  double *x = pgs->x,/* freed only if grid */				
    *xstart= pgs->xstart,						
    *inc=pgs->inc,							
    *y = pgs->supportmin,	/* freed only if grid */		
    *ystart =pgs->supportmax,			
    *yy = NULL,		/* set by Transform2NoGrid */			
    *xx = NULL,		/* dito			   */			
    *incy =pgs->supportcentre;
  if (grid) {								
    cumgridlen[0] = 1;							
    for (d=0; d<tsxdim; d++) {						
      inc[d] = loc->xgr[d][XSTEP];					
      gridlen[d] = loc->length[d];					
      cumgridlen[d+1] = gridlen[d] * cumgridlen[d];			
									
      start[d] = 0;							
      end[d] = gridlen[d];						
      nx[d] = start[d];							
      x[d] = xstart[d] = loc->xgr[d][XSTART];				
    }									
  }									
  loc->i_col = loc->i_row = 0;
	      

  //  assert(false);

#define FUNCTION FCTN(x, sub, value)
#define FUNCTION_Y NONSTATCOV(x, y, sub, value)

  if (grid) {
    if (ygiven) {							
      STANDARDSTART_Y_SUPPL;
      assert(loc->i_row==0);
      while (true) {
	FUNCTION_Y; 
	loc->i_row ++;
	value += vdim12;
	STANDARDINKREMENT_Y;
	RECYCLE_Y;
	STANDARDINKREMENT; 
      }
    } else {	// grid, y not given
      assert(loc->i_row==0);
      while (true) {
	FUNCTION;
	loc->i_row ++;;  
	value += vdim12;
	STANDARDINKREMENT; 
      } 
    }
  } else { // not a grid
    //   PMI(cov);
    if (trafo) {
      Transform2NoGrid(cov, &xx, &yy);   
      x = xx;
      y = yy;
    } else { x=loc->x; y=loc->y;  }

    if (ygiven) {
      double *y0, *yend;
      y0 = y;
      yend = y + tsxdim * loc->ly;

      for (; loc->i_row<tot; value+=vdim12, x+=tsxdim, y+=tsxdim, loc->i_row++){
	if (y >= yend) y = y0;
	FUNCTION_Y;
      }
    } else {
      for (; loc->i_row < tot; value+=vdim12, x+=tsxdim, loc->i_row++) {
	FUNCTION;
	//printf("x=%f %f %d\n", x[0], x[1], tsxdim, *value);
	//	assert(loc->i_row < 5);
      }
    }
  }

  STANDARD_ENDE;
  if (err != NOERROR) XERR(err); 
}
