/*
Authors
Martin Schlather, schlather@math.uni-mannheim.de

main library for unconditional simulation of random fields

 Copyright (C) 2005 -- 2017 Martin Schlather

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
  
 
#include <stdio.h>  
#include <string.h> 
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "variogramAndCo.h"

model* wheregenuineStatOwn(model *cov) {
  model *sub = cov;
  if (equalsnowGaussMethod(sub) || SUBNR==GAUSSPROC) {
    sub = sub->sub[0];
    while (equalsnowGaussMethod(sub) || SUBNR==GAUSSPROC) sub = sub->sub[0];
  } else if (isnowProcess(sub)) {
    NotProgrammedYet("");
  };

  if (cov->pref[Nothing] == PREF_NONE ||
      (!isnowPosDef(sub) && (!isnowVariogram(sub) || !isXonly(SYSOF(sub))))) {
    // Variogramme sind hier definitiv erlaubt, da durch Addition einer 
    // Konstanten, das Zeug zu einer Kovarianzmatrix gemacht werden kann
    // siehe direct.cc
    //assert(({PMI(cov, "cov matrix"); true;})); //
    ERR("covariance matrix: given model is not a covariance function");
  }
  
  return sub;
} 


#define STANDARDSTART							\
  model *cov = Cov;					\
  assert(cov != NULL); 				\
  if (equalsnowGaussMethod(cov) || COVNR==GAUSSPROC) cov = cov->sub[0];	\
  model *calling = cov;				\
  if (calling->Spgs == NULL) { 		\
    assert(isnowVariogram(calling)); 		\
    calling = cov->calling; /* either interface or process	*/	\
    if (calling != NULL && calling->Spgs == NULL) {			\
      assert(isnowProcess(calling));					\
      calling = calling->calling;					\
    }									\
  }									\
  FINISH_START(calling, cov, false, 0); 

//  printf("OKcc\n\n");PMI0(calling);PMI0(cov);assert(calling != NULL && (equalsnowInterface(calling) || isnowProcess(calling))); printf("OK\n\n");


/* assert(VDIM0 == VDIM1); */	
			

void CovVario(model *Cov, bool is_cov, bool pseudo, double *value) {
  // does not return a matrix, just a vector of values !
  // therefor loc->distances do not make sense
  // STANDARDSTART; //printf("xxxxxx\n");




 model *cov = Cov;					
  assert(cov != NULL); 				
  if (equalsnowGaussMethod(cov) || COVNR==GAUSSPROC) cov = cov->sub[0];	
  model *calling = cov;				
  if (calling->Spgs == NULL) { 		
    assert(isnowVariogram(calling)); 		
    calling = cov->calling; /* either interface or process	*/	
    if (calling != NULL && calling->Spgs == NULL) {			
      assert(isnowProcess(calling));					
      calling = calling->calling;					
    }									
  }
  assert(calling != NULL);

  // PMI(calling); PMI(cov);
  
  if ((calling)->Spgs == NULL) { BUG;	  }

  
  pgs_storage *pgs= (calling)->Spgs;					
  assert(pgs->x != NULL);						
  location_type *loc = Loc(calling);					
  assert(loc != NULL);							
  /*bool grid =loc->grid && !loc->Time && loc->caniso==NULL, why !time ?? 10.3.15*/		
  bool grid =loc->grid && loc->caniso==NULL,				
    trafo = (loc->Time || loc->caniso != NULL) && !loc->distances;	
  long cumgridlen[MAXMPPDIM +1],					
    tot = loc->totalpoints;						
  int d,								
    *gridlen=pgs->gridlen,						
    *start=pgs->start,							
    *end=pgs->end,							
    i_row = 0,								
    *nx=pgs->nx,							
    tsxdim = PREVTOTALXDIM; /* weicht u.U. von logicaldim bei dist=true ab */ 
  double *x = pgs->x,							
    *xstart= pgs->xstart,						
    *inc=pgs->inc,							
    *xx = NULL;		/* dito			   */			
   if (grid) {								
    STANDARDSTART_X;							
  }   								
  /*bool grid =loc->grid && !loc->Time && loc->caniso==NULL, why !time ?? 10.3.15*/		
  bool ygiven = !(false) && (loc->y != NULL || loc->ygr[0] != NULL);/* might be changed !*/ 
  int i_col = 0,							
    vdim0  = (cov)->vdim[0],					
    vdimSq = vdim0 * (cov)->vdim[0], /*not nec. squ!!*/	
    *endy =pgs->endy,							
    *startny =pgs->startny,						
    *ny = pgs->delta;	 						
  long VARIABLE_IS_NOT_USED vdim0tot = vdim0 * tot,			
    vdimSqtot = vdimSq * tot,						
    err = NOERROR;							
  double *ystart =pgs->supportmax,					
     *yy = NULL,		/* set by Transform Loc; freed  */	
    *incy =pgs->supportcentre					;
  



    
  
  long m, n;
  bool stat;
  assert(pgs->supportmin != NULL);
  double 
    *zero = ZERO(cov),
    *y = pgs->supportmin, // vgl. def von *y bei Matrizen
    *C0x = pgs->C0x,
    *C0y = pgs->C0y;
  bool kernel = !isXonly(PREV);
  INCLUDE_VAL;
 
  if (loc->distances) BUG;

  model *genuine = wheregenuineStatOwn(cov);
  bool isvario = isnowVariogram(genuine);
  stat = !kernel && isvario;

  if (is_cov) {
    assert(({/*P M I(cov, "Cov/vario"); */ cov->pref[Nothing] != PREF_NONE &&
	    isnowShape(genuine);})); //
  } else {
    if (cov->pref[Nothing] == PREF_NONE || !isvario) {
      assert(({PMI(cov); true;})); //
      ERR("given model is not a variogram");
    }
    if (stat && isCartesian(PREV)) {   
      if (!isCartesian(OWN)) BUG;
      COV(zero, genuine, C0y);      
    } else {
      if (vdim0 > 1 && !stat) 
	ERR("multivariate variogram only calculable for stationary models");
      NONSTATCOV(zero, zero, genuine, C0y);
    } 
  }
   

#define UNIVAR COV(x, genuine, value)
#define UNIVAR_Y NONSTATCOV(x, y, genuine, value) 

#define MULT 					\
  COV(x, genuine, cross);				\
  for (v = 0; v<vdimSq; v++) Val[v][i_row] = cross[v];

#define MULT_Y				\
  NONSTATCOV(x, y, genuine, cross);	\
  for (v = 0; v<vdimSq; v++) Val[v][i_row] = cross[v];	
 
#define VARIO_UNIVAR	  \
  COV(x, genuine, cross);	  \
  *value = *C0y - *cross; 

#define VARIO_UNIVAR_Y	   	\
  NONSTATCOV(x, y, genuine, cross);		\
  NONSTATCOV(x, x, genuine, C0x);		\
  NONSTATCOV(y, y, genuine, C0y);		\
    *value = 0.5 * (*C0x + *C0y) - *cross;
    
  
#define VARIO_MULT 							\
  COV(x, genuine, cross);						\
  for  (v=0, m=0; m<vdim0; m++) {					\
    for (n=0; n<vdim0; n++, v++) {					\
      Val[v][i_row] =						\
	/* 0.5 * (C0x[v] + C0y[v] - cross[m*vdim0+n] -cross[n*vdim0+m]); */ \
	C0y[v] - 0.5 * (cross[v] + cross[n*vdim0+m]);			\
    }									\
  }				 					

   
#define VARIO_MULT_Y 							\
  NONSTATCOV(x, y, genuine, cross);					\
  NONSTATCOV(x, x, genuine, C0x);					\
  NONSTATCOV(y, y, genuine, C0y);					\
  for  (v=m=0; m<vdim0; m++) {						\
    for (n=0; n<vdim0; n++, v++) {					\
	Val[v][i_row] =						\
	  0.5 * (C0x[v] + C0y[v] - cross[v] -cross[n*vdim0+m]);		\
    }									\
  }								
 
#define PSEUDO_MULT 							\
  COV(x, genuine, cross);						\
  for (v=0; v<vdimSq; v++) Val[v][i_row] = C0y[v] - cross[v];

#define PSEUDO_MULT_Y							\
  NONSTATCOV(x, y, genuine, cross);					\
  NONSTATCOV(x, x, genuine, C0x);					\
  NONSTATCOV(y, y, genuine, C0y);					\
  for (v=0; v<vdimSq; v++) {						\
    Val[v][i_row] = 0.5 * (C0x[v] + C0y[v]) - cross[v];		\
  }								

  // printf("covvario 2 %.50s, line %d %d %ld %ld %ld\n",				
  //     __FILE__, __LINE__, ygiven, y, zero, pgs->supportmin) ;	

  //  printf("grid = %d ygiven=%d kernel=%d trafo=%d\n", grid, ygiven, kernel, trafo);
  
  if (is_cov) {  
    // printf("XX\n"); 
    PERFORM(UNIVAR, MULT, UNIVAR_Y, MULT_Y);
  } else if (pseudo) {
    PERFORM(VARIO_UNIVAR, PSEUDO_MULT, VARIO_UNIVAR_Y, PSEUDO_MULT_Y);
  } else {
    PERFORM(VARIO_UNIVAR, VARIO_MULT, VARIO_UNIVAR_Y, VARIO_MULT_Y);
    /*
if (grid) {					
  if (ygiven || kernel) {  		 
      STANDARDSTART_Y_SUPPL;	 		
      if (vdim0 == 1) { GRIDCYCLES(VARIO_UNIVAR_Y; value+=vdimSq);	
      } else { GRIDCYCLES(VARIO_MULT_Y); }				
    } else { // grid, y not given					
      if (vdim0 == 1) {							
	GRIDCYCLE_X(VARIO_UNIVAR; value+=vdimSq);				
      } else {GRIDCYCLE_X(VARIO_MULT); }				
    }									
} else { // not a grid 						
    if (trafo) {							
      TransformLoc(cov, &xx, &yy, false);				
      x = xx;								
      if (ygiven) y = yy;						
    } else {								
      x=loc->x;								
      if (ygiven) y=loc->y;						
    }	  								
    assert(ygiven xor (y==zero));					
    if (ygiven || kernel) {						
      double *y0, *yend;						
      yend = ygiven ? y + tsxdim * loc->ly : zero;			
      y0 = y;								
      NONGRIDCYCLE(DO_INCREMENTY, DO_RECYCLY, VARIO_UNIVAR_Y; value+=vdimSq,
		   VARIO_MULT_Y);					
    } else { 
      //NONGRIDCYCLE(EMPTY, EMPTY, VARIO_UNIVAR; value+=vdimSq, VARIO_MULT); } 
      //#define NONGRIDCYCLE(INCREMENTY, RECYCLY, FCTN1, FCTN2)	

      if (vdim0 == 1) { 							
	for (; i_row<tot; i_row++, x+=tsxdim ){
	  p rintf("i=%d % ld %d\n", i_row, value, tot); 
	  COV(x, genuine, cross);		       
	  *value = *C0y - *cross;	
	  value+=vdimSq;
	}									
      } else {								
	for (; i_row<tot; i_row++, x+=tsxdim ){		
	  VARIO_MULT;							
	}									
      }
    }

    
    STANDARD_ENDE;							
    if (err != NOERROR) XERR(err);					
  }
*/


  }
} 
 
 
#define swap(x, y) { double swapdummy = x;  x = y;  y = swapdummy; }
#define swapInt(x, y) { int swapdummy = x;  x = y;  y = swapdummy; }

void CovarianceMatrix(model *Cov, double *v) {
  STANDARDSTART;
  
  model * genuine = wheregenuineStatOwn(cov);

  bool dist = loc->distances,
    vdim_closetogether = GLOBAL.general.vdim_close_together;
  long l,n,m, VDIM, NEND, NINCR, MINCR, ENDFORINCR,
    totM1 = tot - 1,							
    vdim0totSq = vdim0tot * tot,
    vdimSqtotSq = vdimSq * tot * tot;				       
  double *C = NULL,
    *x0 = NULL,	/* never free it  */				
    *z = pgs->z;							
 
  double *y = pgs->supportmin;
  if (grid) {
    STANDARDSTART_Y_SUPPL;
  }

  if (ygiven && (loc->x != loc->y || loc->xgr[0] != loc->ygr[0])) {
    // ein Paerchen ist NULL;
    GERR("for the covariance matrix, no y-value may be given");
  }
   
  if (vdim_closetogether) {
    // v-dimensions close together
    VDIM = vdim0;
    NEND = vdimSqtot;
    NINCR = vdim0tot;
    ENDFORINCR = vdim0;
    MINCR = 1;
  } else {
    // values of any single multivariate component close together
    // default in GLOBAL.CovMatrixMulti
    VDIM = 1;
    NEND = vdimSqtotSq;
    NINCR = vdim0totSq;
    ENDFORINCR = vdim0tot;
    MINCR = tot;
  }
  
#define MULTICOV						\
  C = v + VDIM * (i_col + i_row * vdim0tot);			\
  for (l=n=0; n<NEND; n+=NINCR) {					\
    long endfor = n + ENDFORINCR;					\
    for (m=n; m<endfor; m+=MINCR) {					\
      C[m] = z[l++];							\
    }									\
  }									\
									\
  if (i_col != i_row) {					\
    C = v + VDIM * (i_row + i_col * vdim0tot);		\
    for (l=m=0; m<ENDFORINCR; m+=MINCR) {				\
      for (n=m; n<NEND; n+=NINCR) {					\
	C[n] = z[l++];							\
      }									\
    }									\
  }

  
  // printf("enterign %.50s\n", NAME(cov)); if (zaehlerx==1) crash();
  // zaehlerx++;
  
  if (grid) {
    while (true) {
      i_col = i_row;
      for (d=0; d<tsxdim; d++) {
	y[d] = x[d];
	ystart[d] = xstart[d];
	ny[d] = nx[d];
      }
      y[tsxdim] = i_col; 
      while (true) {
	//printf("irow = %d %d dim=%d total=%d nx=%d %d, %d %d\n", i_col, i_row, tsxdim, Gettotalpoints(cov), nx[0], ny[0], end[0], endy[0]);
	NONSTATCOV(x, y, genuine, z);
	//	APMI(genuine); 
	//	printf("nostat %d %d x=%10g %10g y=%10g %10g z=%10g \n", i_row, i_col, *x, x[1], *y, y[1], *z);

	MULTICOV; 
	STANDARDINKREMENT_Y;
	//	APMI(genuine);
	if (d >= tsxdim) break; 
      }
      STANDARDINKREMENT_X;
    }
  } else { // not a grid
    int localdim = tsxdim; // just to control
    if (trafo) {
      localdim = TransformLoc(cov, &xx, false);    
      assert(localdim == tsxdim);

      x0 = xx;
    } else x0 = loc->x;
    x = x0;

    // i_row/col fkt nicht mit parallel!!!!
    // am saubersten waere es, das i als zusaetzliche Koordinate
    // in einem extra koordinatensystem zu uebertragen
    // auch muss nugget zu non-stat nugget umgeschrieben werden,
    // so dass nur fuer nicht-stat. cov-fkten die weitere Koordinate
    // uebertragen werden muss + zusaetzlich die Abfrage "stationaer"
    // ausreicht, um loc->i* auszuschliessen. Oder aber ein Flag
    // setzen in model, das anzeigt, ob loc->i in einem untermodel
    // vorliegt Auch nugget (messfehler) ist eigentlich nicht-stationaer!
    //#ifdef DO_PARALLEL
    //#pragma omp parallel for num_threads(CORES) if (tot > 20) schedule(dynamic, 10)
    //#endif
    if (dist) {
      double zero3[3] = {0.0};
      for (i_row=0; i_row<tot; i_row++) {
	x = x0 + localdim * i_row;
	for (y=x, i_col=i_row; i_col<tot; i_col++, y+=localdim) {
	  if (i_col==i_row) {
	    zero3[1] = zero3[2] = (double) i_row;
	    COV(zero3, genuine, z); 
	  } else {
	    double x3[3] =  { x0[(i_row * totM1 - (i_row * (i_row + 1)) / 2
				  + i_col -1) * localdim],
			      (double) i_row,
			      (double) i_col };
	    COV(x3, genuine, z);
	  }
	  MULTICOV; 
	}
      }
    } else { // no distances
      int xMem = (localdim + 1)* sizeof(double),
	 xCpy = localdim * sizeof(double);
      double *X = (double*) MALLOC(xMem),
	*Y = (double*) MALLOC(xMem);
      for (i_row=0; i_row<tot; i_row++) {
	x = x0 + localdim * i_row;
	MEMCOPY(X, x, xCpy);
	X[localdim] = i_row;
	
	for (y=x, i_col=i_row; i_col<tot; i_col++, y+=localdim) {
	  MEMCOPY(Y, y, xCpy);
	  Y[localdim] = i_col;
	  NONSTATCOV(X, Y, genuine, z);
	  if (!R_FINITE(z[0])) GERR("model creates non-finte values.");	
	  MULTICOV;
	}
      }
      FREE(X);
      FREE(Y);
    } // no distances
  } // not a grid
  

  if (false) {
    for (m=0; m<27; m++) {
      for (n=0; n<27; n++) {
	//printf("%+2.2f ", v[n * 27 + m]);
      }
      //printf("\n"); 
    }
  }
  //assert(false);   FREE(cov->Spgs->z); printf("zzz\n");

 ErrorHandling: 
  STANDARD_ENDE;
  if (err!=NOERROR) XERR(err); 
   
} // CovarianzMatrix





void CovarianceMatrixCol(model *Cov, int column, double *v) {
  STANDARDSTART;
  
  model * genuine = wheregenuineStatOwn(cov);      
  bool dist = loc->distances,
    vdim_closetogether = GLOBAL.general.vdim_close_together;
  long l,n,m, VDIM, NEND, NINCR, MINCR, ENDFORINCR,
    totM1 = tot - 1;
    //    vdim0totSq = vdim0tot * tot,
  //vdimSqtotSq = vdimSq * tot * tot; 
  double *C = NULL,
    *x0 = NULL,	/* never free it  */				
    *z = pgs->z;							
  
  double *y = pgs->supportmin;
  if (grid) {
    STANDARDSTART_Y_SUPPL;
  }

  if (ygiven && (loc->x != loc->y || loc->xgr[0] != loc->ygr[0])) {
    // ein Paerchen ist NULL;
    GERR("for the covariance matrix, no y-value may be given");
  }
   
  if (vdim_closetogether) {
    // v-dimensions close together
    VDIM = vdim0;
    NEND = vdimSqtot;
    NINCR = vdim0tot;
    ENDFORINCR = vdim0;
    MINCR = 1;
  } else {
    // values of any single multivariate component close together
    // default in GLOBAL.CovMatrixMulti
    VDIM = 1;
    NEND = vdimSqtot; // vdimSqtotSq;
    NINCR = vdim0tot; // vdim0totSq;
    ENDFORINCR = vdim0tot;
    MINCR = tot;
  }

  
#define MULTICOV_COL							\
  C = v + VDIM * i_row;						\
  for (l=m=0; m<ENDFORINCR; m+=MINCR) {					\
    for (n=m; n<NEND; n+=NINCR) {					\
      C[n] = z[l++];							\
    }									\
 }

  // symmetrischer Teil gegenueber CovMatrix fehlt natuerlich.
  
  
  if (grid) {
    BUG;
    i_col = 0;
    while (true) {
      i_row = 0; /// i_col;
      for (d=0; d<tsxdim; d++) {
	y[d] = x[d];
	ystart[d] = xstart[d];
	ny[d] = nx[d];
      }
      while (true) {
	NONSTATCOV(x, y, genuine, z);
	//	APMI(genuine); 
	//	printf("nostat %d %d x=%10g %10g y=%10g %10g z=%10g \n", i_row, i_col, *x, x[1], *y, y[1], *z);

	MULTICOV_COL; 
	STANDARDINKREMENT_Y;
	//	APMI(genuine);
	if (d >= tsxdim) break; 
      }
      STANDARDINKREMENT_X;
    }
  } else { 
    int localdim = tsxdim; // just to control
    if (trafo) {
      localdim = TransformLoc(cov, &xx, false);    
      assert(localdim == tsxdim);

      x0 = xx;
    } else x0 = loc->x;
    x = x0;

    i_row = column;
    if (dist) {
      double zero3[3] = {0.0};
      for (y=x0, i_col=0; i_col<tot;  i_col++, y+=localdim) {
	if (i_row==i_col) {
	  zero3[1] = zero3[2] = (double) i_row;
	  COV(zero3, genuine, z); 
	} else {
	  int i_min = MIN(i_col, i_row),
	    i_max = MAX(i_col, i_row);
	  double x3[3] =  { x0[(i_min * totM1 - (i_min * (i_min + 1)) / 2
				  + i_max -1) * localdim],
			    (double) i_row,
			      (double) i_col };	   
	  COV(x3, genuine, z);
	}
	
	MULTICOV_COL;
      }
    } else {// no distances
      int xBytes = (localdim + 1) * sizeof(double);
      double *X = (double*) MALLOC(xBytes),
	*Y = (double*) MALLOC(xBytes);
      x += localdim * column;
      MEMCOPY(X, x, xBytes);
      X[localdim] = i_row;
      for (y=x0, i_col=0; i_col<tot;  i_col++, y+=localdim) {
	MEMCOPY(Y, y, xBytes);
	Y[localdim] = i_col;
	NONSTATCOV(x, y, genuine, z);
	assert(R_FINITE(z[0]));
  	MULTICOV_COL;
      }
      FREE(X);
      FREE(Y);
    } // no distances
  } // not a grid
  

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
  //    for (j=0; j<tot;j++) printf("%10g ", v[k++]);
  //    printf("\n");  }
  
} // CovarianzMatrixCol



void partial_loc_set_matrix(model *cov, double *x, long lx, bool dist,
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

void partial_loc_set_matrixOZ(model *cov, double *x, long lx, bool dist,
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



void partial_loc_set(model *cov, double *x, long lx, bool dist, bool grid){
  location_type *loc = Loc(cov);
  int err;
  //  bool cartesian = isCartesian(OWNISO(0));
  //  if (!cartesian && loc->ly==0) add_y_zero(loc);
  if ((err = partial_loc_set(loc, x, NULL, // cartesian ? NULL : ZERO(), 
			     lx, 0, //!cartesian,
			     dist, loc->xdimOZ,
			     NULL, grid, false)) 
      != NOERROR) XERR(err);
}


void partial_loc_setOZ(model *cov, double *x, double *y, 
		       long lx, long ly, bool dist, int *xdimOZ){
  // *xdimOZ to distinguish from the previous partial_loc_set definition
  location_type *loc = Loc(cov);
  int err;
  
  //  printf("partial_loc_set dist = %d %d \n", dist, loc->ly);
  
  if ((err = partial_loc_set(loc, x, y, lx, ly, dist, *xdimOZ, 
			     NULL, loc->grid, false)) 
      != NOERROR) XERR(err);
}


void partial_loc_setOZ(model *cov, double *x, long lx, bool dist, int *xdimOZ){
  // *xdimOZ to distinguish from the previous partial_loc_set definition
  location_type *loc = Loc(cov);
  int err;
  //  bool cartesian = isCartesian(OWNISO(0));
  // if (!cartesian && loc->ly==0) add_y_zero(loc);
  
  //  printf("partial_loc_set dist = %d %d \n", dist, loc->ly);
   
  if ((err = partial_loc_set(loc, x, NULL, // cartesian ? NULL : ZERO(),  
			     lx, 0, //!cartesian,
			     dist, *xdimOZ, 
			     NULL, loc->grid, false)) 
      != NOERROR) XERR(err);
  //PMI(cov);
}



void partial_loc_setXY(model *cov, double *x, double *y, long lx, long ly) {
  location_type *loc = Loc(cov);
  int err;

  // if (y == NULL)  crash();
  // assert(y != NULL);
 
  if ((err = partial_loc_set(loc, x, y, lx, ly, false,
			     loc->xdimOZ, NULL, loc->grid, 
			     false)) 
      != NOERROR) XERR(err);
}


void partial_loc_setXY(model *cov, double *x, double *y, long lx) {
  location_type *loc = Loc(cov);
  int err;

  //  PMI(cov);
  
  if ((err = partial_loc_set(loc, x, y, lx, y == NULL ? 0 : lx, false,
			     loc->xdimOZ, NULL, loc->grid, 
			     false)) 
      != NOERROR) XERR(err);
}


void partial_loc_null(model *cov) {
  location_type *loc = Loc(cov);
  loc->lx = loc->ly = 0;
  loc->x = NULL;
  loc->y = NULL;
}


void InverseCovMatrix(model *Cov, double *v, double *det) {
  model *cov = Cov;
  model *genuine = wheregenuineStatOwn(cov);
  location_type *loc = Loc(cov);
  long vdimtot = loc->totalpoints * VDIM0;
  KEY_type *KT = cov->base;    
  assert(VDIM0 == VDIM1);
  DefList[COVNR].covariance(genuine, v);
  if (cov->Ssolve == NULL) SOLVE_STORAGE;
  Ext_setErrorLoc(KT->error_loc);
  //  printf("inverse\n");
  int Exterr = Ext_solvePosDef(v, vdimtot, true, NULL, 0, det, cov->Ssolve);
  if (Exterr != NOERROR){
    Ext_getErrorString(cov->err_msg);
    OnErrorStop(Exterr, cov->err_msg);
  }
}


//////////////////////////////////////////////////////////////////////
// Schnittstellen

#define STANDARDINTERN					\
  if (reg < 0 || reg > MODEL_MAX) XERR(ERRORREGISTER);	\
   assert(currentNrCov != UNSET);		\
   model *cov = KEY()[reg];				\
  if (cov == NULL) { ERR("register not initialised") }	\
  model *truecov = !equalsnowInterface(cov) ?		\
    cov : cov->key == NULL ? cov->sub[0] : cov->key
 
#define STANDARDINTERN_SEXP_BASIC					\
  if (INTEGER(reg)[0] < 0 || INTEGER(reg)[0] > MODEL_MAX) XERR(ERRORREGISTER); \
   assert(currentNrCov != UNSET);				\
   model *cov = KEY()[INTEGER(reg)[0]];				\
  if (cov == NULL) { ERR("register not initialised") }			

#define STANDARDINTERN_SEXP						\
  STANDARDINTERN_SEXP_BASIC;						\
  model VARIABLE_IS_NOT_USED *truecov =  !equalsnowInterface(cov)	? \
    cov : cov->key == NULL ? cov->sub[0] : cov->key;			\
  if (equalsnowGaussMethod(truecov)) truecov = truecov->sub[0]

//  if (cov->pref[Nothing] == PREF_NONE) { PMI(cov); XERR(ERRORINVALIDMODEL) }

SEXP CovLoc(SEXP reg, SEXP x, SEXP y, SEXP xdimOZ, SEXP lx,
	    SEXP result) {
  STANDARDINTERN_SEXP;

  //   PMI(cov);
  if (Loc(cov)->len > 1) BUG;  

  partial_loc_setXY(cov, REAL(x), TYPEOF(y) == NILSXP ? NULL : REAL(y),
		    INTEGER(lx)[0]);
  DefList[MODELNR(truecov)].covariance(truecov, REAL(result));
  //   pmi(cov,0); printf("%d xdimOZ %d \n", Loc(cov)->xdimOZ, INTEGER(xdimOZ)[0]);
  partial_loc_null(cov);
  if (Loc(cov)->xdimOZ != INTEGER(xdimOZ)[0]) BUG;
  return R_NilValue;
}

void CovIntern(int reg, double *x, double *y, long lx, long ly, double *value) {
  STANDARDINTERN;
  partial_loc_setXY(cov, x, y, lx, ly);
  DefList[MODELNR(truecov)].covariance(truecov, value);
  partial_loc_null(cov);
}

 
SEXP CovMatrixLoc(SEXP reg, SEXP x, SEXP dist, SEXP xdimOZ, 
		  SEXP lx, SEXP result) {
  STANDARDINTERN_SEXP;
  partial_loc_set_matrixOZ(cov, REAL(x), INTEGER(lx)[0], LOGICAL(dist)[0], 
			 INTEGER(xdimOZ)); //

  //  PMI(cov, "covmatrixloc");

  DefList[MODELNR(truecov)].covmatrix(truecov, REAL(result));
 
 //  printf("***************************************************\n");
  //  { int i,j,k, tot=Loc(cov)->totalpoints;
  //    printf("covmatrixloc %ld %ld %d\n", DefList[COVNR].covmatrix, covmatrix_select, INTEGER(lx)[0]);
  //  for (k=i=0; i<tot*tot; i+=tot) {
  //   for (j=0; j<tot; j++) printf("%10g ", value[k++]);
  //   printf("\n");  }}


  partial_loc_null(cov);
  
  return(NULL);
}


SEXP CovMatrixIntern(SEXP reg, SEXP x, SEXP dist, SEXP grid,
		     SEXP lx, SEXP result) {
  
  STANDARDINTERN_SEXP;
  partial_loc_set_matrix(cov, REAL(x), INTEGER(lx)[0], LOGICAL(dist)[0],
			 LOGICAL(grid)[0]);
  DefList[MODELNR(truecov)].covmatrix(truecov, REAL(result));
  partial_loc_null(cov);
  return(NULL);
}


SEXP VariogramIntern(SEXP reg) {
  STANDARDINTERN_SEXP;
  //  location_type *loc = Loc(cov);
  SEXP ans;
  int 
    vdim =VDIM0,
    size = Gettotalpoints(cov) * vdim * vdim;
  PROTECT(ans = allocVector(REALSXP, size));
  DefList[MODELNR(truecov)].variogram(truecov, REAL(ans)); 
  UNPROTECT(1);
  return(ans);
}


SEXP PseudovariogramIntern(SEXP reg) {
  STANDARDINTERN_SEXP;
  //  location_type *loc = Loc(cov);
  SEXP ans;
  int 
    vdim =VDIM0,
    size = Gettotalpoints(cov) * vdim * vdim;
  PROTECT(ans = allocVector(REALSXP, size));
  DefList[MODELNR(truecov)].pseudovariogram(truecov, REAL(ans)); 
  UNPROTECT(1);
  return(ans);
}

/*
void PseudovariogramIntern(int reg, double *x, double *y,
			   long lx, long ly, double *value) {
  STANDARDINTERN;
  location_type *loc = Loc(cov);
  partial_loc_setOZ(cov, x, y, lx, ly, false, &(loc->xdimOZ));
  DefList[MODELNR(truecov)].pseudovariogram(truecov, value);
  partial_loc_null(cov);
}

void PseudovariogramIntern(int reg, double *x, double *value) {
  STANDARDINTERN;
  location_type *loc = Loc(cov);
  // bool cartesian = isCartesian(OWNISO(0));
  // if (!cartesian && loc->ly==0) add_y_zero(loc);
  partial_loc_setOZ(cov, x, NULL, // cartesian ? NULL : ZERO(), 
		    1, 0, // !cartesian,
		    false, &(loc->xdimOZ));
  DefList[MODELNR(truecov)].pseudovariogram(truecov, value);
  partial_loc_null(cov);
}

*/
