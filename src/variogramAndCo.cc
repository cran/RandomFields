/*
Authors
Martin Schlather, schlather@math.uni-mannheim.de

main library for unconditional simulation of random fields

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2015 Martin Schlather
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
#include "Operator.h"
#include "variogramAndCo.h"
#include "xport.h"

void genuineStatOwn(cov_model *cov, domain_type *stat, Types *type) {
  cov_model *sub = cov;
  while (isGaussMethod(sub)) sub = sub->sub[0];
  *stat = sub->domown;
  *type = sub->typus;
}

#define STANDARDSTART							\
  assert(cov != NULL);							\
  cov_model *prev = cov->calling;					\
  assert(prev != NULL);							\
  bool ok_ = prev->nr == COVFCTN || prev->nr == VARIOGRAM_CALL ||	\
    prev->nr == PREDICT_CALL ||						\
    CovList[prev->nr].check == check_cov || prev->nr == COVMATRIX ||	\
    CovList[prev->nr].check == check_fctn || prev->nr == DIRECT ||	\
    (prev->nr == GAUSSPROC && prev->calling != NULL &&			\
     CovList[prev->calling->nr].check == check_likelihood) ||		\
    CovList[prev->nr].check == checkS;					\
  if (!ok_) BUG;							\
  FINISH_START(prev, cov, false, 0);					\


/* assert(cov->vdim[0] == cov->vdim[1]); */	
			
//  PMI(prev); print("\n%d %d %d\n", 	prev->nr, PREDICT_CALL, ok_);  


void CovVario(cov_model *cov, bool is_cov, bool pseudo, double *value) { 
   STANDARDSTART;
 domain_type domown;
  Types type;
  long m, n;
  bool stat;
  double 
    *y = ygiven ? pgs->supportmin : ZERO, // vgl. def von *y bei Matrizen
    *C0x = pgs->C0x,
    *C0y = pgs->C0y;
  bool kernel = cov->domprev != XONLY;  
  INCLUDE_VAL;		

  if (loc->distances) BUG;

  // if(isCartesian(prev->isoown) || ygiven)
  // printf("covvario %s, line %d %d %ld %ld %ld\n", __FILE__, __LINE__, ygiven, y, ZERO, pgs->supportmin) ;	

  genuineStatOwn(cov, &domown, &type);
  stat = !kernel && isVariogram(type);

  if (is_cov) {
    assert(({/*PMI(cov, "Cov/vario"); */ cov->pref[Nothing] != PREF_NONE && isShape(type);})); //
  } else {
    if (cov->pref[Nothing] == PREF_NONE || !isVariogram(type)) {
      assert(({PMI(cov); true;})); //
      ERR("given model is not a variogram");
    }
    if (stat && isCartesian(prev->isoown)) {   
      //  PMI(cov);     printf("%f\n", ZERO[0]);

      COV(ZERO, cov, C0y);      
    } else {
      if (vdim0 > 1 && !stat) 
	ERR("multivariate variogram only possible for stationary models");
      NONSTATCOV(ZERO, ZERO, cov, C0y);
    } 
  }
   

#define UNIVAR COV(x, cov, value)
#define UNIVAR_Y NONSTATCOV(x, y, cov, value) 

#define MULT 					\
  COV(x, cov, cross);				\
  for (v = 0; v<vdimSq; v++) Val[v][loc->i_row] = cross[v];

#define MULT_Y					\
  NONSTATCOV(x, y, cov, cross);	\
  for (v = 0; v<vdimSq; v++) Val[v][loc->i_row] = cross[v];	
 
#define VARIO_UNIVAR 			\
  COV(x, cov, cross);			\
  *value = *C0y - *cross;	

#define VARIO_UNIVAR_Y							\
    NONSTATCOV(x, y, cov, cross);					\
    NONSTATCOV(x, x, cov, C0x);						\
    NONSTATCOV(y, y, cov, C0y);						\
    *value = 0.5 * (*C0x + *C0y) - *cross;			       

  
#define VARIO_MULT 							\
  COV(x, cov, cross);							\
  for  (v=0, m=0; m<vdim0; m++) {					\
    for (n=0; n<vdim0; n++, v++) {					\
      Val[v][loc->i_row] =						\
	/* 0.5 * (C0x[v] + C0y[v] - cross[m*vdim0+n] -cross[n*vdim0+m]); */ \
	C0y[v] - 0.5 * (cross[v] + cross[n*vdim0+m]);			\
    }									\
  }									

  //printf("\n\nygiven %d %d\n\n\n", ygiven, kernel);APMI(cov);
  
#define VARIO_MULT_Y 							\
  NONSTATCOV(x, y, cov, cross);						\
  NONSTATCOV(x, x, cov, C0x);						\
  NONSTATCOV(y, y, cov, C0y);						\
  for  (v=m=0; m<vdim0; m++) {						\
    for (n=0; n<vdim0; n++, v++) {					\
	Val[v][loc->i_row] =						\
	  0.5 * (C0x[v] + C0y[v] - cross[v] -cross[n*vdim0+m]);		\
    }									\
  }								
 
#define PSEUDO_MULT 							\
     COV(x, cov, cross);						\
     for (v=0; v<vdimSq; v++) Val[v][loc->i_row] = C0y[v] - cross[v];

#define PSEUDO_MULT_Y							\
  NONSTATCOV(x, y, cov, cross);						\
  NONSTATCOV(x, x, cov, C0x);						\
  NONSTATCOV(y, y, cov, C0y);						\
  for (v=0; v<vdimSq; v++) {						\
    Val[v][loc->i_row] = 0.5 * (C0x[v] + C0y[v]) - cross[v];		\
  }								

  // printf("covvario 2 %s, line %d %d %ld %ld %ld\n",				
  //     __FILE__, __LINE__, ygiven, y, ZERO, pgs->supportmin) ;	
  
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
      if (vdim0 == 1) { GRIDCYCLE_Y(VARIO_UNIVAR_Y; value+=vdimSq);	
      } else { GRIDCYCLE_Y(VARIO_MULT_Y); }				
    } else { // grid, y not given					
      if (vdim0 == 1) {							
	GRIDCYCLE(VARIO_UNIVAR; value+=vdimSq);				
      } else {GRIDCYCLE(VARIO_MULT); }				
    }									
} else { // not a grid 						
    if (trafo) {							
      Transform2NoGrid(cov, &xx, &yy, false);				
      x = xx;								
      if (ygiven) y = yy;						
    } else {								
      x=loc->x;								
      if (ygiven) y=loc->y;						
    }									
    assert(ygiven xor (y==ZERO));					
    if (ygiven || kernel) {						
      double *y0, *yend;						
      yend = ygiven ? y + tsxdim * loc->ly : ZERO;			
      y0 = y;								
      NONGRIDCYCLE(DO_INCREMENTY, DO_RECYCLY, VARIO_UNIVAR_Y; value+=vdimSq,
		   VARIO_MULT_Y);					
    } else { 
      //NONGRIDCYCLE(EMPTY, EMPTY, VARIO_UNIVAR; value+=vdimSq, VARIO_MULT); } 
      //#define NONGRIDCYCLE(INCREMENTY, RECYCLY, FCTN1, FCTN2)	

      if (vdim0 == 1) { 							
	for (; loc->i_row<tot; loc->i_row++, x+=tsxdim ){
	  p rintf("i=%d %ld %d\n", loc->i_row, value, tot); 
	  COV(x, cov, cross);		       
	  *value = *C0y - *cross;	
	  value+=vdimSq;
	}									
      } else {								
	for (; loc->i_row<tot; loc->i_row++, x+=tsxdim ){		
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
 
 

void CovarianceMatrix(cov_model *cov, double *v) {
  domain_type domown;
  Types type;
  genuineStatOwn(cov, &domown, &type);
  if (cov->pref[Nothing] == PREF_NONE ||
      (!isPosDef(type) && (!isVariogram(type) || domown!=XONLY))) {
    // Variogramme sind hier definitiv erlaubt, da durch Addition einer 
    // Konstanten, das Zeug zu einer Kovarianzmatrix gemacht werden kann
    // siehe direct.cc
    //assert(({PMI(cov, "cov matrix"); true;})); //
    ERR("covariance matrix: given model is not a covariance function");
  }
      
  STANDARDSTART;
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
  
#define MULTICOV							\
  C = v + VDIM * (loc->i_col + loc->i_row * vdim0tot);			\
  for (l=n=0; n<NEND; n+=NINCR) {					\
    long endfor = n + ENDFORINCR;					\
    for (m=n; m<endfor; m+=MINCR) {					\
      C[m] = z[l++];							\
    }									\
  }									\
									\
  if (loc->i_col != loc->i_row) {					\
    C = v + VDIM * (loc->i_row + loc->i_col * vdim0tot);		\
    for (l=m=0; m<ENDFORINCR; m+=MINCR) {				\
      for (n=m; n<NEND; n+=NINCR) {					\
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
	//	APMI(cov); 
	//	printf("nostat %d %d x=%f %f y=%f %f z=%f \n", loc->i_row, loc->i_col, *x, x[1], *y, y[1], *z);

	MULTICOV; 
	loc->i_row++;
	STANDARDINKREMENT_Y;
	//	APMI(cov);
	if (d >= tsxdim) break; 
      }
      loc->i_col++;  
      STANDARDINKREMENT;
    }
  } else { 
    int localdim = tsxdim; // just to control
    if (trafo) {
      localdim = Transform2NoGrid(cov, &xx, false);    
      assert(localdim == tsxdim);

      x0 = xx;
    } else x0 = loc->x;
    x = x0;

    for (loc->i_col=0; loc->i_col<tot; loc->i_col++, x+=localdim) {
      for (y=x, loc->i_row=loc->i_col; loc->i_row<tot; 
	   loc->i_row++, y+=localdim) {

	if (dist) {	  	  
	  // here x and y are ignored	  
	  
	  if (loc->i_row==loc->i_col) {
	    COV(ZERO, cov, z); // !! Achtung; wegen der Mixed models
	    //                    koennen hier verschiedene Werte auftreten
	    //                   aufgrund von loc->i_col und loc->i_row !!
	  } else {
	    COV(x0 + (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2
		      + loc->i_row -1) * localdim , cov, z);


	    //	 PMI(cov->calling->calling->calling);
	    //		    int idx = (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2  + loc->i_row -1) * localdim;
	    //  printf("loc=%d %f %f %f\n", idx, x[idx], x[idx +1], z);

	  }

	  	 
	  if (false) {
	    if (loc->i_col < 10 && loc->i_row<10) 
	      printf("cm %d %d totM1=%d %ld %f %f \n", (int) loc->i_row, (int)loc->i_col, (int) totM1, (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2 + loc->i_row -1) * localdim, x0[ (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2 + loc->i_row -1) * localdim] , *z); // else assert(false);
	  }

	} else {
	  NONSTATCOV(x, y, cov, z);
	  assert(R_FINITE(z[0]));
	}
	// printf("%s\n", NAME(cov));
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




int ptrStart(int *ptr, int *selected, int nsel, long tot, int vdim) {
  long v,i, step, start, newstart;
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


void ptrNext(int *ptr, int *selected, int nsel, long tot, int vdim, int* min) {
  long v,
    step = tot;
  int
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


void partial_loc_set_matrix(cov_model *cov, double *x, long lx, bool dist,
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

void partial_loc_set_matrixOZ(cov_model *cov, double *x, long lx, bool dist,
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



void partial_loc_set(cov_model *cov, double *x, long lx, bool dist, bool grid){
  location_type *loc = Loc(cov);
  int err;
  //  bool cartesian = isCartesian(cov->isoown);
  //  if (!cartesian && loc->ly==0) add_y_zero(loc);
  if ((err = partial_loc_set(loc, x, NULL, // cartesian ? NULL : ZERO, 
			     lx, 0, //!cartesian,
			     dist, loc->xdimOZ,
			     NULL, grid, false)) 
      != NOERROR) XERR(err);
}


void partial_loc_setOZ(cov_model *cov, double *x, double *y, 
		       long lx, long ly, bool dist, int *xdimOZ){
  // *xdimOZ to distinguish from the previous partial_loc_set definition
  location_type *loc = Loc(cov);
  int err;
  
  //  printf("partial_loc_set dist = %d %d \n", dist, loc->ly);
  
  if ((err = partial_loc_set(loc, x, y, lx, ly, dist, *xdimOZ, 
			     NULL, loc->grid, false)) 
      != NOERROR) XERR(err);
}


void partial_loc_setOZ(cov_model *cov, double *x, long lx, bool dist, int *xdimOZ){
  // *xdimOZ to distinguish from the previous partial_loc_set definition
  location_type *loc = Loc(cov);
  int err;
  //  bool cartesian = isCartesian(cov->isoown);
  // if (!cartesian && loc->ly==0) add_y_zero(loc);
  
  //  printf("partial_loc_set dist = %d %d \n", dist, loc->ly);
   
  if ((err = partial_loc_set(loc, x, NULL, // cartesian ? NULL : ZERO,  
			     lx, 0, //!cartesian,
			     dist, *xdimOZ, 
			     NULL, loc->grid, false)) 
      != NOERROR) XERR(err);
  //PMI(cov);
}



void partial_loc_setXY(cov_model *cov, double *x, double *y, long lx, long ly) {
  location_type *loc = Loc(cov);
  int err;

  // if (y == NULL)  crash();
  // assert(y != NULL);
 
  if ((err = partial_loc_set(loc, x, y, lx, ly, false,
			     loc->xdimOZ, NULL, loc->grid, 
			     false)) 
      != NOERROR) XERR(err);
}


void partial_loc_setXY(cov_model *cov, double *x, double *y, long lx) {
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


void InverseCovMatrix(cov_model *cov, double *v, double *det) {
  location_type *loc = Loc(cov);
  long vdimtot = loc->totalpoints * cov->vdim[0];  
  assert(cov->vdim[0] == cov->vdim[1]);
  CovList[cov->nr].covariance(cov, v);
  if (cov->Ssolve == NULL) SOLVE_STORAGE;
  Ext_setErrorLoc(ERROR_LOC);
  //  printf("inverse\n");
  int Exterr = Ext_solvePosDef_(v, vdimtot, true, NULL, 0, det, cov->Ssolve, 
			 &(GLOBAL.solve), PL);
  if (Exterr != NOERROR){
    Ext_getErrorString(ERRORSTRING);
    ErrorStop(Exterr);
  }
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
 
#define STANDARDINTERN_SEXP_BASIC						\
  if (INTEGER(reg)[0] < 0 || INTEGER(reg)[0] > MODEL_MAX) XERR(ERRORREGISTER); \
  if (currentNrCov==-1) InitModelList();				\
  cov_model *cov = KEY[INTEGER(reg)[0]];				\
  if (cov == NULL) { ERR("register not initialised") }			

#define STANDARDINTERN_SEXP						\
  STANDARDINTERN_SEXP_BASIC;						\
  cov_model VARIABLE_IS_NOT_USED *truecov =  !isInterface(cov)	?	\
    cov : cov->key == NULL ? cov->sub[0] : cov->key

//  if (cov->pref[Nothing] == PREF_NONE) { PMI(cov); XERR(ERRORINVALIDMODEL) }

SEXP Delete_y(SEXP reg) {
  STANDARDINTERN_SEXP_BASIC; 
  int d;
  location_type *loc = Loc(cov);
  if (loc->y != NULL) {
    if (loc->y != loc->x) UNCONDFREE(loc->y);
    loc->y = NULL;
  }
  if (loc->ygr[0] != NULL) {
    if (loc->ygr[0] != loc->xgr[0]) UNCONDFREE(loc->ygr[0]);
    for (d=0; d<MAXSIMUDIM; d++) loc->ygr[d] = NULL;
  }
  loc->ly = 0;
  return NULL;
}

SEXP CovLoc(SEXP reg, SEXP x, SEXP y, SEXP xdimOZ, SEXP lx,
	    SEXP result) {
  STANDARDINTERN_SEXP;

  //   PMI(cov);
  if (Loc(cov)->len > 1) BUG;  

  partial_loc_setXY(cov, REAL(x), TYPEOF(y) == NILSXP ? NULL : REAL(y),
		    INTEGER(lx)[0]);
  CovList[truecov->nr].covariance(truecov, REAL(result));
  // PMI(cov, "covloc"); printf("xdimOZ %d \n", INTEGER(xdimOZ)[0]);
  partial_loc_null(cov);
  if (Loc(cov)->xdimOZ != INTEGER(xdimOZ)[0]) BUG;
  return NULL;
}

void CovIntern(int reg, double *x, double *y, long lx, long ly, double *value) {
  STANDARDINTERN;
  partial_loc_setXY(cov, x, y, lx, ly);
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


SEXP VariogramIntern(SEXP reg) {
  STANDARDINTERN_SEXP;
  //  location_type *loc = Loc(cov);
  SEXP ans;
  int 
    vdim =cov->vdim[0],
    size = Gettotalpoints(cov) * vdim * vdim;
  PROTECT(ans = allocVector(REALSXP, size));
  CovList[truecov->nr].variogram(truecov, REAL(ans)); 
  UNPROTECT(1);
  return(ans);
}


void PseudovariogramIntern(int reg, double *x, double *y,
			   long lx, long ly, double *value) {
  STANDARDINTERN;
  location_type *loc = Loc(cov);
  partial_loc_setOZ(cov, x, y, lx, ly, false, &(loc->xdimOZ));
  CovList[truecov->nr].pseudovariogram(truecov, value);
  partial_loc_null(cov);
}

void PseudovariogramIntern(int reg, double *x, double *value) {
  STANDARDINTERN;
  location_type *loc = Loc(cov);
  // bool cartesian = isCartesian(cov->isoown);
  // if (!cartesian && loc->ly==0) add_y_zero(loc);
  partial_loc_setOZ(cov, x, NULL, // cartesian ? NULL : ZERO, 
		    1, 0, // !cartesian,
		    false, &(loc->xdimOZ));
  CovList[truecov->nr].pseudovariogram(truecov, value);
  partial_loc_null(cov);
}

