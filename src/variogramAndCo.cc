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
#include "Operator.h"
#include "variogramAndCo.h"

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
  FINISH_START(prev);							\
  assert(cov->vdim[0] == cov->vdim[1]);	
			



void CovVario(cov_model *cov, bool is_cov, bool pseudo, double *value) { 
  STANDARDSTART;
  domain_type domown;
  Types type;
  long m, n;
  bool stat;
  double 
    *y = ygiven ? pgs->supportmin : ZERO,
    *C0x = pgs->C0x,
    *C0y = pgs->C0y;
  bool kernel = cov->domown != XONLY;  
  INCLUDE_VAL;				       

  // if(isCartesian(prev->isoown) || ygiven);


  genuineStatOwn(cov, &domown, &type);
  stat = !kernel && isVariogram(type);

  if (is_cov) {
    assert(({/*PMI(cov, "Cov/vario"); */ cov->pref[Nothing] != PREF_NONE && isShape(type);})); //
  } else {
    if (cov->pref[Nothing] == PREF_NONE || !isVariogram(type)) {
      assert(({PMI(cov, "cov/Vario"); true;})); //
      ERR("given model is not a variogram");
    }
    if (stat && isCartesian(prev->isoown)) {
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
 

  if (is_cov) {
    PERFORM(UNIVAR, MULT, UNIVAR_Y, MULT_Y);
  } else if (pseudo) {
    PERFORM(VARIO_UNIVAR, PSEUDO_MULT, VARIO_UNIVAR_Y, PSEUDO_MULT_Y);
  } else {
    PERFORM(VARIO_UNIVAR, VARIO_MULT, VARIO_UNIVAR_Y, VARIO_MULT_Y);
  }
} 
 
 

void CovarianceMatrix(cov_model *cov, double *v, int *nonzeros) {
   domain_type domown;
  Types type;
  genuineStatOwn(cov, &domown, &type);
  if (cov->pref[Nothing] == PREF_NONE ||
      (!isPosDef(type) && (!isVariogram(type) || domown!=XONLY))) {
    // Variogramme sind hier definitiv erlaubt, da durch Addition einer 
    // Konstanten, das Zeug zu einer Kovarianzmatrix gemacht werden kann
    // siehe direct.cc
    //assert(({PMI(cov, "cov matrix"); true;})); //
    error("covariance matrix: given model is not a covariance function");
  }
  long l,n,m, VDIM, NEND, NINCR, MINCR, ENDFORINCR;
  int nz = 0;
  double *C = NULL;
  bool vdim_closetogether = GLOBAL.general.vdim_close_together;
  //
  
  STANDARDSTART;
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
      nz += z[l] != 0.0;						\
      C[m] = z[l++];							\
    }									\
  }									\
									\
  if (loc->i_col != loc->i_row) {					\
    C = v + VDIM * (loc->i_row + loc->i_col * vdim0tot);		\
    for (l=m=0; m<ENDFORINCR; m+=MINCR) {				\
      for (n=m; n<NEND; n+=NINCR) {					\
	nz += z[l] != 0.0;						\
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
	    //		    int idx = (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2  + loc->i_row -1) * tsxdim;
	    //  printf("loc=%d %f %f %f\n", idx, x[idx], x[idx +1], z);

	  }

	  	 
	  if (false) {
	    if (loc->i_col < 10 && loc->i_row<10) 
	      printf("cm %d %d totM1=%d %ld %f %f \n", (int) loc->i_row, (int)loc->i_col, (int) totM1, (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2 + loc->i_row -1) * tsxdim, x0[ (loc->i_col * totM1 - (loc->i_col * (loc->i_col + 1)) / 2 + loc->i_row -1) * tsxdim] , *z); // else assert(false);
	  }

	} else {
	  
	  //	  printf("%d %s\n", cov->gatternr, NICK(cov));
	  //assert(cov->gatternr == 5);

	  // PMI(cov->calling);

	  NONSTATCOV(x, y, cov, z);
	  //		  printf("x=%f %f y=%f %f z=%f\n", x[0], x[1], y[0], y[1], z[0]);

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

  *nonzeros = nz;
  
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

void SelectedCovMatrix(cov_model *cov,
		       int *selected /* nr! */, int nsel, double *v,
		       int *nonzeros) {
  bool vdim_close_together = GLOBAL.general.vdim_close_together;
  domain_type domown;
  Types type;
  int l,vrow, vcol, 
    nz = 0,
    oldCMC = NA_INTEGER, 
    oldCMR = NA_INTEGER;
  //   vdimselrow, vdimselcol, selrow, selcol, 

  genuineStatOwn(cov, &domown, &type);
  if (cov->pref[Nothing] == PREF_NONE || !isPosDef(type)) {
    assert(({PMI(cov, "cov matrix"); true;})); //
    error("cov. matrix: given model is not a covariance function");
  }
  if (vdim_close_together) {
    NotProgrammedYet("vdim_closetogether for selected access on covariance matrices");
  }
  
    
  STANDARDSTART;
  int *ptrcol = pgs->ptrcol, 
    *ptrrow = pgs->ptrrow;
  double *y = pgs->supportmin;

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
  
  loc->i_col = ptrStart(ptrcol, selected, nsel, tot, vdim0);
  
  while (loc->i_col < tot) {
    //  print("col=%d\n", loc->i_col);
    if (grid) split(loc->i_col, tsxdim, cumgridlen, inc, x);
    else x = x0 + tsxdim * loc->i_col; 
    MEMCOPY(ptrrow, ptrcol, sizeof(int) * vdim0);
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
      
      for (vcol=l=0; vcol<vdim0; vcol++) {
	//
	//	print(">> %d %d sel=%d tot=%d i=%d %d j=%d\n",
	//	      vcol, ptrcol[vcol], selected[ptrcol[vcol]],
	//	      tot, loc->i_col, selected[ptrcol[vcol]] % tot == loc->i_col,
	//     loc->i_row);
	
	if (ptrcol[vcol] >= 0 && 
	    selected[ptrcol[vcol]] % tot == loc->i_col) {
	  for (vrow=0; vrow<vdim0; vrow++) {
	    //	    print("vcol=%d, vrow=%d %d %d j=%d %d z=%f, %d (%d %d)\n",
	    //		  vcol, vrow, ptrrow[vrow],
	    //	   	  selected[ptrrow[vrow]], loc->i_row, 
	    //		  selected[ptrrow[vrow]] % tot == loc->i_row,
	    //	  z[l], l,
	    //		  ptrcol[vcol] * nsel + ptrrow[vrow],
	    //	  ptrcol[vcol] + ptrrow[vrow] * nsel );
	    if (ptrrow[vrow] >= 0 && 
		selected[ptrrow[vrow]] % tot == loc->i_row) {	      
	      nz += z[l] != 0.0;
	      v[ptrcol[vcol] * nsel + ptrrow[vrow]] = 
		v[ptrcol[vcol] + ptrrow[vrow] * nsel] = z[l++];
	      //   print("%d %d %d %d \n", loc->i_col, loc->i_row,
	      //          ptrcol[vcol], ptrrow[vrow]);
	    }
	  } // for vrow
	}
      } // for vcol
      ptrNext(ptrrow, selected, nsel, tot, vdim0, &loc->i_row);
    } // while loc->i_row < tot
    ptrNext(ptrcol, selected, nsel, tot, vdim0, &loc->i_col);  
    // assert( loc->i_col < 5);
  } // while loc->i_col < tot
  
  if (cov->domown==XONLY && isVariogram(cov->typus) && 
      !isPosDef(cov->typus)) {
    double first=v[0];
    nz = 0;
    for (l=0; l<vdimSqtotSq; l++) {
      v[l] = first - v[l];
      nz += v[l] != 0.0;
    }
  }

  *nonzeros = nz;
  
 ErrorHandling:
  STANDARD_ENDE; // 	ErrorHandling:
  
  loc->i_col=oldCMC;
  loc->i_row=oldCMR; 
  
  if (err!=NOERROR) XERR(err); 
    
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
  bool cartesian = isCartesian(cov->isoown);
  if (!cartesian && loc->ly==0) add_y_zero(loc);
  if ((err = partial_loc_set(loc, x, cartesian ? NULL : ZERO, 
			     lx, !cartesian, dist, loc->xdimOZ,
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
  bool cartesian = isCartesian(cov->isoown);
  if (!cartesian && loc->ly==0) add_y_zero(loc);
  
  //  printf("partial_loc_set dist = %d %d \n", dist, loc->ly);
   
  if ((err = partial_loc_set(loc, x, cartesian ? NULL : ZERO,  
			     lx, !cartesian, dist, *xdimOZ, 
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


void InverseCovMatrix(cov_model *cov, double *v) {
  location_type *loc = Loc(cov);
  long vdimtot = loc->totalpoints * cov->vdim[0];  
  assert(cov->vdim[0] == cov->vdim[1]);
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
  cov_model *cov = KEY[INTEGER(reg)[0]];				\
  if (cov == NULL) { ERR("register not initialised") }			\
  cov_model VARIABLE_IS_NOT_USED *truecov =  !isInterface(cov)	?	\
    cov : cov->key == NULL ? cov->sub[0] : cov->key

//  if (cov->pref[Nothing] == PREF_NONE) { PMI(cov); XERR(ERRORINVALIDMODEL) }

SEXP Delete_y(SEXP reg) {
  STANDARDINTERN_SEXP; 
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
			  SEXP result, SEXP nonzeros) {
  STANDARDINTERN_SEXP;
  partial_loc_set_matrixOZ(cov, REAL(x), INTEGER(lx)[0], LOGICAL(dist)[0], 
			   INTEGER(xdimOZ));
  CovList[truecov->nr].selectedcovmatrix(truecov, INTEGER(selected), 
					 INTEGER(nsel)[0], REAL(result),
					 INTEGER(nonzeros));
  partial_loc_null(cov);
  if (Loc(cov)->xdimOZ != INTEGER(xdimOZ)[0]) BUG;
  
  return NULL;
}


SEXP CovMatrixSelected(SEXP reg, SEXP selected, SEXP nsel,
		       SEXP result, SEXP nonzeros) {
  
  STANDARDINTERN_SEXP;
  CovList[truecov->nr].selectedcovmatrix(truecov, INTEGER(selected), 
					 INTEGER(nsel)[0], REAL(result),
					 INTEGER(nonzeros));
  
  return NULL;
}

/*

  .C("setListElements",Reg, nteger(1), as.integer(1),
  PACKAGE="RandomFields");

*/      

void CovIntern(int reg, double *x, double *y, long lx, long ly, double *value) {
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
  
  // PMI(cov);


  partial_loc_setXY(cov, x, cartesian ? NULL : ZERO, 1, !cartesian);
  CovList[truecov->nr].covariance(truecov, value);
  partial_loc_null(cov);
}


SEXP CovMatrixLoc(SEXP reg, SEXP x, SEXP dist, SEXP xdimOZ,
		  SEXP lx, SEXP result, SEXP nonzeros) {
  STANDARDINTERN_SEXP;
  partial_loc_set_matrixOZ(cov, REAL(x), INTEGER(lx)[0], LOGICAL(dist)[0], 
			 INTEGER(xdimOZ)); //

  //  PMI(cov, "covmatrixloc");

  CovList[truecov->nr].covmatrix(truecov, REAL(result), INTEGER(nonzeros));
 
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
		     SEXP lx, SEXP result, SEXP nonzeros) {
  
  STANDARDINTERN_SEXP;
  //PMI(truecov);
  partial_loc_set_matrix(cov, REAL(x), INTEGER(lx)[0], LOGICAL(dist)[0],
			 LOGICAL(grid)[0]);
  //PMI(cov);  
  //PMI(truecov);
  //assert(cov->xdimprev == 2);
  CovList[truecov->nr].covmatrix(truecov, REAL(result), INTEGER(nonzeros));
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
  bool cartesian = isCartesian(cov->isoown);
  if (!cartesian && loc->ly==0) add_y_zero(loc);
  partial_loc_setOZ(cov, x, cartesian ? NULL : ZERO, 
		    1, !cartesian, false, &(loc->xdimOZ));
  CovList[truecov->nr].pseudovariogram(truecov, value);
  partial_loc_null(cov);
}

