#ifndef variogramAndCo_H
#define variogramAndCo_H 1


#define FINISH_START(whereSpgs)						\
  assert(whereSpgs->Spgs != NULL);					\
  pgs_storage *pgs= whereSpgs->Spgs;					\
  assert(pgs->x != NULL);						\
  location_type *loc = Loc(cov);					\
  assert(loc != NULL);							\
  bool grid =loc->grid && !loc->Time && loc->caniso==NULL,		\
    trafo = (loc->Time || loc->caniso != NULL) && !loc->distances,	\
    VARIABLE_IS_NOT_USED dist = loc->distances,				\
    VARIABLE_IS_NOT_USED isposdef = isPosDef(cov),			\
    ygiven = loc->y != NULL || loc->ygr[0] != NULL;/* might be changed !*/ \
  long cumgridlen[MAXMPPDIM +1],					\
    err = NOERROR;							\
  int d,								\
    tsdim VARIABLE_IS_NOT_USED = loc->timespacedim,			\
    vdim0 VARIABLE_IS_NOT_USED = cov->vdim[0],				\
    VARIABLE_IS_NOT_USED vdimSq =cov->vdim[0]*cov->vdim[1], /*not nec. squ!!*/ \
    *gridlen=pgs->gridlen,						\
    *end=pgs->end,							\
    *endy VARIABLE_IS_NOT_USED =pgs->endy,				\
    *start=pgs->start,							\
    *startny VARIABLE_IS_NOT_USED =pgs->startny,			\
    *nx=pgs->nx,							\
    *ny VARIABLE_IS_NOT_USED = pgs->delta,				\
    VARIABLE_IS_NOT_USED xdimOZ = loc->xdimOZ, /* weicht u.U. von tsdim bei dist=true ab */ \
    tsxdim = cov->xdimprev; /* weicht u.U. von tsdim bei dist=true ab */ \
  long tot = loc->totalpoints,						\
    VARIABLE_IS_NOT_USED totM1 = tot - 1,				\
    VARIABLE_IS_NOT_USED vdim0P1 = vdim0 + 1,				\
    vdim0tot = vdim0 * tot,						\
    VARIABLE_IS_NOT_USED vdimSqtot = vdimSq * tot,			\
    VARIABLE_IS_NOT_USED vdim0totSq = vdim0tot * tot,			\
    VARIABLE_IS_NOT_USED vdimSqtotSq = vdimSq * tot * tot;		\
  double *x = pgs->x,							\
    *xstart= pgs->xstart,						\
    *inc=pgs->inc,							\
    *ystart VARIABLE_IS_NOT_USED =pgs->supportmax,			\
    *yy = NULL,		/* set by Transform2NoGrid; freed */		\
    *xx = NULL,		/* dito			   */			\
    *x0 VARIABLE_IS_NOT_USED = NULL,	/* never free it  */		\
    *incy VARIABLE_IS_NOT_USED =pgs->supportcentre,			\
    *z VARIABLE_IS_NOT_USED = pgs->z;					\
  if (grid) {								\
    cumgridlen[0] = 1;							\
    for (d=0; d<tsxdim; d++) {						\
      inc[d] = loc->xgr[d][XSTEP];					\
      gridlen[d] = loc->xgr[d][XLENGTH];				\
      cumgridlen[d+1] = gridlen[d] * cumgridlen[d];			\
									\
      start[d] = 0;							\
      end[d] = gridlen[d];						\
      nx[d] = start[d];							\
      x[d] = xstart[d] = loc->xgr[d][XSTART];				\
    }									\
  }									\
  loc->i_row = 0;							\
  loc->i_col = I_COL_NA


#define INCLUDE_VAL							\
  long i, ii, v, j;							\
  double **Val= pgs->Val,						\
    *cross=pgs->cross;							\
  if (vdim0 > 1) {							\
    for (v = i = 0; i<vdimSqtot ; i+=vdim0tot) {			\
      for (ii=i, j=0; j<vdim0; ii+=tot, v++, j++) {			\
	Val[v] = value + ii;						\
      }									\
    }									\
  }									\
  if (kernel && !ygiven && PL > 0) {					\
    char wrn[200]; 				\
    sprintf(wrn, "'%s' is called with a single variable only, although it is used as a kernel. So, the second variable is set to zero, here.\n", NICK(cov)); \
    warning(wrn);							\
  }
 


#define STANDARD_START(NX, START, END, X, XSTART, INC)	\
  int 						\
    *nx = NX,						\
    *start = START,					\
    *end = END;						\
  double						\
    *x = X,						\
    *xstart = XSTART,					\
    *inc = INC;						\
							\
  for (d=0; d<tsxdim; d++) {				\
    nx[d] = start[d];					\
    x[d] = xstart[d];					\
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
  if (d >= tsxdim) break		
  // end define StandardInkrement


				
// Achtung hier *nicht* incy, starty
#define STANDARDINKREMENT_Y	\
  d = 0;	\
  ny[d]++;			\
  y[d] += incy[d]; 	\
   while (ny[d] >= endy[d]) {	\
     ny[d] = startny[d];	\
     y[d] = ystart[d];	\
     if (++d >= tsxdim) break;	\
    ny[d]++;			\
    y[d] += incy[d]; 		\
   }				
  // end define StandardInkrement_Y


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
								

#define RECYCLE_Y					\
  if (d>tsxdim)						\
    for (d=0; d<tsxdim; d++) {				\
      y[d] = ystart[d] = loc->ygr[d][XSTART];		\
      ny[d] = startny[d];				\
    }
				

// loc->grid bei SelectedCovMatrix notwendig
#define STANDARD_ENDE				\
  FREE(xx);			\
  FREE(yy);			\
  loc->i_col = loc->i_row = I_COL_NA;


#define GRIDCYCLE_Y(SOME_FCTN)			\
    while (true) {				\
      SOME_FCTN;				\
      loc->i_row++;				\
      if (ygiven) { STANDARDINKREMENT_Y; RECYCLE_Y;}	\
      STANDARDINKREMENT;			\
    }					

#define GRIDCYCLE(SOME_FCTN)			\
    while (true) {				\
      SOME_FCTN;				\
      loc->i_row++;				\
      STANDARDINKREMENT;			\
    }
    
#define DO_INCREMENTY , y+=tsxdim
#define DO_RECYCLY if (y >= yend) y = y0;
#define EMPTY 
#define NONGRIDCYCLE(INCREMENTY, RECYCLY, FCTN1, FCTN2)			\
  if (vdim0 == 1) {							\
    for (; loc->i_row<tot; loc->i_row++, x+=tsxdim INCREMENTY){		\
      RECYCLY			\
      FCTN1;								\
    }									\
  } else {								\
    for (; loc->i_row<tot; loc->i_row++, x+=tsxdim INCREMENTY){		\
      RECYCLY								\
      FCTN2;								\
    }									\
  }

#define PERFORM(UNIVAR_FCTN, MULTIVAR_FCTN, UNIVAR_FCTN_Y, MULTIVAR_FCTN_Y)  \
if (grid) {					\
  if (ygiven || kernel) { 		\
      STANDARDSTART_Y_SUPPL;			\
      if (vdim0 == 1) { GRIDCYCLE_Y(UNIVAR_FCTN_Y; value+=vdimSq);	\
      } else { GRIDCYCLE_Y(MULTIVAR_FCTN_Y); }				\
    } else { /* grid, y not given */					\
      if (vdim0 == 1) {							\
	GRIDCYCLE(UNIVAR_FCTN; value+=vdimSq);				\
      } else {GRIDCYCLE(MULTIVAR_FCTN); }				\
    }									\
  } else { /* not a grid */						\
    if (trafo) {							\
      Transform2NoGrid(cov, &xx, &yy);					\
      x = xx;								\
      if (ygiven) y = yy;						\
    } else {								\
      x=loc->x;								\
      if (ygiven) y=loc->y;						\
    }									\
    assert(ygiven xor (y==ZERO));					\
    if (ygiven || kernel) {						\
      double *y0, *yend;						\
      yend = ygiven ? y + tsxdim * loc->ly : ZERO;			\
      y0 = y;								\
      NONGRIDCYCLE(DO_INCREMENTY, DO_RECYCLY, UNIVAR_FCTN_Y; value+=vdimSq,\
		   MULTIVAR_FCTN_Y);					\
    } else { NONGRIDCYCLE(EMPTY, EMPTY, UNIVAR_FCTN; value+=vdimSq, MULTIVAR_FCTN); } \
    STANDARD_ENDE;							\
    if (err != NOERROR) XERR(err);					\
  }




void CovVario(cov_model *cov, bool is_cov, bool pseudo, 
	      double *value);
void CovIntern(int reg, double *x, double *value);
void CovIntern(int reg, double *x, double *y, long lx, long ly, double *value);
void PseudovariogramIntern(int reg, double *x, double *y,
			   long lx, long ly, double *value);
void PseudovariogramIntern(int reg, double *x, double *value);
void CovarianceMatrix(cov_model *cov,double *COV, int *nonzeros);
void InverseCovMatrix(cov_model *cov, double *v);
void SelectedCovMatrix(cov_model *cov, int *selected /* nr! */, int nsel, 
		       double *v, int *nonzeros);
void StandardCovariance(cov_model *cov, double *v);
void StandardCovMatrix(cov_model *cov, double *v, int *nonzeros);
void StandardInverseCovMatrix(cov_model *cov, double *v);
void StandardVariogram(cov_model *cov, double *v);
void StandardPseudoVariogram(cov_model *cov, double *v);
void StandardSelectedCovMatrix(cov_model *cov, int *sel, int nsel, double *v,
			       int* nonzeros);


#endif

