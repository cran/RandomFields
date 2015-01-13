#ifndef variogramAndCo_H
#define variogramAndCo_H 1


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
  if (xx != NULL) free(xx);			\
  if (yy != NULL) free(yy);			\
   loc->i_col = loc->i_row = 0;
  // standard ende


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

