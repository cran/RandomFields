#ifndef extremes_H
#define extremes 1
#include "RFsimu.h"

extern double EXTREMES_STANDARDMAX;

#ifdef CPP2C 
EXTERN void InitMaxStableRF(double *x, double *y, double *z, int *dim, int *lx, 
			    int *grid, int *covnr, double *ParamList, int *nParam, 
			    int *method, 
			    int *distr,
			    int *keyNr,
			    int *error);

EXTERN void DoMaxStableRF(int *keyNr, double *res, int *error); 


#else
extern void InitMaxStableRF(double *x, double *y, double *z, int *dim, int *lx, 
			    int *grid, int *covnr, double *ParamList, int *nParam, 
			    int *method, 
			    int *distr,
			    int *keyNr,
			    int *error);

extern void DoMaxStableRF(int *keyNr, double *res, int *error);
#endif /* CPP2C */

#endif /* extremes */
