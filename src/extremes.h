#ifndef extremes_H
#define extremes 1
#include "RFsimu.h"

extern Real EXTREMES_STANDARDMAX;

#ifdef CPP2C 
EXTERN void InitMaxStableRF(Real *x, Real *y, Real *z, int *dim, int *lx, 
			    int *grid, int *covnr, Real *ParamList, int *nParam, 
			    int *method, 
			    int *distr,
			    int *keyNr,
			    int *error);

EXTERN void DoMaxStableRF(int *keyNr, Real *res, int *error); 


#else
extern void InitMaxStableRF(Real *x, Real *y, Real *z, int *dim, int *lx, 
			    int *grid, int *covnr, Real *ParamList, int *nParam, 
			    int *method, 
			    int *distr,
			    int *keyNr,
			    int *error);

extern void DoMaxStableRF(int *keyNr, Real *res, int *error);
#endif /* CPP2C */

#endif /* extremes */
