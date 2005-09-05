
#ifndef AUXILIARY_H
#define AUXILIARY_H 1

#include "GSLvsR.h"

EXTERN void vectordist(double *v, int *dim, double *dist, int *diag); 
EXTERN void sleepMilli(int *milli);

EXTERN void I0ML0(double *x, int *n);
double I0mL0(double x);
EXTERN double struve(double x, double nu,  double factor_sign, bool expscaled);
EXTERN void StruveH(double *x, double *nu);
EXTERN void StruveL(double *x, double *nu, int * expScaled);

EXTERN void ordering(double *d, int len, int dim, int *pos);
EXTERN void Ordering(double *d, int *len, int *dim, int *pos);

bool is_diag(double *aniso, int dim);

#endif /* AUXILIARY_H */




