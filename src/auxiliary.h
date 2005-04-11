
#ifndef AUXILIARY_H
#define AUXILIARY_H 1

#include "GSLvsR.h"

EXTERN void xyInchToPixels(double *xy, double *Pixels);
EXTERN void RandomPermutation(double *x, int n, double *y 
				  );
EXTERN double quantile(double *X, int lb, double p);
EXTERN void Rquantile(double *X, int *lb, double *p, double *res);
EXTERN void pid(int *i);
EXTERN void hostname(char **h, int *i);
EXTERN void orderdouble(double *d,int *pos, int start, int end);
EXTERN void quicksortdouble(double *d, int start, int end);
EXTERN void vectordist(double *v, int *dim, double *dist, int *diag); 


EXTERN void I0ML0(double *x, int *n);
double I0mL0(double x);
EXTERN double struve(double x, double nu,  double factor_sign, bool expscaled);
EXTERN void StruveH(double *x, double *nu);
EXTERN void StruveL(double *x, double *nu, int * expScaled);
EXTERN int InitRandom();

EXTERN void ordering(double *d, int len, int dim, int *pos);
EXTERN void Ordering(double *d, int *len, int *dim, int *pos);

#endif /* AUXILIARY_H */




