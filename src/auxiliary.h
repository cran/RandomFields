

#ifndef AUXILIARY_H
#define AUXILIARY_H 1
 
#include "RFsimu.h"

#ifdef CPP2C 
#ifdef RF_GSL
extern "C" void RandomPermutation(double *x,int n,double *y,gsl_rng *RANDOM);
#else 
extern "C" void RandomPermutation(double *x,int n,double *y);
#endif
extern "C" double quantile(double *X, int lb,double p);
extern "C" void Rquantile(double *X, int *lb,double *p,double *res);
extern "C" void pid(int * i);
extern "C" void hostname(char **h, int *i);
extern "C" void orderdouble(double *d,int *pos, int start, int end);
extern "C" void gauss(int *n, double *G);
#else
extern double quantile(double *X, int lb,double p);
extern void pid(int * i);
extern void hostname(char **h, int *i);
extern void orderdouble(double *d,int *pos, int start, int end);
extern void gauss(int *n, double *G);
#endif
#endif /* AUXILIARY_H */




