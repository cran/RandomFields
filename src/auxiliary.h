
#ifndef AUXILIARY_H
#define AUXILIARY_H 1

#include "basic.h"

extern "C" void vectordist(double *v, int *dim, double *dist, int *diag); 
extern "C" void sleepMilli(int *milli);

extern "C" void I0ML0(double *x, int *n);
double I0mL0(double x);
extern "C" double struve(double x, double nu,  double factor_sign, bool expscaled);
extern "C" void StruveH(double *x, double *nu);
extern "C" void StruveL(double *x, double *nu, int * expScaled);

extern "C" void ordering(double *d, int len, int dim, int *pos);
extern "C" void Ordering(double *d, int *len, int *dim, int *pos);

bool is_diag(double *aniso, int dim);

#endif /* AUXILIARY_H */




