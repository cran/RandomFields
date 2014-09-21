
#ifndef AUXILIARY_H
#define AUXILIARY_H 1

#include "basic.h"

double I0mL0(double x);
extern "C" {
  SEXP vectordist(SEXP V, SEXP diag); 
  void sleepMilli(int *milli);
  void I0ML0(double *x, int *n);
  double struve(double x, double nu,  double factor_sign, bool expscaled);
  void StruveH(double *x, double *nu);
  void StruveL(double *x, double *nu, int * expScaled);
  void ordering(double *d, int len, int dim, int *pos);
  void Ordering(double *d, int *len, int *dim, int *pos);
}
bool is_diag(double *aniso, int dim);

#endif /* AUXILIARY_H */




