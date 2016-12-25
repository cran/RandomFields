/* SpherModels.h */
#ifndef SPHERMODELS_H
#define SPHERMODELS_H 1


// declaration of functions
void SinePower(double *x, cov_model *cov, double *v);
void rangeSinePower(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range);

void Multiquad(double *x, cov_model *cov, double *v);
void rangeMultiquad(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range);
/*
  void kappa_choquet(int i, cov_model *cov, int *nr, int *nc);
  void Choquet(double *x, cov_model *cov, double *v);
  void rangeChoquet(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range);
  int checkChoquet(cov_model *cov);
*/


#endif /* SpherModels.h */
