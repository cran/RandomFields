/* SpherModels.h */
#ifndef SPHERMODELS_H
#define SPHERMODELS_H 1

#ifdef __cplusplus
extern "C"
{
#endif

  // declaration of functions
  void SinePower(double *x, cov_model *cov, double *v);
  void rangeSinePower(cov_model *cov, range_type *range);

  void Eq16(double *x, cov_model *cov, double *v);
  void rangeEq16(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range);

  
#ifdef __cplusplus
}
#endif


#endif /* SpherModels.h */
