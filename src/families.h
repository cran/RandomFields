#ifndef Families_H
#define Families_H 1

void arcsqrtD(double *x, cov_model *cov, double *v);
void arcsqrtDlog(double *x, cov_model *cov, double *v);
void arcsqrtDinverse(double *v, cov_model *cov, double *left, double *right);
void arcsqrtP(double *x, cov_model *cov, double *v);
void arcsqrtQ(double *x, cov_model *cov, double *v);
void arcsqrtR(double *x, cov_model *cov, double *v);
int check_arcsqrt_distr(cov_model *cov);
int init_arcsqrt(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s);
void do_arcsqrt(cov_model *cov, double *v);
void range_arcsqrt(cov_model *cov, range_type* range);


void determD(double *x, cov_model *cov, double *v);
void determDlog(double *x, cov_model *cov, double *v);
void determDinverse(double *v, cov_model *cov, double *left, double *right);
void determP(double *x, cov_model *cov, double *v);
void determP2sided(double *x, double *y, cov_model *cov, double *v);
void determQ(double *x, cov_model *cov, double *v);
void determR(double *x, cov_model *cov, double *v);
void determR2sided(double *x, double *y, cov_model *cov, double *v);
void kappa_determ(int i, cov_model *cov, int *nr, int *nc);
int check_determ(cov_model *cov);
int init_determ(cov_model *cov, gen_storage *s);
void do_determ(cov_model *cov, double *);
void range_determ(cov_model *cov, range_type *range);


void distrD(double *x, cov_model *cov, double *v);
void distrDlog(double *x, cov_model *cov, double *v);
void distrDinverse(double *v, cov_model *cov, double *left, double *right);
void distrP(double *x, cov_model *cov, double *v);
void distrP2sided(double *x, double *y, cov_model *cov, double *v);
void distrQ(double *x, cov_model *cov, double *v);
void distrR(double *x, cov_model *cov, double *v);
void distrR2sided(double *x, double *y, cov_model *cov, double *v);
int check_distr(cov_model *cov);
int init_distr(cov_model *cov, gen_storage *s);
void do_distr_do(cov_model *cov, double *);
void range_distr(cov_model *cov, range_type *range);



void sphericD(double *x, cov_model *cov, double *v);
void sphericDlog(double *x, cov_model *cov, double *v);
void sphericDinverse(double *v, cov_model *cov, double *left, double *right);
void sphericP(double *x, cov_model *cov, double *v);
void sphericQ(double *x, cov_model *cov, double *v);
void sphericR(double *x, cov_model *cov, double *v);
int check_RRspheric(cov_model *cov);
int init_RRspheric(cov_model *cov, gen_storage *s);
void do_RRspheric(cov_model *cov, double *);
void range_RRspheric(cov_model *cov, range_type *range);


void gaussD(double *x, cov_model *cov, double *v);
void gaussDlog(double *x, cov_model *cov, double *v);
void gaussDinverse(double *v, cov_model *cov, double *left, double *right);
void gaussP(double *x, cov_model *cov, double *v);
void gaussP2sided(double *x, double *y, cov_model *cov, double *v);
void gaussQ(double *x, cov_model *cov, double *v);
void gaussR(double *x, cov_model *cov, double *v);
void gaussR2sided(double *x, double *y, cov_model *cov, double *v);
void kappa_gauss_distr(int i, cov_model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc);
int check_gauss_distr(cov_model *cov);
int init_gauss_distr(cov_model *cov, gen_storage *s);
void do_gauss_distr(cov_model *cov, double *v);
void range_gauss_distr(cov_model *cov, range_type *range);


void locD(double *x, cov_model *cov, double *v);
void locDlog(double *x, cov_model *cov, double *v);
void locDinverse(double *v, cov_model *cov, double *left, double *right);
void locP(double *x, cov_model *cov, double *v);
void locP2sided(double *x, double *y, cov_model *cov, double *v);
void locQ(double *x, cov_model *cov, double *v);
void locR(double *x, cov_model *cov, double *v);
void locR2sided(double *x, double *y, cov_model *cov, double *v);
void kappa_loc(int i, cov_model *cov, int *nr, int *nc);
int check_loc(cov_model *cov);
int init_loc(cov_model *cov, gen_storage *s);
void do_loc(cov_model *cov, double *v);
void range_loc(cov_model *cov, range_type *range);

void rectangularD(double *x, cov_model *cov, double *v);
void rectangularDlog(double *x, cov_model *cov, double *v);
void rectangularDinverse(double *v, cov_model *cov, double *left,double *right);
void rectangularP(double *x, cov_model *cov, double *v);
void rectangularP2sided(double *x, double *y, cov_model *cov, double *v);
void rectangularQ(double *x, cov_model *cov, double *v);
void rectangularR(double *x, cov_model *cov, double *v); 
void rectangularR2sided(double *x, double *y, cov_model *cov, double *v); 
//void kappa_rectangular(int i, cov_model *cov, int *nr, int *nc);
int check_rectangular(cov_model *cov);
void do_rectangular(cov_model *cov, double *v);
int init_rectangular(cov_model *cov, gen_storage *s);
void range_rectangular(cov_model *cov, range_type *range);


void setParamD(double *x, cov_model *cov, double *v); 
void setParamDlog(double *x, cov_model *cov, double *v); 
void setParamDinverse(double *v, cov_model *cov, double *left, double *right);
void setParamP(double *x, cov_model *cov, double *v);   
void setParamP2sided(double *x, double *y, cov_model *cov, double *v);   
void setParamQ(double *x, cov_model *cov, double *v);   
void setParamR(double *x, cov_model *cov, double *v);   
void setParamR2sided(double *x, double *y, cov_model *cov, double *v);   
int check_setParam(cov_model *cov);
void range_setParam(cov_model *cov, range_type *range);
int init_setParam(cov_model *cov, gen_storage *s);
void do_setParam(cov_model *cov, double *v);




void unifD(double *x, cov_model *cov, double *v);
void unifDlog(double *x, cov_model *cov, double *v);
void unifDinverse(double *v, cov_model *cov, double *left, double *right);
void unifP(double *x, cov_model *cov, double *v);
void unifP2sided(double *x, double *y, cov_model *cov, double *v);
void unifQ(double *x, cov_model *cov, double *v);
void unifR(double *x, cov_model *cov, double *v); 
void unifR2sided(double *x, double *y, cov_model *cov, double *v); 
void kappa_unif(int i, cov_model *cov, int *nr, int *nc);
int check_unif(cov_model *cov);
void do_unif(cov_model *cov, double *v);
int init_unif(cov_model *cov, gen_storage *s);
void range_unif(cov_model *cov, range_type *range);

#endif /* Families_H*/
