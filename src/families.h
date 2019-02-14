

/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 -- 2017 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/



#ifndef Families_H
#define Families_H 1

#define ARCSQRT_SCALE 0
void arcsqrtD(double *x, model *cov, double *v);
void arcsqrtDlog(double *x, model *cov, double *v);
void arcsqrtDinverse(double *v, model *cov, double *left, double *right);
void arcsqrtP(double *x, model *cov, double *v);
void arcsqrtQ(double *x, model *cov, double *v);
void arcsqrtR(double *x, model *cov, double *v);
int check_arcsqrt_distr(model *cov);
int init_arcsqrt(model *cov, gen_storage VARIABLE_IS_NOT_USED *s);
void do_arcsqrt(model *cov, double *v);
void range_arcsqrt(model *cov, range_type* range);


void determD(double *x, model *cov, double *v);
void determDlog(double *x, model *cov, double *v);
void determDinverse(double *v, model *cov, double *left, double *right);
void determP(double *x, model *cov, double *v);
void determP2sided(double *x, double *y, model *cov, double *v);
void determQ(double *x, model *cov, double *v);
void determR(double *x, model *cov, double *v);
void determR2sided(double *x, double *y, model *cov, double *v);
void kappa_determ(int i, model *cov, int *nr, int *nc);
int check_determ(model *cov);
int init_determ(model *cov, gen_storage *s);
void do_determ(model *cov, double *);
void range_determ(model *cov, range_type *range);


#define DISTR_NAME 0
#define DISTR_NROW 1
#define DISTR_NCOL 2
#define DISTR_DX 3
#define DISTR_PX 4
#define DISTR_QX 5
#define DISTR_RX 6
#define DISTR_ENV 7
#define DISTR_LAST DISTR_ENV
void distrD(double *x, model *cov, double *v);
void distrDlog(double *x, model *cov, double *v);
void distrDinverse(double *v, model *cov, double *left, double *right);
void distrP(double *x, model *cov, double *v);
void distrP2sided(double *x, double *y, model *cov, double *v);
void distrQ(double *x, model *cov, double *v);
void distrR(double *x, model *cov, double *v);
void distrR2sided(double *x, double *y, model *cov, double *v);
int check_distr(model *cov);
int init_distr(model *cov, gen_storage *s);
void do_distr_do(model *cov, double *);
void range_distr(model *cov, range_type *range);



#define SPHERIC_SPACEDIM 0
#define SPHERIC_BALLDIM 1
void sphericD(double *x, model *cov, double *v);
void sphericDlog(double *x, model *cov, double *v);
void sphericDinverse(double *v, model *cov, double *left, double *right);
void sphericP(double *x, model *cov, double *v);
void sphericQ(double *x, model *cov, double *v);
void sphericR(double *x, model *cov, double *v);
int check_RRspheric(model *cov);
int init_RRspheric(model *cov, gen_storage *s);
void do_RRspheric(model *cov, double *);
void range_RRspheric(model *cov, range_type *range);


#define GAUSS_DISTR_MEAN 0
#define GAUSS_DISTR_SD 1
#define GAUSS_DISTR_LOG 2
void gaussD(double *x, model *cov, double *v);
void gaussDlog(double *x, model *cov, double *v);
void gaussDinverse(double *v, model *cov, double *left, double *right);
void gaussP(double *x, model *cov, double *v);
void gaussP2sided(double *x, double *y, model *cov, double *v);
void gaussQ(double *x, model *cov, double *v);
void gaussR(double *x, model *cov, double *v);
void gaussR2sided(double *x, double *y, model *cov, double *v);
void kappa_gauss_distr(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc);
int check_gauss_distr(model *cov);
int init_gauss_distr(model *cov, gen_storage *s);
void do_gauss_distr(model *cov, double *v);
void range_gauss_distr(model *cov, range_type *range);


#define LOC_LOC 0
#define LOC_SCALE 1
#define LOC_POWER 2
void locD(double *x, model *cov, double *v);
void locDlog(double *x, model *cov, double *v);
void locDinverse(double *v, model *cov, double *left, double *right);
void locP(double *x, model *cov, double *v);
void locP2sided(double *x, double *y, model *cov, double *v);
void locQ(double *x, model *cov, double *v);
void locR(double *x, model *cov, double *v);
void locR2sided(double *x, double *y, model *cov, double *v);
void kappa_loc(int i, model *cov, int *nr, int *nc);
int check_loc(model *cov);
int init_loc(model *cov, gen_storage *s);
void do_loc(model *cov, double *v);
void range_loc(model *cov, range_type *range);

#define RECT_SAFETY 0
#define RECT_MINSTEPLENGTH 1
#define RECT_MAXSTEPS 2
#define RECT_PARTS 3
#define RECT_MAXIT 4
#define RECT_INNERMIN 5
#define RECT_OUTERMAX 6
#define RECT_MCMC_N 7
#define RECT_NORMED 8 // unsused
#define RECT_APPROX 9
#define RECT_ONESIDED 10
void rectangularD(double *x, model *cov, double *v);
void rectangularDlog(double *x, model *cov, double *v);
void rectangularDinverse(double *v, model *cov, double *left,double *right);
void rectangularP(double *x, model *cov, double *v);
void rectangularP2sided(double *x, double *y, model *cov, double *v);
void rectangularQ(double *x, model *cov, double *v);
void rectangularR(double *x, model *cov, double *v); 
void rectangularR2sided(double *x, double *y, model *cov, double *v); 
//void kappa_rectangular(int i, model *cov, int *nr, int *nc);
int check_rectangular(model *cov);
void do_rectangular(model *cov, double *v);
int init_rectangular(model *cov, gen_storage *s);
void range_rectangular(model *cov, range_type *range);


#define SETPARAM_LOCAL 0
#define SET_PERFORMDO 0
//#define SETPARAM_VARIANT 1
//#define SETPARAM_FROM 1
//#define SETPARAM_SET 2
void setParamD(double *x, model *cov, double *v); 
void setParamDlog(double *x, model *cov, double *v); 
void setParamDinverse(double *v, model *cov, double *left, double *right);
void setParamP(double *x, model *cov, double *v);   
void setParamP2sided(double *x, double *y, model *cov, double *v);   
void setParamQ(double *x, model *cov, double *v);   
void setParamR(double *x, model *cov, double *v);   
void setParamR2sided(double *x, double *y, model *cov, double *v);   
int check_setParam(model *cov);
void range_setParam(model *cov, range_type *range);
int init_setParam(model *cov, gen_storage *s);
void do_setParam(model *cov, double *v);


#define UNIF_MIN 0
#define UNIF_MAX 1
#define UNIF_NORMED 2
void unifD(double *x, model *cov, double *v);
void unifDlog(double *x, model *cov, double *v);
void unifDinverse(double *v, model *cov, double *left, double *right);
void unifP(double *x, model *cov, double *v);
void unifP2sided(double *x, double *y, model *cov, double *v);
void unifQ(double *x, model *cov, double *v);
void unifR(double *x, model *cov, double *v); 
void unifR2sided(double *x, double *y, model *cov, double *v); 
void kappa_unif(int i, model *cov, int *nr, int *nc);
int check_unif(model *cov);
void do_unif(model *cov, double *v);
int init_unif(model *cov, gen_storage *s);
void range_unif(model *cov, range_type *range);


void kappa_mcmc(int i, model *cov, int *nr, int *nc);
void mcmcD(double *x, model *cov, double *v);
void mcmcDlog(double *x, model *cov, double *v);
void mcmcDinverse(double *v, model *cov, double *left,double *right);
void mcmcP(double *x, model *cov, double *v);
void mcmcP2sided(double *x, double *y, model *cov, double *v);
void mcmcQ(double *x, model *cov, double *v);
void mcmcR(double *x, model *cov, double *v); 
void mcmcR2sided(double *x, double *y, model *cov, double *v); 
//void kappa_mcmc(int i, model *cov, int *nr, int *nc);
int check_mcmc(model *cov);
void do_mcmc(model *cov, double *v);
int init_mcmc(model *cov, gen_storage *s);
void range_mcmc(model *cov, range_type *range);

void randomSign(double *x, model *cov, double *v);
void lograndomSign(double *x, model *cov, double *v, double *Sign); 
void randomSignInverse(double *x, model *cov, double *v);
void randomSignNonstatInverse(double *v, model *cov, double *x, double *y);
int check_randomSign(model *cov);
void range_randomSign(model *cov, range_type *range);
int init_randomSign(model *cov, gen_storage *s);
int struct_randomSign(model *cov, model **newmodel);
void do_randomSign(model *cov, gen_storage *s);



#endif /* Families_H*/
