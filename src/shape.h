

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



#ifndef RandomShape_H
#define RandomShape_H 1

#define IS_IS 1

//-----------------------------------------------------------------
// MPP operator 


void kappamppplus(int i, model *cov, int *nr, int *nc);
int checkmppplus(model *cov);
void rangempplus(model *cov, range_type *range);
int struct_mppplus(model *plus, model **newmodel);
int init_mppplus(model *cov, gen_storage *s);
void do_mppplus(model *cov, gen_storage *s);
void mppplus(double *x, model *cov, double *v); 
//void logmppplus(double *x, model *cov, double *v, double *Sign); 


//-----------------------------------------------------------------
//    Huetchen.cc 

#define PGS_FCT 0 // never change
#define PGS_LOC 1

#define ZHOU_RATIO 0
#define ZHOU_FLATHULL 1
#define ZHOU_INFTY_SMALL 2
#define ZHOU_NORMED 3
#define ZHOU_ISOTROPIC 4
void mcmc_pgs(double *x, model *cov, double *v); 
void logmcmc_pgs(double *x, model *cov, double *v, double *Sign); 
int check_mcmc_pgs(model *cov);
int struct_mcmc_pgs(model *cov, model **newmodel);
int init_mcmc_pgs(model *cov, gen_storage *S);  
void do_mcmc_pgs(model *cov, gen_storage *S);
void range_mcmc_pgs(model *cov, range_type *range);

void Zhou(double *x, model *cov, double *v); 
void logZhou(double *x, model *cov, double *v, double *Sign); 
int check_Zhou(model *cov);
int struct_Zhou(model *cov, model **newmodel);
int init_Zhou(model *cov, gen_storage *S);  
void do_Zhou(model *cov, gen_storage *S);
void range_Zhou(model *cov, range_type *range);

void Ballani(double *x, model *cov, double *v); 
void logBallani(double *x, model *cov, double *v, double *Sign); 
int check_Ballani(model *cov);
int struct_Ballani(model *cov, model **newmodel);
int init_Ballani(model *cov, gen_storage *S);  
void do_Ballani(model *cov, gen_storage *S);
void range_Ballani(model *cov, range_type *range);

void standard_shape(double *x, model *cov, double *v); 
void logstandard_shape(double *x, model *cov, double *v, double *Sign); 
int check_standard_shape(model *cov);
int struct_standard_shape(model *cov, model **newmodel);
int init_standard_shape(model *cov, gen_storage *S);  
void do_standard_shape(model *cov, gen_storage *S);

void stationary_shape(double *x, model *cov, double *v); 
void logstationary_shape(double *x, model *cov, double *v, double *Sign); 
int check_stationary_shape(model *cov);
int struct_stationary_shape(model *cov, model **newmodel);
int init_stationary_shape(model *cov, gen_storage *S);  
void do_stationary_shape(model *cov, gen_storage *S);


//-----------------------------------------------------------------
// trend.cc
#define TREND_MEAN 0
void trend(double *x, model *cov, double *v);
void trend_nonstat(double *x, double *y, model *cov, double *v);
bool allowedItrend(model *cov);
bool settrend(model *cov);
int checktrend(model *cov);
void rangetrend(model *cov, range_type* ra);
int init_trend(model *cov, gen_storage *s);
void do_trend(model *cov, gen_storage *s);
void GetInternalMean(model *cov, int vdim, double *mean);
void kappatrend(int i, model *cov, int *nr, int *nc);
Types Typetrend(Types required, model *cov, isotropy_type required_iso);


#define COVARIATE_C 0
#define COVARIATE_X 1
#define COVARIATE_RAW 2
#define COVARIATE_ADDNA 3
#define COVARIATE_FACTOR 4
void covariate(double *x, model *cov, double *v);
int checkcovariate(model *cov);
void rangecovariate(model VARIABLE_IS_NOT_USED *cov, range_type *range);
void kappa_covariate(int i, model *cov, int *nr, int *nc);


#define FIXCOV_M COVARIATE_C // never change
#define FIXCOV_X COVARIATE_X  // never change
#define FIXCOV_RAW COVARIATE_RAW  // never change
void fixStat(double *x, model *cov, double *v);
void fix(double *x, double *y, model *cov, double *v);
int checkfix(model *cov);
void rangefix(model VARIABLE_IS_NOT_USED *cov, range_type *range);
void kappa_fix(int i, model *cov, int *nr, int *nc);
bool setfix(model *cov);
Types Typefix(Types required, model *cov, isotropy_type requ_iso);
bool allowedIfix(model *cov);
//  bool allowedDfix(model *cov);



//-----------------------------------------------------------------
// shapes

#define BALL_RADIUS 1.0 // nicht mehr aendern!!
void ball(double *x, model *cov, double *v);
int init_ball(model *cov, gen_storage *s);
void do_ball(model *cov, gen_storage *s); 
int struct_ball(model *cov, model **newmodel);
void Inverseball(double *x, model *cov, double *v);


void kappa_EAxxA(int i, model *cov, int *nr, int *nc);
void EAxxA(double *x, model *cov, double *v);
int checkEAxxA(model *cov);
void rangeEAxxA(model *cov, range_type* ra);
void minmaxEigenEAxxA(model *cov, double *mm);

void kappa_EtAxxA(int i, model *cov, int *nr, int *nc);
void EtAxxA(double *x, model *cov, double *v);
int checkEtAxxA(model *cov);
void rangeEtAxxA(model *cov, range_type* ra);
void minmaxEigenEtAxxA(model *cov, double *mm);


#define POLYGON_BETA 0
#define POLYGON_SAFETY 1
void Polygon(double *x, model *cov, double *v);
int check_polygon(model *cov);
void range_polygon(model *cov, range_type *range);
int init_polygon(model *cov, gen_storage *s);
int struct_polygon(model *cov, model **newmodel);
void Inversepolygon(double *x, model *cov, double *v);
void InversepolygonNonstat(double *v, model *cov, double *x, double *y);

void kappa_rotat(int i, model *cov, int *nr, int *nc);
void rotat(double *x, model *cov, double *v);
int checkrotat(model *cov);
void rangerotat(model *cov, range_type* ra);
void minmaxEigenrotat(model *cov, double *mm);

void kappa_Rotat(int i, model *cov, int *nr, int *nc);
void Rotat(double *x, model *cov, double *v);
int checkRotat(model *cov);
void rangeRotat(model *cov, range_type* ra);


void kappa_rational(int i, model *cov, int *nr, int *nc);
void minmaxEigenrational(model *cov, double *mm);
void rational(double *x, model *cov, double *v);
int checkrational(model *cov);
void rangerational(model *cov, range_type* ra);


#define TRUNC_RADIUS 0
void truncsupport(double *x, model *cov, double *v);
int checktruncsupport(model *cov);
void truncsupportInverse(double *x, model *cov, double *v);
void rangetruncsupport(model *cov, range_type *range);
int struct_truncsupport(model *cov, model **newmodel);
int init_truncsupport(model *cov, gen_storage *s);
void do_truncsupport(model *cov, gen_storage *s);



//-----------------------------------------------------------------
// PointShape
int FillInPts(model *key, model *shape);




#endif /* RandomShape_H */
 
