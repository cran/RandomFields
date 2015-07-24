

/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 Martin Schlather

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

//-----------------------------------------------------------------
// MPP operator


void kappamppplus(int i, cov_model *cov, int *nr, int *nc);
int checkmppplus(cov_model *cov);
void rangempplus(cov_model *cov, range_type *range);
int struct_mppplus(cov_model *plus, cov_model **newmodel);
int init_mppplus(cov_model *cov, gen_storage *s);
void do_mppplus(cov_model *cov, gen_storage *s);
void mppplus(double *x, cov_model *cov, double *v); 
//void logmppplus(double *x, cov_model *cov, double *v, double *Sign); 


void mcmc_pgs(double *x, cov_model *cov, double *v); 
void logmcmc_pgs(double *x, cov_model *cov, double *v, double *Sign); 
int check_mcmc_pgs(cov_model *cov);
int struct_mcmc_pgs(cov_model *cov, cov_model **newmodel);
int init_mcmc_pgs(cov_model *cov, gen_storage *S);  
void do_mcmc_pgs(cov_model *cov, gen_storage *S);
void range_mcmc_pgs(cov_model *cov, range_type *range);

void pts_given_shape(double *x, cov_model *cov, double *v); 
void logpts_given_shape(double *x, cov_model *cov, double *v, double *Sign); 
int check_pts_given_shape(cov_model *cov);
int struct_pts_given_shape(cov_model *cov, cov_model **newmodel);
int init_pts_given_shape(cov_model *cov, gen_storage *S);  
void do_pts_given_shape(cov_model *cov, gen_storage *S);
void range_pts_given_shape(cov_model *cov, range_type *range);

void standard_shape(double *x, cov_model *cov, double *v); 
void logstandard_shape(double *x, cov_model *cov, double *v, double *Sign); 
int check_standard_shape(cov_model *cov);
int struct_standard_shape(cov_model *cov, cov_model **newmodel);
int init_standard_shape(cov_model *cov, gen_storage *S);  
void do_standard_shape(cov_model *cov, gen_storage *S);

void stationary_shape(double *x, cov_model *cov, double *v); 
void logstationary_shape(double *x, cov_model *cov, double *v, double *Sign); 
int check_stationary_shape(cov_model *cov);
int struct_stationary_shape(cov_model *cov, cov_model **newmodel);
int init_stationary_shape(cov_model *cov, gen_storage *S);  
void do_stationary_shape(cov_model *cov, gen_storage *S);


//-----------------------------------------------------------------
// Trends
void mixed(double *x,  cov_model *cov, double *v);
void MLEmixed(double *x,  cov_model *cov, double *v);
void mixed_nonstat(double *x,  double *y, cov_model *cov, double *v);
void MLEmixed_nonstat(double *x,  double *y, cov_model *cov, double *v);
void kappamixed(int i, cov_model *cov, int *nr, int *nc);
int  checkmixed(cov_model *cov);
void rangemixed(cov_model *cov, range_type* ra);
int  initmixed(cov_model *cov, gen_storage *s);
void domixed(cov_model *cov, gen_storage *s);
void covmatrix_mixed(cov_model *cov, double *v);
char iscovmatrix_mixed(cov_model *cov);

void trend(double *x, cov_model *cov, double *v);
void trend_nonstat(double *x, double *y, cov_model *cov, double *v);
int checktrend(cov_model *cov);
int checktrendproc(cov_model *cov);
void rangetrend(cov_model *cov, range_type* ra);
int init_trend(cov_model *cov, gen_storage *s);
void do_trend(cov_model *cov, gen_storage *s);
void GetInternalMean(cov_model *cov, int vdim, double *mean);
void kappatrend(int i, cov_model *cov, int *nr, int *nc);


//-----------------------------------------------------------------
// shapes

void kappa_ave(int i, cov_model *cov, int *nr, int *nc);
void ave(double *x, cov_model *cov, double *v);
int checkave(cov_model *cov);
void rangeave(cov_model *cov, range_type* ra);
//void sd_standard(mpp_storage *s, cov_model *cov);
int structAve(cov_model *cov, cov_model **newmodel);

int check_shapeave(cov_model *cov);
int init_shapeave(cov_model *cov, gen_storage *s);
void do_shapeave(cov_model *cov, gen_storage *s);
void logshapeave(double *x, cov_model *cov, double *v, double *Sign);

void ball(double *x, cov_model *cov, double *v);
int init_ball(cov_model *cov, gen_storage *s);
void do_ball(cov_model *cov, gen_storage *s); 
int struct_ball(cov_model *cov, cov_model **newmodel);
void Inverseball(double *x, cov_model *cov, double *v);

void kappa_EAxxA(int i, cov_model *cov, int *nr, int *nc);
void EAxxA(double *x, cov_model *cov, double *v);
int checkEAxxA(cov_model *cov);
void rangeEAxxA(cov_model *cov, range_type* ra);
void minmaxEigenEAxxA(cov_model *cov, double *mm);

void kappa_EtAxxA(int i, cov_model *cov, int *nr, int *nc);
void EtAxxA(double *x, cov_model *cov, double *v);
int checkEtAxxA(cov_model *cov);
void rangeEtAxxA(cov_model *cov, range_type* ra);
void minmaxEigenEtAxxA(cov_model *cov, double *mm);

void Polygon(double *x, cov_model *cov, double *v);
int check_polygon(cov_model *cov);
void range_polygon(cov_model *cov, range_type *range);
int init_polygon(cov_model *cov, gen_storage *s);
int struct_polygon(cov_model *cov, cov_model **newmodel);
void Inversepolygon(double *x, cov_model *cov, double *v);
void InversepolygonNonstat(double *v, cov_model *cov, double *x, double *y);

void kappa_rotat(int i, cov_model *cov, int *nr, int *nc);
void rotat(double *x, cov_model *cov, double *v);
int checkrotat(cov_model *cov);
void rangerotat(cov_model *cov, range_type* ra);
void minmaxEigenrotat(cov_model *cov, double *mm);

void kappa_Rotat(int i, cov_model *cov, int *nr, int *nc);
void Rotat(double *x, cov_model *cov, double *v);
int checkRotat(cov_model *cov);
void rangeRotat(cov_model *cov, range_type* ra);


void kappa_rational(int i, cov_model *cov, int *nr, int *nc);
void minmaxEigenrational(cov_model *cov, double *mm);
void rational(double *x, cov_model *cov, double *v);
int checkrational(cov_model *cov);
void rangerational(cov_model *cov, range_type* ra);

void kappa_stp(int i, cov_model *cov, int *nr, int *nc);
void stp(double *x,  double *y, cov_model *cov, double *v);
int checkstp(cov_model *cov);
void rangestp(cov_model *cov, range_type* ra);
int structStp(cov_model *cov, cov_model **newmodel);

int check_shapestp(cov_model *cov);
int init_shapestp(cov_model *cov, gen_storage *s);
void do_shapestp(cov_model *cov, gen_storage *s);
void logshapestp(double *x, double *u, cov_model *cov, double *v, double *);

void truncsupport(double *x, cov_model *cov, double *v);
int checktruncsupport(cov_model *cov);
void truncsupportInverse(double *x, cov_model *cov, double *v);
void rangetruncsupport(cov_model *cov, range_type *range);
int struct_truncsupport(cov_model *cov, cov_model **newmodel);
int init_truncsupport(cov_model *cov, gen_storage *s);
void do_truncsupport(cov_model *cov, gen_storage *s);

void randomSign(double *x, cov_model *cov, double *v);
void lograndomSign(double *x, cov_model *cov, double *v, double *Sign); 
void randomSignInverse(double *x, cov_model *cov, double *v);
void randomSignNonstatInverse(double *v, cov_model *cov, double *x, double *y);
int check_randomSign(cov_model *cov);
void range_randomSign(cov_model *cov, range_type *range);
int init_randomSign(cov_model *cov, gen_storage *s);
int struct_randomSign(cov_model *cov, cov_model **newmodel);
void do_randomSign(cov_model *cov, gen_storage *s);





//-----------------------------------------------------------------
// Gauss Methods
void kappa_ce(int i, cov_model *cov, int *nr, int *nc);
int check_ce(cov_model *cov); 
void range_ce(cov_model *cov, range_type *range);  
int init_circ_embed(cov_model *cov, gen_storage *s);
void do_circ_embed(cov_model *cov, gen_storage *s);

int struct_ce_approx(cov_model *cov, cov_model **newmodel);
int init_ce_approx(cov_model *cov, gen_storage *S);
void do_ce_approx(cov_model *cov, gen_storage *S);

void kappa_localproc(int i, cov_model *cov, int *nr, int *nc);

int check_co_proc(cov_model *cov);
int check_co_intern(cov_model *cov);
void range_co_proc(cov_model *cov, range_type* ra);

void distrD(double *x, cov_model *cov, double *v);
void distrP(double *x, cov_model *cov, double *v);
void distrQ(double *x, cov_model *cov, double *v);
void distrR(double *x, cov_model *cov, double *v);
void kappa_distr(int i, cov_model *cov, int *nr, int *nc);
int check_distr(cov_model *cov);
int init_distr(cov_model *cov, gen_storage *s);
void do_distr(cov_model *cov, gen_storage *s);
void range_distr(cov_model *cov, range_type *range);
void range_intrinCE(cov_model *cov, range_type* ra);
int check_local_proc(cov_model *cov);

int check_directGauss(cov_model *cov);
void range_direct(cov_model *cov, range_type *range);
int init_directGauss(cov_model *cov, gen_storage *s);
void do_directGauss(cov_model *cov, gen_storage *s);


int check_hyperplane(cov_model *cov);
int check_hyperplane_intern(cov_model *cov);
void range_hyperplane(cov_model *cov, range_type *range);
int struct_hyperplane(cov_model *cov, cov_model **newmodel);
int init_hyperplane(cov_model *cov, gen_storage *s);
void do_hyperplane(cov_model *cov, gen_storage *s);


int check_nugget_proc(cov_model *cov);
void range_nugget_proc(cov_model *cov, range_type *range);
int init_nugget(cov_model *cov, gen_storage *s);
void do_nugget(cov_model *cov, gen_storage *s );
int struct_nugget(cov_model *cov, cov_model **newmodel);


int check_randomcoin(cov_model *cov);
void range_randomcoin(cov_model *cov, range_type *range);
//void coin(double *x, cov_model *cov, double *v);
int init_randomcoin(cov_model *cov, gen_storage *s);
int struct_randomcoin(cov_model *cov, cov_model **newmodel);
//void coinInverse(double *x, cov_model *cov, double *v);
int struct_randomcoin_extract(cov_model *cov, cov_model **newmodel);



int check_sequential(cov_model *cov) ;
void range_sequential(cov_model *cov, range_type *range) ;
int init_sequential(cov_model *cov, gen_storage *s);
void do_sequential(cov_model *cov, gen_storage *s);


int check_spectral(cov_model *cov) ;
void range_spectral(cov_model *cov, range_type *range) ;
int init_spectral(cov_model *cov, gen_storage *s);
void do_spectral(cov_model *cov, gen_storage *s );
int struct_spectral(cov_model *plus, cov_model **newmodel);


int check_specificGauss(cov_model *cov);
int struct_specificGauss(cov_model *cov, cov_model **newmodel);
int init_specificGauss(cov_model *cov, gen_storage *S);
void do_specificGauss(cov_model *cov, gen_storage *S);
void range_specificGauss(cov_model *cov, range_type *range);


void tbm_kappasproc(int i, cov_model *cov, int *nr, int *nc);
int checktbmproc(cov_model *cov);
void rangetbmproc(cov_model *cov, range_type* ra);
int struct_tbmproc(cov_model *cov, cov_model **newmodel);
int init_tbmproc(cov_model *cov, gen_storage *s);
void do_tbmproc(cov_model *cov, gen_storage *s);
//void tbm2num(double *x, cov_model *cov, double *v);


//-----------------------------------------------------------------
// MPP Methods

int init_BRorig(cov_model *cov, gen_storage *s);
void do_BRorig(cov_model *cov, gen_storage *s);

int init_BRshifted(cov_model *cov, gen_storage *s);
void do_BRshifted(cov_model *cov, gen_storage *s);

int check_BRmixed(cov_model *cov);
int init_BRmixed(cov_model *cov, gen_storage *s);
void do_BRmixed(cov_model *cov, gen_storage *s);
void range_BRmixed(cov_model *cov, range_type *range);
void kappaBRmixed(int i, cov_model *cov, int *nr, int *nc);
 
int initBRuser (cov_model *cov, gen_storage *S);
int structBRuser(cov_model *cov, cov_model **newmodel);
int structBRintern(cov_model *cov, cov_model **newmodel);



//-----------------------------------------------------------------
// Processes

// checkSproc by checkS
int structSproc(cov_model *cov, cov_model **newmodel);
int initSproc(cov_model *cov, gen_storage *s);
void doSproc(cov_model *cov, gen_storage *s);

int checkplusproc(cov_model *cov);
 int structplusproc(cov_model *cov, cov_model **newmodel);
int initplusproc(cov_model *cov, gen_storage *s);
void doplusproc(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) ;
 
int checkmultproc(cov_model *cov);
int structmultproc(cov_model *cov, cov_model **newmodel);
void rangemultproc(cov_model *cov, range_type *range);
int initmultproc(cov_model *cov, gen_storage *s);
void domultproc(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s) ;
 
void kappa_binaryprocess(int i, cov_model *cov, int *nr, int *nc);
void binary(double *x, cov_model *cov, double *v);
int checkbinaryprocess(cov_model *cov);
void rangebinaryprocess(cov_model *cov, range_type *range);
int struct_binaryprocess(cov_model *cov, cov_model **newmodel);
int init_binaryprocess(cov_model *cov, gen_storage *s);
void do_binaryprocess(cov_model *cov, gen_storage *s);

void BrownResnick(double *x, cov_model *cov, double *v);
int checkBrownResnickProc(cov_model *cov);
int structBrownResnick(cov_model *cov, cov_model **newmodel);
int initBrownResnick(cov_model *cov, gen_storage *s);
void doBrownResnick(cov_model *cov, gen_storage *s);
void loglikelihoodBR(double *data, cov_model *cov, double *v);

void kappaGProc(int i, cov_model *cov, int *nr, int *nc);
int checkgaussprocess(cov_model *cov);
void rangegaussprocess(cov_model *cov, range_type *range);
int struct_gaussprocess(cov_model *cov, cov_model **newmodel);
int init_gaussprocess(cov_model *cov, gen_storage *s);
void do_gaussprocess(cov_model *cov, gen_storage *s);
void gaussprocessDlog(double *x, cov_model *cov, double *v);
int struct_gauss_logli(cov_model *cov);
SEXP gauss_linearpart(SEXP model_reg, SEXP Set);
void gauss_predict(cov_model *predict, cov_model *Cov, double *v);
void gauss_trend(cov_model *predict, cov_model *cov, double *v, int set);


int struct_extractdollar(cov_model *cov, cov_model **newmodel);

int check_poisson(cov_model *cov);
int struct_poisson(cov_model *cov, cov_model **newmodel);
void range_poisson(cov_model *cov, range_type *range);
int init_poisson(cov_model *cov, gen_storage *S);

void extremalgaussian(double *x, cov_model *cov, double *v);
int check_schlather(cov_model *cov);
int struct_schlather(cov_model *cov, cov_model **newmodel);
void loglikelihoodSchlather(double *data, cov_model *cov, double *v);

int struct_smith(cov_model *cov, cov_model **newmodel);
int check_smith(cov_model *cov);
int struct_smith_pts(cov_model **Key, cov_model *shape, cov_model *calling,
		     int tsdim, int vdim);


int checkchisqprocess(cov_model *cov) ;
void rangechisqprocess(cov_model *cov, range_type *range) ;
int struct_chisqprocess(cov_model *cov, cov_model **newmodel) ;
int init_chisqprocess(cov_model *cov, gen_storage *s) ;
void do_chisqprocess(cov_model *cov, gen_storage *s);

 
void rangetprocess(cov_model *cov, range_type *range) ;
void do_tprocess(cov_model *cov, gen_storage *s);

int init_opitzprocess(cov_model *cov, gen_storage *s);
void range_opitz(cov_model *cov, range_type *range);

void dompp(cov_model *cov, gen_storage *S); 
int init_mpp(cov_model *cov, gen_storage *S);
void range_mpp(cov_model *cov, range_type *range);
int SetGEVetc(cov_model *cov, int type);


//-----------------------------------------------------------------
// Interfaces
void simulate(double *x, cov_model *cov, double *v);
int check_simulate(cov_model *cov); 
void range_simulate(cov_model *cov, range_type *range);  
int struct_simulate(cov_model *cov, cov_model **newmodel);
int init_simulate(cov_model *cov, gen_storage *S);
// void do_simulate(cov_model *cov, gen_storage *S);
void range_simulate(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range);

void density(double *x, cov_model *cov, double *v);
int check_density(cov_model *cov); 
void range_density(cov_model *cov, range_type *range);  
int struct_density(cov_model *cov, cov_model **newmodel);
int init_density(cov_model *cov, gen_storage *S);
// void do_density(cov_model *cov, gen_storage *S);
void range_density(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range);


void kappalikelihood(int i, cov_model VARIABLE_IS_NOT_USED *cov, 
		     int *nr, int *nc);
void likelihood(double *data, cov_model *cov, double *v);
int check_likelihood(cov_model *cov);
int struct_likelihood(cov_model *cov, cov_model **newmodel);
void range_likelihood(cov_model *cov, range_type* range);

void linearpart(double *data, cov_model *cov, double *v);
int check_linearpart(cov_model *cov);
int struct_linearpart(cov_model *cov, cov_model **newmodel);
//void range_linearpart(cov_model *cov, range_type* range);


void predict(double VARIABLE_IS_NOT_USED *x, cov_model *cov, double *v);
int check_predict(cov_model *predict);
int struct_predict(cov_model *cov, cov_model VARIABLE_IS_NOT_USED  **newmodel);
void range_predict(cov_model VARIABLE_IS_NOT_USED *cov, range_type* range);


int check_cov_intern(cov_model *cov, Types type, bool close, bool kernel);
void Cov(double *x, cov_model *cov, double *value) ;
int check_cov(cov_model *cov) ;
int struct_cov(cov_model *cov, cov_model **newmodel);

void FctnIntern(cov_model *cov, cov_model *covVdim, cov_model *sub, double *value, bool ignore_y);
void Fctn(double *x, cov_model *cov, double *value);
int check_fctn(cov_model *cov);

void CovMatrix(double *x, cov_model *cov, double *value);
int check_covmatrix(cov_model *cov) ;

void EvalDistr(double *x, cov_model *cov, double *v);
void kappa_EvalDistr(int i, cov_model *cov, int *nr, int *nc);
int check_EvalDistr(cov_model *cov); 
void range_EvalDistr(cov_model *cov, range_type *range);  
int struct_EvalDistr(cov_model *cov, cov_model **newmodel);
int init_EvalDistr(cov_model *cov, gen_storage *S);
// void do_EvalDistr(cov_model *cov, gen_storage *S);

void RFget(double *x, cov_model *cov, double *v);
int SearchParam(cov_model *cov, get_storage *s) ;
int check_RFget(cov_model *cov) ;
void range_RFget(cov_model *cov, range_type* range);
int struct_RFget(cov_model *cov, cov_model **newmodel);


void Variogram(double *x, cov_model *cov, double *value) ;
void Pseudovariogram(double *x, cov_model *cov, double *value) ;
int check_vario(cov_model *cov);
int struct_variogram(cov_model *cov, cov_model **newmodel);

void Dummy(double *x, cov_model *cov, double *value);
int check_dummy(cov_model *cov);
int struct_dummy(cov_model *cov, cov_model **newmodel);


int checkTrendEval(cov_model *cov);
int init_TrendEval(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s);
void do_TrendEval(cov_model *cov, gen_storage *s);
void range_TrendEval(cov_model  *cov, range_type *range);

//-----------------------------------------------------------------
// unsorted


#endif /* RandomShape_H */
 
