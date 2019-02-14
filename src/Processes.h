
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



#ifndef RF_PROCESSES_H
#define RF_PROCESSES_H 1


//-----------------------------------------------------------------
// Gauss Methods
#define GAUSS_BOXCOX 0
#define COMMON_GAUSS 0


#define CE_FORCE (COMMON_GAUSS + 1)
#define CE_MMIN  (COMMON_GAUSS + 2)
#define CE_STRATEGY  (COMMON_GAUSS + 3)
#define CE_MAXGB  (COMMON_GAUSS + 4)
#define CE_MAXMEM  (COMMON_GAUSS + 5)
#define CE_TOLIM  (COMMON_GAUSS + 6)
#define CE_TOLRE  (COMMON_GAUSS + 7)
#define CE_TRIALS (COMMON_GAUSS +  8)
#define CE_USEPRIMES  (COMMON_GAUSS + 9)
#define CE_DEPENDENT  (COMMON_GAUSS + 10)
#define CE_APPROXSTEP  (COMMON_GAUSS + 11)
#define CE_APPROXMAXGRID  (COMMON_GAUSS + 12) // alway last of CE_*
#define CE_LAST CE_APPROXMAXGRID 
void kappa_ce(int i, model *cov, int *nr, int *nc);
int check_ce(model *cov); 
void range_ce(model *cov, range_type *range);  
int init_circ_embed(model *cov, gen_storage *s);
void do_circ_embed(model *cov, gen_storage *s);


void Xkappa_ce(int i, model *cov, int *nr, int *nc);
int Xcheck_ce(model *cov); 
void Xrange_ce(model *cov, range_type *range);  
int Xinit_circ_embed(model *cov, gen_storage *s);
void Xdo_circ_embed(model *cov, gen_storage *s);
int Xstruct_ce_approx(model *cov, model **newmodel);
int Xinit_ce_approx(model *cov, gen_storage *S);
void Xdo_ce_approx(model *cov, gen_storage *S);
void Xkappa_localproc(int i, model *cov, int *nr, int *nc);
int Xcheck_co_proc(model *cov);
int Xcheck_co_intern(model *cov);
void Xrange_co_proc(model *cov, range_type* ra);


int struct_ce_approx(model *cov, model **newmodel);
int init_ce_approx(model *cov, gen_storage *S);
void do_ce_approx(model *cov, gen_storage *S);

void kappa_localproc(int i, model *cov, int *nr, int *nc);


int check_co_proc(model *cov);
int check_co_intern(model *cov);
void range_co_proc(model *cov, range_type* ra);

void distrD(double *x, model *cov, double *v);
void distrP(double *x, model *cov, double *v);
void distrQ(double *x, model *cov, double *v);
void distrR(double *x, model *cov, double *v);
void kappa_distr(int i, model *cov, int *nr, int *nc);
int check_distr(model *cov);
int init_distr(model *cov, gen_storage *s);
void do_distr(model *cov, gen_storage *s);
void range_distr(model *cov, range_type *range);
void range_intrinCE(model *cov, range_type* ra);
int check_local_proc(model *cov);

int check_directGauss(model *cov);
void range_direct(model *cov, range_type *range);
int init_directGauss(model *cov, gen_storage *s);
void do_directGauss(model *cov, gen_storage *s);


int check_hyperplane(model *cov);
int check_hyperplane_intern(model *cov);
void range_hyperplane(model *cov, range_type *range);
int struct_hyperplane(model *cov, model **newmodel);
int init_hyperplane(model *cov, gen_storage *s);
void do_hyperplane(model *cov, gen_storage *s);


int check_nugget_proc(model *cov);
void range_nugget_proc(model *cov, range_type *range);
int init_nugget(model *cov, gen_storage *s);
void do_nugget(model *cov, gen_storage *s );
int struct_nugget(model *cov, model **newmodel);


#define COIN_COV 0
#define COIN_SHAPE 1
int check_randomcoin(model *cov);
void range_randomcoin(model *cov, range_type *range);
//void coin(double *x, model *cov, double *v);
int init_randomcoin(model *cov, gen_storage *s);
int struct_randomcoin(model *cov, model **newmodel);
//void coinInverse(double *x, model *cov, double *v);
int struct_randomcoin_extract(model *cov, model **newmodel);


int checkNormedBrProc(model *cov);
int structNormedBrProc(model *cov, model **newmodel);
int initNormedBrProc(model *cov, gen_storage *s);
 void doNormedBrProc(model *cov, gen_storage *s);



int check_sequential(model *cov) ;
void range_sequential(model *cov, range_type *range) ;
int init_sequential(model *cov, gen_storage *s);
void do_sequential(model *cov, gen_storage *s);


int check_spectral(model *cov) ;
void range_spectral(model *cov, range_type *range) ;
int init_spectral(model *cov, gen_storage *s);
void do_spectral(model *cov, gen_storage *s );
int struct_spectral(model *plus, model **newmodel);


int check_specificGauss(model *cov);
int struct_specificGauss(model *cov, model **newmodel);
int init_specificGauss(model *cov, gen_storage *S);
void do_specificGauss(model *cov, gen_storage *S);
void range_specificGauss(model *cov, range_type *range);


#define TBM_FULLDIM (COMMON_GAUSS + 1)
#define TBM_TBMDIM (COMMON_GAUSS + 2)
#define TBM_LAYERS (COMMON_GAUSS + 3)
void tbm_kappasproc(int i, model *cov, int *nr, int *nc);
int checktbmproc(model *cov);
void rangetbmproc(model *cov, range_type* ra);
int struct_tbmproc(model *cov, model **newmodel);
int init_tbmproc(model *cov, gen_storage *s);
void do_tbmproc(model *cov, gen_storage *s);



#define GAUSSPROC_STATONLY  (COMMON_GAUSS + 1)
#define GAUSSPROC_LAST GAUSSPROC_STATONLY
void kappaGProc(int i, model *cov, int *nr, int *nc);
int checkgaussprocess(model *cov);
void rangegaussprocess(model *cov, range_type *range);
int struct_gaussprocess(model *cov, model **newmodel);
int init_gaussprocess(model *cov, gen_storage *s);
void do_gaussprocess(model *cov, gen_storage *s);
void gaussprocessDlog(double *x, model *cov, double *v);
int struct_gauss_logli(model *cov);
SEXP gauss_linearpart(SEXP model_reg, SEXP Set);
void gauss_predict(model *predict, model *Cov, double *v);
void gauss_trend(model *predict, model *cov, double *v, int set);
int struct_gaussmethod(model *cov, model **newmodel);


void boxcox_trafo(double boxcox[], int vdim, double *res, long pts, int repet);	
void boxcox_inverse(double boxcox[], int vdim, double *res, int pts, 
		    int repet);	
#define SAVE_GAUSS_TRAFO    
#define BOXCOX_TRAFO boxcox_trafo(boxcox, res, Gettotalpoints(cov) * VDIM0); 
#define BOXCOX_INVERSE_PARAM(GAUSS_BOXCOX)				\
  boxcox_inverse(P(GAUSS_BOXCOX), VDIM0, res, Gettotalpoints(cov), 1)
#define BOXCOX_INVERSE BOXCOX_INVERSE_PARAM(GAUSS_BOXCOX)
int kappaBoxCoxParam(model *cov, int BC);


#define  GAUSS_COMMON_RANGE_PARAM(GAUSS_BOXCOX)	\
  range->min[GAUSS_BOXCOX] = RF_NEGINF;		\
  range->max[GAUSS_BOXCOX] = RF_INF;		\
  range->pmin[GAUSS_BOXCOX] = 0;		\
  range->pmax[GAUSS_BOXCOX] = 2;		\
  range->openmin[GAUSS_BOXCOX] = false;		\
  range->openmax[GAUSS_BOXCOX] = false
#define  GAUSS_COMMON_RANGE GAUSS_COMMON_RANGE_PARAM(GAUSS_BOXCOX)


#define FRAME_ASSERT_GAUSS assert(hasGaussMethodFrame(cov));
#define FRAME_ASSERT_GAUSS_INTERFACE assert(hasGaussMethodFrame(cov) || (hasInterfaceFrame(cov) && cov->calling != NULL && cov->calling->calling == NULL));

//-----------------------------------------------------------------
// max-stable processes
#define GEV_XI 0
#define GEV_MU 1
#define GEV_S 2
#define LAST_MAXSTABLE GEV_S

#define MPP_SHAPE 0
#define MPP_TCF 1

void extremalgaussian(double *x, model *cov, double *v);
int check_schlather(model *cov);
int struct_schlather(model *cov, model **newmodel);
void loglikelihoodSchlather(double *data, model *cov, double *v);

int struct_smith(model *cov, model **newmodel);
int check_smith(model *cov);
int struct_smith_pts(model **Key, model *shape, model *calling,
		     int logicaldim, int vdim);

void BrownResnick(double *x, model *cov, double *v);
int checkBrownResnickProc(model *cov);
int structBrownResnick(model *cov, model **newmodel);
int initBrownResnick(model *cov, gen_storage *s);
void doBrownResnick(model *cov, gen_storage *s);
void finaldoBrownResnick(model *cov, double *res, int n, gen_storage *s);
void loglikelihoodBR(double *data, model *cov, double *v);


int init_opitzprocess(model *cov, gen_storage *s);
void range_opitz(model *cov, range_type *range);


//-----------------------------------------------------------------
// max-stable methods
int SetGEVetc(model *cov);

int init_BRorig(model *cov, gen_storage *s);
void do_BRorig(model *cov, gen_storage *s);
void finalmaxstable(model *cov, double *res, int n, gen_storage *s);

int init_BRshifted(model *cov, gen_storage *s);
void do_BRshifted(model *cov, gen_storage *s);

int check_BRmixed(model *cov);
int init_BRmixed(model *cov, gen_storage *s);
void do_BRmixed(model *cov, gen_storage *s);
void range_BRmixed(model *cov, range_type *range);
void kappaBRmixed(int i, model *cov, int *nr, int *nc);
 
int initBRuser (model *cov, gen_storage *S);
int structBRuser(model *cov, model **newmodel);
int structBRintern(model *cov, model **newmodel);

void kappabrnormed(int i, model *cov, int *nr, int *nc);
int check_brnormed(model *cov);
int struct_brnormed(model *cov, model **newmodel);
int init_brnormed(model *cov, gen_storage *s);
void do_brnormed(model *cov, gen_storage *s);
void range_brnormed(model  *cov, range_type *range); 


//-----------------------------------------------------------------
// Specific Processes

int structSproc(model *cov, model **newmodel);
int initSproc(model *cov, gen_storage *s);
void doSproc(model *cov, gen_storage *s);

int checkplusproc(model *cov);
 int structplusproc(model *cov, model **newmodel);
int initplusproc(model *cov, gen_storage *s);
void doplusproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) ;
 
int checkmultproc(model *cov);
int structmultproc(model *cov, model **newmodel);
void rangemultproc(model *cov, range_type *range);
int initmultproc(model *cov, gen_storage *s);
void domultproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) ;
 


// other processes
void kappa_binaryprocess(int i, model *cov, int *nr, int *nc);
void binary(double *x, model *cov, double *v);
int checkbinaryprocess(model *cov);
void rangebinaryprocess(model *cov, range_type *range);
int struct_binaryprocess(model *cov, model **newmodel);
int init_binaryprocess(model *cov, gen_storage *s);
void do_binaryprocess(model *cov, gen_storage *s);



int check_poisson(model *cov);
int struct_poisson(model *cov, model **newmodel);
void range_poisson(model *cov, range_type *range);
int init_poisson(model *cov, gen_storage *S);


int checkchisqprocess(model *cov) ;
void rangechisqprocess(model *cov, range_type *range) ;
int struct_chisqprocess(model *cov, model **newmodel) ;
int init_chisqprocess(model *cov, gen_storage *s) ;
void do_chisqprocess(model *cov, gen_storage *s);

 
void rangetprocess(model *cov, range_type *range) ;
void do_tprocess(model *cov, gen_storage *s);



void dompp(model *cov, gen_storage *S); 
int init_mpp(model *cov, gen_storage *S);
void range_mpp(model *cov, range_type *range);


int checktrafoproc(model *cov);
int structtrafoproc(model  *cov, model **newmodel);
int inittrafoproc(model *cov, gen_storage *s);
void dotrafoproc(model *cov, gen_storage *s);

int checkprodproc(model *cov);
int structprodproc(model  *cov, model **newmodel);
int initprodproc(model *cov, gen_storage *s);
void doprodproc(model *cov, gen_storage *s);

int structMproc(model *cov, model **newmodel);
int initMproc(model *cov, gen_storage *s);
void doMproc(model *cov, gen_storage *s);

int checkTrendproc(model *cov);
int init_Trendproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s);
void do_Trendproc(model *cov, gen_storage *s);

int checkvar2covproc(model *cov);
int structvar2covproc(model *cov, model **newmodel);
int initvar2covproc(model *cov, gen_storage *s);
void dovar2covproc(model *cov, gen_storage *s);

  
#endif
