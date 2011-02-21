#ifndef Operators_H
#define Operators_H 1

#include "primitive.h"

void iso2iso(double *x, cov_model *cov, double *v) ;
void spaceDiso2iso(double *x, cov_model *cov, double *v) ;
void spiso2spiso(double *x, cov_model *cov, double *v) ;
void spaceDspiso2spiso(double *x, cov_model *cov, double *v);
void spacetime2iso(double *x, cov_model *cov, double *v) ;
void spaceDspiso2iso(double *x, cov_model *cov, double *v);
void Stat2iso(double *x, cov_model *cov, double *v) ;
void Stat2spacetime(double *x, cov_model *cov, double *v) ;
void Stat2Stat(double *x, cov_model *cov, double *v);
void Nonstat2iso(double *x, double *y, cov_model *cov, double *v) ;
void Nonstat2spacetime(double *x, double *y, cov_model *cov, double *v) ;
void Nonstat2Stat(double *x, double *y, cov_model *cov, double *v) ;
void D_2(double *x, cov_model *cov, double *v);
void DD_2(double *x, cov_model *cov, double *v);
void take2(cov_model *cov);
int check2(cov_model *cov) ;
void range2(cov_model *cov, range_arraytype* ra);
int init2(method_type *meth);
void do2(method_type *meth, res_type *res);
int initspectral2(cov_model *cov) ;
void spectral2(cov_model *cov, spectral_storage *s, double *e) ;
void inv2(double *x, cov_model *cov, double *v);
void mppinit2(mpp_storage *s, cov_model *cov);
void coin2(mpp_storage *s, cov_model *cov);
void sd2(mpp_storage *s, cov_model *cov);
res_type mppgetiso2iso(double *x, cov_model *cov, mpp_storage *s) ;
res_type mppgetspaceDiso2iso(double *x, cov_model *cov, mpp_storage *s) ;
res_type mppgetspiso2spiso(double *x, cov_model *cov, mpp_storage *s) ;
res_type mppgetspaceDspiso2spiso(double *x, cov_model *cov, mpp_storage *s);
res_type mppgetspacetime2iso(double *x, cov_model *cov, mpp_storage *s) ;
res_type mppgetspaceDspiso2iso(double *x, cov_model *cov, mpp_storage *s);
res_type mppgetStat2iso(double *x, cov_model *cov, mpp_storage *s) ;
res_type mppgetStat2spacetime(double *x, cov_model *cov, mpp_storage *s) ;
res_type mppgetStat2Stat(double *x, cov_model *cov, mpp_storage *s);
res_type mppgetNonstat2iso(double *x, double *y, cov_model *cov, mpp_storage *s) ;
res_type mppgetNonstat2spacetime(double *x, double *y, cov_model *cov,
			       mpp_storage *s);
res_type mppgetNonstat2Stat(double *x, double *y, cov_model *cov, mpp_storage *s) ;

void plusStat(double *x, cov_model *cov, double *v);
void MLEplusStat(double *x, cov_model *cov, double *v);
void plusNonStat(double *x,  double *y, cov_model *cov, double *v);
void MLEplusNonStat(double *x,  double *y, cov_model *cov, double *v);
int checkplus(cov_model *cov) ;
void rangeplus(cov_model *cov, range_arraytype* ra);
void Dplus(double *x, cov_model *cov, double *v);
void DDplus(double *x, cov_model *cov, double *v);
int initspectralplus(cov_model *cov) ;
void spectralplus(cov_model *cov, spectral_storage *s, double *e);
int initplus(method_type *meth);
void doplus(method_type *meth, res_type *res);

void malStat(double *x, cov_model *cov, double *v);
void malNonStat(double *x,  double *y, cov_model *cov, double *v);
int checkmal(cov_model *cov) ;
void Dmal(double *x, cov_model *cov, double *v);
void rangemal(cov_model *cov, range_arraytype* ra);
int initmal(method_type *meth);
void domal(method_type *meth, res_type *res);

void kappaM(int i, cov_model *cov, int *nr, int *nc);
void Mstat(double *x, cov_model *cov, double *v);
void Mnonstat(double *x, double *y, cov_model *cov, double *v);
int checkM(cov_model *cov) ;
void rangeM(cov_model *cov, range_arraytype* ra);

void IdStat(double *x, cov_model *cov, double *v) ;
void IdNonStat(double *x, double *y, cov_model *cov, double *v);
int checkId(cov_model *cov) ;
void DId(double *x, cov_model *cov, double *v);
void DDId(double *x, cov_model *cov, double *v);
void TBM2Id(double *x, cov_model *cov, double *v);
int initspectralId(cov_model *cov);
void spectralId(cov_model *cov, spectral_storage *s, double *e);
void coinitId(cov_model *cov, localinfotype *li) ;
void ieinitId(cov_model *cov, localinfotype *li) ;
void rangeId(cov_model *cov, range_arraytype* ra); 


void kappaS(int i, cov_model *cov, int *nr, int *nc);
void Siso(double *x, cov_model *cov, double *v);
void Sstat(double *x, cov_model *cov, double *v);
void DS(double *x, cov_model *cov, double *v);
void DDS(double *x, cov_model *cov, double *v);
void tbm2S(double *x, cov_model *cov, double *v);
void Snonstat(double *x, double *y, cov_model *cov, double *v);  
int checkS(cov_model *cov);
void rangeS(cov_model *cov, range_arraytype* ra);
int initspectralS(cov_model *cov) ;
void spectralS(cov_model *cov, spectral_storage *s, double *e) ;
extern int SSTAT, SNONSTAT;
int initS(method_type *meth);
void doS(method_type *meth, res_type *res);
void coinitS(cov_model *cov, localinfotype *li);
void ieinitS(cov_model *cov, localinfotype *li);
void invS(double *x, cov_model *cov, double *v);
void nablaS(double *x, cov_model *cov, double *v);
void hessS(double *x, cov_model *cov, double *v);

void tbm2(double *x, cov_model *cov, double *v);
int checktbm2(cov_model *cov);
void tbm2num(double *x, cov_model *cov, double *v);
void rangetbm2(cov_model *cov, range_arraytype* ra);

void tbm3(double *x, cov_model *cov, double *v);
int checktbm3(cov_model *cov) ;
void rangetbm3(cov_model *cov, range_arraytype* ra);



//-----------------------------------------------------------------
void coin_ave(mpp_storage *s, cov_model *cov);

void mppinit_ave1(mpp_storage *s, cov_model *cov);
void kappa_ave1(int i, cov_model *cov, int *nr, int *nc);
void ave1(double *x, cov_model *cov, double *v) ;
int checkave1(cov_model *cov) ;
void rangeave1(cov_model *cov, range_arraytype* ra);
void ave1_logg(double *u, cov_model *cov, int dim, double *logg, double *sign);
double ave1_f(double *u, cov_model *cov, int dim) ;
res_type mppget_ave1(double *x, cov_model *cov, mpp_storage *s);
void sd_ave_stp(mpp_storage *s, cov_model *cov);
void sd_standard(mpp_storage *s, cov_model *cov);

void mppinit_ave2(mpp_storage *s, cov_model *cov);
void kappa_ave2(int i, cov_model *cov, int *nr, int *nc);
void ave2(double *x, cov_model *cov, double *v) ;
int checkave2(cov_model *cov) ;
void rangeave2(cov_model *cov, range_arraytype* ra);
res_type mppget_ave2(double *x, cov_model *cov, mpp_storage *s) ;


void co(double *x, cov_model *cov, double *v);
bool alternativeparam_co(cov_model *cov);
void range_co(cov_model *cov, range_arraytype* ra);
int check_co(cov_model *cov) ;


void kappa_cox1(int i, cov_model *cov, int *nr, int *nc);
void cox1(double *x, cov_model *cov, double *v) ;
void rangecox1(cov_model *cov, range_arraytype* ra);
int checkcox1(cov_model *cov);
int initspectralcox1(cov_model *cov) ;
void spectralcox1(cov_model *cov, spectral_storage *s, double *e) ; 
void coxnabla(double *x, cov_model *cov, double *v);
void coxhess(double *x, cov_model *cov, double *v);


void lp(double *x, cov_model *cov, double *v);
int checklp(cov_model *cov);
void rangelp(cov_model *cov, range_arraytype* ra);


void kappa_EAxxA(int i, cov_model *cov, int *nr, int *nc);
void EAxxA(double *x, cov_model *cov, double *v);
int checkEAxxA(cov_model *cov);
void rangeEAxxA(cov_model *cov, range_arraytype* ra);
double minEigenEAxxA(cov_model *cov);

void kappa_EtAxxA(int i, cov_model *cov, int *nr, int *nc);
void EtAxxA(double *x, cov_model *cov, double *v);
int checkEtAxxA(cov_model *cov);
void rangeEtAxxA(cov_model *cov, range_arraytype* ra);
double minEigenEtAxxA(cov_model *cov);

void kappa_rotat(int i, cov_model *cov, int *nr, int *nc);
void rotat(double *x, cov_model *cov, double *v);
int checkrotat(cov_model *cov);
void rangerotat(cov_model *cov, range_arraytype* ra);
double minEigenrotat(cov_model *cov);

void kappa_Rotat(int i, cov_model *cov, int *nr, int *nc);
void Rotat(double *x, cov_model *cov, double *v);
int checkRotat(cov_model *cov);
void rangeRotat(cov_model *cov, range_arraytype* ra);

void kappa_stp(int i, cov_model *cov, int *nr, int *nc);
void stp(double *x,  double *y, cov_model *cov, double *v) ;
int checkstp(cov_model *cov);
void rangestp(cov_model *cov, range_arraytype* ra);
void mppinit_stp(mpp_storage *s, cov_model *cov) ;
void coin_stp(mpp_storage *s, cov_model *cov);
res_type mppget_stp(double *x, double *u, cov_model *cov, mpp_storage *s);

void kappa_rational(int i, cov_model *cov, int *nr, int *nc);
double minEigenrational(cov_model *cov);
void rational(double *x, cov_model *cov, double *v);
int checkrational(cov_model *cov);
void rangerational(cov_model *cov, range_arraytype* ra);

void Stein(double *x, cov_model *cov, double *v);
bool alternativeparam_Stein(cov_model *cov);
int check_Stein(cov_model *cov);
void range_Stein(cov_model *cov, range_arraytype* ra);


void MaStein(double *x, cov_model *cov, double *v);
int check_MaStein(cov_model *cov);
void range_MaStein(cov_model *cov, range_arraytype* ra);

void kappaNonStWM(int i, cov_model *cov, int *nr, int *nc);
void NonStWMQ(double *x, double *y,  double srqtQ, cov_model *cov, double *v);
void NonStWM(double *x, double *y, cov_model *cov, double *v);
int checkNonStWM(cov_model *cov) ; 
void rangeNonStWM(cov_model *cov, range_arraytype* ra); 
double DrawLogMixNonStWM(cov_model *cov, mpp_storage *s);
double LogMixWeightNonStWM(double *x, cov_model *cov, mpp_storage *s);


/* nsst */
/* Tilmann Gneiting's space time models, part I */
void nsst(double *x, cov_model *cov, double *v);
void Dnsst(double *x, cov_model *cov, double *v);
void TBM2nsst(double *x, cov_model *cov, double *v);
int checknsst(cov_model *cov);
void rangensst(cov_model *cov, range_arraytype* ra);


/* undefined model -- for technical reasons only */
void rangeundefined(cov_model *cov, range_arraytype* ra);
int checkundefined(cov_model *cov);
int checkOK(cov_model *cov);

int initstandard(method_type *meth);
void dostandard(method_type *meth, res_type *res);

void kappashift(int i, cov_model *cov, int *nr, int *nc);
void shift(double *x, cov_model *cov, double *v);
int checkshift(cov_model *cov);
void rangeshift(cov_model *cov, range_arraytype* ra);

void vector(double *x, cov_model *cov, double *v);
void vectorAniso(double *x, cov_model *cov, double *v);
int checkvector(cov_model *cov);
void rangevector(cov_model *cov, range_arraytype* ra);

void div(double *x, cov_model *cov, double *v);
void curl(double *x, cov_model *cov, double *v);
int checkdivcurl(cov_model *cov);
void rangedivcurl(cov_model *cov, range_arraytype* ra);


void kappaqam(int i, cov_model *cov, int *nr, int *nc);
void qam(double *x, cov_model *cov, double *v);
int checkqam(cov_model *cov) ;
void rangeqam(cov_model *cov, range_arraytype* ra);

void kappamqam(int i, cov_model *cov, int *nr, int *nc);
void mqam(double *x, cov_model *cov, double *v);
int checkmqam(cov_model *cov) ;
void rangemqam(cov_model *cov, range_arraytype* ra);


void Exp(double *x, cov_model *cov, double *v);
void DExp(double *x, cov_model *cov, double *v);
void DDExp(double *x, cov_model *cov, double *v);
int checkExp(cov_model *cov) ;
void rangeExp(cov_model *cov, range_arraytype* ra); 


void Pow(double *x, cov_model *cov, double *v);
void DPow(double *x, cov_model *cov, double *v);
void DDPow(double *x, cov_model *cov, double *v);
int checkPow(cov_model *cov) ;
void rangePow(cov_model *cov, range_arraytype* ra); 
void invPow(double *x, cov_model *cov, double *v);


void ma1(double *x, cov_model *cov, double *v);
int checkma1(cov_model *cov);
void rangema1(cov_model *cov, range_arraytype* ra);
void ma2(double *x, cov_model *cov, double *v);
void rangema2(cov_model *cov, range_arraytype* ra);
int checkma2(cov_model *cov);


////////////////////////////////////////////////////////////////////
//   covariance models only used in MLE
////////////////////////////////////////////////////////////////////

//void Hatstat(double *x, cov_model *cov, double *v);
//void Hatnonstat(double *x, double *y, cov_model *cov, double *v);
//int checkHat(cov_model *cov);
//void rangeHat(cov_model *cov, range_arraytype* ra);

//void Outstat(double *x, cov_model *cov, double *v);
//void Outnonstat(double *x, double *y, cov_model *cov, double *v);
//int checkOut(cov_model *cov);
//void rangeOut(cov_model *cov, range_arraytype* ra);


//void splitstat(double *x, cov_model *cov, double *v);
//void splitnonstat(double *x, double *y, cov_model *cov, double *v);
//int checksplit(cov_model *cov);
///void rangesplit(cov_model *cov, range_arraytype* ra);


void iso2iso_MLE(double *x, cov_model *cov, double *v) ;
void spiso2spiso_MLE(double *x, cov_model *cov, double *v) ;
void spacetime2iso_MLE(double *x, cov_model *cov, double *v) ;
void Stat2iso_MLE(double *x, cov_model *cov, double *v) ;
void Nonstat2iso_MLE(double *x, double *y, cov_model *cov, double *v) ;
void Stat2spacetime_MLE(double *x, cov_model *cov, double *v) ;
void Nonstat2spacetime_MLE(double *x, double *y, cov_model *cov, double *v) ;
void Stat2Stat_MLE(double *x, cov_model *cov, double *v) ;
void Nonstat2Stat_MLE(double *x, double *y, cov_model *cov, double *v) ;
void Siso_MLE(double *x, cov_model *cov, double *v);
void Sstat_MLE(double *x, cov_model *cov, double *v);
void Snonstat_MLE(double *x, double *y, cov_model *cov, double *v);
void DS_MLE(double *x, cov_model *cov, double *v);
void DDS_MLE(double *x, cov_model *cov, double *v);


void natsc(double *x, cov_model *cov, double *v);
void Dnatsc(double *x, cov_model *cov, double *v);
void DDnatsc(double *x, cov_model *cov, double *v);
void invnatsc(double *x, cov_model *cov, double *v) ;
void coinitnatsc(cov_model *cov, localinfotype *li) ;
void ieinitnatsc(cov_model *cov, localinfotype *li) ;
void tbm2natsc(double *x, cov_model *cov, double *v);
int checknatsc(cov_model *cov) ;
int initspectralnatsc(cov_model *cov) ;
void spectralnatsc(cov_model *cov, spectral_storage *s, double *e) ;
void rangenatsc(cov_model *cov, range_arraytype* ra);
int initnatsc(method_type *meth);
void donatsc(method_type *meth, res_type *res);




void mixed(double *x,  cov_model *cov, double *v);
void MLEmixed(double *x,  cov_model *cov, double *v);
void mixed_nonstat(double *x,  double *y, cov_model *cov, double *v);
void MLEmixed_nonstat(double *x,  double *y, cov_model *cov, double *v);
int  checkmixed(cov_model *cov) ;
void rangemixed(cov_model *cov, range_arraytype* ra);
int  initmixed(method_type *meth);
void domixed(method_type *meth, res_type *res);


void X(double *x, cov_model *cov, double *v);
void Xnonstat(double *x, double *y, cov_model *cov, double *v);
int checkX(cov_model *cov);
void rangeX(cov_model *cov, range_arraytype* ra);


#endif /* Operators_H */
 
