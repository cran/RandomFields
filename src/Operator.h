#ifndef Operators_H
#define Operators_H 1

#include "primitive.h"
#include "shape_processes.h"
#include "Families.h"
#include "SpherModels.h"
#include "Coordinate_systems.h"


void D_2(double *x, cov_model *cov, double *v);
void DD_2(double *x, cov_model *cov, double *v);
//int check2(cov_model *cov);
void inverse2(double *x, cov_model *cov, double *v);
void nonstatinverse2(double *v, cov_model *cov, double *left, double *right);
void nonstat_loginverse2(double *v, cov_model *cov, double *x, double *y);
//int struct2(cov_model *cov, cov_model **newmodel);
int struct2(cov_model *cov, cov_model **newmodel);
int init2(cov_model *cov, gen_storage *s);
void do2(cov_model *cov, gen_storage *s);
void dorandom2(cov_model *cov, double *v);

void plusStat(double *x, cov_model *cov, double *v);
void plusNonStat(double *x,  double *y, cov_model *cov, double *v);
void covmatrix_plus(cov_model *cov, double *v, int*);
char iscovmatrix_plus(cov_model *cov);
int checkplus(cov_model *cov);
bool Typeplus(Types required, cov_model *cov, int depth);
void Dplus(double *x, cov_model *cov, double *v);
void DDplus(double *x, cov_model *cov, double *v);
void spectralplus(cov_model *cov, gen_storage *s, double *e);
int structplus(cov_model *cov, cov_model **newmodel);
int initplus(cov_model *cov, gen_storage *s);
void doplus(cov_model *cov, gen_storage *s);


void malStat(double *x, cov_model *cov, double *v);
void malNonStat(double *x,  double *y, cov_model *cov, double *v);
void logmalStat(double *x, cov_model *cov, double *v, double *Sign);
void logmalNonStat(double *x, double *y, cov_model *cov, double *v, 
		   double *Sign);
int checkmal(cov_model *cov);
bool Typemal(Types required, cov_model *cov, int depth);
void Dmal(double *x, cov_model *cov, double *v);
int initmal(cov_model *cov, gen_storage *s);
void domal(cov_model *cov, gen_storage *s);

void kappaM(int i, cov_model *cov, int *nr, int *nc);
void Mstat(double *x, cov_model *cov, double *v);
void Mnonstat(double *x, double *y, cov_model *cov, double *v);
int checkM(cov_model *cov);
sortsofparam paramtype_M(int k, int row, int col);
void rangeM(cov_model *cov, range_type* ra);

void kappaSchur(int i, cov_model *cov, int *nr, int *nc);
void Schurstat(double *x, cov_model *cov, double *v);
void Schurnonstat(double *x, double *y, cov_model *cov, double *v);
int checkSchur(cov_model *cov);
void rangeSchur(cov_model *cov, range_type* ra);

void IdStat(double *x, cov_model *cov, double *v);
void IdNonStat(double *x, double *y, cov_model *cov, double *v);
void IdInverse(double *x, cov_model *cov, double *v);
int checkId(cov_model *cov);
void DId(double *x, cov_model *cov, double *v);
void DDId(double *x, cov_model *cov, double *v);
void TBM2Id(double *x, cov_model *cov, double *v);
int initId(cov_model *cov, gen_storage *s);
void spectralId(cov_model *cov, gen_storage *s, double *e);
void coinitId(cov_model *cov, localinfotype *li);
void ieinitId(cov_model *cov, localinfotype *li);
void rangeId(cov_model *cov, range_type* ra); 


void kappaS(int i, cov_model *cov, int *nr, int *nc);
void Siso(double *x, cov_model *cov, double *v);
void logSiso(double *x, cov_model *cov, double *v, double * Sign);
void Sstat(double *x, cov_model *cov, double *v);
void covmatrixS(cov_model *cov, double *v, int*);
char iscovmatrixS(cov_model *cov);
void logSstat(double *x, cov_model *cov, double *v, double * Sign);
void DS(double *x, cov_model *cov, double *v);
void DDS(double *x, cov_model *cov, double *v);
void D3S(double *x, cov_model *cov, double *v);
void D4S(double *x, cov_model *cov, double *v);
void tbm2S(double *x, cov_model *cov, double *v);
void Snonstat(double *x, double *y, cov_model *cov, double *v);  
void logSnonstat(double *x, double *y, cov_model *cov, double *v, double *);  
int checkS(cov_model *cov);
bool TypeS(Types required, cov_model *cov, int depth);
void rangeS(cov_model *cov, range_type* ra);
void spectralS(cov_model *cov, gen_storage *s, double *e);
extern int SSTAT, SNONSTAT;
void coinitS(cov_model *cov, localinfotype *li);
void ieinitS(cov_model *cov, localinfotype *li);
void inverseS(double *x, cov_model *cov, double *v);
void nonstatinverseS(double *x, cov_model *cov, double *left, double*right);
void nonstat_loginverseS(double *v, cov_model *cov, double *x, double *y);
void nablaS(double *x, cov_model *cov, double *v);
void hessS(double *x, cov_model *cov, double *v);
int structS(cov_model *cov, cov_model **newmodel);
int initS(cov_model *cov, gen_storage *s);
void doS(cov_model *cov, gen_storage *s);





void binary(double *x, cov_model *cov, double *v);
int checkbinary(cov_model *cov);
void rangebinary(cov_model *cov, range_type *range);

void brownresnick(double *x, cov_model *cov, double *v);
int checkbrownresnick(cov_model *cov);
void Dbrownresnick(double *x, cov_model *cov, double *v);
void DDbrownresnick(double *x, cov_model *cov, double *v);
void D3brownresnick(double *x, cov_model *cov, double *v);
int struct_brownresnick(cov_model *cov, cov_model **newmodel);
int init_brownresnick(cov_model *cov, gen_storage *s);
void do_brownresnick(cov_model *cov, gen_storage *s);

void BR2BG(double *x, cov_model *cov, double *v);
int check_BR2BG(cov_model *cov);

void BR2EG(double *x, cov_model *cov, double *v);
int check_BR2EG(cov_model *cov);


void kappa_cox(int i, cov_model *cov, int *nr, int *nc);
void cox(double *x, cov_model *cov, double *v);
void rangecox(cov_model *cov, range_type* ra);
int checkcox(cov_model *cov);
int initcox(cov_model *cov, gen_storage *s);
void spectralcox(cov_model *cov, gen_storage *s, double *e); 
void coxnabla(double *x, cov_model *cov, double *v);
void coxhess(double *x, cov_model *cov, double *v);

void extremalgaussian(double *x, cov_model *cov, double *v);
int check_extremalgaussian(cov_model *cov);

void binaryGauss(double *x, cov_model *cov, double *v);
int check_binaryGauss(cov_model *cov);

void extrgauss(double *x, cov_model *cov, double *v);
  int check_extrgauss(cov_model *cov);

void kappaNonStWM(int i, cov_model *cov, int *nr, int *nc);
void NonStWMQ(double *x, double *y,  double srqtQ, cov_model *cov, double *v);
void NonStWM(double *x, double *y, cov_model *cov, double *v);
int checkNonStWM(cov_model *cov); 
sortsofparam paramtype_nonstWM(int k, int row, int col);
void rangeNonStWM(cov_model *cov, range_type* ra); 
void DrawMixNonStWM(cov_model *cov, double *random);
double LogMixWeightNonStWM(double *x, double logV, cov_model *cov);


void lp(double *x, cov_model *cov, double *v);
int checklp(cov_model *cov);
void rangelp(cov_model *cov, range_type* ra);

void MaStein(double *x, cov_model *cov, double *v);
int check_MaStein(cov_model *cov);
void range_MaStein(cov_model *cov, range_type* ra);

/* nsst */
/* Tilmann Gneiting's space time models, part I */
void nsst(double *x, cov_model *cov, double *v);
void Dnsst(double *x, cov_model *cov, double *v);
void TBM2nsst(double *x, cov_model *cov, double *v);
int checknsst(cov_model *cov);
sortsofparam paramtype_nsst(int k, int row, int col);
void rangensst(cov_model *cov, range_type* ra);

void nonstatgennsst(double *x, double *y, cov_model *cov, double *v);
int checkgennsst(cov_model *cov);

void kappa_gennsst_intern(int i, cov_model *cov, int *nr, int *nc);
void nonstatgennsst_intern(double *x, double *y, cov_model *cov, double *v);
int checkgennsst_intern(cov_model *cov);
void range_gennsst_intern(cov_model *cov, range_type* range);

/* undefined model -- for technical reasons only */
//void rangeundefined(cov_model *cov, range_type* ra);
//int checkundefined(cov_model *cov);
int checkOK(cov_model *cov);

void kappashift(int i, cov_model *cov, int *nr, int *nc);
void shift(double *x, cov_model *cov, double *v);
int checkshift(cov_model *cov);
void rangeshift(cov_model *cov, range_type* ra);

void vector(double *x, cov_model *cov, double *v);
void vectorAniso(double *x, cov_model *cov, double *v);
int checkvector(cov_model *cov);
void rangevector(cov_model *cov, range_type* ra);

void div(double *x, cov_model *cov, double *v);
void curl(double *x, cov_model *cov, double *v);
int checkdivcurl(cov_model *cov);
//void rangedivcurl(cov_model *cov, range_type* ra);


void Exp(double *x, cov_model *cov, double *v);
void DExp(double *x, cov_model *cov, double *v);
void DDExp(double *x, cov_model *cov, double *v);
int checkExp(cov_model *cov);
void rangeExp(cov_model *cov, range_type *range);

void ma1(double *x, cov_model *cov, double *v);
int checkma1(cov_model *cov);
void rangema1(cov_model *cov, range_type* ra);
void ma2(double *x, cov_model *cov, double *v);
int checkma2(cov_model *cov);


////////////////////////////////////////////////////////////////////
//   covariance models only used in MLE
////////////////////////////////////////////////////////////////////

//void Hatstat(double *x, cov_model *cov, double *v);
//void Hatnonstat(double *x, double *y, cov_model *cov, double *v);
//int checkHat(cov_model *cov);
//void rangeHat(cov_model *cov, range_type* ra);

//void Outstat(double *x, cov_model *cov, double *v);
//void Outnonstat(double *x, double *y, cov_model *cov, double *v);
//int checkOut(cov_model *cov);
//void rangeOut(cov_model *cov, range_type* ra);


//void splitstat(double *x, cov_model *cov, double *v);
//void splitnonstat(double *x, double *y, cov_model *cov, double *v);
//int checksplit(cov_model *cov);
///void rangesplit(cov_model *cov, range_type* ra);

/*
void iso2iso_MLE(double *x, cov_model *cov, double *v);
void spiso2spiso_MLE(double *x, cov_model *cov, double *v);
void spacetime2iso_MLE(double *x, cov_model *cov, double *v);
void Stat2iso_MLE(double *x, cov_model *cov, double *v);
void Nonstat2iso_MLE(double *x, double *y, cov_model *cov, double *v);
void Stat2spacetime_MLE(double *x, cov_model *cov, double *v);
void Nonstat2spacetime_MLE(double *x, double *y, cov_model *cov, double *v);
void Stat2Stat_MLE(double *x, cov_model *cov, double *v);
void Nonstat2Stat_MLE(double *x, double *y, cov_model *cov, double *v);
void Siso_MLE(double *x, cov_model *cov, double *v);
void Sstat_MLE(double *x, cov_model *cov, double *v);
void DS_MLE(double *x, cov_model *cov, double *v);
void DDS_MLE(double *x, cov_model *cov, double *v);
*/

void natsc(double *x, cov_model *cov, double *v);
void Dnatsc(double *x, cov_model *cov, double *v);
void DDnatsc(double *x, cov_model *cov, double *v);
void Inversenatsc(double *x, cov_model *cov, double *v);
void coinitnatsc(cov_model *cov, localinfotype *li);
void ieinitnatsc(cov_model *cov, localinfotype *li);
void tbm2natsc(double *x, cov_model *cov, double *v);
int checknatsc(cov_model *cov);
int initnatsc(cov_model *cov);
void spectralnatsc(cov_model *cov, gen_storage *s, double *e);
int initnatsc(cov_model *cov, gen_storage *s);
void donatsc(cov_model *cov, gen_storage *s);

void NullModel(double *x, cov_model *cov, double *v);
bool TypeNullModel(Types required, cov_model *cov, int depth);
void rangeNullModel(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range);


void Pow(double *x, cov_model *cov, double *v);
void DPow(double *x, cov_model *cov, double *v);
void DDPow(double *x, cov_model *cov, double *v);
int checkPow(cov_model *cov);
void rangePow(cov_model *cov, range_type* ra); 
void InversePow(double *x, cov_model *cov, double *v);

void shapePow(double *x, cov_model *cov, double *v);
void DshapePow(double *x, cov_model *cov, double *v);
void DDshapePow(double *x, cov_model *cov, double *v);
void InverseShapePow(double *x, cov_model *cov, double *v);



void PowSstat(double *x, cov_model *cov, double *v);
void logPowSstat(double *x, cov_model *cov, double *v, double *Sign);
void PowSnonstat(double *x, double *y, cov_model *cov, double *v);
void logPowSnonstat(double *x, double *y, cov_model *cov, double *v, 
		 double *Sign);
void inversePowS(double *x, cov_model *cov, double *v) ;
int checkPowS(cov_model *cov) ;
bool TypePowS(Types required, cov_model *cov, int depth) ;
void rangePowS(cov_model *cov, range_type* range);
void PowScaleToLoc(cov_model *to, cov_model *from, int VARIABLE_IS_NOT_USED depth) ;
int structPowS(cov_model *cov, cov_model **newmodel) ;
int initPowS(cov_model *cov, gen_storage *s);
void doPowS(cov_model *cov, gen_storage *s);



void kappaqam(int i, cov_model *cov, int *nr, int *nc);
void qam(double *x, cov_model *cov, double *v);
int checkqam(cov_model *cov);
sortsofparam paramtype_qam(int k, int row, int col);
void rangeqam(cov_model *cov, range_type* ra);

void kappamqam(int i, cov_model *cov, int *nr, int *nc);
void mqam(double *x, cov_model *cov, double *v);
int checkmqam(cov_model *cov);
void rangemqam(cov_model *cov, range_type* ra);


void select(double *x, cov_model *cov, double *v);
void covmatrix_select(cov_model *cov, double *v, int*);
char iscovmatrix_select(cov_model *cov);
int checkselect(cov_model *cov);
void rangeselect(cov_model *cov, range_type *range);

void tbm(double *x, cov_model *cov, double *v);
void Dtbm(double *x, cov_model *cov, double *v);
int checktbmop(cov_model *cov);
void rangetbm_common(cov_model *cov, range_type *range, bool tbmop);
void rangetbmop(cov_model *cov, range_type* ra);


void co(double *x, cov_model *cov, double *v);
int check_co(cov_model *cov);
bool alternativeparam_co(cov_model *cov);
void range_co(cov_model *cov, range_type* ra);


void Stein(double *x, cov_model *cov, double *v);
int check_Stein(cov_model *cov);
bool alternativeparam_Stein(cov_model *cov);
void range_Stein(cov_model *cov, range_type* ra);

void strokorb(double *x, cov_model *cov, double *v);
int checkstrokorb(cov_model *cov);
int init_strokorb(cov_model *cov, gen_storage *s);
void do_strokorb(cov_model *cov, gen_storage *s);


//void strokorbBall(double *x, cov_model *cov, double *v);
int checkstrokorbBall(cov_model *cov);
int struct_strokorbBall(cov_model *cov, cov_model **newmodel);  
//int init_strokorbBall(cov_model *cov, gen_storage *s);
//void do_strokorbBall(cov_model *cov, gen_storage *s);
void rangestrokorbball(cov_model  VARIABLE_IS_NOT_USED *cov, range_type *range);


void strokorbBallInner(double *x, cov_model *cov, double *v);
int check_strokorbBallInner(cov_model *cov);
void range_strokorbBallInner(cov_model *cov, range_type *range);
int init_strokorbBallInner(cov_model *cov, gen_storage *s);
void do_strokorbBallInner(cov_model *cov, gen_storage *s);


void strokorbPoly(double *x, cov_model *cov, double *v);
int checkstrokorbPoly(cov_model *cov);
int struct_strokorbPoly(cov_model *cov, cov_model **newmodel); 


void mult_inverse(double *x, cov_model *cov, double *v);
void mult_inverseNonstat(double *x, double *y, cov_model *cov, double *v);
int checkmult_inverse(cov_model *cov);


void addSetParam(cov_model **newmodel, cov_model * remote, 
		 param_set_fct set,  bool performdo, int variant);
void addSetDistr(cov_model **newmodel, cov_model * remote, 
		 param_set_fct set,  bool performdo, int variant);
void setparamStat(double *x, cov_model *cov, double *v);
void setparamNonStat(double *x,  double *y, cov_model *cov, double *v);
void covmatrix_setparam(cov_model *cov, double *v, int *nonzeros);
char iscovmatrix_setparam(cov_model *cov);
int checksetparam(cov_model *cov);
void range_setparam(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range);
void Inverse_setparam(double *v, cov_model *cov, double *x);
void NonstatInverse_setparam(double *v, cov_model *cov, double *x, double *y);
void LogNonstatInverse_setparam(double *v, cov_model *cov, double *x, double *y);
bool Typesetparam(Types required, cov_model *cov, int depth);
void Dsetparam(double *x, cov_model *cov, double *v);
void DDsetparam(double *x, cov_model *cov, double *v);
void D3setparam(double *x, cov_model *cov, double *v);
void D4setparam(double *x, cov_model *cov, double *v);
void spectralsetparam(cov_model *cov, gen_storage *s, double *e);
int initsetparam(cov_model *cov, gen_storage *s);
void dosetparam(cov_model *cov, gen_storage *s);

void kappaScatter(int i, cov_model *cov, int *nr, int *nc);
void Scatter(double *xx, cov_model *cov, double *v);
int checkScatter(cov_model *cov);
void rangeScatter(cov_model *cov, range_type *range); 
int struct_scatter(cov_model *cov, cov_model **newmodel);
int init_scatter(cov_model *cov, gen_storage *s);


void oesting(double *x, cov_model *cov, double *v);
void Doesting(double *x, cov_model *cov, double *v);
void DDoesting(double *x, cov_model *cov, double *v); 
int checkoesting(cov_model *cov);
int initoesting(cov_model *cov, gen_storage VARIABLE_IS_NOT_USED *s);
void rangeoesting(cov_model *cov, range_type *range);

void idcoord(double *x, cov_model *cov, double *v);
int checkidcoord(cov_model *cov);

void trafo(double *x, cov_model *cov, double *v);  
void logtrafo(double *x, cov_model *cov, double *v, 
		     double *Sign);  
void nonstattrafo(double *x, double *y, cov_model *cov, double *v);  
void lognonstattrafo(double *x, double *y, cov_model *cov, double *v, 
		     double *Sign);  
int checktrafo(cov_model *cov);
void rangetrafo(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range);
bool Typetrafo(Types required, cov_model *cov, int depth);


void covariate(double *x, cov_model *cov, double *v);
int checkcovariate(cov_model *cov);
void rangecovariate(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range);
void kappa_covariate(int i, cov_model *cov, int *nr, int *nc);


void nonstatprod(double *x, double *y, cov_model *cov, double *v);
int checkprod(cov_model *cov);
int checkprodproc(cov_model *cov);
int structprodproc(cov_model  *cov, cov_model **newmodel);
int initprodproc(cov_model *cov, gen_storage *s);
void doprodproc(cov_model *cov, gen_storage *s);

void nonstatsum(double *x, double *y, cov_model *cov, double *v);
int checksum(cov_model *cov);

#endif /* Operators_H */
 
