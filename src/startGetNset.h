
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

#ifndef RFstartGetNset_H
#define RFstartGetNset_H 1


// extern int etc, if other are enabled to use these functions?
int IncludePrim(const char *name, Types type, int kappas,  
		domain_type domain, isotropy_type isotropy,
		checkfct check, rangefct range,
		int maxdim, 
		ext_bool finiterange, monotone_type monotonicity);
int IncludePrim(const char *name, Types type, int kappas,  size_fct kappasize,
		domain_type domain, isotropy_type isotropy,
		checkfct check, rangefct range,
		int maxdim, 
		ext_bool finiterange, monotone_type monotonicity);
int IncludePrim(const char *name, Types type, int kappas,  
		domain_type domain, isotropy_type isotropy,
		checkfct check, rangefct range,
		int vdim, int maxdim,
		ext_bool finiterange, monotone_type monotonicity);
int IncludePrim(const char *name,  Types type, int kappas,  size_fct kappasize,
		domain_type domain, isotropy_type isotropy,
		checkfct check, rangefct range,
		int vdim, int maxdim,
		ext_bool finiterange, monotone_type monotonicity);

int IncludePrim(const char *name,Types type,  int kappas, 
		 domain_type domain, isotropy_type isotropy,	
		 checkfct check, rangefct range, pref_type pref,
		 int vdim, int maxdim,
		ext_bool finiterange, monotone_type monotonicity);
int IncludePrim(const char *name, Types type, int kappas, size_fct kappasize,
		 domain_type domain, isotropy_type isotropy,	
		 checkfct check, rangefct range, pref_type pref,
		int vdim, int maxdim,
		ext_bool finiterange, monotone_type monotonicity);

void make_internal();


int IncludeModel(const char *name, Types type, int minsub, int maxsub,
		 int kappas, size_fct kappasize,
		 domain_type domain, isotropy_type isotropy,
		 checkfct check, rangefct range, pref_type pref, 
		 int internal, int vdim, int maxdim, ext_bool finiterange,
		 monotone_type monotonicity);
int IncludeScalar(const char *name, Types type, int minsub, int maxsub, 
		 int kappas, domain_type domain, isotropy_type isotropy,
		 checkfct check, rangefct range, pref_type pref, 
		 int maxdim, ext_bool finiterange, monotone_type monotonicity);
int CopyModel(const char *name, int which);
int CopyModel(const char *name, int which, Types type);
int CopyModel(const char *name, int which, checkfct check);
void nickname(const char *nick);
void AddVariant(Types type, isotropy_type iso);
void addCov(covfct cf, covfct D, covfct inverse);
void addCov(covfct cf, covfct D, covfct D2, covfct inverse);
void addCov(covfct cf, covfct D, covfct D2, nonstat_inv inverse);
void addCov(covfct cf, covfct D, covfct D2, 
	    covfct inverse, nonstat_inv nonstat_inverse);
void addCov(covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse);
void addCov(covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse,  nonstat_inv nonstat_inverse);
void addCov( int F_derivs, covfct cf, covfct D, covfct inverse);
void addCov(int F_derivs, covfct cf, covfct D, covfct D2, covfct inverse,
	    nonstat_inv nonstatinverse);
void addCov( int F_derivs, covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse);
void addCov(nonstat_covfct cf);
void addCov( int F_derivs , nonstat_covfct cf);
void addCov(aux_covfct auxcf);
void addCov(covfct distrD, covfct logdistrD, nonstat_inv Dinverse,
	    covfct distrP, nonstat_covfct distrP2sided,
	    covfct distrQ, covfct distrR, nonstat_covfct distrR2sided);

void addlogCov(logfct log, nonstat_logfct nonstatlog, 
	       nonstat_inv lognonstat_inverse);
void addlogCov(logfct log);
void addlogCov(nonstat_logfct nonstatlog);


void addlogD(covfct logD);

void nablahess(covfct nabla, covfct hess);
void addSpaceD(covfct spaceD);
void addLocal(getlocalparam coinit, getlocalparam ieinit);
void addCallLocal(altlocalparam alt);
//void addTBM(covfct tbm2, covfct spaceD);
int addTBM(covfct tbm2);
void addTBM(covfct tbm2, initfct Init, spectral_do spectral);
void addTBM(initfct Init, spectral_do spectral);
void addHyper(hyper_pp_fct hyper_pp);
void addSpecific(int cov);
void addSpecific(int cov, bool copy);
void RandomShape(int, structfct Struct, initfct init, dofct Do,
		 do_random_fct DoRandom, 
		 bool average, bool randomcoin, bool specific);
void RandomShape(int, structfct Struct, initfct Init, dofct Do, 
		 bool average, bool randomcoin, bool specific);
void RandomShape(int, structfct Struct, initfct Init, do_random_fct DoRandom);
void RandomShape(int, structfct Struct, initfct init, dofct Do);
void RandomShape(int maxmoments, structfct Struct, initfct Init,
		 dofct Do, finaldofct FinalDo);
void RandomShape(structfct Struct, bool average);
void RandomShape(structfct Struct);
void RandomShape(int, initfct init, dofct Do, bool average);
void RandomShape(int, initfct init, dofct Do);
void addSpecial(minmaxfct minmaxeigen);
void addGaussMixture(draw_random drawmix, log_mixdens logmixdens);

int addFurtherCov(covfct cf, covfct D);
int addFurtherCov(covfct cf, covfct D, covfct D2);
int addFurtherCov(nonstat_covfct cf);
void add_sortof(sortof_fct sortof);
void change_sortof(int i, sortsofparam sort);
void change_typeof(int i, Types type);
void change_typeof(int i, Types type, const char *names[]);
bool allowsRandomParam(model *cov, int i);
bool allowsShapeParam(model *cov, int i);
bool IsParamTrend(model *cov, int i);


void addReturns(//return_fct Covariance, ext_bool_ret_fct isCovariance, 
		return_covmat CovMatrix, ext_bool_ret_fct isCovMatrix
		//tworeturns_fct InverseCovMatrix,
		//ext_bool_ret_fct isInverseCovMatrix,
		//return_fct Variogram, ext_bool_ret_fct isVariogram,
		//return_fct PseudoVariogram, ext_bool_ret_fct isPseudoVariogram
		);


void kappanames(const char* n1, SEXPTYPE t1);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15);
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15, const char* n16, SEXPTYPE t16);

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15, const char* n16, SEXPTYPE t16,
		const char* n17, SEXPTYPE t17, const char* n18, SEXPTYPE t18
		);

void subnames(const char* n1);
void subnames(const char* n1, const char* n2);
void subnames(const char* n1, const char* n2, const char* n3);
void subnames(const char* n1, const char* n2, const char* n3, const char* n4);
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5);
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5, const char* n6);

int checkMissing(model *cov);

void ScaleOne(double *x, model VARIABLE_IS_NOT_USED *cov, double *v);
int struct_failed(model *cov, model **atom);
int init_failed(model *cov, gen_storage *s);
int init_statiso(model *cov, gen_storage *s);
void do_failed(model *cov, gen_storage *s);
void do_statiso(model *cov, gen_storage VARIABLE_IS_NOT_USED *s);
int checkOK(model *cov);
void ErrInverse(double *v, model *cov, double *x);
void InverseIsotropic(double *U, model *cov, double *inverse);

int structOK(model *cov, model **newmodel);
int initOK(model *cov, gen_storage *s);
void doOK(model *cov, gen_storage *s);
void do_random_ok(model *cov, double *v);
void do_random_failed(model *cov, double *v);

					
void setptwise(ptwise_type pt);

void includeStandardMath();


void Taylor(double c, double pow);
void TailTaylor( double t, double tpow, double texp, double texppow);
void Taylor(double c, double pow, double c1, double pow1);
void Taylor(double c, double pow, double c1, double pow1, 
	    double c2, double pow2);

int checkNotOK(model VARIABLE_IS_NOT_USED *cov);

void setDI(allowedD_fct D, allowedI_fct I, setDI_fct f);
void addTypeFct(type_fct TypeFct);

bool allowedDstandard(model *cov);
bool allowedIstandard(model *cov);
bool allowedPrevModelI(model *cov);

#endif
