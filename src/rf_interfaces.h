
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

#ifndef RF_INTERFACE_H
#define RF_INTERFACE_H 1




#define XLIST_X 0
#define XLIST_Y 1
#define XLIST_T 2
#define XLIST_GRID 3
#define XLIST_SPATIALDIM 4
#define XLIST_TIME 5
#define XLIST_DIST 6
#define XLIST_RESTOT 7
#define XLIST_L 8
#define XLIST_UNITS 9
#define XLIST_NEWUNITS 10

model *InitIntern(int cR, SEXP Model, SEXP x, bool NA_OK);

void simulate(double *x, model *cov, double *v);
int check_simulate(model *cov); 
void range_simulate(model *cov, range_type *range);  
int struct_simulate(model *cov, model **newmodel);
int init_simulate(model *cov, gen_storage *S);
// void do_simulate(model *cov, gen_storage *S);
void range_simulate(model VARIABLE_IS_NOT_USED *cov, range_type* range);

void density(double *x, model *cov, double *v);
int check_density(model *cov); 
void range_density(model *cov, range_type *range);  
int struct_density(model *cov, model **newmodel);
int init_density(model *cov, gen_storage *S);
// void do_density(model *cov, gen_storage *S);
void range_density(model VARIABLE_IS_NOT_USED *cov, range_type* range);



#define LIKELIHOOD_DATA 0
#define LIKELIHOOD_NA_VAR 1
#define LIKELIHOOD_BETASSEPARATE 2
#define LIKELIHOOD_IGNORETREND 3
#define LIKELIHOOD_LAST LIKELIHOOD_IGNORETREND
#define LIKELI_NA_INTEGER false
#define LIKELI_EXCLUDE_TREND true
void kappalikelihood(int i, model VARIABLE_IS_NOT_USED *cov, 
		     int *nr, int *nc);
void likelihood(double *data, model *cov, double *v);
int check_likelihood(model *cov);
int struct_likelihood(model *cov, model **newmodel);
void range_likelihood(model *cov, range_type* range);

void linearpart(double *data, model *cov, double *v);
int check_linearpart(model *cov);
int struct_linearpart(model *cov, model **newmodel);
//void range_linearpart(model *cov, range_type* range);


void predict(double VARIABLE_IS_NOT_USED *x, model *cov, double *v);
int check_predict(model *predict);
int struct_predict(model *cov, model VARIABLE_IS_NOT_USED  **newmodel);
void range_predict(model VARIABLE_IS_NOT_USED *cov, range_type* range);


void Cov(double *x, model *cov, double *value);
int check_cov(model *cov);
int struct_cov(model *cov, model **newmodel);
int init_cov(model *cov, gen_storage *s);


void FctnIntern(model *cov, model *covVdim, model *sub,
		double *value, bool ignore_y);
void FctnExtern(model *cov, model *covVdim, model *sub,
		double *value, bool ignore_y);
void Fctn(double *x, model *cov, double *value);
int check_fctn(model *cov);
//int check_fct_intern(model *cov, Types type, bool close, bool kernel,
//		     int rows, int cols, Types frame);

void CovMatrix(double *x, model *cov, double *value);
int check_covmatrix(model *cov) ;

void EvalDistr(double *x, model *cov, double *v);
void kappa_EvalDistr(int i, model *cov, int *nr, int *nc);
int check_EvalDistr(model *cov); 
void range_EvalDistr(model *cov, range_type *range);  
int struct_EvalDistr(model *cov, model **newmodel);
int init_EvalDistr(model *cov, gen_storage *S);
// void do_EvalDistr(model *cov, gen_storage *S);

void RFget(double *x, model *cov, double *v);
int SearchParam(model *cov, get_storage *s) ;
int check_RFget(model *cov) ;
void range_RFget(model *cov, range_type* range);
int struct_RFget(model *cov, model **newmodel);


void Variogram(double *x, model *cov, double *value) ;
void Pseudovariogram(double *x, model *cov, double *value) ;
int check_vario(model *cov);
int struct_variogram(model *cov, model **newmodel);

void Dummy(double *x, model *cov, double *value);
int check_dummy(model *cov);
int struct_dummy(model *cov, model **newmodel);


//-----------------------------------------------------------------
// unsorted

#endif
