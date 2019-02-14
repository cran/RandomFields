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


#ifndef RF_INIT_H
#define RF_INIT_H 1

#include "AutoRandomFields.h"
#include "kleinkram.h"

extern int GENERALISEDCAUCHY, STABLE, BROWNIAN, CAUCHY,  GENNSST_INTERN,
  GAUSS, NUGGET, PLUS, MLEPLUS, BALL, MULT, CONST, BIND,
  DISTRIBUTION,  DETERM_DISTR, GAUSS_DISTR, SETPARAM, SET_DISTR, PROD, ANGLE,
  COVMATRIX, RFGET, COVFCTN,LIKELIHOOD_CALL, LINEARPART_CALL, MATRIX,
  PREDICT_CALL,
  DOLLAR, LASTDOLLAR, 
  POWER_DOLLAR,  MLE_ENDCOV,  SPLIT, SCATTER, MCMC_PGS,
  MCMC,
  FIRST_PLANE, LAST_PLANE, 
  FIRST_TRAFO, LAST_TRAFO, EARTHKM2CART, EARTHMILES2CART,
  EARTHKM2GNOMONIC, EARTHMILES2GNOMONIC,
  EARTHKM2ORTHOGRAPHIC, EARTHMILES2ORTHOGRAPHIC, 
  NATSC_USER, NATSC_INTERN, TREND,
  LOC, SET_DISTR, SCALESPHERICAL, TRUNCSUPPORT, BROWNRESNICK, 
  STROKORB_MONO, STROKORB_BALL_INNER, POLYGON, RECTANGULAR,
  IDCOORD, MULT_INVERSE,
  SHAPESTP, SHAPEAVE, SPHERICAL, UNIF, MPPPLUS, ZHOU, BALLANI,
  STATIONARY_SHAPE, STANDARD_SHAPE, TRAFO, TRAFOPROC, PROJMODEL, COVARIATE,
// DENSITY, 
  VARIOGRAM_CALL, SIMULATE,
  BRORIGINAL_USER, BRMIXED_USER, BRSHIFTED_USER, BRNORMED,  
  ARCSQRT_DISTR,SHAPEPOW, POW, SCATTER, GNEITING_INTERN,
  SCALEMODEL, BUBBLE,
  
  BRORIGINAL_INTERN, BRMIXED_INTERN, BRSHIFTED_INTERN,
  MISSING_COV, NULL_MODEL, TBM_OP, USER, 
  DOLLAR_PROC, PLUS_PROC,  VARIOGRAM2COV,
  MPPPLUS_PROC, MULT_PROC, // TREND_PROC,
  BINARYPROC, BROWNRESNICKPROC, GAUSSPROC, POISSONPROC,
  SCHLATHERPROC, SMITHPROC, CHI2PROC, TPROC, EXTREMALTPROC, TREND_PROC,
  PROD_PROC,
  NUGGET_USER,  NUGGET_INTERN,
  CIRCEMBED, CUTOFF, STEIN, SPECTRAL_PROC_USER, SPECTRAL_PROC_INTERN, 
  DIRECT, SEQUENTIAL, 
  AVERAGE_USER, AVERAGE_INTERN, HYPERPLANE_USER, HYPERPLANE_INTERN,
  RANDOMCOIN_USER, CE_CUTOFFPROC_USER, CE_CUTOFFPROC_INTERN, 
  CE_INTRINPROC_USER, CE_INTRINPROC_INTERN, TBM_PROC_USER, TBM_PROC_INTERN, 
  SPECIFIC, SELECTNR, M_PROC,
  BRSHIFTED, BRMIXED, BRORIGINAL, EXTREMALGAUSSIAN, RANDOMSIGN,
  TBM2NR, VECTOR;

#define FIRSTGATTER 0
#define ISO2ISO 1
#define SP2SP 2
#define SP2ISO 3
#define S2ISO 4
#define S2S 5
#define SId 6
#define S2SP 7
#define E2EIso 8
#define E2E 9
#define E2SphIso 10
#define E2Sph 11
#define Sph2SphIso 12
#define Sph2Sph 13
#define LASTGATTER Sph2Sph


typedef enum ptwise_type {pt_posdef, pt_indef, pt_negdef, 
			  pt_zero, pt_paramdep, pt_submodeldep,
			  pt_undefined, pt_unknown, pt_mismatch}
  ptwise_type;
#define last_pt_definite  pt_zero


typedef int pref_type[Forbidden+1];
extern int currentNrCov;
extern pref_type PREF_ALL, PREF_NOTHING, PREF_TREND, PREF_AUX, PREF_MATHDEF;
extern double ONE;
extern char *FREEVARIABLE;
extern name_type STANDARDPARAM, STANDARDSUB;
extern const char *METHOD_NAMES[Forbidden+1],  
  *CAT_TYPE_NAMES[OtherType + 1],
  *REG_NAMES[MODEL_MAX+1],
  **LIST_OF_NAMES[nNamesOfNames],
  *POSITIVITY_NAMES[pt_mismatch + 1];


void includeOtherModels();
void includeCovModels();



#endif /* RF_INIT_H */
