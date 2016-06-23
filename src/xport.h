

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


#ifndef RFxport_H
#define RFxport_H 1

#include <Options_utils.h>

typedef void (*solve_DELETE_type)(solve_storage **S); 
typedef void (*solve_NULL_type)(solve_storage *x); 
typedef int (*solvePosDef_type)(double *, int, bool, double *, int, double *, 
		solve_storage *);

typedef int (*sqrtPosDef_type)(double*, int, solve_storage*);
typedef int (*sqrtRHS_type)(solve_storage*, double*, double*);

typedef void (*getErrorString_type)(char errorstring[MAXERRORSTRING]);
typedef void (*setErrorLoc_type)(char errorloc[nErrorLoc]);
typedef double (*I0mL0_type)(double x);
typedef int (*invertMatrix_type)(double *M, int size);
typedef void (*getUtilsParam_type)(utilsparam**); 
typedef void (*attachRFoptions_type)(const char **, int,  const char ***, int *,
			  setparameterfct, finalsetparameterfct, 
			  getparameterfct);
typedef void (*detachRFoptions_type)(const char **, int);
typedef void (*relaxUnknownRFoption_type)(bool);

extern solve_DELETE_type Ext_solve_DELETE;
extern solve_NULL_type Ext_solve_NULL;
extern solvePosDef_type Ext_solvePosDef;
extern sqrtPosDef_type Ext_sqrtPosDef;
extern sqrtRHS_type Ext_sqrtRHS;
extern getErrorString_type Ext_getErrorString;
extern setErrorLoc_type Ext_setErrorLoc;
extern I0mL0_type Ext_I0mL0;
extern invertMatrix_type Ext_invertMatrix;
extern getUtilsParam_type Ext_getUtilsParam;
extern attachRFoptions_type Ext_attachRFoptions;
extern detachRFoptions_type Ext_detachRFoptions;
extern relaxUnknownRFoption_type Ext_relaxUnknownRFoption;

void includeXport();
extern utilsparam* GLOBAL_UTILS;
#endif
