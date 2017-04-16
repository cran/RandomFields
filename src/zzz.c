
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- 2017 Martin Schlather, Reinhard Furrer

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

#include <R_ext/Rdynload.h>
#include "RandomFields.h"
#include "auxiliary2.h"


static R_NativePrimitiveArgType one_int[] = { INTSXP },
  intdouble[] = {INTSXP, REALSXP }, 
  intintdouble[] = {INTSXP, INTSXP, REALSXP }, 
  charint[] = { STRSXP, INTSXP}, 
  inttwochar[] = {INTSXP, STRSXP, STRSXP }, 
  two_int[] = { INTSXP, INTSXP },
  three_int[] = { INTSXP, INTSXP, INTSXP},
  attr_arg[] = { INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,  INTSXP,
		 INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,  INTSXP };
  //  static R_NativeArgStyle argin[] = {R_ARG_IN},
  //    argout[] = {R_ARG_OUT},
  //   hostarg[] = {R_ARG_OUT, R_ARG_OUT};
static const R_CMethodDef cMethods[]  = {
  {"GetAttr", (DL_FUNC) &GetAttr, 14, attr_arg},
  {"GetModelName", (DL_FUNC) &GetModelName, 3, inttwochar},
  {"GetModelNr", (DL_FUNC) &GetModelNr, 2, charint},
  {"GetCurrentNrOfModels", (DL_FUNC) &GetCurrentNrOfModels, 2, two_int},
  {"GetNrParameters", (DL_FUNC) &GetNrParameters, 2, two_int},
  {"PrintModelList", (DL_FUNC) &PrintModelList, 3, three_int},
  {"PutValuesAtNA", (DL_FUNC) &PutValuesAtNAnoInit, 2, intdouble},
  {"PutValuesAtNAnoInit", (DL_FUNC) &PutValuesAtNAnoInit, 2, intdouble},
  {"expliciteDollarMLE", (DL_FUNC) &expliciteDollarMLE, 2, intdouble},
  {"MultiDimRange", (DL_FUNC) &MultiDimRange, 3, intintdouble},
  {"GetModelRegister", (DL_FUNC) &GetModelRegister, 2, charint},
   {"ResetWarnings", (DL_FUNC) &ResetWarnings, 1, one_int},
  {"NoCurrentRegister", (DL_FUNC) &NoCurrentRegister, 0},
  {"GetCurrentRegister", (DL_FUNC) &GetCurrentRegister, 1, one_int},
  {"PutGlblVar", (DL_FUNC) &PutGlblVar, 2, intdouble},
  {"DeleteKey", (DL_FUNC) &DeleteKey, 1, one_int},
  {"detachRFoptionsRandomFields", (DL_FUNC) &detachRFoptionsRandomFields, 0},
  {"attachRFoptionsRandomFields", (DL_FUNC) &attachRFoptionsRandomFields, 0},
  {"RelaxUnknownRFoption", (DL_FUNC) &RelaxUnknownRFoption, 1, one_int},
  // {"attachRFoptionsUtils", (DL_FUNC) &attachRFoptionsUtils, 0, NULL, NULL},
  // {"detachRFoptionsUtils", (DL_FUNC) &detachRFoptionsUtils, 0, NULL, NULL},
  {NULL, NULL, 0, NULL}
};



#define CALLDEF_DO(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  CALLDEF_DO(GetParameterNames, 1),
  CALLDEF_DO(GetSubNames, 1),
  CALLDEF_DO(GetAllModelNames, 0),
  CALLDEF_DO(GetCoordSystem, 3),
  CALLDEF_DO(GetExtModelInfo, 4),
  CALLDEF_DO(GetModel, 6),
  CALLDEF_DO(Init, 4),
  CALLDEF_DO(EvaluateModel, 2),
  CALLDEF_DO(GetProcessType, 2),
  CALLDEF_DO(empiricalvariogram, 10),
  CALLDEF_DO(empvarioXT, 13),
  CALLDEF_DO(fftVario3D, 14),
  CALLDEF_DO(boxcounting, 5),
  CALLDEF_DO(detrendedfluc, 5),
  CALLDEF_DO(periodogram, 6),
  CALLDEF_DO(minmax, 5), 
  CALLDEF_DO(VariogramIntern, 1),
  CALLDEF_DO(CovLoc, 6),
  //CALLDEF_DO(GetNAPositions, 5),   
  CALLDEF_DO(SetAndGetModelInfo, 10),   
  CALLDEF_DO(SetAndGetModelLikeli, 3), 
  CALLDEF_DO(Take2ndAtNaOf1st, 8),
  CALLDEF_DO(countelements, 3), 
  CALLDEF_DO(countneighbours, 5), 
  CALLDEF_DO(getelements, 5), 
  CALLDEF_DO(getneighbours, 5), 
  CALLDEF_DO(set_boxcox, 1), 
  CALLDEF_DO(get_boxcox, 0),
  //  CALLDEF_DO(BoxCox_inverse, 2), 
  CALLDEF_DO(BoxCox_trafo, 4),
  CALLDEF_DO(get_logli_residuals, 1),
  CALLDEF_DO(get_likeliinfo, 1),
  CALLDEF_DO(simple_residuals, 1), 
  CALLDEF_DO(get_linearpart, 2),
  CALLDEF_DO(vectordist, 2),
//  CALLDEF_DO(),
  {NULL, NULL, 0}
};


void R_init_RandomFields(DllInfo  *dll) {
  R_registerRoutines(dll, cMethods, callMethods, NULL, // .Fortran
		     NULL // extended
		     );
  R_useDynamicSymbols(dll, FALSE);
}


void R_unload_RandomFields(DllInfo *info) {
  /* Release resources. */
}

