
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
#include "def.h"
#include "RandomFields.h"
#include "auxiliary2.h"
#include <Basic_utils.h>

#define none 0
static R_NativePrimitiveArgType
  one_int[] = { INTSXP },
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

#define CDEF(name, n, type) {#name, (DL_FUNC) &name, n, type}
static const R_CMethodDef cMethods[]  = {  
  CDEF(GetAttr, 14, attr_arg),
  CDEF(GetModelName, 3, inttwochar),
  CDEF(GetModelNr, 2, charint),
  CDEF(GetCurrentNrOfModels, 1, one_int),
  CDEF(GetNrParameters, 2, two_int),
  CDEF(PrintModelList, 3, three_int),
  CDEF(PutValuesAtNA, 2, intdouble),
  CDEF(PutValuesAtNAnoInit, 2, intdouble),
  CDEF(expliciteDollarMLE, 2, intdouble),
  CDEF(MultiDimRange, 3, intintdouble),
  CDEF(GetModelRegister, 2, charint),
  CDEF(ResetWarnings, 1, one_int),
  CDEF(NoCurrentRegister, 0, none),
  CDEF(GetCurrentRegister, 1, one_int),
  CDEF(PutGlblVar, 2, intdouble),
  CDEF(attachRFoptionsRandomFields, 1, one_int),
  CDEF(detachRFoptionsRandomFields, 0, none),
  CDEF(RelaxUnknownRFoptions, 1, one_int),
  // CDEF(attachRFoptionsUtils, 0, NULL, NULL},
  // CDEF(detachRFoptionsUtils, 0, NULL, NULL},
  {NULL, NULL, 0, none}
};



#define CALLDEF_DO(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  CALLDEF_DO(GetParameterNames, 1),
  CALLDEF_DO(GetSubNames, 1),
  CALLDEF_DO(GetAllModelNames, 1),
  CALLDEF_DO(GetCoordSystem, 3),
  CALLDEF_DO(GetModelInfo, 5),
  CALLDEF_DO(GetModel, 7),
  CALLDEF_DO(Init, 4),
  CALLDEF_DO(EvaluateModel, 2),
  CALLDEF_DO(GetProcessType, 2),
  CALLDEF_DO(empirical, 10),
  CALLDEF_DO(empvarioXT, 13),
  CALLDEF_DO(fftVario3D, 14),
  CALLDEF_DO(boxcounting, 5),
  CALLDEF_DO(detrendedfluc, 5),
  CALLDEF_DO(periodogram, 6),
  CALLDEF_DO(minmax, 5), 
  CALLDEF_DO(VariogramIntern, 1),
  CALLDEF_DO(CovLoc, 6),
  CALLDEF_DO(GetNAPositions, 7),   
  CALLDEF_DO(SetAndGetModelInfo, 10),   
  CALLDEF_DO(SetAndGetModelLikelihood, 4), 
  CALLDEF_DO(Take2ndAtNaOf1st, 8),
  CALLDEF_DO(countelements, 3), 
  CALLDEF_DO(countneighbours, 6), 
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
  CALLDEF_DO(maintainers_machine, 0),
//  CALLDEF_DO(),
  {NULL, NULL, 0}
};


// otherwise clang does not recognize it
#ifdef __cplusplus
extern "C" {
#endif 

void R_init_RandomFields(DllInfo  *dll) {
  R_registerRoutines(dll, cMethods, callMethods, NULL, // .Fortran
		     NULL // extended
		     );
  R_useDynamicSymbols(dll, FALSE); // OK
}


void R_unload_RandomFields(DllInfo *info) {
  // just to avoid warning from compiler
  //int *x = (int*) info; x = x + 0;
  // if (0) { PRINTF("%ld\n", (long int) info);}
  /* Release resources. */
}

 
#ifdef __cplusplus
}
#endif

