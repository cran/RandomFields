/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by Cholesky or SVD decomposition

 Copyright (C) 2015 -- Martin Schlather, 

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

#include "RF.h"

solve_DELETE_type Ext_solve_DELETE = NULL;
solve_NULL_type Ext_solve_NULL = NULL;
solvePosDef__type Ext_solvePosDef_ = NULL;
getErrorString_type Ext_getErrorString = NULL;
setErrorLoc_type Ext_setErrorLoc = NULL;
I0mL0_type Ext_I0mL0 = NULL;
invertMatrix_type Ext_invertMatrix = NULL;


#include <R_ext/Rdynload.h>

#define pkg "RandomFieldsUtils"
#define CALL(what) Ext_##what = (what##_type) R_GetCCallable(pkg, #what)

void includeXport() {

  CALL(solve_DELETE); 
  CALL(solve_NULL); 
  CALL(solvePosDef_); 
  CALL(getErrorString); 
  CALL(setErrorLoc);
  CALL(I0mL0);
  CALL(invertMatrix);
 /*
   CALL(= 
   (_type) R_GetCCallable("RandomFieldsUtils", "");
  */

} // export C


