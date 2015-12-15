

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


typedef void (*solve_DELETE_type)(solve_storage **S); 
typedef void (*solve_NULL_type)(solve_storage *x); 
typedef int (*solvePosDef__type)(double *, int, bool, double *, int, double *, 
		solve_storage *, solve_param *,  int);

typedef int (*sqrt__type)(double*, int, solve_storage*, solve_param*, int);
typedef int (*sqrt_RHS__type)(solve_storage*, double*, double*);

typedef void (*getErrorString_type)(char errorstring[MAXERRORSTRING]);
typedef void (*setErrorLoc_type)(char errorloc[nErrorLoc]);
typedef double (*I0mL0_type)(double x);
typedef int (*invertMatrix_type)(double *M, int size);


extern solve_DELETE_type Ext_solve_DELETE;
extern solve_NULL_type Ext_solve_NULL;
extern solvePosDef__type Ext_solvePosDef_;
extern sqrt__type Ext_sqrt_;
extern sqrt_RHS__type Ext_sqrt_RHS_;
extern getErrorString_type Ext_getErrorString;
extern setErrorLoc_type Ext_setErrorLoc;
extern I0mL0_type Ext_I0mL0;
extern invertMatrix_type Ext_invertMatrix;

void includeXport();

#endif
