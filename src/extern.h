 /* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 simulation of max-stable random fields

 Copyright (C) 2001 -- 2017 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#ifndef randomfields_extern_H
#define randomfields_exterm_H 1

#include "basic.h"
#include "Options.h"

typedef char errorstring_type[MAXERRORSTRING];
typedef char errorloc_type[nErrorLoc];

typedef int Rint;
typedef unsigned int Uint;

extern globalparam GLOBAL;
extern int PL, CORES, NS;
extern char pch[nr_modes];


extern bool 
  allowdistance0[nr_modes] ,
  skipchecks[nr_modes],
  ce_force[nr_modes] ,
  ce_dependent[nr_modes],
  sp_grid[nr_modes] ,
  fit_reoptimise[nr_modes] ,
  fit_ratiotest_approx[nr_modes],
  fit_cross_refit[nr_modes] 
  ;

extern int locmaxn[nr_modes] ,
  ce_trials[nr_modes] ,
  spectral_lines[nr_modes] ,
  tbm_lines[nr_modes] ,
  mpp_n_estim_E[nr_modes],
  hyper_superpos[nr_modes],
  maxclique[nr_modes] ,
  fit_critical[nr_modes] ,
  fit_ncrit[nr_modes] ,
  fit_split[nr_modes] ,
  sequ_back_initial[nr_modes]
  ;

extern usr_bool exactness[nr_modes] ;

extern double 
ce_tolRe[nr_modes] ,
  ce_tolIm[nr_modes],
  ce_approx_step[nr_modes] ,
  svd_tol[nr_modes]  ,
  nugget_tol[nr_modes] ,
  tbm_linefactor[nr_modes] ,
  mpp_intensity[nr_modes] ,
  mpp_zero[nr_modes] ,
  max_max_gauss[nr_modes] ,
  fit_pgtol[nr_modes]     ,
  fit_pgtol_recall[nr_modes],
  fit_factr[nr_modes]       ,
  fit_factr_recall[nr_modes]
   ;



//extern utilsparam GLOBAL;
//extern const char *basic[basicN];
//extern const char * InversionNames[nr_InversionMethods];
//extern const char * solve[solveN];
// extern bool ToFalse[1];
//




#ifndef R_PRINTLEVEL


#ifdef Long
#undef Long
#endif

#ifdef Ulong
#undef Ulong
#endif


#define R_PRINTLEVEL 1
#define C_PRINTLEVEL 1


#define LENMSG 1000
#define LENERRMSG 1000
#define MAXERRORSTRING 1000
#define nErrorLoc 1000
typedef char errorstring_type[MAXERRORSTRING];
typedef char errorloc_type[nErrorLoc];

typedef int Rint;
typedef unsigned int Uint;


#ifdef SCHLATHERS_MACHINE
#define INITCORES 4
#define DO_TESTS true
#else
#define INITCORES 1
#define DO_TESTS false
#endif

#endif


#endif
