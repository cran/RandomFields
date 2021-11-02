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

#include <Basic_utils.h>
#include "AutoRandomFields.h"
#include "basic.h"





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


#ifdef SCHLATHERS_MACHINE
#define INITCORES 4
#define DO_TESTS true
#else
#define INITCORES 1
#define DO_TESTS false
#endif

#endif


char pch[nr_modes] =         {'\0',  '\0',  '\0',  '.',   '.',   '.',   '.'};
#include "Options.h"


bool 
  allowdistance0[nr_modes] = {true, false, false, false, false, false, false},
  skipchecks[nr_modes] =  {true, false, false, false, false, false, false},
  ce_force[nr_modes] =       {true, true, true, false, false, false, false},
  ce_dependent[nr_modes] =   {true, true, true, false, false, false, false},
  sp_grid[nr_modes] =           {true, true, true, true, true, false, false},
  fit_reoptimise[nr_modes] = {false, false, false, false, true, true, true},
  fit_ratiotest_approx[nr_modes]={true, true, true, true, false, false, false},
  fit_cross_refit[nr_modes] = {false, false, false, false, true, true, true}
  ;
int locmaxn[nr_modes] =      {3000, 4000, 6000, 8000, 9000, 10000000, 10000000},
  ce_trials[nr_modes] =      {1,    2,    3,    3,    4,    5,   6},
  spectral_lines[nr_modes] = {300, 500, 1500,  2500, 3500, 5000, 10000},
  tbm_lines[nr_modes] =      {30,   40,   50,   60,   80,   100,   200},
  mpp_n_estim_E[nr_modes] =  {10000, 10000, 20000,50000,80000,100000,1000000},
  hyper_superpos[nr_modes] = {200, 300, 450, 700, 1000, 2000, 5000},
  maxclique[nr_modes] =      {1000, 2000, 3500, 5000, 6000, 7000, 10000000},
  fit_critical[nr_modes] =   {-1,  -1,   -1,   0,   1,   2,   3}, 
  fit_ncrit[nr_modes] =      { 2,   4,    5,   5,   5,  10,  20}, 
  //  fit_optim_var[nr_modes] =  { 2,   2,    2,   2,   1,   1,   1} ,
  fit_split[nr_modes] =      {10000,  5,  4, 4, 4, 3, 3},
  sequ_back_initial[nr_modes]={3,   3,    5,  10,  20,  30,  40}
  ;

usr_bool exactness[nr_modes] = {False, False, False, Nan, True, True, True};

double 
  ce_tolRe[nr_modes] =       {-1e2,  -1e1,  -1e-3, -1e-7,  -1e-14, 0, 0},
  ce_tolIm[nr_modes] =       {1e2,  1e1,  1e-1,  1e-3,   1e-7, 0, 0},
  ce_approx_step[nr_modes] = {1.0,  1.0,   RF_NA, RF_NA, RF_NA, 0.0, 0.0},
#ifdef SCHLATHERS_MACHINE
  svd_tol[nr_modes] =     {1e-5, 1e-6, 1e-7, 1e-8, 1e-10, 1e-12, 0},
#else
  svd_tol[nr_modes] =     {0, 0, 0, 0, 0, 0, 0},
#endif
  nugget_tol[nr_modes] =     {1e-8, 1e-8, 1e-12, 0,     0,     0,     0},
  tbm_linefactor[nr_modes] = {1.5,  1.5,  1.7,   2.0,  3.0,  5.0,  6.0},
  mpp_intensity[nr_modes] =  {50, 50, 80, 100, 200, 500, 1000},
  mpp_zero[nr_modes] =       {1e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7},
  max_max_gauss[nr_modes] =  {2,    2,    3,    4,    5,   6,   6},
  fit_pgtol[nr_modes]       ={10, 1e-2, 1e-3, 1e-4, 1e-4, 1e-6, 0},
  fit_pgtol_recall[nr_modes]={10, 1e-1, 1e-2, 1e-3, 1e-3, 1e-4, 0},
  fit_factr[nr_modes]       ={1e14, 1e13, 1e12, 1e11, 1e11, 1e9, 1e7},
  fit_factr_recall[nr_modes]={1e15, 1e14, 1e13, 1e12, 1e12, 1e10, 1e8}
   ;


#ifndef INITCORES
  #ifdef SCHLATHERS_MACHINE
    #define INITCORES 4
    #define DO_TESTS true
  #else
    #define INITCORES 1
    #define DO_TESTS false
  #endif
#endif


int PL = 1,
  CORES = INITCORES, // FINITE 
  NS = NAT_SCALE;
globalparam GLOBAL = {
  general_START,
  gauss_START,
  krige_START, 
  ce_START,
  spectral_START,
   tbm_START,
  direct_START,
  sequ_START,
  ave_START,
  nugget_START,
  mpp_START,

  hyper_START,
  extreme_START,
  br_START,
  distr_START,
  fit_START, 
  
  empvario_START, 
  gui_START,  
  graphics_START,
  register_START,
  internal_START,
  coords_START, 
  special_START
  
};

char ERRMSG[LENERRMSG], MSG[LENERRMSG],
#ifndef ERROR_LOC
  ERROR_LOC[LENERRMSG],
#endif  
#ifndef ERRORSTRING  
  ERRORSTRING[LENERRMSG],
#endif  
  MSG2[LENERRMSG];
