

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





#ifndef GSL_VS_R_H
#define GSL_VS_R_H 1

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <errno.h>
#include <R_ext/Complex.h>

#ifdef SCHLATHERS_MACHINE
#undef SCHLATHERS_MACHINE
#endif

#ifdef RANDOMFIELDS_DEBUGGING
#undef RANDOMFIELDS_DEBUGGING
#endif

#ifdef showfree
#undef showfree
#endif

#ifdef DOPRINT
#undef DOPRINT
#endif



#define MAXNRCOVFCTS 300
#define MAXUNITSCHAR 10
#define MAXINVERSIONS 2


// Formerly in <R_ext/Applic.h>LinkedTo: 
void fft_factor_(int n, int *pmaxf, int *pmaxp);
Rboolean fft_work_(double *a, double *b, int nseg, int n, int nspn,
		  int isn, double *work, int *iwork);/* TRUE: success */




//
// 
// 
// 



#endif /* GSL_VS_R_H */


