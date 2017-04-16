
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


#ifndef AUXILIARY22_H
#define AUXILIARY22_H 1

#include "basic.h"
#include "AutoRandomFields.h"

//double I0mL0(double x);

#ifdef __cplusplus
extern "C" {
#endif 
  SEXP vectordist(SEXP V, SEXP diag); 
#ifdef __cplusplus
}
#endif


#endif /* AUXILIARY22_H */




