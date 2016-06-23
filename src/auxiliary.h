
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


#ifndef AUXILIARY_H
#define AUXILIARY_H 1

#include "basic.h"
#include "AutoRandomFields.h"

//double I0mL0(double x);
extern "C" {
  SEXP vectordist(SEXP V, SEXP diag); 
  void sleepMilli(int *milli);
  //  void I0ML0(double *x, int *n);
  //double struve(double x, double nu,  double factor_Sign, bool expscaled);
  // void StruveH(double *x, double *nu);
  // void StruveL(double *x, double *nu, int * expScaled);
  void ordering(double *d, int len, int dim, int *pos);
  void Ordering(double *d, int *len, int *dim, int *pos);
}
bool is_diag(double *aniso, int dim);

void Abbreviate(char *old, char *abbr);


#endif /* AUXILIARY_H */




