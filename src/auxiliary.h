
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


#ifndef AUXILIARY_H
#define AUXILIARY_H 1

#include "auxiliary2.h"

bool is_diag(double *aniso, int dim);

void Abbreviate(char *old, char *abbr);


double getMinimalAbsEigenValue(double *Aniso, int dim);
double getDet(double *Aniso, int dim); // arbitrary matrix


int addressbits(void *addr);
double SurfaceSphere(int d, double r);
double VolumeBall(int d, double r);

double intpow(double x, int p);


void analyse_matrix(double *aniso, int row, int col,
		    bool *diag, bool *quasidiag, int *idx,
		    bool *semiseparatelast, bool *separatelast);
			
double *EinheitsMatrix(int dim);

double incomplete_gamma(double start, double end, double s);


SEXP Mat(double* V, int row, int col);
SEXP Mat_t(double* V, int row, int col);



#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))



#define ATAN2 atan2
#define COSH cosh
#define SINH sinh
#define TANH tanh
#define FMIN fmin2
#define FMAX fmax2

#ifdef EVEN_NOT_USED_ON_SCHLATHERS_MACHINE
#define EXP2 exp2
#define LOGB logb
#define LOG_BASE2 log2
#define CBRT cbrt
#define FMOD fmod
#define REMAINDER remainder
#define FDIM fdim

#else
//#define EXPM1(X) (EXP(X) - 1.0)
#define LOGB(X) FLOOR(LOG(X) *  INVLOG2) //  EXP(X) - 1.0
#define EXP2(X) POW(2.0, X)
#define LOG_BASE2(X) (LOG(X) * INVLOG2)
#define CBRT(X) POW(X, ONETHIRD)
#define FMOD(X, Y) ((double) ((X) - (Y) * (int) ((X) / (Y))))
#define REMAINDER(X, Y) ((double) ((X) - (Y) * fround((X) / (Y), 0))) // ROUND(
#define FDIM(X, Y) FMAX((double) (X) - (Y), 0.0)
#endif

#define ERRSYS(N) ERR1("Function '%.50s' is not available on this system.", N)

#endif /* AUXILIARY_H */




