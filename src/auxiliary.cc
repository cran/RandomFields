//#define DEBUG 1


/* 
 Authors
 Martin Schlather, schlath@hsu-hh.de

 Collection of auxiliary functions

 Copyright (C) 2001 -- 2006 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#include <math.h>
#include <unistd.h>
#include <assert.h>
#include "auxiliary.h"

bool is_diag(double *aniso, int dim) {
  int diag = dim + 1, size = dim * dim, i;
  bool notdiag=false;
  for (i=0; i<size; i++) {
    if (i % diag != 0 && (notdiag = aniso[i] != 0.0)) break;
  }
  return !notdiag;
};

void I0ML0(double *x, int *n) {
  int i;
  for (i=0; i<*n; i++) x[i] = I0mL0(x[i]);
} 

double I0mL0(double x){
  /* Bessel I_0 - Struve L_0 for non-negative arguments x */
  /* see J. MacLeod, Chebyshev expansions for modified {S}truve and 
     related functions, Mathematics of Computation, 60, 735-747, 1993 */
  static double g2[24] = {
	0.52468736791485599138e0,
	-0.35612460699650586196e0,
	0.20487202864009927687e0,
	-0.10418640520402693629e0,
	0.4634211095548429228e-1,
	-0.1790587192403498630e-1,
	0.597968695481143177e-2,
	-0.171777547693565429e-2,
	0.42204654469171422e-3,
	-0.8796178522094125e-4,
	0.1535434234869223e-4,
	-0.219780769584743e-5,
	0.24820683936666e-6,
	-0.2032706035607e-7,
	0.90984198421e-9,
	0.2561793929e-10,
	-0.710609790e-11,
	0.32716960e-12,
	0.2300215e-13,
	-0.292109e-14,
	-0.3566e-16,
	0.1832e-16,
	-0.10e-18,
	-0.11e-18
  };
  static double g3[24] = {
	2.00326510241160643125e0,
	0.195206851576492081e-2,
	0.38239523569908328e-3,
	0.7534280817054436e-4,
	0.1495957655897078e-4,
	0.299940531210557e-5,
	0.60769604822459e-6,
	0.12399495544506e-6,
	0.2523262552649e-7,
	0.504634857332e-8,
	0.97913236230e-9,
	0.18389115241e-9,
	0.3376309278e-10,
	0.611179703e-11,
	0.108472972e-11,
	0.18861271e-12,
	0.3280345e-13,
	0.565647e-14,
	0.93300e-15,
	0.15881e-15,
	0.2791e-16,
	0.389e-17,
	0.70e-18,
	0.16e-18
  };
  double r, x2, ac;
  int i;
  
  if (x < 0.0) {return RF_NAN;}
  if (x < 16.0) {
    r = 0.5 * g2[0];
    ac = acos((6.0 * x - 40.0) / (x + 40.0));
    for (i=1; i<24; i++) {
      r += g2[i] * cos(i * ac);
    }
  } else {
    r = 0.5 * g3[0];
    x2 = x * x;
    ac = acos((800.0 - x2) / (288.0 + x2));
    for (i=1; i<24; i++) {
      r += g3[i] * cos(i * ac);
    }
    r *= T_PI /* 2/pi */ / x;
  }
  return r;
}

void StruveH(double *x, double *nu) {*x=struve(*x, *nu, -1.0, false);}
void StruveL(double *x, double *nu, int * expScaled) 
{
  *x=struve(*x, *nu, 1.0, (bool) *expScaled);
}
double struve(double x, double nu, double factor_sign, bool expscaled)
{ 
  double sign, res, epsilon=1e-20;
  double dummy, logx, x1, x2;
  if ((x==0.0) && (nu>-1.0)) return 0.0;
  if (x<=0) return RF_NAN; // not programmed yet
  logx = log(0.5 * x);
  x1=1.5;   
  x2=nu+1.5;   
  sign=1.0;
  if (x2 > 0.0) { 
    dummy = (nu + 1.0) * logx - lgammafn(x1) - lgammafn(x2);
    if (expscaled) dummy -= x;
    res = exp(dummy);
  } else {
    if ( (double) ((int) (x1-0.5)) != x1-0.5 ) return RF_NAN;
    res=pow(0.5 * x, nu + 1.0) / (gammafn(x1) * gammafn(x2));
    if (expscaled) res *= exp(-x);
    if ((dummy= res) <0) {
      dummy = -dummy;
      sign = -1.0;
    }
    dummy = log(dummy);
  }
  logx *= 2.0;
  do {
    if (x2<0) { sign = -sign; }    
    dummy += logx - log(x1) - log(fabs(x2));
    res +=  sign * exp(dummy);
    x1 += 1.0;
    x2 += 1.0;
    sign = factor_sign * sign; 
  } while (exp(dummy) > fabs(res) * epsilon);
  return res;
}

void vectordist(double *v, int *Dim, double *dist, int *diag){
  int m, n, d, l, dim, r, lr, dr, add;
  add = (*diag==0) ? 1 : 0;
  l = Dim[0];
  dim = Dim[1] * l;
  lr = (l * (l + 1 - 2 * add)) / 2;
  for (r=0, m=0; m<l; m++) { // if add==1 loop is one to large, but
      // but doesn't matter
    for (n=m+add; n<l; n++, r++) {
      for (d=0, dr=0; d<dim; d+=l, dr+=lr) {
	  dist[r + dr] = v[n + d] - v[m + d];
      }
    }
  }
} 

static int ORDERDIM;
static double *ORDERD;
bool smaller(int i, int j)
{
  double *x, *y;
  int d;
  x = ORDERD + i * ORDERDIM;
  y = ORDERD + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++)
     if (x[d] != y[d]) return x[d] < y[d];
  return false;
}

bool greater(int i, int j)
{
  double *x, *y;
  int d;
  x = ORDERD + i * ORDERDIM;
  y = ORDERD + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++)
    if (x[d] != y[d]) return x[d] > y[d];
  return false;
}

void order(int *pos, int start, int end) {
  int randpos, pivot, left, right, pivotpos, swap;
  
  if( start < end ) {   
    //GetRNGstate();randpos = start + (int) (UNIFORM_RANDOM * (end-start+1)); PutRNGstate();
    randpos = (int) (0.5 * (start + end));
    pivot = pos[randpos];
    pos[randpos] = pos[start];
    pos[start] = pivot;
    
    pivotpos=start; 
    left = start;
    right=end+1;   
    while (left < right) {      
      while (++left < right && smaller(pos[left], pivot)) pivotpos++;
      while (--right > left && greater(pos[right], pivot));      
      if (left < right) {
	swap=pos[left]; pos[left]=pos[right]; pos[right]=swap;
	pivotpos++;
      }
    }
    pos[start] = pos[pivotpos];
    pos[pivotpos] = pivot;
    order(pos, start, pivotpos-1 );
    order(pos, pivotpos + 1, end );
  }
}

void ordering(double *d, int len, int dim, int *pos) 
{
  /* quicksort algorithm, slightly modified:
     does not sort the data, but d[pos] will be ordered 
     NOTE: pos must have the values 0,1,2,...,start-end !
     (orderdouble is a kind of sorting pos according to
     the variable d)
  */  
  int i;
  for (i=0; i<len;i++) pos[i]=i;
  ORDERD = d;
  ORDERDIM = dim;
  order(pos, 0, len-1);
}


void Ordering(double *d, int *len, int *dim, int *pos) 
{
  int i;
  for (i=0; i<*len; i++) pos[i]=i;
  ORDERD = d;
  ORDERDIM = *dim;
  order(pos, 0, *len-1 );
}

