//#define DEBUG 1


/* 
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Collection of auxiliary functions

 Copyright (C) 2001 -- 2011 Martin Schlather, 

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
 
#include "auxiliary.h"
//#include <curses.h>
#include "RandomFields.h"
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


void distInt(int *X, int*N, int *Genes, double *dist) {
    int i,j, k, di, diff, *x, *y, ve, ho, endfor,
	n = *N,
	nP1 = n + 1,
	genes = *Genes;
 
  x = y = X;
  for (j=0, i=0;  j<n;  i += nP1, j++, y += genes) {
    dist[i] = 0.0;
    endfor = i + (n - j);
    for (ve = i + 1, ho = i + n, x = y + genes; 
         ve < endfor; 
	 ve++, ho += n) {
      for (di=0.0, k=0; k<genes; k++, x++) {
	diff = *x - y[k];
	di += diff * diff;
      }
      dist[ve] = dist[ho] = sqrt((double) di);
    }
  }
}


SEXP GetChar(SEXP N, SEXP Choice, SEXP Shorter, SEXP Beep, SEXP Show) {  
  int i, j, ret,
      start = RF_NAN,
      milli = 500,
      len = LENGTH(Choice), 
      n = INTEGER(N)[0],
      nline = 100;
  bool shorter = LOGICAL(Shorter)[0],
      piep = LOGICAL(Beep)[0],
      show = LOGICAL(Show)[0];
  char choice[256], *zk,
      backspace = 127,
      eof = 'e'; // crl-D, ^D
  SEXP Zk;
  
  zk = (char*) malloc(sizeof(char) * (n+1));
  if (len > 256) len = 256;  
  for (i=0; i<len; i++) {
      choice[i] = CHAR(STRING_ELT(Choice, i))[0];
  }
 
  ret = system("/bin/stty cbreak -echo iuclc"); /* or "stty raw" */
  if (ret < 0) error("GetChar failed.");
  for (j=0; j<n; ) {
   if (j % nline == 0) {
       start = j;
       if (show) {
	   if (j != 0) PRINTF("\n");
	   int m = n-j;
	   if (m > nline) m = nline;
	   for (i=0; i<m; i++) if (i % 20 == 19) PRINTF(":"); else PRINTF(".");
	   for (i=0; i<m; i++) PRINTF("\b");
       }
   }

    zk[j] = getchar();
    if (zk[j] == eof && shorter) break;
    if (zk[j] == backspace && j>start) {
      j--;
      PRINTF("\b \b");
      continue;
    }
    for (i=0; i<len; i++)
      if (zk[j] == choice[i]) break;
    if (i < len) {
	if (show) PRINTF("%c", zk[j], n-1-j);
        j++; 
    } else {
	if (piep) PRINTF("beep does not work.\n"); // beep();
    }
  }
//  beep();
  if (show) {
      PRINTF("\n");
      sleepMilli(&milli);
  }
  ret = system("/bin/stty -cbreak echo -iuclc");
  if (ret < 0) error("GetChar failed.");
  zk[j] = '\0';
  PROTECT(Zk = allocVector(STRSXP, 1));
  SET_STRING_ELT(Zk, 0, mkChar(zk));
  UNPROTECT(1);      
  return Zk;
}

void strcopyN(char *dest, const char *src, int n) {
  n--;
  strncpy(dest, src, n);
  dest[n] = '\0';
}

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
  double sign, value, epsilon=1e-20;
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
    value = exp(dummy);
  } else {
    if ( (double) ((int) (x1-0.5)) != x1-0.5 ) return RF_NAN;
    value=pow(0.5 * x, nu + 1.0) / (gammafn(x1) * gammafn(x2));
    if (expscaled) value *= exp(-x);
    if ((dummy= value) <0) {
      dummy = -dummy;
      sign = -1.0;
    }
    dummy = log(dummy);
  }
  logx *= 2.0;
  do {
    if (x2<0) { sign = -sign; }    
    dummy += logx - log(x1) - log(fabs(x2));
    value +=  sign * exp(dummy);
    x1 += 1.0;
    x2 += 1.0;
    sign = factor_sign * sign; 
  } while (exp(dummy) > fabs(value) * epsilon);
  return value;
}

void vectordist(double *v, int *Dim, double *Dist, int *diag){
  int d, dim, dr;
  double *v1, *v2, *end;
  bool notdiag = (*diag==0);
  dim = Dim[0];
  end = v + Dim[1] * dim; 

//  printf("%d %d %f %f\n", dim , Dim[0], v, end);

  for (dr=0, v1=v; v1<end; v1+=dim) { // if add==1 loop is one to large, but
      // but doesn't matter
    v2 = v1;
    if (notdiag) {
       v2 += dim;
    }
    for (; v2<end; ) {
      for (d=0; d<dim; v2++) {
	Dist[dr++] = v1[d++] - *v2;
      }
    }
  }
} 

static int ORDERDIM;
static double *ORDERD;
static int *ORDERDINT;

typedef bool (*vergleich)(int, int);
vergleich SMALLER=NULL, GREATER=NULL;

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

bool smallerInt(int i, int j)
{
  int *x, *y;
  int d;
  x = ORDERDINT + i * ORDERDIM;
  y = ORDERDINT + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++)
     if (x[d] != y[d]) return x[d] < y[d];
  return false;
}



bool greaterInt(int i, int j)
{
  int *x, *y;
  int d;
  x = ORDERDINT + i * ORDERDIM;
  y = ORDERDINT + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++)
    if (x[d] != y[d]) return x[d] > y[d];
  return false;
}

int is_positive_definite(double *C, int dim) {
    int err,
	bytes = sizeof(double) * dim * dim;
  double *test;
  test = (double*) malloc(bytes);
  memcpy(test, C, bytes);
  F77_CALL(dpofa)(test, &dim, &dim, &err); 
  free(test);
  return(err == 0);
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
      while (++left < right && SMALLER(pos[left], pivot)) pivotpos++;
      while (--right > left && GREATER(pos[right], pivot));      
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
  SMALLER = smaller;
  GREATER = greater;
  order(pos, 0, len-1);
}


void Ordering(double *d, int *len, int *dim, int *pos) 
{
  int i;
  for (i=0; i<*len; i++) pos[i]=i;
  ORDERD = d;
  ORDERDIM = *dim;
  SMALLER = smaller;
  GREATER = greater;
  order(pos, 0, *len-1 );
}


void orderingInt(int *d, int len, int dim, int *pos) 
{
  /* quicksort algorithm, slightly modified:
     does not sort the data, but d[pos] will be ordered 
     NOTE: pos must have the values 0,1,2,...,start-end !
     (orderdouble is a kind of sorting pos according to
     the variable d)
  */  
  int i;
  for (i=0; i<len;i++) pos[i]=i;
  ORDERDINT = d;
  ORDERDIM = dim;
  SMALLER = smallerInt;
  GREATER = greaterInt;
  order(pos, 0, len-1);
}


int xMatch(char *name, char **list, unsigned int llen)  {
  unsigned int ln, nr;
  // == -1 if no matching function is found
  // == -2 if multiple matching fctns are found, without one matching exactly
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)
  // see also GetModelNr !

  nr=0;
  ln=strlen(name);
  while ((nr < llen) && strncmp(name, list[nr], ln)) nr++;
  if (nr < llen) { 
    // a matching method is found. Are there other methods that match?
    unsigned int j; 
    if (ln == strlen(list[nr])) return nr; // exact matching
    j = nr + 1; // if two or more methods have the same name the last one is
    //          taken; stupid here, but nice in GetCovNr
    while (j < llen && strncmp(name, list[j], ln)) j++;
    if (j < llen) {
      if (ln == strlen(list[j])) return j; //exact matching 
      else return -2; // multiple matching
    }
    return nr;
  } else return -1; // unmatched
}


