/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 wrapper for Havard Rue's GMRF library

 Copyright (C) 2001 -- 2011 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

// was passiert wenn geo-coef nicht verfuegbar?
// wie wird *problem geloescht ??



#include "RF.h"


typedef struct markov_storage {
// typedef cannot be moved to RF.h since
// it requires GMRFLib/GMRFLib.h and the inclusion
// in RF.h leads to conflicts with blas.h definitions
// used in circulant.cc
int rowcol[2];
} markov_storage;


void Markov_destruct(void **S) 
{ 
  if (*S!=NULL) {
   free(*S);   
    *S = NULL;
  }
}


void SetParamMarkov(int *action,int *neighbours, double *precision,
		    int *cyclic, int *maxmem) {
   markov_param *lp = &(GLOBAL.markov);
  if (*action) {
    lp->neighbours = *neighbours; 
    if (lp->neighbours < 2) {
	if (PL>0) PRINTF("minimal neighbourhood is 2");
	lp->neighbours = 2;
    } 
    if (lp->neighbours > 3) {
	if (PL>0) PRINTF("maximal neighbourhood is 3");
	lp->neighbours = 3;
    }
    lp->precision = *precision;
    lp->cyclic = *cyclic;
    lp->maxmem = *maxmem;
  } else {
    *neighbours = lp->neighbours; 
    *precision = lp->precision;
    *cyclic = lp->cyclic;
    *maxmem = lp->maxmem;
  }
}



int init_markov(method_type *meth){ 
  return ERRORFAILED;
}

    
 
void do_markov(method_type *meth, res_type *res ) {
  error("Markov not programmed anymore");
}

