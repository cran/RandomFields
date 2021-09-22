/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields 

 Copyright (C) 2001 -- 2017 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
RO
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "RF.h"



SEXP countelements(SEXP Idx, SEXP N, SEXP Totparts) {
  int i, *boxes,
    *idx = INTEGER(Idx),
    nbox = INTEGER(Totparts)[0],
    n = INTEGER(N)[0]; 
  SEXP Boxes;

  PROTECT(Boxes = allocVector(INTSXP, nbox));
  boxes = INTEGER(Boxes);
  for (i=0; i<nbox; i++) boxes[i]= 0;

  for (i=0; i<n; i++) {
    boxes[idx[i]]++;
  }

  UNPROTECT(1);
  return Boxes;
}


SEXP countneighbours(SEXP Xdim, SEXP Parts, SEXP Squarelength,
		     SEXP Cumgridlen, SEXP Boxes, SEXP MAXN) {
  int d,  totcumlen, relstart, y, e,
    dim = INTEGER(Xdim)[0],
    *nb = (int*) MALLOC(sizeof(int) * dim), 
    *loc = (int*) MALLOC(sizeof(int) * dim), 
    x = 0,
    sum = 0,
    sl = INTEGER(Squarelength)[0],
      boundary = (sl - 1) / 2,
    //   total = cumgridlen[dim],
    maxn = INTEGER(MAXN)[0],
    *parts = INTEGER(Parts),
    *cumgridlen = INTEGER(Cumgridlen),
    nboxes = length(Boxes),
    *boxes = INTEGER(Boxes);
  SEXP Neighbours;
  PROTECT(Neighbours = allocVector(INTSXP, nboxes));
  int *neighbours = INTEGER(Neighbours);
 
  totcumlen = 0;
  for (d=0; d<dim; d++) {
    loc[d] = -boundary; 
    nb[d] = 0;
    totcumlen += cumgridlen[d];
    //print("%d %d %d\n", d, totcumlen, cumgridlen[d]);
  }

  relstart = totcumlen * boundary;

  d = 0;
  while(d < dim) {
    y = x - relstart;
    //    print("%d\n", x);
    neighbours[x] = 0; 
    sum = e = 0;
    while(e < dim) {
      bool inside = true;
      int j;
      for (j=0; j<dim; j++) {
	double abs = loc[j] + nb[j];
	if (abs < 0 || abs>=parts[j]) {inside=false; break;}
      }
      if (inside) {
	//     	print("inside %d (%d, %d)\n", y, loc[0], loc[1]);
	sum += boxes[y];
	neighbours[x]++;
      }
      e = 0;	
      loc[e]++;
      y++;
      while (loc[e] > boundary) {
	loc[e] = -boundary; 
	y -= cumgridlen[e] * sl;
	if (++e >= dim) break;
	loc[e]++;
	y += cumgridlen[e];
      }
      //print("e=%d\n", e);
    }
    //     print("sum=%d maxn=%d (%d %d) %d %d parts=%d  %d; %d %dl nei=%d\n",
    // 	    sum, maxn, nb[0], nb[1], totcumlen, relstart, parts[0], parts[1],
	      // 	   cumgridlen[0], cumgridlen[1], neighbours[x]);
     //     assert(false);
    if (sum > maxn) {
      Neighbours = NILSXP;
      goto ErrorHandling;
    }

    d = 0;					
    nb[d]++;
    x++;
    while (nb[d] >= parts[d]) {
      nb[d] = 0;
      if (++d >= dim) break;
      nb[d]++;
    }
    // print("d=%d (%d %d %d) [%d %d %d]\n", d, nb[0], nb[1], nb[2], 
    //	   parts[0], parts[1], parts[2]);
    //  assert(false);
  }

 ErrorHandling:
  UNPROTECT(1);
  FREE(nb);
  FREE(loc);
  return Neighbours;
}


SEXP getelements(SEXP Idx, SEXP Xdim, SEXP N, SEXP Cumgridlen, SEXP Boxes) {
  int i, err = NOERROR, 
    *idx = INTEGER(Idx),
    dim = INTEGER(Xdim)[0],
    n = INTEGER(N)[0],
    *cumgridlen = INTEGER(Cumgridlen),
    *boxes = INTEGER(Boxes),
    total = cumgridlen[dim],
    *count = NULL,
    **elm = NULL;
  SEXP subel = R_NilValue;
   
  if ((elm = (int **) MALLOC(sizeof(int*) * total)) == NULL ||
      (count = (int*) MALLOC(sizeof(int) * total)) == NULL) {
    err = ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  for (i=0; i<total; i++) elm[i] = NULL;
 
  for (i=0; i<total; i++) {
    //   print("%d %d l=%d \n", i, total, boxes[i]); 
    if ((elm[i] = (int *) MALLOC(sizeof(int) * boxes[i])) == NULL) {
      err = ERRORMEMORYALLOCATION;
      goto ErrorHandling;
    }
    count[i] = 0;
  }

  for (i=0; i<n; i++) {
    int k=idx[i];
    elm[k][(count[k])++] = i + 1;
  }

  PROTECT(subel = allocVector(VECSXP, total));
  
  for (i=0; i<total; i++) {
    //  print("%d %d %d\n", i, total, boxes[i]);
    SEXP el;
    int k,
      end = boxes[i];
   
    PROTECT(el = allocVector(INTSXP, boxes[i]));

    for (k=0; k<end; k++) INTEGER(el)[k] = elm[i][k];
    
    SET_VECTOR_ELT(subel, i, el);
    UNPROTECT(1);
  }

  UNPROTECT(1);
 
 ErrorHandling :
  if (elm != NULL) {
    for (i=0; i<total; i++) FREE(elm[i]);
    UNCONDFREE(elm);
  }
  FREE(count);
 
  if (err!=NOERROR) XERR(err);

  return subel; 
}



SEXP getneighbours(SEXP Xdim, SEXP Parts, SEXP Squarelength, 
		   SEXP Cumgridlen, SEXP Neighbours) {
  int i, d, sum, totcumlen, relstart, x, y, e,
    err = NOERROR,
    dim = INTEGER(Xdim)[0],    
    *nb = (int*) MALLOC(sizeof(int) * dim), 
    *loc = (int*) MALLOC(sizeof(int) * dim), 
    *parts = INTEGER(Parts),
    sl = INTEGER(Squarelength)[0],
    *cumgridlen = INTEGER(Cumgridlen),
    *neighbours = INTEGER(Neighbours),
    boundary = (sl - 1) / 2,
    total = cumgridlen[dim],
    ** neighb = NULL;
  SEXP subnei = R_NilValue;

  if ( (neighb = (int **) MALLOC(sizeof(int*) * total)) == NULL) {
    err =  ERRORMEMORYALLOCATION;
    goto ErrorHandling;
  }
  // print("OK %d\n", total);

  for (i=0; i<total; i++) neighb[i] = NULL;
 
  for (i=0; i<total; i++) {
    //    print("%d %d l=%d %d\n", i, total, neighbours[i]); 
    // R_CheckUserInterrupt();
    if ((neighb[i] = (int *) MALLOC(sizeof(int) * neighbours[i])) == NULL) {
      err =  ERRORMEMORYALLOCATION;
      goto ErrorHandling;
    }
  }

  x = 0;
  totcumlen = 0;
  for (d=0; d<dim; d++) {
    loc[d] = -boundary; 
    nb[d] = 0;
    totcumlen += cumgridlen[d];
  }

  relstart = totcumlen * boundary;

  d = 0;
  while(d < dim) {
    y = x - relstart;
    sum = e = 0;
    while(e < dim) {
      bool inside = true;
      int j;
      for (j=0; j<dim; j++) {
	double abs = loc[j] + nb[j];
	if (abs < 0 || abs>=parts[j]) {inside=false; break;}
      }
      if (inside) {
	neighb[x][sum] = y + 1; 
	sum++;
      }

      e = 0;	
      loc[e]++;
      y++;
      while (loc[e] > boundary) {
	loc[e] = -boundary; 
	y -= cumgridlen[e] * sl;
	if (++e >= dim) break;
	loc[e]++;
	y += cumgridlen[e];
      }
    }

    d = 0;					
    nb[d]++;
    x++;
    while (nb[d] >= parts[d]) {
      nb[d] = 0;
      if (++d >= dim) break;
      nb[d]++;
    }
  }

  // print("heres\n");

  PROTECT(subnei = allocVector(VECSXP, total));
  
  for (i=0; i<total; i++) {
    //    print("%d %d %d %d\n", i, total, neighbours[i]);
    R_CheckUserInterrupt();
    SEXP nei;
    PROTECT(nei = allocVector(INTSXP, neighbours[i]));

    int k,
      end = neighbours[i];
    for (k=0; k<end; k++) INTEGER(nei)[k] = neighb[i][k];
    
    SET_VECTOR_ELT(subnei, i, nei);
    UNPROTECT(1);
  }

  UNPROTECT(1);
 

 ErrorHandling :
  FREE(loc);
  FREE(nb);
  if (neighb != NULL) {
    for (i=0; i<total; i++) FREE(neighb[i]);
    UNCONDFREE(neighb);
  }
  if (err!=NOERROR) XERR(err);

  return subnei;
}

