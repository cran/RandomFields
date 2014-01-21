/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part and error messages

 Copyright (C) 2011 -- 2014 Martin Schlather
tree
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


#include <math.h>
#include "RandomFields.h"

#define G 0
#define M 1
#define F 2
#define S 3

#define vor 0
#define links 1
#define oben 2

#define ARGS L, B, H, lambda, states, pi, doublet, maxdepth, cum_p

/*

// R : Achtung reverse ordering of indices!

states[4][4][3] = {
  // G           M          F          S (vorwaerts, links, oben)
  { { 0, 0, 0}, {0, 0, 0}, {14,16, 0}, {0, 0, 7} }, // G
  { {12, 0,19}, {0, 0,20}, {11, 0, 0}, {9, 0,10} }, // M
  { {13,15,17}, {0, 0, 0}, { 0, 0, 0}, {3, 5, 1} }, // F
  { { 8, 0,18}, {0, 0, 0}, { 4, 6, 2}, {0, 0, 0} }  // S
},


// R : 21 Zeilen, 2 Spalten

  change[21][2] = { {-1, -1},
		    {S, F}, {F, S}, {S, F}, {F, S}, {S, F}, {F, S}, // 6
		    {S, M}, {S, M}, // 8
		    {S, M}, {S, M}, // 10
		    {F, G}, {G, G}, // 12
		    {G, F}, {F, G}, {G, F}, {F, G}, // 16
		    {G, F}, {G, S}, {G, G}, {G, M}
  };

  lambda[0] = 0;
  lambda[1] = lW;
  lambda[2] = lW;
  lambda[3] = lWM;
  lambda[4] = lWP;
  lambda[5] = lW;
  lambda[6] = lW;
  
  lambda[7] = lS;
  lambda[8] = lS;
  lambda[9] = lS;
  lambda[10] = lS;
  lambda[11] = lS;
  lambda[12] = a * lS;
  
  lambda[13] = lD;
  lambda[14] = lD;
  lambda[15] = lD;
  lambda[16] = lD;
  
  lambda[17] = lG; 
  lambda[18] = lG; 
  lambda[19] = lG; 
  lambda[20] = lG;
 
*/


void checkS(int *Z, int L, int B, int H, bool stop){
  int i, 
    LBH = L * B * H, 
    s=0;
  for (i=0; i<LBH; s+=Z[i++] == S);
  //printf(" s=%d\n",  s);
  if (stop) assert(s == 10000);
}

void checkS(int *Z, int L, int B, int H){
  checkS(Z, L, B, H, true);
}

void sort(int *a, int n) {
  int i,j;
  // for (i=0; i<n; i++) printf("i=%d %d %d\n", i, n, a[i]);
  for (i=0; i<n; i++) {
    for (j=i; j<n; j++) {
      if (a[j] < a[i]) {
	int dummy = a[i];
	a[i] = a[j];
	a[j] = dummy;
      }
    }
  }
  //  for (i=0; i<n; i++) printf("sort: i=%d %d\n", i, a[i]);
}

void update_all(int *idx_array, int n, int maxdepth, double* pi,
		double* cum_p) { // 0.65
  int d, newbase, zeiger, idx, i,
    base = (1 << maxdepth) - 1;
   sort(idx_array, n);
   for (zeiger=-1, i=0; i<n; i++) {
    idx = idx_array[i];
    if (zeiger>=0 && idx_array[zeiger] == idx) continue;
    //   printf("%d %d %d %d \n", i, zeiger+1, idx_array[zeiger], idx);
    idx_array[++zeiger] = idx;   
    double *ppi = pi + 3 * idx;
    cum_p[base + idx] = ppi[0] + ppi[1] + ppi[2];
  }
  for (d=maxdepth-1; d>=0; d--) {
    newbase = (base + 1) / 2 - 1;
    n = zeiger;
    for (zeiger=-1, i=0; i<=n; i++) {
      idx = idx_array[i] / 2;
      if (zeiger>=0 && idx_array[zeiger] == idx) continue;
      //   printf("%d %d %d %d\n", d, i, zeiger+1, idx);
      idx_array[++zeiger] = idx;    
      int base_idx = base + 2*idx;
      cum_p[newbase + idx] = cum_p[base_idx] + cum_p[base_idx + 1];
    }
    base = newbase;
  }
  // assert(false);
}


// WK'en fuer die einzelnen Doppelpixel
void set_dblt_pi(int *Z, int i, int j, int k, int p, int L, int B, int VARIABLE_IS_NOT_USED H, 
		 double* lambda, int *states, double* pi, char* doublet,
		 int VARIABLE_IS_NOT_USED maxdepth, double VARIABLE_IS_NOT_USED * cum_p, bool VARIABLE_IS_NOT_USED update, int *Idx) {
  assert(k < H - 1 || p < 2);
  int  idx2, 
    im = i % L,
    jm = j % B,
    jL = jm * L, 
    LB = L * B,
    kLB = k * LB,
    idx = kLB + jL + im,
    p_idx = p + 3 * idx;

  idx2 = (p==0) ? kLB + jL + ((i + 1) % L)
    : (p==1) ?  kLB + ((j + 1) % B) * L + im
    : idx + LB;

  //  if (update) 
  //   printf("%d %d Z=%d %d %d states=%d\n", idx, idx2, Z[idx], Z[idx2], p,
  //	   states[ Z[idx] * 12 +  Z[ idx2 ] * 3 + p]);

  doublet[p_idx] = states[ Z[idx] * 12 +  Z[ idx2 ] * 3 + p];
  pi[p_idx] = lambda[(int) doublet[p_idx]];
  *Idx = idx; 
}



// WK'en fuer die einzelnen Doppelpixel
void set_dblt_pi_all(int *Z, int i, int j, int k, int L, int B, int VARIABLE_IS_NOT_USED H, 
		     double* lambda, int *states, double* pi, char* doublet,
		     int VARIABLE_IS_NOT_USED maxdepth, double VARIABLE_IS_NOT_USED *cum_p, bool VARIABLE_IS_NOT_USED update, int **Idx) {
  int  vidx[3], 
    im = i % L,
    jm = j % B,
    jL = jm * L, 
    LB = L * B,
    kLB = k * LB,
    idx = kLB + jL + im,
    p_idx = 3 * idx;

  vidx[0] = kLB + jL + ((i + 1) % L);
  vidx[1] =  kLB + ((j + 1) % B) * L + im;
  vidx[2] = idx + LB;

  //  if (update) 
  //   printf("%d %d Z=%d %d %d states=%d\n", idx, idx2, Z[idx], Z[idx2], p,
  //	   states[ Z[idx] * 12 +  Z[ idx2 ] * 3 + p]);

  int p,
    endfor = k<H-1 ? 3 : 2;
  for (p=0; p<endfor; p++, p_idx++) {
    doublet[p_idx] = states[ Z[idx] * 12 +  Z[ vidx[p] ] * 3 + p];
    pi[p_idx] = lambda[(int) doublet[p_idx]];
    **Idx = idx; 
    (*Idx)++;
  }
}


void random_doublet(int *i, int *j, int *k, int *p, int *idx,
		    int L, int B, int VARIABLE_IS_NOT_USED H, int maxdepth,
		    double VARIABLE_IS_NOT_USED *cum_p, double* pi) {
  int d, g,
    LB = L * B;
    //size = LB * H;
  double *ppi,
    prob = unif_rand() * cum_p[0];  

  for (g=0, d=1; d<=maxdepth; d++) {
    g *= 2;
    double cum = cum_p[ (1 << d) - 1 + g];
    if (prob > cum) {
      prob -= cum;
      g++;
      //assert(prob <= cum_p[d][g]);
    }
  }
  *idx = g;
  *p = 0;
  ppi = pi + 3 * g;
    
  // printf("prob=%f %d %d %f\n", prob, *p, g, cum_p[maxdepth][g]);
  while (prob > ppi[*p]) {
    // printf("prob=%f %f %d %d %f\n", prob, ppi[*p], *p, g, cum_p[d-1][g]);
    prob -= ppi[(*p)++];
  }
  assert(*p <= 2);
  *k = g / LB;
  assert(*k < H);
  g -= *k * LB;
  *j = g / L;
  assert(*j < B);
  *i = g - *j * L; 
  assert(*i < L);
}




void duene(int *Z, int L, int B, int H, int maxdepth, 
	   int *states, int *change, double *lambda, int l_lambda,
	   char* doublet, double* pi, double* cum_p,
	   int maxtime, bool periodic) {
  int i, j, k, p, idx, m, status, 
    idx_array[6 * 3],
    *Idx,
    //newstatus = 0,
    LB = L * B,
    //LBH = LB * H,
    //H1 = H - 1,
    B1 = B - 1,
    L1 = L - 1;
  double time = 0;
  bool ShareAtZero = false;

  while (time < maxtime) {
    time++;

    //  printf("time =%d %d\n", time, maxtime);

    Idx = idx_array;
    random_doublet(&i, &j, &k, &p, &idx, L, B, H, maxdepth, cum_p, pi);//0.44
    status = (int) doublet[p + 3 * idx];
    if (false) {
      if (status <= 0 || change[status] < 0 || change[status + l_lambda] < 0 ||
	  change[status] >20 || change[status + l_lambda] > 20)
	PRINTF("status=%d %d %d %d %d -- lp1=%d change %d %d \n", 
	       status, p,k,j,i,
	       l_lambda,
	       change[status], change[status + l_lambda]);
      //    printf("i=%d %d %d p=%d Z[idx]=%d, %d -> %d, %d \n", 
      //   i,j,k,p, Z[idx], status, change[status], change[status + l_lambda]);
    }

    Z[idx] = change[status];   
    for (m=0; m<3; m++) {
      //printf("%d %d %d k=%d %d\n", m, i, j, k, p);
      
      //  block kostet 0.935, davon 0.65 fuer update_cum_p
      set_dblt_pi(Z, (i+L1) % L, j, k, 0, ARGS, true, Idx); Idx++;
      set_dblt_pi(Z, i, (j+B1) % B, k, 1, ARGS, true, Idx); Idx++; 
      if (k>0) {
	set_dblt_pi(Z, i, j, k-1, 2, ARGS, true, Idx); 
	Idx++; 
      }

      /*
      set_dblt_pi(Z, i, j, k, 0, ARGS, true, Idx);
      set_dblt_pi(Z, i, j, k, 1, ARGS, true, Idx); 
      if (k<H1) {
	set_dblt_pi(Z, i, j, k, 2, ARGS, true, Idx); 
	Idx++;
      }
      */
      set_dblt_pi_all(Z, i, j, k, ARGS, true, &Idx);
      
     if (p==0) {
	if (m==0) {
	  i = (i + 1) % L;
	  ShareAtZero = i == 0 && !periodic && change[status + l_lambda] == S;
	  Z[k * LB + j * L + i] = ShareAtZero ? F : change[status + l_lambda];
	} else {
	  if (m != 1 || !ShareAtZero) break;
	  int idx2;
	  while (true) {	    
	    k = 1 + (int) (unif_rand() * (H-1));
	    j = (int) (unif_rand() * B);
	    idx2 = k * LB + j * L + 0;
	    assert(idx2 < LB * H);
	    if (Z[idx2] == F) break;
	  }
	  Z[idx2] = S;
	}
      } else {
	if (m==1) break;	
	if (p==1) {
	  j = (j + 1) % B;
	  Z[k * LB + j * L + i] = change[status + l_lambda];
	} else { // p = 2
	  k++;
	  Z[idx + LB] = change[status + l_lambda]; // passt das ? idx vs idx2 ?
	}
      }
    } // for m

    update_all(idx_array, Idx - idx_array, maxdepth, pi, cum_p);

    // checkS(Z, L, B, H); // assert(false);


  } // while time
}


void duenen(int *Z, int *LL, int *BB, int *HH,
	    int *states, int *change, double *lambda, int *l_lambda,
	    int *maxtime, int *continuing, int *Lperiodic) { 

  static double* cum_p = NULL;
  static double* pi = NULL;
  static char* doublet = NULL;
  static int maxdepth = -1;

  int i, j, k, p, size, d, Idx,
    L = *LL,
    B = *BB,
    H = *HH,
    H1 = H - 1,
    LBH = L * B * H,
    LBH3 = 3 * LBH;
  bool cont = (bool) *continuing;
  
  if (!cont) {
    if (pi != NULL) { free(pi); pi=NULL; }
    if (doublet != NULL) { free(doublet); doublet=NULL; }
    if (cum_p != NULL) {free(cum_p); cum_p=NULL;}
    maxdepth = (int) (log((double)LBH) / log(2.) + 0.5);
    assert(pow(2.0, maxdepth == LBH));
    pi = (double*) malloc(sizeof(double) * LBH * 3);
    doublet = (char*) malloc(sizeof(char) * LBH * 3);
    cum_p = (double *) malloc(sizeof(double) * 2 * LBH);
      
    for (i=0; i<LBH3; i++) { 
      pi[i] = 0.0;
      doublet[i] = 0;
    }
			      
    for (i=0; i<L; i++) {
      for (j=0; j<B; j++) {
	for (k=0; k<H1; k++) {
	  for (p=0; p<=2; p++) {
	    set_dblt_pi(Z, i, j, k, p, ARGS, false, &Idx);
	  }
	}
      }
    }
    
    int newbase, 
      base = LBH - 1;
    for (i=0; i<LBH; i++) {
      double *ppi = pi + 3 * i;
      cum_p[base + i] = ppi[0] + ppi[1] + ppi[2];
    }

    size = LBH;
    for (d=maxdepth-1; d>=0; d--) {
      size /= 2;
      newbase = (base + 1) / 2 - 1;
      for (i=0; i<size; i++) 
	cum_p[newbase + i] = cum_p[base + 2*i] + cum_p[base + 2*i+1];
      base = newbase;
    }
  }
  
 
  GetRNGstate();  
  // int totalS = 0; for (i=0; i<LBH; i++) totalS += (Z[i] == S); printf("totalS = %d\n", totalS);
  duene(Z, L, B, H, maxdepth,
	states, change, lambda, *l_lambda,
	doublet, pi, cum_p, *maxtime, (bool) *Lperiodic);
  // printf("here\n");

  PutRNGstate();
 }
