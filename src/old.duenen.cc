/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part and error messages

 Copyright (C) 2011 -- 2013 Martin Schlather
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

// must be the same as in duenen.R
#define L 200
#define B 200
#define H 25
#define LB  (L * B)
#define G 0
#define M 1
#define F 2
#define S 3

void GetLH(int *l, int *b, int *h){
  *h = H;
  *b = B;
  *l = L;
}

static int  
states[4][4][3] = {
  // G           M          F          S (vorwaerts, links, oben)
  { { 0, 0, 0}, {0, 0, 0}, {14,16, 0}, {0, 0, 7} }, // G
  { {12, 0,19}, {0, 0,20}, {11, 0, 0}, {9, 0,10} }, // M
  { {13,15,17}, {0, 0, 0}, { 0, 0, 0}, {3, 5, 1} }, // F
  { { 8, 0,18}, {0, 0, 0}, { 4, 6, 2}, {0, 0, 0} }  // S
},

  change[21][2] = { {-1, -1},
		    {S, F}, {F, S}, {S, F}, {F, S}, {S, F}, {F, S}, // 6
		    {S, M}, {S, M}, // 8
		    {S, M}, {S, M}, // 10
		    {F, G}, {G, G}, // 12
		    {G, F}, {F, G}, {G, F}, {F, G}, // 16
		    {G, F}, {G, S}, {G, G}, {G, M}
  };

static double pi[3][H][B][L], pi3[H][B][L], pi3H[B][L], pi3HB[L], pi3HBL;
static char doublet[3][H][B][L];


void checkdoublets(int  VARIABLE_IS_NOT_USED *Z, double  VARIABLE_IS_NOT_USED *lambda,
		   int  VARIABLE_IS_NOT_USED starti, int endi, int startj, int endj, 
		   int startk, int endk) {
  int i, j, k, p,
    H1 = H - 1;
  double dummy;
    for (i=0; i<endi; i++) {
      for (j=startj; j<endj; j++) {
	for (k=startk; k<endk; k++) {
	  dummy = pi3[k][j][i];
	  for(p=0; p<=2; p++) {
	    assert(pi[p][k][j][i] == lambda[ doublet[p][k][j][i] ]);
	    dummy -= pi[p][k][j][i];
	  }
	  assert(fabs(dummy) < 10^-10);
	} // k
	dummy = pi3H[j][i];
	for (k=0; k<H1; k++) dummy -= pi3[k][j][i];
	assert(fabs(dummy) < 10^-10);
      } // j
      dummy = pi3HB[i];
      for (j=0; j<B; j++) dummy -= pi3H[j][i];      
      assert(fabs(dummy) < 10^-10);
    } // i
    dummy = pi3HBL;
    for (i=0; i<L; i++) dummy -= pi3HB[i];
    assert(fabs(dummy) < 10^-10);
}





// WK'en fuer die einzelnen Doppelpixel
void calculatePi(int *Z, double *lambda,
		 int starti, int endi,
		 int startj, int endj,
		 int startk, int endk, bool pr) {
  // assert(false);

  int i, j, k, p, idx, jm, im, kLB, jL,  idx2,
    //    L1 = L - 1,
    H1 = H - 1;

  // print("cpi %d %d %d %d %d %d\n", 
  //	 starti, endi,startj, endj, startk, endk);
 assert(startk >=0 && endk<=H);


  for (i=starti; i<endi; i++) {
    im = i % L;
    for (j=startj; j<endj; j++) {
      jm = j % B;
      jL = jm * L;
      for (k=startk; k<endk; k++) {
	kLB = k * LB;
	assert(k < H1);
	
	if (pr) PRINTF("c-pi %d %d %d; %d %d; %d %d; %d %d\n", i, j, k, 
		      starti, endi, startj, endj, startk, endk);

	idx = kLB + jL + im;
	
	//	print("%d %d %d idx=%d; km=%d %d %d kLB=%d %d %d\n", k, i, j, idx, k, im, jm, kLB, im, jL);
	// print(" %d %d %d %d\n",
	//		 (int)Z[idx],  (int)Z[idx+1], (int)Z[idx + L], (int)Z[idx + LB]);
	
	/*	  doublet[0][k][jm][im] = (im >= L1) ? 0
		  : states[ (int) Z[idx] ][ (int) Z[idx + 1] ][0];
		  doublet[1][k][jm][im] = (jm >= L1) ? 0
		  : states[ (int) Z[idx] ][ (int) Z[idx + L] ][1];
		  doublet[2][k][jm][im] = (k >= H1) ? 0 
		  : states[ (int) Z[idx] ][ (int) Z[idx + LB] ][2];
	*/
	
	idx2 = kLB + jL + ((i + 1) % L);
	assert(idx2 <= LB * H);

	if (pr) PRINTF("idx=%d %d %d %d\n", idx, Z[idx],   idx2,  Z[idx2]);

	doublet[0][k][jm][im] = states[ (int) Z[idx] ][ (int) Z[ idx2 ] ][0];
	//	print("here 1 %d\n", kLB + ((j + 1) % B) * L + im);

	idx2 = kLB + ((j + 1) % B) * L + im;
	assert(idx2 <= LB * H);
	doublet[1][k][jm][im] = states[ (int) Z[idx] ][ (int) Z[idx2] ][1];
	//	print("here 2 %d %d %d Z=%d\n", k + 1, L, ((k + 1) % L) * LB + jL + im, (int) Z[((k + 1) % H) * LB + jL + im]);
	if (k < H1) {
	  idx2 = (k + 1) * LB + jL + im;
	  assert(idx2 <= LB * H);
	  doublet[2][k][jm][im] = states[ (int) Z[idx] ][ (int) Z[idx2] ][2];
	} else {
	  doublet[2][k][jm][im] = 0 ;
	}
	//print("here 3\n");
	
	//doublet[0][k][jm][im] =doublet[1][k][jm][im] =  doublet[2][k][jm][im] =0;

	pi3[k][jm][im] = 0.0;
	//print("here\n");
	for(p=0; p<=2; p++) {
	  pi[p][k][jm][im] = lambda[ (int) doublet[p][k][jm][im] ];
	  pi3[k][jm][im] += pi[p][k][jm][im];
	}
      } // k
      pi3H[jm][im] = 0.0;
      for (k=0; k<H1; k++) pi3H[jm][im] += pi3[k][jm][im];
    } // j
    pi3HB[im] = 0.0;
    for (j=0; j<B; j++) {
      pi3HB[im] += pi3H[j][im];      
      // print("%d %d %f\n", i, j,  pi3H[j][i]); 
    }
  } // i
  pi3HBL = 0.0;
  for (i=0; i<L; i++) (pi3HBL) += pi3HB[i];

  //  assert(false);
} 
 

void calculatePi(int *Z, double *lambda,
		 int starti, int endi,
		 int startj, int endj,
		 int startk, int endk) {
  calculatePi(Z, lambda, starti, endi, startj, endj, startk, endk, false);
}

void duene(int *Z,
	   double lS, double lW, double lWM, double lWP, double lD, double lG, 
	   double a, double MaxTime) {
  int i,j,k,p, starti, startj, startk, endi, endj, endk, status, idx, 
    //LBH = LB * H,
    H1 = H - 1,
    B1 = B - 1,
    L1 = L - 1;
  double choose, time, lambda[21];

  pi3HBL = 0.0;
  for (i=0; i<L; i++) {
    pi3HB[i] = 0.0;
    for (j=0; j<B; j++) {
      pi3H[j][i] = 0.0;
      for (k=0; k<H; k++) {
	pi3[k][j][i] = 0.0;
	for (p=0; p<=2; p++) {
	  pi[p][k][j][i] = 0.0;
	}
      }
    }
  }
    
  //  print("%d\n", (int) Z[0]); 
  // assert(false);
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
 
  //  print("init\n");
  // assert(false);

  
  calculatePi(Z, lambda, 0, L, 0, L, 0, H - 1);
  //  checkdoublets(Z, lambda, 0, L, 0, B, 0, H-1);
  // assert(false);
  time = rexp(1.0);// / pi3HBL);

  PRINTF("time %f %f %f\n", time, MaxTime, pi3HBL);

  int zaehler = 0;

  while (time < MaxTime) {
    //  
    //
    // choose doublet
    zaehler++;
    
    choose = unif_rand() * pi3HBL;
    //   PRINTF("t=%f %f %f\n", time, pi3HBL, choose);
    for (i=0; i < L && pi3HB[i] < choose; i++) choose -= pi3HB[i];
    if (i<L && pi3HB[i] < choose)
      PRINTF("** 3HB %f %f \n", pi3HB[i], choose);
    for (j=0; j < B && pi3H[j][i] < choose; j++) choose -= pi3H[j][i];
    if (i<L && j<B && pi3H[j][i] < choose) 
      PRINTF("** 3H %f %f \n", pi3H[j][i], choose);
    for (k=0; k < H && pi3[k][j][i] < choose; k++) choose -= pi3[k][j][i];
    if (i<L && j<B && k<H1 && pi3[k][j][i] < choose) 
      PRINTF("** 3 %f %f \n", pi3[k][j][i], choose);
    for (p=0; p < 2 && pi[p][k][j][i] < choose; p++) choose -= pi[p][k][j][i];
    if (i<L && j<B && k<H && p<2 && pi[p][k][j][i] < choose)
      PRINTF("** %f %f \n", pi[p][k][j][i], choose);

    //     PRINTF("***&** %f %f \n", pi[p][k][j][i], choose);

    status = doublet[p][k][j][i];
    if (status <= 0 || change[status][0] < 0 || change[status][1] < 0 ||
	change[status][0] >20 || change[status][1] > 20)
      PRINTF("status=%d %d %d %d %d -- change %d %d \n", status, p,k,j,i,
	     change[status][0], change[status][1]);
    //  
    assert(status != 0);
  
    starti = i + L1; // ~= i - 1
    startj = j + B1;
    startk = k > 0 ? k-1 : 0;
    endi = i + L + 1; // ~= i + 1
    endj = j + B + 1; // 
    endk = k == H1 ? H1 : k + 1; //
    assert(endk<=H1);
    idx = k * LB + j * L + i;
    //    print("change %d %d\n", change[status][0], change[status][0]);
    //  print("%d %d %d %d status=%d\n",  p, k, j, i, status);
    assert(change[status][0] >= 0 && change[status][1] >= 0 && 
	   change[status][0] <= 20 && change[status][1] <=20);
  // 

 
    if (false) {
      PRINTF("chosen i=%d j=%d k=%d p=%d; pi=%f status=%d %d->%d %d->%d idx=%d %d %d\n", //
	    i, j, k, p, pi[p][k][j][i],
	     status,  Z[idx], change[status][0],
	     Z[idx + (int) (0.1 + pow((double)L, p))],change[status][1], idx, L, LB);
    }

    //   assert((Z[idx] <= M) + (Z[idx + (int) (0.1+ pow(L, p))] <= M)
    //	   == 
    //	   (change[status][0]<=M) + ( change[status][1] <=M));

    assert(idx < LB * H);
    int oldZ = Z[idx];
    Z[idx] =  change[status][0];
    int newstatus = change[status][1];
    if (p==0) {  
      assert(i < L);
      assert(k * LB + j * L + (i + 1) % L < LB * H);
      if (newstatus == S && (i + 1) % L == 0) {
 
	if (false) {
	  int i0, t=0, T=L*B*H;
	  for (i0=0; i0<T; i0++) t += Z[i0]==S; 
	  PRINTF("\nstart t=%d zaehler=%d oldZ=%d %d change %d %d status=%d \n",
		 t, zaehler, oldZ,
		 Z[k * LB + j * L + (i + 1) % L], 
		 change[status][0], change[status][1], status);
	}
 
 	int k0, j0, idx0;
	Z[k * LB + j * L + 0] = F;
	while (true) {
	  k0 = 1 + (int) (unif_rand() * (H-1));
	  j0 = (int) (unif_rand() * B);
	  idx0 = k0 * LB + j0 * L + 0;
	  assert(idx0 < LB * H);
	  if (Z[idx0] == F) break;
	}
	Z[idx0] = S;
 
	//	if (zaehler == 13628) { printf("%d %d %d\n", L, j0, k0);  }
	
 	calculatePi(Z, lambda, L-1, L + 1, j0 + B1, j0+B+1, k0-1, 
		    k0 == H1 ? H1 : k0 + 1, zaehler==13628);
	//	checkdoublets(Z, lambda, 0, L, 0, B, 0, H-1);

	if (false) {
	  int ii, t=0, T=L*B*H;
	  for (ii=0; ii<T; ii++) t += Z[ii]==S; 
	  //printf("t=%d zaehler=%d\n", t, zaehler);
	  // assert(t == totalS); 
	}
 
      } else {
	Z[k * LB + j * L + (i + 1) % L] = newstatus;
      }
      endi++;
    } else if (p==1) {
      assert(j < B);
      assert(k * LB + ((j + 1) % B) * L + i <  LB * H);
      Z[k * LB + ((j + 1) % B) * L + i] = newstatus;
      endj++;
    } else  {
      assert(k < H1);
      assert(idx + LB < LB * H);
      Z[idx + LB] = newstatus;
      if (endk < H1) endk++;
    }

    if (false) {
      int ii, t=0, T=L*B*H;
      for (ii=0; ii<T; ii++) t += Z[ii]==S; 
      //printf("t=%d zaehler=%d\n", t, zaehler);
      // assert(t == totalS); 
    }
 
    //  if (zaehler == 13628)
    //  print("here p=%d %d %d; %d %d; %d %d \n", p, starti, endi, startj, endj,
    //	  startk, endk);

    //if (doublet[0][16][85][199] == 4 && Z[16 * LB + 85 * L + 199] == S &&
    //	Z[16 * LB + 85 * L] == S ) {
    // printf("z=%d %f\n", zaehler, pi[0][16][85][199]);
    //  assert(false); // 13628
    // }


    calculatePi(Z, lambda, starti, endi, startj, endj, startk, endk);
    //  checkdoublets(Z, lambda, 0, L, 0, B, 0, H-1);

    // print("hereB\n");

  
    time += rexp(1.0); // / pi3HBL);

    //if (k==20 && j==0 && i==43)
    //   checkdoublets(Z, lambda, doublet, 
    //		0, L, 0, L, 0, H,
    //		pi, pi3, pi3H, pi3HB, &pi3HBL);

  } // while
  
  // Rand
}


void Xduenen(int *ZZ, 
	   double *lambdas, double *lambdaw, double *lambdawm,
	   double *lambdawp, double *lambdad, double *lambdag, 
	    double *a, double *maxtime) { 
  /*
  int i;
  double min=1e20, max=-1e20;
  for (i=0; i<L*L*H; i++) {
    if (ZZ[i] < min) {min = ZZ[i]; print("%f %f \n", min , max); } //
    if (ZZ[i] > max) {max = ZZ[i]; print("%f %f \n", min , max); }//
  }
  assert(min==0 && max==3);
  */
 
  GetRNGstate();  
  // int i, total = L * B * H, totalS = 0; for (i=0; i<total; i++) totalS += (ZZ[i] == S); printf("totalS = %d\n", totalS);
  duene(ZZ, *lambdas, *lambdaw,*lambdawm,
	*lambdawp,*lambdad, *lambdag,  *a,
	 *maxtime);
  // printf("here\n");

  PutRNGstate();
 }
