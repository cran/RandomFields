/*
 Authors  Martin Schlather, Martin.Schlather@uni-bayreuth.de 

 Simulation of a random field by circular embedding
 (see Wood and Chan, or Dietrich and Newsam for the theory )

 Copyright (C) 2000 Martin Schlather, 

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


#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include "RFsimu.h"
#include <assert.h>

// z coordinate run the fastest in values, x the slowest   

//#define debug_tools 1
#define TOOLS_MEMORYERROR 1
#define TOOLS_XERROR 2
#define TOOLS_BIN_ERROR 3

double EV_TIMES_LESS = 2.0; // 0:always new method; inf:always usual method
int EV_SYSTEM_LARGER = 100; // 0:always new method; inf:always usual method

Real variogram(Real a, Real b) {register Real dummy;dummy=a-b;return dummy*dummy;}
Real Efunction(Real a, Real b) {return a + b;} 
/*  see Schlather (2001), Bernoulli, 7(1), 99-117
    Schlather (????), Math. Nachr, conditionally accepted
    Schlather, Ribeiro, Diggle (in prep.) 
*/
Real kmm(Real a, Real b) {return a * b;}  /* see Stoyan, Kendall, & Mecke, 
					     or Stoyan & Stoyan (both Wiley) 
					  */

int LOW = -1;

void empiricalvariogram(Real *x, Real *y, Real *z, int *dim, int *lx, 
			Real *values, int *repet, int *grid,
			Real *bin, int *nbin, int *charact, 
			Real *res, int *n)
/*    x,y,z  : vector of coordinates, 
 *    dim    : dimension,
 *    lx     : length of x,y, and z
 *    values : (univariate) data for the points
 *    repet  : number of repetitions (the calculated emprical variogram will 
 *             be an overall average
 *    grid   : (boolean) if true lx must be 3 [ x==c(start,end,step) ]
 *    bin    : specification of the bins 
 *             sequence is checked whether it is strictly isotone
 *    nbin   : number of bins. i.e. length(bin)-1
 *    character: 1 : function E
 *               2 : Stoyan's kmm function
 *               default : semivariogram
 *    res    : empirical variogram
 *    n      : number of pairs for each bein
 */
{ 
  int i,ii,j,d,halfnbin,gridpoints[MAXDIM],dimM1,EV_METHOD,error;   
  long totallength,totallengthrepet, segment, factor; 
  Real (*characteristic)(Real, Real);
  Real  *xx[MAXDIM],maxdist[MAXDIM],dd,*BinSq; 

#ifdef debug_tools
  bool DEBUG=true;
  Real Speicher[1000];
  int SPI[100][100][100],SPIX[100][100][100],zaehler,sumzaehler=0;
#endif  

  BinSq = NULL;
  xx[0]=x; xx[1]=y; xx[2]=z;

  if (xx[0]==NULL) {error=TOOLS_XERROR; goto ErrorHandling;}
  for (i=0; i<*nbin;i++) {
    if (bin[i]>=bin[i+1])  {error=TOOLS_BIN_ERROR; goto ErrorHandling;}
  }


  dimM1 = *dim-1;
  halfnbin = *nbin / 2;

  switch (*charact) {
  case 0 :
    characteristic = variogram;
    factor = 2 * *repet;
    break;
  case 1 : 
    characteristic = Efunction;
    factor = 2 * *repet;
    break;
  case 2  : 
    characteristic = kmm;
    factor = *repet;
    break;
  default : assert(false);
  }
   
  //  if ((n = (int*) malloc(sizeof(int)* (*nbin + 1)))==NULL) { 
  //    error=TOOLS_MEMORYERROR; goto ErrorHandling;
  //  }
  if ((BinSq = (Real *) malloc(sizeof(Real)* (*nbin + 1)))==NULL) {
    error=TOOLS_MEMORYERROR; goto ErrorHandling; 
  }
  for (i=0;i<*nbin;i++){res[i]=0.0;n[i]=0;}
  for (i=0;i<=*nbin;i++){if (bin[i]>0) BinSq[i]=bin[i] * bin[i]; 
  else BinSq[i]=bin[i];
  // PRINTF(" %d %f %f \n",i,bin[i],BinSq[i]);
  }


  //////////////////////////////////// GRID ////////////////////////////////////
  if (*grid) {
    int d1,d2,head,tail,swapstart;
    Real p1[MAXDIM],p2[MAXDIM],distSq,dx[MAXDIM]; 
    int delta[MAXDIM],deltaTwo[MAXDIM],maxtail,low,cur;    
    long indextail[MAXDIM],indexhead[MAXDIM],valuePtail,valuePointhead,
      segmentbase[MAXDIM],DeltaTwoSegmentbase[MAXDIM],
      DeltaFourSegmentbase[MAXDIM],
      SegmentbaseTwo[MAXDIM],SegmentbaseFour[MAXDIM];
    Real psq[MAXDIM],p[MAXDIM],maxbinsquare;
    bool forbid;    

    // does not work!! :
    // GetGridSize(x,y,z,dim,&(gridpoints[0]),&(gridpoints[1]),&(gridpoints[2])); 
    // totallength = gridpoints[0] * gridpoints[1] * gridpoints[2];
    //
    // instead :
    for (dd=0.0,totallength=1,i=0; i<=dimM1; i++) {
      gridpoints[i] = (int) ((xx[i][XEND]-xx[i][XSTART])/xx[i][XSTEP]+1.5);
      totallength *= gridpoints[i];
      maxdist[i]=(gridpoints[i]-1)*xx[i][2]; dd += maxdist[i]*maxdist[i];
    }
#ifdef debug_tools
    if (DEBUG) { EV_METHOD=0; } else 
#endif  
      EV_METHOD = (int) ( ((EV_TIMES_LESS * BinSq[*nbin]) < dd) ||
			  (EV_SYSTEM_LARGER<totallength)); // ? large system ?    
    for (i=0;i<=dimM1;i++) { dx[i] = xx[i][2] * 0.99999999; }
    totallengthrepet = totallength * *repet;
    segmentbase[0]=1; // here x runs the fastest;  segmentbase[dimM1]=1 else!
    SegmentbaseTwo[0]=segmentbase[0] << 1;//SegmentbaseTwo[i-1]=segmentbase[i-1] <<1;
    SegmentbaseFour[0]=segmentbase[0] << 2;//SegmentbaseFour[i-1]=segmentbase[i-1] <<2;
    for (i=1;i<=dimM1;i++) { //for (i=dimM1;i>0;i--)
      segmentbase[i]=segmentbase[i-1] * gridpoints[i-1];//segmentbase[i-1]=segmentbase[i] * gridpoints[i];
      SegmentbaseTwo[i]=segmentbase[i] <<1;//SegmentbaseTwo[i-1]=segmentbase[i-1] <<1;
      SegmentbaseFour[i]=segmentbase[i] <<2;//SegmentbaseFour[i-1]=segmentbase[i-1] <<2;
    } 
    
    switch (EV_METHOD) {
    case 0 :
      for (i=0;i<=dimM1;i++) {indexhead[i]=0; p1[i]=0.0;}

#ifdef debug_tools
      for (i=dimM1; i<MAXDIM; i++ ) {indexhead[i]=indextail[i]=0;} 
#endif  
       
      // loop through all pair of points, except (head,tail) for head=tail, 
      // which is treated separately at the end, if necessary (not relevant for 
      // variogram itself, but for function E and V)
      for (head=0;head<totallength;) {
	for (i=0;i<=dimM1;i++) {indextail[i]=0; p2[i]=0.0;}
	for (tail=0;tail<head;) {
	  distSq=0.0;
	  for (i=0;i<=dimM1;i++) {
	    Real dx;
	    dx = p1[i]-p2[i];
	    distSq += dx * dx;
	  }
	  //	  distSq = sqrt(distSq); 26.2.
  
#ifdef debug_tools
	  if ((abs(indexhead[0]-indextail[0])<99) && 
	      (abs(indexhead[1]-indextail[1])<99) &&
	      (abs(indexhead[2]-indextail[2])<99))
	    (SPI[abs(indexhead[0]-indextail[0])][abs(indexhead[1]-indextail[1])]
	     [abs(indexhead[2]-indextail[2])])++;
#endif  
	  
	  if ((distSq>BinSq[0]) && (distSq<=BinSq[*nbin])) { 
	    /* search which bin distSq in */
	    register int up, cur, low;
	    low=0; up= *nbin; /* 21.2.01, *nbin-1 */  cur= halfnbin;
	    while (low!=up) { 
	      if (distSq> BinSq[cur]) {low=cur;} else {up=cur-1;} /* ( * ; * ] */ 
	      cur=(up+low+1)/2; 
	    }
#ifdef debug_tools
	    if (low==LOW) {
	      segment=0; 
	      PRINTF(" value %d %f %d %f %f low=%d\n",
		     head+segment,values[head+segment],
		     tail+segment,values[tail+segment], res[low],low ); }
#endif 
	    for (segment=0;segment<totallengthrepet;segment += totallength) {	
	      res[low] += 
		characteristic(values[head+segment],values[tail+segment]);
	    } 	
	    //for (i=0;i<*nbin;i++){PRINTF("%f ",res[i]);}
	    assert(low<*nbin);
	    n[low]++;

#ifdef debug_tools
	  if ((abs(indexhead[0]-indextail[0])<99) && 
	      (abs(indexhead[1]-indextail[1])<99) &&
	      (abs(indexhead[2]-indextail[2])<99))
	    (SPIX[abs(indexhead[0]-indextail[0])][abs(indexhead[1]-indextail[1])]
	     [abs(indexhead[2]-indextail[2])])=low;
#endif  

	  }	
	  tail++;
	  d2=dimM1; //d2=0;
	  indextail[d2]++; p2[d2]+=dx[d2];
	  while (indextail[d2]>=gridpoints[d2]) { 
	    indextail[d2]=0; p2[d2]=0; 
	    d2--; //d2++;
	    assert(d2>=0); // assert(d2<dim);	
	    indextail[d2]++; p2[d2]+=dx[d2];
	  }
	}   
	head++;
	if (head<totallength) {
	  d1=dimM1;  // d1=0;  // if x should run the slowest
	  indexhead[d1]++; p1[d1]+=dx[d1];      
	  while (indexhead[d1]>=gridpoints[d1]) { 
	    indexhead[d1]=0; p1[d1]=0; 
	    d1--;   // d1++;
	    assert(d1>=0); 
	    indexhead[d1]++; p1[d1]+=dx[d1];
	  }
	}
      }
      if(GENERAL_PRINTLEVEL>=8) {
	for(i=0;i<*nbin;i++){
	  PRINTF(" %d:%f(%f,%d)",i,res[i]/(Real)(factor * n[i]),res[i],n[i]);
	}
	PRINTF("  XXX\n");}      
#ifdef debug_tools
      if (DEBUG) 
	for (i=0;i<*nbin;i++) Speicher[i]=res[i]/(Real)(factor * n[i]);
      else 
#endif  
	break;

      /************************************************************************/
      /************************************************************************/

    case 1 :
      /* idea: 
	 1) fix a vector `delta' of positive components
	 2) search for any pair of points whose distSqance vector dv equals
	    |dv[i]|=delta
	    and some the f(value1,value2) up.
	 
	realisation :
	1) one preferred direction, namely dimM1-direction (of 0..dimM1)
	2) for the minor directions, swap the values of delta from positive 
	   to negative (all 2^{dimM1} posibilities)   --> indexhead
	   EXCEPTION :  delta[dimM1]==0. Then only the minor directions up to 
                        dimM1-2 are swapped!
                        (otherwise values will be counted twice as indextail is 
                         running through all the positions with
                         coordinate[dimM1]==0)
           if delta[dimM1-1] is also 0 then only minor direction up to dimM1-3 
           are swapped, etc.
        3) indextail running through all positions with coordinate[dimM1]==0. 
	   This is not very effective if the grid is small in comparison to the
	   right boundary of the bins!    
      */
    
      for (i=0;i<*nbin;i++){res[i]=0.0;n[i]=0;}
      
#ifdef debug_tools
      for (i=dimM1;i<MAXDIM;i++) {delta[i]=0;} 
#endif  

      maxbinsquare = BinSq[*nbin];
      for (i=0;i<=dimM1;i++) {delta[i]=0; psq[i]=p[i]=0.0;}
      for (;;) { // delta
	do {
	  d=0;  // d=dimM;  // if x should run the lowest
	  delta[d]++;  // this excludes distSq==0 !	
#ifdef debug_tools
	  zaehler = 0;
	  if (GENERAL_PRINTLEVEL>=10) 
	    PRINTF("\n {%d %d %d}",delta[0],delta[1],delta[2]);
#endif  
	  p[d]+=dx[d]; psq[d]=p[d]*p[d]; //////////
	  for (distSq=0.0,i=0;i<=dimM1;i++) {distSq+=psq[i];}
	  //PRINTF(" dsq %f ",distSq);
	  while ((distSq>maxbinsquare)  || (p[d]>maxdist[d])) {/* (*) */ 
#ifdef debug_tools
	    if (GENERAL_PRINTLEVEL>=10) 
	      PRINTF("failure: %f %f   %f %f %d\n",
		     distSq,maxbinsquare,p[d],maxdist[d],d);
#endif  
	    delta[d]=0; psq[d]=p[d]=0.0; 
	    d++;   // d--;
	    if (d<=dimM1) { //d>=0
	      delta[d]++; p[d]+=dx[d]; psq[d]=p[d]*p[d];
	      for (distSq=0.0,i=0;i<=dimM1;i++) {distSq+=psq[i];}
	      //  PRINTF(" xdsq %f ",distSq);
	    } else break;	
	  }
	} while ((distSq<=BinSq[0]) && (d<=dimM1));
	if (d>dimM1) break;//(d<0)
	// distSq = sqrt(sum);
	// PRINTF("distSq %f %f %f\n",distSq,BinSq[*nbin],BinSq[0]);
	assert((distSq<=BinSq[*nbin]) && (distSq>BinSq[0]));
	{
	  register int up;
	  low=0; up= *nbin; /* */ cur= halfnbin;
	  while(low!=up){
	    if (distSq> BinSq[cur]) {low=cur;} else {up=cur-1;} /* ( * ; * ] */ 
	    cur=(up+low+1)/2; }
	}	

#ifdef debug_tools
	if (GENERAL_PRINTLEVEL>=10){
	  PRINTF(" [%d %d] ",SPI[delta[0]][delta[1]][delta[2]],
		 SPIX[delta[0]][delta[1]][delta[2]]);
	  PRINTF("distSq: %f bin:%f %f %d  ",distSq,BinSq[low],BinSq[low+1],low);
	}
#endif  
	
	i=0; //i=dimM1;
	while ((i<dimM1) && (delta[i]==0)) i++; //while ((i>0) && (delta[i]==0)) {i--;} 
	swapstart=i+1; //swapstart=i-1;
	
	for(i=0;i<=dimM1;i++) {
	  deltaTwo[i]=delta[i]<<1;
	  DeltaTwoSegmentbase[i] = delta[i] * SegmentbaseTwo[i];
	  DeltaFourSegmentbase[i]= delta[i] * SegmentbaseFour[i];
	}

	for (i=0;i<=dimM1;i++) {indextail[i]=0;} 
	valuePtail = valuePointhead =0;
	
	maxtail=gridpoints[0]-delta[0]; //maxtail=gridpoints[dimM1]-delta[dimM1];
     
	for(;;) { // indextail
	  for (i=0,valuePointhead=0; i<=dimM1;i++) {
	    indexhead[i]=indextail[i]+delta[i]; /* Do not use -delta[i], as last 
						   dimension must always be 
						   +delta[] */
	    valuePointhead += indexhead[i] * segmentbase[i];
	  }  // one point is given by indextail, the other by indexhead	
	  for (;;) {// indexhead
	    forbid = false; // in case dim=1 ! (???, 21.2.01)
	    for (i=0; i<=dimM1; i++) { // for (i=1; ... ) should be enough...
	      if (forbid=(indexhead[i]<0)||(indexhead[i]>=gridpoints[i])){break;}
	    }
	    if (!forbid) { 	  
	      long segtail,seghead;
	      for (tail=0; tail<maxtail; tail++) {	
#ifdef debug_tools
		zaehler++;
#endif  	
		
		segtail = valuePtail+tail;
		seghead = valuePointhead+tail;
#ifdef debug_tools
		if (low==LOW) {
		  segment=0; 
		  PRINTF("\ncase1 value %d %d  %f %f %f {t%d,%d} {h%d,%d} %d",
			 segtail,seghead,values[segtail+segment],
			 values[seghead+segment],res[low],
			 indextail[0],indextail[1],indexhead[0],indexhead[1],tail); }
#endif  
		for (segment=0;segment<totallengthrepet;segment+=totallength) {	  
		  res[low] += characteristic(values[segtail+segment],
					     values[seghead+segment]);
		} 	
		assert(low<*nbin);
		n[low]++;
	      }	
	    }
	    
	    d=swapstart;  
	    if (d<=dimM1) {//(d>=0)
	      indexhead[d]-=deltaTwo[d];  
	      valuePointhead-= DeltaTwoSegmentbase[d];
	      while ( (indexhead[d]<indextail[d]-delta[d]) || (delta[d]==0)) { 
		if (delta[d]!=0) {
		  // so, indexhead[d]<indextail[d]-delta[d]; 
		  // i.e.indexhead[d]=indextail[d]-3*delta[d]
		  valuePointhead += DeltaFourSegmentbase[d]; 
		  indexhead[d]=indextail[d]+delta[d]; 
		}
		d++; //d--; 
		if (d<=dimM1) {//(d>=0)
		  indexhead[d]-=deltaTwo[d]; 
		  valuePointhead-=DeltaTwoSegmentbase[d];
		} else break;
	      }
	    } else break; 
	    if (d>dimM1) break;//(d<0)
	  } //for(;;), indexhead
	  
	  d=1;  // d=dimM1-1;  // if x should run the fastest
	  if (d<=dimM1) {//(d>=0)    // not satisfied for one dimensional data!
	    indextail[d]++;  valuePtail+=segmentbase[d];
	    while (indextail[d]>=gridpoints[d]) { 
	      valuePtail -= indextail[d] * segmentbase[d];
	      indextail[d]=0;
	      d++; //d--; 
	      if (d<=dimM1) {//(d>=0) 
		indextail[d]++; valuePtail+=segmentbase[d];
	      } else break;
	    }
	  } else break;
	  if (d>dimM1) break; // (d<0)
	} // for(;;), indextail
#ifdef debug_tools
	sumzaehler = sumzaehler + zaehler;
        if(GENERAL_PRINTLEVEL>=10) PRINTF(" z=%d %d",zaehler,sumzaehler);
#endif
      }  //for(;;), delta
      if(GENERAL_PRINTLEVEL>=8) {
	PRINTF("\n");
	for(i=0;i<*nbin;i++){
	  PRINTF(" %d:%f(%f,%d)  ",i,res[i]/(Real)(factor * n[i]),res[i],n[i]);
	}
	PRINTF(" YYY\n");
      }
#ifdef debug_tools
      if (DEBUG) {
	for (i=0;i<*nbin;i++) {
	  if (fabs(Speicher[i] - ( res[i]/(Real)(factor * n[i])))/Speicher[i] >
	      0.000000001 ) {
	    PRINTF("Inconsistency: %d %.15f %.15f %.15f\n",
		   i,Speicher[i],res[i]/(Real)(factor * n[i]),
		   Speicher[i] - res[i]/(Real)(factor * n[i]));
	    LOW = i;
	    if (GENERAL_PRINTLEVEL<8) {
	      GENERAL_PRINTLEVEL=10;
	      empiricalvariogram(x,y,z,dim,lx,values,repet,grid,bin,nbin,charact,
				 res, n);
	    } else {
	      PRINTF("\nBin ");
	      for (i=0;i<=*nbin;i++) PRINTF("%d:%f ",i,BinSq[i]);
	      PRINTF("\nSPEZ\n");
	      for (i=0;i<=dimM1;i++) 
		PRINTF(" %f %f %f\n",xx[i][XSTART],xx[i][XEND],xx[i][XSTEP]);
	      PRINTF("\n totallength=%d, %d %d %d \n",
		     totallength,gridpoints[0],gridpoints[1],gridpoints[2]);
	      assert(false);
	    }
	  }	  
	}
      }
#endif
      break;
    default : assert(false);
    }
  } else {
    ////////////////////////////////////  ARBITRARY /////////////////////////////
    totallength = *lx;
    totallengthrepet = totallength * *repet;
    for (i=0;i<totallength;i++){ // to have a better performance for large 
      //                            data sets, group the data first into blocks
      for (j=0;j<i;j++){
        Real distSq;
      	for (distSq=0.0,d=0; d<=dimM1;d++) {
	  register Real dx;
	  dx=xx[d][i]-xx[d][j];
	  distSq += dx * dx;
	}		
	// see also above
	//distSq = sqrt(distSq); 26.2.
	if ((distSq>BinSq[0]) && (distSq<=BinSq[*nbin])) { 
	  register int up, cur,low;
	  low=0; up= *nbin; cur= halfnbin;
	  while (low!=up) {  
	    if (distSq> BinSq[cur]) {low=cur;} else {up=cur-1;} // ( * ; * ]  
	    cur=(up+low+1)/2; 
	  } 	
	  for (segment=0;segment<totallengthrepet;segment += totallength) {	  
	    res[low] += characteristic(values[i+segment],values[j+segment]);
	  } 	
	  assert(low<*nbin);
	  n[low]++;
	}
	/* */
      }
    }
  }
  
  // ii is used in switch by this value: !!!!
  ii=0;while( (ii<=*nbin) && (bin[ii]<0) ){ii++;}
  ii--;

  for (i=0;i<*nbin;i++){
    if (n[i]>0) { res[i]/= (Real) (factor * n[i]); } 
    else if (i!=ii) {res[i]= RF_NAN;}
  }

 
  switch (*charact) {
  case 0 :
    if ((ii<*nbin)&&(ii>=0)){
      res[ii] *= ((Real) (factor * n[ii]));
      n[ii] += totallength;
      res[ii] /=  (Real) (factor * n[ii]);
    }    
    break;
  case 1 : // E, calculating E(0) correctly, taking into account that the
    //        bin including 0 may include further points
    if ((ii<*nbin)&&(ii>=0)){
      long j; Real sum;
      res[ii] *= ((Real) (factor * n[ii]));
      for (sum=0.0, j=0;j<totallengthrepet;j++) { sum += values[j]; }
      n[ii] += totallength;
      res[ii] = (res[ii] + 2.0 * sum) / (Real) (factor * n[ii]);
      printf(" sum=%f %f %d\n", sum,res[ii],ii);
    }
    break;
  case 2 : // kmm, calculating kmm(0} and dividing all the values by mean^2
    register Real mean,square;
    for (j=0,mean=square=0.0;j<totallengthrepet;j++) { 
      mean += values[j]; square += values[j]*values[j];
    }
    mean /= (Real) totallengthrepet;
    square /= (Real) totallengthrepet;

    if ((ii<*nbin)&&(ii>=0)){
      res[ii] *= ((Real) (factor * n[ii]));
      n[ii] += totallength;     
      res[ii]= (res[ii] + square) / (Real) (factor * n[ii]);
    }
    mean *= mean;  for (j=0;j<*nbin;j++) { res[j] /= mean; }
    break;
  }
  //  if (n!=NULL) free(n); 
  if (BinSq!=NULL) free(BinSq);
  return;
  
 ErrorHandling:
  PRINTF("Error: ");
  switch (error) {
  case TOOLS_MEMORYERROR :  
    PRINTF("Memory alloc failed in empiricalvariogram.\n"); break;
  case TOOLS_XERROR :  
    PRINTF("The x coordinate may not be NULL.\n"); break;
  case TOOLS_BIN_ERROR :
    PRINTF("Bin components not an increasing sequence.\n"); break;
  default : assert(false);
  }
  //if (n!=NULL) free(n); 
  if (BinSq!=NULL) free(BinSq);
#ifdef RF_GSL
  assert(false); 
#else
  for (i=0;i<*nbin;i++){res[i]=RF_NAN;} 
#endif
} 
	 
       
// A simple function for calculating the empirical variogram
// only equidistant bins, only non-grids allowed
// obsolete soon?
void binnedvariogram( Real *x, Real *y, Real *z, int *lx, Real *step,
		 Real *res, int *n, int *lr)  
{ 
  int i,j,index;
  Real stepinv;
 //  int *n;
 
 //  n = (int*) malloc(sizeof(int)* *lr);
  stepinv = 1/ *step;
  for (i=0;i<*lr;i++){res[i]=0;n[i]=0;}
  for (i=0;i<*lx;i++){ /* if (i % 128 ==0) PRINTF("%d\n",i); */
    //    PRINTF("III %d\n",i);
    for (j=0;j<i;j++){
      register Real d,dx;
      dx=x[i]-x[j];
      d=y[i]-y[j];
      d=dx*dx + d*d;
      if (d>0) {
	index = (int) floor(sqrt(d)*stepinv)+1;
	//if (index>*lr) {PRINTF(" %d > %d\n",index,*lr);exit(0);}
	if (index<*lr) {
	  //	PRINTF("index %d\n",index);
	  d=z[i]-z[j];
	  res[index] += d*d;
	  n[index]++;
	}
      }
      /* PRINTF(" d=%f, index=%d, v=%f, res=%f, n=%d\n", sqrt(d),
	 index,v,res[index],n[index]);
      */ 
    }
  }
  for (i=1;i<*lr;i++){
    // PRINTF(" %d ",n[i]);
    if (n[i]>0) {res[i]/=( (Real) (2*n[i]));} else {res[i]=-1.0;}
  }
  //  free(n);
  //PRINTF("vario OK\n");
}
	 
         





