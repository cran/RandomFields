// Version 27.6.00
//#define DEBUG 1

#define DD if (0)

#include "math.h" 
#include "RFsimu.h"
#include "auxiliary.h"
#include <sys/timeb.h>
#include <unistd.h>



#define TWOPI 6.2831853071795864769252867666
 
void pid(int *i)  {
#ifdef Unix  
  *i = getpid();
#else
  PRINTF("pid not programmed yet");
  *i = 0;
#endif
}
void hostname(char **h, int *i){
#ifdef Unix  
  gethostname(*h,*i);
#else
  PRINTF("pid not programmed yet");
  *i = 0;
#endif
}  


void RandomPermutation(double *x,int n,double *y END_WITH_GSLRNG) {
  /* calculates a random permutation of the n-ector x and returns it in y */
  double *pos;
  int *sort,j;
  
  sort = (int*) malloc(sizeof(int)*n);
  pos = (double*) malloc(sizeof(double)*n);
  for (j=0;j<n;j++) {  
    pos[j]=UNIFORM_RANDOM;
     sort[j]=j; 
  }
  orderdouble(pos,sort,0,n-1);
  for (j=0;j<n;j++) {  y[j] = x[sort[j]]; }
}


void orderdouble(double *d, int *pos, int start, int end END_WITH_GSLRNG)
     /* quicksort algorithm, slightly modified:
        does not sort the data, but d[pos] will be ordered 
	NOTE: pos must have the values 0,1,2,...,start-end !
	(orderdouble is a kind of sorting pos according to
	the variable d)
     */
{
  int randpos, pivot, left, right, pivotpos, swap;
  double Dpivot;

  if( start < end )
  {   
    randpos = start + (int) (UNIFORM_RANDOM * (end-start+1));
    pivot = pos[randpos];
    pos[randpos] = pos[start];
    pos[start] = pivot;
    Dpivot = d[pivot];
   
    pivotpos=start; 
    left = start;
    right=end+1;   

    while (left < right) {      
      while (++left < right && d[pos[left]] < Dpivot) pivotpos++;
      while (--right > left && d[pos[right]] > Dpivot) ;      
      if (left < right) {
        swap=pos[left]; pos[left]=pos[right]; pos[right]=swap;
        pivotpos++;
      }
    }
    pos[start] = pos[pivotpos];
    pos[pivotpos] = pivot;
    orderdouble( d, pos, start, pivotpos-1);
    orderdouble( d, pos, pivotpos + 1, end);
  }
}

#define EPSILON_QUANTILE 0.00000001
 //  different definition than in R !!!!
double quantile(double *X,int lb,double p){
  int *pos,i; 
  double result;
  pos = (int *) malloc(sizeof(int) * lb);
  for (i=0;i<lb;i++) { pos[i] = i; /* needed for orderdouble */  }
  orderdouble(X,pos,0,lb-1);  
  while (!(X[pos[lb-1]]<1.8E308) && (lb>0)) {lb--;}
  
  result = X[pos[(int) (p * (double) lb + EPSILON_QUANTILE)]];
  DD PRINTF("\n pos=%d [%f,%20.20f,%d] ",
	    (int) (p * (double) lb + EPSILON_QUANTILE),
	    p,(p * (double) lb+ EPSILON_QUANTILE),lb);
  DD {for (i=0;i<lb;i++) {PRINTF("X[pos[%d]]=%f  ",i, X[pos[i]]);}}
  free(pos);
  return result;
} 

void gauss(int *n,double *G) 
     /* generates n iid standard Gaussinan variables and returns them in G */
{        
  double UU,VV,sqrttwo;
  int halfn,i;
 
#ifdef RF_GSL
  struct timeb tp; 
  if (RANDOM==NULL) {    
    ftime(&tp);
    gsl_rng_default_seed =  tp.time * 257 + tp.millitm; 
    //gsl_rng_default_seed = 0;    
    RANDOM = gsl_rng_alloc(RANDOMNUMBERGENERATOR);         
  }
#endif

  sqrttwo  = sqrt(2);

  halfn = (*n/2)*2;
  for (i=0;i<halfn;) {
    UU= sqrttwo * sqrt(-log(1.0 - UNIFORM_RANDOM));
    VV = TWOPI * UNIFORM_RANDOM;
    G[i++] = UU * sin(VV);
    G[i++] = UU * cos(VV);
  }
  if (halfn< *n) {
    UU= sqrttwo * sqrt(-log(1.0 - UNIFORM_RANDOM));
    VV = TWOPI * UNIFORM_RANDOM;
    G[i++] = UU * sin(VV);
  }
}
  
