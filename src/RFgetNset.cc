

/* 
 Authors
 Martin Schlather, martin.schlather@cu.lu

 library for simulation of random fields -- get's, set's, and print's

 Copyright (C) 2001 -- 2004 Martin Schlather, 

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
#include <assert.h>
#include <string.h>
#include "RFsimu.h"
#include "MPPFcts.h"
#include <R_ext/Applic.h>


void GetParameterIndices(int *variance, int *scale, 
			 int *kappa, int *lastkappa, int *invscale,
   			 int *aniso)
{
  *variance=VARIANCE; *scale=SCALE; *kappa=KAPPA;
  *lastkappa=LASTKAPPA; *invscale=INVSCALE; *aniso=ANISO;
}


void GetrfParameters(int *covmaxchar, int *methodmaxchar, 
		     int *distrmaxchar,
		     int *covnr, int *methodnr, int *distrnr,
		     int *maxdim, int *maxmodels) {
  if (currentNrCov==-1) InitModelList();
  *covmaxchar=COVMAXCHAR;
  *methodmaxchar=METHODMAXCHAR;
  *distrmaxchar=DISTRMAXCHAR;
  *covnr=currentNrCov;
  *methodnr=(int) Nothing;
  *distrnr=DISTRNR;
  *maxdim=MAXDIM;
  *maxmodels=MAXKEYS;
}


void GetDistrName(int *nr,char **name){
  if ((*nr<0) ||(*nr>=DISTRNR)) {strncpy(*name,"",DISTRMAXCHAR); return;}
  strncpy(*name, DISTRNAMES[*nr], DISTRMAXCHAR);
}

void GetDistrNr(char **name, int *n, int *nr) 
{
  unsigned int ln, v, nn;
  nn = (unsigned int) *n;
  // == -1 if no matching function is found
  // == -2 if multiple matching fnts are found, without one matching exactly
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)

  for (v=0; v<nn; v++) {
    nr[v]=0;
    ln=strlen(name[v]);
    while ( (nr[v] < DISTRNR) && strncmp(name[v], DISTRNAMES[nr[v]], ln)) {
      (nr[v])++;
    }
    if (nr[v]<DISTRNR) { 
      // a matching function is found. Are there other functions that match?
      int j; 
      bool exactmatching,multiplematching;
      exactmatching=(ln==strlen(DISTRNAMES[nr[v]]));
      multiplematching=false;
      j=nr[v]+1; // if two or more distributions have the same name 
      //            the last one is taken 
      while (j<DISTRNR) {
	while ( (j<DISTRNR) && strncmp(name[v],DISTRNAMES[j],ln)) {j++;}
	if (j<DISTRNR) {
	  if (ln==strlen(DISTRNAMES[j])) {nr[v]=j; exactmatching=true;} 
	  else {multiplematching=true;}
	}
	j++;
      } 
      if (!exactmatching && multiplematching) {nr[v]=-2;}
    } else nr[v]=-1;
  }
}

void SetParamCE( int *action, int *force, Real *tolRe, Real *tolIm,
		 int *trials, 
		 int *mmin, // mmin must be vector of length MAXDIM!!
		 int *userfft, int *strategy,
		 ce_param *CE, char *name)
{
  int d;
  switch(*action) {
  case 0 :
    CE->force=(bool)*force;

    CE->tol_re=*tolRe;
   if (CE->tol_re>0) {
      CE->tol_re=0;
      if (GENERAL_PRINTLEVEL>0)
	PRINTF("\nWARNING! %s.tol_re had been positive.\n",name);
    }

    CE->tol_im=*tolIm; 
    if (CE->tol_im<0) {
      CE->tol_im=0; 
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("\nWARNING! %s.tol_im had been neagtive.\n",name);
    }

    CE->trials=*trials; 
    if (CE->trials<1) {
      CE->trials=1;
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("\nWARNING! %s.trials had been less than 1\n",name);
    }

    for (d=0; d<MAXDIM; d++) {
      CE->mmin[d] = mmin[d];
      if (CE->mmin[d]>0) {
	CE->mmin[d] = 
	  1 << (1+(int) (log(((double) CE->mmin[d])-0.1)/log(2.0))); 

	if ((CE->mmin[d]!=mmin[d]) && (GENERAL_PRINTLEVEL>0))
	  PRINTF("\nWARNING! %s.mmin[%d] set to %d.\n",name,d,CE->mmin[d]);
      }
    }
 
    CE->userfft=(bool)*userfft;

    if (*strategy>LASTSTRATEGY) {  
      if (GENERAL_PRINTLEVEL>0) 
	PRINTF("\nWARNING! %s_STRATEGY not set\n",name);
    } else CE->strategy=(char) *strategy;      
    break;
  case 1 :
    *force=CE->force;
    *tolRe=CE->tol_re;
    *tolIm=CE->tol_im;
    *trials=CE->trials;
    for (d=0; d<MAXDIM; d++) {mmin[d]=CE->mmin[d];}   
    *userfft= CE->userfft;
    *strategy=(int) CE->strategy;
    if (GetNotPrint) break;
  case 2 :
    PRINTF("\n%s\n==========\nforce=%d\ntol_re=%e\ntol_im=%e\ntrials=%d\nuserfft=%d\nstrategy=%d\n",name,
	   CE->force,CE->tol_re,CE->tol_im, CE->trials, CE->userfft,
	   (int) CE->strategy);
    for (d=0; d<MAXDIM; d++) PRINTF("%d ", CE->mmin[d]);
    PRINTF("\n");
    break;
  default : PRINTF(" unknown action\n");
  }
}

void GetModelName(int *nr,char **name){
  if (currentNrCov==-1) InitModelList();
  if ((*nr<0) ||(*nr>=currentNrCov)) {strncpy(*name,"",COVMAXCHAR); return;}
  strncpy(*name,CovList[*nr].name,COVMAXCHAR);
}

void GetNrParameters(int *covnr, int *n, int* kappas) {
  int v;
  if (currentNrCov==-1) InitModelList();
  for (v=0; v<*n; v++) {
    if ((covnr[v]<0) ||(covnr[v]>=currentNrCov)) {kappas[v]=-1;}
    else kappas[v] = CovList[covnr[v]].kappas;
  }
}
void GetModelNr(char **name, int *n, int *nr) 
{
  unsigned int ln, v, nn;
  nn = (unsigned int) *n;
  // == -1 if no matching function is found
  // == -2 if multiple matching fnts are found, without one matching exactly
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)
  if (currentNrCov==-1) InitModelList();

  for (v=0; v<nn; v++) {
    nr[v]=0;
    ln=strlen(name[v]);
    while ( (nr[v]<currentNrCov) && strncmp(name[v],CovList[nr[v]].name,ln)) {
      (nr[v])++;
    }
    if (nr[v]<currentNrCov) { 
      // a matching function is found. Are there other functions that match?
      int j; 
      bool exactmatching,multiplematching;
      exactmatching=(ln==strlen(CovList[nr[v]].name));
      multiplematching=false;
      j=nr[v]+1; // if two or more covariance functions have the same name 
      //            the last one is taken 
      while (j<currentNrCov) {
	while ( (j<currentNrCov) && strncmp(name[v],CovList[j].name,ln)) {j++;}
	if (j<currentNrCov) {
	  if (ln==strlen(CovList[j].name)) {nr[v]=j; exactmatching=true;} 
	  else {multiplematching=true;}
	}
	j++;
      } 
      if (!exactmatching && multiplematching) {nr[v]=-2;}
    } else nr[v]=-1;
  }
}

void GetMethodNr(char **name, int *n, int *nr) 
{
  unsigned int ln, v, nn;
  nn = (unsigned int) *n;
  // == -1 if no matching method is found
  // == -2 if multiple matching methods are found, without one matching exactly
  for (v=0; v<nn; v++) {
    nr[v]=0;
    ln=strlen(name[v]);
    while ( (nr[v]<(int) Forbidden) && strncmp(name[v],METHODNAMES[nr[v]],ln)) {
      (nr[v])++;
    }
    if (nr[v]<(int) Forbidden) { 
      // a matching method is found. Are there other methods that match?
      int j; 
      bool exactmatching,multiplematching;
      exactmatching=(ln==strlen(METHODNAMES[nr[v]]));
      multiplematching=false;
      j=nr[v]+1; // if two or more methods have the same name the last one is
      //          taken; stupid here, but nice in GetCovNr
      while (j<(int) Forbidden) {
	while ( (j<(int) Forbidden) && strncmp(name[v],METHODNAMES[j],ln)) j++;
	if (j<(int) Forbidden) {
	  if (ln==strlen(METHODNAMES[j])) {nr[v]=j; exactmatching=true;} 
	  else {multiplematching=true;}
	}
	j++;
      } 
      if (!exactmatching && multiplematching) {nr[v]=-2;}
    } else nr[v]=-1;
  }
}

void GetMethodName(int *nr, char **name) 
{   
  if ((*nr<0) ||(*nr>=(int) Forbidden)) {
    strncpy(*name,"",METHODMAXCHAR); 
    return;
  }
  strncpy(*name,METHODNAMES[*nr],METHODMAXCHAR);
}

void PrintMethods()
{ 
  int i;
  PRINTF("\n\n  Methods for generating Gaussian random fields\n  =============================================\n\n");  
  for (i=0;i<(int) Nothing;i++) { PRINTF("  * %s\n",METHODNAMES[i]); }
  PRINTF("\n\n  Methods for non-Gaussian random fields\n  ======================================\n");
  for (i=1+(int) Nothing;i<(int) Forbidden; i++) { 
    PRINTF("  * %s\n",METHODNAMES[i]); }
  PRINTF("\n");
}

void PrintCovDetailed(int i)
{
  PRINTF("%5d: %s cov=%d %d %d, tbm=%d %d %d, spec=%d,\n        abf=%d hyp=%d, oth=%d %d\n",
	 i,CovList[i].name,CovList[i].cov,CovList[i].cov_loc,
	 CovList[i].square_factor,
	 CovList[i].cov_tbm2,CovList[i].cov_tbm3,CovList[i].tbm_method,
	 CovList[i].spectral,CovList[i].add_mpp_scl,CovList[i].hyperplane,
	 CovList[i].initother,CovList[i].other);
}

void PrintModelListDetailed()
{
  int i;
  if (CovList==NULL) {PRINTF("There are no functions available!\n");} 
  else {
    PRINTF("\nList of covariance functions:\n=============================\n");
    for (i=0;i<currentNrCov;i++) PrintCovDetailed(i);
  }
}

void GetRange(int *nr, int *dim, int *index, Real *range, int *lrange){
  // never change Real without crosschecking with fcts in RFCovFcts.cc!
  // index is increased by one except index is the largest value possible
  int i;
  if (currentNrCov==-1) InitModelList();
  if ((*nr<0) || (*nr>=currentNrCov) || 
      (*lrange != CovList[*nr].kappas * 4) ) {
    for (i=0; i<*lrange; i++) range[i]=RF_NAN;
    *index = -100;
    return;
  }
  assert(CovList[*nr].range != NULL);
  CovList[*nr].range(*dim, index, range);
  // index>0 : further parts of the range are missing
  // index=-1: no further part of the range
  // index=-2: dimension not valid
  //assert(false);
}

extern int IncludeModel(char *name, int kappas, checkfct check,
			int isotropic, bool variogram, methodfct method,
			rangefct range) 
{  
  int i;
  if (currentNrCov==-1) InitModelList(); assert(CovList!=NULL);
  assert(currentNrCov>=0);
  if (currentNrCov>=MAXNRCOVFCTS) 
    {PRINTF("Error. Model list full.\n");return -1;}
   strncpy(CovList[currentNrCov].name, name, COVMAXCHAR-1);
  CovList[currentNrCov].name[COVMAXCHAR-1]='\0';
  if (strlen(name)>=COVMAXCHAR) {
    PRINTF("Warning! Covariance name is truncated to `%s'",
	   CovList[currentNrCov].name);
  }
  CovList[currentNrCov].kappas = kappas;
  CovList[currentNrCov].check = check;
  CovList[currentNrCov].isotropic = isotropic;
  CovList[currentNrCov].cov=NULL;
  CovList[currentNrCov].naturalscale=NULL;
  CovList[currentNrCov].cov_loc=NULL;
  CovList[currentNrCov].cov_tbm2=NULL;
  CovList[currentNrCov].cov_tbm3=NULL;
  CovList[currentNrCov].ableitung=NULL,
  CovList[currentNrCov].spectral=NULL;
  CovList[currentNrCov].add_mpp_scl=NULL;  
  CovList[currentNrCov].add_mpp_rnd=NULL;  
  CovList[currentNrCov].hyperplane=NULL;
  CovList[currentNrCov].initother=NULL;
  CovList[currentNrCov].other=NULL;
  CovList[currentNrCov].square_factor=NULL;
  CovList[currentNrCov].tbm_method = Forbidden;
  CovList[currentNrCov].isotropic = isotropic;
  CovList[currentNrCov].variogram = variogram;
  CovList[currentNrCov].even = true; 
  CovList[currentNrCov].range= range;
  assert((CovList[currentNrCov].method=method)!=NULL); 
   for (i=0; i<MAXDIM; i++) {
     CovList[currentNrCov].odd[i] = false;
   }
  if (kappas>0) assert(check!=NULL);
  currentNrCov++;
  return(currentNrCov-1);
}

extern void addodd(int nr, int dim) {
  assert((dim>=0) && (dim<MAXDIM));
  assert((nr>=0) && (nr<currentNrCov));
  assert(CovList[nr].isotropic==ANISOTROPIC);
  PRINTF("addodd not implemented yet. Sorry.");
  assert(false);
  CovList[nr].even = false;
  CovList[nr].odd[dim] = true;
}
			
//extern void addSimu(int nr, SimulationType r1, SimulationType r2,
//		    SimulationType r3)
//  // see RFsimu.h for comments
//{ 
//  int i;
//  assert((nr>=0) && (nr<currentNrCov));
//  CovList[nr].first[0]=r1;
//  CovList[nr].first[1]=r2;
//  for (i=2; i<MAXDIM; i++) CovList[nr].first[i]=r3;
//}

extern void addCov(int nr, covfct cov, natscalefct naturalscale, 
		   covfct cov_loc, scalefct square_factor)
{
  assert((nr>=0) && (nr<currentNrCov));
  CovList[nr].cov=cov;
  CovList[nr].naturalscale=naturalscale;
  CovList[nr].cov_loc=cov_loc;
  CovList[nr].square_factor=square_factor;
  if ((square_factor!=NULL) || (cov_loc!=NULL) )
    assert( (square_factor!=NULL) && (cov_loc!=NULL) );
}

extern void addTBM(int nr, isofct cov_tbm2, isofct cov_tbm3,
		   isofct ableitung,
		   SimulationType tbm_method,
		   randommeasure spectral) {
  // must be called always AFTER addCov !!!!
  assert((nr>=0) && (nr<currentNrCov));
  CovList[nr].cov_tbm2=cov_tbm2;
  CovList[nr].cov_tbm3=cov_tbm3;
  CovList[nr].ableitung=ableitung,
  CovList[nr].tbm_method = tbm_method;
  CovList[nr].spectral=spectral;

  if ( (ableitung!=NULL) || (tbm_method!=Nothing)) 
    assert((ableitung!=NULL) && (CovList[nr].cov!=NULL)
	   && (tbm_method!=Nothing));
}
	
extern void addOther(int nr, initmppfct add_mpp_scl, MPPRandom add_mpp_rnd, 
		     generalSimuInit initother, generalSimuMethod other){
  assert((nr>=0) && (nr<currentNrCov));
  CovList[nr].add_mpp_scl=add_mpp_scl;  
  CovList[nr].add_mpp_rnd=add_mpp_rnd;  
  CovList[nr].initother=initother;
  CovList[nr].other=other;
  CovList[nr].hyperplane=NULL;
 if ((add_mpp_scl!=NULL) || (add_mpp_rnd!=NULL)) 
    assert((add_mpp_scl!=NULL) && (add_mpp_rnd!=NULL));
  if ((other!=NULL) || (initother!=NULL)) 
    assert((other!=NULL) && (initother!=NULL));
}

void PrintModelList()
{
  int i, tbm2;
  char percent[]="%";
  char empty[]="";
  char header[]="Circ local TBM2 TBM3 sp dir add hyp oth\n";
  char coded[3][2]={"-","X","+"};
  char firstcolumn[20],line[80];
  if (currentNrCov==-1) {
    InitModelList(); 
    if (GENERAL_PRINTLEVEL>5) 
      PRINTF("List of covariance functions initiated.\n"); 
  }
  if (CovList==NULL) {PRINTF("There are no functions available!\n");} 
  else {
    sprintf(firstcolumn,"%s%ds ",percent,COVMAXCHAR);
    sprintf(line,"%s3s  %s3s   %s3s  %s3s  %s2s %s2s  %s2s  %s2s  %s2s\n",
	    percent,percent,percent,percent,percent,percent,percent,percent,
	    percent);
    PRINTF("\n\n");
    PRINTF(firstcolumn,empty); PRINTF("       List of models\n");
    PRINTF(firstcolumn,empty); PRINTF("       ==============\n");
    PRINTF(firstcolumn,empty); PRINTF("[See also PrintMethodList()]\n\n");
    PRINTF(firstcolumn,empty);PRINTF(header,empty);   
    for (i=0;i<currentNrCov;i++) {
      PRINTF(firstcolumn,CovList[i].name);
      if (strcmp(CovList[i].name,"nugget")==0) {
	PRINTF(line,
	       coded[2], coded[false], coded[2],     coded[2],     coded[2],
	       coded[2], coded[2], coded[false], coded[false]);
      } else {
	if (!(tbm2 = CovList[i].cov_tbm2!=NULL))
	  tbm2 = 2 * (CovList[i].ableitung != NULL);
	PRINTF(line,
	       coded[(CovList[i].cov!=NULL) && (CovList[i].cov_loc==NULL)],
	       coded[CovList[i].cov_loc!=NULL],
	       coded[tbm2],
	       coded[CovList[i].cov_tbm3!=NULL],
	       coded[CovList[i].spectral!=NULL], 
	       coded[(CovList[i].cov!=NULL) && (CovList[i].cov_loc==NULL)],
	       // nugget is left out
	       coded[CovList[i].add_mpp_scl!=NULL],
	       coded[CovList[i].hyperplane!=NULL],
	       coded[CovList[i].initother!=NULL]
	       );
      }
    }
    PRINTF(firstcolumn,empty);PRINTF(header,empty);   
  }
}

void GetModelList(int* idx) {
  int i, j, methods;
  if (currentNrCov==-1) {
    InitModelList(); 
    if (GENERAL_PRINTLEVEL>5) 
      PRINTF("List of covariance functions initiated.\n"); 
  }
  if (CovList==NULL) return;
  assert(Nothing==10);
  for (j=i=0; i<currentNrCov; i++) {
    idx[j++] = (CovList[i].cov!=NULL) && (CovList[i].cov_loc==NULL);
    idx[j++] = CovList[i].cov_loc!=NULL;
    idx[j++] = CovList[i].cov_tbm2!=NULL ||  CovList[i].ableitung!=NULL;
    idx[j++] = CovList[i].cov_tbm3!=NULL;
    idx[j++] = CovList[i].spectral!=NULL;
    idx[j++] = (CovList[i].cov!=NULL) && (CovList[i].cov_loc==NULL);
    // nugget is left out;
    idx[j++] = CovList[i].add_mpp_scl!=NULL;
    idx[j++] = CovList[i].hyperplane!=NULL;
    idx[j++] = CovList[i].initother!=NULL;
  }
  return;
}

void GetNaturalScaling(int *covnr, Real *q, int *naturalscaling,
		       Real *natscale, int *error)
{
  //  q: kappas only

  // values of naturalscaling:
  //#define NATSCALE_EXACT 1   
  //#define NATSCALE_APPROX 2
  //#define NATSCALE_MLE 3 /* check mleRF when changing !! */
  // +10 if numerical is allowed
  static int oldcovnr = -99;
  static Real oldp[TOTAL_PARAM], p[TOTAL_PARAM];
  static Real OldNatScale;
  int TEN, endfor, k;
  bool numeric;  
  *error = 0;
  TEN = 10;

  endfor = KAPPA + CovList[*covnr].kappas;
  for (k=KAPPA; k<endfor; k++) p[k]=q[k-KAPPA];
  if (*naturalscaling) {	 
    if (CovList[*covnr].isotropic!=FULLISOTROPIC) {
      *error=ERRORANISOTROPIC; return;}
    if (!(numeric=*naturalscaling>TEN)) {TEN=0;}
    *natscale=0.0;
    if (CovList[*covnr].naturalscale!=NULL) { 
      *natscale = CovList[*covnr].naturalscale(p,*naturalscaling-TEN);
      if (*natscale!=0.0) return;
    }
    if (numeric) {
      Real x,newx,yold,y,newy;
      covfct cov;
      int parami, wave,i;
      if ((cov=CovList[*covnr].cov)==NULL) {*error=ERRORNOTDEFINED;return;}
      if ((cov=CovList[*covnr].cov)==nugget) {*error=ERRORRESCALING; return;}
      
      // already calculated ?
      parami=KAPPA; // do not compare mean,variance, etc.
      if (oldcovnr==*covnr) {
	for (; parami<endfor; parami++) {
	  if (oldp[parami]!=p[parami]) {
	    break;
	  }
	}
	if (parami==endfor) {
	  *natscale=OldNatScale; 
	  return;
	}
      }   


      // the other parameters need not to be copied as checked that 
      // they are identical to previous ones
      for (;parami<endfor;parami++) {oldp[parami]=p[parami];}
      oldp[VARIANCE] = oldp[SCALE] = oldp[INVSCALE] = 1.0; 
      oldcovnr = -99; /* if error occurs, the next call will realise that 
			 the previous result cannot be used */

      /* **************************************
	 Now, find a quick and good solution for NewInvScale --
      */
      wave  = 0;
      x = 1.0; 
      if ( (yold=cov(&x,oldp,1)) > 0.05) {
	Real leftx;
	x *= 2.0;
	while ( (y=cov(&x,oldp,1)) > 0.05) {  
	  if (yold<y){ wave++;  if (wave>10) {*error=ERRORWAVING; return;} }
	  yold = y;
	  x *= 2.0;
	  if (x>1E30) {*error=ERRORRESCALING; return;} // or use a separate ERROR
	} 
	leftx = x * 0.5;
	for (i=0; i<3 /* good choice?? */ ;i++) {          
	  if (y==yold) {*error=ERRORWAVING; return;} // should never appear
	  newx = x + (x-leftx)/(y-yold)*(0.05-y);
	  if ( (newy=cov(&newx,oldp,1)) >0.05) {
	    leftx = newx;
	    yold  = newy;
	  } else {
	    x = newx;
	    y = newy;
	  }
	}
	if (y==yold)  {*error=ERRORWAVING; return;} // should never appear
	*natscale = 1.0 / ( x + (x-leftx)/(y-yold)*(0.05-y) );
      } else {
	Real rightx;
	x *= 0.5;
	while ( (y=cov(&x,oldp,1)) < 0.05) {  
	  if (yold>y){ wave++;  if (wave>10) {*error=ERRORWAVING; return;} }
	  yold = y;
	  x *= 0.5;
	  if (x<1E-30) {*error=ERRORRESCALING; return;} //or use a separate ERROR
	}    
	rightx = x * 2.0;
	for (i=0; i<3 /* good choice?? */ ;i++) {          
	  if (y==yold) {*error=ERRORWAVING; return;} // should never appear
	  newx = x + (x-rightx)/(y-yold)*(0.05-y);
	  if ( (newy=cov(&newx,oldp,1)) <0.05) {
	    rightx = newx;
	    yold   = newy;
	  } else {
	    x = newx;
	    y = newy;
	  }
	}
	if (y==yold)  {*error=ERRORWAVING; return;} // should never appear
	*natscale = 1.0 / ( x + (x-rightx)/(y-yold)*(0.05-y) );
      }
      oldcovnr=*covnr;  /* result is reliable -> stored for next call */
      OldNatScale = *natscale;
    } else {*error=ERRORRESCALING; }
  } else {
    *natscale = 1.0;
  }
}
void SetParam(int *action, int *storing,int *printlevel,int *naturalscaling,
	      char **pch) {
  switch (*action) {
  case 0 :
    GENERAL_STORING= (bool) *storing;
    GENERAL_PRINTLEVEL=*printlevel;
    GENERAL_NATURALSCALING=*naturalscaling;
    if ((strlen(*pch)>1) && (GENERAL_PRINTLEVEL>0)) 
	PRINTF("\n`pch' has more than one character -- first character taken only\n");
    strncpy(GENERAL_PCH, *pch, 1);
    break;
  case 1 :
    *storing =  GENERAL_STORING;
    *printlevel = GENERAL_PRINTLEVEL;
    *naturalscaling = GENERAL_NATURALSCALING;
    strcpy(*pch, GENERAL_PCH);
   if (GetNotPrint) break;
 case 2 : 
    PRINTF("\nGeneral Parameters\n==================\nstoring=%d\nprint level=%d\nnatural scaling=%d\npch=%s\n",
	   GENERAL_STORING,GENERAL_PRINTLEVEL,GENERAL_NATURALSCALING,*pch);
    break;
  default : PRINTF(" unknown action\n"); 
  }
}

void StoreTrend(int *keyNr, int *modus, char **trend, 
		int *lt, Real *lambda, int *ll, int *error)
{
  unsigned long len;
  key_type *key;
  if (currentNrCov==-1) InitModelList(); assert(CovList!=NULL); 
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*error=ERRORREGISTER;goto ErrorHandling;}
  key = &(KEY[*keyNr]);
  if (key->LinearTrend!=NULL) free(key->LinearTrend);
  key->LinearTrend = NULL;
  if (key->TrendFunction!=NULL) free(key->TrendFunction);
  key->TrendFunction = NULL;
  key->lTrendFct = 0;
  key->lLinTrend = 0;

  switch(*modus){
  case 0 : break;
  case 3 : case 1 : 
    key->LinearTrend = (Real *) malloc(len = sizeof(Real) * *ll);
    memcpy(key->LinearTrend, lambda, len);
    if (*modus==1) break; 
    key->lLinTrend = *ll;
  case 2 : 
    key->TrendFunction = (char *) malloc(len = sizeof(char) * *lt);
    memcpy(key->TrendFunction, *trend, len);
    key->lTrendFct = *lt;
    break; 
  default: key->TrendModus = -1; *error=ERRORUNSPECIFIED; goto ErrorHandling;
  }
  key->TrendModus = *modus;
  *error = 0;
  return;

 ErrorHandling:
  return;
}


void GetTrueDim(bool anisotropy, int timespacedim, Real* param, 
	        int *TrueDim, bool *no_last_comp,
		int *start_aniso, int *index_dim)
  // output: start_aniso (where the non-zero positions start, anisotropy)
  //         index_dim   (dimension nr of the non-zero positions, aniso)
  //         no_last_comp(last dim of aniso matrix. is identically zero, i.e,
  //                      effectively no time component)
  //         TrueDim     The reduced dimension, including the time direction
{ 
  int i,j;
  long unsigned int endfor, startfor, k;
  assert(timespacedim>0);
  if (anisotropy) {
    /* check whether the dimension can be reduced */
    /* TODO:this algorithm only detects zonal anisotropies -- could be improved*/
    for (j=-1, i=0; i<timespacedim; i++) {
      startfor = i * timespacedim + ANISO;
      endfor= startfor + timespacedim;
      for (k=startfor; k<endfor; k++) {
	if (param[k]!=0.0) break;
      }
      if (k<endfor) {
	start_aniso[++j]=startfor; /* where the non-zero positions*/
 	/* of the anisotropy matrix start */
 	index_dim[j] = i;
      }
    }
    *no_last_comp = k>=endfor; // last i is time, if k=endfor then all 
    //                            components zero
    *TrueDim = j+1;
  } else { 
    *no_last_comp = true;
    *TrueDim = timespacedim; // == spatialdim
  }
}


void GetTrendLengths(int *keyNr, int *modus, int *lt, int *ll, int *lx)
{
  key_type *key;
  if (currentNrCov==-1) {*modus=-1;goto ErrorHandling;}
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*modus=-2; goto ErrorHandling;}
  key = &(KEY[*keyNr]);
  if (!key->active) {*modus=-3; goto ErrorHandling;}
  *modus = key->TrendModus;
  *lt = key->lTrendFct;
  *ll = key->lLinTrend;
  *lx = key->totalpoints;
  return;
 ErrorHandling:
  return;
}


void GetTrend(int *keyNr, char **trend, Real *lambda, Real *x, int *error)
{
  key_type *key;
  if (currentNrCov==-1) {*error=-1;goto ErrorHandling;}
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*error=-2; goto ErrorHandling;}
  key = &(KEY[*keyNr]);
  if (!key->active) {*error=-3; goto ErrorHandling;}
  switch(key->TrendModus){
  case 0 : break;
  case 3 : case 1 : 
    memcpy(lambda, key->LinearTrend, key->lLinTrend);
    if (key->TrendModus==1) break; 
  case 2:
    strcpy(*trend, key->TrendFunction);
    break; 
  default: *error=-4; goto ErrorHandling;
  }
  *error = 0;
  return;

 ErrorHandling:
  return;
}


void GetCoordinates(int *keyNr, Real *x, int *error)
{
  key_type *key;
  if (currentNrCov==-1) {*error=-1;goto ErrorHandling;}
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*error=-2; goto ErrorHandling;}
  key = &(KEY[*keyNr]);
  if (!key->active) {*error=-3; goto ErrorHandling;}
  if (key->TrendModus==0) {*error=-4; goto ErrorHandling;}
 
  *error = ERRORNOTPROGRAMMED;  // TODO !!

  return;

 ErrorHandling:
  return;
}


void GetKeyInfo(int *keyNr, int *total, int *lengths, int *spatialdim, 
		int *timespacedim, int *grid, int *distr, int *maxdim)
{
  int d;
  key_type *key;
  if (*maxdim<MAXDIM) { *total=-3; *maxdim=MAXDIM; return;} 
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*total=-1; return;} 
  key = &(KEY[*keyNr]);  
  if (!key->active) {*total=-2;} 
  else { 
    *total=key->totalpoints;
    *spatialdim = key->spatialdim;
    *timespacedim = key->timespacedim;
    for (d=0; d<MAXDIM; d++) lengths[d]=key->length[d];
    *grid = (int) key->grid;
    *distr = (int) key->distribution; 
    *maxdim = MAXDIM;
   }
}


void GetCornersOfGrid(key_type  *key, int Stimespacedim, int* start_aniso, 
		      Real *param, Real *sxx){
  // returning sequence (#=2^MAXDIM) vectors of dimension s->simutimespacedim
  Real sx[ZWEIHOCHMAXDIM * MAXDIM];
  int index,i,k,g,endforM1;
  long endfor,indexx;
  char j[MAXDIM];      
  endfor = 1 << key->timespacedim;
  endforM1 = endfor -1;
  for (i=0; i<key->timespacedim; i++) j[i]=0;
  for (index=i=0; i<endfor; i++) {
    for (k=0; k<key->timespacedim; k++) {
      sx[index] = key->x[k][XSTART];
      if (j[k]) 
	sx[index] += key->x[k][XSTEP] * (double) (key->length[k]-1);
      index++;
    }
    k=0; (j[k])++;
    if (i!=endforM1)
      while(j[k]>=2) { assert(j[k]==2); j[k]=0;  k++; (j[k])++;}	
  }
  for (index=indexx=i=0; i<endfor; 
       i++, index+=key->timespacedim, indexx+=Stimespacedim) {
    for (k=0; k<Stimespacedim; k++) {
      register Real dummy;
      dummy = 0.0;
      for (g=0; g<key->timespacedim; g++) 
	dummy += param[start_aniso[k]+g] * sx[index+g];
      sxx[indexx + k] = dummy;
    }
  }
}

void GetCornersOfElement(key_type  *key, int Stimespacedim, int* start_aniso, 
		      Real *param, Real *sxx){
  Real sx[ZWEIHOCHMAXDIM * MAXDIM];
  int index,i,k,g,endforM1;
  long endfor,indexx;
  char j[MAXDIM];      
  endfor = 1 << key->timespacedim;
  endforM1 = endfor -1;
  for (i=0; i<key->timespacedim; i++) j[i]=0;
  for (index=i=0; i<endfor; i++) {
    for (k=0; k<key->timespacedim; k++) {
      if (j[k]) sx[index] = key->x[k][XSTEP]; else sx[index]=0.0;
      index++;
    }
    k=0; (j[k])++;
    if (i!=endforM1)
      while(j[k]>=2) { assert(j[k]==2); j[k]=0;  k++; (j[k])++;}	
  }
  for (index=indexx=i=0; i<endfor; 
       i++, index+=key->timespacedim, indexx+=Stimespacedim) {
    for (k=0; k<Stimespacedim; k++) {
      register Real dummy;
      dummy = 0.0;
      for (g=0; g<key->timespacedim; g++) 
	dummy += param[start_aniso[k]+g] * sx[index+g];
      sxx[indexx + k] = dummy;
    }
  }
}

void GetRangeCornerDistances(key_type *key, Real *sxx, int Stimespacedim,
			int Ssimuspatialdim, Real *min, Real *max) {
  long endfor,indexx,index;
  int d,i,k;
  endfor = 1 << key->timespacedim;
  *min = RF_INF;
  *max = RF_NEGINF;
  for (index=i=0; i<endfor; i++, index+=Stimespacedim)
    for (indexx=k=0; k<i; k++, indexx+=Stimespacedim) {
      register double sum = 0.0, dummy;
      for (d=0; d<Ssimuspatialdim; d++) {
	dummy = sxx[index+d] - sxx[indexx+d];
	sum += dummy * dummy;
      }
      if ((sum > 0) && (sum < *min)) *min = sum; 
      if ((sum > 0) && (sum > *max)) *max = sum; 
    }
  *min = sqrt(*min);
  *max = sqrt(*max);
}

void niceFFTnumber(int *n) {
  *n = (int) NiceFFTNumber( (unsigned long) *n);
}
bool ok_n(int n, int *f, int nf) // taken from fourier.c of R
{
  int i;
  for (i = 0; i < nf; i++)
    while(n % f[i] == 0) if ((n = n / f[i]) == 1) return TRUE;
  return n == 1;
}
int nextn(int n, int *f, int nf) // taken from fourier.c of R
{
  while(!ok_n(n, f, nf)) n++;
  return n;
}

bool HOMEMADE_NICEFFT=true;
unsigned long NiceFFTNumber(unsigned long n) {
  unsigned long i,ii,j,jj,l,ll,min, m=1, f[4]={2,3,5,7};
  if (HOMEMADE_NICEFFT) {
    if (n<=1) return n;
    for (i=0; i<4; i++) 
      while ( ((n % f[i])==0) && (n>10000)) { m*=f[i]; n/=f[i]; }
    if (n>10000) {
      while (n>10000) {m*=10; n/=10;}
      n++;
    }
    min = 10000000;
    for (i=0, ii=1; i<=14; i++, ii<<=1) { // 2^14>10.000
      if (ii>=n) { if (ii<min) min=ii; break; }
      for (j=0, jj=ii; j<=9; j++, jj*=3) {// 3^9>10.000
	if (jj>=n) { if (jj<min) min=jj; break; }
	
	//for (k=0, kk=jj; k<=6; k++, kk*=5) {// 5^6>10.000
	//if (kk>=n) { if (kk<min) min=kk; break; }
	//for (l=0, ll=kk; l<=5; l++, ll*=7) {// 7^5>10.000
	
	// instead of (if 7 is included)
	for (l=0, ll=jj; l<=6; l++, ll*=5) {// 5^5>10.000
	  	  
	  if (ll>=n) { 
	    if (ll<min) min=ll;
	    break;
	  }
	  //}
	}
      }
    }
    return m*min;
  } else { // not HOMEMADE_NICEFFT
#define F_NUMBERS 3
    int f[F_NUMBERS]={2,3,5}; 
    return nextn(n, f, F_NUMBERS);
  }
}

int eigenvalues(Real *C, int dim, double *ev)
{
  // SVD decomposition
  double *V,*e, *U, *G, *D;
  long longrow = dim, job=11;
  int error, d;
  error=0;

  D = V = e = U = G = NULL;
  if ((U =(double *) malloc(sizeof(double) * longrow * longrow))==NULL){
    error=ERRORMEMORYALLOCATION;goto ErrorHandling;
  }
  if ((G = (double *)  malloc(sizeof(double) * longrow))==NULL){
    error=ERRORMEMORYALLOCATION;goto ErrorHandling;
  } 
  if ((D = (double *)  malloc(sizeof(double) * longrow  * longrow))==NULL){
    error=ERRORMEMORYALLOCATION;goto ErrorHandling;
  } 
  if ((V =(double *) malloc(sizeof(double) * longrow * longrow))==NULL){
    error=ERRORMEMORYALLOCATION;goto ErrorHandling;
  }
  if ((e = (double *) malloc(sizeof(double) * longrow))==NULL){
    error=ERRORMEMORYALLOCATION;goto ErrorHandling;
  }
  
  for (d=dim * dim -1; d>=0; d--) D[d]=C[d]; // dsvdc destroys the 
  // input matrix !!!!!!!!!!!!!!!!!!!!
  { 
    int row,jobint; 
    row=longrow; jobint=job;
    F77_CALL(dsvdc)(D,&row,&row,&row,ev,e,U,&row,V,&row,G,&jobint,
		    &error);
  }

 ErrorHandling:
  if (D!=NULL) free(D);
  if (V!=NULL) free(V);
  if (e!=NULL) free(e);
  if (U!=NULL) free(U);
  if (G!=NULL) free(G);
  return error;
}

