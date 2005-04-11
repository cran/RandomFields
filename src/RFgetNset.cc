

/* 
 Authors
 Martin Schlather, schlath@hsu-hh.de

 library for simulation of random fields -- get's, set's, and print's

 Copyright (C) 2001 -- 2004 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
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

extern "C" {
#include <R.h>
#include <Rdefines.h>
}
#include <R_ext/Applic.h>
#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "RFsimu.h"
#include "RFCovFcts.h"
#include "MPPFcts.h"


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

void SetParamDecision( int *action, int *stationary_only, int *exactness)
{
  int NA=-1;
  switch(*action) {
  case 0 :
     DECISION_PARAM.stationary_only = (*stationary_only==NA) ?
       DECISION_CASESPEC : *stationary_only ? DECISION_TRUE : DECISION_FALSE;
     DECISION_PARAM.exactness = (*exactness==NA) ?
       DECISION_CASESPEC : *exactness ? DECISION_TRUE : DECISION_FALSE;
    break;
  case 1 :
    *stationary_only = DECISION_PARAM.stationary_only==DECISION_TRUE ? true :
      DECISION_PARAM.stationary_only==DECISION_FALSE ? false : NA;
    *exactness = DECISION_PARAM.exactness==DECISION_TRUE ? true :
      DECISION_PARAM.exactness==DECISION_FALSE ? false : NA;
    if (GetNotPrint) break;
  case 2 :
    PRINTF("\nDECISION_PARAM\n==================\nstationary_only=%d\nexactness=%d\n",
	   DECISION_PARAM.stationary_only, DECISION_PARAM.exactness);
    break;
  default : PRINTF(" unknown action\n");
  }
}

void SetParamCE( int *action, int *force, double *tolRe, double *tolIm,
		 int *trials, 
		 int *mmin, // mmin must be vector of length MAXDIM!!
		 int *userfft, int *strategy, double *maxmem,
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

    CE->maxmem = *maxmem;

    break;
  case 1 :
    *force=CE->force;
    *tolRe=CE->tol_re;
    *tolIm=CE->tol_im;
    *trials=CE->trials;
    for (d=0; d<MAXDIM; d++) {mmin[d]=CE->mmin[d];}   
    *userfft= CE->userfft;
    *strategy = (int) CE->strategy;
    *maxmem = CE->maxmem;
    if (GetNotPrint) break;
  case 2 :
    PRINTF("\n%s\n==========\nforce=%d\ntol_re=%e\ntol_im=%e\ntrials=%d\nuserfft=%d\nstrategy=%d\nmaxmem=%e\n",name,
	   CE->force,CE->tol_re,CE->tol_im, CE->trials, CE->userfft,
	   (int) CE->strategy, CE->maxmem);
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
  // == -2 if multiple matching fctns are found, without one matching exactly
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
  PRINTF("%5d: %s cov=%d, intr=%d, cut=%d tbm=%d %d %d, spec=%d,\n        abf=%d hyp=%d, oth=%d %d\n",
	 i,CovList[i].name,CovList[i].cov,
	 CovList[i].intrinsic_strategy,CovList[i].cutoff_strategy,
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

void GetRange(int *nr, int *dim, int *index, double *range, int *lrange){
  // never change double without crosschecking with fcts in RFCovFcts.cc!
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
			int isotropic, bool variogram, infofct info,
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
  for (i=0; i<SimulationTypeLength; i++)
    CovList[currentNrCov].implemented[i] = NOT_IMPLEMENTED;
  CovList[currentNrCov].intrinsic_strategy=NULL;
  CovList[currentNrCov].cutoff_strategy=NULL;
  CovList[currentNrCov].derivative=NULL;
  CovList[currentNrCov].secondderivt=NULL;

  CovList[currentNrCov].cov_tbm2=NULL;
  CovList[currentNrCov].cov_tbm3=NULL;
  CovList[currentNrCov].derivative=NULL,
  CovList[currentNrCov].spectral=NULL;
  CovList[currentNrCov].add_mpp_scl=NULL;  
  CovList[currentNrCov].add_mpp_rnd=NULL;  
  CovList[currentNrCov].hyperplane=NULL;
  CovList[currentNrCov].initother=NULL;
  CovList[currentNrCov].other=NULL;
  CovList[currentNrCov].tbm_method = Forbidden;
  CovList[currentNrCov].isotropic = isotropic;
  CovList[currentNrCov].variogram = variogram;
  CovList[currentNrCov].even = true; 
  CovList[currentNrCov].range= range;
  assert((CovList[currentNrCov].info=info)!=NULL); 
   for (i=0; i<MAXDIM; i++) {
     CovList[currentNrCov].odd[i] = false;
   }
  if (kappas>0) assert(check!=NULL);
  currentNrCov++;
  return(currentNrCov-1);
}

extern void addodd(int nr, int dim) {
  PRINTF("addodd not implemented yet. Sorry.");
  assert(false);

  assert((dim>=0) && (dim<MAXDIM));
  assert((nr>=0) && (nr<currentNrCov));
  assert(CovList[nr].isotropic==ANISOTROPIC);
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

extern void addCov(int nr, covfct cov, isofct derivative,
		   natscalefct naturalscale)
{
  assert((nr>=0) && (nr<currentNrCov) && cov!=NULL);
  CovList[nr].cov=cov;
  CovList[nr].derivative=derivative;
  CovList[nr].naturalscale=naturalscale;
  CovList[nr].implemented[CircEmbed] =
    CovList[nr].implemented[Direct] = !CovList[nr].variogram;
}

extern void addLocal(int nr, IntrinsicStrategyFct intrinsic_strategy, 
		     CutoffStrategyFct cutoff_strategy, 
		     isofct secondderiv)
{
  assert((nr>=0) && (nr<currentNrCov) && CovList[nr].derivative!=NULL);
  CovList[nr].intrinsic_strategy=intrinsic_strategy;
  CovList[nr].cutoff_strategy=cutoff_strategy;
  CovList[nr].implemented[CircEmbedCutoff] = cutoff_strategy!=NULL;
  CovList[nr].secondderivt=secondderiv;
  if ((CovList[nr].implemented[CircEmbedIntrinsic] = intrinsic_strategy!=NULL))
    assert(secondderiv!=NULL);
}

extern void addTBM(int nr, isofct cov_tbm2, isofct cov_tbm3,
		   SimulationType tbm_method,
		   randommeasure spectral) {
  // must be called always AFTER addCov !!!!
  assert((nr>=0) && (nr<currentNrCov));
  CovList[nr].cov_tbm2=cov_tbm2;
  CovList[nr].implemented[TBM2] = (cov_tbm2!=NULL) ? IMPLEMENTED : 
    (CovList[nr].derivative!=NULL) ? NUM_APPROX : NOT_IMPLEMENTED;
  CovList[nr].cov_tbm3=cov_tbm3;
  CovList[nr].tbm_method = tbm_method;
  if (cov_tbm3!=NULL || cov_tbm2!=NULL || CovList[nr].derivative!=NULL || 
      tbm_method!=Nothing){
    assert(CovList[nr].derivative!=NULL && tbm_method!=Nothing);
  }
  CovList[nr].implemented[TBM3] = CovList[nr].cov_tbm3!=NULL;
  CovList[nr].spectral=spectral;
  CovList[nr].implemented[SpectralTBM] = spectral!=NULL;
}
	
extern void addOther(int nr, initmppfct add_mpp_scl, MPPRandom add_mpp_rnd, 
		     hyper_pp_fct hyper_pp,
		     generalSimuInit initother, generalSimuMethod other){
  assert((nr>=0) && (nr<currentNrCov));
  CovList[nr].add_mpp_scl=add_mpp_scl;  
  CovList[nr].add_mpp_rnd=add_mpp_rnd;  
  if (add_mpp_scl!=NULL || add_mpp_rnd!=NULL)
    assert(add_mpp_scl!=NULL && add_mpp_rnd!=NULL);
  CovList[nr].implemented[AdditiveMpp] = add_mpp_scl!=NULL;
  CovList[nr].hyperplane=hyper_pp;
  CovList[nr].implemented[Hyperplane] = hyper_pp!=NULL;
  CovList[nr].initother=initother;
  CovList[nr].other=other;
  if ((other!=NULL) || (initother!=NULL)) 
    assert((other!=NULL) && (initother!=NULL));
}

void PrintModelList()
{
  int i, last_method, m;
  char percent[]="%";
  char header[]="Circ cut intr TBM2 TBM3 spec drct nugg add hyp oth\n";
  char coded[4][2]={"-", "X", "+", "o"};
  char firstcolumn[20], empty[5][5]={"", " ", "  ", "   ", "    "};
  int empty_idx[(int) Nothing] = {1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 0};

  last_method = (int) Nothing;    
  if (currentNrCov==-1) {
    InitModelList(); 
    if (GENERAL_PRINTLEVEL>5) 
      PRINTF("List of covariance functions initiated.\n"); 
  }
  if (CovList==NULL) {PRINTF("There are no functions available!\n");} 
  else {
    sprintf(firstcolumn,"%s%ds ",percent, COVMAXCHAR);
    PRINTF("\n\n");
    PRINTF(firstcolumn,empty[0]); PRINTF("       List of models\n");
    PRINTF(firstcolumn,empty[0]); PRINTF("       ==============\n");
    PRINTF(firstcolumn,empty[0]); PRINTF("[See also PrintMethodList()]\n\n");
    PRINTF(firstcolumn,empty[0]);PRINTF(header,empty[0]);  
    for (i=0;i<currentNrCov;i++) {
      PRINTF(firstcolumn, CovList[i].name);
      for (m=(int) CircEmbed; m<(int) Nothing; m++)
	PRINTF("%3s%s", coded[(int) CovList[i].implemented[m]], 
	       empty[empty_idx[m]]);
      PRINTF("\n");
    }
    PRINTF(firstcolumn,empty[0]);PRINTF(header,empty[0]);  
    PRINTF("Legend:");
    PRINTF("'-': method not available\n");
    PRINTF("'X': method available for at least some parameters\n");
    PRINTF("'+': parts are evaluated only approximatively\n");
    PRINTF("'o': given method is ignored and an alternative one is used\n");
  }
}

void GetModelList(int* idx) {
  int i, j, m;
  if (currentNrCov==-1) {
    InitModelList(); 
    if (GENERAL_PRINTLEVEL>5) 
      PRINTF("List of covariance functions initiated.\n"); 
  }
  if (CovList==NULL) return;
  for (j=i=0; i<currentNrCov; i++)
    for (m=(int) CircEmbed; m<(int) Nothing; m++)
      idx[j++] = CovList[i].implemented[m];
  return;
}

void GetNaturalScaling(int *covnr, double *q, int *naturalscaling,
		       double *natscale, int *error)
{
  //  q: kappas only

  // values of naturalscaling:
  //#define NATSCALE_EXACT 1   
  //#define NATSCALE_APPROX 2
  //#define NATSCALE_MLE 3 /* check mleRF when changing !! */
  // +10 if numerical is allowed
  static int oldcovnr = -99;
  static double oldp[TOTAL_PARAM], p[TOTAL_PARAM];
  static double OldNatScale;
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
      double x,newx,yold,y,newy;
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
	double leftx;
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
	double rightx;
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
		int *lt, double *lambda, int *ll, int *error)
{
  unsigned long len;
  key_type *key;
  if (currentNrCov==-1) InitModelList(); assert(CovList!=NULL); 
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*error=ERRORREGISTER;goto ErrorHandling;}

  key = &(KEY[*keyNr]);
  DeleteKeyTrend(key);
  switch(*modus){
      case TREND_MEAN : //0
	key->mean = *lambda;
	assert(*lt==0 && *ll==1);
	break;
      case TREND_PARAM_FCT : case TREND_FCT: // 3, 1
	key->LinearTrend = (double *) malloc(len = sizeof(double) * *ll);
	memcpy(key->LinearTrend, lambda, len);
	if (*modus==1) break; 
	key->lLinTrend = *ll;
      case TREND_LINEAR : //2
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


void GetTrueDim(bool anisotropy, int timespacedim, double* param, 
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
  endfor = k = (int) RF_NAN;
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
  // currently unused, since AddTrend in RFsimu.cc is being programmed
  assert(false);

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


void GetTrend(int *keyNr, char **trend, double *lambda, double *x, int *error)
{
  // currently unused, since AddTrend in RFsimu.cc is being programmed
  assert(false);
  
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


void GetCoordinates(int *keyNr, double *x, int *error)
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

SEXP keynum(double* V, int n) {
  int i;
  SEXP dummy;
  PROTECT(dummy=allocVector(REALSXP, n));
  for (i=0; i<n; i++) REAL(dummy)[i] = V[i];
  UNPROTECT(1);
  return dummy;
}

SEXP GetExtKeyInfo(SEXP keynr, SEXP Ignoreactive) { // extended
#define ninfo 13
  char *infonames[ninfo] = 
    {"active", "grid", "anisotropy", "time.comp",
     "distrib.", "model",  "mean", "TBMline", "spatial.dim", 
     "timespacedim", "spatial.pts", "total.pts", "trend.mode", 
    };
#define nmodelinfo 5 // mit op
  char *modelinfo[nmodelinfo] = {"op", "name", "method", "param", "processed"};
  int knr, ni, i, actninfo, subi, k;
  bool ignore_active;
  SEXP info, namevec, model, submodel[MAXCOV], namemodelvec[2]; 
  key_type *key;

  knr = INTEGER(keynr)[0];
  ignore_active = LOGICAL(Ignoreactive)[0];
  if ((knr<0) || (knr>=MAXKEYS)) {
    return allocVector(VECSXP, 0);
  }
  key = &(KEY[knr]);

  actninfo = (key->active || ignore_active) ? ninfo : 1;
  PROTECT(info=allocVector(VECSXP, actninfo));
  PROTECT(namevec = allocVector(STRSXP, actninfo));
  for (i=0; i<actninfo; i++) SET_VECTOR_ELT(namevec, i, mkChar(infonames[i]));
  setAttrib(info, R_NamesSymbol, namevec);
  UNPROTECT(1);

  ni = 0;
  SET_VECTOR_ELT(info, ni++, ScalarLogical(key->active));
  if (!key->active && !ignore_active) {UNPROTECT(1); return info;}
  SET_VECTOR_ELT(info, ni++, ScalarLogical(key->grid)); 
  SET_VECTOR_ELT(info, ni++, ScalarLogical(key->anisotropy)); 
  SET_VECTOR_ELT(info, ni++, ScalarLogical(key->Time)); 
  SET_VECTOR_ELT(info, ni++, mkString(DISTRNAMES[key->distribution]));

  PROTECT(model = allocVector(VECSXP, key->ncov));
  for (i=0; i<=1; i++) {
    int j;
    PROTECT(namemodelvec[i] = allocVector(STRSXP, nmodelinfo - (i==0)));
    for (j=0,k=(i==0); k<nmodelinfo; k++, j++) 
      SET_VECTOR_ELT(namemodelvec[i], j, mkChar(modelinfo[k]));
  }
  for (i=0; i<key->ncov; i++) {
    PROTECT(submodel[i] = allocVector(VECSXP, nmodelinfo - (i==0)));
    setAttrib(submodel[i], R_NamesSymbol, namemodelvec[i>0]);

    subi = 0;
    if (i>0) SET_VECTOR_ELT(submodel[i], subi++,mkString(OP_SIGN[key->op[i-1]]));
    SET_VECTOR_ELT(submodel[i], subi++, mkString(CovList[key->covnr[i]].name));
    SET_VECTOR_ELT(submodel[i], subi++, mkString(METHODNAMES[key->method[i]]));
    SET_VECTOR_ELT(submodel[i], subi++, keynum(key->param[i], key->totalparam));
    SET_VECTOR_ELT(submodel[i], subi++, ScalarLogical(!key->left[i]));
    assert(subi==nmodelinfo - (i==0));
 
    SET_VECTOR_ELT(model, i, submodel[i]);
    UNPROTECT(1);
  }
  SET_VECTOR_ELT(info, ni++, model);
  UNPROTECT(3); // model + namemodelvec

  SET_VECTOR_ELT(info, ni++, ScalarReal(key->mean));
  SET_VECTOR_ELT(info, ni++, mkString(METHODNAMES[key->tbm_method]));
  SET_VECTOR_ELT(info, ni++, ScalarInteger(key->spatialdim));
  SET_VECTOR_ELT(info, ni++, ScalarInteger(key->timespacedim));
  SET_VECTOR_ELT(info, ni++, ScalarInteger((int) key->spatialtotalpoints));
  SET_VECTOR_ELT(info, ni++, ScalarInteger((int) key->totalpoints));
  SET_VECTOR_ELT(info, ni++, ScalarInteger(key->TrendModus));

  assert(ni==ninfo);
  UNPROTECT(1);
  return info;
}


void GetCornersOfGrid(key_type  *key, int Stimespacedim, int* start_aniso, 
		      double *param, double *sxx){
  // returning sequence (#=2^MAXDIM) vectors of dimension s->simutimespacedim
  double sx[ZWEIHOCHMAXDIM * MAXDIM];
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
      register double dummy;
      dummy = 0.0;
      for (g=0; g<key->timespacedim; g++) 
	dummy += param[start_aniso[k]+g] * sx[index+g];
      sxx[indexx + k] = dummy;
    }
  }
}

void GetCornersOfElement(key_type  *key, int Stimespacedim, int* start_aniso, 
		      double *param, double *sxx){
  double sx[ZWEIHOCHMAXDIM * MAXDIM];
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
      register double dummy;
      dummy = 0.0;
      for (g=0; g<key->timespacedim; g++) 
	dummy += param[start_aniso[k]+g] * sx[index+g];
      sxx[indexx + k] = dummy;
    }
  }
}

void GetRangeCornerDistances(key_type *key, double *sxx, int Stimespacedim,
			int Ssimuspatialdim, double *min, double *max) {
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

int eigenvalues(double *C, int dim, double *ev)
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

