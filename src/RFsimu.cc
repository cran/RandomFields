
/* 
 Authors
 Martin Schlather, Martin.Schlather@uni-bayreuth.de

 library for unconditional simulation of stationary and isotropic random fields

 Copyright (C) 2001 Martin Schlather, 

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


/*

  PRINTING LEVELS
  ===============
  error messages: 1
  forcation : 2
  minor tracing information : 3--5
  large debugging information: >10

*/


#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include <sys/timeb.h>
#include <assert.h>
#include <string.h>
#include "RFsimu.h"
#include "RFCovFcts.h"
#include "MPPFcts.h"

#include <unistd.h>

static cov_fct *CovList=NULL;
#ifdef RF_GSL
static gsl_rng *RF_RNG_STORED=NULL;
#endif

key_type KEY[MAXKEYS]; //static key_type KEY[MAXKEYS];
int currentNrCov=-1;
int GENERAL_STORING=true; 
/* true: intermediate results are stored: might be rather memory consuming,
         but the simulation can (depending on the method chosen) be much faster
	 when done the second time with exactly the same parameters
	 do not forget to call DeleteAllKeys when all the simulation are done.
   false: intermediate results are not stored if SimulateRF is called or
         stored only until DoSimulateRF is called.
	 minimum memory consumption, but might be rather slow if many
	 simulations are performed with exactly the same parameters
	 to be safe call DeleteAllKeys when all the simulations are done
	 Warning! init_foo may depend on GENERAL_STORING! This may cause problems
	 if GENERAL_STORING is changed between init_foo and do_foo from false to
	 true and if do_foo is called more than once! (SimulateRF is safe in
	 this case)

	 indifferent behaviour for the simulation methods if parameters are 
	 changed after each simulation.
	 In case, only MAXKEYS different parameter sets are used, but after each 
	 step the parameter set is changed, use different keynr for each 
	 parametere set, and STORING==true, to get fast simulation 
*/

int GENERAL_PRINTLEVEL=1;
/* 0:no messages; 
   1:error messages; 
   2: hints when algorithm is at a forcation; 
   >=3: more and more debugging information
*/

int GENERAL_NATURALSCALING=0; 
/* is called PracticalRange in R (see Chiles&Delfiner?!) 
   Note: RFparameters() allows only 0 and 1 as values!

   has an effect only if cov (and not only cov.local, e.g. Brownian Motion) is 
   defined
   0 : using the covariance function as defined
   1 : rescaling of cov fctn such that cov(1)~=0.05, if rescaling function 
       does not exist then failure, e.g. for Bessel model
   2 : exact or approximate value (e.g. expPLUScirc)
   3 : MLE (special needs taken into account, long memory covariance functions
       have too long tails, which are shortened (so threshold 0.05 is changed to
       higher values in this case))
  +10: if any of the above fails : numerical evaluation!
   else : using rescaling if rescaling function exists, otherwise without 
          rescaling
*/        

char ERRORSTRING_OK[MAXERRORSTRING],ERRORSTRING_WRONG[MAXERRORSTRING];

char METHODNAMES[][METHODMAXCHAR]={"circulant embedding",
				   "local CE (not implem.)",
                                   "TBM2 (2dim. turning bands)", 
				   "TBM3",
				   "spectral TBM (2dim)",
                                   "direct matrix decomposition",
				   "nugget",
                                   "add.MPP",
                                   "hyperplanes (not implem.)",
                                   "other methods(not implem.)",
				   "nothing", 
				   "max.MPP",
				   "forbidden"};

void GetParameterIndices(int *mean, int *variance, int *nugget, int *scale, 
			 int *kappa, int *lastkappa, int *invscale, int *sill) 
{
  *mean=MEAN; *variance=VARIANCE; *nugget=NUGGET; *scale=SCALE; *kappa=KAPPA;
  *lastkappa=LASTKAPPA; *invscale=INVSCALE; *sill=SILL;
}

void GetrfParameters(int *covmaxchar,int *methodmaxchar,
			    int *covnr, int *methodnr) {
  if (currentNrCov==-1) InitModelList();
  *covmaxchar=COVMAXCHAR;
  *methodmaxchar=METHODMAXCHAR;
  *covnr=currentNrCov;
  *methodnr=(int) Nothing;
}

void GetModelName(int *nr,char **name){
  if (currentNrCov==-1) InitModelList();
  if ((*nr<0) ||(*nr>=currentNrCov)) {strncpy(*name,"",COVMAXCHAR); return;}
  strncpy(*name,CovList[*nr].name,COVMAXCHAR);
}

void GetNrParameters(int *covnr,int* kappas) {
  if (currentNrCov==-1) InitModelList();
  if ((*covnr<0) ||(*covnr>=currentNrCov)) {*kappas=-1;}
  else *kappas = CovList[*covnr].kappas;
}

void GetModelNr(char **name, int *nr) 
{
  unsigned int ln;
  // == -1 if no matching function is found
  // == -2 if multiple matching functions are found, without one matching exactly
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)
  if (currentNrCov==-1) InitModelList();

  *nr=0;
  ln=strlen(*name);
  while ( (*nr<currentNrCov) && strncmp(*name,CovList[*nr].name,ln)) {
    (*nr)++;
  }
  if (*nr<currentNrCov) { 
    // a matching function is found. Are there other functions that match?
    int j; 
    bool exactmatching,multiplematching;
    exactmatching=(ln==strlen(CovList[*nr].name));
    multiplematching=false;
    j=*nr+1; // if two or more covariance functions have the same name the last 
    //          one is taken 
    while (j<currentNrCov) {
      while ( (j<currentNrCov) && strncmp(*name,CovList[j].name,ln)) {j++;}
      if (j<currentNrCov) {
	if (ln==strlen(CovList[j].name)) {*nr=j; exactmatching=true;} 
	else {multiplematching=true;}
      }
      j++;
    } 
    if (!exactmatching && multiplematching) {*nr=-2;}
  } else *nr=-1;
}

void GetMethodNr(char **name, int *nr) 
{
  unsigned int ln;
  // == -1 if no matching method is found
  // == -2 if multiple matching methods are found, without one matching exactly
  *nr=0;
  ln=strlen(*name);
  while ( (*nr<(int) Forbidden) && strncmp(*name,METHODNAMES[*nr],ln)) {
    (*nr)++;
  }
  if (*nr<(int) Forbidden) { 
    // a matching method is found. Are there other methods that match?
    int j; 
    bool exactmatching,multiplematching;
    exactmatching=(ln==strlen(METHODNAMES[*nr]));
    multiplematching=false;
    j=*nr+1; // if two or more methods have the same name the last one is
    //          taken; stupid here, but nice in GetCovNr
    while (j<(int) Forbidden) {
      while ( (j<(int) Forbidden) && strncmp(*name,METHODNAMES[j],ln)) {j++;}
      if (j<(int) Forbidden) {
	if (ln==strlen(METHODNAMES[j])) {*nr=j; exactmatching=true;} 
	else {multiplematching=true;}
      }
      j++;
    } 
    if (!exactmatching && multiplematching) {*nr=-2;}
  } else *nr=-1;
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
  PRINTF("%5d: %s cov=%d %d %d %d, tbm=%d %d %d, spec=%d,\n        abf=%d hyp=%d, oth=%d %d, f1=%d f2=%d f3=%d\n",
	 i,CovList[i].name,CovList[i].cov,CovList[i].cov_loc,
	 CovList[i].extension_factor,CovList[i].square_factor,
	 CovList[i].cov_tbm2,CovList[i].cov_tbm3,CovList[i].tbm_method,
	 CovList[i].spectral,CovList[i].add_mpp_scl,CovList[i].hyperplane,
	 CovList[i].initother,CovList[i].other,
	 CovList[i].first[0],CovList[i].first[1],CovList[i].first[2]);
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

void InternalIncludeModel(char *name,  
			   covfct cov,   
			   scalefct naturalscale,		   
			   covfct cov_loc, 
			   parameterfct extension_factor,
			   scalefct square_factor,
			   covfct cov_tbm2, covfct cov_tbm3, 
			   SimulationType tbm_method,
			   randommeasure spectral, 
			   MPPScales add_mpp_scl, 
			   MPPRandom add_mpp_rnd, 
			   generalSimuInit initother,
			   generalSimuMethod other,
			   int kappas,
			   SimulationType r1,SimulationType r2,SimulationType r3,
			   checkfct check
			   )
{  
  if (currentNrCov==-1) InitModelList();  assert(CovList!=NULL);
  if (currentNrCov>=MAXNRCOVFCTS) {PRINTF("Error. List full.\n");return;}
  strncpy(CovList[currentNrCov].name,name,COVMAXCHAR-1);
  if (strlen(name)>=COVMAXCHAR) {
    PRINTF("Warning! Covariance name is truncated to `%s'",
	   CovList[currentNrCov].name);
  }
  CovList[currentNrCov].cov=cov;
  CovList[currentNrCov].naturalscale=naturalscale;
  CovList[currentNrCov].cov_loc=cov_loc;
  CovList[currentNrCov].cov_tbm2=cov_tbm2;
  CovList[currentNrCov].cov_tbm3=cov_tbm3;
  CovList[currentNrCov].spectral=spectral;
  CovList[currentNrCov].add_mpp_scl=add_mpp_scl;  
  CovList[currentNrCov].add_mpp_rnd=add_mpp_rnd;  
  CovList[currentNrCov].hyperplane=NULL;
  CovList[currentNrCov].initother=initother;
  CovList[currentNrCov].other=other;
  CovList[currentNrCov].first[0]=r1;
  CovList[currentNrCov].first[1]=r2;
  CovList[currentNrCov].first[2]=r3;
  CovList[currentNrCov].other=other;
  CovList[currentNrCov].extension_factor=extension_factor;
  CovList[currentNrCov].square_factor=square_factor;
  CovList[currentNrCov].tbm_method      =tbm_method;
  CovList[currentNrCov].kappas = kappas;
  CovList[currentNrCov].check = check;
  currentNrCov++;
}

void IncludeModel(char *name, covfct cov, scalefct naturalscale,   
		   covfct cov_loc, parameterfct extension_factor,
		   scalefct square_factor,
		   covfct cov_tbm2, covfct cov_tbm3, SimulationType tbm_method,
		   randommeasure spectral, 
		   MPPScales add_mpp_scl, 
		   MPPRandom add_mpp_rnd, 
		   generalSimuInit initother, generalSimuMethod other,
		   int kappas,
		   SimulationType r1,SimulationType r2,SimulationType r3,
		   checkfct check
		   )
{
  InternalIncludeModel(name,cov,naturalscale,cov_loc,extension_factor,
			square_factor,
			cov_tbm2, cov_tbm3, tbm_method,
			spectral, 
			add_mpp_scl, add_mpp_rnd,
			initother,other, kappas,
			r1,r2,r3,check
			);
}


void IncludeModel(char *name, covfct cov, scalefct naturalscale,	
		   covfct cov_loc, parameterfct extension_factor, 
		   scalefct square_factor,
		   covfct cov_tbm2, covfct cov_tbm3, SimulationType tbm_method,
		   randommeasure spectral, 
		   int kappas,
		   SimulationType r1,SimulationType r2,SimulationType r3,
		   checkfct check
		   )
{
  InternalIncludeModel(name,cov,naturalscale,cov_loc,extension_factor,
			square_factor,
			cov_tbm2,cov_tbm3,tbm_method,spectral, 
			NULL,NULL,NULL,NULL,kappas,
			r1,r2,r3,check
			);
}

void IncludeModel(char *name, covfct cov,scalefct naturalscale,  
		   covfct cov_tbm2, covfct cov_tbm3, SimulationType tbm_method,
		   randommeasure spectral, 
		   int kappas, // how many additional parameters??
		   SimulationType r1,SimulationType r2,
		   SimulationType r3,//if a d dimen. field is to be generated
		   //                  which method is the prefered one?
		   checkfct check
		   )
{
  InternalIncludeModel(name,cov,naturalscale,NULL,NULL,NULL,
			cov_tbm2,cov_tbm3,tbm_method,spectral, 
			NULL,NULL,NULL,NULL,kappas,
			r1,r2,r3,check
			);
}



void InitModelList()
{
  int i,d;

#ifdef RF_GSL
  StoreSeed();
#endif
  for (i=0;i<MAXKEYS;i++) {
    KEY[i].active = false;
    KEY[i].S=NULL;
    KEY[i].destruct=NULL;
    KEY[i].SX=NULL;
    KEY[i].destructX=NULL;
    for (d=0; d<MAXDIM; d++) KEY[i].x[d]=NULL;
    KEY[i].cov = NULL;
    KEY[i].naturalscaling =-1;
  }
  if (CovList!=NULL) {
    PRINTF("List of covariance functions looks already initiated.\n"); 
    return;
  }
  assert(currentNrCov=-1);
  CovList = (cov_fct*) malloc(sizeof(cov_fct) * MAXNRCOVFCTS);
  currentNrCov=0;
  IncludeModel("bessel",Bessel,NULL,NULL,NULL,Nothing,spectralBessel,
		1,CircEmbed,SpectralTBM,CircEmbed,checkBessel);
  IncludeModel("cauchy",Cauchy,ScaleCauchy,TBM2Cauchy,TBM3Cauchy,CircEmbed,NULL,
		1,CircEmbed,TBM3,TBM3,checkCauchy);
  IncludeModel("cauchytbm",Cauchytbm,NULL,NULL,TBM3Cauchytbm,CircEmbed,NULL,
		3,CircEmbed,TBM3,TBM3,checkCauchytbm); 
  IncludeModel("circular",circular,Scalecircular, 
		NULL, //covfct cov_loc,    
		NULL, //     parameterfct extension_factor,   
		NULL, //     scalefct square_factor,
		NULL, //covfct cov_tbm2,    
		NULL, //covfct cov_tbm3,    
		Nothing, //   SimulationType tbm_method,
		NULL, //randommeasure spectral, 
		circular_init, //MPPScales,
		circularMpp, //MPPRandom, 
		NULL, //generalSimuInit initother,   
		NULL, //generalSimuMethod other,
		0,    //   int kappas,		      
		CircEmbed, //SimulationType r1,   
		CircEmbed, //SimulationType r2,   
		Forbidden, //SimulationType r3,
		NULL);     
  IncludeModel("cubic",cubic,Scalecubic,NULL,TBM3cubic,CircEmbed,NULL,
		0,CircEmbed,CircEmbed,CircEmbed,NULL); 
  IncludeModel("cone", 
		NULL, //covfct cov, 
		NULL, //scalefct naturalscale,   
		NULL, //covfct cov_loc,    
		NULL, //     parameterfct extension_factor,   
		NULL, //     scalefct square_factor,
		NULL, //covfct cov_tbm2,    
		NULL, //covfct cov_tbm3,    
		Nothing, //   SimulationType tbm_method,
		NULL , //randommeasure spectral, 
		cone_init, //MPPScale,
		cone, //MPPRandom, 
		NULL, //generalSimuInit initother,    
		NULL, //generalSimuMethod other,
		3,    //   int kappas,		      
		AdditiveMpp, //SimulationType r1,   
		AdditiveMpp, //SimulationType r2,   
		AdditiveMpp, //SimulationType r3,
		checkcone //   checkfct check
		);
  IncludeModel("exponential", exponential,Scaleexponential,NULL,
		TBM3exponential,CircEmbed,NULL,
		0,CircEmbed,CircEmbed,CircEmbed,NULL);
  // deleted as not of general interest
  //  IncludeModel("expPLUScirc",expPLUScirc,ScaleexpPLUScirc,
  //		NULL,NULL,Nothing,NULL,
  //		2,CircEmbed,CircEmbed,Forbidden,checkexpPLUScirc);  
  IncludeModel("gauss", 
		Gauss, //covfct cov, 
		ScaleGauss, //scalefct naturalscale,   
		NULL, //covfct cov_loc,    
		NULL, //     parameterfct extension_factor,   
		NULL, //     scalefct square_factor,
		NULL, //covfct cov_tbm2,    
		TBM3Gauss, //covfct cov_tbm3,    
		CircEmbed, //   SimulationType tbm_method,
		spectralGauss, //randommeasure spectral, 
		gaussmpp_init, //MPPScales,
		gaussmpp, //MPPRandom, 
		NULL, //generalSimuInit initother,    
		NULL, //generalSimuMethod other,
		0,    //   int kappas,		      
		CircEmbed, //SimulationType r1,   
		SpectralTBM, //SimulationType r2,   
		CircEmbed, //SimulationType r3,
		NULL //   checkfct check
		);     
  IncludeModel("gencauchy",generalisedCauchy, ScalegeneralisedCauchy,
		NULL,TBM3generalisedCauchy,CircEmbed,NULL,
		2,CircEmbed,TBM3,TBM3,checkgeneralisedCauchy);
  IncludeModel("gengneiting",genGneiting,NULL,NULL,
	       TBM3genGneiting, CircEmbed,NULL,
		2,CircEmbed,CircEmbed,CircEmbed,checkgenGneiting);
  IncludeModel("gneiting",Gneiting,ScaleGneiting,NULL,
		TBM3Gneiting,CircEmbed,NULL,
		0,CircEmbed,CircEmbed,CircEmbed,NULL);
  IncludeModel("gneitingdiff",Gneitingdiff,NULL,NULL,
		 TBM3Gneitingdiff, CircEmbed,NULL,
		2,CircEmbed,CircEmbed,CircEmbed,checkGneitingdiff);
  IncludeModel("holeeffect",holeeffect,Scaleholeeffect,NULL,
		TBM3holeeffect,CircEmbed,NULL,  
    		1,CircEmbed,CircEmbed,CircEmbed,checkholeeffect);
  IncludeModel("hyperbolic",hyperbolic,NULL,NULL,
		TBM3hyperbolic,CircEmbed,NULL,   
		3,CircEmbed,CircEmbed,CircEmbed,checkhyperbolic);
  IncludeModel("nugget",nugget,NULL,NULL,NULL,Nothing,NULL,
		0,Nugget,Nugget,Nugget,NULL);
  IncludeModel("penta",penta,Scalepenta,NULL,TBM3penta,CircEmbed,NULL,
		0,CircEmbed,CircEmbed,CircEmbed,NULL); 
  IncludeModel("power",power,Scalepower,TBM2power,
		TBM3power,CircEmbed,NULL,
		1,CircEmbed,CircEmbed,CircEmbed,checkpower);
  IncludeModel("qexponential", qexponential,Scaleqexponential,NULL,
		NULL,Nothing,NULL,
		1,CircEmbed,CircEmbed,CircEmbed,checkqexponential);
  IncludeModel("spherical",spherical,Scalespherical,
		NULL, //covfct cov_loc,    
		NULL, //     parameterfct extension_factor,   
		NULL, //     scalefct square_factor,
		TBM2spherical,//covfct cov_tbm2,
		TBM3spherical,//covfct cov_tbm3,
		CircEmbed,//   SimulationType tbm_method,
		NULL, //randommeasure spectral, 
		spherical_init, //MPPScales,
		sphericalMpp, //MPPRandom, 
		NULL, //generalSimuInit initother,    
		NULL, //generalSimuMethod other,
		0,    //   int kappas,		      
		CircEmbed, //SimulationType r1,   
		CircEmbed, //SimulationType r2,   
		CircEmbed, //SimulationType r3,
		NULL);     
  IncludeModel("stable",stable,Scalestable,NULL,TBM3stable,CircEmbed,NULL,
		1,CircEmbed,CircEmbed,CircEmbed,checkstable);
  IncludeModel("wave",wave,Scalewave,NULL,NULL,Nothing,spectralwave,
		0,CircEmbed,CircEmbed,CircEmbed,NULL);
  IncludeModel("whittlematern",WhittleMatern,ScaleWhittleMatern,NULL,
		TBM3WhittleMatern,CircEmbed,spectralWhittleMatern,
		1,CircEmbed,SpectralTBM,CircEmbed,checkWhittleMatern);

   // see /home/martin/article/C/RFCovBrownian.cc
}

void PrintModelList()
{
  int i;
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
	PRINTF(line,
	       coded[(CovList[i].cov!=NULL) && (CovList[i].cov_loc==NULL)],
	       coded[CovList[i].cov_loc!=NULL],
	       coded[CovList[i].cov_tbm2!=NULL],
	       coded[CovList[i].cov_tbm3!=NULL],
	       coded[CovList[i].spectral!=NULL], 
	       coded[(CovList[i].cov!=NULL) && (CovList[i].cov_loc==NULL)],
	       coded[CovList[i].add_mpp_scl!=NULL],
	       coded[CovList[i].hyperplane!=NULL],
	       coded[CovList[i].initother!=NULL]
	       );
      }
      /*
      if (i==18) {
	PRINTF(firstcolumn,empty);PRINTF(header,empty); 
	PRINTF("\n\
**** The following functions are part of current research work.       ****\n\ 
**** Please check with Tilmann Gneiting, tilmann@stat.washington.edu, ****\n\ 
**** before using them!                                               ****\n");
      }
      */
    }
    PRINTF(firstcolumn,empty);PRINTF(header,empty);   
  }
}


int InternalGetGridSize(Real *x[MAXDIM], int *dim, int lx[MAXDIM])
{
  int d;
  for (d=0; d<*dim; d++) {
    if ((x[d][XEND]<=x[d][XSTART]) || (x[d][XSTEP]<=0))
      return ERRORCOORDINATES;
    lx[d] = 1 + (int)((x[d][XEND]-x[d][XSTART])/x[d][XSTEP]); 
  }
  for (d=*dim; d<MAXDIM; d++)
    lx[d]=1;
  return NOERROR;
}

void GetGridSize(Real *x,Real *y,Real *z,int *dim,int *lx,int *ly,int *lz)
{ 
  // should now be consistent with `seq' of R
  // if not, please report!!!!!!!!!!!! 
  Real *xx[MAXDIM];
  int lxx[MAXDIM];
  xx[0]=x; xx[1]=y; xx[2]=z;
  if (InternalGetGridSize(xx,dim,lxx)) {
    *lx = *ly = *lz = 0;
  } else {
    *lx = lxx[0];  *ly = lxx[1];  *lz=lxx[2];
  }
}

void ErrorMessage(SimulationType m, int error) 
{
  char MS[50], EM[100+2*MAXERRORSTRING];
  switch(m) {
  case CircEmbed : strcpy(MS,"circulant embedding"); break;
  case CircEmbedLocal : strcpy(MS,"local circulant embedding"); break;
  case TBM2 : strcpy(MS,"2-dim. TBM"); break;
  case TBM3: strcpy(MS,"3-dim. TBM"); break;
  case SpectralTBM: strcpy(MS,"spectral TBM"); break;
  case Direct : strcpy(MS,"direct Gaussian (decomposition of cov. matrix)");
    break;
  case Nugget: strcpy(MS,"nugget effect modelling"); break;
  case AdditiveMpp: strcpy(MS,"additive MPP (random coins)"); break;
  case Hyperplane : strcpy(MS,"hyperplane tessellation"); break;
  case Special: strcpy(MS,"special procedure"); break;
  case Nothing: 
    if (error>0) strcpy(MS,"Error");
    else strcpy(MS,"Message"); 
    break;
  case MaxMpp: strcpy(MS,"max. MPP (Boolean functions)"); break;
  case Forbidden: strcpy(MS,"forbidden"); break;
  default : assert(false);
  }
  switch (error) {
  case USEOLDINIT : strcpy(EM,"Using stored initialization"); break;
  case NOERROR : strcpy(EM,"fine"); break;

  case ERRORNOTDEFINED :       
    strcpy(EM,"undefined for this model");break;
  case ERRORMETHODNOTALLOWED : 
    strcpy(EM,"not allowed for the specified parameters or points");break;
  case ERRORNOTPROGRAMMED :    
    strcpy(EM,"not programmed yet. Sorry"); break;
  case ERRORCOVNOTALLOWED :    
    strcpy(EM,"model not allowed for specified dimension");break;
  case ERRORFAILED: 
    strcpy(EM,"algorithm failed");break;
  case ERRORMEMORYALLOCATION: 
    strcpy(EM,"memory allocation error");break;
  case ERRORNOTINITIALIZED: 
    strcpy(EM,"not initialized or RFparameter()$Storage==FALSE");break;
  case ERRORKEYNROUTOFRANGE: 
    strcpy(EM,"wrong key number");break;
  case ERRORDECOMPOSITION:
    strcpy(EM,"matrix decomposition failed");break;
  case ERRORPRECISION: 
    strcpy(EM,"required precision not attained");break;
  case ERRORRESCALING:
    strcpy(EM,"rescaling not defined");break;
  case ERRORFOURIER: 
    strcpy(EM,"fft factorization failed");break;
  case ERRORCOVFAILED: 
    sprintf(EM,
	    "model and method only valid for %s. Got %s.\n",
	    ERRORSTRING_OK,ERRORSTRING_WRONG);
    break;
  case ERRORREGISTER: 
    sprintf(EM,"Register number out of [0,%d]",MAXKEYS-1);break;
  case ERRORCOORDINATES: 
    strcpy(EM,"coordinates are not given or invalid grid specification");break;
  case ERRORNEGATIVEVAR: 
    strcpy(EM,"Variance or nugget non-positive");break;
  case ERRORPARAMNUMBER: 
    strcpy(EM,"number of parameters not correct");break;
  case ERRORDIM: 
    sprintf(EM,"dimension specification not in [1,%d]",MAXDIM);break;
  case ERRORNEGATIVESCALE :  
    strcpy(EM,"scale parameter must be positive");break;
  case ERRORWAVING :
    strcpy(EM,"Rescaling not possible (waving)");break;
  case ERRORDUMMY : 
    strcpy(EM,"none"); break;
   case ERRORWRONGINIT: 
    strcpy(EM,"Wrong initialisation");break;
  case ERRORNN:
    strcpy(EM,"The number of points on the line is too large. Check your
parameters and make sure that none of the locations are given twice");break;
   //
   // extremes:
  case ERRORSILLNULL : 
    strcpy(EM,"Vanishing sill not allowed");break;     
    //    case : strcpy(EM,"");break;
   //
   // Poisson:
  case ERRORVARMEAN :
    strcpy(EM,"Variance differs from mean");break;
  default : assert(false);
  }
  if (m==Nothing) { PRINTF(" %s: %s.\n",MS,EM); }
  else  PRINTF("Method `%s': %s.\n",MS,EM);
}


void DeleteKey(int *keyNr)
{  
  int d;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {return;}
  if (GENERAL_PRINTLEVEL>=4) { PRINTF("deleting stored parameters...\n");}
  for (d=0; d<MAXDIM; d++) {
    if (KEY[*keyNr].x[d]!=NULL) {free(KEY[*keyNr].x[d]);KEY[*keyNr].x[d]=NULL;}
  }
  if (KEY[*keyNr].destruct!=NULL) KEY[*keyNr].destruct(&(KEY[*keyNr].S));
  KEY[*keyNr].destruct=NULL;
  if (KEY[*keyNr].destructX!=NULL) KEY[*keyNr].destructX(&(KEY[*keyNr].SX));
  KEY[*keyNr].destructX=NULL;
  KEY[*keyNr].active=false;
}

void DeleteAllKeys() {
  int i;
  for (i=0;i<MAXKEYS;i++) {DeleteKey(&i);}
}


void printKEY(int keyNr) {
  int i,actparam;  
  if ((keyNr<0) || (keyNr>=MAXKEYS)) {PRINTF("\nkeynr out of range!\n");return;}
  PRINTF("Nr=%d, grid=%d, active=%d, nug.incl=%d, method=%d tbm_method=%d\n",
	 keyNr,
	 KEY[keyNr].grid, KEY[keyNr].active, KEY[keyNr].nuggetincluded,
	 KEY[keyNr].method, KEY[keyNr].tbm_method);
  PRINTF("^x=%d ^y=%d ^z=%d ^S=%d ^destruct=%d \n", 
	 KEY[keyNr].x[0], KEY[keyNr].x[1], KEY[keyNr].x[2],
	 KEY[keyNr].S,KEY[keyNr].destruct);
  PRINTF("mean=%f, variance=%f, nugget=%f, scale=%f, sill=%f,  invscale=%f",
	 KEY[keyNr].param[MEAN],
	 KEY[keyNr].param[VARIANCE],
	 KEY[keyNr].param[NUGGET],
	 KEY[keyNr].param[SCALE],
	 KEY[keyNr].param[SILL],
	 KEY[keyNr].param[INVSCALE]);
  if (KEY[keyNr].cov) {
    PRINTF(", #kappa=%d",KEY[keyNr].cov->kappas);
    actparam = KAPPA+KEY[keyNr].cov->kappas;
    for (i=KAPPA;i<actparam;i++) { 
      PRINTF(", kappa%d=%f",i-KAPPA,KEY[keyNr].param[i]);
    }
  } 
  PRINTF("\ncovnr=%d, lx=%d, dim=%d, totallength=%d\n",
	 KEY[keyNr].covnr,KEY[keyNr].lx,KEY[keyNr].dim,KEY[keyNr].totallength);
}

void printAllKeys() {
  int i;
  for (i=0;i<MAXKEYS;i++) {printKEY(i);PRINTF("\n");}
}


void GetNaturalScaling(int *covnr, Real *p, int *actparam, int *naturalscaling,
		       Real *natscale, int *error)
{
  // values of naturalscaling:
  //#define NATSCALE_EXACT 1   
  //#define NATSCALE_APPROX 2
  //#define NATSCALE_MLE 3 /* check mleRF when changing !! */
  // +10 if numerical is allowed
  static int oldcovnr = -99;
  static Real oldp[TOTAL_PARAM];
  static Real OldNatScale;
  int TEN;
  bool numeric;  
  *error = 0;
  TEN = 10;
  if (*naturalscaling) {	 
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
      
      //printf("numeric\n");
      // already calculated ?
      parami=KAPPA; // do not compare mean,variance, etc.
      if (oldcovnr==*covnr) {
	for (; parami<*actparam; parami++) {
	  if (oldp[parami]!=p[parami]) {
	    break;
	  }
	}
	if (parami==*actparam) {
	  *natscale=OldNatScale; 
	  return;
	}
      }   

      // the other parameters need not to be copied as checked that 
      // they are identical to previous ones
      for (;parami<*actparam;parami++) {oldp[parami]=p[parami];}
      oldp[VARIANCE] = oldp[SCALE] = oldp[INVSCALE] = 1.0; 
      oldp[NUGGET] = oldp[MEAN] = 0.0;
      oldcovnr = -99; /* if error occurs, the next call will realise that 
			 the previous result cannot be used */

      /* **************************************
	 Now, find a quick and good solution for NewInvScale --
      */
      wave  = 0;
      x = 1.0; 
      if ( (yold=cov(x,oldp)) > 0.05) {
	Real leftx;
	x *= 2.0;
	while ( (y=cov(x,oldp)) > 0.05) {  
	  if (yold<y){ wave++;  if (wave>10) {*error=ERRORWAVING; return;} }
	  yold = y;
	  x *= 2.0;
	  if (x>1E30) {*error=ERRORRESCALING; return;} // or use a separate ERROR
	} 
	leftx = x * 0.5;
	for (i=0; i<3 /* good choice?? */ ;i++) {          
	  if (y==yold) {*error=ERRORWAVING; return;} // should never appear
	  newx = x + (x-leftx)/(y-yold)*(0.05-y);
	  if ( (newy=cov(newx,oldp)) >0.05) {
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
	while ( (y=cov(x,oldp)) < 0.05) {  
	  if (yold>y){ wave++;  if (wave>10) {*error=ERRORWAVING; return;} }
	  yold = y;
	  x *= 0.5;
	  if (x<1E-30) {*error=ERRORRESCALING; return;} //or use a separate ERROR
	}    
	rightx = x * 2.0;
	for (i=0; i<3 /* good choice?? */ ;i++) {          
	  if (y==yold) {*error=ERRORWAVING; return;} // should never appear
	  newx = x + (x-rightx)/(y-yold)*(0.05-y);
	  if ( (newy=cov(newx,oldp)) <0.05) {
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
     } else *error=ERRORRESCALING; 
  } else {
    *natscale = 1.0;
  }
}

int CheckAndRescale(int naturalscaling, int covnr, Real *p, int np, int dim, 
		    Real *newscale) 
{
  int actparam,error;
  Real natscale;

  if (currentNrCov==-1) InitModelList();assert(CovList!=NULL); 
  if ((covnr>=currentNrCov) || (covnr<0)) return ERRORNOTDEFINED;
  if ((dim<1) || (dim>MAXDIM)) return ERRORDIM;
  if (CovList[covnr].first[dim-1]==Forbidden ) return ERRORCOVNOTALLOWED;
  actparam = KAPPA + CovList[covnr].kappas;
  if (np!=actparam) {
    return ERRORPARAMNUMBER;}
  if (p[SCALE] <= 0.0) return ERRORNEGATIVESCALE;
  if ((p[VARIANCE]<0.0) || (p[NUGGET]<0.0)) return ERRORNEGATIVEVAR;

  if (CovList[covnr].check!=NULL) {    
    key_type key; int error,i;
    key.method = Nothing;
    key.dim = dim; 
    for (i=0; i<actparam; i++) {key.param[i]=p[i];}
    if (error=(CovList[covnr].check(&key))) return error;
  }
  GetNaturalScaling(&covnr,p,&actparam,&naturalscaling,&natscale,&error);
  if (error) return  error;
  *newscale = p[SCALE] * natscale;
  return 0;
}

int CheckCovariance(int covnr,int naturalscaling, Real *p, int np, Real *param, 
		    int dim, covfct *cov, bool *variogram)
{
  // called by Covariance and Variogram, only
  int endfor,i,error;
  Real newscale;
  // must be the very first !!
  if (error=(CheckAndRescale(naturalscaling,covnr,p,np,dim,&newscale))) 
    return error; 
  if ((*cov=CovList[covnr].cov)==NULL) return ERRORNOTDEFINED;
  endfor = KAPPA + CovList[covnr].kappas;
  for (i=0; i<endfor; i++) {param[i]=p[i];}
  param[SILL] = 1.0;
  *variogram = ( (*cov)(0,param)==0);   
  param[SILL] = param[NUGGET]+param[VARIANCE]; 
  param[SCALE] = newscale;
  param[INVSCALE] =  1.0 / param[SCALE];
  return 0;
}
 
// slow if only a single value is calculated -- due to lots of checks beforehand
// define an additional parameter CHECK? or an UnCheckedCovariance function?
void CovarianceNatSc(Real *x,int *lx,int *covnr,Real *p, int *np, int *dim,
		     Real *result, int *naturalscaling)
{
  int error,i;
  bool variogram;
  Real param[TOTAL_PARAM];
  covfct cov;
  if (error=(CheckCovariance(*covnr,*naturalscaling,p,*np,param,
			    *dim,&cov,&variogram))) {
    if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,error);
    for (i=0;i<*lx;i++) { result[i] = RF_NAN; }
    return;	
  }
  if (variogram) {
    for (i=0;i<*lx;i++) {  result[i] = RF_NAN;  }
  } else {
    for (i=0;i<*lx;i++) {  result[i] = cov(x[i],param);  }
  }
}

void Covariance(Real *x,int *lx,int *covnr,Real *p, int *np, int *dim,
		Real *result)
{
  CovarianceNatSc(x,lx,covnr,p,np,dim,result,&GENERAL_NATURALSCALING);
}
  
void VariogramNatSc(Real *x,int *lx,int *covnr,Real *p,int *np, int *dim,
		    Real *result, int *naturalscaling)
{
  int error,i;
  bool variogram;
  Real param[TOTAL_PARAM];
  covfct cov;
  // only if variogram corresponds to covariance matrix !!!
  if (error=(CheckCovariance(*covnr,*naturalscaling,p,*np,param,
			    *dim,&cov,&variogram)))  {
     if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,error);    
    for (i=0;i<*lx;i++) { result[i] = RF_NAN;}
    return;
  }
  if (variogram) {
    for (i=0;i<*lx;i++) {result[i] = cov(x[i],param);}
  }  else {
    for (i=0;i<*lx;i++) { result[i] = param[SILL]-cov(x[i],param);}
  }
}
void Variogram(Real *x,int *lx,int *covnr,Real *p,int *np, int *dim,
	       Real *result)
{
  VariogramNatSc(x,lx,covnr,p,np,dim,result,&GENERAL_NATURALSCALING);
}
  
void CovarianceMatrixNatSc(Real *dist,int *lx,int *covnr,Real *p,int *np, 
			   int *dim, Real *result, int *naturalscaling)
{
  int i,j,zaehler,index,jindex,error;
   bool variogram;
   Real param[TOTAL_PARAM];
   covfct cov;
    
   if (error=(CheckCovariance(*covnr,*naturalscaling,p,*np,param,
			    *dim,&cov,&variogram))) {
     if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,error); 
     for(i=0,index=0;i<*lx;i++) {  
       for(j=0;j<*lx;j++) { result[index] = RF_NAN; index++; }
     }
     return;
   }
 
   long lxq, lxqhalf;
   lxq = *lx * *lx;
   lxqhalf = lxq/2;

  index=zaehler=0;
  for (i=0;i<*lx;i++) { 
    for(j=0,jindex=i; j<i; j++,jindex+=*lx) {
      result[index+j] = result[jindex];
      assert(index+j<lxq);
      assert(jindex<lxq);
    }
    result[index + i] = param[SILL];
    if (variogram) {
      for(j=index+i+1,index+=*lx; j<index; j++){ result[j] = RF_NAN;}
     } else {
      for(j=index+i+1,index+=*lx; j<index; j++){
	assert(zaehler<lxqhalf); ////////////////
	assert(j< lxq); ////////////////
	result[j] = cov(dist[zaehler++],param);
      }
    }
  } 
}
void CovarianceMatrix(Real *dist,int *lx,int *covnr,Real *p,int *np, int *dim,
		      Real *result)
{
  CovarianceMatrixNatSc(dist,lx,covnr,p,np,dim,result,&GENERAL_NATURALSCALING);
}

void VariogramMatrix(Real *dist,int *lx,int *covnr,Real *p,int *np, Real *result)
{
}
  
void SetParam(int *action, int *storing,int *printlevel,int *naturalscaling) {
  switch (*action) {
  case 0 :
    GENERAL_STORING= (bool) *storing;
    GENERAL_PRINTLEVEL=*printlevel;
    GENERAL_NATURALSCALING=*naturalscaling;
    break;
  case 1 :
    *storing =  GENERAL_STORING;
    *printlevel = GENERAL_PRINTLEVEL;
    *naturalscaling = GENERAL_NATURALSCALING;
   if (GetNotPrint) break;
 case 2 : 
    PRINTF("\nGeneral Parameters\n==================\nstoring=%d\nprint level=%d\nnatural scaling=%d\n",
	   GENERAL_STORING,GENERAL_PRINTLEVEL,GENERAL_NATURALSCALING);
    break;
  default : PRINTF(" unknown action\n"); 
  }
}


#define DIFF_PARAMLIST 1
#define DIFF_COORDINATES 2
#define DIFF_OTHERS 3
#define DIFF_NOTACTIVE 4
#define DIFF_COVNR 5
#define DIFF_GRID 6
#define DIFF_METHOD 7
void InitSimulateRF(Real *x, Real *y, Real *z, int *dim, int *lx, int *grid, 
	            int *covnr, Real *ParamList, int *nParam,
		    int *method, int *distr,
		    int *keyNr , int *error)
{ // grid  : boolean
  // keyNr : label where intermediate results are to be stored;
  //         can be chosen freely between 0 and MAXKEYS-1
  bool user_defined,finished;
  int actparam,i,result,act_number,d;
  long totalBytes;
  Real newscale, *coord[MAXDIM];
  
  coord[0] = x;
  coord[1] = y;
  coord[2] = z; 


  // preference lists, distinguished by grid==true/false and dimension
  // lists must end with Nothing!
  SimulationType pg1[]={CircEmbed, CircEmbedLocal, Direct, AdditiveMpp, Nothing}; 
  SimulationType pg2[]={CircEmbed, CircEmbedLocal, TBM2, SpectralTBM, 
			 AdditiveMpp, TBM3, Direct, Nothing};
  SimulationType pg3[]={CircEmbed, CircEmbedLocal, TBM3, AdditiveMpp, 
			Direct, Nothing};
  SimulationType png1[]={Direct, AdditiveMpp, Nothing};
  SimulationType png2[]={TBM2, SpectralTBM, AdditiveMpp, TBM3,  Direct, 
			 Nothing};
  SimulationType png3[]={TBM3, AdditiveMpp, Direct, Nothing};
  SimulationType *preference_list,first_method;

 
  // check parameters
  // CheckAndRescale must be at the beginning !!

  if (*error=(CheckAndRescale(GENERAL_NATURALSCALING,*covnr,ParamList,*nParam,
			     *dim,&newscale))) {
    goto FatalError;
  }

  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*error=ERRORREGISTER;goto FatalError;}
  if (coord[0]==NULL) {*error=ERRORCOORDINATES;goto FatalError;}

  // bytes to be stored for each vector of coordinates:
  totalBytes=sizeof(Real) * *lx;  
  actparam = KAPPA+CovList[*covnr].kappas;

  result=0;
  if (GENERAL_PRINTLEVEL>=3) {PRINTF("comparing with stored parameters...\n");}

  if (KEY[*keyNr].active) {// if not active, all the parameters have to be 
    //                        stored again
 
    KEY[*keyNr].active = false;
    assert(KEY[*keyNr].x[0]!=NULL);
    
    // if lx, dim, or grid is different from previous run, there is a severe 
    // difference; thus everything has to be stored from scratch; i.e., 
    // act as if active==false
    if ((KEY[*keyNr].lx==*lx) && 
	(KEY[*keyNr].dim==*dim) && 
	(KEY[*keyNr].grid==(bool)*grid) && 
	(KEY[*keyNr].distribution==*distr)) {

      if (KEY[*keyNr].covnr!=*covnr) {
	result=DIFF_COVNR; 
	KEY[*keyNr].covnr=*covnr;
      }
      
     { 
	Real dummy;
	dummy = ParamList[SCALE];
	ParamList[SCALE] = newscale;
	
	for (i=0;i<actparam;i++) {
	  if (KEY[*keyNr].param[i]!=ParamList[i]) {result=DIFF_PARAMLIST; break;}
	}
	if (result) {
	  for (i=0;i<actparam;i++) { KEY[*keyNr].param[i]=ParamList[i]; }
	}
	ParamList[SCALE]=dummy;
      }

      // check coordinates
      for (d=0; d<*dim; d++) {
	if ((KEY[*keyNr].x[d]!=NULL) && 
	    (memcmp(coord[d],KEY[*keyNr].x[d],totalBytes))) {
	  result=DIFF_COORDINATES;
	  free(KEY[*keyNr].x[d]); KEY[*keyNr].x[d]=NULL;
	  if ((KEY[*keyNr].x[d]=(Real*)malloc(totalBytes))==NULL){
	    *error=ERRORMEMORYALLOCATION;goto FatalError;
	  }
	  memcpy(KEY[*keyNr].x[d],coord[d],totalBytes);
	}
      }
      // should be excluded as here it is asserted that the dimensions
      // *dim and KEY[*keyNr].dim equal
      for (d=*dim; d<MAXDIM; d++) {assert(KEY[*keyNr].x[d]==NULL);} 
      goto ende;
    } else {result=DIFF_OTHERS; } // ...(KEY[*keyNr].grid==(bool)*grid))
  } else {
    result=DIFF_NOTACTIVE;
  }

  // save all the parameters if key has not been active

  KEY[*keyNr].lx=*lx;
  KEY[*keyNr].dim=*dim;
  KEY[*keyNr].grid=(bool) *grid;
  KEY[*keyNr].covnr=*covnr; 
  KEY[*keyNr].distribution=*distr;

  for (i=0;i<actparam;i++) { KEY[*keyNr].param[i]=ParamList[i]; }
  KEY[*keyNr].param[SCALE] = newscale;

   for (d=0; d<MAXDIM; d++) {
    if (KEY[*keyNr].x[d]!=NULL) {
      free(KEY[*keyNr].x[d]);
      KEY[*keyNr].x[d]=NULL;
    }
  }

  for (d=0; d<KEY[*keyNr].dim; d++) {
    if ((KEY[*keyNr].x[d] = (Real*) malloc(totalBytes))==NULL) {
      *error=ERRORMEMORYALLOCATION; goto FatalError;
    }
    memcpy(KEY[*keyNr].x[d],coord[d],totalBytes);
  }

   if (GENERAL_PRINTLEVEL>=5)  PRINTF("before ");
 ende:
  if (GENERAL_PRINTLEVEL>=5)  PRINTF("label `ende'\n");
  if (user_defined = (*method>=0)) {
    if (*method >= (int) Nothing) {*error=ERRORNOTDEFINED; goto FatalError;}
    if ((KEY[*keyNr].method!=(SimulationType)*method) && 
	((KEY[*keyNr].method!=Nugget) || (ParamList[VARIANCE]!=0))) {	
      result=DIFF_METHOD; 
      if (ParamList[VARIANCE]==0) {
	KEY[*keyNr].method=Nugget;
      } else {
	KEY[*keyNr].method=(SimulationType) *method;
      }
    }
  } else {
    if (result) {
      assert(CovList[*covnr].first[KEY[*keyNr].dim-1]<Nothing);
      KEY[*keyNr].method=CovList[*covnr].first[KEY[*keyNr].dim-1];
    }	
  }

   if (!result) { // detected that all relevant parameters agree, 
    //              no initialization is necessary
     KEY[*keyNr].active = true;
     *error=USEOLDINIT; 
     if (GENERAL_PRINTLEVEL>=2) ErrorMessage(Nothing,*error);
     return; 
  }

 
  if (GENERAL_PRINTLEVEL>=5) { 
    printf("disagreement at %d\n",result);
    printKEY(*keyNr);
  }

  // Delete stored, intermediate result if there is any
  if (KEY[*keyNr].destruct!=NULL) {
    KEY[*keyNr].destruct(&(KEY[*keyNr].S));
    KEY[*keyNr].destruct = NULL;
  }
  if (KEY[*keyNr].destructX!=NULL) {
    KEY[*keyNr].destructX(&(KEY[*keyNr].SX));
    KEY[*keyNr].destructX = NULL;
  }
  
  // minor settings, usually overtaken from previous run, if all the 
  // other parameters are identical
  KEY[*keyNr].param[INVSCALE] = 1.0/KEY[*keyNr].param[SCALE];
  KEY[*keyNr].param[SILL] = 
    KEY[*keyNr].param[VARIANCE] + KEY[*keyNr].param[NUGGET];
  if (TBM_METHOD==Nothing) {KEY[*keyNr].tbm_method=CovList[*covnr].tbm_method;}
  else {KEY[*keyNr].tbm_method=TBM_METHOD;}

  if (*grid) { 
    if ((*lx!=3) || 
	(InternalGetGridSize(coord,&KEY[*keyNr].dim,KEY[*keyNr].length))) {
      *error=ERRORCOORDINATES; goto FatalError;
    }
    KEY[*keyNr].totallength = 1;
    for (d=0; d<*dim; d++) {
      KEY[*keyNr].totallength *= KEY[*keyNr].length[d];
    }
    switch (KEY[*keyNr].dim) {
    case 1 : preference_list = pg1; break;
    case 2 : preference_list = pg2; break;
    case 3 : preference_list = pg3; break;
    default : assert(false);
    }
  } else {
    KEY[*keyNr].totallength=KEY[*keyNr].length[0]=*lx; 
    for (d=1; d<*dim; d++) KEY[*keyNr].length[d]=0; 
    switch (KEY[*keyNr].dim) {
    case 1 : preference_list = png1; break;
    case 2 : preference_list = png2; break;
    case 3 : preference_list = png3; break;
    default : assert(false);
    }
  }

  KEY[*keyNr].cov=&CovList[KEY[*keyNr].covnr]; // pointer to the complete struct, 
  //                                              not to `cov' within struct

  // (hidden) pure nugget effect ??
  if (KEY[*keyNr].param[VARIANCE]==0) { KEY[*keyNr].method=Nugget;  } 


  if (GENERAL_PRINTLEVEL>=5) { printKEY(*keyNr);}


  first_method=KEY[*keyNr].method;
  finished=false; 
  while (!finished) {
   //  PrintCovDetailed(KEY[*keyNr].covnr); 
    if (GENERAL_PRINTLEVEL>=5) {ErrorMessage(KEY[*keyNr].method,ERRORDUMMY);}

    assert((KEY[*keyNr].destruct==NULL) && (KEY[*keyNr].S==NULL));

    if ((KEY[*keyNr].cov->check==NULL) ||
	((*error=KEY[*keyNr].cov->check(&KEY[*keyNr]))==0)) {
      switch (KEY[*keyNr].method) {
      case CircEmbed : 
	if (KEY[*keyNr].cov->cov==NULL) {*error=ERRORNOTDEFINED;} 
	else  {
	  *error=init_circulantembedding(&KEY[*keyNr]);
	  KEY[*keyNr].nuggetincluded = true;
	}
	break;
	
      case CircEmbedLocal : // not programmed yet
	if (KEY[*keyNr].cov->cov_loc==NULL) {*error=ERRORNOTDEFINED;} 
	else {
	  *error=ERRORNOTPROGRAMMED;
	  // *error = init_circ_embed_local(&KEY[*keyNr]);
	  // see /home/martin/article/C/RFce_local.cc
	  KEY[*keyNr].nuggetincluded = true; 
	}
	break;
	
      case TBM2 :
	if (KEY[*keyNr].cov->cov_tbm2==NULL) {*error=ERRORNOTDEFINED;} 
	else if (KEY[*keyNr].dim==2) { 
	  *error=init_turningbands(&KEY[*keyNr]);
	  KEY[*keyNr].nuggetincluded = false;
	}
	else if (KEY[*keyNr].dim==1) {
	  *error=ERRORNOTPROGRAMMED;
	}
	else { *error=ERRORMETHODNOTALLOWED; }
	break;
	
      case TBM3 :
	if (KEY[*keyNr].cov->cov_tbm3==NULL) {*error=ERRORNOTDEFINED;} 
	else if ((KEY[*keyNr].dim==2) || (KEY[*keyNr].dim==3)) { 
	  *error=init_turningbands(&KEY[*keyNr]);
	  KEY[*keyNr].nuggetincluded = false;
	}
	else if (KEY[*keyNr].dim==1) {
	  *error=ERRORNOTPROGRAMMED; 
	}
	else { *error=ERRORMETHODNOTALLOWED; }
	break;
	
      case SpectralTBM :
	if (KEY[*keyNr].cov->spectral==NULL) {*error=ERRORNOTDEFINED;}     
	else *error=init_simulatespectral(&KEY[*keyNr]);
	KEY[*keyNr].nuggetincluded = false;
	break;
	
      case Direct :
	if (KEY[*keyNr].cov->cov==NULL) {*error=ERRORNOTDEFINED;} 
	else *error=init_directGauss(&KEY[*keyNr]);
	KEY[*keyNr].nuggetincluded = true;
	break;
	
      case Nugget : 
	KEY[*keyNr].nuggetincluded = false;
	if (KEY[*keyNr].param[VARIANCE]==0.0) {*error=0;} 
	else {*error=ERRORMETHODNOTALLOWED;}
	break;
	
      case AdditiveMpp :// not programmed yet
	if (KEY[*keyNr].cov->add_mpp_scl==NULL) {*error=ERRORNOTDEFINED;} 
	else { 	  
	  init_mpp(&KEY[*keyNr]);
	  KEY[*keyNr].nuggetincluded = false;
	}
	break;
	
      case Hyperplane : // not programmed yet
	if (KEY[*keyNr].cov->hyperplane==NULL) {*error=ERRORNOTDEFINED;} 
	else {
	  // Abhaengigkeiten:
	  //    Funktion, die "stationaeren" punktprozess of [0,2 pi] x RR_+ liefert
	  //    dichte der punkte [KAPPA]
	  //    Wahl der Verteilungsfunktion [KAPPAX]
	  //                   ((0<alpha<2:heavy tailed), 2:z.B.uniform, 
	  //                   3:exponentiell, 4:gauss, etc.
	  //    Weiterer Parameter der Verteilungsfunktion [KAPPAXX] 
	  //                   (stable: schiefe; uniform[0,KAPPAXX], ...)
	  *error=ERRORNOTPROGRAMMED;
	  KEY[*keyNr].nuggetincluded = false;
	  // change preference list once it is programmed !
	}
	break;
	
      case Special :
	if (KEY[*keyNr].cov->other==NULL) { *error=ERRORNOTDEFINED; }
	else {
	  *error=(KEY[*keyNr].cov->initother)(&KEY[*keyNr]);      
	  KEY[*keyNr].nuggetincluded = false; // !!
	  // change preference list once it is programmed !
	}
	break;
	
      default:
	if (KEY[*keyNr].method==MaxMpp) ErrorMessage(MaxMpp,ERRORWRONGINIT);
	else assert(false);
	break;
      }
    } // if check covariance is ok.
    finished= user_defined || !(*error);  

    // the is a first_method given by programmer when defining a cov model;
    // this method is tried first; afterwards the standard list of methods is 
    // checked; however the first_method should not be reconsidered!
    if (!finished) {
      if (GENERAL_PRINTLEVEL>1) ErrorMessage(KEY[*keyNr].method,*error);
      if (first_method==KEY[*keyNr].method) { 
	act_number=0; // now start with the list; nonetheless, first_method 
	//               could have been the method nr==0, this is corrected 
	//               in the next step
      } else {
	act_number++; // we are in the middle of the standard list, just try 
	//               the next one
      }
      KEY[*keyNr].method=preference_list[act_number];
      if (KEY[*keyNr].method==first_method) {
	// the new(next or first) method out of the standard list might be 
	// the first_method
	act_number++;
	KEY[*keyNr].method=preference_list[act_number];
      }
      // end of the list?
      finished=finished || (KEY[*keyNr].method==Nothing);
    } //if !finished
  } // while !finished

  if (*error) { 
    if (user_defined) {
      if (GENERAL_PRINTLEVEL>0) ErrorMessage(KEY[*keyNr].method,*error);
    } else {
      *error=ERRORFAILED; 
      if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,*error);      
    }
  } else 
    if (!user_defined && GENERAL_PRINTLEVEL>1) 
      ErrorMessage(KEY[*keyNr].method,0);
  KEY[*keyNr].active = !*error;
  return;
  
  FatalError:
  //PRINTF("** Fatal Error occured **\n");
   if (GENERAL_PRINTLEVEL>0) ErrorMessage(Nothing,*error); 
  DeleteKey(keyNr);  
}

void GetKeyInfo(int *keyNr,int *total, int *lengths, int *dim, int *grid,
		int *distr)
{
  int d;
  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {*total=-1;} 
  else if (!KEY[*keyNr].active) {*total=-2;} 
  else {
    *total=KEY[*keyNr].totallength;
    *dim = KEY[*keyNr].dim;
    for (d=0; d<MAXDIM; d++) lengths[d]=KEY[*keyNr].length[d];
    *grid = (int) KEY[*keyNr].grid;
    *distr = (int) KEY[*keyNr].distribution;
  }
}
    
void RNGtest() {
  unsigned long  z0,z1;
  double UU;
#ifndef RF_GSL
  GetRNGstate();
#endif
  for (z0=z1=0; z1<100;) { 
    z0++; if (z0==1000000000) {z0=0;z1++; PRINTF("* z1=%d\n",z1);}
    UU=UNIFORM_RANDOM;
    if (UU==0.0) PRINTF("NULL: %d %d\n",z1,z0); else
      if (UU==0.5) PRINTF("0.5: %d %d\n",z1,z0); else
	if ((UU>0.5) && (UU<=0.500000001)) PRINTF("   ~0.5: %.15f\n",UU);
  }
#ifndef RF_GSL
  PutRNGstate();
#endif
  return;
}

void DoPoissonRF(int *keyNr, Real *res, int *error)
{
}

void DoSimulateRF(int *keyNr, Real *res, int *error)
{
  Real mean;
  long nx,endfor;

  *error=0; 

  if ((*keyNr<0) || (*keyNr>=MAXKEYS)) {
    *error=ERRORKEYNROUTOFRANGE; goto ErrorHandling;
  }
  if (!KEY[*keyNr].active) {*error=ERRORNOTINITIALIZED; goto ErrorHandling;}
  //  if (KEY[*keyNr].distribution!=DISTR_GAUSS) {
  //    *error=ERRORWRONGINIT; 
  //    goto ErrorHandling;
  //  }

#ifndef RF_GSL
     GetRNGstate();
#endif
  switch (KEY[*keyNr].method) {
  case CircEmbed : 
    do_circulantembedding(&KEY[*keyNr],res END_WITH_RANDOM); break; 
  case CircEmbedLocal : 
    do_circ_embed_local(&KEY[*keyNr],res END_WITH_RANDOM); break; 
  case TBM2 : 
    do_turningbands(&KEY[*keyNr],res END_WITH_RANDOM); break; 
  case TBM3 : 
    do_turningbands(&KEY[*keyNr],res END_WITH_RANDOM); break;   
  case SpectralTBM : 
    do_simulatespectral(&KEY[*keyNr],res END_WITH_RANDOM); break;
  case Direct : 
    do_directGauss(&KEY[*keyNr],res END_WITH_RANDOM);  break;    
  case Nugget : 
    *error=0;{int i;for(i=0;i<KEY[*keyNr].totallength;i++){res[i]=0.0;}}break;
  case AdditiveMpp: 
  case MaxMpp: /* **************** */
    do_addmpp(&KEY[*keyNr],res END_WITH_RANDOM); break;
  case Hyperplane : 
    *error=ERRORNOTPROGRAMMED; break;
  case Special : 
    KEY[*keyNr].cov->other(&KEY[*keyNr],res END_WITH_RANDOM); break;
  default: 	
    assert(false);
  }
  if (*error) goto ErrorHandling;
  mean =  KEY[*keyNr].param[MEAN];
  if ((!KEY[*keyNr].nuggetincluded) && (KEY[*keyNr].param[NUGGET]>0)) { //nugget
    Real sqrttwonugget,UU,VV;
    sqrttwonugget = sqrt(KEY[*keyNr].param[NUGGET]*2.0);  
    for (nx=0,endfor=2*(KEY[*keyNr].totallength/2); nx<endfor; nx++) {
      UU = sqrttwonugget * sqrt(-log(1.0 - UNIFORM_RANDOM));
      VV = TWOPI * UNIFORM_RANDOM; 
      res[nx] +=  UU * sin(VV) + mean;
      nx++;
      res[nx] += UU * cos(VV)+mean;
    }
    if (endfor<KEY[*keyNr].totallength) {
      do UU=UNIFORM_RANDOM; while (UU==0.0); 
      UU = sqrttwonugget * sqrt(-log(1.0 - UNIFORM_RANDOM));
      VV = TWOPI * UNIFORM_RANDOM; 
      res[endfor] += UU * sin(VV)+ mean;
    }
  } else { // no nugget 
    if (mean!=0.0) {
      for (nx=0,endfor=KEY[*keyNr].totallength; nx<endfor; nx++) res[nx] += mean;
    }
  }
#ifndef RF_GSL
    PutRNGstate();
#endif
  return;
 ErrorHandling:
  KEY[*keyNr].active=false;  
  if (GENERAL_PRINTLEVEL>0) ErrorMessage(KEY[*keyNr].method,*error); 
  return;
}

void SimulateRF(Real *x, Real *y, Real *z, int *dim,int *lx, int *grid, 
		int *covnr, Real *ParamList, int *nParam,
		int *method, int *distr,
		int *keyNr,
		Real *res, int *error)
{
  // grid   : boolean; if true then x,y,z are supposed to have three values: 
  //                      start, end, step
  //                   if false then (x,y,z) is a list of points where the 
  //                      value of the random field should be simulated
  // x,y,z  : see grid; if z==NULL it is assumed that the field is at most 
  //                    2 dimensional
  //                    if y==z==NULL, the random field is 1 dimensional
  // covnr  : nr of covariance function by getCovFct,
  // ParamList: vector with 5 to 8 elements
  //          see RFsimu.h, the definition of MEAN..KAPPAXX, for the meaning 
  //          of the components
  // method : if <0 then programme chooses automatically for the method
  //          if ==-1 it starts with the preferred ("first") method for the 
  //                  specific covariance fct
  //            else starts search with previous successful method (if any)
  //          otherwise see list `SimulationType' in RFsimu.h
  // keyNR  : label where intermediate results should be stored
  // res : array/grid of simulation results

  InitSimulateRF(x,y,z,dim,lx,grid,covnr,ParamList,nParam,method,
		 distr,
		 keyNr,error);
  if (error<=0) {
    DoSimulateRF(keyNr,res,error);
  } else {    
    *error=ERRORFAILED;
    if (GENERAL_PRINTLEVEL>0) PRINTF(" ** All methods have failed. **\n");
  } 
}
  

void GetRandomSize(int *l)
{
#ifdef RF_GSL
  struct timeb tp;   
  ftime(&tp);
  if (RANDOM==NULL) {    
    gsl_rng_default_seed =  tp.time * 257 + tp.millitm; 
    RANDOM = gsl_rng_alloc(RANDOMNUMBERGENERATOR);      
  }
  *l=gsl_rng_size(RANDOM);
#else 
  *l = -1;
#endif
}


#ifdef RF_GSL
void StoreSeed()
{ // gsl_library has no getrandomseed-function. so we restart the generator with 
  //a new seed
  struct timeb tp;   
  ftime(&tp);
  if (RANDOM==NULL) {    
    gsl_rng_default_seed =  tp.time * 257 + tp.millitm;
    RANDOM = gsl_rng_alloc(RANDOMNUMBERGENERATOR);     
  }
  if (RF_RNG_STORED!=NULL) gsl_rng_free(RF_RNG_STORED);
  RF_RNG_STORED=gsl_rng_clone(RANDOM);
}
void RestoreSeed()
{
  gsl_rng_free(RANDOM);
  RANDOM=gsl_rng_clone(RF_RNG_STORED);
}

void XGetSeed(int *seed,int *l)
{ // gsl_library has no getrandomseed-function. so we restart the generator
  // with a new seed
  struct timeb tp;   
  ftime(&tp);
  if (RANDOM==NULL) {    
    gsl_rng_default_seed =  tp.time * 257 + tp.millitm; 
    RANDOM = gsl_rng_alloc(RANDOMNUMBERGENERATOR);      
  }
  *l = gsl_rng_size(RANDOM);
  if (seed!=NULL) free(seed);
  seed = (int*) malloc(sizeof(int) * *l);
  memcpy(seed,RANDOM->state,*l);
}

void SetSeed(int *seed, int *l, int *error)
{ // gsl_library has no getrandomseed-function. so we restart the generator 
  // with a new seed
  gsl_rng_free(RANDOM);
  RANDOM = gsl_rng_alloc(RANDOMNUMBERGENERATOR);
  if (*l<gsl_rng_size(RANDOM)) {*error=1; return;}
  memcpy(RANDOM->state,seed,gsl_rng_size(RANDOM));
  *error=0;
}

void GetSeed(int *seed, int *l, int *error)
{
  struct timeb tp;   
  ftime(&tp);
  if (RANDOM==NULL) {    
    gsl_rng_default_seed =  tp.time * 257 + tp.millitm; 
    RANDOM = gsl_rng_alloc(RANDOMNUMBERGENERATOR);      
  }
  if (*l<gsl_rng_size(RANDOM)) {*error=1; return;}
  else (*l=gsl_rng_size(RANDOM));
  //gsl_rng_print_state(RANDOM);
  memcpy(seed,RANDOM->state,*l);
  *error=0;
}

#endif   /* of RF_GSL */
