
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part and error messages

 Copyright (C) 2017 -- 2018 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.i

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#include "def.h"
#include "operator.h"


// todo: Folgendes kann nur mit extremen Aufwand fuer \pkg{parallel}
// sicher gemacht werden
#ifdef DO_PARALLEL
// hier
#else  // not DO_PARALLEL
char ERRMSG[LENERRMSG], MSG[LENERRMSG], MSG2[LENERRMSG];
#endif

#define SCPY(A,B) STRNCPY(A,B,LENERRMSG);

void errorMSG(int err, errorstring_type errorstring, KEY_type *KT,
	      char* M, int len) {
  char m[LENERRMSG];
  //if (err >= ERRORM && err <= ERRORMEND) err = ERRORM;
 
  
  switch (err) {
  case NOERROR : SCPY(m,"none"); break;
  case NOERROR_REPEAT : SCPY(m,"none; looking for further covariances applicable to the same method");break;
  case NOERROR_ENDOFLIST : SCPY(m,"none; end of list");break;
  case XERRORNOTDEFINED :       
    SCPY(m,"specified method undefined for the given model or no simulation method found for the given model");break;
  case ERRORM: case XERRORCHANGESYSTEM:
    SCPY(m, errorstring);
    break;
    //  case XERRORTRAFOSYSTEM:
    //SCPY(m,"coordinate systems do not match without (internal) transformation");
    //   break;
  case XERRORCARTESIAN :
    SCPY(m,"only cartesian system allowed (currently)"); break;
   case XERRORVDIMNOTPROGRAMMEDYET :    
    SCPY(m,"multivariate version not programmed yet. Sorry."); break;
  case XERRORTYPECONSISTENCY :
     SCPY(m,"incorrect choice of submodel (type inconsistency)"); break;
  case XERRORNOTINITIALIZED: 
    SCPY(m,"not initialized or storing=FALSE");break; // OK
  case XERRORSVD:
    SCPY(m,"SVD decomposition failed");break;
  case XERRORDECOMPOSITION:
    SCPY(m,"covariance function does not seem to be (strictly) positive definite");break;
  case XERRORNOMULTIVARIATE :
    SCPY(m, "multivariate models not allowed (yet)"); 
    break;
  case XERROR_MATRIX_SQUARE :
    SCPY(m, "square matrix expected"); break;
  case XERROR_MATRIX_VDIM :
    SCPY(m, "size of matrix is not a multiple of the multivariate dimension"); break;
  case XERROR_MATRIX_POSDEF :
    SCPY(m, "matrix does not seem to be strictly positive definite"); break;
    //  case XERROR_MATRIX_ :   SCPY(m, ""); break;
  case XERRORCEDIM: 
    //
    //    { printf("error dimension\n"); model *cov; crash(cov); }
    SPRINTF(m,"dimension specification not in [1,%d] or dimension of coordinates larger than that the supposed spatio-temporal process",
	    MAXCEDIM);break;
  case XERRORSORTOF :
   SCPY(m, "argument 'whichparam' out of range."); break;
    //  case XERROR_MATRIX_ :   SCPY(m, ""); break;
  case ERRORDIM: 
    //
    //    { printf("error dimension\n"); model *cov; crash(cov); }
    SPRINTF(m,"dimension specification not correct or dimension of coordinates larger than that the supposed spatio-temporal process");break;
  case XERRORWAVING :
    SCPY(m,"Rescaling not possible (waving or large nugget effect?)");break;
  case XERRORRESCALING:
    SCPY(m,"practical range not defined");
    break;
  case XERRORANISO_T :
    SPRINTF(m, "'%.50s' may not be given at the same time with '%.50s' or '%.50s'", 
	    DefList[DOLLAR].kappanames[DANISO], 
	    DefList[DOLLAR].kappanames[DAUSER], 
	    DefList[DOLLAR].kappanames[DPROJ]);
    break;
  case XERRORANISO:
    SCPY(m,"anisotropy parameter not allowed"); break; 
 case XERRORANISO_DET:
    SCPY(m,"determinant of the anisotropy matrix could not be calculated.");
    break; 
  case XERRORANISO_SQUARE:
    SCPY(m,"only square anisotropic matrices allowed");
    break; 
  case XERRORANISO_INV:
    SCPY(m,"inversion of anisotropy matrix failed");
    break; 
  case XERRORFOURIER:
    SCPY(m,"fourier transformation has failed");
    break;
  case XERRORTOOMANYLOC:
    SCPY(m,"fourier transformation has failed");
    break;
  case XERRORNONSTATSCALE:
    SCPY(m,"unstationary scale not allowed yet for RMS con-struct-ions");
    break;
    
  case XERRORNOSTATMATCH : 
    SCPY(m,"no matching assumption found for the domains");
    break;
  case XERRORUNKNOWNMETHOD:
    SCPY(m,"Unknown method in this context or unallowed mixture of methods"); 
    break;
  case XERRORWRONGDIM:
    SCPY(m,"wrong dimension"); break;
   case XERRORUNKOWNSXPTYPE:
    SCPY(m, "parameter value of unknown SXP type");
    break;
  case XERROROUTOFMETHODLIST:
    char restrictive[100], info[500];
    SPRINTF(restrictive, "Are the %.50s() too restrictive?", RFOPTIONS);
    SPRINTF(info, "\n You get (more) internal information if you set %.10s(%.50s=%d) before running your code.",
	    RFOPTIONS, 
	    "cPrintlevel",
	    PL_DETAILSUSER);
    
    SPRINTF(m, 
	    "Running out of list of methods. %.100s%.50s",
	    GLOBAL_UTILS->basic.skipchecks
	    ? "Did you try an invalid parameter combination?" : restrictive,
	    PL <= 2 ? info : "" );
    break;
  case XERRORWRONGISO: 
    SCPY(m, "unexpected appearance of a rather general function. Should the model be isotropic, but isn't?")
    break;
  case XERRORWRONGDOM:
    SCPY(m, "unexpected appearance of a kernal function. Should the model be stationary, but isn't?");
    break;
  case XERRORKERNEL:
    SCPY(m, "the mapping 'earth -> cartesian' keeps definiteness only if it is used as a kernel.");
    break;
  case XERRORBADVDIM: 
    SCPY(m, "m-dimensionality could not be detected");
    break;

  case XERRORNOTCARTESIAN:
       SCPY(m, "Currently only cartesian coordinate systems are allowed");
    break;
  case XERRORODDCOORDTRAFO:
    SCPY(m, "coordinate transformation not possible or not programmed yet");
    break;
  case XERRORFULLCOORDNEEDED:
    SCPY(m, "full coordinate system needed");
    break;
  case XERRORSUBMETHODFAILED:
    SPRINTF(m, "no good submethods exist");
    break;
  case  XERRORSTATVARIO:
    SCPY(m, 
	   "negative definite function expected depending on 1 variable only");
    break;
  case  XERRORREDUCED:
    SCPY(m, "calling system not unreduced");
    break;
   case XERRORNOVARIOGRAM:
    SCPY(m, "Variogram model not allowed in this context");
    break;
  case XERRORNORMALMIXTURE:
    SCPY(m, "only normal mixtures as first submodel allowed (Gneiting, 2002)");
    break;
  case XERRORMAXDIMMETH:
    SCPY(m, "maximal dimension of variables for the method exceeded");
    break;
  case XERRORPREVDOLLAR:
    SCPY(m, "method may not be initialised by preceding initS");
    break;
  case XERRORSPECTRAL: 
    SCPY(m, "submodel does not have spectral representation");
    break;    
  case XERRORTBMCOMBI: 
    SCPY(m, "the given combination of 'fulldim' and 'reduceddim' is not possible yet.");
    break;    
 case XERRORMAXVDIM: 
   SPRINTF(m, "maximum of multivariate components (%d) exceeded. If necessary change the value od 'MAXVDIM' in 'RandomFields/src/AutoRandomFields.cc' and reinstall the package. Note that starting with MAXVDIM > 10^3 the 'RandomFields' needs a lot of memory.", MAXVDIM);
    break;

  case XERRORINVALIDMODEL : // gauss distribution, no method
    SCPY(m, "Invalid covariance model: did you wrongly use an auxiliary function to construct the model?");
    break;    
  case XERRORODDMODEL : // gauss distribution, no method
    SCPY(m, "Odd covariance model: the use of auxiliary functions and/or your choice of the parameters lead to a covariance model for which no simulation methods exist.");
    break;    
  case XERRORDIAMETERNOTGIVEN:
    SCPY(m, "Diameter must always be given");
    break;
  case XERRORPREFNONE:
    SCPY(m, "the simulation method does not allow for the given model.");
    break;
  case XERRORPREFNONECOV:
    SCPY(m, "the given model does not allow for calculation the covariance values.");
    break;
 case XERRORRANDOMKAPPA:
    SCPY(m, "Only (sub)models with deterministic parameters allowed.");
    break;
   
    //    case : SCPY(m,"");break;
    //
    // Poisson:
  case XERRORUNKNOWNMAXTYPE :
    SCPY(m, "unknown type of max-stable process");
    break;
 
  case XERRORATOMP :
    SCPY(m, "p must be given everywhere or nowhere");
    break;
   
  case XERRORKRIGETOL :
    SCPY(m,"sigma must be at most KRIGE_TOLERANCE");
    break;


  case MSGLOCAL_OK :
    SCPY(m,"fine");
    break;
  case XMSGLOCAL_JUSTTRY :
    SCPY(m,
	   "unclear whether algorithm will work for specified parameters");
    break;
  case XMSGLOCAL_NUMOK :
    SCPY(m,"fine. Algorithm should work for specified parameters");
    break;
  case XMSGLOCAL_ENDOFLIST :
    SCPY(m,"end of list for variants of the algorithm");
    break;
  case XMSGLOCAL_SIGNPHI :
    SCPY(m,"wrong sign of covariance function at break point");
    break;
  case XMSGLOCAL_SIGNPHIFST :
    SCPY(m, "wrong sign of 1st derivative of the covariance function at the break point");
    break;
  case XMSGLOCAL_SIGNPHISND :
    SCPY(m, "wrong sign of 2nd derivative of the covariance function at the break point");
    break;
 case XMSGLOCAL_NOPOSITIVEROOT :
    SCPY(m, "algorithm faild to find positive radius");
    break;
  case XMSGLOCAL_INITINTRINSIC :
    SCPY(m,"one of a2, b or a0+phi(0) has wrong sign");
    break;
  case XERRORUNSPECIFIED :
    SCPY(m,"(unspecified)");
    break;
  case ERRORNOTPROGRAMMEDYET :    // imported from Utils 
    SCPY(m,"Not programmed yet in RandomFields Version 3. Sorry."); break;
  case ERRORFAILED:               // imported from Utils 
   SCPY(m,"algorithm failed (partially)");break;
  case ERRORMEMORYALLOCATION:     // imported from Utils 
    SCPY(m, "memory allocation error -- too much space demanded or non-positive number of bytes requested"); 
    break;
    // case ERRORDUMMY : SCPY(m,"none (dummy)"); break; 
  case ERRORREGISTER: 
    SCPY(m, "register number out of range or model definition empty");
    break;

  case TOOLS_DIM :
    SCPY(m,"dimension exceed max dimension of empirical variogram estimation");
     break;
   case TOOLS_XERROR :  
    SCPY(m,"The x coordinate may not be NULL.\n");
     break;
  case TOOLS_BIN_ERROR :
    SCPY(m,"Bin components not an increasing sequence.\n");
     break;
  case TOOLS_UNKNOWN_CHAR:
    SCPY(m,"unknown type of second order characteristic");
     break;
  case TOOLS_METHOD:
    SCPY(m,"unknown method in empirical variogram calculation. Please contact author.");
    break;
  case XMSGLOCAL_FAILED :
     SCPY(m,"creation of the cut-off model failed.");
     break;   
  case XMSGLOCAL_WRONGRADII :
     SCPY(m,"Creation the cut-off model failed.");
     break;   
  default : 
     PRINTF(" error=%d\n", err); 
     // crash();
     BUG;
  }

  char m2[LENERRMSG];
  if (KT != NULL) SPRINTF(m2, "%.100s %.800s", KT->error_loc, m);
  else SCPY(m2, m);
  
  if (STRLEN(m) > (unsigned int) len && len > 6) {    
    //  printf("%s %d %d\n", m, STRLEN(m), len);
    m2[len-2] = m2[len-3] = m2[len-4] = '.';
    m2[len-5] = ' ';
    m2[len-1] ='\0';
    // printf("out %s %d %d\n", m, STRLEN(m), len);
  }
  STRNCPY(M, m2, MAXERRORSTRING);

  // printf("!!!! err=%d M=%s m2=%s m=%s\n", err, M, m2, m);
  
  if (PL >= PL_ERRORS) { 
    PRINTF("err code %d [%s]\n", err, m2);
  }
}



void errorMSG(int err, char* m, errorstring_type err_msg) {
  errorMSG(err, err_msg, NULL, m, MAXERRORSTRING);}

void errorMSG(int err, errorstring_type err_msg, KEY_type *KT, char* m) {
  errorMSG(err, err_msg, KT, m, MAXERRORSTRING); }


//errorstring_type err_msg = "error could not be identified in multicore modus. Recompile without OMP to see the error msg.";
errorstring_type err_msg = "error could not be identified.";
void errorMSG(int err, char* m) {
  errorMSG(err, err_msg, NULL, m, MAXERRORSTRING);
}

void errorMSG(int err, KEY_type *KT, char* m) {
  errorMSG(err, err_msg, KT, m, MAXERRORSTRING);
}


void OnErrorStop(int Err, errorstring_type errmsg) {
  if (Err == NOERROR) return;
  char m[MAXERRORSTRING];
  errorMSG(Err, errmsg, NULL, m, MAXERRORSTRING);
  RFERROR(m);
}

void OnErrorStop(int Err, errorstring_type errmsg, KEY_type *KT) {
  if (Err == NOERROR) return;
  char m[MAXERRORSTRING];
  errorMSG(Err, errmsg, KT, m, MAXERRORSTRING);
  RFERROR(m);
}

void OnErrorStop(int Err, model *cov) {
  if (Err == NOERROR) return;
  char m[MAXERRORSTRING];
  errorMSG(Err, cov->err_msg, cov->base, m, MAXERRORSTRING);
  RFERROR(m);
}

