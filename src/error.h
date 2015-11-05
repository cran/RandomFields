


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 Martin Schlather

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


#ifndef RFerror_H
#define RFerror_H 1

extern char ERRMSG[LENERRMSG], MSG[LENERRMSG], BUG_MSG[250], MSG2[LENERRMSG],
  ERRORSTRING[MAXERRORSTRING], ERROR_LOC[nErrorLoc];
extern int ERRORMODELNUMBER;
void errorMSG(int error, char* EM);
void FinalErrorMSG(int err, char* m);
void ErrorStop(int err);

// Error codes & messages
// positive values are error codes
// negative values are messages

// #define NOERROR 0                 
// #define ERRORMEMORYALLOCATION 1 /* NEVER CHANGE VALUE; MALLOC returned NULL pointer */



// LESS THAN 10 IS RESERVED TO BE JOINTLY USED BY RFUTILS


#define ERRORMEND 11         /* a single error message -- und alles dazwischen */
#define ERRORCARTESIAN 12

#define ERRORTYPECONSISTENCY 14
#define ERRORWRONGVDIM 15
#define ERRORBADVDIM 16
#define ERRORNOTCARTESIAN 17
#define ERRORODDCOORDTRAFO 18
#define ERRORVDIMNOTPROGRAMMEDYET 19

#define ERRORNOTDEFINED 20       /* the specification for the  covariance and 
				    method is not given/known, e.g. TBM2 for 
				    many covariance functions  */
#define ERRORUNKOWNSXPTYPE 22   /* parameter value of unknown SXP type */
#define ERRORSTATVARIO 23      /* only stat or variogram allowed as submodel */
#define ERRORANISO_T 24         /* anisoT' may not be given at the same
				time with 'Aniso' or 'proj' */
//#define ERRORVARIOGRAMONLY 25    /* attempt to get the covariance whereas only 				    the variogram is defined */
#define ERRORDECOMPOSITION 26   /* direct.cc */
#define ERRORPREFNONE 27
#define ERRORPREFNONECOV 28

#define ERRORSUBMETHODFAILED 30 /* aufsplitten der Modelle bringt nichts */
#define ERRORMAXDIMMETH 31     /* max dimension of method exceeded */
#define ERRORPREVDOLLAR 32      /* one of the prev method was dollar,
				   but cannot be taken into account */
#define ERRORINVALIDMODEL 33 /* no gauss method available, including 'Nothing'*/
#define ERRORODDMODEL 34     /* no gauss method available */
#define ERROROUTOFMETHODLIST 35 /* no further method we can try is available */

#define ERRORREGISTER 36
#define ERRORWRONGISO 37
#define ERRORKERNEL 38

#define ERRORWRONGDIM 40        /* (tbm) dimension couldn't be reduced to be 
                                   appropriate for the respective method */
#define ERRORTOOMANYLINES 41    /* Hyperplane tesselation: estimated simulated
				   lines passes given threshold*/
#define ERRORNOVARIOGRAM 43     /* variogram models not allowed in 
				   ci and direct */
#define ERRORSPECTRAL 44
#define ERRORTBMCOMBI 45



// for getnaturalscaling !!
#define ERRORWAVING 50           /* numerical determination of practical range
				    fails */
#define ERRORRESCALING 51        /* attempt to rescale where not possible */

 

#define NOERROR_REPEAT 97
#define NOERROR_ENDOFLIST 98


#define ERRORDUMMY 99            /* no error, only for tracing purposes */


// *** severe errors
#define ERRORNORMALMIXTURE 101    /* normal mixture required */
#define ERRORNOMULTIVARIATE 102   /* vdim > 1 not allowed */
#define ERROR_MATRIX_SQUARE 103
#define ERROR_MATRIX_VDIM 104
#define ERROR_MATRIX_POSDEF 105
#define MATRIX_NOT_CHECK_YET -999  /* must be a negative value not used somewhere else ! */

//#define ERROR_MATRIX_ 10

#define ERRORNOTINITIALIZED 107   /* key.active==false in DoSimulateRF; is only 
				    checked there !!! */
#define ERRORDIM 119              /* dim<1 or dim>MAXDIM */
#define ERRORANISO 127            /* anisotropic call not allowed */
#define ERRORDIAMETERNOTGIVEN 128
#define ERRORUNKNOWNMETHOD 131  

#define ERRORNOSTATMATCH 134     /* no matching ISO/STAT found in check2 */


// **** extremes
#define ERRORTREND 203           /* extremes: no trend allowed */
#define ERRORUNKNOWNMAXTYPE 204 

//#define ERRORONLYATOMS 250   /* Poisson distribution */
#define ERRORANISODOLLARNOTALLOWED_YET 251
//#define ERRORATOMNOTPROGRAMMED_YET 252
#define ERRORATOMP 253


// *** kriging
#define ERRORKRIGETOL 350

// messages and strategies for local methods
// ENDOFLIST also used for any kind of error appearing in getparam
// never change order of the message numbers up to ENDOFLIST
#define MSGLOCAL_OK 400 /* do not change ordering! (the less the better) */
#define MSGLOCAL_NUMOK 401
#define MSGLOCAL_JUSTTRY 402
#define MSGLOCAL_ENDOFLIST 403
#define MSGLOCAL_SIGNPHI 404
#define MSGLOCAL_SIGNPHIFST 405
#define MSGLOCAL_SIGNPHISND 406
#define MSGLOCAL_INITINTRINSIC 407
#define MSGLOCAL_FAILED 408
#define MSGLOCAL_NOPOSITIVEROOT 409
#define ERRORUNSPECIFIED 999  

/* do not use numbers 800 -- 900 : reserved to MPP package */



#endif

