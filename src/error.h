#ifndef RFerror_H
#define RFerror_H 1


// Error codes & messages
// positive values are error codes
// negative values are messages
#define NOERROR 0                 
#define ERRORNOTDEFINED 1        /* the specification for the  covariance and 
				    method is not given/known, e.g. TBM2 for 
				    many covariance functions  */
#define ERRORNOTPROGRAMMED 3     
#define ERRORCOVNOTALLOWED 4     /* wrong dimension for the specified covariance 
				    function */
#define ERRORFAILED 5            /* method didn't work for the specified 
				    parameters */
#define ERRORDECOMPOSITION 9     /* direct simulation; matrix inversion failed */
#define ERRORPRECISION 10        /* direct simulation; matrix inversion didn't 
				    reach the desired presicion */
#define ERRORCOVFAILED 14        /* covariance model not defined for specified 
				    parameters and method 
				    -- NOTE: ERRORSTRING_OK and 
				    ERRORSTRING_WRONG must be set*/
#define ERRORNN 15              /* nn in TBM too large */

// for getnaturalscaling !!
#define ERRORWAVING 16           /* numerical determination of practical range
				    fails */
#define ERRORRESCALING 17        /* attempt to rescale where not possible */

 
#define ERRORM 20                /* a single error message */
#define ERRORMEND 30         /* a single error message -- und alles dazwischen */

#define ERRORISOTROPICMETHOD 31  /* TBM might be called only by isotropic
                                    functions */
#define ERRORTIMESEPARATE 32     /* last line of matrix 
				  */
#define ERRORWRONGDIM 33        /* (tbm) dimension couldn't be reduced to be 
                                   appropriate for the respective method */
#define ERRORANISOMIX 34         /* diff. anisotropies for multiplication 
				    covariance and tbm */
#define ERRORTOOMANYPOINTS 37   /* too many points to determine smallest 
				   distance between within TBM! */
#define ERRORTIMENOTALLOWED 38  /* time component not allowed for the specified 
                                   method */
#define ERRORCOVNUMERICAL 39    /* covariance function of derivative cannot
				   be calculated exactly, e.g. in TBM2 -- 
				   NOTE: ERRORSTRING_OK and 
				   ERRORSTRING_WRONG must be set */
#define ERRORCEMAXMEMORY 40       /* in circulant embedding: circulant vector too
				   large  -- NOTE: ERRORSTRING_OK and 
				    ERRORSTRING_WRONG must be set */
#define ERRORTOOMANYLINES 41    /* Hyperplane tesselation: estimated simulated
				   lines passes given threshold*/
#define ERRORLAYERS 42           /* value of TBM*.layers does not match 
				   the situation */
#define ERRORMETHODEXCLUDED 43  /* method explicitely excluded by other 
				   parameter; here: stationary_only=TRUE */
#define ERRORTBMCENTER 45       /* given values for TBM center contain NAs */
#define ERRORTBMPOINTS 46       /* given number of points is too small */
#define ERRORHYPERMETHOD 48         /* no of announced submodels invalid */
#define ERRORHYPERNOTALLOWED 49 /* most methods do not allow for hypermodels 
				   yet */
#define ERRORISO    50         /* isotropic model not allowed -- internal warning */ 
#define ERRORMSG 52            /* parenthesis as operator iff isohypermodel 
				  -- NOTE: ERRORSTRING_OK and 
				  ERRORSTRING_WRONG must be set*/
#define ERRORTBMPOINTSZERO 53  /* linesimustep and linesimufactor are zero 
				   and TBM_POINTS<4 */
#define ERRORFULLRANK 54       /* anisotropy matrix does not have full rank,
				  but needed in circ_embed_local */
#define ERRORTIMECOMPONENT 55  /* some hyper models need vanishing temporal
				  component in the anisotropy matrix*/
#define ERRORLOWRANK 56   
#define ERRORLOWRANKTBM 57       /* tbm does not allow that the matrix `a' of 
                                   a product matrix are all zero
				 */ 
#define ERRORSTATVARIO 58       /* only stat or variogram allowed as submodel */
#define ERRORONLYISOTROPIC 59   /* only isotropic random fields -- only
				   for Markov implementation */
#define ERRORMARKOVNOTINCLUDED 60 /* see includeMarkov.h for details */
#define ERRORONLYREGULARGRID 61 /* markov */
#define ERRORMARKOVPARAMETER 62 /*parameter set not available for GMRF method */ 
#define ERRORMARKOVMAXMEMORY 63   /* too many points for markov
				     -- NOTE: ERRORSTRING_OK and 
				    ERRORSTRING_WRONG must be set */
#define ERRORTIMECOMPLETESEPARATE 64    /* last line and row of matrix 
				  */
#define ERRORLASTGRID 65       /* time or grid */
#define ERRORTRIVIALTIMEDIRECTION 66       /* time or grid */
#define ERRORMETROPOLIS     67     /* Metropolis search alg. failed */
#define ERRORLOGMIX 68           /* log mix not defined */
#define ERRORNOVARIOGRAM 69     /* variogram models not allowed in 
				   ci and direct */
#define ERRORDIAG 70            /* diagonal matrices needed */
#define ERRORSUBMETHODFAILED 71 /* aufsplitten der Modelle bringt nichts */
#define ERRORUNITVAR  72        /* unit variance required */
#define ERRORMAXDIMMETH 73     /* max dimension of method exceeded */
#define ERRORPREVDOLLAR 74      /* one of the prev method was dollar,
				   but cannot be taken into account */
#define ERRORNOGENUINEMETHOD 75  /* the leave of a method tree is not 
				  randomcoin, see extremes.cc */
#define ERRORANYLEFT 96         /* some covariances have not been simulated yet
				  -- this error code may be used only in #
				  InitSimulation */
#define NOERROR_REPEAT 97
#define NOERROR_ENDOFLIST 98

#define ERRORDUMMY 99            /* no error, only for tracing purposes */

#define ERRORCURRENTLY 100    /* currently a fucntionality is not available */


// *** severe errors
#define ERRORNORMALMIXTURE 101    /* normal mixture required */
#define ERRORNOMULTIVARIATE 102   /* vdim > 1 not allowed */
#define ERRORMULTIMISMATCH 104    /* another multivariate dimension expected */
#define ERRORAUXVECTOR 105    /* another multivariate dimension expected */
#define ERRORMEMORYALLOCATION 106 /* malloc returned NULL pointer */
#define ERRORNOTINITIALIZED 107   /* key.active==false in DoSimulateRF; is only 
				    checked there !!! */
#define ERRORKEYNROUTOFRANGE 108   /* wrong keynr */
#define ERRORVARIOGRAMONLY 112    /* attempt to get the covariance whereas only 
				    the variogram is defined */
#define ERRORFOURIER 113          /* specific error in Fourier transformation */
#define ERRORREGISTER 115         /* wrong register number */
#define ERRORCOORDINATES 116      /* coordinates are not given or invalid grid 
				    specification */
#define ERRORNEGATIVEVAR 117
#define ERRORDIM 119              /* dim<1 or dim>MAXDIM */
#define ERRORNEGATIVESCALE 121
#define ERRORWRONGINIT 122        /* dosimulate and initsimulate do not match */
#define ERRORANISOTROPIC 124      /* trying to call anisotropic fct with 
				    isotropic parameters */
#define ERRORTIMENOTANISO 126      /* if time is used then, anisotropic must
				    be true */
#define ERRORANISO 127            /* anisotropic call not allowed */
#define ERRORCIRCNONPOS 128       /* circ embed: fourier trafo not pos def */
#define ERRORUNKNOWNMETHOD 131  

#define ERRORWITHOUTTIME 132     /* spatially isotropic covariance functions
				   must always be called with a time 
				   component */
#define ERRORCOVNROUTOFRANGE 133 /* not allowed covariance number */
#define ERRORNOSTATMATCH 134     /* no matching ISO/STAT found in check2 */
#define ERRORNOTANISO 135        /* make sure that SPACEISOTROPIC and ANISOTROPIC
                                   covariance function are called by anisotropy
				   definition */
#define ERRORDIMMISMATCH 136    /* logical dim does not match physical dim in the
                                   covariance function call */
#define ERROROUTOFMETHODLIST 144 /* no further method we can try is available */
#define ERRORHYPERNR 147         /* no of announced submodels invalid */
#define ERRORHANGING 148         /* hanging covariance (#), but
				    method can not deal with that */
#define ERRORNCOVOUTOFRANGE 150 /* number of cov. functions not in 0/1..MAXDIM*/
#define ERRORCOVUNKNOWN 154    

#define ERRORMATRIX_VECTOR 160 /* must be vector */
#define ERRORMATRIX_Z 161 /* wrong values for Z*/


// **** extremes
#define ERRORSILLNULL 201        /* extremes : sill must be positive */
#define ERRORPAIRS 202          /* extremes : paired not allowed */
#define ERRORTREND 203           /* extremes: no trend allowed */


// *** poisson
#define ERRORVARMEAN 301         /* Poisson distr.: mean<>variance */
/* do not use numbers 800 -- 900 : reserved to MPP package */

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
#define ERRORUNSPECIFIED 999  

#endif
