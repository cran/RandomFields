#ifndef RFerror_H
#define RFerror_H 1


// Error codes & messages
// positive values are error codes
// negative values are messages
#define USEOLDINIT -1            /* not an error, just a message !*/
#define NOERROR 0                 
#define ERRORNOTDEFINED 1        /* the specification for the  covariance and 
				    method is not given/known, e.g. TBM2 for 
				    many covariance functions  */
#define ERRORMETHODNOTALLOWED 2  /* wrong method for the specified covariance 
				    function or grid */
#define ERRORNOTPROGRAMMED 3     
#define ERRORCOVNOTALLOWED 4     /* wrong dimension for the specified covariance 
				    function */
#define ERRORFAILED 5            /* method didn't work for the specified 
				    parameters */
#define ERRORMEMORYALLOCATION 6  /* malloc returned NULL pointer */
#define ERRORNOTINITIALIZED 7    /* key.active==false in DoSimulateRF; is only 
				    checked there !!! */
#define ERRORKEYNROUTOFRANGE 8   /* wrong keynr */
#define ERRORDECOMPOSITION 9     /* direct simulation; matrix inversion failed */
#define ERRORPRECISION 10        /* direct simulation; matrix inversion didn't 
				    reach the desired presicion */
#define ERRORRESCALING 11        /* attempt to rescale where not possible */
#define ERRORVARIOGRAMONLY 12    /* attempt to get the covariance whereas only 
				    the variogram is defined */
#define ERRORFOURIER 13          /* specific error in Fourier transformation */
#define ERRORCOVFAILED 14        /* covariance model not defined for specified 
				    parameters */
#define ERRORREGISTER 15         /* wrong register number */
#define ERRORCOORDINATES 16      /* coordinates are not given or invalid grid 
				    specification */
#define ERRORNEGATIVEVAR 17
#define ERRORPARAMNUMBER 18
#define ERRORDIM 19              /* dim<1 or dim>MAXDIM */
#define ERRORNEGATIVESCALE 20
#define ERRORWAVING 21           /* numerical determination of practical range
				    fails */
#define ERRORWRONGINIT 22        /* dosimulate and initsimulate do not match */
#define ERRORNN 23               /* nn in TBM too large */
#define ERRORANISOTROPIC 24      /* trying to call anisotropic fct with 
				    isotropic parameters */
#define ERRORMETHODMIX 25        /* diff. methods for multiplication 
				    covariance */
#define ERRORTIMENOTANISO 26      /* if time is used then, anisotropic must
				    be true */
#define ERRORNOMULTIPLICATION 27 /* method not allowed for multiplicative 
				    covariance functions*/
#define ERRORANISOMIX 28         /* diff. anisotropies for multiplication 
				    covariance and tbm */
#define ERRORISOTROPICMETHOD 29  /* TBM might be called only by isotropic
                                    functions */
#define ERRORTIMESEPARATE 30     /* last line of matrix 
				  */
#define ERRORUNKNOWNMETHOD 31  

#define ERRORWITHOUTTIME 32     /* spatially isotropic covariance functions
				   must always be called with a time 
				   component */
#define ERRORCOVNROUTOFRANGE 33 /* not allowed covariance number */
#define ERRORWRONGDIM 34        /* (tbm) dimension couldn't be reduced to be 
                                   appropriate for the respective method */
#define ERRORNOTANISO 35        /* make sure that SPACEISOTROPIC and ANISOTROPIC
                                   covariance function are called by anisotropy
				   definition */
#define ERRORDIMMISMATCH 36     /* logical dim does not match physical dim in the
                                   covariance function call */
#define ERRORTOOMANYPOINTS 37   /* two many points to determine smallest 
				   distance between within TBM! */
#define ERRORTIMENOTALLOWED 38  /* time component not allowed for the specified 
                                   method */
#define ERRORCOVNUMERICAL 39    /* covariance function of derivative cannot
				   be calculated exactly, e.g. in TBM2 */
#define ERRORMAXMEMORY 40       /* in circulant embedding: circulant vector too
				   large */
#define ERRORTOOMANYLINES 41    /* Hyperplane tesselation: estimated simulated
				   lines passes given threshold*/
#define ERRORANYLEFT 42         /* some covariances have not been simulated yet
				 */
#define ERRORMETHODEXCLUDED 43  /* method explicitely excluded by other 
				   parameter; here: stationary_only=TRUE */
#define ERROROUTOFMETHODLIST 44 /* no further method we can try is available */
#define ERRORTBMCENTER 45       /* given values for TBM center contain NAs */
#define ERRORTBMPOINTS 46       /* given number of points is too small */
#define ERRORHYPERNR 47         /* no of announced submodels invalid */
#define ERRORHYPERMETHOD 48         /* no of announced submodels invalid */
#define ERRORHYPERNOTALLOWED 49 /* most methods do not allow for hypermodels 
				   yet */
#define ERRORNCOVOUTOFRANGE 50 /* number of cov. functions not in 0/1..MAXDIM*/
#define ERRORPARENTHESIS 51    /* parenthesis as operator iff isohypermodel */
#define ERRORMSG 52            /* parenthesis as operator iff isohypermodel */
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
#define ERRORHYPERNOTISO 58      /* only ISOTROPIC covariance models allowed
				    for submodels of hypermodels */

#define NOERROR_REPEAT 97
#define NOERROR_ENDOFLIST 98

#define ERRORDUMMY 99            /* no error, only for tracing purposes */

#define ERRORSILLNULL 101        /* extremes : sill must be positive */
#define ERRORPAIRS 102          /* extremes : paired not allowed */
#define ERRORTREND 103           /* extremes: no trend allowed */

#define ERRORVARMEAN 201         /* Poisson distr.: mean<>variance */
/* do not use numbers 800 -- 900 : reserved to MPP package */

// messages and strategies for local methods
// ENDOFLIST also used for any kind of error appearing in getparam
// never change order of the message numbers up to ENDOFLIST
#define MSGLOCAL_OK 300
#define MSGLOCAL_NUMOK 301
#define MSGLOCAL_JUSTTRY 302
#define MSGLOCAL_ENDOFLIST 303
#define MSGLOCAL_SIGNPHI 304
#define MSGLOCAL_SIGNPHIFST 305
#define MSGLOCAL_SIGNPHISND 306
#define MSGLOCAL_INITINTRINSIC 307
#define ERRORUNSPECIFIED 999  

#endif
