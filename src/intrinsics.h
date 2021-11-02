
#ifndef randomfields_initrinsics_H
#define randomfields_initrinsics_H 1

#include <inttypes.h> // uintptr_t

// PKG_CXXFLAGS =  $(SHLIB_OPENMP_CXXFLAGS) -mavx -msse -msse2 -msse3 -msse4 -msse4a -march=core-avx2

#ifdef asdfsdaf

#ifdef __MMX__
#define MMX __MMX__
#endif

#ifdef __SSE__
#define SSE __SSE__
#else
#endif


#ifdef  __SSE2__
#define SSE2 __SSE2__
#endif
#ifdef  __SSE3__
#define SSE3 __SSE3__
#endif
#ifdef  __SSSE3__
#define SSSE3 __SSSE3__
#endif
#ifdef  __SSE4A__
#define SSE4A __SSE4A__
#endif
#if defined __SSE41__ || defined __SS42__
#define SSE412 1
#endif

//
#ifdef __AVX__
#define AVX 1
#endif
#ifdef __AVX2__
#define AVX2 1
#else
#endif

#endif


#if defined (AVX512)
#define SSEBITS 512
#define SSEMODE 30
#elif defined (SSE)
#define SSEBITS 256
#define SSEMODE 20
#elif defined (SSE)
#define SSEBITS 128
#define SSEMODE 10
#else
#define SSEBITS 64
#define SSEMODE 0
#endif

#ifndef WIN32
// #define FMA_AVAILABLE __FMA__
#endif


#if __GNUC__ > 4 ||							\
  (__GNUC__ == 4 && (__GNUC_MINOR__ > 9 ||				\
		     (__GNUC_MINOR__ == 9 &&  __GNUC_PATCHLEVEL__ >= 1)))
//#define OpenMP4 1
#endif


//#define ALIGNED __declspec(align(SSEBITS/8))


// ??   __LP64__


/* 
until version 0.1.x there was a distinction between 32 Bit machines and 64 bit machines; 64 bit machines could be run as 32 Bit machine, but with a loss of about 10% of speed

 ...


*/

#ifdef AVX512
//#include <immintrin.h>
#endif


#if defined AVX || defined AVX2
#include <x86intrin.h>
#else 

#ifdef MMX
//#include <mmintrin.h>
#endif

#ifdef SSE
//#include <xmmintrin.h>
#endif

#ifdef SSE2
  #include <emmintrin.h>
#endif

#ifdef SSE3
//#include <pmmintrin.h>
#endif

#ifdef SSSE3
//#include <tmmintrin.h>
#endif

#ifdef SSE4A
//#include <ammintrin.h>
#endif

#ifdef SSE412
//#include <smmintrin.h>
#endif

#endif


#if defined AVX
#define BytesPerBlock 32
#define BlockType __m256i ALIGNED
#define Double __m256d
#define MAXDOUBLE _mm256_max_pd
#define LOAD _mm256_load_si256
#define EXPDOUBLE _mm256_exp_pd
#define ADDDOUBLE  _mm256_add_pd
#define MULTDOUBLE _mm256_mul_pd 
#define LOADuDOUBLE _mm256_loadu_pd
#define LOADDOUBLE _mm256_load_pd
#define STOREuDOUBLE _mm256_storeu_pd
#define ZERODOUBLE _mm256_setzero_pd()
#elif defined SSE2
#define BytesPerBlock 16
#define BlockType __m128i ALIGNED
#define Double __m128d
#define MAXDOUBLE _mm_max_pd
#define LOAD _mm_load_si128
#define EXPDOUBLE _mm_exp_pd
#define ADDDOUBLE  _mm_add_pd
#define MULTDOUBLE _mm_mul_pd 
#define LOADuDOUBLE _mm_loadu_pd
#define LOADDOUBLE _mm_load_pd
#define STOREuDOUBLE _mm_storeu_pd
#define ZERODOUBLE _mm_setzero_pd()
#else
#define BytesPerBlock 8
#endif
#define algn_general(X)  ((1L + (uintptr_t) (((uintptr_t) X - 1L) / BytesPerBlock)) * BytesPerBlock)
double inline *algn(double *X) { return (double *) algn_general(X); }
#define ALIGNED __attribute__ ((aligned (BytesPerBlock)))
#define doubles (BytesPerBlock / 8)





#endif


