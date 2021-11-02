/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2021 -- 2021 Martin Schlather

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

#ifndef miraculix_IntrinsicsBase_H
#define miraculix_IntrinsicsBase_H 1

#ifdef _OPENMP
  #ifdef SCHLATHERS_MACHINE
  #else
    #define DO_PARALLEL 1
  #endif
#else
  #ifdef DO_PARALLEL
    #undef DO_PARALLEL
  #endif
#endif

//#ifdef WIN32
//#ifdef DO_PARALLEL
//#undef DO_PARALLEL // make a comment to get parallel (part 1, see also part 2)
//#endif
//#endif

#if defined AVX512
#undef AVX512
#endif

#if defined __AVX2__
#define AVX2 1
#else
#ifdef AVX2
#undef AVX2
#endif
#endif

#if defined __AVX__ 
#define AVX 1
#else
#ifdef AVX2
#undef AVX2
#endif
#endif

#if defined __SSE41__ 
#define SSE41 1
#else
#ifdef SSE41
#undef SSE41
#endif
#endif

#if defined  __SSSE3__  // needs SSE2 as well
#define SSSE3 1
#else
#ifdef SSSE3
#undef SSSE3
#endif
#endif

#if defined  __SSE3__
#define SSE3 1
#else
#ifdef SSE3
#undef SSE3
#endif
#endif

#if defined  __SSE2__ 
#define SSE2 1 
#else
#ifdef SSE2
#undef SSE2
#endif
#endif



#endif
