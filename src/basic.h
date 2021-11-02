

/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 -- 2017 Martin Schlather

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





#ifndef GSL_VS_R_H
#define GSL_VS_R_H 1

#ifdef SCHLATHERS_MACHINE
#undef SCHLATHERS_MACHINE
#endif

#ifdef RANDOMFIELDS_DEBUGGING
#undef RANDOMFIELDS_DEBUGGING
#endif

#ifdef showfree
#undef showfree
#endif

#ifdef DOPRINT
#undef DOPRINT
#endif



#define MAXNRCOVFCTS 300
#define MAXUNITSCHAR 10
#define MAXINVERSIONS 2
#define DISTMAXSTEPS 1000


#ifdef __GNUC__
#define HIDE_UNUSED_VARIABLE 1
#endif

#if __GNUC__ >= 7
#define FALLTHROUGH_OK __attribute__ ((fallthrough));
#else
#define FALLTHROUGH_OK   
#endif      


#endif /* GSL_VS_R_H */
