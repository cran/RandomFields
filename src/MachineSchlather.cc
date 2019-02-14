
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives (spectral measures, 
 tbm operators)

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2017 Martin Schlather
 

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

#ifdef XXSCHLATHERS_MACHINE
#include "RF.h"

void inline COPYALLSYSTEMS(system_type *to, system_type *from, bool keepnr) {
  assert_sys(to); assert_sys(from);					
  int nr_ = MISMATCH;						       
  if (keepnr) nr_= SYSMODEL(to);					
  MEMCOPY(to, from, sizeof(Systems_type));			       
  if (keepnr) set_nr(to, nr_);
}    

#endif
