/* 
 Authors
 Martin Schlather, martin.schlather@cu.lu

 all around the nugget effect -- needs special treatment 

 Copyright (C) 2001 -- 2004 Martin Schlather, 

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


#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "RFsimu.h"



int init_hyperplane(key_type *key, int m){ return ERRORNOTPROGRAMMED; }
void do_hyperplane(key_type *key, int m, Real *res ) {
  assert(false);
 {
    register long i,endfor;
    endfor = key->totalpoints;
    for (i=0;i<endfor;i++) { res[i]=0.0; }
  }
}
int init_special(key_type *key, int m){ 
  return ERRORNOTPROGRAMMED; 
   //  calls key->initother

  // do_special calls key->other      
  /* 
     if (m==0) KEY[*keyNr].cov->other(key,res ); 
      else { 
	KEY[*keyNr].cov->other(key,part_result ); 
	for (i=0, i<KEY[*keyNr].totalpoints; i++) res[i] += part_result[i];
      } 
  */
}

void do_special(key_type *key, int m, Real *res ){
  assert(false);
 {
    register long i,endfor;
    endfor = key->totalpoints;
    for (i=0;i<endfor;i++) { res[i]=0.0; }
  }
}

 
