/*
 Authors 
 Felix Ballani

 Copyright (C) 2014 -- 2017 Felix Ballani

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


#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "RF.h"
#include "cubicsolver.h"

#define ERRORNOTCUBIC 501

//Constants for the cubic equation solver
// the function solve the equation ax^3 + bx^2 + cx + d = 0. It returns the smallest positive real root or real root
// returns 0 if the equation is trully cubic.
int cubicsolver(double a, double b, double c, double d, double roots[][2]) {

    //the equation is quadratic, but not cubic. The first coefficient must be non 0.
    if (a == 0) FAILED1("a=%10e NOT OK", a);

    b /= a;
    c /= a;
    d /= a;
    double disc, q, r, dum1, s, t, term1, r13;

    q = (3.0*c - (b*b))/9.0;
    r = -(27.0*d) + b*(9.0*c - 2.0*(b*b));
    r /= 54.0;
    disc = q*q*q + r*r;
    roots[0][1] = 0; //The first root is always real.
    term1 = (b/3.0);
    if (disc > 0) { // one root real, two are complex
        s = r + SQRT(disc);
        if (s < 0) s = -POW(-s, (1.0/3.0)); else s = POW(s, (1.0/3.0));
	t = r - SQRT(disc);
        if (t < 0) t = -POW(-t, (1.0/3.0)); else t = POW(t, (1.0/3.0));
        roots[0][0] = -term1 + s + t;
        term1 += (s + t)/2.0;
        roots[2][0] = -term1;
        roots[1][0] = -term1;
        term1 = SQRT(3.0)*(-t + s)/2;
        roots[1][1] = term1;
        roots[2][1] = -term1;
        return NOERROR;
    }
    // End if (disc > 0)
    // The remaining options are all real
    roots[1][1] = 0;
    roots[2][1] = 0;
    if (disc == 0){ // All roots real, at least two are equal.      
      if (r < 0) r13 = -POW(-r,(1.0/3.0));
      else {r13 = POW(r,(1.0/3.0));}
      roots[0][0] = -term1 + 2.0*r13;
      roots[1][0] = -(r13 + term1);
      roots[2][0] = -(r13 + term1);
      return 0;
    } // End if (disc == 0)
    // Only option left is that all roots are real and unequal (to get here, q < 0)
    q = -q;
    dum1 = q*q*q;
    dum1 = ACOS(r/SQRT(dum1));
    r13 = 2.0*SQRT(q);
    roots[0][0]= -term1 + r13*COS(dum1/3.0);
    roots[1][0] = -term1 + r13*COS((dum1 + 2.0*M_PI)/3.0);
    roots[2][0] = -term1 + r13*COS((dum1 + 4.0*M_PI)/3.0);
    return NOERROR;
}
