/*
 Authors 
 Christoph Berreth, cberreth@mail.uni-mannheim.de

 library for simulation of random fields on spheres 
 -- init part and error messages
 
 Copyright (C) 2014 -- 2014 Christoph Berreth
               2015 -- 2017 Martin Schlather

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

/* header files */
#include "RF.h"    /* for struct cov_model, macro P0,
		    enthält u.A. basic.h für vordefinierte Konstanten*/

#include "SpherModels.h"   /* muss nach RF.h stehen, da cov_model 
			      enthalten. */



/* function to initialize the sperical models */
void  SpherModelsinit( ){
   

  /* ------ sinepower model ------ */
 
  pref_type
    psinepower = {0, 0, 0,  0, 0, 5, 0, 0, 0, 0,  0, 0,  0,  5};
  //             CE CO CI TBM Sp di sq Ma av  n mpp Hy spf any
  IncludePrim("sinepower", PosDefType, 1, NULL,
	      XONLY, SPHERICAL_ISOTROPIC, checkOK, rangeSinePower,
	      psinepower, SCALAR, 2, false, MONOTONE);
  kappanames("alpha", REALSXP);
  addCov(SinePower, NULL, NULL);


  /* ------ model of the multiquadratic family ------ */
 
  pref_type
    pmultiquad = {0, 0, 0,  0, 0, 5, 0, 0, 0, 0,  0, 0,  0,  5};
  //        CE CO CI TBM Sp di sq Ma av  n mpp Hy spf any
  IncludePrim("multiquad", PosDefType, 2, NULL,
	      XONLY, SPHERICAL_ISOTROPIC, checkOK, rangeMultiquad,
	      pmultiquad, SCALAR, 2, false, MONOTONE);
  kappanames("delta", REALSXP, "tau", REALSXP);
  addCov(Multiquad, NULL, NULL);


 

  /* ------ Choquet's representation ------ */
  /*
  pref_type
    pchoquet = {0, 0, 0,  0, 0, 5, 0, 0, 0, 0,  0, 0,  0,  5};
  //            CE CO CI TBM Sp di sq Ma av  n mpp Hy spf any
  IncludeModel("choquet", PosDefType, 0, 0, 1, kappa_choquet,
	      XONLY, ISOTROPIC, checkChoquet, rangeChoquet,
	       pchoquet, false, 1, INFDIM, true, MONOTONE);
              // IncludeModel: see startGetNset.cc
  kappanames("b", REALSXP);
  addCov(Choquet, NULL, NULL);
  //add_paramtype(paramtype_choquet);
  */

 

}


