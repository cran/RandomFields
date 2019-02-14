
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2017 -- 2017 Martin Schlather

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





#ifndef RF_MACHINE_OTHERS
#define RF_MACHINE_OTHERS 1

#include "Machine.h"

#define LOC_NOT_INITIALISED SERR("locations not initialised.")
#define LNRC(I, rc) LNRC_(I, rc)
#define STOPAFTER(COND, DO)

#define PREVSYSOF(cov) (cov)->prev
#define GATTERSYSOF(cov) (cov)->gatter
#define SYSOF(cov) (cov)->own


#define LASTSYSTEM(sys) LASTi((sys)[0])
#define RESETLAST(sys) LASTi((sys)[0])=0
#define CUMXOHNE(sys, s) ((s==0) ? 0 : CUMXMITi((sys)[(s)-1]))
#define CUMXMIT(sys, s) CUMXMITi((sys)[s])
#define TOTALXDIM(sys) CUMXMITi((sys)[LASTSYSTEM(sys)])
#define XDIM(sys, s) XDIMi((sys)[s]) // OK
#define LOGDIM(sys, s) LOGDIMi((sys)[s])
#define MAXDIM(sys, s) MAXDIMi((sys)[s])


#define ISO(sys, s) ISOi((sys)[s])
#define DOM(sys, s) DOMi((sys)[s])
#define SYSTYPE(sys, s) TYPEi((sys)[s])
#define COPYSYS(to, from) MEMCOPY(to, from, sizeof(system_type))
#define COPYALLSYSTEMS(to, from, keepnr) {		\
    int nr_;						\
    if (keepnr) nr_ = SYSMODEL(to);			\
    MEMCOPY(to, from, sizeof(Systems_type));		\
    if (keepnr) { set_nr(to, nr_); }			\
    /* muss hier mit SET_NR statt set_nr gearbeitet werden ?? */	\
  }
#define COPYALLSYSTEMS_COND(to, from, keepnr) { 			\
    bool gatter_set_ = LASTi(((cov)->gatter)[0]) >= 0;			\
    int nr_ = MISMATCH; if (gatter_set_ && keepnr) nr_= SYSMODEL(to);	\
    MEMCOPY(to, from, sizeof(Systems_type));				\
    if (gatter_set_ && keepnr) { set_nr(to, nr_); }			\
    /* muss hier mit SET_NR statt set_nr gearbeitet werden ?? */	\
  }
#define ANYOWNDIM OWNTOTALXDIM
#define ANYDIMOF(cov) OWNTOTALXDIM
#define SYSMODEL(sys) NRi((sys)[0])



#define SET_IDX(Cov, IDX) (GLOBAL.general.set % (Cov)->nrow[IDX])

#define CHECK(C,L,X,type,D,I,V,R) check2X(C,L,X,type,D,I,V,R)
#define CHECK_THROUGHOUT(C,P,type,D,I,V,R) check2Xthroughout(C,P,type,D,I,V,R)
#define CHECK_NO_TRAFO(C,T,X,type,D,I,V,R) check2Xnotrafo(C,T,X,type,D,I,V,R)


#define STRUCT(Cov, NM) DefList[FIRSTGATTER].Struct(Cov, NM)

#define PARAM(Cov, IDX) ((double *) (Cov)->px[IDX])
#define PARAMINT(Cov, IDX) ((int *) (Cov)->px[IDX])
#define PARAMCHAR(Cov, IDX) ((char **) (Cov)->px[IDX])
#define PARAMVEC(Cov, IDX) ((sexp_type *) (Cov)->px[IDX])->sexp
#define PARAMENV(Cov, IDX) ((sexp_type *) (Cov)->px[IDX])
#define PARAMLANG(Cov, IDX) ((sexp_type *) (Cov)->px[IDX])
#define PARAMLIST(Cov, IDX) ((listoftype *) (Cov)->px[IDX])
#define LPARAM(Cov, IDX) ((double *) (PARAMLIST(Cov, IDX)->lpx[SET_IDX(Cov, IDX)]))
#define LPARAMINT(Cov, IDX) ((int *) (PARAMLIST(Cov, IDX)->lpx[SET_IDX(Cov, IDX)]))
//#define LISTLIST(Cov, IDX) ((listoftype *) (Cov)->px[IDX])

#define PARAM0(Cov, IDX) PARAM(Cov, IDX)[0]
#define PARAM0INT(Cov, IDX) PARAMINT(Cov, IDX)[0]
#define LPARAM0(Cov, IDX) LPARAM(Cov, IDX)[0]
#define LPARAM0INT(Cov, IDX) LPARAMINT(Cov, IDX)[0]

#define PCOPY(TO, FROM, IDX) 						\
  MEMCOPY((TO)->px[IDX], (FROM)->px[IDX],				\
	  ((FROM)->nrow[IDX]) * ((FROM)->ncol[IDX]) *			\
	  (DefList[MODELNR(FROM)].kappatype[IDX]==REALSXP ? sizeof(double) : \
	   DefList[MODELNR(FROM)].kappatype[IDX]==INTSXP ? sizeof(int) : -1))

#define QcovALLOC(cov, nr) {						\
    (cov)->qlen = (int) (nr);						\
    assert((cov)->q == NULL);						\
    if (((cov)->q = (double*) CALLOC((int) (nr), sizeof(double))) == NULL) \
      RFERROR("memory allocation error for local memory");		\
}


#define B2D(X) (double) (X)
#define B2I(X) (int) (X)
#define E2D(X) (X)
#define E2I(X) (int) (X)
#define E2B(X) (bool) (X)
#define I2D(X) (double) (X)

#define RootOf(cov) ((cov)->root)
#define KEYtypeOf(cov)((cov)->base)
#define KT_NAOK_RANGE KT->naok_range
#define setKT_NAOK_RANGE(KT,X) KT->naok_range = X
#define NAOK_RANGE KEYtypeOf(cov)->naok_range
#define set_NAOK_RANGE(X) KEYtypeOf(cov)->naok_range = X

#endif
