
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

#ifndef RFMachines_H
#define RFMachines_H 1

#define LASTi(sys) (sys).last
#define CUMXMITi(sys) (sys).cumxdim
#define XDIMi(sys) (sys).xdim // OK
#define LOGDIMi(sys) (sys).logicaldim
#define MAXDIMi(sys) (sys).maxdim
#define ISOi(sys) (sys).iso
#define DOMi(sys) (sys).dom
#define TYPEi(sys) (sys).type
#define NRi(sys) (sys).nr
#define MODELNR(cov) SYSMODEL(SYSOF(cov))
#define isSetLastSystem(sys) (LASTi((sys)[0]) != UNSET)

#define ANYDIM ANYDIMOF(cov)
#define SYSTEMS(s) (LASTSYSTEM(s) + 1)


#define CALLING SYSOF(cov->calling)
#define CALLINGNR SYSMODEL(CALLING)

#define PREV PREVSYSOF(cov)
#define PREVCUMXOHNE(s) CUMXOHNE(PREV, s)
#define PREVCUMXMIT(s) CUMXMIT(PREV, s)
#define PREVXDIM(s) XDIM(PREV, s)
#define PREVISO(s) ISO(PREV, s)
#define PREV_INITIALISEDXX isSetLastSystem(PREVSYSOF(cov)) //
#define PREV_INITIALISED (isSetLastSystem(PREVSYSOF(cov)) && XDIM(PREV, 0) != UNSET) // siehe settrafo in op...cc
#define CONDPREVISO(s) (PREV_INITIALISED ? PREVISO(s) : ISO_MISMATCH)
#define CONDPREVDOM(s) (PREV_INITIALISED ? PREVDOM(s) : DOMAIN_MISMATCH)
#define PREVDOM(s) DOM(PREV, s)
#define PREVTYPE(s) SYSTYPE(PREV, s)
#define PREVLOGDIM(s) LOGDIM(PREV, s)
#define PREVTOTALXDIM TOTALXDIM(PREV)
#define PREVLASTSYSTEM LASTSYSTEM(PREV)
#define PREVSYSTEMS SYSTEMS(PREV)
#define PREVRESETLAST RESETLAST(PREV)
#define PREVLASTSYS(s) LASTSYS(PREV, s)
#define TRAFONR SYSMODEL(PREV)

#define GATTER GATTERSYSOF(cov)
#define GATTERCUMXOHNE(s) CUMXOHNE(GATTER, s)
#define GATTERCUMXMIT(s) CUMXMIT(GATTER, s)
#define GATTERXDIM(s) XDIM(GATTER, s)
#define GATTERISO(s) ISO(GATTER, s)
#define GATTERDOM(s) DOM(GATTER, s)
#define GATTERTYPE(s) SYSTYPE(GATTER, s)
#define GATTERLOGDIM(s) LOGDIM(GATTER, s)
#define GATTERTOTALXDIM TOTALXDIM(GATTER)
#define GATTERLASTSYSTEM LASTSYSTEM(GATTER)
#define GATTERSYSTEMS SYSTEMS(GATTER)
#define GATTERNR SYSMODEL(GATTER)

#define OWN SYSOF(cov) 
#define OWNCUMXOHNE(s) CUMXOHNE(OWN, s)
#define OWNCUMXMIT(s) CUMXMIT(OWN, s)
#define OWNXDIM(s) XDIM(OWN, s)
#define OWNISO(s) ISO(OWN, s)
#define OWNDOM(s) DOM(OWN, s)
#define CONDOWNISO(s) (PREV_INITIALISED ? OWNISO(s) : ISO_MISMATCH)
#define CONDOWNDOM(s) (PREV_INITIALISED ? OWNDOM(s) : DOMAIN_MISMATCH)
#define OWNTYPE(s) SYSTYPE(OWN, s)
#define OWNLOGDIM(s) LOGDIM(OWN, s)
#define OWNMAXDIM(s) MAXDIM(OWN, s)
#define OWNTOTALXDIM TOTALXDIM(OWN)
#define OWNLASTSYSTEM LASTSYSTEM(OWN)
#define OWNSYSTEMS SYSTEMS(OWN)
#define OWNRESETLAST RESETLAST(OWN)
#define COVNR (SYSMODEL(OWN))
#define VDIM0 cov->vdim[0]
#define VDIM1 cov->vdim[1]

#define SUB SYSOF(sub)
#define SUBCUMXOHNE(s) CUMXOHNE(SUB, s)
#define SUBCUMXMIT(s) CUMXMIT(SUB, s)
#define SUBXDIM(s) XDIM(SUB, s)
#define SUBISO(s) ISO(SUB, s)
#define SUBDOM(s) DOM(SUB, s)
#define SUBTYPE(s) SYSTYPE(SUB, s)
#define SUBLOGDIM(s) LOGDIM(SUB, s)
#define SUBMAXDIM(s) MAXDIM(SUB, s)
#define SUBTOTALXDIM TOTALXDIM(SUB)
#define SUBLASTSYSTEM LASTSYSTEM(SUB)
#define SUBSYSTEMS SYSTEMS(SUB)
#define SUBNR SYSMODEL(SUB)

#define NEXT SYSOF(next)
#define NEXTCUMXOHNE(s) CUMXOHNE(NEXT, s)
#define NEXTCUMXMIT(s) CUMXMIT(NEXT, s)
#define NEXTXDIM(s) XDIM(NEXT, s)
#define NEXTISO(s) ISO(NEXT, s)
#define NEXTDOM(s) DOM(NEXT, s)
#define NEXTTYPE(s) SYSTYPE(NEXT, s)
#define NEXTLOGDIM(s) LOGDIM(NEXT, s)
#define NEXTMAXDIM(s) MAXDIM(NEXT, s)
#define NEXTTOTALXDIM TOTALXDIM(NEXT)
#define NEXTLASTSYSTEM LASTSYSTEM(NEXT)
#define NEXTSYSTEMS SYSTEMS(NEXT)
#define NEXTNR SYSMODEL(NEXT)


#define DEFSYS(Cov) DefList[MODELNR(Cov)].systems[(Cov)->variant == UNSET \
						  ? 0 : (Cov)->variant]
#define DEF DEFSYS(cov)
#define DEFISO(s) ISO(DEF, s)
#define DEFDOM(s) DOM(DEF, s)
#define DEFTYPE(s) SYSTYPE(DEF, s)
#define DEFMAXDIM(s) MAXDIM(DEF, s)
#define DEFSYSTEMS SYSTEMS(DEF)
#define DEFLASTSYSTEM LASTSYSTEM(DEF)


#define SET_OUT_OF(D) (D)->lpx[GLOBAL.general.set]
#define NROW_OUT_OF(D) (D)->nrow[GLOBAL.general.set]
#define NCOL_OUT_OF(D) (D)->ncol[GLOBAL.general.set]

#define P(IDX) PARAM(cov, IDX) 
#define PINT(IDX) PARAMINT(cov, IDX)
#define PCHAR(IDX) PARAMCHAR(cov, IDX)
#define P0(IDX) PARAM0(cov, IDX) 
#define P0INT(IDX) PARAM0INT(cov, IDX)
#define P0CHAR(IDX) PARAM0CHAR(cov, IDX)
#define PVEC(IDX) PARAMVEC(cov, IDX) 
#define PENV(IDX) PARAMENV(cov, IDX) 
#define PLANG(IDX) PARAMLANG(cov, IDX) 
#define PLIST(IDX) PARAMLIST(cov, IDX)
#define PARAMSEXP(Cov, IDX) ((sexp_type *) (Cov)->px[IDX]) /* kein assert! */
#define PSEXP(IDX) PARAMSEXP(cov, IDX)
#define NROW(IDX) cov->nrow[IDX]
#define NCOL(IDX) cov->ncol[IDX]
#define LP(IDX) LPARAM(cov, IDX) 
#define LPINT(IDX) LPARAMINT(cov, IDX)
#define LP0(IDX) LPARAM0(cov, IDX) 
#define LP0INT(IDX) LPARAM0INT(cov, IDX)
#define LNRC_(I, rc) PLIST(I)->rc[SET_IDX(cov, I)] /* see also Machine*.h */
#define LNROW(IDX) LNRC(IDX, nrow)
#define LNCOL(IDX) LNRC(IDX, ncol)
#define LPARAMNROW(Cov, IDX) PARAMLIST(Cov, IDX)->nrow[SET_IDX(Cov, IDX)]
#define LPARAMNCOL(Cov, IDX) PARAMLIST(Cov, IDX)->ncol[SET_IDX(Cov, IDX)]

//#define KAPPASUBisNULL(COV, IDX) COV->kappasub[IDX] == NULL
//#define KSUBisNULL(IDX) KAPPASUBisNULL(cov, IDX)


#define PARAMFREE(Cov, IDX) if ((Cov)->px[IDX] == NULL) { } else {	\
  if (DefList[MODELNR(Cov)].kappatype[IDX] < LISTOF) {		\
    UNCONDFREE((Cov)->px[IDX]);				\
  } else {						\
    LIST_DELETE( (listoftype **) ((Cov)->px + IDX));	\
  }							\
  (Cov)->ncol[IDX] = (Cov)->nrow[IDX] = 0;			\
  }


#define PARAMALLOC(Cov, IDX, ROW, COL) {					\
    int _PARAMsize;							\
    switch(DefList[MODELNR(Cov)].kappatype[IDX]) {			\
    case REALSXP : _PARAMsize = sizeof(double); break;			\
    case INTSXP : _PARAMsize = sizeof(int); break;			\
    default : 						\
      if ((Cov)->kappasub[IDX]!=NULL && MODELNR((Cov)->kappasub[IDX])==DISTRIBUTION) { \
	  ERR("argument value recognized as distribution family although it should not. Maybe the error is caused by a non-existing variable."); \
      } else BUG;							\
    }									\
    assert((Cov)->px[IDX]==NULL);					\
    (Cov)->nrow[IDX] = ROW; (Cov)->ncol[IDX] = COL;			\
    if (((Cov)->px[IDX] =						\
	 (double*) CALLOC((ROW) * (COL), _PARAMsize)) == NULL) {	\
      XERR(ERRORMEMORYALLOCATION)					\
    }									\
  }
#define PALLOC(IDX, ROW, COL) PARAMALLOC(cov, IDX, ROW, COL)

#define PFREE(IDX) PARAMFREE(cov, IDX)
#define PARAMisNULL(Cov, IDX) ((Cov)->px[IDX] == NULL)
#define PisNULL(IDX) PARAMisNULL(cov, IDX)
#define PARAMtoNULL(Cov, IDX) (Cov)->px[IDX] = NULL
#define PtoNULL(IDX) PARAMtoNULL(cov, IDX)

#define QALLOC(nr) QcovALLOC(cov, nr)
#define QALLOC1(x1) if (cov->q == NULL) { QALLOC(1); cov->q[0] = x1; }
#define QALLOC2(x1,x2) if (cov->q == NULL) { QALLOC(2); cov->q[0] = x1; cov->q[1] = x2; }
#define QALLOC3(x1,x2,x3) if (cov->q == NULL) { QALLOC(3); cov->q[0] = x1; cov->q[1] = x2; cov->q[2] = x3;}  

#define QFREE					\
  if ((cov)->q == NULL) { } else {		\
    UNCONDFREE((cov)->q);			\
    (cov)->qlen = 0;				\
  }



#define set_xdim(S, I, V) DEBUG(set_xdim_intern(S, I, V))
#define set_xdim_blank(S, I, V) XDIMi((S)[I])=V // OK
#define set_logdim(S, I, V) LOGDIMi((S)[I])=V
#define set_nr(S, V) NRi((S)[0])=V
#define set_trafo(V) set_nr(PREV, V)
#define set_nri(S, I, V) NRi((S)[I])=V
#define set_maxdim(S,I,V) MAXDIMi((S)[I])=V
#define set_last(S,I,V) LASTi((S)[I])=V
#define set_iso(S,I,V) ISOi((S)[I])=V
#define set_dom(S,I,V) DOMi((S)[I])=V
#define set_type(S,I,V) TYPEi((S)[I])=V
#define set_cumxmit(S,I,V) CUMXMITi((S)[I])=V


void set_system_type(system_type *sys, Types type);
void set_system_domain(system_type *sys, domain_type dom);
void set_system(system_type * sys, int idx, int logicaldim, int maxdim, 
		int xdim, Types type, domain_type dom, isotropy_type iso);
void set_system(system_type * sys, int idx, int logicaldim, int maxdim, 
		int xdim, Types type, domain_type dom, isotropy_type iso, 
		bool check_unset);

void set_both_systems(model *cov, int dim, Types type);

sortsofparam SortOf(model *cov,int k, int row, int col, sort_origin origin);
int total_logicaldim(system_type * sys);
void set_xdim_intern(system_type *sys, int i, int value);



#define LocLoc(loc) ((loc) != NULL ? (loc)[GLOBAL.general.set % (loc)[0]->len]\
		     : NULL)
#define Loc(cov) ((cov)->ownloc!=NULL ? LocLoc((cov)->ownloc) : \
		  LocLoc((cov)->prevloc))
#define PLoc(cov) ((cov)->ownloc != NULL ? (cov)->ownloc : (cov)->prevloc)
#define PrevLoc(cov)						\
  ((cov)->prevloc == NULL ? NULL :				\
   (cov)->prevloc[GLOBAL.general.set % (cov)->prevloc[0]->len])
#define OwnLoc(cov) ((cov)->ownloc == NULL ? NULL : \
		     (cov)->ownloc[GLOBAL.general.set % (cov)->ownloc[0]->len])
#define Getgrid(cov) (PLoc(cov) != NULL ? Loc(cov)->grid : false)
#define GetTime(cov) (PLoc(cov) != NULL ? Loc(cov)->Time : false)
#define GetLoctsdim(cov) (PLoc(cov) != NULL ? Loc(cov)->timespacedim : 0)
#define GetLocspatialdim(cov) (PLoc(cov) != NULL ? Loc(cov)->spatialdim : 0)
#define Gettotalpoints(cov) (PLoc(cov) != NULL ? Loc(cov)->totalpoints : 0)
#define Getspatialpoints(cov) (PLoc(cov) != NULL ? Loc(cov)->spatialtotalpoints\
			       : 0)
#define Getcaniso(cov) (PLoc(cov) != NULL ? Loc(cov)->caniso : NULL)
#define Getxgr(cov) (PLoc(cov) != NULL ? Loc(cov)->xgr : NULL)
#define Getx(cov) (PLoc(cov) != NULL ? Loc(cov)->x : NULL)
#define DistancesGiven(cov) (PLoc(cov) != NULL ? Loc(cov)->distances : false)
#define GET_LOC_SETS(cov) (PLoc(cov) != NULL ? PLoc(cov)[0]->len : 0)





#define DEBUG(A) A

//////////////////////////////////////////////////////////////////////
// CHECKING 
//////////////////////////////////////////////////////////////////////
#define ASSERT_GATTER(Cov) assert(TrafoOK(Cov, __FILE__, __LINE__))
#define ASSERT_CHECKED(Cov) assert((Cov)->checked)
#define ASSERT_LOC_GIVEN \
  if (loc != NULL) { } else {PMI(cov); LOC_NOT_INITIALISED;}
#define ASSERT_CARTESIAN  if (isCartesian(OWN)) {} else RETURN_ERR(ERRORCARTESIAN)
#define ASSERT_QUASIONESYSTEM if (QuasiOneSystem(cov)) { } else BUG;
#define ASSERT_ONESYSTEM {						\
    if (OWNLASTSYSTEM != 0 && (OWNLASTSYSTEM != 1 ||			\
			       !equalsIsotropic(OWNISO(0))		\
			       /*  ||  !equalsIsotropic(OWNISO(1)) */ )) BUG;}
#define ASSERT_UNREDUCED {assert(OWNXDIM(0) == OWNLOGDIM(0) && PREVXDIM(0) == OWNLOGDIM(0));}
#define ASSERT_GAUSS_METHOD(METHOD) DEBUGINFO;				\
  if(!hasGaussMethodFrame(cov) || cov->method != (METHOD))	{ \
    SERR4("Gaussian field for '%.50s' only possible with '%.50s' as method. Got frame '%.50s' and method '%.50s'.", \
	 NICK(cov),							\
	  DefList[(METHOD) == RandomCoin ? RANDOMCOIN_USER	:	\
		 gaussmethod[METHOD] - 				\
		 DefList[gaussmethod[METHOD]].internal].nick,		\
	  TYPE_NAMES[cov->frame],					\
	  gaussmethod[cov->method] <= 0 ? "MISMATCH"	:\
	  DefList[cov->method == RandomCoin ? RANDOMCOIN_USER	:	\
		  gaussmethod[cov->method] -				\
		  DefList[gaussmethod[cov->method]].internal].nick	\
	  ) }



#endif
