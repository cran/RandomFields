/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- 2017 Martin Schlather

This program is freeall software; you can redistribute it and/or
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

 
#include "def.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "extern.h"
#include "questions.h"
#include "operator.h"
#include "Processes.h"
#include "startGetNset.h"

//#define xdebug 1
//#define ydebug 1

//  printf("L=%d E=%d cov:err=%d err_lev=%d %.50s\n", L, X, cov->err,  cov->err_level, NAME(cov));

//int Zaehler = 0;

#define LEVEL_UPDATE(L,X)			\
  if ((X) == NOERROR) {				\
    cov->err_level=0;				\
    cov->err = NOERROR;				\
  } else if (L > cov->err_level) {		\
    cov->err_level=L;				\
    cov->err = X;				\
  }				       

#define NOERROR_RETURN(L) {	\
  if (L > cov->err_level) {	\
    cov->err_level=L;		\
    cov->err = NOERROR;		\
  }				\
  return NOERROR; }

#define GOTO(X) { err = X; goto ErrorHandling; }
// #define GOTO(X) {printf("error re-turned at line %d\n", __LINE__); err = X; goto ErrorHandling;  }
#ifdef xdebug
#define RETURN(L, X) { LEVEL_UPDATE(L,X);  if (cov->err != NOERROR) {PRINTF("error %d L=%d in '%.50s'.\n", cov->err, L, NAME(cov));} return cov->err; }
#define ERR_RETURN(L, X) { FERR(X); LEVEL_UPDATE(L,ERRORM);PRINTF("error in '%.50s': %.50s\n", NAME(cov), ERRORM); return cov->err; }
#define ERR_RETURN1(L,X,Y) { FERR1(X,Y); LEVEL_UPDATE(L,ERRORM); PRINTF("error in '%.50s': %.50s\n", NAME(cov), ERRORM);return cov->err;}
#define ERR_RETURN2(L,X,Y,Z) {FERR2(X,Y,Z); LEVEL_UPDATE(L,ERRORM); PRINTF("error in '%.50s': %.50s\n", NAME(cov), ERRORM); return cov->err;}
#else 
#define RETURN(L, X) { LEVEL_UPDATE(L,X);  return cov->err; }
#define ERR_RETURN(L, X) { FERR(X); LEVEL_UPDATE(L,ERRORM); return cov->err; }
#define ERR_RETURN1(L,X,Y) { FERR1(X,Y); LEVEL_UPDATE(L,ERRORM); return cov->err;}
#define ERR_RETURN2(L,X,Y,Z) {FERR2(X,Y,Z); LEVEL_UPDATE(L,ERRORM); return cov->err;}
#endif


#define LEVEL_DOM 40
#define LEVEL_CHECK 43
#define LEVEL_LAST 48
#define CONT(L) { levels[dom] = L; continue;} 
#define LEVEL(L, X) {						\
    errs[dom] = X;						\
    levels[dom] = L;							\
    errorMSG(errs[dom], cov->err_msg, KT, ERRSTR[dom]);			\
    if (PL > PL_ERRORS) { PRINTF("check level %d failed for %.50s (%d), %.50s %.50s.\n", levels[dom], NAME(cov), errs[dom], KT->error_loc, ERRSTR[dom]);} \
    continue;								\
  }

//if (false) { printf("check level %d failed for %.50s: %.50s (%.50s)\n", levels[dom], NAME(cov), cov->err_msg, KT->error_loc);}


#define ShortList 5
#define ShortN 8
char kurz_[ShortList][ShortN + 1] = {""};
const char *Short(int i, const char * x) {
  assert(i < ShortList);
  strcopyN(kurz_[i], x, ShortN);
  kurz_[i][ShortN] = '\0';
  return kurz_[i];
}


int SetXdimLogdim(model *cov, isotropy_type *Iso, int lastsystem);
int set_own_domain_and_check(model *cov, int vdim0, int vdim1,
			     ext_bool coord_trafo);				
int change_coord_system(model *cov, ext_bool coordinate_trafo,
			errorstring_type ERRSTR);
int SetGatterNr(model *cov, errorstring_type ERRSTR);
int checkDims(model *cov, int vdim0, int vdim1, errorstring_type ERRSTR);
int check_within_range(model *cov, bool NAOK, errorstring_type ERRSTR);
int checkkappas(model *cov, bool errornull, errorstring_type ERRSTR);

int check2passtype(model *cov, system_type* s, Types type,
		   int vdim0, int vdim1, Types frame){
  assert(LASTi(s[0]) >= 0);
  COPYALLSYSTEMS(PREV, s, false);
  ASSERT_ONESYSTEM;
  set_type(PREV, 0, type);
  int vdim = DefList[SYSMODEL(s)].vdim;
  return check2X(cov, vdim == PARAM_DEP ? vdim0 : vdim,
		    vdim == PARAM_DEP ? vdim1 : vdim, 
		 frame, false);
}


int check2passframe(model *cov, system_type* s, int vdim0, int vdim1, 
		 Types frame){
  assert(LASTi(s[0]) >= 0);
  COPYALLSYSTEMS(PREV, s, false);
  return check2X(cov, vdim0, vdim1, frame, false);
} 

int check2passTF(model *cov, system_type* s, Types type, int vdim,
		 Types frame) {
  assert(LASTi(s[0]) >= 0);
  COPYALLSYSTEMS(PREV, s, false);
  assert(LASTi(PREV[0]) >= 0);
  ASSERT_ONESYSTEM; 
  set_type(PREVSYSOF(cov), 0, type);
  return check2X(cov, vdim, vdim, frame, false);
}


int CheckPos2Neg(model *cov, int vdim, Types frame, int ntype,
		 domain_type maxdom) {
  /// genutzt von $

  // function must be preceeded by COPYALLSYSTEMS
#define ntype3 3
  int Err=ERRORFAILED;
  // statselect[nsel]={STATIONARY, VARIOGRAM, COVARIANCE, GEN_VARIOGRAM};
  Types typeselect[ntype3] = {PosDefType, VariogramType, NegDefType};
  assert(LASTi(CALLING[0]) >= 0);
  COPYALLSYSTEMS(PREV, CALLING, false);
  
  if (isAnySpherical(PREVISO(0))) ntype = 1;// Variogram does not make sense on  the sphere
  if (isAnyIsotropic(PREVISO(0))) maxdom = XONLY;
  
  for (int t=0; t<ntype; t++) {
    for (int s=XONLY; s<=maxdom; s++) {
      //printf("1 statt 0\n");
      set_system_type(PREV, typeselect[t]); 
      set_system_domain(PREV, (domain_type) s);
      if ((Err = check2X(cov, vdim, vdim, frame, true)) == NOERROR) return Err; 
    }
  }
  //printf("maxdom=%d\n", maxdom);  BUG;

  return Err;
}




// served by CHECK_NO_TRAFO // err 
int check2Xnotrafo(model *cov, int logicaldim, int xdimprev,
		   Types type, domain_type dom, isotropy_type iso,
		   int vdim, Types frame) {
  return check2X(cov, logicaldim, xdimprev, type, dom, iso,
		 vdim, vdim, frame, false);
}


// served by  CHECK // err 
int check2X(model *cov, int logicaldim, int xdimprev,
	    Types type, domain_type dom, isotropy_type iso,
	    int vdim, Types frame) { 
  // genutzt von dirct -> rms

  //  printf("check2x: type=%.50s\n", TYPE_NAMES[type]);
  
 return check2X(cov, logicaldim, xdimprev, type, dom, iso,
		 vdim, vdim, frame, true);
}

// served by  // err 
int check2X(model *cov, int logicaldim, int xdimprev, Types type, 
	    domain_type dom, isotropy_type iso, int *vdim, Types frame) {
  // genutzt von Simulate
  return check2X(cov, logicaldim, xdimprev, type, dom, iso, 
		 vdim[0], vdim[1], frame, true);
}

// served by CHECK_THROUGHOUT // err 
int check2Xthroughout(model *cov, model *calling, 
		      Types type,  domain_type dom, isotropy_type iso,
		      int vdim, Types frame) {
  assert(LASTi(SYSOF(calling)[0]) >= 0);
  COPYALLSYSTEMS(PREV, SYSOF(calling), false);
  set_system_type(PREV, type);
  int last = PREVLASTSYSTEM;
  if (dom != KEEPCOPY_DOM) for (int j=0; j<=last; j++) set_dom(PREV, j, dom);
  if (iso != KEEPCOPY_ISO) for (int j=0; j<=last; j++) set_iso(PREV, j, iso);
  return check2X(cov, vdim, vdim, frame, true);
}


// served by macro CHECK // err
int check2X(model *cov, int logicaldim, int xdimprev,
	    Types type, domain_type dom, isotropy_type iso,
	    int vdim0, int vdim1, Types frame, bool coord_trafo) {
  // genutzt von RMS

  //  PMI(cov->calling);
  
  set_system(PREV, 0, logicaldim, UNSET, xdimprev, type, dom,
	     equalsSpaceIsotropic(iso) ?
#if MAXSYSTEMS == 1    
	     DOUBLEISOTROPIC
#else
	     ISOTROPIC
#endif
	     : (equalsUnreduced(iso) && cov->calling != NULL)//neu, apr 15
	     ? ISO(CALLING, 0)
	     : iso
	     );

  //  printf("check2xx %.50s\n", NAME(cov));
  
#if MAXSYSTEMS > 1
  if (equalsSpaceIsotropic(iso))
    set_system(PREV, 1, 1, UNSET, 1, SameAsPrevType, dom, ISOTROPIC);
#endif

  return check2X(cov, vdim0, vdim1, frame, coord_trafo);
}


int basic_asserts(model *cov, Types frame, bool *coord_trafo) {
    // erst bei check unten
  assert(cov != NULL);
 int 
    prev_lastsystem = PREVLASTSYSTEM;
  model *calling = cov->calling;
  defn //*P = Cov List + prev -> nr, //  nicht gatternr
    *C = DefList + COVNR;        //  nicht gatternr
  KEY_type *KT = cov->base;
 
  cov->checked = false;
  assert(cov != NULL);
  //if (cov->base == NULL) {printf("basic %d %.50s\n", cov->calling==NULL,NAME(cov)); crash(); BUG;}
  //  if (cov->base != NULL) crash();
  //  PMI(cov);
  // if (cov->base == NULL) BUG;
  assert(cov->base != NULL);
  SPRINTF(KT->error_loc, "'%.50s'", NICK(cov));
  if (PL >= PL_COV_STRUCTURE) { 
    if (calling == NULL) { PRINTF("\n"); }
    LPRINT("%s\n", KT->error_loc); 
  }


  if (false) {
    //    PMI0(cov);
    PRINTF("InternalCov.cc: frame %d %.50s %d %d\n", frame,TYPE_NAMES[frame],equalsInterface(frame), //
	 //  equalsProcess(frame) ||
	 equalsGaussMethod(frame) ||
	 equalsBrMethod(frame) ||
	 equalsEvaluation(frame) ||
	 equalsInterface(frame) ||
	 equalsSmith(frame) ||
	 equalsSchlather(frame) ||
	 equalsRandom(frame) ||
	 equalsTrend(frame) ||
	 equalsNormedProcess(frame) ||
	 equalsPoissonGauss(frame) ||
	 equalsPoisson(frame)
	 );
  }
  
  
  assert(//equalsVariogramType(frame) ||
	 // equalsPosDefType(frame) ||
	 //equalsProcess(frame) ||
	 equalsGaussMethod(frame) ||
	 equalsBrMethod(frame) ||
	 equalsEvaluation(frame) ||
	 equalsInterface(frame) ||
	 equalsSmith(frame) ||
	 equalsSchlather(frame) ||
	 equalsRandom(frame) ||
	 equalsTrend(frame) ||
	 equalsNormedProcess(frame) ||
	 equalsPoissonGauss(frame) ||
	 equalsPoisson(frame) ||
	 equalsLikelihood(frame)
	 );

  // basic input checks
  //  printf("basic %.50s %.50s(%d)\n", NAME(cov), TYPE_NAMES[frame], frame);
  
  if (isManifold(frame) || isBad(frame)) {
    ERR_RETURN(1, "frame undefined");
  }

   //  printf("'%.50s' variant = %d %d %d is.interface=%d\n", NAME(cov), cov->variant, frame, BadType, isInterface(cov));
 
  
  if (calling != NULL && isInterface(cov)) {
    // PMI(cov);    crash();
    ERR_RETURN1(2, "'%.50s' may be used only as top model", NAME(cov));
  }

  for (int s=0; s<=prev_lastsystem; s++) {
    Types curtype = PREVTYPE(s);
    isotropy_type curiso = PREVISO(s);
    assert(!equalsIsoMismatch(curiso));

    if (isManifold(curtype) || isBad(curtype)) {
      ERR_RETURN2(3, "type '%.50s' not allowed for %.50s", TYPE_NAMES[curtype], NAME(cov));
    }
    if (equalsVariogram(curtype) && isAnySpherical(curiso)) {
      ERR_RETURN(4, "variograms do not make sense on spheres");
    }
    if (equalsKernel(PREVDOM(s)) &&
	(isAnyIsotropic(PREVISO(s)) || isSpaceIsotropic(PREVISO(s)) 
	 || equalsTrend(curtype))) {
      //      printf("%.50s %.50s %.50s\n", DOMAIN_NAMES[PREVDOM(s)],  ISO_NAMES[PREVISO(s)], TYPE_NAMES[curtype]);
      //PMI(cov);
      //APMI0(cov);
      RETURN(5, ERRORFAILED);
    }
    if (PREVXDIM(s) < 1) ERR_RETURN(5, "dimension less than 1"); 
     
    if (false && PL > PL_STRUCTURE) {
      LPRINT("#%d[%s -> %s] requ: %s, %s; has: %s, %s\n", s,
	     calling == NULL ? "NULL" : NICK(calling), NICK(cov),
	     Short(0, DOMAIN_NAMES[PREVDOM(s)]), Short(1, ISO_NAMES[curiso]),
	     PREVDOM(s) == DOM(C->systems[0], 0) 
	     ? "~" : Short(2, DOMAIN_NAMES[DOM(C->systems[0], 0) ]),  
	     curiso == ISO(C->systems[0],0)
	     ? "~" : Short(3,ISO_NAMES[ISO(C->systems[0], 0)])); 
    }
  }

  // assert(GLOBAL.coords.allow_earth2cart);
  
  if (calling != NULL && isEarth(PREVISO(0))) {
    *coord_trafo = COVNR == TRAFO || 
      (*coord_trafo && GLOBAL.coords.allow_earth2cart && !isAnyDollar(calling));
  } else {
    *coord_trafo = false;
    for (int s=1; s<=prev_lastsystem; s++) 
      if (isEarth(PREVISO(s))) BUG; // only the first system may be earth
  }

  NOERROR_RETURN(7);
}


/*
check2X : anything concerning OWN, except DOMain (and xdim, logdim if no TRAFO)
  PREV has already been set
  COPYALLSYSTEMS(OWN, DEF, true);
  switch owniso 
    case SubModelI :
      SetXdimLogdim: sets alwaus iso;  logdim & xdim if prev=EARTH & owniso=CART 
    case not SubModelI: 
      case Unreduced :  
        COPYALLSYSTEMS(OWN, PREV, prev)
      case isParamDepI : 
        setDI : all setting of OWN that are parameter specific
      case not is ParamDepI :
        none (i.e. own == def)
	set_own_domain_and_check  
    setDI
    c hange_coord_system
      setgatter_but_nr (trying if possible)
      sets TRAFO for change of coord. system (inkl. Gatter.xdim,logd,dom,type)
    sets OWNDOMoop
    C->check
    TypeConsistency
    checkDims (checks vdims and sets max(x)dim)
    check_within_range
    SetGatterNr
*/



int check2Xintern(model *cov, int vdim0, int vdim1, Types frame,
		  bool coord_trafo) {

  //PMI(cov);
  
  if (COVNR == DOLLAR && equalsProcess(PREVTYPE(0))) BUG; // crash();

  
  //  printf("C HECK2X%.50s %.50s(%d) call=%ld root=%ld %ld\n", NAME(cov), TYPE_NAMES[PREVTYPE(0)], PREVTYPE(0), (long) cov->calling, (long) cov->root, (long) cov->base);
  defn //*P = Cov List + prev -> nr, //  nicht gatternr
    *C = DefList + COVNR;        //  nicht gatternr
  int 
    Err = ERRORFAILED,
    prev_lastsystem = PREVLASTSYSTEM,
    variants = C->variants;
  model *calling = cov->calling;
  bool left[MAXVARIANTS];
  for (int i=0; i<variants; left[i++] = true);
  set_trafo(UNSET);

  if (equalsLikelihood(frame) && cov->calling != NULL) {
    if (TOTALXDIM(PREVSYSOF(cov->calling)) != TOTALXDIM(SYSOF(cov->calling)))
      cov->frame = EvaluationType; // wann passiert das?? frage am 26.2.19
  }
  

  assert(vdim0 != 0 && vdim1 != 0);
  if ((Err=basic_asserts(cov, frame, &coord_trafo)) != NOERROR) return Err;
  cov->frame = frame;  
  if (calling != NULL) cov->prevloc = PLoc(calling);
 
  // collect possibilities for coordinate systems
  // note: always try to keep coordinate system, but always try to 
  // change to a simpler representation within the coordinate system

  // set / try out anything concerning OWN, except DOMain

  //  printf("entering check2X %.50s (%d) coordtrafo = %d (%d %d)\n", NAME(cov), cov->zaehler, coord_trafo, C->Iallowed != NULL, !cov->IallowedDone);

  if (C->Iallowed != NULL && !cov->IallowedDone) {
    // printf("get allowed %.50s\n", NAME(cov));
    cov->IallowedDone = true; // muss zwingend vor dem Aufruf stehen
    bool *I = cov->allowedI;
    if (C->Iallowed(cov)) {
      for(int i =FIRST_ISOUSER; i<= LAST_ISOUSER; I[i++] = false);
      I[PREVISO(0)] = true;
      //      printf("XXX XXX XXX\n");
      //BUG;
    } else {
      int i =(int) FIRST_ISOUSER;
      for( ; i<=(int) LAST_ISOUSER; i++) if (I[i]) break;

      // PMI(cov);      printI(cov);
      
      if (i > LAST_ISOUSER) BUG;
      //
#ifdef xdebug
      printI(cov); //
#endif
    }
  } else if (C->Iallowed != NULL) {
    // printf("from former: ");
    // printI(cov);
    // BUG;
  }

  if (C->Dallowed != NULL && !cov->DallowedDone) {
    cov->DallowedDone = true;
    assert(FIRST_DOMAIN == XONLY);
    if (C->Dallowed(cov)) {
      for(int i=FIRST_DOMAIN ; i<=LAST_DOMAINUSER; cov->allowedD[i++] = false);
      cov->allowedD[PREVDOM(0)] = true;
      //printf("XXX XXX XXX DDD\n");
      //BUG;
    } else {
      int i = (int)FIRST_DOMAIN;
      for (; i <= (int) LAST_DOMAINUSER; i++) if (cov->allowedD[i]) break;

      //     printD(cov);
      
      if (i > LAST_DOMAINUSER) BUG;
      if (PREVDOM(0) < i || (i > XONLY && isAnyIsotropic(PREVISO(0)))) {
	//	printf("here\n"); TREE0(cov); APMI(cov);
	RETURN(10, CERRORWRONGDOM);
      }
#ifdef xdebug
      printD(cov);//
#endif
    }

    if (C->Dallowed != NULL) {
      //intf("from former: "); printD(cov);
    }
    cov->variant = 0;
    domain_type dom = DEFDOM(0);
    if (dom <= LAST_DOMAINUSER && PREVDOM(0) < dom) RETURN(11, ERRORWRONGDOM);
  }
#ifdef xdebug
  TREE0(cov); //
  //if (cov->calling != NULL) crash();
#endif
#ifdef ydebug 
  //PMI0(cov->calling); //
  // if (cov->Snugget !=NULL) printf("spatial nugget = %d\n", cov->Snugget->spatialnugget);
  //PMI0(cov);
#endif
  
 
  for (ext_bool ct=falsch; ct<=(ext_bool) coord_trafo;
       ct=ct == falsch ? wahr : BothOK) { // c*oord-t*rafo
    //intf("ct=%d of %d\n", ct, coord_trafo);
    for (cov->variant=0; cov->variant < variants; cov->variant++) {
      //
#ifdef ydebug
     PRINTF("variant=%d %d\n", cov->variant, ct);
#endif     
     //  PRINTF("C HECKct%.50s %.50s(%d)\n", NAME(cov), TYPE_NAMES[PREVTYPE(0)], PREVTYPE(0));
      if (C->TypeFct == NULL &&
	  isBad(TypeConsistency(PREVTYPE(0), DEFTYPE(0)))) {
	//	
#ifdef ydebug
	PRINTF("variant=%d prev=%.50s %.50s %d %d\n", cov->variant,
	       TYPE_NAMES[PREVTYPE(0)], TYPE_NAMES[DEFTYPE(0)],
	       isTrend(DEFTYPE(0)), TypeConsistency(PREVTYPE(0), DEFTYPE(0)));
#endif     
	//     	PMI(cov->calling);
	Err = ERRORTYPECONSISTENCY;
	continue;
      }
      if (!left[cov->variant]) continue; // never happens when ct=false
      isotropy_type owniso = DEFISO(0); // OWNISO(0); ///
      //printf("A here ct = %d %.50s;  %d %.50s --> %d %.50s\n", ct, NAME(cov), PREVISO(0),ISO_NAMES[PREVISO(0)],  owniso, ISO_NAMES[owniso]); 
      if (isFixed(owniso) && !atleastSpecialised(owniso, PREVISO(0))) {
	// short cut for most cases
	Err = ERRORFAILED;
	//printf("X\n");	
	continue;
      }
      
      int n = prev_lastsystem;
      //
#ifdef ydebug
      PRINTF(" ------------ variant=%d owniso='%.50s' in '%.50s' [n=%d] \n", cov->variant, ISO_NAMES[owniso], NAME(cov), n);
#endif     

      if (equalsSubModelI(owniso)) {
	//	printf("%.50s :: SubmodelI\n", NAME(cov));	
	assert(C->Iallowed != NULL);
	COPYALLSYSTEMS(OWN, DEF, true);
	set_xdim(OWN, 0, PREVXDIM(0));
	set_logdim(OWN, 0, PREVLOGDIM(0));

	if (C->setDI != NULL) { // at least settbm does not set OWNISO
	  //
	  // bisschen doppelt gemoppelter Code. Siehe SetXdimLogdim & set_own...
	  if (!C->setDI(cov)) BUG;
	  owniso = OWNISO(0);
	  //	  printf("\n\nsSETDI !!! %d %.50s\n\n", owniso, ISO_NAMES[owniso]);
	  
	  if (!equalsSubModelI(owniso)) {
	    isotropy_type iso[MAXSYSTEMS + 1];
	    if (!cov->allowedI[owniso]) {
	      // printf("%.50s not allowed in DI case\n", ISO_NAMES[owniso]);
	      continue;
	    }
	    iso[0] = owniso;
	    // YYY	    
	    if ((Err = SetXdimLogdim(cov, iso, prev_lastsystem)) == NOERROR) {
#if MAXSYSTEMS > 1
	      BUG;	  // call join_systems(...) here !
#endif
	      
	      // XXXXXX	      
	      if ((Err = set_own_domain_and_check(cov, vdim0, vdim1, ct))
		  == NOERROR) {
		return Err;
	      }	else {
		assert(false);
	      continue;
	    }  
	    } else {
	      assert(false);
	      continue;
	    }
	  }
	}
	
	if (equalsSubModelI(owniso)) {
	  //	  printf("%.50s ::|| SubmodelI\n", NAME(cov));	
	//	  printf("XAA %d %.50s\n", COVNR, NAME(cov));	
	  isotropy_type iso[MAXSYSTEMS],
	    start[MAXSYSTEMS] = { ISO_MISMATCH }, // rest initialised with 0
	    end[MAXSYSTEMS] = { ISO_MISMATCH };
	  assert(n < MAXSYSTEMS);
	  for (int s=0; s<=n; s++) {	  
	    isotropy_type 
	      previso = PREVISO(s),
	      newiso = CoordinateSystemOf(previso);
	    
	    if (ct == wahr && isEarth(newiso) && s==0) newiso = CARTESIAN_COORD;
	    // printf("hiere %d %d PREVISO = %.50s; %.50s; %.50s ######\n", s, n, ISO_NAMES[previso], ISO_NAMES[newiso], ISO_NAMES[OWNISO(0)]);
	    switch (newiso) {
	    case CARTESIAN_COORD : 
	      start[s] = FIRST_CARTESIAN;
	      //  isUnreducedCartesian(previso)
	      // ? (isotropy_type) (LAST_R EDUCEDxdim_CART + 1)
	      // : FIRST_CARTESIAN;
	      end[s] = CARTESIAN_COORD;
	      break;
	    case SPHERICAL_COORD :
	      // hilfskonstruktion, solange keine echten multiplen
	      // koordinatensysteme programmiert sind
	      start[s] =  OWNXDIM(0) > 2 ? SPHERICAL_SYMMETRIC :FIRST_SPHERICAL;
	      //isUnreducedSpherical(previso)
	      //  ? SPHERICAL_SYMMETRIC
	      //  : FIRST_SPHERICAL;
	      end[s] = LAST_TRUESPHERICAL;
	      break;
	    case EARTH_COORD :
	      //	      printf("xdim = %d\n", OWNXDIM(0));
	      //  start[s] = OWNXDIM(0) > 2 ? EARTH_SYMMETRIC : FIRST_EARTH;
	      start[s] = OWNXDIM(0) > 2 ? SPHERICAL_SYMMETRIC : FIRST_SPHERICAL;
	      
	      //isUnreducedEarth(previso)
	      // ? EARTH_SYMMETRIC 
	      //  : FIRST_EARTH;
	      end[s] = LAST_EARTH;
	      break;
	      //	  case  LOGCART_COORDS : 
	      //	    iso[s] = start[s] = isUnreducedLogCart(previso)
	      //	      ? LOGCART_SYMMETRIC
	      //	      : FIRST_LOGCART;
	      //	    end[s] = LAST_LOGCART;
	      //	    break;
	    default: BUG;
	    }
	    iso[s] = start[s]; 
	  }

	  Err = XERRORCHANGESYSTEM;
	  while (true) {
	    //
	    //
#ifdef ydebug
	    PRINTF("while !!!!!!!!!!!!!!!!!!! %.50s %d %.50s; [... %.50s] max-ct=%d ct=%d allwd=%d\n", NAME(cov), n, ISO_NAMES[iso[0]], ISO_NAMES[end[0]], coord_trafo, ct, cov->allowedI[iso[0]]); printI(cov); //
#endif
	    //printI(cov);
	    // for (int s=0; s<=n; s++) set_iso(OWN, s, iso[s]);
	    //printf("%.50s %d\n", ISO_NAMES[iso[0]], cov->allowedI[iso[0]]);
	    if (cov->allowedI[iso[0]] &&
		// YYY		
		(Err = SetXdimLogdim(cov, iso, prev_lastsystem)) == NOERROR) {
#if MAXSYSTEMS > 1
	      BUG;	  // call join_systems(...) here !
#endif
	      //
	      //
#ifdef ydebug
	      PRINTF(" ++++++++++++++= xxxxxx ct=%d [%.50s]\n", ct, ISO_NAMES[iso[0]]);
#endif
	      // XXXX
	      Err = set_own_domain_and_check(cov, vdim0, vdim1, ct);
	      //
#ifdef ydebug
	      PRINTF("back from set_own_dom %d %.50s err = %d\n", COVNR, NAME(cov), Err);
#endif
	      if (Err == NOERROR) return Err;

	      //
#ifdef ydebug
	      PRINTF("A not ok err=%d, '%.50s'\n", Err, cov->err_msg);
#endif	      
	    }
	    int s = 0;			
	    iso[s] = (isotropy_type) (iso[s] + 1);
	    //  printf("iso: s=%d %d %d %d !!!!!!!!!\n", s, iso[s], end[s], n);
	    assert(n < MAXSYSTEMS);
	    n = MIN(n, MAXSYSTEMS - 1);
	    while (s>=0 && s<=n && iso[s] > end[s]) {
	      assert(s < MAXSYSTEMS);
	      iso[s] = start[s];		
	      if (++s > n) break;
	      assert(s < MAXSYSTEMS);
	      iso[s] = (isotropy_type) (iso[s] + 1);
	    }
	    if (s > n) break;
	  } //  while true 
	} // still issubmodeli 
	left[cov->variant] = Err==CERRORCHANGESYSTEM;//indifferent when ct=true
	//printf("left = %d %d\n", cov->variant, left[cov->variant]);

      } else { // if not issubmodeli
	if (ct == wahr) continue;
        bool def_Unreduced = isPrevModelI(C) || equalsUnreduced(C);

	//	printf("def_unred %d %d; def_unred=%d %d %.50s\n",  isPrevModelI(C) , equalsUnreduced(C), def_Unreduced, hasFullXdim(owniso), NAME(cov));
	
	if (!def_Unreduced && hasFullXdim(owniso)) { // isUnred not eqUnred!
	  // E.g. Def = owniso = "Cart system"

	  // printf("prev->own %.50s\n", NAME(cov));
	  
	  COPYALLSYSTEMS(OWN, PREV, true);
	  set_iso(OWN, 0, owniso);
	  int s,
	    nn = PREVLASTSYSTEM;
	  for (s=0; s<=nn; s++)
	    if (!hasFullXdim(PREVISO(s))) break;
	  if (s <= nn) RETURN(13, ERRORREDUCED);
	  //	  print("CCC  variant %d %d\n", cov->variant, OWNISO(0));
	  // XXXX
	  if ((Err = set_own_domain_and_check(cov, vdim0, vdim1, BothOK))
	      == NOERROR) return Err;
	} else {
	  if (def_Unreduced) { // essentially == Unreduced
	    COPYALLSYSTEMS(OWN, PREV, true);
	  } else {
	    COPYALLSYSTEMS(OWN, DEF, true);

	    //	if (cov->allowedI[iso[0]] &&
	    isotropy_type iso[MAXSYSTEMS];
	    iso[0] = OWNISO(0);
	    if ((Err = SetXdimLogdim(cov, iso, 0)) != NOERROR) continue;
	    set_iso(OWN, 0, equalsPrevModelI(owniso) ? PREVISO(0) : owniso);
	    assert(equalsParamDepI(OWNISO(0)) || isFixed(OWNISO(0)));
	  }
	  if ((Err = set_own_domain_and_check(cov, vdim0, vdim1, BothOK))
	      ==NOERROR) return Err;    // HIER
	}	  
      } // !issubmodeli
    } // for cov->variant
    if (cov->variant >= variants) {
      assert(Err != NOERROR);
      cov->variant = UNSET;
    }
  } // for ct
  RETURN(14, Err);

  // ErrorHandling:
  //  SYSTEM_NULL(PREV, 0);
  //  cov->IallowedDone = cov->DallowedDone = false;
  //  RETURN_ERR(Err);
}


int SetXdimLogdim(model *cov, isotropy_type *Iso, int lastsystem) {
  //  printf("entering setxdimlogdim ");
  //  printf("SetX %.50s %.50s(%d)\n", NAME(cov), TYPE_NAMES[OWNTYPE(0)], OWNTYPE(0));
  // guessing what could be the right value for OWN.
  for (int s=0, sneu=0; s<=lastsystem; s++, sneu++) {
    isotropy_type iso = Iso[s];
    set_iso(OWN, sneu, iso);
    if (isCartesian(PREVISO(s))) { // Cartesian, no trafo      
      set_logdim(OWN, sneu, PREVLOGDIM(s));
      // printf("%.50s %.50s %d \n", NAME(cov), ISO_NAMES[iso], isAnyIsotropic(iso));
      if (isAnyIsotropic(iso)) {set_xdim(OWN, sneu, 1);}
      else if (equalsUnreduced(iso)) { set_xdim(OWN, sneu, PREVXDIM(s));}
      else if (equalsSpaceIsotropic(iso)) {
#if MAXSYSTEMS == 1	
	if (PREVXDIM(s) < 2) {
	  ERR_RETURN2(20, "'%.50s' not possible in %.50s", ISO_NAMES[iso], NAME(cov));
	}
	set_iso(OWN, sneu, DOUBLEISOTROPIC);
	assert(PREVLOGDIM(s) == PREVXDIM(s));
	set_xdim(OWN, sneu, 2);
#else	  
	if (s==0 && (PREVXDIM(s) > 1)) {
	  set_iso(OWN, sneu, ISOTROPIC);
	  assert(PREVLOGDIM(s) == PREVXDIM(s));
	  set_logdim(OWN, sneu, PREVLOGDIM(s) - 1);
	  set_xdim(OWN, sneu, 1);
	  (sneu)++; assert(sneu < MAXSYSTEMS);
	  set_iso(OWN, sneu, ISOTROPIC);
	  set_logdim(OWN, sneu, 1);
	  set_xdim(OWN, sneu, 1);
	} else ERR_RETURN(21, "split in spatial and temporal part not possible");
#endif
      } else set_xdim(OWN, sneu, PREVXDIM(s));
    } else if (isAnySpherical(PREVISO(s))) {
      if (isCartesian(iso)) { // Earth to cartesian, genuine trafo
	//         Spherical is senseless -- must be stopped elsewhere       
	assert(s == 0);
	set_logdim(OWN, sneu, 3);
	switch(iso) {
	case ISOTROPIC: set_xdim(OWN, sneu, 1); break;
	case DOUBLEISOTROPIC: ERR_RETURN(22, "not allowed"); break;
	case VECTORISOTROPIC: case SYMMETRIC: case CARTESIAN_COORD:
	  set_xdim(OWN, sneu, 3);
	  break;
	default: BUG;
	}
      } else { // stays within Earth or Spherical
	assert(s == 0);
	set_logdim(OWN, sneu, PREVLOGDIM(s));
	set_xdim(OWN, sneu, isAnyIsotropic(iso) ? 1 : PREVXDIM(s)); 
      } // not earth->cart
    } else { 
      BUG }
  } // end for

  //  printf("setxdim ok\n");
  NOERROR_RETURN(23);
}


int check_rec(model *cov) {
  defn *C = DefList + COVNR;
  // printf(".");
  if (!TrafoOK(cov, __FILE__, __LINE__) ||
      (cov->err_level >= LEVEL_DOM && cov->err_level <= LEVEL_LAST)) {
    // APMI0(cov);
    return false;
  }
  for (int i=0; i<cov->nsub; i++)
    if (!check_rec(cov->sub[i])) return false;
  for (int i=0; i<C->kappas; i++)
    if (cov->kappasub[i] != NULL && !check_rec(cov->kappasub[i])) return false;
  return true;
}



int check2X(model *cov, int vdim0, int vdim1, Types frame,
	    bool coord_trafo) {  
  //   int level = cov->err_level;
  //  printf("initial level%.50s\n", NAME(cov));
 
    
  int Err = check2Xintern(cov, vdim0, vdim1, frame, coord_trafo);
  model *calling = cov->calling;
  //  printf("out cov=%.50s (z=%d L=%d E=%d) calling=%.50s(z=%d, L=%d E=%d) err=%d \n", NAME(cov), cov->zaehler, cov->err_level,  cov->err,  cov->calling == NULL ? "NONE" : NAME(cov->calling),cov->calling == NULL ? -999 : cov->calling->zaehler, cov->calling == NULL ? -999 : cov->calling->err_level, cov->calling == NULL ? -999 : cov->calling->err, Err);

  if (Err >= ERRORM && Err <= ERRORMEND && calling != NULL &&
      calling->err_level < LEVEL_CHECK)
    STRCPY(calling->err_msg, cov->err_msg);
  
  // printf("\t\t\t    calling=%.50s(z=%d, L=%d E=%d) \n",NAME(cov->calling),cov->calling->zaehler, cov->calling->err_level, cov->calling->err);
  // PMIE(calling); PMIE(cov);

  assert(Err != NOERROR || check_rec(cov));
  
  return Err;
}

int set_own_domain_and_check(model *cov, int vdim0, int vdim1,
			     ext_bool coord_trafo) {
  // printf("entring set_own_domain %.50s\n", NAME(cov));
  //printf("set_own %.50s %.50s(%d)\n", NAME(cov), TYPE_NAMES[OWNTYPE(0)], OWNTYPE(0));
  KEY_type *KT = cov->base;
  defn *C = DefList + COVNR;
  if (coord_trafo == false && isEarth(PREVISO(0)) && isCartesian(OWNISO(0)))
    RETURN(30, ERRORWRONGISO);
  
  if (isManifold(OWNTYPE(0))) set_type(OWN, 0, PREVTYPE(0));

  isotropy_type iso = OWNISO(0);
  //printf("here\n");
  if ( C->setDI != NULL && (!isFixed(iso) || !isFixed(OWNDOM(0)))
       && !C->setDI(cov) ) { // muss i.A. ein 2. Mal
    // aufgerufen werden (1. Mal in check2x, da beim 1. Aufruf u.U.
    // noch nicht alle Daten vorhanden
    // Und insbesondere bei ParamDepI, wenn noch gar nichts gesetzt ist
    BUG;    
  }
  if (isFixed(iso) && iso != OWNISO(0)) BUG;
  
  // printf("B here %.50s;  %d %.50s --> %d %.50s\n", NAME(cov), PREVISO(0),ISO_NAMES[PREVISO(0)],  OWNISO(0), ISO_NAMES[OWNISO(0)]); 
  if (C->TypeFct == NULL) {
    if (isBad(TypeConsistency(PREVTYPE(0), cov, PREVISO(0)))) {
      //      printf("bad type = %d\n", BadType);
      //APMI0(cov);     
      RETURN(31, ERRORTYPECONSISTENCY);
    }
  } else {
    if (OWNISO(0) <= LAST_ISOUSER &&
	CoordinateSystemOf(OWNISO(0)) == CoordinateSystemOf(PREVISO(0)) && 
	!atleastSpecialised(OWNISO(0), PREVISO(0))) {
      //
      if (false) {
	PMI0(cov->calling);//
	printf("wrong iso: own=%.50s prev=%.50s coordtrafo=%d\n", ISO_NAMES[OWNISO(0)], ISO_NAMES[PREVISO(0)], coord_trafo);//
	//
	printI(cov->calling);//
	printI(cov);      //	
	//APMI0(cov);//
      }
      RETURN(31, ERRORWRONGISO);
    }
  }

  // Function itself sets OWNDOM and special check for spherical coord
  
  // but also calls
  // setDI, cov->allowedD
  // set_system_domain
  // c hange_coord_system sets everything concerning TRAFO and GATTER,
  //                     except GATTERNR
  //                 prev -> gatter (incl gatterxdim, gatterlogdim, gatteriso)
  //                 and prev -> own (incl. gatter: xdim,iso,dom,type,logdim)
  // check_kappas
  // C->check
  // setgatternr
  // check_range

 

  model *calling = cov->calling;
  domain_type
    first_dom = LAST_DOMAIN,
    last_dom = FIRST_DOMAIN,
    defdom = DEFDOM(0); 
  bool skipchecks = GLOBAL_UTILS->basic.skipchecks;
  
  //printf("gonna up here (%.50s) %.50s \n", NAME(cov), DOMAIN_NAMES[defdom]);

  switch (defdom) {
  case DOMAIN_MISMATCH : 
    if (!isRandom(cov)) BUG;    
    first_dom =last_dom = XONLY;
    break;
  case PREVMODEL_D :
    first_dom = last_dom = PREVDOM(0);
    break;
  case PARAMDEP_D : {
    domain_type owndom = OWNDOM(0);
    if (!equalsParamDepD(owndom)) first_dom = last_dom = owndom;
    else {
      PMI(cov->calling); //
      BUG;
    }
  }
    break;
  case KEEPCOPY_DOM : BUG; break;
  case SUBMODEL_D :
    first_dom = XONLY;
    last_dom = KERNEL;
    break;
  case XONLY : case KERNEL : first_dom = last_dom = defdom;
    // printf("hier!!\n");
    break;
  default : BUG;
  }
  //printf("dom %.50s %.50s\n", DOMAIN_NAMES[first_dom], DOMAIN_NAMES[last_dom]);
  //  PMI0(cov);
  if (isAnyIsotropic(OWNISO(0)) || isSpaceIsotropic(OWNISO(0))) {
    last_dom = XONLY;
  } 

  //printf("dom %.50s %.50s\n", DOMAIN_NAMES[first_dom], DOMAIN_NAMES[last_dom]);
  

  if (last_dom > PREVDOM(0)) last_dom = PREVDOM(0);
  
  // if (isSubModelD(last_dom) && KERNEL <= GATTERDOM(0)) last_dom = KERNEL;
   
  //  printf("\nfirst=%d last=%d def=%d %.50s %d\n", first_dom, last_dom, defdom, NAME(cov), !STRCMP(NAME(cov), "prod"));
  //assert(STRCMP(NAME(cov), "prod"));

  if (first_dom > last_dom) {
    if (calling == NULL) BUG;
    if (PL >= PL_ERRORS) {
      LPRINT("(%s dom.start=%d, end=%d)\n", NAME(cov), first_dom, last_dom);
    }
    
    //PMI0(cov);
    RETURN(32, ERRORNOSTATMATCH);
  }

  //
#ifdef ydebug
  PRINTF("cov=%s first/last/prev = %d %d %d; %d %d %d\n",
	 NAME(cov), first_dom, last_dom,
	 PREVDOM(0),
	 cov->DallowedDone,
	 C->Dallowed == NULL ? 9999 : cov->allowedD[0],
	 C->Dallowed == NULL ? 9999 :  cov->allowedD[1]);
#endif

  int dom,
    errs[LAST_DOMAINUSER + 1],
    levels[LAST_DOMAINUSER + 1],
    save_level = cov->err_level,
    save_err = cov->err;
  errorstring_type ERRSTR[LAST_DOMAINUSER + 1] = {""},
    save_errmsg;
  bool save_errorm = save_err >= ERRORM && save_err <= ERRORMEND,
    save_copy = cov->err_level >= LEVEL_CHECK && save_errorm;
  if (save_copy) STRCPY(save_errmsg, cov->err_msg);
  
  // printf("!! %s %d err_lev=%d err=%d '%s'\n", NAME(cov), cov->zaehler, cov->err_level, cov->err , cov->err_msg);

 
  for (dom = first_dom; dom <= last_dom; dom++) {
    // printf("cov=%s dom=%d\n", NAME(cov), dom);
    if (dom > PREVDOM(0) || (cov->DallowedDone && !cov->allowedD[dom]))
      LEVEL(LEVEL_DOM, CERRORNOSTATMATCH);       
    DEBUG(set_system_domain(OWN, (domain_type) dom));
    if ((errs[dom] = change_coord_system(cov, coord_trafo, ERRSTR[dom]))
	!= NOERROR) CONT(41);
  
    // FROM HERE ON, GATTER SHOULD ALWAYS BE WELL-DEFINED
    setdefault(cov, vdim0, vdim1);            
    if ((errs[dom] = checkkappas(cov, C->primitive, ERRSTR[dom])) != NOERROR)
      CONT(42);
    //    printf("@@ %s %d  err_lev=%d err=%d '%s'\n", NAME(cov), cov->zaehler, cov->err_level, cov->err , cov->err_msg);
    //

    errs[dom] = C->check(cov);   
    if (errs[dom] != NOERROR) {
      if (errs[dom] >= ERRORM && errs[dom] <= ERRORMEND)
	STRCPY(ERRSTR[dom], cov->err_msg);

      // A      PMI(cov->calling->calling);

      // printf("YY %s %d err_lev=%d err=%d '%s'\n", NAME(cov), cov->zaehler, cov->err_level, cov->err , cov->err_msg);
     
      CONT(LEVEL_CHECK); // 43
    }
    
    // Zaehler++;
    // printf("?? %s %d z=%d err_lev=%d err=%d '%s'\n", NAME(cov), cov->zaehler, Zaehler, cov->err_level, cov->err , cov->err_msg);
    //    assert(Zaehler != 86);

    //    printf(">>> %s %s %s\n", NAME(cov), ISO_NAMES[PREVISO(0)], ISO_NAMES[OWNISO(0)]);

    
    
    if (isBad(TypeConsistency(PREVTYPE(0), cov, PREVISO(0)))) {
      LEVEL(44, ERRORTYPECONSISTENCY);
    }
    if (cov->monotone == PARAM_DEP) BUG;
    if ((errs[dom] = checkDims(cov, vdim0, vdim1, ERRSTR[dom])) != NOERROR)
      CONT(45);
    if (!skipchecks &&
	(errs[dom]=check_within_range(cov, NAOK_RANGE, ERRSTR[dom])) != NOERROR)
      CONT(46);
    if ((errs[dom] = SetGatterNr(cov,  ERRSTR[dom])) != NOERROR) CONT(47);
    assert(errs[dom] == NOERROR);
    levels[dom] = LEVEL_LAST;
    break;
  }

  if (dom > last_dom) { // some error occured
    for (dom = first_dom; dom <= last_dom; dom++) { //could be a void loop
      //  printf("%s dom=%d [%d %d] \n", NAME(cov), dom, first_dom, last_dom);
      //    printf(">> %s %d err_lev=%d err=%d dom=%d level=%d error=%d '%s'\n", NAME(cov), cov->zaehler, cov->err_level, cov->err , dom, levels[dom], errs[dom], ERRSTR[dom]);
      assert(errs[dom] != NOERROR);
      if (levels[dom] > cov->err_level) {
	cov->err = errs[dom];
	cov->err_level = levels[dom];
	STRCPY(cov->err_msg, ERRSTR[dom]);
	// printf("up %s %d err_lev=%d err=%d '%s'\n", NAME(cov), cov->zaehler, cov->err_level, cov->err , cov->err_msg);
      }
    }
    
    if (cov->err_level == save_level) {//old err level is not passed, so restore
      cov->err = save_err;
      if (save_copy) STRCPY(cov->err_msg, save_errmsg);
      // printf("RE %s %d err_lev=%d err=%d '%s'\n", NAME(cov), cov->zaehler, cov->err_level, cov->err , cov->err_msg);
    }
    if (PL > PL_COV_STRUCTURE) { // && calling == NULL) {
       LPRINT("error in %s: end check2x level=%d err=%d %s\n",
	     NAME(cov), cov->err_level, cov->err, cov->err_msg);
      //APMI(cov);     
    }
    return cov->err; // #define RETURN
  }
  
  if (PL >= PL_COV_STRUCTURE) {
    LPRINT("Continuing '%s' (no err):\n", Nick(cov));
  }
  SPRINTF(KT->error_loc, "'%.50s'",
	  calling == NULL ? "parsing the model" : Nick(calling));

  GATTER_STORAGE;
  if (isDollar(cov) && cov->finiterange == wahr && 
      isAnySpherical(ISO(PREVSYSOF(cov->sub[0]), 0))) {
    // erst abfragbar, nachdem checked=true gesetzt wurde !
    double range;
    INVERSE(ZERO(cov), cov, &range);
    if (!(isSpherical(OWNISO(0)) && range <= M_PI) ||
	(isEarth(OWNISO(0)) && range <= 180)){
      ERR_RETURN2(49,
		  "model '%.50s' does not match required coordinate system '%.50s'",
		  NICK(cov->sub[0]),COORD_SYS_NAMES[GetCoordSystem(OWNISO(0))]);
    }
  }
  
  assert(!equalsUnreduced(ISO(C->systems[cov->variant], 0)) ||
	 OWNISO(0) == PREVISO(0));

  //printf("CHECKED = true %s", NAME(cov));
  cov->checked = true;

  ASSERT_GATTER(cov); 
  //  printf("set_own ENDE %s %s(%d)\n", NAME(cov), TYPE_NAMES[OWNTYPE(0)], OWNTYPE(0));
  RETURN(50, NOERROR);// nicht NOERROR_RETURN(50) !! Da diese fkt die letzte ist
}


#define RET_ENR(NR) return NR;
#define RERR(X) { STRCPY(ERRSTR, X); return ERRORM; }
#define RERR1(X,A) {SPRINTF(ERRSTR,X,A); return ERRORM; }
#define RERR2(X,A,B) {SPRINTF(ERRSTR,X,A,B); return ERRORM; }
#define RERR3(X,A,B,C) {SPRINTF(ERRSTR,X,A,B,C); return ERRORM; }
#define RERR4(X,A,B,C,D) {SPRINTF(ERRSTR,X,A,B,C,D); return ERRORM; }
#define RERR5(X,A,B,C,D,E) {SPRINTF(ERRSTR,X,A,B,C,D,E); return ERRORM; }
#define RERR6(X,A,B,C,D,E,F) {SPRINTF(ERRSTR,X,A,B,C,D,E,F); return ERRORM; }
#define RERR7(X,A,B,C,D,E,F,G) {SPRINTF(ERRSTR,X,A,B,C,D,E,F,G); return ERRORM;}

#define NRERR(NR,X) { STRCPY(ERRSTR, X); return NR; }


int setgatter_but_nr(model *cov, int minp, int maxp, int mino, int maxo,
		     errorstring_type ERRSTR){
  // hier nur verschiebung TRAFinnerhalb der coordinatensysteme
  // hier wird gatter gesetzt
  // minp-revious .. maxo-wn
			  
  defn *C = DefList + COVNR;        //  nicht gatternr
  int
    sp = minp, // s-tart-p-revious
    so = mino, 
    ep = sp, // e-nd-p-revious
    eo = so;
  assert(maxp == 1 && maxo == 1);
 
  bool def_Unreduced = isPrevModelI(C) || equalsUnreduced(C);

  if (def_Unreduced) {
    COPYALLSYSTEMS(GATTER, PREV, false);
    set_nr(GATTER, UNSET);    
    return NOERROR;
  } else if (equalsSubModelI(OWNISO(0))) {
    COPYALLSYSTEMS(GATTER, PREV, false);
    set_nr(GATTER, UNSET);
    return NOERROR; // 14.11.17 -- OK?
  } else if (equalsParamDepI(OWNISO(0)) || equalsParamDepD(OWNDOM(0))) {
    COPYALLSYSTEMS(GATTER, PREV, false);
    set_nr(GATTER, UNSET);
  }
  
  // first only fill missing logdims and missing xdims, if there are
  // still unsolved issues. Should happen only rarely
  while (sp < maxp && so < maxo) { // ist also keine echte Schleife
    if (equalsUnreduced(OWNISO(so))) {
      if (!equalsCoordinateSystem(PREVISO(sp)))
	NRERR(ERRORCHANGESYSTEM, "unexpected reduced system of coordinates");
 
      while ((ep == sp || (ep < maxp && isSameAsPrev(OWNTYPE(ep)))) && 
	     (eo == so || (eo < maxo && isSameAsPrev(OWNTYPE(eo))))) {  
	if (OWNLOGDIM(eo) == UNSET) set_logdim(OWN, eo, PREVLOGDIM(ep));
	else if (OWNLOGDIM(eo) != PREVLOGDIM(ep))
	  NRERR(ERRORCHANGESYSTEM, "mismatch of logical dimensions");
	if (OWNXDIM(eo) == UNSET) { set_xdim(OWN, eo, PREVXDIM(ep)); }
	else if (OWNXDIM(eo) != PREVXDIM(ep))
	  NRERR(ERRORCHANGESYSTEM, "mismatch of dimensions");
	eo++;
	ep++;
      }
      if (ep == maxp && eo == maxo) {
	sp = ep;
	so = eo;
	break;
      }
     
      if (ep == maxp || eo == maxo ||
	  (isSameAsPrev(OWNTYPE(ep)) xor isSameAsPrev(OWNTYPE(eo)))
	  ) NRERR(ERRORCHANGESYSTEM, "coordinate systems cannot be matched");
      sp = ep;
      so = eo;
      continue;
    }
    
    int 
      xdims = 0,
      logdims = 0,
      unsetxdims = 0,
      unsetlogdims = 0;
    isotropy_type c1 = CoordinateSystemOf(PREVISO(sp));
    while (ep < maxp && c1 == CoordinateSystemOf(PREVISO(ep)) &&
	   (ep == sp || isSameAsPrev(PREVTYPE(ep) ))) {
      xdims += PREVXDIM(ep);         assert(PREVXDIM(ep) > 0);
      logdims += PREVLOGDIM(ep);     assert(PREVLOGDIM(ep) > 0);
      ep++;
    } 

    //    printf("c1=%s\n", ISO_NAMES[c1]);
    isotropy_type
      c0 = CoordinateSystemOf(OWNISO(eo)),
      c2 = EssentialCoordinateSystemOf(PREVISO(sp));

    if (false)
      printf("%s %s %d %d so=%d %d c1=%d c0=%d c2=%d %s ep=%d\n", //
	     ISO_NAMES[c0], ISO_NAMES[c2], eo, maxo,
	     so,isSameAsPrev(OWNTYPE(eo)), c1, c0, c2, ISO_NAMES[c1], ep);
    
    while (eo < maxo && (eo == so || isSameAsPrev(OWNTYPE(eo))) && (c1==c0 || c2==c0)) {
      //printf("eo=%d %d\n", eo, so);
      if (OWNXDIM(eo) != UNSET) {
	//	printf("eo=%d %d\n", eo, OWNXDIM(eo));
	xdims -= OWNXDIM(eo);
	assert(OWNXDIM(eo) > 0);
      } else unsetxdims++;
      if (OWNLOGDIM(eo) != UNSET) {
	logdims -= OWNLOGDIM(eo);
 	assert(OWNLOGDIM(eo) > 0);
      } else unsetlogdims++;
      eo++;
      if (eo < maxo) c0 = CoordinateSystemOf(OWNISO(eo));
    }

    if (eo == so) { // no machting possible without genuine coord trafo

      //APMI(cov);
      
      NRERR(ERRORCHANGESYSTEM,
	    "Non-matching coordinate systems -- coord trafo needed");
    }

    //    PMI(cov);
    // printf("%d %d  %d %d  %d  %d %d: %d %d\n", xdims < 0, unsetxdims > 1, logdims < 0, unsetlogdims > 1, (unsetxdims == 1 && xdims == 0), /* (unsetxdims == 0 && xdims > 0), passiert von z.B. Cart -> iso */ (unsetlogdims == 1 && logdims == 0), (unsetlogdims == 0 && logdims > 0), unsetlogdims, logdims);
    
    if (xdims < 0 || unsetxdims > 1 ||
	logdims < 0 || unsetlogdims > 1 || 
	(unsetxdims == 1 && xdims == 0) ||
	// (unsetxdims == 0 && xdims > 0) || passiert von z.B. Cart -> iso
	(unsetlogdims == 1 && logdims == 0) ||
	(unsetlogdims == 0 && logdims > 0)
	) {

      assert(false);
      
      NRERR(ERRORCHANGESYSTEM,
	    "non-matching dimensions of the coordinate systems");      
    }
    //    printf("ccs here ...\n");
   
    if (unsetxdims == 1) 
      for (int i=so; i<eo; i++)
	if (OWNXDIM(i) == UNSET) {
	  set_xdim(OWN, i, xdims);
	  break;
	}
    if (unsetlogdims == 1) 
      for (int i=so; i<eo; i++)
	if (OWNLOGDIM(i) == UNSET) {
	  set_logdim(OWN, i, logdims);
	  break;
	}    
    sp = ep;
    so = eo;
  } // end while
  //  printf("ccs here !!\n");
 

  if ((sp < maxp) xor (so < maxo))
    NRERR(ERRORCHANGESYSTEM, "non-matching coordinate systems"); 


  //////// jetzt werden Gatter gesetzt
  int curprev = minp,
    curlogdim = PREVLOGDIM(curprev);
  
  // if no real coordinate system transformation, 
  // it still my happen that a coordinate system is split or several are
  // unified. This must be fixed in GATTER
  COPYALLSYSTEMS(GATTER, OWN, false);
  set_nr(GATTER, UNSET);
  
  //PSYS(cov); 
  //   printf("ccs here !!\n");
 
  for (int i=mino; i<maxo; i++) {                 
    set_xdim(GATTER, i, PREVXDIM(curprev));
    set_iso(GATTER, i, PREVISO(curprev));
    set_dom(GATTER, i, PREVDOM(curprev));
    set_type(GATTER, i, PREVTYPE(curprev));
    curlogdim -= OWNLOGDIM(i);      
    bool fullxdim = hasFullXdim(OWNISO(i));
    if (curlogdim != 0) {
      if (!fullxdim) NRERR(ERRORCHANGESYSTEM, "coordinate change does not work.");
      while (curlogdim < 0) {
	curprev++;
	if (curprev > PREVLASTSYSTEM) BUG;
	if (PREVISO(curprev) != GATTERISO(i) ||
	    PREVDOM(curprev) != GATTERDOM(i) ||
	    PREVTYPE(curprev) != GATTERTYPE(i))
	  NRERR(ERRORCHANGESYSTEM, "coordinate change does not work");
	curlogdim += PREVLOGDIM(curprev);
      }
    } //  kein else
    if (curlogdim == 0) {
      curprev++;
      if (curprev > PREVLASTSYSTEM) {
	assert(eo != -9999);
	if (i < eo - 1) BUG;
      } else curlogdim = PREVLOGDIM(curprev);	
    }
    //    PSYS(cov); PMI(cov->calling);
    if (fullxdim) {
      if (GATTERXDIM(i) == UNSET) {set_xdim(GATTER, i, GATTERLOGDIM(i));}
      else if (GATTERXDIM(i) != GATTERLOGDIM(i)) {
	//
	//	printf("i=%d %d %d\n", i, GATTERXDIM(i), GATTERLOGDIM(i)); // 0, 1, 3
	//	APMI(cov->calling);
	BUG;
      }
    } else {
      if (GATTERXDIM(i) != PREVXDIM(i)) BUG;
    }
    set_cumxmit(GATTER, i, GATTERXDIM(i));
    //if (i>0) set_cumxmit(GATTER, i, CUMXMIT(GATTER, i) + CUMXMIT(GATTER, i-1)); 
  }
  //  PSYS(cov); 
  set_trafo(UNSET);
  return NOERROR;
}



int change_coord_system(model *cov, ext_bool coordinate_trafo,
			errorstring_type ERRSTR) {
  // ACHTUNG!!! FUNKTION DARF NUR IN InternalCov.cc VERWENDET WERDEN !!
  // hier wird das grosse Rad gedreht: Wechsel zwischen den Systemen
  // function itself sets TRAFO, 
  bool checkerror = true;
  int err = CERRORCHANGESYSTEM;

  /* version until end august 2018
     int err = setgatter_but_nr(cov, 0, PREVSYSTEMS, 0, OWNSYSTEMS);
     if (err == NOERROR || !coordinate_trafo) RET_ENR(err);
  */

  if (coordinate_trafo != wahr) { // version since august 2018
    if (!(isEarth(PREVISO(0)) && isCartesian(OWNISO(0))))// kann nur BothOK sein
      //               wegen abfrage in set_own_domain_and_check
      err = setgatter_but_nr(cov, 0, PREVSYSTEMS, 0, OWNSYSTEMS, ERRSTR);
    if (coordinate_trafo == falsch || err != CERRORCHANGESYSTEM) {
      //PMI(cov);
      //printf("hier fehler %d err=%d is %d %d\n", coordinate_trafo, err, isEarth(PREVISO(0)) ,  isCartesian(OWNISO(0)) );
      RET_ENR(err);
    }
  }
  
  if (isnowNegDef(cov) && equalsXonly(PREVDOM(0))) RET_ENR(ERRORWRONGDOM);

  COPYALLSYSTEMS_COND((cov)->gatter, OWN, true); // nicht GATTER, da unset
  set_trafo(UNSET); 
   
  isotropy_type 
    isogatter = GATTERISO(0),
    isoprev = PREVISO(0);
  int 
    err2,
    logdimprev = PREVLOGDIM(0),
    xdimprev = PREVXDIM(0);
  
  switch(isoprev) {
  case EARTH_COORD : case EARTH_SYMMETRIC: {
    int additional_xdim = xdimprev - 2, // time and/or height
      additional_logdim = logdimprev -2;
    assert(isogatter >= 0);
    if (isCartesian(isogatter)) {
      if (xdimprev != logdimprev) BUG;
      if (STRCMP(GLOBAL.coords.newunits[0], UNITS_NAMES[units_km]) == 0){
	set_trafo(equalsGnomonic(isogatter) ? EARTHKM2GNOMONIC
		  : equalsOrthographic(isogatter) ?  EARTHKM2ORTHOGRAPHIC
		  : EARTHKM2CART);
      } else if (STRCMP(GLOBAL.coords.newunits[0], 
			UNITS_NAMES[units_miles]) == 0) {
	set_trafo(equalsGnomonic(isogatter) ? EARTHMILES2GNOMONIC  
		  : equalsOrthographic(isogatter) ?  EARTHMILES2ORTHOGRAPHIC
		  : EARTHMILES2CART);
      } else {
	RERR4("only units '%.50s' and '%.50s' are allowed. Got '%.50s' (user's '%.50s').",
	      UNITS_NAMES[units_km], UNITS_NAMES[units_miles], 
	      GLOBAL.coords.newunits[0], GLOBAL.coords.curunits[0]);
      }
      if (isEarthProjection(isogatter)) {
	set_xdim(GATTER, 0, 2 + additional_xdim);
	set_logdim(GATTER, 0, 2 + additional_logdim);
      } else if (isCartesian(isogatter)) {
	set_xdim(GATTER, 0, 3 + additional_xdim);
	set_logdim(GATTER, 0, 3 + additional_logdim);
      } else {
	BUG;
      }
      set_dom(GATTER, 0, PREVDOM(0));
      set_type(GATTER, 0, PREVTYPE(0));
    } else {
      BUG;
      RET_ENR(ERRORODDCOORDTRAFO); // NotProgrammedYet("");
    }
  }
    break;
  case EARTH_ISOTROPIC : case SPHERICAL_ISOTROPIC :
    RET_ENR(ERRORWRONGISO);
  default:    
    RET_ENR(ERRORODDCOORDTRAFO); // NotProgrammedYet("");
  }
  

  if ((err2 = checkEarth(cov)) != NOERROR) {
    if (err == NOERROR || (err == ERRORNOSTATMATCH && !checkerror))
      err = err2;
    //continue;
  }

  //ERR("XDIM != UNSET beruecksichtigen"); // irgendwas fehlt noch
  // zusaetzlich abgleichen mit Coord trafo, um insb. spaceiso->iso
  // hinzubekommen. Muss gelockert werden, dass gleiches ISO

  
  //assert(*newlogdim == 4);
  RET_ENR(NOERROR);
}



int SetGatterNr(model *cov, errorstring_type ERRSTR) {
  // ACHTUNG!!! FUNKTION DARF NUR IN InternalCov.cc VERWENDET WERDEN !!
  int n = OWNSYSTEMS;
  assert(n == GATTERSYSTEMS);
  for (int s=0; s<n; s++) {
    domain_type owndom = OWNDOM(s),
      gdom = GATTERDOM(s);
    isotropy_type 
      owniso = OWNISO(s),
      giso = GATTERISO(s);
    int nr;
    
    if (gdom < owndom) 
      RERR2("Cannot call more complex models ('%.50s') from simpler ones ('%.50s')",
	    DOMAIN_NAMES[(int) owndom], DOMAIN_NAMES[(int) gdom]);

    //////////////
    // cartesian
    //////////////
    if (isCartesian(owniso)) {
      if (!isGenuineCartesian(owniso)) {	
	//RERR1("'%.50s' is not genuinely cartesian", ISO_NAMES[(int) owniso]);
      }
      // bool isoOK = isoprev == isonext || 
      // (isoprev > isonext && isonext <= CARTESIAN_COORD);
      //if (!isoOK) 
      if ((owniso > giso && isGenuineCartesian(giso)))
	RERR2("cannot call more complex models ('%.50s') from simpler ones ('%.50s')",
	      ISO_NAMES[(int) owniso], ISO_NAMES[(int) giso]);
      
      if (equalsXonly(owndom) &&
	  (equalsKernel(gdom) || hasFullCartesian(giso))) {	
	switch (owniso) {
	case ISOTROPIC :
	  nr = S2ISO;
	  break;
	case DOUBLEISOTROPIC :
	  nr = S2SP;
	  break;
	case VECTORISOTROPIC: case SYMMETRIC: 
	case CARTESIAN_COORD: case GNOMONIC_PROJ: case ORTHOGRAPHIC_PROJ:
	  nr = S2S;
	  break;
	case UNREDUCED:
	  nr = SId;
	  break;
	default: BUG;
	}
      } else {    
	if (equalsXonly(gdom)) { 
	  switch(giso) {
	  case ISOTROPIC :
	    nr = ISO2ISO; 
	    break;
	  case DOUBLEISOTROPIC :
	    nr = (equalsIsotropic(owniso)) ? SP2ISO : SP2SP;
	    break;
	  default: 
	    PRINTF("SetGatter prev=%s; %s\n          next=%s %s\n",
		   DOMAIN_NAMES[gdom], ISO_NAMES[giso], 
		   DOMAIN_NAMES[owndom], ISO_NAMES[owniso]);
	    BUG;
	  }
	} else { // K ERNEL gdom und owndom
	  nr = SId;
	}
      }    
    }


    else if (isAnySpherical(owniso)) {
      if(!isAnySpherical(giso)) BUG; //to do?: back from cartesian 
      //                                  not possible currently
      if (equalsUnreduced(giso)) {
	if (giso != owniso || gdom != owndom)
	  RERR("unclear transformation of earth/spherical coordinates");
	nr = isEarth(giso) ? E2E : Sph2Sph;
      }
      if (equalsXonly(owndom)) {     
	if (isAnySphericalIso(owniso)) {	  
	  if ( (!isAnySphericalIso(giso) && !equalsKernel(gdom)) || 
	       (isAnySphericalIso(giso) && !equalsXonly(gdom))) {
	    //   A 	    PMI(cov);
	    RERR("impossible spherical trafo");
	  }
	  
	  if (isEarth(giso)) nr = isSpherical(owniso) ? E2SphIso : E2EIso;
	  else nr = isSpherical(owniso) ? Sph2SphIso : ISO_MISMATCH; 	
	} else { // next: xonly, not iso; e.g. for shape functions
	  assert(!isAnySphericalIso(owniso));
	  if (equalsKernel(gdom)) RERR("earth trafo not possible");// mathemati-
	  // cally likely not correct
	  assert(equalsXonly(gdom));
	  if (isAnySphericalIso(giso))  
	    RERR("mathematically not correct coordinate transformation"); // mathematically likely not correct. to do
	  if (isEarth(giso)) nr = (isSpherical(owniso) ? E2Sph : E2E);
	  else nr =  (isSpherical(owniso) ? Sph2Sph : ISO_MISMATCH);	
	}
      } else { // owndom == K ERNEL
	assert(!isAnySphericalIso(owniso));
	assert(equalsKernel(gdom));
	assert(!isAnySphericalIso(giso));
	
	if (isEarth(giso)) nr = (isSpherical(owniso) ? E2Sph : E2E);
	else nr =  (isSpherical(owniso) ? Sph2Sph : ISO_MISMATCH);
      }
    }
    
    
    else if (isCylinder(giso) || isCylinder(owniso)) {
      RERR("general spherical coordinates not programmed yet");
    }
    
  
    else BUG;


    set_nri(GATTER, s, nr);

    
  } // for system
  RET_ENR(NOERROR);
}


#define PERR(X) { SPRINTF(ERRSTR, "'%.50s' : %.50s", param_name, X); return ERRORM;}
#define PERRX(ERR, X) { errorstring_type msg_1; errorMSG(ERR, msg_1);	\
    SPRINTF(ERRSTR, "'%.50s' : %.50s (%.50s)", param_name, X, msg_1); return ERRORM;}

int checkkappas(model *cov, bool errornull, errorstring_type ERRSTR){
  defn *C = DefList + COVNR; // nicht gatternr

  int i,nr, nc,
    kappas= C->kappas,
    *ncol = cov->ncol,
    *nrow = cov->nrow;
  char param_name[PARAMMAXCHAR]; // used in PERR

  for (i=0; i<kappas; i++) {
    model *ks = cov->kappasub[i];     
    if (SortOf(cov, i, 0, 0, original) == DONOTVERIFYPARAM) {
      if (ks != NULL) BUG;
      continue;
    }
    STRCPY(param_name, 
	   cov->ownkappanames != NULL && cov->ownkappanames[i]!=NULL 
	   ? cov->ownkappanames[i]
	   : C->kappanames[i]);
    
    if (ks != NULL) {
      if (isRandom(ks)) { // allgemeiner TypeConsistent(RandomType, ... )??
	if (!allowsRandomParam(cov,i)) PERR("argument must be deterministic");
	cov->randomkappa = true;
	int err, len;       
	
	nr = nrow[i];
	nc = ncol[i];
	if (nr==SIZE_NOT_DETERMINED || nc==SIZE_NOT_DETERMINED) 
	  C->kappasize(i, cov, &nr, &nc); 
	if (nr==SIZE_NOT_DETERMINED || nc==SIZE_NOT_DETERMINED) {
	  int d;
	  for (d=9; d>=1; d--) { // to do -- 9 sieht nicht gut aus.
	    //	    printf("check_r %s : vdim = %d\n", NAME(ks), d);
	    err = CHECK_R(ks, d);
	    if (err != NOERROR && ks->vdim[0] > 0 && ks->vdim[0] != d) {
	      //    printf("zweites check_r %s : vdim = %d\n", NAME(ks), d);
	      err = CHECK_R(ks, ks->vdim[0]);
	      d = 0;
	    }
	    if (err == NOERROR) {
	      nr = ks->vdim[0];
	      nc = ks->vdim[1];
	      len = nr * nc;
	      if (MODELNR(ks) == DISTRIBUTION) {
		if (PARAMINT(ks, DISTR_NROW) != NULL) 
		  nr = PARAM0INT(ks, DISTR_NROW);
		if (PARAMINT(ks, DISTR_NCOL) != NULL) 
		  nc = PARAM0INT(ks, DISTR_NCOL);
	      }
	      break;
	    }
	  }
	  if (err != NOERROR) RET_ENR(err);
	} else {
	  //if (nr==SIZE_NOT_DETERMINED || nc==SIZE_NOT_DETERMINED)
	  // PERR("size of random parameter could not be determined -- please give the size explicitely");	
	  len = nr * nc;	

	  //	  PMI(ks);
	  // printf("drittes check_r %s : vdim = %d\n", NAME(ks), len);

	  if ((err = CHECK_R(ks, len)) != NOERROR) {
	    PERRX(err, "random parameter not well defined");
	  }
	}
	
	if ( (ks->vdim[0] != nr || ks->vdim[1] != nc) &&
	     (ks->vdim[0] != len || ks->vdim[1] != 1) )
	  PERR("required size of random parameter does not match the model");
     
	if (cov->Sgen == NULL) NEW_STORAGE(gen);
	

  	if (PisNULL(i)) {
	  PALLOC(i, nr, nc);
	}

	//printf("init_random %s : vdim = %d\n", NAME(ks), 5);
	if (!cov->initialised && 
	    (err = INIT_RANDOM(ks, 0, cov->Sgen, P(i))) != NOERROR) { 
	  // wird spaeter
	  // gegebememfalls nochmals initialisiert mit richtigen moments
	  //X4
	  PERRX(err, "random parameter cannot be initialized");	  
	}	
      } else { // not random, e.g. Aniso
	if (!allowsShapeParam(cov, i)) 
	  PERR("argument may not be an arbitrary function");
	// no automatic check possible
	// if (!ks->checked) BUG; // fails
	// ks->checked = false; // fails as well
      }
    } // end ks != NULL (function given)

    if (PisNULL(i)) {
      if (errornull) { PERR("unset"); }
      else continue;    
    }
    
 
    C->kappasize(i, cov, &nr, &nc);   
    if ( (nc < 1 && nc != SIZE_NOT_DETERMINED) || 
	 (nr < 1 && nr != SIZE_NOT_DETERMINED)) { 
      BUG;
    }
    
    if (nc == 1 && ncol[i] != nc) {
      if (nr == 1 && nrow[i] != 1) PERR("must be a scalar");

      //      PERR("must be a vector, not a matrix");
    }
    if (nc > 1 && ncol[i] == 1) PERR("parameter must be a (larger) matrix");
    
    if ((nc > 0 && ncol[i] != nc) || (nr > 0 && nrow[i] != nr)) {
      
      // nc==0, nr==0 is coded as SIZE_NOT_DETERMINED
      char info[255], info2[255];
      SPRINTF(info, "not of the required size: (%d, %d) instead of (",
	      nrow[i], ncol[i]);
      if (nr!=SIZE_NOT_DETERMINED) SPRINTF(info2, "%.50s%d, ", info, nr);
      else SPRINTF(info2, "%.50sundetermined, ", info);
      if (nc!=SIZE_NOT_DETERMINED) SPRINTF(info, "%.50s%d)", info2, nc);
      else SPRINTF(info, "%.50sundetermined)", info2);
      // printf("info = %s\n", info);
      //      crash();      APMI0(cov);
      PERR(info);
    }
    // if nc==0 (nr==0) then this is undetermined.
  } // for i < kappas

  RET_ENR(NOERROR);
}


int checkkappas(model *cov){
  return checkkappas(cov, true, cov->err_msg);// #define RETURN
}


int checkkappas(model *cov, bool on){
  return checkkappas(cov, on, cov->err_msg);// #define RETURN
}




int checkDims(model *cov, int vdim0, int vdim1, errorstring_type ERRSTR) {
  model *calling = cov->calling;
  int n = OWNSYSTEMS;
  system_type *def = DEF;
  assert(n == SYSTEMS(def));
  for (int s=0; s<n; s++) {
    if (MAXDIM(def, s) >= 0 && MAXDIM(OWN, s) > MAXDIM(def, s)) {
      set_maxdim(OWN, s, MAXDIM(def, s));
      assert(MAXDIM(def, s) > 0);
    }
  }
  if (VDIM0 <= 0 || VDIM1 <= 0) RET_ENR(ERRORBADVDIM);
  
  if ((vdim0 > 0 && VDIM0 != vdim0) || (vdim1 > 0 && VDIM1 != vdim1))
    RERR7("multivariate dimension (of submodel '%.50s'), which is %d x %d, does not match the specification of '%.50s', which is %d x %d and is required by '%.50s'",
	  NICK(cov), VDIM0, VDIM1, NAME(cov), vdim0, vdim1, 
	  calling == NULL ? "-- none --" :  NAME(calling));
 
  RET_ENR(NOERROR);
}


int check_within_range(model *cov, bool NAOK, errorstring_type ERRSTR) {
  // sets also maxdim and finiterange in cov !

  defn *C = DefList + COVNR; //nicht gatternr
  int len, 
    err = NOERROR,
    k = 0, 
    i = 0,   // compiler dummy values
    kappas = C->kappas;
  range_type range;
  char Msg[255];
  SEXPTYPE *type = C->kappatype;
  rangefct getrange = C->range;
  double min, max,
    value= RF_NA;

  if (GLOBAL_UTILS->basic.skipchecks) RET_ENR(NOERROR);

  getrange(cov, &range); 

  if (!maxdim_ok(cov)) {
    int s = -maxdim_ok(cov);
    RERR3("Max. dimension in '%.50s' is %d. Got %d", NAME(cov), OWNMAXDIM(s), OWNLOGDIM(s));
  }

  for (i=0; i<kappas; i++) {
    if (cov->kappasub[i] != NULL) continue;

    if (!STRCMP(C->kappanames[i], FREEVARIABLE)) {
      // i.e. equal
      if (PisNULL(i)) continue;
    }
    if (type[i] >= LISTOF) {
      // to simplify things -- otherwise in simu.cc
      // code must also be changed.
      assert(range.min[i] == RF_NEGINF && range.max[i]==RF_INF);
      continue;
    }
    // full range !!
    
    len = cov->ncol[i] * cov->nrow[i];
    min = range.min[i];
    max = range.max[i];
    
    err = CERRORM; 
    for (k=0; k<len; k++) {
      if (type[i] == REALSXP) value = P(i)[k];
      else if (type[i] == INTSXP)
	value = PINT(i)[k] == NA_INTEGER ? RF_NA : (double) PINT(i)[k];
      else if (isRObject(type[i]) || type[i]==STRSXP) continue;
      else { STRCPY(Msg, "is not finite"); goto ErrorHandling; }
      if (ISNAN(value)) {
	if (NAOK) continue;
	else { STRCPY(Msg, "is not finite"); goto ErrorHandling; }
      }
      
      if (range.openmin[i] && value <= min) { 
	addmsg(value, ">", min, Msg);
	goto ErrorHandling;
      } else if (value < min) {
	addmsg(value, ">=", min, Msg); 
        goto ErrorHandling;
      }
      if (range.openmax[i] && value >= max) { 
	addmsg(value, "<", max, Msg); 
        goto ErrorHandling;
      } else if (value > max) { 
	addmsg(value, "<=", max, Msg);
	goto ErrorHandling;
      }
    }
  } // kappas

  RET_ENR(NOERROR);

 ErrorHandling:
  if (PL>PL_ERRORS) {
    PRINTF("error in check_witih_range (%s): %s %s(%d) err=%d ('%s' does not hold.)\n",
	   __FILE__, C->name, C->kappanames[i], i, err, Msg);
  }
  //crash();
   /*
     if (err == ERRORUNSPECIFIED)
     RERR4("%.50s[%d] = %.50s does not hold (logicaldim(0)=%d).",
 	     C->kappanames[i], k+1, Msg, OWNLOGDIM(0)); // + r eturn
   */
   RERR3("%.50s[%d]=%.50s does not hold.", C->kappanames[i], k+1, Msg);
}



int check_recursive_range(model *cov, bool NAOK) {  // FOR MLE ONLY !!
    /*
      called by "Set And GetModelInfo"
     */
  int i, err,
    kappa = DefList[COVNR].kappas;
  KEY_type *KT = cov->base;
  
  SPRINTF(KT->error_loc, "'%.50s'", NICK(cov));
  if ((err = check_within_range(cov, NAOK, cov->err_msg)) != NOERROR)
    RET_ENR(err);
  for (i=0; i<kappa; i++) {
     if (cov->kappasub[i] != NULL &&
	 (err = check_recursive_range(cov->kappasub[i], NAOK))
	 != NOERROR) RET_ENR(err);
  }
  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL && 
	(err = check_recursive_range(cov->sub[i], NAOK))
	!= NOERROR) RET_ENR(err);
  }
  RETURN_NOERROR;
}
