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



#ifndef RF_MACHINE_DEBUG
#define RF_MACHINE_DEBUG 1

#undef COV_DELETE
#undef NICK
#undef INIT_RANDOM
#undef CHECK
#undef CHECK_NO_TRAFO
#undef CHECK_VDIM
#undef CHECKPD2ND


#undef CHECK_GEN
#undef CHECK_ONLY
#undef CHECK_NOPASS
#undef CHECK_PASSTF
#undef CHECK_PASSFRAME
#undef CHECK_PASSTYPE
#undef CHECK_R

#undef INIT
#undef REINIT
#undef INIT_RANDOM
#undef STRUCT
#undef DEBUG


/*
//#define CHECKFOR(C) //
void checkfor(model *cov);
#ifdef CHECKFOR
#undef CHECKFOR
#endif
#ifdef CHECKFOR_BACK
#undef CHECKFOR_BACK
#endif

#define CHECKFOR(C)							\
  PRINTF("CHECKFOR: '%.50s', line %d", __FILE__, __LINE__); checkfor(C)
#define CHECKFOR_BACK(C)						\
  PRINTF("CHECKFORBACK: '%.50s', line %d", __FILE__, __LINE__);checkfor(C)
*/

#define DEBUG(A) {							\
  /* //  printf("at '%.50s', line %d, %.50s\n", __FILE__, __LINE__, #A);    */ \
  A;									\
  /*  PRINTF("done (debug): '%.50s', line %d, %.50s \n", __FILE__, __LINE__, #A);*/ \
  PRINTF("done (debug): '%.50s', line %d\n", __FILE__, __LINE__);\
  }


#define COV_DELETE(Cov) {					\
    PRINTF("\nCOV_DELETE: '%.50s', line %d", __FILE__, __LINE__);	\
    COV_DELETE_(Cov);}

#define NICK(Cov) (DefList[MODELNR(Cov)].nick)

#define XX(C) assert((C)->simu.expected_number_simu >= 0 ||		\
		     __extension__({DOPRINTF("Start.\n"); false;}))
#define YY(C) assert((C)->simu.expected_number_simu >= 0 ||		\
		     __extension__ ({DOPRINTF("End.\n"); false;}))

#define LLPRINT(SIGN, cov, Z)						\
    if (!leading_spaces(cov, SIGN)) { } else 				\
      PRINTF("(%.50s, %.50s, line %d : %.50s)\n", Z, __FILE__, __LINE__, NAME(cov)) 

#define CHECKSIGN "_"
#define INITSIGN "_"
#define STRUCTSIGN "_"

#define ANYSTART(C, MESS, SIGN)  __extension__({	\
  LLPRINT(SIGN, C, MESS);			\
  XX(C)
    
#define CHECKSTART(C) ANYSTART(C,"CHECK", CHECKSIGN)

#define ANYEND(C,MESS, SIGN) 						\
  YY(C);								\
  if (_x==NOERROR){LLPRINT(SIGN, C, " DONE");}				\
  else { LLPRINT(SIGN, C, "#MESS FAILED"); errorstring_type msg_0;		\
    errorMSG(_x, cov->err_msg, cov->base, msg_0); PRINTF("%.50s\n", msg_0);} \
  _x;})

#define CHECKEND(C) ANYEND(C,"CHECK", CHECKSIGN)

#define CHECK(C,T,X,type,D,I,V,R) CHECKSTART(C)				\
  assert((type)!=RandomType);						\
  int _x = check2X(C,T,X,type,D,I,V,R);					\
  CHECKEND(C)

#define CHECK_NO_TRAFO(C,T,X,type,D,I,V,R) CHECKSTART(C)		\
      assert((type)!=RandomType);					\
      int _x = check2Xnotrafo(C,T,X,type,D,I,V,R);			\
      CHECKEND(C)

#define CHECK_VDIM(C,T,X,type,D,I,V0,V1,R) CHECKSTART(C)		\
      int _x = check2X(C,T,X,type,D,I,V0,V1,R,true);		\
      CHECKEND(C)

#define CHECKPD2ND(C,V,R) CHECKSTART(C)				\
      int _x = CheckPD2ND(C,V,R);YY(C);				\
      CHECKEND(C)

#define CHECK_GEN(C,V0, V1, R, CT) CHECKSTART(C)	\
  int _x = check2X(C,V0, V1, R, CT);		\
  CHECKEND(C)

#define CHECK_ONLY(C) CHECKSTART(C)	\
  int _x = check2X(C,(C)->vdim[0],(C)->vdim[1],(C)->frame,false);	\
  CHECKEND(C)

#define CHECK_NOPASS(C) CHECKSTART(C)				\
  int _x = check2passframe(C, OWN, VDIM0, VDIM1, cov->frame);		\
  CHECKEND(C)

#define CHECK_PASSTF(C,T,V,R) CHECKSTART(C)	\
  int _x = check2passTF(C, OWN, T, V, R);		\
  CHECKEND(C)

#define CHECK_PASSFRAME(C,R) CHECKSTART(C)	\
  int _x = check2passframe(C, OWN, VDIM0, VDIM1, R);		\
  CHECKEND(C)

#define CHECK_PASSTYPE(C,T) CHECKSTART(C)	\
  int _x = check2passtype(C, OWN, T, VDIM0, VDIM1, cov->frame);		\
  CHECKEND(C)

#define CHECK_R(C, vdim) CHECKSTART(C)	\
  int _x = check2X(C, vdim, vdim, RandomType, KERNEL, CARTESIAN_COORD, \
		   vdim, 1, RandomType, true);			       \
  CHECKEND(C)			       


#define INIT(C, Moments, S) ANYSTART(C,"INIT", INITSIGN);		\
      int _x = INIT_intern(C, Moments, S);			\
      YY(C);							\
      if (_x==NOERROR){LLPRINT(STRUCTSIGN, C, "INIT DONE");}	\
      else {LLPRINT(STRUCTSIGN, C, "INIT FAILED");}		\
      _x;})

#define INITEND(C)	ANYEND(C,"INIT", INITSIGN)

#define REINIT(C, Moments, S) ANYSTART(C,"REINIT", INITSIGN);	\
  int _x = REINIT_intern(C, Moments, S);	\
  INITEND(C)


#define INIT_RANDOM(C, Moments, S, P) ANYSTART(C,"INITRANDOM", INITSIGN); \
      int _x = INIT_RANDOM_intern(C, Moments, S, P);\
      INITEND(C)

#undef STRUCT
#define STRUCT(C, NM) ANYSTART(C,"STRUCT", STRUCTSIGN);		 \
	   ASSERT_GATTER(C);					 \
	   int _x = DefList[FIRSTGATTER].Struct(C, NM);		 \
	   INITEND(C)

   /*  printf("//%ld %ld %d %d %d idx=%d\n", (TO)->px[IDX], (FROM)->px[IDX], \
       (FROM)->nrow[IDX], (FROM)->ncol[IDX],				\
       DefList[MODELNR(Cov)].kappatype[IDX]==REALSXP ? sizeof(double) :	\
       DefList[MODELNR(Cov)].kappatype[IDX]==INTSXP ? sizeof(int) :	\
       -1, IDX);							\
		
   */

#ifndef XXX_ABC_FGT
#undef DEBUGINFOERR 
#define DEBUGINFOERR {						\
    errorstring_type dummy_; STRCPY(dummy_, cov->err_msg); 	\
    SPRINTF(cov->err_msg, "%.50s (%.50s, line %d)\n", dummy_, __FILE__, __LINE__); \
    PRINTF("note: %.50s\n", cov->err_msg);					\
  }
#endif



#endif
