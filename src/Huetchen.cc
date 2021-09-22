/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 simulation of max-stable random fields

 Copyright (C) 2001 -- 2017 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#include <stdio.h>
#include "def.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "shape.h"
#include "operator.h"
#include "Processes.h"
#include "primitive.h" 
#include "families.h"

 
#define POISSON_INTENSITY 0
#define RANDOMCOIN_INTENSITY (COMMON_GAUSS + 1) 



#define ZHOU_VOXEL 0
#define ZHOU_CORNER 1
#define ZHOU_SIDES (ZHOU_CORNER + 1)

#define ZHOU_2D_SIDE1 ZHOU_SIDES
#define ZHOU_2D_SIDE2 (ZHOU_SIDES + 1)

#define ZHOU_3D_SIDE1 ZHOU_SIDES
#define ZHOU_3D_SIDE2 (ZHOU_SIDES + 1)
#define ZHOU_3D_SIDE3 (ZHOU_SIDES + 2)
#define ZHOU_3D_EDGE1 (ZHOU_SIDES + 3)
#define ZHOU_3D_EDGE2 (ZHOU_SIDES + 4)
#define ZHOU_3D_EDGE3 (ZHOU_SIDES + 5)


void closest(double *x, model *cov, double *delta) {
  location_type *loc = Loc(cov);
  int dim = OWNXDIM(0);
  assert(Getgrid(cov));

  for (int d=0; d<dim; d++) {
    double
      start = loc->xgr[d][XSTART],
      step = loc->xgr[d][XSTEP],
      idx = ROUND((x[d] - start) / step),
      maxidx = loc->xgr[d][XLENGTH] - 1.0;
    idx = idx < 0.0 ? 0.0 : idx > maxidx ? maxidx : idx;
    delta[d] = x[d] - (start + idx * step);
  }
}


void Zhou(double *x, model *cov, // pointshapetype
		     double *v) {
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  double w;
  TALLOC_X1(delta, OWNXDIM(0));
  NONSTATCOV(x, cov->q, shape, v);
  closest(cov->q, cov, delta);
  VTLG_D(delta, pts, &w);
  *v /= w;
  //printf("x=%f %f v=%f\n", x[0], x[1], *v);
 // PMI(cov);
  
  /// if (*v > (1  + 1e-12) * pts->mpp.unnormedmass) printf("violating v=%10e\n", *v);
  assert(*v < (1 + 1e-12) * pts->mpp.unnormedmass);
    END_TALLOC_X1;
}

void logZhou(double *x, model *cov, double *v, double *Sign) { 
 model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  double w;
  TALLOC_X1(delta, OWNXDIM(0));
  LOGNONSTATCOV(x, cov->q, shape, v, Sign);
  closest(cov->q, cov, delta);
  VTLG_DLOG(delta, pts, &w);

  *v -= w; // + LOG(pts->mpp.unnormedmass);
  assert(*v - cov->Spgs->log_density < 1e-12);
  //printf("violating v=%10e > 0\n", *v);

  END_TALLOC_X1;
}

int check_Zhou(model *cov) {
  ASSERT_ONESYSTEM;
  ASSERT_UNREDUCED;
 model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  location_type *loc = Loc(cov);
  int err,
    dim = OWNLOGDIM(0);
 
  ASSERT_CARTESIAN;
  if (loc->Time) SERR("Time component not allowed yet"); // todo
   
  kdefault(cov, ZHOU_RATIO, GLOBAL.extreme.density_ratio); 
  kdefault(cov, ZHOU_FLATHULL, (int) GLOBAL.extreme.flathull);
  kdefault(cov, ZHOU_INFTY_SMALL, P0INT(ZHOU_FLATHULL) != (int) False);
  kdefault(cov, ZHOU_NORMED, true);
  kdefault(cov, ZHOU_ISOTROPIC, true);

  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);
  
  if (cov->q == NULL) QALLOC(dim);
  
  assert(!isProcess(cov));
  assert(hasSmithFrame(cov));
 
  // isotropy is the property -- user should not define unisotropic functions!
  // so !isotropic only for special, internal definitions, e.g. 
  // strokorbPoly
  // int err2 = NOERROR;
  //  if (P0INT(ZHOU_ISOTROPIC)) 
  //    err2 = CHECK(shape, dim, dim, ShapeType, XONLY, ISOTROPIC,
  //		 SCALAR, frame);


  // but it must be treated in any case as if it was non-isotropic
  if ((err = CHECK(shape, dim, dim, ShapeType, XONLY, CARTESIAN_COORD,
		   SCALAR, cov->frame)) != NOERROR) {
    if (P0INT(ZHOU_ISOTROPIC)) BUG;
    XERR(err);
    RETURN_ERR(err);
  }

  setbackward(cov, shape);

  //if (shape->randomkappa) RETURN_ERR(ERRORNOTPROGRAMMED); // Dichte muss nicht normiert sein und stattdessen durch das mittlere Volumen dividiert werden.

  if (pts != NULL) {
    if ((err = CHECK_R(pts, dim)) != NOERROR) RETURN_ERR(err);
  }
 
  EXTRA_STORAGE;
  RETURN_NOERROR;
  
}



int struct_Zhou(model *cov, model **newmodel){
  BUG; // should not be called at all! 2.2.19

  
  model *shape = cov->sub[PGS_FCT];
  //  location_type *loc = Loc(cov);
  int err = NOERROR;

  ASSERT_NEWMODEL_NULL;
  if (cov->Spgs != NULL)  pgs_DELETE(&(cov->Spgs), cov);

  if (!hasSmithFrame(shape)) ILLEGAL_FRAME;
  if (cov->sub[PGS_LOC] == NULL) {
     if ((err = STRUCT(shape, cov->sub + PGS_LOC)) != NOERROR)
      RETURN_ERR(err);

     if (cov->sub[PGS_LOC] == NULL) 
      SERR1("no intensity found for '%.50s'", NICK(shape));

    assert(cov->sub[PGS_LOC]->calling != NULL);
  }
  RETURN_NOERROR;
}



int calculate_mass_maxstable(model *cov) {// = pointshapetype

  pgs_storage *pgs = cov->Spgs;
  location_type *loc = Loc(cov);
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  double voxels, value_orig,
    *single = pgs->single, // out
    *total = pgs->total,   // out
    *halfstepvector = pgs->halfstepvector, //out
    *x = pgs->x;
  assert(x != NULL);
  int d,
    dim = OWNLOGDIM(0);
  usr_bool flathull = (usr_bool) P0INT(ZHOU_FLATHULL);
 
  if (hasPoissonFrame(shape)){
    BUG;
    // to do : felix's paper; derzeit einfach wie max-stable
    //   if (pts->nr !=    SERR("currently, only simple domain fields are allowed that are based on Poisson processes");

  } // else {

  VTLG_D(ZERO(pts), pts, &value_orig);

  //  PMI(cov, "calculate max");
  //  printf("%.50s %d %d %10g\n",NICK(cov), ZHOU_FLATHULL, P0INT(ZHOU_FLATHULL), value_orig );// APMI(cov); assert(false);
  
  // VOXEL 
  for (d=0; d<dim; d++) halfstepvector[d] = 0.5 * loc->xgr[d][XSTEP];

  if (flathull == Nan) {
    // flathull == im Kerngebiet wird konstant simuliert; ausserhalb dann
    // abfallend
    if (loc->grid) {
      double v;
      VTLG_D(halfstepvector, pts, &v);
      double threshold = 
	value_orig == RF_INF//&& cov->p[ZHOU_RATIO][0] == 0.0 
	? RF_INF
	: value_orig * P0(ZHOU_RATIO);
      pgs->flathull = (usr_bool) (threshold < v && !cov->randomkappa); // sicherheitshalber

      //printf("thres %10g %10g %10g\n", threshold, v, value_orig);

    } else { // not grid
      BUG;
      // long term project: to do
      pgs->flathull = True;
      /* Bei letzterem und nicht Gitter wird dann wird eine maximale 
	 Kaestchengroesse bestimmt so dass gerade noch beruehrt wird
	 und der Rest ausserhalb aller Kaestchen auf einen konstanten 
	 wert gesetzt, etwas niedriger als der hoechste Grenzwert. Was 
	 "hoch" bedeutet ist hier unklar. 
      */
    }      
  } else pgs->flathull = flathull;

  if (pgs->flathull == True) {
    if (P0INT(ZHOU_INFTY_SMALL)) {
      SERR2("'%.50s' and '%.50s' may not be positive at the same time",
	    KNAME(ZHOU_FLATHULL), KNAME(ZHOU_INFTY_SMALL));
    }
    single[ZHOU_VOXEL] = value_orig;
    for (d=0; d<dim; d++) single[ZHOU_VOXEL] *= 2.0 * halfstepvector[d];
  } else {
    //    APMI(cov);
    //       printf("two sided %10g %10g %10g\n", halfstepvector[0], halfstepvector[1], halfstepvector[2]); //assert(false);
    VTLG_P2SIDED(NULL, halfstepvector, pts, single + ZHOU_VOXEL);
    //     printf("two sided done %10e %10g %.50s\n", single[ZHOU_VOXEL], halfstepvector[0], NICK(pts)); //assert(false);
  }
  voxels = 1.0;
  for (d=0; d<dim; d++) voxels *= loc->xgr[d][XLENGTH] - 1.0;  
  total[ZHOU_VOXEL] = single[ZHOU_VOXEL] * voxels;
  
  // CORNER
  for (d=0; d<dim; d++) x[d] = RF_INF;  
  VTLG_P2SIDED(NULL, x, pts, single + ZHOU_CORNER);

  //  printf("%d %10g\n", !P0INT(ZHOU_NORMED), single[ZHOU_CORNER]);

  assert(!P0INT(ZHOU_NORMED) || single[ZHOU_CORNER] == 1.0); // nur wenn normiert
  total[ZHOU_CORNER] = single[ZHOU_CORNER] + total[ZHOU_VOXEL];

  if (dim >= 2) {
    for (d=0; d<dim; d++) {
      int nr = 1;
      if (pgs->flathull == True) for (int i=0; i<dim; i++) x[i] = 0.0;
      else for (int i=0; i<dim; i++) x[i] = halfstepvector[i];
      x[d] = RF_INF; // !! changed
      VTLG_P2SIDED(NULL, x, pts, single + ZHOU_SIDES + d);

      for (int i=0; i<dim; i++) {
	if (R_FINITE(x[i])) {
	  if (pgs->flathull == True) single[ZHOU_SIDES + d] *= loc->xgr[i][XSTEP]; 
	  nr *= loc->xgr[i][XLENGTH] - 1.0;
	}
      }
      //
      //PMI(cov);
      //  printf("  single d=%d %d x=(%10g %10g %10g) val=%10g %d\n", d, ZHOU_SIDES + d, 
      //     x[0], x[1], x[2], single[ZHOU_SIDES + d], nr); assert(false);

      total[ZHOU_SIDES + d] =
	total[ZHOU_SIDES + d - 1] + nr * single[ZHOU_SIDES + d];
    }

    if (dim == 3) { // Kanten
      for (d=0; d<dim; d++) {
	int nr = 1.0;
	for (int i=0; i<dim; i++) x[i] = RF_INF;
	x[d] = pgs->flathull == True? 0.0 : halfstepvector[d];
	VTLG_P2SIDED(NULL, x, pts, single + ZHOU_SIDES + dim + d);
	if (pgs->flathull == True)
	  single[ZHOU_SIDES + dim + d] *= loc->xgr[d][XSTEP]; 
	nr *= loc->xgr[d][XLENGTH] - 1.0;
	
	//	printf("  single d=%d %d x=(%10g %10g %10g) val=%10g %d\n", d, ZHOU_SIDES + dim + d,  x[0], x[1], x[2], single[ZHOU_SIDES + dim + d], nr); //assert(false);

	total[ZHOU_SIDES + dim + d] = 
	  total[ZHOU_SIDES + dim + d - 1] + nr * single[ZHOU_SIDES + dim + d];
      }
    } else if (dim > 3) BUG;
  }

  pgs->totalmass = total[pgs->size - 1];
  if (!R_FINITE(pgs->totalmass)) {
    // int i; for (i=0; i<pgs->size; i++) printf("i=%d total=%10g\n", i, total[i]);
    SERR("the total intensity mass is not finite");
  }

  RETURN_NOERROR;
}


model *prunecov(model *newmodel, model *cov) {
  //printf("names= %.50s %.50s\n", NICK(newmodel->calling), NICK(cov->calling));
  model
    **sub = NULL,
    *next = NULL,
    *calling = cov->calling;
  if (calling == newmodel->calling) {
    return newmodel;
  }
  if (calling == NULL) BUG;
  prunecov(newmodel, calling);
  if (calling->key == cov) sub = &(newmodel->key);
  else if (calling->sub[0] == cov) sub = newmodel->sub + 0;
  else if (calling->sub[1] == cov) sub = newmodel->sub + 1;
  else BUG;
  next = *sub;
  *sub = NULL;
  COV_DELETE(&newmodel, cov);
  return next;
}

int complete_copy(model **newmodel, model *covIn) {
  model *cov = NULL,
    *start = covIn;
  int err;
  while (start->calling != NULL) start = start->calling;
 if (!equalsnowInterface(start)) BUG; // nur zur Sicherheit, derzeit; kann erweitert werden
  if (start == covIn) BUG;
  cov = start->key != NULL ? start->key : start->sub[0];
  if (!equalsnowProcess(cov)) BUG;
  
  if ((err = covcpy(newmodel, cov)) != NOERROR) RETURN_ERR(err);
  SET_CALLING(*newmodel, covIn);
  Types frame = cov->frame;
  ASSERT_ONESYSTEM;
  if ((err = CHECK(*newmodel, OWNLOGDIM(0), PREVXDIM(0), OWNTYPE(0),
		   PREVDOM(0), PREVISO(0), cov->vdim, frame))
      //      if ((err = CHECK(*newmodel, cov->tsdim, cov->xdimcov, cov->typus, 
      //	   cov->domcov, cov->isocov, cov->vdim, frame)) 
      != NOERROR)  RETURN_ERR(err);  
   if ((err = STRUCT(*newmodel, NULL)) != NOERROR) { RETURN_ERR(err); }
  if (!(*newmodel)->initialised) {
   if ((err =  CHECK(*newmodel, OWNLOGDIM(0), PREVXDIM(0), OWNTYPE(0), 
		      PREVDOM(0), PREVISO(0), cov->vdim, frame)) 
       //    if ((err =  CHECK(*newmodel, cov->tsdim, cov->xdimcov, cov->typus, 
       //	      cov->domcov, cov->isocov, cov->vdim, frame)) 
	!= NOERROR) RETURN_ERR(err);  
 
    NEW_COV_STORAGE(*newmodel, gen);
   //APMI(*newmodel);
    if ((err = INIT(*newmodel, 0, covIn->Sgen)) != NOERROR) {
      //APMI(covIn); // !!! ?? hier weitermachen
      RETURN_ERR(err); 
    }
  }
 
  SET_CALLING(*newmodel, start);
  *newmodel = prunecov(*newmodel, covIn);
  SET_CALLING_NULL(*newmodel, start);
  RETURN_NOERROR;
}

int init_Zhou(model *cov, gen_storage *S) {// = pointshapetype
  ASSERT_ONESYSTEM;
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  defn *Cshape = DefList + MODELNR(shape);
  location_type *loc = Loc(cov);
  int 
    err = NOERROR,
    dim = XDIM(PREVSYSOF(shape), 0);
  pgs_storage *pgs = cov->Spgs;
  bool 
    grid = Loc(cov)->grid,
    pgsnull = pgs == NULL;
  if (!grid) { // todo
    SERR("non-grid not programmed yet"); // meiste da; fehlt noch
    // welche Vektoren ich wirklich brauche
  }

  assert(hasSmithFrame(cov));
  assert(cov->sub[PGS_LOC] != NULL);
  assert(dim == XDIM(PREVSYSOF(cov->sub[PGS_LOC]), 0));
  assert(dim == Loc(cov)->timespacedim);
  
  if (Cshape->inverse == ErrInverse)//&&Cshape->nonstat_inverse==Err CovNonstat)
    SERR1("support of the model is unknown. Use '%.50s' to determine the support",
	  DefList[TRUNCSUPPORT].nick);
 
  if (pgsnull) {
    if ((err = alloc_pgs(cov)) != NOERROR) RETURN_ERR(err); // insb pgs->x
    pgs = cov->Spgs;
    
    if ((pgs->v = (double*) CALLOC(dim, sizeof(double))) == NULL ||
	(pgs->y = (double*) CALLOC(dim, sizeof(double))) == NULL
	) RETURN_ERR(ERRORMEMORYALLOCATION);

    pgs->zhou_c = 1.0;
    pgs->sum_zhou_c = pgs->sq_zhou_c = 0.0;
    pgs->n_zhou_c = 0;
  }
  pgs->size = loc->totalpoints;
  if (grid) pgs->size = intpow(2, dim);
 
  // selbst wenn zufaelliger Shape: 1x laufen lassen, ob 
  // Fehler auftauchen. Unter "Do" lassen sie sich nicht mehr
  // schoen abfangen.


  if ((err = INIT(pts, 0, S)) != NOERROR) RETURN_ERR(err);// , 0, S statt ,1,S
  if ((err = INIT(shape, cov->mpp.moments, S)) != NOERROR)  RETURN_ERR(err);
  assert(cov->mpp.moments >= 1);
  assert(shape->mpp.moments >= 1);

  
  if (pgsnull &&
      ((pgs->single = (double*) CALLOC(pgs->size, sizeof(double)))==NULL ||
       (pgs->total = (double*) CALLOC(pgs->size, sizeof(double))) ==NULL ||
       (pgs->halfstepvector = (double*) CALLOC(dim, sizeof(double)))==NULL
       ))
    RETURN_ERR(ERRORMEMORYALLOCATION);
  
  cov->mpp.maxheights[0] = 1.0; 
  if (P0INT(ZHOU_NORMED)) {
    if (R_FINITE(pts->mpp.unnormedmass)) {
      
      //      printf("masses shape=%10g pts=%10g\n", shape->mpp.mMplus[1], pts->mpp.unnormedmass);
      
      //APMI(cov);
    } else {
      pts->mpp.unnormedmass = shape->mpp.mMplus[1];
      if (!R_FINITE(pts->mpp.unnormedmass))
	SERR("not possible to be simulated at the moment");
      //	cov->mpp.maxheights[0] =  shape->mpp.maxheights[0] * pts->mpp.maxheights[0];
    }
  } else {
    if (R_FINITE(pts->mpp.maxheights[0])) {
      // eigentlich ein 1.0 / pts->mpp.minheight
      
      BUG;
      
      cov->mpp.maxheights[0] =
	shape->mpp.maxheights[0] / pts->mpp.maxheights[0];
    } else {
      BUG;
    }
  }
  
  // f_raw,rectangular must be always >= the shape function.
  // this is the fundamental condition
  //  APMI(cov);
  if (!R_FINITE(cov->mpp.maxheights[0])) BUG;
  
  assert(cov->randomkappa || P0INT(ZHOU_NORMED));
  assert(cov->randomkappa == shape->randomkappa);
  if ((pgs->estimated_zhou_c =cov->randomkappa)) {
    pgs->zhou_c = RF_NA;
    if (pgs->cov == NULL) {
      model *start=cov->calling;
      if (start == NULL) BUG;	
      //  APMI(start);
      while (start->calling != NULL && MODELNR(start) != ZHOU)
	start = start->calling;	
      if (MODELNR(start) != ZHOU) {
	if ((err=complete_copy(&(pgs->cov), cov)) != NOERROR) RETURN_ERR(err);
	SET_CALLING(pgs->cov, cov->calling);
	pgs->cov->Spgs->cov = cov;
      }
    }
    pgs->old_zhou = 0.0;
  } else {
    if ((err = calculate_mass_maxstable(cov)) != NOERROR) RETURN_ERR(err);
    pgs->zhou_c =
      pts->mpp.unnormedmass * pgs->totalmass / shape->mpp.mMplus[1];
    //     printf("%10g %10g %10g %10g\n", pgs->zhou_c, pts->mpp.unnormedmass, pgs->totalmass, shape->mpp.mMplus[1]);assert(false);
  }
  

  if (DefList[MODELNR(shape)].nonstat_inverse == ErrInverseNonstat) {
    if (MODELNR(pts) != RECTANGULAR) {
      //       PMI(pts);
      // APMI(shape);
      warning("Inverse of shape function cannot be determined. Simulation speed  might be heavily decreased.");
    }
  }
  
  for (int i=0; i<=cov->mpp.moments; i++) {
    cov->mpp.mM[i] = shape->mpp.mM[i] * pts->mpp.mMplus[0];
    cov->mpp.mMplus[i] = shape->mpp.mMplus[i] * pts->mpp.mMplus[0];
  }
  
  ReturnOtherField(cov, shape); // see also do for fieldreturn
 
  RETURN_ERR(err);
}

int DrawCathegory(int size, double *single, 
		  double *total, 
		  bool calculate_elements, int *elmts) { 

  int i;
  double mass;
  mass = UNIFORM_RANDOM * total[size - 1];
  if (calculate_elements) {
    //printf("mass=%10g\n", mass);
    
    for (i = 0; mass > total[i]; i++);
    if (i > 0) mass -= total[i-1];
    *elmts = FLOOR(mass / single[i]);
    //printf("draw: i=%d tot=%10g mass=%10g sing=%10g elm=%d\n", i, total[i], mass, single[i], *elmts);
    return i;
  } else {
    // falls nicht calculate_elements, so sind letztere unwichtig,
    // d.h. zu single wird nur der Faktor 1.0 multipliziert.
    // dies ist aber der Fall nur falls !grid
    return CeilIndex(mass, total, size);
  }
}

#define SET_COV					 \
  pgs = cov->Spgs;				 \
  shape = cov->sub[PGS_FCT];			 \
  pts = cov->sub[PGS_LOC];			 \
  loc = Loc(cov);				 \
  single = pgs->single;				 \
  total = pgs->total;				 \
  halfstepvector = pgs->halfstepvector;		 \
  x = pgs->x;					 \
  v = pgs->v;					 \
  dim = XDIM(PREVSYSOF(shape), 0);		 \
  flathull = (bool) pgs->flathull;			 \
  assert(pgs != NULL);				 \
  assert(!cov->randomkappa ||			 \
	 (pgs->cov!=NULL && pgs->cov->Spgs!=NULL \
	  && pgs->cov->Spgs->cov==cov));

#define SWAP_PGS				\
  cov = cov->Spgs->cov;				\
  SET_COV

 
void do_pgs_maxstable(model *cov, gen_storage *S) {// = pointshapetype
  pgs_storage *pgs = NULL;
  model *shape = NULL,
    *pts = NULL;
  location_type *loc = NULL;
  int d, elmts, mcmc, dim, err;
  double 
    *single = NULL, // in
    *total = NULL,   // in
    *halfstepvector = NULL, // in
    *x = NULL, // dummy
    *v = NULL; // dummy
  bool flathull;
  assert(cov->calling != NULL);

  if (!cov->randomkappa) {
    SET_COV;
  } else {
    // assert(PARAM0(shape->sub[0], POWSCALE) == PARAM0(pts, LOC_SCALE));
    double cmaxDmu;
    SWAP_PGS;    
    for (mcmc=0; mcmc<GLOBAL.extreme.mcmc_zhou; mcmc++) {
      if ((err = REINIT(cov, cov->mpp.moments, S)) != NOERROR) BUG;
      DO(shape, S);     
      if (calculate_mass_maxstable(cov) != NOERROR) 
	ERR("unexpected error in 'do_Zhou' (maxstable)");

      cmaxDmu = pgs->totalmass / shape->mpp.mMplus[1];
      if (pgs->n_zhou_c < GLOBAL.extreme.max_n_zhou) {
	pgs->n_zhou_c++;
	pgs->sum_zhou_c += cmaxDmu;
	pgs->sq_zhou_c += cmaxDmu * cmaxDmu;
	// sumArea += getArea( shape->sub[0]->Spolygon->P);
	//if (pgs->n_zhou_c % 1000 == 0) printf("mean area = %10g\n", sumArea / (double) pgs->n_zhou_c);
      }

      if (false)
      printf("//cmaxDmu %4.1f; mass=%4.1f m1=%4.1e beta=%4.1f (%4.6f, %4.6f) (%4.6f, %4.6f) %4.3f %10e\n",	     
	     cmaxDmu, (double) pgs->totalmass, shape->mpp.mMplus[1], 
	     PARAM0(shape->sub[0], 0), 
	     shape->sub[0]->Spolygon->P->box0[0],
	     shape->sub[0]->Spolygon->P->box0[1],
	     shape->sub[0]->Spolygon->P->box1[0],
	     shape->sub[0]->Spolygon->P->box1[1],
	     101 * 101 *  shape->sub[0]->Spolygon->P->box1[0]* shape->sub[0]->Spolygon->P->box1[1] * 4,
	     101 * 101 *  shape->sub[0]->Spolygon->P->box1[0]* shape->sub[0]->Spolygon->P->box1[1] * 4 -  pgs->totalmass);

      //assert( shape->sub[0] != NULL);
      //assert(shape->sub[0]->Spolygon != NULL);
      //assert(shape->sub[0]	->Spolygon->P != NULL); 


      //assert(false);
      //    if (!(shape->mpp.mMplus[1] > 0.0001 || pgs->totalmass < 100)) APMI(cov);
      //    assert( pgs->totalmass <= 101 * 101 *  shape->sub[0]->Spolygon->P->box1[0]* shape->sub[0]->Spolygon->P->box1[1] * 4 + 1e-12);

      //PMI(cov);
      //assert( shape->sub[0]->Spolygon->P->box0[0] == PARAM(pts, 0)[0]);

      double old_zhou = pgs->cov->Spgs->old_zhou;
      if (old_zhou < cmaxDmu || UNIFORM_RANDOM * old_zhou < cmaxDmu) {
	//printf(". %10g %10g\n", old_zhou, cmaxDmu);
	pgs->old_zhou = cmaxDmu;
	SWAP_PGS;
      } //else printf("* %10g %10g\n", old_zhou, cmaxDmu);
    } // mcmc
    
    SWAP_PGS;
    model *calling=cov->calling;		
    assert(calling != NULL);				
    if (calling->key != NULL) calling->key = cov;		
    else if (calling->sub[0] != NULL) calling->sub[0] = cov;	
    else if (calling->sub[1] != NULL) calling->sub[1] = cov;	
    else ERR("structure mismatch");			

    // printf("\n");
    //PMI(pts);
    // printf("%10g %10g\n", PARAM(pts, 1)[0] *PARAM(pts, 1)[1] * 4 , pts->mpp.mMplus[0]);
    
    //assert(FABS(PARAM(pts, 1)[0] *PARAM(pts, 1)[1] * 4 - pts->mpp.mMplus[0]) < 1e-14 *  pts->mpp.mMplus[0]);

    //   PMI(save_pts);
    //APMI(pts);

    DORANDOM(pts, cov->q); // cov->q ist hier nur dummy    

    // if (!newshape) printf("failed\n");
    //printf("old %10g %10g %10g\n", (double) pgs->old_zhou, (double) pgs->totalmass, PARAM0(shape->sub[0], POWSCALE));

    // assert(pgs->old_zhou == pgs->totalmass);
    // if (pgs->old_zhou > 1e6 && PARAM0(shape->sub[0], POWSCALE)>1) APMI(cov);  //    assert(shape->sub[0]->nr == POWER_DOLLAR && pts->nr == LOC);
    // assert(PARAM0(shape->sub[0], POWSCALE) == PARAM0(pts, LOC_SCALE));


    //
    if (false)
      printf("//zhou %4.1f %4.3f %ld mean=%4.3f, shape=%.50s %10g \n", (double) pgs->sum_zhou_c, pgs->totalmass, pgs->n_zhou_c, 
	     (double) pgs->sum_zhou_c / (double) pgs->n_zhou_c,
	     NICK(shape->sub[0]), PARAM0(shape->sub[0], POWSCALE)
	     ); // assert(pgs->n_zhou_c <= 500);
  }

  int idx = DrawCathegory(pgs->size, single, total, loc->grid, &elmts);
 
  if (flathull) for (d=0; d<dim; d++) x[d] = 0.0; 
  else for (d=0; d<dim; d++) x[d] = halfstepvector[d];

  assert(halfstepvector[0] != 0.0);
  assert(*x != 0.0);

  if (dim <= 2) {
    switch (idx) {
    case ZHOU_VOXEL : 
      break;
    case ZHOU_CORNER :
      for (d=0; d<dim; d++) x[d] = RF_INF; 
      break;
    case ZHOU_2D_SIDE1 : x[0] = RF_INF;
      break;
    case ZHOU_2D_SIDE2 : x[1] = RF_INF;
      break;
    default : BUG;
    }
  } else if (dim == 3) {
    switch (idx) {
    case ZHOU_VOXEL : 
      break;
    case ZHOU_CORNER :
      for (d=0; d<dim; d++) x[d] = RF_INF; 
      break;
    case ZHOU_3D_SIDE1 : x[0] = RF_INF;
      break;
    case ZHOU_3D_SIDE2 : x[1] = RF_INF;
      break;
    case ZHOU_3D_SIDE3 : x[2] = RF_INF;
      break;
    case ZHOU_3D_EDGE1 : x[1] = x[2] = RF_INF;
      break;
    case ZHOU_3D_EDGE2 : x[0] = x[2] = RF_INF;
      break;
    case ZHOU_3D_EDGE3 : x[0] = x[1] = RF_INF;
      break;
    default : BUG;
    }    
  } else BUG;
 
  //APMI(cov);
  //printf("flathull = %d %ld %.50s i=%d %d \n", flathull, x, NICK(pts), i, ZHOU_VOXEL);
 
  if (flathull) { 
    // to do: programmierfehler!
    // frage: lohnt sich 'flathull=TRUE' ueberhaupt noch hinsichtlich Aufwand
    //        und Rechenzeitersparnis?
    // Bsp das nicht funktioniert fuer flathull=FLATHULL_UNDETERMINED
    //    x _ seq(1, 10, 0.01)
    //      z _ RFsimulate(RPsmith(RMgauss(), xi=0), x, x, grid=T, print=20)
    //      Print(EXP(-EXP(range(-z[[1]]))))
    //      plot(z); X11(); hist(z[[1]], freq=FALSE, 50)
    //			curve(EXP(-x-2) * EXP(-EXP(-x-2)), -5, 2, add=TRUE)
    if (idx != ZHOU_VOXEL) VTLG_R(x, pts, v); 
    //  printf("v=%10g x=%10g %d %d %d \n", v[0],  x[0], i, ZHOU_VOXEL,  ZHOU_CORNER);
    for (d = 0; d<dim; d++)
      if (x[d] == 0.0) 
	v[d] = loc->xgr[d][XSTEP] * UNIFORM_RANDOM - halfstepvector[d];
  } else { 
    //    printf("here %10g %10g %10g\n", x[0], x[1], x[2]);
    // double xx[3]; /* print */
    //  xx[0] = xx[1] = 0.5;   xx[2] = RF_INF;  bad
    // xx[0] = 0.5;  xx[1] = xx[2] = RF_INF; badd
    // xx[0] = xx[1] = xx[2] = RF_INF; ok
    //xx[0] = 0.5;  xx[1] = xx[2] = RF_INF;
    
    VTLG_R2SIDED(NULL, x, pts, v); 
  }

  // v[0] = (2 * UNIFORM_RANDOM - 1) * PARAM0(shape->sub[0], POWSCALE) ;
  //  v[1] = (2 * UNIFORM_RANDOM - 1) * PARAM0(shape->sub[0], POWSCALE);
  // v[2] = (2 * UNIFORM_RANDOM - 1) * PARAM0(shape->sub[0], POWSCALE);

  //  if (v[0] * v[0] + v[1] * v[1] +  v[2] *  v[2] > 0.1001 * 0.1001 * 2)
  // printf("i=%d %10g %10g %10g v=%10g %10g %10g\n", i, x[0], x[1], x[2], v[0], v[1], v[2]);

  //printf("i=%d %10g %10g v=%10g %10g\n", i, x[0], x[1], v[0], v[1]);

  for (d=0; d<dim; d++) {
    cov->q[d] = loc->xgr[d][XSTART] + v[d];
    // 
    // printf("q=%10g %10g %10g %d\n", cov->q[d], loc->xgr[d][XSTART], v[d], flathull);
//APMI(pts->calling);
    if (R_FINITE(x[d])) {      
      int 
	len = (int) loc->xgr[d][XLENGTH] - 1,
	nr = elmts % len;
      elmts /= len;
      cov->q[d] += loc->xgr[d][XSTEP] * (nr + (v[d] > 0.0));	
      // simuliert von VTLF_R2SIDED ist die abgeschnittene Dichte symmetrisch
      // um 0; in Ursprung muss die Dichtekurve aufgeschnitten werden und 
      // versetzt wieder zusammengesetzt werden, so dass in der Mitte i.a.
      // der kleinste Werte ist und aussen die groessten, ehemaligen Werte
      // in der Naehe des Ursprungs. Dies bewirkt gerade (v>0) und "+ v"
    } else {
      if (v[d] > 0.0) 
	cov->q[d] += loc->xgr[d][XSTEP] * (loc->xgr[d][XLENGTH] - 1.0);	
      // dichte function der pts muss "zerschnitten" werden; der negative
      // Teil wird vor den Beginn des Gitters gesetzt, der positive Teil
      // ans Ende des Gitters.  Somit immer "+v", nie "-v";
      // loc->xgr[0][XSTEP] * (loc->xgr[0][XLENGTH] - 1.0) schiebt den 
      // Punkt ans Ende des Gitters.       
      x[d] = v[d];
    }
    //   printf("D=%d %10g %10g %d\n", d, cov->q[d], v[d], elmts);

  }
 
  if (flathull) for (d = 0; d<dim; d++) if (x[d] == 0.0)  v[d] = 0.0;
 
  pgs->log_density = LOG(pts->mpp.unnormedmass);
  if (!R_FINITE(pgs->log_density)) {  BUG; }
}



void do_Zhou(model *cov, gen_storage *S) {// = pointshapetype
  // muss zu allerst stehen, da cov sich aendern kann!
  assert(hasSmithFrame(cov));
  do_pgs_maxstable(cov, S); 

  // do_gs_maxstable might break links of the cov structure:
  model *calling=cov->calling;		
  assert(calling != NULL);				
  if (calling->key != NULL) cov = calling->key;		
  else if (calling->sub[0] != NULL) cov = calling->sub[0];	
  else if (calling->sub[1] != NULL) cov = calling->sub[1];	
  else ERR("structure mismatch");			
  
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  pgs_storage *pgs = cov->Spgs;
  int d,
    dim = XDIM(PREVSYSOF(shape), 0);
  double 
    eps,
    *x = pgs->x,
    *y = pgs->y;

  eps = pgs->currentthreshold;
  assert(R_finite(eps));

  if (!R_FINITE(pgs->log_density)) {
    BUG;
  }
  if (cov->loggiven)  eps += pgs->log_density;
  else eps *= EXP(pgs->log_density);
  
  if (cov->loggiven) { NONSTATLOGINVERSE(&eps, shape, x, y); }
  else NONSTATINVERSE(&eps, shape, x, y);

  //  printf("%d %f %f--%f %s\n",cov->loggiven, eps, x[0], y[0], NAME(shape));

  if (ISNAN(x[0]) || x[0] > y[0]) {
    //double eps_neu = eps / cov->mpp.maxheights[0]; // warum ?? 29.12.2013
    double eps_neu = eps; //  29.12.2013
    if (cov->loggiven) { BUG; }
    else NONSTATINVERSE_D(&eps_neu, pts, x, y);    
    if (ISNAN(x[0])  || x[0] > y[0]) BUG;
  }

  for (d=0; d<dim; d++) { 
    pgs->supportmin[d] = cov->q[d] - y[d] ; // 4 * for debugging...
    pgs->supportmax[d] = cov->q[d] - x[d] ;
    if (ISNAN(pgs->supportmin[d]) || ISNAN(pgs->supportmax[d]) ||
	pgs->supportmin[d] > pgs->supportmax[d]) {
      BUG;
    }
  }

  // printf("q=%f %f sup=(%f %f)-(%f %f)\n", cov->q[0], cov->q[1], pgs->supportmin[0],pgs->supportmin[1] , pgs->supportmax[0],pgs->supportmax[1]);
  
  cov->fieldreturn = shape->fieldreturn;
}


void range_Zhou(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[ZHOU_RATIO] = 0;
  range->max[ZHOU_RATIO] = 1;
  range->pmin[ZHOU_RATIO] = 0;
  range->pmax[ZHOU_RATIO] = 1;
  range->openmin[ZHOU_RATIO] = false;
  range->openmax[ZHOU_RATIO] = false; 

  booleanRange(ZHOU_FLATHULL);
  booleanRange(ZHOU_INFTY_SMALL);
  booleanRange(ZHOU_NORMED);
  booleanRange(ZHOU_ISOTROPIC);

}





#define STAT_SHAPE_FCT 0
void stationary_shape(double *x, model *cov, double *v) { // = pointshapetype
  model *shape = cov->sub[STAT_SHAPE_FCT];
  FCTN(x, shape, v);
  //  if (cov->q[0] < 1) {
  //   printf("x=%10g q=%10g v=%10g\n", x[0], cov->q[0], v[0]);
  //   assert(false);
  // } 
  //  APMI(cov);
}

void logstationary_shape(double *x, model *cov, double *v, double *Sign) { 
  model *shape = cov->sub[STAT_SHAPE_FCT];
  LOGCOV(x, shape, v, Sign);
}

int check_stationary_shape(model *cov) {
  model *shape = cov->sub[STAT_SHAPE_FCT];
  int err, 
    dim = OWNLOGDIM(0); 

   ASSERT_CARTESIAN;
   ASSERT_UNREDUCED;
  
  if ((err = CHECK(shape, dim, dim, ProcessType, XONLY, CARTESIAN_COORD, 
		   SCALAR, cov->frame
		   )) != NOERROR)  RETURN_ERR(err);
  setbackward(cov, shape);

  RETURN_NOERROR;
}

int struct_stationary_shape(model *cov, model **newmodel){
  model *shape = cov->sub[STAT_SHAPE_FCT];
  //  location_type *loc = Loc(cov);

  ASSERT_NEWMODEL_NULL;

  if (!hasPoissonFrame(shape) && !hasSchlatherFrame(shape)) ILLEGAL_FRAME;

  //printf("here\n");

  RETURN_NOERROR;
}


int init_stationary_shape(model *cov, gen_storage *S) {  
  ASSERT_ONESYSTEM;
  model *shape = cov->sub[STAT_SHAPE_FCT];
  int d, i,
    err = NOERROR,
    dim = XDIM(PREVSYSOF(shape), 0);

  assert(dim == Loc(cov)->timespacedim);
  if ((err = alloc_pgs(cov)) != NOERROR) RETURN_ERR(err);
  pgs_storage *pgs = cov->Spgs;

  // selbst wenn zufaelliger Shape: 1x laufen lassen, ob 
  // Fehler auftauchen. Unter "Do" lassen sie sich nicht mehr
  // schoen abfangen.
   
  if ((err = INIT(shape, 1, S)) != NOERROR) {
    RETURN_ERR(err); //gatter?
  }

  if (shape->mpp.moments >= 1) {
    for (i=0; i<=cov->mpp.moments; i++) {
      cov->mpp.mM[i] = shape->mpp.mM[i];
      cov->mpp.mMplus[i] = shape->mpp.mMplus[i];
    }
    pgs->zhou_c = 1.0 / cov->mpp.mMplus[1]; // passt fuer binary, und auch fuer 
    if (!R_FINITE(pgs->zhou_c))
      SERR1("max height of '%.50s' not finite", NICK(shape));
    pgs->estimated_zhou_c = false;
  } else {
    assert(hasSchlatherFrame(cov)); // meaning any other process with
    assert(shape->mpp.maxheights[0] == 1);
    // unknown properties.
    pgs->zhou_c = pgs->sq_zhou_c = pgs->sum_zhou_c = 0.0;
    pgs->estimated_zhou_c = true;
    pgs->n_zhou_c = 0;
    pgs->logmean = false;
  }

  if (!isProcess(shape)) // war: if (cov->randomkappa) 
    SERR("shapes must be a (stationary) process in stationary modelling -- pls contact author");
  
  pgs->log_density = 0; 

  for (d=0; d<dim; d++) {
    pgs->supportmin[d] = RF_NEGINF; // 4 * for debugging...
    pgs->supportmax[d] = RF_INF;
  }

  cov->mpp.maxheights[0] = shape->mpp.maxheights[0];
  ReturnOtherField(cov, shape);
  if (cov->fieldreturn == falsch) BUG;
  
  RETURN_NOERROR;
}


void do_stationary_shape(model *cov, gen_storage *S) {
  pgs_storage *pgs = cov->Spgs;
  model *shape = cov->sub[STAT_SHAPE_FCT]; 
  DO(shape, S);
//printf("%ld %ld\n", cov->rf, shape->rf);
//printf("%10g\n", shape->rf[0]);
// APMI0(cov);

  if (pgs->estimated_zhou_c) {
    BUG;
    //    pgs->zhou_c += max;
    //    pgs->n_zhou_c++;
  }
  
  cov->mpp.maxheights[0] = shape->mpp.maxheights[0];
  assert(shape->fieldreturn == wahr);
}


void standard_shape(double *x, model *cov, double *v) { // = pointshapetype,
  // poisson & smith
  model *shape = cov->sub[PGS_FCT];
  NONSTATCOV(x, cov->q, shape, v);
}

void logstandard_shape(double *x, model *cov, double *v, double *Sign) { 
  model *shape = cov->sub[PGS_FCT];
  LOGNONSTATCOV(x, cov->q, shape, v, Sign);
}

int check_standard_shape(model *cov) {
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  int err, 
    dim =  ANYDIM;
  Types frame;
  
  ASSERT_CARTESIAN;
  ASSERT_UNREDUCED;

  if (cov->q == NULL) QALLOC(dim);
  
  if (hasPoissonFrame(cov)) {
    frame = PoissonType;
  } else if (hasSmithFrame(cov)) {
    frame = cov->frame;
  } else ILLEGAL_FRAME;


  if ((err = CHECK(shape, dim, dim, ShapeType, XONLY, CARTESIAN_COORD, 
		   SCALAR, frame)) != NOERROR)  RETURN_ERR(err);
  setbackward(cov, shape);
  
  if (shape->randomkappa)
    SERR1("random shapes for '%.50s' not allowed yet", NICK(cov));

  if (pts != NULL) {
    if ((err = CHECK_R(pts, dim)) != NOERROR) RETURN_ERR(err);
  }

  RETURN_NOERROR;
}



int struct_standard_shape(model *cov, model **newmodel){
  model *shape = cov->sub[PGS_FCT];

  // if (newmodel != NULL) crash();
  
  ASSERT_NEWMODEL_NULL;

  if (!hasPoissonFrame(shape) && !hasSmithFrame(shape)) ILLEGAL_FRAME;
  
  RETURN_NOERROR;
}
 



int init_standard_shape(model *cov, gen_storage *S) {  
  ASSERT_ONESYSTEM;
  model *shape = cov->sub[PGS_FCT];
  //  location_type *loc = Loc(cov);

  assert(cov->sub[PGS_LOC] != NULL);
  location_type *loc = Loc(cov);
  int d,
    dim = XDIM(PREVSYSOF(shape), 0),
     err = NOERROR;
  pgs_storage *pgs = cov->Spgs;

  assert(XDIM(PREVSYSOF(shape), 0) == Loc(cov)->timespacedim);

  if (pgs == NULL) {
    if ((err = alloc_pgs(cov)) != NOERROR) RETURN_ERR(err);
    pgs = cov->Spgs;
    if ((pgs->localmin = (double*) CALLOC(dim, sizeof(double))) == NULL ||
	(pgs->localmax = (double*) CALLOC(dim, sizeof(double))) == NULL ||
	(pgs->minmean = (double*) CALLOC(dim, sizeof(double))) == NULL ||
	(pgs->maxmean = (double*) CALLOC(dim, sizeof(double))) == NULL)
      RETURN_ERR(ERRORMEMORYALLOCATION);
  }

 
  // selbst wenn zufaelliger Shape: 1x laufen lassen, ob 
  // Fehler auftauchen. Unter "Do" lassen sie sich nicht mehr
  // schoen abfangen.
  // FRAME_ASSE RT(Poisson || cov->frame = = M axStableType);

  //  PMI(cov->calling);
  //if (cov->mpp.moments != hasSmithFrame(cov)) crash();
  
  assert(cov->mpp.moments == hasSmithFrame(cov));
  
  if ((err = INIT(shape, cov->mpp.moments, S)) != NOERROR) RETURN_ERR(err); //gatter?
 
  //APMI(cov);

  model *u = cov->sub[PGS_LOC]; 
  assert(MODELNR(u) == UNIF);
  double 
    *min = PARAM(u, UNIF_MIN),
    *max = PARAM(u, UNIF_MAX),
    *x = pgs->minmean, // !! wird gespeichert
    *y = pgs->maxmean; // !!
  
  // try out whether it works:
  NONSTATINVERSE(ZERO(shape), shape, x, y);

  //  printf("x=%10g %10g %10g %10g\n", x[0], x[1], y[0], y[1]); BUG;
  
  if (ISNAN(x[0]) || x[0] > y[0])
    SERR1("inverse of '%.50s' unknown", NICK(shape));
  GetDiameter(loc, pgs->localmin, pgs->localmax, 
	      pgs->supportcentre); // last is only a dummy
  pgs->totalmass = 1.0;
  for (d=0; d<dim; d++) {
    min[d] = pgs->localmin[d] - y[d];
    max[d] = pgs->localmax[d] - x[d];
    if (!(R_FINITE(min[d]) && R_FINITE(max[d])))
      SERR1("simulation window does not have compact support. Should '%.50s' be used?", DefList[TRUNCSUPPORT].nick);
    pgs->totalmass *= max[d] - min[d];
  }
  
  if (hasPoissonFrame(cov)) {
    pgs->log_density = 0.0;
  } else {    
    pgs->log_density = 0.0;
    pgs->zhou_c = pgs->totalmass / shape->mpp.mMplus[1];
    cov->mpp.maxheights[0] = shape->mpp.maxheights[0];
    pgs->estimated_zhou_c = cov->randomkappa;
    if (pgs->estimated_zhou_c) SERR("random shapes in standard approach not coded yet -- please contact author");
  }

  ReturnOtherField(cov, shape);

  RETURN_NOERROR;
}


void do_standard_shape(model *cov, gen_storage *S) {
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC]; 
  assert(cov->sub[PGS_LOC] != NULL);
  pgs_storage *pgs = cov->Spgs;
  double *x = pgs->x,
    *y = pgs->xstart;
  int d,
    dim = XDIM(PREVSYSOF(shape), 0);

  DO(shape, S);
  DORANDOM(pts, cov->q);

  ///  assert(shape->fieldreturn != falsch);  3.4.2018. Warum?


  //PMI(shape, "standard");
  NONSTATINVERSE(ZERO(shape), shape, x, y);
  if (ISNAN(x[0])|| x[0] > y[0]) BUG;

  for (d=0; d<dim; d++) {
    pgs->supportmin[d] = cov->q[d] - y[d]; // 4 * for debugging...
    pgs->supportmax[d] = cov->q[d] - x[d];

    assert(pgs->supportmin[d] != NA_INTEGER && 
	   pgs->supportmax[d] != NA_INTEGER);
    
    assert(pgs->supportmin[d] <= pgs->supportmax[d]);
  }  

  pgs->log_density = 0.0;

}



 






void mcmc_pgs(double VARIABLE_IS_NOT_USED *x, model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *v) { 
}

void logmcmc_pgs(double *x, model *cov, double *v, double *Sign) { 
  model *shape = cov->sub[PGS_FCT];
  LOGNONSTATCOV(x, cov->q, shape, v, Sign);

  //  printf("x = %10g %10g log=%10g\n", *x, cov->q[0], *v, *Sign);
}
int check_mcmc_pgs(model *cov) { // = pointshapetype
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  location_type *loc = Loc(cov);
  int err, 
    dim = OWNLOGDIM(0);
  Types frame;

  ASSERT_CARTESIAN;
  ASSERT_UNREDUCED;
 if (loc->Time) SERR("Time component not allowed yet"); // todo
   
  kdefault(cov, ZHOU_RATIO, GLOBAL.extreme.density_ratio); 
  kdefault(cov, ZHOU_FLATHULL, (int) GLOBAL.extreme.flathull);
  kdefault(cov, ZHOU_INFTY_SMALL, P0INT(ZHOU_FLATHULL) != (int) False);
  kdefault(cov, ZHOU_NORMED, true);
  kdefault(cov, ZHOU_ISOTROPIC, true);

  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);

  if (cov->q == NULL) QALLOC(dim);
  
  if (hasGaussMethodFrame(cov))
    frame = isGaussMethod(shape) || equalsBernoulliProcess(shape)
      ? GaussMethodType : cov->frame ; 
  else if (hasPoissonFrame(cov)) frame = PoissonType;
  //  else if (hasMaxStableFrame(cov)) frame = cov->frame;
  else ILLEGAL_FRAME;
  

  //    printf("here %d\n");
  
  // isotropy is the property -- user should not define unisotropic functions!
  // so !isotropic only for special, internal definitions, e.g. 
  // strokorbPoly
  // int err2 = NOERROR;
  //  if (P0INT(ZHOU_ISOTROPIC)) 
  //    err2 = CHECK(shape, dim, dim, ShapeType, XONLY, ISOTROPIC,
  //		 SCALAR, frame);

  //  printf("here2\n");
  // PMI(shape);
 // but it must be treated in any case as if it was non-isotropic
  if ((err = CHECK(shape, dim, dim, ShapeType, XONLY, CARTESIAN_COORD,
		   SCALAR, frame)) != NOERROR) {
    if (P0INT(ZHOU_ISOTROPIC)) BUG;
    XERR(err);
    RETURN_ERR(err);
  }
  // if (err2 != NOERROR) SERR("isotropic shape function expected.");

  setbackward(cov, shape);

  //if (shape->randomkappa) RETURN_ERR(ERRORNOTPROGRAMMED); // Dichte muss nicht normiert sein und stattdessen durch das mittlere Volumen dividiert werden.

  if (pts != NULL) {
    if ((err = CHECK_R(pts, dim)) != NOERROR) RETURN_ERR(err);
  }

  //  printf("done\n");
 
  RETURN_NOERROR;
}



int struct_mcmc_pgs(model VARIABLE_IS_NOT_USED *cov, model VARIABLE_IS_NOT_USED **newmodel){
  //  model *shape = cov->sub[PGS_FCT];
  //  location_type *loc = Loc(cov);
  //  int err = NOERROR;
  ASSERT_NEWMODEL_NULL;
 RETURN_NOERROR;
}


int init_mcmc_pgs(model *cov, gen_storage VARIABLE_IS_NOT_USED *S) {  
  ASSERT_ONESYSTEM;
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  int i,
    err = NOERROR;

  for (i=0; i<=cov->mpp.moments; i++) {
    cov->mpp.mM[i] = shape->mpp.mM[i] * pts->mpp.mMplus[0];
    cov->mpp.mMplus[i] = shape->mpp.mMplus[i] * pts->mpp.mMplus[0];
  }
  cov->mpp.maxheights[0] = RF_NAN; // 2.2.19 -- zur sicherheit

  ReturnOtherField(cov, shape);

  RETURN_ERR(err);
}


void do_mcmc_pgs(model  VARIABLE_IS_NOT_USED *cov, gen_storage  VARIABLE_IS_NOT_USED *S) {
}

void range_mcmc_pgs(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  
  range->min[ZHOU_RATIO] = 0;
  range->max[ZHOU_RATIO] = 1;
  range->pmin[ZHOU_RATIO] = 0;
  range->pmax[ZHOU_RATIO] = 1;
  range->openmin[ZHOU_RATIO] = false;
  range->openmax[ZHOU_RATIO] = false; 

  booleanRange(ZHOU_FLATHULL);
  booleanRange(ZHOU_INFTY_SMALL);
  booleanRange(ZHOU_NORMED);
  booleanRange(ZHOU_ISOTROPIC);

}







/// ballani




int calculate_mass_gauss(model *cov) { // = pointshapetype
  /* 
     erste Idee bei Gitter --- diese ist aber viel zu kompliziert,
     insbesondere beim Auswerten der Dichtefunktion bei der Simulation:
     jeweils von der Seite wird mit konstanter, groesserer 
     Schrittweite (newstep) ins zentrum gegangen. Das Verfahren waere
     schoen wenn nicht in der mitte etwas Chaos gegeben wuerde, sofern
     die Schrittweite nicht aufgeht (!addingup[d]):
     (a) die Massen um die Punkte am naechsten zum Zentrum 
         ueberlappen sich sehr (oldmass > reference)
         dann weiss nicht so richtig was ich tun soll; hier wird einfach
         alles so gelassen wie es ist
     (b) die Massen ueberlappen sich nicht genug. Dann werden 
	 zusaetzliche Punkte aufgenommen (und zwar fuer alle 2^d bzw.
         2^d -1 Moeglichkeiten einzeln, ob die Punkte dem Zentrum
         am naechsten in einer Richtung sind (fixed[d]=true) oder nicht.

     Zweite Idee (ab V.2.1.31): die Dichten werden auf einem "beliebigen"
     Gitter (pgs->xgr) definiert. Dies vereinfacht die Berechnung der 
     Dichten. Ob dieses Gitter auf den zu simulierenden Punkten liegen
     muss oder nicht sei dahingestellt (auch wenn Felix' Resultate
     zeigen, dass dies aus theoretischen Gruenden so sein muss). Was
     letztendlich die kleinste reale Rechenzeit verursacht ist nicht
     klar. 
     Auch sind die Gewichte voellig unklar. Insbesondere wo/wieviele
     Nullen auftreten und die Schwere der Randeffekte bei ueberall 
     konstantem Gewicht.

     Um hier den ersten Algorithmus einfach zu halten werden konstante
     Gewichte genommen und (dafuer das Gebiet eher etwas groesser 
     gemacht)
     
     TODO!

  */  

  pgs_storage *pgs = cov->Spgs;
  location_type *loc = Loc(cov); 
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  double zx, zy, 
    **xgr = pgs->xgr, // out
    *v = pgs->v,             // dummy
    *x = pgs->x,             // dummy  
    *y = pgs->y;             // dummy  

  int d, 
    dim = OWNLOGDIM(0);

  if (loc->grid) {
    // int size = pgs->size;
    COV(ZERO(shape), shape, v);
    v[0] *= intpow(0.5, dim);
    NONSTATINVERSE_D(v, shape, x, y); 
    if (ISNAN(x[0]) || x[0] > y[0]) 
      SERR1("inverse function of '%.50s' unknown", NICK(shape));

  
    // is the underlying point distribution uniform? Then
    // no safety division by SQRT(dim) seems to be necessary // todo
    // unif might be modified by $. So a check on the equality
    // of three values is done
    VTLG_D(ZERO(pts), pts, v);
    VTLG_D(x, pts, &zx);
    VTLG_D(y, pts, &zy);
 
    //printf("> x=%10g %10g;   %10g %10g %10g\n", x[0], x[1], zx, zy, v[0]);

    for (d=0; d<dim; d++) {
      y[d] -= x[d]; // !! this line must be before next line
      // and after VTLG_D(y, pts, &zy); !!
      assert(y[d] > 0);
    }
    if (zx != v[0] || zy != v[0] ||  true) {
      for (d=0; d<dim; d++) {
	y[d] /= SQRT((double) dim); // falls nicht runif
      }
    }
    
    //    PMI(shape);    printf("%10g %10g\n", v[0], x[0]);
    //     printf("! x=%10e %d\n", x[0], (int) -0.5); assert(false);

    // um Pythagoras fuer dim > 1
    pgs->totalmass = 1;
    for (d=0; d<dim; d++) { // Berechnung der ausgeduennten Reihe
      if (loc->xgr[d][XLENGTH] > 1) {
	double range = loc->xgr[d][XSTEP] * (loc->xgr[d][XLENGTH] - 1.0);
	xgr[d][XLENGTH] = CEIL(range / y[d] + 1.0);
	if (xgr[d][XLENGTH] >= loc->xgr[d][XLENGTH]) {
	  BUG;
	  MEMCOPYX(xgr[d], loc->xgr[d], 3 * sizeof(double));
	} else {
	  xgr[d][XSTART] =
	    loc->xgr[d][XSTART] - 0.5 *((xgr[d][XLENGTH] - 1.0) * y[d] - range);
	  xgr[d][XSTEP] = y[d];
	}
	pgs->totalmass *= xgr[d][XLENGTH];
      
      // printf("%d %10g %10g %10g \n", d,xgr[d][XSTART], xgr[d][XSTEP], xgr[d][XLENGTH] ); 
	assert(xgr[d][XSTEP] >= 0);
      } else {
	int i; assert(XSTART == 0 && XLENGTH == 2);
	for (i=XSTART; XLENGTH; i++) xgr[d][i] = loc->xgr[d][i];
      }

    }
    //for (d=0; d<=1; d++)
    //    printf("%d %10g %10g %10g \n", d,xgr[d][XSTART], xgr[d][XSTEP], xgr[d][XLENGTH] );    assert(false);

  } else { // not grid
    pgs->totalmass = loc->totalpoints;
  }

  RETURN_NOERROR;
}



double gauss_eps = 1e-10;
void do_pgs_gauss(model *cov, gen_storage *S) {// = pointshapetype
  pgs_storage *pgs = cov->Spgs;
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  location_type *loc = Loc(cov);
  //  window_info *w = &(S->window);
  long i;
  int d, 
    *min = pgs->min, // dummy
    *max = pgs->max, // dummy
    *pos = pgs->pos, // dummy
    dim = XDIM(PREVSYSOF(shape), 0);
  double value,
    total = 0.0,
    *x = pgs->x, // dummy or (future) Static local
    *y = pgs->y, // dummy or (future) Static local
    //    *supportmin = pgs->supportmin,
    //   *supportmax = pgs->supportmax,
    // *pgstotal = pgs->total, // in
    //*single = pgs->single, // in
    **xgr = pgs->xgr, // in
    *v = pgs->v;    
  bool 
    grid = Loc(cov)->grid;


  if (cov->randomkappa) {
    DO(shape, S); 
    DORANDOM(pts, cov->q);  // cov->q nur dummy. Wird ueberschrieben
    if (hasPoissonGaussFrame(cov) || !grid) { 
      if (calculate_mass_gauss(cov) != NOERROR) 
	ERR("unexpected error in 'do_Zhou' (maxstable)");
    } else BUG;
  }

  // printf("name = %.50s\n", DefList[pts->nr].name);
  //  APMI(pts);
  
  VTLG_R(NULL, pts, v); 
  i = UNIFORM_RANDOM * pgs->totalmass;
  
  // to determine the points that contribute to the value of the 
  // random field by the shape function centered at cov->q
  //printf("cat=%d el=%d\n", i, elmts);
  if (loc->grid) {
    NONSTATINVERSE_D(&gauss_eps, pts, x, y);
    if (ISNAN(x[0]) || x[0] > y[0]) BUG;

    //    printf("extremes %10g %10g %10g\n", eps, x[0], y[0]);

    for (d=0; d<dim; d++) {
      //printf("e1 %10g %10g\n", x[d], y[d]);
      int which = i % (int) xgr[d][XLENGTH];
      i /= xgr[d][XLENGTH];
      cov->q[d] = xgr[d][XSTART] + xgr[d][XSTEP] * which + v[d];
      // to determine the points that contribute to the value of the 
      // density function at cov->q
      min[d] = CEIL((cov->q[d] - y[d] - xgr[d][XSTART]) / xgr[d][XSTEP]);
      max[d] = (cov->q[d] - x[d] - xgr[d][XSTART]) / xgr[d][XSTEP];

      // printf("%d %d delat=%10g %10g y=%10g %10g step=%10g\n", min[d], max[d], 
      //     cov->q[d] - y[d] - xgr[d][XSTART],
      //     cov->q[d] - x[d] - xgr[d][XSTART],
      //     y[d], x[d], xgr[d][XSTEP]);

      assert(xgr[d][XSTEP] > 0);

      if (min[d] < 0) min[d] = 0;
      if (max[d] >= xgr[d][XLENGTH]) max[d] = xgr[d][XLENGTH] -1;
      if (min[d] > max[d]) {
	// no point contributed to cov->q -- this might be possible if 
	// the simulating grid goes beyond the actual grid and 	

	//printf("redo d=%d x=%10g y=%10g q=%10g i=%d d=%d range=%d %d\n",
	// d, x[d], y[d], cov->q[d], i, d, min[d], max[d]);
	do_pgs_gauss(cov, S);
	pgs->log_density = RF_INF;
	return;
      } // else printf("ok q=%10g i=%d\n", cov->q[d], i, d, min[d], max[d]);
      pos[d] = min[d];
      y[d] = x[d] = cov->q[d] - (xgr[d][XSTART] + pos[d] * xgr[d][XSTEP]);
    }
   
    while (true) {
      VTLG_D(y, pts, &value); 
      total += value;
      
     //printf("y=%10g v=%10g q=%10g val=%10g tot=%10g min=%d max=%d pos=%d xgr.len=%d\n", 
      //     y[0], v[0], cov->q[0], value, total, min[0], max[0], pos[0],
      //(int) loc->xgr[0][XLENGTH]);
      
      d = 0;
      while(d < dim && pos[d] == max[d]) {
	pos[d] = min[d];
	y[d] = x[d];
	d++;
      }
      if (d >= dim) break;
      pos[d] ++;
      y[d] -= xgr[d][XSTEP]; // !! nicht '+' !!
    }
  } else { // not grid
    if (loc->timespacedim != dim) BUG;
    double *xx = loc->x + dim * i;
    for (d=0; d<dim; d++) cov->q[d] = v[d] + xx[d];
    // todo : mit Unterteilung des Feldes viel schneller
    xx = loc->x;
    long endfor = loc->totalpoints;
    for (i=0; i<endfor; i++, xx+=dim) {
      for (d=0; d<dim; d++) y[d] = cov->q[d] - xx[d];
      VTLG_D(y, pts, &value);
      total += value;
    }
  }
  pgs->log_density = LOG(total / pgs->totalmass);

  //printf("total=%10g %10g grid=%d dim=%d\n", total, pgs->totalmass,loc->grid,dim);

  assert(R_FINITE(pgs->log_density));
  // assert(false);

} // do_pgs_gauss


void Ballani(double *x, model *cov, double *v){
  model *shape = cov->sub[PGS_FCT];
  //    *pts = cov->sub[PGS_LOC];
  NONSTATCOV(x, cov->q, shape, v);

  //  TALLOC_X1(delta, OWNXDIM(0));
  //  closest(cov->q, cov, delta);
  //  double w;
  //  VTLG_D(delta, pts, &w);
  //  *v /= w;
  //  assert(*v < 1e-12 * pts->mpp.unnormedmass);
  // if (*v > 1e-12 * pts->mpp.unnormedmass) printf("violating v=%10e\n", *v);
  //   END_TALLOC_X1;
}

void logBallani(double *x, model *cov, double *v, double *Sign){
  model *shape = cov->sub[PGS_FCT];
    //    *pts = cov->sub[PGS_LOC];
  LOGNONSTATCOV(x, cov->q, shape, v, Sign);
}

int check_Ballani(model *cov){
  RETURN_ERR(ERRORNOTPROGRAMMEDYET);
  
  ASSERT_ONESYSTEM;
  ASSERT_UNREDUCED;
 model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  location_type *loc = Loc(cov);
  int err,
    dim = OWNLOGDIM(0);
  
  ASSERT_CARTESIAN;
  if (loc->Time) SERR("Time component not allowed yet"); // todo
   
  //  kdefault(cov, ZHOU_RATIO, GLOBAL.extreme.density_ratio); 
  //  kdefault(cov, ZHOU_FLATHULL, (int) GLOBAL.extreme.flathull);
  //  kdefault(cov, ZHOU_INFTY_SMALL, P0INT(ZHOU_FLATHULL) != (int) False);
  //  kdefault(cov, ZHOU_NORMED, true);
  //  kdefault(cov, ZHOU_ISOTROPIC, true);

  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err); 
  if (cov->q == NULL) QALLOC(dim);
  
  assert(!isProcess(cov));
  assert(hasPoissonGaussFrame(cov));

  // but it must be treated in any case as if it was non-isotropic
  if ((err = CHECK(shape, dim, dim, ShapeType, XONLY, CARTESIAN_COORD,
		   SCALAR,
		   isNormedProcess(shape) ? NormedProcessType : cov->frame)
       ) != NOERROR) {
    if (P0INT(ZHOU_ISOTROPIC)) BUG;
    XERR(err);
    RETURN_ERR(err);
  }

  setbackward(cov, shape);

  //if (shape->randomkappa) RETURN_ERR(ERRORNOTPROGRAMMED); // Dichte muss nicht normiert sein und stattdessen durch das mittlere Volumen dividiert werden.

  if (pts != NULL) {
    if ((err = CHECK_R(pts, dim)) != NOERROR) RETURN_ERR(err);
  }
  
  EXTRA_STORAGE;
  RETURN_NOERROR;
}
int struct_Ballani(model *cov, model **newmodel){
  BUG; // should not be called at all! 2.2.19

  
  RETURN_ERR(ERRORNOTPROGRAMMEDYET);
  model *shape = cov->sub[PGS_FCT];
  //  location_type *loc = Loc(cov);
  int err = NOERROR;

  ASSERT_NEWMODEL_NULL;
  if (cov->Spgs != NULL)  pgs_DELETE(&(cov->Spgs), cov);

  if (!hasPoissonGaussFrame(shape)) ILLEGAL_FRAME;
  if (cov->sub[PGS_LOC] == NULL) {
    if ((err = STRUCT(shape, cov->sub + PGS_LOC)) != NOERROR)
      RETURN_ERR(err);
  

    if (cov->sub[PGS_LOC] == NULL) 
      SERR1("no intensity found for '%.50s'", NICK(shape));

    assert(cov->sub[PGS_LOC]->calling != NULL);
  }
  RETURN_NOERROR;


}
int init_Ballani(model *cov, gen_storage *S){
  RETURN_ERR(ERRORNOTPROGRAMMEDYET);
  ASSERT_ONESYSTEM;
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  defn *Cshape = DefList + MODELNR(shape);
 location_type *loc = Loc(cov);
  int d, i,
    err = NOERROR,
     dim = XDIM(PREVSYSOF(shape), 0);
  pgs_storage *pgs = cov->Spgs;
  bool 
    grid = Loc(cov)->grid,
    pgsnull = pgs == NULL;
  if (!grid) { // todo
    SERR("non-grid not programmed yet"); // meiste da; fehlt noch
    // welche Vektoren ich wirklich brauche
  }

  assert(hasPoissonGaussFrame(cov));
  assert(cov->sub[PGS_LOC] != NULL);
  assert(dim == XDIM(PREVSYSOF(cov->sub[PGS_LOC]), 0));
  assert(dim == Loc(cov)->timespacedim);
  
  if (Cshape->inverse == ErrInverse)//&&Cshape->nonstat_inverse==Err CovNonstat)
    SERR1("support of the model is unknown. Use '%.50s' to determine the support",
	  DefList[TRUNCSUPPORT].nick);
 
  if (pgsnull) {
    if ((err = alloc_pgs(cov)) != NOERROR) RETURN_ERR(err);
    pgs = cov->Spgs;
    
    if ((pgs->v = (double*) CALLOC(dim, sizeof(double))) == NULL ||
	(pgs->y = (double*) CALLOC(dim, sizeof(double))) == NULL
	) RETURN_ERR(ERRORMEMORYALLOCATION);
  }
  pgs->size = loc->totalpoints;
  if (grid) pgs->size = intpow(2, dim);
 
  // selbst wenn zufaelliger Shape: 1x laufen lassen, ob 
  // Fehler auftauchen. Unter "Do" lassen sie sich nicht mehr
  // schoen abfangen.


  if ((err = INIT(pts, 0, S)) != NOERROR) RETURN_ERR(err);// , 0, S statt ,1,S
  if ((err = INIT(shape, cov->mpp.moments, S)) != NOERROR)  RETURN_ERR(err);
  assert(cov->mpp.moments >= 1);
  assert(shape->mpp.moments >= 1);

  
  if (pgsnull &&
      ((pgs->xgr[0] = (double*) CALLOC(dim * 3, sizeof(double))) == NULL ||
       (pgs->pos = (int *) CALLOC(dim, sizeof(int))) ==NULL ||
       (pgs->min = (int*) CALLOC(dim, sizeof(int))) == NULL ||
       (pgs->max = (int*) CALLOC(dim, sizeof(int))) == NULL
       )) RETURN_ERR(ERRORMEMORYALLOCATION);
  for (i=3, d=1; d<dim; d++, i+=3) pgs->xgr[d] = pgs->xgr[0] + i;
  // printf("%d\n", pgs->xgr[1] - pgs->xgr[0]); assert(false);
  if ((err = calculate_mass_gauss(cov)) != NOERROR) RETURN_ERR(err);


  if (DefList[MODELNR(shape)].nonstat_inverse == ErrInverseNonstat) {
    if (MODELNR(pts) != RECTANGULAR) {
      //       PMI(pts);
      // APMI(shape);
      warning("Inverse of shape function cannot be determined. Simulation speed  might be heavily decreased.");
    }
  }
  
  for (i=0; i<=cov->mpp.moments; i++) {
    cov->mpp.mM[i] = shape->mpp.mM[i] * pts->mpp.mMplus[0];
    cov->mpp.mMplus[i] = shape->mpp.mMplus[i] * pts->mpp.mMplus[0];
  }
  
  ReturnOtherField(cov, shape); // see also do for fieldreturn
 
  RETURN_ERR(err);
}


void do_Ballani(model *cov, gen_storage *S) {// = pointshapetype
  assert(hasPoissonGaussFrame(cov));  
  do_pgs_gauss(cov, S);

  // do_gs_maxstable might break links of the cov structure:
  model *calling=cov->calling;		
  assert(calling != NULL);				
  if (calling->key != NULL) cov = calling->key;		
  else if (calling->sub[0] != NULL) cov = calling->sub[0];	
  else if (calling->sub[1] != NULL) cov = calling->sub[1];	
  else ERR("structure mismatch");			
  
  model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  pgs_storage *pgs = cov->Spgs;
  int d,
    dim = XDIM(PREVSYSOF(shape), 0);
  double 
    eps,
    *x = pgs->x,
    *y = pgs->y;

  eps = GLOBAL.mpp.about_zero * EXP(pgs->log_density);

  if (cov->loggiven) { NONSTATLOGINVERSE(&eps, shape, x, y); }
  else NONSTATINVERSE(&eps, shape, x, y);

  // warum fkt obiges nicht ??


  if (ISNAN(x[0]) || x[0] > y[0]) {
    //double eps_neu = eps / cov->mpp.maxheights[0]; // warum ?? 29.12.2013
    double eps_neu = eps; //  29.12.2013
    if (cov->loggiven) { BUG; }
    else NONSTATINVERSE_D(&eps_neu, pts, x, y);    
    if (ISNAN(x[0])  || x[0] > y[0]) BUG;
  }


  for (d=0; d<dim; d++) {
    pgs->supportmin[d] = cov->q[d] - y[d] ; // 4 * for debugging...
    pgs->supportmax[d] = cov->q[d] - x[d] ;
    if (ISNAN(pgs->supportmin[d]) || ISNAN(pgs->supportmax[d]) ||
	pgs->supportmin[d] > pgs->supportmax[d]) {
      BUG;
    }
  }
  cov->fieldreturn = shape->fieldreturn;
}


void range_Ballani(model VARIABLE_IS_NOT_USED *cov,
		   range_type VARIABLE_IS_NOT_USED *range){ BUG; }
