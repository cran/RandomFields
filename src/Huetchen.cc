 /* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 simulation of max-stable random fields

 Copyright (C) 2001 -- 2014 Martin Schlather, 

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


#include <math.h>
#include <stdio.h>
#include "RF.h"
#include "Covariance.h"

 
#define POISSON_INTENSITY 0
#define RANDOMCOIN_INTENSITY (COMMON_GAUSS + 1) 



#define PGS_VOXEL 0
#define PGS_CORNER 1
#define PGS_SIDES PGS_CORNER + 1

#define PGS_2D_SIDE1 PGS_SIDES
#define PGS_2D_SIDE2 PGS_SIDES + 1

#define PGS_3D_SIDE1 PGS_SIDES
#define PGS_3D_SIDE2 PGS_SIDES + 1
#define PGS_3D_SIDE3 PGS_SIDES + 2
#define PGS_3D_EDGE1 PGS_SIDES + 3
#define PGS_3D_EDGE2 PGS_SIDES + 4
#define PGS_3D_EDGE3 PGS_SIDES + 5


int addUnifModel(cov_model *cov, double radius, cov_model **newmodel) {
  addModel(newmodel, UNIF, cov);
  kdefault(*newmodel, UNIF_MIN, -radius);
  kdefault(*newmodel, UNIF_MAX, radius);
  return NOERROR;
}





void pts_given_shape(double *x, cov_model *cov, double *v) { 
  cov_model *shape = cov->sub[PGS_FCT];
  NONSTATCOV(x, cov->q, shape, v);

  //printf("x = %f %f %v\n", *x, cov->q[0], *v);

}

void logpts_given_shape(double *x, cov_model *cov, double *v, double *sign) { 
  cov_model *shape = cov->sub[PGS_FCT];
  LOGNONSTATCOV(x, cov->q, shape, v, sign);

  //  printf("x = %f %f log=%f\n", *x, cov->q[0], *v, *sign);
}

int check_pts_given_shape(cov_model *cov) {
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  location_type *loc = Loc(cov);
  int err, role,
    dim = cov->tsdim; 

  if (loc->Time) SERR("Time component not allowed yet"); // todo
   
  kdefault(cov, PGS_RATIO, GLOBAL.extreme.density_ratio); 
  kdefault(cov, PGS_FLAT, GLOBAL.extreme.flat);
  kdefault(cov, PGS_INFTY_SMALL, !PINT(PGS_FLAT));
  kdefault(cov, PGS_NORMED, true);
  kdefault(cov, PGS_ISOTROPIC, true);


  //PMI(cov, "checking");
 
  if ((err = checkkappas(cov)) != NOERROR) return err;

  if (cov->q == NULL) {
    // zwingend CALLOC !
    if ((cov->q  = (double*) CALLOC(sizeof(double), dim)) == NULL)
      return ERRORMEMORYALLOCATION;
    cov->qlen = dim;
  }
  
  if (cov->xdimprev != cov->xdimown || cov->xdimprev != cov->tsdim) 
    return ERRORDIM;
  
  if (cov->role == ROLE_GAUSS) {
    role = isShape(shape) ? cov->role  //Marco, 29.5.13
      : isGaussProcess(shape) ? ROLE_GAUSS
      : shape->nr == BINARYPROC ? ROLE_GAUSS
      : ROLE_UNDEFINED; 
    ASSERT_ROLE_DEFINED(shape);
  } else if (hasPoissonRole(cov)) {
    role = ROLE_POISSON;
  } else if (hasMaxStableRole(cov)) {
    role = ROLE_MAXSTABLE;
  } else ILLEGAL_ROLE;


  //    printf("here %d\n");
  
  // isotropy is the property -- user should not define unisotropic functions!
  // so !isotropic only for special, internal definitions, e.g. 
  // strokorbPoly
  // int err2 = NOERROR;
  //  if (P0INT(PGS_ISOTROPIC)) 
  //    err2 = CHECK(shape, dim, dim, ShapeType, XONLY, ISOTROPIC,
  //		 SCALAR, role);

  //  printf("here2\n");
  // PMI(shape);
 // but it must be treated in any case as if it was non-isotropic
  if ((err = CHECK(shape, dim, dim, ShapeType, XONLY, CARTESIAN_COORD,
		   SCALAR, role)) != NOERROR) {
    if (P0INT(PGS_ISOTROPIC)) BUG;
    XERR(err);
    return err;
  }
  // if (err2 != NOERROR) SERR("isotropic shape function expected.");

  setbackward(cov, shape);

  //if (!shape->deterministic) return ERRORNOTPROGRAMMED; // Dichte muss nicht normiert sein und stattdessen durch das mittlere Volumen dividiert werden.

  if (pts != NULL) {

    //printf("pts %s %s\n", CovList[pts->nr + 1].nick, CovList[pts->nr + 1].name);
    //PMI(pts);

    //      assert(dim == 2); print("x");
    // printf("dim=%d\n", dim);
    // APMI(pts);

    if ((err = CHECK_R(pts, dim)) != NOERROR) return err;
  }

  //  printf("done\n");
 
  return NOERROR;
}



int struct_pts_given_shape(cov_model *cov, cov_model **newmodel){
  cov_model *shape = cov->sub[PGS_FCT];
  //  location_type *loc = Loc(cov);
  int err = NOERROR;

  ASSERT_NEWMODEL_NULL;
  if (cov->Spgs != NULL)  PGS_DELETE(&(cov->Spgs));

  if (shape->role != ROLE_POISSON && shape->role != ROLE_MAXSTABLE)
    ILLEGAL_ROLE;
  if (cov->sub[PGS_LOC] == NULL) {
    // COV_DELETE(cov->sub + PGS_LOC);
    // assert(cov->sub[PGS_LOC] == NULL);
    if ((err = STRUCT(shape, cov->sub + PGS_LOC)) != NOERROR)
      return err;
  

    if (cov->sub[PGS_LOC] == NULL) 
      SERR1("no intensity found for '%s'", NICK(shape));

    assert(cov->sub[PGS_LOC]->calling != NULL);

    //   printf("here %s\n", NICK(cov->sub[PGS_LOC]));
    
    //if ((err = CHECK_R(cov->sub[PGS_LOC], cov->tsdim)) != NOERROR) return err;

  }
  return NOERROR;
}




int calculate_mass_gauss(cov_model *cov) {
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
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  double zx, zy, 
    **xgr = pgs->xgr, // out
    *v = pgs->v,             // dummy
    *x = pgs->x,             // dummy  
    *y = pgs->y;             // dummy  

  int d, 
    dim = cov->tsdim;

  if (loc->grid) {
    // int size = pgs->size;
    COV(ZERO, shape, v);
    v[0] *= intpow(0.5, dim);
    NONSTATINVERSE_D(v, shape, x, y); 
    if (ISNAN(x[0]) || x[0] > y[0]) 
      SERR1("inverse function of '%s' unknown", NICK(shape));

  
    // is the underlying point distribution uniform? Then
    // no safety division by sqrt(dim) seems to be necessary // todo
    // unif might be modified by $. So a check on the equality
    // of three values is done
    VTLG_D(ZERO, pts, v);
    VTLG_D(x, pts, &zx);
    VTLG_D(y, pts, &zy);
 
    //printf("> x=%f %f;   %f %f %f\n", x[0], x[1], zx, zy, v[0]);

    for (d=0; d<dim; d++) {
      y[d] -= x[d]; // !! this line must be before next line
      // and after VTLG_D(y, pts, &zy); !!
      assert(y[d] > 0);
    }
    if (zx != v[0] || zy != v[0] ||  true) {
      for (d=0; d<dim; d++) {
	y[d] /= sqrt((double) dim); // falls nicht runif
      }
    }
    
    //    PMI(shape);    printf("%f %f\n", v[0], x[0]);
    //     printf("! x=%e %d\n", x[0], (int) -0.5); assert(false);

    // um Pythagoras fuer dim > 1
    pgs->totalmass = 1;
    for (d=0; d<dim; d++) { // Berechnung der ausgeduennten Reihe
      if (loc->xgr[d][XLENGTH] > 1) {
	double range = loc->xgr[d][XSTEP] * (loc->xgr[d][XLENGTH] - 1.0);
	xgr[d][XLENGTH] = ceil(range / y[d] + 1.0);
	if (xgr[d][XLENGTH] >= loc->xgr[d][XLENGTH]) {
	  BUG;
	  memcpy(xgr[d], loc->xgr[d], 3 * sizeof(double));
	} else {
	  xgr[d][XSTART] =
	    loc->xgr[d][XSTART] - 0.5 *((xgr[d][XLENGTH] - 1.0) * y[d] - range);
	  xgr[d][XSTEP] = y[d];
	}
	pgs->totalmass *= xgr[d][XLENGTH];
      
      // printf("%d %f %f %f \n", d,xgr[d][XSTART], xgr[d][XSTEP], xgr[d][XLENGTH] ); 
	assert(xgr[d][XSTEP] >= 0);
      } else {
	int i; assert(XSTART == 0 && XLENGTH == 2);
	for (i=XSTART; XLENGTH; i++) xgr[d][i] = loc->xgr[d][i];
      }

    }
    //for (d=0; d<=1; d++)
    //    printf("%d %f %f %f \n", d,xgr[d][XSTART], xgr[d][XSTEP], xgr[d][XLENGTH] );    assert(false);

  } else { // not grid
    pgs->totalmass = loc->totalpoints;
  }


  //  PMI(cov, "ccccc"); //  assert(false);

  return NOERROR;
}

int calculate_mass_maxstable(cov_model *cov) {

  pgs_storage *pgs = cov->Spgs;
  location_type *loc = Loc(cov);
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  double voxels, value_orig,
    *single = pgs->single, // out
    *total = pgs->total,   // out
    *halfstepvector = pgs->halfstepvector, //out
    *x = pgs->x;
  //  bool flat = false;
  int d,
    dim = cov->tsdim,
    flat =P0INT(PGS_FLAT);
 
  if (shape->role == ROLE_POISSON) {
    BUG;
    // to do : felix's paper; derzeit einfach wie max-stable
    //   if (pts->nr !=    SERR("currently, only simple domain fields are allowed that are based on Poisson processes");

  } // else {

  assert(dim <= MAXSIMUDIM);
  VTLG_D(ZERO, pts, &value_orig);

  //  PMI(cov, "calculate max");
  //  printf("%s %d %d %f\n",NICK(cov), PGS_FLAT, P0INT(PGS_FLAT), value_orig );// APMI(cov); assert(false);
  
  // VOXEL 
  for (d=0; d<dim; d++) halfstepvector[d] = 0.5 * loc->xgr[d][XSTEP];     
  if (flat == FLAT_UNDETERMINED) {
    // flat == im Kerngebiet wird konstant simuliert; ausserhalb dann
    // abfallend
    if (loc->grid) {
      double v;
      VTLG_D(halfstepvector, pts, &v);
      double threshold = 
	value_orig == RF_INF//&& cov->p[PGS_RATIO][0] == 0.0 
	? RF_INF
	: value_orig * P0(PGS_RATIO);
      pgs->flat =  threshold < v && cov->deterministic; // sicherheitshalber

      //printf("thres %f %f %f\n", threshold, v, value_orig);

    } else { // not grid
      BUG;
      // long term project: to do
      pgs->flat = true;
      /* Bei letzterem und nicht Gitter wird dann wird eine maximale 
	 Kaestchengroesse bestimmt so dass gerade noch beruehrt wird
	 und der Rest ausserhalb aller Kaestchen auf einen konstanten 
	 wert gesetzt, etwas niedriger als der hoechste Grenzwert. Was 
	 "hoch" bedeutet ist hier unklar. 
      */
    }      
  } else pgs->flat = flat;

  if (pgs->flat) {
    if (P0INT(PGS_INFTY_SMALL)) {
      //PMI(cov);
      SERR2("'%s' and '%s' may not be positive at the same time",
	    KNAME(PGS_FLAT), KNAME(PGS_INFTY_SMALL));
    }
    single[PGS_VOXEL] = value_orig;
    for (d=0; d<dim; d++) single[PGS_VOXEL] *= 2.0 * halfstepvector[d];
  } else {
    //    APMI(cov);
    //       printf("two sided %f %f %f\n", halfstepvector[0], halfstepvector[1], halfstepvector[2]); //assert(false);
    VTLG_P2SIDED(NULL, halfstepvector, pts, single + PGS_VOXEL);
    //     printf("two sided done %e %f %s\n", single[PGS_VOXEL], halfstepvector[0], NICK(pts)); //assert(false);
  }
  voxels = 1.0;
  for (d=0; d<dim; d++) voxels *= loc->xgr[d][XLENGTH] - 1.0;  
  total[PGS_VOXEL] = single[PGS_VOXEL] * voxels;
  
  // CORNER
  for (d=0; d<dim; d++) x[d] = RF_INF;  
  VTLG_P2SIDED(NULL, x, pts, single + PGS_CORNER);

  //  printf("%d %f\n", !P0INT(PGS_NORMED), single[PGS_CORNER]);

  assert(!P0INT(PGS_NORMED) || single[PGS_CORNER] == 1.0); // nur wenn normiert
  total[PGS_CORNER] = single[PGS_CORNER] + total[PGS_VOXEL];

  //  printf("ier\n");
  //  PMI(cov);

  if (dim >= 2) {
    for (d=0; d<dim; d++) {
      int i, nr = 1;
      if (pgs->flat) for (i=0; i<dim; i++) x[i] = 0.0;
      else for (i=0; i<dim; i++) x[i] = halfstepvector[i];
      x[d] = RF_INF; // !! changed
      VTLG_P2SIDED(NULL, x, pts, single + PGS_SIDES + d);

      for (i=0; i<dim; i++) {
	if (R_FINITE(x[i])) {
	  if (pgs->flat) single[PGS_SIDES + d] *= loc->xgr[i][XSTEP]; 
	  nr *= loc->xgr[i][XLENGTH] - 1.0;
	}
      }
      //
      //PMI(cov);
      //  printf("  single d=%d %d x=(%f %f %f) val=%f %d\n", d, PGS_SIDES + d, 
      //     x[0], x[1], x[2], single[PGS_SIDES + d], nr); assert(false);

      total[PGS_SIDES + d] =
	total[PGS_SIDES + d - 1] + nr * single[PGS_SIDES + d];
    }

    if (dim == 3) { // Kanten
      for (d=0; d<dim; d++) {
	int i, nr = 1.0;
	for (i=0; i<dim; i++) x[i] = RF_INF;
	x[d] = pgs->flat ? 0.0 : halfstepvector[d];
	VTLG_P2SIDED(NULL, x, pts, single + PGS_SIDES + dim + d);
	if (pgs->flat) single[PGS_SIDES + dim + d] *= loc->xgr[d][XSTEP]; 
	nr *= loc->xgr[d][XLENGTH] - 1.0;
	
	//	printf("  single d=%d %d x=(%f %f %f) val=%f %d\n", d, PGS_SIDES + dim + d,  x[0], x[1], x[2], single[PGS_SIDES + dim + d], nr); //assert(false);

	total[PGS_SIDES + dim + d] = 
	  total[PGS_SIDES + dim + d - 1] + nr * single[PGS_SIDES + dim + d];
      }
    } else if (dim > 3) BUG;
  }
  //	APMI(cov);

  pgs->totalmass = total[pgs->size - 1];
  if (!R_FINITE(pgs->totalmass)) {
    // int i; for (i=0; i<pgs->size; i++) printf("i=%d total=%f\n", i, total[i]);
    SERR("the total intensity mass is not finite");
  }

  return NOERROR;
}


cov_model *prunecov(cov_model *newmodel, cov_model *cov) {
  //printf("names= %s %s\n", NICK(newmodel->calling), NICK(cov->calling));
  cov_model
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
  COV_DELETE(&newmodel);
  return next;
}

int complete_copy(cov_model **newmodel, cov_model *cov) {
  cov_model *prev = NULL,
    *start = cov;
  int err;
  while (start->calling != NULL) start = start->calling;
  if (start->typus != InterfaceType) BUG; // nur zur Sicherheit, derzeit; kann erweitert werden
  if (start == cov) BUG;
  prev = start->key != NULL ? start->key : start->sub[0];
  if (prev->typus != ProcessType) BUG;
  
  if ((err = covcpy(newmodel, prev)) != NOERROR) return err;
  (*newmodel)->calling = cov;
  int role = prev->role ;// role_of_process(prev->nr);
  if ((err = CHECK(*newmodel, prev->tsdim, prev->xdimprev, prev->typus, 
		   prev->domprev, prev->isoprev, prev->vdim2, role)) 
      != NOERROR)  return err;  
   if ((err = STRUCT(*newmodel, NULL)) != NOERROR) { return err; }
  if (!(*newmodel)->initialised) {
    if ((err =  CHECK(*newmodel, prev->tsdim, prev->xdimprev, prev->typus, 
		      prev->domprev, prev->isoprev, prev->vdim2, role)) 
	!= NOERROR) return err;  
 
    NEW_COV_STORAGE(*newmodel, stor, STORAGE, gen_storage);
   //APMI(*newmodel);
    if ((err = INIT(*newmodel, 0, cov->stor)) != NOERROR) {
      //APMI(cov); // !!! ?? hier weitermachen
      return err; 
    }
  }
 
  ///APMI(cov);
  (*newmodel)->calling = start;
  *newmodel = prunecov(*newmodel, cov);
  (*newmodel)->calling = NULL;
  return NOERROR;
}

int init_pts_given_shape(cov_model *cov, gen_storage *S) {  
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  cov_fct *Cshape = CovList + shape->nr;
  location_type *loc = Loc(cov);
  int d, i,
    err = NOERROR,
    dim = shape->xdimprev;
  pgs_storage *pgs = cov->Spgs;
  bool 
    grid = Loc(cov)->grid,
    pgsnull = pgs == NULL;
  //PMI(cov);

  assert(cov->sub[PGS_LOC] != NULL);
  assert(dim == cov->sub[PGS_LOC]->xdimprev);
  assert(dim == Loc(cov)->timespacedim);
  
  if (Cshape->inverse == ErrCov) // && Cshape->nonstat_inverse == ErrCovNonstat)
    SERR1("support of the model is unknown. Use '%s' to determine the support",
	  CovList[TRUNCSUPPORT].nick);
 
  if (pgsnull) {
    if ((err = alloc_pgs(cov)) != NOERROR) return err;
    pgs = cov->Spgs;
    
    if ((pgs->v = (double*) CALLOC(dim, sizeof(double))) == NULL ||
	(pgs->y = (double*) CALLOC(dim, sizeof(double))) == NULL
      ) return ERRORMEMORYALLOCATION;

    pgs->zhou_c = 1.0;
    pgs->sum_zhou_c = pgs->sq_zhou_c = 0.0;
    pgs->n_zhou_c = 0;
  }
  
  // PMI(cov, "hi");
  assert(pts == cov->sub[PGS_LOC]);

  // selbst wenn zufaelliger Shape: 1x laufen lassen, ob 
  // Fehler auftauchen. Unter "Do" lassen sie sich nicht mehr
  // schoen abfangen.

  //  printf("A moments=%d\n", cov->mpp.moments);
  //  APMI(cov);
  if ((err = INIT(shape, cov->mpp.moments, S)) != NOERROR) return err; //gatter?
  // APMI(shape);


  assert(cov->mpp.moments >= 1);
  assert(shape->mpp.moments >= 1);

  //  printf("AB moments=%d\n", cov->mpp.moments);
  //APMI(pts->calling);
  if ((err = INIT(pts, 0, S)) != NOERROR) { // , 0, S statt ,1,S
    return err; 
  }
  //     printf("AC moments=%d\n", cov->mpp.moments);

  if (!grid) { // todo
    SERR("non-grid not programmed yet"); // meiste da; fehlt noch
    // welche Vektoren ich wirklich brauche
  }

  pgs->size =  grid ? intpow(2, dim) : loc->totalpoints;
  if (cov->role == ROLE_POISSON_GAUSS) {
    if (pgsnull &&
	(
	 (pgs->xgr[0] = (double*) CALLOC(dim * 3, sizeof(double))) == NULL ||
	 (pgs->pos = (int *) CALLOC(dim, sizeof(int))) ==NULL ||
	 (pgs->min = (int*) CALLOC(dim, sizeof(int))) == NULL ||
	 (pgs->max = (int*) CALLOC(dim, sizeof(int))) == NULL
	 )) return ERRORMEMORYALLOCATION;
    for (i=3, d=1; d<dim; d++, i+=3) pgs->xgr[d] = pgs->xgr[0] + i;
    // printf("%d\n", pgs->xgr[1] - pgs->xgr[0]); assert(false);
    if ((err = calculate_mass_gauss(cov)) != NOERROR) return err;
  } else if (hasMaxStableRole(cov)) {
    if (pgsnull &&
	(
	(pgs->single = (double*) CALLOC(pgs->size, sizeof(double)))==NULL ||
	(pgs->total = (double*) CALLOC(pgs->size, sizeof(double))) ==NULL ||
	(pgs->halfstepvector = (double*) CALLOC(dim, sizeof(double)))==NULL
	 ))
      return ERRORMEMORYALLOCATION;

    // APMI(cov);

    assert(!cov->deterministic || P0INT(PGS_NORMED));
    
   
    //         printf("init %f zhou=%f\n", pgs->totalmass, cov->mpp.zhou_c);
    if (P0INT(PGS_NORMED)) {
      if (R_FINITE(pts->mpp.unnormedmass)) { 
	// APMI(pts);

	// TO DO !! Richtig nur falls Huetchen derministisch ist !!!

	cov->mpp.maxheights[0] = pts->mpp.unnormedmass / shape->mpp.mMplus[1]; // korrekt, da pts * mass >= shape  !!

	//APMI(cov);

	//	
      } else {
	//PMI(pts->calling);
	//	
	// BUG; // not programmed yet
	// to do: funktioniert nur bei ausgewaehlten funktionen. Z.B. Gauss
	cov->mpp.maxheights[0] = shape->mpp.maxheights[0] * pts->mpp.maxheights[0];
	//     APMI(cov);
      }
    } else {
      if (R_FINITE(pts->mpp.maxheights[0])) {
	// eigentlich ein 1.0 / pts->mpp.minheight
	cov->mpp.maxheights[0] = pts->mpp.maxheights[0] * shape->mpp.maxheights[0];
      } else {
	BUG;
      }
    }			   
 
      // f_raw,rectangular must be always >= the shape function.
    // this is the fundamental condition
    if (!R_FINITE(cov->mpp.maxheights[0])) BUG;

    if ((cov->deterministic = shape->deterministic)) {
      if ((err = calculate_mass_maxstable(cov)) != NOERROR) return err;
      pgs->zhou_c = pgs->totalmass / shape->mpp.mMplus[1];

      //     printf("%f %f %f\n", pgs->zhou_c, pgs->totalmass, shape->mpp.mMplus[1]);assert(false);

    } else { // currently both cases are identical ... 
      pgs->zhou_c = RF_NA;
      if (pgs->cov == NULL) {
	cov_model *start=cov->calling;
	if (start == NULL) BUG;	
	//	APMI(start);
	while (start->calling != NULL && start->nr != PTS_GIVEN_SHAPE)
	  start = start->calling;	
	if (start->nr != PTS_GIVEN_SHAPE) {
	  if ((err = complete_copy(&(pgs->cov), cov)) != NOERROR) return err;
	  pgs->cov->calling = cov->calling;
	  pgs->cov->Spgs->cov = cov;
	  //PMI(cov);	  APMI(pgs->cov);
	} 
      }
      pgs->old_zhou = 0.0;
    }
    pgs->estimated_zhou_c = !cov->deterministic;


    //    PMI(cov);    assert(!pgs->estimated_zhou_c);

  } else {
    // APMI(cov);
    BUG;
  }
  


  if (CovList[shape->nr].nonstat_inverse == ErrInverseNonstat) {
    if (pts->nr != RECTANGULAR) {
      //       PMI(pts);
      // APMI(shape);
      warning("Inverse of shape function cannot be determined. Simulation speed  might be heavily decreased.");
    }
  }
  
 
  for (i=0; i<=cov->mpp.moments; i++) {
    //    printf("%d %f %f %d\n", i, pts->mpp.mM[i], shape->mpp.mMplus[i], cov->mpp.moments);

    //assert(pts->mpp.mMplus[0] == 1.0);
    
    cov->mpp.mM[i] = shape->mpp.mM[i] * pts->mpp.mMplus[0];
    cov->mpp.mMplus[i] = shape->mpp.mMplus[i] * pts->mpp.mMplus[0];
  }

  
   
  cov->rf = shape->rf;
  cov->origrf = false;
 
  // APMI(cov);

  return err;
}

int DrawCathegory(int size, double *single, 
		  double *total, 
		  bool calculate_elements, int *elmts) { 

  int i;
  double mass;
  mass = UNIFORM_RANDOM * total[size - 1];
  if (calculate_elements) {
    //printf("mass=%f\n", mass);
    
    for (i = 0; mass > total[i]; i++);
    if (i > 0) mass -= total[i-1];
    *elmts = floor(mass / single[i]);
    //printf("draw: i=%d tot=%f mass=%f sing=%f elm=%d\n", i, total[i], mass, single[i], *elmts);
    return i;
  } else {
    // falls nicht calculate_elements, so sind letztere unwichtig,
    // d.h. zu single wird nur der Faktor 1.0 multipliziert.
    // dies ist aber der Fall nur falls !grid
    return CeilIndex(mass, total, size);
  }
}

static double gauss_eps = 1e-10;
void do_pgs_gauss(cov_model *cov, gen_storage *S) {
  pgs_storage *pgs = cov->Spgs;
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  location_type *loc = Loc(cov);
  //  window_info *w = &(S->window);
  long i;
  int d, 
    *min = pgs->min, // dummy
    *max = pgs->max, // dummy
    *pos = pgs->pos, // dummy
     dim = shape->xdimprev;
  double value,
    total = 0.0,
    *x = pgs->x, // dummy or (future) static local
    *y = pgs->y, // dummy or (future) static local
    //    *supportmin = pgs->supportmin,
    //   *supportmax = pgs->supportmax,
    // *pgstotal = pgs->total, // in
    //*single = pgs->single, // in
    **xgr = pgs->xgr, // in
    *v = pgs->v;    
  bool 
    grid = Loc(cov)->grid;


  if (!cov->deterministic) {
    DO(shape, S); 
    DORANDOM(pts, cov->q);  // cov->q nur dummy. Wird ueberschrieben
    if (cov->role == ROLE_POISSON_GAUSS || !grid) { 
      if (calculate_mass_gauss(cov) != NOERROR) 
	error("unexpected error in 'do_pts_given_shape' (maxstable)");
    } else BUG;
  }

  // printf("name = %s\n", CovList[pts->nr].name);
  //  APMI(pts);
  
  VTLG_R(NULL, pts, v); 
  i = UNIFORM_RANDOM * pgs->totalmass;
  
  // to determine the points that contribute to the value of the 
  // random field by the shape function centered at cov->q
  //printf("cat=%d el=%d\n", i, elmts);
  if (loc->grid) {
    NONSTATINVERSE_D(&gauss_eps, pts, x, y);
    if (ISNAN(x[0]) || x[0] > y[0]) BUG;

    //    printf("extremes %f %f %f\n", eps, x[0], y[0]);

    for (d=0; d<dim; d++) {
      //printf("e1 %f %f\n", x[d], y[d]);
      int which = i % (int) xgr[d][XLENGTH];
      i /= xgr[d][XLENGTH];
      cov->q[d] = xgr[d][XSTART] + xgr[d][XSTEP] * which + v[d];
      // to determine the points that contribute to the value of the 
      // density function at cov->q
      min[d] = ceil((cov->q[d] - y[d] - xgr[d][XSTART]) / xgr[d][XSTEP]);
      max[d] = (cov->q[d] - x[d] - xgr[d][XSTART]) / xgr[d][XSTEP];

      // printf("%d %d delat=%f %f y=%f %f step=%f\n", min[d], max[d], 
      //     cov->q[d] - y[d] - xgr[d][XSTART],
      //     cov->q[d] - x[d] - xgr[d][XSTART],
      //     y[d], x[d], xgr[d][XSTEP]);

      // PMI(cov);
      assert(xgr[d][XSTEP] > 0);

      if (min[d] < 0) min[d] = 0;
      if (max[d] >= xgr[d][XLENGTH]) max[d] = xgr[d][XLENGTH] -1;
      if (min[d] > max[d]) {
	// no point contributed to cov->q -- this might be possible if 
	// the simulating grid goes beyond the actual grid and 	

	//printf("redo d=%d x=%f y=%f q=%f i=%d d=%d range=%d %d\n",
	// d, x[d], y[d], cov->q[d], i, d, min[d], max[d]);
	do_pgs_gauss(cov, S);
	pgs->log_density = RF_INF;
	return;
      } // else printf("ok q=%f i=%d\n", cov->q[d], i, d, min[d], max[d]);
      pos[d] = min[d];
      y[d] = x[d] = cov->q[d] - (xgr[d][XSTART] + pos[d] * xgr[d][XSTEP]);
    }
   
    while (true) {
      VTLG_D(y, pts, &value); 
      total += value;
      
     //printf("y=%f v=%f q=%f val=%f tot=%f min=%d max=%d pos=%d xgr.len=%d\n", 
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
  pgs->log_density = log(total / pgs->totalmass);

  //printf("total=%f %f grid=%d dim=%d\n", total, pgs->totalmass,loc->grid,dim);

  assert(R_FINITE(pgs->log_density));
  // assert(false);

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
  dim = shape->xdimprev;			 \
  flat = pgs->flat;				 \
  assert(pgs != NULL);				 \
  assert(cov->deterministic ||			 \
	 (pgs->cov!=NULL && pgs->cov->Spgs!=NULL \
	  && pgs->cov->Spgs->cov==cov));

#define SWAP_PGS				\
  cov = cov->Spgs->cov;				\
  SET_COV

 
//static double sumArea = 0.0;
void do_pgs_maxstable(cov_model *cov, gen_storage *S) {
  pgs_storage *pgs = NULL;
  cov_model *shape = NULL,
    *pts = NULL;
  location_type *loc = NULL;
  int i, d, elmts, mcmc, dim, err;
  double 
    *single = NULL, // in
    *total = NULL,   // in
    *halfstepvector = NULL, // in
    *x = NULL, // dummy
    *v = NULL; // dummy
  bool flat;
  assert(cov->calling != NULL);

  if (cov->deterministic) {
    SET_COV;
  } else {
    // assert(PARAM0(shape->sub[0], POWSCALE) == PARAM0(pts, LOC_SCALE));
    double cmaxDmu;
    SWAP_PGS;    
    for (mcmc=0; mcmc<GLOBAL.extreme.mcmc_zhou; mcmc++) {
      if ((err = REINIT(cov, cov->mpp.moments, S)) != NOERROR) BUG;
      DO(shape, S);     
      if (calculate_mass_maxstable(cov) != NOERROR) 
	error("unexpected error in 'do_pts_given_shape' (maxstable)");

      cmaxDmu = pgs->totalmass / shape->mpp.mMplus[1];
      if (pgs->n_zhou_c < GLOBAL.extreme.max_n_zhou) {
	pgs->n_zhou_c++;
	pgs->sum_zhou_c += cmaxDmu;
	pgs->sq_zhou_c += cmaxDmu * cmaxDmu;
	// sumArea += getArea( shape->sub[0]->Spolygon->P);
	//if (pgs->n_zhou_c % 1000 == 0) printf("mean area = %f\n", sumArea / (double) pgs->n_zhou_c);
      }

      //    PMI(cov);
      // PMI(pts);
      
      if (false)
      printf("//cmaxDmu %4.1f; mass=%4.1f m1=%4.1e beta=%4.1f (%4.6f, %4.6f) (%4.6f, %4.6f) %4.3f %e\n",	     
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
	//printf(". %f %f\n", old_zhou, cmaxDmu);
	pgs->old_zhou = cmaxDmu;
	SWAP_PGS;
      } //else printf("* %f %f\n", old_zhou, cmaxDmu);
    } // mcmc
    
    SWAP_PGS;
    cov_model *prev=cov->calling;		
    assert(prev != NULL);				
    if (prev->key != NULL) prev->key = cov;		
    else if (prev->sub[0] != NULL) prev->sub[0] = cov;	
    else if (prev->sub[1] != NULL) prev->sub[1] = cov;	
    else error("structure mismatch");			

    // APMI(cov);


    // printf("\n");
    //PMI(pts);
    // printf("%f %f\n", PARAM(pts, 1)[0] *PARAM(pts, 1)[1] * 4 , pts->mpp.mMplus[0]);
    
    //assert(fabs(PARAM(pts, 1)[0] *PARAM(pts, 1)[1] * 4 - pts->mpp.mMplus[0]) < 1e-14 *  pts->mpp.mMplus[0]);

    //   PMI(save_pts);
    //APMI(pts);

    DORANDOM(pts, cov->q); // cov->q ist hier nur dummy    

    // if (!newshape) printf("failed\n");
    //printf("old %f %f %f\n", (double) pgs->old_zhou, (double) pgs->totalmass, PARAM0(shape->sub[0], POWSCALE));

    // assert(pgs->old_zhou == pgs->totalmass);
    // if (pgs->old_zhou > 1e6 && PARAM0(shape->sub[0], POWSCALE)>1) APMI(cov);  //    assert(shape->sub[0]->nr == POWER_DOLLAR && pts->nr == LOC);
    // assert(PARAM0(shape->sub[0], POWSCALE) == PARAM0(pts, LOC_SCALE));


    //
    if (false)
      printf("//zhou %4.1f %4.3f %ld mean=%4.3f, shape=%s %f \n", (double) pgs->sum_zhou_c, pgs->totalmass, pgs->n_zhou_c, 
	     (double) pgs->sum_zhou_c / (double) pgs->n_zhou_c,
	     NICK(shape->sub[0]), PARAM0(shape->sub[0], POWSCALE)
	     ); // assert(pgs->n_zhou_c <= 500);
  }

  i = DrawCathegory(pgs->size, single, total, loc->grid, &elmts);
 
  if (flat) for (d=0; d<dim; d++) x[d] = 0.0; 
  else for (d=0; d<dim; d++) x[d] = halfstepvector[d];

  if (dim <= 2) {
    switch (i) {
    case PGS_VOXEL : 
      break;
    case PGS_CORNER :
      for (d=0; d<dim; d++) x[d] = RF_INF; 
      break;
    case PGS_2D_SIDE1 : x[0] = RF_INF;
      break;
    case PGS_2D_SIDE2 : x[1] = RF_INF;
      break;
    default : BUG;
    }
  } else if (dim == 3) {
    switch (i) {
    case PGS_VOXEL : 
      break;
    case PGS_CORNER :
      for (d=0; d<dim; d++) x[d] = RF_INF; 
      break;
    case PGS_3D_SIDE1 : x[0] = RF_INF;
      break;
    case PGS_3D_SIDE2 : x[1] = RF_INF;
      break;
    case PGS_3D_SIDE3 : x[2] = RF_INF;
      break;
    case PGS_3D_EDGE1 : x[1] = x[2] = RF_INF;
      break;
    case PGS_3D_EDGE2 : x[0] = x[2] = RF_INF;
      break;
    case PGS_3D_EDGE3 : x[0] = x[1] = RF_INF;
      break;
    default : BUG;
    }    
  } else BUG;
 
  //APMI(cov);
  //printf("flat = %d %ld %s i=%d %d \n", flat, x, NICK(pts), i, PGS_VOXEL);
 
  if (flat) { 
    // to do: programmierfehler!
    // frage: lohnt sich 'flat=TRUE' ueberhaupt noch hinsichtlich Aufwand
    //        und Rechenzeitersparnis?
    // Bsp das nicht funktioniert fuer flat=FLAT_UNDETERMINED
    //    x _ seq(1, 10, 0.01)
    //      z _ RFsimulate(RPsmith(RMgauss(), xi=0), x, x, grid=T, print=20)
    //      Print(exp(-exp(range(-z[[1]]))))
    //      plot(z); X11(); hist(z[[1]], freq=FALSE, 50)
    //			curve(exp(-x-2) * exp(-exp(-x-2)), -5, 2, add=TRUE)
    if (i != PGS_VOXEL) VTLG_R(x, pts, v); 
    //  printf("v=%f x=%f %d %d %d \n", v[0],  x[0], i, PGS_VOXEL,  PGS_CORNER);
    for (d = 0; d<dim; d++)
      if (x[d] == 0.0) 
	v[d] = loc->xgr[d][XSTEP] * UNIFORM_RANDOM - halfstepvector[d];
  } else { 
    //    printf("here %f %f %f\n", x[0], x[1], x[2]);
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
  // printf("i=%d %f %f %f v=%f %f %f\n", i, x[0], x[1], x[2], v[0], v[1], v[2]);

  // APMI(cov);
 
  //printf("i=%d %f %f v=%f %f\n", i, x[0], x[1], v[0], v[1]);

  for (d=0; d<dim; d++) {
    cov->q[d] = loc->xgr[d][XSTART] + v[d];
    // 
    // printf("q=%f %f %f %d\n", cov->q[d], loc->xgr[d][XSTART], v[d], flat);
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
    //   printf("D=%d %f %f %d\n", d, cov->q[d], v[d], elmts);

  }
 
  if (flat) for (d = 0; d<dim; d++) if (x[d] == 0.0)  v[d] = 0.0;
 
  if (CovList[pts->nr].logD != ErrCov) {
    // printf("log!\n");
    VTLG_DLOG(v, pts, &(pgs->log_density)); // hier muss allgemein die nicht-normierte Dichte stehen und stattdessen durch das Mittel des Volumens geteilt  werden
   

    
    //printf("log density v=%f, %f, %f %e\n", v[0],v[1],v[2], pgs->log_density);

    //printf("log_dens %f\n", pgs->log_density);
    //  pgs->log_density -= log(pgs->totalmass); // ueberprueft: normierung
    // zu einer WK'Dichte.
  } else {
    double density;
    VTLG_D(v, pts, &density);

    //  printf("dens %f\n", density);
      //  printf("density %e v=%f %f %f\n", density, v[0], v[1], v[2]);

    pgs->log_density = log(density);
  }

 
  if (!R_FINITE(pgs->log_density)) {
    //PMI(cov);
    BUG;
  }

  //assert(pgs->log_density == 0);

  //  PMI(cov);
  
}



void do_pts_given_shape(cov_model *cov, gen_storage *S) {
  // muss zu allerst stehen, da cov sich aendern kann!
  if (cov->role == ROLE_POISSON_GAUSS) {
    do_pgs_gauss(cov, S);
  } else if (hasMaxStableRole(cov)) { // todo: sauber trennen, wann 
    do_pgs_maxstable(cov, S); 
  } else {
    PMI(cov); BUG; //
  }

  // do_gs_maxstable might break links of the cov structure:
  cov_model *prev=cov->calling;		
  assert(prev != NULL);				
  if (prev->key != NULL) cov = prev->key;		
  else if (prev->sub[0] != NULL) cov = prev->sub[0];	
  else if (prev->sub[1] != NULL) cov = prev->sub[1];	
  else error("structure mismatch");			
  
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  pgs_storage *pgs = cov->Spgs;
  int d,
    dim = shape->xdimprev;
  double 
    eps,
    *x = pgs->x,
    *y = pgs->y;

  if (cov->role == ROLE_POISSON_GAUSS) {
    eps = GLOBAL.mpp.about_zero * exp(pgs->log_density);
  } else if (hasMaxStableRole(cov)) { // todo: sauber trennen, wann 
    // max-stable, smith etc
    eps = pgs->currentthreshold;

    //printf("eps= %f\n", eps);

    if (!R_FINITE(pgs->log_density)) {
      //PMI(cov->Spgs->cov);   PMI(cov);
      BUG;
    }
    if (cov->loggiven)  eps += pgs->log_density;
    else eps *= exp(pgs->log_density);
  } else BUG;

  
  //printf("here %f %d\n", eps, cov->loggiven);
  // APMI(cov);
   
  if (cov->loggiven) { NONSTATLOGINVERSE(&eps, shape, x, y); }
  else NONSTATINVERSE(&eps, shape, x, y);

  // assert(false);
  // warum fkt obiges nicht bei gauss??
  //  PMI(shape); 
  //printf("eps=%f x=%f %d %f %f\n", eps, *x, cov->loggiven,
  //
  //	 pgs->currentthreshold, pgs->log_density);


  if (ISNAN(x[0]) || x[0] > y[0]) {
    //double eps_neu = eps / cov->mpp.maxheights[0]; // warum ?? 29.12.2013
    double eps_neu = eps; //  29.12.2013
    if (cov->loggiven) { BUG; }
    else NONSTATINVERSE_D(&eps_neu, pts, x, y); 
    
    if (ISNAN(x[0])  || x[0] > y[0]) BUG;
 
    //printf("eps=%f %f %f x=(%e %e %e) y=(%e %e %e) q=(%f %f %f)\n",
    //     eps, pgs->currentthreshold,  exp(pgs->log_density),
    //     x[0], x[1], x[2], y[0], y[1], y[2],
    //     cov->q[0], cov->q[1], cov->q[2]);
  }

  //  printf("Xx=%f %d %d %f %f\n", *x, ISNAN(x[0]),cov->loggiven, pgs->currentthreshold, pgs->log_density);


  //APMI(cov);
  
  for (d=0; d<dim; d++) {
 
    pgs->supportmin[d] = cov->q[d] - y[d]; // 4 * for debugging...
    pgs->supportmax[d] = cov->q[d] - x[d];

    //     printf("do d=%d q=%f x=%f y=%f eps=%e thr=%e %f %f\n",
    //	  	  d, cov->q[d], x[d], y[d], eps, pgs->currentthreshold, pgs->supportmin[d], pgs->supportmax[d]);

  if (ISNAN(pgs->supportmin[d]) || ISNAN(pgs->supportmax[d]) ||
	pgs->supportmin[d] > pgs->supportmax[d]) {
      //      printf("do d=%d q=%f min=%f max=%f eps=%e thr=%e supp=(%e, %e)\n",    d, cov->q[d], x[d], y[d], eps, pgs->currentthreshold,  pgs->supportmin[d], pgs->supportmax[d]);
      // APMI(shape);

// BUG;
    }
  
    //assert(pgs->supportmin[d] <= pgs->supportmax[d]);
  }
  // assert(false);
  

  //  BUG;
  //   APMI(cov);

  cov->fieldreturn = shape->fieldreturn;
  cov->origrf = false;
}


void range_pts_given_shape(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[PGS_RATIO] = 0;
  range->max[PGS_RATIO] = 1;
  range->pmin[PGS_RATIO] = 0;
  range->pmax[PGS_RATIO] = 1;
  range->openmin[PGS_RATIO] = false;
  range->openmax[PGS_RATIO] = false; 

  range->min[PGS_FLAT] = -1;
  range->max[PGS_FLAT] = 1;
  range->pmin[PGS_FLAT] = -1;
  range->pmax[PGS_FLAT] = 1;
  range->openmin[PGS_FLAT] = false;
  range->openmax[PGS_FLAT] = false; 

  range->min[PGS_INFTY_SMALL] = 0;
  range->max[PGS_INFTY_SMALL] = 1;
  range->pmin[PGS_INFTY_SMALL] = 0;
  range->pmax[PGS_INFTY_SMALL] = 1;
  range->openmin[PGS_INFTY_SMALL] = false;
  range->openmax[PGS_INFTY_SMALL] = false; 

  range->min[PGS_NORMED] = 0;
  range->max[PGS_NORMED] = 1;
  range->pmin[PGS_NORMED] = 0;
  range->pmax[PGS_NORMED] = 1;
  range->openmin[PGS_NORMED] = false;
  range->openmax[PGS_NORMED] = false; 

  range->min[PGS_ISOTROPIC] = 0;
  range->max[PGS_ISOTROPIC] = 1;
  range->pmin[PGS_ISOTROPIC] = 0;
  range->pmax[PGS_ISOTROPIC] = 1;
  range->openmin[PGS_ISOTROPIC] = false;
  range->openmax[PGS_ISOTROPIC] = false; 
}



void standard_shape(double *x, cov_model *cov, double *v) { 
  cov_model *shape = cov->sub[PGS_FCT];
  NONSTATCOV(x, cov->q, shape, v);
}

void logstandard_shape(double *x, cov_model *cov, double *v, double *sign) { 
  cov_model *shape = cov->sub[PGS_FCT];
  LOGNONSTATCOV(x, cov->q, shape, v, sign);
}

int check_standard_shape(cov_model *cov) {
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  int err, role,
    dim = cov->tsdim; 
  

  if (cov->q == NULL) {
    // zwingend CALLOC !
    if ((cov->q  = (double*) CALLOC(sizeof(double), dim)) == NULL)
      return ERRORMEMORYALLOCATION;
    cov->qlen = dim;
  }

  if (cov->xdimprev != cov->xdimown || cov->xdimprev != cov->tsdim) 
    return ERRORDIM;
  
  if (hasPoissonRole(cov)) {
    role = ROLE_POISSON;
  } else if (hasMaxStableRole(cov)) {
    role = ROLE_MAXSTABLE;
  } else ILLEGAL_ROLE;


  if ((err = CHECK(shape, dim, dim, ShapeType, XONLY, CARTESIAN_COORD, 
		     SCALAR, role)) != NOERROR)  return err;
  setbackward(cov, shape);
  
  if (!shape->deterministic) 
    SERR1("random shapes for '%s' not allowed yet", NICK(cov));

  if (pts != NULL) {
    if ((err = CHECK_R(pts, dim)) != NOERROR) return err;
  }

  return NOERROR;
}

int struct_standard_shape(cov_model *cov, cov_model **newmodel){
  cov_model *shape = cov->sub[PGS_FCT];
  //  int    err = NOERROR;
  //dim = shape->xdimprev;

  // printf("ttt\n");

  ASSERT_NEWMODEL_NULL;


  if (shape->role != ROLE_POISSON && shape->role != ROLE_MAXSTABLE)
    ILLEGAL_ROLE;

  //   printf("here %s\n", NICK(cov->sub[PGS_LOC]));
  
  //  APMI(cov);
       
  //if ((err = CHECK_R(cov->sub[PGS_LOC], cov->tsdim)) != NOERROR) return err;
  
  return NOERROR;
}
 

int init_standard_shape(cov_model *cov, gen_storage *S) {  
  cov_model *shape = cov->sub[PGS_FCT];
  //  location_type *loc = Loc(cov);

  //APMI(cov);

  assert(cov->sub[PGS_LOC] != NULL);
  location_type *loc = Loc(cov);
  int d,
    dim = shape->xdimprev,
    err = NOERROR;
  pgs_storage *pgs = cov->Spgs;

  //APMI(cov);

  assert(shape->xdimprev == Loc(cov)->timespacedim);

  if (pgs == NULL) {
    if ((err = alloc_pgs(cov)) != NOERROR) return err;
    pgs = cov->Spgs;
    if ((pgs->localmin = (double*) CALLOC(dim, sizeof(double))) == NULL ||
	(pgs->localmax = (double*) CALLOC(dim, sizeof(double))) == NULL ||
	(pgs->minmean = (double*) CALLOC(dim, sizeof(double))) == NULL ||
	(pgs->maxmean = (double*) CALLOC(dim, sizeof(double))) == NULL)
      return ERRORMEMORYALLOCATION;
  }

 
  // selbst wenn zufaelliger Shape: 1x laufen lassen, ob 
  // Fehler auftauchen. Unter "Do" lassen sie sich nicht mehr
  // schoen abfangen.
  // ROLE_ASSERT(ROLE_POISSON || cov->role == ROLE_MAXSTABLE);

  //PMI(cov->calling, -1);
  
  assert(cov->mpp.moments == hasMaxStableRole(cov));
  if ((err = INIT(shape, cov->mpp.moments, S)) != NOERROR) return err; //gatter?

  cov_model *u = cov->sub[PGS_LOC]; 
  //PMI(cov);
  assert(u->nr == UNIF);
  double 
    *min = PARAM(u, UNIF_MIN),
    *max = PARAM(u, UNIF_MAX),
    *x = pgs->minmean, // !! wird gespeichert
    *y = pgs->maxmean; // !!
  
  // try out whether it works:
  NONSTATINVERSE(ZERO, shape, x, y);
  if (ISNAN(x[0]) || x[0] > y[0])
    SERR1("inverse of '%s' unknown", NICK(shape));
  GetDiameter(loc, pgs->localmin, pgs->localmax, 
	      pgs->supportcentre); // last is only a dummy
  pgs->totalmass = 1.0;
  for (d=0; d<dim; d++) {
    min[d] = pgs->localmin[d] - y[d];
    max[d] = pgs->localmax[d] - x[d];
    if (!(R_FINITE(min[d]) && R_FINITE(max[d])))
      SERR1("simulation window does not have compact support. Should '%s' be used?", CovList[TRUNCSUPPORT].nick);
    pgs->totalmass *= max[d] - min[d];
  }
  
  if (cov->role == ROLE_POISSON) {
    pgs->log_density = 0.0;
  } else {    
    pgs->log_density = 0.0;
    pgs->zhou_c = pgs->totalmass / shape->mpp.mMplus[1];
    cov->mpp.maxheights[0] = shape->mpp.maxheights[0];
    pgs->estimated_zhou_c = !cov->deterministic;
    if (pgs->estimated_zhou_c) SERR("random shapes in standard approach not coded yet -- please contact author");
  }


  cov->rf = shape->rf;
  cov->origrf = false;
  cov->fieldreturn = shape->fieldreturn;
 
  // APMI(cov);

  return NOERROR;
}


void do_standard_shape(cov_model *cov, gen_storage *S) {
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC]; 
  assert(cov->sub[PGS_LOC] != NULL);
  pgs_storage *pgs = cov->Spgs;
  double *x = pgs->x,
    *y = pgs->xstart;
  int d,
    dim = shape->xdimprev;

  DO(shape, S);
  DORANDOM(pts, cov->q);


  //printf("q=%f\n", cov->q[0]);
  // APMI(cov->calling);

  assert(!shape->fieldreturn);


  //PMI(shape, "standard");
  NONSTATINVERSE(ZERO, shape, x, y);
  if (ISNAN(x[0])|| x[0] > y[0]) BUG;

  for (d=0; d<dim; d++) {


    //    printf("do d=%d q=%f x=%f y=%f thr=%e\n",
    //	  d, cov->q[d], x[d], y[d], pgs->currentthreshold);
    pgs->supportmin[d] = cov->q[d] - y[d]; // 4 * for debugging...
    pgs->supportmax[d] = cov->q[d] - x[d];

    assert(pgs->supportmin[d] != NA_INTEGER && 
	   pgs->supportmax[d] != NA_INTEGER);
    
    assert(pgs->supportmin[d] <= pgs->supportmax[d]);
  }  

  pgs->log_density = 0.0;

}



 






#define STAT_SHAPE_FCT 0
void stationary_shape(double *x, cov_model *cov, double *v) { 
  cov_model *shape = cov->sub[STAT_SHAPE_FCT];
  FCTN(x, shape, v);
  //  if (cov->q[0] < 1) {
  //   printf("x=%f q=%f v=%f\n", x[0], cov->q[0], v[0]);
  //   assert(false);
  // } 
  //  APMI(cov);
}

void logstationary_shape(double *x, cov_model *cov, double *v, double *sign) { 
  cov_model *shape = cov->sub[STAT_SHAPE_FCT];
  LOGCOV(x, shape, v, sign);
}

int check_stationary_shape(cov_model *cov) {
  cov_model *shape = cov->sub[STAT_SHAPE_FCT];
  int err, role,
    dim = cov->tsdim; 
  
  if (cov->xdimprev != cov->xdimown || cov->xdimprev != cov->tsdim) 
    return ERRORDIM;
  
  if (cov->role == ROLE_GAUSS) {
    role = isGaussProcess(shape) ? ROLE_GAUSS
      : shape->nr == BINARYPROC ? ROLE_GAUSS
      : ROLE_UNDEFINED; 
    ASSERT_ROLE_DEFINED(shape);
  } else if (hasPoissonRole(cov)) {
    role = ROLE_POISSON;
  } else if (hasMaxStableRole(cov)) {
    role = ROLE_MAXSTABLE;
  } else ILLEGAL_ROLE;


  //printf("here\n");
  
  if ((err = CHECK(shape, dim, dim, ProcessType, XONLY, CARTESIAN_COORD, 
		     SCALAR, role)) != NOERROR)  return err;
  setbackward(cov, shape);

  return NOERROR;
}

int struct_stationary_shape(cov_model *cov, cov_model **newmodel){
  cov_model *shape = cov->sub[STAT_SHAPE_FCT];
  //  location_type *loc = Loc(cov);

  ASSERT_NEWMODEL_NULL;

  if (shape->role != ROLE_POISSON && shape->role != ROLE_MAXSTABLE)
    ILLEGAL_ROLE;

  //printf("here\n");

  return NOERROR;
}


int init_stationary_shape(cov_model *cov, gen_storage *S) {  
  cov_model *shape = cov->sub[STAT_SHAPE_FCT];
  int d, i,
    err = NOERROR,
    dim = shape->xdimprev;
 
  //PMI(cov);

  assert(dim == Loc(cov)->timespacedim);
  if ((err = alloc_pgs(cov)) != NOERROR) return err;
  pgs_storage *pgs = cov->Spgs;

  // selbst wenn zufaelliger Shape: 1x laufen lassen, ob 
  // Fehler auftauchen. Unter "Do" lassen sie sich nicht mehr
  // schoen abfangen.

  
  assert(cov->mpp.moments >= 1);
  if ((err = INIT(shape, 1, S)) != NOERROR) return err; //gatter?
  assert(shape->mpp.moments >= 1);
  for (i=0; i<=cov->mpp.moments; i++) {
    cov->mpp.mM[i] = shape->mpp.mM[i];
    cov->mpp.mMplus[i] = shape->mpp.mMplus[i];
  }

   
  pgs->zhou_c = 1.0 / cov->mpp.mMplus[1]; // passt fuer binary, und auch fuer 
  if (!R_FINITE(pgs->zhou_c))
    SERR1("max height of '%s' not finite", NICK(shape));
  pgs->estimated_zhou_c = false; 
  if (!cov->deterministic) SERR("not deterministic shapes in stationary modelling -- please contact author");
  
  pgs->log_density = 0; 
   
  //  printf("dims %d %d\n", cov->xdimown, shape->xdimprev);

  for (d=0; d<dim; d++) {
    pgs->supportmin[d] = RF_NEGINF; // 4 * for debugging...
    pgs->supportmax[d] = RF_INF;
  }

  cov->mpp.maxheights[0] = shape->mpp.maxheights[0];
  cov->rf = shape->rf;
  cov->origrf = false;
  cov->fieldreturn = shape->fieldreturn;
  if (!cov->fieldreturn) BUG;


  // APMI(cov);

  return NOERROR;
}


void do_stationary_shape(cov_model *cov, gen_storage *S) {
  cov_model *shape = cov->sub[STAT_SHAPE_FCT]; 
  DO(shape, S);
  cov->mpp.maxheights[0] = shape->mpp.maxheights[0];
  assert(shape->fieldreturn);
}


