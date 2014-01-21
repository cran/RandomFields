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


#define PGS_RATIO 0
#define PGS_FLAT 1

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


int addUnifModel(cov_model VARIABLE_IS_NOT_USED *cov, double radius, cov_model **newmodel) {
  addModel(newmodel, UNIF);
  kdefault(*newmodel, UNIF_MIN, -radius);
  kdefault(*newmodel, UNIF_MAX, radius);
  return NOERROR;
}


int alloc_pgs(cov_model *cov) {
  return alloc_pgs(cov, cov->xdimown);
}

int alloc_pgs(cov_model *cov, int dim) { // all what is necessary for dompp
  pgs_storage *pgs = NULL;

  if (cov->Spgs != NULL) PGS_DELETE(&(cov->Spgs));
  if ((pgs = cov->Spgs = (pgs_storage*) MALLOC(sizeof(pgs_storage))) == NULL) 
    return ERRORMEMORYALLOCATION;
  PGS_NULL(pgs);

  //printf("here\n");
  
  if ((pgs->supportmin = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->supportmax = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->supportcentre = (double*) CALLOC(dim, sizeof(double))) == NULL || 
      
      (pgs->gridlen = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->end = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->start = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->delta = (int*) CALLOC(dim, sizeof(int))) == NULL ||
      (pgs->nx = (int*) CALLOC(dim, sizeof(int))) == NULL ||
   
      (pgs->xstart = (double*) CALLOC(dim, sizeof(double))) == NULL || 
      (pgs->x = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->inc = (double*) CALLOC(dim, sizeof(double))) == NULL)
    return ERRORMEMORYALLOCATION;

  // printf("here end\n");
  return NOERROR;
}




void pts_given_shape(double *x, cov_model *cov, double *v) { 
  cov_model *shape = cov->sub[PGS_FCT];
  NONSTATCOV(x, cov->q, shape, v);
  //  if (cov->q[0] < 1) {
  //   printf("x=%f q=%f v=%f\n", x[0], cov->q[0], v[0]);
  //   assert(false);
  // } 
  //  APMI(cov);
}

void logpts_given_shape(double *x, cov_model *cov, double *v, double *sign) { 
  cov_model *shape = cov->sub[PGS_FCT];
  LOGNONSTATCOV(x, cov->q, shape, v, sign);
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


  //printf("here\n");
  
  // isotropy is the property
  if ((err = CHECK(shape, dim, dim, ShapeType, XONLY, ISOTROPIC,
		     SCALAR, role)) != NOERROR)  return err;
  // but it must be treated as if it was non-isotropic
  if ((err = CHECK(shape, dim, dim, ShapeType, XONLY, CARTESIAN_COORD,
		     SCALAR, role)) != NOERROR)  BUG;
  setbackward(cov, shape);

  //if (!shape->deterministic) return ERRORNOTPROGRAMMED; // Dichte muss nicht normiert sein und stattdessen durch das mittlere Volumen dividiert werden.

  if (pts != NULL) {

    //printf("pts %s %s\n", CovList[pts->nr + 1].nick, CovList[pts->nr + 1].name);
    //PMI(pts);

    if ((err = CHECK_R(pts, dim)) != NOERROR) return err;
  }

  //  printf("done\n");
 
  return NOERROR;
}



int struct_pts_given_shape(cov_model *cov, cov_model **newmodel){
  cov_model *shape = cov->sub[PGS_FCT];
  //  location_type *loc = Loc(cov);
  int err = NOERROR;

  if (newmodel != NULL) BUG;
  if (cov->Spgs != NULL) free(cov->Spgs);

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
    //   *single = pgs->single,   // out
    //  *total = pgs->total,     // out
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
    if (ISNA(x[0]) || x[0] > y[0]) 
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
	y[d] /= sqrt((double)dim); // falls nicht runif
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
	  BUG; assert(false);
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
  double voxels,
    *single = pgs->single, // out
    *total = pgs->total,   // out
    *halfstepvector = pgs->halfstepvector, //out
    *x = pgs->x;
  //  bool flat = false;
  int d,
    dim = cov->tsdim;
 
  if (shape->role == ROLE_POISSON) {
    BUG;
    // to do : felix's paper; derzeit einfach wie max-stable
    //   if (pts->nr !=    SERR("currently, only simple domain fields are allowed that are based on Poisson processes");

  } // else {

  for (d=0; d<dim; d++) halfstepvector[d] = 0.0; 

  //  printf("hier calculate_mass\n");
  VTLG_D(halfstepvector, pts, &(pgs->value_orig));

  //  PMI(cov, "calculate max");
  
  int flat =P0INT(PGS_FLAT);
  //  printf("%s %d %d %f\n",NICK(cov), PGS_FLAT, ((int*)cov->p[PGS_FLAT])[0], pgs->value_orig );// APMI(cov); assert(false);
    
  if (flat == FLAT_UNDETERMINED) {
    // flat == im Kerngebiet wird konstant simuliert; ausserhalb dann
    // abfallend
    if (loc->grid) {
      double v;
      for (d=0; d<dim; d++) halfstepvector[d] = 0.5 * loc->xgr[d][XSTEP];     
      VTLG_D(halfstepvector, pts, &v);
      double threshold = 
	pgs->value_orig == RF_INF//&& cov->p[PGS_RATIO][0] == 0.0 
	? RF_INF
	: pgs->value_orig * P0(PGS_RATIO);
      pgs->flat =  threshold < v;
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
    single[PGS_VOXEL] = pgs->value_orig;
    for (d=0; d<dim; d++) single[PGS_VOXEL] *= loc->xgr[d][XSTEP];
  } else {
    //    APMI(cov);
    //       printf("two sided %f %f %f\n", halfstepvector[0], halfstepvector[1], halfstepvector[2]); //assert(false);
    VTLG_P2SIDED(NULL, halfstepvector, pts, single + PGS_VOXEL);
    //   printf("two sided done %f %s\n", single[PGS_VOXEL], NICK(pts)); assert(false);
  }
  voxels = 1.0;
  for (d=0; d<dim; d++) voxels *= loc->xgr[d][XLENGTH] - 1.0;  
  total[PGS_VOXEL] = single[PGS_VOXEL] * voxels;
  
  for (d=0; d<dim; d++) x[d] = RF_INF;  
  VTLG_P2SIDED(NULL, x, pts, single + PGS_CORNER);
  assert(single[PGS_CORNER] == 1.0); // nur wenn normiert
  total[PGS_CORNER] = 1.0 + total[PGS_VOXEL];

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

  if (!R_FINITE(pgs->totalmass = total[pgs->size - 1])) {
    // int i; for (i=0; i<pgs->size; i++) printf("i=%d total=%f\n", i, total[i]);
    SERR("the total intensity mass is not finite");
  }

  return NOERROR;
}



int init_pts_given_shape(cov_model *cov, storage *S) {  
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  cov_fct *Cshape = CovList + shape->nr;
  location_type *loc = Loc(cov);
  int d, i,
    err = NOERROR,
    dim = shape->xdimprev;
  bool grid = Loc(cov)->grid;

  //PMI(cov);

  assert(cov->sub[PGS_LOC] != NULL);
  assert(dim == cov->sub[PGS_LOC]->xdimprev);
  assert(dim == Loc(cov)->timespacedim);
  
  if (Cshape->inverse == ErrCov)
    SERR1("support of the model is unknown. Use '%s' to determine the support",
	  CovList[TRUNCSUPPORT].nick);
 
  if ((err = alloc_pgs(cov)) != NOERROR) return err;
  pgs_storage *pgs = cov->Spgs;

  if ((pgs->v = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->y = (double*) CALLOC(dim, sizeof(double))) == NULL
      ) return ERRORMEMORYALLOCATION;
  
  
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
  if ((err = INIT(pts, 1, S)) != NOERROR) {
    return err; 
  }
  //   printf("AC moments=%d\n", cov->mpp.moments);

  if (!grid) { // todo
    SERR("non-grid not programmed yet"); // meiste da; fehlt noch
    // welche Vektoren ich wirklich brauche
  }

  pgs->size =  grid ? intpow(2, dim) : loc->totalpoints;
  if (cov->role == ROLE_POISSON_GAUSS) {
    if (
	(pgs->xgr[0] = (double*) CALLOC(dim * 3, sizeof(double))) == NULL ||
	(pgs->pos = (int *) CALLOC(dim, sizeof(int))) ==NULL ||
	(pgs->min = (int*) CALLOC(dim, sizeof(int))) == NULL ||
	(pgs->max = (int*) CALLOC(dim, sizeof(int))) == NULL
	) return ERRORMEMORYALLOCATION;
    for (i=3, d=1; d<dim; d++, i+=3) pgs->xgr[d] = pgs->xgr[0] + i;
    // printf("%d\n", pgs->xgr[1] - pgs->xgr[0]); assert(false);
   if ((err = calculate_mass_gauss(cov)) != NOERROR) return err;
  } else if (hasMaxStableRole(cov)) {
    if (
	(pgs->single = (double*) CALLOC(pgs->size, sizeof(double)))==NULL ||
	(pgs->total = (double*) CALLOC(pgs->size, sizeof(double))) ==NULL ||
	(pgs->halfstepvector = (double*) CALLOC(dim, sizeof(double)))==NULL
	)
      return ERRORMEMORYALLOCATION;

    // APMI(cov);

    if ((err = calculate_mass_maxstable(cov)) != NOERROR) return err;
    cov->mpp.log_zhou_c = log(pgs->totalmass); 

    //  printf("zhou %f %f\n",  cov->mpp.log_zhou_c, pgs->totalmass); APMI(cov); assert(false); // x,x,x:zhou 2.161274 8.682194; x,x,0: zhou 1.177634 3.246682


    if (shape->deterministic) {
      cov->mpp.maxheight = pts->mpp.maxheight * shape->mpp.maxheight;
      if (!R_FINITE(cov->mpp.maxheight)) BUG;
    } else { // currently both cases are identical ...
      cov->mpp.maxheight = pts->mpp.maxheight * shape->mpp.maxheight;
      if (!R_FINITE(cov->mpp.maxheight)) {
	//	APMI(cov);
	BUG;
      }
    }



  } else {
    // APMI(cov);
    BUG;
  }


  if (CovList[shape->nr].nonstat_inverse == ErrInverseNonstat) {
    if (pts->nr != RECTANGULAR) {
      // APMI(shape);
      warning("Inverse of shape function cannot be determined. Simulation speed  might be heavily decreased.");
    }
  }
  
 
  for (i=0; i<=cov->mpp.moments; i++) {
    //    printf("%d %f %f %d\n", i, shape->mpp.M[i], shape->mpp.Mplus[i], cov->mpp.moments);

    cov->mpp.M[i] = pts->mpp.M[i];
    cov->mpp.Mplus[i] = pts->mpp.Mplus[i];
  }
   
  cov->rf = shape->rf;
  cov->origrf = false;
 

  //APMI(cov);

  // assert(false);

  //  printf("init done %d\n", err);

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
void do_pgs_gauss(cov_model *cov, storage *S) {
  pgs_storage *pgs = cov->Spgs;
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  location_type *loc = Loc(cov);
  //  window_info *w = &(S->window);
  int i, d, 
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
    int err;
    DO(shape, S); 
    // APMI(pts);
    if ((err = INIT(pts, 1, S)) != NOERROR) {
      XERR(err); 
    }
    DORANDOM(pts, cov->q);  // cov->q nur dummy. Wird ueberschrieben
    if (cov->role == ROLE_POISSON_GAUSS || !grid) { 
      if (calculate_mass_gauss(cov) != NOERROR) 
	error("unexpected error in 'do_pts_given_shape' (maxstable)");
    } else if (cov->role == ROLE_MAXSTABLE) { 
      BUG;
      if (calculate_mass_maxstable(cov) != NOERROR)
	error("Unexpected error in 'do_pts_given_shape' (maxstable)");  ;
      cov->mpp.log_zhou_c = log(pgs->totalmass);
      cov->mpp.maxheight = pts->mpp.maxheight * shape->mpp.maxheight;
      if (!R_FINITE(cov->mpp.maxheight)) BUG;
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
    if (ISNA(x[0]) || x[0] > y[0]) BUG;

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
    for (i=0; i<loc->totalpoints; i++, xx+=dim) {
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


//static int counter[2]= {0,0};

void do_pgs_maxstable(cov_model *cov, storage *S) {
  pgs_storage *pgs = cov->Spgs;
  cov_model *shape = cov->sub[PGS_FCT],
    *pts = cov->sub[PGS_LOC];
  location_type *loc = Loc(cov);
  //window_info *w = &(S->window);
  int i, d, elmts,
    dim = shape->xdimprev;
  double 
    *single = pgs->single, // in
    *total = pgs->total,   // in
    *halfstepvector = pgs->halfstepvector, // in
    *x = pgs->x, // dummy
    *v = pgs->v; // dummy
  bool
    flat = pgs->flat;


  if (!cov->deterministic) {
    int err;
    DO(shape, S); 
    if ((err = INIT(pts, 1, S)) != NOERROR) {
      XERR(err); 
    }
    DORANDOM(pts, cov->q);  // cov->q nur dummy. Wird ueberschrieben
     // CovList[shape->nr].Do(shape, S); // nicht gatternr
    if (calculate_mass_maxstable(cov) != NOERROR) 
      error("unexpected error in 'do_pts_given_shape' (maxstable)");
  }

  i = DrawCathegory(pgs->size, single, total, loc->grid, &elmts);

  //   pgs:single 0.000000 1.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.020262 
  //     pgs:total  0.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 21.261911 
   //  printf("i=%d total=%f grid=%d elm=%d sin=%f, size=%d\n", i, *total, loc->grid, elmts, *single, pgs->size);
 
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

  if (flat) { 
    if (i != PGS_VOXEL) VTLG_R(x, pts, v); 
    //  printf("v=%f x=%f %d %d %d \n", v[0],  x[0], i, PGS_VOXEL,  PGS_CORNER);
    for (d = 0; d<dim; d++)
      if (x[d] == 0.0) 
	v[d] = loc->xgr[d][XSTEP] * UNIFORM_RANDOM - halfstepvector[d];
  } else { 
    //  printf("here\n");
    VTLG_R2SIDED(NULL, x, pts, v); 
  }

  //printf("i=%d %f %f %f v=%f %f %f\n", i, x[0], x[1], x[2], v[0], v[1], v[2]);

  for (d=0; d<dim; d++) {
    cov->q[d] = loc->xgr[d][XSTART] + v[d];
    // printf("q=%f %f %f %d\n", cov->q[d], loc->xgr[d][XSTART], v[d], flat);
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


  //   if (counter[0]==0 && counter[1]==0) PMI(cov);
 
  if (flat) for (d = 0; d<dim; d++) if (x[d] == 0.0)  v[d] = 0.0;
 
  if (CovList[pts->nr].logD != ErrLogCov) {
    // printf("log!\n");
    double sign;
    VTLG_DLOG(v, pts, &(pgs->log_density), &sign); // hier muss allgemein die nicht-normierte Dichte stehen und stattdessen durch das Mittel des Volumens geteilt  werden
   

    
    //      printf("log density v=%f, %f, %f %e\n", v[0],v[1],v[2], pgs->log_density);

     //printf("log_dens %f\n", pgs->log_density);
    //  pgs->log_density -= log(pgs->totalmass); // ueberprueft: normierung
    // zu einer WK'Dichte.
  } else {
    double density;
    VTLG_D(v, pts, &density);

    //    printf("density %e v=%f %f %f\n", density, v[0], v[1], v[2]);

    pgs->log_density = log(density);
  }

  //counter[cov->q[1] > 0.5]++; 
 //printf("counter %d %d cat =%d v>0:%d q=(%f %f %f)\n", 
 //	counter[0], counter[1], i, v[1]>0, cov->q[0], cov->q[1], cov->q[2]
 //	);
 //assert(counter[0] < 40);

  //  PMI(pts);
  
  if (!R_FINITE(pgs->log_density)) {
    //PMI(cov);
    BUG;
  }
  
}




void do_pts_given_shape(cov_model *cov, storage *S) {
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
    do_pgs_gauss(cov, S);
    eps = GLOBAL.mpp.about_zero * exp(pgs->log_density);
  } else if (hasMaxStableRole(cov)) { // todo: sauber trennen, wann 
    // max-stable, smith etc
    do_pgs_maxstable(cov, S); 
    eps = pgs->currentthreshold;
    if (cov->loggiven)  eps += pgs->log_density;
    else eps *= exp(pgs->log_density);
  } else {
    PMI(cov); BUG; //
  }

  //  APMI(shape);
  
  NONSTATINVERSE(&eps, shape, x, y);
  if (ISNA(x[0]) || x[0] > y[0]) {
    assert(!cov->loggiven);
    //double eps_neu = eps / cov->mpp.maxheight; // warum ?? 29.12.2013
    double eps_neu = eps; //  29.12.2013
    NONSTATINVERSE_D(&eps_neu, pts, x, y); 
    if (ISNA(x[0]) || x[0] > y[0]) BUG;
 
      //printf("eps=%f %f %f x=(%e %e %e) y=(%e %e %e) q=(%f %f %f)\n",
      //     eps, pgs->currentthreshold,  exp(pgs->log_density),
      //     x[0], x[1], x[2], y[0], y[1], y[2],
      //     cov->q[0], cov->q[1], cov->q[2]);
  }
 
  //  assert(!R_FINITE(y[0]) || y[0] < 55.281);

  //APMI(cov);
  
  for (d=0; d<dim; d++) {

    //     printf("do d=%d q=%f x=%f y=%f eps=%e thr=%e\n",
    //	  d, cov->q[d], x[d], y[d], eps, pgs->currentthreshold);

    pgs->supportmin[d] = cov->q[d] - 10 * y[d]; // 4 * for debugging...
    pgs->supportmax[d] = cov->q[d] - 10 * x[d];

    if (ISNA(pgs->supportmin[d]) || ISNA(pgs->supportmax[d]) ||
	pgs->supportmin[d] > pgs->supportmax[d]) {
      //printf("do d=%d q=%f x=%f y=%f eps=%e thr=%e supp=(%e, %e)\n",    d, cov->q[d], x[d], y[d], eps, pgs->currentthreshold,  pgs->supportmin[d], pgs->supportmax[d]);

      BUG;
    }
  
    //assert(pgs->supportmin[d] <= pgs->supportmax[d]);
  }

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

  if (newmodel != NULL) BUG;


  if (shape->role != ROLE_POISSON && shape->role != ROLE_MAXSTABLE)
    ILLEGAL_ROLE;

  //   printf("here %s\n", NICK(cov->sub[PGS_LOC]));
  
  //  APMI(cov);
       
  //if ((err = CHECK_R(cov->sub[PGS_LOC], cov->tsdim)) != NOERROR) return err;
  
  return NOERROR;
}
 

int init_standard_shape(cov_model *cov, storage *S) {  
  cov_model *shape = cov->sub[PGS_FCT];
  //  location_type *loc = Loc(cov);

  //APMI(cov);

  assert(cov->sub[PGS_LOC] != NULL);
  location_type *loc = Loc(cov);
  int d,
    dim = shape->xdimprev,
    err = NOERROR;
 
  //APMI(cov);

  assert(shape->xdimprev == Loc(cov)->timespacedim);

  if ((err = alloc_pgs(cov)) != NOERROR) return err;
  pgs_storage *pgs = cov->Spgs;
  if ((pgs->localmin = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->localmax = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->minmean = (double*) CALLOC(dim, sizeof(double))) == NULL ||
      (pgs->maxmean = (double*) CALLOC(dim, sizeof(double))) == NULL)
    return ERRORMEMORYALLOCATION;

 
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
  if (ISNA(x[0]) || x[0] > y[0])
    SERR1("inverse of '%s' unknown", NICK(shape));
  GetDiameter(loc, pgs->localmin, pgs->localmax, 
	      pgs->supportcentre); // last is only a dummy
  pgs->totalmass = 1.0;
  for (d=0; d<dim; d++) {
    min[d] = pgs->localmin[d] - y[d];
    max[d] = pgs->localmax[d] - x[d];
    if (!(R_FINITE(min[d]) && R_FINITE(max[d])))
      SERR("simulation window does not have compact support. Should 'RMtruncsupport' be used?");
    pgs->totalmass *= max[d] - min[d];
  }
  
  if (cov->role == ROLE_POISSON) {
    pgs->log_density = 0.0;
  } else {    
    pgs->log_density = 0.0;
    cov->mpp.log_zhou_c = log(pgs->totalmass); // todo: logzhou gehoert in pgs
    cov->mpp.maxheight = shape->mpp.maxheight;
  }

  cov->rf = shape->rf;
  cov->origrf = false;
  cov->fieldreturn = shape->fieldreturn;

  // APMI(cov);

  return NOERROR;
}


void do_standard_shape(cov_model *cov, storage *S) {
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
  if (ISNA(x[0])|| x[0] > y[0]) BUG;

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

  if (newmodel != NULL) BUG;


  if (shape->role != ROLE_POISSON && shape->role != ROLE_MAXSTABLE)
    ILLEGAL_ROLE;

  //printf("here\n");

  return NOERROR;
}


int init_stationary_shape(cov_model *cov, storage *S) {  
  cov_model *shape = cov->sub[STAT_SHAPE_FCT];
   int d,
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
 
  cov->mpp.log_zhou_c = 0.0;
  if (!R_FINITE(cov->mpp.log_zhou_c))
    SERR1("max height of '%s' not finite", NICK(shape));
  
  pgs->log_density = 0; 
   
  //  printf("dims %d %d\n", cov->xdimown, shape->xdimprev);

  for (d=0; d<dim; d++) {
    pgs->supportmin[d] = RF_NEGINF; // 4 * for debugging...
    pgs->supportmax[d] = RF_INF;
  }

  cov->mpp.maxheight = shape->mpp.maxheight;
  cov->rf = shape->rf;
  cov->origrf = false;
  cov->fieldreturn = shape->fieldreturn;
  if (!cov->fieldreturn) BUG;

  return NOERROR;
}


void do_stationary_shape(cov_model *cov, storage *S) {
  cov_model *shape = cov->sub[STAT_SHAPE_FCT]; 
  DO(shape, S);
  cov->mpp.maxheight = shape->mpp.maxheight;
  assert(shape->fieldreturn);
}


