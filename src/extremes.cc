 /* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 simulation of max-stable random fields

 Copyright (C) 2001 -- 2013 Martin Schlather, 

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


double GetTotal(cov_model* cov) {
  cov_fct *C = CovList + cov->nr;
  double v = 1.0;
  int i,
    nsub = cov->nsub,
    kappas = C->kappas;							

  if (C->Type == RandomType) {
    if (cov->total_n >= 0) {
      int VARIABLE_IS_NOT_USED old_n = cov->total_n;
      if (cov->total_n < GLOBAL.distr.repetitions) {
	double *r = (double*) MALLOC(cov->xdimown * sizeof(double));
	while(cov->total_n < GLOBAL.distr.repetitions) {
	  DORANDOM(cov, r);
	  assert(cov->total_n > cov->total_n);
	}
	free(r);
      }
      v *= cov->total_sum / (double) cov->total_n;
    }
  } else {
    for (i=0; i<nsub; i++) {
      v *= GetTotal(cov->sub[i]);
    }    
    for (i=0; i<kappas; i++) {						
      if (cov->kappasub[i] != NULL) 
	v *= GetTotal(cov->kappasub[i]);  
    }
  }
  return v;
}


int addStandard(cov_model **Cov) {
  cov_model *cov = (*Cov)->calling,
    *next = *Cov;
  int i, err,
    dim = next->xdimprev,
    vdim = next->vdim,
    role = next->role;     

  addModel(Cov, STANDARD_SHAPE);
  cov_model *key = *Cov;
  SetLoc2NewLoc(key, Loc(cov));
  assert(key->calling == cov);
  assert(key->sub[PGS_LOC] == NULL && key->sub[PGS_FCT] != NULL);

  for (i=0; i<2; i++) {
    if ((err = CHECK(key, dim, dim, PointShapeType, XONLY, NO_ROTAT_INV,
		     vdim, role)) != NOERROR) return err;   
    if (i==0) {
      if (hasPoissonRole(cov)) {
	addModel(key->sub + PGS_LOC, UNIF);
	key->sub[PGS_LOC]->calling = cov;
      } else {
	if ((err = STRUCT(key, key->sub + PGS_LOC)) != NOERROR) return err;
      }

    }

  }
  return NOERROR;
}


int addPGS(cov_model **Cov) {
  // for m3 & random coin
  cov_model *shape = *Cov;
  assert(shape->calling != NULL);
  //domain_type dom = shape->domprev;
  //  isotropy_type iso = shape->isoprev;
  int err,
    dim = shape->xdimprev,
    vdim = shape->vdim,
    role = shape->role;


  //PMI(shape);
  
  assert(dim == shape->tsdim);
  assert(vdim == 1);
    

  // most models split into a shape function and location distribution
  // given the shape function, see oesting, schlather, chen
  //
  // for anisotropic models also the reverse modelling could be of interest:
  // first the location_type aere modelled. Then conditional on the locations
  // the shape functions are modelled. Um diese Option zu ermoeglichen muesste
  // noch ein Schalter programmiert werden, oder es ist implizit dadurch
  // gegeben, dass das Modell eben als ganzen gegeben ist.
  //
  // Or the model is given as a whole. This is the first if-condition:
  //      if (isPointShape(shape)) return NOERROR;
  // else: model is given by the shape and then conditional on the 
  // location

  
  addModel(Cov, PTS_GIVEN_SHAPE);
  assert((*Cov)->sub[PGS_LOC] == NULL && (*Cov)->sub[PGS_FCT] != NULL);
  if ((err = CHECK(*Cov, dim, dim, PointShapeType, XONLY, NO_ROTAT_INV,
		     vdim, role)) != NOERROR) return err; 
  if ((err = STRUCT((*Cov), (*Cov)->sub + PGS_FCT)) != NOERROR) return err;
  if ((err = CHECK(*Cov, dim, dim, PointShapeType, XONLY, NO_ROTAT_INV, 
		     vdim, role)) != NOERROR) return err;

  //   APMI(*Cov);

  return NOERROR;
}


  void range_mpp(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
    //   assert(cov != NULL);
  range->min[GEV_XI] = RF_NEGINF;
  range->max[GEV_XI] = RF_INF;
  range->pmin[GEV_XI] = -10;
  range->pmax[GEV_XI] = 10;
  range->openmin[GEV_XI] = true;
  range->openmax[GEV_XI] = true; 

  range->min[GEV_MU] = RF_NEGINF;
  range->max[GEV_MU] = RF_INF;
  range->pmin[GEV_MU] = -1000;
  range->pmax[GEV_MU] = 1000;
  range->openmin[GEV_MU] = true;
  range->openmax[GEV_MU] = true; 

  range->min[GEV_S] = 0;
  range->max[GEV_S] = RF_INF;
  range->pmin[GEV_S] = 0.001;
  range->pmax[GEV_S] = 1000;
  range->openmin[GEV_S] = true;
  range->openmax[GEV_S] = true; 
}


int init_mpp(cov_model *cov, storage *S) {
  cov_model *sub = cov->key != NULL ? cov->key :
    cov->sub[0] != NULL ? cov->sub[0] : cov->sub[1];
  if (sub == NULL) SERR("substructure could be detected");
  location_type *loc = Loc(cov);
  //window_info *w = &(S->window);
  int d, err,
    dim = cov->tsdim,
    role = cov->role, 
    maxstable = hasMaxStableRole(cov);
  pgs_storage *pgs;


  if (!maxstable && !hasPoissonRole(cov)) ILLEGAL_ROLE;
  if (!isPointShape(sub)) SERR1("%s is not shape/point function", NICK(sub));
  if (loc->distances) return ERRORFAILED;
  
  //   printf("%d \n", 1 + (role == ROLE_POISSON_GAUSS));APMI(sub->calling);
  //  print("nr =%d gatter=%d \n", sub->nr, sub->gatternr);
  //  assert(sub->gatternr >= 0);
  //  print("Init=%ld %ld\n", CovList[sub->gatternr].Init, init2);
  //  print("S=%ld\n", S);
  //
  // print("end list\n");
 
  if ((err = INIT(sub, maxstable ? 1 : role == ROLE_POISSON ? 0 : 2, S)) 
      != NOERROR)  return err;
  pgs = sub->Spgs;

  //assert(false);

  //PMI(cov);


  assert(pgs != NULL);
   
  for (d=0; d<dim; d++)
    pgs->supportmin[d] = pgs->supportmax[d] = pgs->supportcentre[d] = RF_NAN;
  GetDiameter(loc, pgs->supportmin, pgs->supportmax, pgs->supportcentre); 
  
  if (maxstable) {      
    if (ISNA(sub->mpp.Mplus[1]) || !R_FINITE(sub->mpp.Mplus[1])
	|| sub->mpp.Mplus[1] <= 0.0 ) {
      // A
      //      PMI(sub);
      //      crash();
      SERR("integral of positive part of submodel unkown");
    }

    //PMI(sub, "sub");

    if (!R_FINITE(sub->mpp.log_zhou_c) )
      SERR("maximal height of submodel unkown");
  } else if (role == ROLE_POISSON_GAUSS) {
    if (ISNA(sub->mpp.M[2]) || !R_FINITE(sub->mpp.M[2] || sub->mpp.M[2] <=0.0)){
      SERR("second moment unkown");
    }
  } // else role == ROLE_POISSON -- nothing to do
 
  if ((err = FieldReturn(cov)) != NOERROR) { // must be later than INIT !
    return err;  
  }
  
  cov->simu.active = err == NOERROR;

  return err;
}



void dompp(cov_model *cov, storage *s) {
  assert(cov->simu.active);
  location_type *loc = Loc(cov);	
  if (loc->caniso != NULL) BUG;

  // window_info *w = &(s->window);
  cov_model  
    *sub = (cov->key != NULL) ? cov->key : 
    cov->sub[0] != NULL ? cov->sub[0] : cov->sub[1];
  if (sub == NULL) ERR("substructure could be detected");
  /*  mpp_model randomcoin = C->randomcoin; */	
  //
  pgs_storage *pgs = sub->Spgs;

  assert(pgs != NULL);
 
  if (!isPointShape(sub)) BUG;

  //  PMI(cov);
  
  int i, k, d, 
    *gridlen = pgs->gridlen,
    *end = pgs->end,
    *start = pgs->start,
    *delta = pgs->delta,
    *nx = pgs->nx,
    dim = cov->tsdim,
    spatial = loc->totalpoints,
    every = GLOBAL.general.every, //
    nthreshold = (every>0) ? every : MAXINT,	
    //    covrole = cov->role,
    subrole = sub->role,
    deltathresh = nthreshold;			

  double logdens,  factor, summand, 
    sign[2], value[2], logthreshold, 
    logM1 = RF_NAN,
    logM2 = RF_NAN,
    gumbel = RF_NAN,
    *xstart = pgs->xstart, // nach DO gesetzt
    *x = pgs->x,           // nach DO gesetzt
    *inc = pgs->inc,
    threshold = RF_NAN,
    *rf=NULL, 
    poisson=0.0;
  res_type
    *res = cov->rf;
  bool fieldreturn,
    maxstable = hasMaxStableRole(cov),
    compoundpoisson = subrole == ROLE_POISSON,
    poissongauss = subrole == ROLE_POISSON_GAUSS,
    simugrid = loc->grid,
    // AVERAGE braucht Sonderbehandlung:
     randomcoin = cov->method == RandomCoin //only for ROLE_POISSON_GAUSS
    ;
  long zaehler, n, cumgridlen[MAXMPPDIM +1],
    Total_n = maxstable ? -1 : rpois(pgs->intensity),			
    total_pts = loc->totalpoints, 
    vdimtot = total_pts * cov->vdim,
    control = 0;

  assert(maxstable || Total_n > 0);
    
  cov_model *key = cov->key;
  ext_bool loggiven = key == NULL ? cov->sub[0]->loggiven : key->loggiven;
  
  //  assert(covrole != ROLE_POISSON_GAUSS || subrole == ROLE_POISSON_GAUSS);
  //  assert(loc->caniso == NULL);
  // assert(covrole != ROLE_POISSON_GAUSS || sub->mpp.moments >= 2);
  // assert(covrole == ROLE_POISSON_GAUSS || covrole == ROLE_POISSON || 
  //	 sub->mpp.moments >= 1);
  
  if (simugrid) {
    cumgridlen[0] = 1;
    for (d=0; d<dim; d++) {
      inc[d] = loc->xgr[d][XSTEP];
      gridlen[d] =  loc->length[d]; 
      cumgridlen[d+1] = loc->length[d] * cumgridlen[d];

      //printf("set d=%d inc=%f gridlen=%d cum=%d\n", d, inc[d], gridlen[d],
      //     (int) cumgridlen[d]);
    }
  }

  if (maxstable) {
    if (loggiven) {
      for (i=0; i<vdimtot; i++) res[i] = R_NegInf;
    }
    else for (i=0; i<vdimtot; i++) res[i] = 0.0;
    if (sub->mpp.moments < 1 || !R_FINITE(sub->mpp.Mplus[1]) || 
      sub->mpp.Mplus[1] == 0.0)
      ERR("integral of positive part of the submodel is unknown or not finite");
    logM1 = log(sub->mpp.Mplus[1]);
    threshold = logthreshold = RF_INF;
    pgs->globalmin = loggiven ? RF_NEGINF : 0;
    pgs->currentthreshold =
      loggiven ? pgs->globalmin - threshold : pgs->globalmin / threshold;
  } else if (compoundpoisson){
    for (i=0; i<vdimtot; i++) res[i] = 0.0;
  } else {
    for (i=0; i<vdimtot; i++) res[i] = 0.0;
    logM2 = log(sub->mpp.M[2]);
    pgs->currentthreshold = 1e-8;
   }
 
 
  for (n=0; ; n++) {    
    //     printf("n=%d tot=%d contr=%d: %f < %f; res[0]=%f\n", (int) n, (int) Total_n, (int) control, res[control], threshold, res[0]);// assert(n <= 10);


    //    assert(res[0] < 1e-4);
    //    assert(res[0] < 80);

    if (!maxstable && (n >= Total_n)) break; 
    
    // for (d=0; d<dim; d++) info->min[d] = info->max[d] = 0.0;

    // printf(".\n");

    DO(sub, s);   
    //    printf("here!!!\n");

    //PMI(sub);

    fieldreturn = sub->fieldreturn;
    if (subrole == ROLE_BROWNRESNICK) {
      assert(sub->SBR != NULL);
      n += sub->SBR->hatnumber;
    }
    
    logdens = pgs->log_density;

    //    PMI(sub, -1);

    if (!R_FINITE(logdens)) {
      // PMI(sub, -1);
      BUG;
    //    get_shape = CovList[sub->gatternr].cov;
    //   get_logshape = CovList[sub->gatternr].log;
    }
   
    if (compoundpoisson) {
      summand = 0.0;
    } else if (poissongauss) {
      summand = -0.5 * (logdens + logM2);
    } else { // max-stable field    
      //PMI(cov);
      
      assert(subrole == ROLE_SMITH || subrole == ROLE_SCHLATHER ||
	     subrole == ROLE_BROWNRESNICK);

      if (n >= GLOBAL.extreme.maxpoints) {
	PRINTF("'maxpoints' (= %d) exceeded. Simulation is likely to be a rough approximation only. Consider increasing 'maxpoints'.\n", 
	       GLOBAL.extreme.maxpoints
		);
	//	printf("break1 %d %d\n",  n, ((int*) cov->p[MAXSTABLE_POINTS])[0]);
	break;
      }
      poisson += rexp(1.0); 
      gumbel = -log(poisson);
      // MARCO, habe ich am 17.11.13 geaendert
      summand = gumbel + sub->mpp.log_zhou_c - logdens - logM1;
      
      //   threshold = sub->loggiven ? GumbelThreshold : FrechetThreshold;
      threshold = logthreshold =
	(res_type) (gumbel + sub->mpp.log_zhou_c + sub->mpp.maxheight - logM1); //Gumbel
      if (!loggiven) {
	threshold = exp(threshold); // Frechet
      }

 
      // printf("for thres=%e %d res=%e log.zhou=%f(%f) logden=%f logM1=%e gumbel=%f loggiven=%d summand=%f\n", threshold, (int) control, res[control], sub->mpp.log_zhou_c, sub->mpp.maxheight, logdens, logM1, gumbel, loggiven, summand);

      //   PMI(sub);

      assert(R_FINITE(threshold ));

// {   int i; for (i=0; i<total_pts; i++) printf("%e ", res[i]); }

      //APMI(cov);
      while ((control<total_pts) && (res[control]>=threshold)) {
	//   print("control: %d %f %f\n", control, res[control], threshold);
	control++;
      }
      if (control >= total_pts) {
	//printf("break2\n");
	break;
      }
    
      if (n % GLOBAL.extreme.check_every == 0) {
	pgs->globalmin = res[control];
	for (i=control + 1; i<total_pts; i++) 
	  if (res[i] < pgs->globalmin) pgs->globalmin = res[i];
      }
 
       
      pgs->currentthreshold =
	loggiven ? pgs->globalmin - gumbel : pgs->globalmin * exp(-gumbel);
     

      //printf("loggiven =%d %e %f cur.thr=%e\n",
      //     loggiven, pgs->globalmin, gumbel,  pgs->currentthreshold);

    }
    factor = exp(summand);

  //
  //printf("factor %4.4f sumd=%4.4f logdens=%4.4f logM2=%4.4f th=%4.4f %d<%d \n", factor, summand, logdens, logM2, threshold, control, total_pts);  //  assert(false);
    
  //   printf("dim=%d loggiven=%d maxstable=%d simugrid=%d randomcoin=%d\n",
  //	 dim, loggiven, maxstable, simugrid, randomcoin);// assert(false);

    //printf("A=%d\n", n);
    zaehler = 0;
    if (simugrid) {

      for (d=0; d<dim; d++) {	
	double dummy = ceil((pgs->supportmin[d] - loc->xgr[d][XSTART]) / 
			    inc[d]);
	start[d] = (dummy < 0) ? 0 : (int) dummy;
	dummy = 1.0 + ((pgs->supportmax[d] - loc->xgr[d][XSTART]) / inc[d]);
	end[d] = (dummy >= gridlen[d]) ? gridlen[d] : (int) dummy;

	//printf("window [%f %f] [%d %d]\n", pgs->supportmin[d], pgs->supportmax[d], start[d], end[d]);


	if (start[d] >= end[d]) { // 
	  //PMI(cov);
	  //  printf("broken %d %d %d supp=%f %f, loc.start=%f %f\n", d,
	  //	 start[d], end[d],
	  //	 pgs->supportmin[d], pgs->supportmax[d],
	  //	 loc->xgr[d][XSTART], inc[d]); 
	  break;
	}
	delta[d] = (end[d] - start[d]) * cumgridlen[d];
	nx[d] = start[d];
	zaehler += start[d] * cumgridlen[d];
	x[d] = xstart[d] = loc->xgr[d][XSTART] + (double) start[d] * inc[d]; 
	
	if (false) 
	   printf("d=%d start=%d, end=%d xstart=%f %f pgs=[%4.3f, %4.3f] xgr=%f %f %f inc=%3.2f\n", //
	  	 d, start[d], end[d], xstart[d], pgs->xstart[d], pgs->supportmin[d], pgs->supportmax[d],
	 	 loc->xgr[d][XSTART], loc->xgr[d][XSTEP], loc->xgr[d][XLENGTH],
	  	 inc[d]);
	// assert(false);
	
	}

      //     APMI(cov);
      if (d < dim) continue;
    }

    //  printf("simugrid %d\n", simugrid);
    //  assert(false);
    //    APMI(cov);
 

#define STANDARDINKREMENT_ZAEHLER		\
  d = 0;			\
  nx[d]++;			\
  x[d] += inc[d];		\
  zaehler += cumgridlen[d];	\
  while (nx[d] >= end[d]) {	\
    nx[d] = start[d];		\
    x[d] = xstart[d];		\
    zaehler -= delta[d]; /* delta is positive */	\
    if (++d >= dim) break;	\
    nx[d]++;			\
    x[d] += inc[d];		\
    zaehler += cumgridlen[d];					\
  }								\
  /* printf("nx=%d %d %d\n", nx[0], nx[1], zaehler); // */	\
  if (d >= dim) break;		
    // end define StandardInkrement

 
#define GET SHAPE(x, sub, value); /*printf("%f : %f fctr=%f %ld \n", x[0], *value, factor, zaehler); // */
#define LOGGET LOGSHAPE(x, sub, value, sign);
#define LOGGETPOS LOGGET if (sign[0] > 0)

#define RFMAX(OP) {				\
      rf = sub->rf;				\
      for (i=0; i<total_pts; i++) {		\
	*value = (res_type) (OP);		\
	if (res[i] < *value) res[i]=*value;	\
      }						\
    }

#define RFSUM(OP) {				\
      rf = sub->rf;				\
      for (i=0; i<total_pts; i++) {		\
	*value = (res_type) (OP);		\
	if (res[i] < *value) res[i] += OP;	\
      }						\
    }

    //  gridlcoations/not;
#define MAXAPPLY(G,OP)	{					\
      if (simugrid) {						\
	while (true) {						\
	  G {							\
	    OP;							\
	    if (res[zaehler] < *value) res[zaehler]=*value;	\
	    STANDARDINKREMENT_ZAEHLER;				\
	  }							\
	}							\
      } else {							\
	x=loc->x;						\
	for(; zaehler<spatial; zaehler++, x += dim) {		\
	  G {							\
	    OP;							\
	    if (res[zaehler] < *value) res[zaehler] = *value;	\
	  }							\
	}							\
      }								\
    }
 
    // printf("%d %f %f res=%f %ld\n", n, *value, OP, res[zaehler], zaehler);
#define SUMAPPLY(G,OP) {						\
      if (simugrid) {							\
	while (true) {G; res[zaehler]+=OP; STANDARDINKREMENT_ZAEHLER;}	\
      } else {								\
	x=loc->x;						\
	for(; zaehler<spatial; zaehler++, x += dim) {			\
	  G; res[zaehler] += OP;					\
	}								\
      }									\
    }

#define AVEAPPLY(G,OP0,OP1) {						\
      warning(" * timescale ?!");					\
      if (simugrid) {							\
	while (true) {							\
	  int zt;  double zw,inct, segt, ct, st, A, cit, sit;		\
	  inct = sub->q[AVERAGE_YFREQ] * loc->xgr[dim][XSTEP];		\
	  cit = cos(inct);						\
	  sit = sin(inct);						\
	  G;								\
	  segt = OP1;							\
	  A = OP0;							\
	  ct = A * cos(segt);						\
	  st = A * sin(segt);						\
	  for (zt = zaehler; zt < total_pts; zt += spatial) {		\
	    res[zt] += (res_type) ct;					\
	    zw = ct * cit - st * sit;					\
	    st = ct * sit + st * cit;					\
	    ct = zw;							\
	    res[zaehler] += zw;						\
	  }								\
	  STANDARDINKREMENT_ZAEHLER;					\
	}								\
      } else {  /* not a grid */					\
        x=loc->x;							\
	/* the following algorithm can greatly be improved ! */		\
	/* but for ease, just the stupid algorithm	        */	\
	for (; zaehler<spatial; zaehler++, x += dim) {			\
	  G; res[zaehler] += OP0 * cos(OP1);				\
	}								\
      }									\
    }	

    // *atomdependent rfbased/not; (( recurrent/dissipative: min/max gesteuert))
    if (poissongauss  && !randomcoin) { // AVERAGE
      if (sub->loggiven) {
	AVEAPPLY(LOGGET, sign[0] * exp(value[0] + summand), 
		 sign[1] * exp(value[1]));
      } else AVEAPPLY(GET, value[0] * factor, value[1]);
    } else {
      assert(subrole == ROLE_SMITH || subrole == ROLE_SCHLATHER ||
	     subrole ==ROLE_BROWNRESNICK || subrole==ROLE_POISSON || 
	     subrole==ROLE_POISSON_GAUSS);

      if (fieldreturn) { // atomdependent rfbased/not;     
	// todo : !simugrid! && Bereiche einengen!!     
	// note ! no sign may be given currently when fieldreturn and loggiven!!
	if (sub->loggiven) {
	  if (maxstable) {RFMAX(rf[i] + summand);}
	  else RFSUM(exp(rf[i] + summand));	  
	} else {
	  if (maxstable) {RFMAX(rf[i] * factor);}
	  else RFSUM(rf[i] * factor);
	}
      } else { // not field return
	if (loggiven == true) {
	  if (maxstable) {MAXAPPLY(LOGGETPOS, (*value) += summand);}
	  else SUMAPPLY(LOGGET, sign[0] * exp(*value + summand));
	} else { // not loggiven
	  if (sub->loggiven){ 
	    if (maxstable) {
 	      MAXAPPLY(LOGGETPOS, *value=sign[0] * exp(*value+summand));
	    } else SUMAPPLY(LOGGET, sign[0] * exp(*value + summand));
	  } else {

 	    //
	    //
	    // printf("factor %f %f \n", factor, *value);
	    assert(R_FINITE(factor));
	    if (maxstable) { MAXAPPLY(GET, (*value) *= factor); }
	    else SUMAPPLY(GET, *value * factor);
	  }
	}
      }
    }
  
    // printf("B=%d\n", n);
    if (n >= nthreshold) {
      if (maxstable) {
	LPRINT("%d-th position: value=%f threshold=%f gu=%f %f %d (%d %d %d)\n",
	       control, (double) res[control], (double) threshold,
	       gumbel, logdens, loggiven,
	       n, nthreshold, deltathresh); 
   	nthreshold += deltathresh;
      } else {
	PRINTF("n=%d Total=%d\n", (int) (n / 1000), (int)(Total_n));
      }
    }

    R_CheckUserInterrupt();
    
    //printf("n=%d %f %f %f 0:%f\n", n, res[0], res[1], res[2], res[3]);

    //    assert(res[0] < 80);

  }// for(n,...    --     while control < total_pts

  //printf("n=%d %f %f %f 3::%f\n", n, res[0], res[1], res[2], res[3]);

  // Trafo to margins

  // printf("%d frechet %d \n", sub->loggiven, ev_p->Frechet);
 
  
  double finalc;
  finalc = GetTotal(sub);
  if (fabs(finalc - 1.0) > 1e-15) {
    assert(false);
    if (maxstable) {
      if (loggiven) {
	finalc = log(finalc);
	for (k=0; k<total_pts; k++) res[k] += finalc;
      } else for (k=0; k<total_pts; k++) res[k] *= finalc;
    } else { // gauss
      BUG; // to do
    }
  }
  

  if (maxstable){

    double xi = cov->p[GEV_XI][0];
    if (loggiven) { // 3.6.13 sub->loggiven geaendert zu loggiven
 
      //  printf("here %e %f\n", xi, res[0]);
      //      APMI(cov->calling);

      if (xi != 0.0) 
	for (k=0; k<total_pts; k++) res[k] = (exp(xi * res[k]) - 1.0) / xi;
    } else {
      if (xi==0.0) for (k=0; k<total_pts; k++) res[k] = log(res[k]);
      else if (xi == 1.0) for (k=0; k<total_pts; k++) res[k] -= 1.0;
      else for (k=0; k<total_pts; k++) res[k] = pow(res[k], xi) - 1.0;
    }
    
    // muss das vorletzte sein !
    if (cov->kappasub[GEV_S] != NULL) {
      error("'s' cannot be chosen as a function yet.");
      // check whether s is positive !!
    } else {
      double ss = xi == 0 ? cov->p[GEV_S][0] : (cov->p[GEV_S][0] / xi);
      if (ss != 1.0) for (k=0; k<total_pts; k++) res[k] *= ss;
    }  

    // muss das letzte sein !
    if (cov->kappasub[GEV_MU] != NULL) {
      error("mu cannot be chosen as a function yet.");
    } else {
      double mu = cov->p[GEV_MU][0];
      if (mu != 0.0) for (k=0; k<total_pts; k++) res[k] += mu;
    }
  } 
  //printf("n=%d %f %f %f 0:%f\n", n, res[0], res[1], res[2], res[3]);
  // for (k=0; k<total_pts; k++) printf("%f ", res[k]); printf("\n");

  return;
} // end dompp



//void calculate_integral(int d, double R, cov_model *cov,
//			integral_type *integral) {
//  error("calculate_integral not programmed yet.\n");
//}




////////////////////////////////////////////////////////////////////
// Poisson

int check_poisson(cov_model *cov) {
  cov_model *next=cov->sub[0],
    *key = cov->key,
    *sub = key == NULL ? next : key;
  int err,
    dim = cov->tsdim; // taken[MAX DIM],
  mpp_param *gp  = &(GLOBAL.mpp);
  Types type = sub == key ? PointShapeType : ShapeType; 

  //print("eee\n");
  
  cov->role = ROLE_POISSON;

  kdefault(cov, POISSON_INTENSITY, gp->intensity[dim]);
  if ((err = checkkappas(cov, true)) != NOERROR) return err;

  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown)
    return ERRORDIM;
  
  if ((err = CHECK(sub, dim, dim, type, XONLY, NO_ROTAT_INV,
		     SUBMODEL_DEP, cov->role)) != NOERROR) return err;
  setbackward(cov, sub);
 
  return NOERROR;
}


void range_poisson(cov_model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[POISSON_INTENSITY] = RF_NEGINF;
  range->max[POISSON_INTENSITY] = RF_INF;
  range->pmin[POISSON_INTENSITY] = -10;
  range->pmax[POISSON_INTENSITY] = 10;
  range->openmin[POISSON_INTENSITY] = true;
  range->openmax[POISSON_INTENSITY] = true; 
}


int struct_poisson(cov_model *cov, cov_model **newmodel){
  cov_model *next=cov->sub[0];
  location_type *loc = Loc(cov);
 
  if (newmodel != NULL) SERR("unexpected call of struct_poisson"); /// ?????
  if (cov->role != ROLE_POISSON)
    SERR1("'%s' not called as random coin", NICK(cov));  
  if (cov->key != NULL) COV_DELETE(&(cov->key));

  if (loc->Time || (loc->grid && loc->caniso != NULL)) {
    Transform2NoGrid(cov, false, GRIDEXPAND_AVOID);
  }

  if (!isPointShape(next)) {
    int err;
    if ((err = covcpy(&(cov->key), next)) != NOERROR) return err;
    if ((err = addStandard(&(cov->key))) != NOERROR) return err;
  }
 
  return NOERROR;
}


int init_poisson(cov_model *cov, storage *S) {
  //  location_type *loc = Loc(cov);
  // mpp_storage *s = &(S->mpp);
  cov_model *key=cov->key;
  int err;
 
  if ((err = init_mpp(cov, S)) != NOERROR) {
    return err;    
  }

  //APMI(cov);
  if (!isPointShape(key)) SERR("no definition of a shape function found");
  
  key->Spgs->intensity = key->Spgs->totalmass * cov->p[POISSON_INTENSITY][0];

  return err;
}




//////////////////////////////////////////////////////////////////////
//           Random Coin
//////////////////////////////////////////////////////////////////////


/*
  void coin(double *x, cov_model *cov, double *v){
  cov_model
  *shape = cov->sub[COIN_SHAPE],
  *next = cov->sub[COIN_COV];
  if (shape == NULL) COV(x, next, v);
  else {
  int i, d,
  vdim = cov->vdim;
  for (i=0; i<vdim; i++) v[i] = RF_NAN;
  for (d=0; d<cov->xdimown; d++) {
  if (x[d] != 0.0) {
  return;
  }
  }
  v[0] = 1.0; // only in the univariate case well defined
  }
  }


  void coinInverse(double *v, cov_model *cov, double *x){
  cov_model
  *shape = cov->sub[COIN_SHAPE],
  *next = cov->sub[COIN_COV];
  if (shape == NULL) INVERSE(v, next, x);
  else {
  int i, 
  vdim = cov->vdim;
  for (i=0; i<vdim; i++) x[i] = RF_NAN;
  }
  }
*/

int check_randomcoin(cov_model *cov) {
  cov_model 
    *next = cov->sub[COIN_COV],
    *shape = cov->sub[COIN_SHAPE],
    *key = cov->key;
  int err,
    dim = cov->tsdim; // taken[MAX DIM],
  mpp_param *gp  = &(GLOBAL.mpp);
  //extremes_param *ep = &(GLOBAL.extreme);

  //PMI(cov->calling, "check_randomcoin!");
  
  //PMI(cov, "check random coins");
  SERR("Currently not programmed yet.");
 
  ROLE_ASSERT(ROLE_POISSON_GAUSS || (cov->role==ROLE_GAUSS && cov->key!=NULL));
  ASSERT_ONE_SUBMODEL(cov);

  if (cov->sub[COIN_COV] != NULL && cov->sub[COIN_SHAPE]!=NULL)
    SERR("either 'tcf' or 'shape' must be given");

  if ((err = check_common_gauss(cov)) != NOERROR) return err;
  kdefault(cov, RANDOMCOIN_INTENSITY, gp->intensity[dim]);
  if ((err = checkkappas(cov, true)) != NOERROR) {
    //AERR(err);
    return err;
  }
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown)  {
    //AERR(1);
    return ERRORDIM;
  }


  if (key == NULL) {

    //printf("SHAPE %ld\n", shape);

    if (shape == NULL) {      
      if ((err = CHECK(next, dim, dim, PosDefType,
			 XONLY, SYMMETRIC,
			 SCALAR, ROLE_COV)) != NOERROR)  {
	return err;
      }      
      if ((next->pref[Average] == PREF_NONE) + 
	  (next->pref[RandomCoin] == PREF_NONE) !=1){

	//APMI(next);

	//assert(false);
	return ERRORPREFNONE;
      }
    } else { // shape != NULL
      if (next != NULL) SERR("phi and shape may not be given at the same time");
    
      if ((err = CHECK(shape, dim, dim, ShapeType, XONLY, NO_ROTAT_INV, 
			 SCALAR, ROLE_POISSON)) != NOERROR)  {
	return err;
      }     
      //   APMI(cov);
    }
    //print("ok\n");
  } else { // key != NULL
    // note: key calls first function to simulate the points
    // if this is true then AverageIntern/RandomCoinIntern calls Average/RandomCoin 
    cov_model *intern = cov;
    while (intern != NULL && intern->nr!=AVERAGE_INTERN  && isAnyDollar(intern))
      intern = intern->key != NULL ? intern->key : intern->sub[0];
   if (intern == NULL) {
      //APMI(key);
      BUG;
    } //else if (intern != cov) {
    //      printf("Here\n");
    //      APMI(cov);
    //      paramcpy(intern, cov, true, false);
    //    }
    
    //CHECK(key, dim, dim, PosDefType, XONLY, SYMMETRIC, 
    //		   SUBMODEL_DEP, ROLE_COV)
    
    if ((err = cov->role == ROLE_BASE || cov->role == ROLE_GAUSS
	 ? CHECK(key, dim, dim, ProcessType, XONLY, NO_ROTAT_INV, 
    		   SUBMODEL_DEP, ROLE_POISSON_GAUSS)
	 : cov->role == ROLE_POISSON_GAUSS 
	 ? CHECK(key, dim,  dim, PointShapeType, XONLY, NO_ROTAT_INV, 
		   SUBMODEL_DEP, ROLE_POISSON_GAUSS)
	 : ERRORFAILED
	 ) != NOERROR) {
      // APMI(cov)
      return err;
    }
    // APMI(cov)
  }
 
  cov_model *sub = key != NULL ? key : shape != NULL ? shape : next;
  setbackward(cov, sub);
 
  //  APMI(cov);

  return NOERROR;
}

void range_randomcoin(cov_model *cov, range_type *range) {
  range_common_gauss(cov, range);

  range->min[RANDOMCOIN_INTENSITY] = RF_NEGINF;
  range->max[RANDOMCOIN_INTENSITY] = RF_INF;
  range->pmin[RANDOMCOIN_INTENSITY] = -10;
  range->pmax[RANDOMCOIN_INTENSITY] = 10;
  range->openmin[RANDOMCOIN_INTENSITY] = true;
  range->openmax[RANDOMCOIN_INTENSITY] = true; 
}

 
int struct_randomcoin(cov_model *cov, cov_model **newmodel){
  cov_model 
    *next = cov->sub[COIN_COV],
    *shape = cov->sub[COIN_SHAPE];
  location_type *loc = Loc(cov);
  int err,
    dim = cov->tsdim; // taken[MAX DIM],

  // APMI(cov);

  ROLE_ASSERT(ROLE_POISSON_GAUSS);

  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if (loc->Time || (loc->grid && loc->caniso != NULL)) {
    Transform2NoGrid(cov, true, GRIDEXPAND_AVOID);
    SetLoc2NewLoc(next == NULL ? shape : next, Loc(cov)); 
  }

  if (newmodel != NULL) {
    //    printf("error in struct random coin %ld\n", (long int) (newmodel));
    //    PMI(cov->calling);
    SERR("unexpected call of stuct_randomcoin"); /// ?????
  }

  if (shape != NULL) {
    if ((err = covcpy(&(cov->key), shape)) > NOERROR) {
      return err;
    }
    if ((err = addPGS(&(cov->key))) != NOERROR) return err;
   
    // APMI(cov);
    return NOERROR;
  } else { // shape == NULL, i.e. covariance given
    if (next == NULL) BUG;
    // next not NULL
    if (next->pref[Average] == PREF_NONE && next->pref[RandomCoin]==PREF_NONE) {
      // if (next->nr > LASTDOLLAR) AERR(ERRORPREFNONE);
      return ERRORPREFNONE;
    }

    //    printf("%d \n", dim);
    //    PMI(next, "randomcoins");
    
    if ((err = CHECK(next, dim,  dim, PosDefType, XONLY, ISOTROPIC, 
		       SCALAR, ROLE_POISSON_GAUSS))
	!= NOERROR) {
      // APMI(cov)
      return err;
    }

    if ((err = STRUCT(next, &(cov->key))) > NOERROR) return err;
    if (cov->key == NULL)
      SERR("no structural information for random coins given");
    
    cov->key->calling = cov;
    if ( cov->pref[Average] == PREF_NONE ) {
      if (cov->key->nr != RANDOMSIGN) addModel(&(cov->key), RANDOMSIGN);
      assert(!isPointShape(cov->key));
      if ((err = addPGS(&(cov->key))) != NOERROR) return err;
    }
    
    //APMI(cov);
    return NOERROR;
  }
}


int init_randomcoin(cov_model *cov, storage *S) {
  cov_model
    *covshape = cov->sub[ cov->sub[COIN_SHAPE] != NULL ? COIN_SHAPE : COIN_COV],
    *key = cov->key,
    *sub = key == NULL ? covshape : key; 
  char name[] = "Poisson-Gauss";
  location_type *loc = Loc(cov);
  //mpp_storage *s = &(S->mpp);  
  int 
    err = NOERROR;

  sprintf(ERROR_LOC, "%s process: ", name);

  if (cov->role != ROLE_POISSON_GAUSS)  {
    // PMI(cov->calling, "init_randomcoin");
    ILLEGAL_ROLE;
  }

  cov->method = covshape->pref[Average] == PREF_NONE ? RandomCoin : Average;
    
  if (cov->method == Average &&  loc->caniso != NULL) {
    bool diag, quasidiag, semiseparatelast, separatelast;
    int idx[MAXMPPDIM];
    assert(loc->timespacedim <= MAXMPPDIM);
    analyse_matrix(loc->caniso, loc->cani_nrow, loc->cani_ncol, &diag,
		   &quasidiag, idx, &semiseparatelast, &separatelast);
    if (!separatelast) SERR("not a model where time is separated");
  }
 
  if ((err = init_mpp(cov, S)) != NOERROR) {
    //PMI(cov);    AERR(err);
    return err;    
  }

  sub->Spgs->intensity = 
    sub->Spgs->totalmass * cov->p[RANDOMCOIN_INTENSITY][0]; 

  //  PMI(cov, "randmcoin");
  assert(sub->mpp.moments >= 2);
  if (!R_FINITE(sub->Spgs->totalmass) || !R_FINITE(sub->mpp.M[2]))
    SERR("Moments of submodels not known");

  sub->role = ROLE_POISSON_GAUSS;
  return NOERROR;
}



void do_randomcoin(cov_model *cov, storage *s) {   
  assert(cov->simu.active);
  bool loggauss = (bool) ((int*) cov->p[LOG_GAUSS])[0];
  location_type *loc = Loc(cov);
  int i;
  double *res = cov->rf;

  dompp(cov, cov->stor==NULL ? s : cov->stor);// letzteres falls shape gegeben
  if (loggauss) {
    int vdimtot = loc->totalpoints * cov->vdim;
    for (i=0; i<vdimtot; i++) res[i] = exp(res[i]);
  }
}


////////////////////////////////////////////////////////////////////
// Schlather

void extremalgaussian(double *x, cov_model *cov, double *v) {
  // schlather process
  cov_model *next = cov->sub[0];
  
  if (cov->role == ROLE_SCHLATHER) COV(x, next, v)
  else {
    COV(x, next, v);
    *v = 1.0 - sqrt(0.5 * (1.0 - *v));
  }
}

//#define MPP_SHAPE 0
//#define MPP_TCF 1

int SetGEVetc(cov_model *cov, int role) {
  int err;
  extremes_param *ep = &(GLOBAL.extreme);
  
  if (cov->role != ROLE_COV) cov->role = role;

  if (cov->sub[MPP_TCF] != NULL && cov->sub[MPP_SHAPE]!=NULL)
    SERR("either 'tcf' or a shape definition must be given");
  if ((err = checkkappas(cov, false)) != NOERROR) return err;
  kdefault(cov, GEV_XI, ep->GEV_xi);
  kdefault(cov, GEV_S, 1.0);
  kdefault(cov, GEV_MU, cov->p[GEV_XI][0] / cov->p[GEV_S][0]);
  if ((err = checkkappas(cov, true)) != NOERROR) return err;
  if (cov->tsdim != cov->xdimprev || cov->tsdim != cov->xdimown) {
    return ERRORDIM;
  }
  return NOERROR;
}


int check_schlather(cov_model *cov) {

  //print("check schlather\n");

  cov_model
    *key = cov->key,
    *next = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE:MPP_TCF];
  int err,
    dim = cov->tsdim; // taken[MAX DIM],
  // mpp_param *gp  = &(GLOBAL.mpp);
  double v;

  //  printf("A \n");
 
  ASSERT_ONE_SUBMODEL(cov);
  ///  printf("B A \n");
  if ((err = SetGEVetc(cov, ROLE_SCHLATHER)) != NOERROR) return err;
  //   printf("CA \n");

  cov_model *sub = cov->key==NULL ? next : key;
  //        printf("GA \n");
       
  if (key == NULL) {
    // printf("key=NULL\n");
    int  role = isNegDef(sub) ? ROLE_COV  //Marco, 29.5.13
      : isShape(sub) ? ROLE_SCHLATHER
      : isGaussProcess(sub) ? ROLE_GAUSS
      : isBernoulliProcess(sub) ? ROLE_BERNOULLI
      : ROLE_UNDEFINED;    
    ASSERT_ROLE_DEFINED(sub);

    if ((err = isPosDef(next)
	 ? CHECK(next, dim, dim, PosDefType, XONLY, ISOTROPIC,
		   SCALAR, role)
	 : CHECK(next, dim, dim, ProcessType, XONLY, NO_ROTAT_INV, 
		   SCALAR, role)) != NOERROR)  return err;

    setbackward(cov, sub);

    if (role==ROLE_COV) {
      if (next->pref[Nothing] == PREF_NONE) return ERRORPREFNONE;
      COV(ZERO, next, &v);
      if (v != 1.0) 
	SERR("extremalgaussian requires a correlation function as submodel.");
    }
    
  } else { // key != NULL
    //printf("key!=NULL\n");
   if ( (err = CHECK(key, dim,  dim, PointShapeType, XONLY, NO_ROTAT_INV,  
		      SUBMODEL_DEP, ROLE_SCHLATHER)) != NOERROR) {
     //APMI(cov);
      return err;
    }
    setbackward(cov, sub);
  }
  //print("end check schlather\n");

  return NOERROR;
}


int struct_schlather(cov_model *cov, cov_model **newmodel){
  cov_model
    *sub = cov->sub[cov->sub[MPP_SHAPE] != NULL ? MPP_SHAPE:MPP_TCF];
  int err, role, ErrNoInit;

  if (cov->role != ROLE_SCHLATHER) BUG;
  if (newmodel != NULL)
    SERR1("unexpected structure request for '%s'", NICK(cov));
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if (cov->sub[MPP_TCF] != NULL) {
    if ((err = STRUCT(sub, &(cov->key))) >  NOERROR) return err;
    cov->key->calling = cov;
  } else {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) return err;
  }
  assert(cov->key->calling == cov);
 
  if (cov->key->nr != GAUSSPROC && !isBernoulliProcess(cov->key)) {
    if (isNegDef(cov->key)) addModel(&(cov->key), GAUSSPROC); 
    else {
      if (isGaussProcess(cov->key)) {
	SERR("invalid model specification");
      } else SERR("schlather process currently only allowed for gaussian processes and binary gaussian processes");
    }
  }
  assert(cov->key->calling == cov);
  
  role = cov->key->nr == GAUSSPROC ? ROLE_GAUSS   //Marco, 29.5.13
    : isBernoulliProcess(cov->key) ? ROLE_BERNOULLI
    : ROLE_UNDEFINED;

  //PMI(cov->key, "role ");
  ASSERT_ROLE_DEFINED(cov->key);

  if ((err = CHECK(cov->key, cov->tsdim, cov->xdimown, ProcessType, 
		     cov->domown,
		     cov->isoown, cov->vdim, role)) != NOERROR) return err;  

  // APMI(cov->key);
 
  if ((ErrNoInit = STRUCT(cov->key, NULL)) > NOERROR) return ErrNoInit;
  
  addModel(&(cov->key), STATIONARY_SHAPE);
  
  if ((err = CHECK(cov->key, cov->tsdim, cov->xdimown, PointShapeType, 
		     cov->domown,
		     cov->isoown, cov->vdim, ROLE_SCHLATHER)) != NOERROR)
    return err;
  
  return ErrNoInit;
}


////////////////////////////////////////////////////////////////////
// Smith


int check_smith(cov_model *cov) {
  cov_model 
    *shape = cov->sub[MPP_SHAPE],
    *TCF =  cov->sub[MPP_TCF],
    *next = shape != NULL ? shape : TCF,
    *key = cov->key,
    *sub = (key==NULL) ? next : key;
  ;
  int err, role,
    dim = cov->tsdim; // taken[MAX DIM],
  //Types type;

  ASSERT_ONE_SUBMODEL(cov);
  
  if ((err = SetGEVetc(cov, ROLE_MAXSTABLE)) != NOERROR) return err;

  if (key != NULL) {
    if ((err = CHECK(key, dim,  dim, PointShapeType, XONLY, NO_ROTAT_INV,
		       SUBMODEL_DEP, ROLE_SMITH)) != NOERROR) return err;
  } else { // key == NULL
    if (next == cov->sub[MPP_TCF]) {

      //      PMI(cov);

      if ((err = CHECK(next, dim, dim, TcfType, XONLY, ISOTROPIC, 
			 SCALAR, ROLE_SMITH)) != NOERROR) {
	return err;
      }

      if ((dim == 1 && next->rese_derivs < 1) ||
	  (dim >= 2 && dim <= 3 && next->rese_derivs < 2) ||
	  dim > 3) {

	//printf("dim = %d\n", dim);
        //APMI(next);
	SERR("submodel does not have enough derivatives (programmed).");
      }

    } else { // shape
      role = // (isNegDef(sub) && !isPosDef(sub)) ?  ROLE_COV  //Marco, 29.5.13
	isShape(sub) ? ROLE_MAXSTABLE
	: isPointShape(sub) ? ROLE_SMITH
	: isGaussProcess(sub) ? ROLE_GAUSS
	: isBernoulliProcess(sub) ? ROLE_BERNOULLI
	: ROLE_UNDEFINED; 
      ASSERT_ROLE_DEFINED(sub);
 
      if ((err = CHECK(next, dim, dim, ShapeType, XONLY, NO_ROTAT_INV, 
		       SCALAR, role)) != NOERROR) {
	return err;
      }
      if (next->full_derivs < 0) 
	SERR1("'%s' requires a an explicit submodel.", NICK(cov));
    }
  } 
  
  setbackward(cov, next);

  return NOERROR;
}

 


void param_set_identical(cov_model *to, cov_model *from, int depth) {
  if (depth <= 0) return;
  assert(to != NULL && from != NULL);
  assert(strcmp(NICK(from), NICK(to)) == 0); // minimal check

  cov_fct *C = CovList + from->nr;
  //*TO = CovList + to->nr;
  int i;

  assert(from->qlen == to->qlen);
  if (from->q != NULL) MEMCOPY(to->q, from->q, sizeof(double) * from->qlen);

  for (i=0; i<MAXPARAM; i++) {
    int n;
    assert(from->ncol[i] == to->ncol[i]);
    assert(from->nrow[i] == to->nrow[i]);
    assert(C->kappatype[i]==CovList[to->nr].kappatype[i]);
    if (C->kappatype[i]==REALSXP) {
      n = sizeof(double);
    } else if (C->kappatype[i]==INTSXP) {
      n = sizeof(int);
    } else n = -1;
    if (n > 0) {
      n *= from->nrow[i] * from->ncol[i];
      MEMCOPY(to->p[i], from->p[i], n);
    }
  }
 
  for (i=0; i<MAXSUB; i++) {
    if (from->sub[i] != NULL) 
      param_set_identical(to->sub[i], from->sub[i], depth - 1);
  }
 
}


int PointShapeLocations(cov_model *key, cov_model *shape) {
  int err,
    pgs = key->nr;
  if (key->sub[PGS_LOC] != NULL) return NOERROR;
  if (pgs == PTS_GIVEN_SHAPE) {   
    if ((err = covcpy(key->sub + PGS_LOC, shape)) != NOERROR) return err;
    addModel(key->sub + PGS_LOC, SETPARAM);
    cov_model *set = key->sub[PGS_LOC];
    assert(SETPARAM_TO == 0);
    if (set->Sset != NULL) SET_DELETE(&(set->Sset));
    set->Sset = (set_storage *) MALLOC(sizeof(set_storage));
    SET_NULL(set->Sset);
    set_storage *S = set->Sset;
    S->from = key->sub[PGS_FCT];
    S->set = param_set_identical;    
    S->variant = MAXINT;
  } else if (pgs == STANDARD_SHAPE) {
    assert(key != NULL);
    if ((err = STRUCT(shape, key->sub + PGS_LOC)) != NOERROR) return err;
  } else BUG;
  return NOERROR;
}

int addLocations(cov_model **Key, cov_model *shape, int dim, int vdim) {
#define PGS_N 2
  int err, i,
    pgs[PGS_N] = {PTS_GIVEN_SHAPE, STANDARD_SHAPE};
  
  assert(*Key == NULL);
  for (i=0; i<PGS_N; i++) {  
    if (*Key != NULL) COV_DELETE(Key);
    covcpy(Key, shape);
    addModel(Key, pgs[i]);
    assert((*Key)->sub[PGS_LOC] == NULL && (*Key)->sub[PGS_FCT] != NULL);
    if ((err = PointShapeLocations(*Key, shape)) != NOERROR) continue;

    if ((err = CHECK(*Key, dim, dim, PointShapeType, XONLY, NO_ROTAT_INV,
		     vdim, ROLE_MAXSTABLE)) != NOERROR) continue; 
    (*Key)->stor = (storage *) MALLOC(sizeof(storage)); 
    STORAGE_NULL((*Key)->stor);
    if ((err = INIT(*Key, 1, (*Key)->stor)) != NOERROR) continue;
  }
     
  return err;
}

int addPointShape(cov_model **Key, cov_model *shape, cov_model *pts, 
		  int dim, int vdim) {
#define PGS_N 2
  int err, i,
    pgs[PGS_N] = {PTS_GIVEN_SHAPE, STANDARD_SHAPE};
  
  assert(*Key == NULL);
  for (i=0; i<PGS_N; i++) {  
    if (*Key != NULL) COV_DELETE(Key);
    covcpy(Key, shape);
    addModel(Key, pgs[i]);
    covcpy((*Key)->sub + PGS_LOC, pts);
    (*Key)->sub[PGS_LOC]->calling = *Key;
    assert((*Key)->sub[PGS_LOC] != NULL && (*Key)->sub[PGS_FCT] != NULL);
    if ((err = CHECK(*Key, dim, dim, PointShapeType, XONLY, NO_ROTAT_INV,
		     vdim, ROLE_MAXSTABLE)) != NOERROR) continue; 
    (*Key)->stor = (storage *) MALLOC(sizeof(storage)); 
    STORAGE_NULL((*Key)->stor);
    if ((err = INIT(*Key, 1, (*Key)->stor)) != NOERROR) continue;
  }
    
  return err;
}

int struct_smith(cov_model *cov,  cov_model **newmodel){
  cov_model
    *tmp_shape = NULL,
    *shape = cov->sub[MPP_SHAPE],
    *tcf =  cov->sub[MPP_TCF],
    *dummy = NULL,
    *tcf_shape = NULL,
    *sub = shape != NULL ? shape : tcf;
  location_type *loc = Loc(cov);
  int err = NOERROR;
  if (newmodel != NULL) SERR("unexpected call of struct_smith");

  if (cov->role != ROLE_SMITH) BUG;  

  if (loc->Time || (loc->grid && loc->caniso != NULL)) {
    Transform2NoGrid(cov, false, GRIDEXPAND_AVOID);
    SetLoc2NewLoc(sub, Loc(cov));
  }

  if (cov->key != NULL) COV_DELETE(&(cov->key));
  assert(cov->key == NULL);
  
  if (tcf != NULL) {
    // to do: ausbauen 
    if ((err = covcpy(&tcf_shape, sub)) != NOERROR) return err;   
    addModel(&tcf_shape, STROKORB_MONO);

    if ((err = CHECK(tcf_shape, tcf->tsdim, tcf->xdimprev, ShapeType, 
		     tcf->domprev, tcf->isoprev, tcf->vdim, 
		     ROLE_MAXSTABLE)) != NOERROR) goto ErrorHandling;
    tmp_shape = tcf_shape; 
  } else {
    tmp_shape = shape;
  }

  
  if ((err = STRUCT(tmp_shape, &(cov->key))) == NOERROR && cov->key != NULL) {
    // APMI(cov);
    
    cov->key->calling = cov;
    if (TypeConsistency(PointShapeType, cov->key)) {
      if ((err = PointShapeLocations(cov->key, tmp_shape)) != NOERROR) 
	goto ErrorHandling;
    } else {
      dummy = cov->key;
      cov->key = NULL;
      if (TypeConsistency(RandomType, dummy)) { // random locations given
	// so, it must be of pgs type (or standard):
	err = addPointShape(&(cov->key), tmp_shape, dummy, 
			    cov->tsdim, cov->vdim);
	if (cov->key == NULL) BUG;
	cov->key->calling = cov;
	// APMI(cov);
      } else if (TypeConsistency(ShapeType, cov->key)) {
	// suche nach geeigneten locationen
	if ((err = addLocations(&(cov->key), dummy, cov->tsdim, cov->vdim))
	    != NOERROR) goto ErrorHandling; 
      } else BUG;
    } // !TypeConsistency(PointShapeType, cov->key)
  } else { // STRUCT does not return anything
    if ((err = addLocations(&(cov->key), tmp_shape, cov->tsdim, cov->vdim))
	!= NOERROR) goto ErrorHandling; 
  }

  err = NOERROR;

 ErrorHandling:
  if (tcf_shape != NULL) COV_DELETE(&tcf_shape);
  if (dummy != NULL) COV_DELETE(&dummy);
 
  return err;
}
