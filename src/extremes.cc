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


#include <Rmath.h>
#include <stdio.h>
#include "questions.h"
#include "Processes.h"
#include "shape.h"
#include "operator.h"

 
#define POISSON_INTENSITY 0
#define RANDOMCOIN_INTENSITY (COMMON_GAUSS + 1) 
#define RANDOMCOIN_METHOD (COMMON_GAUSS + 2) 

/* not used 
double GetTotal(model* cov) {
  defn *C = DefList + COVNR;
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
	FREE(r);
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
*/




  void range_mpp(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
    //   assert(cov != NULL);
  range->min[GEV_XI] = RF_NEGINF;
  range->max[GEV_XI] = RF_INF;
  range->pmin[GEV_XI] = - 10;
  range->pmax[GEV_XI] = 10;
  range->openmin[GEV_XI] = true;
  range->openmax[GEV_XI] = true; 

  range->min[GEV_MU] = RF_NEGINF;
  range->max[GEV_MU] = RF_INF;
  range->pmin[GEV_MU] = - 1000;
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


int init_mpp(model *cov, gen_storage *S) { // cov ist process
  model *sub = cov->key;
  if (sub == NULL) sub = cov->sub[cov->sub[0] == NULL];
  if (sub == NULL) SERR("substructure could be detected");
  location_type *loc = Loc(cov);
  //window_info *w = &(S->window);
  int err;
  bool
    poisson = hasPoissonFrame(sub),
    maxstable = hasMaxStableFrame(sub);
  pgs_storage *pgs;


  if (!maxstable && !poisson) ILLEGAL_FRAME;
  if (!equalsnowPointShape(sub))
    SERR1("%.50s is not shape/point function", NICK(sub));
  
  if (loc->distances) RETURN_ERR(ERRORFAILED);

  if ((err = INIT(sub, maxstable ? 1 : hasPoissonFrame(sub) ? 0 : 2, S)) 
      != NOERROR)  RETURN_ERR(err);
  pgs = sub->Spgs;

  assert(pgs != NULL);
   
  GetDiameter(loc, pgs->supportmin, pgs->supportmax, pgs->supportcentre); 
  
  if (maxstable) {      
    if (!R_FINITE(sub->mpp.mMplus[1]) || sub->mpp.mMplus[1] <= 0.0) {
      SERR("integral of positive part of submodel unkown");
    }

    if (!R_FINITE(pgs->zhou_c) && !R_FINITE(pgs->sum_zhou_c))
      SERR("maximal height of submodel unkown or the set of locations exceeds possibilities of the programme");
  } else if (hasPoissonGaussFrame(sub)) {
    if (ISNAN(sub->mpp.mM[2]) || !R_FINITE(sub->mpp.mM[2] || sub->mpp.mM[2] <=0.0)){
      SERR("second moment unkown");
    }
  } else assert(hasPoissonFrame(sub)); 
 
  if ((err = ReturnOwnField(cov)) != NOERROR) { // must be later than INIT !
    RETURN_ERR(err);  
  }
  
  cov->initialised = cov->simu.active = err == NOERROR;

  RETURN_ERR(err);
}

#define SET_SUB	{							\
    sub = cov->key;  /* if changed, change do_max_pgs! */		\
    if (sub == NULL) sub = cov->sub[cov->sub[0] == NULL];		\
    if (sub == NULL) ERR("substructure could not be detected");		\
    pgs = sub->Spgs;							\
    assert(pgs != NULL);						\
    gridlen = pgs->gridlen;						\
    end = pgs->end;							\
    start = pgs->start;							\
    delta = pgs->delta;							\
    nx = pgs->nx;							\
    xstart = pgs->xstart; /* nach DO gesetzt */				\
    x = pgs->x;           /* nach DO gesetzt */				\
    inc = pgs->inc;}

 

void dompp(model *cov, gen_storage *s, double *simuxi) { // cov ist process
  assert(cov->simu.active);
  location_type *loc = Loc(cov);	
  if (Getcaniso(cov) != NULL) BUG;

  model *sub = NULL,
    *key = cov->key;
  pgs_storage *pgs = NULL;
  long spatial = Getspatialpoints(cov);
  int  err,
    *gridlen = NULL,
    *end = NULL,
    *start = NULL,
    *delta = NULL,
    *nx = NULL,
    dim = ANYOWNDIM,
     every = GLOBAL.general.every, //
    nthreshold = (every>0) ? every : MAXINT,	
    deltathresh = nthreshold;			
 
  double logdens,  factor, summand, Sign[2], value[2], logthreshold, 
    *xstart = NULL, // nach DO gesetzt
    *x = NULL,           // nach DO gesetzt
    *inc = NULL,
    *rf = NULL,
    //M1 = RF_NA,
    logM2 = RF_NA,
    gumbel = RF_NA,
    frechet = RF_NA,
    threshold = RF_NA,
    poisson=0.0,
    Minimum = GLOBAL.extreme.min_shape_gumbel, // -10 // TO DO; -1e150
    *res = cov->rf;
  coord_type xgr = Getxgr(cov);

  assert(res != NULL);
  ext_bool fieldreturn, loggiven;
  bool   
    maxstable = hasMaxStableFrame(key),
    simugrid = Getgrid(cov),
    partialfield = false,
    // AVERAGE braucht Sonderbehandlung:
     randomcoin = cov->method == RandomCoin //only for FRAME_ POISSON_ GAUSS
    ;
  long zaehler, nn, cumgridlen[MAXMPPDIM +1], Total_n,			
    total = Gettotalpoints(cov), 
    vdim = VDIM0,
    vdimtot = total * vdim,
    control = 0;
  assert(VDIM0 == VDIM1);

 
  if (vdim != 1)
    ERR("Poisson point process based methods only work in the univariate case");

  SET_SUB;
  if (!equalsnowPointShape(sub)) BUG;					       
  bool compoundpoisson = hasPoissonFrame(sub),
    poissongauss = hasPoissonGaussFrame(sub);

  
  Total_n = -1;
  if (!maxstable) Total_n = rpois(pgs->intensity);

  // printf("pg = %d %10g %d\n", poissongauss, pgs->intensity, maxstable);
  // A  PMI0(sub); 

  assert(maxstable || Total_n > 0);
  loggiven = key == NULL ? cov->sub[0]->loggiven : key->loggiven;
  if (!loggiven) Minimum =EXP(Minimum);
   
  if (simugrid) {
    cumgridlen[0] = 1;
    for (int d=0; d<dim; d++) {
      inc[d] = xgr[d][XSTEP];
      gridlen[d] = (int) xgr[d][XLENGTH]; 
      cumgridlen[d+1] = gridlen[d] * cumgridlen[d];

      // printf("set d=%d inc=%10g gridlen=%d cum=%d\n", d, inc[d], gridlen[d], (int) cumgridlen[d]);
    }
  }
  
  if (maxstable) {
    if (loggiven) for (int i=0; i<vdimtot; i++) res[i] = RF_NEGINF;
    else for (int i=0; i<vdimtot; i++) res[i] = 0.0;
    if (sub->mpp.moments < 1 || !R_FINITE(sub->mpp.mMplus[1]) || 
      sub->mpp.mMplus[1] == 0.0)
      ERR("integral of positive part of the submodel is unknown or not finite");
    //    M1 = sub->mpp.mMplus[1];
    threshold = logthreshold = 1e15; //RF_INF;
    pgs->globalmin = Minimum;
    pgs->currentthreshold =
      loggiven ? pgs->globalmin - threshold : pgs->globalmin / threshold;
    //    printf("pgs->current = %10g loggiven=%d %10g %10g\n", pgs->currentthreshold, loggiven, threshold, pgs->globalmin);
    //    BUG;
    
    if (simuxi != NULL) {
      RFERROR("extremal t not programmed yet");
      // to do : RPstudent, extremal t process
    }
  } else if (compoundpoisson){
    assert(simuxi == NULL);
    for (int i=0; i<vdimtot; i++) res[i] = 0.0;
  } else {
    assert(simuxi == NULL);
    for (int i=0; i<vdimtot; res[i++] = 0.0);
    logM2 =LOG(sub->mpp.mM[2]);
    pgs->currentthreshold = GLOBAL.mpp.about_zero;
  }

  for(nn=1; ; nn++) {
    //printf("%d %d tot=%d\n", (int) nn, (int) control, (int) Total_n); assert(nn < 1000);
    //       if (n % 1000 == 0)
    // printf("n=%d tot=%d  contr=%d: %10g < %10g; res[0]=%10g\n", (int) n, (int) Total_n,  (int) control, res[control], threshold, res[0]);
    // assert(n <= 10);
    

    if (!maxstable && (nn > Total_n)) break; 
    
    // for (d=0; d<dim; d++) info->min[d] = info->max[d] = 0.0;
    
    if (sub->randomkappa) {
      if ((err = REINIT(sub, sub->mpp.moments, s)) != NOERROR) BUG;
      DO(sub, s);
      SET_SUB;
    } else DO(sub, s);
  
    //    printf("%10g %d\n", sub->rf[0], sub->fieldreturn);
    //APMI0(sub);
   
    fieldreturn = sub->fieldreturn;
    // MARCO: hattest Du dies reingenommen?
    //   if (subf rame == BrMethodType) {
    //      assert(sub->SBR != NULL);
    //      n += sub->SBR->hatnumber;
    //    }
    
    //printf("log_den %10g\n",  pgs->log_density);
    logdens = pgs->log_density;
    
    if (!R_FINITE(logdens)) {
      BUG;
    //    get_shape = DefList[sub->gatternr].cov;
    //   get_logshape = DefList[sub->gatternr].log;
    }
 
    if (compoundpoisson) {
      summand = 0.0;
    } else if (poissongauss) {
      summand = -0.5 * (logdens + logM2);
    } else { // max-stable field          
      assert(hasMaxStableFrame(sub));

      if (nn > GLOBAL.extreme.maxpoints) {
	PRINTF("'maxpoints' (= %d) exceeded. Simulation is likely to be a rough approximation only. Consider increasing 'maxpoints'.\n", 
	       GLOBAL.extreme.maxpoints
		);
       	break;
      } 
      poisson += rexp(1.0);
      frechet =POW(poisson, - 1.0 / pgs->alpha);
      
      gumbel = -LOG(poisson);
      threshold = logthreshold = (double) (gumbel+LOG(sub->mpp.maxheights[0])); 
      if (!loggiven) threshold = EXP(threshold); // Frechet
 
      summand = gumbel - logdens; //@MARCO: summand ist eine rein
      // technische Variable, die die Anzahl der Berechnungen ein
      // winziges bisschen reduziert im Falle des Zhou-Ansatzes.

      //printf("logdens = %10g\n", logdens);BUG;


      //printf("for thres=%10e %d res=%10e zhou=%10g(%10g) logden=%10g  gumbel=%10g loggiven=%d summand=%10g\n", threshold, (int) control, res[control], pgs->zhou_c, sub->mpp.maxheights[0], logdens, gumbel, loggiven, summand);// assert(false);

      //        assert(R_FINITE(threshold ));

// {   int i; for (i=0; i<total; i++) printf("%10e ", res[i]); }


      while ((control<total) && (res[control]>=threshold)) {
	//	printf("control: %d %10g %10g\n",  (int) control, res[control], threshold); //assert(false);
	control++;
      }
      if (control >= total) {
	break;
      }
      
       
       //     printf("n=%d every =%d %d %d\n", (int) n, GLOBAL.extreme.check_every, PL, PL_STRUCTURE);
      if (nn % GLOBAL.extreme.check_every == 0) {
	double min = res[control];
	for (int i=control + 1; i<total; i++)  {
	  //printf("i=%d control=%d total=%d %10g\n", i, control, total,  res[i]);
	  min = res[i] < min ? res[i] : min;
	} 
	  
	if (PL >= PL_STRUCTURE) {
	  PRINTF("control: %ld %10g %10g glob=%10g n=%ld every=%d %10g log.max=%10g\n", // %ld OK
		 control, res[control], threshold, min,
		 nn, GLOBAL.extreme.check_every, sub->mpp.maxheights[0],
		 LOG(sub->mpp.maxheights[0]));	 
	}
	pgs->globalmin = min < Minimum ? Minimum : min;     
      }
 	
      pgs->currentthreshold = loggiven
	? pgs->globalmin - gumbel
	: pgs->globalmin / frechet;
      //   printf("xxx = %10e %10e %10e Min=%10e \n",pgs->globalmin,  gumbel, pgs->currentthreshold, Minimum);
      
    } // maxstable
    factor =EXP(summand);
    
  // printf("factor %4.4f sumd=%4.4f logdens=%4.4f logM2=%4.4f th=%4.4f %d<%d \n", factor, summand, logdens, logM2, threshold, (int) control, (int) total);  //  assert(false);
    
  //   printf("dim=%d loggiven=%d maxstable=%d simugrid=%d randomcoin=%d\n", dim, loggiven, maxstable, simugrid, randomcoin);// assert(false);
    
    //printf("A=%d\n", n);
    zaehler = 0;
    if (simugrid) {
      partialfield = false;
      int d;
      for (d=0; d<dim; d++) {	
	double
	  step = inc[d] == 0.0 ? 1.0 : inc[d], 
	  dummy = CEIL((pgs->supportmin[d] - xgr[d][XSTART]) / step
		       - 1e-8);	
	start[d] = dummy < 0 ? 0 : dummy > gridlen[d] ? gridlen[d] : (int)dummy;
	partialfield |= start[d] > 0;
	dummy = TRUNC(1.00000001 + 
		      ((pgs->supportmax[d] - xgr[d][XSTART]) / step ));
	end[d] = dummy < 0 ? 0 : dummy > gridlen[d] ? gridlen[d] : (int) dummy;
	partialfield |= end[d] < gridlen[d];

	// 	printf("window [%10g %10g] [%d %d] %d %d\n", pgs->supportmin[d], pgs->supportmax[d], start[d], end[d], gridlen[d], partialfield);  //assert(n < 150); //APMI(cov);


	if (start[d] >= end[d]) { // 
	  //	  printf("broken %d %d %d supp=%10g %10g, loc.start=%10g %10g\n", d, start[d], end[d], pgs->supportmin[d], pgs->supportmax[d], xgr[d][XSTART], inc[d]); 
	  break;
	}
	delta[d] = (end[d] - start[d]) * cumgridlen[d];
	nx[d] = start[d];
	zaehler += start[d] * cumgridlen[d];
	x[d] = xstart[d] = xgr[d][XSTART] + (double) start[d] * inc[d]; 

	if (false || zaehler < 0) {
	  PRINTF("d=%d start=%d, end=%d xstart=%10g %10g pgs=[%4.3f, %4.3f] xgr=%10g %10g %10g inc=%3.2f\n", 
	  	 d, start[d], end[d], xstart[d], pgs->xstart[d], pgs->supportmin[d], pgs->supportmax[d],
	 	 xgr[d][XSTART], xgr[d][XSTEP], xgr[d][XLENGTH],
	  	 inc[d]);
	   // 	  assert(false);
	}
      }

      if (d < dim) {
	//printf("cont: d=%d\n", d);
	continue;
      }
    }


#define STANDARDINKREMENT_ZAEHLER			\
    int d = 0;						\
    nx[d]++;						\
    x[d] += inc[d];					\
    zaehler += cumgridlen[d];				\
    while (nx[d] >= end[d]) {				\
      nx[d] = start[d];					\
      x[d] = xstart[d];					\
      zaehler -= delta[d]; /* delta is positive */	\
      if (++d >= dim) break;				\
      nx[d]++;						\
      x[d] += inc[d];					\
      zaehler += cumgridlen[d];				\
    }							\
    if (d >= dim) { break;}	
    // end define StandardInkrement

 
#define GET SHAPE(x, sub, value); /* // printf("%10g : %10g fctr=%10g %ld \n", x[0], *value, factor, zaehler); // */
#define LOGGET LOGSHAPE(x, sub, value, Sign);
#define LOGGETPOS LOGGET if (Sign[0] > 0)
#define GETRF *value = (double) (rf[zaehler]);


#define RFMAX(OP) {				\
      rf = sub->rf;				\
      for (int i=0; i<total; i++) {		\
	double val = (double) (OP);			\
	res[i] = res[i] < val ? val : res[i];		\
      }						\
    }

#define RFSUM(OP) {				\
      rf = sub->rf;				\
      for (int i=0; i<total; i++) {		\
	double val = (double) (OP);		\
	/*	if (res[i] < *value) res[i] += OP;  */	\
	res[i] += val;				\
      }						\
    }



#define MAXAPPLY(GET, OP)						\
      assert(zaehler >=0 && zaehler<total);				\
      /* if (zaehler<0 || zaehler>=total) { PMI(cov);  printf("//z=%d n=%d\n", zaehler, n); } */ \
      GET {								\
	OP;								\
	/*  assert( {if (*value > sub->mpp.maxheights[0] *EXP(gumbel) * (1+1e-15)) {  printf("//r=%10g %10g delta=%10e %d %10e\n", *value /EXP(gumbel), sub->mpp.maxheights[0], *value /EXP(gumbel) - sub->mpp.maxheights[0], loggiven, factor); };true;}); */ \
	/*assert(loggiven || *value <= sub->mpp.maxheights[0] *EXP(gumbel) * (1+1e-15)); */ \
	res[zaehler] = res[zaehler] < *value ? *value : res[zaehler];	\
      }								       
   
#define SUMAPPLY(GET,OP) GET; res[zaehler] += OP;

#define APPLY(GET, OP, WHAT)	{				\
      if (simugrid) {						\
	while(true) {						\
	  WHAT(GET, OP);					\
	  STANDARDINKREMENT_ZAEHLER;				\
	}							\
      } else {							\
	for(x=loc->x; zaehler<spatial; zaehler++, x += dim) {	\
	  WHAT(GET, OP);					\
	}							\
      }								\
    }
 

#define OG_BASIC(OP)						\
    assert(zaehler >=0 && zaehler < total);			\
    int jj, idx, index=0;					\
    for (jj=0; jj<dim; jj++) {					\
      idx = ROUND((x[jj] - ogstart[jj]) / ogstep[jj]);		\
      if (idx < 0 || idx >= pgs->own_grid_len[jj]) break;	\
      index += pgs->own_grid_cumsum[jj] * idx;			\
    }								\
    if (jj >= dim) {  OP;  } else if (false) printf("%d %10g %10g %10g len=%10g %d %d\n", jj, x[jj],  ogstart[jj], ogstep[jj], pgs->own_grid_len[jj], pgs->own_grid_cumsum[jj], idx);  /* // */  \
   
  				

#define OG_MAXAPPLY(NONE, OP)						\
    rf = sub->rf;							\
    OG_BASIC(*value = rf[index];					\
    OP; res[zaehler] = res[zaehler] < *value ? *value : res[zaehler];)
    
#define OG_SUMAPPLY(NONE, OP)				\
    rf = sub->rf;					\
    OG_BASIC(*value = rf[index]; res[zaehler] += OP;) 
    
    //   OWNGRIDAPPLY(OP, OG_SUMAPPLY_BASIC)
#define OWNGRIDAPPLY(NONE, OP, WHAT) {		\
      double *ogstart = pgs->own_grid_start,	\
	*ogstep = pgs->own_grid_step;		\
      APPLY(NONE, OP, WHAT);			\
    }


#define	ALL_APPLY(APPLY, MGET, MAPPLY, SGET, SAPPLY, NONLOGGET)		\
    if (loggiven == true) {						\
      if (maxstable) APPLY(MGET, (*value) += summand, MAPPLY)		\
	else APPLY(SGET, Sign[0] * EXP(*value + summand), SAPPLY);	\
    } else { /* not loggiven */						\
      if (sub->loggiven){/* the realized loggiven  */			\
	if (maxstable)							\
	  APPLY(MGET, *value=Sign[0] *EXP(*value+summand), MAPPLY)	\
	else APPLY(SGET, Sign[0] *EXP(*value + summand), SAPPLY);	\
      } else {								\
	assert(R_FINITE(factor));					\
	if (maxstable) APPLY(NONLOGGET, (*value) *= factor, MAPPLY)	\
	  else APPLY(NONLOGGET, *value * factor, SAPPLY);		\
      }									\
    }


#define AVEAPPLY(G,OP0,OP1) {						\
      warning(" * timescale ?!");					\
      if (simugrid) {							\
	while (true) {							\
	  double zw,inct, segt, ct, st, A, cit, sit;			\
	  inct = sub->q[AVERAGE_YFREQ] * xgr[dim][XSTEP];		\
	  cit = COS(inct);						\
	  sit = SIN(inct);						\
	  G;								\
	  segt = OP1;							\
	  A = OP0;							\
	  ct = A * COS(segt);						\
	  st = A * SIN(segt);						\
	  for (long zt = zaehler; zt < total; zt += spatial) {		\
	    res[zt] += (double) ct;					\
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
	  G; res[zaehler] += OP0 * COS(OP1);				\
	}								\
      }									\
    }	


    // TO DO: OMP: use parallel and SIMD max-operators
    
    // *atomdependent rfbased/not; (( recurrent/dissipative: min/max gesteuert))

    //  assert(maxstable);

       
    if (poissongauss  && !randomcoin) { // AVERAGE
      if (sub->loggiven) {
	AVEAPPLY(LOGGET, Sign[0] *EXP(value[0] + summand), 
		 Sign[1] *EXP(value[1]));
      } else AVEAPPLY(GET, value[0] * factor, value[1]);
    } else {
      // printf("here\n");
      assert(hasMaxStableFrame(sub) || hasPoissonFrame(sub)); // reicht das?
      if (fieldreturn == Huetchenownsize) { // atomdependent rfbased/not;   
	ALL_APPLY(OWNGRIDAPPLY, NULL, OG_MAXAPPLY, NULL, OG_SUMAPPLY, NULL);
      } else if (fieldreturn == wahr) { // extended boolean
	// todo : !simugrid! && Bereiche einengen!!     
	// note ! no Sign may be given currently when fieldreturn and loggiven!!
	if (partialfield) {
	  //!! Windows
	  if (!simugrid) BUG;
	  ALL_APPLY(APPLY, GETRF, MAXAPPLY, GETRF, SUMAPPLY, GETRF);
	} else { // !partialfield
	  if (sub->loggiven) {       
	    if (maxstable) RFMAX(rf[i] + summand)
	    else RFSUM(EXP(rf[i] + summand));	  
	  } else { // Linux
	    if (maxstable) RFMAX(rf[i] * factor)
	    else RFSUM(rf[i] * factor);
	  }
	}
      } else if (fieldreturn == falsch) { // not field return
	//	printf("kein feld\n");
	ALL_APPLY(APPLY, LOGGETPOS, MAXAPPLY, LOGGET, SUMAPPLY, GET);
      } else BUG; // neither Huetchenownsize, nor true nor false
    }
  
    if (nn > nthreshold) {
      if (maxstable) {
	LPRINT("%ld%.50s pos.: value=%1.3f threshold=%1.3f gu=%1.3f %1.3f %d (%ld %d)\n", // %ld ok
	       control, TH(control), (double) res[control], (double) threshold,
	       gumbel, logdens, loggiven,
	       nn, nthreshold); //, deltathresh); 
   	nthreshold += deltathresh;
      } else {
	PRINTF("n=%d Total=%d\n", (int) nn, (int)(Total_n));
      }
    }

    R_CheckUserInterrupt();
 
  }// for(n,...    --     while control < total


  // formerly (before 19 Aug 2017: contents of 
  // void finalmaxstable(model *cov, double *res, int n)
  // had been here --- much better to have it at the very end --
  // much less addional evaluations necessary to get c_zhou!


  // for (k=0; k<total; k++) printf("%10g ", res[k]); printf("\n");
  if (PL >= PL_DETAILSUSER) {
    PRINTF("number of shape functions used: %ld\n", nn);
  }

  return;
} // end dompp


void dompp(model *cov, gen_storage *s) {
  dompp(cov, s, NULL);
}


// 4805920678128 121



void finalmaxstable(model *cov, double *Res, int n, gen_storage *s) {
  
  // cov is process
  model *key = cov->key,
    *sub = cov->key;/* if changed, change do_max_pgs! */
  if (sub == NULL) sub = cov->sub[cov->sub[0] == NULL];			
  if (sub == NULL) ERR("substructure could be detected");		
  long
    total = Gettotalpoints(cov), 
    vdim = VDIM0,
    vdimtot = total * vdim;  
  double meansq, var,
     mean = RF_NA, 
    n_min = RF_INF,      
    eps = GLOBAL.extreme.eps_zhou,
    *endres = Res + n * vdimtot;
  ext_bool loggiven = key == NULL ? cov->sub[0]->loggiven : key->loggiven;
  pgs_storage *pgs = sub->Spgs;							
  assert(pgs != NULL);						

  
  if (pgs->estimated_zhou_c) {
    if (PL > 5) {
      PRINTF("zhou_c: %ld %d",pgs->n_zhou_c, GLOBAL.extreme.min_n_zhou);
    }
    
    while (pgs->n_zhou_c < GLOBAL.extreme.min_n_zhou) {
      int err;
      if ((err = REINIT(sub, sub->mpp.moments, s)) != NOERROR) BUG;
      DO(sub, s);   
    }
    
    while (true) { // n_min depends on the estimation precision
      mean = (double) (pgs->sum_zhou_c / (long double) pgs->n_zhou_c);
      meansq = mean * mean;
      var = (double) (pgs->sq_zhou_c / (long double) pgs->n_zhou_c - meansq);
      n_min = var / (meansq * eps * eps);
      
      if (PL > 5) {
	PRINTF(" n=%ld, min=%10g %10g mean=%10g zhou=%10g eps=%10g\n",
	       pgs->n_zhou_c,
	       n_min, (double) pgs->sum_zhou_c, mean, pgs->zhou_c, eps);
      }

      if (n_min >= GLOBAL.extreme.max_n_zhou || pgs->n_zhou_c >= n_min) break;

      int endfor=(n_min - pgs->n_zhou_c); 
      endfor = endfor < 10 ? 10 : endfor > 50 ? 50 : endfor;
      for (int k=0; k<endfor; k++) {
	int err;
	if ((err = REINIT(sub, sub->mpp.moments, s)) != NOERROR) BUG;
	DO(sub, s);
	// SET_SUB missing ??
      }
    }
    //
  } else {
    mean = pgs->zhou_c;
  }
    
    // mean /= M1;
    
    //printf(" n=%ld, min=%10g, mean=%10g zhou=%10g eps=%10g %10g max=%10g\n", pgs->n_zhou_c, n_min, (double) pgs->sum_zhou_c, mean, pgs->zhou_c, eps, sub->mpp.maxheights[0]);
  
  
    //  printf("loggiven = %d %10g %10g\n", loggiven, mean, pgs->zhou_c);

  double xi = P0(GEV_XI);
  if (loggiven && !pgs->logmean) mean =LOG(mean);  

  // sieht noch umstaendlich aus, aber jede variable kann ein anderes
  // xi irgendwann haben.
  for (double *res = Res; res<endres; res+=vdimtot) {
    if (pgs->alpha != 1.0) { // default = 1 in NULL.cc; s. opitz
      if (loggiven) for (int k=0; k<total; k++) res[k] *= pgs->alpha; 
      else for (int k=0; k<total; k++) res[k] =POW(res[k], pgs->alpha);      
    }
    //    printf("xi=%10g %d logmean=%d %10g\n", xi, loggiven, pgs->logmean, mean);
    
    if (loggiven) { // 3.6.13 sub->loggiven geaendert zu loggiven
      for (int k=0; k<total; k++) res[k] += mean;
      if (xi != 0.0) 
	for (int k=0; k<total; k++) res[k] =EXP(xi * res[k]) - 1.0;
      // @MARCO unten ein s/xi statt nur s wie erwartet
      //       
    } else {
      for (int k=0; k<total; k++) res[k] *= mean;  
      if (FABS(xi)<1e-15) for (int k=0; k<total; k++) res[k] = LOG(res[k]);
      else if (FABS(xi - 1.0) < 1e-15) for (int k=0; k<total; k++) res[k] -=1.0;
      else for (int k=0; k<total; k++) res[k] = POW(res[k], xi) - 1.0;
    }
  
    
    // muss das vorletzte sein !
    if (cov->kappasub[GEV_S] != NULL) {
      ERR1("'%.50s' cannot be chosen as a function yet.", KNAME(GEV_S));
      // check whether s is positive !!
    } else {
      double ss = xi == 0 ? P0(GEV_S) : (P0(GEV_S) / xi);
      //    printf("s/xi=%10g\n", ss);
      if (FABS(ss - 1.0) > 1e-15) for (int k=0; k<total; k++) res[k] *= ss;
    }
    
    
    // muss das letzte sein !
    if (cov->kappasub[GEV_MU] != NULL) {
      ERR1("'%.50s' cannot be chosen as a function yet.", KNAME(GEV_MU));
    } else {
      double mu = P0(GEV_MU);
      //   printf("mu=%10g xi=%10g %10g\n", mu, xi, s);
      //      printf("mu=%10g\n", mu);
      if (FABS(mu) > 1e-15) for (int k=0; k<total; k++) res[k] += mu;
    }
  }
} // finalmaxstable


//void calculate_integral(int d, double R, model *cov,
//			integral_type *integral) {
// 
//}




////////////////////////////////////////////////////////////////////
// Poisson


int addStandardPoisson(model **Cov) { // a submodel of poisson only
  model *cov = (*Cov)->calling, // cov is poisson process
    *next = *Cov;
  ASSERT_ONESYSTEM;
  int  err,
     dim = XDIM(PREVSYSOF(next), 0),
    vdim = next->vdim[0];
  Types
    frame = next->frame;     
  assert(VDIM0 == VDIM1);

  addModel(Cov, STANDARD_SHAPE, cov);
  model *key = *Cov;
  SetLoc2NewLoc(key, PLoc(cov));
  assert(key->calling == cov);
  assert(key->sub[PGS_LOC] == NULL && key->sub[PGS_FCT] != NULL);

#define checkstandardnext\
  if ((err = CHECK(key, dim, dim, PointShapeType, XONLY,	\
		   CoordinateSystemOf(OWNISO(0)),		\
		   vdim, frame)) != NOERROR) RETURN_ERR(err)

  checkstandardnext;
  if (!CallingSet(*Cov)) BUG;
  if (hasPoissonFrame(next)) {
    addModel(key, PGS_LOC, UNIF);
    assert(key->nsub == 2); // ???
    PARAMALLOC(key->sub[PGS_LOC], UNIF_MIN, dim, 1);
    PARAMALLOC(key->sub[PGS_LOC], UNIF_MAX, dim, 1);
  } else {
    // TREE(cov);
    // printf("structh %ld %.50s\n", key->sub + PGS_LOC, TYPE_NAMES[cov->frame]);
    
    if ((err = STRUCT(key, key->sub + PGS_LOC)) != NOERROR) RETURN_ERR(err);
    SET_CALLING(key->sub[PGS_LOC], key);
  }
  if (!CallingSet(*Cov)) BUG;
  checkstandardnext;
  
  RETURN_NOERROR;
}

int check_poisson(model *cov) {
  model *next=cov->sub[0],
    *key = cov->key,
    *sub = key == NULL ? next : key;
  int err,
    dim = ANYDIM; // taken[MAX DIM],
  mpp_param *gp  = &(GLOBAL.mpp);
  Types type = sub == key ? PointShapeType : ShapeType; 

  //print("eee\n");
  
  kdefault(cov, POISSON_INTENSITY, gp->intensity[dim]);
  if ((err = checkkappas(cov, true)) != NOERROR) RETURN_ERR(err);

  if ((err = CHECK(sub, dim, dim, type, XONLY, CoordinateSystemOf(OWNISO(0)),
		     SUBMODEL_DEP, PoissonType)) != NOERROR) RETURN_ERR(err);
  setbackward(cov, sub);
 
  RETURN_NOERROR;
}


void range_poisson(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[POISSON_INTENSITY] = RF_NEGINF;
  range->max[POISSON_INTENSITY] = RF_INF;
  range->pmin[POISSON_INTENSITY] = - 10;
  range->pmax[POISSON_INTENSITY] = 10;
  range->openmin[POISSON_INTENSITY] = true;
  range->openmax[POISSON_INTENSITY] = true; 
}


int struct_poisson(model *cov, model **newmodel){
  model *next=cov->sub[0];
  location_type *loc = Loc(cov);

  // if (newmodel != NULL) crash();
  ASSERT_NEWMODEL_NULL;
  assert(isnowProcess(cov));

  if (cov->key != NULL) COV_DELETE(&(cov->key));

  if (loc->Time || (loc->grid && loc->caniso != NULL)) {
    TransformLoc(cov, false, GRIDEXPAND_AVOID, false);
    SetLoc2NewLoc(next, PLoc(cov)); // passt das?
   }

  if (!equalsnowPointShape(next)) {
    int err;
    if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
    if ((err = addStandardPoisson(&(cov->key))) != NOERROR) RETURN_ERR(err);
  }
 
  RETURN_NOERROR;
}


int init_poisson(model *cov, gen_storage *S) {
  //  location_type *loc = Loc(cov);
  // mpp_storage *s = &(S->mpp);
  model *key=cov->key;
  int err;
 
  if ((err = init_mpp(cov, S)) != NOERROR) {
    RETURN_ERR(err);    
  }

  if (!equalsnowPointShape(key)) SERR("no definition of a shape function found");
  
  key->Spgs->intensity = key->Spgs->totalmass
    * P0(POISSON_INTENSITY);

  cov->simu.active = cov->initialised = err == NOERROR;

  RETURN_ERR(err);
}



////////////////////////////////////////////////////////////////////
// Schlather

void extremalgaussian(double *x, model *cov, double *v) {
  // schlather process
  model *next = cov->sub[0]; 
  COV(x, next, v);
  if (hasGaussMethodFrame(next)) *v = 1.0 - SQRT(0.5 * (1.0 - *v));
}
	

//#define MPP_SHAPE 0
//#define MPP_TCF 1

int SetGEVetc(model *cov) { // cov is process
  int err;
  extremes_param *ep = &(GLOBAL.extreme);

  if (cov->sub[MPP_TCF] != NULL && cov->sub[MPP_SHAPE]!=NULL)
    SERR2("either '%.50s' or '%.50s' must be given",
	  SNAME(MPP_TCF), SNAME(MPP_SHAPE));
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);
  kdefault(cov, GEV_XI, ep->GEV_xi);
  kdefault(cov, GEV_S, P0(GEV_XI) == 0.0 ? 1.0 : FABS(P0(GEV_XI)));
  kdefault(cov, GEV_MU, P0(GEV_XI) == 0.0 ? 0.0 : 1.0);

  // print("xi s mu %10g %10g %10g\n",  P0(GEV_XI), P0(GEV_S), P0(GEV_MU)); assert(false);

  if ((err = checkkappas(cov, true)) != NOERROR) RETURN_ERR(err);
  RETURN_NOERROR;
}


int check_schlather(model *cov) {
  //print("check schlather\n");
  if ((cov->sub[MPP_TCF] != NULL) xor (cov->sub[MPP_SHAPE] == NULL))
    SERR("two submodels given instead of one.");
 model
    *key = cov->key,
    *next = cov->sub[cov->sub[MPP_TCF] != NULL ? MPP_TCF : MPP_SHAPE];
  int err,
    dim = ANYDIM; // taken[MAX DIM],
  // mpp_param *gp  = &(GLOBAL.mpp);
  double v;
  bool schlather = DefList[COVNR].Init == init_mpp; // else opitz

  //  printf("A \n");
 
  ASSERT_ONE_SUBMODEL(cov);
  ///  printf("B A \n");
  if ((err = SetGEVetc(cov)) != NOERROR) RETURN_ERR(err);
  //   printf("CA \n");

  model *sub = cov->key==NULL ? next : key;
  //        printf("GA \n");
       
  if (key == NULL) {
    //    printf("key=NULL\n");
    if (equalsBernoulliProcess(sub) && !schlather)
      SERR1("'%.50s' does not work with Bernoulli processes", NAME(cov));
    Types frame = isProcess(sub) || isPointShape(sub) ? SchlatherType :
      EvaluationType;
      //isGaussMethod(sub) ? GaussMethodType
      // : isPointShape(sub) && schlather ? SchlatherType
      //: equalsBernoulliProcess(sub) && schlather ? NormedProcessType // egal
      //: isProcess(sub) ? SchlatherType
      // : EvaluationType;

    //  printf(">>> %.50s %d %.50s\n", TYPE_NAMES[frame], isProcess(sub), NAME(sub));
    //

    if (isProcess(next))
      err = CHECK(next, dim, dim, ProcessType, XONLY, 
		  CoordinateSystemOf(OWNISO(0)), SCALAR, frame);
    else
      err = CHECK(next, dim,dim, PosDefType, XONLY, IsotropicOf(OWNISO(0)),
		  SCALAR, frame);
    if (err != NOERROR) RETURN_ERR(err);

    if (sub->vdim[0] != 1) SERR("only univariate processes are allowed");
    assert(VDIM0 == VDIM1);

    setbackward(cov, sub);

    if (hasEvaluationFrame(next)) {
      if (next->pref[Nothing] == PREF_NONE) RETURN_ERR(ERRORPREFNONE);
      COV(ZERO(next), next, &v);
      if (v != 1.0) 
	SERR2("a correlation function is required as submodel, but '%.50s' has variance %10g.", NICK(next), v);
    }
    
  } else { // key != NULL
    //printf("key!=NULL\n");
   if ( (err = CHECK(key, dim,  dim, PointShapeType, XONLY, 
		     CoordinateSystemOf(OWNISO(0)),  
		     SUBMODEL_DEP, SchlatherType)) != NOERROR) {
      RETURN_ERR(err);
    }
    setbackward(cov, sub);
  }
  //print("end check schlather\n");

  RETURN_NOERROR;
}


int struct_schlather(model *cov, model **newmodel){
  model
    *sub = cov->sub[cov->sub[MPP_TCF] != NULL ? MPP_TCF : MPP_SHAPE];
  int err,  ErrNoInit;
  bool schlather = DefList[COVNR].Init == init_mpp; // else opitz
  //  assert(isnowSchlather(cov));
  ASSERT_NEWMODEL_NULL;
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  if (cov->sub[MPP_TCF] != NULL) {
    if ((err = STRUCT(sub, &(cov->key))) >  NOERROR) RETURN_ERR(err);
    SET_CALLING(cov->key, cov);
  } else {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) RETURN_ERR(err);
  }
  assert(cov->key->calling == cov);
 
  if (MODELNR(cov->key) != GAUSSPROC && !equalsBernoulliProcess(cov->key) &&
      MODELNR(cov->key) != BRNORMED
      ) {
    if (isnowVariogram(cov->key)) addModel(&(cov->key), GAUSSPROC); 
    else {
      if (isGaussMethod(cov->key)) {
	SERR("invalid model specification");
      } else SERR2("'%.50s' currently only allowed for gaussian processes %.50s", 
		   NICK(cov), schlather ? "and binary gaussian processes" : "");
    }
  }
  assert(cov->key->calling == cov);
  
  Types frame = SchlatherType;
  //    MODELNR(cov->key) == GAUSSPROC ? GaussMethodType   //Marco, 29.5.13
  //  : equalsBernoulliProcess(cov->key) ? NormedProcessType // egal
  //  : SchlatherType;

  
  // if ((err = CHECK(cov->key, cov->tsdim, cov->xdimown, ProcessType, 
  //		     OWNDOM(0),
  //		     OWNISO(0), cov->vdim, frame)) != NOERROR) RETURN_ERR(err);  

   if ((err = CHECK_PASSTF(cov->key, ProcessType, VDIM0, frame)) != NOERROR)
    RETURN_ERR(err);  
 
  if ((ErrNoInit = STRUCT(cov->key, NULL)) > NOERROR)//err
    return ErrNoInit; 
  
  addModel(&(cov->key), STATIONARY_SHAPE);
  
  // if ((err = CHECK(cov->key, cov->tsdim, cov->xdimown, PointShapeType, 
  //		     OWNDOM(0),
  //		     OWNISO(0), cov->vdim, SchlatherType)) != NOERROR)
  //RETURN_ERR(err);

   if ((err = CHECK_PASSTF(cov->key, PointShapeType, VDIM0,
			   SchlatherType)) != NOERROR) RETURN_ERR(err);
  
  return ErrNoInit;
}

#define OPITZ_ALPHA (LAST_MAXSTABLE + 1)
int init_opitzprocess(model *cov, gen_storage *S) {
  int err;
  if ( (err = init_mpp(cov, S)) != NOERROR) RETURN_ERR(err);
  
  double alpha = P0(OPITZ_ALPHA);
  model *key = cov->key;
  assert(key != NULL);
  assert(key->mpp.moments == 1);
  pgs_storage *pgs = key->Spgs;
  assert(pgs != NULL);
  key->mpp.mMplus[1] = INVSQRTTWOPI *POW(2, 0.5 * alpha - 0.5) * 
    gammafn(0.5 * alpha + 0.5); 
  pgs->zhou_c = 1.0 / key->mpp.mMplus[1];
  pgs->alpha = alpha;

  cov->simu.active = cov->initialised = true;
  RETURN_NOERROR;
}


void range_opitz(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
    //   assert(cov != NULL);
  range_mpp(cov, range);

  range->min[OPITZ_ALPHA] = 0;
  range->max[OPITZ_ALPHA] = RF_INF;
  range->pmin[OPITZ_ALPHA] = 0.1;
  range->pmax[OPITZ_ALPHA] = 5;
  range->openmin[OPITZ_ALPHA] = true;
  range->openmax[OPITZ_ALPHA] = true; 
}



////////////////////////////////////////////////////////////////////
// Smith

// siegburg 7:47, 7:54; 8:19, Bus 611 8:32; 8:41 Immenburg



void param_set_identical(model *to, model *from, int depth) {
  assert(depth >= 0);
  assert(to != NULL && from != NULL);
  assert(STRCMP(NICK(from), NICK(to)) == 0); // minimal check

  int i;

  assert(from->qlen == to->qlen);
  if (from->q != NULL) MEMCOPY(to->q, from->q, sizeof(double) * from->qlen);

  for (i=0; i<MAXPARAM; i++) { PCOPY(to, from, i); }

  if (depth > 0) {
    for (i=0; i<MAXSUB; i++) {
      if (from->sub[i] != NULL) 
	param_set_identical(to->sub[i], from->sub[i], depth - 1);
    }
  }
}

int FillInPts(model *cov,  // struct of shape
			model *shape // shape itself
			) {
  // fuegt im max-stabilen Fall diegleiche Funktion als
  // intensitaetsfunktion fuer die Locationen ein; random scale
  // wird mitbeachtet
  //
  // fuegt im coin fall die ge-powerte Shape-Fkt als
  // intensitaetsfunktion fuer die Locationen ein; random scale
  // wird mitbeachtet
  //

  // smith, randomcoin

  assert(cov != NULL);

  ASSERT_ONESYSTEM;

  //  printf("fx \n");
  int err, i,
    dim = XDIM(SYSOF(shape), 0),
    nr = COVNR;

  
  if (cov->sub[PGS_LOC] != NULL) RETURN_NOERROR;
  if ((err = covcpy(cov->sub + PGS_FCT, shape)) != NOERROR) RETURN_ERR(err);
 
  if (nr == ZHOU) {   
    //printf("here %ld %d %d %d\n", cov->sub[PGS_LOC], ScaleOnly(shape),
    //	   !shape->randomkappa, shape->sub[0]->randomkappa);
    cov->nsub = 2;
    if (cov->sub[PGS_LOC] != NULL) BUG;
    bool randomscale = ScaleOnly(shape) && shape->randomkappa &&
      !shape->sub[0]->randomkappa;
    
    if ((err = covcpyWithoutRandomParam(cov->sub + PGS_LOC, 
					randomscale ? shape->sub[0] : shape))
	  != NOERROR) RETURN_ERR(err);
    if (hasPoissonGaussFrame(shape) //|| hasMaxStableFrame(shape)
	) {
      assert(dim > 0 && dim <= MAXMPPDIM);
      addModel(cov, PGS_LOC, POW);
      assert(cov->nsub == 2);
      kdefault(cov->sub[PGS_LOC], POW_ALPHA, GLOBAL.mpp.shape_power);
      addModel(cov, PGS_LOC, SCATTER);
      model *scatter = cov->sub[PGS_LOC];

      // PMI0(scatter); printf("sm=%d %d %d\n", SCATTER_MAX, dim,SCATTER_STEP);
 
      PARAMALLOC(scatter, SCATTER_MAX, dim, 1);
      for (i=0; i<dim; i++) 
	PARAMINT(scatter, SCATTER_MAX)[i] = GLOBAL.mpp.scatter_max[i];
      PARAMALLOC(scatter, SCATTER_STEP, dim, 1);
      for (i=0; i<dim; i++) 
	PARAM(scatter, SCATTER_STEP)[i] = GLOBAL.mpp.scatter_step[i];
      BUG;
      addModel(cov, PGS_FCT, RANDOMSIGN); 
    } else if (!hasSmithFrame(shape))  BUG; 
    
    if (!randomscale && shape->randomkappa) {	
      addSetDistr(cov->sub + PGS_LOC, cov->sub[PGS_FCT], 
		  param_set_identical, true, MAXINT);
    }
    
    addModel(cov, PGS_LOC, RECTANGULAR);

    if (randomscale) {
      addModel(cov, PGS_LOC, LOC);
      addSetDistr(cov->sub + PGS_LOC, shape, ScaleDollarToLoc, true, 0);
    } 
  } else if (nr == STANDARD_SHAPE) {
    assert(cov != NULL && shape != NULL);
    if ((err = STRUCT(shape,
		      cov->sub + PGS_LOC)) != NOERROR) RETURN_ERR(err);
    SET_CALLING(cov->sub[PGS_LOC], cov);
    if (cov->sub[PGS_LOC] == NULL)
      SERR("simple enlarging of the simulation window does not work");
  } else if (nr == BALLANI) {
    RETURN_ERR(ERRORNOTPROGRAMMEDYET);
  } else BUG;

  RETURN_NOERROR;
}

int addPGS(model **Key, // struct of shape
		  model *shape,// shape itself
		  model *pts, 
		  int dim, int vdim, Types frame) {
  // SEE ALSO addPGSLocal in RMS.cc !!!!!!
  /// random coin und smith

  // versucht die automatische Anpassung einer PointShapeType-Funktion;
  // derzeit wird 
  // * PGS und
  // * STANDARD_SHAPE (weiss nicht mehr wieso -> coins?)
  // versucht; inklusive Init
  //
  // ruft t.w.  FillInPts auf
  
#define specific (nPOISSON_SCATTER - 1)
  bool maxstable = hasMaxStableFrame(shape);
  int i,
    err = NOERROR,
    method = GLOBAL.extreme.scatter_method,
    pgs[specific] = {maxstable ? ZHOU : BALLANI, STANDARD_SHAPE}; 
 #define infoN 200
  char info[nPOISSON_SCATTER - 1][LENERRMSG];
  assert(shape != NULL);
  assert(*Key == NULL);

  // to do: strokorbball: raeumlich waere Standard/Uniform besser;
  // praeferenzen programmieren?
  for (i=0; i<specific; i++) {
    if (method != POISSON_SCATTER_ANY && i != method) continue;
    if (i > 0) {
      errorMSG(err, info[i-1]);
      //     XERR(err);  // eigentlich muss das hier weg
    }
    //    if (i > 0) XERR(err); assert(i ==0);
    if (*Key != NULL) COV_DELETE(Key);
    assert(shape->calling != NULL);
    addModel(Key, pgs[i], shape->calling);
    
    if (pts == NULL) {
      if ((err = FillInPts(*Key, shape)) != NOERROR) continue;
    } else {
      if ((err = covcpy((*Key)->sub + PGS_FCT, shape)) != NOERROR ||
	  (err = covcpy((*Key)->sub + PGS_LOC, pts)) != NOERROR
 	  ) {
	model *cov = *Key;
	RETURN_ERR(err);
      }
      // wechselseitiger kreuzverweis 
      Ssetcpy((*Key)->sub[PGS_FCT], (*Key)->sub[PGS_LOC], shape, pts);
      Ssetcpy((*Key)->sub[PGS_LOC], (*Key)->sub[PGS_FCT], pts, shape);
    }
 
    model *cov = *Key;
    SET_CALLING(cov, shape->calling);
    SET_CALLING(cov->sub[PGS_FCT], cov);
    SET_CALLING(cov->sub[PGS_LOC], cov);

    assert(cov->sub[PGS_LOC] != NULL && cov->sub[PGS_FCT] != NULL);
    cov->nsub = 2;

    //    printf(">>>>>>> KEY() start %d %.50s \n", i, NICK((*Key)));
    if ((err = CHECK(*Key, dim, dim, PointShapeType, XONLY, 
		     CoordinateSystemOf(ISO(SYSOF(shape), 0)),
		     vdim, frame)) != NOERROR) {
      //        printf(">>>>>>> KEY() break %d %.50s %d\n", i, Nick(*Key), dim);  
    
      //  XERR(err);
      continue; 
    }
    NEW_COV_STORAGE(cov, gen);

    if ((err = INIT(cov, 1, cov->Sgen)) == NOERROR) break;
  } // for i_pgs

  model *cov = *Key;
  if (err != NOERROR) {
    SERR("error occured when creating the point-shape fctn");
    //   errorstring_type xx = "error occured when creating the point-shape fctn";
    // errorMSG(err, xx, cov->base, info[i-1]);
    // SERR2("no method found to simulate the given model:\n\tpgs: %.50s\n\tstandard: %.50s)", info[0], info[1]);
  }
  RETURN_NOERROR;
}


int struct_ppp_pts(model **Key, // = *cov->key
		   model *shape, 
		   model *calling, // = process
		   int timespacedim, int vdim, Types frame) {
  // MAIN FUNCTION FOR CREATING THE POINT-SHAP_OBJECT FOR
  // SMITH AND RANDOM COIN
  // Ruft S TRUCT(shape, Key) auf und verarbeitet weiter, je nachdem
  // was nun in Key steht: pgs, shape, random, nichts.
  model
    *dummy = NULL;
  int err = NOERROR,
    logdim = LOGDIM(PREVSYSOF(shape), 0),
    xdim = XDIM(PREVSYSOF(shape), 0);


#ifdef RANDOMFIELDS_DEBUGGING
  model *cov = shape;
#endif
  err = STRUCT(shape, Key); // get the random location corresponding to shape;
  //                            either zhou paper or ballani

  if (err == NOERROR && *Key != NULL) {
    // indeed there is a corresponding random location function
    // could be both smith model and random coin

    SET_CALLING(*Key, calling);

    Types type = BadType;
    if (isPointShape(*Key)) type = PointShapeType;
    else if ((err = CHECK_R(*Key, logdim)) == NOERROR) type = RandomType;

    if (!equalsRandom(type)) {
      if ((err = CHECK(*Key, logdim, xdim, type, 
		       DOM(PREVSYSOF(shape), 0), ISO(PREVSYSOF(shape), 0),
		       shape->vdim, frame)) != NOERROR) goto ErrorHandling;
      type = type == BadType ? ShapeType : type;
    }
 
    if (type == RandomType) { // random locations given for the shape fct;
      // so, it must be of pgs type (or standard):
      // could be both smith model and random coin
#ifdef RANDOMFIELDS_DEBUGGING
      model *cov = *Key;
#endif
      dummy = *Key;
      *Key = NULL;     
      if ((err = addPGS(Key, shape, dummy, timespacedim, vdim, frame))
	  != NOERROR) goto ErrorHandling;      
      if (*Key == NULL) BUG;
      SET_CALLING(*Key, calling);
    } else {
      model *cov = *Key; ASSERT_ONESYSTEM;
      if (type==PointShapeType) { // happens with struct_strokorbPoly
	if ((err = FillInPts(*Key, shape)) != NOERROR) goto ErrorHandling;
      } else { // (type == ShapeType)
	assert(CALLINGNR != SMITHPROC); // i.e. random coin
	// i.e. it delivers the (non-random) coin for a given
	// covariance function
	dummy = *Key;
	*Key = NULL;
 	  // suche nach geeigneten locationen
	if ((err = addPGS(Key, dummy, NULL, timespacedim, vdim, frame))
	    != NOERROR) goto ErrorHandling; 
 	//printf("e ddd\n");
      }
    } // ! randomtype
  } else { // S-TRUCT does not return anything
    int err2;
    //assert(false);
    // XERR(err);
    if (*Key != NULL) {
      SET_CALLING(*Key, calling);// make sure that the locations are not deleted
      COV_DELETE(Key);
    }
  
    if ((err2 = addPGS(Key, shape, NULL, timespacedim, vdim, frame))!=NOERROR){
      if (err == NOERROR) err = err2;
      goto ErrorHandling; 
    }
    err = NOERROR;
  }
 ErrorHandling:
  if (dummy != NULL) COV_DELETE(&dummy);
  model *cov = *Key;
  RETURN_ERR(err);
}


int check_smith(model *cov) {
  model 
    *shape = cov->sub[MPP_SHAPE],
    *TCF =  cov->sub[MPP_TCF],
    *next = shape != NULL ? shape : TCF,
    *key = cov->key;
    // *sub = (key==NULL) ? next : key;
  int err, 
    dim = ANYDIM; // taken[MAX DIM],
  

  ASSERT_ONE_SUBMODEL(cov);
  if ((err = SetGEVetc(cov)) != NOERROR) RETURN_ERR(err);
 
  if (key != NULL) {
    if ((err = CHECK(key, dim,  dim, PointShapeType, XONLY,
		     CoordinateSystemOf(OWNISO(0)), 
		       SUBMODEL_DEP, SmithType)) != NOERROR) RETURN_ERR(err);
  } else { // key == NULL
    if (next == TCF) {
      if ((err = CHECK(next, dim, dim, TcfType, XONLY, ISOTROPIC, 
			 SCALAR, SmithType)) != NOERROR) {
	RETURN_ERR(err);
      }

      if ((dim == 1 && next->rese_derivs < 1) ||
	  (dim >= 2 && dim <= 3 && next->rese_derivs < 2) ||
	  dim > 3) {

	//printf("dim = %d\n", dim);
	SERR("submodel does not have enough derivatives (programmed).");
      }

    } else { // shape -- geaendert
      Types frame = SmithType;
      //	isPointShape(sub) || is Shape(sub) ? SmithType
      //: isGaussMethod(sub) ? GaussMethodType
      //: equalsBernoulliProcess(sub) ?  NormedProcessType // egal
	    //: BadType; 
      
      if ((err = CHECK(next, dim, dim, ShapeType, XONLY, 
		       CoordinateSystemOf(OWNISO(0)),
		       SCALAR, frame)) != NOERROR) {
	RETURN_ERR(err);
      }
      if (next->full_derivs < 0) 
	SERR1("'%.50s' requires an explicit submodel.", NICK(cov));
    }
  } 
  
  setbackward(cov, next);

  RETURN_NOERROR;
}

int struct_smith(model *cov,  model **newmodel){
  model
    *tmp_shape = NULL,
    *shape = cov->sub[MPP_SHAPE],
    *tcf =  cov->sub[MPP_TCF],
    *tcf_shape = NULL,
    *sub = shape != NULL ? shape : tcf;
  location_type *loc = Loc(cov);
  int err = NOERROR,
    logdim = LOGDIM(PREVSYSOF(sub), 0),
    xdim = XDIM(PREVSYSOF(sub), 0);

  assert(hasSmithFrame(sub));
  if (loc->Time || (loc->grid && loc->caniso != NULL)) {
    TransformLoc(cov, false, GRIDEXPAND_AVOID, false);
    SetLoc2NewLoc(sub, PLoc(cov));
  }
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  ASSERT_NEWMODEL_NULL;
  
  if (tcf != NULL) {
    // to do: ausbauen 
    if ((err = covcpy(&tcf_shape, sub)) != NOERROR) goto ErrorHandling;   
    addModel(&tcf_shape, STROKORB_MONO);

    if ((err = CHECK(tcf_shape, logdim, xdim, ShapeType, 
		     DOM(PREVSYSOF(tcf), 0), ISO(PREVSYSOF(tcf), 0), tcf->vdim, 
		     SmithType)) != NOERROR) goto ErrorHandling;
   tmp_shape = tcf_shape;
   // SET_CALLING_NULL tcf_shape->calling ; ??
  } else {
    tmp_shape = shape;
  }
  // 

  if ((err = struct_ppp_pts(&(cov->key), tmp_shape, cov,
			    ANYDIM, VDIM0, SmithType))
      != NOERROR) goto ErrorHandling;

  err = NOERROR;

 ErrorHandling:
  if (tcf_shape != NULL && tmp_shape != NULL) COV_DELETE(&tmp_shape);
  RETURN_ERR(err);
}




typedef double (*logDfct)(double *data, double gamma);

double HueslerReisslogD(double *data, double gamma) {  
  double 
    g = SQRT(2.0 * gamma),
    logy2y1 =LOG(data[1] / data[0]);
  return -pnorm(0.5 * g + logy2y1 / g, 0.0, 1.0, true, false) / data[0] 
    -pnorm(0.5 * g - logy2y1 / g, 0.0, 1.0, true, false) / data[1];
}


double schlatherlogD(double *data, double gamma) {  
  double 
    sum = data[0] + data[1],
    prod = data[0] * data[1];
  return 0.5 * sum / prod * 
    (1.0 + SQRT(1.0 - 2.0 * (2.0 - gamma) * prod / (sum * sum)));
}

#define LL_AUTO 0
#define LL_FULL 1
#define LL_COMPOSITE 2
#define LL_SELECTION 3
#define LL_UNKNOWN ERR("unknown value for 'likelihood' -- please contact author");


void loglikelihoodMaxstable(double *data,
			    model *cov, // = process
			    logDfct logD, double *v) {
  // DARF JA NICHT DURCH MEHRERE CORES AUFGERUFEN WERDEN!!
  model *sub = cov;
  while (isnowProcess(sub))
    sub = sub->key != NULL ? sub->key : sub->sub[0];

  if (cov->q == NULL) {
    location_type *loc = Loc(cov);
    long len = loc->totalpoints;    
    QALLOC(len);
    if (loc->grid || loc->Time) TransformLoc(sub, false, True, false);   
  }

  ASSERT_ONESYSTEM;
  location_type *loc = Loc(cov);
  int dim = OWNTOTALXDIM;
  int len = loc->totalpoints;    
    
  if (data != NULL) {
    double 
      xi = P0(GEV_XI),
      mu = P0(GEV_MU),
      s = P0(GEV_S);    
    if (xi == 0) {
      for (int i=0; i<len; i++) cov->q[i] =EXP((data[i] - mu) / s);
    } else {
      double 
	xiinv = 1.0 / xi,
	xis = xi / s;      
      for (int i=0; i<len; i++)
	cov->q[i] =POW(1.0 + xis * (data[i] - mu), xiinv);
    }
  }

  switch(GLOBAL.fit.likelihood) {
  case LL_AUTO: case LL_COMPOSITE: {
    double z[2], gamma, gamma0, x[MAXMPPDIM + 1], y[MAXMPPDIM + 1],
      *xx = loc->x;    
    COV(ZERO(sub), sub, &gamma0);
    for (int i=0; i<len; i++, xx += dim) {
      MEMCOPY(x, xx, dim);
      x[dim] = i;
      double *yy = xx + dim;
      for (int j=i+1; j<len; j++, yy += dim) {
	MEMCOPY(y, yy, dim);
	y[dim] = j;
	NONSTATCOV(x, y, sub, &gamma);
	z[0] = cov->q[i];
	z[1] = cov->q[j];
	*v += logD(z, gamma0 - gamma);
      }
    }
  }
    break;
  case LL_FULL: ERR("full MLE not available for Brown Resnick");
    break;
  case LL_SELECTION: NotProgrammedYet("LL_SELECTION");
  default : LL_UNKNOWN;
  }
}


void loglikelihoodBR(double *data, model *cov, double *v) {
  loglikelihoodMaxstable(data, cov, HueslerReisslogD, v);
}

void loglikelihoodSchlather(double *data, model *cov, double *v) {
  loglikelihoodMaxstable(data, cov, schlatherlogD, v);
}






//////////////////////////////////////////////////////////////////////
//           Random Coin
//////////////////////////////////////////////////////////////////////



int check_randomcoin(model *cov) {
  model 
    *pdf = cov->sub[COIN_COV],
    *shape = cov->sub[COIN_SHAPE],
    *next = shape != NULL ? shape : pdf,
    *key = cov->key;
  int err,
    dim = ANYDIM; // taken[MAX DIM],
  mpp_param *gp  = &(GLOBAL.mpp);
  //extremes_param *ep = &(GLOBAL.extreme);

  SERR("'random coin' method does not work for the current version");

  INTERNAL;

  ASSERT_ONE_SUBMODEL(cov);
  if ((!hasPoissonGaussFrame(next) && !hasGaussMethodFrame(next) &&
       !hasProcessFrame(next))
      || cov->key!=NULL) ILLEGAL_FRAME;
  kdefault(cov, RANDOMCOIN_INTENSITY, gp->intensity[dim]);
  kdefault(cov, RANDOMCOIN_METHOD, 0);
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);
 
  if (key != NULL) {
    // note: key calls first function to simulate the points
    // if this is true then AverageIntern/RandomCoinIntern calls Average/RandomCoin  
    model *intern = cov;
    while (intern != NULL && MODELNR(intern) != AVERAGE_INTERN &&
	   isAnyDollar(intern))
      intern = intern->key != NULL ? intern->key : intern->sub[0];
    if (intern == NULL) {
       BUG;
    }

    
    if (hasGaussMethodFrame(next))
      err = CHECK(key, dim, dim, ProcessType, XONLY, 
		  CoordinateSystemOf(OWNISO(0)),
		  SUBMODEL_DEP, PoissonGaussType);
    else if (hasPoissonGaussFrame(next))
      err = CHECK(key, dim,  dim, PointShapeType, XONLY,
		  CoordinateSystemOf(OWNISO(0)),
		  SUBMODEL_DEP, PoissonGaussType);
    else err = ERRORFAILED;
    if (err != NOERROR) RETURN_ERR(err);

  } else { // key == NULL
    
    if (next == pdf) {      // s ymmetric: insbesondere also univariate
      if ((err = CHECK(next, dim, dim, PosDefType, XONLY,
		       SYMMETRIC, SCALAR, PoissonGaussType)) != NOERROR)  {
	RETURN_ERR(err);
      }      
      if ((pdf->pref[Average] == PREF_NONE) + 
	  (pdf->pref[RandomCoin] == PREF_NONE) !=1){

	//assert(false);
	RETURN_ERR(ERRORPREFNONE);
      }
    } else { // shape != NULL
      if ((err = CHECK(next, dim, dim, ShapeType, XONLY, 
		       CoordinateSystemOf(OWNISO(0)),
		       SCALAR, PoissonType)) != NOERROR)  { //POISS_GAUSS ??
	RETURN_ERR(err);
      }
    }    
  }
 
  setbackward(cov, key != NULL ? key : next);
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);


  RETURN_NOERROR;
}

 
void range_randomcoin(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  GAUSS_COMMON_RANGE;

  range->min[RANDOMCOIN_INTENSITY] = RF_NEGINF;
  range->max[RANDOMCOIN_INTENSITY] = RF_INF;
  range->pmin[RANDOMCOIN_INTENSITY] = - 10;
  range->pmax[RANDOMCOIN_INTENSITY] = 10;
  range->openmin[RANDOMCOIN_INTENSITY] = true;
  range->openmax[RANDOMCOIN_INTENSITY] = true; 

  range->min[RANDOMCOIN_METHOD] = 0;
  range->max[RANDOMCOIN_METHOD] = 1;
  range->pmin[RANDOMCOIN_METHOD] = 0;
  range->pmax[RANDOMCOIN_METHOD] = 1;
  range->openmin[RANDOMCOIN_METHOD] = false;
  range->openmax[RANDOMCOIN_METHOD] = false; 
}

 
int struct_randomcoin(model *cov, model **newmodel){
  model *tmp_shape = NULL,
    *pdf = cov->sub[COIN_COV],
    *shape = cov->sub[COIN_SHAPE];
  location_type *loc = Loc(cov);
  int err,
    dim = ANYDIM; // taken[MAX DIM],

  if (loc->Time || (loc->grid && loc->caniso != NULL)) {
    TransformLoc(cov, true, GRIDEXPAND_AVOID, false);
    SetLoc2NewLoc(pdf == NULL ? shape : pdf, PLoc(cov)); 
  }
  if (cov->key != NULL) COV_DELETE(&(cov->key));
  ASSERT_NEWMODEL_NULL;

 
  if (pdf != NULL) {// s ymmetric: insbesondere also univariate
    if ((err = CHECK(pdf, dim,  dim, PosDefType, XONLY,
		     SYMMETRIC, SCALAR, PoissonGaussType)) != NOERROR) {
      //sicherstellen, dass pdf von der richtigen Form, insb. frame_poisson_gauss
      RETURN_ERR(err);
    }
    if (pdf->pref[Average] == PREF_NONE && pdf->pref[RandomCoin]==PREF_NONE) {
      // if (pdf->nr > LASTDOLLAR) AERR(ERRORPREFNONE);
      RETURN_ERR(ERRORPREFNONE);
    }
    if ((err = STRUCT(pdf, &tmp_shape)) != NOERROR) 
      goto ErrorHandling; //&(cov->key)
    if (tmp_shape == NULL)
      GERR("no structural information for random coins given");

    SET_CALLING(tmp_shape, cov);

    if ((err = CHECK(tmp_shape, dim, dim, ShapeType, XONLY,
		     CoordinateSystemOf(OWNISO(0)),
		     SCALAR, PoissonGaussType)) != NOERROR) {
      goto ErrorHandling;
    }    
  } else {
    tmp_shape = shape;
    //  if ((err = covcpy(&(cov->key), shape)) > NOERROR) {
    //      RETURN_ERR(err);
    //    }
  } 
 
  // if ((err = STRUCT(cov, NULL)) != NOERROR) RETURN_ERR(err); ?????

  SERR("Sorry, 'random coin' does not work currently.");

  assert(false); // todo
  switch(P0INT(RANDOMCOIN_METHOD)) {
  case 0: {
    // Alternativ?!
    // if ((err = addPGS(&(cov->key), shape, NULL, timespacedim, vdim)) != NOERROR)
    if ((err = struct_ppp_pts(&(cov->key), tmp_shape, cov, dim, VDIM0,
			      PoissonGaussType)) != NOERROR) goto ErrorHandling;
   break;
  }
  case 1: {
 
    if ((err = covcpy(&(cov->key), shape)) != NOERROR) goto ErrorHandling;
    addModel(&(cov->key), RANDOMSIGN, cov); 
    addModel(&(cov->key), MCMC_PGS, cov);
    model *key = cov->key;
    if ((err = covcpy(key->sub + PGS_LOC, shape)) != NOERROR)
	  goto ErrorHandling;
    addModel(key, PGS_LOC, POW);
    kdefault(key->sub[PGS_LOC], POW_ALPHA, GLOBAL.mpp.shape_power);
    addModel(key, PGS_LOC, SCATTER);
    model *scatter = key->sub[PGS_LOC];
    PARAMALLOC(scatter, SCATTER_MAX, dim, 1);
    for (int i=0; i<dim; i++) PARAM(scatter, SCATTER_MAX)[i] = NA_INTEGER;
    PARAMALLOC(scatter, SCATTER_STEP, dim, 1);
    for (int i=0; i<dim; i++) PARAM(scatter, SCATTER_STEP)[i] = RF_NA;
    addModel(&(cov->key), MCMC, cov);
    break;
  }
  default: BUG;
  }
 
  if ((err = CHECK(cov->key, dim, dim, PointShapeType, XONLY,
		   CoordinateSystemOf(OWNISO(0)),
		   VDIM0, PoissonGaussType)) != NOERROR) goto ErrorHandling;

 ErrorHandling:
  if (pdf != NULL && tmp_shape != NULL) COV_DELETE(&tmp_shape);
  RETURN_ERR(err);
}


int init_randomcoin(model *cov, gen_storage *S) {
  model
    *covshape = cov->sub[ cov->sub[COIN_SHAPE] != NULL ? COIN_SHAPE : COIN_COV],
    *key = cov->key,
    *sub = key == NULL ? covshape : key; 
  char name[] = "Poisson-Gauss";
  location_type *loc = Loc(cov);
  //mpp_storage *s = &(S->mpp);  
  int 
    err = NOERROR;
  KEY_type *KT = cov->base;

  SPRINTF(KT->error_loc, "%.50s process", name);

  assert(hasPoissonGaussFrame(sub));

  cov->method = covshape->pref[Average] == PREF_NONE ? RandomCoin : Average;
    
  if (cov->method == Average && loc->caniso != NULL) {
    bool diag, quasidiag, semiseparatelast, separatelast;
    int idx[MAXMPPDIM];
    assert(loc->timespacedim <= MAXMPPDIM);
    analyse_matrix(loc->caniso, loc->cani_nrow, loc->cani_ncol, &diag,
		   &quasidiag, idx, &semiseparatelast, &separatelast);
    if (!separatelast) SERR("not a model where time is separated");
  }
 
  if ((err = init_mpp(cov, S)) != NOERROR) {
    RETURN_ERR(err);    
  }

  sub->Spgs->intensity = 
    sub->Spgs->totalmass * P0(RANDOMCOIN_INTENSITY); 
  sub->Spgs->log_density =LOG(P0(RANDOMCOIN_INTENSITY));

  assert(sub->mpp.moments >= 2);
  if (!R_FINITE(sub->Spgs->totalmass) || !R_FINITE(sub->mpp.mM[2]))
    SERR("Moments of submodels not known");

  RETURN_NOERROR;
}



void do_randomcoin(model *cov, gen_storage *s) {   
  assert(cov->simu.active);
  SAVE_GAUSS_TRAFO;
  double *res = cov->rf;
  dompp(cov, cov->Sgen==NULL ? s : cov->Sgen);// letzteres falls shape gegeben
  BOXCOX_INVERSE;
}

