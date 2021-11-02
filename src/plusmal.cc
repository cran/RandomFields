/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of auxiliary correlation functions 

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


#include "def.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "questions.h"
#include "operator.h"
#include "Processes.h"
#include "startGetNset.h"
#include "variogramAndCo.h"
#include "primitive.others.h"
#include "extern.h"






int checkplusmal(model *cov) {
  model *sub;
  int  err, 
    vdim[2] = {1, 1},
    last = OWNLASTSYSTEM;
  bool plus = DefList[COVNR].check == checkplus,
    fullXdim = hasLikelihoodFrame(cov),
    trend = equalsnowTrend(cov),
    top = cov->calling == NULL || !isnowShape(cov->calling);
  Types covtype = OWNTYPE(0),
    frame = fullXdim ? EvaluationType : cov->frame; // see also below exception:
  // LikelihoodType is passed downwards if sub is plus!
  
  /*
  domain_type covdom = trend ? XONLY : OWNDOM(0);
  int trendiso = U pgradeToCoordinateSystem(OWNISO(0));
  if (trendiso == ISO_MISMATCH) trendiso = OWNISO(0);
  assert(trendiso != ISO_MISMATCH);
  int coviso = trend ? trendiso : OWNISO(0);
  */

  bool *conform = cov->Splus == NULL ? NULL : cov->Splus->conform;
 
  int possibilities = 1 +
    (int) (!trend && top && (isNegDef(covtype) || isProcess(covtype))
	   && hasAnyEvaluationFrame(cov)
	   );

  //PMI(cov);  printf("hasEval = %d\n", hasAnyEvaluationFrame(cov));

  //  printf("top = %d %d %d shape=%d %s\n", top, cov->calling == NULL, cov->calling != NULL && !isnowShape(cov->calling), cov->calling != NULL && isnowShape(cov->calling), cov->calling == NULL ? "Niente" : NAME(cov->calling));

  cov->matrix_indep_of_x = true;
  //  printf("checking plus A\n");
  for (int i=0; i<cov->nsub; i++) { 
    Types type = covtype;
    sub = cov->sub[i];
    if (sub == NULL) 
      SERR("+ or *: named submodels are not given in a sequence!");

    if (plus) {
      if (top && equalsnowMathDef(sub) && isNegDef(type)) continue;
    } else {
      if (equalsVariogram(type)) type=PosDefType;
    }
  
    err = CERRORTYPECONSISTENCY;


    //printf("+ %s %d\n", NAME(sub), possibilities);
    
    for (int j=0; j<possibilities; j++) {// nur trend als abweichender typus erlaubt
      //
      if (false && COVNR != PLUS) {
	//printf("mal j = %d %.50s %.50s %d trend=%d top=%d \n", j, TYPE_NAMES[type], NAME(sub), possibilities, trend, top);
	//PMI0(cov); PMI0(sub);
      }
      COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);
      for (int s=0; s<=last; s++) {	
	if (j != 0) set_type(PREVSYSOF(sub), s, TrendType);
	if (equalsTrend(type)) {
	  //if (!equalsCoordinateSystem(OWNISO(s))) RETURN_ERR(ERRORWRONGISO);
	  set_dom(PREVSYSOF(sub), s, XONLY);
 	} else if (false && isNegDef(type)) {
	  if (!isAnyIsotropic(PREVISO(s)) && !isSpaceIsotropic(PREVISO(s))
	      && !equalsKernel(PREVDOM(s))) RETURN_ERR(ERRORWRONGISO);
	} else {
	  set_dom(PREVSYSOF(sub), s, OWNDOM(s));
	}
	if (fullXdim && equalsTrend(type) && !equalsCoordinateSystem(OWNISO(s)))
	  RETURN_ERR(ERRORWRONGISO);
	set_iso(PREVSYSOF(sub), s, OWNISO(s));
      }

      //     A     PMI(cov);

      // incorrect :

      //   TREE0(cov);
      //    printf("!!plus %d %.50s,%d=%d %d\n", plus, NAME(sub), SUBNR,  PLUS, fullXdim);
      
      if ((err = CHECK_GEN(sub, SUBMODEL_DEP, SUBMODEL_DEP,
	        	   plus && (SUBNR == PLUS) && fullXdim ? LikelihoodType
	          : type == TrendType ? TrendType : frame, true))
	  == NOERROR) break;

      // 26.2.19
      if ((!isNegDef(type) && !isProcess(type)) || isTrend(type)) break;
 
      //APMI(cov);
      //BUG;
      
      type = TrendType;     
    }

    // printf("mult err = %d %d\n", err, cov->err);

    
    if (err != NOERROR) RETURN_ERR(err);

    if (i == 0) {
      setbackward(cov, sub);
      vdim[0] = cov->vdim[0];
      vdim[1] = cov->vdim[1];      
    } else {      
      updatepref(cov, sub);
      //      VDIM0 = sub->vdim[0];
      // VDIM1 = sub->vdim[1];
      cov->randomkappa |= sub->randomkappa;

      bool vdim0ok = vdim[0] == cov->vdim[0];
      if ( !vdim0ok||
	   (vdim[1] != cov->vdim[1] && vdim[1]!=1 && cov->vdim[1] != 1 &&
	    !equalsnowTrend(sub) && !equalsnowTrend(cov->sub[0]))) {
	//	printf("i=%d\n", i);	APMI(cov);
	SERR4("multivariate dimensionality is different in the submodels (%.50s is %d-variate; %.50s is %d-variate)", NICK(cov->sub[0]), cov->vdim[vdim0ok], NICK(sub), vdim[vdim0ok]);
      }
    }
    
    if (false)
      for(int j=0; j<2; j++) {
      if (vdim[j] == 1) {
	if (cov->vdim[j] != 1) vdim[j] = cov->vdim[j];
      } else {
	if (cov->vdim[j] != 1 && cov->vdim[j] != vdim[j])
	  SERR4("multivariate dimensionality is different in the submodels (%.50s is %d-variate; %.50s is %d-variate)", NICK(cov->sub[0]), cov->vdim[j], NICK(sub), vdim[j]);
      }
    }
    /*
    if (i == 0) vdim[0] = VDIM0;
    if (vdim[1] == 1 && VDIM1 != 1) vdim[1] = VDIM1;
    if (i > 0 && (vdim[0] != VDIM0 ||
		  (VDIM1 != 1 && VDIM1 != vdim[1])))
      SERR4("multivariate dimensionality is different in the submodels (%.50s is %d-variate; %.50s is %d-variate)", NICK(cov->sub[0]), cov->sub[0]->vdim[0],
	    NICK(sub), sub->vdim[0]);
    */		  

    cov->matrix_indep_of_x &= sub->matrix_indep_of_x;
    if (conform != NULL)
      conform[i] = !isBad(TypeConsistency(covtype, SUBTYPE(0)));
  } // i, nsub

  VDIM0 = vdim[0];
  VDIM1 = vdim[1];

  // !! incorrect  !!
  // cov->semiseparatelast = false; 
  //cov->separatelast = false;

  //  printf("plus checked\n");

  RETURN_NOERROR;
}





// see private/old.code/operator.cc for some code including variable locations
void select(double *x, model *cov, double *v) {
  int len,
    *element = PINT(SELECT_SUBNR);
  model *sub = cov->sub[*element];
  assert(VDIM0 == VDIM1);
  if (*element >= cov->nsub) ERR("select: element out of range");
  COV(x, sub, v);
  if ( (len = cov->nrow[SELECT_SUBNR]) > 1) {
    int i, m,
      vdim = VDIM0,
      vsq = vdim * vdim;
    TALLOC_XX1(z, vsq);
    for (i=1; i<len; i++) {
      sub = cov->sub[element[i]];
      COV(x, sub, z);
      for (m=0; m<vsq; m++) v[m] += z[m]; 
    }
    END_TALLOC_XX1;
 }
}
  

void covmatrix_select(model *cov, double *v) {
  int len = cov->nrow[SELECT_SUBNR];
  
  if (len == 1) {
    int element = P0INT(SELECT_SUBNR);
    model *next = cov->sub[element];
    if (element >= cov->nsub) ERR("select: element out of range");
    DefList[NEXTNR].covmatrix(next, v);
  }  else StandardCovMatrix(cov, v);
}

char iscovmatrix_select(model VARIABLE_IS_NOT_USED *cov) {  return 2; }

int checkselect(model *cov) {
  int err;

  assert(cov->Splus == NULL);

  if (!isCartesian(OWNISO(0))) BUG; // RETURN_ERR(ERRORNOTCARTESIAN);
  kdefault(cov, SELECT_SUBNR, 0);
  if ((err = checkplus(cov)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);

  EXTRA_STORAGE;
  RETURN_NOERROR;
}

bool allowedDselect(model *cov) { return allowedDplus(cov); }
bool allowedIselect(model *cov) { return allowedIplus(cov); }


void rangeselect(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  range->min[SELECT_SUBNR] = 0;
  range->max[SELECT_SUBNR] = MAXSUB-1;
  range->pmin[SELECT_SUBNR] = 0;
  range->pmax[SELECT_SUBNR] = MAXSUB-1;
  range->openmin[SELECT_SUBNR] = false;
  range->openmax[SELECT_SUBNR] = false;
}



void plusStat(double *x, model *cov, double *v){
  model *sub;
  int i, m,
    nsub=cov->nsub,
    vdim = VDIM0,
    vsq = vdim * vdim;
  bool *conform = cov->Splus->conform;
  TALLOC_XX1(z, vsq);

  for (m=0; m<vsq; m++) v[m] = 0.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];    

    if (conform[i]) {
      //printf("i=%d\n", i);
      //PMI(cov);
      
      FCTN(x, sub, z);

      //     TREE0(cov);
      //     PMI(cov);
      //  printf("i=%d m=%d %.50s %10g\n", i,  m, NAME(sub), z[0]);
      
      if (sub->vdim[0] == 1) for (m=0; m<vsq; m++) v[m] += z[0]; 
      else for (m=0; m<vsq; m++) v[m] += z[m];
      //printf("%% %d\n", i);
    }
    // if (!R_FINITE(v[0]) || !R_FINITE(v[1]) || !R_FINITE(v[2]) || !R_FINITE(v[3])) { PMI(sub); printf("\n\nplus i=%d m=%d x=%10g %10g %10g\n", i, m, x[0], v[0], z[0]);    printf("(%4.3f, %4.3f; %4.3e %4.3e %4.3e %4.3e)\t", x[0], x[1], v[0], v[1], v[2], v[3]);BUG; }
    //APMI(cov->calling);
   // 
  }
  END_TALLOC_XX1;
}

void plusNonStat(double *x, double *y, model *cov, double *v){
  model *sub;
  int i, m,
    nsub=cov->nsub,
    vsq = VDIM0 * VDIM1;
  bool *conform = cov->Splus->conform;
  assert(VDIM0 == VDIM1);
  TALLOC_XX1(z, vsq);
  for (m=0; m<vsq; m++) v[m] = 0.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    if (conform[i]) {
      //   if (Type Consistency(cov->typus, sub->typus)) {
      NONSTATCOV(x, y, sub, z);
      if (sub->vdim[0] == 1) for (m=0; m<vsq; m++) v[m] += z[0]; 
      else for (m=0; m<vsq; m++) v[m] += z[m]; 
    }
  }
  END_TALLOC_XX1;
}

void Dplus(double *x, model *cov, double *v){
  model *sub;
  int n = cov->nsub, i,
    vsq = VDIM0 * VDIM1;
  bool *conform = cov->Splus->conform;
  TALLOC_XX1(z, vsq);
  for (int m=0; m<vsq; m++) v[m] = 0.0;
  for (i=0; i<n; i++) { 
    sub = cov->sub[i];
    if (conform[i]) {
      Abl1(x, sub, z);
      if (sub->vdim[0] == 1) for (int m=0; m<vsq; m++) v[m] += z[0]; 
      else for (int m=0; m<vsq; m++) v[m] += z[m]; 
    }
  }
  END_TALLOC_XX1;
}

void DDplus(double *x, model *cov, double *v){
  model *sub;
  int n = cov->nsub, i,
    vsq = VDIM0 * VDIM1;
  bool *conform = cov->Splus->conform;
  TALLOC_XX1(z, vsq);
  for (int m=0; m<vsq; m++) v[m] = 0.0;
  for (i=0; i<n; i++) { 
    sub = cov->sub[i];
   if (conform[i]) {
     Abl2(x, sub, z);
      if (sub->vdim[0] == 1) for (int m=0; m<vsq; m++) v[m] += z[0]; 
      else for (int m=0; m<vsq; m++) v[m] += z[m]; 
    }
  }
  END_TALLOC_XX1;
}


int checkplus(model *cov) {
  int err, i;
  ONCE_NEW_STORAGE(plus);
  bool *conform = cov->Splus->conform;
  if ((err = checkplusmal(cov)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  if (OWNDOM(0) == DOMAIN_MISMATCH) RETURN_ERR(ERRORNOVARIOGRAM);
  if (cov->nsub == 0) cov->pref[SpectralTBM] = PREF_NONE;

  if (isnowPosDef(cov) && isXonly(OWN)) cov->logspeed = 0.0;
  else if (isnowVariogram(cov) && isXonly(OWN)) {
   cov->logspeed = 0.0;
    for (i=0; i<cov->nsub; i++) {
      model *sub = cov->sub[i];
      if (conform[i]) {
	if (ISNAN(sub->logspeed)) {
	  cov->logspeed = RF_NA;
	  break;
	} else cov->logspeed += sub->logspeed;
      }
    }
  } else cov->logspeed = RF_NA;

  EXTRA_STORAGE;
  RETURN_NOERROR;

  // spectral mit "+" funktioniert, falls alle varianzen gleich sind,
  // d.h. nachfolgend DOLLAR die Varianzen gleich sind oder DOLLAR nicht
  // auftritt; dann zufaellige Auswahl unter "+"
}


bool allowedDplus(model *cov) {
  //assert(COVNR == MULT);
  model **sub = cov->Splus != NULL && cov->Splus->keys_given
    ? cov->Splus->keys : cov->sub;
 
  bool *D = cov->allowedD;
  int nsub = MAXSUB, // cov->nsub,
    j=0,
    idx = (int) FIRST_DOMAIN;

  while (true) {
    while (j<nsub && sub[j] == NULL) j++;
    if (j >= nsub) return allowedDtrue(cov);
    if (!allowedD(sub[j])) break;
    j++;
  }
  MEMCOPY(D, sub[j]->allowedD, sizeof(allowedD_type));

  while (idx <= (int) LAST_DOMAINUSER && !D[idx]) idx++;
  assert(idx <= (int) LAST_DOMAINUSER);
  if (idx == (int) LAST_DOMAINUSER) return false;
  for (j++; j<nsub; j++) {
    if (sub[j] == NULL) continue;

    if (allowedD(sub[j])) continue;
    bool *subD = sub[j]->allowedD;
    int subDidx = FIRST_DOMAIN;
    while (subDidx <= (int) LAST_DOMAINUSER && !subD[subDidx]) subDidx++;
    assert(subDidx <= (int) LAST_DOMAINUSER);
    for (; idx < subDidx; D[idx++] = false);
    for (int i = idx; i<=(int) LAST_DOMAINUSER; i++) D[i] |= subD[i];
    if (idx == (int) LAST_DOMAINUSER) return false;
  }
  return false;
}


bool allowedIplus(model *cov) {
  bool *I = cov->allowedI;
  if (COVNR == PLUS && hasLikelihoodFrame(cov)) {
    for(int i=FIRST_ISOUSER; i<= LAST_ISOUSER; i++) I[i] = false;
    I[CARTESIAN_COORD] = I[SPHERICAL_COORD] = I[EARTH_COORD] = true;
    return false;
  }
  //ssert(COVNR == MULT);
  model **Sub = cov->Splus != NULL  && cov->Splus->keys_given
    ? cov->Splus->keys : cov->sub;
  //  bool *I = cov->allowedI;
  int nsub = cov->nsub,
    z = 0;

  model *sub[MAXSUB];
  for (int i=0; z<nsub; i++) if (Sub[i]!=NULL) sub[z++] = Sub[i];
  assert(z == nsub);
  bool allowed = allowedIsubs(cov, sub, z);
  if (COVNR == PLUS) {
    I[CARTESIAN_COORD] = I[SPHERICAL_COORD] = I[EARTH_COORD] = true;
  }
  return allowed;
}


//int Zaehler = 0;

Types Typeplus(Types required, model *cov, isotropy_type required_iso) {
  //  assert(false);
  bool allowed = isShape(required) || isTrend(required) ||
    equalsRandom(required);
    //||required==ProcessType ||required==GaussMethodType; not yet allowed;to do
 
  if (!allowed) return BadType;

  // printf("here %.50s %d\n", NAME(cov), ++Zaehler);
 
  if (isManifold(required)) { // nue 6.8.17
    BUG;
    int last = PREVLASTSYSTEM;
    for (int s=1; s<=last; s++)
      if (!isSameAsPrev(PREVTYPE(s))) {
	// to do : not programmed yet
	//	printf("FAILED! %.50s %d\n", NAME(cov), Zaehler);
	return BadType;
      }
    
    //   printf("RETURNING %.50s %d\n", NAME(cov), Zaehler);
    return TypeConsistency(PREVTYPE(0), cov, required_iso);
  }
  
  for (int i=0; i<cov->nsub; i++) {
    //   print("i=%d %.50s\n", i, NAME(cov->sub[i]));
    if (!isBad(TypeConsistency(required, cov->sub[i], required_iso)))
      return required;
    //   print("i not ok\n");
  }
  
  // printf("FAILED %.50s %d\n", NAME(cov), Zaehler);
  return BadType;
}

void spectralplus(model *cov, gen_storage *s, double *e){
  assert(VDIM0 == 1);
  int nr;
  double dummy;
  spec_properties *cs = &(s->spec);
  double *sd_cum = cs->sub_sd_cum;

  nr = cov->nsub - 1;
  dummy = UNIFORM_RANDOM * sd_cum[nr];
  if (ISNAN(dummy)) BUG;
  while (nr > 0 && dummy <= sd_cum[nr - 1]) nr--;
  model *sub = cov->sub[nr];
  SPECTRAL(sub, s, e);  // nicht gatternr
}


int structplus(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  int m, err;  
  switch(cov->frame) {
  case EvaluationType :  RETURN_NOERROR; // VariogramType vor 11.1.19
  case GaussMethodType : 
    if (isnowProcess(cov)) {
      assert(COVNR != PLUS_PROC);
      assert(COVNR == PLUS);
      //COVNR = PLUS_PROC;
      BUG;
      //return structplusproc(cov, newmodel); // kein S-TRUCT !!
    }
    if (cov->Splus != NULL && cov->Splus->keys_given) BUG;
    for (m=0; m<cov->nsub; m++) {
      model *sub = cov->sub[m];
      if ((err = STRUCT(sub, newmodel))  > NOERROR) {
	RETURN_ERR(err);
      }
    }
    break;
  default :
    SERR2("frame '%.50s' not allowed for '%.50s'", TYPE_NAMES[cov->frame],
	  NICK(cov));
  }
  RETURN_NOERROR;
}


int initplus(model *cov, gen_storage *s){
  int i, err,
    vdim = VDIM0;
  if (VDIM0 != VDIM1) BUG; // ??

  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = RF_NA;

  if (hasGaussMethodFrame(cov)) {
    spec_properties *cs = &(s->spec);
    double *sd_cum = cs->sub_sd_cum;
 
    if (VDIM0 == 1) {
      for (i=0; i<cov->nsub; i++) {
	model *sub = cov->Splus == NULL || !cov->Splus->keys_given
	  ? cov->sub[i] : cov->Splus->keys[i];
	
	assert(sub != NULL);
	if (sub->pref[Nothing] > PREF_NONE) { // to do ??
	  // for spectral plus only
	  COV(ZERO(sub), sub, sd_cum + i);
	  if (i>0) sd_cum[i] += sd_cum[i-1];
	}
	cov->sub[i]->Sgen = (gen_storage *) MALLOC(sizeof(gen_storage));
	if ((err = INIT(sub, cov->mpp.moments, s)) != NOERROR) {
	  //  AERR(err);
	  RETURN_ERR(err);
	}
	sub->simu.active = true;
      }
    } 
 
    cov->fieldreturn = (ext_bool) (cov->Splus!=NULL && cov->Splus->keys_given);
    cov->origrf = false;
    if (cov->fieldreturn) cov->rf = cov->Splus->keys[0]->rf;
     
    RETURN_NOERROR;
  }

  else if (hasAnyEvaluationFrame(cov)) {
    RETURN_NOERROR;
  }

  RETURN_ERR(ERRORFAILED);
}


void doplus(model *cov, gen_storage *s) {
  int i;
  
  if (hasGaussMethodFrame(cov) && cov->method==SpectralTBM) {
    ERR("error in doplus with spectral");
  }
  
  for (i=0; i<cov->nsub; i++) {
    model *sub = cov->Splus!=NULL && cov->Splus->keys_given
      ? cov->Splus->keys[i] : cov->sub[i];
    DO(sub, s);
  }
}




void covmatrix_plus(model *cov, double *v) {
  location_type *loc = Loc(cov);
  //  defn *C = DefList + COVNR; // nicht gatternr
  long i, 
    totalpoints = loc->totalpoints,
    vdimtot = totalpoints * VDIM0,
    vdimtotSq = vdimtot * vdimtot;
  int
    nsub = cov->nsub;
  bool is = iscovmatrix_plus(cov) >= 2;

  if (is) {
    TALLOC_X1(mem, vdimtotSq);// da sowieso hoffnungslos
    // Bei der Berechnung der einzelnen Kovarianzmatrizen, muss darauf
    // geachtet werden, dass u.U. Koord-Trafos stattfinden. Somit
    // wird selectnr aufgerufen
    if (mem != NULL) {      
      int j;
      if (PisNULL(SELECT_SUBNR)) PALLOC(SELECT_SUBNR, 1, 1); // alloc scalar
      P(SELECT_SUBNR)[0] = 0;
      DefList[SELECTNR].covmatrix(cov, v);
      for (i=1; i<nsub; i++) {
	if (Loc(cov->sub[i])->totalpoints != totalpoints) BUG;
	P(SELECT_SUBNR)[0] = i;
	DefList[SELECTNR].covmatrix(cov, mem);
	for (j=0; j<vdimtotSq; j++) v[j] += mem[j];
      }
      END_TALLOC_X1;
      return;
    }
    END_TALLOC_X1;
  }
  StandardCovMatrix(cov, v);
}

char iscovmatrix_plus(model *cov) {
  char max=2, is; // 0
  int i,
    nsub = cov->nsub;
  for (i=0; i<nsub; i++) {
    model *sub = cov->sub[i];
    is = DefList[SUBNR].is_covmatrix(sub);
    if (is < max) max = is; // > 
  }
  return max;
}


void malStat(double *x, model *cov, double *v){
  model *sub;
  int i, m,
    nsub=cov->nsub,
    vdim = VDIM0,
    vsq = vdim * vdim;
  TALLOC_XX1(z, vsq);
  
  //  assert(x[0] >= 0.0 || OWNXDIM(0) > 1);
  for (m=0; m<vsq; m++) v[m] = 1.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    COV(x, sub, z);
    if (sub->vdim[0] == 1) for (m=0; m<vsq; m++) v[m] *= z[0];
    else for (m=0; m<vsq; m++) v[m] *= z[m];  
  }
  END_TALLOC_XX1;
}

void logmalStat(double *x, model *cov, double *v, double *Sign){
  model *sub;
  int i, m,
    nsub=cov->nsub,
    vdim = VDIM0,
    vsq = vdim * vdim;
  TALLOC_XX1(z, vsq);
  TALLOC_XX2(zSign, vsq);

  assert(VDIM0 == VDIM1); 
  assert(x[0] >= 0.0 || OWNTOTALXDIM > 1);
  for (m=0; m<vsq; m++) {v[m] = 0.0; Sign[m]=1.0;}
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    LOGCOV(x, sub, z, zSign);
    if (sub->vdim[0] == 1) {
      for (m=0; m<vsq; m++) {
	v[m] += z[0]; 
	Sign[m] *= zSign[0];
      }
    } else {
      for (m=0; m<vsq; m++) {
	v[m] += z[m]; 
	Sign[m] *= zSign[m];
      }
    }
  }
  END_TALLOC_XX1;
  END_TALLOC_XX2;
}

void malNonStat(double *x, double *y, model *cov, double *v){
  model *sub;
  int i, m, nsub=cov->nsub,
    vdim = VDIM0,
    vsq = vdim * vdim;
  TALLOC_XX1(z, vsq);

  assert(VDIM0 == VDIM1);
  assert(VDIM0 == VDIM1);

  for (m=0; m<vsq; m++) v[m] = 1.0;
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    NONSTATCOV(x, y, sub, z);
    if (sub->vdim[0] == 1) for (m=0; m<vsq; m++) v[m] *= z[0];
    else for (m=0; m<vsq; m++) v[m] *= z[m];  
  }
  END_TALLOC_XX1;
}

void logmalNonStat(double *x, double *y, model *cov, double *v, 
		   double *Sign){
  model *sub;
  int i, m, nsub=cov->nsub,
    vdim = VDIM0,
    vsq = vdim * vdim;
  assert(VDIM0 == VDIM1);
  TALLOC_XX1(z, vsq);
  TALLOC_XX2(zSign, vsq);
  for (m=0; m<vsq; m++) {v[m] = 0.0; Sign[m]=1.0;}
  for(i=0; i<nsub; i++) {
    sub = cov->sub[i];
    LOGNONSTATCOV(x, y, sub, z, zSign);
    if (sub->vdim[0] == 1) {
      for (m=0; m<vsq; m++) {
	v[m] += z[0]; 
	Sign[m] *= zSign[0];
      }
    } else {
      for (m=0; m<vsq; m++) {
	v[m] += z[m]; 
	Sign[m] *= zSign[m];
      }
    }
  }
  END_TALLOC_XX1;
  END_TALLOC_XX2;
}

void Dmal(double *x, model *cov, double *v){
  model *sub;
  int n = cov->nsub, i,
    vsq = VDIM0 * VDIM1;
  TALLOC_XXX1(c, vsq * n);
  TALLOC_XXX2(d, vsq * n);

  for (i=0; i<n; i++) {
    sub = cov->sub[i];
    COV(x, sub, c + i * vsq);
    Abl1(x, sub, d + i * vsq);
  }
  *v = 0.0;
  for (i=0; i<n; i++) {
    double *zw = d + i *vsq;
    int j;
    for (j=0; j<n; j++) {
      double *cc = c + j * vsq;
      if (j!=i) {
	for (int k=0; k<vsq; k++) zw[j] *= cc[j]; 
      }
    }
    for (int k=0; k<vsq; k++) v[k] += zw[k];
  }
  END_TALLOC_XXX1;
  END_TALLOC_XXX2;
}


int checkmal(model *cov) {
  model *next1 = cov->sub[0];
  model *next2 = cov->sub[1];
  int i, err,
    nsub = cov->nsub;

  if (next2 == NULL) next2 = next1;
  if ((err = checkplusmal(cov)) != NOERROR) {
    RETURN_ERR(err);
  }

  bool ok = OWNDOM(0) != DOMAIN_MISMATCH &&
    (equalsnowTrend(cov) || equalsnowRandom(cov) ||
    (isnowShape(cov) && (!isnowNegDef(cov) || isnowPosDef(cov) )));
  if (!ok) RETURN_ERR(ERRORNOVARIOGRAM);
   
  // to do istcftype und allgemeiner typ zulassen
  if (equalsnowTrend(cov)) {
    ok = false;
    for (i=0; i<nsub; i++)
      if ((ok = MODELNR(cov->sub[i]) == CONST || MODELNR(cov->sub[i]) == BIND))
	break;
    if (!ok) SERR2("misuse as a trend function. At least one factor must be a constant (including 'NA') or a vector built with '%.50s(...)' or '%.50s(...).",
		  DefList[BIND].name,  DefList[BIND].nick);
  }
  cov->logspeed = isXonly(OWN) ? 0 : RF_NA;
      

  if (OWNTOTALXDIM >= 2) cov->pref[TBM] = PREF_NONE;
  if (OWNTOTALXDIM==2 && cov->nsub == 2 && 
      isAnyDollar(next1) && isAnyDollar(next2)) {
    double *aniso1 = PARAM(next1, DANISO),
      *aniso2 = PARAM(next2, DANISO);
    if (aniso1 != NULL && aniso2 != NULL) {
      if (aniso1[0] == 0.0 && next1->ncol[DANISO] == 1) {
	cov->pref[TBM] = next2->pref[TBM];
      } else if (aniso2[0] == 0.0 && next2->ncol[DANISO] == 1) {
	cov->pref[TBM] = next1->pref[TBM];
      }
    }
  }
  
  if (cov->ptwise_definite <= last_pt_definite) {
    cov->ptwise_definite = next1->ptwise_definite;
    if (cov->ptwise_definite != pt_zero) {
      for (i=1; i<cov->nsub; i++) {
	model *sub = cov->sub[i];
	if (sub->ptwise_definite == pt_zero) {
	  cov->ptwise_definite = pt_zero;
	  break;
	}
	if (sub->ptwise_definite != pt_posdef) {
	  if (sub->ptwise_definite == pt_negdef) {
	    cov->ptwise_definite = 
	      cov->ptwise_definite == pt_posdef ? pt_negdef
	      : cov->ptwise_definite == pt_negdef ? pt_posdef
	      : pt_indef;
	  } else { // sub = indef
	    cov->ptwise_definite = sub->ptwise_definite;
	    break;
	  }
	}
      }
    }
  }
	  
  EXTRA_STORAGE;
  RETURN_NOERROR;
}



Types Typemal(Types required, model *cov, isotropy_type required_iso){
  bool allowed = isShape(required) || isTrend(required) ||
    equalsRandom(required);
    //||required==ProcessType ||required==GaussMethodType; not yet allowed;to do
  if (!allowed) return BadType;
  for (int i=0; i<cov->nsub; i++) {
    if (isBad(TypeConsistency(required, cov->sub[i], required_iso)))
      return BadType;
  }
  return required;
}



int initmal(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
//  int err;
//  RETURN_ERR(err);
  RETURN_ERR(ERRORFAILED);
  int 
    vdim = VDIM0;
  if (VDIM0 != VDIM1) BUG;

  for (int i=0; i<vdim; i++) cov->mpp.maxheights[i] = RF_NA;
}
void domal(model VARIABLE_IS_NOT_USED *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  BUG;
}


//////////////////////////////////////////////////////////////////////
// mpp-plus


#define PLUS_P 0 // parameter
int  CheckAndSetP(model *cov){   
  int n,
    nsub = cov->nsub;
  double 
    cump = 0.0;
  if (PisNULL(PLUS_P)) {
    assert(cov->nrow[PLUS_P] == 0 && cov->ncol[PLUS_P] == 0);
    PALLOC(PLUS_P, nsub, 1);
    for (n=0; n<nsub; n++) P(PLUS_P)[n] = 1.0 / (double) nsub;    
  } else {
    cump = 0.0;
    for(n = 0; n<nsub; n++) {
      cump += P(PLUS_P)[n];
      if (cump > 1.0 && n+1<nsub) RETURN_ERR(ERRORATOMP);
    }
    if (cump != 1.0) {
      if (nsub == 1) {
	WARN0("the p-values do not sum up to 1.\nHere only one p-value is given which must be 1.0");
	P(PLUS_P)[0] = 1.0;
      } else {
	if (cump < 1.0 && P(PLUS_P)[nsub-1] == 0) {
	  WARN1("The value of the last component of '%.50s' is increased.",
		KNAME(PLUS_P)); 

	}
	else SERR1("The components of '%.50s' do not sum up to 1.", KNAME(PLUS_P));
	P(PLUS_P)[nsub-1] = 1.0 - (cump - P(PLUS_P)[nsub-1]);
      }  
    }
  }
  RETURN_NOERROR;
}

void kappamppplus(int i, model *cov, int *nr, int *nc){
  *nr = cov->nsub;
  *nc = i < DefList[COVNR].kappas ? 1 : -1;
}

void mppplus(double *x, model *cov, double *v) { 
  int i, n,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  model *sub;
  TALLOC_XX1(z, vdimSq);

  if (hasAnyEvaluationFrame(cov)) {  
    for (i=0; i<vdimSq; i++) v[i] = 0.0;
    for (n=0; n<cov->nsub; n++, sub++) {
      sub = cov->sub[n];
      COV(x, sub, z); // urspruenglich : covlist[sub].cov(x, cov, z); ?!
      for (i=0; i<vdimSq; i++) v[i] += P(PLUS_P)[n] * z[i];
    }
  } else {
    assert(hasPoissonFrame(cov));
    BUG;
  }
  END_TALLOC_XX1;
}

int checkmppplus(model *cov) { 
  ASSERT_ONESYSTEM;
 int err, 
    size = 1;
  
  SERR("the current version does not support RMmppplus\n");
  set_maxdim(OWN, 0, MAXMPPDIM);

  if ((err = checkplusmal(cov)) != NOERROR) {
    RETURN_ERR(err);
  }

  if ((err = CheckAndSetP(cov)) != NOERROR) RETURN_ERR(err);

  if (cov->q == NULL) QALLOC(size);
  EXTRA_STORAGE;
  RETURN_NOERROR;
}


void rangempplus(model VARIABLE_IS_NOT_USED *cov, range_type *range){ 
  range->min[PLUS_P] = 0.0;
  range->max[PLUS_P] = 1.0;
  range->pmin[PLUS_P] = 0.0;
  range->pmax[PLUS_P] = 1.0;
  range->openmin[PLUS_P] = false;
  range->openmax[PLUS_P] = false;
}

int struct_mppplus(model *cov, model **newmodel){
  int m,
    //nsub = cov->nsub,
    err = NOERROR;
  
  if (!hasMaxStableFrame(cov) && !hasPoissonFrame(cov)) {
    SERR("method is not based on Poisson point process");
  }


  RETURN_ERR(ERRORNOTPROGRAMMEDYET);

  
  // Ausnahme: mppplus wird separat behandelt:
  // if (nr == MPPPLUS) return S TRUCT(shape, NULL);
 
  ASSERT_NEWMODEL_NOT_NULL;
  NEW_STORAGE_WITH_SAVE(plus);
  plus_storage *s = cov->Splus;

  for (m=0; m<cov->nsub; m++) {
    model *sub = cov->sub[m];          
    if (s->keys[m] != NULL) COV_DELETE(s->keys + m, cov);    
    if ((err = covcpy(s->keys + m, sub)) != NOERROR) RETURN_ERR(err);
    if ((err = addShapeFct(s->keys + m)) != NOERROR) RETURN_ERR(err);
    SET_CALLING(s->keys[m], cov);
  }
  RETURN_NOERROR;
}


int init_mppplus(model *cov, gen_storage *S) {
  model  *sub;
  double M2[MAXMPPVDIM], M2plus[MAXMPPVDIM], Eplus[MAXMPPVDIM], 
    maxheight[MAXMPPVDIM];
  int i,n, err,
    vdim = VDIM0;
  if (VDIM0 != VDIM1) BUG;
  if (vdim > MAXMPPVDIM) BUG;
  pgs_storage *pgs = NULL;
  for (i=0; i<vdim; i++) {
    maxheight[i] = RF_NEGINF;
    M2[i] = M2plus[i] = Eplus[i] = 0.0;
  }
    
  NEW_STORAGE_WITH_SAVE(pgs);
  pgs = cov->Spgs;
  pgs->totalmass = 0.0;
  
 // loggiven could depend on the realasation; this is determined here
  // in case the different submodels have different values for loggiven
  // formerly loggiven got the value SUBMODEL_DEP. Now, it gets the
  // value false
  cov->loggiven = wahr; 
  assert(cov->nsub > 0);
  for (n=0; n<cov->nsub; n++) {
    sub = cov->sub[n];
    //if (!sub->mpp.loc_done) 
    //SERR1("submodel %d of '++': the location of the shapes is not defined", 
    // n);
    if ((err = INIT(sub, cov->mpp.moments, S)) != NOERROR) RETURN_ERR(err);
    //if (!sub->mpp.loc_done) SERR("locations are not initialised");

    if (cov->loggiven != falsch)      // at least 1 non-loggiven
      cov->loggiven = sub->loggiven;  //  leads to non-loggiven
    if (n==0) cov->fieldreturn = sub->fieldreturn;
    else if (cov->fieldreturn != sub->fieldreturn) 
      cov->fieldreturn = (ext_bool) SUBMODEL_DEP;

     pgs->totalmass += sub->Spgs->totalmass * P(PLUS_P)[n];
    for (i=0; i<vdim; i++)   
      if (cov->mpp.maxheights[i] > maxheight[i]) 
	maxheight[i] = cov->mpp.maxheights[i];
    //  loggiven &= cov->loggiven;

    // Achtung cov->mpp.mM2 und Eplus koennten nicht direkt gesetzt
    // werden, da sie vom DefList[SUBNR].Init ueberschrieben werden !!
    
    if (cov->mpp.moments >= 1) {
      int nmP1 = sub->mpp.moments + 1;
      for (i=0; i<vdim; i++) {
	int idx = i * nmP1;
	Eplus[i] += PARAM0(sub, PLUS_P) * sub->mpp.mMplus[idx + 1]; 
      }
      if (cov->mpp.moments >= 2) {
	for (i=0; i<vdim; i++) {
	  int idx = i * nmP1;
	  M2[i] += PARAM0(sub, PLUS_P)  * sub->mpp.mM[idx + 2];
	  M2plus[i] += PARAM0(sub, PLUS_P) * sub->mpp.mM[idx + 2];
	}
      }
    }
    //assert(sub->mpp.loc_done);
  }
  for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = maxheight[i];
  //cov->mpp.refsd = RF_NA;
  //  cov->mpp.refradius = RF_NA;

  
  if (cov->mpp.moments >= 1) {
    int nmP1 = cov->mpp.moments + 1;
    for (i=0; i<vdim; i++) {
      int idx = i * nmP1;
      cov->mpp.mMplus[idx + 1] = Eplus[i];
      cov->mpp.mM[idx + 1] = RF_NA;
    }
    if (cov->mpp.moments >= 2) {
      for (i=0; i<vdim; i++) {
	 int idx = i * nmP1;
	 cov->mpp.mM[idx + 2] = M2[i];
	 cov->mpp.mMplus[idx + 2] = M2plus[i];
      }
    }
  }
  
  //  cov->loggiven = loggiven;
  // cov->fieldreturn = fieldreturn;
  //cov->mpp.loc_done = true;       
  cov->origrf = false;
  cov->rf = NULL;

  RETURN_NOERROR;
}

void do_mppplus(model *cov, gen_storage *s) {
  model *sub;
  double subselect = UNIFORM_RANDOM;
  int subnr,
    vdim = VDIM0;
  assert(VDIM0 == VDIM1);
  for (subnr=0; (subselect -= PARAM0(cov->sub[subnr], PLUS_P)) > 0; subnr++);
  cov->q[0] = (double) subnr; // wofuer???
  sub = cov->sub[subnr];
  
  DefList[SUBNR].Do(sub, s);  // nicht gatternr
  for (int i=0; i<vdim; i++) cov->mpp.maxheights[i] = sub->mpp.maxheights[i];
  cov->fieldreturn = sub->fieldreturn;
  cov->loggiven = sub->loggiven;
}

//////////////////////////////////////////////////////////////////////
// PROCESSES
//////////////////////////////////////////////////////////////////////


int checkplusmalproc(model *cov) {
  model *sub;
  int err;

  
  assert(cov->Splus != NULL);
  assert(cov->Splus->keys_given);
  
  for (int i=0; i<cov->nsub; i++) {

    sub = cov->Splus->keys[i];
     
    if (sub == NULL) 
      SERR("named submodels are not given in a sequence.");

    Types type = isTrend(sub) ? ProcessType : OWNTYPE(0); // sollte ok sein.

    
    assert(equalsCoordinateSystem(OWNISO(0)));
    
    if ((err = CHECK_THROUGHOUT(sub, cov, type, KEEPCOPY_DOM, KEEPCOPY_ISO,
				SUBMODEL_DEP, cov->frame)) != NOERROR) {
     //   if ((err= CHECK(sub, dim, xdim, type, dom, iso, SUBMODEL_DEP, frame))
     //	!= NOERROR) {
     
      
      if ((cov->calling == NULL || cov->calling->calling == NULL) &&
	  isSymmetric(OWNISO(0)) && isVariogram(OWNTYPE(0))) {
	err = CHECK_THROUGHOUT(sub, cov, type, KEEPCOPY_DOM, 
			       CoordinateSystemOf(OWNISO(0)),
			       SUBMODEL_DEP, cov->frame);
      }
      // APMI(sub);
      if (err != NOERROR) RETURN_ERR(err);
    }
   
    if (!isnowProcess(sub) && !equalsnowTrend(sub))
      RETURN_ERR(ERRORTYPECONSISTENCY);
    if (i==0) {
      VDIM0=sub->vdim[0];  // to do: inkonsistent mit vorigen Zeilen !!
      VDIM1=sub->vdim[1];  // to do: inkonsistent mit vorigen Zeilen !!
    } else {
      if (VDIM0 != sub->vdim[0] || VDIM1 != sub->vdim[1]) {
	SERR("multivariate dimensionality must be equal in the submodels");
      }
    }
  }

  RETURN_NOERROR;
}



int checkplusproc(model *cov) {
  int err;

  if ((err = checkplusmalproc(cov)) != NOERROR) {
    RETURN_ERR(err);
  }

  EXTRA_STORAGE;
  RETURN_NOERROR;
}


int structplusmalproc(model *cov, model VARIABLE_IS_NOT_USED**newmodel){
  int err,
    dim =  PREVXDIM(0);
  //  bool plus = COVNR == PLUS_PROC ;
 
  //printf("struct plusmalproc\n"); PMI(cov->calling);
  
  switch(cov->frame) {
  case GaussMethodType : {
    location_type *loc = Loc(cov);
    ONCE_NEW_STORAGE(plus);
    plus_storage *s = cov->Splus;
    s->keys_given = true;
    for (int m=0; m<cov->nsub; m++) {
      model *sub = cov->sub[m];
      bool trend = isnowTrend(sub);
      if (s->keys[m] != NULL) COV_DELETE(s->keys + m, cov);
      if ((err =  covcpy(s->keys + m, sub)) != NOERROR) {
	RETURN_ERR(err);
      }
      assert(s->keys[m] != NULL);
      assert(s->keys[m]->calling == cov);
      
      if (PL >= PL_STRUCTURE) {
	LPRINT("plus: trying initialisation of submodel #%d (%.50s).\n", m+1, 
	       NICK(sub));
      }

      addModel(s->keys + m, trend ? TREND_PROC : GAUSSPROC);
     
      //printf("isTrend = %d %s\n", isTrend(sub), NAME(sub));

     if (trend && sub->Spgs == NULL) {
	if ((err = alloc_cov(sub, dim, sub->vdim[0], sub->vdim[1])) != NOERROR)
	  RETURN_ERR(err);
      }

      SET_CALLING(s->keys[m], cov);

      
#if MAXSYSTEMS > 1    
      COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);
      int last = OWNLASTSYSTEM;
      for (int ss=0; ss<=last; ss++) {
      	BUG;
      }
#endif
 
      if ((err = CHECK(s->keys[m], loc->timespacedim, loc->timespacedim,
		       trend ? ProcessType: OWNTYPE(0), XONLY,
		       PREVISO(0), 
		       cov->vdim,
		       GaussMethodType)) != NOERROR) {
	RETURN_ERR(err);
      }

     
      if ((err = STRUCT(s->keys[m], NULL))  > NOERROR) RETURN_ERR(err);
      
    }
    
    RETURN_NOERROR;
  }
    
  default :
    SERR2("frame '%.50s' not allowed for '%.50s'", TYPE_NAMES[cov->frame],
	  NICK(cov));
  }

  RETURN_NOERROR;
}


int structplusproc(model *cov, model **newmodel){
  assert(COVNR == PLUS_PROC);
  return structplusmalproc(cov, newmodel);
}


int structmultproc(model *cov, model **newmodel){
  assert(DefList[COVNR].check == checkmultproc);
  return structplusmalproc(cov, newmodel);
}

int initplusmalproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  int  err,
    vdim = VDIM0;
  assert(VDIM0 == VDIM1);
  bool plus = COVNR == PLUS_PROC ;

  for (int i=0; i<vdim; i++) cov->mpp.maxheights[i] = RF_NA;
  if (cov->Splus == NULL || !cov->Splus->keys_given) BUG;

  if (hasGaussMethodFrame(cov)) {
 
    for (int i=0; i<cov->nsub; i++) {
      model *sub = cov->Splus != NULL && cov->Splus->keys_given
	? cov->Splus->keys[i] : cov->sub[i];
      if (!plus && 
	  (SUBNR == CONST 
	   //|| DefList[sub[0]->nr].check == checkconstant ||
	   // (isDollar(sub) && DefList[sub->sub[0]->nr].check == checkconstant)
	   ))
	continue;
      assert(cov->sub[i]->Sgen==NULL);
      cov->sub[i]->Sgen = (gen_storage *) MALLOC(sizeof(gen_storage));
      if ((err = INIT(sub, 0, cov->sub[i]->Sgen)) != NOERROR) {
	RETURN_ERR(err);
      }
      sub->simu.active = true;
    }
    cov->simu.active = true;
    RETURN_NOERROR;
  }
    
  else {
    BUG;
  }

  RETURN_ERR(ERRORFAILED);
}

 
int initplusproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  int err;
  if ((err = initplusmalproc(cov, s)) != NOERROR) RETURN_ERR(err);

  if (hasGaussMethodFrame(cov)) {
    cov->fieldreturn = (ext_bool) (cov->Splus!=NULL && cov->Splus->keys_given);
    cov->origrf = false;
    assert(cov->fieldreturn);
    if (cov->fieldreturn) cov->rf = cov->Splus->keys[0]->rf;
     
    RETURN_NOERROR;
  }
     
  else {
    BUG;
  }

  RETURN_ERR(ERRORFAILED);
}


void doplusproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  int m, i,
    total = Loc(cov)->totalpoints * VDIM0;
  double *res = cov->rf;
  assert(cov->rf == cov->Splus->keys[0]->rf);
  assert(VDIM0 == VDIM1);

  if (hasGaussMethodFrame(cov) && cov->method==SpectralTBM) {
    ERR("error in doplus with spectral");
  }
  assert(cov->Splus != NULL && cov->Splus->keys_given);

  for (m=0; m<cov->nsub; m++) {
    model *key = cov->Splus->keys[m],
      *sub = cov->sub[m];
    double *keyrf = key->rf;
    DO(key, sub->Sgen);
    if (m > 0) for(i=0; i<total; i++) res[i] += keyrf[i];
  }
  return;
}



#define MULTPROC_COPIES 0
int checkmultproc(model *cov) {
  int err;
  kdefault(cov, MULTPROC_COPIES, GLOBAL.special.multcopies);
  if ((err = checkplusmalproc(cov)) != NOERROR) {
    RETURN_ERR(err);
  }

  EXTRA_STORAGE;
  RETURN_NOERROR;
}


int initmultproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  int  err;

  if ((err = initplusmalproc(cov, s)) != NOERROR) {
    BUG;
   RETURN_ERR(err);
  }

  if (hasGaussMethodFrame(cov)) {
    ReturnOwnField(cov);
    RETURN_NOERROR;
  }
     
  else {
    BUG;
  }

  RETURN_ERR(ERRORFAILED);
}




void domultproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  double *res = cov->rf;
  assert(VDIM0 == VDIM1);
  int  m, i, c, idx,
    vdim = VDIM0,
   //  vdimSq= vdim * vdim,
    total = Loc(cov)->totalpoints,
   totalvdim = total * vdim,
   copies = GLOBAL.special.multcopies,
   factors = 0;

  if (hasGaussMethodFrame(cov) && cov->method==SpectralTBM) {
    ERR("error in do_mult with spectral");
  }
  
  if (cov->nsub==2) {
    int m0 = MODELNR(cov->sub[0]),
      m1 = MODELNR(cov->sub[1]);      
    if (((m0 == PROD_PROC) xor (idx = m1==PROD_PROC))
	&& m0 != CONST &&  m1 != CONST) {
   // koennte noch allgemeiner gemacht werden in verbindung mit CONST; todo
      copies = 1;
      // PMI(cov);
      //printf("isx=%d\n", idx);
      assert(cov->sub[idx]->qlen > PRODPROC_RANDOM);
      cov->sub[idx]->q[PRODPROC_RANDOM] = (double) false;
    }
  }

  assert(cov->Splus != NULL && cov->Splus->keys_given);
  TALLOC_X1(z, totalvdim); // da eh hoffnungslos

  SAVE_GAUSS_TRAFO;
  for (c=0; c<copies; c++) {
    for(i=0; i<totalvdim; res[i++] = 1.0);
    for (m=0; m<cov->nsub; m++) {

      if (PL > PL_RECURSIVE) {
	PRINTF("\rcopies=%d sub=%d\n", c, m);
      }

      model *key = cov->Splus->keys[m],
	*sub = cov->sub[m];
      double *keyrf = key->rf;
      if (SUBNR == CONST) {
	double 
	  cc = isnowTrend(sub) ? PARAM0(sub, CONST_C) 
	  : SQRT(PARAM0(sub, CONST_C));
	for(i=0; i<totalvdim; i++) res[i] *= cc;
	//printf("cc=%10g\n", cc);
      }
      /* else {
	bool dollar = isDollar(sub);
	model *Const = dollar ? Const->sub[0] : sub;
	if (DefList[Check->nr].check == checkconstant) {
	  if (dollar) {
	    double var = PARAM0(sub, DVAR);
	    bool random = false;
	    if (sub->kappasub[DVAR] != NULL) {
	      if (random = is Random(sub->kappasub[DVAR])) {
		Do(sub->kappasub[DVAR], sub->Sgen);
	      } else {
		TA LLOC_EXTRA2(VarMem, totalvdim);
		F ct n(NULL, sub->kappasub[DVAR], VarMem);
	      }
	    }
	    if (var != 1.0) {
	      double sd = SQRT(var);
	      for(i=0; i<totalvdim; i++) res[i] *= sd;
	    }
	  } else {
	  TA LLOC_EXT RA3(ConstMem, vdimSq);
	  }
	  } */
      else {	  
	factors ++;
	DO(key, sub->Sgen);
	for(i=0; i<totalvdim; i++) res[i] *= keyrf[i];
      }
    }
    if (factors == 1) return; // no Error, just exit
    if (c == 0) res = z;
    else for(i=0; i<totalvdim; i++) cov->rf[i] += res[i];
  }

  double f;
  f = 1 / SQRT((double) copies);
  //printf("f=%10g\n", f);
  for(i=0; i<totalvdim; i++) cov->rf[i] *= f;

  END_TALLOC_X1;

}



void rangemultproc(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[MULTPROC_COPIES] = 1.0;
  range->max[MULTPROC_COPIES] = RF_INF;
  range->pmin[MULTPROC_COPIES] = 1.0;
  range->pmax[MULTPROC_COPIES] = 1000;
  range->openmin[MULTPROC_COPIES] = false;
  range->openmax[MULTPROC_COPIES] = true;
}


// $power

void PowSstat(double *x, model *cov, double *v){
  logPowSstat(x, cov, v, NULL);
}

void logPowSstat(double *x, model *cov, double *v, double *Sign){
  model *next = cov->sub[POW_SUB];
  double 
    factor,
    var = P0(POWVAR),
    scale =P0(POWSCALE), 
    p = P0(POWPOWER),
    invscale = 1.0 / scale;
  int i,
    vdim = VDIM0,
    vdimSq = vdim *vdim,
    xdimown = OWNXDIM(0);
  assert(VDIM0 == VDIM1);
  TALLOC_X1(z, xdimown);

  for (i=0; i < xdimown; i++) z[i] = invscale * x[i];
  if (Sign==NULL) {
    COV(z, next, v);
    factor = var * POW(scale, p);
    for (i=0; i<vdimSq; i++) v[i] *= factor; 
  } else {
    LOGCOV(z, next, v, Sign);
    factor = LOG(var) + p * LOG(scale);
    for (i=0; i<vdimSq; i++) v[i] += factor; 
  }
  END_TALLOC_X1;
}

void PowSnonstat(double *x, double *y, model *cov, double *v){
  logPowSnonstat(x, y, cov, v, NULL);
}

void logPowSnonstat(double *x, double *y, model *cov, double *v, 
		 double *Sign){
  model *next = cov->sub[POW_SUB];
  double 
    factor,
    var = P0(POWVAR),
    scale =P0(POWSCALE), 
    p = P0(POWPOWER),
     invscale = 1.0 / scale;
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim,
    xdimown = OWNXDIM(0);
  assert(VDIM0 == VDIM1);
  TALLOC_X1(z1, xdimown);
  TALLOC_X2(z2, xdimown);

  for (i=0; i<xdimown; i++) {
    z1[i] = invscale * x[i];
    z2[i] = invscale * y[i];
  }

  if (Sign==NULL) {
    NONSTATCOV(z1, z2, next, v);
    factor = var * POW(scale, p);
    for (i=0; i<vdimSq; i++) v[i] *= factor; 
  } else {
    LOGNONSTATCOV(z1, z2, next, v, Sign);
    factor = LOG(var) + p * LOG(scale);
    for (i=0; i<vdimSq; i++) v[i] += factor; 
  }
  END_TALLOC_X1;
  END_TALLOC_X2;
}

 
void inversePowS(double *x, model *cov, double *v) {
  model *next = cov->sub[POW_SUB];
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  double y, 
    scale =P0(POWSCALE),
    p = P0(POWPOWER),
    var = P0(POWVAR);
 assert(VDIM0 == VDIM1);

  y= *x / (var * POW(scale, p)); // inversion, so variance becomes scale
  if (DefList[NEXTNR].inverse == ErrInverse) BUG;
  INVERSE(&y, next, v);
 
  for (i=0; i<vdimSq; i++) v[i] *= scale; 
}


int TaylorPowS(model *cov) {
  if (VDIM0 != 1) SERR("Taylor only known in the unvariate case");
  model 
    *next = cov->sub[POW_SUB];
  int i;
  double scale = PisNULL(POWSCALE) ? 1.0 : P0(POWSCALE);
  cov->taylorN = next->taylorN;  
  for (i=0; i<cov->taylorN; i++) {
    cov->taylor[i][TaylorPow] = next->taylor[i][TaylorPow];
    cov->taylor[i][TaylorConst] = next->taylor[i][TaylorConst] *
      P0(POWVAR) * POW(scale, P0(POWPOWER) - next->taylor[i][TaylorPow]);   
  }
  
  cov->tailN = next->tailN;  
  for (i=0; i<cov->tailN; i++) {
    cov->tail[i][TaylorPow] = next->tail[i][TaylorPow];
    cov->tail[i][TaylorExpPow] = next->tail[i][TaylorExpPow];
    cov->tail[i][TaylorConst] = next->tail[i][TaylorConst] *
      P0(POWVAR) * POW(scale, P0(POWPOWER) - next->tail[i][TaylorPow]);   
    cov->tail[i][TaylorExpConst] = next->tail[i][TaylorExpConst] *
      POW(scale, -next->tail[i][TaylorExpPow]);
  }
  RETURN_NOERROR;
}


int checkPowS(model *cov) {
  // hier kommt unerwartet  ein scale == nan rein ?!!
  model 
    *next = cov->sub[POW_SUB];
  int err, 
    dim = OWNLOGDIM(0),
    xdimown = OWNXDIM(0),
    xdimNeu = xdimown;

  if (!isCartesian(OWNISO(0))) RETURN_ERR(ERRORNOTCARTESIAN);
    
  kdefault(cov, POWVAR, 1.0);
  kdefault(cov, POWSCALE, 1.0);
  kdefault(cov, POWPOWER, 0.0);
  if ((err = checkkappas(cov)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  ASSERT_ONESYSTEM;
  if ((err = CHECK(next, dim, xdimNeu, OWNTYPE(0), OWNDOM(0),
		   OWNISO(0), SUBMODEL_DEP, cov->frame)) != NOERROR) {
    RETURN_ERR(err);
  }

  setbackward(cov, next);
  if ((err = TaylorPowS(cov)) != NOERROR) RETURN_ERR(err);

  EXTRA_STORAGE;
  RETURN_NOERROR;
}



Types TypePowS(Types required, model *cov, isotropy_type required_iso){
 bool allowed = isShape(required) || isTrend(required) ||
   equalsRandom(required);
    //||required==ProcessType ||required==GaussMethodType; not yet allowed;to do
  if (!allowed) return BadType;

  model *next = cov->sub[0];
  return TypeConsistency(required, next, required_iso);
}


void rangePowS(model *cov, range_type* range){
  range->min[POWVAR] = 0.0;
  range->max[POWVAR] = RF_INF;
  range->pmin[POWVAR] = 0.0;
  range->pmax[POWVAR] = 100000;
  range->openmin[POWVAR] = false;
  range->openmax[POWVAR] = true;

  range->min[POWSCALE] = 0.0;
  range->max[POWSCALE] = RF_INF;
  range->pmin[POWSCALE] = 0.0001;
  range->pmax[POWSCALE] = 10000;
  range->openmin[POWSCALE] = true;
  range->openmax[POWSCALE] = true;

   int dim = OWNLOGDIM(0);
 range->min[POWPOWER] = RF_NEGINF;
  range->max[POWPOWER] = RF_INF;
  range->pmin[POWPOWER] = -dim;
  range->pmax[POWPOWER] = +dim;
  range->openmin[POWPOWER] = true;
  range->openmax[POWPOWER] = true;
 }



void PowScaleToLoc(model *to, model *from, int VARIABLE_IS_NOT_USED depth) {
  assert(!PARAMisNULL(to, LOC_SCALE) && !PARAMisNULL(from, POWSCALE));
  PARAM(to, LOC_SCALE)[0] = PARAM0(from, POWSCALE);
}

int structPowS(model *cov, model **newmodel) {
  model
    *next = cov->sub[POW_SUB],
    *scale = cov->kappasub[POWSCALE];
  int err; 

  if (next->randomkappa) SERR("random shapes not programmed yet");

  switch (cov->frame) {
  case SmithType :  case GaussMethodType :
    ASSERT_NEWMODEL_NOT_NULL;
    
    if ((err = STRUCT(next, newmodel)) > NOERROR) RETURN_ERR(err);
    
    addModel(newmodel, POWER_DOLLAR, cov);
    if (!PisNULL(POWVAR)) kdefault(*newmodel, POWVAR, P0(POWVAR));
    if (!PisNULL(POWSCALE)) kdefault(*newmodel, POWSCALE, P0(POWSCALE));
    if (!PisNULL(POWPOWER)) kdefault(*newmodel, POWPOWER, P0(POWPOWER));
    
    break;
  case BrMethodType : case SchlatherType : {
    ASSERT_NEWMODEL_NOT_NULL;
    
    if ((err = STRUCT(next, newmodel)) > NOERROR) RETURN_ERR(err);
    
    if (!isnowRandom(scale)) SERR("unstationary scale not allowed yet");
    addModel(newmodel, LOC, cov);
    addSetDistr(newmodel, scale, PowScaleToLoc, true, MAXINT);
  }
    break;
  default :
    SERR2("'%.50s': changes in scale/variance not programmed yet for '%.50s'", 
	  NICK(cov), TYPE_NAMES[cov->frame]);      
  }  
   
  RETURN_NOERROR;
}




int initPowS(model *cov, gen_storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!
  model *next = cov->sub[POW_SUB],
    *Var = cov->kappasub[POWVAR],
    *Scale = cov->kappasub[POWSCALE];
  int 
    vdim = VDIM0,
    nm = cov->mpp.moments,
    nmvdim = (nm + 1) * vdim,
    err = NOERROR;
  bool 
    maxstable = hasMaxStableFrame(cov);// Realisationsweise 
  assert(VDIM0 == VDIM1);


  assert(cov->key == NULL || ({PMI(cov);false;}));//   
  
  if (maxstable) { // !! ohne maxstable selbst !!
    double
      var[MAXMPPVDIM],  
      p = P0(POWPOWER),
      scale = P0(POWSCALE);
    int 
      intp = (int) p,
      dim = OWNLOGDIM(0);


    // Achtung I-NIT_RANDOM ueberschreibt mpp.* !!
    if (Var != NULL) {
      int nm_neu = nm == 0 && !maxstable ? 1 : nm;
      if ((err = INIT_RANDOM(Var, nm_neu, s, P(POWVAR))) != NOERROR) 
	RETURN_ERR(err);
      int nmP1 = Var->mpp.moments + 1;
      for (int i=0; i<vdim; i++) {
	int idx = i * nmP1;
	var[i] = maxstable ? P0(POWVAR) : Var->mpp.mM[idx + 1];
      }
    } else for (int i=0; i<vdim; i++) var[i++] = P0(POWVAR);
 
    if (Scale != NULL) {
      if (p != intp)
	SERR1("random scale can be initialised only for integer values of '%.50s'",
	     KNAME(POWPOWER));
      if (dim + intp < 0) SERR("negative power cannot be calculated yet");
      int idx_s = maxstable ? nm : dim + nm + intp < 1 ? 1 : dim + nm + intp;
      if ((err = INIT_RANDOM(Scale, idx_s, s, P(POWSCALE))) != NOERROR)
	RETURN_ERR(err);
      assert(Scale->mpp.moments == 1);
      scale = maxstable ? P0(DSCALE) : Scale->mpp.mM[1];
    }
    if ((err = INIT(next, nm, s)) != NOERROR) RETURN_ERR(err);


    for (int i=0; i < nmvdim; i++) {
      cov->mpp.mM[i] = next->mpp.mM[i];
      cov->mpp.mMplus[i] = next->mpp.mMplus[i];
   }


    if (Var != NULL && !maxstable) {
      int j,
	nmP1 = Var->mpp.moments + 1;
      for (j=0; j<vdim; j++) {
	int idx = j * nmP1;
	for (int i=0; i <= nm; i++) {
	  cov->mpp.mM[i] *= Var->mpp.mM[idx + i];
	  cov->mpp.mMplus[i] *= Var->mpp.mMplus[idx + i];
	}
      }
    } else {
      int j, k;
      double pow_var;
      for (k=j=0; j<vdim; j++) { 
	pow_var = 1.0;
	for (int i=0; i<=nm; i++, pow_var *= var[j], k++) {	
	  cov->mpp.mM[k] *= pow_var;
	  cov->mpp.mMplus[k] *= pow_var;
	}
      }
    }

    if (Scale != NULL && !maxstable) {
      if (dim + nm * intp < 0 || dim + intp * nm > Scale->mpp.moments) 
	SERR("moments cannot be calculated");
      assert(Scale->vdim[0] == 1 && Scale->vdim[1] == 1 );
      for (int i=0; i <= nm; i++) {
	int idx = dim + i * intp;
	cov->mpp.mM[i] *= Scale->mpp.mM[idx];
	cov->mpp.mMplus[i] *= Scale->mpp.mMplus[idx];
      }
    } else {
      int j,k;
      double 
	pow_scale,
	pow_s = POW(scale, dim),
	pow_p = POW(scale, p);
      for (k=j=0; j<vdim; j++) { 
	pow_scale = pow_s;
	for (int i=0; i <= nm; i++, pow_scale *= pow_p, k++) {
	  cov->mpp.mM[k] *= pow_scale;
	  cov->mpp.mMplus[k] *= pow_scale;
	}
      }
    }

    bool unnormed = R_FINITE(next->mpp.unnormedmass);
    if (unnormed) {
      if (vdim > 1) BUG;
      cov->mpp.unnormedmass =
	next->mpp.unnormedmass * (var[0] * POW(scale, p + dim));
    } else cov->mpp.unnormedmass = RF_NAN;
    
    if (Scale != NULL || isnowRandom(next)) {
      if (unnormed) for (int i=0; i<vdim; i++) cov->mpp.maxheights[i] = RF_INF;
      else RETURN_ERR(ERRORNOTPROGRAMMEDYET);
    } else {
      double pp = POW(scale, p);
      for (int i=0; i<vdim; i++)   
	cov->mpp.maxheights[i] = next->mpp.maxheights[i] * (var[i] * pp);
    }
  }
  

  else if (hasGaussMethodFrame(cov)) {
    if ((err=INIT(next, 0, s)) != NOERROR) RETURN_ERR(err);
  }

  //  else if (cov->frame = = Any Type) {
  //    if ((err=INIT(next, 0, s)) != NOERROR) RETURN_ERR(err);    
  //  }

  else SERR("Initiation of scale and/or variance failed");

 
  if ((err = TaylorPowS(cov)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;
}


void doPowS(model *cov, gen_storage *s){
 
  if (hasMaxStableFrame(cov)) {
    model *next = cov->sub[POW_SUB];
         
    DO(next, s);// nicht gatternr
    double factor = P0(POWVAR) * POW(P0(POWSCALE), P0(POWPOWER));
    int i,
      vdim = VDIM0;
    assert(VDIM0 == VDIM1);
    for (i=0; i<vdim; i++)   
      cov->mpp.maxheights[i] = next->mpp.maxheights[i] * factor;
    return;
  } 
 
  BUG;
}
 
