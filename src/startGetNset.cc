/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 (library for simulation of random fields)

 Copyright (C) 2001 -- 2017 Martin Schlather, 

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

#include <stdio.h>  
#include <string.h>
#include "def.h"
#include <Basic_utils.h>
#include <Rdefines.h>
#include <R_ext/Linpack.h>
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>

#include "extern.h"
#include "questions.h"
#include "primitive.h"
#include "operator.h"
#include "variogramAndCo.h"
#include "Coordinate_systems.h"
#include "shape.h"
#include "Processes.h"
#include "init.h"
#include "RF.h"
#include "QMath.h"



int checkOK(model VARIABLE_IS_NOT_USED *cov){
   RETURN_NOERROR;
}

int checkMissing(model *cov){
  if (cov->calling == NULL) ERR0("missing may not be called by the user");
  model *calling=cov->calling;
  ERR1("'%.50s' does have not enough submodels", NICK(calling));
  RETURN_ERR(ERRORFAILED); // damit compiler keine Warnung bringt
}

int checkNotOK(model VARIABLE_IS_NOT_USED *cov){
  RETURN_ERR(ERRORFAILED);
}

void ScaleOne(double *x, model VARIABLE_IS_NOT_USED *cov, double *v){ 
  *v = *x <= 0.05 ? 1.0 : RF_NA;
} 


void addkappa(int i, const char *n, SEXPTYPE t) {
  defn *C = DefList + currentNrCov - 1;
  assert(n[0] != '\0' && 
	 (n[0] != ONEARGUMENT_NAME || n[1] != '\0') // reserved for standard parameter names
	 && (n[0] != DEFAULT_SUBNAME || 
	     (n[1] != '\0' && (n[1] < '1' ||  n[1] > '9')))
	 // reserved for standard submodel names
	 );
  assert(i < C->kappas);
  // assert(STRCMP(n, ELEMENT) || C->check == checkfix || C->check == checkmixed);
  strcopyN(C->kappanames[i], n, PARAMMAXCHAR);
  C->kappatype[i] = t;
  assert(STRCMP(n, FREEVARIABLE) || C->internal);
  // if (t >= LISTOF) assert(STRCMP(C->kappanames[0], ELEMENT) == 0 
  //			  ||  C->check == checkselect
  //			  );
}



void kappanames(const char* n1, SEXPTYPE t1) {
  assert({defn *C = DefList + currentNrCov - 1; C->kappas == 1;});
  addkappa(0, n1, t1);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2) {
  assert({defn *C = DefList + currentNrCov - 1; C->kappas == 2;});
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3) {
  assert({defn *C = DefList + currentNrCov - 1; C->kappas == 3;});
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 4;});
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
  addkappa(3, n4, t4);
}
void kappanamesX(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5) {
  addkappa(0, n1, t1);
  addkappa(1, n2, t2);
  addkappa(2, n3, t3);
  addkappa(3, n4, t4);
  addkappa(4, n5, t5);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 5;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 6;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
		const char* n7, SEXPTYPE t7) {
  assert({ defn *C = DefList + currentNrCov - 1; C->kappas == 7;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
  addkappa(6, n7, t7);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
		const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 8;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
  addkappa(6, n7, t7);
  addkappa(7, n8, t8);
}
void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 9;});
 kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
  addkappa(6, n7, t7);
  addkappa(7, n8, t8);
  addkappa(8, n9, t9);
}
void kappanamesX(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10) {
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5);
  addkappa(5, n6, t6);
  addkappa(6, n7, t7);
  addkappa(7, n8, t8);
  addkappa(8, n9, t9);
  addkappa(9, n10, t10);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 10;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11) {
  assert({ defn *C = DefList + currentNrCov - 1; C->kappas == 11;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 12;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
  addkappa(11, n12, t12);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13) {
   assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 13;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
  addkappa(11, n12, t12);
  addkappa(12, n13, t13);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 14;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
  addkappa(11, n12, t12);
  addkappa(12, n13, t13);
  addkappa(13, n14, t14);
}

void kappanamesX(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		 const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15) {
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	     n9, t9, n10, t10);
  addkappa(10, n11, t11);
  addkappa(11, n12, t12);
  addkappa(12, n13, t13);
  addkappa(13, n14, t14);
  addkappa(14, n15, t15);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 15;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	      n9, t9, n10, t10, n11, t11, n12, t12, n13, t13, n14, t14, n15, 
	      t15);  
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15, const char* n16, SEXPTYPE t16) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 16;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	      n9, t9, n10, t10, n11, t11, n12, t12, n13, t13, n14, t14, n15, 
	      t15);  
  addkappa(15, n16, t16);
}

void kappanames(const char* n1, SEXPTYPE t1, const char* n2, SEXPTYPE t2,
		const char* n3, SEXPTYPE t3, const char* n4, SEXPTYPE t4,
		const char* n5, SEXPTYPE t5, const char* n6, SEXPTYPE t6,
	       	const char* n7, SEXPTYPE t7, const char* n8, SEXPTYPE t8,
		const char* n9, SEXPTYPE t9, const char* n10, SEXPTYPE t10,
		const char* n11, SEXPTYPE t11, const char* n12, SEXPTYPE t12,
		const char* n13, SEXPTYPE t13, const char* n14, SEXPTYPE t14, 
		const char* n15, SEXPTYPE t15, const char* n16, SEXPTYPE t16,
		const char* n17, SEXPTYPE t17, const char* n18, SEXPTYPE t18
		) {
  assert({ defn *C = DefList + currentNrCov - 1;C->kappas == 18;});
  kappanamesX(n1, t1, n2, t2, n3, t3, n4, t4, n5, t5, n6, t6, n7, t7, n8, t8,
	      n9, t9, n10, t10, n11, t11, n12, t12, n13, t13, n14, t14, n15, 
	      t15);  
  addkappa(15, n16, t16);
  addkappa(16, n17, t17);
  addkappa(17, n18, t18);
}


void change_sortof(int i, sortsofparam sort) {
  defn *C = DefList + currentNrCov - 1;
  assert(i >= 0 && i < C->kappas);
  C->sortof_tab[i] = sort;
}


void change_typeof(int i, Types type, const char *names[]) {
  defn *C = DefList + currentNrCov - 1;
  assert(i >= 0 && i < C->kappas);
  assert(type == RandomType ||  // parameter might be random
	 type == ShapeType ||   // parameter must be determinist
	 type == RandomOrShapeType ||
	 type == MixedInputType ||
	 type == CharInputType ||
	 (type >= NN1 && type <= NN4));
  C->kappaParamType[i] = type;
  C->kappaParamTypeNames[i] = names;
}


void change_typeof(int i, Types type) {
  assert(type != MixedInputType);
  change_typeof(i, type, NULL);
}


void add_sortof(sortof_fct sortof) {
  defn *C = DefList + currentNrCov - 1;
  C->sortof = sortof;
}

void addsub(int i, const char *n) {
  defn *C = DefList + currentNrCov - 1;
  int j;
  assert(n[0] != ONEARGUMENT_NAME && n[0] != DEFAULT_SUBNAME);
  assert(i < MAXSUB);

  strcopyN(C->subnames[i], n, PARAMMAXCHAR);
  C->subintern[i] = false;
  for (j=0; j<C->kappas; j++) {
    assert(C->kappanames[j] != NULL);
    if ((C->subintern[i] = STRCMP(C->kappanames[j], C->subnames[i]) == 0))
      break;
  }
}

void subnames(const char* n1) {
  addsub(0, n1);
}
void subnames(const char* n1, const char* n2) {
  addsub(0, n1);
  addsub(1, n2);
}
void subnames(const char* n1, const char* n2, const char* n3) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
}
void subnames(const char* n1, const char* n2, const char* n3, const char* n4) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
  addsub(3, n4);
}
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
  addsub(3, n4);
  addsub(4, n5);
}
void subnames(const char* n1, const char* n2, const char* n3, const char* n4,
	      const char* n5, const char* n6) {
  addsub(0, n1);
  addsub(1, n2);
  addsub(2, n3);
  addsub(3, n4);
  addsub(4, n5);
  addsub(5, n6);
}


int xxx(int x) {return (int) POW(10, (double) x);}

void ErrCovX(double VARIABLE_IS_NOT_USED *x, model *cov,
	     double VARIABLE_IS_NOT_USED *v, const char *name) {
  // 
  PRINTF("\nErrCov.%s %s [%d] trafo=%d gatter=%d :\n", name,
	 NICK(cov), COVNR, TRAFONR, GATTERNR); // ok
  if (PL >= PL_ERRORS){
    PMI(cov);//
    crash();
  }
  ERR0("unallowed or undefined call of function");
}
void ErrCov(double *x, model *cov, double *v) { ErrCovX(x, cov, v, "Cov");}
void ErrD(double *x, model *cov, double *v) { ErrCovX(x, cov, v, "D");}
void ErrRnd(double *x, model *cov, double *v) { ErrCovX(x, cov, v, "Rnd");}
void ErrInverse(double *x, model *cov, double *v) { ErrCovX(x, cov, v, "Inv");}

void ErrLogCov(double VARIABLE_IS_NOT_USED *x, model *cov, 
	       double VARIABLE_IS_NOT_USED *v, 
	       double VARIABLE_IS_NOT_USED *Sign) {
  // 
  PRINTF("\nErrLogCov %s:\n", NICK(cov));
  if (PL >=  PL_ERRORS) {
    PMI(cov);//
    crash(); 
  }
  ERR0("unallowed or undefined call of function (log)");
}
void ErrCovNonstatX(double VARIABLE_IS_NOT_USED *x, 
		   double VARIABLE_IS_NOT_USED *y, model *cov, 
		    double VARIABLE_IS_NOT_USED *v, const char *name) {
  PRINTF("\nErrCovNonstat.%s %s: (%d)\n", name, NICK(cov), COVNR);
  if (PL >= PL_ERRORS) {
    PMI(cov->calling); //
    crash();
  }
  ERR1("unallowed or undefined call of '%.50s' as a kernel", NAME(cov));
}
void ErrCovNonstat(double *x, double *y, model *cov, double *v) { ErrCovNonstatX(x, y, cov, v, "Cov");}
void ErrDNonstat(double *x, double *y, model *cov, double *v) { ErrCovNonstatX(x,  y,cov, v, "D");}
void ErrRndNonstat(double *x, double *y, model *cov, double *v) { ErrCovNonstatX(x, y, cov, v, "Rnd");}
void ErrLogCovNonstat(double VARIABLE_IS_NOT_USED *x, 
		      double VARIABLE_IS_NOT_USED *y, model *cov,
		      double VARIABLE_IS_NOT_USED *v, 
		      double VARIABLE_IS_NOT_USED *Sign) {
  PRINTF("\nErrLogCovNonstat %s: (%d)\n", NICK(cov), COVNR);
  if (PL >= PL_ERRORS) {
    PMI(cov->calling);//
    crash();
  }
  ERR1("unallowed or undefined call of '%.50s' (log) as a kernel", NAME(cov));
}
void Errspectral(model *cov,
		 gen_storage VARIABLE_IS_NOT_USED *s, 
		 double VARIABLE_IS_NOT_USED *e) {
#ifdef SCHLATHERS_MACHINE
  defn *C = DefList + COVNR;
  PRINTF("\nErrspectral %s: (%d); %ld %ld\n",
	 NAME(cov), COVNR, (long) C->spectral , (long) Errspectral);
#else
  PRINTF("\nErrspectral %s: (%d)\n", NICK(cov), COVNR);
#endif
  if (PL >= PL_ERRORS) {
    PMI(cov->calling);//
    crash();
  }
  ERR0("unallowed or undefined call of spectral function");
} 

void ErrInverseNonstat(double VARIABLE_IS_NOT_USED *v, model *cov,
		       double *x, double *y) {
  x[0] = y[0] = RF_NAN;
  return;

  PRINTF("\nErrInverseNonstat %s: (%d)\n", NICK(cov), COVNR);
  if (PL >= PL_ERRORS) {
    PMI(cov->calling); //
    crash();
  }
  ERR0("unallowed or undefined call of non-domain function (inverse)");
}

char InternalName[]="-";
void kappasize1(int VARIABLE_IS_NOT_USED i,
		model VARIABLE_IS_NOT_USED *cov, int *nrow, int *ncol) {
  *nrow = *ncol = 1;
}


void rangeOK(model VARIABLE_IS_NOT_USED *cov, 
	     range_type VARIABLE_IS_NOT_USED *range) { 
  assert(DefList[COVNR].kappas == 0);
}

int structOK(model VARIABLE_IS_NOT_USED *cov, model VARIABLE_IS_NOT_USED **newmodel){ RETURN_NOERROR; }

int struct_statiso(model *cov, model **newmodel) {
  defn *C = DefList + COVNR;

  ASSERT_NEWMODEL_NOT_NULL;
  
  if (hasSmithFrame(cov) || hasPoissonFrame(cov)) {
    int i,
      vdim = VDIM0;
    for (i=0; i<vdim; i++) cov->mpp.maxheights[i] = 1.0;
     if (C->finiterange == wahr) {
      return addUnifModel(cov, 1.0, newmodel);
    } else SERR2("The function '%.50s' has inifinite support use '%.50s' to truncate the support.", NICK(cov), DefList[TRUNCSUPPORT].nick);
     RETURN_NOERROR;
  } else ILLEGAL_FRAME_STRUCT;
  
  RETURN_NOERROR;
}

int struct_failed(model *cov, model VARIABLE_IS_NOT_USED **newmodel) {
  //  PMI(cov);
  //  printf("%ld %ld %ld\n", DefList[COVNR].Struct, DefList[DISTRIBUTION].Struct, struct_failed);  crash();
  SERR4("initialization failed --  model '%.50s' (%d) does not fit (yet) the properties required by '%.50s'. %.50s",
	NICK(cov), COVNR, 
	cov->calling == NULL ? "<null>" : NICK(cov->calling),
	TRAFONR == MISMATCH || TRAFONR == UNSET ? "" : 
	"NOTE THAT THE ERROR CAN ALSO BE CAUSED BY COORDINATE TRANSFORMATION\n"
	); 
}

int initOK(model *cov, gen_storage *s) {
  defn *C = DefList + COVNR; // nicht gatternr  
  int i, err = NOERROR,
    nk = C->kappas;
  bool random = false;
  for (i=0; i<nk; i++) {
    model *ks = cov->kappasub[i];
    if (ks != NULL) {
      if (isRandom((Types) C->kappaParamType[i])) {
	random = true;
	if ((err = INIT(ks, cov->mpp.moments, s)) != NOERROR) RETURN_ERR(err);
      } else {
	SERR2("%.50s : parameter %.50s is not of random type", 
	      NICK(cov), C->kappanames[i]);
      }
    }
  }
  if (random) SERR("'initOK' not programmed yet for 'random'");
  RETURN_ERR(err); 
}

int init_failed(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  if (PL >= PL_ERRORS) { PRINTF("init failed cov=%s:\n", NICK(cov));}
  SERR1("'%.50s' cannot be initialised", NAME(cov));
  // SERR("Init failed. Compound Poisson fields are essentially only programmed for simple and isotropic shape functions (not kernels)");
  RETURN_ERR(ERRORFAILED);
}

int init_statiso(model *cov, gen_storage *s) {
  // only domain and isotropic models
  defn *C = DefList + COVNR; // nicht gatternr
  int err;

  if (cov->finiterange == wahr) {
    if  (C->finiterange == wahr) {
      // ?? to do ??
      // cov->mpp.refradius = 1.0; // jetzt ueber Inverse abgefangen
    }
  }
  
  if ((err = initOK(cov, s)) == NOERROR) RETURN_ERR(err);
 
  if (hasPoissonFrame(cov)) {
      RETURN_NOERROR;
  }


  if (PL >= PL_ERRORS){ PRINTF("init failed cov=%s:\n", NICK(cov));}

  //  PMI(cov);
  
  SERR("Call of init: Compound Poisson fields are essentially only programmed for domain and isotropic functions");

 
  RETURN_NOERROR;
}


void doOK(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  if (!cov->randomkappa) return;
  defn *C = DefList + COVNR; 
  int kappas = C->kappas;
  for (int i=0; i<kappas; i++) {
    model *sub = cov->kappasub[i];
    if (isnowRandom(sub)) {
      assert(!isnowProcess(sub));
      assert(!PisNULL(i));
      DORANDOM(sub, P(i));
    } else if (sub->randomkappa) {
      XERR(ERRORNOTPROGRAMMEDYET);
      // muesste ja irgendwie auf kappa uebertragen werden
    }
  }  
}

void do_failed(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  if (PL >= PL_ERRORS) {PRINTF("do failed for %s:\n", NICK(cov));}
  ERR0("call of do: compound Poisson fields are essentially only programmed for isotropic shape functions (not kernels)");
}

void do_random_failed(model *cov, double VARIABLE_IS_NOT_USED *v) {
  if (PL >= PL_ERRORS) {PRINTF("do_random failed for %s:\n", NICK(cov));}
  ERR0("Call of do: Compound Poisson fields are essentially only programmed for isotropic shape functions (not kernels)");
}
void do_random_ok(model VARIABLE_IS_NOT_USED *cov, double VARIABLE_IS_NOT_USED *v) {
}

void do_statiso(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) {
  
  if (hasPoissonFrame(cov) || hasMaxStableFrame(cov)) return;

  if (PL >= PL_ERRORS) {
    // PMI(cov);
    PRINTF("do_statosp failed for '%s' and frame='%s':\n", 
				NICK(cov), TYPE_NAMES[cov->frame]);
    if (PL >= PL_ERRORS) ERR0("Call of do_statiso: compound Poisson fields are essentially only programmed for isotropic shape functions (not kernels)");
  }
}

void nickname(const char *name, int nr, int type) {
  static int badname = UNSET; // ok, since only used in the initialisation
  char dummy[MAXCHAR];
  defn *C = DefList + nr; // nicht gatternr
  int sl = STRLEN(CAT_TYPE_NAMES[type]);  
  strcopyN(dummy, name, MAXCHAR-sl);
  //printf("%s %s\n",  CAT_TYPE_NAMES[type], dummy);
  SPRINTF(C->nick, "%.3s%.*s", CAT_TYPE_NAMES[type], MAXCHAR-4, dummy);
  STRCPY(CovNickNames[nr], C->nick);

  if ((int) STRLEN(name) >= (int) MAXCHAR - sl) {
    badname = nr;
  } else {
    if (badname >= 0 && badname != nr) {
      PRINTF("Warning! Nick name is truncated to '%s'.\n", 
	     DefList[badname].nick);
    }
    badname = UNSET;
  }
}

void insert_name(int curNrCov, const char *name, int type) {
  defn *C = DefList + curNrCov;
  char dummy[MAXCHAR];
  strcopyN(dummy, name, MAXCHAR);
  STRCPY(CovNames[curNrCov], dummy);
  STRCPY(C->name, dummy);
  if (STRLEN(name)>=MAXCHAR) {
    PRINTF("Warning! Covariance name is truncated to '%s'.\n", C->name);
  }
  assert(STRCMP(InternalName, name));
  nickname(name, curNrCov, type);
}

void StandardCovariance(model *cov, double *v){
  CovVario(cov, true, false, v); 
}
void StandardCovMatrix(model *cov, double *v) { 
  model *calling  = cov->calling;
  int
    dim = GetLoctsdim(cov),
    vdim = VDIM0;
   if (calling == NULL ||
      (!equalsnowInterface(calling) && !isnowProcess(calling)))
    calling = cov;
  if (calling->Spgs == NULL && alloc_cov(calling, dim, vdim, vdim) != NOERROR)
    XERR(ERRORMEMORYALLOCATION);
  
  CovarianceMatrix(cov, v); 
}
void StandardInverseCovMatrix(model *cov, double *inverse, double *det) { 
  InverseCovMatrix(cov, inverse, det); 
}
void StandardVariogram(model *cov, double *v) {
  //assert(false);
  CovVario(cov, false, false, v); 
}
void StandardPseudoVariogram(model *cov, double *v) {
  CovVario(cov, false, true, v); 
}

void StandardInverseNonstat(double *v, model *cov,
			    double *left, double *right) {
  assert(DefList[COVNR].inverse != NULL);
  double x;
  int d,
    dim = PREVLOGDIM(0); // !!und nicht OWNXDIM, xdimprev !!
  if (!equal_coordinate_systems(PREV, OWN)) BUG; // indeed,
  // some kind of reverse of the #-fctns should be implemented
   
  INVERSE(v, cov, &x);

  for (d=0; d<dim; d++) {
    left[d] = -x;
    right[d] = x;
  }
}

void StandardLogInverseNonstat(double *v, model *cov,
			    double *left, double *right) {
  assert(DefList[COVNR].inverse != NULL);
  double x, 
    w = EXP(*v);
  int d,
    dim = PREVLOGDIM(0); // !!und nicht OWNXDIM, xdimprev !!
  if (!equal_coordinate_systems(PREV, OWN)) BUG; // indeed,
  // some kind of reverse of the #-fctns should be implemented
 
  INVERSE(&w, cov, &x);
  for (d=0; d<dim; d++) {
    //    printf("inverse %d %10g %10g\n", d, *v, x);
    left[d] = -x;
    right[d] = x;
  }
}  

void StandardNonLogDistrD(double *x, model *cov, double *D) {
  VTLG_DLOG(x, cov, D);
  *D = EXP(*D);
}


char isTrue(model VARIABLE_IS_NOT_USED *cov) {ERR0("isTrue may not be used anymore"); return true;}
char isFalse(model VARIABLE_IS_NOT_USED *cov) {return false;}

void InverseFiniteRange(double VARIABLE_IS_NOT_USED *x, model VARIABLE_IS_NOT_USED *cov, double *v){ *v = 1.0; }

void InverseIsotropic(double *U, model *cov, double *inverse){
  // in case of vdim > 1, just take the first component.
  // Of course also mea over all components might be passible as well

#define step 2.0
#define tol 1e-8
#define max_it 30
  
  if (VDIM0 != VDIM1 || OWNLASTSYSTEM != 0) BUG;
  int vdim = VDIM0;
  double v[MAXVDIM * MAXVDIM], wert[MAXVDIM * MAXVDIM];
  if (vdim > MAXVDIM) BUG;
 
  double left, right, middle, leftinverse,
    x = 0.0,
    u = *U;
  bool greater;
  int i;
  

  COV(&x, cov, v);
  greater = *v > u;
  *wert = *v;
  x = step;  
  for (i=0; i<max_it; i++) {
    leftinverse = *wert;
    COV(&x, cov, wert);
    if (greater xor (*wert >= u)) break;
    x *= step;
  }
  if (i >= max_it) {
    *inverse = FABS(*v - u) <= FABS(*wert - u) ? 0.0 : RF_INF;
    return;
  }
  *inverse = *wert;
  right = x;
  left = x == step ? 0 : x / step;
  for (i=0; i<max_it; i++) {
    middle = 0.5 * (left + right);
    COV(&middle, cov, wert);    
    if (greater xor (*wert >= u)) {
      right = middle;
   } else {
      left = middle;
      leftinverse = *wert;
    }
  }
 
  *inverse = leftinverse == u ? left : right;

  //  if (PL > 1) {PMI(cov); PRINTF("inverse: %10g -> %10g\n", u, *inverse);}

}

void InverseIsoMon(double VARIABLE_IS_NOT_USED *x, model *cov, double VARIABLE_IS_NOT_USED *v){
  searchInverse(DefList[COVNR].cov, cov, 1.0, *x, 0.001);
}


bool allowedDstandard(model *cov) {
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  if (allowedD(sub)) return true;
  MEMCOPY(cov->allowedD, sub->allowedD, sizeof(allowedD_type));
  return false;
}



bool allowedIstandard(model *cov) {
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  if (allowedI(sub)) return true;
  MEMCOPY(cov->allowedI, sub->allowedI, sizeof(allowedI_type));
  return false;
}


bool allowedPrevModelI(model *cov) {
  defn *C = DefList + COVNR;
  int nsub = cov->nsub,
    z = 0;
  model *sub[MAXPARAM + MAXSUB];
  for (int i=0; z<nsub; i++) if (cov->sub[i]!=NULL) sub[z++] = cov->sub[i];
  for (int i=0; i<C->kappas; i++)
    if (cov->kappasub[i] != NULL) sub[z++] = cov->kappasub[i];

  bool unspecified = allowedIsubs(cov, sub, z);
  if (isMathDef(C) && (C->cov==MathCos || C->cov==MathSin ||  C->cov==MathTan)
      ) {
    assert(!unspecified);
    for(int i=FIRST_EARTH; i<=LAST_EARTH; cov->allowedI[i++] = false);
  } // else printf("not ok: %s\n", NAME(cov));
 

  return unspecified;
}


void createmodel(const char *name, Types type, int kappas, size_fct kappasize,	
		 domain_type domain, isotropy_type isotropy,
		 checkfct check, rangefct range, int vdim, pref_shorttype pref,
		 int maxdim, ext_bool finiterange, monotone_type monotone) {
  defn *C = DefList + currentNrCov;
  bool 
    stat_iso = equalsXonly(domain) && equalsIsotropic(isotropy);


  // printf("%s %d\n", name, currentNrCov); 

  if (PL >= PL_DETAILS) {
    PRINTF("%d %s vdim=%d statiso=%d iso=%d type=%d\n",
	   currentNrCov, name, vdim, stat_iso, isotropy, type);
  }
  
  assert(currentNrCov != UNSET); 
  assert(DefList!=NULL && currentNrCov>=0);
  assert(currentNrCov < MAXNRCOVFCTS);
    // ERR1("maximum number of covariance models reached. Last included  model is '%.50s'.", DefList[MAXNRCOVFCTS-1].name);
  C->TypeFct = NULL;
  assert(type >=0 && type <= OtherType);


  assert((isotropy >= 0 && (isotropy != ISO_MISMATCH || equalsRandom(type))) 
	 || isMathDef(type));
  assert(!equalsParamDepD(domain) || equalsParamDepI(isotropy));


  for (int i=0; i<MAXVARIANTS; i++) SYSTEM_NULL(C->systems[i], MAXSYSTEMS);

#if MAXSYSTEMS == 1
  set_system(C->systems[0], 0, maxdim==1 ? 1 : UNSET, maxdim, 
	       equalsSpaceIsotropic(isotropy) ? 2 : maxdim == 1 ? 1 : UNSET,
	       type, domain, 
	       isotropy, false);
#else
    set_system(C->systems[0], 0, maxdim==1 ? 1 : UNSET, maxdim, 
	       equalsSpaceIsotropic(isotropy) || maxdim == 1 ? 1 : UNSET,
	       type, domain, 
	       equalsSpaceIsotropic(isotropy) ? ISOTROPIC : isotropy, false);
    if (equalsSpaceIsotropic(isotropy))
      set_system(C->systems[0], 1, 1, 1, 1, SameAsPrevType, domain, ISOTROPIC, 
		 false);
#endif
  set_nr(C->systems[0], currentNrCov);
  C->variants = 1;


  if ((finiterange == wahr && isPosDef(type) && vdim == SCALAR) || 
      monotone == COMPLETELY_MON) {
    assert(LASTSYSTEM(C->systems[0]) == 0);
    set_system(C->systems[C->variants], 0, 2, 2, 2,
	       PosDefType, domain, SPHERICAL_ISOTROPIC, false);
    set_nr(C->systems[C->variants], currentNrCov);    
    //   printf("name %s\n", name); PSYS(C->systems[C->variants]);  assert(STRCMP(C->name, "exp") != 0);

    C->variants++;
   }
  insert_name(currentNrCov, name, type);
  C->Dallowed = equalsSubModelD(domain) ? allowedDstandard : NULL;
  C->Iallowed = equalsSubModelI(isotropy) ? allowedIstandard :
    equalsPrevModelI(isotropy) ? allowedPrevModelI : NULL;

  //if (domain == PARAM_DEP) printf("%s\n", C->nick);

  // nachfolgendes unbedingt lassen, da sonst vor Aufruf von Interface
  // reduzierende Gatter gesetzt werden koennten
  assert(!isInterface(type) ||
	 (equalsXonly(domain) && equalsUnreduced(isotropy)));

  C->kappas = kappas;
  assert(kappas >= 0 && kappas <= MAXPARAM);
  C->minsub = C->maxsub = 0;
  C->vdim = vdim;
  
  C->maxmoments = isShape(type) ? 0 : MISMATCH;
  for (int i=0; i<kappas; i++) {
    SPRINTF(C->kappanames[i], "%c%d", ONEARGUMENT_NAME, i); // default (repeated twice)
    C->kappatype[i] = REALSXP;
  }
  C->kappasize = (kappasize == NULL) ? kappasize1 : kappasize;
  C->sortof = NULL;
  if (isProcess(type)) {
    // UNBEDEINGT BIS ZUM SCHLUSS LAUFEN LASSEN, DA VERDECKTE PARAMTER
    for (int i=0; i<MAXPARAM; i++) C->sortof_tab[i] = FORBIDDENPARAM;
  } else {
    for (int i=0; i<MAXPARAM; i++) C->sortof_tab[i] = ANYPARAM;
  }
  if (isMathDef(type)) {
    for (int i=0; i<MAXPARAM; i++) {
      C->kappaParamType[i] = ShapeType;//29.12.17
      C->kappaParamTypeNames[i] = NULL;
    }
  } else {
    for (int i=0; i<MAXPARAM; i++) {
      C->kappaParamType[i] = RandomType;
      C->kappaParamTypeNames[i] = NULL;
    }
  }

  if (kappas==0) {
    assert(range == NULL);
    C->range=rangeOK; 
  } else {
    assert(range != NULL);
    C->range= range;
  }
  // printf("check\n");
  C->check = check;
  assert(check != NULL);
  for (int i=0; i<Forbidden; C->implemented[i++] = NOT_IMPLEMENTED);
  C->internal = false;
  C->Specific = isProcess(type) || isInterface(type) ? MISMATCH : UNSET;
  assert(finiterange != Paramdep || check==checktbmop);
  C->finiterange = finiterange;
  assert(monotone != MON_MISMATCH || 
	 (monotone == MON_MISMATCH && (type==RandomType || type== PointShapeType
				       || type== InterfaceType
				       || isProcess(type)
				       )));
  // assert(monotone != PARAM_DEP);
  C->Monotone = monotone;
  C->ptwise_definite = !isShape(type) && type != MathDefType ? pt_mismatch
    : isTcf(type) || isBernstein(monotone)
    || (isVariogram(type) && isMonotone(monotone) && C->vdim == 1)
    ? pt_posdef : pt_unknown;


  MEMCOPY(C->pref, pref, sizeof(pref_shorttype));
 
  C->cov = ErrCov;
  C->logD = C->D = C->D2 = C->D3 = C->D4 =C->tbm2 = C->nabla = C->hess = ErrD;
  C->random = ErrRnd;
  C->nonstat_inverse =  C->nonstat_loginverse = C->nonstat_inverse_D =
    ErrInverseNonstat;
  // printf("log\n");
  //C->density = MISMATCH;
  C->log = ErrLogCov;
  C->F_derivs = C->RS_derivs = isProcess(type) ? 0 :
    isInterface(type) ? MISMATCH : UNSET;
  C->nonstat_cov  = ErrCovNonstat;
  C->nonstat_D = ErrDNonstat;
  C->nonstat_random = ErrRndNonstat;
  C->nonstatlog = ErrLogCovNonstat;
  C->aux_cov = NULL;
  C->coinit = C->ieinit = NULL;
  C->alternative = NULL;

  C->spectral=Errspectral;

  C->drawmix = NULL;
  C->logmixdens = NULL;

  //  if (isVariogram(type) || isShape(type) || 
  //    type == MathDefType) C->logD = standard_likelihood;

  C->inverse = ErrInverse;
  C->Struct = struct_failed;
  C->Init = init_failed;
  C->Do = do_failed;
  C->FinalDo = NULL;

  //  printf("do\n");
  
  if (LASTSYSTEM(C->systems[0]) == 0) {
    if (finiterange == wahr) C->inverse = InverseFiniteRange;
    else if (stat_iso) C->inverse = InverseIsotropic; 
    if  (stat_iso) {
      C->Struct = struct_statiso;
      C->Init = init_statiso;
      C->Do = do_statiso;
    }
  }  

  C->DoRandom = do_random_failed;
  C->minmaxeigenvalue = NULL;

  C->hyperplane=NULL;  
  C->primitive = true;
  C->covariance = StandardCovariance;
  C->covmatrix = StandardCovMatrix;
  C->inversecovmatrix = StandardInverseCovMatrix;
  C->variogram = StandardVariogram;
  C->pseudovariogram = StandardPseudoVariogram;
  C->is_covmatrix= isFalse;

  C->TaylorN = C->TailN = isShape(type) ? UNSET : MISMATCH;
  C->setDI = NULL;
  
  currentNrCov++;
  //  printf("end %s\n", name);
  
}


bool isDummyInit(initfct Init) {
  return Init == init_statiso || Init == init_failed;
}


int CopyModel(const char *name, int which) {
  MEMCOPYX(DefList + currentNrCov, DefList + which, sizeof(defn)); 
  int type = SYSTYPE(DefList[which].systems[0], 0);
  assert(type <= OtherType);
  insert_name(currentNrCov, name, type);
  currentNrCov++;
  return currentNrCov - 1;
}


int CopyModel(const char *name, int which, Types type) {
  CopyModel(name, which);
  int nr = currentNrCov - 1;
  assert(DefList[nr].variants == 1);
  set_type(DefList[nr].systems[0], 0, type);
  return nr;
}


int CopyModel(const char *name, int which, checkfct check) {  
  CopyModel(name, which);
  int nr = currentNrCov - 1;  
  DefList[nr].check = check;
  return nr;
}

/*
  int CopyModel(const char *name, int which, checkfct check) {
  CopyModel(name, which);
  int nr = currentNrCov - 1;
  DefList[nr].check = check;
  }
*/

void nickname(const char *name) {
  int nr = currentNrCov - 1;
  int type = SYSTYPE(DefList[nr].systems[0], 0);
  assert(type <= MathDefType);
  nickname(name, nr, type);
}

int IncludePrim(const char *name, Types type,  int kappas, 
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, int maxdim, 
		ext_bool finiterange, monotone_type monotonicity) {  
  createmodel(name, type, kappas, NULL, domain, isotropy, check, range,
	      SCALAR, PREF_ALL, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}
int IncludePrim(const char *name,Types type,  int kappas, size_fct kappasize,
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, int maxdim, 
		ext_bool finiterange, monotone_type monotonicity) {  
  createmodel(name, type, kappas, kappasize, domain, isotropy, check, range,
	      SCALAR, PREF_ALL, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}

int IncludePrim(const char *name, Types type,  int kappas, 
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, int vdim, 
		int maxdim, ext_bool finiterange, monotone_type monotonicity) {
  createmodel(name, type, kappas, NULL, domain, isotropy, check, range,
	      vdim, PREF_ALL, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}
int IncludePrim(const char *name,Types type,  int kappas, size_fct kappasize,
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, int vdim, 
		int maxdim, ext_bool finiterange, monotone_type monotonicity) {
  createmodel(name, type, kappas, kappasize, domain, isotropy, check, range,
	      vdim, PREF_ALL, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}

int IncludePrim(const char *name,Types type,  int kappas, 
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, pref_type pref,
		int vdim, int maxdim, ext_bool finiterange, 
		monotone_type monotonicity) {  
  createmodel(name, type, kappas, NULL, domain, isotropy, check, range,
	      vdim, pref, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}

int IncludePrim(const char *name,Types type,  int kappas, size_fct kappasize,
		domain_type domain, isotropy_type isotropy,	
		checkfct check, rangefct range, pref_type pref,
		int vdim, int maxdim, ext_bool finiterange, 
		monotone_type monotonicity) { 
  createmodel(name, type, kappas, kappasize, domain, isotropy, check, range,
	      vdim, pref, maxdim, finiterange, monotonicity);
  return currentNrCov - 1;
}


void make_internal() {  
  int nr = currentNrCov - 1;  
  defn *C = DefList + nr; // nicht gatternr
  C->internal = true; 
}



// extern ?!
int IncludeModel(const char *name, Types type, int minsub, int maxsub,
		 int kappas, size_fct kappasize,
		 domain_type domain, isotropy_type isotropy,
		 checkfct check, rangefct range, pref_type pref, 
		 int internal,
		 int vdim, int maxdim, ext_bool finiterange,
		 monotone_type monotonicity) {  
  createmodel(name, type, kappas, kappasize, domain, isotropy, check, range,
	      vdim, pref, maxdim, finiterange, monotonicity);
  //    assert(maxsub > 0); // check deleted 25. nov 2008 due to nugget 
  // printf("name = %s\n", name);
  assert(!equalsPrevModelI(isotropy) || maxsub != 0 || isMathDef(type)
	 || STRCMP("trend", name) == 0
	 || STRCMP("U", name) == 0 
	 || STRCMP("declare", name) == 0 
	 || STRCMP("whittle", name) == 0 
	 || STRCMP("matern", name) == 0 
	 || STRCMP("id", name) == 0
	 || STRCMP("constant", name) == 0);

  assert(maxsub >= minsub && maxsub <= MAXSUB);
  assert(check != checkOK || maxsub==0);
  int i, 
    nr = currentNrCov - 1;  
  defn *C = DefList + nr; // nicht gatternr
  C->minsub = minsub;
  C->maxsub = maxsub;  

  /*
    if (maxsub == 0) printf("not primitive: %s \n", name); //
    not primitive: biWM 
    not primitive: constant 
    not primitive: epsC 
    not primitive: matern 
    not primitive: nugget 
    not primitive: trend 
    not primitive: whittle 
    (and others?)
  */

  assert(minsub <= maxsub);
  if (false && PL>=PL_SUBDETAILS && maxsub == 0) {
    PRINTF("note: %s has no submodels\n", C->name);
  }
  C->primitive = false; // muss falls sein, sonst wird kappacheck 
  // aufgerufen
  C->internal = internal;


  //  if (internal == 1) printf("internal: %s %s\n", C->nick, C->name);

  if (maxsub <= 2) {
    if (maxsub >= 1) {
      addsub(0, "phi");
      if (maxsub >= 2) {
	addsub(1, "psi");
      }
    }
  } else {
    for (i=0; i<maxsub; i++) {      
      SPRINTF(C->subnames[i], "C%d", i); // default (repeated twice)
      C->subintern[i] = false;
    }
  }
  return nr;
}

int IncludeScalar(const char *name, Types type, int minsub, int maxsub,
		 int kappas, 
		 domain_type domain, isotropy_type isotropy,
		 checkfct check, rangefct range, pref_type pref, 
		 int maxdim, ext_bool finiterange, 
		 monotone_type monotonicity) {
  return
    IncludeModel(name, type, minsub, maxsub, kappas,
		 NULL, domain, isotropy, check, 
		 range, pref, false,
		 SCALAR, maxdim, finiterange, monotonicity);
}

bool addvariantOK(Types type, isotropy_type iso) {
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
   if (isManifold(SYSTYPE(C->systems[0], 0))) return false;
  if (C->variants >= MAXVARIANTS) return false;
  if (!(nr == BIND)
      && equalsSubModelD(DOM(C->systems[0], 0))) return false;
  if (LASTSYSTEM(C->systems[0]) > 0) BUG;
   if (isPrevModelI(C) || equalsPrevModelI(iso)) {
    if (C->check != checkconstant && 
	SYSTYPE(C->systems[0], 0) != MathDefType &&
	C->check != checkcovariate &&
	C->check != checkMatern &&
	C->check != checkWM &&
	C->check != checkdeclare
	) return false; 
  }
  int n = SYSTEMS(C->systems[C->variants - 1]);
  for (int s=0; s<n; s++) {
    isotropy_type sysiso = ISO(C->systems[C->variants - 1], s);
 if (equal_coordinate_system(sysiso, iso, true)) {
       if (sysiso > iso && C->check != checkpower) return false;  // see check2x
       if (sysiso != iso ||  // see check2x
	 !isBad(TypeConsistency(type, SYSTYPE(C->systems[C->variants - 1],
						  s))))
	return false;
    }
   }
  Types systype = SYSTYPE(C->systems[0], 0);
   if (!isNegDef(systype) &&
     systype != type &&  
     !equalsShape(systype) &&
     !isProcess(systype) &&
     !isMathDef(C) &&
     C->check != checktrend) 
    return false;// see also 
   if (isAnySphericalIso(iso) &&
      ((C->finiterange == wahr && isPosDef(type) && C->vdim == SCALAR)
       || C->Monotone == COMPLETELY_MON)) return false;
     return true;
}

void AddVariant(Types type, isotropy_type iso) {
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
  
  //  printf("addvariant %s\n", C->name);  PSYS(C->systems[0]);
 
  assert(addvariantOK(type, iso));
  set_system(C->systems[C->variants], 0,
	     LOGDIM(C->systems[0], 0),
	     MAXDIM(C->systems[0], 0),
	     XDIM(C->systems[0], 0),
	     type,
	     DOM(C->systems[0], 0),
	     iso, false);
  //  printf("addvariant %s\n", C->name);
  //  PSYS(C->systems[0]);
  // PSYS(C->systems[C->variants]);  assert(STRCMP(C->name, "exp") != 0);
  set_nr(C->systems[C->variants], nr);
  C->variants++; 
}


#define IMPLEMENTED_CE							\
  ( (anyVariant(isPosDef, C) || anyVariant(isManifold, C)) &&		\
    !equalsKernel(DOM(C->systems[0], 0)) )

#define IMPLEMENTED_SEQUENTIAL (C->vdim <= 1 &&	IMPLEMENTED_CE)

void addCov(int F_derivs, covfct cf, covfct D, covfct D2, covfct inverse,
	    nonstat_inv nonstat_inverse) {
  int nr = currentNrCov - 1;

  
  assert(nr>=0 && nr<currentNrCov && cf!=NULL);
  defn *C = DefList + nr; // nicht gatternr
  bool stat_iso = isIsotropicXonly(C->systems[0]);

  //  printf("addcov %d\n", stat_iso); if (nr > 114) BUG;;
 
  C->cov = cf;  
  if (C->RS_derivs < 0) C->RS_derivs = 0;
  assert(C->nonstat_cov == ErrCovNonstat);

  if (D != NULL) {
    assert(cf != NULL);
    if (C->cov!=NULL && C->RS_derivs < 1) C->RS_derivs = 1;
    C->D=D;
    
    assert(isSpaceIsotropic(C->systems[0])|| isSubModelI(C) ||
	   isPrevModelI(C) || isParamDepI(C));
    // C->implemented[TBM2] = NUM_APPROX;
    C->implemented[TBM] = true; 
  }
  if (D2 != NULL) {
    assert(D != NULL);
    C->D2 = D2;
   if (C->cov!=NULL && C->D != NULL && C->RS_derivs < 2)  C->RS_derivs = 2;
  }
  if (inverse != NULL) C->inverse = inverse;
  else if (isMonotone(C->Monotone) && isIsotropic(C->systems[0]) && 
	   C->inverse==ErrInverse)
    C->inverse = InverseIsoMon;
  if (stat_iso && C->inverse != ErrInverse) 
    C->nonstat_loginverse = StandardLogInverseNonstat;

  //  if (nonstat_inverse != NULL) printf("nonstat %s\n", C->nick);

  C->nonstat_inverse = nonstat_inverse!=NULL ? nonstat_inverse :
    stat_iso && inverse != NULL ? StandardInverseNonstat : ErrInverseNonstat;
  C->implemented[Direct] = cf != NULL;

  // here * posdef=0 mani=1 xonly=0 prevD=0
  //  printf("here %s posdef=%d mani=%d xonly=%d prevD=%d\n", C->name, anyVariant(isPosDef, C), anyVariant(isManifold, C),  isXonly(DOM(C->systems[0], 0)), isPrevModelD(DOM(C->systems[0], 0)) );
  
  C->implemented[CircEmbed] = cf != NULL && IMPLEMENTED_CE;

  // printf("%s %d %d\n", C->nick, C->implemented[CircEmbed], C->pref[CircEmbed]);

  C->implemented[Sequential] = IMPLEMENTED_SEQUENTIAL;
  C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
  assert(C->F_derivs <= C->RS_derivs);
  assert(C->F_derivs >= 0);
}

void addCov(covfct cf, covfct D, covfct D2, covfct inverse) {
  addCov(UNSET, cf, D, D2, inverse, NULL);
}

void addCov(covfct cf, covfct D, covfct D2, nonstat_inv inverse) {
  addCov(UNSET, cf, D, D2, NULL, inverse);
}

void addCov(covfct cf, covfct D, covfct D2, 
	    covfct inverse, nonstat_inv nonstat_inverse) {
  addCov(UNSET, cf, D, D2, inverse, nonstat_inverse);
}

void addCov(int F_derivs, covfct cf, covfct D, covfct inverse) {
  addCov(F_derivs, cf, D, NULL, inverse, NULL);
}


void addCov(covfct cf, covfct D, covfct inverse) {
  addCov(UNSET, cf, D, inverse);
}


void addCov(int F_derivs, covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse, nonstat_inv nonstat_inverse) {
  int nr = currentNrCov - 1;
  assert(nr>=0 && nr<currentNrCov);
  addCov(UNSET, cf, D, D2, inverse, nonstat_inverse);
  defn *C = DefList + nr;
  C->D3 = D3;
  assert(C->RS_derivs == 2 && D3 != NULL);
  if (D4 == NULL) {
    C->RS_derivs = 3;
  } else {
    C->RS_derivs = 4;
    C->D4 = D4;
  }
  C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
  assert(C->F_derivs <= C->RS_derivs);
}


void addCov(covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse) {
  addCov(-1, cf, D, D2, D3, D4, inverse, NULL);
}

void addCov(covfct cf, covfct D, covfct D2, covfct D3, covfct D4,
	    covfct inverse, nonstat_inv nonstat_inverse) {
  addCov(-1, cf, D, D2, D3, D4, inverse, nonstat_inverse);
}



void addCov(int F_derivs, nonstat_covfct cf) {
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov) && cf!=NULL);
  defn *C = DefList + nr; // nicht gatternr

  C->nonstat_cov = cf;
  C->implemented[Direct] = true;
  C->implemented[CircEmbed] = IMPLEMENTED_CE;
  C->implemented[Sequential] = IMPLEMENTED_SEQUENTIAL;

  if (C->RS_derivs < 0) {
    C->RS_derivs = 0;
    C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
    assert(C->F_derivs <= C->RS_derivs);
  }
}
void addCov(nonstat_covfct cf) {
  addCov(-1, cf);
}

void addCov(aux_covfct auxcf){
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov) && auxcf!=NULL);
  defn *C = DefList + nr; // nicht gatternr
  assert(C->cov == ErrCov && C->nonstat_cov==ErrCovNonstat);
  C->aux_cov = auxcf;
}

void addCov(covfct distrD, covfct logdistrD, nonstat_inv Dinverse,
	    covfct distrP, nonstat_covfct distrP2sided,
	    covfct distrQ, covfct distrR, nonstat_covfct distrR2sided) {
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
  assert(nr>=0 && nr<currentNrCov);
  assert(equalsDomMismatch(DOM(C->systems[0], 0)));
  assert(equalsIsoMismatch(ISO(C->systems[0], 0)));
  assert(distrD != NULL && logdistrD!=NULL && Dinverse != NULL &&
	 distrP != NULL && distrQ != NULL &&
	 distrR != NULL);
  assert(C->Struct != struct_failed && C->Init !=init_failed);

  C->RS_derivs = 1;
  C->F_derivs = 0;

  C->cov = distrP;
  C->D = distrD; 
  C->logD = logdistrD;
  C->nonstat_inverse_D = Dinverse;
  C->inverse = distrQ;
  C->random = distrR;
  if (distrP2sided != NULL || distrR2sided !=NULL) {
    assert(distrP2sided != NULL && distrR2sided !=NULL);
    C->nonstat_cov = distrP2sided;
    C->nonstat_random = distrR2sided;
    set_dom(C->systems[0], 0, PREVMODEL_D);
  } else {
    assert(distrP2sided == NULL && distrR2sided ==NULL);
    set_dom(C->systems[0], 0, XONLY);
  }
  set_iso(C->systems[0], 0, CARTESIAN_COORD);
}


void addlogD(covfct logdistrD) {
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
  assert(nr>=0 && nr<currentNrCov);
  assert(isProcess(SYSTYPE(C->systems[0], 0)));
  assert(equalsXonly(DOM(C->systems[0], 0)));
  assert(equalsUnreduced(ISO(C->systems[0], 0)));
  assert(logdistrD != NULL);
  assert(C->Struct != struct_failed);
  assert( C->Init !=init_failed);

  C->RS_derivs = 1;
  C->F_derivs = 0;

  C->logD = logdistrD;
  C->D = StandardNonLogDistrD;
}

/*
void addlogD(covfct logD) {
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov));
  defn *C = DefList + nr; // nicht gatternr
  assert(logD != NULL);
  C->logD = logD;
}
*/

void addlogCov(logfct log, nonstat_logfct nonstatlog, 
	       nonstat_inv nonstat_loginverse) {
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov));
  defn *C = DefList + nr; // nicht gatternr
  if (log != NULL) C->log = log;
  else assert(nonstatlog != NULL);
  if (nonstatlog != NULL) C->nonstatlog = nonstatlog;
  if (nonstat_loginverse!=NULL)  C->nonstat_loginverse = nonstat_loginverse;
}
 

void addlogCov(logfct log) {
  addlogCov(log, NULL, NULL);
}

void addlogCov(nonstat_logfct nonstatlog) {
  addlogCov(NULL, nonstatlog, NULL);
}


 
int addFurtherCov(int F_derivs, covfct cf, covfct D, covfct D2) {
  assert(currentNrCov > 0);
  defn *C = DefList + currentNrCov;
  MEMCOPY(C, C - 1, sizeof(defn));
  assert(C->vdim == SCALAR || C->vdim == SUBMODEL_DEP ||
	 C->vdim == PARAM_DEP || D == NULL);

  strcopyN(CovNames[currentNrCov], InternalName, MAXCHAR);
  C->name[0] = InternalName[0];
  strcopyN(C->name + 1, DefList[currentNrCov-1].name, MAXCHAR - 1);
  if (cf != NULL) {
    C->cov = cf;
    C->RS_derivs = 0;
  }
  if (D != NULL) {
    assert(cf != NULL);
    C->D = D;
    C->RS_derivs = 1;

    assert(isSpaceIsotropic(C->systems[0]) || isPrevModelI(C) ||isSubModelI(C));
 
    //C->implemented[TBM2] = NUM_APPROX;
    C->implemented[TBM] = IMPLEMENTED; 
  }
  if (D2 != NULL) {
    assert(D != NULL);
    C->D2 = D2;
    C->RS_derivs = 2;
  }
  C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
  assert(C->F_derivs <= C->RS_derivs);
  C->internal = true; // addCov is used without previous call of IncludeModel
  currentNrCov++;
  return currentNrCov - 1;
}
int addFurtherCov(covfct cf, covfct D, covfct D2) {
  return addFurtherCov(-1, cf, D, D2);
}

int addFurtherCov(covfct cf, covfct D) {
  return addFurtherCov(cf, D, NULL);
}
int addFurtherCov(int F_derivs, covfct cf, covfct D) {
  return addFurtherCov(F_derivs, cf, D, NULL);
}
int addFurtherCov(int F_derivs, nonstat_covfct cf) {
  assert(currentNrCov > 0);
  defn *C = DefList + currentNrCov;
  MEMCOPY(C, C - 1, sizeof(defn));
  strcopyN(CovNames[currentNrCov], InternalName, MAXCHAR);
  C->name[0] = InternalName[0];
  strcopyN(C->name + 1, DefList[currentNrCov-1].name, MAXCHAR - 1);
  C->RS_derivs = MISMATCH;
  if (cf != NULL) {
    C->nonstat_cov = cf;
    C->RS_derivs = 0;
  }
  C->F_derivs = F_derivs >= 0 ?  F_derivs : C->RS_derivs;
  assert(C->F_derivs <= C->RS_derivs);
  C->D = ErrD;
  C->internal = true; // addCov is used without previous call of IncludeModel
  currentNrCov++;
  return currentNrCov - 1;
}

int addFurtherCov(nonstat_covfct cf) {
  return addFurtherCov(-1, cf);
}


void addTypeFct(type_fct TypeFct) {
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
  C->TypeFct = TypeFct;
  //  printf("Typefct: %s\n", C->name);
}


//int addFurtherCov(nonstat_covfct cf, covfct D, covfct D2) {
//  return addFurtherCov(-1, cf, D, D2);
//}

void nablahess(covfct nabla, covfct hess) {
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr

  assert((nr>=0) && (nr<currentNrCov));
  assert(C->vdim == SCALAR || C->vdim == SUBMODEL_DEP || 
	 C->vdim == PARAM_DEP || nabla==NULL);
  assert(C->cov != NULL && nabla!=NULL && hess != NULL);
  
  C->nabla=nabla;    
  C->hess = hess;

}


void addLocal(getlocalparam coinit, getlocalparam ieinit) {
  int nr = currentNrCov - 1;
  assert(nr>=0 && nr < currentNrCov) ;
  defn *C = DefList + nr; // nicht gatternr
  int *pref = C->pref;
  assert(C->D!=ErrD);
  if ((C->implemented[CircEmbedIntrinsic] = ieinit != NULL)) {
    assert(C->D2 != NULL);
    C->ieinit = ieinit;
    if (pref[CircEmbedIntrinsic] == PREF_NONE)
      pref[CircEmbedIntrinsic] = PREF_BEST;
  }
  if ((C->implemented[CircEmbedCutoff] = coinit != NULL)) {
    C->coinit = coinit;
    if (pref[CircEmbedCutoff] == PREF_NONE) pref[CircEmbedCutoff] = PREF_BEST;
    if (pref[CircEmbedIntrinsic] > 2)
      pref[CircEmbedIntrinsic] = 2;
  }
}

void addCallLocal(altlocalparam alt) {
  int nr = currentNrCov - 1;
  assert(nr>=0 && nr < currentNrCov) ;
  defn *C = DefList + nr; // nicht gatternr
  C->alternative = alt;
}
   
#define ASSERT_TBM  \
  /* // printf("assert_tbm %s %d %d %d %d\n", C->name, isPrevModelD(DOM(C->systems[0], 0)), anyVariant(isPosDef, C), anyVariant(isManifold, C), isXonly(DOM(C->systems[0], 0))); */ \
  assert(equalsPrevModelD(DOM(C->systems[0], 0)) ||			\
	 equalsSubModelD(DOM(C->systems[0], 0)) ||			\
	 equalsParamDepD(DOM(C->systems[0], 0)) ||			\
 	 ((anyVariant(isPosDef, C) || anyVariant(isManifold, C))	\
	  && equalsXonly(DOM(C->systems[0], 0))));			\
  /* // printf("%d %d %d\n", C->vdim == SCALAR, C->vdim == SUBMODEL_DEP,	C->vdim == PARAM_DEP) ;	*/ \
  assert(C->vdim == SCALAR || C->vdim == SUBMODEL_DEP ||		\
	 C->vdim == PARAM_DEP)					

int addTBM(covfct tbm2) { 
  // must be called always AFTER addCov !!!!
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov));
  defn *C = DefList + nr; // nicht gatternr
  ASSERT_TBM; 
  C->tbm2=tbm2;
  if (tbm2 != NULL) {
    // addTBM is called from the other addTBM's -- so tbm2 might
    // be NULL
    assert(isSpaceIsotropic(C->systems[0]) ||isPrevModelI(C) || isSubModelI(C));
    assert(C->D != ErrD);
    C->implemented[TBM] = IMPLEMENTED;
    if (C->pref[TBM] == PREF_NONE) C->pref[TBM] = PREF_BEST;
  }
  // IMPLEMENTED must imply the NUM_APPROX to simplify the choice
  // between TBM2 and Tbm2Num
  return nr;
}

void addTBM(covfct tbm2, initfct Init, spectral_do spectralDo) {
  int nr = addTBM(tbm2);
  defn *C = DefList + nr; // nicht gatternr
  ASSERT_TBM;
  C->spectral=spectralDo;
  C->Init=Init;
  C->implemented[SpectralTBM] = true;
  if (C->pref[SpectralTBM] == PREF_NONE) C->pref[SpectralTBM] = PREF_BEST;
  }

void addTBM(initfct Init, spectral_do spectralDo) {
  addTBM((covfct) NULL, Init, spectralDo);
}
	

void addSpecific(int cov, bool copy) {
  int nr = currentNrCov - 1;
  assert((nr>=0) && (nr<currentNrCov));
  defn *X = DefList + cov,
    *C = DefList + nr; // nicht gatternr
  assert(//(SYSTYPE(DefList[nr].systems[0], 0) == ProcessType &&
	  // DefList[nr].kappas == X->kappas) || 
	 (DefList[nr].kappas == X->kappas || 
	 ( X->check == checkmal && DefList[nr].check==checkmultproc)  ||
	  ( X->check == checktrend && DefList[nr].check==checkTrendproc)) );
  // to do ... und die namen sollten auch gleich sein...

  if (copy) {
    int kappas = X->kappas;
    if (kappas == C->kappas) {
      for (int i=0; i<kappas; i++){
	STRCPY(C->kappanames[i], X->kappanames[i]);
	C->kappatype[i] = X->kappatype[i];
	C->sortof_tab[i] = X->sortof_tab[i];
	C->kappaParamType[i] = X->kappaParamType[i];
	C->kappaParamTypeNames[i] = X->kappaParamTypeNames[i];
      }
    }
    
    int maxsub = X->maxsub;
    if (maxsub == C->maxsub) {
      for (int i=0; i<maxsub; i++) {
	C->subintern[i] = X->subintern[i];
	STRCPY(C->subnames[i], X->subnames[i]);
      }
    }
  } else {
    make_internal();    
  }

  nickname(X->nick + STRLEN(CAT_TYPE_NAMES[SYSTYPE(X->systems[0], 0)]));

  while (true) {
    assert(X->Specific == MISMATCH ||X->Specific == UNSET || X->Specific == nr);
    X->Specific = nr;
    if (X->pref[Specific] == PREF_NONE) X->pref[Specific] = PREF_BEST;
    X->implemented[Specific] = IMPLEMENTED;
    //    printf("addspecific done %s\n", X->name);
    X++;
    //    printf("addspecific next %s\n", X->name);
    if (X->name[0] != InternalName[0]) break;
  }
}

void addSpecific(int cov) { addSpecific(cov, true); }

	
void addHyper(hyper_pp_fct hyper_pp) {
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
  assert((nr>=0) && (nr<currentNrCov));
  C->hyperplane=hyper_pp;
  C->implemented[Hyperplane] = hyper_pp!=NULL;
  if (C->pref[Hyperplane] == PREF_NONE) C->pref[Hyperplane] = PREF_BEST;
}
		   
//void addSpecialMeth(initstandard initspecial, dometh special)  {
///  int nr = currentNrCov - 1;
//  defn *C = DefList + nr; // nicht gatternr
//  C->initspecial=initspecial;
//  C->special=special;
//  if ((special!=NULL) || (initspecial!=NULL)) 
//    assert((special!=NULL) && (initspecial!=NULL));
//  C->implemented[Special] = true;
//}


void RandomShape(int maxmoments, structfct Struct, initfct Init, 
		 dofct Do, do_random_fct DoRandom, 
		 bool average, bool randomcoin, bool self_specific){
 // rein auf shape-Function bezogen, ohne Kenntnis von irgendeinem
  // Kovarianzmodell
  //
  // init and do der elementaren shape-Funktion
  // insbesondere notwendig wenn die elementare shape-Function zufaellig ist.
  // Erst dann koennen werte abgefragt werden
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr

  assert(Init != NULL && Do != NULL && DoRandom != NULL && Struct!=NULL);

  C->maxmoments = maxmoments;
  C->Struct = Struct;
  C->Do = Do;
  C->DoRandom = DoRandom;
  C->Init = Init;
  C->implemented[Average] = average; 
  C->implemented[RandomCoin] = randomcoin; 
  C->implemented[Specific] = self_specific;
  if (self_specific) {
    if (C->Specific == MISMATCH ||  C->Specific == UNSET) C->Specific = nr;
    else {
      assert( C->name[0] == InternalName[0]);
    } 
  } else {
    Types type = SYSTYPE(DefList[nr].systems[0], 0);
    C->Specific = isProcess(type) || isInterface(type) ? MISMATCH : UNSET;
  }
}

void RandomShape(int maxmoments, structfct Struct, initfct Init, dofct Do,
		 bool average, bool randomcoin, bool self_specific){
 RandomShape(maxmoments, Struct, Init, Do, do_random_failed,
	     average, randomcoin, self_specific);
}

void RandomShape(int maxmoments, structfct Struct, initfct Init,
		 dofct Do){
  RandomShape(maxmoments, Struct, Init, Do, do_random_failed, 
	      false, false, false);
}

void RandomShape(int maxmoments, structfct Struct, initfct Init,
		 dofct Do, finaldofct FinalDo){
  RandomShape(maxmoments, Struct, Init, Do, do_random_failed, 
	      false, false, false);
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
  C->FinalDo = FinalDo;
}

void RandomShape(int maxmoments, structfct Struct, initfct Init, 
		 do_random_fct DoRandom){
  RandomShape(maxmoments, Struct, Init, do_failed, DoRandom, 
	      false, false, false);
}

void RandomShape(int maxmoments, initfct Init, dofct Do, 
		 bool average){
  RandomShape(maxmoments, structOK, Init, Do, do_random_failed,
	      average, !average, false);
}

void RandomShape(int maxmoments, initfct init, dofct Do){
   RandomShape(maxmoments, structOK, init, Do, do_random_failed, 
 	      false, true, false); 
}

void RandomShape(structfct Struct, bool average){
   RandomShape(-1, Struct, initOK, doOK,  do_random_failed, 
	      average, !average, false);
}

void RandomShape(structfct Struct){
  RandomShape(-1, Struct, initOK, doOK, do_random_failed, false, true, false);
}


void addSpecial(minmaxfct minmaxeigen){
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
  C->minmaxeigenvalue = minmaxeigen;
}

void addGaussMixture(draw_random drawmix,
		     log_mixdens logmixdens) {
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
  assert(drawmix != NULL && logmixdens != NULL);
  C->drawmix = drawmix;
  C->logmixdens = logmixdens;
}

void addReturns(//return_fct Covariance, ext_bool_ret_fct isCovariance, 
		return_covmat CovMatrix, ext_bool_ret_fct isCovMatrix
		//tworeturns_fct InverseCovMatrix,
		//		ext_bool_ret_fct isInverseCovMatrix,
		//return_fct Variogram, ext_bool_ret_fct isVariogram,
		//return_fct PseudoVariogram, ext_bool_ret_fct isPseudoVariogram
		) {
  int nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr
  /*
  if (Covariance!=NULL) {
    C->covariance=Covariance;
    C->is_covariance = isCovariance==NULL ? isFalse : isCovariance;
  } else assert(isCovariance==NULL);
  */
  if (CovMatrix!=NULL) {
    C->covmatrix=CovMatrix;
    C->is_covmatrix = isCovMatrix==NULL ? isFalse : isCovMatrix;
  } else assert(isCovMatrix==NULL);
  /*
  if (InverseCovMatrix!=NULL) {
    C->inversecovmatrix=InverseCovMatrix;
    C->is_inversecovmatrix =
      isInverseCovMatrix==NULL ? isFalse : isInverseCovMatrix;
  } else assert(isInverseCovMatrix==NULL);
  if (Variogram!=NULL) {
    C->variogram=Variogram;
    C->is_variogram = isVariogram==NULL ? isFalse : isVariogram;
  } else assert(isVariogram==NULL);
  if (PseudoVariogram!=NULL) {
    C->pseudovariogram=PseudoVariogram;
    C->is_pseudovariogram =
      isPseudoVariogram==NULL ? isFalse : isPseudoVariogram;
  } else assert(isPseudoVariogram==NULL);
  */
}



void TailTaylor(double t, double tpow, double texpc, double texppow) {
  int
    nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr

  if (C->finiterange == wahr) {
    assert(t == 0 && tpow==0 && texpc==0 && texppow==0);
  }


  // TaylorConst * x^TaylorPow * EXP(-TaylorExpConst * x^TaylorExpPow)
  C->TailN = 0;
  C->Tail[C->TailN][TaylorConst] = t;
  C->Tail[C->TailN][TaylorPow] = tpow;
  C->Tail[C->TailN][TaylorExpConst] = texpc;
  C->Tail[C->TailN][TaylorExpPow] = texppow;
  C->TailN++;

  assert(C->TailN <= MAXTAYLOR);
}


void Taylor(double c, double pow) {
  int 
    nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr

  // TaylorConst * x^TaylorPow

  C->TaylorN = 0;
  Types type = SYSTYPE(C->systems[0], 0);
  if (isPosDef(type) || isManifold(type)) {
    assert(pow != 0.0);
    C->Taylor[C->TaylorN][TaylorConst] = 1.0;
    C->Taylor[C->TaylorN][TaylorPow] = 0.0;    
    C->TaylorN++;
  }

  C->Taylor[C->TaylorN][TaylorConst] = c;
  C->Taylor[C->TaylorN][TaylorPow] = pow;
  C->TaylorN++;
  assert(C->TaylorN <= MAXTAYLOR);

  if (C->finiterange == wahr) TailTaylor(0, 0, 0, 0);

}


void Taylor(double c, double pow, double c1, double pow1) {
  int  
    nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr

  Taylor(c, pow);
  C->Taylor[C->TaylorN][TaylorConst] = c1;
  C->Taylor[C->TaylorN][TaylorPow] = pow1;
  C->TaylorN++;
  assert(C->TaylorN <= MAXTAYLOR);
}

void Taylor(double c, double pow, double c1, double pow1, 
	    double c2, double pow2) {
  int 
    nr = currentNrCov - 1;
  defn *C = DefList + nr; // nicht gatternr

  Taylor(c, pow, c1, pow1);
  C->Taylor[C->TaylorN][TaylorConst] = c2;
  C->Taylor[C->TaylorN][TaylorPow] = pow2;
  C->TaylorN++;

  assert(C->TaylorN <= MAXTAYLOR);

}


void setptwise(ptwise_type pt) {
  int
    nr = currentNrCov - 1;
  DefList[nr].ptwise_definite = pt;
}


void setDI(allowedD_fct D, allowedI_fct I, setDI_fct f) {
  int nr = currentNrCov - 1;
  assert(DefList[nr].Dallowed == NULL ||
	 DefList[nr].Dallowed == allowedDstandard);
  assert(DefList[nr].Iallowed == NULL ||
	 DefList[nr].Iallowed == allowedIstandard ||
	 DefList[nr].Iallowed == allowedPrevModelI);
  assert(D != NULL || I != NULL);
  if (D!=NULL) DefList[nr].Dallowed = D;
  if (I!=NULL) DefList[nr].Iallowed = I;
  DefList[nr].setDI = f;
}

