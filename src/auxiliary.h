
#ifndef AUXILIARY_H
#define AUXILIARY_H 1

#include "basic.h"
#include "AutoRandomFields.h"

double I0mL0(double x);
extern "C" {
  SEXP vectordist(SEXP V, SEXP diag); 
  void sleepMilli(int *milli);
  void I0ML0(double *x, int *n);
  double struve(double x, double nu,  double factor_sign, bool expscaled);
  void StruveH(double *x, double *nu);
  void StruveL(double *x, double *nu, int * expScaled);
  void ordering(double *d, int len, int dim, int *pos);
  void Ordering(double *d, int *len, int *dim, int *pos);
}
bool is_diag(double *aniso, int dim);


#define INT Integer(el, name, 0)
#define LOG Logical(el, name, 0)
#define NUM Real(el, name, 0)
#define CHR Char(el, name)
#define STR(X, N)  strcopyN(X, CHAR(STRING_ELT(el, 0)), N);
#define POS0INT NonNegInteger(el, name) /* better: non-negative */
#define POS0NUM NonNegReal(el, name)
#define NEG0NUM NonPosReal(el, name)
#define POSINT PositiveInteger(el, name) /* better: non-negative */
#define POSNUM PositiveReal(el, name)



double Real(SEXP p, char *name, int idx);
void Real(SEXP el,  char *name, double *vec, int maxn) ;
int Integer(SEXP p, char *name, int idx, bool nulltoNA) ;
int Integer(SEXP p, char *name, int idx);
void Integer(SEXP el, char *name, int *vec, int maxn) ;
void Integer2(SEXP el, char *name, int *vec) ;
bool Logical(SEXP p, char *name, int idx);
char Char(SEXP el, char *name) ;
void String(SEXP el, char *name, char names[MAXUNITS][MAXCHAR]) ;
double NonNegInteger(SEXP el, char *name) ;
double NonNegReal(SEXP el, char *name) ;
double NonPosReal(SEXP el, char *name) ;
double PositiveInteger(SEXP el, char *name) ;
double PositiveReal(SEXP el, char *name) ;
int GetName(SEXP el, char *name, const char * List[], int n) ;
 int GetName(SEXP el, char *name, const char * List[], int n,
	    int defaultvalue) ;





#endif /* AUXILIARY_H */




