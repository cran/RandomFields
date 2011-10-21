#ifndef MPPstandard
#define MPPstandard 1


#define STANDARDDEFINITIONEN(COV)	\
  globalparam *gp = meth->gp;\
  location_type *loc = meth->loc;\
  mpp_storage *s = (mpp_storage*) meth->S;\
  cov_model *cov = COV; \
  cov_fct *C = CovList + cov->nr;\
  mpp_model randomcoin = C->randomcoin;\
  \
  long zaehler= 0,\
    cumgridlen[MAXMPPDIM +1],\
    total = loc->totalpoints;\
  int d, gridlen[MAXMPPDIM], n,\
    every = gp->general.every,\
    nthreshold = (every>0) ? every : MAXINT,\
    deltathresh = nthreshold;\
  res_type factor; \
  double Radius,  \
    inc[MAXMPPDIM],  \
    *min = s->min, \
    *u = s->u\
 
#define STANDARDINIT\
  for (zaehler=0; zaehler<total; res[zaehler++]=0.0);\
  if (meth->type <= TypeDiag && loc->grid) {\
    cumgridlen[0] = 1;\
    for (d=0; d<dim; d++) {\
      inc[d] = meth->grani[d][XSTEP];\
      gridlen[d] =  loc->length[d];\
      cumgridlen[d+1] = gridlen[d] * cumgridlen[d];\
    }\
  }

#define STANDARDGRID\
  int end[MAXMPPDIM], start[MAXMPPDIM], delta[MAXMPPDIM], nx[MAXMPPDIM];\
  double xstart[MAXMPPDIM], x[MAXMPPDIM];\
  bool inside=true;\
  \
  zaehler = 0;\
  for (d=0; d<dim; d++) {\
    start[d] = (int) ((u[d] - Radius - min[d]) / inc[d]);\
    if (start[d] < 0) start[d]=0;\
    end[d] = 1 + (int) ((u[d] + Radius - min[d]) / inc[d]);\
    if (end[d] > gridlen[d]) end[d] = gridlen[d];\
    if (start[d] >= end[d]) inside = false;\
    delta[d] = (end[d] - start[d]) * cumgridlen[d];\
    nx[d] = start[d];\
    zaehler += start[d] * cumgridlen[d];\
    x[d] = xstart[d] = min[d] + (double) start[d] * inc[d]; \
  }

#define STANDARDINKREMENT\
  d = 0;\
  nx[d]++;\
  x[d] += inc[d];\
  zaehler += cumgridlen[d];\
  \
  while (nx[d] >= end[d]) {\
    nx[d] = start[d];\
    x[d] = xstart[d];\
    zaehler -= delta[d]; /* delta is positive */ \
    if (++d >= dim) break;\
    nx[d]++;\
    x[d] += inc[d];\
    zaehler += cumgridlen[d];\
  }\
  if (d >= dim) break;\

#define STANDARDUSER\
  if (n >= nthreshold) {\
 /*   PRINTF("\b\b\b\b%d%%", (int) (n * 100 / ntot)); */ \
    PRINTF("%7d %3d%%\n", n, (int) (n * 100 / ntot)); \
    nthreshold += deltathresh;\
  } \
 R_CheckUserInterrupt();\

#endif
