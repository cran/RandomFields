#ifndef MppFcts_H
#define MppFcts_H 1
extern double MPP_APPROXZERO;

void cone_init(mpp_storage *s, int v);
void cone(mpp_storage *s, int v, double *min, double *max,
	  mppmodel *model );
int checkcone(double *param, int timespacedim, SimulationType method);
SimulationType methodcone(int spacedim, bool grid);
void rangecone(int spatialdim, int *index, double* range);
void infocone(double *p, int *maxdim, int *CEbadlybehaved);

void gaussmpp_init(mpp_storage *s, int v);
void gaussmpp(mpp_storage *s, int v, double *min, double *max,
	  mppmodel *model );

void circular_init(mpp_storage *s, int v);
void circularMpp(mpp_storage *s, int v, double *min, double *max,
	  mppmodel *model );

void spherical_init(mpp_storage *s, int v);
void sphericalMpp(mpp_storage *s, int v, double *min, double *max,
	  mppmodel *model );

#endif /* MppFcts_H */
 
