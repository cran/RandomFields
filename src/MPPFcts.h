#ifndef MppFcts_H
#define MppFcts_H 1
extern double MPP_APPROXZERO;

void cone_init(mpp_storage *s, int dim, param_type param);
void cone(mpp_storage *s, double *min, double *max,
	  mppmodel *model );
int checkcone(double *param, int timespacedim, SimulationType method);
SimulationType methodcone(int spacedim, bool grid);
void rangecone(int spatialdim, int *index, double* range);
void infocone(double *p, int *maxdim, int *CEbadlybehaved);

void gaussmpp_init(mpp_storage *s, int dim, param_type param);
void gaussmpp(mpp_storage *s, double *min, double *max,
	  mppmodel *model );

void circular_init(mpp_storage *s, int dim, param_type param);
void circularMpp(mpp_storage *s, double *min, double *max,
	  mppmodel *model );

void spherical_init(mpp_storage *s, int dim, param_type param);
void sphericalMpp(mpp_storage *s, double *min, double *max,
	  mppmodel *model );

#endif /* MppFcts_H */
 
