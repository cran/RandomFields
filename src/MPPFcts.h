#ifndef MppFcts_H
#define MppFcts_H 1
extern Real MPP_APPROXZERO;

void cone_init(mpp_storage *s, int v);
void cone(mpp_storage *s, int v, Real *min, Real *max,
	  mppmodel *model );
int checkcone(Real *param, int timespacedim, SimulationType method);
SimulationType methodcone(int spacedim, bool grid);
void rangecone(int spatialdim, int *index, Real* range);

void gaussmpp_init(mpp_storage *s, int v);
void gaussmpp(mpp_storage *s, int v, Real *min, Real *max,
	  mppmodel *model );

void circular_init(mpp_storage *s, int v);
void circularMpp(mpp_storage *s, int v, Real *min, Real *max,
	  mppmodel *model );

void spherical_init(mpp_storage *s, int v);
void sphericalMpp(mpp_storage *s, int v, Real *min, Real *max,
	  mppmodel *model );

#endif /* MppFcts_H */
 
