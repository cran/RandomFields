#ifndef MppFcts_H
#define MppFcts_H 1
extern Real MPP_APPROXZERO;

void cone_init(key_type *key);
void cone(mpp_storage *bs, Real *min, Real *max,
	  mppmodel *model END_WITH_RANDOM);
int checkcone(key_type *key);

void gaussmpp_init(key_type *key);
void gaussmpp(mpp_storage *bs, Real *min, Real *max,
	  mppmodel *model END_WITH_RANDOM);

void circular_init(key_type *key);
void circularMpp(mpp_storage *bs, Real *min, Real *max,
	  mppmodel *model END_WITH_RANDOM);

void spherical_init(key_type *key);
void sphericalMpp(mpp_storage *bs, Real *min, Real *max,
	  mppmodel *model END_WITH_RANDOM);

#endif /* MppFcts_H */
 
