#ifndef hana_public_H
#define hana_public_H 1
#include "GSLvsR.h"

EXTERN void addownfunctions(); // to be called by .C("addownfunctions") 
//                           at the very beginning -- after library(RandomFields)

EXTERN void SetMinDiam(double *x); 

#endif /* hana_public_H */
