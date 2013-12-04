#ifndef CONVEX_H
#define CONVEX_H 1


/*******************************************************************************
 * convhull2D.h
 ****/
 

int ccw(double **P, int i, int j, int k);
int cmpl(const void *a, const void *b);
int cmph(const void *a, const void *b);
int make_chain(double** V, int n, int (*cmp)(const void*, const void*));
int ch2d(double **P, int n);

#endif /* CONVEX_H */
