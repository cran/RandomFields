


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 Felix Ballani

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


#ifndef POLYGON_H
#define POLYGON_H 1


/*******************************************************************************
 * PoissonPolygon.h
 ****/

// Definitions
//#define PI        3.141592653589793238462643383279

// Structs
typedef struct vertex
{
	double x[2];				// vertex coordinates
} vertex;

typedef struct edge
{
	double u[2];				// normal vector
	double p;						// distance from the origin
} edge;

typedef struct polygon
{
	int n;							// number of vertices = number of edges
        vertex *v;		// array of vertices
        edge *e;			// array of edges
	double box0[2];			// coordinates of the left lower vertex of smallest box containing the polygon
	double box1[2];			// coordinates of the right upper vertex of smallest box containing the polygon
} polygon;

// Functions
double scProd(const double *x, const double *y);
int compareAngles(const void * a, const void * b);
void rTriangle(double *phi);
int rPoissonPolygon(struct polygon *P, double lambda);
int rPoissonPolygon2(struct polygon *P, double lambda, bool);
void freePolygon(struct polygon *P);
int isInside(struct polygon *P, double *x);
double getArea(struct polygon *P);


#endif /* POLYGON_H */
