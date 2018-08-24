/************************************************************************

vectors.h---the header for vectors.c

This tile declares the struct vector and prototypes functions to deal
with them, such as adding, subtracting, cross and dot products, rotating,
and the magnitude of a vector.

***********************************************************************/

#pragma once
#include <math.h>
	
/* define the vector subscripts for the x, y, and z components
   as x, y, and z are popular variable neames, Letterman has 
   torn the V off of his sweatshirt and prefixed it to stand
   for vector */

#define VECSIZE 3

#define Vx 0
#define Vy 1
#define Vz 2

#define Vr 0
#define Vtheta 1
#define Vphi 2

#ifndef pi
	#define pi 3.14159265359
#endif

/* macro defn's  */

#define v_copy(a,b) \
	{b[Vx]=a[Vx]; \
	b[Vy]=a[Vy]; \
	b[Vz]=a[Vz];}
	
#define v_mag(a) (sqrt(a[Vx]*a[Vx]+a[Vy]*a[Vy]+a[Vz]*a[Vz]))

#define v_inver(a) \
	{a[Vx]*=-1;\
	a[Vy]*=-1;\
	a[Vz]*=-1;}
	
#define v_add(a, b) \
	{a[Vx]+=b[Vx]; \
	a[Vy]+=b[Vy]; \
	a[Vz]+=b[Vz];}
	
#define v_sub(a, b) \
	{a[Vx]-=b[Vx]; \
	a[Vy]-=b[Vy]; \
	a[Vz]-=b[Vz];}


typedef double vector[VECSIZE];


/* here are ansi prototypes */

void v_parallel(double, vector, vector);

void v_cross(vector, vector, vector);
double v_dot(vector, vector);

void v_sph2cart(vector);
void v_cart2sph(vector);
void v_cyl2cart(vector);
void v_cart2cyl(vector);
void v_rotate (vector, double, double);

