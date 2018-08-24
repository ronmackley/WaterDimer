/******************************************************************

vector_prime.c --- the source file to work with vectors.

******************************************************************/

#include "vector.h"

void v_parallel(double a, vector v, vector result)
{
	/* returns a vector of magnitude a parallel with vector v */

	double Mv, coeff;
	
	Mv=v_mag(v);
	coeff=sqrt(a/Mv);
	result[Vx]=coeff*(v[Vx]/Mv);
	result[Vy]=coeff*(v[Vx]/Mv);
	result[Vz]=coeff*(v[Vx]/Mv);
	
}

void v_cross(vector vec1, vector vec2, vector result)
{

	result[Vx]=vec1[Vy]*vec2[Vz]-vec1[Vz]*vec2[Vy];
	result[Vy]=vec1[Vx]*vec2[Vz]-vec1[Vz]*vec2[Vx];
	result[Vz]=vec1[Vx]*vec2[Vy]-vec1[Vy]*vec2[Vx];
	
}


double v_dot(vector vec1, vector vec2)

{
	double result;
	
	result=vec1[Vx]*vec2[Vx]+vec1[Vy]*vec2[Vy]+vec1[Vz]*vec2[Vz];
	
	return(result);
}

void v_cart2sph(vector thevec)
{

	/*
		returns:     x=r
		             y=theta	angle in xy-plane
		             z=phi		angle from z-axis
	*/
	vector result;
			
	result[Vr]=v_mag(thevec);
	result[Vtheta]=atan2(thevec[Vy], thevec[Vx]);
	result[Vphi]=acos(thevec[Vz]/v_mag(thevec));
	v_copy(result, thevec);
	
}

void v_sph2cart(vector thevec)
{

	/*
		when called,   x=r
		               y=theta	angle in xy-plane
		               z=phi	angle from z-axis
	*/
	vector result;
		
	result[Vx]=thevec[Vr]*sin(thevec[Vphi])*cos(thevec[Vtheta]);
	result[Vy]=thevec[Vr]*sin(thevec[Vphi])*sin(thevec[Vtheta]);
	result[Vz]=thevec[Vr]*cos(thevec[Vphi]);
	v_copy(result, thevec);
		
}

void v_cyl2cart(vector thevec)
{

	/*
		in thevec,     r    =x
		               theta=y	angle in xy-plane
		               z	=z
	*/
	vector result;
	
	result[Vx]=thevec[Vr]*cos(thevec[Vtheta]);
	result[Vy]=thevec[Vr]*sin(thevec[Vtheta]);
	result[Vz]=thevec[Vz];
	v_copy(result, thevec);

}

void v_cart2cyl(vector thevec)
{
	
	/*
		when returned, x=r
		               y=theta	angle in xy-plane
		               z=z
	*/
	vector result;
	
	result[Vr]=sqrt(thevec[Vx]*thevec[Vx]+thevec[Vy]*thevec[Vy]);
	result[Vtheta]=atan2(thevec[Vy], thevec[Vx]);
	result[Vz]=thevec[Vz];
	v_copy(result, thevec);

}


void v_rotate (vector thevec, double theta, double phi)
{ 
	v_cart2sph(thevec);
	thevec[Vtheta]+=theta; 
	thevec[Vphi]+=phi;
	v_sph2cart(thevec);
} 
