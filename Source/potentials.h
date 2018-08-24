/***************************************************************************

	potentials.h----the header for potentials.c

	
	this file contains all the definitions for the potential functions
	for the forces between the various atoms for TIP4P water simulation
**************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "stdDefines.h"
#include "vector.h"
				
#pragma once

/*constants from _Physical Chemistry_ by Atkins, 4th ed.*/

#define e 1.6602177 /*Coulombs  */
#define e2 2.566971139 /* e squared */

#define N0 6.02214e23
#define FourPiEpsilon0 1.11265e-10

#define A 2.5104 /*kJ Å^12 mol^-1  */
#define B 2552.24 /*kJ Å^6 mol^-1  */

#define theGradName "tip4p.gra"
#define ENDOFGRAMARK '@'

/* data structures */
/*
typedef dm_el;
typedef struct
{
	short	ab;
	short	abc;
	short	abc_bcd;
} del_el;
*/

/* the del_el is a wonderful creation of mine to allow me to control the 
	constraints on a molecule by deliberately specifying what the delta value 
	used to find the derivative will be.  In so doing, I can set the delta for
	those bond lengths etc. which I don't want to change to FALSE and they will
	remain constant as the gradient routine checks the value of del and finds
	the gradient only if it is NOT FALSE */

/*  prototypes */

double electro(vector, double, vector, double);
double leonard_jones (vector, vector);
double energy_tip(double *, vector *, short);
double gradient_tip(zm_el *, del_el *, vector *, double *, bonds *, short, double *);
short makeGradTable(del_el *, FILE *);
