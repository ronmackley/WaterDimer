/************************************************************

	zmat.h
	
	the header file for all code which copes with
	z-matricies.
	
************************************************************/

#pragma once

#include <stdio.h>
#include "vector.h"
#include "potentials.h"
					
#define COMMENTMARK ';' 	/*this is the char which appears at the */
							/* beginning of a line to delimit comments*/
#define ENDOFMATRIX '@'
#define LINELEN 80			/* the maximum line length*/

#define Hid 'H'
#define Mid 'M'
#define Oid 'O'

#define qH 0.52
#define qM -1.04
#define qDefault 0.0

// typedef struct
// {
// 	short 	theCount;
// 	short	toWhere[8];
// } bonds;

/*
typedef struct zm_el
{
	short	a;
	short	b;
	double	ab;
	short	c;
	double	abc;
	short	d;
	double	abc_bcd;
} zm_el;
*/

/* ansi prototypes follow */

void notComment(FILE *, char *);
short countAtoms(FILE *);
short readZmat(FILE *, zm_el *, bonds *, vector *, short);
void makeCartesian(zm_el *, bonds *, vector *, short);
short makeZmat(zm_el *, bonds *, vector *, short);
double bondAng(vector *, short, short, short);
double dihedralAng(vector *, short, short, short, short);
void addBond(bonds *, short, short);
void makeQtable(zm_el *, double *, short);
short isBonded(bonds, short);