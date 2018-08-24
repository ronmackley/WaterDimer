/************************************************

	stdDefines.h---miscelnaeous handy things
	
*************************************************/

#pragma once

#define TRUE 1
#define FALSE 0
#define LINELEN 80

#define mark_of_the_Beast 666

#ifndef pi
	#define pi 3.14159265359
#endif

#define deg2rad 0.017453293

#define rad2deg 57.29577951

typedef struct
{
	short 	theCount;
	short	toWhere[8];
} bonds;

typedef struct
{
	short	ab;
	short	abc;
	short	abc_bcd;
} del_el;

typedef struct
{
	short	a;
	short	b;
	double	ab;
	short	c;
	double	abc;
	short	d;
	double	abc_bcd;
} zm_el;
