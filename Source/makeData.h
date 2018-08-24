/************************************************************

	makeData.h
	
	header file for the makeData program.

*************************************************************/

#pragma once

#include "zmat.h"		
#include "vector.h"
		
#define theName "tip4p.zma"
		
/* ansi prototypes */

int main(void);
void printMatrix(zm_el *, short);
void printCarts(vector *, short);