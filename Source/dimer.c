/******************************************************

	dimer.c
	
	This is the shell of the minimization program.  It 
	reads the z-matrix into memory, makes cartesians
	and everything else needed and then minimizes.

*********************************************************/

#include <stdlib.h>
#include "dimer.h"
#include "potentials.h"

int main(void)
{
	FILE	*outfile;

	FILE	*theFile;
	FILE	*theGraFile;
	short	numAtoms;
	zm_el	*theMatrix;
	double	*theResults;
	vector	*theCarts;
	bonds	*theBonds;
	del_el	*theGradTable;
	double	*theQs;
	
	double	theEnergy;
	short	gradCount;
	
/* Old Mac Think C library functions. Assume that everyting will go to stdout */
/*
	cshow(stdout);
	
	cecho2file("zmat.out", 0, stdout);
	cecho2printer(stdout);
*/
	numAtoms=countAtoms(theFile);

	printf("there are %hd atoms.\n\n", numAtoms);

	theMatrix=(zm_el *)(calloc((size_t)(numAtoms),sizeof(zm_el)));
	theCarts=(vector *)calloc((size_t)numAtoms, sizeof(vector));
	theBonds=(bonds *)calloc((size_t)numAtoms, sizeof(bonds));
	theGradTable=(del_el *)calloc((size_t)(numAtoms-1), sizeof(del_el));
	theQs=(double *)calloc((size_t)numAtoms, sizeof(double));
	theResults=(double *)(calloc((size_t)(gradCount),sizeof(double)));
		
	gradCount=makeGradTable(theGradTable, theGraFile);

	readZmat(theFile, theMatrix, theBonds, theCarts, numAtoms);
	printMatrix(theMatrix, numAtoms);
	
	makeCartesian(theMatrix, theBonds, theCarts, numAtoms);	
	printCarts(theCarts, numAtoms);

	makeQtable(theMatrix, theQs, numAtoms);
		
	theEnergy=energy_tip(theQs, theCarts, numAtoms);
	printf("The energy is ..... %3.2E\n", theEnergy);

/*	exit(0);*/

	printf("\n\n------COMMENCING MINIMIZATION-------\n\n");
	
	theEnergy=Newtons(theMatrix, theGradTable, theCarts, theQs, theBonds,  
			numAtoms,  gradCount,  theResults);

	printf("The minimized energy is ..... %3.2E\n\n", theEnergy);
	printMatrix(theMatrix, numAtoms);
	printCarts(theCarts, numAtoms);
	return 0;
}

void printMatrix(zm_el *theMatrix, short numAtoms)
{
	short index;
	
	for (index=0;index<numAtoms;index++)
	{
		printf("%hd\t%c\t%hd\t%3.2e\t%hd\t%3.2e\t%hd\t%3.2e\n", index,
			theMatrix[index].a,
			theMatrix[index].b,
			theMatrix[index].ab,
			theMatrix[index].c,
			theMatrix[index].abc,
			theMatrix[index].d,
			theMatrix[index].abc_bcd);
	}
	printf ("\n");	
}

void printCarts(vector *theCarts, short numAtoms)
{
	short index;
	
	for (index=0;index<numAtoms;index++)
	{
		printf("%hd\t%3.2e\t%3.2e\t%3.2e\n", index,
			theCarts[index][Vx],
			theCarts[index][Vy],
			theCarts[index][Vz]);
	}
	printf ("\n");	
}
