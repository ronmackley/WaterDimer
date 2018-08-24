/******************************************************


	makeData.c---uses the z-matrix routines etc to
	calculate a table of energies for various con-
	figurations of atoms.  It varies the length of
	the H-bond and the dihedral angle along the h-
	bond axis and holds everything else constant as
	per the users input

*********************************************************/

#include <time.h>

#include "makeData.h"

#define DIV 20

int main(void)
{
	FILE	*outFile;
	FILE	*theFile;
	short	numAtoms;
	zm_el	*theMatrix;
	double	*theResults;
	vector	*theCarts;
	bonds	*theBonds;
	del_el	*theGradTable;
	double	*theQs;
	
	char 	FileName[80];
	struct	tm *date;
	
	time_t	now;
	double	theEnergy, t1, t2, results[DIV+1], x[DIV+1], y[DIV+1];
	double	S1, S2, S3, S4;
	short	i1, i2;
	
//	cshow(stdout);	/*show the console window--Mac Specific*/
	
	theMatrix=(zm_el *)(calloc((size_t)(numAtoms),sizeof(zm_el)));
	theCarts=(vector *)calloc((size_t)numAtoms, sizeof(vector));
	theBonds=(bonds *)calloc((size_t)numAtoms, sizeof(bonds));
	theQs=(double *)calloc((size_t)numAtoms, sizeof(double));

	numAtoms=countAtoms(theFile);
	readZmat(theFile, theMatrix, theBonds, theCarts, numAtoms);

	makeCartesian(theMatrix, theBonds, theCarts, numAtoms);	
	makeQtable(theMatrix, theQs, numAtoms);

/*make a unique file name based on the date */

	now=time(NULL);
	date=localtime(&now);
	strftime(FileName, 25, "outfile%m%d%y%H%M%S", date);

	outFile=fopen(FileName, "w+");

/*make the arrays of the steps for x and y */

/* normal ranges 0.0=4.0 Å for the bond length 
							and 0-2π for the angle*/
							
	for (i1=0,t1=.01;i1<=DIV&&t1<=4.2;i1++,t1+=0.2)
		{
			x[i1]=t1;
		}
	for (i1=0,t1=0.0;i1<=DIV&&t1<=6.28;i1++,t1+=0.314)
		{
			y[i1]=t1;
		}


	printf("H2a-O1-m1:  ");
	scanf("%lf", &S1);
	printf("H2a-O1-m1-H1b:  ");
	scanf("%lf", &S2);
	printf("O2-H2a-O1:  ");
	scanf("%lf", &S3);
	printf("m2-O2-H2a-O1:  ");
	scanf("%lf", &S4);

	theMatrix[4].abc=S1*deg2rad;		/*H2a-O1-m1*/
	theMatrix[4].abc_bcd=S2*deg2rad;	/*H2a-O1-m1-H1b*/
	theMatrix[5].abc=S3*deg2rad;		/*O2-H2a-O1*/
	theMatrix[6].abc_bcd=S4*deg2rad;	/*m2-O2-H2a-O1*/

	fprintf(outFile, "Serial Number:  %s\n\n", FileName);
	
	fprintf(outFile, "H2a-O1-m1\t%3.2f\n", theMatrix[4].abc*rad2deg);
	fprintf(outFile, "H2a-O1-m1-H1b\t%3.2f\n", theMatrix[4].abc_bcd*rad2deg);
	fprintf(outFile, "O2-H2a-O1\t%3.2f\n", theMatrix[5].abc*rad2deg);
	fprintf(outFile, "m2-O2-H2a-O1\t%3.2f\n", theMatrix[6].abc_bcd*rad2deg);


	fprintf(outFile, "\n\n\t");
	
	for (i1=0;i1<=DIV;i1++)
		fprintf(outFile, "%3.2f\t", y[i1]);
	fprintf(outFile, "\n");
	
	for (i1=0;i1<=DIV;i1++)
	{
		theMatrix[4].ab=x[i1];
		fprintf(outFile, "%4.3f\t", x[i1]);
		for (i2=0;i2<DIV;i2++)
		{
			theMatrix[5].abc_bcd=y[i2];
			makeCartesian(theMatrix, theBonds, theCarts, numAtoms);	
			theEnergy=energy_tip(theQs, theCarts, numAtoms);
			fprintf(outFile, "%4.3e\t", theEnergy);
		}
		fprintf(outFile, "\n");
		printf("Done a Row!!!\n");
	}
					
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
