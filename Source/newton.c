/*************************************************************************

	newton.c
	
	C code for newtons methods stuff

*************************************************************************/

#include "newton.h"

double Newtons(zm_el *theMatrix, del_el *theGradTable, vector *theCarts,
		double *theQs, bonds *theBonds, short numAtoms, short gradCount, double *theResults)
{
	double epsilon, theEnergy, lastEnergy;
	short done=FALSE, itercnt=0;
		
	lastEnergy=energy_tip(theQs, theCarts, numAtoms);

	printf("energy\t\t\tlast energy\t\tquotient\n------\t\t\t-----------\t\t------\n");
	while(done==FALSE)
	{
		itercnt++;
		theEnergy=gradient_tip(theMatrix, theGradTable, theCarts, theQs,
				theBonds, numAtoms, theResults);
		epsilon=fabs(theEnergy/lastEnergy);
		if((epsilon<1e-4)||(itercnt>25))
			done=TRUE;
		update(theResults, theMatrix, theGradTable, theEnergy, 
						numAtoms, gradCount);
		makeCartesian(theMatrix, theBonds, theCarts, numAtoms);	
		printf("%5.4e\t\t%5.4e\t\t%5.4e\n", theEnergy, lastEnergy, epsilon);
		lastEnergy=theEnergy;
	}
	printf ("Minimization done in %hi iterations\n\n", itercnt);
	return lastEnergy;
}

void update(double *theResults, zm_el *theMatrix, del_el *theGradTable,
				double theEnergy, short numAtoms, short gradCount)
{
	short index, gCount=0;
	double	*dummy;
	
/* commented out as it's only used for debugging */
/*	printMatrix(theMatrix, numAtoms);*/

	for(index=1;index<numAtoms;index++)
	{
		if (theGradTable[index].ab==TRUE)
		{
			theMatrix[index].ab-=theMatrix[index].ab*
			pow(10.0, -fabs(log10(fabs(theResults[gCount])))/100.0)
			*(theResults[gCount]/fabs(theResults[gCount]));
			gCount++;
		}
		if (theGradTable[index].abc==TRUE)
		{
			theMatrix[index].abc-=theMatrix[index].abc*
			pow(10.0, -fabs(log10(fabs(theResults[gCount])))/1000.0)
			*(theResults[gCount]/fabs(theResults[gCount]));
			gCount++;
		}
		if (theGradTable[index].abc_bcd==TRUE)
		{
			theMatrix[index].abc_bcd-=theMatrix[index].abc_bcd*
			pow(10.0, -fabs(log10(fabs(theResults[gCount])))/1000.0)
			*(theResults[gCount]/fabs(theResults[gCount]));
			gCount++;
		}
	}
	
/*see comment above*/
/*	printMatrix(theMatrix, numAtoms);*/
}	

double correctAngle(double theAngle)
{
	double temp;
	short done=FALSE;
	
	temp=theAngle;
	if(temp>pi)
		{
			done=FALSE;
			while (done==FALSE)
			{
				temp-=pi;
				if(temp<pi) done=TRUE;
			}
		}
	else if(temp<-pi)
		{
			done=FALSE;
			while (done==FALSE)
			{
				temp+=pi;
				if(temp>-pi) done=TRUE;
			}
		}
	return temp;
}
