/***************************************************************************

potentials.c

***************************************************************************/

#include "potentials.h"
#include "zmat.h"

#define del(x)  x*1e-6


double electro(vector r1, double q1, vector r2, double q2)
{
	double result, mag, temp2;
	vector tempv;
	
	v_copy(r1, tempv);
	v_sub(tempv, r2);
	mag=v_mag(tempv);
	mag*=1e-10;
	
	result=q1*q2*e2/mag;
	
	/*returns answer in coulombs, div by 4*pi*epsilon0 to get joules*/
	result /= FourPiEpsilon0;

/* commented out.  only used in debugging
	printf ("%3.2e\t%3.2e\n",result, mag);*/
	return(result);
}

double leonard_jones (vector r1, vector r2)
{
	double 	result, mag;
	vector 	tempv;
		
	v_copy(r1, tempv);
	v_sub(tempv, r2);
	mag=v_mag(tempv);
	
	result = (A/pow(mag,12.0))-(B/pow(mag,6.0)); 
	/*returns kJ/mol div by N0 then multiply by 1000 to get joules*/
	
	result /= N0;
	result *=1000;

/* commented out.  only used in debugging
	printf ("\t\t%3.2e\t%3.2e\n",result, mag);*/
	return(result);
}

double energy_tip(double *theQs, vector *theCarts, short numAtoms)

/*I screwed up.  This has been hacked to be TIP4P dimer specific.  IE
	atoms on water 1 interract only with those on atom 2 and vice versa
	maybe a future version will be more general and thus applicable to
	a wider variety of applications*/

{
	short 	i,j, O1=1, O2=5;
	double 	result=0.0, qi, qj;
	
	/*numbers are obscure.  #define atom names and put them in instead*/
	
	result+=electro(theCarts[0], theQs[0], theCarts[4], theQs[4]);
	result+=electro(theCarts[0], theQs[0], theCarts[6], theQs[6]);
	result+=electro(theCarts[0], theQs[0], theCarts[7], theQs[7]);
	result+=electro(theCarts[2], theQs[2], theCarts[4], theQs[4]);
	result+=electro(theCarts[2], theQs[2], theCarts[6], theQs[6]);
	result+=electro(theCarts[2], theQs[2], theCarts[7], theQs[7]);
	result+=electro(theCarts[3], theQs[3], theCarts[4], theQs[4]);
	result+=electro(theCarts[3], theQs[3], theCarts[6], theQs[6]);
	result+=electro(theCarts[3], theQs[3], theCarts[7], theQs[7]);

	result+=leonard_jones(theCarts[1], theCarts[5]);
	
	return result;
}


double gradient_tip(zm_el *theMatrix, del_el *theDeltas, vector *theCarts,
		 double *theQs, bonds *theBonds, short numAtoms, double *theResults)
{
	double 	temp1, temp2, delta;
	short 	index, index2=0;
	vector 	*newCarts;
	
	newCarts=(vector *)calloc((size_t)numAtoms, sizeof(vector));

	temp1=energy_tip(theQs, theCarts, numAtoms);

	for (index=1;index<numAtoms;index++)
	{
		if(theDeltas[index-1].ab==TRUE)
		{
			delta=del(theMatrix[index].ab);
			theMatrix[index].ab+=delta;
			makeCartesian(theMatrix, theBonds, newCarts, numAtoms);	
			temp2=(energy_tip(theQs, newCarts, numAtoms));
			theMatrix[index].ab-=delta;
			theResults[index2++]=(temp2-temp1)/delta;
/*			printf("%4.3e\t%4.3e\t%4.3e\n", temp1, temp2, delta);
		*/}
		if(theDeltas[index-1].abc==TRUE)
		{
			delta=del(theMatrix[index].abc);
			theMatrix[index].abc+=delta;
			makeCartesian(theMatrix, theBonds, newCarts, numAtoms);	
			temp2=(energy_tip(theQs, newCarts, numAtoms));
			theMatrix[index].abc-=delta;
			theResults[index2++]=(temp2-temp1)/delta;
		/*	printf("%4.3e\t%4.3e\t%4.3e\n", temp1, temp2, delta);
		*/}
		if(theDeltas[index-1].abc_bcd==TRUE)
		{
			delta=del(theMatrix[index].abc_bcd);
			theMatrix[index].abc_bcd+=delta;
			makeCartesian(theMatrix, theBonds, newCarts, numAtoms);	
			temp2=(energy_tip(theQs, newCarts, numAtoms));
			theMatrix[index].abc_bcd-=delta;
			theResults[index2++]=(temp2-temp1)/delta;
		/*	printf("%4.3e\t%4.3e\t%4.3e\n", temp1, temp2, delta);
		*/}
	}
	free(newCarts);
	return temp1;
}

short makeGradTable(del_el *theGradTable, FILE *theFile)
{
	
	struct 
	{
		short ab, abc, abc_bcd;
	} temp_el;
	
	short 	index=0, count=0, done=FALSE;
	char 	*theLine;
	
	theFile=fopen(theGradName, "r");

	theLine=(char *)malloc(LINELEN+1);

	while(done==FALSE)
	{
		theLine=fgets(theLine, LINELEN, theFile);
		if (theLine[0]==ENDOFGRAMARK)
			done=TRUE;
		else
		{
			sscanf(theLine, "%hd%hd%hd",
					&temp_el.ab,&temp_el.abc,&temp_el.abc_bcd);
			if (temp_el.ab==TRUE)
			{
				theGradTable[index].ab=TRUE;
				count++;
			}
			else theGradTable[index].abc=FALSE;
			if (temp_el.abc==TRUE)
			{
				theGradTable[index].abc=TRUE;
				count++;
			}
			else theGradTable[index].abc=FALSE;
			if (temp_el.abc_bcd==TRUE)
			{
				theGradTable[index].abc_bcd=TRUE;
				count++;
			}
			else theGradTable[index].abc_bcd=FALSE;
			index++;
		}
	}
	return count;
}

