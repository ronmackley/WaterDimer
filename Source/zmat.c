/****************************************************************

	zmat.c---this is the code to read a zmatrix from a file
	and build a list of cartesian coordinates from it.
	Also ,I will, from the cartesian coords, build a 
	zmatrix, if I ever figure out how.
		
*****************************************************************/

#include "zmat.h"
#include "dimer.h"

void notComment(FILE *fp, char *theLine)

/*
	reads and reads and reads until it finds a line which
	does not begin with COMMENTMARK, then it returns that
	line in theLine
*/

{
	int 	flag=FALSE;
	char 	*temp;
	
	while(flag==FALSE)
	{
		temp=fgets(theLine, LINELEN, fp);	
		if(theLine[0]!=COMMENTMARK)
			flag=TRUE;
	}
}


short countAtoms(FILE *fp)
{
	char 	*theLine;
	short 	numAtoms=0, done=FALSE;
	
	theLine=(char *)malloc(LINELEN+1);
	fp=fopen(theName,"r");
	while((done==FALSE)&&!feof(fp))
	{
		notComment(fp, theLine);
		/* terminate the z-matrix with a @ aka ENDOFMATRIX */
		if(theLine[0]!=ENDOFMATRIX)
			numAtoms++;
		else done=TRUE;
	}
	free(theLine);
	return numAtoms;
}

short readZmat(FILE *fp, zm_el *theMatrix, bonds *theBonds, 
						vector *theCarts, short numAtoms)

{
	char 	*theLine;
	short 	i;
		
	theLine=(char *)malloc(LINELEN+1);

	fp=fopen(theName,"r");

	for (i=0;i<numAtoms;i++)
	{
		notComment(fp, theLine);
		sscanf(theLine, "%hd%hd%Lf%hd%Lf%hd%Lf",
				&theMatrix[i].a,
				&theMatrix[i].b,
				&theMatrix[i].ab, 
				&theMatrix[i].c,
				&theMatrix[i].abc, 
				&theMatrix[i].d,
				&theMatrix[i].abc_bcd);
		
		/*theMatrix[i].ab*=1e-10; convert form A to m*/
		theMatrix[i].abc*=deg2rad;
		theMatrix[i].abc_bcd*=deg2rad;
		
		addBond(theBonds, i, theMatrix[i].b);
	}
	while(!feof(fp))
	{
		short 	a, b;
		
		notComment(fp, theLine);
		sscanf(theLine, "%hd%hd", &a, &b);
		addBond(theBonds, a, b);
	
/*
	alas, this should take care of any additional bonds being specified
  	the bond count is critical, tho, as any lines not comments after
  	numAtoms lines are treated to be additional linkages.
*/
	
	}
	fclose(fp);
	
	free(theLine);
	return i;
}

void makeCartesian(zm_el *theMatrix, bonds *theBonds,
					 vector *theCarts, short numAtoms)

/*
	this procedure is the crux of the program. Without it, this progrm
	would be useless.  It was tested on specially designed tetrahedral
	and octahedral and some asymetric sample z-matricies and worked on 
	those systems.  It should therefore work with any system, as there 
	was nothing special about the test systems, however, THIS HAS NOT
	BEEN VERIFIED and may be why this program doesn't quite work properly
*/

{
	short 			index;
	const vector	theOrigin={0.0,0.0,0.0};
	
/* set up the first atom at the origin */

	v_copy(theOrigin, theCarts[0]);

/* the second lies along the abritrarily chosen z-axis @ ab*/

	theCarts[1][Vx]=0.0;
	theCarts[1][Vy]=0.0;
	theCarts[1][Vz]=theMatrix[1].ab;
	
/*
 	the third bond lies at (ab, abc) in polar coordinates, assuming the origin
   	at b.  Convert to cart so that it lies in the arbitrarly chosen xz plane.
*/

	theCarts[2][Vr]=theMatrix[2].ab;
	theCarts[2][Vtheta]=0.0;
	theCarts[2][Vphi]=theMatrix[2].abc;

	v_sph2cart(theCarts[2]);
	v_add(theCarts[2], theCarts[theMatrix[2].b]);
	
/* 	everything from here on in copes with zmatrix rows with dihedral angles */

	for (index=3;index<numAtoms;index++)
	{
		vector axis, cd, tc, td;
		double cdTheta;
		
/*	get a vector for the axis for future reference */
		v_copy(theCarts[theMatrix[index].c], axis);
		v_sub(axis, theCarts[theMatrix[index].b]);

/*
	and put it in spherical coords this will give me theta and phi
	describing the orientation of the bc axis relative to the "absolute
	coords of the rest of the atom 
*/

		v_cart2sph(axis);

/*	get the cd vector so that i can get the dehedral wrt. the xz plane */

		v_copy(theCarts[theMatrix[index].d], td);
		v_copy(theCarts[theMatrix[index].c], tc);
		
		v_cart2sph(tc);
		v_cart2sph(td);
		
		tc[Vtheta]-=axis[Vtheta];
		td[Vphi]-=axis[Vphi];

		v_sph2cart(tc);
		v_sph2cart(td);
		
		v_copy(td, cd);
		v_sub(cd, tc);

/*	find the orientation wrt the xz plane via an impromptu dot product */ 

		cdTheta=acos(cd[Vx]/(v_mag(cd)));

/* set up the atom for positionin g in spherical coords */
		
		theCarts[index][Vr]=theMatrix[index].ab;
		theCarts[index][Vtheta]=theMatrix[index].abc_bcd;
		theCarts[index][Vphi]=theMatrix[index].abc;		

/*
	I believe this should give me cartesian coordinates for the atom assuming
	the origin lies at atom b and an axis (z I think) lies along the bc bond

	now it must be translated.  axis contains the orientation of the bond axis
	from where we derive the position of the new atom.  The vector cb is used
	as the Z-axis.  So, add the angle of the axis and we will now get 
	orientation of a wrt. the absolute coords. 
*/

		theCarts[index][Vtheta]-=cdTheta;

		theCarts[index][Vtheta]-=axis[Vtheta];
		theCarts[index][Vphi]-=axis[Vphi];
		
/*
	now we put it into cartesian coords and add the position of b to it and
	we have it.  the position of a derived from its z-matrix entry.
*/

		v_sph2cart(theCarts[index]);
		v_add(theCarts[index], theCarts[theMatrix[index].b]);
		
		
		
/* and it's done!!!! */
	}

}

short makeZmat(zm_el *theMatrix, bonds *theBonds,
					 vector *theCarts, short numAtoms)
/*
	as of now this function will take the order of the current zmatrix
	and compute all of the new bond lengths and angles based on the
	current configuration.  To make it of general utility, I must ba able
	to compute a Zmatrix given any first atom.  Maybe someday I will do it.
	This has never been tested. 
*/

{
	short	index;
	vector	Scratch;


	v_copy(theCarts[1], Scratch);
	v_sub(Scratch, theCarts[theMatrix[1].b]);
	theMatrix[1].ab=v_mag(Scratch);

	v_copy(theCarts[2], Scratch);
	v_sub(Scratch, theCarts[theMatrix[2].b]);
	theMatrix[2].ab=v_mag(Scratch);
	theMatrix[2].abc=bondAng(theCarts, 2, theMatrix[2].b, theMatrix[2].c);
	
	for(index=3;index<numAtoms;index++)
	{	
		v_copy(theCarts[index], Scratch);
 		v_sub(Scratch, theCarts[theMatrix[index].b]);
		theMatrix[index].ab=v_mag(Scratch);
		theMatrix[index].abc=bondAng(theCarts, index,
										theMatrix[index].b,
										theMatrix[index].c);
		theMatrix[index].abc_bcd=dihedralAng(theCarts, index,
										theMatrix[index].b,
										theMatrix[index].c,
										theMatrix[index].d);
	}
	return 0;		
}


void makeQtable(zm_el *theMatrix, double *theQs, short numAtoms)
{
	short	index;
	
	for (index=0;index<numAtoms;index++)
	{
		switch (theMatrix[index].a)
		{
			case Hid:
			{
				theQs[index]=qH;
				break;
			}
			case Mid:
			{
				theQs[index]=qM;
				break;
			}
			default:
			{
				theQs[index]=qDefault;
				break;
			}
		}
	}
}


double bondAng(vector *theCarts, short a, short b, short c)
{
	vector 	ba, bc;
	double 	result;
	
	v_copy(theCarts[a], ba);
	v_copy(theCarts[c], bc);
	
	v_sub(ba, theCarts[b]);
	v_sub(bc, theCarts[b]);
	result=acos(v_dot(ba, bc)/(v_mag(ba)*v_mag(bc)));
	return result;
}

double dihedralAng(vector *theCarts, short a, short b, short c, short d)
{
	vector	ab, bc, cb, dc, abc, bcd;
	double	result, q1, q2;
	
/* cross ab into bc, and cd into bc and get the angle between the two 
   resultant vectors to get the dihedral angle  */
   
	v_copy(theCarts[a], ab);
	v_sub(ab, theCarts[b]);

	v_copy(theCarts[d], dc);
	v_sub(dc, theCarts[c]);

	v_copy(theCarts[b], bc);
	v_sub(bc, theCarts[c]);

	v_copy(bc, cb);
	v_inver(cb);
	
	v_cross(ab, bc, abc);
	v_cross(dc, cb, bcd);
	
	q1=v_mag(abc);
	q2=v_mag(bcd);
	
	result=acos(v_dot(abc, bcd)/(q1*q2));
	
	return result;
}
	

void addBond(bonds *theBonds, short a, short b)
{

/*
	this procedure is useless right now as I have no need for
	a bond list.  When mzkeZmat works properly, one will probably
	be needed.  this procedure is known to work, tho.
*/

	theBonds[a].toWhere[theBonds[a].theCount]=b;
	theBonds[b].toWhere[theBonds[b].theCount]=a;
	/* tell each atom what was bonded to it */
	theBonds[a].theCount++;
	theBonds[b].theCount++;
	/*increment the bond count for each atom */
}

short isBonded(bonds theBondRec, short serialNo)
{
	short i, flag=FALSE;
	
	for(i=0; i<theBondRec.theCount;i++)
	{
		if(theBondRec.toWhere[i]==serialNo)
		{
			flag=TRUE;
			break;
		}
	}
	return flag;
}