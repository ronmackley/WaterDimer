/*************************************************************************

	newton.h
	
	Header file for newtons methods stuff

*************************************************************************/

#include <stdio.h>

//#include "types.h"
#include "stdDefines.h"
#include "potentials.h"
#include "vector.h"
#include "zmat.h"
//#include "potentials.h"
				
double Newtons(zm_el *, del_el *, vector *, double *, bonds *,
		 short , short , double *);
void update(double *, zm_el *, del_el *, double, short , short );
double correctAngle(double);

