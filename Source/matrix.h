/****************************************************

	matrix.h---the header for matrix routines
	
******************************************************/

void sqTranspose(short size, double **source, double **dest);
void sqAdd(short size, double **souce1, double **souce2, double **dest);
void sqMul(short size, double **source1, double **source2, double **dest);
void vecMatMult(short size, double *vec, double **mat, double *dest);
void matInv(short n, double *a, double **y, double *indx);
void LUDCMP(double *a, short n, short np, double *indx, short d);
void LUBKSB(double *a, short n, short np, double *indx, double *b);