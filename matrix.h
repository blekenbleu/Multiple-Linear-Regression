/*
*Functions for manipulating mathematical matrices
*/

#ifndef _MATRIX_
#define _MATRIX_

typedef struct
{
  int rows;
  int cols;
  double *data;
} Matrix;

//My matrix functions
extern Matrix initMatrix(int rows, int cols);
extern Matrix loadMatrix(int rows, int cols, double *data);
extern void printMatrix(Matrix m);
extern Matrix addMatrix(Matrix a, Matrix b);
extern Matrix multiMatrix(Matrix a, Matrix b);
extern Matrix transMatrix(Matrix m);	// used only to generate xtrans from x
extern int squareMatrix(Matrix m, double square[25][25]);
extern Matrix inverseMatrix(Matrix m);
extern Matrix genCoefficients(Matrix x, Matrix y, Matrix xtrans);
extern Matrix calcResiduals(Matrix x, Matrix y, Matrix Yhat, double * errStdDev);
extern Matrix stdErr(Matrix x, double errStdDev, Matrix xtrans);

//Inverse matrix functions from: http://scanftree.com/programs/c/c-program-to-find-the-inverse-of-the-matrix/
extern double determinant(double[25][25], double);
extern void cofactors(double[25][25], double);
extern void trans(double[25][25], double[25][25], double);

#endif
