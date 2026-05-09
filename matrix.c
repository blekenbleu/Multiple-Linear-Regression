#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "matrix.h"

Matrix initMatrix(int rows, int cols)
{
  Matrix result = { 0 };
  result.rows = rows;
  result.cols = cols;
  size_t size = result.rows;
  size *= result.cols;
  result.data = (double *)calloc(size, sizeof(double));

  return result;
}

Matrix loadMatrix(int rows, int cols, double * elements)
{
  Matrix result = { 0 };
  result.rows = rows;
  result.cols = cols;
  size_t size = result.rows;
  size *= result.cols;
  result.data = (double *)calloc(size, sizeof(double));

  for(int i = 0; i < rows * cols; i++) {
    result.data[i] = elements[i];
  }

  return result;
}

void printMatrix(Matrix m)
{
  printf("\nRows: %d\tColumns: %d\n", (int)m.rows, (int)m.cols);
  for(int i = 0; i < m.rows; i++) {
    printf("| ");
    for(int j = 0; j < m.cols; j++) {
      printf("%4.2lf ", m.data[i*m.cols+j]);
    }
    printf(" |\n");
  }
}

Matrix addMatrix(Matrix a, Matrix b)
{
  if(a.rows != b.rows || a.cols != b.cols) {
    printf("Can't add matrices of different dimensions");
    return a;
  }

  else {
    Matrix result = initMatrix(a.rows, a.cols);

    for(int i = 0; i < a.rows * a.cols; i++) {
      result.data[i] = a.data[i] + b.data[i];
    }

    return result;
  }

}

Matrix multiMatrix(Matrix a, Matrix b)
{
  if(a.cols != b.rows) {
    printf("Can't multiply matrices");
    return a;
  }

  else {
    Matrix result = initMatrix(a.rows, b.cols);
    
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            double sum = 0;
            for (int k = 0; k < a.cols; k++)
                sum = sum + a.data[i * a.cols + k] * b.data[k * b.cols + j];
            result.data[i * b.cols + j] = sum;
        }
    }

    return result;
  }
}

Matrix transMatrix(Matrix m)
{
  Matrix result = initMatrix(m.cols, m.rows);
  int k = 0;
  for(int j = 0; j < m.rows; j++) {
    for(int i = 0; i < m.cols; i++) {
      result.data[i*m.rows + j] = m.data[k];
      k++;
    }
  }
  return result;
}

int squareMatrix(Matrix m, double square[25][25])
{
  if(m.rows != m.cols || m.rows > 24) {
    printf("Matrix isn't square!");
    return 0;
  }

  else {

    int k = 0;
    for(int i = 0; i < m.rows; i++) {
      for(int j = 0; j < m.rows; j++) {
        square[i][j] = m.data[k];
        k++;
      }
    }
    return m.rows;
  }
}

Matrix inverseMatrix(Matrix m)
{
  double squareTemp [25][25];
  memset(squareTemp, 0, 625 * sizeof(double)); // 25 * 25
  Matrix result = initMatrix(m.rows,m.rows);
  int n = squareMatrix(m, squareTemp);
  double d = determinant(squareTemp, n);
 
	//printf("\nThe determinant is: %.0f", d);
 
	if (d == 0) {
	  printf("\nMATRIX IS NOT INVERSIBLE\n"); 
  }
  else {
	  cofactors(squareTemp, n);
  }
  
  int k = 0;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      result.data[k] = squareTemp[i][j];
      k++;
    }
  }

  return result;

}

Matrix genCoefficients(Matrix x, Matrix y)
{
  Matrix B = initMatrix(x.cols, x.cols);

  B = multiMatrix((inverseMatrix(multiMatrix(transMatrix(x), x))),(multiMatrix(transMatrix(x),y)));

  return B;
}

Matrix calcResiduals(Matrix x, Matrix y, Matrix B, double * errStdDev)
{
  int i;
  double sum = 0;
  Matrix result = initMatrix(y.rows, 1);
  Matrix Yhat = initMatrix(y.rows, 1);
  Yhat = multiMatrix(x, B);

  for(i = 0; i < y.rows; i++) {
    result.data[i] = y.data[i] - Yhat.data[i];
  }

  for(i = 0; i < result.rows; i++) {
    sum += result.data[i] * result.data[i];
  }

  *errStdDev = sqrt(sum / (result.rows - (x.cols - 1) - 1));

  free(Yhat.data);

  return result;
}

Matrix stdErr(Matrix x, double errStdDev)
{
  int i;
  Matrix result = initMatrix(x.cols, x.cols);

  result = inverseMatrix(multiMatrix(transMatrix(x),x));

  for(i = 0; i < x.cols * x.cols; i++) {
    result.data[i] *= errStdDev * errStdDev;
  }

  return result;
}

//Inverse matrix functions from : http://scanftree.com/programs/c/c-program-to-find-the-inverse-of-the-matrix/

double determinant(double a[25][25], double k) {
	double s = 1, det = 0, b[25][25] = { 0 };
	int i, j, m, n, c;
	if (k == 1) {
		return (a[0][0]);
	} else {
		det = 0;
		for (c = 0; c < k; c++) {
			m = 0;
			n = 0;
			for (i = 0; i < k; i++) {
				for (j = 0; j < k; j++) {
					b[i][j] = 0;
					if (i != 0 && j != c) {
						b[m][n] = a[i][j];
						if (n < (k - 2))
						       n++; else {
							n = 0;
							m++;
						}
					}
				}
			}
			det = det + s * (a[0][c] * determinant(b, k - 1));
			s = -1 * s;
		}
	}
	return (det);
}
 
void cofactors(double num[25][25], double f) {
	double b[25][25] = { 0 }, fac[25][25] = { 0 };
	int p, q, m, n, i, j;
	for (q = 0; q < f; q++) {
		for (p = 0; p < f; p++) {
			m = 0;
			n = 0;
			for (i = 0; i < f; i++) {
				for (j = 0; j < f; j++) {
					b[i][j] = 0;
					if (i != q && j != p) {
						b[m][n] = num[i][j];
						if (n < (f - 2))
						       n++; else {
							n = 0;
							m++;
						}
					}
				}
			}
			fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
		}
	}
	trans(num, fac, f);
}
 
void trans(double num[25][25], double fac[25][25], double r) {
	int i, j = 0;
	double b[25][25] = { 0 }, inv[25][25] = { 0 }, d;
	for (i = 0; i < r; i++) {
		for (j = 0; j < r; j++) {
			b[i][j] = fac[j][i];
		}
	}
	d = determinant(num, r);
	inv[i][j] = 0;
	for (i = 0; i < r; i++) {
		for (j = 0; j < r; j++) {
			inv[i][j] = b[i][j] / d;
		}
	}

  for (i = 0; i < r; i++) {
		for (j = 0; j < r; j++) {
			num[i][j] = inv[i][j];
		}
	}
}