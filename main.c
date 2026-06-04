#define _CRT_SECURE_NO_DEPRECATE			// fopen(), fscanf() warnings
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "matrix.h"
#include "dist.h"
#include "t_test.h"

static double mean(Matrix m, unsigned int column)
{
	double sum = 0;
	column %= m.cols;
	double *d = m.data + column++;

	for (int i = 0; i < m.rows; i++)
	{
		sum += *d;
		d += column;
	}
	return sum / m.rows;
}

static Metrics regress(Matrix x, Matrix y, double *coefficientMetrics)
{
  int i, j;
  double variance = 0;
  double yMean = mean(y, 0);
  Matrix xtrans = transMatrix(x);
  Matrix B = genCoefficients(x, y, xtrans); // B.rows = x.cols; B.cols = y.cols = 1
  Matrix Yhat = multiMatrix(x, B);	//  y estimates
  // variance = sum((y-Yhat)**2)
  Matrix residuals = calcResiduals(x, y, Yhat, &variance);
  Matrix stdErrMatrix = stdErr(x, variance, xtrans);
  double *m = stdErrMatrix.data;
  printf("stdErrmatrix\n");
  for (i = 0; i < stdErrMatrix.rows; i++)
  {
	printf("%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
			m[0], m[1], m[2], m[3], m[4], m[5], m[6]);
	m += 7;
  }

  // predicted values:  x.rows by B.cols = 1
  Metrics modelMetrics = { 0 };
  
  /*
   * Generate Model Metrics and store in order (for later printing)
   */

  // Regression Sum of Squares:  explained / estimated variance
  //Model Sum of Squares
  modelMetrics.SSR = 0;
  for(i = 0; i < y.rows; i++)
	modelMetrics.SSR += (Yhat.data[i] - yMean) * (Yhat.data[i] - yMean);

  //Model mean square; recall that one x column is NOT independent variables
  double MSR = modelMetrics.SSR / (x.cols - 1);	// Mean Square Regression

  //Residuals Sum of Squares (RSS):  Unexplained Variance
  modelMetrics.RSS = 0;
  for(i = 0; i < y.rows; i++)
	modelMetrics.RSS += (y.data[i] - Yhat.data[i]) * (y.data[i] - Yhat.data[i]);

  //Standard Error of Estimate https://www.socscistatistics.com/tests/errorofestimate/
  modelMetrics.SEE = (modelMetrics.RSS / (x.rows - 1));

  //Residuals mean square = modelMetrics.RSS / dof
  modelMetrics.RMS = modelMetrics.RSS / (x.rows - x.cols);

  // https://www.numberanalytics.com/blog/ultimate-f-test-regression-guide#calculating-the-f-test
  // F-test statistic = mean square Regression / Residuals mean square
  modelMetrics.F = MSR / modelMetrics.RMS;

  // F-test p-value
  modelMetrics.p_value = gsl_cdf_fdist_P(modelMetrics.F, x.cols - 1, x.rows - x.cols);

  //Total Sum of Squares (TSS)
  modelMetrics.TSS = 0;
  for(i = 0; i < y.rows; i++)
	modelMetrics.TSS += (y.data[i] - yMean) * (y.data[i] - yMean);

  //R squared
  modelMetrics.R2 = 1 - (modelMetrics.RSS / modelMetrics.TSS);

  //Adjusted R squared
  modelMetrics.AR2 = 1 - (1 - modelMetrics.R2) * (double)(x.rows - 1) / (x.rows - x.cols);

  // Model estimate
  modelMetrics.Mest = mean(Yhat, 0);

  // Model t-value:  estimate - actual / standard error
  modelMetrics.Mtv = (modelMetrics.Mest - yMean) / modelMetrics.SEE;

  //Root mean square error = sqrt(RSS / x.rows)
  modelMetrics.RMSE = sqrt(modelMetrics.RSS / x.rows);

  /*
  * Generate Coefficient Metrics and store in order (for later printing)
  */

  i = 0; //Track metric in the Coefficent Metric array

  //Non constant variables
  printf("stdErrmatrix for t-value\n");
  double d;
  for(j = 1; j < x.cols; j++) {	// independent variable
	// Model coefficient
	coefficientMetrics[i] = B.data[j];
	i++;
	//Standard error for regression coefficients 
	coefficientMetrics[i] = sqrt(d = stdErrMatrix.data[j * x.cols + j]);
    printf("[% d, % d]: % .4f ", j, j, d);
	i++;
	//t test statistic
	coefficientMetrics[i] = coefficientMetrics[i-2] / (coefficientMetrics[i-1]);
	i++;
	//t test p-value
	coefficientMetrics[i] = critical_value(x.rows - 1);
	i++;
	//Confidence inteval lower
	coefficientMetrics[i] = coefficientMetrics[i-4] - 1.96 * coefficientMetrics[i-3];
	i++;
	//Confidence inteval upper
	coefficientMetrics[i] = coefficientMetrics[i-5] + 1.96 * coefficientMetrics[i-4];
	i++;
  }
  printf("\n");

  //Constant
  j = 0;
  //Coefficient
  coefficientMetrics[i] = B.data[j];
  i++;
  //Standard error
  coefficientMetrics[i] = sqrt(stdErrMatrix.data[j * x.cols + j]);
  i++;
   //t test statistic
  coefficientMetrics[i] = coefficientMetrics[i-2] / (coefficientMetrics[i-1]);
  i++;
  //t test p-value
  coefficientMetrics[i] = critical_value(x.rows - 1);
  i++;
  //Confidence inteval lower
  coefficientMetrics[i] = coefficientMetrics[i-4] - 1.96 * coefficientMetrics[i-3];
  i++;
  //Confidence inteval upper
  coefficientMetrics[i] = coefficientMetrics[i-5] + 1.96 * coefficientMetrics[i-4];
  i++;

  free(B.data);
  free(xtrans.data);
  free(residuals.data);
  free(stdErrMatrix.data);
  free(Yhat.data);
  return modelMetrics;
}

static void loadXY(Matrix *data, Matrix *x, Matrix *y)
{
  size_t size = x->rows = y->rows = data->rows;
  x->cols = data->cols;
  y->cols = 1;
  size *= x->cols;
  x->data = (double *)calloc(size, sizeof(double));
  y->data = (double *)calloc(y->rows, sizeof(double));
  double foo, *data_pt = data->data, *x_pt = x->data, *y_pt = y->data;

  if (NULL != y_pt && NULL != x_pt)
  for (int j = 0; j < x->rows; j++)
  {
	*x_pt++ = foo = 1.0;		// insert constant coefficients
	*y_pt++ = foo = *data_pt++;
	for (double *z = x_pt + data->cols - 1; x_pt < z; x_pt++)
		*x_pt = foo = *data_pt++;
  }
}

static int readData(FILE *fp, char *varNames[10], int columns, Matrix *x, Matrix *y)
{
  int i, size = 100;
  double tempDouble;
  Matrix read = {0};

  read.data = malloc(sizeof(double) * size);
  if(read.data == NULL)
	return 0;

  for (i = 0; fscanf(fp,"%lf,", &tempDouble) != EOF; i++) {
	if(i >= size - 1) {
	  size += 100;
	  double *more = realloc(read.data, size * sizeof(double));
	  if (NULL == more)
	  {
		  free(read.data);
		  exit(1);
	  } else {
		read.data = more;
		read.data[i] = tempDouble;
	  }
	} else read.data[i] = tempDouble;
  }

  read.rows = i / columns;		// rows
  read.cols = columns;

  loadXY(&read, x, y);
  free(read.data);
  return 1;
}

static void vnprint(char *varName)
{
	printf("\n ");
	size_t l = strlen(varName);
	for(int j = 0; j < 9; j++)
	  printf("%c", j < l ? varName[j] : ' '); 
}

static void printModel(char *varNames[10], Matrix x, Matrix y)
{
  double *coefficientMetrics = (double *)malloc(sizeof(double) * 6 * x.cols);
  if (NULL == coefficientMetrics)
	  return;

  Metrics modelMetrics = regress(x, y, coefficientMetrics);

  /*
   * Print equation
   */
  printf("\nRegression Model Equation:  significant Coef. if |t-value| > critial\n%s = %.2lf",
	varNames[0], coefficientMetrics[x.cols * 6 - 6]);
  for(int i = 1, j = 0; i < x.cols; i++, j += 6)
	printf(" %+.3lf %s",coefficientMetrics[j], varNames[i]);
  printf(";  RMS error = %.3lf", modelMetrics.RMSE);

  /*
   * Print model metrics
   */
  printf("\n\n Source    |  Sum of       dof  Mean      0.05 significance, %d observations", x.rows);
  printf("\n           |  Squares           Squares                  F(%d, %d) = %.0lf",
    x.cols - 1, x.rows - x.cols, modelMetrics.F);
  printf("\n-----------+------------------------------               F p-value     =  %6.4lf", modelMetrics.p_value);
  printf("\n Model     |  %10.3lf %5d  %10.3lf       significant if 0.05 > F p-value",
   modelMetrics.SSR, x.cols - 1, modelMetrics.SSR/(x.cols - 1));
  printf("\n Residuals |  %10.3lf %5d  %10.3lf               R-squared     =  %6.3lf",
   modelMetrics.RSS, x.rows - x.cols, modelMetrics.RMS, modelMetrics.R2);
  printf("\n-----------+------------------------------               Adj R-squared =  %6.3lf", modelMetrics.AR2);
  printf("\n Total     |  %10.3lf %5d  %10.3lf               Root RSS      =  %6.3lf",
   modelMetrics.TSS, x.rows - 1, modelMetrics.TSS/(x.rows - 1), sqrt(modelMetrics.RSS));

  /*
   * Print coefficient metrics
   */
  char *s, **vn = varNames;
  printf(s = "\n-----------+--------------------------------------------------------------------");
  vnprint(*vn++);
  printf(" |      Coef.   Std. Err.   t-value critical      [95%% Conf. Interval]%s", s);
  for (int k = 0; k < x.cols; k++)
  {
	vnprint(k == x.cols - 1 ? "Intercept" : *vn++);
	int j = k * 6;
	printf(" |%11.3lf %11.3lf %9.2lf  %4.3lf  %13.3lf  %10.3lf",
	  coefficientMetrics[j],  coefficientMetrics[j+1], coefficientMetrics[j+2],
	  coefficientMetrics[j+3],  coefficientMetrics[j+4],  coefficientMetrics[j+5]);
  }
  printf(s);
}

int main(int argc, char **argv)
{
//char *fin = (1 == argc) ? "../../../data/health_data.txt" : argv[1];
  char *fin = (1 == argc) ? "../../../data/Before_redefineG.txt" : argv[1];

  FILE *text = fopen(fin, "r");
  if(text == NULL)
	return - printf("Unable to open data file '%s'.", fin);

  char *varNames[10] = {""};
  char varString[100];
  int j, columns, gp = 0;
  int i = fscanf(text, "%[^\n]", varString);
  varString[99] = '\0';
  size_t l = strlen(varString);

  if (7 < l && 0 == strncmp("# rows ", varString, 7)) // intercept gnuplot data
  {
	char *endptr;
	gp = strtol(7 + varString, &endptr, 10);
	fgets(varString, 98, text);				// pop that first line
	// next line should be parm names:
	// xG,yG,dxR,dyR,dxB,dyB,xG2,yG2,xG3,yG3
	i = fscanf(text, "%[^\n]", varString);
	varString[99] = '\0';
	l = strlen(varString);
  }

  for(columns = i = j = 0; i <= l; i++)
	if(varString[i] == ',' || varString[i] == '\0')
	{
		varString[i] = '\0';	
  		varNames[columns++] = j + varString;
		j = 1 + i;		// next varNames
	}


  Matrix x = {0}, y = {0};
  if (0 < gp)
	gpdata(text, varNames, columns, gp, &x, &y);
  else readData (text, varNames, columns, &x, &y);

  if (NULL == x.data)
  {
	printf("******************************************************************************\n");
	printf("*                                                                            *\n");
	printf("* Title: Multiple Linear Regression                                          *\n");
	printf("* Description:  This program takes an inputted data file and performs        *\n");
	printf("*               multiple linear regression analysis on the data.             *\n");
	printf("* Author:       Oscar Zealley                                                *\n");
	printf("* Instructions: Put your data in a .txt file in the same directory as        *\n");
	printf("*               this program. Data must be formatted like this:              *\n");
	printf("*                                                                            *\n");
	printf("*  Dependent Variable, Independent Variable 1, Indendent Variable 2,...      *\n");
	printf("*                   4,                      5,                    8,...      *\n");
	printf("*                   7,                     12,                    5,...      *\n");
	printf("*                   .                       .                     .          *\n");
	printf("*                   .                       .                     .          *\n");
	printf("*                                                                            *\n");
	printf("* NB: Max number of variables is 10.                                         *\n");
	printf("*                                                                            *\n");
	printf("* Press enter for sample data statistics.                                    *\n");
	printf("*                                                                            *\n");
	printf("******************************************************************************\n");

	char response = getchar();
	  return -1;
  }

  printModel(varNames, x, y);
  
  free(x.data);
  free(y.data);
  return 0;
}
