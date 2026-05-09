#define _CRT_SECURE_NO_DEPRECATE        // fopen(), fscanf() warnings
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "matrix.h"
#include "dist.h"

double * readData(FILE * fp, int * numVarPtr, int * sampleSizePtr, char variableNames[10][25]);
Matrix loadX(int numVar, int sampleSize, double * data);
Matrix loadY(int numVar, int sampleSize, double * data);
void printModel(char varNames[10][25], double modelMetrics[17], double * coefficientMetrics, int numVar);

static void regress(Matrix x, Matrix y, int numVar, int sampleSize, double modelMetrics[15], double * coefficientMetrics)
{
  int i, j;
  double errStdDev = 0;
  double * errStdDevPtr = &errStdDev;
  double sum = 0;
  double yMean = 0;
  for(i = 0; i < y.rows; i++) {
    sum+= y.data[i];
  }
  yMean = sum / y.rows;
  sum = 0;

  Matrix B = genCoefficients(x, y);
  Matrix residuals = calcResiduals(x, y, B, errStdDevPtr);
  Matrix stdErrMatrix = stdErr(x, errStdDev);
  Matrix Yhat = initMatrix(y.rows, 1);
  Yhat = multiMatrix(x, B);
  
  /*
  * Generate Model Metrics and store in order (for later printing)
  */

  // Number of observations
  modelMetrics[0] = sampleSize;

  // Degrees of freedom lower
  modelMetrics[1] = numVar - 1;

  //Degrees of freedom upper
  modelMetrics[2] = sampleSize - numVar;

  //Model Sum of Squares (ESS)
  for(i = 0; i < y.rows; i++) {
    sum += (Yhat.data[i] - yMean) * (Yhat.data[i] - yMean);
  }
  modelMetrics[4] = sum;
  sum = 0;

  //Degrees of freedom lower
  modelMetrics[5] = numVar - 1;

  //Model mean square
  modelMetrics[6] = modelMetrics[4] / modelMetrics[5];

  //Residuals Sum of Squares (SSR)
  for(i = 0; i < y.rows; i++) {
    sum += (y.data[i] - Yhat.data[i]) * (y.data[i] - Yhat.data[i]);
  }
  modelMetrics[8] = sum;
  sum = 0;

  //Degrees of freedom upper
  modelMetrics[9] = sampleSize - numVar;

  //Residuals mean square
  modelMetrics[10] = modelMetrics[8] / modelMetrics[9];

  //Total Sum of Squares (TSS)
  for(i = 0; i < y.rows; i++) {
    sum += (y.data[i] - yMean) * (y.data[i] - yMean);
  }
  modelMetrics[13] = sum;
  sum = 0;

  //Corrected degrees of freedom
  modelMetrics[14] = sampleSize - 1;

  //Total Mean Square
  modelMetrics[15] = modelMetrics[13] / modelMetrics[14];

  // F-test statistic
  modelMetrics[3] = modelMetrics[6] / modelMetrics[10];

  // F-test p value
  modelMetrics[7] = gsl_cdf_fdist_P(modelMetrics[3], modelMetrics[1], modelMetrics[2]);

  //R squared
  modelMetrics[11] = 1 - (modelMetrics[8] / modelMetrics[13]);

  //Adjusted R squared
  modelMetrics[12] = 1 - (1 - modelMetrics[11]) * ((modelMetrics[0] - 1) / (modelMetrics[0] - modelMetrics[5] - 1));

  //Root MSE
  modelMetrics[16] = modelMetrics[8] / modelMetrics[0] + 1;

  /*
  * Generate Coefficient Metrics and store in order (for later printing)
  */

  i = 0; //Track metric in the Coefficent Metric array
  j = 1; //Track independent variable 

  //Non constant variables
  for(; j < numVar; j++) {
    //Coefficient
    coefficientMetrics[i] = B.data[j];
    i++;
    //Standard error
    coefficientMetrics[i] = sqrt(stdErrMatrix.data[j * numVar + j]);
    i++;
    //t test statistic
    coefficientMetrics[i] = coefficientMetrics[i-2] / (coefficientMetrics[i-1]);
    i++;
    //t test p value
    coefficientMetrics[i] = 0; //PLACEHOLDER - FIND FUNCTION TO CALCUALTE T DIST
    i++;
    //Confidence inteval lower
    coefficientMetrics[i] = coefficientMetrics[i-4] - 1.96 * coefficientMetrics[i-3];
    i++;
    //Confidence inteval upper
    coefficientMetrics[i] = coefficientMetrics[i-5] + 1.96 * coefficientMetrics[i-4];
    i++;
  }

  //Constant
  j = 0;
  //Coefficient
  coefficientMetrics[i] = B.data[j];
  i++;
  //Standard error
  coefficientMetrics[i] = sqrt(stdErrMatrix.data[j * numVar + j]);
  i++;
   //t test statistic
  coefficientMetrics[i] = coefficientMetrics[i-2] / (coefficientMetrics[i-1]);
  i++;
  //t test p value
  coefficientMetrics[i] = 0; //PLACEHOLDER - FIND FUNCTION TO CALCUALTE T DIST
  i++;
  //Confidence inteval lower
  coefficientMetrics[i] = coefficientMetrics[i-4] - 1.96 * coefficientMetrics[i-3];
  i++;
  //Confidence inteval upper
  coefficientMetrics[i] = coefficientMetrics[i-5] + 1.96 * coefficientMetrics[i-4];
  i++;

  free(B.data);
  free(residuals.data);
  free(stdErrMatrix.data);
  free(Yhat.data);
}

int main(int argc, char **argv)
{
  char response = '0', *fin;
  char varNames[10][25];
  double modelMetrics[17];
  double * coefficientMetrics;
  int numVar = 0;
  int * numVarPtr = &numVar;
  int sampleSize = 0;
  int * sampleSizePtr = &sampleSize;

//FILE* text = fopen(fin = (1 == argc) ? "../../../data/health_data.txt" : argv[1], "r");
  FILE* text = fopen(fin = (1 == argc) ? "../../../data/Before_redefine.gp" : argv[1], "r");
  if(text == NULL) {
    printf("Unable to open data file '%s'.", fin);
    return 1;
  }

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
  
  response = getchar();

  double * data = readData(text, numVarPtr, sampleSizePtr, varNames);
  coefficientMetrics = (double*)malloc(sizeof(double) * 6 * (numVar));
  Matrix x = loadX(numVar, sampleSize, data);
  Matrix y = loadY(numVar, sampleSize, data);
  free(data);

  regress(x, y, numVar, sampleSize, modelMetrics, coefficientMetrics);
 
  printModel(varNames, modelMetrics, coefficientMetrics, numVar);
  
  free(x.data);
  free(y.data);
  return 0;
}

double * readData(FILE * fp, int * numVarPtr, int * sampleSizePtr, char variableNames[10][25])
{
  char varString[100];
  int i = fscanf(fp, "%[^\n]", varString);
  varString[99] = '\0';

  int j, ctr;
  size_t l = strlen(varString);
  j = 0; 
  ctr = 0;
  for(i = 0; i <= l; i++){
    if(varString[i] == ','||varString[i] == '\0') {
        variableNames[ctr][j]='\0';
        ctr++;
        j = 0;
    }
    else {
      variableNames[ctr][j]=varString[i];
      j++;
    }
  }
  int size = 100;
  double tempDouble;
  double * data = malloc(sizeof(double) * size);
  if(data == NULL) exit(1);
  i = 0;
  while(fscanf(fp,"%lf,", &tempDouble) != EOF) {
    if(i >= size - 1) {
      double* more;
      size += 100;
      if (NULL == (more = realloc(data, size * sizeof(double))))
      {
          free(data);
          exit(1);
      }
      else data = more;
    }
    if(NULL != data && i == 0 || (i % (ctr + 1)) == 0 ) {
      data[i] = 1.0;
      i++;
    }
    if (NULL != data)
      data[i] = tempDouble;
    i++;
  }
  *sampleSizePtr = i/(ctr + 1);
  *numVarPtr = ctr;

  return data;
}

Matrix loadX(int numVar, int sampleSize, double * data)
{
    Matrix result = { 0 };
  result.rows = sampleSize;
  result.cols = numVar;
  size_t size = result.rows;
  size *= result.cols;
  result.data = (double *)calloc(size, sizeof(double));
  int i,j;
  for(i = 0, j = 0; i < (numVar + 1) * sampleSize; i++) {
    if((i - 1) % (numVar + 1) != 0) {
      result.data[j] = data[i];
      j++;
    }
  }

  return result;
}

Matrix loadY(int numVar, int sampleSize, double * data)
{
    Matrix result = { 0 };
  result.rows = sampleSize;
  result.cols = 1;
  result.data = (double *)calloc(result.rows, sizeof(double));
  int i,j;
  for(i = 0, j = 0; i < (numVar + 1) * sampleSize; i++) {
    if((i - 1) % (numVar + 1) == 0) {
      result.data[j] = data[i];
      j++;
    }
  }

  return result;
}

void printModel(char varNames[10][25], double modelMetrics[17], double * coefficientMetrics, int numVar)
{
  int i = 0, j = 0, k;

  /*
  *Print equation
  */
  printf("\nRegression Equation:\n%s = %.2lf ", varNames[0], coefficientMetrics[numVar * 6 - 6]);
  for(i = 0; i < numVar - 1; i++) {
    printf("%+.2lf %s ",coefficientMetrics[6 *i], varNames[i+1]);
  }
  printf("\n");
  i = 0;

  /*
  *Print model metrics
  */
  printf("\n  Source |       SS       df       MS                  Number of obs = %7d", (int)modelMetrics[0]);
  printf("\n---------+------------------------------               F(%3d,%6d) =  %6.5g", (int)modelMetrics[1], (int)modelMetrics[2], modelMetrics[3]);
  printf("\n   Model |  %10.9g %5d  %10.9g               Prob > F      =  %6.4lf", modelMetrics[4], (int)modelMetrics[5], modelMetrics[6], modelMetrics[7]);
  printf("\nResidual |  %10.9g %5d  %10.9g               R-squared     =  %6.4lf", modelMetrics[8], (int)modelMetrics[9], modelMetrics[10], modelMetrics[11]);
  printf("\n---------+------------------------------               Adj R-squared =  %6.4lf", modelMetrics[12]);
  printf("\n   Total |  %10.9g %5d  %10.9g               Root MSE      =  %6.5g\n", modelMetrics[13], (int)modelMetrics[14], modelMetrics[15], modelMetrics[16]);

  /*
  *Print coefficient metrics
  */
  printf("\n------------------------------------------------------------------------------");
  printf("\n ");
  for(j = 0; j < 8; j++) {
    if(varNames[i][j] == '\0') {
      for( ; j < 8; j++) printf(" ");
    }
    else printf("%c",varNames[i][j]); 
  }
  i++;
  printf("|      Coef.   Std. Err.       t     P>|t|       [95%% Conf. Interval]");
  printf("\n---------+--------------------------------------------------------------------");

  for (k = 0; k < (int)modelMetrics[1]; k++) {
    printf("\n ");
    for(j = 0; j < 8; j++) {
      if(varNames[i][j] == '\0') {
        for( ; j < 8; j++) printf(" ");
      }
      else printf("%c",varNames[i][j]); 
    }
  i++;
  printf("|%11.7g %11.7g %9.5g   %4.3lf  %13.7g  %10.7g", coefficientMetrics[k*6],  coefficientMetrics[k*6 + 1],  coefficientMetrics[k*6 + 2],  coefficientMetrics[k*6 +3],  coefficientMetrics[k*6 + 4],  coefficientMetrics[k*6 +5]);
  }
  k++;
  printf("\n Const   |%11.7g %11.7g %9.5g   %4.3lf  %13.7g  %10.7g", coefficientMetrics[numVar * 6 - 6],  coefficientMetrics[numVar * 6 - 5],  coefficientMetrics[numVar * 6 - 4],  coefficientMetrics[numVar * 6 - 3],  coefficientMetrics[numVar * 6 - 2],  coefficientMetrics[numVar * 6 - 1]);
  printf("\n------------------------------------------------------------------------------\n");
}
