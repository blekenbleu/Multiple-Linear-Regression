#define _CRT_SECURE_NO_DEPRECATE            // fopen(), fscanf() warnings
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"

// process files from CAcorr

int gpdata(FILE *fp, char **varNames, int columns, int gp, Matrix *x, Matrix *y)
{
	if (10 != columns)
	{
		printf("gpdata():  expected 10 columns; got %d\n", columns);
		return 0;
	}
	size_t size = x->rows = y->rows = gp;
	x->cols = 7;
	y->cols = 1;
	char buffer[256] = { '\0' };
	double* dd[4] = {
		(double*)calloc(gp, sizeof(double)),
		(double*)calloc(gp, sizeof(double)),
		(double*)calloc(gp, sizeof(double)),
		(double*)calloc(gp, sizeof(double))
	};
	double  *dp[4] = { dd[0],dd[1], dd[2], dd[3] },
			*xp = x->data = (double*)calloc(gp * 7, sizeof(double)),
			d0, d1, d2, d3, d4, d5;

	if (NULL == xp || NULL == dd[0] || NULL == dd[1] || NULL == dd[2] || NULL == dd[3])
	{
		printf("gpdata():  calloc() fail\n");
		return 0;
	}

	// rearrange varNames
	char *depNames[4] = { varNames[2], varNames[3], varNames[4], varNames[5] };
	varNames[2] = varNames[1];
	varNames[1] = varNames[0];
	varNames[3] = varNames[6];
	varNames[4] = varNames[7];
	varNames[5] = varNames[8];
	varNames[6] = varNames[9];
	varNames[0] = depNames[0];

	y->data = dd[0];
	int i, rc = 10;
	for (i = 0; !feof(fp) && i < gp && 10 == rc; i++)
	{
		*xp++ = 1.0;
		int rc = fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
			&d0, &d1, dp[0]++, dp[1]++, dp[2]++, dp[3]++, &d2, &d3, &d4, &d5);
		*xp++ = d0;
		*xp++ = d1;
		*xp++ = d2;
		*xp++ = d3;
		*xp++ = d4;
		*xp++ = d5;
	}

	return gp;
}

