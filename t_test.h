// https://statisticsbyjim.com/hypothesis-testing/how-to-find-p-value/
#include <math.h>

double t_value(double sample_mean, double null_value, double stddev, int sample_size)
{
	return (sample_mean - null_value) / (stddev / sqrt(sample_size));
}

// two-tailed significance 0.05 https://statisticsbyjim.com/hypothesis-testing/t-distribution-table/
float t_table[] = { 12.71f, 4.303f, 3.182f, 2.776f, 2.571f, 2.447f, 2.365f, 2.306f, 2.262f,
	2.228f, 2.201f, 2.179f, 2.160f, 2.145f, 2.131f, 2.120f, 2.110f, 2.101f, 2.093f, 2.086f,
	2.080f, 2.074f, 2.069f, 2.064f, 2.060f, 2.056f, 2.052f, 2.048f, 2.045f, 2.042f, 2.021f,
	2.000f, 1.990f, 1.984f, 1.962f, 1.960f };

float critical_value(unsigned int df)
{
	if (30 >= df)
		return t_table[df - 1];
	else if (df >= 1000)
		return t_table[34];
	else if (df >= 100)
		return t_table[33] + (df - 100) * (t_table[34] - t_table[33]) / 900;
	else if (df >= 80)
		return t_table[33] + (df - 80) * (t_table[33] - t_table[32]) / 20;
	else if (df >= 60)
		return t_table[33] + (df - 60) * (t_table[32] - t_table[31]) / 20;
	else if (df >= 40)
		return t_table[33] + (df - 40) * (t_table[31] - t_table[30]) / 20;
	else
		return t_table[29] + (df - 30) * (t_table[30] - t_table[29]) / 10;
}
