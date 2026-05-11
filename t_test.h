// https://statisticsbyjim.com/hypothesis-testing/how-to-find-p-value/
#include <math.h>

double t_value(double sample_mean, double null_value, double stddev, int sample_size)
{
	return (sample_mean - null_value) / (stddev / sqrt(sample_size));
}
