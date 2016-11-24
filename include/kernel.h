#include <math.h>

//general wrapper
//euclidean time=0.....beta
//program assumes that beta is in units of the step in euclidean time
double kernel_function(double omega, double beta, double euclidean_time);

//set of kernels
double kernel_conductivity(double omega, double beta, double euclidean_time);

double kernel_DOS(double omega, double beta, double euclidean_time);
