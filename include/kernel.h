#include <math.h>

//general wrapper
//euclidean time=0.....beta
//program assumes that beta is in units of the step in euclidean time
double kernel_function(double omega, int beta, int euclidean_time);

//set of kernels
double kernel_conductivity(double omega, int beta, int euclidean_time);

double kernel_DOS(double omega, int beta, int euclidean_time);

//lattice version of kernel_DOS (takes into account time discretization)
double kernel_lattice_DOS(double omega, int  beta, int euclidean_time);

//lattice exponent with discrete time
double lattice_exp (double omega, int euclidean_time);

double my_pow(double x, int b);
