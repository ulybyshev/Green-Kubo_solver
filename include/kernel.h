#include <math.h>
#include "math_functions.h"

//general wrapper
//euclidean time=0.....beta
//program assumes that beta is in units of the step in euclidean time
double kernel_function(double omega, int beta, int euclidean_time);

//set of kernels
double kernel_conductivity(double omega, int beta, int euclidean_time);

double kernel_DOS_even(double omega, int beta, int euclidean_time);

//lattice version of kernel_DOS (takes into account time discretization)
double kernel_lattice_DOS_even(double omega, int  beta, int euclidean_time);

//verions of the kernel for DOS in case of non-symmtrical correlator
double kernel_DOS_odd(double omega, int beta, int euclidean_time);

//lattice version of kernel_DOS (takes into account time discretization)
double kernel_lattice_DOS_odd(double omega, int  beta, int euclidean_time);
