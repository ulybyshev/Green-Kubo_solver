#include "kernel.h"
#include "basic_structures.h"
#include "constants.h"

//here new kernels can be added, they are choosen by kernel_switcher


//  \pi is already defined as PI constant
 

//general wrapper
//euclidean time=0.....beta
//program assumes that beta is in units of the step in euclidean time
double kernel_function(double omega, double beta, double euclidean_time)
{

    switch (kernel_switcher)
    {
	case 0: 
	    return kernel_conductivity(omega, beta, euclidean_time);
	case 1:
	    return kernel_DOS(omega, beta, euclidean_time);
	default:
	    return kernel_conductivity(omega, beta, euclidean_time);
    }
}


//set of kernels
//to extract conductivity from current-current correlator
double kernel_conductivity(double omega, double beta, double euclidean_time)
{
    return (omega/PI)*(exp(omega*(-euclidean_time))+exp(omega*(euclidean_time-beta)))/(1.0-exp(-omega*beta));
}

//to extract Density of States from fermionic propagator
double kernel_DOS(double omega, double beta, double euclidean_time)
{
    return (exp(omega*(-euclidean_time))+exp(omega*(euclidean_time-beta)))/(1.0+exp(-omega*beta));
}
