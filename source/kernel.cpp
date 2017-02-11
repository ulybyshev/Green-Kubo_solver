#include "kernel.h"
#include "basic_structures.h"
#include "constants.h"

//here new kernels can be added, they are choosen by kernel_switcher


//  \pi is already defined as PI constant
 

//general wrapper
//euclidean time=0.....beta

//IMPORTANT: about dimensions !
//program assumes that beta is in units of the step in euclidean time and omega in units of inverse step in Euclidean time
//the same is for dimension of resulting spectral function and input correlator:  energy unit for it is defined by the the step in Euclidean time
double kernel_function(double omega, int beta, int euclidean_time)
{

    switch (kernel_switcher)
    {
	case 0: 
	    return kernel_conductivity(omega, beta, euclidean_time);
	case 1:
	    return kernel_DOS_even(omega, beta, euclidean_time);
	case 2:
	    return kernel_lattice_DOS_even(omega, beta, euclidean_time);
	case 3:
	    return kernel_meson(omega, beta, euclidean_time);
	default:
	    return kernel_conductivity(omega, beta, euclidean_time);
    }
}


//set of kernels
//to extract conductivity from current-current correlator (written using assumption that the energy unit for current-current correlator is the inverse step in Euclidean time)
double kernel_conductivity(double omega, int  beta,  int euclidean_time)
{
    return (omega/PI)*(exp(omega*(-(double)euclidean_time))+exp(omega*((double)euclidean_time-(double)beta)))/(1.0-exp(-omega*(double)beta));
}


//to extract spectral function for meson: the same as for conductivity but with tanh(omega*beta/2) regularization at zero omega
//it means that the kernel is cosh(owega*(tau-beta/2))/cosh(omega*beta/2),
// but the spectral function should be multiplied in the end by sinh(omega*beta/2)/cosh(omega*beta/2) = sinh (omega/2)/cosh (omega/2) if omega is in units of temperature
double kernel_meson(double omega, int  beta,  int euclidean_time)
{
    return   (exp(omega*(-(double)euclidean_time))+exp(omega*((double)euclidean_time-(double)beta)))/(1.0+exp(-omega*(double)beta));
}


//to extract Density of States from fermionic propagator
//written using assumption that the energy unit for density of states is he inverse step in Euclidean time
//and normalization condition for the spectral function is \int_0^\infty \frac {\rho(\omega) d\omega } = 1/2
//  (we assume the \rho (-omega) = \rho (\omega)) so the full integral  \int_{-infty}^\infty \frac {\rho(\omega) d\omega } = 1
double kernel_DOS_even(double omega, int beta, int  euclidean_time)
{
    return (exp(omega*(-(double)euclidean_time))+exp(omega*((double)euclidean_time-(double)beta)))/(1.0+exp(-omega*(double)beta));
}

//for non-symmetrical case
double kernel_DOS_odd(double omega, int beta, int euclidean_time)
{
    return (exp(omega*(-(double)euclidean_time))-exp(omega*((double)euclidean_time-(double)beta)))/(1.0+exp(-omega*(double)beta));
}

//lattice version of kernel_DOS (takes into account time discretization)
double kernel_lattice_DOS_even(double omega, int beta, int  euclidean_time)
{
    double lat_exp_beta=lattice_exp(omega, beta);
    double lat_exp_beta1=lattice_exp(-omega, beta);
    return  (lattice_exp(-omega, euclidean_time)* (1.0+lat_exp_beta) + lattice_exp(omega, euclidean_time)* (1.0+lat_exp_beta1)) / ( ( 1.0+ lat_exp_beta1 ) *(1.0+ lat_exp_beta ) );
}

//lattice version of kernel_DOS (takes into account time discretization)
//odd version of the kernel for non-symmetrical case
double kernel_lattice_DOS_odd(double omega, int beta, int  euclidean_time)
{
    double lat_exp_beta=lattice_exp(omega, beta);
    double lat_exp_beta1=lattice_exp(-omega, beta);
    return (lattice_exp(-omega, euclidean_time)* (1.0+lat_exp_beta) - lattice_exp(omega, euclidean_time)* (1.0+lat_exp_beta1))/ ( ( 1.0+ lat_exp_beta1 ) *(1.0+ lat_exp_beta ) );
}

