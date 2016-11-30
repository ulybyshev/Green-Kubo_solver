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
	    return kernel_DOS(omega, beta, euclidean_time);
	case 2:
	    return kernel_lattice_DOS(omega, beta, euclidean_time);
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

//to extract Density of States from fermionic propagator
//written using assumption that the energy unit for density of states is he inverse step in Euclidean time
//and normalization condition for the spectral function is \int_0^\infty \frac {\rho(\omega) d\omega } = 1/2
//  (we assume the \rho (-omega) = \rho (\omega)) so the full integral  \int_{-infty}^\infty \frac {\rho(\omega) d\omega } = 1
double kernel_DOS(double omega, int beta, int  euclidean_time)
{
    return (exp(omega*(-(double)euclidean_time))+exp(omega*((double)euclidean_time-(double)beta)))/(1.0+exp(-omega*(double)beta));
}

//lattice version of kernel_DOS (takes into account time discretization)
double kernel_lattice_DOS(double omega, int beta, int  euclidean_time)
{
    return (lattice_exp(omega,(-euclidean_time))+lattice_exp(omega,(euclidean_time-beta)))/(1.0+lattice_exp(omega,-beta));
}



//lattice exponent with discrete time
double lattice_exp(double omega, int euclidean_time)
{
    if(euclidean_time>0.0)
    {
	if (omega>-1.0)
	    return my_pow(1.0+omega, euclidean_time);
	else
	    return 0.0;
    }
    else
    {
	if(omega<1.0)
{
//	    printf("%.5le\t%d\t%.5le\t%.5le\n",omega,euclidean_time,exp(omega*(double)euclidean_time) , my_pow(1.0-omega, -euclidean_time));
	    return my_pow(1.0-omega, -euclidean_time);//exp(omega*euclidean_time);//return pow(1.0-omega, euclidean_time);

}
	else
	    return 0.0;
    }
}


double my_pow(double x, int b)
{
    int i=0;
    double res=1.0;
    for(i=0;i<b;i++)
	res*=x;
    return res;
}
