#include "constants.h"


//for lattice version of the kernel functions
double lattice_exp(double omega, int euclidean_time);
//for lattice exponent
double my_pow(double x, int b);



//simple kernel for time points (works in cases of flag_model==1)
double kernel_points(int i, correlator* pC,double omega);

//kernel at zero omega (to exclude zeros of denominator from the calculation)
double kernel0_points(int i, correlator* pC);

//kernel in case of averaging over intervals (workd s in kase of flag_model==2)
double kernel_intervals(int i, correlator* pC, double omega);

//kernel at zero omega (to exclude zeros of denominator from the calculation)
double kernel0_intervals(int i, correlator* pC);

//general wrapper
double kernel(int i, correlator* pC, double omega);

double kernel0(int i, correlator* pC);


//functions for integration
//simple kernel integration (integer parameter t should NOT be equal to 0)
double kernel_int (double x, void * params);


//omega*kernel integration (integer parameter t should NOT be equal to 0)
//we need  it to define real center of the delta-function
double omega_kernel_int (double x, void * params);


double delta(double omega, gsl_vector* Q, correlator* pC);//omega_center is defined in calculation of vector Q
double delta0(gsl_vector* Q, correlator* pC);//omega_center is defined in calculation of vector Q

//simple kernel integration (integer parameter t should NOT be equal to 0)
double W_int (double x, void * params);


double W_function (w_parameters params, double omega);


//integration to obtain R vector
void R_integration(calc_structures* pA,correlator* pC);

//integration to obtain omega_R vector
void omega_R_integration(calc_structures* pA,correlator* pC);

//calculation of Q vector
void calculate_Q(gsl_vector* Q, calc_structures* pA, correlator* pC, double center);


void calculate_rho(gsl_vector* Q, correlator* pC, double* rho, double* rho_stat_err);


void delta_characteristics_calculation(double* start, double* stop, double* center_real, gsl_vector* Q, double center, correlator* pC);

double delta_width_calculation(gsl_vector* Q, double center, correlator* pC);

double delta_sq_int (double x, void* params);


