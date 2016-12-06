#include "constants.h"

//computation of relative error of spectral function for given regularization constant lambda
//lambda is global variable
double relative_error_computation(correlator* pC, calc_structures* pA);
double cov_reg_lambda_definition(correlator* pC, calc_structures* pA, int* flag_limit, FILE* general_log);
double svd_reg_lambda_definition(correlator* pC, calc_structures* pA, int* flag_limit, FILE* general_log);

//final computation of the spectral function
void delta_rho_calculation_and_output(correlator * pC, calc_structures* pA, FILE* file_out_excl, int flag_mode=0, spectral_functions* pRho=NULL, int lambda_count=0);


