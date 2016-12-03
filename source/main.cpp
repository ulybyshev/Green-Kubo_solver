#include "constants.h"
#include "parser_cmd_line.h"
#include "parser_const_file.h"
#include "input_output.h"
#include "math_functions.h"
#include "spectral_func.h"


int main(int argc, char ** argv)
{
    flag_log_output=false;
    special_flag_log_output=false;
    correlator C;

    if(parse_cmd_line(argc, argv))
    {
	printf("command line is ok\n");
    }
    else
    {
	printf("something is wrong with command line options\n");
	return 0;
    }


    FILE* parameters_file;
    parameters_file=fopen(parameters_filename,"r");
    parse_const_file(parameters_file, &C);
    fclose(parameters_file);

    FILE* general_log;
    general_log=fopen_control("general_log.txt","w");
    print_parameters(general_log, &C);
    
    

//files for input-output operations
  FILE* file_in_current;
  FILE* file_in_matrix;


//input and conversion of correlator and covariance matrix
    file_in_current=fopen(correlator_filename,"r");
    file_in_matrix=fopen(cov_matrix_filename,"r");

    input_correlator_matrix(file_in_current, file_in_matrix, &C);

    fclose(file_in_current);
    fclose(file_in_matrix);
    
//now data for correlator and covariance matrix are ready
    calc_structures A(C.N_valid_points);

//integration to obtain R vector
    R_integration( &A, &C);
    omega_R_integration( &A, &C);

//calculation of integrals in W matrix

    fprintf(general_log,"\n W matrix calculation started\nN_center=%d\n", A.N_center);fflush(general_log);
    int count_center=0;
    for(count_center=0; count_center<A.N_center; count_center++)
    {
	fprintf(general_log,"count_center=%d\n", count_center);fflush(general_log);
	W_integration(A.W[count_center], &C,  A.center[count_center]/(2.0*C.length));
    }
    fprintf(general_log,"\n W matrix calculation finished\n");fflush(general_log);

    double lambda_final;
    int flag_limit;
    if(flag_lambda_regularization==-1)
    {
	lambda_final=cov_reg_lambda_definition(&C, &A,  &flag_limit, general_log);
	fprintf(general_log,"\n\nFinal lambda=%.15le\nfinal_flag_limit=%d\n", lambda_final, flag_limit);fflush(general_log);
        lambda=lambda_final;
    }
    
    special_flag_log_output=true;
    
    //now the final calculation of spectral function is launched for the found value of lambda
    FILE* file_out_rho;
    //clean file before output
    file_out_rho=fopen_control("rho_basic.txt", "w");  
    delta_rho_calculation_and_output(&C, &A, file_out_rho, 0);
    fclose(file_out_rho); 

    printf("Ok\n");

    fclose(general_log);
    return 0;
}
