#include "constants.h"
#include "parser_cmd_line.h"
#include "parser_const_file.h"
#include "input_output.h"
#include "math_functions.h"


int main(int argc, char ** argv)
{

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
  FILE* file_out;
  FILE* file_out_excl;


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
    int count_center=0, i;
    for(count_center=0; count_center<A.N_center; count_center++)
    {
	W_integration(A.W[count_center], &C,  A.center[count_center]/(2.0*C.length));
    }




//now the calculation of conductivity is launched
{
    //clean file before output
    file_out_excl=fopen_control("rho_excl.txt", "w");  
    fclose(file_out_excl);
    file_out=fopen_control("rho.txt", "w");  
    fclose(file_out);
    

    double B, B0;
    
    gsl_vector * Q;
    gsl_vector* Q_real;
    gsl_vector* Q_initial;
    double center, omega;
    
    Q_initial=gsl_vector_calloc(C.N_valid_points);
  
  
//    for(center=center_start; center<=center_stop; center+=center_delta)
    for(count_center=0; count_center<A.N_center; count_center++)
    {
	double center_real=A.center[count_center]/(2.0*C.length);
	char file_name[1000];
	sprintf(file_name,"delta_function_c=%3.3leT.txt", A.center[count_center]);
        Q=gsl_vector_calloc(C.N_valid_points);
	Q_real=gsl_vector_calloc(C.N_valid_points);

	calculate_Q(Q, &A, &C, count_center);

	if(count_center==0)
	{
    	    for(i=0;i<C.N_valid_points;i++)
    		gsl_vector_set(Q_initial, i, gsl_vector_get(Q,i));
	}
	if(count_center==0)
	{
    	    for(i=0;i<C.N_valid_points;i++)
    		gsl_vector_set(Q_real, i, gsl_vector_get(Q,i));
	}
	else
	{
    	    if(flag_exclude_delta==1 && count_center>=count_start_exclude)
    	    {
    		FILE* log_exclusion;
    		log_exclusion=fopen_log("log_exclusion.txt","a", center_real);
    		B=delta0(Q, &C);
    		B0=delta0(Q_initial, &C);
    		fprintf(log_exclusion,"B=%.15le\nC0=%.15le\n",B,B0);fflush(log_exclusion);
        
    		for(i=0;i<C.N_valid_points;i++)
        	    gsl_vector_set(Q_real, i, gsl_vector_get(Q,i)-(B/B0)*gsl_vector_get(Q_initial,i));
        
    		//normalization
    		double new_norma=0.0;
    		for(i=0;i<C.N_valid_points;i++)
    		{
        	    new_norma += gsl_vector_get(Q_real, i) * gsl_vector_get(A.R,i);
    		}
    		fprintf(log_exclusion, "norma_old=%.15le\n", new_norma);fflush(log_exclusion);
        
    		for(i=0;i<C.N_valid_points;i++)
    		{
        	    gsl_vector_set(Q_real,i,gsl_vector_get(Q_real, i)/new_norma );
    		}
    		new_norma=0.0; 
    		for(i=0;i<C.N_valid_points;i++)
    		{
        	    new_norma += gsl_vector_get(Q_real, i) * gsl_vector_get(A.R,i);
    		}
    		fprintf(log_exclusion, "norma_new=%.15le\n", new_norma);fflush(log_exclusion);
        
    		fclose(log_exclusion);     
    	    }
    	    else
    	    {
    		for(i=0;i<C.N_valid_points;i++)
    		    gsl_vector_set(Q_real, i, gsl_vector_get(Q,i));
    	    }
	}
	//delta function output
	file_out=fopen_control(file_name,"w");
	for(omega=0;omega<omega_plot_limit/(2.0*C.length);omega+=omega_plot_delta/(2.0*C.length))
	{
    	    fprintf(file_out,"%.15le\t%.15le\t%.15le\n", omega*2.0*C.length, delta(omega,Q, &C), delta(omega, Q_real, &C));
    	    fflush(file_out);
	}
	fclose(file_out);

    
	//output of dimensionless spectral function
	double rho, rho_stat_err, width;
	double rho_real, rho_stat_err_real, width_real;
	//values to really characterize delta functions 
	double start, stop, center1, width1;
	double real_start, real_stop, real_center1, real_width1;


	file_out=fopen_control("rho.txt", "a");
	file_out_excl=fopen_control("rho_excl.txt", "a");    
	calculate_rho(Q, &C, &rho, &rho_stat_err);
	width=delta_width_calculation(Q, center_real, &C);

	delta_characteristics_calculation(&start, &stop, &center1, Q, center_real, &C);
	width1=(stop-start)/2.0;

	calculate_rho(Q_real,  &C, &rho_real, &rho_stat_err_real);
	width_real=delta_width_calculation(Q_real, center_real, &C);

	delta_characteristics_calculation(&real_start, &real_stop, &real_center1, Q_real, center_real, &C);
	real_width1=(real_stop-real_start)/2.0;


	fprintf(file_out,"%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\n", A.center[count_center], width*2.0*C.length, center1, width1, rho, rho_stat_err, start, stop);
	fprintf(file_out_excl,"%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\n", A.center[count_center], width_real*2.0*C.length, real_center1, real_width1,  rho_real, rho_stat_err_real, real_start, real_stop);


	fclose(file_out);
	fclose(file_out_excl);


	gsl_vector_free(Q);
	gsl_vector_free(Q_real);

    }

gsl_vector_free(Q_initial);
  
 
}

 

    printf("Ok\n");

    fclose(general_log);
    return 0;
}
