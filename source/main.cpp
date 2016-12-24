#include "constants.h"
#include "parser_cmd_line.h"
#include "parser_const_file.h"
#include "input_output.h"
#include "math_functions.h"
#include "spectral_func.h"

void print_test(correlator *pC) {
  int i, j;
   for(i=0;i<num_jack_samples;i++) {
     for(j=0;j<pC[i].N_valid_points;j++) {
       printf("DEBUG %d %d %e %e\n",i,j,pC[i].corr_full[j],pC[i].error_full[j]);
     }
   }
}


int main(int argc, char ** argv)
{
    int i,j;
    flag_log_output=false;
    special_flag_log_output=false;
    correlator tempC;
    correlator *C;
    initial_data_description D;
    char file_name[1000];

    if(parse_cmd_line(argc, argv))
    {
	printf("command line is ok\n");
    }
    else
    {
	printf("something is wrong with command line options\n");
	return 0;
    }

    set_default_values();

    if(flag_constants_file)
    {
	FILE* parameters_file;
	parameters_file=fopen(parameters_filename,"r");
	parse_const_file(parameters_file, &tempC);
	fclose(parameters_file);
    }
    else
    {
	tempC.format(Nt_2, Nt_2);
	for(i=0;i<tempC.N_valid_points; i++)
	    tempC.points_numbers[i]=i+1;
    }


    FILE* general_log;
    general_log=fopen_control("general_log.txt","w");
    print_parameters(general_log, &tempC);        

    //files for input-output operations
    FILE* file_in_current;
  
    file_in_current=fopen(correlator_filename,"r");
    
    if(flag_jackknife) 
    {
      input_raw_data(file_in_current); 
      if(n_conf<4)
	{
	    fprintf(general_log,"ERROR: too small amount of configurations (<4). Exit.\n"); fflush(general_log);
	    return 0;
	}
      if(flag_tune_blocking)
      {
            if(!input_data_analysis(&D))
    	    {
        	fprintf(general_log,"WARNING: too small amount of configurations to make appropriate data binning\n");fflush(general_log);
    	    }
    	    else
    	    {
    		fprintf(general_log,"Binning was successful\n");fflush(general_log);
    		//output for correlation lengths
    		{
		    int count_time;
		    for(count_time=0; count_time<D.N_histories;count_time++)
		    {
			fprintf(general_log,"autocorrelation at time %d is equal to %.15le\n", D.times[count_time], D.corr_lengths[count_time]);fflush(general_log);
		    }
		    fprintf(general_log, "maximum autocorrelation time is  %.15le\n",D.largest_corr_length);fflush(general_log);
		}
    	    }
      }
      fprintf(general_log,"Number of bins in data blocking \t %d\n", num_jack_samples );fflush(general_log);
      if(num_jack_samples==N_BINS_MINIMUM && flag_tune_blocking)
      {
    	fprintf(general_log,"WARNING:not enough amount of data or data is too autocorrelated\n");fflush(general_log);
      }
    }
    else {
      FILE* file_in_matrix=fopen(cov_matrix_filename,"r");
      input_correlator_matrix(file_in_current, file_in_matrix, &tempC);
      fclose(file_in_matrix);
    }
    fclose(file_in_current);

    if(flag_jackknife) {
      C=(correlator *)calloc(num_jack_samples,sizeof(correlator));
      for(i=0;i<num_jack_samples;i++) {
	C[i].format(tempC.N_full_points,tempC.N_valid_points);
	for(j=0;j<tempC.N_valid_points;j++) {
	  C[i].points_numbers[j]=tempC.points_numbers[j];
	}
	C[i].construct_intervals();
	get_jack_sample(&C[i], i+1);
      }
    }

    
//now data for correlator and covariance matrix are ready
    calc_structures A(tempC.N_valid_points);
        
//integration to obtain R vector
    if(flag_jackknife) {
      R_integration( &A, &C[0]);
      omega_R_integration( &A, &C[0]);
    }
    else {
      R_integration( &A, &tempC);
      omega_R_integration( &A, &tempC);
    }

//calculation of integrals in W matrix

    fprintf(general_log,"\n W matrix calculation started\nN_center=%d\n", A.N_center);fflush(general_log);
    int count_center=0;
    for(count_center=0; count_center<A.N_center; count_center++)
    {
	fprintf(general_log,"count_center=%d\n", count_center);fflush(general_log);
	if(flag_jackknife)
	  W_integration(A.W[count_center], &C[0],  A.center[count_center]/(2.0*C[0].length));
	else
	  W_integration(A.W[count_center], &tempC,  A.center[count_center]/(2.0*tempC.length));
    }
    fprintf(general_log,"\n W matrix calculation finished\n");fflush(general_log);

    double lambda_final;
    int flag_limit;
    if(!flag_jackknife) {
      
      if(flag_lambda_regularization<0) {
	if(flag_lambda_regularization==-1)
	  lambda_final=cov_reg_lambda_definition(&tempC, &A,  &flag_limit, general_log);
	else
	  lambda_final=svd_reg_lambda_definition(&tempC, &A,  &flag_limit, general_log);
	
	fprintf(general_log,"\n\nFinal lambda=%.15le\nfinal_flag_limit=%d\n", lambda_final, flag_limit);fflush(general_log);
	lambda=lambda_final;
      }    
    
      //special_flag_log_output=true;
      
      //now the final calculation of spectral function is launched for the found value of lambda
      FILE* file_out_rho;
      //clean file before output
      file_out_rho=fopen_control("rho_basic.txt", "w");
      delta_rho_calculation_and_output(&tempC, &A, file_out_rho, 0);
      fclose(file_out_rho); 
      if(flag_lambda_regularization<0) {
	special_flag_log_output=false;
	
	
	//calculation of specral function for several values of lambda in the vicinity of the found one
	spectral_functions Rho(Number_lambda_points, lambda_final, &A);
	fprintf(general_log,"\n\nStudy of regularization influence\n");fflush(general_log);
	int count_lambda=0;
	for(count_lambda=0;count_lambda<Rho.N_lambda; count_lambda++) {
	  fprintf(general_log,"lambda=%.15le\n",Rho.lambda_array[count_lambda]);fflush(general_log);
	  char rho_filename[1000];
	  sprintf(rho_filename,"rho_lambda_%.15le.txt", Rho.lambda_array[count_lambda]);
	  file_out_rho=fopen_control(rho_filename, "w");
	  lambda=Rho.lambda_array[count_lambda];
	  delta_rho_calculation_and_output(&tempC, &A, file_out_rho, 1, &Rho, count_lambda);
	  fclose(file_out_rho);
	}
	//final output of dependence of spectral function on regularization
	FILE* rho_lambda;
	rho_lambda=fopen_control("rho_lambda.txt","w");
	fprintf(rho_lambda,"#value of lambda\t");fflush(rho_lambda);
	for(count_center=0;count_center<A.N_center;count_center++) {
	  fprintf(rho_lambda,"Rho for c=%.3le T\tstat error\t",A.center[count_center]);fflush(rho_lambda);
	}
	fprintf(rho_lambda,"\n");fflush(rho_lambda);
	for(count_lambda=0;count_lambda<Rho.N_lambda; count_lambda++) {
	  fprintf(rho_lambda,"%.15le\t", Rho.lambda_array[count_lambda] );fflush(rho_lambda);
	  
	  for(count_center=0;count_center<A.N_center;count_center++) {
	    fprintf(rho_lambda,"%.15le\t%.15le\t",Rho.rho_array[count_lambda][count_center], Rho.rho_err_array[count_lambda][count_center]);fflush(rho_lambda);
	  }
	  fprintf(rho_lambda,"\n");fflush(rho_lambda);
	}
	fclose(rho_lambda);
      }
      printf("Ok\n");
    }//if(!flag_jackknife)
    else {
      
      if(flag_lambda_regularization<0) {
	if(flag_lambda_regularization==-1)
	  lambda_final=cov_reg_lambda_definition_jack(C, &A,  &flag_limit, general_log);
	else
	  lambda_final=svd_reg_lambda_definition_jack(C, &A,  &flag_limit, general_log);
	
	fprintf(general_log,"\n\nFinal lambda=%.15le\nfinal_flag_limit=%d\n", lambda_final, flag_limit);fflush(general_log);
	lambda=lambda_final;
      }
    
      special_flag_log_output=false;
      flag_log_output=false;
      
      //now the final calculation of spectral function is launched for the found value of lambda
      delta_rho_calculation_and_output_jack(C, &A, 0);
      //return 0;

      if(flag_lambda_regularization<0) {
	special_flag_log_output=false;
	
	//calculation of specral function for several values of lambda in the vicinity of the found one
	spectral_functions *Rho=(spectral_functions *)calloc(num_jack_samples,sizeof(spectral_functions));
	spectral_functions Rho_avg(Number_lambda_points, lambda_final, &A);
	
	for(i=0;i<num_jack_samples;i++) {
	  Rho[i].format(Number_lambda_points, lambda_final, &A);
	}
	
	fprintf(general_log,"\n\nStudy of regularization influence\n");fflush(general_log);
	int count_lambda=0;
	for(count_lambda=0;count_lambda<Rho[0].N_lambda; count_lambda++) {
	  fprintf(general_log,"lambda=%.15le\n",Rho[0].lambda_array[count_lambda]);fflush(general_log);
	  lambda=Rho[0].lambda_array[count_lambda];
	  delta_rho_calculation_and_output_jack(C, &A, 1, Rho, count_lambda, &Rho_avg);
	}

	//final output of dependence of spectral function on regularization
	FILE* rho_lambda;
	
	  sprintf(file_name,"rho_lambda.txt");
	  rho_lambda=fopen_control(file_name,"w");
	  fprintf(rho_lambda,"#value of lambda\t");fflush(rho_lambda);
	  for(count_center=0;count_center<A.N_center;count_center++) {
	    fprintf(rho_lambda,"Rho for c=%.3le T\tstat error\t",A.center[count_center]);fflush(rho_lambda);
	  }
	  fprintf(rho_lambda,"\n");fflush(rho_lambda);
	  for(count_lambda=0;count_lambda<Rho_avg.N_lambda; count_lambda++) {
	    fprintf(rho_lambda,"%.15le\t", Rho_avg.lambda_array[count_lambda] );fflush(rho_lambda);
	    
	    for(count_center=0;count_center<A.N_center;count_center++) {
	      fprintf(rho_lambda,"%.15le\t%.15le\t",Rho_avg.rho_array[count_lambda][count_center], Rho_avg.rho_err_array[count_lambda][count_center]);fflush(rho_lambda);
	    }
	    fprintf(rho_lambda,"\n");fflush(rho_lambda);
	  }
	  fclose(rho_lambda);
	
	free(Rho);
      }
      printf("Ok\n");   
    }
    
    fclose(general_log);
    if(flag_jackknife) {
      free(raw_data);
      free(C);
    }
    return 0;
}
