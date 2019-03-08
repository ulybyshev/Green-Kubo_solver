#include "spectral_func.h"
#include "math_functions.h"
#include "input_output.h"
#include "kernel.h"

//high-level routines for calculation of spectral function, resolution functions and defining regularization constants

//computation of relative error of spectral function for given regularization constant lambda
//lambda is global variable
double relative_error_computation(correlator* pC, calc_structures* pA)
{
    gsl_vector * Q;
    int count_center;

    double average_rho=0.0;
    double average_error=0.0;
    for(count_center=0; count_center<pA->N_center; count_center++)
    {
        Q=gsl_vector_calloc(pC->N_valid_points);
	double rho, rho_stat_err;
	calculate_Q(Q, pA, pC, count_center);
	calculate_rho(Q, pC, &rho, &rho_stat_err);
	average_rho+=fabs(rho);
	average_error+=rho_stat_err;
	gsl_vector_free(Q);
    }
    average_rho=average_rho/(double)(pA->N_center);
    average_error=average_error/(double)(pA->N_center);

    return average_error/average_rho;
}


//defines lambda to satisfy the realtive error in case of covariance matrix regularization
//flag_limit=-1 if algorithm stopped at smallest regularization, flag_lambda=1 if algorithm reached largest regularization
double cov_reg_lambda_definition(correlator* pC, calc_structures* pA, int* flag_limit, FILE* general_log)
{

    fprintf(general_log,"Lambda definition for covariance matrix regularization is started\n"); fflush(general_log);

//full scan of the lambda interval
//    lambda =1.0-accuracy/gsl_matrix_get(pC->S,0,0);
    lambda =1.0-accuracy;
    double lambda_limit=1.0-2.0*pow(10.0, (double)limit_power);
    double error;
    int count_lambda=0;
    *flag_limit=1;
    while(lambda>lambda_limit)
    {
	error=relative_error_computation(pC, pA);
	fprintf(general_log,"Lambda=%.15le\t relative error=%.15le\n",  lambda, error); fflush(general_log);
	if (error<relative_error)
	{
	    if (count_lambda==0)
		*flag_limit=-1;
	    else
		*flag_limit=0;
	    break;
	}
	lambda=1.0-(1.0-lambda)*10.0;
	count_lambda++;
    }
    if(lambda<lambda_limit)
    {
	lambda=lambda_limit;
    }

    fprintf(general_log,"flag_limit=%d\n",*flag_limit); fflush(general_log);
    if ((*flag_limit)!=0)
	return lambda;

//now exact computation of lambda
    double lambda1=lambda;
    double lambda2=1.0-0.1*(1.0-lambda1);

    fprintf(general_log,"\n****\nlambda1=%.15le\t lambda2=%.15le\n",  lambda1, lambda2); fflush(general_log);
    
    count_lambda=0;
    while(1)
    {
	lambda=(lambda2+lambda1)*0.5;
	error=relative_error_computation(pC, pA);
	if(error<relative_error)
	{
	    lambda1=lambda;
	}
	else
	{
	    lambda2=lambda;
	}
	fprintf(general_log,"lambda1=%.15le\t lambda2=%.15le\n",  lambda1, lambda2); fflush(general_log);
	fprintf(general_log,"lambda=%.15le\t relative error=%.15le\n",  lambda, error); fflush(general_log);
	if(fabs(error-relative_error)<0.1*relative_error || count_lambda>1000)
	    break;
	count_lambda++;
    }

    return lambda;
}


//defines lambda to satisfy the realtive error in case of some kind of SVD regularization
//flag_limit=-1 if algorithm stopped at smallest regularization, flag_lambda=1 if algorithm reached largest regularization
double svd_reg_lambda_definition(correlator* pC, calc_structures* pA, int* flag_limit, FILE* general_log)
{

    fprintf(general_log,"Lambda definition for SVD regularization is started\n"); fflush(general_log);

//full scan of the lambda interval
//    lambda =accuracy/gsl_matrix_get(pC->S,0,0);
    lambda =accuracy;
    double lambda_limit=2.0*pow(10.0, (double)limit_power);
    double error;
    int count_lambda=0;
    *flag_limit=1;
    while(lambda<lambda_limit)
    {
	error=relative_error_computation(pC, pA);
	fprintf(general_log,"Lambda=%.15le\t relative error=%.15le\n",  lambda, error); fflush(general_log);
	if (error<relative_error)
	{
	    if (count_lambda==0)
		*flag_limit=-1;
	    else
		*flag_limit=0;
	    break;
	}
	lambda=lambda*10.0;
	count_lambda++;
    }
    if(lambda>lambda_limit)
    {
	lambda=lambda_limit;
    }

    fprintf(general_log,"flag_limit=%d\n",*flag_limit); fflush(general_log);
    if ((*flag_limit)!=0)
	return lambda;

//now exact computation of lambda
    double lambda1=lambda;
    double lambda2=0.1*lambda1;

    fprintf(general_log,"\n****\nlambda1=%.15le\t lambda2=%.15le\n",  lambda1, lambda2); fflush(general_log);
    
    count_lambda=0;
    while(1)
    {
	lambda=(lambda2+lambda1)*0.5;
	error=relative_error_computation(pC, pA);
	if(error<relative_error)
	{
	    lambda1=lambda;
	}
	else
	{
	    lambda2=lambda;
	}
	fprintf(general_log,"lambda1=%.15le\t lambda2=%.15le\n",  lambda1, lambda2); fflush(general_log);
	fprintf(general_log,"lambda=%.15le\t relative error=%.15le\n",  lambda, error); fflush(general_log);
	if(fabs(error-relative_error)<0.1*relative_error || count_lambda>1000)
	    break;
	count_lambda++;
    }

    return lambda;
}





//flag_mode==0 - usual calculation (includes output of  the spectral function and resolution functions)
//flag_mode==1 - only output of spectral function is performed, alongside with filling arrays for spectral function
void delta_rho_calculation_and_output(correlator * pC, calc_structures* pA, FILE* file_out_excl, int flag_mode, spectral_functions* pRho, int lambda_count)
{
    FILE* file_out;
    int count_center, i;

    double B, B0;
    
    gsl_vector * Q;
    gsl_vector* Q_real;
    gsl_vector* Q_initial;
    double omega;
    
    Q_initial=gsl_vector_calloc(pC->N_valid_points);
  
    FILE* file_rho_fin;

    if(!flag_mode)
    {
	file_rho_fin=fopen_control("rho_final.txt","w");
	fprintf(file_rho_fin, "#center\t spectral_function\t error\t resolution_function_peak_position\n");fflush(file_rho_fin);   
	fprintf(file_out_excl,"#center   resolution_func_width_from_D  resolution_func_peak_pos  resolution_func_width_at_half_peak  rho  rho_err  resolution_func_start  resolution_func_stop\n");
    }
    else
    {
	fprintf(file_out_excl, "#center\t spectral_function\t error\n");fflush(file_out_excl);   
    }
    for(count_center=0; count_center<pA->N_center; count_center++)
    {
      
	double center_real=pA->center[count_center]/(2.0*pC->length);
	char file_name[1000];
	sprintf(file_name,"delta_function_c=%3.3leT.txt", pA->center[count_center]);
        Q=gsl_vector_calloc(pC->N_valid_points);
	Q_real=gsl_vector_calloc(pC->N_valid_points);

	calculate_Q(Q, pA, pC, count_center);
	
	if(count_center==0)
	{
    	    for(i=0;i<pC->N_valid_points;i++)
    		gsl_vector_set(Q_initial, i, gsl_vector_get(Q,i));
	}
	if(count_center==0)
	{
    	    for(i=0;i<pC->N_valid_points;i++)
    		gsl_vector_set(Q_real, i, gsl_vector_get(Q,i));
	}
	else
	{
    	    if(flag_exclude_delta==1 && count_center>=count_start_exclude)
    	    {
    		FILE* log_exclusion;
    		LOG_FILE_OPERATION(log_exclusion=fopen_log("log_exclusion.txt","a", center_real);)
    		B=delta0(Q, pC);
    		B0=delta0(Q_initial, pC);
    		LOG_FILE_OPERATION(fprintf(log_exclusion,"B=%.15le\nC0=%.15le\n",B,B0);fflush(log_exclusion);)
        
    		for(i=0;i<pC->N_valid_points;i++)
        	    gsl_vector_set(Q_real, i, gsl_vector_get(Q,i)-(B/B0)*gsl_vector_get(Q_initial,i));
        
    		//normalization
    		double new_norma=0.0;
    		for(i=0;i<pC->N_valid_points;i++)
    		{
        	    new_norma += gsl_vector_get(Q_real, i) * gsl_vector_get(pA->R,i);
    		}
    		LOG_FILE_OPERATION(fprintf(log_exclusion, "norma_old=%.15le\n", new_norma);fflush(log_exclusion);)
        
    		for(i=0;i<pC->N_valid_points;i++)
    		{
        	    gsl_vector_set(Q_real,i,gsl_vector_get(Q_real, i)/new_norma );
    		}
		
    		new_norma=0.0; 
    		for(i=0;i<pC->N_valid_points;i++)
    		{
        	    new_norma += gsl_vector_get(Q_real, i) * gsl_vector_get(pA->R,i);
    		}
		
    		LOG_FILE_OPERATION(fprintf(log_exclusion, "norma_new=%.15le\n", new_norma);fflush(log_exclusion);)
        
    		LOG_FILE_OPERATION(fclose(log_exclusion);)     
    	    }
    	    else
    	    {
    		for(i=0;i<pC->N_valid_points;i++)
    		    gsl_vector_set(Q_real, i, gsl_vector_get(Q,i));
    	    }
	}
	if(flag_mode==0)
	{
	    //delta function output
	    file_out=fopen_control(file_name,"w");
	    fprintf(file_out,"#omega(in units of T)\t resolution_function\n");
	    for(omega=0;omega<omega_plot_limit/(2.0*pC->length);omega+=omega_plot_delta/(2.0*pC->length))
	    {
		double f_d;
		if (flag_exclude_delta==1 && count_center>=count_start_exclude)
		    f_d=delta(omega,Q_real, pC);
		else
		    f_d=delta(omega, Q, pC);
		    
    		fprintf(file_out,"%.5le\t%.5le\n", omega*2.0*pC->length, f_d);
    		fflush(file_out);
	    }
	    fclose(file_out);
	    
	}
	
	//output of dimensionless spectral function
	double rho, rho_stat_err, width;
	double rho_real, rho_stat_err_real, width_real;
	//values to characterize delta functions really
	double start, stop, center1, width1;
	double real_start, real_stop, real_center1, real_width1;


	LOG_FILE_OPERATION(file_out=fopen_control("rho_wo_excl.txt", "a");)
	calculate_rho(Q, pC, &rho, &rho_stat_err);
	width=delta_width_calculation(Q, center_real, pC);
	
	delta_characteristics_calculation(&start, &stop, &center1, Q, center_real, pC);
	width1=(stop-start)/2.0;

	calculate_rho(Q_real,  pC, &rho_real, &rho_stat_err_real);
	width_real=delta_width_calculation(Q_real, center_real, pC);
	
	delta_characteristics_calculation(&real_start, &real_stop, &real_center1, Q_real, center_real, pC);
	real_width1=(real_stop-real_start)/2.0;

	if(flag_mode==1)
	{
		pRho->rho_array[lambda_count][count_center]=rho;
		pRho->rho_err_array[lambda_count][count_center]=rho_stat_err;
	}

	LOG_FILE_OPERATION(fprintf(file_out,"%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\n", pA->center[count_center], width*2.0*pC->length, center1, width1, rho, rho_stat_err, start, stop);)
	
	if(!flag_mode)
	    fprintf(file_out_excl,"%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\n", pA->center[count_center], width_real*2.0*pC->length, real_center1, real_width1,  rho_real, rho_stat_err_real, real_start, real_stop);
	else
	    fprintf(file_out_excl,"%.15le\t%.15le\t%.15le\n",pA->center[count_center],rho_real, rho_stat_err_real);

	
	if(!flag_mode)
	{
	    fprintf(file_rho_fin,"%.15le\t%.15le\t%.15le\t%.15le\n",pA->center[count_center],rho_real, rho_stat_err_real, real_center1); fflush(file_rho_fin);
	}
	
	LOG_FILE_OPERATION(fclose(file_out);)

	gsl_vector_free(Q);
	gsl_vector_free(Q_real);

    }

gsl_vector_free(Q_initial);
	if(!flag_mode)
	{
	    fclose(file_rho_fin);
	} 
}

double relative_error_computation_jack(correlator* pC, calc_structures* pA)
{
    gsl_vector * Q;
    int i, count_center;
    double total_average_rho=0.0;
    double total_average_error=0.0;
    double rho_stat_err;
    double *rho=(double *)calloc(num_jack_samples,sizeof(double));
    double jack_avg, jack_err;

    for(count_center=0; count_center<pA->N_center; count_center++) {
      jack_avg=jack_err=0.0;
      for(i=0;i<num_jack_samples;i++) {
        Q=gsl_vector_calloc(pC[i].N_valid_points);
	calculate_Q(Q, pA, &pC[i], count_center);
	calculate_rho(Q, &pC[i], &rho[i], &rho_stat_err);
	total_average_rho+=fabs(rho[i]);
	jack_avg+=fabs(rho[i]);
	gsl_vector_free(Q);
      }
      jack_avg=jack_avg/(double)num_jack_samples;
      for(i=0;i<num_jack_samples;i++) {
	jack_err+=(rho[i]-jack_avg)*(rho[i]-jack_avg);
      }
      jack_err=sqrt(jack_err/((double)num_jack_samples*(num_jack_samples-1.0)));
      total_average_error+=jack_err;
    }
    total_average_rho=total_average_rho/(double)(pA->N_center*num_jack_samples);
    total_average_error=total_average_error/(double)(pA->N_center);

    free(rho);
    
    return total_average_error/total_average_rho;
}


//defines lambda to satisfy the realtive error in case of covariance matrix regularization
//flag_limit=-1 if algorithm stopped at smallest regularization, flag_lambda=1 if algorithm reached largest regularization
double cov_reg_lambda_definition_jack(correlator* pC, calc_structures* pA, int* flag_limit, FILE* general_log)
{

    fprintf(general_log,"Lambda definition for covariance matrix regularization is started\n"); fflush(general_log);

//full scan of the lambda interval
//    lambda =1.0-accuracy/gsl_matrix_get(pC[0].S,0,0);
    lambda =1.0-accuracy;
    double lambda_limit=1.0-2.0*pow(10.0, (double)limit_power);
    double error;
    int count_lambda=0;
    *flag_limit=1;
    while(lambda>lambda_limit)
    {
	error=relative_error_computation_jack(pC, pA);
	fprintf(general_log,"Lambda=%.15le\t relative error=%.15le\n",  lambda, error); fflush(general_log);
	if (error<relative_error)
	{
	    if (count_lambda==0)
		*flag_limit=-1;
	    else
		*flag_limit=0;
	    break;
	}
	lambda=1.0-(1.0-lambda)*10.0;
	count_lambda++;
    }
    if(lambda<lambda_limit)
    {
	lambda=lambda_limit;
    }

    fprintf(general_log,"flag_limit=%d\n",*flag_limit); fflush(general_log);
    if ((*flag_limit)!=0)
	return lambda;

//now exact computation of lambda
    double lambda1=lambda;
    double lambda2=1.0-0.1*(1.0-lambda1);

    fprintf(general_log,"\n****\nlambda1=%.15le\t lambda2=%.15le\n",  lambda1, lambda2); fflush(general_log);
    
    count_lambda=0;
    while(1)
    {
	lambda=(lambda2+lambda1)*0.5;
	error=relative_error_computation_jack(pC, pA);
	if(error<relative_error)
	{
	    lambda1=lambda;
	}
	else
	{
	    lambda2=lambda;
	}
	fprintf(general_log,"lambda1=%.15le\t lambda2=%.15le\n",  lambda1, lambda2); fflush(general_log);
	fprintf(general_log,"lambda=%.15le\t relative error=%.15le\n",  lambda, error); fflush(general_log);
	if(fabs(error-relative_error)<0.1*relative_error || count_lambda>1000)
	    break;
	count_lambda++;
    }

    return lambda;
}


//defines lambda to satisfy the realtive error in case of some kind of SVD regularization
//flag_limit=-1 if algorithm stopped at smallest regularization, flag_lambda=1 if algorithm reached largest regularization
double svd_reg_lambda_definition_jack(correlator* pC, calc_structures* pA, int* flag_limit, FILE* general_log)
{

    fprintf(general_log,"Lambda definition for SVD regularization is started\n"); fflush(general_log);

//full scan of the lambda interval
//    lambda =accuracy/gsl_matrix_get(pC[0].S,0,0);
    lambda =accuracy;
    double lambda_limit=2.0*pow(10.0, (double)limit_power);
    double error;
    int count_lambda=0;
    *flag_limit=1;
    while(lambda<lambda_limit)
    {
	error=relative_error_computation_jack(pC, pA);
	fprintf(general_log,"Lambda=%.15le\t relative error=%.15le\n",  lambda, error); fflush(general_log);
	if (error<relative_error)
	{
	    if (count_lambda==0)
		*flag_limit=-1;
	    else
		*flag_limit=0;
	    break;
	}
	lambda=lambda*10.0;
	count_lambda++;
    }
    if(lambda>lambda_limit)
    {
	lambda=lambda_limit;
    }

    fprintf(general_log,"flag_limit=%d\n",*flag_limit); fflush(general_log);
    if ((*flag_limit)!=0)
	return lambda;

//now exact computation of lambda
    double lambda1=lambda;
    double lambda2=0.1*lambda1;

    fprintf(general_log,"\n****\nlambda1=%.15le\t lambda2=%.15le\n",  lambda1, lambda2); fflush(general_log);
    
    count_lambda=0;
    while(1)
    {
	lambda=(lambda2+lambda1)*0.5;
	error=relative_error_computation_jack(pC, pA);
	if(error<relative_error)
	{
	    lambda1=lambda;
	}
	else
	{
	    lambda2=lambda;
	}
	fprintf(general_log,"lambda1=%.15le\t lambda2=%.15le\n",  lambda1, lambda2); fflush(general_log);
	fprintf(general_log,"lambda=%.15le\t relative error=%.15le\n",  lambda, error); fflush(general_log);
	if(fabs(error-relative_error)<0.1*relative_error || count_lambda>1000)
	    break;
	count_lambda++;
    }

    return lambda;
}


//flag_mode==0 - usual calculation (includes output of  the spectral function and resolution functions)
//flag_mode==1 - only output of spectral function is performed, alongside with filling arrays for spectral function
//here pC and pRho are pointers to an array of correlators and spectral functions respectively
void delta_rho_calculation_and_output_jack(correlator * pC, calc_structures* pA, int flag_mode, spectral_functions* pRho, int lambda_count, spectral_functions* Rho_avg)
{
    FILE* file_out;
    FILE* file_out_excl;
    FILE* file_out_avg;
    int count_center, i, jack;

    double B, B0;
    
    gsl_vector * Q;
    gsl_vector* Q_real;
    gsl_vector* Q_initial;
    double omega;
    double center_real;
    char file_name[1000];
    
    Q_initial=gsl_vector_calloc(pC[0].N_valid_points);
    Q=gsl_vector_calloc(pC[0].N_valid_points);
    Q_real=gsl_vector_calloc(pC[0].N_valid_points);

    double *jack_values=(double *)calloc(num_jack_samples*pA->N_center,sizeof(double));
    
    double* real_center1_array=(double*) calloc(pA->N_center, sizeof(double));
        
    if(flag_mode==0) {
      sprintf(file_name,"rho_final.txt");
    }
    else {
      sprintf(file_name,"rho_lambda_%.15le_avg.txt",pRho[0].lambda_array[lambda_count]);
    }
    file_out_avg=fopen_control(file_name,"w");
    fprintf(file_out_avg, "#center\t spectral_function\t error\t resolution_function_peak_position\n");fflush(file_out_avg);   
    
    for(jack=0;jack<num_jack_samples;jack++) {
      
      if(flag_mode==0) {
	sprintf(file_name,"rho_basic_%d.txt",jack+1);
	file_out_excl=fopen_control(file_name,"w");
	fprintf(file_out_excl,"#center   resolution_func_width_from_D  resolution_func_peak_pos  resolution_func_width_at_half_peak  rho  rho_err  resolution_func_start  resolution_func_stop\n");
      }
            
      for(count_center=0; count_center<pA->N_center; count_center++) {
	
	  center_real=pA->center[count_center]/(2.0*pC[0].length);
	  	  	  
	  calculate_Q(Q, pA, &pC[0], count_center);
	  	  
	  if(count_center==0) {
	    for(i=0;i<pC[0].N_valid_points;i++)
	      gsl_vector_set(Q_initial, i, gsl_vector_get(Q,i));
	  }
	  if(count_center==0) {
	    for(i=0;i<pC[0].N_valid_points;i++)
	      gsl_vector_set(Q_real, i, gsl_vector_get(Q,i));
	  }
	  else {
	    if(flag_exclude_delta==1 && count_center>=count_start_exclude) {
	      FILE* log_exclusion;
	      
	      LOG_FILE_OPERATION(log_exclusion=fopen_log("log_exclusion.txt","a", center_real);)
	      B=delta0(Q, &pC[0]);
	      B0=delta0(Q_initial, &pC[0]);
	      LOG_FILE_OPERATION(fprintf(log_exclusion,"B=%.15le\nC0=%.15le\n",B,B0);fflush(log_exclusion);)
	      
	      for(i=0;i<pC[0].N_valid_points;i++)
		gsl_vector_set(Q_real, i, gsl_vector_get(Q,i)-(B/B0)*gsl_vector_get(Q_initial,i));
	      
	      //normalization
	      double new_norma=0.0;
	      for(i=0;i<pC[0].N_valid_points;i++) {
		new_norma += gsl_vector_get(Q_real, i) * gsl_vector_get(pA->R,i);
	      }
	      LOG_FILE_OPERATION(fprintf(log_exclusion, "norma_old=%.15le\n", new_norma);fflush(log_exclusion);)
        
	      for(i=0;i<pC[0].N_valid_points;i++) {
		gsl_vector_set(Q_real,i,gsl_vector_get(Q_real, i)/new_norma );
	      }
	      new_norma=0.0; 
	      for(i=0;i<pC[0].N_valid_points;i++) {
		new_norma += gsl_vector_get(Q_real, i) * gsl_vector_get(pA->R,i);
	      }
	      LOG_FILE_OPERATION(fprintf(log_exclusion, "norma_new=%.15le\n", new_norma);fflush(log_exclusion);)
	    
	      LOG_FILE_OPERATION(fclose(log_exclusion);)     
	    }
	    else {
	      for(i=0;i<pC[0].N_valid_points;i++)
		gsl_vector_set(Q_real, i, gsl_vector_get(Q,i));
	    }
        
	  }
	  if(flag_mode==0) {
	    //delta function output
	    if(jack==0) {
	      sprintf(file_name,"delta_function_c=%3.3leT.txt", pA->center[count_center]);
	      file_out=fopen_control(file_name,"w");
	      fprintf(file_out,"#omega(in units of T)\t resolution_function\n");
	      for(omega=0;omega<omega_plot_limit/(2.0*pC[0].length);omega+=omega_plot_delta/(2.0*pC[0].length)) 
	      {
		double f_d;
		if (flag_exclude_delta==1 && count_center>=count_start_exclude)
		    f_d=delta(omega,Q_real, pC);
		else
		    f_d=delta(omega, Q, pC);
		    
		fprintf(file_out,"%.5le\t%.5le\n", omega*2.0*pC[0].length, f_d);
		fflush(file_out);
	      }
	      fclose(file_out);
	    }
	  }
	  
	//output of dimensionless spectral function
	double rho, rho_stat_err, width;
	double rho_real, rho_stat_err_real, width_real;
	//values to characterize delta functions really
	double start, stop, center1, width1;
	double real_start, real_stop, real_center1, real_width1;
	
	sprintf(file_name,"rho_wo_excl_%d.txt",jack+1);
	LOG_FILE_OPERATION(file_out=fopen_control(file_name, "a");)
	
	calculate_rho(Q, &pC[jack], &rho, &rho_stat_err);
	width=delta_width_calculation(Q, center_real, &pC[jack]);
		
	delta_characteristics_calculation(&start, &stop, &center1, Q, center_real, &pC[jack]);
	width1=(stop-start)/2.0;
	
	calculate_rho(Q_real,  &pC[jack], &rho_real, &rho_stat_err_real);
	width_real=delta_width_calculation(Q_real, center_real, &pC[jack]);
	      
	delta_characteristics_calculation(&real_start, &real_stop, &real_center1, Q_real, center_real, &pC[jack]);
	real_width1=(real_stop-real_start)/2.0;
		
	if(flag_mode==1) {
	  pRho[jack].rho_array[lambda_count][count_center]=rho;
	  pRho[jack].rho_err_array[lambda_count][count_center]=rho_stat_err;
	}
      
	LOG_FILE_OPERATION(fprintf(file_out,"%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\n", pA->center[count_center], width*2.0*pC[0].length, center1, width1, rho, rho_stat_err, start, stop);)
	
	  if(flag_mode==0) {
	    fprintf(file_out_excl,"%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\n", pA->center[count_center], width_real*2.0*pC[0].length, real_center1, real_width1,  rho_real, rho_stat_err_real, real_start, real_stop); }
	
	jack_values[count_center+jack*pA->N_center]=rho_real;
	
	real_center1_array[count_center]=real_center1;
      
	LOG_FILE_OPERATION(fclose(file_out);)

      }
      if(flag_mode==0) {
	fclose(file_out_excl); }
    }
    gsl_vector_free(Q);
    gsl_vector_free(Q_real);
    gsl_vector_free(Q_initial);

    //now compute average rho and it's jack knife error
    double avg, err;
    for(count_center=0; count_center<pA->N_center; count_center++) {
      avg=err=0.0;
      for(jack=0;jack<num_jack_samples;jack++) {
	avg+=jack_values[count_center+jack*pA->N_center];
      }
      avg=avg/(double)num_jack_samples;
      for(jack=0;jack<num_jack_samples;jack++) {
	err+=(jack_values[count_center+jack*pA->N_center]-avg)*(jack_values[count_center+jack*pA->N_center]-avg);
      }
      err=sqrt(err/((double)num_jack_samples*(num_jack_samples-1.0)));
      if(flag_mode==1) {
	Rho_avg->rho_array[lambda_count][count_center]=avg;
	Rho_avg->rho_err_array[lambda_count][count_center]=err;
      }
      fprintf(file_out_avg,"%.15le\t%.15le\t%.15le\t%.15le\n", pA->center[count_center], avg, err, real_center1_array[count_center] );
    }

    free(jack_values);
    free(real_center1_array);
    fclose(file_out_avg);
}
