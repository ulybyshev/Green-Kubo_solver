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
    lambda =1.0-accuracy/gsl_matrix_get(pC->S,0,0);
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
    lambda =accuracy/gsl_matrix_get(pC->S,0,0);
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
	    for(omega=0;omega<omega_plot_limit/(2.0*pC->length);omega+=omega_plot_delta/(2.0*pC->length))
	    {
    		fprintf(file_out,"%.15le\t%.15le\t%.15le\n", omega*2.0*pC->length, delta(omega,Q, pC), delta(omega, Q_real, pC));
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
	fprintf(file_out_excl,"%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\n", pA->center[count_center], width_real*2.0*pC->length, real_center1, real_width1,  rho_real, rho_stat_err_real, real_start, real_stop);


	LOG_FILE_OPERATION(fclose(file_out);)


	gsl_vector_free(Q);
	gsl_vector_free(Q_real);

    }

gsl_vector_free(Q_initial);
  
 
}



