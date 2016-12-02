#include "math_functions.h"
#include "input_output.h"
#include "kernel.h"


//for lattice version of the kernel functions
double lattice_exp(double omega, int euclidean_time)
{//we assume that euclidean time is >=0 
    if(euclidean_time>=0)
    {
	if (omega>-1.0)
	    return my_pow(1.0+omega, euclidean_time);
	else
	    return 0.0;
    }
    else
    {
	return 0.0;
    }
}
//for lattice exponent
double my_pow(double x, int b)
{
    int i=0;
    double res=1.0;
    for(i=0;i<b;i++)
	res*=x;
    return res;
}




//simple kernel for time points (works in cases of flag_model==1)
double kernel_points(int i, correlator* pC,double omega)
{
    int real_t;
    real_t=pC->points_numbers[i-1];
    if(omega*pC->length<accuracy)
	return kernel_function(omega+accuracy/pC->length, 2*pC->N_full_points, real_t);
//	1.0/(PI*pC->length);
    else
	return kernel_function(omega, 2*pC->N_full_points, real_t);
//	 (omega/PI)*(exp(omega*(-(double)real_t))+exp(omega*((double)real_t-2.0*pC->length)))/(1.0-exp(-2.0*omega*pC->length));
}
//kernel at zero omega (to exclude zeros of denominator from the calculation)
double kernel0_points(int i, correlator* pC)
{
  return kernel_function(accuracy/pC->length, 2*pC->N_full_points, i);
}

//kernel in case of averaging over intervals (works in case of flag_model==2)
double kernel_intervals(int i, correlator* pC, double omega)
{
    int interval_length=pC->interval_numbers[i-1]->size;  
    int real_t, j;
    double result=0.0;
    for(j=0;j<interval_length;j++) 
    {
      real_t=pC->interval_numbers[i-1]->times[j];
      if(omega*pC->length<accuracy) 
      {
	result+=kernel_function(omega+accuracy/pC->length, 2*pC->N_full_points, real_t);
      }
      else 
      {
	result+=kernel_function(omega, 2*pC->N_full_points, real_t);;
      }
    }
    return result/((double)interval_length);
}
//kernel at zero omega (to exclude zeros of denominator from the calculation)
double kernel0_intervals(int i, correlator* pC)
{
  return kernel_function(accuracy/pC->length, 2*pC->N_full_points, i);
}


//general wrapper
double kernel(int i, correlator* pC, double omega)
{
    if(flag_model==1 || flag_model==0)
	return kernel_points(i, pC, omega);
    else if(flag_model==2)
	return kernel_intervals(i, pC, omega);
    else
	return 0;
}
double kernel0(int i, correlator* pC)
{
    if(flag_model==1 || flag_model==0)
	return kernel0_points(i, pC);
    else if(flag_model==2)
	return kernel0_intervals(i, pC);
    else
	return 0;
}





double W_function (w_parameters params, double omega)
{
  return (omega-params.omega_center)*(omega-params.omega_center)*kernel(params.i,params.pC,omega)*kernel(params.j,params.pC, omega);
}

double delta(double omega, gsl_vector* Q, correlator* pC)//omega_center is defined in calculation of vector Q
{
  int i;
  double result=0.0;
  for(i=0;i<Q->size;i++)
  {
    result+=gsl_vector_get(Q,i)*kernel(i+1,pC, omega);
  }
  return result;
}
//at zero omega (to exclude zeros in denominator)
double delta0(gsl_vector* Q, correlator* pC)//omega_center is defined in calculation of vector Q
{
  int i;
  double result=0.0;
  for(i=0;i<Q->size;i++)
  {
    result+=gsl_vector_get(Q,i)*kernel0(i+1, pC);
  }
  return result;
}



double delta_sq_int (double x, void* params){
  d_parameters  Q_and_center=* (d_parameters*) params;
  double res_delta=delta(x, Q_and_center.Q, Q_and_center.pC);
  double result=res_delta*res_delta*(x-Q_and_center.center)*(x-Q_and_center.center);
  return result;
}


//functions for integration
//simple kernel integration (integer parameter t should NOT be equal to 0)
double kernel_int (double x, void * params) {
  kernel_parameters buffer=* (kernel_parameters*) params;
  return kernel(buffer.t,buffer.pC,x);
}


//omega*kernel integration (integer parameter t should NOT be equal to 0)
//we need  it to define real center of the delta-function
double omega_kernel_int (double x, void * params) {
  kernel_parameters buffer=* (kernel_parameters*) params;
  return  x*kernel(buffer.t,buffer.pC,x);
}



//simple kernel integration (integer parameter t should NOT be equal to 0)
double W_int (double x, void * params) {
  w_parameters buffer = *(w_parameters *) params;
  return  W_function(buffer,x);
}




//integration to obtain R vector
void R_integration(calc_structures* pA,correlator* pC)
{
    double int_result=2.0, int_error=1.0;

    int t;
    FILE* file_out;
    LOG_FILE_OPERATION(file_out=fopen_control("R_control.txt","w");)
    
    {
	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (N_int_steps);
	gsl_function F;
	F.function = &kernel_int;
	kernel_parameters buffer;
	buffer.pC=pC;
	for(t=1;t<=pC->N_valid_points;t++)
	{
	    buffer.t=t;
	    F.params = &buffer;
	    if(kernel_switcher!=2)
	    	    gsl_integration_qagiu (&F, 0.0, 0.0, accuracy, N_int_steps, w, &int_result, &int_error); 
	    else
	    {
	    	    gsl_integration_qag (&F, 0.0, 1.0, 0.0, accuracy, N_int_steps,GSL_INTEG_GAUSS41, w, &int_result, &int_error); 
	    }
	    gsl_vector_set(pA->R,t-1,int_result);
	    
	    LOG_FILE_OPERATION(fprintf(file_out,"%d\t%.15le\t%.15le\n", t, int_result, int_error); fflush(file_out);)
	}
	gsl_integration_workspace_free (w);
    }
    LOG_FILE_OPERATION(fclose(file_out);)
}
///////

//integration to obtain omega_R vector
void omega_R_integration(calc_structures* pA,correlator* pC)
{
    double int_result, int_error;

    int t;
    FILE* file_out;
    LOG_FILE_OPERATION(file_out=fopen_control("omega_R_control.txt","w");)
  
    {
	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (N_int_steps);
	gsl_function F;
	F.function = &omega_kernel_int;
	kernel_parameters buffer;
	buffer.pC=pC;
	for(t=1;t<=pC->N_valid_points;t++)
	{
	    buffer.t=t;
	    F.params = &buffer;
	    if(kernel_switcher!=2)
	    	gsl_integration_qagiu (&F, 0.0, 0.0, accuracy, N_int_steps, w, &int_result, &int_error); 
	    else
		gsl_integration_qag (&F, 0.0, 1.0, 0.0, accuracy, N_int_steps,GSL_INTEG_GAUSS41, w, &int_result, &int_error); 
    	    gsl_vector_set(pA->omega_R,t-1,int_result);
	    LOG_FILE_OPERATION(fprintf(file_out,"%d\t%.15le\t%.15le\n", t, int_result, int_error); fflush(file_out);)
	}
	gsl_integration_workspace_free (w);

    }
    LOG_FILE_OPERATION(fclose(file_out);)
}
///////

void W_integration(gsl_matrix * W, correlator* pC, double center)
{
    int i,j;
     double int_result, int_error;
    FILE* file_out;
    //Integration
  LOG_FILE_OPERATION(file_out=fopen_log("W_calc_control.txt","a", center);)
  gsl_integration_workspace * w1  = gsl_integration_workspace_alloc (N_int_steps);
  gsl_function F;
  F.function = &W_int;
  w_parameters buffer;
  buffer.omega_center=center;  
  buffer.pC=pC;
  for(i=1;i<=pC->N_valid_points;i++)
    for(j=1;j<=pC->N_valid_points;j++)
    {
      buffer.i=i;
      buffer.j=j;
      F.params = &buffer;
      if(kernel_switcher!=2)
            gsl_integration_qagiu (&F, 0.0, 0.0, accuracy, N_int_steps, w1, &int_result, &int_error);
	else
	    gsl_integration_qag (&F, 0.0, 1.0, 0.0, accuracy, N_int_steps,GSL_INTEG_GAUSS41, w1, &int_result, &int_error); 
      gsl_matrix_set(W,i-1,j-1,int_result);
      LOG_FILE_OPERATION(fprintf(file_out,"%d\t%d\t%.15le\t%.15le\n", i,j, int_result, int_error); fflush(file_out);)
    }
  gsl_integration_workspace_free (w1);
  LOG_FILE_OPERATION(fclose(file_out);)
}


void calculate_Q(gsl_vector* Q, calc_structures* pA, correlator* pC, int center_count)
{
  int i,j,k;
  double center;
  FILE* file_out;
  double par_double;
  double int_result, int_error;

  gsl_matrix * W, *W_1;  //W is the full W matrix which takes into account regularization
  center=pA->center[center_count]/(2.0*pC->length);

//full W matrix(with regularization) and inverse W matrix
 W=gsl_matrix_alloc(pC->N_valid_points, pC->N_valid_points);
 W_1=gsl_matrix_alloc(pC->N_valid_points, pC->N_valid_points);

{
 
  //creating local copy of W matrix
  for(i=1;i<=pC->N_valid_points;i++)
    for(j=1;j<=pC->N_valid_points;j++)
    {
        gsl_matrix_set(W,i-1,j-1,gsl_matrix_get(pA->W[center_count], i-1, j-1));
    }
  
  //Regularization of the 1st type
  if (flag_lambda_regularization==1)
  {
    for(i=1;i<=pC->N_valid_points;i++)
    for(j=1;j<=pC->N_valid_points;j++)
    {
        par_double=lambda*gsl_matrix_get(W, i-1, j-1)+(1.0-lambda) * gsl_matrix_get(pC->S, i-1, j-1);
        gsl_matrix_set(W,i-1,j-1,par_double);
    }
  }
  
  //Save
  LOG_FILE_OPERATION(file_out=fopen_log("W_matrix.txt","a", center);)
  for(i=1;i<=pC->N_valid_points;i++){
    for(j=1;j<=pC->N_valid_points;j++)
      LOG_FILE_OPERATION(fprintf(file_out,"%.5le\t", gsl_matrix_get(W,i-1,j-1)); fflush(file_out); )
    LOG_FILE_OPERATION(fprintf(file_out,"\n"); fflush(file_out); )
  }
  LOG_FILE_OPERATION(fclose(file_out);)
}


//CALCULATION Q = W^(-1) pA->R
//SVD decomposition is used  
if (flag_lambda_regularization==2)
{
  
  gsl_matrix * Uk, *Vk;
  gsl_vector * workT;
  gsl_vector* Wk;//W eigenvalues

  Uk=gsl_matrix_calloc(pC->N_valid_points,pC->N_valid_points);
  Vk=gsl_matrix_calloc(pC->N_valid_points,pC->N_valid_points);
  Wk=gsl_vector_calloc(pC->N_valid_points);
  workT=gsl_vector_calloc(pC->N_valid_points);

  for(i=0;i<pC->N_valid_points;i++)
    for(j=0;j<pC->N_valid_points;j++)
      gsl_matrix_set(Uk, j,i, gsl_matrix_get(W, i,j));

  gsl_linalg_SV_decomp(Uk,Vk,Wk,workT);
  
  //Different regularization type
  double sv0 = gsl_vector_get(Wk,0);
  for(j=1;j<pC->N_valid_points;j++){
    par_double = gsl_vector_get(Wk,j);
    if(par_double < sv0*lambda)
      par_double = 0;
    gsl_vector_set(Wk,j,par_double);
  }
  
  gsl_linalg_SV_solve(Uk,Vk,Wk,pA->R,Q);

  //Check W^(-1)
  double sum=0.0;
  for(i=0;i<pC->N_valid_points;i++){
    par_double = 0.0;
    for(j=0;j<pC->N_valid_points;j++)
      par_double += gsl_matrix_get(W,i,j)*gsl_vector_get(Q, j);
    sum+= fabs(par_double / gsl_vector_get(pA->R, i) - 1.0);
  }    
  LOG_FILE_OPERATION(file_out=fopen_log("Discrepancy.txt","a", center);)
  LOG_FILE_OPERATION(fprintf(file_out,"Discrepancy:%.15le\n",sum);fflush(file_out);)
  LOG_FILE_OPERATION(fclose(file_out);)
    
  LOG_FILE_OPERATION(file_out=fopen_log("Q_vector.txt","a", center);)
  //Normalize Q = Q/(pA->R*Q)
  double norma=0.0;
  for(i=0;i<pC->N_valid_points;i++)
    norma +=  gsl_vector_get(Q, i)*gsl_vector_get(pA->R, i);
  LOG_FILE_OPERATION(fprintf(file_out, "norma_old=%.15le\n", norma);fflush(file_out);)
  
  for(i=0;i<pC->N_valid_points;i++){
    par_double = gsl_vector_get(Q, i);
    gsl_vector_set(Q,i,par_double/norma);
  }
  
  par_double = 0.0;
  for(i=0;i<pC->N_valid_points;i++)
    par_double+=gsl_vector_get(Q, i)*gsl_vector_get(pA->R, i);
  LOG_FILE_OPERATION(fprintf(file_out, "norma_new=%.15le\n", par_double);fflush(file_out);)   
  
  //Save Q vector
  for(i=0;i<pC->N_valid_points;i++)
    LOG_FILE_OPERATION(fprintf(file_out,"%.15le\n",gsl_vector_get(Q, i));fflush(file_out);)
    
  LOG_FILE_OPERATION(fclose(file_out); )
  gsl_vector_free(Wk);
  gsl_matrix_free(Uk);
  gsl_matrix_free(Vk);    
  gsl_vector_free(workT); 
  }
else  
{
 //calculation of inverse W matrix (stored in W_1) 
//SVD decomposition is used
{
  gsl_matrix * Uk, *Vk;
  gsl_vector * workT;
  gsl_vector* Wk;//W eigenvalues

    Wk=gsl_vector_calloc(pC->N_valid_points);
    Uk=gsl_matrix_calloc(pC->N_valid_points,pC->N_valid_points);
    Vk=gsl_matrix_calloc(pC->N_valid_points,pC->N_valid_points);
    workT=gsl_vector_calloc(pC->N_valid_points);
    
    for(i=0;i<pC->N_valid_points;i++)
    for(j=0;j<pC->N_valid_points;j++)
    {
      par_double=gsl_matrix_get(W, i,j);
      gsl_matrix_set(Uk, j,i, par_double);
    }
    
    gsl_linalg_SV_decomp(Uk,Vk,Wk,workT);
    
    LOG_FILE_OPERATION(file_out=fopen_log("W_1_matrix.txt","a", center);)//output of inverse kernel matrix
    for(i=0;i<pC->N_valid_points;i++)
    {
      for(k=0;k<pC->N_valid_points;k++)
      {
        par_double=0.0;
        for(j=0;j<pC->N_valid_points;j++)
        {
          par_double+=gsl_matrix_get(Vk,k,j)*gsl_matrix_get(Uk,i,j)*(1.0/gsl_vector_get(Wk,j));
        }
        gsl_matrix_set(W_1, i,k, par_double);
        LOG_FILE_OPERATION(fprintf(file_out,"%.5le\t",gsl_matrix_get(W_1, i,k));fflush(file_out);)
      }
      LOG_FILE_OPERATION(fprintf(file_out,"\n");fflush(file_out);)
    }
    LOG_FILE_OPERATION(fclose(file_out);)
    

    LOG_FILE_OPERATION(file_out=fopen_log("Wk_vector.txt","a", center);)//output of eigenvalues of kernel matrix
    for(i=0;i<pC->N_valid_points;i++)
    {
      LOG_FILE_OPERATION(fprintf(file_out,"%.15le\n",gsl_vector_get(Wk, i));fflush(file_out);)
    }
    LOG_FILE_OPERATION(fclose(file_out);)

    LOG_FILE_OPERATION(file_out=fopen_log( "W_inv_control.txt","a", center);)//for control - output of K^(-1)*K
    for(i=0;i<pC->N_valid_points;i++)
    {
      for(k=0;k<pC->N_valid_points;k++)
      {
        par_double=0.0;
        for(j=0;j<pC->N_valid_points;j++)
        {
          par_double+=gsl_matrix_get(W,i,j)*gsl_matrix_get(W_1,j,k);
        }
      LOG_FILE_OPERATION(fprintf(file_out,"%.5le\t",par_double);fflush(file_out);)
      }
      LOG_FILE_OPERATION(fprintf(file_out,"\n");fflush(file_out);)
    }
    LOG_FILE_OPERATION(fclose(file_out);)
 }//end of W_1 calculation 
 
 //vector Q calculation
  {
    double norma;
  
  par_double=0.0;
  for(i=0;i<pC->N_valid_points;i++)
  {
    for(j=0;j<pC->N_valid_points;j++)
    {
      par_double+=gsl_vector_get(pA->R, i) * gsl_vector_get(pA->R, j) * gsl_matrix_get(W_1,i,j);
    }
  }
  norma=par_double;
  for(i=0;i<pC->N_valid_points;i++)
  {
    par_double=0.0;
    for(j=0;j<pC->N_valid_points;j++)
    {
      par_double+= gsl_vector_get(pA->R, j) * gsl_matrix_get(W_1,i,j);
    }
    gsl_vector_set(Q,i,par_double/norma);
  }
 
  }

//vector Q output
  LOG_FILE_OPERATION(file_out=fopen_log("Q_vector.txt","a", center);)

//normalization
  par_double=0.0;
  for(i=0;i<pC->N_valid_points;i++)
  {
    par_double += gsl_vector_get(Q, i) * gsl_vector_get(pA->R,i);
  }
  LOG_FILE_OPERATION(fprintf(file_out, "norma_old=%.15le\n", par_double);fflush(file_out);)
  
  for(i=0;i<pC->N_valid_points;i++)
  {
    gsl_vector_set(Q,i,gsl_vector_get(Q, i)/par_double);
  }
  for(i=0;i<pC->N_valid_points;i++)
  {
     LOG_FILE_OPERATION(fprintf(file_out,"%.15le\n",gsl_vector_get(Q, i));fflush(file_out);)
  }
 ////
  par_double=0.0;
  for(i=0;i<pC->N_valid_points;i++)
  {
    par_double += gsl_vector_get(Q, i) * gsl_vector_get(pA->R,i);
  }
  LOG_FILE_OPERATION(fprintf(file_out, "norma_new=%.15le\n", par_double);fflush(file_out);)

  LOG_FILE_OPERATION(fclose(file_out);)



}
  gsl_matrix_free(W_1);
  gsl_matrix_free(W);
 

}


void calculate_rho(gsl_vector* Q, correlator* pC, double* rho, double* rho_stat_err)
{
    int i;
    double res, res_error;
    res=0.0;
    res_error=0.0;
    for(i=0;i<pC->N_valid_points;i++)
    {
	res+=pC->corr[i]*gsl_vector_get(Q, i);
	res_error+=pC->error[i]*pC->error[i]*gsl_vector_get(Q, i)*gsl_vector_get(Q, i);
    }
    res_error=sqrt(res_error);
    *rho=res;
    *rho_stat_err=res_error;
}

double delta_width_calculation(gsl_vector* Q, double center, correlator* pC)
{
  double int_result, int_error;
  gsl_function F;
  gsl_integration_workspace * w1  = gsl_integration_workspace_alloc (N_int_steps);
  FILE* file_out;

  LOG_FILE_OPERATION(file_out=fopen_log("delta_width_control.txt","a", center);)

  F.function = &delta_sq_int;
  d_parameters buffer;
  buffer.Q=Q;
  buffer.center=center;  
  buffer.pC=pC;
    F.params = &buffer;
    if(kernel_switcher!=2)
	gsl_integration_qagiu (&F, 0.0, 0.0, accuracy, N_int_steps, w1, &int_result, &int_error);
    else
	gsl_integration_qag (&F, 0.0, 1.0, 0.0, accuracy, N_int_steps,GSL_INTEG_GAUSS41, w1, &int_result, &int_error); 
    LOG_FILE_OPERATION(fprintf(file_out,"%.15le\t%.15le\t%.15le\n",center, int_result, int_error); fflush(file_out);)
   gsl_integration_workspace_free (w1);
  LOG_FILE_OPERATION(fclose(file_out);)
  return sqrt(int_result);
}


//calculation of real chracteristics of delta function
//returns all values in temperature units
void delta_characteristics_calculation(double* start, double* stop, double* center_real, gsl_vector* Q, double center, correlator* pC)
{
  //delta function output
  FILE* file_out;
  double omega,omega_max, function_max;
  unsigned long int count=0;

  LOG_FILE_OPERATION(file_out=fopen_log("real_delta_width_control.txt","a", center);)
  for(omega=omega_plot_delta/(2.0*pC->length);omega<omega_plot_limit/(2.0*pC->length);omega+=omega_plot_delta/(2.0*pC->length))
  {
    if(count==0)
      {
        function_max= delta(omega,Q, pC);
        omega_max=omega;
      }
    else
      {
        if(function_max<delta(omega, Q, pC))
        {
          function_max=delta(omega,Q, pC);
          omega_max=omega;
        }
      }
    count++;  
  }
  LOG_FILE_OPERATION(fprintf(file_out,"omega_max=%.15le\n",omega_max*(2.0*pC->length)); fflush(file_out);)
  (*start)=0.0;
 


  //start calculation
  for(omega=omega_max;omega>omega_plot_delta/(2.0*pC->length);omega-=omega_plot_delta/(2.0*pC->length))
  {
    if(delta(omega,Q, pC)<function_max/2.0)
     {
       (*start)=omega;
       break;
     }
  }
  LOG_FILE_OPERATION(fprintf(file_out,"omega_start=%.15le\n",(*start)*(2.0*pC->length)); fflush(file_out);)

//stop calculation
  for(omega=omega_max;omega<omega_plot_limit/(2.0*pC->length);omega+=omega_plot_delta/(2.0*pC->length))
  {
    if(delta(omega,Q, pC)<function_max/2.0)
     {
       (*stop)=omega;
       break;
     }
  }
  LOG_FILE_OPERATION(fprintf(file_out,"omega_stop=%.15le\n",(*stop)*(2.0*pC->length)); fflush(file_out);)
  (*center_real)=omega_max;;//(  (*start)   +   (*stop))/2.0;
  
  (*start)=(*start) * 2.0*pC->length;
  (*stop)=(*stop) * 2.0*pC->length;
  (*center_real)=(*center_real) * 2.0*pC->length;

}




