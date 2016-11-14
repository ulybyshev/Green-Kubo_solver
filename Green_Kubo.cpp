//to run the program:   ./program  Nt  file_current.txt  file_matrix.txt path_for_output  file_with_parameters.txt
// file for currents:
// 1    current_correlator    current_correlator_error 
// 2 ....
// 3 ....
// 4 ....
// Nt ....
//N_points = HALF of physical timeslices  in real simulations,
// current-current correlator should be symmetrized before input
// value of correlator at zero euclidean time should be excluded from input
//errors  from current-current correlator file are taken into account in calculation of statistical errors of spectral function 
//covariance matrix should be calculated for symmetrized correlator (size of the matrix thus is also reduced)
//input for matrix is formatted as
//C_{11} .......  C_{1 Nt}
//..........................
//C_{Nt 1} .......  C_{Nt Nt}
//(without any indexes)


//file with parameters (comments lines alternate with parameter values). All dimensional quantities are in units of temperature
//double accuracy=1e-8;
//
//long int N_int_steps=1000000;
//
//double omega_plot_delta=0.01;
//
//double omega_plot_limit=30.0;
//
//double lambda=1.0;
//
//double center_start=3.0;
//
//double center_stop=18.01;
//
//double center_delta=3.0;
//
//int flag_exclude_delta_function (0 if not exclude - then the next parameter is not taken into account, 1 - exclude starting from the count_start iteration in center)
//
//int count_start - iteration in change the center to start exclude initial delta function (interations start from 0)
//
//int flag model
//      //(0 - all timslices are taken into account, 1 - partially, then the list of valid timeslices starts immediately after in the format:
//N_{valid_timeslices}
//i1
//i2
//.....
//i_{N_{valid}}  
//(sorted ascending, contact term in correlator is under number 0, so valid timeslices start from 1 and go till Nt)

//output - all dimensional quantities in units of temperature

//output - delta functions
//omega(in units of temperature)    function   function_after_exclude

//output rho
// center(as in was in input)  resolution_width(as it was in D calculation)  center(Maximum)  width((start-stop)/2, start and stop are at half-peak)  rho  rho_stat_err

//output rho_exclude - the same for spectral functions with exclusion at zero 
// center(as in was in input)  resolution_width(as it was in D calculation)  center(Maximum)  width((start-stop)/2, start and stop are at half-peak)  rho  rho_stat_err  start stop

//TODO
//rewrite the output of resolution function width for non-symmetrical error (center, low, high), including separate processing of the case when low=zero


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>

#define PI 3.1415926

#define SKIP_COMMENT_LINE(file) { while(getc(file)!='\n'); }
#define SKIP_REMAINING_CHARS(file) SKIP_COMMENT_LINE(file)
#define SKIP_OPTION(file) { SKIP_COMMENT_LINE(file) SKIP_REMAINING_CHARS(file) }

void parse_option(FILE* file, const char* pattern, void* ptr)
{
    SKIP_COMMENT_LINE(file)
     fscanf(file, pattern, ptr);
    SKIP_REMAINING_CHARS(file)
}

struct interval
{  
  int size;
  int *times;
};

int Nt_2;//actually  it is equal to HALF of real timeslices (current-current correlator should be symmetrized before further usage)
double dNt_2;

double accuracy=1e-8;

long int N_int_steps=1000000;

double omega_plot_delta=0.01;
double omega_plot_limit=30.0;

double lambda=1.0;

//in units of temperature
double center_start=3.0;
double center_stop=18.01;
double center_delta=3.0;

char directory[1000];

int flag_model;
int N_valid_points;
int* points_numbers;
int N_intervals=0;
interval *interval_numbers;
int Nt_2_pre;
double dNt_2_pre;

int flag_exclude_delta;
int count_start_exclude;


//in all parameters structures  we have the number of timeslice in array, NOT real time
struct w_parameters
{
  double omega_center;
  int i;
  int j;
};

struct d_parameters
{
  double center;
  gsl_vector* Q;
};
//new idea - regularization by the free conductivity (including delta-function - but with  zero-omega delta function there 
//could be a problem), then it will be possible to work with fourier-transformed kernel
// in case of this kernel the cut-off of high frequencies (\nu) looks more natural

//real width = width at half heigth of global maximum 
//we calculate start and stop and
//real center = middle between start and stop

//fourier transformed kernel
/*
double kernel (int i, double omega)
{
    return  1.0/((double)(i*i)+omega*omega);
}
*/

//in this kernel egularization at low frequncies is made by tanh
/*
double kernel (int i, double omega)
{
    return  (exp(omega*(-(double)i))+exp(omega*((double)i-2.0*dNt_2_pre)))/(1.0+exp(-2.0*omega*dNt_2_pre));
}*/

void construct_intervals(int *points, int num_points, interval *intervals) {

  int i,j,k,count_interval,interval_length;

  count_interval=num_points-1;
  for(i=(num_points-1);i>=1;i--) {
    if((points[i]-1)==points[i-1]) { //single "red" point
      intervals[count_interval].size=1;
      intervals[count_interval].times=(int *)calloc(1,sizeof(int));
      intervals[count_interval].times[0]=points[i];
    }
    else { //contains "black points"
      interval_length=(points[i]-points[i-1]);
      intervals[count_interval].size=interval_length;
      intervals[count_interval].times=(int *)calloc(interval_length,sizeof(int));
      for(j=points[i],k=(interval_length-1);j>points[i-1];j--,k--) {
	intervals[count_interval].times[k]=j;
      }
    }
    count_interval--;
  }
  intervals[0].size=1;
  intervals[0].times=(int *)calloc(1,sizeof(int));
  intervals[0].times[0]=points[0];
  
}

//regularization is made by omega 
    //i is the number of timeslice in array+1, NOT real time
double kernel(int i, double omega)
{
    int interval_length=interval_numbers[i-1].size;   
    int real_t, j;
    double result=0.0;
    for(j=0;j<interval_length;j++) {
      //real_t=points_numbers[i-1];
      real_t=interval_numbers[i-1].times[j];
      if(omega*dNt_2_pre<accuracy) {
	//return cosh(omega*(dNt_2-(double)real_t))/(PI*dNt_2_pre);
	result+=cosh(omega*(dNt_2-(double)real_t))/(PI*dNt_2_pre);
      }
      else {
	//return  (omega/PI)*(exp(omega*(-(double)real_t))+exp(omega*((double)real_t-2.0*dNt_2_pre)))/(1.0-exp(-2.0*omega*dNt_2_pre));
	result+=(omega/PI)*(exp(omega*(-(double)real_t))+exp(omega*((double)real_t-2.0*dNt_2_pre)))/(1.0-exp(-2.0*omega*dNt_2_pre));
      }
    }
    return result/((double)interval_length);
}
//kernel at zero omega (to exclude zeros of denominator from the calculation)
double kernel0(int i)
{
  int interval_length=interval_numbers[i-1].size;
  return 1.0/(((double)interval_length)*PI*dNt_2_pre);
}


double W_function (w_parameters params, double omega)
{
  return (omega-params.omega_center)*(omega-params.omega_center)*kernel(params.i,omega)*kernel(params.j, omega);
}

double delta(double omega, gsl_vector* Q)//omega_center is defined in calculation of vector Q
{
  int i;
  double result=0.0;
  for(i=0;i<Q->size;i++)
  {
    result+=gsl_vector_get(Q,i)*kernel(i+1, omega);
  }
  return result;
}
//at zero omega (to exclude zeros in denominator)
double delta0(gsl_vector* Q)//omega_center is defined in calculation of vector Q
{
  int i;
  double result=0.0;
  for(i=0;i<Q->size;i++)
  {
    result+=gsl_vector_get(Q,i)*kernel0(i+1);
  }
  return result;
}



double delta_sq_int (double x, void* params){
  d_parameters  Q_and_center=* (d_parameters*) params;
  double res_delta=delta(x, Q_and_center.Q);
  double result=res_delta*res_delta*(x-Q_and_center.center)*(x-Q_and_center.center);
  return result;
}

//functions for integration
//simple kernel integration (integer parameter t should NOT be equal to 0)
double kernel_int (double x, void * params) {
  int t = *(int *) params;
  double result = kernel(t,x);
  return result;
}



//omega*kernel integration (integer parameter t should NOT be equal to 0)
//we need  it to define real center of the delta-function
double omega_kernel_int (double x, void * params) {
  int t = *(int *) params;
  double result = x*kernel(t,x);
  return result;
}



//simple kernel integration (integer parameter t should NOT be equal to 0)
double W_int (double x, void * params) {
  w_parameters buffer = *(w_parameters *) params;
  double result = W_function(buffer,x);
  return result;
}

FILE* fopen_log(const char* name, const char* aim, double parameter)
{
  FILE* file_out;
  char final_name[1000];
  sprintf(final_name, "%s/%s", directory, name);
  file_out=fopen(final_name, aim);
  fprintf(file_out,"\n\n\n******************************\n");
  fprintf(file_out,"center/T=%.15le\n", parameter*2.0*dNt_2_pre);
  fprintf(file_out,"******************************\n");

  fflush(file_out);
  return file_out;
}

FILE* fopen_control(const char* name, const char* aim)
{
  FILE* file_out;
  char final_name[1000];
  sprintf(final_name, "%s/%s", directory, name);
  file_out=fopen(final_name, aim);
  
  return file_out;
}

void calculate_Q(gsl_vector* Q, gsl_vector* R, gsl_matrix* S, double center)
{
  int i,j,k,t;
  
  FILE* file_out;
  int par_int;
  double par_double;
  double int_result, int_error;
  double omega;

  gsl_matrix * W;

  gsl_matrix* W_1;
  gsl_vector* Wk;//W eigenvalues

//integration to obtain W matrix
 W=gsl_matrix_alloc(Nt_2, Nt_2);

 file_out=fopen_log("W_calc_control.txt","a", center);
 
{
  gsl_integration_workspace * w1  = gsl_integration_workspace_alloc (N_int_steps);
 
  gsl_function F;
  F.function = &W_int;
  w_parameters buffer;
  buffer.omega_center=center;  
  for(i=1;i<=Nt_2;i++)
  for(j=1;j<=Nt_2;j++)
  {
    buffer.i=i;
    buffer.j=j;
    F.params = &buffer;
    gsl_integration_qagiu (&F, 0.0, accuracy, accuracy, N_int_steps, w1, &int_result, &int_error);
    gsl_matrix_set(W,i-1,j-1,int_result);
    fprintf(file_out,"%d\t%d\t%.15le\t%.15le\n", i,j, int_result, int_error); fflush(file_out);
  }
  gsl_integration_workspace_free (w1);
  
  //regularization
  for(i=1;i<=Nt_2;i++)
  for(j=1;j<=Nt_2;j++)
  {
      if(i==j)
      {
        par_double=gsl_matrix_get(W, i-1, j-1);
        par_double=(1.0-lambda) * gsl_matrix_get(S, i-1, j-1)+ lambda*par_double;
        gsl_matrix_set(W,i-1,j-1,par_double);
      } 
  }



}
  fclose(file_out);

  file_out=fopen_log("W_matrix.txt","a", center);
  for(i=1;i<=Nt_2;i++)
  {
    for(j=1;j<=Nt_2;j++)
    {
     fprintf(file_out,"%.5le\t", gsl_matrix_get(W,i-1,j-1)); fflush(file_out); 
    }
     fprintf(file_out,"\n"); fflush(file_out); 
  }
  fclose(file_out);
//calculation of inverse W matrix (stored in W_1) 
//SVD decomposition is used
  
    W_1=gsl_matrix_calloc(Nt_2,Nt_2);
    Wk=gsl_vector_calloc(Nt_2);
{
    gsl_matrix * Uk, *Vk;
    gsl_vector * workT;

    Uk=gsl_matrix_calloc(Nt_2,Nt_2);
    Vk=gsl_matrix_calloc(Nt_2,Nt_2);
    workT=gsl_vector_calloc(Nt_2);
    
    for(i=0;i<Nt_2;i++)
    for(j=0;j<Nt_2;j++)
    {
      par_double=gsl_matrix_get(W, i,j);
      gsl_matrix_set(Uk, j,i, par_double);
    }
    
    gsl_linalg_SV_decomp(Uk,Vk,Wk,workT);
    
    file_out=fopen_log("W_1_matrix.txt","a", center);//output of inverse kernel matrix
    for(i=0;i<Nt_2;i++)
    {
      for(k=0;k<Nt_2;k++)
      {
        par_double=0.0;
        for(j=0;j<Nt_2;j++)
        {
          par_double+=gsl_matrix_get(Vk,k,j)*gsl_matrix_get(Uk,i,j)*(1.0/gsl_vector_get(Wk,j));
        }
        gsl_matrix_set(W_1, i,k, par_double);
        fprintf(file_out,"%.5le\t",gsl_matrix_get(W_1, i,k));fflush(file_out);
      }
      fprintf(file_out,"\n");fflush(file_out);
    }
    fclose(file_out);
    

    file_out=fopen_log("Wk_vector.txt","a", center);//output of eigenvalues of kernel matrix
    for(i=0;i<Nt_2;i++)
    {
      fprintf(file_out,"%.15le\n",gsl_vector_get(Wk, i));fflush(file_out);
    }
    fclose(file_out);

    file_out=fopen_log( "W_inv_control.txt","a", center);//for control - output of K^(-1)*K
    for(i=0;i<Nt_2;i++)
    {
      for(k=0;k<Nt_2;k++)
      {
        par_double=0.0;
        for(j=0;j<Nt_2;j++)
        {
          par_double+=gsl_matrix_get(W,i,j)*gsl_matrix_get(W_1,j,k);
        }
      fprintf(file_out,"%.5le\t",par_double);fflush(file_out);
      }
      fprintf(file_out,"\n");fflush(file_out);
    }
    fclose(file_out);
 }//end of W_1 calculation 

 //vector Q calculation
  {
    double norma;
  
  par_double=0.0;
  for(i=0;i<Nt_2;i++)
  {
    for(j=0;j<Nt_2;j++)
    {
      par_double+=gsl_vector_get(R, i) * gsl_vector_get(R, j) * gsl_matrix_get(W_1,i,j);
    }
  }
  norma=par_double;
  for(i=0;i<Nt_2;i++)
  {
    par_double=0.0;
    for(j=0;j<Nt_2;j++)
    {
      par_double+= gsl_vector_get(R, j) * gsl_matrix_get(W_1,i,j);
    }
    gsl_vector_set(Q,i,par_double/norma);
  }
 
  }

//vector Q output
  file_out=fopen_log("Q_vector.txt","a", center);

//normalization
  par_double=0.0;
  for(i=0;i<Nt_2;i++)
  {
    par_double += gsl_vector_get(Q, i) * gsl_vector_get(R,i);
  }
  fprintf(file_out, "norma_old=%.15le\n", par_double);fflush(file_out);
  
  for(i=0;i<Nt_2;i++)
  {
    gsl_vector_set(Q,i,gsl_vector_get(Q, i)/par_double);
  }
  for(i=0;i<Nt_2;i++)
  {
     fprintf(file_out,"%.15le\n",gsl_vector_get(Q, i));fflush(file_out);
  }
 ////
  par_double=0.0;
  for(i=0;i<Nt_2;i++)
  {
    par_double += gsl_vector_get(Q, i) * gsl_vector_get(R,i);
  }
  fprintf(file_out, "norma_new=%.15le\n", par_double);fflush(file_out);

  fclose(file_out); 




  gsl_matrix_free(W_1);
  gsl_vector_free(Wk);
  gsl_matrix_free(W);
}

void calculate_rho(gsl_vector* Q,  double* current, double* error, double* rho, double* rho_stat_err)
{
  int i,j,interval_length;
  double res, res_error;
  double avg,err_avg;
  res=0.0;
  res_error=0.0;
 for(i=0;i<Nt_2;i++)
  {
    interval_length=interval_numbers[i].size;
    for(j=0;j<interval_length;j++) {
      avg+=current[interval_numbers[i].times[j]];
      err_avg+=error[interval_numbers[i].times[j]]*error[interval_numbers[i].times[j]];
    }
    err_avg=sqrt(err_avg)/((double)interval_length);
    avg=avg/((double)interval_length);
    res+=avg*gsl_vector_get(Q, i);
    res_error+=err_avg*err_avg*gsl_vector_get(Q, i)*gsl_vector_get(Q, i);
    //res+=current[i]*gsl_vector_get(Q, i);
    //res_error+=error[i]*error[i]*gsl_vector_get(Q, i)*gsl_vector_get(Q, i);
  }
  res_error=sqrt(res_error);
  *rho=res;
  *rho_stat_err=res_error;
}

double delta_width_calculation(gsl_vector* Q, double center)
{
  double int_result, int_error;
  gsl_function F;
  gsl_integration_workspace * w1  = gsl_integration_workspace_alloc (N_int_steps);
  FILE* file_out;

  file_out=fopen_log("delta_width_control.txt","a", center);

  F.function = &delta_sq_int;
  d_parameters buffer;
  buffer.Q=Q;
  buffer.center=center;  
    F.params = &buffer;
    gsl_integration_qagiu (&F, 0.0, accuracy, accuracy, N_int_steps, w1, &int_result, &int_error);
    fprintf(file_out,"%.15le\t%.15le\t%.15le\n",center, int_result, int_error); fflush(file_out);
   gsl_integration_workspace_free (w1);
  fclose(file_out);
  return sqrt(int_result);
}


//calculation of real chracteristics of delta function
//returns all values in temperature units
void delta_characteristics_calculation(double* start, double* stop, double* center_real, gsl_vector* Q, double center)
{
  //delta function output
  FILE* file_out;
  double omega,omega_max, function, function_max;
  unsigned long int count=0;

  file_out=fopen_log("real_delta_width_control.txt","a", center);
  for(omega=omega_plot_delta/(2.0*dNt_2_pre);omega<omega_plot_limit/(2.0*dNt_2_pre);omega+=omega_plot_delta/(2.0*dNt_2_pre))
  {
    if(count==0)
      {
        function_max= delta(omega,Q);
        omega_max=omega;
      }
    else
      {
        if(function_max<delta(omega, Q))
        {
          function_max=delta(omega,Q);
          omega_max=omega;
        }
      }
    count++;  
  }
  fprintf(file_out,"omega_max=%.15le\n",omega_max*(2.0*dNt_2_pre)); fflush(file_out);
  (*start)=0.0;
 


  //start calculation
  for(omega=omega_max;omega>omega_plot_delta/(2.0*dNt_2_pre);omega-=omega_plot_delta/(2.0*dNt_2_pre))
  {
    if(delta(omega,Q)<function_max/2.0)
     {
       (*start)=omega;
       break;
     }
  }
  fprintf(file_out,"omega_start=%.15le\n",(*start)*(2.0*dNt_2_pre)); fflush(file_out);

//stop calculation
  for(omega=omega_max;omega<omega_plot_limit/(2.0*dNt_2_pre);omega+=omega_plot_delta/(2.0*dNt_2_pre))
  {
    if(delta(omega,Q)<function_max/2.0)
     {
       (*stop)=omega;
       break;
     }
  }
  fprintf(file_out,"omega_stop=%.15le\n",(*stop)*(2.0*dNt_2_pre)); fflush(file_out);
  (*center_real)=omega_max;;//(  (*start)   +   (*stop))/2.0;
  
  (*start)=(*start) * 2.0*dNt_2_pre;
  (*stop)=(*stop) * 2.0*dNt_2_pre;
  (*center_real)=(*center_real) * 2.0*dNt_2_pre;

}




int main(int np, char** p)
{
  double * current, *error;
  double * current_pre, *error_pre;//preliminary input
  int i,j,k,t;
  FILE* file_in_current;
  FILE* file_in_matrix;
  FILE* file_out;
  FILE* file_out_excl;
  FILE* file_const;
  int par_int;
  double par_double;
  double int_result, int_error;
  double center;
  double omega;
  
  
  gsl_vector * R;
  gsl_vector * omega_R; //for real center definition

  gsl_vector * Q;

  gsl_vector * Q_initial;

  gsl_matrix * S;//covariance matrix
  gsl_matrix * S_pre;//covariance matrix (preliminary input)

  file_in_current=fopen(p[2],"r");
  file_in_matrix=fopen(p[3],"r");
  sprintf(directory,"%s", p[4]);
  Nt_2=atoi(p[1]);
  dNt_2=(double) Nt_2;
//reading parameters file
  file_const=fopen(p[5],"r");
  
//double accuracy=1e-8;
parse_option(file_const,"%le",&accuracy);
//long int N_int_steps=1000000;
parse_option(file_const,"%ld",&N_int_steps);
//double omega_plot_delta=0.01;
parse_option(file_const,"%le",&omega_plot_delta);
//double omega_plot_limit=30.0;
parse_option(file_const,"%le",&omega_plot_limit);
//double lambda=1.0;
parse_option(file_const,"%le",&lambda);
//double center_start=3.0;
parse_option(file_const,"%le",&center_start);
//double center_stop=18.01;
parse_option(file_const,"%le",&center_stop);
//double center_delta=3.0;
parse_option(file_const,"%le",&center_delta);
//int flag_exclude
parse_option(file_const,"%d",&flag_exclude_delta);
//int count_start
parse_option(file_const,"%d",&count_start_exclude);
//int flag_model
parse_option(file_const,"%d",&flag_model);
if(flag_model==1)
{
  fscanf(file_const, "%d", &N_valid_points);
  points_numbers=(int*) calloc(N_valid_points, sizeof(int));
  for(i=0;i<N_valid_points; i++)
    fscanf(file_const, "%d", &(points_numbers[i]));
}
else
{
  N_valid_points=Nt_2;
  points_numbers=(int*) calloc(N_valid_points, sizeof(int));
  for(i=0;i<N_valid_points; i++)
    points_numbers[i]=i+1;
}
 N_intervals=N_valid_points;
 interval_numbers=(interval *)calloc(N_intervals, sizeof(interval));
 
 construct_intervals(points_numbers,N_valid_points,interval_numbers);

  fclose(file_const);

file_out=fopen_control("parameters.txt","w");
fprintf(file_out,"number of points (=half of physical number of timeslices)=%d\n", Nt_2);

fprintf(file_out,"file with current correlator:\n %s\n", p[2]);
fprintf(file_out,"file with covariance matrix:\n %s\n", p[3]);
fprintf(file_out,"path for output:\n %s\n", p[4]);
fprintf(file_out,"file with parameters:\n %s\n", p[5]);




//double accuracy=1e-8;
fprintf(file_out,"accuracy=%.15le\n",accuracy);
//long int N_int_steps=1000000;
fprintf(file_out,"Number os steps in numerical interation=%ld\n",N_int_steps);
//double omega_plot_delta=0.01;
fprintf(file_out,"Step in plots of resolution function(omega)=%.15le\n",omega_plot_delta);
//double omega_plot_limit=30.0;
fprintf(file_out,"Limit in plots of resolution function(omega)=%.15le\n",omega_plot_limit);
//double lambda=1.0;
fprintf(file_out,"Lambda(regularization)=%.15le\n(=1 - without regularization)\n",lambda);
//double center_start=3.0;
fprintf(file_out,"Center of resolution function starts from = %.15le\n",center_start);
//double center_stop=18.01;
fprintf(file_out,"Center of resolution function stops at %.15le\n",center_stop);
//double center_delta=3.0;
fprintf(file_out,"Step in center of resolution function = %.15le\n",center_delta);
//flag_exclude_delta;
fprintf(file_out,"TExclude or not (1 or 0) initial delta finction %d\n",flag_exclude_delta);
//count_start_exclude
fprintf(file_out,"Number of iteration for center to start exclusion of initial delta-function = %d\n",count_start_exclude);
//model_flag;
fprintf(file_out,"Take into account all (0)  or not all (1) timeslices = %d\n",flag_model);
if(flag_model==1)
{
  fprintf(file_out,"Number of timeslices taken into account = %d\n",N_valid_points);
  for(i=0;i<N_valid_points;i++)
  {
    fprintf(file_out,"%d\n",points_numbers[i]);
  }
}
fclose(file_out);


  current_pre=(double*) calloc(Nt_2, sizeof(double));
  error_pre=(double*) calloc(Nt_2, sizeof(double));

  S_pre=gsl_matrix_calloc(Nt_2, Nt_2);

  for(t=1;t<=Nt_2;t++)
  {
    fscanf(file_in_current, "%d", &par_int);
    fscanf(file_in_current, "%le", &par_double);
    current_pre[t-1]=par_double;
    
    fscanf(file_in_current, "%le", &par_double);
    error_pre[t-1]=par_double;
  }
  fclose(file_in_current);

  file_out=fopen_control("current_control_pre.txt","w");
  for(t=0;t<Nt_2;t++)
  {
    fprintf(file_out,"%d\t%.15le\t%.15le\n", t, current_pre[t], error_pre[t]);
  }
  fclose(file_out);

  for(i=1;i<=Nt_2;i++)
  for(t=1;t<=Nt_2;t++)
  {
    fscanf(file_in_matrix, "%le", &par_double);
    gsl_matrix_set(S_pre, i-1, t-1, par_double);

  }
  fclose(file_in_matrix);

  file_out=fopen_control("cov_matrix_control_pre.txt","w");
  for(i=0;i<Nt_2;i++){
  for(t=0;t<Nt_2;t++)
  {
    fprintf(file_out,"%.15le\t", gsl_matrix_get(S_pre, i,t));
  }
  fprintf(file_out,"\n");
  }
  fclose(file_out);
  //If we are going to do intervals: do we need full data set to construct new covariance matrix?
//conversion to real arrays taken into account only certain timeslices
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////
{ 
  int count1, count2;
  Nt_2_pre=Nt_2;
  dNt_2_pre=(double) Nt_2_pre;
  if(flag_model==1)
  {
    Nt_2=N_valid_points;
    dNt_2=(double)Nt_2;
  } 
  current=(double*) calloc(Nt_2, sizeof(double));
  error=(double*) calloc(Nt_2, sizeof(double));

  S=gsl_matrix_calloc(Nt_2, Nt_2);

  for(count1=0;count1<Nt_2;count1++)
  {
    current[count1]=current_pre[points_numbers[count1]-1];
    error[count1]=error_pre[points_numbers[count1]-1];
  }
  fclose(file_in_current);

  file_out=fopen_control("current_control_fin.txt","w");
  for(t=0;t<Nt_2;t++)
  {
    fprintf(file_out,"%d\t%.15le\t%.15le\n", points_numbers[t], current[t], error[t]);
  }
  fclose(file_out);

  double result=0.0;
  for(count1=0;count1<Nt_2;count1++)
  for(count2=0;count2<Nt_2;count2++)
  {
    for(i=0;i<interval_numbers[count1].size;i++) {
      for(j=0;j<interval_numbers[count2].size;j++) {
	par_double=gsl_matrix_get(S_pre,interval_numbers[count1].times[i]-1,interval_numbers[count2].times[j]-1);
	result+=par_double;
      }
    }
    par_double=result/((double)interval_numbers[count1].size*interval_numbers[count2].size);
    gsl_matrix_set(S,count1,count2,par_double);
    result=0.0;
    //par_double=gsl_matrix_get(S_pre,points_numbers[count1]-1, points_numbers[count2]-1);
    //gsl_matrix_set(S, count1, count2, par_double);
  }
  fclose(file_in_matrix);

  file_out=fopen_control("cov_matrix_control.txt","w");
  for(i=0;i<Nt_2;i++){
  for(t=0;t<Nt_2;t++)
  {
    fprintf(file_out,"%.15le\t", gsl_matrix_get(S, i,t));
  }
  fprintf(file_out,"\n");
  }
  fclose(file_out);
}
//end of conversion
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////



//integration to obtain R vector
  file_out=fopen_control("R_control.txt","w");
  
  R=gsl_vector_alloc(Nt_2);
{
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (N_int_steps);
  gsl_function F;
  F.function = &kernel_int;
  for(t=1;t<=Nt_2;t++)
  {
    F.params = &t;
    gsl_integration_qagiu (&F, 0.0, accuracy, accuracy, N_int_steps, w, &int_result, &int_error); 
    gsl_vector_set(R,t-1,int_result);
    fprintf(file_out,"%d\t%.15le\t%.15le\n", t, int_result, int_error); fflush(file_out);
  }
  gsl_integration_workspace_free (w);

}
  fclose(file_out);
///////
//integration to obtain omega_R vector
  file_out=fopen_control("omega_R_control.txt","w");
  
  omega_R=gsl_vector_alloc(Nt_2);
{
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (N_int_steps);
  gsl_function F;
  F.function = &omega_kernel_int;
  for(t=1;t<=Nt_2;t++)
  {
    F.params = &t;
    gsl_integration_qagiu (&F, 0.0, accuracy, accuracy, N_int_steps, w, &int_result, &int_error); 
    gsl_vector_set(omega_R,t-1,int_result);
    fprintf(file_out,"%d\t%.15le\t%.15le\n", t, int_result, int_error); fflush(file_out);
  }
  gsl_integration_workspace_free (w);

}
  fclose(file_out);
///////


  {
  int count_center=0;
  double C, C0;
  gsl_vector* Q_real;
  Q_initial=gsl_vector_calloc(Nt_2);
  for(center=center_start; center<=center_stop; center+=center_delta)
  {
    double center_real=center/(2.0*dNt_2_pre);
    double center_calculated=0.0;
    double center_calculated_real=0.0;
    char file_name[1000];
    sprintf(file_name,"delta_function_c=%3.3leT.txt", center);
    Q=gsl_vector_calloc(Nt_2);
    Q_real=gsl_vector_calloc(Nt_2);
    calculate_Q(Q, R, S, center_real);
    if(count_center==0)
    {
      for(i=0;i<Nt_2;i++)
        gsl_vector_set(Q_initial, i, gsl_vector_get(Q,i));
    }
    if(count_center==0)
    {
      for(i=0;i<Nt_2;i++)
        gsl_vector_set(Q_real, i, gsl_vector_get(Q,i));
    }
    else
    {
      if(flag_exclude_delta==1 && count_center>=count_start_exclude)
      {
        FILE* log_exclusion;
        log_exclusion=fopen_log("log_exclusion.txt","a", center_real);
        C=delta0(Q);
        C0=delta0(Q_initial);
        fprintf(log_exclusion,"C=%.15le\nC0=%.15le\n",C,C0);fflush(log_exclusion);
        
        for(i=0;i<Nt_2;i++)
          gsl_vector_set(Q_real, i, gsl_vector_get(Q,i)-(C/C0)*gsl_vector_get(Q_initial,i));
        
        //normalization
        double new_norma=0.0;
        for(i=0;i<Nt_2;i++)
        {
          new_norma += gsl_vector_get(Q_real, i) * gsl_vector_get(R,i);
        }
        fprintf(log_exclusion, "norma_old=%.15le\n", new_norma);fflush(log_exclusion);
        
        for(i=0;i<Nt_2;i++)
        {
          gsl_vector_set(Q_real,i,gsl_vector_get(Q_real, i)/new_norma );
        }
        new_norma=0.0; 
        for(i=0;i<Nt_2;i++)
        {
          new_norma += gsl_vector_get(Q_real, i) * gsl_vector_get(R,i);
        }
        fprintf(log_exclusion, "norma_new=%.15le\n", new_norma);fflush(log_exclusion);
        
        fclose(log_exclusion);     
      }
      else
      {
        for(i=0;i<Nt_2;i++)
        gsl_vector_set(Q_real, i, gsl_vector_get(Q,i));
      }
    }
    //delta function output
    file_out=fopen_control(file_name,"w");
    for(omega=0;omega<omega_plot_limit/(2.0*dNt_2_pre);omega+=omega_plot_delta/(2.0*dNt_2_pre))
    {
      fprintf(file_out,"%.15le\t%.15le\t%.15le\n", omega*2.0*dNt_2_pre, delta(omega,Q), delta(omega, Q_real));
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
    //calculate_rho(Q, current, error, &rho, &rho_stat_err);
    calculate_rho(Q, current_pre, error_pre, &rho, &rho_stat_err);
    width=delta_width_calculation(Q, center_real);

    delta_characteristics_calculation(&start, &stop, &center1, Q, center_real);
    width1=(stop-start)/2.0;

    //calculate_rho(Q_real, current, error, &rho_real, &rho_stat_err_real);
    calculate_rho(Q_real, current_pre, error_pre, &rho_real, &rho_stat_err_real);
    width_real=delta_width_calculation(Q_real, center_real);

    delta_characteristics_calculation(&real_start, &real_stop, &real_center1, Q_real, center_real);
    real_width1=(real_stop-real_start)/2.0;


    fprintf(file_out,"%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\n", center, width*2.0*dNt_2_pre, center1, width1, rho, rho_stat_err);
    fprintf(file_out_excl,"%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\t%.15le\n", center, width_real*2.0*dNt_2_pre, real_center1, real_width1,  rho_real, rho_stat_err_real, real_start, real_stop);
    
    fclose(file_out);
    fclose(file_out_excl);
    gsl_vector_free(Q);
    gsl_vector_free(Q_real);
    count_center++;
  }
  gsl_vector_free(Q_initial);
  }
  gsl_vector_free(R);
  gsl_vector_free(omega_R);
  gsl_matrix_free(S);
  gsl_matrix_free(S_pre);
  
  free(current);
  free(error);

  free(current_pre);
  free(error_pre);

  if(flag_model==1)
    free(points_numbers);
  return 0;
}

