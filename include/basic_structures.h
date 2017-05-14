#ifndef   __GKS_BASIC
#define   __GKS_BASIC

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>

#define PI 3.1415926

class interval
{
public:  
  int size;
  int *times;
  interval(int size_interval=1);
  void format(int size_interval);
  ~interval();
};

class correlator
{
public:
    int N_full_points; //equal to number of points in input file for correlator
    int N_valid_points;
    double length;
//all points    
    double* corr_full;
    double* error_full;
    gsl_matrix* S_full;

//only valid points
    double* corr;
    double* error;
    gsl_matrix* S;
    

    int* points_numbers;
    int* points_stop_numbers;

    interval **interval_numbers;    
    
    correlator(int N_full=1, int N_valid=1);
    void format(int N_full,int N_valid);
    ~correlator();
    void construct_intervals();
};


class calc_structures
{
public:
    int N_t_points;//number of valid on correlator - equal to the size of W matrix
    double length;
    gsl_vector * R;
    gsl_vector * omega_R; 
    
    int N_center;//number of resolution functions
    double* center;//arrays for values of the centersof resolution functions
    
    gsl_matrix** W;//array of pointers for the "integral" part of W matrix

    calc_structures(int N_t);//values of center_start and center_stop are global
    ~calc_structures();
};

class spectral_functions
{
public:
    calc_structures* pCS;
    int N_lambda;
    double* lambda_array;
    double** rho_array;
    double** rho_err_array;
    
    spectral_functions(int N_lambda_input, double lambda_base, calc_structures* pCS);
    void format(int N_lambda_input,double lambda_base, calc_structures* pCS);
    ~spectral_functions();
    
};

class initial_data_description
{
public:
    int N_histories;
    double* corr_lengths;
    int* times;
    double largest_corr_length;
    initial_data_description(int N_histories_in=1);
    void format(int N_histories_in=1);
    ~initial_data_description();
    double largest_corr_length_calc();

};


//in all parameters structures  we have the number of timeslice in array, NOT real time
struct kernel_parameters
{
    int t;//number of point in time (if flag_model==1) or number of interval if flag_model==2
    correlator* pC;
};

struct w_parameters
{
  double omega_center;
  int i;
  int j;
  correlator* pC;
};

struct d_parameters
{
  double center;
  gsl_vector* Q;
  correlator* pC;
};



#endif
