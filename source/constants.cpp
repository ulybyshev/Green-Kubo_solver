#include "constants.h"
//real definition of global constants and default values

/////////////////////////////
//constants from command line options
char output_directory[1024];
const char default_directory[]="./";

char correlator_filename[1024];
char cov_matrix_filename[1024];
char parameters_filename[1024];

int Nt_2;//number of points in correlator
double dNt_2;//the same in double precision

////////////////////////////////
//parameters from constants file
int kernel_switcher;
//==0 for conductitvity kernel
//==1 for Density of States kernel


//relative accuracy of numerical integration
double accuracy=1e-8;

//maxium number of steps in numberical integration
long int N_int_steps=1000000;

//all dimensional quantities are in units of temperature. Correlator is assumed to be symmetrized so Nt_2 = corresponds to the half of inverse temperature

//step in plots for resolution function
double omega_plot_delta=0.01;

//maximal value for plots for resolution functions
double omega_plot_limit=30.0;

//==0 if  L-regularization is not introduced
//==1 if regularization through addition of covariance matrix is introduced  (1-\lambda) S_{ij}
//==2 if regularization through neglecting all eigenvalues of W kernel less than \lambda
int flag_lambda_regularization=0;

//regularization constant
double lambda=1.0;

//initial position of the center of resolution functions
double center_start=0.00;

//final  position of the center of resolution functions
double center_stop=18.01;

//the step in positioning of the center of resolution functions
double center_delta=3.0;

//flag_model==0 if we don't introduce regulariztion with neglecting of some points of correaltor
//==1 if we neglect some points
//==2 if we take into account averages over intervals
int flag_model=0;

//number of points which is really taken into account (in case of intervals one interval is counted as one point)
//int N_valid_points;

//array of point numbers which are taken into account
//int* points_numbers;

//int N_intervals;
//interval *interval_numbers;

//int Nt_2_pre;
//double dNt_2_pre;

//==0 if we impose additional requirement on resolution functions to be zero at zero frequency
int flag_exclude_delta=0;

//number of resolution function from which we start to impose this additional requirement
int count_start_exclude;

