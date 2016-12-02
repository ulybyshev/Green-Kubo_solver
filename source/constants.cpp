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
//==2 for Density of States kernel taking into account discrete time

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
//==1 if regularization through addition of covariance matrix is introduced  (1-\lambda) S_{ij}    (and -1 if lambda should be choosen automatically)
//==2 if regularization through neglecting all eigenvalues of W kernel less than \lambda  (and -2 if lambda should be choosen automatically)
int flag_lambda_regularization=0;

//regularization constant
double lambda=1.0;

double relative_error=0.1;//if some value of relative error is ordered and the program should set proper lambda value automatically

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

//==0 if we impose additional requirement on resolution functions to be zero at zero frequency
int flag_exclude_delta=0;

//number of resolution function from which we start to impose this additional requirement
int count_start_exclude;

//global constant for performing the output in various log files
bool flag_log_output;
//for the most important logs
bool special_flag_log_output;

int limit_power=-1;//corresponding to the minimal value of lambda=0.9 = 1-10 ^{-1} for regularization with covariance matrix

