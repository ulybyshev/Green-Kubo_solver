#include "constants.h"
//real definition of global constants and default values

/////////////////////////////
//constants from command line options
char output_directory[1024];
const char default_directory[]="./";

char correlator_filename[1024];
char cov_matrix_filename[1024];
char parameters_filename[1024];

bool flag_tune_blocking;
bool flag_constants_file;

int Nt_2;//number of points in correlator
double dNt_2;//the same in double precision

int flag_jackknife;
//==0 proceeds in normal way
//==1 computes analytic continuation for jackknife samples
int num_jack_samples;

////////////////////////////////
//parameters from constants file
int kernel_switcher;
//==0 for conductitvity kernel
//==1 for Density of States kernel
//==2 for Density of States kernel taking into account discrete time

//relative accuracy of numerical integration
double accuracy;

//maxium number of steps in numberical integration
long int N_int_steps;

//all dimensional quantities are in units of temperature. Correlator is assumed to be symmetrized so Nt_2 = corresponds to the half of inverse temperature

//step in plots for resolution function
double omega_plot_delta;

//maximal value for plots for resolution functions
double omega_plot_limit;

//==0 if  L-regularization is not introduced
//==1 if regularization through addition of covariance matrix is introduced  (1-\lambda) S_{ij}    (and -1 if lambda should be choosen automatically)
//==2 if regularization through neglecting all eigenvalues of W kernel less than \lambda  (and -2 if lambda should be choosen automatically) (Truncated SVD decomposition)
//==3 if regularization through Tikhonov filtering in SVD (and -3 if lambda should be choosen automatically) (minimize |K x-y| +lambda |x| (suppress fluctuations in spectral funtions)
//==4 if alternative regularization through Tikhonov filtering in SVD (and -4 if lambda should be choosen automatically) (minimize  |K x-y| +lambda |Kx| (sppress fluctuations in correlator obtained from calculated spetral function)
int flag_lambda_regularization;

//regularization constant
double lambda;

double relative_error;//if some value of relative error is ordered and the program should set proper lambda value automatically

//initial position of the center of resolution functions
double center_start;

//final  position of the center of resolution functions
double center_stop;

//the step in positioning of the center of resolution functions
double center_delta;

//flag_model==0 if we don't introduce regulariztion with neglecting of some points of correaltor
//==1 if we neglect some points
//==2 if we take into account averages over intervals
int flag_model;

//==0 if we impose additional requirement on resolution functions to be zero at zero frequency
int flag_exclude_delta;

//number of resolution function from which we start to impose this additional requirement
int count_start_exclude;

//global constant for performing the output in various log files
bool flag_log_output;
//for the most important logs
bool special_flag_log_output;

int limit_power=-1;//corresponding to the minimal value of lambda=0.9 = 1-10 ^{-1} for regularization with covariance matrix

int Number_lambda_points=15;//number of point in plots for lambda (1/3 of them for smaller regularization and 2/3 for larger)

int n_conf; //number of configurations
double dn_conf;

double *raw_data;


//presence of covariance matrix in input data
bool flag_covariance_matrix_input=true;
//presence of the errors in correlator input (makes sence only in regime without data blocking)
bool flag_error_corr_input=true;
