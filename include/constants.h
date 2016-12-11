//declaration of global constants and default values 
//detailed description of variables see in constants.txt
#ifndef   __GKS_CONSTANTS
#define   __GKS_CONSTANTS
#include "basic_structures.h"

//constants from command line options
extern char output_directory[1024];
extern const char default_directory[];

extern char correlator_filename[1024];
extern char cov_matrix_filename[1024];
extern char parameters_filename[1024];

extern int Nt_2;//number of points in correlator
extern double dNt_2;//the same in double precision

extern int n_conf;
extern double dn_conf;

extern int kernel_switcher;

extern int flag_jackknife;

extern int num_jack_samples;

extern double accuracy;

extern long int N_int_steps;

extern double omega_plot_delta;

extern double omega_plot_limit;

extern int flag_lambda_regularization;

extern double lambda;

extern double relative_error;

extern double center_start;

extern double center_stop;

extern double center_delta;

extern int flag_model;

extern int flag_exclude_delta;

extern int count_start_exclude;

extern bool flag_log_output;

extern bool special_flag_log_output;


extern int limit_power;//corresponding to the minimal value of lambda=0.9 = 1-10 ^{-1} for regularization with covariance matrix
extern int Number_lambda_points;//number of point in plots for lambda

extern double *raw_data;

#endif
