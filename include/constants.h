//declaration of global constants and default values 
//detailed description of variables see in constants.txt
#ifndef   __GKS_CONSTANTS
#define   __GKS_CONSTANTS
#include "basic_structures.h"

//options for constants file
#define KERNEL_SWITCHER_OPTION  "kernel_switcher"
#define ACCURACY_OPTION  "accuracy"

#define INTEGRAL_NUMBER_STEPS_OPTION   "N_int_steps"
#define RES_FUNCTION_PLOT_STEP_OPTION   "delta_plot_step"
#define RES_FUNCTION_PLOT_LIMIT_OPTION   "delta_plot_limit"
#define RES_FUNCTION_CENTER_START_OPTION   "center_start"
#define RES_FUNCTION_CENTER_STOP_OPTION   "center_stop"
#define RES_FUNCTION_CENTER_DELTA_OPTION   "center_delta"
#define REGULARIZATION_OPTION                      "lambda"
#define EXCLUDE_CORR_POINTS_OPTION   "flag_exclude_corr"
#define RES_FUNCTION_ZERO_AT_ZERO_OMEGA_OPTION   "flag_force_zero"

//default values

#define  DEFAULT_FLAG_REGULARIZATION	-4
#define  DEFAULT_RELATIVE_ERROR    	0.05
#define  DEFAULT_LAMBDA 		1.0e-06
#define  DEFAULT_FLAG_MODEL   		0
#define  DEFAULT_FLAG_EXCLUDE_DELTA	0
#define  DEFAULT_COUNT_START_EXCLUDE	0
#define  DEFAULT_KERNEL_SWITCHER	1
#define  DEFAULT_ACCURACY 		1.0e-12
#define  DEFAULT_N_INT_STEPS 		1000000
#define  DEFAULT_OMEGA_PLOT_DELTA	0.5
#define  DEFAULT_CENTER_START 		0.0
#define  DEFAULT_CENTER_DELTA 		1.0
#define  DEFAULT_FLAG_EXCLUDE_DELTA	0

#define N_BINS_MINIMUM			5
#define N_BINS_MAXIMUM			20
#define N_CONF_IN_BIN			50


//constants from command line options
extern char output_directory[1024];
extern const char default_directory[];

extern char correlator_filename[1024];
extern char cov_matrix_filename[1024];
extern char parameters_filename[1024];

extern bool flag_constants_file;
extern bool flag_tune_blocking;

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

extern bool flag_covariance_matrix_input;
extern bool flag_error_corr_input;


#endif
