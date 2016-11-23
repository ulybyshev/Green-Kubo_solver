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


extern double accuracy;

extern long int N_int_steps;

extern double omega_plot_delta;

extern double omega_plot_limit;

extern int flag_lambda_regularization;

extern double lambda;

extern double center_start;

extern double center_stop;

extern double center_delta;

extern int flag_model;

//extern int N_valid_points;

//extern int* points_numbers;

//extern int N_intervals;
//extern interval *interval_numbers;

//extern int Nt_2_pre;
//extern double dNt_2_pre;

extern int flag_exclude_delta;

extern int count_start_exclude;



#endif
