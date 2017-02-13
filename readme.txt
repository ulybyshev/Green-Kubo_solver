The programs performs analytical continuation from Euclidean to real time through the solution of Green-Kubo relation. 
Example calculation of Density of States for Square Hubbard model from Monte Carlo data is included.

____________________________________________________________________________________________________________________
COMPILATION

1) requires GNU GSL library: path to the library should be  inserted in  LIB and INCLUDE_DIR variables in Makefile
2) make all


____________________________________________________________________________________________________________________
COMMAND LINE OPTIONS


without -a option: enters the regime with  error estimation through data blocking, needs the full set of data (correlators for each monte-carlo (MC) configuration)
    in this regime:
    -b Nb  (Nb=number of data blocks) - if not provided by the user then the number of data blocks is set up automatically 
    -c filename.txt  - file with full set of raw data: full-time Euclidean correlator for each MC configuration
	format:
	0    correlator_real_part   correlator_imaginary_part
	1 ....
	2 ....
	Nt_full-1 ......
	(empty line)
	......(data for the next configuration)

with -a option: enters the regime without data blocking, which works with the average correlator and covariance matrix:
    -c  filename.txt   - file with correlator  (it is assuned that the correlator is symmetrical with respect to the half of Euclidean time)
	format:
	1    correlator    correlator_error 
	2 ....
	3 ....
	4 ....
	Nt ....
	
	The input of errors is optional: if -e option exists, the format is the following
	1    correlator
	2 ....
	3 ....
	4 ....
	Nt ....
	
	But in this case the error estimation and the automated procedures for choosing regularization don't work
	
    -m  filename.txt -file with covariance matrix   -   the input of covariance matrix is optional, can work without it, but then the regularization with covariance matrix doesn't work
	format:
	C_{11} .......  C_{1 Nt}
        ..........................
	C_{Nt 1} .......  C_{Nt Nt}
	(without any indexes)

joint parameters for both regimes:
-t  Nt  (Nt  = HALF of physical timeslices  in real simulations)
-o path_for_output_data
-p filename.txt  - file with parameters (if needed: without it  all parameters are equal to their default values)

____________________________________________________________________________________________________________________
FILE WITH PARAMETERS (formatting)


    parameter_name
    parameter_value (or values in case of multiple ones)

    parameter names:
    
  "kernel_switcher"
    possible values:  
    ==0 for conductitvity kernel
    ==1 for Density of States kernel	
    ==2 for Density of States kernel taking into account discrete time
    ==3 the same as for conductivity but with cosh(omega*beta/2) in denominator
 
  
  "accuracy"
    floating point number: relative accuracy of numerical integration

  "N_int_steps"
    integer: maxium number of steps in numberical integration

  "delta_plot_step"
    floating point number: the step size in the plots for resolution functions
  
  "delta_plot_limit"
    floating point number: maximal value of frequency in the plots for resolution functions
    
  "center_start"
    floating point number: starting value of the center of resolution function
  
  "center_stop"
    floating point number: maximal value of the center of resolution function
  
  "center_delta"
    floating point number: the step size for the center of resolution functions 
  
  "lambda"
    format:
    integer_number     floating_point_parameter
	integer_number:  
	==0 without regularization (thus no floating point parameter)
	==+-1 regularization by addition of covariance matrix:  (1-\lambda) S_{ij} 
	    +1 floating_point_parameter=lambda
	    -1 floating_point_parameter=average relative error  (and  lambda is choosen automatically)
	==+-2 regularization by neglecting all eigenvalues of W kernel less than \lambda  (Truncated SVD decomposition)
	    +2 floating_point_parameter=lambda
	    -2 floating_point_parameter=average relative error  (and  lambda is choosen automatically)
	==3 regularization by Tikhonov filtering in SVD decomposition  (minimize |K x-y| +lambda |x| (suppress fluctuations in spectral funtions)
	    +3 floating_point_parameter=lambda
	    -3 floating_point_parameter=average relative error  (and  lambda is choosen automatically)
	==4 alternative regularization through Tikhonov filtering in SVD  decomposition   (minimize  |K x-y| +lambda |Kx| (suppress fluctuations in correlator obtained from calculated spectral function)
	    +4 floating_point_parameter=lambda
	    -4 floating_point_parameter=average relative error  (and  lambda is choosen automatically)

  "flag_exclude_corr" additional regularization by neglecting some points in correlator or by averaging over some intervals in correlator
    format:  
    ==0 no additional regularization
    ==1 additional regularization is introduced: all points except the listed ones are neglected
	format:
	    1   N_{valid_timeslices}    i1   i2  .....  i_{N_{valid}}  
	    (numbers should be sorted in ascending order, contact term in correlator is under number 0, so valid timeslices start from 1 and proceed till Nt)
    ==2 the same as in the previous case, but correlator is averaged over intervals between listed points

example of the file with parameters

kernel_switcher
1
accuracy
1.0e-12
N_int_steps
1000000
delta_plot_step
0.1
delta_plot_limit
80.0
center_start
0.0
center_stop
40.0
center_delta
0.5
lambda
-4  0.1
flag_exclude_corr
1  10  1 2 5 10 15 20 25 30 35 40

all dimensional parameteres are in the units of temperature


Default values of parameters:
regularization type:				-4
relative error:   				0.05
value of regularization constant lambda:	1.0e-06
kernel type					1
accuracy 					1.0e-12
number of steps in numerical integration 	1000000
delta_plot_step					0.25
delta_plot_limit				Nt *2.0
center_start 					0.0
center_delta 					1.0
center_stop=					Nt

____________________________________________________________________________________________________________________
OUTPUT

all dimensional quantities are in the units of temperature

1)	correlator_control_pre_N.txt
	Preliminary values of correlator for N-th bin (before averaging over intervals or neglecting some data points)
	Format:
	#time  correlator_Re   error_Re

2)	correlator_control_intervals_N.txt   of  correlator_control_fin_N.txt 
	Final values of correlator for N-th bin (after averaging over intervals or neglecting some data points, if this type of regularization was ordered in parameters file)
	Format:
	
	correlator_control_fin_N.txt 
	#time  correlator_Re   error_Re

	correlator_control_intervals_N.txt 
	#interval_number	average_correlator_Re  error_Re

3)	cov_matrix_control_pre_N.txt and cov_matrix_control_fin.txt: 
	covariance matrix for correlators from correlator_control_pre_N.txt and correlator_control_fin_N.txt (or correlator_control_intervals_N.txt )
	


4)      correlataion_study_timeN.txt
	input data for calculation of autocorrelation length for N-th time slice.
	Format:
	#binsize   bin_error/initial_error
	
	bin_error is statistical error calculated after data binning, initial error is ordinary statistical error. This ratio, in the limit of large binsize, is equal to
	autocorrelation length.
	 
5)	delta_function_c=_VALUE_T.txt
	resolution function with center at frequency = _VALUE_ (in units of temperature)
	Format:
	#frequency(in units of T)  resolution_function
	
	resolution functions are saved only for the final value of regularization parameter lambda (when it's tuned to get some fixed value of relative error)

6)	rho_final.txt
	final output for spectral function (at final value of regularization parameter lambda)
	Format:
	#center_of_resolution_function   spectral_function    error
	
7)      rho_basic.txt
	the same as rho_final.txt but includes more technical information about resolution functions
	In case of data blocking this file appears for each block of data separately for the final value of regularization parameter lambda
	Format:
	#center_of_resolution_function  resolution_width(as it was in the  calculation of D functional) \\
	position_of_the_maximum_of_resolution_function  width=(start-stop)/2  spectral_function  error  start  stop

	start and stop points are half-max points before and after maximum of resolution function
	

8)	rho_lambda_VALUE_avg.txt
	spectral function (at some value of regularization parameter lambda in the vicinity of the final one)
	Format:
	#center_of_resolution_function   spectral_function    error


9)      rho_lambda.txt
	dependence of recolution function on regularization parameter lambda
	Format:
	#value_of_lambda	rho_for_c=0.000e+00 T	rho_error  ...... (for all positions of the center of resolution functions)
	
	
	
