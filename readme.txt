The programs performs analytical continuation from Euclidean to real time through the solution of Green-Kubo relation. 
Example calculation of Density of States for Square Hubbard model from Monte Carlo data is included.


to run the program:
command line options:
without -a option: enters the regime with  error estimation through the data blocking, needs full set of data (correlators for each monte-carlo (MC) configuration)
    in this regime:
    -b Nb  (Nb=number of data blocks) - if doesnt't exist then number of data blocks is set up automatically (10 is preferrable value) 
    -c filename.txt  - file with full set of raw data: full-time Euclidean correlator for each MC configuration
	format:
	0    correlator_real_part   correlator_imaginary_part
	1 ....
	2 ....
	Nt_full-1 ......
	(empty line)
	......(data for the next configuration)

without data blocking (works with only average data and covariance matrix) (with -a option):
    -t  Nt  (Nt  = HALF of physical timeslices  in real simulations)
    -c  filename.txt   - file with correlator  (it is assuned that the correlator is symmetrical with respect to the half of Euclidean time)
	format:
	1    correlator    correlator_error 
	2 ....
	3 ....
	4 ....
	Nt ....
    -m  filename.txt -file with covariance matrix
	format:
	C_{11} .......  C_{1 Nt}
        ..........................
	C_{Nt 1} .......  C_{Nt Nt}
	(without any indexes)

joint parameters for bith regimes:
-o path_for_output_data
-p filename.txt  - file with parameters (if needed: without it all parameters are equal to default values)
    format:


    parameter_name
    parameter_value (or values in case of multiple ones)

    parameter names:
  KERNEL_SWITCHER_OPTION  "kernel_switcher"
  ACCURACY_OPTION  "accuracy"
  INTEGRAL_NUMBER_STEPS_OPTION   "N_int_steps"
  RES_FUNCTION_PLOT_STEP_OPTION   "delta_plot_step"
  RES_FUNCTION_PLOT_LIMIT_OPTION   "delta_plot_limit"
  RES_FUNCTION_CENTER_START_OPTION   "center_start"
  RES_FUNCTION_CENTER_STOP_OPTION   "center_stop"
  RES_FUNCTION_CENTER_DELTA_OPTION   "center_delta"
  REGULARIZATION_OPTION                      "lambda"
  EXCLUDE_CORR_POINTS_OPTION   "flag_exclude_corr"
  RES_FUNCTION_ZERO_AT_ZERO_OMEGA_OPTION   "flag_force_zero"



    Example:
    
    kernel_switcher
    1
    accuracy
    1e-8
    N_int_steps
    1000000
    delta_plot_step
    0.01
    delta_plot_limit
    30.0
    center_start
    0.0
    center_stop
    18.01
    center_delta
    3.0
    lambda
    1  0.9999
    flag_exclude_corr   (numbers should be sorted ascending, contact term in correlator is under number 0, so valid timeslices start from 1 and proceed till Nt)
    1   N_{valid_timeslices}    i1   i2  .....  i_{N_{valid}}  
    flag_force_zero (0 if not exclude - then the next parameter is not taken into account, 1 - exclude starting from the count_start iteration in center) and count_start
    1  2
 

More detailes about setting up the parameters
 
 lambda regularization: the following set up is possible:
  [0] - without regularization, [1 lambda] - reg. with covariance matrix   [2 lambda] - we neglect eigenvalues of W lees than lambda,
   [-1 rel_error]  or [-2 rel_error] - when program tries to set up the lambda to fit the stated value of relative error
   



 
output - all dimensional quantities in units of temperature

output - delta functions
omega(in units of temperature)    function   function_forced_to_be_zero_at_omwga_zero

output rho
center(as in was in input)  resolution_width(as it was in D calculation)  center(Maximum)  width((start-stop)/2, start and stop are at half-peak)  rho  rho_stat_err start stop

output rho_exclude - the same for resolution functions forced to be zero at zero frequency
center(as in was in input)  resolution_width(as it was in D calculation)  center(Maximum)  width((start-stop)/2, start and stop are at half-peak)  rho  rho_stat_err  start  stop

