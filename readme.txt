The programs performs analytical continuation from Euclidean to real time through the solution of Green-Kubo relation. 
Example calculation os Density of States for Square Hubbard model from Monte Carlo data is included.


to run the program:   ./program  -t Nt -c  file_current.txt  -m file_matrix.txt -o path_for_output -p  file_with_parameters.txt

  file for currents:
  1    correlator    correlator_error 
  2 ....
  3 ....
  4 ....
  Nt ....

Nts = HALF of physical timeslices  in real simulations,
correlator should be symmetrized before input
value of correlator at zero euclidean time should be excluded from input
errors  from  correlator file are taken into account in calculation of statistical errors of spectral function 
covariance matrix should be calculated for symmetrized correlator (size of the matrix thus is also reduced)

 input for matrix is formatted as
 C_{11} .......  C_{1 Nt}
 ..........................
 C_{Nt 1} .......  C_{Nt Nt}
 (without any indexes)


 file with parameters (comments lines alternate with parameter values). All dimensional quantities are in units of temperature
 int kernel_switcher (=0 for conductivity kernel, ==1 for Density of States kernel) 
 1
 double accuracy
 1e-8
 long int N_int_steps
 1000000
 double omega_plot_delta
 0.01
 double omega_plot_limit
 30.0
 double center_start
 0.0
 double center_stop
 18.01
 double center_delta
 3.0
 lambda regularization
 1  0.9999
 int flag model   (numbers sorted ascending, contact term in correlator is under number 0, so valid timeslices start from 1 and go till Nt)
 1   N_{valid_timeslices}    i1   i2  .....  i_{N_{valid}}  
 int flag_exclude_delta_function (0 if not exclude - then the next parameter is not taken into account, 1 - exclude starting from the count_start iteration in center) and count_start
 1  2
 
 
 
output - all dimensional quantities in units of temperature

output - delta functions
omega(in units of temperature)    function   function_forced_to_be_zero_at_omwga_zero

output rho
center(as in was in input)  resolution_width(as it was in D calculation)  center(Maximum)  width((start-stop)/2, start and stop are at half-peak)  rho  rho_stat_err start stop

output rho_exclude - the same for resolution functions forced to be zero at zero frequency
center(as in was in input)  resolution_width(as it was in D calculation)  center(Maximum)  width((start-stop)/2, start and stop are at half-peak)  rho  rho_stat_err  start  stop

