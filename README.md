[![DOI](https://zenodo.org/badge/50623371.svg)](https://zenodo.org/badge/latestdoi/50623371)

# Solver for Green-Kubo relations

The programs performs the analytic continuation from Euclidean to real time through the solution of the Green-Kubo relation. 
The modified Backus-Gilbert method is empoyed to regularize this ill-posed problem. 
The solver can work with several kernels in the Green-Kubo relations including kernel for Density of States (DOS) and conductivity (see below). 
An example calculation for the Density of States using Monte Carlo data for the square lattice Hubbard model is included.


## Compilation

The program requires the GNU GSL library. The path to the library should be inserted in LIB and INCLUDE_DIR variables in Makefile.
Once the path is set up, compilation is performed by typing

```
$ make all
```

## Command line options

The program can work in two distinct regimes: 1) errors of the spectral function are estimated through data binning 
(in this regime it needs the full Monte Carlo time series for the data)
and 2) errors of the spectral function are estimated on the basis of the average value of the Euclidean correlator and its errors (ideally, 
in this regime the covariance matrix of the correlator should be provided).

-a : If this option is present the program enters the 2nd regime. If it's absent, the 1st regime with data binning is activated.

### Common options for both regimes:

-t  Nt : Nt is positive integer number equal to half the number of Euclidean timeslices in the initial correlator. 
-o path : With this option the path where the output data is placed can be defined.
-p filename.txt :  If one wants to tune the internal parameters of the algorithm, the path to 
the file which contains these parameters can be defined here.  Without this option all parameters are equal to their default values. 
The format of the file with parameters is decribed below.

### Additional options in the regime with data binning:

-b Nb : Nb is a positive integer which sets the number of bins for the data binning procedure. If this number is not provided 
by the user then the number of bins is set up automatically on the basis of autocorrelation length measurement.
-c filename.txt : file with full Monte Carlo time series for the Euclidean correlator. File format:

```
	0    		Real_Part   Imaginary_Part
	1    		....
	2    		....
	....................
	2 Nt-1		....
	(empty line)
	......(data for the next configuration in the same format)
```
	
-i : This options turns off the imaginary part of the correlator. In this case the input file is organized in two columns.

Depending on the kernel, the program performs symmetrization or antisymmetrization of the input data with respect to \beta/2 where \beta is the inverse temperature. 
 If kernel_switcher=5 (see below), it  performs antisymmetrization, in all other cases  symmetrization is performed. 


### Additional options in the regime without data binning:

-c  filename.txt : File with correlator  (it is assuned that the correlator is symmetrical/antisymmetrical with respect to \beta/2).
```
	1    correlator    correlator_error 
	2 ....
	3 ....
	4 ....
	Nt ....
```
-e : If this option exists, the input of errors is canceled (only two columns in the file with correlator data).
In this case the error estimation and the automated procedures for choosing the regularization parameter don't work.
-m  filename.txt : File with covariance matrix.   The input of covariance matrix is optional. The program can work without it 
but then the regularization with covariance matrix doesn't work. File format:

```
	C_{11} ..........  C_{1 Nt}
        ..........................
	C_{Nt 1} .......  C_{Nt Nt}
```

## File with parameters

In this file it's possible to tune some internal parameters of the algorithm, e.g. to change the range of spectral function, choose the regularization 
or the kernel in Green-Kubo relation. Generally, the file consists of the following two-string expressions:
```
    parameter_name
    parameter_value (or several values if needed)
```
All dimensional parameteres are assumed to be in units of temperature.

### The list of all possible parameters.

"kernel_switcher" : Changes the kernel in Green-Kubo relations. Possible values:  
  * 0 : kernel for conductivity
  * 1 : kernel for Density of States (DOS is assumed to be symmetrical with respect to zero and the correlator is symmetrical with respect to half of inverse temperature)
  * 2 : kernel for Density of States kernel taking into account discrete Euclidean time
  * 4 : lattice version of the kernel for conductitvity
  * 5 : odd kernel for Density of States (is used to compute antisymmetrical woth respect to zero part of DOS, the correlator is assumed to be antisymmetrical with respect to half of inverse temperature)
 
  
"accuracy" : Parameter value is some floating point number. It defines the relative accuracy of numerical integration.

"N_int_steps" : Parameter value is some integer number. It defines maxium number of steps in numerical integration.

"delta_plot_step" : Parameter value is floating point number. It defines the step size in the plots for resolution functions.
  
"delta_plot_limit" : Parameter value is floating point number. It defines the maximal value of frequency in the plots for resolution functions.
    
"center_start" : Parameter value is floating point number. It defines the minimal value of the center of resolution function.
  
"center_stop" : Parameter value is floating point number. It defines the maximal value of the center of resolution function.
  
"center_delta" : Parameter value is floating point number. It defines the step size for the center of resolution functions.
  
"lambda" : With the help of this parameter it's possible to tune the regularization procedure.
  General format of the second string:
```
    integer_number     floating_point_parameter
```
  Possible values of the integer number (it defines the type of regularization):  
  *	0 : without regularization (thus no floating point parameter)
  *    +-1 : regularization by addition of covariance matrix:  (1-\lambda) S_{ij}.
  *    +-2 : regularization by neglecting all eigenvalues of W-kernel  which are less than \lambda (Truncated SVD decomposition).
  *    +-3 : regularization by Tikhonov filtering in SVD decomposition.
  *    +-4 : alternative variant of Tikhonov filtering in SVD  decomposition (more smooth).
  In all cases if the integer parameter is positive then the floating point parameter=\lambda. If it's negaive then the floating point parameter=average relative error
  and  lambda is choosen automatically.


"flag_exclude_corr" : Defines the additional regularization by neglecting some points in correlator or by averaging over some intervals in Euclidean time.
  Format of the second string: 
```
    integer_number   additional_parameters....
```
  Possible values of the integer number:
  *	0 : no additional regularization (thus no further parameters)
  *	1 : additional regularization is introduced: all points except the listed ones are neglected. 
  Thus the  additional parameters in the second string define the list of the points which are taken into account.
```
	    1   N    i_1   i_2  .....  i_N 
```
  N is the number of valid timeslices.
  Numbers should be sorted in ascending order, contact term in correlator is under number 0, so valid timeslices can start from 1 and end at i_N=Nt.
  *	2 : the correlator is averaged over intervals between listed points.
  The additional parameters define these intervals:
```
	    2   N   start_1   stop_1  .....  start_N  stop_N
```
  N is the number of intervals. Averaging is performed in each  inteval  over points t=start_i...stop_i  (including the first and the last point).


### Example of the file with parameters

```
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
-3  0.1
flag_exclude_corr
1  10  1 2 5 10 15 20 25 30 35 40
```


### Default values of parameters.
```
regularization type:				-3
relative error:   				0.05
value of regularization constant (\lambda):	1.0e-06
kernel type					1
accuracy 					1.0e-12
number of steps in numerical integration 	1000000
delta_plot_step					0.25
delta_plot_limit				Nt *2.0
center_start 					0.0
center_delta 					1.0
center_stop					Nt
```

## Output

In the output all dimensional quantities are in the units of temperature.
Below we present the list of files with a brief description of their formats.

1. "correlator_control_pre_N.txt"

   Preliminary values of correlator for N-th bin (before averaging over intervals or neglecting some data points).

   Format:
```
	time  correlator   error
```
2. "correlator_control_intervals_N.txt"  and  "correlator_control_fin_N.txt" 

   Final values of correlator for N-th bin (after averaging over intervals or neglecting some data points, if this type of regularization was ordered in parameters file).
 
   Format:

```
	time(or interval_number)  correlator(or average correlator)   error
```

3. "cov_matrix_control_pre_N.txt" and "cov_matrix_control_fin.txt" 

   These files contain covariance matrixes for correlators from "correlator_control_pre_N.txt" and "correlator_control_fin_N.txt" (or "correlator_control_intervals_N.txt").

   Format:
```
	C_{11} ..........  C_{1 Nt}
        ..........................
	C_{Nt 1} .......  C_{Nt Nt}
```
   
4. "correlataion_study_timeN.txt"

   The data produced during the calculation of autocorrelation length for N-th time slice.
   
   Format:
```   
	binsize   bin_error/initial_error
```

   "bin_error" is the statistical error calculated after data binning, "initial_error" is the ordinary statistical error. This ratio, in the limit of large binsize, is equal to the
    autocorrelation length.
	 
5. "delta_function_c=_VALUE_T.txt"

   Resolution function with center at frequency = _VALUE_ (in units of temperature).

   Format:
```
	frequency(in units of T)  resolution_function
```

   Resolution functions are saved only for the final value of regularization parameter lambda (when it's tuned to get some fixed value of relative error).
   

6. "rho_final.txt"
   
   Final output for spectral function (at final value of regularization parameter lambda).
   
   Format:
```
	center_of_resolution_function   spectral_function    error
```
	
7. "rho_basic.txt"

   This file is the same as "rho_final.txt" but it includes more technical information about resolution functions (like width, etc.).
   In case of data blocking this file appears for each block of data separately for the final value of regularization parameter ambda
   
   Format:
```   
    center_of_resolution_function  resolution_width(as it was in the  calculation of D functional)  position_of_the_maximum_of_resolution_function  width=(start-stop)/2  spectral_function  error  start  stop
```


   Start and stop points are the half-max points before and after the maximum of resolution function.

8. "rho_lambda_VALUE_avg.txt"

   Spectral function at some value of regularization parameter lambda in the vicinity of the final one.

   Format:
```
	center_of_resolution_function   spectral_function    error
```

9. "rho_lambda.txt"

   Dependence of resolution function on regularization parameter lambda.

   Format:

```   
	value_of_lambda	rho_for_c=0.000e+00 T	rho_error  ...... (for all positions of the center of resolution functions)
```	
	
	

