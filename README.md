# Solver for Green-Kubo relations

The programs performs analytical continuation from Euclidean to real time through the solution of Green-Kubo relation. 
The modified Backus-Gilbert method is empoyed to regularize this ill-posed problem. 
The solver can work with several kernels in Green-Kubo relations including kernel for Density of States (DOS) and conductivity (see below). 
Example calculation of Density of States for Square Hubbard model from Monte Carlo data is included.


## Compilation

The program requires GNU GSL library. The path to the library should be inserted in LIB and INCLUDE_DIR variables in Makefile.
Once the path is set up, compilation can be made just by typing

```
$ make all
```

## Command line options

The program can work in two distinct regimes: 1) errors of the spectral function are estimated through the data binning 
(in this regime it needs the full Monte Carlo ensemble of the data)
and 2) errors of the spectral function are estimated on the basis of average value of the Euclidean correlator and its errors (ideally, 
in this regime the covariance matrix of the correlator should be provided).

-a : If this option is present the program enters the 2nd regime. If it's absent, the 1st regime with data binning is activated.

### Common options for both regimes:

-t  Nt : Nt is positive integer number  equal to the half number of Euclidean timeslices in initial correlator. 
-o path : With this option the path where to place the output data can be defined.
-p filename.txt :  If one wants to tune the internal parameters of the algorithm, the path to 
the file which contains these parameters can be defined here.  Without this option all parameters are equal to their default values. 
The format of the file with parameters is decribed below.

### Additional options in the regime with data binning:

-b Nb : Nb is positive integer number, it sets the number of bins for the data binning procedure. If this number is not provided 
by the user then the number of bins is set up automatically on the basis of autocorrelation length measurement.
-c filename.txt : file with full ensemble of full-time Euclidean correlators. File format:

```
	0    		Real_Part   Imaginary_Part
	1    		....
	2    		....
	2 Nt-1 	....
	(empty line)
	......(data for the next configuration)
```
	
-i : This options turns off the imaginary part of correlator. In this case the input file is organized in two columns.

	depending on the kernel, performs symmetrization/antisymmetrization of the input data:
	if kernel_switcher=5(see below) - performs antisymmetrization
	else - symmetrization
	

### Additional options in the regime without data binning:

-c  filename.txt   - file with correlator  (it is assuned that the correlator is symmetrical/antisymmetrical with respect to the half of Euclidean time)
	format:
```
	1    correlator    correlator_error 
	2 ....
	3 ....
	4 ....
	Nt ....
```	

	The input of errors is optional: if -e option exists, the format is the following
	But in this case the error estimation and the automated procedures for choosing regularization don't work
	
    -m  filename.txt -file with covariance matrix   -   the input of covariance matrix is optional, can work without it, but then the regularization with covariance matrix doesn't work
	format:
```
	C_{11} .......  C_{1 Nt}
        ..........................
	C_{Nt 1} .......  C_{Nt Nt}
```
	(without any indexes)


