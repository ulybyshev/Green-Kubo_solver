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

### Additional options in the regime with data binning
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
	
	-i option turns off imaginary part. In this case input is organized in two columns

	depending on the kernel, performs symmetrization/antisymmetrization of the input data:
	if kernel_switcher=5(see below) - performs antisymmetrization
	else - symmetrization
	


with -a option: enters the regime without data blocking, which works with the average correlator and covariance matrix:
    -c  filename.txt   - file with correlator  (it is assuned that the correlator is symmetrical/antisymmetrical with respect to the half of Euclidean time)
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





```
$ git clone https://github.com/Anny-Moon/Simple_PDB_parser
$ cd Simple_PDB_parser/
$ sh compile_script
$ ./pdb-reader
```
## Example: 5dn7
```
$ ./pdb-reader 5dn7
```
The output is:
```
Protein: 5dn7
The first CA atom has number 345.
Missing atoms from 361 to 365 (5 atoms).
Missing atoms from 558 to 562 (5 atoms).
The last CA atom has number  594.
Number of CA atoms in the model: 250.
But there is data only for 240 of them.

***************************************
*      Maps of missing atoms          *
*  atom: .         missing atom: 0    *
***************************************
Percentage: string length = 100 chars.

......00.............................................................................00.............

Actual: string length = number of atoms in model.

................00000................................................................................................................................................................................................00000................................

*****************************************
If you want to rewrite dat-file with only one segment, call
./pdb-reader (without arguments) for instructions.
```

The program created file `results/xyz_5dn7.dat` with 240 lines.

Let's concider the maps. They represent the positions and amounts of missing atoms, where `0` stands for missing atom and `.` for atom from the pdb-file.
The first map depicture the percentege picture (so not very precise becouse of rounding). The second one is
an actual picture: each char is one char.

Now I want to rewrite the dat-file with one segment in oder to not have missing data there. There are 3 segments
there. I will take the largest one, between two missing parts. I can do this by calling the program again but
now give the number of the first atom (or any missing atom before the first) in segment as the **second**
argument.

Like this:
```
$ ./pdb-reader 5dn7 365
```
or the same:
```
$ ./pdb-reader 5dn7 366
```
And the output is:
```
Protein: 5dn7
Number of CA atoms: 192.
```
what means that now in `results/xyz_5dn7.dat` we have 192 lines wich correspond to 192 atoms between the 
firs and the second missing parts.

If you want to take the first segment (from the first atom to the first missing atom) you always can put
second argument is egual to `0`:

```
$ ./pdb-reader 5dn7 0
```
### Note
If you want to have your outpun in file, not on the screen you always can do this:
```
$ ./pdb-reader 5dn7 > create_and_put_everything_here.dat
```

## Bugs
Since the autor doesn't know all possible formats of notation in PDB, there should be cases when the program will give absurd result. These cases will be fixed when revealed.

___
Anna Sinelnikova

Uppsala, 2017