#include  <stdio.h>
#include <stdlib.h>
#include  "basic_structures.h"

//open file in directory for logs
FILE* fopen_control(const char* name, const char* aim);

FILE* fopen_log(const char* name, const char* aim, double parameter);

//print all parameters in file
bool print_parameters(FILE* file_out, correlator* pC);

//input of correlator and covariance matrix 
//and simultaneous conversion from initial data to final data taking into account that we neglect some points or treat only averages over intervals
bool input_correlator_matrix(FILE* file_in_current, FILE* file_in_matrix, correlator* pC);
