#include  <stdio.h>
#include <stdlib.h>
#include  "basic_structures.h"

#define LOG_FILE_OPERATION(operation)   {if(flag_log_output){operation;}}   
#define SPECIAL_LOG_FILE_OPERATION(operation)   {if(special_flag_log_output||flag_log_output){operation;}}   

//open file in directory for logs
FILE* fopen_control(const char* name, const char* aim);

FILE* fopen_log(const char* name, const char* aim, double parameter);

//print all parameters in file
bool print_parameters(FILE* file_out, correlator* pC);

//input of correlator and covariance matrix 
//and simultaneous conversion from initial data to final data taking into account that we neglect some points or treat only averages over intervals
bool input_correlator_matrix(FILE* file_in_current, FILE* file_in_matrix, correlator* pC);
void input_raw_data(FILE* file_in_current);
void get_jack_sample(correlator *C_jack, int jack_sample);
