#include "constants.h"
#include "parser_cmd_line.h"

bool option_exists (char** start, char** stop, const std::string& option)
{
    return std::find(start, stop, option)!=stop;
}

char* option_parameter (char** start, char** stop, const std::string& option)
{
    char ** parameter = std::find(start, stop, option);
    if (parameter != stop && ++parameter != stop)
    {
        return *parameter;
    }
    return 0;
}

bool parse_cmd_line(const int& argc, char ** argv)
{


  //jackknife samples
  if(option_exists(argv, argv+argc, "-a")) 
  {
    flag_jackknife=0;
    num_jack_samples=1;
    flag_tune_blocking=false;
  }
  else if(option_exists(argv, argv+argc, "-b")) 
  {
    flag_jackknife=1;
    char *num_jack_samples_string = option_parameter(argv, argv+argc, "-b");
    if(!num_jack_samples_string) 
    {
      return false;
    }
    else 
    {
      num_jack_samples=atoi(num_jack_samples_string);
      
      if(num_jack_samples<2) 
      {
	printf("Invalid number of jack knife samples!\n");
	return false;
      }
    }
  }
  else 
  {
    flag_jackknife=1;
    flag_tune_blocking=true;
  }
  
//output directory name
    if(option_exists (argv, argv+argc, "-o"))
    {
        sprintf(output_directory,"%s", option_parameter(argv, argv + argc, "-o"));
    }
    else
    {
    	sprintf(output_directory,"%s",default_directory);
    }

//correlator/raw_data file  name
    if(option_exists (argv, argv+argc, "-c"))
    {
        sprintf(correlator_filename,"%s",option_parameter(argv, argv + argc, "-c"));
    }
    else
    {
	return false;
    }

//existence of errors in correlator
    if(option_exists (argv, argv+argc, "-e"))
    {
	flag_error_corr_input=false;
    }
    else
    {
	flag_error_corr_input=true;
    }


//existence imaginary part data in raw data input
    if(option_exists (argv, argv+argc, "-i"))
    {
	flag_imag_part_input=false;
    }
    else
    {
	flag_imag_part_input=true;
    }
    

    if(!flag_jackknife) {
      //covariance matrix file  name
      if(option_exists (argv, argv+argc, "-m")) {
	sprintf(cov_matrix_filename,"%s",option_parameter(argv, argv + argc, "-m"));
	flag_covariance_matrix_input=true;
      }
      else {
	flag_covariance_matrix_input=false;
      }
    }
    

//file with parameters
    if(option_exists (argv, argv+argc, "-p"))
    {
        sprintf(parameters_filename,"%s", option_parameter(argv, argv + argc, "-p"));
        flag_constants_file=true;
    }
    else
    {
	flag_constants_file=false;
    }
    
//number of correlator points
    if(option_exists (argv, argv+argc, "-t"))
    {
	char* Nt_string  = option_parameter(argv, argv + argc, "-t");
        if(!Nt_string)
	{
	    return false;
	}
	else
	{
	    Nt_2=atoi(Nt_string);
	}
    }
    else
    {
	return false;
    }
    dNt_2=(double) Nt_2;

    return true;
} 

