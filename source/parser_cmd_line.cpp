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
//output directory name
    if(option_exists (argv, argv+argc, "-o"))
    {
        sprintf(output_directory,"%s", option_parameter(argv, argv + argc, "-o"));
    }
    else
    {
    	sprintf(output_directory,"%s",default_directory);
    }

//correlator file  name
    if(option_exists (argv, argv+argc, "-c"))
    {
        sprintf(correlator_filename,"%s",option_parameter(argv, argv + argc, "-c"));
    }
    else
    {
	return false;
    }

//covariance matrix file  name
    if(option_exists (argv, argv+argc, "-m"))
    {
        sprintf(cov_matrix_filename,"%s",option_parameter(argv, argv + argc, "-m"));
    }
    else
    {
	return false;
    }
    

//file with parameters
    if(option_exists (argv, argv+argc, "-p"))
    {
        sprintf(parameters_filename,"%s", option_parameter(argv, argv + argc, "-p"));
    }
    else
    {
	return false;
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

    //jackknife samples
    if(option_exists(argv, argv+argc, "-a")) {
      flag_jackknife=0;
      num_jack_samples=1;
    }
    else if(option_exists(argv, argv+argc, "-b")) {
      flag_jackknife=1;
      char *num_jack_samples_string = option_parameter(argv, argv+argc, "-t");
      if(!num_jack_samples_string) {
	return false;
      }
      else {
	num_jack_samples=atoi(num_jack_samples_string);
      }
    }
    else {
      return false;
    }

    return true;
} 

