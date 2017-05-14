#include "constants.h"
#include "parser_const_file.h"

void get_description(FILE* file, char* description)
{
    char in;
    int count=0;
    while(!feof(file))
    {
	in=getc(file);
        if(in!='\n')
        {
	    description[count]=in;
	    count++;
        }
        else
	    break;
    }
}


void parse_option(FILE* file, const char* pattern, const char* option, void* ptr)
{
    char *description;
    while(!feof(file))
    {
        description=(char*)calloc(1000,sizeof(char));
	get_description(file, description); 
	if(strcmp(description,option)==0)
	{
    	    fscanf(file, pattern, ptr);
    	    free(description);
    	    break;
	}
	free(description);
    }
    fseek(file, 0, SEEK_SET);

}

void parse_lambda_option(FILE* file_const, const char* option)
{
    char *description;
  
    while(!feof(file_const))
    {
        description=(char*)calloc(1000,sizeof(char));
	get_description(file_const, description); 
	if(strcmp(description,option)==0)
	{
	    fscanf(file_const, "%d", &flag_lambda_regularization);
    	    if(!flag_lambda_regularization)
	    {SKIP_REMAINING_CHARS(file_const)}
    	    else
	    {
		if(flag_lambda_regularization>0)
		    fscanf(file_const, "%le", &lambda);
		else
		    fscanf(file_const, "%le", &relative_error);
		SKIP_REMAINING_CHARS(file_const)
    	    }
	    
    	    free(description);
    	    break;
	}
	free(description);
    }
    fseek(file_const, 0, SEEK_SET);

}

void parse_exclusion_option(FILE* file_const, const char* option, correlator* pC)
{
    int i;
    char *description;
    bool is_option=false;

    while(!feof(file_const))
    {
        description=(char*)calloc(1000,sizeof(char));
	get_description(file_const, description); 
	if(strcmp(description,option)==0)
	{
	    is_option=true;
	    
    	    fscanf(file_const, "%d", &flag_model);
    	    if(!flag_model)
	    {
	  	  
		printf("Nt_2=%d\n",Nt_2);
	  	pC->format(Nt_2, Nt_2);
		for(i=0;i<pC->N_valid_points; i++)
		    pC->points_numbers[i]=i+1;
		SKIP_REMAINING_CHARS(file_const)
	    }
    	    else
	    {
		int n_valid;
		fscanf(file_const, "%d", &n_valid);
		pC->format(Nt_2, n_valid);

		if(flag_model==1)
		{
			for(i=0;i<n_valid; i++)
			{
			    fscanf(file_const, "%d", &(pC->points_numbers[i]));
			    pC->points_stop_numbers[i]=pC->points_numbers[i]; 
			}
		    SKIP_REMAINING_CHARS(file_const)
		}
		else
		{
			for(i=0;i<n_valid; i++)
			{
			    fscanf(file_const, "%d%d", &(pC->points_numbers[i]), &(pC->points_stop_numbers[i]));
			}
		}
	    }
        
    	    free(description);
    	    break;
	}
	free(description);
    }
    if(!is_option)
    {
	pC->format(Nt_2, Nt_2);
	for(i=0;i<pC->N_valid_points; i++)
	{
	    pC->points_numbers[i]=i+1;
	    pC->points_stop_numbers[i]=pC->points_numbers[i]; 
	}
    }

    fseek(file_const, 0, SEEK_SET);

}

void parse_force_zero_option(FILE* file_const, const char* option )
{
    char *description;
    bool is_option=false;
    
    while(!feof(file_const))
    {
        description=(char*)calloc(1000,sizeof(char));
	get_description(file_const, description); 
	if(strcmp(description,option)==0)
	{
	    is_option=true;
    	    fscanf(file_const, "%d",&flag_exclude_delta);
    	    if(!flag_exclude_delta)
		{SKIP_REMAINING_CHARS(file_const)}
    	    else
	    {
		fscanf(file_const,"%d",&count_start_exclude);
		SKIP_REMAINING_CHARS(file_const)
	    }
	    free(description);
	    break;
	}
	free(description);
    }
    if(!is_option)
    {
	flag_exclude_delta=0;
    }
    fseek(file_const, 0, SEEK_SET);

}

bool parse_const_file(FILE* file_const, correlator* pC)
{
    parse_option(file_const,"%d",KERNEL_SWITCHER_OPTION,&kernel_switcher);
    parse_option(file_const,"%le",ACCURACY_OPTION ,&accuracy);
    parse_option(file_const,"%ld",INTEGRAL_NUMBER_STEPS_OPTION ,&N_int_steps);
    parse_option(file_const,"%le",RES_FUNCTION_PLOT_STEP_OPTION ,&omega_plot_delta);
    parse_option(file_const,"%le",RES_FUNCTION_PLOT_LIMIT_OPTION ,&omega_plot_limit);
    parse_option(file_const,"%le",RES_FUNCTION_CENTER_START_OPTION ,&center_start);
    parse_option(file_const,"%le",RES_FUNCTION_CENTER_STOP_OPTION ,&center_stop);
    parse_option(file_const,"%le",RES_FUNCTION_CENTER_DELTA_OPTION ,&center_delta);
    parse_lambda_option(file_const, REGULARIZATION_OPTION);
    parse_exclusion_option(file_const, EXCLUDE_CORR_POINTS_OPTION , pC);
    parse_force_zero_option(file_const, RES_FUNCTION_ZERO_AT_ZERO_OMEGA_OPTION );

//setup intervals
    pC->construct_intervals();

    return true;
}

void set_default_values()
{
  flag_lambda_regularization=DEFAULT_FLAG_REGULARIZATION;
  relative_error=DEFAULT_RELATIVE_ERROR;
  lambda=DEFAULT_LAMBDA;
  flag_model=DEFAULT_FLAG_MODEL;
  flag_exclude_delta=DEFAULT_FLAG_EXCLUDE_DELTA;
  count_start_exclude=DEFAULT_COUNT_START_EXCLUDE;
  kernel_switcher=DEFAULT_KERNEL_SWITCHER;
  accuracy=DEFAULT_ACCURACY;
  N_int_steps=DEFAULT_N_INT_STEPS ;
  omega_plot_delta=DEFAULT_OMEGA_PLOT_DELTA;
  omega_plot_limit= (double)Nt_2 *2.0;
  center_start=DEFAULT_CENTER_START;
  center_stop=(double)Nt_2;
  center_delta=DEFAULT_CENTER_DELTA;
    
}

