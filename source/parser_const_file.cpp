#include "constants.h"
#include "parser_const_file.h"

void parse_option(FILE* file, const char* pattern, void* ptr)
{
    SKIP_COMMENT_LINE(file)
     fscanf(file, pattern, ptr);
    SKIP_REMAINING_CHARS(file)
}

void parse_lambda_option(FILE* file_const)
{
    SKIP_COMMENT_LINE(file_const)
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
}

void parse_exclusion_option(FILE* file_const, correlator* pC)
{
    SKIP_COMMENT_LINE(file_const)
     fscanf(file_const, "%d", &flag_model);
    if(!flag_model)
    {
    
    
printf("Nt_2=%d\n",Nt_2);
    
	int i;
	pC->format(Nt_2, Nt_2);
	for(i=0;i<pC->N_valid_points; i++)
	    pC->points_numbers[i]=i+1;
	SKIP_REMAINING_CHARS(file_const)
    }
    else
    {
	int i;
	int n_valid;
	fscanf(file_const, "%d", &n_valid);
	pC->format(Nt_2, n_valid);
	for(i=0;i<pC->N_valid_points; i++)
	    fscanf(file_const, "%d", &(pC->points_numbers[i]));
	SKIP_REMAINING_CHARS(file_const)
    }
}

void parse_force_zero_option(FILE* file_const)
{
    SKIP_COMMENT_LINE(file_const)
     fscanf(file_const, "%d",&flag_exclude_delta);
    if(!flag_exclude_delta)
    {SKIP_REMAINING_CHARS(file_const)}
    else
    {
	fscanf(file_const,"%d",&count_start_exclude);
	SKIP_REMAINING_CHARS(file_const)
    }
}

bool parse_const_file(FILE* file_const, correlator* pC)
{
    parse_option(file_const,"%d",&kernel_switcher);
    parse_option(file_const,"%le",&accuracy);
    parse_option(file_const,"%ld",&N_int_steps);
    parse_option(file_const,"%le",&omega_plot_delta);
    parse_option(file_const,"%le",&omega_plot_limit);
    parse_option(file_const,"%le",&center_start);
    parse_option(file_const,"%le",&center_stop);
    parse_option(file_const,"%le",&center_delta);

    parse_lambda_option(file_const);

    parse_exclusion_option(file_const, pC);

    parse_force_zero_option(file_const);

//setup intervals
    pC->construct_intervals();

    return true;
}


