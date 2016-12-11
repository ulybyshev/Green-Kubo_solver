#include "constants.h"
#include "parser_const_file.h"

void parse_option(FILE* file, const char* pattern, const char* option, void* ptr)
{
    int i;
    char *description=(char*)calloc(1000,sizeof(char));
    char in;
    int count=0;
    
    while(1) {
      in=getc(file);
      if(in!='\n') {
	description[count]=in;
	count++;
      }
      else
	break;
    }
    if(strcmp(description,option)==0) {
      fscanf(file, pattern, ptr);
      SKIP_REMAINING_CHARS(file)
    }
    else {
      ungetc('\n',file);
      for(i=count-1;i>=0;i--)
	ungetc(description[i],file);
    }
    free(description);
}

void parse_lambda_option(FILE* file_const)
{
    int i;
    char *description=(char*)calloc(1000,sizeof(char));
    char in;
    int count=0;
  
    while(1) {
      in=getc(file_const);
      if(in!='\n') {
	description[count]=in;
	count++;
      }
      else
	break;
    }
    if(strcmp(description,"lambda")==0) {
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
    else {
      ungetc('\n',file_const);
      for(i=count-1;i>=0;i--)
	ungetc(description[i],file_const);
    }
    free(description);
}

void parse_exclusion_option(FILE* file_const, correlator* pC)
{
    int i;
    char *description=(char*)calloc(1000,sizeof(char));
    char in;
    int count=0;

    while(1) {
      in=getc(file_const);
      if(in!='\n') {
	description[count]=in;
	count++;
      }
      else
	break;
    }
    if(strcmp(description,"flag_exclude_corr")==0) {
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
	  for(i=0;i<pC->N_valid_points; i++)
	    fscanf(file_const, "%d", &(pC->points_numbers[i]));
	  SKIP_REMAINING_CHARS(file_const)
	    }
    }
    else {
      ungetc('\n',file_const);
      for(i=count-1;i>=0;i--)
	ungetc(description[i],file_const);
    }
    free(description);
}

void parse_force_zero_option(FILE* file_const)
{
    int i;
    char *description=(char*)calloc(1000,sizeof(char));
    char in;
    int count=0;

    while(1) {
      in=getc(file_const);
      if(in!='\n') {
	description[count]=in;
	count++;
      }
      else
	break;
    }
    if(strcmp(description,"flag_force_zero")==0) {   
      fscanf(file_const, "%d",&flag_exclude_delta);
      if(!flag_exclude_delta)
	{SKIP_REMAINING_CHARS(file_const)}
      else
	{
	  fscanf(file_const,"%d",&count_start_exclude);
	  SKIP_REMAINING_CHARS(file_const)
	    }
    }
    else {
      ungetc('\n',file_const);
      for(i=count-1;i>=0;i--)
	ungetc(description[i],file_const);
    }
    free(description);
}

bool parse_const_file(FILE* file_const, correlator* pC)
{
    parse_option(file_const,"%d","int kernel switcher",&kernel_switcher);
    parse_option(file_const,"%le","accuracy",&accuracy);
    parse_option(file_const,"%ld","N_int_steps",&N_int_steps);
    parse_option(file_const,"%le","omega_plot_delta",&omega_plot_delta);
    parse_option(file_const,"%le","omega_plot_limit",&omega_plot_limit);
    parse_option(file_const,"%le","center_start",&center_start);
    parse_option(file_const,"%le","center_stop",&center_stop);
    parse_option(file_const,"%le","center_delta",&center_delta);

    parse_lambda_option(file_const);

    parse_exclusion_option(file_const, pC);

    parse_force_zero_option(file_const);

//setup intervals
    pC->construct_intervals();

    return true;
}

void set_default_values() {

  flag_lambda_regularization=-4;
  relative_error=0.1;
  lambda=1.0e-06;
  flag_model=0;
  flag_exclude_delta=0;
  count_start_exclude=0;
  kernel_switcher=1;
  accuracy=1.0e-12;
  N_int_steps=1000000;
  omega_plot_delta=1.0;
  omega_plot_limit=Nt_2;
  center_start=0.0;
  center_stop=Nt_2/2;
  center_delta=1.0;
  flag_exclude_delta=0;
    
}

