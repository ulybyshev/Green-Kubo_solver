#include "constants.h"
#include "input_output.h"

FILE* fopen_control(const char* name, const char* aim)
{
  FILE* file_out;
  char* final_name=(char*) malloc(sizeof(output_directory)+sizeof(name)+10);
  sprintf(final_name, "%s/%s", output_directory, name);
  file_out=fopen(final_name, aim);
  free(final_name);
  return file_out;
}


FILE* fopen_log(const char* name, const char* aim, double parameter)
{
  FILE* file_out;
  char final_name[1000];
  sprintf(final_name, "%s/%s", output_directory, name);
  file_out=fopen(final_name, aim);
  fprintf(file_out,"\n\n\n******************************\n");
  fprintf(file_out,"center/T=%.15le\n", parameter*2.0*dNt_2);
  fprintf(file_out,"******************************\n");

  fflush(file_out);
  return file_out;
}

bool print_parameters(FILE* file_out, correlator* pC)
{
    fprintf(file_out,"***********************Analytical continuation**********************\n");

    fprintf(file_out,"number of points in symmetrized correlator:%d\n", pC->N_full_points);

    fprintf(file_out,"file with current correlator:\n %s\n", correlator_filename);

    fprintf(file_out,"file with covariance matrix:\n %s\n", cov_matrix_filename);

    fprintf(file_out,"path for output:\n %s\n", output_directory);

    fprintf(file_out,"file with parameters:\n %s\n", parameters_filename);

    fprintf(file_out,"kernel switcher=: %d\n", kernel_switcher);

    fprintf(file_out,"accuracy=%.15le\n",accuracy);
    fprintf(file_out,"Number os steps in numerical interation=%ld\n",N_int_steps);
    fprintf(file_out,"Step in plots of resolution function(omega)=%.15le\n",omega_plot_delta);
    fprintf(file_out,"Limit in plots of resolution function(omega)=%.15le\n",omega_plot_limit);
    fprintf(file_out,"Center of resolution function starts from = %.15le\n",center_start);
    fprintf(file_out,"Center of resolution function stops at %.15le\n",center_stop);
    fprintf(file_out,"Step in center of resolution function = %.15le\n",center_delta);

    fprintf(file_out,"flag_lambda_regularization=%d\n",flag_lambda_regularization);
    if(flag_lambda_regularization>0)
	fprintf(file_out,"Lambda(regularization)=%.15le\n",lambda);
    else
	fprintf(file_out,"ordered relative error=%.15le\n",relative_error);

    fprintf(file_out,"Take into account all (0)  or not all (1) timeslices or average over intervals (2) = %d\n",flag_model);
    if(flag_model)
    {
	int i;
	fprintf(file_out,"Number of timeslices taken into account = %d\n",pC->N_valid_points);
	for(i=0;i<pC->N_valid_points;i++)
	{
	    fprintf(file_out,"%d\n",pC->points_numbers[i]);
	}
    }

    fprintf(file_out,"Force or not (1 or 0) resolution functions to be zero at zero frequency: %d\n",flag_exclude_delta);
    fprintf(file_out,"Number of resolution function  to start enforcement of condition f(0)=0 : %d\n",count_start_exclude);

    if(flag_model==2)
    {	
	int i,j;
	fprintf(file_out,"\n\nN_intervals %d\n",pC->N_valid_points);
	for(i=0;i<pC->N_valid_points;i++) 
	{
	    fprintf(file_out,"******************\n");
	    fprintf(file_out,"Interval %d\n",i);
	    fprintf(file_out,"******************\n");
	    for(j=0;j<pC->interval_numbers[i]->size;j++)
		fprintf(file_out,"%d\n",pC->interval_numbers[i]->times[j]);
	}
    }
    fflush(file_out);
    return true;
}



bool input_correlator_matrix(FILE* file_in_current, FILE* file_in_matrix, correlator* pC)
{
    int t,i, j,par_int;
    double par_double;
    FILE* file_out;


    for(t=1;t<=pC->N_full_points;t++)
    {
	fscanf(file_in_current, "%d", &par_int);
	fscanf(file_in_current, "%le", &(pC->corr_full[t-1]));
	fscanf(file_in_current, "%le", &(pC->error_full[t-1]));
    }

    file_out=fopen_control("correlator_control_pre.txt","w");
    for(t=0;t<pC->N_full_points;t++)
    {
	fprintf(file_out,"%d\t%.15le\t%.15le\n", t+1, pC->corr_full[t], pC->error_full[t]);fflush(file_out);
    }
    fclose(file_out);

    for(i=1;i<=pC->N_full_points;i++)
    for(t=1;t<=pC->N_full_points;t++)
     {
	fscanf(file_in_matrix, "%le", &par_double);
	gsl_matrix_set(pC->S_full, i-1, t-1, par_double);
    }

  file_out=fopen_control("cov_matrix_control_pre.txt","w");
  for(i=0;i<pC->N_full_points;i++){
  for(t=0;t<pC->N_full_points;t++)
  {
    fprintf(file_out,"%.15le\t", gsl_matrix_get(pC->S_full, i,t));fflush(file_out);
  }
  fprintf(file_out,"\n");
  }
  fclose(file_out);
  
  
  //conversion to real arrays taken into account only certain timeslices or intervals
///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////

if(flag_model==2)//intervals
{ 
  int count1, count2;
  double result=0.0;
  double temp_err=0.0;
  for(count1=0;count1<pC->N_valid_points;count1++)
  {
    for(i=0;i<pC->interval_numbers[count1]->size;i++) {
      result+=pC->corr_full[pC->interval_numbers[count1]->times[i]-1];
      temp_err+=pC->error_full[pC->interval_numbers[count1]->times[i]-1]*pC->error_full[pC->interval_numbers[count1]->times[i]-1];
    }
    pC->corr[count1]=result/((double)(pC->interval_numbers[count1]->size));
    pC->error[count1]=sqrt(temp_err)/((double)(pC->interval_numbers[count1]->size));
    result=0.0;
    temp_err=0.0; 
  }

  file_out=fopen_control("correlator_control_intervals.txt","w");
  for(t=0;t<pC->N_valid_points;t++)
  {
    fprintf(file_out,"%d\t%d\t%.15le\t%.15le\n" ,pC->points_numbers[t]-pC->interval_numbers[t]->size +1, pC->points_numbers[t], pC->corr[t], pC->error[t]);fflush(file_out);
  }
  fclose(file_out);

//TODO check covariance matrix calculation
  result=0.0;
  for(count1=0;count1<pC->N_valid_points;count1++)
  for(count2=0;count2<pC->N_valid_points;count2++)
  {
    for(i=0;i<pC->interval_numbers[count1]->size;i++) {
      for(j=0;j<pC->interval_numbers[count2]->size;j++) {
	par_double=gsl_matrix_get(pC->S_full,pC->interval_numbers[count1]->times[i]-1,pC->interval_numbers[count2]->times[j]-1);
	result+=par_double;
      }
    }
    par_double=result/((double)(pC->interval_numbers[count1]->size*pC->interval_numbers[count2]->size));
    gsl_matrix_set(pC->S,count1,count2,par_double);
    result=0.0;
  }

  file_out=fopen_control("cov_matrix_control.txt","w");
  for(i=0;i<pC->N_valid_points;i++){
  for(t=0;t<pC->N_valid_points;t++)
  {
    fprintf(file_out,"%.15le\t", gsl_matrix_get(pC->S, i,t));fflush(file_out);
  }
  fprintf(file_out,"\n");
  }
  fclose(file_out);
}
else//just neglecting points (the case when we save full correlator is also here)
{ 

  int count1, count2;

  for(count1=0;count1<pC->N_valid_points;count1++)
  {
	pC->corr[count1]=pC->corr_full[pC->points_numbers[count1]-1];
	pC->error[count1]=pC->error_full[pC->points_numbers[count1]-1];
  }


  file_out=fopen_control("correlator_control_fin.txt","w");
  for(t=0;t<pC->N_valid_points;t++)
  {
    fprintf(file_out,"%d\t%.15le\t%.15le\n", pC->points_numbers[t], pC->corr[t], pC->error[t]); fflush(file_out);
  }
  fclose(file_out);


  for(count1=0;count1<pC->N_valid_points;count1++)
  for(count2=0;count2<pC->N_valid_points;count2++)
  {
    par_double=gsl_matrix_get(pC->S_full,pC->points_numbers[count1]-1, pC->points_numbers[count2]-1);
    gsl_matrix_set(pC->S, count1, count2, par_double);
  }

  file_out=fopen_control("cov_matrix_control.txt","w");
  for(i=0;i<pC->N_valid_points;i++){
  for(t=0;t<pC->N_valid_points;t++)
  {
    fprintf(file_out,"%.15le\t", gsl_matrix_get(pC->S, i,t));fflush(file_out);
  }
  fprintf(file_out,"\n");
  }
  fclose(file_out);
}
//end of conversion
///////////////////////////////////////////
    return true;
}




