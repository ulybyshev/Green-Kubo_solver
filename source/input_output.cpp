#include "constants.h"
#include "input_output.h"
#include "sample_sequence.h"

#define MAXT 1000
#define MAXCONF 2000
#define SKIP_LINE(file) { while(getc(file)!='\n'); }
#define SKIP_REMAINING_CHARS(file) SKIP_LINE(file)
#define SKIP_BARRIER(file) { SKIP_LINE(file); SKIP_LINE(file); }

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

    fprintf(file_out,"flag_data_blocking=%d\n",flag_jackknife);

    fprintf(file_out,"flag_tune_blocking=%d\n",flag_tune_blocking);

    fprintf(file_out,"number of points in symmetrized correlator:%d\n", pC->N_full_points);

    fprintf(file_out,"file with current correlator:\n %s\n", correlator_filename);

    if(flag_covariance_matrix_input)
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
    if(flag_model==1)
    {
	int i;
	fprintf(file_out,"Number of timeslices taken into account = %d\n",pC->N_valid_points);
	for(i=0;i<pC->N_valid_points;i++)
	{
	    fprintf(file_out,"%d\n",pC->points_numbers[i]);
	}
    }

    if(flag_model==2)
    {
	int i;
	fprintf(file_out,"Number of intervals taken into account = %d\n",pC->N_valid_points);
	for(i=0;i<pC->N_valid_points;i++)
	{
	    fprintf(file_out,"%d\t%d\n",pC->interval_numbers[i]->times[0], pC->interval_numbers[i]->times[pC->interval_numbers[i]->size-1]);
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
    int t,i, j;
    double par_double;
    FILE* file_out;

    for(t=1;t<=pC->N_full_points;t++)
    {
		fscanf(file_in_current, "%le", &par_double);
		fscanf(file_in_current, "%le", &(pC->corr_full[t-1]));
		if(flag_error_corr_input)
			fscanf(file_in_current, "%le", &(pC->error_full[t-1]));
		else
			pC->error_full[t-1]=0.0;
    }

    file_out=fopen_control("correlator_control_pre.txt","w");
    fprintf(file_out,"#time\t correlator_Re\t  error_Re\n");
    for(t=0;t<pC->N_full_points;t++)
    {
		fprintf(file_out,"%d\t%.15le\t%.15le\n", t+1, pC->corr_full[t], pC->error_full[t]);
		fflush(file_out);
    }
    fclose(file_out);

    for(i=1;i<=pC->N_full_points;i++)
    	for(t=1;t<=pC->N_full_points;t++)
    	{
			if(flag_covariance_matrix_input)
	    		fscanf(file_in_matrix, "%le", &par_double);
			else
	    		par_double=0.0;
			gsl_matrix_set(pC->S_full, i-1, t-1, par_double);
    	}

    if(flag_covariance_matrix_input)
    {
		file_out=fopen_control("cov_matrix_control_pre.txt","w");
		for(i=0;i<pC->N_full_points;i++)
		{
			for(t=0;t<pC->N_full_points;t++)
			{
	    		fprintf(file_out,"%.15le\t", gsl_matrix_get(pC->S_full, i,t));fflush(file_out);
			}
			fprintf(file_out,"\n");
		}
		fclose(file_out);
    }
  
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
  fprintf(file_out,"#interval_number\t average_correlator_Re\t  error_Re\n");
  for(t=0;t<pC->N_valid_points;t++)
  {
    fprintf(file_out,"%d\t%d\t%.15le\t%.15le\n" ,pC->interval_numbers[t]->times[0],pC->interval_numbers[t]->times[pC->interval_numbers[t]->size-1], pC->corr[t], pC->error[t]);fflush(file_out);
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

if(flag_covariance_matrix_input)
{
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
  fprintf(file_out,"#time\t correlator_Re\t  error_Re\n");
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
  
  if(flag_covariance_matrix_input)
  {	
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
}
//end of conversion
///////////////////////////////////////////
    return true;
}

int input_raw_data(FILE* file_in_current) 
{
	int t_count,conf_count,t;
	double re_g,im_g;
	int check;
  
	double *temp_g=(double *)calloc(MAXT*MAXCONF,sizeof(double));
	if(temp_g==NULL)
	{
    	fprintf(stderr,"error allocating memory for raw data!\n");
    	exit(1);
  	};
    
  	conf_count=0;
  	int Nt=2*Nt_2;
  	
  	int ts = 0, t0;
  	char tstr[512];
  	
  	while(fgets(tstr, 512, file_in_current)!=NULL && ts<MAXT*MAXCONF)
  	{
  		//printf("%s",tstr);
  		sscanf(tstr, "%d%lf", &t0, &re_g);
  		printf("\t%02i\t%+2.4E\n", t0, re_g);
  		temp_g[ts] = re_g;
  		ts++;
	};
	printf("\n");
	conf_count = ts/Nt;
	if(ts%Nt!=0)
	{
		fprintf(stderr, "The number of lines %i does not match the integer number of configurations with Nt=%i", ts, Nt);	fflush(stdout);
		fprintf(stdout, "The number of lines %i does not match the integer number of configurations with Nt=%i", ts, Nt);	fflush(stderr);
	};
	printf("\n\t >> \x1b[1;35m %i configurations in total \x1b[0m\n\n", conf_count);
  
  //set parameters
  n_conf=conf_count;
  dn_conf=(double)n_conf;
  //allocate memory for raw data
  raw_data=(double *)calloc(Nt_2*n_conf,sizeof(double));

  //fold data
  for(conf_count=0;conf_count<n_conf;conf_count++) {
    for(t_count=1;t_count<Nt_2;t_count++) {
    
	if(kernel_switcher!=5)
	{//symmetrical case
    	    raw_data[t_count-1+conf_count*Nt_2]=(temp_g[t_count+conf_count*Nt]+temp_g[(Nt-t_count)+conf_count*Nt])/2;
    	}
    	else
    	{//antisymmetrical case
    	    raw_data[t_count-1+conf_count*Nt_2]=(temp_g[t_count+conf_count*Nt]-temp_g[(Nt-t_count)+conf_count*Nt])/2;
    	}
      
    }
    if(kernel_switcher!=5)
	raw_data[Nt_2-1+conf_count*Nt_2]=temp_g[Nt_2+conf_count*Nt];
    else
	raw_data[Nt_2-1+conf_count*Nt_2]=0.0;
    
    
  }


 //output of full average and error 
{ 
  int i,j;
  double *avg=(double *)calloc(Nt_2,sizeof(double));
  double *err=(double *)calloc(Nt_2,sizeof(double));
  for(i=0;i<Nt_2;i++) {
    avg[i]=0.0;
    err[i]=0.0;
  }
  
   for(i=1;i<=Nt_2;i++) {
    for(j=0;j<=n_conf;j++) {
      avg[i-1]+=raw_data[(i-1)+(j-1)*Nt_2];
    }
    avg[i-1]=avg[i-1]/((double)n_conf);
  }
  
  for(i=1;i<=Nt_2;i++) {
    for(j=0;j<=n_conf;j++) {
      err[i-1]+=(raw_data[(i-1)+(j-1)*Nt_2]-avg[i-1])*(raw_data[(i-1)+(j-1)*Nt_2]-avg[i-1]);
    }
    err[i-1]=sqrt(err[i-1]/((double)n_conf*(n_conf-1.0)));
  }
  FILE* file_out=fopen_control("correlator_control_full.txt","w");
  fprintf(file_out,"#time\t correlator_Re\t  error_Re\n");

  for(i=0;i<=Nt_2-1;i++) {
    fprintf(file_out,"%d\t%.15le\t%.15le\n", i+1, avg[i], err[i]);fflush(file_out);
  }
  fclose(file_out);
  free(avg);
  free(err);
  
  printf("Bye from input_raw_data!!!\n"); fflush(stdout);
}

  
  free(temp_g);
  return n_conf;
}

void get_jack_sample(correlator *C_jack, int jack_sample) {

  int i,j,k,first_conf,last_conf, t;
  int num_configs=n_conf/num_jack_samples;
  FILE *file_out;
  char file_name[1000];
  
  double *avg=(double *)calloc(Nt_2,sizeof(double));
  double *err=(double *)calloc(Nt_2,sizeof(double));
  double *cov=(double *)calloc(Nt_2*Nt_2,sizeof(double));
  
  for(i=0;i<Nt_2;i++) {
    avg[i]=0.0;
    err[i]=0.0;
  }
  
  if(jack_sample==num_jack_samples) {
    first_conf=(num_jack_samples-1)*num_configs+1;
    last_conf=n_conf;
    num_configs=last_conf-first_conf+1;
  }
  else {
    first_conf=(jack_sample-1)*num_configs+1;
    last_conf=jack_sample*num_configs;
  }
  
  //construct average and error for jackknife sample
  for(i=1;i<=Nt_2;i++) {
    for(j=first_conf;j<=last_conf;j++) {
      avg[i-1]+=raw_data[(i-1)+(j-1)*Nt_2];
    }
    avg[i-1]=avg[i-1]/((double)num_configs);
    C_jack->corr_full[i-1]=avg[i-1];
  }

  for(i=1;i<=Nt_2;i++) {
    for(j=first_conf;j<=last_conf;j++) {
      err[i-1]+=(raw_data[(i-1)+(j-1)*Nt_2]-avg[i-1])*(raw_data[(i-1)+(j-1)*Nt_2]-avg[i-1]);
      
    }
    err[i-1]=sqrt(err[i-1]/((double)num_configs*(num_configs-1.0)));
    C_jack->error_full[i-1]=err[i-1];
  }

  //output jack sample average and error
  sprintf(file_name,"correlator_control_pre_%d.txt",jack_sample);
  file_out=fopen_control(file_name,"w");
  fprintf(file_out,"#time\t correlator_Re\t  error_Re\n");

  for(t=0;t<C_jack->N_full_points;t++) {
    fprintf(file_out,"%d\t%.15le\t%.15le\n", t+1, C_jack->corr_full[t], C_jack->error_full[t]);fflush(file_out);
  }
  fclose(file_out);
  
  //construct covariance matrix
  for(i=1;i<=Nt_2;i++)
    for(j=1;j<=Nt_2;j++)
    {
      for(k=first_conf;k<=last_conf;k++) {
	
	cov[(i-1)+(j-1)*Nt_2]+=raw_data[(i-1)+(k-1)*Nt_2]*raw_data[(j-1)+(k-1)*Nt_2];
      }
      
      cov[(i-1)+(j-1)*Nt_2]=cov[(i-1)+(j-1)*Nt_2]/(double)num_configs - avg[i-1]*avg[j-1];
      gsl_matrix_set(C_jack->S_full, i-1, j-1, cov[(i-1)+(j-1)*Nt_2]);
    }

  //output jack sample covariance matrix
  sprintf(file_name,"cov_matrix_control_pre_%d.txt",jack_sample);
  file_out=fopen_control(file_name,"w");
  for(i=0;i<C_jack->N_full_points;i++){
    for(t=0;t<C_jack->N_full_points;t++) {
      fprintf(file_out,"%.15le\t", gsl_matrix_get(C_jack->S_full, i,t));fflush(file_out);
    }
    fprintf(file_out,"\n");
  }
  fclose(file_out);
  
  free(avg);
  free(err);
  free(cov);

  //compute for cases of intervals or just neglecting points (TBD)
  if(flag_model==2) {

    int count1, count2;
    double result=0.0;
    double temp_err=0.0;
    double par_double;
    for(count1=0;count1<C_jack->N_valid_points;count1++) {
      for(i=0;i<C_jack->interval_numbers[count1]->size;i++) {
	result+=C_jack->corr_full[C_jack->interval_numbers[count1]->times[i]-1];
	temp_err+=C_jack->error_full[C_jack->interval_numbers[count1]->times[i]-1]*C_jack->error_full[C_jack->interval_numbers[count1]->times[i]-1];
      }
      C_jack->corr[count1]=result/((double)(C_jack->interval_numbers[count1]->size));
      C_jack->error[count1]=sqrt(temp_err)/((double)(C_jack->interval_numbers[count1]->size));
      result=0.0;
      temp_err=0.0; 
    }

    sprintf(file_name,"correlator_control_intervals_%d.txt",jack_sample);
    file_out=fopen_control(file_name,"w");
    fprintf(file_out,"#interval_number\t average_correlator_Re\t  error_Re\n");

    for(t=0;t<C_jack->N_valid_points;t++) {
      fprintf(file_out,"%d\t%d\t%.15le\t%.15le\n",C_jack->interval_numbers[t]->times[0] ,C_jack->interval_numbers[t]->times[C_jack->interval_numbers[t]->size-1], C_jack->corr[t], C_jack->error[t]);fflush(file_out);
    }
    fclose(file_out);
    
    //TODO check covariance matrix calculation
    result=0.0;
    for(count1=0;count1<C_jack->N_valid_points;count1++)
      for(count2=0;count2<C_jack->N_valid_points;count2++) {
	for(i=0;i<C_jack->interval_numbers[count1]->size;i++) {
	  for(j=0;j<C_jack->interval_numbers[count2]->size;j++) {
	    par_double=gsl_matrix_get(C_jack->S_full,C_jack->interval_numbers[count1]->times[i]-1,C_jack->interval_numbers[count2]->times[j]-1);
	    result+=par_double;
	  }
	}
	par_double=result/((double)(C_jack->interval_numbers[count1]->size*C_jack->interval_numbers[count2]->size));
	gsl_matrix_set(C_jack->S,count1,count2,par_double);
	result=0.0;
      }

    sprintf(file_name,"cov_matrix_control_intervals_%d.txt",jack_sample);
    file_out=fopen_control(file_name,"w");
    for(i=0;i<C_jack->N_valid_points;i++){
      for(t=0;t<C_jack->N_valid_points;t++) {
	fprintf(file_out,"%.15le\t", gsl_matrix_get(C_jack->S, i,t));fflush(file_out);
      }
      fprintf(file_out,"\n");
    }
    fclose(file_out);
    
  }
  else { //just neglecting points (the case when we save full correlator is also here)  

    int count1, count2;
    double par_double;
    
    for(count1=0;count1<C_jack->N_valid_points;count1++) {
      C_jack->corr[count1]=C_jack->corr_full[C_jack->points_numbers[count1]-1];
      C_jack->error[count1]=C_jack->error_full[C_jack->points_numbers[count1]-1];
    }

    sprintf(file_name,"correlator_control_fin_%d.txt",jack_sample);
    file_out=fopen_control(file_name,"w");
    fprintf(file_out,"#time\t correlator_Re\t  error_Re\n");
    for(t=0;t<C_jack->N_valid_points;t++) {
      fprintf(file_out,"%d\t%.15le\t%.15le\n", C_jack->points_numbers[t], C_jack->corr[t], C_jack->error[t]); fflush(file_out);
    }
    fclose(file_out);
    
    
    for(count1=0;count1<C_jack->N_valid_points;count1++)
      for(count2=0;count2<C_jack->N_valid_points;count2++) {
	par_double=gsl_matrix_get(C_jack->S_full,C_jack->points_numbers[count1]-1, C_jack->points_numbers[count2]-1);
	gsl_matrix_set(C_jack->S, count1, count2, par_double);
      }

    sprintf(file_name,"cov_matrix_control_fin_%d.txt",jack_sample);
    file_out=fopen_control(file_name,"w");
    for(i=0;i<C_jack->N_valid_points;i++){
      for(t=0;t<C_jack->N_valid_points;t++) {
	fprintf(file_out,"%.15le\t", gsl_matrix_get(C_jack->S, i,t));fflush(file_out);
      }
      fprintf(file_out,"\n");
    }
    fclose(file_out);
  }
  //end of conversion
  ///////////////////////////////////////////
  
    
}

//defines automatically binsize for blocking on the basis of autocorrelation length calculation
bool input_data_analysis(initial_data_description* D)
{
    int N_histories, delta, delta_output, count_history;
    bool result=true;

    N_histories=Nt_2;
    delta=1;

    if(Nt_2<5)
    {
	delta_output=1;
    }
    else
    {
	delta_output=Nt_2/5;
    }


    
    D->format(N_histories);
    for(count_history=0; count_history<N_histories; count_history++)
    {
	D->times[count_history]=delta*count_history+1;
    }
	
    sample_sequence<double>* corr_history;
    corr_history=new sample_sequence<double>[N_histories];
    
    //formatting the arrays
    
    for(count_history=0; count_history<N_histories; count_history++)
    {
	corr_history[count_history].format(n_conf);
    }
    
    //filling the arrays and computation statistical characteristics
    for(count_history=0; count_history<N_histories; count_history++)
    {
	int count_conf, cur_time;
	cur_time=D->times[count_history]-1;
	for(count_conf=0;count_conf<n_conf; count_conf++)
	{
	    corr_history[count_history].data[count_conf]=raw_data[cur_time+count_conf*Nt_2];
	}	
	
	FILE* file_output;
	char filename[1024];
	
	int corr_res;
	if(cur_time%delta_output==0)
	{
	    sprintf(filename, "correlation_study_time%d.txt", D->times[count_history]);
	    file_output=fopen_control(filename, "w");
	    corr_res=corr_history[count_history].autocorrelation_calc(1, file_output);
	}
	else
	{
	    corr_res=corr_history[count_history].autocorrelation_calc(0);
	}	
	
	if(!corr_res)
	{
	    result=false;
	    if(cur_time%delta_output==0)
	    {
		fclose(file_output);
	    }
	    break;    
	}
	else
	{
	    D->corr_lengths[count_history]=corr_history[count_history].corr_length;
	}
	if(cur_time%delta_output==0)
	{
	    fclose(file_output);
	}
    }
    if(result)
    {
	double c_r=D->largest_corr_length_calc();
	if (c_r<1.0)
	    c_r=1.0;
	double n_c=c_r*(double)N_CONF_IN_BIN;
	num_jack_samples=(int) ( (double)n_conf/n_c ); 
	if (num_jack_samples<N_BINS_MINIMUM)
	{
	    num_jack_samples=N_BINS_MINIMUM;
	}
	else if(num_jack_samples>N_BINS_MAXIMUM)
	{
	    num_jack_samples=N_BINS_MAXIMUM;
	}
    }
    else
    {    
	if (n_conf<2*N_BINS_MINIMUM)
	{
	    num_jack_samples=n_conf/2;
	}
	else
	{
	    num_jack_samples=N_BINS_MINIMUM;
	}
    }


 //output of full average and error 
{ 
  int i,j;
  double *avg=(double *)calloc(Nt_2,sizeof(double));
  double *err=(double *)calloc(Nt_2,sizeof(double));
  for(i=0;i<Nt_2;i++) {
    avg[i]=0.0;
    err[i]=0.0;
  }
  
   for(i=1;i<=Nt_2;i++) {
    for(j=0;j<=n_conf;j++) {
      avg[i-1]+=raw_data[(i-1)+(j-1)*Nt_2];
    }
    avg[i-1]=avg[i-1]/((double)n_conf);
  }
  
  for(i=1;i<=Nt_2;i++) {
    for(j=0;j<=n_conf;j++) {
      err[i-1]+=(raw_data[(i-1)+(j-1)*Nt_2]-avg[i-1])*(raw_data[(i-1)+(j-1)*Nt_2]-avg[i-1]);
    }
    err[i-1]=sqrt(err[i-1]/((double)n_conf*(n_conf-1.0)));
  }
  FILE* file_out=fopen_control("correlator_control_full_autocorr.txt","w");
  fprintf(file_out,"#time\t correlator_Re\t  error_Re (taking into account autocorrelation)\n");

  for(i=0;i<=Nt_2-1;i++) {
    fprintf(file_out,"%d\t%.15le\t%.15le\n", i+1, avg[i], err[i]*sqrt(1.0+2.0*D->corr_lengths[i]));fflush(file_out);
  }
  fclose(file_out);
  free(avg);
  free(err);
}

    
    delete[]corr_history;
    return result;
}

