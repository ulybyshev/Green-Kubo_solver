#include "basic_structures.h"
#include "constants.h"

interval::interval(int size_interval)
{
    size=size_interval;
    times=(int*) calloc(size_interval, sizeof(int));    
}
void interval::format(int size_interval)
{
    size=size_interval;
    free(times);
    times=(int*) calloc(size_interval, sizeof(int));    
}
interval::~interval()
{
    free(times);
}



correlator::correlator(int N_full, int N_valid)
{
    N_full_points=N_full;
    N_valid_points=N_valid;
    length=(double) N_full;
    
    corr_full=(double*) calloc(N_full_points, sizeof(double));
    error_full=(double*) calloc(N_full_points, sizeof(double));
    S_full=gsl_matrix_alloc(N_full_points, N_full_points);

    corr=(double*) calloc(N_valid_points, sizeof(double));
    error=(double*) calloc(N_valid_points, sizeof(double));
    S=gsl_matrix_alloc(N_valid_points, N_valid_points);
    
    points_numbers=(int*)calloc(N_valid_points, sizeof(int));

    interval_numbers=new interval*[N_valid_points];
    int i;
    for(i=0;i<N_valid_points;i++)
	interval_numbers[i]=new interval;
    
}
void correlator::format(int N_full,int N_valid)
{
    free(corr_full);
    free(corr);
    
    free(error);
    free(error_full);
    
    gsl_matrix_free(S_full);
    gsl_matrix_free(S);
    
    free(points_numbers);

    int i;
    for(i=0;i<N_valid_points;i++)
	delete interval_numbers[i];
    delete[] interval_numbers;

    
    N_full_points=N_full;
    N_valid_points=N_valid;
    length=(double) N_full;
    
    corr_full=(double*) calloc(N_full_points, sizeof(double));
    error_full=(double*) calloc(N_full_points, sizeof(double));
    S_full=gsl_matrix_alloc(N_full_points, N_full_points);

    corr=(double*) calloc(N_valid_points, sizeof(double));
    error=(double*) calloc(N_valid_points, sizeof(double));
    S=gsl_matrix_alloc(N_valid_points, N_valid_points);
    
    points_numbers=(int*)calloc(N_valid_points, sizeof(int));

    interval_numbers=new interval*[N_valid_points];
    for(i=0;i<N_valid_points;i++)
	interval_numbers[i]=new interval;

}


correlator:: ~correlator()
{
    free(corr_full);
    free(corr);
    
    free(error);
    free(error_full);
    
    gsl_matrix_free(S_full);
    gsl_matrix_free(S);
    
    free(points_numbers);
    int i;
    for(i=0;i<N_valid_points;i++)
	delete interval_numbers[i];

    delete[] interval_numbers;

}
    
    

void correlator::construct_intervals()//int *points, int num_points, interval *intervals) 
{

  int i,j,k,count_interval,interval_length;

  count_interval=N_valid_points-1;
  for(i=(N_valid_points-1);i>=1;i--) 
  {
	if((points_numbers[i]-1)==points_numbers[i-1]) 
	{ //single "red" point
    	    interval_numbers[count_interval]->format(1);
    	    interval_numbers[count_interval]->times[0]=points_numbers[i];
	}
	else
	{ //contains "black points"
    	    interval_length=(points_numbers[i]-points_numbers[i-1]);
    	    interval_numbers[count_interval]->format(interval_length);
    	    for(j=points_numbers[i],k=(interval_length-1);j>points_numbers[i-1];j--,k--)
    	    {
		interval_numbers[count_interval]->times[k]=j;
    	    }
	}
	count_interval--;
  }
  interval_numbers[0]->format(1);
  interval_numbers[0]->times[0]=points_numbers[0];
}


calc_structures::calc_structures(int N_t)
{
    N_t_points=N_t;
    length=(double)N_t;
    R=gsl_vector_alloc(N_t);
    omega_R=gsl_vector_alloc(N_t);
    
    int i;
    N_center=(int) ((center_stop-center_start)/center_delta);
    center=(double*) calloc(N_center, sizeof(double));

    center[0]=center_start;
    for(i=1; i<N_center; i++)
	center[i]=center[i-1]+center_delta;

    W=(gsl_matrix**) calloc(N_center, sizeof(gsl_matrix*));
    for(i=0; i<N_center; i++)
	W[i]=gsl_matrix_alloc(N_t_points, N_t_points);    
}    
calc_structures::~calc_structures()
{
    gsl_vector_free(R);
    gsl_vector_free(omega_R);
    free(center);
    int i;
    for(i=0; i<N_center; i++)
	gsl_matrix_free(W[i]);
    free(W);
    
}


