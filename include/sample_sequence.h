#include <math.h>
#include <stdio.h>

//class for computing errors, averages, autocorrelation, etc. fir sample sequences
template <class element>
class sample_sequence
{
public:
    element* data;
    element average;
    element average2;//<a^2>
    element error;
    element corr_length;
    int N_elements;
    sample_sequence(int N_elements_in=1);
    ~sample_sequence();
    void format(int N_elements_in=1);
    void error_average_calc();
    element binning_error(int N_binning);
    bool autocorrelation_calc(int flag_output=0, FILE* file_output=NULL);
};



//class for computing errors, averages, autocorrelation, etc. fir sample sequences
template <class element>
sample_sequence<element>::sample_sequence(int N_elements_in)
    {
	int i=0;
	N_elements=N_elements_in;
	data=new element[N_elements_in];
	for(i=0;i<N_elements;i++)
	    data[i]=0.0;
    }
template <class element>
sample_sequence<element>:: ~sample_sequence()
    {
	delete[] data;
    }
    
template <class element>
void sample_sequence<element>::format(int N_elements_in)
    {
    	int i;
    	delete[] data;
    	N_elements=N_elements_in;
	data=new element[N_elements_in];
	for(i=0;i<N_elements;i++)
	    data[i]=0.0;
    }
    
template <class element>
void sample_sequence<element>::error_average_calc()
    {
	int i;
	average=0.0;
	for(i=0;i<N_elements;i++)
	{
	    average=average+data[i];	
	}
	
	average=average/((double)N_elements);

	average2=0.0;
	for(i=0;i<N_elements;i++)
	{
	    average2=average2+fabs(data[i])*fabs(data[i]);	
	}
	
	average2=average2/((double)N_elements);

	if(N_elements!=1)
	{
	    error=0.0;
	    for(i=0;i<N_elements;i++)
	    {
		error=error+(data[i]-average)*(data[i]-average);
	    }
	    error=sqrt(fabs(error)/((double)(N_elements*(N_elements-1))));
	}
	else
	    error=0.0;
    }
template <class element>
element sample_sequence<element>::binning_error(int N_binning)
    {
	int N_new_data;
	int i,k;
	element result;
	sample_sequence<element> new_data;
	N_new_data=this->N_elements/N_binning;

	new_data.format(N_new_data);
	
	for(i=0;i<N_new_data; i++)
	{
	    new_data.data[i]=0.0;
	    for(k=0;k<N_binning;k++)
	    {
		new_data.data[i]=new_data.data[i]+this->data[i*N_binning+k];
	    }
	    new_data.data[i]=new_data.data[i]/(double)N_binning;
	}
	new_data.error_average_calc();
	result=new_data.error;
	return result;
    }

template <class element>
bool sample_sequence<element>::autocorrelation_calc(int flag_output, FILE* file_output)
{
    if(N_elements<30)
	return false;
	
    this->error_average_calc();
    
    int N_points=N_elements/10-1;
    
    element* x;
    element* y;
    
    x=new element[N_points];
    y=new element[N_points];
    

    fprintf(file_output, "#binsize\tbin_error\n"); fflush(file_output);
    
    int count_points;
    for(count_points=2;count_points<=N_points; count_points++)
    {
	y[count_points-2]=binning_error(count_points+2)/error;
	x[count_points-2]=1.0/(element)(count_points);
    }
    
    element x2=0.0; 
    element x_av=0.0; 
    element y2=0.0; 
    element y_av=0.0;
    element xy=0.0; 
    
    for(count_points=2;count_points<=N_points; count_points++)
    {
	x2=x2+x[count_points-2]*x[count_points-2];
	y2=y2+y[count_points-2]*y[count_points-2];
	
	xy=xy+x[count_points-2]*y[count_points-2];
	
	x_av=x_av+x[count_points-2];
	y_av=y_av+y[count_points-2];
    }
    
    x2=x2/(double)(N_points-1);
    y2=y2/(double)(N_points-1);
    xy=xy/(double)(N_points-1);
    x_av=x_av/(double)(N_points-1);
    y_av=y_av/(double)(N_points-1);

    corr_length=(x2*y_av-x_av*xy)/(x2-x_av*x_av);

    if(flag_output)
    {
	for(count_points=2;count_points<=N_points; count_points++)
	{
	    fprintf(file_output, "%d\t%.15le\n", count_points, y[count_points-2]); fflush(file_output);
	} 
    }

    delete[] x;
    delete[] y;
    return true;
}

