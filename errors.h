class observ_data//class for history of observable values
{
public:
    complex* data;
    complex average;
    complex average2;//<a^2>
    complex error;
    int N_conf;
    observ_data(int N_conf_in=1)
    {
	int i=0;
	N_conf=N_conf_in;
	data=new complex[N_conf_in];
	for(i=0;i<N_conf;i++)
	    data[i]=complex(0.0,0.0);
    }
    ~observ_data()
    {
	delete[] data;
    }
    void format(int N_conf_in=1)
    {
    	int i;
    	delete[] data;
    	N_conf=N_conf_in;
	data=new complex[N_conf_in];
	for(i=0;i<N_conf;i++)
	    data[i]=complex(0.0,0.0);

    }
    void error_average_calc()
    {
	int i;
	average=complex(0.0,0.0);
	for(i=0;i<N_conf;i++)
	{
	    average=average+data[i];	
	}
	
	average=average/((double)N_conf);

	average2=complex(0.0,0.0);
	for(i=0;i<N_conf;i++)
	{
	    average2=average2+complex(data[i].re*data[i].re,data[i].im*data[i].im);	
	}
	
	average2=average2/((double)N_conf);

if(N_conf!=1)
{
	error=complex(0.0,0.0);
	for(i=0;i<N_conf;i++)
	{
	    error.re=error.re+(data[i].re-average.re)*(data[i].re-average.re);
	    error.im=error.im+(data[i].im-average.im)*(data[i].im-average.im);
	}
	
	error=complex(sqrt(error.re/((double)(N_conf*(N_conf-1)))),sqrt(error.im/((double)(N_conf*(N_conf-1)))) );
}
else
	error=complex(0.0,0.0 );


    }

    complex binning_error(int N_binning)
    {
	int N_new_data;
	int i,k;
	complex result;
	observ_data new_data;
	N_new_data=this->N_conf/N_binning;

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
};




//we have 2d massive G[Npoints][Nvalues], this function makes covariance matrix calculation
void cov_matrix_calculation(int Npoints, int Nvalues,double** cov_matrix, double** G)
{
    int count_point;
    int count_value;
    int count_m1, count_m2;
    double res;
    double* av_values;
    
    av_values=(double*) calloc(Npoints, sizeof(double));
        
    for(count_m1=0;count_m1<Npoints;count_m1++)
    for(count_m2=0;count_m2<Npoints;count_m2++)
    {
	cov_matrix[count_m1][count_m2]=0.0;
    }
    
    for(count_point=0;count_point<Npoints;count_point++)
    {
	for(count_value=0;count_value<Nvalues;count_value++)
	{
	av_values[count_point]+=G[count_point][count_value];
	}
	av_values[count_point]=av_values[count_point]/(double)Nvalues;
    }

    for(count_m1=0;count_m1<Npoints;count_m1++)
    for(count_m2=0;count_m2<Npoints;count_m2++)
    {
	for(count_value=0;count_value<Nvalues;count_value++)
	{
	cov_matrix[count_m1][count_m2]+=G[count_m1][count_value] * G[count_m2][count_value];
	}
	cov_matrix[count_m1][count_m2]=cov_matrix[count_m1][count_m2]/(double)Nvalues - av_values[count_m1]* av_values[count_m2];
    }

    free(av_values);
}