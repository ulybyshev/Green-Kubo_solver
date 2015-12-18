//stable 17
#include <math.h>
class complex
{
public:
	double re;
	double im;

	complex (double par1=0.0, double par2=0.0)
	{
		re=par1;
		im=par2;
	}
	complex operator = (complex A)
	{
		 re=A.re;
		 im=A.im;
		return complex(re,im);
	}

	complex operator = (double A)
	{
		 re=A;
		 im=0.0;
		return complex(re,im);
	}

	bool operator ==(complex A)
	{
		if ( re==A.re &&  im==A.im)
			return true;
		else
			return false;
	}

	complex operator-()
	{
		return complex(-re,-im);
	}

	complex erm()
	{
		return complex (re,-im);
	}

	complex operator +(complex A)
	{
		return complex( re+A.re, im+A.im);
	}

	complex operator -(complex A)
		{
		return complex( re-A.re, im-A.im);
	}

	complex operator *(complex A)
	{
		return complex( re*A.re- im*A.im, re*A.im+ im*A.re);
	}

	complex operator /(complex A)
	{
		return complex(( re*A.re+ im*A.im)/(A.re*A.re+A.im*A.im),(- re*A.im+ im*A.re)/(A.re*A.re+A.im*A.im));
	}


	friend complex operator +(double par,complex A);
	friend complex operator -(double par,complex A);
	friend complex operator *(double par,complex A);
	friend complex operator /(double par,complex A);
	friend complex operator +(complex A,double par);
	friend complex operator -(complex A,double par);
	friend complex operator *(complex A,double par);
	friend complex operator /(complex A,double par);

};

struct complex_vector
{
    int N;
    double* re;
    double* im;
    
};


complex operator +(double par,complex A)
{
	return complex(par+A.re,A.im);
}
complex operator -(double par,complex A)
{
	return complex(par-A.re,A.im);
}
complex operator *(double par,complex A)
{
	return complex(par*A.re,par*A.im);
}
complex operator /(double par,complex A)
{
	return complex((par*A.re)/(A.re*A.re+A.im*A.im),(-par*A.im)/(A.re*A.re+A.im*A.im));
}

complex operator +(complex A,double par)
{
	return complex(par+A.re,A.im);
}
complex operator -(complex A,double par)
{
	return complex(A.re-par,A.im);
}
complex operator *(complex A,double par)
{
	return complex(par*A.re,par*A.im);
}
complex operator /(complex A,double par)
{
	return complex(A.re/par,A.im/par);
}
complex exp(complex x)
{
	return exp(x.re)*(cos(x.im)+sin(x.im)*complex(0.0,1.0));
}
double fabs(complex A)
{
	return sqrt(A.re*A.re+A.im*A.im);
}
double real(complex A)
{
  return A.re;
}

double phase (complex A)
{
    double result=0.0;
    if(A.re>0)
    {
    	result=atan(A.im/A.re);
    }
    else
    {
    if(A.im>0)
    {
	result=4.0*atan(1.0)+atan(A.im/A.re);
    }
    else
    {
	result=-4.0*atan(1.0)+atan(A.im/A.re);
    }}
    return result;
}


