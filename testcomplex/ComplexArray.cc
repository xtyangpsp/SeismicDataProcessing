#include "ComplexArray.h"

ComplexArray::ComplexArray()
{
	nsamp=0;
	data=NULL;
}
ComplexArray::ComplexArray(vector<Complex64> &d)
{
	nsamp=d.size();
	data=new FortranComplex64[nsamp];
	for(int i=0;i<nsamp;i++)
	{
		data[i].real=d[i].real();
		data[i].imag=d[i].imag();
	}
}
ComplexArray::ComplexArray(vector<Complex32> &d)
{
	nsamp=d.size();
	data=new FortranComplex64[nsamp];
	for(int i=0;i<nsamp;i++)
	{
		data[i].real=d[i].real();
		data[i].imag=d[i].imag();
	}
}
ComplexArray::ComplexArray(int n, FortranComplex32 *d)
{
	nsamp=n;
	data=new FortranComplex64[nsamp];
	for(int i=0;i<nsamp;i++)
	{
		data[i].real=d[i].real;
		data[i].imag=d[i].imag;
	}
}
ComplexArray::ComplexArray(int n, FortranComplex64 *d)
{
	nsamp=n;
	data=new FortranComplex64[nsamp];
	for(int i=0;i<nsamp;i++)
	{
		data[i].real=d[i].real;
		data[i].imag=d[i].imag;
	}
}
ComplexArray::ComplexArray(int n, float *d)
{
	nsamp=n;
	data=new FortranComplex64[nsamp];
	for(int i=0;i<nsamp;i++)
	{
		data[i].real=d[i];
		data[i].imag=0.0;
	}
}
ComplexArray::ComplexArray(int n, double *d)
{
	nsamp=n;
	data=new FortranComplex64[nsamp];
	for(int i=0;i<nsamp;i++)
	{
		data[i].real=d[i];
		data[i].imag=0.0;
	}
}
ComplexArray::ComplexArray(int n)
{
	nsamp=n;
	data=new FortranComplex64[nsamp];
	for(int i=0;i<nsamp;i++)
	{
		data[i].real=0.0;
		data[i].imag=0.0;
	}
}

ComplexArray::ComplexArray(vector<double> mag,vector<double> phase)
{
	nsamp=mag.size();
	if(nsamp==phase.size())
	{
		data=new FortranComplex64[nsamp];
		for(int i=0;i<nsamp;i++)
		{
			Complex64 temp=polar(mag[i],phase[i]);
			data[i].real=temp.real();
			data[i].imag=temp.imag();
		}
	}
	else
	{
		cout<<"Length of magnitude vector and phase vector doesn't match"<<endl;
		throw SeisppError("ComplexArray::ComplexArray(vector<double> mag,vector<double> phase): Length of magnitude vector and phase vector doesn't match");
	}
}
ComplexArray::ComplexArray(const ComplexArray &parent)
{
	nsamp=parent.nsamp;
	data=new FortranComplex64[nsamp];
	for(int i=0;i<nsamp;i++)
	{
		data[i].real=parent.data[i].real;
		data[i].imag=parent.data[i].imag;
	}
}
ComplexArray& ComplexArray::operator=(const ComplexArray &parent)
{
	nsamp=parent.nsamp;
	delete [] data;
	data=new FortranComplex64[nsamp];
	for(int i=0;i<nsamp;i++)
	{
		data[i].real=parent.data[i].real;
		data[i].imag=parent.data[i].imag;
	}
	return *this;
}
ComplexArray::~ComplexArray()
{
	delete[] data;
}
double *ComplexArray::ptr()
{
	return reinterpret_cast<double*>(&data[0].real);
}
double *ComplexArray::ptr(int sample)
{
	return reinterpret_cast<double*>(&data[sample].real);
}
Complex64 ComplexArray::operator[](int sample)
{
	return *reinterpret_cast<Complex64*>(&data[sample].real);
}
ComplexArray& ComplexArray::operator +=(const ComplexArray& other)
{
	int n;
	if(nsamp>other.nsamp)
		n=other.nsamp;
	else
		n=nsamp;
	for(int i=0;i<n;i++)
	{
		data[i].real+=other.data[i].real;
		data[i].imag+=other.data[i].imag;
	}
	return *this;
}
ComplexArray& ComplexArray::operator -=(const ComplexArray& other)
{
	int n;
	if(nsamp>other.nsamp)
		n=other.nsamp;
	else
		n=nsamp;
	for(int i=0;i<n;i++)
	{
		data[i].real-=other.data[i].real;
		data[i].imag-=other.data[i].imag;
	}
	return *this;
}
ComplexArray ComplexArray::operator +(const ComplexArray& other)
{
	ComplexArray result(*this);
	int n;
	if(nsamp>other.nsamp)
		n=other.nsamp;
	else
		n=nsamp;
	for(int i=0;i<n;i++)
	{
		result.data[i].real=data[i].real+other.data[i].real;
		result.data[i].imag=data[i].imag+other.data[i].imag;
	}
	return result;
}
ComplexArray ComplexArray::operator -(const ComplexArray& other)
{
	ComplexArray result(*this);
	int n;
	if(nsamp>other.nsamp)
		n=other.nsamp;
	else
		n=nsamp;
	for(int i=0;i<n;i++)
	{
		result.data[i].real=data[i].real-other.data[i].real;
		result.data[i].imag=data[i].imag-other.data[i].imag;
	}
	return result;
}
ComplexArray ComplexArray::operator *(const ComplexArray& other)
{
	ComplexArray result(*this);
	int n;
	if(nsamp>other.nsamp)
		n=other.nsamp;
	else
		n=nsamp;
	for(int i=0;i<n;i++)
	{
		Complex64 temp1(data[i].real,data[i].imag);
		Complex64 temp2(other.data[i].real,other.data[i].imag);
		temp1*=temp2;
		result.data[i].real=temp1.real();
		result.data[i].imag=temp1.imag();
	}
	return result;
}
ComplexArray ComplexArray::operator /(const ComplexArray& other)
{
	ComplexArray result;
	int n;
	if(nsamp>other.nsamp)
	{
		n=other.nsamp;
		result=other;
	}
	else
	{
		n=nsamp;
		result=*this;
	}
	for(int i=0;i<n;i++)
	{
		Complex64 temp1(data[i].real,data[i].imag);
		Complex64 temp2(other.data[i].real,other.data[i].imag);
		temp1/=temp2;
		result.data[i].real=temp1.real();
		result.data[i].imag=temp1.imag();
		if(gsl_isnan(temp1.real()))
			result.data[i].real=0;
		else
			result.data[i].real=temp1.real();
		if(gsl_isnan(temp1.imag()))
			result.data[i].imag=0;
		else
			result.data[i].imag=temp1.imag();
	}
	return result;
}
void ComplexArray::conj()
{
	for(int i=0;i<nsamp;i++)
		data[i].imag=-data[i].imag;
}
vector<double> ComplexArray::abs()
{
	vector<double> result;
    result.reserve(nsamp);
	for(int i=0;i<nsamp;i++)
		result.push_back(sqrt((double)data[i].real*data[i].real+data[i].imag*data[i].imag));
	return result;
}
double ComplexArray::rms()
{
	double result=0;
	for(int i=0;i<nsamp;i++)
		result+=((double)data[i].real*data[i].real+data[i].imag*data[i].imag)/nsamp/nsamp;
	return sqrt(result);
}
double ComplexArray::norm2()
{
	double result=0;
	for(int i=0;i<nsamp;i++)
		result+=((double)data[i].real*data[i].real+data[i].imag*data[i].imag);
	return sqrt(result);
}
vector<double> ComplexArray::phase()
{
	vector<double> result;
    result.reserve(nsamp);
	for(int i=0;i<nsamp;i++)
		result.push_back(atan2((double)data[i].imag,(double)data[i].real));
	return result;
}
int ComplexArray::size()
{
	return nsamp;
}
