#ifndef __DECON_OPERATOR_H__
#define __DECON_OPERATOR_H__
#include <vector>
#include <cmath>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_cblas.h>
#include "dpss.h"
#include "seispp.h"
#include "Metadata.h"
#include "ComplexArray.h"

using namespace std;
using namespace SEISPP;

enum IterDeconType {XCORR,WATER_LEVEL,LEAST_SQ,MULTI_TAPER};

double *rickerwavelet(float fpeak, float dt, int n);
double *gaussian(float sigma, float dt, int n);
double power(double* data,int len);
ComplexArray gaussianFilter(ComplexArray x,float gwidthFactor,float dt);

class CoreDeconOperator
{
public:
	virtual void changeparameter(Metadata &md)=0;
	vector<double> getresult(){return result;}
	virtual ~CoreDeconOperator(){}
protected:
	/* data is set as a pointer to improve efficiency */
	vector<double> *data;
	vector<double> *wavelet;
	vector<double> result;
};

class SimpleDecon: public CoreDeconOperator	
{
public:
	int load(vector<double> &wavelet,vector<double> &data);
	~SimpleDecon(){}
protected:
	//should be defined in the main function by, metadata.put
	unsigned int nfft;
	int sample_shift;
private:
	/* Apply the deconvolution. This will be called in load method. */
	virtual int apply()=0;
};

class ArrayDecon: public CoreDeconOperator
{
public:
	int load(vector<double> &beam);
	int loaddata(vector<double> &data);
	~ArrayDecon(){}
protected:
	ComplexArray invbeam;
	//should be defined in the main function by, metadata.put
	unsigned int nfft;
	int sample_shift;
private:
	/* Apply the deconvolution. This will be called in loaddata method. */
	virtual int apply()=0;
	/* Calculate the inverse operator and store it. This will be called in load method. */
	virtual int getinverse()=0;
};

class SimpleLeastSquareDecon: public SimpleDecon
{
public:
	SimpleLeastSquareDecon(const SimpleLeastSquareDecon &parent);
	SimpleLeastSquareDecon(Metadata &md);
	SimpleLeastSquareDecon(Metadata &md,vector<double> &wavelet,vector<double> &data);
	void changeparameter(Metadata &md);
	~SimpleLeastSquareDecon(){delete shapingwavelet;}
private:
	int read_metadata(Metadata &md);
	int apply();
	float damp;
	ComplexArray *shapingwavelet;
};

class SimpleWaterLevelDecon: public SimpleDecon
{
public:
	SimpleWaterLevelDecon(const SimpleWaterLevelDecon &parent);
	SimpleWaterLevelDecon(Metadata &md);
	SimpleWaterLevelDecon(Metadata &md,vector<double> &wavelet,vector<double> &data);
	void changeparameter(Metadata &md);
	~SimpleWaterLevelDecon(){delete shapingwavelet;}
private:
	int read_metadata(Metadata &md);
	int apply();
	float wlv;
	ComplexArray *shapingwavelet;
};

class SimpleMultiTaperDecon: public SimpleDecon
{
public:
	SimpleMultiTaperDecon(const SimpleMultiTaperDecon &parent);
	SimpleMultiTaperDecon(Metadata &md);
	SimpleMultiTaperDecon(Metadata &md,vector<double> &noise,vector<double> &wavelet,vector<double> &data);
	void changeparameter(Metadata &md);
	int loadnoise(vector<double> &noise);	
	~SimpleMultiTaperDecon(){delete shapingwavelet;delete [] tapers;}
private:
	int read_metadata(Metadata &md);
	int apply();
	ComplexArray *shapingwavelet;
	vector<double> *noise;
	float nw,damp;
	int seql;
	unsigned int taperlen;
	double *tapers;
};

class SimpleGeneralIterDecon: public SimpleDecon
{
public:
	int numberiter(){return(icount);}
	double epsilon(){return(eps);}
	SimpleGeneralIterDecon(const SimpleGeneralIterDecon &parent);
	SimpleGeneralIterDecon(Metadata &md);
	SimpleGeneralIterDecon(Metadata &md,vector<double> &wavelet,vector<double> &data);
	SimpleGeneralIterDecon(Metadata &md,vector<double> &noise,vector<double> &wavelet,vector<double> &data);
	void changeparameter(Metadata &md);
	int loadnoise(vector<double> &noise);
	~SimpleGeneralIterDecon(){delete shapingwavelet;delete invshapingwavelet;delete [] tapers;}
private:
	int read_metadata(Metadata &md);
	int apply();
	ComplexArray *shapingwavelet;
	ComplexArray *invshapingwavelet;
	int itermax;
	float tol;
	float tolinstep;
	
	int icount;
	double eps;
	
	IterDeconType type;
	vector<double> *noise;
	float nw,damp;
	int seql;
	unsigned int taperlen;
	double *tapers;
};

//here begins different ArrayDecon

class ArrayLeastSquareDecon: public ArrayDecon
{
public:
	ArrayLeastSquareDecon(const ArrayLeastSquareDecon &parent);
	ArrayLeastSquareDecon(Metadata &md);
	ArrayLeastSquareDecon(Metadata &md,vector<double> &beam,vector<double> &data);
	void changeparameter(Metadata &md);
	~ArrayLeastSquareDecon(){delete shapingwavelet;}
private:
	int read_metadata(Metadata &md);
	int apply();
	int getinverse();
	float damp;
	ComplexArray *shapingwavelet;
};

class ArrayWaterLevelDecon: public ArrayDecon
{
public:
	ArrayWaterLevelDecon(const ArrayWaterLevelDecon &parent);
	ArrayWaterLevelDecon(Metadata &md);
	ArrayWaterLevelDecon(Metadata &md,vector<double> &beam,vector<double> &data);
	void changeparameter(Metadata &md);
	~ArrayWaterLevelDecon(){delete shapingwavelet;}
private:
	int read_metadata(Metadata &md);
	int apply();
	int getinverse();
	float wlv;
	ComplexArray *shapingwavelet;
};

class ArrayMultiTaperDecon: public ArrayDecon
{
public:
	ArrayMultiTaperDecon(const ArrayMultiTaperDecon &parent);
	ArrayMultiTaperDecon(Metadata &md);
	ArrayMultiTaperDecon(Metadata &md,vector<double> &noise,vector<double> &beam,vector<double> &data);
	void changeparameter(Metadata &md);
	int loadnoise(vector<double> &noise);	
	~ArrayMultiTaperDecon(){delete shapingwavelet;delete [] tapers;}
private:
	int read_metadata(Metadata &md);
	int apply();
	int getinverse();
	ComplexArray *shapingwavelet;
	vector<double> *noise;
	float nw,damp;
	int seql;
	//should be defined in the main function best when taperlen=data->size()
	unsigned int taperlen;
	double *tapers;
};

class ArrayGeneralIterDecon: public ArrayDecon
{
public:
	int numberiter(){return(icount);}
	double epsilon(){return(eps);}
	ArrayGeneralIterDecon(const ArrayGeneralIterDecon &parent);
	ArrayGeneralIterDecon(Metadata &md);
	ArrayGeneralIterDecon(Metadata &md,vector<double> &beam,vector<double> &data);
	ArrayGeneralIterDecon(Metadata &md,vector<double> &noise,vector<double> &beam,vector<double> &data);
	void changeparameter(Metadata &md);
	int loadnoise(vector<double> &noise);
	~ArrayGeneralIterDecon(){delete shapingwavelet;delete invshapingwavelet;delete [] tapers;}
private:
	int read_metadata(Metadata &md);
	int apply();
	int getinverse();
	ComplexArray *shapingwavelet;
	ComplexArray *invshapingwavelet;
	int itermax;
	float tol;
	float tolinstep;
	
	int icount;
	double eps;
	
	IterDeconType type;
	vector<double> *noise;
	float nw,damp;
	int seql;
	//should be defined in the main function best when taperlen=data->size()
	unsigned int taperlen;
	double *tapers;
};
#endif
