#include "DeconOperator.h"

using namespace std;
using namespace SEISPP;
//SimpleDecon
int SimpleDecon::load(vector<double> &w,vector<double> &d)
{
	wavelet=&w;
	data=&d;
	result.clear();
	if(this->apply())
		return 1;
	else
		return 0;
}
//ArrayDecon
int ArrayDecon::load(vector<double> &beam)
{
	wavelet=&beam;
	if(this->getinverse())
		return 1;
	else
		return 0;
}
int ArrayDecon::loaddata(vector<double> &d)
{
	data=&d;
	result.clear();
	if(this->apply())
		return 1;
	else
		return 0;
}
//SimpleLeastSquareDecon
SimpleLeastSquareDecon::SimpleLeastSquareDecon(const SimpleLeastSquareDecon &parent)
{
	data=parent.data;
	wavelet=parent.wavelet;
	if(parent.shapingwavelet!=NULL)
		shapingwavelet=new ComplexArray(*(parent.shapingwavelet));
	else
		shapingwavelet=NULL;
	damp=parent.damp;
	nfft=parent.nfft;
    sample_shift=parent.sample_shift;
	result=parent.result;
}
int SimpleLeastSquareDecon::read_metadata(Metadata &md)
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	try{
		damp=md.get_double("damping_factor");
		nfft=md.get_int("operator_nfft");
		sample_shift=md.get_int("sample_shift");
		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
	}catch(SeisppError& err)
	{
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleLeastSquareDecon::read_metadata"<<endl;
		exit(-1);
	}
	try{
		string wavelettype=md.get_string("shaping_wavelet_type");
		if(wavelettype=="gaussian")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=gaussian(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else if(wavelettype=="ricker")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=rickerwavelet(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else
		{
			cout<<"WARNING: Using the raw output of deconvolution"<<endl;
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,1.0);
		}
	}catch(SeisppError& err)
	{
		cerr<<err.what()<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}catch(...)
	{
		cerr<<"Unknown error from SimpleLeastSquareDecon::read_metadata"<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
SimpleLeastSquareDecon::SimpleLeastSquareDecon(Metadata &md)
{
	data=NULL;
	wavelet=NULL;
	shapingwavelet=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleLeastSquareDecon constructor"<<endl;
		exit(-1);
	}
}
void SimpleLeastSquareDecon::changeparameter(Metadata &md)
{
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot change the decon parameters ";
		throw err;
	}catch(...)
	{
		throw SeisppError("Error at SimpleLeastSquareDecon::changeparameter");
	}
}
SimpleLeastSquareDecon::SimpleLeastSquareDecon(Metadata &md,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleLeastSquareDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
int SimpleLeastSquareDecon::apply()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	//apply fft to the input trace data
	ComplexArray d_fft(nfft,*data);
	gsl_fft_complex_forward(d_fft.ptr(), 1, nfft, wavetable, workspace);
	
	//apply fft to wavelet
	ComplexArray b_fft(nfft,*wavelet);
	gsl_fft_complex_forward(b_fft.ptr(), 1, nfft, wavetable, workspace);
	
	//deconvolution: RF=conj(B).*D./(conj(B).*B+damp) 
	b_fft.conj();
	ComplexArray rf_fft,conj_b_fft(b_fft);
	b_fft.conj();
	
    double b_rms=b_fft.rms();
	rf_fft=(conj_b_fft*d_fft)/(conj_b_fft*b_fft+b_rms*damp);
	
	//apply shaping wavelet
	rf_fft=*shapingwavelet*rf_fft;
	
	//ifft gets result
	gsl_fft_complex_inverse(rf_fft.ptr(), 1, nfft, wavetable, workspace);
	if(sample_shift>0)
    {
        for(int k=sample_shift;k>0;k--)
            result.push_back(rf_fft[nfft-k].real());
        for(int k=0;k<data->size()-sample_shift;k++)
            result.push_back(rf_fft[k].real());
    }
    else
	{
        for(int k=0;k<data->size();k++)
            result.push_back(rf_fft[k].real());
    }
	
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
//SimpleWaterLevelDecon
SimpleWaterLevelDecon::SimpleWaterLevelDecon(const SimpleWaterLevelDecon &parent)
{
	data=parent.data;
	wavelet=parent.wavelet;
	if(parent.shapingwavelet!=NULL)
		shapingwavelet=new ComplexArray(*(parent.shapingwavelet));
	else
		shapingwavelet=NULL;
	wlv=parent.wlv;
	nfft=parent.nfft;
    sample_shift=parent.sample_shift;
	result=parent.result;
}
int SimpleWaterLevelDecon::read_metadata(Metadata &md)
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	try{
		wlv=md.get_double("damping_factor");
		nfft=md.get_int("operator_nfft");
		sample_shift=md.get_int("sample_shift");
		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
	}catch(SeisppError& err)
	{
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleWaterLevelDecon::read_metadata"<<endl;
		exit(-1);
	}
	try{
		string wavelettype=md.get_string("shaping_wavelet_type");
		if(wavelettype=="gaussian")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=gaussian(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else if(wavelettype=="ricker")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=rickerwavelet(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else
		{
			cout<<"WARNING: Using the raw output of deconvolution"<<endl;
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,1.0);
		}
	}catch(SeisppError& err)
	{
		cerr<<err.what()<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}catch(...)
	{
		cerr<<"Unknown error from SimpleWaterLevelDecon::read_metadata"<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
SimpleWaterLevelDecon::SimpleWaterLevelDecon(Metadata &md)
{
	data=NULL;
	wavelet=NULL;
	shapingwavelet=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleWaterLevelDecon constructor"<<endl;
		exit(-1);
	}
}
void SimpleWaterLevelDecon::changeparameter(Metadata &md)
{
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot change the decon parameters ";
		throw err;
	}catch(...)
	{
		throw SeisppError("Error at SimpleWaterLevelDecon::changeparameter");
	}
}
SimpleWaterLevelDecon::SimpleWaterLevelDecon(Metadata &md,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleWaterLevelDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
int SimpleWaterLevelDecon::apply()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	//apply fft to the input trace data
	ComplexArray d_fft(nfft,*data);
	gsl_fft_complex_forward(d_fft.ptr(), 1, nfft, wavetable, workspace);
	
	//apply fft to wavelet
	ComplexArray b_fft(nfft,*wavelet);
	gsl_fft_complex_forward(b_fft.ptr(), 1, nfft, wavetable, workspace);
	
    double b_rms=b_fft.rms();
    
    //water level
    for(int i=0;i<nfft;i++)
    {
        if(abs(b_fft[i])<b_rms*wlv)
        {
            *b_fft.ptr(i)=*b_fft.ptr(i)/abs(b_fft[i])*b_rms*wlv;
			*(b_fft.ptr(i)+1)=*(b_fft.ptr(i)+1)/abs(b_fft[i])*b_rms*wlv;
        }
    }
    
	//deconvolution: RF=D./(B+wlv)
	ComplexArray rf_fft;
	rf_fft=d_fft/b_fft;
	
	//apply shaping wavelet
	rf_fft=*shapingwavelet*rf_fft;
	
	//ifft gets result
	gsl_fft_complex_inverse(rf_fft.ptr(), 1, nfft, wavetable, workspace);
	if(sample_shift>0)
    {
        for(int k=sample_shift;k>0;k--)
            result.push_back(rf_fft[nfft-k].real());
        for(int k=0;k<data->size()-sample_shift;k++)
            result.push_back(rf_fft[k].real());
    }
    else
	{
        for(int k=0;k<data->size();k++)
            result.push_back(rf_fft[k].real());
    }
	
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
//SimpleMultiTaperDecon
SimpleMultiTaperDecon::SimpleMultiTaperDecon(const SimpleMultiTaperDecon &parent)
{
	data=parent.data;
	wavelet=parent.wavelet;
	noise=parent.noise;
	if(parent.shapingwavelet!=NULL)
		shapingwavelet=new ComplexArray(*(parent.shapingwavelet));
	else
		shapingwavelet=NULL;
	nw=parent.nw;
	seql=parent.seql;
	taperlen=parent.taperlen;
	if(parent.tapers!=NULL)
	{
		int sequ=nw*2-2;
		int nseq=sequ-seql+1;
		tapers=new double[nseq*taperlen];
		cblas_dcopy(nseq*taperlen,parent.tapers,1,tapers,1);
	}
	else
		tapers=NULL;
	damp=parent.damp;
	nfft=parent.nfft;
    sample_shift=parent.sample_shift;
	result=parent.result;
}
int SimpleMultiTaperDecon::read_metadata(Metadata &md)
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	try{
		damp=md.get_double("damping_factor");
		nfft=md.get_int("operator_nfft");
		sample_shift=md.get_int("sample_shift");
		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
		nw=md.get_double("time_bandwidth_product");
		nw=static_cast<int>(nw*2)/2.0;
		seql=md.get_int("lower_dpss");
		if(seql>2*nw-2)
			seql=static_cast<int>(2*nw-2);
        taperlen=md.get_int("taper_length");
        int sequ=nw*2-2;
        int nseq=sequ-seql+1;
        if(tapers!=NULL)
            delete [] tapers;
        tapers=new double[nseq*taperlen];
        dpss_calc(taperlen, nw, seql, sequ, tapers);
	}catch(SeisppError& err)
	{
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleMultiTaperDecon::read_metadata"<<endl;
		exit(-1);
	}
	try{
		string wavelettype=md.get_string("shaping_wavelet_type");
		if(wavelettype=="gaussian")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=gaussian(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else if(wavelettype=="ricker")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=rickerwavelet(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else
		{
			cout<<"WARNING: Using the raw output of deconvolution"<<endl;
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,1.0);
		}
	}catch(SeisppError& err)
	{
		cerr<<err.what()<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}catch(...)
	{
		cerr<<"Unknown error from SimpleMultiTaperDecon::read_metadata"<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
SimpleMultiTaperDecon::SimpleMultiTaperDecon(Metadata &md)
{
	data=NULL;
	wavelet=NULL;
	shapingwavelet=NULL;
	noise=NULL;
	tapers=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleMultiTaperDecon constructor"<<endl;
		exit(-1);
	}
}
void SimpleMultiTaperDecon::changeparameter(Metadata &md)
{
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot change the decon parameters ";
		throw err;
	}catch(...)
	{
		throw SeisppError("Error at SimpleMultiTaperDecon::changeparameter");
	}
}
int SimpleMultiTaperDecon::loadnoise(vector<double> &n)
{
	noise=&n;
	return 0;
}
SimpleMultiTaperDecon::SimpleMultiTaperDecon(Metadata &md,vector<double> &n,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	tapers=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleMultiTaperDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	noise=&n;
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
int SimpleMultiTaperDecon::apply()
{
	if(noise==NULL)
	{
		cerr<<"SimpleMultiTaperDecon::apply: noise data is empty."<<endl;
		throw SeisppError("SimpleMultiTaperDecon::apply: noise data is empty.");
		return 1;
	}
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	int sequ=nw*2-2;
	int nseq=sequ-seql+1;

	//apply fft to the input trace data
	ComplexArray od_fft(nfft,*data);
	gsl_fft_complex_forward(od_fft.ptr(), 1, nfft, wavetable, workspace);
	//apply tapers to the data
	ComplexArray tdata(nfft*nseq);
	ComplexArray d_fft(nfft);
	for(int i=0;i<nseq;i++)
	{
		for(int k=0;k<taperlen;k++)
		{
			*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*data->at(k);
		}
		gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
	}
	for(int i=0;i<nfft;i++)
	{
		double sumr=0,sumi=0;
		for(int k=0;k<nseq;k++)
		{
			sumr+=tdata[k*nfft+i].real();
			sumi+=tdata[k*nfft+i].imag();
		}
		*d_fft.ptr(i)=sumr/nseq;
		*(d_fft.ptr(i)+1)=sumi/nseq;
	}
	double odmax=abs(od_fft[0]),dmax=abs(d_fft[0]);
	for(int i=1;i<nfft;i++)
	{
		if(abs(od_fft[i])>odmax)
			odmax=abs(od_fft[i]);
		if(abs(d_fft[i])>dmax)
			dmax=abs(d_fft[i]);
	}
	d_fft=d_fft*(odmax/dmax);
	
	//apply fft to beam
	ComplexArray ob_fft(nfft,*wavelet);
	gsl_fft_complex_forward(ob_fft.ptr(), 1, nfft, wavetable, workspace);
	//apply tapers to the beam
	ComplexArray b_fft(nfft);
	int waveletlen;
	if(wavelet->size()<taperlen)
		waveletlen=wavelet->size();
	else
		waveletlen=taperlen;
	for(int i=0;i<nseq;i++)
	{
		for(int k=0;k<nfft;k++)
		{
			if(k<waveletlen)
				*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*wavelet->at(k);
			else
				*tdata.ptr(i*nfft+k)=0.0;
			*(tdata.ptr(i*nfft+k)+1)=0.0;
		}
		gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
	}
	for(int i=0;i<nfft;i++)
	{
		double sumr=0,sumi=0;
		for(int k=0;k<nseq;k++)
		{
			sumr+=tdata[k*nfft+i].real();
			sumi+=tdata[k*nfft+i].imag();
		}
		*b_fft.ptr(i)=sumr/nseq;
		*(b_fft.ptr(i)+1)=sumi/nseq;
	}
	double obmax=abs(ob_fft[0]),bmax=abs(b_fft[0]);
	for(int i=1;i<nfft;i++)
	{
		if(abs(ob_fft[i])>obmax)
			obmax=abs(ob_fft[i]);
		if(abs(b_fft[i])>bmax)
			bmax=abs(b_fft[i]);
	}
	b_fft=b_fft*(obmax/bmax);
	
	//apply fft to noise
	ComplexArray on_fft(nfft,*noise);
	gsl_fft_complex_forward(on_fft.ptr(), 1, nfft, wavetable, workspace);
	//apply tapers to the noise
	ComplexArray n_fft(nfft);
	int noiselen;
	if(noise->size()<taperlen)
		noiselen=noise->size();
	else
		noiselen=taperlen;
	for(int i=0;i<nseq;i++)
	{
		for(int k=0;k<nfft;k++)
		{
			if(k<noiselen)
				*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*noise->at(k);
			else
				*tdata.ptr(i*nfft+k)=0.0;
			*(tdata.ptr(i*nfft+k)+1)=0.0;
		}
		gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
	}
	for(int i=0;i<nfft;i++)
	{
		double sumr=0,sumi=0;
		for(int k=0;k<nseq;k++)
		{
			sumr+=tdata[k*nfft+i].real();
			sumi+=tdata[k*nfft+i].imag();
		}
		*n_fft.ptr(i)=sumr/nseq;
		*(n_fft.ptr(i)+1)=sumi/nseq;
	}
	double onmax=abs(on_fft[0]),nmax=abs(n_fft[0]);
	for(int i=1;i<nfft;i++)
	{
		if(abs(on_fft[i])>onmax)
			onmax=abs(on_fft[i]);
		if(abs(n_fft[i])>nmax)
			nmax=abs(n_fft[i]);
	}
	n_fft=n_fft*(onmax/nmax);
	
	//deconvolution: RF=conj(B).*D./(conj(B).*B+abs(S)) 
	b_fft.conj();
	ComplexArray rf_fft,conj_b_fft(b_fft);
	b_fft.conj();
	
	n_fft=n_fft*damp;
	
	rf_fft=(conj_b_fft*d_fft)/(conj_b_fft*b_fft+n_fft.abs());
	
	//apply shaping wavelet
	rf_fft=*shapingwavelet*rf_fft;
	
	//ifft gets result
	gsl_fft_complex_inverse(rf_fft.ptr(), 1, nfft, wavetable, workspace);
    if(sample_shift>0)
    {
        for(int k=sample_shift;k>0;k--)
            result.push_back(rf_fft[nfft-k].real());
        for(int k=0;k<data->size()-sample_shift;k++)
            result.push_back(rf_fft[k].real());
    }
    else
	{
        for(int k=0;k<data->size();k++)
            result.push_back(rf_fft[k].real());
    }
	
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
//SimpleGeneralIterDecon
SimpleGeneralIterDecon::SimpleGeneralIterDecon(const SimpleGeneralIterDecon &parent)
{
	data=parent.data;
	wavelet=parent.wavelet;
	if(parent.shapingwavelet!=NULL)
		shapingwavelet=new ComplexArray(*(parent.shapingwavelet));
	else
		shapingwavelet=NULL;
	if(parent.invshapingwavelet!=NULL)
		invshapingwavelet=new ComplexArray(*(parent.invshapingwavelet));
	else
		invshapingwavelet=NULL;
	itermax=parent.itermax;
	tol=parent.tol;
	tolinstep=parent.tolinstep;
	type=parent.type;
	nfft=parent.nfft;
	sample_shift=parent.sample_shift;
	result=parent.result;
	icount=parent.icount;
	eps=parent.eps;
	if(type!=XCORR)
		damp=parent.damp;
	if(type==MULTI_TAPER)
	{
		nw=parent.nw;
		seql=parent.seql;
		taperlen=parent.taperlen;
		if(parent.tapers!=NULL)
		{
			int sequ=nw*2-2;
			int nseq=sequ-seql+1;
			tapers=new double[nseq*taperlen];
			cblas_dcopy(nseq*taperlen,parent.tapers,1,tapers,1);
		}
		else
			tapers=NULL;
		noise=parent.noise;
	}
	else
	{
		noise=NULL;
		tapers=NULL;
	}
}
int SimpleGeneralIterDecon::read_metadata(Metadata &md)
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	string itertype=md.get_string("iterative_decon_type");
	if(itertype=="multi_taper")
		type=MULTI_TAPER;
	else if(itertype=="water_level")
		type=WATER_LEVEL;
	else if(itertype=="least_square")
		type=LEAST_SQ;
	else if(itertype=="xcorr")
		type=XCORR;
	else
	{
		cout<<"WARNING: iterative_decon_type not defined, set to xcorr by default."<<endl;
		type=XCORR;
	}
	try{
		itermax=md.get_int("max_iteration_time");
		tol=md.get_double("tolerance_of_misfit");
		tolinstep=md.get_double("tolerance_of_misfit_in_step");
		nfft=md.get_int("operator_nfft");
		sample_shift=md.get_int("sample_shift");
		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
		if(type!=XCORR)
			damp=md.get_double("damping_factor");
		if(type==MULTI_TAPER)
		{
			nw=md.get_double("time_bandwidth_product");
			nw=static_cast<int>(nw*2)/2.0;
			seql=md.get_int("lower_dpss");
			if(seql>2*nw-2)
				seql=static_cast<int>(2*nw-2);
			taperlen=md.get_int("taper_length");
			int sequ=nw*2-2;
			int nseq=sequ-seql+1;
			if(tapers!=NULL)
				delete [] tapers;
			tapers=new double[nseq*taperlen];
			dpss_calc(taperlen, nw, seql, sequ, tapers);
		}
	}catch(SeisppError& err)
	{
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleGeneralIterDecon::read_metadata"<<endl;
		exit(-1);
	}
	try{
		string wavelettype=md.get_string("shaping_wavelet_type");
		if(wavelettype=="gaussian")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			double *r=gaussian(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
			if(type!=XCORR)
			{
				float invfpeak=md.get_double("shaping_wavelet_frequency_for_inverse");
				double *ir=gaussian(invfpeak,dt,nfft);
				if(invshapingwavelet!=NULL)
					delete invshapingwavelet;
				invshapingwavelet=new ComplexArray(nfft,ir);
				gsl_fft_complex_forward(invshapingwavelet->ptr(), 1, nfft, wavetable, workspace);
				delete [] ir;
			}
		}
		else if(wavelettype=="ricker")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=rickerwavelet(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
			if(type!=XCORR)
			{
				float invfpeak=md.get_double("shaping_wavelet_frequency_for_inverse");
				double *ir=rickerwavelet(invfpeak,dt,nfft);
				if(invshapingwavelet!=NULL)
					delete invshapingwavelet;
				invshapingwavelet=new ComplexArray(nfft,ir);
				gsl_fft_complex_forward(invshapingwavelet->ptr(), 1, nfft, wavetable, workspace);
				delete [] ir;
			}
		}
		else
		{
			cout<<"WARNING: Using the raw output of deconvolution"<<endl;
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,1.0);
			if(type!=XCORR)
			{
				float invfpeak=md.get_double("shaping_wavelet_frequency_for_inverse");
				float dt=md.get_double("shaping_wavelet_dt");
				double *ir=rickerwavelet(invfpeak,dt,nfft);
				if(invshapingwavelet!=NULL)
					delete invshapingwavelet;
				invshapingwavelet=new ComplexArray(nfft,ir);
				gsl_fft_complex_forward(invshapingwavelet->ptr(), 1, nfft, wavetable, workspace);
				delete [] ir;
			}
		}
	}catch(SeisppError& err)
	{
		if(type!=XCORR)
		{
			throw err;
			return 1;
		}
		cerr<<err.what()<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}catch(...)
	{
		if(type!=XCORR)
		{
			throw SeisppError("Error from SimpleGeneralIterDecon::read_metadata");
			return 1;
		}
		cerr<<"Unknown error from SimpleGeneralIterDecon::read_metadata"<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
SimpleGeneralIterDecon::SimpleGeneralIterDecon(Metadata &md)
{
	data=NULL;
	wavelet=NULL;
	shapingwavelet=NULL;
	invshapingwavelet=NULL;
	noise=NULL;
	tapers=NULL;
	icount=-1;
	eps=-1.0;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleGeneralIterDecon constructor"<<endl;
		exit(-1);
	}
}
void SimpleGeneralIterDecon::changeparameter(Metadata &md)
{
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot change the decon parameters ";
		throw err;
	}catch(...)
	{
		throw SeisppError("Error at SimpleGeneralIterDecon::changeparameter");
	}
}
int SimpleGeneralIterDecon::loadnoise(vector<double> &n)
{
	noise=&n;
	return 0;
}
SimpleGeneralIterDecon::SimpleGeneralIterDecon(Metadata &md,vector<double> &n,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	invshapingwavelet=NULL;
	tapers=NULL;
	icount=-1;
	eps=-1.0;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleGeneralIterDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	if(type==MULTI_TAPER)
		noise=&n;
	else
		noise=NULL;
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
SimpleGeneralIterDecon::SimpleGeneralIterDecon(Metadata &md,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	invshapingwavelet=NULL;
	tapers=NULL;
	noise=NULL;
	icount=-1;
	eps=-1.0;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from SimpleGeneralIterDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	if(type==MULTI_TAPER)
	{
		cout<<"WARNING: no noise read in for decon, set decon type to least_square"<<endl;
		type=LEAST_SQ;
	}
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
int SimpleGeneralIterDecon::apply()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	ComplexArray b_fft(nfft,*wavelet),n_fft(nfft),d;
	if(type==MULTI_TAPER && noise!=NULL)
	{
		int sequ=nw*2-2;
		int nseq=sequ-seql+1;
				
		//apply fft to the input trace data
		ComplexArray od_fft(nfft,*data);
		gsl_fft_complex_forward(od_fft.ptr(), 1, nfft, wavetable, workspace);
//        od_fft=od_fft**shapingwavelet;
        ComplexArray od(nfft,*data);
//        gsl_fft_complex_inverse(od.ptr(), 1, nfft, wavetable, workspace);
		//apply tapers to the data
		ComplexArray tdata(nfft*nseq);
		ComplexArray d_fft(nfft);
		for(int i=0;i<nseq;i++)
		{
			for(int k=0;k<taperlen;k++)
			{
				*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]**od.ptr(k);
			}
			gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
		}
		for(int i=0;i<nfft;i++)
		{
			double sumr=0,sumi=0;
			for(int k=0;k<nseq;k++)
			{
				sumr+=tdata[k*nfft+i].real();
				sumi+=tdata[k*nfft+i].imag();
			}
			*d_fft.ptr(i)=sumr/nseq;
			*(d_fft.ptr(i)+1)=sumi/nseq;
		}
		double odmax=abs(od_fft[0]),dmax=abs(d_fft[0]);
		for(int i=1;i<nfft;i++)
		{
			if(abs(od_fft[i])>odmax)
				odmax=abs(od_fft[i]);
			if(abs(d_fft[i])>dmax)
				dmax=abs(d_fft[i]);
		}
		d_fft=d_fft*(odmax/dmax);
		d=d_fft;
		gsl_fft_complex_inverse(d.ptr(), 1, nfft, wavetable, workspace);
	
		//apply fft to beam
		ComplexArray ob_fft(nfft,*wavelet);
		gsl_fft_complex_forward(ob_fft.ptr(), 1, nfft, wavetable, workspace);
		//apply tapers to the beam
		int waveletlen;
		if(wavelet->size()<taperlen)
			waveletlen=wavelet->size();
		else
			waveletlen=taperlen;
		for(int i=0;i<nseq;i++)
		{
			for(int k=0;k<nfft;k++)
			{
				if(k<waveletlen)
					*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*wavelet->at(k);
				else
					*tdata.ptr(i*nfft+k)=0.0;
				*(tdata.ptr(i*nfft+k)+1)=0.0;
			}
			gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
		}
		for(int i=0;i<nfft;i++)
		{
			double sumr=0,sumi=0;
			for(int k=0;k<nseq;k++)
			{
				sumr+=tdata[k*nfft+i].real();
				sumi+=tdata[k*nfft+i].imag();
			}
			*b_fft.ptr(i)=sumr/nseq;
			*(b_fft.ptr(i)+1)=sumi/nseq;
		}
		double obmax=abs(ob_fft[0]),bmax=abs(b_fft[0]);
		for(int i=1;i<nfft;i++)
		{
			if(abs(ob_fft[i])>obmax)
				obmax=abs(ob_fft[i]);
			if(abs(b_fft[i])>bmax)
				bmax=abs(b_fft[i]);
		}
		b_fft=b_fft*(obmax/bmax);
	
		//apply fft to noise
		ComplexArray on_fft(nfft,*noise);
		gsl_fft_complex_forward(on_fft.ptr(), 1, nfft, wavetable, workspace);
		//apply tapers to the noise
		int noiselen;
		if(noise->size()<taperlen)
			noiselen=noise->size();
		else
			noiselen=taperlen;
		for(int i=0;i<nseq;i++)
		{
			for(int k=0;k<nfft;k++)
			{
				if(k<noiselen)
					*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*noise->at(k);
				else
					*tdata.ptr(i*nfft+k)=0.0;
				*(tdata.ptr(i*nfft+k)+1)=0.0;
			}
			gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
		}
		for(int i=0;i<nfft;i++)
		{
			double sumr=0,sumi=0;
			for(int k=0;k<nseq;k++)
			{
				sumr+=tdata[k*nfft+i].real();
				sumi+=tdata[k*nfft+i].imag();
			}
			*n_fft.ptr(i)=sumr/nseq;
			*(n_fft.ptr(i)+1)=sumi/nseq;
		}
		double onmax=abs(on_fft[0]),nmax=abs(n_fft[0]);
		for(int i=1;i<nfft;i++)
		{
			if(abs(on_fft[i])>onmax)
				onmax=abs(on_fft[i]);
			if(abs(n_fft[i])>nmax)
				nmax=abs(n_fft[i]);
		}
		n_fft=n_fft*(onmax/nmax);
	}
	else
	{
		//convert input trace data to ComplexArray
		ComplexArray d_fft(nfft,*data);
		
		//gsl_fft_complex_forward(d_fft.ptr(), 1, nfft, wavetable, workspace);
		//d_fft=gaussianFilter(d_fft,2.5,0.1);
		//gsl_fft_complex_inverse(d_fft.ptr(), 1, nfft, wavetable, workspace);
		
		d=d_fft;

		//apply fft to beam
		gsl_fft_complex_forward(b_fft.ptr(), 1, nfft, wavetable, workspace);
		
		//b_fft=gaussianFilter(b_fft,2.5,0.1);
	}
	
	//construct inverse wavelet
    double b_rms=b_fft.rms();
	ComplexArray iw;
	if(type==XCORR)
	{
		b_fft.conj();
		iw=b_fft;
		b_fft.conj();
	}
	else if(type==LEAST_SQ)
	{
		b_fft.conj();
		ComplexArray conj_b_fft(b_fft);
		b_fft.conj();
		iw=(conj_b_fft**invshapingwavelet)/(conj_b_fft*b_fft+b_rms*damp);
	}
	else if(type==WATER_LEVEL)
	{
		//water level
		for(int i=0;i<nfft;i++)
		{
			if(abs(b_fft[i])<b_rms*damp)
			{
				*b_fft.ptr(i)=*b_fft.ptr(i)/abs(b_fft[i])*b_rms*damp;
				*(b_fft.ptr(i)+1)=*(b_fft.ptr(i)+1)/abs(b_fft[i])*b_rms*damp;
			}
		}
        
		iw=*invshapingwavelet/b_fft;
	}
	else if(type==MULTI_TAPER && noise!=NULL)
	{
		b_fft.conj();
		ComplexArray conj_b_fft(b_fft);
		b_fft.conj();
		n_fft=n_fft*damp;
		iw=(conj_b_fft**invshapingwavelet)/(conj_b_fft*b_fft+n_fft.abs());
	}
	else if(type==MULTI_TAPER && noise==NULL)
	{
		cout<<"WARNING: no noise read in for decon, set decon type to least_square"<<endl;
		b_fft.conj();
		ComplexArray conj_b_fft(b_fft);
		b_fft.conj();
		iw=(conj_b_fft**invshapingwavelet)/(conj_b_fft*b_fft+b_rms*damp);
	}
	else
	{
		throw SeisppError("Unexpected error at SimpleGeneralIterDecon::apply");
		return 1;
	}

	//deconvolution
	ComplexArray s(d);
	ComplexArray refl(nfft,0.0);
	double *z=new double[nfft];
	cblas_dcopy(nfft,d.ptr(),2,z,1);
	double dpower=power(z,nfft);
	double prevpower = dpower;
	double improvement=100;
	int gind=-1;
	vector<int> fakegind;
	int i;
	icount=0;
	for(i=0;i<itermax;i++)
	{
		ComplexArray originrefl(refl);
		ComplexArray origins(s);
		ComplexArray g;
		ComplexArray s_fft(s);
		gsl_fft_complex_forward(s_fft.ptr(), 1, nfft, wavetable, workspace);
		g=iw*s_fft;
		gsl_fft_complex_inverse(g.ptr(), 1, nfft, wavetable, workspace);
		
		for(int j=0;j<fakegind.size();j++)
		{
			g[fakegind[j]]=0;
		}
		if(gind!=-1)
		{
			g[gind]=0;
		}
		
		double gmax=g[0].real();
		gind=0;
		for(int j=0;j<nfft-sample_shift;j++)
		{
			if(abs(g[j].real())>abs(gmax))
			{
				gmax=g[j].real();
				gind=j;
			}
		}
		ComplexArray temp(nfft,0.0);
		*temp.ptr(gind)=gmax;
		ComplexArray temp_fft(temp);
		gsl_fft_complex_forward(temp_fft.ptr(), 1, nfft, wavetable, workspace);
		//apply shaping wavelet for output
//		temp_fft=temp_fft**shapingwavelet;
//		temp=temp_fft;
//		gsl_fft_complex_inverse(temp.ptr(), 1, nfft, wavetable, workspace);
		ComplexArray adtemp;
		adtemp=b_fft*temp_fft;
		gsl_fft_complex_inverse(adtemp.ptr(), 1, nfft, wavetable, workspace);
		double nume=0.0,denome=0.0;
		for(int j=0;j<nfft;j++)
		{
			nume+=adtemp[j].real()*s[j].real();
			denome+=adtemp[j].real()*adtemp[j].real();
		}
		//have to set imag part of temp to zero here?
//		for(int j=0;j<nfft;j++)
//		{
//			*temp.ptr(j)=temp[j].real()*nume/denome;
//			*(temp.ptr(j)+1)=0.0;
//		}
		*temp.ptr(gind)=temp[gind].real()*nume/denome;
		refl+=temp;
		ComplexArray refl_fft(refl);
		gsl_fft_complex_forward(refl_fft.ptr(), 1, nfft, wavetable, workspace);
		s=b_fft*refl_fft;
		gsl_fft_complex_inverse(s.ptr(), 1, nfft, wavetable, workspace);
		for(int j=0;j<nfft;j++)
		{
			*(s.ptr(j)+1)=0.0;
		}
		s=d-s;
		cblas_dcopy(nfft,s.ptr(),2,z,1);
		double spower=power(z,nfft);
		improvement = 100.0*(prevpower-spower)/dpower;
		eps=100*spower/dpower;
		if(eps<=tol)
		{
			break;
		}
		if(improvement<=tolinstep)
		{
			refl=originrefl;
			s=origins;
			fakegind.push_back(gind);
			icount--;
		}
		else if(i%(itermax/5)==0)
		{
			fakegind.clear();
		}
		prevpower = spower;
		icount++;
	}
	gsl_fft_complex_forward(refl.ptr(), 1, nfft, wavetable, workspace);
	refl=refl**shapingwavelet;
	gsl_fft_complex_inverse(refl.ptr(), 1, nfft, wavetable, workspace);
	cout<<"iter time: "<<icount<<endl;

	//save result
	if(sample_shift>0)
	{
		for(int k=sample_shift;k>0;k--)
			result.push_back(refl[nfft-k].real());
		for(int k=0;k<data->size()-sample_shift;k++)
			result.push_back(refl[k].real());
	}
	else
	{
		for(int k=0;k<data->size();k++)
			result.push_back(refl[k].real());
	}
	delete [] z;
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}


//ArrayLeastSquareDecon
ArrayLeastSquareDecon::ArrayLeastSquareDecon(const ArrayLeastSquareDecon &parent)
{
	data=parent.data;
	wavelet=parent.wavelet;
	if(parent.shapingwavelet!=NULL)
		shapingwavelet=new ComplexArray(*(parent.shapingwavelet));
	else
		shapingwavelet=NULL;
	damp=parent.damp;
	nfft=parent.nfft;
	sample_shift=parent.sample_shift;
	invbeam=parent.invbeam;
	result=parent.result;
}
int ArrayLeastSquareDecon::read_metadata(Metadata &md)
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	try{
		damp=md.get_double("damping_factor");
		nfft=md.get_int("operator_nfft");
		sample_shift=md.get_int("sample_shift");
		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
	}catch(SeisppError& err)
	{
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayLeastSquareDecon::read_metadata"<<endl;
		exit(-1);
	}
	try{
		string wavelettype=md.get_string("shaping_wavelet_type");
		if(wavelettype=="gaussian")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=gaussian(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else if(wavelettype=="ricker")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=rickerwavelet(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else
		{
			cout<<"WARNING: Using the raw output of deconvolution"<<endl;
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,1.0);
		}
	}catch(SeisppError& err)
	{
		cerr<<err.what()<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}catch(...)
	{
		cerr<<"Unknown error from ArrayLeastSquareDecon::read_metadata"<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
ArrayLeastSquareDecon::ArrayLeastSquareDecon(Metadata &md)
{
	data=NULL;
	wavelet=NULL;
	shapingwavelet=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayLeastSquareDecon constructor"<<endl;
		exit(-1);
	}
}
void ArrayLeastSquareDecon::changeparameter(Metadata &md)
{
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot change the decon parameters ";
		throw err;
	}catch(...)
	{
		throw SeisppError("Error at ArrayLeastSquareDecon::changeparameter");
	}
}
ArrayLeastSquareDecon::ArrayLeastSquareDecon(Metadata &md,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayLeastSquareDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	if(this->getinverse())
	{
		cerr<<"Failed in getting inverse beam"<<endl;
		throw SeisppError("Failed in getting inverse beam");
	}
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
int ArrayLeastSquareDecon::getinverse()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	//apply fft to wavelet
	ComplexArray b_fft(nfft,*wavelet);
	gsl_fft_complex_forward(b_fft.ptr(), 1, nfft, wavetable, workspace);
	
	//deconvolution: RF=conj(B)./(conj(B).*B+damp) 
	b_fft.conj();
	ComplexArray conj_b_fft(b_fft);
	b_fft.conj();
	
	double b_rms=b_fft.rms();
	invbeam=(conj_b_fft**shapingwavelet)/(conj_b_fft*b_fft+b_rms*damp);
	
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
int ArrayLeastSquareDecon::apply()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	//apply fft to the input trace data
	ComplexArray d_fft(nfft,*data);
	gsl_fft_complex_forward(d_fft.ptr(), 1, nfft, wavetable, workspace);

	ComplexArray rf_fft;
	
	rf_fft=invbeam*d_fft;
	
	//ifft gets result
	gsl_fft_complex_inverse(rf_fft.ptr(), 1, nfft, wavetable, workspace);
	if(sample_shift>0)
	{
		for(int k=sample_shift;k>0;k--)
			result.push_back(rf_fft[nfft-k].real());
		for(int k=0;k<data->size()-sample_shift;k++)
			result.push_back(rf_fft[k].real());
	}
	else
	{
		for(int k=0;k<data->size();k++)
			result.push_back(rf_fft[k].real());
	}
	
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
//ArrayWaterLevelDecon
ArrayWaterLevelDecon::ArrayWaterLevelDecon(const ArrayWaterLevelDecon &parent)
{
	data=parent.data;
	wavelet=parent.wavelet;
	if(parent.shapingwavelet!=NULL)
		shapingwavelet=new ComplexArray(*(parent.shapingwavelet));
	else
		shapingwavelet=NULL;
	wlv=parent.wlv;
	nfft=parent.nfft;
	sample_shift=parent.sample_shift;
	invbeam=parent.invbeam;
	result=parent.result;
}
int ArrayWaterLevelDecon::read_metadata(Metadata &md)
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	try{
		wlv=md.get_double("damping_factor");
		nfft=md.get_int("operator_nfft");
		sample_shift=md.get_int("sample_shift");
		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
	}catch(SeisppError& err)
	{
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayWaterLevelDecon::read_metadata"<<endl;
		exit(-1);
	}
	try{
		string wavelettype=md.get_string("shaping_wavelet_type");
		if(wavelettype=="gaussian")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=gaussian(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else if(wavelettype=="ricker")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=rickerwavelet(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else
		{
			cout<<"WARNING: Using the raw output of deconvolution"<<endl;
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,1.0);
		}
	}catch(SeisppError& err)
	{
		cerr<<err.what()<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}catch(...)
	{
		cerr<<"Unknown error from ArrayWaterLevelDecon::read_metadata"<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
ArrayWaterLevelDecon::ArrayWaterLevelDecon(Metadata &md)
{
	data=NULL;
	wavelet=NULL;
	shapingwavelet=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayWaterLevelDecon constructor"<<endl;
		exit(-1);
	}
}
void ArrayWaterLevelDecon::changeparameter(Metadata &md)
{
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot change the decon parameters ";
		throw err;
	}catch(...)
	{
		throw SeisppError("Error at ArrayWaterLevelDecon::changeparameter");
	}
}
ArrayWaterLevelDecon::ArrayWaterLevelDecon(Metadata &md,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayWaterLevelDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	if(this->getinverse())
	{
		cerr<<"Failed in getting inverse beam"<<endl;
		throw SeisppError("Failed in getting inverse beam");
	}
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
int ArrayWaterLevelDecon::getinverse()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	//apply fft to wavelet
	ComplexArray b_fft(nfft,*wavelet);
	gsl_fft_complex_forward(b_fft.ptr(), 1, nfft, wavetable, workspace);
	
	double b_rms=b_fft.rms();
    
	//water level
	for(int i=0;i<nfft;i++)
	{
		if(abs(b_fft[i])<b_rms*wlv)
		{
			*b_fft.ptr(i)=*b_fft.ptr(i)/abs(b_fft[i])*b_rms*wlv;
			*(b_fft.ptr(i)+1)=*(b_fft.ptr(i)+1)/abs(b_fft[i])*b_rms*wlv;
		}
	}
    
	//deconvolution: RF=1/(B+wlv)
	invbeam=*shapingwavelet/b_fft;
	
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
int ArrayWaterLevelDecon::apply()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	//apply fft to the input trace data
	ComplexArray d_fft(nfft,*data);
	gsl_fft_complex_forward(d_fft.ptr(), 1, nfft, wavetable, workspace);
	
	//deconvolution: RF=D./(B+wlv)
	ComplexArray rf_fft;
	rf_fft=invbeam*d_fft;
	
	//ifft gets result
	gsl_fft_complex_inverse(rf_fft.ptr(), 1, nfft, wavetable, workspace);
	if(sample_shift>0)
	{
		for(int k=sample_shift;k>0;k--)
			result.push_back(rf_fft[nfft-k].real());
		for(int k=0;k<data->size()-sample_shift;k++)
			result.push_back(rf_fft[k].real());
	}
	else
	{
		for(int k=0;k<data->size();k++)
			result.push_back(rf_fft[k].real());
	}
	
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
//ArrayMultiTaperDecon
ArrayMultiTaperDecon::ArrayMultiTaperDecon(const ArrayMultiTaperDecon &parent)
{
	data=parent.data;
	wavelet=parent.wavelet;
	noise=parent.noise;
	if(parent.shapingwavelet!=NULL)
		shapingwavelet=new ComplexArray(*(parent.shapingwavelet));
	else
		shapingwavelet=NULL;
	nw=parent.nw;
	seql=parent.seql;
	taperlen=parent.taperlen;
	if(parent.tapers!=NULL)
	{
		int sequ=nw*2-2;
		int nseq=sequ-seql+1;
		tapers=new double[nseq*taperlen];
		cblas_dcopy(nseq*taperlen,parent.tapers,1,tapers,1);
	}
	else
		tapers=NULL;
	damp=parent.damp;
	nfft=parent.nfft;
	sample_shift=parent.sample_shift;
	invbeam=parent.invbeam;
	result=parent.result;
}
int ArrayMultiTaperDecon::read_metadata(Metadata &md)
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	try{
		damp=md.get_double("damping_factor");
		nfft=md.get_int("operator_nfft");
		sample_shift=md.get_int("sample_shift");
		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
		nw=md.get_double("time_bandwidth_product");
		nw=static_cast<int>(nw*2)/2.0;
		seql=md.get_int("lower_dpss");
		if(seql>2*nw-2)
			seql=static_cast<int>(2*nw-2);
		taperlen=md.get_int("taper_length");
		int sequ=nw*2-2;
		int nseq=sequ-seql+1;
		if(tapers!=NULL)
			delete [] tapers;
		tapers=new double[nseq*taperlen];
		dpss_calc(taperlen, nw, seql, sequ, tapers);
	}catch(LAPACK_ERROR &err)
	{
		cerr<<err.getmsg()<<endl;
		throw SeisppError("Failed in generating tapers");
	}catch(SeisppError& err)
	{
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayMultiTaperDecon::read_metadata"<<endl;
		exit(-1);
	}
	try{
		string wavelettype=md.get_string("shaping_wavelet_type");
		if(wavelettype=="gaussian")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=gaussian(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else if(wavelettype=="ricker")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=rickerwavelet(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
		}
		else
		{
			cout<<"WARNING: Using the raw output of deconvolution"<<endl;
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,1.0);
		}
	}catch(SeisppError& err)
	{
		cerr<<err.what()<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}catch(...)
	{
		cerr<<"Unknown error from ArrayMultiTaperDecon::read_metadata"<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
ArrayMultiTaperDecon::ArrayMultiTaperDecon(Metadata &md)
{
	data=NULL;
	wavelet=NULL;
	shapingwavelet=NULL;
	noise=NULL;
	tapers=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayMultiTaperDecon constructor"<<endl;
		exit(-1);
	}
}
void ArrayMultiTaperDecon::changeparameter(Metadata &md)
{
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot change the decon parameters ";
		throw err;
	}catch(...)
	{
		throw SeisppError("Error at ArrayMultiTaperDecon::changeparameter");
	}
}
int ArrayMultiTaperDecon::loadnoise(vector<double> &n)
{
	noise=&n;
	return 0;
}
ArrayMultiTaperDecon::ArrayMultiTaperDecon(Metadata &md,vector<double> &n,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	tapers=NULL;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayMultiTaperDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	noise=&n;
	if(this->getinverse())
	{
		cerr<<"Failed in getting inverse beam"<<endl;
		throw SeisppError("Failed in getting inverse beam");
	}
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
int ArrayMultiTaperDecon::getinverse()
{
	if(noise==NULL)
	{
		cerr<<"ArrayMultiTaperDecon::getinverse: noise data is empty."<<endl;
		throw SeisppError("ArrayMultiTaperDecon::getinverse: noise data is empty.");
		return 1;
	}
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	int sequ=nw*2-2;
	int nseq=sequ-seql+1;
	//apply fft to beam
	ComplexArray tdata(nfft*nseq);
	ComplexArray ob_fft(nfft,*wavelet);
	gsl_fft_complex_forward(ob_fft.ptr(), 1, nfft, wavetable, workspace);
	//apply tapers to the beam
	ComplexArray b_fft(nfft);
	int waveletlen;
	if(wavelet->size()<taperlen)
		waveletlen=wavelet->size();
	else
		waveletlen=taperlen;
	for(int i=0;i<nseq;i++)
	{
		for(int k=0;k<waveletlen;k++)
		{
			*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*wavelet->at(k);
		}
		gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
	}
	for(int i=0;i<nfft;i++)
	{
		double sumr=0,sumi=0;
		for(int k=0;k<nseq;k++)
		{
			sumr+=tdata[k*nfft+i].real();
			sumi+=tdata[k*nfft+i].imag();
		}
		*b_fft.ptr(i)=sumr/nseq;
		*(b_fft.ptr(i)+1)=sumi/nseq;
	}
	double obmax=abs(ob_fft[0]),bmax=abs(b_fft[0]);
	for(int i=1;i<nfft;i++)
	{
		if(abs(ob_fft[i])>obmax)
			obmax=abs(ob_fft[i]);
		if(abs(b_fft[i])>bmax)
			bmax=abs(b_fft[i]);
	}
	b_fft=b_fft*(obmax/bmax);
	
	//apply fft to noise
	ComplexArray on_fft(nfft,*noise);
	gsl_fft_complex_forward(on_fft.ptr(), 1, nfft, wavetable, workspace);
	//apply tapers to the noise
	ComplexArray n_fft(nfft);
	int noiselen;
	if(noise->size()<taperlen)
		noiselen=noise->size();
	else
		noiselen=taperlen;
	for(int i=0;i<nseq;i++)
	{
		for(int k=0;k<nfft;k++)
		{
			if(k<noiselen)
				*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*noise->at(k);
			else
				*tdata.ptr(i*nfft+k)=0.0;
			*(tdata.ptr(i*nfft+k)+1)=0.0;
		}
		gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
	}
	for(int i=0;i<nfft;i++)
	{
		double sumr=0,sumi=0;
		for(int k=0;k<nseq;k++)
		{
			sumr+=tdata[k*nfft+i].real();
			sumi+=tdata[k*nfft+i].imag();
		}
		*n_fft.ptr(i)=sumr/nseq;
		*(n_fft.ptr(i)+1)=sumi/nseq;
	}
	double onmax=abs(on_fft[0]),nmax=abs(n_fft[0]);
	for(int i=1;i<nfft;i++)
	{
		if(abs(on_fft[i])>onmax)
			onmax=abs(on_fft[i]);
		if(abs(n_fft[i])>nmax)
			nmax=abs(n_fft[i]);
	}
	n_fft=n_fft*(onmax/nmax);
	
	//deconvolution: RF=conj(B)./(conj(B).*B+abs(S)) 
	b_fft.conj();
	ComplexArray conj_b_fft(b_fft);
	b_fft.conj();
	
	n_fft=n_fft*damp;
	
	invbeam=(conj_b_fft**shapingwavelet)/(conj_b_fft*b_fft+n_fft.abs());
	
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;

}
int ArrayMultiTaperDecon::apply()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	int sequ=nw*2-2;
	int nseq=sequ-seql+1;
	//apply fft to the input trace data
	ComplexArray od_fft(nfft,*data);
	gsl_fft_complex_forward(od_fft.ptr(), 1, nfft, wavetable, workspace);
	//apply tapers to the data
	ComplexArray tdata(nfft*nseq);
	ComplexArray d_fft(nfft);
	int datalen;
	if(data->size()<taperlen)
		datalen=data->size();
	else
		datalen=taperlen;
	for(int i=0;i<nseq;i++)
	{
		for(int k=0;k<datalen;k++)
		{
			*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*data->at(k);
		}
		gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
	}
	for(int i=0;i<nfft;i++)
	{
		double sumr=0,sumi=0;
		for(int k=0;k<nseq;k++)
		{
			sumr+=tdata[k*nfft+i].real();
			sumi+=tdata[k*nfft+i].imag();
		}
		*d_fft.ptr(i)=sumr/nseq;
		*(d_fft.ptr(i)+1)=sumi/nseq;
	}
	double odmax=abs(od_fft[0]),dmax=abs(d_fft[0]);
	for(int i=1;i<nfft;i++)
	{
		if(abs(od_fft[i])>odmax)
			odmax=abs(od_fft[i]);
		if(abs(d_fft[i])>dmax)
			dmax=abs(d_fft[i]);
	}
	d_fft=d_fft*(odmax/dmax);
	
	//deconvolution: RF=conj(B).*D./(conj(B).*B+abs(S))
	ComplexArray rf_fft;
	
	rf_fft=invbeam*d_fft;
	
	//ifft gets result
	gsl_fft_complex_inverse(rf_fft.ptr(), 1, nfft, wavetable, workspace);
	if(sample_shift>0)
	{
		for(int k=sample_shift;k>0;k--)
			result.push_back(rf_fft[nfft-k].real());
		for(int k=0;k<data->size()-sample_shift;k++)
			result.push_back(rf_fft[k].real());
	}
	else
	{
		for(int k=0;k<data->size();k++)
			result.push_back(rf_fft[k].real());
	}
	
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
//ArrayGeneralIterDecon
ArrayGeneralIterDecon::ArrayGeneralIterDecon(const ArrayGeneralIterDecon &parent)
{
	data=parent.data;
	wavelet=parent.wavelet;
	if(parent.shapingwavelet!=NULL)
		shapingwavelet=new ComplexArray(*(parent.shapingwavelet));
	else
		shapingwavelet=NULL;
	if(parent.invshapingwavelet!=NULL)
		invshapingwavelet=new ComplexArray(*(parent.invshapingwavelet));
	else
		invshapingwavelet=NULL;
	itermax=parent.itermax;
	tol=parent.tol;
	type=parent.type;
	tolinstep=parent.tolinstep;
	nfft=parent.nfft;
	sample_shift=parent.sample_shift;
	invbeam=parent.invbeam;
	result=parent.result;
	icount=parent.icount;
	eps=parent.eps;
	if(type!=XCORR)
		damp=parent.damp;
	if(type==MULTI_TAPER)
	{
		nw=parent.nw;
		seql=parent.seql;
		taperlen=parent.taperlen;
		if(parent.tapers!=NULL)
		{
			int sequ=nw*2-2;
			int nseq=sequ-seql+1;
			tapers=new double[nseq*taperlen];
			cblas_dcopy(nseq*taperlen,parent.tapers,1,tapers,1);
		}
		else
			tapers=NULL;
		noise=parent.noise;
	}
	else
	{
		noise=NULL;
		tapers=NULL;
	}
}
int ArrayGeneralIterDecon::read_metadata(Metadata &md)
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	string itertype=md.get_string("iterative_decon_type");
	if(itertype=="multi_taper")
		type=MULTI_TAPER;
	else if(itertype=="water_level")
		type=WATER_LEVEL;
	else if(itertype=="least_square")
		type=LEAST_SQ;
	else if(itertype=="xcorr")
		type=XCORR;
	else
	{
		cout<<"WARNING: iterative_decon_type not defined, set to xcorr by default."<<endl;
		type=XCORR;
	}
	try{
		itermax=md.get_int("max_iteration_time");
		tol=md.get_double("tolerance_of_misfit");
		tolinstep=md.get_double("tolerance_of_misfit_in_step");
		nfft=md.get_int("operator_nfft");
		sample_shift=md.get_int("sample_shift");
		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
		if(type!=XCORR)
			damp=md.get_double("damping_factor");
		if(type==MULTI_TAPER)
		{
			nw=md.get_double("time_bandwidth_product");
			nw=static_cast<int>(nw*2)/2.0;
			seql=md.get_int("lower_dpss");
			if(seql>2*nw-2)
				seql=static_cast<int>(2*nw-2);
			taperlen=md.get_int("taper_length");
			int sequ=nw*2-2;
			int nseq=sequ-seql+1;
			if(tapers!=NULL)
				delete [] tapers;
			tapers=new double[nseq*taperlen];
			dpss_calc(taperlen, nw, seql, sequ, tapers);
		}
	}catch(LAPACK_ERROR &err)
	{
		cerr<<err.getmsg()<<endl;
		throw SeisppError("Failed in generating tapers");
	}catch(SeisppError& err)
	{
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayGeneralIterDecon::read_metadata"<<endl;
		exit(-1);
	}
	try{
		string wavelettype=md.get_string("shaping_wavelet_type");
		if(wavelettype=="gaussian")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			double *r=gaussian(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
			if(type!=XCORR)
			{
				float invfpeak=md.get_double("shaping_wavelet_frequency_for_inverse");
				double *ir=gaussian(invfpeak,dt,nfft);
				if(invshapingwavelet!=NULL)
					delete invshapingwavelet;
				invshapingwavelet=new ComplexArray(nfft,ir);
				gsl_fft_complex_forward(invshapingwavelet->ptr(), 1, nfft, wavetable, workspace);
				delete [] ir;
			}
		}
		else if(wavelettype=="ricker")
		{
			float fpeak=md.get_double("shaping_wavelet_frequency");
			float dt=md.get_double("shaping_wavelet_dt");
			//construct wavelet and fft
			double *r=rickerwavelet(fpeak,dt,nfft);
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,r);
			gsl_fft_complex_forward(shapingwavelet->ptr(), 1, nfft, wavetable, workspace);
			delete [] r;
			if(type!=XCORR)
			{
				float invfpeak=md.get_double("shaping_wavelet_frequency_for_inverse");
				double *ir=rickerwavelet(invfpeak,dt,nfft);
				if(invshapingwavelet!=NULL)
					delete invshapingwavelet;
				invshapingwavelet=new ComplexArray(nfft,ir);
				gsl_fft_complex_forward(invshapingwavelet->ptr(), 1, nfft, wavetable, workspace);
				delete [] ir;
			}
		}
		else
		{
			cout<<"WARNING: Using the raw output of deconvolution"<<endl;
			if(shapingwavelet!=NULL)
				delete shapingwavelet;
			shapingwavelet=new ComplexArray(nfft,1.0);
			if(type!=XCORR)
			{
				float invfpeak=md.get_double("shaping_wavelet_frequency_for_inverse");
				float dt=md.get_double("shaping_wavelet_dt");
				double *ir=rickerwavelet(invfpeak,dt,nfft);
				if(invshapingwavelet!=NULL)
					delete invshapingwavelet;
				invshapingwavelet=new ComplexArray(nfft,ir);
				gsl_fft_complex_forward(invshapingwavelet->ptr(), 1, nfft, wavetable, workspace);
				delete [] ir;
			}
		}
	}catch(SeisppError& err)
	{
		if(type!=XCORR)
		{
			throw err;
			return 1;
		}
		cerr<<err.what()<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}catch(...)
	{
		if(type!=XCORR)
		{
			throw SeisppError("Error from ArrayGeneralIterDecon::read_metadata");
			return 1;
		}
		cerr<<"Unknown error from ArrayGeneralIterDecon::read_metadata"<<endl<<"set wavelet_type to spike"<<endl;
		cout<<"WARNING: Using the raw output of deconvolution"<<endl;
		if(shapingwavelet!=NULL)
			delete shapingwavelet;
		shapingwavelet=new ComplexArray(nfft,1.0);
	}
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
ArrayGeneralIterDecon::ArrayGeneralIterDecon(Metadata &md)
{
	data=NULL;
	wavelet=NULL;
	shapingwavelet=NULL;
	invshapingwavelet=NULL;
	noise=NULL;
	tapers=NULL;
	icount=-1;
	eps=-1.0;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayGeneralIterDecon constructor"<<endl;
		exit(-1);
	}
}
void ArrayGeneralIterDecon::changeparameter(Metadata &md)
{
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot change the decon parameters ";
		throw err;
	}catch(...)
	{
		throw SeisppError("Error at ArrayGeneralIterDecon::changeparameter");
	}
}
int ArrayGeneralIterDecon::loadnoise(vector<double> &n)
{
	noise=&n;
	return 0;
}
ArrayGeneralIterDecon::ArrayGeneralIterDecon(Metadata &md,vector<double> &n,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	invshapingwavelet=NULL;
	tapers=NULL;
	icount=-1;
	eps=-1.0;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayGeneralIterDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	if(type==MULTI_TAPER)
		noise=&n;
	else
		noise=NULL;
	if(this->getinverse())
	{
		cerr<<"Failed in getting inverse beam"<<endl;
		throw SeisppError("Failed in getting inverse beam");
	}
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
ArrayGeneralIterDecon::ArrayGeneralIterDecon(Metadata &md,vector<double> &w,vector<double> &d)
{
	shapingwavelet=NULL;
	invshapingwavelet=NULL;
	tapers=NULL;
	noise=NULL;
	icount=-1;
	eps=-1.0;
	try{
		this->read_metadata(md);
	}catch(SeisppError& err)
	{
		err.message+=" Cannot construct decon operator ";
		throw err;
	}catch(...)
	{
		cerr<<"Unknown error from ArrayGeneralIterDecon constructor"<<endl;
		exit(-1);
	}
	wavelet=&w;
	data=&d;
	if(type==MULTI_TAPER)
	{
		cout<<"WARNING: no noise read in for decon, set decon type to least_square"<<endl;
		type=LEAST_SQ;
	}
	if(this->getinverse())
	{
		cerr<<"Failed in getting inverse beam"<<endl;
		throw SeisppError("Failed in getting inverse beam");
	}
	if(this->apply())
	{
		cerr<<"Failed in doing deconvolution"<<endl;
		throw SeisppError("Failed in doing deconvolution");
	}
}
int ArrayGeneralIterDecon::getinverse()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);

	ComplexArray b_fft(nfft,*wavelet);
	ComplexArray n_fft(nfft);
	if(type==MULTI_TAPER && noise!=NULL)
	{
		int sequ=nw*2-2;
		int nseq=sequ-seql+1;
		//apply fft to beam
		ComplexArray tdata(nfft*nseq);
		ComplexArray ob_fft(nfft,*wavelet);
		gsl_fft_complex_forward(ob_fft.ptr(), 1, nfft, wavetable, workspace);
        //ob_fft=ob_fft**shapingwavelet;
        //ComplexArray ob(ob_fft);
        //gsl_fft_complex_inverse(ob.ptr(), 1, nfft, wavetable, workspace);
		//apply tapers to the beam
		int waveletlen;
		if(wavelet->size()<taperlen)
			waveletlen=wavelet->size();
		else
			waveletlen=taperlen;
		for(int i=0;i<nseq;i++)
		{
			for(int k=0;k<waveletlen;k++)
			{
				*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*wavelet->at(k);
			}
			gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
		}
		for(int i=0;i<nfft;i++)
		{
			double sumr=0,sumi=0;
			for(int k=0;k<nseq;k++)
			{
				sumr+=tdata[k*nfft+i].real();
				sumi+=tdata[k*nfft+i].imag();
			}
			*b_fft.ptr(i)=sumr/nseq;
			*(b_fft.ptr(i)+1)=sumi/nseq;
		}
		double obmax=abs(ob_fft[0]),bmax=abs(b_fft[0]);
		for(int i=1;i<nfft;i++)
		{
			if(abs(ob_fft[i])>obmax)
				obmax=abs(ob_fft[i]);
			if(abs(b_fft[i])>bmax)
				bmax=abs(b_fft[i]);
		}
		b_fft=b_fft*(obmax/bmax);
	
		//apply fft to noise
		ComplexArray on_fft(nfft,*noise);
		gsl_fft_complex_forward(on_fft.ptr(), 1, nfft, wavetable, workspace);
        //on_fft=on_fft**shapingwavelet;
        //ComplexArray on(on_fft);
        //gsl_fft_complex_inverse(on.ptr(), 1, nfft, wavetable, workspace);
		//apply tapers to the noise
		int noiselen;
		if(noise->size()<taperlen)
			noiselen=noise->size();
		else
			noiselen=taperlen;
		for(int i=0;i<nseq;i++)
		{
			for(int k=0;k<nfft;k++)
			{
				if(k<noiselen)
					*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*noise->at(k);
				else
					*tdata.ptr(i*nfft+k)=0.0;
				*(tdata.ptr(i*nfft+k)+1)=0.0;
			}
			gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
		}
		for(int i=0;i<nfft;i++)
		{
			double sumr=0,sumi=0;
			for(int k=0;k<nseq;k++)
			{
				sumr+=tdata[k*nfft+i].real();
				sumi+=tdata[k*nfft+i].imag();
			}
			*n_fft.ptr(i)=sumr/nseq;
			*(n_fft.ptr(i)+1)=sumi/nseq;
		}
		double onmax=abs(on_fft[0]),nmax=abs(n_fft[0]);
		for(int i=1;i<nfft;i++)
		{
			if(abs(on_fft[i])>onmax)
				onmax=abs(on_fft[i]);
			if(abs(n_fft[i])>nmax)
				nmax=abs(n_fft[i]);
		}
		n_fft=n_fft*(onmax/nmax);
	}
	else
	{
		//apply fft to beam
		gsl_fft_complex_forward(b_fft.ptr(), 1, nfft, wavetable, workspace);
        //b_fft=b_fft**shapingwavelet;
	}
	//construct inverse wavelet
	double b_rms=b_fft.rms();
	if(type==XCORR)
	{
		b_fft.conj();
		invbeam=b_fft;
	}
	else if(type==LEAST_SQ)
	{
		b_fft.conj();
		ComplexArray conj_b_fft(b_fft);
		b_fft.conj();
		invbeam=(conj_b_fft**invshapingwavelet)/(conj_b_fft*b_fft+b_rms*damp);
	}
	else if(type==WATER_LEVEL)
	{
		//water level
		for(int i=0;i<nfft;i++)
		{
			if(abs(b_fft[i])<b_rms*damp)
			{
				*b_fft.ptr(i)=*b_fft.ptr(i)/abs(b_fft[i])*b_rms*damp;
				*(b_fft.ptr(i)+1)=*(b_fft.ptr(i)+1)/abs(b_fft[i])*b_rms*damp;
			}
		}
		invbeam=*invshapingwavelet/b_fft;
	}
	else if(type==MULTI_TAPER && noise!=NULL)
	{
		b_fft.conj();
		ComplexArray conj_b_fft(b_fft);
		b_fft.conj();
		n_fft=n_fft*damp;
		invbeam=(conj_b_fft**invshapingwavelet)/(conj_b_fft*b_fft+n_fft.abs());
	}
	else if(type==MULTI_TAPER && noise==NULL)
	{
		cout<<"WARNING: no noise read in for decon, set decon type to least_square"<<endl;
		b_fft.conj();
		ComplexArray conj_b_fft(b_fft);
		b_fft.conj();
		invbeam=(conj_b_fft**invshapingwavelet)/(conj_b_fft*b_fft+b_rms*damp);
	}
	else
	{
		throw SeisppError("Unexpected error at ArrayGeneralIterDecon::apply");
		return 1;
	}

	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}
int ArrayGeneralIterDecon::apply()
{
	//allocate space for fft
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (nfft);
	workspace = gsl_fft_complex_workspace_alloc (nfft);
	
	//convert input trace data to ComplexArray
	ComplexArray d(nfft,*data);
//    gsl_fft_complex_forward(d.ptr(), 1, nfft, wavetable, workspace);
//    d=d**shapingwavelet;
//    gsl_fft_complex_inverse(d.ptr(), 1, nfft, wavetable, workspace);
	ComplexArray b_fft(nfft,*wavelet);
	if(type==MULTI_TAPER)
	{
		int sequ=nw*2-2;
		int nseq=sequ-seql+1;
		
		//apply fft to the input trace data
		ComplexArray od_fft(nfft,*data);
		gsl_fft_complex_forward(od_fft.ptr(), 1, nfft, wavetable, workspace);
//        od_fft=od_fft**shapingwavelet;
		ComplexArray od(nfft,*data);
//        gsl_fft_complex_inverse(od.ptr(), 1, nfft, wavetable, workspace);
		//apply tapers to the data
		ComplexArray tdata(nfft*nseq);
		ComplexArray d_fft(nfft);
		int datalen;
		if(data->size()<taperlen)
			datalen=data->size();
		else
			datalen=taperlen;
		for(int i=0;i<nseq;i++)
		{
			for(int k=0;k<datalen;k++)
			{
				*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]**od.ptr(k);
			}
			gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
		}
		for(int i=0;i<nfft;i++)
		{
			double sumr=0,sumi=0;
			for(int k=0;k<nseq;k++)
			{
				sumr+=tdata[k*nfft+i].real();
				sumi+=tdata[k*nfft+i].imag();
			}
			*d_fft.ptr(i)=sumr/nseq;
			*(d_fft.ptr(i)+1)=sumi/nseq;
		}
		double odmax=abs(od_fft[0]),dmax=abs(d_fft[0]);
		for(int i=1;i<nfft;i++)
		{
			if(abs(od_fft[i])>odmax)
				odmax=abs(od_fft[i]);
			if(abs(d_fft[i])>dmax)
				dmax=abs(d_fft[i]);
		}
		d_fft=d_fft*(odmax/dmax);
		d=d_fft;
		gsl_fft_complex_inverse(d.ptr(), 1, nfft, wavetable, workspace);
	
		//apply fft to beam
		ComplexArray ob_fft(nfft,*wavelet);
		gsl_fft_complex_forward(ob_fft.ptr(), 1, nfft, wavetable, workspace);
        //ob_fft=ob_fft**shapingwavelet;
        //ComplexArray ob(ob_fft);
        //gsl_fft_complex_inverse(ob.ptr(), 1, nfft, wavetable, workspace);
		//apply tapers to the beam
		int waveletlen;
		if(wavelet->size()<taperlen)
			waveletlen=wavelet->size();
		else
			waveletlen=taperlen;
		for(int i=0;i<nseq;i++)
		{
			for(int k=0;k<nfft;k++)
			{
				if(k<waveletlen)
					*tdata.ptr(i*nfft+k)=tapers[i*taperlen+k]*wavelet->at(k);
				else
					*tdata.ptr(i*nfft+k)=0.0;
				*(tdata.ptr(i*nfft+k)+1)=0.0;
			}
			gsl_fft_complex_forward(tdata.ptr(i*nfft), 1, nfft, wavetable, workspace);
		}
		for(int i=0;i<nfft;i++)
		{
			double sumr=0,sumi=0;
			for(int k=0;k<nseq;k++)
			{
				sumr+=tdata[k*nfft+i].real();
				sumi+=tdata[k*nfft+i].imag();
			}
			*b_fft.ptr(i)=sumr/nseq;
			*(b_fft.ptr(i)+1)=sumi/nseq;
		}
		double obmax=abs(ob_fft[0]),bmax=abs(b_fft[0]);
		for(int i=1;i<nfft;i++)
		{
			if(abs(ob_fft[i])>obmax)
				obmax=abs(ob_fft[i]);
			if(abs(b_fft[i])>bmax)
				bmax=abs(b_fft[i]);
		}
		b_fft=b_fft*(obmax/bmax);
	}
	else
	{
		//apply fft to beam
		gsl_fft_complex_forward(b_fft.ptr(), 1, nfft, wavetable, workspace);
        //b_fft=b_fft**shapingwavelet;
	}

	//deconvolution
	ComplexArray s(d);
	ComplexArray refl(nfft,0.0);
	double *z=new double[nfft];
	cblas_dcopy(nfft,d.ptr(),2,z,1);
	double dpower=power(z,nfft);
	double prevpower = dpower;
	double improvement=100;
	int gind=-1;
	vector<int> fakegind;
	int i;
	icount=0;
	for(i=0;i<itermax;i++)
	{
		ComplexArray originrefl(refl);
		ComplexArray origins(s);
		ComplexArray g;
		ComplexArray s_fft(s);
		gsl_fft_complex_forward(s_fft.ptr(), 1, nfft, wavetable, workspace);
		g=invbeam*s_fft;
		gsl_fft_complex_inverse(g.ptr(), 1, nfft, wavetable, workspace);
		
		for(int j=0;j<fakegind.size();j++)
		{
			g[fakegind[j]]=0;
		}
		if(gind!=-1)
		{
			g[gind]=0;
		}
		
		double gmax=g[0].real();
		gind=0;
		for(int j=0;j<nfft-sample_shift;j++)
		{
			if(abs(g[j].real())>abs(gmax))
			{
				gmax=g[j].real();
				gind=j;
			}
		}
		ComplexArray temp(nfft,0.0);
		*temp.ptr(gind)=gmax;
		ComplexArray temp_fft(temp);
		gsl_fft_complex_forward(temp_fft.ptr(), 1, nfft, wavetable, workspace);
		//apply shaping wavelet for output
//		temp_fft=temp_fft**shapingwavelet;
//		temp=temp_fft;
//		gsl_fft_complex_inverse(temp.ptr(), 1, nfft, wavetable, workspace);
		ComplexArray adtemp;
		adtemp=b_fft*temp_fft;
		gsl_fft_complex_inverse(adtemp.ptr(), 1, nfft, wavetable, workspace);
		double nume=0.0,denome=0.0;
		for(int j=0;j<nfft;j++)
		{
			nume+=adtemp[j].real()*s[j].real();
			denome+=adtemp[j].real()*adtemp[j].real();
		}
		//have to set imag part of temp to zero here?
//		for(int j=0;j<nfft;j++)
//		{
//			*temp.ptr(j)=temp[j].real()*nume/denome;
//			*(temp.ptr(j)+1)=0.0;
//		}
		*temp.ptr(gind)=temp[gind].real()*nume/denome;
		refl+=temp;
		ComplexArray refl_fft(refl);
		gsl_fft_complex_forward(refl_fft.ptr(), 1, nfft, wavetable, workspace);
		s=b_fft*refl_fft;
		gsl_fft_complex_inverse(s.ptr(), 1, nfft, wavetable, workspace);
		for(int j=0;j<nfft;j++)
		{
			*(s.ptr(j)+1)=0.0;
		}
		s=d-s;
		cblas_dcopy(nfft,s.ptr(),2,z,1);
		double spower=power(z,nfft);
		improvement = 100.0*(prevpower-spower)/dpower;
		eps=100*spower/dpower;
		if(eps<=tol)
		{
			break;
		}
		if(improvement<=tolinstep)
		{
			refl=originrefl;
			s=origins;
			fakegind.push_back(gind);
			icount--;
		}
		else if(i%(itermax/5)==0)
		{
			fakegind.clear();
		}
		prevpower = spower;
		icount++;
	}
	gsl_fft_complex_forward(refl.ptr(), 1, nfft, wavetable, workspace);
	refl=refl**shapingwavelet;
	gsl_fft_complex_inverse(refl.ptr(), 1, nfft, wavetable, workspace);
	cout<<"iter time: "<<icount<<endl;

	//save result
	if(sample_shift>0)
	{
		for(int k=sample_shift;k>0;k--)
			result.push_back(refl[nfft-k].real());
		for(int k=0;k<data->size()-sample_shift;k++)
			result.push_back(refl[k].real());
	}
	else
	{
		for(int k=0;k<data->size();k++)
			result.push_back(refl[k].real());
	}
	delete [] z;
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return 0;
}


double *rickerwavelet(float fpeak, float dt, int n)
{
	//modified from http://toto-share.com/2011/05/cc-ricker-wavelets-code/
	int i, k;
	int nw;
	int nc;
	double nw1, alpha, beta;
 
	nw1 = 2.2/fpeak/dt;
	nw = 2*floor(nw1/2)+1;
	nc = floor(nw/2);
	double *wricker=new double[n];

	for (i=0; i<nw; i++)
	{
		k = i+1;
		alpha = (nc-k+1)*fpeak*dt*M_PI;
		beta = pow(alpha, 2.0);
		wricker[i] = (1 - (beta*2)) * exp(-beta);
	}
	for (i=nw;i<n;i++)
		wricker[i]=0;
	double *ricker=new double[n];
	if(fpeak==0)
	{
		ricker[0]=1.0;
		for(int i=1;i<n;i++)
			ricker[i]=0.0;
	}
	else
	{
		for(i=0;i<nc;i++)
			ricker[n-nc+i]=wricker[i];
		for(i=nc;i<nw;i++)
			ricker[i-nc]=wricker[i];
		for(i=nw-nc;i<n-nc;i++)
			ricker[i]=0;
	}
	delete [] wricker;
    return(ricker);
}

double *gaussian(float sigma, float dt, int n)
{
	double total=dt*(n-1);
	double tover2=total/2.0;
	double *t=new double[n];
	if(n%2)
	{
		for(int i=0;i<=tover2/dt;i++)
			t[i]=i*dt;
		for(int i=tover2/dt+1;i<n;i++)
			t[i]=-tover2+(i-tover2/dt-1)*dt;
	}
	else
	{
		for(int i=0;i<=n/2;i++)
			t[i]=i*dt;
		for(int i=n/2+1;i<n;i++)
			t[i]=-(n-i)*dt;
	}
	double *gw=new double[n];
	if(sigma==0)
	{
		gw[0]=1.0;
		for(int i=1;i<n;i++)
			gw[i]=0.0;
	}
	else
	{
		for(int i=0;i<n;i++)
			gw[i]=exp(-(t[i]/sigma)*(t[i]/sigma));
	}
	delete [] t;
	return(gw);
}
double power(double* data,int len) 
{
        double power=0;
        for (int i=0; i<len; i++) 
        {
            power += data[i]*data[i];
        }
        return power;
}
/** Modified from SOD: convolve a function with a unit-area Gaussian filter.
     *   G(w) = exp(-w^2 / (4 a^2))
     *  The 1D gaussian is: f(x) = 1/(2*PI*sigma) e^(-x^2/(q * sigma^2))
     *  and the impluse response is: g(x) = 1/(2*PI)e^(-sigma^2 * u^2 / 2)
     *
     * If gwidthFactor is zero, does not filter.
     */
ComplexArray gaussianFilter(ComplexArray x,
                            float gwidthFactor,
                            float dt) 
{
        // gwidthFactor of zero means no filter
        if (gwidthFactor == 0) {
            return x;
        }
        //float[] forward = forwardFFT(x);
        double df = 1/(x.size() * dt);
        double d_omega = 2*M_PI*df;
        double gwidth = 4*gwidthFactor*gwidthFactor;
        double gauss;
        double omega;

        // Handle the nyquist frequency
        omega = M_PI/dt; // eliminate 2 / 2
        gauss = exp(-omega*omega / gwidth);
        *x.ptr(x.size()/2) *= gauss;
        *(x.ptr(x.size()/2)+1) *= gauss;

        for (int i=1; i<x.size()/2; i++) {
            omega = i*d_omega;
            gauss = exp(-omega*omega / gwidth);
            *x.ptr(i) *= gauss;
            *(x.ptr(i)+1) *= gauss;
            *x.ptr(x.size()-i) *= gauss;
            *(x.ptr(x.size()-i)+1) *= gauss;
        }
        return x;
    }
