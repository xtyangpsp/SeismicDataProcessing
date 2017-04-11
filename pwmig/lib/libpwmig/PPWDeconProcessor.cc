#include <algorithm>
#include "coords.h"
#include "TimeWindow.h"
#include "PPWDeconProcessor.h"
#include "VectorStatistics.h"
/* Needed only for debug */
#include "SeismicPlot.h"
using namespace std;
using namespace SEISPP;
namespace SEISPP{
PWIndexPosition::PWIndexPosition()
{
    iu=0;
    lag0=0;
    time0=0.0;
    amp=0.0;
    SNR=0.0;
    for(int k=0;k<3;++k) v[k]=0.0;
}
PWIndexPosition::PWIndexPosition(int iuin, int lag0in, 
	double t0in, double ampin, double SNRin, double *vin)
{
	iu=iuin;
	lag0=lag0in;
        time0=t0in;
	amp=ampin;
        SNR=SNRin;
	for(int k=0;k<3;++k)v[k]=vin[k];
}
PWIndexPosition::PWIndexPosition(const PWIndexPosition& parent)
{
	iu=parent.iu;
	lag0=parent.lag0;
	time0=parent.time0;
	amp=parent.amp;
        SNR=parent.SNR;
	for(int k=0;k<3;++k)v[k]=parent.v[k];
}
PWIndexPosition& PWIndexPosition::operator=(const PWIndexPosition& parent)
{
    if(this != &parent)
    {
	iu=parent.iu;
	lag0=parent.lag0;
	time0=parent.time0;
	amp=parent.amp;
        SNR=parent.SNR;
	for(int k=0;k<3;++k)v[k]=parent.v[k];
    }
    return *this;
}

/*! Primary constructor for this object.

Construction of this beast fills in data attributes through 
a combination of initialization and computation.  Be aware raw
data (draw) is simply copied.  This means it is important that the
caller guarantee sufficient padding is present in draw to allow for
all feasible lags of any plane wave in the deconvolution or downstream
problems are likely. 
*/
PWDataMember::PWDataMember(ThreeComponentSeismogram& draw,
	TimeSeries& winp,
	    SlownessVector& u0in,
	        vector<SlownessVector>& rayp,
		    double wt, double lat0, double lon0)
		: ThreeComponentSeismogram(draw), u0(u0in), wavelet(winp)
{
	weight=wt;
	try {
                /* We need these to compute moveout */
		lat=get_double("sta_lat");
		lon=get_double("sta_lon");
		elev=get_double("sta_elev");
                /* Assume lat,lon stored in Metadata always in degrees */
                lat=rad(lat);
                lon=rad(lon);
		/* Intentionally inline this.  There is a procedure
		in pwstack that does this, but this calculation is so
		simple it seems preferable here to just inline it */
		const double Rearth(6378.164);
		double dnorth,deast,delta,azimuth;
		dist(lat0,lon0,lat,lon,&delta, &azimuth);
		delta *= Rearth;
		deast=delta*sin(azimuth);
		dnorth=delta*cos(azimuth);
		/* Here we compute lag in samples for each slowness vector */
		int np=rayp.size();
		lag.reserve(np);
		int i;
		for(i=0;i<np;++i)
		{
			//double ux=u0.ux+rayp[i].ux;
			//double uy=u0.uy+rayp[i].uy;
			double ux=rayp[i].ux;
			double uy=rayp[i].uy;
			double tlag=ux*deast+uy*dnorth;
                        /* lags are negative to get moveout */
			lag.push_back(nint(-tlag/dt));
		}
	}
	catch (MetadataGetError mderr)
	{
		throw mderr;
	}
}

PWDataMember::PWDataMember(const PWDataMember& parent)
  : ThreeComponentSeismogram(dynamic_cast<const ThreeComponentSeismogram&>(parent)),
    u0(parent.u0), lag(parent.lag),wavelet(parent.wavelet)
{
	weight=parent.weight;
	lat=parent.lat;
	lon=parent.lon;
	elev=parent.elev;
}
/* Uses trick described at:  
http://www.fredosaurus.com/notes-cpp/oop-overloading/overloadassign2.html
also in a more coherent/expanded form in 
http://www.acm.org/crossroads/xrds1-4/ovp.html
This essentially invokes the operator= function of the class 
from which this is derived.   In this case this is a ThreeComponentSeismogram.
*/
PWDataMember& PWDataMember::operator=(const PWDataMember& parent)
{
	if(this != &parent)
	{
		this->ThreeComponentSeismogram::operator=(parent);
		/* Now the actual members of this class */
		weight=parent.weight;
		lat=parent.lat;
		lon=parent.lon;
		elev=parent.elev;
		u0=parent.u0;
		wavelet=parent.wavelet;
		lag=parent.lag;
	}
	return *this;
}
dvector PWDataMember::polarization(int iu,double t0)
{
    /* compute lag formula.  Intentionally do not check index errors for speed */
	int rlag0=this->sample_number(t0)-lag[iu]-wavelet.sample_number(0.0);	
        //DEBUG
        cout << "In polarization:  sn(t0)="<<this->sample_number(t0)
            <<" lag[iu]="<<lag[iu]
            <<" wavelet.sample_number(0.0)="<<wavelet.sample_number(0.0)<<endl;
        cout << "neighbors"<<endl;
        int id0=this->sample_number(t0)-lag[iu];
        for(int idbg=id0-10;idbg<=(id0+10); ++idbg)
        {
            cout << idbg <<" ";
            for(int k=0;k<3;++k) cout << u(k,idbg)<<" ";
            cout <<endl;
        }
	double *baseptr=u.get_address(0,rlag0);
	int i;
        dvector result(3);
        result.zero();
        for(i=0;i<3;++i,baseptr++)
        {
            /* Note 3 interval on baseptr because data are in component order */
            result(i)=ddot(wavelet.ns,&(wavelet.s[0]),1,baseptr,3);
        }
        return(result);
}

/* remove algorithm has evolved. This version takes the 3vector
   data a as input and computes a scale factor using dot products
   with the wavelet stored internally.  The method assumes the
   wavelet has been normalized to have L2 norm of unity or this
   will have the amplitude really messsed up */
int PWDataMember::remove(int iu, double t0, double *a)
{
	/* Note this method intentionally does NOT check for bounds
	for efficiency.  User must guarantee this to avoid chaos */
	int rlag0=this->sample_number(t0)-lag[iu]-wavelet.sample_number(0.0);	
	double *baseptr=u.get_address(0,rlag0);
        for(int i=0;i<3;++i,baseptr++)
        {
            double sclthis;
            sclthis=-(this->weight)*a[i]; // Each data member is weighted by weight
            daxpy(wavelet.ns,sclthis,&(wavelet.s[0]),1,baseptr,3);
        }
        return(rlag0);
}
PWStackMember::PWStackMember(ThreeComponentSeismogram& din,
	SlownessVector& sv, double noise_estimate) 
            : ThreeComponentSeismogram(din)
{
	try{
		slow=sv;
		/*We initialize amp vector to 0s to inialize the container.
		Note ns is inherited from BasicTimeSeries*/
		amp.reserve(ns);
		for(int i=0;i<ns;++i)amp.push_back(0.0);
		ampmax=0.0;
		lagmax=0;
                noise_rms=noise_estimate;
	} catch (...) {throw;};
}
PWStackMember::PWStackMember(const PWStackMember& parent) 
    : ThreeComponentSeismogram(parent)
{
	slow=parent.slow;
	/* These are the private attributes */
	ampmax=parent.ampmax;
	lagmax=parent.lagmax;
	amp=parent.amp;
        noise_rms=parent.noise_rms;
}
/* Use same trick as PWDataMember to copy ThreeComponentSeismogram
attributes */
PWStackMember& PWStackMember::operator=(const PWStackMember& parent)
{
	if(this != &parent)
	{
		this->ThreeComponentSeismogram::operator=(parent);
		slow=parent.slow;
		/* These are the private attributes */
		ampmax=parent.ampmax;
		lagmax=parent.lagmax;
		amp=parent.amp;
                noise_rms=parent.noise_rms;
	}
	return *this;
}
double PWStackMember::maxamp()
{
	double *baseptr=u.get_address(0,0);
	/* Always recompute the amplitudes before searching for
	the maximum */
	for(int i=0;i<ns;++i)
	{
		/* Intentionally do not use push_back in 
		setting amp components for efficiency. */
		amp[i]=0.0;
		for(int k=0;k<3;++k,baseptr++)
				amp[i]+=(*(baseptr))*(*(baseptr));
		amp[i]=sqrt(amp[i]);
	}
	/* This duplicates the PeakAmplitude methods in seispp but
	differences in the algorithm require duplication of small
	code fragment */
	vector<double>::iterator amptr;
	amptr=max_element(amp.begin(),amp.end());
        /* Use the STL distance algorithm to get position of lagmax */
	lagmax=distance(amp.begin(),amptr);
        ampmax=*amptr;
	return(*amptr);
}

PPWDeconProcessor::PPWDeconProcessor(ThreeComponentEnsemble& threecd,
        ThreeComponentEnsemble& noised,
            SlownessVector& u0in,
	        vector<SlownessVector>& ulist,
	            vector<TimeSeries>& w,
		        Metadata& md) 
		: u0(u0in), parent_data(threecd), noise_data(noised),
                    wavelet(w),slow(ulist)
                        , plot_handle(md)  // Debug only 
{
    try {
        /* First fetch any control parameters */
    	maxiteration=md.get_int("maximum_iterations");
        MinimumNoiseWindow=md.get_double("minimum_noise_window_length");
        SNRfloor=md.get_double("snr_convergence_floor");
        UseSNR=md.get_bool("use_snr");
        if(parent_data.member.size() != noise_data.member.size())
        {
            cerr << "PPWDeconProcessor::constructor: "
                <<"Warning:  data and noise ensemble sizes differ."<<endl
                <<"Data size="<<parent_data.member.size() 
                <<" and noise ensemble size="
                <<noise_data.member.size()<<endl;
        }
    	maxiteration=md.get_int("maximum_iterations");
	nu=ulist.size();
        int ndens=threecd.member.size();
	if(ndens<=0) throw SeisppError(string("PPWDeconProcessor::constructor:  ")
			+ "ensemble passed has no data.");
	dt_stack=parent_data.member[0].dt;
	aperture=md.get_double("pseudostation_array_aperture");
	cutoff=md.get_double("pseudostation_array_cutoff");
        snravgsize=md.get_int("snr_averaging_length");
        // We need to normalize the wavelets to have unit L2 norm
        vector<TimeSeries>::iterator wptr;
        for(wptr=wavelet.begin();wptr!=wavelet.end();++wptr)
        {
            double scale=dnrm2(wptr->ns,&(wptr->s[0]),1);
            dscal(wptr->ns,1.0/scale,&(wptr->s[0]),1);
            //DEBUG
            cout << "Variable bandwidth constructor debug:  scale="<<scale<<endl;
        }
	stalat.reserve(ndens);
	stalon.reserve(ndens);
	vector<ThreeComponentSeismogram>::iterator dptr;
	for(dptr=parent_data.member.begin();
		dptr!=parent_data.member.end();++dptr)
	{
		double lat,lon;
		lat=dptr->get_double("sta_lat");
		lon=dptr->get_double("sta_lon");
		stalat.push_back(rad(lat));
		stalon.push_back(rad(lon));
	}
	/* These are initialized to useless values.  They are set when
	the compute method is called */
	ns_stack=-1;
	t0_stack=0.0;
	start_offset.clear();
	stacksize=-1;
	current_stack.clear();
	d.clear();
	nd=-1;
    } catch (...) {throw;};
}
/* This is a near duplicate of the above.  There is, no doubt, a clever
way to utilize the above constructor as part of this one, but in my 
view it is better to suffer this duplication than it is to utilize 
a devious trick and make the code obscure.  HOWEVER, be sure in any 
revisions that errors found in the above propagate to here.  The only
difference here is we build the wavelet ensemble as copies of a single
master */
PPWDeconProcessor::PPWDeconProcessor(ThreeComponentEnsemble& threecd,
        ThreeComponentEnsemble& noised,
	    SlownessVector& u0in,
		vector<SlownessVector>& ulist,
		    TimeSeries& w,
			Metadata& md) 
		: u0(u0in), parent_data(threecd), noise_data(noised),
                    slow(ulist)
                        , plot_handle(md)  // Debug only 
{
    try {
        /* First fetch any control parameters */
    	maxiteration=md.get_int("maximum_iterations");
        MinimumNoiseWindow=md.get_double("minimum_noise_window_length");
        SNRfloor=md.get_double("snr_convergence_floor");
        UseSNR=md.get_bool("use_snr");
	nu=ulist.size();
        int ndens=threecd.member.size();
	if(ndens<=0) throw SeisppError(string("PPWDeconProcessor::constructor:  ")
			+ "ensemble passed has no data.");
	dt_stack=parent_data.member[0].dt;
	aperture=md.get_double("pseudostation_array_aperture");
	cutoff=md.get_double("pseudostation_array_cutoff");
        snravgsize=md.get_int("snr_averaging_length");
        // We need to normalize the wavelets to have unit L2 norm
        // Since w is passed by reference have to make copy
        TimeSeries wtmp(w);
        double wscale=dnrm2(wtmp.ns,&(wtmp.s[0]),1);
        dscal(wtmp.ns,1.0/wscale,&(wtmp.s[0]),1);
	stalat.reserve(ndens);
	stalon.reserve(ndens);
	vector<ThreeComponentSeismogram>::iterator dptr;
	stalat.reserve(ndens);
	stalon.reserve(ndens);
	wavelet.reserve(ndens);
	for(dptr=parent_data.member.begin();
		dptr!=parent_data.member.end();++dptr)
	{
		double lat,lon;
		lat=dptr->get_double("sta_lat");
		lon=dptr->get_double("sta_lon");
		stalat.push_back(rad(lat));
		stalon.push_back(rad(lon));
		/* This sets all wavelet members to a single waveform*/
		wavelet.push_back(wtmp);
	}
	/* These are initialized to useless values.  They are set when
	the compute method is called */
	ns_stack=-1;
	t0_stack=0.0;
	start_offset.clear();
	stacksize=-1;
	current_stack.clear();
	d.clear();
	nd=-1;
    } catch (...) {throw;};
}
PPWDeconProcessor::PPWDeconProcessor(const PPWDeconProcessor& parent)
{
	u0=parent.u0;
	maxiteration=parent.maxiteration;
	nd=parent.nd;
	nu=parent.nu;
	ns_stack=parent.ns_stack;
	t0_stack=parent.t0_stack;
	dt_stack=parent.dt_stack;
	aperture=parent.aperture;
	cutoff=parent.cutoff;
	stalat=parent.stalat;
	stalon=parent.stalon;
	start_offset=parent.start_offset;
	stacksize=parent.stacksize;
	d=parent.d;
	parent_data=parent.parent_data;
	wavelet=parent.wavelet;
	slow=parent.slow;
	current_stack=parent.current_stack;
	snravgsize=parent.snravgsize;
}
PPWDeconProcessor& PPWDeconProcessor::operator=(const PPWDeconProcessor& parent)
{
    if(this!=&parent)
    {
	u0=parent.u0;
	maxiteration=parent.maxiteration;
	nd=parent.nd;
	nu=parent.nu;
	ns_stack=parent.ns_stack;
	t0_stack=parent.t0_stack;
	dt_stack=parent.dt_stack;
	aperture=parent.aperture;
	cutoff=parent.cutoff;
	stalat=parent.stalat;
	stalon=parent.stalon;
	start_offset=parent.start_offset;
	stacksize=parent.stacksize;
	d=parent.d;
	parent_data=parent.parent_data;
	wavelet=parent.wavelet;
	slow=parent.slow;
	current_stack=parent.current_stack;
	snravgsize=parent.snravgsize;
    }
    return *this;
}
PWIndexPosition PPWDeconProcessor::maxamp()
{
	vector<PWStackMember>::iterator mptr;
	PWIndexPosition result;
	mptr=current_stack.begin();
	result.amp=mptr->maxamp();
        result.SNR=mptr->SNR();
	result.lag0=mptr->lag_at_max();
	result.iu=0;
	++mptr;
        double testamp,testsnr;
        /* Warning:  this loop assumes current_stack vector has more than
           one element */
	for(int i=1;mptr!=current_stack.end();++mptr,++i)
	{
		testamp=mptr->maxamp();
                testsnr=mptr->SNR();
//DEBUG
//cout << "In PPWDeconProcessor::maxamp:  maxamp for stack"<<i<<"="<<testamp<<endl<<" with SNR="<<testsnr<<endl;
                if(UseSNR)
                {
                    if(testsnr>result.SNR)
                    {
			result.amp=testamp;
                        result.SNR=testsnr;
			result.lag0=mptr->lag_at_max();
                        result.time0=mptr->time_at_max();
			result.iu=i;
                    }
                }
                else if(testamp>result.amp)
		{
		    result.amp=testamp;
                    result.SNR=testsnr;
	            result.lag0=mptr->lag_at_max();
                    result.time0=mptr->time_at_max();
		    result.iu=i;
		}
	}
	for(int k=0;k<3;++k) result.v[k]
	 		= current_stack[result.iu].u(k,result.lag0);
        return(result);
}
void PPWDeconProcessor::stack()
{
	vector<PWDataMember>::iterator dptr;
	vector<PWStackMember>::iterator sptr;
	vector<SlownessVector>::iterator uptr;
	vector<int>::iterator tlptr;
	double *drp,*srp;
	int ntosum,i,j,total_lag;
	for(i=0,uptr=slow.begin(),sptr=current_stack.begin();
		sptr!=current_stack.end();++i,++uptr,++sptr)
	{
		sptr->u.zero();
		srp=sptr->u.get_address(0,0);
		
		for(dptr=d.begin(),tlptr=start_offset.begin();
			dptr!=d.end();++dptr,++tlptr)
		{
			total_lag=dptr->lag[i]+(*tlptr);
                        srp=sptr->u.get_address(0,total_lag);
			drp=dptr->u.get_address(0,0);
                        ntosum=3*(dptr->ns);
//DEBUG
/*
if((total_lag+dptr->ns)>(sptr->ns) ){

  cerr << "Coding error:  "
    << "ntosum="<<ntosum<<" but stack ns ="<< sptr->ns<<endl;
  exit(-1);
}
*/
			for(j=0;j<ntosum;++j,srp++,drp++) *srp+=(*drp);
		}
//cout << "Stack member number "<<i<<" has peak amplitude="<<Linf(sptr->u)<<endl;
	}
}
/* Private method to appraise if solution has converged.  Intentionally 
made a method because it is not yet clear how to define convergence.
Stubbed initally as we'll just go on count for initial testing.*/
bool PPWDeconProcessor::converged(int current_count,
        PWIndexPosition& ip)
{
    if(current_count>=maxiteration)
        return(true);
    else if((ip.SNR)<SNRfloor)
        return(true);
    /* This is a test to stop when snr reaches a floor.  Tests show this
       is essential when the noise is colored by significant low velocity
       components */
    recent_snr.push_back(ip.SNR);
    if(current_count>=snravgsize)
    {
        recent_snr.pop_front();
        list<double>::iterator qptr;
        double snravg;
        for(snravg=0.0,qptr=recent_snr.begin();qptr!=recent_snr.end();++qptr)
        {
            snravg+=(*qptr);
        }
        snravg /= static_cast<double>(snravgsize);
        if(current_count==snravgsize) last_snravg=snravg;
        if(snravg>last_snravg) return(true);
        last_snravg=snravg;
    }
    return(false);
}
/* Note caller must test for empty ensemble indicating no data for
   this pseudostation point */
auto_ptr<ThreeComponentEnsemble> PPWDeconProcessor::compute(double pslat,
					double pslon)
{
	try {
		int i,j;
		/* First we have to form the pseudostation gather */
		d.clear();
		double sumwt;
		int npd=parent_data.member.size();
		for(i=0,sumwt=0.0;i<npd;++i)
		{
			double wt;
			double lat,lon;
			wt=pseudostation_weight(stalat[i],stalon[i],pslat,pslon,
				aperture, cutoff);
			if(wt>0.0) 
			{
			  d.push_back(PWDataMember(parent_data.member[i],
			  	wavelet[i],u0,slow,wt,pslat,pslon));
			  sumwt += wt;
			}
		}
                nd=d.size();
                if(nd==0)
                {
                    auto_ptr<ThreeComponentEnsemble> 
                            nodata(new ThreeComponentEnsemble());
                    return(nodata);
                }

                /* First normalize the weights so sum of wt is 1 */
                for(i=0;i<nd;++i)
                {
                    d[i].weight/=sumwt;
                }

                /* We normalize each PWDataMember by weights so we do not have 
                   to normalized the stack repeatedly. */
                for(i=0;i<nd;++i){
                    dscal(3*d[i].ns,d[i].weight,d[i].u.get_address(0,0),1);
                }

		/* This section does two things.  First, it computes a size
		that allows blind sums without error checking.  Second, it initializes
		the seismic data in each element of the container to all zeros. */

		int maxnt0;
		double maxtend;
		vector<int>::iterator lptr;
		vector<PWDataMember>::iterator dptr;
		for(i=0,dptr=d.begin();dptr!=d.end();++i,++dptr)
		{
		    int nt0tmp=dptr->sample_number(0.0);
		    double tendtmp=dptr->endtime();
		    for(j=0,lptr=dptr->lag.begin();
		    	lptr!=dptr->lag.end();++j,++lptr)
		    {
		    	int nt0test=nt0tmp - (*lptr);
			double tendtest=tendtmp+((*lptr)*dt_stack);
			if( (i==0) && (j==0) )
			{
			    maxnt0=nt0test;
			    maxtend=tendtmp;
			}
			else
			{
			    if(maxnt0<nt0test) 
			    	maxnt0=nt0test;
			    if(maxtend<tendtest)
			        maxtend=tendtest;
			}
		    }
		}
                /* Next we need to compute noise estimates.  These
                   are created here as a vector of numbers of length
                   equal to size of the vector of slowness (vectors).  
                   */
                vector<double> noise_levels
                        = compute_noise_estimates(pslat,pslon);
		/* We add npad on top and bottom to be sure we don't overflow 
		the bounds with computed stack indices.  Do NOT make this too 
		large as it adds an inefficieny. */
		const int npad(2);  
		t0_stack=-(maxnt0+npad)*dt_stack;
		ns_stack=nint((maxtend-t0_stack)/dt_stack);
		ns_stack+=npad;
		current_stack.clear();
		current_stack.reserve(nu);
		ThreeComponentSeismogram stack_pattern(ns_stack);
		stack_pattern.live=true;
		stack_pattern.t0=t0_stack;
		stack_pattern.dt=dt_stack;
		stack_pattern.tref=relative;
		// This is not necessary because the stack method always 0s before stack
		//stack_pattern.zero();
		vector<SlownessVector>::iterator uptr;
		for(uptr=slow.begin(),i=0;uptr!=slow.end();++uptr,++i)
		{
			current_stack.push_back(PWStackMember(stack_pattern,
				*uptr,noise_levels[i]));
                        //DEBUG
                        //cout << "slowness "<<i<<" t0="<<current_stack[i].t0<<endl;
		}
		/* Compute lags in a parallel vector for efficiency */
		start_offset.clear();
		start_offset.reserve(nd);
		for(dptr=d.begin();dptr!=d.end();++dptr)
		{
			int thislag;
			thislag=current_stack[0].sample_number(dptr->t0);
			start_offset.push_back(thislag);
		}

		/* NOTE DURING DEVELOPMENT:  may want to compute coherence
		estimaes.  How to most efficiently do this is less than clear. 
		Probably need a way to access residuals after the iterative
		method finishes and compare it to original data.  This may 
		require caching either raw data with nonzero weights or initial 		stacks.  */
		
		auto_ptr<ThreeComponentEnsemble> result
			= auto_ptr<ThreeComponentEnsemble> (new 
				ThreeComponentEnsemble(nu,ns_stack));
		/* Prep this gather */
		vector<ThreeComponentSeismogram>::iterator mptr;
		for(mptr=result->member.begin(),i=0;
		    mptr!=result->member.end();++mptr,++i)
		{
			/* These need to be set. Defaults are ok for other
			attributes of this object*/
			mptr->ns=ns_stack;
			mptr->t0=t0_stack;
			mptr->dt=dt_stack;
			mptr->tref=relative;
			mptr->u.zero();
			/* These are needed in output for pwmig.  Posting these
			in this way is needed for interface to save results
			with the PwmigFileHandle object. */
			mptr->put("ux0",u0.ux);
			mptr->put("uy0",u0.uy);
			mptr->put("ux",slow[i].ux);
			mptr->put("uy",slow[i].uy);
		}
		/* Main processing loop */
		int itercount(0);
                PWIndexPosition ip;
                //DEBUG
                cout << "Input data is being plotted"<<endl;
                plot_current_data();
		do {
			this->stack();
			ip=this->maxamp();
                        PWStackMember tmp(current_stack[ip.iu]);
                        /* Peak amplitude is an unstable measurement.  Use average
                           computed from dot product with individual wavelets.  
                           This sum works because weights are normalized */
                        //DEBUG
                        cout << "Polarization vectors"<<endl;
                        dvector l2amp(3);
			for(i=0;i<nd;++i)
			{
                            dvector l2this;
                            l2this=d[i].polarization(ip.iu,ip.time0);
                            l2amp+=l2this;
                            cout <<d[i].weight<<" "<< l2this(0)<<" "
                            << l2this(1)<<" "
                            << l2this(2)<<endl
                            << "Corrected"  
                            << l2this(0)/d[i].weight<<" "
                            << l2this(1)/d[i].weight<<" "
                            << l2this(2)/d[i].weight<<endl;
			}

                        /* This is needed only to make plotting work.  Remove it when
                           plotting is removed */
                        tmp.put("samprate",1.0/tmp.dt);
                        tmp.put("nsamp",tmp.ns);
                        tmp.put("time",tmp.t0);
                        plot_handle.plot(dynamic_cast<ThreeComponentSeismogram&>(tmp));
                            cout << "DEBUG  itercount="<<itercount<<endl
                                << "lag="<<ip.lag0
                                <<" or time="<<ip.time0
                                <<" iu="<<ip.iu
                                << " which is ux,uy="
                                << slow[ip.iu].ux <<" "<<slow[ip.iu].uy
                                << " peak vector values="
                                << ip.v[0]<<", "
                                << ip.v[1]<<", "
                                << ip.v[2]
                                << " L2 vector values="
                                << l2amp(0)<<", "
                                << l2amp(1)<<", "
                                << l2amp(2)<<endl
                                <<"SNR estimate for this point is "<<ip.SNR<<endl
                                <<"Refreshing plot to remove this component"
                                <<endl;
			for(int k=0;k<3;++k) ip.v[k]=l2amp(k);
			for(int k=0;k<3;++k) 
                        {
			    result->member[ip.iu].u(k,ip.lag0)+=ip.v[k];
                        }
			for(i=0;i<nd;++i)
			{
				d[i].remove(ip.iu,ip.time0,ip.v);
			}
                        plot_current_data();
			++itercount;
		} while(!converged(itercount,ip));
	} catch(...) {throw;};
}
void PPWDeconProcessor::plot_current_data()
{
    ThreeComponentEnsemble dtoplot;
    for(int i=0;i<d.size();++i)
    {
        ThreeComponentSeismogram d3c(dynamic_cast<ThreeComponentSeismogram&>
                (d[i]));
        double wt=d[i].weight;
        wt=1.0/wt;
        int n3c=3*(d3c.ns);
        dscal(n3c,wt,d3c.u.get_address(0,0),1);
        dtoplot.member.push_back(d3c);
    }
    plot_handle.plot(dtoplot);
}
/* Private method to compute noise estimates.  Algorithm uses assumes noise
   on each component is a zero mean process.  Returned estimate uses median
   absolute three-component amplitude. */
vector<double> PPWDeconProcessor::compute_noise_estimates
        (double pslat, double pslon)
{
    const string base_error("PPWDeconProcessor::compute_noise_estimates:  ");
    try {
        double dnorth,deast;  // used repeatedly so declare here
        vector<ThreeComponentSeismogram> dnoise_used;
        /* We know the size this should be from the working data ensemble
           stored as d which is the private variable nd*/
        dnoise_used.reserve(nd);
        vector<double> wt;
        wt.reserve(nd);
        /* This loop repeats calculations in compute method but the overhead
           compared to storing the weights as baggage is trivial so I intentionally
           recompute the weights here */
        double sumwt;
        int nnd=noise_data.member.size();
        int i,j;
        for(i=0,sumwt=0.0;i<nnd;++i)
        {
            double thiswt;
            double lat,lon;
            double delta,azimuth;
            const double Rearth(6378.164);
            lat = noise_data.member[i].get_double("sta_lat");
            lon = noise_data.member[i].get_double("sta_lon");
            lat=rad(lat);
            lon=rad(lon);
            thiswt=pseudostation_weight(lat,lon,pslat,pslon,
                    aperture,cutoff);
            if(thiswt>0.0)
            {
                /* We need these below and a fast and convenient
                   way to store them is in the Metadata area.  
                   Beware this will load random values for these
                   attributes in the master ensemble.  We depend on
                   them being copied.  Minor inefficiency. Further 
                   since this is a duplicate of the constructor or
                   PWDataMembers this probably should be a function
                   but I'm going to inline it */
                dist(pslat,pslon,lat,lon,&delta,&azimuth);
                delta *= Rearth;
                deast=delta*sin(azimuth);
                dnorth=delta*cos(azimuth);
                noise_data.member[i].put("dnorth",dnorth);
                noise_data.member[i].put("deast",deast);
                dnoise_used.push_back(noise_data.member[i]);
                wt.push_back(thiswt);
                sumwt += thiswt;
            }
        }
        /* This could get verbose, but this requires a warning*/
        if(dnoise_used.size() != nd)
        {
            cerr<< "PPWDeconProcessor::compute_noise_estimates (WARNING):  "
                << "noise and data size mismatch."<<endl
                <<"Pseudostation at latitude="<<deg(pslat)
                <<" and longitude="<<deg(pslon)<<endl
                <<"data ensemble has size="<<nd
                <<" but noise ensemble has size="<<dnoise_used.size()<<endl;
        }
        if(sumwt<=0.0)
        {
            stringstream errss;
            errss<< base_error
                << "sum of pseudostations weight is 0.0"<<endl;
            throw SeisppError(errss.str());
        }
        /* normalize weights.  Although it looks different this is the 
           same as what is used in the compute method */
        for(i=0;i<wt.size();++i) wt[i]/=sumwt;
        /* now scale data so weighted sum are normalized */
        for(i=0;i<wt.size();++i)
        {
            dscal(3*dnoise_used[i].ns,wt[i],dnoise_used[i].u.get_address(0,0),1);
        }
        /* The noise has to be handled different from the data because
           we require the data are assumed to be in an arrival time 
           reference frame.  For the noise we have to apply the additional
           full moveout that includes u0.  We also assume here that
           the noise was prefltered by the deconvolution operator and
           any subsequent post-filtering or the estimates made by
           this method will be very misleading to just plain wrong */
        vector<double> t0values;  // hold t0 values 
        vector<double> tendvalues;  // hold end time values
        t0values.reserve(dnoise_used.size()); 
        tendvalues.reserve(dnoise_used.size()); 
        /* Intentionally use subscripting for clarity.  Assume this is
           not a performance bottleneck here */
        double moveout,uxt,uyt;;
        dmatrix noiselags(dnoise_used.size(),slow.size());
        vector<double> tlags;
        tlags.reserve(dnoise_used.size()*slow.size());
        for(i=0;i<dnoise_used.size();++i)
        {
            dnorth=dnoise_used[i].get_double("dnorth");
            deast=dnoise_used[i].get_double("deast");
            t0values.push_back(dnoise_used[i].t0);
            tendvalues.push_back(dnoise_used[i].endtime());
            for(j=0;j<slow.size();++j)
            {
                uxt=slow[j].ux+u0.ux;
                uyt=slow[j].uy+u0.uy;
                moveout=uxt*deast+uyt*dnorth;
                noiselags(i,j)=moveout;
                tlags.push_back(moveout);
            }
        }
        vector<double>::iterator maxmoveout,minmoveout;
        vector<double>::iterator maxt0,mint0;
        vector<double>::iterator maxtend,mintend;
        int lagmax,lagmin;
        maxmoveout=max_element(tlags.begin(),tlags.end());
        minmoveout=min_element(tlags.begin(),tlags.end());
        lagmax=nint((*maxmoveout)/dt_stack);
        lagmin=nint((*minmoveout)/dt_stack);
        maxt0=max_element(t0values.begin(),t0values.end());
        mint0=min_element(t0values.begin(),t0values.end());
        maxtend=max_element(tendvalues.begin(),tendvalues.end());
        mintend=min_element(tendvalues.begin(),tendvalues.end());
        /* We now define the noise analysis window by a simple
           method that assumes the data are close to coming from a
           fixed time window.  i.e. if the data have a lot of gaps 
           this is will cause trouble.  
           */
        TimeWindow process_window((*maxt0)-(*minmoveout),
                (*mintend)-(*maxmoveout));
        if(process_window.length()<MinimumNoiseWindow)
        {
            stringstream errss;
            errss<< base_error
                << "Noise Window length is too short"<<endl
                << "Computed moveout range="<<*minmoveout<<" to "
                <<*maxmoveout<<endl
                <<"With input data this yields a time window of length "
                << process_window.length()<<endl
                <<"This is less than input cutoff parameter of "
                << MinimumNoiseWindow<<" seconds"<<endl;
            throw SeisppError(errss.str());
        }

        /* Check for gaps and if they exist just zero them for now */
        bool noise_data_have_gaps(false);
        for(i=0;i<dnoise_used.size();++i)
        {
            if(dnoise_used[i].is_gap(process_window))
            {
                noise_data_have_gaps=true;
                dnoise_used[i].zero_gaps();
            }
        }
        if(noise_data_have_gaps)
        {
            cerr<< "PPWDeconProcessor::compute_noise_estimates (WARNING):  "
                << "noise data have marked gaps."<<endl
                <<"Pseudostation at latitude="<<deg(pslat)
                <<" and longitude="<<deg(pslon)<<endl
                <<"Result may contain convergence artifacts."<<endl;
        }
        int ntostack=static_cast<int>(process_window.length()/dt_stack);
        int n3c=3*ntostack;
        vector<double> result;
        int nu=slow.size();
        result.reserve(nu);
        double *noise_stack=new double[n3c];
        vector<double> work; // holds vector amplitudes 
        work.reserve(ntostack);
        /* This factor times the mad is sigma for a normal distribution */
        const double mad_normal_multiplier(1.4826);
        for(j=0;j<nu;++j)
        {
            for(i=0;i<n3c;++i) noise_stack[i]=0.0;
            for(i=0;i<dnoise_used.size();++i)
            {
                double t=process_window.start+noiselags(i,j);
                int istart=dnoise_used[i].sample_number(t);
                if(istart<0) 
                {
                    delete [] noise_stack;
                    throw SeisppError(base_error
                        + "Coding error - computed sample number of noise segment lt 0");
                }
                double *ptr=dnoise_used[i].u.get_address(0,istart);
                for(int k=0;k<n3c;++k,ptr++) noise_stack[k]+=(*ptr);
            }
            /* compute vector amplitudes and store these in the work vector */
            work.clear();
            double atmp;
            for(i=0;i<n3c;i+=3)
            {
                atmp=0.0;
                for(int k=0;k<3;++k)atmp+=noise_stack[i+k]*noise_stack[i+k];
                atmp=sqrt(atmp);
                work.push_back(atmp);
            }
            VectorStatistics<double> stackstats(work);
            result.push_back(mad_normal_multiplier*stackstats.median());
        }
        /* Apply ceiling a floor based on upper and lower quartile.  May eventually
           want to generalize this to allow using a median and/or apply a multiplier
           to the quartiles. */
        VectorStatistics<double> rmstats(result);
        double rmslow=rmstats.q1_4();
        double rmshigh=rmstats.q3_4();
//DEBUG
cout << "Applying upper and lower bounds="<<rmshigh<<" " <<rmslow<<endl
    << "Total range of all components="<<rmstats.range()<<endl;;
        vector<double>::iterator rmsptr;
        for(rmsptr=result.begin();rmsptr!=result.end();++rmsptr)
        {
            if((*rmsptr)>rmshigh) *rmsptr=rmshigh;
            if((*rmsptr)<rmslow) *rmsptr=rmslow;
        }

        delete [] noise_stack;
        return(result);
    }
    catch (...){throw;};
}

/* Helper functions */
double pseudostation_weight(double lat,double lon,
	double pslat, double pslon, double aperture, double cutoff)
{
	const double Rearth(6378.164);
	double dnorth,deast,delta,azimuth;
	dist(pslat,pslon,lat,lon, &delta, &azimuth);
	delta *= Rearth;
	if(delta>cutoff) 
		return(-1.0);
	else
		return(exp(-(delta*delta)/(2.0*aperture*aperture)));
}
}  // end SEISPP namespace encapsulation
