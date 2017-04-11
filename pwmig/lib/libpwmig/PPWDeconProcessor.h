#ifndef _PPWSDECONPROCESSOR_H_
#define _PPWSDECONPROCESSOR_H_
#include <vector>
#include <list>
#include "slowness.h"
#include "TimeSeries.h"
#include "seispp.h"
/* This is needed only for debug.  Needs to eventually be removed */
#include "SeismicPlot.h"

namespace SEISPP {
using namespace std;
using namespace SEISPP;
/*! \brief Data structure for PPWDeconProcessor.

This is a data object used for ensemble members in the plane wave 
deconvolution processing object.  It is an extension of a ThreeComponentSeismogram.
It is used to avoid excessive metadata fetches to improve efficiency.
*/
class PWDataMember: public ThreeComponentSeismogram
{
public:
	/*! Weight used to scale original data to produce this object. */
	double weight;
	/*! Location of receiver (extracted from inherited Metadata ) */
	double lat,lon,elev;
	/*! Incident wavefield slowness vector */
	SlownessVector u0;  
	TimeSeries wavelet;
	/* These hold parallel vector of lags for this this plane wave
	stack for each slowness vector */
	vector<int> lag;
	PWDataMember(ThreeComponentSeismogram& draw,
	    TimeSeries& w,
	        SlownessVector& u0in,
			vector<SlownessVector>& rayp,
				double weight,
					double lat0,
						double lon0);
	/*! Standard copy constructor */
	PWDataMember(const PWDataMember& parent);
	/*! Standard assignment operator */
	PWDataMember& operator=(const PWDataMember& parent);
        /*! Return polarization vector at a give lag.

          This is a companion to remove and the arguments are similar.
          This method results a 3 vector of particle motion at a 
          specified lag using a dot product of the wavelet with each 
          three component seismograms. Lag is computed by this formula:
            rlag0=this->sample_number(t0)-lag[iu]-wavelet.sample_number(0.0);
          \param iu integer of slowness component to use to enter lag array.
          \param t0 is time reference to which the lag is applied.
          */
        dvector polarization(int iu,double t0);
	/*! \brief Remove an estimated plane wave component 

          Subtracts a*wavelet shifted by lag computed from time t0
          and plane wave shift lag[iu] from each data member.
          This is the key recursion of this method.  Be warned
          this method does not do any error checking for speed.
          */
	int remove(int iu, double t0,double *a);
};
/*! Stacked trace used in iterative, plane wave method.

For efficiency the iterative plane wave deconvolution processing
object uses a container of these data objects to hold the current
suite of plane wave stacks.  This is a useful object as it 
contains embedded methods that simplify the plane wave algorithm. */
class PWStackMember : public ThreeComponentSeismogram
{
public:
	/*! Relative slowness vector of this stack. Note this is relative
         to a base slowness stored in processors as u0*/
	SlownessVector slow;
        /*! Contains estimate of noise rms level.

          The iteration sequence is terminated by checking against
          noise level of the stacks.  This is stored as a public
          attribute because it has to be altered during procesing. */
        double noise_rms; 
	/*! Primary constructor. 

	\param din parent 3c seismogram
	\param sv slowness vector of this stack member. */
	PWStackMember(ThreeComponentSeismogram& din,SlownessVector& sv,
                double noise_estimate);
	/*! Standard copy constructor. */
	PWStackMember(const PWStackMember& parent);
	/*! Standard assignment operator. */
	PWStackMember& operator=(const PWStackMember& parent);

	/*! Return the maximum amplitude of this member.

	This is a core method for the iterative method.  Returns
	the maximum amplitude of this 3d (stack) seismogram.  This ALWAYS 
	initiates a recalculation of the amplitudes from the current
	data.  Calls to this method must ALWAY precede calls to 
	the lag_at_max method or the return of that method will 
	be wrong.  The reason is that this method initates the 
	recalculation of the amplitude while lag_at_max simply 
	returns a value. */
	double maxamp();
	/*! Returns the lag in samples of the maximum amplitude in
	the 3c data attached to this object. */
	int lag_at_max(){return lagmax;};
        double time_at_max(){return (this->time(lagmax));};
        double SNR(){return(ampmax/noise_rms);};
private:
	double ampmax;
	int lagmax;
	/* This holds an array of 3-c amplitude values.  Update on 
	each call to maxamp() method. */
	vector<double> amp;
};
class PWIndexPosition
{
public:
	PWIndexPosition();
	PWIndexPosition(int iuin, int lag0in, double t0in, 
                double ampin, double SNRin, double *vin);
	PWIndexPosition(const PWIndexPosition& parent);
	PWIndexPosition& operator=(const PWIndexPosition& parent);
        /*! Index to slowness grid position */
	int iu;
        /*! Lag in samples for this plane wave component at origin */
	int lag0;
        /*! Relative time corresponding to lag0 */
        double time0;
        /*! L2 of amplitude at this point. */
	double amp;
        /*! Signal to noise ratio estimate at this point */
        double SNR;
        /*! Vector values at this point. */
	double v[3];
};
class PPWDeconProcessor
{
public:
	/*! Slowness vector of incident wavefield */
	SlownessVector u0;
	/*! Construct from a single gather.

	Use this constructor for a single event gather from a single
	source.  In that case we can assume the actual output wavelet
	is the same for every data member. 

	\param threecd is the input event gather.
	\param threecn is a similar ensemble to threed but containing 
            a window of noise data.  The noise data are not
            expected to be in an arrival time reference frame so 
            a larger time window is normally needed to handle
            moveout.   Further is is assumed the noise sample
            given has been passed through the same spiking filter and
            any other filters applied to the data ensemble 
	\param u0in is the incident wave slowness vector.
	\param w is the actual output used in the iterative method.
	\param md is a Metadata object used as a back door way to 
		pass an indeterminate set of control parameters.  This
		is not especially elegant, but useful for an experimental
		algorithm to avoid constantly changing the argument list. 
	\exception will throw a SeisppError object if there are problems.
	*/
	PPWDeconProcessor(ThreeComponentEnsemble& threecd,
                ThreeComponentEnsemble& threecn,
		    SlownessVector& u0in,
			vector<SlownessVector>& ulist, 
				TimeSeries& w,
					Metadata& md);
	/*! Constructor for common region gathers.

	This constructor is designed to allow this algorithm to work
	on ensembles constructed by stacking data from a common source area.
	(telecluser or equivalent).  In that case we take
	a vector of ensembles, stack them without moveouts, stack the
	actual outputs, and result looks just like a single event gather.
	The BIG complication of this is that it makes the actual output dependent
	on each station.  As a result the PWmember function has an actual
	output member.  (see above).  For single events these are all 
	the same. Use this constructor for composites and the single
	event constructor otherwise.
	\param threecd is the input event gather.
	\param threecn is a similar ensemble to threed but containing 
            a window of noise data.  The noise data are not
            expected to be in an arrival time reference frame so 
            a larger time window is normally needed to handle
            moveout.   Further is is assumed the noise sample
            given has been passed through the same spiking filter and
            any other filters applied to the data ensemble 
	\param u0in is the incident wave slowness vector.
	\param w a vector of actual output wavelets for each ensemble
		member.  This MUST be a parallel vector to the member
		vector container of threecd or chaos will result.
	\param md is a Metadata object used as a back door way to 
		pass an indeterminate set of control parameters.  This
		is not especially elegant, but useful for an experimental
		algorithm to avoid constantly changing the argument list. 
	\exception will throw a SeisppError object if there are problems.
	*/
	PPWDeconProcessor(ThreeComponentEnsemble& threecd,
            ThreeComponentEnsemble& threecn,
		SlownessVector& u0in,
			vector<SlownessVector>& ulist,
				vector<TimeSeries>& w,
					Metadata& md);
	/*! Standard copy constructor. */
	PPWDeconProcessor(const PPWDeconProcessor& parent);
	/*! Standard assignment operator.  */
	PPWDeconProcessor& operator=(const PPWDeconProcessor& parent);
				
	/*! Compute method.

	This is the primary method for this object.  It runs the
	plane wave, iterative, deconvolution algorithm and returns
	the result as a  ThreeComponentEnsemble of deconvolved data.  
	This method does all the work and can run for a lon time. 
	Each output member has key parameters stored in it's Metadata 
	area including:  gridd, ix1, ix2, ux0, uy0, ux, uy, and elev.  

	\return auto_ptr to the deconvolved data.  
	\param pslat pseudostation latitude of point to construct this
		deconvolution (radians).
	\param pslon pseudostation longitude of point to construct 
		deconvolution (radians).  
	\exception can throw a SeisppError if there are unreconciable 
		errors.
	*/
	auto_ptr<ThreeComponentEnsemble> compute(double pslat, double pslon);
private:
        /* DEBUG only */
        SeismicPlot plot_handle;
        void plot_current_data();
        /// end DEBUG additons
	int nd; /* d.size() cached for efficiency */
	int nu;  /* slow.size() cached for efficiency */
	int ns_stack;  /* stack ns */
	double t0_stack;  /* Relative time t0 for plane wave stacks */
	double dt_stack;  /* sample interval of stack data */
	double aperture, cutoff;  /* pseudostation sigma and cutoff */
	/* Cache these for efficiency */
	vector<double> stalat,stalon;
	/* This is an nu length vector parallel to slow containing offset
	value for an internal buffer used to accumulate plane wave 
	components.  */
	vector<int> start_offset;
	/* This is used to determine stack buffer size.  It is set
	by constructor as a conservative estimate to allow any 
	feasible plane wave shifts (scaled by maximum aperture and
	range of slowness vectors.  Note because the stack is 3c the
	total stack size is 3*stacksize, but held in a dmatrix object
	created and then destroyed  when the stack method is called*/
	int stacksize;
	/* This object uses a vector of slowness values to be more general
	than the rectangular slowness grid used in pwmig.  This is
	intentional to make this more general. */
	vector<SlownessVector> slow;
	vector<PWDataMember> d;
	vector<TimeSeries> wavelet;
	ThreeComponentEnsemble parent_data;
        ThreeComponentEnsemble noise_data;
	/* Plane wave stack components (length nu) */
	vector<PWStackMember> current_stack;
        /* cutoff parameter.  Will throw an exception if noise time window is 
           less than this length when allowing for moveout */
        double MinimumNoiseWindow; 
        /* if true use SNR instead of absolute amplitude to find max in
           each iteration */
        bool UseSNR;
        /* This is a critical convergence variable.  Stop when computed max snr
           falls below this value */
        double SNRfloor;
        /* Limit on number of iterative cycles */
	int maxiteration;  
        /* Used to compute a running average of snr values over past 
           queuesize iterations*/
        int snravgsize;
        list<double> recent_snr;
        /*previous average snr computed from recent_snr*/
        double last_snravg;
	/*! workhorse internal method to reform stack in place.  i.e
	no copying is done, just components of current_stack are updated*/
	void stack();
	PWIndexPosition maxamp();
	bool converged(int itercount,PWIndexPosition& ip);
        vector<double>  compute_noise_estimates(double pslat, double pslon);
};
double pseudostation_weight(double lat, double lon,
	double pslat, double pslon, double aperture, double cutoff);
}  /* End SEISPP Namespace encapsulation */
#endif
