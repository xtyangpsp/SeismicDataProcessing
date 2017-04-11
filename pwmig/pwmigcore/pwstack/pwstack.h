#include "SeisppError.h"
#include "ensemble.h"
#include "mute.h"
#include "PfStyleMetadata.h"
#include "PwmigFileHandle.h"
/*  Avoided now with a define in Makefile, but leave here 
for reference. 
#define _FILE_OFFSET_BITS 64
*/
using namespace std;
using namespace SEISPP;

class DepthDependentAperture
{
public:
	int npoints;
	double *t;
	double *aperture;
	double *cutoff;
	DepthDependentAperture()
		{npoints=0;t=NULL;aperture=NULL;cutoff=NULL;};
        /*! Construct using linear segments specified by parallel vectors*/
	DepthDependentAperture(int npt){
		npoints=npt;
		t=new double[npt];
		aperture=new double[npt];
		cutoff=new double[npt];
	};
        /*! Use Fresnel zone formula to construct.  
          *
          * A useful depth dependent aperture formula is to use the
          * Fresnel zone size for vertical incident S in a constant
          * velocity medium.  This constructor implements that formula.
          * The aperture is parameterized internally by two parallel 
          * vectors just like the manual constructors.  The dtau and ntau
          * variable set the uniform grid interval and number of points
          * respectively.
          *
          *\param vs is S wave velocity.
          *\param vp is P wave velocity.
          *\param period is wave period to compute the Fresnel zone size.
          *     Be aware this is effectively a width parameter because at
          *     zero lag the width becomes this value and increases with 
          *     lag.
          *\param dtau is lag sample interval for formula
          *\param ntau is the number of lags to compute (Note that the 
          *      interpolation method continues the value of (ntau-1)*dtau 
          *      to infinity.
          *\param cutoffmultiplier is a scale factor used to compute cutoff
          *     parameter.  For this constructor cutoff is computed as this
          *     number times the aperture.
          *\param echoresult if set true will print computed aperture sizes
          *     using cout. 
         */
        DepthDependentAperture(double vs, double vp,double period,
            double dtau, int ntau,double cutoff_multiplier,bool echoresult);
        /*! Parameter file driven constructor.   
          
          Uses a pf format file through a nonantelope interface called
          a PfStyleMetadata object.   

          \param md is the object created from a pf file.
          \param tag is an &Arr tag bounding he parameters for this object.
          */
	DepthDependentAperture(PfStyleMetadata& md,string tag);
        /*! Copy constructor. */
	DepthDependentAperture(const DepthDependentAperture& parent);
	~DepthDependentAperture(){
		if(t!=NULL)delete [] t; 
		if(aperture!=NULL)delete [] aperture;
		if(cutoff!=NULL)delete [] cutoff;
	};
	double get_aperture(double tnow);
	double get_cutoff(double tnow);
        DepthDependentAperture& operator=(const DepthDependentAperture& parent);
};

/* Special internal object for this program used to return
two different coherence estimate: (1) component by component
estimates and (2) vector sum coherence (scalar). 
*/
class Coharray
{
public:
	double t0;
	double dt;
	int ns;
	double winlen;  // width of averaging window in time units
	dmatrix compcoh;
	vector<double> coh;
	Coharray(double t0i,double dti, int nsi,double win)
	{
		t0=t0i;
		dt=dti;
		ns=nsi;
		winlen=win;
		coh.reserve(ns);
		compcoh=dmatrix(3,ns);
		compcoh.zero();
	};
	/* Standard copy constructor.*/
	Coharray(const Coharray& parent);
	Coharray& operator=(const Coharray& other);
};
void geographic_to_dne(double lat0, double lon0, double lat, double lon, 
		double *dnorth, double *deast);
void geographic_to_dne(double lat0, double lon0,
        double lat, double lon, double *dnorth, double *deast);
int compute_pseudostation_weights(int nsta, double *dnorth, double *deast,
                double aperture, double cutoff, double *w);
void compute_pwmoveout(int nsta,double *deast, double *dnorth,
                double ux, double uy, double *moveout);
int pwstack_ensemble(ThreeComponentEnsemble& indata,
	RectangularSlownessGrid& ugrid,
	TopMute& mute,
	TopMute& stackmute,
	int stack_count_cutoff,
	double tstart,
	double tend,
	DepthDependentAperture& aperture,
	double aperture_taper_length,
	double centroid_cutoff,
	double dtcoh,
	double cohtwin,
	MetadataList& mdlcopy,
	PwmigFileHandle& dfh,
	PwmigFileHandle& coh3cfh,
	PwmigFileHandle& cohfh);
string virtual_station_name(int ix1, int ix2);
Coharray compute_stack_coherence(vector<dmatrix> d,
	dmatrix& wgt,ThreeComponentSeismogram& stack,
	double dtcoh,double cohavgwin,TopMute& mute);
void save_coh(Coharray& coh, ThreeComponentSeismogram& stack,
	PwmigFileHandle& pfhcoh3c, PwmigFileHandle& pfhcoh);
//DEBUG
#ifdef MATLABDEBUG
extern MatlabProcessor mp;
#endif
