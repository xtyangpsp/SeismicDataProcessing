#include <vector>
#include "Hypocenter.h"
#include "SyntheticSeismogram.h"
#include "ThreeComponentSeismogram.h"
using namespace std;
using namespace SEISPP;
/* This defines a format name for a single velocity.  Currently this 
      is the only name allowed */
const string default_model_format("glpvmodels");
class LayeredModel
{
public:
    string name;
    vector<float> alpha;
    vector<float> beta;
    vector<float> rho;
    vector<float> dz;
};
/*! \brief Synthethic seismogram generator for station variable layered models.

  The purpose of this object is to provide a way to generate synthetic P to S 
response functions with constant velocity layered models but with the model being
specified independently for individual seismic stations.   The internal algorithm 
is a version of plane wave synthetics I obtained originally from Hersh Gilbert
at Purdue University, but which seems to be a descendent of fortran subroutine 
originally written by Brian Kennett or at least based on his book.   This object is
effectively a C++ interface to that propagator matrix code.  
  Models are read from files.   The files are assumed to be defined by unix path 
file names with the file name at the end of the path being the station name.   Each 
file is assumed a simple ascii file with four columns of data: (1) P-velocity,
(2) S-velocity, (3) density, and (4) thickness.   Velocities are assumed to be in km/s,
density in gm/cc, and thicknesses in km.   (I take no responsibiity for the bastard units).  
  Note it is allowed to have an empty list of file names in which case the default will 
always be used.  This is an easy way to use a common model for all data.   In most cases,
however, that will likely also require calling the set_default_model method.
  The synthetics created by this program are always filtered on output by a Gaussian
filter to a specified pulse with (in time using a gaussian sigma parameter).  Be warned this
filter is not properly normalized so the absolute amplitudes are distorted from the original
Green's function computed by the propagator matrix code.  
  The output can be in one of two forms.  The default output is simulated receiver functions.
Receiver functions are not the true impulse response because true receiver functions are 
produced by deconvolution of the radial component with the vertical station by station.  
When run this way the output has the vertical component zeroed and the horizontal components
defined by spliting the radial component into its components in cardinal directions.   
(ThreeC seismograms in SEISPP have 1=East, 2=North, and 3=Up when set cardinal).   The 
user can also compute the true impulse response function in which case the vertical will 
have a train of pulses and the horizonals will have radial date in cardinal components as
with the RF case.  
*/
class StaVariableLayeredSynthetic : public SyntheticSeismogram
{
public:
    /*! \brief construct from a list of file names.

      Assumes the list is a set of strings that define valid file names.
      The assumption is that the names a relative or absolute unix path 
      names with the leaf (right) being the file name that is the name
      of a seismic station for which a model is to be defined.   These
      names are used internally as keys to index the entire suite of models.
      Each file is assumed to contain ascii data with four columns:
      p-velocity, s-velocity, density, thickness
      \param filelist is an stl list of file names.   The name structure
        is important.  It is blindly assumed file names are unix file
        names with the leaf of the tree being a station name.
      \param tsig is time sigma.  The output seismograms are filtered by 
        a Gaussian filter (in the frequency domain) to build pulses of this
        width measure in time. 
      \param wlev is a water level used for frequency domain waterlevel 
        deconvolution to simulate the radial being a standard receiver
        function computed by deconvolution with the vertical component.
      \param MakeRF is a boolean that controls the type of output.  When
        true output is plane wave radial and vertical response to model.  
        When false the vertical component will be zero and the radial will
        be a simulated receiver function produced by deconvolving the 
        radial with the vertical (waterlevel).

       \exception throws a SeisppError exception mostly related to reading
       files or inconsistencies in files.
       */
   StaVariableLayeredSynthetic(list<string> filelist, double tsig,
           double wlev=0.01,bool MakeRF=false);
   /*! \brief Construct from a single file.

     This constructor builds the object from a single large file.   
     At present this is a single format that is defaulted, but 
     the interface allows alternatives if one wished to do the
     coding to support such an alternative.  The default format
     is best viewed as the same as concatenating the single station
     files used by the list driven constructor adding a header line
     before each station model with two arguments:  station_name and
     number of layers (lines to follow).

     An important reason this constuctor was added is that the list
     driven version was ridiculously inefficient due to the large 
     number of tiny file it required. 

     \param modfile is the name of the file containing all the models.
      \param tsig is time sigma.  The output seismograms are filtered by 
        a Gaussian filter (in the frequency domain) to build pulses of this
        width measure in time. 
      \param wlev is a water level used for frequency domain waterlevel 
        deconvolution to simulate the radial being a standard receiver
        function computed by deconvolution with the vertical component.
      \param MakeRF is a boolean that controls the type of output.  When
        true output is plane wave radial and vertical response to model.  
        When false the vertical component will be zero and the radial will
        be a simulated receiver function produced by deconvolving the 
        radial with the vertical (waterlevel).
      \param format - optional format specification.   Currently will
       throw an exception if it is not the default.

       \exception throws a SeisppError exception mostly related to reading
       files or inconsistencies in files.

       */
   StaVariableLayeredSynthetic(string modfile, double tsig,
           double wlev=0.01,bool MakeRF=false,string format=default_model_format);
   /*! Compute a scalar time series - not implemented. */
   TimeSeries ComputeScalar(int nsamp, double dt, 
           Hypocenter& hypo,double rlat,double rlon,double relev,
           string options);
   /*! Compute a scalar time series - not implemented. */
   TimeSeries ComputeScalar(const TimeSeries& parent, Hypocenter& hypo,
           double rlat, double rlon, double relev, string type);
   /* Compute a 3C seismogram using a pattern seismogram.

      This is the primary method for computing synthetics with this 
   algorithm.   It uses the parent seismogram and builds a clone 
   of parent aimed as a simulation of parent.   It assumes the data
   can be cast to a relative reference frame with the zero of the 
   seismogram set to the P wave arrival time.  Output data are rotated
   into cardinal coordinates in standard seispp 3c coordinates.   

  \param parent is the 3c trace to clone.
  \param hypo is the hypocenter object used to get source coordinates.
  \param rlat is station (receiver) latitude in radians
  \param rlon is station (receiver) longitude in radians
  \param relev is station (receiver) elevation in km
  \param options is an option string (currently ignored).

  */

   ThreeComponentSeismogram Compute3C(const ThreeComponentSeismogram& parent,
               Hypocenter& hypo, double rlat, double rlon, double relev,
               string options);
   /* Compute a 3C seismogram with minimal input. 

      This is the secondaryethod for computing synthetics with this 
   algorithm.   It builds a 3c seismogram of a specified sample interval
   and length.   sta name is passed in an ugly fashion through the options
   string.  (see options description)  It is actually implemented by 
   building a simple pattern and calling the Compute3C method that
   clones a parent.

  \param nsamp is the number of (vector) samples for output.
  \param dt is the sample interval of the output in seconds.
  \param hypo is the hypocenter object used to get source coordinates.
  \param rlat is station (receiver) latitude in radians
  \param rlon is station (receiver) longitude in radians
  \param relev is station (receiver) elevation in km
  \param options is used to define the station name and (optionally) 
    a time offset for sample 0.   String should contain at least a
    station name.  Optionally append a time (in seconds) separated by
    white space to define t0 of the output. 

  */
   ThreeComponentSeismogram Compute3C(int nsamp, double dt, Hypocenter& hypo,
               double rlat, double rlon, double relev,string options);
   /*! Set default velocity model.

     This processing object uses a default crustal model that is hard wired
     into the code.  It is essentially a one layer iaspei91 continental crust.  
     This little method allows the user to change this default. */
   void set_default_model(LayeredModel& defmod){default_model=defmod;};
   /*! Return the velocity model associated with a particular station name.
     
     This method does the job of finding a velocity model by sta.  It silently
   returns default if staname is not found.  Hence this routine never fails.
   To be explicit to get default model use staname="default" and unless you
   (foolishly) have a station named default you will get what you want.  
   \param staname  is the name of the station for which a model is requested.*/
   LayeredModel get_model(string staname);
   /*! Update a model for a station.

     I may be necessary at times to change a velocity model for one or more
     stations.  This does that.  If that station did not have a model before
     it will be added. 

     \param sta is the station name used as a key for the model passed.
     \param stamod is the model assigned to this station 
     */
   void set_model(string sta, LayeredModel& stamod);
private:
   /* Models are stored in this associative array keyed by station name */
   map<string,LayeredModel> mods;
   /* Stations without a model defined default to this one */
   LayeredModel default_model;   
   /* Set true if sythetics are to be converted to RF with a waterlevel
      deconvolution operator */
   bool ConvertToRF;
   float tsigma,wlevel;
};
