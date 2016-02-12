#include <vector>
#include <string>
#include "gclgrid.h"
/* This is duplicated from StaVariableSynthethic.h.   Duplicated to avoid
   the baggage of handling yet another include file */
class LayeredModel
{
public:
    string name;
    vector<float> alpha;
    vector<float> beta;
    vector<float> rho;
    vector<float> dz;
};


/* This object is the core of this program.   We abstract it as a 
   GCLvectorfield3D object assembled from isopach layers.

   The vector is used to hold the physical properties and methods are 
   used to make it easier to extract the quantity of interest. */
class GriddedBasinModel : public GCLvectorfield3d
{
  public:
    /*! \brief Construct constant layer property model.

       For initial work this is the primary constructor.  It builds
       a basin model with variable thickness of units but constant
       physical properties.   It is riven by parallel list of file names
       that define the GCLscalarfield isopachs and list of P, S, and
       density values to be assigned to each isopach.   
       The isopach files are assumed to be output of isopach2gcl with the
       layer thickness at each node stored as the field variable. These
       are used to define dz values.   The list is also ASSUMED to be given
       from the top downward with the first surface top at sea level.   

       \param isopachs is a vector of GClscalarfield objects to be 
          assembled to build the full model.  This version is not
          robust and blindly assumes all elements are in the same
          coordinate system and of identical size.  Sanity checks 
          will cause an abort if a size mismatch is found.   
       \param vp contains a parallel (to isofiles) list of P wave velocities.
       \param vs contains a parallel (to isofiles) list of S wave velocities.
       \param rho contains a parallel (to isofiles) list of densities.
       \param nm is a name that should be assigned to the grid geometry created
     
    \exception SeisppError objects can be thrown for a variety of problems.
    This constructor is not bullet proof but I've tried to make it reasonably
    robust..
     */
    GriddedBasinModel(vector<GCLscalarfield> s,
          vector<double> vp, vector<double> vs, vector<double>rho,string nm);
    /*! Construct from a file. */
    GriddedBasinModel(string fname);
    /* If we need that sophistication we could eventually put a 
       constructor here to build a fully 3D basin model */
    /* Since this is all about velocity models this is the main 
       method.  Returns a constant velocity layered model at the
       requested lat, lon value (specify in radians please) */
    LayeredModel fetch_model(double lat0, double lon0);
    void save(string fname, string dir, 
            string format=default_output_format);
  private:
        string model_name;
};
