#include <sstream>
#include "stock.h"
#include "perf.h"
#include "PfStyleMetadata.h"
#include "pwstack.h"
using namespace std;
//
// Pf-based constrructor for DepthDependentAperture object.
//
//  Arguments:
//      md - PfStyleMetadata object created from a pf file.
//	tag - name identifying the Tbl defining the designed input
//
//  This routine throws an exception with a message string if the required
// information from the pf cannot be retrieved.
// This should probably normally be fatal.  A more complex return
// to the error was deemed unnesseary since this is a specialized
// object for pwstack
//
//  Author:  Gary L. Pavlis
//

DepthDependentAperture::DepthDependentAperture(PfStyleMetadata& md,
        string tag) 
{
    try {
        list<string> tlist=md.get_tbl(tag);
        npoints=tlist.size();
        t=new double[npoints];
        aperture=new double[npoints];
        cutoff=new double[npoints];
        list<string>::iterator tptr;
        int i;
        for(tptr=tlist.begin(),i=0;tptr!=tlist.end();++tptr,++i)
        {
            stringstream instr(tptr->c_str());
            instr >> t[i];
            instr >> aperture[i];
            instr >> cutoff[i];
        }
    }catch(...){throw;};
}
/* Fresnel zone constructor.  See include file for parameter definitions. */
DepthDependentAperture::DepthDependentAperture(double vs, 
        double vp,
            double period,
                double dtau, 
                    int ntau,
                        double cutoff_multiplier,
                            bool echoresult)
{
    npoints = ntau;
    t=new double[npoints];
    aperture=new double[npoints];
    cutoff=new double[npoints];
    int i;
    for(i=0;i<npoints;++i) t[i]=static_cast<double>(i)*dtau;
    double lagtots=vp/(vp-vs);  // multiplier from S-P time to S time
    double ts;
    if(echoresult) cout << "Fresnel zone aperture calculation used"<<endl
                        << "lag   aperture(km)    cutoff(km)"<<endl;
    for(i=0;i<npoints;++i)
    {
        ts=lagtots*t[i];
        double term1=ts+period/2.0;
        aperture[i]=vs*sqrt(term1*term1 - ts*ts);
        cutoff[i]=cutoff_multiplier*aperture[i];
        if(echoresult)cout<<t[i]<<"  "<<aperture[i]<<"  "<<cutoff[i]<<endl;
    }
}

DepthDependentAperture::DepthDependentAperture(const DepthDependentAperture& a)
{
    npoints=a.npoints;
    if(npoints>0)
    {
        if(a.t == NULL)
        {
            t=NULL;
        }
        else
        {
            t=new double[npoints];
            dcopy(npoints,a.t,1,t,1);
        }
        if(a.aperture == NULL)
        {
            aperture=NULL;
        }
        else
        {
            aperture=new double[npoints];
            dcopy(npoints,a.aperture,1,aperture,1);
        }
        if(a.cutoff == NULL)
        {
            cutoff=NULL;
        }
        else
        {
            cutoff=new double[npoints];
            dcopy(npoints,a.cutoff,1,cutoff,1);
        }
    }
    else
    {
        t=NULL;
        aperture=NULL;
        cutoff=NULL;
    }
}


double linear_interpolate_irregular(int npts, double *x, double *y,
double xsearch)
{
    int i1,i2;
    double yout;
    //silently return endpoints if outside the range
    if (xsearch<x[0]) return(y[0]);
    if (xsearch>x[npts-1]) return(y[npts-1]);
    for(i2=1;i2<npts;++i2)
        if(x[i2]>xsearch) break;
    i1 = i2-1;
    yout = y[i1] + (xsearch-x[i1])*(y[i2]-y[i1])/(x[i2]-x[i1]);
    return(yout);
}


double DepthDependentAperture::get_aperture(double t0)
{
    return(linear_interpolate_irregular(npoints,t,aperture,t0));
}


double DepthDependentAperture::get_cutoff(double t0)
{
    return(linear_interpolate_irregular(npoints,t,cutoff,t0));
}
DepthDependentAperture& DepthDependentAperture::operator=(const DepthDependentAperture& parent)
{
    if(this!=&parent)
    {
        npoints=parent.npoints;
        if(npoints>0)
        {
            if(parent.t == NULL)
            {
                t=NULL;
            }
            else
            {
                t=new double[npoints];
                dcopy(npoints,parent.t,1,t,1);
            }
            if(parent.aperture == NULL)
            {
                aperture=NULL;
            }
            else
            {
                aperture=new double[npoints];
                dcopy(npoints,parent.aperture,1,aperture,1);
            }
            if(parent.cutoff == NULL)
            {
                cutoff=NULL;
            }
            else
            {
                cutoff=new double[npoints];
                dcopy(npoints,parent.cutoff,1,cutoff,1);
            }
        }
        else
        {
            t=NULL;
            aperture=NULL;
            cutoff=NULL;
        }
    }
    return(*this);
}
