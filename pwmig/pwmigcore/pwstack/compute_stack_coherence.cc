#include <math.h>
#include <vector>
#include "TimeWindow.h"
#include "perf.h"
#include "dmatrix.h"
#include "ThreeComponentSeismogram.h"
#include "pwstack.h"
#include "seispp.h"
#include "SeisppError.h"
using namespace std;
using namespace SEISPP;

double coherence(double ssqr,double ssqd)
{
    // This is an ambiguous choice to avoid a NaN but
    // I make it 0 as a safer choice.  Avoids 0/0
    if(ssqr<=0.0 || ssqd <= 0.0)
        return(0.0);
    else if(ssqr>ssqd)
        return(0.0);
    else
        return(1.0-sqrt(ssqr)/sqrt(ssqd));
}


/*  Computes stack coherence given an input ensemble sent as
a set of 3 matrices (one for each component) and parallel
weight matrix.  This procedure is linked to pwstack making
assumptions that make it far from general.  Key assumptions
are:  (1)  d and wgt have same t0 as stack
(2)  top mute has times consistent with above.

Uses a moving window computation of coherence over a
specified time gate of length ncoh in stack time sampling
units.
Arguments
d - 3 components stored in matrices with constant start time.
wgt - weight array to form stack.
stack - stack used to compute coherence.
dtcoh - sample interval for output coherence array.
cohwinlen - length (in time units) of coherence averaging window.
mute - top mute that is assumed to have been applied to d.
Coherence is forced to 0 in the zero area of the mute zone.
*/
Coharray compute_stack_coherence(vector<dmatrix> d,
dmatrix& wgt,
ThreeComponentSeismogram& stack,
double dtcoh,
double cohwinlen,
TopMute& mute)
{
    int i,j,k;

    // IMPORTANT assume size of d components and stack are always the same
    // Caller must be sure this is true. This was done for
    // efficiency assuming this code will not be used outside
    // of the pwstack program
    int ns=stack.ns;
    int nd=d[0].columns();
    // This is the number of coherence points for output
    int ncoh=static_cast<int>((stack.endtime()-stack.t0)/dtcoh);
    ++ncoh;                                       // standard interval versus points issue
    Coharray result(stack.t0,dtcoh,ncoh,cohwinlen);
    // Note these arrays are in ThreeComponentSeismogram
    // convension with components in rows of the matrix
    dmatrix r(3,ns);                              // accumulated weighted squared residual array
    // holds accumulated weighted sum of squares for norm d
    dmatrix dwtssq(3,ns);
    dmatrix work(3,ns);                           // simplifies algorithm below
    r.zero();   work.zero();  dwtssq.zero();
    double dtover2=dtcoh/2.0;
    double wtd;
    // compute residual matrix
    for(k=0;k<nd;++k)
    {
        for(i=0;i<3;++i)
        {
            // d is a matrix with seismograms in columns
            double *ptr=d[i].get_address(0,k);
            dcopy(ns,ptr,1,work.get_address(i,0),3);
        }
        for(i=0;i<3;++i)
        {
            for(j=0;j<ns;++j)
            {
                wtd=work(i,j)*wgt(j,k);
                dwtssq(i,j)=dwtssq(i,j)+wtd*wtd;
            }
        }
        work=work-stack.u;
        for(i=0;i<3;++i)
        {
            for(j=0;j<ns;++j)
            {
                wtd=work(i,j)*wgt(j,k);
                r(i,j)=r(i,j)+wtd*wtd;
            }
        }
    }
    // Some useful variables best computed only once for efficiency
    double winhalfwidth=cohwinlen/2.0;
    int nwin=static_cast<int>(cohwinlen/stack.dt) + 1;
    int count;
    double rcompsum[3],dcompsum[3],rsum,dsum;
    // Now we need to fill the output arrays.  This is a simple moving
    // average algorithm.  The oddities here are: (1) the top mute zone
    // is zeroed, and (2) the endtime limit has to deal with the overlap
    // variable.  the time gate is shortened when necessary on both ends
    double t;
    int is;
    for(j=0,t=stack.t0;j<ncoh;++j,t=t+dtcoh)
    {
        TimeWindow gate(t-winhalfwidth,t+winhalfwidth);
        if(t<mute.t0e)
        {
            // Don't need to alter compcoh because creation
            // initializes it to 0s
            result.coh.push_back(0.0);
        }
        else
        {
            for(i=0;i<3;++i)
            {
                rcompsum[i]=0.0;
                dcompsum[i]=0.0;
            }
            rsum=0.0;
            dsum=0.0;
            for(count=0,is=stack.sample_number(t);
                count<nwin && is<ns;
                ++is)
            {
                if(is>=0)
                {
                    ++count;
                    for(i=0;i<3;++i)
                    {
                        rcompsum[i]+=r(i,is);
                        dcompsum[i]+=dwtssq(i,is);
                        rsum+=r(i,is);
                        dsum+=dwtssq(i,is);
                    }
                }
            }
            if(count>0)
            {
                for(i=0;i<3;++i)
                    result.compcoh(i,j)=coherence(rcompsum[i],dcompsum[i]);
                result.coh.push_back(coherence(rsum,dsum));
            }
            else
            {
                for(i=0;i<3;++i)result.compcoh(i,j)=0.0;
                result.coh.push_back(0.0);
            }
        }
    }
    return(result);
}
Coharray::Coharray(const Coharray& parent)
{
	t0=parent.t0;
	dt=parent.dt;
	ns=parent.ns;
	winlen=parent.winlen;
	compcoh=parent.compcoh;
	coh=parent.coh;
}
Coharray& Coharray::operator=(const Coharray& parent)
{
	if(this!=&parent)
	{
		t0=parent.t0;
		dt=parent.dt;
		ns=parent.ns;
		winlen=parent.winlen;
		compcoh=parent.compcoh;
		coh=parent.coh;
	}
	return(*this);
}

/*!  Save coherence data as a set of seismic attribute traces
linked by attribute link table.

Arguments -
    coh - result of coherence calculator used in this program.
    stack - stack from which coh was derived.  It is assumed
    pfhcoh3c - special file handle for pwmig to hold 3c data
    pfhcoh - special file handle for pwmig for scalar coh data
*/

void save_coh(Coharray& coh,
  ThreeComponentSeismogram& stack,
    PwmigFileHandle& pfhcoh3c,
      PwmigFileHandle& pfhcoh)
{
    try
    {
        // First clone the stack.  We'll then modify it's contents to match the
        // coherence data and then save the results using dbsave for ThreeComponentSeismogram
        // and TimeSeries data
        ThreeComponentSeismogram coh3ctrace(stack);
	/* We have to put these into the trace objects as the
	ExtractComponent procedure requires them */
	coh3ctrace.put("samprate",1.0/coh.dt);
	coh3ctrace.put("time",coh.t0);
	coh3ctrace.put("nsamp",coh.ns);
        auto_ptr<TimeSeries>cohtrace(ExtractComponent(coh3ctrace,0));
        coh3ctrace.ns=coh.ns;
        cohtrace->ns=coh.ns;
        coh3ctrace.t0=coh.t0;
        cohtrace->t0=coh.t0;
        coh3ctrace.dt=coh.dt;
        cohtrace->dt=coh.dt;
        cohtrace->s=coh.coh;
        coh3ctrace.u=coh.compcoh;
	pfhcoh3c.save(coh3ctrace);
	pfhcoh.save(*cohtrace);;
    }
    catch (...)
    {
        throw;
    }
}
