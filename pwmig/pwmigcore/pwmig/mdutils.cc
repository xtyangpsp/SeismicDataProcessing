#include <vector>
using namespace std;
#include "stock.h"
#include "Metadata.h"
#include "gclgrid.h"
#include "slowness.h"
#include "Hypocenter.h"
#include "ensemble.h"
#include "seispp.h"
using namespace SEISPP;
#include "pwmig.h"

SlownessVector get_stack_slowness(ThreeComponentEnsemble& ensemble)
{
	SlownessVector slow;

	// We extract these from the global metadata area for the ensemle.
	// conceptually correct, but requires more careful implementation and assumes
	// the data are there.
	try {
		slow.ux = ensemble.get_double("plane_wave_ux");
		slow.uy = ensemble.get_double("plane_wave_uy");
	} catch (MetadataError& mde)
	{
		throw(mde);
	}
	return(slow);
}

/*
Hypocenter get_event_hypocenter(ThreeComponentEnsemble& ensemble)
{
	Hypocenter h;

	try {
		h.lat = ensemble.get_double("origin.lat");
		h.lon = ensemble.get_double("origin.lon");
		h.z = ensemble.get_double("origin.depth");
		h.time = ensemble.get_double("origin.time");
	} catch (MetadataError& mde)
	{
		throw mde;
	}
	return(h);
}
*/
