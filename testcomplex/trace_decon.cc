#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <omp.h>
#include "perf.h"
#include "DeconOperator.h"
#include "SeisppKeywords.h"
#include "TimeSeries.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
#include "Metadata.h"
#include "Hypocenter.h"
#include "dbpp.h"
#include "filter++.h"
#include "resample.h"
#include "seispp.h"
#include "stack.h"

#include "SignalToNoise.h"
void usage()
{
	cerr << "Version: 2015.12.01" <<endl;
	cerr << "trace_decon db [-d output_dir -o outputdb -n start_evid -pf pfname]" << endl;
	exit(-1);
}

int main(int argc, char **argv)
{
	cout << "This is a test function" << endl;

}
