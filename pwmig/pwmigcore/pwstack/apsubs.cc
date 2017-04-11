#include <math.h>
#include "coords.h"
#include "stock.h"
#include <iostream>
using namespace std;

#define RADIUS_EARTH 6378.164
/* For plane wave moveout computations a local cartesian coordinate
system is used wrt to a particular origin.  This approximation is
reasonable until the distance of a station from the origin
becomes a significant fraction of the epicentral distance.
In pseudostation stacking we can minimize the impact of this approximation
by continually translating the origin to the pseudostation point and
computing distances to all stations wrt to that point.  Since large
distance stations receive a low weight, time alignments at large
distances become unimportant.  this function is used to computer
vector distances in this context.

Arguments:
lat0, lon0 - origin to compute dnorth deast from
(IN RADIANS)
lat, lon - lat an lon of point to be converted to dnorth deast
dnorth, deast - output dnorth, deast vector

Author:  G Pavlis
Written:  June 2000
*/
void geographic_to_dne(double lat0, double lon0,
double lat, double lon, double *dnorth, double *deast)
{
    double azimuth;
    double d;

    dist(lat0,lon0,lat,lon,&d,&azimuth);
    d *= RADIUS_EARTH;
    *deast = d*sin(azimuth);
    *dnorth = d*cos(azimuth);

}


/* This function computes pseudostation weights for a 2d Gaussian function
to implement pseudostation stacking as described in Neal and Pavlis (1999) grl.

Arguments:
    nsta - number of stations
    dnorth, deast - parallel vectors of length nsta of deast-dnorth coordinates
        relative to a reference geographical position.
    aperture - sigma parameter in Gaussian formula
    cutoff - stations outside this distance (km) are given 0 weight
    w - vector of length nsta (output assumed allocated by caller) to
        hold weights.

The function uses a dnorth-deast approximation.  For every data set I know
of this is not a problem for computing these weights or even for computing
moveout.  The reason is that coherence is lost before this approximation
becomes a problem.

Returns count of nonzero weights.

Author:  G Pavlis
Written:  June 2000
Modified:  Converted to C++ and assimilated into new a working implementation
in March 2003
*/
int compute_pseudostation_weights(int nsta, double *dnorth, double *deast,
double aperture, double cutoff, double *w)
{
    int i;
    double distance, azimuth, denom;
    int count=0;

    denom = 2.0*aperture*aperture;

    for(i=0;i<nsta;++i)
    {
        distance = hypot(dnorth[i],deast[i]);
        if(distance>cutoff)
            w[i] = 0.0;
        else
        {
            w[i] = exp(-distance*distance/denom);
            ++count;
        }
    }
    return(count);
}


/* Computes vector of nsta moveout corrections for slowness vector
 * (ux,uy) using local coordinates stored in parallel deast, dnorth
 * vectors (both assumed length nsta).
 * Assumes moveout has been allocated previusly.
 */
void compute_pwmoveout(int nsta,double *deast, double *dnorth,
double ux, double uy, double *moveout)
{
    int i;
    for(i=0;i<nsta;++i)
    {
        moveout[i] = ux*deast[i] + uy*dnorth[i];
    }
}
