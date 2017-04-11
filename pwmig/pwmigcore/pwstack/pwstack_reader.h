#include <stdio.h>
#include <vector>
#include <string>
#include "SlownessVectorMatrix.h"
using namespace std;
using namespace SEISPP;
using namespace PWMIG;
/*! \brief Abstract base class of pwstack data reader.

This is a virtual base class that allows pwstack to use polymorphism to 
sort out the correct reader.  This provides a handy way to allow multiple
input data formats.  Useful example of that type of abstraction of nothing else.
*/
class pwstack_reader
{
public:
    virtual ThreeComponentEnsemble *read_gather(int id)=0;
};
class PwstackBinaryFileReader : public pwstack_reader
{
public:
    PwstackBinaryFileReader(string fname);
    ~PwstackBinaryFileReader();
    ThreeComponentEnsemble *read_gather(int id);
    SlownessVectorMatrix gather_svm()
    {
        return svm;
    };
    int number_events(){return(ids.size());};
private:
    FILE *fp;
    vector<long> ids;
    vector<long> foffs;
    SlownessVectorMatrix svm;
};
/* Note that the gather header is a bit messy.   This can be
   read as a binary fread with sizeof(PwstackGatherHeader).  
   This is followed by a block of data that defines the
   SlownessVectorMatrix data that will drive the processing */
typedef struct PwstackGatherHeader
{
    int number_members;
    /* This holds event id inherited from database.*/
    long evid;
    /* This holds the sequence number of this gather in the
    file starting with 0 and counting sequentially */
    int sequence_number;
    /* These are event location  (radians) */
    double lon;
    double lat;
    double depth;
    double origin_time;
    /* These are the dimension of the SlownessVectorMatrix stored
       as part 2 of the gather header */
    int svmrows, svmcolumns;
} PwStackGatherHeader;
typedef struct PwstackTraceHeader
{
    char sta[8];
    double time;
    double endtime;
    int nsamp;
    double samprate;
    double atime;
    /* These hold receiver location (lat lon in radians internally)*/
    double lon;
    double lat;
    double elev;
} PwstackTraceHeader;
