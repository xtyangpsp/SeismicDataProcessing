#include <stdio.h>
#include <vector>
#include <sstream>
#include "seispp.h"
#include "SeisppError.h"
#include "Metadata.h"
#include "pwstack_reader.h"
using namespace std;
using namespace SEISPP;
/* Helpers.

   This one takes the gather header structure and returns a Metadata 
   object that can be passed to the 3C constructor.  A bit inefficient, but
   a stable solution */
Metadata PGHtoMD(PwStackGatherHeader& gh)
{
    /* The names here are frozen and a maintenance issue is matching
       these to get routines in pwstack code. */
    Metadata result;
    result.put("evid",static_cast<long>(gh.evid));
    /* assumed lon and lat are radians and depth is km */
    result.put("lon",gh.lon);
    result.put("lat",gh.lat);
    result.put("time",gh.origin_time);
    result.put("depth",gh.depth);
    return(result);
}
/* comparable routine for trace headers */
Metadata THtoMD(PwstackTraceHeader& th)
{
    Metadata result;
    result.put("sta",th.sta);
    result.put("time",th.time);
    result.put("endtime",th.endtime);
    result.put("nsamp",static_cast<int>(th.nsamp));
    result.put("samprate",th.samprate);
    /* assumed this is relative time - maintenance issue warning*/
    result.put("atime",th.atime);
    result.put("site.lon",th.lon);
    result.put("site.lat",th.lat);
    result.put("site.elev",th.elev);
    /* Data read are assumed to have been oriented to cardinal
       directions.  We need to set these to connect with 
       constructor in seispp library. */
    result.put("components_are_cardinal",true);
    return(result);
}
PwstackBinaryFileReader::PwstackBinaryFileReader(string fname)
{
    const string base_error("PwstackBinaryFileReader constructor:  ");
    fp=fopen(fname.c_str(),"r");
    if(fp==NULL) throw SeisppError(base_error
            + "fopen failed for input file="+fname);
    /* Load the index vectors.  First read a directory seek position */
    long diroffset;
    if(fread(&diroffset,sizeof(long),1,fp)!=1) 
        throw SeisppError(base_error+"fread failed reading diroffset");
    int nevents;
    if(fread(&nevents,sizeof(int),1,fp)!=1)
        throw SeisppError(base_error+"attempt to read nevents of index failed");
    if(fseek(fp,diroffset,SEEK_SET))
        throw SeisppError(base_error+"fseek error to diroffset value read");
    long *idin,*foffin;
    idin = new long[nevents];
    foffin = new long[nevents];
    if(fread(idin,sizeof(long),nevents,fp)!=nevents)
        throw SeisppError(base_error+"attempt to read evid vector in index failed");
    if(fread(foffin,sizeof(long),nevents,fp)!=nevents)
        throw SeisppError(base_error+"attempt to read file offset vector in index failed");
    /* A potential memory lead above if either of the two throws are
       executed.  Did this intentionally because in this use if any of these
       errors are thown pwstack will abort so there is no risk of a memory 
       leak. */
    int i;
    ids.reserve(nevents);
    for(i=0;i<nevents;++i) ids.push_back(idin[i]);
    foffs.reserve(nevents);
    for(i=0;i<nevents;++i) foffs.push_back(foffin[i]);
    delete [] idin;
    delete [] foffin;
}
PwstackBinaryFileReader::~PwstackBinaryFileReader()
{
    fclose(fp);
}
SlownessVectorMatrix build_svm_from_buffer(double *buf, int nrow, int ncol)
{
    SlownessVectorMatrix svm(nrow,ncol);
    /* Assume the data are arranged in fortran order with 2 vectors in
       each slot.  */
    int i,j,k;
    for(k=0,i=0;i<nrow;++i)
        for(j=0;j<ncol;++j,k+=2)
        {
            /* Intentionally do not test for index range since 
               svm is constructed above and there is no way it can
               get out of bounds in this procedure */
            SlownessVector sv(buf[k],buf[k+1]);
            svm(i,j)=sv;
        }
    return svm;
}


/* id is the vector index for foffs */
ThreeComponentEnsemble *PwstackBinaryFileReader::read_gather(int id)
{
    ThreeComponentEnsemble *result;
    const string base_error("PwstackBinaryFileReader::read_gather method:  ");
    if(fseek(fp,foffs[id],SEEK_SET))
    {
        stringstream ss;
        ss << base_error <<" fseek to offset "<<foffs[id]<<" failed for "
            << "evid="<<ids[id]<<endl;
        throw SeisppError(ss.str());
    }
    PwstackGatherHeader gh;
    if(fread(&gh,sizeof(PwstackGatherHeader),1,fp)!=1)
    {
        stringstream ss;
        ss << base_error << "fread error while attempting to read "
            <<"gather header for evid="<<ids[id];
        throw SeisppError(ss.str());
    }
    if(id!=gh.sequence_number)
    {
        stringstream ss;
        ss << base_error << "sequence number  mismatch.  "<<endl
            <<"gather header sequence number="<<gh.sequence_number <<endl
            <<"Expected sequence number "<<id<<" for evid="<<ids[id]<<endl;
        throw SeisppError(ss.str());
    }
    /* Now read the block of data containing the matrix of 
       slowness vectors.   The procedure called after fread 
       converts this to a SlownessVectorMatrix object. */
    int gridsize=gh.svmrows*gh.svmcolumns*2; //2 for ux,uy values
    double *svmbuf=new double [gridsize];
    if(fread(svmbuf,sizeof(double),gridsize,fp)!=gridsize)
    {
        delete [] svmbuf;
        stringstream ss;
        ss << base_error << "fread error while reading slowness vector "
            << "matrix section of gather header"<<endl
            << "Error while reading data for evid="<<gh.evid<<endl;
        throw SeisppError(ss.str());
    }
    svm=build_svm_from_buffer(svmbuf, gh.svmrows, gh.svmcolumns);
    delete [] svmbuf;

    try {
        Metadata ghmd=PGHtoMD(gh);
        result=new ThreeComponentEnsemble(ghmd,gh.number_members);
        for(int i=0;i<gh.number_members;++i)
        {
            PwstackTraceHeader th;
            if(fread(&th,sizeof(PwstackTraceHeader),1,fp)!=1)
            {
                stringstream ss;
                ss << base_error <<" fread failed reading trace header for "
                    <<i<<"th member of ensemble of size"<<gh.number_members;
                // Temporarily remove this delete.  Generating a double free 
                // error
                //delete result;
                throw SeisppError(ss.str());
            }
            Metadata trmd=THtoMD(th);
            ThreeComponentSeismogram seis(trmd,false);
            /* Set required parameters.  Some of these may reset attributes 
               set by constructor.  Intentional to make this more stable*/
            seis.dt=1.0/th.samprate;
            seis.t0=th.time;
            seis.ns=th.nsamp;
            seis.tref=relative;  // We assume arrival time reference frame
            seis.live=true;
            /* temporary pointer just to make this less obscure.
               Used to point to first byte of u matrix memory block*/
            double *uptr=seis.u.get_address(0,0);
            int nsamp3c=3*th.nsamp;
            if(fread(uptr,sizeof(double),nsamp3c,fp)!=nsamp3c)
            {
                stringstream ss;
                ss << base_error <<" fread failed reading 3C data samples for "
                    <<i<<"th member of ensemble of size"<<gh.number_members
                    <<"This is data for station name "<<th.sta;
                delete result;
                throw SeisppError(ss.str());
            }
            result->member.push_back(seis);
        }
        return result;
    }catch(...)
    {
        delete result;
        throw;
    };
}
