#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <float.h>
#include <string.h>
#include <sstream>
#include <algorithm>
#include "PwmigFileHandle.h"
using namespace std;
/* This appears necessary on quarray to allow large file support */
#define _FILE_OFFSET_BITS 64
/* used to sort records by gridid */
struct comp_gridid : public binary_function<PwmigFileRecord,PwmigFileRecord,bool> {
    bool operator()(PwmigFileRecord x, PwmigFileRecord y)
    {
	if(x.gridid<y.gridid)
		return true;
	else
		return false;
    }
};
/* Function used to construct common name for slowness vector data file */
string build_slowness_file_name(string base)
{
    string name;
    name=base+SlownessFileExtension+"."+dfile_ext;
    return name;
}
/* Read mode constructor */
PwmigFileHandle::PwmigFileHandle(string fname, bool smode, bool cohmode)
{
	const string base_error("PwmigFileHandle constructor:  ");

	/* Save this root name in private area */
    rootname=fname;
    /* assume output will be good by default */
    delete_files_on_exit=false;
    readmode=true;

	/* dfile_ext and hdr_ext come from the PwmigFileHandle.h include file */
	string dfile=fname+"."+dfile_ext;
	string hfile=fname+"."+hdr_ext;
	scalar_mode=smode;

	datafd=open(dfile.c_str(),O_RDONLY);
	if(datafd<0)
		throw SeisppError(base_error
			 + string("Cannot open input data file=")
			 + dfile);
	hdrfd=open(hfile.c_str(),O_RDONLY);
	if(hdrfd<0)	
	{
		close(datafd);
		throw SeisppError(base_error
			 + string("Cannot open input header data file=")
			 + dfile);
	}
	/* First read the global data block */	
	ssize_t read_count;
	read_count=read(hdrfd,static_cast<void *>(&filehdr),sizeof(PwmigFileGlobals));
	if(read_count!=sizeof(PwmigFileGlobals))
			throw SeisppError(base_error
			 + string("fread error on file = ")+hfile
			 + string("\nError on first read attempt to read globals parameters") );
	PwmigFileRecord thisrec;
	/* Load up all the header data */
	while(read(hdrfd,static_cast<void *>(&thisrec),sizeof(PwmigFileRecord))>0)
	{
		recs.push_back(thisrec);
	}
	if(recs.size()<=0) throw SeisppError(base_error
		+ string("fread error reading header data from file=")
		+hfile);
	/* Sort the recs list by gridid so we can read in data for each plane
		wave component together */
	sort(recs.begin(),recs.end(),comp_gridid());
	current_record=recs.begin();
        if(cohmode)
            svmfp=NULL;
        else
        {
    	/* In read mode we need this slowness data file open */
            string slownessfname=build_slowness_file_name(rootname);
            svmfp=fopen(slownessfname.c_str(),"r");
            if(svmfp==NULL)
                throw SeisppError(base_error + "fopen failed for file="
                    + fname);
        }
}
/* Write mode constructor */
PwmigFileHandle::PwmigFileHandle(string fname, bool smode,RectangularSlownessGrid& ug)
{
	const string base_error("PwmigFileHandle constructor:  ");
	readmode=false;

    /* Save this root name in private area */
    rootname=fname;
    /* assume output will be good by default */
    delete_files_on_exit=false;
	/* dfile_ext and hdr_ext come from the PwmigFileHandle.h include file */
	string dfile=fname+"."+dfile_ext;
	string hfile=fname+"."+hdr_ext;
	scalar_mode=smode;
	/* Copy in needed ug data.   A bit complicated by need to require ug geometry to 
	have a zero delta u position - error thrown to abort pwmig in that situation. */
	strncpy(this->filehdr.gridname,ug.name.c_str(),16);  //16 is the size of name - maintenance warning
	/* uxlow and uylow are always negative so this formula is a bit odd.   
	Computation is for C indexing starting at 0 */
	filehdr.i0=SEISPP::nint(-ug.uxlow/ug.dux);
	filehdr.j0=SEISPP::nint(-ug.uylow/ug.duy);
	/* using float eps is a simple test because slowness is of order 1 */
	double dux0test,duy0test;
	dux0test=ug.uxlow+ug.dux*((double)(filehdr.i0));
	duy0test=ug.uylow+ug.duy*((double)(filehdr.j0));
	if( (fabs(dux0test)>FLT_EPSILON) 
		    || (fabs(duy0test)>FLT_EPSILON) )
	{
		throw SeisppError(base_error
		  + "Unacceptable Slowness_Grid_Definition parameters\n"
		  + "Grid must hit 0,0 slowness vector delta - ulow+i*du must = 0 for some i\n");
	}
	
	mode_t fumask=0775;
	/* In output mode all we do is open the file in write mode */
	datafd=open(dfile.c_str(),
			O_WRONLY | O_CREAT | O_TRUNC, fumask);
	if(datafd<0)
		throw SeisppError(base_error
			 + string("Cannot open output data file=")
			 + dfile);
	hdrfd=open(hfile.c_str(),
			O_WRONLY | O_CREAT | O_TRUNC,fumask);
	if(hdrfd<0)
	{
		close(datafd);
		throw SeisppError(base_error
			 + string("Cannot open output header data file=") );
	}
}
PwmigFileHandle::~PwmigFileHandle()
{
    const string base_error("PwmigFileHandle destructor:  ");
	close(datafd);
	if(!readmode)
	{
            if(delete_files_on_exit)
            {
                string dfile=rootname+"."+dfile_ext;
                if(remove(dfile.c_str()))
                    throw SeisppError(base_error + "stdio C function remove failed while "
                            + "trying to remove file="+dfile);
            }
            else
            {
                string message=base_error + "error dumping output header data";
		ssize_t test;
		test=write(hdrfd,static_cast<void *>(&filehdr),sizeof(PwmigFileGlobals) );
		if(test!=sizeof(PwmigFileGlobals) )
		{
			close(hdrfd);
                        ostringstream ess;
                        ess << message <<endl
                            << "Tried to write PwmigFileGlobals block of size "
                            << sizeof(PwmigFileGlobals) <<endl
                            << "binary write function returned count of "
                            << test << "  Mismatch is fatal error."<<endl;
			throw SeisppError(ess.str());
		}
		for(current_record=recs.begin();current_record!=recs.end();++current_record)
		{
			test=write(hdrfd,static_cast<void *>(&(*current_record)),
				sizeof(PwmigFileRecord) );
			if(test!=sizeof(PwmigFileRecord))
			{
				close(hdrfd);
                        ostringstream ess;
                        ess << message <<endl
                            << "Tried to write PwmigFileRecord block of size "
                            << sizeof(PwmigFileRecord) <<endl
                            << "binary write function returned count of "
                            << test << "  Mismatch is fatal error."<<endl;
				throw SeisppError(message);
			}
		}
            }
	}
	else if(svmfp!=NULL)
	/* In read mode the slowness data file is kept open until destruction here */
		fclose(svmfp);
}
template <class T> void PwmigFileRecord_load(T& ts, 
	PwmigFileRecord& hdr) 
{
	try {
		hdr.gridid=ts.get_int("gridid");
		hdr.ix1=ts.get_int("ix1");
		hdr.ix2=ts.get_int("ix2");
		hdr.ux0=ts.get_double("ux0");
		hdr.uy0=ts.get_double("uy0");
		hdr.ux=ts.get_double("ux");
		hdr.uy=ts.get_double("uy");
		hdr.elev=ts.get_double("elev");
		hdr.nsamp=ts.ns;
		hdr.dt=ts.dt;
		hdr.t0=ts.t0;
	} catch (MetadataGetError mderr)
	{
		throw mderr;
	}
}
/* complicated error function needed for errors from lseek */
string buildlseek_message(off_t foff)
{
	char buf[256];
	ostringstream ss(buf);
	ss<<"PwmigFileHandle::save method:  lseek error"<<endl
		<<"lseek returned off_t="<<foff <<endl;
	switch (errno)
	{
	case EBADF:
		ss <<"errno = EBADF - not an open file descriptor";
		break;
	case ESPIPE:
		ss <<"errno=ESPIPE - filedes is a pipe, socket of FIFO";
		break;
	case EINVAL:
		ss<<"errno=EINVAL - whence arguement invalid";
		break;
	case EOVERFLOW:
		ss << "errno=EOVERFLOW - file offset exceeds word length"
		  << endl<<"Probably a large file error";
		break;
	default:
		ss << "Unknown value for errno - "
			<< "program bug has probably overwritten data segment";
	}
	return(ss.str());
}
void PwmigFileHandle::save(TimeSeries& ts)
{
	PwmigFileRecord hdr;
	try {
		PwmigFileRecord_load<TimeSeries>(ts,hdr);
	} catch (MetadataGetError mderr)
	{
		throw mderr;
	}
	/* Load the current file position as foff for read method */
	hdr.foff=lseek(datafd,0,SEEK_CUR);
	if(hdr.foff<0)
	{
		string errmess=buildlseek_message(hdr.foff);
		throw SeisppError(errmess);
	}
	/* Assume a seek to end is not necessary and all saves will be appends */
	ssize_t test;
	void *sptr=static_cast<void *>(&(ts.s[0]));
	test=write(datafd,sptr,(ts.ns)*sizeof(double));
	if(test!=(ts.ns)*sizeof(double) )
        {
            ostringstream ess;
            ess << "PwmigFileHandle::save(scalar data):  write error" <<endl
                << "Tried to write "<< (ts.ns)*sizeof(double) << "bytes" <<endl
                << "write function return count = "<<test<<endl
                << "Mismatch is a fatal error"<<endl;
		throw SeisppError(ess.str());
        }
	recs.push_back(hdr);
}
void PwmigFileHandle::save(ThreeComponentSeismogram& tcs)
{
	PwmigFileRecord hdr;
	try {
		PwmigFileRecord_load<ThreeComponentSeismogram>(tcs,hdr);
	} catch (MetadataGetError mderr)
	{
		throw mderr;
	}
	/* Load the current file position as foff for read method */
	hdr.foff=lseek(datafd,0,SEEK_CUR);
	if(hdr.foff<0)
	{
		string errmess=buildlseek_message(hdr.foff);
		throw SeisppError(errmess);
	}
	recs.push_back(hdr);
	/* similar to above, but here we use the fact that in this implementation
	the data in a 3c seismogram is a contiguous block 3xnsamp long */
	ssize_t test;
	double *sptr=tcs.u.get_address(0,0);
	test=write(datafd,static_cast<void *>(sptr),
		3*sizeof(double)*tcs.ns);
	if(test!=(3*sizeof(double)*tcs.ns) )
        {
            ostringstream ess;
            ess << "PwmigFileHandle::save(3Cdata):  write error" <<endl
                << "Tried to write "<< 3*(tcs.ns)*sizeof(double) << "bytes" <<endl
                << "write function return count = "<<test<<endl
                << "Mismatch is a fatal error"<<endl;
		throw SeisppError(ess.str());
        }
}
/* This is the reciprocal of PwmigFileRecord_load above */
template <class T> void PwmigFileRecord_load_Metadata(PwmigFileRecord& hdr, T& d) 
{
	d.put("gridid",hdr.gridid);
	d.put("ix1",hdr.ix1);
	d.put("ix2",hdr.ix2);
	d.put("ux0",hdr.ux0);
	d.put("uy0",hdr.uy0);
	d.put("ux",hdr.ux);
	d.put("uy",hdr.uy);
	d.put("elev",hdr.elev);
	d.put("nsamp",hdr.nsamp);
	d.ns=hdr.nsamp;
	d.t0=hdr.t0;
	d.dt=hdr.dt;
}
/* similar for global parameters */
template <class T> void PwmigFileHandle_load_global(PwmigFileGlobals& g,T *d)
{
	d->put("evid",g.evid);
	d->put("origin.lat",g.slat);
	d->put("origin.lon",g.slon);
	d->put("origin.depth",g.sdepth);
	d->put("origin.time",g.stime);
}
/* companion to immediately below.  Loads on seismogram from file handle record */
ThreeComponentSeismogram *load_3c_seis(PwmigFileRecord& hdr, int fd)
{
	ThreeComponentSeismogram *seis=new ThreeComponentSeismogram(hdr.nsamp);
	PwmigFileRecord_load_Metadata<ThreeComponentSeismogram>(hdr,*seis);
	lseek(fd,hdr.foff,SEEK_SET);
	ssize_t test;
	test=read(fd,static_cast<void *>(seis->u.get_address(0,0)),
			(3*hdr.nsamp)*sizeof(double));
	if(test!=(3*hdr.nsamp)*sizeof(double) )
        {
            ostringstream ess;
            ess << "PwmigFileHandle load_3c_seis procedure:  read error"
                <<endl
                << "Tried to read "<< (3*hdr.nsamp)*sizeof(double) << " bytes"
                <<endl
                << "read function return value = "<<test
                << endl
                << "Mismatch is a fatal error"<<endl;
		throw SeisppError(ess.str());
        }
	return seis;
}
TimeSeries *load_seis(PwmigFileRecord& hdr, int fd)
{
	TimeSeries *seis=new TimeSeries(hdr.nsamp);
	PwmigFileRecord_load_Metadata<TimeSeries>(hdr,*seis);
	lseek(fd,hdr.foff,SEEK_SET);
	ssize_t test;
	double *ptr;
	ptr=&(seis->s[0]);
	test=read(fd,static_cast<void *>(ptr),(hdr.nsamp)*sizeof(double));
	if(test!=(hdr.nsamp)*sizeof(double) )
        {
            ostringstream ess;
            ess << "PwmigFileHandle load_seis procedure:  read error"
                <<endl
                << "Tried to read "<< (hdr.nsamp)*sizeof(double) << " bytes"
                <<endl
                << "read function return value = "<<test
                << endl
                << "Mismatch is a fatal error"<<endl;
		throw SeisppError(ess.str());
        }
	return seis;
}
ThreeComponentEnsemble *PwmigFileHandle::load_next_3ce()
{
	if(scalar_mode) throw SeisppError(
		string("PwmigFileHandle::load_next:  attempt to load ")
		+ string("three component data from a scalar file\n")
		+ string("Coding error") );
	if(current_record==recs.end()) return(NULL);
	ThreeComponentEnsemble *ens=new ThreeComponentEnsemble();
	PwmigFileHandle_load_global<ThreeComponentEnsemble>(filehdr,ens);
	int count;
	int current_gridid;
	double ux,uy;
	current_gridid=current_record->gridid;
	ThreeComponentSeismogram *tcs;
	count=0;
	do {
		// load the ensemble metadata on the first pass
		if(count==0)
		{
			ens->put("gridid",current_gridid);
			ens->put("ux",current_record->ux);
			ens->put("uy",current_record->uy);
			ens->put("ux0",current_record->ux0);
			ens->put("uy0",current_record->uy0);
		}
		try {
			tcs=load_3c_seis(*current_record,datafd);
		} catch (SeisppError serr)
		{
			serr.log_error();
			cerr << "Problems reading data for gridid="<<current_gridid<<endl
				<< "Grid index values="<<current_record->ix1
				<< " " << current_record->ix2<<endl;
			delete tcs;
		}
		tcs->live=true;
		tcs->tref=relative;
		ens->member.push_back(*tcs);
		delete tcs;
		++current_record;
		++count;
	} while( (current_record->gridid == current_gridid)
		&& (current_record!=recs.end()) );
	return(ens);
}
TimeSeriesEnsemble *PwmigFileHandle::load_next_tse()
{
	if(!scalar_mode) throw SeisppError(
		string("PwmigFileHandle::load_next:  attempt to load ")
		+ string("scalar data from a three component file\n")
		+ string("Coding error") );
	if(current_record==recs.end()) return(NULL);
	TimeSeriesEnsemble *ens=new TimeSeriesEnsemble();
	PwmigFileHandle_load_global<TimeSeriesEnsemble>(filehdr,ens);
	int count;
	int current_gridid;
	double ux,uy;
	current_gridid=current_record->gridid;
	TimeSeries *ts;
	count=0;
	do {
		// load the ensemble metadata on the first pass
		if(count==0)
		{
			ens->put("gridid",current_gridid);
			ens->put("ux",current_record->ux);
			ens->put("uy",current_record->uy);
			ens->put("ux0",current_record->ux0);
			ens->put("uy0",current_record->uy0);
		}
		try {
			ts=load_seis(*current_record,datafd);
		} catch (SeisppError serr)
		{
			serr.log_error();
			cerr << "Problems reading data for gridid="<<current_gridid<<endl
				<< "Grid index values="<<current_record->ix1
				<< " " << current_record->ix2<<endl;
			delete ts;
		}
		ts->live=true;
		ts->tref=relative;
		ens->member.push_back(*ts);
		delete ts;
		++current_record;
		++count;
	} while( (current_record->gridid == current_gridid)
		&& (current_record!=recs.end()) );
	return(ens);
}
/* This was added Jan 2015 to save SlownessVector data - a 4D array */
void PwmigFileHandle::save_slowness_vectors(SlownessVectorMatrix& u0,
                RectangularSlownessGrid& ugrid)
{
    const string base_error("PwmigFileHandle::save_slowness_vector:  ");
    /* This method is kind of showhorned into this object, but this is 
       a more maintainable approach than carrying around a file name and
       another set of read/write procedures.   First we create a file name
       from the base name and then try to open such a file in write mode. */
    int i,j,k,l;
    string fname=build_slowness_file_name(rootname);
    FILE *fp=fopen(fname.c_str(),"w");
    if(fp==NULL)
        throw SeisppError(base_error + "fopen failed for file="
                + fname);
    /* First record records sizes */
    int dims[4];
    dims[0]=ugrid.nux;
    dims[1]=ugrid.nuy;
    dims[2]=u0.rows();
    dims[3]=u0.columns();
    if(fwrite(dims,sizeof(int),4,fp)!=4)
    {
        fclose(fp);
        throw SeisppError(base_error + "fwrite failed writing first data block"
                + "(dimensions)");
    }
    /* Now we write the slowness vectors out 2 doubles at a time.
       Order should be clear from the for loops */
    double ubuf[2];
    for(i=0;i<ugrid.nux;++i)
        for(j=0;j<ugrid.nuy;++j)
{
//DEBUG
//cout << "ugrid i,j="<<i<<","<<j<<endl;
            for(k=0;k<u0.rows();++k)
                for(l=0;l<u0.columns();++l)
                {
                    SlownessVector uincident=u0(k,l);
                    SlownessVector du=ugrid.slow(i,j);
                    /* 0 is ux, 1 is uy */
                    ubuf[0]=uincident.ux + du.ux;
                    ubuf[1]=uincident.uy + du.uy;
//cout << k <<","<<l<<" "<<uincident.ux<<" "<<uincident.uy<<" "<<du.ux<<" "<<du.uy<<endl;
                    if(fwrite(ubuf,sizeof(double),2,fp)!=2)
                    {
                        fclose(fp);
                        throw SeisppError(base_error
                                + "fwrite error saving slowness vector data");
                    }
                }
//DEBUG - remove when finished
}
    fclose(fp);
}
SlownessVectorMatrix PwmigFileHandle::plane_wave_slowness_vectors(int iux,
        int iuy)
{
    const string base_error("PwmigFileHandle::plane_wave_slowness_vectors(int,int) mehthod:  ");
    if(svmfp==NULL)
    {
        cerr << base_error
             << "coding error.   svmfp is NULL.  This method was "
             << "probably called for a handl pointing at coherence files"
             <<endl
             << "Exiting - FATAL ERROR"<<endl;
         exit(-1);
    }
    int psgn1, psgn2;   // pseudostation grid dimensions
    int ugn1, ugn2;   // slowness grid dimensions
    size_t count(0);
    rewind(svmfp);
    count += fread((void *)(&ugn1),sizeof(int),1,svmfp);
    count += fread((void *)(&ugn2),sizeof(int),1,svmfp);
    count += fread((void *)(&psgn1),sizeof(int),1,svmfp);
    count += fread((void *)(&psgn2),sizeof(int),1,svmfp);
    if(count!=4)
    {
        throw SeisppError(base_error
            + "fread error while reading array dimensions at head of file");
    }
    int gridid=ugn1*iux + iuy;
    try {
	   return (this->plane_wave_slowness_vectors(gridid));
    }catch(...){throw;};
}
/* This is the main read routine for the slowness grid data file.   Others are 
just wrappers on this one */
SlownessVectorMatrix PwmigFileHandle::plane_wave_slowness_vectors(int gridid)
{
   string base_error("PwmigFileHandle::plane_wave_slowness_vectors:  ");
    if(svmfp==NULL)
    {
        cerr << base_error
            << "coding error.   svmfp is NULL.  This method was "
            << "probably called for a handl pointing at coherence files"
            <<endl
            << "Exiting - FATAL ERROR"<<endl;
        exit(-1);
    }
    string fname=build_slowness_file_name(rootname);  //not needed but useful for errors
    int psgn1, psgn2;   // pseudostation grid dimensions
    int ugn1, ugn2;   // slowness grid dimensions
    size_t count(0);
    rewind(svmfp);
    count += fread((void *)(&ugn1),sizeof(int),1,svmfp);
    count += fread((void *)(&ugn2),sizeof(int),1,svmfp);
    count += fread((void *)(&psgn1),sizeof(int),1,svmfp);
    count += fread((void *)(&psgn2),sizeof(int),1,svmfp);
    if(count!=4)
    {
        throw SeisppError(base_error
            + "fread error while reading array dimensions at head of file");
    }
   /* The gridid is the count (starting at 0) working through the slowness grid in C order 
      (column index runs fastest).  The slowness data file is stored in the order
      pwmig needs it, which is the gridid components are the slowest varying index.  
      We can then compute the file offset to the start of the needed data by the following
      simple formula.*/
   long pwblocksize=2*sizeof(double)*(psgn1)*(psgn2);
    
   long foff=gridid*pwblocksize + 4*sizeof(int);   // this is foff from file start
   if(fseek(svmfp,foff,SEEK_SET)<0) 
   {
       stringstream ss;
       ss << base_error << "fseek failed of file "<<fname<<endl
           << "Attempted to seek to file offset="<<foff<<endl;
       throw SeisppError(ss.str());
   }
   SlownessVectorMatrix result(psgn1,psgn2);
   int i,j;
   count=0;
   double u[2];  // read vectors one at a time into this array
   for(i=0;i<psgn1;++i)
       for(j=0;j<psgn2;++j)
       {
           count = fread((void *)u,sizeof(double),2,svmfp);
           if(count != 2)
           {
               stringstream ss;
               ss << base_error << "fread error reading block for gridid="<<gridid<<endl
                   << "Reading data at pseudostation grid index ("<<i<<","<<j<<")"<<endl;
               throw SeisppError(ss.str());
           }
           SlownessVector slow(u[0],u[1]);
           result(i,j)=slow;
       }
   return result;
}
SlownessVectorMatrix PwmigFileHandle::incident_wave_slowness_vectors()
{
    try {
        return(this->plane_wave_slowness_vectors(this->filehdr.i0,
                    this->filehdr.j0));
    }catch(...){throw;};
}

