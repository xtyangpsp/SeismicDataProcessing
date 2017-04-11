#include "coords.h"
#include "SeisppError.h"
#include "SEEDStaGeom.h"
using namespace SEISPP;
SEEDStaGeom::SEEDStaGeom()
{
    net=string("");
    sta=string("");
    lat=-9999.;
    lon=-9999.;
    elev=-9999.;
}
SEEDStaGeom::SEEDStaGeom(Dbptr db, string netin) : net(netin)
{
    char s[20];
    if(dbgetv(db,0,"sta",s,"lat",&lat,"lon",&lon,"elev",&elev,NULL)
            ==dbINVALID)
        throw SeisppError(string("SEEDStaGeom constructor:  ")
                + "dbgetv error");
    sta=string(s);
}
SEEDStaGeom::SEEDStaGeom(const SEEDStaGeom& parent) 
    : net(parent.net),sta(parent.sta)

{
    lat=parent.lat;
    lon=parent.lon;
    elev=parent.elev;
}
SEEDStaGeom&  SEEDStaGeom::operator=(const SEEDStaGeom& parent) 
{
    if(this!=&parent)
    {
        lat=parent.lat;
        lon=parent.lon;
        elev=parent.elev;
        sta=parent.sta;
        net=parent.net;
    }
    return *this;
}
bool SEEDStaGeom::operator==(const SEEDStaGeom& other)
{
    if(net!=other.net) return false;
    if(sta!=other.sta) return false;
    /* Assume lat and lon are in degrees.  Start with a very rough
       test for difference larger than nominally 1 km */
    if(fabs(lat-other.lat)>0.01) return false;
    if(fabs(lon-other.lon)>0.01) return false;
    /* use 10 m of elevation for test*/
    if(fabs(elev-other.elev)>0.01) return false;
    /* A more precise test to 100 m */
    const double dxmax(0.1);
    /* If we are this close a great circle path is not needed. */
    const double Rearth(6378.164);
    double dx,dy;
    dy=lat-other.lat;
    dy=Rearth*rad(dy);
    dx=lon-other.lon;
    dx=cos(rad(lat))*Rearth*rad(dx);
    double delta=hypot(dx,dy);
    if(delta>dxmax) return false;
    return true;
}
bool SEEDStaGeom::operator!=(const SEEDStaGeom& other)
{
    return(!((*this)==other));
}

ostream& operator<<(ostream& os,SEEDStaGeom& ssg)
{
    os << "Net code="<<ssg.net<<" sta code="<<ssg.sta<<endl
        <<"lat,lon,elev="<<ssg.lat<<" "<<ssg.lon<<" "<<ssg.elev<<endl;
    return(os);
}


