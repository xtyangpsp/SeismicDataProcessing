#include <string>
#include "db.h"
using namespace std;
class SEEDStaGeom
{
public:
    string net;
    string sta;
    /*Here these are stored in degrees */
    double lat,lon;
    double elev;  /* elev in km */ 
    /* Need this explicitly it seems */
    SEEDStaGeom();
    /* Primary contructor 
       db must point to a row in a site table
       netin - network code that will be assigned */
    SEEDStaGeom(Dbptr db, string netin);
    SEEDStaGeom(const SEEDStaGeom& parent);
    SEEDStaGeom& operator=(const SEEDStaGeom& parent);
    bool operator==(const SEEDStaGeom& other);
    bool operator!=(const SEEDStaGeom& other);
    friend ostream& operator<<(ostream& os,SEEDStaGeom& ssg);
};
