#!/bin/csh
set jday=$1
set db=oiink
foreach sta (`cat rtsta.list`)
  echo Running dbpick for $sta
  dbpick -sc "${sta}:VM." -ts 2013${jday}:00:00:00.0 -tw 86400 \
   -nostarttalk $db
end
