#!/bin/csh
set jday=$1
set year=2015
set db=oiink
foreach sta (`cat rtsta.list`)
  echo Running dbpick for $sta
#  dbpick -sc "${sta}:[VL][HM]." -ts $year${jday}:00:00:00.0 -tw 86400 \
  dbpick -sc "${sta}:*" -ts $year${jday}:00:00:00.0 -tw 86400 \
   -nostarttalk $db
end
