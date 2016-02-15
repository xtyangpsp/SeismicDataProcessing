#!/bin/csh
set masterdir=/export/home/field/oiink_continuous/wf_XO_rt
set sta=$1
set year=2013
cd ~/oiink/db/$year
foreach day (???)
echo Copying $sta data for day $year/$day
mkdir $masterdir/$year/$day
cp $day/$sta* $masterdir/$year/$day
end
