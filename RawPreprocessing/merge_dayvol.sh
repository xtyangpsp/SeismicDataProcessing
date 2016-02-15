#!/bin/bash
masterdir="/export/home/field/oiink_continuous/wf_TA"
origdir="/export/home/field/oiink_continuous/tmp_othernetwork/TA/TA_late_2012_dayvol"
cd $origdir
for year in 2012
do
  for yearday in $year/???
  do
    echo "Moving data from:" $origdir/$yearday "to" $masterdir/$yearday
    if [ ! -d $masterdir/$yearday ]; then
    mkdir $masterdir/$yearday
    echo "making directory: "$masterdir/$yearday
    fi
    mv $yearday/* $masterdir/$yearday
  done
done
