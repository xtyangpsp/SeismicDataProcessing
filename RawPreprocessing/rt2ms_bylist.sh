#!/bin/bash
# this script runs the rt2ms and miniseed2days program following the file list(zip file/das list).
master_dir=/export/home/field/stationdownloads
miniseed_mst=$1 # for rt2ms
#dayvol_mst=$2
file_list=$2   # for rt2ms
par=$3
#daslist=$3    #for miniseed2days

logfile=rt2ms_by_list.log

echo 'log file' > $logfile

#### for rt2ms
for zipfile in `cat $file_list`

###for miniseed2days. use the following before 'do'
#/bin/ls -d $miniseed_mst/* > das.list

#for das in `cat $daslist`
do
  echo 'run rt2ms for: '$zipfile
  das=`echo $zipfile | sed -e 's/\./ /g' | awk '{print $2}'`
  ##echo hello1
  echo Creat miniseedfiles to: $master_dir/$miniseed_mst/$das
  rt2ms -p $par -L -f $zipfile -o $master_dir/$miniseed_mst/$das

  #echo 'run miniseed2days for msd files of: ' $das >> $logfile
 # miniseed2days -v -u -S $dayvol_mst/$das $miniseed_mst/$das
  
  #zipfile=$das
  echo 'Done: '$zipfile
  echo $zipfile '----->' `date` >> $logfile
done

