#!/bin/bash
masterdir="/export/home/field/oiink_continuous/wf_XO"
origdir="/export/home/field/stationdownloads/"$1
cd $origdir
flist=$2
for sta in `cat $flist`
#echo $sta
do

cd ./$sta

for year in 2015 # 2014
do
/bin/ls -d $year/* > temp.d  
for yearday in `cat temp.d` 
  do
    echo "Moving data from:" $origdir/$sta/$yearday "to" $masterdir/$yearday
    if [ ! -d $masterdir/$year ]; then mkdir $masterdir/$year; fi
    if [ ! -d $masterdir/$yearday ]; then
    mkdir $masterdir/$yearday
    echo "making directory: "$masterdir/$yearday
    fi
    #echo Move $yearday/* $masterdir/$yearday
    mv -v $yearday/* $masterdir/$yearday
  done
  rm temp.d
done
cd ..
done
