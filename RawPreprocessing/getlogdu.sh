#!/bin/bash
# this script is used to get disk-usage from rt log files for each station.
# Xiaotao Yang on Jan 2, 2014

if [ $# -lt 1 ]
	then
		echo "USAGE: getlogdu.sh logdirlistfile"
		exit 1
fi

dirlist=$1

echo "STA  DISK#  USED(Mb)  AVAILABLE(Mb)  TOTAL(Mb)  USED(%)"

for dir in `cat $dirlist`
do
	cat $dir/log | grep DISK | grep USED | tail -2 | \
		awk '{printf "%5s  %5s  %8d  %8d  %8d  %6d\n", $4,$8, $10/1024,$12/1024, \
		$14/1024, 100 - 100*$12/$14}' | sed -e 's/://g'
done

