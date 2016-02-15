#!/bin/bash
#Description:
#
#USAGE: excludelist.sh mstlist excludelist
#
#Xiaotao Yang (xtyang@indiana.edu)	Tue Apr 2 15:34:11 EDT 2013	Indiana University
# 1001 E 10TH ST., BLOOMINGTON, IN 47405
mstlist=$1
excludelist=$2

for flist in `cat $mstlist`
do
	num=`cat $excludelist | grep -c $flist`
	if [ $num == 0 ]; then
		echo $flist
	fi
done
