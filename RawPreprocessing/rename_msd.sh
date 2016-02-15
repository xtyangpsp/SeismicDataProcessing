#!/bin/bash
/bin/ls *.msd > flist
for f in `cat flist`
do
	ff=`echo $f | sed -e 's/Mar2015/May2_2015/g' `
	mv -v $f $ff
done
rm flist
