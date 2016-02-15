#!/bin/bash
#Description:
# Get das list for the corresponding *.ZIP/*.tar data set. And match the das list with corresponding stations.
# Use output: das_station.list.
#
#USAGE example:  ./match_das_sta.sh -p oiink.par -f data_file_list
# 			-p: parameter file includes das and station list, I recommend to use the same parameter file as rt2ms.
#			-f: name list file includes all of the data packages.
#
#Xiaotao Yang (xtyang@indiana.edu)	
#Tue Oct 2-13 15:34:34 EDT 2012
#Indiana University
#1001 E 10TH ST., BLOOMINGTON, IN 47405
#

if [ $# == 0 ]||[ $1 == "-h" ]
	then
	  echo '************'
echo 'USAGE example:  ./match_das_sta.sh -p oiink.par -f data_file_list'
echo ' 		-p: parameter file includes das and station list, I recommend to use the same parameter file as rt2ms.'
echo '		-f: name list file includes all of the data packages.'
#echo 'Output file: das_station.list'
echo ' '
echo 'Xiaotao Yang (xtyang@indiana.edu)	'
echo 'Tue Oct 2-13 15:34:34 EDT 2012'
echo 'Indiana University'
echo '1001 E 10TH ST., BLOOMINGTON, IN 47405'
	  echo '************'
	  exit 1
fi


#fileout='das_station.list'

if [ $1 == '-p' ]
	then
	  parfile=$2
	  if [ $3 == '-f' ]
	  	then
	  	  listfile=$4
	  else
	    echo 'Please specify the name list file! ERROR!'
	    exit 1
	  fi
elif [$1 == '-f' ]
	then
	  listfile=$2
	  if [ $3 == '-p' ]
	  	then
	  	  parfile=$4
	  else
	    echo 'Please specify the parameter file! ERROR!'
	    exit 1
	  fi
else
  echo 'ERROR! Use option -h to get usage information!'
  exit 1
fi

sed -e 's/\./ /g' $listfile | awk '{print $2}' |sort > das_list.temp
i=0
j=0
temp_das='++++'

#echo 'Matching das and stations for the data files ...'

echo '>> DAS and station pairs with data:' # > $fileout

echo DAS     Station #>> $fileout

for das in `cat das_list.temp`
do
#echo $das   $temp_das
    if [ $das != $temp_das ]
	then
      		#echo $das	`grep -iw -m 1 $das $parfile |awk '{print $5}'|sed -e 's/;/ /g'` >> $fileout
      		echo $das	`grep -iw -m 1 $das $parfile |awk '{print $5}'|sed -e 's/;/ /g'` 
      		((i=i+1))
      		temp_das=$das
    fi
    ((j=j+1))
done    

# das list not included in the data list.
echo ' ' # >> $fileout
echo '>> Stations have no data' #>> $fileout

awk '{if ( NR >1 ) print $1}' $parfile > das_all.temp
temp_das='++++'


#echo 'Searching sations have no data ...'

k=0
p=0
for das in `cat das_all.temp | sed -e 's/\;/ /g'`
do
  if [ $das != $temp_das ]
    then
  #echo $das
  count=`grep -c $das das_list.temp`
  #echo $count
  #exit 1
    if [ $count == 0 ]
      then
        echo $das	`grep -iw -m 1 $das $parfile |awk '{print $5}'|sed -e 's/;/ /g'` #>>$fileout
        ((k=k+1))
    fi
    ((p=p+1))
  fi
  temp_das=$das
done

# summary
echo ' ' #>> $fileout
#echo '+++++++++++++++++++++++ SUMMARY ++++++++++++++++++++++++++ ' >> $fileout
echo '+++++++++++++++++++++++ SUMMARY ++++++++++++++++++++++++++ ' 

#echo 'TOTAL NUMBER of stations in par file: 	'$p >> $fileout
echo 'TOTAL NUMBER of stations in par file: 	'$p

#echo 'DATA FILES in the list: 		'$j >>$fileout
echo 'DATA FILES in the list: 		'$j

#echo 'Number of stations HAVE data: 		'$i >> $fileout
echo 'Number of stations HAVE data: 		'$i

#echo 'Number of stations have NO data: 	'$k >> $fileout
echo 'Number of stations have NO data: 	'$k

#echo '>> Output file:  '$fileout
rm das_all.temp das_list.temp

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
