#!/usr/bin/env python
PROG_VERSION = "2009.005"
PROG_NAME = "CATHUNT"

#program to
# find all "LOG:" lines in a rt130 antelope log file and concatinate
# all logs into one file by station

import os.path
import pickle
import re
from string import split, join, replace
from sys import exit
from commands import getoutput


def mkloglist(arg,dirname,names):
 for name in names:
  if name == "log":
   arg.append(dirname+"/"+ name)
 arg.sort()

#find year from path
year = os.path.split(os.path.abspath('.'))
#year = year[-1]
# hack by glp for dec 2013
year = 2013
print "working on year: ", year

#place to put new log files
putdir = "stationlogs"
if os.path.isdir(putdir):
 print "Writing log files in dir: ", putdir
else:
 print "dir \"%s\" does not exist! bye" % putdir
 exit()

#get todays date NOT to process
julout = getoutput("julday")
re_julday = re.compile( r'Julian Date\s*(\d*)\s*(\d*)')
daymatch = re_julday.search(julout)
if daymatch:
 today = daymatch.group(1)
 toyear = daymatch.group(2)
 print "Will not process files on or after today: %s:%s" % (toyear, today)
 if int(toyear) > int(year):
  today = "367"
else:
 print "Can not determine julday! bye."
 exit()


#if has been run before open filelist from pickle
try:
 fpp = open( putdir+"/.cathunt.txt",'r') 
 p = pickle.Unpickler( fpp)
 donefilelist = p.load()
 fpp.close()
except Exception, e:
 donefilelist = []

#walk the tree finding logs to be filed
loglist = []
os.path.walk('.',mkloglist,loglist)

print "found %d antelope log files to process" % (len(loglist))
print

# loop on each log file
for daylog in loglist:
 if daylog in donefilelist:
  print "skipping " +daylog +": has already been processed"
  continue


 parts  = split(daylog,"/" )
 day =parts[1]
 station = parts[2]
 putpath = putdir + '/' + station + '.log'

 if day >= today:
  print "skipping " +daylog+": on or after today: " + today
  continue

 fpin = open(daylog,'r')

 if os.path.isfile(putpath) == False:
  #need to make new file
  fpout = open(putpath, 'w')
  
  timestr = year[2:] + ":" + day + ":" + replace( split( fpin.readline())[1] ,'.',':',1 )[:-1]

  header = "State of Health  %s   ST: %s\n"% \
		(timestr, station)

  print putpath + "Does not exist creating new file with header:\n   " + header[:-1]
  fpout.write(header)
  fpout.close()

 fpout = open(putpath, 'a')

 fpin.seek(0)


 print "Appending " + daylog + " to " +putpath

 count = 0
 for line in fpin:
  where = line.find("LOG:")
  if where != -1:
   fpout.write("%s\n"%line[where+5:].strip())
   count += 1

 print "Wrote %d \"LOG:\" lines" % count 

 fpin.close()
 fpout.close()
 donefilelist.append(daylog)
 print

fpp = open( putdir+"/.cathunt.txt",'w+') 
p = pickle.Pickler( fpp)
p.dump(donefilelist)
fpp.close()
