#!/bin/csh
#
# This is a master build script to construct a full dbmaster set of tables for the oiink 
# experiment.   It is designed to make it easy to rebuild dbmaster with the addition of
# new stations.   It is not user fiendly, however, so beware the following:
#  1.  Make sure you remove the dbmaster directory here or move it away before running this
#      script.  It will create the directory and write into a new set each time.   
#  2.  All the dataless SEED files used to build the master are in the dataless directory.
#      If new network codes are added please make sure they are placed in the appropriate
#      subdirectory there.  These files must end in ".dataless" for this script to work 
#      properly.
#  3.  The script tacitly assumes each of the dataless files are othogonal.  i.e. trouble
#      is guaranteed if you have two dataless files with duplicate or conflicting information
#      on the same station.
#  4.  The oiink database is used as the starting database.  Note that the first step is to
#      dbcp the oiink dbmaster to the merged database directory. The concept is we will run
#      this script every time we update the oiink dbmaster.
#
#  IMPORTANT final note:  this is the working copy.  To install dbcp the database in dbmaster here
#  to the work area.   Wise to pause stop orb2db processes before doing that though.   
#
# ~~~~~~
# modified by Xiaotao Yang on Jan 28, 2013 to seperately process the TA dbmaster by firstly merge TA with oiink 
#and then merge IU & NM into it.
# modified by XTY on Jan 29, 2013: merge databases for IU, NM, and TA all.
# modified by glp on June 7, 2013 changing confusing (to me) arrangment of files coming from other places.
# oiink master is not below this working directory and we doing have the dbcp that created problems 
# creating a "portable" dbmaster database.
#
# CAUTION:
# Before running this script, make sure you have had the databases for each net under netdb. Usually this 
#could be done by running ./build_netdb.csh.
# NOTE on N4 collision issue with TA:
# See the readme file in netdb/N4 directory for details. This script generates a readme file on this too.

set outdir=dbmaster
set outdb=oiink
set oiink_original_dir=dbm_oiink_only
set netdbdir=netdb   #=TAdbmaster/TA4OIINK
#
# clear the output
#
\rm -r $outdir
mkdir $outdir
#
# First we copy the oiink master 
# IMPORTANT:  this assumes db name is $outdb
#
cp -r ${oiink_original_dir}/*  $outdir
#

#generates N4 & TA collision readme file.
echo "Be ware of the following issue if you see those stations in the database" > $outdir/DO_NOT_REMOVE_N4TAISSUE.readme
cat $netdbdir/N4/DO_NOT_REMOVE.readme >> $outdir/DO_NOT_REMOVE_N4TAISSUE.readme
# now crack each network's set of dataless files 
#
foreach net (IU NM TA N4)
  dbmerge -v $netdbdir/$net/$net $outdir/$outdb
end
