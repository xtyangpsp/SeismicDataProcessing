#!/bin/csh
# This script is used to build the databases for other networks by 
#running seed2db. TA db is treated separately.
#CAUTION:
# Before running this script, make sure there are no existing databases 
#under netdb/net, or say, netdb has ONLY TA sub- directory.
# added N4 network, Dec 16, 2014 by Xiaotao Yang

set outdir_master=netdb
#
#

# now crack each network's set of dataless files 
# copied and modified from GLP's script.
foreach net (IU NM N4)
  \rm -rf $outdir_master}/${net}
  mkdir ${outdir_master}/${net}

  foreach f (dataless/${net}/*.dataless)
    seed2db -respdir response -stagedir response/stage $f $outdir_master/$net/$net
  end
end

