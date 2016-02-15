#!/bin/csh
foreach net (IU NM TA XO XO_rt)
foreach year (2011 2012 2013)
miniseed2db wf_$net/$year/??? oiinkc
end
end
