#!/bin/csh
foreach net (XO)
foreach year (2014)
miniseed2db wf_$net/$year/??? oiinkc_XO_2014
end
end
