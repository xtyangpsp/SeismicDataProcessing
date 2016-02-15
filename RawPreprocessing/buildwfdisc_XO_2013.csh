#!/bin/csh
foreach net (XO)
foreach year (2013)
miniseed2db wf_$net/$year/??? oiinkc_XO_2013
end
end
