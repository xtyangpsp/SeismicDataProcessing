#!/bin/csh
foreach net (XO)
foreach year (2015)
miniseed2db wf_$net/$year/??? oiinkc_XO_2015
end
end
