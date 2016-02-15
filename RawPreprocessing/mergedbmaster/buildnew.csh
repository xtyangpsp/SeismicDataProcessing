#!/bin/csh
\rm oiink.*
\rm -r response nom_response
dbbuild -b oiink OIINK_dbbuild_batch
