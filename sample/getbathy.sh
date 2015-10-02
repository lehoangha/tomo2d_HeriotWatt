#!/bin/csh -f
# getbathy.sh

awk 'NR==2{for(i=1;i<=NF;i++) print $i}' smesh.dat > tmp.x
awk 'NR==3{for(i=1;i<=NF;i++) print $i}' smesh.dat > tmp.t
paste tmp.x tmp.t > tmp.xt
