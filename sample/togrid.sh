#!/bin/csh -f
# togrid.sh

set GMTsurface = /usr/local/bin/surface

set xyz=$1 Xmin=$2 Xmax=$3 Zmin=$4 Zmax=$5 dx=$6 dz=$7 grd=$8

set R = -R$Xmin/$Xmax/$Zmin/$Zmax
set I = -I$dx/$dz
blockmean $xyz $R $I -V  > tmp.bmean
$GMTsurface tmp.bmean -G$grd $I $R -V 
