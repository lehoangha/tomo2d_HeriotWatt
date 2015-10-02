#!/bin/csh -f
# make_fig.sh
set path = (~/src/tomo $path)

set ps = invex.ps
set cpt = cpt.gray
cat > $cpt <<EOF
-10 255 255 255 -5 200 200 200
-5 200 200 200 0 170 170 170
0 170 170 170 5 100 100 100
5 100 100 100 10 0 0 0
EOF

set X=40 Z=10

set Wdot = -W5/0/0/0to Msrc = "-Sc0.05 -W3 -G255"

set ex1 = out.w1.type1/out.type1.smesh.2.1
set ex2 = out.w100.type1/out.type1.smesh.2.1
set ex1r = out.w1.type1/out.type1.refl.2.1
set ex2r = out.w100.type1/out.type1.refl.2.1

#goto skip
tt_forward -Msmesh.dat -Itmp.v0 > /dev/null
./togrid.sh tmp.v0 0 $X -1.0 $Z 1 0.25 tmp.grd0
tt_forward -M$ex1 -Itmp.v1 > /dev/null
./togrid.sh tmp.v1 0 $X -1.0 $Z 1 0.25 tmp.grd1
tt_forward -M$ex2 -Itmp.v2 > /dev/null
./togrid.sh tmp.v2 0 $X -1.0 $Z 1 0.25 tmp.grd2

grdmath tmp.grd1 tmp.grd0 - = tmp.diff
grdmath tmp.diff tmp.grd0 / = tmp.diff2
grdmath tmp.diff2 100 x = tmp.par1
grdmath tmp.grd2 tmp.grd0 - = tmp.diff
grdmath tmp.diff tmp.grd0 / = tmp.diff2
grdmath tmp.diff2 100 x = tmp.par2

skip:

psbasemap -JX2.5/-1.4 -R0/$X/-1/$Z -Ba5:"x [km]":/a2:"z [km]":WeSn -X1.5 -Y7 -K -P > $ps
(echo "0 0\n $X 0") | psxy -JX -R -W2 -O -K >> $ps
(cat bathy.dat; echo "40 10"; echo "0 10") \
| psxy -JX -R -W5 -L -G170 -O -K >> $ps
(echo "0 10"; cat refl.dat.cos; echo "40 10") \
| psxy -JX -R -W5 -L -G255 -O -K >> $ps
psxy refl.dat -JX -R $Wdot -O -K >> $ps
pstext -JX -R -O -K -N >> $ps <<EOF
0 -2.3 14 0 4 5 (a) True and Starting Models
EOF

psscale -C$cpt -D2.9/-3.5/3.5/0.2h -B:"Velocity Perturbation [%]": -O -K >> $ps

psbasemap -JX -R -X3.3 -Ba5:"x [km]":/a2:"z [km]":WeSn -O -K >> $ps
(echo "0 0\n $X 0") | psxy -JX -R -W2 -O -K >> $ps
psxy bathy.dat -JX -R -W5 -O -K >> $ps
psxy refl.dat.cos -JX -R -W5 -O -K >> $ps
psxy ray.dat.type1 -JX -R -W1/0 -M -O -K >> $ps
psxy src.dat.type1 -JX -R $Msrc -O -K >> $ps

pstext -JX -R -O -K -N >> $ps <<EOF
0 -2.3 14 0 4 5 (b) Ray Paths
EOF

psbasemap -JX -R -X-3.3 -Y-2.6 -B -O -K >> $ps
grdimage tmp.par1 -JX -R -C$cpt -O -K -B >> $ps
grdcontour tmp.par1 -JX -R -C1 -L0.1/10 -O -K >> $ps
grdcontour tmp.par1 -JX -R -C1 -L-10/-0.1 -O -K >> $ps
(cat bathy.dat; echo "40 -1"; echo "0 -1") \
| psxy -JX -R -W5 -L -G255 -O -K >> $ps
(echo "0 0\n $X 0") | psxy -JX -R -W2 -O -K >> $ps
(echo "0 10"; cat $ex1r; echo "40 10") \
| psxy -JX -R -W5 -L -G255 -O -K >> $ps
psxy refl.dat -JX -R $Wdot -O -K >> $ps

pstext -JX -R -O -K -N >> $ps <<EOF
0 -2.3 14 0 4 5 (c) Recovery (w=1)
EOF

psbasemap -JX -R -X3.3 -B -O -K >> $ps
grdimage tmp.par2 -JX -R -C$cpt -O -K -B >> $ps
(cat bathy.dat; echo "40 -1"; echo "0 -1") \
| psxy -JX -R -W5 -L -G255 -O -K >> $ps
(echo "0 0\n $X 0") | psxy -JX -R -W2 -O -K >> $ps
(echo "0 10"; cat $ex2r; echo "40 10") \
| psxy -JX -R -W5 -L -G255 -O -K >> $ps
psxy refl.dat -JX -R $Wdot -O -K >> $ps

pstext -JX -R -O -N >> $ps <<EOF
0 -2.3 14 0 4 5 (d) Recovery (w=100)
EOF

