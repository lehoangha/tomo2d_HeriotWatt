#!/bin/csh -f
# make_data.sh
set path = (~/src/tomo $path)

set X=40 Z=10

set id = type1
set geom = geom.dat.$id
./getbathy.sh
echo "1 20 39" > tmp.xsrc
(cat tmp.xsrc ; cat tmp.xt) | \
awk 'NR==1{ns=NF; split($0,xsrc);i=1}NR>1 && $1==xsrc[i] {printf "%s ", $2; i++} END{printf "\n"}' > tmp.zsrc

cat tmp.xsrc > tmp.in
cat tmp.zsrc >> tmp.in
cat >> tmp.in <<EOF
1 39 2
0
2 39
2 39
EOF

./make_geom.awk tmp.in > $geom

# make observed traveltimes 
prep:
set smesh = smesh.dat
set refl = refl.dat.cos
tt_forward -M$smesh -G$geom -F$refl \
           -V1 -N5/10/1/8/1e-4/1e-5 -Ccputime.fine.$id \
           -Rray.dat.$id -Ssrc.dat.$id > obs_tt.dat.$id



