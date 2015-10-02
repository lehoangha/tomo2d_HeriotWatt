#!/bin/csh -f 
# prep_model.sh - 
set path = (~/src/tomo $path)

set X=40 Z=10 Zsf=2.5
set vwater = 1.5

echo $X $Zsf | awk '{srand(1); for(x=0;x<=$1+0.01;x+=0.1){print x,$2+5*(rand()-0.5);}}' \
| filter1d -Fg5 -E | sample1d -I2 > bathy.dat
#cp ../invex/bathy.dat bathy.dat

# make a Zelt model first
cat > zelt.dat <<EOF
1     0.0    $X
0     0.0   0.0
        0     0
1      $X
0     1.5
        0
1      $X
0     1.5
        0
EOF

awk 'BEGIN{printf "2 ";}{printf "%.1f ", $1}END{printf "\n";}' bathy.dat >> zelt.dat
awk 'BEGIN{printf "0 ";}{printf "%.2f ", $2}END{printf "\n";}' bathy.dat >> zelt.dat
awk 'BEGIN{printf "  "; }{printf "0 "; }END{printf "\n";}' bathy.dat >> zelt.dat

cat >> zelt.dat <<EOF
2     $X
0     3.0
        0
2     $X
0     5.5
        0
3     0.0    $X
0     3.0   3.0
        0     0
3     $X
0     5.5
        0
3     $X
0     6.5
        0
4     0.0    $X
0     3.5    3.5
        0      0
4     $X
0     6.5
        0
4     $X
0     7.0
        0
5     0.0  $X
0     9.0  9.0
        0    0
5     $X
0     7.0
        0
5     $X
0     7.0
        0
6     0       $X
0     15.0  15.0
EOF

# make a rough velocity mesh
set dX=1.0 dZ=0.05
echo "$dZ" | awk '{for (z=0; z<=8; ){ dz=$1+sqrt(0.02*z); dz*=1; print z; z+=dz}}' \
> tmp.zs
gen_smesh -Czelt.dat/2 -F5/refl.dat \
          -Q$vwater -R$vwater -E$dX -Ztmp.zs > smesh.dat
awk '{print $1, $2+0.5*sin(2*3.14*$1/40)}' refl.dat > refl.dat.cos






