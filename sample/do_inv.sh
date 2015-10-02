#!/bin/csh -f
# do_inv.sh
set path = (~/src/tomo $path)

if ($#argv != 2) then
    echo "needs two arg";
    exit;
endif

set id = $1
set weight = $2
set X=40 Z=10
set data = obs_tt.dat.$id
set init_mesh = smesh.dat
set init_refl = refl.dat
set vcorr = vcorr.dat
set dcorr = dcorr.dat
set N = -N5/10/1/8/1e-4/1e-5

set Lht=2 Lhb=5 Lvt=0.5 Lvb=1 LhR=5
cat > $vcorr <<EOF
2 2 
0.0 40.0
0.0 0.0
0.0 10.0
$Lht $Lhb
$Lht $Lhb
$Lvt $Lvb
$Lvt $Lvb
EOF

cat > $dcorr <<EOF
0.0 $LhR
40  $LhR
EOF

set dir = out.w$2.$1
if (! -e $dir) then
    mkdir $dir
endif

tt_inverse -M$init_mesh -F$init_refl -G$data $N \
           -L$dir/log.all.$id -O$dir/out.$id -V1 -W$weight -Q1e-3 \
           -I2 -SV20 -SD2 -CV$vcorr -CD$dcorr 




