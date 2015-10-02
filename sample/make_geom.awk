
# make_geom.awk
# inputs: line 1: sxs
#         line 2: szs
#         line 3: rx_begin rx_end drx
#         line 4: rz
#         line 5: range_min range_max (for refraction)
#         line 6: range_min range_max (for reflection)

awk '
NR == 1 { 
	nsrc = NF;
	for (i=1; i<=nsrc; i++) sx[i] = $i;
}
NR == 2 { 
	if (NF != nsrc){ print "error in line 2" > "/dev/tty"; }
	for (i=1; i<=nsrc; i++) sz[i] = $i;
}
NR == 3 {
	i=0;
	for (x=$1; x<=$2; x+=$3) rx[++i] = x;
	nrcv = i;
}
NR == 4 {
	for (i=1; i<=nrcv; i++) rz[i] = $1;
}
NR == 5 { rmin0 = $1; rmax0 = $2; }
NR == 6 { rmin1 = $1; rmax1 = $2; }
END {
	print nsrc;
	for (i=1; i<=nsrc; i++){ 
                nr0 = 0; nr1 = 0; 
                for (j=1; j<=nrcv; j++){ 
			if (abs(rx[j]-sx[i])>=rmin0 && 
			    abs(rx[j]-sx[i])<=rmax0) nr0++;
			if (abs(rx[j]-sx[i])>=rmin1 && 
			    abs(rx[j]-sx[i])<=rmax1) nr1++;
		}
		print "s", sx[i], sz[i], nr0+nr1;
		for (j=1; j<=nrcv; j++){ 
			if (abs(rx[j]-sx[i])>=rmin0 && 
			    abs(rx[j]-sx[i])<=rmax0){
				print "r", rx[j], rz[j], 0, 0, 0;\
			}
		}
		for (j=1; j<=nrcv; j++){ 
			if (abs(rx[j]-sx[i])>=rmin1 && 
			    abs(rx[j]-sx[i])<=rmax1){
				print "r", rx[j], rz[j], 1, 0, 0;\
			}
		}
	}
}

function abs(x) {
	return x >= 0 ? x : -x;
}
' $*
