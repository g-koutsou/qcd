for z in $(seq 0 9)
do
    for y in $(seq 0 9)
    do
	for x in $(seq 0 9)
	do
	    cat <<EOF > rms.ini
<processors_txyz>1 1 1 1</processors_txyz>
<lattice_txyz>20 10 10 10</lattice_txyz>
<source_pos_txyz>0 $x $y $z</source_pos_txyz>
<alpha_gauss> 5 1 5 </alpha_gauss>
<nsmear_gauss>1</nsmear_gauss>
<alpha_APE>0.5</alpha_APE>
<nsmear_APE>5</nsmear_APE>
<cfg_name>/home/koutsou/data/Confs/Nf0_b5p72_L10T20/conf_Nf0_b5p72_L10T20.0080000.scidac</cfg_name>
EOF
	    echo $x $y $z $(mpiexec -np 1 ./rms.exe rms.ini | grep 'Nsmear =   0,')
	done
    done
done