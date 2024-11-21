etamins=(5 6 7 8)
njob=10
for etamin in "${etamins[@]}" 
do
    etamax=$(($etamin+1))
    sed -e "s/&ymin&/${etamin}/g" -e "s/&ymax&/${etamax}/g" -e "s/&njobs&/$njob/g"  run-lund-lxplus.sub > run-lund-lxplus.tmp
    condor_submit run-lund-lxplus.tmp
    rm run-lund-lxplus.tmp
done
