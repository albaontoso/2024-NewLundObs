etamins=(5 6 7 8 9)
njob=10
showers=("panglobal -beta 0.0" "panglobal -beta 0.5" "panlocal -beta 0.5" "panlocal -beta 0.0")
for etamin in "${etamins[@]}" 
do
    etamax=$(($etamin+1))
    for shower in "${showers[@]}"
    do
	showername=$(echo ${shower%%' '*})
	beta=${shower#*-beta}
        sed -e "s/&ymin&/${etamin}/g" -e "s/&ymax&/${etamax}/g" -e "s/&njobs&/$njob/g" -e "s/&shower&/${showername}/g" -e "s/&betaps&/${beta}/g"  run-lund-lxplus.sub > run-lund-lxplus.tmp
	condor_submit run-lund-lxplus.tmp
        rm run-lund-lxplus.tmp
    done
done
