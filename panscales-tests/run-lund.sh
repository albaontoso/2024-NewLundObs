nevents=(1)

for nevent in "${nevents[@]}" 
do
    sed -e "s/&nev&/${nevent}/g" run-lund-lxplus.sub > run-lund-lxplus.tmp
    condor_submit run-lund-lxplus.tmp
    rm run-lund-lxplus.tmp
done
