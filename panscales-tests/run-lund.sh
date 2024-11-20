nevents=(1)
nruns=(10)
for nevent in "${nevents[@]}" 
do
    nrun="${nruns[iruns]}"
    sed -e "s/&nev&/${nevent}/g" -e "s/&njobs&/$nrun/g" run-lund-lxplus.sub > run-lund-lxplus.tmp
    condor_submit run-lund-lxplus.tmp
    rm run-lund-lxplus.tmp
    iruns=$iruns+1
done
