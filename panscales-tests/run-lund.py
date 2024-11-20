#!/usr/bin/env python3
import os 
import sys
import numpy as np

print(sys.version)
print(sys.executable)

base     = "/afs/cern.ch/work/a/asotoont/private/2024-NewLundObs/panscales-tests/"
outdir   = "results-lxplus"
prg      = "/afs/cern.ch/work/a/asotoont/private/2024-NewLundObs/panscales-tests/build-double/lund-analysis"

os.system("mkdir -p "+base+outdir)

print("Base= "+base)
print("Results will be copied in "+base+outdir)

nev        = sys.argv[1]
seed       = sys.argv[2]

print("Input arguments: nev = "+nev+"; seed = "+seed)

basic_settings = f"-shower panglobal -beta 0.0 -physical-coupling -no-spin-corr -match-process -process ee2qq -rts 1000"
extra          = " -nev "+str(nev)
outname        ="pg00"+"_seed_"+str(seed)+".dat"
    
command=prg+" -out "+outname+extra
command+=" "+basic_settings+" -rseq "+seed

print("Command = "+command)
print("Load compiler and run")

os.system(command)
os.system("mv "+outname+" "+base+outdir)
