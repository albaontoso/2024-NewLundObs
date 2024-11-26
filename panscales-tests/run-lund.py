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
etamin     = sys.argv[3]
etamax     = sys.argv[4]
shower     = sys.argv[5]
betaps     = sys.argv[6]
print("Input arguments: shower = "+shower+betaps+"; nev = "+nev+"; seed = "+seed+"; etamin = "+etamin+"; etamax = "+etamax)

basic_settings = f"-nloops 0 -no-spin-corr -process ee2qq -rts 10000"
extra          = " -shower "+shower+" -beta "+str(betaps)+" -nev "+str(nev)+" -eta-min "+str(etamin)+" -eta-max "+str(etamax)
outname        =shower+"_beta"+str(betaps)+"_etamin_"+str(etamin)+"_etamax_"+str(etamax)+"_seed_"+str(seed)+".dat"
    
command=prg+" -out "+outname+extra
command+=" "+basic_settings+" -rseq "+seed

print("Command = "+command)
print("Load compiler and run")

os.system(command)
os.system("mv "+outname+" "+base+outdir)
