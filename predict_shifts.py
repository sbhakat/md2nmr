#! /usr/bin/env python

import sys, os
import subprocess as sub
import cPickle as pickle
import MDAnalysis as mda
import ShiftPred as sp
import md2nmr

simPath         = '../..'
runName         = '04_prod01_protein'
expDatFile      = simPath + '/' + 'common_files/IL6_S2_exp.dat'
simulations     = {}
skipFrames      = 1000
rerunMD2NMR     = False
rerunShiftPred  = True
shiftPredMethod = 'sparta+'
 
#indices     = [6, 12, 13, 14 ,15, 21, 22, 23, 26, 31]
indices     = [12]
forcefields = ['charmm27']
#forcefields = ['amber03',
#               'amber03-star',
#               'amber99sb',
#               'amber99sb-star',
#               'amber99sb-ildn',
#               'amber99sb-star-ildn',
#               'amber99sbnmr1-ildn',
#               'charmm22star',
#               'charmm27',
#               'oplsaa',
#               'gromos54a7']


print "Creating simulation objects"

simulations = {}
for ff in forcefields:
    for index in indices: 

        simName = "{:s}_2IL6_{:d}".format(ff, index)
        print simName
        simulations[simName] = md2nmr.md2nmr(runName, path=simPath+'/'+simName, rerun=rerunMD2NMR, verbose=False)
 

print ""

print "Predicting chemical shifts"
shiftPredictions = {}
for key in simulations.keys():

    pickleFilename = "{}/{}/{}_shiftpred_{}.dat".format(simPath, simName, runName, shiftPredMethod)

    if not rerunShiftPred and os.path.isfile(pickleFilename):
        print "unpickling average shifts"
        # unpickle average shifts
        loadFile = open(pickleFilename, 'rb')
        shiftPredictions[key] = pickle.load(loadFile)
        loadFile.close()  

    else:
        shiftPredictions[key] = sp.ShiftPred(simulations[key].universe, method=shiftPredMethod)
        shiftPredictions[key].predict(skip=skipFrames)
        shiftPredictions[key].universe = None # unset universe, so we can pickle

        # pickle average shifts
        print "pickling average shifts"
        dumpFile = open(pickleFilename, 'wb')
        pickle.dump(shiftPredictions[key], dumpFile, protocol=pickle.HIGHEST_PROTOCOL)
        dumpFile.close()  



