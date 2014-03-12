#! /usr/bin/env python

import sys, os
from mpi4py import MPI
import subprocess as sub
import cPickle as pickle
import MDAnalysis as mda
import ShiftPred as sp
import md2nmr

simPath         = '/local/jubio/oschill/iGRASP/IL-6/IL-6_ffcomp'
runName         = '04_prod01_protein'
expDatFile      = simPath + '/' + 'common_files/IL6_S2_exp.dat'
simulations     = {}
skipFrames      = 10
rerunMD2NMR     = False
rerunShiftPred  = True
shiftPredMethod = 'sparta+'

comm        = MPI.COMM_WORLD
rank        = comm.Get_rank()
numprocs    = comm.Get_size()
root        = 0 

# root distributes work
if rank == root:
    indices     = [6, 12, 13, 14 ,15, 21, 22, 23, 26, 31]
    forcefields = ['amber03',
                    'amber03-star',
                    'amber99sb',
                    'amber99sb-star',
                    'amber99sb-ildn',
                    'amber99sb-star-ildn',
                    'amber99sbnmr1-ildn',
                    'charmm22star',
                    'charmm27',
                    'oplsaa',
                    'gromos54a7']
   #indices     = [15, 21]
   #forcefields = ['charmm27']


    # create work items
    workItems = []
    for ff in forcefields:
        for index in indices:
            simName = "{:s}_2IL6_{:d}".format(ff, index)
            workItems.append(simName)

    # compute workload per process
    totalWorkload = len(workItems)
    workload  = totalWorkload / numprocs
    extraWork = totalWorkload % numprocs

    # loop through processes
    for proc in range(numprocs):
        currentWorkload = workload
        if proc < extraWork:
            currentWorkload += 1

        # pack work items per process
        work = []
        for w in range(currentWorkload):
            work.append(workItems.pop())

        # send work items
        if proc != 0:
            comm.send(work, dest=proc, tag=123)
        else:
            mywork = work


# receive work from root
if rank != root:
    work = comm.recv(source=0, tag=123)
else:
    work = mywork


# report work
for proc in range(numprocs):
    if rank == proc:
        for w in work:
            print rank, w
        #print ""
        sys.stdout.flush()
    comm.Barrier()


# Creating simulation objects
for simName in work:
    simulations[simName] = md2nmr.md2nmr(runName, path=simPath+'/'+simName, rerun=rerunMD2NMR, verbose=False)


# Predicting chemical shifts
shiftPredictions = {}
for i, key in enumerate(simulations.keys()):

    simName = key
    pickleFilename = "{}/{}/{}_shiftpred_{}_skip{}.dat".format(simPath, simName, runName, shiftPredMethod, skipFrames)
    #print pickleFilename

    if not rerunShiftPred and os.path.isfile(pickleFilename):
        #print "unpickling average shifts"
        # unpickle average shifts
        loadFile = open(pickleFilename, 'rb')
        shiftPredictions[key] = pickle.load(loadFile)
        loadFile.close()  

    else:
        shiftPredictions[key] = sp.ShiftPred(simulations[key].universe, method=shiftPredMethod)
        shiftPredictions[key].predict(skip=skipFrames)
        shiftPredictions[key].universe = None # unset universe, so we can pickle

        # pickle average shifts
        #print "pickling average shifts"
        dumpFile = open(pickleFilename, 'wb')
        pickle.dump(shiftPredictions[key], dumpFile, protocol=pickle.HIGHEST_PROTOCOL)
        dumpFile.close()  

    print "rank {}, {}/{} done".format(rank, i+1, len(simulations.keys()))
 

print "Process {} done.".format(rank) 












