#! /bin/bash

import sys, os
from mpi4py import MPI
import md2nmr

rerun       = True
simPath     = '../..'
runName     = '04_prod01_protein'
simulations = {}
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


# do work
for simName in work:
    simulations[simName] = md2nmr.md2nmr(runName, path=simPath+'/'+simName, rerun=rerun, verbose=False)
    simulations[simName].compute_order_parameters()


print "Process {} done.".format(rank)
