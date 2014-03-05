#! /bin/bash

import sys, os
import md2nmr

rerun       = False
simPath     = '../..'
runName     = '04_prod01_protein'
simulations = {}
 
indices     = [6, 12, 13, 14 ,15, 21, 22, 23, 26, 31]
forcefields = ['amber03']
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
        simulations[simName] = md2nmr.md2nmr(runName, path=simPath+'/'+simName, rerun=rerun, verbose=False)

print ""

print "Loading Data into simulation objects"
for simName in simulations.keys():
    print simName
    simulations[simName].compute_order_parameters()



