#! /bin/bash

import sys, os
import md2nmr
import matplotlib.pyplot as plt
import numpy as np

rerun       = False
simPath     = '../..'
runName     = '04_prod01_protein'
expDatFile  = simPath + '/' + 'common_files/IL6_S2_exp.dat'
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
S2 = {}
for simName in simulations.keys():
    print simName
    simulations[simName].compute_order_parameters()
    S2[simName] = simulations[simName].NHcorrelations[0].S2


# load experimental data
resids   = simulations[simName].resids
resnames = simulations[simName].resnames
simulations[simName].NHcorrelations[0].read_exp_S2(expDatFile, resids, resnames)
S2exp    = simulations[simName].NHcorrelations[0].S2_exp



# compute average S2 per force field
averages = {}
for ff in forcefields:
    averages[ff] = np.zeros_like(S2exp)

    nsims = 0
    for index in indices:
        simName = "{:s}_2IL6_{:d}".format(ff, index)
        averages[ff] += S2[simName]
        nsims += 1

    averages[ff] /= nsims
    


def plot_average(average, S2exp):

    plt.figure()
    plt.plot(average, 'b')
    plt.plot(S2exp, 'r')
    plt.show()

