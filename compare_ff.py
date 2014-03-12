#! /usr/bin/env python

import sys, os, time
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import ShiftPred as sp
import md2nmr

# ============================================================================ #

simPath         = '/local/jubio/oschill/iGRASP/IL-6/IL-6_ffcomp'
runName         = '04_prod01_protein'
expS2File       = simPath + '/' + 'common_files/IL6_S2_exp.dat'
skipFrames      = 10
rerunMD2NMR     = False
rerunShiftPred  = False
shiftPredMethod = 'sparta+'
 

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


#indices     = [6, 12]
#forcefields = ['amber03',
#               'amber03-star']

if __name__ is "__main__":

    from compare_ff import *
    main()

# ============================================================================ #

def main():
    cff = Compare_ff(simPath, runName, skipFrames, indices, forcefields,
                     rerunMD2NMR, rerunShiftPred, shiftPredMethod)

    return cff

# ============================================================================ #

class Compare_ff:

    def __init__(self, simPath,
                       runName,
                       skipFrames,
                       indices,
                       forcefields,
                       rerunMD2NMR=False,
                       rerunShiftPred=False,
                       shiftPredMethod='sparta+'):

        self.simPath         = simPath
        self.runName         = runName
        self.skipFrames      = skipFrames
        self.indices         = indices
        self.forcefields     = forcefields
        self.rerunMD2NMR     = rerunMD2NMR
        self.rerunShiftPred  = rerunShiftPred
        self.shiftPredMethod = shiftPredMethod

        self.simNames         = []
        self.simulations      = {}
        self.shiftPredictions = {}

        self.S2               = {}
        self.S2resids         = []
        self.S2resnames       = []
        self.S2exp            = []

        # generate simNames
        for ff in forcefields:
            for index in indices:
                self.simNames.append("{:s}_2IL6_{:d}".format(ff, index))

# ==================================== #

    def load_S2_predictions(self, expS2File):

        self.simulations = {}
        self.S2 = {}

        starttime = time.time()
        messageTemplate = "loading S2 predictions [{:3d}/{:3d}], ETA: {:4.1f} sec."
        message         = messageTemplate.format(0, len(self.simNames), 0)
        sys.stdout.write(message)

        for n, simName in enumerate(self.simNames):

            self.simulations[simName] = md2nmr.md2nmr(self.runName,
                                        path=self.simPath+'/'+simName,
                                        rerun=self.rerunMD2NMR, verbose=False)

            self.simulations[simName].compute_order_parameters()
            self.S2[simName] = self.simulations[simName].NHcorrelations[0].S2 

            eta = (time.time() - starttime) / (n+1) * (len(self.simNames) - n+1)
            backspace = len(message) * '\b'
            message = messageTemplate.format(n+1, len(self.simNames), eta)
            sys.stdout.write(backspace + message)

        sys.stdout.write('\n')

        # load experimental data
        simName         = self.simulations.keys()[0]
        self.S2resids   = self.simulations[simName].resids
        self.S2resnames = self.simulations[simName].resnames
        self.simulations[simName].NHcorrelations[0].read_exp_S2(expS2File,
                                  self.S2resids, self.S2resnames)
        self.S2exp = self.simulations[simName].NHcorrelations[0].S2_exp 

# ==================================== #

    def load_shift_predictions(self):

        self.shiftPredictions = {}

        starttime = time.time()
        messageTemplate = "loading shift predictions [{:3d}/{:3d}], ETA: {:4.1f} sec."
        message         = messageTemplate.format(0, len(self.simNames), 0)
        sys.stdout.write(message) 

        for n, simName in enumerate(self.simNames):
            pickleFilename = "{}/{}/{}_shiftpred_{}_skip{}.dat".format(self.simPath, 
                              simName, self.runName, self.shiftPredMethod, self.skipFrames)

            try:
                # unpickle shift predictions
                loadFile = open(pickleFilename, 'rb')
                self.shiftPredictions[simName] = pickle.load(loadFile)
                loadFile.close()   
            except IOError:
                print "IOError upon opening {}".format(pickleFilename)
                sys.exit(1)

            eta = (time.time() - starttime) / (n+1) * (len(self.simNames) - n+1)
            backspace = len(message) * '\b'
            message = messageTemplate.format(n+1, len(self.simNames), eta)
            sys.stdout.write(backspace + message) 

        sys.stdout.write('\n')

# ==================================== #

    def rank_S2_predictions(self):

        norms = {}
        means = {}
        SDs   = {}

        for simName in self.S2.keys():

            s2 = self.S2[simName]
            diff = s2 - self.S2exp

            norms[simName] = np.linalg.norm(diff[np.invert(np.isnan(diff))])
            means[simName] = diff[np.invert(np.isnan(diff))].mean()
            SDs[simName]   = diff[np.invert(np.isnan(diff))].std()

        return (norms, means, SDs)

# ==================================== #

    def rank_shift_predictions(self):
        pass

# ============================================================================ #
