#! /usr/bin/env python
# TODO: Complete ranking of chemical shift data, especially boolean selection

import sys, os, time, csv
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import ShiftPred as sp
import md2nmr

# ============================================================================ #

#simPath         = '/local/jubio/oschill/iGRASP/IL-6/IL-6_ffcomp' # jubio
simPath         = '/home/oliver/SiSc/Courses/Thesis/iGRASP/IL-6/IL-6_ffcomp/analysis/MDs' # IdeaPad
runName         = '04_prod01_protein'
expS2File       = simPath + '/' + 'common_files/IL6_S2_exp.dat'
expCSFile       = simPath + '/' + 'common_files/IL6_CS_101801_clean.txt'
skipFrames      = 99
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
        self.expShifts        = {}

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
        sys.stdout.flush()

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
            sys.stdout.flush()

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
        sys.stdout.flush()

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
            sys.stdout.flush()

        sys.stdout.write('\n')

# ==================================== #

    def read_experimental_shifts(self, expFile):
        """
        Read experimentally determinded chemical shifts from expFile
        """

        # read data
        with open(expFile, 'r') as expF:
            reader = csv.reader(expF)
            header = [field.strip() for field in reader.next()]
            data   = [row for row in reader]

        # initialize data structures
        self.expShifts = {}
        self.expShifts[header[0]] = np.zeros(len(data), dtype=np.int64) # resids
        self.expShifts[header[1]] = []                                  # resnames (one letter code)
        for field in header[2:]:
            self.expShifts[field] = np.zeros(len(data), dtype=np.float64)

        # parse shift data
        for i, line in enumerate(data):

            # first resids and resnames
            self.expShifts[header[0]][i] = int(line[0].strip()) # resid
            self.expShifts[header[1]].append(line[1].strip())   # resname

            # then the shift data
            for j, field in enumerate(line[2:]):

                try: 
                    self.expShifts[header[j+2]][i] = float(field)

                except ValueError:
                    if field.isspace() or len(field) == 0:
                        self.expShifts[header[j+2]][i] = np.nan
                    else:
                        raise

        
# ==================================== #

    def rank_S2_predictions(self, forcefields, indices):

        norms = {}
        means = {}
        SDs   = {}

        for simName in self.S2.keys():

            s2 = self.S2[simName]
            diff = s2 - self.S2exp

            norms[simName] = np.linalg.norm(diff[np.invert(np.isnan(diff))])
            means[simName] = diff[np.invert(np.isnan(diff))].mean()
            SDs[simName]   = diff[np.invert(np.isnan(diff))].std()

        data = means
        
        fig = plt.figure()
        i = 0
        for ff in forcefields:
            x = []; y = []
            for index in indices:
                simName = "{:s}_2IL6_{:d}".format(ff, index)
                x.append(i)
                y.append(data[simName])
                i += 1
            plt.plot(x,y, 'o', 'MarkerSize', 2)

        plt.show()

# ==================================== #

    def rank_shift_predictions(self):
        """
        Rank simulations according to best aggreement with experimental chemical shifts
        """

        simName = self.simNames[0]

        # loop through all elements for which there are chemical shift predictions
        for element in self.shiftPredictions[simName].averageShifts.keys():

            # construct boolean array for comparison
            resids    = self.shiftPredictions[simName].averageShifts[element].resids
            expResids = self.expShifts['Resid']

            preBool = np.zeros_like(resids,    dtype=np.bool)
            expBool = np.zeros_like(expResids, dtype=np.bool)

            # for predicted shifts
            for ID in expResids:
                preBool = np.logical_or(resids==ID, preBool)

            # for experimental shifts
            for ID in resids:
                expBool = np.logical_or(expResids==ID, expBool)

            print expResids[expBool]
            print resids[preBool]

            break

# ============================================================================ #
