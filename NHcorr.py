#! /usr/bin/env python

import sys, os
import cPickle as pickle
import threading
import multiprocessing
import numpy as np

# ============================================================================ #

class NHcorr: 
    """Compute the NH correlation function from NH vector data and estimate
    convergence and S2 order parameters"""

    def __init__(self, NHvectors, name, path, dt, t_lag, t0=0.0, S2expFile=None, rerun=False, verbose=False):
        """NHvectors is the time series of NH bond vectors as NHvectors[nframes, nvectors, 3]
        dt is the timestep
        t0 is the starting time
        t_lag is the maximum lag time"""

        self.NHvectors = NHvectors
        self.nframes   = NHvectors.shape[0]
        self.nvectors  = NHvectors.shape[1]

        self.dt        = dt      # timestep
        self.t0        = t0      # starting time
        self.t_lag     = t_lag   # maximum lag time
        self.f0        = int(t0 / dt)          # starting frame number
        self.nf        = int((t_lag-t0)/dt)  # number of frames for correlation function

        self.name       = name  # trajectory name
        self.path       = path  # trajectory path
        self.NHcorrFile = self.path + '/' + self.name + '_NHcorr_f0-{0}_nf-{1}.dat'.format(self.f0, self.nf)
        self.S2File     = self.path + '/' + self.name + '_S2.dat'
        self.S2expFile  = S2expFile # experimental S2 values
        self.rerun      = rerun     # do not load intermediate results from pickled files
        self.verbose    = verbose   # if True, write verbose output to sys.stdout
 
        self.corr      = np.zeros([self.nvectors, self.nf]) # correlation function per NH vector
        self.corr_std  = np.zeros([self.nvectors, self.nf]) # standard deviations of correlation functions
        self.corr_conf = np.zeros([self.nvectors, self.nf]) # standard error of the mean per NH vector

        self.S2         = None    # NMR order parameters per NH vector
        self.S2_std     = None    # standard deviation of each order parameter
        self.S2_conf    = None    # standard error of the mean of each order parameter 
        self.S2_exp     = None    # experimental S2 values
        self.S2_exp_err = None    # errors of experimental S2 values
        self.S2_diff    = None    # difference of estimated and experimental S2 values

        if self.f0 < 0:
            sys.stderr.write("starting time for correlation function must be >= 0\n")
            sys.exit(1)

        if self.f0 + 2*self.nf > self.nframes:
            sys.stderr.write("starting time + 2 * lag time to use for correlation function must be <= total time.\n")
            sys.exit(1)


# ============================================================================ #

    def compute_NH_correlation(self):
        """compute NH correlation function for each NH vector"""

        if os.path.isfile(self.NHcorrFile) and not self.rerun:
            self.vout("Loading NH correlations from file.\n")

            # unpickle data
            loadFile = open(self.NHcorrFile, 'rb')
            (self.corr, self.corr_std, self.corr_conf) = pickle.load(loadFile)
            loadFile.close()  

        else:

            self.vout('Computing NH correlations [  0.0%]') 

            for vector in range(self.nvectors):

                message = "[{0:6.1%}]".format(1.0*(vector+1)/self.nvectors)
                self.vout(len(message)*'\b' + message)
                sys.stdout.flush() 

                self.single_NH_correlation(vector)


            self.vout('\n')

            # pickle trajectory NH correlations
            dumpFile = open(self.NHcorrFile, 'wb')
            pickle.dump((self.corr, self.corr_std, self.corr_conf), dumpFile, protocol=pickle.HIGHEST_PROTOCOL)
            dumpFile.close()   


# ============================================================================ #

    def compute_NH_correlation_threaded(self):
        """compute NH correlation function for each NH vector"""

        sys.stderr.write("compute_NH_correlation_threaded does not work.\n")
        sys.stderr.write("ToDo: Retrieve result from Child Processes.\n")
        sys.stderr.write("Maybe through Queue or Pool.\n")
        sys.exit(1)

        nCores = multiprocessing.cpu_count()
        vector = 0

        self.vout('Computing NH correlations [  0.0%]') 
        
        while vector < self.nvectors:

            # start as many threads as cores
            threadID = 0
            threads = []
            while threadID < nCores and vector < self.nvectors:
                #threads.append( multiprocessing.Process(target=single_NH_correlation, args=(self.NHvectors, self.f0, self.nf, vector)) )
                threads.append( multiprocessing.Process(target=self.single_NH_correlation, args=(vector, )) )
                threads[-1].start()
                threadID += 1
                vector += 1

            # wait for threads to finish
            for i in range(threadID):
                threads[-1].join()

            message = "[{0:6.1%}]".format(1.0*(vector)/self.nvectors)
            self.vout(len(message)*'\b' + message)
            sys.stdout.flush()  


        self.vout('\n') 
        

# ============================================================================ #

    def single_NH_correlation(self, vector):
        """compute NH correlation function of a single NH vector
        vector is the vector index"""
            

        # second order legendre polynomial of NH vector dotproducts P2[i,j] <=> P2[t,t+tau]
        P2 = np.polynomial.legendre.legval( np.dot(self.NHvectors[self.f0:self.f0+2*self.nf,vector,:],
                                                   self.NHvectors[self.f0:self.f0+2*self.nf,vector,:].T),
                                                   [0,0,1])

        # compute the correlation function for each lag time
        for frame in range(self.nf):
            d = np.diagonal(P2, frame)
            self.corr     [vector, frame] = d.mean()
            self.corr_std [vector, frame] = d.std()
            self.corr_conf[vector, frame] = self.corr_std[vector, frame] / d.shape[0]**0.5 
 
# ============================================================================ #

    def compute_S2(self):
        """Compute S2 order parameters for each NH bond"""

        if self.corr == None:
            self.compute_NH_correlation()

        self.vout("Computing S2 order parameters\n")

        NMRanalysisBondLength = 1.02
        NMRanalysisBondLengthProper = 1.04
        #scalingFactor = (NMRanalysisBondLength / self.NHbondlength.mean())**6
        scalingFactor = (NMRanalysisBondLength / NMRanalysisBondLengthProper)**6

        self.S2      = scalingFactor * self.corr.mean(1)
        self.S2_std  = self.corr.std(1)
        self.S2_conf = self.S2_std / self.corr.shape[1]**0.5 

# ============================================================================ #

    def compare_S2_exp(self, expDatFile, resids, resnames):
        """Compare S2 order parameter estimates to experimental data
        stored in file expDatFile"""

        self.read_exp_S2(expDatFile, resids, resnames)

        self.S2_diff = self.S2_exp - self.S2

# ============================================================================ #

    def read_exp_S2(self, expDatFile, resids, resnames):
        """read experimental S2 order parameters from file"""

        if not os.path.isfile(expDatFile):
            sys.stderr.write("File {0} with experimental S2 parameters does not exist.\n")
            sys.exit(1)

        expDatFp = open(expDatFile, 'r')
        expDat   = expDatFp.readlines()
        expDatFp.close()

        self.S2_exp     = np.zeros_like(self.S2)
        self.S2_exp_err = np.zeros_like(self.S2)

        pos = 0
        for line in expDat[1:]:
            fields = line.split()

            resname   = fields[0]
            resid     = int(fields[1])

            # if residues match
            if pos < len(resids) and resids[pos] == resid and resnames[pos] == resname:
                if len(fields) == 4:
                    self.S2_exp[pos]     = float(fields[2])
                    self.S2_exp_err[pos] = float(fields[3])
                else:
                    self.S2_exp[pos] = np.nan
                pos += 1

# ============================================================================ #

    def vout(self, message):
        """verbose output
        write to sys.stdout if verbose output has been requested"""

        if self.verbose:
            sys.stdout.write(message)
 
# ============================================================================ #
