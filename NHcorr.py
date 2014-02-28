#! /usr/bin/env python

import sys
import threading
import multiprocessing
import numpy as np

# ============================================================================ #

class NHcorr: 
    """Compute the NH correlation function from NH vector data and estimate
    convergence and S2 order parameters"""

    def __init__(self, NHvectors, dt, t_lag, t0=0.0):
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

        self.corr      = np.zeros([self.nvectors, self.nf]) # correlation function per NH vector
        self.corr_std  = np.zeros([self.nvectors, self.nf]) # standard deviations of correlation functions
        self.corr_conf = np.zeros([self.nvectors, self.nf]) # standard error of the mean per NH vector

        self.S2        = None    # NMR order parameters per NH vector
        self.S2_std    = None    # standard deviation of each order parameter
        self.S2_conf   = None    # standard error of the mean of each order parameter 

        if self.f0 < 0:
            print "starting time for correlation function must be >= 0"
            sys.exit(1)

        if self.f0 + 2*self.nf > self.nframes:
            print "starting time + 2 * lag time to use for correlation function must be <= total time."
            sys.exit(1)


# ============================================================================ #

    def compute_NH_correlation(self, verbose=False):
        """compute NH correlation function for each NH vector"""

        if verbose:
            sys.stdout.write('Computing NH correlations [  0.0%]') 

        for vector in range(self.nvectors):

            if verbose:
                message = "[{0:6.1%}]".format(1.0*(vector+1)/self.nvectors)
                sys.stdout.write(len(message)*'\b' + message)
                sys.stdout.flush() 

            self.single_NH_correlation(vector)

        
        if verbose:
            sys.stdout.write('\n')

# ============================================================================ #

    def compute_NH_correlation_threaded(self, verbose=False):
        """compute NH correlation function for each NH vector"""

        print "Does not work."
        print "ToDo: Retrieve result from Child Processes"
        print "Maybe through Queue or Pool"
        return

        nCores = multiprocessing.cpu_count()
        vector = 0

        if verbose:
            sys.stdout.write('Computing NH correlations [  0.0%]') 
        
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

            if verbose:
                message = "[{0:6.1%}]".format(1.0*(vector)/self.nvectors)
                sys.stdout.write(len(message)*'\b' + message)
                sys.stdout.flush()  


        if verbose:
            sys.stdout.write('\n') 
        

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
            self.compute_NH_correlation_threaded()

        NMRanalysisBondLengths = 1.02
        NMRanalysisBondLengthsProper = 1.04
        #scalingFactor = (NMRanalysisBondLengths / self.NHbondlengths.mean())**6
        scalingFactor = (NMRanalysisBondLengths / NMRanalysisBondLengthsProper)**6

        self.S2      = scalingFactor * self.corr.mean(1)
        self.S2_std  = self.corr.std(1)
        self.S2_conf = self.S2_std / self.corr.shape[1]**0.5 

# ============================================================================ #


#def single_NH_correlation(NHvectors, f0, nf, vector):
#    """compute NH correlation function of a single NH vector
#    vector is the vector index"""
#
#    corr      = np.zeros([nf]) # correlation function
#    corr_std  = np.zeros([nf]) # standard deviations of correlation function
#    corr_conf = np.zeros([nf]) # standard error of the mean 
#        
#
#    # second order legendre polynomial of NH vector dotproducts P2[i,j] <=> P2[t,t+tau]
#    P2 = np.polynomial.legendre.legval( np.dot(NHvectors[f0:f0+2*nf,vector,:],
#                                               NHvectors[f0:f0+2*nf,vector,:].T),
#                                               [0,0,1])
#
#    # compute the correlation function for each lag time
#    for frame in range(nf):
#        d = np.diagonal(P2, frame)
#        corr     [frame] = d.mean()
#        corr_std [frame] = d.std()
#        corr_conf[frame] = corr_std[frame] / d.shape[0]**0.5  


