#! /usr/bin/env python

import sys, os
import cPickle as pickle
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.align as align

# ============================================================================ #

class md2nmr:
    """md2nmr objects read trajectory and toplogy data and compute NMR observables
    (S2 order parameters and chemical shifts)
    The results can be stored and retrieved"""

    def __init__(self, name, topology=None, trajectory=None, path='.', rerun=False):

        self.name            = name         # name of the trajectory (will be used as basename of output files)
        self.topology        = topology     # pdb or gro file
        self.trajectory      = trajectory   # xtc or trr file 
        self.path            = path         # path of the trajectory data (analysis results will be written here)
        self.rerun           = rerun        # if True, do the whole analysis again and don't load anything precomputed results from disk

        self.NselectionText  = 'name N and not resname PRO and not resid 20'
        self.HselectionText  = 'name H and not resname PRO and not resid 20' 

        self.universe      = None           # create MDAnalysis universe
        self.NHvectors     = None           # NH vectors per timeframe
        self.NHbondlengths = None           # NH bondlengths per timeframe
        self.N             = None           # N atom coordinates per timeframe
        self.H             = None           # H atom coordinates per timeframe
        self.time          = None           # time in ps 
        self.resids        = None           # residue IDs
        self.resnames      = None           # redidue names

        if self.topology == None:
            self.topology = self.name + '.pdb'

        if self.trajectory == None:
            self.trajectory = self.name + '.xtc' 

        self.name = "_".join(self.name.split()) # substitute underscores for whitespace
        self.topology   = self.path + '/' + self.topology
        self.trajectory = self.path + '/' + self.trajectory
        self.alignedTrj = self.path + '/' + self.name + '_aligned.xtc'  # aligned trajectory

        self.trajectoryAligned = False  # True if trajectory is already aligned

        if not os.path.isfile(self.topology):
            print "Cannot locate file: ", self.topology
            sys.exit(0) 

        if not os.path.isfile(self.trajectory):
            print "Cannot locate file: ", self.trajectory 
            sys.exit(0) 

        self.create_universe()

# ============================================================================ #

    def create_universe(self):
        """Create a universe from the data in the given trajectory and topology"""

        print self.name
        print "Creating universe for {0}".format(self.name)
        if os.path.isfile(self.alignedTrj) and not self.rerun:
            self.universe = mda.Universe(self.topology, self.alignedTrj)
            self.trajectoryAligned = True
        else:
            self.universe = mda.Universe(self.topology, self.trajectory)
            self.trajectoryAligned = False

# ============================================================================ #

    def compute_order_parameters(self):
        """Compute correlation functions and order parameters"""

        self.rms_fit()
        self.compute_NH_vectors()
#        self.compute_NH_correlation() 
#        self.compute_S2()
 
# ============================================================================ #

    def rms_fit(self):
        """Superimpose whole trajectory on first frame"""

        if self.trajectoryAligned and not self.rerun:
            sys.stdout.write('Trajectory of {0} already fitted\n'.format(self.name))
        else:

            sys.stdout.write('Fitting trajectory of {0}\n'.format(self.name))
            self.universe.trajectory.rewind()
            self.alignedTrj = align.rms_fit_trj(self.universe,
                                             self.universe,
                                             select='backbone',
                                             filename=self.alignedTrj)

            # load fittet trajectories
            print "Reading fittet trajectory of {0}".format(self.name)
            self.universe = mda.Universe(self.topology, self.alignedTrj) 

            self.trajectoryAligned = True
                                                         

# ============================================================================ #

    def compute_NH_vectors(self):
        """Compute NH bond vectors"""

        NHvectorFile = self.path + '/' + self.name + '_NHvectors.dat'

        if os.path.isfile(NHvectorFile) and not self.rerun:
            # load NH bond vectors etc. from file
            sys.stdout.write('Loading NH vectors for {0}\n'.format(self.name))

            # unpickle data
            loadFile = open(NHvectorFile, 'rb')
            (self.NHvectors, self.NHbondlengths) = pickle.load(loadFile)
            loadFile.close()

        else:

            # compute NH bond vectors

            sys.stdout.write('Computing NH vectors for {0}         '.format(self.name))

            self.universe.trajectory.rewind()
            nframes = self.universe.trajectory.numframes

            self.N = self.universe.selectAtoms(self.NselectionText)
            self.NHvectors     = np.zeros([nframes, self.N.numberOfAtoms(), 3])
            self.NHbondlengths = np.zeros([nframes, self.N.numberOfAtoms()])

            self.resids        = self.N.resids()
            self.resnames      = self.N.resnames()

            # loop through timesteps
            for ts in self.universe.trajectory:

                message = "[{0:6.1%}]".format(1.0*ts.frame/nframes)
                sys.stdout.write(len(message)*'\b' + message)
                sys.stdout.flush() 

                self.N = self.universe.selectAtoms(self.NselectionText)
                self.H = self.universe.selectAtoms(self.HselectionText)

                self.NHvectors[ts.frame-1,:,:] = self.H.coordinates() - self.N.coordinates()

                for atomIndex in range(self.NHvectors.shape[1]):
                    norm = np.linalg.norm(self.NHvectors[ts.frame-1, atomIndex])
                    self.NHvectors[ts.frame-1, atomIndex]    /= norm
                    self.NHbondlengths[ts.frame-1, atomIndex] = norm

            sys.stdout.write('\n')

            # pickle data
            dumpFile = open(NHvectorFile, 'wb')
            pickle.dump((self.NHvectors, self.NHbondlengths), dumpFile, protocol=pickle.HIGHEST_PROTOCOL)
            dumpFile.close()

# ============================================================================ #
    
    def compute_NH_correlation(self, f=0.5, start=0, stop=None):
        """ Compute the angular correlation function of NH bond vectors
        f is the fraction of the data series that is used as maximum lag time"""

        sys.stdout.write('Computing NH correlations for {0}         '.format(self.name))

        if not f >= 0 and f <= 1:
            print "ERROR: maximum lag time as a fraction of data series length must be between 0 and 1"
            sys.exit(1)

        Nframes  = self.NHvectors.shape[0]
        nvectors = self.NHvectors.shape[1]

        if stop == None:
            nframes  = int(round(self.NHvectors.shape[0] * f))
            stop     = nframes
        else:
            nframes = int(round((stop - start) * f))

        self.corr      = np.zeros([nvectors, nframes]) # correlation functions
        self.corr_std  = np.zeros([nvectors, nframes]) # standard deviation
        self.corr_conf = np.zeros([nvectors, nframes]) # standard error of the mean

        for vector in range(nvectors):

            message = "[{0:6.1%}]".format(1.0*(vector+1)/nvectors)
            sys.stdout.write(len(message)*'\b' + message)
            sys.stdout.flush()

            # second order legendre polynomial of NH vector dotproducts P2[i,j] <=> P2[t,t+tau]
            P2 = np.polynomial.legendre.legval( np.dot(self.NHvectors[start:2*stop,vector,:],
                                                       self.NHvectors[start:2*stop,vector,:].T),
                                                       [0,0,1])

            # compute the correlation function for each lag time
            for frame in range(nframes):
                d = np.diagonal(P2, frame)
                self.corr     [vector, frame] = d.mean()
                self.corr_std [vector, frame] = d.std()
                self.corr_conf[vector, frame] = self.corr_std[vector, frame] / d.shape[0]**0.5

        sys.stdout.write('\n')

# ============================================================================ # 
