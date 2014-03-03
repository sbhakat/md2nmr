#! /usr/bin/env python

import sys, os
import cPickle as pickle
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.align as align
import NHcorr
from md2nmr import *

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
        self.trjmetafile     = ""           # Filename of trajectory information that is not contained in fitted trajectory
        self.NHvectorFile    = ""           # Filename of pickled NH vectors
        self.rerun           = rerun        # if True, do the whole analysis again and don't load anything precomputed results from disk

        self.NselectionText  = 'name N and not resname PRO and not resid 20'
        self.HselectionText  = 'name H and not resname PRO and not resid 20' 

        self.universe       = None           # create MDAnalysis universe
        self.NHvectors      = None           # NH vectors per timeframe
        self.NHbondlengths  = None           # NH bondlengths per timeframe
        self.N              = None           # N atom coordinates per timeframe
        self.H              = None           # H atom coordinates per timeframe
        self.time           = None           # time in ps 
        self.trjMetaData    = TrjMetaData()  # Meta data for trajectory that is not stored in aligned trajectory
        self.resids         = None           # residue IDs
        self.resnames       = None           # redidue names
        self.NHcorrelations = []             # list of NH correlation functions and S2 order parameters for different starting and lag times

        if self.topology == None:
            self.topology = self.name + '.pdb'

        if self.trajectory == None:
            self.trajectory = self.name + '.xtc' 

        self.name         = "_".join(self.name.split()) # substitute underscores for whitespace
        self.topology     = self.path + '/' + self.topology
        self.trajectory   = self.path + '/' + self.trajectory
        self.alignedTrj   = self.path + '/' + self.name + '_aligned.xtc'   # aligned trajectory
        self.trjmetafile  = self.path + '/' + self.name + '_metadata.dat'  # metadata not contained in aligned trajectory
        self.NHvectorFile = self.path + '/' + self.name + '_NHvectors.dat' # pickled NHvectors

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

        print "Creating universe for {0}".format(self.name)
        if os.path.isfile(self.alignedTrj) and not self.rerun:
            self.universe = mda.Universe(self.topology, self.alignedTrj)
            self.trajectoryAligned = True

            # unpickle trajectory meta data
            loadFile = open(self.trjmetafile, 'rb')
            self.trjMetaData = pickle.load(loadFile)
            loadFile.close() 

        else:
            self.universe = mda.Universe(self.topology, self.trajectory)
            self.trajectoryAligned = False
            self.trjMetaData.dt        = self.universe.trajectory.dt
            self.trjMetaData.totaltime = self.universe.trajectory.totaltime
            self.trjMetaData.numframes = self.universe.trajectory.numframes
            self.trjMetaData.units     = self.universe.trajectory.units

# ============================================================================ #

    def compute_order_parameters(self):
        """Compute correlation functions and order parameters"""

        self.rms_fit()
        self.compute_NH_vectors()
        self.analyze_NH(verbose=True)
 
# ============================================================================ #

    def rms_fit(self):
        """Superimpose whole trajectory on first frame"""

        if self.trajectoryAligned and not self.rerun:
            sys.stdout.write('Trajectory of {0} already fitted\n'.format(self.name))
        else:

            sys.stdout.write('Fitting trajectory of {0}\n'.format(self.name))
            self.universe.trajectory.rewind()
            align.rms_fit_trj(self.universe,
                              self.universe,
                              select='backbone',
                              filename=self.alignedTrj)

            # pickle trajectory meta data
            dumpFile = open(self.trjmetafile, 'wb')
            pickle.dump(self.trjMetaData, dumpFile, protocol=pickle.HIGHEST_PROTOCOL)
            dumpFile.close() 

            # load fittet trajectory
            print "Reading fittet trajectory of {0}".format(self.name)
            self.universe = mda.Universe(self.topology, self.alignedTrj) 

            self.trajectoryAligned = True
                                                         

# ============================================================================ #

    def compute_NH_vectors(self):
        """Compute NH bond vectors"""

        self.N = self.universe.selectAtoms(self.NselectionText)

        self.resids        = self.N.resids()
        self.resnames      = self.N.resnames()
        #self.trjMetaData.resids  = self.N.resids()
        #self.trjMetaData.renames = self.N.resnames()
 

        if os.path.isfile(self.NHvectorFile) and not self.rerun:
            # load NH bond vectors etc. from file
            sys.stdout.write('Loading NH vectors for {0}\n'.format(self.name))

            # unpickle data
            loadFile = open(self.NHvectorFile, 'rb')
            (self.NHvectors, self.NHbondlengths) = pickle.load(loadFile)
            loadFile.close()

        else:

            # compute NH bond vectors

            sys.stdout.write('Computing NH vectors for {0}         '.format(self.name))

            self.universe.trajectory.rewind()
            nframes = self.universe.trajectory.numframes

            self.NHvectors     = np.zeros([nframes, self.N.numberOfAtoms(), 3])
            self.NHbondlengths = np.zeros([nframes, self.N.numberOfAtoms()])

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
            dumpFile = open(self.NHvectorFile, 'wb')
            pickle.dump((self.NHvectors, self.NHbondlengths), dumpFile, protocol=pickle.HIGHEST_PROTOCOL)
            dumpFile.close()

# ============================================================================ #
    
    def analyze_NH(self, verbose=False):
        """Compute (or load precomputed) NH correlation functions.
        Then compute S2 order parameters)"""

        self.NHcorrelations = []

        t_lag = self.trjMetaData.totaltime/2
        
        self.analyze_NH_block(t_lag, verbose=verbose)
   
# ============================================================================ #

    def analyze_NH_block(self, t_lag, t0=0.0, verbose=False):
        """Compute (or load precomputed) NH correlation function for given lag and starting times.
        Then compute S2 order parameters"""

        NHc = NHcorr.NHcorr(self.NHvectors, self.name, self.path, self.trjMetaData.dt, t_lag, t0=t0, rerun=self.rerun)
        NHc.compute_NH_correlation(verbose=verbose)
        NHc.compute_S2()

        self.NHcorrelations.append(NHc)

# ============================================================================ # 

class TrjMetaData:
    """Class containing Meta Data of trajectory
    i.e. timestep, totaltime, etc"""

    def __init__(self, dt=0.0, totaltime=0.0, numframes=0, units={}, resids=None, resnames=None):
        self.dt        = dt        # timestep in ps
        self.totaltime = totaltime # total trajectory time in ps
        self.numframes = numframes # total number of frames
        self.units     = units     # dictionary: {'length': 'nm', 'time': 'ps'}

        #self.resids    = resids    # residue IDs
        #self.resnames  = resnames  # residue names
      
# ============================================================================ # 
