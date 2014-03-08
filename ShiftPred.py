#! /usr/bin/env python
#TODO:
# write_pdb() write only 'protein'
# sparta()

import sys, os
import string, random
import subprocess as sub
import cPickle as pickle
import numpy as np
import MDAnalysis as mda

# ============================================================================ #

class ShiftPred:  
    """Chemical shift predictor for MD trajectories with the help of 3rd party programs"""

    def __init__(self, universe, method="shiftx2"):

        self.universe = universe  # MDAnalysis trajectory
        self.method   = method    # 3rd party program to use for prediction
        
        self.PDBpath     = ""           # path at which to store temporary pdb files
        self.PDBbasename = "shiftpred"  # basename of temporary pdb files
        self.oneFile     = True         # write all frames to single or multiple pdb files
        self.PDBswritten = []           # list of pdb files written, so they can be removed later
        self.skip        = 0            # Number of frames to skip for average chemical shift prediction


        self.methods  = ['shiftx2'] # list of implemented methods

# ============================================================================ #

    def __del__(self):

        self.rm_pdb()

# ============================================================================ #
    
    def predict(self, method=None, skip=0):
        """Predict chemical shifts
        method    specifies the method to use, if None, fall back to method given at initialization
        skip      number of frames to skip for average
        """

        self.skip = int(skip)

        if method:
            self.method = method

        if not self.method in self.methods:
            print "method {} unknown.".format(self.method)
            sys.exit(1)

        if method == "shiftx2":
            self.shiftx2()

# ============================================================================ #

    def shiftx2(self):
        """Use shiftx2 to predict average chemical shifts"""

        shiftFile = '/tmp/{:s}.shiftx2'.format(random_string())

        # write frames as models in one pdb file
        PDBfiles = self.write_pdb(oneFile=True)

        # predict average chemical shifts
        try:

            cmd = ['shiftx2.py',
                   '--nmr',
                   '--infile={:s}'.format(PDBfiles[-1]),
                   '--outfile={:s}'.format(shiftFile),
                   '--outformat=CSV',
                   '--atoms=BACKBONE',
                   '--ph=7.0',
                   '--temperature=298',]
                   #'--shiftx1']

            sub.Popen(cmd).communicate()

        except OSError:
            print "shiftx2.py executable not in PATH"
            sys.exit(1)

        # read average chemical shifts

        # remove pdb file and shift file
        self.rm_pdb()

        try:
            os.remove(shiftFile)
        except OSError:
            pass
        
# ============================================================================ #

    def sparta(self):
        """Use SPARTA+ to predict average chemical shifts"""

        shiftFile = '/tmp/{:s}.shiftx2'.format(random_string())

        # write frames as models in one pdb file
        PDBfiles = self.write_pdb(oneFile=False) 

        # predict chemical shifts
        try:

            cmd = ['sparta+',
                   '-out {}'.format(),
                   '-in']

            cmd += PDBfiles

            sub.Popen(cmd).communicate()

        except OSError:
            self.rm_pdb()
            os.remove(shiftFile)
            print "sparta+ executable not in PATH"
            sys.exit(1)

        # read average chemical shifts 

        # remove pdb file and shift file
        self.rm_pdb()

        try:
            os.remove(shiftFile)
        except OSError:
            pass 

# ============================================================================ #

    def write_pdb(self, basename=None, PDBpath=None, oneFile=True):
        """Write the trajectory in self.universe as PDB files
        basename        base file name without .pdb extension
                        a number is added for multiple pdb files
        PDBpath         path at which to write trajectories
        oneFile         if True, write one file with multiple models
                        else write multiple files with one model each
        """

        rndm             = random_string(size=5)
        self.PDBbasename = "shiftPred_" + rndm  # basename of temporary pdb files
        self.PDBpath     = "/tmp"               # path at which to store temporary pdb files
        self.oneFile     = oneFile              # write all frames to single or multiple pdb files 
        PDBfilenames     = []                   # PDB filenames written during the current call of the function

        if basename:
            self.PDBbasename = basename + '_' + rndm  # basename of temporary pdb files
        if PDBpath:
            self.PDBpath     = PDBpath   # path at which to store temporary pdb files

        if self.oneFile:

            # create one writer for whole trajectory if all frames go into same file
            PDBfilename = "{}/{}.pdb".format(self.PDBpath, self.PDBbasename)
            writer = mda.Writer(PDBfilename, multiframe=True)
            self.PDBswritten.append(PDBfilename)
            PDBfilenames.append(PDBfilename)
            
            # loop through frames
            for nframe, ts in enumerate(self.universe.trajectory):

                # skip frames
                if nframe % (self.skip + 1 ) == 0:
                    writer.write(self.universe)

            writer.close_trajectory() 

        else:

            # loop through frames
            for nframe, ts in enumerate(self.universe.trajectory):

                # skip frames
                if nframe % (self.skip + 1 ) == 0:

                    # create one writer per frame if each gets its one file
                    nwrite = nframe / (self.skip + 1)
                    PDBfilename = "{}/{}_{:0>5d}.pdb".format(self.PDBpath, self.PDBbasename, nwrite)
                    writer = mda.Writer(PDBfilename)
                    self.PDBswritten.append(PDBfilename)
                    PDBfilenames.append(PDBfilename)

                    writer.write(self.universe)

                    # close each writer
                    writer.close_trajectory()

        return PDBfilenames

# ============================================================================ #

    def rm_pdb(self):
        """Remove all PDB files written by write_pdb()"""

        while True:
            try:
                os.remove(self.PDBswritten.pop())
            except OSError:
                pass
            except IndexError:
                break

# ============================================================================ #

def random_string(size=20):
    randomString = ''.join(random.choice(string.ascii_letters + string.digits) for x in range(size))
    return randomString

# ============================================================================ #
