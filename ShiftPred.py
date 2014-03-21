#! /usr/bin/env python

import sys, os
import string, random
import subprocess as sub
import cPickle as pickle
import numpy as np
import MDAnalysis as mda

# ============================================================================ #

class ShiftPred:  
    """Chemical shift predictor for MD trajectories with the help of 3rd party programs"""

    def __init__(self, universe, method="sparta+"):

        self.universe = universe  # MDAnalysis trajectory
        self.method   = method    # 3rd party program to use for prediction
        
        self.PDBpath      = ""           # path at which to store temporary pdb files
        self.PDBbasename  = "shiftpred"  # basename of temporary pdb files
        self.oneFile      = True         # write all frames to single or multiple pdb files
        self.FilesWritten = []           # list of files written, so they can be removed later
        self.skip         = 0            # Number of frames to skip for average chemical shift prediction
        self.verbose      = False

        self.methods      = ['shiftx2', 'sparta+'] # list of implemented methods

        self.averageShifts = {}       # dictionary with AverageShifts objects per atom of which shifts have been predicted

# ==================================== #

    def __del__(self):

        self.rm_files()

# ==================================== #
    
    def predict(self, method=None, skip=None, verbose=None):
        """Predict chemical shifts
        method    specifies the method to use, if None, fall back to method given at initialization
        skip      number of frames to skip for average
        """

        if verbose:
            self.verbose = verbose

        if skip:
            self.skip = int(skip)

        if method:
            self.method = method

        if not self.method in self.methods:
            print "method {} unknown.".format(self.method)
            sys.exit(1)

        if self.method == "shiftx2":
            self.shiftx2()
        elif self.method == "sparta+":
            self.averageShifts = self.spartaplus()

# ==================================== #

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
        self.rm_files()

        try:
            os.remove(shiftFile)
        except OSError:
            pass
        
# ==================================== #

    def spartaplus(self):
        """Use SPARTA+ to predict average chemical shifts"""

        chunkSize = 100 # maximum number of frames to predict at once

        predictionFiles = []
        structureFiles  = []
        shiftsPerFrame  = {}  # Per atom a list of ChemShift object, one per frame

        # write frames as models in one pdb file
        PDBfiles = self.write_pdb(oneFile=False) 

        fileIndex  = 0
        while fileIndex < len(PDBfiles):

            startIndex = fileIndex
            endIndex   = fileIndex + chunkSize

            if endIndex > len(PDBfiles):
                endIndex = len(PDBfiles)

            chunk = PDBfiles[startIndex:endIndex]

            fileIndex = endIndex

            predictionFiles = []
            structureFiles  = [] 

            # predict chemical shifts
            try:

                cmd = ['sparta+',
                       '-in']
                       #'-out {}'.format(predictionFile),
                       #'-outS {}'.format(structureFile),

                cmd += chunk

                cwd = os.getcwd()

                # execute sparta+ and wait for completion
                spartaOutFilename = '/tmp/sparta_{:s}.out'.format(random_string())
                spartaOut         = open(spartaOutFilename, 'w')

                if self.verbose:
                    spartaProc        = sub.Popen(cmd)
                else:
                    spartaProc        = sub.Popen(cmd, stdout=spartaOut, stderr=sub.STDOUT)

                # store output filenames
                if len(chunk) > 1:
                    for PDBfile in chunk:
                        basename               = os.path.basename(PDBfile)
                        (filename, extension)  = os.path.splitext(basename)
                        predictionFiles.append("{}/{}_pred.tab".format(os.getcwd(), filename))
                        structureFiles.append("{}/{}_struct.tab".format(os.getcwd(), filename)) 
                else:
                    predictionFiles.append("pred.tab")
                    structureFiles.append("struct.tab")


                (stdoutdata, stderrdata) = spartaProc.communicate()
                if self.verbose:
                    print stdoutdata
                    print stderrdata

                spartaOut.close()


            # if Popen failed, complain
            except OSError:
                spartaOut.close()
                self.rm_files()
                print "sparta+ executable not in PATH"
                sys.exit(1)


            # read chemical shifts
            for filename in predictionFiles:

                shiftsThisFrame = self.spartaplus_read_shifts(filename)

                for key in shiftsThisFrame.keys():

                    # initialize list if not done yet
                    if not shiftsPerFrame.has_key(key):
                        shiftsPerFrame[key] = []

                    # add shifts of this frame to list
                    shiftsPerFrame[key].append(shiftsThisFrame[key])

            # remove pdb file and prediction and structure files
            self.rm_files(chunk)
            self.rm_files(predictionFiles)
            self.rm_files(structureFiles) 
            

#        self.rm_files()
#        sys.exit(0)


#        # predict chemical shifts
#        try:
#
#            cmd = ['sparta+',
#                   '-in']
#                   #'-out {}'.format(predictionFile),
#                   #'-outS {}'.format(structureFile),
#
#            cmd += PDBfiles
#
#            cwd = os.getcwd()
#
#            # execute spara+ and wait for completion
#            spartaOutFilename = '/tmp/sparta_{:s}.out'.format(random_string())
#            spartaOut         = open(spartaOutFilename, 'w')
#
#            if self.verbose:
#                spartaProc        = sub.Popen(cmd)
#            else:
#                spartaProc        = sub.Popen(cmd, stdout=spartaOut, stderr=sub.STDOUT)
#
#            # store output filenames
#            for PDBfile in PDBfiles:
#                basename               = os.path.basename(PDBfile)
#                (filename, extension)  = os.path.splitext(basename)
#                predictionFiles.append("{}/{}_pred.tab".format(os.getcwd(), filename))
#                structureFiles.append("{}/{}_struct.tab".format(os.getcwd(), filename)) 
#
#            (stdoutdata, stderrdata) = spartaProc.communicate()
#            if self.verbose:
#                print stdoutdata
#                print stderrdata
#
#
#        # if Popen failed, complain
#        except OSError:
#            self.rm_files()
#            print "sparta+ executable not in PATH"
#            sys.exit(1)
#
#
#        # read chemical shifts
#        shiftsPerFrame = {}  # Per atom a list of ChemShift object, one per frame
#        for filename in predictionFiles:
#
#            shiftsThisFrame = self.spartaplus_read_shifts(filename)
#
#            for key in shiftsThisFrame.keys():
#
#                # initialize list if not done yet
#                if not shiftsPerFrame.has_key(key):
#                    shiftsPerFrame[key] = []
#
#                # add shifts of this frame to list
#                shiftsPerFrame[key].append(shiftsThisFrame[key])


        # average chemical shifts over frames
        averageShifts = {}
        for key in shiftsPerFrame.keys():
            averageShifts[key] = AverageShifts(shiftsPerFrame[key])

#        # remove pdb file and prediction and structure files
#        self.rm_files()
#        self.rm_files(predictionFiles)
#        self.rm_files(structureFiles)

        return averageShifts

# ==================================== #

    def spartaplus_read_shifts(self, filename):
        """read chemical shifts from sparta plus prediction file"""

        shiftFile = open(filename, 'r')
        shiftData = shiftFile.readlines()
        shiftFile.close()

        sequence = ''

        shiftStart = 0
        firstResid = 0
        for l, line in enumerate(shiftData):
            fields = line.split()

            # read ID of first residue
            if len(fields) == 3 and fields[0] == 'DATA' and fields[1] == 'FIRST_RESID':
                firstResid = int(fields[2])

            # read sequence
            if len(fields) > 1 and fields[0] == 'DATA' and fields[1] == 'SEQUENCE':
                subsequence = ''.join(fields[2:])
                sequence += subsequence.upper()

            # get line number prediction table start
            if len(fields) > 0 and fields[0] == 'FORMAT':
                shiftStart = l + 2
                break

        # read ID of last residue
        lastResid = int(shiftData[-1].split()[0])

        shiftsPerElement = {}
        shiftsPerElement["N" ] = ChemShifts(sequence, firstResid, 'N')
        shiftsPerElement["HN"] = ChemShifts(sequence, firstResid, 'HN')
        shiftsPerElement["C" ] = ChemShifts(sequence, firstResid, 'C')
        shiftsPerElement["CA"] = ChemShifts(sequence, firstResid, 'CA')
        shiftsPerElement["CB"] = ChemShifts(sequence, firstResid, 'CB')
        shiftsPerElement["HA"] = ChemShifts(sequence, firstResid, 'HA')

        # read chemical shifts
        for line in shiftData[shiftStart:]:

            fields  = line.split()

            resid   = int(fields[0])
            resname = fields[1].upper()
            element = fields[2].upper()
            shift   = float(fields[4]) 

            # strip trailing index
            while element[-1].isdigit():
                element = element[:-1]

            # set chemical shift
            shiftsPerElement[element].set_shift(shift, resid, resname)

        return shiftsPerElement

# ==================================== #

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
        #self.PDBpath     = "/tmp"               # path at which to store temporary pdb files
        self.PDBpath     = os.getcwd()          # path at which to store temporary pdb files
        self.oneFile     = oneFile              # write all frames to single or multiple pdb files 
        PDBfilenames     = []                   # PDB filenames written during the current call of the function

        if basename:
            self.PDBbasename = basename + '_' + rndm  # basename of temporary pdb files
        if PDBpath:
            self.PDBpath     = PDBpath   # path at which to store temporary pdb files

        if self.oneFile:

            # rewind and select only protein atoms
            self.universe.trajectory.rewind()
            protein = self.universe.selectAtoms('protein')

            # create one writer for whole trajectory if all frames go into same file
            PDBfilename = "{}/{}.pdb".format(self.PDBpath, self.PDBbasename)
            writer = mda.Writer(PDBfilename, multiframe=True)
            self.FilesWritten.append(PDBfilename)
            PDBfilenames.append(PDBfilename)
            
            # loop through frames
            for nframe, ts in enumerate(self.universe.trajectory):

                # skip frames
                if nframe % (self.skip + 1 ) == 0:
                    #writer.write(self.universe)
                    writer.write(protein)

            writer.close_trajectory() 

        else:

            # rewind and select only protein atoms
            self.universe.trajectory.rewind()
            protein = self.universe.selectAtoms('protein')

            # loop through frames
            for nframe, ts in enumerate(self.universe.trajectory):

                # skip frames
                if nframe % (self.skip + 1 ) == 0:

                    # create one writer per frame if each gets its one file
                    nwrite = nframe / (self.skip + 1)
                    PDBfilename = "{}/{}_{:0>5d}.pdb".format(self.PDBpath, self.PDBbasename, nwrite)
                    writer = mda.Writer(PDBfilename)
                    self.FilesWritten.append(PDBfilename)
                    PDBfilenames.append(PDBfilename)

                    #writer.write(self.universe)
                    writer.write(protein)

                    # close each writer
                    writer.close_trajectory()

        return PDBfilenames

# ==================================== #

    def rm_files(self, filelist=None):
        """Remove all files in the given list
        If no list is given, all PDB files written by write_pdb() are deleted"""

        if filelist:
            removelist = filelist
        else:
            removelist = self.FilesWritten

        while True:
            try:
                os.remove(removelist.pop())
            except OSError:
                pass
            except IndexError:
                break

# ============================================================================ #

class ChemShifts:
    """
    Hold the chemical shifts of one atom for a single structure for all residues
    """

    def __init__(self, sequence, firstResid, atom):

        self.possibleAtoms = ['N', 'HN', 'C', 'CA', 'HA', 'CB']

        self.firstResid    = firstResid
        self.sequence      = sequence.upper()   # amino acid sequence as one letter codes in one string
        self.atom          = atom               # one atom from possible atoms for which the shifts are stored

        self.nShifts       = len(self.sequence)                     # number of chemical shifts for the given atom
        self.shifts        = np.zeros(self.nShifts) * np.nan        # shifts in ppm
        self.resids        = np.zeros(self.nShifts, dtype=np.int64) # residue IDs of chemical shifts
        self.resn          = self.sequence                          # residue names as one   letter codes
        self.resname       = []                                     # residue names as three letter codes

        # amino acid dictionaries

        self.aa3to1 = {
                       "Ala": "A",
                       "Cys": "C",
                       "Asp": "D",
                       "Glu": "E",
                       "Phe": "F",
                       "Gly": "G",
                       "His": "H",
                       "Ile": "I",
                       "Lys": "K",
                       "Leu": "L",
                       "Met": "M",
                       "Asn": "N",
                       "Pro": "P",       
                       "Gln": "Q",
                       "Arg": "R",
                       "Ser": "S",
                       "Thr": "T",
                       "Val": "V",
                       "Trp": "W",
                       "Tyr": "Y"}

        self.aa1to3 = {
                       "A": "Ala",
                       "C": "Cys",
                       "D": "Asp",
                       "E": "Glu",
                       "F": "Phe",
                       "G": "Gly",
                       "H": "His",
                       "I": "Ile",
                       "K": "Lys",
                       "L": "Leu",
                       "M": "Met",
                       "N": "Asn",
                       "P": "Pro",
                       "Q": "Gln",
                       "R": "Arg",
                       "S": "Ser",
                       "T": "Thr",
                       "V": "Val",
                       "W": "Trp",
                       "Y": "Tyr"}

        if self.atom not in self.possibleAtoms:
            print "Chemical shifts for {} not implemented".format(self.atom)
            sys.exit(1)

        resid = self.firstResid
        for letter in self.resn:
            self.resids[resid-firstResid] = resid
            self.resname.append(self.aa1to3[letter])
            resid += 1

# ==================================== #

#    def set_length(self):
#        """Determine number of shifts to store from sequence
#        no N  shift for: first res, Pro
#        no HN shift for: first res, Pro
#        no CB shift for: Gly
#        no C  shift for: last res"""
#
#        # count number of prolines and glycines in sequence
#        nPro = self.sequence.count('P')
#        nGly = self.sequence.count('G')
#
#        self.nShifts = len(self.sequence)
#
#        # determine proper number of shifts for element
#        if self.atom == 'N' or self.atom == 'HN':
#            self.nShifts -= 1    # first residue
#            self.nShifts -= nPro # prolines
#        elif self.atom == 'CB':
#            self.nShifts -= nGly # glycines
#        elif self.atom == 'C':
#            self.nShifts -= 1    # last residue
#        
#        # initialize data structures
#        self.shifts  = np.zeros(self.nShifts, dtype=np.float64)
#        self.resids  = np.zeros(self.nShifts, dtype=np.int64)
#        self.nextPos = 0

# ==================================== #

    def set_shift(self, shift, resid, resname):
        """set next chemical shift"""

#        # set resname
#        if len(resname) == 1:
#            self.resn   .append(resname.upper())
#            self.resname.append(self.aa1to3[resname])
#
#        elif len(resname) == 3:
#            self.resn   .append(self.aa3to1[resname])
#            self.resname.append(resname)
#
#        else:
#            print "Residue name unknown {}.".format(resname)
#            sys.exit(1)

        # set resid and shift
        position = resid - self.firstResid
        self.resids[position] = resid
        self.shifts[position] = shift 


# ============================================================================ #

class AverageShifts:
    """Compute average chemical shifts for the given element"""

    def __init__(self, shiftList):

        self.possibleAtoms = ['N', 'HN', 'C', 'CA', 'HA', 'CB']

        self.sequence      = ""          # amino acid sequence as one letter codes in one string
        self.atom          = ""          # one atom from possible atoms for which the shifts are stored

        self.shiftList     = shiftList   # list of ChemShift objects
        self.nShifts       = 0           # number of chemical shifts for the given atom
        self.shifts        = np.zeros(0) # average shifts in ppm
        self.resids        = np.zeros(0) # residue IDs of chemical shifts
        self.resn          = []          # residue names as one   letter codes
        self.resname       = []          # residue names as three letter codes

        self.average()

# ==================================== #

    def average(self):
        """Compute average of shifts"""

        # copy relevant data
        try:
            representative = self.shiftList[0]

            self.sequence = representative.sequence
            self.atom     = representative.atom
            self.nShifts  = representative.nShifts
            self.shifts   = np.zeros_like(representative.shifts)
            self.resids   = representative.resids
            self.resn     = representative.resn
            self.resname  = representative.resname

        except IndexError:
            print "shiftList is empty"
            sys.exit(1)

        except AttributeError:
            print "Whatever is in shiftList, its not a ChemShift object"
            sys.exit(1)

        # compute average
        for frame in self.shiftList:
            self.shifts += frame.shifts

        self.shifts /= len(self.shiftList)
 
# ============================================================================ #

def random_string(size=20):
    randomString = ''.join(random.choice(string.ascii_letters + string.digits) for x in range(size))
    return randomString

# ============================================================================ #
