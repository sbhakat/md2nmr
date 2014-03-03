#! /usr/bin/env python

from md2nmr import *
from NHcorr import *

name = '04_prod01_protein'
#path = '/home/oliver/externalDisk/MD_BACKUP/IL-6/IL-6_ffcomp/amber99sb-star_2IL6_13'
path = '/home/oliver/SiSc/Courses/Thesis/iGRASP/remote_fs/gpufs/externalDisk/MD_BACKUP/IL-6/IL-6_ffcomp/amber99sb-star_2IL6_13'

md = md2nmr(name, path=path)




