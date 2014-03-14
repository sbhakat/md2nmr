#! /usr/bin/env python

from md2nmr import *
from NHcorr import *
import ShiftPred as sp

name = '04_prod01_protein'
#path = '/home/oliver/externalDisk/MD_BACKUP/IL-6/IL-6_ffcomp'
#path = '/home/oliver/SiSc/Courses/Thesis/iGRASP/remote_fs/gpufs/externalDisk/MD_BACKUP/IL-6/IL-6_ffcomp'
path = '/local/jubio/oschill/iGRASP/IL-6/IL-6_ffcomp'
#path = '/home/oliver/SiSc/Courses/Thesis/iGRASP/remote_fs/jubiofs/iGRASP/IL-6/IL-6_ffcomp'

#simName = 'amber99sb-star_2IL6_13'
simName = 'gromos54a7_2IL6_23'

md = md2nmr(name, path=path+'/'+simName, rerun=False, verbose=True)
md.compute_order_parameters()


s = sp.ShiftPred(md.universe, method='sparta+')
s.skip = 1000



