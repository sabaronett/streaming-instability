#!/usr/bin/env python3
#=======================================================================
# vpmax.py
#
# Returns the maximum radial and vertical particle velocities among
# particle output (.pout) data.
#
# Author: Stanley A. Baronett
# Created: 2022-06-09
#=======================================================================

import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np

cases = ['AB', 'BA']
Pis = ['0.01', '0.02', '0.05', '0.10']
res = 2048
dirs = ['vpx', 'vpz']
tab = [['Run', 'v_p,x,max [H/T]', 'v_p,z,max [H/T]']]

for case in cases:
    for Pi in Pis:
        run = f'{case:s}/{Pi:s}'
        output = 'run'+f'/{res:n}/SI.pout.00001.dat'
        print(f'Loading {run:s} particle data...')
        time, pdata = athena_read.particles(output)
        row = [run, None, None]

        for i in range(2):
            dir = dirs[i]
            print(f'  Finding the maximum {dir:s}...')
            vpmax = np.amax(pdata[dir])
            row[i+1] = f'{vpmax:0.16e}'

        print(f'  Done.')
        tab.append(row)

fname = 'vpmaxs.txt'
print(f'Saving to {fname:s}...')
np.savetxt(fname, tab)
print(f'Done.')
